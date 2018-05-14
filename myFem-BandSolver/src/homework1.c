
#include"fem.h"

femDiffusionProblem *theProb;

int compareX (void const *a, void const *b){
	int *a1 = (int *) a;
	int *b1 = (int *) b;
	double *X = theProb->mesh->X;
	return  X[*a1] - X[*b1];	
}
int compareY (void const *a, void const *b){
	int *a1 = (int *) a;
	int *b1 = (int *) b;
	double *Y = theProb->mesh->Y;
	return  Y[*a1] - Y[*b1];
}

void femDiffusionRenumber(femDiffusionProblem *theProblem, femRenumType renumType)
{
    int i;    
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < theProblem->mesh->nNode; i++) 
                theProblem->number[i] = i;
            break;
// 
// A modifier :-)
// debut
//
        case FEM_XNUM :
			for (i = 0; i < theProblem->mesh->nNode; i++) 
                theProblem->number[i] = i;
			theProb = theProblem;
			qsort(theProblem->number, theProblem->size, sizeof(int), compareX);
			for (i = 0; i < theProblem->mesh->nNode; i++) 
                theProblem->number[theProblem->number[i]] = i;
			break;			
        case FEM_YNUM : 
            for (i = 0; i < theProblem->mesh->nNode; i++) 
                theProblem->number[i] = i;
			theProb = theProblem;
			qsort(theProblem->number, theProblem->size, sizeof(int), compareY);
			for (i = 0; i < theProblem->mesh->nNode; i++) 
                theProblem->number[theProblem->number[i]] = i;
            break;           
// 
// end
//

        default : Error("Unexpected renumbering option"); }
}

double minhArray(int *x, int n) 
{
    int myMinh = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMinh = fmin(myMinh,x[i]);
    return myMinh;
}

double maxArray(int *x, int n)
{
    int myMax = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMax = fmax(myMax,x[i]);
    return myMax;
}

int femDiffusionComputeBand(femDiffusionProblem *theProblem)
{	
    femMesh *theMesh = theProblem->mesh;    
	int nLoc = theMesh->nLocalNode;
	int i,j,numEdge[nLoc], locMinh, locMax, locDiff;
	int myBand = 0;		
	for(i = 0; i < theMesh->nElem; i++){
		for(j = 0; j < nLoc; j++){
			numEdge[j] = theProblem->number[theMesh->elem[i*nLoc + j]];			
		}
		locMinh = minhArray(numEdge, nLoc);
		locMax = maxArray(numEdge, nLoc);
		locDiff = locMax - locMinh;
		if(locDiff > myBand){
			myBand = locDiff;
		}
	}		
    return(myBand+1);
}

void femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc)
{
    int i,j;
    for (i = 0; i < nLoc; i++) { 
        int myRow = map[i];
        for(j = 0; j < nLoc; j++) {
            mySolver->R[myRow] -= Aloc[i*nLoc+j]*Uloc[j]; }
        mySolver->R[myRow] += Bloc[i]; }
	
}

void femIterativeSolverConstrain(femIterativeSolver* mySolver, int myNode, double myValue) 
{   
    mySolver->R[myNode] = myValue; 
	mySolver->S[myNode] = myValue; 
}

double *femIterativeSolverEliminate(femIterativeSolver *mySolver)
{
    mySolver->iter++;	
	
    double error = 0.0; int i;
	double alphaDen = 0.0;
	double newError = 0.0;
    for (i=0; i < mySolver->size; i++) {
        error += (mySolver->R[i])*(mySolver->R[i]); // alphaNum = error
		//alphaDen += (mySolver->R[i])*(mySolver->S[i]);
        		
	}
	if (mySolver->iter == 1) {			
        for (i=0; i < mySolver->size; i++) {
            mySolver->X[i] = 0; 
            mySolver->D[i] = mySolver->R[i];
		}
	}
else{
	for (i=0; i < mySolver->size; i++) {
		alphaDen += mySolver->R[i] * mySolver->S[i]; }
	double alpha = - error / alphaDen;
	int j;
	double betaNum = 0.0;
	for(j = 0; j < mySolver->size; j++){
		mySolver->R[j] = mySolver->R[j] + alpha * mySolver->S[j];
		betaNum += mySolver->R[j] * mySolver->R[j]; // rk+1 * rk+1
		mySolver->X[j] = alpha * mySolver->D[j];
	}
	double beta = betaNum / error;
	int k;
	for(k = 0; k < mySolver->size; k++){
		mySolver->D[k] = mySolver->R[k] + beta * mySolver->D[k];		
	}
}
    mySolver->error = sqrt(error);
    return(mySolver->X);
}
