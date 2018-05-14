//authors: Alexandre Gobeaux and Gildas Mulders
#include"fem.h"

femDiffusionProblem *theProb;

int compareX (void const *a, void const *b){
	int *a1 = (int *) a;
	int *b1 = (int *) b;
	double *X = theProb->mesh->X;
	if( X[*a1] == X[*b1] ){
		return 0;
	}
	else if( X[*a1] < X[*b1] ){
		return -1;
	}
	else {
		return 1;
	}	 
}
int compareY (void const *a, void const *b){
	int *a1 = (int *) a;
	int *b1 = (int *) b;
	double *Y = theProb->mesh->Y;	
	if( Y[*a1] == Y[*b1] ){
		return 0;
	}
	else if( Y[*a1] < Y[*b1] ){
		return -1;
	}
	else {
		return 1;
	}	
}

void femDiffusionRenumber(femDiffusionProblem *theProblem, femRenumType renumType)
{
    int i;    
	int nbrIndex[theProblem->size];
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
			for( i = 0; i < theProblem->size; i++ ) 
                nbrIndex[i] = i;
			theProb = theProblem;
			qsort(nbrIndex, theProblem->size, sizeof(int), compareX);
			for( i = 0; i < theProblem->size; i++ ) 
                theProblem->number[nbrIndex[i]] = i;
			break;			
        case FEM_YNUM : 
            for( i = 0; i < theProblem->size; i++ ) 
                nbrIndex[i] = i;
			theProb = theProblem;
			qsort(nbrIndex, theProblem->size, sizeof(int), compareY);
			for( i = 0; i < theProblem->size; i++ ) 
				theProblem->number[nbrIndex[i]] = i;
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
    for( i=1 ;i < n; i++ ) 
        myMinh = fmin(myMinh,x[i]);
    return myMinh;
}

double maxArray(int *x, int n)
{
    int myMax = x[0];
    int i;
    for( i=1 ;i < n; i++ ) 
        myMax = fmax(myMax,x[i]);
    return myMax;
}

int femDiffusionComputeBand(femDiffusionProblem *theProblem)
{	
    femMesh *theMesh = theProblem->mesh;    
	int nLoc = theMesh->nLocalNode;
	int i,j,numEdge[nLoc], locMinh, locMax, locDiff;
	int myBand = 0;		
	for( i = 0; i < theMesh->nElem; i++ ){
		for( j = 0; j < nLoc; j++ ){
			numEdge[j] = theProblem->number[theMesh->elem[i*nLoc + j]];			
		}
		locMinh = minhArray(numEdge, nLoc);
		locMax = maxArray(numEdge, nLoc);
		locDiff = locMax - locMinh;
		if( locDiff > myBand ){
			myBand = locDiff;
		}
	}		
    return myBand+1;
}