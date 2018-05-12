

#include"fem.h"

double *theGlobalCoord;

int compare(const void *nodeOne, const void *nodeTwo) 
{
    int *iOne = (int *)nodeOne;
    int *iTwo = (int *)nodeTwo;
    double diff = theGlobalCoord[*iOne] - theGlobalCoord[*iTwo];
    if (diff < 0)    return  1;
    if (diff > 0)    return -1;
    return  0;  
}

void femDiffusionRenumber(femDiffusionProblem *theProblem, femRenumType renumType)
{
    int i, *inverse;
    
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < theProblem->mesh->nNode; i++) 
                theProblem->number[i] = i;
            break;
        case FEM_XNUM : 
            inverse = malloc(sizeof(int)*theProblem->mesh->nNode);
            for (i = 0; i < theProblem->mesh->nNode; i++) 
                inverse[i] = i; 
            theGlobalCoord = theProblem->mesh->X;
            qsort(inverse, theProblem->mesh->nNode, sizeof(int), compare);
            for (i = 0; i < theProblem->mesh->nNode; i++)
                theProblem->number[inverse[i]] = i;
            free(inverse);  
            break;
        case FEM_YNUM : 
            inverse = malloc(sizeof(int)*theProblem->mesh->nNode);
            for (i = 0; i < theProblem->mesh->nNode; i++) 
                inverse[i] = i; 
            theGlobalCoord = theProblem->mesh->Y;
            qsort(inverse, theProblem->mesh->nNode, sizeof(int), compare);
            for (i = 0; i < theProblem->mesh->nNode; i++)
                theProblem->number[inverse[i]] = i;
            free(inverse);  
            break;
        default : Error("Unexpected renumbering option"); }
}

int femDiffusionComputeBand(femDiffusionProblem *theProblem)
{
    femMesh *theMesh = theProblem->mesh;
    int iElem,j,myMax,myMin,myBand,map[4];
    int nLocal = theMesh->nLocalNode;
    myBand = 0;
    for(iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; ++j) 
            map[j] = theProblem->number[theMesh->elem[iElem*nLocal+j]];
        myMin = map[0];
        myMax = map[0];
        for (j=1; j < nLocal; j++) {
            myMax = fmax(map[j],myMax);
            myMin = fmin(map[j],myMin); }
        if (myBand < (myMax - myMin)) myBand = myMax - myMin; }         
    return(++myBand);
}

void femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, 
                                double *Bloc, double *Uloc, int *map, int nLoc)
{
    int i,j,myRow;
    if (mySolver->iter == 0) {
        for (i = 0; i < nLoc; i++) { 
            myRow = map[i];
            mySolver->R[myRow] -= Bloc[i];
            for(j = 0; j < nLoc; j++) {
                mySolver->R[myRow] += Aloc[i*nLoc+j]*Uloc[j];}}}
    
    for (i = 0; i < nLoc; i++) {
        myRow = map[i];
        for(j = 0; j < nLoc; j++) {
            int myCol = map[j];
            mySolver->S[myRow] += Aloc[i*nLoc+j] * mySolver->D[myCol]; }}
}

void femIterativeSolverConstrain(femIterativeSolver* mySolver, int myNode, double myValue) 
{

//
// Pour conserver les conditions essentielles homogenes,
// il suffit d'annuler les deux increments pour le residu et la direction :-)

// Observer qu'il est impossible de resoudre le probleme avec des conditions ess. non-homogenes
// tel que c'etait implemente dans le devoir :-)
// Il faudrait modifier tres legerement la fonction femDiffusionCompute 
// pour que cela soit possible :-)
//

        mySolver->R[myNode] = myValue;
        mySolver->S[myNode] = myValue; 
}

double *femIterativeSolverEliminate(femIterativeSolver *mySolver)
{


// R   = vecteur R (residu)
// X   = dX
// D   = vecteur D (direction)
// S   = vecteur A D (direction conjugee)

// Demarrage de l'algorithme
// La direction initiale est le residu et on initialise l'increment a zero.

    mySolver->iter++;
    double error = 0.0; int i;
    for (i=0; i < mySolver->size; i++) {
        error += (mySolver->R[i])*(mySolver->R[i]); }
    
  
    if (mySolver->iter == 1) {
        for (i=0; i < mySolver->size; i++) {
            mySolver->X[i] = 0; 
            mySolver->D[i] = mySolver->R[i]; }}           
    else {    
        double denAlpha = 0.0;
        for (i=0; i < mySolver->size; i++) {
            denAlpha += mySolver->R[i] * mySolver->S[i]; }

        double alpha = -error/denAlpha;
    
        double numBeta = 0.0;
        for (i=0; i < mySolver->size; i++) {            
            mySolver->R[i] = mySolver->R[i] + alpha * mySolver->S[i];
            numBeta += mySolver->R[i] * mySolver->R[i]; }
            
        double beta = numBeta/error;
        for (i=0; i < mySolver->size; i++) {
            mySolver->X[i] = alpha * mySolver->D[i];
            mySolver->D[i] = mySolver->R[i] + beta * mySolver->D[i]; 
            mySolver->S[i] = 0.0; }}
     
    mySolver->error = sqrt(error);
    return(mySolver->X);
}



