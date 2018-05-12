
#include"fem.h"

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
        case FEM_YNUM : 
            for (i = 0; i < theProblem->mesh->nNode; i++) 
                theProblem->number[i] = i;
            break;            
// 
// end
//

        default : Error("Unexpected renumbering option"); }
}

int femDiffusionComputeBand(femDiffusionProblem *theProblem)
{
    femMesh *theMesh = theProblem->mesh;
    int myBand = theMesh->nNode;
    return(myBand);
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
}

double *femIterativeSolverEliminate(femIterativeSolver *mySolver)
{
    mySolver->iter++;
    double error = 0.0; int i;
    for (i=0; i < mySolver->size; i++) {
        error += (mySolver->R[i])*(mySolver->R[i]);
        mySolver->R[i] = mySolver->R[i]/5;  }
    mySolver->error = sqrt(error);
    return(mySolver->R);
}

