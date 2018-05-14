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