/*
 *  fem.c
 *  Library for MECA1120 : Finite Elements for dummies
 *  
 *  Basé sur le travail de :
 *  Copyright (C) 2017 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

static const double _gaussTri3Xsi[3]     = { 0.166666666666667, 0.666666666666667, 0.166666666666667};
static const double _gaussTri3Eta[3]     = { 0.166666666666667, 0.166666666666667, 0.666666666666667};
static const double _gaussTri3Weight[3]  = { 0.166666666666667, 0.166666666666667, 0.166666666666667};


femIntegration *femIntegrationCreate(int n, femElementType type)
{
    femIntegration *theRule = malloc(sizeof(femIntegration));
    if (type == FEM_TRIANGLE && n == 3) {
        theRule->n      = 3;
        theRule->xsi    = _gaussTri3Xsi;
        theRule->eta    = _gaussTri3Eta;
        theRule->weight = _gaussTri3Weight; }
    else Error("Cannot create such an integration rule !");
    return theRule; 
}

void femIntegrationFree(femIntegration *theRule)
{
    free(theRule);
}

void _p1c0_x(double *xsi, double *eta) 
{
    xsi[0] =  0.0;  eta[0] =  0.0;
    xsi[1] =  1.0;  eta[1] =  0.0;
    xsi[2] =  0.0;  eta[2] =  1.0;
}

void _p1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = 1 - xsi - eta;  
    phi[1] = xsi;
    phi[2] = eta;
}

void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] = -1.0;  
    dphidxsi[1] =  1.0;
    dphidxsi[2] =  0.0;
    dphideta[0] = -1.0;  
    dphideta[1] =  0.0;
    dphideta[2] =  1.0;

}


femDiscrete *femDiscreteCreate(int n, femElementType type)
{
    femDiscrete *theSpace = malloc(sizeof(femDiscrete));
    if (type == FEM_TRIANGLE && n == 3) {
        theSpace->n       = 3;
        theSpace->x2      = _p1c0_x;
        theSpace->phi2    = _p1c0_phi;
        theSpace->dphi2dx = _p1c0_dphidx; }
    else Error("Cannot create such a discrete space !");
    return theSpace; 
}

void femDiscreteFree(femDiscrete *theSpace)
{
    free(theSpace);
}

void femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta)
{
    mySpace->x2(xsi,eta);
}

void femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi)
{
    mySpace->phi2(xsi,eta,phi);
}

void femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta)
{
    mySpace->dphi2dx(xsi,eta,dphidxsi,dphideta);
}

void femDiscretePrint(femDiscrete *mySpace)
{
    int i,j;
    int n = mySpace->n;
    double xsi[4], eta[4], phi[4], dphidxsi[4], dphideta[4];
    
    femDiscreteXsi2(mySpace,xsi,eta);
    for (i=0; i < n; i++) {
        
        femDiscretePhi2(mySpace,xsi[i],eta[i],phi);
        femDiscreteDphi2(mySpace,xsi[i],eta[i],dphidxsi,dphideta);

        for (j=0; j < n; j++)  {
            printf("(xsi=%+.1f,eta=%+.1f) : ",xsi[i],eta[i]);
            printf(" phi(%d)=%+.1f",j,phi[j]);  
            printf("   dphidxsi(%d)=%+.1f",j,dphidxsi[j]);  
            printf("   dphideta(%d)=%+.1f \n",j,dphideta[j]);  }
        printf(" \n"); }
}



femMesh *femMeshRead(const char *filename)
{
    femMesh *theMesh = malloc(sizeof(femMesh));

    int i,trash,*elem;
    
    FILE* file = fopen(filename,"r");
    if (file == NULL) Error("No mesh file !");

    ErrorScan(fscanf(file, "Number of nodes %d \n", &theMesh->nNode));
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        ErrorScan(fscanf(file,"%d : %le %le \n",&trash,&theMesh->X[i],&theMesh->Y[i])); }
    
    char str[256]; if (fgets(str, sizeof(str), file) == NULL) Error("Corrupted mesh file !");

    if (!strncmp(str,"Number of triangles",19))  { 
        ErrorScan(sscanf(str,"Number of triangles %d \n", &theMesh->nElem));
        theMesh->elem = malloc(sizeof(int)*3*theMesh->nElem);
        theMesh->nLocalNode = 3;
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*3]);
            ErrorScan(fscanf(file,"%d : %d %d %d\n", &trash,&elem[0],&elem[1],&elem[2])); }}
    else{
		printf("Corrupted mesh file, only triangles are allowed");
		fclose(file);
		return NULL;
	}
  
    fclose(file);
    return theMesh;
}

void femMeshFree(femMesh *theMesh)
{
    free(theMesh->X);
    free(theMesh->Y);
    free(theMesh->elem);
    free(theMesh);
}

void femMeshWrite(const femMesh *theMesh, const char *filename)
{
    int i,*elem;
    
    FILE* file = fopen(filename,"w");
    
    fprintf(file, "Number of nodes %d \n", theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        fprintf(file,"%6d : %14.7e %14.7e \n",i,theMesh->X[i],theMesh->Y[i]); }
    
    if (theMesh->nLocalNode == 4) {
        fprintf(file, "Number of quads %d \n", theMesh->nElem);  
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*4]);
            fprintf(file,"%6d : %6d %6d %6d %6d \n", i,elem[0],elem[1],elem[2],elem[3]);   }}
    else if (theMesh->nLocalNode == 3) {
        fprintf(file, "Number of triangles %d \n", theMesh->nElem);  
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*3]);
            fprintf(file,"%6d : %6d %6d %6d \n", i,elem[0],elem[1],elem[2]);   }}    
    fclose(file);
}
                                                                                                          
femEdges *femEdgesCreate(femMesh *theMesh)
{
    femEdges *theEdges = malloc(sizeof(femEdges));
    int nLoc = theMesh->nLocalNode;
    int i,j,n = theMesh->nElem * nLoc;
    femEdge* edges = malloc(n * sizeof(femEdge));
    theEdges->mesh  = theMesh;
    theEdges->edges = edges;
    theEdges->nEdge = n;
    theEdges->nBoundary = n;
    
    for (i = 0; i < theMesh->nElem; i++) {
        int *elem = &(theMesh->elem[i*nLoc]);
        for (j = 0; j < nLoc; j++) {
            int id = i * nLoc + j;
            edges[id].elem[0] = i;
            edges[id].elem[1] = -1;
            edges[id].node[0] = elem[j];
            edges[id].node[1] = elem[(j + 1) % nLoc]; }}

    qsort(theEdges->edges, theEdges->nEdge, sizeof(femEdge), femEdgesCompare);

    int index = 0;          
    int nBoundary = 0;
    
    for (i=0; i < theEdges->nEdge; i++) {
      if (i == theEdges->nEdge - 1 || femEdgesCompare(&edges[i],&edges[i+1]) != 0) {
              edges[index] = edges[i];
              nBoundary++; }
      else {  edges[index] = edges[i];
              edges[index].elem[1] = edges[i+1].elem[0];
              i = i+1;}
      index++; }
      
    theEdges->edges = realloc(edges, index * sizeof(femEdge));
    theEdges->nEdge = index;
    theEdges->nBoundary = nBoundary;
    return theEdges;
}

void femEdgesPrint(femEdges *theEdges)
{
    int i;    
    for (i = 0; i < theEdges->nEdge; ++i) {
        printf("%6d : %4d %4d : %4d %4d \n",i,
               theEdges->edges[i].node[0],theEdges->edges[i].node[1],
               theEdges->edges[i].elem[0],theEdges->edges[i].elem[1]); }
}

void femEdgesFree(femEdges *theEdges)
{
    free(theEdges->edges);
    free(theEdges);
}

int femEdgesCompare(const void *edgeOne, const void *edgeTwo)
{
    int *nodeOne = ((femEdge*) edgeOne)->node;
    int *nodeTwo = ((femEdge*) edgeTwo)->node;  
    int  diffMin = fmin(nodeOne[0],nodeOne[1]) - fmin(nodeTwo[0],nodeTwo[1]);
    int  diffMax = fmax(nodeOne[0],nodeOne[1]) - fmax(nodeTwo[0],nodeTwo[1]);
    
    if (diffMin < 0)    return  1;
    if (diffMin > 0)    return -1;
    if (diffMax < 0)    return  1;
    if (diffMax > 0)    return -1; 
                        return  0;
}

femBandSystem *femBandSystemCreate(int size, int band)
{
    femBandSystem *myBandSystem = malloc(sizeof(femBandSystem));
    myBandSystem->B = malloc(sizeof(double)*size*(band+1));
    myBandSystem->A = malloc(sizeof(double*)*size);        
    myBandSystem->size = size;
    myBandSystem->band = band;
    myBandSystem->A[0] = myBandSystem->B + size;
    int i;
    for (i=1 ; i < size ; i++) 
        myBandSystem->A[i] = myBandSystem->A[i-1] + band - 1;
    femBandSystemInit(myBandSystem);
    return(myBandSystem);
}
 
void femBandSystemFree(femBandSystem *myBandSystem)
{
    free(myBandSystem->B);
    free(myBandSystem->A); 
    free(myBandSystem);
}
 
void femBandSystemInit(femBandSystem *myBandSystem)
{
    int i;
    int size = myBandSystem->size;
    int band = myBandSystem->band;
    for (i=0 ; i < size*(band+1) ; i++) 
        myBandSystem->B[i] = 0;        
}
   
void femBandSystemPrintInfos(femBandSystem *myBand)
{
    int size = myBand->size;
    int band = myBand->band;
    printf(" \n");
    printf("    Banded Gaussian elimination \n");
    printf("    Storage informations \n");
    printf("    Matrix size      : %8d\n",size);
    printf("    Matrix band      : %8d\n",band);
    printf("    Bytes required   : %8d\n",(int)sizeof(double)*size*(band+1));     
}

void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    int i,j;
    for (i = 0; i < nLoc; i++) { 
        int myRow = map[i];
        for(j = 0; j < nLoc; j++) {
            int myCol = map[j];
            if (myCol >= myRow)  myBandSystem->A[myRow][myCol] += Aloc[i*nLoc+j]; }
        myBandSystem->B[myRow] += Bloc[i]; }
}
 
void femBandSystemConstrain(femBandSystem *myBand, int myNode, double myValue) 
{
    double  **A, *B;
    int     i, size, band, ifirst, iend;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    ifirst = fmax(0,myNode - band + 1);
    iend   = myNode;
    for (i=ifirst; i < iend; i++) {
        B[i] -= myValue * A[i][myNode];
        A[i][myNode] = 0; }

    ifirst = myNode+1;
    iend = fmin(myNode + band,size);
    for (i=ifirst; i < iend; i++) {
        B[i] -= myValue * A[myNode][i];
        A[myNode][i] = 0; }
        
    A[myNode][myNode] = 1;
    B[myNode] = myValue;
}
 
double  *femBandSystemEliminate(femBandSystem *myBand)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    /* Incomplete Cholesky factorization */ 

    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-4 ) {
            Error("Cannot eleminate with such a pivot"); }
        jend = fmin(k + band,size);
        for (i = k+1 ; i <  jend; i++) {
            factor = A[k][i] / A[k][k];
            for (j = i ; j < jend; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
        
    /* Back-substitution */

    for (i = (size-1); i >= 0 ; i--) {
        factor = 0;
        jend = fmin(i + band,size);
        for (j = i+1 ; j < jend; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }

    return(myBand->B);
}

//Checks whether the grain at position (x,y) is in the triangle iElem
int elemContains(double x, double y, femMesh *theMesh, int iElem){    
    int jacobian[3];
    double x2, x3, y2, y3;
    int tab[3] = {1,2,0};
    for(int i = 0; i < 3; i++){

        x2 = theMesh->X[theMesh->elem[iElem*3+i]];
        y2 = theMesh->Y[theMesh->elem[iElem*3+i]];
        x3 = theMesh->X[theMesh->elem[iElem*3+tab[i]]];
        y3 = theMesh->Y[theMesh->elem[iElem*3+tab[i]]];
        jacobian[i] = (((x2 - x) * (y3 - y) - (x3 - x) * (y2 - y)) >= 0);

    }
    if(jacobian[0] == jacobian[1] && jacobian[1] == jacobian[2]){
        return 1;
    }
    return 0;
}

double femMin(double *x, int n) 
{
    double myMin = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMin = fmin(myMin,x[i]);
    return myMin;
}

double femMax(double *x, int n)
{
    double myMax = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMax = fmax(myMax,x[i]);
    return myMax;
}

void femError(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);                                                 
}

void femErrorScan(int test, int line, char *file)                                  
{ 
    if (test >= 0)  return;
    
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in fscanf or fgets in %s at line %d : \n", file, line);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");   
    exit(69);                                       
}


void femWarning(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Warning in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");                                              
}
