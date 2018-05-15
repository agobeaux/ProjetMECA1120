/*
 *  fem.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
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
 
void femBandSystemPrint(femBandSystem *myBand)
{
    double  **A, *B;
    int     i, j, band, size;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    for (i=0; i < size; i++) {
        for (j=i; j < i+band; j++)
            if (A[i][j] == 0) printf("         ");   
            else              printf(" %+.1e",A[i][j]);
        printf(" :  %+.1e \n",B[i]); }
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


double femBandSystemGet(femBandSystem* myBandSystem, int myRow, int myCol)
{
    double value = 0;
    if (myCol >= myRow && myCol < myRow+myBandSystem->band)  value = myBandSystem->A[myRow][myCol]; 
    return(value);
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

femCouetteProblem *femCouetteCreate(const char *filename,  femRenumType renumType)
{
    int i,band;
    femCouetteProblem *theProblem = malloc(sizeof(femCouetteProblem));
    theProblem->mesh  = femMeshRead(filename);
    if(theProblem->mesh == NULL){
		return NULL;
	}      
    theProblem->edges = femEdgesCreate(theProblem->mesh); 
    if (theProblem->mesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE); 
    }
    theProblem->size = theProblem->mesh->nNode;
    theProblem->number  = malloc(sizeof(int)*theProblem->size);
    femCouetteRenumber(theProblem,renumType);
    band = femCouetteComputeBand(theProblem);
    theProblem->systemX = femBandSystemCreate(theProblem->size, band);
    theProblem->systemY = femBandSystemCreate(theProblem->size, band);
    theProblem->soluceX = malloc(sizeof(double)*theProblem->size);
    theProblem->soluceY = malloc(sizeof(double)*theProblem->size);
    for (i = 0; i < theProblem->size; i++)
    {      
        theProblem->soluceX[i] = 0;
        theProblem->soluceY[i] = 0;
    }
 
    return theProblem;
}

void femSoluceInit(femCouetteProblem *theProblem)
{
    for (int i = 0; i < theProblem->size; i++)
    {      
        theProblem->soluceX[i] = 0;
        theProblem->soluceY[i] = 0;
    }
}

void femCouetteFree(femCouetteProblem *theProblem)
{
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    femEdgesFree(theProblem->edges);
    femMeshFree(theProblem->mesh);
    femBandSystemFree(theProblem->systemX);
    femBandSystemFree(theProblem->systemY);
    free(theProblem->number);
    free(theProblem->soluceX);
    free(theProblem->soluceY);
    free(theProblem);
}
    

void femCouetteMeshLocal(const femCouetteProblem *theProblem, const int iElem, int *map, double *x, double *y)
{
    femMesh *theMesh = theProblem->mesh;
    int j,nLocal = theMesh->nLocalNode;
    
    for (j=0; j < nLocal; ++j) {
        map[j] = theMesh->elem[iElem*nLocal+j];
        x[j]   = theMesh->X[map[j]];
        y[j]   = theMesh->Y[map[j]]; 
        //u[j]   = theProblem->soluce[map[j]];
        map[j] = theProblem->number[map[j]];
        }   
}

void femCouetteCompute(femCouetteProblem *theProblem, femGrains *theGrains, double mu, double gamma, double vExt, int systIsY)
{
    femMesh *theMesh = theProblem->mesh;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;
    femBandSystem *theSystem;
    if(systIsY)
    {
        theSystem = theProblem->systemY;
    }
    else
    {
        theSystem = theProblem->systemX;    
    }
    femEdges *theEdges = theProblem->edges;
    int *number = theProblem->number;
       
    if (theSpace->n > 4) Error("Unexpected discrete space size !"); 
    
    double Xloc[4],Yloc[4],phi[4],dphidxsi[4],phiGrains[3], dphideta[4],dphidx[4],dphidy[4],Aloc[16],Bloc[4],xGrains, yGrains, xsiGrains, etaGrains, vGrains;
    int iEdge,iElem,iInteg,i,j,map[4], iGrains;
   
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (i = 0; i < theSpace->n; i++)      Bloc[i] = 0;
        for (i = 0; i < (theSpace->n)*(theSpace->n); i++) Aloc[i] = 0;
        femCouetteMeshLocal(theProblem,iElem,map,Xloc,Yloc);  
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0; 
            double dydeta = 0;
            for (i = 0; i < theSpace->n; i++) {    
                dxdxsi += Xloc[i]*dphidxsi[i];       
                dxdeta += Xloc[i]*dphideta[i];   
                dydxsi += Yloc[i]*dphidxsi[i];   
                dydeta += Yloc[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    Aloc[i*(theSpace->n)+j] += mu*(dphidx[i] * dphidx[j] 
                                            + dphidy[i] * dphidy[j]) * jac * weight; 
                }
            }                                                                                            
        }        
        for(iGrains = 0; iGrains < theGrains->n; iGrains++)
        {            
            xGrains = theGrains->x[iGrains];
            yGrains = theGrains->y[iGrains]; 
            if(elemContains(xGrains,yGrains,theMesh,iElem)==1)
            {
                xsiGrains = -(Xloc[0] * (Yloc[2] - yGrains) + Xloc[2] * (yGrains - Yloc[0]) + xGrains * (Yloc[0] - Yloc[2]))/(Xloc[0] * (Yloc[1] - Yloc[2]) + Xloc[1] * (Yloc[2] - Yloc[0]) + Xloc[2] * (Yloc[0] - Yloc[1]));
                etaGrains = (Xloc[0] * (Yloc[1] - yGrains) + Xloc[1] * (yGrains - Yloc[0]) + xGrains * (Yloc[0] - Yloc[1]))/(Xloc[0] * (Yloc[1] - Yloc[2]) + Xloc[1] * (Yloc[2] - Yloc[0]) + Xloc[2] * (Yloc[0] - Yloc[1]));
                femDiscretePhi2(theSpace,xsiGrains,etaGrains,phiGrains);
                if(systIsY)
                {
                    vGrains = theGrains->vy[iGrains];
                }
                else
                {
                    vGrains = theGrains->vx[iGrains];
                }
                for (i = 0; i < theSpace->n; i++) 
                { 
                    Bloc[i] += gamma*phiGrains[i]*vGrains; 
                    for(j = 0; j < theSpace->n; j++) 
                    {
                        Aloc[i*(theSpace->n)+j] += gamma*phiGrains[i]*phiGrains[j];
                    }
                }
            }
        }
        femBandSystemAssemble(theSystem, Aloc,Bloc,map,theSpace->n);
    } 

    for (iEdge= 0; iEdge < theEdges->nEdge; iEdge++) {      
        if (theEdges->edges[iEdge].elem[1] < 0) {  
            for (i = 0; i < 2; i++) 
            {
                double v;
                int iNode = theEdges->edges[iEdge].node[i];
                double xloc = theMesh->X[iNode];
                double yloc = theMesh->Y[iNode];                
                double Rin = 0.4;
                double Rout = 2.0;  
                double value = 0.0;              
                double NormeCarree = xloc*xloc + yloc*yloc;
                if((Rout*Rout - NormeCarree) < (NormeCarree - Rin*Rin))
                {
                    value = vExt;
                }
                if(systIsY)
                {
                    v = -value*(xloc/sqrt(NormeCarree));
                }
                else
                {
                    v = value*(yloc/sqrt(NormeCarree));
                }     
                femBandSystemConstrain(theSystem, number[iNode],v);                  
            }
        }
    }
  
    double *soluce = femBandSystemEliminate(theSystem);
    if(systIsY){
        for (i = 0; i < theProblem->mesh->nNode; i++)        
        theProblem->soluceY[i] += soluce[number[i]];
    }
    else{
        for (i = 0; i < theProblem->mesh->nNode; i++)        
        theProblem->soluceX[i] += soluce[number[i]];
    }  
}

double *femCouetteNorme(femCouetteProblem *theProblem)
{
    double X, Y;
    double *soluce = malloc(theProblem->mesh->nNode*sizeof(double));
    for(int i = 0; i < theProblem->mesh->nNode; i++)
    {
        X = theProblem->soluceX[i];
        Y = theProblem->soluceY[i];
        soluce[i] = sqrt(X*X+Y*Y);
    }
    return soluce;
}

femGrains *femGrainsCreateSimple(int n, double r, double m, double radiusIn, double radiusOut, femMesh *theMesh, double gamma)
{
    int i,nContact = n*(n-1)/2;
    
    femGrains *theGrains = malloc(sizeof(femGrains));
    theGrains->n = n;
    theGrains->radiusIn = radiusIn;
    theGrains->radiusOut = radiusOut;
    theGrains->gravity[0] =  0.0;
    theGrains->gravity[1] = -9.81;
    theGrains->gamma = gamma;
    
       
    theGrains->x  = malloc(n*sizeof(double));
    theGrains->y  = malloc(n*sizeof(double));
    theGrains->vx = malloc(n*sizeof(double));
    theGrains->vy = malloc(n*sizeof(double));
    theGrains->r  = malloc(n*sizeof(double));
    theGrains->m  = malloc(n*sizeof(double));       
    theGrains->dvBoundary = malloc(n * sizeof(double));
    theGrains->dvContacts = malloc(nContact * sizeof(double));
   
    for(i = 0; i < n; i++) {
        theGrains->r[i] = r;
        theGrains->m[i] = m;
        theGrains->x[i] = (i%5) * r * 2.5 - 5 * r + 1e-8; 
        theGrains->y[i] = (i/5) * r * 2.5 + 2 * r + radiusIn;       
        theGrains->vx[i] = 0.0;
        theGrains->vy[i] = 0.0; 
        theGrains->dvBoundary[i] = 0.0; }
    for(i = 0; i < nContact; i++)  
        theGrains->dvContacts[i] = 0.0;

  
    return theGrains;
}

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

void femGrainsFree(femGrains *theGrains)
{
    free(theGrains->x);
    free(theGrains->y);
    free(theGrains->vx);
    free(theGrains->vy);
    free(theGrains->r);
    free(theGrains->m);
    free(theGrains->dvBoundary);
    free(theGrains->dvContacts);
    free(theGrains);
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
