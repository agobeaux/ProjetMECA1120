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
       
    if (theSpace->n > 3) Error("Unexpected discrete space size !"); 
    
    double Xloc[3],Yloc[3],phi[3],dphidxsi[3],phiGrains[3], dphideta[3],dphidx[3],dphidy[3],Aloc[9],Bloc[3],xGrains, yGrains, xsiGrains, etaGrains, vGrains;
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

void femCouetteRenumber(femCouetteProblem *theProblem, femRenumType renumType)
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

int femCouetteComputeBand(femCouetteProblem *theProblem)
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

double femGrainsContactIterate(femGrains *myGrains, double dt, int iter)  
{
    int i,j,iContact;    
    int n = myGrains->n;
    
    
    double *x          = myGrains->x;
    double *y          = myGrains->y;
    double *m          = myGrains->m;
    double *r          = myGrains->r;
    double *vy         = myGrains->vy;
    double *vx         = myGrains->vx;
    double *dvBoundary = myGrains->dvBoundary;
    double *dvContacts = myGrains->dvContacts;
    double rIn         = myGrains->radiusIn;
    double rOut        = myGrains->radiusOut;
    
    double error = 0.0;
    double gap,rr, rx, ry, nx, ny, vn, dv, dvIn, dvOut;

    error = 0;
    iContact = 0;
    for(i = 0; i < n; i++) {
        for(j = i+1; j < n; j++) { 
            rx = (x[j]-x[i]);
            ry = (y[j]-y[i]);
            rr = sqrt(rx*rx+ry*ry);
            nx = rx/rr;
            ny = ry/rr;
            if (iter == 0) {
                  dv = dvContacts[iContact]; }
            else {
                  vn = (vx[i]-vx[j])*nx + (vy[i]-vy[j])*ny ;
                  gap = rr - (r[i]+r[j]);
                  dv = fmax(0.0, vn + dvContacts[iContact] - gap/dt);
                  dv = dv - dvContacts[iContact];                      
                  dvContacts[iContact] += dv; 
                  error = fmax(fabs(dv),error); }
            vx[i] -= dv * nx * m[j] / ( m[i] + m[j] );
            vy[i] -= dv * ny * m[j] / ( m[i] + m[j] );
            vx[j] += dv * nx * m[i] / ( m[i] + m[j] );
            vy[j] += dv * ny * m[i] / ( m[i] + m[j] );                  
            iContact++; }}

    for(i = 0; i < n; i++) {
        rr = sqrt(x[i]*x[i]+y[i]*y[i]);      
        nx = x[i]/rr;
        ny = y[i]/rr;
        if (iter == 0) {
            dv = dvBoundary[i]; }
        else {
            vn = vx[i]*nx + vy[i]*ny ;      
            gap = rOut - rr - r[i];
            dvOut = fmax(0.0, vn + dvBoundary[i] - gap/dt);
            gap = rr - rIn - r[i];
            dvIn  = fmax(0.0,-vn - dvBoundary[i] - gap/dt);
            dv = dvOut - dvIn - dvBoundary[i]; 
            dvBoundary[i] += dv; 
            error = fmax(fabs(dv),error); }
        vx[i] -= dv * nx;
        vy[i] -= dv * ny; }  
    return error;
}

double fluidSpeed(double xGrains, double yGrains, femCouetteProblem *theProblem, int systIsY)
{
    femMesh *theMesh = theProblem->mesh;    
    femBandSystem *theSystem;
    double *soluce;
    if(systIsY)
    {
        theSystem = theProblem->systemY;
        soluce = theProblem->soluceY;
        
    }
    else
    {
        theSystem = theProblem->systemX; 
        soluce = theProblem->soluceX;   
    }    
    femDiscrete *theSpace = theProblem->space;
    double x[3],y[3],phiGrains[3], xsiGrains, etaGrains, speed = 0.0;
    int i,iElem,map[3], iGrains, flag = 1;  
    for(iElem = 0; iElem < theMesh->nElem && flag; iElem++){
        femCouetteMeshLocal(theProblem,iElem,map,x,y);
        if(elemContains(xGrains,yGrains,theMesh,iElem)==1)
        {
            xsiGrains = -(x[0] * (y[2] - yGrains) + x[2] * (yGrains - y[0]) + xGrains * (y[0] - y[2]))/(x[0] * (y[1] - y[2]) + x[1] * (y[2] - y[0]) + x[2] * (y[0] - y[1]));
            etaGrains = (x[0] * (y[1] - yGrains) + x[1] * (yGrains - y[0]) + xGrains * (y[0] - y[1]))/(x[0] * (y[1] - y[2]) + x[1] * (y[2] - y[0]) + x[2] * (y[0] - y[1]));
            femDiscretePhi2(theSpace,xsiGrains,etaGrains,phiGrains);
            for (i = 0; i < theSpace->n; i++) 
            { 
                speed += soluce[theMesh->elem[iElem*3+i]]*phiGrains[i]; 
            }
            flag = 0;
        }        
    }
    return speed;
}

void femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax, femCouetteProblem *theProblem)
{
    int i;    
    int n = myGrains->n;
    
    double *x          = myGrains->x;
    double *y          = myGrains->y;
    double *m          = myGrains->m;
    double *vy         = myGrains->vy;
    double *vx         = myGrains->vx;
    double gamma       = myGrains->gamma;
    double gx          = myGrains->gravity[0];
    double gy          = myGrains->gravity[1];

// 
// -1- Calcul des nouvelles vitesses des grains sur base de la gravit� et de la trainee
//

    for(i = 0; i < n; i++) {
        double xGrains = x[i];
        double yGrains = y[i];
        double u = fluidSpeed(xGrains,yGrains,theProblem, 0);
        double v = fluidSpeed(xGrains,yGrains,theProblem, 1);
        double fx = - gamma * (vx[i]-u);
        double fy = m[i] * gy - gamma * (vy[i]-v);
        //printf("u = %f\n", u);
        //printf("v = %f\n", v);
        vx[i] += fx * dt / m[i];
        vy[i] += fy * dt / m[i];  }

//
// -2- Correction des vitesses pour tenir compte des contacts        
//       

    int iter = 0;
    double error;
           
    do {
        error = femGrainsContactIterate(myGrains,dt,iter);
        iter++; }
    while ((error > tol/dt && iter < iterMax) || iter == 1);
    printf("iterations = %4d : error = %14.7e \n",iter-1,error);
 
//  
// -3- Calcul des nouvelles positions sans penetrations de points entre eux
//

    for (i = 0; i < n; ++i) {
        x[i] += vx[i] * dt;
        y[i] += vy[i] * dt; }
}
