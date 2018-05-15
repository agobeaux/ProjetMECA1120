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