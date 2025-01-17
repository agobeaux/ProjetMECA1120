

#include "fem.h"

# ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename)
{
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->mesh  = femMeshRead(filename);   
    femMeshClean(theProblem->mesh);        
    theProblem->edges = femEdgesCreate(theProblem->mesh);
    //femNeighbours(theProblem->mesh, theProblem->edges); // rajouté ici
    if (theProblem->mesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4,FEM_QUAD); }
    else if (theProblem->mesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE); }
    theProblem->systemX = femFullSystemCreate(theProblem->mesh->nNode);
    theProblem->systemY = femFullSystemCreate(theProblem->mesh->nNode);
    return theProblem;
}



# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{
    femFullSystemFree(theProblem->systemX);
    femFullSystemFree(theProblem->systemY);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    femEdgesFree(theProblem->edges);
    femMeshFree(theProblem->mesh);
    free(theProblem);
}
    
# endif
# ifndef NOMESHLOCAL

void femMeshLocal(const femMesh *theMesh, const int iElem, int *map, double *x, double *y)
{
    int j,nLocal = theMesh->nLocalNode;
    
    for (j=0; j < nLocal; ++j) {
        map[j] = theMesh->elem[iElem*nLocal+j];
        x[j]   = theMesh->X[map[j]];
        y[j]   = theMesh->Y[map[j]]; }   
}

# endif
# ifndef NOPOISSONSOLVE

void femPoissonSolve(femPoissonProblem *theProblem, femGrains *theGrains, double mu, double gamma, double vExt, int systIsY)
{
    femMesh *theMesh = theProblem->mesh;
    femEdges *theEdges = theProblem->edges;
    femFullSystem *theSystem;
    if(systIsY)
    {
        theSystem = theProblem->systemY;
    }
    else
    {
        theSystem = theProblem->systemX;    
    }
    
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;

    if (theSpace->n > 4) Error("Unexpected discrete space size !");  
    double x[4],y[4],phi[4],phiGrains[3],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4], xGrains, yGrains, xsiGrains, etaGrains, vGrains;
    int iElem,iInteg,iEdge,i,j,map[4], iGrains;    

    for (iElem = 0; iElem < theMesh->nElem; iElem++) 
    {
        femMeshLocal(theMesh,iElem,map,x,y);         
        for (iInteg=0; iInteg < theRule->n; iInteg++) 
        {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0; 
            double dydeta = 0;
            double xloc = 0;
            double yloc = 0;
            for (i = 0; i < theSpace->n; i++) 
            {  
                xloc   += x[i]*phi[i];  
                yloc   += y[i]*phi[i];  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; 
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < theSpace->n; i++) 
            {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; 
            }            
            for (i = 0; i < theSpace->n; i++) 
            { 
                for(j = 0; j < theSpace->n; j++) 
                {
                    theSystem->A[map[i]][map[j]] += mu*(dphidx[i] * dphidx[j] 
                     + dphidy[i] * dphidy[j]) * jac * weight; 
                }
            }
        }
        for(iGrains = 0; iGrains < theGrains->n; iGrains++)
        {            
            xGrains = theGrains->x[iGrains];
            yGrains = theGrains->y[iGrains]; 
            //if(elemContains(xGrains,yGrains,theMesh,iElem)==1)
            //{
                xsiGrains = -(x[0] * (y[2] - yGrains) + x[2] * (yGrains - y[0]) + xGrains * (y[0] - y[2]))/(x[0] * (y[1] - y[2]) + x[1] * (y[2] - y[0]) + x[2] * (y[0] - y[1]));
                etaGrains = (x[0] * (y[1] - yGrains) + x[1] * (yGrains - y[0]) + xGrains * (y[0] - y[1]))/(x[0] * (y[1] - y[2]) + x[1] * (y[2] - y[0]) + x[2] * (y[0] - y[1]));
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
                    theSystem->B[map[i]] += gamma*phiGrains[i]*vGrains; 
                    for(j = 0; j < theSpace->n; j++) 
                    {
                        theSystem->A[map[j]][map[i]] += gamma*phiGrains[i]*phiGrains[j];
                    }
                }
            //}
        }
    }    
    for (iEdge= 0; iEdge < theEdges->nEdge; iEdge++) 
    {      
        if (theEdges->edges[iEdge].elem[1] < 0) 
        {  
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
                
                femFullSystemConstrain(theSystem,iNode,v);  
            }
        }
    }

    femFullSystemEliminate(theSystem);
}


# endif

# ifndef NOCONTACTITERATE

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

# endif
# ifndef NOUPDATE

double fluidSpeed(double xGrains, double yGrains, femPoissonProblem *theProblem, int systIsY)
{
    femMesh *theMesh = theProblem->mesh;    
    femFullSystem *theSystem;
    if(systIsY)
    {
        theSystem = theProblem->systemY;
    }
    else
    {
        theSystem = theProblem->systemX;    
    }    
    femDiscrete *theSpace = theProblem->space;
    double x[3],y[3],phiGrains[3], xsiGrains, etaGrains, speed = 0.0;
    int i,iElem,map[3], iGrains, flag = 1;  
    for(iElem = 0; iElem < theMesh->nElem && flag; iElem++){
        femMeshLocal(theMesh,iElem,map,x,y);
        if(elemContains(xGrains,yGrains,theMesh,iElem)==1)
        {
            xsiGrains = -(x[0] * (y[2] - yGrains) + x[2] * (yGrains - y[0]) + xGrains * (y[0] - y[2]))/(x[0] * (y[1] - y[2]) + x[1] * (y[2] - y[0]) + x[2] * (y[0] - y[1]));
            etaGrains = (x[0] * (y[1] - yGrains) + x[1] * (yGrains - y[0]) + xGrains * (y[0] - y[1]))/(x[0] * (y[1] - y[2]) + x[1] * (y[2] - y[0]) + x[2] * (y[0] - y[1]));
            femDiscretePhi2(theSpace,xsiGrains,etaGrains,phiGrains);
            for (i = 0; i < theSpace->n; i++) 
            { 
                speed += theSystem->B[map[i]]*phiGrains[i]; 
            }
            flag = 0;
        }        
    }
    return speed;
}

void femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax, femPoissonProblem *theProblem)
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

# endif

