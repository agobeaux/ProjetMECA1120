
/*
 *  fem.h
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2017 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#ifndef _FEM_H_
#define _FEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define ErrorScan(a)   femErrorScan(a,__LINE__,__FILE__)
#define Error(a)       femError(a,__LINE__,__FILE__)
#define Warning(a)     femWarning(a,  __LINE__, __FILE__)
#define FALSE 0 
#define TRUE  1



typedef enum {FEM_TRIANGLE,FEM_QUAD} femElementType;
typedef enum {FEM_NO,FEM_XNUM,FEM_YNUM} femRenumType;

typedef struct {
    int *elem;
    double *X;
    double *Y;
    int nElem;
    int nNode;
    int nLocalNode;
} femMesh;

typedef struct {
    int elem[2];
    int node[2];
} femEdge;

typedef struct {
    femMesh *mesh;
    femEdge *edges;
    int nEdge;
    int nBoundary;
} femEdges;

typedef struct {
    int n;
    void (*x2)(double *xsi, double *eta);
    void (*phi2)(double xsi, double eta, double *phi);
    void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
} femDiscrete;
    
typedef struct {
    int n;
    const double *xsi;
    const double *eta;
    const double *weight;
} femIntegration;


typedef struct 
{
    double *B;
    double **A;        
    int size;
    int band;        
} femBandSystem;

typedef struct {
    femMesh *mesh;
    femEdges *edges;
    femDiscrete *space;
    femIntegration *rule;
    femBandSystem *systemX;
    femBandSystem *systemY;
    int size;
    int *number;
    double *soluceX;
    double *soluceY;
} femCouetteProblem;

typedef struct {
    int n;
    double radiusIn;
    double radiusOut;
    double gravity[2];
    double gamma;
    double *x;
    double *y;
    double *vx;
    double *vy;
    double *r;
    double *m;
    double *dvBoundary;
    double *dvContacts;
} femGrains;

femGrains  *femGrainsCreateSimple(int n, double r, double m, double radiusIn, double radiusOut, femMesh *theMesh, double gamma);
void        femGrainsFree(femGrains *myGrains);
double      fluidSpeed(double xGrains, double yGrains, femCouetteProblem *theProblem, int systIsY);
void        femGrainsUpdate(femGrains *myGrains, double dt, double tol, double iterMax, femCouetteProblem *theProblem);
double      femGrainsContactIterate(femGrains *myGrains, double dt, int iter);
int         elemContains(double x, double y, femMesh *theMesh, int iElem);

femIntegration      *femIntegrationCreate(int n, femElementType type);
void                 femIntegrationFree(femIntegration *theRule);

femMesh             *femMeshRead(const char *filename);
void                 femMeshWrite(const femMesh* myMesh, const char *filename);
void                 femMeshFree(femMesh *theMesh);

femEdges*            femEdgesCreate(femMesh *theMesh);
void                 femEdgesFree(femEdges *theEdges);
void                 femEdgesPrint(femEdges *theEdges);
int                  femEdgesCompare(const void *edgeOne, const void *edgeTwo);

femDiscrete*         femDiscreteCreate(int n, femElementType type);
void                 femDiscreteFree(femDiscrete* mySpace);
void                 femDiscretePrint(femDiscrete* mySpace);
void                 femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta);
void                 femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi);
void                 femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);

femBandSystem*       femBandSystemCreate(int size, int band);
void                 femBandSystemFree(femBandSystem* myBandSystem);
void                 femBandSystemInit(femBandSystem *myBand);
void                 femBandSystemPrintInfos(femBandSystem *myBand);
double*              femBandSystemEliminate(femBandSystem *myBand);
void                 femBandSystemConstrain(femBandSystem *myBand, int myNode, double myValue);
void                 femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc);
 
femCouetteProblem*   femCouetteCreate(const char *filename, femRenumType renumType);
void                 femSoluceInit(femCouetteProblem *theProblem);
void                 femCouetteFree(femCouetteProblem *theProblem);
void                 femCouetteMeshLocal(const femCouetteProblem *theProblem, const int i, int *map, double *x, double *y);
void                 femCouetteCompute(femCouetteProblem *theProblem, femGrains *theGrains, double mu, double gamma, double vExt, int systIsY);
double*              femCouetteNorme(femCouetteProblem *theProblem);
void                 femCouetteRenumber(femCouetteProblem *theProblem, femRenumType renumType);
int                  femCouetteComputeBand(femCouetteProblem *theProblem);

int                  compare(const void *nodeOne, const void *nodeTwo) ;
double               femMin(double *x, int n);
double               femMax(double *x, int n);
void                 femError(char *text, int line, char *file);
void                 femErrorScan(int test, int line, char *file);
void                 femWarning(char *text, int line, char *file);

#endif