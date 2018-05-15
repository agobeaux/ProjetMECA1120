/*
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-EPL : Vincent Legat    
 *  All rights reserved.
 *
 */


#include "glfem.h"
#include <time.h>



int main(void)
{
	double mu = 1.0;
    double gamma = 1;
    double vExt = 4.0;

    //GRAINS
    int    n = 10;
    double radius    = 0.1;
    double mass      = 0.1;
    double radiusIn  = 0.4;
    double radiusOut = 2.0;    
    double dt      = 1e-1;
    double tEnd    = 8.0;
    double tol     = 1e-6;
    double t       = 0;
    double iterMax = 100; 
    
    femSolverType solverType = FEM_BAND;
    femRenumType  renumType  = FEM_YNUM;
    char meshFileName[] = "../data/meshMedium.txt";  
    
    // Pour Windows, remplacer l'argument :
    // ("../data/rect_quad_1601.txt") 
    // par :
    // ("..\\data\\rect_quad_1601.txt") 
    //
    // Sorry for the inconvenience :-)
    // On réfléchit pour rendre cela plus transparent dans les homeworks suivants :-)
    // Be patient !

    
    
    femDiffusionProblem* theProblem = femDiffusionCreate(meshFileName,solverType,renumType);
    clock_t tic = clock();
    femDiffusionCompute(theProblem);  
    femSolverPrintInfos(theProblem->solver);
    printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
    printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
    fflush(stdout);
    

    int option = 1;    
    femSolverType newSolverType = solverType;
    femRenumType  newRenumType  = renumType;

    GLFWwindow* window = glfemInit("MECA1120 : FEM PROJECT ");
    glfwMakeContextCurrent(window);

    do 
    {
        int w,h;
        char theMessage[256];
        sprintf(theMessage, "Max : %.4f ",femMax(theProblem->soluce,theProblem->size));
        glfwGetFramebufferSize(window,&w,&h);        
        glfemReshapeWindows(theProblem->mesh,w,h);
        glfemPlotField(theProblem->mesh,theProblem->soluce);   

        glColor3f(0.0,0.0,0.0); glfemDrawMessage(20,460,theMessage);        
        for (i=0 ;i < theGrains->n; i++) {     
            glColor3f(1,0,0); 
            glfemDrawDisk(theGrains->x[i],theGrains->y[i],theGrains->r[i]); 
        }         
        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusOut);
        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusIn);         
        sprintf(theMessage,"Time = %g sec",t);
        glColor3f(1,0,0); glfemDrawMessage(20,460,theMessage);  

        
        glfwSwapBuffers(window);
        glfwPollEvents();
        
        if (t < tEnd && theRunningMode == 1) {
            printf("Time = %4g : ",t);
  //
  // A decommenter pour pouvoir progresser pas par pas
  //          printf("press CR to compute the next time step >>");
  //          char c= getchar();
  //
            femGrainsUpdate(theGrains,dt,tol,iterMax, theProblem);
            femFullSystemInit2(theProblem->systemX, theProblem->systemY);
            femPoissonSolve(theProblem, theGrains, mu, gamma, vExt, 1);
            femPoissonSolve(theProblem, theGrains, mu, gamma, vExt, 0);
            t += dt;
            
            femDiffusionFree(theProblem);
            theProblem = femDiffusionCreate(meshFileName,solverType,renumType);
            theProblem = femDiffusionCreate(meshFileName,solverType,renumType);
            do {
				femDiffusionCompute(theProblem);
				femSolverPrintInfos(theProblem->solver);
				theProblem = femDiffusionCreate(meshFileName,solverType,renumType);
		}
		while ( glfwGetTime()-currentTime < theVelocityFactor ) {
          if (glfwGetKey(window,'R') == GLFW_PRESS) 
                theRunningMode = 1; 
          if (glfwGetKey(window,'S') == GLFW_PRESS) 
                theRunningMode = 0; 
		}
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed
               
    glfwTerminate(); 
    femDiffusionFree(theProblem);
    femGrainsFree(theGrains);
    exit(EXIT_SUCCESS);
    
    

}

