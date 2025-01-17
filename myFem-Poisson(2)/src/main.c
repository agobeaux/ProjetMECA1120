/*
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 */

#include "glfem.h"



int main(void)
{ 
    double mu = 1.0;
    double gamma = 1;
    double vExt = 4.0;

    //GRAINS
    int    n = 0;
    double radius    = 0.1;
    double mass      = 0.1;
    double radiusIn  = 0.4;
    double radiusOut = 2.0;    
    double dt      = 1e-1;
    double tEnd    = 8.0;
    double tol     = 1e-6;
    double t       = 0;
    double iterMax = 100;    
 
    femPoissonProblem* theProblem = femPoissonCreate("../data/meshMedium.txt");
    femGrains* theGrains = femGrainsCreateSimple(n,radius,mass,radiusIn,radiusOut, theProblem->mesh, gamma);
    
    // Pour Windows, remplacer l'argument :
    // ("../data/triangles_166.txt") 
    // par :
    // ("..\\data\\triangles_166.txt") 
    //
    // Sorry for the inconvenience :-)
    // On r�fl�chit pour rendre cela plus transparent dans les homeworks suivants :-)
    // Be patient !
    
     
    printf("Number of elements    : %4d\n", theProblem->mesh->nElem);
    printf("Number of local nodes : %4d\n", theProblem->mesh->nLocalNode);
    printf("Number of segments    : %4d\n", theProblem->edges->nBoundary);
    printf("Number of unknowns    : %4d\n", theProblem->systemY->size);

    femPoissonSolve(theProblem, theGrains, mu, gamma, vExt, 1);
    femPoissonSolve(theProblem, theGrains, mu, gamma, vExt, 0);

        
 
    printf("Maximum value : %.4f\n", femMax(theProblem->systemY->B,theProblem->systemY->size));
    fflush(stdout);
    
    char theMessage[256];
    sprintf(theMessage, "Max : %.4f", femMax(theProblem->systemY->B,theProblem->systemY->size));
    
    GLFWwindow* window = glfemInit("MECA1120 : FEM PROJECT ");
    glfwMakeContextCurrent(window);
    int theRunningMode = 1.0;
    float theVelocityFactor = 0.25;
    do {
        int i,w,h;
        double currentTime = glfwGetTime();

        glfwGetFramebufferSize(window,&w,&h);
        glColor3f(1.0,0.0,0.0);
        glfemPlotSolver(theProblem->systemX,theProblem->systemX->size,w,h); 
        /*if(theRunningMode){
            glfemReshapeWindows(theProblem->mesh,w,h);
            glfemPlotField(theProblem->mesh,theProblem->systemY->B); 

           
        glColor3f(1.0,0.0,0.0); glfemDrawMessage(20,460,theMessage);  
        for (i=0 ;i < theGrains->n; i++) {     
            glColor3f(1,0,0); 
            glfemDrawDisk(theGrains->x[i],theGrains->y[i],theGrains->r[i]); 
        }       
        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusOut);
        glColor3f(0,0,0); glfemDrawCircle(0,0,radiusIn);         
        sprintf(theMessage,"Time = %g sec",t);
        glColor3f(1,0,0); glfemDrawMessage(20,460,theMessage);             



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
            t += dt; }
                    }
        else{
            glColor3f(1.0,0.0,0.0);
            glfemPlotSolver(theProblem->systemX,theProblem->systemX->size,w,h); 
        }
         
        while ( glfwGetTime()-currentTime < theVelocityFactor ) {
          if (glfwGetKey(window,'R') == GLFW_PRESS) 
                theRunningMode = 1; 
          if (glfwGetKey(window,'S') == GLFW_PRESS) 
                theRunningMode = 0; }
        glfwSwapBuffers(window);
        glfwPollEvents();*/
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             !glfwWindowShouldClose(window));
           
    // Check if the ESC key was pressed or the window was closed
               
    glfwTerminate(); 
    femPoissonFree(theProblem);
    femGrainsFree(theGrains);
    exit(EXIT_SUCCESS);
    
    

}

