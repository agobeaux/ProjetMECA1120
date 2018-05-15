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
    int    n = 100;
    double radius    = 0.1;
    double mass      = 0.1;
    double radiusIn  = 0.4;
    double radiusOut = 2.0;    
    double dt      = 1e-1;
    double tEnd    = 8.0;
    double tol     = 1e-6;
    double t       = 0;
    double iterMax = 100; 
    
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

    
    
    femCouetteProblem* theProblem = femCouetteCreate(meshFileName,renumType);
    femGrains* theGrains = femGrainsCreateSimple(n,radius,mass,radiusIn,radiusOut, theProblem->mesh, gamma);
    clock_t tic = clock();
    femCouetteCompute(theProblem, theGrains, mu, gamma, vExt, 1);
    femCouetteCompute(theProblem, theGrains, mu, gamma, vExt, 0);
    femBandSystemPrintInfos(theProblem->systemX);
    printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
    printf("    Maximum value : %.4f\n", femMax(femCouetteNorme(theProblem),theProblem->size));
    fflush(stdout);
    
    GLFWwindow* window = glfemInit("MECA1120 : FEM PROJECT ");
    glfwMakeContextCurrent(window);

    int theRunningMode = 1.0;
    float theVelocityFactor = 0.25;

    do 
    {
        int i,w,h;
        char theMessage[256];
        double currentTime = glfwGetTime();        
        glfwGetFramebufferSize(window,&w,&h);        
        glfemReshapeWindows(theProblem->mesh,w,h);
        glfemPlotField(theProblem->mesh,femCouetteNorme(theProblem));   
       
        for (i=0 ;i < theGrains->n; i++) {     
            
            glColor3f(255,255,0); 
            glfemDrawDisk(theGrains->x[i],theGrains->y[i],theGrains->r[i]); 
            double xNodes[2] = {theGrains->x[i]+theGrains->r[i]/2, theGrains->x[i]-(2/3)*theGrains->r[i]};
            double yNodes[2] = {theGrains->y[i]+theGrains->r[i]/2, theGrains->y[i]+(2/3)*theGrains->r[i]};
            glColor3f(1,0,0);
            glfemDrawNodes(xNodes, yNodes, 2, theGrains->r[i]);
            
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
            t += dt;  
            femSoluceInit(theProblem); 
            femBandSystemInit(theProblem->systemX);
            femBandSystemInit(theProblem->systemY);
			femCouetteCompute(theProblem, theGrains, mu, gamma, vExt,1);
			femCouetteCompute(theProblem, theGrains, mu, gamma, vExt,0);
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
    femCouetteFree(theProblem);
    femGrainsFree(theGrains);
    exit(EXIT_SUCCESS);
    
    

}

