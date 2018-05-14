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
             
        glfwSwapBuffers(window);
        glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed
               
    glfwTerminate(); 
    femDiffusionFree(theProblem);
    exit(EXIT_SUCCESS);
    
    

}

