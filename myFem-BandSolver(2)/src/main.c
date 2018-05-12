<<<<<<< HEAD
/*
*  main.c
*  Library for MECA1120 : Finite Elements for dummies
*
*  Copyright (C) 2018 UCL-EPL : Vincent Legat    
*  All rights reserved.
*
*/
=======
ï»¿/*
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2018 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 */
>>>>>>> 892f4706b2c091bd45fab6dc95d4a32aa05e66b1


#include "glfem.h"
#include <time.h>



int main(void)
{
    printf("\n\n    V : Plot results \n");
    printf("    S : Spy matrix \n");
    printf("    X-Y-N : Renumbering along x - along y - No renumbering \n");

    femSolverType solverType = FEM_ITER;
    femRenumType  renumType  = FEM_YNUM;
<<<<<<< HEAD
    char meshFileName[] = "../data/meshMedium.txt";  

// Pour Windows, remplacer l'argument :
// ("../data/rect_quad_1601.txt") 
// par :
// ("..\\data\\rect_quad_1601.txt") 
//
// Sorry for the inconvenience :-)
// On réfléchit pour rendre cela plus transparent dans les homeworks suivants :-)
// Be patient !
=======
    char meshFileName[] = "../data/rect_quad_1601.txt";

    // Pour Windows, remplacer l'argument :
    // ("../data/rect_quad_1601.txt")
    // par :
    // ("..\\data\\rect_quad_1601.txt")
    //
    // Sorry for the inconvenience :-)
    // On rÃ©flÃ©chit pour rendre cela plus transparent dans les homeworks suivants :-)
    // Be patient !
>>>>>>> 892f4706b2c091bd45fab6dc95d4a32aa05e66b1



    femDiffusionProblem* theProblem = femDiffusionCreate(meshFileName,solverType,renumType);
    clock_t tic = clock();
    femDiffusionCompute(theProblem);
    int testConvergence;
    do {
<<<<<<< HEAD
        femDiffusionCompute(theProblem);  
        femSolverPrintInfos(theProblem->solver); 
        testConvergence = femSolverConverged(theProblem->solver); 
=======
        femDiffusionCompute(theProblem);
        femSolverPrintInfos(theProblem->solver);
        testConvergence = femSolverConverged(theProblem->solver);
>>>>>>> 892f4706b2c091bd45fab6dc95d4a32aa05e66b1
    }
    while ( testConvergence == 0);
    printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
    printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
    fflush(stdout);


<<<<<<< HEAD
    int option = 1;     
=======
    int option = 1;
>>>>>>> 892f4706b2c091bd45fab6dc95d4a32aa05e66b1
    femRenumType  newRenumType  = renumType;

    GLFWwindow* window = glfemInit("MECA1120 : homework 5 ");
    glfwMakeContextCurrent(window);

<<<<<<< HEAD
    do 
=======
    do
>>>>>>> 892f4706b2c091bd45fab6dc95d4a32aa05e66b1
    {
        int w,h;
        char theMessage[256];
        sprintf(theMessage, "Max : %.4f ",femMax(theProblem->soluce,theProblem->size));
        glfwGetFramebufferSize(window,&w,&h);

        if (option == 1) {
            glfemReshapeWindows(theProblem->mesh,w,h);
<<<<<<< HEAD
            glfemPlotField(theProblem->mesh,theProblem->soluce);   
        }
        else {
            glColor3f(1.0,0.0,0.0);
            glfemPlotSolver(theProblem->solver,theProblem->size,w,h); 
        }
        glColor3f(0.0,0.0,0.0); glfemDrawMessage(20,460,theMessage);              

        if (renumType != newRenumType) {             
=======
            glfemPlotField(theProblem->mesh,theProblem->soluce);
        }
        else {
            glColor3f(1.0,0.0,0.0);
            glfemPlotSolver(theProblem->solver,theProblem->size,w,h);
        }
        glColor3f(0.0,0.0,0.0);
        glfemDrawMessage(20,460,theMessage);

        if (renumType != newRenumType) {
>>>>>>> 892f4706b2c091bd45fab6dc95d4a32aa05e66b1
            renumType = newRenumType;
            femDiffusionFree(theProblem);
            theProblem = femDiffusionCreate(meshFileName,solverType,renumType);
            clock_t tic = clock();
            do {
<<<<<<< HEAD
                femDiffusionCompute(theProblem);  
                femSolverPrintInfos(theProblem->solver); 
                testConvergence = femSolverConverged(theProblem->solver); 
=======
                femDiffusionCompute(theProblem);
                femSolverPrintInfos(theProblem->solver);
                testConvergence = femSolverConverged(theProblem->solver);
>>>>>>> 892f4706b2c091bd45fab6dc95d4a32aa05e66b1
            }
            while ( testConvergence == 0);
            if (testConvergence == -1)  printf("    Iterative solver stopped afer a maximum number of iterations\n");
            printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
            switch (renumType) {
<<<<<<< HEAD
                case FEM_XNUM : printf("    Renumbering along the x-direction\n"); break;
                case FEM_YNUM : printf("    Renumbering along the y-direction\n"); break;
                default : break; 
            }
            printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
            fflush(stdout); 
        }
        if (glfwGetKey(window,'V') == GLFW_PRESS)   option = 1;
        if (glfwGetKey(window,'S') == GLFW_PRESS)   option = 0; 
        if (glfwGetKey(window,'X') == GLFW_PRESS)   newRenumType  = FEM_XNUM; 
        if (glfwGetKey(window,'Y') == GLFW_PRESS)   newRenumType  = FEM_YNUM; 
        if (glfwGetKey(window,'N') == GLFW_PRESS)   newRenumType  = FEM_NO; 

        glfwSwapBuffers(window);
        glfwPollEvents();
    } 
    while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) != 1 );

// Check if the ESC key was pressed or the window was closed

    glfwTerminate(); 
=======
            case FEM_XNUM :
                printf("    Renumbering along the x-direction\n");
                break;
            case FEM_YNUM :
                printf("    Renumbering along the y-direction\n");
                break;
            default :
                break;
            }
            printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
            fflush(stdout);
        }
        if (glfwGetKey(window,'V') == GLFW_PRESS)   option = 1;
        if (glfwGetKey(window,'S') == GLFW_PRESS)   option = 0;
        if (glfwGetKey(window,'X') == GLFW_PRESS)   newRenumType  = FEM_XNUM;
        if (glfwGetKey(window,'Y') == GLFW_PRESS)   newRenumType  = FEM_YNUM;
        if (glfwGetKey(window,'N') == GLFW_PRESS)   newRenumType  = FEM_NO;

        glfwSwapBuffers(window);
        glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );

    // Check if the ESC key was pressed or the window was closed

    glfwTerminate();
>>>>>>> 892f4706b2c091bd45fab6dc95d4a32aa05e66b1
    femDiffusionFree(theProblem);
    exit(EXIT_SUCCESS);


<<<<<<< HEAD

}
=======
>>>>>>> 892f4706b2c091bd45fab6dc95d4a32aa05e66b1

}
