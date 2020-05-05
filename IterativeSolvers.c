#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <string.h>



/* function used to print the data to a text file */
int print_current_data(int step, double** laplace_phi, double** f, double** phi, double* error, int NX, int NY, char method[]) {
    int i, j;


    char filename1[40] = "step_";
    itoa(step, filename1 + 4, 10);
    strcat(filename1, "_");
    strcat(filename1, method);
    strcat(filename1, "_laplace_phi_data.txt");

    FILE* fpointer1 = fopen(filename1, "w");
    for (j = 0; j < NY; j++) {
        for (i = 0; i < NX; i++) {
            fprintf(fpointer1, "%.20lf ", laplace_phi[j][i]);
        }
        fprintf(fpointer1, "\n");
    }
    fclose(fpointer1);



    char filename2[25] = "step_";
    itoa(step, filename2 + 4, 10);
    strcat(filename2, "_");
    strcat(filename2, method);
    strcat(filename2, "_f_data.txt");

    FILE* fpointer2 = fopen(filename2, "w");
    for (j = 0; j < NY; j++) {
        for (i = 0; i < NX; i++) {
            fprintf(fpointer2, "%.20lf ", f[j][i]);
        }
        fprintf(fpointer2, "\n");
    }
    fclose(fpointer2);



    char filename3[30] = "step_";
    itoa(step, filename3 + 4, 10);
    strcat(filename3, "_");
    strcat(filename3, method);
    strcat(filename3, "_phi_data.txt");

    FILE* fpointer3 = fopen(filename3, "w");
    for (j = 0; j < NY; j++) {
        for (i = 0; i < NX; i++) {
            fprintf(fpointer3, "%.20lf ", phi[j][i]);
        }
        fprintf(fpointer3, "\n");
    }
    fclose(fpointer3);



    char filename4[30] = "step_";
    itoa(step, filename4 + 4, 10);
    strcat(filename4, "_");
    strcat(filename4, method);
    strcat(filename4, "_error_data.txt");

    FILE* fpointer4 = fopen(filename4, "w");
    for (j = 0; j < step; j++) {
        fprintf(fpointer4, "%.20lf ", error[j]);
    }
    fclose(fpointer4);


    return 1;
}



/* ==========================================================================
======================= 3 Step Gauss-Seidel Solver ==========================
========================================================================== */

int GS_3step(double** f, double** phi, int Nx, int Ny, double epsilon) {

    /* Initilaizations */
    int i, j;
    double nx = Nx;
    double ny = Ny;
    double D_x = 1 / nx;
    double D_y = 1 / ny;

    double  f_norm;
    double integral;
    double lambda = pow(D_x, -2);
    double RHS;    
    double laplace_phi_minus_f_norm;
    
    int step = 1;
    int max_num_steps = 3;

    double** laplace_phi = (double**)calloc(Ny, sizeof(double*));
    for (j = 0; j < Ny; j++) {
        laplace_phi[j] = (double*)calloc(Nx, sizeof(double));
    }


    /* Solving -------------------------------------------------------------------------------- */

    f_norm = 0; /* compute f_norm */
    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            f_norm = f_norm + pow(f[j][i], 2);
        }
    }

    f_norm = sqrt(f_norm);


    if (Nx > 1 && Ny > 1) {

        do {

            if (Nx > 2 && Ny > 2) {

                for (j = 1; j < Ny - 1; j++) {
                    for (i = 1; i < Nx - 1; i++) {
                        phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 4 - f[j][i] / (4 * lambda);

                    }
                }
                

                //update left boundary values
                for (j = 1; j < Ny - 1; j++) {
                    i = 0;
                    phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

                }

                //update right boundary values
                for (j = 1; j < Ny - 1; j++) {
                    i = Nx - 1;
                    phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

                }

                //update bottom boundary values 
                for (i = 1; i < Nx - 1; i++) {
                    j = 0;
                    phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
                }

                //update top boundary values
                for (i = 1; i < Nx - 1; i++) {
                    j = Ny - 1;
                    phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) / 3 - f[j][i] / (3 * lambda);
                }

            }


            //update corner points
            phi[0][0] = (phi[1][0] + phi[0][1]) / 2 - f[0][0] / (2 * lambda);
            phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) / 2 - f[0][Nx - 1] / (2 * lambda);
            phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) / 2 - f[Ny - 1][0] / (2 * lambda);
            phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) / 2 - f[Ny - 1][Nx - 1] / (2 * lambda);

            
            if (Nx > 2 && Ny > 2) {

                //compute laplace_p matrix
                for (j = 1; j < Ny - 1; j++) {
                    for (i = 1; i < Nx - 1; i++) {
                        laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 4 * lambda * phi[j][i];
                    }
                }
                for (j = 1; j < Ny - 1; j++) {
                    i = 0;
                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
                }
                for (j = 1; j < Ny - 1; j++) {
                    i = Nx - 1;
                    laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
                }
                for (i = 1; i < Nx - 1; i++) {
                    j = 0;
                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
                }
                for (i = 1; i < Nx - 1; i++) {
                    j = Ny - 1;
                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) * lambda - phi[j][i] * (3 * lambda);
                }

            }

            laplace_phi[0][0] = (phi[1][0] + phi[0][1]) * lambda - phi[0][0] * (2 * lambda);
            laplace_phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) * lambda - phi[0][Nx - 1] * (2 * lambda);
            laplace_phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) * lambda - phi[Ny - 1][0] * (2 * lambda);
            laplace_phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) * lambda - phi[Ny - 1][Nx - 1] * (2 * lambda);



            //compute the norm
            laplace_phi_minus_f_norm = 0;

            for (j = 0; j < Ny; j++) {
                for (i = 0; i < Nx; i++) {
                    laplace_phi_minus_f_norm = laplace_phi_minus_f_norm + pow((laplace_phi[j][i] - f[j][i]), 2);
                }
            }

            

            laplace_phi_minus_f_norm = sqrt(laplace_phi_minus_f_norm);


            if (f_norm == 0) {
                RHS = epsilon;
            }
            else {
                RHS = epsilon * f_norm;
            }


            /* break out of the loop if the max number of steps has been reached */
            printf("\tstep %d\n", step);

            if (step == max_num_steps) {
                break;
            }


           step += 1;

        } while (laplace_phi_minus_f_norm > RHS);
    }


    /* Impose condition that the integral over the domain is equal to zero */
    integral = 0;
    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            integral += phi[j][i] * D_x * D_y;
        }
    }

    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            phi[j][i] = phi[j][i] - integral / (Nx * Ny);
        }
    }


    for (i = 0; i < nx; i++) {
        free(laplace_phi[i]);
    }

    free(laplace_phi);

    printf("\tGS done\n");

    return 0;

}

/* ====================================================================================================================
=================================== MG_recusion =======================================================================
==================================================================================================================== */


int MG_recursion(double** f, double** phi, int Nx, int Ny, double epsilon) {

    printf("Start MG_recursion: %dx%d\n", Nx, Ny);

        /* Initializing variables */
        int finished;
        int i, j;

        double nx = Nx;
        double ny = Ny;

        double half_nx = ceil(nx / 2);
        double half_ny = ceil(ny / 2);
        int half_Nx = half_nx;
        int half_Ny = half_ny;

        double** phi2 = (double**)calloc(half_ny, sizeof(double*));
        double** f2 = (double**)calloc(half_ny, sizeof(double*));

        for (i = 0; i < half_Nx; i++) {
            phi2[i] = (double*)calloc(half_nx, sizeof(double));
            f2[i] = (double*)calloc(half_nx, sizeof(double));

        }
        
        /* ============================================ Recursive code ============================================*/

        /* Initial GS solving */
        printf("1st GS_3step\n");
        finished = GS_3step(f, phi, Nx, Ny, epsilon);
        printf("done\n");

        /* Shrink phi & f by taking average of a point and its surounding 4 points and placing it into phi2 & f2 (adjusted for edges & corners) */
        if (Nx > 2 && Ny > 2) { /* For phi & f sizes greater than 2x2 */

            /* For interior points */
            for (j = 1; j < (half_Ny - 1); j++){
                for(i = 1; i < (half_Nx - 1);  i++) {

                    phi2[j][i] = (phi[j * 2][i * 2] + phi[j * 2][i * 2 + 1] + phi[j * 2 + 1][i * 2] + phi[j * 2][i * 2 - 1] + phi[j * 2 - 1][i * 2]) / 5;
                    f2[j][i] = (f[j * 2][i * 2] + f[j * 2][i * 2 + 1] + f[j * 2 + 1][i * 2] + f[j * 2][i * 2 - 1] + f[j * 2 - 1][i * 2]) / 5;
                
                }
            }

            /* For the left boundary */
            for (j = 1; j < (half_Ny - 1); j++) {
                i = 0;
                phi2[j][i] = (phi[j * 2][i * 2] + phi[j * 2][i * 2 + 1] + phi[j * 2 + 1][i * 2] + phi[j * 2 - 1][i * 2]) / 4;
                f2[j][i] = (f[j * 2][i * 2] + f[j * 2][i * 2 + 1] + f[j * 2 + 1][i * 2] + f[j * 2 - 1][i * 2]) / 4;

            }

            /* For the bottom boundary */
            for (i = 1; j < (half_Nx - 1); j++) {
                j = 0;
                phi2[j][i] = (phi[j * 2][i * 2] + phi[j * 2 + 1][i * 2] + phi[j * 2][i * 2 + 1] + phi[j * 2][i * 2 - 1]) / 4;
                f2[j][i] = (f[j * 2][i * 2] + f[j * 2 + 1][i * 2] + f[j * 2][i * 2 + 1] + f[j * 2][i * 2 - 1]) / 4;

            }

            /* ODD */
            if (nx / half_nx != 2 && ny / half_ny != 2) { /* For odd-sized square meshes (7x7, 23x23, etc.) */
                /* For the right boundary */
                for (j = 1; j < (half_Ny - 1); j++) {
                    i = half_Nx - 1;
                    phi2[j][i] = (phi[j * 2][i * 2] + phi[j * 2][i * 2 + 1] + phi[j * 2 + 1][i * 2] + phi[j * 2 - 1][i * 2]) / 4;
                    f2[j][i] = (f[j * 2][i * 2] + f[j * 2][i * 2 + 1] + f[j * 2 + 1][i * 2] + f[j * 2 - 1][i * 2]) / 4;

                }

                /* For the top boundary */
                for (i = 1; j < (half_Nx - 1); j++) {
                    j = half_Ny - 1;
                    phi2[j][i] = (phi[j * 2][i * 2] + phi[j * 2 + 1][i * 2] + phi[j * 2][i * 2 + 1] + phi[j * 2][i * 2 - 1]) / 4;
                    f2[j][i] = (f[j * 2][i * 2] + f[j * 2 + 1][i * 2] + f[j * 2][i * 2 + 1] + f[j * 2][i * 2 - 1]) / 4;

                }

                /* For the corners */
                phi2[0][0] = (phi[0][0] + phi[0][1] + phi[1][0]) / 3; /* bottom left */
                phi2[0][half_Nx - 1] = (phi[0][(half_Nx - 1) * 2] + phi[0][((half_Nx - 1) * 2) - 1] + phi[1][(half_Nx - 1) * 2]) / 3; /* bottom right */
                phi2[half_Ny - 1][0] = (phi[(half_Ny - 1) * 2][0] + phi[((half_Ny - 1) * 2) - 1][0] + phi[(half_Ny - 1) * 2][1]) / 3; /* top left */
                phi2[half_Ny - 1][half_Nx - 1] = (phi[(half_Ny - 1) * 2][(half_Nx - 1) * 2] +  + phi[(half_Ny - 1) * 2][(half_Nx - 1) * 2 - 1] + phi[(half_Ny - 1) * 2 - 1][(half_Nx - 1) * 2]) / 3; /* top right */
                
            }
            /* EVEN */
            else if (nx / half_nx == 2 && ny / half_ny == 2){ /* For even-sized square meshes (6x6, 24x24, etc.) */
                 /* For the right boundary */
                for (j = 1; j < (half_Ny - 1); j++) {
                    i = half_Nx - 1;
                    phi2[j][i] = (phi[j * 2][i * 2] + phi[j * 2][i * 2 + 1] + phi[j * 2 + 1][i * 2] + phi[j * 2][i * 2 - 1] + phi[j * 2 - 1][i * 2]) / 5;
                    f2[j][i] = (f[j * 2][i * 2] + f[j * 2][i * 2 + 1] + f[j * 2 + 1][i * 2] + f[j * 2][i * 2 - 1] + f[j * 2 - 1][i * 2]) / 5;

                }

                /* For the top boundary */
                for (i = 1; j < (half_Nx - 1); j++) {
                    j = half_Ny - 1;
                    phi2[j][i] = (phi[j * 2][i * 2] + phi[j * 2][i * 2 + 1] + phi[j * 2 + 1][i * 2] + phi[j * 2][i * 2 - 1] + phi[j * 2 - 1][i * 2]) / 5;
                    f2[j][i] = (f[j * 2][i * 2] + f[j * 2][i * 2 + 1] + f[j * 2 + 1][i * 2] + f[j * 2][i * 2 - 1] + f[j * 2 - 1][i * 2]) / 5;
                
                }

                /* For the corners */
                phi2[0][0] = (phi[0][0] + phi[0][1] + phi[1][0]) / 3; /* bottom left */
                phi2[0][half_Nx - 1] = (phi[0][(half_Nx - 1) * 2] + phi[0][((half_Nx - 1) * 2) + 1] + phi[0][((half_Nx - 1) * 2) - 1] + phi[1][(half_Nx - 1) * 2]) / 4; /* bottom right */
                phi2[half_Ny - 1][0] = (phi[(half_Ny - 1) * 2][0] + phi[((half_Ny - 1) * 2) + 1][0] + phi[((half_Ny - 1) * 2) - 1][0] + phi[(half_Ny - 1) * 2][1]) / 4; /* top left */
                phi2[half_Ny - 1][half_Nx - 1] = (phi[(half_Ny - 1) * 2][(half_Nx - 1) * 2] + phi[(half_Ny - 1) * 2][(half_Nx - 1) * 2 + 1] + phi[(half_Ny - 1) * 2 + 1][(half_Nx - 1) * 2] + phi[(half_Ny - 1) * 2][(half_Nx - 1) * 2 - 1] + phi[(half_Ny - 1) * 2 - 1][(half_Nx - 1) * 2]) / 5; /* top right */

            }
            else { /* Making sure the mesh is square, and terminating recursion otherwise (also frees memory) */
                printf("\n\n YOU DID NOT INPUT A SQUARE MESH\n\n");

                for (j = 0; j < half_Ny; j++) {
                    free(phi2[j]);
                    free(f2[j]);
                }
                free(phi2);
                free(f2);                

                return 0;
            }

        }
        else if (Nx == 2 && Ny == 2){ /* For phi & f sizes of 2x2 */
            phi2[0][0] = (phi[0][0] + phi[0][1] + phi[1][0]) / 3;
        }

        /* phi & f sizes of 1x1 are not shrunk */

        
        if (Nx > 1 && Ny > 1) {
            finished = MG_recursion(f2, phi2, half_Nx, half_Ny, epsilon);
        }
        printf("Exit MG_recursion: %dx%d\n", Nx, Ny);

        /* Merge phi */
        for (j = 0; j < Ny; j += 2){
            for(i = 0; i < Nx; i += 2) {

                phi[j][i] = phi2[j / 2][i / 2];

            }
        }

        /* Exit GS solving */
        printf("2nd GS_3step\n");
        finished = GS_3step(f, phi, Nx, Ny, epsilon);
        printf("done\n");


        for (j = 0; j < half_Ny; j++) {
            free(phi2[j]);
            free(f2[j]);
        }
        free(phi2);
        free(f2);

    return 0;
}




/* ================================================================================
= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
======================================= MAIN ======================================
= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
================================================================================ */


int main() {

    /* ----------------------------------------------------------------------------------------------------------
    Initializing variables -------------------------------------------------------------------------------- */
    
    char method[] = "SOR";      /* possible methods "PJ", "GS", "SOR", "CG", "MG" */

    int i, j;
    int step = 1;
    int max_num_steps = 2000;   /* Increase this number if the method isn't converging */

    int Nx = 12;
    int Ny = 12;
    double nx = Nx;
    double ny = Ny;
    double D_x = 1 / nx;
    double D_y = 1 / ny;
    double x = D_x / 2;
    double y = D_y / 2;

    double f_norm;
    double lambda = pow(D_x, -2);
    double laplace_phi_minus_f_norm;
    double epsilon = pow(10, -3);
    double RHS;
    double integral;
    double sum;

    double omega = 72;   /* SOR method variables */
    double phi_GS;
    double d_GS;

    int print_now;
    int finished;


    // The following variables are used in conjugate gradient method
    double** A = (double**)calloc(ny*nx, sizeof(double*));
    for (i = 0; i <nx*ny; i++) {
        A[i] = (double*)calloc(nx*ny, sizeof(double));
    }

    double* phi_vec = (double*)calloc(nx*ny, sizeof(double));
    double* f_vec = (double*)calloc(nx*ny, sizeof(double));
    double* r = (double*)calloc(nx*ny, sizeof(double));
    double* d = (double*)calloc(nx*ny, sizeof(double));
    double* epsilon_vec = (double*)calloc(nx*ny, sizeof(double));
    double* rho = (double*)calloc(max_num_steps, sizeof(double));

    double beta;
    double alpha;

    
    /* ----------------------------------------------------------------------------------------------------
    Initializing arrays -------------------------------------------------------------------------------- */

    double* error = (double*)calloc(max_num_steps, sizeof(double));

    double** f = (double**)calloc(ny, sizeof(double*));
    double** phi = (double**)calloc(ny, sizeof(double*));
    double** phi_new = (double**)calloc(ny, sizeof(double*));
    double** laplace_phi = (double**)calloc(ny, sizeof(double*));

    for (i = 0; i < Nx; i++) {
        f[i] = (double*)calloc(nx, sizeof(double));
        phi[i] = (double*)calloc(nx, sizeof(double));
        phi_new[i] = (double*)calloc(nx, sizeof(double));
        laplace_phi[i] = (double*)calloc(nx, sizeof(double));

    }

    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            f[j][i] = exp((-pow(x - 0.75, 2) - pow(y - 0.75, 2)) / pow(0.05, 2)) - exp((-pow(x - 0.25, 2) - pow(y - 0.25, 2)) / pow(0.05, 2));

            x += D_x;
            
        }

        y += D_y;
        x = D_x / 2;
    }

    /*=============================================== SOLVERS =================================================*/

    /* ----------------------------------------------------------------------------------------------------------
    Point Jacobi method -------------------------------------------------------------------------------- */

    if ( strcmp(method, "PJ") == 0 ) {

        f_norm = 0; /* compute f_norm */
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                f_norm = f_norm + pow(f[j][i], 2);
            }
        }

        f_norm = sqrt(f_norm);




        do {


            for (j = 1; j < Ny - 1; j++) {
                for (i = 1; i < Nx - 1; i++) {

                    phi_new[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 4 - f[j][i] / (4 * lambda);

                }

            }

            //update left boundary values
            for (j = 1; j < Ny - 1; j++) {
                i = 0;
                phi_new[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

            }

            //update right boundary values
            for (j = 1; j < Ny - 1; j++) {
                i = Nx - 1;
                phi_new[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

            }

            //update bottom boundary values 
            for (i = 1; i < Nx - 1; i++) {
                j = 0;
                phi_new[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
            }

            //update top boundary values
            for (i = 1; i < Nx - 1; i++) {
                j = Ny - 1;
                phi_new[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) / 3 - f[j][i] / (3 * lambda);
            }

            //update corner points
            phi_new[0][0] = (phi[1][0] + phi[0][1]) / 2 - f[0][0] / (2 * lambda);
            phi_new[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) / 2 - f[0][Nx - 1] / (2 * lambda);
            phi_new[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) / 2 - f[Ny - 1][0] / (2 * lambda);
            phi_new[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) / 2 - f[Ny - 1][Nx - 1] / (2 * lambda);



            //update phi matrix
            for (j = 0; j < Ny; j++) {
                for (i = 0; i < Nx; i++) {
                    phi[j][i] = phi_new[j][i];
                }
            }



            //compute laplace_p matrix
            for (j = 1; j < Ny - 1; j++) {
                for (i = 1; i < Nx - 1; i++) {
                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 4 * lambda * phi[j][i];
                }
            }
            for (j = 1; j < Ny - 1; j++) {
                i = 0;
                laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
            }
            for (j = 1; j < Ny - 1; j++) {
                i = Nx - 1;
                laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
            }
            for (i = 1; i < Nx - 1; i++) {
                j = 0;
                laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
            }
            for (i = 1; i < Nx - 1; i++) {
                j = Ny - 1;
                laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) * lambda - phi[j][i] * (3 * lambda);
            }
            laplace_phi[0][0] = (phi[1][0] + phi[0][1]) * lambda - phi[0][0] * (2 * lambda);
            laplace_phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) * lambda - phi[0][Nx - 1] * (2 * lambda);
            laplace_phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) * lambda - phi[Ny - 1][0] * (2 * lambda);
            laplace_phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) * lambda - phi[Ny - 1][Nx - 1] * (2 * lambda);



            //compute the norm
            laplace_phi_minus_f_norm = 0;

            for (j = 0; j < Ny; j++) {
                for (i = 0; i < Nx; i++) {
                    laplace_phi_minus_f_norm = laplace_phi_minus_f_norm + pow((laplace_phi[j][i] - f[j][i]), 2);
                }
            }

            

            laplace_phi_minus_f_norm = sqrt(laplace_phi_minus_f_norm);


            if (f_norm == 0) {
                RHS = epsilon;
            }
            else {
                RHS = epsilon * f_norm;
            }


            //save the error here and break out of the loop if the max number of steps has been reached
            error[step] = laplace_phi_minus_f_norm;

            if (step == max_num_steps) {
                break;
            }

            step += 1; 

        } while (laplace_phi_minus_f_norm > RHS);


        // Impose condition that the integral over the domain is equal to zero
        integral = 0;
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                integral += phi[j][i] * D_x * D_y;
            }
        }

        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                phi[j][i] = phi[j][i] - integral / (Nx * Ny);
            }
        }

        print_now = print_current_data(step, laplace_phi, f, phi, error, Nx, Ny, method);
        printf("Data was printed for Point Jacobi method");




    }



    


    /* ----------------------------------------------------------------------------------------------------------
    Gauss Seidel method -------------------------------------------------------------------------------- */ 

    if ( strcmp(method, "GS") == 0 ) {
        
        
        f_norm = 0; /* compute f_norm */
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                f_norm = f_norm + pow(f[j][i], 2);
            }
        }

        f_norm = sqrt(f_norm);




        do {


            for (j = 1; j < Ny - 1; j++) {
                for (i = 1; i < Nx - 1; i++) {

                    phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 4 - f[j][i] / (4 * lambda);

                }

            }

            //update left boundary values
            for (j = 1; j < Ny - 1; j++) {
                i = 0;
                phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

            }

            //update right boundary values
            for (j = 1; j < Ny - 1; j++) {
                i = Nx - 1;
                phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

            }

            //update bottom boundary values 
            for (i = 1; i < Nx - 1; i++) {
                j = 0;
                phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
            }

            //update top boundary values
            for (i = 1; i < Nx - 1; i++) {
                j = Ny - 1;
                phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) / 3 - f[j][i] / (3 * lambda);
            }

            //update corner points
            phi[0][0] = (phi[1][0] + phi[0][1]) / 2 - f[0][0] / (2 * lambda);
            phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) / 2 - f[0][Nx - 1] / (2 * lambda);
            phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) / 2 - f[Ny - 1][0] / (2 * lambda);
            phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) / 2 - f[Ny - 1][Nx - 1] / (2 * lambda);



            //compute laplace_p matrix
            for (j = 1; j < Ny - 1; j++) {
                for (i = 1; i < Nx - 1; i++) {
                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 4 * lambda * phi[j][i];
                }
            }
            for (j = 1; j < Ny - 1; j++) {
                i = 0;
                laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
            }
            for (j = 1; j < Ny - 1; j++) {
                i = Nx - 1;
                laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
            }
            for (i = 1; i < Nx - 1; i++) {
                j = 0;
                laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
            }
            for (i = 1; i < Nx - 1; i++) {
                j = Ny - 1;
                laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) * lambda - phi[j][i] * (3 * lambda);
            }
            laplace_phi[0][0] = (phi[1][0] + phi[0][1]) * lambda - phi[0][0] * (2 * lambda);
            laplace_phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) * lambda - phi[0][Nx - 1] * (2 * lambda);
            laplace_phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) * lambda - phi[Ny - 1][0] * (2 * lambda);
            laplace_phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) * lambda - phi[Ny - 1][Nx - 1] * (2 * lambda);



            //compute the norm
            laplace_phi_minus_f_norm = 0;

            for (j = 0; j < Ny; j++) {
                for (i = 0; i < Nx; i++) {
                    laplace_phi_minus_f_norm = laplace_phi_minus_f_norm + pow((laplace_phi[j][i] - f[j][i]), 2);
                }
            }

            

            laplace_phi_minus_f_norm = sqrt(laplace_phi_minus_f_norm);


            if (f_norm == 0) {
                RHS = epsilon;
            }
            else {
                RHS = epsilon * f_norm;
            }


            //save the error here and break out of the loop if the max number of steps has been reached
            error[step] = laplace_phi_minus_f_norm;

            if (step == max_num_steps) {
                break;
            }

            step += 1; 

        } while (laplace_phi_minus_f_norm > RHS);


        // Impose condition that the integral over the domain is equal to zero
        integral = 0;
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                integral += phi[j][i] * D_x * D_y;
            }
        }

        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                phi[j][i] = phi[j][i] - integral / (Nx * Ny);
            }
        }

        print_now = print_current_data(step, laplace_phi, f, phi, error, Nx, Ny, method);
        printf("Data was printed for Gauss-Seidel method");



    }


    /* ----------------------------------------------------------------------------------------------------------
    Successive over-relaxation method -------------------------------------------------------------------------------- */

    if ( strcmp(method, "SOR") == 0 ) {
        
        f_norm = 0; /* compute f_norm */
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                f_norm = f_norm + pow(f[j][i], 2);
            }
        }

        f_norm = sqrt(f_norm);




        do {


            for (j = 1; j < Ny - 1; j++) {
                for (i = 1; i < Nx - 1; i++) {

                    phi_GS = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 4 - f[j][i] / (4 * lambda);
                    d_GS = phi_GS - phi[j][i];
                    phi[j][i] += omega * d_GS;

                }

            }

            //update left boundary values
            for (j = 1; j < Ny - 1; j++) {
                i = 0;
                phi_GS = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
                d_GS = phi_GS - phi[j][i];
                phi[j][i] += omega * d_GS;

            }

            //update right boundary values
            for (j = 1; j < Ny - 1; j++) {
                i = Nx - 1;
                phi_GS = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
                d_GS = phi_GS - phi[j][i];
                phi[j][i] += omega * d_GS;

            }

            //update bottom boundary values 
            for (i = 1; i < Nx - 1; i++) {
                j = 0;
                phi_GS = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
                d_GS = phi_GS - phi[j][i];
                phi[j][i] += omega * d_GS;

            }

            //update top boundary values
            for (i = 1; i < Nx - 1; i++) {
                j = Ny - 1;
                phi_GS = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) / 3 - f[j][i] / (3 * lambda);
                d_GS = phi_GS - phi[j][i];
                phi[j][i] += omega * d_GS;

            }

            //update corner points
            phi_GS = (phi[1][0] + phi[0][1]) / 2 - f[0][0] / (2 * lambda);
            d_GS = phi_GS - phi[0][0];
            phi[0][0] += omega * d_GS;

            phi_GS = (phi[1][Nx - 1] + phi[0][Nx - 2]) / 2 - f[0][Nx - 1] / (2 * lambda);
            d_GS = phi_GS - phi[0][Nx - 1];
            phi[0][Nx - 1] += omega * d_GS;

            phi_GS = (phi[Ny - 1][1] + phi[Ny - 2][0]) / 2 - f[Ny - 1][0] / (2 * lambda);
            d_GS = phi_GS - phi[Ny - 1][0];
            phi[Ny - 1][0] += omega * d_GS;

            phi_GS = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) / 2 - f[Ny - 1][Nx - 1] / (2 * lambda);
            d_GS = phi_GS - phi[Ny - 1][Nx - 1];
            phi[Ny - 1][Nx - 1] += omega * d_GS;



            //compute laplace_p matrix
            for (j = 1; j < Ny - 1; j++) {
                for (i = 1; i < Nx - 1; i++) {
                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 4 * lambda * phi[j][i];
                }
            }
            for (j = 1; j < Ny - 1; j++) {
                i = 0;
                laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
            }
            for (j = 1; j < Ny - 1; j++) {
                i = Nx - 1;
                laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
            }
            for (i = 1; i < Nx - 1; i++) {
                j = 0;
                laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
            }
            for (i = 1; i < Nx - 1; i++) {
                j = Ny - 1;
                laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) * lambda - phi[j][i] * (3 * lambda);
            }
            laplace_phi[0][0] = (phi[1][0] + phi[0][1]) * lambda - phi[0][0] * (2 * lambda);
            laplace_phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) * lambda - phi[0][Nx - 1] * (2 * lambda);
            laplace_phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) * lambda - phi[Ny - 1][0] * (2 * lambda);
            laplace_phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) * lambda - phi[Ny - 1][Nx - 1] * (2 * lambda);



            //compute the norm
            laplace_phi_minus_f_norm = 0;

            for (j = 0; j < Ny; j++) {
                for (i = 0; i < Nx; i++) {
                    laplace_phi_minus_f_norm = laplace_phi_minus_f_norm + pow((laplace_phi[j][i] - f[j][i]), 2);
                }
            }

            

            laplace_phi_minus_f_norm = sqrt(laplace_phi_minus_f_norm);


            if (f_norm == 0) {
                RHS = epsilon;
            }
            else {
                RHS = epsilon * f_norm;
            }


            //save the error here and break out of the loop if the max number of steps has been reached
            error[step] = laplace_phi_minus_f_norm;

            if (step == max_num_steps) {
                break;
            }

            step += 1; 

        } while (laplace_phi_minus_f_norm > RHS);


        // Impose condition that the integral over the domain is equal to zero
        integral = 0;
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                integral += phi[j][i] * D_x * D_y;
            }
        }

        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                phi[j][i] = phi[j][i] - integral / (Nx * Ny);
            }
        }

        print_now = print_current_data(step, laplace_phi, f, phi, error, Nx, Ny, method);
        printf("Data was printed for Successive over-relaxation method");


    }


    /* ----------------------------------------------------------------------------------------------------------
    Conjugate gradient method -------------------------------------------------------------------------------- */

    if ( strcmp(method, "CG") == 0 ) {

        /* compute A matrix */

        //top left point
        A[0][0] = -2*lambda;
        A[0][1] = lambda;
        A[0][Nx] = lambda;

        //top right point
        A[Nx - 1][Nx - 1] = -2 * lambda;
        A[Nx - 1][Nx - 2] = lambda;
        A[Nx - 1][2 * Nx - 1] = lambda;
        
        //bottom left point
        A[(Ny - 1)* Nx][(Ny - 1)* Nx] = -2 * lambda;
        A[(Ny - 1)* Nx][(Ny - 1)* Nx + 1] = lambda;
        A[(Ny - 1)* Nx][(Ny - 2)* Nx] = lambda;
        
        //bottom right point
        A[Ny*Nx - 1][Ny*Nx - 1] = -2 * lambda; 
        A[Ny*Nx - 1][Ny*Nx - 2] = lambda;
        A[Ny*Nx - 1][(Ny-1)*Nx - 1] = lambda;


        //top boudary points
        for (i = 1; i < Nx - 1; i++) {
            A[i][i] = -3 * lambda;
            A[i][i - 1] = lambda;
            A[i][i + 1] = lambda;
            A[i][i + Nx] = lambda; 
        }

        //bottom boundary points 
        for (i = 1; i < Nx - 1; i++) {
            A[Nx*(Ny-1) + i][Nx*(Ny-1) + i] = -3 * lambda;
            A[Nx*(Ny-1) + i][Nx*(Ny-1) + i + 1] = lambda;
            A[Nx*(Ny-1) + i][Nx*(Ny-1) + i - 1] = lambda;
            A[Nx*(Ny-1) + i][Nx*(Ny-1) + i - Nx] = lambda;
        }

        //left boundary points
        for (i = 1; i < Ny - 1; i++) {
            A[i*Nx][i*Nx] = - 3 * lambda;
            A[i*Nx][i*Nx + 1] = lambda;
            A[i*Nx][(i - 1)*Nx] = lambda;
            A[i*Nx][(i + 1)*Nx] = lambda;
        }


        //right boundary points
        for (i = 2; i < Ny; i++) {
            A[i*Nx - 1][i*Nx - 1] = -3 * lambda;
            A[i*Nx - 1][i*Nx - 2] = lambda;
            A[i*Nx - 1][(i - 1)*Nx - 1] = lambda;
            A[i*Nx - 1][(i + 1)*Nx - 1] = lambda;
        }


        //interior points
        for (j = 1; j < Ny - 1; j++) {
            for (i = 1; i < Nx - 1; i++) {
                A[j * Nx + i][j * Nx + i] = -4 * lambda;
                A[j * Nx + i][j * Nx + i + 1] = lambda;
                A[j * Nx + i][j * Nx + i - 1] = lambda;
                A[j * Nx + i][j * Nx + i + Nx] = lambda;
                A[j * Nx + i][j * Nx + i - Nx] = lambda;
            }
        }



        /* compute f_vec vector */

        for (j = 0; j < Ny; j++) {
            for(i = 0; i < Nx; i++) {
                f_vec[i + Nx*j] = f[j][i];
            }
        }

        /*  USED THIS TO PRINT A AND MAKE SURE IT LOOKS GOOD... AND IT LOOKS RIGHT
        //print array to make sure it looks good
        for (j = 0; j < Ny*Nx; j++) {
            for (i = 1; i < Ny*Nx; i++) {
                printf("%f ", A[j][i]);
            }
            printf(" \n");
        }
        //save A to text file to look at it
        char filename1[] = "A_matrix.txt";
        FILE* fpointer1 = fopen(filename1, "w");
        for (j = 0; j < Nx*Ny; j++) {
            for (i = 0; i < Nx*Ny; i++) {
                fprintf(fpointer1, "%.20lf ", A[j][i]);
            }
            fprintf(fpointer1, "\n");
        }
        fclose(fpointer1);
        */


        // Initialization
        for (i = 0; i < Nx*Ny; i++) {
            sum = 0;
            for (j = 0; j < Nx*Ny; j++) {
                sum += A[i][j]*phi_vec[j];
            }

            r[i] = f_vec[i] - sum;
        }
        
        
        for (i = 0; i < Nx*Ny; i++) {
            rho[0] += pow(r[i], 2);
        }
        

        

        do {

            error[step] = sqrt(rho[step]);
            printf("Step %d has error %f \n", step, error[step]);
            step += 1;

            if (step == 1) {

                for (i = 0; i < Nx*Ny; i++) {
                    d[i] = r[i];
                }

            }
            else {
                beta = rho[step - 1] / rho[step - 2];

                for (i = 0; i < Nx*Ny; i++) {
                    d[i] = r[i] + beta * d[i];
                }

            }


            for (i = 0; i < Nx*Ny; i++) {
                epsilon_vec[i] = 0;
                for (j = 0; j < Nx*Ny; j++) {
                    epsilon_vec[i] += A[i][j] * d[j];
                }
            }


            sum = 0;
            for (i = 0; i < Nx*Ny; i++) {
                sum += d[i]*epsilon_vec[i];
            }
            alpha = rho[step - 1] / sum;

            

            for (i = 0; i < Nx*Ny; i++) {
                phi_vec[i] = phi_vec[i] + alpha * d[i];
                r[i] = r[i] - alpha * epsilon_vec[i];
            }

            for (i = 0; i < Nx * Ny; i++) {
                rho[step] += r[i] * r[i];
            }

            
            

        } while ( sqrt(rho[step]) > epsilon * sqrt(rho[0]) );
        


        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                phi[j][i] = phi_vec[j*Nx + i];
            }
        }

        //print final error
        error[step] = sqrt(rho[step]);


        //compute laplace_p matrix (ONLY NECESSARY FOR VISUALIZATION)
        for (j = 1; j < Ny - 1; j++) {
            for (i = 1; i < Nx - 1; i++) {
                laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 4 * lambda * phi[j][i];
            }
        }
        for (j = 1; j < Ny - 1; j++) {
            i = 0;
            laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }
        for (j = 1; j < Ny - 1; j++) {
            i = Nx - 1;
            laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }
        for (i = 1; i < Nx - 1; i++) {
            j = 0;
            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }
        for (i = 1; i < Nx - 1; i++) {
            j = Ny - 1;
            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }
        laplace_phi[0][0] = (phi[1][0] + phi[0][1]) * lambda - phi[0][0] * (2 * lambda);
        laplace_phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) * lambda - phi[0][Nx - 1] * (2 * lambda);
        laplace_phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) * lambda - phi[Ny - 1][0] * (2 * lambda);
        laplace_phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) * lambda - phi[Ny - 1][Nx - 1] * (2 * lambda);



        // Impose condition that the integral over the domain is equal to zero
        integral = 0;
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                integral += phi[j][i] * D_x * D_y;
            }
        }

        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                phi[j][i] = phi[j][i] - integral / (Nx * Ny);
            }
        }

        print_now = print_current_data(step, laplace_phi, f, phi, error, Nx, Ny, method);
        printf("Data was printed for conjugate gradient method");



    }


    /* ----------------------------------------------------------------------------------------------------------
    Multigrid method -------------------------------------------------------------------------------- */

    if ( strcmp(method, "MG") == 0 ) {

        printf("Start MG method\n");

        do {
            
            printf("Iteration step %d\n", step);

            MG_recursion(f, phi, Nx, Ny, epsilon);
            
            /* compute laplace_p matrix */
            for (j = 1; j < Ny - 1; j++) {
                for (i = 1; i < Nx - 1; i++) {
                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 4 * lambda * phi[j][i];
                }
            }
            for (j = 1; j < Ny - 1; j++) {
                i = 0;
                laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
            }
            for (j = 1; j < Ny - 1; j++) {
                i = Nx - 1;
                laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
            }
            for (i = 1; i < Nx - 1; i++) {
                j = 0;
                laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
            }
            for (i = 1; i < Nx - 1; i++) {
                j = Ny - 1;
                laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) * lambda - phi[j][i] * (3 * lambda);
            }

            laplace_phi[0][0] = (phi[1][0] + phi[0][1]) * lambda - phi[0][0] * (2 * lambda);
            laplace_phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) * lambda - phi[0][Nx - 1] * (2 * lambda);
            laplace_phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) * lambda - phi[Ny - 1][0] * (2 * lambda);
            laplace_phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) * lambda - phi[Ny - 1][Nx - 1] * (2 * lambda);


            //compute the norm
            laplace_phi_minus_f_norm = 0;

            for (j = 0; j < Ny; j++) {
                for (i = 0; i < Nx; i++) {
                    laplace_phi_minus_f_norm = laplace_phi_minus_f_norm + pow((laplace_phi[j][i] - f[j][i]), 2);
                }
            }

            
            laplace_phi_minus_f_norm = sqrt(laplace_phi_minus_f_norm);


            if (f_norm == 0) {
                RHS = epsilon;
            }
            else {
                RHS = epsilon * f_norm;
            }


            //save the error here and break out of the loop if the max number of steps has been reached
            error[step] = laplace_phi_minus_f_norm;

            if (step == max_num_steps) {
                break;
            }


            step += 1;

        } while (laplace_phi_minus_f_norm > RHS);
    }

    /* Impose condition that the integral over the domain is equal to zero */
    integral = 0;
    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            integral += phi[j][i] * D_x * D_y;
        }
    }

    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            phi[j][i] = phi[j][i] - integral / (Nx * Ny);
        }
    }

    print_now = print_current_data(step, laplace_phi, f, phi, error, Nx, Ny, method);
    printf("Data was printed for Multigrid method\n");




    /* ----------------------------------------------------------------------------------------------------------
    Freeing variables -------------------------------------------------------------------------------- */

    for (j = 0; j < Ny; j++) {
        free(f[j]);
        free(phi[j]);
        free(phi_new[j]);
        free(laplace_phi[j]);
    }

    free(f);
    free(phi);
    free(phi_new);
    free(laplace_phi);
    free(error);

    printf("end program\n");

    return 0;

}