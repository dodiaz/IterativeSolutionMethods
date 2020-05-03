#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <string.h>



//function used to print the data to a text file
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






int main() {

    /* ----------------------------------------------------------------------------------------------------------
    Initializing variables -------------------------------------------------------------------------------- */
    
    char method[] = "SOR";      //possible methods "PJ", "GS", "SOR", "CG", "MG"

    int i, j;
    int step = 0;
    int max_num_steps = 2000;   //Increase this number if the method isn't converging

    int Nx = 11;
    int Ny = 11;
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

    double omega = 1.3;   /* SOR method variables */
    double phi_GS;
    double d_GS;

    int print_now;

    /* ----------------------------------------------------------------------------------------------------------
    Initializing f -------------------------------------------------------------------------------- */

    double* error = (double*)calloc(max_num_steps, sizeof(double));

    double** f = (double**)calloc(ny, sizeof(double*));
    double** phi = (double**)calloc(ny, sizeof(double*));
    double** phi_new = (double**)calloc(ny, sizeof(double*));
    double** laplace_phi = (double**)calloc(ny, sizeof(double*));

    for (i = 0; i < nx; i++) {
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




    }


    /* ----------------------------------------------------------------------------------------------------------
    Multigrid method -------------------------------------------------------------------------------- */

    if ( strcmp(method, "MG") == 0 ) {
        



    }




    /* ----------------------------------------------------------------------------------------------------------
    Freeing variables -------------------------------------------------------------------------------- */

    for (i = 0; i < nx; i++) {
        free(f[i]);
        free(phi[i]);
        free(phi_new[i]);
        free(laplace_phi[i]);
    }

    free(f);
    free(phi);
    free(phi_new);
    free(laplace_phi);

    
}