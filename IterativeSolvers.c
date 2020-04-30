#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <string.h>



//function used to print the data to a text file
int print_current_data(int step, double** u, double** v, double** p, int NX, int NY) {
    int i, j;


    char filename1[25] = "step_";
    itoa(step, filename1 + 4, 10);
    strcat(filename1, "_u_data.txt");

    FILE* fpointer1 = fopen(filename1, "w");
    for (j = 0; j < NY; j++) {
        for (i = 0; i < NX; i++) {
            fprintf(fpointer1, "%.20lf ", u[j][i]);
        }
        fprintf(fpointer1, "\n");
    }
    fclose(fpointer1);



    char filename2[25] = "step_";
    itoa(step, filename2 + 4, 10);
    strcat(filename2, "_v_data.txt");

    FILE* fpointer2 = fopen(filename2, "w");
    for (j = 0; j < NY; j++) {
        for (i = 0; i < NX; i++) {
            fprintf(fpointer2, "%.20lf ", v[j][i]);
        }
        fprintf(fpointer2, "\n");
    }
    fclose(fpointer2);




    char filename3[25] = "step_";
    itoa(step, filename3 + 4, 10);
    strcat(filename3, "_p_data.txt");

    FILE* fpointer3 = fopen(filename3, "w");
    for (j = 0; j < NY; j++) {
        for (i = 0; i < NX; i++) {
            fprintf(fpointer3, "%.20lf ", p[j][i]);
        }
        fprintf(fpointer3, "\n");
    }
    fclose(fpointer3);

    return 1;
}













int main()
{


}