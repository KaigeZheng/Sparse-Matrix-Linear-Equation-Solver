#include <stdio.h>
#include <stdlib.h>
#include "compute_det.h"
#include "solver.h"
#include "generator.h"

#define MAX_SIZE 1000

#define ErrorOpeningFile -2
#define SuccessReadMatrix 1

int readMatrix(const char *filename, double *matrix, int n, int m) {
    FILE *file = fopen(filename, "r");
    if(file == NULL) {
        printf("Error opening file!\n");
        return ErrorOpeningFile;
    }

    // Skip comments
    char line[256];

    while(fgets(line, sizeof(line), file)) {
        if(line[0] != '%') break;
    }

    int rows, cols, nonzeros;
    sscanf(line, "%d %d %d", &rows, &cols, &nonzeros);

    int i, j;
    double value;

    for(int k = 0; k < nonzeros; ++k) {
        fgets(line, sizeof(line), file);
        int numItems = sscanf(line, "%d %d %lf", &i, &j, &value);
        if(numItems == 2) {
            value = 1.0;
        }
        matrix[(i - 1) * m + (j - 1)] = value;
        if(i != j) {
            matrix[(j - 1) * m + (i - 1)] = value;
        }
    }

    fclose(file);
    return SuccessReadMatrix;
}

int main(){
    double A[MAX_SIZE * MAX_SIZE];
    double b[MAX_SIZE];
    double x[MAX_SIZE];
    int n = 39;
    int m = 39;

    const char *filename = "./data/bcspwr01/bcspwr01.mtx";

    int status = readMatrix(filename, A, n, m);
    generate_random_vector(b, n, 0.0, 10.0);
    
    int cnt = 0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (A[(i * m) + j] != 0) {  
                ++cnt;
                printf("A[%d][%d] = %lf\n", i, j, A[(i * m) + j]);
                //if(i != j) cnt--;
            }
        }
    }

    printf("Total Num(Nonzeros) = %d\n", cnt);

    double det = compute_det_gaussion(A, n);

    printf("det = %lf\n", det);

    for(int i = 0; i < n; ++i) printf("b[%d] = %lf\n", i, b[i]);

    if(det != 0) {
        for(int i = 0; i < n; ++i) {
            x[i] = 0.0;
        }
        jacobi_iteration(A, b, x, n);
        printf("Solution x = [ \n");
        for (int i = 0; i < n; i++) {
            printf("%lf \n", x[i]);
        }
        printf("]\n");
    }
    return 0;
}