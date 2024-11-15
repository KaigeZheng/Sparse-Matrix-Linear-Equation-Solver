#include <stdio.h>
#include <stdlib.h>
#include "compute_det.h"

#define MAX_SIZE 1000

#define ErrorOpeningFIle -1
#define SuccessReadMatrix 1

int readMatrix(const char *filename, double *matrix, int n, int m) {
    FILE *file = fopen(filename, "r");
    if(file == NULL) {
        printf("Error opening file!\n");
        return -1;
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
    double matrix[MAX_SIZE * MAX_SIZE];
    int n = 39;
    int m = 39;

    const char *filename = "./data/bcspwr01.mtx";

    int status = readMatrix(filename, matrix, n, m);
    
    int cnt = 0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (matrix[(i * m) + j] != 0) {  
                ++cnt;
                printf("matrix[%d][%d] = %lf\n", i, j, matrix[(i * m) + j]);
                //if(i != j) cnt--;
            }
        }
    }

    printf("Total Num(Nonzeros) = %d\n", cnt);

    double det = compute_det(matrix, n, DEFAULT, DEFAULT);
    printf("det = %lf", det);
    return 0;
}