#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void generate_random_vector(double *vec, int n, double min, double max) {
    srand(time(NULL));

    for(int i = 0; i < n; ++i) {
        vec[i] = min + (rand() / (RAND_MAX + 1.0)) * (max - min);
    }
}

void generate_random_matrix(double *vec, int m, int n, double min, double max) {
    srand(time(NULL));

    for(int i = 0; i < m; ++i) {
        for(int j = 0; j < n; ++j) {
            matrix[i * n + j] = min + (rand() / (RAND_MAX + 1.0)) * (max - min);
        }
    }
}