#include <stdio.h>
#include <stdlib.h>

#define DEFAULT -1

double compute_det(double* matrix, int n, int expand_row = -1, int expand_col = -1) {
    if(n == 1) {
        return *matrix;
    }

    double det = 0.0;
    double* submatrix = (double*)malloc((n - 1) * (n - 1) * sizeof(double));

    if(expand_row >= 0) {
        for(int col = 0; col < n; ++col) {
            int submatrix_row = 0;
            for(int i = 0; i < n; ++i) {
                if(i == expand_row) continue;
                int submatrix_col = 0;
                for(int j = 0; j < n; ++j) {
                    if(j == col) continue;
                    submatrix[submatrix_row * (n - 1) + submatrix_col] = matrix[i * n + j];
                    ++submatrix_col;
                }
                ++submatrix_row;
            }
            double subdet = compute_det(submatrix, n - 1, DEFAULT, DEFAULT);
            det += (((col + expand_row) % 2 == 0 ? 1: -1) * matrix[expand_row * n + col] * subdet);
        }
    } else if (expand_col >= 0) {
        for(int row = 0; row < n; row++) {
            int submatrix_col = 0;
            for(int j = 0; j < n; ++j) {
                if(j == expand_col) continue;
                int submatrix_row = 0;
                for(int i = 0; i < n; ++i) {
                    if(i == row) continue;
                    submatrix[submatrix_row * (n - 1) + submatrix_col] = matrix[i * n + j];
                    ++submatrix_row;
                }
                ++submatrix_col;
            }
            double subdet = compute_det(submatrix, n-1, DEFAULT, DEFAULT);
            det += (((row + expand_col) % 2 == 0 ? 1: -1) * matrix[row * n + expand_col] * subdet);
        }
    } else {
        // By default, perform the Laplace expansion along the first row.
        for(int col = 0; col < n; ++col) {
            int submatrix_row = 0;
            for(int i = 1; i < n; ++i) {
                int submatrix_col = 0;
                for(int j = 0; j < n; ++j) {
                    if(j == col) continue; // Skip the current column
                    submatrix[submatrix_row * (n - 1) + submatrix_col] = matrix[i * n + j];
                    ++submatrix_col;
                }
                ++submatrix_row;
            }
            double subdet = compute_det(submatrix, n - 1, DEFAULT, DEFAULT);
            det += ((col % 2 == 0 ? 1: -1) * matrix[col] * subdet);
        }
    }

    free(submatrix);
    return det;
}

int main() {
    int n, expand_row, expand_col;
    n = 3;
    double matrix[9] = {1,2,3,4,5,6,7,8,8};
    double det = compute_det(matrix, n, DEFAULT, DEFAULT);
    // The correct result should be 3.000000.
    printf("det = %lf\n", det);
    return 0;
}