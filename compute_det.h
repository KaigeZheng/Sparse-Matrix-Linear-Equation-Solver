#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEFAULT -1
//#define DEBUG 1

double compute_det_laplace(double* matrix, int n, int expand_row, int expand_col);
double compute_det_gaussion(double* matrix, int n);

#ifdef DEBUG
int main() {
    int n, expand_row, expand_col;
    n = 3;
    double matrix[9] = {1,2,3,4,5,6,7,8,8};
    double det = compute_det_gaussion(matrix, n);
    // The correct result should be 3.000000.
    printf("det = %lf\n", det);
    return 0;
}
#endif

double compute_det_laplace(double* matrix, int n, int expand_row = -1, int expand_col = -1) {
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
            if(matrix[expand_row * n + col] == 0) {
                det += 0;
            } else {
                double subdet = compute_det_laplace(submatrix, n - 1, DEFAULT, DEFAULT);
                det += (((col + expand_row) % 2 == 0 ? 1: -1) * matrix[expand_row * n + col] * subdet);
            }
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
            if(matrix[row * n + expand_col]){
                //det += 0;
                continue;
            } else {
                double subdet = compute_det_laplace(submatrix, n-1, DEFAULT, DEFAULT);
                det += (((row + expand_col) % 2 == 0 ? 1: -1) * matrix[row * n + expand_col] * subdet);
            }
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
            if(matrix[col] == 0) {
                det += 0;
            } else {
                double subdet = compute_det_laplace(submatrix, n - 1, DEFAULT, DEFAULT);
                det += ((col % 2 == 0 ? 1: -1) * matrix[col] * subdet);
            }
        }
    }

    free(submatrix);
    return det;
}

double compute_det_gaussion(double* matrix, int n) {
    if(n == 1) {
        return *matrix;
    }

    double det = 1.0;
    int sign = 1;
    double* temp_matrix = (double*)malloc(n * n * sizeof(double));
    for(int i = 0; i < n * n; ++i) {
        temp_matrix[i] = matrix[i];
    }

    for(int i = 0; i < n; ++i) {
        // searching the column pivot
        int max_row = i;
        for(int j = i + 1; j < n; ++j) {
            if(fabs(temp_matrix[j * n + i]) > fabs(temp_matrix[max_row * n + i])) {
                max_row = j;
            }
        }
        // if the pivot is 0, the determinant is 0.
        if(fabs(temp_matrix[max_row * n + i]) < 1e-9) {
            return 0.0;
        }
        // swap rows
        if(max_row != i) {
            for(int k = 0; k < n; ++k) {
                double temp = temp_matrix[i * n + k];
                temp_matrix[i * n + k] = temp_matrix[max_row * n + k];
                temp_matrix[max_row * n + k] = temp;
            }
            sign = -sign;
        }
        //perform elimination
        for(int j = i + 1; j < n; ++j) {
            double factor = temp_matrix[j * n + i] / temp_matrix[i * n + i];
            for(int k = i; k < n; ++k) {
                temp_matrix[j * n + k] -= factor * temp_matrix[i * n + k];
            }
        }
        // update the determinant
        det *= temp_matrix[i * n + i];
    }
    free(temp_matrix);
    return det * sign;
}