#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEBUG -1
#define MemoryError -3
#define SuccessSolved 2
#define FinishSolved 3
#define Solving 0

#define MAX_ITER 10000
#define TOLERANCE 1e-6

int jacobi_iteration(double *A, double *b, double *x, int n);

int main() {
    int n = 3;

    double *A = (double*)malloc(n * n * sizeof(double));
    double b[3] = {20, 33, 36};
    double x[3] = {0, 0, 0};

    A[0 * n + 0] = 8; A[0 * n + 1] = -3; A[0 * n + 2] = 2;
    A[1 * n + 0] = 4; A[1 * n + 1] = 11; A[1 * n + 2] = -1;
    A[2 * n + 0] = 6; A[2 * n + 1] = 3; A[2 * n + 2] = 12;

    jacobi_iteration(A, b, x, n);

    printf("Solution x = [ ");
    for (int i = 0; i < n; i++) {
        printf("%lf ", x[i]);
    }
    printf("]\n");

    free(A);

    return 0;
}

int jacobi_iteration(double *A, double *b, double *x, int n) {

    int return_status = Solving;

    double *x_solution = (double*)malloc(n * sizeof(double)); 

    if (x_solution == NULL) {
        printf("Error Memory Allocation!\n");
        return MemoryError;
    }

    int iter = 0;
    while (iter < MAX_ITER) {
        for (int i = 0; i < n; i++) {
            x_solution[i] = b[i];
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    x_solution[i] -= A[i * n + j] * x[j];
                }
            }
            x_solution[i] /= A[i * (n + 1)];  // div A[i][i]
        }

        double max_error = 0.0;
        for (int i = 0; i < n; i++) {
            max_error = fmax(max_error, fabs(x_solution[i] - x[i]));
        }

        if (max_error < TOLERANCE) {
            return_status = SuccessSolved;
            break;
        }

        // update solution
        for (int i = 0; i < n; i++) {
            x[i] = x_solution[i];
        }

        iter++;

        if(iter == MAX_ITER - 1) {
            return_status = FinishSolved;
        }
    }

    printf("Iteration Finished: Iteartion:%d\n", iter);

    free(x_solution);

    return return_status;
}