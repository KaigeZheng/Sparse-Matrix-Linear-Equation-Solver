#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <ctime>
#include <cmath>

#define SuccessReadMatrix 0
#define ErrorOpeningFile -1

typedef struct {
    double *values;       // 存储非零元素
    int *columns;         // 存储非零元素对应的列索引
    int *rowPointers;     // 存储每行的起始位置
    int numRows;          // 矩阵的行数
    int numCols;          // 矩阵的列数
    int numNonZeros;      // 矩阵中的非零元素个数
} CSRMatrix;

double compute_det(double** matrix, int n);
int readMatrixToArray(const char *filename, double ***matrix, int *numRows, int *numCols);
int convertToCSR(double **matrix, CSRMatrix *csrMatrix, int numRows, int numCols);
void printCSRMatrix(const CSRMatrix *csrMatrix);
void generate_b(double *b, int size);
int jacobiIterationCSR(const CSRMatrix *csrMatrix, const double *b, double *x, double tol, int maxIter);

int main() {
    const char *filename = "./data/msc01440/msc01440.mtx";// bratu3d/bratu3d.mtx
    double **matrix = NULL;
    CSRMatrix csrMatrix;
    int numRows, numCols;

    // Load data into 2-dim array
    int result = readMatrixToArray(filename, &matrix, &numRows, &numCols);
    if (result != SuccessReadMatrix) {
        printf("Error reading matrix!\n");
        return -1;
    }

    // double det = compute_det(matrix, numRows);
    // printf("Matrix Determinant: %f\n", det);

    // if(det == 0) {
    //     printf("A singulary matrix!\n");
    //     return -1;
    // }

    // Convert 2-dim array into CSR format
    result = convertToCSR(matrix, &csrMatrix, numRows, numCols);
    if (result != SuccessReadMatrix) {
        printf("Error converting to CSR format!\n");
        return -1;
    }

    // 输出CSR矩阵
    printCSRMatrix(&csrMatrix);
    
    // 生成向量b
    double *b = (double *)malloc(numRows * sizeof(double));
    generate_b(b, numRows);
    //b[0] = 20; b[1] = 33; b[2] = 36;

    // 初始化x
    double *x = (double *)malloc(numRows * sizeof(double));
    for (int i = 0; i < numRows; ++i) {
        x[i] = 0.0;
    }

    // 计时
    clock_t start_time = clock();
    
    // 雅可比迭代 Input:A, b Output: x
    double tol = 1e-6;  // Tolerance for convergence
    int maxIter = 1000000;  // Maximum number of iterations
    int iter = jacobiIterationCSR(&csrMatrix, b, x, tol, maxIter);

    // 结束计时
    clock_t end_time = clock();
    double duration = double(end_time - start_time) / CLOCKS_PER_SEC;

    printf("Solution vector x:\n");
    for (int i = 0; i < numRows; ++i) {
        printf("%f ", x[i]);
    }
    printf("\n");
    // 验证结果正确性
    std::cout << "Iterations: " << iter << std::endl;
    std::cout << "Time taken: " << duration << " seconds" << std::endl;

    // 清理内存
    for (int i = 0; i < numRows; ++i) {
        free(matrix[i]);
    }

    free(matrix);
    free(csrMatrix.values);
    free(csrMatrix.columns);
    free(csrMatrix.rowPointers);

    return 0;
}

// Load .mtx format data into 2-dim array
int readMatrixToArray(const char *filename, double ***matrix, int *numRows, int *numCols) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error opening file!\n");
        return ErrorOpeningFile;
    }

    // Skip comments
    char line[256];
    while (fgets(line, sizeof(line), file)) {
        if (line[0] != '%') break;
    }

    // 读取矩阵的维度
    int rows, cols, nonzeros;
    sscanf(line, "%d %d %d", &rows, &cols, &nonzeros);

    *numRows = rows;
    *numCols = cols;

    // allocate memory
    *matrix = (double **)malloc(rows * sizeof(double *));
    for (int i = 0; i < rows; ++i) {
        (*matrix)[i] = (double *)calloc(cols, sizeof(double));
    }

    int i, j;
    double value;
    while (fgets(line, sizeof(line), file)) {
        int num = sscanf(line, "%d %d %lf", &i, &j, &value);
        if(num == 2) {
            value = 1.0;
        }
        (*matrix)[i - 1][j - 1] = value;  // convert 1-based index into 0-based index
    }

    fclose(file);
    return SuccessReadMatrix;
}

// Convert 2-dim array into CSR format
int convertToCSR(double **matrix, CSRMatrix *csrMatrix, int numRows, int numCols) {
    int numNonZeros = 0;
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            if (matrix[i][j] != 0) {
                numNonZeros++;
            }
        }
    }

    csrMatrix->numRows = numRows;
    csrMatrix->numCols = numCols;
    csrMatrix->numNonZeros = numNonZeros;

    // 为CSR格式的数组分配内存
    csrMatrix->values = (double *)malloc(numNonZeros * sizeof(double));
    csrMatrix->columns = (int *)malloc(numNonZeros * sizeof(int));
    csrMatrix->rowPointers = (int *)malloc((numRows + 1) * sizeof(int));

    if (!csrMatrix->values || !csrMatrix->columns || !csrMatrix->rowPointers) {
        printf("Memory allocation failed!\n");
        return -2;
    }

    // 初始化rowPointers
    csrMatrix->rowPointers[0] = 0;

    // 填充CSR数组
    int currentNonZero = 0;
    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            if (matrix[i][j] != 0) {
                csrMatrix->values[currentNonZero] = matrix[i][j];
                csrMatrix->columns[currentNonZero] = j;
                currentNonZero++;
            }
        }
        csrMatrix->rowPointers[i + 1] = currentNonZero;  // 更新每行的起始位置
    }

    return SuccessReadMatrix;
}

void printCSRMatrix(const CSRMatrix *csrMatrix) {
    // 输出 Row pointers
    printf("Row pointers: ");
    for (int i = 0; i < csrMatrix->numRows + 1; ++i) {
        printf("%d ", csrMatrix->rowPointers[i]);
    }
    printf("\n");

    // 输出 Values
    printf("Values: ");
    for (int i = 0; i < csrMatrix->numNonZeros; ++i) {
        printf("%f ", csrMatrix->values[i]);
    }
    printf("\n");

    // 输出 Columns
    printf("Columns: ");
    for (int i = 0; i < csrMatrix->numNonZeros; ++i) {
        printf("%d ", csrMatrix->columns[i]);
    }
    printf("\n");
}

double compute_det(double** matrix, int n) {
    // 创建矩阵副本
    double** matrixCopy = new double*[n];
    for (int i = 0; i < n; ++i) {
        matrixCopy[i] = new double[n];
        for (int j = 0; j < n; ++j) {
            matrixCopy[i][j] = matrix[i][j];
        }
    }

    // 高斯消元法计算行列式
    double det = 1.0;
    
    for (int i = 0; i < n; ++i) {
        // 查找当前列最大的元素进行主元选择
        int maxRow = i;
        for (int j = i + 1; j < n; ++j) {
            if (abs(matrixCopy[j][i]) > abs(matrixCopy[maxRow][i])) {
                maxRow = j;
            }
        }

        // 交换当前行和最大行
        if (maxRow != i) {
            for (int j = 0; j < n; ++j) {
                double temp = matrixCopy[i][j];
                matrixCopy[i][j] = matrixCopy[maxRow][j];
                matrixCopy[maxRow][j] = temp;
            }
            det = -det;  // 行交换会改变行列式的符号
        }

        // 如果主元为0，行列式为0
        if (matrixCopy[i][i] == 0) {
            det = 0;
            break;
        }

        // 对当前行进行消元
        for (int j = i + 1; j < n; ++j) {
            double factor = matrixCopy[j][i] / matrixCopy[i][i];
            for (int k = i; k < n; ++k) {
                matrixCopy[j][k] -= factor * matrixCopy[i][k];
            }
        }

        det *= matrixCopy[i][i];
    }

    // Free memory
    for (int i = 0; i < n; ++i) {
        delete[] matrixCopy[i];
    }
    delete[] matrixCopy;

    return det;
}

void generate_b(double *b, int size) {
    for (int i = 0; i < size; ++i) {
        b[i] = 1.0 + (rand() % 1000000) / 100000.0;
        b[i] = 1.0;
    }
}

int jacobiIterationCSR(const CSRMatrix *csrMatrix, const double *b, double *x, double tol, int maxIter) {
    int numRows = csrMatrix->numRows;
    int iter;
    double *x_new = (double *)malloc(numRows * sizeof(double));  // 用于存储下一步迭代的解向量

    for (int i = 0; i < numRows; i++) {
        x_new[i] = 0.0;  // 初始化新解向量
    }

    for (iter = 0; iter < maxIter; iter++) {
        // 逐行计算新解
        for (int i = 0; i < numRows; i++) {
            double sum = 0.0;
            double diag = 0.0; // 存储对角线元素

            // 计算 x(i) 的新值：x(i) = (b(i) - Lx(i) - Ux(i)) / A(i,i)
            for (int j = csrMatrix->rowPointers[i]; j < csrMatrix->rowPointers[i + 1]; j++) {
                int col = csrMatrix->columns[j];
                double value = csrMatrix->values[j];

                if (col < i) {
                    sum += value * x[col];  // Lx(i)
                } else if (col > i) {
                    sum += value * x[col];  // Ux(i)
                } else {
                    diag = value;  // 对角线元素
                }
            }

            // 对角线元素直接从 b(i) 减去上三角和下三角的贡献
            x_new[i] = (b[i] - sum) / diag; // 对角线元素是 values[rowPointers[i]]
        }

        // 计算当前解与上一解的差值，检查是否满足容忍度
        double norm = 0.0;
        for (int i = 0; i < numRows; i++) {
            norm += (x_new[i] - x[i]) * (x_new[i] - x[i]);
        }
        norm = sqrt(norm);

        if (norm < tol) {
            printf("Converged in %d iterations.\n", iter + 1);
            break;
        }

        // 更新 x
        for (int i = 0; i < numRows; i++) {
            x[i] = x_new[i];
        }
    }

    free(x_new);
    return iter;
}
