#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <fstream>
#include <ctime>
#include <random>
#include <vector>

// Read Sparse Matrix(Matrix Market format + pattern matrix)
Eigen::SparseMatrix<double> readSparseMatrix(const std::string &path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << path << std::endl;
        exit(1);
    }

    std::string line;
    
    while (std::getline(file, line)) {
        if (line[0] != '%') {
            break;
        }
    }

    // 读取矩阵维度和非零元素数
    int rows, cols, nonZeros;
    std::stringstream ss(line);
    ss >> rows >> cols >> nonZeros;

    Eigen::SparseMatrix<double> mat(rows, cols);
    std::vector<Eigen::Triplet<double>> triplets;

    // 读取非零元素
    while (std::getline(file, line)) {
        if (line[0] == '%') continue;
        std::stringstream tripletStream(line);
        int row, col;
        tripletStream >> row >> col;
        
        // 矩阵是1-based indexing, 转为0-based
        row--; col--;

        // 将非零元素的值设为1
        triplets.push_back(Eigen::Triplet<double>(row, col, 1.0));
        if (row != col) {
            // 对称矩阵，插入对称位置
            triplets.push_back(Eigen::Triplet<double>(col, row, 1.0));
        }
    }

    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

// 判断矩阵是否奇异
bool isSingular(const Eigen::SparseMatrix<double>& mat) {
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(mat);
    return solver.info() != Eigen::Success;
}

int main() {
    std::string path = "./data/bcspwr01/bcspwr01.mtx";
    std::cout << "Reading matrix from: " << path << std::endl;

    // 读取稀疏矩阵
    Eigen::SparseMatrix<double> A = readSparseMatrix(path);
    
    // 判断矩阵是否奇异
    if (isSingular(A)) {
        std::cout << "Matrix is singular!" << std::endl;
        return -1;
    } else {
        std::cout << "Matrix is non-singular." << std::endl;
    }

    // 随机生成同阶向量b
    int n = A.rows();
    Eigen::VectorXd b = Eigen::VectorXd::Random(n);

    // 计时开始
    clock_t start_time = clock();

    // 使用稀疏LU分解求解Ax = b
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Decomposition failed!" << std::endl;
        return -1;
    }
    Eigen::VectorXd x = solver.solve(b);

    // 计时结束
    clock_t end_time = clock();
    double duration = double(end_time - start_time) / CLOCKS_PER_SEC;
    
    // 输出结果
    std::cout << "Solution x: \n" << x.transpose() << std::endl;
    std::cout << "Time taken: " << duration << " seconds" << std::endl;

    return 0;
}
