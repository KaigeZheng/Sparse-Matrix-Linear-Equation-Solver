#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <fstream>
#include <ctime>
#include <random>
#include <vector>
#include <cmath>
#include <sstream>

// 读取稀疏矩阵（矩阵市场格式 + 图模式矩阵）
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
        double value;
        tripletStream >> row >> col >> value;

        // 矩阵是1-based indexing, 转为0-based
        row--; col--;

        // 将非零元素添加到矩阵
        triplets.push_back(Eigen::Triplet<double>(row, col, value));
    }

    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
}

// 雅可比迭代法求解稀疏线性方程 Ax = b
int jacobiIteration(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b, Eigen::VectorXd& x, double tol = 1e-6, int maxIter = 1000) {
    int n = A.rows();
    Eigen::VectorXd x_new = x;

    // 对角线矩阵D
    Eigen::VectorXd D(n);
    for (int i = 0; i < n; ++i) {
        D(i) = A.coeff(i, i);
    }

    // 剩余矩阵R = A - D
    Eigen::SparseMatrix<double> R = A;
    R.diagonal().setZero();  // 设置对角线为零，即R = A - D

    int iter = 0;
    double error = tol + 1; // 初始误差大于容忍误差

    while (iter < maxIter && error > tol) {
        // x_new = D^(-1) * (b - R * x)
        Eigen::VectorXd r = b - R * x;

        // 更新解向量 x_new
        for (int i = 0; i < n; ++i) {
            if (D(i) != 0) {
                x_new(i) = r(i) / D(i); // 更新解
            }
        }

        // 计算误差
        error = (x_new - x).norm();

        // 更新解向量
        x = x_new;
        ++iter;
    }

    return iter; // 返回迭代次数
}

int main() {
    std::string path = "../data/1138_bus/1138_bus.mtx";
    std::cout << "Reading matrix from: " << path << std::endl;

    // 读取稀疏矩阵
    Eigen::SparseMatrix<double> A = readSparseMatrix(path);

    // 随机生成同阶向量b
    int n = A.rows();
    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);
    //Eigen::VectorXd b(3);
    //b << 20, 33, 36;

    // 初始猜测解x
    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);

    // 计时开始
    clock_t start_time = clock();

    // 使用雅可比迭代法求解Ax = b
    int iter = jacobiIteration(A, b, x);

    // 计时结束
    clock_t end_time = clock();
    double duration = double(end_time - start_time) / CLOCKS_PER_SEC;

    // 输出结果
    std::cout << "Solution x: \n" << x.transpose() << std::endl;
    std::cout << "Iterations: " << iter << std::endl;
    std::cout << "Time taken: " << duration << " seconds" << std::endl;

    return 0;
}

