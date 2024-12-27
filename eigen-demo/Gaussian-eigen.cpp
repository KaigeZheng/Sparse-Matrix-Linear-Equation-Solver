#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <fstream>
#include <ctime>
#include <random>
#include <vector>
#include <cmath>

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
    double value = 1.0;
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

        if(tripletStream >> value){
            
        }

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

// 高斯-赛德尔迭代法求解稀疏线性方程 Ax = b
int gaussSeidelIteration(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b, Eigen::VectorXd& x, double tol = 1e-6, int maxIter = 1000) {
    int n = A.rows();
    Eigen::VectorXd x_new = x;

    // 高斯-赛德尔迭代
    int iter = 0;
    double error = tol + 1; // 初始误差大于容忍误差

    while (iter < maxIter && error > tol) {
        // 逐步更新x的每个元素
        for (int i = 0; i < n; ++i) {
            double sum = b(i);
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum -= A.coeff(i, j) * x_new(j);
                }
            }
            x_new(i) = sum / A.coeff(i, i); // 更新解
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
    std::string path = "./data/bcspwr01/bcspwr01.mtx";
    std::cout << "Reading matrix from: " << path << std::endl;

    // 读取稀疏矩阵
    Eigen::SparseMatrix<double> A = readSparseMatrix(path);

    // 随机生成同阶向量b
    int n = A.rows();
    Eigen::VectorXd b = Eigen::VectorXd::Random(n);

    // 初始猜测解x
    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);

    // 计时开始
    clock_t start_time = clock();

    // 使用高斯-赛德尔迭代法求解Ax = b
    int iter = gaussSeidelIteration(A, b, x);

    // 计时结束
    clock_t end_time = clock();
    double duration = double(end_time - start_time) / CLOCKS_PER_SEC;

    // 输出结果
    std::cout << "Solution x: \n" << x.transpose() << std::endl;
    std::cout << "Iterations: " << iter << std::endl;
    std::cout << "Time taken: " << duration << " seconds" << std::endl;

    return 0;
}
