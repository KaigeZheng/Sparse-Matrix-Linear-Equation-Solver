#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

int main() {
    // 定义矩阵 A 和向量 b
    Matrix3d A;
    Vector3d b, x;

    // 初始化 A 矩阵
    A << 3, -2,  5,
         2,  5, -1,
         1,  4,  3;

    // 初始化 b 向量
    b << 1,
         6,
         3;

    // 求解 Ax = b
    // x = A.colPivHouseholderQr().solve(b);
    x = A.fullPivLu().solve(b);


    // 输出结果
    cout << "解 x = " << endl << x << endl;

    return 0;
}
