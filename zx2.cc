#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
int main()
{
    Matrix2cf X, Y;
    X.real() << 1, 2, 3, 4;
    X.imag() << 1, 2, 3, 4;
    Y.real() << 1, 1, 4, 5;
    Y.imag() << 1, 4, 1, 1;
    cout << "加法：" << endl
         << X + Y << endl;
    cout << "减法：" << endl
         << X - Y << endl;
    cout << "标量乘法：" << endl
         << X * 0.5 << endl;
    cout << "矩阵乘法：" << endl
         << X * Y << endl;
    cout << "共轭：" << endl
         << X.conjugate() << endl;
    cout << "共轭转置：" << endl
         << X.adjoint() << endl;
    cout << "矩阵的迹：" << endl
         << X.trace() << endl;
    cout << "矩阵的各元素之和" << endl
         << X.sum() << endl;
    ptrdiff_t i, j;
    double MaxOfM = X.real().maxCoeff(&i, &j);
    cout << "矩阵实数部分的最大系数，最小系数：" << endl
         << MaxOfM << endl
         << "位置在(" << i << "," << j << ")" << endl;
    cout << "特征值:" << X.eigenvalues() << endl;
    return 0;
}