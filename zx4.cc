#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
int main()
{
    Eigen::Matrix2d mat;
    mat << 1, 2,
        3, 4;
    cout << "Here is mat.sum():       " << mat.sum() << endl;
    cout << "Here is mat.prod():      " << mat.prod() << endl;
    cout << "Here is mat.mean():      " << mat.mean() << endl;
    cout << "Here is mat.minCoeff():  " << mat.minCoeff() << endl;
    cout << "Here is mat.maxCoeff():  " << mat.maxCoeff() << endl;
    cout << "Here is mat.trace():     " << mat.trace() << endl;
    Matrix3f m = Matrix3f::Random();
    std::ptrdiff_t i, j;
    float minOfM = m.minCoeff(&i, &j);
    cout << "Here is the matrix m:\n"
         << m << endl;
    cout << "Its minimum coefficient (" << minOfM
         << ") is at position (" << i << "," << j << ")\n\n";
    RowVector4i v = RowVector4i::Random();
    int maxOfV = v.maxCoeff(&i);
    cout << "Here is the vector v: " << v << endl;
    cout << "Its maximum coefficient (" << maxOfV
         << ") is at position " << i << endl;
MatrixXf x(7,8);

}
