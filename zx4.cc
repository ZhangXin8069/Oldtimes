#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace std;

int main()
{
    Eigen::MatrixXi m(2, 2);
    m << 1, 2, 3, 4;
    cout << m << endl;
    return 0;
}
