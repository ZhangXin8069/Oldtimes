#include <iostream>
#include <complex>
#include "Eigen/Core"
using namespace std;
using namespace Eigen;

class biCGstab
// Stable Biconjugate Gradient Method
{
public:
    biCGstab(int max4loop = 1e2,
             int nu4dim = 44,
             float min4diff = 1e-4,
             MatrixXcd *poit4A = NULL,
             VectorXcd *poit4b = NULL)
    {
        this->max4loop = max4loop;
        this->nu4dim = nu4dim;
        this->min4diff = min4diff;
        this->poit4A = poit4A;
        this->poit4b = poit4b;
    };
    struct values2out
    {
        float diff;
        int loop;
        VectorXcd xi;
    };
    values2out _values2out;
    void Run()
    {
        Config();
        while (loop + 1)
        {
            Calculate();
        };
    };
    void del4biCGstab(biCGstab *poit)
    {
        poit = NULL;
        delete poit;
    };
    void del4values2out(values2out *poit)
    {
        poit = NULL;
        delete poit;
    };

protected:
    int max4loop;
    int nu4dim;
    float min4diff;
    void *poit4A;
    void *poit4b;
    MatrixXcd A;
    VectorXcd b;
    VectorXcd r0;
    VectorXcd pi;
    VectorXcd xi;
    VectorXcd ri;
    VectorXcd vi;
    double alpha;
    double proi;
    double wi;
    double tmp4double;
    float diff;
    int loop;
    complex<double> tmp4complex;
    void Config()
    {
        loop = 0;
        proi = 1.0;
        alpha = 1.0;
        wi = 1.0;
        A = MatrixXcd::Random(nu4dim, nu4dim);
        b = VectorXcd::Random(nu4dim);
        xi = VectorXcd::Random(nu4dim);
        ri = b - A * xi;
        pi = VectorXcd::Zero(nu4dim);
        r0 = ri;
        vi = VectorXcd::Zero(nu4dim);
        cout << "A:" << endl
             << A << "\n*************************************\n"
             << "b:" << endl
             << b
             << "\n*************************************\n"
            //  << "x(0):"
            //  << endl
            //  << xi << endl
            //  << "\n*************************************\n"
            //  << "r(0):"
            //  << endl
            //  << ri << endl
            //  << "\n*************************************\n"
            ;
    };
    void Compare()
    {
        if (loop == max4loop || (A * xi - b).norm() < min4diff)
        {
            cout
                << "R(" << loop << "):" << endl
                << ri << "\n*************************************\n"
                << "X(" << loop << "):" << endl
                << xi << "\n*************************************\n";
            _values2out.loop = loop;
            _values2out.xi = xi;
            _values2out.diff = (A * xi - b).norm();
            loop = -1;
        }
    };
    void Calculate()
    // Refer to https://zh.m.wikipedia.org/wiki/稳定双共轭梯度法.
    {
        loop += 1;
        double &proI = tmp4double;
        proI = proi;

        double &proi1 = proi;
        tmp4complex= r0.transpose() * ri;
        cout <<tmp4complex;
        exit(0);
        // double &beta = tmp4double;
        // beta = (proi1 / proI) * (alpha / wi);

        // VectorXcd &pi1 = pi;
        // pi1 = ri + beta * (pi - wi * vi);

        // VectorXcd &vi1 = vi;
        // vi1 = A * pi1; // dslash,A is special.

        // alpha = proi1 / (r0.adjoint() * vi1);
        // VectorXcd &s = ri;
        // s = ri - alpha * vi1;

        // double &wi1 = wi;
        // tmp4double = (A * s).adjoint() * s;
        // wi1 = tmp4double / ((A * s).adjoint() * (A * s));

        // VectorXcd &xi1 = xi;
        // xi1 = xi + alpha * pi1 + wi1 * s;
        // Compare();
        // VectorXcd &ri1 = ri;
        // ri1 = s - wi1 * (A * s);

        // cout << "w(" << loop << "):" << endl
        //      << wi << "\n*************************************\n"
        //      << "a(" << loop << "):" << endl
        //      << alpha << "\n*************************************\n"
        //      << "r(" << loop << "):" << endl
        //      << ri << "\n*************************************\n"
        //      << "x(" << loop << "):" << endl
        //      << xi << "\n*************************************\n";
    };
};
int main()
{
    biCGstab bi;
    bi.Run();
    biCGstab::values2out values = bi._values2out;
    cout << values.diff << endl;
    bi.del4values2out(&values);
    bi.del4biCGstab(&bi);
    return 0;
};
