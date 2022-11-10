#include <iostream>
#include <Eigen/Core>
using namespace std;
using namespace Eigen;

class biCGstab
{
public:
    biCGstab(int max4loop_ = 1e2,
             int nu4dim_ = 44,
             float min4diff_ = 1e-4,
             MatrixXd *poit4A_ = NULL,
             VectorXd *poit4b_ = NULL)
    {
        max4loop = max4loop_;
        nu4dim = nu4dim_;
        min4diff = min4diff_;
        poit4A = poit4A_;
        poit4b = poit4b_;
    };
    struct values2out
    {
        float diff;
        int loop;
        VectorXd xi;
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

private:
    int max4loop;
    int nu4dim;
    float min4diff;
    void *poit4A;
    void *poit4b;
    MatrixXd A;
    VectorXd b;
    VectorXd r0;
    VectorXd pi;
    VectorXd xi;
    VectorXd ri;
    VectorXd vi;
    double alpha;
    double proi;
    double wi;
    double tmp4double;
    float diff;
    int loop;
    void Config()
    {
        loop = 0;
        proi = 1.0;
        alpha = 1.0;
        wi = 1.0;
        A = MatrixXd::Random(nu4dim, nu4dim);
        b = VectorXd::Random(nu4dim);
        xi = VectorXd::Random(nu4dim);
        ri = b - A * xi;
        pi = VectorXd::Zero(nu4dim);
        r0 = ri;
        vi = VectorXd::Zero(nu4dim);
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
    {
        loop += 1;
        double &proI = tmp4double;
        proI = proi;

        double &proi1 = proi;
        proi1 = r0.adjoint() * ri;

        double &beta = tmp4double;
        beta = (proi1 / proI) * (alpha / wi);

        VectorXd &pi1 = pi;
        pi1 = ri + beta * (pi - wi * vi);

        VectorXd &vi1 = vi;
        vi1 = A * pi1;

        alpha = proi1 / (r0.adjoint() * vi1);
        VectorXd &s = ri;
        s = ri - alpha * vi1;

        double &wi1 = wi;
        tmp4double = (A * s).adjoint() * s;
        wi1 = tmp4double / ((A * s).adjoint() * (A * s));

        VectorXd &xi1 = xi;
        xi1 = xi + alpha * pi1 + wi1 * s;
        Compare();
        VectorXd &ri1 = ri;
        ri1 = s - wi1 * (A * s);
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
    biCGstab bi(1e5, 10, 1e-10);
    bi.Run();
    biCGstab::values2out values  = bi._values2out;
    cout << values.diff << endl;
    bi.del4values2out(&values);
    bi.del4biCGstab(&bi);
}
