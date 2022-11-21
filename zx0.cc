#include "zx_h.h"
class biCGstabX : protected biCGstab
// Going to be Preconditioned Stable Biconjugate Gradient Method...
{
public:
    biCGstabX(int max4loop = 1e5,
              float min4diff = 1e-5,
              int nx = 3,
              int nt = 3,
              int ns = 2,
              int nd = 2,
              double mass = 1,
              MatrixXd *poit4A = NULL,
              VectorXd *poit4b = NULL)
    {
        this->max4loop = max4loop;
        this->nu4dim = nx * nt * ns;
        this->nx = nx;
        this->nt = nt;
        this->ns = ns;
        this->nd = nd;
        this->min4diff = min4diff;
        this->poit4A = poit4A;
        this->poit4b = poit4b;
    };
    using biCGstab::_values2out;
    using biCGstab::del4values2out;
    void Run()
    {
        Config();
        while (loop + 1)
        {
            Calculate();
        };
    };
    void del4biCGstab(biCGstabX *poit)
    {
        poit = NULL;
        delete poit;
    };

private:
    int nx;
    int nt;
    int ns;
    int nd;
    double mass;
    using biCGstab::Calculate;
    using biCGstab::Compare;
    using biCGstab::Config;
    void CalculateX()
    {
        lattice_fermi src(nx, nt, ns);
        lattice_fermi dest(nx, nt, ns);
        lattice_gauge U(nx, nt, ns);
        for (int i = 0; i < U.size; i++)
        {
            U[i] = 1.0;
        };
        for (int i = 0; i < src.size; i++)
        {
            src[i] = pi(i);
        };
        Dslash2(src, dest, U, mass, true);
        for (int i = 0; i < dest.size; i++)
        {
            vi.real()(i) = dest[i].real();
        };
    };
    void Calculate()
    // Refer to https://zh.m.wikipedia.org/wiki/稳定双共轭梯度法
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
        CalculateX();
        // Can't get the value of the A matrix.
        // Still don't use plurals,but chang it later.
        //  vi1 = A * pi1;

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
        //     cout << "w(" << loop << "):" << endl
        //          << wi << "\n*************************************\n"
        //          << "a(" << loop << "):" << endl
        //          << alpha << "\n*************************************\n"
        //          << "r(" << loop << "):" << endl
        //          << ri << "\n*************************************\n"
        //          << "x(" << loop << "):" << endl
        //          << xi << "\n*************************************\n";
        //
    };
};
int main()
{
    biCGstabX bi;
    bi.Run();
    biCGstab::values2out values = bi._values2out;
    cout << values.diff << endl;
    bi.del4values2out(&values);
    bi.del4biCGstab(&bi);
    return 0;
};
