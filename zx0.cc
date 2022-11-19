#include "zx_h.h"
class biCGstabX : protected biCGstab
// Going to be Preconditioned Stable Biconjugate Gradient Method...
{
public:
    using biCGstab::_values2out;
    using biCGstab::biCGstab;
    using biCGstab::del4values2out;
    void del4biCGstab(biCGstabX *poit)
    {
        poit = NULL;
        delete poit;
    };
    void Run()
    {
        Config();
        while (loop + 1)
        {
            Calculate();
        };
    };

private:
    using biCGstab::Calculate;
    using biCGstab::Compare;
    void Config(){
        // gird distance
    int nx = 4;
    int nt = 4;
    int ns = 2;
    int nd = 2;
    double mass = 1;
    lattice_fermi src(nx, nt, ns);
    lattice_fermi ssrc(nx, nt, ns);
    lattice_fermi dest(nx, nt, ns);
    lattice_fermi dest_1(nx, nt, ns);
    lattice_fermi src1(nx, nt, ns);
    lattice_gauge U(nx, nt, ns);
    lattice_propagator prop(nx, nt, ns);
    for (int i = 0; i < U.size; i++)
        U[i] = 1.0;
    for (int i = 0; i < src.size; i++)
        if (i == 0)
            src[i] = 1;
        else
            src[i] = 0;
    for (int i = 0; i < src1.size; i++)
        if (i == 1)
            src1[i] = 1;
        else
            src1[i] = 0;

    Dslash2(src, ssrc, U, mass, true);

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
    void CalculateX(){

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
    biCGstabX bi(1e5, 33, 1e-5);
    bi.Run();
    biCGstab::values2out values = bi._values2out;
    cout << values.diff << endl;
    bi.del4values2out(&values);
    bi.del4biCGstab(&bi);
    return 0;
};
