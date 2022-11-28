#include "zx_h.h"
void Dslash(VectorXcd &value4src, VectorXcd &value4dest, int nx = 3, int nt = 3, int ns = 2, float mass = 1.0)
{
     cout
         << "value4src:" << endl
         << value4src
         << "\n*************************************\n"
         << "value4dest:"
         << endl
         << value4dest << endl
         << "\n#############################\n";
     lattice_fermi src(nx, nt, ns);
     lattice_fermi dest(nx, nt, ns);
     lattice_gauge U(nx, nt, ns);
     for (int i = 0; i < U.size; i++)
     {
          cout << U.size;
          U[i] = 1.0;
     };
     for (int i = 0; i < src.size; i++)
     {
          cout << src.size;
          src[i] = value4src(i);
     };
     Dslash2(src, dest, U, mass, true);
      value4src=value4dest;
     for (int i = 0; i < dest.size; i++)
     {
          cout << dest.size;
          value4dest(i) = dest[i];
     }
     cout
         << "value4src:" << endl
         << value4src
         << "\n*************************************\n"
         << "value4dest:"
         << endl
         << value4dest << endl
         << "\n*************************************\n";
};
int main()
{
     int nu4dim;
     int nx = 4;
     int nt = 4;
     int ns = 2;
     int nd = 2;
     double mass = 1;
     nu4dim = nx * nt * ns;
     VectorXcd bi;
     VectorXcd xi;
     bi = VectorXcd::Random(nu4dim);
     xi = VectorXcd::Random(nu4dim);

     Dslash(bi, xi);
};
