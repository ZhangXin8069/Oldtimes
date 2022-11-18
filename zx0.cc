#include "2D_test.h"
#include "biCGstab.h"
int main()
{
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
    CG(ssrc, dest, U, mass, 1000);
    Dslash2(src1, ssrc, U, mass, true);
    CG(ssrc, dest_1, U, mass, 1000);
    // Dslash2(dest,ssrc,U,mass,false);
    fermi_to_prop(dest, prop, 0);
    // fermi_to_prop(dest_1,prop,1);

    for (int x = 0; x < dest.lat_x; x++)
        for (int t = 0; t < dest.lat_t; t++)
            for (int s = 0; s < dest.lat_spin; s++)
            {
                //    printf("dest=%f\n",prop[(x*prop.lat_t+t)*prop.lat_spin*prop.lat_spin+(0*prop.lat_spin)+s].real());
            }

    // printf("s1=%f\n",s1.real());

    // printf("norm_propagator0=%f\n",norm_2(dest));
    printf("norm_propagator1=%f\n", norm_2(prop));
    // printf("norm_src-propagator=%.10e\n",norm_2(ssrc-src));
    // printf("dslash_1=%f\n",norm_2(ssrc));
    // printf("dslash_2=%f\n",norm_2(dest));
    biCGstab bi(1e5, 33, 1e-5);
    bi.Run();
    biCGstab::values2out values = bi._values2out;
    cout << values.diff << endl;
    bi.del4values2out(&values);
    bi.del4biCGstab(&bi);
    return 0;
};
