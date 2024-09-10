#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <stdlib.h>
//#include <math.h>
//#if defined(HAVE_UNISTD_H)
//#include <unistd.h>
//#endif

#include <rkmethods.h>
#include <solvers/lagrange.h>

void RadauIIA3_b(double theta, double *b);


//Builds a dense output version of the 3 stage RadauIIA method with dense output
void RadauIIA3_dense(RKMethod *method)
{
    method->num_stages = 3;
    method->unique_c = 3;
    method->exp_imp = 1;
    //method->A = m_get(method->num_stages, method->num_stages);
    method->b = malloc(method->num_stages * sizeof(double));
    method->b_theta = malloc(method->num_stages * sizeof(double));
    method->b_theta_deriv = NULL;
    //method->c = v_get(method->num_stages);
    //method->e = v_get(method->num_stages + 1);
    //method->d = v_get(method->num_stages + 1);
    
    method->dense_b = &RadauIIA3_b;
    method->dense_bderiv = NULL;
    
    method->e_order = 4;
    method->e_order_ratio = 4.0 / 3.0;
    method->d_order = 3;
    method->d_order_ratio = 3.0 / 2.0;
    method->localorder = 5;

    //Build the parameters for the method
    static const double A[][3] = {
        { 0.1968154772236604, -0.06553542585019839, 0.02377097434822015 },
        { 0.3944243147390873, 0.2920734116652285, -0.04154875212599793 },
        { 0.3764030627004673, 0.5124858261884216, 0.1111111111111111 }
    };
    //A[0][0] = 0.1968154772236604;
    //A[0][1] = -0.06553542585019839;
    //A[0][2] = 0.02377097434822015;
    //A[1][0] = 0.3944243147390873;
    //A[1][1] = 0.2920734116652285;
    //A[1][2] = -0.04154875212599793;
    //A[2][0] = 0.3764030627004673;
    //A[2][1] = 0.5124858261884216;
    //A[2][2] = 0.1111111111111111;
    method->A = A[0];

    method->dense_b(1.0, method->b);
    method->dense_b(1.0, method->b_theta);

    static const double c[3] = { 0.1550510257216822, 0.6449489742783178, 1.0 };
    //c[0] = 0.1550510257216822;
    //c[1] = 0.6449489742783178;
    //c[2] = 1.0;
    method->c = c;

    static const double e[4] = { -2.762305454748599, 0.379935598252729, -0.091629609865226, 0.2748888295956774 };
    //e[0] = -2.762305454748599;
    //e[1] = 0.379935598252729;
    //e[2] = -0.091629609865226;
    //e[3] = 0.2748888295956774;
    method->e = e;

    /*
    //Two trees 0
    method->d[0] = 0.8097732865210982;
    method->d[1] = -2.701318593546288;
    method->d[2] = 1.616656477429513;
    method->d[3] = 0.2748888295956774;
    */

    //Three trees 0 (Original)
    static const double d[4] = { -0.428298294115368, 0.245039074384917, -0.091629609865226, 0.2748888295956774 };
    //d[0] = -0.428298294115368;
    //d[1] = 0.245039074384917;
    //d[2] = -0.091629609865226;
    //d[3] = 0.2748888295956774;
    method->d = d;

    /*
    //No d_0, two tree 0
    method->d[0] = 1.238071580636466;
    method->d[1] = -2.946357667931205;
    method->d[2] = 1.708286087294739;

    */

    method->w = malloc(method->num_stages * sizeof(double));
    lagrange_weights(method->c, method->num_stages, method->w);
}

//The b(theta) coefficients for RadauIIA2_dense()
void RadauIIA3_b(double theta, double *b)
{
    double t2 = theta * theta;
    double t3 = t2 * theta;

    b[0] = 1.558078204724922 * theta - 1.986947221348443 * t2 + 0.805272079323988 * t3;
    b[1] = -0.891411538058256 * theta + 3.320280554681776 * t2 - 1.916383190435099 * t3;
    b[2] = 0.3333333333333333 * theta - 1.333333333333333 * t2 + 1.111111111111111 * t3;
}
