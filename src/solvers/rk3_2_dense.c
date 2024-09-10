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

void RKDense3_2_b(double theta, double *b);
void RKDense3_2_bderiv(double theta, double *b);


//Builds a dense output RK method of order 3 at each step (order 2 for dense output)
void RKDense3_2(RKMethod* method)
{
    method->num_stages = 3;
    method->unique_c = 3;
    method->exp_imp = 0;
    //method->A = m_get(method->num_stages, method->num_stages);
    method->b = malloc(method->num_stages * sizeof(double));
    method->b_theta = malloc(method->num_stages * sizeof(double));
    method->b_theta_deriv = malloc(method->num_stages * sizeof(double));
    //method->c = v_get(method->num_stages);
    method->dense_b = &RKDense3_2_b;
    method->dense_bderiv = &RKDense3_2_bderiv;
    //method->e = v_get(method->num_stages);
    //method->d = v_get(method->num_stages);
    method->e_order = 3;
    method->e_order_ratio = 3.0 / 2.0;
    //method->e_order_ratio = 1.0/2.0;
    method->d_order = 2;
    method->d_order_ratio = 2.0 / 2.0;
    //method->d_order_ratio = 1.0/2.0;
    //	method->d_max_error = 2.0/3.0;
    method->localorder = 3;

    //Build the coefficients for the method
    const double A[][3] = {
        { 0.0, 0.0, 0.0 },
        { 0.5, 0.0, 0.0 },
        { -1.0, 2.0, 0.0 }
    };
    //A[1][0] = .5;
    //A[2][0] = -1;
    //A[2][1] = 2;
    method->A = A[0];

    method->dense_b(1.0, method->b);
    method->dense_b(1.0, method->b_theta);

    const double c[] = { 0.0, 0.5, 1.0 };
    //c[0] = 0.0;
    //c[1] = .5;
    //c[2] = 1.0;
    method->c = c;

    const double e[] = { 2.0 / 3.0, -4.0 / 3.0, 2.0 / 3.0 };
    //e[0] = 2.0 / 3.0;
    //e[1] = -4.0 / 3.0;
    //e[2] = 2.0 / 3.0;
    method->e = e;

    const double d[] = { 1.0 / 3.0, -2.0 / 3.0, 1.0 / 3.0 };
    //d[0] = 1.0 / 3.0;
    //d[1] = -2.0 / 3.0;
    //d[2] = 1.0 / 3.0;
    method->d = d;

    method->w = NULL;
}

//The b(theta) coefficients for RKDense3_2()
void RKDense3_2_b(double theta, double *b)
{
    b[0] = theta*(1 + theta*(-3.0 / 2.0 + 2.0 / 3.0*theta));
    b[1] = 2 * theta*theta*(1 - 2.0 / 3.0*theta);
    b[2] = theta*theta*(2.0 / 3.0*theta - .5);
}

//The b'(theta) coefficients for RKDense3_2()
void RKDense3_2_bderiv(double theta, double *b)
{
    b[0] = 1.0 - 3.0*theta + 2.0*theta*theta;
    b[1] = 4.0*theta - 4.0*theta*theta;
    b[2] = 2.0*theta*theta - theta;
}