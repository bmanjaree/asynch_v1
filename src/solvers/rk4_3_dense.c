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

void TheRKDense4_3_b(double theta, double *b);

// Constructs a dense output RK method of order 4 at each step (order 3 for dense output)
void TheRKDense4_3(RKMethod* method)
{
    method->num_stages = 4;
    method->unique_c = 4;
    method->exp_imp = 0;
    //method->A = m_get(method->num_stages, method->num_stages);
    //method->b = v_get(method->num_stages);
    method->b_theta = malloc(method->num_stages * sizeof(double));
    method->b_theta_deriv = NULL;
    //method->c = v_get(method->num_stages);
    
    method->dense_b = &TheRKDense4_3_b;
    method->dense_bderiv = NULL;
    
    //method->e = v_get(method->num_stages);
    //method->d = v_get(method->num_stages);
    
    method->e_order = 4;
    method->e_order_ratio = 4.0 / 2.0;
    //method->e_order_ratio = 1.0/2.0;
    method->d_order = 3;
    method->d_order_ratio = 3.0 / 2.0;
    //method->d_order_ratio = 1.0/2.0;
    //	method->d_max_error = 1.0/3.0;
    method->localorder = 4;

    //Build the parameters for the method
    const double A[][4] = {
        { 0.0, 0.0, 0.0, 0.0 },
        { 0.5, 0.0, 0.0, 0.0 },
        { 0.0, 0.5, 0.0, 0.0 },
        { 0.0, 0.0, 1.0, 0.0 }
    };
    //A[1][0] = .5;
    //A[2][1] = .5;
    //A[3][2] = 1;
    method->A = A[0];

    double b[] = { 1.0 / 6.0, 2.0 / 6.0, 2.0 / 6.0, 1.0 / 6.0 };
    //b[0] = 1.0 / 6.0;
    //b[1] = 2.0 / 6.0;
    //b[2] = 2.0 / 6.0;
    //b[3] = 1.0 / 6.0;
    method->b = b;

    method->dense_b(1.0, method->b_theta);

    const double c[] = { 0.0, 0.5 , 0.5, 1.0 };
    //c[0] = 0;
    //c[1] = .5;
    //c[2] = .5;
    //c[3] = 1;
    method->c = c;

    const double e[] = { 2.0 / 3.0, 0.0, -4.0 / 3.0, 2.0 / 3.0 };
    //e[0] = 2.0 / 3.0;
    //e[1] = 0.0;
    //e[2] = -4.0 / 3.0;
    //e[3] = 2.0 / 3.0;
    method->e = e;

    const double d[] = { 2.0 / 9.0, 0.0, -4.0 / 9.0, 2.0 / 9.0 };
    //d[0] = 2.0 / 9.0;
    //d[1] = 0.0;
    //d[2] = -4.0 / 9.0;
    //d[3] = 2.0 / 9.0;
    method->d = d;

    method->w = NULL;
}

//The b(theta) coefficients for TheRKDense4_3()
void TheRKDense4_3_b(double theta, double *b)
{
    b[0] = theta - 3.0 / 2.0*theta*theta + 2.0 / 3.0*theta*theta*theta;
    b[1] = theta*theta - 2.0 / 3.0*theta*theta*theta;
    b[2] = b[1];
    b[3] = -1.0 / 2.0*theta*theta + 2.0 / 3.0*theta*theta*theta;
}
