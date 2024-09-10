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


void DOPRI5_b(double theta, double *);
void DOPRI5_bderiv(double theta, double *b);


static double sq(double x) { return x * x; }


//Builds a dense output version of Dormand-Prince 5(4) with dense order 4
void DOPRI5_dense(RKMethod* method)
{
    method->num_stages = 7;
    method->unique_c = 6;
    method->exp_imp = 0;
    //method->A = m_get(method->num_stages, method->num_stages);
    method->b = malloc(method->num_stages * sizeof(double));
    method->b_theta = malloc(method->num_stages * sizeof(double));
    method->b_theta_deriv = malloc(method->num_stages * sizeof(double));
    //method->c = v_get(method->num_stages);
    method->dense_b = &DOPRI5_b;
    method->dense_bderiv = &DOPRI5_bderiv;
    //method->e = v_get(method->num_stages);
    //method->d = v_get(method->num_stages);
    method->e_order = 5;
    //method->e_order = 6;
    method->e_order_ratio = 6.0 / 5.0;
    //method->e_order_ratio = 1.0/5.0;
    method->d_order = 4;
    //method->d_order = 5;
    method->d_order_ratio = 5.0 / 4.0;
    //method->d_order_ratio = 1.0/4.0;
    //	method->d_max_error = .6510416666666667;
    method->localorder = 5;

    //Build the parameters for the method
    static const double A[][7] = {
        { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        { 1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        { 3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        { 44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0, 0.0 },
        { 19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0, 0.0, 0.0 },
        { 9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0, 0.0, 0.0 },
        { 35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0 }
    };
    //A[1][0] = 1.0 / 5.0;
    //A[2][0] = 3.0 / 40.0;
    //A[2][1] = 9.0 / 40.0;
    //A[3][0] = 44.0 / 45.0;
    //A[3][1] = -56.0 / 15.0;
    //A[3][2] = 32.0 / 9.0;
    //A[4][0] = 19372.0 / 6561.0;
    //A[4][1] = -25360.0 / 2187.0;
    //A[4][2] = 64448.0 / 6561.0;
    //A[4][3] = -212.0 / 729.0;
    //A[5][0] = 9017.0 / 3168.0;
    //A[5][1] = -355.0 / 33.0;
    //A[5][2] = 46732.0 / 5247.0;
    //A[5][3] = 49.0 / 176.0;
    //A[5][4] = -5103.0 / 18656.0;
    //A[6][0] = 35.0 / 384.0;
    //A[6][1] = 0.0;
    //A[6][2] = 500.0 / 1113.0;
    //A[6][3] = 125.0 / 192.0;
    //A[6][4] = -2187.0 / 6784.0;
    //A[6][5] = 11.0 / 84.0;
    method->A = A[0];

    method->dense_b(1.0, method->b);
    method->dense_b(1.0, method->b_theta);
    method->dense_bderiv(1.0, method->b_theta_deriv);

    static const double c[] = { 0.0, 0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0 };
    //c[0] = 0.0;
    //c[1] = .2;
    //c[2] = .3;
    //c[3] = .8;
    //c[4] = 8.0 / 9.0;
    //c[5] = 1.0;
    //c[6] = 1.0;
    method->c = c;

    static const double e[] = { 71.0 / 57600.0, 0.0, -71.0 / 16695.0, 71.0 / 1920.0, -17253.0 / 339200.0, 22.0 / 525.0, -1.0 / 40.0 };
    //e[0] = 71.0 / 57600.0;
    //e[1] = 0.0;
    //e[2] = -71.0 / 16695.0;
    //e[3] = 71.0 / 1920.0;
    //e[4] = -17253.0 / 339200.0;
    //e[5] = 22.0 / 525.0;
    //e[6] = -1.0 / 40.0;
    method->e = e;

    static const double d[] = { .610351562499951, 0.0, -2.105795148247852, 18.310546874999346, -25.185639003536881, 20.749496981890658, -12.378961267605213 };
    //d[0] = .610351562499951;
    //d[1] = 0.0;
    //d[2] = -2.105795148247852;
    //d[3] = 18.310546874999346;
    //d[4] = -25.185639003536881;
    //d[5] = 20.749496981890658;
    //d[6] = -12.378961267605213;
    method->d = d;

    method->w = malloc(method->num_stages * sizeof(double));
    lagrange_weights(method->c, method->num_stages, method->w);
}

//The b(theta) coefficients for DOPRI5_dense()
void DOPRI5_b(double theta, double *b)
{
    /*
    double theta_sq = sq(theta);
    double ntheta2m3 = 3.0 - 2.0*theta;
    double thetam1_sq = sq(theta-1.0);

    b[0] = theta_sq*( ntheta2m3*(35.0/384.0) - thetam1_sq*(5.0/11282082432.0)*(2558722523.0-31403016.0*theta) ) + theta*thetam1_sq;
    b[1] = 0.0;
    b[2] = theta_sq*( ntheta2m3*(500.0/1113.0) + thetam1_sq*(100.0/32700410799.0)*(882725551.0-15701508.0*theta) );
    b[3] = theta_sq*( ntheta2m3*(125.0/192.0) - thetam1_sq*(25.0/1880347072.0)*(443332067.0-31403016.0*theta) );
    b[4] = theta_sq*( ntheta2m3*(-2187.0/6784.0) + thetam1_sq*(32805.0/199316789632.0)*(23143187.0-3489224.0*theta) );
    b[5] = theta_sq*( ntheta2m3*(11.0/84.0) - thetam1_sq*(55.0/822651844.0)*(29972135.0-7076736.0*theta) );
    b[6] = theta_sq*( (theta-1.0) + thetam1_sq*(10.0/29380423.0)*(7414447.0-829305.0*theta) );
    */

    b[0] = sq(theta)*((3.0 - 2.0*theta)*(35.0 / 384.0) - sq(theta - 1.0)*(5.0 / 11282082432.0)*(2558722523.0 - 31403016.0*theta)) + theta*sq(theta - 1.0);
    b[1] = 0.0;
    b[2] = sq(theta)*((3.0 - 2.0*theta)*(500.0 / 1113.0) + sq(theta - 1.0)*(100.0 / 32700410799.0)*(882725551.0 - 15701508.0*theta));
    b[3] = sq(theta)*((3.0 - 2.0*theta)*(125.0 / 192.0) - sq(theta - 1.0)*(25.0 / 1880347072.0)*(443332067.0 - 31403016.0*theta));
    b[4] = sq(theta)*((3.0 - 2.0*theta)*(-2187.0 / 6784.0) + sq(theta - 1.0)*(32805.0 / 199316789632.0)*(23143187.0 - 3489224.0*theta));
    b[5] = sq(theta)*((3.0 - 2.0*theta)*(11.0 / 84.0) - sq(theta - 1.0)*(55.0 / 822651844.0)*(29972135.0 - 7076736.0*theta));
    b[6] = sq(theta)*((theta - 1.0) + sq(theta - 1.0)*(10.0 / 29380423.0)*(7414447.0 - 829305.0*theta));

    /*
    b[0] = sq(theta)*(3.0-2.0*theta)*35.0/384.0 + theta*sq(theta-1.0) - sq(theta)*sq(theta-1.0)*5.0*(2558722523.0-31403016.0*theta)/11282082432.0;
    b[1] = 0.0;
    b[2] = sq(theta)*(3.0-2.0*theta)*500.0/1113.0 + sq(theta)*sq(theta-1.0)*100.0*(882725551.0-15701508.0*theta)/32700410799.0;
    b[3] = sq(theta)*(3.0-2.0*theta)*125.0/192.0 - sq(theta)*sq(theta-1.0)*25.0*(443332067.0-31403016.0*theta)/1880347072.0;
    b[4] = sq(theta)*(3.0-2.0*theta)*-2187.0/6784.0 + sq(theta)*sq(theta-1.0)*32805.0*(23143187.0-3489224.0*theta)/199316789632.0;
    b[5] = sq(theta)*(3.0-2.0*theta)*11.0/84.0 - sq(theta)*sq(theta-1.0)*55.0*(29972135.0-7076736.0*theta)/822651844.0;
    b[6] = sq(theta)*(theta-1.0) + sq(theta)*sq(theta-1.0)*10.0*(7414447.0-829305.0*theta)/29380423.0;
    */
}

//The b'(theta) coefficients for DOPRI5_dense()
void DOPRI5_bderiv(double theta, double *b)
{
    const double prod1 = 6.0*theta*(1.0 - theta);
    const double prod2 = 2.0*theta*(theta - 1.0)*(2.0*theta - 1.0);
    const double prod3 = sq(theta)*sq(theta - 1.0);

    b[0] = prod1*(35.0 / 384.0) + (theta - 1.0)*(3.0*theta - 1.0) - 2.0*theta*(theta - 1.0)*(2.0*theta - 1.0)*(5.0 / 11282082432.0)*(2558722523.0 - 31403016.0*theta) + prod3*(157015080.0 / 11282082432.0);
    b[1] = 0.0;
    b[2] = prod1*(500.0 / 1113.0) + prod2*(100.0 / 32700410799.0)*(882725551.0 - 15701508.0*theta) - prod3*(1570150800.0 / 32700410799.0);
    b[3] = prod1*(125.0 / 192.0) - prod2*(25.0 / 1880347072.0)*(443332067.0 - 31403016.0*theta) + prod3*(785075400.0 / 1880347072.0);
    b[4] = prod1*(-2187.0 / 6784.0) + prod2*(32805.0 / 199316789632.0)*(23143187.0 - 3489224.0*theta) - prod3*(1144640195640.0 / 199316789632.0);
    b[5] = prod1*(11.0 / 84.0) - prod2*(55.0 / 822651844.0)*(29972135.0 - 7076736.0*theta) + prod3*(389220480.0 / 822651844.0);
    b[6] = theta*(3.0*theta - 2.0) + prod2*(10.0 / 29380423.0)*(7414447.0 - 829305.0*theta) - prod3*(8293050.0 / 29380423.0);
}
