#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <memory.h>

#include <blas.h>

// Copies a vector, x, to a vector, y
void dcopy(const double * restrict const x, double * restrict y, unsigned int begin, unsigned int end)
{
    assert(x != NULL);
    assert(y != NULL);
    assert(begin < end);

    for (unsigned int i = begin; i < end; i++)
        y[i] = x[i];
}


// y = a*x + y
void daxpy(double alpha, const double * restrict const x, double * restrict y, unsigned int begin, unsigned int end)
{
    assert(x != NULL);
    assert(y != NULL);
    assert(begin < end);

    for (unsigned int i = begin; i < end; i++)
        y[i] += alpha * x[i];
}


// scales a vector by a constant
void dscal(double val, double * restrict v, unsigned int begin, unsigned int end)
{
    assert(v != NULL);
    assert(begin < end);

    for (unsigned int i = begin; i < end; i++)
        v[i] *= val;
}


//Calculates w = u - v. start is the index of the first entries.
void dsub(const double * restrict const u, const double * restrict const v, double * restrict w, unsigned int begin, unsigned int end)
{
    assert(u != NULL);
    assert(v != NULL);
    assert(w != NULL);
    assert(begin < end);

    for (unsigned int i = begin; i < end; i++)
        w[i] = u[i] - v[i];
}

//Computes the infinity norm of v. v_i is divided first by w_i.
// max(v_i / w_i)
double nrminf2(const double * restrict const v, const double * restrict const w, unsigned int begin, unsigned int end)
{
    assert(v != NULL);
    assert(w != NULL);
    assert(begin < end);

    double max = fabs(v[begin] / w[begin]);
    for (unsigned int i = begin + 1u; i < end; i++)
    {
        double val = fabs(v[i] / w[i]);
        max = (val > max) ? val : max;
    }

    return max;
}


//Computes the infinity norm of the vector v.
double nrminf(const double * restrict const v, unsigned int begin, unsigned int end)
{
    double norm = fabs(v[begin]);
    for (unsigned int i = begin + 1u; i < end; i++)
    {
        double val = fabs(v[i]);
        norm = (norm < val) ? val : norm;
    }
    return norm;
}


