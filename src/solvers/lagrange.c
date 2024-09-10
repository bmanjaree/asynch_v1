#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <stdlib.h>

#include <solvers/lagrange.h>

//Calculate the weights for the Lagrange interpolation polynomial
//Assumes that 0 and each entry of c are node points
//This does not store the weight for 0
void lagrange_weights(const double * const c, unsigned short num_stages, double *w)
{
    assert(w != NULL);

    for (unsigned int i = 0; i < num_stages; i++)
    {
        w[i] = 1.0 / c[i];	//For the 0 node
        for (unsigned int j = 0; j < i; j++)
            w[i] *= 1.0 / (c[i] - c[j]);
        for (unsigned int j = i + 1; j < num_stages; j++)
            w[i] *= 1.0 / (c[i] - c[j]);
    }
}