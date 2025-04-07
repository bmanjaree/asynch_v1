#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <stddef.h>

#include <models/check_consistency.h>


void CheckConsistency_Nonzero_1States(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user)
{
    assert(y != NULL);
    assert(num_dof >= 1);

    if (y[0] < 1e-14)
        y[0] = 1e-14;
}

void CheckConsistency_Nonzero_AllStates_q(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user)
{
    assert(y != NULL);
    assert(num_dof >= 1);

    if (y[0] < 1e-14) 
        y[0] = 1e-14;
    for (unsigned int i = 1; i < num_dof; i++)
        if (y[i] < 0.0)//causing issues !!!
            y[i] = 0.0;
}