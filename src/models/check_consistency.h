#if !defined(ASYNCH_MODEL_CHECK_CONSISTENCY_H)
#define ASYNCH_MODEL_CHECK_CONSISTENCY_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


// Shared consistency checking functions
void CheckConsistency_Nonzero_1States(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,    
    void *user);

void CheckConsistency_Nonzero_AllStates_q(
    double *y, unsigned int num_dof,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user);

#endif //!defined(ASYNCH_MODEL_CHECK_CONSISTENCY_H)

