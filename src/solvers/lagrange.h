#if !defined(ASYNCH_SOLVER_LAGRANGE_H)
#define ASYNCH_SOLVER_LAGRANGE_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

//Calculate the weights for the Lagrange interpolation polynomial
void lagrange_weights(const double * const c, unsigned short num_stages, double *w);

#endif //!defined(ASYNCH_SOLVER_LAGRANGE_H)