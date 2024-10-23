#ifndef MATHMETHODS_H
#define MATHMETHODS_H

#if _MSC_VER > 1000
#pragma once
#define restrict __restrict
#endif // _MSC_VER > 1000


// Copies a vector, x, to a vector, y
void dcopy(const double * restrict const x, double * restrict y, unsigned int begin, unsigned int end);

// Computes w = u - v.
// \param start the index of the first entries.
void dsub(const double * restrict const u, const double * restrict const v, double * restrict w, unsigned int begin, unsigned int end);


/// Computes the infinity norm of the vector v.
double nrminf(const double * restrict const v, unsigned int begin, unsigned int end);

/// Computes the infinity norm of v. v_i is divided first by w_i.
/// max(v_i / w_i)
double nrminf2(const double * restrict const v, const double * restrict const w, unsigned int begin, unsigned int end);

void daxpy(double alpha, const double * restrict const x, double * restrict y, unsigned int begin, unsigned int end//void sv_mlt(double val, VEC v, unsigned int begin);

// scales a vector by a constant
void dscal(double val, double * restrict v, unsigned int begin, unsigned int end);


#endif //MATHMETHODS_H
