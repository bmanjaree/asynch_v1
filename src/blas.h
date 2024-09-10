#ifndef MATHMETHODS_H
#define MATHMETHODS_H

#if _MSC_VER > 1000
#pragma once
#define restrict __restrict
#endif // _MSC_VER > 1000

////BLAS
//void v_add(VEC u,VEC v,VEC w,unsigned int start);
//void m_add(VEC2 A,VEC2 B,VEC2 C);
//void mv_mlt(VEC2 A,VEC u,VEC v);
//void mTv_mlt(VEC2 A,VEC u,VEC v);
//void mm_mlt(VEC2 A,VEC2 B,VEC2 C);
//void mTm_mlt(VEC2 A,VEC2 B,VEC2 C);
//void mmT_mlt(VEC2 A,VEC2 B,VEC2 C);




// Copies a vector, x, to a vector, y
void dcopy(const double * restrict const x, double * restrict y, unsigned int begin, unsigned int end);

// Computes w = u - v.
// \param start the index of the first entries.
void dsub(const double * restrict const u, const double * restrict const v, double * restrict w, unsigned int begin, unsigned int end);

//void sm_mlt(double alpha,VEC2 A,VEC2 B,unsigned int startrow,unsigned int startcol);
//void dipaa(double alpha,VEC2 A,VEC2 B,unsigned int startrow,unsigned int startcol);
//void diag_mm_mlt(VEC2 A,VEC2 B,VEC2 C);
//void diag_mm_mltT(VEC2 A,VEC2 B,VEC2 C);
//void diag_m_add(VEC2 A,VEC2 B,VEC2 C);
//void diag_invert(VEC2 A,VEC2 B);
//void diag_mv_mlt(VEC2 A,VEC u,VEC v);
//void mTdiag_mlt(VEC2 A,VEC B,VEC2 C);
//void mdiag_mlt(VEC2 A,VEC B,VEC2 C);
//void diagm_mlt(VEC A,VEC2 B,VEC2 C);
//
//void VDVT_mlt(VEC2 V,VEC D,VEC2 A);
//void VDinvVT_mlt(VEC2 V,VEC D,VEC2 A,VEC temp);
//void VsqrtDinvVT_mlt(VEC2 V,VEC D,VEC2 A,VEC temp);
//void VsqrtDVT_mlt(VEC2 V,VEC D,VEC2 A,VEC temp);
//double vTAv(VEC2 A,VEC v);
//
//VEC lagrange_weights(unsigned short num_stages, VEC c);
////double lagrange_bary(double theta,short unsigned int s,VEC c,VEC* Q,VEC w);
//void lagrange_bary(double theta,VEC c,VEC* Z,VEC w,VEC sum);
//double lagrange(double theta,short unsigned int s,VEC c,VEC* Q);
//double norm_inf(VEC v,VEC w,unsigned int start);


/// Computes the infinity norm of the vector v.
double nrminf(const double * restrict const v, unsigned int begin, unsigned int end);

/// Computes the infinity norm of v. v_i is divided first by w_i.
/// max(v_i / w_i)
double nrminf2(const double * restrict const v, const double * restrict const w, unsigned int begin, unsigned int end);


//#pragma omp declare simd notinbranch
//void daxpy_impl(double alpha, double *x, double *y, unsigned int begin, unsigned int end);
//
//void daxpy(double alpha, VEC x, VEC y, unsigned int begin);


void daxpy(double alpha, const double * restrict const x, double * restrict y, unsigned int begin, unsigned int end);
//void sv_mlt(double val, VEC v, unsigned int begin);

// scales a vector by a constant
void dscal(double val, double * restrict v, unsigned int begin, unsigned int end);


#endif //MATHMETHODS_H