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

//void daxpy(double alpha, VEC x, VEC y, unsigned int begin)
//{
//    assert(begin < y.dim);
//    assert(x.dim == y.dim);
//    daxpy_impl(alpha, x.storage, y.storage, begin, x.dim);
//}

//void daxpy_u(double alpha, VEC x, VEC y, unsigned int begin, unsigned int end)
//{
//    assert(begin < y.dim);
//    assert(end <= y.dim);
//    daxpy_impl(alpha, x.storage, y.storage, begin, end);
//}

//void sv_mlt(double val, VEC v, unsigned int begin)
//{
//    assert(begin < v.dim);
//    unsigned int dim = v.dim;
//    for (unsigned int i = begin; i < dim; i++)
//        v.storage[i] *= val;
//}

// scales a vector by a constant
void dscal(double val, double * restrict v, unsigned int begin, unsigned int end)
{
    assert(v != NULL);
    assert(begin < end);

    for (unsigned int i = begin; i < end; i++)
        v[i] *= val;
}

////Calculates w = u + v. start is the index of the first entry of u.
//void v_add(VEC u, VEC v, VEC w, unsigned int start)
//{
//    double* a = u.storage + start;
//    double* b = v.storage;
//    double* c = w.storage;
//
//    for (unsigned int i = 0; i < w.dim; i++)
//        c[i] = a[i] + b[i];
//}
//
////Calculates C = A + B
//void m_add(MAT A, MAT B, MAT C)
//{
//    unsigned int m = C.m;
//    unsigned int n = C.n;
//
//    for (unsigned int i = 0; i < m; i++)
//        for (unsigned int j = 0; j < n; j++)
//            C.me[i][j] = A.me[i][j] + B.me[i][j];
//}
//
////Calculates v = A * u. Assumes dimensions work.
//void mv_mlt(MAT A, VEC u, VEC v)
//{
//    unsigned int m = A.m;
//    unsigned int n = A.n;
//
//    for (unsigned int i = 0; i < m; i++)
//    {
//        v.storage[i] = A.me[i][0] * u.storage[0];
//        for (unsigned int j = 1; j < n; j++)
//            v.storage[i] += A.me[i][j] * u.storage[j];
//    }
//}
//
////Calculates v = A^T * u. Assumes dimensions work.
////!!!! Cache issues? !!!!
//void mTv_mlt(MAT A, VEC u, VEC v)
//{
//    unsigned int m = A.m;
//    unsigned int n = A.n;
//
//    for (unsigned int i = 0; i < n; i++)
//    {
//        v.storage[i] = A.me[0][i] * u.storage[0];
//        for (unsigned int j = 1; j < m; j++)
//            v.storage[i] += A.me[j][i] * u.storage[j];
//    }
//}
//
////Calculates C = A * B. Assumes dimensions work.
//void mm_mlt(MAT A, MAT B, MAT C)
//{
//    unsigned int m = A.m;
//    unsigned int n = B.n;
//    unsigned int inner = A.n;
//
//    for (unsigned int i = 0; i < m; i++)
//    {
//        for (unsigned int j = 0; j < n; j++)
//            C.me[i][j] = 0.0;
//        for (unsigned int k = 0; k < inner; k++)
//        {
//            for (unsigned int j = 0; j < n; j++)
//                C.me[i][j] += A.me[i][k] * B.me[k][j];
//        }
//    }
//
//    /*
//        for(i=0;i<m;i++)
//        {
//            for(j=0;j<n;j++)
//            {
//                C.me[i][j] = 0.0;
//                for(k=0;k<inner;k++)
//                    C.me[i][j] += A.me[i][k] * B.me[k][j];
//            }
//        }
//    */
//}
//
////Calculates C = A^T * B. Assumes dimensions work.
//void mTm_mlt(MAT A, MAT B, MAT C)
//{
//    unsigned int m = A.n;
//    unsigned int n = B.n;
//    unsigned int inner = A.m;
//
//    for (unsigned int i = 0; i < m; i++)
//    {
//        for (unsigned int j = 0; j < n; j++)
//            C.me[i][j] = 0.0;
//        for (unsigned int k = 0; k < inner; k++)
//        {
//            for (unsigned int j = 0; j < n; j++)
//                C.me[i][j] += A.me[k][i] * B.me[k][j];
//        }
//    }
//
//    /*
//        for(i=0;i<m;i++)
//        {
//            for(j=0;j<n;j++)
//            {
//                C.me[i][j] = 0.0;
//                for(k=0;k<inner;k++)
//                    C.me[i][j] += A.me[k][i] * B.me[k][j];
//            }
//        }
//    */
//}
//
////Calculates C = A * B^T. Assumes dimensions work.
//void mmT_mlt(MAT A, MAT B, MAT C)
//{
//    unsigned int m = A.m;
//    unsigned int n = B.m;
//    unsigned int inner = A.n;
//
//    for (unsigned int i = 0; i < m; i++)
//    {
//        for (unsigned int j = 0; j < n; j++)
//        {
//            C.me[i][j] = 0.0;
//            for (unsigned int k = 0; k < inner; k++)
//                C.me[i][j] += A.me[i][k] * B.me[j][k];
//        }
//    }
//
//    /*
//        for(i=0;i<m;i++)
//        {
//            for(j=0;j<n;j++)
//            {
//                C.me[i][j] = 0.0;
//                for(k=0;k<inner;k++)
//                    C.me[i][j] += A.me[i][k] * B.me[j][k];
//            }
//        }
//    */
//}

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

////Calculates B = alpha * A
//void sm_mlt(double alpha, MAT A, MAT B, unsigned int startrow, unsigned int startcol)
//{
//    unsigned int m = A.m;
//    unsigned int n = A.n;
//    for (unsigned int i = 0; i < m; i++)
//        for (unsigned int j = 0; j < n; j++)
//            B.me[i + startrow][j + startcol] = alpha * A.me[i][j];
//}
//
////Calculates B = I + alpha * A
//void dipaa(double alpha, MAT A, MAT B, unsigned int startrow, unsigned int startcol)
//{
//    unsigned int m = A.m;
//    unsigned int n = A.n;
//    for (unsigned int i = 0; i < m; i++)
//    {
//        for (unsigned int j = 0; j < n; j++)
//            B.me[i + startrow][j + startcol] = alpha * A.me[i][j];
//        B.me[i + startrow][i + startcol] += 1.0;
//    }
//}
//
////Calculates A * B = C, assuming A, B, C are diagonal.
//void diag_mm_mlt(MAT A, MAT B, MAT C)
//{
//    unsigned int size = (A.m < A.n) ? A.m : A.n;
//    size = (size < B.n) ? size : B.n;
//    for (unsigned int i = 0; i < size; i++)
//        C.me[i][i] = A.me[i][i] * B.me[i][i];
//}
//
////Calculates A * B^T = C, assuming A, B, C are diagonal.
//void diag_mm_mltT(MAT A, MAT B, MAT C)
//{
//    unsigned int size = (A.m < A.n) ? A.m : A.n;
//    size = (size < B.m) ? size : B.m;
//    for (unsigned int i = 0; i < size; i++)
//        C.me[i][i] = A.me[i][i] * B.me[i][i];
//}
//
////Calculates A + B = C, assuming A, B, C are diagonal.
//void diag_m_add(MAT A, MAT B, MAT C)
//{
//    unsigned int size = (A.m < A.n) ? A.m : A.n;
//    for (unsigned int i = 0; i < size; i++)
//        C.me[i][i] = A.me[i][i] + B.me[i][i];
//}
//
////Calculates B = A^-1, assuming A is diagonal.
//void diag_invert(MAT A, MAT B)
//{
//    unsigned m = A.m;
//    for (unsigned int i = 0; i < m; i++)
//        B.me[i][i] = 1.0 / A.me[i][i];
//}
//
////Calculates v = A * u, assuming A is diagonal.
//void diag_mv_mlt(MAT A, VEC u, VEC v)
//{
//    unsigned int size = (A.m < A.n) ? A.m : A.n;
//    for (unsigned int i = 0; i < size; i++)
//        v.storage[i] = A.me[i][i] * u.storage[i];
//}
//
//
////Calculates C = A^T * B, assuming B is diagonal.
//void mTdiag_mlt(MAT A, VEC B, MAT C)
//{
//    unsigned int m = A.n;
//    unsigned int n = B.dim;
//
//    for (unsigned int  j = 0; j < n; j++)
//        for (unsigned int  i = 0; i < m; i++)
//            C.me[i][j] = A.me[j][i] * B.storage[j];
//}
//
////Calculates C = A * B, assuming B is diagonal
//void mdiag_mlt(MAT A, VEC B, MAT C)
//{
//    unsigned int m = A.m;
//    unsigned int n = B.dim;
//
//    for (unsigned int i = 0; i < m; i++)
//    {
//        for (unsigned int j = 0; j < n; j++)
//            C.me[i][j] = A.me[i][j] * B.storage[j];
//    }
//}
//
////Calculates C = A * B, assuming A is diagonal
//void diagm_mlt(VEC A, MAT B, MAT C)
//{
//    unsigned int m = A.dim;
//    unsigned int n = B.n;
//
//    for (unsigned int i = 0; i < m; i++)
//    {
//        for (unsigned int j = 0; j < n; j++)
//            C.me[i][j] = A.storage[i] * B.me[i][j];
//    }
//}
//
//
////Calculates A = V*D*V^T, where D is diagonal and square
//void VDVT_mlt(MAT V, VEC D, MAT A)
//{
//    unsigned int m = A.m;
//    unsigned int n = A.n;
//    unsigned int inner = D.dim;
//
//    for (unsigned int i = 0; i < m; i++)
//    {
//        for (unsigned int j = 0; j < n; j++)
//        {
//            A.me[i][j] = 0.0;
//            for (unsigned int k = 0; k < inner; k++)
//                A.me[i][j] += V.me[i][k] * V.me[j][k] * D.storage[k];
//        }
//    }
//}
//
////Calculates A = V*D^-1*V^T, where D is diagonal and square
//void VDinvVT_mlt(MAT V, VEC D, MAT A, VEC temp)
//{
//    unsigned int m = A.m;
//    unsigned int n = A.n;
//    unsigned int inner = D.dim;
//
//    for (unsigned int i = 0; i < m; i++)
//        temp.storage[i] = 1.0 / D.storage[i];
//
//    for (unsigned int i = 0; i < m; i++)
//    {
//        for (unsigned int j = 0; j < n; j++)
//        {
//            A.me[i][j] = 0.0;
//            for (unsigned int k = 0; k < inner; k++)
//                A.me[i][j] += V.me[i][k] * V.me[j][k] * temp.storage[k];
//        }
//    }
//}
//
////Calculates A = V*sqrt(D)^-1*V^T, where D is diagonal and square
//void VsqrtDinvVT_mlt(MAT V, VEC D, MAT A, VEC temp)
//{
//    unsigned int i, j, k;
//    unsigned int m = A.m;
//    unsigned int n = A.n;
//    unsigned int inner = D.dim;
//
//    for (i = 0; i < m; i++)
//        temp.storage[i] = 1.0 / sqrt(D.storage[i]);
//
//    for (i = 0; i < m; i++)
//    {
//        for (j = 0; j < n; j++)
//        {
//            A.me[i][j] = 0.0;
//            for (k = 0; k < inner; k++)	A.me[i][j] += V.me[i][k] * V.me[j][k] * temp.storage[k];
//        }
//    }
//}
//
////Calculates A = V*sqrt(D)*V^T, where D is diagonal and square
//void VsqrtDVT_mlt(MAT V, VEC D, MAT A, VEC temp)
//{
//    unsigned int i, j, k;
//    unsigned int m = A.m;
//    unsigned int n = A.n;
//    unsigned int inner = D.dim;
//
//    for (i = 0; i < m; i++)
//        temp.storage[i] = sqrt(D.storage[i]);
//
//    for (i = 0; i < m; i++)
//    {
//        for (j = 0; j < n; j++)
//        {
//            A.me[i][j] = 0.0;
//            for (k = 0; k < inner; k++)	A.me[i][j] += V.me[i][k] * V.me[j][k] * temp.storage[k];
//        }
//    }
//}
//
////Computes v^T*A*v
//double vTAv(MAT A, VEC v)
//{
//    unsigned int i, j;
//    double result = 0.0;
//    unsigned size = A.m;
//
//    for (i = 0; i < size; i++)
//        for (j = 0; j < size; j++)
//            result += v.storage[i] * A.me[i][j] * v.storage[j];
//
//    return result;
//}
//
////Calculate the weights for the Lagrange interpolation polynomial
////Assumes that 0 and each entry of c are node points
////This does not store the weight for 0
//VEC lagrange_weights(unsigned short num_stages, VEC c)
//{
//    unsigned int i, j;
//    VEC w = v_init(num_stages);
//
//    for (i = 0; i < num_stages; i++)
//    {
//        w.storage[i] = 1.0 / v_at(c, i);	//For the 0 node
//        for (j = 0; j < i; j++)
//            w.storage[i] *= 1.0 / (v_at(c, i) - v_at(c, j));
//        for (j = i + 1; j < num_stages; j++)
//            w.storage[i] *= 1.0 / (v_at(c, i) - v_at(c, j));
//    }
//
//    return w;
//}
//
////Evaluates the Lagrange polynomial at theta using the barycentric formula
//double lagrange_bary(double theta,short unsigned int s,VEC c,VEC* Q,VEC w)
//{
//    unsigned int i;
//    double sum = 0.0;
//    double prod = 1.0;
//
//    for(i=0;i<s;i++)
//    {
//        if( fabs(theta - c.storage[i]) < 1e-14)	return Q[i].storage[0];
//        prod *= (theta - c.storage[i]);
//        sum += w.storage[i] / (theta - c.storage[i]) * Q[i].storage[0];
//    }
//
//    return prod * sum;
//}
//*/
//
////Evaluates the Lagrange polynomial at theta using the barycentric formula.
////Assumes 0 and each component of c are the node points. At theta = 0, the poly returns 0.
////This assumes that theta will be larger than 1.
//void lagrange_bary(double theta, VEC c, VEC* Z, VEC w, VEC sum)
//{
//    unsigned short int i;
//    unsigned short int s = w.dim;
//    unsigned short int dim = sum.dim;
//    for (i = 0; i < dim; i++)	sum.storage[i] = 0.0;
//    double prod = theta;	//This is for the value at theta = 0
//
//    for (i = 0; i < s; i++)
//    {
//        prod *= (theta - c.storage[i]);
//        daxpy(w.storage[i] / (theta - c.storage[i]), Z[i], sum, 0);
//    }
//
//    for (i = 0; i < dim; i++)	sum.storage[i] *= prod;
//}
//
////Evaluates the Lagrange polynomial at theta
//double lagrange(double theta, short unsigned int s, VEC c, VEC* Q)
//{
//    unsigned int i, j;
//    double sum = 0.0;
//    double product;
//
//    for (i = 0; i < s; i++)
//    {
//        product = 1.0;
//        for (j = 0; j < i; j++)	product *= (theta - c.storage[j]) / (c.storage[i] - c.storage[j]);
//        for (j = i + 1; j < s; j++)	product *= (theta - c.storage[j]) / (c.storage[i] - c.storage[j]);
//        sum += Q[i].storage[0] * product;
//    }
//
//    return sum;
//}
//
////Computes the infinity norm of v. v_i is divided first by w_i.
//double norm_inf(VEC v, VEC w, unsigned int start)
//{
//    unsigned int i;
//    double max = fabs(v.storage[start] / w.storage[start]);
//    double val;
//    for (i = start + 1; i < v.dim; i++)
//    {
//        val = fabs(v.storage[i] / w.storage[i]);
//        max = (val > max) ? val : max;
//    }
//
//    return max;
//}

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

////Computes the infinity norm of the vector v.
//double vector_norminf(VEC v, unsigned int start)
//{
//    unsigned int i;
//    double norm = fabs(v_at(v, start));
//    double val;
//    for (i = start + 1; i < v.dim; i++)
//    {
//        val = fabs(v_at(v, i));
//        norm = (norm < val) ? val : norm;
//    }
//    return norm;
//}

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


////Prints the vector v to stdout.
//void Print_Vector(VEC v)
//{
//    unsigned int i;
//    for (i = 0; i < v.dim; i++)
//        printf("%.16f ", v.storage[i]);
//    printf("\n");
//}
//
////Prints the matrix A to stdout.
//void Print_Matrix(MAT A)
//{
//    unsigned int i, j;
//    for (i = 0; i < A.m; i++)
//    {
//        for (j = 0; j < A.n; j++)	printf("%.16f ", A.me[i][j]);
//        printf("; \n");
//    }
//}
//
////Prints the vector v in C format to stdout.
//void Print_VectorC(VEC v)
//{
//    unsigned int i;
//    printf("{ %.16f", v.storage[0]);
//    for (i = 1; i < v.dim; i++)
//        printf(", %.16f", v.storage[i]);
//    printf(" }\n");
//}
//
////Prints the matrix A in C format to stdout.
//void Print_MatrixC(MAT A)
//{
//    unsigned int i, j;
//    printf("{ ");
//    for (i = 0; i < A.m; i++)
//    {
//        printf("{ %.16f", A.me[i][0]);
//        for (j = 1; j < A.n; j++)	printf(", %.16f", A.me[i][j]);
//        printf("},\n");
//    }
//    printf("}\n");
//}
//
