#ifndef RKMETHODS_H
#define RKMETHODS_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdbool.h>

#include <structs.h>

//RK methods data
void RKDense3_2(RKMethod* method);
//void RKDense3_2_b(double theta, double *b);
//void RKDense3_2_bderiv(double theta, double *b);

void TheRKDense4_3(RKMethod* method);
//void TheRKDense4_3_b(double theta, double *b);

void DOPRI5_dense(RKMethod* method);
//void DOPRI5_b(double theta, double *b);
//void DOPRI5_bderiv(double theta, double *b);

void RadauIIA3_dense(RKMethod* method);
//void RadauIIA3_b(double theta, double *b);

#endif

