#ifndef ASSIM_H
#define ASSIM_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>
#include <stdlib.h>

#include <asynch_interface.h>
#include <assim/structs.h>

int LSSolveSys(AsynchSolver* asynch, AssimWorkspace* ws, double* q);

void LSResetSys(Link* sys, unsigned int N, GlobalVars* GlobalVars, double t_0, double* backup, unsigned int problem_dim, unsigned int num_forcings, TransData* my_data);

void Print_MATRIX(double** A, unsigned int m, unsigned int n);
void Print_VECTOR(double* v, unsigned int dim);

double*** DownloadGaugeReadings(unsigned int start_time, unsigned int stop_time, unsigned int** id_to_loc, unsigned int N, unsigned int* numlinks, unsigned int** ids, unsigned int** locs, unsigned int** numsteps);
double LSComputeDistance(const double * const d, const double * const q, unsigned int size);




#endif //ASSIM_H
