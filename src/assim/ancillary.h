#if !defined(ASSIM_ANCILLARY_H)
#define ASSIM_ANCILLARY_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <asynch_interface.h>
#include <assim/structs.h>


//void FindAllDischarges(double*** data, double t, unsigned int numlinks, unsigned int* numsteps, double* d);
unsigned int GaugeDownstream(const AsynchSolver* asynch, const unsigned int* obs_locs, unsigned int num_obs, unsigned int** above_gauges, bool **is_above_gauges);
int AdjustDischarges(const AsynchSolver* asynch, const unsigned int* obs_locs, const double * obs, unsigned int num_obs, unsigned int problem_dim, double* x);

void FindUpstreamLinks(const AsynchSolver* const asynch, AssimData* const assim, unsigned int problem_dim, bool trim, double obs_time_step, unsigned int num_steps, unsigned int* obs_locs, unsigned int num_obs);
void FindUpstreamLinks2(const AsynchSolver* const asynch, AssimData* const assim, unsigned int problem_dim, bool trim, double obs_time_step, unsigned int num_steps, unsigned int* obs_locs, unsigned int num_obs);

void CleanUpstreamLinks(const AsynchSolver* asynch);
void FreeUpstreamLinks(const AsynchSolver* asynch);

bool InitAssimData(AssimData* assim, const char* assim_filename);
void FreeAssimData(AssimData* assim);

int GetObservationsIds(const AsynchSolver* asynch, AssimData* assim);
int GetObservationsData(const AssimData* assim, const Lookup * const id_loc_loc, unsigned int N, unsigned int background_time_unix, double* d);

bool ReduceBadDischargeValues(Link* sys, int* assignments, unsigned int N, double* d_full, double* q, unsigned int num_steps, unsigned int* data_locs, unsigned int numdata, double* x_start, unsigned int assim_dim, double limit);

int SnapShot_ModelStates(AsynchSolver* asynch, unsigned int problem_dim);

//int GaugeDataAvailable(AssimData* Assim, unsigned int start_time, unsigned int end_time);

#endif //ASSIM_ANCILLARY_H

