#ifndef ASSIM_MODELS_H
#define ASSIM_MODELS_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdlib.h>

#include <asynch_interface.h>

//#include "assim_ls_methods.h"

//int* Partition_METIS_ByEqs(Link* sys, unsigned int N, Link** leaves, unsigned int numleaves, unsigned int** my_sys, unsigned int* my_N, TransData* my_data, short int *getting);

void Setup_Errors(AsynchSolver* asynch, unsigned int problem_dim);
unsigned int BuildStateShift(AsynchSolver* asynch, unsigned int allstates, unsigned int* data_locs, unsigned int numdata, unsigned int** vareq_shift, unsigned int** inv_vareq_shift);

////Assim model 15
//void SetParamSizes_Assim(GlobalVars* GlobalVars, void* external);
//void ConvertParams_Assim(double *params, unsigned int type, void* external);
//void InitRoutines_Assim(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external);
//void InitRoutines_Model(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external);
////void Precalculations_Assim(Link* link_i, const double * const global_params, const double * const params, unsigned int disk_params, unsigned int params_size, unsigned short int dam, unsigned int type, void* external);
//void Precalculations_Assim(Link* link_i, const double * const global_params, double * const params, unsigned short has_dam, void *user);
//
//int ReadInitData_Assim(
//    const double * const global_params, unsigned int num_global_params,
//    const double * const params, unsigned int num_params,
//    double *y, unsigned int dim,
//    void *user);
//
//void assim_river_rainfall_adjusted_custom(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);

//Common routine for model 254
void SetParamSizes_Assim_254(GlobalVars* GlobalVars, void* external);
void ConvertParams_Assim_254(double *params, unsigned int type, void* external);
void InitRoutines_Assim_254(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external);
//void InitRoutines_Model_254(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external);
void Precalculations_Assim_254(Link* link_i, const double * const global_params, double * const params, unsigned short has_dam, void *user);

void CheckConsistency_Nonzero_Model254(double *y, unsigned int dim, const double * const global_params, unsigned int num_global_params, const double * const params, unsigned int num_params, void *user);
void CheckConsistency_Nonzero_Model252(double *y, unsigned int dim, const double * const global_params, unsigned int num_global_params, const double * const params, unsigned int num_params, void *user);

//Assim model 254
int ReadInitData_Assim_254(
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    double *y, unsigned int dim,
    void *user);

void TopLayerHillslope_extras_assim(
    double t,
    const double * const y_i, unsigned int dim,
    const double * const y_p, unsigned short num_parents, unsigned int max_dim,
    const double * const global_params, const double * const params,
    const double * const forcing_values, const QVSData * const qvs,
    int state, void* user, double *ans);

void Setup_Fitting_Data_Model254(AsynchSolver* asynch, unsigned int* data_locs, unsigned int numdata);

//Assim model 254, q
void InitRoutines_Assim_254_q(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external);
//void InitRoutines_Model_252(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external);

int ReadInitData_Assim_254_q(
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    double *y, unsigned int dim,
    void *user);

void TopLayerHillslope_assim_q(
    double t,
    const double * const y_i, unsigned int dim,
    const double * const y_p, unsigned short num_parents, unsigned int max_dim,
    const double * const global_params, const double * const params,
    const double * const forcing_values, const QVSData * const qvs,
    int state, void* user, double *ans);

void Setup_Fitting_Data_Model254_q(AsynchSolver* asynch, unsigned int* data_locs, unsigned int numdata);

//Assim model 254, q and s_p
void InitRoutines_Assim_254_qsp(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external);

int ReadInitData_Assim_254_qsp(
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    double *y, unsigned int dim,
    void *user);

void TopLayerHillslope_assim_qsp(
    double t,
    const double * const y_i, unsigned int dim,
    const double * const y_p, unsigned short num_parents, unsigned int max_dim,
    const double * const global_params, const double * const params,
    const double * const forcing_values, const QVSData * const qvs,
    int state, void* user, double *ans);

void Setup_Fitting_Data_Model254_qsp(AsynchSolver* asynch, unsigned int* data_locs, unsigned int numdata);

//Assim model 254, q and s_t
void InitRoutines_Assim_254_qst(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external);

int ReadInitData_Assim_254_qst(
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    double *y_0, unsigned int dim,
    void *user);

void TopLayerHillslope_assim_qst(
    double t,
    const double * const y_i, unsigned int dim,
    const double * const y_p, unsigned short num_parents, unsigned int max_dim,
    const double * const global_params, const double * const params,
    const double * const forcing_values, const QVSData * const qvs,
    int state, void* user, double *ans);

void Setup_Fitting_Data_Model254_qst(AsynchSolver* asynch, unsigned int* data_locs, unsigned int numdata);

void CheckConsistency_Nonzero_Model252_st(double *y, unsigned int dim, const double * const global_params, unsigned int num_global_params, const double * const params, unsigned int num_params, void *user);

#endif //ASSIM_MODELS_H