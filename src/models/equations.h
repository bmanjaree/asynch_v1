#if !defined(ASYNCH_MODELS_EQUATIONS_H)
#define ASYNCH_MODELS_EQUATIONS_H


#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <structs.h>
#include <rkmethods.h>

//RHS


void dam_Tetis03(const double * const y_i, unsigned int dim, const double * const global_params, const double * const params, const QVSData * const qvs, int state, void* user, double *ans);



//Data Assimilation Models
void assim_river_rainfall(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void assim_river_rainfall_adjusted(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);

void TopLayerHillslope_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);

void model254(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);

void model400(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void model401(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void model401reservoir(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void model402(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void model403(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void model404(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void model405(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void model406(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);

double get_discharge_from_storage(double storage);
double get_storage_from_discharge(double discharge);
void Tetis03_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void routing_runoff1(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);
void routing_runoff2(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans);



#endif //ASYNCH_MODELS_EQUATIONS_H
