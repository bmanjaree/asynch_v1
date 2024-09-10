#if !defined(ASYNCH_MODEL_H)
#define ASYNCH_MODEL_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <stdio.h>
#include <stdbool.h>

#include <structs_fwd.h>


/// Get a model given its uid
/// 
/// \param_uid Model uid
/// \return A pointer to the model
AsynchModel const * GetModel(unsigned short model_uid);


/// These are the right-hand side functions for the differential equations.
///
/// \param t The current time (typically measured in minutes).
/// \param y_i The vector of the current states of the system at this link [num_dof]. Only states defined by a differential equation are available. This means the indices from diff_start and beyond are available. States defined by algebraic equations must be calculated, if needed for the differential equations.
/// \param num_dof The number of degree of freedom at this link. In other word, the dimension of the vector y_i. 
/// \param y_p The vector of the states of the system of each upstream (parent) link [num_parents][max_num_dof]. Only states defined by a differential equation are available. States defined by algebraic equations must be calculated. Further, only those states listed in dense_indices (defined in SetParamSizes) are available.
/// \param num_parents The number of upstream links (parents) to link i.
/// \param max_num_dof The maximum number of degree of freedom of every links.
/// \param global_params The vector of parameters constant in both space and time.
/// \param params The vector of parameters for link i.
/// \param forcings The array of current forcing values.
/// \param qvs The table of discharge vs storage relationships. This is only available if a dam is present at this link, and the dam_flag is set to 2.
/// \param state The current discontinuity state of the states.
/// \param user A pointer to user specified data.
/// \param ans The vector of function evaluations. Each entry of ans from diff_start (and including diff_start) should be set by this routine.
typedef void (DifferentialFunc) (
    double t,
    const double * const y_i, unsigned int num_dof,
    const double * const y_p, unsigned short num_parents, unsigned int max_num_dof,
    const double * const global_params,
    const double * const params,
    const double * const forcing_values,
    const QVSData * const qvs,
    int state,
    void *user,
    double *ans);

/// These are the right-hand side functions for the algebraic equations.
///
/// \param y The vector of current states.Only the states with index greater than or equal to *diff\_start* are available for use.
/// \param global_params The vector of parameters constant in both space and time.
/// \param params The vector of parameters for this link.
/// \param qvs The table of discharge vs storage relationships.This is only available if a dam is present at this link, and the *dam\_flag* is set to 2.
/// \param state The current discontinuity state of the states.
/// \param user A pointer to user specified data.
/// \param ans The vector of function evaluations.Each entry of *ans* from 0 to *diff\_start* (exclusive)should be set by this routine.
typedef void (AlgebraicFunc)(
    const double * const y_i, unsigned int num_dof,
    const double * const global_params,
    const double * const params,
    const QVSData * const qvs,
    int state,
    void *user,
    double *ans);

/// Jacobian of right-hand side function
typedef void (JacobianFunc)(
    double t,
    const double * const y_i, unsigned int num_dof,
    const double * const y_p, unsigned short num_parents, unsigned int max_num_dof,
    const double * const global_params,
    const double * const params,
    const double * const forcing_values,
    double *ans);

/// RK solver
typedef int (RKSolverFunc)(
    Link* link,
    GlobalVars* globals,
    int* assignments,
    bool print_flag,
    FILE* outputfile,
    ConnData* conninfo,
    Forcing* forcings,
    Workspace* workspace);

/// This routine determines in which discontinuity state the system currently is.
///
/// \param y The vector of current states.Only the states with index greater than or equal to *diff\_start* are available for use.
/// \param global_params The vector of parameters constant in both space and time.
/// \param params The vector of parameters for this link.
/// \param qvs The table of discharge vs storage relationships.This is only available if a dam is present at this link, and only if *dam\_flag* is 2.
/// \param dam The dam flag for this link.If 1, a dam is present at this link.If 0, no dam is present.
typedef int (CheckStateFunc)(
    double *y, unsigned int dim,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,    
    QVSData *qvs,
    bool has_dam,
    void *user);                                                    

/// This routine is called by the integrator to guarantee these constraints are satisfied.
///
/// \param y The vector of current states.Only the states with index greater than or equal to *diff\_start* are available for use.
/// \param params The vector of parameters for this link.
/// \param global_params The vector of parameters constant in both space and time.
typedef void (CheckConsistencyFunc)(
    double *y, unsigned int dim,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,    
    void *user);

// Models function signatures
typedef void (SetParamSizesFunc)(GlobalVars* globals, void* user);
typedef void (ConvertFunc)(double *params, unsigned int type, void* user);
typedef void (RoutinesFunc)(Link*, unsigned int, unsigned int, unsigned short has_dam, void *user);
//typedef void (PrecalculationsFunc)(Link* link_i, double *global_params, double *params, unsigned int disk_params, unsigned int params_size, unsigned short int dam, unsigned int type, void *user);
typedef void (PrecalculationsFunc)(Link* link_i, const double * const gparams, double * const lparams, unsigned short has_dam, void *user);

typedef int (InitializeEqsFunc)(
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    double *y, unsigned int dim,
    void *user);

typedef int* (PartitionFunc)(Link *sys, unsigned int N, Link **leaves, unsigned int num_leaves, Link ***my_sys, unsigned int *my_N, TransData *my_data, short int *getting);

typedef struct AsynchModel
{
    unsigned short uid;                 //!< Unique identifier of the model

    unsigned short dim;                 //!< Dimension of the problem at this link
    unsigned int diff_start;            //!< Starting index of differential variables in solution vectors
    unsigned int no_ini_start;          //!< Starting index of differential variables not read from disk

    unsigned int num_dense;             //!< Number of states where dense output is calculated (usually only discharge is used)
    unsigned int *dense_indices;        //!< List of offsets in state vector where dense output is needed

    unsigned int num_global_params;     //!< Number of global parameters

    bool uses_dam;                      //!< true if this type can use dams, false else
    unsigned int num_params;            //!< The number of params at each link without a dam
    unsigned int num_dam_params_size;   //!< The number of params at each link with a dam
    unsigned int num_disk_params;       //!< Number of parameters to read from disk
    
    unsigned int area_idx;              //!< Index of upstream area (A_i) in params
    unsigned int areah_idx;             //!< Index of hillslope area (A_h) in params
    bool convertarea_flag;              //!< true if hillslope and upstream areas are converted from km^2 to m^2, false if not
    
    unsigned int min_error_tolerances;  //!< The minimum number of error tolerances needed at every link. Used for uniform error tolerances.
    
    unsigned int num_forcings;          //!< The number of forcings

    DifferentialFunc *differential;         //!< Right-hand side function for ODE
    JacobianFunc *jacobian;                 //!< jacobian of right-hand side function
    AlgebraicFunc *algebraic;               //!< Function for algebraic variables
    CheckStateFunc *check_state;            //!< Function to check what "state" the state variables are in (for discontinuities)
    CheckConsistencyFunc *check_consistency; //!< Function to check state variables
    RKSolverFunc *solver;                   //!< RK solver to use
    
    SetParamSizesFunc *set_param_sizes;
    ConvertFunc *convert;
    RoutinesFunc *routines;
    PrecalculationsFunc *precalculations;
    InitializeEqsFunc *initialize_eqs;
    PartitionFunc *partition;
} AsynchModel;


#endif //!defined(ASYNCH_MODEL_H)
