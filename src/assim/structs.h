#if !defined(ASSIM_STRUCTS_H)
#define ASSIM_STRUCTS_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <asynch_interface.h>
#include <structs.h>

#if defined(HAVE_PETSC)
#include <petsc.h>
#endif


#define ASYNCH_MAX_MODEL_LENGTH 256


typedef struct AssimData
{
    char model[ASYNCH_MAX_MODEL_LENGTH];        // The name of the model variant
    char db_filename[ASYNCH_MAX_PATH_LENGTH];   // The path tho the discharge .dbc file
    ConnData conninfo;	        // Query to get link ids with gauges, query to download gauge readings
    unsigned int num_obs;       // Number of observation sites
    unsigned int* obs_locs;     // Link index in the sys[] vector associtated with the site
    unsigned int num_steps;     // Number of time step to use for the optimization
    double obs_time_step;       // Observation time step
    unsigned int max_least_squares_iters;   // Maximum number of LS iterations
    Lookup* id_to_assim;
} AssimData;

typedef struct UpstreamData
{
    unsigned int* fit_states;       //Holds the index in each state vector of the ith sensitivity at this link.
    unsigned int* fit_to_universal; //Holds universal index of the ith sensitivity at this link.
    unsigned int num_fit_states;    //Number of sensitivity at this link
    unsigned int num_upstreams;     //Number of the upstream links
    Link** upstreams;               //List of the upstream links
    unsigned int num_parents;       //Number of the parents links
    Link** parents;                 //List of the parents links
} UpstreamData;

typedef struct AssimWorkspace
{
    double *d_full;
    double *x_start;    //Assimilated initial condition
    double t_b;         //Background time
    double *x_b;        //Background vector
    double *HM_buffer;  //Buffer of the HM Matrix
    Vec rhs;            //Right Hand Side vector
    Vec x;              //Solution of the LS
    Vec B;              //Diagonal of the B Matrix
    Vec R;              //Diagonal of the R Matrix
    Mat HM;
    Mat HTH;
    Mat HMTR;
    KSP ksp;
    Vec invupareas;
    double obs_time_step;
    unsigned int problem_dim;
    unsigned int allstates;
    unsigned int allstates_needed;
    unsigned int num_steps;
    unsigned int *obs_locs;
    unsigned int num_obs;
    unsigned int *above_gauges;
    unsigned int num_above;
    unsigned int assim_dim;
    unsigned int *vareq_shift, *inv_vareq_shift;

    PetscInt *HM_col_indices;   //For inserting HM values
    PetscInt *d_indices;        //For inserting d values
} AssimWorkspace;

#endif //ASSIM_STRUCTS_H
