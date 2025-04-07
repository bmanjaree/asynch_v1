#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif

#include <rksteppers.h>
#include <models/definitions.h>
#include <models/equations.h>
#include <models/check_consistency.h>
#include <models/output_constraints.h>
#include <models/check_state.h>

//Sets the various sizes and flags for the model. This method should set the following fields:
//dim:			The number of unknowns in the differential equations (or the number of ODEs at each link).
//diff_start:		The starting index of the differential unknowns.
//str_flag:		1 if reading rainfall data from a .str file, 0 else.
//binrain_flag:		1 if reading rainfall data from binary files, 0 else.
//uses_dam:		1 if dams are compatible with the model given by model_uid.
//params_size:		The total number of parameters (including precalculations) at links with no dam.
//dam_params_size:	The total number of parameters (including precalculations) at links with a dam.
//area_idx:		The entry in params where the upstream area is stored.
//disk_params:		The number of enries in param that are read from DEM data.
//Currently, this program assumes the same number of differential equations at each link.
//UnivVars* GlobalVars:	Contains the global variables for the system.
void SetParamSizes(GlobalVars* globals, void* external) {
	unsigned short int model_uid = globals->model_uid;
	unsigned int num_global_params;

    //Set dim and start of differential variables
    switch (model_uid)
    {
    //--------------------------------------------------------------------------------------------
    // case 0:	num_global_params = 6;
    //     globals->uses_dam = 0;
    //     globals->num_params = 20;
    //     globals->dam_params_size = 0;
    //     globals->area_idx = 2;
    //     globals->areah_idx = 1;
    //     globals->num_disk_params = 12;
    //     globals->convertarea_flag = 1;
    //     globals->num_forcings = 0;
    //     globals->min_error_tolerances = 1;
    //     break;
	/******************************************************************************************************
	 * Model 100s Routing only models
	 ******************************************************************************************************/
	//--------------------------------------------------------------------------------------------
	case 100:
		num_global_params = 2;
		globals->uses_dam = 0;
		globals->num_params = 4;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 3;
		globals->convertarea_flag = 0;
		globals->num_forcings = 1;
		globals->min_error_tolerances = 1;
		break;
	//--------------------------------------------------------------------------------------------
	case 101:
		num_global_params = 1; //none
		globals->uses_dam = 0;
		globals->num_params = 6;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 5;
		globals->convertarea_flag = 0;
		globals->num_forcings = 1;
		globals->min_error_tolerances = 1;
		break;
	//--------------------------------------------------------------------------------------------
	case 102:
		num_global_params = 3;
		globals->uses_dam = 0;
		globals->num_params = 4;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 3;
		globals->convertarea_flag = 0;
		globals->num_forcings = 1;
		globals->min_error_tolerances = 1;
		break;
	//--------------------------------------------------------------------------------------------
	case 103:
		num_global_params = 1; //none
		globals->uses_dam = 0;
		globals->num_params = 7;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 6;
		globals->convertarea_flag = 0;
		globals->num_forcings = 1;
		globals->min_error_tolerances = 1;
		break;

	/******************************************************************************************************
	 * Model 200s Runoff only models
	 ******************************************************************************************************/
	//--------------------------------------------------------------------------------------------
	case 200:
		num_global_params = 8;//hu,infil,perc,surfvel,subrestime,gwrestime,meltfactor,tempthres
		globals->uses_dam = 0;
		globals->num_params = 3;//ai,li,ah
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 3;
		globals->convertarea_flag = 0;
		globals->num_forcings = 4; //precip, et, temperature,soil temperature
		globals->min_error_tolerances = 7; //as many as states:static,surface,subsurf,gw,snow,surface runoff, subsurface runoff
		break;
	//--------------------------------------------------------------------------------------------
	case 204:
		num_global_params = 1;//
		globals->uses_dam = 0;
		globals->num_params = 11;//ai,li,ah,hu,infil,perc,surfvel,subrestime,gwrestime,meltfactor,tempthres
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 11;
		globals->convertarea_flag = 0;
		globals->num_forcings = 4; //precip, et, temperature,soil temperature
		globals->min_error_tolerances = 7; //as many as states:static,surface,subsurf,gw,snow,surface runoff, subsurface runoff
		break;

	/******************************************************************************************************
	* Model 400s original runoff+routing models
	******************************************************************************************************/
    //--------------------------------------------------------------------------------------------
	case 400://tetis01
		num_global_params = 11;//v0,l1,l2,hu,infil,perc,surfvel,subrestime,gwrestime,meltfactor,tempthres
		globals->uses_dam = 0;
		globals->num_params = 6;//ai,li,ah,invtau,c1,c2
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
		globals->num_disk_params = 3;
		globals->convertarea_flag = 0;
		globals->num_forcings = 5; //precip, et, temperature,soil temperature,discharge
		globals->min_error_tolerances = 6; //as many as states:static,surface,subsurf,gw,channel,snow,
		break;
	//--------------------------------------------------------------------------------------------
	case 404://tetis01
		num_global_params = 1;//none
		globals->uses_dam = 0;
		globals->num_params = 14;
		globals->dam_params_size = 0;
		globals->area_idx = 0;
		globals->areah_idx = 2;
        globals->num_disk_params = 14; //        
		globals->convertarea_flag = 0;
		globals->num_forcings = 5; //precip, et, temperature,soil temperature,discharge
		globals->min_error_tolerances = 6; //link->dim; //as many as states:static,surface,subsurf,gw,channel,snow,
		break;
	//--------------------------------------------------------------------------------------------
	default:
		printf("Error: Invalid model_uid (%u) in SetParamSizes.\n", model_uid);
		MPI_Abort(MPI_COMM_WORLD, 1);
		//--------------------------------------------------------------------------------------------
	}
    
  //Make sure the appropriate number of global parameters are given
	if (globals->num_global_params < num_global_params) {
		printf(
				"\nError: Obtained %u parameters from .gbl file. Expected %u for model model_uid %hu.\n",
				globals->num_global_params, num_global_params, model_uid);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	if (globals->num_global_params > num_global_params)
		printf(
				"\nWarning: Obtained %u parameters from .gbl file. Expected %u for model model_uid %hu.\n",
				globals->num_global_params, num_global_params, model_uid);
}

//Sets the function to be used when writing outputs. This method should set the following field:
//output_constrains_hdf5
//output_constrains_psql
//output_constrains_rec
void SetOutputConstraints(GlobalVars* globals)
{
    unsigned short int model_uid = globals->model_uid;
    //Set dim and start of differential variables
    switch (model_uid)
    {
        default:
            globals->OutputConstrainsHdf5 = NULL;
            globals->OutputConstrainsPsql = NULL;
            globals->OutputConstrainsRec = NULL;
            break;
    }
}

//Performs some unit conversions on the data in params. This takes place immediately after reading in the DEM data,
//so these changes are available in all routines of definetype.c, if params is available. Note that dam data
//and precalculations are not available here.
//VEC* params:		Vector of parameters to convert.
//unsigned int model_uid:	The index of the model.
void ConvertParams(
    double *params,
    unsigned int model_uid,
    void* external)
{	
	// All models for now use same structure so these are redundant
	// if (model_uid == 100 || model_uid == 101){
	// 	params[1] *= 1000;	//L: km -> m
    //     params[2] *= 1e6;	//A_h: km^2 -> m^2
	// }
	//Do same conversion for all
    params[1] *= 1000;		//L_h: km -> m
    params[2] *= 1e6;		//A_h: km^2 -> m^2
}

//Sets the system of ODEs and the Runge-Kutta solver for link. This method MUST set both link->differential
//	and link->solver. The Jacobian of f (link->jacobian) may be set here, if using an
//	implicit solver.
//Link* link: 		The link at which the ODEs and Runge-Kutta solver are selected.
//unsigned int model_uid: 	The index of the model to be set.
//unsigned int exp_imp: 0 if using an explicit solver, 1 if implicit.
//unsigned int dam: 	0 if no dam is present at link, 1 if a dam is present.
void InitRoutines(
    Link* link,
    unsigned int model_uid,
    unsigned int exp_imp,
    unsigned short dam,
    void* external)
{
    //Select appropriate RK Solver for the numerical method (link->solver)
    if ((model_uid == 21 || model_uid == 22 || model_uid == 23 || model_uid == 40 || model_uid == 261 || model_uid == 262) && dam == 1)
        link->solver = &ExplicitRKIndex1SolverDam;
    else if ((model_uid == 21 || model_uid == 22 || model_uid == 23 || model_uid == 40 || model_uid == 261 || model_uid == 262) && dam == 0)
        link->solver = &ExplicitRKIndex1Solver;
    else if (exp_imp == 0)
        link->solver = &ExplicitRKSolver;
    //	else if(link->method->exp_imp == 1)
    //		link->solver = &RadauRKSolver;
    else
        printf("Warning: No solver selected for link ID %u.\n", link->ID);

    //Select the RHS function of the ODE (link->differential, link->jacobian)
    // if (model_uid == 0)
    // {
    //     link->dim = 1;
    //     link->no_ini_start = link->dim;
    //     link->diff_start = 0;

    //     link->num_dense = 1;
    //     link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
    //     link->dense_indices[0] = 0;

    //    // link->differential = &simple_river;
    //    // link->jacobian = &Jsimple;
    //     link->algebraic = NULL;
    //     link->check_state = NULL;
    //     link->check_consistency = &CheckConsistency_Nonzero_1States;
    // }
    /******************************************************************************************************
	 * Model 100s Routing only models
	 ******************************************************************************************************/
	if (model_uid == 100)
    {
        link->dim = 1;
        link->no_ini_start = 1;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &routing_100;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
	else if (model_uid == 101)
    {
        link->dim = 1;
        link->no_ini_start = 1;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &routing_101;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
	else if (model_uid == 102)
    {
        link->dim = 1;
        link->no_ini_start = 1;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &routing_102;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
	else if (model_uid == 103)
    {
        link->dim = 1;
        link->no_ini_start = 1;
        link->diff_start = 0;

        link->num_dense = 1;
        link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
        link->dense_indices[0] = 0;

        link->differential = &routing_103;
        link->algebraic = NULL;
        link->check_state = NULL;
        link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
    }
    /******************************************************************************************************
	* Model 200s Routing only models
	******************************************************************************************************/
	else if (model_uid == 200) //tetis01
	{
		link->dim = 7; //static,surface,interflow,aquifer,snow,surface runoff, subsurface runoff
		link->no_ini_start = 5; //for runoff only tanks need initial starts
		link->diff_start = 0;

		link->num_dense = 1;
		link->dense_indices = (unsigned int*) realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;

		link->differential = &model200;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	} 

	else if (model_uid == 204) //tetis01
	{
		link->dim = 7; //static,surface,interflow,aquifer,snow,surface runoff, subsurface runoff
		link->no_ini_start = 5; //for runoff only tanks need initial starts
		link->diff_start = 0;

		link->num_dense = 1;
		link->dense_indices = (unsigned int*) realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;

		link->differential = &model204;
		link->algebraic = NULL;
		link->check_state = NULL;
		link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	} 

	/******************************************************************************************************
	* Model 400s Original routing+runoff models
	******************************************************************************************************/
	else if (model_uid == 400) //tetis01
	{
		link->dim = 6;
		link->no_ini_start = link->dim;
		link->diff_start = 0;

		link->num_dense = 1;
		link->dense_indices = (unsigned int*) realloc(link->dense_indices,
				link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		//link->dense_indices[1] = 7;

		if (link->has_res) {
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		} else
			link->differential = &model400;
		    link->algebraic = NULL;
		    link->check_state = NULL;
		    link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	} 
    else if (model_uid == 404) //
	{
		link->dim = 6;
		link->no_ini_start = link->dim;
		link->diff_start = 0;

		link->num_dense = 1;
		link->dense_indices = (unsigned int*) realloc(link->dense_indices,
				link->num_dense * sizeof(unsigned int));
		link->dense_indices[0] = 0;
		//link->dense_indices[1] = 7;

		if (link->has_res) {
			link->differential = &TopLayerHillslope_Reservoirs;
			link->solver = &ForcedSolutionSolver;
		} else
			link->differential = &model404;
		    link->algebraic = NULL;
		    link->check_state = NULL;
		    link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
	}
    else
		printf("Warning: No ODE selected for link ID %u.\n", link->ID);
}

//Perform precalculations needed for the differential equation.  These should be stored in params after the DEM
//	data and after the dam data (i.e. params[disk_params] is the first precalcuation, params[params_size]
//	is the first dam datum). This method is run at each link. This method can do nothing, if no precalculations are
//	expected for the model. This assumes the same number of precalculations regardless if there is a dam or not.
//VEC* global_params:		Vector of the global parameters of the system. These are already set and are available for use.
//VEC* params:			Vector of the parameters at a link. The DEM data and dam data are already set (but may be
//					altered here). Only the entries for precalculations need to be set.
//unsigned int disk_params:	The first entry of params that should be set here.
//unsigned int params_size:	First entry of the dam data. Don't change this entry or later unless you want to modify the dam!
//unsigned int model_uid:		The index of the model.

void Precalculations(
    Link* link_i,
    double *global_params, unsigned int num_global_params,
    double *params, unsigned int num_disk_params, unsigned int num_params,
    unsigned short dam,
    unsigned int model_uid,
    void* external)
{
	/******************************************************************************************************
	 * Model 100s Routing only models
	 ******************************************************************************************************/
    if (model_uid == 100)
    {
        //Order of parameters: A_i,L_i,A_h,invtau
        //The numbering is:	0   1   2	3
        //Order of global_params: v_r,lambda_1
        //The numbering is:        0      1
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_0 = global_params[0];
        double lambda_1 = global_params[1];
        vals[3] = 60.0*v_0 / ((1.0 - lambda_1)*L_i);	// [1/min]  invtau
    }
    else if (model_uid == 101)
    {
        //Order of parameters: A_i,L_i,A_h,v_0,lambda_1,invtau
        //The numbering is:	0   1   2	3	4	5
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
		double v_0 = params[3];
        double lambda_1 = params[4];
        vals[5] = 60.0*v_0 / ((1.0 - lambda_1)*L_i);	// [1/min]  invtau
    }
    else if (model_uid == 102)
    {
        //Order of parameters: A_i,L_i,A_h,invtau
        //The numbering is:	0   1   2   3
        //Order of global_params: v_0,lambda_1,lambda_2,
        //The numbering is:        0      1        2
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_0 = global_params[0];
        double lambda_1 = global_params[1];
        double lambda_2 = global_params[2];
        vals[3] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	// [1/min]  invtau
    }
    else if (model_uid == 103)
    {
        //Order of parameters: A_i,L_i,A_h,v_0,lambda_1,lambda_2,invtau
        //The numbering is:	0   1   2   3	4	5	6
        double* vals = params;
        double A_i = params[0];
        double L_i = params[1];
        double A_h = params[2];
        double v_0 = params[3];
        double lambda_1 = params[4];
        double lambda_2 = params[5];
        vals[6] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	// [1/min]  invtau
    }

	/******************************************************************************************************
	* Model 200s Runoff only models
	******************************************************************************************************/
	else if (model_uid == 200) //tetis01 model
	{
		double* vals = params;
		double A_i = params[0]; //upstream area of the hillslope
		double L_i = params[1];	// channel lenght
		double A_h = params[2]; //area of the hillslope
		double Hu = global_params[0]; //max available storage static storage [mm]
		double infiltration = global_params[1]; //infiltration rate [mm/hr]
		double percolation = global_params[2]; //percolation rate [mm/hr]
		double alfa2 = global_params[3]; //surface velocity [m/s]
		double alfa3 = global_params[4]; //linear reserv. coef gravitational storage [days]
		double alfa4 = global_params[5]; //linear reserv. coef aquifer storage [days]
		double melt_factor = global_params[6]; // melting factor in mm/hour/degree
		double temp_thres = global_params[7]; // in celsius degrees
	} 
	else if (model_uid == 204) //spatially varying
	{
		double* vals = params;
		double A_i = params[0]; // //upstream area of the hillslope
		double L_i = params[1];	// channel lenght
		double A_h = params[2]; //area of the hillslope
		double Hu = params[3]; //max available storage static storage [mm]
		double infiltration = params[4]; //infiltration rate [mm/hr]
		double percolation = params[5]; //percolation rate [mm/hr]
		double vsurf = params[6]; //surf velocity [m/s]
		double alfa3 = params[7]; //linear reserv. coef gravitational storage [days]
		double alfa4 = params[8]; //linear reserv. coef aquifer storage [days]
		double melt_factor = params[9]; // melting factor in mm/hour/degree
		double temp_thres = params[10]; // in celsius degrees
	}
	/******************************************************************************************************
	* Model 400s Original combined models
	******************************************************************************************************/
	else if (model_uid == 400) //tetis01 model
	{
		double* vals = params;
		double A_i = params[0]; //upstream area of the hillslope
		double L_i = params[1];	// channel lenght
		double A_h = params[2]; //area of the hillslope
		double v_0 = global_params[0]; //velocity river in channels [m/s]
		double lambda_1 = global_params[1]; //power discharge in routing function
		double lambda_2 = global_params[2]; //power of area in routing function
		double Hu = global_params[3]; //max available storage static storage [mm]
		double infiltration = global_params[4]; //infiltration rate [mm/hr]
		double percolation = global_params[5]; //percolation rate [mm/hr]
		double alfa2 = global_params[6]; //surface velocity [m/s]
		double alfa3 = global_params[7]; //linear reserv. coef gravitational storage [days]
		double alfa4 = global_params[8]; //linear reserv. coef aquifer storage [days]
        double melt_factor = global_params[9]; // melting factor in mm/hour/degree
        double temp_thres = global_params[10]; // in celsius degrees
		vals[3] = 60.0 * v_0 * pow(A_i, lambda_2) / ((1.0 - lambda_1) * L_i);//[1/min]  invtau params[3]
		vals[4] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
		vals[5] = A_h / 60.0;	//  c_2

	} else if (model_uid == 404) //tetis01 model
		{
		double* vals = params;
		double A_i = params[0]; // //upstream area of the hillslope
		double L_i = params[1];	// channel lenght
		double A_h = params[2]; //area of the hillslope
		double v_0 = params[3]; //velocity river in channels [m/s]
		double lambda_1 = params[4]; //power discharge in routing function
		double lambda_2 = params[5]; //power of area in routing function
		double Hu = params[6]; //max available storage static storage [mm]
		double infiltration = params[7]; //infiltration rate [mm/hr]
		double percolation = params[8]; //percolation rate [mm/hr]
		double vsurf = params[9]; //surf velocity [m/s]
		double alfa3 = params[10]; //linear reserv. coef gravitational storage [days]
		double alfa4 = params[11]; //linear reserv. coef aquifer storage [days]
        double melt_factor = params[12]; // melting factor in mm/hour/degree
        double temp_thres = params[13]; // in celsius degrees
	}
}

//Set the initial condition for the differential equations a link. This method will be called once for each link. The differential
//	variables from any .ini, .uini, or database are already set here. Precalculations have also been made previously. All
//	algebraic variables MUST be set here. Other initial conditions may be set here as well, depending on the model.
//VEC* global_params:	The vector of global parameters for the system.
//VEC* params:		The parameters for this link. DEM, dam parameters, and precalculations are available.
//unsigned int dam:	1 if a dam is present at this link, 0 if no dam is present.
//VEC* y_0:		The initial condition vector. Store the initial data here.
//unsigned int model_uid:	The index of the model.
//Returns the state of the solution (use 0 if state discontinuities are not a concern). !!!! Should the return value be an int? !!!!
//This function is only called if the model is initialized with .ini or .uini files.
int ReadInitData(
    double *global_params, unsigned int num_global_params,
    double *params, unsigned int num_params,
    QVSData* qvs,
    unsigned short int dam,
    double *y_0, unsigned int dim,
    unsigned int model_uid,
    unsigned int diff_start, unsigned int no_init_start,
    void* user,
    void* external)
{
    unsigned int state;

    if (model_uid==100 || model_uid==101 || model_uid==102 || model_uid==103)
    {
		return 0;
	}

	else if (model_uid == 400 || model_uid == 404)        //tetis
	{

	} 
    else {
		//If not using algebraic variables, then everything is already set
		return 0;
	}

	return 0;
}


