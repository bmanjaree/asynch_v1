#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#if !defined(_MSC_VER)
#define ASYNCH_SLEEP sleep
#else
#include <windows.h>
#define ASYNCH_SLEEP Sleep
#endif

#include <assert.h>
#include <stdbool.h>
#include <math.h>

//#include <metis.h>

// Internal Asynch stuffs
#include <sort.h>
#include <rkmethods.h>
#include <rksteppers.h>
#include <models/check_consistency.h>
#include <models/check_state.h>
#include <models/equations.h>

#include <assim/structs.h>
#include <assim/ancillary.h>
#include <assim/models.h>


void Setup_Errors(AsynchSolver* asynch, unsigned int problem_dim)
{
    GlobalVars* globals = asynch->globals;
    ErrorData* errors_tol = &asynch->errors_tol;
    unsigned int max_dim = globals->max_dim;

    errors_tol->abstol = realloc(errors_tol->abstol, max_dim * sizeof(double));
    errors_tol->reltol = realloc(errors_tol->reltol, max_dim * sizeof(double));
    errors_tol->abstol_dense = realloc(errors_tol->abstol_dense, max_dim * sizeof(double));
    errors_tol->reltol_dense = realloc(errors_tol->reltol_dense, max_dim * sizeof(double));

    //Setup error
    for (unsigned int i = problem_dim + 1; i < max_dim; i++)
    {
        errors_tol->abstol[i] = errors_tol->abstol[problem_dim];
        errors_tol->reltol[i] = errors_tol->reltol[problem_dim];
        errors_tol->abstol_dense[i] = errors_tol->abstol_dense[problem_dim];
        errors_tol->reltol_dense[i] = errors_tol->reltol_dense[problem_dim];
    }
}

//Creates an array vareq_shift with allstates = sum(assim_dim).
//This is used to remove the sensitivities that are not used.
unsigned int BuildStateShift(AsynchSolver* asynch, unsigned int allstates, unsigned int* obs_locs, unsigned int num_obs, unsigned int** vareq_shift, unsigned int** inv_vareq_shift)
{
    unsigned int i, j, allstates_needed;
    Link *sys = asynch->sys;
    int *assignments = asynch->assignments;
    bool *needed = (bool*)calloc(allstates, sizeof(bool));

    for (i = 0; i < num_obs; i++)
    {
        if (assignments[obs_locs[i]] == asynch->my_rank)
        {
            UpstreamData * updata = (UpstreamData*)sys[obs_locs[i]].user;
            for (j = 0; j < updata->num_fit_states; j++)
                needed[updata->fit_to_universal[j]] = true;
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, needed, allstates, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

    *vareq_shift = (unsigned int*)malloc(allstates * sizeof(unsigned int));
    *inv_vareq_shift = (unsigned int*)calloc(allstates, sizeof(unsigned int));

    //Compute the shifts
    unsigned int shift = 0;
    for (i = 0; i < allstates; i++)
    {
        if (needed[i] == false)
            shift++;
        else
        {
            (*vareq_shift)[i] = i - shift;
            (*inv_vareq_shift)[i - shift] = i;
        }
    }

    /*
    printf("Shift is (%u, %u)\n",allstates,allstates-shift);
    for(i=0;i<allstates;i++)
    printf("%u ",(*vareq_shift)[i]);
    printf("\n");
    printf("Inv shift is (%u)\n",allstates-shift);
    for(i=0;i<allstates;i++)
    printf("%u ",(*inv_vareq_shift)[i]);
    printf("\n");
    */
    allstates_needed = allstates - shift;
    *inv_vareq_shift = (unsigned int*)realloc(*inv_vareq_shift, allstates_needed * sizeof(unsigned int));
    free(needed);
    return allstates_needed;
}

//Data assimilation model (Old Model 315) ************************************************************************************


//void SetParamSizes_Assim(GlobalVars* GlobalVars, void* external)
//{
//    GlobalVars->uses_dam = 0;
//    GlobalVars->num_params = 20;
//    GlobalVars->dam_params_size = 0;
//    GlobalVars->area_idx = 2;
//    GlobalVars->areah_idx = 1;
//    GlobalVars->num_disk_params = 12;
//    GlobalVars->convertarea_flag = 0;
//    GlobalVars->num_forcings = 1;
//    GlobalVars->min_error_tolerances = 4;
//}
//
//
//void ConvertParams_Assim(double *params, unsigned int type, void* external)
//{
//    params[0] *= 1000;	//L: km -> m
//    params[3] *= .001;	//h_b: mm -> m
//    params[4] *= .001;	//h_H: mm -> m
//}
//
//void InitRoutines_Assim(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external)
//{
//    UpstreamData* updata = (UpstreamData*)(link->user);
//    unsigned int i, problem_dim = 2;	//Number of model eqs
//
//    link->dim = problem_dim + problem_dim + (problem_dim - 1)*(problem_dim - 1) //Model eqs + variational eqs from this link
//        + updata->num_upstreams * problem_dim;	                                //Variational eqs from upstreams !!!! Too high? !!!!
//    //for(i=0;i<link->num_parents;i++)
//    //	link->dim += updata->num_upstreams[i] * problem_dim;	//Variational eqs from upstreams !!!! Too high? !!!!
//
//    if (link->dim != updata->dim)
//        printf("[%i]: Warning: calculated the number of equations at link %u twice and got %u and %u.\n", my_rank, link->ID, link->dim, updata->dim);
//    link->no_ini_start = 2;
//    link->diff_start = 0;
//
//    link->num_dense = link->dim - 1;	//Take out s_p
//    link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
//    link->dense_indices[0] = 0;
//    for (i = 1; i < link->num_dense; i++)	link->dense_indices[i] = i + 1;
//
//    link->differential = &assim_river_rainfall_adjusted_custom;
//    link->algebraic = NULL;
//    link->check_state = NULL;
//    link->check_consistency = &CheckConsistency_Nonzero_2States;
//    link->solver = &ExplicitRKSolver;
//}
//
//
//void InitRoutines_Model(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external)
//{
//    UpstreamData* updata = (UpstreamData*)(link->user);
//    unsigned int problem_dim = 2;	//Number of model eqs
//
//    link->dim = 2;
//    link->no_ini_start = 2;
//    link->diff_start = 0;
//
//    link->num_dense = 1;
//    link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
//    link->dense_indices[0] = 0;
//
//    link->differential = &river_rainfall_adjusted;
//    link->algebraic = NULL;
//    link->check_state = NULL;
//    link->check_consistency = &CheckConsistency_Nonzero_2States;
//    link->solver = &ExplicitRKSolver;
//}
//
//void Precalculations_Assim(Link* link_i, const double * const global_params, double * const params, unsigned short has_dam, void *user)
//{
//    //Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//    //The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//    //Order of global_params: v_r,lambda_1,lambda_2,v_h,A_r,RC
//    //The numbering is:        0      1        2     3   4   5
//    //Need to set entries 12-19 of params.
//    double K_T = 1.0;
//    double s_r = 1.0;
//    double rootS_h = pow(params[7], .5);
//    double L = params[0];
//    double A_h = params[1] * 1e6;	//Put into m^2
//    double eta = params[8];
//    double v_r = global_params[0];
//    double lambda_1 = global_params[1];
//    double lambda_2 = global_params[2];
//    double v_h = global_params[3];
//    double A_r = global_params[4];
//    double RC = global_params[5];
//
//    //!!!! Clean this model. You don't really need 20 parameters... !!!!
//    params[12] = 60.0*v_r*pow(params[2] / A_r, lambda_2) / ((1.0 - lambda_1)*params[0]);	//invtau [1/min]
//    params[13] = params[3] / s_r; //epsilon
//    params[14] = v_h*L;	//c_1 [m^2/s]
//    params[15] = params[6] * params[0] * params[3] / 3600.0; //c_2
//    params[16] = (1e-3 / 60.0) * RC;	//c_3
//    params[17] = 60.0*v_h*L / A_h;	//c_4 [1/min], A_h converted above
//    params[18] = K_T / 60.0;
//    params[19] = params[6] / (60.0*s_r);
//
//    //iparams[0] = link_i->location; //!!!! Is this even needed anywhere? !!!!
//}
//
//int ReadInitData_Assim(
//    const double * const global_params, unsigned int num_global_params,
//    const double * const params, unsigned int num_params,
//    double *y_0, unsigned int dim,
//    void *user)
//{
//    //For this type, all initial conditions for variational equation must be set here.
//    //Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
//    //The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
//    //Order of global_params: v_r,lambda_1,lambda_2,Q_r,A_r,RC
//    //The numbering is:        0      1        2     3   4   5
//    unsigned int i;
//    unsigned int offset = 2;
//
//    y_0[offset] = 1.0;
//    y_0[offset + 1] = 1.0;
//    y_0[offset + 2] = 0.0;
//    for (i = offset + 3; i < dim; i++)
//        y_0[i] = 0.0;
//
//    return 0;
//}

////Function for simple river system with data assimilation.
////Calculates the flow using simple parameters, using only the flow q.
////Order of parameters: L_i,A_h,A_i,h_b,h_H,max_inf_rate,K_sat,S_h,eta,b_H,c_H,d_H,invtau,epsilon,c_1,c_2,c_3,c_4,c_5,c_6
////The numbering is:     0   1   2   3   4       5         6    7   8   9   10  11  12    13      14  15  16  17  18  19
////Order of global_params: v_r,lambda_1,lambda_2,v_h,A_r,RC
////The numbering is:        0      1        2     3   4   5
////This uses the units and functions from September 18, 2011 document
////y_i[0] = q, y_i[1] = s, followed by N entries for the variational equation
//void assim_river_rainfall_adjusted_custom(
//    double t,
//    const double * const y_i, unsigned int dim,
//    const double * const y_p, unsigned short num_parents,
//    const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
//{
//    unsigned int i, j;
//    unsigned int offset = 2;		//!!!! This needs to be num_dense, but without variational eqs !!!!
//    unsigned int parent_offset;
//    unsigned int problem_dim = 2;
//    unsigned int all_states = (dim - offset) / problem_dim;
//    double inflow = 0.0;
//    UpstreamData* updata = (UpstreamData*)user;
//
//    double q = y_i[0];
//    double s_p = y_i[1];
//
//    double L = params[0];
//    double invtau = params[12];
//    double c_1 = params[14];
//    double c_3 = params[16];
//    double c_4 = params[17];
//    double lambda_1 = global_params[1];
//
//    double q_to_lambda_1 = pow(q, lambda_1);
//    double q_to_lambda_1_m1 = (q > 1e-12) ? q_to_lambda_1 / q : pow(1e-12, lambda_1 - 1.0);
//    double deriv_qpl = 1.0;
//
//    double q_pl = s_p;
//
//    //Flux equation (y_i[0])
//    ans[0] = -q + c_1 * q_pl;
//    for (i = 0; i < num_parents; i++)
//        inflow += y_p[i][0];
//    ans[0] = invtau * q_to_lambda_1 * (inflow + ans[0]);
//
//    //Ponded water equation (y_i[1])
//    ans[1] = c_3 * forcing_values[0] - c_4 * q_pl;
//    //ans[1] = c_3 * ( max(forcing_values[0] + 20.0*sin(t/5.0),0.0)) - c_4 * q_pl;
//
//    //!!!! Pull if statements out of loops (should just need two cases total) !!!!
//    //!!!! A lot of terms get repeated !!!!
//
//    //Eqs for variational equations
//    for (i = offset; i < dim; i++)	ans[i] = 0.0;
//
//    //s variable from this link
//    ans[offset] = -c_4*deriv_qpl*y_i[offset];
//
//    //q variables from this link
////	if(lambda_1 > 1e-12 && (inflow) > 1e-12)
//    ans[offset + 1] = (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i[offset + 1];
//    //	else
//    //		ans[offset + 1] = -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i[offset + 1];
//
//    //	if(lambda_1 > 1e-12 && (inflow) > 1e-12)
//    ans[offset + 2] = (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i[offset + 2] + invtau*c_1*q_to_lambda_1*deriv_qpl * y_i[offset];
//    //	else
//    //		ans[offset + 2] = -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i[offset + 2] + invtau*c_1*deriv_qpl*y_i[offset];
//
//        //Adjust offset
//    offset += 3;
//
//    //Variables from parents
//    for (i = 0; i < num_parents; i++)
//    {
//        parent_offset = 1 + problem_dim;
//
//        //Get the number of upstreams link for that parent
//        unsigned int num_upstreams = ((UpstreamData *)updata->upstreams[i]->user)->num_upstreams;
//
//        for (j = 0; j < num_upstreams; j++)
//        {
//            ans[offset] = invtau * q_to_lambda_1 * y_p[i][parent_offset];
//            //			if(lambda_1 > 1e-12 && (inflow) > 1e-12)
//            ans[offset] += (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i[offset];
//            //			else
//            //				ans[offset] += -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i[offset];
//
//            ans[offset + 1] = invtau * q_to_lambda_1 * y_p[i][parent_offset + 1];
//            //			if(lambda_1 > 1e-12 && (inflow) > 1e-12)
//            ans[offset + 1] += (lambda_1 * invtau * q_to_lambda_1_m1 * (inflow + c_1*s_p) - (lambda_1 + 1) * invtau * q_to_lambda_1) * y_i[offset + 1];
//            //			else
//            //				ans[offset + 1] += -(lambda_1 + 1.0) * invtau * q_to_lambda_1 * y_i[offset + 1];
//
//            offset += 2;
//            parent_offset += 2;
//        }
//    }
//}



//Data assimilation model (Model 254) ************************************************************************************


void SetParamSizes_Assim_254(GlobalVars* GlobalVars, void* external)
{
    GlobalVars->uses_dam = 0;
    GlobalVars->num_params = 8;
    GlobalVars->dam_params_size = 0;
    GlobalVars->area_idx = 0;
    GlobalVars->areah_idx = 2;
    GlobalVars->num_disk_params = 3;
    GlobalVars->convertarea_flag = 0;
    GlobalVars->num_forcings = 3;
    GlobalVars->min_error_tolerances = 8;
}


void ConvertParams_Assim_254(double *params, unsigned int type, void* external)
{
    params[1] *= 1000;		//L_h: km -> m
    params[2] *= 1e6;		//A_h: km^2 -> m^2
}

void InitRoutines_Assim_254(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external)
{
    UpstreamData* updata = (UpstreamData*)(link->user);
    unsigned int i, problem_dim = 4;	//Number of model eqs

    link->dim = (problem_dim + 3) + 11	    //Model eqs + variational eqs from this link
        + updata->num_upstreams * problem_dim;	//Variational eqs from upstreams
    //for(i=0;i<link->num_parents;i++)
    //	link->dim += updata->num_upstreams[i] * problem_dim;	//Variational eqs from upstreams
    //if(link->dim != updata->dim)
    //	printf("[%i]: Warning: calculated the number of equations at link %u twice and got %u and %u.\n",my_rank,link->ID,link->dim,updata->dim);
    link->no_ini_start = 4;
    link->diff_start = 0;

    link->num_dense = link->dim - 5;	//Only q, q_b, variational eqs
    link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
    link->dense_indices[0] = 0;
    link->dense_indices[1] = 6;
    for (i = 2; i < link->num_dense; i++)	link->dense_indices[i] = i + 5;

    link->differential = &TopLayerHillslope_extras_assim;
    link->algebraic = NULL;
    link->check_state = NULL;
    link->check_consistency = &CheckConsistency_Nonzero_Model254;
    link->solver = &ExplicitRKSolver;
}

//void InitRoutines_Model_254(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external)
//{
//    link->dim = 7;
//    link->no_ini_start = 4;
//    link->diff_start = 0;
//
//    link->num_dense = 2;
//    link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
//    link->dense_indices[0] = 0;
//    link->dense_indices[1] = 6;
//
//    if (link->has_res)
//    {
//        link->differential = &TopLayerHillslope_Reservoirs;
//        link->solver = &ForcedSolutionSolver;
//    }
//    else			
//        link->differential = &TopLayerHillslope_extras;
//    link->algebraic = NULL;
//    link->check_state = NULL;
//    link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
//}

void Precalculations_Assim_254(Link* link_i, const double * const global_params, double * const params, unsigned short has_dam, void *user)
{
    //Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
    //The numbering is:	0   1   2    3     4   5   6   7 
    //Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
    //The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11
    double* vals = params;
    double A_i = params[0];
    double L_i = params[1];
    double A_h = params[2];

    double v_0 = global_params[0];
    double lambda_1 = global_params[1];
    double lambda_2 = global_params[2];
    double v_h = global_params[3];
    double k_i_factor = global_params[5];

    vals[3] = 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
    vals[4] = v_h * L_i / A_h * 60.0;	//[1/min] k_2
    vals[5] = vals[4] * k_i_factor;	//[1/min] k_i
    vals[6] = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
    vals[7] = A_h / 60.0;	//  c_2
}

int ReadInitData_Assim_254(
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    double *y_0, unsigned int dim,
    void *user)
{
    //For this type, all initial conditions for variational equation must be set here.
    unsigned int i;
    unsigned int offset = 7;

    //For this type, the extra states need to be set (4,5,6)
    y_0[4] = 0.0;
    y_0[5] = 0.0;
    y_0[6] = y_0[0];

    y_0[offset++] = 1.0;  //ds_p/ds_p0
    y_0[offset++] = 0.0;  //ds_p/ds_t0

    y_0[offset++] = 0.0;  //ds_t/ds_p0
    y_0[offset++] = 1.0;  //ds_t/ds_t0

    y_0[offset++] = 0.0;  //ds_s/ds_p0
    y_0[offset++] = 0.0;  //ds_s/ds_t0
    y_0[offset++] = 1.0;  //ds_s/ds_s0

    y_0[offset++] = 1.0;  //dq/dq_0
    y_0[offset++] = 0.0;  //dq/ds_p0
    y_0[offset++] = 0.0;  //dq/ds_t0
    y_0[offset++] = 0.0;  //dq/ds_s0

    for (i = offset; i < dim; i++)	y_0[i] = 0.0;	//From upstreams

    return 0;
}

void CheckConsistency_Nonzero_Model254(double *y, unsigned int dim,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user)
{
    unsigned int i, problem_dim = 7;

    if (y[0] < 1e-14)	y[0] = 1e-14;
    for (i = 1; i < problem_dim; i++)
        if (y[i] < 0.0)	y[i] = 0.0;
}

//Function for simple river system with data assimilation.
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: A_i,L_i,A_h | invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2       3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s, followed by N entries for the variational equation
void TopLayerHillslope_extras_assim(
    double t,
    const double * const y_i, unsigned int dim,
    const double * const y_p, unsigned short num_parents, unsigned int max_dim,
    const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    UpstreamData *updata = (UpstreamData*)user;

    double lambda_1 = global_params[1];
    double k_3 = global_params[4];	//[1/min]
    double h_b = global_params[6];	//[m]
    double S_L = global_params[7];	//[m]
    double A = global_params[8];
    double B = global_params[9];
    double exponent = global_params[10];
    double v_B = global_params[11];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]

    double L = params[1];	//[m]
    double A_h = params[2];	//[m^2]
    double invtau = params[3];	//[1/min]
    double k_2 = params[4];	//[1/min]
    double k_i = params[5];	//[1/min]
    double c_1 = params[6];
    double c_2 = params[7];

    double q = y_i[0];		//[m^3/s]
    double s_p = y_i[1];	//[m]
    double s_t = y_i[2];	//[m]
    double s_s = y_i[3];	//[m]
    //double s_precip = y_i[4];	//[m]
    //double V_r = y_i[5];	//[m^3]
    double q_b = y_i[6];	//[m^3/s]

    //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    //A few calculations...
    double q_to_lambda_1 = pow(q, lambda_1);
    double q_to_lambda_1_m1 = (q > 1e-12) ? q_to_lambda_1 / q : pow(1e-12, lambda_1 - 1.0);
    double remaining = 1.0 - s_t / S_L;
    double pow_term = (remaining > 0.0) ? pow(remaining, exponent) : 0.0;
    double pow_term_m1 = (remaining > 1e-12) ? pow_term / remaining : pow(1e-12, exponent - 1.0);
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_sl = k_3 * s_s;	//[m/min]

    //Discharge
    double inflow = 0.0;
    ans[0] = -q + (q_pl + q_sl) * c_2;
    for (unsigned int i = 0, p = 0; i < num_parents; i++, p += updata->parents[i]->dim)
        inflow += y_p[p];
    ans[0] = invtau * q_to_lambda_1 * (inflow + ans[0]);

    //Hillslope
    ans[1] = forcing_values[0] * c_1 - q_pl - q_pt - e_p;
    ans[2] = q_pt - q_ts - e_t;
    ans[3] = q_ts - q_sl - e_s;

    ////Additional states
    //ans[4] = forcing_values[0] * c_1;
    //ans[5] = q_pl;
    //ans[6] = q_sl * A_h - q_b*60.0;
    //for (i = 0; i < num_parents; i++)
    //    ans[6] += y_p[i][6] * 60.0;
    ////ans[6] += k_3*y_p[i][3]*A_h;
    //ans[6] *= v_B / L;


    //!!!! Pull if statements out of loops (should just need two cases total) !!!!
    //!!!! A lot of terms get repeated !!!!

    //Init for variational equations
    unsigned int offset = 7, problem_dim = 7;
    for (unsigned int i = offset; i < dim; i++)
        ans[i] = 0.0;

    //Compute partial derivatives (local variables)
    double dfq_dq = lambda_1 * invtau * q_to_lambda_1_m1 * (-q + c_2*(k_2*s_p + k_3*s_s) + inflow) - invtau * q_to_lambda_1;
    double dfq_dsp = invtau*q_to_lambda_1*c_2*k_2;
    double dfq_dss = invtau*q_to_lambda_1*c_2*k_3;
    double dfsp_dsp = -k_2 - k_t;
    double dfsp_dst = k_2 / S_L*B*exponent*pow_term_m1*s_p;
    double dfst_dsp = k_t;
    double dfst_dst = -dfsp_dst - k_i;
    double dfss_dst = k_i;
    double dfss_dss = -k_3;

    //Compute partial derivatives (upstreams variables)
    double dfq_dupq = invtau*q_to_lambda_1;

    //Hillslope variational eqs
    //!!!! I think these are needed only if changing the appropriate state... !!!!
    ans[offset] = dfsp_dsp * y_i[offset] + dfsp_dst * y_i[offset + 2];	//s_p, s_p
    ans[offset + 1] = dfsp_dsp * y_i[offset + 1] + dfsp_dst * y_i[offset + 3];	//s_p, s_t
    ans[offset + 2] = dfst_dsp * y_i[offset] + dfst_dst * y_i[offset + 2];	//s_t, s_p
    ans[offset + 3] = dfst_dsp * y_i[offset + 1] + dfst_dst * y_i[offset + 3];	//s_t, s_t	
    ans[offset + 4] = dfss_dst * y_i[offset + 2] + dfss_dss * y_i[offset + 4];	//s_s, s_p
    ans[offset + 5] = dfss_dst * y_i[offset + 3] + dfss_dss * y_i[offset + 5];	//s_s, s_t
    ans[offset + 6] = dfss_dss * y_i[offset + 6];	//s_s, s_s

    //Discharge variational eqs from this link
    ans[offset + 7] = dfq_dq * y_i[offset + 7]; //q, q
    ans[offset + 8] = dfq_dq * y_i[offset + 8] + dfq_dsp * y_i[offset + 0] + dfq_dss * y_i[offset + 4]; //q, s_p
    ans[offset + 9] = dfq_dq * y_i[offset + 9] + dfq_dsp * y_i[offset + 1] + dfq_dss * y_i[offset + 5]; //q, s_t
    ans[offset + 10] = dfq_dq * y_i[offset + 10] + dfq_dss * y_i[offset + 6]; //q, s_s

    //Discharge variational eqs from parent links
    unsigned int current_idx = offset + 11, parent_idx;
    for (unsigned int i = 0, p = 0; i < num_parents; i++, p += updata->parents[i]->dim)
    {
        //parent_idx = offset + 8;
        parent_idx = offset + 7;

        //Get the number of upstreams link for that parent
        unsigned int num_upstreams = ((UpstreamData *)updata->upstreams[i]->user)->num_upstreams;

        for (unsigned int j = 0; j < num_upstreams; j++)
        {
            ans[current_idx] = dfq_dupq * y_p[p + parent_idx] + dfq_dq * y_i[current_idx]; //q, upq
            ans[current_idx + 1] = dfq_dupq * y_p[p + parent_idx + 1] + dfq_dq * y_i[current_idx + 1]; //q, ups_p
            ans[current_idx + 2] = dfq_dupq * y_p[p + parent_idx + 2] + dfq_dq * y_i[current_idx + 2]; //q, ups_t
            ans[current_idx + 3] = dfq_dupq * y_p[p + parent_idx + 3] + dfq_dq * y_i[current_idx + 3]; //q, ups_s
            current_idx += 4;
            parent_idx += 4;
        }
    }
}

//I think these values need to be set only for assigned links, not getting.
//!!!! This data is really only needed at the gauges. !!!!
//fit_states[i] holds the index in each state vector of the ith sensitivity at this link.
//fit_to_universal[i] holds universal index of the ith sensitivity at this link.
//These only store those sensitivites used for the fitting.
void Setup_Fitting_Data_Model254(AsynchSolver* asynch, unsigned int* obs_locs, unsigned int num_obs)
{
    unsigned int i, j, k, my_N = asynch->my_N;
    int *assignments = asynch->assignments;
    Link *sys = asynch->sys, **my_sys = asynch->my_sys, *current;
    UpstreamData *updata;

    //Number of states to fit
    unsigned int counter;
    unsigned int problem_dim = 7;	//!!!! Should be allowed to vary by link !!!!
    unsigned int assim_dim = 4;	//!!!! Should be allowed to vary by link !!!!
    //unsigned int num_change_states = 1;	//For q
    unsigned int num_change_states = 2;	//For q and s_p

    //Find links upstreams from gauges
    bool *is_above_gauges;
    unsigned int *above_gauges;
    unsigned int num_above = GaugeDownstream(asynch, obs_locs, num_obs, &above_gauges, &is_above_gauges);

    //Calculate the number of states needed for the fitting
    //unsigned int allstates_needed = num_above;	//For q
    //unsigned int allstates_needed = num_above * 2;	//For q and s_p

    //for(i=0;i<my_N;i++)
    for (i = 0; i < num_obs; i++)
    {
        if (assignments[obs_locs[i]] == asynch->my_rank)
        {
            //current = sys[my_sys[i]];
            current = &sys[obs_locs[i]];
            updata = (UpstreamData*)current->user;

            updata->num_fit_states = 1;
            for (j = 0; j < current->num_parents; j++)
            {
                //Get the number of upstreams link for that parent
                unsigned int num_upstreams = ((UpstreamData *)updata->upstreams[i]->user)->num_upstreams;
                updata->num_fit_states += num_upstreams;
            }
            updata->num_fit_states *= num_change_states;

            updata->fit_states = (unsigned int*)malloc(updata->num_fit_states * sizeof(unsigned int));
            updata->fit_to_universal = (unsigned int*)malloc(updata->num_fit_states * sizeof(unsigned int));
            /*
                        //For q
                        updata->fit_states[0] = problem_dim + 7;	//q, q here
                        updata->fit_to_universal[0] = current->location * assim_dim;
                        counter = 1;
                        for(j=0;j<current->num_parents;j++)
                        {
                            for(k=0;k<updata->num_upstreams[j];k++)
                            {
                                updata->fit_states[counter] = problem_dim + 7 + assim_dim*counter;
                                updata->fit_to_universal[counter] = updata->upstreams[j][k] * assim_dim;
                                counter++;
                            }
                        }
            */

            //For q, s_p
            updata->fit_states[0] = problem_dim + 7;	//q, q here
            updata->fit_to_universal[0] = current->location * assim_dim;
            updata->fit_states[1] = problem_dim + 8;	//q, s_p here
            updata->fit_to_universal[1] = current->location * assim_dim + 1;
            counter = 2;
            for (j = 0; j < current->num_parents; j++)
            {
                unsigned int num_upstreams = ((UpstreamData *)updata->upstreams[i]->user)->num_upstreams;
                Link **upstreams = ((UpstreamData *)updata->upstreams[i]->user)->upstreams;

                for (k = 0; k < num_upstreams; k++)
                {
                    updata->fit_states[counter] = problem_dim + 7 + assim_dim*counter / 2;
                    updata->fit_to_universal[counter] = upstreams[k]->location * assim_dim;
                    updata->fit_states[counter + 1] = problem_dim + 8 + assim_dim*counter / 2;
                    updata->fit_to_universal[counter + 1] = upstreams[k]->location * assim_dim + 1;
                    counter += 2;
                }
            }

        }
    }


    free(is_above_gauges);
    free(above_gauges);
    //return allstates_needed;
}


//Data assimilation model (Model 254 trimmed, q) ************************************************************************************


//For modifying q
void InitRoutines_Assim_254_q(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external)
{
    UpstreamData* updata = (UpstreamData*)(link->user);
    unsigned int i, problem_dim = 4;	//Number of model eqs

    //For q only
    link->dim = problem_dim + 1 	//Model eqs + variational eqs from this link
        + updata->num_upstreams;
    
    //Variational eqs from upstreams
    //for(i=0;i<link->num_parents;i++)
    //    link->dim += ((UpstreamData *)updata->parents[i]->user)->num_upstreams;
    //	link->dim += updata->num_upstreams[i];	//Variational eqs from upstreams
    
    link->no_ini_start = 4;
    link->diff_start = 0;

    link->num_dense = link->dim - 3;	//Only q, variational eqs
    link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
    link->dense_indices[0] = 0;
    for (i = 1; i < link->num_dense; i++)	link->dense_indices[i] = i + 3;

    link->differential = &TopLayerHillslope_assim_q;
    link->algebraic = NULL;
    link->check_state = NULL;
    link->check_consistency = &CheckConsistency_Nonzero_Model252;
    link->solver = &ExplicitRKSolver;
}

//void InitRoutines_Model_252(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external)
//{
//    link->dim = 4;
//    link->no_ini_start = 4;
//    link->diff_start = 0;
//
//    link->num_dense = 1;
//    link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
//    link->dense_indices[0] = 0;
//
//    if (link->has_res)
//    {
//        link->differential = &TopLayerHillslope_Reservoirs;
//        link->solver = &ForcedSolutionSolver;
//    }
//    else			link->differential = &TopLayerHillslope;
//    link->algebraic = NULL;
//    link->check_state = NULL;
//    link->check_consistency = &CheckConsistency_Nonzero_AllStates_q;
//}


int ReadInitData_Assim_254_q(
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    double *y_0, unsigned int dim,
    void *user)
{
    //UpstreamData *updata = (UpstreamData*)user;

    //double lambda_1 = global_params[1];
    //double k_3 = global_params[4];	//[1/min]
    //double h_b = global_params[6];	//[m]
    //double S_L = global_params[7];	//[m]
    //double A = global_params[8];
    //double B = global_params[9];
    //double exponent = global_params[10];

    //double L = params[1];	//[m]
    //double A_h = params[2];	//[m^2]
    //double invtau = params[3];	//[1/min]
    //double k_2 = params[4];	//[1/min]
    //double k_i = params[5];	//[1/min]
    //double c_1 = params[6];
    //double c_2 = params[7];

    //double q = y_0[0];		//[m^3/s]
    //double s_p = y_0[1];	//[m]
    //double s_t = y_0[2];	//[m]
    //double s_s = y_0[3];	//[m]

    //For this type, all initial conditions for variational equation must be set here.
    unsigned int i;
    unsigned int offset = 4;

    y_0[offset++] = 1.0;  //dq/dq_0

    ////A few calculations...
    //double q_to_lambda_1 = pow(q, lambda_1);
    //double q_to_lambda_1_m1 = (q > 1e-12) ? q_to_lambda_1 / q : pow(1e-12, lambda_1 - 1.0);

    ////Discharge
    //double inflow = 0.0;
    //for (unsigned int i = 0; i < updata->num_parents; i++)
    //    inflow += updata->parents[i]->list->head->y_approx[0];

    ////Compute partial derivatives (local variables)
    //double dfq_dq = lambda_1 * invtau * q_to_lambda_1_m1 * (-q + c_2*(k_2*s_p + k_3*s_s) + inflow) - invtau * q_to_lambda_1;

    ////Compute partial derivatives (upstreams variables)
    //double dfq_dupq = invtau*q_to_lambda_1;

    //if (updata->num_upstreams)
    //    y_0[offset++] = dfq_dupq + dfq_dq;

    for (i = offset ; i < dim ; i++)
        y_0[i] = 0.0;

    return 0;
}

void CheckConsistency_Nonzero_Model252(
    double *y, unsigned int dim,
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    void *user)
{
    unsigned int i, problem_dim = 4;

    if (y[0] < 1e-14)	y[0] = 1e-14;
    for (i = 1; i < problem_dim; i++)
        if (y[i] < 0.0)	y[i] = 0.0;
}

//Function for simple river system with data assimilation.
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: A_i,L_i,A_h | invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2       3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10     
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s, followed by N entries for the variational equation
void TopLayerHillslope_assim_q(
    double t,
    const double * const y_i, unsigned int dim,
    const double * const y_p, unsigned short num_parents, unsigned int max_dim,
    const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    UpstreamData *updata = (UpstreamData*)user;

    double lambda_1 = global_params[1];
    double k_3 = global_params[4];	//[1/min]
    double h_b = global_params[6];	//[m]
    double S_L = global_params[7];	//[m]
    double A = global_params[8];
    double B = global_params[9];
    double exponent = global_params[10];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]

    double L = params[1];	//[m]
    double A_h = params[2];	//[m^2]
    double invtau = params[3];	//[1/min]
    double k_2 = params[4];	//[1/min]
    double k_i = params[5];	//[1/min]
    double c_1 = params[6];
    double c_2 = params[7];

    double q = y_i[0];		//[m^3/s]
    double s_p = y_i[1];	//[m]
    double s_t = y_i[2];	//[m]
    double s_s = y_i[3];	//[m]

    //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    //A few calculations...
    double q_to_lambda_1 = pow(q, lambda_1);
    double q_to_lambda_1_m1 = (q > 1e-12) ? q_to_lambda_1 / q : pow(1e-12, lambda_1 - 1.0);
    double remaining = 1.0 - s_t / S_L;
    double pow_term = (remaining > 0.0) ? pow(remaining, exponent) : 0.0;
    double pow_term_m1 = (remaining > 1e-12) ? pow_term / remaining : pow(1e-12, exponent - 1.0);
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_sl = k_3 * s_s;	//[m/min]

    //Discharge
    double inflow = 0.0;
    ans[0] = -q + (q_pl + q_sl) * c_2;
    for (unsigned int i = 0, p = 0; i < num_parents; i++, p += max_dim)
        inflow += y_p[p];
    ans[0] = invtau * q_to_lambda_1 * (inflow + ans[0]);

    //Hillslope
    ans[1] = forcing_values[0] * c_1 - q_pl - q_pt - e_p;
    ans[2] = q_pt - q_ts - e_t;
    ans[3] = q_ts - q_sl - e_s;

    //Init for variational equations
    unsigned int offset = 4, problem_dim = 4;
    //for (unsigned int i = offset; i < dim; i++)
    //    ans[i] = 0.0;	//!!!! Is this needed? !!!!

    //Compute partial derivatives (local variables)
    double dfq_dq = lambda_1 * invtau * q_to_lambda_1_m1 * (-q + c_2*(k_2*s_p + k_3*s_s) + inflow) - invtau * q_to_lambda_1;
    //double dfq_dsp = invtau*q_to_lambda_1*c_2*k_2;
    //double dfq_dss = invtau*q_to_lambda_1*c_2*k_3;
    //double dfsp_dsp = -k_2 - k_t;
    //double dfsp_dst = k_2 / S_L*B*exponent*pow_term_m1*s_p;
    //double dfst_dsp = k_t;
    //double dfst_dst = -dfsp_dst - k_i;
    //double dfss_dst = k_i;
    //double dfss_dss = -k_3;

    //Compute partial derivatives (upstreams variables)
    double dfq_dupq = invtau*q_to_lambda_1;

    const unsigned int num_variational_eqs = 1;

    //Discharge variational eqs from this link
    ans[offset] = dfq_dq * y_i[offset]; //q, q

/* TEST SAM
    //Discharge variational eqs from parent links
    unsigned int current_idx = offset + 1;

    //if (updata->num_upstreams)
    //{
    //    for (unsigned int i = 0; i < num_parents; i++)
    //    {
    //        unsigned int num_upstreams = ((UpstreamData *)updata->parents[i]->user)->num_upstreams;
    //        ans[current_idx++] = dfq_dupq * y_p[i][offset] + dfq_dq * y_i[current_idx]; //q, upq

    //        for (unsigned int j = 0; j < num_upstreams; j++)
    //        {
    //            unsigned int parent_idx = j + offset + 1;

    //            assert(current_idx < ans.dim);
    //            assert(parent_idx < y_p[i].dim);
    //            assert(current_idx < y_i.dim);
    //            ans[current_idx++] = dfq_dupq * y_p[i][parent_idx] + dfq_dq * y_i[current_idx]; //q, upq
    //        }
    //    }

    //    assert(current_idx == ans.dim);
    //}


    unsigned int j = 0, p = 0;
    // For every upstream links
    for (unsigned int i = 0, j = 0; i < updata->num_upstreams; i++, j++)
    {
        unsigned int np = p + 1;

        // If switch to next parent
        if (np < updata->num_parents && updata->upstreams[i] == updata->parents[np])
        {        
            p++;
            j = 0;
        }

        unsigned int current_idx = offset + i + 1;
        unsigned int parent_idx = offset + j;

        assert(current_idx < dim);
        assert(parent_idx < dim_parents[p]);
        assert(current_idx < dim);
        ans[current_idx] = dfq_dupq * y_p[p + parent_idx] + dfq_dq * y_i[current_idx]; //q, upq
    }
*/

    //Discharge variational eqs from parent links
    unsigned int current_idx = offset + num_variational_eqs, parent_idx;
    for (unsigned int i = 0, p = 0; i < num_parents; i++, p += max_dim)
    {
        parent_idx = offset;
        unsigned int num_upstreams2 = ((UpstreamData *)updata->parents[i]->user)->num_upstreams;
        unsigned int num_upstreams = updata->parents[i]->dim - problem_dim;
        for (unsigned int j = 0; j < num_upstreams; j++)
        {
            ans[current_idx] = dfq_dupq * y_p[p + parent_idx] + dfq_dq * y_i[current_idx]; //q, upq
            current_idx += 1;
            parent_idx += 1;
        }
    }
}

//I think these values need to be set only for assigned links, not getting.
//!!!! This data is really only needed at the gauges. !!!!
//fit_states[i] holds the index in each state vector of the ith sensitivity at this link.
//fit_to_universal[i] holds universal index of the ith sensitivity at this link.
//These only store those sensitivites used for the fitting.
void Setup_Fitting_Data_Model254_q(AsynchSolver* asynch, unsigned int* obs_locs, unsigned int num_obs)
{
    unsigned int my_N = asynch->my_N;
    int *assignments = asynch->assignments;
    Link *sys = asynch->sys, **my_sys = asynch->my_sys, *current;
    UpstreamData *updata;

    //Number of states to fit
    unsigned int counter;
    unsigned int problem_dim = 4;	//!!!! Should be allowed to vary by link !!!!
    unsigned int assim_dim = 4;	//!!!! Should be allowed to vary by link !!!!
    unsigned int num_change_states = 1;	//For q
    //unsigned int num_change_states = 2;	//For q and s_p

    ////Find links upstreams from gauges
    //short int* is_above_gauges;
    //unsigned int* above_gauges;
    //unsigned int num_above = GaugeDownstream(asynch,&above_gauges,&is_above_gauges,obs_locs,num_obs);

    //Calculate the number of states needed for the fitting
    //unsigned int allstates_needed = num_above * 2;	//For q and s_p
    //unsigned int allstates_needed = num_above;	//For q

    for (unsigned int i = 0; i < num_obs; i++)
    {
        if (assignments[obs_locs[i]] == asynch->my_rank)
        {
            current = &sys[obs_locs[i]];
            updata = (UpstreamData*)current->user;

            updata->num_fit_states = 1 + updata->num_upstreams;
            //for(j=0;j<current->num_parents;j++)
            //	updata->num_fit_states += updata->num_upstreams[j];
            updata->num_fit_states *= num_change_states;

            updata->fit_states = (unsigned int*)malloc(updata->num_fit_states * sizeof(unsigned int));
            updata->fit_to_universal = (unsigned int*)malloc(updata->num_fit_states * sizeof(unsigned int));

            //For q
            updata->fit_states[0] = problem_dim;	//q, q here
            updata->fit_to_universal[0] = current->location * assim_dim;
            counter = 1;
            //for(j=0;j<current->num_parents;j++)
            //{
            unsigned int num_upstreams = ((UpstreamData *)updata)->num_upstreams;
            Link **upstreams = ((UpstreamData *)updata)->upstreams;

            for (unsigned int k = 0; k < num_upstreams; k++)
            {
                updata->fit_states[counter] = problem_dim + counter;
                updata->fit_to_universal[counter] = upstreams[k]->location * assim_dim;
                counter++;
            }
            //}


        }
    }

    //free(is_above_gauges);
    //free(above_gauges);
    //return allstates_needed;
}




//Data assimilation model (Model 254 trimmed, q and s_p) ************************************************************************************


//For modifying q and s_p
void InitRoutines_Assim_254_qsp(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external)
{
    UpstreamData* updata = (UpstreamData*)(link->user);
    unsigned int i, problem_dim = 4;	//Number of model eqs

    //For q only
    link->dim = problem_dim + 5         //Model eqs + variational eqs from this link
        + updata->num_upstreams * 2;	    //Variational eqs from upstreams
    //for(i=0;i<link->num_parents;i++)
    //	link->dim += updata->num_upstreams[i] * 2;	
    link->no_ini_start = 4;
    link->diff_start = 0;

    link->num_dense = link->dim - 3;	//Only q, variational eqs
    link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
    link->dense_indices[0] = 0;
    for (i = 1; i < link->num_dense; i++)	link->dense_indices[i] = i + 3;

    link->differential = &TopLayerHillslope_assim_qsp;
    link->algebraic = NULL;
    link->check_state = NULL;
    link->check_consistency = &CheckConsistency_Nonzero_Model252;
    link->solver = &ExplicitRKSolver;
}


int ReadInitData_Assim_254_qsp(
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    double *y_0, unsigned int dim,
    void *user)
{
    //For this type, all initial conditions for variational equation must be set here.
    unsigned int i;
    unsigned int offset = 4;

    y_0[offset++] = 1.0;  //ds_p/ds_p0

    y_0[offset++] = 0.0;  //ds_t/ds_p0

    y_0[offset++] = 0.0;  //ds_s/ds_p0

    y_0[offset++] = 1.0;  //dq/dq_0
    y_0[offset++] = 0.0;  //dq/ds_p0

    for (i = offset; i < dim; i++)
        y_0[i] = 0.0;	//From upstreams

    return 0;
}


//Function for simple river system with data assimilation.
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: A_i,L_i,A_h | invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2       3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10     
//This uses the units and functions from September 18, 2011 document
//y_i[0] = q, y_i[1] = s, followed by N entries for the variational equation
void TopLayerHillslope_assim_qsp(
    double t,
    const double * const y_i, unsigned int dim,
    const double * const y_p, unsigned short num_parents, unsigned int max_dim,
    const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned int i, j;

    UpstreamData *updata = (UpstreamData*)user;

    double lambda_1 = global_params[1];
    double k_3 = global_params[4];	//[1/min]
    double h_b = global_params[6];	//[m]
    double S_L = global_params[7];	//[m]
    double A = global_params[8];
    double B = global_params[9];
    double exponent = global_params[10];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]

    double L = params[1];	//[m]
    double A_h = params[2];	//[m^2]
    double invtau = params[3];	//[1/min]
    double k_2 = params[4];	//[1/min]
    double k_i = params[5];	//[1/min]
    double c_1 = params[6];
    double c_2 = params[7];

    double q = y_i[0];		//[m^3/s]
    double s_p = y_i[1];	//[m]
    double s_t = y_i[2];	//[m]
    double s_s = y_i[3];	//[m]

    //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    //A few calculations...
    double q_to_lambda_1 = pow(q, lambda_1);
    double q_to_lambda_1_m1 = (q > 1e-12) ? q_to_lambda_1 / q : pow(1e-12, lambda_1 - 1.0);
    double remaining = 1.0 - s_t / S_L;
    double pow_term = (remaining > 0.0) ? pow(remaining, exponent) : 0.0;
    double pow_term_m1 = (remaining > 1e-12) ? pow_term / remaining : pow(1e-12, exponent - 1.0);
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_sl = k_3 * s_s;	//[m/min]

    //Discharge
    double inflow = 0.0;
    ans[0] = -q + (q_pl + q_sl) * c_2;
    for (unsigned int i = 0, p = 0; i < num_parents; i++, p += updata->parents[i]->dim)
        inflow += y_p[p];
    ans[0] = invtau * q_to_lambda_1 * (inflow + ans[0]);

    //Hillslope
    ans[1] = forcing_values[0] * c_1 - q_pl - q_pt - e_p;
    ans[2] = q_pt - q_ts - e_t;
    ans[3] = q_ts - q_sl - e_s;


    //Init for variational equations
    unsigned int offset = 4, problem_dim = 4;
    for (i = offset; i < dim; i++)
        ans[i] = 0.0;

    //Compute partial derivatives (local variables)
    double dfq_dq = lambda_1 * invtau * q_to_lambda_1_m1 * (-q + c_2*(k_2*s_p + k_3*s_s) + inflow) - invtau * q_to_lambda_1;
    double dfq_dsp = invtau*q_to_lambda_1*c_2*k_2;
    double dfq_dss = invtau*q_to_lambda_1*c_2*k_3;
    double dfsp_dsp = -k_2 - k_t;
    double dfsp_dst = k_2 / S_L*B*exponent*pow_term_m1*s_p;
    double dfst_dsp = k_t;
    double dfst_dst = -dfsp_dst - k_i;
    double dfss_dst = k_i;
    double dfss_dss = -k_3;

    //Compute partial derivatives (upstreams variables)
    double dfq_dupq = invtau*q_to_lambda_1;

    const unsigned int num_variational_eqs = 5;

    //Hillslope variational eqs
    ans[offset] = dfsp_dsp * y_i[offset] + dfsp_dst * y_i[offset + 1];	//s_p, s_p
    ans[offset + 1] = dfst_dsp * y_i[offset] + dfst_dst * y_i[offset + 1];	//s_t, s_p
    ans[offset + 2] = dfss_dst * y_i[offset + 1] + dfss_dss * y_i[offset + 2];	//s_s, s_p

    //Discharge variational eqs from this link
    ans[offset + 3] = dfq_dq * y_i[offset + 3]; //q, q
    ans[offset + 4] = dfq_dq * y_i[offset + 4] + dfq_dsp * y_i[offset + 0] + dfq_dss * y_i[offset + 2]; //q, s_p

    //Discharge variational eqs from parent links
    unsigned int current_idx = offset + num_variational_eqs, parent_idx;
    for (i = 0; i < num_parents; i++)
    {
        parent_idx = offset + 3;
        unsigned int num_upstreams = ((UpstreamData *)updata->upstreams[i]->user)->num_upstreams;
        for (j = 0; j < num_upstreams; j++)
        {
            // TODO
            //ans[current_idx] = dfq_dupq * y_p[i][parent_idx] + dfq_dq * y_i[current_idx]; //q, upq
            //ans[current_idx + 1] = dfq_dupq * y_p[i][parent_idx + 1] + dfq_dq * y_i[current_idx + 1]; //q, ups_p
            current_idx += 2;
            parent_idx += 2;
        }
    }
}

//I think these values need to be set only for assigned links, not getting.
//fit_states[i] holds the index in each state vector of the ith sensitivity at this link.
//fit_to_universal[i] holds universal index of the ith sensitivity at this link.
//These only store those sensitivites used for the fitting.
void Setup_Fitting_Data_Model254_qsp(AsynchSolver* asynch, unsigned int* obs_locs, unsigned int num_obs)
{
    unsigned int i, j, k, my_N = asynch->my_N;
    int *assignments = asynch->assignments;
    Link *sys = asynch->sys, **my_sys = asynch->my_sys, *current;
    UpstreamData *updata;

    //Number of states to fit
    unsigned int counter;
    unsigned int problem_dim = 4;	//!!!! Should be allowed to vary by link !!!!
    unsigned int assim_dim = 4; 	//!!!! Should be allowed to vary by link !!!!
    unsigned int num_change_states = 2;	//For q and s_p

    //Find links upstreams from gauges
    bool *is_above_gauges;
    unsigned int *above_gauges;
    unsigned int num_above = GaugeDownstream(asynch, obs_locs, num_obs, &above_gauges, &is_above_gauges);

    //Calculate the number of states needed for the fitting
    //unsigned int allstates_needed = num_above * 2;	//For q and s_p

    for (i = 0; i < num_obs; i++)
    {
        if (assignments[obs_locs[i]] == asynch->my_rank)
        {
            current = &sys[obs_locs[i]];
            updata = (UpstreamData*)current->user;

            updata->num_fit_states = 1 + updata->num_upstreams;
            //for(j=0;j<current->num_parents;j++)
            //	updata->num_fit_states += updata->num_upstreams[j];
            updata->num_fit_states *= num_change_states;

            updata->fit_states = (unsigned int*)malloc(updata->num_fit_states * sizeof(unsigned int));
            updata->fit_to_universal = (unsigned int*)malloc(updata->num_fit_states * sizeof(unsigned int));

            //For q, s_p
            updata->fit_states[0] = problem_dim + 3;	//q, q here
            updata->fit_to_universal[0] = current->location * assim_dim;
            updata->fit_states[1] = problem_dim + 4;	//q, s_p here
            updata->fit_to_universal[1] = current->location * assim_dim + 1;
            counter = 2;
            for (j = 0; j < current->num_parents; j++)
            {
                unsigned int num_upstreams = ((UpstreamData *)updata->upstreams[i]->user)->num_upstreams;
                Link **upstreams = ((UpstreamData *)updata->upstreams[i]->user)->upstreams;

                for (k = 0; k < num_upstreams; k++)
                {
                    updata->fit_states[counter] = problem_dim + 3 + counter;
                    updata->fit_to_universal[counter] = upstreams[k]->location * assim_dim;
                    updata->fit_states[counter + 1] = problem_dim + 4 + counter;
                    updata->fit_to_universal[counter + 1] = upstreams[k]->location * assim_dim + 1;
                    counter += 2;
                }
            }
            /*
            printf("ID = %u\n",current->ID);
            for(k=0;k<current->num_parents;k++)
                printf("upstreams = %u\n",updata->num_upstreams[k]);
            for(k=0;k<updata->num_fit_states;k++)
            {
            printf("%u ",updata->fit_to_universal[k]);
            }
            printf("\n++++++\n");
            */
        }

    }


    free(is_above_gauges);
    free(above_gauges);
    //return allstates_needed;
}



//Data assimilation model (Model 254, q and s_t) ************************************************************************************


//For modifying q and s_t
void InitRoutines_Assim_254_qst(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* external)
{
    UpstreamData* updata = (UpstreamData*)(link->user);
    unsigned int i, problem_dim = 4;	//Number of model eqs

    link->dim = problem_dim + 5	    //Model eqs + variational eqs from this link
        + updata->num_upstreams * 2;  //Variational eqs from upstreams
    //for(i=0;i<link->num_parents;i++)
    //	link->dim += updata->num_upstreams[i] * 2;	//Variational eqs from upstreams
    link->no_ini_start = 4;
    link->diff_start = 0;

    link->num_dense = link->dim - 3;	//Only q, variational eqs	!!!! Do all of the variational eqs really need to be passed down? !!!!
    link->dense_indices = (unsigned int*)realloc(link->dense_indices, link->num_dense * sizeof(unsigned int));
    link->dense_indices[0] = 0;
    for (i = 1; i < link->num_dense; i++)	link->dense_indices[i] = i + 3;

    link->differential = &TopLayerHillslope_assim_qst;
    link->algebraic = NULL;
    link->check_state = NULL;
    link->check_consistency = &CheckConsistency_Nonzero_Model252_st;
    link->solver = &ExplicitRKSolver;
}


int ReadInitData_Assim_254_qst(
    const double * const global_params, unsigned int num_global_params,
    const double * const params, unsigned int num_params,
    double *y_0, unsigned int dim,
    void *user)
{
    //For this type, all initial conditions for variational equation must be set here.
    unsigned int i;
    unsigned int offset = 4;

    y_0[offset++] = 0.0;  //ds_p/ds_t0

    y_0[offset++] = 1.0;  //ds_t/ds_t0

    y_0[offset++] = 0.0;  //ds_s/ds_t0

    y_0[offset++] = 1.0;  //dq/dq_0
    y_0[offset++] = 0.0;  //dq/ds_t0

    for (i = offset; i < dim; i++)
        y_0[i] = 0.0;	//From upstreams

    return 0;
}


//Function for river system with data assimilation. Uses model 252/254.
//Calculates the flow using simple parameters, using only the flow q.
//Order of parameters: A_i,L_i,A_h | invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2       3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10     
//y_i[0] = q, y_i[1] = s_p, y_i[2] = s_t, y_i[3] = s_s followed by N entries for the variational equation
//!!!! Note: this actually works out to be the same as the function for qsp, I think... !!!!
void TopLayerHillslope_assim_qst(
    double t,
    const double * const y_i, unsigned int dim,
    const double * const y_p, unsigned short num_parents, unsigned int max_dim,
    const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned int i, j;

    UpstreamData *updata = (UpstreamData*)user;

    double lambda_1 = global_params[1];
    double k_3 = global_params[4];	//[1/min]
    double h_b = global_params[6];	//[m]
    double S_L = global_params[7];	//[m]
    double A = global_params[8];
    double B = global_params[9];
    double exponent = global_params[10];
    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));	//[mm/month] -> [m/min]

    double L = params[1];	//[m]
    double A_h = params[2];	//[m^2]
    double invtau = params[3];	//[1/min]
    double k_2 = params[4];	//[1/min]
    double k_i = params[5];	//[1/min]
    double c_1 = params[6];
    double c_2 = params[7];

    double q = y_i[0];		//[m^3/s]
    double s_p = y_i[1];	//[m]
    double s_t = y_i[2];	//[m]
    double s_s = y_i[3];	//[m]

    //Evaporation
    double e_p, e_t, e_s;
    double Corr = s_p + s_t / S_L + s_s / (h_b - S_L);
    if (e_pot > 0.0 && Corr > 1e-12)
    {
        e_p = s_p * e_pot / Corr;
        e_t = s_t / S_L * e_pot / Corr;
        e_s = s_s / (h_b - S_L) * e_pot / Corr;
    }
    else
    {
        e_p = 0.0;
        e_t = 0.0;
        e_s = 0.0;
    }

    //A few calculations...
    double q_to_lambda_1 = pow(q, lambda_1);
    double q_to_lambda_1_m1 = (q > 1e-12) ? q_to_lambda_1 / q : pow(1e-12, lambda_1 - 1.0);
    double remaining = 1.0 - s_t / S_L;
    double pow_term = (remaining > 0.0) ? pow(remaining, exponent) : 0.0;
    double pow_term_m1 = (remaining > 1e-12) ? pow_term / remaining : pow(1e-12, exponent - 1.0);
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_sl = k_3 * s_s;	//[m/min]

    //Discharge
    double inflow = 0.0;
    ans[0] = -q + (q_pl + q_sl) * c_2;
    for (unsigned int i = 0, p = 0; i < num_parents; i++, p += updata->parents[i]->dim)
        inflow += y_p[p];
    ans[0] = invtau * q_to_lambda_1 * (inflow + ans[0]);

    //Hillslope
    ans[1] = forcing_values[0] * c_1 - q_pl - q_pt - e_p;
    ans[2] = q_pt - q_ts - e_t;
    ans[3] = q_ts - q_sl - e_s;


    //Init for variational equations
    unsigned int offset = 4, problem_dim = 4;
    for (i = offset; i < dim; i++)	ans[i] = 0.0;

    //Compute partial derivatives (local variables)
    double dfq_dq = lambda_1 * invtau * q_to_lambda_1_m1 * (-q + c_2*(k_2*s_p + k_3*s_s) + inflow) - invtau * q_to_lambda_1;
    double dfq_dsp = invtau*q_to_lambda_1*c_2*k_2;
    double dfq_dss = invtau*q_to_lambda_1*c_2*k_3;
    double dfsp_dsp = -k_2 - k_t;
    double dfsp_dst = k_2 / S_L*B*exponent*pow_term_m1*s_p;
    double dfst_dsp = k_t;
    double dfst_dst = -dfsp_dst - k_i;
    double dfss_dst = k_i;
    double dfss_dss = -k_3;

    //Compute partial derivatives (upstreams variables)
    double dfq_dupq = invtau*q_to_lambda_1;

    //Hillslope variational eqs
    ans[offset] = dfsp_dsp * y_i[offset] + dfsp_dst * y_i[offset + 1];	//s_p, s_t
    ans[offset + 1] = dfst_dsp * y_i[offset] + dfst_dst * y_i[offset + 1];	//s_t, s_t
    ans[offset + 2] = dfss_dst * y_i[offset + 1] + dfss_dss * y_i[offset + 2];	//s_s, s_t

    //Discharge variational eqs from this link
    ans[offset + 3] = dfq_dq * y_i[offset + 3]; //q, q
    ans[offset + 4] = dfq_dq * y_i[offset + 4] + dfq_dsp * y_i[offset + 0] + dfq_dss * y_i[offset + 2]; //q, s_t

    //Discharge variational eqs from parent links
    unsigned int current_idx = offset + 5, parent_idx;
    for (i = 0; i < num_parents; i++)
    {
        parent_idx = offset + 3;
        unsigned int num_upstreams = ((UpstreamData *)updata->upstreams[i]->user)->num_upstreams;
        for (j = 0; j < num_upstreams; j++)
        {
            // TODO
            //ans[current_idx] = dfq_dupq * y_p[i][parent_idx] + dfq_dq * y_i[current_idx]; //q, upq
            //ans[current_idx + 1] = dfq_dupq * y_p[i][parent_idx + 1] + dfq_dq * y_i[current_idx + 1]; //q, ups_t
            current_idx += 2;
            parent_idx += 2;
        }
    }
}

//I think these values need to be set only for assigned links, not getting.
//fit_states[i] holds the index in each state vector of the ith sensitivity at this link.
//fit_to_universal[i] holds universal index of the ith sensitivity at this link.
//These only store those sensitivites used for the fitting.
void Setup_Fitting_Data_Model254_qst(AsynchSolver* asynch, unsigned int* obs_locs, unsigned int num_obs)
{
    unsigned int i, j, k, my_N = asynch->my_N;
    int *assignments = asynch->assignments;
    Link *sys = asynch->sys, **my_sys = asynch->my_sys, *current;
    UpstreamData *updata;

    //Number of states to fit
    unsigned int counter;
    unsigned int problem_dim = 4;	//!!!! Should be allowed to vary by link !!!!
    unsigned int assim_dim = 4;	//!!!! Should be allowed to vary by link !!!!
    unsigned int num_change_states = 2;	//For q and s_t

    //Find links upstreams from gauges
    bool *is_above_gauges;
    unsigned int *above_gauges;
    unsigned int num_above = GaugeDownstream(asynch, obs_locs, num_obs, &above_gauges, &is_above_gauges);

    for (i = 0; i < num_obs; i++)
    {
        if (assignments[obs_locs[i]] == asynch->my_rank)
        {
            current = &sys[obs_locs[i]];
            updata = (UpstreamData*)current->user;

            updata->num_fit_states = 1 + updata->num_upstreams;
            //for(j=0;j<current->num_parents;j++)
            //	updata->num_fit_states += updata->num_upstreams[j];
            updata->num_fit_states *= num_change_states;

            updata->fit_states = (unsigned int*)malloc(updata->num_fit_states * sizeof(unsigned int));
            updata->fit_to_universal = (unsigned int*)malloc(updata->num_fit_states * sizeof(unsigned int));

            //For q, s_t
            updata->fit_states[0] = problem_dim + 3;	//q, q here
            updata->fit_to_universal[0] = current->location * assim_dim;
            updata->fit_states[1] = problem_dim + 4;	//q, s_t here
            updata->fit_to_universal[1] = current->location * assim_dim + 2;
            counter = 2;
            for (j = 0; j < current->num_parents; j++)
            {
                unsigned int num_upstreams = ((UpstreamData *)updata->upstreams[i]->user)->num_upstreams;
                Link **upstreams = ((UpstreamData *)updata->upstreams[i]->user)->upstreams;

                for (k = 0; k < num_upstreams; k++)
                {
                    updata->fit_states[counter] = problem_dim + 3 + counter;
                    updata->fit_to_universal[counter] = upstreams[k]->location * assim_dim;
                    updata->fit_states[counter + 1] = problem_dim + 4 + counter;
                    updata->fit_to_universal[counter + 1] = upstreams[k]->location * assim_dim + 2;
                    counter += 2;
                }
            }
            /*
            printf("ID = %u\n",current->ID);
            //for(k=0;k<current->num_parents;k++)
            //	printf("upstreams = %u\n",updata->num_upstreams[k]);
            for(k=0;k<updata->num_fit_states;k++)
            {
            printf("%u %u\n",updata->fit_states[k],updata->fit_to_universal[k]);
            }
            printf("\n++++++\n");
            */
        }

    }


    free(is_above_gauges);
    free(above_gauges);
    //return allstates_needed;
}

void CheckConsistency_Nonzero_Model252_st(double *y, unsigned int dim, const double * const global_params, unsigned int num_global_params, const double * const params, unsigned int num_params, void *user)
{
    unsigned int i, problem_dim = 4;

    if (y[0] < 1e-14)	y[0] = 1e-14;
    if (y[1] > global_params[7])		y[1] = global_params[7];
    for (i = 1; i < problem_dim; i++)
        if (y[i] < 0.0)	y[i] = 0.0;
}
