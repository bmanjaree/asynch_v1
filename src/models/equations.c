/******************************************************************************************************
 * Includes
 ******************************************************************************************************/
#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <minmax.h>
#include <models/equations.h>

/******************************************************************************************************
 * Descriptions of  different variables in the models (remove later?)
 ******************************************************************************************************/
//extern int flaggy;
//These are the right-hand side functions for the ODEs. For each of these functions:
//double t: The current time
//double *y_i: The approximate value of the solution to the ODEs at the current link at time t
//double *y_p: y_p[j] has the approximate values for the immediately upstream link j to the current link at time t
//unsigned short num_parents: The number of upstream links (parents) to link i
//double *global_params: The global parameters
//RainData* rain: The rain fall values for link i
//double *params: The parameters for link i
//int state: The current state of the system
//double *ans (set by method, assumed that space is allocated): The value returned by the right-hand side function
// double sq(double x) { return x * x; } #why is this here?

/******************************************************************************************************
 * Model 100s Routing only models
 ******************************************************************************************************/

//Type 100
//Non-linear routing equation with two parameter velocity equation globally
//Order of parameters: A_i, L_i, A_h, invtau
//The numbering is:	0   1   2   3
//Order of global_params: v_0,lambda_1 
//The numbering is:        0      1
void routing_100(double t, \
    const double * const y_i, \
    unsigned int dim, \
    const double * const y_p, \
    unsigned short num_parents, \
    unsigned int max_dim, \
    const double * const global_params, \
    const double * const params, \
    const double * const forcing_values, \
    const QVSData * const qvs, \
    int state, \
    void* user, \
    double *ans)
{
    unsigned short i;

    // Spatially variant parameters
	double L_i = params[1]; //Length of stream km -> m converted in definitions.c
	double A_h = params[2]; //hillslope area km2 -> m2 converted in definitions.c
    double invtau = params[3]; // 60.0*v_0 / ((1.0 - lambda_1)*L_i) [1/min]  invtau

    // Global Parameters
    double lambda_1 = global_params[1]; //discharge exponent
    
    //Discharge and Runoff forcing
    double q = y_i[0];		                                        // [m^3/s]
    double runoff = forcing_values[0]*(0.001 / 60.0);               //(mm/hr ->m/min)

    //Nonlinear routing equation
    ans[0] = -q + (runoff * A_h / 60.0); //m3/min to m3/s
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];
}

//Type 101
//Non-linear routing equation with two parameter velocity equation spatially variant
//Order of parameters: A_i, L_i, A_h, v_0, lambda_1, invtau
//The numbering is:	0   1   2   3   4   5
void routing_101(double t, \
    const double * const y_i, \
    unsigned int dim, \
    const double * const y_p, \
    unsigned short num_parents, \
    unsigned int max_dim, \
    const double * const global_params, \
    const double * const params, \
    const double * const forcing_values, \
    const QVSData * const qvs, \
    int state, \
    void* user, \
    double *ans)
{
    unsigned short i;

    // Spatially variant parameters
	double L_i = params[1]; //Length of stream km -> m converted in definitions.c
	double A_h = params[2]; //hillslope area km2 -> m2 converted in definitions.c
    double lambda_1 = params[4]; //discharge exponent
    double invtau = params[5]; // 60.0*v_0 / ((1.0 - lambda_1)*L_i) [1/min]  invtau

    //Discharge and Runoff forcing
    double q = y_i[0];		                                        // [m^3/s]
    double runoff = forcing_values[0]*(0.001 / 60.0);               //(mm/hr ->m/min)

    //Nonlinear routing equation
    ans[0] = -q + (runoff * A_h / 60.0); //m3/min to m3/s
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];
}


//Type 102
//Non-linear routing equation with three drainage area dependent parameter velocity equation globally
//Order of parameters: A_i,L_i,A_h,invtua
//The numbering is:	0   1   2   3
//Order of global_params: v_r,lambda_1,lambda_2
//The numbering is:        0      1        2
void routing_102(double t, \
    const double * const y_i, \
    unsigned int dim, \
    const double * const y_p, \
    unsigned short num_parents, \
    unsigned int max_dim, \
    const double * const global_params, \
    const double * const params, \
    const double * const forcing_values, \
    const QVSData * const qvs, \
    int state, \
    void* user, \
    double *ans)
{
    unsigned short i;

    // Spatially variant parameters
	double L_i = params[1]; //Length of stream km -> m converted in definitions.c
	double A_h = params[2]; //hillslope area km2 -> m2 converted in definitions.c
    double invtau = params[3]; // 60.0*v_0*pow(A_i,lambda_2) / ((1.0 - lambda_1)*L_i) [1/min]  invtau

    // Global Parameters
    double lambda_1 = global_params[1]; //discharge exponent
    
    //Discharge and Runoff forcing
    double q = y_i[0];		                                        // [m^3/s]
    double runoff = forcing_values[0]*(0.001 / 60.0);               //(mm/hr ->m/min)

    //Nonlinear routing equation
    ans[0] = -q + (runoff * A_h / 60.0); //m3/min to m3/s
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];
}

//Type 103
//Non-linear routing equation with three drainage area dependent parameter velocity equation spatially variant
//Order of parameters: A_i,L_i,A_h,v_0,lambda_1,lambda_2,invtua
//The numbering is:	0   1   2   3   4   5   6
void routing_103(double t, \
    const double * const y_i, \
    unsigned int dim, \
    const double * const y_p, \
    unsigned short num_parents, \
    unsigned int max_dim, \
    const double * const global_params, \
    const double * const params, \
    const double * const forcing_values, \
    const QVSData * const qvs, \
    int state, \
    void* user, \
    double *ans)
{
    unsigned short i;

    // Spatially variant parameters
	double L_i = params[1]; //Length of stream km -> m converted in definitions.c
	double A_h = params[2]; //hillslope area km2 -> m2 converted in definitions.c
    double lambda_1 = params[4]; //discharge exponent
    double invtau = params[6]; // 60.0*v_0 / ((1.0 - lambda_1)*L_i) [1/min]  invtau
    
    //Discharge and Runoff forcing
    double q = y_i[0];		                                        // [m^3/s]
    double runoff = forcing_values[0]*(0.001 / 60.0);               //(mm/hr ->m/min)

    //Nonlinear routing equation
    ans[0] = -q + (runoff * A_h / 60.0); //m3/min to m3/s
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];
}


/******************************************************************************************************
 * Model 200s Now changed to runoff only models
 ******************************************************************************************************/



/******************************************************************************************************
 * Model 400s Original Combined TETIS Runoff Routing models
 ******************************************************************************************************/
//Type 253 , 255,256 Reservoir modeling structure that we will change later
//Contains 3 layers on hillslope: ponded, top layer, soil
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2     3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10
void TopLayerHillslope_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{	
	if(forcing_values[2] >0){
		ans[0] = forcing_values[2];
	}
    if(forcing_values[2] <=0){
		unsigned short i;
		for (i = 0; i<num_parents; i++)
			ans[0] += y_p[i * dim];
	}
    ans[1] = 0.0;
    ans[2] = 0.0;
    ans[3] = 0.0;
}


//Type 400
//Tetis model structure for runoff generation + normal routing (NO stream order based velocity)
//Four layers
//Global parameters:
//The numbering is:    0   1   2   3      4      5   6   7...8

//y_i: vector with model states to be resolved by the solver
//dim:scalar with number of dimensions (states?) of the model
//y_p
//num_parents: number of tributary links
//max_dim
//global_params: global parameters applied to all hillslopes. See Precalculations in definitions.c
//params: distributed parameters per hillslope. see Precalculations in definitions.c
//forcing_values:
//qvs: unused
//state: unused
//user: unused
//ans: unused
void model400(double t, \
		const double * const y_i, \
		unsigned int dim, \
		const double * const y_p, \
		unsigned short num_parents, \
		unsigned int max_dim, \
		const double * const global_params, \
		const double * const params, \
		const double * const forcing_values, \
		const QVSData * const qvs, \
		int state, \
		void* user, \
		double *ans)
{

	 	unsigned short i; //auxiliary variable for loops
	    double L = params[1];   // Length of the channel [m]
	    double A_h = params[2]; //Area of the hillslopes [m^2]
	    double c_1 = params[4]; //factor .converts [mm/hr] to [m/min]
	    double rainfall = forcing_values[0] * c_1; //rainfall. from [mm/hr] to [m/min]
	    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));//potential et[mm/month] -> [m/min]
        double temperature = forcing_values[2]; //daily temperature in Celsius
        double temp_thres=global_params[10]; // celsius degrees
        double melt_factor = global_params[9] *(1/(24*60.0)) *(1/1000.0); // mm/day/degree to m/min/degree
        double frozen_ground = forcing_values[3]; // 1 if ground is frozen, 0 if not frozen 
        double x1 =0;

        //states
        unsigned int STATE_DISCHARGE=0;
        unsigned int STATE_STATIC= 1;
        unsigned int STATE_SURFACE=2;
        unsigned int STATE_SUBSURF=3;
        unsigned int STATE_GW = 4;
        unsigned int STATE_SNOW = 5;

        //INITIAL VALUES
        double h5 = y_i[STATE_SNOW];//snow storage [m]
		double h1 = y_i[STATE_STATIC]; //static storage [m]
		double h2 = y_i[STATE_SURFACE];//water in the hillslope surface [m]
		double h3 = y_i[STATE_SUBSURF]; //water in the gravitational storage in the upper part of soil [m]
		double h4 = y_i[STATE_GW]; //water in the aquifer storage [m]
	    double q = y_i[STATE_DISCHARGE];      //[m^3/s]

        //snow storage
        //temperature =0 is the flag for no forcing the variable. no snow process
        if(temperature==0){
            x1 = rainfall;
            ans[STATE_SNOW]=0;
        }
        else{
            if(temperature>=temp_thres){
                double snowmelt = min(h5,temperature * melt_factor); // in [m]
                ans[STATE_SNOW]=-snowmelt; //melting outs of snow storage
                x1 = rainfall + snowmelt; // in [m]
               // printf("temp > th: %f\n", temperature);
               // printf("snowmelt : %f\n", snowmelt);
            }
            if(temperature != 0 & temperature <temp_thres){
                ans[STATE_SNOW]=rainfall; //all precipitation is stored in the snow storage
                x1=0;
                //printf("temp < th: %f\n", temperature);
            }
        }
        

		//static storage
		double Hu = global_params[3]/1000; //max available storage in static tank [mm] to [m]
		double x2 = max(0,x1 + h1 - Hu ); //excedance flow to the second storage [m] [m/min] check units
        //if ground is frozen, x1 goes directly to the surface
        //therefore nothing is diverted to static tank
        if(frozen_ground == 1){
            x2 = x1;
        }
            
		double d1 = x1 - x2; // the input to static tank [m/min]
		double out1 = min(e_pot, h1); //evaporation from the static tank. it cannot evaporate more than h1 [m]
		//double out1 = (e_pot > h1) ? e_pot : 0.0;
		ans[STATE_STATIC] = d1 - out1; //differential equation of static storage


		//surface storage tank
		double infiltration = global_params[4]*c_1; //infiltration rate [m/min]
         if(frozen_ground == 1){
            infiltration = 0;
        }
		double x3 = min(x2, infiltration); //water that infiltrates to gravitational storage [m/min]
		double d2 = x2 - x3; // the input to surface storage [m] check units
        double alfa2 =global_params[6]; //velocity in m/s
        double w = alfa2 * L / A_h  * 60; // [1/min]
        w = min(1,w); //water can take less than 1 min (dt) to leave surface
        double out2 =0;
        out2  = h2 * w; //direct runoff [m/min]
		ans[STATE_SURFACE] = d2 - out2; //differential equation of surface storage


		// SUBSURFACE storage
		double percolation = global_params[5]*c_1; // percolation rate to aquifer [m/min]
		double x4 = min(x3,percolation); //water that percolates to aquifer storage [m/min]
		double d3 = x3 - x4; // input to gravitational storage [m/min]
		double alfa3 = global_params[7]* 24*60; //residence time [days] to [min].
        double out3=0;
        if(alfa3>=1)
		    out3 = h3/alfa3; //interflow [m/min]
		ans[STATE_SUBSURF] = d3 - out3; //differential equation for gravitational storage

		//aquifer storage
		double x5 = 0;//water loss to deeper aquifer [m]
		double d4 = x4 - x5;
		double alfa4 = global_params[8]* 24*60; //residence time [days] to [min].
        double out4=0;
        if(alfa4>=1)
		    out4 = h4/alfa4 ; //base flow [m/min]
		ans[STATE_GW] = d4 - out4; //differential equation for aquifer storage

		//channel storage
		double lambda_1 = global_params[1];
	    double invtau = params[3];// 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
	   	double c_2 = params[5];// = A_h / 60.0;	//  c_2

	    ans[STATE_DISCHARGE] = -q + (out2 + out3 + out4) * c_2; //[m/min] to [m3/s]
	    for (i = 0; i < num_parents; i++)
	        ans[STATE_DISCHARGE] += y_p[i * dim + STATE_DISCHARGE];
	    ans[STATE_DISCHARGE] = invtau * pow(q, lambda_1) * ans[STATE_DISCHARGE];    // discharge[0]
}

//Type 404
//Tetis model structure for runoff generation + normal routing (NO stream order based velocity)
//Four layers
//Global parameters:
//The numbering is:    0   1   2   3      4      5   6   7...8

//y_i: vector with model states to be resolved by the solver
//dim:scalar with number of dimensions (states?) of the model
//y_p
//num_parents: number of tributary links
//max_dim
//global_params: global parameters applied to all hillslopes. See Precalculations in definitions.c
//params: distributed parameters per hillslope. see Precalculations in definitions.c
//forcing_values:
//qvs: unused
//state: unused
//user: unused
//ans: unused
void model404(double t, \
		const double * const y_i, \
		unsigned int dim, \
		const double * const y_p, \
		unsigned short num_parents, \
		unsigned int max_dim, \
		const double * const global_params, \
		const double * const params, \
		const double * const forcing_values, \
		const QVSData * const qvs, \
		int state, \
		void* user, \
		double *ans)
{

	 	unsigned short i; //auxiliary variable for loops

        double c_1 = (0.001 / 60.0);		//(mm/hr->m/min)  c_1
        double A_i = params[0]; //drainage area in km2
	    double L = params[1];   // Length of the channel [m]
	    double A_h = params[2]; //Area of the hillslopes [m^2] 
        double c_2 = A_h / 60.0;	//  c_2
        double v_0 = params[3];
        double lambda_1 = params[4];
        double lambda_2 = params[5];
        double Hu = params[6]/1000; //[m]
        double infiltration = params[7]*c_1; //infiltration rate [m/min]
		double percolation = params[8]*c_1; // percolation rate to aquifer [m/min]
        double alfa2 =params[9]; //velocity in m/s 
		double alfa3 = params[10]* 24*60; //residence time [days] to [min].
		double alfa4 = params[11]* 24*60; //residence time [days] to [min].
		double melt_factor = params[12] *(1/(24*60.0)) *(1/1000.0); // mm/day/degree to m/min/degree
        double temp_thres= params[13]; // celsius degrees
        double invtau = 60.0 * v_0 * pow(A_i, lambda_2) / ((1.0 - lambda_1) * L);//[1/min]  invtau params[3] !!!here


        double x1 =0;


        //forcings
	    double rainfall = forcing_values[0] * c_1; //rainfall. from [mm/hr] to [m/min]
	    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));//potential et[mm/month] -> [m/min]
        double temperature = forcing_values[2]; //daily temperature in Celsius
        double frozen_ground = forcing_values[3]; // 1 if ground is frozen, 0 if not frozen 
        
        //states
        unsigned int STATE_DISCHARGE=0;
        unsigned int STATE_STATIC= 1;
        unsigned int STATE_SURFACE=2;
        unsigned int STATE_SUBSURF=3;
        unsigned int STATE_GW = 4;
        unsigned int STATE_SNOW = 5;

        //INITIAL VALUES
        double h5 = y_i[STATE_SNOW];//snow storage [m]
		double h1 = y_i[STATE_STATIC]; //static storage [m]
		double h2 = y_i[STATE_SURFACE];//water in the hillslope surface [m]
		double h3 = y_i[STATE_SUBSURF]; //water in the gravitational storage in the upper part of soil [m]
		double h4 = y_i[STATE_GW]; //water in the aquifer storage [m]
	    double q = y_i[STATE_DISCHARGE];      //[m^3/s]

        //snow storage
        //temperature =0 is the flag for no forcing the variable. no snow process
        if(temperature==0){
            x1 = rainfall;
            ans[STATE_SNOW]=0;
        }
        else{
            if(temperature>=temp_thres){
                double snowmelt = min(h5,temperature * melt_factor); // in [m]
                ans[STATE_SNOW]=-snowmelt; //melting outs of snow storage
                x1 = rainfall + snowmelt; // in [m]
               // printf("temp > th: %f\n", temperature);
               // printf("snowmelt : %f\n", snowmelt);
            }
            if(temperature != 0 & temperature <temp_thres){
                ans[STATE_SNOW]=rainfall; //all precipitation is stored in the snow storage
                x1=0;
                //printf("temp < th: %f\n", temperature);
            }
        }
        

		//static storage
		 //max available storage in static tank [mm] to [m]
		double x2 = max(0,x1 + h1 - Hu ); //excedance flow to the second storage [m] [m/min] check units
        //if ground is frozen, x1 goes directly to the surface
        //therefore nothing is diverted to static tank
        if(frozen_ground == 1){
            x2 = x1;
        }
            
		double d1 = x1 - x2; // the input to static tank [m/min]
		double out1 = min(e_pot, h1); //evaporation from the static tank. it cannot evaporate more than h1 [m]
		//double out1 = (e_pot > h1) ? e_pot : 0.0;
		ans[STATE_STATIC] = d1 - out1; //differential equation of static storage


		//surface storage tank
         if(frozen_ground == 1){
            infiltration = 0;
        }
		double x3 = min(x2, infiltration); //water that infiltrates to gravitational storage [m/min]
		double d2 = x2 - x3; // the input to surface storage [m] check units
        //double alfa2 =global_params[6]; //velocity in m/s
        
        
        //double alfa2 = (pow(h2,2./3.) * pow(slope,0.5))/manning; //m/s

        double w = alfa2 * L / A_h  * 60; // [1/min] !!!here
        w = min(1,w); //water can take less than 1 min (dt) to leave surface
        double out2 =0;
        out2  = h2 * w; //direct runoff [m/min]
		ans[STATE_SURFACE] = d2 - out2; //differential equation of surface storage


		// SUBSURFACE storage
		double x4 = min(x3,percolation); //water that percolates to aquifer storage [m/min]
		double d3 = x3 - x4; // input to gravitational storage [m/min]
        double out3=0;
        if(alfa3>=1)
		    out3 = h3/alfa3; //interflow [m/min]
		ans[STATE_SUBSURF] = d3 - out3; //differential equation for gravitational storage

		//aquifer storage
		double x5 = 0;//water loss to deeper aquifer [m]
		double d4 = x4 - x5;
        double out4=0;
        if(alfa4>=1)
		    out4 = h4/alfa4 ; //base flow [m/min]
		ans[STATE_GW] = d4 - out4; //differential equation for aquifer storage

		//channel storage

		//double lambda_1 = params[1];
	    //double invtau = params[3];// 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
	   	//double c_2 = params[5];// = A_h / 60.0;	//  c_2

	    ans[STATE_DISCHARGE] = -q + (out2 + out3 + out4) * c_2; //[m/min] to [m3/s]
	    for (i = 0; i < num_parents; i++)
	        ans[STATE_DISCHARGE] += y_p[i * dim + STATE_DISCHARGE];
	    ans[STATE_DISCHARGE] = invtau * pow(q, lambda_1) * ans[STATE_DISCHARGE];    // discharge[0]
}
