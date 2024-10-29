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



double sq(double x) { return x * x; }


//Type 253 , 255,256
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


//Type 254
//Contains 3 layers on hillslope: ponded, top layer, soil. Also has 3 extra states: total precip, total runoff, base flow
//Order of parameters: A_i,L_i,A_h,invtau,k_2,k_i,c_1,c_2
//The numbering is:	0   1   2     3    4   5   6   7
//Order of global_params: v_0,lambda_1,lambda_2,v_h,k_3,k_I_factor,h_b,S_L,A,B,exponent,v_B
//The numbering is:        0      1        2     3   4     5        6   7  8 9  10       11
void model254(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

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
                                //double h_r = params[3];	//[m]
    double invtau = params[3];	//[1/min]
    double k_2 = params[4];	//[1/min]
    double k_i = params[5];	//[1/min]
    double c_1 = params[6];
    double c_2 = params[7];

    double q =   y_i[0];		//[m^3/s]
    double s_p = y_i[1];	//[m]
    double s_t = y_i[2];	//[m]
    double s_s = y_i[3];	//[m]
                            //double s_precip = y_i[4];	//[m]
                            //double V_r = y_i[5];	//[m^3]
    double q_b = max(0.001,y_i[6]);	//[m^3/s]

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

    double pow_term = (1.0 - s_t / S_L > 0.0) ? pow(1.0 - s_t / S_L, exponent) : 0.0;
    double k_t = (A + B * pow_term) * k_2;

    //Fluxes
    double q_pl = k_2 * s_p;
    double q_pt = k_t * s_p;
    double q_ts = k_i * s_t;
    double q_sl = k_3 * s_s;	//[m/min]

                                //Discharge
    ans[0] = -q + (q_pl + q_sl) * c_2;
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];

    //Hillslope
    ans[1] = forcing_values[0] * c_1 - q_pl - q_pt - e_p;
    ans[2] = q_pt - q_ts - e_t;
    ans[3] = q_ts - q_sl - e_s;

    //Additional states

    ans[4] = forcing_values[0] * c_1;
    ans[5] = q_pl;


    ans[6] = q_sl * A_h - q_b*60.0;
    for (i = 0; i<num_parents; i++)
        ans[6] += y_p[i * dim + 6] * 60.0;
    //ans[6] += k_3*y_p[i].ve[3]*A_h;
    ans[6] *= v_B / L;
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

        // if (forcing_values[0]>1 && ratio<1) {
        //     printf("time: %f\n", t);
        //     printf(" rain in mm/hour: %f\n", forcing_values[0]);
        //     printf(" area hill, area basin, area ratio: %f %f %f\n", A_h,A_i,ratio);
        //     MPI_Abort(MPI_COMM_WORLD, 1);
        // }

}
//Type 401
//model 400 + tracking of total runoff, surface runoff, subsurface runoff and groundwater runoff
void model401(double t, \
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
        double A_i = params[0]; //Area of the basin [km^2]
         A_i*=1e6; //[m^2]
	    double c_1 = params[4]; //factor .converts [mm/hr] to [m/min]
	    double rainfall = forcing_values[0] * c_1; //rainfall. from [mm/hr] to [m/min]
	    double e_pot = forcing_values[1] * (1e-3 / (30.0*24.0*60.0));//potential et[mm/month] -> [m/min]
		double temperature = forcing_values[2]; //daily temperature in Celsius
        double temp_thres=global_params[10]; // celsius degrees
        double melt_factor = global_params[9] *(1/(24*60.0)) *(1/1000.0); // mm/day/degree to m/min/degree
        double frozen_ground = forcing_values[3]; // 1 if ground is frozen, 0 if not frozen 
        double x1 =0;

        // 9 states
        // i need to put the fluxes on top because cant print more than state7. bug
        //y0=q discharge[m3/s]
        //y1 = basin rainfall and snowmelt [m3/hour]
        //y2 = basin surface runoff [m3/hour]
        //y3= basin subsurface runoff [m3/hour]
        //y4 = basin gw rounoff [m3/hour]
        //y5= h1 static storage[m]
        //y6= h2 water hill surface[m]
        //y7 = h3 water upper soil [m]
        //y8 = h4  water lower soil [m]
        //y9 = h5 snow storage [m]
        unsigned int STATE_DISCHARGE=0;
        unsigned int STATE_CUMRAINFALL=1;
        unsigned int STATE_CUMSURF = 2;
        unsigned int STATE_CUMSUB=3;
        unsigned int STATE_CUMGW = 4;
        unsigned int STATE_STATIC=5;
        unsigned int STATE_SURFACE=6;
        unsigned int STATE_SUBSURF =7;
        unsigned int STATE_GW = 8;
        unsigned int STATE_SNOW = 9;

        #ifndef minf
        #define minf(a,b) ((a) < (b) ? (a) : (b))
        #endif

        #ifndef maxf
        #define maxf(a,b) ((a) > (b) ? (a) : (b))
        #endif

        //INITIAL VALUES
        double h5 = y_i[STATE_SNOW];//snow storage [m]
        double basin_rainfall = y_i[STATE_CUMRAINFALL]; //[m3/hour] see last lines
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
                double snowmelt = min(h5,temperature * melt_factor); // in [m/min]
                ans[STATE_SNOW]=-snowmelt; //melting outs of snow storage
                x1 = rainfall + snowmelt; // in [m/min]
                //printf("temp > th: %f\n", temperature);
                //printf("snowmelt : %f\n", snowmelt);
            }
            if(temperature != 0 & temperature <temp_thres){
                ans[STATE_SNOW]=rainfall; //all precipitation is stored in the snow storage
                x1=0;
                //printf("temp < th: %f\n", temperature);
            }
        }

		//static storage
		double Hu = global_params[3]/1000; //max available storage in static tank [mm] to [m]
        double aux1 = x1 + h1 - Hu;
		double x2 = max(0,aux1); //excedance flow to the second storage [m] [m/min] check units
        //double x2 = (aux1>0)? aux1: 0.0;
        //if ground is frozen, x1 goes directly to the surface
        //therefore nothing is diverted to static tank
        if(frozen_ground == 1){
            x2 = x1;
        }
		double d1 = x1 - x2; // the input to static tank [m/min]
		double out1 = min(e_pot, h1); //evaporation from the static tank. it cannot evaporate more than h1 [m]
		//double out1 = (e_pot > h1) ? e_pot : 0.0;
		ans[STATE_STATIC] = d1 - out1; //differential equation of static storage
        //printf("t %f\n",t);
        //printf("d1 %f\n",d1);
        //printf("x1 and x2 %f %f\n",x1,x2);
        //printf("out1 %f\n", out1);
        //printf(" rain in mm/hour: %f\n", forcing_values[0]);
        //printf("h1 %f\n ",h1);
        //MPI_Abort(MPI_COMM_WORLD, 1);

		//surface storage tank
		double infiltration = global_params[4]*c_1; //infiltration rate [m/min]
        if(frozen_ground == 1){
            infiltration = 0;
        }
		double x3 = min(x2, infiltration); //water that infiltrates to gravitational storage [m/min]
		double d2 = x2 - x3; // the input to surface storage [m] check units
		//double alfa2 = global_params[6]* 24*60; //residence time [days] to [min].
		double alfa2 =global_params[6]; //velocity in m/s
        double w = alfa2 * L / A_h  * 60; // [1/min]
        w = min(w,1); // water can take less than 1 min to
        double out2 =0;
        out2  = h2 * w; //direct runoff [m/min]
		ans[STATE_SURFACE] = d2 - out2; //differential equation of surface storage
        double surface_runoff = y_i[STATE_CUMSURF]; //[m3/hour]

		// SUBSURFACE storage
		double percolation = global_params[5]*c_1; // percolation rate to aquifer [m/min]
		double x4 = min(x3,percolation); //water that percolates to aquifer storage [m/min]
		double d3 = x3 - x4; // input to gravitational storage [m/min]
		double alfa3 = global_params[7]* 24*60; //residence time [days] to [min].
		double out3=0;
        if(alfa3>=1)
		    out3 = h3/alfa3; //interflow [m/min]
		ans[STATE_SUBSURF] = d3 - out3; //differential equation for gravitational storage
        double subsurface_runoff = y_i[STATE_CUMSUB];  //[m3/hour]

		//aquifer storage
		double x5 = 0;//water loss to deeper aquifer [m]
		double d4 = x4 - x5;
		double alfa4 = global_params[8]* 24*60; //residence time [days] to [min].
		double out4=0;
        if(alfa4>=1)
		    out4 = h4/alfa4 ; //base flow [m/min]
		ans[STATE_GW] = d4 - out4; //differential equation for aquifer storage
        double groundwater_runoff = y_i[STATE_CUMGW]; //[m3/hour]

		//channel storage
		double lambda_1 = global_params[1];
	    double invtau = params[3];// 60.0*v_0*pow(A_i, lambda_2) / ((1.0 - lambda_1)*L_i);	//[1/min]  invtau
	   	double c_2 = params[5];// = A_h / 60.0;	//  c_2

	    ans[STATE_DISCHARGE] = -q + (out2 + out3 + out4) * c_2; //[m/min] to [m3/s]

        double aux = forcing_values[0] *(1/1000.0) * A_h;//[mm/h] to [m3/h] 
        ans[STATE_CUMRAINFALL] = -basin_rainfall + aux; //[m3/hour] 
        ans[STATE_CUMSURF] = -surface_runoff + out2 * 60.0 * A_h; //[m/min] to [m3/hour]
        ans[STATE_CUMSUB] = -subsurface_runoff + out3 *60.0 * A_h ; //[m/min] to [m3/hour]
        ans[STATE_CUMGW] = -groundwater_runoff + out4 *60.0 *A_h ; //[m/min] to [m3/hour]

	    for (i = 0; i < num_parents; i++){
            ans[STATE_DISCHARGE] += y_p[i * dim + STATE_DISCHARGE];
            ans[STATE_CUMRAINFALL] += y_p[i * dim +STATE_CUMRAINFALL];
            ans[STATE_CUMSURF] += y_p[i * dim +STATE_CUMSURF];
            ans[STATE_CUMSUB] += y_p[i * dim +STATE_CUMSUB];
            ans[STATE_CUMGW] += y_p[i * dim +STATE_CUMGW];
        }
	        
	    ans[STATE_DISCHARGE] = invtau * pow(q, lambda_1) * ans[STATE_DISCHARGE];    // discharge[0]
}
void model401reservoir(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{	
    unsigned int FORCING_RESERVOIR = 4;
    unsigned int STATE_DISCHARGE=0;
	if(forcing_values[FORCING_RESERVOIR] >=0){
		ans[STATE_DISCHARGE] = forcing_values[FORCING_RESERVOIR];
	}
}

//Type 402
void model402(double t, \
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
        unsigned int STATE_DAM_STORAGE=6;

        //INITIAL VALUES
        double h5 = y_i[STATE_SNOW];//snow storage [m]
		double h1 = y_i[STATE_STATIC]; //static storage [m]
		double h2 = y_i[STATE_SURFACE];//water in the hillslope surface [m]
		double h3 = y_i[STATE_SUBSURF]; //water in the gravitational storage in the upper part of soil [m]
		double h4 = y_i[STATE_GW]; //water in the aquifer storage [m]
	    double q = y_i[STATE_DISCHARGE];      //[m^3/s]
        double dam_storage =y_i[STATE_DAM_STORAGE] ;//m^3

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
        //control surface storage to not become too small. if not, calculation takes longer
        if(h2<1E-3)
            h2 = 1E-2; // this wont make the model run faster or slower
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

        if(state==0){ //no dams
            ans[STATE_DISCHARGE] = -q + (out2 + out3 + out4) * c_2; //[m/min] to [m3/s]
	        for (i = 0; i < num_parents; i++)
	            ans[STATE_DISCHARGE] += y_p[i * dim + STATE_DISCHARGE];
	        ans[STATE_DISCHARGE] = invtau * pow(q, lambda_1) * ans[STATE_DISCHARGE];    // discharge[0]
            //this lines were intented to keep storage at zero on no dam links, but is making the
           dam_storage=1e-1; // this is the model making run slower
           ans[STATE_DAM_STORAGE] = -dam_storage; // this is the model making run slower

        }
        if(state ==1){// dams
            double dam_input = 0;//m3
            double dam_output=0;//m3
            double dam_outflow=0; //m3s-1
            int debug = 0;
            if(debug) printf("time: %f\n", t);
            //inflow from upstream links
            for (i = 0; i < num_parents; i++)
	            dam_input += y_p[i * dim + STATE_DISCHARGE] * 60.0; //m3/s to m3
            if(debug) printf("inflow from upstream in m3: %f\n", dam_input);
            // storage inputs
            if(debug) printf("initial dam storage in m3: %f\n", dam_storage);
            //dam_storage= dam_storage + dam_input;
            if(debug) printf("step1 dam storage in m3: %f\n", dam_storage + dam_input);

            //calculate dam outflow
            dam_outflow = get_discharge_from_storage(dam_storage + dam_input); //m3s-1
            if(debug) printf("dam outflow in m3/s: %f\n", dam_outflow);

            ans[STATE_DISCHARGE] = dam_outflow - q;
            
            //storage output
            //dam_output =get_storage_from_discharge(dam_outflow);
            dam_output = dam_outflow * 60.0; //m3s-1 to m3
            if(debug) printf("dam output storage in m3: %f\n", dam_output);

            //dam_storage =dam_storage - dam_output;
            if(debug) printf("step2 dam storage in m3: %f\n", dam_storage + dam_input - dam_output);

            //dam_storage = max(dam_storage,0);
            //if(debug) printf("step3 dam storage in m3: %f\n", dam_storage);
            //ans[STATE_DAM_STORAGE] = dam_storage;
            ans[STATE_DAM_STORAGE] = dam_input - dam_output;

            if(debug) printf("\n");

            
        }
}

double get_discharge_from_storage(double storage){
    double discharge=0;
    if(storage<125000) //discharge = 1m3s-1
        discharge = (1-0)/(125000-0.0)*storage;
     if(storage>=125000)
        discharge = ((100.0-1.0)/(225000.0-125000.0))*(storage - 125000.0)+1.0;
    return discharge;
}

double get_storage_from_discharge(double discharge){
    double storage=0;
    if(discharge<1)
        storage = (125000-0)/(1-0)*discharge ;
    else
        storage = (225000-125000)/(100.0-1.0)*(discharge-1.0) + 125000; 
    return storage;
}

//Type 403
//model tetis with dams using specific individual ratings 
void model403(double t, \
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
        unsigned int STATE_DAM_STORAGE=6;

        //INITIAL VALUES
        double h5 = y_i[STATE_SNOW];//snow storage [m]
		double h1 = y_i[STATE_STATIC]; //static storage [m]
		double h2 = y_i[STATE_SURFACE];//water in the hillslope surface [m]
		double h3 = y_i[STATE_SUBSURF]; //water in the gravitational storage in the upper part of soil [m]
		double h4 = y_i[STATE_GW]; //water in the aquifer storage [m]
	    double q = y_i[STATE_DISCHARGE];      //[m^3/s]
        double dam_storage =y_i[STATE_DAM_STORAGE] ;//m^3

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

        //printf("state: %d\n", state);
        if(state == 0){ //no dams
            ans[STATE_DISCHARGE] = -q + (out2 + out3 + out4) * c_2; //[m/min] to [m3/s]
	        for (i = 0; i < num_parents; i++)
	            ans[STATE_DISCHARGE] += y_p[i * dim + STATE_DISCHARGE];
	        ans[STATE_DISCHARGE] = invtau * pow(q, lambda_1) * ans[STATE_DISCHARGE];    // discharge[0]
            dam_storage=1e-1;
            ans[STATE_DAM_STORAGE] = -dam_storage;
        }
        //state is the array index corresponding to the current storage in the qvs 
        if(state == 1 ){// dams
            int debug = 0;
            
            if(debug) printf("state: %d\n", state);

            double dam_input = 0;//m3
            double dam_output=0;//m3
            double dam_outflow=0; //m3s-1
            
            
            if(debug) printf("time: %f\n", t);
            //inflow from upstream links
            for (i = 0; i < num_parents; i++)
	            dam_input += y_p[i * dim + STATE_DISCHARGE] * 60.0; //m3/s to m3
            if(debug) printf("inflow from upstream in m3: %f\n", dam_input);
            // storage inputs
            if(debug) printf("initial dam storage in m3: %f\n", dam_storage);
            //dam_storage= dam_storage + dam_input;
            if(debug) printf("step1 dam storage in m3: %f\n", dam_storage + dam_input);

            //if(debug) printf("state: %u\n", state);
            unsigned int max_storage_pond_index = (int)qvs->n_values - 1;
            double max_storage_pond = qvs->points[max_storage_pond_index][0];
            if (dam_storage + dam_input >= max_storage_pond)
            {
                if(debug) printf("storage >= max.sto.pond\n");
                //S_max = qvs->points[qvs->n_values - 1][0];
                double q_max = qvs->points[max_storage_pond_index][1];
                if(debug) printf("qmax : %f\n", q_max);
                dam_outflow = q_max;//m3s-1
            }
            else if (dam_storage + dam_input < max_storage_pond)
            {
                if(debug) printf("storage within rating range\n");

                //S = (y_i[1] < 0.0) ? 0.0 : y_i[1];
                double S = (dam_storage + dam_input < 0.0) ? 0.0 : dam_storage + dam_input;
                if(debug) printf("storage  : %f\n", S);

                //find discharge value in rating
                unsigned int ii;
                for (ii = 0; ii < qvs->n_values - 1; ii++){
				   // if(debug)printf("model403 storage qvs points: %f\n", qvs->points[ii][0]);
                    if (qvs->points[ii][0] <= S
						&& S < qvs->points[ii + 1][0])    
					break;
                }
                double q2 = qvs->points[ii + 1][1];
                double q1 = qvs->points[ii][1];
                double S2 = qvs->points[ii + 1][0];
                double S1 = qvs->points[ii][0];
                dam_outflow = (q2 - q1) / (S2 - S1) * (S - S1) + q1;//m3s-1
                
            }
            
            //calculate dam outflow
            //dam_outflow = get_discharge_from_storage(dam_storage + dam_input); //m3s-1
            if(debug) printf("dam outflow in m3/s: %f\n", dam_outflow);

            ans[STATE_DISCHARGE] = dam_outflow - q;
            
            //storage output
            //dam_output =get_storage_from_discharge(dam_outflow);
            dam_output = dam_outflow * 60.0; //m3s-1 to m3
            if(debug) printf("dam output storage in m3: %f\n", dam_output);

        
            if(debug) printf("step2 dam storage in m3: %f\n", dam_storage + dam_input - dam_output);

            //dam_storage = max(dam_storage,0);
            //if(debug) printf("step3 dam storage in m3: %f\n", dam_storage);
            //ans[STATE_DAM_STORAGE] = dam_storage;
            ans[STATE_DAM_STORAGE] = dam_input - dam_output;

            if(debug) printf("\n");

            
        }

        // if (forcing_values[0]>1 && ratio<1) {
        //     printf("time: %f\n", t);
        //     printf(" rain in mm/hour: %f\n", forcing_values[0]);
        //     printf(" area hill, area basin, area ratio: %f %f %f\n", A_h,A_i,ratio);
        //     MPI_Abort(MPI_COMM_WORLD, 1);
        // }
}

//Type 402
void Tetis03_Reservoirs(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    ans[0] = forcing_values[3];
    ans[1] = 0.0;
    ans[2] = 0.0;
    ans[3] = 0.0;
    ans[4] = 0.0;
    ans[5] = 0.0;
    ans[6] = 0.0;
    ans[7] = 0.0;
    ans[8] = 0.0;
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
        double invtau = 60.0 * v_0 * pow(A_i, lambda_2) / ((1.0 - lambda_1) * L);//[1/min]  invtau params[3]


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

        double w = alfa2 * L / A_h  * 60; // [1/min]
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

        // if (forcing_values[0]>1 && ratio<1) {
        //     printf("time: %f\n", t);
        //     printf(" rain in mm/hour: %f\n", forcing_values[0]);
        //     printf(" area hill, area basin, area ratio: %f %f %f\n", A_h,A_i,ratio);
        //     MPI_Abort(MPI_COMM_WORLD, 1);
        // }

}

//Type 405
//model tetis with dams using specific individual ratings, but using the flowchart algorithm . 
void model405(double t, \
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
        double factor_low_threshold = global_params[11];
        double factor_high_threshold = global_params[12];
        double factor_operating_window = global_params[13];

        //states
        unsigned int STATE_DISCHARGE=0;
        unsigned int STATE_STATIC= 1;
        unsigned int STATE_SURFACE=2;
        unsigned int STATE_SUBSURF=3;
        unsigned int STATE_GW = 4;
        unsigned int STATE_SNOW = 5;
        unsigned int STATE_DAM_STORAGE=6;

        //INITIAL VALUES
        double h5 = y_i[STATE_SNOW];//snow storage [m]
		double h1 = y_i[STATE_STATIC]; //static storage [m]
		double h2 = y_i[STATE_SURFACE];//water in the hillslope surface [m]
		double h3 = y_i[STATE_SUBSURF]; //water in the gravitational storage in the upper part of soil [m]
		double h4 = y_i[STATE_GW]; //water in the aquifer storage [m]
	    double q = y_i[STATE_DISCHARGE];      //[m^3/s]
        double dam_storage =y_i[STATE_DAM_STORAGE] ;//m^3

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

        //printf("state: %d\n", state);
        if(state == 0){ //no dams
            ans[STATE_DISCHARGE] = -q + (out2 + out3 + out4) * c_2; //[m/min] to [m3/s]
	        for (i = 0; i < num_parents; i++)
	            ans[STATE_DISCHARGE] += y_p[i * dim + STATE_DISCHARGE];
	        ans[STATE_DISCHARGE] = invtau * pow(q, lambda_1) * ans[STATE_DISCHARGE];    // discharge[0]
            dam_storage=1e-1;
            ans[STATE_DAM_STORAGE] = -dam_storage;
        }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    // DAM MODEL
        //state is the array index corresponding to the current storage in the qvs 
        if(state == 1 ){// dams
            int debug = 0;
            
            if(debug) printf("state: %d\n", state);

            double dam_input = 0;//m3
            double dam_output=0;//m3
            double dam_outflow=0; //m3s-1
            double rating_flow=0; //m3s-1
            double output_mode=2; //defines reservoir operation scheme. 0 = passive. 1 = conservation pools. 2=incremental outflow increases with storage            
            

            if(debug) printf("time: %f\n", t);
            //inflow from upstream links
            for (i = 0; i < num_parents; i++)
	            dam_input += y_p[i * dim + STATE_DISCHARGE] * 60.0; //m3/s to m3
            dam_input += ((out2 + out3 + out4) * c_2) * 60; //adds hillslope runoff to the channel input
            if(debug) printf("inflow from upstream in m3: %f\n", dam_input);
            // storage inputs
            if(debug) printf("initial dam storage in m3: %f\n", dam_storage);
            //dam_storage= dam_storage + dam_input;
            if(debug) printf("step1 dam storage in m3: %f\n", dam_storage + dam_input);

            //if(debug) printf("state: %u\n", state);
            unsigned int max_storage_pond_index = (int)qvs->n_values - 1;
            double max_storage_pond = qvs->points[max_storage_pond_index][0];
            // find the upper and lower storage threshold for this link
            double low_storage_threshold = max_storage_pond * factor_low_threshold;
            double high_storage_threshold = max_storage_pond * factor_high_threshold;
            //if(dam_storage + dam_input >0) debug=1; //only print results if there is water in the pond
            if(debug) printf("max.sto.pond %f\n",max_storage_pond);
            if(debug) printf("low_storage_factor and high_storage_factor: %f and %f\n", factor_low_threshold,factor_high_threshold);

            if(debug) printf("low_storage_threshold and high_storage_threshold: %f and %f\n", low_storage_threshold,high_storage_threshold);
            if (dam_storage + dam_input >= max_storage_pond)
            {
                if(debug) printf("storage >= max.sto.pond\n");
                //S_max = qvs->points[qvs->n_values - 1][0];
                double q_max = qvs->points[max_storage_pond_index][1];
                if(debug) printf("qmax : %f\n", q_max);
                dam_outflow = q_max;//m3s-1
            }
            else if (dam_storage + dam_input < max_storage_pond)
            {
                if(debug) printf("storage within rating range\n");
                double S = (dam_storage + dam_input < 0.0) ? 0.0 : dam_storage + dam_input;
                if(debug) printf("storage  : %f\n", S);
                //find discharge value in rating
                unsigned int ii;
                for (ii = 0; ii < qvs->n_values - 1; ii++){
				   // if(debug)printf("model403 storage qvs points: %f\n", qvs->points[ii][0]);
                    if (qvs->points[ii][0] <= S
						&& S < qvs->points[ii + 1][0])    
					break;
                }
                double q2 = qvs->points[ii + 1][1];
                double q1 = qvs->points[ii][1];
                double S2 = qvs->points[ii + 1][0];
                double S1 = qvs->points[ii][0];
                //dam_outflow = (q2 - q1) / (S2 - S1) * (S - S1) + q1;//m3s-1
                //rating_flow is Qrating in the flowchart
                rating_flow = (q2 - q1) / (S2 - S1) * (S - S1) + q1;//m3s-1

                if(output_mode == 0){ //passive
                    dam_outflow = rating_flow;
                }

                if(output_mode==1){
                    if(dam_storage + dam_input <= low_storage_threshold){
                        if(debug) printf("storage <= low_storage_threshold\n");
                    //is thhe pond storage too low, is it an emergency and we
                    //need to hurry to build storage
                        if(dam_storage + dam_input <= low_storage_threshold*0.8){
                            if(debug) printf("storage <= low_storage_threshold*.8\n");
                            dam_outflow = 1e-6;
                        }
                        if(dam_storage + dam_input > low_storage_threshold*0.8){
                            if(debug) printf("storage > low_storage_threshold*0.8\n");

                            //find ideal flow
                            double deltaS = low_storage_threshold-(dam_storage + dam_input);
                            double deltaQ = deltaS / 60.0;
                            double inflow = dam_input / 60.0;
                            double Qideal = inflow - deltaQ;
                            if(Qideal <=rating_flow){
                                dam_outflow = Qideal;
                            }
                            if(Qideal >rating_flow){
                                dam_outflow = rating_flow;
                            }   
                        }
                    } 
                    if(dam_storage + dam_input > low_storage_threshold){
                    //if the pond storage is greater than threslow but less than threshigh
                        if(dam_storage + dam_input < high_storage_threshold){
                            if(debug) printf("storage < high_storage_threshold\n");
                            double deltaS = (dam_storage + dam_input) - low_storage_threshold;
                            if(debug) printf("delS  : %f\n", deltaS);
                            double deltaQ = deltaS / 60.0;
                            if(debug) printf("delQ  : %f\n", deltaQ);
                            double inflow = dam_input / 60.0;
                            if(debug) printf("inflow  : %f\n", inflow);
                            double Qideal = factor_operating_window * (inflow + deltaQ);
                            if(debug) printf("Qideal  : %f\n", Qideal);
                            if(Qideal <=rating_flow){
                                dam_outflow = Qideal;
                            }
                            if(Qideal >rating_flow){
                                dam_outflow = rating_flow;
                            }
                        }
                        if(dam_storage + dam_input >= high_storage_threshold){
                            if(debug) printf("storage >= high_storage_threshold\n");
                            double deltaS = (dam_storage + dam_input) - low_storage_threshold;
                            double deltaQ = deltaS / 60.0;
                            double inflow = dam_input / 60.0;
                            double Qideal = inflow + deltaQ;
                            if(Qideal <=rating_flow){
                                dam_outflow = Qideal;
                            }
                            if(Qideal >rating_flow){
                                dam_outflow = rating_flow;
                            }
                        }
                    }
                }
                if(output_mode==2){
                    if(dam_storage + dam_input <= (max_storage_pond * 0.05)){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.05) && (dam_storage + dam_input <= (max_storage_pond * 0.1))){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.1) && (dam_storage + dam_input <= (max_storage_pond * 0.15))){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.15) && (dam_storage + dam_input <= (max_storage_pond * 0.2))){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.2) && (dam_storage + dam_input <= (max_storage_pond * 0.25))){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.25) && (dam_storage + dam_input <= (max_storage_pond * 0.3))){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.3) && (dam_storage + dam_input <= (max_storage_pond * 0.35))){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.35) && (dam_storage + dam_input <= (max_storage_pond * 0.4))){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.4) && (dam_storage + dam_input <= (max_storage_pond * 0.45))){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.45) && (dam_storage + dam_input <= (max_storage_pond * 0.5))){
                        dam_outflow = rating_flow * 0.05;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.5) && (dam_storage + dam_input <= (max_storage_pond * 0.55))){
                        dam_outflow = rating_flow * 0.1;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.55) && (dam_storage + dam_input <= (max_storage_pond * 0.6))){
                        dam_outflow = rating_flow * 0.3;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.6) && (dam_storage + dam_input <= (max_storage_pond * 0.65))){
                        dam_outflow = rating_flow * 0.3;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.65) && (dam_storage + dam_input <= (max_storage_pond * 0.7))){
                        dam_outflow = rating_flow * 0.5;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.7) && (dam_storage + dam_input <= (max_storage_pond * 0.75))){
                        dam_outflow = rating_flow * 0.5;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.75) && (dam_storage + dam_input <= (max_storage_pond * 0.8))){
                        dam_outflow = rating_flow * 0.7;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.8) && (dam_storage + dam_input <= (max_storage_pond * 0.85))){
                        dam_outflow = rating_flow * 0.7;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.85) && (dam_storage + dam_input <= (max_storage_pond * 0.9))){
                        dam_outflow = rating_flow * 0.9;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.9) && (dam_storage + dam_input <= (max_storage_pond * 0.95))){
                        dam_outflow = rating_flow * 0.9;
                    }
                    if(dam_storage + dam_input > (max_storage_pond * 0.95)){
                        dam_outflow = rating_flow;
                    }                
                }
            
            //calculate dam outflow
            //dam_outflow = get_discharge_from_storage(dam_storage + dam_input); //m3s-1
            if(debug) printf("dam outflow in m3/s: %f\n", dam_outflow);

            ans[STATE_DISCHARGE] = dam_outflow - q;
            
            //storage output
            //dam_output =get_storage_from_discharge(dam_outflow);
            dam_output = dam_outflow * 60.0; //m3s-1 to m3
            if(debug) printf("dam output storage in m3: %f\n", dam_output);
            if(debug) printf("step2 dam storage in m3: %f\n", dam_storage + dam_input - dam_output);

            //dam_storage = max(dam_storage,0);
            //if(debug) printf("step3 dam storage in m3: %f\n", dam_storage);
            //ans[STATE_DAM_STORAGE] = dam_storage;
            ans[STATE_DAM_STORAGE] = dam_input - dam_output;

            if(debug) printf("\n");

            
        }

        // if (forcing_values[0]>1 && ratio<1) {
        //     printf("time: %f\n", t);
        //     printf(" rain in mm/hour: %f\n", forcing_values[0]);
        //     printf(" area hill, area basin, area ratio: %f %f %f\n", A_h,A_i,ratio);
        //     MPI_Abort(MPI_COMM_WORLD, 1);
        // }
    }
}


//Type 193
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
//The numbering is:	0   1   2   3  4    5    6   7
//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g
//The numbering is:        0      1        2     3  4   5
void routing_runoff2(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double lambda_1 = global_params[1];

	double A_i = params[0];
	double L_i = params[1];
	double A_h = params[2];
	double invtau = params[3];
	//printf("invtau: %f\n", invtau);
    double q = y_i[0];		                                        // [m^3/s]

    double runoff = forcing_values[0]*(0.001 / 60.0);               //(mm/hr ->m/min)


    //in this model, all the water is routed downstream
    ans[0] = -q + (runoff * A_h / 60.0); //m3/min to m3/s
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    
	
    //Hillslope
    //ans[1] = q_rp - q_pl;

    //Sub-surface
   // ans[2] = q_ra - q_al - e_a;

    //Accumulated precip
    //ans[3] = q_rp + q_ra;


}

//Type 194
//Order of parameters: A_i,L_i,A_h,k2,k3,invtau,c_1,c_2
//The numbering is:	0   1   2   3  4    5    6   7
//Order of global_params: v_r,lambda_1,lambda_2,RC,v_h,v_g
//The numbering is:        0      1        2     3  4   5
void routing_runoff1(double t, const double * const y_i, unsigned int dim, const double * const y_p, unsigned short num_parents, unsigned int max_dim, const double * const global_params, const double * const params, const double * const forcing_values, const QVSData * const qvs, int state, void* user, double *ans)
{
    unsigned short i;

    double lambda_1 = global_params[1];

	double A_i = params[0];
	double L_i = params[1];
	double A_h = params[2];
	double invtau = params[3];
	//printf("invtau: %f\n", invtau);
    double q = y_i[0];		                                        // [m^3/s]

    double runoff = forcing_values[0]*(0.001 / 60.0);               //(mm/hr ->m/min)



    ans[0] = -q + (runoff * A_h / 60.0); //m3/min to m3/s
    for (i = 0; i<num_parents; i++)
        ans[0] += y_p[i * dim];
    ans[0] = invtau * pow(q, lambda_1) * ans[0];
	
    //Hillslope
    //ans[1] = q_rp - q_pl;

    //Sub-surface
   // ans[2] = q_ra - q_al - e_a;

    //Accumulated precip
    //ans[3] = q_rp + q_ra;


}


