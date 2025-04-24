#include <stdio.h>
#include <math.h>


/***********************************************************************************************************************************
* Put equations for soil temperature here
* *********************************************************************************************************************************/
#include <stdio.h>
#include <math.h>

// Soil Temperature calculation at daily scale 
// Requires: 
//      previous daily air temperature [c]
//      Previous daily ground temperature [c]
//      Depth of snow [m]

double soiltemp(double Tair,double Tz,double Ds) 
{   
    //Parameters from Rankinen et al 2002.
    double Cs = 1e6;    // J/m3/C
    double Kt = 0.516;   // W/m/C
    double Cice = 8.93e6;    // J/m3/C
    double fs = -2.7;    // m-1
    double Zs = 3.5 / 100;    // cm to m, middle of 0-7 cm (level 1)
    double delta_t = 3600*24;    // 1 day in seconds+
    double CA = Cs + Cice; // 
    double f = delta_t * Kt / (CA * pow(2 * Zs,2)); //intermediate factors
    
    //Soil temperature change
    double T_star = Tz + f * (Tair - Tz);
    double Tz_out = T_star * exp(-fs * Ds);
    
    return Tz_out;
}