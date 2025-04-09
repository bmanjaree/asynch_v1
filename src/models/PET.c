#include <stdio.h>
#include <math.h>
#define M_PI 3.14159265358979323846

/***********************************************************************************************************************************
* Put equations for PET estimation based on other forcings here.
* *********************************************************************************************************************************/

//Hamon PET estimation for daily average temperature using the CBM model to compute the number of 
// Requires the temperature input as daily average; latitude of hillslope; and day of year at time step.
double HamonPET(double temperature,\
                double latitude, \
                double doy) 
{
    double PET = 0;
    if(temperature > 0){
        
        //Calculate saturated water vapor
        double esat = 6.108 * exp((17.26939 * temperature)/(temperature + 237.3)); //saturation vapor pressure
        double Wt = 216.7 * esat / (temperature+273.3); //Saturated water vapor [g/m3]
        
        //Calculate the daily daylight hours (per 12 hours) using CBM model 
        double theta = 0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (doy - 186)));
        double phi = asin(0.39795 * cos(theta));
        double p = 0.8333; //the daylength coefficient
        double D = (24 - (24/M_PI)*acos((sin(p*M_PI/180) + sin(latitude*M_PI/180)*sin(phi))/(cos(latitude*M_PI/180)*cos(phi))))/12;
        
        //For artic cases
        if (isnan(D)) {
            D = 0;
            if (phi > 0 && latitude > 0) {
                D = 2;
            } else if (phi < 0 && latitude < 0) {
                D = 2;
            }
        }
        
        //Calculate daily Hamon PET in m/min
        PET = 1.6169 * pow(10,-6) * pow(D,2) * Wt * 60 / 1000;
    }
    
    return PET;
}