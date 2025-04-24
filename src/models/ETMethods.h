#if !defined(ETMETHODS_H)
#define ETMETHODS_H


#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


double HamonPET(double temperature, double latitude, double doy); 
double ETactual(double Emax, double s, double sw, double ss);


#endif //ETMETHODS_H