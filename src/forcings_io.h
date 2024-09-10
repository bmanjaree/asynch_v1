#ifndef RAINFALL_H
#define RAINFALL_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

//

#include "structs.h"


int Create_Rain_Data_Par(
    Link *sys, unsigned int N,
    Link **my_sys, unsigned int my_N,
    const GlobalVars * const globals,
    int* assignments,
    char strfilename[],
    unsigned int first, unsigned int last,
    double t_0, double increment,
    Forcing* forcing, const Lookup * const id_to_loc, unsigned int max_files, unsigned int forcing_idx);

int Create_Rain_Data_GZ(
    Link *sys, unsigned int N,
    Link **my_sys,unsigned int my_N,
    const GlobalVars * const globals, int* assignments, char strfilename[], unsigned int first, unsigned int last, double t_0, double increment, Forcing* forcing, const Lookup * const id_to_loc, unsigned int max_files, unsigned int forcing_idx);

int Create_Rain_Data_Grid(
    Link *sys, unsigned int N,
    Link **my_sys, unsigned int my_N,
    const GlobalVars * const globals, int* assignments, char strfilename[], unsigned int first, unsigned int last, double t_0, double increment, Forcing* forcing, const Lookup * const id_to_loc, unsigned int max_files, unsigned int forcing_idx);

int Create_Rain_Database(
    Link *sys, unsigned int N,
    Link **my_sys, unsigned int my_N,
    const GlobalVars * const globals, int* assignments, ConnData *conninfo, unsigned int first, unsigned int last, Forcing* forcing, const Lookup * const id_to_loc, double maxtime, unsigned int forcing_idx);

//void SetRain0(Link* sys, unsigned int my_N, double maxtime, unsigned int* my_sys, const GlobalVars * const globals, Forcing* forcing, unsigned int forcing_idx);

double CreateForcing_Monthly(
    Link *sys, unsigned int N,
    Link **my_sys, unsigned int my_N,
    const GlobalVars * const globals, TimeSerie* global_forcings, unsigned int forcing_idx, struct tm *current_time, time_t first_time, time_t last_time, double t_0);

int Create_Rain_Database_Irregular(
    Link *sys, unsigned int N,
    Link **my_sys, unsigned int my_N,
    const GlobalVars * const globals, int* assignments, ConnData *conninfo, unsigned int first, unsigned int last, Forcing* forcing, const Lookup * const id_to_loc, double maxtime, unsigned int forcing_idx);

#endif

