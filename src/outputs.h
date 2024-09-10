#ifndef MODELOUTPUTS_H
#define MODELOUTPUTS_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <structs.h>
#include <data_types.h>

extern int np;
extern int my_rank;

/// Set the built-in ouput functions
/// LinkID, Time, TimeI,  State0State1, .., State7
void SetDefaultOutputFunctions(char* outputname, Output *output, unsigned int* states_used, unsigned int* num_states_used);
void SetPeakflowOutputFunctions(char* outputname, PeakflowOutputCallback **peak_output);

/// Returns the number of bytes for the given type
short int GetByteSize(enum AsynchTypes type);

/// Returns the format specifier (see printf specification) for the given type 
const char* GetSpecifier(enum AsynchTypes type);

unsigned int CalcTotalOutputSize(GlobalVars* GlobalVars);
bool AreOutputsSet(GlobalVars* GlobalVars);
int PeakflowOutputsSet(GlobalVars* GlobalVars);

/// Built-in Output functions
int Output_LinkID(unsigned int id, double t, double *y, unsigned int num_dof);
double Output_Time(unsigned int id, double t, double *y, unsigned int num_dof);
float Output_State0(unsigned int id, double t, double *y, unsigned int num_dof);
float Output_State1(unsigned int id, double t, double *y, unsigned int num_dof);
float Output_State2(unsigned int id, double t, double *y, unsigned int num_dof);
float Output_State3(unsigned int id, double t, double *y, unsigned int num_dof);
float Output_State4(unsigned int id, double t, double *y, unsigned int num_dof);
float Output_State5(unsigned int id, double t, double *y, unsigned int num_dof);
float Output_State6(unsigned int id, double t, double *y, unsigned int num_dof);
float Output_State7(unsigned int id, double t, double *y, unsigned int num_dof);
int Output_Time_Int(unsigned int id, double t, double *y, unsigned int num_dof);

// Built-in Peakflow output functions
void OutputPeakflow_Classic_Format(unsigned int ID, double peak_time, double *peak_value, double *params, double *global_params, double conversion, unsigned int area_idx, void* user, char* buffer);
void OutputPeakflow_Forecast_Format(unsigned int ID, double peak_time, double *peak_value, double *params, double *global_params, double conversion, unsigned int area_idx, void* user, char* buffer);

#endif

