#ifndef SYSTEM_H
#define SYSTEM_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <structs.h>

extern int my_rank;
extern int np;

/// Creates a list to hold the data for an ODE.
///
/// \param t0: the initial time.
/// \param y0: Vector of the initial data [num_dof]
/// \param num_dof: the number of degree of freedom of the ODE.
/// \param num_dense_dof
/// \param num_stages: the number of stages in the RKMethod.
/// \param list_length: the maximum number of steps to store in the list.
void Init_List(RKSolutionList* list, double t0, double *y0, unsigned int num_dof, unsigned int num_dense_dof, unsigned short int num_stages, unsigned int list_length);

//Destructors
void Destroy_Link(Link* link_i, int rkd_flag, Forcing* forcings, GlobalVars* GlobalVars);
void Destroy_ForcingData(TimeSerie* forcing_buff);
void Destroy_RKMethod(RKMethod* method);
void Destroy_ErrorData(ErrorData* error);
void Destroy_List(RKSolutionList* list);
void Destroy_UnivVars(GlobalVars* GlobalVars);

//Discontinuity list
unsigned int Insert_Discontinuity(double time, unsigned int start, unsigned int end, unsigned int* count, unsigned int size, double* array, unsigned int id);
void Insert_SendDiscontinuity(double time, unsigned int order, unsigned int* count, unsigned int size, double* array, unsigned int*order_array, unsigned int id);

//RKSolution list
RKSolutionNode* New_Step(RKSolutionList* list);
void Undo_Step(RKSolutionList* list);
void Remove_Head_Node(RKSolutionList* list);

//Workspace methods
void Create_Workspace(Workspace *workspace, unsigned int dim, unsigned short num_stages, unsigned short max_parents);
void Destroy_Workspace(Workspace* workspace, unsigned short int s, unsigned short int max_parents);

#endif