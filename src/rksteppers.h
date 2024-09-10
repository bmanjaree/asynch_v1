#if !defined(RKSTEPPERS_H)
#define RKSTEPPERS_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdbool.h>

#include <structs.h>

//Copies contents of the vectors full_k with dim entries into the vectors k with num_dense entries.
void store_k(const double * const full_k, unsigned int num_dof, double *k, unsigned int num_stages, const unsigned int * const dense_indices, unsigned int num_dense);

double InitialStepSize(double t, Link* link_i, const GlobalVars * const globals, Workspace* workspace);

// Steppers methods
int ExplicitRKSolver(Link* link_i, GlobalVars* globals, int* assignments, bool print_flag, FILE* outputfile, ConnData* conninfo, Forcing* forcings, Workspace* workspace);
int ExplicitRKIndex1SolverDam(Link* link_i, GlobalVars* globals, int* assignments, bool print_flag, FILE* outputfile, ConnData* conninfo, Forcing* forcings, Workspace* workspace);
int ExplicitRKIndex1Solver(Link* link_i, GlobalVars* globals, int* assignments, bool print_flag, FILE* outputfile, ConnData* conninfo, Forcing* forcings, Workspace* workspace);
int ExplicitRKSolverDiscont(Link* link_i, GlobalVars* globals, int* assignments, bool print_flag, FILE* outputfile, ConnData* conninfo, Forcing* forcings, Workspace* workspace);

//Forced solution methods
int ForcedSolutionSolver(Link* link_i, GlobalVars* globals, int* assignments, bool print_flag, FILE* outputfile, ConnData* conninfo, Forcing* forcings, Workspace* workspace);

#endif // !defined(RKSTEPPERS_H)
