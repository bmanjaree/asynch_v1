#ifndef FORCINGS_H
#define FORCINGS_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <structs.h>
#include <globals.h>

void Forcing_Init(Forcing* forcing);
void Forcing_Free(Forcing* forcing);

unsigned int PassesOther(Forcing* forcing, double maxtime, ConnData* conninfo);
unsigned int PassesBinaryFiles(Forcing* forcing, double maxtime, ConnData* conninfo);
unsigned int PassesIrregularBinaryFiles(Forcing* forcing, double maxtime, ConnData* conninfo);
unsigned int PassesDatabase(Forcing* forcing, double maxtime, ConnData* conninfo);
unsigned int PassesRecurring(Forcing* forcing, double maxtime, ConnData* conninfo);
unsigned int PassesDatabase_Irregular(Forcing* forcing, double maxtime, ConnData* conninfo);

double NextForcingOther(Link* sys, unsigned int N, Link **my_sys, unsigned int my_N, int* assignments, const GlobalVars * const globals, Forcing* forcing, ConnData* db_connections, const Lookup * const id_to_loc, unsigned int forcing_idx);
double NextForcingBinaryFiles(Link* sys, unsigned int N, Link **my_sys, unsigned int my_N, int* assignments, const GlobalVars * const globals, Forcing* forcing, ConnData* db_connections, const Lookup * const id_to_loc, unsigned int forcing_idx);
double NextForcingIrregularBinaryFiles(Link* sys, unsigned int N, Link **my_sys, unsigned int my_N, int* assignments, const GlobalVars * const globals, Forcing* forcing, ConnData* db_connections, const Lookup * const id_to_loc, unsigned int forcing_idx);
double NextForcingGZBinaryFiles(Link* sys, unsigned int N, Link **my_sys, unsigned int my_N, int* assignments, const GlobalVars * const globals, Forcing* forcing, ConnData* db_connections, const Lookup * const id_to_loc, unsigned int forcing_idx);
double NextForcingGridCell(Link* sys, unsigned int N, Link **my_sys, unsigned int my_N, int* assignments, const GlobalVars * const globals, Forcing* forcing, ConnData* db_connections, const Lookup * const id_to_loc, unsigned int forcing_idx);
double NextForcingDatabase(Link* sys, unsigned int N, Link **my_sys, unsigned int my_N, int* assignments, const GlobalVars * const globals, Forcing* forcing, ConnData* db_connections, const Lookup * const id_to_loc, unsigned int forcing_idx);
double NextForcingRecurring(Link* sys, unsigned int N, Link **my_sys, unsigned int my_N, int* assignments, const GlobalVars * const globals, Forcing* forcing, ConnData* db_connections, const Lookup * const id_to_loc, unsigned int forcing_idx);
double NextForcingDatabase_Irregular(Link* sys, unsigned int N, Link **my_sys, unsigned int my_N, int* assignments, const GlobalVars * const globals, Forcing* forcing, ConnData* db_connections, const Lookup * const id_to_loc, unsigned int forcing_idx);

#endif

