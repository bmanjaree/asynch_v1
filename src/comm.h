#if !defined(ASYNCH_COMM_H)
#define ASYNCH_COMM_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#include <structs.h>
#include <globals.h>

#include <system.h>
#include <sort.h>
#include <libpq_fwd.h>

#define ASYNCH_MAX_NUMBER_OF_PROCESS 256

//MPI Related Methods
void Transfer_Data(TransData* my_data,Link* sys,int* assignments,GlobalVars* GlobalVars);
void Transfer_Data_Finish(TransData* my_data,Link* sys,int* assignments,GlobalVars* GlobalVars);
void Exchange_InitState_At_Forced(Link* system, unsigned int N, int* assignments, short int* getting, unsigned int* res_list, unsigned int res_size, const Lookup * const id_to_loc, GlobalVars* globals);
TransData* Initialize_TransData();
void Flush_TransData(TransData* data);
void TransData_Free(TransData* data);

#endif //!defined(ASYNCH_COMM_H)
