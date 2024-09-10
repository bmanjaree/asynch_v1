#if !defined(ASYNCH_DB_H)
#define ASYNCH_DB_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#include <structs.h>
#include <libpq_fwd.h>

//PostgreSQL Database Methods
void ConnData_Init(ConnData* const conn, const char* connstring);
void ConnData_Free(ConnData* const conninfo);

int ConnectPGDB(ConnData* conninfo);
void DisconnectPGDB(ConnData* conninfo);

int CheckResError(PGresult* res,char* event);
int CheckResState(PGresult* res,short int error_code);
void CheckConnConnection(ConnData* conninfo);

#endif //ASYNCH_DB_H
