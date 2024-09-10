#if !defined(ASYNCH_CONSTANTS_H)
#define ASYNCH_CONSTANTS_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

//Constants
#define ASYNCH_MAX_DB_CONNECTIONS 20

#define ASYNCH_DB_LOC_TOPO 0
#define ASYNCH_DB_LOC_PARAMS 1
#define ASYNCH_DB_LOC_INIT 2
#define ASYNCH_DB_LOC_QVS 3
#define ASYNCH_DB_LOC_RSV 4
#define ASYNCH_DB_LOC_HYDROSAVE 5
#define ASYNCH_DB_LOC_PEAKSAVE 6
#define ASYNCH_DB_LOC_HYDRO_OUTPUT 7
#define ASYNCH_DB_LOC_PEAK_OUTPUT 8
#define ASYNCH_DB_LOC_SNAPSHOT_OUTPUT 9
#define ASYNCH_DB_LOC_FORCING_START 10

#define ASYNCH_MAX_QUERIES 5

#define ASYNCH_MAX_NUM_FORCINGS 12
#define ASYNCH_MAX_TIMESTAMP_LENGTH 12
#define ASYNCH_MAX_PATH_LENGTH 1024

#define ASYNCH_MAX_LINE_LENGTH 1024
#define ASYNCH_MAX_SYMBOL_LENGTH 64
#define ASYNCH_MAX_QUERY_LENGTH 16384
#define ASYNCH_MAX_CONNSTRING_LENGTH 1024

#define ASYNCH_MAX_SOLVER_STAGES 8      //!< Maximum number of stages in RK solvers

#define ASYNCH_MAX_DIM 256              //!< Maximum number of Degree of Freedom in the outputs

#define ASYNCH_LINK_MAX_PARENTS 8

#endif //ASYNCH_CONSTANTS_H