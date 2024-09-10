#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_POSTGRESQL)
#include <libpq-fe.h>
#endif

#include <comm.h>
#include <minmax.h>

#if defined(HAVE_POSTGRESQL)


//Disables postgresql notices
void Quiet(void *arg, const char *message)
{
    return;
}

//Create a ConnData object
void ConnData_Init(ConnData* const conndata, const char* connstring)
{
    memset(conndata, 0, sizeof(ConnData));
    strcpy(conndata->connectinfo, connstring);
}

//Destroy a ConnData object
void ConnData_Free(ConnData* const conninfo)
{
    if (conninfo)
    {
        //if(my_rank == 0 && conninfo->conn != NULL)	PQfinish(conninfo->conn);
        if (conninfo->conn && PQstatus(conninfo->conn) == CONNECTION_OK)
            PQfinish(conninfo->conn);
        //free(conninfo->query);
        for (unsigned int i = 0; i < conninfo->num_queries; i++)
            free(conninfo->queries[i]);
        //if(conninfo->queries)	free(conninfo->queries);
        //free(conninfo->connectinfo);
        //free(conninfo);
    }
}

//Connect to the database with information stored in connectinfo
int ConnectPGDB(ConnData* conninfo)
{
    conninfo->conn = PQconnectdb(conninfo->connectinfo);
    if (PQstatus(conninfo->conn) == CONNECTION_BAD)
    {
        printf("[%i]: Error: Unable to connect to the database.\n", my_rank);
        return 1;
    }
    PQsetNoticeProcessor(conninfo->conn, Quiet, NULL);	//Disable annoying notices
    return 0;
}

//Disconnect from the database
void DisconnectPGDB(ConnData* conninfo)
{
    if (PQstatus(conninfo->conn) == CONNECTION_OK)
    {
        PQfinish(conninfo->conn);
        conninfo->conn = NULL;
    }
}

//Check if an error related to an sql query occurred.
int CheckResError(PGresult* res, char* event)
{
    short int status = PQresultStatus(res);
    if (!(status == PGRES_COMMAND_OK || status == PGRES_TUPLES_OK))
    {
        printf("[%i]: SQL error encountered while %s. %hi\n", my_rank, event, status);
        printf("[%i]: %s\n", my_rank, PQresultErrorMessage(res));

        return 1;
    }

    return 0;
}

//Check if a query returned a certain value
int CheckResState(PGresult* res, short int error_code)
{
    short int status = PQresultStatus(res);
    if (status == error_code)	return 0;
    else
    {
        printf("[%i]: Error: did not get error code %hi. Got %hi.\n", my_rank, error_code, status);
        return 1;
    }
}

//Check if connection to SQL database is still good
void CheckConnConnection(ConnData* conninfo)
{
    if (PQstatus(conninfo->conn) == CONNECTION_BAD)
    {
        printf("[%i]: Connection to database lost. Attempting to reconnect...\n", my_rank);
        PQreset(conninfo->conn);
        printf("[%i]: Connection reestablished.\n", my_rank);
    }
}

#else

void ConnData_Init(ConnData* const conndata, const char* connstring) {}
void ConnData_Free(ConnData* const conninfo) {}

void CheckConnConnection(ConnData* conninfo) {}
int CheckResState(PGresult* res, short int error_code) {}
int CheckResError(PGresult* res, char* event) {}
void DisconnectPGDB(ConnData* conninfo) {}

void SwitchDB(ConnData* conninfo, char connectinfo[]) {}

int ConnectPGDB(ConnData* conninfo) {}

#endif



