#if !defined(ASYNCH_ADVANCE_H)
#define ASYNCH_ADVANCE_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include <stdio.h>
#include <stdbool.h>

#include <globals.h>
#include <structs.h>

void Advance(
    Link *sys, unsigned int N,
    Link **my_sys, unsigned int my_N,
    GlobalVars* globals,
    int* assignments, short int* getting, unsigned int* res_list, unsigned int res_size, const Lookup * const id_to_loc,
    Workspace* workspace,
    Forcing* forcings,
    ConnData* db_connections,
    TransData* my_data,
    int print_level,
    FILE* outputfile);

#endif //ASYNCH_ADVANCE_H
