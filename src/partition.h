#ifndef PARTITION_H
#define PARTITION_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <structs.h>
#include <system.h>
#include <sort.h>
//#include <metis.h>

extern int my_rank;
extern int np;

int* Partition_System_By_Leaves(Link *sys, unsigned int N, Link **leaves, unsigned int numleaves, Link ***my_sys, unsigned int *my_N, TransData *my_data, short int *getting);
int* Partition_System_By_Leaves_2(Link *sys, unsigned int N, Link **leaves, unsigned int numleaves, Link ***my_sys, unsigned int * my_N, TransData *my_data, short int *getting);

#if defined(HAVE_METIS)

int* Partition_METIS_ByEqs(Link* sys, unsigned int N, Link** leaves, unsigned int numleaves, Link** my_sys, unsigned int* my_N, TransData* my_data, short int *getting);

#endif

#endif
