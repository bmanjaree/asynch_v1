#ifndef IO_H
#define IO_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include <stdbool.h>
#include <stdio.h>

#include <structs.h>

extern int np;
extern int my_rank;

void OutputFunc_Init(unsigned short hydros_loc_flag, unsigned short peaks_loc_flag, unsigned short dump_loc_flag, OutputFunc* output_func);
void WriteValue(FILE* outputfile, const char* specifier, char* data_storage, short int data_type, char* delim);
unsigned int WriteStep(Output *output, unsigned int num_outputs, FILE* outputfile, unsigned int id, double t, double *y, unsigned int num_dof, long int* pos_offset);
unsigned int CatBinaryToString(char* submission, const char* specifier, void* data_storage, short int data_type, char* delim);

#endif
