#ifndef ASYNCH_INTERFACE_PY_H
#define ASYNCH_INTERFACE_PY_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#ifdef _DEBUG
#undef _DEBUG
#include <Python.h>
#define _DEBUG
#else
#include <Python.h>
#endif

#include <asynch_interface.h>

typedef struct PythonInterface
{
	PyObject* lib;
	void (*SetParamSizes_func)(GlobalVars*,PyObject*);
	void (*Routines_func)(Link*,unsigned int,unsigned int,unsigned short int,void*);
} PythonInterface;

void C_inc_ref(PyObject* obj);
unsigned int* Allocate_CUINT_Array(unsigned int n);
void Free_PythonInterface(AsynchSolver* asynch);

AsynchSolver* Asynch_Init_py(int numprocs,int* ranks);
void Asynch_Set_System_State_py(AsynchSolver* asynch,double t_0,double* values);

int Asynch_Custom_Model_py(AsynchSolver* asynch, AsynchModel* model, PyObject* lib);

void SetParamSizes_py(GlobalVars* GlobalVars,void* user);
void InitRoutines_py(Link* link,unsigned int type,unsigned int exp_imp,unsigned short int dam,void* user);
void Asynch_Copy_Local_OutputUser_Data_py(AsynchSolver* asynch,unsigned int location,PyObject* source,unsigned int size);

#endif

