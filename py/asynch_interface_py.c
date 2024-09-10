#include "asynch_interface_py.h"


// Global variables
int my_rank = 0;
int np = 0;

//Utilities for mixing C and Python **************************************

void C_inc_ref(PyObject* obj)
{
    Py_INCREF(obj);
}

unsigned int* Allocate_CUINT_Array(unsigned int n)
{
    return (unsigned int*)malloc(n * sizeof(unsigned int));
}

void Free_PythonInterface(AsynchSolver* asynch)
{
    free((PythonInterface*)(asynch->ExternalInterface));
}


//Routines for the Python interface **************************************

AsynchSolver* Asynch_Init_py(int numprocs, int* ranks)
{
    //!!!! Assumes MPI_COMM_WORLD. This should create a new communicator and pass that along. !!!!
    AsynchSolver* asynch = Asynch_Init(MPI_COMM_WORLD, false);

    asynch->ExternalInterface = malloc(sizeof(PythonInterface));

    return asynch;
}

void Asynch_Set_System_State_py(AsynchSolver* asynch, double t_0, double* values)
{
    Asynch_Set_System_State(asynch, t_0, values);
}

//Custom Model routines *******************************************


int Asynch_Custom_Model_py(AsynchSolver *asynch, AsynchModel *model, PyObject *lib)
{
    //Setup the python interface data
    PythonInterface* external = (PythonInterface*)asynch->ExternalInterface;
    external->lib = lib;
    external->SetParamSizes_func = model->set_param_sizes;
    external->Routines_func = model->routines;

    //Set the functions
    int val = Asynch_Custom_Model(asynch, model);

    return val;
}

//Called by solvers
void SetParamSizes_py(GlobalVars* GlobalVars, void* user)
{
    PythonInterface* external = (PythonInterface*)user;
    external->SetParamSizes_func(GlobalVars, external->lib);
}


void InitRoutines_py(Link* link, unsigned int type, unsigned int exp_imp, unsigned short int dam, void* user)
{
    PythonInterface* external = (PythonInterface*)user;
    external->Routines_func(link, type, exp_imp, dam, external->lib);
}

void Asynch_Copy_Local_OutputUser_Data_py(AsynchSolver* asynch, unsigned int location, PyObject* source, unsigned int size)
{
    Asynch_Copy_Local_OutputUser_Data(asynch, location, (void*)source, size);
}
