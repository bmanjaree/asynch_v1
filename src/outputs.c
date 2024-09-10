#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#include <asynch_interface.h>
#include <outputs.h>

void SetDefaultOutputFunctions(char* outputname, Output *output, unsigned int* states_used, unsigned int* num_states_used)
{
    assert(outputname != NULL);

    if (strcmp(outputname, "Time") == 0)
    {
        output->type = ASYNCH_DOUBLE;
        output->callback.out_double = &Output_Time;
    }
    else if (strcmp(outputname, "State0") == 0)
    {
        states_used[(*num_states_used)++] = 0;
        output->type = ASYNCH_FLOAT;
        output->callback.out_float = &Output_State0;
    }
    else if (strcmp(outputname, "State1") == 0)
    {
        states_used[(*num_states_used)++] = 1;
        output->type = ASYNCH_FLOAT;
        output->callback.out_float = &Output_State1;
    }
    else if (strcmp(outputname, "State2") == 0)
    {
        states_used[(*num_states_used)++] = 2;
        output->type = ASYNCH_FLOAT;
        output->callback.out_float = &Output_State2;
    }
    else if (strcmp(outputname, "State3") == 0)
    {
        states_used[(*num_states_used)++] = 3;
        output->type = ASYNCH_FLOAT;
        output->callback.out_float = &Output_State3;
    }
    else if (strcmp(outputname, "State4") == 0)
    {
        states_used[(*num_states_used)++] = 4;
        output->type = ASYNCH_FLOAT;
        output->callback.out_float = &Output_State4;
    }
    else if (strcmp(outputname, "State5") == 0)
    {
        states_used[(*num_states_used)++] = 5;
        output->type = ASYNCH_FLOAT;
        output->callback.out_float = &Output_State5;
    }
    else if (strcmp(outputname, "State6") == 0)
    {
        states_used[(*num_states_used)++] = 6;
        output->type = ASYNCH_FLOAT;
        output->callback.out_float = &Output_State6;
    }
    else if (strcmp(outputname, "State7") == 0)
    {
        states_used[(*num_states_used)++] = 7;
        output->type = ASYNCH_FLOAT;
        output->callback.out_float = &Output_State7;
    }
    else if (strcmp(outputname, "TimeI") == 0)
    {
        output->type = ASYNCH_INT;
        output->callback.out_int = &Output_Time_Int;
    }
    else if(strcmp(outputname,"LinkID") == 0)
    {
        output->type = ASYNCH_INT;
        output->callback.out_int = &Output_LinkID;
    }
    else	//Must be a user defined output
    {
        output->type = -1;
        return;
        //printf("[%i]: Error setting output: bad output function name %s.\n",my_rank,outputname);
        //MPI_Abort(MPI_COMM_WORLD,1);
    }

    output->size = GetByteSize(output->type);
    output->specifier = GetSpecifier(output->type);
}

void SetPeakflowOutputFunctions(char* outputname, PeakflowOutputCallback **peak_output)
{
    if (strcmp(outputname, "Classic") == 0)
        *peak_output = &OutputPeakflow_Classic_Format;
    else if (strcmp(outputname, "Forecast") == 0)
        *peak_output = &OutputPeakflow_Forecast_Format;
    else	//Must be a user defined output
    {
        *peak_output = NULL;
        //printf("[%i]: Bad peakflow function %s.\n",my_rank,outputname);
        //MPI_Abort(MPI_COMM_WORLD,1);
    }
}

//ASYNCH_CHAR = 0
//ASYNCH_SHORT = 1
//ASYNCH_INT = 2
//ASYNCH_FLOAT = 3
//ASYNCH_DOUBLE = 4
short int GetByteSize(enum AsynchTypes type)
{
    static const short int data_sizes[] = { sizeof(char), sizeof(short int), sizeof(int), sizeof(float), sizeof(double) };
    return data_sizes[type];
}

const char* GetSpecifier(enum AsynchTypes type)
{
    static const char *data_specifiers[5] = { "%c", "%hi", "%i", "%.6e", "%.12e" };
    //memcpy(specifier, data_specifiers[type], strlen(data_specifiers[type]) + 1);

    return data_specifiers[type];
}

unsigned int CalcTotalOutputSize(GlobalVars* GlobalVars)
{
    unsigned total = 0;
    for (unsigned int i = 0; i < GlobalVars->num_outputs; i++)
        total += GlobalVars->outputs[i].size;

    return total;
}

//Returns 0 if an output is not set, 1 if all outputs are good to go.
bool AreOutputsSet(GlobalVars* GlobalVars)
{
    unsigned int i;

    for (i = 0; i < GlobalVars->num_outputs; i++)
    {
        switch (GlobalVars->outputs[i].type)
        {
        case ASYNCH_BAD_TYPE:
            return false;
        case ASYNCH_INT:
            if (!(GlobalVars->outputs[i].callback.out_int))
                return false;
            break;
        case ASYNCH_DOUBLE:
            if (!(GlobalVars->outputs[i].callback.out_double))
                return false;
            break;
        case ASYNCH_FLOAT:
            if (!(GlobalVars->outputs[i].callback.out_float))
                return false;
            break;
        default:
            printf("[%i]: Error: Bad output type (%i) encountered while checking if outputs set.\n", my_rank, GlobalVars->outputs[i].type);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    return true;
}

//Returns 0 if peakflow output is not set, 1 if it is.
int PeakflowOutputsSet(GlobalVars* GlobalVars)
{
    if (GlobalVars->peakflow_function_name && !GlobalVars->peakflow_output)	return 0;
    return 1;
}


//Output functions ***************************************************************************************************

int Output_LinkID(unsigned int id, double t, double *y, unsigned int num_dof)
{
    return id;
}

double Output_Time(unsigned int id, double t, double *y, unsigned int num_dof)
{
    return t;
}

float Output_State0(unsigned int id, double t, double *y, unsigned int num_dof)
{
    return (float)y[0];
}

float Output_State1(unsigned int id, double t, double *y, unsigned int num_dof)
{
    return (float)y[1];
}

float Output_State2(unsigned int id, double t, double *y, unsigned int num_dof)
{
    return (float)y[2];
}

float Output_State3(unsigned int id, double t, double *y, unsigned int num_dof)
{
    return (float)y[3];
}

float Output_State4(unsigned int id, double t, double *y, unsigned int num_dof)
{
    return (float)y[4];
}

float Output_State5(unsigned int id, double t, double *y, unsigned int num_dof)
{
    return (float)y[5];
}

float Output_State6(unsigned int id, double t, double *y, unsigned int num_dof)
{
    return (float)y[6];
}

float Output_State7(unsigned int id, double t, double *y, unsigned int num_dof)
{
    return (float)y[7];
}

int Output_Time_Int(unsigned int id, double t, double *y, unsigned int num_dof)
{
    return (int)(round(t) + 0.1);
}


//Peakflow output functions***********************************************************************************

void OutputPeakflow_Classic_Format(unsigned int ID, double peak_time, double *peak_value, double *params, double *global_params, double conversion, unsigned int area_idx, void* user, char* buffer)
{
    sprintf(buffer, "%u %.4f %.8f %.8f\n", ID, conversion * params[area_idx], peak_time, peak_value[0]);
}

void OutputPeakflow_Forecast_Format(unsigned int ID, double peak_time, double *peak_value, double *params, double *global_params, double conversion, unsigned int area_idx, void* user, char* buffer)
{
    unsigned int offset = *(unsigned int*)user;
    //sprintf(buffer,"%u %u %.12e %u NULL\n",ID,offset + (unsigned int)(peak_time*60 + .1),peak_value.ve[0],offset);
    sprintf(buffer, "%u %u %.12e %u\n", ID, offset + (unsigned int)(peak_time * 60 + .1), peak_value[0], offset);
}

