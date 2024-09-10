# API for Python

Most of the documentation about the user interface routines applies to the Python interface. Naturally, some specifics concerning Python come into play.

The module `asynch_interface` for Python can be found in the directory `py` in the ASYNCH root directory. The module can be loaded a in Python script with a command such as:

```Python
from py.asynch_interface import asynchsolver
```

The routines in the Python interface can be found in Section [user interface routines](user_interface_routines.md). These routines are predominantly the same as those from the C interface, with the exception of the initialization and finalization routines (Sections [solver_initialization]() and [asynch_free](), respectively).

To demonstrate the usage of these routines in Python, consider the following code, which is analogous to the program given in Section [interface usage]():

```Python
#Include the ASYNCH Python interface
from py.asynch_interface import asynchsolver
import sys

#Prepare system
asynch = asynchsolver()
asynch.Parse_GBL(sys.argv[1])
asynch.Load_System()

#Prepare outputs
asynch.Prepare_Temp_Files()
asynch.Write_Current_Step()
asynch.Prepare_Peakflow_Output()
asynch.Prepare_Output()

#Advance solver
asynch.Advance(1)

#Take a snapshot
asynch.Take_System_Snapshot(None)

#Create output files
asynch.Create_Output(None)
asynch.Create_Peakflows_Output()

#Cleanup
asynch.Delete_Temporary_Files()
```

The contents of the source code `asynchdist.py` is essentially the program above. The module `sys` is loaded to read the filename of the global file from the command line. An instance of the `asynchsolver` class is created and stored as `asynch`. From this object, all the routines from the interface are called. As is traditional in Python, clean up of the `asynchsolver` instance is performed automatically.

Creating custom models and outputs with the Python interface is similar to that of the C interface. See Sections [custom outputs]() and [custom models](). Any routines passed into the Python interface written in Python must be decorated as an appropriate function type. The available data types are:

| Type | Description |
| --- | --- |
| `ASYNCH_F_DATATYPE` | ODE |
| `ASYNCH_RKSOLVER_DATATYPE` | Runge-Kutta Solver |
| `ASYNCH_CONSISTENCY_DATATYPE` | System Consistency |
| `ASYNCH_ALG_DATATYPE` | Algebraic Equations |
| `ASYNCH_STATECHECK_DATATYPE` | State Check |
| `ASYNCH_OUTPUT_INT_DATATYPE` | Output Time Series with Integer Values |
| `ASYNCH_OUTPUT_DOUBLE_DATATYPE` | Output Time Series with Double Precision Values |
| `ASYNCH_PEAKOUTPUT_DATATYPE` | Output Peakflow Routine |
| `ASYNCH_SETPARAMSIZES_DATATYPE` | SetParmSizes Routine for Custom Models |
| `ASYNCH_CONVERT_DATATYPE` | ConvertParams Routine for Custom Models |
| `ASYNCH_ROUTINES_DATATYPE` | InitRoutines Routine for Custom Models |
| `ASYNCH_PRECALCULATIONS_DATATYPE` | Precalculations Routine for Custom Models |
| `ASYNCH_INITIALIZEEQS_DATATYPE` | ReadInitData Routine for Custom Models |

Decorating a routine in Python can be done with the `@` symbol. For instance, the following is the definition of a function for converting units of state parameters:

```Python
@ASYNCH_CONVERT_DATATYPE
def ConvertParams_MyModel(params,model_type,lib_ptr):
	params.contents.ve[1] *= 1000
	params.contents.ve[2] *= 1e6
```

A function decorated with an ASYNCH function type cannot be called as a regular Python function. It should only be called from a routine written in C.

The routines for creating a custom model have one additional argument available. This is a library object to the ASYNCH C routines. Such routines may be used by and given to Python interface functions. For instance, consider this example routine used for the `ReadInitData` routine:

```Python
@ASYNCH_ROUTINES_DATATYPE
def InitRoutines_MyModel(link_p,model_type,\
    exp_imp,dam,lib_p):
	lib = cast(lib_p,py_object).value
	link = link_p.contents

	if link.res:
	  link.f = \
	    LinearHillslope_Reservoirs_MyModel
	  link.RKSolver = cast(\
	    lib.ForcedSolutionSolver,
	    ASYNCH_RKSOLVER_DATATYPE)
	else:
	  link.f = LinearHillslope_MyModel
	  link.RKSolver = cast(\
	  lib.ExplicitRKSolver,
	    ASYNCH_RKSOLVER_DATATYPE)
	  link.alg = cast(None,ASYNCH_ALG_DATATYPE)
	  link.state_check = cast(None,\
			   ASYNCH_STATECHECK_DATATYPE)
	  link.CheckConsistency = cast(\
	    lib.CheckConsistency_Nonzero_AllStates_q,
	    ASYNCH_CONSISTENCY_DATATYPE)
```

In this sample, the Runge-Kutta methods are set to functions defined in the ASYNCH library. Also note that if a routine does not need to be set (here, the algebraic and state check routines), then the routine is set to `None` casted as the appropriate function data type. The `SetParamSizes` routine for creating custom models requires creating an array. This can be done with a call to the function `Allocate_CUINT_Array` from the ASYNCH library. It requires only one argument, the size of the array, and its return value can be directly set to the member `dense_indices`.
