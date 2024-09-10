C API
=====

ASYNCH comes with a collection of routines to access specific details of the underlying solver structure. This allows a user to construct his or her own programs, while making calls to the solver in particular ways. This is useful for modifying how the solvers behave, and for producing more specialized outputs. A user can also specify his or her own models in a separate module. Currently, interface routines exist for the C/C++ and Python programming languages.

Typical Interface Usage
-----------------------

In this section, we provide an overview of how to use the API routines for running simulations. Regardless of the language, the general order of calling API routines is the same. Further routines may be added to increase flexibility.

A basic program that uses the ASYNCH solvers to perform calculations and takes advantage of all the features of a global file is the following (written in C):

.. code-block:: c

  //Make sure we finalize MPI
  void asynch_onexit(void)
  {
      int flag;
      MPI_Finalized(&flag);
      if (!flag)
          MPI_Finalize();
  }

  int main(int argc,char* argv[])
  {
    int res;
    char *global_filename = argv[1];

    //Initialize MPI stuff
    res = MPI_Init(&argc, &argv);
    if (res == MPI_SUCCESS)
        atexit(asynch_onexit);
    else
    {
        print_err("Failed to initialize MPI");
        exit(EXIT_FAILURE);
    }

    //Init asynch object and the river network
    AsynchSolver asynch;
    Asynch_Init(&asynch, MPI_COMM_WORLD);
    Asynch_Parse_GBL(&asynch, global_filename);
    Asynch_Load_Network(&asynch);
    Asynch_Partition_Network(&asynch);
    Asynch_Load_Network_Parameters(&asynch, 0);
    Asynch_Load_Dams(&asynch);
    Asynch_Load_Numerical_Error_Data(&asynch);
    Asynch_Initialize_Model(&asynch);
    Asynch_Load_Initial_Conditions(&asynch);
    Asynch_Load_Forcings(&asynch);
    Asynch_Load_Save_Lists(&asynch);
    Asynch_Finalize_Network(&asynch);
    Asynch_Calculate_Step_Sizes(&asynch);

    //Prepare output files
    Asynch_Prepare_Temp_Files(&asynch);
    Asynch_Write_Current_Step(&asynch);
    Asynch_Prepare_Peakflow_Output(&asynch);
    Asynch_Prepare_Output(&asynch);

    //Perform the calculations
    Asynch_Advance(&asynch, 1);

    //Create output files
    Asynch_Take_System_Snapshot(&asynch, NULL);
    Asynch_Create_Output(&asynch, NULL);
    Asynch_Create_Peakflows_Output(&asynch);

    //Clean up
    Asynch_Delete_Temporary_Files(&asynch);
    Asynch_Free(&asynch);

    return EXIT_SUCCESS;
  }

Details of each function call can be found in Section [sec: user interface routines].

The program begins by initializing the asynchsolver object. The MPI communicator consisting of all processes (:code:`MPI_COMM_WORLD`) is used for the calculations. Next, the global file specified as a command line argument to the program is parsed. Based upon the information specified by the global file, the different components of the network are constructed.

Next, the outputs for the program are initialized. The files for holding calculation results while the program runs is initialized (i.e. the temporary files). The initial values of the states are written to disk. Any initializations needed for the output peakflow and output time series sources are done.

With a call to :code:`Asynch_Advance`, the calculations are performed. Any time series results are temporarily stored in the temporary files.

Once the calculations are complete, any necessary outputs are created (snapshot, output time series, output peakflow data). This could include writing to files or database tables, depending upon the options selected in the global file.

Lastly, clean up routines are called. The temporary files are deleted. The asynchsolver object is also deleted from memory with a call to . Note that this routine does not need to be called if using the interface functions from a language that supports automatic garbage collection.

Of course, outputting more information might be useful (timing results, command line parameter checking, results printed to screen, etc), and additional features may be needed (outputting data to multiple sources, creating peakflow data over intervals of time, etc). However, this is the basic structure needed to perform simulations. The source files *asynchdist.c* and *asynchdist.py* are essentially duplicates of the above program, but with information printed to screen.

User Interface Routines
-----------------------

In this section, routines for operating the solver are described. These routines can be used to create an instance of an ASYNCH solver, and manipulate properties such as total simulation time, when data output occurs, etc... Creation of custom outputs is discussed in Section [sec: custom outputs] and creation of custom models is discussed in Section [sec: custom-models].

Solver Constructor / Destructor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: Asynch_Init

.. doxygenfunction:: Asynch_Free

Solver Initialization
~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: Asynch_Parse_GBL

.. doxygenfunction:: Asynch_Load_Network

.. doxygenfunction:: Asynch_Partition_Network

.. doxygenfunction:: Asynch_Load_Network_Parameters

.. doxygenfunction:: Asynch_Load_Dams

.. doxygenfunction:: Asynch_Load_Numerical_Error_Data

.. doxygenfunction:: Asynch_Initialize_Model

.. doxygenfunction:: Asynch_Load_Initial_Conditions

.. doxygenfunction:: Asynch_Load_Forcings

.. doxygenfunction:: Asynch_Load_Save_Lists

.. doxygenfunction:: Asynch_Finalize_Network

Integration
~~~~~~~~~~~

.. doxygenfunction:: Asynch_Calculate_Step_Sizes

.. doxygenfunction:: Asynch_Advance

Timeseries Output
~~~~~~~~~~~~~~~~~

.. doxygenfunction:: Asynch_Prepare_Output
.. doxygenfunction:: Asynch_Create_Output

.. doxygenfunction:: Asynch_Prepare_Peakflow_Output
.. doxygenfunction:: Asynch_Create_Peakflows_Output

.. doxygenfunction:: Asynch_Write_Current_Step

Temporary files
~~~~~~~~~~~~~~~

.. doxygenfunction:: Asynch_Prepare_Temp_Files

.. doxygenfunction:: Asynch_Delete_Temporary_Files

.. doxygenfunction:: Asynch_Set_Temp_Files

Snapshots
~~~~~~~~~

.. doxygenfunction:: Asynch_Take_System_Snapshot

.. doxygenfunction:: Asynch_Get_Snapshot_Output_Name
.. doxygenfunction:: Asynch_Set_Snapshot_Output_Name

Getters and Setters
~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: Asynch_Get_Model_Type
.. doxygenfunction:: Asynch_Set_Model_Type

.. doxygenfunction:: Asynch_Get_Total_Simulation_Duration
.. doxygenfunction:: Asynch_Set_Total_Simulation_Duration

.. doxygenfunction:: Asynch_Get_Num_Links
.. doxygenfunction:: Asynch_Get_Links

.. doxygenfunction:: Asynch_Get_Num_Links_Proc
.. doxygenfunction:: Asynch_Get_Links_Proc

Database
~~~~~~~~

.. doxygenfunction:: Asynch_Set_Database_Connection

Forcing
~~~~~~~

.. doxygenfunction:: Asynch_Activate_Forcing
.. doxygenfunction:: Asynch_Deactivate_Forcing

.. doxygenfunction:: Asynch_Get_First_Forcing_Timestamp
.. doxygenfunction:: Asynch_Set_First_Forcing_Timestamp

.. doxygenfunction:: Asynch_Get_Last_Forcing_Timestamp
.. doxygenfunction:: Asynch_Set_Last_Forcing_Timestamp

.. doxygenfunction:: Asynch_Set_Forcing_DB_Starttime

State Variables Getters and Setters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: Asynch_Set_Init_File
.. doxygenfunction:: Asynch_Set_System_State

State Variables Getters and Setters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



Custom Outputs
~~~~~~~~~~~~~~

Time series outputs can be customized so as to produce results specific to a particular model. For
instance, fluxes that are used internally to the model, or values that are not needed at all for
computing state solutions, may be outputted.

To create a custom output, specify the name of your new output in the global file (see Section [sec:
time series output]). Next, in your program after reading the global file but before performing
simulations, call the routine *Asynch_Set_Output* to set your output. Then perform any calculations
or further modifications as usual.

A custom routine must be provided by the user for determining how the output is calculated. The
specifications for this function are provided below.

Similarly, custom peakflow outputs can be created. These are set with a call to the routine
*Asynch_Set_Peakflow_Output*.

The routines for setting custom outputs are described below.

.. doxygenfunction:: Asynch_Set_Output_Int
.. doxygenfunction:: Asynch_Set_Output_Double
