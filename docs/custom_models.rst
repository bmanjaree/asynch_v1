Custom Models
=============

An extendible collection of models is built-in ASYNCH. The evaluation of the differential and algebraic equations occurs in the source code *problems.c*, while the definition of the models (i.e. number of parameters, precalculations, etc.) is set in *definetype.c*. New models can be added here by modifying those two source files, plus adding needed function declarations to *problems.h*.

Every built in model is given a unique id known as the *model type*. This nonnegative integer value is used to identify the model throughout the initialization process. The model type is specified in the global file used to initialize ASYNCH. User defined models are possible, which can be created outside ASYNCH’s built-in collection of models.

Model Definition
----------------

The definition of every model is given in *definetype.c*. This module consists of five routines used to initialize each model. A description of the contents of each of these routines is given below.

SetParamSizes
~~~~~~~~~~~~~

This routine defines several integer values needed to describe a model.

+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Name              | Description                                                                                                                                                                                                                                                                                                                                                                                  |
+===================+==============================================================================================================================================================================================================================================================================================================================================================================================+
| dim               | The number of states modeled by the differential and algebraic equations. This parameter is also the number of such equations at a single link.                                                                                                                                                                                                                                              |
+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| template_flag     | A flag to specify if the model uses an XML parser for evaluating the differential equations. 0 indicates no parser, 1 indicates a parser is used.                                                                                                                                                                                                                                            |
+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| assim_flag        | A flag to specify if the model uses a data assimilation scheme. 0 indicates no data assimilation, 1 indicates data assimilation is used. This feature will be removed in future versions.                                                                                                                                                                                                    |
+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| diff_start        | The index in the equation-value vectors where the differential equations begin. All equations before this index are assumed to be algebraic. If all equations are differential, then this value should be 0.                                                                                                                                                                                 |
+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| no_ini_start      | The index in the state vectors corresponding to the first state not requiring initial conditions specified by the initial states specification in a global file. These states are generally initialized by other parameters of the model in the function *InitRoutines*. If all states require initialization through the global file (typical), then *no\_ini\_start* should be set to dim. |
+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| num_global_params | The number of global parameters specified in a global file for this model. If the number provided by the global file is less than expected, an error occurs. If more parameters are given than expected, a warning is given. The extra parameters are accessible. Although providing more parameters than needed is not recommended, it can be useful for testing.                           |
+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| uses_dam          | A flag to indicate if the model uses dams. 0 indicates no, 1 indicates yes. If dams are not available for this model, a value of 0 is expected for the *dam flag* in the global file.                                                                                                                                                                                                        |
+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| params_size       | The total number of local parameters available at each link. This includes all parameters read from the link parameters of the global file as well as all precalculations (specified in the *Precalculations*).                                                                                                                                                                              |
+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| iparams_size      | The number of integer valued parameters at each link. This may be removed in future versions.                                                                                                                                                                                                                                                                                                |
+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| dam_params_size   | The number of additional parameters at links with a dam. These parameters are included at the end of the vector of parameters at each link with a dam.                                                                                                                                                                                                                                       |
+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| area_idx          | The index in the parameter vector of the upstream area parameter. This para eter is used frequently with peakflow data.                                                                                                                                                                                                                                                                      |
+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| areah_idx         | The index in the parameter vector of the hillslope area. This parameter is frequently used with peakflow data.                                                                                                                                                                                                                                                                               |
+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| disk_params       | The number of local parameters available at each link read from a parameter file or database table. The *params\_size* minus the *disk\_params* is the number of recalculated parameters plus any dam parameters at the link.                                                                                                                                                                |
+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| num_dense         | The number of states passed down from one link to another. This number cannot be larger than dim. If equal to 0, then the links are totally disconnected.                                                                                                                                                                                                                                    |
+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| convertarea_flag  | Flag to indicate whether the model converts the parameter with index area\_idx from :math:`km^2` to :math:`m^2`. This can be needed for peakflow output data. The flag is set to 1 if the units are converted, 0 if not.                                                                                                                                                                     |
+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| num_forcings      | The number of forcings for the model. If a global file specifies less than this number of forcings, an error occurs. If more than this number of forcings is specified, a warning is given.                                                                                                                                                                                                  |
+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| dense_indices     | An array containing the indices in the state vectors that are passed from one link to another. This array must contain *num\_dense* indices.                                                                                                                                                                                                                                                 |
+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


ConvertParams
~~~~~~~~~~~~~

This routine allows for unit conversions on the local parameters at each link. These conversions occur immediately after the parameters are loaded into memory, and thus will be in place for all calculations. This feature is useful for when a data source provides values with units different than those expected by the model.

InitRoutines
~~~~~~~~~~~~

This routine specifies routines associated with the model. In this routine, the following arguments are available.

+---------+---------------------------------------------------------------------------------------------------------------------------+
| Name    | Description                                                                                                               |
+=========+===========================================================================================================================+
| link    | The current link where the routines are to be set.                                                                        |
+---------+---------------------------------------------------------------------------------------------------------------------------+
| type    | The model index.                                                                                                          |
+---------+---------------------------------------------------------------------------------------------------------------------------+
| exp_imp | A flag to determine if an implicit or explicit RK method is to be used. 0 if the method is explicit, 1 if it is implicit. |
+---------+---------------------------------------------------------------------------------------------------------------------------+
| dam     | A flag for whether a dam is present at this link. 0 if no dam is present, 1 if a dam is present.                          |
+---------+---------------------------------------------------------------------------------------------------------------------------+

The following routines must be set at each link.

+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Name             | Description                                                                                                                                                                                                                                  |
+==================+==============================================================================================================================================================================================================================================+
| RKSolver         | The routine for the numerical integrator.                                                                                                                                                                                                    |
+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| f                | The routine to evaluate the differential equations of the model.                                                                                                                                                                             |
+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| alg              | The routine to evaluate the algebraic equations of the model.                                                                                                                                                                                |
+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Jacobian         | The routine to evaluate the Jacobian of the system of differential equations. This must be set if an implicit RK method is used.                                                                                                             |
+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| state_check      | The routine to check in what discontinuity state the system is. The number of the discontinuity state is determined by the model.                                                                                                            |
+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| CheckConsistency | This routine alters the state vectors to be consistent with constraints of the system. **Notice: these constraints MUST exist in the exact solution of the equations for the link (for example, nonnegative solutions to a linear system).** |
+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Precalculations
~~~~~~~~~~~~~~~

This routine allows computations that are static in time and independent of state to be performed. The results are stored with the link parameters. This feature can be used to prevent redundant computations. The following information is accessible in this routine:

+---------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Name          | Description                                                                                                                                                           |
+===============+=======================================================================================================================================================================+
| link_i        | The current link where precalculations are performed.                                                                                                                 |
+---------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| global_params | The parameters which are constant in space and time.                                                                                                                  |
+---------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| params        | The parameters for this link. Results from this routine will be stored in this vector. Other parameters from a database or parameter file (.prm) are accessible here. |
+---------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| iparams       | The integer parameters for this link.                                                                                                                                 |
+---------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| disk_params   | The first entry of params that should be set for this location.                                                                                                       |
+---------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| params_size   | The first entry for dam parameters. These are only accessible if the dam flag is set.                                                                                 |
+---------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| dam           | The flag to indicate if a dam is present at the current link. If dam is 1, then a dam model is present here. If dam is 0, then a dam model is not present.            |
+---------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| type          | The index of the model.                                                                                                                                               |
+---------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Before exiting, all entries in params from index *disk_params* up to (but not including) *params_size* should be set.

ReadInitData
~~~~~~~~~~~~

This routine sets any initial conditions which are *not* determined through the *Initial Conditions* section of the global file (.gbl) (Section :ref:`Initial States`). Generally, this is to set the initial conditions for unknowns in models determined by algebraic equations, or those ODEs which have hardcoded initial conditions. The *ReadInitData* routine sets the initial conditions link by link. The following information is available in this routine:

+---------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Name          | Description                                                                                                                                                                                                                                                                                      |
+===============+==================================================================================================================================================================================================================================================================================================+
| global_params | The parameters which are constant in space and time.                                                                                                                                                                                                                                             |
+---------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| params        | The parameters for this link.                                                                                                                                                                                                                                                                    |
+---------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| iparams       | Integral parameters for this link.                                                                                                                                                                                                                                                               |
+---------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| qvs           | Discharge vs storage table. This information is available only if a dam is present at this link.                                                                                                                                                                                                 |
+---------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| dam           | The flag to indicate if a dam is present at the current link. If *dam* is 1, then a dam model is present here. If *dam* is 0, then a dam model is not present.                                                                                                                                   |
+---------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| y_0           | The vector of initial values. All indices between *diff\_start* (inclusive) and *no\_ini\_start* (exclusive) are set. These values were determined from the initial conditions specified in the global file. Both *diff\_start* and *no\_ini\_start* are defined by the routine *SetParamSizes*. |
+---------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| type          | The index of the model.                                                                                                                                                                                                                                                                          |
+---------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

The return value *ReadInitData* is the discontinuity state of the system, based upon the initial value vector *y\_0*.

Model Equations Definition
--------------------------

The equations for the model are defined in the file *problems.c*. Each set of built-in equations requires a routine to be defined here. Further, the differential and algebraic equations for a model must be defined in separate routines (although the routine for the differential equations may call the function for the algebraic equations). As is typical in C, any routines created in *problems.c* should be declared in *problems.h*. The routines defined here should be attached to each model in the *InitRoutines* method in *definetype.c*.

Differential Equations
~~~~~~~~~~~~~~~~~~~~~~

Every model must have a set of differential equations. The equations defined in this routine are for a single link only. Mathematically, the form of these equations should appear as

.. math::

  \frac{d y_s}{dt} &= f_s(...) \\
  \frac{d y_{s+1}}{dt} &= f_{s+1}(...) \\
  &\vdots \\
  \frac{d y_{dim}}{dt} &= f_{dim}(...)

where :math:`s` is *diff\_start*. Note that the index of the first state determined by a differential equation is *diff\_start* (or :math:`s` here). Thus, these states should appear after any states determined through algebraic equations in state and equation-value vectors. When the differential equation routine is called, the rate of change of each of the state variables :math:`y_i` is the expected output. Thus, this routine should evaluate all of the functions on the right of the equations. Examples of differential equations used for ASYNCH can be found in Section :ref:`Built-in Models`.

.. doxygentypedef:: DifferentialFunc

It is worth noting that only states from the upstream links are available in this routine. Dependence upon further upstream links breaks the underlying tree structure.

Algebraic Equations
~~~~~~~~~~~~~~~~~~~

Some models may have a set of algebraic equations. The equations defined in this routine are for a single link only. Mathematically, the form of these equations should appear as

.. math::

  y_0 &= g_0(...) \\
  y_1 &= g_1(...) \\
  &\vdots \\
  y_{s-1} &= g_{s-1}(...)

where :math:`s` is *diff\_start*. Note that the index of the first state determined by an algebraic equation is 0. Thus, these states should appear before any states determined through differential equations in state and equation-value vectors. When this routine is called, the expected output is the evaluation of the right side function. Support for algebraic equations is limited to explicit equations of the state variables. This means none of the states :math:`y_0`, ..., :math:`y_{s-1}` are available for use in this routine. Only the states defined through differential equations are available (:math:`y_s`, ..., :math:`y_{dim}`). Examples of models with algebraic equations can be found in Section :ref:`Built-in Models`.

.. doxygentypedef:: AlgebraicFunc

It is worth noting that only states from this link are available in this routine.

State Check
~~~~~~~~~~~

Some models may include discontinuities in the states of the system. This routine determines in which discontinuity state the system currently is. The return value is the integer representing the current discontinuity state.

.. doxygentypedef:: CheckStateFunc

System Consistency
~~~~~~~~~~~~~~~~~~

For many models, the equations describing the differential and algebraic system states come with built-in constraints. Common examples include non-negative values or maximum state values. These constraints may not necessarily be satisfied due to numerical errors. A routine for system consistency is called by the integrator to guarantee these constraints are satisfied.

.. doxygentypedef:: CheckConsistencyFunc

.. note::

  The solutions to the algebraic and differential equations MUST support these constraints. For instance, an equation with an exponential decaying solution has a minimum value for the solution. However, such an equation has no limit on the maximum value of its solution. Thus, a consistency routine can be created to impose the minimum value, but not a maximum value.

The values of states derived through algebraic equations are not available in the consistency routine. This is done for efficiency, as the algebraic states may not be needed to check consistency.
