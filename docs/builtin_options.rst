Built-in Options
================

Through a global file, many options are selected, including which model to use, what time series to output, and how to calculate peakflows. Each of these can be customized by the user. However, several precreated options do exist. This section provides a list of these options.

Built-In Output Time Series
---------------------------

Figure [fig: built-in output time series] contains the names and a description of built-in output time series. These outputs are defined in the source file *modeloutputs.c*. Up to seven states can be outputted with the built-in output time series. In addition to these, users can create their own time series outputs. See Section :ref:`Custom Outputs`.

.. table:: Built-in output time series

  +---------------+--------------------------------------------+
  | Output Name   | Description                                |
  +===============+============================================+
  | Time          | Simulation time                            |
  +---------------+--------------------------------------------+
  | TimeI         | Simulation time, truncated to an integer   |
  +---------------+--------------------------------------------+
  | State0        | State 0 of the model                       |
  +---------------+--------------------------------------------+
  | State1        | State 1 of the model                       |
  +---------------+--------------------------------------------+
  | ...           |                                            |
  +---------------+--------------------------------------------+
  | State6        | State 6 of the model                       |
  +---------------+--------------------------------------------+

Built-In Peakflow Functions
---------------------------

Two built-in peakflow functions exist: *Classic* and *Forecast*. The two are described in Figure [fig: built-in peakflow functions]. The peak discharges are the largest values obtained in the state with index 0 in the state vectors. The time to peak for the *Classic* function is given in simulation time. For *Forecast*, the time to peak is measured in unix time. The time period output is a parameter that can be altered by user programs to provide additional output information.

.. table:: Built-in peakflow functions

  +-----------------+--------------------------------------------------------+
  | Function Name   | Outputs                                                |
  +=================+========================================================+
  | Classic         | Link ID, upstream area, time to peak, peak discharge   |
  +-----------------+--------------------------------------------------------+
  | Forecast        | Link ID, time to peak, peak discharge, time period     |
  +-----------------+--------------------------------------------------------+

Built-In Runge-Kutta Methods
----------------------------

The ASYNCH solver is based upon using Runge-Kutta methods at the link level. These methods are selected either in the input global file or in a Runge-Kutta data file by the *RK index* in Table :ref:`rk-methods`.

.. _rk-methods:

.. table:: Built-in RK methods

  +------------+-------------------------------+-----------------------------+
  | RK index   | Name                          | Local order / Dense order   |
  +============+===============================+=============================+
  | 0          | Kutta’s Method                | 3 / 2                       |
  +------------+-------------------------------+-----------------------------+
  | 1          | The RK Method                 | 4 / 3                       |
  +------------+-------------------------------+-----------------------------+
  | 2          | Dormand and Prince’s Method   | 5 / 4                       |
  +------------+-------------------------------+-----------------------------+
  | 3          | RadauII 3A                    | 3 / 2                       |
  +------------+-------------------------------+-----------------------------+

The application of these methods is done through the *RKSolver* routine in the *UnivVars* structure. This is set with a call to the *InitRoutines* method. See Section [sec: initroutines]. Several choices exist for the *RKSolver*. They are given in Table  :ref:`rk-solvers`. Some solvers are only appropriate if the model uses ODEs, while others support DAEs. Similarly, some methods support discontinuity states, while others do not. Currently, only one method is equipped to handle stiff ODEs. Certainly, the routine *ExplicitRKIndex1SolverDam* could be used to solve any problem. However, using a more appropriate solver is significantly more efficient.

.. _rk-solvers:

.. table:: Built-in RK solvers

  +-----------------------------+--------+-------------------+---------+
  | Name                        | DAEs   | Discontinuities   | Stiff   |
  +=============================+========+===================+=========+
  | ExplicitRKSolver            | No     | No                | No      |
  +-----------------------------+--------+-------------------+---------+
  | ExplicitRKIndex1SolverDam   | Yes    | Yes               | No      |
  +-----------------------------+--------+-------------------+---------+
  | ExplicitRKIndex1Solver      | Yes    | No                | No      |
  +-----------------------------+--------+-------------------+---------+
  | ExplicitRKSolverDiscont     | No     | Yes               | No      |
  +-----------------------------+--------+-------------------+---------+
  | RadauRKSolver               | No     | No                | Yes     |
  +-----------------------------+--------+-------------------+---------+
