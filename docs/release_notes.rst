Releases Notes
==============

ASYNCH release notes provide information on the features and improvements in each release. This page includes release notes for major releases and minor (bugfix) releases. If you are upgrading from an earlier version of ASYNCH, you will find essential information in the Breaking Changes associated with the relevant release notes.

Version 1.4
-----------

Breaking Changes
~~~~~~~~~~~~~~~~

The signature of the function that implements the differential equations has an extra parameter ``max_num_dim``. This is required because data assimilation uses a variable number of equations at links, the maximum number of degree of freedom beeing known at runtime.

.. code-block:: c

  typedef void (DifferentialFunc) (
    double t,
    const double * const y_i, unsigned int num_dof,
    const double * const y_p, unsigned short num_parents, unsigned int max_num_dof,
    const double * const global_params,
    const double * const params,
    const double * const forcing_values,
    const QVSData * const qvs,
    int state,
    void *user,
    double *ans);

New Features
~~~~~~~~~~~~

This release introduces :ref:`Data Assimilation`.

An additional, more compact HDF5 format for the time series is available, see option ``6`` of :ref:`Time Series Location`.

Three new models are available:

Constant Runoff Model ``195``:
  * based on ``190`` with precipitation forcing replaced by runoff and infiltration forcings to describe spatial variability of such processes, and including an additional state for rainfall accumulation.

Top Layer Model ``256``:
  * based on ``254``;
  * has a new state variable ``ans[7]`` : accumulated evapotranspiration;
  * has one more toplayer storage to link parameter ``k_tl = global_params[12]``;
  * has a flux component from toplayer storage to channel that represents the interflow :math:`q_{tl} = k_{tl} s_t`.

Top Layer Model ``257``:
  * based on ``256`` but channel velocity is spatialized as a function of the Horton order;
  * the Horton order is an additional local parameter.

Version 1.3
-----------

This release is the result of a loong run profiling and optimization work.

Memory footprint improvements, better data structure reduce the memory usage typically by a factor two. If you have memory available you may want to increase the number of buffers, see :ref:`Buffer Sizes`.

Performance improvements, ASYNCH runs about 30% faster.

Version 1.2
-----------

Breaking Changes
~~~~~~~~~~~~~~~~

The global file structure has changed. The time of the simulation is now given in absolute time, see :ref:`Simulation period`.

The snapshots in HDF5 format has changed, see :ref:`Ini HDF5 Files`.

New Features
~~~~~~~~~~~~

Time series outputs can be written in HDF5 format, see :ref:`Time Series Location`.
