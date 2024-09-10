#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <mpi.h>

#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif

#if defined(HAVE_POSTGRESQL)
#include <libpq-fe.h>
#endif

#if defined(HAVE_PETSC)
#include <petsc.h>
#endif

#if !defined(_MSC_VER)
#define ASYNCH_SLEEP(seconds) sleep(seconds)
#else
#include <windows.h>
#define ASYNCH_SLEEP(seconds) Sleep((seconds) * 1000)
#endif

#include <asynch_interface.h>

//Internal stuffs
#include <blas.h>
#include <comm.h>
#include <system.h>

#include <assim/models.h>
#include <assim/structs.h>
#include <assim/linear_least_squares.h>


/// This computes the least squares fit assuming the background and analysis difference is linear in the innovations. 
/// HM is (num_obs*steps_to_use) X allstates_needed
/// HM_els is 1 X allstates (i.e. a 1D array)
/// HM_buffer is 1 X allstates_needed
/// \param asynch The asynch solver instance
/// \param asynch The assimilation workspace
int LSSolveSys(AsynchSolver* asynch, AssimWorkspace* ws, double* q)
{
    //Unpack ptr
    unsigned int N = asynch->N;
    Link* sys = asynch->sys;
    GlobalVars* globals = asynch->globals;
    int* assignments = asynch->assignments;
    short int* getting = asynch->getting;
    unsigned int *obs_locs = ws->obs_locs, assim_dim = ws->assim_dim;
    unsigned int problem_dim = ws->problem_dim, allstates = ws->allstates;
    int *HM_col_indices = ws->HM_col_indices, *d_indices = ws->d_indices;
    double t_b = ws->t_b;
    unsigned int allstates_needed = ws->allstates_needed;
    double /**RHS_els,*/*x_start = ws->x_start, *HM_buffer = ws->HM_buffer;
    AsynchModel* custom_model = asynch->model;
    unsigned int *vareq_shift = ws->vareq_shift, *inv_vareq_shift = ws->inv_vareq_shift;

    double start = MPI_Wtime();

    unsigned int num_total_obs = ws->num_steps * ws->num_obs;

    //Build a vector structure for d
    //TODO Optimize this variable out
    Vec d;
    VecCreateSeq(MPI_COMM_SELF, num_total_obs, &d);	//!!!! This needs to be fixed !!!!
    VecSet(d, 0.0);
    VecSetValues(d, num_total_obs, d_indices, ws->d_full, INSERT_VALUES);
    VecAssemblyBegin(d);
    VecAssemblyEnd(d);

    //Initialize the system
    LSResetSys(sys, N, globals, t_b, x_start, assim_dim, globals->num_forcings, asynch->my_data);

    for (unsigned int i = 0; i < N; i++)
        if (assignments[i] == asynch->my_rank || getting[i])
            custom_model->initialize_eqs(
                globals->global_params, globals->num_global_params,
                sys[i].params, globals->num_params,
                sys[i].my->list.head->y_approx, sys[i].dim,
                sys[i].user);

    for (unsigned int i = 0; i < asynch->globals->num_forcings; i++)
    {
        //TODO Recurring and binary files may need this too
        if (asynch->forcings[i].flag == 3)
        {
            //printf("Setting to %u and %u\n",asynch->forcings[i]->first_file,asynch->forcings[i]->last_file);
            Asynch_Set_Forcing_State(asynch, i, t_b, asynch->forcings[i].first_file, asynch->forcings[i].last_file);
        }
    }

    //Advance the system and extract the HM matrix
    //HM here holds the values of M that are needed
    //!!!! Start at i=1? For i = 0, I don't think we need to set anything... !!!!    
    for (unsigned int i = 0; i < ws->num_steps; i++)
    {
        globals->t = 0.0;

        if (i > 0)
        {
            // Adjust the end of the simulation
            globals->maxtime = t_b + i * ws->obs_time_step;

            MPI_Barrier(MPI_COMM_WORLD);
            double start = MPI_Wtime();

            //Advance the simuation to globals->maxtime
            Asynch_Advance(asynch, 0);

            MPI_Barrier(MPI_COMM_WORLD);
            double stop = MPI_Wtime();

            if (asynch->my_rank == 0)
                printf("Time for advance to time %f: %.0f\n", globals->maxtime, stop - start);
        }

        //Build HM
        for (unsigned int j = 0; j < ws->num_obs; j++)
        {
            //Assumes only discharges
            //TODO Generalize this
            Link *current = &sys[obs_locs[j]];
            int owner = assignments[obs_locs[j]];
            bool is_my_link = (owner == asynch->my_rank);

            UpstreamData *updata = (UpstreamData*)(current->user);
            memset(HM_buffer, 0, allstates_needed * sizeof(double));

            //From my link
            if (is_my_link)
            {
                //Pull out needed data
                for (unsigned int n = 0; n < updata->num_fit_states; n++)
                {
                    if (asynch->verbose)
                        printf("ID = %u | Loading %e (from %u) into spot %u\n",
                            current->ID,
                            current->my->list.tail->y_approx[updata->fit_states[n]],
                            updata->fit_states[n],
                            vareq_shift[updata->fit_to_universal[n]]);

                    assert(updata->fit_states[n] < current->dim);
                    HM_buffer[vareq_shift[updata->fit_to_universal[n]]] = current->my->list.tail->y_approx[updata->fit_states[n]];
                }

                //Extract calculationed q's (Just needed for testing. Maybe...)
                q[i * ws->num_obs + j] = current->my->list.tail->y_approx[0];
            }

            //MPI_Bcast(HM_buffer, allstates_needed, MPI_DOUBLE, owner, MPI_COMM_WORLD);	//!!!! Only proc 0 needs this !!!!
            if (asynch->my_rank == 0)
            {
                MPI_Reduce(MPI_IN_PLACE, HM_buffer, allstates_needed, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

#if !defined(NDEBUG)
                unsigned int k;
                for (k = 0; k < allstates_needed; k++)
                    if (HM_buffer[k] != 0.)
                        break;

                assert(k < allstates_needed);
#endif

                unsigned int row_idx = i * ws->num_obs + j;
                MatSetValues(ws->HM, 1, &row_idx, allstates_needed, HM_col_indices, HM_buffer, INSERT_VALUES);
            }
            else
                MPI_Reduce(HM_buffer, NULL, allstates_needed, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


            MPI_Bcast(&(q[i * ws->num_obs + j]), 1, MPI_DOUBLE, owner, MPI_COMM_WORLD);

        }
    }

    double stop = MPI_Wtime();

    if (asynch->my_rank == 0)
        printf("Time for advance to time %f: %.0f\n", globals->maxtime, stop - start);

    if (asynch->my_rank == 0)
    {
        //Assemble the HM matrix
        MatAssemblyBegin(ws->HM, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(ws->HM, MAT_FINAL_ASSEMBLY);

        if (asynch->verbose)
        {
            printf("Matrix HM\n");
            MatView(ws->HM, PETSC_VIEWER_STDOUT_SELF);
        }

        start = MPI_Wtime();

        //Calculate innovations
        double *buffer = NULL;
        VecGetArray(d, &buffer);
        for (unsigned int i = 0; i < num_total_obs; i++)
            buffer[i] = buffer[i] - q[i];
        VecRestoreArray(d, &buffer);

        //Build the linear system \f$ A x = rhs \f$
        //HMTR is allstates_needed x (num_obs*max_or_steps)
        //HM is (num_obs*max_or_steps) x allstates_needed

        /// \f$ HMTR = H(y_0)^T R \f$
        /// HMTR is a temporary variable used for rhs and A computation
        MatTranspose(ws->HM, MAT_REUSE_MATRIX, &ws->HMTR);
        MatDiagonalScale(ws->HMTR, NULL, ws->R);

        /// \f$ A = B + H(y_0)^T R H(y_0) \f$
        MatMatMult(ws->HMTR, ws->HM, MAT_REUSE_MATRIX, PETSC_DEFAULT, &ws->HTH);
        MatDiagonalSet(ws->HTH, ws->B, ADD_VALUES);

        /// \f$ rhs = H(y_0)^T R \alpha(y_0^b) \f$
        MatMult(ws->HMTR, d, ws->rhs);
        MatAssemblyBegin(ws->HTH, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(ws->HTH, MAT_FINAL_ASSEMBLY);

        if (asynch->verbose)
        {
            printf("Matrix HTH\n");
            MatView(ws->HTH, PETSC_VIEWER_STDOUT_SELF);
        }

        //MPI_Barrier(MPI_COMM_WORLD);
        double stop = MPI_Wtime();
        if (asynch->my_rank == 0)
            printf("Time for matrix computations: %.0f\n", stop - start);

        //Compute analysis
        //MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();

        /// \f$ x = y_0 - y_0^b \f$
        KSPSetOperators(ws->ksp, ws->HTH, ws->HTH);     //Maybe not actually necessary
        KSPSolve(ws->ksp, ws->rhs, ws->x);
        KSPConvergedReason reason;
        KSPGetConvergedReason(ws->ksp, &reason);

        if (asynch->my_rank == 0)
            printf("Converged reason: %s\n", KSPConvergedReasons[reason]);

        //MPI_Barrier(MPI_COMM_WORLD);
        stop = MPI_Wtime();
        if (asynch->my_rank == 0)
            printf("Time for inversion: %.0f\n", stop - start);

        if (asynch->verbose)
        {
            printf("Solution x\n");
            VecView(ws->x, PETSC_VIEWER_STDOUT_SELF);
        }

        //Copy new solution to x_start
        VecGetArray(ws->x, &buffer);
        for (unsigned int i = 0; i < allstates_needed; i++)	//!!!! I think this is right... !!!!
        {
            x_start[inv_vareq_shift[i]] += buffer[i];
            //x_start[above_gauges[i]*assim_dim] += x_els[i];	//!!!! To skip hillslope !!!!
            //for(j=0;j<assim_dim;j++)
            //	x_start[above_gauges[i]*assim_dim+j] += x_els[i*assim_dim+j];
        }
        VecRestoreArray(ws->x, &buffer);
    }

    //Send solution to everyone
    MPI_Bcast(x_start, allstates, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (asynch->verbose && asynch->my_rank == 0)
    {
        //unsigned int idxm[num_obs*steps_to_use];
        //double temp_matptr[(num_obs*steps_to_use*allstates_needed > allstates_needed*allstates_needed) ? num_obs*steps_to_use*allstates_needed : allstates_needed*allstates_needed];
        //for(i=0;i<num_obs*steps_to_use;i++)
        //    idxm[i] = i;

        //printf("x_start\n");
        //for(i=0;i<allstates;i++)
        //    printf("%.15e ",x_start[i]);
        //printf("\n");

        double* buffer;
        printf("difference (x)\n");
        VecGetArray(ws->x, &buffer);
        for (unsigned int i = 0; i < allstates_needed; i++)
            printf("%.2e ", buffer[i]);
        printf("\n");
        VecRestoreArray(ws->x, &buffer);

        printf("d\n");
        VecGetArray(d, &buffer);
        for (unsigned int i = 0; i < num_total_obs; i++)
            printf("%.2e ", buffer[i]);
        printf("\n");
        VecRestoreArray(d, &buffer);

        //printf("HM\n");
        //MatGetValues(*HM,num_obs*max_or_steps,idxm,allstates_needed,cols_allstates_needed,temp_matptr);
        //for(i=0;i<num_obs*max_or_steps;i++)
        //{
           // for(j=0;j<allstates_needed;j++)
              //  printf("%.15e ",temp_matptr[i*allstates_needed + j]);
           // printf(";\n");
        //}

        //printf("HTH\n");
        //MatGetValues(*HTH,allstates_needed,cols_allstates_needed,allstates_needed,cols_allstates_needed,temp_matptr);
        //for(i=0;i<allstates_needed;i++)
        //{
           // for(j=0;j<allstates_needed;j++)
              //  printf("%.15e ",temp_matptr[i*allstates_needed + j]);
           // printf(";\n");
        //}
    }

    //Clean up
    VecDestroy(&d);

    MPI_Barrier(MPI_COMM_WORLD);
    stop = MPI_Wtime();
    if (asynch->my_rank == 0)
        printf("Total time for linear least squares fit: %.0f\n", stop - start);

    return 0;
}

// COmput least square distance 
double LSComputeDistance(const double * const d, const double * const q, unsigned int size)
{
    unsigned int i;
    double result = 0.0;

    for (i = 0; i < size; i++)
        result += (d[i] - q[i]) * (d[i] - q[i]);

    //return pow(result, 0.5);
    return result;
}


void LSResetSys(Link* sys, unsigned int N, GlobalVars* globals, double t_0, double* x_start, unsigned int problem_dim, unsigned int num_forcings, TransData* my_data)
{
    unsigned i, j, k, l;
    Link* current;

    Flush_TransData(my_data);

    for (i = 0; i < N; i++)
    {
        current = &sys[i];
        if (current->my != NULL)
        {
            while (current->current_iterations > 1)
            {
                Remove_Head_Node(&current->my->list);
                (current->current_iterations)--;
            }
            current->my->list.head->t = t_0;
            current->last_t = t_0;
            current->steps_on_diff_proc = 1;
            current->iters_removed = 0;
            current->rejected = 0;
            if (current->num_parents == 0)
                current->ready = 1;
            else
                current->ready = 0;
            for (j = 0; j < problem_dim; j++)
                current->my->list.head->y_approx[j] = x_start[i*problem_dim + j];
            //v_copy(backup[i],current->my->list.head->y_approx);

            //Reset the next_save time
            if (current->save_flag)
            {
                current->next_save = t_0;		//!!!! This forces the print times to match up with the assimilation times !!!!
                                                //current->disk_iterations = 1;
            }

            //Reset peak flow information
            current->peak_time = t_0;
            dcopy(current->my->list.head->y_approx, current->peak_value, 0, current->dim);

            //Set hydrograph scale
            //current->Q_TM = backup[i]->ve[0];

            //Reset current state
            if (current->check_state != NULL)
                current->state = current->check_state(
                    current->my->list.head->y_approx, current->dim,
                    globals->global_params, globals->num_global_params,
                    current->params, globals->num_params,
                    current->qvs, current->has_dam, NULL);
            current->my->list.head->state = current->state;

            //Set forcings
            if (current->my->forcing_data)
            {
                for (k = 0; k < num_forcings; k++)
                {
                    if (current->my->forcing_values[k])
                    {
                        //Find the right index in forcings
                        for (l = 0; l < current->my->forcing_data[k].num_points - 1; l++)
                            if (current->my->forcing_data[k].data[l].time <= t_0 && t_0 < current->my->forcing_data[k].data[l + 1].time)	break;
                        double rainfall_buffer = current->my->forcing_data[k].data[l].value;
                        current->my->forcing_values[k] = rainfall_buffer;
                        current->my->forcing_indices[k] = l;

                        //Find and set the new change in data
                        for (j = l + 1; j < current->my->forcing_data[k].num_points; j++)
                        {
                            if (fabs(current->my->forcing_data[k].data[j].value - rainfall_buffer) > 1e-12)
                            {
                                current->my->forcing_change_times[k] = current->my->forcing_data[k].data[j].time;
                                break;
                            }
                        }
                        if (j == current->my->forcing_data[k].num_points)
                            current->my->forcing_change_times[k] = current->my->forcing_data[k].data[j - 1].time;

                        //Select new step size
                        //current->h = InitialStepSize(current->last_t,current,globals,workspace);
                    }
                }
            }
        }
    }
}

void Print_MATRIX(double** A, unsigned int m, unsigned int n)
{
    unsigned int i, j;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)	printf("%.2e ", A[i][j]);
        printf(";\n");
    }
    printf("\n");
}

void Print_VECTOR(double* v, unsigned int dim)
{
    unsigned int i;

    for (i = 0; i < dim; i++)
        printf("[%d]: %.2e\n", i, v[i]);
    printf(";\n");
}

