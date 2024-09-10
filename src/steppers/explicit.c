#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <math.h>
#include <memory.h>

#include <minmax.h>
#include <system.h>
#include <blas.h>
#include <io.h>
#include <rksteppers.h>


//Computes one step of a method to solve the ODE at a link. Assumes parents have enough computed solutions.
//Link* link_i: the link to apply a numerical method to.
//Returns 1 if the step was successfully taken, 0 if the step was rejected.
int ExplicitRKSolver(Link* link_i, GlobalVars* globals, int* assignments, bool print_flag, FILE* outputfile, ConnData* conninfo, Forcing* forcings, Workspace* workspace)
{
    unsigned int idx;
    
    RKSolutionNode *curr_node[ASYNCH_LINK_MAX_PARENTS], *node, *new_node;
    double t_needed, current_theta;

    //Some variables to make things easier to read
    double *y_0 = link_i->my->list.tail->y_approx;
    double h = link_i->h;
    double t = link_i->my->list.tail->t;
    const double * const A = link_i->method->A;
    double *b = link_i->method->b;
    const double * const c = link_i->method->c;
    const double * const e = link_i->method->e;
    const double * const d = link_i->method->d;
    unsigned int num_stages = link_i->method->num_stages;
    RKMethod* meth = link_i->method;
    ErrorData* error = link_i->my->error_data;
    unsigned int dim = link_i->dim;
    unsigned int num_dense = link_i->num_dense;
    unsigned int* dense_indices = link_i->dense_indices;
    double *temp = workspace->temp;
    double *sum = workspace->sum;
    double **temp_k = workspace->temp_k_slices;  

	// adlz
	// printf("+Using ExplicitRKSolver...\n");

    //Get the approximate solutions from each parent
    for (unsigned int i = 0; i < link_i->num_parents; i++)
    {
        Link* curr_parent = link_i->parents[i];
        curr_node[i] = curr_parent->my->list.head;
        for (unsigned int j = 0; j < num_stages; j++)
        {
            //Find the needed value of t and corresponding node for y_p
            //Assuming everything needed is already calculated
            //This assumes c_s is the biggest. If not, one extra step may not get freed.
            t_needed = min(t + c[j] * h, curr_parent->last_t);

            //Find the corresponding theta value and approximate solution
            while (t_needed > curr_node[i]->t)
                curr_node[i] = curr_node[i]->next;
            if (curr_node[i] != curr_parent->my->list.head)
                curr_node[i] = curr_node[i]->prev;

            double dt = curr_node[i]->next->t - curr_node[i]->t;
            current_theta = (t_needed - curr_node[i]->t) / dt;
            curr_parent->method->dense_b(current_theta, curr_parent->method->b_theta);

            //[num_stages][max_parents][max_dim] -> [max_dim]
            double *parent_approx = workspace->stages_parents_approx
                + j * globals->max_parents * globals->max_dim
                + i * globals->max_dim;

            for (unsigned int m = 0; m < curr_parent->num_dense; m++)
            {
                idx = curr_parent->dense_indices[m];
                double approx = curr_node[i]->y_approx[idx];

                for (int l = 0; l < curr_parent->method->num_stages; l++)
                    approx += dt * curr_parent->method->b_theta[l] *
                        curr_node[i]->next->k[l * curr_parent->num_dense + m];

                parent_approx[idx] = approx;
            }

            link_i->check_consistency(parent_approx, curr_parent->dim, globals->global_params, globals->num_global_params, curr_parent->params, link_i->num_params, curr_parent->user);
        }
    }

    //Do the RK method to get the next approximation

    //Setup variables for the new data
    new_node = New_Step(&link_i->my->list);
    new_node->t = t + h;
    //k = new_node->k;
    double *new_y = new_node->y_approx;

    //Compute the k's
    for (unsigned int i = 0; i < num_stages; i++)
    {
        //v_copy_n(y_0, sum, link_i->dim);
		/* printf("+begin: ");
		for (unsigned int i = 0; i < link_i->dim; i++)
			printf("%f, ", y_0[i]);
		printf(" (y_0)\n"); */
        memcpy(sum, y_0, link_i->dim * sizeof(double));
        for (unsigned int j = 0; j < i; j++)
        {
            double alpha = h * A[i * num_stages + j];
            daxpy(alpha, temp_k[j], sum, 0, link_i->dim);
        }

        link_i->check_consistency(sum, link_i->dim, globals->global_params, globals->num_global_params, link_i->params, link_i->num_params, link_i->user);

        //[num_stages][max_parents][max_dim]
        double *y_p = workspace->stages_parents_approx + i * globals->max_parents * globals->max_dim;

        double dt = c[i] * h;

		/* printf("+Pre diff: ");
		for (unsigned int i = 0; i < link_i->dim; i++)
			printf("%f, ", sum[i]);
		printf(" (sum)\n"); */

        link_i->differential(
            t + dt,
            sum, link_i->dim,
            y_p, link_i->num_parents, globals->max_dim,
            globals->global_params,
            link_i->params,
            link_i->my->forcing_values,
            link_i->qvs,
            link_i->state,
            link_i->user,
            temp_k[i]);

		/* printf("+Pos diff: ");
		for (unsigned int i = 0; i < link_i->dim; i++)
			printf("%f, ", temp_k[i]);
		printf(" (temp_k)\n"); */
    }

    //Build the solution
    dcopy(y_0, new_y, 0, link_i->dim);
    for (unsigned int i = 0; i < num_stages; i++)
        daxpy(h * b[i], temp_k[i], new_y, 0, link_i->dim);

    // Check constistency
    link_i->check_consistency(new_y, link_i->dim, globals->global_params, globals->num_global_params, link_i->params, link_i->num_params, link_i->user);

    new_node->state = link_i->state;
    

    //Error estimation and step size selection

    //Check the error of y_1 (in inf norm) to determine if the step can be accepted
    double err_1;
    dcopy(temp_k[0], sum, 0, link_i->dim);
    dscal(h * e[0], sum, 0, link_i->dim);
    for (unsigned int i = 1; i < num_stages; i++)
        daxpy(h * e[i], temp_k[i], sum, 0, link_i->dim);

    //Build SC_i
    for (unsigned int i = 0; i < dim; i++)
        temp[i] = max(fabs(new_y[i]), fabs(y_0[i])) * error->reltol[i] + error->abstol[i];

    err_1 = nrminf2(sum, temp, 0, link_i->dim);
	/* printf("+err_1: from ");
	for (unsigned int i = 0; i < dim; i++)
		printf("[%f, %f], ", sum[i], temp[i]);
	printf("got %f.\n", err_1); */

    double value_1 = pow(1.0 / err_1, 1.0 / meth->e_order);

    //Check the dense error (in inf norm) to determine if the step can be accepted
    double err_d;
    dcopy(temp_k[0], sum, 0, link_i->dim);
    dscal(h * d[0], sum, 0, link_i->dim);

    for (unsigned int i = 1; i < num_stages; i++)
        daxpy(h * d[i], temp_k[i], sum, 0, link_i->dim);

    for (unsigned int i = 0; i < dim; i++)
        temp[i] = max(fabs(new_y[i]), fabs(y_0[i])) * error->reltol_dense[i] + error->abstol_dense[i];

    err_d = nrminf2(sum, temp, 0, link_i->dim);
	/* adlz
	printf("+err_d: from ");
	for (unsigned int i = 0; i < dim; i++)
		printf("[%f, %f], ", sum[i], temp[i]);
	printf("got %f.\n", err_d); */

    double value_d = pow(1.0 / err_d, 1.0 / meth->d_order);

    //Determine a new step size for the next step
    double step_1 = h * min(error->facmax, max(error->facmin, error->fac * value_1));
    double step_d = h * min(error->facmax, max(error->facmin, error->fac * value_d));
    link_i->h = min(step_1, step_d);

    if (err_1 < 1.0 && err_d < 1.0)
    {
        //Check if a discontinuity has been stepped on
        if (link_i->discont_count > 0 && (t + h) >= link_i->discont[link_i->discont_start])
        {
            (link_i->discont_count)--;
            link_i->discont_start = (link_i->discont_start + 1) % globals->discont_size;
            link_i->h = InitialStepSize(link_i->last_t, link_i, globals, workspace);
        }

        //Save the new data
		// adlz
		// printf("+Saving the data as (%f < 1.0) and (%f < 1.0)...\n", err_1, err_d);
        link_i->last_t = t + h;
        link_i->current_iterations++;
        store_k(workspace->temp_k, globals->max_dim, new_node->k, num_stages, dense_indices, num_dense);

        //Check if new data should be written to disk
        if (print_flag)
        {
            //while(t <= link_i->next_save && link_i->next_save <= link_i->last_t)
            while (t <= link_i->next_save && (link_i->next_save < link_i->last_t || fabs(link_i->next_save - link_i->last_t) / link_i->next_save < 1e-12))
            {
                if (link_i->disk_iterations == link_i->expected_file_vals)
                {
                    printf("[%i]: Warning: Too many steps computed for link id %u. Expected no more than %u. No more values will be stored for this link.\n", my_rank, link_i->ID, link_i->expected_file_vals);
                    break;
                }
                (link_i->disk_iterations)++;
                node = link_i->my->list.tail->prev;
                current_theta = (link_i->next_save - t) / h;
                link_i->method->dense_b(current_theta, link_i->method->b_theta);
                for (unsigned int m = 0; m < num_dense; m++)
                {
                    idx = dense_indices[m];
                    double approx = node->y_approx[idx];
                    for (unsigned int l = 0; l < link_i->method->num_stages; l++)
                        approx += h * link_i->method->b_theta[l] * node->next->k[l * num_dense + m];

                    sum[idx] = approx;
                }

                link_i->check_consistency(sum, link_i->dim, globals->global_params, globals->num_global_params, link_i->params, link_i->num_params, link_i->user);

                //Write to a file
                WriteStep(globals->outputs, globals->num_outputs, outputfile, link_i->ID, link_i->next_save, sum, link_i->dim, &(link_i->pos_offset));

                link_i->next_save += link_i->print_time;
            }
        }

        //Check if this is a peak value
        if (link_i->peak_flag && (new_y[0] > link_i->peak_value[0]))
        {
            dcopy(new_y, link_i->peak_value, 0, link_i->dim);
            link_i->peak_time = link_i->last_t;
        }

        //Check if the newest step is on a change in rainfall
        short int propagated = 0;	//Set to 1 when last_t has been propagated
        for (unsigned int j = 0; j < globals->num_forcings; j++)
        {
            if (forcings[j].active && (link_i->my->forcing_data[j].num_points > 0) && (fabs(link_i->last_t - link_i->my->forcing_change_times[j]) < 1e-8))
            {
                //Propagate the discontinuity to downstream links
                if (!propagated)
                {
                    propagated = 1;
                    Link* next = link_i->child;
                    Link* prev = link_i;
                    for (unsigned int i = 0; i < globals->max_localorder && next != NULL; i++)
                    {
                        if (assignments[next->location] == my_rank && i < next->method->localorder)
                        {
                            //Insert the time into the discontinuity list
                            next->discont_end = Insert_Discontinuity(link_i->my->forcing_change_times[j], next->discont_start, next->discont_end, &(next->discont_count), globals->discont_size, next->discont, next->ID);
                        }
                        else if (next != NULL && assignments[next->location] != my_rank)
                        {
                            //Store the time to send to another process
                            Insert_SendDiscontinuity(link_i->my->forcing_change_times[j], i, &(prev->discont_send_count), globals->discont_size, prev->discont_send, prev->discont_order_send, prev->ID);
                            break;
                        }

                        prev = next;
                        next = next->child;
                    }
                }

                //Find the right index in rainfall
                //for(l=1;l<link_i->my->forcing_data[j].n_times;l++)
                unsigned int l;
                for (l = link_i->my->forcing_indices[j] + 1; l < link_i->my->forcing_data[j].num_points; l++)
                    if (fabs(link_i->my->forcing_change_times[j] - link_i->my->forcing_data[j].data[l].time) < 1e-8)
                        break;
                link_i->my->forcing_indices[j] = l;

                double forcing_buffer = link_i->my->forcing_data[j].data[l].value;
                link_i->my->forcing_values[j] = forcing_buffer;

                //Find and set the new change in rainfall
                unsigned int i;
                for (i = l + 1; i < link_i->my->forcing_data[j].num_points; i++)
                {
                    //if(link_i->my->forcing_data[j].rainfall[i].value != forcing_buffer)
                    if (fabs(link_i->my->forcing_data[j].data[i].value - forcing_buffer) > 1e-8)
                    {
                        link_i->my->forcing_change_times[j] = link_i->my->forcing_data[j].data[i].time;
                        break;
                    }
                }
                if (i == link_i->my->forcing_data[j].num_points)
                    link_i->my->forcing_change_times[j] = link_i->my->forcing_data[j].data[i - 1].time;
            }
        }

        //Select new step size, if forcings changed
        if (propagated)
            link_i->h = InitialStepSize(link_i->last_t, link_i, globals, workspace);

        //Free up parents' old data
        for (unsigned int i = 0; i < link_i->num_parents; i++)
        {
            Link *curr_parent = link_i->parents[i];
            while (curr_parent->my->list.head != curr_node[i])
            {
                Remove_Head_Node(&curr_parent->my->list);
                curr_parent->current_iterations--;
                curr_parent->iters_removed++;
            }
        }

        return 1;
    }
    else
    {
		// adlz
		// printf("+Trashing the data as (%f > 1.0) or (%f > 1.0)...\n", err_1, err_d);

        //Trash the data from the failed step
        Undo_Step(&link_i->my->list);

        return 0;
    }
}
