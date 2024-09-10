#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <memory.h>
#include <math.h>

#include <minmax.h>
#include <system.h>
#include <blas.h>
#include <io.h>
#include <rksteppers.h>

//Computes a solution for a location with forced system states.
//Computes a solution at either the last time of the upstream link, or the time when a change in the system state occurs.
//Should return 1.
int ForcedSolutionSolver(Link* link_i, GlobalVars* globals, int* assignments, bool print_flag, FILE* outputfile, ConnData* conninfo, Forcing* forcings, Workspace* workspace)
{
    unsigned int i, j, l;
    double *new_y;
    RKSolutionNode *curr_node[ASYNCH_LINK_MAX_PARENTS], *new_node;
    Link* currentp;
    double t_needed;
    short int change_value = 0;

    //Some variables to make things easier to read
    double *y_0 = link_i->my->list.tail->y_approx;
    double t = link_i->my->list.tail->t;
    double h = link_i->h;
    unsigned int num_stages = link_i->method->num_stages;
    double *params = link_i->params;
    RKMethod* meth = link_i->method;
    ErrorData* error = link_i->my->error_data;
    const unsigned int dim = link_i->dim;
    unsigned int num_dense = link_i->num_dense;
    unsigned int* dense_indices = link_i->dense_indices;
    unsigned int num_outputs = globals->num_outputs;
    double **temp_k = workspace->temp_k_slices;

    //Find the next time to step on
    t_needed = globals->maxtime;
    for (i = 0; i < globals->num_forcings; i++)
    {
        if (forcings[i].active && (link_i->my->forcing_data[i].num_points > 0))
        {
            if (t_needed > link_i->my->forcing_change_times[i])
            {
                t_needed = link_i->my->forcing_change_times[i];
                change_value = 1;
            }
        }
    }

    for (i = 0; i < link_i->num_parents; i++)
    {
        if (t_needed > link_i->parents[i]->last_t)
        {
            t_needed = link_i->parents[i]->last_t;
            change_value = 0;
        }
    }

    h = t_needed - t;

    //Setup the current nodes at each parent (for deleting data)
    for (i = 0; i < link_i->num_parents; i++)
    {
        currentp = link_i->parents[i];
        curr_node[i] = currentp->my->list.head;

        //Find the corresponding theta value and approximate solution
        while (t_needed > curr_node[i]->t)
            curr_node[i] = curr_node[i]->next;

        if (curr_node[i] != currentp->my->list.head)	curr_node[i] = curr_node[i]->prev;
    }

    //Setup variables for the new data
    new_node = New_Step(&link_i->my->list);
    new_node->t = t + h;
    new_y = new_node->y_approx;

    //Compute the k's
    memset(workspace->temp_k, 0, num_stages * dim * sizeof(double));

    //Build the solution
    if (change_value)
    {
        //Check if the newest step is on a change in rainfall
        short int propagated = 0;	//Set to 1 when last_t has been propagated
        for (j = 0; j < globals->num_forcings; j++)
        {
            if (forcings[j].active && (link_i->my->forcing_data[j].num_points > 0) && (fabs(new_node->t - link_i->my->forcing_change_times[j]) < 1e-8))
            {
                //Propagate the discontinuity to downstream links
                if (!propagated)
                {
                    propagated = 1;
                    Link* next = link_i->child;
                    Link* prev = link_i;
                    for (i = 0; i < globals->max_localorder && next != NULL; i++)
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
                //for(l=1;l<link_i->my->forcing_data[j]->n_times;l++)
                for (l = link_i->my->forcing_indices[j] + 1; l < link_i->my->forcing_data[j].num_points; l++)
                    if (fabs(link_i->my->forcing_change_times[j] - link_i->my->forcing_data[j].data[l].time) < 1e-8)
                        break;
                link_i->my->forcing_indices[j] = l;

                double forcing_buffer = link_i->my->forcing_data[j].data[l].value;
                link_i->my->forcing_values[j] = forcing_buffer;

                //Find and set the new change in rainfall
                for (i = l + 1; i < link_i->my->forcing_data[j].num_points; i++)
                {
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
    }
    link_i->differential(
        t + h,
        y_0, dim,
        NULL, 0, 0,
        globals->global_params, params, link_i->my->forcing_values, link_i->qvs, link_i->state, link_i->user, new_y);
    if (link_i->check_state)
        new_node->state = link_i->check_state(new_y, dim, globals->global_params, globals->num_global_params, link_i->params, link_i->num_params, link_i->qvs, link_i->has_dam, link_i->user);

    //Set stepsize
    link_i->h = h;

    //Ignore propagated discontinuities
    link_i->discont_count = 0;
    link_i->discont_start = 0;
    link_i->discont_end = globals->discont_size - 1;

    //Save the new data
    link_i->last_t = t + h;
    link_i->current_iterations++;
    store_k(workspace->temp_k, globals->max_dim, new_node->k, num_stages, dense_indices, num_dense);

    //Check if new data should be written to disk
    if (print_flag)
    {
        while (t <= link_i->next_save && (link_i->next_save < link_i->last_t || fabs(link_i->next_save - link_i->last_t) / link_i->next_save < 1e-12))
        {
            if (link_i->disk_iterations == link_i->expected_file_vals)
            {
                printf("[%i]: Warning: Too many steps computed for link id %u. Expected no more than %u. No more values will be stored for this link.\n", my_rank, link_i->ID, link_i->expected_file_vals);
                break;
            }
            (link_i->disk_iterations)++;

            //Write to a file
            if (change_value && fabs((link_i->next_save - link_i->last_t) / link_i->next_save) < 1e-12)
                WriteStep(globals->outputs, globals->num_outputs, outputfile, link_i->ID, link_i->next_save, new_y, dim, &link_i->pos_offset);
            else
                WriteStep(globals->outputs, globals->num_outputs, outputfile, link_i->ID, link_i->next_save, y_0, dim, &link_i->pos_offset);
            link_i->next_save += link_i->print_time;
        }
    }

    //Check if this is a max discharge
    if (link_i->peak_flag && (new_y[0] > link_i->peak_value[0]))
    {
        dcopy(new_y, link_i->peak_value, 0, link_i->dim);
        link_i->peak_time = link_i->last_t;
    }


    //Check if the newest step is on a change in rainfall
    if (!change_value)
    {
        short int propagated = 0;	//Set to 1 when last_t has been propagated
        for (j = 0; j < globals->num_forcings; j++)
        {
            if (forcings[j].active && (link_i->my->forcing_data[j].num_points > 0) && (fabs(link_i->last_t - link_i->my->forcing_change_times[j]) < 1e-8))
            {
                //Propagate the discontinuity to downstream links
                if (!propagated)
                {
                    propagated = 1;
                    Link* next = link_i->child;
                    Link* prev = link_i;
                    for (i = 0; i < globals->max_localorder && next != NULL; i++)
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
                //for(l=1;l<link_i->my->forcing_data[j]->n_times;l++)
                for (l = link_i->my->forcing_indices[j] + 1; l < link_i->my->forcing_data[j].num_points; l++)
                    if (fabs(link_i->my->forcing_change_times[j] - link_i->my->forcing_data[j].data[l].time) < 1e-8)	break;
                link_i->my->forcing_indices[j] = l;

                double forcing_buffer = link_i->my->forcing_data[j].data[l].value;
                link_i->my->forcing_values[j] = forcing_buffer;

                //Find and set the new change in rainfall
                for (i = l + 1; i < link_i->my->forcing_data[j].num_points; i++)
                {
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
    }

    //Free up parents' old data
    for (i = 0; i < link_i->num_parents; i++)
    {
        currentp = link_i->parents[i];
        while (currentp->my->list.head != curr_node[i])
        {
            Remove_Head_Node(&currentp->my->list);
            currentp->current_iterations--;
            currentp->iters_removed++;
        }
    }

    //if(link_i->ID == 2)
    //printf("%f %f %u\n",link_i->last_t,link_i->my->forcing_values[2],link_i->my->forcing_indices[2]);

    return 1;
}
