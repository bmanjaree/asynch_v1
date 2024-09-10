#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <math.h>

#include <minmax.h>
#include <system.h>
#include <blas.h>
#include <io.h>
#include <rksteppers.h>

//Computes one step of a method to solve the ODE at a link. Assumes parents have enough computed solutions.
//This method should be used for links that have a is_dam.
//Link* link_i: the link to apply a numerical method to.
//VEC* sum: some space for temporary calculations. Should have same dimension as number of equations.
//VEC* temp: same as sum.
//Returns 1 if the step was successfully taken, 0 if the step was rejected.
int ExplicitRKSolverDiscont(Link* link_i, GlobalVars* GlobalVars, int* assignments, bool print_flag, FILE* outputfile, ConnData* conninfo, Forcing* forcings, Workspace* workspace)
{
    unsigned int i, j, l, m, idx;
    VEC new_y;
    RKSolutionNode *curr_node[ASYNCH_LINK_MAX_PARENTS], *node, *new_node;
    Link* currentp;
    double t_needed, timediff, current_theta;

    //Some variables to make things easier to read
    VEC y_0 = link_i->my->list.tail->y_approx;
    double h = link_i->h;
    double t = link_i->my->list.tail->t;
    VEC2 A = link_i->method->A;
    VEC b = link_i->method->b;
    VEC c = link_i->method->c;
    unsigned int s = link_i->method->s;
    VEC params = link_i->params;
    VEC e = link_i->method->e;
    VEC d = link_i->method->d;
    RKMethod* meth = link_i->method;
    ErrorData* error = &link_i->my->error_data;
    const unsigned int dim = link_i->dim;
    unsigned int num_dense = link_i->num_dense;
    unsigned int* dense_indices = link_i->dense_indices;
    unsigned int num_outputs = GlobalVars->num_outputs;
    unsigned int* print_indices = GlobalVars->print_indices;
    VEC temp = workspace->temp;
    VEC sum = workspace->sum;
    VEC2 temp_k = workspace->temp_k;

    //Get the approximate solutions from each parent
    for (i = 0; i < link_i->num_parents; i++)
    {
        currentp = link_i->parents[i];
        curr_node[i] = currentp->my->list.head;
        for (j = 0; j < s; j++)
        {
            //Find the needed value of t and corresponding node for y_p
            //Assuming everything needed is already calculated
            //This assumes c_s is the biggest. If not, one extra step may not get freed.
            t_needed = min(t + v_at(c, j) * h, currentp->last_t);
            //t_needed = t + v_at(c, j)*h;

            //Find the corresponding theta value and approximate solution
            while (t_needed > curr_node[i]->t)
                curr_node[i] = curr_node[i]->next;
            if (curr_node[i] != currentp->my->list.head)	curr_node[i] = curr_node[i]->prev;

            timediff = curr_node[i]->next->t - curr_node[i]->t;
            current_theta = (t_needed - curr_node[i]->t) / timediff;
            currentp->method->dense_b(current_theta, currentp->method->b_theta);

            VEC parent_approx = v3_slice2(workspace->temp_parent_approx, j, i);

            //v_copy(curr_node[i]->y_approx,temp_parent_approx[j][i]);
            //for(l=0;l<currentp->method->s;l++)
            //	daxpy(timediff*currentp->method->b_theta.ve[l],curr_node[i]->next->k[l],temp_parent_approx[j][i],1);
            for (m = 0; m < currentp->num_dense; m++)
            {
                idx = currentp->dense_indices[m];
                double approx = v_at(curr_node[i]->y_approx, idx);
                for (l = 0; l < currentp->method->s; l++)
                    approx += timediff * v_at(currentp->method->b_theta, l) * v2_at(curr_node[i]->next->k, l, m);

                v_set(parent_approx, idx, approx);
            }

            //Build the algebraic variable
            currentp->check_consistency(parent_approx, currentp->params, GlobalVars->global_params);
        }
    }

    //Do the RK method to get the next approximation

    //Setup variables for the new data
    new_node = New_Step(&link_i->my->list);
    new_node->t = t + h;
    //k = new_node->k;
    new_y = new_node->y_approx;

    //Compute the k's
    for (i = 0; i < s; i++)
    {
        v_copy_n(y_0, sum, link_i->dim);
        for (j = 0; j < i; j++)
            daxpy(h * v2_at(A, i, j), v2_slice(temp_k, j), sum, 0, link_i->dim);
        link_i->check_consistency(sum, params, GlobalVars->global_params);
        link_i->differential(t + v_at(c, i) * h, sum, v3_slice(workspace->temp_parent_approx, i), GlobalVars->global_params, link_i->forcing_values, link_i->qvs, params, link_i->state, link_i->user, v2_slice(temp_k, i));
    }

    //Build the solution
    v_copy_n(y_0, new_y, link_i->dim);
    for (i = 0; i < s; i++)	daxpy_u(h*v_at(b, i), v2_slice(temp_k, i), new_y, 0, link_i->dim);
    link_i->check_consistency(new_y, params, GlobalVars->global_params);
    new_node->state = link_i->check_state(new_y, GlobalVars->global_params, params, link_i->qvs, link_i->is_dam);

    //Error estimation and step size selection

    //Check the error of y_1 (in inf norm) to determine if the step can be accepted
    double err_1;
    v_copy_n(v2_slice(temp_k, 0), sum, link_i->dim);
    dscal(h * v_at(e, 0), sum, 0, link_i->dim);
    for (i = 1; i < s; i++)	daxpy_u(h*v_at(e, i), v2_slice(temp_k, i), sum, 0, link_i->dim);

    //Build SC_i
    for (i = 0; i < dim; i++)
        v_set(temp, i, max(fabs(v_at(new_y, i)), fabs(v_at(y_0, i)) * error->reltol[i] + error->abstol[i]));

    //err_1 = norm_inf(sum,temp,meth->e_order_ratio,1);
    err_1 = nrminf2(sum, temp, 0, link_i->dim);
    double value_1 = pow(1.0 / err_1, 1.0 / meth->e_order);


    //Check the dense error (in inf norm) to determine if the step can be accepted
    double err_d;
    v_copy_n(v2_slice(temp_k, 0), sum, link_i->dim);
    dscal(h * v_at(d, 0), sum, 0, link_i->dim);
    for (i = 1; i < s; i++)	daxpy_u(h*v_at(d, i), v2_slice(temp_k, i), sum, 0, link_i->dim);
    for (i = 0; i < dim; i++)
        v_set(temp, i, max(fabs(v_at(new_y, i)), fabs(v_at(y_0, i)) * error->reltol_dense[i] + error->abstol_dense[i]));

    //err_d = norm_inf(sum,temp,meth->d_order_ratio,1);
    err_d = nrminf2(sum, temp, 0, link_i->dim);
    double value_d = pow(1.0 / err_d, 1.0 / meth->d_order);

    //Determine a new step size for the next step
    double step_1 = h*min(error->facmax, max(error->facmin, error->fac * value_1));
    double step_d = h*min(error->facmax, max(error->facmin, error->fac * value_d));
    link_i->h = min(step_1, step_d);
    //link_i->h = step_1;

    if (err_1 < 1.0 && err_d < 1.0)
        //if(err_1 < 1.0)
    {
        //Check for issues with state discontinuities
        if (link_i->rejected == 2)
        {
            int old_state = link_i->state;
            link_i->state = link_i->check_state(new_y, GlobalVars->global_params, params, link_i->qvs, link_i->is_dam);

            if (old_state != link_i->state)
            {
                //Pass the time to the next link
                double xh = t + h;
                Link* next = link_i->child;
                Link* prev = link_i;

                for (i = 1; i < GlobalVars->max_localorder && next != NULL; i++)
                {
                    if (assignments[next->location] == my_rank && i < next->method->localorder)
                    {
                        //Insert the time into the discontinuity list
                        next->discont_end = Insert_Discontinuity(xh, next->discont_start, next->discont_end, &(next->discont_count), GlobalVars->discont_size, next->discont, next->ID);
                    }
                    else if (assignments[next->location] != my_rank)
                    {
                        //Store the time to send to another process
                        Insert_SendDiscontinuity(xh, i, &(prev->discont_send_count), GlobalVars->discont_size, prev->discont_send, prev->discont_order_send, prev->ID);
                        break;
                    }

                    prev = next;
                    next = next->child;
                }

                link_i->h = InitialStepSize(link_i->last_t, link_i, GlobalVars, workspace);
            }
        }
        else
        {
            //Check for discontinuities
            unsigned int sign_l, sign_h, sign_m;
            //sign_l = link_i->check_state(y_0,GlobalVars->global_params,params,link_i->qvs,link_i->is_dam);
            sign_l = link_i->state;
            sign_h = link_i->check_state(new_y, GlobalVars->global_params, params, link_i->qvs, link_i->is_dam);
            new_node->state = sign_h;
            if (sign_l != sign_h)
            {
                double xl = t;
                double xh = t + h;
                double xm = (xl + xh) / 2.0;

                //Form the derivative of the interpolant
                meth->dense_bderiv(1.0, meth->b_theta_deriv);
                v_zero(sum);
                for (i = 0; i < s; i++)
                    daxpy_u(h * v_at(meth->b_theta_deriv, i), v2_slice(temp_k, i), sum, 0, link_i->dim);

                //Get the approximate solutions from each parent
                for (i = 0; i < link_i->num_parents; i++)
                {
                    currentp = link_i->parents[i];
                    curr_node[i] = currentp->my->list.head;

                    //Find the needed value of t and corresponding node for y_p
                    //Assuming everything needed is already calculated
                    t_needed = min(t + h, currentp->last_t);
                    //t_needed = t + h;

                    //Find the corresponding theta value and approximate solution
                    while (t_needed > curr_node[i]->t)
                        curr_node[i] = curr_node[i]->next;
                    if (curr_node[i] != currentp->my->list.head)	curr_node[i] = curr_node[i]->prev;

                    timediff = curr_node[i]->next->t - curr_node[i]->t;
                    current_theta = (t_needed - curr_node[i]->t) / timediff;
                    currentp->method->dense_b(current_theta, currentp->method->b_theta);

                    VEC parent_approx = v3_slice2(workspace->temp_parent_approx, 0, i);

                    for (m = 0; m < currentp->num_dense; m++)
                    {
                        idx = currentp->dense_indices[m];
                        double approx = v_at(curr_node[i]->y_approx, idx);
                        for (l = 0; l < currentp->method->s; l++)
                            approx += timediff * v_at(currentp->method->b_theta, l) * v2_at(curr_node[i]->next->k, l, m);

                        v_set(parent_approx, idx, approx);
                    }
                }

                //Exact derivative at time t + h
                link_i->differential(t + h, new_y, v3_slice(workspace->temp_parent_approx, 0), GlobalVars->global_params, link_i->forcing_values, link_i->qvs, params, link_i->state, link_i->user, temp);

                //Form the defect for a stopping criteria
                dsub(sum, temp, temp, 0);
                double norm_defect = vector_norminf(temp, 0);
                double prev_error;
                double curr_error = (xh - xl)*norm_defect;

                //Locate the discontinuity
                do
                {
                    prev_error = curr_error;

                    //Build the solution at time xm
                    v_copy_n(y_0, sum, link_i->dim);
                    current_theta = (xm - t) / h;
                    meth->dense_b(current_theta, meth->b_theta);
                    for (i = 0; i < s; i++)
                        daxpy_u(h * v_at(meth->b_theta, i), v2_slice(temp_k, i), sum, 0, link_i->dim);
                    sign_m = link_i->check_state(sum, GlobalVars->global_params, params, link_i->qvs, link_i->is_dam);

                    if (sign_l != sign_m)
                    {
                        sign_h = sign_m;
                        xh = xm;
                    }
                    else if (sign_m != sign_h)
                    {
                        sign_l = sign_m;
                        xl = xm;
                    }
                    else
                        printf("Uh oh....\n");

                    xm = (xl + xh) / 2.0;
                    curr_error = (xh - xl)*norm_defect;

                } while (curr_error > GlobalVars->discont_tol  &&  fabs(curr_error - prev_error) > 1e-13);

                //Prepare to redo the step
                if (curr_error <= GlobalVars->discont_tol)
                {
                    //Prepare to step on the time of the discontinuity next iteration
                    link_i->h = xh - t;
                    Undo_Step(&link_i->my->list);

                    return 2;
                }
                else
                {
                    //printf("[%i]: Notice: Discontinuity error has converged. Rejecting current step.\n",my_rank);
                    //printf("[%i]: ID = %u t = %f old h = %e current error = %.16e  previous error = %.16e\n",my_rank,link_i->ID,t,h,curr_error,prev_error);

                    //Cut the step size and redo the iteration
                    link_i->h = h * 0.5;
                    Undo_Step(&link_i->my->list);

                    return 0;
                }
            }
        }

        //Check if a propagated discontinuity has been stepped on
        if (link_i->discont_count > 0 && (t + h) >= link_i->discont[link_i->discont_start])
        {
            (link_i->discont_count)--;
            link_i->discont_start = (link_i->discont_start + 1) % GlobalVars->discont_size;
            link_i->h = InitialStepSize(link_i->last_t, link_i, GlobalVars, workspace);
        }

        //Save the new data
        link_i->last_t = t + h;
        link_i->current_iterations++;
        store_k(temp_k, link_i->dim, new_node->k, s, dense_indices, num_dense);

        //Check if new data should be written to disk
        if (print_flag)
        {
            //while(t <= link_i->next_save && link_i->next_save <= link_i->last_t)
            while (t <= link_i->next_save && (link_i->next_save < link_i->last_t || fabs(link_i->next_save - link_i->last_t) / link_i->next_save < 1e-12))
            {
                /*
                //Don't write anything if using data assimilation and at a time when data is available
                if(GlobalVars->assim_flag)
                {
                double rounded = 1e-13*rint(1e13*(GlobalVars->maxtime - link_i->next_save));
                if(rounded < 1e-13 && -rounded < 1e-13)		break;
                }
                */
                if (link_i->disk_iterations == link_i->expected_file_vals)
                {
                    printf("[%i]: Warning: Too many steps computed for link id %u. Expected no more than %u. No more values will be stored for this link.\n", my_rank, link_i->ID, link_i->expected_file_vals);
                    break;
                }
                (link_i->disk_iterations)++;
                node = link_i->my->list.tail->prev;
                current_theta = (link_i->next_save - t) / h;
                link_i->method->dense_b(current_theta, link_i->method->b_theta);

                //v_copy(node->y_approx,sum);
                //for(l=0;l<link_i->method->s;l++)
                //	daxpy(h * link_i->method->b_theta.ve[l],node->next->k[l],sum,0);
                for (m = 0; m < num_dense; m++)
                {
                    idx = dense_indices[m];
                    double approx = v_at(node->y_approx, idx);
                    for (l = 0; l < link_i->method->s; l++)
                        approx += h * v_at(link_i->method->b_theta, l) * v2_at(node->next->k, l, m);

                    v_set(sum, idx, approx);
                }
                link_i->check_consistency(sum, params, GlobalVars->global_params);

                //Write to a file
                //fsetpos(outputfile,&(link_i->pos));
                //WriteStep(link_i->next_save,sum,GlobalVars,params,link_i->state,outputfile,link_i->output_user,&(link_i->pos));
                WriteStep(outputfile, link_i->ID, link_i->next_save, sum, GlobalVars, params, link_i->state, link_i->output_user, &(link_i->pos_offset));
                //fgetpos(outputfile,&(link_i->pos));
                /*
                fsetpos(outputfile,&(link_i->pos));
                fwrite(&(link_i->next_save),sizeof(double),1,outputfile);
                for(j=0;j<num_print;j++)
                fwrite(&(sum.ve[print_indices[j]]),sizeof(double),1,outputfile);
                fgetpos(outputfile,&(link_i->pos));
                */
                /*
                #ifdef PRINT2DATABASE
                sprintf(conninfo->query,"%i,%u,%.12e,%.12e\n",link_i->ID,(unsigned int)(link_i->next_save * 60.0) + conninfo->time_offset,suv_at(m, 0)/link_i->Q_TM,suv_at(m, 0));
                unsigned int length = strlen(conninfo->query);
                (conninfo->submission_content) += length;
                while(conninfo->submission_content > conninfo->submission_size)
                {
                //Allocate more space
                (conninfo->submission_size) *= 2;
                conninfo->submission = realloc(conninfo->submission,conninfo->submission_size);
                }
                strcat(conninfo->submission,conninfo->query);
                #else //Write to a file
                fsetpos(outputfile,&(link_i->pos));
                fwrite(&(link_i->next_save),sizeof(double),1,outputfile);
                for(j=0;j<num_print;j++)
                fwrite(&(sum.ve[print_indices[j]]),sizeof(double),1,outputfile);
                fgetpos(outputfile,&(link_i->pos));
                #endif
                */
                link_i->next_save += link_i->print_time;
            }
        }

        //Check if this is a max discharge
        if (link_i->peak_flag && (v_at(new_y, 0) > v_at(link_i->peak_value, 0)))
        {
            v_copy_n(new_y, link_i->peak_value, link_i->dim);
            link_i->peak_time = link_i->last_t;
        }

        //Check if the newest step is on a change in rainfall
        short int propagated = 0;	//Set to 1 when last_t has been propagated
        for (j = 0; j < GlobalVars->num_forcings; j++)
        {
            if (forcings[j].active && link_i->forcing_data[j] && fabs(link_i->last_t - link_i->forcing_change_times[j]) < 1e-8)
            {
                //Propagate the discontinuity to downstream links
                if (!propagated)
                {
                    propagated = 1;
                    Link* next = link_i->child;
                    Link* prev = link_i;
                    for (i = 0; i < GlobalVars->max_localorder && next != NULL; i++)
                    {
                        if (assignments[next->location] == my_rank && i < next->method->localorder)
                        {
                            //Insert the time into the discontinuity list
                            next->discont_end = Insert_Discontinuity(link_i->forcing_change_times[j], next->discont_start, next->discont_end, &(next->discont_count), GlobalVars->discont_size, next->discont, next->ID);
                        }
                        else if (next != NULL && assignments[next->location] != my_rank)
                        {
                            //Store the time to send to another process
                            Insert_SendDiscontinuity(link_i->forcing_change_times[j], i, &(prev->discont_send_count), GlobalVars->discont_size, prev->discont_send, prev->discont_order_send, prev->ID);
                            break;
                        }

                        prev = next;
                        next = next->child;
                    }
                }

                //Find the right index in rainfall
                //for(l=1;l<link_i->forcing_data[j]->n_times;l++)
                for (l = link_i->forcing_indices[j] + 1; l < link_i->forcing_data[j]->nrows; l++)
                    if (fabs(link_i->forcing_change_times[j] - link_i->forcing_data[j]->data[l][0]) < 1e-8)	break;
                link_i->forcing_indices[j] = l;

                double forcing_buffer = link_i->forcing_data[j]->data[l][1];
                v_set(link_i->forcing_values, j, forcing_buffer);

                //Find and set the new change in rainfall
                for (i = l + 1; i < link_i->forcing_data[j]->nrows; i++)
                {
                    //if(link_i->forcing_data[j]->rainfall[i][1] != forcing_buffer)
                    if (fabs(link_i->forcing_data[j]->data[i][1] - forcing_buffer) < 1e-8)
                    {
                        link_i->forcing_change_times[j] = link_i->forcing_data[j]->data[i][0];
                        break;
                    }
                }
                if (i == link_i->forcing_data[j]->nrows)
                    link_i->forcing_change_times[j] = link_i->forcing_data[j]->data[i - 1][0];
            }
        }

        //Select new step size, if forcings changed
        if (propagated)	link_i->h = InitialStepSize(link_i->last_t, link_i, GlobalVars, workspace);

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

        return 1;
    }
    else
    {
        //Trash the data from the failed step
        Undo_Step(&link_i->my->list);

        return 0;
    }
}