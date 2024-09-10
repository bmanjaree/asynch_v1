

#if defined(ASYNCH_HAVE_ASSIM_SOLVER)

//Computes one step of a method to solve the ODE at a link. Assumes parents have enough computed solutions.
//This is for use with data assimilation. It assumes the vectors are a lot smaller than they really are.
//Link* link_i: the link to apply a numerical method to.
//Returns 1 if the step was successfully taken, 0 if the step was rejected.
int ExplicitRKSolver_DataAssim(Link* link_i, UnivVars* GlobalVars, int* assignments, bool print_flag, FILE* outputfile, ConnData* conninfo, Forcing* forcings, Workspace* workspace)
{
    unsigned int i, j, l, m, idx;
    VEC* new_y;
    RKSolutionNode *curr_node[link_i->num_parents], *node, *new_node;
    Link* currentp;
    double t_needed, timediff, current_theta;

    //Some variables to make things easier to read
    VEC* y_0 = link_i->my->list.tail->y_approx;
    double h = link_i->h;
    double t = link_i->my->list.tail->t;
    VEC2* A = link_i->method->A;
    VEC* b = link_i->method->b;
    VEC* c = link_i->method->c;
    unsigned int s = link_i->method->s;
    VEC* params = link_i->params;
    VEC* e = link_i->method->e;
    VEC* d = link_i->method->d;
    RKMethod* meth = link_i->method;
    ErrorData* error = link_i->error_data;
    //const unsigned int dim = GlobalVars->dim;
    //unsigned int num_dense = GlobalVars->num_dense;
    unsigned int* dense_indices = GlobalVars->dense_indices;
    unsigned int num_print = GlobalVars->num_print;
    VEC* temp = workspace->temp;
    VEC* sum = workspace->sum;
    VEC*** temp_parent_approx = workspace->temp_parent_approx;
    VEC** temp_k = workspace->temp_k;
    unsigned int problem_dim = GlobalVars->problem_dim;

    //Calculate total number of upstream links
    unsigned int practical_dim, num_dense_parents[link_i->num_parents], num_upstream_links = 0;
    for (i = 0; i<link_i->num_parents; i++)
    {
        num_upstream_links += link_i->numupstream[i];
        num_dense_parents[i] = (link_i->numupstream[i] * problem_dim + 1) * problem_dim;	//!!!! This is actually an upper bound. The extra problem_dim assumes all states are needed. !!!!
    }
    practical_dim = 2 * problem_dim + (problem_dim - 1)*(problem_dim - 1) + problem_dim * num_upstream_links;

    //Get the approximate solutions from each parent
    for (i = 0; i<link_i->num_parents; i++)
    {
        currentp = link_i->parents[i];
        curr_node[i] = currentp->my->list.head;

        for (j = 0; j<s; j++)
        {
            //Find the needed value of t and corresponding node for y_p
            //Assuming everything needed is already calculated
            //This assumes c_s is the biggest. If not, one extra step may not get freed.
            t_needed = min(t + v_at(c, j)*h, currentp->last_t);
            //t_needed = t + v_at(c, j)*h;

            //Find the corresponding theta value and approximate solution
            while (t_needed > curr_node[i]->t)
                curr_node[i] = curr_node[i]->next;
            if (curr_node[i] != currentp->my->list.head)	curr_node[i] = curr_node[i]->prev;

            timediff = curr_node[i]->next->t - curr_node[i]->t;
            current_theta = (t_needed - curr_node[i]->t) / timediff;
            currentp->method->dense_b(current_theta, currentp->method->b_theta);

            //for(m=0;m<num_dense;m++)
            for (m = 0; m<num_dense_parents[i]; m++)
            {
                idx = dense_indices[m];
                temp_parent_approx[j][i].ve[idx] = curr_node[i]->y_approx.ve[idx];
                for (l = 0; l<currentp->method->s; l++)
                    temp_parent_approx[j][i].ve[idx] += timediff*currentp->method->b_theta.ve[l] * curr_node[i]->next->k[l].ve[m];
            }
            link_i->check_consistency(temp_parent_approx[j][i], params, GlobalVars->global_params);
            //printf("ID = %u parent = %u num_dense_parents = %u\n",link_i->ID,currentp->ID,num_dense_parents[i]);
            //Print_Vector(temp_parent_approx[j][i]);
        }
        //getchar();
    }

    //Do the RK method to get the next approximation

    //Setup variables for the new data
    new_node = New_Step(&link_i->my->list);
    new_node->t = t + h;
    new_y = new_node->y_approx;

    //Compute the k's
    for (i = 0; i<s; i++)
    {
        v_copy_n(y_0, sum, practical_dim);
        for (j = 0; j<i; j++)
            daxpy_u(h*A.me[i][j], v2_slice(temp_k, j), sum, 0, practical_dim);
        link_i->check_consistency(sum, params, GlobalVars->global_params);
        link_i->f(t + v_at(c, i) * h, sum, temp_parent_approx[i], link_i->num_parents, GlobalVars->global_params, link_i->forcing_values, link_i->qvs, params, link_i->state, link_i->user, v2_slice(temp_k, i));
    }

    //Build the solution
    v_copy_n(y_0, new_y, practical_dim);
    for (i = 0; i<s; i++)	daxpy_u(h*v_at(b, i), v2_slice(temp_k, i), new_y, 0, practical_dim);
    link_i->check_consistency(new_y, params, GlobalVars->global_params);

    //Error estimation and step size selection

    //Check the error of y_1 (in inf norm) to determine if the step can be accepted
    double err_1;
    v_copy_n(v2_slice(temp_k, 0), sum, practical_dim);
    sv_mlt_u(h * v_at(e, 0), sum, 0, practical_dim);
    for (i = 1; i<s; i++)	daxpy_u(h*v_at(e, i), v2_slice(temp_k, i), sum, 0, practical_dim);

    //Build SC_i
    //for(i=0;i<dim;i++)
    for (i = 0; i<practical_dim; i++)
        temv_at(p, i) = max(fabs(new_v_at(y, i)), fabs(y_v_at(0, i))) * error->reltov_at(l, i) + error->abstov_at(l, i);

    err_1 = norm_inf_u(sum, temp, 0, practical_dim);
    double value_1 = pow(1.0 / err_1, 1.0 / meth->e_order);

    //Check the dense error (in inf norm) to determine if the step can be accepted
    double err_d;
    v_copy_n(v2_slice(temp_k, 0), sum, practical_dim);
    sv_mlt_u(h * v_at(d, 0), sum, 0, practical_dim);
    for (i = 1; i<s; i++)	daxpy_u(h*v_at(d, i), v2_slice(temp_k, i), sum, 0, practical_dim);

    //for(i=0;i<dim;i++)
    for (i = 0; i<practical_dim; i++)
        temv_at(p, i) = max(fabs(new_v_at(y, i)), fabs(y_v_at(0, i))) * error->reltol_densv_at(e, i) + error->abstol_densv_at(e, i);

    err_d = norm_inf_u(sum, temp, 0, practical_dim);
    double value_d = pow(1.0 / err_d, 1.0 / meth->d_order);

    //Determine a new step size for the next step
    double step_1 = h*min(error->facmax, max(error->facmin, error->fac * value_1));
    double step_d = h*min(error->facmax, max(error->facmin, error->fac * value_d));
    link_i->h = min(step_1, step_d);

    //printf("ID = %u t_0 = %e 1 = %e d = %e\n",link_i->ID,t,step_1,step_d);
    //Print_Vector(new_y);
    //getchar();

    if (err_1 < 1.0 && err_d < 1.0)
    {
        //Check if a discontinuity has been stepped on
        if (link_i->discont_count > 0 && (t + h) >= link_i->discont[link_i->discont_start])
        {
            (link_i->discont_count)--;
            link_i->discont_start = (link_i->discont_start + 1) % GlobalVars->discont_size;
            link_i->h = InitialStepSize(link_i->last_t, link_i, GlobalVars, workspace);
        }

        //Save the new data
        link_i->last_t = t + h;
        link_i->current_iterations++;
        store_k(temp_k, new_node->k, s, dense_indices, num_dense);

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
                node = link_i->my->list.tail->prev;
                current_theta = (link_i->next_save - t) / h;
                link_i->method->dense_b(current_theta, link_i->method->b_theta);
                //for(m=0;m<num_dense;m++)
                for (m = 0; m<problem_dim; m++)	//!!!! This is wasteful if not all states are printed. Also forbids printing variational eq states. !!!!
                {
                    idx = dense_indices[m];
                    sum.ve[idx] = node->y_approx.ve[idx];
                    for (l = 0; l<link_i->method->s; l++)
                        sum.ve[idx] += h*link_i->method->b_theta.ve[l] * node->next->k[l].ve[m];
                }
                link_i->check_consistency(sum, params, GlobalVars->global_params);

                //Write to a file
                WriteStep(link_i->next_save, sum, GlobalVars, params, link_i->state, outputfile, link_i->output_user, &(link_i->pos_offset));

                link_i->next_save += link_i->print_time;
            }
        }

        //Check if this is a max discharge
        if (link_i->peak_flag && (new_v_at(y, 0) > link_i->peak_valuv_at(e, 0)))
        {
            v_copy_n(new_y, link_i->peak_value, practical_dim);
            link_i->peak_time = link_i->last_t;
        }

        //Check if the newest step is on a change in rainfall
        short int propagated = 0;	//Set to 1 when last_t has been propagated
        for (j = 0; j<GlobalVars->num_forcings; j++)
        {
            if (forcings[j].active && link_i->forcing_buff[j] && (fabs(link_i->last_t - link_i->forcing_change_times[j]) < 1e-8))
            {
                //Propagate the discontinuity to downstream links
                if (!propagated)
                {
                    propagated = 1;
                    Link* next = link_i->child;
                    Link* prev = link_i;
                    for (i = 0; i<GlobalVars->max_localorder && next != NULL; i++)
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
                for (l = link_i->forcing_indices[j] + 1; l<link_i->forcing_buff[j]->n_times; l++)
                    if (fabs(link_i->forcing_change_times[j] - link_i->forcing_buff[j]->rainfall[l][0]) < 1e-8)	break;
                link_i->forcing_indices[j] = l;

                double forcing_buffer = link_i->forcing_buff[j]->rainfall[l][1];
                link_i->forcing_values[j] = forcing_buffer;

                //Find and set the new change in rainfall
                for (i = l + 1; i<link_i->forcing_buff[j]->n_times; i++)
                {
                    if (fabs(link_i->forcing_buff[j]->rainfall[i][1] - forcing_buffer) > 1e-8)
                    {
                        link_i->forcing_change_times[j] = link_i->forcing_buff[j]->rainfall[i][0];
                        break;
                    }
                }
                if (i == link_i->forcing_buff[j]->n_times)
                    link_i->forcing_change_times[j] = link_i->forcing_buff[j]->rainfall[i - 1][0];
            }
        }

        //Select new step size, if forcings changed
        if (propagated)	link_i->h = InitialStepSize(link_i->last_t, link_i, GlobalVars, workspace);

        //Free up parents' old data
        for (i = 0; i<link_i->num_parents; i++)
        {
            currentp = link_i->parents[i];
            while (currentp->my->list.head != curr_node[i])
            {
                Remove_Head_Node(currentp->list);
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

#endif // defined(ASYNCH_HAVE_ASSIM_SOLVER)
