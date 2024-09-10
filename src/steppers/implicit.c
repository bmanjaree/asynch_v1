#if defined(ASYNCH_HAVE_IMPLICIT_SOLVER)

//Computes one step of an implicit RK method to solve the ODE at a link. Assumes parents have enough computed solutions.
//Link* link_i: the link to apply a numerical method to.
//VEC* sum: some space for temporary calculations. Should have same dimension as number of equations.
//VEC* temp: same as sum.
//Returns 1 if the step was successfully taken, 0 if the step was rejected.
int RadauRKSolver(Link* link_i, UnivVars* GlobalVars, int* assignments, bool print_flag, FILE* outputfile, ConnData* conninfo, Forcing* forcings, Workspace* workspace)
{
    unsigned int i, j, l, m;
    VEC** k;
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
    const unsigned int dim = link_i->dim;
    VEC* temp = workspace->temp;
    VEC* sum = workspace->sum;
    VEC*** temp_parent_approx = workspace->temp_parent_approx;

    //Get the approximate solutions from each parent
    for (i = 0; i<link_i->num_parents; i++)
    {
        currentp = link_i->parents[i];
        curr_node[i] = currentp->my->list.head;

        //Find the needed value of t_0 and corresponding node for y_p for the extra stage
        t_needed = t;

        //Find the corresponding theta value and approximate solution
        while (t_needed > curr_node[i]->t)
            curr_node[i] = curr_node[i]->next;
        if (curr_node[i] != currentp->my->list.head)	curr_node[i] = curr_node[i]->prev;

        timediff = curr_node[i]->next->t - curr_node[i]->t;
        current_theta = (t_needed - curr_node[i]->t) / timediff;
        currentp->method->dense_b(current_theta, currentp->method->b_theta);

        v_copy(curr_node[i]->y_approx, temp_parent_approx[s][i]);
        for (l = 0; l<currentp->method->s; l++)
            daxpy(timediff*currentp->method->b_theta.ve[l], curr_node[i]->next->k[l], temp_parent_approx[s][i], 0);


        //currentp = link_i->parents[i];
        //curr_node[i] = currentp->my->list.head;
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

            v_copy(curr_node[i]->y_approx, temp_parent_approx[j][i]);
            for (l = 0; l<currentp->method->s; l++)
                daxpy(timediff*currentp->method->b_theta.ve[l], curr_node[i]->next->k[l], temp_parent_approx[j][i], 0);
        }
    }

    //Do the RK method to get the next approximation

    //Allocate space. !!!! Do this elsewhere !!!!
    int info;
    int* ipiv = workspace->ipiv;
    VEC* RHS = workspace->RHS;
    //	VEC2* CoefMat = workspace->CoefMat;
    VEC2* CoefMat = link_i->CoefMat;
    //	VEC2* JMatrix = workspace->JMatrix;
    VEC2* tempmat = workspace->JMatrix;
    VEC** Z_i = workspace->Z_i;
    VEC** old_Z = link_i->Z_i;
    VEC* err = workspace->err;
    unsigned int iters = 10;
    double normDelta = -1.0;
    double normDeltam1 = -1.0;
    double Theta, eta, error_est;
    unsigned short int abort = 0;
    unsigned short int newt = 0;
    unsigned short int newtmax = 20;

    //Construct a guess for the nonlinear solver
    if (link_i->h_old < 0.0)
        for (i = 0; i<s; i++)	for (j = 0; j<dim; j++)	Z_i[i].ve[j] = 0.0;
    else
    {
        for (i = 0; i<s; i++)
        {
            lagrange_bary(1.0 + h / link_i->h_old*v_at(c, i), c, old_Z, meth->w, Z_i[i]);
            v_add(Z_i[i], link_i->sol_diff, Z_i[i], 0);
            if (Z_i[i].ve[0] + y_v_at(0, 0) < 0.0)	Z_i[i].ve[0] = 0.0 - y_v_at(0, 0);
            if (dim > 1 && Z_i[i].ve[1] + y_0.ve[1] < 0.0)	Z_i[i].ve[1] = -y_0.ve[1];
            if (dim > 2 && Z_i[i].ve[2] + y_0.ve[2] < 0.0)	Z_i[i].ve[2] = -y_0.ve[2];
            if (dim > 2 && Z_i[i].ve[2] + y_0.ve[2] > 1.0)	Z_i[i].ve[2] = 1.0 - y_0.ve[2] - 0.0;
            if (dim > 3 && Z_i[i].ve[3] + y_0.ve[3] < 0.0)	Z_i[i].ve[3] = -y_0.ve[3];
            if (dim > 3 && Z_i[i].ve[3] + y_0.ve[3] > 1.0 - Z_i[i].ve[2] - y_0.ve[2])	Z_i[i].ve[3] = 1.0 - Z_i[i].ve[2] - y_0.ve[2] - y_0.ve[3];
        }
    }

    //Build the coefficient matrix, if necessary
    if (link_i->compute_J == 1)
        link_i->jacobian(t, y_0, temp_parent_approx[s], link_i->num_parents, GlobalVars, link_i->forcing_values, params, link_i->JMatrix);

    if (link_i->compute_LU == 1)
    {
        for (i = 0; i<s; i++) for (j = 0; j<s; j++)
        {
            for (l = 0; l<dim; l++) for (m = 0; m<dim; m++)
                CoefMat.me[i*dim + l][j*dim + m] = -h*A.me[i][j] * link_i->JMatrix.me[l][m];
        }
        for (i = 0; i<s*dim; i++)	CoefMat.me[i][i] += 1.0;
        info = clapack_dgetrf(CblasRowMajor, s*dim, s*dim, CoefMat->array, CoefMat->m, ipiv);
    }

    //Perform the simplified Newton iterations
    do
    {
        //Initialize RHS
        for (j = 0; j<s; j++) for (l = 0; l<dim; l++)
            RHS.ve[l + j*dim] = -Z_i[j].ve[l];

        //Evaluate the guess in the system created by the numerical method
        for (j = 0; j<s; j++)
        {
            v_add(y_0, Z_i[j], sum, 0);
            if (suv_at(m, 0) < 0.0)	suv_at(m, 0) = 0.0;	//A negative discharge is clearly a poor approximation
            if (dim > 1 && sum.ve[1] < 0.0)	sum.ve[1] = 0.0;
            if (dim > 2 && sum.ve[2] < 0.0)	sum.ve[2] = 0.0;
            if (dim > 2 && sum.ve[2] > 1.0)	sum.ve[2] = 1.0;
            if (dim > 3 && sum.ve[3] < 0.0)	sum.ve[3] = 0.0;
            if (dim > 3 && sum.ve[3] > 1.0 - sum.ve[2])	sum.ve[3] = 1.0 - sum.ve[2];
            link_i->f(t + v_at(c, j)*h, sum, temp_parent_approx[j], link_i->num_parents, GlobalVars, link_i->forcing_values, link_i->qvs, params, link_i->state, link_i->user, temp);
            for (l = 0; l<s; l++)
                daxpy(h*A.me[l][j], temp, RHS, l*dim);
        }

        //Make sure RHS has valid values
        for (j = 0; j<s; j++)
        {
            if (RHS.ve[j*dim] + (y_v_at(0, 0) + Z_i[j].ve[0]) < 0.0)	RHS.ve[j*dim] = -(y_v_at(0, 0) + Z_i[j].ve[0]);
            if (dim > 1 && RHS.ve[j*dim + 1] + (y_0.ve[1] + Z_i[j].ve[1]) < 0.0)	RHS.ve[j*dim + 1] = -(y_0.ve[1] + Z_i[j].ve[1]);
            if (dim > 2 && RHS.ve[j*dim + 2] + (y_0.ve[2] + Z_i[j].ve[2]) < 0.0)	RHS.ve[j*dim + 2] = -(y_0.ve[2] + Z_i[j].ve[2]);
            if (dim > 2 && RHS.ve[j*dim + 2] + (y_0.ve[2] + Z_i[j].ve[2]) > 1.0)	RHS.ve[j*dim + 2] = 1.0 - (y_0.ve[2] + Z_i[j].ve[2]);
            if (dim > 3 && RHS.ve[j*dim + 3] + (y_0.ve[3] + Z_i[j].ve[3]) < 0.0)	RHS.ve[j*dim + 3] = -(y_0.ve[3] + Z_i[j].ve[3]);
            if (dim > 3 && RHS.ve[j*dim + 3] + (y_0.ve[3] + Z_i[j].ve[3]) > 1.0 - RHS.ve[j*dim + 2] - (y_0.ve[2] + Z_i[j].ve[2]))	RHS.ve[j*dim + 3] = 1.0 - RHS.ve[j*dim + 2] - (y_0.ve[2] + Z_i[j].ve[2]) - (y_0.ve[3] + Z_i[j].ve[3]);
        }

        //Solve the system
        info = clapack_dgetrs(CblasRowMajor, 111, s*dim, 1, CoefMat->array, CoefMat->m, ipiv, RHS.ve, s*dim);

        //Form the solution to the iteration
        for (j = 0; j<s; j++)	v_add(RHS, Z_i[j], Z_i[j], j*dim);

        //Make sure Z_i has valid values
        for (j = 0; j<s; j++)
        {
            if (Z_i[j].ve[0] + y_v_at(0, 0) < 0.0)	Z_i[j].ve[0] = 0.0 - y_v_at(0, 0);
            if (dim > 1 && Z_i[j].ve[1] + y_0.ve[1] < 0.0)	Z_i[j].ve[1] = -y_0.ve[1];
            if (dim > 2 && Z_i[j].ve[2] + y_0.ve[2] < 0.0)	Z_i[j].ve[2] = -y_0.ve[2];
            if (dim > 2 && Z_i[j].ve[2] + y_0.ve[2] > 1.0)	Z_i[j].ve[2] = 1.0 - y_0.ve[2] - 0.0;
            if (dim > 3 && Z_i[j].ve[3] + y_0.ve[3] < 0.0)	Z_i[j].ve[3] = -y_0.ve[3];
            if (dim > 3 && Z_i[j].ve[3] + y_0.ve[3] > 1.0 - Z_i[j].ve[2] - y_0.ve[2])	Z_i[j].ve[3] = 1.0 - Z_i[j].ve[2] - y_0.ve[2] - y_0.ve[3];
        }

        //Check if solution is accurate enough
        normDeltam1 = normDelta;
        normDelta = vector_norminf(RHS);
        Theta = normDelta / normDeltam1;
        if (normDeltam1 < 0.0)	eta = link_i->last_eta;
        else			eta = Theta / (1.0 - Theta);

        //Check that the error estimation is roughly correct
        //error_est = pow(Theta,imax-i) / (1.0 - Theta) * normDelta;

        //Check abort conditions
        newt++;
        if (Theta >= 1.0 || newt == newtmax)	abort = 1;
        else			abort = 0;
    } while (eta*normDelta > 1e-10 && abort == 0);

    //Check that the nonlinear solver converged
    if (abort)
    {
        link_i->h = link_i->h / 2;
        if (link_i->ID == 10)	printf("Solver diverged. ID: %u Time: %f New Step: %.16f\n", link_i->ID, link_i->last_t, link_i->h);
        link_i->compute_J = 1;
        link_i->compute_LU = 1;
        return 0;
    }
    else
        link_i->last_eta = pow(eta, .8);

    //Check if jacobian should be recomputed
    if (normDeltam1 < 0.0 || Theta < 1e-3)
        link_i->compute_J = 0;
    else
        link_i->compute_J = 1;

    //Setup variables for the new data
    new_node = New_Step(&link_i->my->list);
    new_node->t = t + h;
    k = new_node->k;
    new_y = new_node->y_approx;

    //Build the solution
    v_add(y_0, Z_i[s - 1], new_y, 0);
    if (new_v_at(y, 0) < 0.0)	new_v_at(y, 0) = 1e-15;	//A negative discharge is clearly a poor approximation
    if (dim > 1 && new_y.ve[1] < 0.0)	new_y.ve[1] = 0.0;
    if (dim > 2 && new_y.ve[2] < 0.0)	new_y.ve[2] = 0.0;
    if (dim > 2 && new_y.ve[2] > 1.0)	new_y.ve[2] = 1.0 - 1e-15;
    if (dim > 3 && new_y.ve[3] < 0.0)	new_y.ve[3] = 1e-15;
    if (dim > 3 && sum.ve[3] > 1.0 - sum.ve[2])	sum.ve[3] = 1.0 - sum.ve[2];

    //Error estimation and step size selection

    //Compute LU factorization of (I - h e[s] JMatrix)
    dipaa(-h*e.ve[s], link_i->JMatrix, tempmat, 0, 0);
    info = clapack_dgetrf(CblasRowMajor, dim, dim, tempmat->array, dim, ipiv);

    //Build err for y_1
    v_copy(Z_i[0], sum);
    sv_mlt(v_at(e, 0), sum, 0);
    for (i = 1; i<s; i++)	daxpy(v_at(e, i), Z_i[i], sum, 0);
    v_copy(sum, err);
    link_i->f(t, y_0, temp_parent_approx[s], link_i->num_parents, GlobalVars, link_i->forcing_values, link_i->qvs, params, link_i->state, link_i->user, temp);	//Compute the extra k
    daxpy(h*e.ve[s], temp, err, 0);	//Extra k
    info = clapack_dgetrs(CblasRowMajor, 111, dim, 1, tempmat->array, dim, ipiv, err.ve, dim);

    if (link_i->rejected == 0)	//If the previous step was rejected
    {
        v_add(y_0, err, temp, 0);
        v_copy(sum, err);
        link_i->f(t, temp, temp_parent_approx[s], link_i->num_parents, GlobalVars, link_i->forcing_values, link_i->qvs, params, link_i->state, link_i->user, temp);
        daxpy(h*e.ve[s], temp, err, 0);
        info = clapack_dgetrs(CblasRowMajor, 111, dim, 1, tempmat->array, dim, ipiv, err.ve, dim);
    }

    //Build SC_i for y_1
    for (i = 0; i<dim; i++)
        temv_at(p, i) = max(fabs(new_v_at(y, i)), fabs(y_v_at(0, i))) * error->reltov_at(l, i) + error->abstov_at(l, i);

    //err_1 = norm_inf(sum,temp,meth->e_order_ratio,0);
    double err_1 = norm_inf(err, temp, 1.0, 0);
    double value_1 = pow(1.0 / err_1, 1.0 / meth->e_order);

    //Build err for dense output
    //	v_copy(Z_i[0],sum);
    //	sv_mlt(v_at(d, 0),sum,0);
    //	for(i=1;i<s;i++)	daxpy(v_at(d, i),Z_i[i],sum,0);
    //	v_copy(sum,err);
    //	link_i->f(t,y_0,temp_parent_approx[s],link_i->num_parents,GlobalVars,link_i->forcing_values,link_i->qvs,params,link_i->user,temp);	//Compute the extra k
    //	daxpy(h*d.ve[s],temp,err,0);	//Extra k
    //	info = clapack_dgetrs(CblasRowMajor,111,dim,1,JMatrix->array,dim,ipiv,err.ve,dim);	//This assumes d.ve[s] == e.ve[s]

    //	if(link_i->rejected == 0)	//If the previous step was rejected
    //	{
    //		v_add(y_0,err,temp,0);
    //		v_copy(sum,err);
    //		link_i->f(t,temp,temp_parent_approx[s],link_i->num_parents,GlobalVars,link_i->forcing_values,link_i->qvs,params,link_i->user,temp);
    //		daxpy(h*d.ve[s],temp,err,0);
    //		info = clapack_dgetrs(CblasRowMajor,111,dim,1,JMatrix->array,dim,ipiv,err.ve,dim);
    //	}

    //Build SC_i for dense
    //	for(i=0;i<dim;i++)
    //		temv_at(p, i) = max(fabs(new_v_at(y, i)),fabs(y_v_at(0, i))) * error->reltol_densv_at(e, i) + error->abstol_densv_at(e, i);

    //	double err_d = norm_inf(err,temp,1.0,0);
    //	double value_d = pow(1.0/err_d,1.0/meth->d_order);

    //Determine a new step size for the next step
    //double step_1 = h*min(error->facmax,max(error->facmin,error->fac * value_1));
    double step_1 = error->fac * (2 * newtmax + 1) / (2 * newtmax + newt) * h * value_1;
    double step_log;
    if (link_i->h_old > 0.0)
        step_log = error->fac * h * value_1 * h / link_i->h_old * value_1 / link_i->value_old;
    else
        step_log = 1e5;
    //	double step_d = h*min(error->facmax,max(error->facmin,error->fac * value_d));
    //	link_i->h = min(step_1,step_d);
    //	link_i->h = step_1;

    //	if(err_1 < 1.0 && err_d < 1.0)
    if (err_1 < 1.0)
    {
        //Check if an LU factorization should take place in the next step
        if (link_i->compute_J == 1 || h > step_1 || step_1 > 1.2*h)
        {
            link_i->h = min(step_1, step_log);
            link_i->compute_LU = 1;
        }
        else
            link_i->compute_LU = 0;

        //Save the new data
        for (j = 0; j<s; j++)
        {
            v_add(y_0, Z_i[j], sum, 0);
            if (suv_at(m, 0) < 0.0)	suv_at(m, 0) = 1e-10;	//A negative discharge is clearly a poor approximation
            if (dim > 1 && sum.ve[1] < 0.0)	sum.ve[1] = 0.0;
            if (dim > 2 && sum.ve[2] < 0.0)	sum.ve[2] = 0.0;
            if (dim > 2 && sum.ve[2] > 1.0)	sum.ve[2] = 1.0 - 1e-10;
            if (dim > 3 && sum.ve[3] < 0.0)	sum.ve[3] = 1e-10;
            if (dim > 3 && sum.ve[3] > 1.0 - sum.ve[2])	sum.ve[3] = 1.0 - sum.ve[2];
            link_i->f(t + v_at(c, j)*h, sum, temp_parent_approx[j], link_i->num_parents, GlobalVars, link_i->forcing_values, link_i->qvs, params, link_i->state, link_i->user, k[j]);
        }
        link_i->last_t = t + h;
        link_i->current_iterations++;
        link_i->h_old = h;
        link_i->value_old = value_1;
        v_sub(y_0, new_y, link_i->sol_diff, 0);
        for (i = 0; i<s; i++)	v_copy(Z_i[i], old_Z[i]);

        //Check if new data should be written to disk
        if (print_flag)
        {
            while (t <= link_i->next_save && link_i->next_save <= link_i->last_t)
            {
                //Don't write anything if using data assimilation and at a time when data is available
                if (GlobalVars->assim_flag)
                {
                    double rounded = 1e-13*rint(1e13*(GlobalVars->maxtime - link_i->next_save));
                    if (rounded < 1e-13 && -rounded < 1e-13)		break;
                }

                (link_i->disk_iterations)++;
                node = link_i->my->list.tail->prev;
                current_theta = (link_i->next_save - t) / h;
                link_i->method->dense_b(current_theta, link_i->method->b_theta);

                v_copy(node->y_approx, sum);
                for (l = 0; l<link_i->method->s; l++)
                    daxpy(h * link_i->method->b_theta.ve[l], node->next->k[l], sum, 0);
                if (suv_at(m, 0) < 0.0)	suv_at(m, 0) = 0.0;
                if (dim > 1 && sum.ve[1] < 0.0)	sum.ve[1] = 0.0;
                if (dim > 2 && new_y.ve[2] < 0.0)	new_y.ve[2] = 0.0;
                if (dim > 2 && new_y.ve[2] > 1.0)	new_y.ve[2] = 1.0;
                if (dim > 3 && new_y.ve[3] < 0.0)	new_y.ve[3] = 0.0;
                if (dim > 3 && sum.ve[3] > 1.0 - sum.ve[2])	sum.ve[3] = 1.0 - sum.ve[2];

                fsetpos(outputfile, &(link_i->pos));
                fwrite(&(link_i->next_save), sizeof(double), 1, outputfile);
                for (j = 0; j<num_print; j++)
                    fwrite(&(sum.ve[print_indices[j]]), sizeof(double), 1, outputfile);
                fgetpos(outputfile, &(link_i->pos));

                link_i->next_save += link_i->print_time;
            }
        }

        //Check if this is a max discharge
        if (new_v_at(y, 0) > link_i->peak_valuv_at(e, 0))
        {
            v_copy(new_y, link_i->peak_value);
            link_i->peak_time = link_i->last_t;
        }

        //Check if the newest step is on a change in rainfall
        short int propagated = 0;	//Set to 1 when last_t has been propagated
        for (j = 0; j<GlobalVars->num_forcings; j++)
        {
            if (link_i->forcing_buff[j] && fabs(link_i->last_t - link_i->forcing_change_times[j]) < 1e-14)
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
                //for(l=1;l<link_i->forcing_buff[j]->n_times;l++)
                for (l = link_i->forcing_indices[j] + 1; l<link_i->forcing_buff[j]->n_times; l++)
                    if (fabs(link_i->forcing_change_times[j] - link_i->forcing_buff[j]->rainfall[l][0]) < 1e-8)	break;
                link_i->forcing_indices[j] = l;

                double forcing_buffer = link_i->forcing_buff[j]->rainfall[l][1];
                link_i->forcing_values[j] = forcing_buffer;

                //Find and set the new change in rainfall
                for (i = l + 1; i<link_i->forcing_buff[j]->n_times; i++)
                {
                    if (link_i->forcing_buff[j]->rainfall[i][1] != forcing_buffer)
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

        link_i->h = min(step_1, step_log);
        link_i->compute_J = 1;
        link_i->compute_LU = 1;

        return 0;
    }

}

#endif // defined(ASYNCH_HAVE_IMPLICIT_SOLVER)