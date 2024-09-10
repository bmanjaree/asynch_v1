#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

//#include <stdio.h>
//#include <stdlib.h>
#include <math.h>
//#include <string.h>

//#include <io.h>
#include <minmax.h>
#include <rkmethods.h>
#include <blas.h>


//Copies contents of the vectors full_k [num_stages][max_dim] into the vector k [num_stages][num_dense]
void store_k(const double * const full_k, unsigned int num_dof, double *k, unsigned int num_stages, const unsigned int * const dense_indices, unsigned int num_dense)
{
    unsigned int i, j;

    // TODO optimize this loop
    for (i = 0; i < num_stages; i++)
    {
        for (j = 0; j < num_dense; j++)
            k[i * num_dense + j] = full_k[i * num_dof + dense_indices[j]];
    }
}


double InitialStepSize(double t, Link* link_i, const GlobalVars * const globals, Workspace* workspace)
{
    unsigned int start = link_i->diff_start;
    double *y_0 = link_i->my->list.tail->y_approx;
    double t_0 = link_i->my->list.tail->t;
    double *fy0 = workspace->temp;
    double *fy1 = workspace->sum;
    double *y_1 = workspace->temp2;
    double *SC = workspace->temp3;
    unsigned int p = link_i->method->localorder;
    RKSolutionNode* curr_node;
    double d0, d1, d2, h0, h1, largest, timediff, current_theta;
    unsigned int idx;
    unsigned int dim = link_i->dim;
    unsigned int num_dense = link_i->num_dense;
    unsigned int* dense_indices = link_i->dense_indices;
    double *parents_approx = workspace->parents_approx;
    ErrorData* error = link_i->my->error_data;
	bool from_reservoir;

    //Build SC for this link
    for (unsigned int i = 0; i < dim; i++)
        SC[i] = fabs(y_0[i]) * error->reltol[i] + error->abstol[i];

    //Grab parents data
    for (unsigned int i = 0; i < link_i->num_parents; i++)
    {
        Link *curr_parent = link_i->parents[i];
        curr_node = curr_parent->my->list.head;

        unsigned int num_stages = curr_parent->method->num_stages;
        unsigned int num_dense = curr_parent->num_dense;

        double *curr_parent_approx = workspace->parents_approx + i * globals->max_dim;

        if (fabs(curr_parent->my->list.tail->t - t) < 1e-10)	//!!!! Ugh... !!!!
            dcopy(curr_parent->my->list.tail->y_approx, curr_parent_approx, 0, curr_parent->dim);
        else
        {
            //Find the corresponding theta value and approximate solution
            while (t > curr_node->t)
                curr_node = curr_node->next;
            if (curr_node != curr_parent->my->list.head)
                curr_node = curr_node->prev;

            timediff = curr_node->next->t - curr_node->t;
            current_theta = (t - curr_node->t) / timediff;

            curr_parent->method->dense_b(current_theta, curr_parent->method->b_theta);

            // !!!! Note: this varies with num_print. Consider doing a linear interpolation. !!!!
            for (unsigned int m = 0; m < num_dense; m++)
            {
                idx = curr_parent->dense_indices[m];
                curr_parent_approx[idx] = curr_node->y_approx[idx];

                for (unsigned int l = 0; l < num_stages; l++)
                    curr_parent_approx[idx] += timediff * curr_parent->method->b_theta[l] * curr_node->next->k[l * num_dense + m];
            }
            link_i->check_consistency(curr_parent_approx, curr_parent->dim, globals->global_params, globals->num_global_params, curr_parent->params, link_i->num_params, curr_parent->user);

            if (link_i->algebraic)
                link_i->algebraic(curr_parent_approx, curr_parent->dim, globals->global_params, curr_parent->params, curr_parent->qvs, curr_parent->state, curr_parent->user, curr_parent_approx);
        }
    }

    //Step a
    //d0 = vector_norminf(y0,start);
    d0 = nrminf2(y_0, SC, start, link_i->dim);
    link_i->differential(
        t_0,
        y_0, link_i->dim,
        parents_approx, link_i->num_parents, globals->max_dim,
        globals->global_params, link_i->params, link_i->my->forcing_values, link_i->qvs, link_i->state, link_i->user, fy0);
    d1 = nrminf2(fy0, SC, start, link_i->dim);

    //Step b
    if (d0 < 1e-5 || d1 < 1e-5)
        h0 = 1e-6;
    else
        h0 = 0.01*(d0 / d1);

    //Step c
    //Note: This assumes the parents discharge is the same. It also assume no change in rain or state
    dcopy(y_0, y_1, 0, link_i->dim);
    daxpy(h0, fy0, y_1, start, link_i->dim);
    link_i->check_consistency(y_1, link_i->dim, globals->global_params, globals->num_global_params, link_i->params, link_i->num_params, link_i->user);
    link_i->differential(
        t_0 + h0,
        y_1, link_i->dim,
        parents_approx, link_i->num_parents, globals->max_dim,
        globals->global_params, link_i->params, link_i->my->forcing_values, link_i->qvs, link_i->state, link_i->user, fy1);

    //Step d
    dsub(fy1, fy0, fy1, start, link_i->dim);
    d2 = nrminf2(fy1, SC, start, link_i->dim);

    //Step e
    largest = max(d1, d2);
    h1 = (largest < 1e-1) ? max(1e-6, h0*1e-3) : pow(1e-2 / largest, 1.0 / (p + 1.0));

    //Step f
    // adlzanchetta: if has a parent reservoir, ignores h0
	from_reservoir = false;
	for (unsigned int i = 0; i < link_i->num_parents; i++)
	{
		Link *curr_parent = link_i->parents[i];
		if (curr_parent->has_res) {
			from_reservoir = true;
		}
	}
	if (!from_reservoir) {
		h1 = min(100.0*h0, h1);
	}

    //Make sure returned value is not too small. This can create roundoff problems when time gets larger
    //return h1;
    //frexp(h1,&exponent);
    //if(exponent < -29)	return 1e-10;
    if (y_0[0] < 1e-6)
        return 1e-2;
    else		
        return h1;
}
