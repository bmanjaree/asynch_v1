#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#if !defined(_MSC_VER)
#define ASYNCH_SLEEP sleep
#else
#include <windows.h>
#define ASYNCH_SLEEP Sleep
#endif

#if defined(HAVE_POSTGRESQL)
#include <libpq-fe.h>
#endif

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>


// Internal stuffs
#include <blas.h>
#include <db.h>
#include <comm.h>
#include <riversys.h>
#include <sort.h>

#include <assim/ancillary.h>


// For older version of OpenMPI
#if !defined(MPI_C_BOOL)
#define MPI_C_BOOL MPI_CHAR
#endif


//Finds the link ids upstreams from every link in obs_locs. If trim is 1, then only links which can affect the links in obs_locs (assuming a constant channel velocity) are used.
void FindUpstreamLinks(const AsynchSolver * const asynch, AssimData* const assim, unsigned int problem_dim, bool trim, double obs_time_step, unsigned int num_steps, unsigned int* obs_locs, unsigned int num_obs)
{
    Link *sys = asynch->sys, *current;
    unsigned int N = asynch->N, parentsval, leaves_size = 0, i, j;
    Lookup *id_to_loc = asynch->id_to_loc;
    int *assignments = asynch->assignments;
    GlobalVars *globals = asynch->globals;

    short int* getting = asynch->getting;
    UpstreamData* updata;

    //For every links
    for (i = 0; i < N; i++) {
        //Allocate UpstreamData
        UpstreamData *upstreams = malloc(sizeof(UpstreamData));
        memset(upstreams, 0, sizeof(UpstreamData));

        //Copy the vector of parents list
        upstreams->num_parents = sys[i].num_parents;
        if (upstreams->num_parents > 0)
        {
            upstreams->parents = malloc(upstreams->num_parents * sizeof(Link*));
            memcpy(upstreams->parents, sys[i].parents, upstreams->num_parents * sizeof(Link*));
        }

        sys[i].user = upstreams;
    }

    //Find leaves
    Link **leaves = (Link**)malloc(N * sizeof(Link*));
    for (i = 0; i < N; i++)
        if (sys[i].num_parents == 0)
            leaves[leaves_size++] = &sys[i];

    //Set the number of upstreams parents to 0 for leaves
    for (i = 0; i < leaves_size; i++)
        ((UpstreamData*)leaves[i]->user)->num_upstreams = 0;

    //Count upstreams links
    for (i = 0; i < leaves_size; i++)
    {
        //Traverse the tree downstream the leave
        for (current = leaves[i]->child; current != NULL; current = current->child)
        {
            parentsval = 0;
            for (j = 0; j < current->num_parents; j++)
            {
                Link *parent = current->parents[j];
                UpstreamData *parent_updata = (UpstreamData*)parent->user;
                //parentsval += (temp_numupstream[current->parents[j]->location] > 0);
                parentsval += (parent_updata->num_upstreams > 0) || (parent->num_parents == 0);
            }

            if (parentsval == current->num_parents)	//All parents have upstreams set
            {
                //temp_numupstream[current->location] = 1;
                ((UpstreamData*)current->user)->num_upstreams = current->num_parents;
                for (j = 0; j < current->num_parents; j++)
                    //temp_numupstream[current->location] += temp_numupstream[current->parents[j]->location];
                    ((UpstreamData*)current->user)->num_upstreams += ((UpstreamData*)current->parents[j]->user)->num_upstreams;
            }
            else
                break;
        }
    }

    //Set the upstreams links
    //unsigned int** temp_upstream = (unsigned int**) malloc(N*sizeof(unsigned int*));	//temp_upstream[i] is a list of all links upstreams from link i
    //for (i = 0; i < N; i++)
    //    temp_upstream[i] = (unsigned int*)malloc(temp_numupstream[sys[i].location] * sizeof(unsigned int));

    for (i = 0; i < N; i++)
    {
        int num_upstreams = ((UpstreamData*)sys[i].user)->num_upstreams;
        if (num_upstreams > 0)
            ((UpstreamData*)sys[i].user)->upstreams = calloc(num_upstreams, sizeof(Link*));
    }

    unsigned int* counter = (unsigned int*)calloc(N, sizeof(unsigned int));

    unsigned int stack_size = leaves_size;
    Link **stack = (Link**)calloc(N, sizeof(Link*));
    for (i = 0; i < leaves_size; i++)
        stack[i] = leaves[i];

    while (stack_size > 0)
    {
        unsigned int l;

        current = stack[stack_size - 1];
        updata = (UpstreamData*)current->user;

        l = current->location;

        ////Add this link to its own upstreams list
        //updata->upstreams[counter[l]] = current;
        //counter[l]++;

        //Add each parents' upstreams list
        unsigned int count = 0;
        for (i = 0; i < current->num_parents; i++)
        {
            Link *parent = current->parents[i];
            if (parent->has_res)
                continue;

            assert(count < updata->num_upstreams);
            updata->upstreams[count] = parent;
            count++;

            UpstreamData *parent_updata = ((UpstreamData*)parent->user);
            for (j = 0; j < parent_updata->num_upstreams; j++)
            {
                assert(count < updata->num_upstreams);
                updata->upstreams[count] = parent_updata->upstreams[j];
                count++;
            }
        }

        assert(count == updata->num_upstreams);
        counter[l] = count;

        stack_size--;

        //If every parent of current's child has an upstreams list determined, add it to the stack
        if (current->child != NULL)
        {
            parentsval = 0;
            for (i = 0; i < current->child->num_parents; i++)
            {
                Link *p = current->child->parents[i];
                parentsval += (counter[p->location] > 0) || (p->num_parents == 0);
            }

            if (parentsval == current->child->num_parents)
            {
                stack[stack_size] = current->child;
                stack_size++;
            }
        }
    }

    //	//Move the data from temp_upstream into the child upstreams
    //	short int* used = (short int*) calloc(N,sizeof(short int));	//1 if temp_upstream[i] was used, 0 if not
    //	for(i=0;i<N;i++)
    //	{
    //		//if(assignments[i] == my_rank || getting[i])
    //		{
    //			sys[i].user = malloc(sizeof(UpstreamData));
    //			updata = (UpstreamData*) (sys[i].user);
    //
    //			updata->upstreams = (unsigned int**) malloc(sys[i].num_parents * sizeof(unsigned int*));
    //			updata->num_upstreams = (unsigned int*) malloc(sys[i].num_parents * sizeof(unsigned int));
    //			updata->fit_states = NULL;
    //			updata->fit_to_universal = NULL;
    //			for(j=0;j<sys[i].num_parents;j++)
    //			{
    //				updata->upstreams[j] = temp_upstream[sys[i].parents[j]->location];
    //				updata->num_upstreams[j] = temp_numupstream[sys[i].parents[j]->location];
    //				used[sys[i].parents[j]->location] = 1;
    //			}
    ///*
    //			//Calculate the dimension at each link
    //			updata->dim = problem_dim + problem_dim + (problem_dim-1)*(problem_dim-1);	//Model eqs + variational eqs from this link
    //			for(j=0;j<sys[i].num_parents;j++)
    //				updata->dim += updata->num_upstreams[j] * problem_dim;	//Variational eqs from upstreams !!!! Too high? !!!!
    //*/
    //		}
    //	}


        //Cleanup
    //for (i = 0; i < N; i++)
    //{
    //    if (!used[i])	free(temp_upstream[i]);
    //}
    //free(temp_upstream);
    //free(temp_numupstream);
    free(counter);
    free(leaves);
    free(stack);
    //free(used);

    //Remove extra links from the upstreams lists
    //!!!! This needs to be generalized for the case of many outlets !!!!
    if (trim)	//!!!! Does this work correctly for the entire state? !!!!
    {
        unsigned int loc;
        PGresult *res;
        double *distance = (double*)calloc(N, sizeof(double));
        double influence_radius;
        double speed = 3.0 * 60.0;	//In m/min

        //!!!! Hard coding right now. Blah... !!!!
        unsigned int outlet = asynch->globals->outletlink;
        unsigned int n = 0;
        //unsigned int outlet = 434478;  //Turkey River above French Hollow
        //unsigned int outlet = 434514;	//Turkey River at Garber
        //unsigned int outlet = 307864;  //Half Squaw Creek
        //unsigned int outlet = 292254;	//Squaw Creek at Ames

        if (my_rank == 0)
        {
            ConnectPGDB(&assim->conninfo);

            char buffer[16];
            snprintf(buffer, 16, "%d", outlet);
            const char *paramValues[1];
            paramValues[0] = buffer;

            res = PQexecParams(assim->conninfo.conn, assim->conninfo.queries[2], 1, NULL, paramValues, NULL, NULL, 0);
            if (CheckResError(res, "getting list of distances to outlet"))
                MPI_Abort(MPI_COMM_WORLD, 1);
            n = PQntuples(res);
            if (n != N)
            {
                printf("Error: got a different number of links for the distances to outlet than links in network. (%u vs %u)\n", i, N);
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            else
            {
                int i_link_id = PQfnumber(res, "link_id");
                int i_distance = PQfnumber(res, "distance");

                //Sort the data
                for (i = 0; i < N; i++)
                {
                    if (!PQgetisnull(res, i, i_link_id) && PQgetlength(res, i, i_link_id) > 0)
                    {
                        char *ptr = PQgetvalue(res, i, 0);
                        //int link_id = be32toh(*((uint32_t *)ptr));
                        int link_id = atoi(ptr);

                        loc = find_link_by_idtoloc(link_id, id_to_loc, N);
                        if (loc >= N)
                        {
                            i = loc;
                            break;
                        }

                        if (!PQgetisnull(res, i, i_distance) && PQgetlength(res, i, i_distance) > 0)
                        {
                            ptr = PQgetvalue(res, i, i_distance);
                            //distance[loc] = (double) (be64toh(*((uint64_t *)ptr)));
                            distance[loc] = atof(ptr);
                        }
                    }
                }
            }

            //Clean up db connection
            PQclear(res);
            DisconnectPGDB(&assim->conninfo);

            //printf("!!!! Loading distances for test basin !!!!\n");
            //distance[0] = 0;distance[1] = 1;distance[2] = 2;distance[3] = 2;distance[4] = 3;distance[5] = 1;distance[6] = 3;distance[7] = 4;distance[8] = 2;distance[9] = 3;distance[10] = 3;
        }
        //speed = 0.021;

        MPI_Bcast(distance, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        //Calculate the radius of influence
        influence_radius = speed * num_steps * obs_time_step;

        //For each link, find the gauges influenced. Then, check which upstreams links influence one of those gauges.
        unsigned int *influenced_gauges = (unsigned int*)calloc(num_obs, sizeof(unsigned int));
        unsigned int num_influenced, drop;
        double difference;
        for (i = 0; i < N; i++)
        {
            current = &sys[i];
            updata = (UpstreamData*)current->user;

            //Find influenced gauges
            num_influenced = 0;
            for (j = 0; j < num_obs; j++)
            {
                difference = distance[i] - distance[obs_locs[j]];
                if (-1e-12 < difference && difference < influence_radius)	//!!!! Does this work if the network is disconnected? !!!!
                {
                    influenced_gauges[num_influenced++] = obs_locs[j];
                    //if(sys[obs_locs[j]]->ID == 399711)
                    //printf("%u\n",current->ID);
                }
            }

            //Find which upstreams links are also influenced
            //for (j = 0; j < current->num_parents; j++)
            //{
            unsigned int *num_upstreams = &updata->num_upstreams;
            Link **upstreams = updata->upstreams;

            drop = 0;
            for (unsigned int l = 0; l < *num_upstreams; l++)
            {
                unsigned int m;

                for (m = 0; m < num_influenced; m++)
                {
                    if (distance[upstreams[l]->location] - distance[influenced_gauges[m]] < influence_radius)
                        break;
                }

                if (m == num_influenced)	//Upstream link does not influence any gauged location
                    drop++;
                //else
                //    upstreams[l - drop] = upstreams[l];
            }
            //*num_upstreams -= drop;
            //if (drop)
            //{
            //    if (*num_upstreams)
            //        upstreams = (Link**)realloc(upstreams, *num_upstreams * sizeof(Link*));
            //    else {
            //        free(upstreams);
            //        upstreams = NULL;
            //    }
            //}
            //}
        }

        //For each gauges
        stack_size = 0;
        Link **stack = (Link**)calloc(N, sizeof(Link*));
        for (i = 0; i < num_obs; i++)
        {
            Link *gauge = &sys[obs_locs[i]];
            double d = distance[obs_locs[i]];

            // Visit from gauge to source
            while (stack_size > 0)
            {
                Link *current = stack[stack_size - 1];
                UpstreamData* updata = (UpstreamData*)current->user;

                // Pop from the stack
                stack_size--;

                double difference = distance[current->location] - distance[gauge->location];
                //if (-1e-12 < difference && difference < influence_radius)	//!!!! Does this work if the network is disconnected? !!!!

                // If upstream and in influence radius
                if (difference >= 0. && difference < influence_radius)



                    for (unsigned int i = 0; i < updata->num_parents; i++)
                    {
                        stack[stack_size] = &sys[i];
                        stack_size++;
                    }
            }


        }



        //For each link, find the gauges influenced. Then, check which upstreams links influence one of those gauges.
        //Link **stack = (Link**)calloc(N, sizeof(Link*));
        stack_size = 0;
        for (i = 0; i < N; i++)
            if (sys[i].num_parents == 0)
            {
                stack[stack_size] = &sys[i];
                stack_size++;
            }

        // Visit from source to outlet
        unsigned int *visits = (unsigned int *)calloc(N, sizeof(unsigned int));
        while (stack_size > 0)
        {
            Link *current = stack[stack_size - 1];
            UpstreamData* updata = (UpstreamData*)current->user;

            // Pop from the stack
            stack_size--;

            // Increment visit counter of child
            if (current->child)
                visits[current->child->location]++;

            //Find gauges that influence the current link
            unsigned int num_influenced = 0;
            for (j = 0; j < num_obs; j++)
            {
                double difference = distance[current->location] - distance[obs_locs[j]];
                //if (-1e-12 < difference && difference < influence_radius)	//!!!! Does this work if the network is disconnected? !!!!

                // If upstream and in influence radius
                if (difference >= 0. && difference < influence_radius)
                    //influenced_gauges[num_influenced++] = obs_locs[j];
                    num_influenced++;
            }

            if (num_influenced == 0)
            {
                updata->num_upstreams = 0;

                free(updata->upstreams);
                updata->upstreams = NULL;
            }
            else
            {
                //Compute the number of upstreams link
                unsigned int num_upstreams = updata->num_parents;
                for (unsigned int l = 0; l < updata->num_parents; l++)
                {
                    Link *parent = current->parents[l];
                    UpstreamData *parent_updata = (UpstreamData*)parent->user;

                    num_upstreams += parent_updata->num_upstreams;
                }

                assert(num_upstreams <= N);

                //If this link neeeds to be fixed
                if (updata->num_upstreams != num_upstreams)
                {
                    updata->num_upstreams = num_upstreams;

                    //Reallocate the data
                    if (num_upstreams > 0)
                    {
                        updata->upstreams = (Link**)realloc(updata->upstreams, num_upstreams * sizeof(Link*));
                        assert(updata->upstreams != NULL);
                    }
                    else
                    {
                        free(updata->upstreams);
                        updata->upstreams = NULL;
                    }

                    //Add each parents' upstreams list
                    unsigned int count = 0;
                    for (i = 0; i < current->num_parents; i++)
                    {
                        Link *parent = current->parents[i];
                        if (parent->has_res)
                            continue;

                        assert(count < updata->num_upstreams);
                        updata->upstreams[count] = parent;
                        count++;

                        UpstreamData *parent_updata = ((UpstreamData*)parent->user);
                        for (j = 0; j < parent_updata->num_upstreams; j++)
                        {
                            assert(count < updata->num_upstreams);
                            updata->upstreams[count] = parent_updata->upstreams[j];
                            count++;
                        }
                    }

                    assert(count == updata->num_upstreams);
                }
            }

            if (current->child && visits[current->child->location] == current->child->num_parents)
            {
                stack[stack_size] = current->child;
                stack_size++;
            }
        }

        for (i = 0; i < N; i++)
        {
            current = &sys[i];
            updata = (UpstreamData*)current->user;

            unsigned int num_upstreams = updata->num_parents;

            for (unsigned int l = 0; l < updata->num_parents; l++)
            {
                Link *parent = current->parents[l];
                UpstreamData *parent_updata = (UpstreamData*)parent->user;

                num_upstreams += parent_updata->num_upstreams;
            }

            if (num_upstreams != updata->num_upstreams)
                printf("Olala");

            if (num_upstreams > 0 && updata->upstreams == NULL)
                printf("Olala");
        }


        /*
        for(i=0;i<N;i++)
        {
        current = sys[i];
        updata = (UpstreamData*) current->user;
        printf("ID = %u loc = %u\n",current->ID,i);
        for(j=0;j<current->num_parents;j++)
        {
        printf("upstreams: %u\n",updata->num_upstreams[j]);
        for(l=0;l<updata->num_upstreams[j];l++)
            printf("%u ",updata->upstreams[j][l]);
        printf("\n");
        }
        printf("+++++++\n");
        }
        */

        //Clean up
        free(stack);
        free(visits);
        free(distance);
    }
}


//Finds the link ids upstreams from every link in obs_locs. If trim is 1, then only links which can affect the links in obs_locs (assuming a constant channel velocity) are used.
void FindUpstreamLinks2(const AsynchSolver * const asynch, AssimData* const assim, unsigned int problem_dim, bool trim, double obs_time_step, unsigned int num_steps, unsigned int* obs_locs, unsigned int num_obs)
{
    Link *sys = asynch->sys;
    unsigned int N = asynch->N;
    Lookup *id_to_loc = asynch->id_to_loc;
    int *assignments = asynch->assignments;
    GlobalVars *globals = asynch->globals;

    short int* getting = asynch->getting;
    //UpstreamData* updata;

    //For every links
    for (unsigned int i = 0; i < N; i++)
    {
        //Allocate UpstreamData
        UpstreamData *updata = malloc(sizeof(UpstreamData));
        memset(updata, 0, sizeof(UpstreamData));

        //Copy the vector of parents list
        updata->num_parents = sys[i].num_parents;
        if (updata->num_parents > 0)
        {
            updata->parents = malloc(updata->num_parents * sizeof(Link*));
            memcpy(updata->parents, sys[i].parents, updata->num_parents * sizeof(Link*));
        }

        assert(updata != NULL);
        sys[i].user = updata;
    }

    //Remove extra links from the upstreams lists
    double *distance = (double*)calloc(N, sizeof(double));

    //!!!! Hard coding right now. Blah... !!!!
    unsigned int outlet = asynch->globals->outletlink;
    unsigned int n = 0;
    //unsigned int outlet = 434478;  //Turkey River above French Hollow
    //unsigned int outlet = 434514;	//Turkey River at Garber
    //unsigned int outlet = 307864;  //Half Squaw Creek
    //unsigned int outlet = 292254;	//Squaw Creek at Ames

    if (my_rank == 0)
    {
        ConnectPGDB(&assim->conninfo);

        char buffer[16];
        snprintf(buffer, 16, "%d", outlet);
        const char *paramValues[1];
        paramValues[0] = buffer;

        PGresult *res = PQexecParams(assim->conninfo.conn, assim->conninfo.queries[2], 1, NULL, paramValues, NULL, NULL, 0);
        if (CheckResError(res, "getting list of distances to outlet"))
            MPI_Abort(MPI_COMM_WORLD, 1);
        n = PQntuples(res);
        if (n != N)
        {
            printf("Error: got a different number of links for the distances to outlet than links in network. (%u vs %u)\n", n, N);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        else
        {
            int i_link_id = PQfnumber(res, "link_id");
            int i_distance = PQfnumber(res, "distance");

            //Sort the data
            for (unsigned int i = 0; i < N; i++)
            {
                if (!PQgetisnull(res, i, i_link_id) && PQgetlength(res, i, i_link_id) > 0)
                {
                    char *ptr = PQgetvalue(res, i, 0);
                    //int link_id = be32toh(*((uint32_t *)ptr));
                    int link_id = atoi(ptr);

                    unsigned int loc = find_link_by_idtoloc(link_id, id_to_loc, N);
                    if (loc >= N)
                    {
                        i = loc;
                        break;
                    }

                    if (!PQgetisnull(res, i, i_distance) && PQgetlength(res, i, i_distance) > 0)
                    {
                        ptr = PQgetvalue(res, i, i_distance);
                        //distance[loc] = (double) (be64toh(*((uint64_t *)ptr)));
                        distance[loc] = atof(ptr);
                    }
                }
            }
        }

        //Clean up db connection
        PQclear(res);
        DisconnectPGDB(&assim->conninfo);

        //printf("!!!! Loading distances for test basin !!!!\n");
        //distance[0] = 0;distance[1] = 1;distance[2] = 2;distance[3] = 2;distance[4] = 3;distance[5] = 1;distance[6] = 3;distance[7] = 4;distance[8] = 2;distance[9] = 3;distance[10] = 3;
    }
    //speed = 0.021;

    MPI_Bcast(distance, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //Calculate the radius of influence
    double speed = 3.0 * 60.0;	//In m/min
    double influence_radius = speed * num_steps * obs_time_step;

    //For each gauges
    Link **stack = (Link**)calloc(N, sizeof(Link*));
    Link **sources = (Link**)calloc(N, sizeof(Link*));
    unsigned int *visits = (unsigned int *)calloc(N, sizeof(unsigned int));
    for (unsigned int i = 0; i < num_obs; i++)
    {
        Link *gauge = &sys[obs_locs[i]];
        double d = distance[obs_locs[i]];

        unsigned int stack_size = 0;
        unsigned int source_size = 0;
        memset(visits, 0, N * sizeof(unsigned int));

        stack[stack_size++] = gauge;

        unsigned int count_phase1 = 0;

        // Visit from gauge to source
        while (stack_size > 0)
        {
            Link *current = stack[stack_size - 1];
            UpstreamData* updata = (UpstreamData*)current->user;

            // Pop from the stack
            stack_size--;

            count_phase1++;

            double difference = distance[current->location] - distance[gauge->location];
            assert(difference >= 0.);
            //if (-1e-12 < difference && difference < influence_radius)	//!!!! Does this work if the network is disconnected? !!!!

            //Check if this link is an upstream gauge
            bool is_gauge = false;
            for (unsigned int j = 0; j < num_obs; j++)
                if ((i != j) && (current->location == obs_locs[j]))
                {
                    is_gauge = true;
                    break;
                }

            // If upstream and in influence radius
            if (!is_gauge && (difference < influence_radius))
            {
                //assert(updata->num_upstreams == 0);
                updata->num_upstreams = 1;

                // This is a source
                if (current->num_parents == 0)
                    sources[source_size++] = current;
                else
                    //Push the parents links to the queue
                    for (unsigned int j = 0; j < current->num_parents; j++)
                        stack[stack_size++] = current->parents[j];
            }
            else
                //Else this is a source (technically a source of the domaine of influence)
                sources[source_size++] = current;
        }

        unsigned int count_phase2 = 0;

        // Visit from source to outlet
        while (source_size > 0)
        {
            Link *current = sources[source_size - 1];
            UpstreamData* updata = (UpstreamData*)current->user;

            // Pop from the source stack
            source_size--;

            count_phase2++;

            // Increment visit counter of child
            if (current->child)
                visits[current->child->location]++;

            //Compute the number of upstreams link
            unsigned int num_upstreams = updata->num_parents;
            for (unsigned int l = 0; l < updata->num_parents; l++)
            {
                Link *parent = current->parents[l];
                UpstreamData *parent_updata = (UpstreamData*)parent->user;

                num_upstreams += parent_updata->num_upstreams;
            }

            updata->num_upstreams = num_upstreams;

            //Reallocate the data
            if (num_upstreams > 0)
            {
                updata->upstreams = (Link**)calloc(num_upstreams, sizeof(Link*));
                assert(updata->upstreams != NULL);

                //Add each parents' upstreams list
                unsigned int count = 0;
                for (unsigned int i = 0; i < current->num_parents; i++)
                {
                    Link *parent = current->parents[i];
                    if (parent->has_res)
                        continue;

                    assert(count < updata->num_upstreams);
                    updata->upstreams[count] = parent;
                    count++;

                    UpstreamData *parent_updata = ((UpstreamData*)parent->user);
                    for (unsigned int j = 0; j < parent_updata->num_upstreams; j++)
                    {
                        assert(count < updata->num_upstreams);
                        updata->upstreams[count] = parent_updata->upstreams[j];
                        count++;
                    }
                }

                assert(count == updata->num_upstreams);
            }

            assert(!(updata->num_upstreams > 0 && updata->upstreams == NULL));

            if (current != gauge && current->child && visits[current->child->location] == current->child->num_parents)
                sources[source_size++] = current->child;
        }

        assert(count_phase1 == count_phase2);
    }

    //Clean up
    free(stack);
    free(sources);
    free(visits);
    free(distance);
}


//Deletes any unneeded upstreams link information. Use this after system partitioning is determined.
void CleanUpstreamLinks(const AsynchSolver* asynch)
{
    Link* sys = asynch->sys;
    unsigned int N = asynch->N, i;
    int* assignments = asynch->assignments;
    short int* getting = asynch->getting;

    for (i = 0; i < N; i++)
    {
        if (assignments[i] != my_rank && !getting[i])
        {
            UpstreamData *data = (UpstreamData*)(sys[i].user);
            if (data)
            {
                //for (j = 0; j < sys[i].num_parents; j++)
                //    if (data->upstreams[j])	free(data->upstreams[j]);
                if (data->fit_states)
                    free(data->fit_states);
                if (data->fit_to_universal)
                    free(data->fit_to_universal);
                if (data->upstreams)
                    free(data->upstreams);
                if (data->parents)
                    free(data->parents);
                free(data);
                sys[i].user = NULL;
            }
        }
    }
}

void FreeUpstreamLinks(const AsynchSolver* asynch)
{
    Link* sys = asynch->sys;
    unsigned int N = asynch->N, i;

    for (i = 0; i < N; i++)
    {
        UpstreamData *data = (UpstreamData*)(sys[i].user);
        if (data)
        {
            //for (j = 0; j < sys[i].num_parents; j++)
            //    if (data->upstreams[j])	free(data->upstreams[j]);
            if (data->fit_states)
                free(data->fit_states);
            if (data->fit_to_universal)
                free(data->fit_to_universal);
            if (data->upstreams)
                free(data->upstreams);
            if (data->parents)
                free(data->parents);
            free(data);
            sys[i].user = NULL;
        }
    }
}

////Read into memory the times and discharges stored in a .dat file.
//double*** ReadSolution(char filename[], unsigned int** id_to_loc, unsigned int N, unsigned int* numlinks, unsigned int** ids, unsigned int** locs, unsigned int** numsteps)
//{
//    unsigned int i, j, k, dim;
//    double*** data = NULL;
//
//    if (my_rank == 0)
//    {
//        FILE* file = fopen(filename, "r");
//        if (!file)
//        {
//            printf("Error reading true solution. File %s not found.\n", filename);
//            abort();
//        }
//
//        //Setup data structures
//        fscanf(file, "%u%u", numlinks, &dim);
//        *ids = (unsigned int*)malloc(*numlinks * sizeof(unsigned int));
//        *locs = (unsigned int*)malloc(*numlinks * sizeof(unsigned int));
//        *numsteps = (unsigned int*)calloc(*numlinks, sizeof(unsigned int));
//        data = (double***)malloc(*numlinks * sizeof(double**));
//
//        //Read in the file
//        for (i = 0; i < *numlinks; i++)
//        {
//            fscanf(file, "%u %u", &((*ids)[i]), &((*numsteps)[i]));
//            data[i] = (double**)malloc((*numsteps)[i] * sizeof(double*));
//            for (j = 0; j < (*numsteps)[i]; j++)
//            {
//                data[i][j] = (double*)malloc(2 * sizeof(double));
//                fscanf(file, "%lf %lf", &(data[i][j][0]), &(data[i][j][1]));	//time and discharge
//                for (k = 2; k < dim; k++)	fscanf(file, "%*f");
//            }
//        }
//
//        //Find locations from ids
//        for (i = 0; i < *numlinks; i++)
//            (*locs)[i] = find_link_by_idtoloc((*ids)[i], id_to_loc, N);
//
//        //Cleanup
//        fclose(file);
//    }
//
//    /*
//        //Send data to all procs
//        MPI_Bcast(numlinks,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
//        if(my_rank != 0)
//        {
//            *ids = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
//            *locs = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
//            *numsteps = (unsigned int*) malloc(*numlinks*sizeof(unsigned int));
//            data = (double***) malloc(*numlinks*sizeof(double**));
//        }
//        MPI_Bcast(*ids,*numlinks,MPI_UNSIGNED,0,MPI_COMM_WORLD);
//        MPI_Bcast(*locs,*numlinks,MPI_UNSIGNED,0,MPI_COMM_WORLD);
//        MPI_Bcast(*numsteps,*numlinks,MPI_UNSIGNED,0,MPI_COMM_WORLD);
//
//        for(i=0;i<*numlinks;i++)
//        {
//            if(my_rank != 0)
//                data[i] = (double**) malloc((*numsteps)[i]*sizeof(double*));
//            for(j=0;j<(*numsteps)[i];j++)
//            {
//                if(my_rank != 0)	data[i][j] = (double*) malloc(2*sizeof(double));
//                MPI_Bcast(data[i][j],2,MPI_DOUBLE,0,MPI_COMM_WORLD);
//            }
//        }
//    */
//
//    MPI_Bcast(numlinks, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//    if (my_rank != 0)
//    {
//        *numsteps = NULL;
//        *ids = (unsigned int*)malloc(*numlinks * sizeof(unsigned int));
//        *locs = (unsigned int*)malloc(*numlinks * sizeof(unsigned int));
//    }
//    MPI_Bcast(*ids, *numlinks, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//    MPI_Bcast(*locs, *numlinks, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//
//    return data;
//}


////Returns the discharges for all links with observations in data at time t.
//void FindAllDischarges(double*** data, double t, unsigned int numlinks, unsigned int* numsteps, double* d)
//{
//    unsigned int i, j, min, max;
//
//    if (my_rank == 0)
//    {
//        for (j = 0; j < numlinks; j++)
//        {
//            i = numsteps[j] / 2;
//            max = numsteps[j];
//            min = 0;
//            if (t > data[j][numsteps[j] - 1][0])
//                i = max - 1;
//            else
//            {
//                while (t*.99999 > data[j][i][0] || t*1.00001 < data[j][i][0])
//                {
//                    if (data[j][i][0] < t)	min = i + 1;
//                    else			max = i;
//                    i = (max + min) / 2;
//
//                    if (min >= max)
//                    {
//                        printf("Time %f not found for data %u.\n", t, j);
//                        break;
//                    }
//                }
//            }
//
//            d[j] = data[j][i][1];
//        }
//    }
//
//    MPI_Bcast(d, numlinks, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//}

//Gives a list of all locations with a downstream gauge, sorted by location.
//This allocates an appropriate amount of space for above_gauges and is_above_gauges.
//is_above_gauges[i] == 1 if sys[i] is above a gauge. 0 if not.
//above_gauges[i] is the list of link locations above a gauge.
//Returns the total number of gauges above a gauge.
unsigned int GaugeDownstream(const AsynchSolver* asynch, const unsigned int* obs_locs, unsigned int num_obs, unsigned int** above_gauges, bool **is_above_gauges)
{
    if (!above_gauges || !is_above_gauges)	return 0;
    unsigned int N = asynch->N;
    unsigned int i, j;
    unsigned int num_above = 0;
    Link *current, *next, *sys = asynch->sys;
    *is_above_gauges = (bool*)calloc(N, sizeof(bool));

    if (my_rank == 0)
    {
        //Set the location with gauges
        for (i = 0; i < num_obs; i++)
            (*is_above_gauges)[obs_locs[i]] = 1;
        num_above = num_obs;

        //Trickle down from each external link, until a gauge is reached, marking all links
        for (i = 0; i < N; i++)
        {
            // For all source links
            if (sys[i].num_parents == 0)
            {
                current = &sys[i];
                next = current;

                //See if there is a gauge
                while (next && (*is_above_gauges)[next->location] == 0)	next = next->child;

                if (next)	//Gauge found
                {
                    for (; current != next; current = current->child)
                    {
                        (*is_above_gauges)[current->location] = 1;
                        num_above++;
                    }
                }
            }
        }

        //Setup above_gauges
        *above_gauges = (unsigned int*)malloc(num_above * sizeof(unsigned int));
        j = 0;
        for (i = 0; i < N; i++)
        {
            if ((*is_above_gauges)[i])
                (*above_gauges)[j++] = i;
        }
    }

    MPI_Bcast(*is_above_gauges, N, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&num_above, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    if (my_rank != 0)
        *above_gauges = (unsigned int*)malloc(num_above * sizeof(unsigned int));
    MPI_Bcast(*above_gauges, num_above, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    return num_above;
}

//Scales discharge upstreams by area from the gauges
//!!!! Going to need areas when calling this with multiple procs !!!!
int AdjustDischarges(const AsynchSolver* asynch, const unsigned int* obs_locs, const double * obs, unsigned int num_obs, unsigned int problem_dim, double* x)
{
    //Unpack
    GlobalVars* globals = asynch->globals;
    unsigned int N = asynch->N;
    Lookup *id_to_loc = asynch->id_to_loc;
    int *assignments = asynch->assignments;
    Link *sys = asynch->sys;
    unsigned int area_idx = globals->area_idx;

    unsigned int i, j, num_upstreams;
    unsigned int *upstreams = (unsigned int*)malloc(N * sizeof(unsigned int));
    short int *locs_set = (short int*)calloc(N, sizeof(short int));	//1 if the discharge is to be changed, 0 if not
    double *locs_newq = (double*)malloc(N * sizeof(double));	//The new value for the discharge at location i
    Link* current;
    double ratio;

    //Get the upstreams area for each link
    double* upareas = (double*)calloc(N, sizeof(double));
    for (i = 0; i < N; i++)
    {
        if (assignments[i] == my_rank)
            upareas[i] = sys[i].params[globals->area_idx];
    }

    MPI_Allreduce(MPI_IN_PLACE, upareas, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    //Initialize the new discharges to a negative value
    for (i = 0; i < N; i++)
        locs_newq[i] = -1.0;

    //Store gauge readings
    for (i = 0; i < num_obs; i++)
    {
        locs_set[obs_locs[i]] = 1;
        locs_newq[obs_locs[i]] = obs[i];
    }

    //Perform the trickle down one last time. This is for links with no upstreams or downstream gauges.
    for (i = 0; i < N; i++)
    {
        if (!locs_set[i])
        {
            current = &sys[i];
            num_upstreams = 0;
            while (current && !(locs_set[current->location]))
            {
                upstreams[num_upstreams++] = current->location;
                current = current->child;
            }

            if (current)	//Gauge downstream
            {
                ratio = locs_newq[current->location] / upareas[current->location];
                //ratio = locs_newq[current->location] / current->params->ve[area_idx];
                for (j = 0; j < num_upstreams; j++)
                {
                    locs_set[upstreams[j]] = 1;
                    locs_newq[upstreams[j]] = upareas[upstreams[j]] * ratio;
                    //locs_newq[upstreams[j]] = sys[upstreams[j]]->params->ve[area_idx] * ratio;
                }
            }
        }
    }

    //Set the discharges downstream from each gauge. This is for locations with no downstream gauges.
    unsigned int *counter = (unsigned int*)malloc(N * sizeof(unsigned int));

    //This follows each gauge downstream until a link is found with locs_set or an outlet.
    //counter is set to help track the closest gauge.
    for (i = 0; i < N; i++)
        counter[i] = N;

    for (i = 0; i < num_obs; i++)
    {
        unsigned int loc = obs_locs[i];
        unsigned int prev_loc = loc;
        current = sys[loc].child;
        counter[loc] = 0;

        if (current)
        {
            unsigned int curr_loc = current->location;
            while (!locs_set[curr_loc])
            {
                if (counter[curr_loc] > counter[prev_loc] + 1)
                {
                    counter[curr_loc] = counter[prev_loc] + 1;
                    locs_newq[curr_loc] = upareas[curr_loc] * locs_newq[loc] / upareas[loc];
                    //locs_newq[curr_loc] = current->params->ve[area_idx] * locs_newq[loc] / sys[loc]->params->ve[area_idx];
                }
                prev_loc = curr_loc;
                current = current->child;
                if (current)	curr_loc = current->location;
                else		break;
            }
        }
    }

    for (i = 0; i < N; i++)
        if (counter[i] < N)	locs_set[i] = 1;

    free(counter);

    //Set the determined discharge. If a link's discharge was not determined by the above process, then it lies in a totally ungauged basin.
    for (i = 0; i < N; i++)
    {
        if (locs_set[i])
            x[problem_dim*i] = locs_newq[i];
    }

    //Clean up
    free(locs_set);
    free(upstreams);
    free(locs_newq);
    free(upareas);
    return 0;
}

//Reads a .das file
bool InitAssimData(AssimData* assim, const char* assim_filename)
{
    //unsigned int N = asynch->N, string_size = asynch->globals->string_size;
    //Lookup *id_to_loc = asynch->id_to_loc;
    int errorcode, valsread;
    FILE* inputfile = NULL;
    char end_char;
    char line_buffer[ASYNCH_MAX_LINE_LENGTH];
    const unsigned int line_buffer_len = ASYNCH_MAX_LINE_LENGTH;

    MPI_Barrier(MPI_COMM_WORLD);

    memset(assim, 0, sizeof(AssimData));

    if (my_rank == 0)
    {
        //Open file
        inputfile = fopen(assim_filename, "r");
        errorcode = 0;
        if (!inputfile)
        {
            printf("[%i]: Error opening assimilation file %s.\n", my_rank, assim_filename);
            errorcode = 1;
        }
    }

    //Check if the assimilation file was openned
    MPI_Bcast(&errorcode, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (errorcode)	return false;

    //Read the .model variant
    ReadLineFromTextFile(inputfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%s", assim->model);
    if (ReadLineError(valsread, 1, "model variant"))	return false;

    //Read the .dbc file
    ReadLineFromTextFile(inputfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%s %lf", assim->db_filename, &assim->obs_time_step);
    if (ReadLineError(valsread, 2, "assimilation dbc filename and time resolution of observations"))	return false;
    ReadDBC(assim->db_filename, &assim->conninfo);

    //Read the num of observation (assimilation window)
    ReadLineFromTextFile(inputfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%u", &(assim->num_steps));
    if (ReadLineError(valsread, 1, "number of observations to use"))	return false;

    //Read the max least squares iterations
    ReadLineFromTextFile(inputfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%u", &(assim->max_least_squares_iters));
    if (ReadLineError(valsread, 1, "least squares iterations"))	return false;

    //Read ending mark
    ReadLineFromTextFile(inputfile, line_buffer, line_buffer_len);
    valsread = sscanf(line_buffer, "%c", &end_char);
    if (ReadLineError(valsread, 1, "ending mark"))	return false;

    //Clean up
    if (my_rank == 0)
        fclose(inputfile);
    MPI_Barrier(MPI_COMM_WORLD);
    if (end_char != '#')
    {
        if (my_rank == 0)	printf("Error: Ending mark not seen in %s.\n", assim_filename);
        return false;
    }

    return true;
}

//Trashes an AssimData object
void FreeAssimData(AssimData* assim)
{
    if (assim->id_to_assim)
    {
        free(assim->id_to_assim);
    }
    if (assim->obs_locs)
        free(assim->obs_locs);
}


// Assign the list of observation in AssimData (obs_locs, num_obs)
int GetObservationsIds(const AsynchSolver* asynch, AssimData* assim)
{
    int errorcode = 0;
    unsigned int i, dropped, *gauged_ids = NULL, N = asynch->N;
    Lookup *id_to_loc = asynch->id_to_loc;
    char query[ASYNCH_MAX_QUERY_LENGTH];
    PGresult *res;

    //Get link ids with gauges
    if (my_rank == 0)
    {
        ConnectPGDB(&assim->conninfo);
        sprintf(query, assim->conninfo.queries[0]);
        res = PQexec(assim->conninfo.conn, query);
        if (CheckResError(res, "getting list of gauge link ids"))
            errorcode = 1;
        else
        {
            dropped = 0;
            assim->num_obs = PQntuples(res);
            assim->obs_locs = (unsigned int*)malloc(assim->num_obs * sizeof(unsigned int));
            gauged_ids = (unsigned int*)malloc(assim->num_obs * sizeof(unsigned int));
            for (i = 0; i < assim->num_obs; i++)
            {
                gauged_ids[i - dropped] = atoi(PQgetvalue(res, i, 0));
                assim->obs_locs[i - dropped] = find_link_by_idtoloc(gauged_ids[i], id_to_loc, N);
                if (assim->obs_locs[i - dropped] > N)
                {
                    printf("Warning: Ignoring gauge at link id %u. No link with this id is in the network.\n", gauged_ids[i - dropped]);
                    dropped++;
                }
            }

            //Resize the obs_locs array
            assim->num_obs -= dropped;
            assim->obs_locs = (unsigned int*)realloc(assim->obs_locs, assim->num_obs * sizeof(unsigned int));
        }

        PQclear(res);
        DisconnectPGDB(&assim->conninfo);
    }

    //Give the locations to the other procs
    MPI_Bcast(&errorcode, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!errorcode)
    {
        MPI_Bcast(&(assim->num_obs), 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        if (my_rank != 0)
        {
            assim->obs_locs = (unsigned int*)malloc(assim->num_obs * sizeof(unsigned int));
            gauged_ids = (unsigned int*)malloc(assim->num_obs * sizeof(unsigned int));
        }
        MPI_Bcast(assim->obs_locs, assim->num_obs, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        MPI_Bcast(gauged_ids, assim->num_obs, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        //Build id_to_assim
        assim->id_to_assim = malloc(assim->num_obs * sizeof(Lookup));	//!!!! Should really have an array (should REALLY be a hash table) !!!!
        for (i = 0; i < assim->num_obs; i++)
        {
            assim->id_to_assim[i].id = gauged_ids[i];
            assim->id_to_assim[i].loc = i;
        }
        merge_sort_by_ids(assim->id_to_assim, assim->num_obs);
    }

    if (gauged_ids)	free(gauged_ids);
    return errorcode;
}

//Download gauge data and store them in d. d orders the data by time, then location.
//Data are obtained from times starting at background_time_unix and ending at steps_to_use*inc+background_time_unix.
//Returns 0 if everything went as planned. Returns 1 if some data was not available. Returns -1 if an error occurred.
int GetObservationsData(const AssimData* assim, const Lookup * const id_loc_loc, unsigned int N, unsigned int background_time_unix, double* d)
{
    unsigned int i, n, end_time_unix, *obs_locs = assim->obs_locs, idx, inc_secs = (unsigned int)(assim->obs_time_step*60.0 + 1e-3), unix_t, id, num_obs = assim->num_obs;
    unsigned int num_steps = assim->num_steps, current_time;
    int errorcode = 0;
    PGresult *res;
    const char * const query = assim->conninfo.query;

    //Reset d
    for (i = 0; i < num_steps*num_obs; i++)	d[i] = -1.0;

    //Download the data
    if (my_rank == 0)
    {
        if (ConnectPGDB(&assim->conninfo))
            errorcode = -1;
        else
        {
            end_time_unix = background_time_unix + num_steps * 60 * assim->obs_time_step;
            sprintf(query, assim->conninfo.queries[1], background_time_unix, end_time_unix);	//Assumes data sorted by time, then link id (!!!! No, is it the reverse? !!!!)
            res = PQexec(assim->conninfo.conn, query);
            if (CheckResError(res, "downloading gauge data") == 0)
            {
                //Unpack the data
                n = PQntuples(res);
                for (i = 0; i < n;)
                {
                    id = atoi(PQgetvalue(res, i, 0));
                    /*
                    printf("Unpacking data for id %u...\n",id);
                    int stopper = 1;
                    while(id == 204046 && stopper == 1)
                    {
                    sleep(5);
                    }
                    */
                    idx = find_link_by_idtoloc(id, assim->id_to_assim, num_obs);	//!!!! What if no data for a link is available? Doing nothing keeps the values the same. Or zero... !!!!
                    if (idx < num_obs)
                    {
                        while (i < n && atoi(PQgetvalue(res, i, 0)) == id)
                        {
                            //Check if the first time needed is available
                            unix_t = atoi(PQgetvalue(res, i, 1));
                            if (background_time_unix == unix_t)
                            {
                                //d[idx*steps_to_use] = atof(PQgetvalue(res,i,2));
                                d[idx] = atof(PQgetvalue(res, i, 2));
                                i++;
                            }
                            else	//For now, just take from the next available reading
                            {
                                d[idx] = atof(PQgetvalue(res, i, 2));
                                //printf("!!!! Error: didn't get a gauge value for id = %u at time %u. !!!!\n",id,background_time_unix);
                                //MPI_Abort(MPI_COMM_WORLD,1);
                            }

                            //i++;
                            current_time = 1;
                            unix_t = background_time_unix + inc_secs;
                            //unix_t += inc_secs;
                            //for(;i<n && atoi(PQgetvalue(res,i,0)) == id;i++)
                            while (i < n && atoi(PQgetvalue(res, i, 0)) == id)
                            {
                                if (atoi(PQgetvalue(res, i, 1)) == unix_t)
                                {
                                    //d[idx*steps_to_use+current_time] = atof(PQgetvalue(res,i,2));
                                    d[idx + current_time*num_obs] = atof(PQgetvalue(res, i, 2));
                                    current_time++;
                                    unix_t += inc_secs;
                                    i++;
                                }
                                else if (atoi(PQgetvalue(res, i, 1)) > unix_t)	//Not available!
                                {
                                    //Just use the last observation
                                    //d[idx*steps_to_use+current_time] = atof(PQgetvalue(res,i-1,2));
                                    //d[idx+current_time*num_obs] = atof(PQgetvalue(res,i-1,2));
                                    d[idx + current_time*num_obs] = d[idx + (current_time - 1)*num_obs];
                                    current_time++;
                                    unix_t += inc_secs;
                                    errorcode = 1;
                                }
                                else	//The time step for this location does not match time_inc
                                {
                                    i++;	//This does NOT use the closest gauge reading. It basically ignores all off-timestep observations
                                }
                            }
                        }
                    }
                    else	//ID not a gauged location. So skip it.
                    {
                        i++;
                        while (i < n && atoi(PQgetvalue(res, i, 0)) == id)	i++;
                    }
                }
            }
            else
                errorcode = -1;

            //Clean up
            PQclear(res);
            DisconnectPGDB(&assim->conninfo);
        }
    }

    //Send the data to everyone
    MPI_Bcast(&errorcode, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (errorcode > -1)
        MPI_Bcast(d, num_steps*num_obs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (my_rank == 0)
    {
        printf("At time %u, got:\n", background_time_unix);
        for (i = 0; i < num_steps*num_obs; i++)
            printf("%f ", d[i]);
        printf("\n+++++++++++++++++++++++\n");
    }
    return errorcode;
}


//Checks which dicharge values are having problems converging. Everything upstreams is cut in half.
//!!!! This is a little sloppy with overlapping domains of influence. Could reduce by more than 1/2. Also broadcasting static information. !!!!
bool ReduceBadDischargeValues(Link* sys, int* assignments, unsigned int N, double* d_full, double* q, unsigned int num_steps, unsigned int* obs_locs, unsigned int num_obs, double* x_start, unsigned int assim_dim, double limit)
{
    bool res = false;

    unsigned int i, j, k, num_links;
    Link* current;
    double new_diff, old_diff, factor = 0.8;
    UpstreamData* updata;

    for (i = 0; i < num_obs; i++)
    {
        if (d_full[i] < 0.0)
        {
            if (my_rank == 0)
                printf("!!!! Skipping scaling of link %u. Observation is %f. !!!!\n", obs_locs[i], d_full[i]);
            continue;
        }

        current = &sys[obs_locs[i]];

        //Check how far off each gauge is
        //old_diff = d_full[i] - q[i];
        old_diff = fabs(q[i] - d_full[i]);
        for (j = 1; j < num_steps; j++)
        {
            //new_diff = d_full[i+j*num_obs] - q[i+j*num_obs];
            new_diff = fabs(q[i + j*num_obs] - d_full[i + j*num_obs]);
            if (old_diff < new_diff)	//Getting worse...
                old_diff = new_diff;
            else
                break;
        }

        //if(j == steps_to_use && old_diff > limit)	//Got progressively worse, and ended badly. Make a change!
        if (j == num_steps)
        {
            res = true;

            //Calculate the scaling factor
            factor = d_full[i + (num_steps - 1)*num_obs] / q[i + (num_steps - 1)*num_obs];

            //Reduce q at the gauge
            x_start[obs_locs[i] * assim_dim] *= factor;

            if (assignments[obs_locs[i]] == my_rank)
            {
                printf("[%i] !!!! Reducing discharges above link %u. Factor is %f. !!!!\n", my_rank, current->ID, factor);
                updata = (UpstreamData*)current->user;

                //Sum up the total number of links upstreams from current
                num_links = updata->num_upstreams;
                //for (j = 0; j < current->num_parents; j++)
                //    num_links += updata->num_upstreams[j];
                MPI_Bcast(&num_links, 1, MPI_UNSIGNED, my_rank, MPI_COMM_WORLD);

                //for (j = 0; j < current->num_parents; j++)
                //{
                //    unsigned int num_upstreams = ((UpstreamData *)updata->upstreams[i]->user)->num_upstreams;
                //    Link **upstreams = ((UpstreamData *)current->user)->upstreams;

                //    for (k = 0; k < num_upstreams; k++)
                //    {
                //        MPI_Bcast(&(updata->upstreams[j][k]), 1, MPI_UNSIGNED, my_rank, MPI_COMM_WORLD);
                //        x_start[updata->upstreams[j][k] * assim_dim] *= factor;
                //    }
                //}

                unsigned int *upstream_location = malloc(num_links * sizeof(unsigned int));

                for (k = 0; k < num_links; k++)
                {
                    upstream_location[k] = updata->upstreams[k]->location;
                    x_start[upstream_location[k] * assim_dim] *= factor;
                }

                MPI_Bcast(upstream_location, num_links, MPI_UNSIGNED, my_rank, MPI_COMM_WORLD);

                free(upstream_location);
            }
            else
            {
                MPI_Bcast(&num_links, 1, MPI_UNSIGNED, assignments[obs_locs[i]], MPI_COMM_WORLD);

                unsigned int *upstream_location = malloc(num_links * sizeof(unsigned int));

                //for (j = 0; j < num_links; j++)
                //{
                //    MPI_Bcast(&loc, 1, MPI_UNSIGNED, assignments[obs_locs[i]], MPI_COMM_WORLD);
                //    x_start[loc * assim_dim] *= factor;
                //}

                MPI_Bcast(&upstream_location, num_links, MPI_UNSIGNED, assignments[obs_locs[i]], MPI_COMM_WORLD);

                for (k = 0; k < num_links; k++)
                    x_start[upstream_location[k] * assim_dim] *= factor;
            }
        }
    }

    return res;
}


//This creates a snapshot of the model states only (i.e. no variational equation states are saved)
int SnapShot_ModelStates(AsynchSolver* asynch, unsigned int problem_dim)
{
    Link *sys = asynch->sys;
    GlobalVars *globals = asynch->globals;
    int *assignments = asynch->assignments;
    unsigned int i, j, N = asynch->N;
    FILE* output;
    double buffer[ASYNCH_MAX_DIM];

    if (my_rank == 0)	//Creating the file
    {
        output = fopen(globals->dump_loc_filename, "w");
        if (output == NULL)
        {
            printf("[%i]: Error opening file %s.\n", my_rank, globals->dump_loc_filename);
            i = 1;
        }
        else	i = 0;
        MPI_Bcast(&i, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        if (i)	return 1;

        fprintf(output, "%hu\n%u\n0.0\n\n", globals->model_uid, N);

        for (i = 0; i < N; i++)
        {
            if (assignments[i] != 0)
                MPI_Recv(buffer, problem_dim, MPI_DOUBLE, assignments[i], i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            else
                for (j = 0; j < problem_dim; j++)	buffer[j] = sys[i].my->list.tail->y_approx[j];

            fprintf(output, "%u\n", sys[i].ID);
            for (j = 0; j < problem_dim; j++)	fprintf(output, "%.6e ", buffer[j]);
            fprintf(output, "\n");
        }

        fclose(output);
    }
    else			//Sending data to proc 0
    {
        //Check for error
        MPI_Bcast(&i, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        if (i)	return 1;

        for (i = 0; i < N; i++)
        {
            if (assignments[i] == my_rank)
            {
                for (j = 0; j < problem_dim; j++)	buffer[j] = sys[i].my->list.tail->y_approx[j];
                MPI_Send(buffer, problem_dim, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}

//int GaugeDataAvailable(AssimData* assim, unsigned int start_time, unsigned int end_time)
//{
//    int isnull = 1, num_obs, numvalues, expect;
//    PGresult* res;
//    char* query = NULL;
//
//    if (my_rank == 0)
//    {
//        query = (char*)malloc(1024 * sizeof(char));
//        ConnectPGDB(&assim->conninfo);
//        num_obs = assim->num_obs;
//
//        //Get number of observations available
//        sprintf(query, assim->conninfo->queries[2], start_time, end_time);
//        res = PQexec(assim->conninfo->conn, query);
//        if (CheckResError(res, "checking for new observation data"))	goto error;
//        numvalues = atoi(PQgetvalue(res, 0, 0));
//        PQclear(res);
//
//        //Is a good number of observations available? Waiting for 80%.
//        expect = (assim->num_steps * num_obs * 8) / 10;
//        if (numvalues >= expect)	isnull = 0;
//
//        printf("%i/%i new observations available. Expecting at least %i.\n", numvalues, assim->num_steps * num_obs, expect);
//
//    error:
//        DisconnectPGDB(assim->conninfo);
//    }
//
//    MPI_Bcast(&isnull, 1, MPI_INT, 0, MPI_COMM_WORLD);
//    if (query)	free(query);
//    return isnull;
//}

