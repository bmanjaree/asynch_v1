#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include <mpi.h>

#if defined(HAVE_METIS)
#include <metis.h>
#endif

#include <partition.h>

//Partitions a river system by first partitioning the leaves.
//Link** sys: The river system.
//int N: The number of links in the river system.
//Link** leaves: The list of leaves to the system.
//int numleaves: Number of leaves in the system.
//int** my_sys (set by this method): An array that will contain the location in the system of each link assigned to this process.
//int* my_N (set by this method): Will be the number of links assigned to this process.
//TransData* my_data (this is assumed to be allocated already, but contents will be set here): Information about the links that
//				each process will communicate information about will be stored in my_data.
//short int *getting (set by this method, but assumed space is allocated with N entries): getting[i] will equal 1 if this process needs
//				to receive information about link i, 0 otherwise.
//Returns an array of N integers that take values from 0 to np-1. The i-th entry contains which process the link in location i of
//				the system is assigned to.
int* Partition_System_By_Leaves(Link *sys, unsigned int N, Link **leaves, unsigned int numleaves, Link ***my_sys, unsigned int *my_N, TransData *my_data, short int *getting)
{
    //Assign leaves
    unsigned int i, start_index, end_index, extras;
    int j;
    unsigned int nodes_per_proc = numleaves / np;	//Number of leaves assigned to each process (except the last)

    start_index = nodes_per_proc * my_rank;
    if (my_rank == np - 1)
        end_index = numleaves;
    else
        end_index = nodes_per_proc * (my_rank + 1);
    *my_N = end_index - start_index;
    
    unsigned int my_max_nodes = N - numleaves + nodes_per_proc;
    
    //Pointer to this processes links (in sys)
    *my_sys = malloc(my_max_nodes * sizeof(Link *));
    memset(*my_sys, 0, my_max_nodes * sizeof(Link *));
    
    for (i = start_index; i < end_index; i++)
        (*my_sys)[i - start_index] = leaves[i];

    //Initialize assignments
    int* assignments = malloc(N * sizeof(int));
    for (i = 0; i < N; i++)
        assignments[i] = -1;

    //Calculate and store the assignments for the leaves
    extras = numleaves % np;	//This is how many extra leaves we have (they are assigned to proc np above). We assign them to np-1.
    end_index = numleaves - extras;
    for (i = 0; i < end_index; i++)
        assignments[leaves[i]->location] = i / nodes_per_proc;
    for (i = end_index; i < numleaves; i++)
        assignments[leaves[i]->location] = np - 1;

    //Initialize getting
    for (i = 0; i < N; i++)	getting[i] = 0;

    //Assign the rest of the links
    //I do this by starting at each leaf, and assigning the leaf's child to its left parent's (0) process.
    Link *current, *q;
    unsigned int curr_idx;
    for (i = 0; i < numleaves; i++)
    {
        current = leaves[i]->child;
        if (current != NULL)
        {
            curr_idx = current->location;
            q = leaves[i];
            while (assignments[curr_idx] == -1)
            {
                assignments[curr_idx] = assignments[current->parents[0]->location];
                if (assignments[curr_idx] == my_rank)	//If this node is assigned to this process
                {
                    (*my_sys)[*my_N] = &sys[curr_idx];
                    (*my_N)++;
                }
                q = current;
                current = current->child;
                if (current != NULL)	curr_idx = current->location;
                else			break;
            }

            if (current != NULL)
            {
                if (assignments[curr_idx] != assignments[q->location])
                {
                    int q_proc = assignments[q->location];
                    int curr_proc = assignments[curr_idx];

                    //Check if this process will receive data
                    if (my_rank == curr_proc)
                    {
                        //my_data->receive_data[q_proc][my_data->receive_size[q_proc]] = q;
                        (my_data->receive_size[q_proc])++;
                        getting[q->location] = 1;
                    }

                    //Check if this process will send data
                    if (my_rank == q_proc)
                    {
                        //my_data->send_data[curr_proc][my_data->send_size[curr_proc]] = q;
                        (my_data->send_size[curr_proc])++;
                    }
                }
            }
        }
    }

    //Reorder my_sys so that the links with lower numbering are towards the beginning
    merge_sort_by_distance(*my_sys, *my_N);

    //Allocate space in my_data for recieving and sending
    for (j = 0; j < np; j++)
    {
        my_data->receive_data[j] = (Link**)malloc(my_data->receive_size[j] * sizeof(Link*));
        my_data->send_data[j] = (Link**)malloc(my_data->send_size[j] * sizeof(Link*));
    }
    int* current_receive_size = (int*)calloc(np, sizeof(int));
    int* current_send_size = (int*)calloc(np, sizeof(int));
    for (i = 0; i < N; i++)	assignments[i] = -1;

    //Recalculate and store the assignments for the leaves
    extras = numleaves % np;
    end_index = numleaves - extras;
    for (i = 0; i < end_index; i++)
        assignments[leaves[i]->location] = i / nodes_per_proc;
    for (i = end_index; i < numleaves; i++)
        assignments[leaves[i]->location] = np - 1;

    //Assign links to receive_data and send_data
    for (i = 0; i < numleaves; i++)
    {
        current = leaves[i]->child;
        if (current != NULL)
        {
            curr_idx = current->location;
            q = leaves[i];

            while (assignments[curr_idx] == -1)
            {
                assignments[curr_idx] = assignments[current->parents[0]->location];
                q = current;
                current = current->child;
                if (current != NULL)	curr_idx = current->location;
                else			break;
            }

            if (current != NULL)
            {
                if (assignments[curr_idx] != assignments[q->location])
                {
                    int q_proc = assignments[q->location];
                    int curr_proc = assignments[curr_idx];

                    //Check if this process will receive data
                    if (my_rank == curr_proc)
                    {
                        my_data->receive_data[q_proc][current_receive_size[q_proc]] = q;
                        current_receive_size[q_proc]++;
                    }

                    //Check if this process will send data
                    if (my_rank == q_proc)
                    {
                        my_data->send_data[curr_proc][current_send_size[curr_proc]] = q;
                        current_send_size[curr_proc]++;
                    }
                }
            }
        }
    }

    //Clean up
    free(current_receive_size);
    free(current_send_size);

    //Resize down
    *my_sys = realloc(*my_sys, *my_N * sizeof(Link *));

    return assignments;
}


//Partitions a river system by first partitioning the leaves. Makes adjustments so all leaves are assigned to their child's process.
//Link** sys: The river system.
//int N: The number of links in the river system.
//Link** leaves: The list of leaves to the system.
//int numleaves: Number of leaves in the system.
//int** my_sys (set by this method): An array that will contain the location in the system of each link assigned to this process.
//int* my_N (set by this method): Will be the number of links assigned to this process.
//TransData* my_data (this is assumed to be allocated already, but contents will be set here): Information about the links that
//				each process will communicate information about will be stored in my_data.
//short int *getting (set by this method, but assumed space is allocated with N entries): getting[i] will equal 1 if this process needs
//				to receive information about link i, 0 otherwise.
//Returns an array of N integers that take values from 0 to np-1. The i-th entry contains which process the link in location i of
//				the system is assigned to.
int* Partition_System_By_Leaves_2(Link *sys, unsigned int N, Link **leaves, unsigned int numleaves, Link ***my_sys, unsigned int * my_N, TransData *my_data, short int *getting)
{
    //Assign leaves
    unsigned int i, j, k, l, start_index, end_index, extras;
    int ii;
    unsigned int nodes_per_proc = numleaves / np;	//Number of leaves assigned to each process (except the last)

    start_index = nodes_per_proc * my_rank;
    if (my_rank == np - 1)		
        end_index = numleaves;
    else				
        end_index = nodes_per_proc * (my_rank + 1);
    *my_N = end_index - start_index;

    unsigned my_max_nodes = N - numleaves + nodes_per_proc;

    //Pointer to this processes links (in sys)
    *my_sys = malloc(my_max_nodes * sizeof(Link *));
    memset(*my_sys, 0, my_max_nodes * sizeof(Link *));

    for (i = 0; i < my_max_nodes; i++)
        (*my_sys)[i] = NULL;
    for (i = start_index; i < end_index; i++)
        (*my_sys)[i - start_index] = leaves[i];

    //Initialize assignments
    int* assignments = (int*)malloc(N * sizeof(int));
    for (i = 0; i < N; i++)	assignments[i] = -1;

    //Calculate and store the assignments for the leaves
    extras = numleaves % np;	//This is how many extra leaves we have (they are assigned to proc np above). We assign them to np-1.
    end_index = numleaves - extras;
    for (i = 0; i < end_index; i++)
        assignments[leaves[i]->location] = i / nodes_per_proc;
    for (i = end_index; i < numleaves; i++)
        assignments[leaves[i]->location] = np - 1;

    //Initialize getting
    for (i = 0; i < N; i++)	getting[i] = 0;

    //Assign the rest of the links
    //I do this by starting at each leaf, and assigning the leaf's child to its left parent's (0) process.
    Link *current, *q;
    int curr_idx;
    for (i = 0; i < numleaves; i++)
    {
        current = leaves[i]->child;
        if (current != NULL)
        {
            curr_idx = current->location;
            q = leaves[i];

            while (assignments[curr_idx] == -1)
            {
                assignments[curr_idx] = assignments[current->parents[0]->location];
                if (assignments[curr_idx] == my_rank)	//If this node is assigned to this process
                {
                    (*my_sys)[*my_N] = &sys[curr_idx];
                    (*my_N)++;
                }

                //Check if any parents are leaves. If so, assign them to this process
                for (j = 0; j < current->num_parents; j++)
                {
                    if (current->parents[j]->num_parents == 0 && assignments[current->parents[j]->location] != assignments[curr_idx])
                    {
                        //Remove the parent from the my_sys it is currently in
                        if (assignments[current->parents[j]->location] == my_rank)
                        {
                            for (k = 0; k < *my_N; k++)
                            {
                                if ((*my_sys)[k] == current->parents[j])
                                {
                                    for (l = k; l < *my_N - 1; l++)
                                        (*my_sys)[l] = (*my_sys)[l + 1];
                                    break;
                                }
                            }
                            (*my_N)--;
                        }

                        //Assign the parent
                        assignments[current->parents[j]->location] = assignments[current->location];
                        if (assignments[current->parents[j]->location] == my_rank)
                        {
                            (*my_sys)[*my_N] = current->parents[j];
                            (*my_N)++;
                        }
                    }
                }

                q = current;
                current = current->child;
                if (current != NULL)	curr_idx = current->location;
                else			break;
            }

            if (current != NULL)
            {
                if (assignments[curr_idx] != assignments[q->location])
                {
                    int q_proc = assignments[q->location];
                    int curr_proc = assignments[curr_idx];

                    //Check if this process will receive data
                    if (my_rank == curr_proc)
                    {
                        (my_data->receive_size[q_proc])++;
                        getting[q->location] = 1;
                    }

                    //Check if this process will send data
                    if (my_rank == q_proc)	(my_data->send_size[curr_proc])++;
                }
            }
        }
    }

    //Reorder my_sys so that the links with lower numbering are towards the beginning
    merge_sort_by_distance(*my_sys, *my_N);

    //Allocate space in my_data for recieving and sending
    for (ii = 0; ii < np; ii++)
    {
        my_data->receive_data[ii] = (Link**)malloc(my_data->receive_size[ii] * sizeof(Link*));
        my_data->send_data[ii] = (Link**)malloc(my_data->send_size[ii] * sizeof(Link*));
    }
    int* current_receive_size = (int*)calloc(np, sizeof(int));
    int* current_send_size = (int*)calloc(np, sizeof(int));
    for (i = 0; i < N; i++)	assignments[i] = -1;

    //Calculate and store the assignments for the leaves
    extras = numleaves % np;	//This is how many extra leaves we have (they are assigned to proc np above). We assign them to np-1.
    end_index = numleaves - extras;
    for (i = 0; i < end_index; i++)
        assignments[leaves[i]->location] = i / nodes_per_proc;
    for (i = end_index; i < numleaves; i++)
        assignments[leaves[i]->location] = np - 1;

    for (i = 0; i < numleaves; i++)
    {
        current = leaves[i]->child;
        if (current != NULL)
        {
            curr_idx = current->location;
            q = leaves[i];

            while (assignments[curr_idx] == -1)
            {
                assignments[curr_idx] = assignments[current->parents[0]->location];

                //Check if any parents are leaves. If so, assign them to this process
                for (j = 0; j < current->num_parents; j++)
                {
                    if (current->parents[j]->num_parents == 0 && assignments[current->parents[j]->location] != assignments[curr_idx])
                    {
                        //Assign the parent
                        assignments[current->parents[j]->location] = assignments[current->location];
                    }
                }

                q = current;
                current = current->child;
                if (current != NULL)	curr_idx = current->location;
                else			break;
            }

            if (current != NULL)
            {
                if (assignments[curr_idx] != assignments[q->location])
                {
                    int q_proc = assignments[q->location];
                    int curr_proc = assignments[curr_idx];

                    //Check if this process will receive data
                    if (my_rank == curr_proc)
                    {
                        my_data->receive_data[q_proc][current_receive_size[q_proc]] = q;
                        current_receive_size[q_proc]++;
                    }

                    //Check if this process will send data
                    if (my_rank == q_proc)
                    {
                        my_data->send_data[curr_proc][current_send_size[curr_proc]] = q;
                        current_send_size[curr_proc]++;
                    }
                }
            }
        }
    }

    //Clean up
    free(current_receive_size);
    free(current_send_size);

    //Resize down
    *my_sys = realloc(*my_sys, *my_N * sizeof(Link *));

    return assignments;
}


#if defined(HAVE_METIS)

int* Partition_METIS_ByEqs(Link* sys, unsigned int N, Link** leaves, unsigned int numleaves, Link** my_sys, unsigned int* my_N, TransData* my_data, short int *getting)
{
    unsigned int loc, retval;
    unsigned int nodes_per_proc = numleaves / np;	//Number of leaves assigned to each process (except the last)
    Link* current;

    unsigned int start_index = nodes_per_proc * my_rank;
    unsigned int end_index;
    if (my_rank == np - 1)		
        end_index = numleaves;
    else
        end_index = nodes_per_proc * (my_rank + 1);
    //	*my_N = end_index - start_index;
    *my_N = 0;
    unsigned int my_max_nodes = N - numleaves + nodes_per_proc;
    *my_sys = malloc(my_max_nodes * sizeof(unsigned int));	//The indices of this processes links (in sys)
    for (unsigned int i = 0; i < my_max_nodes; i++)
        my_sys[i] = NULL;
    for (unsigned int i = start_index; i < end_index; i++)
        my_sys[i - start_index] = leaves[i];
    for (unsigned int i = 0; i < N; i++)
        getting[i] = 0;

    //Initialize assignments
    int* assignments = malloc(N * sizeof(int));
    for (unsigned int i = 0; i < N; i++)
        assignments[i] = -1;

    //Form the graph to partition
    idx_t* xadj = malloc((N + 1) * sizeof(idx_t));
    idx_t* adjncy = malloc(2 * (N - 1) * sizeof(idx_t));
    idx_t index = 0;

    for (unsigned int i = 0; i < N; i++)
    {
        xadj[i] = index;
        current = &sys[i];
        if (current->child != NULL)
        {
            adjncy[index] = current->child->location;
            index++;
        }
        for (unsigned int j = 0; j < current->num_parents; j++)
        {
            adjncy[index] = current->parents[j]->location;
            index++;
        }
    }
    xadj[N] = 2 * (N - 1);

    //Partition the system
    idx_t nverts = N;
    idx_t parts = np;
    idx_t ncon = 1;
    idx_t objval;
    idx_t* partitions = calloc(N, sizeof(idx_t));
    idx_t* vwgt = malloc(N * sizeof(idx_t));

    for (unsigned int i = 0; i < N; i++)
        vwgt[i] = sys[i].dim;
    if (np != 1)
    {
        retval = METIS_PartGraphKway(&nverts, &ncon, xadj, adjncy, vwgt, vwgt, NULL, &parts, NULL, NULL, NULL, &objval, partitions);
        //retval = METIS_PartGraphKway(&nverts,&ncon,xadj,adjncy,NULL,NULL,NULL,&parts,NULL,NULL,NULL,&objval,partitions);
        if (retval != METIS_OK)
        {
            printf("Error: METIS returned error code %i.\n", retval);
            return NULL;
        }
    }

    *my_N = 0;
    for (unsigned int i = 0; i < N; i++)
    {
        assignments[i] = partitions[i];	//!!!! Just use assignments? !!!!
        if (partitions[i] == my_rank)
        {
            my_sys[*my_N] = &sys[i];
            (*my_N)++;
        }
    }

    //Set the getting array and determine number of sending and receiving links
    for (unsigned int i = 0; i < *my_N; i++)
    {
        //Receiving
        for (unsigned int j = 0; j < my_sys[i]->num_parents; j++)
        {
            loc = my_sys[i]->parents[j]->location;
            if (assignments[loc] != my_rank)
            {
                getting[loc] = 1;
                my_data->receive_size[assignments[loc]]++;
            }
        }

        //Sending
        if (my_sys[i]->child != NULL)
        {
            loc = my_sys[i]->child->location;
            if (assignments[loc] != my_rank)
                my_data->send_size[assignments[loc]]++;
        }
    }

    //Reorder my_sys so that the links with lower numbering are towards the beginning
    merge_sort_by_distance(my_sys, *my_N);

    //Allocate space in my_data for recieving and sending
    for (int j = 0; j < np; j++)
    {
        my_data->receive_data[j] = malloc(my_data->receive_size[j] * sizeof(Link*));
        my_data->send_data[j] = malloc(my_data->send_size[j] * sizeof(Link*));
    }

    //Set the receive_data and send_data arrays
    int* current_receive_size = (int*)calloc(np, sizeof(int));
    int* current_send_size = (int*)calloc(np, sizeof(int));
    for (unsigned int i = 0; i < *my_N; i++)
    {
        //Receiving
        for (unsigned int j = 0; j < my_sys[i]->num_parents; j++)
        {
            loc = my_sys[i]->parents[j]->location;
            if (assignments[loc] != my_rank)
            {
                my_data->receive_data[assignments[loc]][current_receive_size[assignments[loc]]] = &sys[loc];
                current_receive_size[assignments[loc]]++;
            }
        }

        //Sending
        if (my_sys[i]->child != NULL)
        {
            loc = my_sys[i]->child->location;
            if (assignments[loc] != my_rank)
            {
                my_data->send_data[assignments[loc]][current_send_size[assignments[loc]]] = my_sys[i];
                current_send_size[assignments[loc]]++;
            }
        }
    }

    //Clean up
    free(current_receive_size);
    free(current_send_size);
    free(xadj);
    free(adjncy);
    free(partitions);
    free(vwgt);
    *my_sys = realloc(*my_sys, *my_N * sizeof(Link *));

    return assignments;
}


int* Partition_METIS_Traditional(Link* sys, unsigned int N, Link** leaves, unsigned int numleaves, Link** my_sys, unsigned int* my_N, TransData* my_data, short int *getting)
{
    unsigned int loc, retval;
    unsigned int nodes_per_proc = numleaves / np;	//Number of leaves assigned to each process (except the last)
    Link* current;

    unsigned int start_index = nodes_per_proc * my_rank;
    unsigned int end_index;
    if (my_rank == np - 1)	
        end_index = numleaves;
    else
        end_index = nodes_per_proc * (my_rank + 1);
    //	*my_N = end_index - start_index;
    *my_N = 0;
    unsigned int my_max_nodes = N - numleaves + nodes_per_proc;
    for (unsigned int i = 0; i < my_max_nodes; i++)
        my_sys[i] = NULL;
    for (unsigned int i = start_index; i < end_index; i++)
        my_sys[i - start_index] = leaves[i];
    for (unsigned int i = 0; i < N; i++)
        getting[i] = 0;

    //Initialize assignments
    int* assignments = malloc(N * sizeof(int));
    for (unsigned int i = 0; i < N; i++)
        assignments[i] = -1;

    //Form the graph to partition
    idx_t* xadj = malloc((N + 1) * sizeof(idx_t));
    idx_t* adjncy = malloc(2 * (N - 1) * sizeof(idx_t));
    idx_t index = 0;

    for (unsigned int i = 0; i < N; i++)
    {
        xadj[i] = index;
        current = &sys[i];
        if (current->child != NULL)
        {
            adjncy[index] = current->child->location;
            index++;
        }
        for (unsigned int j = 0; j < current->num_parents; j++)
        {
            adjncy[index] = current->parents[j]->location;
            index++;
        }
    }
    xadj[N] = 2 * (N - 1);

    //Partition the system
    idx_t nverts = N;
    idx_t parts = np;
    idx_t ncon = 1;
    idx_t objval;
    idx_t* partitions = calloc(N, sizeof(idx_t));
    if (np != 1)
    {
        retval = METIS_PartGraphKway(&nverts, &ncon, xadj, adjncy, NULL, NULL, NULL, &parts, NULL, NULL, NULL, &objval, partitions);
        if (retval != METIS_OK)
        {
            printf("Error: METIS returned error code %i.\n", retval);
            return NULL;
        }
    }

    *my_N = 0;
    for (unsigned int i = 0; i < N; i++)
    {
        assignments[i] = partitions[i];	//!!!! Just use assignments? !!!!
        if (partitions[i] == my_rank)
        {
            my_sys[*my_N] = &sys[i];
            (*my_N)++;
        }
    }

    //Set the getting array and determine number of sending and receiving links
    for (unsigned int i = 0; i < *my_N; i++)
    {
        //Receiving
        for (unsigned int j = 0; j < my_sys[i]->num_parents; j++)
        {
            loc = my_sys[i]->parents[j]->location;
            if (assignments[loc] != my_rank)
            {
                getting[loc] = 1;
                my_data->receive_size[assignments[loc]]++;
            }
        }

        //Sending
        if (my_sys[i]->child != NULL)
        {
            loc = my_sys[i]->child->location;
            if (assignments[loc] != my_rank)
                my_data->send_size[assignments[loc]]++;
        }
    }

    //Reorder my_sys so that the links with lower numbering are towards the beginning
    merge_sort_by_distance(my_sys, *my_N);

    //Allocate space in my_data for recieving and sending
    for (int j = 0; j < np; j++)
    {
        my_data->receive_data[j] = (Link**)malloc(my_data->receive_size[j] * sizeof(Link*));
        my_data->send_data[j] = (Link**)malloc(my_data->send_size[j] * sizeof(Link*));
    }

    //Set the receive_data and send_data arrays
    int* current_receive_size = (int*)calloc(np, sizeof(int));
    int* current_send_size = (int*)calloc(np, sizeof(int));
    for (unsigned int i = 0; i < *my_N; i++)
    {
        //Receiving
        for (unsigned int j = 0; j < my_sys[i]->num_parents; j++)
        {
            loc = my_sys[i]->parents[j]->location;
            if (assignments[loc] != my_rank)
            {
                my_data->receive_data[assignments[loc]][current_receive_size[assignments[loc]]] = &sys[loc];
                current_receive_size[assignments[loc]]++;
            }
        }

        //Sending
        if (my_sys[i]->child != NULL)
        {
            loc = my_sys[i]->child->location;
            if (assignments[loc] != my_rank)
            {
                my_data->send_data[assignments[loc]][current_send_size[assignments[loc]]] = my_sys[i];
                current_send_size[assignments[loc]]++;
            }
        }
    }

    //Clean up
    free(current_receive_size);
    free(current_send_size);
    free(xadj);
    free(adjncy);
    free(partitions);
    *my_sys = realloc(my_sys, *my_N * sizeof(Link *));

    return assignments;
}


#endif //defined(HAVE_METIS)
