#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>

#include <misc.h>

static double sq(double val)
{
    return val * val;
}

//Calculates the order of each link in sys
//Assumes order and complete have N spaces already reserved.
//This could be fairly easily incorporated into riversys.c to be more efficient, if needed.
void CalcHortonOrder(Link** sys,unsigned int N,unsigned int* order,unsigned short int* complete)
{
	unsigned int i,j,loc,parentsval;
	Link *current,*root = NULL;

	Link** stack = malloc(N * sizeof(Link*));
	int stack_size = 0;
	Link** leaves = malloc(N * sizeof(Link*));
	unsigned int leaves_size = 0;
	unsigned short int num_parents;

	for(i=0;i<N;i++)
	{
		order[i] = 0;
		complete[i] = 0;
	}

	//Find the root
	for(i=0;i<N;i++)
	{
		if(sys[i]->child == NULL)
		{
			root = sys[i];
			break;
		}
	}

	//Find the leaves
	//Note: this assumes only one root. If there are multiple roots, this needs a for loop.
	stack[0] = root;
	stack_size = 1;
	while(stack_size > 0)
	{
		current = stack[stack_size-1];	//Top of stack
		num_parents = current->num_parents;

		if(num_parents == 0)
		{
			stack_size--;
			leaves[leaves_size] = current;
			leaves_size++;
		}
		else
		{
			//If current is not a leaf, replace it with it's parents
			for(i=0;i<num_parents;i++)
			{
				stack[stack_size - 1 + i] = current->parents[num_parents - 1 - i];
				stack[stack_size - 1 + i]->child = current;
			}
			stack_size += num_parents - 1;
		}
	}

	//Calculate order
	for(i=0;i<leaves_size;i++)
	{
		current = leaves[i];
		order[current->location] = 1;
		//complete[current->location] = 1;	//Exclude order 1 links

		while(current->child != NULL)
		{
			current = current->child;
			loc = current->location;
			num_parents = current->num_parents;
			for(j=0;j<num_parents;j++)
				order[loc] = (order[loc] > order[current->parents[j]->location]) ? order[loc] : order[current->parents[j]->location];

			//Check if current is a complete ordered link
			parentsval = 0;
			for(j=0;j<num_parents;j++)
				if(order[loc] == order[current->parents[j]->location])	parentsval++;
			if(parentsval == num_parents)
			{
				order[loc]++;
				complete[loc] = 1;
			}
		}
	}

	//for(i=0;i<N;i++)
	//	printf("id = %u  order = %u  complete = %hu\n",sys[i]->ID,order[i],complete[i]);

	free(stack);
	free(leaves);
}

//Generate a .str file for use with complete ordered links
void CreateStrComplete(Link** sys,unsigned int N)
{
	unsigned int i,j,changes,rain_order;
	FILE* output = fopen("Cedar30Complete9.str","w");
	if(!output)
	{
		printf("Error creating file in CreateStrComplete.\n");
		return;
	}
	fprintf(output,"%u\n\n",N);

	//Get the order of each link
	printf("Calculating orders...\n");
	unsigned int* order = malloc(N*sizeof(unsigned int));
	unsigned short int* complete = malloc(N*sizeof(unsigned short int));
	CalcHortonOrder(sys,N,order,complete);

	//Create .str file
	changes = 100;
	rain_order = 9;
	printf("Creating .str file...\n");

	for(i=0;i<N;i++)
	{
		if(complete[i] == 1 && order[i] == rain_order)
		{
			fprintf(output,"%u\n%u\n",sys[i]->ID,changes);
			for(j=0;j<changes-1;j++)
				fprintf(output,"%f %f\n",20.0*j,5.0*j);
			for(j=0;j<changes-1;j++)
				fprintf(output,"%f %f\n",20.0*(changes+j),5.0*j);
			for(j=0;j<changes-1;j++)
				fprintf(output,"%f %f\n",20.0*(2*changes+j),5.0*j);
			fprintf(output,"%f %f\n",20.0*3*changes,0.0);
		}
		else
		{
			fprintf(output,"%u\n%u\n",sys[i]->ID,1);
			fprintf(output,"%f %f\n",0.0,0.0);
		}
		fprintf(output,"\n");
	}

	free(order);
	free(complete);
	fclose(output);
}

//Generates a graph for use by the standalone programs from METIS
void CreateGraph(Link** sys,unsigned int N)
{
	unsigned int i,j;
	int offset = sys[0]->ID + 2;
	FILE* output = fopen("Cedar30.gra","w");
	if(!output)
	{
		printf("Error opening file for graph.\n");
		return;
	}

	fprintf(output,"%u %u\n",N,N-1);

	for(i=0;i<N;i++)
	{
		if(sys[i]->child != NULL)	fprintf(output,"%u ",sys[i]->child->location + offset);
		for(j=0;j<sys[i]->num_parents;j++)
			fprintf(output,"%u ",sys[i]->parents[j]->location + offset);
		fprintf(output,"\n");
	}

	fclose(output);
}


//Generates a graph for use by the standalone programs from METIS but takes rainfall into account
/*
void CreateGraphRain(Link** sys,unsigned int N,UnivVars* GlobalVars)
{
	unsigned int i,j,holder;
	char filename[256];
	int offset = -sys[0]->ID + 2;
	unsigned int* total = calloc(N,sizeof(unsigned int));
	float* last = malloc(N*sizeof(float));
	float rainfall_buffer;
	FILE* rainfile;
	FILE* output = fopen("Cedar30Complete1.gra","w");
	if(!output)
	{
		printf("Error opening file for graph.\n");
		return;
	}

	printf("Calculating the number of rainfall changes...\n");

	if(GlobalVars->rain_flag == 2)
	{
		//Calculate the number of rainfall changes that occur
		unsigned int totalfiles = GlobalVars->last_file - GlobalVars->first_file + 1;
		for(i=0;i<totalfiles;i++)
		{
			sprintf(filename,"%srain%i",GlobalVars->rain_filename,GlobalVars->first_file+i);
			rainfile = fopen(filename,"r");
			if(!rainfile)
			{
				printf("Error opening file %s.\n",filename);
				return;
			}

			for(j=0;j<N;j++)
			{
				//Read in the storm data for this link
				fread(&rainfall_buffer,sizeof(float),1,rainfile);

				//This assumes different endianness
				holder = *(unsigned int*) &rainfall_buffer;
				holder = (((holder & 0x0000ffff)<<16) | ((holder & 0xffff0000)>>16));
				holder = (((holder & 0x00ff00ff)<<8) | ((holder & 0xff00ff00)>>8));
				rainfall_buffer = *(float*) &holder;

				if(i == 0 || !(.999 * rainfall_buffer <= last[j] && last[j] <= 1.001 * rainfall_buffer) )
				{
					last[j] = rainfall_buffer;
					total[j]++;
				}
			}

			fclose(rainfile);
		}
	}
	else if(GlobalVars->rain_flag == 1)
	{
		unsigned int id,numtimes;
		rainfile = fopen(GlobalVars->rain_filename,"r");
		if(!rainfile)
		{
			printf("Error opening file %s.\n",GlobalVars->rain_filename);
			return;
		}

		fscanf(rainfile,"%*u");	//N

		for(i=0;i<N;i++)
		{
			fscanf(rainfile,"%u\n%u",&id,&numtimes);
			total[i] = numtimes;
			for(j=0;j<numtimes;j++)
				fscanf(rainfile,"%*f %*f");
		}
	}
	else
	{
		printf("Invalid rain flag (%hu) in CreateGraphRain.\n",GlobalVars->rain_flag);
		return;
	}

	printf("Writing graph to file...\n");

	//Make the graph
	fprintf(output,"%u %u %u%u%u\n",N,N-1,0,1,0);
	for(i=0;i<N;i++)
	{
		fprintf(output,"%u ",total[i]);
		if(sys[i]->child != NULL)	fprintf(output,"%u ",sys[i]->child->location + offset);
		for(j=0;j<sys[i]->num_parents;j++)
			fprintf(output,"%u ",sys[i]->parents[j]->location + offset);
		fprintf(output,"\n");
	}

	fclose(output);
	free(total);
	free(last);
}
*/

/*
int* Partition_System_File(Link** sys,unsigned int N,Link** leaves,unsigned int numleaves,unsigned int** my_sys,unsigned int* my_N,unsigned int* my_max_nodes,TransData* my_data,short int *getting)
{
	unsigned int i,j,start_index,end_index,extras,partition,loc;
	unsigned int nodes_per_proc = numleaves/np;	//Number of leaves assigned to each process (except the last)
	char filename[256];

	start_index = nodes_per_proc * my_rank;
	if(my_rank == np-1)		end_index = numleaves;
	else				end_index = nodes_per_proc * (my_rank + 1);
//	*my_N = end_index - start_index;
	*my_N = 0;
	*my_max_nodes = N - numleaves + nodes_per_proc;
	*my_sys = (unsigned int*) malloc(*my_max_nodes * sizeof(unsigned int));	//The indices of this processes links (in sys)
	for(i=0;i<*my_max_nodes;i++)	(*my_sys)[i] = -1;
	for(i=start_index;i<end_index;i++)	(*my_sys)[i-start_index] = leaves[i]->location;
	for(i=0;i<N;i++)	getting[i] = 0;

	//Initialize assignments
	int* assignments = (int*) malloc(N*sizeof(int));
	for(i=0;i<N;i++)	assignments[i] = -1;

	//Open a file
	sprintf(filename,"testgraph.gra.part.%i",np);
	FILE* graph = fopen(filename,"r");
	if(!graph)	printf("Error reading graph file %s.\n",filename);

	//Read in the partitioning data
	//Note: The links should be listed in the graph file in the same order as in the .rvr file
	for(i=0;i<N;i++)
	{
		fscanf(graph,"%i",&(assignments[i]));
		if(assignments[i] == my_rank)
		{
			(*my_sys)[*my_N] = i;
			(*my_N)++;
		}
	}

	//Close file
	fclose(graph);

	//Set the getting array and determine number of sending and receiving links
	for(i=0;i<*my_N;i++)
	{
		//Receiving
		for(j=0;j<sys[(*my_sys)[i]]->num_parents;j++)
		{
			loc = sys[(*my_sys)[i]]->parents[j]->location;
			if(assignments[loc] != my_rank)
			{
				getting[loc] = 1;
				my_data->receive_size[assignments[loc]]++;
			}
		}

		//Sending
		if(sys[(*my_sys)[i]]->child != NULL)
		{
			loc = sys[(*my_sys)[i]]->child->location;
			if(assignments[loc] != my_rank)
				my_data->send_size[assignments[loc]]++;
		}
	}

	//Reorder my_sys so that the links with lower numbering are towards the beginning
	merge_sort_by_distance(sys,*my_sys,*my_N);

	//Allocate space in my_data for recieving and sending
	for(j=0;j<np;j++)
	{
		my_data->receive_data[j] = (Link**) malloc(my_data->receive_size[j] * sizeof(Link*));
		my_data->send_data[j] = (Link**) malloc(my_data->send_size[j] * sizeof(Link*));
	}

	//Set the receive_data and send_data arrays
	int* current_receive_size = (int*) calloc(np,sizeof(int));
	int* current_send_size = (int*) calloc(np,sizeof(int));
	for(i=0;i<*my_N;i++)
	{
		//Receiving
		for(j=0;j<sys[(*my_sys)[i]]->num_parents;j++)
		{
			loc = sys[(*my_sys)[i]]->parents[j]->location;
			if(assignments[loc] != my_rank)
			{
				my_data->receive_data[assignments[loc]][current_receive_size[assignments[loc]]] = sys[loc];
				current_receive_size[assignments[loc]]++;
			}
		}

		//Sending
		if(sys[(*my_sys)[i]]->child != NULL)
		{
			loc = sys[(*my_sys)[i]]->child->location;
			if(assignments[loc] != my_rank)
			{
				my_data->send_data[assignments[loc]][current_send_size[assignments[loc]]] = sys[(*my_sys)[i]];
				current_send_size[assignments[loc]]++;
			}
		}
	}

	free(current_receive_size);
	free(current_send_size);

	return assignments;
}
*/
