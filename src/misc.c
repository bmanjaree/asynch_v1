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



