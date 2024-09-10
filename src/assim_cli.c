#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <mpi.h>

#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif

#if defined(HAVE_POSTGRESQL)
#include <libpq-fe.h>
#endif

#if defined(HAVE_PETSC)
#include <petsc.h>
#endif

#if !defined(_MSC_VER)
#define ASYNCH_SLEEP(seconds) sleep(seconds)
#else
#include <windows.h>
#define ASYNCH_SLEEP(seconds) Sleep((seconds) * 1000)
#endif

#include <asynch_interface.h>
#include <assim/linear_least_squares.h>
#include <assim/models.h>
#include <assim/ancillary.h>
#include <optparse.h>

// Global variables
bool verbose = false;
int my_rank;
int np;

////Output functions
//int Output_Linkid(double t, VEC y_i, VEC global_params, VEC params, int state, void* user);
//void Set_Output_User_LinkID(asynchsolver* asynch);
//
//int Output_Timestamp(double t, VEC y_i, VEC global_params, VEC params, int state, void* user);
//void Set_Output_User_Timestamp(asynchsolver* asynch);


//Print to stdout only for process of rank 0
int print_out(const char* format, ...)
{
    int res = 0;
    if (my_rank == 0)
    {
        va_list args;
        va_start(args, format);
        res = vprintf(format, args);
        va_end(args);
    }

    return res;
}


//Print to stderr only for process of rank 0
int print_err(const char* format, ...)
{
    int res = 0;
    if (my_rank == 0)
    {
        va_list args;
        va_start(args, format);
        res = vfprintf(stderr, format, args);
        va_end(args);
    }

    return res;
}


//Make sure we finalize MPI
void asynch_onexit(void)
{
    int flag;
    MPI_Finalized(&flag);
    if (!flag)
        MPI_Finalize();
}


int main(int argc, char* argv[])
{
    int res;

    //Initialize MPI stuff
    res = MPI_Init(&argc, &argv);
    if (res == MPI_SUCCESS)
        atexit(asynch_onexit);
    else
    {
        print_err("Failed to initialize MPI");
        exit(EXIT_FAILURE);
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    //PETSC
    PetscInitialize(&argc, &argv, NULL, NULL);

    //Command line options
    bool debug = false;
    bool help = false;
    bool version = false;

    //Parse command line
    struct optparse options;
    optparse_init(&options, argv);
    struct optparse_long longopts[] = {
        { "debug", 'd', OPTPARSE_NONE },
        { "help", 'h', OPTPARSE_NONE },
        { "version", 'v', OPTPARSE_NONE },
        { "verbose", 'w', OPTPARSE_NONE },
        { 0 }
    };
    int option;
    while ((option = optparse_long(&options, longopts, NULL)) != -1) {
        switch (option) {
        case 'd':
            debug = true;
            break;
        case 'h':
            help = true;
            break;
        case 'v':
            version = true;
            break;
        case 'w':
            verbose = true;
            break;
        case '?':
            print_err("%s: %s\n", argv[0], options.errmsg);
            exit(EXIT_FAILURE);
        }
    }

    if (version) print_out("This is %s\n", PACKAGE_STRING);
    if (help)
    {
        print_out("Usage: assim <global file> <das file>\n", PACKAGE_STRING);
        print_out(
            "  -d [--debug]   : Wait for the user input at the begining of the program (useful" \
            "                   for attaching a debugger)\n" \
            "  -w [--verbose] : Print debugging information to stdout\n" \
            "  -v [--version] : Print the current version of ASYNCH\n");
        exit(EXIT_SUCCESS);
    }
    if (version || help) exit(EXIT_SUCCESS);

    //Parse remaining arguments
    char *global_filename = optparse_arg(&options);
    if (global_filename == NULL && my_rank == 0)
    {
        print_err("Command line parameter required:  A universal variable file (.gbl).\n");
        exit(EXIT_FAILURE);
    }

    char *assim_filename = optparse_arg(&options);
    if (assim_filename == NULL && my_rank == 0)
    {
        print_err("Command line parameter required:  An assim file (.das).\n");
        exit(EXIT_FAILURE);
    }

    //Disable stdout buffering
    setvbuf(stdout, NULL, _IONBF, 0);

    if (debug)
    {
        //Disable stdout buffering
        setvbuf(stdout, NULL, _IONBF, 0);

        //When the program first starts to execute, at the very beginning of our program, we 
        //ask the user to type some sort of input to simply stall the application until start your
        //"Attach to Process" and you can attach to all the different threads in your program.
        if (my_rank == 0)
        {
            printf("You may now attach the debugger then press enter.\n");
            //fflush(stdout);
            int ch = getchar();
        }

        MPI_Barrier(MPI_COMM_WORLD); // All threads will wait here until you give thread 0 an input
    }

    //Declare variables
    unsigned int i, j/*, k,l,m,n,*/;
    //RKMethod** AllMethods;
    //char additional[16];	//For output filename
    //srand(time(NULL));	//!!!! Is this needed? !!!!

    //Init asynch object and the river network
    AsynchSolver *asynch = Asynch_Init(MPI_COMM_WORLD, verbose);

    //Create assim structure
    AssimData assim;
    //Read data assimilation file
    InitAssimData(&assim, assim_filename);

    //Model 254, full
    AsynchModel model_254_assim;
    memset(&model_254_assim, 0, sizeof(AsynchModel));
    model_254_assim.dim = 4;
    model_254_assim.set_param_sizes = SetParamSizes_Assim_254;
    model_254_assim.convert = ConvertParams_Assim_254;
    model_254_assim.routines = InitRoutines_Assim_254;
    model_254_assim.precalculations = Precalculations_Assim_254;
    model_254_assim.initialize_eqs = ReadInitData_Assim_254;

    //Model 254, q
    AsynchModel model_254_assim_q;
    memset(&model_254_assim_q, 0, sizeof(AsynchModel));
    model_254_assim_q.dim = 4;
    model_254_assim_q.set_param_sizes = SetParamSizes_Assim_254;
    model_254_assim_q.convert = ConvertParams_Assim_254;
    model_254_assim_q.routines = InitRoutines_Assim_254_q;
    model_254_assim_q.precalculations = Precalculations_Assim_254;
    model_254_assim_q.initialize_eqs = ReadInitData_Assim_254_q;

    //Model 254, q and s_p
    AsynchModel model_254_assim_qsp;
    memset(&model_254_assim_qsp, 0, sizeof(AsynchModel));
    model_254_assim_qsp.dim = 4;
    model_254_assim_qsp.set_param_sizes = SetParamSizes_Assim_254;
    model_254_assim_qsp.convert = ConvertParams_Assim_254;
    model_254_assim_qsp.routines = InitRoutines_Assim_254_qsp;
    model_254_assim_qsp.precalculations = Precalculations_Assim_254;
    model_254_assim_qsp.initialize_eqs = ReadInitData_Assim_254_qsp;

    //Model 254, q and s_t
    AsynchModel model_254_assim_qst;
    memset(&model_254_assim_qst, 0, sizeof(AsynchModel));
    model_254_assim_qst.dim = 4;
    model_254_assim_qst.set_param_sizes = SetParamSizes_Assim_254;
    model_254_assim_qst.convert = ConvertParams_Assim_254;
    model_254_assim_qst.routines = InitRoutines_Assim_254_qst;
    model_254_assim_qst.precalculations = Precalculations_Assim_254;
    model_254_assim_qst.initialize_eqs = ReadInitData_Assim_254_qst;

    if (strcmp(assim.model, "254") == 0)
        Asynch_Custom_Model(asynch, &model_254_assim);
    else if (strcmp(assim.model, "254_q") == 0)
        Asynch_Custom_Model(asynch, &model_254_assim_q);
    else if (strcmp(assim.model, "254_qsp") == 0)
        Asynch_Custom_Model(asynch, &model_254_assim_qsp);
    else if (strcmp(assim.model, "254_qst") == 0)
        Asynch_Custom_Model(asynch, &model_254_assim_qst);
    else /* default: */
    {
        print_err("Invalid model variant (expected 254, 254_q, 254_qsp or 254_qst an got %s", assim.model);
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    //For model 254
    //TODO Generalize this
    unsigned int problem_dim = 4;
    unsigned int assim_dim = 4;

    //Read global file
    if (my_rank == 0)	printf("Reading global file...\n");
    Asynch_Parse_GBL(asynch, global_filename);
    if (my_rank == 0)	printf("Loading network...\n");
    Asynch_Load_Network(asynch);

    //Find the gauged locations
    if (GetObservationsIds(asynch, &assim))
        MPI_Abort(MPI_COMM_WORLD, 1);

    //Finds the link ids upstreams from every gauged locations
    const bool trim = true;
    FindUpstreamLinks2(asynch, &assim, problem_dim, trim, assim.obs_time_step, assim.num_steps, assim.obs_locs, assim.num_obs);

    print_out("Partitioning network...\n");
    Asynch_Partition_Network(asynch);
    //CleanUpstreamLinks(asynch);
    print_out("Loading parameters...\n");
    Asynch_Load_Network_Parameters(asynch);
    print_out("Reading dam and reservoir data...\n");
    Asynch_Load_Dams(asynch);
    print_out("Setting up numerical error data...\n");
    Asynch_Load_Numerical_Error_Data(asynch);
    print_out("Initializing model...\n");
    Asynch_Initialize_Model(asynch);
    Setup_Errors(asynch, problem_dim);
    print_out("Loading initial conditions...\n");
    Asynch_Load_Initial_Conditions(asynch);
    print_out("Loading forcings...\n");
    Asynch_Load_Forcings(asynch);
    print_out("Loading output data information...\n");
    Asynch_Load_Save_Lists(asynch);
    print_out("Finalizing network...\n");
    Asynch_Finalize_Network(asynch);
    print_out("Calculating initial step sizes...\n");
    Asynch_Calculate_Step_Sizes(asynch);

    //Pull data from asynch
    Link **my_sys = asynch->my_sys;
    Lookup *id_to_loc = asynch->id_to_loc;
    unsigned int my_N = asynch->my_N, N = asynch->N, num_forcings = asynch->globals->num_forcings;
    int *assignments = asynch->assignments;
    Link* sys = asynch->sys;
    short int *getting = asynch->getting;
    GlobalVars *globals = asynch->globals;
    AsynchModel* custom_model = asynch->model;

    //Initialize choices
    unsigned int num_total_obs = assim.num_steps * assim.num_obs;
    time_t begin_time = Asynch_Get_Begin_Timestamp(asynch);
    time_t end_time = Asynch_Get_End_Timestamp(asynch);
    double duration = Asynch_Get_Total_Simulation_Duration(asynch);
    double t_b = 0.0;
    unsigned int allstates = assim_dim * N;

    // Allocate background
    double  *x_b = calloc(allstates, sizeof(double));
    for (i = 0; i < N; i++)
    {
        if (assignments[i] == my_rank)
        {
            for (j = 0; j < assim_dim; j++)
                x_b[i*assim_dim + j] = sys[i].my->list.tail->y_approx[j];	//!!!! Need to be able to specify which states are used !!!!
        }
        //MPI_Bcast(&(x_b[i*assim_dim]), assim_dim, MPI_DOUBLE, assignments[i], MPI_COMM_WORLD);
    }

    int mpi_res = MPI_Allreduce(MPI_IN_PLACE, x_b, allstates, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // Vector of observations (single step)
    double* d = calloc(assim.num_obs, sizeof(double));

    // Vector of observations (multiple steps)
    double* d_full = calloc(num_total_obs, sizeof(double));

    //Values used to start asynch solver in tao solvers
    double* x_start = calloc(allstates, sizeof(double));	//Values used to start asynch solver in tao solvers

    //Call model specific data assimilation routines
    if (strcmp(assim.model, "254") == 0)
        //For Model 254
        Setup_Fitting_Data_Model254(asynch, assim.obs_locs, assim.num_obs);
    else if (strcmp(assim.model, "254_q") == 0)
        //For Model 254 trim, q
        Setup_Fitting_Data_Model254_q(asynch, assim.obs_locs, assim.num_obs);
    else if (strcmp(assim.model, "254_qsp") == 0)
        //For Model 254 trim, q and s_p
        Setup_Fitting_Data_Model254_qsp(asynch, assim.obs_locs, assim.num_obs);
    else if (strcmp(assim.model, "254_qst") == 0)
        //For Model 254 trim, q and s_t
        Setup_Fitting_Data_Model254_qst(asynch, assim.obs_locs, assim.num_obs);

    //Find locations unaffected by gauges
    unsigned int *vareq_shift, *inv_vareq_shift;
    unsigned int allstates_needed = BuildStateShift(asynch, allstates, assim.obs_locs, assim.num_obs, &vareq_shift, &inv_vareq_shift);

    printf("allstates_needed: %u allstates: %u\n", allstates_needed, allstates);

    /*
        //!!!! Assuming q and s_p are changing !!!!
        unsigned int curr_idx = 0;
        Vec B,R;
        if(my_rank == 0)
        {
            VecCreateSeq(MPI_COMM_SELF,allstates_needed,&B);
            for(i=0;i<N;i++)
            {
                if(links_needed[i])
                {
                    VecSetValue(B,curr_idx++,1.0,INSERT_VALUES);
                    VecSetValue(B,curr_idx++,1e2,INSERT_VALUES);
                }
            }
            VecAssemblyBegin(B);
            VecAssemblyEnd(B);

            VecCreateSeq(MPI_COMM_SELF,num_obs*steps_to_use,&R);
            for(i=0;i<num_obs;i++)
            {
                for(j=0;j<steps_to_use;j++)
                    VecSetValue(R,j*steps_to_use + i,invupareas[obs_locs[i]] * 10.0,INSERT_VALUES);
            }
            VecAssemblyBegin(R);
            VecAssemblyEnd(R);
        }
        free(links_needed);
    */


    //Prep PetSC
    if (my_rank == 0)
        printf("\nPrepping PetSc...\n");
    
    AssimWorkspace ws;

    //For linear least squares
    int *HM_col_indices = NULL;
    if (my_rank == 0)
    {
        VecCreateSeq(MPI_COMM_SELF, allstates_needed, &ws.rhs);
        VecCreateSeq(MPI_COMM_SELF, allstates_needed, &ws.x);
        MatCreateSeqDense(MPI_COMM_SELF, num_total_obs, allstates_needed, NULL, &ws.HM);
        MatCreateSeqDense(MPI_COMM_SELF, allstates_needed, allstates_needed, NULL, &ws.HTH);
        MatCreateSeqDense(MPI_COMM_SELF, allstates_needed, num_total_obs, NULL, &ws.HMTR);
        HM_col_indices = (int*)malloc(allstates_needed * sizeof(int));
        for (i = 0; i < allstates_needed; i++)
            HM_col_indices[i] = i;
        KSPCreate(MPI_COMM_SELF, &ws.ksp);
        KSPSetOperators(ws.ksp, ws.HTH, ws.HTH);
        //KSPSetTolerances(ws.ksp, 1.e-12, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
        KSPSetFromOptions(ws.ksp);	// This is used to override the solver setting from the command line
    }

    int* d_indices = (int*)malloc(num_total_obs * sizeof(int));
    for (i = 0; i < num_total_obs; i++)
        d_indices[i] = i;

    //Transfer upstreams areas to all procs
    double* inv_upareas = (double*)calloc(N, sizeof(double));
    UpstreamData* updata;
    for (i = 0; i < N; i++)
    {
        if (assignments[i] == my_rank)
            inv_upareas[i] = 1.0 / (sys[i].params[globals->area_idx] * 1e3);
    }

    MPI_Allreduce(MPI_IN_PLACE, inv_upareas, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    //Links needed for fitting
    bool *links_needed = (bool*)calloc(N, sizeof(bool));
    for (i = 0; i < assim.num_obs; i++)
    {
        if (assignments[assim.obs_locs[i]] == my_rank)
        {
            Link *current = &sys[assim.obs_locs[i]];
            updata = (UpstreamData*)(current->user);
            links_needed[current->location] = true;
            //for(j=0;j<current->num_parents;j++)
            //{
            //	for(k=0;k<updata->num_upstreams[j];k++)
            //		my_links_needed[updata->upstreams[j][k]] = 1;
            //}

            for (j = 0; j < updata->num_upstreams; j++)
            {
                links_needed[updata->upstreams[j]->location] = true;
            }
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, links_needed, N, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

    //Build weight matrices
    //!!!! Assuming only q is changing !!!!
    unsigned int curr_idx = 0;
    if (my_rank == 0)
    {
        VecCreateSeq(MPI_COMM_SELF, allstates_needed, &ws.B);
        for (i = 0; i < N; i++)
        {
            if (links_needed[i])
                //VecSetValue(ws.B, curr_idx++, inv_upareas[i], INSERT_VALUES);
                VecSetValue(ws.B, curr_idx++, 1.0, INSERT_VALUES);
        }
        VecAssemblyBegin(ws.B);
        VecAssemblyEnd(ws.B);

        VecCreateSeq(MPI_COMM_SELF, num_total_obs, &ws.R);
        for (i = 0; i < assim.num_obs; i++)
        {
            for (j = 0; j < assim.num_steps; j++)
                //VecSetValue(ws.R, j*num_obs + i, inv_upareas[obs_locs[i]] * 10.0, INSERT_VALUES);
                VecSetValue(ws.R, j * assim.num_obs + i, 1.0, INSERT_VALUES);
        }
        VecAssemblyBegin(ws.R);
        VecAssemblyEnd(ws.R);

        if (verbose)
        {
            printf("Weighting Matrix B (diagonal)\n");
            VecView(ws.B, PETSC_VIEWER_STDOUT_SELF);
            printf("Weighting Matrix R (diagonal)\n");
            VecView(ws.R, PETSC_VIEWER_STDOUT_SELF);
        }
    }
    free(links_needed);

    ws.HM_buffer = (double*)calloc(allstates_needed, sizeof(double));
    ws.HM_col_indices = HM_col_indices;
    ws.d_indices = d_indices;
    ws.d_full = d_full;
    ws.x_start = x_start;
    ws.problem_dim = problem_dim;
    ws.assim_dim = assim_dim;
    ws.allstates = allstates;
    ws.allstates_needed = allstates_needed;
    ws.vareq_shift = vareq_shift;
    ws.inv_vareq_shift = inv_vareq_shift;
    ws.obs_time_step = assim.obs_time_step;
    ws.num_steps = assim.num_steps;
    ws.obs_locs = assim.obs_locs;
    ws.num_obs = assim.num_obs;
    ws.t_b = t_b;
    ws.x_b = x_b;

    //Print out some information
    unsigned int my_eqs = 0, total_eqs;
    for (i = 0; i < my_N; i++)	my_eqs += my_sys[i]->dim;
    MPI_Reduce(&my_eqs, &total_eqs, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

    printf("[%i]: Good to go with %u links (%u eqs).\n", my_rank, my_N, my_eqs);
    if (my_rank == 0)
    {
        ASYNCH_SLEEP(1);
        printf("\nNetwork has a total of %u links and %u equations.\n\n", N, total_eqs);
        printf("Making calculations...\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    unsigned int forcing_idx_rain = 0;
    unsigned int forcing_idx_tep = 1;

    unsigned int begin_assim_window = asynch->forcings[forcing_idx_rain].first_file;
    unsigned int end_assim_window = asynch->forcings[forcing_idx_rain].first_file + (unsigned int)(assim.num_steps * assim.obs_time_step * 60.0);

    //Reset each link
    //Asynch_Set_System_State(asynch, 0.0, backup);
    //Set_Output_User_forecastparams(asynch, forecast_time_unix, simulation_time_with_data_secs);
    //Set_Output_PeakflowUser_Offset(asynch,forecast_time_unix);
    //Asynch_Write_Current_Step(asynch);
    Asynch_Set_Forcing_State(asynch, forcing_idx_rain, 0.0, begin_assim_window, end_assim_window);

    for (i = 0; i < asynch->globals->num_forcings; i++)	//Set any other database forcings to begin at first_file
    {
        if (asynch->forcings[i].flag == 3)
            Asynch_Set_Forcing_State(asynch, i, 0.0, begin_assim_window, end_assim_window);
    }

    //Make sure all buffer flushing is done
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

    //Start the analysis
    unsigned int max_least_squares_iters = assim.max_least_squares_iters;
    double *analysis = calloc(allstates, sizeof(double));   //!!!! Should be removed !!!!
    double *q = calloc(num_total_obs, sizeof(double));

    //Get the observations
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if (my_rank == 0)
            printf("Downloading observations...\n");

        double start = MPI_Wtime();

        while (GetObservationsData(&assim, id_to_loc, N, begin_assim_window, d_full) == -1)
        {
            if (my_rank == 0)	printf("Error downloading observations. Retrying...\n");
            ASYNCH_SLEEP(5);
            MPI_Barrier(MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        double stop = MPI_Wtime();

        if (my_rank == 0)
            printf("Time to get new discharges: %.0f\n", stop - start);
    }

    if (verbose && my_rank == 0)
    {
        printf("d_full\n");
        Print_VECTOR(d_full, num_total_obs);
        printf("\n");
    }

    ////Scale by upstream area
    //printf("Scaling init discharges by upstream area...\n");
    //AdjustDischarges(asynch, obs_locs, d_full, num_obs, assim_dim, x_b);

    //Set the initial guess from background
    memcpy(x_start, x_b, allstates * sizeof(double));
    //Copy in states that won't be affected by the optimization solver
    memcpy(analysis, x_b, allstates * sizeof(double));

    //Calculate the analysis
    bool try_again = false;
    do
    {
        int iterations = 0;
        double error = 1.0, prev_error = -1.0;
        for (j = 0; (j < max_least_squares_iters) && (error > 1e-2); j++)
        {
            iterations++;

            //Solve the Linear sytems
            LSSolveSys(asynch, &ws, q);

            //Compute the distance between observations and model
            error = LSComputeDistance(d_full, q, num_total_obs);
            if (prev_error >= 0.0)
            {
                //Check if the distance is decreasing
                double diff = prev_error - error;
                if (error > prev_error)
                {
                    if (my_rank == 0)
                    {
                        printf("!!!! LS error got worse. Breaking... !!!!\n");
                        printf("Errors are %f and %f\n", error, prev_error);
                    }

                    //Go back to previous solution
                    for (i = 0; i < allstates; i++)
                        x_start[i] = analysis[i];

                    break;
                }
                if (my_rank == 0)
                    printf("Difference is %f (%f vs %f)\n", diff, error, prev_error);
            }

            prev_error = error;
            for (i = 0; i < N; i++)
            {
                //TODO Add test for negative discharge
                analysis[i] = x_start[i] > 1.e-14 ? x_start[i] : 1.e-14;
                analysis[i + 1] = x_start[i + 1] > 0. ? x_start[i + 1] : 0.;
                analysis[i + 2] = x_start[i + 2] > 0. ? x_start[i + 2] : 0.;
                analysis[i + 3] = x_start[i + 3] > 0. ? x_start[i + 3] : 0.;
            }
        }
        if (my_rank == 0)
            printf("Total iterations = %i\n", iterations);

        //Check if the numerical scheme is having convergence issues
        if (!try_again)
            //Where dicharge values are having problems converging, upstreams discharge are cut in half
            try_again = ReduceBadDischargeValues(sys, assignments, N, d_full, q, assim.num_steps, assim.obs_locs, assim.num_obs, x_start, assim_dim, 1.0);	//!!!! Not sure what to use for the limit... !!!!
        else
            try_again = false;
    } while (try_again);

    if (verbose && my_rank == 0)
    {
        //if(k == 10)
        {
            printf("x_b\n");
            Print_VECTOR(x_b, allstates);
            printf("\n");

            printf("analysis [%d - %d]\n", begin_assim_window, end_assim_window);
            Print_VECTOR(analysis, allstates);
        }
    }

    free(analysis);
    free(q);

    double stop = MPI_Wtime();
    print_out("\nTime for calculations: %f. All done!\n", stop - start);

    //Prepare snapshots
    for (i = 0; i < N; i++)
    {
        Link *current = &asynch->sys[i];
        if (current->my != NULL)
        {
            double *y = current->my->list.tail->y_approx;
            for (j = 0; j < problem_dim; j++)
                y[j] = x_start[i * problem_dim + j];

            // Idem CheckConsistency_Nonzero_4States
            if (y[0] < 1e-14)	y[0] = 1e-14;
            if (y[1] < 1e-20)	y[1] = 0.0;
            if (y[2] < 1e-20)	y[2] = 0.0;
            if (y[3] < 1e-20)	y[3] = 0.0;
        }
    }

    //Make a snaphsot
    print_out("Making snapshot\n");
    Asynch_Take_System_Snapshot(asynch, NULL);

    //Clean up
    print_out("Cleaning up\n");
    free(d);
    free(d_full);
    free(x_start);
    free(vareq_shift);
    free(inv_vareq_shift);

    if (my_rank == 0)
    {
        MatDestroy(&ws.HM);
        MatDestroy(&ws.HTH);
        VecDestroy(&ws.rhs);
        VecDestroy(&ws.x);
        VecDestroy(&ws.R);
        VecDestroy(&ws.B);
        MatDestroy(&ws.HMTR);
        KSPDestroy(&ws.ksp);
    }
    free(HM_col_indices);
    free(d_indices);

    free(inv_upareas);
    FreeAssimData(&assim);

    //Petsc clean up
    PetscFinalize();

    //Asynch clean up
    FreeUpstreamLinks(asynch);
    Asynch_Delete_Temporary_Files(asynch);
    Asynch_Free(asynch);

    return res;
}
