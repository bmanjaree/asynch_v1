from asynch_py import ASYNCH_LIBRARY_LOCATION
from ctypes import *
from mpi4py import MPI
import numpy
import sys

#ASYNCH_LIBRARY_LOCATION = '/home/ssma/NewAsynchVersion/libs/libasynch_py.so'

#Predefined values *******************************************************
ASYNCH_BAD_TYPE = -1
ASYNCH_CHAR = 0
ASYNCH_SHORT = 1
ASYNCH_INT = 2
ASYNCH_FLOAT = 3
ASYNCH_DOUBLE = 4

ASYNCH_MAX_QUERY_SIZE = 1024
ASYNCH_MAX_DB_CONNECTIONS = 16
ASYNCH_DB_LOC_TOPO = 0
ASYNCH_DB_LOC_PARAMS = 1
ASYNCH_DB_LOC_INIT = 2
ASYNCH_DB_LOC_QVS = 3
ASYNCH_DB_LOC_RSV = 4
ASYNCH_DB_LOC_HYDROSAVE = 5
ASYNCH_DB_LOC_PEAKSAVE = 6
ASYNCH_DB_LOC_HYDRO_OUTPUT = 7
ASYNCH_DB_LOC_PEAK_OUTPUT = 8
ASYNCH_DB_LOC_SNAPSHOT_OUTPUT = 9
ASYNCH_DB_LOC_FORCING_START = 10

#Data structures **********************************************************
class FILE(Structure): pass
class RKMethod(Structure): pass
class RKSolutionList(Structure): pass
class ErrorData(Structure): pass
class ConnData(Structure): pass
class Forcing(Structure): pass
class TempStorage(Structure): pass
class io(Structure): pass
class Link(Structure): pass
class UnivVars(Structure): pass
class ForcingData(Structure): pass

class VEC(Structure):
	_fields_ = [("ve",POINTER(c_double)),("dim",c_uint)]

class MAT(Structure):
	_fields_ = [("array",POINTER(c_double)),("me",POINTER(POINTER(c_double))),("m",c_uint),("n",c_uint)]

class QVSData(Structure):
	_fields_ = [("points",POINTER(POINTER(c_double))),("n_values",c_uint)]

ASYNCH_F_DATATYPE = CFUNCTYPE(None,c_double,POINTER(VEC),POINTER(POINTER(VEC)),c_ushort,POINTER(VEC),POINTER(c_double),POINTER(QVSData),POINTER(VEC),c_int,c_void_p,POINTER(VEC))
ASYNCH_RKSOLVER_DATATYPE = CFUNCTYPE(c_int,POINTER(Link),POINTER(UnivVars),POINTER(c_int),c_short,POINTER(FILE),POINTER(ConnData),POINTER(POINTER(Forcing)),POINTER(TempStorage))
ASYNCH_CONSISTENCY_DATATYPE = CFUNCTYPE(None,POINTER(VEC),POINTER(VEC),POINTER(VEC))
ASYNCH_ALG_DATATYPE = CFUNCTYPE(None,POINTER(VEC),POINTER(VEC),POINTER(VEC),POINTER(QVSData),c_int,c_void_p,POINTER(VEC))
ASYNCH_STATECHECK_DATATYPE = CFUNCTYPE(c_int,POINTER(VEC),POINTER(VEC),POINTER(VEC),POINTER(QVSData),c_uint)

ASYNCH_OUTPUT_DATATYPE = CFUNCTYPE(None,c_double,POINTER(VEC),POINTER(VEC),POINTER(VEC),c_int,c_void_p)
ASYNCH_OUTPUT_INT_DATATYPE = CFUNCTYPE(c_int,c_double,POINTER(VEC),POINTER(VEC),POINTER(VEC),c_int,c_void_p)
ASYNCH_OUTPUT_DOUBLE_DATATYPE = CFUNCTYPE(c_double,c_double,POINTER(VEC),POINTER(VEC),POINTER(VEC),c_int,c_void_p)

ASYNCH_PEAKOUTPUT_DATATYPE = CFUNCTYPE(None,c_uint,c_double,POINTER(VEC),POINTER(VEC),POINTER(VEC),c_double,c_uint,c_void_p,POINTER(c_char))

class UnivVars(Structure):
	_fields_ = [("type",c_ushort),("method",c_ushort),("max_s",c_ushort),("max_parents",c_ushort),
			("iter_limit",c_int),("max_transfer_steps",c_int),("maxtime",c_double),("t_0",c_double),
			("discont_size",c_uint),("max_localorder",c_uint),("uses_dam",c_ushort),("global_params",POINTER(VEC)),
			("params_size",c_uint),("dam_params_size",c_uint),("disk_params",c_uint),("area_idx",c_uint),
			("areah_idx",c_uint),("rain_filename",c_char_p),("init_filename",c_char_p),("rvr_filename",c_char_p),
			("prm_filename",c_char_p),("init_flag",c_ushort),("rvr_flag",c_ushort),("prm_flag",c_ushort),
			("output_flag",c_ushort),("temp_filename",c_char_p),("dam_filename",c_char_p),("print_time",c_double),
			("print_par_flag",c_ushort),("dam_flag",c_ushort),("hydrosave_flag",c_ushort),("peaksave_flag",c_ushort),
			("hydrosave_filename",c_char_p),("peaksave_filename",c_char_p),("peakfilename",c_char_p),("max_dim",c_uint),
			("outletlink",c_uint),("string_size",c_uint),("query_size",c_uint),("rkd_flag",c_short),
			("convertarea_flag",c_ushort),("discont_tol",c_double),("min_error_tolerances",c_uint),("num_forcings",c_uint),
			("hydros_loc_flag",c_ushort),("peaks_loc_flag",c_ushort),("dump_loc_flag",c_ushort),("res_flag",c_ushort),
			("hydros_loc_filename",c_char_p),("peaks_loc_filename",c_char_p),("dump_loc_filename",c_char_p),("rsv_filename",c_char_p),
			("init_timestamp",c_uint),("res_forcing_idx",c_ushort),("num_states_for_printing",c_uint),("num_print",c_uint),
			("print_indices",POINTER(c_uint)),
			("outputs_d",POINTER(CFUNCTYPE(c_double,c_double,POINTER(VEC),POINTER(VEC),POINTER(VEC),c_int,c_void_p))),
			("outputs_i",POINTER(CFUNCTYPE(c_int,c_double,POINTER(VEC),POINTER(VEC),POINTER(VEC),c_int,c_void_p))),
			("output_names",POINTER(c_char_p)),("output_specifiers",POINTER(c_char_p)),("output_types",POINTER(c_short)),("output_sizes",POINTER(c_short)),
			("output_data",POINTER(io)),("peakflow_function_name",c_char_p),
			("peakflow_output",CFUNCTYPE(None,c_uint,c_double,POINTER(VEC),POINTER(VEC),POINTER(VEC),c_double,c_uint,c_void_p,c_char_p)),
			("hydro_table",c_char_p),("peak_table",c_char_p),("dump_table",c_char_p)]

class Link(Structure):
	_fields_ = [("method",POINTER(RKMethod)),("list",POINTER(RKSolutionList)),("errorinfo",POINTER(ErrorData)),("params",POINTER(VEC)),
			("f",ASYNCH_F_DATATYPE),
			("alg",ASYNCH_ALG_DATATYPE),
			("state_check",ASYNCH_STATECHECK_DATATYPE),
			("Jacobian",CFUNCTYPE(None,c_double,POINTER(VEC),POINTER(POINTER(VEC)),c_ushort,POINTER(VEC),POINTER(c_double),POINTER(VEC),POINTER(MAT))),
			("RKSolver",ASYNCH_RKSOLVER_DATATYPE),
			("CheckConsistency",ASYNCH_CONSISTENCY_DATATYPE),
			("h",c_double),("last_t",c_double),("print_time",c_double),("next_save",c_double),
			("ID",c_uint),("location",c_uint),("ready",c_short),("numparents",c_ushort),
			("disk_iterations",c_int),("peak_time",c_double),("peak_value",POINTER(VEC)),("parents",POINTER(POINTER(Link))),
			("c",POINTER(Link)),("current_iterations",c_int),("steps_on_diff_proc",c_int),("iters_removed",c_int),
			("distance",c_uint),("rejected",c_int),("save_flag",c_ushort),("peak_flag",c_ushort),
			("qvs",POINTER(QVSData)),("pos_offset",c_long),
			("expected_file_vals",c_uint),("dam",c_ushort),("res",c_ushort),
			("dim",c_uint),("diff_start",c_uint),("no_ini_start",c_uint),("num_dense",c_uint),("dense_indices",POINTER(c_uint)),
			("forcing_buff",POINTER(POINTER(ForcingData))),
			("forcing_change_times",POINTER(c_double)),("forcing_values",POINTER(c_double)),("forcing_indices",POINTER(c_uint)),
			("output_user",c_void_p),("peakoutput_user",c_void_p),
			("user",c_void_p),
			("last_eta",c_double),
			("JMatrix",POINTER(MAT)),("CoefMat",POINTER(MAT)),("Z_i",POINTER(POINTER(VEC))),("sol_diff",POINTER(VEC)),
			("h_old",c_double),("value_old",c_double),("compute_J",c_short),("compute_LU",c_short),
			("state",c_int),("discont",POINTER(c_double)),("discont_count",c_uint),("discont_start",c_uint),
			("discont_end",c_uint),("discont_send_count",c_uint),("discont_send",POINTER(c_double)),("discont_order_send",POINTER(c_uint))]

ASYNCH_SETPARAMSIZES_DATATYPE = CFUNCTYPE(None,POINTER(UnivVars),c_void_p)
ASYNCH_CONVERT_DATATYPE = CFUNCTYPE(None,POINTER(VEC),c_uint,c_void_p)
ASYNCH_ROUTINES_DATATYPE = CFUNCTYPE(None,POINTER(Link),c_uint,c_uint,c_ushort,c_void_p)
ASYNCH_PRECALCULATIONS_DATATYPE = CFUNCTYPE(None,POINTER(Link),POINTER(VEC),POINTER(VEC),c_uint,c_uint,c_ushort,c_uint,c_void_p)
ASYNCH_INITIALIZEEQS_DATATYPE = CFUNCTYPE(c_int,POINTER(VEC),POINTER(VEC),POINTER(QVSData),c_ushort,POINTER(VEC),c_uint,c_uint,c_uint,c_void_p,c_void_p)

# asynchsolver class with interface functions **********************************************************

class asynchsolver:
	def __init__(self):
		#self.lib = pydll.LoadLibrary('./libs/libasynch_py.so')
		self.lib = cdll.LoadLibrary(ASYNCH_LIBRARY_LOCATION)
		self.comm = MPI.COMM_WORLD
		self.np = self.comm.Get_size()
		self.my_rank = self.comm.Get_rank()
		ranks = range(0,self.np)
		self.asynch_obj = self.lib.Asynch_Init_py(self.np,(c_int * self.np)(*ranks))
		self.tempfiles_exist = False
		#ranks = numpy.zeros((self.np,1),numpy.dtype('i4'))
		#for i in range(self.np):	ranks[i] = i
		#self.asynch_obj = self.lib.Asynch_Init_py(self.np,ranks.ctypes.data)

		#Set return types of functions. Default is c_int.
		self.lib.Asynch_Parse_GBL.restype = None
		self.lib.Asynch_Load_Network.restype = None
		self.lib.Asynch_Partition_Network.restype = None
		self.lib.Asynch_Load_Network_Parameters.restype = None
		self.lib.Asynch_Load_Numerical_Error_Data.restype = None
		self.lib.Asynch_Initialize_Model.restype = None
		self.lib.Asynch_Load_Initial_Conditions.restype = None
		self.lib.Asynch_Load_Forcings.restype = None
		self.lib.Asynch_Load_Dams.restype = None
		self.lib.Asynch_Load_Save_Lists.restype = None
		self.lib.Asynch_Finalize_Network.restype = None
		self.lib.Asynch_Calculate_Step_Sizes.restype = None
		self.lib.Asynch_Free.restype = None
		self.lib.Asynch_Advance.restype = None
		self.lib.Asynch_Prepare_Output.restype = None
		self.lib.Asynch_Prepare_Temp_Files.restype = None
		self.lib.Asynch_Prepare_Peakflow_Output.restype = None
		self.lib.Asynch_Set_Database_Connection.restype = None
		self.lib.Asynch_Get_Total_Simulation_Time.restype = c_double
		self.lib.Asynch_Set_Total_Simulation_Time.restype = None
		self.lib.Asynch_Get_Last_Rainfall_Timestamp.restype = c_uint
		self.lib.Asynch_Set_Last_Rainfall_Timestamp.restype = None
		self.lib.Asynch_Get_First_Rainfall_Timestamp.restype = c_uint
		self.lib.Asynch_Set_First_Rainfall_Timestamp.restype = None
		self.lib.Asynch_Set_RainDB_Starttime.restype = None
		self.lib.Asynch_Set_Init_File.restype = None
		self.lib.Asynch_Get_Number_Links.restype = c_uint
		self.lib.Asynch_Get_Local_Number_Links.restype = c_uint
		self.lib.Asynch_Set_System_State_py.restype = None
		self.lib.Asynch_Reset_Peakflow_Data.restype = None
		self.lib.Asynch_Get_Local_LinkID.restype = c_uint
		self.lib.Asynch_Get_Init_Timestamp.restype = c_uint
		self.lib.Asynch_Copy_Local_OutputUser_Data.restype = None
		self.lib.Asynch_Set_Size_Local_OutputUser_Data.restype = None
		self.lib.Asynch_Get_Size_Global_Parameters.restype = c_uint

		#Functions created specifically for the Python interface
		self.lib.C_inc_ref.restype = None
		self.lib.Allocate_CUINT_Array.restype = POINTER(c_uint)
		self.lib.Free_PythonInterface.restype = None
		self.lib.SetParamSizes_py.restype = None
		self.lib.InitRoutines_py.restype = None
		self.lib.Asynch_Copy_Local_OutputUser_Data_py.restype = None

	def __del__(self):
		if self.tempfiles_exist == True:
			self.Delete_Temporary_Files()

		self.lib.Free_PythonInterface(self.asynch_obj)
		self.lib.Asynch_Free(self.asynch_obj)

	def Custom_Model(self,SetParamSizes,Convert,Routines,Precalculations,InitializeEqs):
		return self.lib.Asynch_Custom_Model_py(self.asynch_obj,SetParamSizes,Convert,Routines,Precalculations,InitializeEqs,py_object(self.lib))

	#Routines to initialize the system
	def Parse_GBL(self,gbl_filename):
		self.lib.Asynch_Parse_GBL(self.asynch_obj,gbl_filename)

	def Load_Network(self):
		self.lib.Asynch_Load_Network(self.asynch_obj)

	def Partition_Network(self):
		self.lib.Asynch_Partition_Network(self.asynch_obj)

	def Load_Network_Parameters(self,load_all):
		if load_all == True:
			c_load_all = 1
		else:
			c_load_all = 0
		self.lib.Asynch_Load_Network_Parameters(self.asynch_obj,c_load_all)

	def Load_Numerical_Error_Data(self):
		self.lib.Asynch_Load_Numerical_Error_Data(self.asynch_obj)

	def Initialize_Model(self):
		self.lib.Asynch_Initialize_Model(self.asynch_obj)

	def Load_Initial_Conditions(self):
		self.lib.Asynch_Load_Initial_Conditions(self.asynch_obj)

	def Load_Forcings(self):
		self.lib.Asynch_Load_Forcings(self.asynch_obj)

	def Load_Dams(self):
		self.lib.Asynch_Load_Dams(self.asynch_obj)

	def Load_Save_Lists(self):
		self.lib.Asynch_Load_Save_Lists(self.asynch_obj)

	def Finalize_Network(self):
		self.lib.Asynch_Finalize_Network(self.asynch_obj)

	def Calculate_Step_Sizes(self):
		self.lib.Asynch_Calculate_Step_Sizes(self.asynch_obj)

	#Forcing routines
	def Activate_Forcing(self,idx):
		return self.lib.Asynch_Activate_Forcing(self.asynch_obj,idx)

	def Deactivate_Forcing(self,idx):
		return self.lib.Asynch_Deactivate_Forcing(self.asynch_obj,idx)

	#Advance solver
	def Advance(self,print_flag):
		if print_flag == True:
			c_print_flag = 1
		else:
			c_print_flag = 0
		self.lib.Asynch_Advance(self.asynch_obj,c_print_flag)

	#Data file routines
	def Prepare_Output(self):
		self.lib.Asynch_Prepare_Output(self.asynch_obj)

	def Prepare_Temp_Files(self):
		self.tempfiles_exist = True
		self.lib.Asynch_Prepare_Temp_Files(self.asynch_obj)

	def Prepare_Peakflow_Output(self):
		self.lib.Asynch_Prepare_Peakflow_Output(self.asynch_obj)

	def Create_Output(self,additional_out):
		return self.lib.Asynch_Create_Output(self.asynch_obj,additional_out)

	def Create_Peakflows_Output(self):
		return self.lib.Asynch_Create_Peakflows_Output(self.asynch_obj)

	def Delete_Temporary_Files(self):
		self.tempfiles_exist = False
		return self.lib.Asynch_Delete_Temporary_Files(self.asynch_obj)

	def Write_Current_Step(self):
		return self.lib.Asynch_Write_Current_Step(self.asynch_obj)

	#Snapshot
	def Take_System_Snapshot(self,name):
		return self.lib.Asynch_Take_System_Snapshot(self.asynch_obj,name)

	#Set and get routines
	def Set_Database_Connection(self,database_info,conn_idx):
		self.lib.Asynch_Set_Database_Connection(self.asynch,database_info,conn_idx)

	def Get_Total_Simulation_Time(self):
		return self.lib.Asynch_Get_Total_Simulation_Time(self.asynch_obj)

	def Set_Total_Simulation_Time(self,new_time):
		self.lib.Asynch_Set_Total_Simulation_Time(self.asynch_obj,c_double(new_time));

	def Get_Last_Rainfall_Timestamp(self,forcing_idx):
		return self.lib.Asynch_Get_Last_Rainfall_Timestamp(self.asynch_obj,forcing_idx)

	def Set_Last_Rainfall_Timestamp(self,epoch_timestamp,forcing_idx):
		self.lib.Asynch_Set_Last_Rainfall_Timestamp(self.asynch_obj,epoch_timestamp,forcing_idx)

	def Get_First_Rainfall_Timestamp(self,forcing_idx):
		return self.lib.Asynch_Get_First_Rainfall_Timestamp(self.asynch_obj,forcing_idx)

	def Set_First_Rainfall_Timestamp(self,epoch_timestamp,forcing_idx):
		self.lib.Asynch_Set_First_Rainfall_Timestamp(self.asynch_obj,epoch_timestamp,forcing_idx)

	def Set_RainDB_Starttime(self,epoch_timestamp,forcing_idx):
		self.lib.Asynch_Set_RainDB_Starttime(self.asynch_obj,epoch_timestamp,forcing_idx)

	def Set_Init_File(self,filename):
		self.lib.Asynch_Set_Init_File(self.asynch_obj,filename)

	def Get_Number_Links(self):
		return self.lib.Asynch_Get_Number_Links(self.asynch_obj)

	def Get_Local_Number_Links(self):
		return self.lib.Asynch_Get_Local_Number_Links(self.asynch_obj)

	def Get_Local_LinkID(self,location):
		return self.lib.Asynch_Get_Local_LinkID(self.asynch_obj,location)

	def Set_Init_Timestamp(self,epoch_timestamp):
		return self.lib.Asynch_Set_Init_Timestamp(self.asynch_obj,epoch_timestamp)

	def Get_Init_Timestamp(self):
		return self.lib.Asynch_Get_Init_Timestamp(self.asynch_obj)

	def Get_Size_Global_Parameters(self):
		return self.lib.Asynch_Get_Size_Global_Parameters(self.asynch_obj)

	def Get_Global_Parameters(self):
		n = self.lib.Asynch_Get_Size_Global_Parameters(self.asynch_obj)
		c_array_type = (c_double*(n))
		arr = c_array_type()
		ret_val = self.lib.Asynch_Get_Global_Parameters(self.asynch_obj,arr)
		return [list(arr), ret_val]

	def Set_Global_Parameters(self,gparams):
		c_array_type = (c_double*len(gparams))
		arr = c_array_type()
		for i in range(len(gparams)):
			arr[i] = gparams[i]
		return self.lib.Asynch_Set_Global_Parameters(self.asynch_obj,arr,len(gparams))

	#Probably not the most efficient. This currently assumes every proc has space for every link.
	def Set_System_State(self,t_0,values):
		c_array_type = ( c_double * (len(values[0]) * len(values)) )
		arr = c_array_type()
		for i in range(len(values)):
			for j in range(len(values[i])):
				arr[i*len(values[i]) + j] = values[i][j]
		self.lib.Asynch_Set_System_State_py(self.asynch_obj,c_double(t_0),arr)

	def Reset_Peakflow_Data(self):
		self.lib.Asynch_Reset_Peakflow_Data(self.asynch_obj)

	def Set_Forcing_State(self,idx,t_0,first_file,last_file):
		return self.lib.Asynch_Set_Forcing_State(self.asynch_obj,idx,c_double(t_0),first_file,last_file)

	def Set_Temp_Files(self,set_time,set_value,output_idx):
		return self.lib.Asynch_Set_Temp_Files(self.asynch_obj,c_double(set_time),byref(set_value),output_idx)

	def Reset_Temp_Files(self,set_time):
		return self.lib.Asynch_Reset_Temp_Files(self.asynch_obj,c_double(set_time))

	def Get_Peakflow_Output_Name(self):
		peakflowname = 1024*'\0'
		ret_value = self.lib.Asynch_Get_Peakflow_Output_Name(self.asynch_obj,c_char_p(peakflowname))
		i = 0
		for i in range(0,1024):
			if peakflowname[i] == '\0':	break
		true_peakflowname = peakflowname[0:i]
		return [true_peakflowname, ret_value]

	def Set_Peakflow_Output_Name(self,peakflowname):
		return self.lib.Asynch_Set_Peakflow_Output_Name(self.asynch_obj,c_char_p(peakflowname))

	def Get_Snapshot_Output_Name(self):
		filename = 1024*'\0'
		ret_value = self.lib.Asynch_Get_Snapshot_Output_Name(self.asynch_obj,filename)
		i = 0
		for i in range(0,1024):
			if filename[i] == '\0':	break
		true_filename = filename[0:i]
		return [true_filename, ret_value]

	def Set_Snapshot_Output_Name(self,filename):
		return self.lib.Asynch_Get_Snapshot_Output_Name(asynch,filename)

	#Routines for output
	def Set_Output(self,name,data_type,func,used_states_list):
		#Note: This only allows Python functions for the output.
		if used_states_list == None:
			num_states = 0
		else:
			num_states = len(used_states_list)
		used_states = numpy.array(used_states_list)
		return self.lib.Asynch_Set_Output(self.asynch_obj,name,data_type,cast(func,ASYNCH_OUTPUT_DATATYPE),used_states.ctypes.data,num_states)

	def Check_Output(self,name):
		return self.lib.Asynch_Check_Output(self.asynch_obj,name)

	def Check_Peakflow_Output(self,name):
		return self.lib.Asynch_Check_Peakflow_Output(self.asynch_obj,name)

	def Set_Peakflow_Output(self,name,func):
		return self.lib.Asynch_Set_Peakflow_Output(self.asynch_obj,name,func)

	def Create_OutputUser_Data(self,data_size):
		return self.lib.Asynch_Create_OutputUser_Data(self.asynch_obj,data_size)

	def Free_OutputUser_Data(self):
		return self.lib.Asynch_Free_OutputUser_Data(self.asynch_obj)

	def Copy_Local_OutputUser_Data(self,location,source):
		self.lib.Asynch_Copy_Local_OutputUser_Data.restype = None
		self.lib.Asynch_Copy_Local_OutputUser_Data_py(self.asynch_obj,location,py_object(source),sys.getsizeof(source))

	def Set_Size_Local_OutputUser_Data(self,location,size):
		self.lib.Asynch_Set_Size_Local_OutputUser_Data.restype = None
		self.lib.Asynch_Set_Size_Local_OutputUser_Data(self.asynch_obj,location,size)


