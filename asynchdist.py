from asynch_py.asynch_interface import asynchsolver
import sys

#Parse command line arguments
numargs = len(sys.argv)
if numargs != 2:
	print 'Need an input .gbl file'
	sys.exit(1)

#Prepare system
asynch = asynchsolver()
print "Reading global file..."
asynch.Parse_GBL(sys.argv[1])
print "Loading network..."
asynch.Load_Network()
print "Partitioning network..."
asynch.Partition_Network()
print "Loading parameters..."
asynch.Load_Network_Parameters(False)
print "Reading dam and reservoir data..."
asynch.Load_Dams()
print "Setting up numerical error data..."
asynch.Load_Numerical_Error_Data()
print "Initializing model..."
asynch.Initialize_Model()
print "Loading initial conditions..."
asynch.Load_Initial_Conditions()
print "Loading forcings..."
asynch.Load_Forcings()
print "Loading output data information..."
asynch.Load_Save_Lists()
print "Finalizing network..."
asynch.Finalize_Network()
print "Calculating initial step sizes..."
asynch.Calculate_Step_Sizes()

N = asynch.Get_Number_Links()
print 'I see',N,'links.'


#Prepare outputs
asynch.Prepare_Temp_Files()
asynch.Write_Current_Step()
asynch.Prepare_Peakflow_Output()
asynch.Prepare_Output()

#Advance solver
asynch.Advance(True)
print 'Calculations done!'

#Take a snapshot
asynch.Take_System_Snapshot(None)

#Create output files
asynch.Create_Output(None)
asynch.Create_Peakflows_Output()

