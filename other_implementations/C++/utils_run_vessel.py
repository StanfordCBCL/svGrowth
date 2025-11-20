from mpi4py import MPI
import numpy as np
import pickle
import os
from ctypes import cdll
from ctypes import c_char_p
from ctypes import c_int
from ctypes import c_double
lib = cdll.LoadLibrary('./libgnr.so')
run = lib.run
run.argtypes = [c_char_p, c_char_p, c_int, c_int, c_int, c_int, c_double,c_double,c_double, c_double, c_double, c_double,c_double,c_double,c_double,c_double,c_double,c_double,c_double,c_double]
#run(dir_prefix.encode('ascii'), filename.encode('ascii'),restart_arg,iter_arg,gnr_arg,num_days,step_size,sigma_arg,tauw_arg,anysm_arg,tevg_arg,F0, F1, F2, F3, F4,F5, F6, F7, F8)


comm = MPI.COMM_WORLD
size = comm.Get_size() # new: gives number of ranks in comm
rank = comm.Get_rank()

input_file = open("input_array.dat","rb")
input_array = pickle.load(input_file)
input_file.close()

numData = len(input_array)
numDataPerRank = numData//size
numRemainder = numData%size

for i in range(numDataPerRank*rank, numDataPerRank*(rank+1)):
    inputData = input_array[i]
    prefix = inputData[2].encode('ascii')
    suffix = ('python_'+str(inputData[0])+'_'+str(inputData[1])).encode('ascii')
    run(prefix, suffix, inputData[3], inputData[4], inputData[5],inputData[6], inputData[7], inputData[8], inputData[9], inputData[10], inputData[11], inputData[12], inputData[13],inputData[14],inputData[15],inputData[16],inputData[17],inputData[18], inputData[19], inputData[20])

if rank < numRemainder:
    i = numDataPerRank*size + rank
    inputData = input_array[i]
    prefix = inputData[2].encode('ascii')
    suffix = ('python_'+str(inputData[0])+'_'+str(inputData[1])).encode('ascii')
    run(prefix, suffix, inputData[3], inputData[4], inputData[5],inputData[6], inputData[7], inputData[8], inputData[9], inputData[10], inputData[11], inputData[12], inputData[13],inputData[14],inputData[15],inputData[16],inputData[17],inputData[18], inputData[19], inputData[20])
