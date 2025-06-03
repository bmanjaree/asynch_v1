#Use to compile
#!/bin/bash
# export C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/local/hdf5/gcc/1.14.4/include
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/hdf5/gcc/1.14.4/lib64

module purge
module load intel-mpi/gcc/2021.15
module load hdf5/gcc/1.14.4        
module load netcdf/gcc/hdf5-1.14.4/4.9.2


autoreconf --install
mkdir build
cd build
make clean
../configure CFLAGS="-g -O0 -Wno-format-security -DDEBUG" --without-postgresql
make
