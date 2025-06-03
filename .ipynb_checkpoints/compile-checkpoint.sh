#Use to compile
#!/bin/bash
# export C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/local/hdf5/gcc/1.14.4/include
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/hdf5/gcc/1.14.4/lib64
autoreconf --install
mkdir build
cd build
make clean
../configure CFLAGS="-g -O0 -Wno-format-security" --without-postgresql
make
