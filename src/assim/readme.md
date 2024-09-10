# Assim


## Building PETSC on Argon

Assuming that we are using the Intel Compiler:


```
module load parallel_studio/2017.1
```


For the purpose of Asynch/Assim, the MPI version of PETSC is not required.


```
$ ./configure PETSC_ARCH=arch-linux-gnu-intel --prefix=$HOME/.local --with-cc=icc --with-fc=0 --with-blas-lapack-dir=/opt/apps/parallel_studio/2017.1/compilers_and_libraries_2017.1.132/linux/mkl/lib/intel64
$ make PETSC_DIR=/Users/sdebionne/src/petsc-3.7.6 PETSC_ARCH=arch-linux-gnu-intel all
$ make PETSC_DIR=/Dedicated/IFC/.argon PETSC_ARCH="" install
$ make PETSC_DIR=/Dedicated/IFC/.argon PETSC_ARCH="" test
``