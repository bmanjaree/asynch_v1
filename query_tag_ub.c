#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int tag_ub;
    int flag;
    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub, &flag);

    if (flag) {
        printf("MPI_TAG_UB = %d\n", tag_ub);
    } else {
        printf("MPI_TAG_UB attribute not found\n");
    }

    MPI_Finalize();
    return 0;
}
