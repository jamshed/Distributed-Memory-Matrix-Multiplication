#include<cstdio>
#include<mpi.h>

int main(int argc, char *argv[])
{
    //printf("Hello MPI world!\n");

    int p, procRank;

    MPI_Init(&argc, &argv);

    printf("Hello MPI world!\n");

    MPI_Finalize();

    return 0;
}
