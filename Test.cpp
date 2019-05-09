#include<cstdio>
#include<mpi.h>

int main(int argc, char *argv[])
{
    int procRank, val = 7, recvBuff;
    MPI_Status status;

    MPI_Init(&argc, &argv);

    //MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    if(procRank == 0)
    {
        MPI_Send(&val, 1, MPI_INT, 1, MPI_ANY_TAG, MPI_COMM_WORLD);
        printf("Process %d sent data %d.\n", procRank, val);
    }
    else
    {
        MPI_Recv(&recvBuff, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        printf("Process %d received data %d.\n", procRank, recvBuff);
    }

    MPI_Finalize();

    return 0;
}
