#include<cstdio>
#include<mpi.h>

int main(int argc, char *argv[])
{
    int procRank, val = 7, recvBuff;
    MPI_Status status;
    MPI_Request req[2];

    MPI_Init(&argc, &argv);

    //MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    if(procRank == 0)
    {
        MPI_Isend(&val, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &req[0]);
        printf("Process %d sent data %d.\n", procRank, val);
    }
    else if(procRank == 1)
    {
        MPI_Irecv(&recvBuff, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &req[1]);

        for(int i = 0; i < 1000000; ++i);

        MPI_Wait(&req[1], &status);
        printf("Process %d received data %d.\n", procRank, recvBuff);
    }
    else
        printf("Nothing to do for process %d.\n", procRank);

    MPI_Finalize();

    return 0;
}
