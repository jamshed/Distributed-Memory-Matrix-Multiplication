#include<cstdio>
#include<mpi.h>

using namespace std;

int main(int argc, char *argv[])
{
    int processRank, data = 77, count;
    MPI_Status status;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);

    if(processRank == 0)
    {
        MPI_Send(&data, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        printf("Process %d sent data %d.\n", processRank, data);
    }
    else
    {
        MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_INT, &count);

        int *recvBuff = new int[count];
        MPI_Recv(recvBuff, count, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        printf("Process %d received data %d.\n", processRank, *recvBuff);
    }

    MPI_Finalize();

    return 0;
}
