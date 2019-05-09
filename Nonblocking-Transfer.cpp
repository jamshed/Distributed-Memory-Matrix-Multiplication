#include<cstdio>
#include<mpi.h>

using namespace std;

int main(int argc, char *argv[])
{
    int processRank, data = 7;
    MPI_Status status;
    MPI_Request request;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);

    if(processRank == 0)
    {
        MPI_Isend(&data, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &request);

        printf("Nonblocking send of data %d from process %d initiated.\n", data, processRank);

        // tasks that do not modify data

        MPI_Wait(&request, &status);

        printf("Nonblocking send at process %d completed.\n", processRank);
    }
    else if(processRank == 1)
    {
        int recvBuff;
        MPI_Irecv(&recvBuff, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &request);

        printf("Nonblocking receive at process %d initiated.\n", processRank);

        // tasks that do not modify recvBuff

        MPI_Wait(&request, &status);

        printf("Nonblocking receive of data %d at process %d completed.\n", recvBuff, processRank);
    }

    return 0;
}
