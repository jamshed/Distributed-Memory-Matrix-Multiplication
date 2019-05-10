#include<cstdio>
#include<mpi.h>

#define SIZE 20
#define RECEIVER 0

int main(int argc, char *argv[])
{
    char A[3][SIZE] = {"Stony", "Brook", "University"}, B[3][SIZE];
    int procRank, commSize;
    MPI_Request sendReq[3], recvReq[3];

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    if(procRank != RECEIVER)
    {
        int memIdx = procRank - (procRank > RECEIVER);

        printf("Process %d initiated sending message %s to process %d.\n", procRank, &A[memIdx][0], RECEIVER);
        MPI_Isend(&A[memIdx][0], SIZE, MPI_CHAR, RECEIVER, 0, MPI_COMM_WORLD, &sendReq[memIdx]);
    }
    else
    {
        for(int i = 0; i < 4; ++i)
            if(i != RECEIVER)
            {
                int memIdx = i - (i > RECEIVER);

                MPI_Irecv(&B[memIdx][0], SIZE, MPI_CHAR, i, MPI_ANY_TAG, MPI_COMM_WORLD, &recvReq[memIdx]);
                printf("Process %d initiated receiving message from process %d.\n", RECEIVER, i);
            }

        for(int i = 0; i < 4; ++i)
            if(i != RECEIVER)
            {
                int memIdx = i - (i > RECEIVER);

                MPI_Wait(&recvReq[memIdx], MPI_STATUS_IGNORE);
                printf("Process %d received message %s from process %d.\n", RECEIVER, &B[memIdx][0], i);
            }
    }

    MPI_Finalize();

    return 0;
}
