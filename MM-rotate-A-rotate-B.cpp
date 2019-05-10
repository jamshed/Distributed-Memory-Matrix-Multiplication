#include<cstdio>
#include<cstring>
#include<cmath>
#include<cstdlib>
#include<mpi.h>

using namespace std;

#define MAX_DIM 1025
#define P 4
#define UB 10
#define MASTER 0

int M1[MAX_DIM][MAX_DIM], M2[MAX_DIM][MAX_DIM];

int dim;

inline int get_1D_index(int x, int y)
{
    return x * dim + y;
}

inline void get_2D_indices(int idx, int &x, int &y)
{
    x = idx / dim;
    y = idx % dim;
}

void parse_input(int argc, char *argv[], int &n)
{
    for(int i = 0; i < argc; ++i)
        if(!strcmp(argv[i], "-n"))
        {
            n = atoi(argv[i + 1]);
            i++;
        }
}


void get_random_square_matrix(int A[MAX_DIM][MAX_DIM], int n)
{
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            A[i][j] = rand() % UB;
}

void print_matrix(const char *message, int A[MAX_DIM][MAX_DIM], int n)
{
    puts(message);

    for(int i = 0; i < n; ++i)
    {
        printf("\n");

        for(int j = 0; j < n; ++j)
            printf("%d ", A[i][j]);
    }

    printf("\n");
}

void distribute(int A[MAX_DIM][MAX_DIM], int B[MAX_DIM][MAX_DIM], int *a, int *b, int n, int p)
{
    int processRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);

    if(processRank == MASTER)
    {
        int **sndBuf = new int *[p];
        int blockDim = n / dim;
        MPI_Request *sndReq = new MPI_Request[p];

        dim = sqrt(p);
        for(int i = 0; i < dim; ++i)
            for(int j = 0; j < dim; ++j)
            {
                int receiver = get_1D_index(i, j);

                sndBuf[receiver] = new int[n * n / p];

                int idx = 0;
                for(int k = 0; k < blockDim; ++k)
                    for(int l = 0; l < blockDim; ++l)
                        sndBuf[receiver][idx++] = A[i * blockDim + k][j * blockDim + l];

                MPI_Isend(sndBuf[receiver], n * n / p, MPI_INT, receiver, 0, MPI_COMM_WORLD, &sndReq[receiver]);
            }
    }
}

void distribute(int M[MAX_DIM][MAX_DIM], int *rcvBuf, int n, int p)
{
    int *sndBuf = new int[n * n];
    int blockDim, processRank, idx = 0;

    dim = sqrt(p);
    blockDim = n / dim;

    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);

    if(processRank == MASTER)
        for(int i = 0; i < dim; ++i)
            for(int j = 0; j < dim; ++j)
                for(int k = 0; k < blockDim; ++k)
                    for(int l = 0; l < blockDim; ++l)
                        sndBuf[idx++] = M[i * blockDim + k][j * blockDim + l];

    MPI_Scatter(sndBuf, n * n / p, MPI_INT, rcvBuf, n * n / p, MPI_INT, MASTER, MPI_COMM_WORLD);

    delete sndBuf;
}

int main(int argc, char *argv[])
{
    int n, p, processRank;

    srand(time(NULL));

    parse_input(argc, argv, n);

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);

    // printf("Process count = %d, rank = %d\n", p, processRank);

    if(processRank == MASTER)
    {
        get_random_square_matrix(M1, n);
        get_random_square_matrix(M2, n);

        print_matrix("\nA", M1, n);
        //print_matrix("\nB", M2, n);

        dim = sqrt(p);
    }

    int *a = new int[n * n  / p];

    distribute(M1, a, n, p);

    if(processRank == MASTER)
        for(int i = 0; i < n * n / p; ++i)
            printf("%d\t", a[i]);

    delete a;
    MPI_Finalize();

    return 0;
}
