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

int M1[MAX_DIM][MAX_DIM], M2[MAX_DIM][MAX_DIM], M3[MAX_DIM][MAX_DIM];
int dim;


// Map the 2-d index (x, y) into a 1-d index, where the upper-bound for the dimensions of the 2-d indices is d.

inline int get_1D_index(int x, int y, int d = dim)
{
    return x * d + y;
}

// Map the 1-d index idx into a 2-d index (x, y), where the upper-bound for the dimensions of the 2-d indices is d.

inline void get_2D_indices(int idx, int &x, int &y, int d = dim)
{
    x = idx / d;
    y = idx % d;
}

// Get the input matrices dimension n from command line arguments.

void parse_input(int argc, char *argv[], int &n)
{
    for(int i = 0; i < argc; ++i)
        if(!strcmp(argv[i], "-n"))
        {
            n = atoi(argv[i + 1]);
            i++;
        }
}


// Generate a random square matrix of dimension n in A, each entry within (-UB, +UB).

void get_random_square_matrix(int A[MAX_DIM][MAX_DIM], int n)
{
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
        {
            A[i][j] = rand() % UB;
            if(rand() % 2)
                A[i][j] = -A[i][j];
        }
}

// Print the n-dimensional square matrix A, preceded by a message.

void print_matrix(const char *message, int A[MAX_DIM][MAX_DIM], int n)
{
    printf("\n%s\n", message);

    for(int i = 0; i < n; ++i)
    {
        printf("\n");

        for(int j = 0; j < n; ++j)
            printf("%d\t", A[i][j]);
    }

    printf("\n");
}

void print_array(const char *message, int *A, int n)
{
    printf("\n%s\n", message);

    for(int i = 0; i < n; ++i)
        printf("%d\t", A[i]);

    printf("\n");
}


// Initial distribution of the n / sqrt(p) - dimensional blocks for the n-dimensional matrix M to p processes;
// This process has rank 'processorRank', and will receive its corresponding block into a buffer 'recvBuf'.

void distribute_blocks(int M[MAX_DIM][MAX_DIM], int *recvBuf, int n, int p, int processRank)
{
    int *sendBuf = new int[n * n];
    int blockDim, blockSize, idx = 0;

    dim = sqrt(p);
    blockDim = n / dim;
    blockSize = blockDim * blockDim;

    if(processRank == MASTER)
        for(int i = 0; i < dim; ++i)
            for(int j = 0; j < dim; ++j)
                for(int k = 0; k < blockDim; ++k)
                    for(int l = 0; l < blockDim; ++l)
                        sendBuf[idx++] = M[i * blockDim + k][j * blockDim + l];

    MPI_Scatter(sendBuf, blockSize, MPI_INT, recvBuf, blockSize, MPI_INT, MASTER, MPI_COMM_WORLD);

    delete sendBuf;
}


// Shifts the data stored in:
//      - buffer 'a' for this process to buffer 'a' of the process 'leftRecvr';
//      - buffer 'b' for this process to buffer 'b' of the process 'upRecvr'.
// Also receives the data stored in:
//      - buffer 'a' of the process 'rightSendr' to buffer 'a' of this process;
//      - buffer 'b' of the process 'downSendr' to buffer 'b' of this process.
// Use the buffers 'temp_a' and 'temp_b' as temporary storage to avoid data overlap.

inline void shift_blocks(int *a, int *b, int *temp_a, int *temp_b, int blockSize,
                         int leftRecvr, int rightSendr, int upRecvr, int downSendr)
{
    MPI_Request sendReq[2], recvReq[2];

    MPI_Isend(a, blockSize, MPI_INT, leftRecvr, 0, MPI_COMM_WORLD, &sendReq[0]),
    MPI_Irecv(temp_a, blockSize, MPI_INT, rightSendr, MPI_ANY_TAG, MPI_COMM_WORLD, &recvReq[0]);

    MPI_Isend(b, blockSize, MPI_INT, upRecvr, 0, MPI_COMM_WORLD, &sendReq[1]),
    MPI_Irecv(temp_b, blockSize, MPI_INT, downSendr, MPI_ANY_TAG, MPI_COMM_WORLD, &recvReq[1]);


    MPI_Wait(&sendReq[0], MPI_STATUS_IGNORE), MPI_Wait(&recvReq[0], MPI_STATUS_IGNORE);
    memcpy(a, temp_a, blockSize * sizeof(int)),

    MPI_Wait(&sendReq[1], MPI_STATUS_IGNORE), MPI_Wait(&recvReq[1], MPI_STATUS_IGNORE);
    memcpy(b, temp_b, blockSize * sizeof(int));
}


// Initial block alignment for the n-dimensional input matrices into the buffers a and b for p processes.

void align_initially(int *a, int *b, int *c, int n, int p, int processRank)
{
    dim = sqrt(p);
    int blockSize = n * n / p, processRow, processCol, leftRecvr, upRecvr, rightSendr, downSendr;
    int *temp_a = new int[blockSize], *temp_b = new int[blockSize];

    memset(c, 0, blockSize * sizeof(int));


    get_2D_indices(processRank, processRow, processCol);

    leftRecvr = get_1D_index(processRow, (processCol - processRow + dim) % dim),
    rightSendr = get_1D_index(processRow, (processCol + processRow) % dim);

    upRecvr = get_1D_index((processRow - processCol + dim) % dim, processCol);
    downSendr = get_1D_index((processRow + processCol) % dim, processCol);


    shift_blocks(a, b, temp_a, temp_b, blockSize, leftRecvr, rightSendr, upRecvr, downSendr);

    delete temp_a, delete temp_b;
}


void multiply_matrices(int *A, int *B, int *C, int n)
{
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            for(int k = 0; k < n; ++k)
                C[get_1D_index(i, j, n)] += A[get_1D_index(i, k, n)] * B[get_1D_index(k, j, n)];
}

void iterate_and_multiply(int *a, int *b, int *c, int n, int p, int processRank)
{
    dim = sqrt(p);
    int blockDim = n / dim, blockSize = blockDim * blockDim;
    int *temp_a = new int[blockSize], *temp_b = new int[blockSize];

    for(int i = 1; i <= dim; ++i)
    {
        multiply_matrices(a, b, c, blockDim);

        if(i < dim)
        {
            int processRow, processCol;

            get_2D_indices(processRank, processRow, processCol);

            int leftRecvr = get_1D_index(processRow, (processCol - 1 + dim) % dim),
                rightSendr = get_1D_index(processRow, (processCol + 1) % dim),
                upRecvr = get_1D_index((processRow - 1 + dim) % dim, processCol),
                downSendr = get_1D_index((processRow + 1) % dim, processCol);

            shift_blocks(a, b, temp_a, temp_b, blockSize, leftRecvr, rightSendr, upRecvr, downSendr);
        }
    }

    delete temp_a, delete temp_b;
}

void collect_blocks(int M[MAX_DIM][MAX_DIM], int *sendBuf, int n, int p, int processRank)
{
    int *recvBuf = new int[n * n];
    int blockDim, blockSize, idx = 0;

    dim = sqrt(p);
    blockDim = n / dim;
    blockSize = blockDim * blockDim;

    //if(processRank == MASTER)
    //    printf("P_dim = %d, Block_dim = %d, Block_size = %d\n", dim, blockDim, blockSize);

    MPI_Gather(sendBuf, blockSize, MPI_INT, recvBuf, blockSize, MPI_INT, MASTER, MPI_COMM_WORLD);

    if(processRank == MASTER)
        for(int i = 0; i < dim; ++i)
            for(int j = 0; j < dim; ++j)
                for(int k = 0; k < blockDim; ++k)
                    for(int l = 0; l < blockDim; ++l)
                        M[i * blockDim + k][j * blockDim + l] = recvBuf[idx++];

    //if(processRank == MASTER)
    //    printf("Collected %d items.\n", idx);

    delete recvBuf;
}

bool verify_correctness(int A[MAX_DIM][MAX_DIM], int B[MAX_DIM][MAX_DIM], int C[MAX_DIM][MAX_DIM],
                        int R[MAX_DIM][MAX_DIM], int n)
{
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
        {
            C[i][j] = 0;
            for(int k = 0; k < n; ++k)
                C[i][j] += A[i][k] * B[k][j];

            if(R[i][j] != C[i][j])
                return false;
        }

    /*for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            if(C[i][j] != R[i][j])
                return false;
*/
    return true;
}

void MM_rotate_A_rotate_B(int n, int p, int processRank)
{
    int blockSize = n * n / p;
    int *a = new int[blockSize], *b = new int[blockSize], *c = new int[blockSize];
    dim = sqrt(p);

    distribute_blocks(M1, a, n, p, processRank);
    distribute_blocks(M2, b, n, p, processRank);

    align_initially(a, b, c, n, p, processRank);

    iterate_and_multiply(a, b, c, n, p, processRank);

    collect_blocks(M3, c, n, p, processRank);

    delete a, delete b, delete c;
}

void iterate_and_multiply(int *a, int *b, int *c, int n, int p, int processorRank, MPI_Comm columnGroup)
{

}

void MM_rotate_A_broadcast_B(int n, int p, int processRank)
{
    int blockSize = n * n / p, processRow, processCol;
    int *a = new int[blockSize], *b = new int[blockSize], *c = new int[blockSize];
    MPI_Comm columnGroup;
    dim = sqrt(p);

    distribute_blocks(M1, a, n, p, processRank);
    distribute_blocks(M2, b, n, p, processRank);

    memset(c, 0, blockSize * sizeof(int));

    get_2D_indices(processRank, processRow, processCol);
    MPI_Comm_split(MPI_COMM_WORLD, processCol, processRow, &columnGroup);



    delete a, delete b, delete c;
}

int main(int argc, char *argv[])
{
    int n, p, processRank;

    srand(time(NULL));

    parse_input(argc, argv, n);

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);

    dim = sqrt(p);

    if(processRank == MASTER)
        get_random_square_matrix(M1, n), get_random_square_matrix(M2, n);

    MM_rotate_A_rotate_B(n, p, processRank);

    /*
    int *a = new int[n * n  / p], *b = new int[n * n / p], *c = new int[n * n / p];

    distribute_blocks(M1, a, n, p, processRank);
    distribute_blocks(M2, b, n, p, processRank);


    align_initially(a, b, c, n, p, processRank);

    //puts("initial distribution done.");

    iterate_and_multiply(a, b, c, n, p, processRank);

    //printf("Iteration of matrix multiplications done for process %d.\n", processRank);

    collect_blocks(M3, c, n, p, processRank);

    delete a, delete b, delete c;
    */

    if(processRank == MASTER)
    {
        //print_matrix("\nA", M1, n);
        //print_matrix("\nB", M2, n);

        int C[MAX_DIM][MAX_DIM];
        printf("%s matrix multiplication for n = %d.\n", verify_correctness(M1, M2, C, M3, n) ? "Correct" : "Incorrect", n);

        //print_matrix("\nCorrect result:", C, n);
        //print_matrix("\n Result from parallel computation:", M3, n);
    }

    MPI_Finalize();

    return 0;
}
