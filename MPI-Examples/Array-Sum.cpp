#include<cstdio>
#include<ctime>
#include<cstdlib>
#include<mpi.h>

using namespace std;

#define UB 10

void generate_random_array(int *A, int n)
{
    for(int i = 0; i < n; ++i)
        A[i] = rand() % UB;
}

void print_array(int *A, int n)
{
    for(int i = 0; i < n; ++i)
        printf("%d ", A[i]);

    puts("\n");
}

int main(int argc, char *argv[])
{
    int processRank, commSize, n, distCount, *a, *b, *c, *d, *e, *f;

    srand(time(NULL));

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);

    n = atoi(argv[1]);
    printf("Rank = %d, n  = %d\n", processRank, n);

    distCount = n / commSize;

    d = new int[distCount];
    e = new int[distCount];
    f = new int[distCount];

    if(processRank == 0)
    {
        a = new int[n];
        b = new int[n];
        c = new int[n];

        generate_random_array(a, n), generate_random_array(b, n);

        print_array(a, n);
        print_array(b, n);
    }

    MPI_Scatter(a, distCount, MPI_INT, d, distCount, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(b, distCount, MPI_INT, e, distCount, MPI_INT, 0, MPI_COMM_WORLD);

    for(int i = 0; i < distCount; ++i)
        f[i] = d[i] + e[i];

    MPI_Gather(f, distCount, MPI_INT, c, distCount, MPI_INT, 0, MPI_COMM_WORLD);

    if(processRank == 0)
        print_array(c, n);

    MPI_Finalize();

    return 0;
}
