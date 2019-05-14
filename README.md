# Distributed-Memory-Matrix-Multiplication

## Compile
----------

To work with input matrices of sizes upto 2^10 x 2^10:

`mpiicc Distributed-Memory-Matrix-Multiplication.cpp -o matmul`

To work with input matrix of sizes more than 2^10 x 2^10:

`mpiicc Distributed-Memory-Matrix-Multiplication.cpp -o matmul -mcmodel medium -shared-intel`.

Note that at the *Distributed-Memory-Matrix-Multiplication.cpp*, you need to modify the macro `MAX_DIM` to suit your input matrices dimension. Specifically, modify your dimension at line 12, `#define MAX_DIM 1024`.
Still, there might be segmentation faults present while running at local machines. Please use computing clusters while working with such dimensional matrices.


## Run
------

`mpiexec -np <processing_element_count> ./matmul -k <matrix_dimension_parameter> -a <algorithm_id>`

- The `k` parameter is used to set the dimensions of the input matrices (2^k x 2^k).

- The `a` parameter value should be in {0: MM-rotate-A-rotate-B, 1: MM-rotate-A-broadcast-B, 2: MM-broadcast-A-broadcast-B}.

- After running the executable *matmul*, there will be an output file *output.txt* generated at the same directory, with relevant informations (correctness here) depending on the source files.


For example, to run the executable *matmul* on 4 processes and using the MM-broadcast-A-broadcast-B algorithm for 16 x 16 = 2^4 x 2^4 sized input matrices, execute the following.

`mpiexec -np 4 ./matmul -k 4 -a 2`

The output file will be *output.txt*.
