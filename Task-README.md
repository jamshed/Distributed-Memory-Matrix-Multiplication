# Tasks

## Compile
----------

`mpiicc <source-file> -o <executable> -mcmodel medium -shared-intel`.

Here, set `<source-file>` to either *Task-bc.cpp*, *Task-def.cpp*, *Task-bc-2.cpp*, or *Task-def-2.cpp*. Set `<executable>` as desired.


## Run
------

Edit the file *batch-job.slurm* with required parameter settings such as the number of compute nodes, the number of processes per compute node, maximum time allocation etc. Conclude the file with a statement of the following format:

`ibrun ./<executable> -k <processing_element_count> -a <algorithm_id>`


- The `k` parameter is used to set the dimensions of the input matrices (2^k x 2^k).

- The `a` parameter value should be in {0: MM-rotate-A-rotate-B, 1: MM-rotate-A-broadcast-B, 2: MM-broadcast-A-broadcast-B}.


Submit the batch job for execution using:

`sbatch batch-job.slurm`


After running the batch job, there will be an output file *output.txt* generated at the same directory, with relevant time measurements.


For example, to run the executable *def* on 4 nodes with 16 processes/node and using the MM-broadcast-A-broadcast-B algorithm for 512 x 512 = 2^9 x 2^9 sized input matrices, set the `--nodes=4` and `--ntasks-per-node=1` at *batch-job.slurm*, with the concluding statement:

`ibrun ./def -k 9 -a 2`

The output file will be *output.txt*.
