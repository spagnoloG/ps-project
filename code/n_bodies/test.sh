#!/bin/sh

set -x -e -u -o pipefail

VERSIONS_DEFAULT=(0 1)
VERSIONS_PTHREAD=(0 1 2)

N_BODIES=(100 1000 10000)
N_STEPS=(100 1000 10000)
N_THREADS=(1 2 4 8 16 32)
N_NODES_MPI=(1 2)

# Generate CSV
echo "filename, n_bodies, n_steps, n_threads, version, time, n_nodes" > results.csv


## SEQUENTIAL ##
N_BODIES=(100 1000)
# Loop through STEPS
for STEPS in "${N_STEPS[@]}"; do
    # Loop through BODIES
    for BODIES in "${N_BODIES[@]}"; do
        # Loop through VERSIONS
        for VERSION in "${VERSIONS_DEFAULT[@]}"; do
            time=$(make run_default_remote n_bodies=$BODIES n_steps=$STEPS version=$VERSION | grep CPU | awk '{ print $2 }' )
            echo "sequential, $BODIES, $STEPS, 1, $VERSION, $time, 1" >> results.csv
        done
    done
done

## PTHREAD ##
N_BODIES=(100 1000 10000)
# Loop through STEPS
for STEPS in "${N_STEPS[@]}"; do
    # Loop through BODIES
    for BODIES in "${N_BODIES[@]}"; do
        # Loop through VERSIONS
        for VERSION in "${VERSIONS_PTHREAD[@]}"; do
            # Loop through THREADS
            for THREADS in "${N_THREADS[@]}"; do
                time=$(make run_pthread_remote n_bodies=$BODIES n_steps=$STEPS version=$VERSION n_threads=$THREADS | grep CPU | awk '{ print $2 }')
                echo "pthread, $BODIES, $STEPS, $THREADS, $VERSION, $time, 1" >> results.csv
            done
        done
    done
done

## OpenMP ##

# Loop through STEPS
for STEPS in "${N_STEPS[@]}"; do
    # Loop through BODIES
    for BODIES in "${N_BODIES[@]}"; do
        # Loop through THREADS
        for THREADS in "${N_THREADS[@]}"; do
            time=$(make run_openmp_remote n_bodies=$BODIES n_steps=$STEPS n_threads=$THREADS | grep CPU | awk '{ print $2 }')
            echo "openmp, $BODIES, $STEPS, $THREADS, 0, $time, 1" >> results.csv
        done
    done
done

