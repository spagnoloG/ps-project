#!/bin/sh

set -x -e -u -o pipefail

VERSIONS_DEFAULT=(0 1)
VERSIONS_PTHREAD=(0 1 2)
VERSIONS_MPI=(0 1)
VERSIONS_MP=(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14)

N_STEPS=(100 1000)
N_THREADS=(8 16 24)

# Generate CSV
#echo "filename, n_bodies, n_steps, n_threads, version, time, n_procs" > results.csv
#
#
### SEQUENTIAL ##
#N_BODIES=(10000)
### Loop through STEPS
#for STEPS in "${N_STEPS[@]}"; do
#    # Loop through BODIES
#    for BODIES in "${N_BODIES[@]}"; do
#        # Loop through VERSIONS
#        for VERSION in "${VERSIONS_DEFAULT[@]}"; do
#            time=$(make run_default_remote n_bodies=$BODIES n_steps=$STEPS version=$VERSION | grep CPU | awk '{ print $2 }' )
#            echo "sequential, $BODIES, $STEPS, 1, $VERSION, $time, 1" >> results.csv
#        done
#    done
#done
#
### PTHREAD ##
N_BODIES=(100 1000 10000)
## Loop through STEPS
#N_THREADS=(8 16 24)
#for STEPS in "${N_STEPS[@]}"; do
#    # Loop through BODIES
#    for BODIES in "${N_BODIES[@]}"; do
#        # Loop through VERSIONS
#        for VERSION in "${VERSIONS_PTHREAD[@]}"; do
#            # Loop through THREADS
#            for THREADS in "${N_THREADS[@]}"; do
#                #if [ $BODIES -gt 1000 ] && [ $THREADS -lt 4 ]; then
#                #    continue
#                #fi
#                #if [ $BODIES -gt 10000 ] && [ $THREADS -lt 8 ]; then
#                #    continue
#                #fi
#                #if [ $STEPS -gt 10000 ] && [ $THREADS -lt 16 ]; then
#                #    continue
#                #fi 
#                time=$(make run_pthread_remote n_bodies=$BODIES n_steps=$STEPS version=$VERSION n_threads=$THREADS | grep CPU | awk '{ print $2 }')
#                echo "pthread, $BODIES, $STEPS, $THREADS, $VERSION, $time, 1" >> results.csv
#            done
#        done
#    done
#done

## OpenMP ##
# Loop through STEPS
#for STEPS in "${N_STEPS[@]}"; do
#    # Loop through BODIES
#    for BODIES in "${N_BODIES[@]}"; do
#        # Loop through VERSIONS
#        for VERSION in "${VERSIONS_MP[@]}"; do
#            # Loop through THREADS
#            for THREADS in "${N_THREADS[@]}"; do
#                time=$(make run_mp_remote n_bodies=$BODIES n_steps=$STEPS n_threads=$THREADS version=$VERSION | grep CPU | awk '{ print $2 }')
#                echo "openmp, $BODIES, $STEPS, $THREADS, $VERSION, $time, 1" >> results.csv
#            done
#        done
#    done
#done
#
### MPI ##
N_PROCS=(1 2 4 8 16 24)
N_THREADS=(1 2 4 8 16 24)
N_STEPS=(100)
N_BODIES=(10000)
# Loop through N_STEPS
for STEPS in "${N_STEPS[@]}"; do
    # Loop through N_BODIES
    for BODIES in "${N_BODIES[@]}"; do
        # Loop through N_PROCS
        for PROC in "${N_PROCS[@]}"; do
            # Loop through N_NODEs
             for THREADS in "${N_THREADS[@]}"; do
                # Check if sum of threads and procs is equal to 24
                if [ $((PROC + (PROC * THREADS))) -gt 24 ]; then
                    continue
                fi
                # Loop through THREADS
                for VERSION in "${VERSIONS_MPI[@]}"; do
                    time=$(make run_mpi_remote n_bodies=$BODIES n_steps=$STEPS version=$VERSION n_proc=$PROC n_threads=$THREADS | grep CPU | grep CPU | cut --delimiter=":" --fields=2 | head -n 1 | cut --delimiter="s" --fields=1)
                    echo "mpi, $BODIES, $STEPS, $THREADS, $VERSION, $time, $PROC" >> results.csv
                done
            done
        done
    done
done
