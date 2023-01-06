CC=gcc 
CC_MPI=mpicc
CFLAGS=-g -Wall -lm  -ffast-math -O3 -march=native
MPI_FLAGS=-fopenmp
PTHREAD_CFLAGS=-pthread
MAIN=n_bodies
MAIN_PTHREAD=n_bodies_pthread
MAIN_MP=n_bodies_openmp
MAIN_MPI=n_bodies_mpi
USERNAME=$(shell ssh nsc "whoami")
DESTINATION=/ceph/grid/home/$(USERNAME)/n_bodies

all: $(MAIN) $(MAIN_PTHREAD)

$(MAIN): $(MAIN).c
	$(CC) $(CFLAGS) $(MAIN).c -o $(MAIN).o

$(MAIN_PTHREAD): $(MAIN_PTHREAD).c
	$(CC) $(CFLAGS) $(PTHREAD_CFLAGS) $(MAIN_PTHREAD).c -o $(MAIN_PTHREAD).o

compile_default:
	@test -n "$(version)" || (echo "Usage: make compile_default version=0"; exit 1)

	@echo "Compiling default version $(version)"

	@$(CC) $(CFLAGS) $(MAIN)$(version).c -o $(MAIN)$(version).o

compile_pthread:
	@test -n "$(version)" || (echo "Usage: make compile_pthread version=0"; exit 1)
	@echo "Compiling pthread version $(version)"

	@$(CC) $(CFLAGS) $(PTHREAD_CFLAGS) $(MAIN_PTHREAD)$(version).c -o $(MAIN_PTHREAD)$(version).o

compile_mpi:
	@test -n "$(n_threads)" || (echo "Please specify n_threads" && exit 1)
	$(CC_MPI) $(CFLAGS) $(MPI_FLAGS) -DOMP_NUM_THREADS=$(n_threads) $(MAIN_MPI).c -o $(MAIN_MPI).o

compile_mp:
	@test -n "$(n_threads)" || (echo "Please specify n_threads" && exit 1)
	$(CC) $(CFLAGS) $(CFLAGS) $(MPI_FLAGS) -DOMP_NUM_THREADS=$(n_threads) $(MAIN_MP).c -o $(MAIN_MP).o

clean: 
	rm *.o || true
	rm *.txt || true

run_local_main:
	@test -n "$(n_bodies)" || (echo "n_bodies is not set"; exit 1)
	@test -n "$(n_steps)" || (echo "n_steps is not set"; exit 1)

	./$(MAIN).o $(n_bodies) $(n_steps)

run_local_main_pthread:
	@test -n "$(n_bodies)" || (echo "n_bodies is not set"; exit 1)
	@test -n "$(n_steps)" || (echo "n_steps is not set"; exit 1)
	@test -n "$(n_threads)" || (echo "n_threads is not set"; exit 1)

	./$(MAIN_PTHREAD).o $(n_bodies) $(n_steps) $(n_threads)

run_local_main_mp:
	@test -n "$(n_bodies)" || (echo "n_bodies is not set"; exit 1)
	@test -n "$(n_steps)" || (echo "n_steps is not set"; exit 1)

	./$(MAIN_MP).o $(n_bodies) $(n_steps)

run_default_remote:
	@test -n "$(n_bodies)" || (echo "n_bodies is not set"; exit 1)
	@test -n "$(n_steps)" || (echo "n_steps is not set"; exit 1)
	@test -n "$(version)" || (echo "version is not set"; exit 1)

	ssh nsc "mkdir -p $(DESTINATION)"

	@echo "Copying files to nsc"
	rsync -a --progress ./Makefile ./n_bodies$(version).c nsc:$(DESTINATION)

	@echo "Compiling on nsc"
	ssh nsc "cd $(DESTINATION) && make clean compile_default version=$(version)"

	@echo "Running on nsc"

	ssh nsc "cd $(DESTINATION) && srun --reservation=fri ./$(MAIN)$(version).o $(n_bodies) $(n_steps)"
	
	@echo "Copying files from nsc"
	rsync -a --progress nsc:$(DESTINATION)/n_bodies$(version).txt ./

run_mp_remote:
	@test -n "$(n_bodies)" || (echo "n_bodies is not set"; exit 1)
	@test -n "$(n_steps)" || (echo "n_steps is not set"; exit 1)
	@test -n "$(n_threads)" || (echo "n_threads is not set"; exit 1)

	ssh nsc "mkdir -p $(DESTINATION)"
	
	@echo "Copying files to nsc"
	rsync -a --progress ./Makefile ./n_bodies_openmp.c nsc:$(DESTINATION)

	@echo "Compiling on nsc"
	ssh nsc "cd $(DESTINATION) && make clean compile_mp n_threads=$(n_threads)"

	@echo "Running on nsc"
	ssh nsc "cd $(DESTINATION) && srun --reservation=fri --cpus-per-task=$(n_threads) --mem-per-cpu=10G ./$(MAIN_MP).o $(n_bodies) $(n_steps)"
	
	@echo "Copying the results back"
	rsync -a --progress nsc:"$(DESTINATION)/output_mp.txt" ./

run_pthread_remote:
	@test -n "$(n_bodies)" || (echo "n_bodies is not set"; exit 1)
	@test -n "$(n_threads)" || (echo "n_threads is not set"; exit 1)
	@test -n "$(n_steps)" || (echo "n_steps is not set"; exit 1)
	@test -n "$(version)" || (echo "version is not set"; exit 1)
	
	ssh nsc "mkdir -p $(DESTINATION)"
	
	@echo "Copying files to nsc"
	rsync -a --progress ./Makefile ./n_bodies_pthread$(version).c  nsc:"$(DESTINATION)"

	@echo "Compiling and running on nsc"
	ssh nsc "cd $(DESTINATION) && make clean compile_pthread version=$(version)"
	ssh nsc "cd $(DESTINATION) && srun --reservation=fri --cpus-per-task=$(n_threads) --mem-per-cpu=10G ./$(MAIN_PTHREAD)$(version).o $(n_bodies) $(n_steps) $(n_threads)"
	
	@echo "Copying results back to local machine"
	rsync -a --progress nsc:"$(DESTINATION)/output_pthread$(version).txt" ./


run_mpi_remote:
	@test -n "$(n_bodies)" || (echo "n_bodies is not set"; exit 1)
	@test -n "$(n_steps)" || (echo "n_steps is not set"; exit 1)
	@test -n "$(n_tasks)" || (echo "n_tasks is not set"; exit 1)
	@test -n "$(n_nodes)" || (echo "n_nodes is not set"; exit 1)
	@test -n "$(n_cpus)" || (echo "n_cpus is not set"; exit 1)
	@test -n "$(n_threads)" || (echo "n_threads is not set"; exit 1)

	ssh nsc "mkdir -p $(DESTINATION)"

	@echo "Copying files to nsc"
	rsync -a --progress ./Makefile ./n_bodies_mpi.c nsc:"$(DESTINATION)"

	@echo "Compiling and running on nsc"
	ssh nsc 'cd $(DESTINATION); make clean;module load OpenMPI/4.1.0-GCC-10.2.0; make compile_mpi n_threads=$(n_threads);'
	ssh nsc 'module load OpenMPI/4.1.0-GCC-10.2.0; srun --mpi=pmix --reservation=fri --ntasks=$(n_tasks) -N$(n_nodes) --cpus-per-task=$(n_cpus) $(DESTINATION)/$(MAIN_MPI).o $(n_bodies) $(n_steps)'

	@echo "Copying results back to local machine"
	rsync -a --progress nsc:"$(DESTINATION)/output_mpi.txt" ./
