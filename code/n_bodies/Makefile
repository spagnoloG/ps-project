CC=gcc
CC_MPI=mpicc
CFLAGS=-g -Wall -lm  -ffast-math -O3 -march=native
MPI_FLAGS=-fopenmp
PTHREAD_CFLAGS=-pthread
MAIN=n_bodies
MAIN_PTHREAD=n_bodies_pthread
MAIN_MP=n_bodies_openmp
MAIN_MPI=n_bodies_mpi
DESTINATION=/home/devel/tests

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
	@test -n "$(version)" || (echo "Please specify version"; exit 1)
	@test -n "$(n_threads)" || (echo "Please specify n_threads" && exit 1)

	$(CC_MPI) $(CFLAGS) $(MPI_FLAGS) -DOMP_NUM_THREADS=$(n_threads) $(MAIN_MPI)$(version).c -o $(MAIN_MPI)$(version).o

compile_mp:
	@test -n "$(version)" || (echo "Please specify version"; exit 1)
	@test -n "$(n_threads)" || (echo "Please specify n_threads" && exit 1)

	$(CC) $(CFLAGS) $(CFLAGS) $(MPI_FLAGS) -DOMP_NUM_THREADS=$(n_threads) $(MAIN_MP)$(version).c -o $(MAIN_MP)$(version).o

clean: 
	rm *.o || true
	#rm *.txt || true

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

	ssh hsrv_devel "mkdir -p $(DESTINATION)"

	@echo "Copying files to hsrv_devel"
	rsync -a --progress ./Makefile ./n_bodies$(version).c hsrv_devel:$(DESTINATION)

	@echo "Compiling on hsrv_devel"
	ssh hsrv_devel "cd $(DESTINATION) && make clean compile_default version=$(version)"

	@echo "Running on hsrv_devel"

	ssh hsrv_devel "cd $(DESTINATION) &&  ./$(MAIN)$(version).o $(n_bodies) $(n_steps)"
	
	#@echo "Copying files from hsrv_devel"
	#rsync -a --progress hsrv_devel:$(DESTINATION)/n_bodies$(version).txt ./

run_mp_remote:
	@test -n "$(n_bodies)" || (echo "n_bodies is not set"; exit 1)
	@test -n "$(n_steps)" || (echo "n_steps is not set"; exit 1)
	@test -n "$(n_threads)" || (echo "n_threads is not set"; exit 1)
	@test -n "$(version)" || (echo "version is not set"; exit 1)

	ssh hsrv_devel "mkdir -p $(DESTINATION)"
	
	@echo "Copying files to hsrv_devel"
	rsync -a --progress ./Makefile ./n_bodies_openmp$(version).c hsrv_devel:$(DESTINATION)

	@echo "Compiling on hsrv_devel"
	ssh hsrv_devel "cd $(DESTINATION) && make clean compile_mp n_threads=$(n_threads) version=$(version)"

	@echo "Running on hsrv_devel"
	ssh hsrv_devel "cd $(DESTINATION) && OMP_NUM_THREADS=$(n_threads) ./$(MAIN_MP)$(version).o $(n_bodies) $(n_steps)"
	
	@echo "Copying the results back"
	#rsync -a --progress hsrv_devel:"$(DESTINATION)/output_mp.txt" ./

run_pthread_remote:
	@test -n "$(n_bodies)" || (echo "n_bodies is not set"; exit 1)
	@test -n "$(n_threads)" || (echo "n_threads is not set"; exit 1)
	@test -n "$(n_steps)" || (echo "n_steps is not set"; exit 1)
	@test -n "$(version)" || (echo "version is not set"; exit 1)
	
	ssh hsrv_devel "mkdir -p $(DESTINATION)"
	
	@echo "Copying files to hsrv_devel"
	rsync -a --progress ./Makefile ./n_bodies_pthread$(version).c  hsrv_devel:"$(DESTINATION)"

	@echo "Compiling and running on hsrv_devel"
	ssh hsrv_devel "cd $(DESTINATION) && make clean compile_pthread version=$(version)"
	ssh hsrv_devel "cd $(DESTINATION) &&  ./$(MAIN_PTHREAD)$(version).o $(n_bodies) $(n_steps) $(n_threads)"
	
	#@echo "Copying results back to local machine"
	#rsync -a --progress hsrv_devel:"$(DESTINATION)/output_pthread$(version).txt" ./


run_mpi_remote:
	@test -n "$(n_bodies)" || (echo "n_bodies is not set"; exit 1)
	@test -n "$(n_steps)" || (echo "n_steps is not set"; exit 1)
	@test -n "$(n_proc)" || (echo "n_proc is not set"; exit 1)
	#@test -n "$(n_nodes)" || (echo "n_nodes is not set"; exit 1)
	@test -n "$(n_threads)" || (echo "n_threads is not set"; exit 1)
	@test -n "$(version)" || (echo "version is not set"; exit 1)

	ssh hsrv_devel "mkdir -p $(DESTINATION)"

	@echo "Copying files to hsrv_devel"
	rsync -a --progress ./Makefile ./n_bodies_mpi$(version).c hsrv_devel:"$(DESTINATION)"

	@echo "Compiling and running on hsrv_devel"
	ssh hsrv_devel 'cd $(DESTINATION); make clean; make compile_mpi n_threads=$(n_threads) version=$(version);'
	ssh hsrv_devel 'cd $(DESTINATION); OMP_NUM_THREADS=$(n_threads) mpirun -np $(n_proc) $(DESTINATION)/$(MAIN_MPI)$(version).o $(n_bodies) $(n_steps)'

	#@echo "Copying results back to local machine"
	#rsync -a --progress hsrv_devel:"$(DESTINATION)/output_mpi.txt" ./
