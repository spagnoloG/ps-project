CC=gcc
CFLAGS=-g -Wall -lm
PTHREAD_CFLAGS=-pthread -ffast-math -O3 -march=native
MAIN=n_bodies
MAIN_PTHREAD=n_bodies_pthread
USERNAME=$(shell ssh nsc "whoami")
DESTINATION=/ceph/grid/home/$(USERNAME)/n_bodies

all: $(MAIN) $(MAIN_PTHREAD)

$(MAIN): $(MAIN).c
	$(CC) $(CFLAGS) $(MAIN).c -o $(MAIN).o

$(MAIN_PTHREAD): $(MAIN_PTHREAD).c
	$(CC) $(CFLAGS) $(PTHREAD_CFLAGS) $(MAIN_PTHREAD).c -o $(MAIN_PTHREAD).o

clean: 
	@test -f $(MAIN).o && rm $(MAIN).o || true
	@test -f $(MAIN_PTHREAD).o && rm $(MAIN_PTHREAD).o || true

run_local_main:
	./$(MAIN).o

run_local_main_pthread:
	./$(MAIN_PTHREAD).o

run_pthread_remote:
	@test -n "$(n_bodies)" || (echo "n_bodies is not set"; exit 1)
	@test -n "$(n_threads)" || (echo "n_threads is not set"; exit 1)
	@test -n "$(n_steps)" || (echo "n_steps is not set"; exit 1)
	
	ssh nsc "mkdir -p $(DESTINATION)"
	
	@echo "Copying files to nsc"
	rsync -a --progress ./Makefile nsc:"$(DESTINATION)"
	rsync -a --progress ./n_bodies.c nsc:"$(DESTINATION)"
	rsync -a --progress ./n_bodies_pthread.c nsc:"$(DESTINATION)"

	@echo "Compiling and running on nsc"
	ssh nsc "cd $(DESTINATION) && make clean"
	ssh nsc "cd $(DESTINATION) && make"
	ssh nsc "cd $(DESTINATION) && srun --reservation=fri --cpus-per-task=$(n_threads) --mem-per-cpu=10G ./$(MAIN_PTHREAD).o $(n_bodies) $(n_steps) $(n_threads)"
	
	@echo "Copying results back to local machine"
	rsync -a --progress nsc:"$(DESTINATION)/output_pthread.txt" ./output_pthread.txt
