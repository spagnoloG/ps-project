CC=gcc
CFLAGS=-g -Wall -lm
PTHREAD_CFLAGS=-pthread
MAIN=n_bodies
MAIN_PTHREAD=n_bodies_pthread

all: $(MAIN) $(MAIN_PTHREAD)

$(MAIN): $(MAIN).c
	$(CC) $(CFLAGS) $(MAIN).c -o $(MAIN).o

$(MAIN_PTHREAD): $(MAIN_PTHREAD).c
	$(CC) $(CFLAGS) -pthread $(MAIN_PTHREAD).c -o $(MAIN_PTHREAD).o

clean: 
	rm $(MAIN).o $(MAIN_PTHREAD).o

run_local_main:
	./$(MAIN).o

run_local_main_pthread:
	./$(MAIN_PTHREAD).o



