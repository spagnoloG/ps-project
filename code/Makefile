CC=gcc
CFLAGS=-g -Wall -lm

all: main 

server: main.c
	$(CC) $(CFLAGS) main.c -o main 

clean: 
	rm main
run:
	./main
