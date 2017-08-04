CC = gcc
cflags = -Wall

main:  graph.o main.c
	$(CC) $@.c $< -o $@ $(cflags)
	
graph.o: graph.c graph.h
	$(CC) -c graph.c -o $@ $(cflags) 
