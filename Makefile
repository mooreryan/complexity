CFLAGS=-Wall -g -std=c11 -lmpfr -lgmp

so:
	gcc -Wall -c complexity.c -o complexity.o
	gcc -lmpfr -lgmp -shared -o complexity.so complexity.o
