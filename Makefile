CC=clang++
CFLAGS=-g -Wall -std=c++11 -O0
LIBS=-fopenmp
PROG=TreeNome
SRC=$(wildcard *.C)
OBJ=$(patsubst %.C, obj/%.o, $(SRC))

$(PROG): $(OBJ)
	$(CC) $(CFLAGS) $(LIBS) -o $@ $^

obj/%.o: %.C
	$(CC) $(CFLAGS) $(LIBS) -o $@ -c $<

.PHONY: run clean

run:
	./$(PROG) -f ./fastq/frag_2.fastq

clean:
	rm $(OBJ) $(PROG)
