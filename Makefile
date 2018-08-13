CC=g++
CFLAGS=-g -Wall -std=c++11 -O0
LIBS=-fopenmp
INC=-I./includes -I./includes/tclap
PROG=TreeNome
SRC=$(wildcard src/*.C)
OBJ=$(patsubst src/%.C, obj/%.o, $(SRC))

$(PROG): $(OBJ)
	$(CC) $(INC) $(CFLAGS) $(LIBS) -o $@ $^

obj/%.o: src/%.C
	$(CC) $(INC) $(CFLAGS) $(LIBS) -o $@ -c $<

.PHONY: run clean

run:
	./$(PROG) -f ./fastq/frag_2.fastq

clean:
	rm $(OBJ) $(PROG)
