CC=g++ #clang++
CFLAGS=-g -Wall -std=c++11 -O0
LIBS=-fopenmp
INC=-I./includes -I./includes/tclap
PROG=TreeNome
SRC=$(wildcard src/*.C)
OBJ=$(patsubst src/%.C, obj/%.o, $(SRC))
OPTIM=-O3

$(PROG): $(OBJ)
	$(CC) $(INC) $(CFLAGS) $(LIBS) $(OPTIM) -o $@ $^

obj/%.o: src/%.C
	$(CC) $(INC) $(CFLAGS) $(LIBS) $(OPTIM) -o $@ -c $<

.PHONY: run clean

run:
	./$(PROG) -f ./fastq/example.fastq -s ./data_out/test.gno

clean:
	rm $(OBJ) $(PROG)
