CC=g++ #clang++
CFLAGS=-g -Wall -std=c++11 -O0
LIBS=-fopenmp
INC=-I./includes -I./includes/tclap
PROG=TreeNome
SRC=$(wildcard src/*.C)
OBJ=$(patsubst src/%.C, obj/%.o, $(SRC))
#OPTIM?=0
#OPTIMFLAGS=-fsave-optimization-record -foptimization-record-file=./
#ifeq($(OPTIM), 1)
	#$(CC) $(INC) $(CFLAGS) $(LIBS) -o $@ -c $<
#else
	#$(CC) $(INC) $(OPTIMFLAGS) $(CFLAGS) $(LIBS) -o $@ -c $<
#endif

$(PROG): $(OBJ)
	$(CC) $(INC) $(CFLAGS) $(LIBS) -o $@ $^

obj/%.o: src/%.C
	$(CC) $(INC) $(CFLAGS) $(LIBS) -o $@ -c $<

.PHONY: run clean

run:
	./$(PROG) -f ./fastq/frag_2.fastq

clean:
	rm $(OBJ) $(PROG)
