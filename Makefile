CC=clang++
LIBS=-g -Wall -std=c++11 -fopenmp=libomp
PROG=TreeNome
SRC=$(wildcard *.C)
OBJ=$(patsubst %.C, obj/%.o, $(SRC))

$(PROG): $(OBJ)
	$(CC) $(LIBS) -o $@ $^

obj/%.o: %.C
	$(CC) $(LIBS) -o $@ -c $<

.PHONY: run clean

run:
	./$(PROG) -f ./fastq/frag_2.fastq

clean:
	rm $(OBJ) $(PROG)
