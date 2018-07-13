CC=clang++
CFLAGS=-g -Wall -std=c++11 -O0
LIBS=-pthread
PROG=TreeNome
SRC=$(wildcard *.C)
OBJ=$(patsubst %.C, obj/%.o, $(SRC))

all: $(OBJ) $(PROG)
	
$(PROG): $(OBJ)
	$(CC) $(CFLAGS) $(LIBS) -o $@ $?

$(OBJ): obj/%.o: %.C
	$(CC) $(CFLAGS) $(LIBS) -o $@ -c $<

.PHONY: run clean

run:
	./$(PROG) -f ./fastq/frag_2.fastq

clean:
	rm $(OBJ) $(PROG)
