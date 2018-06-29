CC=clang++-3.8
LIBS=-g -Wall -std=c++11
PROG=TreeNome
SRC=$(wildcard *.C)
OBJ=$(patsubst %.C, obj/%.o, $(SRC))

all: $(OBJ) $(PROG)
	
$(PROG): $(OBJ)
	$(CC) $(LIBS) -o $@ $?

$(OBJ): obj/%.o: %.C
	$(CC) $(LIBS) -o $@ -c $<

.PHONY: clean

clean:
	rm $(OBJ) $(PROG)
