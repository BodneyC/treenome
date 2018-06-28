CC=g++
LIBS=-g -Wall
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
