CC=g++

all: obj exec

obj: 
	$(CC) -std=c++11 -c src/DNATranslator.cpp src/FASTAParse.cpp src/SequenceAlignment.cpp src/main.cpp 

exec: 
	$(CC) -lpthread -o test main.o DNATranslator.o FASTAParse.o SequenceAlignment.o 

.PHONY: clean

clean: