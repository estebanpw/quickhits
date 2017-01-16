CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall -DVERBOSE
BIN=../bin

all: quickhits

quickhits: quickhits.c
	$(CC) $(CFLAGS) binaryTreeFunctions.c commonFunctions.c -lm quickhits.c -lpthread -lm -o $(BIN)/quickhits

clean:
	rm -rf $(BIN)/quickhits
