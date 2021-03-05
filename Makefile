
CFLAGS=-O3
CC=g++

all : test_sa

test_sa : needleman_wunsch_sa.h smith_waterman_sa.h utils.h test_sa.cpp
	$(CC)  $(CFLAGS) test_sa.cpp -o $@



clean : 
	rm -f test_sa