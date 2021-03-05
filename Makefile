
CFLAGS=-O3
CC=g++

all : test_sw_sa

test_sw_sa : smith_waterman_sa.h test_sw_sa.cpp
	$(CC)  $(CFLAGS) test_sw_sa.cpp -o $@

clean : 
	rm -f test_sw_sa