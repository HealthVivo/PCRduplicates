
#CC=gcc -Wall
CC=gcc -D_GNU_SOURCE
CFLAGS=-c -Wall
SAMTOOLS=parsebam/samtools-0.1.18/
HAPCUT=parsebam/

all:    
	$(MAKE) -C parsebam/samtools-0.1.18 all
	$(MAKE) -C parsebam hairs
	$(MAKE) PCR

PCR: PCRdups.c 
	$(CC) -I$(SAMTOOLS) -I$(HAPCUT) -g -O2 parsebam/bamread.o parsebam/hapfragments.o parsebam/hashtable.o parsebam/readfasta.o parsebam/readvariant.o -o extract_duplicates PCRdups.c -L$(HAPCUT) -L$(SAMTOOLS) -lbam -lm -lz

clean:
	$(MAKE) -C parsebam/samtools-0.1.18 clean
	$(MAKE) -C parsebam clean
	rm -f extract_duplicates
