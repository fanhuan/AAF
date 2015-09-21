CC=gcc
CFLAGS=-g -W -Wall -O2 -D_FILE_OFFSET_BITS=64
GLIBS=-lm
GENERIC_SRC=string.h vector.h file_reader.h hashset.h

all: kmer_count kmer_countx kmer_merge

kmer_count: $(GENERIC_SRC) file_reader.c kmer_count.c
	$(CC) $(CFLAGS) -o kmer_count file_reader.c kmer_count.c $(GLIBS)

kmer_countx: $(GENERIC_SRC) file_reader.c kmer_countx.c
	$(CC) $(CFLAGS) -o kmer_countx file_reader.c kmer_countx.c $(GLIBS)

kmer_merge: $(GENERIC_SRC) file_reader.c filter.c
	$(CC) $(CFLAGS) -o kmer_merge file_reader.c filter.c $(GLIBS)

clean:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out kmer_count kmer_countx kmer_merge

clear:
	rm -f *.o *.gcda *.gcno *.gcov gmon.out
