CC        = cc
CFLAGS    = -O3 -fomit-frame-pointer
LIBS      = -lm

all:            fitch_kmerX consense treedist fitch_kmerX_long
phylip.o:       phylip.h
dist.o:         phylip.h dist.h
dist_long.o:         phylip.h dist.h
fitch_kmerX.o:  fitch_kmerX.c dist.h phylip.h
fitch_kmerX_long.o:  fitch_kmerX_long.c dist.h phylip.h

fitch_kmerX:    fitch_kmerX.o dist.o phylip.o
	$(CC) $(CFLAGS) fitch_kmerX.o dist.o phylip.o $(LIBS) -o fitch_kmerX
fitch_kmerX_long:    fitch_kmerX_long.o dist_long.o phylip.o
	$(CC) $(CFLAGS) fitch_kmerX_long.o dist_long.o phylip.o $(LIBS) -o fitch_kmerX_long
cons.o:         cons.c cons.h phylip.h 
consense.o:     consense.c cons.h phylip.h

consense:       consense.o phylip.o cons.o
	$(CC) $(CFLAGS) consense.o phylip.o cons.o $(LIBS) -o consense
treedist.o:     treedist.c cons.h phylip.h

treedist:     treedist.o phylip.o cons.o
	$(CC) $(CFLAGS) treedist.o cons.o phylip.o $(LIBS) -o treedist

clean:
	rm -rf *o fitch_kmerX fitch_kmerX_long consense treedist

