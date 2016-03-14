CC=gcc -std=c99  #for snow leopard
#C89=cc
#C99=c99
#CC = $(C99)
CFLAGS_MAIN=-O3
CFLAGS_IIO =-O3 -D_GNU_SOURCE -DNDEBUG

OMPFLAGS=-fopenmp

all: iio.o xmalloc.o TV_ROF_Box.o  utils.o tvl1Occflow

xmalloc.o: xmalloc.h xmalloc.c
	$(CC) $(CFLAGS_IIO) $(OMPFLAGS) -c xmalloc.c

iio.o: iio.c
	$(CC) $(CFLAGS_IIO) -c iio.c

TV_ROF_Box.o: constants.h TV_ROF_Box.h TV_ROF_Box.c operators.h operators.c
	$(CC) $(CFLAGS_IIO) $(OMPFLAGS) -c operators.c TV_ROF_Box.c

utils.o: constants.h utils.h utils.c pyramid.h pyramid.c warp.h warp.c
	$(CC) $(CFLAGS_IIO) $(OMPFLAGS) -c utils.c pyramid.c warp.c

tvl1Occflow: constants.h main.c tvl1OccFlow_lib.c bicubic_interpolation.c solvers.c operators.c  iio.o TV_ROF_Box.o
	$(CC) $(CFLAGS_MAIN) $(OMPFLAGS) -o tvl1Occflow  bicubic_interpolation.c tvl1OccFlow_lib.c  main.c solvers.c xmalloc.o iio.o utils.o  operators.o warp.o pyramid.o TV_ROF_Box.o -lm -lpng -ljpeg -ltiff


clean:
	rm -f *.o tvl1Occflow
