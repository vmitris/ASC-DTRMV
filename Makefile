CC = gcc

CFLAGS = -Wall -g -O5
CFLAGS1 = /opt/tools/libraries/atlas/3.10.1-opteron-gcc-4.4.6/lib
CFLAGS2 = /opt/tools/libraries/atlas/3.10.1-opteron-gcc-4.4.6/include

CFLAGS3 = /opt/tools/libraries/atlas/3.10.1-nehalem-gcc-4.4.6/lib
CFLAGS4 = /opt/tools/libraries/atlas/3.10.1-nehalem-gcc-4.4.6/include

CFLAGS5 = /opt/tools/libraries/atlas/3.10.1-quad-gcc-4.4.6/lib
CFLAGS6 = /opt/tools/libraries/atlas/3.10.1-quad-gcc-4.4.6/include

build:
	./run.sh

build1:

	$(CC) $(CFLAGS) -o dtrmv -L $(CFLAGS1) dtrmv.c -I $(CFLAGS2) -l cblas -l atlas -lm

build2:

	$(CC) $(CFLAGS) -o dtrmv -L $(CFLAGS3) dtrmv.c -I $(CFLAGS4) -l cblas -l atlas -lm

build3:

	$(CC) $(CFLAGS) -o dtrmv -L $(CFLAGS5) dtrmv.c -I $(CFLAGS6) -l cblas -l atlas -lm

clean:
	rm -fr dtrmv