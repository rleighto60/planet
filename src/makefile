#
CC = gcc -O3 -Wall -fomit-frame-pointer
LIBS = -lregina -lm -ldl

OFILES = planet.o parse.o

../bin/planet:	$(OFILES)
	$(CC) -o ../bin/planet $(OFILES) $(LIBS)

planet.o:  planet.c planet.h
	$(CC) -c -o planet.o planet.c

parse.o:    parse.c
	$(CC) -c -o parse.o parse.c

clean:
	rm -f *.o
