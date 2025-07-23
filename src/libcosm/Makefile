CC=cc
CFLAGS=-std=c99 -O3
CLIBS=-L. -lcosm -lm
AR=ar

#ALL=libcosm.a distance formation_rate growth_rate \
#	mass_func power_spec sigma volume

ALL=libcosm.a distance

libcosm.a: cosmology.o
	$(AR) rv libcosm.a cosmology.o
	ranlib libcosm.a

distance: distance.o libcosm.a
	$(CC) $(CFLAGS) $^ -o $@ $(CLIBS)

.c.o: 
	$(CC) $(CFLAGS) -c $*.c -o $*.o

.PHONY: clean distclean

clean:
	\rm -rf *.o

distclean: clean
	\rm -rf $(ALL)
	\rm -rf *.o *~
