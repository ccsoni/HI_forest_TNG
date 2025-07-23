#CC=gcc
#CFLAGS=-O3 -mavx2 -ffast-math -funroll-loops -std=c99 -fopenmp
CC=icx
CFLAGS=-O3 -qopenmp -funroll-loops -ffast-math

HDF5_ROOT=/home/kohji/tmp/hdf5-1.8.22/hdf5
CFLAGS+=-I$(HDF5_ROOT)/include
HDF5LIB=-L$(HDF5_ROOT)/lib -lhdf5_hl -lhdf5 -lz

MAP_OBJ = map_main.o map_2D.o tng_input_data.o
MAP_DEP := $(MAP_OBJ)

HI_FOREST_OBJ = HI_forest.o los_1D.o tng_input_data.o graph.o
HI_FOREST_DEP := $(HI_FOREST_OBJ)

READ_SPECTRUM_OBJ = read_spectrum.o los_1D.o graph.o
READ_SPECTRUM_DEP := $(READ_SPECTRUM_OBJ)

#CFLAGS+=-D__MASS_WEIGHTED_MAP__

ALL = map_2D HI_forest read_spectrum

all: $(ALL)

map_2D: $(MAP_DEP)
	$(CC) $(CFLAGS) -o $@ $(MAP_OBJ) $(HDF5LIB) -lm

HI_forest: $(HI_FOREST_DEP)
	$(CC) $(CFLAGS) -o $@ $(HI_FOREST_OBJ) $(HDF5LIB) -Llibcosm -lcosm -lm

read_spectrum: $(READ_SPECTRUM_DEP)
	$(CC) $(CFLAGS) -o $@ $(READ_SPECTRUM_OBJ) -Llibcosm -lcosm -lm

clean:
	\rm -rf *.o

distclean: clean
	\rm -rf *~
	\rm -rf $(ALL)
