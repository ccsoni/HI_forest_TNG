#include <stdio.h>
#include <stdlib.h>
#include "hdf5.h"
#include "tng_data.h"
#include "map_2D.h"

int main(int argc, char **argv)
{

  hid_t file;
  struct tng_header head;
  struct tng_gas *cell;
  struct tng_cdm *cdm;

  struct map_2D _map;

  if(argc != 2) {
    fprintf(stderr, "Usage : %s prefix \n", argv[0]);
    exit(EXIT_FAILURE);
  }

  static char filename[128];
  sprintf(filename, "%s.0.hdf5", argv[1]);

  file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  read_header(&head, file);
  show_header(&head);

  H5Fclose(file);

  int dm_io_flag, gas_io_flag;

  dm_io_flag = gas_io_flag = 0;

  //  dm_io_flag |= _DM_POS;
  gas_io_flag |= (_GAS_POS | _GAS_MASS | _GAS_FHI | _GAS_UENE | _GAS_FELE | _GAS_DENS);

  _map.slice_pos = head.boxsize/2.0;
  zero_out_map(&_map);

  for(int ichunk=0;ichunk < head.nchunk;ichunk++) {
    sprintf(filename, "%s.%d.hdf5", argv[1], ichunk);
    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);    
    read_header(&head, file);

    printf("ichunk = %d\n", ichunk);
  
    cell = (struct tng_gas *) malloc(sizeof(struct tng_gas)*head.np_local[0]);
    cdm  = (struct tng_cdm *) malloc(sizeof(struct tng_cdm)*head.np_local[1]);

    read_gas_cell(cell, file, &head, gas_io_flag);
    read_cdm(cdm, file, &head, dm_io_flag);    
    
    accumulate_gas_to_map(&_map, cell, &head);

    free(cdm);
    free(cell);

    H5Fclose(file);
  }

  normalize_map(&_map, &head);

  char output_filename[256];
  sprintf(output_filename, "%s_map.dat", argv[1]);
  output_map(output_filename, &_map, &head);
  
  exit(0);  

}
