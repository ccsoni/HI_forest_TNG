#include <stdio.h>
#include <stdlib.h>
#include "hdf5.h"
#include "tng_data.h"
#include "los_1D.h"
#include "libcosm/cosmology.h"

int main(int argc, char **argv)
{

  hid_t file;
  struct tng_header head;
  struct tng_gas *cell;
  struct tng_cdm *cdm;

  struct cosm_model lcdm;

  struct los_1D los;
  struct line_segment *seg;
  struct vector lseg;
  struct point start_pnt;
  int nseg;

  if(argc != 8) {
    fprintf(stderr, "Usage : %s prefix x_0 y_0 z_0 dx dy dz,\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  static char filename[128];
  sprintf(filename, "%s.0.hdf5", argv[1]);

  file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  read_header(&head, file);
  show_header(&head);

  init_cosm(&lcdm, head.omegam, head.omegab, head.omegav, head.hubble,
	    0.9667, 0.8159, 1);

  // number of segments 
  nseg = 8;
  seg = (struct line_segment *)malloc(sizeof(struct line_segment)*nseg);

  // starting point of the line segment (here, the origin of the simulatin box)
  start_pnt.xpos = atof(argv[2]);
  start_pnt.ypos = atof(argv[3]);
  start_pnt.zpos = atof(argv[4]);
  
  // unit vector along the line segment 
  lseg.x = atof(argv[5]);
  lseg.y = atof(argv[6]);
  lseg.z = atof(argv[7]);
  normalize(&lseg);

  int ok = check_line_segment_setting(start_pnt, lseg, nseg, &head);
  if(ok != 0) fprintf(stderr, "# Invalid setting for the line segment setting.\n");

  setup_line_segment(seg, nseg, start_pnt, lseg, &head, lcdm);

  H5Fclose(file);

  zero_out_los(&los);
  setup_los(&los, seg, nseg, &head, lcdm);
  
  int dm_io_flag, gas_io_flag;

  dm_io_flag = gas_io_flag = 0;

  //  dm_io_flag |= _DM_POS;
  gas_io_flag |= (_GAS_POS| _GAS_MASS | _GAS_DENS | _GAS_VEL | _GAS_FHI | _GAS_UENE | _GAS_FELE );

  for(int ichunk=0;ichunk<head.nchunk;ichunk++){
    sprintf(filename, "%s.%d.hdf5", argv[1], ichunk);
    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);    
    read_header(&head, file);

    printf("# ichunk = %d / %d \n", ichunk, head.nchunk);fflush(stdout);
    cell = (struct tng_gas *) malloc(sizeof(struct tng_gas)*head.np_local[0]);

    read_gas_cell(cell, file, &head, gas_io_flag);
    printf("# data read\n");fflush(stdout);
    calc_los_tau(cell, &los, seg, nseg, lcdm, &head);
    printf("# optical depth computed\n");fflush(stdout);

    free(cell);
  }

  for(int iz=0;iz<NMESH_LOS;iz++) {
    los.tau[iz]  *= HI21_K;
    los.vel[iz]  /= (los.column_mass[iz]+1.0e-33);
    los.tmpr[iz] /= (los.column_mass[iz]+1.0e-33);
    los.fHI[iz]  /= (los.column_mass[iz]+1.0e-33);
    //los.column_mass[iz] /= (los.path_length[iz]*MPROTON+1.0e-33);
    los.NHI[iz] *= XHYDROG/MPROTON;
    los.dens[iz] /= (los.count[iz]+1.0e-33);
  }

  char output_los_filename[256];

  sprintf(output_los_filename, "%s_los.dat", argv[1]);
  output_los_data(&los, output_los_filename, 1);  
}
