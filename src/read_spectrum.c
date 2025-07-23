#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include "hdf5.h"
#include "tng_data.h"
#include "los_1D.h"
#include "libcosm/cosmology.h"

int main(int argc, char **argv)
{
  struct los_1D *los;
  int nlos;

  if(argc != 3){
    fprintf(stderr, "Usage: %s <data file> <output_prefix>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  input_los_data_header(argv[1], &nlos);

  printf("# of line-of-sights: %d\n", nlos);

  los = (struct los_1D *) malloc(nlos*sizeof(struct los_1D));

  input_los_data(los, argv[1], nlos);

  FILE *fp_spec, *fp_los;
  static char spec_name[256], los_name[256];

  sprintf(spec_name, "%s.spec", argv[2]);
  sprintf(los_name, "%s.los", argv[2]);
  
  fp_spec = fopen(spec_name, "w");
  fp_los  = fopen(los_name, "w");

  for(int inu=0;inu<NMESH_LOS;inu++) {
    double nu = los[0].nu_min + (double)(inu+0.5)*los[0].delta_nu;
    double zred = los[0].zred_min + (double)(inu+0.5)*los[0].dzred;

    fprintf(fp_spec,"%14.6e %14.6e\n", nu, los[0].tau[inu]);

    fprintf(fp_los,"%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %d\n",
	   zred, los[0].tmpr[inu], los[0].fHI[inu], los[0].vel[inu], los[0].column_mass[inu], los[0].NHI[inu], los[0].count[inu]);

  }

  fflush(fp_spec);
  fflush(fp_los);
  fclose(fp_spec);
  fclose(fp_los);
}
