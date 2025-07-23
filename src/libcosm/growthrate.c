#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cosmology.h"

int main(int argc,char **argv)
{
  struct cosm_model lcdm,scdm;
  double z,D1,delta;

  init_cosm(&lcdm, 0.2, 0.0, 0.8, 0.7, 1.0, 1.0, 1);
  init_cosm(&scdm, 1.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1);

  for(z=10.0;z>0.0;z-=0.05){
    D1 = d_growthrate_dt(lcdm,z);
    delta = Delta_c(lcdm,z);
    /*      printf("%f %e %f\n",z,D1,delta);*/
  }

  for(z=10.0;z>0.0;z-=0.05){
    D1 = d_growthrate_dt(scdm,z);
    delta = Delta_c(scdm,z);
    printf("%f %e %f\n",z,D1,delta);
  }


  exit(0);
}
