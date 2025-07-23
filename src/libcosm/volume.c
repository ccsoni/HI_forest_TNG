#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cosmology.h"

int main(int argc,char **argv)
{
  struct cosm_model lcdm,scdm;
  double z,D1,D2;

  init_cosm(&lcdm, 0.3, 0.0, 0.7, 0.7, 1.0, 1.0, 1);
  init_cosm(&scdm, 1.0, 0.0, 0.0, 0.7, 1.0, 1.0, 1);

  for(z=0.e0;z<5.0;z+=0.05){
    D1 = dVoverdz(lcdm,z)/pow(3.e3*MPC/lcdm.hubble,3.0);
    D2 = comoving_volume(lcdm,0.00,z)/pow(3.e3*MPC/lcdm.hubble,3.0);
    printf("%f %e %e\n",z,D1,D2);
  }

  for(z=0.e0;z<5.0;z+=0.05){
    D1 = dVoverdz(scdm,z)/pow(3.e3*MPC/scdm.hubble,3.0);
    D2 = comoving_volume(scdm,0.00,z)/pow(3.e3*MPC/scdm.hubble,3.0);
    printf("%f %e %e\n",z,D1,D2);
  }


  exit(0);
}
