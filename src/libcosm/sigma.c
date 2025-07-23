#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cosmology.h"

int main(int argc,char **argv)
{
  struct cosm_model lcdm,scdm;
  double mass,D1;

  init_cosm(&lcdm, 0.2, 0.0, 0.8, 0.7, 1.0, 1.0, 1);
  init_cosm(&scdm, 1.0, 0.0, 0.0, 0.7, 1.0, 1.0, 1);

  for(mass=1.e9;mass<1.e20;mass*=1.5){
    D1 = sigmaM(lcdm,mass);
    printf("%e %f\n",mass/lcdm.omega0*lcdm.hubble,D1);
  }

  /*   for(mass=1.e9;mass<1.e20;mass*=1.5){
       D1 = sigmaM(scdm,mass);
       printf("%e %f\n",mass/scdm.omega0*scdm.hubble,D1);
       }*/

  exit(0);
}

