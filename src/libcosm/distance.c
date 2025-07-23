#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cosmology.h"

int main(int argc,char **argv)
{
  struct cosm_model lcdm,scdm;


  init_cosm(&lcdm, 0.3, 0.0, 0.7, 0.7, 1.0, 1.0, 1);
  init_cosm(&scdm, 0.3, 0.0, 0.0, 0.7, 1.0, 1.0, 1);

  for(double z=0.0;z<10.0;z+=0.05){
    double D1 = angular_distance(lcdm, z)/MPC;
    double D2 = luminosity_distance(lcdm, z)/MPC;
    double D3 = comoving_distance(lcdm, z)/MPC;
    printf("%f %e %e %e\n",z,D1,D2,D3);
  }

  for(double z=0.0;z<10.0;z+=0.05){
    double D1 = angular_distance(scdm, z)/MPC;
    double D2 = luminosity_distance(scdm, z)/MPC;
    /*      printf("%f %e \n",z,D1);*/
  }


  exit(0);
}
