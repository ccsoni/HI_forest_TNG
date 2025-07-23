#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cosmology.h"

int main(int argc,char **argv)
{
  struct cosm_model lcdm,scdm;
  double k,D1;

  init_cosm(&lcdm, 0.357, 0.0, 0.0, 0.7, 1.0, 1.0, 1);
  init_cosm(&scdm, 1.0, 0.0, 0.0, 0.7, 1.0, 1.0, 1);

  /*   for(k=1.e-3;k<1.e1;k*=1.5){
       D1 = power_spec(lcdm,k);
       printf("%f %f\n",k,D1);
       }*/

  for(k=3.e-4;k<1.e1;k*=1.2){
    D1 = power_spec(scdm,k);
    printf("%f %f\n",k,D1);
  }


  exit(0);
}
