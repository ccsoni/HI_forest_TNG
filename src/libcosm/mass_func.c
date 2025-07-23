#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cosmology.h"

int main(int argc,char **argv)
{
  struct cosm_model lcdm,scdm;
  double mass,D1;

  init_cosm(&lcdm, 0.3, 0.0, 0.7, 0.7, 1.0, 1.0, 1);
  init_cosm(&scdm, 1.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1);

  for(mass=1.e8;mass<1.e16;mass*=1.1){
    D1 = psmassfunc(lcdm,mass,0.0);
    printf("%e %e\n",mass,D1);
  }

  /*   for(mass=3.e8;mass<3.e16;mass*=1.1){
       D1 = psmassfunc(scdm,mass,10.0);
       printf("%e %e\n",mass,D1);
       }*/


  exit(0);
}
