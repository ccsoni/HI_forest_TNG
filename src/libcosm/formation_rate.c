#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cosmology.h"

double ps_total_rate(struct cosm_model,double,double);
double ps_formation_rate(struct cosm_model,double,double);
double ps_destruction_rate(struct cosm_model,double,double);

int main(int argc,char **argv)
{
  struct cosm_model lcdm,scdm;
  double mass,total_rate,form_rate,psmass;
  double redshift;

  /* init_cosm(struct cosm_model*,double omega0,double omegab,double lambda0,
     double hubble,double pkindex,double sigma8,int DM);*/
  init_cosm(&scdm, 1.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1);
  redshift = 0.0;

  for(mass=1.e8;mass<1.e16;mass*=1.1){
    psmass = psmassfunc(scdm,mass,redshift);
    total_rate = ps_total_rate(scdm,mass,redshift);
    form_rate = ps_formation_rate(scdm,mass,redshift);
    printf("%e\t%e\t%e\t%e\n",log10(mass),log10(fabs(total_rate))
	   ,log10(form_rate),log10(psmass));
  }

  exit(0);
}

double ps_total_rate(struct cosm_model cosm,double Mass,double z)
{
  double total_rate;
   
  double sigma,delta,delta_0;
  double dndM,D1,dD1_dz,scale,k,Hubble,dtdz;

  scale = 1.e0/(1.e0+z);
  k = cosm.omega0+cosm.lambda0-1.e0;
  Hubble = H0*cosm.hubble*
    sqrt(cosm.omega0/pow(scale,3.0)-k/pow(scale,2.0)+cosm.lambda0);
  /*    dtdz = (1.e0+scale)*(1.e0+scale)/(scale*Hubble); */
  dtdz = scale/Hubble;

  sigma = sigmaM(cosm,Mass);
  delta = delta_c(cosm,z);
  delta_0 = delta_c(cosm,0.e0);
  D1 = growthrate(cosm,z);
  dD1_dz = d_growthrate_dt(cosm,z)*dtdz;

  dndM = psmassfunc(cosm,Mass,z);
  total_rate = dndM*(delta_0/D1/(sigma*sigma)-D1/delta_0)*
    (delta_0*dD1_dz/(D1*D1));

  return total_rate;
}

double ps_formation_rate(struct cosm_model cosm,double Mass,double z)
{
  double M_form;
  double sigma_form,sigma;
  double dndM,D1,dD1_dz,scale,k,Hubble,dtdz;
  double delta,delta_0;
  double formation_rate;

  M_form = Mass/2.e0;

  scale = 1.e0/(1.e0+z);
  k = cosm.omega0+cosm.lambda0-1.e0;
  Hubble = H0*cosm.hubble*
    sqrt(cosm.omega0/pow(scale,3.0)-k/pow(scale,2.0)+cosm.lambda0);
  dtdz = scale/Hubble;
   
  delta = delta_c(cosm,z);
  delta_0 = delta_c(cosm,0.e0);
  D1 = growthrate(cosm,z);
  dD1_dz = d_growthrate_dt(cosm,z)*dtdz;
  sigma_form = sigmaM(cosm,M_form);
  sigma      = sigmaM(cosm,Mass);
  dndM = psmassfunc(cosm,Mass,z);

  formation_rate = sqrt(2.e0/PI)*dndM*(delta_0*dD1_dz/(D1*D1))
    /sqrt(sigma_form*sigma_form-sigma*sigma);

  return formation_rate;
}
