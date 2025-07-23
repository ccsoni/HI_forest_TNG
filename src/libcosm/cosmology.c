#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cosmology.h"

/* prototypes */
void init_cosm(struct cosm_model*, double, double, double, double, double, double, int);
double growthrate(struct cosm_model,double);
double d_growthrate_dt(struct cosm_model,double);
double power_spec(struct cosm_model,double);
double delta_c(struct cosm_model,double);
double Delta_c(struct cosm_model,double);
double Transfer(struct cosm_model,double);
double sigmaR(struct cosm_model,double);
double sigmaM(struct cosm_model,double);
double psmassfunc(struct cosm_model,double,double);
double cumulative_psmassfunc(struct cosm_model,double,double);
double d_kai(struct cosm_model,double);
double kai_distance(struct cosm_model,double);
double comoving_distance(struct cosm_model,double);
double angular_distance(struct cosm_model,double);
double luminosity_distance(struct cosm_model,double);
double dVoverdz(struct cosm_model,double);
double comoving_volume(struct cosm_model,double,double);

double cosm_trapzd(double (*func)(struct cosm_model,double), struct cosm_model,double,double,int);
void   cosm_polint(double*,double*,int,double,double*,double*);
/* -------- */

void init_cosm(struct cosm_model *cosm, double om, double omb, double lmd, double h
	       ,double pkindx, double sigm8, int dm)
{
  double sigma2;
  double kmin,kmax;
  double delta,kr,k;
  int nbin,p;

  cosm->omega_m  = om;
  cosm->omega_b  = omb;
  cosm->omega_v = lmd;
  cosm->hubble  = h;
  cosm->pkindex = pkindx;
  cosm->sigma8  = sigm8;
  cosm->DM      = dm;

  cosm->shapeparam = (cosm->omega_m)*(cosm->hubble)
    *exp(-(cosm->omega_b)*(1.e0+sqrt(2.e0*(cosm->hubble))/(cosm->omega_m)));

  kmin = 1.e-5;kmax = 1.e3;
  nbin = 1000;

  delta  = (log(kmax)-log(kmin))/(double)nbin;
  sigma2 = 0.e0;

  for(p=1;p<=nbin;p++){
    k       = kmin*exp(((float)p-0.5)*delta);
    kr      = k*8.0/(cosm->hubble);
    sigma2 += Transfer(*cosm,k)*pow(k,(cosm->pkindex))*pow(k,3.e0)
      *9.e0*pow((sin(kr)-kr*cos(kr))/pow(kr,3.e0),2.e0);
  }   

  sigma2 *= delta;

  cosm->pk_amplitude = 2.e0*PI*PI/sigma2*(cosm->sigma8)*(cosm->sigma8);
  cosm->cobe_amplitude = 1.67e-5/2.7*PI*sqrt(1.2)*4.e0*pow(2997.9/(cosm->hubble),2.e0)
    *pow((cosm->omega_m),-0.77);
  cosm->cobe_amplitude = (cosm->cobe_amplitude)*(cosm->cobe_amplitude);
}

double growthrate(struct cosm_model cosm,double z)
{
  double D1;
  double x,x0,D10;

  if(cosm.omega_m>0.99 && cosm.omega_v<0.01) D1 = 1.e0/(1.e0+z);
  if(cosm.omega_m<1.e0 && cosm.omega_v<0.01) {
    x  = (1.e0/cosm.omega_m-1.e0)/(1.e0+z);
    x0 = (1.e0/cosm.omega_m-1.e0);
    D10 = 1.0 + 3.0/x0 
      + 3.0*sqrt((1.0+x0)/pow(x0,3.e0))*log(sqrt(1.0+x0)-sqrt(x0));
    D1 = 1.0 + 3.0/x 
      + 3.0*sqrt((1.0+x)/pow(x,3.e0))*log(sqrt(1.0+x)-sqrt(x));
    D1 = D1/D10;
  }
  if(cosm.omega_m<1.e0 && cosm.omega_v>0.01){
    x  = pow(2.0*(1.0/cosm.omega_m-1.0),1.0/3.0)/(1.0+z);
    x0 = pow(2.0*(1.0/cosm.omega_m-1.0),1.0/3.0);
    D10 = (0.2*x0+0.044*pow(x0,1.37))/(1.0+0.328*pow(x0,1.37));
    D1 = (0.2*x+0.044*pow(x,1.37))/(1.0+0.328*pow(x,1.37));
    D1 = D1/D10;
  }
   
  return D1;
}

double d_growthrate_dt(struct cosm_model cosm,double z)
{
  double D1,Hubble,scale;
  double x,x0,D10,k;

  scale = 1.e0/(1.e0+z);
  k = cosm.omega_m+cosm.omega_v-1.e0;
  Hubble = H0*cosm.hubble*
    sqrt(cosm.omega_m/pow(scale,3.0)-k/pow(scale,2.0)+cosm.omega_v);
  if(cosm.omega_m>0.99 && cosm.omega_v<0.01) {
    D1 = scale*Hubble;
  }
  if(cosm.omega_m<1.e0 && cosm.omega_v<0.01) {
    x  = (1.e0/cosm.omega_m-1.e0)/(1.e0+z);
    x0 = (1.e0/cosm.omega_m-1.e0);
    D10 = 1.0 + 3.0/x0 
      + 3.0*sqrt((1.0+x0)/pow(x0,3.e0))*log(sqrt(1.0+x0)-sqrt(x0));
    D1 = -3.0/(x*x)
      +1.5*sqrt(x*x*x/(1.e0+x))*(-2.0*x-3.0)/(x*x*x*x)*log(sqrt(1.e0+x)-sqrt(x))
      +1.5*sqrt((1.0+x)/(x*x*x))*(1.0/sqrt(1.0+x)-1.0/sqrt(x))/(sqrt(1.0+x)-sqrt(x));
    D1 = D1*(1.e0/cosm.omega_m-1.e0)*scale*Hubble;
    D1 = D1/D10;
  }
  if(cosm.omega_m<1.e0 && cosm.omega_v>0.01){
    x  = pow(2.0*(1.0/cosm.omega_m-1.0),1.0/3.0)/(1.0+z);
    x0 = pow(2.0*(1.0/cosm.omega_m-1.0),1.0/3.0);
    D10 = (0.2*x0+0.044*pow(x0,1.37))/(1.0+0.328*pow(x0,1.37));
    D1 = ((0.2+0.044*1.37*pow(x,0.37))*(1.0+0.328*pow(x,1.37))
	  -(0.2*x+0.044*pow(x,1.37))*0.328*1.37*pow(x,0.37))
      /(1.0+0.328*pow(x,1.37))/(1.0+0.328*pow(x,1.37));
    D1 = D1*pow(2.0*(1.0/cosm.omega_m-1.0),1.0/3.0)*scale*Hubble;
    D1 = D1/D10;
  }
   
  return D1;
}

double power_spec(struct cosm_model cosm, double k)
{
  double power_spectrum;

  power_spectrum = cosm.pk_amplitude*pow(k,cosm.pkindex)*Transfer(cosm,k);

  return power_spectrum;
}

double delta_c(struct cosm_model cosm, double zvir)
{
  double delta,omegavir,etavir,wvir;
  double Dt;

  Dt = growthrate(cosm,zvir);

  if(cosm.omega_m>0.99 && cosm.omega_v<0.01) delta = 1.69;
  if(cosm.omega_m<1.e0 && cosm.omega_v<0.01){
    /*       omegavir = cosm.omega_m*pow(1.e0+zvir,3.0) */
    /* 	 /(cosm.omega_m*pow(1.e0+zvir,3.0) */
    /* 	 +(1.e0-cosm.omega_m-cosm.omega_v)*pow(1.e0+zvir,2.0)+cosm.omega_v); */
    omegavir = cosm.omega_m;
    etavir = acosh(2.0/omegavir - 1.0);
    delta = 
      1.5*((3.0*sinh(etavir)*(sinh(etavir)-etavir))/pow(cosh(etavir)-1.0,2.e0)
	   -2.0)*(1.0+pow(2*PI/(sinh(etavir)-etavir),0.6667));
  }
  if(cosm.omega_m<1.e0 && cosm.omega_v>0.01){
    /*       omegavir = cosm.omega_m*pow(1.e0+zvir,3.0) */
    /*          /(cosm.omega_m*pow(1.e0+zvir,3.0)  */
    /*          +(1.e0-cosm.omega_m-cosm.omega_v)*pow(1.e0+zvir,2.0)+cosm.omega_v); */
    omegavir = cosm.omega_m;
    wvir = 1.0/omegavir-1.0;
    delta = 3.e0*pow(12.0*PI,0.6667)/20.0*(1.0+0.0123*log10(omegavir));
  }

  return delta/Dt;
}

double Delta_c(struct cosm_model cosm, double zvir)
{
  double Delta;
  double omegavir,etavir,wvir;

  if(cosm.omega_m>0.99 && cosm.omega_v<0.01) Delta = 18.e0*PI*PI;
  if(cosm.omega_m<1.e0 && cosm.omega_v<0.01) {
    omegavir = cosm.omega_m*pow(1.e0+zvir,3.0)
      /(cosm.omega_m*pow(1.e0+zvir,3.0) 
	+(1.e0-cosm.omega_m-cosm.omega_v)*pow(1.e0+zvir,2.0)+cosm.omega_v);
    etavir = acosh(2.0/omegavir - 1.0);
    Delta = 4.0*PI*PI*pow(cosh(etavir)-1.0,3.e0)/pow(sinh(etavir)-etavir,2.e0);
  }
  if(cosm.omega_m<1.e0 && cosm.omega_v>0.01){
    omegavir = cosm.omega_m*pow(1.e0+zvir,3.0)
      /(cosm.omega_m*pow(1.e0+zvir,3.0) 
	+(1.e0-cosm.omega_m-cosm.omega_v)*pow(1.e0+zvir,2.0)+cosm.omega_v);
    wvir = 1.e0/omegavir-1.e0;
    Delta = 18.e0*PI*PI*(1.e0+0.4093*pow(wvir,0.9052));
  }

  return Delta;
}

double Transfer(struct cosm_model cosm, double k)
{
  double q,T,k_FS;
   
  if(cosm.DM==1){ /* CDM */
    q=k/(cosm.shapeparam*cosm.hubble);

    T = pow(log(1.e0+2.34*q)/(2.34*q),2.e0)
      /sqrt(1.e0+3.89*q+pow(16.1*q,2.e0)+pow(5.46*q,3.e0)+pow(6.71*q,4.e0));
  }else if(cosm.DM==2){ /* HDM with m=30eV */
    double k_FS = 0.16;
    T = exp(-4.61*pow(k/k_FS,1.5));
  }
   
  return T;
}

double sigmaR(struct cosm_model cosm,double r)
{
  double sigma;
  double kmin,kmax;
  double delta,kr,k;
  int    nbin,p;

  kmin = 1.e-5;kmax = 1.e3;
  nbin = 1000;

  delta = (log(kmax)-log(kmin))/(double)nbin;
  sigma = 0.e0;

  for(p=1;p<=nbin;p++){
    k     = kmin*exp(((float)p-0.5)*delta);
    kr    = k*r;
    sigma += Transfer(cosm,k)*pow(k,cosm.pkindex)*pow(k,3.e0)
      *9.e0*pow((sin(kr)-kr*cos(kr))/pow(kr,3.e0),2.e0);
  }

  sigma = sqrt(cosm.pk_amplitude*sigma*delta/(2.e0*PI*PI));
   
  return sigma;
}

double sigmaM(struct cosm_model cosm,double Mass)
{
  double sigma,rTH;

  rTH = pow((3.e0*Mass)/(4.e0*PI*cosm.omega_m*pow(cosm.hubble,2.e0)*RHOC_HI2),1.0/3.0);
  sigma = sigmaR(cosm,rTH);

  return sigma;
}

double psmassfunc(struct cosm_model cosm, double Mass, double z)
{
  double nps;

  double sigma,dsdm;
  double kmin,kmax;
  double delta_k,kr,k,r;
  int nbin,p;

  r = pow((3.e0*Mass)/(4.e0*PI*cosm.omega_m*pow(cosm.hubble,2.e0)*RHOC_HI2),1.0/3.0);

  kmin = 1.e-5;kmax = 1.e3;
  nbin = 1000;

  delta_k = (log(kmax)-log(kmin))/(double)nbin;
  sigma   = 0.e0;
  dsdm    = 0.e0;

  for(p=1;p<=nbin;p++){
    k     = kmin*exp(((float)p-0.5)*delta_k);
    kr    = k*r;
    sigma += Transfer(cosm,k)*pow(k,cosm.pkindex)*pow(k,3.e0)
      *9.e0*pow((sin(kr)-kr*cos(kr))/pow(kr,3.e0),2.e0);
    dsdm  += Transfer(cosm,k)*pow(k,cosm.pkindex)*pow(k,4.e0)
      *18.e0*((sin(kr)-kr*cos(kr))/pow(kr,3.e0)
	      *(pow(kr,2.e0)*sin(kr)-3.e0*(sin(kr)-kr*cos(kr)))/pow(kr,4.0));
  }

  sigma = sqrt(cosm.pk_amplitude*sigma*delta_k/(2.e0*PI*PI));
  /*   sigma = sigma*growthrate(cosm,z);*/
  dsdm  = cosm.pk_amplitude*dsdm*delta_k*0.5/pow(PI,2.0)*r/Mass/3.e0;
  /*    dsdm  = dsdm/2.e0/sigma*pow(growthrate(cosm,z),2.0); */
  dsdm  = dsdm/2.e0/sigma;

  nps   = sqrt(2.0/PI)*delta_c(cosm,z)*(-dsdm/(sigma*sigma))*RHOC_HI2/Mass
    *cosm.omega_m*(cosm.hubble*cosm.hubble)
    *exp(-(delta_c(cosm,z)*delta_c(cosm,z))/2.0/(sigma*sigma));

  return nps;
}
   

double cumulative_psmassfunc(struct cosm_model cosm,double Mass,double z)
{
  double cumps;

  double Mmax,M,delta_M;
  int nbin,p;

  Mmax = 1.e17;
  nbin = 500;
  delta_M = (log(Mmax)-log(Mass))/(double)nbin;

  cumps = 0.e0;

  for(p=1;p<=nbin;p++){
    M = Mass*exp(((float)p-0.5)*delta_M);
    cumps += psmassfunc(cosm,M,z)*M;
  }

  cumps *= delta_M;

  return cumps;
}

double d_kai(struct cosm_model cosm, double z)
{
  double delta_kai;
  double om,ov,zp1;

  om  = cosm.omega_m;
  ov  = cosm.omega_v;
  zp1 = 1.e0+z;

  delta_kai = 1.0/sqrt(om*zp1*zp1*zp1+(1.e0-om-ov)*zp1*zp1+ov);

  return delta_kai;
}

double kai_distance(struct cosm_model cosm,double z)
{
  double kai;
  double k0,sqk0;
  double ss,dss,h[21],s[21],eps;
  int    j,k = 5;

  if(z<1.e-5){
    kai = z;
    return kai;
  }

  k0   = cosm.omega_m+cosm.omega_v-1.0;
  sqk0 = sqrt(fabs(k0));
  eps  = 1.e-6;
   
  h[1] = 1.e0;

  for(j=1;j<=20;j++){
    s[j] = cosm_trapzd(d_kai,cosm,0.e0,z,j);
    if(j>=k){
      cosm_polint(&h[j-k],&s[j-k],k,0.e0,&ss,&dss);
      if(fabs(dss)<eps*fabs(ss)){
	kai = ss;
	return kai;
      }
    }
    s[j+1]=s[j];
    h[j+1]=0.25*h[j];
  }

  fprintf(stderr,"too many steps in kai_distance \n");
}

double comoving_distance(struct cosm_model cosm, double z)
{
  double ss,k0,sqk0;
  double r_z;

  ss = kai_distance(cosm,z);

  k0   = cosm.omega_m+cosm.omega_v-1.0;
  sqk0 = sqrt(fabs(k0));

  if(k0<0.0 && fabs(k0)>0.01) {
    r_z = sinh(sqk0*ss)/sqk0*CSPEED/(H0*cosm.hubble);
  }else if(fabs(k0)<0.01) {
    r_z = ss*(CSPEED)/(H0*cosm.hubble);
  }else if(k0>0.0 && fabs(k0)>0.01) {
    r_z = sin(sqk0*ss)/sqk0*CSPEED/(H0*cosm.hubble);
  }
  return r_z;
}

double angular_distance(struct cosm_model cosm, double z)
{
  double d_A;
  double k0,sqk0;
  double ss;

  k0   = cosm.omega_m+cosm.omega_v-1.0;
  sqk0 = sqrt(fabs(k0));
  ss   = kai_distance(cosm,z);

  if(k0<0.0 && fabs(k0)>0.01){
    d_A = sinh(sqk0*ss)/sqk0*CSPEED/(H0*cosm.hubble)/(1.0+z);
  }else if(fabs(k0)<0.01){
    d_A = ss*(CSPEED)/(H0*cosm.hubble)/(1.0+z);
  }else if(k0>0.0 && fabs(k0)>0.01) {
    d_A = sin(sqk0*ss)/sqk0*CSPEED/(H0*cosm.hubble)/(1.0+z);
  }
  return d_A;

}

double luminosity_distance(struct cosm_model cosm,double z)
{
  double d_L;

  d_L = angular_distance(cosm,z)*(1.e0+z)*(1.e0+z);

  return d_L;
}

double dVoverdz(struct cosm_model cosm,double z)
{
  double dvdz,dA;
  double oneplusz,oneplusz2,oneplusz3;
  double om,ov,h;

  om = cosm.omega_m;
  ov = cosm.omega_v;
  h  = cosm.hubble;

  oneplusz  = 1.e0+z;
  oneplusz2 = oneplusz*oneplusz;
  oneplusz3 = oneplusz2*oneplusz;

  dA   = angular_distance(cosm,z);
  dvdz = CSPEED/(H0*h)*oneplusz2*dA*dA/sqrt(om*oneplusz3+(1.e0-om-ov)*oneplusz2+ov);

  return dvdz;
}

double comoving_volume(struct cosm_model cosm, double z1, double z2)
{
  double ss,dss,h[21],s[21];
  int i,k = 5;
  double eps = 1.0e-6;

  h[1] = 1.e0;
  for(i=1;i<=20;i++){
    s[i] = cosm_trapzd(dVoverdz,cosm,z1,z2,i);
    if(i>=k){
      cosm_polint(&h[i-k],&s[i-k],k,0.e0,&ss,&dss);
      if(fabs(dss)<=eps*fabs(ss)) return 2.e0*ss;
    }
    s[i+1]=s[i];
    h[i+1]=0.25*h[i];
  }
  fprintf(stderr,"too many steps in comoving_volume \n");
}

#define FUNC(x,y) ((*func)(x,y))

double cosm_trapzd(double (*func)(struct cosm_model,double), struct cosm_model cosm, double a, double b, int n)
{
  double x,tnm,sum,del;
  static double s;
  int it,j;

  if (n == 1) {
    return (s=0.5*(b-a)*(FUNC(cosm,a)+FUNC(cosm,b)));
  } else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(cosm,x);
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
}
#undef FUNC

void cosm_polint(double *xa,double *ya,int n,double x,double *y,double *dy)
{
  int i,m,ns=1;
  double den,dif,dift,ho,hp,w;
  double *c,*d;
   
  double *ctmp = (double *)malloc(n*sizeof(double));
  double *dtmp = (double *)malloc(n*sizeof(double));
  c = ctmp-1;
  d = dtmp-1;
   
  dif=fabs(x-xa[1]);

  for (i=1;i<=n;i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[(ns--)];
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      den=ho-hp;
      if (den == 0.0){
	fprintf(stderr,"Error in routine cosm_polint\n");
	exit(1);
      }
      den=w/den;
      d[i]=hp*den;
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
   
  free((void *)ctmp);
  free((void *)dtmp);
}

double cosm_erfc(double x)
{
  double t,z,ans;
   
  z=fabs(x);
  t=1.0/(1.0+0.5*z);
  ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));
  return x >= 0.0 ? ans : 2.0-ans;
}   
