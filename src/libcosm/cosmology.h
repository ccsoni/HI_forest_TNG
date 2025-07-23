/*
------------------------------------------------------------
        Copyright (c) Kohji Yoshikawa, 1999

  All Rights Reserved. Unpublished rights reserved under the
  copyright laws of Japan
------------------------------------------------------------
HISTORY

  last updated 99/06/23

  Cosmology.h
------------------------------------------------------------
*/

#ifdef __cplusplus
extern "C"
{
#endif

#ifndef _COSMOLOGY_
#define _COSMOLOGY_

#define RHOC_HI2  (2.77536627e11)   /* in units of Msolar/Mpc^3 */
#define T_CMB     (2.73)            /* in units of Kelvin */
#define MPC       (3.086e24)        /* in units of cm */
#define KPC       (3.086e21)        /* in units of cm */
#define SOLARMASS (1.989e+33)       /* in units of g */
#define CSPEED    (2.99792458e10)   /* in units of c.g.s */
#define H0        (3.240779e-18)    /* in units of c.g.s */
#define GNEWTON   (6.67259e-8)	    /* in units of c.g.s */

struct cosm_model{
  double omega_m, omega_b, omega_v, hubble;
  double sigma8, bias, shapeparam;
  double pk_amplitude, cobe_amplitude;
  double pkindex;
  
  int DM; /* 1 for CDM, 2 for HDM with m=30eV */
};

#endif /* _COSMOLOGY_ */

#ifndef PI
#define PI (3.14159265359)
#endif /* PI */

extern void init_cosm(struct cosm_model*, double, double, double, double, double, double, int);
extern double growthrate(struct cosm_model,double);
extern double d_growthrate_dt(struct cosm_model,double);
extern double power_spec(struct cosm_model,double);
extern double delta_c(struct cosm_model,double);
extern double Delta_c(struct cosm_model,double);
extern double Transfer(struct cosm_model,double);
extern double sigmaR(struct cosm_model,double);
extern double sigmaM(struct cosm_model,double);
extern double psmassfunc(struct cosm_model,double,double);
extern double cumulative_psmassfunc(struct cosm_model,double,double);
extern double d_kai(struct cosm_model,double);
extern double kai_distance(struct cosm_model,double);
extern double comoving_distance(struct cosm_model,double);
extern double angular_distance(struct cosm_model,double);
extern double luminosity_distance(struct cosm_model,double);
extern double dVoverdz(struct cosm_model,double);
extern double comoving_volume(struct cosm_model,double,double);

extern double cosm_trapzd(double (*func)(struct cosm_model,double), struct cosm_model,double,double,int);
extern void   cosm_polint(double*,double*,int,double,double*,double*);
extern double cosm_erfc(double);

#ifdef __cplusplus
}
#endif
