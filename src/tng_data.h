#pragma once
#include "hdf5.h"

#define _DM_POS     (1)
#define _DM_VEL     (1<<1)
#define _DM_DENS    (1<<2)
#define _DM_ID      (1<<3)

#define _GAS_POS    (1)
#define _GAS_VEL    (1<<1)
#define _GAS_DENS   (1<<2)
#define _GAS_MASS   (1<<3)
#define _GAS_FELE   (1<<4)
#define _GAS_FHI    (1<<5)
#define _GAS_UENE   (1<<6)
#define _GAS_ZMETAL (1<<7)
#define _GAS_ID     (1<<8)

struct tng_header {
  double boxsize;  // BoxSize in units of ckpc/h
  double hubble;   // HubbleParam
  int32_t np_local[6]; // NumPart_ThisFile
  uint64_t np_total[6]; // NumPart_Total
  int32_t nchunk; // Number of chunk files for this snapshot
  double mcdm;   // mass of a single CDM particles in units of 10^10 Msun/h
  double omegam; // Omega0
  double omegab; // OmegaBaryon
  double omegav; // OmegaLambda
  double zred;   // Redshift
  double tnow;   // Time
  double lunit;  // UnitLength_in_cm
  double munit;  // UnitMass_in_g
  double vunit;  // UnitVelocity_in_cm_per_s for peculiar velocity
};

struct tng_gas {
  float xpos, ypos, zpos; // in units of [ckpc/h]
  float xvel, yvel, zvel; // in units of [km*sqrt(a)/s]
  float dens;             // in units of (10^10 Msun/h) / ((ckpc/h)^3)
  float vol;              // defined as mass/dens [(ckpc/h)^3]
  float mass; // gas mass in units of 10^10 Msun/h
  float felec; // fraction of ne relative to total hydrogen number density
  float fHI;   // fraction of neutral hydrogen mass (density) 
  float uene; // specific thermal energy in units of  (km/s)^2
  float tmpr; // temperature in units of Kelvin
  float zmetal[10]; // metal abundances of H, He, C, N, O, Ne, Mg, Si, Fe and Z
  uint64_t ID; // particle ID
};

struct tng_cdm {
  float xpos, ypos, zpos; // in units of [ckpc/h]
  float xvel, yvel, zvel; // in units of [km*sqrt(a)/s]
  float dens; //  cdm density in units of 10^10 Msun/h / (ckpc/h)^3
  uint64_t ID; // particle ID
};

void read_cdm(struct tng_cdm *, hid_t, struct tng_header *, int);
void read_gas_cell(struct tng_gas *, hid_t, struct tng_header *, int);
void read_gas_cell_fast(struct tng_gas *, hid_t, struct tng_header *, int);
void read_header(struct tng_header *, hid_t);
void show_header(struct tng_header *);

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUART(x) ((x)*(x)*(x)*(x))
#define NORM2(x,y,z) (SQR(x)+SQR(y)+SQR(z))
#define MIN(x,y) ((x)>(y) ? (y) : (x))
#define MAX(x,y) ((x)<(y) ? (y) : (x))

#define NU_0 (1.4204057e9) // rest frequency of HI 21 cm line in units of Hz
#define KBOLTZ (1.3806488e-16) // Boltzmann constant in the cgs unit
#define XHYDROG (0.755)
#define MPROTON (1.672623e-24) // proton mass in the cgs unit
#define GAMMA_EOS  (5.0/3.0)
#define PI (3.14159265359)
