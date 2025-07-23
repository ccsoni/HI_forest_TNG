

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
  double vunit;  // UnitVelocity_in_cm_per_s
};

struct tng_gas {
  float xpos, ypos, zpos; // in units of [ckpc/h]
  float xvel, yvel, zvel; // in units of [km*sqrt(a)/s]
  float dens;             // in units of 10^10 Msun/h / (ckpc/h)^3
  float mass; // gas mass in units of 10^10 Msun
  float felec; // fraction of ne relative to total hydrogen number density
  float fHI;   // fraction of neutral hydrogen mass (density) 
  float uene; // specific thermal energy in units of  (km/s)^2
  float zmetal[10]; // metal abundances of H, He, C, N, O, Ne, Mg, Si, Fe and Z
  uint64_t ID; // particle ID
};

struct tng_cdm {
  float xpos, ypos, zpos; // in units of [ckpc/h]
  float xvel, yvel, zvel; // in units of [km*sqrt(a)/s]
  float dens; //  cdm density in units of 10^10 Msun/h / (ckpc/h)^3
  uint64_t ID; // particle ID
};
