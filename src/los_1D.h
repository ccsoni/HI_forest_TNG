#pragma once
#define NMESH_LOS (1024)
//#define NMESH_LOS (512)

#include "libcosm/cosmology.h"
#include "graph.h"

struct los_1D{ // line-of-sight profile in redshift space
  int count[NMESH_LOS];
  double tau[NMESH_LOS];
  double vel[NMESH_LOS];
  double tmpr[NMESH_LOS];
  double fHI[NMESH_LOS];
  double NHI[NMESH_LOS];
  double zred[NMESH_LOS];
  double dens[NMESH_LOS];
  //  double path_length[NMESH_LOS];
  double column_mass[NMESH_LOS];
  double zred_min, zred_max, dzred;
  double nu_min, nu_max, delta_nu;  // in units of Hz
};

struct line_segment{
  double delta_x, delta_y; //displacement between adjuscent line segments along x- and y-direction
  double length;  // legth between two end points
  double dist_from_origin; // distance between the near end of this line segment and the origin
  double zred_min, zred_max;
  struct point p0, p1;
};

#define HI21_A10 (2.86888e-15)
#define HI21_T10 (0.0681685261)
#define HI21_K (3.0*SQR(CSPEED)/(32.0*PI*SQR(NU_0))*HI21_A10*XHYDROG/MPROTON*HI21_T10)

void zero_out_los(struct los_1D *);
void setup_los(struct los_1D *, struct line_segment*, int, struct tng_header *, struct cosm_model);
void calc_los_tau(struct tng_gas *, struct los_1D *, struct line_segment *, int, struct cosm_model, struct tng_header *);
void output_los_data(struct los_1D *, char *, int);
void input_los_data_header(char *, int*);
void input_los_data(struct los_1D*, char *, int);
double kernel_line_integral(double, double);
void setup_line_segment(struct line_segment *, int, struct point, struct vector, struct tng_header*, struct cosm_model);
int check_line_segment_setting(struct point, struct vector, int, struct tng_header *);
