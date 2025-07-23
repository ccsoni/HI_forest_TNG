#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "tng_data.h"
#include "los_1D.h"
#include "libcosm/cosmology.h"

FILE *fp_pnt;

// return the redshift of the position correspoinding to (x,y,z) = (*,*,zpos)
// zred_max is the possible maximum redshift required to use the bisection method
double redshift_space(double zpos, struct tng_header *header,
		      struct cosm_model cosm, double zred_max)
{
  //zpos is in the simulation unit [ckpc/h]
  double zred_min = header->zred;
  double zpos_in_mpc = zpos/header->hubble/1.0e3;

  // find the zred using bisection-method
  double zred_mid = 0.5*(zred_max+zred_min);

  double rs_min = comoving_distance(cosm, zred_min)/MPC;
  double rs_max = comoving_distance(cosm, zred_max)/MPC;

  do {
    double rs_mid = comoving_distance(cosm, zred_mid)/MPC;
    if(rs_mid > rs_min + zpos_in_mpc) {
      zred_max = zred_mid;
      zred_mid = 0.5*(zred_mid+zred_min);
    }else{
      zred_min = zred_mid;
      zred_mid = 0.5*(zred_max+zred_mid);
    }
  }while(fabs(zred_max-zred_min)/zred_min > 1.0e-6);
  
  double zred_true = 0.5*(zred_max+zred_min);

  return zred_true;
  
}

void setup_line_segment(struct line_segment *seg, int nseg,
			struct point start_pnt, struct vector lseg,
			struct tng_header *header, struct cosm_model cosm)
{

  struct surface surfx_up, surfy_up, surfz_up;
  struct surface surfx_lo, surfy_lo, surfz_lo;

  // x=35000 ckpc
  surfx_up.cent.xpos = header->boxsize;
  surfx_up.cent.ypos = 0.0;
  surfx_up.cent.zpos = 0.0;
  surfx_up.norm_vect.x = 1.0;
  surfx_up.norm_vect.y = 0.0;
  surfx_up.norm_vect.z = 0.0;

  // x=0 ckpc
  surfx_lo.cent.xpos = 0.0;
  surfx_lo.cent.ypos = 0.0;
  surfx_lo.cent.zpos = 0.0;
  surfx_lo.norm_vect.x = 1.0;
  surfx_lo.norm_vect.y = 0.0;
  surfx_lo.norm_vect.z = 0.0;

  // y=35000 ckpc
  surfy_up.cent.xpos = 0.0;
  surfy_up.cent.ypos = header->boxsize;
  surfy_up.cent.zpos = 0.0;
  surfy_up.norm_vect.x = 0.0;
  surfy_up.norm_vect.y = 1.0;
  surfy_up.norm_vect.z = 0.0;  

  // y=0 ckpc
  surfy_lo.cent.xpos = 0.0;
  surfy_lo.cent.ypos = 0.0;
  surfy_lo.cent.zpos = 0.0;
  surfy_lo.norm_vect.x = 0.0;
  surfy_lo.norm_vect.y = 1.0;
  surfy_lo.norm_vect.z = 0.0;  

  // z=35000 ckpc 
  surfz_up.cent.xpos = 0.0;
  surfz_up.cent.ypos = 0.0;
  surfz_up.cent.zpos = header->boxsize;
  surfz_up.norm_vect.x = 0.0;
  surfz_up.norm_vect.y = 0.0;
  surfz_up.norm_vect.z = 1.0;

  // z=0 ckpc 
  surfz_lo.cent.xpos = 0.0;
  surfz_lo.cent.ypos = 0.0;
  surfz_lo.cent.zpos = 0.0;
  surfz_lo.norm_vect.x = 0.0;
  surfz_lo.norm_vect.y = 0.0;
  surfz_lo.norm_vect.z = 1.0;
  
  //  seg = (struct line_segment *) malloc(sizeof(struct line_segment)*nseg);

  double distance_from_origin = 0.0;  
  seg[0].p0 = start_pnt;
  seg[0].zred_min = header->zred;
  seg[0].dist_from_origin = 0.0;

  double x_abs, y_abs, z_abs;

  x_abs = fabs(lseg.x);
  y_abs = fabs(lseg.y);
  z_abs = fabs(lseg.z);

  double maximum_comp = fmax(x_abs, fmax(y_abs, z_abs));

  // surface through which the line segment first penetrate
  struct surface surf;  
#if 0  // this code implicitly assumes that maximum_comp has the positive direction.
  if(maximum_comp == x_abs) {
    surf = surfx_up;
  }else if(maximum_comp == y_abs) {
    surf = surfy_up;
  }else{
    surf = surfz_up;
  }
#else
  if(maximum_comp == x_abs) {
    if(lseg.x > 0.0) surf = surfx_up;
    else surf = surfx_lo;
  }else if(maximum_comp == y_abs) {
    if(lseg.y > 0.0) surf = surfy_up;
    else surf = surfy_lo;
  }else{
    if(lseg.z > 0.0) surf = surfz_up;
    else surf = surfz_lo;
  }
#endif

  for(int il=0;il<nseg;il++)  {
    struct line ll;
    ll.pnt =  seg[il].p0;
    ll.vect = lseg;

    double affine = compute_crosspoint_line_surface(&seg[il].p1, &ll, &surf);

    struct vector v01 = relative_vector(&seg[il].p0, &seg[il].p1);

    seg[il].delta_x = fabs(v01.x);
    seg[il].delta_y = fabs(v01.y);
    seg[il].length = sqrt(dot_product(&v01, &v01));

    distance_from_origin += seg[il].length;

    seg[il].zred_max = redshift_space(distance_from_origin,
				      header, cosm, header->zred+0.5);

    printf("# zred = %14.6e -- %14.6e\n", seg[il].zred_min, seg[il].zred_max);

    if(il+1<nseg) {
      seg[il+1].dist_from_origin = distance_from_origin;
      seg[il+1].zred_min = seg[il].zred_max;

      seg[il+1].p0.xpos = fmod(seg[il].p1.xpos, header->boxsize);
      seg[il+1].p0.ypos = fmod(seg[il].p1.ypos, header->boxsize);
      seg[il+1].p0.zpos = fmod(seg[il].p1.zpos, header->boxsize);
    }

  }
  
}

int check_line_segment_setting(struct point start_pnt, struct vector lseg,
			       int nseg, struct tng_header *header)
{
  double x_abs, y_abs, z_abs;

  x_abs = fabs(lseg.x);
  y_abs = fabs(lseg.y);
  z_abs = fabs(lseg.z);

  double maximum_comp = fmax(x_abs, fmax(y_abs, z_abs));

  int end_pnt_location;  // 0 for inside [0:boxsize)x[0:boxsize) otherwise 1
  int start_pnt_flag = 0;

  if (x_abs == maximum_comp) {

    if(lseg.x>0.0 && start_pnt.xpos!=0.0) start_pnt_flag |= 1;
    if(lseg.x<0.0 && start_pnt.xpos!=header->boxsize) start_pnt_flag |= 0b10;
    
    //double yend = start_pnt.ypos + nseg*lseg.y/lseg.x*header->boxsize;
    //double zend = start_pnt.zpos + nseg*lseg.z/lseg.x*header->boxsize;
    double yend = start_pnt.ypos + nseg*lseg.y/x_abs*header->boxsize;
    double zend = start_pnt.zpos + nseg*lseg.z/x_abs*header->boxsize;

    if((yend>0.0 && yend<header->boxsize) &&
       (zend>0.0 && zend<header->boxsize)) end_pnt_location = 0;
    else end_pnt_location = 1;
    
  }else if(y_abs == maximum_comp) {

    if(lseg.y>0.0 && start_pnt.ypos!=0.0) start_pnt_flag |= 0b100;
    if(lseg.y<0.0 && start_pnt.ypos!=header->boxsize) start_pnt_flag |= 0b1000;

    //double xend = start_pnt.xpos + nseg*lseg.x/lseg.y*header->boxsize;
    //double zend = start_pnt.zpos + nseg*lseg.z/lseg.y*header->boxsize;
    double xend = start_pnt.xpos + nseg*lseg.x/y_abs*header->boxsize;
    double zend = start_pnt.zpos + nseg*lseg.z/y_abs*header->boxsize;

    if((xend>0.0 && xend<header->boxsize) &&
       (zend>0.0 && zend<header->boxsize)) end_pnt_location = 0;
    else end_pnt_location = 1;

  }else{
    
    if(lseg.z>0.0 && start_pnt.zpos!=0.0) start_pnt_flag |= 0b10000;
    if(lseg.z<0.0 && start_pnt.zpos!=header->boxsize) start_pnt_flag |= 0b100000;

    //double xend = start_pnt.xpos + nseg*lseg.x/lseg.z*header->boxsize;
    //double yend = start_pnt.ypos + nseg*lseg.y/lseg.z*header->boxsize;
    double xend = start_pnt.xpos + nseg*lseg.x/z_abs*header->boxsize;
    double yend = start_pnt.ypos + nseg*lseg.y/z_abs*header->boxsize;
    
    if((xend>0.0 && xend<header->boxsize) &&
       (yend>0.0 && yend<header->boxsize)) end_pnt_location = 0;
    else end_pnt_location = 1;

  }

  return (end_pnt_location+start_pnt_flag);
}

void input_los_data_header(char *filename, int *nlos)
{
  FILE *fp = fopen(filename, "r");

  fread(nlos, sizeof(int), 1, fp);

  fclose(fp);
}

void input_los_data(struct los_1D *los, char *filename, int nlos)
{
  FILE *fp = fopen(filename, "r");

  int nlos_;
  
  fread(&nlos_, sizeof(int), 1, fp);

  if(nlos<nlos_) {
    fprintf(stderr,"nlos is too small in input_los_data()\n");
    exit(EXIT_FAILURE);
  }

  fread(los, sizeof(struct los_1D), nlos_, fp);

  printf("# %14.6e\n",los->nu_min);
  printf("# %14.6e\n",los->zred_min);
  for(int32_t inu=0;inu<NMESH_LOS;inu++) {
    printf("%14.6e %14.6e\n",los->nu_min + (inu+0.5)*los->delta_nu, los->tau[inu]);
  }
			    

  fclose(fp);
}

void output_los_data(struct los_1D *los, char *filename, int nlos)
{
  FILE *fp = fopen(filename, "w");

  fwrite(&nlos, sizeof(int), 1, fp);
  fwrite(los, sizeof(struct los_1D), nlos, fp);

  fclose(fp);

}

void zero_out_los(struct los_1D *los)
{
  for(int iz=0;iz<NMESH_LOS;iz++) {
    los->count[iz]=0.0;
    los->column_mass[iz]=0.0;
    //    los->path_length[iz]=0.0;
    los->tau[iz]=0.0;
    los->vel[iz]=0.0;
    los->tmpr[iz]=0.0;
    los->fHI[iz]=0.0;
    los->NHI[iz]=0.0;
    los->dens[iz]=0.0;
  }
}

double profile(double nu, double dnu, double T, double zred) // return averaged profile nu-dnu/2 < \nu < nu+dnu/2
{

  // frequency width in units of Hz
  double thermal_dnu = NU_0*sqrt(2.0*KBOLTZ*T/(MPROTON*SQR(CSPEED)));

  long double x_lo = (nu-0.5*dnu-NU_0/(1.0+zred))/(thermal_dnu/(1.0+zred));
  long double x_hi = (nu+0.5*dnu-NU_0/(1.0+zred))/(thermal_dnu/(1.0+zred));

  double phi = 0.5/(1.0+zred)*(erfl(x_hi)-erfl(x_lo))/dnu;

  //double phi = 1.0/(sqrt(PI)*thermal_dnu)*exp(-SQR((nu*(1.0+zred)-NU_0)/thermal_dnu));

  return phi;

}

void accum_tau(struct los_1D *los, double column, double fHI, double tmpr,
	       double zred, struct tng_header *header, int il)
{
  //  double dens_conv = header->munit/CUBE(header->lunit);

  // density * neutral fraction  
  //  double dens = cell->dens*CUBE(1.0+zred)*cell->fHI*dens_conv;

  // convert to physical path length in cgs unit  
  //  path *= (1.0/(1.0+header->zred)*header->lunit);
  
  double delta_nu = NU_0*sqrt(2.0*KBOLTZ*tmpr/(MPROTON*SQR(CSPEED)));
  if(il==2) fprintf(fp_pnt," %14.6e \n", HI21_K*column*fHI/fmax(tmpr,T_CMB*(1.0+zred))/delta_nu);

  for(int inu=0;inu<NMESH_LOS;inu++) {
    double nu = (los->nu_min + (inu+0.5)*los->delta_nu); // in units of Hz
    los->tau[inu] += column*fHI*profile(nu, los->delta_nu, tmpr, zred)/fmax(tmpr, T_CMB*(1.0+zred));
  }
  
}

void setup_los(struct los_1D *los, struct line_segment *seg, int nseg,
	       struct tng_header *header, struct cosm_model cosm)
{
  los->zred_min = seg[0].zred_min;
  los->zred_max = seg[nseg-1].zred_max;
  los->dzred = (los->zred_max-los->zred_min)/NMESH_LOS;

  los->nu_min = NU_0/(1.0 + los->zred_max);
  los->nu_max = NU_0/(1.0 + los->zred_min);
  los->delta_nu = (NU_0/(1.0 + los->zred_min) - NU_0/(1.0 + los->zred_max))/NMESH_LOS;

  printf("# zred_min = %12.4e : zred_max = %12.4e\n",
	 los->zred_min, los->zred_max);
  printf("# delta_nu = %12.4e [Hz]\n", los->delta_nu);
}

double kernel_line_integral(double b_impact, double hsm)
{
  double d  = b_impact/hsm;
  double d2 = SQR(d);
  double d4 = SQR(d2);

  double column;

  double zo = sqrt(1.0-d2);
  double zi = sqrt(1.0-4.0*d2);

  if(d<1.0e-4) {
    column = 6.0/PI;
  }  else if(1.0e-4<d && d<0.5){
    column = 32.0/PI*((((-3.0*d4)-12.0*d2)*asinh(sqrt(1.0-d2)/d)+sqrt(1.0-d2)*(13.0*d2+2.0))/8.0-
		      (((-24.0*d4)-96.0*d2)*asinh((sqrt(1.0-4*d2)/d)/2.0)+sqrt(1.0-4.0*d2)*(58.0*d2+15.0))/64.0)
      -(sqrt(1.0-4.0*d2)*(46.0*d2-11.0)-72.0*d4*asinh((sqrt(1.0-4*d2))/(2.0*d)))/(2.0*PI);
    //    column = 32.0/PI*((((-3.0*d4)-12.0*d2)*asinh(zo/d)+zo*(13.0*d2+2.0))/8.0-
    //    		      (((-24.0*d4)-96.0*d2)*asinh((zi/d)/2.0)+zi*(58.0*d2+15.0))/64.0)
    //      -(zi*(46.0*d2-11.0)-72.0*d4*asinh(zi/(2.0*d)))/(2.0*PI);
  }else if(d>0.5 && d<1.0) {
    column = 4.0/PI*(((-3.0*d4)-12.0*d2)*asinh(sqrt(1.0-d2)/d)+sqrt(1.0-d2)*(13.0*d2+2.0));
    //column = 4.0/PI*(((-3.0*d4)-12.0*d2)*asinh(zo/d)+zo*(13.0*d2+2.0));
  }else{
    column = 0.0;
  }

  return column;
}


void calc_los_tau(struct tng_gas *cell, struct los_1D *los,
		  struct line_segment *seg, int nseg,
		  struct cosm_model cosm, struct tng_header *header)
{
  struct timespec start, end;
  double column_unit_conv
    = header->munit/SQR(header->lunit)*SQR(1.0+header->zred);

  clock_gettime(CLOCK_MONOTONIC, &start);

//#pragma omp parallel for schedule(auto)  
  for(int i=0;i<header->np_local[0];i++) {

#if 1
    if(i%100000 == 0) {
      clock_gettime(CLOCK_MONOTONIC, &end);
      float elapsed_time = end.tv_sec-start.tv_sec+1.0e-9*(end.tv_nsec-start.tv_nsec);
      float rate = (float)i/elapsed_time;
      printf("\33[2K\r %8d / %d (%.2f \%) : %.1f cells/sec : ETR = %.1f [sec]",
             i, header->np_local[0],
             (float)i/header->np_local[0]*100.0, rate,
             (float)(header->np_local[0]-i)/rate);
      fflush(stdout);
    }
#endif

    double vol = cell[i].mass/cell[i].dens;   // comoving volume in (ckpc/h)^3
    double rad = cbrt(3.0*vol/(4.0*PI));      // comoving radius in ckpc/h
    double rad2 = rad*rad;
    double meanmu = 4.0/(1.0+3.0*XHYDROG+4.0*XHYDROG*cell[i].felec)*MPROTON;
    double thermal_dnu = NU_0*sqrt(2.0*KBOLTZ*cell[i].tmpr/(MPROTON*SQR(CSPEED)));

#pragma omp parallel for schedule(auto)
    for(int iseg=0;iseg<nseg;iseg++) {
      struct point pcell;
      pcell.xpos = cell[i].xpos;
      pcell.ypos = cell[i].ypos;
      pcell.zpos = cell[i].zpos;
      
      struct line lseg = connecting_line(&seg[iseg].p0, &seg[iseg].p1);

      // comoving distance in units of kpc/h between gas cell and line segments
      double dist_to_line = compute_distance_point_line(&pcell, &lseg);

      // comoving length in unit of kpc/h of the segment through which the LOS intersects with the kernel sphere
      double path = 2.0*sqrt(rad2-SQR(dist_to_line));

      if(rad > dist_to_line) { // if the line segments intersects the gas cell
	double dist_to_orig = seg[iseg].dist_from_origin
	  + sqrt(NORM2(pcell.xpos-seg[iseg].p0.xpos,
		       pcell.ypos-seg[iseg].p0.ypos,
		       pcell.zpos-seg[iseg].p0.zpos));

	double zred = redshift_space(dist_to_orig, header, cosm, los->zred_max);

	// dot producet of gas cell's velocity and unit vector along the LOS.
	double vpec = cell[i].xvel*lseg.vect.x
	  + cell[i].yvel*lseg.vect.y + cell[i].zvel*lseg.vect.z;
	vpec *= header->vunit;
#if 1	
	zred += vpec/CSPEED*(1.0+header->zred);
#endif

	int iz = (int)(floor((zred-los->zred_min)/los->dzred));

	// physical column density \rho * (path length) in cgs unit
#if 0
	double W = kernel_line_integral(dist_to_line, rad);
	double column
	  = W*cell[i].mass/rad2*header->munit/SQR(header->lunit)*SQR(1.0+header->zred);
#else
	double column = cell[i].dens*path*header->munit/SQR(header->lunit)*SQR(1.0+header->zred);
#endif

#pragma omp critical
	{
	  if(iz >= 0 && iz < NMESH_LOS) {
	    los->count[iz] += 1;
	    los->column_mass[iz] += column;
	    //	  los->path_length[iz] += path;
	    los->NHI[iz] += column*cell[i].fHI;
	    los->vel[iz] += column*vpec;
	    los->tmpr[iz] += column*cell[i].tmpr;
	    los->fHI[iz] += column*cell[i].fHI;
	    los->dens[iz] += cell[i].dens;
	  }

	
	  for(int inu=0;inu<NMESH_LOS;inu++) {
	    double nu_lo = los->nu_min + inu*los->delta_nu;
	    double nu_hi = los->nu_min + (inu+1)*los->delta_nu;

	    long double x_lo = (nu_lo-NU_0/(1.0+zred))/(thermal_dnu/(1.0+zred));
	    long double x_hi = (nu_hi-NU_0/(1.0+zred))/(thermal_dnu/(1.0+zred));

	    double phi = 0.5*(erfl(x_hi)-erfl(x_lo))/(los->delta_nu*(1.0+zred));
	    los->tau[inu] += column*cell[i].fHI*phi/fmax(cell[i].tmpr, T_CMB*(1.0+zred));	  
	  }
	}
      }
    }
    
  }

#if 1
  clock_gettime(CLOCK_MONOTONIC, &end);
  printf("\33[2K\r# elapsed time: %12.4e [sec]\n",
         end.tv_sec-start.tv_sec+1.0e-9*(end.tv_nsec-start.tv_nsec));
  fflush(stdout);
#endif
  
}
