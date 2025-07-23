#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "hdf5.h"
#include "tng_data.h"
#include "libcosm/cosmology.h"

void read_cdm(struct tng_cdm *cdm, hid_t file, struct tng_header *head, int io_flag)
{
  hid_t grp_top, grp_p1;
  hid_t dataset;

  grp_top = H5Gopen(file, "/", H5P_DEFAULT);
  grp_p1 = H5Gopen(file, "PartType1", H5P_DEFAULT);

  float *data3d = (float *)malloc(sizeof(float)*head->np_local[1]*3);
  float *data1d = (float *)malloc(sizeof(float)*head->np_local[1]);
  uint64_t *idata = (uint64_t *)malloc(sizeof(uint64_t)*head->np_local[1]);

  if(io_flag & _DM_POS) {
    dataset = H5Dopen(grp_p1, "Coordinates", H5P_DEFAULT);
    H5Dread(dataset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data3d);
#pragma omp parallel for schedule(auto)
    for(int i=0;i<head->np_local[1];i++) {      
      cdm[i].xpos = data3d[3*i+0];
      cdm[i].ypos = data3d[3*i+1];
      cdm[i].zpos = data3d[3*i+2];
    }
    H5Dclose(dataset);
  }

  if(io_flag & _DM_VEL) {
    dataset = H5Dopen(grp_p1, "Velocities", H5P_DEFAULT);
    H5Dread(dataset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data3d);
#pragma omp parallel for schedule(auto)
    for(int i=0;i<head->np_local[1];i++) {
      cdm[i].xvel = data3d[3*i+0];
      cdm[i].yvel = data3d[3*i+1];
      cdm[i].zvel = data3d[3*i+2];
    }
    H5Dclose(dataset);
  }

  if(io_flag & _DM_DENS) {  
    dataset = H5Dopen(grp_p1, "SubfindDMDensity", H5P_DEFAULT);
    H5Dread(dataset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data1d);
#pragma omp parallel for schedule(auto)
    for(int i=0;i<head->np_local[1];i++) {
      cdm[i].dens = data3d[i];
    }
    H5Dclose(dataset);  
  }

  if(io_flag & _DM_ID) {  
    dataset = H5Dopen(grp_p1, "ParticleIDs", H5P_DEFAULT);
    H5Dread(dataset, H5T_STD_U64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, idata);
#pragma omp parallel for schedule(auto)
    for(int i=0;i<head->np_local[1];i++) {      
      cdm[i].ID = idata[i];
    }
    H5Dclose(dataset);      
  }
  
  H5Gclose(grp_p1);
  H5Gclose(grp_top);

  free(data3d);
  free(data1d);
  free(idata);
}

void read_gas_cell(struct tng_gas *cell, hid_t file, struct tng_header *head, int io_flag)
{
  hid_t grp_top, grp_p0;
  hid_t dataset;

  grp_top = H5Gopen(file, "/", H5P_DEFAULT);
  grp_p0 = H5Gopen(file, "PartType0", H5P_DEFAULT);

  float *data3d = (float *)malloc(sizeof(float)*head->np_local[0]*3);
  float *data1d = (float *)malloc(sizeof(float)*head->np_local[0]);
  float *data10d = (float *)malloc(sizeof(float)*head->np_local[0]*10);
  uint64_t *idata = (uint64_t *)malloc(sizeof(uint64_t)*head->np_local[0]);

  if(io_flag & _GAS_POS) {
    dataset = H5Dopen(grp_p0, "Coordinates", H5P_DEFAULT);
    H5Dread(dataset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data3d);
#pragma omp parallel for schedule(auto)
    for(int i=0;i<head->np_local[0];i++) {      
      cell[i].xpos = data3d[3*i+0];
      cell[i].ypos = data3d[3*i+1];
      cell[i].zpos = data3d[3*i+2];
    }
    H5Dclose(dataset);
  }

  if(io_flag & _GAS_VEL) {  
    dataset = H5Dopen(grp_p0, "Velocities", H5P_DEFAULT);
    H5Dread(dataset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data3d);
#pragma omp parallel for schedule(auto)
    for(int i=0;i<head->np_local[0];i++) {      
      cell[i].xvel = data3d[3*i+0];
      cell[i].yvel = data3d[3*i+1];
      cell[i].zvel = data3d[3*i+2];
    }
    H5Dclose(dataset);
  }

  if(io_flag & _GAS_DENS) {  
    dataset = H5Dopen(grp_p0, "Density", H5P_DEFAULT);
    H5Dread(dataset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data1d);
#pragma omp parallel for schedule(auto)
    for(int i=0;i<head->np_local[0];i++) {      
      cell[i].dens = data1d[i];
    }
    H5Dclose(dataset);
  }

  if(io_flag & _GAS_FELE) {
    dataset = H5Dopen(grp_p0, "ElectronAbundance", H5P_DEFAULT);
    H5Dread(dataset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data1d);
#pragma omp parallel for schedule(auto)
    for(int i=0;i<head->np_local[0];i++) {
      cell[i].felec = data1d[i];
    }
    H5Dclose(dataset);
  }

  if(io_flag & _GAS_UENE) {
    dataset = H5Dopen(grp_p0, "InternalEnergy", H5P_DEFAULT);
    H5Dread(dataset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data1d);
#pragma omp parallel for schedule(auto)
    for(int i=0;i<head->np_local[0];i++) {
      cell[i].uene = data1d[i];
    }
    H5Dclose(dataset);
  }

  if((io_flag & _GAS_UENE) && (io_flag & _GAS_FELE)) {
    //double conv_fact = SQR(head->vunit);
    double conv_fact = 1.0e10;
#pragma omp parallel for schedule(auto)
    for(int i=0;i<head->np_local[0];i++) {
      double meanmu = 4.0/(1.0+3.0*XHYDROG+4.0*XHYDROG*cell[i].felec)*MPROTON;
      cell[i].tmpr = cell[i].uene*(GAMMA_EOS-1.0)*meanmu*conv_fact/KBOLTZ;
    }
  }

  if(io_flag & _GAS_MASS) {
    dataset = H5Dopen(grp_p0, "Masses", H5P_DEFAULT);
    H5Dread(dataset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data1d);
#pragma omp parallel for schedule(auto)
    for(int i=0;i<head->np_local[0];i++) {
      cell[i].mass = data1d[i];
    }
    H5Dclose(dataset);
  }

  if(io_flag & _GAS_FHI) {
    dataset = H5Dopen(grp_p0, "NeutralHydrogenAbundance", H5P_DEFAULT);
    H5Dread(dataset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data1d);
#pragma omp parallel for schedule(auto)
    for(int i=0;i<head->np_local[0];i++) {
      cell[i].fHI = data1d[i];
    }
    H5Dclose(dataset);
  }

  if((io_flag & _GAS_DENS) && (io_flag & _GAS_MASS)) {
#pragma omp parallel for schedule(auto)
    for(int i=0;i<head->np_local[0];i++) {
      cell[i].vol = cell[i].mass/cell[i].dens;
    }
  }

  if(io_flag & _GAS_ID) {
    dataset = H5Dopen(grp_p0, "ParticleIDs", H5P_DEFAULT);
    H5Dread(dataset, H5T_STD_U64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, idata);
#pragma omp parallel for schedule(auto)
    for(int i=0;i<head->np_local[0];i++) {
      cell[i].ID = idata[i];
    }
    H5Dclose(dataset);
  }

  if(io_flag & _GAS_ZMETAL) {
    dataset = H5Dopen(grp_p0, "GFM_Metals", H5P_DEFAULT);
    H5Dread(dataset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data10d);
#pragma omp parallel for schedule(auto)
    for(int i=0;i<head->np_local[0];i++) {
      cell[i].zmetal[0] = data10d[10*i];
      cell[i].zmetal[1] = data10d[10*i+1];
      cell[i].zmetal[2] = data10d[10*i+2];
      cell[i].zmetal[3] = data10d[10*i+3];
      cell[i].zmetal[4] = data10d[10*i+4];
      cell[i].zmetal[5] = data10d[10*i+5];
      cell[i].zmetal[6] = data10d[10*i+6];
      cell[i].zmetal[7] = data10d[10*i+7];
      cell[i].zmetal[8] = data10d[10*i+8];
      cell[i].zmetal[9] = data10d[10*i+9];
    }
    H5Dclose(dataset);  
  }

  H5Gclose(grp_p0);
  H5Gclose(grp_top);


  free(idata);
  free(data10d);
  free(data1d);
  free(data3d);
}


void read_gas_cell_fast(struct tng_gas *cell, hid_t file, struct tng_header *head, int io_flag)
{
  hid_t grp_top, grp_p0;
  hid_t dataset;

  grp_top = H5Gopen(file, "/", H5P_DEFAULT);
  grp_p0 = H5Gopen(file, "PartType0", H5P_DEFAULT);

  float *data3d = (float *)malloc(sizeof(float)*head->np_local[0]*3);
  float *data1d = (float *)malloc(sizeof(float)*head->np_local[0]);
  float *data10d = (float *)malloc(sizeof(float)*head->np_local[0]*10);
  uint64_t *idata = (uint64_t *)malloc(sizeof(uint64_t)*head->np_local[0]);

#pragma omp sections
  {

#pragma omp section
    if(io_flag & _GAS_POS) {
      float *data3d = (float *)malloc(sizeof(float)*head->np_local[0]*3);
      dataset = H5Dopen(grp_p0, "Coordinates", H5P_DEFAULT);
      H5Dread(dataset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data3d);
      for(int i=0;i<head->np_local[0];i++) {      
	cell[i].xpos = data3d[3*i+0];
	cell[i].ypos = data3d[3*i+1];
	cell[i].zpos = data3d[3*i+2];
      }
      H5Dclose(dataset);
      free(data3d);
    }

#pragma omp section
    if(io_flag & _GAS_VEL) {
      float *data3d = (float *)malloc(sizeof(float)*head->np_local[0]*3);
      dataset = H5Dopen(grp_p0, "Velocities", H5P_DEFAULT);
      H5Dread(dataset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data3d);
      for(int i=0;i<head->np_local[0];i++) {
	cell[i].xvel = data3d[3*i+0];
	cell[i].yvel = data3d[3*i+1];
	cell[i].zvel = data3d[3*i+2];
      }
      H5Dclose(dataset);
      free(data3d);
    }
#pragma omp section
    if(io_flag & _GAS_DENS) {
      float *data1d = (float *)malloc(sizeof(float)*head->np_local[0]);
      dataset = H5Dopen(grp_p0, "Density", H5P_DEFAULT);
      H5Dread(dataset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data1d);
      for(int i=0;i<head->np_local[0];i++) {      
	cell[i].dens = data1d[i];
      }
      H5Dclose(dataset);
      free(data1d);
    }
#pragma omp section
    if(io_flag & _GAS_FELE) {
      float *data1d = (float *)malloc(sizeof(float)*head->np_local[0]);
      dataset = H5Dopen(grp_p0, "ElectronAbundance", H5P_DEFAULT);
      H5Dread(dataset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data1d);
      for(int i=0;i<head->np_local[0];i++) {
	cell[i].felec = data1d[i];
      }
      H5Dclose(dataset);
      free(data1d);
    }
#pragma omp section
    if(io_flag & _GAS_UENE) {
      float *data1d = (float *)malloc(sizeof(float)*head->np_local[0]);
      dataset = H5Dopen(grp_p0, "InternalEnergy", H5P_DEFAULT);
      H5Dread(dataset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data1d);
      for(int i=0;i<head->np_local[0];i++) {
	cell[i].uene = data1d[i];
      }
      H5Dclose(dataset);
      free(data1d);
    }
#pragma omp section
    if(io_flag & _GAS_MASS) {
      float *data1d = (float *)malloc(sizeof(float)*head->np_local[0]);
      dataset = H5Dopen(grp_p0, "Masses", H5P_DEFAULT);
      H5Dread(dataset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data1d);
      for(int i=0;i<head->np_local[0];i++) {
	cell[i].mass = data1d[i];
      }
      H5Dclose(dataset);
      free(data1d);
    }
#pragma omp section
    if(io_flag & _GAS_FHI) {
      float *data1d = (float *)malloc(sizeof(float)*head->np_local[0]);
      dataset = H5Dopen(grp_p0, "NeutralHydrogenAbundance", H5P_DEFAULT);
      H5Dread(dataset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data1d);
      for(int i=0;i<head->np_local[0];i++) {
	cell[i].fHI = data1d[i];
      }
      H5Dclose(dataset);
      free(data1d);
    }
#pragma omp section
    if(io_flag & _GAS_ID) {
      float *data1d = (float *)malloc(sizeof(float)*head->np_local[0]);
      dataset = H5Dopen(grp_p0, "ParticleIDs", H5P_DEFAULT);
      H5Dread(dataset, H5T_STD_U64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, idata);
      for(int i=0;i<head->np_local[0];i++) {
	cell[i].ID = idata[i];
      }
      H5Dclose(dataset);
      free(data1d);
    }
#pragma omp section
    if(io_flag & _GAS_ZMETAL) {
      float *data10d = (float *)malloc(sizeof(float)*head->np_local[0]*10);
      dataset = H5Dopen(grp_p0, "GFM_Metals", H5P_DEFAULT);
      H5Dread(dataset, H5T_IEEE_F32LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data10d);
      for(int i=0;i<head->np_local[0];i++) {
	cell[i].zmetal[0] = data10d[10*i];
	cell[i].zmetal[1] = data10d[10*i+1];
	cell[i].zmetal[2] = data10d[10*i+2];
	cell[i].zmetal[3] = data10d[10*i+3];
	cell[i].zmetal[4] = data10d[10*i+4];
	cell[i].zmetal[5] = data10d[10*i+5];
	cell[i].zmetal[6] = data10d[10*i+6];
	cell[i].zmetal[7] = data10d[10*i+7];
	cell[i].zmetal[8] = data10d[10*i+8];
	cell[i].zmetal[9] = data10d[10*i+9];
      }
      H5Dclose(dataset);
      free(data10d);
    }
    
  } // omp sections

  if((io_flag & _GAS_UENE) && (io_flag & _GAS_FELE)) {
    //double conv_fact = SQR(head->vunit);
    double conv_fact = 1.0e10;
    for(int i=0;i<head->np_local[0];i++) {
      double meanmu = 4.0/(1.0+3.0*XHYDROG+4.0*XHYDROG*cell[i].felec)*MPROTON;
      cell[i].tmpr = cell[i].uene*(GAMMA_EOS-1.0)*meanmu*conv_fact/KBOLTZ;
    }
  }

  if((io_flag & _GAS_DENS) && (io_flag & _GAS_MASS)) {
#pragma omp parallel for schedule(auto)
    for(int i=0;i<head->np_local[0];i++) {
      cell[i].vol = cell[i].mass/cell[i].dens;
    }
  }



  H5Gclose(grp_p0);
  H5Gclose(grp_top);


  free(idata);
  free(data10d);
  free(data1d);
  free(data3d);
}


void read_header(struct tng_header *header, hid_t file)
{
  hid_t grp_top, grp_header;
  hid_t attr;

  grp_top = H5Gopen(file, "/", H5P_DEFAULT);
  grp_header = H5Gopen(file, "Header", H5P_DEFAULT);

  attr = H5Aopen(grp_header, "BoxSize", H5P_DEFAULT);
  H5Aread(attr, H5T_IEEE_F64LE, &header->boxsize);
  H5Aclose(attr);

  static double mass_tbl[6];
  attr = H5Aopen(grp_header, "MassTable", H5P_DEFAULT);
  H5Aread(attr, H5T_IEEE_F64LE, mass_tbl);
  H5Aclose(attr);
  header->mcdm = mass_tbl[1];

  attr = H5Aopen(grp_header, "NumFilesPerSnapshot", H5P_DEFAULT);
  H5Aread(attr, H5T_STD_I32LE, &header->nchunk);
  H5Aclose(attr);

  attr = H5Aopen(grp_header, "HubbleParam", H5P_DEFAULT);  
  H5Aread(attr, H5T_IEEE_F64LE, &header->hubble);
  H5Aclose(attr);

  attr = H5Aopen(grp_header, "NumPart_ThisFile", H5P_DEFAULT);  
  H5Aread(attr, H5T_STD_I32LE, &header->np_local);
  H5Aclose(attr);

  static uint32_t np_total[6];
  static uint32_t hword[6];  
  attr = H5Aopen(grp_header, "NumPart_Total", H5P_DEFAULT);  
  H5Aread(attr, H5T_STD_U32LE, np_total);
  H5Aclose(attr);

  attr = H5Aopen(grp_header, "NumPart_Total_HighWord", H5P_DEFAULT);
  H5Aread(attr, H5T_STD_U32LE, hword);
  H5Aclose(attr);

  for(int i=0;i<6;i++) {
    header->np_total[i] = hword[i]*(uint64_t)(1L<<31) + np_total[i];
  }

  attr = H5Aopen(grp_header, "Redshift", H5P_DEFAULT);  
  H5Aread(attr, H5T_IEEE_F64LE, &header->zred);
  H5Aclose(attr);

  attr = H5Aopen(grp_header, "Time", H5P_DEFAULT);  
  H5Aread(attr, H5T_IEEE_F64LE, &header->tnow);
  H5Aclose(attr);

#if 0
  attr = H5Aopen(grp_header, "UnitLength_in_cm", H5P_DEFAULT);  
  H5Aread(attr, H5T_IEEE_F64LE, &header->lunit);
  H5Aclose(attr);    

  attr = H5Aopen(grp_header, "UnitMass_in_g", H5P_DEFAULT);  
  H5Aread(attr, H5T_IEEE_F64LE, &header->munit);
  H5Aclose(attr);    

  attr = H5Aopen(grp_header, "UnitVelocity_in_cm_per_s", H5P_DEFAULT);  
  H5Aread(attr, H5T_IEEE_F64LE, &header->vunit);
  H5Aclose(attr);
#else
  header->lunit = KPC/header->hubble;                    // comoving kpc/h
  header->munit = SOLARMASS*1.0e10/header->hubble;       // 10^10 Msolar/h
  header->vunit = 1.0e5*sqrt(1.0/(1.0+header->zred));    // sqrt(a)*km/sec
#endif

  H5Gclose(grp_header);
  H5Gclose(grp_top);
}

void show_header(struct tng_header *head)
{

  fprintf(stdout,
	  "# boxsize                            : %10.4e [cMpc/h]\n",
	  head->boxsize/1000.0);
  fprintf(stdout,
	  "# number of chunks for this snapshot : %d\n", head->nchunk);
  fprintf(stdout,
	  "# number of local CDM particles      : %d\n", head->np_local[1]);
  fprintf(stdout,
	  "# number of local gas cells          : %d\n", head->np_local[0]);
  fprintf(stdout,
	  "# number of total CDM particles      : %d\n", head->np_total[1]);
  fprintf(stdout,
	  "# number of total gas cells          : %d\n", head->np_total[0]);
  fprintf(stdout,
	  "# redshift                           : %10.4e\n", head->zred);
  fprintf(stdout,
	  "# Hubble parameter                   : %10.4e\n", head->hubble);
  fprintf(stdout,
	  "# lunit                              : %10.4e\n", head->lunit);
  fprintf(stdout,
	  "# munit                              : %10.4e\n", head->munit);
  fprintf(stdout,
	  "# vunit                              : %10.4e\n", head->vunit);  
  fflush(stdout);

}
