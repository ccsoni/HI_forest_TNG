#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tng_data.h"
#include "map_2D.h"

#define __MASS_WEIGHTED_MAP__

void zero_out_map(struct map_2D *_map)
{

#pragma omp parallel for schedule(auto)
  for(int ix=0;ix<NMESH_MAP;ix++) {
    for(int iy=0;iy<NMESH_MAP;iy++) {
      _map->cdm_dens[ix][iy] = 0.0;
      _map->gas_dens[ix][iy] = 0.0;
      _map->tmpr[ix][iy] = 0.0;
      _map->fHI[ix][iy] = 0.0;
      _map->vol[ix][iy] = 0.0;
    }
  }
}

void output_map(char *filename, struct map_2D *_map, struct tng_header *header)
{
  FILE *fp;
  fp = fopen(filename, "w");

  float delta_x;
  delta_x = header->boxsize/NMESH_MAP;

  for(int ix=0;ix<NMESH_MAP;ix++) {
    for(int iy=0;iy<NMESH_MAP;iy++) {
      fprintf(fp, "%12.4e %12.4e %12.4e %12.4e %12.4e %12.4e \n",
	      ix*delta_x, iy*delta_x,
	      _map->cdm_dens[ix][iy], _map->gas_dens[ix][iy],
	      _map->tmpr[ix][iy],_map->fHI[ix][iy]);
    }
  }
}

void accumulate_cdm_to_map(struct map_2D *_map, struct tng_cdm *cdm,
			   struct tng_header *header)
{
  float delta_x;

  delta_x = header->boxsize/NMESH_MAP;

  for(int i=0;i<header->np_local[1];i++) {
    if(fabs(cdm[i].zpos-_map->slice_pos) < 0.5*delta_x) {
      int ix = (int)(cdm[i].xpos/delta_x);
      int iy = (int)(cdm[i].ypos/delta_x);

      _map->cdm_dens[ix][iy] += header->mcdm;
    }
  }
}

void accumulate_gas_to_map(struct map_2D *_map, struct tng_gas *cell,
			   struct tng_header *header)
{
  float delta_x;

  delta_x = header->boxsize/NMESH_MAP;

  for(int i=0;i<header->np_local[0];i++) {
    if(fabs(cell[i].zpos-_map->slice_pos) < 0.5*delta_x) {
      int ix = (int)(cell[i].xpos/delta_x);
      int iy = (int)(cell[i].ypos/delta_x);

#ifdef __MASS_WEIGHTED_MAP__
      _map->gas_dens[ix][iy] += cell[i].mass;
      _map->tmpr[ix][iy] += cell[i].tmpr*cell[i].mass;
      _map->fHI[ix][iy] += cell[i].fHI*cell[i].mass;
#else  // volume weighted
      _map->vol[ix][iy] += cell[i].vol;
      _map->gas_dens[ix][iy] += cell[i].mass;
      _map->tmpr[ix][iy] += cell[i].tmpr*cell[i].vol;
      _map->fHI[ix][iy] += cell[i].fHI*cell[i].vol;
#endif
    }
      
  }
  
}

void normalize_map(struct map_2D *_map, struct tng_header *header)
{

  double delta_x = header->boxsize/NMESH_MAP;
  double vol = CUBE(delta_x);

#pragma omp parallel for schedule(auto)    
  for(int ix=0;ix<NMESH_MAP;ix++) {
    for(int iy=0;iy<NMESH_MAP;iy++) {
#ifdef __MASS_WEIGHTED_MAP__
      _map->tmpr[ix][iy] /= (_map->gas_dens[ix][iy]+1.0e-30);
      _map->fHI[ix][iy] /= (_map->gas_dens[ix][iy]+1.0e-30);  
      _map->cdm_dens[ix][iy] /= vol;
      _map->gas_dens[ix][iy] /= vol;
#else // volume weighted
      _map->tmpr[ix][iy] /= (_map->vol[ix][iy]+1.0e-30);
      _map->fHI[ix][iy] /= (_map->vol[ix][iy]+1.0e-30);
      _map->cdm_dens[ix][iy] /= vol;
      _map->gas_dens[ix][iy] /= (_map->vol[ix][iy]+1.0e-30);
#endif
    }
  }
}
