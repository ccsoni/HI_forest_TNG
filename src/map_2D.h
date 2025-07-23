#pragma once
#define NMESH_MAP (256)

struct map_2D{
  float cdm_dens[NMESH_MAP][NMESH_MAP];  
  float gas_dens[NMESH_MAP][NMESH_MAP];
  float tmpr[NMESH_MAP][NMESH_MAP];
  float fHI[NMESH_MAP][NMESH_MAP];
  float vol[NMESH_MAP][NMESH_MAP];  // for computing volume weighted values
  float xmin, xmax;
  float ymin, ymax;
  float slice_pos;
};

void zero_out_map(struct map_2D *);
void accumulate_cdm_to_map(struct map_2D *, struct tng_cdm *, struct tng_header *);
void accumulate_gas_to_map(struct map_2D *, struct tng_gas *, struct tng_header *);
void normalize_map(struct map_2D *, struct tng_header *);
void output_map(char *, struct map_2D *, struct tng_header *);
