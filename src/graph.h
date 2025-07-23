#pragma once

typedef double REAL;
typedef unsigned int UINT;

struct point
{
  UINT indx;
  REAL xpos, ypos, zpos;
};

struct vector
{
  REAL x,y,z;
};

struct line
{
  struct point pnt;
  struct vector vect;
};

struct surface
{
  struct point cent;
  struct vector norm_vect;
};

/* prototypes */
struct vector relative_vector(struct point*, struct point*);
struct line connecting_line(struct point*, struct point*);
REAL dot_product(struct vector*, struct vector*);
struct vector cross_product(struct vector*, struct vector*);
void normalize(struct vector*);
void compute_bisector_surface(struct surface*, struct point*, struct point*);
void compute_crosspoint_three_surface(struct point*, struct surface*, struct surface*, struct surface*);
void compute_three_point_surface(struct surface*, struct point*, struct point*, struct point*);
REAL compute_distance_point_surface(struct point*, struct surface*);
REAL compute_distance_point_line(struct point*, struct line*);
REAL compute_crosspoint_line_surface(struct point*, struct line*, struct surface*);
