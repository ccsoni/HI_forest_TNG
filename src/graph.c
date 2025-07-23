#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include "tng_data.h"
#include "graph.h"

struct vector relative_vector(struct point *p0, struct point *p1)
{
  struct vector v;
  v.x = p1->xpos - p0->xpos;
  v.y = p1->ypos - p0->ypos;
  v.z = p1->zpos - p0->zpos;

  return v;
}

struct line connecting_line(struct point *p0, struct point *p1)
{
  struct line connect_line;
  connect_line.vect = relative_vector(p0, p1);
  
  normalize(&connect_line.vect);
  connect_line.pnt = *p0;

  return connect_line;
}

REAL dot_product(struct vector *v1, struct vector *v2)
{
  REAL dot;

  dot = v1->x*v2->x + v1->y*v2->y + v1->z*v2->z;

  return dot;
}

struct vector cross_product(struct vector *v1, struct vector *v2)
{
  struct vector vect;

  vect.x = v1->y*v2->z - v1->z*v2->y;
  vect.y = v1->z*v2->x - v1->x*v2->z;
  vect.z = v1->x*v2->y - v1->y*v2->x;

  return vect;
}

void normalize(struct vector *v)
{
  REAL norm;

  norm = sqrt(NORM2(v->x, v->y, v->z))+FLT_MIN;
  v->x /= norm;
  v->y /= norm;
  v->z /= norm;

}

void compute_bisector_surface(struct surface *surf, struct point *p1, struct point *p2)
{
  REAL dx, dy, dz, dr2;

  surf->cent.xpos = 0.5*(p1->xpos+p2->xpos);
  surf->cent.ypos = 0.5*(p1->ypos+p2->ypos);
  surf->cent.zpos = 0.5*(p1->zpos+p2->zpos);

  dx = p1->xpos - p2->xpos;
  dy = p1->ypos - p2->ypos;
  dz = p1->zpos - p2->zpos;
  dr2 = NORM2(dx, dy, dz)+FLT_MIN;

  surf->norm_vect.x = dx/sqrt(dr2);
  surf->norm_vect.y = dy/sqrt(dr2);
  surf->norm_vect.z = dz/sqrt(dr2);
}

void compute_three_point_surface(struct surface *s, 
                                 struct point *p0, 
                                 struct point *p1, 
                                 struct point *p2)
{
  struct vector v01, v02;

  v01.x = p1->xpos - p0->xpos;
  v01.y = p1->ypos - p0->ypos;
  v01.z = p1->zpos - p0->zpos;

  v02.x = p2->xpos - p0->xpos;
  v02.y = p2->ypos - p0->ypos;
  v02.z = p2->zpos - p0->zpos;

  s->norm_vect = cross_product(&v01, &v02);

  normalize(&s->norm_vect);
  s->cent = *p1;
}

void compute_crosspoint_three_surface(struct point *pnt, 
                                      struct surface *s0, 
                                      struct surface *s1, 
                                      struct surface *s2)
{
  REAL det; /* determinant of the three normal vectors */
  REAL sub_det0, sub_det1, sub_det2;

  sub_det0 = s1->norm_vect.y*s2->norm_vect.z - s1->norm_vect.z*s2->norm_vect.y;
  sub_det1 = s1->norm_vect.z*s2->norm_vect.x - s1->norm_vect.x*s2->norm_vect.z;
  sub_det2 = s1->norm_vect.x*s2->norm_vect.y - s1->norm_vect.y*s2->norm_vect.x;

  det = s0->norm_vect.x*sub_det0 
    +   s0->norm_vect.y*sub_det1 
    +   s0->norm_vect.z*sub_det2;

  REAL pn0, pn1, pn2;

  pn0 = s0->norm_vect.x*s0->cent.xpos 
    +   s0->norm_vect.y*s0->cent.ypos
    +   s0->norm_vect.z*s0->cent.zpos;

  pn1 = s1->norm_vect.x*s1->cent.xpos 
    +   s1->norm_vect.y*s1->cent.ypos
    +   s1->norm_vect.z*s1->cent.zpos;

  pn2 = s2->norm_vect.x*s2->cent.xpos 
    +   s2->norm_vect.y*s2->cent.ypos
    +   s2->norm_vect.z*s2->cent.zpos;

  struct vector n12, n20, n01;
  n12 = cross_product(&s1->norm_vect, &s2->norm_vect);
  n20 = cross_product(&s2->norm_vect, &s0->norm_vect);
  n01 = cross_product(&s0->norm_vect, &s1->norm_vect);

  pnt->xpos = (pn0*n12.x + pn1*n20.x + pn2*n01.x)/det;
  pnt->ypos = (pn0*n12.y + pn1*n20.y + pn2*n01.y)/det;
  pnt->zpos = (pn0*n12.z + pn1*n20.z + pn2*n01.z)/det;

}

REAL compute_distance_point_surface(struct point *pnt,
                                    struct surface *s)
{
  struct line l;
  struct point p;
  REAL dist;

  l.pnt = *pnt;
  l.vect = s->norm_vect;

  dist = compute_crosspoint_line_surface(&p,&l,s);

  return fabs(dist);
}

REAL compute_distance_point_line(struct point *pnt, struct line *l)
{
  struct vector rel_v,proj_v;
  REAL proj_length;

  rel_v.x = l->pnt.xpos - pnt->xpos;
  rel_v.y = l->pnt.ypos - pnt->ypos;
  rel_v.z = l->pnt.zpos - pnt->zpos;

  proj_length = dot_product(&rel_v, &l->vect);

  proj_v.x = proj_length*l->vect.x;
  proj_v.y = proj_length*l->vect.y;
  proj_v.z = proj_length*l->vect.z;

  struct vector v;

  v.x = rel_v.x - proj_v.x;
  v.y = rel_v.y - proj_v.y;
  v.z = rel_v.z - proj_v.z;

  pnt->xpos += v.x;
  pnt->ypos += v.y;
  pnt->zpos += v.z;

  return sqrt(NORM2(v.x, v.y, v.z));
  
}

REAL compute_crosspoint_line_surface(struct point *p, 
                                     struct line *l, struct surface *s)
{
  struct vector v;

  REAL affine, denom;
  
  v.x = s->cent.xpos - l->pnt.xpos;
  v.y = s->cent.ypos - l->pnt.ypos;
  v.z = s->cent.zpos - l->pnt.zpos;

  denom = dot_product(&(s->norm_vect), &(l->vect));
  assert(fabs(denom) > 1.0e-33);

  affine = 
    dot_product(&(s->norm_vect), &v)/dot_product(&(s->norm_vect), &(l->vect));

  p->xpos = l->pnt.xpos + affine*l->vect.x;
  p->ypos = l->pnt.ypos + affine*l->vect.y;
  p->zpos = l->pnt.zpos + affine*l->vect.z;

  return affine;
}

