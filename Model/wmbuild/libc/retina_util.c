/*****************************************************************************/
/*                                                                           */
/*   retina_util.c                                                           */
/*   wyeth bair                                                              */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "my_util.h"
#include "misc_util.h"
#include "carray_util.h"
#include "iarray_util.h"
#include "farray_util.h"
#include "func_util.h"
#include "myrand_util.h"

#define RXCENTER  2000.0 // When these are 0.0, there is an artifact
#define RYCENTER  2000.0

struct cell_struct{ // For modeling coordinates of photoreceptors
  float x,y;
  struct cell_struct *prev,*next;
};

//
//  Will hold pointers to hex neighbors
//
struct h_cell_struct{ // For modeling coordinates of photoreceptors
  float x,y;
  int n;                           // Number of neighbors
  int ui;                          // Unit index (if we're keeping this)
  int cid;
  int edge_flag;                   // Used to extract unique edges
  struct h_cell_struct **neighbor; // Pointers to neighbors [n][]
};

/**************************************-**************************************/
/*                                                                           */
/*                          RETUTIL_GET_CONE_DIAM_UM                         */
/*                                                                           */
/*  Human data from Tyler (1985) JOSA 2:393--398, Fig 5, outer segment.      */
/*                                                                           */
/*****************************************************************************/
float retutil_get_cone_diam_um(ecc,spec)
     float ecc;
     char spec[];
{
  float radum;  // Cone radius in microns

  if (strcmp(spec,"Tyler85_OD")==0)
    radum = 1.4 * pow((ecc + 0.2),0.2);
  else
    radum = 1.4 * pow((ecc + 0.2),0.2);  // Use Tyler85

  return radum;
}
/**************************************-**************************************/
/*                                                                           */
/*                            RETUTIL_GET_CONE_DENSITY                       */
/*                                                                           */
/*  Compute cone density as a function of distance from fovea in mm.         */
/*  Cone density is given in cones/mm^2.                                     */
/*                                                                           */
/*  Given an average sphere diameter of 21.45 mm, there are 5.34 degr/mm.    */
/*  (See table 1 in Packer et al. 1990).                                     */
/*  Note also, retinal area 800mm^2, 55% of sphere covered by retina.        */
/*  Cone size varies from 2.3--10um from fovea to periphery.                 */
/*                                                                           */
/*  This was another set of equations for density:                           */
/*    x = -(0.45 + r/(r+0.1) * 0.45);                                        */
/*    dens = 28000.0 * pow(r,x);                                             */
/*    if (r-0.023 >= 0.0)                                                    */
/*    fprintf(fout,"%.3f %.3f\n",r-0.023,dens);                              */
/*                                                                           */
/*****************************************************************************/
float retutil_get_cone_density(c,r)
     float c;  // some constant, e.g., 0.95
     float r;  // distance from center (mm)
{
  float a,b,dens,d0,d10;

  d0 = 210000.0; // Density at center
  d10 = 3500.0;  // Density at 10mm from center
  b = d10*pow(10.0,c)/(d0-d10);
  a = b*d0;
  dens = a/(b+pow(r,c));

  return dens;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MINIMIZE_MASK_LENGTH                           */
/*                                                                           */
/*  Reduce the mask length from the right until the first point is found     */
/*  which exceeds 'tolerance' times the maximum of the absolute value of     */
/*  the mask.                                                                */
/*                                                                           */
/*****************************************************************************/
void minimize_mask_length(mask,nm,tolerance)
     float **mask;
     int *nm;
     float tolerance;
{
  int n;
  float eps,min,max,*t;

  n = *nm;

  get_min_max_farray(*mask,n,&min,&max);
  if (min < 0.0)
    eps = tolerance * mymaxf(-min,max);
  else
    eps = tolerance * max;

  while ((fabs((*mask)[n-1]) < eps) && (n > 1))
    n--;

  if (n == 1)
    exit_error("MINIMIZE_MASK_LENGTH","Mask shrunk to length 1");
  if (n == *nm)
    printf("  *** WARNING: last value in mask exceeds tolerance.\n");

  if (n < *nm){ /*** Shorten the mask. ***/
    t = copy_farray(*mask,n);
    myfree(*mask);
    *mask = t;
    *nm = n;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            RETUTIL_GET_CONE_FILTER                        */
/*                                                                           */
/*  Choose the kernel that will represent the cone impulse responses.        */
/*                                                                           */
/*****************************************************************************/
float *retutil_get_cone_filter(tau_r,tau_d,tau_p,phi,rn)
     float tau_r,tau_d,tau_p,phi;
     int *rn;
{
  int n;
  float *mask;

  /*** These are other possible functions. ***/
  /* mask = alpha_farray(0.0,0.2,1.0,nm1);*/
  /* mask = maxwell_farray(0.0,4.0,1.0,nm1);*/ /* Has area 1.0 */
  /* mask = func_lognormal_farray(0.0,10.0,1.0,nm1);*/
  /* mask = func_photoreceptor_farray(0.0,1.0,1.0,nm1);*/

  n = 500;

  mask = func_primate_photocurrent_farray(0.0,1.0,tau_r,tau_d,tau_p,phi,n);

  append_farray_plot("zzz.photo.mask.pl","MASK",mask,n,1);

  minimize_mask_length(&mask,&n,0.005); /* 0.5% precision */
  norm_area_farray(mask,n,1.0);

  *rn = n;
  return mask;
}
/**************************************-**************************************/
/*                                                                           */
/*                             RETUTIL_INIT_RETINA                           */
/*                                                                           */
/*****************************************************************************/
void retutil_init_retina(flag,n,dist,rcx,rcy,rcn,rb,rbn)
     int flag,n;
     float dist,**rcx,**rcy;
     int *rcn;
     struct cell_struct **rb;
     int *rbn;
{
  int cn,bn;
  float *cx,*cy;
  struct cell_struct *b,*p;

  // Set number of photoreceptors in retina.  Make storage for coords
  cx = (float *)myalloc(n*sizeof(float));
  cy = (float *)myalloc(n*sizeof(float));

  if (flag == 1){ // Start with triangle around the center
    cx[0] = RXCENTER;
    cy[0] = RYCENTER;
    cx[1] = cx[0] + dist;
    cy[1] = cy[0];
    cx[2] = cx[0] + dist/2.0;
    cy[2] = cy[0] + dist*sqrt(3.0)/2.0;
    cn = 3;
    // Set up border ring, b points to the original beginning
    b = (struct cell_struct *)myalloc(sizeof(struct cell_struct));
    b->x = cx[0]; b->y = cy[0];
    b->next = (struct cell_struct *)myalloc(sizeof(struct cell_struct));
    p = b->next;
    p->x = cx[1]; p->y = cy[1]; p->prev = b;
    p->next = (struct cell_struct *)myalloc(sizeof(struct cell_struct));
    p->next->prev = p;
    p = p->next;
    p->x = cx[2]; p->y = cy[2];
    p->next = b;
    b->prev = p;
    bn = 3;
  }else if (flag == 2){ /*** Make square at center. ***/
    cx[0] = RXCENTER; cy[0] = RYCENTER;
    cx[1] = cx[0] + dist; cy[1] = cy[1];
    cx[2] = cx[1]; cy[2] = cy[1] + dist;
    cx[3] = cx[0];  cy[3] = cy[2];
    cn = 4;
    /*** Set up border ring, b points to the original beginning. ***/
    b = (struct cell_struct *)myalloc(sizeof(struct cell_struct));
    b->x = cx[0]; b->y = cy[0];
    b->next = (struct cell_struct *)myalloc(sizeof(struct cell_struct));
    b->next->prev = b;
    p = b->next;
    p->x = cx[1]; p->y = cy[1];
    p->next = (struct cell_struct *)myalloc(sizeof(struct cell_struct));
    p->next->prev = p;
    p = p->next;
    p->x = cx[2]; p->y = cy[2];
    p->next = (struct cell_struct *)myalloc(sizeof(struct cell_struct));
    p->next->prev = p;
    p = p->next;
    p->x = cx[3]; p->y = cy[3];
    p->next = b;
    b->prev = p;
    bn = 4;
  }else{
    cn = bn = 0;
    b = NULL;
    exit_error("INIT_RETINA","Unknown flag value");
  }

  *rcx = cx; *rcy = cy; *rcn = cn; *rb = b; *rbn = bn;
}
/**************************************-**************************************/
/*                                                                           */
/*                               RETUTIL_GET_ANGLE                           */
/*                                                                           */
/*  a,b,c - points along the border.                                         */
/*                                                                           */
/*****************************************************************************/
float retutil_get_angle(a1,a2,b1,b2,c1,c2)
     float a1,a2,b1,b2,c1,c2;
{
  float x,y,u1,u2,v1,v2,theta;

  u1 = b1-a1; u2 = b2-a2;
  v1 = c1-b1; v2 = c2-b2;
  x =  u1*v1 + u2*v2;
  y = -u2*v1 + u1*v2;

  theta = 180.0/M_PI * (atan2(y,x)+ M_PI);
  return theta;
}
/**************************************-**************************************/
/*                                                                           */
/*                              PROCESS_NEXT_EDGE                            */
/*                                                                           */
/*                                                                           */
/*                                           ^                               */
/*        (c1,c2) *                         /.\                              */
/*                 \                      __ .                               */
/*                  ^                     vo .                               */
/*                   \  _                    .                               */
/*                    \ u                    .                               */
/*                     \                     .                               */
/*                      \                    .                               */
/*                       \          _        .                               */
/*                        \         v        .                               */
/*                 (b1,b2) *-<---------------* (a1,a2)                       */
/*                                                                           */
/*  a,b,c - points along the border.                                         */
/*  u,v,vo - vectors, (vo is orthogonal to v)                                */
/*  r,ro - vectors, r is projection of u on v, ro is u on vo.                */
/*                                                                           */
/*****************************************************************************/
void process_next_edge(rb,bn,cx,cy,cn,dist,dnoise,pseed,k)
     struct cell_struct **rb;
     int *bn;
     float *cx,*cy;
     int *cn;
     float dist,dnoise;
     int *pseed;
     int k;
{
  struct cell_struct *p,*t;
  int flag;
  /*
  float a1,a2,b1,b2,c1,c2;
  float u1,u2,v1,v2,vo1,vo2,w1,w2,wo1,wo2,r1,r2,ro1,ro2;
  float m,x,y,mr2,mu2,mv2,mvo2,nx,ny,xnew,ynew,tdist2,distsq;
  float dab2,dnew;
  */

  //
  //  WYETH - I made these double because floats were different on
  //  MacBook Pro and linux boxes
  //
  double a1,a2,b1,b2,c1,c2;
  double u1,u2,v1,v2,vo1,vo2,w1,w2,wo1,wo2,r1,r2,ro1,ro2;
  double m,x,y,mr2,mu2,mv2,mvo2,nx,ny,xnew,ynew,tdist2,distsq;
  double dab2,dnew;


  p = *rb;
  distsq = dist*dist;
  a1 = p->x; a2 = p->y;
  b1 = p->next->x; b2 = p->next->y;
  c1 = p->next->next->x; c2 = p->next->next->y;
  v1 = b1-a1; v2 = b2-a2;
  u1 = c1-b1; u2 = c2-b2;
  w1 = c1-a1; w2 = c2-a2;
  vo1 = v2; vo2 = -v1;
  wo1 = w2; wo2 = -w1;

  /*** projection_2d(u1,u2,v1,v2,&r1,&r2);  ***  Get proj. u on v. ***/
  mv2 = (v1*v1 + v2*v2); /* Squared magnitude of v */
  m = (u1*v1 + u2*v2)/mv2;
  r1 = m*v1;
  r2 = m*v2;

  /*** projection_2d(u1,u2,vo1,vo2,&ro1,&ro2); ** Get proj. u on orth. v. ***/
  mvo2 = (vo1*vo1 + vo2*vo2); /* Squared magnitude of v */
  m = (u1*vo1 + u2*vo2)/mvo2;
  ro1 = m*vo1;
  ro2 = m*vo2;

  mr2 = r1*r1 + r2*r2; /* Mag. of r, squared */
  mu2 = u1*u1 + u2*u2; /* Mag. of u, squared */

  /* Compute distance between a and b, squared */
  dab2 = mv2/4.0;

  if (distsq <= dab2)
    dnew = 0.0;
  else
    dnew = sqrt(distsq - dab2);

  m = sqrt(mvo2);
  vo1 /= m;
  vo2 /= m;
  xnew =  dnew*vo1 + (a1 + b1)/2.0;
  ynew =  dnew*vo2 + (a2 + b2)/2.0;

  tdist2 = (xnew-c1)*(xnew-c1) + (ynew-c2)*(ynew-c2);

  /*** Case 0-close, 1-parallelogram add, 2-simple add ***/
  if ((ro1*vo1 < 0.0)||(ro2*vo2 < 0.0)||(tdist2 >= distsq)) /*** > 180 ***/
    flag = 2;
  else if (mr2 < 0.25*mu2) /*** Within 60--120 range ***/
    flag = 1;
  else if ((r1*v1 > 0.0)||(r2*v2 > 0.0)) /*** Within 120--90 range ***/
    flag = 2;
  else /*** Within 0--60 range ***/
    flag = 0;

  if (tdist2 < distsq)
    if (100.0 > retutil_get_angle(a1,a2,c1,c2,p->next->next->next->x,
				  p->next->next->next->y)){
      flag = 0;
      /*printf(" ***THETA_AHEAD < 100.0,  %.2f\n",sqrt(tdist2));*/
    }

  if (flag==1){
    if (dnoise > 0.0){
      nx = (myrand_util_ran2(pseed)-0.5)*dnoise;
      ny = (myrand_util_ran2(pseed)-0.5)*dnoise;
    }else
      nx = ny = 0.0;
    /* Compute distance between a and c, squared */
    dab2 = w1*w1 + w2*w2;
    if (distsq <= dab2/4.0)
      dnew = 0.0;
    else
      dnew = sqrt(distsq - dab2/4.0);
    m = sqrt(dab2);
    wo1 /= m;
    wo2 /= m;
    xnew =  dnew*wo1 + (a1 + c1)/2.0;
    ynew =  dnew*wo2 + (a2 + c2)/2.0;
    x = cx[*cn] = xnew + nx;
    y = cy[*cn] = ynew + ny;
    *cn += 1;
    /*** Replace border item "b" with new x,y values. ***/
    t = p->next;
    t->x = x; t->y = y;
    *rb = t; /*** Move current position to new cell. ***/
  }else if (flag==2){
    if (dnoise > 0.0){
      nx = (myrand_util_ran2(pseed)-0.5)*dnoise;
      ny = (myrand_util_ran2(pseed)-0.5)*dnoise;
    }else
      nx = ny = 0.0;
    x = cx[*cn] = nx + xnew;
    y = cy[*cn] = ny + ynew;
    *cn += 1;

    /*printf("a: %.2f,%.2f  b: %.2f,%.2f  c: %.2f,%.2f\n",a1,a2,b1,b2,c1,c2);
      printf("new: %.2f,%.2f\n",x,y);*/

    /*** Insert into border array. ***/
    t = (struct cell_struct *)myalloc(sizeof(struct cell_struct));
    t->x = x; t->y = y;
    t->next = p->next;
    t->prev = p;
    p->next->prev = t;
    p->next = t;
    *bn += 1;
    *rb = t; /*** Move current position to new cell. ***/
  }else{
    t = p->next;
    p->next = t->next;
    p->next->prev = p;
    myfree(t);
    *bn -= 1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            RETUTIL_SET_NEIGHBORS                          */
/*                                                                           */
/*****************************************************************************/
void retutil_set_neighbors(clist,n,xi,yi,nw,nh,ci,wideflag)
     struct h_cell_struct *clist;
     int n;
     int xi;  // Horizontal coord
     int yi;  // Vertical coord
     int nw;  // Number in wide row
     int nh;  // Number of rows
     int ci;  // Cell index
     int wideflag;  // Is current index in wide row
{
  int k;
  struct h_cell_struct *hc;
  int f0,f1,f2,f3,f4,f5;

  //
  //        2 3
  //      0  *  1
  //        5 4
  //

  // All flags to 1
  f0 = f1 = f2 = f3 = f4 = f5 = 1;

  // Take care of the ends of each row
  if (wideflag == 1){
    if (xi == 0)
      f0 = f2 = f5 = 0;
    else if (xi == (nw-1))
      f1 = f3 = f4 = 0;
  }else{
    if (xi == 0)
      f0 = 0;
    else if (xi == (nw-2))
      f1 = 0;
  }

  // Take care of top and bottom
  if (yi == 0)
    f4 = f5 = 0;
  else if (yi == (nh-1))
    f2 = f3 = 0;

  hc = &(clist[ci]);

  hc->n = f0 + f1 + f2 + f3 + f4 + f5;

  hc->neighbor = (struct h_cell_struct **)
    myalloc(hc->n*sizeof(struct h_cell_struct *));

  k = 0;
  if (f0 == 1){
    hc->neighbor[k] = &(clist[ci-1]);
    k += 1;
  }
  if (f1 == 1){
    hc->neighbor[k] = &(clist[ci+1]);
    k += 1;
  }
  if (f2 == 1){
    hc->neighbor[k] = &(clist[ci+(nw-1)]);
    k += 1;
  }
  if (f3 == 1){
    hc->neighbor[k] = &(clist[ci+nw]);
    k += 1;
  }
  if (f4 == 1){
    hc->neighbor[k] = &(clist[ci-(nw-1)]);
    k += 1;
  }
  if (f5 == 1){
    hc->neighbor[k] = &(clist[ci-nw]);
    k += 1;
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                        RETUTIL_GET_CONE_COORDS_REGJIT                     */
/*                                                                           */
/*  Start with a regular hex array, and then jitter the coordinates.         */
/*                                                                           */
/*  Coordinates returned are in pixels centered around (0,0).                */
/*                                                                           */
/*****************************************************************************/
void retutil_get_cone_coords_regjit(n,dens,uperpix,dnoise,seed,rn,rclist,
				    rnw,rnh)
     int n;              // Number of cones to attempt
     float dens;         // cones / mm^2
     float uperpix;      // microns per pixel
     float dnoise;
     int seed;
     int *rn;            // Number of cones, will typically be > n
     struct h_cell_struct **rclist; // Structure of linked cones
     int *rnw;           // Width of wide row, in cones
     int *rnh;           // Height - Total number of rows
{
  int i,j,k;
  int nw,nh,ncone,nhalf,ymid,wideflag;
  float *cx,*cy,cdist,edist,side_len_mm,dy,yval,rx,ry,d,noise_amp;
  struct h_cell_struct *clist;

  side_len_mm = sqrt(n/dens);     // Length of side of squre region (mm)
  edist = sqrt(1.0/(sqrt(3.0)/2.0*dens));  // Dist between cones (mm)
  cdist = edist * 1000.0/uperpix; // Estimated distance in pixels
  dy = sqrt(3.0)/2.0 * cdist;     // Distance between rows

  nhalf = (int)(side_len_mm / edist / 2);  // Number of cones to right of mid
  nw    = nhalf * 2 + 1;                   // Number of cones across middle


  //
  //  Determine 'ncone' (total cones) and 'nh' (number of rows)
  //
  ncone = nw;   // Add in the first (central) row, which is a wide row
  nh = 1;       // Count how many rows we have added
  wideflag = 1; // Indicate that the most recently added row is a wide one
  while(ncone < n){
    if (wideflag == 1){  // WYETH BUG: was always narrow ((nh-1)/2)%1 == 0){
      ncone += 2 * (nw-1);  // Two new narrow rows
    }else{
      ncone += 2 * nw;      // Two new wide rows
    }
    wideflag = 1 - wideflag;
    nh += 2;
    //printf("%d  %d  %d\n",wideflag,nh,ncone);
  }

  //
  //  Now,
  //    nh     - Height of hex-grid in rows
  //    ncones - Total cones in grid
  //    nw     - Number of cones in a wide row
  //  BUT - THESE WILL NOT ALL BE KEPT BY THE CALLER
  //
  //printf("   n  = %d\n",n);
  //printf("   nh = %d\n",nh);
  //printf("   nw = %d\n",nw);
  //printf("   ncones = %d\n",ncone);
  //printf("   wideflag = %d\n",wideflag);
  //exit(0);

  //
  //  Build regular hex grid as a linked structure, to indicate neighbors
  //
  clist = (struct h_cell_struct *)myalloc(ncone*sizeof(struct h_cell_struct));
  ymid = (nh-1)/2;
  k = 0;


  for(i=0;i<nh;i++){  // For each row; 0 is BOTTOM row

    //printf("i = %d  (nh %d,  ncone %d)  k= %d\n",i,nh,ncone,k);

    yval = (i - ymid)*dy;

    if (wideflag == 1){
      for(j=0;j<nw;j++){
	clist[k].x = (j - nhalf)*cdist;
	clist[k].y = yval;
	retutil_set_neighbors(clist,ncone,i,j,nw,nh,k,wideflag);
	k += 1;
      }
    }else{
      for(j=0;j<(nw-1);j++){
	clist[k].x = ((float)(j - nhalf) + 0.5)*cdist;
	clist[k].y = yval;
	retutil_set_neighbors(clist,ncone,i,j,nw,nh,k,wideflag);
	k += 1;
      }
    }
    wideflag = 1 - wideflag;
  }


  //
  //  Add noise to the cone positions
  //
  if ((dnoise > 0.0) && (dnoise <= 1.0)){

    //  Maximum noise distance is 1/2 way to near cone
    noise_amp = dnoise * cdist/2.0;

    if (seed > 0)
      seed = -seed;

    for(i=0;i<ncone;i++){

      rx = 2.0*myrand_util_ran2(&seed) - 1.0;  // -1 to 1
      ry = 2.0*myrand_util_ran2(&seed) - 1.0;  // -1 to 1
      d = rx*rx + ry*ry;
      while(d >= 1.0){
	rx = 2.0*myrand_util_ran2(&seed) - 1.0;  // -1 to 1
	ry = 2.0*myrand_util_ran2(&seed) - 1.0;  // -1 to 1
	d = rx*rx + ry*ry;
      }

      clist[i].x += rx * noise_amp;
      clist[i].y += ry * noise_amp;
    }
  }else if (dnoise != 0.0){
    exit_error("RETUTIL_GET_CONE_COORDS_REGJIT","Noise value out of range");
  }

  *rn = ncone;
  *rclist = clist;
  *rnw = nw;
  *rnh = nh;
}
/**************************************-**************************************/
/*                                                                           */
/*                           RETUTIL_GET_CONE_COORDS                         */
/*                                                                           */
/*  Coordinates returned are in pixels centered around (0,0).                */
/*                                                                           */
/*****************************************************************************/
void retutil_get_cone_coords(n,dens,uperpix,dnoise,seed,rcx,rcy)
     int n;
     float dens;         // cones / mm^2
     float uperpix;      // microns per pixel
     float dnoise;
     int seed;
     float **rcx,**rcy;
{
  int i;
  int flag,cn,bn,diam;
  float *cx,*cy,cdist,edist;
  struct cell_struct *b;

  if (seed > 0)
    seed = -seed;

  edist = sqrt(1.0/(sqrt(3.0)/2.0*dens));
  cdist = edist * 1000.0/uperpix; // Estimated distance in pixels

  flag = 1; // 1-triangle, 2-square
  retutil_init_retina(flag,n,cdist,&cx,&cy,&cn,&b,&bn);

  diam = (int)(cdist-1.0);
  i = 0;
  while(cn < n){
    process_next_edge(&b,&bn,cx,cy,&cn,cdist,dnoise,&seed,i);
    // cn may or maynot be incremented by 1

    i += 1;
  }
  
  for(i=0;i<n;i++){
    cx[i] -= RXCENTER;
    cy[i] -= RYCENTER;
  }

  // WYETH - SHOULD FREE STORAGE ring? pointed to by 'b'
  
  *rcx = cx; *rcy = cy;
}
/**************************************-**************************************/
/*                                                                           */
/*                           RETUTIL_MOSAIC_REG_S                            */
/*                                                                           */
/*  Make a regular array of S-cones.                                         */
/*                                                                           */
/*****************************************************************************/
void retutil_mosaic_reg_s(betw,clist,n,nw,nh)
     int betw;         // Number of other cones between S-cones
     struct h_cell_struct *clist; // [n] list of cone structures
     int n;            // Number of cones
     int nw;           // Width of wide row, in cones
     int nh;           // Height - Total number of rows
{
  int i,j,k;
  int wide_row,rowi,tw,shift,shv,shift_val;

  // Note, the middle row is wide
  if ((nh - 1)%4 == 0){
    wide_row = 1;  // First row is wide
    tw = nw;
  }else{
    wide_row = 0;  // First row is wide
    tw = nw - 1;
  }

  shift_val = 1 + betw/2;

  //printf("    nw = %d\n",nw);
  //printf("    nh = %d\n",nh);
  //printf("    wide_flag = %d\n",wide_row);
  //printf("    shift_val = %d\n",shift_val);

  k = 0;  // Cone index w/i current row
  rowi = 0;
  shift = 0;  // Don't shift this S-row
  
  for(i=0;i<n;i++){
    if (rowi % (betw+1) == 0){
      // Put S cones in this row
      if ((k+shv) % (betw+1) == 0){
	if (clist[i].cid != -1){  // If this one is being kept
	  clist[i].cid = 0;       // Make an S-cone here
	}
      }
    }
    k += 1;
    if (k == tw){  // New row
      k = 0;
      wide_row = 1 - wide_row;  // Alternate between wide and narrow
      if (wide_row == 1)  // Set 'tn' to be the width of the current row
	tw = nw;
      else
	tw = nw - 1;

      if (rowi % (betw+1) == 0){  // This was an S-row
	shift = 1 - shift;
	if (shift == 1)
	  shv = shift_val;
	else
	  shv = 0;
      }
	
      rowi += 1;
    }
  }
  //exit(0);
}
/**************************************-**************************************/
/*                                                                           */
/*                           RETUTIL_GET_MOSAIC_COLOR                        */
/*                                                                           */
/*****************************************************************************/
void retutil_get_mosaic_color(xn,yn,n,dens,mupix,dnoise,seedx,seedc,prs,lmr,
			      mylogf,
			      rcx,rcy,rcid,rcn)
     int xn,yn;          // Underlying cartesian grid
     int n;              // Number of cones to attempt
     float dens;
     float mupix;
     float dnoise;
     int seedx,seedc;    // Seeds for location and color
     float prs,lmr;      // Prob of S cone; L:M ratio
     char *mylogf;       // Log file
     float **rcx,**rcy;
     int **rcid,*rcn;
{     
  int i,k;
  int *cid,*tid,cn,nn;
  float *cx,*cy,*tx,*ty,prm,rv;
  char tstr[SLEN];

  // The 0th cone should have coords (0,0) after this call
  retutil_get_cone_coords(n,dens,mupix,dnoise,seedx,&cx,&cy);

  if (seedc > 0)
    seedc = -seedc;

  cid = get_zero_iarray(n);

  cn = 0;
  prm = (1.0 - prs) / (lmr + 1.0); // Prob of M-cone
  for(i=0;i<n;i++){
    // Center coords in [0..xn)[0..yn]

    //cx[i] += (float)(xn-1)/2.0 - 1.0;  // WHY DID I HAVE "- 1.0" ???

    // WYETH - CONE CENTERING CHANGED, Mar 16, 2009
    cx[i] += (float)(xn-1)/2.0;
    cy[i] += (float)(yn-1)/2.0;

    //printf("cx,y  %f %f\n",cx[i],cy[i]);
    
    // For each cone within the 'xn' by 'yn' frame, assign a color
    if ((cx[i] >= 0.0) && (cx[i] < (double)(xn-1)) &&
	(cy[i] >= 0.0) && (cy[i] < (double)(yn-1))){
      rv = myrand_util_ran2(&seedc);
      if (rv < prs)
	cid[i] = 0;
      else if (rv < (prs + prm))
	cid[i] = 1;
      else
	cid[i] = 2;
      cn += 1;  // Count how many cones we'll keep
    }else
      cid[i] = -1;
  }
  sprintf(tstr,"    %d of %d cones within grid\n",cn,n);
  mylog(mylogf,tstr);

  // Revise the lists to contain only the cones we're keeping
  tx = (float *)myalloc(cn*sizeof(float));
  ty = (float *)myalloc(cn*sizeof(float));
  tid = (int *)myalloc(cn*sizeof(int));
  k = 0;
  for(i=0;i<n;i++){
    if (cid[i] >= 0){
      tx[k] = cx[i];
      ty[k] = cy[i];
      tid[k] = cid[i];
      k += 1;
    }
  }
  myfree(cx); myfree(cy); myfree(cid);
  *rcn = cn;
  *rcx = tx;
  *rcy = ty;
  *rcid = tid;
}
/**************************************-**************************************/
/*                                                                           */
/*                          RETUTIL_GET_MOSAIC_REGJIT                        */
/*                                                                           */
/*****************************************************************************/
void retutil_get_mosaic_regjit(xn,yn,n,dens,mupix,dnoise,seedx,seedc,prs,lmr,
			       regs,mylogf,
			       rcx,rcy,rcid,rcn,
			       redge0,redge1,redged,redgen,rmow,rmoh)
     int xn,yn;          // Underlying cartesian grid
     int n;              // Number of cones to attempt
     float dens;
     float mupix;
     float dnoise;
     int seedx,seedc;    // Seeds for location and color
     float prs,lmr;      // Prob of S cone; L:M ratio
     int regs;           // Regular S-grid (over-rides 'prs')
     char *mylogf;       // Log file
     float **rcx,**rcy;
     int **rcid,*rcn;
     int **redge0,**redge1;  // [redgen] Unique edge list endpoint unit indices
     double **redged;         // [redgen] Unique edge list distances
     int *redgen;            // Number of unique edges
     int *rmow;          // Width of mosaic - a wide row, the central row
     int *rmoh;          // Height of mosiac - total rows
{     
  int i,j,k;
  int *tid,cn,nn,ne,*edge0,*edge1,edgen,mow,moh;
  float x,y,*tx,*ty,prm,rv,dx,dy;
  double *edgedist;
  char tstr[SLEN];
  struct h_cell_struct *clist,**nlist;

  // The 0th cone should have coords (0,0) after this call
  retutil_get_cone_coords_regjit(n,dens,mupix,dnoise,seedx,&nn,&clist,
				 rmow,rmoh);  // Mosaic width and height


  //
  //  (0) Count which units we keep:  only those w/i the xn,yn stimulus region
  //  (1) Set the CID for each unit
  //  (2) Center each cone
  //  (3) Set the 'edge_flag' = 1 for units to be kept
  //  (4) Set the final index for the units to keep
  //
  if (seedc > 0)
    seedc = -seedc;

  if (regs > 0)  // If we're making a regular S-mosaic, override 'prs' value
    prs = 0.0;

  cn = 0;  // Number of units to keep
  prm = (1.0 - prs) / (lmr + 1.0); // Prob of M-cone
  for(i=0;i<nn;i++){

    // Center coords in [0..xn)[0..yn]
    x = clist[i].x + (float)(xn-1)/2.0;
    y = clist[i].y + (float)(yn-1)/2.0;
    clist[i].x = x;
    clist[i].y = y;

    // For each cone within the 'xn' by 'yn' frame, assign a color
    if ((x >= 0.0) && (x < (double)(xn-1)) &&
	(y >= 0.0) && (y < (double)(yn-1))){
      rv = myrand_util_ran2(&seedc);
      if (rv < prs)
	clist[i].cid = 0;
      else if (rv < (prs + prm))
	clist[i].cid = 1;
      else
	clist[i].cid = 2;

      clist[i].edge_flag = 1;
      clist[i].ui = cn;

      cn += 1;  // Count how many cones we'll keep
    }else{
      clist[i].cid = -1;  // Not keeping
      clist[i].edge_flag = 0;
      clist[i].ui = -1;   // Not keeping
    }
  }
  sprintf(tstr,"    %d of %d cones within grid\n",cn,nn);
  mylog(mylogf,tstr);

  // Make a regular set of S-cones
  if (regs > 0)
    retutil_mosaic_reg_s(regs,clist,nn,*rmow,*rmoh);


  //
  //  Get list of unique edges, with distances
  //
  edge0 = (int *)myalloc(cn*6*sizeof(int));
  edge1 = (int *)myalloc(cn*6*sizeof(int));
  edgedist = (double *)myalloc(cn*6*sizeof(double));
  k = 0;  // Edge counter
  for(i=0;i<nn;i++){  // For each Unit

    if (clist[i].cid != -1){   // If this unit is being kept

      x = clist[i].x;
      y = clist[i].y;

      ne = clist[i].n;             // Number of neigbhors for this unit
      nlist = clist[i].neighbor;   // Neighbor list

      for(j=0;j<ne;j++){
	if (nlist[j]->edge_flag == 1){

	  dx = x - nlist[j]->x;
	  dy = y - nlist[j]->y;
	  edgedist[k] = sqrt(dx*dx + dy*dy);

	  edge0[k] = clist[i].ui;
	  edge1[k] = nlist[j]->ui;
	  k += 1;
	  //  Add to edge list
	}
      }

      clist[i].edge_flag = 0;  // Edges have been counted
    }
  }
  edgen = k;

  sprintf(tstr,"    %d unique edges\n",edgen);
  mylog(mylogf,tstr);


  // Revise the lists to contain only the cones we're keeping
  tx = (float *)myalloc(cn*sizeof(float));
  ty = (float *)myalloc(cn*sizeof(float));
  tid = (int *)myalloc(cn*sizeof(int));
  k = 0;
  for(i=0;i<nn;i++){
    if (clist[i].cid >= 0){
      tx[k]  = clist[i].x;
      ty[k]  = clist[i].y;
      tid[k] = clist[i].cid;
      if (k != clist[i].ui)
	exit_error("(retina_util)","THIS SHOULD NOT HAPPEN - ui bad");
      k += 1;
    }
  }

  //
  //  WYETH - FREE THE 'clist' ???
  //

  *rcn = cn;
  *rcx = tx;
  *rcy = ty;
  *rcid = tid;

  *redge0 = edge0;   // *** WYETH - could copy to shorter arrays, [edgen]
  *redge1 = edge1;
  *redged = edgedist;
  *redgen = edgen;

  // rmoh and rmow are set above
}
/**************************************-**************************************/
/*                                                                           */
/*                            RETUTIL_IRR2SQ_MASK                            */
/*                                                                           */
/*  Return a set of coordinates and weights on a square grid that            */
/*  surround the irregularly placed point (x,y) within that grid.  The       */
/*  weights are drawn from a Gaussian with SD 'sig'.  Points are only        */
/*  considered within 'nsd' SDs of the central point.                        */
/*                                                                           */
/*****************************************************************************/
void retutil_irr2sq_mask(x,y,xn,yn,sig,nsd,rn,rx,ry,rw)
     float  x,y;  // coordinate of irregular location
     int  xn,yn;  // Size of square grid to create
     float  sig;  // SD of Gaussian
     float  nsd;  // Number of SDs to consider weights
     int    *rn;  // Number of coordinates returned
     int   **rx;  // [rn] x-coordinates
     int   **ry;  // [rn] y-coordinates
     float **rw;  // [rn] weight
{
  int i,j,k;
  int n,xi,yi,x0,x1,y0,y1,rad,*cx,*cy;
  float *w,fc,dsq;

  rad = my_rint(nsd * sig);
  if (rad < 1)
    rad = 1;

  xi = my_rint(x);
  yi = my_rint(y);

  x0 = xi - rad;
  x1 = xi + rad;
  y0 = yi - rad;
  y1 = yi + rad;

  if (x0 < 0)
    x0 = 0;
  if (x1 > (xn-1))
    x1 = xn-1;

  if (y0 < 0)
    y0 = 0;
  if (y1 > (yn-1))
    y1 = yn-1;

  n = (x1 - x0 + 1) * (y1 - y0 + 1);
  if (n <= 0)
    exit_error("RETUTIL_IRR2SQ_MASK","Invalid n");

  cx = (int *)myalloc(n*sizeof(int));
  cy = (int *)myalloc(n*sizeof(int));
  w  = (float *)myalloc(n*sizeof(float));

  fc = -0.5/(sig*sig);
  
  k = 0;
  for(i=x0;i<=x1;i++){
    for(j=y0;j<=y1;j++){

      // Distance squared from (x,y) to (i,j)
      dsq = ((float)i-x)*((float)i-x) + ((float)j-y)*((float)j-y);

      cx[k] = i;
      cy[k] = j;
      w[k] = exp(fc * dsq);  // Gaussian function of distance
      k += 1;
    }
  }

  if (k != n)
    exit_error("RETUTIL_IRR2SQ_MASK","Should not happen");

  *rx = cx;
  *ry = cy;
  *rw = w;
  *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                             RETUTIL_HEX_TO_3D                             */
/*                                                                           */
/*****************************************************************************/
float ***retutil_hex_to_3d(hv,cn,tn,cx,cy,rxn,ryn)
     float **hv;  // [cn][tn] Values on hex grid
     int cn;      // Number of cones
     int tn;      // Duration
     float *cx;   // [cn] Hex x-coords on stim grid
     float *cy;   // [cn] Hex y-coords on stim grid
     int *rxn;    // Width of returned array
     int *ryn;    // Height of returned array
{
  int i;
  int ti,nw,nh,nwc,nhc,xi,yi;
  float ***d,xmin,xmax,ymin,ymax,val;
  double dx,dy;

  printf("  RETUTIL_HEX_TO_3D");

  dx = fabs(cx[1] - cx[0]);  // Cone 1 is directly to right
  dy = fabs(cy[2] - cy[0]);  // Cone 2 is in row above

  get_min_max_farray(cx,cn,&xmin,&xmax);
  get_min_max_farray(cy,cn,&ymin,&ymax);


  /*printf("    dx = %f   dy = %f\n",(float)dx,(float)dy);
    printf("    Cone.0  x,y  %f %f\n",cx[0],cy[0]);
    printf("    Cone  x min,max  %f  %f\n",xmin,xmax);
    printf("    Cone  y min,max  %f  %f\n",ymin,ymax);
  */

    
  //  Get width, height of cartesian grid
  nw = my_rint(1.0 + (xmax - xmin) / dx);
  nh = my_rint(1.0 + (ymax - ymin) / dy);

  nwc = 2*nw;
  nhc = 2*nh;

  d = get_zero_3d_farray(nwc,nhc,tn);
    
  //printf("  Cartesian grid:  %d %d\n",nw,nh);

  for(ti=0;ti<tn;ti++){
    for(i=0;i<cn;i++){

      val = hv[i][ti];

      xi = my_rint(2.0 * (cx[i] - xmin) / dx);
      yi = my_rint(2.0 * (cy[i] - ymin) / dy);

      // Store value in four locations.
      d[xi  ][yi  ][ti] = val;
      d[xi  ][yi+1][ti] = val;
      d[xi+1][yi  ][ti] = val;
      d[xi+1][yi+1][ti] = val;
    }
  }

  *rxn = nwc;
  *ryn = nhc;
  
  return d;
}
/**************************************-**************************************/
/*                                                                           */
/*                            RETUTIL_HEX_TO_3D_NEW                          */
/*                                                                           */
/*****************************************************************************/
float ***retutil_hex_to_3d_new(hv,cn,tn,cx,cy,xn,yn)
     float **hv;  // [cn][tn] Values on hex grid
     int cn;      // Number of cones
     int tn;      // Duration
     float *cx;   // [cn] Hex x-coords on stim grid
     float *cy;   // [cn] Hex y-coords on stim grid
     int xn;      // width of square grid
     int yn;      // height of square grid
{
  int i,j,k;
  int ti,xi,yi,*wn,**wx,**wy;
  float **ww,***d,***dw,xmin,xmax,ymin,ymax,*td,*tdw,w,*tval;

  printf("  RETUTIL_HEX_TO_3D_NEW\n");

  get_min_max_farray(cx,cn,&xmin,&xmax);
  get_min_max_farray(cy,cn,&ymin,&ymax);

  printf("    Cone  x min,max  %f  %f\n",xmin,xmax);
  printf("    Cone  y min,max  %f  %f\n",ymin,ymax);

  printf("    xn,yn =   %d , %d\n",xn,yn);

  //
  //  Compute weights for each cone position.
  //
  wn = (int *)myalloc(cn*sizeof(int));
  wx = (int **)myalloc(cn*sizeof(int *));
  wy = (int **)myalloc(cn*sizeof(int *));
  ww = (float **)myalloc(cn*sizeof(float *));
  for(i=0;i<cn;i++){
    retutil_irr2sq_mask(cx[i],cy[i],xn,yn,3.0,3.0,&(wn[i]),&(wx[i]),&(wy[i]),
			&(ww[i]));
  }

  printf("DONE getting weights  xn yn tn %d %d %d\n",xn,yn,tn);

  d  = get_zero_3d_farray(xn,yn,tn);
  dw = get_zero_3d_farray(xn,yn,tn);

  for(i=0;i<cn;i++){

    tval = hv[i];
    for(j=0;j<wn[i];j++){  // For each weighted coord for this cone
      xi = wx[i][j];
      yi = wy[i][j];
      w  = ww[i][j];

      td  =  d[xi][yi];  // Convenient pointers
      tdw = dw[xi][yi];

      for(ti=0;ti<tn;ti++){
	//td[ti] = w * tval[ti];   // Weight x ConeValue
	td[ti] += w * tval[ti];   // Weight x ConeValue
	tdw[ti] += w;
      }
    }
  }

  //write_3d_data("zzz.d.3d",d,xn,yn,tn,4,2,1);  // DEBUG
  //write_3d_data("zzz.w.3d",dw,xn,yn,tn,4,2,1);

  //
  //  Divide by total weight to produce weighted average
  //
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<tn;k++){
	if (dw[i][j][k] > 0.0)
	  d[i][j][k] /= dw[i][j][k];
      }
    }
  }

  //exit_error("  RETUTIL_HEX_TO_3D","Under development ... ");

  // *** WYETH FREE ALL THE WEIGHT STUFF ******
  // *** WYETH FREE ALL THE WEIGHT STUFF ******
  // *** WYETH FREE ALL THE WEIGHT STUFF ******
  
  return d;
}
/**************************************-**************************************/
/*                                                                           */
/*                            RETUTIL_VAL_DIST_PLOT                          */
/*                                                                           */
/*  Compute a list of distance vs. value pairs, given a cone ID.             */
/*                                                                           */
/*****************************************************************************/
void retutil_val_dist_plot(hv,cn,cx,cy,cid,t0,outfile)
     float **hv;  // [cn][tn] Values on hex grid
     int cn;      // Number of cones
     float *cx;   // [cn] Hex x-coords on stim grid
     float *cy;   // [cn] Hex y-coords on stim grid
     int cid;     // Cone ID
     int t0;      // Time
     char *outfile;  // Output file, or NULL
{
  int i;
  float *val,*dist,x,y,x1,y1;
  char pname[SLEN];

  printf("  RETUTIL_VAL_DIST_PLOT\n");

  val = (float *)myalloc(cn*sizeof(float));
  dist = (float *)myalloc(cn*sizeof(float));

  x = cx[cid];
  y = cy[cid];

  for(i=0;i<cn;i++){
    x1 = cx[i];
    y1 = cy[i];
    dist[i] = sqrt((x1-x)*(x1-x) + (y1-y)*(y1-y));
    val[i] = hv[i][t0];
  }

  sprintf(pname,"cone_%d_time_%d",cid,t0);
  append_farray_xy_plot(outfile,dist,val,cn,pname);
}
/**************************************-**************************************/
/*                                                                           */
/*                             RETUTIL_READ_MOSAIC                           */
/*                                                                           */
/*****************************************************************************/
void retutil_read_mosaic(infile,rcx,rcy,rcid,rcn)
     char infile[];
     float **rcx,**rcy;
     int **rcid,*rcn;
{
  FILE *fin,*fopen();
  int i;
  int *cid,cn,ns;
  float *cx,*cy,x,y;
  char c;

  if ((fin = fopen(infile,"r"))==NULL){
    printf("  filename = %s\n",infile);
    exit_error("RETUTIL_READ_MOSAIC","Cannot open file");
  }

  cn = 0;
  while(fscanf(fin,"%f %f %c",&x,&y,&c)!=EOF)
    cn += 1;
  fclose(fin);

  cx = (float *)myalloc(cn*sizeof(float));
  cy = (float *)myalloc(cn*sizeof(float));
  cid = (int *)myalloc(cn*sizeof(int));

  fin = fopen(infile,"r");
  for(i=0;i<cn;i++){
    ns = fscanf(fin,"%f %f %c",&x,&y,&c);
    cx[i] = x;
    cy[i] = y;
    if (c=='S')
      cid[i] = 0;
    else if (c=='M')
      cid[i] = 1;
    else if (c=='L')
      cid[i] = 2;
    else
      exit_error("RETUTIL_READ_MOSAIC","Bad cone type");
  }
  fclose(fin);

  *rcx = cx; *rcy = cy; *rcid = cid; *rcn = cn;
}
/**************************************-**************************************/
/*                                                                           */
/*                       RETUTIL_CUSTOM_MOSAIC_CID_FILE                      */
/*                                                                           */
/*****************************************************************************/
int retutil_custom_mosaic_cid_file(infile,cid,cn)
     char infile[];
     int *cid;        // [n] 0-S, 1-M, 2-L
     int cn;
{
  FILE *fin,*fopen();
  int i;
  int cnt;
  char c;

  if ((fin = fopen(infile,"r"))==NULL){
    printf("  filename = %s\n",infile);
    exit_error("RETUTIL_CUSTOM_MOSAIC_CID_FILE","Cannot open file");
  }
  
  cnt = 0;
  while(fscanf(fin,"%d %c",&i,&c)!=EOF){

    if ((i < 0) || (i >= cn)){
      printf("  Cone index value:  %d\n",i);
      exit_error("RETUTIL_CUSTOM_MOSAIC_CID_FILE","Cannot open file");
    }

    if (c=='S')
      cid[i] = 0;
    else if (c=='M')
      cid[i] = 1;
    else if (c=='L')
      cid[i] = 2;
    else{
      printf("  ConeID = %d, expecting 0,1,2\n",cid[i]);
      exit_error("RETUTIL_CUSTOM_MOSAIC_CID_FILE","Bad cone type");
    }
    
    cnt += 1;
  }
  fclose(fin);

  return cnt;  // The number of cones customised
}
