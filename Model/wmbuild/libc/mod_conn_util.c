/*****************************************************************************/
/*                                                                           */
/*  mod_conn_util.c                                                          */
/*  wyeth bair                                                               */
/*                                                                           */
/*  Establish connections between model units.                               */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "my_util.h"
#include "nr_util.h"
#include "carray_util.h"
#include "iarray_util.h"
#include "farray_util.h"
#include "func_util.h"
#include "myrand_util.h"
#include "data_util.h"
#include "paramfile_util.h"
#include "kernel_util.h"
#include "retina_util.h"
#include "pop_cell_util.h"
#include "mod_util.h"  // For mylog
#include "mod.h" // Data structures
#include "ifc.h" // Data structures
#include "paramfile.h"

// Global
int **radcx = (int **)NULL;
int **radcy;
int radmax = 6;

/**************************************-**************************************/
/*                                                                           */
/*                         MOD_CONN_ZSTRING_GET_UNIQUE_Z                     */
/*                                                                           */
/*****************************************************************************/
int mod_conn_zstring_get_unique_z(zn,zstring,isame)
     int zn;         // Depth of z-dimension
     char *zstring;  // "0","1",..."all", "even", "odd", "first", "last", ...
     int isame;      // Value to return if zstring is "same"
{
  int k;

  if (zn < 1)
    exit_error("MOD_CONN_ZSTRING_GET_UNIQUE_Z","zn < 1");

  if (is_int_string(zstring)){
    k = atoi(zstring);

    if ((k < 0) || (k >= zn))
      exit_error("MOD_CONN_GET_ZFLAG","Integer z-index out of bounds");

  }else if (strcmp(zstring,"first")==0){
    k = 0;
  }else if (strcmp(zstring,"last")==0){
    k = zn-1;
  }else if (strcmp(zstring,"same")==0){
    k = isame;

    if ((k < 0) || (k >= zn))
      exit_error("MOD_CONN_GET_ZFLAG","Integer z-index out of bounds");

  }else if (strcmp(zstring,"all")==0){
    k = -1;
  }else
    exit_error("MOD_CONN_ZSTRING_GET_UNIQUE_Z","Inappropriate z condition");


  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_CONN_ZSTRING_GET_ZFLAG                        */
/*                                                                           */
/*****************************************************************************/
int *mod_conn_zstring_get_zflag(zn,zstring)
     int zn;         // Depth of z-dimension
     char *zstring;  // "0","1",..."all", "even", "odd", "first", "last"
{
  int i;
  int *zf;

  //
  //
  //
  exit_error("MOD_CONN_ZSTRING_GET_ZFLAG","WYETH - NEVER USED");
  //
  //
  //

  if (zn < 1)
    exit_error("MOD_CONN_GET_ZFLAG","zn < 1");

  zf = get_zero_iarray(zn);

  if (is_int_string(zstring)){
    i = atoi(zstring);
    if ((i >= 0) && (i < zn))
      zf[i] = 1;
    else
      exit_error("MOD_CONN_ZSTRING_GET_ZFLAG","Integer z-index out of bounds");
  }else if (strcmp(zstring,"all")){
    for(i=0;i<zn;i++)
      zf[i] = 1;
  }else if (strcmp(zstring,"first")){
    zf[0] = 1;
  }else if (strcmp(zstring,"last")){
    zf[zn-1] = 1;
  }else if (strcmp(zstring,"even")){
    for(i=0;i<zn;i+=2)
      zf[i] = 1;
  }else if (strcmp(zstring,"odd")){
    for(i=1;i<zn;i+=2)
      zf[i] = 1;
  }else
    exit_error("MOD_CONN_ZSTRING_GET_ZFLAG","Unknown z condition");

  return zf;
}
/**************************************-**************************************/
/*                                                                           */
/*                      MOD_CONN_GET_DISTANCE_CELLS_LAYER                    */
/*                                                                           */
/*  Return distance in microns.                                              */
/*                                                                           */
/*****************************************************************************/
float mod_conn_get_distance_cells_layer(c1,c2)
     struct pop_cell *c1,*c2;
{
  float dx,dy,d;
  
  //printf("  Cell_1 cx,cy:  %f %f\n",c1->cx,c1->cy);
  //printf("  Cell_2 cx,cy:  %f %f\n",c2->cx,c2->cy);

  dx = c1->cx - c2->cx;  // These are in um
  dy = c1->cy - c2->cy;

  d = sqrt(dx*dx + dy*dy);

  return d;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_CONN_GET_CIRC_DIST                          */
/*                                                                           */
/*  Return the absolute value of the shortest distance between the two       */
/*  points on the circle of circumference 'cf'.                              */
/*                                                                           */
/*****************************************************************************/
float mod_conn_get_circ_dist(v1,v2,cf)
     float v1,v2;   // Return difference between these values
     float cf;      // Circumference
{
  float dist;

  //  WYETH - COULD USE:   (misc_util.c):   get_circular_diff(a,b,d)
  //  WYETH - COULD USE:   (misc_util.c):   get_circular_diff(a,b,d)
  //  WYETH - COULD USE:   (misc_util.c):   get_circular_diff(a,b,d)


  if (v1 >= v2)
    dist = v1 - v2;
  else
    dist = v2 - v1;

  if (dist > cf/2.0)
    dist = cf - dist;

  return dist;
}
/**************************************-**************************************/
/*                                                                           */
/*                       MOD_CONN_GET_SIGNED_CIRC_DIST                       */
/*                                                                           */
/*  Return the shortest signed-distance from 'ref' to 'v'.                   */
/*                                                                           */
/*****************************************************************************/
float mod_conn_get_signed_circ_dist(ref,v,cf)
     float ref,v;   // reference and value
     float cf;      // Circumference, e.g., 360.0
{
  float diff;

  // Remap 'v' to be within 'cf' in front of (larger than) the ref
  while(v < ref)
    v += cf;
  while(v >= (ref+cf))
    v -= cf;

  diff = v - ref;      // MUST be positive, thanks to remapping above
  if (diff > cf/2.0)
    diff = diff - cf;

  //
  //  returned value in (-cf/2 , cf/2], e.g.,  (-180.0, 180.0] for cf 360.0
  //

  return diff;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_CONN_GET_ORI_DIFF_CELLS                      */
/*                                                                           */
/*****************************************************************************/
float mod_conn_get_ori_diff_cells(c1,c2,oriflag)
     struct pop_cell *c1,*c2;
     int oriflag;  // 1-difference in ori, otherwise difference in dir
{
  float o1,o2,ddir,dori;

  // WYETH ATTRIB
  o1 = pop_cell_attrib_get_f(c1,"ori");
  o2 = pop_cell_attrib_get_f(c2,"ori");
  //o1 = c1->ori;
  //o2 = c2->ori;

  while(o1 >= 360.0)
    o1 -= 360.0;
  while(o1 < 0.0)
    o1 += 360.0;

  while(o2 >= 360.0)
    o2 -= 360.0;
  while(o2 < 0.0)
    o2 += 360.0;

  if (o1 >= o2)
    ddir = o1 - o2;
  else
    ddir = o2 - o1;

  if (ddir > 180.0)
    ddir = 360.0 - ddir;

  if (ddir > 90.0)
    dori = 180.0 - ddir;
  else
    dori = ddir;

  if (oriflag)
    return dori;
  else
    return ddir;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_CONN_GET_ORI_DIFF                          */
/*                                                                           */
/*****************************************************************************/
float mod_conn_get_ori_diff(o1,o2,oriflag)
     float o1,o2;  // Orientations (deg)p
     int oriflag;  // 1-difference in ori, otherwise difference in dir
{
  float ddir,dori;

  while(o1 >= 360.0)
    o1 -= 360.0;
  while(o1 < 0.0)
    o1 += 360.0;

  while(o2 >= 360.0)
    o2 -= 360.0;
  while(o2 < 0.0)
    o2 += 360.0;

  if (o1 >= o2)
    ddir = o1 - o2;
  else
    ddir = o2 - o1;

  if (ddir > 180.0)
    ddir = 360.0 - ddir;

  if (ddir > 90.0)
    dori = 180.0 - ddir;
  else
    dori = ddir;

  if (oriflag)
    return dori;
  else
    return ddir;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_CONN_GET_ORI_DIFF_ONE                        */
/*                                                                           */
/*****************************************************************************/
float mod_conn_get_ori_diff_one(c1,o2,oriflag)
     struct pop_cell *c1;
     float o2;
     int oriflag;  // 1-difference in ori, otherwise difference in dir
                   // WYETH - this use of 'oriflag' was never tested
{
  float o1,ddir,dori;

  // WYETH ATTRIB
  o1 = pop_cell_attrib_get_f(c1,"ori");
  //o1 = c1->ori;

  while(o1 >= 360.0)
    o1 -= 360.0;
  while(o1 < 0.0)
    o1 += 360.0;

  while(o2 >= 360.0)
    o2 -= 360.0;
  while(o2 < 0.0)
    o2 += 360.0;

  if (o1 >= o2)
    ddir = o1 - o2;
  else
    ddir = o2 - o1;

  if (ddir > 180.0)
    ddir = 360.0 - ddir;

  if (ddir > 90.0)
    dori = 180.0 - ddir;
  else
    dori = ddir;

  if (oriflag)
    return dori;
  else
    return ddir;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_CONN_COLOR_ANGLE                           */
/*                                                                           */
/*****************************************************************************/
void mod_conn_color_angle(theta,rr,gg,bb)
     float theta,*rr,*gg,*bb;
{
  int nc;
  float *cr,*cg,*cb,*ca;

  //  ***********
  //  ***********
  //  *** WYETH - THIS IS inefficient, given it is called many times in
  //  a loop when drawing populations of units ***
  //  ***********
  //  ***********

  /*** Specify colors for certain directions ***/
  nc = 5;

  cr = (float *)myalloc(nc*sizeof(float));
  cg = (float *)myalloc(nc*sizeof(float));
  cb = (float *)myalloc(nc*sizeof(float));
  ca = (float *)myalloc(nc*sizeof(float));

  ca[0] =   0.0;  cr[0] = 1.0; cg[0] = 0.0; cb[0] = 0.0;
  ca[1] =  90.0;  cr[1] = 1.0; cg[1] = 1.0; cb[1] = 0.0;
  ca[2] = 180.0;  cr[2] = 0.0; cg[2] = 1.0; cb[2] = 0.0;
  ca[3] = 270.0;  cr[3] = 0.0; cg[3] = 0.0; cb[3] = 1.0;
  ca[4] = 360.0;  cr[4] = 1.0; cg[4] = 0.0; cb[4] = 0.0;

  /*
    i0 = theta/(float)(nc-1);
    i1 = i0 + 1;
    printf("  theta = %f  ca[i0]= %f  ca[i1]= %f\n",theta,ca[i0],ca[i1]);*/

  // WYETH - Added July 2011, for split dir maps ("dir_1")
  while(theta >= 360.0)
    theta -= 360.0;

  *rr = lin_interp_farray(ca,cr,nc,theta);
  *gg = lin_interp_farray(ca,cg,nc,theta);
  *bb = lin_interp_farray(ca,cb,nc,theta);

  myfree(cr); myfree(cg); myfree(cb); myfree(ca);
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_CONN_COLOR_SF                            */
/*                                                                           */
/*****************************************************************************/
void mod_conn_color_sf(sf,rr,gg,bb)
     float sf;            // Spatial freq.  (cyc/deg)
     float *rr,*gg,*bb;   // Return color for color map
{
  int nc;
  float *cr,*cg,*cb,*ca;

  //  ***********
  //  ***********
  //  *** WYETH - THIS IS inefficient, given it is called many times in
  //  a loop when drawing populations of units ***
  //  ***********
  //  ***********

  // Specify colors for certain directions
  nc = 5;

  cr = (float *)myalloc(nc*sizeof(float));
  cg = (float *)myalloc(nc*sizeof(float));
  cb = (float *)myalloc(nc*sizeof(float));
  ca = (float *)myalloc(nc*sizeof(float));

  ca[0] =   0.50;  cr[0] = 1.0; cg[0] = 0.0; cb[0] = 1.0;
  ca[1] =   1.00;  cr[1] = 0.0; cg[1] = 0.0; cb[1] = 1.0;
  ca[2] =   2.00;  cr[2] = 0.0; cg[2] = 1.0; cb[2] = 0.0;
  ca[3] =   3.00;  cr[3] = 1.0; cg[3] = 1.0; cb[3] = 0.0;
  ca[4] =   6.00;  cr[4] = 1.0; cg[4] = 0.0; cb[4] = 0.0;

  /*
    i0 = theta/(float)(nc-1);
    i1 = i0 + 1;
    printf("  theta = %f  ca[i0]= %f  ca[i1]= %f\n",theta,ca[i0],ca[i1]);*/

  // WYETH - Added July 2011, for split dir maps ("dir_1")
  if (sf > 7.0)
    sf = 7.0;
  if (sf < 0.3)
    sf = 0.3;

  *rr = lin_interp_farray(ca,cr,nc,sf);
  *gg = lin_interp_farray(ca,cg,nc,sf);
  *bb = lin_interp_farray(ca,cb,nc,sf);

  myfree(cr); myfree(cg); myfree(cb); myfree(ca);
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_CONN_TEST_COLOR_MAP                          */
/*                                                                           */
/*****************************************************************************/
void mod_conn_test_color_map()
{
  int i,j;
  int xn,yn;
  float **r,**g,**b;

  xn = 360;
  yn = 32;

  r = get_zero_2d_farray(xn,yn);
  g = get_zero_2d_farray(xn,yn);
  b = get_zero_2d_farray(xn,yn);

  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++)
      mod_conn_color_angle((float)i,&(r[i][j]),&(g[i][j]),&(b[i][j]));

  write_ppm_6_rgb_data("zzz_color_test.ppm",xn,yn,r,g,b,0);
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_CONN_WRITE_MAP_TEXT                         */
/*                                                                           */
/*****************************************************************************/
void mod_conn_write_map_text(outfile,map,xn,yn,exitflag)
     char *outfile;
     float **map;
     int xn,yn;
     int exitflag;
{
  char tstr[SLEN];

  if ((strcmp(outfile,"null")==0) || (strcmp(outfile,"NULL")==0))
    return;

  printf("  MOD_CONN_WRITE_MAP_TEXT\n");

  //
  //  WYETH - Could have this query the user for exit after writing the map
  //

  remove_file(outfile);

  printf("    *** Writing map to text file:  %s\n",outfile);

  sprintf(tstr,"#  Map is %d columns by %d rows.\n",xn,yn);
  append_string_to_file(outfile,tstr);

  sprintf(tstr,
	  "#  1st column is written first, starting at lower left corner.\n");
  append_string_to_file(outfile,tstr);

  append_2d_farray(outfile,map,xn,yn,0); // 0-no transpose

  if (exitflag == 1){
    printf("    Exiting after writing map.\n\n");
    exit(0);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_CONN_WRITE_ORI_MAP                         */
/*                                                                           */
/*****************************************************************************/
void mod_conn_write_ori_map(outfile,orimap,xn,yn)
     char outfile[];
     float **orimap;
     int xn,yn;
{
  int i,j;
  float **r,**g,**b;

  printf("  MOD_CONN_WRITE_ORI_MAP\n");

  r = get_zero_2d_farray(xn,yn);
  g = get_zero_2d_farray(xn,yn);
  b = get_zero_2d_farray(xn,yn);
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      mod_conn_color_angle(orimap[i][j],&(r[i][j]),&(g[i][j]),&(b[i][j]));
    }
  }
  write_ppm_6_rgb_data(outfile,xn,yn,r,g,b,0);

  free_2d_farray(r,xn);
  free_2d_farray(g,xn);
  free_2d_farray(b,xn);
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_CONN_GET_ORI_MAP_OLD                        */
/*                                                                           */
/*  'dn' is the distance between pin-wheel centers.  Around the centers,     */
/*  the values change from 0 to 360 degrees.  So, for ori (not direction)    */
/*  the values should be divided by 2.                                       */
/*                                                                           */
/*****************************************************************************/
float **mod_conn_get_ori_map_old(mylogf,xn,yn,dn)
     char *mylogf;
     int xn,yn,dn;
{
  int i,j,k,l;
  int w,h,ix,iy;
  float fc,fx,fy,**orimap;
  char tstr[SLEN];

  mylog(mylogf,"  MOD_CONN_GET_ORI_MAP_OLD\n");

  w = xn/dn;  /* Width in units of pinwheel separations */
  h = yn/dn;  /* Height " */

  sprintf(tstr,"    Map is %d x %d pixels\n",xn,yn);
  mylog(mylogf,tstr);
  sprintf(tstr,"    Contains %d x %d columns\n",w,h);
  mylog(mylogf,tstr);

  fc = (float)(dn-1)/2.0;

  orimap = get_zero_2d_farray(xn,yn);
  for(i=0;i<w;i++){
    for(j=0;j<h;j++){  /* For each pinwheel */
      for(k=0;k<dn;k++){
	for(l=0;l<dn;l++){
	  ix = i*dn + k;
	  iy = j*dn + l;
	  fx = (float)k-fc;
	  fy = (float)l-fc;
	  if (i%2)
	    fx *= -1.0;
	  if (j%2)
	    fy *= -1.0;
	  orimap[ix][iy] = 180.0 + 180.0/M_PI * atan2(fy,fx); /* [0..360] */
	}
      }
    }
  }
  return orimap;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_CONN_GET_ORI_MAP                          */
/*                                                                           */
/*  'dx' and 'dy' are the x- and y-distances between pin-wheel centers.      */
/*  "ori" -  Around the centers, the values change from 0 to 360 degrees.    */
/*           So, for ori (not direction) the values should be divided by 2.  */
/*  "dir_1" -  Values go from 0 to 360 and should *not* be divided by 2.     */
/*             This is a split-dir map where ori-domains are split into      */
/*             opposite directions.                                          */
/*                                                                           */
/*****************************************************************************/
float **mod_conn_get_ori_map(mylogf,xn,yn,dx,dy,xphase,yphase,sphase,maptype)
     char *mylogf;
     int xn,yn;
     float dx,dy;   // Distance between pinwheel centers, in layer coords
     float xphase;  // Phase shift, 0 to 360, for positioning pinwheel ctr
     float yphase;  // Phase shift, 0 to 360, for positioning pinwheel ctr
     float sphase;  // Phase shift, 0 to 360, for spinning around pinwheel
     char *maptype; // 'ori', 'dir_1'
{
  int i,j,k,l;
  int ix,iy;
  //float fx,fy,**orimap,xphs,yphs,pht;
  float **orimap;
  double fx,fy,xphs,yphs,pht;
  char tstr[SLEN];
  
  mylog(mylogf,"  MOD_CONN_GET_ORI_MAP\n");
  
  if ((xphase < 0.0) || (xphase >= 360.0))
    mylog_exit(mylogf,"MOD_CONN_GET_ORI_MAP  xphase not in [0..360)\n");
  if ((yphase < 0.0) || (yphase >= 360.0))
    mylog_exit(mylogf,"MOD_CONN_GET_ORI_MAP  yphase not in [0..360)\n");

  sprintf(tstr,"    Map is %d x %d pixels\n",xn,yn);
  mylog(mylogf,tstr);

  // WYETH - changed May 1, 2009 - old way required perfect number of colums
  // Phase shifts in terms of layer coords
  xphs = (360.0 - xphase)/180.0 * dx; // x phase shift, keep it positive
  yphs = (360.0 - yphase)/180.0 * dy; // y phase shift, keep it positive

  //printf("xphs = %f  dx %f\n",xphs,dx);
  //printf("yphs = %f  dy %f\n",yphs,dy);

  orimap = get_zero_2d_farray(xn,yn);
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){  // For each point

      // Calculate coords on shifted frame
      fx = (double)i+0.5 + xphs;  // Add 0.5 to split margin evenly
      fy = (double)j+0.5 + yphs;

      // Calculate indices of nearest pinwheel center [0...npinwheels)
      ix = (int)(fx/dx);
      iy = (int)(fy/dy);

      // Adjust coords to be relative to closest pinwheel center
      fx = fx - ((double)ix*dx + dx/2.0);
      fy = fy - ((double)iy*dy + dy/2.0);

      if (ix%2)
	fx *= -1.0;
      if (iy%2)
	fy *= -1.0;

      pht = 180.0 + 180.0/M_PI * atan2(fy,fx); // [0..360]
      pht += sphase;  // Rotate around pinwheel
      while (pht < 0.0)
	pht += 360.0;
      while (pht >= 360.0)
	pht -= 360.0;
      orimap[i][j] = (float)pht;

      // WYETH - this catches some rounding differences between float/double
      if (orimap[i][j] >= 360.0)
	orimap[i][j] -= 360.0;

      if (strcmp(maptype,"dir_1")==0){
	//
	//  Make a split direction map, where values go from 0..360, and
	//  values that differ by 180 are next to each other.
	//
	orimap[i][j] /= 2.0;
	if (ix%2){
	  orimap[i][j] = orimap[i][j] + 180.0;  // Opposite direction
	}
      }

    }
  }

  //printf("---------------- dx = %f\n",dx);
  //mod_conn_write_ori_map("zzz.map.ppm",orimap,xn,yn);
  //exit(0);
  
  return orimap;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_CONN_ONODE_MAKE_MAP_ORI                       */
/*                                                                           */
/*****************************************************************************/
void mod_conn_onode_make_map_ori(mylogf,pm,xn,yn)
     char *mylogf;
     struct pop_map *pm;       // Map
     int xn,yn;
{
  int round;
  float phx,phy,php,cpcx,cpcy,ncolx,ncoly;
  char *maptype,*txtfile;

  if (pm->map2 != NULL)
    mylog_exit(mylogf,"MOD_CONN_ONODE_MAKE_MAP_ORI  map2 not NULL");

  ncolx   = onode_getpar_flt_dflt(pm->o,"ncol_x",2);
  ncoly   = onode_getpar_flt_dflt(pm->o,"ncol_y",ncolx);
  phx     = onode_getpar_flt_dflt(pm->o,"phase_x",0.0);
  phy     = onode_getpar_flt_dflt(pm->o,"phase_y",0.0);
  php     = onode_getpar_flt_dflt(pm->o,"phase_p",0.0);
  round   = onode_getpar_int_dflt(pm->o,"round_to_unit",0);
  txtfile = onode_getpar_chr_dflt(pm->o,"write_as_text",NULL);

  cpcx  = (float)xn/ncolx;
  cpcy  = (float)yn/ncoly;

  if (round){
    cpcx = (float)my_rint(cpcx);
    cpcy = (float)my_rint(cpcy);
  }

  maptype = onode_getpar_chr_exit(pm->o,"type");
  pm->map2 = mod_conn_get_ori_map(mylogf,xn,yn,cpcx,cpcy,phx,phy,php,maptype);
  pm->xn = xn;
  pm->yn = yn;
  myfree(maptype);

  if (txtfile != NULL){
    mod_conn_write_map_text(txtfile,pm->map2,xn,yn,1);  // 1-exit flag
    myfree(txtfile);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                        MOD_CONN_SET_LAYER_ORI_MAP                         */
/*                                                                           */
/*  Set 'ori' values for all cells in layer according to ori map.            */
/*                                                                           */
/*****************************************************************************/
void mod_conn_set_layer_ori_map(mylogf,lay,pm,oriflag)
     char *mylogf;
     struct pop_layer *lay;
     struct pop_map *pm;
     int oriflag;   // If 1 convert direction to orientation
{
  int i,j,k;
  int xn,yn,zn,xi,yi,ai_ori;
  float **map,x0,y0;
  struct pop_cell *pc;

  mylog(mylogf,"  MOD_CONN_SET_LAYER_ORI_MAP\n");

  //printf("    layer  %s\n",lay->name);

  ai_ori =  pop_cell_attrib_index_f(lay,"ori"); // Attribute index

  xn = pm->xn;
  yn = pm->yn;
  map = pm->map2;

  x0 = lay->x0;
  y0 = lay->y0;

  zn = lay->zn;
  for(i=0;i<lay->xn;i++){
    for(j=0;j<lay->yn;j++){
      if (j%2 == 0)
	xi = my_rint(x0 + lay->xf * (float)i);
      else
	xi = my_rint(x0 + lay->xf * (float)i + lay->oddxoff);
      yi = my_rint(y0 + lay->yf * (float)j);
	
      if ((xi < 0) || (xi >= xn) || (yi < 0) || (yi >= yn)){
	printf("xf = %f\n",lay->xf);
	printf("xi yi   %d %d\n",xi,yi);
	mylog_exit(mylogf,"MOD_CONN_SET_LAYER_ORI_MAP  index off of map\n");
      }

      for(k=0;k<zn;k++){  // For each cell in the first layer
	pc = &(lay->c[i][j][k]);

	if (oriflag == 1){
	  /*
	  if (map[xi][yi] >= 360.0)
	    printf("HERE WYETH -   ----------------- 360  %f\n",map[xi][yi]);
	  */
	  pc->attrib_f[ai_ori] = map[xi][yi] / 2.0;
	}else{
	  pc->attrib_f[ai_ori] = map[xi][yi];
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_CONN_GET_SF_W                             */
/*                                                                           */
/*  Assume the box is 1 x 1 unit, and dx,dy is a point on/within the box.    */
/*                                                                           */
/*    01__11       0    1                                                    */
/*     |  |                                                                  */
/*    00--10       1    0                                                    */
/*                                                                           */
/*****************************************************************************/
float mod_conn_get_sf_w(dx,dy,cflag,gamma,gam_thresh)
     float dx,dy;  // Distance in x,y from corner
     int cflag;
     float gamma;        // To apply a nonlinear remapping
     float gam_thresh;   // Keep a linear mapping below this threshold [0..1]
{
  float dx1,dy1,d00,d01,d10,d11,w,x,y;

  dx1 = 1.0 - dx;
  dy1 = 1.0 - dy;

  //
  //  Think of these as weights for the value represented at the four corners.
  //  These weights are ONE MINUS the distance from that corner.
  //
  d00 = 1.0 - sqrt(dx *dx  + dy *dy );
  d01 = 1.0 - sqrt(dx *dx  + dy1*dy1);
  d10 = 1.0 - sqrt(dx1*dx1 + dy *dy );
  d11 = 1.0 - sqrt(dx1*dx1 + dy1*dy1);

  // 
  //  Some weight could end up negative, because the maximum distance from
  //  a corner is SQRT(2).  Set any such weights to zero.
  // 
  if (d00 < 0.0) d00 = 0.0;
  if (d01 < 0.0) d01 = 0.0;
  if (d10 < 0.0) d10 = 0.0;
  if (d11 < 0.0) d11 = 0.0;

  if (cflag == 0)
    w = (d00 + d11) / (d00 + d01 + d10 + d11);
  else
    w = (d01 + d10) / (d00 + d01 + d10 + d11);

  //
  //  Apply a nonlinearity, if the exponent 'gamma' is not 1.0
  //
  if (gamma != 1.0){

    // w = pow(w,gamma);   THIS WAS THE OLD WAY - simple gamma correction

    if (w <= gam_thresh)
      ;  // Let this be a linear mapping from 0 up to a
    else{
      y = w - gam_thresh;
      x = pow(y/(1.0-gam_thresh),gamma);
      w = gam_thresh + x*(1-gam_thresh);  // x=0 --> a,  x=1  ==> 1
    }
  }

  return w;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_CONN_GET_SF_MAP                           */
/*                                                                           */
/*  Copied from ..._GET_ORI_MAP and modified to find distances between       */
/*  pinwheels.                                                               */
/*                                                                           */
/*****************************************************************************/
float **mod_conn_get_sf_map(mylogf,xn,yn,dx,dy,xphase,yphase,sphase,
			    sfmin,sfmax,gamma,gamthr)
     char *mylogf;
     int xn,yn;
     float dx,dy;   // Distance between pinwheel centers, in layer coords
     float xphase;  // Phase shift, 0 to 360, for positioning pinwheel ctr
     float yphase;  // Phase shift, 0 to 360, for positioning pinwheel ctr
     float sphase;  // Phase shift, 0 to 360, for spinning around pinwheel
     float sfmin;
     float sfmax;
     float gamma;
     float gamthr;
{
  int i,j;
  int ix,iy;
  float **sfmap,dsf,ddx,ddy;
  double fx,fy,xphs,yphs,pht;
  char tstr[SLEN];
  
  mylog(mylogf,"  MOD_CONN_GET_SF_MAP\n");

  dsf = sfmax - sfmin;
  
  if ((xphase < 0.0) || (xphase >= 360.0))
    mylog_exit(mylogf,"MOD_CONN_GET_SF_MAP  xphase not in [0..360)\n");
  if ((yphase < 0.0) || (yphase >= 360.0))
    mylog_exit(mylogf,"MOD_CONN_GET_SF_MAP  yphase not in [0..360)\n");

  sprintf(tstr,"    Map is %d x %d pixels\n",xn,yn);
  mylog(mylogf,tstr);

  // WYETH - changed May 1, 2009 - old way required perfect number of colums
  // Phase shifts in terms of layer coords
  xphs = (360.0 - xphase)/180.0 * dx; // x phase shift, keep it positive
  yphs = (360.0 - yphase)/180.0 * dy; // y phase shift, keep it positive

  //printf("xphs = %f  dx %f\n",xphs,dx);
  //printf("yphs = %f  dy %f\n",yphs,dy);

  sfmap = get_zero_2d_farray(xn,yn);
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){  // For each point

      // Calculate coords on shifted frame
      fx = (double)i+0.5 + xphs;  // Add 0.5 to split margin evenly
      fy = (double)j+0.5 + yphs;

      // Calculate indices of nearest pinwheel center [0...npinwheels)
      ix = (int)(fx/dx);
      iy = (int)(fy/dy);

      // Adjust coords to be relative to closest pinwheel center
      fx = fx - ((double)ix*dx + dx/2.0);
      fy = fy - ((double)iy*dy + dy/2.0);

      if (fx < 0.0){
	fx = dx + fx;   // fx is now positive offset w/i a box
	ix -= 1;        // Use coord of pinwheel to lower left
      }
      if (fy < 0.0){
	fy = dy + fy;   // fy is not positive offset w/i a box
	iy -= 1;        // Use coord of pinwheel to lower left
      }

      ddx = fx/dx;  // [0..1]
      ddy = fy/dy;  // [0..1]

      if ((ix+iy)%2 == 0){
	sfmap[i][j] = sfmin + dsf * mod_conn_get_sf_w(ddx,ddy,0,gamma,gamthr);
      }else{
	sfmap[i][j] = sfmin + dsf * mod_conn_get_sf_w(ddx,ddy,1,gamma,gamthr);
      }

      //
      //  WYETH LARS - Consider adding some noise (spread) to the SF values
      //  here.
      //

    }
  }

  /*
  {
    float x[100],y[100];

    for(i=0;i<100;i++){
      x[i] = (float)i / 99.0;
      //ddx = ddy = x[i];
      ddy = 0.0;
      ddx = x[i];
      y[i] = mod_conn_get_sf_w(ddx,ddy,1,gamma,gamthr);
    }
    append_farray_plot("zzz.gamtest.pl","y",y,100,1);

    exit(0);
  }
  */


  //printf("---------------- dx = %f\n",dx);
  //mod_conn_write_ori_map("zzz.map.ppm",sfmap,xn,yn);
  //exit(0);
   
  return sfmap;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_CONN_ONODE_MAKE_MAP_SF                       */
/*                                                                           */
/*****************************************************************************/
void mod_conn_onode_make_map_sf(mylogf,pm,xn,yn)
     char *mylogf;
     struct pop_map *pm;       // Map
     int xn,yn;
{
  int round;
  float phx,phy,php,cpcx,cpcy,ncolx,ncoly,sfmin,sfmax,gamma,gammt;
  char *txtfile;

  if (pm->map2 != NULL)
    mylog_exit(mylogf,"MOD_CONN_ONODE_MAKE_MAP_SF  map2 not NULL");

  ncolx   = onode_getpar_flt_dflt(pm->o,"ncol_x",2);
  ncoly   = onode_getpar_flt_dflt(pm->o,"ncol_y",ncolx);
  phx     = onode_getpar_flt_dflt(pm->o,"phase_x",0.0);
  phy     = onode_getpar_flt_dflt(pm->o,"phase_y",0.0);
  php     = onode_getpar_flt_dflt(pm->o,"phase_p",0.0);
  sfmin   = onode_getpar_flt_dflt(pm->o,"sf_min",1.0);
  sfmax   = onode_getpar_flt_dflt(pm->o,"sf_max",2.0);
  gamma   = onode_getpar_flt_dflt(pm->o,"gamma",1.0);
  gammt   = onode_getpar_flt_dflt(pm->o,"gam_thresh",0.0);
  txtfile = onode_getpar_chr_dflt(pm->o,"write_as_text",NULL);

  round = onode_getpar_int_dflt(pm->o,"round_to_unit",0);

  cpcx  = (float)xn/ncolx;
  cpcy  = (float)yn/ncoly;
  if (round){
    cpcx = (float)my_rint(cpcx);
    cpcy = (float)my_rint(cpcy);
  }

  pm->map2 = mod_conn_get_sf_map(mylogf,xn,yn,cpcx,cpcy,phx,phy,php,
				 sfmin,sfmax,gamma,gammt);
  pm->xn = xn;
  pm->yn = yn;

  //
  //  Renormalize the distribution of SFs
  //
  if (onode_item(pm->o,"renorm")==1){
    int i,j,k;
    int nn,*uisf,*sf_ilist,nsf,nu,*ucnt,**mapi,rerule,sfi,seed;
    float *norm_targ,dx,*r,*w,sf_lo,sf_hi,sf_new,sf_delta;
    char *refile;
    struct onode *renorm;

    renorm = onode_get_item_child(pm->o,"renorm");
    rerule = onode_getpar_int_dflt(pm->o,"renorm_rule",1);
    seed   = onode_getpar_int_dflt(pm->o,"renorm_seed",1777);
    refile = onode_getpar_chr_dflt(pm->o,"renorm_plot",NULL);

    if (seed > 0.0)
      seed *= -1.0;

    nn = renorm->nval;
    norm_targ = get_farray(nn);
    for(i=0;i<nn;i++){
      norm_targ[i] = atof(onode_get_nth_val(renorm,i));
    }
    //printf("  RENORM:  %d values\n",renorm->nval);

    //
    //  Get a list of all SF values in the map, multiply by 100, round to
    //  integer.  Also, keep a 2D 'mapi' of these values to look them up
    //  later.
    //
    nsf = xn*yn;
    sf_ilist = (int *)myalloc(nsf*sizeof(int));
    mapi = get_2d_iarray(xn,yn);
    k = 0;
    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){
	mapi[i][j] = my_rint(100.0 * pm->map2[i][j]);
	sf_ilist[k] = mapi[i][j];
	k += 1;
      }
    }

    //  Get the list of unique (int) SF values
    get_unique_count_iarray(sf_ilist,nsf,&uisf,&ucnt,&nu);
    w = i2farray(ucnt,nu);

    /*
    printf("%d unique values\n",nu);
    for(i=0;i<nu;i++)
      printf(" %d occurs %d times:  %f\n",uisf[i],ucnt[i],w[i]);
    */

    dx = (sfmax - sfmin) / (float)nn;
    r = distrib_box_match(norm_targ,nn,sfmin,dx,w,nu,refile);

    /*
    for(i=0;i<nu;i++)
      printf("  Right side value %d =   %f\n",i,r[i]);
    */

    //
    //  Now, change the SF values in the map.
    //
    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){
	sfi = mapi[i][j];
	k = get_index_search_iarray(uisf,nu,sfi);

	//  Find the lower and upper SF bounds for this region:
	//  *** WYETH - inefficient to keep recomputing these sf values...
	if (k > 0)
	  sf_lo = r[k-1];
	else
	  sf_lo = sfmin;
	sf_hi = r[k];
	sf_delta = sf_hi - sf_lo;

	if (rerule == 0){
	  sf_new = (sf_lo + sf_hi)/2.0;
	}else if (rerule == 1){
	  //
	  //  Pick uniform random
	  //
	  sf_new = sf_lo + sf_delta * myrand_util_ran2(&seed);
	}else{
	  mylog_exit(mylogf,"MOD_CONN_ONODE_MAKE_MAP_SF  Bad 'renorm_rule'");
	}
	pm->map2[i][j] = sf_new;
      }
    }

    myfree(uisf);
    myfree(sf_ilist);
    myfree(ucnt);
    free_2d_iarray(mapi,xn);
    myfree(norm_targ);
    myfree(r);
    myfree(w);
  }


  if (txtfile != NULL){
    mod_conn_write_map_text(txtfile,pm->map2,xn,yn,1);  // 1-exit flag
    myfree(txtfile);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_CONN_SET_LAYER_SF_MAP                         */
/*                                                                           */
/*  Set 'sf' values for all cells in layer according to sf map.              */
/*                                                                           */
/*****************************************************************************/
void mod_conn_set_layer_sf_map(mylogf,lay,pm)
     char *mylogf;
     struct pop_layer *lay;
     struct pop_map *pm;
{
  int i,j,k;
  int xn,yn,zn,xi,yi,ai_sf;
  float **map,x0,y0;
  struct pop_cell *pc;

  mylog(mylogf,"  MOD_CONN_SET_LAYER_SF_MAP\n");

  //printf("    layer  %s\n",lay->name);

  ai_sf =  pop_cell_attrib_index_f(lay,"sf"); // Attribute index

  //printf("    ai_SF attrib index  = %d\n",ai_sf);


  xn = pm->xn;
  yn = pm->yn;
  map = pm->map2;

  x0 = lay->x0;
  y0 = lay->y0;

  zn = lay->zn;
  for(i=0;i<lay->xn;i++){
    for(j=0;j<lay->yn;j++){
      if (j%2 == 0)
	xi = my_rint(x0 + lay->xf * (float)i);
      else
	xi = my_rint(x0 + lay->xf * (float)i + lay->oddxoff);
      yi = my_rint(y0 + lay->yf * (float)j);

      if ((xi < 0) || (xi >= xn) || (yi < 0) || (yi >= yn)){
	printf("xf = %f\n",lay->xf);
	printf("xi yi   %d %d\n",xi,yi);
	mylog_exit(mylogf,"MOD_CONN_SET_LAYER_SF_MAP  index off of map\n");
      }

      for(k=0;k<zn;k++){  // For each cell in the first layer
	pc = &(lay->c[i][j][k]);
	pc->attrib_f[ai_sf] = map[xi][yi];
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_CONN_ONODE_MAKE_MAP_OD                        */
/*                                                                           */
/*****************************************************************************/
void mod_conn_onode_make_map_od(mylogf,pm,xn,yn)
     char *mylogf;
     struct pop_map *pm;       // Map
     int xn,yn;
{
  int i,j;
  int round;
  float s,phx,cpcx,ncolx;
  char tstr[SLEN],*txtfile;

  mylog(mylogf,"  MOD_CONN_ONODE_MAKE_MAP_OD\n");

  if (pm->map2 != NULL)
    mylog_exit(mylogf,"MOD_CONN_ONODE_MAKE_MAP_OD  map2 not NULL");

  ncolx   = onode_getpar_flt_dflt(pm->o,"ncol_x",2);
  phx     = onode_getpar_flt_dflt(pm->o,"phase_x",0.0);
  round   = onode_getpar_int_dflt(pm->o,"round_to_unit",0);
  txtfile = onode_getpar_chr_dflt(pm->o,"write_as_text",NULL);

  while(phx < 0.0) // Make phx be positive
    phx += 360.0;

  cpcx  = (float)xn/ncolx;
  if (round)
    cpcx = (float)my_rint(cpcx);

  sprintf(tstr,"    Cells per OD column in x-axis:  %.4f\n",cpcx);
  mylog(mylogf,tstr);

  pm->map2 = get_zero_2d_farray(xn,yn);

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){  // For each point

      s = sin(M_PI * (float)i/cpcx + M_PI*phx/180.0);
      if (s > 0.0)
	pm->map2[i][j] = 1.0;
      else
	pm->map2[i][j] = -1.0;
    }
  }
  pm->xn = xn;
  pm->yn = yn;

  if (txtfile != NULL){
    mod_conn_write_map_text(txtfile,pm->map2,xn,yn,1);  // 1-exit flag
    myfree(txtfile);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_CONN_SET_LAYER_OD_MAP                        */
/*                                                                           */
/*  Set OD values for all cells in layer according to od map.                */
/*                                                                           */
/*****************************************************************************/
void mod_conn_set_layer_od_map(mylogf,lay,pm)
     char *mylogf;
     struct pop_layer *lay;
     struct pop_map *pm;
{
  int i,j,k;
  int xn,yn,zn,xi,yi,ai_ocdom;
  float **map,x0,y0;
  struct pop_cell *pc;

  mylog(mylogf,"  MOD_CONN_SET_LAYER_OD_MAP\n");

  xn = pm->xn;
  yn = pm->yn;
  map = pm->map2;

  x0 = lay->x0;
  y0 = lay->y0;

  ai_ocdom =  pop_cell_attrib_index_f(lay,"ocdom"); // Attribute index


  /*** WYETH - maps are no longer w.r.t. stim grid
  if (lay->area != NULL){
    x0 += (float)(lay->area->x0);  // Adjust to area origin
    y0 += (float)(lay->area->y0);
    }***/

  zn = lay->zn;
  for(i=0;i<lay->xn;i++){
    for(j=0;j<lay->yn;j++){
      if (j%2 == 0)
	xi = my_rint(x0 + lay->xf * (float)i);
      else
	xi = my_rint(x0 + lay->xf * (float)i + lay->oddxoff);
      yi = my_rint(y0 + lay->yf * (float)j);
	
      if ((xi < 0) || (xi >= xn) || (yi < 0) || (yi >= yn)){
	printf("xf = %f\n",lay->xf);
	printf("xi yi   %d %d\n",xi,yi);
	mylog_exit(mylogf,"MOD_CONN_SET_LAYER_OD_MAP  index off of map\n");
      }

      for(k=0;k<zn;k++){  // For each cell in the first layer
	pc = &(lay->c[i][j][k]);

	// WYETH ATTRIB
	pc->attrib_f[ai_ocdom] = map[xi][yi];
	//pc->ocdom = map[xi][yi];
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_CONN_ONODE_MAKE_MAP_RFX                       */
/*                                                                           */
/*****************************************************************************/
void mod_conn_onode_make_map_rfx(mylogf,pm,xn,yn,xf)
     char *mylogf;
     struct pop_map *pm;       // Map
     int xn,yn;
     float xf;
{
  int i,j;
  int round;
  float xx,s,phx,cpcx,ncolx,f,fcyc,offset;
  char ggstr[SLEN];

  mylog(mylogf,"  MOD_CONN_ONODE_MAKE_MAP_RFX\n");
  sprintf(ggstr,"    xn,yn = %d %d   xf = %f\n",xn,yn,xf);
  mylog(mylogf,ggstr);

  if (pm->map2 != NULL)
    mylog_exit(mylogf,"MOD_CONN_ONODE_MAKE_MAP_RFX  map2 not NULL");

  ncolx = onode_getpar_flt_dflt(pm->o,"ncol_x",2);
  phx   = onode_getpar_flt_dflt(pm->o,"phase_x",0.0);
  while(phx < 0.0) // Make phx be positive
    phx += 360.0;
  cpcx  = (float)xn/ncolx;  // Cells per column

  round = onode_getpar_int_dflt(pm->o,"round_to_unit",0);
  if (round)
    cpcx = (float)my_rint(cpcx);

  pm->map2 = get_zero_2d_farray(xn,yn);

  for(i=0;i<xn;i++){

    f = ((float)i/cpcx + phx/180.0);
    fcyc = f - (int)f;    // Fraction through cycle
    offset = (fcyc - 0.5) * (float)cpcx;

    //printf("offset = %f\n",offset);

    // OLD WYETH THIS WAS WRONG?
    // xx = 0.5 * ((float)i + offset); // 0.5 counteracts expansion

    xx = xf * ((float)i + offset);
    for(j=0;j<yn;j++)  // For each point
      pm->map2[i][j] = xx;
  }

  /*** Dump a row from the map
  {
    float *dd;

    dd = get_column_2d_farray(pm->map2,xn,yn,0);
    append_farray_plot("zzz.rfx.pl","rfx",dd,xn,1);

    exit(0); 
    }***/

  pm->xn = xn;
  pm->yn = yn;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_CONN_SET_LAYER_RFX_MAP                        */
/*                                                                           */
/*  Set 'rfx' values for all cells in layer according to map.                */
/*                                                                           */
/*****************************************************************************/
void mod_conn_set_layer_rfx_map(mylogf,lay,pm)
     char *mylogf;
     struct pop_layer *lay;
     struct pop_map *pm;
{
  int i,j,k;
  int xn,yn,zn,xi,yi;
  float **map,x0,y0;
  struct pop_cell *pc;

  mylog(mylogf,"  MOD_CONN_SET_LAYER_RFX_MAP\n");

  xn = pm->xn;
  yn = pm->yn;
  map = pm->map2;

  x0 = lay->x0;
  y0 = lay->y0;

  zn = lay->zn;
  for(i=0;i<lay->xn;i++){
    for(j=0;j<lay->yn;j++){
      if (j%2 == 0)
	xi = my_rint(x0 + lay->xf * (float)i);
      else
	xi = my_rint(x0 + lay->xf * (float)i + lay->oddxoff);
      yi = my_rint(y0 + lay->yf * (float)j);
	
      if ((xi < 0) || (xi >= xn) || (yi < 0) || (yi >= yn)){
	printf("xf = %f\n",lay->xf);
	printf("xi yi   %d %d\n",xi,yi);
	mylog_exit(mylogf,"MOD_CONN_SET_LAYER_RFX_MAP  index off of map\n");
      }

      for(k=0;k<zn;k++){  // For each cell in the first layer
	pc = &(lay->c[i][j][k]);
	pc->rfx = map[xi][yi] + lay->area->x0;
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                       MOD_CONN_SET_LAYER_DIR_MAP_FLAG                     */
/*                                                                           */
/*  Set 'ori' values to represent direction values according to 'dtype'.     */
/*                                                                           */
/*****************************************************************************/
void mod_conn_set_layer_dir_map_flag(mylogf,lay,dtype)
     char *mylogf;
     struct pop_layer *lay;
     char *dtype;            // 'ori', 'ori_invert', 'ori_alt_z'
{
  int i,j,k;
  int ai_ori;
  float tori;
  char ggstr[SLEN];
  struct pop_cell *pc;

  mylog(mylogf,"  MOD_CONN_SET_LAYER_DIR_MAP_FLAG\n");

  ai_ori =  pop_cell_attrib_index_f(lay,"ori"); // Attribute index

  for(i=0;i<lay->xn;i++){
    for(j=0;j<lay->yn;j++){
      for(k=0;k<lay->zn;k++){  // For each cell in the first layer

	pc = &(lay->c[i][j][k]);

	if (strcmp("ori",dtype) == 0){  // Do not change ori


	  /*
  printf("(mod_conn_ut) WYETH HERE -------------------------- NOT CHANGE\n");
	  tori = pc->attrib_f[ai_ori];
	  if (tori == 0.0)
	    pc->attrib_f[ai_ori] = 180.0;
  printf("(mod_conn_ut) WYETH HERE ------REMOVED------ NOT CHANGE\n");
	  */
	  ;  // No change
	}else if (strcmp("ori_invert",dtype) == 0){  // Set dir opposite of ori
	  tori = pc->attrib_f[ai_ori];

	  tori += 180.0;
	  while (tori >= 360.0)  // Won't happen if:  0 <= ori < 180
	    tori -= 360.0;

	  pc->attrib_f[ai_ori] = tori;

	}else if (strcmp("ori_alt_z",dtype) == 0){ // Alternate even/odd z
	  if (k%2 == 1){
	    tori = pc->attrib_f[ai_ori];
	    
	    tori += 180.0;
	    while (tori >= 360.0)  // Won't happen if:  0 <= ori < 180
	      tori -= 360.0;

	    pc->attrib_f[ai_ori] = tori;
	  }
	}else{
	  sprintf(ggstr,"  *** map_dir_flag %s\n",dtype);
	  mylog_exit(mylogf,"MOD_CONN_SET_LAYER_DIR_MAP_FLAG  Unknown value");
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_CONN_PICK_PAIRED_COORD                        */
/*                                                                           */
/*  Pick coordinates to pair with the given coordinates based on the offset  */
/*  and jitter.                                                              */
/*                                                                           */
/*****************************************************************************/
void mod_conn_pick_paired_coord(x1,y1,theta,dr,sdpara,sdorth,pseed,rx,ry)
     int x1,y1;            // Original coords
     float theta;          // Direction (degrees)
     float dr;             // Offset Distance
     float sdpara,sdorth;  // Noise SD
     int *pseed;           // pointer to seed for jitter
     int *rx,*ry;          // New connections, paired w/ orig
{
  int xi,yi;
  float x0,y0,x,y;
  float a,b,c,d;

  a =  cos(M_PI/180.0*theta); b = sin(M_PI/180.0*theta); // rotation matrix
  c = -sin(M_PI/180.0*theta); d = cos(M_PI/180.0*theta);

  x0 = dr + sdpara * nr_util_gasdev(pseed); // Noise along length of vector
  y0 = sdorth * nr_util_gasdev(pseed);      // Noise orthogonal to vector
  
  x = a*x0 + c*y0;  // Rotated coord
  y = b*x0 + d*y0;
  
  xi = my_rint(x) + x1;
  yi = my_rint(y) + y1;
  
  *rx = xi;
  *ry = yi;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_CONN_WRITE_CONN_TXT                         */
/*                                                                           */
/*****************************************************************************/
void mod_conn_write_conn_txt(mylogf,pc,outfile)
     char *mylogf;
     struct pconn *pc;
     char *outfile;
{
  FILE *fopen(),*fout;
  int i,j,k,l;
  struct pconn_el *pce;

  printf("  MOD_CONN_WRITE_CONN_TXT\n");

  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** File %s\n",outfile);
    mylog_exit(mylogf,"MOD_CONN_WRITE_CONN_TEXT  Cannot open outfile");
  }
  
  fprintf(fout,"FTYPE %d\n",pc->ftype);
  fprintf(fout,"NAME1 %s\n",pc->name1);
  fprintf(fout,"NAME2 %s\n",pc->name2);
  fprintf(fout,"SYNTYPE %s\n",pc->syntype);
  fprintf(fout,"%d %d %d\n",pc->xn,pc->yn,pc->zn);

  for(i=0;i<pc->xn;i++){
    for(j=0;j<pc->yn;j++){
      for(k=0;k<pc->zn;k++){
	pce = &(pc->c[i][j][k]);
	fprintf(fout,"N %d\n",pce->nsyn);
	for(l=0;l<pce->nsyn;l++){
	  fprintf(fout,"%d %d %d %.2f\n",pce->xi[l],pce->yi[l],pce->zi[l],
		  pce->w[l]);
	}
      }
    }
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_CONN_WRITE_CONN_T1                          */
/*                                                                           */
/*  Write coordinates and weight only.  All synapses must have same type.    */
/*                                                                           */
/*****************************************************************************/
void mod_conn_write_conn_t1(mylogf,lay1,lay2,syntype,outfile)
     char *mylogf;
     struct pop_layer *lay1;  // From here
     struct pop_layer *lay2;  // to here
     int syntype;             // All connections must be of this type
     char outfile[];
{
  FILE *fopen(),*fout;
  int i,j,k;
  int nsyn,ftype,n,xn,yn,zn;
  struct pop_cell *c;
  struct pop_syn *s,*t;
  char tstr[SLEN];

  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** File %s\n",outfile);
    exit_error("MOD_CONN_WRITE_CONN_T1","Cannot open file");
  }
  printf("    Writing %s to %s connections to %s\n",lay1->name,lay2->name,
	 outfile);

  i = 1;
  fwrite((char *)&i,sizeof(int),1,fout);  // For endian

  strcpy(tstr,"wmpop");
  carray_nchar_write(tstr,fout);      // WM Population
  strcpy(tstr,"conn");
  carray_nchar_write(tstr,fout);      // type of saved object
  ftype = 1;
  fwrite(&ftype,sizeof(int),1,fout);  // format of saved object

  carray_nchar_write(lay1->name,fout);
  carray_nchar_write(lay2->name,fout);
  fwrite(&(lay1->xn),sizeof(int),1,fout);
  fwrite(&(lay1->yn),sizeof(int),1,fout);
  fwrite(&(lay1->zn),sizeof(int),1,fout);

  fwrite(&syntype,sizeof(int),1,fout);  // Synapse type

  xn = lay1->xn;
  yn = lay1->yn;
  zn = lay1->zn;

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	c = &(lay1->c[i][j][k]);
	nsyn = pop_cell_count_syn_from_cell_to_layer(lay2->name,c);
	fwrite(&nsyn,sizeof(int),1,fout); // Number of synapses

	s = pop_cell_get_next_syn_layer_name(c->out,0,lay2->name);
	n = 0;
	while(s != NULL){
	  if (s->stype != syntype)
	    exit_error("MOD_CONN_WRITE_CONN_T1","Synapse type error");

	  fwrite(&(s->post->layx),sizeof(int),1,fout);  // x,y,z
	  fwrite(&(s->post->layy),sizeof(int),1,fout);
	  fwrite(&(s->post->layz),sizeof(int),1,fout);
	  fwrite(&(s->w),sizeof(float),1,fout);         // weight
	  
	  n += 1;
	  t = s->pre_next;
	  s = pop_cell_get_next_syn_layer_name(t,0,lay2->name);
	}
	if (n != nsyn)
	  exit_error("MOD_CONN_WRITE_CONN_T1","Synapse count error");
      }
    }
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_CONN_READ_CONN_T1_FP                         */
/*                                                                           */
/*  Write coordinates and weight only.  All synapses must have same type.    */
/*                                                                           */
/*****************************************************************************/
struct pconn *mod_conn_read_conn_t1_fp(mylogf,fin,revflag)
     char *mylogf;
     FILE *fin;
     int revflag;
{
  int i,j,k,l;
  int nread,nsyn,xn,yn,zn;
  struct pconn *pc;
  struct pconn_el *pce;

  pc = (struct pconn *)myalloc(sizeof(struct pconn));
  pc->ftype = 1;

  /*** Layer names ***/
  carray_nchar_read(&(pc->name1),fin,revflag);
  carray_nchar_read(&(pc->name2),fin,revflag);

  /*** Layer1 size ***/
  nread = revflag_fread(&xn,sizeof(int),1,fin,revflag);
  nread = revflag_fread(&yn,sizeof(int),1,fin,revflag);
  nread = revflag_fread(&zn,sizeof(int),1,fin,revflag);

  pc->xn = xn;
  pc->yn = yn;
  pc->zn = zn;

  /*** Synapse type ***/
  nread = revflag_fread(&(pc->syntype),sizeof(int),1,fin,revflag);

  pc->c = (struct pconn_el ***)myalloc(xn*sizeof(struct pconn_el **));
  for(i=0;i<xn;i++){
    pc->c[i] = (struct pconn_el **)myalloc(yn*sizeof(struct pconn_el *));
    for(j=0;j<yn;j++){
      pc->c[i][j] = (struct pconn_el *)myalloc(zn*sizeof(struct pconn_el));
      for(k=0;k<zn;k++){
	
	/*** Number of synapses ***/
	nread = revflag_fread(&nsyn,sizeof(int),1,fin,revflag);

	pce = &(pc->c[i][j][k]);
	pce->nsyn = nsyn;
	pce->xi = (int *)myalloc(nsyn*sizeof(int));
	pce->yi = (int *)myalloc(nsyn*sizeof(int));
	pce->zi = (int *)myalloc(nsyn*sizeof(int));
	pce->w = (float *)myalloc(nsyn*sizeof(float));
	
	for(l=0;l<nsyn;l++){
	  nread = revflag_fread(&(pce->xi[l]),sizeof(int),1,fin,revflag);
	  nread = revflag_fread(&(pce->yi[l]),sizeof(int),1,fin,revflag);
	  nread = revflag_fread(&(pce->zi[l]),sizeof(int),1,fin,revflag);
	  nread = revflag_fread(&(pce->w[l]),sizeof(float),1,fin,revflag);
	}
      }
    }
  }
  return pc;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_CONN_READ_CONN_T1                          */
/*                                                                           */
/*  Write coordinates and weight only.  All synapses must have same type.    */
/*                                                                           */
/*****************************************************************************/
void mod_conn_read_conn_t1(mylogf,lay1,lay2,infile,axonv,axondt,inindex)
     char *mylogf;
     struct pop_layer *lay1;  // From here
     struct pop_layer *lay2;  // to here
     char infile[];
     float axonv;             // Axon velocity, if > 0.0 (m/s)
     float axondt;            // Constant delay (s)
     short inindex;           // Index into population 'inlist'
{
  FILE *fopen(),*fin;
  int i,j,k,l;
  int revflag,nread,ftype,syntype;
  int nsyn,xn,yn,zn,xn2,yn2,zn2,xi,yi,zi;
  float w,dist,tdelay;
  double spum;
  char *tstr,*name1,*name2,ts[SLEN];
  struct pop_cell *c1,*c2;

  if (axonv > 0.0)
    spum = 1.0 / (axonv * 1000000.0);  // Seconds per Micron
  else
    spum = 0.0;

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** File %s\n",infile);
    exit_error("MOD_CONN_READ_CONN_T1","Cannot open file");
  }

  sprintf(ts,"    Reading %s to %s connections from %s\n",lay1->name,
	  lay2->name,infile);
  mylog(mylogf,ts);

  revflag = 0; // Endian check
  nread = revflag_fread((char *)&i,sizeof(int),1,fin,revflag);
  if (i != 1)
    revflag = 1;

  carray_nchar_read(&tstr,fin,revflag);
  if (strcmp(tstr,"wmpop")!=0)
    exit_error("MOD_CONN_READ_CONN_T1","Expecting wmpop");
  myfree(tstr);

  carray_nchar_read(&tstr,fin,revflag);
  if (strcmp(tstr,"conn")!=0)
    exit_error("MOD_CONN_READ_CONN_T1","Expecting conn");
  myfree(tstr);

  nread = revflag_fread((char *)&ftype,sizeof(int),1,fin,revflag);
  if (ftype != 1)
    exit_error("MOD_CONN_READ_CONN_T1","Expecting ftype 1");

  // Layer names
  carray_nchar_read(&name1,fin,revflag);
  carray_nchar_read(&name2,fin,revflag);
  if (strcmp(name1,lay1->name)!=0)
    printf("      layer1 names differ:  %s  %s\n",name1,lay1->name);
  if (strcmp(name2,lay2->name)!=0)
    printf("      layer2 names differ:  %s  %s\n",name2,lay2->name);
  myfree(name1);
  myfree(name2);
  
  // Layer1 size
  nread = revflag_fread(&xn,sizeof(int),1,fin,revflag);
  nread = revflag_fread(&yn,sizeof(int),1,fin,revflag);
  nread = revflag_fread(&zn,sizeof(int),1,fin,revflag);
  if ((xn != lay1->xn)||(yn != lay1->yn)||(zn != lay1->zn))
    exit_error("MOD_CONN_READ_CONN_T1","Layer1 size differs");

  // Synapse type
  nread = revflag_fread(&syntype,sizeof(int),1,fin,revflag);

  xn2 = lay2->xn;
  yn2 = lay2->yn;
  zn2 = lay2->zn;
  
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	c1 = &(lay1->c[i][j][k]);

	// Number of synapses
	nread = revflag_fread(&nsyn,sizeof(int),1,fin,revflag);

	for(l=0;l<nsyn;l++){
	  nread = revflag_fread(&xi,sizeof(int),1,fin,revflag);
	  nread = revflag_fread(&yi,sizeof(int),1,fin,revflag);
	  nread = revflag_fread(&zi,sizeof(int),1,fin,revflag);
	  nread = revflag_fread(&w,sizeof(float),1,fin,revflag);
	  if ((xi < 0) || (xi >= xn2) ||(yi < 0) || (yi >= yn2) ||
	      (zi < 0) || (zi >= zn2))
	    exit_error("MOD_CONN_READ_CONN_T1","Post-syn index error");
	  
	  c2 = &(lay2->c[xi][yi][zi]);

	  if (syntype == 11)
	    c2->lgn_n += 1;  // 2013 Mar (needed to allow saving LGN input)

	  if (axonv > 0.0){ // Distance dependent delay
	    dist = mod_conn_get_distance_cells_layer(c1,c2); // dist (um)
	    tdelay = axondt + dist*spum;
	  }else
	    tdelay = axondt; // constant delay (s)

	  pop_cell_add_synapse(c1,c2,syntype,w,tdelay,inindex);
	}
      }
    }
  }
  fclose(fin);
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_CONN_WRITE_CONN_T0                         */
/*                                                                           */
/*  Write LGN connections for the specified layer.                           */
/*                                                                           */
/*  HEADER FORMAT                                                            */
/*                                                                           */
/*  1             - (int) for endian                                         */
/*  "wmpop"       - (nstring) int for length, followed by characters         */
/*  "conn"        - (nstring) int for length, followed by characters         */
/*  0             - (int) format of saved object, 0=LGN                      */
/*  <layname>     - (nstring) post-synaptic layer name                       */
/*  <xn><yn><zn>  - (int,int,int) post-syn layer size                        */
/*  <syntype>     - (int) synapse type                                       */
/*                                                                           */
/*  BODY RECORD - One for each cell in 'lay1'                                */
/*                                                                           */
/*  <nsyn>        - (int) number of synapses from LGN                        */
/*  <x1><y1><w1>  - (int,int,float) coords and weight                        */
/*  ...           - REPEAT FOR EACH SYNAPSE                                  */
/*                                                                           */
/*****************************************************************************/
void mod_conn_write_conn_t0(mylogf,lay1,in_i,syntype,outfile)
     char *mylogf;
     struct pop_layer *lay1;  // From here
     int in_i;                // input index
     int syntype;             // All connections must be of this type
     char outfile[];
{
  FILE *fopen(),*fout;
  int i,j,k,si;
  int nsyn,ftype,n,xn,yn,zn,*cx,*cy;
  float w;
  struct pop_cell *c;
  char tstr[SLEN];
  struct pop_lgn_in *ln;

  //
  //  WYETH - MUST UPDATE FOR NEWLGN
  //  WYETH - MUST UPDATE FOR NEWLGN
  //  WYETH - MUST UPDATE FOR NEWLGN
  //
  //mylog_exit(mylogf,"MOD_CONN_WRITE_CONN_T0  Check for lgnnew compatible.");


  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** File %s\n",outfile);
    exit_error("MOD_CONN_WRITE_CONN_T0_LGN","Cannot open file");
  }
  printf("    Writing LGN inputs for %s to %s\n",lay1->name,outfile);

  i = 1;
  fwrite((char *)&i,sizeof(int),1,fout); // For endian
  
  strcpy(tstr,"wmpop");
  carray_nchar_write(tstr,fout);         // WM Population
  strcpy(tstr,"conn");
  carray_nchar_write(tstr,fout);         // type of saved object
  ftype = 0;
  fwrite(&ftype,sizeof(int),1,fout);     // format of saved object, 0=LGN

  carray_nchar_write(lay1->name,fout);
  fwrite(&(lay1->xn),sizeof(int),1,fout);
  fwrite(&(lay1->yn),sizeof(int),1,fout);
  fwrite(&(lay1->zn),sizeof(int),1,fout);

  fwrite(&syntype,sizeof(int),1,fout);  // Synapse type

  // *** WYETH - We should probably change this to write the actual weights
  // *** WYETH - We should probably change this to write the actual weights
  // *** WYETH - We should probably change this to write the actual weights
  w = 1.0;  // All weights are 1 for now

  xn = lay1->xn;
  yn = lay1->yn;
  zn = lay1->zn;

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	c = &(lay1->c[i][j][k]);

	if (c->lgnin == NULL){
	  //
	  //  THIS probably means that there is irregular LGN input
	  //
	  printf("(mod_conn_util) WYETH HERE 001  - Must Debug\n");
	  printf(" in_i = %d\n",in_i);
	  printf(" lgn_n = %d\n",c->lgn_n);

	  exit_error("WYETH HERE","THIS IS NULL");
	}

	if (in_i < c->lgn_n){
	  ln = c->lgnin[in_i];
	}else
	  exit_error("MOD_CONN_WRITE_CONN_TO_LGN","No LGN connections");

	// ON inputs
	nsyn = ln->cn1;
	cx = ln->cx1;
	cy = ln->cy1;

	fwrite(&nsyn,sizeof(int),1,fout); // Number of synapses
	for(si=0;si<nsyn;si++){
	  fwrite(&(cx[si]),sizeof(int),1,fout);  // x,y
	  fwrite(&(cy[si]),sizeof(int),1,fout);
	  fwrite(&w,sizeof(float),1,fout);      // weight
	}

	// OFF inputs
	nsyn = ln->cn0;
	cx = ln->cx0;
	cy = ln->cy0;
	fwrite(&nsyn,sizeof(int),1,fout); // Number of synapses
	for(si=0;si<nsyn;si++){
	  fwrite(&(cx[si]),sizeof(int),1,fout);  // x,y
	  fwrite(&(cy[si]),sizeof(int),1,fout);
	  fwrite(&w,sizeof(float),1,fout);      // weight
	}
      }
    }
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_CONN_READ_CONN_T0                          */
/*                                                                           */
/*  Read LGN inputs from file.                                               */
/*                                                                           */
/*****************************************************************************/
void mod_conn_read_conn_t0(mylogf,lay0,lay1,mech,infile,normw)
     char *mylogf;
     struct pop_layer *lay0;  // LGN layer
     struct pop_layer *lay1;  // Inputs to this layer
     struct pop_mech *mech;   // post-syn receptor mechanism
     char infile[];           // From this file
     float normw;             // Weight normalization *** WYETH NOT IN FILE
{
  FILE *fopen(),*fin;
  int i,j,k,l;
  int revflag,nread,ftype,syntype,nsamp;
  int nsyn,xn,yn,zn,*cx,*cy;
  int *tcx0,*tcy0,tcn0,*tcx1,*tcy1,tcn1;
  float w;
  char *tstr,*name1,ts[SLEN];
  struct pop_cell *c1;
  struct pop_lgn_in *ln;

  //
  //  WYETH - MUST UPDATE FOR NEWLGN
  //  WYETH - MUST UPDATE FOR NEWLGN
  //  WYETH - MUST UPDATE FOR NEWLGN
  //
  //mylog_exit(mylogf,"MOD_CONN_READ_CONN_T0  Check for lgnnew compatible.");



  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** File %s\n",infile);
    exit_error("MOD_CONN_READ_CONN_T0","Cannot open file");
  }
  sprintf(ts,"    Reading LGN inputs for %s from %s\n",lay1->name,infile);
  mylog(mylogf,ts);

  //
  //  See writing routine above for comments on format.
  //

  revflag = 0; // Endian check
  nread = revflag_fread((char *)&i,sizeof(int),1,fin,revflag);
  if (i != 1)
    revflag = 1;

  carray_nchar_read(&tstr,fin,revflag);
  if (strcmp(tstr,"wmpop")!=0)
    exit_error("MOD_CONN_READ_CONN_T0","Expecting wmpop");
  myfree(tstr);

  carray_nchar_read(&tstr,fin,revflag);
  if (strcmp(tstr,"conn")!=0)
    exit_error("MOD_CONN_READ_CONN_T0","Expecting conn");
  myfree(tstr);

  nread = revflag_fread((char *)&ftype,sizeof(int),1,fin,revflag);
  if (ftype != 0)
    exit_error("MOD_CONN_READ_CONN_T0","Expecting ftype 0");

  // Layer names
  carray_nchar_read(&name1,fin,revflag);
  if (strcmp(name1,lay1->name)!=0)
    printf("      layer1 names differ:  %s  %s\n",name1,lay1->name);
  myfree(name1);
  
  // Layer1 size
  nread = revflag_fread(&xn,sizeof(int),1,fin,revflag);
  nread = revflag_fread(&yn,sizeof(int),1,fin,revflag);
  nread = revflag_fread(&zn,sizeof(int),1,fin,revflag);
  if ((xn != lay1->xn)||(yn != lay1->yn)||(zn != lay1->zn))
    exit_error("MOD_CONN_READ_CONN_T0","Layer size differs");

  // Synapse type
  nread = revflag_fread(&syntype,sizeof(int),1,fin,revflag);

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	c1 = &(lay1->c[i][j][k]);

	//c1->cxn = cxn; // Must set grid dimensions for each cell
	//c1->cyn = cyn;

	// ON inputs
	nread = revflag_fread(&nsyn,sizeof(int),1,fin,revflag);
	cx = (int *)myalloc(nsyn*sizeof(int));
	cy = (int *)myalloc(nsyn*sizeof(int));
	for(l=0;l<nsyn;l++){
	  nread = revflag_fread(&(cx[l]),sizeof(int),1,fin,revflag);
	  nread = revflag_fread(&(cy[l]),sizeof(int),1,fin,revflag);
	  nread = revflag_fread(&w,sizeof(float),1,fin,revflag);
	  // *** WYETH - note, all 'w' are set to 1 by '...write' above
	}
	tcx1 = cx;
	tcy1 = cy;
	tcn1 = nsyn;

	// OFF inputs
	nread = revflag_fread(&nsyn,sizeof(int),1,fin,revflag);
	cx = (int *)myalloc(nsyn*sizeof(int));
	cy = (int *)myalloc(nsyn*sizeof(int));
	for(l=0;l<nsyn;l++){
	  nread = revflag_fread(&(cx[l]),sizeof(int),1,fin,revflag);
	  nread = revflag_fread(&(cy[l]),sizeof(int),1,fin,revflag);
	  nread = revflag_fread(&w,sizeof(float),1,fin,revflag);
	  // *** WYETH - note, all 'w' are set to 1 by '...write' above
	}
	tcx0 = cx;
	tcy0 = cy;
	tcn0 = nsyn;

	//
	//  Compute weight
	//
	if (normw > 0.0){
	  nsamp = tcn0 + tcn1;
	  w = normw / (float)nsamp; // WYETH, LGNW
	  //printf("w = %f\n",w);
	}else{
	  w = 1.0;
	}

	// '0' set for 'nposs'
	pop_cell_add_lgn_input(c1,lay0,mech,tcx0,tcy0,tcn0,tcx1,tcy1,tcn1,0,
			       w);
			       //-1.0);
      }
    }
  }
  fclose(fin);
}

/**************************************-**************************************/
/*                                                                           */
/*                             MOD_CONN_READ_CONN                            */
/*                                                                           */
/*  Read a connection file into a struct.                                    */
/*                                                                           */
/*****************************************************************************/
struct pconn *mod_conn_read_conn(mylogf,infile)
     char *mylogf;
     char infile[];
{
  FILE *fopen(),*fin;
  int i;
  int revflag,nread,ftype;
  char *tstr,ts[SLEN];
  struct pconn *pc;

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** File %s\n",infile);
    exit_error("MOD_CONN_READ_CONN","Cannot open file");
  }
  sprintf(ts,"    MOD_CONN_READ_CONN Reading connections from %s\n",infile);
  mylog(mylogf,ts);

  revflag = 0; // Endian check
  nread = revflag_fread((char *)&i,sizeof(int),1,fin,revflag);
  if (i != 1)
    revflag = 1;

  carray_nchar_read(&tstr,fin,revflag);
  if (strcmp(tstr,"wmpop")!=0)
    exit_error("MOD_CONN_READ_CONN","Expecting wmpop");
  myfree(tstr);

  carray_nchar_read(&tstr,fin,revflag);
  if (strcmp(tstr,"conn")!=0)
    exit_error("MOD_CONN_READ_CONN","Expecting conn");
  myfree(tstr);

  nread = revflag_fread((char *)&ftype,sizeof(int),1,fin,revflag);
  if (ftype == 0)
    mylog_exit(mylogf,"    ftype 0 not implemented yet\n");
  else if (ftype == 1)
    pc = mod_conn_read_conn_t1_fp(mylogf,fin,revflag);
  else
    mylog_exit(mylogf,"    Unknown file type, %d\n",ftype);

  fclose(fin);

  return pc;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_CONN_CONN_TO_TEXT                           */
/*                                                                           */
/*****************************************************************************/
void mod_conn_conn_to_text(mylogf,infile,outfile)
     char *mylogf;
     char infile[];
     char outfile[];
{
  struct pconn *pc;

  pc = mod_conn_read_conn(mylogf,infile);
  mod_conn_write_conn_txt(mylogf,pc,outfile);

  /*** WYETH should free 'pc' ***/
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_CONN_IS_LGN_CONNECTED                        */
/*                                                                           */
/*  Return 1 if the lgn ON or OFF cell is connected to the unit.             */
/*                                                                           */
/*****************************************************************************/
int mod_conn_is_lgn_connected(mylogf,c,xi,yi,onflag)
     char *mylogf;
     struct pop_cell *c;  // Post-syn cell
     int xi,yi;           // LGN coordinates
     int onflag;          // 1-ON LGN, 0-OFF LGN
{
  int i;
  int n,*x,*y;
  struct pop_lgn_in *ln;

  // WYETH HERE - assuming zero-th LGN input

  if (c->lgn_n > 0)
    ln = c->lgnin[0];
  else
    exit_error("MOD_CONN_IS_LGN_CONNECTED","No LGN inputs");

  if ((xi < 0) || (xi >= ln->lay->xn) || (yi < 0) || (yi >= ln->lay->yn))
    mylog_exit(mylogf,"MOD_CONN_IS_LGN_CONNECTED  Index error");
  
  if (onflag){
    n = ln->cn1;
    x = ln->cx1;
    y = ln->cy1;
  }else{
    n = ln->cn0;
    x = ln->cx0;
    y = ln->cy0;
  }

  for(i=0;i<n;i++){
    if (x[i] == xi)
      if (y[i] == yi)
	return 1;
  }

  return 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_CONN_ASTAT_LGN_POSTSYN                        */
/*                                                                           */
/*  Compute and write statistics for LGN projections.                        */
/*                                                                           */
/*****************************************************************************/
void mod_conn_astat_lgn_postsyn(mpt,r)
     struct pop_top *mpt;          // Global model structure
     struct response_struct *r;    // Response params
{
  int i,j,k,l;
  int i0,i1,j0,j1,kmax,lmax,maxsyn,xn,yn,ai_ocdom,iocdom;
  int n,cn,*cx,*cy,xi,yi,cnt,onflag;
  int ***lgnx,***lgny,**lgnn,*tx,*ty,eyeflag;
  float umx,umy,*fcx,*fcy,xmu,ymu,xsd,ysd,dxf,dx,dy,dd,dmax;
  char tstr[SLEN3],*lname,nstr[SLEN],outfile[SLEN];
  struct pop_cell *c;
  struct pop_layer *pl;    // post-syn layer to examine
  struct pop_lgn_in *ln;

  printf("  MOD_CONN_ASTAT_LGN_POSTSYN\n");

  // Get a pointer to the post-syn layer, specified in the .rsp file
  lname = paramfile_get_char_param_or_exit(r->ppl,"astat_lgn_proj_layer");
  pl = pop_cell_get_layer_pointer(mpt->lay,mpt->nlay,lname);
  myfree(lname);

  ai_ocdom =  pop_cell_attrib_index_f(pl,"ocdom"); // Attribute index

  // Get any scaling factor for X-axis distance, e.g., for asym. CMF
  dxf = paramfile_get_float_param_default(r->ppl,"astat_scale_factor_x",1.0);

  printf("dxf = %f\n",dxf);

  // Determine whether to process Left or Right (-1 or 1)
  eyeflag = paramfile_get_int_param_default(r->ppl,"astat_lgn_proj_eye",-1);

  // Determine whether to process ON or OFF units
  onflag = paramfile_get_int_param_default(r->ppl,"astat_lgn_proj_onflag",1);
  if (onflag){
    if (eyeflag == -1)
      sprintf(nstr," ON_L %s",pl->name);
    else
      sprintf(nstr," ON_R %s",pl->name);
  }else{
    if (eyeflag == -1)
      sprintf(nstr,"OFF_L %s",pl->name);
    else
      sprintf(nstr,"OFF_R %s",pl->name);
  }
  
  // Make outfile
  sprintf(outfile,"%s.astat",mpt->out_prefix);
  remove_file(outfile);
  printf("    Writing statistics to %s\n",outfile);

  maxsyn = 1000;  // Arbitrary constant, to limit storage for synapses

  //
  //  Build lists of all x and y coordinates that each LGN cell contacts
  //  in the post-syn layer.  Do this by scanning through the input lists
  //  of all cells in the post-syn layer.
  //
  xn = mpt->xn;  // Dimensions of the LGN layer
  yn = mpt->yn;
  lgnx = get_3d_iarray(xn,yn,maxsyn);
  lgny = get_3d_iarray(xn,yn,maxsyn);
  lgnn = get_zero_2d_iarray(xn,yn);

  for(i=0;i<pl->xn;i++){  // For each cell in the post-syn layer
    for(j=0;j<pl->yn;j++){
      for(k=0;k<pl->zn;k++){
	c = &(pl->c[i][j][k]);

	// WYETH - ASSUMING LGNIN[0]
	if (c->lgn_n > 0)
	  ln = c->lgnin[0];
	else
	  exit_error("MOD_CONN_ASTAT_LGN_POSTSYN","No LGN inputs");


	// WYETH ATTRIB
	iocdom = my_rint(c->attrib_f[ai_ocdom]);
	//iocdom = my_rint(c->ocdom);

	if (iocdom == eyeflag){ // -1 or 1 (L or R)
	  if (onflag){
	    cn = ln->cn1;
	    cx = ln->cx1;
	    cy = ln->cy1;
	  }else{
	    cn = ln->cn0;
	    cx = ln->cx0;
	    cy = ln->cy0;
	  }
	  for(l=0;l<cn;l++){
	    xi = cx[l];
	    yi = cy[l];
	    cnt = lgnn[xi][yi];
	    lgnx[xi][yi][cnt] = i;
	    lgny[xi][yi][cnt] = j;
	    lgnn[xi][yi] += 1;
	    if (lgnn[xi][yi] >= maxsyn)
	      exit_error("MOD_CONN_ASTAT_LGN_POSTSYN","Maxsyn exceeded");
	  }
	}
      }
    }
  }

  // Determine window
  i0 = paramfile_get_int_param_or_exit(r->ppl,"astat_lgn_proj_x0");
  j0 = paramfile_get_int_param_or_exit(r->ppl,"astat_lgn_proj_y0");
  i1 = paramfile_get_int_param_or_exit(r->ppl,"astat_lgn_proj_xn");
  j1 = paramfile_get_int_param_or_exit(r->ppl,"astat_lgn_proj_yn");
  i1 += i0 -1;
  j1 += j0 -1;

  // Temporary float storage to compute SD of x- and y-coords
  fcx = get_farray(maxsyn);
  fcy = get_farray(maxsyn);
  umx = pl->area->umx;
  umy = pl->area->umy;

  if (umx != umy)
    exit_error("MOD_CONN_ASTAT_LGN_POSTSYN","umx != umy");
  
  // Find maximum distance between any pair of connections
  for(i=i0;i<=i1;i++){
    for(j=j0;j<=j1;j++){ // For each LGN unit in window
      n = lgnn[i][j];
      dmax = -1.0;         // Max dist, in raw units
      lmax = kmax = 0;     // Indices in synapse list for dmax pair
      if (n > 1){
	tx = lgnx[i][j];
	ty = lgny[i][j];
	for(k=0;k<n;k++){    // For each synapse
	  for(l=k;l<n;l++){  // For each remaining synapse
	    dx = (float)(tx[k] - tx[l]);
	    dy = (float)(ty[k] - ty[l]);
	    if (dxf != 1.0)
	      dx *= dxf;           // Scaling factor, for asym. CMF
	    dd = dx*dx + dy*dy;
	    if (dd > dmax){
	      dmax = dd;
	      lmax = l;
	      kmax = k;
	    }
	  }
	}

	for(k=0;k<n;k++){ // Compute SDs of x and y coords
	  fcx[k] = (float)tx[k] * umx;
	  fcy[k] = (float)ty[k] * umy;
	}
	mean_sdev_farray(fcx,n,&xmu,&xsd);
	mean_sdev_farray(fcy,n,&ymu,&ysd);
	//printf("dxf,xsd = %f %f\n",dxf,xsd);
	if (dxf != 1.0)
	  xsd *= dxf;           // Scaling factor, for asym. CMF

	//printf("dxf,xsd = %f %f\n",dxf,xsd);
	//exit(0);


	// Append stats to outfile
	sprintf(tstr,
	"%3d %3d %s n %3d  p1 %3d %3d  p2 %3d %3d  d %8.2f  sd %8.2f %8.2f\n",
		i,j,nstr,n,tx[kmax],ty[kmax],tx[lmax],ty[lmax],
		sqrt(dmax)*umx,xsd,ysd);
      }else{
	;
	/*
	sprintf(tstr,
	"%3d %3d %s n %3d  p1 %3d %3d  p2 %3d %3d  d %8.2f  sd %8.2f %8.2f\n",
	i,j,nstr,n,0,0,0,0,-1.0,0.0,0.0);*/
      }
      if (n > 1)
	append_string_to_file(outfile,tstr);
    }
  }

  printf("    Statistics columns in outfile are:\n");
  printf("      1.  LGN x index (int)\n");
  printf("      2.  LGN y index (int)\n");
  printf("      3.  ON or OFF (char) and L or R\n");
  printf("      4.  Post-syn. layer (char)\n");
  printf("      5.  'n' (char)\n");
  printf("      6.  Number of connections (int)\n");
  printf("      7.  'p1' (char)\n");
  printf("      8.  x-index in post-syn. layer of p1 (int)\n");
  printf("      9.  y-index in post-syn. layer of p1 (int)\n");
  printf("     10.  'p2' (char)\n");
  printf("     11.  x-index in post-syn. layer of p2 (int)\n");
  printf("     12.  y-index in post-syn. layer of p2 (int)\n");
  printf("     13.  'd' (char)\n");
  printf("     14.  Max. connection dist. microns, p1 to p2 (float)\n");
  printf("     15.  'sd' (char)\n");
  printf("     16.  SD of x-positions, microns, of all synapses (float)\n");
  printf("     17.  SD of y-positions, microns, of all synapses (float)\n");
  
  myfree(fcx); myfree(fcy);
  free_3d_iarray(lgnx,xn,yn,maxsyn);
  free_3d_iarray(lgny,xn,yn,maxsyn);
  free_2d_farray(lgnn,xn);
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_CONN_ASTAT_LMS_INPUT                        */
/*                                                                           */
/*  Compute and write statistics for LMS distribution of inputs.             */
/*                                                                           */
/*****************************************************************************/
void mod_conn_astat_lms_input(mpt,r)
     struct pop_top *mpt;          // Global model structure
     struct response_struct *r;    // Response params
{
  int i;
  int ncell,lms_id,**st;
  char tstr[SLEN],*lname,outfile[SLEN];
  struct pop_syn *s,*t;
  struct pop_cell *c,*cpre;
  struct pop_layer *lay_pre,*lay_post;

  printf("  MOD_CONN_ASTAT_LMS_INPUT\n");

  // Get a pointer to the post-syn layer, specified in the .rsp file
  lname = paramfile_get_char_param_or_exit(r->ppl,"astat_lms_input_post");
  lay_post = pop_cell_get_layer_pointer(mpt->lay,mpt->nlay,lname);
  myfree(lname);

  lname = paramfile_get_char_param_or_exit(r->ppl,"astat_lms_input_pre");
  lay_pre = pop_cell_get_layer_pointer(mpt->lay,mpt->nlay,lname);
  myfree(lname);

  // Make outfile
  sprintf(outfile,"%s.astat",mpt->out_prefix);
  remove_file(outfile);
  printf("    Writing statistics to %s\n",outfile);

  printf("    Connections from %s to %s\n",lay_pre->name,lay_post->name);


  for(i=0;i<r->n;i++){ // For each response requested

    //  Get a pointer to the relevant cell
    if (r->rformat[i] == 0){
      ;
    }else if (r->rformat[i] == 10){
      ;
    }else{
      c = pop_cell_laylist_get_cell_pointer(mpt->lay,mpt->nlay,r->plname[i],
					    r->xi[i],r->yi[i],r->zi[i]);
      if (c == NULL){
	printf("*** MOD_CONN_ASTAT_LMS_INPUT  No cell found for %s\n",
	       r->nd_name[i]);
	exit_error("MOD_CONN_ASTAT_LMS_INPUT","Cell not found for response");
      }

      ncell = pop_cell_count_syn_from_layer_to_cell(lay_pre->name,c);

      st = get_zero_2d_iarray(2,3);  // ON/OFF vs CID 0,1,2

      // Scan pre-syn list for synapses from 'lay'
      s = pop_cell_get_next_syn_layer_name(c->in,1,lay_pre->name);
      while(s != NULL){
	cpre = s->pre;
	lms_id = lay_pre->irr_id[cpre->layx];
	//printf("  Name:  %s   z %d  lms_id:  %d\n",cpre->name,cpre->layz,
	//lms_id);

	st[cpre->layz][lms_id] += 1;

	t = s->post_next; // The post-syn cell of an input syn is this cell
	s = pop_cell_get_next_syn_layer_name(t,1,lay_pre->name);
      }

      sprintf(tstr,"%s  ON_SML %2d %2d %2d  OFF_SML %2d %2d %2d\n",c->name,
	      st[1][0],st[1][1],st[1][2],st[0][0],st[0][1],st[0][2]);
      append_string_to_file(outfile,tstr);
    }
  }

  printf("    Statistics columns in outfile are:\n");
  printf("      1.  Post-syn cell name\n");
  printf("      2.  'ON_SML'\n");
  printf("      3.  Number of ON S inputs (int)\n");
  printf("      4.  Number of ON M inputs (int)\n");
  printf("      5.  Number of ON L inputs (int)\n");
  printf("      6.  'OFF_SML'\n");
  printf("      7.  Number of OFF S inputs (int)\n");
  printf("      8.  Number of OFF M inputs (int)\n");
  printf("      9.  Number of OFF L inputs (int)\n");
  
}
/**************************************-**************************************/
/*                                                                           */
/*                            MAKE_RAD_COORD_ARRAYS                          */
/*                                                                           */
/*  Fill global arrays with coordinates of square rings of various radii.    */
/*                                                                           */
/*****************************************************************************/
void make_rad_coord_arrays()
{
  int i,j,k;

  radcx = (int **)myalloc((radmax+1)*sizeof(int *));
  radcy = (int **)myalloc((radmax+1)*sizeof(int *));

  radcx[0] = (int *)NULL; // unused
  radcy[0] = (int *)NULL;

  for(i=1;i<=radmax;i++){
    radcx[i] = (int *)myalloc((i*8)*sizeof(int));
    radcy[i] = (int *)myalloc((i*8)*sizeof(int));
    k = 0;
    for(j=0;j<(i*2);j++){  // right
      radcx[i][k] = -i;
      radcy[i][k] = 1-i + j;
      k += 1;
    }
    for(j=0;j<(i*2);j++){  // top
      radcx[i][k] = 1-i + j;
      radcy[i][k] = i;
      k += 1;
    }
    for(j=0;j<(i*2);j++){  // left
      radcx[i][k] = i;
      radcy[i][k] = i-1 - j;
      k += 1;
    }
    for(j=0;j<(i*2);j++){  // bottom
      radcx[i][k] = i-1 - j;
      radcy[i][k] = -i;
      k += 1;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            PICK_RAN_NEAREST_COORD                         */
/*                                                                           */
/*****************************************************************************/
void pick_ran_nearest_coord(p,xn,yn,x0,y0,signflag,pseed,rx,ry)
     float **p;
     int xn,yn;
     int x0,y0;
     int signflag;   // -1 or 1, tells whether to look for OFF or ON conn
     int *pseed;
     int *rx,*ry;
{
  int k,kk;
  int done,nrad,irad,xi,yi;
  
  if ((p[x0][y0]*signflag) > 0.0){
    *rx = x0; // Take the center position
    *ry = y0;
  }else{
    done = 0;
    irad = 1;
    while(!done && (irad <= radmax)){
      nrad = irad * 8; // total spots that exist around the square
      // start at 'k', chosen randomly around the square
      k = (int)((float)nrad * myrand_util_ran2(pseed));
      kk = 0; // how far through one cycle of the square
      while(!done && (kk < nrad)){
	xi = x0 + radcx[irad][k];
	yi = y0 + radcy[irad][k];
	if ((xi >= 0) && (xi < xn) && (yi >= 0) && (yi < yn)){
	  if ((p[xi][yi] * signflag) > 0.0){
	    *rx = xi;
	    *ry = yi;
	    done = 1;
	  }else{
	    ;  // Failed because prob 'p' may be zero
	    // printf("%4d %4d - failed  (%d %d)  p %f  (s %d)\n",
	    // xi,yi,x0,y0,p[xi][yi],signflag);
	  }
	}else{
	  ;  // Failed because point is outside [xn,yn] grid
	}
	if (!done){
	  kk += 1;
	  k += 1;
	  if (k == nrad) // Wrap around
	    k = 0;
	}
      }
      if (!done)
	irad += 1; // Go to the next radius
    }
    if (!done){
      *rx = *ry = -1;
      done = 1;
      exit_error("PICK_RAN_NEAREST_COORD","radmax exceeded");
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_CONN_EVEN_OUT                            */
/*                                                                           */
/*  NOTE:  This can cause duplicate sampling points, and does not observe    */
/*         'maxrep'.                                                         */
/*                                                                           */
/*****************************************************************************/
void mod_conn_even_out(mylogf,p,xn,yn,seed,cx,cy,n,maxit,sig,dumpflag)
     char *mylogf;
     float **p;         // [xn][yn] Prob of connection
     int xn,yn,seed;
     int *cx,*cy,n;
     int maxit;         // The maximum number of swaps, -1 for default

     float sig;         // The SD to use for smoothing
           // **** Currently, this is ALWAYS 2.0 ****

     int dumpflag;      // 1-print and write plots
{
  int i,j,k,kk;
  int tn,xi,yi,xlo,ylo,xhi,yhi,sign_lo,sign_hi,*pseed,xmv,ymv,itn;
  float **q,**d,**sm,pmu,psd,smu,ssd,*sqerr;

  //mylog(mylogf,"  MOD_CONN_EVEN_OUT\n");
  //printf("________________ MOD_CONN_EVEN_OUT    maxit = %d\n",maxit);
  //printf("________________ MOD_CONN_EVEN_OUT    sig = %f\n",sig);

  if (maxit <= 0)
    maxit = 50;

  if (radcx == NULL)
    make_rad_coord_arrays();

  pseed = &seed;
  if (seed > 0.0)
    seed *= -1.0;

  single_mean_sdev_2d_farray(p,xn,yn,&pmu,&psd);

  // Establish the connection matrix, +1 for ON, -1 for OFF
  q = get_zero_2d_farray(xn,yn);
  for(i=0;i<n;i++){
    xi = cx[i];
    yi = cy[i];
    if (p[xi][yi] >= 0.0)
      q[xi][yi] += 1.0;
    else
      q[xi][yi] -= 1.0;
  }

  if (dumpflag)
    sqerr = (float *)myalloc(maxit*sizeof(float)); // See how error converge
  else
    sqerr = NULL;
  
  d = (float **)myalloc(xn*sizeof(float *));
  for(i=0;i<xn;i++)
    d[i] = (float *)myalloc(yn*sizeof(float));
  
  itn = 0;
  while(itn < maxit){
    // Smooth the connection array
    sm = smooth_2d_with_gaussian(q,xn,yn,sig,0.05);

    // Match SD to that of prob function
    single_mean_sdev_2d_farray(sm,xn,yn,&smu,&ssd);
    //printf(" ______ before MU = %f  SD = %f\n",smu,ssd);
    multiply_2d_farray(sm,xn,yn,psd/ssd);
    // WYETH Consider whether we should match the mean also:
    //  single_mean_sdev_2d_farray(sm,xn,yn,&smu,&ssd);
    //   add_const_2d_farray(sm,0,xn,0,yn,pmu-smu);
    single_mean_sdev_2d_farray(sm,xn,yn,&smu,&ssd); // Extraneous ??
    //printf(" ________after MU = %f  SD = %f\n",smu,ssd);

    //printf("pmu psd  %f %f\n",pmu,psd);
    //exit_error("Done here","wyeth");

    if ((itn == 0) && dumpflag){
      write_2d_data("zev.p0.2d",p,0,0,xn,yn,4,2,1,0);
      write_2d_data("zev.sm0.2d",sm,0,0,xn,yn,4,2,1,0);
      write_2d_data("zev.q0.2d",q,0,0,xn,yn,4,2,1,0);
    }

    // Subtract 'p' from 'sm', to get *postive* errors
    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){
	if (p[i][j] >= 0.0)
	  d[i][j] = sm[i][j] - p[i][j];
	else
	  d[i][j] = p[i][j] - sm[i][j];
      }
    }

    if (dumpflag)
      sqerr[itn] = sum_square_2d_farray(d,xn,yn);

    min_max_coord_2d_farray(d,xn,yn,&xlo,&ylo,&xhi,&yhi);

    if (dumpflag)
      printf("max at:  %d %d   q(%d,%d) = %f\n",xhi,yhi,xhi,yhi,q[xhi][yhi]);

    /* WYETH REMOVE DEBUG
    if ((xhi == 0) && (yhi == 31)){

    single_mean_sdev_2d_farray(p,xn,yn,&smu,&ssd); // Extraneous ??
    printf(" ________FINAL 'p'  MU = %f  SD = %f\n",smu,ssd);

      write_2d_data("zev.pX.2d",p,0,0,xn,yn,4,2,1,0);
      write_2d_data("zev.smX.2d",sm,0,0,xn,yn,4,2,1,0);
      write_2d_data("zev.qX.2d",q,0,0,xn,yn,4,2,1,0);
      write_2d_data("zev.dX.2d",d,0,0,xn,yn,4,2,1,0);
    }
    */

    
    if (p[xhi][yhi] >= 0.0)  // The largest overshoot is positive
      sign_hi = 1;
    else
      sign_hi = -1;

    if (p[xlo][ylo] >= 0.0)
      sign_lo = 1;
    else
      sign_lo = -1;

    
    // Find coord to move
    pick_ran_nearest_coord(q,xn,yn,xhi,yhi,sign_hi,pseed,&xmv,&ymv);
    if (dumpflag){
      printf("  coord to move = %d %d\n",xmv,ymv);
      printf("   q(%d,%d) = %f\n",xmv,ymv,q[xmv][ymv]);
    }

    // Move it to min (WYETH for now, but perhaps look for empty spot )
    q[xlo][ylo] +=  (float)sign_lo; // add a sampling point
    q[xmv][ymv] += -(float)sign_hi; // subtract a sampling point

    itn += 1;

    if ((itn == maxit) && dumpflag){
      write_2d_data("zev.smf.2d",sm,0,0,xn,yn,4,2,1,0);
      write_2d_data("zev.qf.2d",q,0,0,xn,yn,4,2,1,0);
    }

    free_2d_farray(sm,xn);
  }

  if (dumpflag){
    append_farray_plot("zev.sqerr.pl","sqerr_v_itnum",sqerr,itn,1);
    myfree(sqerr);
  }

  // Update the connection matrices - simply overwrite old connections
  kk = 0;
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      tn = my_rint(fabs(q[i][j]));
      for(k=0;k<tn;k++){
	cx[kk] = i;
	cy[kk] = j;
	kk += 1;
      }
    }
  }
  if (kk != n)
    mylog_exit(mylogf,"MOD_CONN_EVEN_OUT  Connection number mismatch\n");

  free_2d_farray(q,xn);
  free_2d_farray(d,xn);
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_CONN_DIFFUSE_NEAR_N                         */
/*                                                                           */
/*  Return the sum of the distances to the nearest 'nn' points from the      */
/*  (x,y) coordinate.                                                        */
/*                                                                           */
/*****************************************************************************/
float mod_conn_diffuse_near_n(px,py,n,x,y,pi,nn)
     int *px,*py,n;   // point coordinates and number of points
     int x,y;         // reference coordinate
     int pi;          // index of point to exclude, or -1 to 
     int nn;          // Number of nearest neighbors to consider
{
  int i,k;
  float *dist,dx,dy,dsum;

  if (nn < 1)
    exit_error("MOD_CONN_DIFFUSE_NEAR_N","Value of 'nn' must be >= 1");
  else if (nn > (n-1))
    exit_error("MOD_CONN_DIFFUSE_NEAR_N","Value of 'nn' must be <= n-1");

  if ((pi < -1) || (pi >= n)){
    printf("  pi = %d\n",pi);
    exit_error("MOD_CONN_DIFFUSE_NEAR_N","Bad value for 'pi'");
  }

  //
  //  Create a list of the distances from (x,y) to all points except 'pi'
  //
  dist = (float *)myalloc((n-1)*sizeof(float));
  k = 0;
  for(i=0;i<n;i++){
    if (i != pi){
      dx = (float)(px[i] - x);
      dy = (float)(py[i] - y);
      //dist[k] = sqrt(dx*dx + dy*dy);
      //dist[k] = 1.0 / sqrt(dx*dx + dy*dy);
      dist[k] = 1.0 / (dx*dx + dy*dy);
      k += 1;
    }
  }

  //  Sort the list.
  sort_farray(dist,n-1);

  reverse_farray(dist,n-1);

  //  Sum the first 'nn' entries.
  dsum = sum_farray(dist,n-1,0,nn);

  myfree(dist);

  return dsum;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_CONN_DIFFUSE_WORK_1                         */
/*                                                                           */
/*  Given a list of points and one chosen point in that list, find the       */
/*  location next to that point that is furthest from the point's 'nn'       */
/*  nearest neighbors.  If that location is further than the location of     */
/*  the current point, return this location, otherwise return the original   */
/*  location.                                                                */
/*                                                                           */
/*****************************************************************************/
float mod_conn_diffuse_work_1(mylogf,xn,yn,pmap,apcx,apcy,aprad,px,py,pn,nn,pi)
     char *mylogf;
     int xn,yn;         // Size of grid (pix)
     int **pmap;        // [xn][yn] 1 where there are points, 0 elsewhere
     float apcx,apcy;   // Circular aperture center (pix)
     float aprad;       // Circular aperture radius (pix)
     int *px,*py,pn;    // Point coordinates and number of points
     int nn;            // Number of nearest neighbors to consider
     int pi;            // The point to consider moving
{
  int i,j,k;
  int x,y,mini,minj;
  float *dsum,dref,dmin;

  //mylog(mylogf,"  MOD_CONN_DIFFUSE_WORK_1\n");

  dsum = (float *)myalloc(9*sizeof(float));  // Will only ever use 8?

  x = px[pi];  // The coordinate of the point to consider moving
  y = py[pi];

  dref = mod_conn_diffuse_near_n(px,py,pn,x,y,pi,nn);
  dmin = dref;
  mini = minj = -1;

  k = 0;
  for(i=x-1;i<=(x+1);i++){
    for(j=y-1;j<=(y+1);j++){
      if ((i >= 0) && (i < xn) && (j >= 0) && (j < yn)){
	if (pmap[i][j] == 0){
	  dsum[k] = mod_conn_diffuse_near_n(px,py,pn,i,j,pi,nn);
	  if (dsum[k] < dmin){
	    dmin = dsum[k];
	    mini = i;
	    minj = j;
	  }
	  k += 1;
	}else{
	  ; //printf("%d %d ALREADY TAKEN or OFF LIMITS\n",i,j);
	}
      }
    }
  }

  //printf("dref = %f  dmin = %f\n",dref,dmin);
  if (mini != -1){
    px[pi] = mini;  // Move the point to the new coordinate
    py[pi] = minj;

    pmap[x][y] -= 1;  // Update the map
    pmap[mini][minj] += 1;

    //printf("  Point (%d,%d) moved to (%d,%d)\n",x,y,mini,minj);
  }

  //for(i=0;i<k;i++)
  //printf("%4d %f\n",i,dsum[i]);

  myfree(dsum);

  return dmin;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_CONN_EVEN_OUT_DIFFUSE                        */
/*                                                                           */
/*  Use an iterative diffsion algorithm.                                     */
/*                                                                           */
/*****************************************************************************/
void mod_conn_even_out_diffuse(mylogf,xn,yn,cx,cy,n,itmax,apcx,apcy,apdiam,
			       prob,pmin,
			       dumpflag)
     char *mylogf;
     int xn,yn;         // Size of underlying grid
     int *cx,*cy,n;     // Cooridinates and number of points
     int itmax;         // The maximum number of iterations
     float apcx,apcy;   // Aperture to limit diffusion, center
     float apdiam;      // Aperture to limit diffusion, diameter
     float **prob;      // Prob map.
     float pmin;        // Don't allow spots with prob < this value.
     int dumpflag;      // 1-print and write plots
{
  int i,j;
  int nn,x,y,**pmap,itcnt,done;
  float ***fmap,dx,dy,dtot,aprad;

  if (dumpflag == 1)
    mylog(mylogf,"  MOD_CONN_EVEN_OUT_DIFFUSE\n");

  if (itmax <= 0)
    itmax = 20;

  aprad = apdiam/2.0;

  //
  //  Build a map of the current points, which will be updated as points move
  //
  pmap = get_zero_2d_iarray(xn,yn);
  for(i=0;i<n;i++){
    x = cx[i];
    y = cy[i];
    if ((x < 0) || (x >= xn) || (y < 0) || (y >= yn))
      exit_error("MOD_CONN_EVEN_OUT_DIFFUSE","Bad point coordinate");
    pmap[x][y] += 1;
  }

  //  Imprint the aperture on the 'pmap'
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      if (prob[i][j] < pmin)
	pmap[i][j] = 2;
    }
  }

  nn = 6;  // Number of neighbors to consider

  if (dumpflag == 1)
    fmap = get_3d_farray(xn,yn,itmax);

  itcnt = 0;
  done = 0;
  while(done == 0){
    dtot = 0.0;
    for(i=0;i<n;i++){
      dtot += mod_conn_diffuse_work_1(mylogf,xn,yn,pmap,apcx,apcy,aprad,
				      cx,cy,n,nn,i);
    }

    if (dumpflag == 1){
      for(i=0;i<xn;i++){
	for(j=0;j<yn;j++){
	  fmap[i][j][itcnt] = (float)(pmap[i][j]);
	}
      }
    }

    if (dumpflag == 1)
      printf("dtot = %f\n",dtot);

    itcnt += 1;
    if (itcnt >= itmax)
      done = 1;
  }

  free_2d_iarray(pmap,xn);

  if (dumpflag == 1){
    write_3d_data_part("zzz.fmap.3d",fmap,0,xn,0,yn,0,itmax,4,2,1);
    free_3d_farray(fmap,xn,yn,itmax);
    printf("  Exit after dumping output. Set 'dumpflag' to 0 to bypass.\n");
    exit(0);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                        MOD_CONN_PICK_PAIRED_CONNECTIONS                   */
/*                                                                           */
/*  Pick connections that fall around a spatial offset from the given        */
/*  connections.  When no pair is found use -1 for what would have been      */
/*  the coordinates for the paired input.                                    */
/*                                                                           */
/*****************************************************************************/
void mod_conn_pick_paired_connections(cx,cy,cn,xn,yn,theta,dr,sdpara,sdorth,
				      seed,rcx,rcy)
     int *cx,*cy,cn;       // Grid points and number
     int xn,yn;            // Grid size
     float theta;          // Direction (degrees)
     float dr;             // Offset Distance
     float sdpara,sdorth;  // Noise SD
     int seed;             // Seed for noise
     int **rcx,**rcy;      // New connections, paired w/ orig [cn]
{
  int i;
  int *tcx,*tcy,xi,yi;
  float x0,y0,x,y;
  float a,b,c,d;

  // Now pick correlated paired inputs
  tcx = (int *)myalloc(cn*sizeof(int));
  tcy = (int *)myalloc(cn*sizeof(int));

  if (seed > 0)
    seed = -seed;

  a =  cos(M_PI/180.0*theta); b = sin(M_PI/180.0*theta); // rotation matrix
  c = -sin(M_PI/180.0*theta); d = cos(M_PI/180.0*theta);

  for(i=0;i<cn;i++){

    x0 = dr + sdpara * nr_util_gasdev(&seed); // Noise along length of vector
    y0 = sdorth * nr_util_gasdev(&seed);      // Noise orthogonal to vector

    x = a*x0 + c*y0;
    y = b*x0 + d*y0;

    xi = my_rint(x) + cx[i];
    yi = my_rint(y) + cy[i];

    if ((xi >= 0) && (xi < xn) && (yi >= 0) && (yi < yn)){
      tcx[i] = xi;
      tcy[i] = yi;
    }else{
      tcx[i] = -1;   // No pair is made, the original connection is
      tcy[i] = -1;   //   left unpaired, and -1 is the indicator
    }
  }

  *rcx = tcx;
  *rcy = tcy;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_CONN_PHASE_PAIR_MATCH                         */
/*                                                                           */
/*  Return the cell within the 'z' dimension that matches the phase          */
/*  criterion.                                                               */
/*                                                                           */
/*****************************************************************************/
struct pop_cell *mod_conn_phase_pair_match(mylogf,cref,dir,lay_match,ix,iy,
					   ph_shift,ph_dev)
     char *mylogf;
     struct pop_cell *cref;  // reference cell
     float dir;              // target direction preference (e.g. post-syn)
     struct pop_layer *lay_match; // layer in which to search for the match
     int ix,iy;              // coordinates within lay_match
     float ph_shift;         // target phase shift for match (deg)
     float ph_dev;           // plus or minus this amount around ph_shift (deg)
{
  int i;
  int done,ai_phase;
  float ph0,ph1,pha,phb,phd,t,ddir,ori_ref;
  struct pop_cell *c,*cptr;
  struct pop_layer *lay_ref;

  if (ph_dev < 0.0)
    mylog_exit(mylogf,"MOD_CONN_PHASE_PAIR_MATCH  'ph_dev' is negative\n");

  lay_ref = cref->pl;  // layer of ref cell

  ai_phase =  pop_cell_attrib_index_f(lay_ref,"phase"); // Attribute index

  // Compute difference between the directions of the pre- and post-syn cells
  ori_ref = pop_cell_attrib_get_f(cref,"ori");
  ddir = mod_conn_get_circ_dist(dir,ori_ref,360.0);

  ph0 = cref->attrib_f[ai_phase];
  pha = ph_shift - ph_dev;
  phb = ph_shift + ph_dev;
  if (pha < -180.0)  // Ideally, all phase values are in (-180,180]
    pha += 360.0;
  if (phb > 180.0)
    phb -= 360.0;

  if ((pha < -180.0) || (pha > 180.0) || (phb < -180.0) || (phb > 180.0)){
    printf("pha, phb  %f  %f\n",pha,phb);
    mylog_exit(mylogf,"MOD_CONN_PHASE_PAIR_MATCH  Phase range error.\n");
  }

  i = 0;  // Look at all cells
  done = 0;
  cptr = NULL;
  while(!done){
    c = &(lay_match->c[ix][iy][i]);  // pointer to cell to be tested for match

    ph1 = c->attrib_f[ai_phase];  // phase of candidate cell
    phd = mod_conn_get_signed_circ_dist(ph0,ph1,360.0);
    
    // This more complex conditional is designed to handle the case where
    //   the phase range spans the +/- 180 phase boundary.
    if (((pha <= phb) &&
	 (((ddir <= 90.0) && (phd >=  pha) && (phd <=  phb)) ||
	  ((ddir >  90.0) && (phd >= -phb) && (phd <= -pha)))) ||
	((pha > phb) && // the phase range spans the +/-180 boundary
	 (((ddir <= 90.0) && ((phd >=  pha) || (phd <=  phb))) ||
	  ((ddir >  90.0) && ((phd >= -phb) || (phd <= -pha)))))){
      // WYETH - ()'s around || part in 2 previous lines added July 2011

      //printf("NEW one\n");
      /*
	 (((ddir <= 90.0) && (phd >=  pha) || (phd <=  phb)) ||
	  ((ddir >  90.0) && (phd >= -phb) || (phd <= -pha))))){
      */
      
      //printf("%4f Ddir %f  PhDiff %f  pha,b %f %f\n",dir,ddir,phd,pha,phb);
      //printf("%2d %2d phd = %f   ph_cand:  %f  ph_ref %f  [%f   %f %f]\n",
      //i,k,phd,ph1,ph0,ph_shift,pha,phb);
      //printf("  will pair %d with %d\n",iz,k);
      
      cptr = c;
      done = 1;
    }else{
      i += 1;
      if (i >= lay_match->zn)
	done = 1;
    }
  }

  return cptr;
}
/**************************************-**************************************/
/*                                                                           */
/*                       MOD_CONN_CONVERT_SYN_TO_PHASE_PAIR                  */
/*                                                                           */
/*  Convert a regular connection to a paired connection based on the phase   */
/*  parameter of the target population.                                      */
/*                                                                           */
/*****************************************************************************/
int mod_conn_convert_syn_to_phase_pair(mylogf,syn,msi,inindex,dir)
     char *mylogf;
     struct pop_syn *syn;
     struct pop_mech *msi;    // Pointer to SI definition
     short inindex;           // Index into population 'inlist'
     float dir;               // Direction pref. of post-syn DS cell
{
  int ix,iy,circn,cnt;
  float circdt,ph_shift,ph_dev;
  struct pop_cell *cpre,*c;
  struct pop_layer *lpre;
  struct pop_syn *syn1;

  ph_shift = onode_getpar_flt_dflt(msi->o,"ph_shift",90.0);  // Phase shift
  ph_dev   = onode_getpar_flt_dflt(msi->o,"ph_dev",5.0);  // Deviate allowed

  cpre = syn->pre;  // pointer to pre-syn REFERENCE cell
  ix = cpre->layx;  // coords of ref cell
  iy = cpre->layy;  
  lpre = cpre->pl;  // layer of ref cell

  //  Look for one cell that matches along the 'z' dimension
  c = mod_conn_phase_pair_match(mylogf,cpre,dir,lpre,ix,iy,ph_shift,ph_dev);

  if (c == NULL)
    cnt = 0;
  else if (c == cpre){
    // WYETH - WE probably don't want to pair to ourselves???
    mylog_exit(mylogf,"MOD_CONN_CONVERT_SYN_TO_PHASE_PAIR  self-pair.\n");
  }else{

    circn  = msi->wtn;
    circdt = msi->wdt;

    syn1 = pop_cell_add_synapse(c,syn->post,1,1.0,0.0,inindex);
    pop_cell_si_make_001(syn,syn1,msi,circn,circdt);

    cnt = 1;
  }
    
  return cnt;
}
/**************************************-**************************************/
/*                                                                           */
/*                    MOD_CONN_CONVERT_SYN_TO_RFDIST_PAIR                    */
/*                                                                           */
/*  Convert a regular connection to a paired connection based on the         */
/*  distance of the RF centers.                                              */
/*                                                                           */
/*****************************************************************************/
int mod_conn_convert_syn_to_rfdist_pair(mylogf,syn,msi,inindex,pseed)
     char *mylogf;
     struct pop_syn *syn;
     struct pop_mech *msi;    // Pointer to SI definition
     short inindex;           // Index into population 'inlist'
     int *pseed;              // Pointer to randomization seed
{
  int k;
  int x1,y1,x2,y2,cnt;
  float odir,odist,jitp,jito,dr,sdpara,sdorth,dx;
  struct pop_cell *cpre,*c;
  struct pop_layer *lpre;
  struct pop_syn *syn1;

  odir  = onode_getpar_flt_exit(msi->o,"offset_direction");   // Vector dir.
  odist = onode_getpar_flt_exit(msi->o,"offset_distance");    // Vector length
  jitp  = onode_getpar_flt_exit(msi->o,"offset_jitter_para"); // parallel jit.
  jito  = onode_getpar_flt_exit(msi->o,"offset_jitter_orth"); // orthog jitter
  // Note, "offset_jitter_seed" was picked up below, once for entire layer

  cpre = syn->pre;  // pointer to pre-syn REFERENCE cell
  x1 = cpre->layx;  // coords of ref cell
  y1 = cpre->layy;
  k = cpre->layz;
  lpre = cpre->pl;  // layer of ref cell

  /*  WYETH - THIS CONDITION NOT NEEDED ??
  if (lpre->zn > 1)
    mylog_exit(mylogf,"MOD_CONN_CONVERT_TO_RFDIST_PAIR  z-layer >0 ignored\n");
  */

  if (lpre->area == NULL)
    mylog_exit(mylogf,"MOD_CONN_CONVERT_TO_RFDIST_PAIR  Area is null\n");

  //printf("sscale = %f\n",lpre->area->sscale);

  // WYETH - assumes that area has square sscale
  dx = lpre->xf * lpre->area->sscale;  // Degr per unit in this area
  dr = odist / dx;
  sdpara = jitp / dx;
  sdorth = jito / dx;

  // 'pseed' is used to preserve random sequence all calls for layer
  mod_conn_pick_paired_coord(x1,y1,odir,dr,sdpara,sdorth,pseed,&x2,&y2);
  
  cnt = 0;
  //if (x2 != -1){  // WYETH BUG WRONG
  if ((x2 >= 0) && (x2 < lpre->xn) && (y2 >= 0) && (y2 < lpre->yn)){

    c = &(lpre->c[x2][y2][k]);

    // Use weight from 'syn' for new syn
    syn1 = pop_cell_add_synapse(c,syn->post,1,syn->w,0.0,inindex);

    //
    // WYETH - this code seems redundant with 'mod_conn_si_pair_make'
    // could it be called instead?
    //
    if (onode_test_chr(msi->o,"nonlin","mask"))
      pop_cell_si_make_001(syn,syn1,msi,msi->wtn,msi->wdt);
    else if (onode_test_chr(msi->o,"nonlin","symmask"))
      pop_cell_si_make_002(syn,syn1,msi,msi->wtn,msi->wdt,-1.0,0);
    else if (onode_test_chr(msi->o,"nonlin","mult")||
	     onode_test_chr(msi->o,"nonlin","divide")){
      // WYETH-NEW-MULT
      pop_cell_si_make_002(syn,syn1,msi,msi->wtn,msi->wdt,0.0,1);
      syn->post->syn_proc_code = 1;
      syn->si->syn_code = 100002;
      syn1->si->syn_code = 100002;
    }else
      mylogx(mylogf,"MOD_CONN_CONVERT_TO_RFDIST_PAIR","Unknown 'nonlin'");

    cnt += 1;
  }

  return cnt;
}
/**************************************-**************************************/
/*                                                                           */
/*                      MOD_CONN_INPUT_PAIR_MATCH_RFDIST                     */
/*                                                                           */
/*  Return a cell to be paired with 'c' under the conditions of 'ip' using   */
/*  RF distance as a guide.                                                  */
/*                                                                           */
/*****************************************************************************/
void mod_conn_pair_rfdist_get_xy(lf,c,ip,rx,ry)
     char *lf;               // log file name
     struct pop_cell *c;     // Post-syn cell, "reference cell"
     struct input_pair *ip;  // Input pair structure
     int *rx,*ry;            // Return x,y coords for match
{
  int x1,y1;
  float odir,odist,jitp,jito,dx,dr,sdpara,sdorth;
  struct pop_layer *lay;

  odir  = ip->rf_offset_dir;
  odist = ip->rf_offset_dist;
  jitp  = ip->rf_offset_sd_par;
  jito  = ip->rf_offset_sd_orth;

  lay = ip->pop1;   // Seek cell in this layer

  //
  //  Determine (x2,y2), the x,y coords of the "matching cell"
  //
  x1 = c->layx;  // coords of ref cell
  y1 = c->layy;

  //  Be sure there is an 'area' and, FOR NOW, MAKE RESTRICTIONS ON GEOMETRY
  if (lay->area == NULL)
    mylogx(lf,"MOD_CONN_PAIR_RFDIST_GET_XY","Area is null");
  if (lay->xf != lay->yf)  // WYETH - should implement this
    mylogx(lf,"MOD_CONN_PAIR_RFDIST_GET_XY","Layer xf != yf");
  if (lay->area->xf != lay->area->yf)  // WYETH - should implement this
    mylogx(lf,"MOD_CONN_PAIR_RFDIST_GET_XY","Area xf != yf");

  dx = lay->xf * lay->area->sscale;  // Degr per unit in this area
  dr = odist / dx;
  sdpara = jitp / dx;
  sdorth = jito / dx;

  // 'pseed' is used to preserve random sequence all calls for layer
  mod_conn_pick_paired_coord(x1,y1,odir,dr,sdpara,sdorth,&(ip->tseed),rx,ry);

}
/**************************************-**************************************/
/*                                                                           */
/*                      MOD_CONN_INPUT_PAIR_MATCH_RFDIST                     */
/*                                                                           */
/*  Return a cell to be paired with 'c' under the conditions of 'ip' using   */
/*  RF distance as a guide.                                                  */
/*                                                                           */
/*****************************************************************************/
struct pop_cell *mod_conn_input_pair_match_rfdist(lf,c,ip)
     char *lf;               // log file name
     struct pop_cell *c;     // Pre-syn reference cell
     struct input_pair *ip;  // Input pair structure
{
  int x2,y2,z2;
  struct pop_cell *c2;
  struct pop_layer *lay;

  //  Get the x,y coords for the match
  mod_conn_pair_rfdist_get_xy(lf,c,ip,&x2,&y2);

  lay = ip->pop1;   // Seek cell in this layer

  //
  //  Determine 'z2', the z coords of the "matching cell"
  //  The following routine checks for out of bounds conditions.
  //
  z2 = mod_conn_zstring_get_unique_z(lay->zn,ip->z1,c->layz);

  //printf("z2 = %d  Def: %d   ==>%s<==\n",z2,c->layz,ip->z1);


  if ((x2 >= 0) && (x2 < lay->xn) && (y2 >= 0) && (y2 < lay->yn)){
    c2 = &(lay->c[x2][y2][z2]);
  }else
    c2 = NULL;

  return c2;
}
/**************************************-**************************************/
/*                                                                           */
/*                   MOD_CONN_INPUT_PAIR_MATCH_RFDIST_PHASE                  */
/*                                                                           */
/*  Return a cell to be paired with 'c' under the conditions of 'ip' using   */
/*  RF distance AND spatial phase as a guide.                                */
/*                                                                           */
/*****************************************************************************/
struct pop_cell *mod_conn_input_pair_match_rfdist_phase(lf,c,ip,c_post)
     char *lf;                 // log file name
     struct pop_cell *c;       // Pre-syn reference cell
     struct input_pair *ip;    // Input pair structure
     struct pop_cell *c_post;  // Post-syn cell, receives paired input
{
  int x2,y2,z2,ix,iy;
  float dir,ph_shift,ph_dev;
  struct pop_cell *c_match;
  struct pop_layer *lay; //,*lpre;

  //printf("    MOD_CONN_INPUT_PAIR_MATCH_RFDIST_PHASE  version 2\n");

  lay = ip->pop1;   // Seek cell in this layer

  //  Get the x,y coords for the match
  mod_conn_pair_rfdist_get_xy(lf,c,ip,&x2,&y2);

  if ((x2 >= 0) && (x2 < lay->xn) && (y2 >= 0) && (y2 < lay->yn)){

    //  The z-string should be 'all' for now; other options not imp'd
    z2 = mod_conn_zstring_get_unique_z(lay->zn,ip->z1,c->layz);
    if (z2 != -1)  // Note, -1 indicates 'all'
      mylogx(lf,"MOD_CONN_INPUT_PAIR_MATCH_RFDIST_PHASE","z origin not imp'd");

    //
    //  Now, find one z-match based on the appropriate phase criterion
    //

    // direction of post-syn cell is needed to condition on phase
    dir = pop_cell_attrib_get_f(c_post,"ori");
    
    ph_shift  = ip->rf_ph_shift;
    ph_dev    = ip->rf_ph_dev;
    ix = c->layx;  // coords of ref cell
    iy = c->layy;  
    //lpre = c->pl;  // layer of ref cell
    
    //  Look for one cell that matches along the 'z' dimension
    //c_match = mod_conn_phase_pair_match(lf,c,dir,lpre,ix,iy,ph_shift,ph_dev);
    c_match = mod_conn_phase_pair_match(lf,c,dir,lay,ix,iy,ph_shift,ph_dev);
    // Note, c_match may be NULL
  }else{
    c_match = NULL;
  }

  return c_match;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_CONN_SYN_CONVERT_CELL                        */
/*                                                                           */
/*  Convert all connections of 'inindex' from 'lay' to 'c' based on 'msi'    */
/*                                                                           */
/*****************************************************************************/
void mod_conn_syn_convert_cell(mylogf,lay,c,msi,inindex,pseed)
     char *mylogf;
     struct pop_layer *lay;  // pre-syn layer
     struct pop_cell *c;     // Post-syn cell
     struct pop_mech *msi;   // SI definition, e.g. for ds02
     short inindex;          // Index into population 'inlist'
     int *pseed;             // Pointer to randomization seed
{
  int n,cntz,cntm,nnew;
  float ori;
  //char tstr[SLEN];
  struct pop_syn *syn;

  n = 0;
  cntz = cntm = 0;

  syn = pop_cell_get_next_syn_index(c->in,1,inindex,lay->name); // inputs

  while(syn != NULL){
    // The number of new synapses is returned.

    if (strcmp(msi->type,"si_ds02")==0){
      ori = pop_cell_attrib_get_f(c,"ori");
      nnew = mod_conn_convert_syn_to_phase_pair(mylogf,syn,msi,inindex,ori);
    }else if (strcmp(msi->type,"si_pair_rfdist")==0){
      nnew = mod_conn_convert_syn_to_rfdist_pair(mylogf,syn,msi,inindex,pseed);
    }else{
      mylog_cval(mylogf,"*** mech type:",msi->type);
      mylogx(mylogf,"MOD_CONN_SYN_CONVERT_CELL","Unknown mechanism type");
    }

    if (nnew == 0)
      cntz += 1;  // How many had zero matches
    else if (nnew > 1)
      cntm += 1;  // How many had multiple matches

    syn = syn->post_next;
    //syn = pop_cell_get_next_syn_layer_name(syn,1,lay->name);
    syn = pop_cell_get_next_syn_index(syn,1,inindex,lay->name);
    n += 1;
  }
  //
  //  WYETH - should there be a report on the stats for these ??
  //
  /*
  sprintf(tstr,"      %d synapses,",n);
  mylog(mylogf,tstr);
  sprintf(tstr," %d unique match, %d no match, %d multiple match.\n",
	  n-(cntz+cntm),cntz,cntm);
  mylog(mylogf,tstr);
  */
}
/**************************************-**************************************/
/*                                                                           */
/*                       MOD_CONN_INPUT_PAIR_FIND_MATCH                      */
/*                                                                           */
/*  Find a match for 'c' given the constraints in 'ip'.                      */
/*                                                                           */
/*****************************************************************************/
struct pop_cell *mod_conn_input_pair_find_match(lf,c,ip,c_post)
     char *lf;                 // log file name
     struct pop_cell *c;       // "reference" cell (pre-syn)
     struct input_pair *ip;    // Input pair structure
     struct pop_cell *c_post;  // post-syn cell, which will receive pair
{
  struct pop_cell *c2;

  c2 = NULL;

  if (strcmp(ip->distrib,"rf_dist")==0){

    c2 = mod_conn_input_pair_match_rfdist(lf,c,ip);

  }else if (strcmp(ip->distrib,"rf_dist_phase")==0){

    c2 = mod_conn_input_pair_match_rfdist_phase(lf,c,ip,c_post);

  }else if (strcmp(ip->distrib,"phase_diff")==0){

    exit_error("MOD_CONN_INPUT_PAIR_FIND_MATCH","------ WYETH HERE -------");
    //nnew = mod_conn_convert_syn_to_rfdist_pair(mylogf,syn,msi,inindex,pseed);

  }else{
    mylog_cval(lf,"*** <pair_distrib> type:",ip->distrib);
    mylogx(lf,"MOD_CONN_INPUT_PAIR_FIND_MATCH","Unknown 'distrib' type");
  }

  return c2;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_CONN_SI_PAIR_MAKE                          */
/*                                                                           */
/*  Convert all connections of 'inindex' from 'lay' to 'c' based on 'msi'    */
/*                                                                           */
/*****************************************************************************/
void mod_conn_si_pair_make(lf,ip,syn0,syn1)
     char *lf;               // log file name
     struct input_pair *ip;  // Input pair structure
     struct pop_syn *syn0;   // Synapses to be paired
     struct pop_syn *syn1;
{
  struct pop_mech *msi;   // SI definition, e.g. for ds02

  msi = ip->msi;
  
  if (onode_test_chr(msi->o,"nonlin","mask"))
    pop_cell_si_make_001(syn0,syn1,msi,msi->wtn,msi->wdt);
  else if (onode_test_chr(msi->o,"nonlin","symmask"))
    pop_cell_si_make_002(syn0,syn1,msi,msi->wtn,msi->wdt,-1.0,0);
  else if (onode_test_chr(msi->o,"nonlin","mult")||
	   onode_test_chr(msi->o,"nonlin","divide")){
    // WYETH-NEW-MULT, currently, this is the same as 'symmask'
    pop_cell_si_make_002(syn0,syn1,msi,msi->wtn,msi->wdt,0.0,1);
    syn0->post->syn_proc_code = 1;  // Special processing for post-syn cell
    syn0->si->syn_code = 100002;
    syn1->si->syn_code = 100002;
  }else
    mylogx(lf,"MOD_CONN_SI_PAIR_MAKE","Unknown 'nonlin'");

}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_CONN_INPUT_PAIR_CELL                         */
/*                                                                           */
/*  Convert all connections of 'inindex' from 'lay' to 'c' based on 'msi'    */
/*                                                                           */
/*****************************************************************************/
void mod_conn_input_pair_cell(lf,c,io,ip)
     char *lf;               // log file name
     struct pop_cell *c;     // Post-syn cell
     struct onode *io;       // Input onode
     struct input_pair *ip;  // Input pair structure
{
  int n,cnt;
  char tstr[SLEN];
  float w;
  struct pop_syn *syn,*syn1;
  struct pop_cell *c1,*c2;

  n = 0;
  cnt = 0;
  syn = pop_cell_get_next_syn_index(c->in,1,ip->in_index,ip->pop0->name);
  while(syn != NULL){

    c1 = syn->pre;  // Reference cell

    //
    //  Get cell to pair w/ the reference cell
    //
    c2 = mod_conn_input_pair_find_match(lf,c1,ip,c);

    //
    // Create new paired synapse
    //
    if (c2 != NULL){

      //w = 1.0;   // WYETH - HOW TO SET SYN WEIGHT?
      w = syn->w;  // WYETH - Based on the old "rfdist"

      syn1 = pop_cell_add_synapse(c2,syn->post,1,w,0.0,ip->in_index);
      if (ip->msi != NULL){
	mod_conn_si_pair_make(lf,ip,syn,syn1);
      }
      cnt += 1;
    }

    n += 1;
    syn = syn->post_next;
    syn = pop_cell_get_next_syn_index(syn,1,ip->in_index,ip->pop0->name);
  }
  sprintf(tstr,"      %d pairs created from %d original synapses\n",n,cnt);
  mylog(lf,tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_CONN_SYN_CONVERT_LAYER                        */
/*                                                                           */
/*  Convert all connections of 'inindex' from 'lay1' to 'lay2' based on      */
/*  the mech 'msi'.                                                          */
/*                                                                           */
/*****************************************************************************/
void mod_conn_syn_convert_layer(mylogf,lay1,lay2,msi,inindex)
     char *mylogf;
     struct pop_layer *lay1;  // pre-syn layer
     struct pop_layer *lay2;  // post-syn layer
     struct pop_mech *msi;    // SI definition, e.g., for ds02
     short inindex;           // Index into population 'inlist'
{
  //
  //
  //
  //   WYETH - OLD WAY
  //
  //
  int i,j,k;
  int xn,yn,zn,seed;
  struct pop_cell *c;

  mylog(mylogf,"    MOD_CONN_CONVERT_LAYER_TO_DS02\n");

  xn = lay2->xn;
  yn = lay2->yn;
  zn = lay2->zn;

  // For any seed, must get the seed here, and keep a pointer to it so
  // that we have one long random sequence for all cells.
  if (strcmp(msi->type,"si_pair_rfdist")==0){
    seed = onode_getpar_int_exit(msi->o,"offset_jitter_seed"); // jitter seed
    if (seed > 0)
      seed = -seed;
  }else
    seed = -1;  // No seed needed

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	c = &(lay2->c[i][j][k]);
	mod_conn_syn_convert_cell(mylogf,lay1,c,msi,inindex,&seed);
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_CONN_INPUT_PAIR_MAIN                         */
/*                                                                           */
/*  Pair all connections of 'ip->in_index' for the post-synaptic layer.      */
/*  The synaptic interaction 'ip->msi' will be used if not NULL;             */
/*                                                                           */
/*****************************************************************************/
void mod_conn_input_pair_main(mpt,lay_post,io,ip)
     struct pop_top *mpt;         // Global model structure
     struct pop_layer *lay_post;  // Post-synaptic layer
     struct onode *io;            // Input onode
     struct input_pair *ip;       // Input pair structure
{
  int i,j,k;
  int xn,yn,zn;
  struct pop_cell *c;

  mylog(mpt->logf,"    MOD_CONN_INPUT_PAIR_MAIN\n");

  xn = lay_post->xn;
  yn = lay_post->yn;
  zn = lay_post->zn;
  
  if (ip->pair_seed > 0)  // Original seed value, not to be changed
    ip->tseed = -ip->pair_seed;  // Copy, to be used and changed

  // For each cell in "lay_post", pair any appropriate existing synapses.
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	c = &(lay_post->c[i][j][k]);
	mod_conn_input_pair_cell(mpt->logf,c,io,ip);
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_CONN_METH_01                             */
/*                                                                           */
/*  Scan all grid points and make a connection with probability at that      */
/*  point.                                                                   */
/*                                                                           */
/*  - The number of connections can vary from 0 to xn*yn.                    */
/*  - No connection can occur more than once.                                */
/*                                                                           */
/*****************************************************************************/
void mod_conn_meth_01(prob,xn,yn,seed,eps,rcx,rcy,rcn)
     float **prob;          // Connection prob matrix [xn][yn]
     int xn,yn;
     int seed;
     float eps;             // Make no connection if pr < eps
     int **rcx,**rcy,*rcn;  // Return x,y coords, and number of connections
{
  int i,j;
  int *cx,*cy,cn,nn;
  float p;

  if (seed > 0)
    seed *= -1;

  nn = xn*yn;
  cx = (int *)myalloc(nn*sizeof(int));
  cy = (int *)myalloc(nn*sizeof(int));

  cn = 0;
  for(i=0;i<xn;i++){
    for(j=0;i<yn;j++){
      p = prob[i][j];
      if (p > eps)
	if (myrand_util_ran2(&seed) < p){
	  cx[cn] = i;
	  cy[cn] = j;
	  cn += 1;
	}
    }
  }
  *rcx = copy_iarray(cx,cn);
  *rcy = copy_iarray(cy,cn);
  *rcn = cn;

  myfree(cx);
  myfree(cy);
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_CONN_METH_02                             */
/*                                                                           */
/*  For all points with prob > eps (eps >= 0), choose 'n' points such that   */
/*  no point is chosen more than 'maxrep' times.                             */
/*                                                                           */
/*  - The number of connections is fixed at n.                               */
/*                                                                           */
/*  Return the number of points that have prob > eps.                        */
/*                                                                           */
/*****************************************************************************/
int mod_conn_meth_02(mylogf,prob,xn,yn,tflag,n,seed,eps,maxrep,rcx,rcy)
     char *mylogf;
     float **prob;     // Connection prob matrix [xn][yn]
     int xn,yn;
     int **tflag;      // to allow positions to be chosen
     int n;            // Choose this many connections
     int seed;
     float eps;        // Make no connection if pr < eps
     int maxrep;       // 1 or more
     int **rcx,**rcy;  // Return x,y coords in [0..xn-1] and [0..yn-1]
{
  int i,j,k;
  int *tx,*ty,tn,*cx,*cy,nn,done,*flag;
  double *cdf,pmax,f;
  char tstr[SLEN];

  //printf("xn = %d  yn = %d\n",xn,yn);


  if (eps  < 0.0) mylog_exit(mylogf,"MOD_CONN_METH_02  Bad value for eps\n");
  if (maxrep < 1) mylog_exit(mylogf,"MOD_CONN_METH_02  Bad value of maxrep\n");
  if (n < 1)      mylog_exit(mylogf,"MOD_CONN_METH_02  Bad value for n\n");

  if (seed > 0)
    seed *= -1;

  nn = xn*yn;
  tx = (int *)myalloc(nn*sizeof(int));
  ty = (int *)myalloc(nn*sizeof(int));

  // Find all points with pr > eps
  // 'tx' and 'ty' will contain coords from [0..xn-1] and [0..yn-1]
  tn = 0;
  if (tflag == NULL){
    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){
	if (prob[i][j] > eps){  // Only check probability
	  tx[tn] = i;
	  ty[tn] = j;
	  tn += 1;
	}
      }
    }
  }else{
    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){
	if ((prob[i][j] > eps)&&(tflag[i][j])){ // Prob and 'tflag'
	  tx[tn] = i;
	  ty[tn] = j;
	  tn += 1;
	}
      }
    }
  }

  //if (tn < n){
  if ((tn*maxrep) < n){
    sprintf(tstr,"eps= %f  tn= %d maxrep= %d  n= %d\n",eps,tn,maxrep,n);
    mylog(mylogf,tstr);
    
    mylog_exit(mylogf,"MOD_CONN_METH_02  Too few points have pr > eps.\n");
  }

  // Make a float array of cumulative probabilities
  cdf = (double *)myalloc(tn*sizeof(double));
  for(k=0;k<tn;k++){
    i = tx[k];
    j = ty[k];
    if (k==0)
      cdf[k] = (double)prob[i][j];
    else
      cdf[k] = cdf[k-1] + (double)prob[i][j];
  }
  pmax = cdf[tn-1];

  // 'flag' is used to prevent picking same point too many times
  flag = (int *)myalloc(tn*sizeof(int));
  for(i=0;i<tn;i++)
    flag[i] = 0;

  // Choose 'n' connections
  cx = (int *)myalloc(n*sizeof(int));
  cy = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++){
    done = 0;
    while(!done){
      f = pmax * (double)myrand_util_ran2(&seed);
      k = 0;
      while((cdf[k] < f) && (k < (tn-1)))
	k += 1;
      if (k >= tn)
	mylog_exit(mylogf,"MOD_CONN_METH_02  This should not occur.\n");

      if (flag[k] < maxrep){
	cx[i] = tx[k];
	cy[i] = ty[k];
	flag[k] += 1;
	done = 1;
      }
    }
  }
  *rcx = cx;
  *rcy = cy;

  myfree(tx); myfree(ty); myfree(cdf); myfree(flag);

  return tn;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_CONN_GABOR_01                            */
/*                                                                           */
/*  Get connections according to a gabor function.                           */
/*                                                                           */
/*  - Max value of gabor is set to 1                                         */
/*  - No connections will be made where value is <= eps                      */
/*                                                                           */
/*****************************************************************************/
void mod_conn_gabor_01(xn,yn,xc,yc,sdo,sdp,sf,ph,ori,seed,maxrep,eps,n,rcx,rcy)
     int xn,yn;     // get connections to points on grid 'xn' by 'yn'
     float xc,yc;   // center of Gabor on the grid
     float sdo,sdp; // SD of Gaussian (pix) orth and par to ori
     float sf;      // SF of sinusoid (cyc/pix)
     float ph;      // phase of sinusoid (degr)
     float ori;     // orientation
     int seed;      // randomization
     int maxrep;    // maximum times any one point can be chosen (>=1)
     float eps;     // no connections to point w/ lower prob <= eps
     int n;         // number of connections to make
     int **rcx,**rcy;  // Return connection coordinates
{
  float **p;
  int nposs;

  if (eps < 0.0)
    exit_error("MOD_CONN_GABOR_O1","eps must be >=0");

  /*printf("xn,yn,xc,yc,sd,sf,ori,ph = ");
    printf("%d %d %f %f %f %f %f %f\n",xn,yn,xc,yc,sd,sf,ori,ph);*/
  
  p = gabor_2d_space_raw(xn,yn,xc,yc,sdo,sdp,sf,ori,ph);
  make_max_const_2d_farray(p,xn,yn,1.0); // so 'eps' is relative to max

  nposs = mod_conn_meth_02(NULL,p,xn,yn,NULL,n,seed,eps,maxrep,rcx,rcy);
  
  free_2d_farray(p,xn);
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_CONN_GABOR_02                            */
/*                                                                           */
/*  Get connections according to a gabor function.                           */
/*                                                                           */
/*  - Max value of gabor is set to 1                                         */
/*  - No connections will be made where value is <= eps                      */
/*                                                                           */
/*****************************************************************************/
int mod_conn_gabor_02(mylogf,xn,yn,tflag,xc,yc,sdo,sdp,sf,ph,ori,seed,maxrep,
		      eps,n,evflag,rcx1,rcy1,rcn1,rcx0,rcy0,rcn0)
     char *mylogf;  // For cluster
     int xn,yn;     // get connections to points on grid 'xn' by 'yn'
     int **tflag;   // flag for cells position allowed to be picked
     float xc,yc;   // center of Gabor on the grid
     float sdo,sdp; // SD of Gaussian (pix) orth and par to ori
     float sf;      // SF of sinusoid (cyc/pix)
     float ph;      // phase of sinusoid (degr)
     float ori;     // orientation
     int seed;      // randomization
     int maxrep;    // maximum times any one point can be chosen (>=1)
     float eps;     // no connections to point w/ lower prob <= eps
     int n;         // number of connections to make
     int evflag;    // 1-even out connections
     int **rcx1,**rcy1;  // Return ON connection coords [0..xn-1],[0..yn-1]
     int *rcn1;          // Number of ON
     int **rcx0,**rcy0;  // Return OFF connection coordinates
     int *rcn0;          // Number of OFF
{
  int i;
  int *cx,*cy,*cx0,*cy0,*cx1,*cy1,cn0,cn1,k0,k1,nposs;
  float **p,**ppos;

  if (eps < 0.0)
    mylog_exit(mylogf,"MOD_CONN_GABOR_O2  eps must be >=0\n");

  /*printf("xn,yn,xc,yc,sd,sf,ori,ph = ");
    printf("%d %d %f %f %f %f %f %f\n",xn,yn,xc,yc,sd,sf,ori,ph);*/

  p = gabor_2d_space_raw(xn,yn,xc,yc,sdo,sdp,sf,ori,ph);
  ppos = copy_2d_farray(p,xn,yn);
  absolute_value_2d_farray(ppos,xn,yn);

  // Make the largest value be 1, so 'eps' is relative to max=1.0
  make_max_const_2d_farray(ppos,xn,yn,1.0);

  // get connections, return number of possible contacts
  nposs = mod_conn_meth_02(mylogf,ppos,xn,yn,tflag,n,seed,eps,maxrep,&cx,&cy);

  if (evflag){ // Even out the random connections
    mod_conn_even_out(mylogf,p,xn,yn,seed,cx,cy,n,50,2.0,0);
  }

  // Count how many of each type, then extract into separate arrays
  cn0 = cn1 = 0;
  for(i=0;i<n;i++){
    if (p[cx[i]][cy[i]] > 0.0)
      cn1 += 1;
    else
      cn0 += 1;
  }
  cx0 = (int *)myalloc(cn0*sizeof(int));
  cy0 = (int *)myalloc(cn0*sizeof(int));
  cx1 = (int *)myalloc(cn1*sizeof(int));
  cy1 = (int *)myalloc(cn1*sizeof(int));

  k0 = k1 = 0;
  for(i=0;i<n;i++){
    if (p[cx[i]][cy[i]] > 0.0){
      cx1[k1] = cx[i];
      cy1[k1] = cy[i];
      k1 += 1;
    }else{
      cx0[k0] = cx[i];
      cy0[k0] = cy[i];
      k0 += 1;
    }
  }
  *rcx1 = cx1; *rcy1 = cy1; *rcn1 = cn1;
  *rcx0 = cx0; *rcy0 = cy0; *rcn0 = cn0;

  myfree(cx); myfree(cy);
  free_2d_farray(ppos,xn);
  free_2d_farray(p,xn);

  return nposs;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_CONN_GABOR_02_IRR                          */
/*                                                                           */
/*  Get connections according to a gabor function.                           */
/*                                                                           */
/*  - Max value of gabor is set to 1                                         */
/*  - No connections will be made where value is <= eps                      */
/*                                                                           */
/*****************************************************************************/
int mod_conn_gabor_02_irr(mylogf,lay,xn,yn,tflag,xc,yc,sdo,sdp,sf,ph,ori,seed,
			  maxrep,eps,n,evflag,cond,
			  rcx1,rcy1,rcn1,rcx0,rcy0,rcn0)
     char *mylogf;  // For cluster
     struct pop_layer *lay;    // pre-syn layer, has ->irr_... values
     int xn,yn;     // get connections to points on grid 'xn' by 'yn'
     int **tflag;   // flag for cells position allowed to be picked
     float xc,yc;   // center of Gabor on the grid
     float sdo,sdp; // SD of Gaussian (pix) orth and par to ori
     float sf;      // SF of sinusoid (cyc/pix)
     float ph;      // phase of sinusoid (degr)
     float ori;     // orientation
     int seed;      // randomization
     int maxrep;    // maximum times any one point can be chosen (>=1)
     float eps;     // no connections to point w/ lower prob <= eps
     int n;         // number of connections to make
     int evflag;    // 1-even out connections
     char *cond;    // Condition on cell type, "0","1",...
     int **rcx1,**rcy1;  // Return ON connection coords [0..xn-1],[0..yn-1]
     int *rcn1;          // Number of ON
     int **rcx0,**rcy0;  // Return OFF connection coordinates
     int *rcn0;          // Number of OFF
{
  int i;
  int *cx,*cy,*cx0,*cy0,*cx1,*cy1,cn0,cn1,k0,k1,nposs,cnt;
  float *p1,**p,**ppos;

  if (eps < 0.0)
    mylog_exit(mylogf,"MOD_CONN_GABOR_O2_IRR  eps must be >=0\n");

  /*printf("xn,yn,xc,yc,sd,sf,ori,ph = ");
    printf("%d %d %f %f %f %f %f %f\n",xn,yn,xc,yc,sd,sf,ori,ph);*/

  cnt = lay->irr_n;  // Number of units in pre-syn map

  //
  //  Get a 1D array of values on a 2D Gabor map at the set of points given
  //
  p1 = gabor_2d_space_raw_irreg(lay->irr_x,lay->irr_y,cnt,xc,yc,sdo,sdp,sf,
				ori,ph);

  p = get_2d_farray(cnt,1);
  for(i=0;i<cnt;i++)
    p[i][0] = p1[i];

  ppos = copy_2d_farray(p,cnt,1);
  absolute_value_2d_farray(ppos,cnt,1);

  // Make the largest value be 1, so 'eps' is relative to max=1.0
  make_max_const_2d_farray(ppos,cnt,1,1.0);

  // get connections, return number of possible contacts
  nposs = mod_conn_meth_02(mylogf,ppos,cnt,1,tflag,n,seed,eps,maxrep,&cx,&cy);

  if (evflag){ // Even out the random connections
    mod_conn_even_out(mylogf,p,cnt,1,seed,cx,cy,n,50,2.0,0);
  }

  // Count how many of each type, then extract into separate arrays
  cn0 = cn1 = 0;
  for(i=0;i<n;i++){
    if (p[cx[i]][cy[i]] > 0.0)
      cn1 += 1;
    else
      cn0 += 1;
  }
  cx0 = (int *)myalloc(cn0*sizeof(int));
  cy0 = (int *)myalloc(cn0*sizeof(int));
  cx1 = (int *)myalloc(cn1*sizeof(int));
  cy1 = (int *)myalloc(cn1*sizeof(int));

  k0 = k1 = 0;
  for(i=0;i<n;i++){
    if (p[cx[i]][cy[i]] > 0.0){
      cx1[k1] = cx[i];
      cy1[k1] = cy[i];
      k1 += 1;
    }else{
      cx0[k0] = cx[i];
      cy0[k0] = cy[i];
      k0 += 1;
    }
  }
  *rcx1 = cx1; *rcy1 = cy1; *rcn1 = cn1;
  *rcx0 = cx0; *rcy0 = cy0; *rcn0 = cn0;

  myfree(cx); myfree(cy);
  free_2d_farray(ppos,xn);
  free_2d_farray(p,xn);

  return nposs;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_CONN_GABOR_03                            */
/*                                                                           */
/*  Like MOD_CONN_GABOR_02, but pick separately for ON and OFF to balance.   */
/*                                                                           */
/*****************************************************************************/
int mod_conn_gabor_03(mylogf,xn,yn,xc,yc,sdo,sdp,sf,ph,ori,seed,maxrep,eps,n,
		       evflag,rcx1,rcy1,rcn1,rcx0,rcy0,rcn0)
     char *mylogf;  /* For cluster */
     int xn,yn;     /* get connections to points on grid 'xn' by 'yn' */
     float xc,yc;   /* center of Gabor on the grid */
     float sdo,sdp; /* SD of Gaussian (pix) orth and par to ori */
     float sf;      /* SF of sinusoid (cyc/pix) */
     float ph;      /* phase of sinusoid (degr) */
     float ori;     /* orientation */
     int seed;      /* randomization */
     int maxrep;    /* maximum times any one point can be chosen (>=1) */
     float eps;     /* no connections to point w/ lower prob <= eps */
     int n;         /* number of connections to make */
     int evflag;    /* 1-even out connections */
     int **rcx1,**rcy1;  /* Return ON connection coordinates */
     int *rcn1;          /* Number of ON */
     int **rcx0,**rcy0;  /* Return OFF connection coordinates */
     int *rcn0;          /* Number of OFF */
{
  int i;
  int *cx0,*cy0,*cx1,*cy1,no2,k0,k1,nposs0,nposs1,seed2,maxit;
  float **p,**ppos,**pneg;

  seed2 = (seed * 29) / 13;
  no2 = n/2;

  if (eps < 0.0)
    mylog_exit(mylogf,"MOD_CONN_GABOR_O3  eps must be >=0\n");

  /*printf("xn,yn,xc,yc,sd,sf,ori,ph = ");
    printf("%d %d %f %f %f %f %f %f\n",xn,yn,xc,yc,sd,sf,ori,ph);*/

  p = gabor_2d_space_raw(xn,yn,xc,yc,sdo,sdp,sf,ori,ph);
  ppos = copy_2d_farray(p,xn,yn);
  half_wave_rectify_2d_farray(ppos,xn,yn);

  pneg = copy_2d_farray(p,xn,yn);
  multiply_2d_farray(pneg,xn,yn,-1.0);
  half_wave_rectify_2d_farray(pneg,xn,yn);

  // Make the largest value be 1, so 'eps' is relative to max=1.0
  make_max_const_2d_farray(ppos,xn,yn,1.0);
  make_max_const_2d_farray(pneg,xn,yn,1.0);

  // get ON and OFF connections separately, return # of possible contacts
  nposs0 = mod_conn_meth_02(mylogf,pneg,xn,yn,NULL,no2,seed,eps,maxrep,
			    &cx0,&cy0);
  nposs1 = mod_conn_meth_02(mylogf,ppos,xn,yn,NULL,no2,seed,eps,maxrep,
			    &cx1,&cy1);
  if (evflag){ // Even out the random connections
    maxit = n/4;
    mod_conn_even_out(mylogf,pneg,xn,yn,seed,cx0,cy0,no2,maxit,2.0,0);
    mod_conn_even_out(mylogf,ppos,xn,yn,seed,cx1,cy1,no2,maxit,2.0,0);
  }

  *rcx1 = cx1; *rcy1 = cy1; *rcn1 = no2;
  *rcx0 = cx0; *rcy0 = cy0; *rcn0 = no2;

  free_2d_farray(ppos,xn);
  free_2d_farray(pneg,xn);
  free_2d_farray(p,xn);

  return nposs0 + nposs1;  /* Is this an over-estimage? */
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_CONN_GAUSS_02                             */
/*                                                                           */
/*  Get connections according to a Gaussian or disk function.                */
/*                                                                           */
/*  - Max value of gabor is set to 1                                         */
/*  - No connections will be made where value is <= eps                      */
/*                                                                           */
/*****************************************************************************/
int mod_conn_gauss_02(mylogf,xn,yn,shape,xc,yc,sdo,sdp,ori,seed,maxrep,eps,n,
		      evflag,eopar1,itmax,dumpflag,rcx,rcy,rcn)
     char *mylogf;  // Log file
     int xn,yn;     // size of grid
     char *shape;   // 'Gaussian' or 'disk'
     float xc,yc;   // center of Gaussian on the grid
     float sdo,sdp; // SD of Gaussian (pix) orth and par to ori
     float ori;     // orientation
     int seed;      // randomization
     int maxrep;    // maximum times any one point can be chosen (>=1)
     float eps;     // no connections to point w/ lower prob <= eps
     int n;         // number of connections to make
     int evflag;    // 1-even out connections
     float eopar1;  // Even out parameter - e.g., (pix) SD of smoothing
     int itmax;     // Max iterations for even-out algorithm, -1 use default
     int dumpflag;  // [0] 1-dump output for debugging, etc.
     int **rcx,**rcy;  // Return connection coordinates
     int *rcn;         // Number connections
{
  int *cx,*cy,nposs;
  float **p;

  if (eps < 0.0)
    mylog_exit(mylogf,"MOD_CONN_GAUSS_O2  eps must be >=0\n");

  /*printf("xn,yn,xc,yc,sd,ori, = ");
    printf("%d %d %f %f %f %f \n",xn,yn,xc,yc,sd,ori);*/

  if (strcmp(shape,"Gaussian")==0){
    p = gauss_2d_space_raw(xn,yn,xc,yc,sdo,sdp,ori);
  }else if (strcmp(shape,"disk")==0){
    p = disk_2d_space_raw(xn,yn,xc,yc,sdo);  // 'sdo' is Diameter of disk
    eps = 0.5;
  }else
    exit_error("MOD_CONN_GAUSS_02","Bad 'pbflag' value");

  // Make the largest value be 1, so 'eps' is relative to max=1.0
  make_max_const_2d_farray(p,xn,yn,1.0); 

  // get connections, return number of possible contacts
  nposs = mod_conn_meth_02(mylogf,p,xn,yn,NULL,n,seed,eps,maxrep,&cx,&cy);

  if (evflag == 1){ // Even out the random connections
    //mod_conn_even_out(mylogf,p,xn,yn,seed,cx,cy,n,itmax,2.0,0);
    mod_conn_even_out(mylogf,p,xn,yn,seed,cx,cy,n,itmax,eopar1,0);
  }else if (evflag == 2){ // Even out the random connections
    mod_conn_even_out_diffuse(mylogf,xn,yn,cx,cy,n,itmax,xc,yc,sdo,p,eps,
			      dumpflag);
  }else if (evflag != 0)
    mylog_exit(mylogf,"MOD_CONN_GAUSS_O2  Unknown 'evflag' value\n");

  free_2d_farray(p,xn);

  *rcx = cx; *rcy = cy; *rcn = n;

  return nposs;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_CONN_LAYER_GABOR_MAP                        */
/*                                                                           */
/*  Phases are set randomly to [0,90,180,270], controlled by 'phseed'.       */
/*                                                                           */
/*****************************************************************************/
void mod_conn_layer_gabor_map(mylogf,mppl,pl,orimap,xn,yn,conname)
     char *mylogf;
     struct param_pair_list *mppl; // Model parameters
     struct pop_layer *pl;         // layer of cells
     float **orimap;
     int xn,yn;                    // input Grid
     char conname[];
{
  int i,j,k,ic;
  int seed,nsamp,maxrep,nseed,*seedlist,phseed,evflag,nposs,phtype,blncflag;
  int ai_phase,ai_ori;
  float orideg,gabor_sf,gabor_sdo,gabor_sdp,ph,sscale,eps,sdo,sdp,sf;
  float *ranph,x0,y0,xf,yf,xc,yc,ph0,phstep,dph0;
  int *tcx0,*tcy0,tcn0,*tcx1,*tcy1,tcn1;
  char tstr[SLEN];
  struct pop_cell *c;

  mylog(mylogf,"  MOD_CONN_LAYER_GABOR_MAP\n");

  exit_error("WYETH HERE","lgn inputs will have 'null' layer, won't work?");

  ai_phase =  pop_cell_attrib_index_f(pl,"phase"); // Attribute index
  ai_ori   =  pop_cell_attrib_index_f(pl,"ori"); // Attribute index

  sprintf(tstr,"    Layer is %d x %d (x %d) units.\n",pl->xn,pl->yn,pl->zn);
  mylog(mylogf,tstr);

  sprintf(tstr,"%s_sf",conname);
  gabor_sf = paramfile_get_float_param_or_exit(mppl,tstr);
  sprintf(tstr,"%s_sd_orth",conname);
  gabor_sdo = paramfile_get_float_param_or_exit(mppl,tstr);
  sprintf(tstr,"%s_sd_par",conname);
  gabor_sdp = paramfile_get_float_param_or_exit(mppl,tstr);
  sprintf(tstr,"%s_seed",conname);
  seed = paramfile_get_int_param_or_exit(mppl,tstr);
  sprintf(tstr,"%s_phseed",conname);
  phseed = paramfile_get_int_param_or_exit(mppl,tstr);
  sprintf(tstr,"%s_eps",conname);
  eps = paramfile_get_float_param_or_exit(mppl,tstr);
  sprintf(tstr,"%s_nsamp",conname);
  nsamp = paramfile_get_int_param_or_exit(mppl,tstr);
  sprintf(tstr,"%s_maxrep",conname);
  maxrep = paramfile_get_int_param_or_exit(mppl,tstr);
  sprintf(tstr,"%s_balance",conname);
  blncflag = paramfile_get_int_param_default(mppl,tstr,0);

  // Controlling phase for multiple z-dimension neurons
  sprintf(tstr,"%s_phtype",conname);
  phtype = paramfile_get_int_param_default(mppl,tstr,0);
  sprintf(tstr,"%s_ph0",conname);
  ph0 = paramfile_get_float_param_default(mppl,tstr,0.0);
  sprintf(tstr,"%s_phstep",conname);
  phstep = paramfile_get_float_param_default(mppl,tstr,90.0);

  sscale = paramfile_get_float_param_or_exit(mppl,"sscale");

  sdo = gabor_sdo / sscale;  /* Make units pix */
  sdp = gabor_sdp / sscale;  /* Make units pix */
  sf = gabor_sf * sscale;  /* Make units cyc/pix */

  x0 = pl->x0; /* Old way, before onode */
  y0 = pl->y0;
  xf = pl->xf;
  yf = pl->yf;

  // Get a seed for each cell based on 'seed'
  nseed = pl->xn * pl->yn * pl->zn;
  seedlist = get_seeds(seed,1000000,nseed);  /* +1 for phase seed */

  ranph = get_random_floats(phseed,nseed);

  dph0 = 180.0;  /* To alternate initial phase, for phtype = 2*/

  ic = 0; // cell index
  for(i=0;i<pl->xn;i++){
    for(j=0;j<pl->yn;j++){
      orideg = orimap[i][j]/2.0; /* Ori for this location */
      if (phtype == 1)
	ph = ph0;
      else if (phtype == 2) // Checkerboard alternation of 0/180 phase
	if ((i+j)%2 == 0)
	  ph = ph0;
	else
	  ph = ph0 + dph0;
      else
	ph = 90.0 * (float)((int)(4.0 * ranph[ic]));
      xc = x0 + (float)i*xf;     /* Center for this location */
      yc = y0 + (float)j*yf;
      for(k=0;k<pl->zn;k++){     /* For each unit */
	c = &(pl->c[i][j][k]);
	//c->cxn = xn;  WYETH OLD
	//c->cyn = yn;
	// WYETH ATTRIB
	c->attrib_f[ai_ori] = orideg;
	//c->ori = orideg; /* Wyeth new, for gui color ori map, */

	// WYETH ATTRIB
	c->attrib_f[ai_phase] = ph;
	//c->phase = ph;   /*   and useful for synaptic connections */

	if (((i==8) && (j==8)) && strcmp(conname,"pop_ex_conn0")==0){
	  sprintf(tstr,"   WARNING:  evening out connections for 8,8\n");
	  mylog(mylogf,tstr);
	  evflag = 1;
	}else
	  evflag = 0;

	if (blncflag == 1)
	  nposs = mod_conn_gabor_03(mylogf,xn,yn,xc,yc,sdo,sdp,sf,ph,orideg,
				    seedlist[ic],maxrep,eps,nsamp,evflag,
				    &tcx1,&tcy1,&tcn1,&tcx0,&tcy0,&tcn0);
	else
	  nposs = mod_conn_gabor_02(mylogf,xn,yn,NULL,xc,yc,sdo,sdp,sf,ph,
				    orideg,
				    seedlist[ic],maxrep,eps,nsamp,evflag,
				    &tcx1,&tcy1,&tcn1,&tcx0,&tcy0,&tcn0);


	// mech is 'NULL', layer is NULL
	// I put an EXIT ERROR ABOVE to prevent this
	pop_cell_add_lgn_input(c,NULL,NULL,
			       tcx0,tcy0,tcn0,tcx1,tcy1,tcn1,nposs,-1.0);

	//OLD: c->nposs = nposs;

	ic += 1;
	if ((phtype == 1) || (phtype == 2))
	  ph += phstep;
      }
    }
  }
  myfree(seedlist);
  myfree(ranph);
}
/**************************************-**************************************/
/*                                                                           */
/*                  MOD_CONN_ONODE_LAYER_GABOR_MAP_SET_PHASE                 */
/*                                                                           */
/*  Set cell parameters:  phase and sz1                                      */
/*                                                                           */
/*****************************************************************************/
void mod_conn_onode_layer_gabor_map_set_phase(mylogf,io,pl,sscale)
     char *mylogf;
     struct onode *io;           // Input node
     struct pop_layer *pl;       // layer recieving inputs
     float sscale;
{
  int i,j,k,ic;
  int nseed,phseed,phtype,ai_phase,ai_sz1,ai_sf,var_scale;
  float gabor_sdo,gabor_sdp,ph,sdo,sdp,sz1,sf;
  float *ranph,ph0,phstep,dph0;
  struct pop_cell *c;

  mylog(mylogf,"  MOD_CONN_ONODE_LAYER_GABOR_MAP_SET_PHASE\n");

  ai_phase =  pop_cell_attrib_index_f(pl,"phase"); // Attribute index
  ai_sz1   =  pop_cell_attrib_index_f(pl,"sz1");   // Attribute index

  phseed    = onode_getpar_int_exit(io,"phseed");

  if (onode_item(io,"sf")==1){
    var_scale = 0;
    gabor_sdo = onode_getpar_flt_exit(io,"sd_orth");
    gabor_sdp = onode_getpar_flt_exit(io,"sd_par");
    sdo = gabor_sdo / sscale;  // Make units pix
    sdp = gabor_sdp / sscale;  // Make units pix
    if (sdo > sdp)
      sz1 = sdo;
    else
      sz1 = sdp;
  }else{
    var_scale = 1;
    gabor_sdo = onode_getpar_flt_exit(io,"sd_orth_c");  // sdo = c / SF
    gabor_sdp = onode_getpar_flt_exit(io,"sd_par_c");   // sdp = sdo
    ai_sf = pop_cell_attrib_index_f(pl,"sf");   // Attribute index
  }

  // Controlling phase for multiple z-dimension neurons
  phtype    = onode_getpar_int_dflt(io,"phtype",0);
  ph0       = onode_getpar_flt_dflt(io,"ph0",0.0);
  phstep    = onode_getpar_flt_dflt(io,"phstep",90.0);


  // Get a seed for each cell based on 'seed'
  nseed = pl->xn * pl->yn * pl->zn;
  ranph = get_random_floats(phseed,nseed);

  dph0 = 180.0;  // To alternate initial phase, for phtype = 2

  ic = 0; // cell index
  for(i=0;i<pl->xn;i++){
    for(j=0;j<pl->yn;j++){
      //orideg = orimap[i][j]/2.0; // Ori for this location
      if (phtype == 0)  // Random
	ph = 90.0 * (float)((int)(4.0 * ranph[ic]));
      else if (phtype == 1) // Const
	ph = ph0;
      else if (phtype == 2)
	if ((i+j)%2 == 0)
	  ph = ph0;
	else
	  ph = ph0 + dph0;
      else if (phtype == 4){  // Regular w/i 4 boxes
	ph = ph0 + i%2 * 180.0;
	ph += j%2 * 90.0;
      }else
	ph = 90.0 * (float)((int)(4.0 * ranph[ic]));


      for(k=0;k<pl->zn;k++){     // For each unit
	c = &(pl->c[i][j][k]);

	c->attrib_f[ai_phase] = ph;    // Set the phase

	if (var_scale == 1){
	  sf = c->attrib_f[ai_sf] * sscale;  // cyc/pix
	  sdo = gabor_sdo / sf;              // pix
	  if (gabor_sdp > 1.0)
	    sz1 = gabor_sdp * sdo;
	  else
	    sz1 = sdo;
	}

	c->attrib_f[ai_sz1]   = sz1;   // Set the size

	ic += 1;
	if ((phtype == 1) || (phtype == 2))
	  ph += phstep;
	
      }
    }
  }
  myfree(ranph);
}
/**************************************-**************************************/
/*                                                                           */
/*                       MOD_CONN_ONODE_LAYER_GABOR_MAP                      */
/*                                                                           */
/*  Phases are set randomly to [0,90,180,270], controlled by 'phseed'.       */
/*                                                                           */
/*****************************************************************************/
void mod_conn_onode_layer_gabor_map(mylogf,io,plgn,pl,mechr,xn,yn,sscale,
				    inindex)
     char *mylogf;
     struct onode *io;           // Input node
     struct pop_layer *plgn;     // LGN ON-layer (for cell subclass)
     struct pop_layer *pl;       // layer recieving inputs
     struct pop_mech *mechr;     // post-syn receptor mechanism
     int xn,yn;                  // input grid size
     float sscale;               // input grid scale
     short inindex;              // Input index number
{
  int i,j,k,ic;
  int seed,nsamp,nssmp,maxrep,nseed,*seedlist,evflag,nposs,blncflag,evx,evy;
  int **tflag,*tcx0,*tcy0,tcn0,*tcx1,*tcy1,tcn1,ai_phase,ai_ori,ai_sf,var_scale;
  int nsamp_c0;
  float gabor_sf,gabor_sdo,gabor_sdp,ph,eps,sdo,sdp,sf,sz1,xc,yc,ori,nsamp_c;
  float normw,weight,sdo_c0;
  char tstr[SLEN],*lgntype,*tclass,*cond;
  struct pop_cell *c;

  int irr_flag;

  mylog(mylogf,"  MOD_CONN_ONODE_LAYER_GABOR_MAP\n");

  //
  //  1.  MUST SET c->phase AND c->sz1 FOR ALL CELLS
  //
  mod_conn_onode_layer_gabor_map_set_phase(mylogf,io,pl,sscale);

  //
  //  2.  Ideally, only do things that are stored in a .conn file
  //

  sprintf(tstr,"    Layer is %d x %d (x %d) units.\n",pl->xn,pl->yn,pl->zn);
  mylog(mylogf,tstr);

  if (onode_item(io,"sf")==1){
    var_scale = 0;
    gabor_sf  = onode_getpar_flt_exit(io,"sf");
    gabor_sdo = onode_getpar_flt_exit(io,"sd_orth");
    gabor_sdp = onode_getpar_flt_exit(io,"sd_par");
    sdo = gabor_sdo / sscale;  // Make units pix
    sdp = gabor_sdp / sscale;  // Make units pix
    sf  = gabor_sf  * sscale;    // Make units cyc/pix
    nsamp_c = -1.0; // not used
    nsamp_c0 = 0.0; // not used
  }else{
    var_scale = 1;
    gabor_sdo = onode_getpar_flt_exit(io,"sd_orth_c");  // sdo = c / SF
    sdo_c0    = onode_getpar_flt_exit(io,"sd_orth_c0",0.0);  // sdo = c / SF
    gabor_sdp = onode_getpar_flt_exit(io,"sd_par_c");   // sdp = sdo
    nsamp_c   = onode_getpar_flt_dflt(io,"nsamp_c",2.0); // Scale const, var n
    nsamp_c0  = onode_getpar_int_dflt(io,"nsamp_c0",0); // Scale const, var n
  }

  //printf("sdo = %f   var_scale = %d\n",gabor_sdo,var_scale);

  seed      = onode_getpar_int_exit(io,"seed");
  eps       = onode_getpar_flt_exit(io,"eps");
  nsamp     = onode_getpar_int_exit(io,"nsamp");
  maxrep    = onode_getpar_int_exit(io,"maxrep");
  blncflag  = onode_getpar_int_dflt(io,"balance",0);
  lgntype   = onode_getpar_chr_dflt(io,"celltype","a");
  cond      = onode_getpar_chr_dflt(io,"cell_condition","all");

  evx    = onode_getpar_int_dflt(io,"even_x",-1);
  evy    = onode_getpar_int_dflt(io,"even_y",-1);

  normw = onode_getpar_flt_dflt(io,"normw",-1.0);



  //
  //  WYETH VARY these?
  //

  //
  //  3.  Check for irregular geometry in the pre-syn layer
  //
  if (strcmp(plgn->geomt,"irregular")==0){
    irr_flag = 1;
  }else
    irr_flag = 0;


  //
  //  Condition on cell type (e.g., L,M,S type cone classes, or Magno Parvo?)
  //  Make a constraint map 'tflag' based on cell type
  //
  if (irr_flag == 1){
    tflag = popc_layirr_get_conditional_tflag(plgn,cond);
  }else{
    tflag = get_const_2d_iarray(xn,yn,1);
    // WYETH - Sept 2008 - is this needed w/ newlgn ??
    // Presumably this was used to make sparser connections, as if to m-layers
    // which would have fewer cells, for Valeria's SIZE project, but perhaps
    // this won't be used again.
    if (strcmp(lgntype,"a")!=0){
      mylog_exit("MOD_CONN_ONODE_LAYER_GABOR_MAP  CHECK compat. w/ newlgn.");
      for(i=0;i<xn;i++){
	for(j=0;j<yn;j++){
	  tclass = plgn->c[i][j][0].subclass;
	  if (strcmp(tclass,lgntype)!=0){
	    tflag[i][j] = 0;
	  }
	}
      }
    }
  }

  // Get a seed for each cell based on 'seed'
  nseed = pl->xn * pl->yn * pl->zn;
  seedlist = get_seeds(seed,1000000,nseed);  // +1 for phase seed

  ai_phase =  pop_cell_attrib_index_f(pl,"phase"); // Attribute index
  ai_ori   =  pop_cell_attrib_index_f(pl,"ori"); // Attribute index
  if (var_scale == 1)
    ai_sf  =  pop_cell_attrib_index_f(pl,"sf"); // Attribute index
  nssmp = nsamp;

  weight = 1.0;  // Default weight for each LGN input

  ic = 0; // cell index
  for(i=0;i<pl->xn;i++){
    for(j=0;j<pl->yn;j++){
      for(k=0;k<pl->zn;k++){     // For each unit

	c = &(pl->c[i][j][k]);
	xc = c->rfx;  // WYETH FOR ocdom OD col use RFx
	yc = c->rfy;

	ph = c->attrib_f[ai_phase];
	ori = c->attrib_f[ai_ori];
	if (var_scale == 1){
	  sf = c->attrib_f[ai_sf] * sscale;  // cyc/pix
	  sdo = gabor_sdo / sf;                // pix
	  sdo += sdo_c0;
	  sdp = gabor_sdp * sdo;               // pix

	  if (nsamp < 0){
	    nssmp = sdo*sdp * nsamp_c;  // 'nsamp' determined by RF area
	    nssmp += nsamp_c0;

	    if (normw > 0.0){
	      weight = normw / (float)nssmp; // WYETH, LGNW
	    }else{
	      weight = 1.0;
	    }
	  }
	  //if (j == pl->yn/2)
	  //printf("i = %d   sdo,sdp = %f  %f   sf = %f\n",i,sdo,sdp,sf);
	}


	if ((i==evx) && (j==evy)){
	  sprintf(tstr,"    Evening out connections for %d,%d\n",i,j);
	  mylog(mylogf,tstr);
	  evflag = 1;
	}else
	  evflag = 0;

	if (blncflag == 1){
	  nposs = mod_conn_gabor_03(mylogf,xn,yn,xc,yc,sdo,sdp,sf,ph,ori,
				    seedlist[ic],maxrep,eps,nssmp,evflag,
				    &tcx1,&tcy1,&tcn1,&tcx0,&tcy0,&tcn0);
	}else if (irr_flag == 1){
	  nposs = mod_conn_gabor_02_irr(mylogf,plgn,xn,yn,tflag,xc,yc,sdo,sdp,
					sf,ph,ori,seedlist[ic],maxrep,eps,
					nssmp,evflag,cond,
					&tcx1,&tcy1,&tcn1,&tcx0,&tcy0,&tcn0);
	}else{
	  nposs = mod_conn_gabor_02(mylogf,xn,yn,tflag,xc,yc,sdo,sdp,sf,ph,ori,
				    seedlist[ic],maxrep,eps,nssmp,evflag,
				    &tcx1,&tcy1,&tcn1,&tcx0,&tcy0,&tcn0);
	}


	if (irr_flag == 1){
	  pop_cell_add_on_off_input_irr(c,plgn,tcx0,tcy0,tcn0,tcx1,tcy1,tcn1,
					nposs,inindex);
	}else{
	  // 11 Jan 2009
	  // Do not free the connection coord storage, it is attached by ptr
	  pop_cell_add_lgn_input(c,plgn,mechr,tcx0,tcy0,tcn0,tcx1,tcy1,tcn1,
				 nposs,weight);
	  //c->nposs = nposs;
	}

	ic += 1;

	// WYETH TEMP, won't need to check name, once there is layer pointer
	if ((strcmp(plgn->name,"rgc_m")==0)||(strcmp(plgn->name,"lgn_m")==0)){
	  c->lgn_mflag = 1;
	}
      }
    }
  }
  myfree(seedlist);
  myfree(lgntype);
  myfree(cond);
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_CONN_LGN_PAIR                             */
/*                                                                           */
/*  Make paired connections from two LGN populations.                        */
/*                                                                           */
/*****************************************************************************/
void mod_conn_lgn_pair(mpt,io,lpost,lgn1,lgn2,mechr,mechp)
     struct pop_top *mpt;        // Global model structure
     struct onode *io;           // Input node
     struct pop_layer *lpost;    // post-syn layer
     struct pop_layer *lgn1;     // pre-syn LGN layer 1
     struct pop_layer *lgn2;     // pre-syn LGN layer 2
     struct pop_mech *mechr;     // post-syn receptor mechanism
     struct pop_mech *mechp;     // post-syn pairing mechanism
{
  int i,j,k,ic;
  int seed,cseed,nseed,*seedlist,maxrep,evflag,nposs,self;
  int xn,yn,*tx0,*ty0,*tx1,*ty1,tn0,tn1,*tx0p,*ty0p,*tx1p,*ty1p;
  int nsamp,tn,jit_seed,zlay1,zlay2;
  float sscale,vdir,vdist,jit_para,jit_orth,vdistpix,jitp,jito;
  float cdist,sd,minw,normw,prob,xc,yc;
  char tstr[SLEN],*tclass,*mylogf;
  struct pop_cell *c;
  struct onode *disto;   // Distrib node
  struct pop_lgn_in *lgnin;

  mylogf = mpt->logf;

  mylog(mylogf,"  MOD_CONN_LGN_PAIR\n");

  xn     = mpt->xn;
  yn     = mpt->yn;
  sscale = mpt->sscale;

  vdir      = onode_getpar_flt_exit(io,"offset_direction");   // Vector dir.
  vdist     = onode_getpar_flt_exit(io,"offset_distance");    // Vector length
  jit_para  = onode_getpar_flt_exit(io,"offset_jitter_para"); // parallel jit.
  jit_orth  = onode_getpar_flt_exit(io,"offset_jitter_orth"); // orthog jitter
  jit_seed  = onode_getpar_int_exit(io,"offset_jitter_seed"); // jitter seed


  zlay1  = onode_getpar_flt_dflt(io,"pop_layer_z_1",-1); // z-layer
  zlay2  = onode_getpar_flt_dflt(io,"pop_layer_z_2",-1); // z-layer

  if ((zlay1 != -1) || (zlay2 != -1))
    mylogx(mylogf,"MOD_CONN_LGN_PAIR","zlay impl'd only for value -1");

  vdistpix = vdist / sscale;  // vector length (pix)
  jitp = jit_para / sscale;
  jito = jit_orth / sscale;

  // WYETH - THIS SHOULD BE A SEPARATE MODULE
  // WYETH - THIS SHOULD BE A SEPARATE MODULE
  // WYETH - THIS SHOULD BE A SEPARATE MODULE
  disto = onode_child_get_unique(io,"distrib");
  cdist = onode_getpar_flt_exit(disto,"cdist");
  minw  = onode_getpar_flt_dflt(disto,"minw",0.0);
  normw = onode_getpar_flt_dflt(disto,"normw",0.0); // UNUSED HERE??
  seed  = onode_getpar_int_dflt(disto,"seed",1777);
  nsamp = onode_getpar_flt_exit(disto,"n");
  self  = onode_getpar_int_dflt(disto,"self",1);
  maxrep = 1;
  evflag = 0;

  sd = cdist / sscale;  // Gaussian SD (pix)

  // Get a seed for each post-syn cell based on 'seed'
  nseed = lpost->xn * lpost->yn * lpost->zn;
  seedlist = get_seeds(seed,1000000,nseed);

  ic = 0; // cell index
  for(i=0;i<lpost->xn;i++){
    for(j=0;j<lpost->yn;j++){
      for(k=0;k<lpost->zn;k++){     // For each unit

	c = &(lpost->c[i][j][k]);

	cseed = seedlist[ic];

	xc = c->rfx;  // WYETH FOR ocdom OD col use RFx
	yc = c->rfy;

	// Pick OFF coords
	nposs = mod_conn_gauss_02(mylogf,xn,yn,"Gaussian",xc,yc,sd,sd,0.0,
				  cseed,maxrep,minw,nsamp,evflag,2.0,-1,0,
				  &tx0,&ty0,&tn0);
	// Pick ON coords
	nposs = mod_conn_gauss_02(mylogf,xn,yn,"Gaussian",xc,yc,sd,sd,0.0,
				  cseed*27,maxrep,minw,nsamp,evflag,2.0,-1,0,
				  &tx1,&ty1,&tn1);

	// Paired OFF
	mod_conn_pick_paired_connections(tx0,ty0,tn0,xn,yn,vdir,vdistpix,jitp,
					 jito,jit_seed,&tx0p,&ty0p);
	// Paired ON
	mod_conn_pick_paired_connections(tx1,ty1,tn1,xn,yn,vdir,vdistpix,jitp,
					 jito,jit_seed,&tx1p,&ty1p);

	// Do not free the connection coord storage, it is attached by ptr
	pop_cell_add_lgn_input(c,lgn1,mechr,tx0,ty0,tn0,tx1,ty1,tn1,nposs,
			       -1.0);

	lgnin = pop_cell_get_lgn_input_by_lay(c,lgn1);

	pop_cell_add_lgn_pair(lgnin,lgn2,mechp,tx0p,ty0p,tn0,tx1p,ty1p,tn1,0);

	ic += 1;
      }
    }
  }
  myfree(seedlist);
}
/**************************************-**************************************/
/*                                                                           */
/*                               MOD_CONN_RFX_01                             */
/*                                                                           */
/*  Compute a measure of the correlation of the RF masks.                    */
/*                                                                           */
/*****************************************************************************/
float mod_conn_rfx_01(m1,m2,xn,yn)
     float **m1,**m2;  // RF masks [xn][yn]
     int xn,yn;
{
  int i,j;
  float rfx;
  
  rfx = 0.0;
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      rfx += m1[i][j] * m2[i][j];
    }
  }
  return rfx;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_CONN_CLIST_COUNT_SHARED                      */
/*                                                                           */
/*****************************************************************************/
int mod_conn_clist_count_shared(cx1,cy1,cn1,cx2,cy2,cn2)
     int *cx1,*cy1,cn1;  // 1st list of connections
     int *cx2,*cy2,cn2;  // 2nd list of connections
{
  int i,j;
  int n;
  int x,y;

  n = 0;
  for(i=0;i<cn1;i++){
    x = cx1[i];
    y = cy1[i];
    for(j=0;j<cn2;j++){
      if ((x == cx2[j]) && (y == cy2[j]))
	n += 1;
    }
  }
  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_CONN_INPUT_PLAIN                          */
/*                                                                           */
/*  Process 'plain' inputs, i.e., ones that do not relate to 'pop' models.   */
/*                                                                           */
/*  Developed for use with 'mod_me_util' MT models.                          */
/*                                                                           */
/*****************************************************************************/
void mod_conn_input_plain(mylogf,io)
     char *mylogf;
     struct onode *io;              // Input node
{
  struct onode *disto;
  
  mylog(mylogf,"  MOD_CONN_INPUT_PLAIN\n");

  disto = onode_child_get_unique(io,"distrib");

  /*
  minw  = onode_getpar_flt_dflt(disto,"minw",0.0);
  cdist = onode_getpar_flt_dflt(disto,"cdist",0.0);
  prob  = onode_getpar_flt_dflt(disto,"prob",1.0);
  seed  = onode_getpar_int_dflt(disto,"seed",1777);
  csign = onode_getpar_int_dflt(disto,"corr_sign",1);
  normw = onode_getpar_flt_dflt(disto,"normw",0.0);
  self  = onode_getpar_int_dflt(disto,"self",1);
  shape = onode_getpar_chr_dflt(disto,"shape","Gaussian");
  */
  
  exit(0);
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_CONN_CORR_ON_OFF                           */
/*                                                                           */
/*  Connect the inhibitory cells (lay2) to the excitatory cells (lay1)       */
/*  if their RFs are anticorrelated.  The synaptic weight 'w' will range     */
/*  from 1 (if all points in the inhib RF are opposite those of the excit    */
/*  RF) to 0 for pairs of cells with no correlation or positive correlation. */
/*                                                                           */
/*  When the LGN -> CTX connections are sparse, this algorithm may not       */
/*  work so well.                                                            */
/*                                                                           */
/*****************************************************************************/
void mod_conn_corr_on_off(mylogf,lay1,lay2,xn,yn,inindex)
     char *mylogf;
     struct pop_layer *lay1,*lay2;
     int xn,yn;
     short inindex;           // Index into population 'inlist'
{
  int i,ai,aj,ak,bi,bj,bk;
  int n,*x,*y,cnt1,cnt2;
  int xn1,yn1,zn1,xn2,yn2,zn2,nsyn;
  float **cm,tot,w;
  struct pop_cell *pc1,*pc2;
  char tstr[SLEN];
  struct pop_lgn_in *ln1,*ln2;


  mylog(mylogf,"  MOD_CONN_CORR_ON_OFF\n");

  xn1 = lay1->xn;  yn1 = lay1->yn;  zn1 = lay1->zn;
  xn2 = lay2->xn;  yn2 = lay2->yn;  zn2 = lay2->zn;

  cm = get_zero_2d_farray(xn,yn);

  nsyn = 0;
  cnt1 = cnt2 = 0; // Cell number in lay1 and lay2
  for(ai=0;ai<xn1;ai++){
    for(aj=0;aj<yn1;aj++){
      for(ak=0;ak<zn1;ak++){  /* For each cell in the first layer */
	pc1 = &(lay1->c[ai][aj][ak]);

	// WYETH HERE - ASSUMING ONLY ONE LGN INPUT
	if (pc1->lgn_n > 0)
	  ln1 = pc1->lgnin[0];
	else
	  exit_error("MOD_CONN_CORR_ON_OFF","No LGN inputs");


	/* Fill in connection matrix */
	n = ln1->cn1;
	x = ln1->cx1;
	y = ln1->cy1;
	for(i=0;i<n;i++)
	  cm[x[i]][y[i]] += 1.0;  /* ON connection */
	
	n = ln1->cn0;
	x = ln1->cx0;
	y = ln1->cy0;
	for(i=0;i<n;i++)
	  cm[x[i]][y[i]] -= 1.0;  /* OFF connection */

	cnt2 = 0;
	for(bi=0;bi<xn2;bi++){
	  for(bj=0;bj<yn2;bj++){
	    for(bk=0;bk<zn2;bk++){  /* For each cell in the 2nd layer */
	      tot = 0.0;
	      pc2 = &(lay2->c[bi][bj][bk]);

	      // WYETH HERE - ASSUMING ONLY ONE LGN INPUT
	      if (pc2->lgn_n > 0)
		ln2 = pc2->lgnin[0];
	      else
		exit_error("MOD_CONN_CORR_ON_OFF","No LGN inputs");

	      n = ln2->cn1;
	      x = ln2->cx1;
	      y = ln2->cy1;
	      for(i=0;i<n;i++)
		tot += cm[x[i]][y[i]];  /* ON connection */

	      n = ln2->cn0;
	      x = ln2->cx0;
	      y = ln2->cy0;
	      for(i=0;i<n;i++)
		tot -= cm[x[i]][y[i]];  /* OFF connection */

	      w = tot/(float)(ln2->cn1 + ln2->cn0);

	      if (w < 0.0){
		(void)pop_cell_add_synapse(pc2,pc1,2,-w,0.0,inindex);
		nsyn += 1;
		//printf(" %4d %4d  tot = %4f   w_In= %f\n",cnt1,cnt2,tot,-w);
	      }
	      cnt2 += 1;
	    }
	  }
	}

	// Reset connection matrix to 0's
	n = ln1->cn1;
	x = ln1->cx1;
	y = ln1->cy1;
	for(i=0;i<n;i++)
	  cm[x[i]][y[i]] = 0.0;  /* Reset ON connection */
	
	n = ln1->cn0;
	x = ln1->cx0;
	y = ln1->cy0;
	for(i=0;i<n;i++)
	  cm[x[i]][y[i]] = 0.0;  /* Reset OFF connection */

	cnt1 += 1;
      }
    }
  }
  free_2d_farray(cm,xn);

  sprintf(tstr,"    %d synapses from inhibitory to excitatory neurons\n",nsyn);
  mylog(mylogf,tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_CONN_CORR_W_NORM_HIST                        */
/*                                                                           */
/*****************************************************************************/
void mod_conn_corr_w_norm_hist(mylogf,w,xn,yn,zn,minw,phase_flag)
     char *mylogf;
     float ***w;      // Weights (e.g., raw correlation values)
     int xn,yn,zn;
     float minw;      // Don't make connections to wieghts less than this
     int phase_flag;  // -1: anticorr, otherwise correlated
{
  int i,j,k;
  int histn;
  float tw,tot,wmax,*histw,mean,sdev,wcrit,scale;
  
  // Truncate negative (positive) weights
  tot = 0.0;
  wmax = 0.0;
  
  if (phase_flag == -1){
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	for(k=0;k<zn;k++){
	  tw = w[i][j][k];
	  if (tw > 0.0)
	    w[i][j][k] = 0.0;  // positive corr gets 0 syn weight
	  else{
	    w[i][j][k] = -tw;  // negative corr gets positive weight
	    tot -= tw;
	    if (-tw > wmax)
	      wmax = -tw;
	  }
	}
  }else{  // Connect positive correlated RFs
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	for(k=0;k<zn;k++){
	  tw = w[i][j][k];
	  if (tw < 0.0)
	    w[i][j][k] = 0.0;  // negative corr gets 0 syn weight
	  else{
	    w[i][j][k] = tw;   // positive corr gets positive weight
	    tot += tw;
	    if (tw > wmax)
	      wmax = tw;
	  }
	}
  }


  // Make a 1D array of weights, get mean and SD
  histw = (float *)myalloc(xn*yn*zn*sizeof(float));
  histn = 0;
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	tw = w[i][j][k];
	if (tw > 0.0){
	  histw[histn] = tw;
	  histn += 1;
	}
      }
    }
  }
  mean_sdev_farray(histw,histn,&mean,&sdev);
  wcrit = mean + 2.0 * sdev;
  myfree(histw);

  //printf("  MEAN, SD  %f %f\n",mean,sdev);
  
  /*** Now:
       1.  Scale weights so 'wcrit' becomes 1.0
       2.  Set all weights larger than 1.0 to 1.0 ***/
  
  for(i=0;i<xn;i++){            // Rescale weights
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	tw = w[i][j][k];
	if (tw > wcrit)
	  w[i][j][k] = 1.0;  // Limit high values to 1.0
	else if (tw/wcrit < minw)
	  w[i][j][k] = 0.0;  // Allow suppression of small weights
	else
	  w[i][j][k] = tw/wcrit;  // Linearly scale lower values
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_CONN_COMPOSITE_RF_IRR                        */
/*                                                                           */
/*  Return a 2d composit RF map for irregular LGN inputs.                    */
/*                                                                           */
/*****************************************************************************/
float **mod_conn_composite_rf_irr(xn,yn,sig,c,flag)
     int xn,yn;               // size of field
     float sig;               // SD for smoothing, xn,yn units
     struct pop_cell *c;      // Cell for which to build the RF
     int flag;                // 1 - L-ON/OFF is 1/0, M-ON/OFF is 0/1
{
  int i;
  int cid,onflag,wn,*wx,*wy;
  float **rf,*w;
  struct pop_syn *tsyn;
  struct pop_cell *ct;

  rf = get_zero_2d_farray(xn,yn);

  if (flag != 1)  // No other flag options used yet.
    exit_error("MOD_CONN_COMPOSITE_RF_IRR","Flag value not imp'd yet");

  //
  //  Search presyn list for LGN inputs.
  //
  tsyn = c->in;  // Input synapses
  while (tsyn != NULL){
    ct  = tsyn->pre;  // Pre-syn unit

    // *** WYETH should check that this 'ct' UNIT is an LGN INPUT ***
    // *** WYETH should check that this 'ct' UNIT is an LGN INPUT ***
    //  wyeth - search the onodes in the 'inlist' for ones that point to
    //  'lgn_on_off' type of inputs with 'pop_origin ...' layers that are
    //  'irregular' geomt ??? IS THIS THE BEST WAY TO FIGURE THIS OUT??
    //  is there not a query for 'is lgn input?' or 'is irreg input' ????

    cid = my_rint(pop_cell_attrib_get_f(ct,"conetype"));
    onflag = ct->layz;  // z-layer index tells ON/OFF = 1/0

    // Get coords and Gaussian weights on regular grid
    retutil_irr2sq_mask(ct->cx,ct->cy,xn,yn,sig,2.6,&wn,&wx,&wy,&w);

    if (((cid == 2) && (onflag == 1)) ||  // Red ON
	((cid == 1) && (onflag == 0))){   // Green OFF
      for(i=0;i<wn;i++)
	rf[wx[i]][wy[i]] += w[i];  // treat as ON
    }else{
      for(i=0;i<wn;i++)
	rf[wx[i]][wy[i]] -= w[i];  // treat as OFF
    }
    myfree(wx); myfree(wy); myfree(w);

    tsyn = tsyn->post_next;
  }

  return rf;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_CONN_COMPOSITE_RF                         */
/*                                                                           */
/*  Return a 2d composit RF map for the LGN inputs.                          */
/*                                                                           */
/*****************************************************************************/
float **mod_conn_composite_rf(tmap,xn,yn,sig,ln)
     float **tmap;             // [xn][yn] all zeros
     int xn,yn;                // size of field
     float sig;                // SD for smoothing, xn,yn units
     struct pop_lgn_in *ln;    // LGN input
{
  int i;
  int n,*x,*y;
  float **rf;

  n = ln->cn1;
  x = ln->cx1;
  y = ln->cy1;
  for(i=0;i<n;i++)
    tmap[x[i]][y[i]] += 1.0;  // ON connection

  n = ln->cn0;
  x = ln->cx0;
  y = ln->cy0;
  for(i=0;i<n;i++)
    tmap[x[i]][y[i]] -= 1.0;  // OFF connection

  rf = smooth_2d_with_gaussian(tmap,xn,yn,sig,0.05);

  // Reset connection matrix to 0's
  n = ln->cn1;
  x = ln->cx1;
  y = ln->cy1;
  for(i=0;i<n;i++)
    tmap[x[i]][y[i]] = 0.0;  // Reset ON connection
  
  n = ln->cn0;
  x = ln->cx0;
  y = ln->cy0;
  for(i=0;i<n;i++)
    tmap[x[i]][y[i]] = 0.0;  // Reset OFF connection
  
  return rf;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_CONN_CORR_ON_OFF_01                          */
/*                                                                           */
/*  Connect the inhibitory cells (lay2) to the excitatory cells (lay1)       */
/*  if their RFs are anticorrelated.                                         */
/*                                                                           */
/*  *** This algorithm should work better than the previous one for          */
/*      sparse connections.                                                  */
/*                                                                           */
/*****************************************************************************/
void mod_conn_corr_on_off_01(mylogf,lay1,lay2,xn,yn,sig,minw,prob,seed,
			     phase_flag,syntype,cdist,axonv,inindex)
     char *mylogf;
     struct pop_layer *lay1;  // post-syn
     struct pop_layer *lay2;  // pre-syn
     int xn,yn;
     float sig;       // SD for Gaussian weight of RF
     float minw;      // Don't make connections where weight is less than this
     float prob;      // Less than 1.0 for random choice of synapse
     int seed;
     int phase_flag;  // -1: anticorr, otherwise correlated
     int syntype;     // See ifc.h for definitions
     float cdist;     // Cutoff distance, ignored if <= 0.0
     float axonv;     // Axon velocity, m/s, ignored if <= 0.0
     short inindex;   // Index into population 'inlist'
{
  int i,ai,aj,ak,bi,bj,bk;
  int n,*x,*y,xn1,yn1,zn1,xn2,yn2,zn2,nsyn,pflag,irr_flag;
  float **cm,w,tdelay,dist;
  float *****lay1rf,*****lay2rf,***rfx;
  double spum;
  char tstr[SLEN],ggstr[SLEN];
  struct pop_cell *c,*c1,*c2;
  struct pop_layer *lgn_lay;
  struct pop_lgn_in *ln;

  mylog(mylogf,"  MOD_CONN_CORR_ON_OFF_01\n");
  pflag = 0;


  //
  //  Examine the first cell, and establish the number and type of LGN input
  //
  c = &(lay1->c[0][0][0]);
  if (c == NULL)
    exit_error("MOD_CONN_CORR_ON_OFF_01","First unit pointer is null");

  if (c->lgn_n == 0){
    exit_error("MOD_CONN_CORR_ON_OFF_01","No LGN inputs");
  }else if (c->lgn_n > 1){
    sprintf(ggstr,"    There are %d LGN inputs to the first unit\n",c->lgn_n);
    mylog(mylogf,ggstr);
  }

  if (c->lgnin == NULL){    // 'lgnin' is not used for irregular LGN
    irr_flag = 1;
  }else
    irr_flag = 0;


  if (axonv > 0.0)
    spum = 1.0 / (axonv * 1000000.0);  // Seconds per Micron
  else
    spum = 0.0;

  if (pflag){
    printf("laynames = %s %s\n",lay1->name,lay2->name);
    printf("xn,yn = %d %d\n",xn,yn);
    printf("sig  = %f\n",sig);
    printf("minw = %f\n",minw);
    printf("prob = %f\n",prob);
    printf("seed = %d\n",seed);
    printf("phase_flag = %d\n",phase_flag);
    printf("syntype = %d\n",syntype);
  }

  sprintf(tstr,"    Gaussian SD = %.2f grid units\n",sig);
  mylog(mylogf,tstr);
  sprintf(tstr,"    Minimum normalized correlation for connection = %f\n",
	  minw);
  mylog(mylogf,tstr);

  xn1 = lay1->xn;  yn1 = lay1->yn;  zn1 = lay1->zn;
  xn2 = lay2->xn;  yn2 = lay2->yn;  zn2 = lay2->zn;

  cm = get_zero_2d_farray(xn,yn);
  lay1rf = get_5d_farray(xn1,yn1,zn1,-1,-1);
  lay2rf = get_5d_farray(xn2,yn2,zn2,-1,-1);

  if (pflag) printf("HERE 00\n");

  // WYETH - take the first cell, find it's LGN inputs, and then assess what
  // type this is.

  nsyn = 0;
  for(ai=0;ai<xn1;ai++){
    for(aj=0;aj<yn1;aj++){
      for(ak=0;ak<zn1;ak++){  // For each cell in the first layer

	// *****
	// WYETH - BUG FIX  Jan 6, 2012 (bug dates back to at least Aug 2006)
	// *****   (I commented this line out - was wasting memory)
	//cm = get_zero_2d_farray(xn,yn);  // DO NOT KEEP MALLOC'ING

	c = &(lay1->c[ai][aj][ak]);

	if (irr_flag == 1){
	  lay1rf[ai][aj][ak] = mod_conn_composite_rf_irr(xn,yn,sig,c,1);

	  /*** DUMP TEMPLATES
	  if ((ai == 0) && (aj == 0)){
	    if (ak == 0)
	      write_2d_data("z.t.0.2d",lay1rf[ai][aj][ak],0,0,xn,yn,4,2,0,1);
	    if (ak == 1)
	      write_2d_data("z.t.1.2d",lay1rf[ai][aj][ak],0,0,xn,yn,4,2,0,1);
	    if (ak == 2)
	      write_2d_data("z.t.2.2d",lay1rf[ai][aj][ak],0,0,xn,yn,4,2,0,1);
	    if (ak == 3)
	      write_2d_data("z.t.3.2d",lay1rf[ai][aj][ak],0,0,xn,yn,4,2,0,1);
	  }
	  ***/

	}else{
	  if (c->lgn_n > 0)
	    ln = c->lgnin[0];   // *** WYETH - ASSUMING LGNIN[0]
	  else
	    exit_error("MOD_CONN_CORR_ON_OFF_01","No LGN inputs");

	  lay1rf[ai][aj][ak] = mod_conn_composite_rf(cm,xn,yn,sig,ln);
	}
      }
    }
  }
  if (pflag) printf("HERE 01\n");

  for(ai=0;ai<xn2;ai++){
    for(aj=0;aj<yn2;aj++){
      for(ak=0;ak<zn2;ak++){  // For each cell in the second layer
	
	c = &(lay2->c[ai][aj][ak]);

	if (irr_flag == 1){
	  lay2rf[ai][aj][ak] = mod_conn_composite_rf_irr(xn,yn,sig,c,1);
	}else{
	  if (c->lgn_n > 0)
	    ln = c->lgnin[0];  // *** WYETH - ASSUMING LGNIN[0]
	  else
	    exit_error("MOD_CONN_CORR_ON_OFF_01","No LGN inputs");
	  lay2rf[ai][aj][ak] = mod_conn_composite_rf(cm,xn,yn,sig,ln);
	}
      }
    }
  }

  if (pflag){
    printf("WYETH - dumping intermediate 2d data\n");

    write_2d_data("zzz.1.2d",lay1rf[0][0][0],0,0,xn,yn,4,2,0,1);
    write_2d_data("zzz.2.2d",lay2rf[0][0][0],0,0,xn,yn,4,2,0,1);
  }

  if (pflag) printf("HERE 02\n");

  // Store correlation value with each inhib cell
  rfx = get_zero_3d_farray(xn2,yn2,zn2);

  if (prob < 1.0){
    if (seed > 0)
      seed *= -1;
  }

  for(ai=0;ai<xn1;ai++){
    for(aj=0;aj<yn1;aj++){
      for(ak=0;ak<zn1;ak++){  // For each cell in the first layer
	c1 = &(lay1->c[ai][aj][ak]);

	if (pflag) printf("%d %d %d \n",ai,aj,ak);

	for(bi=0;bi<xn2;bi++){
	  for(bj=0;bj<yn2;bj++){
	    for(bk=0;bk<zn2;bk++){  // For each cell in the 2nd layer
	      // Compute correlation of RF masks
	      if (cdist > 0.0){
		c2 = &(lay2->c[bi][bj][bk]);
		dist = mod_conn_get_distance_cells_layer(c1,c2);
		if (dist <= cdist)
		  rfx[bi][bj][bk] = mod_conn_rfx_01(lay1rf[ai][aj][ak],
						    lay2rf[bi][bj][bk],xn,yn);
		else
		  rfx[bi][bj][bk] = 0.0;  // Set weight to zero
	      }else
		rfx[bi][bj][bk] = mod_conn_rfx_01(lay1rf[ai][aj][ak],
						  lay2rf[bi][bj][bk],xn,yn);
	    }
	  }
	}

	//
	//  Put a ceiling of mean + 2*SD on weights, normalize, apply minw
	//
	mod_conn_corr_w_norm_hist(mylogf,rfx,xn2,yn2,zn2,minw,phase_flag);
	
	// Make connections
	tdelay = 0.0;
	for(bi=0;bi<xn2;bi++){
	  for(bj=0;bj<yn2;bj++){
	    for(bk=0;bk<zn2;bk++){  /* For each cell in the 2nd layer */
	      c2 = &(lay2->c[bi][bj][bk]);
	      w = rfx[bi][bj][bk];
	      if (w > 0.0){
		
		if (axonv > 0.0){
		  dist = mod_conn_get_distance_cells_layer(c1,c2);
		  tdelay = dist * spum;
		}
		
		if (prob < 1.0){
		  if (myrand_util_ran2(&seed) < prob){
		    (void)pop_cell_add_synapse(c2,c1,syntype,w,tdelay,inindex);
		    nsyn += 1;
		  }
		}else{
		  (void)pop_cell_add_synapse(c2,c1,syntype,w,tdelay,inindex);
		  nsyn += 1;
		}
	      }
	    }
	  }
	}
	
      }
    }
  }

  free_2d_farray(cm,xn);
  free_3d_farray(rfx,xn2,yn2,zn2);
  free_5d_farray(lay1rf,xn1,yn1,zn1,xn,yn);
  free_5d_farray(lay2rf,xn2,yn2,zn2,xn,yn);

  sprintf(tstr,"    %d synapses from inhibitory to excitatory neurons\n",nsyn);
  mylog(mylogf,tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                        MOD_CONN_CORR_ON_OFF_TEMPLATE                      */
/*                                                                           */
/*  Connect the pre-syn cells (lay2) to the post-syn cells (lay1)            */
/*  if their RFs are anticorrelated.                                         */
/*                                                                           */
/*****************************************************************************/
void mod_conn_corr_on_off_template(mylogf,lay1,lay2,lgn_name,xn,yn,sig,minw,
				   prob,seed,phase_flag,syntype,cdist,axonv,
				   inindex,sf,sdo,sdp,phase,tmpl_fname)
     char *mylogf;
     struct pop_layer *lay1;  // post-syn
     struct pop_layer *lay2;  // pre-syn
     char *lgn_name;          // LGN pop name
     int xn,yn;
     float sig;       // SD for Gaussian weight of RF
     float minw;      // Don't make connections where weight is less than this
     float prob;      // Less than 1.0 for random choice of synapse
     int seed;
     int phase_flag;  // -1: anticorr, otherwise correlated
     int syntype;     // See ifc.h for definitions
     float cdist;     // Cutoff distance, ignored if <= 0.0
     float axonv;     // Axon velocity, m/s, ignored if <= 0.0
     short inindex;           // Index into population 'inlist'
     float sf;        // template SF (cyc/pix)
     float sdo;       // template SD orthogonal (pix)
     float sdp;       // template SD parallel (pix)
     float phase;     // template phase (deg)
     char *tmpl_fname; // Filename to write template, or NULL
{
  int i,j,k,ai,aj,ak,bi,bj,bk;
  int n,*x,*y,ai_ori;
  int xn1,yn1,zn1,xn2,yn2,zn2,nsyn,histn,pflag;
  float **cm,tot,w,wmax,wcrit,*histw,mean,sdev,tdelay,dist,ori,xc,yc;
  float *****lay1rf,*****lay2rf,***rfx,**gabor;
  double spum;
  struct pop_cell *c,*c1,*c2;
  char tstr[SLEN];
  struct pop_lgn_in *ln;

  mylog(mylogf,"  MOD_CONN_CORR_ON_OFF_TEMPLATE\n");

  pflag = 0;

  if (axonv > 0.0)
    spum = 1.0 / (axonv * 1000000.0);  // Seconds per Micron
  else
    spum = 0.0;

  if (pflag){
    printf("laynames = %s %s\n",lay1->name,lay2->name);
    printf("xn,yn = %d %d\n",xn,yn);
    printf("sig  = %f\n",sig);
    printf("minw = %f\n",minw);
    printf("prob = %f\n",prob);
    printf("seed = %d\n",seed);
    printf("phase_flag = %d\n",phase_flag);
    printf("syntype = %d\n",syntype);
  }

  sprintf(tstr,"    Gaussian SD = %.2f grid units\n",sig);
  mylog(mylogf,tstr);
  sprintf(tstr,"    Minimum normalized correlation for connection = %f\n",
	  minw);
  mylog(mylogf,tstr);

  xn1 = lay1->xn;  yn1 = lay1->yn;  zn1 = lay1->zn;
  xn2 = lay2->xn;  yn2 = lay2->yn;  zn2 = lay2->zn;

  cm = get_zero_2d_farray(xn,yn);
  lay1rf = get_5d_farray(xn1,yn1,zn1,-1,-1);
  lay2rf = get_5d_farray(xn2,yn2,zn2,-1,-1);

  ai_ori = pop_cell_attrib_index_f(lay1,"ori"); // Attribute index

  if (pflag) printf("HERE 00\n");

  nsyn = 0;
  for(ai=0;ai<xn1;ai++){
    for(aj=0;aj<yn1;aj++){
      for(ak=0;ak<zn1;ak++){  // For each cell in the pre-syn (first) layer
	
	c = &(lay1->c[ai][aj][ak]);  // Pointer to pre-syn cell

	//
	//  Make the 2D (xn x yn) Gabor template for this cell
	//
	ori = c->attrib_f[ai_ori];
	xc = c->rfx;
	yc = c->rfy;
	//printf("xc,yc = %f %f  i,j %d %d\n",xc,yc,ai,aj);
	//printf("phase = %f  ori=%f  sf=%f corr_sign = %d\n",phase,ori,sf,
	//phase_flag);
	gabor = gabor_2d_space_raw(xn,yn,xc,yc,sdo,sdp,sf,ori,phase);

	if ((ai==0)&&(aj==0)&&(ak==0) && (tmpl_fname != NULL)){
	  write_2d_data(tmpl_fname,gabor,0,0,xn,yn,4,2,1,1);
	}

	lay1rf[ai][aj][ak] = gabor;
      }
    }
  }

  if (pflag) printf("HERE 01\n");

  for(ai=0;ai<xn2;ai++){
    for(aj=0;aj<yn2;aj++){
      for(ak=0;ak<zn2;ak++){  // For each cell in the second layer
	
	// Place connections in 2D farray
	c = &(lay2->c[ai][aj][ak]);

	// WYETH - ASSUMING LGNIN[0]
	if (c->lgn_n > 0){
	  //ln = c->lgnin[0];
	  ln = pop_cell_get_lgn_input_by_name(c,lgn_name);
	  if (ln == NULL){
	    printf("  *** lgn_name:  %s\n",lgn_name);
	    printf("  *** pre_syn layer:  %s\n",lay2->name);
	    exit_error("MOD_CONN_CORR_ON_OFF_TEMPLATE",
		       "No inputs from 'lgn_pname'");
	  }
	}else
	  exit_error("MOD_CONN_CORR_ON_OFF_TEMPLATE","No LGN inputs");


// WYETH - SWITCH TO THIS **** 
//lay2rf[ai][aj][ak] = mod_conn_composite_rf(cm,xn,yn,sig,ln);

	n = ln->cn1;
	x = ln->cx1;
	y = ln->cy1;
	for(i=0;i<n;i++)
	  cm[x[i]][y[i]] += 1.0;  // ON connection
	
	n = ln->cn0;
	x = ln->cx0;
	y = ln->cy0;
	for(i=0;i<n;i++)
	  cm[x[i]][y[i]] -= 1.0;  // OFF connection
	
	lay2rf[ai][aj][ak] = smooth_2d_with_gaussian(cm,xn,yn,sig,0.05);

	// Reset connection matrix to 0's
	n = ln->cn1;
	x = ln->cx1;
	y = ln->cy1;
	for(i=0;i<n;i++)
	  cm[x[i]][y[i]] = 0.0;  // Reset ON connection
	
	n = ln->cn0;
	x = ln->cx0;
	y = ln->cy0;
	for(i=0;i<n;i++)
	  cm[x[i]][y[i]] = 0.0;  // Reset OFF connection

// WYETH - END SWITCH (above)
// WYETH - END SWITCH (above)

      }
    }
  }

  if (pflag){
    printf("WYETH - dumping intermediate 2d data\n");

    write_2d_data("zzz.1.2d",lay1rf[0][0][0],0,0,xn,yn,4,2,0,1);
    write_2d_data("zzz.2.2d",lay2rf[0][0][0],0,0,xn,yn,4,2,0,1);
  }

  if (pflag) printf("HERE 02\n");

  // Store correlation value with each inhib cell
  rfx = get_zero_3d_farray(xn2,yn2,zn2);

  histw = (float *)myalloc(xn2*yn2*zn2*sizeof(float));

  if (prob < 1.0){
    if (seed > 0)
      seed *= -1;
  }

  for(ai=0;ai<xn1;ai++){
    for(aj=0;aj<yn1;aj++){
      for(ak=0;ak<zn1;ak++){  // For each cell in the first layer
	c1 = &(lay1->c[ai][aj][ak]);

	if (pflag) printf("%d %d %d \n",ai,aj,ak);

	for(bi=0;bi<xn2;bi++){
	  for(bj=0;bj<yn2;bj++){
	    for(bk=0;bk<zn2;bk++){  // For each cell in the 2nd layer
	      // Compute correlation of RF masks
	      if (cdist > 0.0){
		c2 = &(lay2->c[bi][bj][bk]);
		dist = mod_conn_get_distance_cells_layer(c1,c2);
		if (dist <= cdist)
		  rfx[bi][bj][bk] = mod_conn_rfx_01(lay1rf[ai][aj][ak],
						    lay2rf[bi][bj][bk],xn,yn);
		else
		  rfx[bi][bj][bk] = 0.0;  // Set weight to zero
	      }else
		rfx[bi][bj][bk] = mod_conn_rfx_01(lay1rf[ai][aj][ak],
						  lay2rf[bi][bj][bk],xn,yn);
	    }
	  }
	}

	mod_conn_corr_w_norm_hist(mylogf,rfx,xn2,yn2,zn2,minw,phase_flag);
	
	// Make connections
	tdelay = 0.0;
	for(bi=0;bi<xn2;bi++){
	  for(bj=0;bj<yn2;bj++){
	    for(bk=0;bk<zn2;bk++){  /* For each cell in the 2nd layer */
	      c2 = &(lay2->c[bi][bj][bk]);
	      w = rfx[bi][bj][bk];
	      if (w > 0.0){
		
		if (axonv > 0.0){
		  dist = mod_conn_get_distance_cells_layer(c1,c2);
		  tdelay = dist * spum;
		}
		
		if (prob < 1.0){
		  if (myrand_util_ran2(&seed) < prob){
		    (void)pop_cell_add_synapse(c2,c1,syntype,w,tdelay,inindex);
		    nsyn += 1;
		  }
		}else{
		  (void)pop_cell_add_synapse(c2,c1,syntype,w,tdelay,inindex);
		  nsyn += 1;
		}
	      }
	    }
	  }
	}
	
	
      }
    }
  }

  myfree(histw);
  free_2d_farray(cm,xn);
  free_3d_farray(rfx,xn2,yn2,zn2);
  free_5d_farray(lay1rf,xn1,yn1,zn1,xn,yn);
  free_5d_farray(lay2rf,xn2,yn2,zn2,xn,yn);

  sprintf(tstr,"    %d synapses from inhibitory to excitatory neurons\n",nsyn);
  mylog(mylogf,tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_CONN_ONODE_CORR_ON_OFF                       */
/*                                                                           */
/*****************************************************************************/
void mod_conn_onode_corr_on_off(mylogf,io,lay1,lay2,syntype,sig,sscale,xn,yn,
				inindex)
     char *mylogf;
     struct onode *io;              // Input node
     struct pop_layer *lay1;        // post-syn layer
     struct pop_layer *lay2;        // pre-syn layer
     int syntype;
     float sig;                     // SD of LGN center
     float sscale;                  // sscale
     int xn,yn;                     // LGN grid
     short inindex;           // Index into population 'inlist'
{
  int seed,csign,alg;
  float minw,prob,axonv,axondt,cdist,normw;
  float sf,sdo,sdp,phase;
  char *tstr,*template,*lgn_name,*tmpl_fname;
  struct onode *disto;   // Distrib node

  disto = onode_child_get_unique(io,"distrib");
  axonv = onode_getpar_flt_dflt(io,"axon_vel",0.0);
  axondt = onode_getpar_flt_dflt(io,"axon_dt",0.0); // (s)

  // Make I-E correlation-based connections from inhib to excit layer
  tstr = onode_getpar_chr_ptr(disto,"file");
  if (tstr != NULL){
    // WYETH THIS IS CHECKED BEFORE CALLING, NO NEED HERE

    mylog_exit(mylogf,"MOD_CONN_ONODE_CORR_ON_OFF  this seems backwards\n");

    mod_conn_read_conn_t1(mylogf,lay1,lay2,tstr,axonv,axondt,inindex);
  }else{
    alg   = onode_getpar_int_dflt(disto,"algorithm",1); // 2008 Sept; was 0
    minw  = onode_getpar_flt_dflt(disto,"minw",0.0);
    cdist = onode_getpar_flt_dflt(disto,"cdist",0.0);
    prob  = onode_getpar_flt_dflt(disto,"prob",1.0);
    seed  = onode_getpar_int_dflt(disto,"seed",1777);
    csign = onode_getpar_int_dflt(disto,"corr_sign",1);
    normw = onode_getpar_flt_dflt(disto,"normw",0.0);

    // NEWLGN - get the population name
    lgn_name = onode_getpar_chr_dflt(disto,"lgn_pname","lgn");

    template = onode_getpar_chr_dflt(disto,"template","self_lgn");

    if (strcmp(template,"Gabor")==0){
      sf    = onode_getpar_flt_exit(disto,"sf");
      sdo   = onode_getpar_flt_exit(disto,"sd_orth");
      sdp   = onode_getpar_flt_exit(disto,"sd_par");
      phase = onode_getpar_flt_exit(disto,"phase");
      tmpl_fname = onode_getpar_chr_null(disto,"write_template");

      //printf("sf = %f\n",sf);
      //printf("sdo = %f\n",sdo);
      //printf("sdp = %f\n",sdp);
      //printf("phase = %f\n",phase);

      sf     *= sscale;  // (cyc/pix)
      sdo    /= sscale;  // (pix)
      sdp    /= sscale;  // (pix)

      mod_conn_corr_on_off_template(mylogf,lay1,lay2,lgn_name,xn,yn,sig/sscale,
				    minw,prob,seed,csign,syntype,cdist,axonv,
				    inindex,sf,sdo,sdp,phase,tmpl_fname);

    }else if (strcmp(template,"self_lgn")==0){

      if (alg == 1){
	// lay1 - post
	// lay2 - pre
	mod_conn_corr_on_off_01(mylogf,lay1,lay2,xn,yn,sig/sscale,minw,
				prob,seed,csign,syntype,cdist,axonv,inindex);
      }else{
	mylog_exit(mylogf,"OLD CORR ALGORITHM, is this intentional?\n");
	mod_conn_corr_on_off(mylogf,lay1,lay2,xn,yn,inindex);
      }
    }else{
      mylog_exit(mylogf,"MOD_CONN_ONODE_CORR_ON_OFF  unknown 'template'\n");
    }
    myfree(template);

    if (normw > 0.0){
      // lay1 is post
      pop_cell_norm_weight_layer_inindex(mylogf,lay1,normw,inindex);
      //pop_cell_norm_weight_layer_layer(mylogf,lay2,lay1,normw,inindex);
    }

    myfree(lgn_name);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_CONN_BINOC_MASK_01                          */
/*                                                                           */
/*  Connect cells from lay2 to those in lay1 on the basis of their ocular    */
/*  dominance and how well their LGN inputs match the mask.                  */
/*                                                                           */
/*  NEWLGN - WYETH - currently, the cells in lay2 are assumed to have        */
/*  LGN inputs that all come from one eye or the other.                      */
/*                                                                           */
/*****************************************************************************/
void mod_conn_binoc_mask_01(mylogf,lay1,lay2,xn,yn,sig,disp_x,disp_y,
			    sf,sdo,sdp,prob,seed,nsamp_l,nsamp_r,syntype,
			    cdist,minw,axonv,inindex)
     char *mylogf;
     struct pop_layer *lay1;  // post-syn
     struct pop_layer *lay2;  // pre-syn
     int xn,yn;
     float sig;     // SD for Gaussian weight of LGN RF center
     float disp_x;  // disparity in x (pix)
     float disp_y;  // disparity in y (pix)
     float sf;      // mask SF (cyc/pix)
     float sdo;     // SD orthogonal to mask (pix)
     float sdp;     // SD parallel to mask (pix)
     float prob;    // Prob. of connection
     int seed;      // Randomization for probabilistic connections
     int nsamp_l;   // Number of connections from left-eye inputs
     int nsamp_r;   // Number of connections from right-eye inputs
     int syntype;   // See ifc.h for definitions
     float cdist;   // Cutoff distance, ignored if <= 0.0
     float minw;    // Don't make connections where weight is less than this
     float axonv;   // Axon velocity, m/s, ignored if <= 0.0
     short inindex;           // Index into population 'inlist'
{
  int i,j,k,ai,aj,ak,bi,bj,bk;
  int n,*x,*y,xn1,yn1,zn1,xn2,yn2,zn2,nsyn,pflag,flag,ai_ocdom,ai_ori1;
  float **cm,tdelay,dist,w,*****lay2rf,***rfx_l,***rfx_r,rv;
  float xc,yc,ori,ph,**maskl,**maskr,ph0;
  double spum;
  char tstr[SLEN];
  struct pop_cell *c,*c1,*c2;
  struct pop_lgn_in *ln;


  mylog(mylogf,"  MOD_CONN_BINOC_MASK_01\n");

  ai_ocdom =  pop_cell_attrib_index_f(lay2,"ocdom"); // Attribute index
  ai_ori1  =  pop_cell_attrib_index_f(lay1,"ori"); // Attribute index

  ph0 = 90.0;
  sprintf(tstr,"    *** Setting mask phase to constant, %.2f\n",ph0);
  mylog(mylogf,tstr);

  pflag = 0;

  if (axonv > 0.0)
    spum = 1.0 / (axonv * 1000000.0);  // Seconds per Micron
  else
    spum = 0.0;

  if (pflag){
    printf("laynames = %s %s\n",lay1->name,lay2->name);
    printf("xn,yn = %d %d\n",xn,yn);
    printf("sig  = %f\n",sig);
    printf("cdist = %f\n",cdist);
    printf("minw = %f\n",minw);
    printf("syntype = %d\n",syntype);
  }

  sprintf(tstr,"    Gaussian SD = %.2f grid units\n",sig);
  mylog(mylogf,tstr);

  xn1 = lay1->xn;  yn1 = lay1->yn;  zn1 = lay1->zn;
  xn2 = lay2->xn;  yn2 = lay2->yn;  zn2 = lay2->zn;

  lay2rf = get_5d_farray(xn2,yn2,zn2,-1,-1);

  if (pflag) printf("HERE 01\n");

  // Fill 'lay2rf' with estimated RFs for all pre-syn cells
  cm = get_zero_2d_farray(xn,yn);
  for(ai=0;ai<xn2;ai++){
    for(aj=0;aj<yn2;aj++){
      for(ak=0;ak<zn2;ak++){  // For each cell in the post-syn layer
	
	// Place connections in 2D farray
	c = &(lay2->c[ai][aj][ak]);

	// WYETH - ASSUMING LGNIN[0]
	if (c->lgn_n > 0)
	  ln = c->lgnin[0];
	else
	  exit_error("MOD_CONN_BINOC_MASK_01","No LGN inputs");


// WYETH - SWITCH TO THIS **** 
//lay2rf[ai][aj][ak] = mod_conn_composite_rf(cm,xn,yn,sig,ln);

	n = ln->cn1;
	x = ln->cx1;
	y = ln->cy1;
	for(i=0;i<n;i++)
	  cm[x[i]][y[i]] += 1.0;  // ON connection
	
	n = ln->cn0;
	x = ln->cx0;
	y = ln->cy0;
	for(i=0;i<n;i++)
	  cm[x[i]][y[i]] -= 1.0;  // OFF connection
	
	lay2rf[ai][aj][ak] = smooth_2d_with_gaussian(cm,xn,yn,sig,0.05);

	// Reset connection matrix to 0's
	n = ln->cn1;
	x = ln->cx1;
	y = ln->cy1;
	for(i=0;i<n;i++)
	  cm[x[i]][y[i]] = 0.0;  // Reset ON connection
	
	n = ln->cn0;
	x = ln->cx0;
	y = ln->cy0;
	for(i=0;i<n;i++)
	  cm[x[i]][y[i]] = 0.0;  // Reset OFF connection

// WYETH - END SWITCH (above)
// WYETH - END SWITCH (above)

      }
    }
  }
  free_2d_farray(cm,xn);

  if (pflag){
    printf("WYETH - dumping intermediate 2d data\n");

    write_2d_data("zzz.2.2d",lay2rf[0][0][0],0,0,xn,yn,4,2,0,1);
  }

  if (pflag) printf("HERE 02\n");

  // Store correlation value with each inhib cell
  rfx_l = get_zero_3d_farray(xn2,yn2,zn2);
  rfx_r = get_zero_3d_farray(xn2,yn2,zn2);

  if (seed > 0)  // For probabilistic connections
    seed *= -1;

  nsyn = 0;
  for(ai=0;ai<xn1;ai++){
    for(aj=0;aj<yn1;aj++){
      for(ak=0;ak<zn1;ak++){  // For each cell in the first layer
	c1 = &(lay1->c[ai][aj][ak]);

	//printf("Making connections for %d %d %d\n",ai,aj,ak);
	
	// Build masks for L and R eyes
	// WYETH ATTRIB
	ori = c1->attrib_f[ai_ori1];
	//ori = c1->ori;
	xc = c1->rfx;
	yc = c1->rfy;
	//printf("xc,yc = %f %f  i,j %d %d\n",xc,yc,ai,aj);
	ph = ph0;
	maskl = gabor_2d_space_raw(xn,yn,xc,yc,sdo,sdp,sf,ori,ph);
	xc += disp_x;
	yc += disp_y;
	maskr = gabor_2d_space_raw(xn,yn,xc,yc,sdo,sdp,sf,ori,ph);

	//printf("sdo %f  sdp %f  sf %f  ori %f  ph %f\n",sdo,sdp,sf,ori,ph);

	if (pflag) printf("%d %d %d \n",ai,aj,ak);

	for(bi=0;bi<xn2;bi++){
	  for(bj=0;bj<yn2;bj++){
	    for(bk=0;bk<zn2;bk++){  // For each cell in the 2nd layer
	      // Compute correlation of RF masks
	      flag = 1;
	      if (cdist > 0.0){
		c2 = &(lay2->c[bi][bj][bk]);

		dist = mod_conn_get_distance_cells_layer(c1,c2);
		if (dist > cdist){
		  flag = 0;        // Do not allow distant connections
		}
	      }

	      rfx_l[bi][bj][bk] = 0.0; // Set weights to zero, in case not
	      rfx_r[bi][bj][bk] = 0.0; //   set below
	      if (flag){

		// WYETH ATTRIB
		if (c2->attrib_f[ai_ocdom] <= 0.0){
		  //if (c2->ocdom <= 0.0){ // Select LEFT mask based on OD
		  rfx_l[bi][bj][bk] = mod_conn_rfx_01(maskl,lay2rf[bi][bj][bk],
						      xn,yn);
		  //printf("HERE left wyeth wyeth %f\n",rfx_l[bi][bj][bk]);
		}else{ // Select RIGHT mask
		  rfx_r[bi][bj][bk] = mod_conn_rfx_01(maskr,lay2rf[bi][bj][bk],
						      xn,yn);
		  //printf("HERE rright wyeth wyeth  %f\n",rfx_r[bi][bj][bk]);
		}
	      }
	    }
	  }
	}
	free_2d_farray(maskl,xn);
	free_2d_farray(maskr,xn);

	
	//
	// Normalize weights so negatives set to zero
	//   and others are scale so 2*SD is w=1, which is also the max
	//
	mod_conn_corr_w_norm_hist(mylogf,rfx_l,xn2,yn2,zn2,minw,1);
	mod_conn_corr_w_norm_hist(mylogf,rfx_r,xn2,yn2,zn2,minw,1);


	// Make connections
	tdelay = 0.0;
	for(bi=0;bi<xn2;bi++){
	  for(bj=0;bj<yn2;bj++){
	    for(bk=0;bk<zn2;bk++){  // For each cell in the 2nd layer
	      c2 = &(lay2->c[bi][bj][bk]);

	      // WYETH ATTRIB
	      if (c2->attrib_f[ai_ocdom] <= 0.0){
		//if (c2->ocdom <= 0.0){
		w = rfx_l[bi][bj][bk];
		//if (bj == yn2/2) printf("  LEFT    w = %f\n",w);
	      }else{
		w = rfx_r[bi][bj][bk];
		//if (bj == yn2/2) printf("    right w = %f\n",w);
	      }

	      if (w > 0.0){

		rv = myrand_util_ran2(&seed);
		if (rv < prob){
		  
		  if (axonv > 0.0){
		    dist = mod_conn_get_distance_cells_layer(c1,c2);
		    tdelay = dist * spum;
		  }
		  (void)pop_cell_add_synapse(c2,c1,syntype,w,tdelay,inindex);
		  nsyn += 1;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  free_3d_farray(rfx_l,xn2,yn2,zn2);
  free_3d_farray(rfx_r,xn2,yn2,zn2);
  free_5d_farray(lay2rf,xn2,yn2,zn2,xn,yn);

  sprintf(tstr,"    %d synapses\n",nsyn);
  mylog(mylogf,tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_CONN_GABOR_MASK_01                          */
/*                                                                           */
/*  Connect cells from lay2 to those in lay1 on the basis of how well        */
/*  their LGN inputs match the mask.                                         */
/*                                                                           */
/*****************************************************************************/
void mod_conn_gabor_mask_01(mylogf,lay1,lay2,xn,yn,sig,
			    sf,sdo,sdp,prob,seed,syntype,
			    cdist,minw,axonv,inindex)
     char *mylogf;
     struct pop_layer *lay1;  // post-syn
     struct pop_layer *lay2;  // pre-syn
     int xn,yn;
     float sig;     // SD for Gaussian weight of LGN RF center
     float sf;      // mask SF (cyc/pix)
     float sdo;     // SD orthogonal to mask (pix)
     float sdp;     // SD parallel to mask (pix)
     float prob;    // Prob. of connection
     int seed;      // Randomization for probabilistic connections
     int syntype;   // See ifc.h for definitions
     float cdist;   // Cutoff distance, ignored if <= 0.0
     float minw;    // Don't make connections where weight is less than this
     float axonv;   // Axon velocity, m/s, ignored if <= 0.0
     short inindex;           // Index into population 'inlist'
{
  int i,j,k,ai,aj,ak,bi,bj,bk;
  int n,*x,*y,xn1,yn1,zn1,xn2,yn2,zn2,nsyn,pflag,flag,ai_ori1;
  float **cm,tdelay,dist,w,*****lay2rf,***rfx_l,rv;
  float xc,yc,ori,ph,**maskl,ph0;
  double spum;
  char tstr[SLEN];
  struct pop_cell *c,*c1,*c2;
  struct pop_lgn_in *ln;

  mylog(mylogf,"  MOD_CONN_GABOR_MASK_01\n");

  ai_ori1  =  pop_cell_attrib_index_f(lay1,"ori"); // Attribute index

  ph0 = 90.0;
  sprintf(tstr,"    *** Setting mask phase to constant, %.2f\n",ph0);
  mylog(mylogf,tstr);

  pflag = 0;

  if (axonv > 0.0)
    spum = 1.0 / (axonv * 1000000.0);  // Seconds per Micron
  else
    spum = 0.0;

  if (pflag){
    printf("laynames = %s %s\n",lay1->name,lay2->name);
    printf("xn,yn = %d %d\n",xn,yn);
    printf("sig  = %f\n",sig);
    printf("cdist = %f\n",cdist);
    printf("minw = %f\n",minw);
    printf("syntype = %d\n",syntype);
  }

  sprintf(tstr,"    Gaussian SD = %.2f grid units\n",sig);
  mylog(mylogf,tstr);

  xn1 = lay1->xn;  yn1 = lay1->yn;  zn1 = lay1->zn;
  xn2 = lay2->xn;  yn2 = lay2->yn;  zn2 = lay2->zn;

  lay2rf = get_5d_farray(xn2,yn2,zn2,-1,-1);

  if (pflag) printf("HERE 01\n");

  // Fill 'lay2rf' with estimated RFs for all pre-syn cells
  cm = get_zero_2d_farray(xn,yn);
  for(ai=0;ai<xn2;ai++){
    for(aj=0;aj<yn2;aj++){
      for(ak=0;ak<zn2;ak++){  // For each cell in the post-syn layer
	
	// Place connections in 2D farray
	c = &(lay2->c[ai][aj][ak]);

	// WYETH - ASSUMING LGNIN[0]
	if (c->lgn_n > 0)
	  ln = c->lgnin[0];
	else
	  exit_error("MOD_CONN_GABOR_MASK_01","No LGN inputs");


// WYETH - SWITCH TO THIS **** 
//lay2rf[ai][aj][ak] = mod_conn_composite_rf(cm,xn,yn,sig,ln);

	n = ln->cn1;
	x = ln->cx1;
	y = ln->cy1;
	for(i=0;i<n;i++)
	  cm[x[i]][y[i]] += 1.0;  // ON connection
	
	n = ln->cn0;
	x = ln->cx0;
	y = ln->cy0;
	for(i=0;i<n;i++)
	  cm[x[i]][y[i]] -= 1.0;  // OFF connection
	
	lay2rf[ai][aj][ak] = smooth_2d_with_gaussian(cm,xn,yn,sig,0.05);

	// Reset connection matrix to 0's
	n = ln->cn1;
	x = ln->cx1;
	y = ln->cy1;
	for(i=0;i<n;i++)
	  cm[x[i]][y[i]] = 0.0;  // Reset ON connection
	
	n = ln->cn0;
	x = ln->cx0;
	y = ln->cy0;
	for(i=0;i<n;i++)
	  cm[x[i]][y[i]] = 0.0;  // Reset OFF connection

// WYETH - END SWITCH (above)
// WYETH - END SWITCH (above)

      }
    }
  }
  free_2d_farray(cm,xn);

  if (pflag){
    printf("WYETH - dumping intermediate 2d data\n");

    write_2d_data("zzz.2.2d",lay2rf[0][0][0],0,0,xn,yn,4,2,0,1);
  }

  if (pflag) printf("HERE 02\n");

  // Store correlation value with each inhib cell
  rfx_l = get_zero_3d_farray(xn2,yn2,zn2);

  if (seed > 0)  // For probabilistic connections
    seed *= -1;

  nsyn = 0;
  for(ai=0;ai<xn1;ai++){
    for(aj=0;aj<yn1;aj++){
      for(ak=0;ak<zn1;ak++){  // For each cell in the first layer
	c1 = &(lay1->c[ai][aj][ak]);

	//printf("Making connections for %d %d %d\n",ai,aj,ak);
	
	// Build masks for L and R eyes
	// WYETH ATTRIB
	ori = c1->attrib_f[ai_ori1];
	//ori = c1->ori;
	xc = c1->rfx;
	yc = c1->rfy;
	//printf("xc,yc = %f %f  i,j %d %d\n",xc,yc,ai,aj);
	ph = ph0;
	maskl = gabor_2d_space_raw(xn,yn,xc,yc,sdo,sdp,sf,ori,ph);

	//printf("sdo %f  sdp %f  sf %f  ori %f  ph %f\n",sdo,sdp,sf,ori,ph);

	if (pflag) printf("%d %d %d \n",ai,aj,ak);

	for(bi=0;bi<xn2;bi++){
	  for(bj=0;bj<yn2;bj++){
	    for(bk=0;bk<zn2;bk++){  // For each cell in the 2nd layer
	      // Compute correlation of RF masks
	      flag = 1;
	      if (cdist > 0.0){
		c2 = &(lay2->c[bi][bj][bk]);

		dist = mod_conn_get_distance_cells_layer(c1,c2);
		if (dist > cdist){
		  flag = 0;        // Do not allow distant connections
		}
	      }

	      rfx_l[bi][bj][bk] = 0.0; // Set weights to zero, in case not
	      if (flag){
		rfx_l[bi][bj][bk] = mod_conn_rfx_01(maskl,lay2rf[bi][bj][bk],
						    xn,yn);
		//printf("HERE left wyeth wyeth %f\n",rfx_l[bi][bj][bk]);
	      }
	    }
	  }
	}
	free_2d_farray(maskl,xn);

	
	//
	// Normalize weights so negatives set to zero
	//   and others are scale so 2*SD is w=1, which is also the max
	//
	mod_conn_corr_w_norm_hist(mylogf,rfx_l,xn2,yn2,zn2,minw,1);


	// Make connections
	tdelay = 0.0;
	for(bi=0;bi<xn2;bi++){
	  for(bj=0;bj<yn2;bj++){
	    for(bk=0;bk<zn2;bk++){  // For each cell in the 2nd layer
	      c2 = &(lay2->c[bi][bj][bk]);

	      w = rfx_l[bi][bj][bk];
	      if (w > 0.0){

		rv = myrand_util_ran2(&seed);
		if (rv < prob){
		  
		  if (axonv > 0.0){
		    dist = mod_conn_get_distance_cells_layer(c1,c2);
		    tdelay = dist * spum;
		  }
		  (void)pop_cell_add_synapse(c2,c1,syntype,w,tdelay,inindex);
		  nsyn += 1;
		}
	      }
	    }
	  }
	}
      }
    }
  }

  free_3d_farray(rfx_l,xn2,yn2,zn2);
  free_5d_farray(lay2rf,xn2,yn2,zn2,xn,yn);

  sprintf(tstr,"    %d synapses\n",nsyn);
  mylog(mylogf,tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_CONN_ONODE_BINOC_MASK                        */
/*                                                                           */
/*****************************************************************************/
void mod_conn_onode_binoc_mask(mylogf,io,lay1,lay2,syntype,sig,sscale,xn,yn,
			       binoc_flag,inindex)
     char *mylogf;
     struct onode *io;              // Input node
     struct pop_layer *lay1;        // post-syn
     struct pop_layer *lay2;        // pre-syn
     int syntype;
     float sig;                     // SD of LGN center
     float sscale;                  // sscale
     int xn,yn;                     // LGN grid
     int binoc_flag;                // 0-monocular, 1-binocular
     int inindex;
{
  int nsamp_l,nsamp_r,seed;
  float axonv,axondt,disp_x,disp_y,sf,sdo,sdp,cdist,minw,prob,normw;
  char *tstr;
  struct onode *disto;   // Distrib node

  disto = onode_child_get_unique(io,"distrib");
  axonv  = onode_getpar_flt_dflt(io,"axon_vel",0.0);
  axondt = onode_getpar_flt_dflt(io,"axon_dt",0.0); // (s)

  // Make I-E correlation-based connections from inhib to excit layer
  tstr = onode_getpar_chr_ptr(disto,"file");
  if (tstr != NULL){
    // WYETH THIS IS CHECKED BEFORE CALLING, NO NEED HERE
    mod_conn_read_conn_t1(mylogf,lay1,lay2,tstr,axonv,axondt,inindex);
  }else{

    if (binoc_flag == 1){
      disp_x = onode_getpar_flt_dflt(disto,"disp_x",0.0);
      disp_y = onode_getpar_flt_dflt(disto,"disp_y",0.0);
    }

    sf    = onode_getpar_flt_exit(disto,"sf");
    sdo   = onode_getpar_flt_exit(disto,"sd_orth");
    sdp   = onode_getpar_flt_exit(disto,"sd_par");
    cdist = onode_getpar_flt_exit(disto,"cdist");
    minw  = onode_getpar_flt_exit(disto,"minw");
    seed  = onode_getpar_int_exit(disto,"seed");
    prob  = onode_getpar_flt_exit(disto,"prob");
    normw = onode_getpar_flt_exit(disto,"normw");

    //nsamp_l  = onode_getpar_int_exit(disto,"nsamp_L");
    //nsamp_r  = onode_getpar_int_exit(disto,"nsamp_R");

    disp_x /= sscale; // Convert parameters to pixels
    disp_y /= sscale;
    sf     *= sscale;
    sdo    /= sscale;
    sdp    /= sscale;

    if (binoc_flag == 1){
      mod_conn_binoc_mask_01(mylogf,lay1,lay2,xn,yn,sig/sscale,
			     disp_x,disp_y,sf,sdo,sdp,prob,seed,nsamp_l,
			     nsamp_r,syntype,cdist,minw,axonv,inindex);
    }else{
      mod_conn_gabor_mask_01(mylogf,lay1,lay2,xn,yn,sig/sscale,
			     sf,sdo,sdp,prob,seed,syntype,cdist,minw,axonv,
			     inindex);
    }

    if (normw > 0.0){
      pop_cell_norm_weight_layer_inindex(mylogf,lay1,normw,inindex);

      //  WYETH - BUG fixed - THESE LAYERS WERE BACKWARDS, fixed 05/Oct/2008
      //pop_cell_norm_weight_layer_layer(mylogf,lay2,lay1,normw);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                      MOD_CONN_LAYER_TO_CELL_ORI_DIST_01                   */
/*                                                                           */
/*  Connect cells in lay1 to cell c2 based on distance and ori criterion.    */
/*                                                                           */
/*****************************************************************************/
void mod_conn_layer_to_cell_ori_dist_01(mylogf,lay1,c2,cdist,cori,minw,wf,
					prob,seed,self,syntype,inindex)
     char *mylogf;
     struct pop_layer *lay1;
     struct pop_cell *c2;
     float cdist;        // Critical distance to limit connections (um)
     float cori;         // Critical orientation difference (degr)
     float minw;         // gaus < minw is ignored, max is 1.0
     float wf;           // multiply weight by this factor
     float prob;         // Probability of making connection
     int seed;           // Randomization seed for random connections
     int self;           // 0-Do not connect to self
     int syntype;        // type of synapse
     short inindex;      // Index into population 'inlist'
{
  int ai,aj,ak;
  int xn1,yn1,zn1,nsyn,ai_ori1;
  float w,dist,dori,rv,c1_ori,c2_ori;
  char tstr[SLEN];
  struct pop_cell *pc1;

  mylog(mylogf,"  MOD_CONN_LAYER_TO_CELL_ORI_DIST_01\n");

  if (seed > 0)
    seed *= -1;

  sprintf(tstr,"    Connecting layer '%s' to cell '%s'\n",lay1->name,
	  c2->name);
  mylog(mylogf,tstr);
  sprintf(tstr,"      Critical distance:  %.4f um\n",cdist);
  mylog(mylogf,tstr);
  sprintf(tstr,"      Critical ori diff:  %.4f degrees\n",cori);
  mylog(mylogf,tstr);

  xn1 = lay1->xn;  yn1 = lay1->yn;  zn1 = lay1->zn;

  c2_ori = pop_cell_attrib_get_f(c2,"ori");
  ai_ori1 =  pop_cell_attrib_index_f(lay1,"ori"); // Attribute index

  nsyn = 0;
  for(ai=0;ai<xn1;ai++){
    for(aj=0;aj<yn1;aj++){
      for(ak=0;ak<zn1;ak++){  // For each cell in the first layer

	rv = myrand_util_ran2(&seed);
	if (rv < prob){
	  
	  pc1 = &(lay1->c[ai][aj][ak]);
	  
	  if ((c2 == pc1) && (self == 0)){
	    ; // Ignoring connections to self
	  }else{

	    c1_ori = pc1->attrib_f[ai_ori1];
	    dori = mod_conn_get_ori_diff(c1_ori,c2_ori,1);

	    // Weight is product of Gaussian in distance and ori
	    w = func_gaussian_one(dori,0.0,cori);
	    dist = mod_conn_get_distance_cells_layer(pc1,c2);
	    w *= func_gaussian_one(dist,0.0,cdist);
	    
	    if (w >= minw){
	      (void)pop_cell_add_synapse(pc1,c2,syntype,wf*w,0.0,inindex);
	      nsyn += 1;
	    }
	  }
	}
      }
    }
  }
  sprintf(tstr,"    %d synapses formed\n",nsyn);
  mylog(mylogf,tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_CONN_ORI_DIST_01                          */
/*                                                                           */
/*  Connect cells in lay1 to cells in lay2 based on distance and ori         */
/*  criterion.                                                               */
/*                                                                           */
/*****************************************************************************/
void mod_conn_ori_dist_01(mylogf,lay1,lay2,cdist,cori,cori_muflag,cori_val,
			  minw,wf,prob,seed,self,syntype,axonv,inindex,shape)
     char *mylogf;
     struct pop_layer *lay1;  // pre-syn
     struct pop_layer *lay2;  // post-syn
     float cdist;        // Critical distance to limit connections (um)
     float cori;         // Critical orientation difference (degr)
     int cori_muflag;    // 0-self, 1-offset, 2-fixed
     float cori_val;     // depends on cori_muflag, offset or fixed value (deg)
     float minw;         // gaus < minw is ignored, max is 1.0
     float wf;           // multiply weight by this factor
     float prob;         // Probability of making connection
     int seed;           // Randomization seed for random connections
     int self;           // 0-Do not connect to self
     int syntype;        // type of synapse
     float axonv;        // Axon AP velocity
     short inindex;      // Index into population 'inlist'
     char *shape;        // 'Gaussian', 'uniform', ...
{
  int ai,aj,ak,bi,bj,bk;
  int xn1,yn1,zn1,xn2,yn2,zn2,nsyn,ai_ori;
  float w,dist,dori,rv,tdelay,tori;
  double spum;
  char tstr[SLEN];
  struct pop_cell *pc1,*pc2;

  mylog(mylogf,"  MOD_CONN_ORI_DIST_01\n");

  ai_ori =  pop_cell_attrib_index_f(lay2,"ori"); // Attribute index

  if (seed > 0)
    seed *= -1;

  sprintf(tstr,"    Connecting layer '%s' to layer '%s'\n",lay1->name,
	  lay2->name);

  mylog(mylogf,tstr);
  sprintf(tstr,"      Critical distance:  %.4f um\n",cdist);
  mylog(mylogf,tstr);
  sprintf(tstr,"      Critical ori diff:  %.4f degrees\n",cori);
  mylog(mylogf,tstr);

  if (axonv > 0.0)
    spum = 1.0 / (axonv * 1000000.0);  // Seconds per Micron
  else
    spum = 0.0;

  xn1 = lay1->xn;  yn1 = lay1->yn;  zn1 = lay1->zn;
  xn2 = lay2->xn;  yn2 = lay2->yn;  zn2 = lay2->zn;

  tdelay = 0.0;
  nsyn = 0;
  for(ai=0;ai<xn1;ai++){
    for(aj=0;aj<yn1;aj++){
      for(ak=0;ak<zn1;ak++){  // For each cell in the pre-syn layer

	pc1 = &(lay1->c[ai][aj][ak]);

	for(bi=0;bi<xn2;bi++){
	  for(bj=0;bj<yn2;bj++){
	    for(bk=0;bk<zn2;bk++){  // For each cell in the post-syn layer

	      rv = myrand_util_ran2(&seed);
	      if (rv < prob){
		
		pc2 = &(lay2->c[bi][bj][bk]);  // Post-syn

		if ((pc2 == pc1) && (self == 0)){
		  ; // printf("  Ignoring self\n");
		}else{
		  // Weight is product of Gaussian in distance and ori

		  if (cori < 0.0){
		    w = 1.0;  // No constraint on orientation
		  }else{
		    if (cori_muflag == 0){
		      //tori = pc2->attrib_f[ai_ori]; // WYETH - UNUSED
		      dori = mod_conn_get_ori_diff_cells(pc1,pc2,1);
		    }else if (cori_muflag == 1){
		      tori = pc2->attrib_f[ai_ori] + cori_val;
		      dori = mod_conn_get_ori_diff_one(pc1,tori,1);
		    }else if (cori_muflag == 2){
		      tori = cori_val;
		      dori = mod_conn_get_ori_diff_one(pc1,cori_val,1);
		    }
		    if (strcmp(shape,"Gaussian")==0){
		      w = func_gaussian_one(dori,0.0,cori);
		    }else if (strcmp(shape,"uniform")==0){
		      if (dori <= cori)
			w = 1.0;
		      else
			w = 0.0;
		    }
		  }

		  if (cdist < 0.0){
		    w *= 1.0;  // No constraint on distance
		  }else{
		    dist = mod_conn_get_distance_cells_layer(pc1,pc2);
		    if (axonv > 0.0){
		      tdelay = dist * spum;
		    }
		    if (strcmp(shape,"Gaussian")==0){
		      w *= func_gaussian_one(dist,0.0,cdist);
		    }else if (strcmp(shape,"uniform")==0){
		      if (dist <= cdist)
			w *= 1.0;
		      else
			w *= 0.0;
		    }
		  }

		  if ((w >= minw) && (w > 0.0)){ // WYETH 2010 MAR
		    (void)pop_cell_add_synapse(pc1,pc2,syntype,wf*w,tdelay,
					       inindex);
		    nsyn += 1;
		  }
		}
	      }


	    }
	  }
	}
      }
    }
  }
  sprintf(tstr,"    %d synapses formed\n",nsyn);
  mylog(mylogf,tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_CONN_ONODE_ORI_DIST_01                       */
/*                                                                           */
/*  Connect cells in lay1 to cells in lay2 based on distance and ori         */
/*  criterion.                                                               */
/*                                                                           */
/*****************************************************************************/
void mod_conn_onode_ori_dist_01(mylogf,io,lay1,lay2,syntype,inindex)
     char *mylogf;
     struct onode *io;              // Input node
     struct pop_layer *lay1,*lay2;  // pre- and post-syn layers
     int syntype;
     short inindex;                 // Index into population 'inlist'
{
  float cdist;        // Critical distance to limit connections (um)
  float cori;         // Critical orientation difference (degr)
  float minw;         // gaus < minw is ignored, max is 1.0
  float wf;           // multiply weight by this factor
  float prob;         // Probability of making connection
  float normw;        // Probability of making connection
  int seed;           // Randomization seed for random connections
  int self;           // 0-Do not connect to self
  float axonv;        // Velocity of spike in axon
  int cori_muflag;    // 0-self, 1-offset, 2-fixed 
  float cori_val;     // could be offset or mean, depending on 'cori_muflag'
  char *shape;        // 'Gaussian' or 'uniform' [Gaussian]

  struct onode *disto;   // Distrib node

  axonv = onode_getpar_flt_dflt(io,"axon_vel",0.0);
  disto = onode_child_get_unique(io,"distrib");
  cdist = onode_getpar_flt_dflt(disto,"cdist",-1.0); // WYETH was 0.0 Sep2008
  cori  = onode_getpar_flt_dflt(disto,"cori",-1.0);  // WYETH was 10.0 Sep2008

  cori_muflag = onode_getpar_int_dflt(disto,"ori_mean_flag",0);
  if (cori_muflag == 1){
    cori_val = onode_getpar_flt_exit(disto,"cori_offset");
  }else if (cori_muflag == 2){
    cori_val = onode_getpar_flt_exit(disto,"cori_mu");
  }

  minw  = onode_getpar_flt_dflt(disto,"minw",0.0);
  normw = onode_getpar_flt_dflt(disto,"normw",0.0);
  seed  = onode_getpar_int_dflt(disto,"seed",23138);
  prob  = onode_getpar_flt_dflt(disto,"prob",1.0);
  self  = onode_getpar_int_dflt(disto,"self",1);
  shape = onode_getpar_chr_dflt(disto,"shape","Gaussian");

  wf = 1.0;
  mod_conn_ori_dist_01(mylogf,lay1,lay2,cdist,cori,cori_muflag,cori_val,
		       minw,wf,prob,seed,self,syntype,axonv,inindex,shape);

  if (normw > 0.0)
    pop_cell_norm_weight_layer_inindex(mylogf,lay2,normw,inindex);
  //pop_cell_norm_weight_layer_layer(mylogf,lay1,lay2,normw);
  
}
/**************************************-**************************************/
/*                                                                           */
/*                        MOD_CONN_CELL_TO_LAYER_GAUSS                       */
/*                                                                           */
/*  Connect cell to cells in lay2 with Gaussian weight.                      */
/*                                                                           */
/*****************************************************************************/
void mod_conn_cell_to_layer_gauss(mylogf,c1,lay2,sd,minw,wf,self,syntype,
				  inindex)
     char *mylogf;
     struct pop_cell *c1;
     struct pop_layer *lay2;
     float sd;                   // SD of Gaussian, cortical distance
     float minw;                 // gaus < minw is ignored, max is 1.0
     float wf;                   // multiply weight by this factor
     int self;                   // 0-Do not connect to self
     int syntype;
     short inindex;              // Index into population 'inlist'
{
  int bi,bj,bk;
  int xn2,yn2,zn2,nsyn;
  float w,dist;
  char tstr[SLEN];
  struct pop_cell *pc2;
  int printflag;

  printflag = 1;
  if (printflag){
    mylog(mylogf,"  MOD_CONN_CELL_TO_LAYER_GAUSS\n");

    sprintf(tstr,"    Connecting cell '%s' to layer '%s'\n",c1->name,
	    lay2->name);
    mylog(mylogf,tstr);
    sprintf(tstr,"      Gaussian SD:  %.4f um\n",sd);
    mylog(mylogf,tstr);
    sprintf(tstr,"      Minimum weight:  %.4f (max is 1.0)\n",minw);
    mylog(mylogf,tstr);
    sprintf(tstr,"      Weight scale:  %.4f\n",wf);
    mylog(mylogf,tstr);
  }

  xn2 = lay2->xn;  yn2 = lay2->yn;  zn2 = lay2->zn;

  nsyn = 0;
  for(bi=0;bi<xn2;bi++){
    for(bj=0;bj<yn2;bj++){
      for(bk=0;bk<zn2;bk++){  /* For each cell in the 2nd layer */
	
	pc2 = &(lay2->c[bi][bj][bk]);
	
	if ((pc2 == c1) && (self == 0)){
	  ; /*printf("  Ignoring self\n");*/
	}else{
	  /* Determine distance between cells */
	  dist = mod_conn_get_distance_cells_layer(c1,pc2);
	  w = func_gaussian_one(dist,0.0,sd);
	  if (w >= minw){
	    (void)pop_cell_add_synapse(c1,pc2,syntype,wf*w,0.0,inindex);
	    nsyn += 1;
	  }
	}
	
      }
    }
  }

  if (printflag){
    sprintf(tstr,"    %d synapses formed\n",nsyn);
    mylog(mylogf,tstr);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                        MOD_CONN_LAYER_TO_CELL_GAUSS                       */
/*                                                                           */
/*  Connect cells in lay1 to cell c2 with Gaussian weight.                   */
/*                                                                           */
/*****************************************************************************/
void mod_conn_layer_to_cell_gauss(mylogf,lay1,c2,sd,minw,wf,self,syntype,
				  inindex)
     char *mylogf;
     struct pop_layer *lay1;
     struct pop_cell *c2;
     float sd;                   // SD of Gaussian, cortical distance
     float minw;                 // gaus < minw is ignored, max is 1.0
     float wf;                   // multiply weight by this factor
     int self;                   // 0-Do not connect to self
     int syntype;
     short inindex;              // Index into population 'inlist'
{
  int ai,aj,ak;
  int xn1,yn1,zn1,nsyn;
  float w,dist;
  char tstr[SLEN];
  struct pop_cell *pc1;
  int printflag;

  printflag = 1;

  if (printflag){
    mylog(mylogf,"  MOD_CONN_LAYER_TO_CELL_GAUSS\n");
    
    sprintf(tstr,"    Connecting layer '%s' to cell '%s'\n",lay1->name,
	    c2->name);
    mylog(mylogf,tstr);
    sprintf(tstr,"      Gaussian SD:  %.4f um\n",sd);
    mylog(mylogf,tstr);
    sprintf(tstr,"      Minimum weight:  %.4f (max is 1.0)\n",minw);
    mylog(mylogf,tstr);
    sprintf(tstr,"      Weight scale:  %.4f\n",wf);
    mylog(mylogf,tstr);
  }

  xn1 = lay1->xn;  yn1 = lay1->yn;  zn1 = lay1->zn;

  nsyn = 0;
  for(ai=0;ai<xn1;ai++){
    for(aj=0;aj<yn1;aj++){
      for(ak=0;ak<zn1;ak++){  /* For each cell in the first layer */
	
	pc1 = &(lay1->c[ai][aj][ak]);

	if ((c2 == pc1) && (self == 0)){
	  ; /*printf("  Ignoring self\n");*/
	}else{
	  /* Determine distance between cells */
	  dist = mod_conn_get_distance_cells_layer(pc1,c2);
	  w = func_gaussian_one(dist,0.0,sd);
	  if (w >= minw){
	    (void)pop_cell_add_synapse(pc1,c2,syntype,wf*w,0.0,inindex);
	    nsyn += 1;
	  }
	}
	
      }
    }
  }
  if (printflag){
    sprintf(tstr,"    %d synapses formed\n",nsyn);
    mylog(mylogf,tstr);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_CONN_LAYERS_GAUSS                           */
/*                                                                           */
/*  Connect cells in lay1 to cells in lay2 with Gaussian weight.             */
/*                                                                           */
/*****************************************************************************/
void mod_conn_layers_gauss(mylogf,lay1,lay2,sd,minw,wf,self,syntype,inindex)
     char *mylogf;
     struct pop_layer *lay1,*lay2;
     float sd;                   // SD of Gaussian, cortical distance
     float minw;                 // gaus < minw is ignored, max is 1.0
     float wf;                   // multiply weight by this factor
     int self;                   // 0-Do not connect to self
     int syntype;
     short inindex;              // Index into population 'inlist'
{
  int ai,aj,ak,bi,bj,bk;
  int xn1,yn1,zn1,xn2,yn2,zn2,nsyn;
  float w,dist;
  char tstr[SLEN];
  struct pop_cell *pc1,*pc2;

  mylog(mylogf,"  MOD_CONN_LAYERS_GAUSS\n");

  sprintf(tstr,"    Connecting layer '%s' to layer '%s'\n",lay1->name,
	  lay2->name);
  mylog(mylogf,tstr);
  sprintf(tstr,"      Gaussian SD:  %.4f um\n",sd);
  mylog(mylogf,tstr);
  sprintf(tstr,"      Minimum weight:  %.4f (max is 1.0)\n",minw);
  mylog(mylogf,tstr);
  sprintf(tstr,"      Weight scale:  %.4f\n",wf);
  mylog(mylogf,tstr);

  xn1 = lay1->xn;  yn1 = lay1->yn;  zn1 = lay1->zn;
  xn2 = lay2->xn;  yn2 = lay2->yn;  zn2 = lay2->zn;

  nsyn = 0;
  for(ai=0;ai<xn1;ai++){
    for(aj=0;aj<yn1;aj++){
      for(ak=0;ak<zn1;ak++){  // For each cell in the first layer

	pc1 = &(lay1->c[ai][aj][ak]);

	for(bi=0;bi<xn2;bi++){
	  for(bj=0;bj<yn2;bj++){
	    for(bk=0;bk<zn2;bk++){  // For each cell in the 2nd layer

	      pc2 = &(lay2->c[bi][bj][bk]);

	      if ((pc2 == pc1) && (self == 0)){
		; // printf("  Ignoring self\n");
	      }else{
		// Determine distance between cells
		dist = mod_conn_get_distance_cells_layer(pc1,pc2);
		w = func_gaussian_one(dist,0.0,sd);
		if (w >= minw){
		  (void)pop_cell_add_synapse(pc1,pc2,syntype,wf*w,0.0,inindex);
		  nsyn += 1;
		}
	      }

	    }
	  }
	}

      }
    }
  }
  sprintf(tstr,"    %d synapses formed\n",nsyn);
  mylog(mylogf,tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_CONN_IMAGE_2D_LGN_GRID                        */
/*                                                                           */
/*****************************************************************************/
float **mod_conn_image_2d_lgn_grid(c,xa,ya,xb,yb,rw,rh)
     struct pop_cell *c;  // Cell
     int xa,ya,xb,yb;     // Corners of grid sub-region to process
     int *rw,*rh;         // Width and Height of returned 2d image
{
  int i,j,k;
  int n,w,h,nw,nh,x1,y1,x2,y2,x,y,nn,wflag,hflag;
  float **d;
  struct pop_lgn_in *ln;

  nn = 5;  // Width of grid box + 1 for margin

  if (xa <= xb){
    x1 = xa;
    x2 = xb;
  }else{
    x1 = xb;
    x2 = xa;
  }
  
  if (ya <= yb){
    y1 = ya;
    y2 = yb;
  }else{
    y1 = yb;
    y2 = ya;
  }

  nw = x2 - x1 + 1;  // Width of region in grid points
  nh = y2 - y1 + 1;  // Height of region in grid points

  w = nn*nw + 1;      // Width of image (pix)
  h = nn*nh + 1;      // Height of image (pix)

  // Keep w and h even
  wflag = hflag = 0;
  if (w%2 == 1){
    w += 1;
    wflag = 1;
  }
  if (h%2 == 1){
    h += 1;
    hflag = 1;
  }

  d = get_const_2d_farray(w,h,0.5);  // Fill w/ gray

  for(i=0;i<w;i+=nn){
    for(j=0;j<h;j++){
      d[i][j] = 0.0;  // Black vertical lines
    }
  }
  for(j=0;j<h;j+=nn){
    for(i=0;i<w;i++){
      d[i][j] = 0.0;  // Black horizontal lines
    }
  }

  // For some reason, postplot seems to want even .2d file
  if (hflag == 1){  // White edges, padding for even
    for(i=0;i<w;i++)
      d[i][h-1] = 1.0;
  }
  if (wflag == 1){  // White edges, padding for even
    for(j=0;j<h;j++)
      d[w-1][j] = 1.0;
  }


  // WYETH - ASSUMING LGNIN[0]
  if (c->lgn_n > 0)
    ln = c->lgnin[0];
  else
    exit_error("MOD_CONN_IMAGE_2D_LGN_GRID","No LGN inputs");

  // ON
  n = ln->cn1;
  for(i=0;i<n;i++){
    x = ln->cx1[i] - x1;
    y = ln->cy1[i] - y1;
    if ((x >= 0) && (x < nw) && (y >= 0) && (y < nh)){ // Point in region
      for(j=1;j<nn;j++){
	for(k=1;k<nn;k++){
	  d[x*nn+j][y*nn+k] = 1.0;
	}
      }
    }
  }

  // OFF
  n = ln->cn0;
  for(i=0;i<n;i++){
    x = ln->cx0[i] - x1;
    y = ln->cy0[i] - y1;
    if ((x >= 0) && (x < nw) && (y >= 0) && (y < nh)){ // Point in region
      for(j=1;j<nn;j++){
	for(k=1;k<nn;k++){
	  d[x*nn+j][y*nn+k] = 0.0;
	}
      }
    }
  }

  *rw = w;
  *rh = h;
  return d;
}
