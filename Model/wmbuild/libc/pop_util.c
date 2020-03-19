/*****************************************************************************/
/*                                                                           */
/*  pop_util.c                                                               */
/*  wyeth bair                                                               */
/*                                                                           */
/*  Solve population of cells simultaneously.                                */
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
#include "spike_util.h"
#include "paramfile_util.h"
#include "sig_util.h"
#include "mod_util.h"
#include "pop_cell_util.h"
#include "mod_conn_util.h"
//#include "mod_mon_util.h"

#include "mod.h" // Data structures
#include "ifc.h" // IFC and POP data structures
#include "paramfile.h"

#define HMAX 0.0005       // Maximum step for Runge-Kutta (sec)

int global_dump_flag = 0;
int global_dump_flag2 = 0;
int globseed;

/*************************************---*************************************/
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
/**************************************-**************************************/
/*                                                                           */
/*                                POP_UTIL_RAN2                              */
/*                                                                           */
/*  To be used within one routine, i.e., when no other routine could inter-  */
/*  fer by trying to reseed.                                                 */
/*                                                                           */
/*  COPIED FROM NumRecInC, Second Edition, page 282.                         */
/*                                                                           */
/*  Long period (>2*10^18) random number generator of L'Ecuyer with          */
/*  Bays-Durham shuffle and added safeguards.  Returns a uniform deviate     */
/*  between 0.0 and 1.0 (exclusive of the endpoint values).  Call with idum  */
/*  a negative integer to initialize; thereafter, do not alter idum between  */
/*  successive deviates in a sequence.  RNMX should approximate the largest  */
/*  floating value that is less than 1.                                      */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float pop_util_ran2(idum)
     int *idum; // Changed from long for consistency on SUNS, DEC ALPHA
{
  int j;
  int k; // Changed from long to int.
  static int idum2=123456789; /*** Changed from long to int. ***/
  static int iy=0; /*** Changed from long to int. ***/
  static int iv[NTAB]; /*** Changed from long to int. ***/
  float temp;

  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/**************************************-**************************************/
/*                                                                           */
/*                               POP_UTIL_GASDEV                             */
/*                                                                           */
/*  Returns a normally distributed deviate with zero mean and unit           */
/*  variance, using ran1(idum) as the source of uniform deviates.            */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float pop_util_gasdev(idum)
     int *idum;
{
  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;
  
  if ((iset == 0) || (*idum <= 0)){  /* Wyeth, reset state if reseeded */
    do {
      v1=2.0*pop_util_ran2(idum)-1.0;
      v2=2.0*pop_util_ran2(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  }else{
    iset=0;
    return gset;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         POP_UTIL_GET_MASK_DIFF_EXP                        */
/*                                                                           */
/*  Return a float array of length '*rnm' to be determined that contains a   */
/*  difference of exponentials mask.  The time units of the returned array   */
/*  are those of the time constants, presumably milliseconds.                */
/*                                                                           */
/*****************************************************************************/
float *pop_util_get_mask_diff_exp(m,name,rnm)
     struct model_struct *m;          // Model
     //struct param_pair_list *mppl;  // Model parameters
     char name[];                     // Name prefix for mask parameters
     int *rnm;                        // Return length of mask
{
  int nm;
  float t1,t2,amp,*mask;
  char tstr[SLEN];

  /// WYETH HERE
  /// WYETH HERE
  /// WYETH HERE  fixing NMDA
  /// WYETH HERE
  /// WYETH HERE

  sprintf(tstr,"%s_tau_f",name);
  t1 = paramfile_get_float_param_or_exit(m->ppl,tstr);  // "simp_ex_tau_f"
  sprintf(tstr,"%s_tau_r",name);
  t2 = paramfile_get_float_param_or_exit(m->ppl,tstr);
  sprintf(tstr,"%s_amp",name);
  amp = paramfile_get_float_param_or_exit(m->ppl,tstr);

  nm = (int)(t1 * 7.0);        // Length of mask is 7 times decay time
  mask = diff_exp_farray(0.0,amp,t1,t2,nm);

  *rnm = nm;
  return mask;
}
/**************************************-**************************************/
/*                                                                           */
/*                      POP_UTIL_GET_MASK_DIFF_EXP_MECH                      */
/*                                                                           */
/*  Return a float array of length '*rnm' to be determined that contains a   */
/*  difference of exponentials mask.  The time units of the returned array   */
/*  are those of the time constants, presumably milliseconds.                */
/*                                                                           */
/*****************************************************************************/
float *pop_util_get_mask_diff_exp_mech(mech,rnm)
     struct pop_mech *mech;  // Mechanism
     int *rnm;               // Return length of mask
{
  int nm;
  float *mask;

  nm = (int)(mech->tau_f * 7.0);    // Length of mask is 7 times decay time
  mask = diff_exp_farray(0.0,mech->amp,mech->tau_f,mech->tau_r,nm);

  *rnm = nm;
  return mask;
}
/**************************************-**************************************/
/*                                                                           */
/*                        POPU_MECH_ONODE_GET_MASK_ALPHA                     */
/*                                                                           */
/*  Return the alpha function mask given the mechanism onode.                */
/*                                                                           */
/*****************************************************************************/
float *popu_mech_onode_get_mask_alpha(o,sampling,rn)
     struct onode *o;    // Mechanism onode
     float sampling;     // samples/s (1000.0 for ms)
     int *rn;            // Return length of array
{
  int maskn;
  float tau,amp,alpha,*mask;
  char *mtype;

  mtype = onode_getpar_chr_exit(o,"type");
  tau   = onode_getpar_flt_exit(o,"tau");
  amp   = onode_getpar_flt_exit(o,"amp");

  if (strcmp(mtype,"psg_alpha")!=0)
    exit_error("POPU_MECH_ONODE_GET_MASK_ALPHA","Wrong mask type");

  maskn = (int)(tau * sampling * 10.0);

  alpha = 1.0/(tau * sampling); // alpha is 1/peaktime
  mask = alpha_farray(0.0,alpha,1.0,maskn);
  make_max_const_farray(mask,maskn,amp);
  
  myfree(mtype);

  *rn = maskn;
  return mask;
}
/**************************************-**************************************/
/*                                                                           */
/*                       POP_UTIL_ONODE_GET_POP_XN_YN_ZN                     */
/*                                                                           */
/*  Return a pointer to the layer with the given name, or NULL.              */
/*                                                                           */
/*****************************************************************************/
void pop_util_onode_get_pop_xn_yn_zn(o,name,rxn,ryn,rzn)
     struct onode *o;    // Top node
     char *name;
     int *rxn,*ryn,*rzn;
{
  int axn,ayn,azn;
  char *aname;
  struct onode *po,*geo,*ao;

  po = onode_get_node_type_item_val(o,"pop","name",name);
  if (po == NULL){
    *rxn = *ryn = *rzn = -1;
    return;
  }
  geo = onode_child_get_unique(po,"geometry");

  if (onode_test_ostr(po,"area")){
    //
    //  If an 'area' is referenced, use it for default 'xn' and 'yn'
    //
    aname = onode_getpar_chr_exit(po,"area");
    ao = onode_get_node_type_item_val(o,"area","name",aname);
    myfree(aname);

    axn = onode_getpar_int_exit(ao,"xn");
    ayn = onode_getpar_int_exit(ao,"yn");
    azn = onode_getpar_int_dflt(ao,"zn",1);

    *rxn = onode_getpar_int_dflt(geo,"xn",axn);
    *ryn = onode_getpar_int_dflt(geo,"yn",ayn);

  }else{
    *rxn = onode_getpar_int_exit(geo,"xn");
    *ryn = onode_getpar_int_exit(geo,"yn");
    azn = 1;
  }

  *rzn = onode_getpar_int_dflt(geo,"zn",azn);
}
/**************************************-**************************************/
/*                                                                           */
/*                          POP_UTIL_GET_AREA_BY_NAME                        */
/*                                                                           */
/*****************************************************************************/
struct pop_area *pop_util_get_area_by_name(mpt,name)
     struct pop_top *mpt;
     char *name;
{
  int i;
  struct pop_area *pa;
  
  pa = NULL;
  for(i=0;i<mpt->narea;i++){
    if (strcmp(mpt->area[i]->name,name)==0){
      if (pa == NULL)
	pa = mpt->area[i];
      else
	exit_error("POP_UTIL_GET_AREA_BY_NAME","Duplicate area names");
    }
  }

  return pa;
}
/**************************************-**************************************/
/*                                                                           */
/*                          POP_UTIL_GET_LAYER_BY_NAME                       */
/*                                                                           */
/*  New for mpt, replaces pop_util_get_layer_for_name above.                 */
/*                                                                           */
/*  Returns NULL if layer not found.  Exits if duplicate names.              */
/*                                                                           */
/*****************************************************************************/
struct pop_layer *pop_util_get_layer_by_name(mpt,name)
     struct pop_top *mpt;
     char *name;
{
  int i;
  struct pop_layer *pl;
  
  pl = NULL;
  for(i=0;i<mpt->nlay;i++){
    if (strcmp(mpt->lay[i]->name,name)==0){
      if (pl == NULL)
	pl = mpt->lay[i];
      else
	exit_error("POP_UTIL_GET_LAYER_BY_NAME","Duplicate population names");
    }
  }

  return pl;
}
/**************************************-**************************************/
/*                                                                           */
/*                         POPU_GET_LAYER_INDEX_BY_NAME                      */
/*                                                                           */
/*  Returns -1 if layer not found.  Exits if duplicate names.                */
/*                                                                           */
/*****************************************************************************/
int popu_get_layer_index_by_name(mpt,name)
     struct pop_top *mpt;
     char *name;
{
  int i,k;
  
  k = -1;
  for(i=0;i<mpt->nlay;i++){
    if (strcmp(mpt->lay[i]->name,name)==0){
      if (k == -1)
	k = i;
      else
	exit_error("POPU_GET_LAYER_INDEX_BY_NAME","Duplicate population name");
    }
  }

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                             POPU_GET_LAY_NAMED                            */
/*                                                                           */
/*  Get the layer named by the value of 'parname' in the object 'o'.         */
/*  Exit if not found.                                                       */
/*                                                                           */
/*****************************************************************************/
struct pop_layer *popu_get_lay_named(mpt,o,parname,callstr)
     struct pop_top *mpt;   // Model structure
     struct onode *o;       // Node in which "parname" should appear
     char *parname;         // Parameter that specifies layer name
     char *callstr;         // String to print on failure and exit
{
  char *lay_name;
  struct pop_layer *lay;

  // Get pointers to two layers of orgin
  lay_name = onode_getpar_chr_exit(o,parname);
  lay = pop_util_get_layer_by_name(mpt,lay_name);
  if (lay == NULL){
    mylog_cval(mpt->logf,"*** Parameter Name",parname);
    mylog_cval(mpt->logf,"*** Caller String",callstr);
    mylog_cval(mpt->logf,"*** Layer Name",lay_name);
    mylog_exit(mpt->logf,"POPU_GET_LAY_NAMED  Layer not found.");
  }
  myfree(lay_name);
  
  return lay;
}
/**************************************-**************************************/
/*                                                                           */
/*                            POPU_GET_MECH_NAMED                            */
/*                                                                           */
/*  Get the mechanism named by the value of 'parname' in the object 'o'.     */
/*  Exit if not found.                                                       */
/*                                                                           */
/*****************************************************************************/
struct pop_mech *popu_get_mech_named(mpt,lay,o,parname,callstr)
     struct pop_top *mpt;   // Model structure
     struct pop_layer *lay; // Layer
     struct onode *o;       // Node in which "parname" should appear
     char *parname;         // Parameter that specifies mechanism name
     char *callstr;         // String to print on failure and exit
{
  char *mech_name;
  struct pop_mech *mech;

  // Get pointers to two layers of orgin
  mech_name = onode_getpar_chr_exit(o,parname);
  mech = pop_cell_get_mech_by_name(lay,mech_name);
  if (mech == NULL){
    mylog_cval(mpt->logf,"*** Parameter Name",parname);
    mylog_cval(mpt->logf,"*** Caller String",callstr);
    mylog_cval(mpt->logf,"*** Mech Name",mech_name);
    mylog_exit(mpt->logf,"POPU_GET_MECH_NAMED  Named mechanism not found.");
  }
  myfree(mech_name);
  
  return mech;
}
/**************************************-**************************************/
/*                                                                           */
/*                              POPU_GET_C_NAMED                             */
/*                                                                           */
/*  Get the cell named by the string in 'parname' in the object 'o'.         */
/*  Exit if not found.                                                       */
/*                                                                           */
/*  The parameter VALUE, which is the cell name, is of the form:             */
/*        <layer>_<x>_<y>_<z>,  e.g.,  ex_9_11_0                             */
/*                                                                           */
/*****************************************************************************/
struct pop_cell *popu_get_c_named(mpt,o,parname,callstr)
     struct pop_top *mpt;   // Model structure
     struct onode *o;       // Node in which "parname" should appear
     char *parname;         // Parameter that specifies layer name
     char *callstr;         // String to print on failure and exit
{
  int x,y,z;
  char *lay_name,*cell_name;
  struct pop_layer *lay;
  struct pop_cell *c;

  cell_name = onode_getpar_chr_exit(o,parname);

  pop_cell_parse_unit_name(cell_name,&lay_name,&x,&y,&z);

  lay = pop_util_get_layer_by_name(mpt,lay_name);
  if (lay == NULL){
    mylog_cval(mpt->logf,"*** Parameter Name",parname);
    mylog_cval(mpt->logf,"*** Caller String",callstr);
    mylog_cval(mpt->logf,"*** Layer Name",lay_name);
    mylog_exit(mpt->logf,"POPU_GET_C_NAMED  Layer not found.");
  }

  c = &(lay->c[x][y][z]);
  if (c == NULL){
    mylog_cval(mpt->logf,"*** cell name",cell_name);
    mylog_exit(mpt->logf,"POPU_GET_C_NAMED  cell not found");
  }
    
  myfree(lay_name);
  myfree(cell_name);

  return c;
}
/**************************************-**************************************/
/*                                                                           */
/*                        POP_UTIL_GET_XN_YN_ZN_BY_NAME                      */
/*                                                                           */
/*  Could also be called, get_pop_dimensions_by_name                         */
/*                                                                           */
/*****************************************************************************/
void pop_util_get_xn_yn_zn_by_name(mpt,name,rxn,ryn,rzn)
     struct pop_top *mpt;
     char *name;
     int *rxn,*ryn,*rzn;
{
  struct pop_layer *pl;

  pl = pop_util_get_layer_by_name(mpt,name);
  
  *rxn = pl->xn;
  *ryn = pl->yn;
  *rzn = pl->zn;
}
/**************************************-**************************************/
/*                                                                           */
/*                           POP_UTIL_GET_MAP_BY_NAME                        */
/*                                                                           */
/*  New for mpt, replaces pop_util_get_layer_for_name above.                 */
/*                                                                           */
/*****************************************************************************/
struct pop_map *pop_util_get_map_by_name(mpt,name)
     struct pop_top *mpt;
     char *name;
{
  int i;
  struct pop_map *pm;
  char *tname;
  
  pm = NULL;
  for(i=0;i<mpt->nmap;i++){
    tname = onode_getpar_chr_ptr_exit(mpt->map[i]->o,"name");
    if (strcmp(tname,name)==0){
      if (pm == NULL)
	pm = mpt->map[i];
      else
	exit_error("POP_UTIL_GET_MAP_BY_NAME","Duplicate map names");
    }
  }

  return pm;
}
/**************************************-**************************************/
/*                                                                           */
/*                             POP_UTIL_MAP_MAKE                             */
/*                                                                           */
/*  Make the 'k'th map.                                                      */
/*                                                                           */
/*****************************************************************************/
void pop_util_map_make(mpt,k)
     struct pop_top *mpt;        // Global model structure
     int k;
{
  char *mylogf,ggstr[SLEN],*aname,*mtype;
  struct onode *mo;       // Map onode
  struct pop_map *mp;     // Map structure
  struct pop_area *area;

  mylogf = mpt->logf;
  mylog(mylogf,"  POP_UTIL_MAP_MAKE\n");

  if ((k < 0) || (k >= mpt->nmap))
    mylogx(mylogf,"POP_UTIL_MAP_MAKE","Map index out of bounds");

  mo = mpt->map[k]->o;
  mp = mpt->map[k];

  if (onode_item(mo,"area")==1){
    aname = onode_getpar_chr_ptr_exit(mo,"area");
    area = pop_util_get_area_by_name(mpt,aname); // NULL if not found

    if (area == NULL){
      sprintf(ggstr,"Could not find 'area' %s for map.\n",aname);
      mylogx(mylogf,ggstr);
    }

    mtype = onode_getpar_chr_ptr_exit(mo,"type");

    if (strcmp(mtype,"ori")==0){
      mod_conn_onode_make_map_ori(mylogf,mp,area->xn,area->yn);
    }else if (strcmp(mtype,"ocdom")==0){
      mod_conn_onode_make_map_od(mylogf,mp,area->xn,area->yn);
    }else if (strcmp(mtype,"sf")==0){
      mod_conn_onode_make_map_sf(mylogf,mp,area->xn,area->yn);
    }else if (strcmp(mtype,"rfx")==0){
      mod_conn_onode_make_map_rfx(mylogf,mp,area->xn,area->yn,area->xf);
    }else{
      sprintf(ggstr,"POP_UTIL_MAP_MAKE:  Unknown map type,  %s\n",mtype);
      mylogx(mylogf,ggstr);
    }
  }else{
    sprintf(ggstr,"*** WARNING: Map should contain 'area' parameter.\n");
    mylog(mylogf,ggstr);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              POPU_WARP_APPLY                              */
/*                                                                           */
/*****************************************************************************/
float **popu_warp_apply(map,xn,yn,wx,wy)
     float **map;      // [nx][ny] original and warped maps
     int xn,yn;        // Size of map
     float **wx,**wy;  // [nx][ny] coordinates of source
{
  int i,j;
  int xi,yi;
  float **wmap;

  wmap = get_2d_farray(xn,yn);

  for (i=0;i<xn;i++){
    for(j=0;j<yn;j++){

      xi = my_rint(wx[i][j]);      // Round source coords to nearest int
      yi = my_rint(wy[i][j]);

      while(xi < 0)	xi += xn;  // Wrap around in [xn,yn] in case the
      while(xi >= xn)   xi -= xn;  //   warp goes outside the map
      while(yi < 0)     yi += yn;
      while(yi >= yn)   yi -= yn;

      wmap[i][j] = map[xi][yi];
    }
  }

  return wmap;
}
/**************************************-**************************************/
/*                                                                           */
/*                              POPU_WARP_REMAP                              */
/*                                                                           */
/*  Return the transformation of 'x' via the function 'y'.                   */
/*  This is largely just a linear interpolation routine for 'y'.             */
/*                                                                           */
/*****************************************************************************/
float popu_warp_remap(y,n,x,dx)
     float *y;    // [n] Mapping array
     int n;       // Length of mapping array
     float x;     // Value to be remapped
     float dx;    // x-axis units per mapping array index
{
  int i;
  float v,r;

  i = (int)(x/dx);   // Index into mapping array
  if (x < 0.0){
    v = y[0]; // Shouldn't happen, but if so, use lowest value in mapping array
  }else if (i >= (n-1)){  // We are at or beyond the end, so extrapolate
    v = y[n-1]*x/((n-1)*dx);
  }else{
    r = x - ((float)i*dx);  // Remainder
    v = y[i] + (r/dx)*(y[i+1]-y[i]);  // Linear interpolation
  }
  return v;
}
/**************************************-**************************************/
/*                                                                           */
/*                            POPU_WARP_MAP_COORDS                           */
/*                                                                           */
/*  Given a set of warps, perfor                                             */
/*                                                                           */
/*****************************************************************************/
void popu_warp_map_coords(xn,yn,wpar,nwarps,fns,n,xfact,rwx,rwy)
     int xn,yn;       // Size of 2D map to be warped
     float **wpar;    // [nwarps] Parameters for each warp
     int nwarps;      // Number of warps to apply
     float **fns;     // [nwarps][n] Arrays for remapping radial distance
     int n;           // Length of mapping arrays
     float xfact;     // x-axis units per map array index
     float ***rwx,***rwy;  // [xn][yn] x and y source coordinates
{
  int i,j,w;
  int xi,yi;
  float cx,cy,dx,dy,rr,x,y,r,theta,**wx,**wy;

  wx = get_2d_farray(xn,yn);
  wy = get_2d_farray(xn,yn);

  for (i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      x = (float)i;
      y = (float)j;

      for(w=0;w<nwarps;w++){  // Transform x,y once for each warp

	cx = wpar[w][4];  // Coordinate for center of this warp
	cy = wpar[w][5];

	dx = x - cx;  // Compute polar coorinates of (x,y)
	dy = y - cy;
	r = sqrt(dx*dx + dy*dy);
	theta = atan2f(dy,dx);

        rr = popu_warp_remap(fns[w],n,r,xfact); // Remap radial distance 'r'

	x = cx + rr * cos(theta); // Convert back to Cartesian coords.
	y = cy + rr * sin(theta);
      }
      wx[i][j] = x; // Save the floating point source to assign to (i,j)
      wy[i][j] = y;
    }
  }
  *rwx = wx;  // Return the coordinates needed to perform the warping
  *rwy = wy;
}
/**************************************-**************************************/
/*                                                                           */
/*                          POPU_WARP_GET_MAP_ARRAYS                         */
/*                                                                           */
/*****************************************************************************/
float **popu_warp_get_map_arrays(nwarps,n,params,xfact,dumpfile)
     int nwarps;        // Number of warps
     int n;             // Length of each mapping array
     float **params;    // Warp parameters [nwarps][6]
     float xfact;       // x-axis units per map array index
     char *dumpfile;    // If not NULL, name of file to write map arrays
{
  int i,j;
  float a,sf,ph,sigma,**fns,x,ff,*t;
  char tstr[SLEN];

  if (dumpfile != NULL)
    remove_file(dumpfile);

  t = (float *)myalloc(n*sizeof(float));
  fns = get_2d_farray(nwarps,n);
  for (i=0;i<nwarps;i++){
    a     = params[i][0];
    sf    = params[i][1];   // cyc/unit
    ph    = params[i][2];   // Radians
    sigma = params[i][3];

    ff = sf * 2.0 * M_PI;   // Convert SF from cyc/unit to rad/unit

    fns[i][0] = 0.0;
    for (j=1;j<n;j++){
      x = (float)j * xfact;
      t[j] = a*cos(ff*x + ph) * exp(-x*x/(2.0*sigma*sigma));
      fns[i][j] = fns[i][j-1] + 1.0 + t[j];
    }

    if (dumpfile != NULL){
      sprintf(tstr,"warp_%d",i);
      append_farray_plot("zz.fmap.pl",tstr,fns[i],n,1);
      sprintf(tstr,"gabor_%d",i);
      append_farray_plot("zz.fmap.pl",tstr,t,n,1);
    }
  }

  myfree(t);

  return fns;
}
/**************************************-**************************************/
/*                                                                           */
/*                            POPU_WARP_GEN_PARAMS                           */
/*                                                                           */
/*****************************************************************************/
float **popu_warp_gen_params(mylogf,wo,xn,yn,rnw,rnpar,rfmlen,rxfact)
     char *mylogf;       // Log file
     struct onode *wo;   // Warp onode
     int xn,yn;          // Size of array to be warped
     int *rnw;           // Return the number of warps
     int *rnpar;         // Return the number of parameters per warp
     int *rfmlen;        // Length of mapping arrays
     float *rxfact;      // Scale factor, x-units per mapping array index
{
  int i;
  int npar,nw,seed;
  int xctype;  // 0-contraction, 1-expansion, 2-alternate, 3-random
  float a,sigma,ph,sf,cx,cy,**params,marg,a_max,a_min,sd_min,sd_max;
  char ggstr[SLEN];

  //  Read params from <warp> object
  nw     = onode_getpar_int_exit(wo,"n");     // Number of warps
  seed   = onode_getpar_int_exit(wo,"seed");  // Randomization seed
  a_min  = onode_getpar_flt_exit(wo,"ampl_min");
  a_max  = onode_getpar_flt_exit(wo,"ampl_max");
  sd_min = onode_getpar_flt_exit(wo,"sd_min");
  sd_max = onode_getpar_flt_exit(wo,"sd_max");
  marg   = onode_getpar_int_exit(wo,"margin");  // Avoid warps at map edges
  xctype = onode_getpar_int_exit(wo,"expcon_type");

  //  Set other constants
  *rfmlen = (int)(1 + sqrt(xn*xn + yn*yn));  // Length of remapping arrays
  *rxfact = 1.0;

  npar = 6;
  params = get_2d_farray(nw,npar);

  if (seed > 0)
    seed = -seed;

  sprintf(ggstr,"                a      sf    ph     sigma     cx      cy\n");
  mylog(mylogf,ggstr);

  for (i=0;i<nw;i++){
    a     =  a_min + ( a_max -  a_min)*myrand_util_ran2(&seed);
    sigma = sd_min + (sd_max - sd_min)*myrand_util_ran2(&seed);

    if (xctype == 0)        // Contraction only
      ph = 0.0;
    else if (xctype == 1)   // Expansion only
      ph = M_PI;
    else if (xctype == 2){  // Alternate expansion and contraction
      if (i%2 == 0)
	ph = M_PI;
      else
	ph = 0.0;
    }else if (xctype == 3){  // Randomly choose expansion or contraction
      if (myrand_util_ran2(&seed) < 0.5){
	ph = 0.0;
      }else{
	ph = M_PI;
      }
    }

    cx = marg + (float)(xn-2*marg)*myrand_util_ran2(&seed);
    cy = marg + (float)(yn-2*marg)*myrand_util_ran2(&seed);

    // WYETH - this '3.2' could be made an external param
    sf = 1.0 / (3.2 * sigma);  // Scale SF inversely with Gaussian SD

    params[i][0] = a;
    params[i][1] = sf;
    params[i][2] = ph;
    params[i][3] = sigma;
    params[i][4] = cx;
    params[i][5] = cy;

    sprintf(ggstr,"    Warp[%2d]  %.4f %.4f %.2f %9.4f %7.2f %7.2f\n",i,
	    params[i][0],params[i][1],params[i][2],params[i][3],params[i][4],
	    params[i][5]);
    mylog(mylogf,ggstr);
  }

  *rnw = nw;
  *rnpar = npar;
  return params;
}
/**************************************-**************************************/
/*                                                                           */
/*                            POPU_WARP_MAP_SIMPLE                           */
/*                                                                           */
/*  Return a simple map with lines in it for testing the warping.            */
/*                                                                           */
/*****************************************************************************/
float **popu_warp_map_simple(nx,ny,linesep)
     int nx,ny;    // Size of returned 2D map
     int linesep;  // Separation between lines (pix)
{
  int i,j;
  float **map;

  map = get_2d_farray(ny,nx);

  for (i=0;i<nx;i++){
    for (j=0;j<ny;j++){
      if ((i % linesep == 0) || (j % linesep == 0))
	map[j][i] = 1.0;
      else
	map[j][i] = 0.0;
    }
  }
  return map;
}
/**************************************-**************************************/
/*                                                                           */
/*                               POPU_WARP_MAP                               */
/*                                                                           */
/*  Apply a warping to the maps.                                             */
/*                                                                           */
/*  *** WYETH ---- for now, ASSUME ALL MAPS WILL BE WARPED                   */
/*  *** WYETH ---- for now, ASSUME ALL MAPS WILL BE WARPED                   */
/*                                                                           */
/*****************************************************************************/
void popu_warp_map(mpt,wo)
     struct pop_top *mpt;    // Global model structure
     struct onode *wo;       // Warp onode
{
  int i;
  int xn,yn,nw,fmlen,npar;
  float xfact,**wpar,**fmap,**wx,**wy,**wmap;
  char *mylogf,ggstr[SLEN],*dumpfile,*txtfile;
  struct pop_map *mp;     // Map structure

  mylogf = mpt->logf;
  mylog(mylogf,"  POPU_WARP_MAP\n");

  if (mpt->nmap <= 0)
    return;

  //  Get size from first map
  xn = mpt->map[0]->xn;
  yn = mpt->map[0]->yn;

  //  Get a simple map for testing
  //map = map_simple(xn,yn,8);

  //  Choose parameters for each warp
  wpar = popu_warp_gen_params(mylogf,wo,xn,yn,&nw,&npar,&fmlen,&xfact);
  
  //  Create a mapping array for each warp
  dumpfile = NULL;
  fmap = popu_warp_get_map_arrays(nw,fmlen,wpar,xfact,dumpfile);

  //  Get the overall coordinates that perform the set of warps
  popu_warp_map_coords(xn,yn,wpar,nw,fmap,fmlen,xfact,&wx,&wy);


  //
  //  Use the warping coords to perform the warp on all maps
  //
  for(i=0;i<mpt->nmap;i++){

    mp = mpt->map[i];  // Pointer to the 'i'th map

    if ((mp->xn != xn) || (mp->yn != yn))
      mylogx(mylogf,"POPU_WARP_MAP  All maps are not the same size");

    wmap = popu_warp_apply(mp->map2,xn,yn,wx,wy); // Get a warped map

    free_2d_farray(mp->map2,xn);  // Free the old map and replace it
    mp->map2 = wmap;

    sprintf(ggstr,"    Warping map %d  xn,yn  %d %d\n",i,xn,yn);
    mylog(mylogf,ggstr);

    txtfile = onode_getpar_chr_dflt(mp->o,"write_as_text_warp",NULL);
    if (txtfile != NULL){
      mod_conn_write_map_text(txtfile,mp->map2,xn,yn,1);  // 1-exit flag
      myfree(txtfile);
    }
  }


  free_2d_farray(fmap,nw);
  free_2d_farray(wpar,nw);
  free_2d_farray(wx,xn);
  free_2d_farray(wy,xn);
}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_UTIL_CHECK_LCLASS                          */
/*                                                                           */
/*  Return 1 if 'ltype' is a type that belongs to class 'lclass'.            */
/*                                                                           */
/*****************************************************************************/
int pop_util_check_lclass(ltype,lclass)
     char *ltype;
     char *lclass;
{
  int flag;

  flag = 0;

  if (strcmp(lclass,"lgn0")==0){
    //if ((strcmp(ltype,"lgn")==0) || (strcmp(ltype,"retina0_gc0")==0))
    if (strcmp(ltype,"lgn")==0)
      flag = 1;
  }else if (strcmp(lclass,"rgci")==0){
    // NOTE - not all 'retina0_rgc0' will be 'rgci' WYETH - IS THIS OK?
    if (strcmp(ltype,"retina0_rgc0")==0)
      flag = 1;
  }else{
    exit_error("POP_UTIL_CHECK_LCLASS","Unknown layer class");
  }

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                           POP_UTIL_LCLASS_IS_LGN0                         */
/*                                                                           */
/*  Return 1 if the layer is of the class LGN0, which has the following      */
/*  properties:                                                              */
/*    - Storage for cells is created only for those returning output.        */
/*    - Having zn=2 implies index 0 has OFF and 1 has ON cells               */
/*    - Having zn=4 implies index 0/1 are OFF/ON for left eye, 2/3 right.    */
/*                                                                           */
/*  NOTES                                                                    */
/*  - A layer class may contain many layer types.                            */
/*  - Layer class is currently assigned in the moo files.                    */
/*                                                                           */
/*****************************************************************************/
int pop_util_lclass_is_lgn0(lay)
     struct pop_layer *lay;
{
  int flag;
  char *lt;

  lt = lay->laytype;
  if (lt == NULL)
    exit_error("POP_UTIL_LCLASS_IS_LGN0","Null layer type");

  flag = pop_util_check_lclass(lt,"lgn0");

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                           POP_UTIL_LCLASS_IS_RGCI                         */
/*                                                                           */
/*  Return 1 if the layer is of the class RGCI, meaning it is an irregular   */
/*  array of units that maps to RGCs in the mesh model.                      */
/*    - Storage is created for all cells.                                    */
/*    - Responses are precomputed before the cortical layers run.            */
/*    - runflag = 0                                                          */
/*                                                                           */
/*  NOTES                                                                    */
/*  - A layer class may contain many layer types.                            */
/*  - Layer class is currently assigned in the moo files.                    */
/*                                                                           */
/*****************************************************************************/
int pop_util_lclass_is_rgci(lay)
     struct pop_layer *lay;
{
  int flag;
  char *lt;

  flag = 0;

  lt = lay->laytype;
  if (lt == NULL)
    exit_error("POP_UTIL_LCLASS_IS_RGCI","Null layer type");

  if (strcmp(lt,"retina0_gc0")==0){
    if (lay->geomt == NULL)
      exit_error("POP_UTIL_LCLASS_IS_RGCI","geomt is NULL");
    else
      if (strcmp(lay->geomt,"irregular")==0){
	flag = 1;
      }
  }
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                            POPU_FREE_POP_CELL                             */
/*                                                                           */
/*****************************************************************************/
void popu_free_pop_cell(c,logf)
     struct pop_cell *c;
     char *logf;
{
  if (c->name     != NULL)  myfree(c->name);
  if (c->subclass != NULL)  myfree(c->subclass);

  if (c->attrib_f != NULL)  myfree(c->attrib_f);

  if ((c->nout > 0) || (c->nin > 0)){
    pop_cell_remove_all_syn_from_cell(logf,c);

    //printf(" nout = %d\n",c->nout);
    //exit_error("POPU_FREE_POP_CELL","Free not impl'd yet for:  out");
  }
 //myfree(c);
}
/**************************************-**************************************/
/*                                                                           */
/*                              POPU_FREE_POP_LAYER                          */
/*                                                                           */
/*****************************************************************************/
void popu_free_pop_layer(t,logf)
     struct pop_layer *t;
     char *logf;
{
  int i,j,k;
  int xn,yn,zn;
  struct pop_cell *tc;

  xn = t->xn;
  yn = t->yn;
  zn = t->zn;

  //printf("HERE a  NAME = %s\n",t->name);

  myfree(t->name);
  myfree(t->laytype);

  if (t->icon != NULL){
    exit_error("POPU_FREE_POP_LAYER","Free not impl'd yet for:  icon");
  }

  if (t->geomt != NULL)  myfree(t->geomt);

  if (t->irr_n > 0){
    myfree(t->irr_x);
    myfree(t->irr_y);
    myfree(t->irr_id);
  }

  if (t->layseed_n > 0)
    myfree(t->layseed);

  for(i=0;i<t->nmech;i++){
    exit_error("POPU_FREE_POP_LAYER","Free not impl'd yet for:  mech");
  }

  //printf("HERE f  xn,yn,zn %d %d %d\n",xn,yn,zn);

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	tc = &(t->c[i][j][k]);
	popu_free_pop_cell(tc,logf);
      }
      myfree(t->c[i][j]);
    }
    myfree(t->c[i]);
  }
  if (t->c != NULL)  // WYETH - do we need to chck this??
    myfree(t->c);


  if (t->attrib_fn > 0){
    free_2d_carray(t->attrib_fname,t->attrib_fn);
  }

  if (t->ifcp != NULL)
    exit_error("POPU_FREE_POP_LAYER","Free not impl'd yet for:  ifcp");

  if (t->poissp != NULL)
    exit_error("POPU_FREE_POP_LAYER","Free not impl'd yet for:  poissp");

  if (t->nmda_p != NULL)
    exit_error("POPU_FREE_POP_LAYER","Free not impl'd yet for:  nmda_p");

  if (t->csav_n > 0)
    exit_error("POPU_FREE_POP_LAYER","Free not impl'd yet for:  csav_");

  if (t->cnt != NULL)
    exit_error("POPU_FREE_POP_LAYER","Free not impl'd yet for:  lgn...");

  myfree(t);
}  
/**************************************-**************************************/
/*                                                                           */
/*                             POPU_FREE_POP_AREA                            */
/*                                                                           */
/*****************************************************************************/
void popu_free_pop_area(t,logf)
     struct pop_area *t;
     char *logf;
{
  myfree(t->name);
  myfree(t);
}  
/**************************************-**************************************/
/*                                                                           */
/*                              POPU_FREE_POP_MAP                            */
/*                                                                           */
/*****************************************************************************/
void popu_free_pop_map(t,logf)
     struct pop_map *t;
     char *logf;
{
  // Assuming 't->o' is a pointer.  Is this correct?

  if (t->map2 != NULL){
    free_2d_farray(t->map2,t->xn);
  }    

  myfree(t);
}  
/**************************************-**************************************/
/*                                                                           */
/*                              POPU_FREE_POP_TOP                            */
/*                                                                           */
/*****************************************************************************/
void popu_free_pop_top(t)
     struct pop_top *t;
{
  int i;
  int n,xn,yn,tn;
  char *logf;

  if (t->logf != NULL)
    logf = strdup(t->logf);
  else
    logf = NULL;

  mylog(logf,"  POPU_FREE_POP_TOP\n");

  if (t->name       != NULL)  myfree(t->name);
  if (t->logf       != NULL)  myfree(t->logf);
  if (t->out_prefix != NULL)  myfree(t->out_prefix);

  xn = t->xn;
  yn = t->yn;
  tn = t->tn;

  n = t->narea;
  for(i=0;i<n;i++)
    popu_free_pop_area(t->area[i],logf);


  n = t->nlay;
  for(i=0;i<n;i++){
    //printf("    Layer %d\n",i);
    popu_free_pop_layer(t->lay[i],logf);
  }

  n = t->nmap;
  for(i=0;i<n;i++)
    popu_free_pop_map(t->map[i],logf);

  myfree(t);

  free(logf);
}
/**************************************-**************************************/
/*                                                                           */
/*                            POPU_MAKE_INIT_CELLS                           */
/*                                                                           */
/*****************************************************************************/
struct pop_cell ***popu_make_init_cells(xn,yn,zn,pl)
     int xn,yn,zn;          // Dimensions of layer
     struct pop_layer *pl;  // Layer pointer
{
  int i,j,k;
  struct pop_cell ***c,*cp;

  c = (struct pop_cell ***)myalloc(xn*sizeof(struct pop_cell **));
  for(i=0;i<xn;i++){
    c[i] = (struct pop_cell **)myalloc(yn*sizeof(struct pop_cell *));
    for(j=0;j<yn;j++){
      c[i][j] = (struct pop_cell *)myalloc(zn*sizeof(struct pop_cell));
      for(k=0;k<zn;k++){
	cp = &(c[i][j][k]);

	pop_cell_init_cell(cp);  // Set zero and NULL values

	cp->rfx  = (float)i;
	cp->rfy  = (float)j;

	cp->layx = i;
	cp->layy = j;
	cp->layz = k;

	cp->pl   = pl;  // WYETH - Added Sep 5, 2015

	// Position in tissue, (0,0) is center
	//cp->cx = cx[i];
	//cp->cy = cy[i];
      }
    }
  }

  return c;
}
/**************************************-**************************************/
/*                                                                           */
/*                             POPU_MAKE_POP_LAYER                           */
/*                                                                           */
/*****************************************************************************/
struct pop_layer *popu_make_pop_layer(pname)
     char *pname;
{
  struct pop_layer *pl;

  pl = (struct pop_layer *)myalloc(sizeof(struct pop_layer));

  if (pname == NULL)
    pl->name = NULL;
  else
    pl->name = strdup(pname);

  pl->runflag = 0;

  pl->lgn_in_flag  = 0;  // Assume no LGN input, until INPUTs are read
  pl->lgn_in_n     = 0;  // Assume no LGN input
  pl->lgn_gain     = 0;  // Assume no LGN gain signal
  pl->lgn_gain_tau = 0.010;
  pl->lgn_gain_ca  = 1.0;
  pl->lgn_gain_cb  = 1.0;

  pl->laytype = NULL;
  pl->icon = NULL;
  pl->area = NULL;

  pl->ninlist = 0;
  pl->inlist = NULL;

  pl->geomt = NULL;
  pl->xn = -1;
  pl->yn = -1;
  pl->zn = -1;

  pl->x0  = pl->y0  = 0.0;
  pl->xf  = pl->yf  = 0.0;
  pl->umx = pl->umy = 0.0;

  pl->oddxoff = 0;

  pl->irr_n = 0;
  pl->irr_x = NULL;
  pl->irr_y = NULL;
  pl->irr_id = NULL;


  pl->save_name_list = NULL;

  pl->layseed_n = -1;
  pl->layseed = NULL;

  pl->nmech = 0;
  pl->mech = NULL;

  pl->c = NULL;

  pl->attrib_fn = 0;
  pl->attrib_fname = NULL;

  pl->ifcp = NULL;
  pl->poissp = NULL;
  pl->nmda_p = NULL;

  pl->csav_n = 0;
  pl->csav_datid = NULL;
  pl->csav_source = NULL;
  pl->csav_samp = NULL;

  pl->cnt = NULL;
  pl->s = NULL;
  pl->rl = NULL;
  pl->rr = NULL;
  pl->rgl = NULL;      // Gain signals
  pl->rgr = NULL;

  return pl;
}
/**************************************-**************************************/
/*                                                                           */
/*                           POPU_MAKE_LAYER_IRREG                           */
/*                                                                           */
/*****************************************************************************/
struct pop_layer *popu_make_layer_irreg(pname,n,zn,cx,cy,rfx,rfy)
     char *pname;      // name for layer
     int n;            // Number of locations
     int zn;           // Number of unit per location
     float *cx,*cy;    // position of units in tissue (0,0) is center
     float *rfx,*rfy;  // rf center position in [0..xn) [0..yn)
{
  int i,j,k;
  struct pop_layer *pl;
  struct pop_cell ***c,*cp;

  pl = popu_make_pop_layer(pname);
  pl->laytype = strdup("irreg");

  pl->xn = n;
  pl->yn = 1;
  pl->zn = zn;
  
  c = (struct pop_cell ***)myalloc(pl->xn*sizeof(struct pop_cell **));
  for(i=0;i<pl->xn;i++){
    c[i] = (struct pop_cell **)myalloc(pl->yn*sizeof(struct pop_cell *));
    for(j=0;j<pl->yn;j++){
      c[i][j] = (struct pop_cell *)myalloc(pl->zn*sizeof(struct pop_cell));

      for(k=0;k<pl->zn;k++){
	cp = &(c[i][j][k]);
	pop_cell_init_cell(cp);
	
	// Position in tissue, (0,0) is center
	cp->cx = cx[i];
	cp->cy = cy[i];
	
	// RF center in stimulus grid [0..xn][0..yn]
	cp->rfx = rfx[i];
	cp->rfy = rfy[i];

	// WYETH - Should we store any attribs here? such as cone type??

	//  I added this, because I found it was needed below for ..._CART,
	//  so we should be sure it is needed here?
	cp->layx = i;  // WYETH - added Apr 4, 2013
	cp->layy = j;  // WYETH - added Apr 4, 2013
	cp->layz = k;  // WYETH - added Apr 4, 2013

      }
    }
  }

  pl->c = c;

  return pl;
}
/**************************************-**************************************/
/*                                                                           */
/*                            POPU_MAKE_LAYER_CART                           */
/*                                                                           */
/*  Cells in Cartesian coords.                                               */
/*                                                                           */
/*****************************************************************************/
struct pop_layer *popu_make_layer_cart(pname,xn,yn,zn)
     char *pname;        // Population name
     int xn,yn,zn;       // Dimensions of layer
{
  struct pop_layer *pl;

  pl = popu_make_pop_layer(pname);
  pl->laytype = strdup("default");

  pl->xn = xn;
  pl->yn = yn;
  pl->zn = zn;

  pl->c = popu_make_init_cells(xn,yn,zn,pl);

  return pl;
}
/**************************************-**************************************/
/*                                                                           */
/*                          POPU_MAKE_LAYER_TEMPLATE                         */
/*                                                                           */
/*  Create a template layer, for use in building .mar files.                 */
/*  The idea here is that one set of input connections ('fx' and 'fy') is    */
/*  specified for each 'z' layer.  Thus, there is no variation in the input  */
/*  pattern (from 'lpre') across units in (x,y) here.                        */
/*                                                                           */
/*****************************************************************************/
struct pop_layer *popu_make_layer_template(name,xn,yn,zn,pa,lpre,fx,fy,fw,fn)
     char *name;
     int xn,yn,zn;  // Each 'z' index will have a single set of connections
     struct pop_area *pa;
     struct pop_layer *lpre;
     int   **fx;  // [zn][fn] x-coord for synaptic connection, in 'lpre'
     int   **fy;  // [zn][fn] y-coord for synaptic connection, in 'lpre'
     float **fw;  // [zn][fn] weights for synaptic connection
     int    *fn;  // [zn]     Number of 'inputs'
{
  int i;
  int zi,stype,inindex,pxn,pyn;
  float w,tdelay;
  struct pop_layer *pl;
  struct onode *ot,*o;
  struct pop_cell *post,*pre;
  struct pop_syn *tsyn;

  pl = popu_make_layer_cart(name,xn,yn,zn);
  pl->geomt = strdup("default");
  pl->area = pa;
  pl->x0 = pl->y0 = 0;
  pl->xf = pl->yf = 1.0;
  pl->ninlist = 1;
  pl->inlist = (struct onode **)myalloc(sizeof(struct onode *));

  // Create an onode for the input
  ot = paramfile_onode_create_onode();
  ot->otype = strdup("input");
  o = paramfile_onode_get_init_onode("item",NULL,"pop_origin",lpre->name,
				     1,NULL,NULL,1);
  onode_insert_child_at_end(ot,o);

  // *** WYETH - I ASSUME this 'template' type means that only the connections
  //   for the first unit need to be stored.
  o = paramfile_onode_get_init_onode("item",NULL,"type","template",
				     1,NULL,NULL,1);
  onode_insert_child_at_end(ot,o);
  pl->inlist[0] = ot;

  //
  //  Write connections for first unit
  //
  stype = 1;
  tdelay = 0.0;
  inindex = 0;
  w = 1.0;

  pxn = lpre->xn;
  pyn = lpre->yn;

  for(zi=0;zi<zn;zi++){
    post = &(pl->c[0][0][zi]);

    for(i=0;i<fn[zi];i++){
      //printf("fx[i] fy[i] = %d  %d    %f\n",fx[zi][i],fy[zi][i],fw[zi][i]);
      pre = &(lpre->c[fx[zi][i]][fy[zi][i]][0]);
      tsyn = pop_cell_add_synapse(pre,post,stype,fw[zi][i],tdelay,inindex);
    }
  }

  return pl;
}
/**************************************-**************************************/
/*                                                                           */
/*                              POPU_MAKE_POP_TOP                            */
/*                                                                           */
/*****************************************************************************/
struct pop_top *popu_make_pop_top(logf,sscale,tscale,xn,yn,tn)
     char *logf;
     float sscale;
     float tscale;
     int xn,yn,tn;
{
  struct pop_top *mpt;

  // Basic global params, determine duration of sim in seconds
  mpt = (struct pop_top *)myalloc(sizeof(struct pop_top));

  mpt->name = NULL;
  
  if  (logf == NULL){
    mpt->logf = NULL;
  }else
    mpt->logf = strdup(logf);
  
  mpt->tscale = tscale;
  mpt->sscale = sscale;
  mpt->xn = xn;
  mpt->yn = yn;
  mpt->tn = tn;
  mpt->x1 = 0.0;
  mpt->x2 = mpt->tscale * (float)mpt->tn;
  mpt->out_prefix = NULL;

  mpt->binoc = 0;     // Assume monocular left eye, unless 'lgn_on_r' below
  mpt->retflag = 0;   // Assume DOG filter

  mpt->narea = 0;
  mpt->area = NULL;

  mpt->nlay = 0;
  mpt->lay = NULL;

  mpt->nmap = 0;
  mpt->map = NULL;

  mpt->conn_rw = NULL;
  mpt->conn_read_dir = NULL;
  mpt->conn_write_dir = NULL;

  return mpt;
}
/**************************************-**************************************/
/*                                                                           */
/*                           POPU_MAKE_AREA_DEFAULT                          */
/*                                                                           */
/*****************************************************************************/
struct pop_area *popu_make_area_default(name,xn,yn,sscale,um)
     char *name;
     int xn,yn;
     float sscale;
     float um;
{
  struct pop_area *pa;

  pa = (struct pop_area *)myalloc(sizeof(struct pop_area));

  pa->name = strdup(name);
  pa->x0   = 0;
  pa->y0   = 0;

  pa->xf   = 1.0;
  pa->yf   = 1.0;

  pa->xn   = xn;
  pa->yn   = yn;
  pa->umx  = um;
  pa->umy  = um;

  pa->sscale = sscale;

  return pa;
}
/**************************************-**************************************/
/*                                                                           */
/*                               POPU_MAKE_AREA                              */
/*                                                                           */
/*****************************************************************************/
struct pop_area *popu_make_area(logf,ao,sscale)
     char *logf;
     struct onode *ao;  // <area> onode
     float sscale;      // e.g., mpt->sscale
{
  struct pop_area *pa;

  pa = (struct pop_area *)myalloc(sizeof(struct pop_area));

  pa->name = onode_getpar_chr_exit(ao,"name");
  pa->x0   = onode_getpar_flt_exit(ao,"x0");
  pa->y0   = onode_getpar_flt_exit(ao,"y0");

  pa->xf   = onode_getpar_flt_dflt(ao,"xf",1.0);
  pa->yf   = onode_getpar_flt_dflt(ao,"yf",1.0);

  pa->xn   = onode_getpar_int_exit(ao,"xn");
  pa->yn   = onode_getpar_int_exit(ao,"yn");
  pa->umx  = onode_getpar_flt_dflt(ao,"umx",0.0);
  pa->umy  = onode_getpar_flt_dflt(ao,"umy",pa->umx);

  pa->sscale = sscale;

  return pa;
}
/**************************************-**************************************/
/*                                                                           */
/*                              POPU_INIT_GEOM                               */
/*                                                                           */
/*****************************************************************************/
void popu_init_geom(logf,o,po,pl,retflag)
     char *logf;
     struct onode *o;    // Head onode
     struct onode *po;   // Population onode
     struct pop_layer *pl;
     int retflag;             // e.g., mpt->retflag
{
  struct onode *geo;

  geo = onode_child_get_unique(po,"geometry");

  pl->geomt = onode_getpar_chr_dflt(geo,"type","default");

  //
  //  Number of cells in layer.  For "irregular", caller must configure later.
  //
  if (pop_util_lclass_is_lgn0(pl) ||  // e.g., "lgn" or square retina0_gc0
      ((strcmp(pl->laytype,"retina0_gc0")==0) && (retflag==100))){
    pl->xn = onode_getpar_int_exit(o,"xn"); // Inherit top level xn, yn
    pl->yn = onode_getpar_int_exit(o,"yn");
    //printf("    %s getting top level xn,yn\n",pl->name);
  }else if (pl->area != NULL){
    pl->xn = onode_getpar_int_dflt(geo,"xn",pl->area->xn); // May inherit area
    pl->yn = onode_getpar_int_dflt(geo,"yn",pl->area->yn);
    //printf("    %s getting area xn,yn\n",pl->name);
  }else{
    pl->xn = onode_getpar_int_exit(geo,"xn");
    pl->yn = onode_getpar_int_exit(geo,"yn");
    //printf("    %s getting own xn,yn\n",pl->name);
  }
  pl->zn = onode_getpar_int_dflt(geo,"zn",1);

  //
  //  Origin and scaling of layer wrt overall 'xn' 'yn'
  //
  pl->x0      = onode_getpar_flt_dflt(geo,"x0",0.0);
  pl->y0      = onode_getpar_flt_dflt(geo,"y0",0.0);
  pl->xf      = onode_getpar_flt_dflt(geo,"xf",1.0);
  pl->yf      = onode_getpar_flt_dflt(geo,"yf",1.0);
  pl->oddxoff = onode_getpar_flt_dflt(geo,"oddxoff",0.0);

}
/**************************************-**************************************/
/*                                                                           */
/*                            POPU_INPUT_GET_ORIGIN                          */
/*                                                                           */
/*  Return the pop (layer) of origin for the input index 'iin'.              */
/*                                                                           */
/*****************************************************************************/
struct pop_layer *popu_input_get_origin(mpt,lay,iin)
     struct pop_top *mpt;     // Model top
     struct pop_layer *lay;   // Layer post-syn for input
     int iin;                 // input index
{
  struct pop_layer *l;
  struct onode *o;

  o = lay->inlist[iin];
  l = popu_get_lay_named(mpt,o,"pop_origin","POPU_INPUT_GET_ORIGIN");

  return l;
}
/**************************************-**************************************/
/*                                                                           */
/*                            POPU_INPUT_IS_PAIRED                           */
/*                                                                           */
/*  Return '1' if the input is paired, '0' otherwise.                        */
/*                                                                           */
/*****************************************************************************/
int popu_input_is_paired(lay,iin)
     struct pop_layer *lay;
     int iin;
{
  int flag;
  struct onode *o;

  // WYETH - THIS MAY NOT BE COMPLETE - WHAT ABOUT OTHER PAIRING TYPES???
  // WYETH - THIS MAY NOT BE COMPLETE - WHAT ABOUT OTHER PAIRING TYPES???
  // WYETH - THIS MAY NOT BE COMPLETE - WHAT ABOUT OTHER PAIRING TYPES???

  flag = 0;

  if (lay->inlist != NULL){
    o = onode_child_get_unique(lay->inlist[iin],"input_pair");
    if (o != NULL)
      flag = 1;

    // WYETH - added for 'syn_int ds02' which makes a pairing
    //o = onode_child_get_unique(lay->inlist[iin],"syn_int");
    o = onode_get_item_child(lay->inlist[iin],"syn_int");
    if (o != NULL)
      flag = 1;
  }
  
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                              POPU_INPUT_IS_BG                             */
/*                                                                           */
/*  Return '1' if the input is paired, '0' otherwise.                        */
/*                                                                           */
/*****************************************************************************/
int popu_input_is_bg(lay,iin)
     struct pop_layer *lay;
     int iin;
{
  int flag;
  struct onode *o;

  // WYETH - THIS MAY NOT BE COMPLETE - WHAT ABOUT OTHER PAIRING TYPES???
  // WYETH - THIS MAY NOT BE COMPLETE - WHAT ABOUT OTHER PAIRING TYPES???
  // WYETH - THIS MAY NOT BE COMPLETE - WHAT ABOUT OTHER PAIRING TYPES???

  flag = 0;

  o = lay->inlist[iin];
  flag = onode_test_chr(o,"type","bg");


  //printf("HERE WYETH -  iin = %d\n",iin);
  //printf("HERE WYETH -  layname = %s\n",lay->name);

  if (o != NULL)
    flag = 1;
  
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                          POPU_INPUT_GET_PAIRED_LAY                        */
/*                                                                           */
/*  Return a pointer to the paired layer, or NULL if no pairing.             */
/*                                                                           */
/*****************************************************************************/
struct pop_layer *popu_input_get_paired_lay(mpt,lay,iin)
     struct pop_top *mpt;     // Model top
     struct pop_layer *lay;
     int iin;
{
  struct onode *po,*io;
  struct pop_layer *lay2;

  lay2 = NULL;

  if (popu_input_is_paired(lay,iin)){
    
    io = lay->inlist[iin];   // <input> onode
    po = onode_child_get_unique(io,"input_pair");  // <input_pair>

    if (po == NULL){
      // WYETH - assuming default is 'pop_origin'
      // WYETH - this handles 'syn_int ds02'
      lay2 = popu_get_lay_named(mpt,io,"pop_origin",
				"POPU_INPUT_GET_PAIRED_LAY");
    }else{
      if (onode_get_item_child(po,"pop_origin_2") != NULL)
	lay2 = popu_get_lay_named(mpt,po,"pop_origin_2",
				  "POPU_INPUT_GET_PAIRED_LAY");
      else // WYETH - assuming default is 'pop_origin'
	lay2 = popu_get_lay_named(mpt,io,"pop_origin",
				  "POPU_INPUT_GET_PAIRED_LAY");
    }
    
  }  
  return lay2;
}
/**************************************-**************************************/
/*                                                                           */
/*                           POP_UTIL_GET_NMDA_SHAPE                         */
/*                                                                           */
/*****************************************************************************/
float *pop_util_get_nmda_shape(alph1,alph2,alph3,amp1,amp2,amp3,n)
     float alph1,alph2,alph3,amp1,amp2,amp3;
     int n;
{
  float *g,*gt;

  g  = alpha_farray(0.0,alph1,amp1,n);
  gt = alpha_farray(0.0,alph2,amp2,n);
  add_to_farray(g,gt,n);
  myfree(gt);

  gt = alpha_farray(0.0,alph3,amp3,n);
  add_to_farray(g,gt,n);
  myfree(gt);

  return g;
}
/**************************************-**************************************/
/*                                                                           */
/*                          POP_UTIL_COMPUTE_NMDA_TABLE                      */
/*                                                                           */
/*  Values from Jahr and Stevens (1990) table 1 and Eqn 4a.                  */
/*                                                                           */
/*****************************************************************************/
void pop_util_compute_nmda_table(ppmv,v1,v2,rg,rn)
     int ppmv;    /* Number of points in table per mV */
     int v1,v2;   /* Table goes from v1 to v2 in mV */
     float **rg;  /* Return the table */
     int *rn;     /* Number of entries in table */
{
  int i;
  int nmv,n;
  float *g,*x,C,a1,a2,b1,b2,A,B1,B2,v,t1,t2;

  nmv = v2 - v1 + 1; /* Number of mV covered in table */
  n = nmv * ppmv;  /* number of points in table */
  g = get_zero_farray(n);
  x = get_zero_farray(n);
  for(i=0;i<n;i++){
    v = v1 + (float)i/(float)ppmv; /* Voltage (mV) */
    C = 100.0; /* Mg conc (uM) */
    a1 = exp(-0.016 * v - 2.91);
    a2 = C * exp(-0.045 * v - 6.97); 
    b1 = exp(0.009 * v + 1.22);
    b2 = exp(0.017 * v + 0.96);
    A = exp(-2.847);
    B1 = exp(-0.693);
    B2 = exp(-3.101);

    t1 = (a1 + a2)*(a1*B1 + a2*B2);
    t2 = A*a1*(b1 + B1) + A*a2*(b2 + B2);

    g[i] = 1.0/(1.0 + t1/t2);
    x[i] = v;
  }
  /*append_farray_xy_plot("zzz.NMDA.table.pl",x,g,n,"g(v)");*/

  myfree(x);

  *rg = g; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                     POP_UTIL_NMDA_MATCH_INTCURR_WAVEFORM                  */
/*                                                                           */
/*  Compute integrated current to compare EPSCs for AMPA and NMDA at V.      */
/*                                                                           */
/*****************************************************************************/
void pop_util_nmda_match_intcurr_waveform(myid,mylogf,gwav,n,gov,gov_n,gov_v0,
					  gov_ppmv,nmda_frac,v,v_rev,
					  alph1,alph2,alph3,ramp1,ramp2,ramp3,
					  outprefix,name)
     int myid;         // For cluster computing
     char *mylogf;     // For cluster computing
     float *gwav;      // waveform to match, typically AMPA
     int n;            // length of gwav
     float *gov;       // table of NMDA voltage dependence g(v)
     int gov_n;        // length of 'gov'
     float gov_v0;     // Initial voltage in table (mV)
     int gov_ppmv;     // points per mV in 'gov'
     float nmda_frac;  // fraction of total integrated current for NMDA
     float v;          // voltage to match at
     float v_rev;      // reversal potential for both
     float alph1,alph2,alph3;
     float *ramp1,*ramp2,*ramp3;   // These must contain starting values
     char outprefix[];
     char name[];
{
  int k;
  int nnmda;
  float *gnmda,ampa_curr,nmda_curr,ff,rescale;
  char outfile[SLEN],tstr[SLEN];

  mylog(mylogf,"  POP_UTIL_NMDA_MATCH_INTCURR_WAVEFORM\n");
  sprintf(tstr,"    Matching at %.2f mV\n",v);
  mylog(mylogf,tstr);

  nnmda = (int)(1.0/alph3 * 10.0);

  gnmda = pop_util_get_nmda_shape(alph1,alph2,alph3,*ramp1,*ramp2,*ramp3,
				  nnmda);

  /*
  if (myid == -1){
    sprintf(outfile,"%s.%s.NMDA.pl",outprefix,name);
    append_farray_plot(outfile,"AMPA_PSC",gwav,n,1);
    append_farray_plot(outfile,"NMDA_PSC_unscaled",gnmda,nnmda,1);
    }*/

  ampa_curr = sum_farray(gwav,n,0,n) * (v_rev - v);

  k = (v - gov_v0) * gov_ppmv;
  if ((k < 0) || (k >= gov_n))
    mylog_exit(mylogf,
	       "POP_UTIL_NMDA_MATCH_INTCURR_WAVEFORM  NMDA table index");
  ff = gov[k];
  sprintf(tstr,"    NMDA voltage-dependent factor at %.2f mV is %f\n",v,ff);
  mylog(mylogf,tstr);

  nmda_curr = ff * sum_farray(gnmda,nnmda,0,nnmda) * (v_rev - v);

  sprintf(tstr,"    Integrated current AMPA = %f\n",ampa_curr);
  mylog(mylogf,tstr);
  sprintf(tstr,"    Integrated current NMDA = %f\n",nmda_curr);
  mylog(mylogf,tstr);
  sprintf(tstr,"    NMDA fraction target is %f\n",nmda_frac);
  mylog(mylogf,tstr);

  rescale = ampa_curr/nmda_curr * nmda_frac/(1.0 - nmda_frac);
  sprintf(tstr,"      rescaling NMDA amplitude by %f\n",rescale);
  mylog(mylogf,tstr);
  *ramp1 *= rescale;
  *ramp2 *= rescale;
  *ramp3 *= rescale;

  myfree(gnmda);
  gnmda = pop_util_get_nmda_shape(alph1,alph2,alph3,*ramp1,*ramp2,*ramp3,
				  nnmda);

  nmda_curr = ff * sum_farray(gnmda,nnmda,0,nnmda) * (v_rev - v);
  sprintf(tstr,"    Integrated current NMDA = %f\n",nmda_curr);
  mylog(mylogf,tstr);
  sprintf(tstr,"    NMDA / (AMPA + NMDA) = %f\n",
	  nmda_curr/(ampa_curr+nmda_curr));
  mylog(mylogf,tstr);

  multiply_farray(gnmda,nnmda,ff);
  /*
    if (myid == -1)
    append_farray_plot(outfile,"NMDA_PSC_NEW_x_g(V)",gnmda,nnmda,1);*/

  myfree(gnmda);
}
/**************************************-**************************************/
/*                                                                           */
/*                   POP_UTIL_COMPUTE_INTCURR_AMPA_NMDA_LAYER                */
/*                                                                           */
/*  Compute integrated current to compare EPSCs for AMPA and NMDA at V.      */
/*                                                                           */
/*****************************************************************************/
void pop_util_compute_intcurr_ampa_nmda_layer(myid,mylogf,lay,v,outprefix)
     int myid;
     char *mylogf;
     struct pop_layer *lay;
     float v;
     char outprefix[];
{
  int nampa;
  float *gampa;
  char tstr[SLEN];
  struct pop_cell *c;
  struct nmda_param *tp;
  
  mylog(mylogf,"  POP_UTIL_COMPUTE_INTCURR_AMPA_NMDA_LAYER\n");

  c = &(lay->c[0][0][0]);
  tp = lay->nmda_p;

  nampa = (int)(c->ex1tau*1000.0 * 10.0); /* length in msec */
  gampa = alpha_farray(0.0,1.0/(c->ex1tau*1000.0),c->ex1amp,nampa);
  sprintf(tstr,"      AMPA alpha-func tau= %.4f  ampl= %.4f\n",c->ex1tau,
	  c->ex1amp);
  mylog(mylogf,tstr);

  pop_util_nmda_match_intcurr_waveform(myid,mylogf,gampa,nampa,tp->nmda_gv,
				       tp->nmda_n,
				       tp->nmda_v0f,tp->nmda_ppmv,
				       tp->nmda_frac,v,c->cifcp->v_ex,
				       tp->nmda_alpha_1,tp->nmda_alpha_2,
				       tp->nmda_alpha_3,&(tp->nmda_amp_1),
				       &(tp->nmda_amp_2),&(tp->nmda_amp_3),
				       outprefix,"Layer");

  myfree(gampa);
}
/**************************************-**************************************/
/*                                                                           */
/*                    POP_UTIL_COMPUTE_INTCURR_AMPA_NMDA_LGN                 */
/*                                                                           */
/*****************************************************************************/
void pop_util_compute_intcurr_ampa_nmda_lgn(myid,mylogf,lay,m,v,outprefix)
     int myid;
     char *mylogf;
     struct pop_layer *lay;
     struct model_struct *m;  /* Model parameters. */
     float v;
     char outprefix[];
{
  int gtn;
  float *gt;
  char name[SLEN];
  struct pop_cell *c;
  struct pop_mech *mech;
  struct nmda_param *tp;

  mylog(mylogf,"  POP_UTIL_COMPUTE_INTCURR_AMPA_NMDA_LGN\n");

  c = &(lay->c[0][0][0]);
  tp = lay->nmda_p;

  // WYETH - this wants to find the "lgn2ctx" EPSG function
  if (m->ppl != NULL){
    strcpy(name,"lgn2ctx");
    gt = pop_util_get_mask_diff_exp(m,name,&gtn);
  }else{
    mech = pop_cell_get_mech_by_name(lay,"lgn_ex");
    gt = pop_util_get_mask_diff_exp_mech(mech,&gtn);
  }

  pop_util_nmda_match_intcurr_waveform(myid,mylogf,gt,gtn,tp->nmda_gv,
				       tp->nmda_n,
				       tp->nmda_v0f,tp->nmda_ppmv,
				       tp->lgn_nmda_frac,v,c->cifcp->v_ex,
				       tp->nmda_alpha_1,tp->nmda_alpha_2,
				       tp->nmda_alpha_3,
				       &(tp->lgn_nmda_amp_1),
				       &(tp->lgn_nmda_amp_2),
				       &(tp->lgn_nmda_amp_3),
				       outprefix,"LGN");
  myfree(gt);
}
/**************************************-**************************************/
/*                                                                           */
/*                      POP_UTIL_NMDA_PREP_CONSTANTS_LAYER                   */
/*                                                                           */
/*  These constants are set to keep 'pop_input_02' clean.                    */
/*                                                                           */
/*****************************************************************************/
void pop_util_nmda_prep_constants_layer(lay)
     struct pop_layer *lay;
{
  float alph1,alph2,alph3;
  struct nmda_param *tp;

  tp = lay->nmda_p;

  // Note, alpha is 1/tau
  alph1 = 1000.0 * tp->nmda_alpha_1; // convert to 1/sec
  alph2 = 1000.0 * tp->nmda_alpha_2; // convert to 1/sec
  alph3 = 1000.0 * tp->nmda_alpha_3; // convert to 1/sec

  tp->nmda1_aetau = tp->nmda_amp_1 * M_E * alph1;
  tp->nmda2_aetau = tp->nmda_amp_2 * M_E * alph2;
  tp->nmda3_aetau = tp->nmda_amp_3 * M_E * alph3;

  tp->nmda1_tau2inv = alph1 * alph1;
  tp->nmda2_tau2inv = alph2 * alph2;
  tp->nmda3_tau2inv = alph3 * alph3;

  tp->nmda1_totau = 2.0 * alph1;
  tp->nmda2_totau = 2.0 * alph2;
  tp->nmda3_totau = 2.0 * alph3;
}
/**************************************-**************************************/
/*                                                                           */
/*                               POPU_FLAG_RESP                              */
/*                                                                           */
/*  Return the number of new marks in the 'cmap' array.                      */
/*                                                                           */
/*****************************************************************************/
int popu_flag_resp(r,lay,cmap)
     struct response_struct *r;     // Response data
     struct pop_layer *lay;         // layer in which responses are req'd
     int ***cmap;                   // set to '1' for each response req'd
{
  int i,j;
  int xi,yi,zi,new,tot;

  // For each response
  new = tot = 0;
  for(i=0;i<r->n;i++){

    // If this response in in our 'lay'
    if (strcmp(r->plname[i],lay->name)==0){
      //printf("MATCHED  RESP: %s    LAY: %s\n",r->plname[i],lay->name);

      xi = r->xi[i];
      yi = r->yi[i];
      zi = r->zi[i];
      
      if (cmap[xi][yi][zi] == 0){
	new += 1;
	cmap[xi][yi][zi] = 1;
      }
      tot += 1;
    }
  }

  return new;
  //printf("    %d new marks of %d total responses in layer %s\n",new,tot,
  //lay->name);
}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_UTIL_RESPONSE_LINK                         */
/*                                                                           */
/*****************************************************************************/
void pop_util_response_link(r,laylist,nlay,mylogf)
     struct response_struct *r;     // Response data
     struct pop_layer **laylist;    // List of layers
     int nlay;                      // Number of layers
     char *mylogf;
{
  int i;
  int nmdaflag;
  struct pop_cell *c;
  struct pop_syn *syn;
  struct pop_layer *pl;
  char tstr[SLEN];

  mylog(mylogf,"  POP_UTIL_RESPONSE_LINK\n");

  for(i=0;i<r->n;i++){ // For each response requested

    //
    //  Get a pointer to the relevant cell or synapse
    //
    if (r->rformat[i] == 0){
      sprintf(tstr,"*** POP_UTIL_RESPONSE_LINK  Format 0 not impl'd yet");
      mylog_exit(mylogf,tstr);
    }else if (r->rformat[i] == 10){
      //
      //  Synapse
      //
      syn = pop_cell_laylist_get_syn_pointer(laylist,nlay,r->plname[i],
					     r->xi[i],r->yi[i],r->zi[i],
					     r->plname1[i],
					     r->xi1[i],r->yi1[i],r->zi1[i],0);
      if (syn == NULL){
	sprintf(tstr,"*** POP_UTIL_RESPONSE_LINK  No synapse found for %s\n",
		r->nd_name[i]);
	mylog_exit(mylogf,tstr);
      }
      pop_cell_syn_set_save(mylogf,syn,r->datid[i]);
    }else{
      //
      //  Cell
      //
      c = pop_cell_laylist_get_cell_pointer(laylist,nlay,r->plname[i],
					    r->xi[i],r->yi[i],r->zi[i]);
      if (c == NULL){
	sprintf(tstr,"*** POP_UTIL_RESPONSE_LINK  No cell found for %s\n",
		r->nd_name[i]);
	mylog_exit(mylogf,tstr);
      }

      if (c->pl->ifcp != NULL){  // Note, LGN pops don't have ifcp ?
	if (c->pl->ifcp->g_tran > 0.0)
	  nmdaflag = 1;
	else
	  nmdaflag = 0;
      }

      pop_cell_set_save(mylogf,c,r->datid[i],nmdaflag);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         POP_UTIL_RESPONSE_HANDOVER                        */
/*                                                                           */
/*  Response Handover                                                        */
/*                                                                           */
/*****************************************************************************/
void pop_util_response_handover(m,r,laylist,nlay,mylogf)
     struct model_struct *m;        // Model parameter pair list
     struct response_struct *r;     // Response data
     struct pop_layer **laylist;    // List of layers
     int nlay;                      // Number of layers
     char *mylogf;
{
  int i;
  int flag,din,dn,csi,csav_flag;
  float *gtot,*di,*d;
  struct pop_layer *pl;
  struct pop_cell *c;
  struct pop_syn *syn;
  char tstr[SLEN];

  for(i=0;i<r->n;i++){ // For each response requested

    //  Get a pointer to the relevant cell
    if (r->rformat[i] == 0){
      sprintf(tstr,"*** POP_UTIL_RESPONSE_HANDOVER  Format 0 not impl'd yet");
      mylog_exit(mylogf,tstr);
    }else if (r->rformat[i] == 10){
      syn = pop_cell_laylist_get_syn_pointer(laylist,nlay,r->plname[i],
					     r->xi[i],r->yi[i],r->zi[i],
					     r->plname1[i],
					     r->xi1[i],r->yi1[i],r->zi1[i],0);

      if (syn == NULL){
	sprintf(tstr,"*** POP_UTIL_RESPONSE_HANDOVER  No syn found for %s\n",
		r->nd_name[i]);
	mylog_exit(mylogf,tstr);
      }
    }else{
      c = pop_cell_laylist_get_cell_pointer(laylist,nlay,r->plname[i],
					    r->xi[i],r->yi[i],r->zi[i]);
      if (c == NULL){
	sprintf(tstr,"*** POP_UTIL_RESPONSE_HANDOVER  No cell found for %s\n",
		r->nd_name[i]);
	mylog_exit(mylogf,tstr);
      }
    }

    //  Handle requests for all valid response entitities

    if (strcmp(r->datid[i],"spikes")==0){
      if (!c->savs)
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER Spikes not saved\n");

      /*** WYETH - must check sampling units ?? ***/

      //printf("-------(pop_util)--RESP HANDOVER name= %s  nspikes= %d  r->gtsi= %d\n",
      //r->plname[i],c->ns,r->gtsi);

      mod_util_resp_store_s(r,i,c->s,c->ns,1,mylogf);

    }else if (strcmp(r->datid[i],"lgn_gx")==0){
      if (!c->savg)
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER g not saved\n");
      if (r->samp[i] != c->samp0) /* Should do resampling here? */
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER sampl mismatch\n");

      if (c->pl->nmda_p != NULL){
	gtot = add_farrays(c->gtx0,c->gtn0,c->n0);
	mod_util_resp_store_f(r,i,gtot,c->n0,0,mylogf); /* 0-Don't copy */
      }else{
	mod_util_resp_store_f(r,i,c->gtx0,c->n0,1,mylogf); /* 1-Copy */
      }
    }else if (strcmp(r->datid[i],"lgn_gain")==0){
      if (c->ggain != NULL)
	mod_util_resp_store_f(r,i,c->ggain,c->n0,1,mylogf); // 1-Copy
      else
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER no gain signal\n");

    }else if (strcmp(r->datid[i],"lgn_ga")==0){
      if (!c->savg)
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER g not saved\n");
      if (r->samp[i] != c->samp0) /* Should do resampling here? */
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER sampl mismatch\n");

      mod_util_resp_store_f(r,i,c->gtx0,c->n0,1,mylogf);

    }else if (strcmp(r->datid[i],"lgn_ia")==0){
      if ((!c->savg) || (!c->savv))
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER g or v not saved\n");
      if (r->samp[i] != c->samp0) /* Should do resampling here? */
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER sampl mismatch\n");

      /* Compute derived current */
      di = mod_util_resp_derive_current(c->gtx0,c->samp0,(float *)NULL,c->n0,
					1.0,c->vm,c->vmt,c->vn,
					c->cifcp->v_th_x,
					c->cifcp->v_ex,
					mylogf,&din);
      /*
      if (din != c->n0){
	printf("  din = %d   c->n0 = %d\n",din,c->n0);
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER din error\n");
	}*/

      mod_util_resp_store_f(r,i,di,c->n0,1,mylogf);
      myfree(di);

    }else if (strcmp(r->datid[i],"in_ii")==0){
      if ((!c->savg) || (!c->savv))
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER g or v not saved\n");

      // Compute derived current
      di = mod_util_resp_derive_current(c->gtin1,0.0,c->gtex1t,c->gtex1n,
					0.001,c->vm,c->vmt,c->vn,
					c->cifcp->v_th_x,
					c->cifcp->v_in,
					mylogf,&din);
      /*
	if (din != c->gtex1n){
	printf("  din = %d   c->n0 = %d\n",din,c->n0);
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER din error\n");
	}*/


      mod_util_resp_store_f_reg(r,i,c->gtex1t,di,din,c->samp0,mylogf);

      myfree(di);

    }else if (strcmp(r->datid[i],"ex_ix")==0){
      if ((!c->savg) || (!c->savv))
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER g or v not saved\n");

      // Compute derived current
      di = mod_util_resp_derive_current(c->gtex1,0.0,c->gtex1t,c->gtex1n,
					0.001,c->vm,c->vmt,c->vn,
					c->cifcp->v_th_x,
					c->cifcp->v_ex,
					mylogf,&din);
      mod_util_resp_store_f_reg(r,i,c->gtex1t,di,din,c->samp0,mylogf);
      myfree(di);
    }else if (strcmp(r->datid[i],"lgn_gn")==0){
      if (!c->savg)
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER g not saved\n");
      if (r->samp[i] != c->samp0) /* Should do resampling here? */
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER sampl mismatch\n");

      if (c->pl->nmda_p == NULL){
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER  No NMDA\n");
      }

      mod_util_resp_store_f(r,i,c->gtn0,c->n0,1,mylogf);

    }else if (strcmp(r->datid[i],"ex_gx")==0){
      if (!c->savg)
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER g not saved\n");
      
      mod_util_resp_store_f_reg(r,i,c->gtex1t,c->gtex1,c->gtex1n,c->samp0,
				mylogf);
    }else if (strcmp(r->datid[i],"in_gi")==0){
      if (!c->savg)
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER g not saved\n");

      mod_util_resp_store_f_reg(r,i,c->gtex1t,c->gtin1,c->gtex1n,c->samp0,
				mylogf);
    }else if (strcmp(r->datid[i],"gad")==0){
      if (!c->savg)
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER g not saved\n");
      if (r->samp[i] != c->samp0) // Should do resampling here?
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER sampl mismatch\n");


      if (c->gta0 == NULL)  // WYETH - I added this, but is it needed?
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER gta0 is NULL\n");

      mod_util_resp_store_f(r,i,c->gta0,c->n0,1,mylogf);

      /***
      // wygad
      // wygad
      {
	float *tz;
	// wygad
	tz = get_zero_farray(c->n0);
	tz[4] = 1.0;
	mod_util_resp_store_f(r,i,tz,c->n0,1,mylogf);
      }
      // wygad
      // wygad
      ***/

    }else if ((strcmp(r->datid[i],"vm")==0)||
	      (strcmp(r->datid[i],"rate")==0)){ // For probabilistic spike gen
      if (strcmp(r->rtype[i],"f")!=0)
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER Vm not 'f'\n");
      if (!c->savv)
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER Vm not saved\n");

      // Sampling is sec for Vm

      if (c->vm == NULL)
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER  vm is NULL\n");

      if (c->vmt != NULL)
	mod_util_resp_store_f_reg(r,i,c->vmt,c->vm,c->vn,1.0,mylogf);
      else{ // LGN layers only
	mod_util_resp_store_f(r,i,c->vm,c->vn,1,mylogf);
      }

    }else if (strcmp(r->datid[i],"gx")==0){   // LGN layers only
      mod_util_resp_store_f(r,i,c->gtx0,c->n0,1,mylogf); // 1-Copy

    }else if (strcmp(r->datid[i],"gi")==0){   // LGN layers only
      mod_util_resp_store_f(r,i,c->gti0,c->n0,1,mylogf); // 1-Copy

    }else if (strcmp(r->datid[i],"vd")==0){
      if (!c->savd)
	mylog_exit(mylogf,"*** POP_UTIL_RESPONSE_HANDOVER Vd not saved\n");

      // Sampling is sec for Vd
      mod_util_resp_store_f_reg(r,i,c->vmt,c->vmd,c->vn,1.0,mylogf);
    }else if ((strcmp(r->datid[i],"si_mask")==0) ||
	      (strcmp(r->datid[i],"si_mask1")==0) ||
	      (strcmp(r->datid[i],"si_mask2")==0)){
      //printf("pop_util.c    HERE wyeth \n");
      d = pop_cell_si_get_data_ptr(mylogf,syn,r->datid[i],&dn);
      
      // Data is regular, sampled at the SI circbuf dt, which is
      // set when the SI is created
      mod_util_resp_store_f(r,i,d,dn,1,mylogf);

    }else{
      //
      //  Check for this datid in csavs
      //
      // WYETH - COULD USE HERE:  int popc_csav_get_csi(c,datid)
      //
      csav_flag = 0;
      if (c->csav_cnt != NULL){  // If this cell has 'csav' data
	csi = search_2d_carray(c->pl->csav_datid,r->datid[i],c->pl->csav_n);
	//printf("ID: %s  N: %d\n",r->datid[i],c->pl->csav_n);
	if (csi >= 0){  // We found the 'datid' in the layer's csav list
	  //printf(" found, csi = %d\n",csi);
	  d  = c->csav_f[csi];   // Pointer to the data
	  dn = c->csav_cnt[csi];
	  if (d != NULL){
	    if (r->samp[i] != c->pl->csav_samp[csi])
	      mylogx(mylogf,"POP_UTIL_RESPONSE_HANDOVER",
		     "Sampling mismatch for csav\n");
	    mod_util_resp_store_f(r,i,d,dn,1,mylogf);  // 1-Copy
	    csav_flag = 1;
	  }else{
	    mylogx(mylogf,"POP_UTIL_RESPONSE_HANDOVER","data is NULL");
	  }
	}else{
	  ; // Not found
	}
      }

      if (csav_flag == 0){
	sprintf(tstr,"  datid:  %s\n",r->datid[i]);
	mylog(mylogf,tstr);
	mylogx(mylogf,"POP_UTIL_RESPONSE_HANDOVER","Unknown data ID");
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_UTIL_CELL_SPIKE_INC                        */
/*                                                                           */
/*****************************************************************************/
void pop_util_cell_spike_inc(c,m)
     struct pop_cell *c;
     int m; // Increase storage by this factor
{
  int i;
  int n,nn;
  float *s,*t;

  n = c->maxscnt;
  nn = m*n;
  printf("  Spike storage limit being increased from %d to %d\n",n,nn);

  s = (float *)myalloc(nn*sizeof(float));
  t = c->s;
  for(i=0;i<n;i++)
    s[i] = t[i];
  myfree(c->s);
  c->s = s;
  c->maxscnt = nn;
}
/**************************************-**************************************/
/*                                                                           */
/*                           POP_UTIL_ADD_SPIKE_TO_CB                        */
/*                                                                           */
/*****************************************************************************/
void pop_util_add_spike_to_cb(syn,t,amp)
     struct pop_syn *syn;  // Synapse, for access to post-syn cell & syn type
     float t;              // Time to add spike
     float amp;            // Spike amplitude
{
  int j,k;
  int n;
  struct pop_circbuf *cb;

  // Choose the appropriate circular buffer, "cb"
  if (syn->stype == 1){
    cb = syn->post->ex1s;
  }else if (syn->stype == 2){
    cb = syn->post->in1s;
  }else{
    cb = NULL;
    exit_error("POP_UTIL_ADD_SPIKE_TO_CB","Unknown synapse type");
  }

  // Determine the index 'j' in CB where spike should go
  n = cb->n;
  if (t < cb->t){
    printf("syn->post  %s\n",syn->post->name);
    printf("syn->pre   %s\n",syn->pre->name);
    printf("t = %f\n",t);
    printf("cb->t = %f\n",cb->t);
    exit_error("POP_UTIL_ADD_SPIKE_TO_CB","Spike arrives in past");
  }
  k = (int)((t - cb->t)/cb->dt); // Number of bins from origin
  if (k >= n){
    printf("  cb->t:              %f \n",cb->t);
    printf("  cb->dt:             %f \n",cb->dt);
    printf("  cb->n:              %d \n",cb->n);
    printf("  Time to add spike:  %f \n",t);
    exit_error("POP_UTIL_ADD_SPIKE_TO_CB","Spike arrives in future");
  }
  j = cb->i + k; // Wrap around if necessary
  if (j >= n)
    j -= n;

  // Add the spike to the chosen buffer
  cb->d[j] += amp;  // The value added controls strength of input
}
/**************************************-**************************************/
/*                                                                           */
/*                          POP_UTIL_CIRCBUF_VAL_AT_T                        */
/*                                                                           */
/*  Return the value at time 't' from the circbuf, or 0.0 if the time is     */
/*  in the future beyond the end of the circbuf data.                        */
/*                                                                           */
/*****************************************************************************/
float pop_util_circbuf_val_at_t(cb,t)
     struct pop_circbuf *cb;  // circular buffer
     double t;                // time at which to sample value
{
  int cbn,nmask,npast,dti,ti0;
  float fval;
  double dt;

  cbn  = cb->n;         // Length of circ buf
  nmask = cb->maskn;    // Number of points in mask (future)
  npast = cbn - nmask;  // Number of points in past

  //  Compute time index within circbuf
  //    dt    - Time step from current origin to 't'
  //    dti   - 'dt' expressed in terms of cb elements
  //    ti0   - cb index of 't'
  //
  dt  = t - cb->t;   // Time between current origin and new mask
  if (dt < 0.0){

    dti = (int)(-0.5 + dt/cb->dt);  // Round negative value to nearest int

    if (dti < -npast){
      printf("  negative dt  %.6f   npast = %d\n",dt,npast);
      exit_error("POP_UTIL_CIRCBUF_VAL_AT_T","Negative dt, not allowed");
    }
  }else{
    dti = (int)(0.5 + dt/cb->dt);  // Round a non-negative value to nearest int
  }

  if (dti >= nmask){  // After any data in circbuf, thus 0.0

    fval = 0.0;  // We have jumped the entire buffer, return zero

  }else{
    ti0 = cb->i + dti;
    if (ti0 >= cbn)  // Wrap around if needed
      ti0 -= cbn;

    if (ti0 < 0)     // Could wrap backwards
      ti0 += cbn;

    fval = cb->d[ti0];
  }

  return fval;
}
/**************************************-**************************************/
/*                                                                           */
/*                       POP_UTIL_CIRCBUF_VAL_AT_T_INTERP                    */
/*                                                                           */
/*  Return the value at time 't' from the circbuf, or 0.0 if the time is     */
/*  in the future beyond the end of the circbuf data.                        */
/*                                                                           */
/*****************************************************************************/
float pop_util_circbuf_val_at_t_interp(cb,t)
     struct pop_circbuf *cb;  // circular buffer
     double t;                // time at which to sample value
{
  int cbn,nmask,npast,dti,ti0,ti1;
  float fval0,fval1,fval,tfrac;
  double dt,cbdt;

  cbn  = cb->n;         // Length of circ buf
  nmask = cb->maskn;    // Number of points in mask (future)
  npast = cbn - nmask;  // Number of points in past
  cbdt = cb->dt;        // Circbuf DT

  //  Compute time index within circbuf
  //    dt    - Time step from current origin to 't'
  //    dti   - 'dt' expressed in terms of cb elements
  //    ti0   - cb index of 't'
  //
  dt = t - cb->t;   // Time between current origin and new mask
  if (dt < 0.0){

    exit_error("POP_UTIL_CIRCBUF_VAL_AT_T_INTERP","Negative Not imp'd");
    // WYETH - FIX THIS
    // WYETH - FIX THIS
    // WYETH - FIX THIS
    // WYETH - FIX THIS
    // WYETH - FIX THIS


    dti = (int)(-0.5 + dt/cbdt);  // Round negative value to nearest int

    if (dti < -npast){
      printf("  negative dt  %.6f   npast = %d\n",dt,npast);
      exit_error("POP_UTIL_CIRCBUF_VAL_AT_T_INTERP","Negative dt not allowed");
    }
  }else{
  //dti = (int)(0.5 + dt/cb->dt);  // Round a non-negative value to nearest int
    dti = (int)(dt/cbdt);  // Index at or before 't'
    tfrac = (dt - dti*cbdt)/cbdt;      // Fractional part of cb->dt
  }

  if (dti >= nmask){  // After any data in circbuf, thus 0.0

    fval = 0.0;  // We have jumped the entire buffer, return zero

    //printf("  dt = %f   dti = %d  nmask = %d\n",dt,dti,nmask);

  }else{
    ti0 = cb->i + dti;
    if (ti0 >= cbn)  // Wrap around if needed
      ti0 -= cbn;

    if (ti0 < 0)     // Could wrap backwards
      ti0 += cbn;

    fval0 = cb->d[ti0];

    if (dti == (nmask-1)){  // If we were at the end of the mask
      fval1 = 0.0;
    }else{
      ti1 = ti0 + 1;
      if (ti1 >= cbn)  // Wrap around if needed
	ti1 = 0;
      fval1 = cb->d[ti1];
    }

    fval = fval0 + tfrac * (fval1-fval0);  // Interpolate

    /*
    if (fval != 0.0)
      printf("__________fval = %f\n",fval);
    */

  }

  return fval;
}
/**************************************-**************************************/
/*                                                                           */
/*                        POP_UTIL_CIRCBUF_MASK_UPDATE_T                     */
/*                                                                           */
/*  All masks written *must* be of the same length, 'maskn', otherwise the   */
/*  erasing may fail.                                                        */
/*                                                                           */
/*  To be TRUE at ALL times:                                                 */
/*  - cb->i  is the location at which the most recent mask was added         */
/*  - maskn <= cbn                                                           */
/*  - if maskn < cbn, then there should be a segment equal to the length     */
/*    difference that sits right *before* cb->i that represents the past.    */
/*  - The current future is only valid from cb->i to cb->i + maskn - 1.      */
/*                                                                           */
/*****************************************************************************/
void pop_util_circbuf_mask_update_t(cb,t,mask,maskn,sav,rsavti)
     struct pop_circbuf *cb;  // circular buffer
     double t;                // time to start writing data
     float *mask;             // mask data
     int    maskn;            // mask duration
     float *sav;              // save old values to long-term storage
     int   *rsavti;           // The last time at which data is saved
{
  int i,k;
  int cbn,cbi,dti,ti0,ti1;
  int erase_i0,erase_i1;
  float *cbd;
  double cbt,cbdt,dt;
  int savi;

  if (cb == NULL)
    exit_error("POP_UTIL_CIRCBUF_MASK_UPDATE_T","Cirbuf is NULL");

  cbn  = cb->n;    // Length of circ buf
  cbi  = cb->i;    // Current index
  cbt  = cb->t;    // Current time
  cbdt = cb->dt;   // DT
  cbd  = cb->d;    // buffer data

  //
  //  1.  Compute time at which new mask will be written
  //         dt    - Time step from current origin to new mask
  //         dti   - 'dt' expressed in terms of cb elements
  //         ti0   - cb index of origin of new mask
  //
  dt  = t - cbt;   // Time between current origin and new mask
  if (dt < 0.0){
    printf("  t %f   cbt %f\n",t,cbt);
    exit_error("POP_UTIL_CIRCBUF_MASK_UPDATE_T","Negative dt, not allowed");
  }
  dti = (int)(0.5 + dt/cbdt);  // Round a non-negative value to nearest int

  //printf("  pop_util:  dti = %d  (cbn %d)\n",dti,cbn);


  //
  //   1.5  Save for long term storage from 'cbi' to 'savt1'
  //
  *rsavti = -1;  // Assume no saving
  if (sav != NULL){
    savi = (int)(0.5 + t/cbdt);
    //printf("savi = %d      t= %f\n",savi,t);
  }

  if (dti >= cbn){
    ti0 = 0;  // We have jumped the entire buffer, essentially we reset

    for(i=0;i<maskn;i++){     // Write the mask into the beginning
      cbd[i] = mask[i];
      if (sav != NULL) sav[savi+i] = cbd[i];
    }
    for(i=maskn;i<cbn;i++)   // Clear the rest of the buffer
      cbd[i] = 0.0;

    if (sav != NULL)
      *rsavti = savi + maskn-1;

    // Now, we only need to do Step 4, below.

  }else{
    ti0 = cbi + dti;
    if (ti0 >= cbn)  // Wrap around if needed
      ti0 -= cbn;


    //
    //  2.  Erase segment in future, under assumption that 'maskn' is constant
    //
    erase_i0 = cbi + maskn;
    if (erase_i0 >= cbn)
      erase_i0 -= cbn;

    erase_i1 = erase_i0 + dti - 1;
    if (erase_i1 >= cbn)
      erase_i1 -= cbn;

    if (erase_i0 <= erase_i1){
      for(i=erase_i0;i<=erase_i1;i++)   // No wrap-around during erase
	cbd[i] = 0.0;
    }else{
      for(i=erase_i0;i<cbn;i++)   // Wrap-around during erase, up to end
	cbd[i] = 0.0;
      for(i=0;i<=erase_i1;i++)    // from beginning onwards
	cbd[i] = 0.0;
    }

    
    //
    //  3.  Add the mask into the buffer
    //
    ti1 = ti0 + maskn - 1;
    if (ti1 >= cbn)
      ti1 -= cbn;

    k = 0;
    if (ti0 <= ti1){
      for(i=ti0;i<=ti1;i++){   // No wrap-around during erase
	cbd[i] += mask[k];
	if (sav != NULL) sav[savi+k] = cbd[i];
	k += 1;
      }
    }else{
      for(i=ti0;i<cbn;i++){   // Wrap-around during erase, up to end
	cbd[i] += mask[k];
	if (sav != NULL) sav[savi+k] = cbd[i];
	k += 1;
      }
      for(i=0;i<=ti1;i++){    // from beginning onwards
	cbd[i] += mask[k];
	if (sav != NULL) sav[savi+k] = cbd[i];
	k += 1;
      }
    }

    if (sav != NULL){
      if (k > 0)
	*rsavti = savi + k-1;
    }
  }


  //
  //  4.  Update current values
  //
  cb->t = cbt + (double)dti*cbdt;   // Time moves only in 'cbdt' increments
  cb->i = ti0;

}
/**************************************-**************************************/
/*                                                                           */
/*                         POP_UTIL_PROCESS_SID_DS02                         */
/*                                                                           */
/*****************************************************************************/
void pop_util_process_sid_ds02(syn,t,amp)
     struct pop_syn *syn;  // Synapse
     float t;              // Time of spike at post-syn
     float amp;            // Spike amplitude
{
  int i,j,k;
  int ci,n,maskn,dsavn,savti;
  float w,*mask,*dsav;
  struct pop_mech *msi;
  struct pop_circbuf *cb;

  // Helpful pointers
  ci = syn->si->ci;               // Cell index: 0-send spike,  1-do mask
  msi = syn->si->siu->s001->msi;  // SI Def
  cb = syn->si->siu->s001->wt;    // circular buffer
  n = cb->n;                      // Length of circbuf

  // Calculate the index 'j' into the CB, set to -1 as over-run flag
  k = (int)((t - cb->t)/cb->dt); // Number of bins from origin

  //if (1){
  if (k < 0){
    printf("*** *** WARNING:\n");
    printf("  Pre-syn cell: %s\n",syn->pre->name);
    printf("  Post-syn cell: %s\n",syn->post->name);
    printf("  tlast: %f  alast: %f  tdelay: %f\n",syn->tlast,syn->alast,
	   syn->tdelay);
    printf("  ci = %d\n",ci);
    printf("  k = %d  cb->dt %f   cb->i %d\n",k,cb->dt,cb->i);
    printf("*** *** WARNING  spike is too early:  t=%f  cb->t = %f\n",t,cb->t);
    //exit_error("POP_UTIL_PROCESS_SID_DS02","Spike in past");
    printf("  Spike set to current time\n");
    printf("CONTINUING...\n");
    k = 0;  // Ignore it if spikes are too early
  }

  if (k >= n)
    j = -1;    // Signifies that we have jumped beyond the buffer
  else{
    j = cb->i + k; // Get cb index for this spike time, wrapping if needed
    if (j >= n)
      j -= n;
  }

  if (ci == 0){ // Get the weight and send spike through
    
    if (msi->si_comp == 2){  //  1-weight
      if (j < 0){
	w = 1.0;   // Beyond the CB, thus, no weight modification
      }else{
	w = 1.0 - cb->d[j];  // Apply the mask to the spike amplitude
	if (w < 0.0)
	  w = 0.0;           // Cannot have a negative weight
      }
    }else{  // default
      if (j < 0){
	w = 0.0;     // Beyond the CB, thus, no weight modification
      }else{
	w = cb->d[j];  // Apply the mask to the spike amplitude
	if (w < 0.0)
	  w = 0.0;           // Cannot have a negative weight
      }
    }

    // w is the multiplying factor
    pop_util_add_spike_to_cb(syn,t,w*amp);

  }else if (ci == 1){ // Modify weight mask in the circular buffer

    mask  = msi->mask1;  // Helpful pointers
    maskn = msi->mask1n;
    //d = cb->d;

    if (syn->si->siu->s001->sisv != NULL){
      dsav = pop_cell_si_get_data_ptr(NULL,syn,"si_mask",&dsavn);
      //ti = (int)(t/cb->dt); // Spike time in storage units
    }else
      dsav = NULL;


    //
    //  Write Mask into CircBuff, and possible long-term save
    //

    //printf("Pre-Post cell: %s %s\n",syn->pre->name,syn->post->name);

    // WYETH - 'dsav' must be long enough to account for 'maskn' after 't'
    // WYETH - 'dsav' must be long enough to account for 'maskn' after 't'
    // WYETH - 'dsav' must be long enough to account for 'maskn' after 't'
    pop_util_circbuf_mask_update_t(cb,t,mask,maskn,dsav,&savti);

    if (savti != -1){
      if (dsav != NULL){  // Update the number of saved data values
	pop_cell_si_sav_set_n(NULL,syn,"si_mask",savti);
      }
    }

  }else
    exit_error("POP_UTIL_PROCESS_SID_DS02","Bad ci value");

}
/**************************************-**************************************/
/*                                                                           */
/*                        POP_UTIL_PROCESS_SID_SYMMASK                       */
/*                                                                           */
/*  Uses 'pop_si002'.                                                        */
/*                                                                           */
/*****************************************************************************/
void pop_util_process_sid_symmask(syn,t,amp)
     struct pop_syn *syn;  // Synapse
     float t;              // Time of spike at post-syn
     float amp;            // Spike amplitude
{
  int i,k;
  int ti,ci,di,n,maskn,dsavn,jw,jm;
  float w,*mask,*d,*dsav;
  struct pop_mech *msi;
  struct pop_circbuf *cbm,*cbw;

  //
  //  ci 0  - lay mask in wt1, get wieght from wt2
  //  ci 1  - lay mask in wt2, get wieght from wt1
  //

  // Helpful pointers
  ci = syn->si->ci;               // Cell index: 0-send spike,  1-do mask
  msi = syn->si->siu->s002->msi;  // SI Def

  if (ci == 0){
    cbm = syn->si->siu->s002->wt1;   // circular buffer to lay mask
    cbw = syn->si->siu->s002->wt2;   // circular buffer for spike weight
  }else if (ci == 1){
    cbm = syn->si->siu->s002->wt2;   // circular buffer to lay mask
    cbw = syn->si->siu->s002->wt1;   // circular buffer for spike weight
  }else{
    exit_error("POP_UTIL_PROCESS_SID_SYMMASK","Bad ci value");
  }
  n = cbm->n;                        // Length of circbuf [Same for both]

  //
  //  Calculate index 'jm' into CB for laying mask, -1 as over-run flag
  //
  k = (int)((t - cbm->t)/cbm->dt); // Number of bins from origin
  if (k < 0){
    printf("k = %d  cbm->dt %f\n",k,cbm->dt);
    printf("*** WARNING  spike is too early:  t=%f  cbm->t = %f\n",t,cbm->t);
    exit_error("POP_UTIL_PROCESS_SID_SYMMASK","Spike in past");
    k = 0;  // Ignore it if spikes are too early
  }
  if (k >= n)
    jm = -1;    // Signifies that we have jumped beyond the buffer
  else{
    jm = cbm->i + k; // Get cb index for this spike time, wrapping if needed
    if (jm >= n)
      jm -= n;
  }

  //
  //  Calculate index 'jw' into CB for getting weight, -1 as over-run flag
  //
  k = (int)((t - cbw->t)/cbw->dt); // Number of bins from origin
  if (k < 0){
    printf("*** WARNING  spike is too early:  t=%f  cbw->t = %f\n",t,cbw->t);
    exit_error("POP_UTIL_PROCESS_SID_SYMMASK","Spike in past");
    k = 0;  // Ignore it if spikes are too early
  }
  if (k >= n)
    jw = -1;    // Signifies that we have jumped beyond the buffer
  else{
    jw = cbw->i + k; // Get cb index for this spike time, wrapping if needed
    if (jw >= n)
      jw -= n;
  }


  //
  //  Get the weight from 'cbw' and send spike through
  //
  if (msi->si_comp == 2){  //  1 - weight  ("one minus the weight")
    if (jw < 0){
      w = 1.0;   // Beyond the CB, thus, no weight modification
    }else{
      w = 1.0 - cbw->d[jw];  // Apply the mask to the spike amplitude
      if (w < 0.0)
	w = 0.0;           // Cannot have a negative weight
    }
    //printf("HERE WYETH    1 - mask           w = %f\n",w);
  }else{  // default
    if (jw < 0){
      w = 0.0;     // Beyond the CB, thus, no weight modification
    }else{
      w = cbw->d[jw];  // Apply the mask to the spike amplitude
      if (w < 0.0)
	w = 0.0;           // Cannot have a negative weight
    }
    //printf("HERE WYETH  default ************ w = %f\n",w);
  }
  pop_util_add_spike_to_cb(syn,t,w*amp); // w is the multiplying factor


  //
  //  Update the mask in 'cbm'
  //
  mask  = msi->mask1;  // Helpful pointers
  maskn = msi->mask1n;
  d = cbm->d;

  if (syn->si->siu->s002->sisv != NULL){

    //printf("HERE 0000 1\n");

    if (ci == 0){
      //printf("pop_util.c    HERE wyeth 22 \n");
      dsav = pop_cell_si_get_data_ptr(NULL,syn,"si_mask1",&dsavn);
    }else
      dsav = pop_cell_si_get_data_ptr(NULL,syn,"si_mask2",&dsavn);
    ti = (int)(t/cbm->dt); // Spike time in storage units

    //printf("----------SHOULD save this spike\n");

  }else
    dsav = NULL;
  
  // If we jumped over the entire buffer, start at beginning of CB
  if (jm < 0){
    for(i=0;i<maskn;i++){
      d[i] = mask[i];
      if (dsav != NULL){ //  Bug fixed here 26/Sept/2008 - this added
	dsav[ti + i] = d[i];   // Save the mask data for response
      }
    }
    while(i < n){
      d[i] = 0.0;  // Use zeros for the last (two) places
      if (dsav != NULL){ // Bug fixed here 26/Sept/2008 - this added
	dsav[ti + i] = d[i];   // Save the mask data for response
      }
      i += 1;
    }
    cbm->i = 0;   // Circ buf index is zero
    cbm->t = t;   // Time at circ buf index
  }else{
    i = cbm->i; // Current index in CB
    while(i != jm){
      d[i] = 0.0; // Clear
      i += 1;
      if (i >= n)
	i -= n;
    }
    // Update CB index to index of current spike
    di = jm;
    for(i=0;i<maskn;i++){
      d[di] += mask[i];
      
      if (dsav != NULL){  // Save the mask data for response
	dsav[ti + i] = d[di];  
      }
      
      di += 1;
      if (di >= n)
	di -= n;
    }
    
    cbm->i = jm;
    cbm->t = t;
  }

  // Update the number of saved data values
  if (dsav != NULL){
    //printf("   ti+maskn = %d\n",ti+maskn);
    if (ci == 0)
      pop_cell_si_sav_set_n(NULL,syn,"si_mask1",ti+maskn);
    else
      pop_cell_si_sav_set_n(NULL,syn,"si_mask2",ti+maskn);
  }
  
}
/**************************************-**************************************/
/*                                                                           */
/*                          POP_UTIL_PROCESS_SID_MULT                        */
/*                                                                           */
/*  Uses 'pop_si002'.                                                        */
/*                                                                           */
/*****************************************************************************/
void pop_util_process_sid_mult(syn,t,amp)
     struct pop_syn *syn;  // Synapse
     float t;              // Time of spike at post-syn
     float amp;            // Spike amplitude
{
  int ci,maskn,dsavn,savti;
  float *mask,*dsav;
  struct pop_mech *msi;
  struct pop_circbuf *cbm;

  ci = syn->si->ci;               // Cell index: 0,1
  msi = syn->si->siu->s002->msi;  // SI Def

  //  For cell index ('ci') 0,1 put mask in circ buf wt1,wt2
  if (ci == 0){
    cbm = syn->si->siu->s002->wt1;   // circular buffer to lay mask
  }else if (ci == 1){
    cbm = syn->si->siu->s002->wt2;   // circular buffer to lay mask
  }else{
    exit_error("POP_UTIL_PROCESS_SID_MULT","Bad ci value");
  }

  mask  = msi->mask1;  // Helpful pointers
  maskn = msi->mask1n;

  //
  //  Prepare for saving if masks are requested in .rsp file
  //
  if (syn->si->siu->s002->sisv != NULL){

    // WYETH - it appears that both 'si_mask1' and 'si_mask2' must be saved
    // if either is???  Thus, in .rsp files, I must request both masks, or
    // neither??? 

    //  To save the SI mask data, we need a pointer to the destination.
    if (ci == 0){
      dsav = pop_cell_si_get_data_ptr(NULL,syn,"si_mask1",&dsavn);
    }else{
      dsav = pop_cell_si_get_data_ptr(NULL,syn,"si_mask2",&dsavn);
    }
  }else{
    dsav = NULL;
  }

  //
  //  Write Mask into CircBuff, and possible long-term save
  //
  // WYETH - 'dsav' must be long enough to account for 'maskn' after 't'
  // WYETH - 'dsav' must be long enough to account for 'maskn' after 't'
  // WYETH - 'dsav' must be long enough to account for 'maskn' after 't'
  pop_util_circbuf_mask_update_t(cbm,t,mask,maskn,dsav,&savti);

  if (savti != -1){
    if (dsav != NULL){  // Update the number of saved data values
      if (ci == 0)
	pop_cell_si_sav_set_n(NULL,syn,"si_mask1",savti);
      else
	pop_cell_si_sav_set_n(NULL,syn,"si_mask2",savti);
    }
  }

  // NOTE:  Do not add any spikes to post-syn cells.  This will be handled
  // in routines such as 'pop_input_01'.
}
/**************************************-**************************************/
/*                                                                           */
/*                           POP_UTIL_PROCESS_SI001                          */
/*                                                                           */
/*****************************************************************************/
void pop_util_process_si001(syn,msi,tgen,amp)
     struct pop_syn *syn;
     struct pop_mech *msi;
     float tgen;           // Time the spike was generated
     float amp;            // Spike amplitude (e.g., after any syn depr)
{
  if (strcmp(msi->type,"si_ds02")==0){
    pop_util_process_sid_ds02(syn,tgen,amp);
  }else if ((strcmp(msi->type,"si_pair_rfdist")==0) ||
	    (strcmp(msi->type,"si_nonlin")==0)){
    if (onode_test_chr(msi->o,"nonlin","mask")){
      pop_util_process_sid_ds02(syn,tgen,amp);
    }else if (onode_test_chr(msi->o,"nonlin","symmask")){
      pop_util_process_sid_symmask(syn,tgen,amp);
    }else if (onode_test_chr(msi->o,"nonlin","mult")||
	      onode_test_chr(msi->o,"nonlin","divide")){
      pop_util_process_sid_mult(syn,tgen,amp);
    }else{
      printf("***  nonlin value not imp'd yet\n");
      exit(0);
    }
  }else{
    exit_error("POP_UTIL_PROCESS_SI001","Unknown SID name");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          POP_UTIL_PROCESS_SI_SPIKE                        */
/*                                                                           */
/*****************************************************************************/
void pop_util_process_si_spike(syn,tgen,amp)
     struct pop_syn *syn;
     float tgen;           // Time the spike was generated
     float amp;            // Spike amplitude (e.g., after any syn depr)
{
  struct pop_mech *msi;

  if (syn->si->siui == 1){

    msi = syn->si->siu->s001->msi;
    pop_util_process_si001(syn,msi,tgen,amp);

  }else if (syn->si->siui == 2){

    msi = syn->si->siu->s002->msi;
    pop_util_process_si001(syn,msi,tgen,amp);  // NOTE, using this for 002

  }else{
    exit_error("POP_UTIL_PROCESS_SI_SPIKE","Unknown SI Union Index");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                       POP_UTIL_PROCESS_SPIKE_POSTSYN                      */
/*                                                                           */
/*  Spikes in [t , t+dt) are put in the bin starting at time "t".            */
/*                                                                           */
/*****************************************************************************/
void pop_util_process_spike_postsyn(syn,tgen,sdf,sdtau)
     struct pop_syn *syn;
     float tgen;           // Time the spike was generated
     float sdf,sdtau;      // For synaptic depression
{
  float t,w,dt,anew;

  // Add appropriate delay and jitter to spike time
  t = tgen + 0.0015;  // Add at least 1ms to any spike
  t += syn->tdelay;   // For example, distance propagation delay
  
  // DEBUG
  /*if (strcmp(syn->pre->name,"in_surr_13_10_0")==0){
    printf("syn->tdelay = %f\n",syn->tdelay);
    }*/

  // Synaptic Depression
  if (sdf < 1.0){
    dt = t - syn->tlast;
    anew = 1.0 - (1.0 - syn->alast*sdf) * exp(-dt/sdtau);
    if (( anew < 0.0) || (anew > 1.0)){
      printf("anew %f  dt %f  (%f - %f)\n",anew,dt,t,syn->tlast);
    }
    syn->tlast = t;
    syn->alast = anew;

    // WYETH HERE
    //printf("%s to %s   anew = %f\n",syn->pre->name,syn->post->name,anew);
  }else{
    anew = 1.0;

    //
    //  *** WYETH - added for debugging, should REMOVE 2010 Feb 09
    //  *** WYETH - added for debugging, should REMOVE
    //  *** WYETH - added for debugging, should REMOVE
    //
    /***
    if (t < syn->tlast){
      printf(" *** WARNING: t     = %f\n",t);
      printf("              tlast = %f\n",syn->tlast);
      //exit_error("POP_UTIL_PROCESS_SPIKE_POSTSYN","spike times reversed");
    }
    if (t-syn->tlast < 0.00245)
      printf("  dt:  %f\n",t-syn->tlast);
      
    syn->tlast = t;
    ***/
  }

  w = syn->w; // Weight

  /*
if (strcmp(syn->pre->name,"exs_11_11_0")==0){
  printf("****************SPIKING 11 11 at time %f\n",tgen);
}else if (strcmp(syn->pre->name,"exs_17_11_0")==0){
  printf("****************SPIKn   17 11 at time %f\n",tgen);
}
  */

  // Handle any special Synaptic Interaction (SI), or do regular add
  if (syn->si != NULL){
    //printf("------------------HERE si\n");
    pop_util_process_si_spike(syn,t,(anew*w));
  }else{
    //printf("------------------HERE not si\n");
    pop_util_add_spike_to_cb(syn,t,(anew*w));
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                          POP_UTIL_PAIR_GET_PARAMS                         */
/*                                                                           */
/*****************************************************************************/
void pop_util_pair_get_params(o,rdir,rdist,rjitp,rjito,rseed)
     struct onode *o;
     float *rdir,*rdist,*rjitp,*rjito;
     int *rseed;
{
  int seed;
  float dir,dist,jitp,jito;

  //WYETH HERE;  UNUSED
  //WYETH HERE;  UNUSED

  dir   = onode_getpar_flt_exit(o,"offset_direction");   // Vector dir.
  dist  = onode_getpar_flt_exit(o,"offset_distance");    // Vector length
  jitp  = onode_getpar_flt_exit(o,"offset_jitter_para"); // parallel jit.
  jito  = onode_getpar_flt_exit(o,"offset_jitter_orth"); // orthog jitter
  seed  = onode_getpar_int_exit(o,"offset_jitter_seed"); // jitter seed
}
/**************************************-**************************************/
/*                                                                           */
/*                          POP_UTIL_MECH_CONFIG_MASK                        */
/*                                                                           */
/*  Configure mask parameters within a mechanism.                            */
/*                                                                           */
/*  Mask type                                                                */
/*    1 - sinusoid, down, up, down                                           */
/*    2 - maxwell, up to peak fast, down slowly, max 1 default.              */
/*    3 - sinusoid, positive lobe                                            */
/*    4 - half Gaussian                                                      */
/*    5 - alpha-function - has long tail.                                    */
/*                                                                           */
/*****************************************************************************/
void pop_util_mech_config_mask(mylogf,m)
     char *mylogf;
     struct pop_mech *m;           // Mechanism
{
  int i;
  int wtn,masktype,maskn,npast;
  float *mask,mc,tt,sigma,wdt,dtms,dt,maskamp,alpha;
  char *mask_comp,*nonlin;

  mylog(mylogf,"  POP_UTIL_MECH_CONFIG_MASK\n");

  dt       = onode_getpar_flt_exit(m->o,"dt");
  masktype = onode_getpar_int_exit(m->o,"masktype");
  maskamp  = onode_getpar_flt_dflt(m->o,"mask_amp",1.0);

  // ******************************************************************
  // WYETH - 31 March 2010: Partially fixed problem here: if I changed
  // the value below (wdt = 0.001) to 0.0002, then I used to get:
  //   ...
  //   negative dt  -0.000800   npast = 2
  //   *** POP_UTIL_CIRCBUF_VAL_AT_T:  Negative dt, not allowed.  Exiting.
  //
  // And then, I might have fixed the above, but 'mod_mon_..' doesn't
  //   plot the mask properly if it is not at 1 msec dt
  //
  // Thus, I believe the problem is that some parts of the code (mod_mon_)
  //   may be assuming a 1ms time scale of the "mask..." circbuf data. 
  //
  wdt = 0.001; // *** SEE NOTES ABOVE  // circular buffer time resolution

  //dtms = dt * 1000.0;  // dt in msec  WYETH - OLD

  dtms = my_rint(dt/wdt);  // dt circbuf time units

  npast = 0.002/wdt;  // Number of points "in the past", 2ms worth of past

  /*
  if (masktype == 1)
    wtn = (int)(2.0*dtms + 2.0);   // circular buffer length in time units
  else if (masktype == 2)
    wtn = (int)(4.0*dtms + 2.0);
  else if (masktype == 3)
    wtn = (int)(2.0*dtms + 2.0);
  else if (masktype == 4)
    wtn = (int)(3.0*dtms + 2.0);
  else if (masktype == 5)
    wtn = (int)(7.0*dtms + 2.0);
  */

  if (masktype == 1)
    wtn = (int)(2.0*dtms);   // circ buffer length in time units
  else if (masktype == 2)
    wtn = (int)(4.0*dtms);
  else if (masktype == 3)
    wtn = (int)(2.0*dtms);
  else if (masktype == 4)
    wtn = (int)(3.0*dtms);
  else if (masktype == 5)
    wtn = (int)(7.0*dtms);

  wtn += npast;

  m->wdt = wdt;
  m->wtn = wtn;

  //maskn = wtn - 2; // Mask is two units shorter than CircBuff
  maskn = wtn - npast; // Mask is two units shorter than CircBuff

  if (masktype == 1){ // SINEWAVE, peaks in middle of mask
    // Build the mask: it has amplitude 1, peaking in the middle at
    // time 'dt' and has negative lobes on either side.  It is over-all
    // negative, thus suppressive.
    mask = (float *)myalloc(maskn*sizeof(float));
    mc = 2.0*M_PI * 3.0/(4.0*dtms);
    for(i=0;i<maskn;i++){
      tt = (float)i;
      mask[i] = -maskamp*sin(mc*tt);
    }
  }else if (masktype == 2){ // MAXWELL
    sigma = dtms/sqrt(2.0);
    mask = maxwell_farray(0.0,sigma,1.0,maskn);
    make_max_const_farray(mask,maskn,maskamp);
  }else if (masktype == 3){ // SINE Positive lobe
    mask = (float *)myalloc(maskn*sizeof(float));
    mc = M_PI /(2.0*dtms);
    for(i=0;i<maskn;i++){
      tt = (float)i;
      mask[i] = maskamp * sin(mc*tt);
    }
  }else if (masktype == 4){ // Half Guassian
    mask = (float *)myalloc(maskn*sizeof(float));
    for(i=0;i<maskn;i++){
      tt = (float)i;
      mask[i] = maskamp * func_gaussian_one(tt,0.0,dtms);
    }
  }else if (masktype == 5){ // ALPHA
    alpha = 1.0/dtms; // alpha is 1/peaktime
    mask = alpha_farray(0.0,alpha,1.0,maskn);
    make_max_const_farray(mask,maskn,maskamp);
  }else
    mylog_exit(mylogf,"POP_UTIL_MECH_CONFIG_MASK  Bad masktype value.");

  m->mask1n = maskn;
  m->mask1 = mask;

  /*
  append_farray_plot("zzz.mask.pl","mask",mask,maskn,1);
  printf("mask dtms = %f  masktype = %d\n",dtms,masktype);
  //exit_error("WYETH------","done");
  */


  //
  //  Type of mask computation
  //
  nonlin = onode_getpar_chr_dflt(m->o,"nonlin","NULL"); // non-linearity
  if (strcmp(nonlin,"mult")==0){
    m->si_comp = 1;
  }else if (strcmp(nonlin,"divide")==0){
    m->si_comp = 2;
  }else{
    mask_comp  = onode_getpar_chr_dflt(m->o,"mask_comp","prod");
    if (strcmp(mask_comp,"prod")==0)
      m->si_comp = 1;
    else if (strcmp(mask_comp,"prod1m")==0){
      m->si_comp = 2;
    }else
      exit_error("POP_UTIL_MECH_CONFIG_MASK","Unknown 'mask_comp' value");
    myfree(mask_comp);
  }
  myfree(nonlin);

  //
  //  Constants
  //
  m->a = onode_getpar_flt_dflt(m->o,"c_a",0.0); // additive constant
  m->b = onode_getpar_flt_dflt(m->o,"c_b",1.0); // additive constant
}
/**************************************-**************************************/
/*                                                                           */
/*                             POP_UTIL_INPUT_PAIR                           */
/*                                                                           */
/*  Create and return a structure to hold parameters and pointers for        */
/*  making a paired input.                                                   */
/*                                                                           */
/*****************************************************************************/
struct input_pair *pop_util_input_pair(mpt,lay,io,inindex)
     struct pop_top *mpt;     // Model top
     struct pop_layer *lay;   // Layer containing this <input>
     struct onode *io;        // <input> which contains the <input_pair>
     int inindex;             // Input index (unique reference)
{
  char *lf;
  struct onode *od,*oip;
  struct input_pair *p;

  lf = mpt->logf;

  mylog(lf,"  POP_UTIL_INPUT_PAIR\n");

  oip = onode_child_get_unique(io,"input_pair");  // Onode for <input_pair>
  if (oip == NULL)
    mylogx(lf,"POP_UTIL_INPUT_PAIR","Cannot find <input_pair>");
    

  // 'p' will be filled and returned
  p = (struct input_pair *)myalloc(sizeof(struct input_pair));
  p->pair_seed = -1;   // Default value
  p->tseed = -1;       // Default value
  p->in_index = inindex;
  p->distrib = NULL;
  p->z1 = NULL;

  //
  //  Pairing rules, distribution parameters
  //
  od = onode_child_get_unique(oip,"pair_distrib");
  if (od != NULL){
    p->distrib = onode_getpar_chr_exit(od,"type");
    if ((strcmp(p->distrib,"rf_dist")==0) ||
	(strcmp(p->distrib,"rf_dist_phase")==0)){
      p->rf_offset_dir         = onode_getpar_flt_exit(od,"direction");
      p->rf_offset_dist        = onode_getpar_flt_exit(od,"distance");
      p->rf_offset_sd_par  = onode_getpar_flt_exit(od,"jitter_para");
      p->rf_offset_sd_orth = onode_getpar_flt_exit(od,"jitter_orth");
      p->pair_seed        = onode_getpar_int_exit(od,"jitter_seed");
      if (strcmp(p->distrib,"rf_dist_phase")==0){
	p->rf_ph_shift  = onode_getpar_flt_dflt(od,"ph_shift",90.0);
	p->rf_ph_dev    = onode_getpar_flt_dflt(od,"ph_dev",5.0);
      }else{
	p->rf_ph_shift  =  0.0;
	p->rf_ph_dev    = -1.0;
      }
      
    }else
      mylogx(lf,"POP_UTIL_INPUT_PAIR","Unknown pair_distrib type");
  }else{
    mylogx(lf,"POP_UTIL_INPUT_PAIR","Cannot find <pair_distrib>");
  }

  
  //
  //  Population names and z-layer specification
  //
  p->pop0 = popu_get_lay_named(mpt,io,"pop_origin","POP_UTIL_INPUT PAIR pop0");
  if (onode_test_ostr(oip,"pop_origin_2")){
    p->pop1 = popu_get_lay_named(mpt,oip,"pop_origin_2",
				 "POP_UTIL_INPUT PAIR pop1");
  }else{
    p->pop1 = p->pop0;
  }
  //p->z0 = onode_getpar_chr_dflt(oip,"pop_origin_1_z","all"); WYETH DELETE?
  p->z1 = onode_getpar_chr_dflt(oip,"pop_origin_2_z","all");

  //
  //  Synaptic Interaction
  //
  if (onode_test_ostr(oip,"syn_int")){
    if (onode_test_chr(oip,"syn_int","NULL") == 0){
      p->msi = popu_get_mech_named(mpt,lay,oip,"syn_int","MOD_POP_INPUT");
    }else
      p->msi = NULL;
  }

  /*
 printf("z0,1  %s %s\n",p->z0,p->z1);
  if (p->msi != NULL)
    printf("Mech name,type:  %s %s\n",p->msi->name,p->msi->type);
  else
    printf("Mech is NULL\n");
  exit(0);
  */

  return p;
}
/**************************************-**************************************/
/*                                                                           */
/*                        POP_UTIL_MAKE_CELL_STYLE_DS01                      */
/*                                                                           */
/*  NOT/NEVER USED??                                                         */
/*                                                                           */
/*****************************************************************************/
void pop_util_make_cell_style_ds01(mylogf,mppl,laylist,nlay,c)
     char *mylogf;
     struct param_pair_list *mppl; /* Model parameter pair list */
     struct pop_layer **laylist;   /* List of all layers */
     int nlay;                     /* Number of layers */
     struct pop_cell *c;           /* Cell to modify */
{
  int i,j;
  struct pop_layer *pl;

  mylog(mylogf,"  POP_UTIL_MAKE_CELL_STYLE_DS01\n");

  if (c->subclass != NULL)
    myfree(c->subclass);
  c->subclass = strdup("ds01");  /* For sub-cortical response processing */

  // Clear any existing inputs to this cell
  for(i=0;i<nlay;i++){
    pl = laylist[i];
    pop_cell_remove_all_syn_from_layer_to_cell(mylogf,pl,c);
  }

  // Create paired synaptic interactions from LGN
  {
    int *cx,*cy,seed,seed0,nsamp,maxrep,nposs,evflag,xn,yn,*tcx,*tcy,tcn;
    int *seeds;
    float sd_orth,sd_par,sdo,sdp,eps,sscale,xc,yc,orideg,eopar;
    char tstr[SLEN3],conname[SLEN];
    float dr,sd_noise;
    struct pop_lgn_in *ln;
    int tcn0,*tcx0,*tcy0,tcn1,*tcx1,*tcy1;

    strcpy(conname,"pop_ex_conn0"); /*** WYETH - use Gabor params for EX ***/
    sprintf(tstr,"%s_sd_orth",conname);
    sd_orth = paramfile_get_float_param_or_exit(mppl,tstr);
    sprintf(tstr,"%s_sd_par",conname);
    sd_par = paramfile_get_float_param_or_exit(mppl,tstr);
    sprintf(tstr,"%s_seed",conname);
    seed = paramfile_get_int_param_or_exit(mppl,tstr);
    sprintf(tstr,"%s_eps",conname);
    eps = paramfile_get_float_param_or_exit(mppl,tstr);
    sprintf(tstr,"%s_nsamp",conname);
    nsamp = paramfile_get_int_param_or_exit(mppl,tstr);
    sprintf(tstr,"%s_maxrep",conname);
    maxrep = paramfile_get_int_param_or_exit(mppl,tstr);

    sscale = paramfile_get_float_param_or_exit(mppl,"sscale");
    
    sdo = sd_orth / sscale;  /* Make units pix */
    sdp = sd_par / sscale;  /* Make units pix */

    i = c->layx;
    j = c->layy;
    xc = c->pl->x0 + (float)i*c->pl->xf;     /* Center for this location */
    yc = c->pl->y0 + (float)j*c->pl->yf;

    if (c->lgn_n > 0)
      ln = c->lgnin[0];
    else
      exit_error("POP_UTIL_MAKE_CELL_STYLE_DS01","No LGN inputs");

    xn = ln->lay->xn;  // LGN grid
    yn = ln->lay->yn;


    printf(" USE ATTRIB TO GET ORI\n");
    //orideg = c->ori;
    printf(" USE ATTRIB TO GET ORI\n");

    // Pick ON and OFF inputs independently, under a 2D Gaussian
    seeds = get_seeds(seed,100000,4);
    evflag = 1; // Even out connections
    eopar = 2.0;  // SD for even-out, used to be fixed in 'mod_conn_gauss_02'
    nposs = mod_conn_gauss_02(mylogf,xn,yn,xc,yc,sdo,sdp,orideg,seeds[0],
			      maxrep,eps,nsamp/2,evflag,eopar,0,0,
			      &tcx1,&tcy1,&tcn1);
    nposs = mod_conn_gauss_02(mylogf,xn,yn,xc,yc,sdo,sdp,orideg,seeds[1],
			      maxrep,eps,nsamp/2,evflag,eopar,0,0,
			      &tcx0,&tcy0,&tcn0);

    //
    // WYETH - In theory, this won't work, because if there are not already
    // LGN inputs, we'll get stopped several lines above, and if there are,
    // then we'll be adding new ones to the list.  A little thought can fix
    // this.
    exit_error("POP_UTIL_MAKE_CELL_STYLE_DS01","WYETH - BROKEN");
    //
    // *************
    // *************
    // *************
    //
    

    // NEEDS mechr - for post-syn receptor
    {
      struct pop_mech *mechr;  // WYETH - put here for compiler
      pop_cell_add_lgn_input(c,ln->lay,mechr,tcx0,tcy0,tcn0,tcx1,tcy1,tcn1,
			     nposs,-1.0);
    }

    
    // create matches for the inputs
    dr = 4.0;
    sd_noise = 0.0;
    mod_conn_pick_paired_connections(tcx1,tcy1,tcn1,xn,yn,orideg,dr,
				     sd_noise,seeds[2],&tcx,&tcy);

    append_iarray(&(ln->cx1),ln->cn1,tcx,ln->cn1);
    append_iarray(&(ln->cy1),ln->cn1,tcy,ln->cn1);
    myfree(tcx); myfree(tcy);
    ln->cn1 = 2*ln->cn1;

    mod_conn_pick_paired_connections(tcx0,tcy0,tcn0,xn,yn,orideg,dr,
				     sd_noise,seeds[3],&tcx,&tcy);
    append_iarray(&(ln->cx0),ln->cn0,tcx,ln->cn0);
    append_iarray(&(ln->cy0),ln->cn0,tcy,ln->cn0);
    myfree(tcx); myfree(tcy);
    ln->cn0 = 2*ln->cn0;

    myfree(seeds);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             POP_UTIL_CUSTOMIZE                            */
/*                                                                           */
/*  09May16 - changed to work w/ moo files.                                  */
/*                                                                           */
/*****************************************************************************/
void pop_util_customize(m,mpt)
     struct model_struct *m;     // Model params
     struct pop_top *mpt;        // Top
{
  int i,k;
  int x,y,z,layi,nlay,sii,cflag;
  char tstr[SLEN],*layname,*mylogf,*ctype;
  float dsdph,dsdev;
  struct param_pair *tp;
  struct pop_layer *l0,*l1;
  struct pop_cell *c0,*c1;
  struct pop_mech *msi;
  struct pop_layer **laylist;
  struct onode *t;

  mylogf = mpt->logf;

  mylog(mylogf,"  POP_UTIL_CUSTOMIZE\n");

  laylist = mpt->lay; // List of all layers
  nlay = mpt->nlay;   // Number of layers

  cflag = 0;  // Was any 'customize' requested?

  t = onode_get_next_type(m->o->o,"customize");

  //tp = mppl->c;
  while(t != NULL){
    cflag = 1;   // There is some customization

    ctype = onode_getpar_chr_exit(t,"type");  // Type of customization
    sprintf(tstr,"    %s\n",ctype);
    mylog(mylogf,tstr);

    if (strcmp(ctype,"synapse_delete")==0){
      
      if (onode_test_ostr(t,"unit_pre")){
	
	c0 = popu_get_c_named(mpt,t,"unit_pre","POP_UTIL_CUSTOMIZE u_pre");
	
	if (onode_test_ostr(t,"unit_post")){
	  c1 = popu_get_c_named(mpt,t,"unit_post","POP_UTIL_CUSTOMIZE u_post");


	  if (onode_test_ostr(t,"si_index")){  // Get the SI index
	    sii = onode_getpar_int_exit(t,"si_index");
	    printf("sii = %d\n",sii);
	  }

	  //
	  exit_error("POP_UTIL_CUSTOMIZE","Not Imp'd yet");
	  //

	}else if (onode_test_ostr(t,"layer_post")){
	  l1 = popu_get_lay_named(mpt,t,"layer_post","POP_UTIL_CUSTOMIZE lpo");

	  if (onode_test_ostr(t,"si_index")){  // Get the SI index
	    sii = onode_getpar_int_exit(t,"si_index");
	    pop_cell_remove_all_syn_from_cell_to_layer_si(mylogf,c0,l1,sii);
	  }else{
	    pop_cell_remove_all_syn_from_cell_to_layer(mylogf,c0,l1);
	  }
	}
	
      }else if (onode_test_ostr(t,"layer_pre")){
	
	l0 = popu_get_lay_named(mpt,t,"layer_pre","POP_UTIL_CUSTOMIZE l_pre");

	if (onode_test_ostr(t,"unit_post")){
	  c1 = popu_get_c_named(mpt,t,"unit_post","POP_UTIL_CUSTOMIZE u_post");
	  
	  //
	  exit_error("POP_UTIL_CUSTOMIZE","Not tested yet");
	  //

	  pop_cell_remove_all_syn_from_layer_to_cell(mylogf,l0,c1);
	}
      }else{
	mylog_exit(mylogf,"POP_UTIL_CUSTOMIZE  synapse_delete bad pre-syn\n");
      }
    }else if (strcmp(ctype,"synapse_adjust")==0){

      float synw;

      if (onode_test_ostr(t,"unit_pre")){
	
	c0 = popu_get_c_named(mpt,t,"unit_pre","POP_UTIL_CUSTOMIZE u_pre");
	
	if (onode_test_ostr(t,"unit_post")){
	  c1 = popu_get_c_named(mpt,t,"unit_post","POP_UTIL_CUSTOMIZE u_post");


	  if (onode_test_ostr(t,"si_index")){  // Get the SI index
	    sii = onode_getpar_int_exit(t,"si_index");
	    printf("sii = %d\n",sii);
	  }

	  //
	  exit_error("POP_UTIL_CUSTOMIZE","Not Imp'd yet");
	  //

	}else if (onode_test_ostr(t,"layer_post")){
	  l1 = popu_get_lay_named(mpt,t,"layer_post","POP_UTIL_CUSTOMIZE lpo");

	  synw = onode_getpar_flt_exit(t,"syn_weight");

	  pop_cell_adjust_all_syn_from_cell_to_layer(mylogf,c0,l1,synw);

	  printf("HERE WYETH w = %f\n",synw);
	  //exit(0);
	}
	
      }else if (onode_test_ostr(t,"layer_pre")){
	
	l0 = popu_get_lay_named(mpt,t,"layer_pre","POP_UTIL_CUSTOMIZE l_pre");

	if (onode_test_ostr(t,"unit_post")){
	  c1 = popu_get_c_named(mpt,t,"unit_post","POP_UTIL_CUSTOMIZE u_post");
	  
	  //
	  exit_error("POP_UTIL_CUSTOMIZE","Not tested yet");
	  //

	  pop_cell_remove_all_syn_from_layer_to_cell(mylogf,l0,c1);
	}
      }else{
	mylog_exit(mylogf,"POP_UTIL_CUSTOMIZE  synapse_delete bad pre-syn\n");
      }
    }else if (strcmp(ctype,"synapse_add")==0){

      float synw,tdelay;
      int syntype,inindex;

      if (onode_test_ostr(t,"unit_pre")){
	c0 = popu_get_c_named(mpt,t,"unit_pre","POP_UTIL_CUSTOMIZE u_pre");
	
	if (onode_test_ostr(t,"unit_post")){
	  c1 = popu_get_c_named(mpt,t,"unit_post","POP_UTIL_CUSTOMIZE u_post");

	  synw = onode_getpar_flt_exit(t,"syn_weight");
	  syntype = onode_getpar_flt_exit(t,"syn_type");
	  tdelay = onode_getpar_flt_dflt(t,"syn_tdelay",0.0);

	  inindex = -2;  // input index flag for customize, -2
	  pop_cell_add_synapse(c0,c1,syntype,synw,tdelay,inindex);

	  /**
	  if (onode_test_ostr(t,"si_index")){  // Get the SI index
	    sii = onode_getpar_int_exit(t,"si_index");
	    printf("sii = %d\n",sii);
	    }**/

	}else if (onode_test_ostr(t,"layer_post")){
	  exit_error("POP_UTIL_CUSTOMIZE","Not Imp'd yet");
	}
	
      }else if (onode_test_ostr(t,"layer_pre")){
	  exit_error("POP_UTIL_CUSTOMIZE","Not Imp'd yet");
      }else{
	mylog_exit(mylogf,"POP_UTIL_CUSTOMIZE  synapse_add bad pre-syn\n");
      }
    }else{
      mylog_exit(mylogf,"POP_UTIL_CUSTOMIZE  Unknown type of customize\n");
    }

    myfree(ctype);
    
    t = onode_get_next_type(t->next,"customize");
  }
  if (cflag == 1)
    mylog(mylogf,"    Done customize\n");
  else
    mylog(mylogf,"    No customization\n");

  
  /***

  // OLD Examples
  //   CUSTOMIZE unit ex_9_8_0 synapse_delete from_layer in
  //   CUSTOMIZE unit ex_9_8_0 synapse_delete to_layer in
  //   CUSTOMIZE unit ex_9_8_0 synapse_copy layer_to_unit in ex_8_8_0
  //   CUSTOMIZE unit ex_8_8_0 set_param bg_ex_rate 100.0
  //   CUSTOMIZE unit ex_9_8_0 set_param bg_in_rate 100.0
  //   CUSTOMIZE unit ex_8_8_0 bg_corr_ex 0.10 ex_9_8_0

     else if (strcmp(slist[2],"synapse_copy")==0){
      if (strcmp(slist[3],"layer_to_unit")==0){

	// Get pointer to layer
	k = pop_cell_get_layer_index(laylist,nlay,slist[4]);
	if (k < 0){
	  mylog_exit(mylogf,"POP_UTIL_CUSTOMIZE  layer name not found");
	}
	pl = laylist[k];

	pop_cell_parse_unit_name(slist[5],&layname,&x,&y,&z);
	layi = pop_cell_get_layer_index(laylist,nlay,layname);
	if (layi < 0){
	  mylog_exit(mylogf,"POP_UTIL_CUSTOMIZE  layer name not found");
	}
	c1 = &(laylist[layi]->c[x][y][z]);

	// Copy synapses of type pl -> c1 to cell 'c'
	pop_cell_copy_syn_layer_from_cell_to_cell(mylogf,pl,c1,c,0);
      }else{
	mylog_exit(mylogf,"POP_UTIL_CUSTOMIZE  synapse_copy type error");
      }
    }else if (strcmp(slist[2],"set_param")==0){
      if (strcmp(slist[3],"bg_ex_rate")==0){
	c->bg_x_rate = atof(slist[4]);
	sprintf(tstr,"      bg_ex_rate  %f\n",c->bg_x_rate);
	mylog(mylogf,tstr);
      }else if (strcmp(slist[3],"bg_in_rate")==0){
	c->bg_i_rate = atof(slist[4]);
	sprintf(tstr,"      bg_in_rate  %f\n",c->bg_i_rate);
	mylog(mylogf,tstr);
      }else{
	mylog_exit(mylogf,"POP_UTIL_CUSTOMIZE  set_param name error");
      }
    }else if (strcmp(slist[2],"bg_corr_ex")==0){
      pop_cell_parse_unit_name(slist[4],&layname,&x,&y,&z);
      layi = pop_cell_get_layer_index(laylist,nlay,layname);
      if (layi < 0){
	mylog_exit(mylogf,"POP_UTIL_CUSTOMIZE  layer name not found");
      }
      c1 = &(laylist[layi]->c[x][y][z]);

      if ((c->bg_x_corr_c != NULL)||(c1->bg_x_corr_c != NULL))
	mylog_exit(mylogf,"POP_UTIL_CUSTOMIZE  bg_corr: cell listed twice");
      
      c->bg_x_corr_f = atof(slist[3]);
      c1->bg_x_corr_f = atof(slist[3]);
      c->bg_x_corr_c = c1;
      c1->bg_x_corr_c = c;
    }else if (strcmp(slist[2],"bg_corr_in")==0){
      pop_cell_parse_unit_name(slist[4],&layname,&x,&y,&z);
      layi = pop_cell_get_layer_index(laylist,nlay,layname);
      if (layi < 0){
	mylog_exit(mylogf,"POP_UTIL_CUSTOMIZE  layer name not found");
      }
      c1 = &(laylist[layi]->c[x][y][z]);

      if ((c->bg_i_corr_c != NULL)||(c1->bg_i_corr_c != NULL))
	mylog_exit(mylogf,"POP_UTIL_CUSTOMIZE  bg_corr: cell listed twice");
      
      c->bg_i_corr_f = atof(slist[3]);
      c1->bg_i_corr_f = atof(slist[3]);
      c->bg_i_corr_c = c1;
      c1->bg_i_corr_c = c;
    }else if (strcmp(slist[2],"style")==0){
      if (strcmp(slist[3],"ds01")==0){
	pop_util_make_cell_style_ds01(mylogf,mppl,laylist,nlay,c);
      }else if (strcmp(slist[3],"ds02")==0){
	pl = pop_cell_get_layer_pointer(laylist,nlay,"ex");
	if (pl == NULL)
	  mylog_exit(mylogf,"POP_UTIL_CUSTOMIZE  Layer 'ex' not found.\n");

	msi = pop_cell_get_mech_by_name(pl,"ds02");
	// sid = pop_cell_layer_get_sid_pointer(pl,"ds02");

	if (msi == NULL){
	  float dsdt1;
	  int dsmask;

	  //mylog_exit(mylogf,"POP_UTIL_CUSTOMIZE  No SI Def for 'ds02'.\n");
	  mylog(mylogf,"  POP_UTIL_CUSTOMIZE  Defining SI for 'ds02'.\n");
	  mylog(mylogf,"  POP_UTIL_CUSTOMIZE  Using 'cmplx_conn_ds_dt1'.\n");
	  dsdt1 = paramfile_get_float_param_or_exit(mppl,"cmplx_conn_ds_dt1");
	  dsmask = paramfile_get_int_param_or_exit(mppl,"cmplx_conn_ds_mask");


	  exit_error("POP_UTIL_CUSTOMIZE","Changed format here");

	  //msi = mod_conn_def_si_ds02(mylogf,dsdt1,dsmask,1.0);
	  pop_cell_layer_add_mech(mylogf,c->pl,msi);
	}

	exit_error("POP_UTIL_CUSTOMIZE","Changed format for ds02 here");
	
	//
	//  WYETH - the phase shift and deviation need to be stored in 'msi'
	//  But, is this used any more?

	dsdph = paramfile_get_float_param_default(mppl,
						  "cmplx_conn_ds_ph_shift",
						  90.0);
	dsdev = paramfile_get_float_param_default(mppl,
						  "cmplx_conn_ds_ph_dev",5.0);

	mod_conn_syn_convert_cell(mylogf,pl,c,msi,-1);  // WYETH, was no inindx
      }else if (strcmp(slist[3],"gauss_output_1")==0){
	//  4 - output layer
	//  5 - SD
	//   6 - min weight
	//   7 - weight factor
	//   8 - autapse flag
	//   9 - synapse type

	float cdist,minw,wf;
	int self,syntype;

	pl = pop_cell_get_layer_pointer(laylist,nlay,slist[4]);
	cdist   = atof(slist[5]);
	minw    = atof(slist[6]);
	wf      = atof(slist[7]);
	self    = atoi(slist[8]);
	syntype = atoi(slist[9]);
	mod_conn_cell_to_layer_gauss(mylogf,c,pl,cdist,minw,wf,self,syntype,
				     -1);
      }else if (strcmp(slist[3],"gauss_input_od1")==0){
	// 4 - input layer
	//   5 - crit dist
	//   6 - crit ori
	//   7 - min weight
	//   8 - norm weight 
	//   9 - prob
	//   10 - seed
	//  11 - autapse flag
	//  12 - synapse type

	float cdist,cori,minw,wf,prob;
	int seed,self,syntype;

	pl = pop_cell_get_layer_pointer(laylist,nlay,slist[4]);
	cdist   = atof(slist[5]);
	cori    = atof(slist[6]);
	minw    = atof(slist[7]);
	wf      = atof(slist[8]);
	prob    = atof(slist[9]);
	seed    = atoi(slist[10]);
	self    = atoi(slist[11]);
	syntype = atoi(slist[12]);
	mod_conn_layer_to_cell_ori_dist_01(mylogf,pl,c,cdist,cori,minw,wf,
					   prob,seed,self,syntype,-1);
      }else{
	sprintf(tstr,"  *** style:  %s\n",slist[3]);
	mylog(mylogf,tstr);
	mylog_exit(mylogf,"POP_UTIL_CUSTOMIZE  Unknown style\n");
      }

    }else if (strcmp(slist[2],"isotropic_demo")==0){
      float wsd,wmin,wf;
      int sflag,syntype;
      struct pop_layer *plex,*plin;

      plex = pop_cell_get_layer_pointer(laylist,nlay,"ex");
      plin = pop_cell_get_layer_pointer(laylist,nlay,"in");

      if (strcmp("in",c->pl->name)==0){
	wf = 0.2;
	wsd =  70.0;
	wmin = 0.01;
	syntype = 2;
      }else{
	wf = 0.5;
	wsd = 140.0;
	wmin = 0.01;
	syntype = 1;
      }
      sflag = 0; // do not connect to self

      mod_conn_cell_to_layer_gauss(mylogf,c,plex,wsd,wmin,wf,sflag,syntype,-1);
      mod_conn_cell_to_layer_gauss(mylogf,c,plin,wsd,wmin,wf,sflag,syntype,-1);

      mod_conn_layer_to_cell_gauss(mylogf,plex,c,140.0,wmin,0.5,sflag,1,-1);
      mod_conn_layer_to_cell_gauss(mylogf,plin,c, 70.0,wmin,0.2,sflag,2,-1);

    
      // Clean up, and look for next CUSTOMIZE
      free_2d_carray(slist,ns);
      
      tp = tp->next;
      
  ***/
}
/**************************************-**************************************/
/*                                                                           */
/*                             POP_UTIL_ADD_NOISE                            */
/*                                                                           */
/*  gfg - Gaussian filtered Gaussian noise                                   */
/*  gshot - Poisson shot, shot is Gaussian shaped.                           */
/*                                                                           */
/*****************************************************************************/
void pop_util_add_noise(data,n,seed,ntype,nmu,nsd,ntsd,tscale)
     float *data;    /* Add noise to this array [n] */
     int n,seed;
     char ntype[];   /* Type of noise */
     float nmu,nsd;  /* Mean and SD of noise */
     float ntsd;     /* Temporal SD for smoothing noise */
     float tscale;   /* model tscale */
{
  int tseed,ns;
  float *sf,*noise,*snoise,rate,sampling,amp;

  /*printf("  POP_UTIL_ADD_NOISE\n");
    printf("    type  %s\n",ntype);
    printf("    mu    %f\n",nmu);
    printf("    sd    %f\n",nsd);
    printf("    tsd   %f\n",ntsd);*/

  if (strcmp(ntype,"gfg")==0){
    noise = gaussian_corr_noise(ntsd,nmu,nsd,n,seed);
    add_to_farray(data,noise,n);
    myfree(noise);
  }else if (strcmp(ntype,"gshot")==0){
    rate = nmu; /* Hz */
    amp = ntsd;
    tseed = seed;
    sampling = 1.0/tscale;
    make_poisson_float_spikes(&sf,&ns,0,n,sampling,rate,0.0,0.0,1,&tseed);
    noise = expand_spike_farray(sf,ns,0,n);
    snoise = smooth_with_gaussian(noise,n,nsd,0.01);
    multiply_farray(snoise,n,amp * sqrt(2.0*M_PI)*nsd);
    myfree(noise); myfree(sf);
    noise = snoise;
    /*append_farray_plot("zzz.noise","noise",noise,n,1);*/
    add_to_farray(data,noise,n);
    myfree(noise);
  }else{
    printf("ntype = %s\n",ntype);
    exit_error("POP_UTIL_ADD_NOISE","Unknown noise type");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                      POP_UTIL_SET_FLAG_GRID_INPUT_LAYER                   */
/*                                                                           */
/*  Set to '1' all points in 'flag' that have connections from the named     */
/*  lgn layer 'lgn_name' to the post-synaptic layer 'lay'.                   */
/*                                                                           */
/*  NOTES                                                                    */
/*    This handles regular and paired LGN inputs                             */
/*                                                                           */
/*****************************************************************************/
void pop_util_set_flag_grid_input_layer(mylogf,flag,xn,yn,onflag,lay,lgn_name)
     char *mylogf;
     int **flag,xn,yn;          // Flag and grid dimensions
     int onflag;                // 0-off connections, 1-on connections
     struct pop_layer *lay;     // Post-syn layer
     char *lgn_name;            // Name of pre-syn LGN layer
{
  int i,j,k,l;
  int n,*xi,*yi;
  struct pop_cell *c;
  struct pop_lgn_in *ln,*lp;

  for(i=0;i<lay->xn;i++){
    for(j=0;j<lay->yn;j++){
      for(k=0;k<lay->zn;k++){ // For each cell
	//printf("  cell - %d %d %d\n",i,j,k);
	c = &(lay->c[i][j][k]);

	for(l=0;l<c->lgn_n;l++){ // For each set of LGN inputs for this cell

	  ln = c->lgnin[l];

	  //
	  //  Regular inputs
	  //
	  if (strcmp(ln->lay->name,lgn_name)==0){ // Only specified LGN pop

	    if ((ln->lay->xn != xn) || (ln->lay->yn != yn)){ // seems wasteful
	      mylogx(mylogf,"POP_UTIL_SET_FLAG_GRID_INPUT_LAYER","Bad xn,yn");
	    }

	    if (onflag == 0){          // connections from OFF cells
	      n  = ln->cn0;
	      xi = ln->cx0;
	      yi = ln->cy0;
	    }else if (onflag == 1){    // connections from ON cells
	      n  = ln->cn1;
	      xi = ln->cx1;
	      yi = ln->cy1;
	    }else
	      mylogx(mylogf,"POP_UTIL_SET_FLAG_GRID_INPUT_LAYER","Bad onflag");
	    
	    for(l=0;l<n;l++){
	      flag[xi[l]][yi[l]] = 1;
	    }
	  }

	  //
	  //  Paired inputs
	  //
	  if (ln->cpair != NULL){
	    lp = ln->cpair;
	    if (strcmp(lp->lay->name,lgn_name)==0){

	      if ((lp->lay->xn != xn) || (lp->lay->yn != yn)){
		mylogx(mylogf,"POP_UTIL_SET_FLAG_GRID_INPUT_LAYER",
		       "Bad xn,yn");
	      }
	    
	      if (onflag == 0){          // connections from OFF cells
		n  = lp->cn0;
		xi = lp->cx0;
		yi = lp->cy0;
	      }else if (onflag == 1){    // connections from ON cells
		n  = lp->cn1;
		xi = lp->cx1;
		yi = lp->cy1;
	      }else
		mylogx(mylogf,"POP_UTIL_SET_FLAG_GRID_INPUT_LAYER",
		       "Bad onflag");
	      
	      for(l=0;l<n;l++){
		if (xi[l] != -1)  // Coords -1 indicates NO PAIRED CONNECTION
		  flag[xi[l]][yi[l]] = 1;
	      }
	    }
	  }
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                      POP_UTIL_EXPAND_DS01_PAIRED_SPIKES                   */
/*                                                                           */
/*  Return a float array of length 'tn' which sums the interactions from     */
/*  all paired spikes.  The ith input connection is paired with the          */
/*  cn/2 + ith connection.                                                   */
/*                                                                           */
/*  WYETH - OLD WAY OF PAIRING lgn inputs?                                   */
/*                                                                           */
/*****************************************************************************/
float *pop_util_expand_ds01_paired_spikes(cx,cy,cn,xn,yn,s,cnt,tn,cw)
     int *cx,*cy,cn;           /* Coordinates of input connections [cn] */
     int xn,yn;                /* size of input grid [xn][yn] */
     float ***s;               /* Spike trains for grid points [xn][yn] */
     int **cnt;                /* number of spikes [xn][yn] */
     int tn;                   /* gt[tn] */
     float *cw;                /* spike weight for connection [cn], or NULL */
{
  int i,j,k;
  int xi,yi,t,tcnt;
  float *ts,*x;

  int nmask,npair;
  float *mask,*ma,mc,mu,tt,mcent;

  /*** WYETH - cw is ignored ***/

  /*** WYETH - masks may ultimately be customized for DX's ***/
  /*** WYETH - SHOULD BE DONE BY LOOKUP TABLES, NO MATH HERE ***/
  nmask = 40; /* Assuming this is msec */
  mask = (float *)myalloc(nmask*sizeof(float));
  mc = 2.0*M_PI*1.5/(float)nmask;
  mcent = (float)nmask/2.0;

  /*** SINEWAVE, peaks in middle of mask. ***/
  for(i=0;i<nmask;i++){
    tt = (float)i;
    mask[i] = -2.0*sin(mc*tt);
  }
  
  /*** DOG, could use this instead of sinusoid ***/
  /*
  for(i=0;i<nmask;i++){
    tt = (float)i;
    mask[i] = (float)(func_gaussian(tt,mcent,nmask/8.0) - 
		      func_gaussian(tt,mcent,nmask/6.0));
  }
  make_max_const_farray(mask,nmask,2.0);
  printf("SUM OF MASK = %f\n",sum_farray(mask,nmask,0,nmask));*/

  /*append_farray_plot("zz.DS01_MASK.pl","mask",mask,nmask,1);*/

  x = get_zero_farray(tn);
  ma = (float *)myalloc(tn*sizeof(float));

  npair = cn/2;
  for(i=0;i<npair;i++){  /* For each input pair */
    xi = cx[i];
    yi = cy[i];
    if ((xi<0)||(xi>=xn)||(yi<0)||(yi>=yn))
      exit_error("POP_UTIL_EXPAND_GRID_SPIKES","invalid grid coordinate");
    
    /*** Set up 'ma' mask array ***/
    for(j=0;j<tn;j++)
      ma[j] = 0.0; /* Re-zero */

    ts = s[xi][yi];
    tcnt = cnt[xi][yi];

    for(j=0;j<tcnt;j++){
      t = my_rint(ts[j]);
      if ((t >= 0)&&(t < (tn-nmask))){
	for(k=0;k<nmask;k++){
	  ma[t+k] += mask[k];
	}
      }
    }

    /*append_farray_plot("zz.DS01_MASK.pl","ma",ma,tn,1);*/

    xi = cx[npair+i];
    yi = cy[npair+i];
    if ((xi<0)||(xi>=xn)||(yi<0)||(yi>=yn))
      exit_error("POP_UTIL_EXPAND_GRID_SPIKES","invalid grid coordinate");

    ts = s[xi][yi];
    tcnt = cnt[xi][yi];
    for(j=0;j<tcnt;j++){
      t = my_rint(ts[j]);
      if ((t >= 0)&&(t < tn)){
	if (ma[t] > 0.0)  /* Only used positive weights */
	  x[t] += ma[t];  /* WYETH - would apply 'w' here, if used */
      }
    }
    /*append_farray_plot("zz.DS01_MASK.pl","xt",x,tn,1);*/
    /*exit(0);*/
  }

  myfree(ma);
  
  return x;
}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_UTIL_EXPAND_SPIKES                         */
/*                                                                           */
/*  Return a float array of length 'tn' which sums all spikes from 's'.      */
/*                                                                           */
/*****************************************************************************/
float *pop_util_expand_spikes(s,n,tn,w,sdf,sdtau)
     float *s;      // Spike train [n]
     int n;         // spike count
     int tn;        // length of expanded array
     float w;       // spike weight
     float sdf;     // Synaptic depression, frac reduction after spike
     float sdtau;   // Synaptic depression, recovery time const
{
  float *x;

  x = get_zero_farray(tn);   // The expanded spike arrays, to be returned

  spikeu_expand_spikes_wsd_accum(s,n,x,tn,w,sdf,sdtau);

  return x;
}
/**************************************-**************************************/
/*                                                                           */
/*                         POP_UTIL_EXPAND_GRID_SPIKES                       */
/*                                                                           */
/*  Return a float array of length 'tn' which sums all spikes from the       */
/*  grid coordinates (cx,cy).                                                */
/*                                                                           */
/*****************************************************************************/
float *pop_util_expand_grid_spikes(cx,cy,cn,xn,yn,s,cnt,tn,cw,sdf,sdtau)
     int *cx,*cy,cn;     // Coordinates of input connections [cn]
     int xn,yn;          // size of input grid [xn][yn]
     float ***s;         // Spike trains for grid points [xn][yn]
     int **cnt;          // number of spikes [xn][yn]
     int tn;             // gt[tn]
     float *cw;          // spike weight for connection [cn], or NULL
     float sdf;          // Synaptic depression, frac reduction after spike
     float sdtau;        // Synaptic depression, recovery time const
{
  int i,j;
  int xi,yi,wflag;
  float *tx,*x,w;

  x = get_zero_farray(tn);

  if (cw == (float *)NULL)
    wflag = 0;
  else
    wflag = 1;
  
  for(i=0;i<cn;i++){  // For each input contributing spikes
    xi = cx[i];
    yi = cy[i];
    if (wflag)
      w = cw[i];
    else
      w = 1.0;

    if ((xi<0)||(xi>=xn)||(yi<0)||(yi>=yn))
      exit_error("POP_UTIL_EXPAND_GRID_SPIKES","invalid grid coordinate");

    tx = pop_util_expand_spikes(s[xi][yi],cnt[xi][yi],tn,w,sdf,sdtau);

    for(j=0;j<tn;j++)
      x[j] += tx[j];

    myfree(tx);
  }

  return x;
}
/**************************************-**************************************/
/*                                                                           */
/*                          POP_UTIL_SUM_GT_SPK_GRID                         */
/*                                                                           */
/*  1. Sum the spike trains from the specified grid points, with optional    */
/*     weighting.                                                            */
/*  2. Convolve with the named DOE.                                          */
/*  3. Add result into 'gt'                                                  */
/*                                                                           */
/*****************************************************************************/
void pop_util_sum_gt_spk_grid(lay,cx,cy,cn,cw,gw,xn,yn,s,cnt,m,name,mechr,
				     gt,gnt,
				     tn,dumpname,subclass)
     struct pop_layer *lay;    // layer (for NMDA params)
     int *cx,*cy,cn;           // Coordinates of input connections [cn]
     float *cw;                // Weight for connection [cn], maybe NULL
     float gw;                 // Single, global weight
     int xn,yn;                // size of input grid [xn][yn]
     float ***s;               // Spike trains for grid points [xn][yn]
     int **cnt;                // number of spikes [xn][yn]
     struct model_struct *m;   // Model parameter pair list.
     char name[];              // Name for PSG parameters (mppl way)
     struct pop_mech *mechr;   // Name for PSG parameters (onode way)
     float *gt;                // conductance to add to
     float *gnt;               // conductance to add NMDA to
     int tn;                   // gt[tn]
     char dumpname[];
     char *subclass;           // NULL or ds01 - how to process spikes
{
  int nm,nn;
  float *mask,*xsp,*g,*nmda_mask,sdf,sdtau;
  char tstr[SLEN];
  struct pop_mech *mech;
  struct nmda_param *tp;

  if (m->ppl == NULL){
    sdf = mechr->sdf;
    sdtau = mechr->sdtau;
    mask = pop_util_get_mask_diff_exp_mech(mechr,&nm);
  }else{
    sprintf(tstr,"%s_sdf",name);
    sdf = paramfile_get_float_param_or_exit(m->ppl,tstr);
    sprintf(tstr,"%s_sdtau",name);
    sdtau = paramfile_get_float_param_or_exit(m->ppl,tstr);
    mask = pop_util_get_mask_diff_exp(m,name,&nm);
  }

  //printf("  name = %s\n",name);
  //printf("  lay name = %s\n",lay->name);

  // If time consts are in msec, so is returned mask
  if (dumpname[0] != '\0')
    append_farray_plot(dumpname,"LGN_EX_PSC",mask,nm,1);

  // Get expanded, weighted sum of spikes, w/ synaptic depression
  if (subclass == NULL)
    xsp = pop_util_expand_grid_spikes(cx,cy,cn,xn,yn,s,cnt,tn,cw,sdf,sdtau);
  else if (strcmp(subclass,"ds01")==0)
    xsp = pop_util_expand_ds01_paired_spikes(cx,cy,cn,xn,yn,s,cnt,tn,cw);
  else
    exit_error("POP_UTIL_SUM_GT_SPK_GRID","Unknown subclass");

  if ((gw >= 0.0) && (gw != 1.0)){ // WYETH LGNW - adding a weight to LGN inputs
    // WYETH - I added the condition to block negative weights because it
    //  seems that w = -1.0 is set some places, possibly as a flag?
    multiply_farray(xsp,tn,gw); // Multiply all input by a global weight
  }

  if (dumpname[0] != '\0')
    append_farray_plot(dumpname,"sum_spikes",xsp,tn,1);

  // Convolve spikes w/ mask
  g = convolve_with_mask_causal(xsp,tn,mask,nm);
  if (dumpname[0] != '\0')
    append_farray_plot(dumpname,"g_conv_AMPA",g,tn,1);

  // Add result to 'gt' which is passed in
  add_to_farray(gt,g,tn);
  myfree(g);
  myfree(mask);

  // NMDA
  tp = lay->nmda_p;
  if (gnt != (float *)NULL){
    nn = (int)(1.0/tp->nmda_alpha_3 * 6.0);
    if (nn > tn){
      //printf("  *** Shortening NMDA mask from %d to ",nn);
      nn = tn-1;
      //printf("length of data minus 1, %d\n",nn);
    }
    
    //nmda_mask = pop_util_get_nmda_shape_layer(lay,nn);
    nmda_mask = pop_util_get_nmda_shape(tp->nmda_alpha_1,tp->nmda_alpha_2,
					tp->nmda_alpha_3,tp->lgn_nmda_amp_1,
					tp->lgn_nmda_amp_2,
					tp->lgn_nmda_amp_3,nn);


    if (dumpname[0] != '\0')
      append_farray_plot(dumpname,"NMDA_PSC",nmda_mask,nn,1);
    g = convolve_with_mask_causal(xsp,tn,nmda_mask,nn);
    if (dumpname[0] != '\0')
      append_farray_plot(dumpname,"g_conv_NMDA",g,tn,1);
    add_to_farray(gnt,g,tn);

    myfree(g);
    myfree(nmda_mask);
  }

  myfree(xsp);
}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_UTIL_PAIR_PROCESS                          */
/*                                                                           */
/*****************************************************************************/
void pop_util_pair_process(logf,mechp,x1,x2,n,mask,nm,mask_nmda,nn,gt,gnt,
			   dumpname)
     char *logf;              // logfile
     struct pop_mech *mechp;  // Pairing mechanism
     float *x1,*x2;           // expanded spike trains [n]
     int n;                   // length of expanded data
     float *mask;             // AMPA mask [nm]
     int nm;                  //
     float *mask_nmda;        // NMDA mask [nn]
     int nn;                  //
     float *gt;               // Add AMPA result to this array
     float *gnt;              // Add NMDA result to this array (NULL if none)
     char *dumpname;          // Dumpfile name, NULL for none
{
  float *g1,*g2,*g,*gn;

  // Convolve spikes w/ mask
  g1 = convolve_with_mask_causal(x1,n,mask,nm);
  g2 = convolve_with_mask_causal(x2,n,mask,nm);

  if (strcmp(mechp->type,"multiply")==0){

    g = multiply_farrays(g1,g2,n);

    if (mechp->normv != 1.0)
      multiply_farray(g,n,mechp->normv);  // Renormalize

  }else if (strcmp(mechp->type,"divide")==0){
    //  g1/(a + b*g2),  typically, a=1, b=1

    multiply_farray(g2,n,mechp->b);
    add_const_farray(g2,n,mechp->a);

    g = divide_farrays(g1,g2,n);

    // No post-divide normalization
    
  }else
    mylogx(logf,"POP_UTIL_PAIR_PROCESS","Unknown mech type");


  add_to_farray(gt,g,n); // Add result to 'gt' which is passed in

  if (dumpname[0] != '\0'){
    append_farray_plot(dumpname,"g1_conv_AMPA",g1,n,1);
    append_farray_plot(dumpname,"g2_conv_AMPA",g2,n,1);
    printf("NORMw = %f\n",mechp->normv);
    append_farray_plot(dumpname,"gMULT",g,n,1);
    append_farray_plot(dumpname,"gt_TOT",gt,n,1);
  }


  if (gnt != (float *)NULL){ // NMDA

    gn = convolve_with_mask_causal(g,n,mask_nmda,nn);
    //if (dumpname[0] != '\0')
    //append_farray_plot(dumpname,"g_conv_NMDA",g,n,1);
    add_to_farray(gnt,gn,n);

    myfree(gn);
  }

  myfree(g1);
  myfree(g2);
  myfree(g);
}
/**************************************-**************************************/
/*                                                                           */
/*                         POP_UTIL_GT_SUM_GRID_PAIR                         */
/*                                                                           */
/*  1. Sum pair-wise interactions of spikes from the specified grid points.  */
/*  2. Add result into 'gt'                                                  */
/*                                                                           */
/*****************************************************************************/
void pop_util_gt_sum_grid_pair(logf,lay,lgn_in1,lgn_in2,cw,xn,yn,
			       gt,gnt,tn,dumpname)
     char *logf;
     struct pop_layer *lay;       // Post-syn layer, used for NMDA params
     struct pop_lgn_in *lgn_in1;  // LGN inputs, original
     struct pop_lgn_in *lgn_in2;  // LGN inputs, paired
     float *cw;                   // Weight for connection [cn], maybe NULL
     int xn,yn;                   // size of input grid [xn][yn]
     float *gt;                   // Add AMPA to this conductance [tn]
     float *gnt;                  // Add NMDA to this conductance [tn]
     int tn;                      // length of gt and gtn
     char dumpname[];
{
  int i;
  int nm,nn,cn;
  int *cx_1,*cy_1,*cx_2,*cy_2,xi,yi,**n_1,**n_2,n1,n2;
  float ***s_1,***s_2,*s1,*s2;
  float *mask,*mask_nmda,sdf,sdtau,*x1,*x2,w;
  char tstr[SLEN];
  struct pop_mech *mechr,*mechp;
  struct nmda_param *tp;

  mechr = lgn_in1->mech;  // Post-syn receptor mechanism
  mechp = lgn_in2->mech;  // Pairing mechanism

  //
  //  Post-syn masks
  //
  sdf   = mechr->sdf;
  sdtau = mechr->sdtau;
  mask  = pop_util_get_mask_diff_exp_mech(mechr,&nm);

  if (dumpname[0] != '\0') // If time consts are in msec, so is returned mask
    append_farray_plot(dumpname,"LGN_EX_PSC",mask,nm,1);

  if (gnt != (float *)NULL){
    tp = lay->nmda_p;

    nn = (int)(1.0/tp->nmda_alpha_3 * 6.0);
    if (nn > tn){
      mylog(logf,"*** Shortening NMDA mask");
      nn = tn-1;
    }

    mask_nmda = pop_util_get_nmda_shape(tp->nmda_alpha_1,tp->nmda_alpha_2,
					tp->nmda_alpha_3,tp->lgn_nmda_amp_1,
					tp->lgn_nmda_amp_2,
					tp->lgn_nmda_amp_3,nn);
    if (dumpname[0] != '\0')
      append_farray_plot(dumpname,"NMDA_PSC",mask_nmda,nn,1);
  }else{
    mask_nmda = NULL;
    nn = -1;
  }


  //
  //  Pre-compute any normalization constants for the pair interaction
  //
  if (strcmp(mechp->type,"multiply")==0){

    mechp->normv = ((double)sum_farray(mask,nm,0,nm) /
		    (double)sum_square_farray(mask,nm));
    //printf("WYETH HERE  NORM v  = %f\n",mechp->normv);
    //exit(0);
  }


  //
  //  ON cells
  //
  cx_1 = lgn_in1->cx1;
  cy_1 = lgn_in1->cy1;
  cx_2 = lgn_in2->cx1;
  cy_2 = lgn_in2->cy1;
  s_1  = lgn_in1->lay->s[1];       // spikes for LGN regular
  n_1  = lgn_in1->lay->cnt[1];     // counts
  s_2  = lgn_in2->lay->s[1];       // spikes for LGN regular
  n_2  = lgn_in2->lay->cnt[1];     // counts

  cn = lgn_in1->cn1;  // Number of connections
  if (lgn_in2->cn1 != cn)
    mylogx(logf,"POP_UTIL_GT_SUM_GRID_PAIR","Paired connection count error");

  for(i=0;i<cn;i++){  // For each pair

    if (cw != NULL)   // Connection weight (is this used?)
      w = cw[i];
    else
      w = 1.0;

    xi = cx_1[i];  // Temporary, *** will be REUSED below
    yi = cy_1[i];
    s1 = s_1[xi][yi];       // spikes for LGN regular
    n1 = n_1[xi][yi];

    xi = cx_2[i];  // -1 if no match
    if (xi >= 0){
      yi = cy_2[i];
      s2 = s_2[xi][yi];    // spikes for LGN pair
      n2 = n_2[xi][yi];    // counts

      //  Now we have the two spike trains:  s1[n1], and s2[n2]

      // Get expanded spike trains, weighted, w/  synaptic depression
      x1 = pop_util_expand_spikes(s1,n1,tn,w,sdf,sdtau);
      x2 = pop_util_expand_spikes(s2,n2,tn,w,sdf,sdtau);
      if (dumpname[0] != '\0'){
	append_farray_plot(dumpname,"expanded_train_1",x1,tn,1);
	append_farray_plot(dumpname,"expanded_train_2",x2,tn,1);
      }
      
      // Compute pair interaction, add to 'gt' and 'gnt'
      pop_util_pair_process(logf,mechp,x1,x2,tn,mask,nm,mask_nmda,nn,gt,gnt,
			    dumpname);
      
      myfree(x1);
      myfree(x2);
    }
  }


  //
  //  OFF cells
  //
  cx_1 = lgn_in1->cx0;
  cy_1 = lgn_in1->cy0;
  cx_2 = lgn_in2->cx0;
  cy_2 = lgn_in2->cy0;
  s_1  = lgn_in1->lay->s[0];       // spikes for LGN regular
  n_1  = lgn_in1->lay->cnt[0];     // counts
  s_2  = lgn_in2->lay->s[0];       // spikes for LGN regular
  n_2  = lgn_in2->lay->cnt[0];     // counts

  cn = lgn_in1->cn0;  // Number of connections
  if (lgn_in2->cn0 != cn)
    mylogx(logf,"POP_UTIL_GT_SUM_GRID_PAIR","Paired connection count error");

  for(i=0;i<cn;i++){  // For each pair

    if (cw != NULL)   // Connection weight (is this used?)
      w = cw[i];
    else
      w = 1.0;

    xi = cx_1[i];
    yi = cy_1[i];
    s1 = s_1[xi][yi];       // spikes for LGN regular
    n1 = n_1[xi][yi];

    xi = cx_2[i];  // -1 if no match
    if (xi >= 0){
      yi = cy_2[i];
      s2 = s_2[xi][yi];    // spikes for LGN pair
      n2 = n_2[xi][yi];    // counts

      // Get expanded spike trains, weighted, w/  synaptic depression
      x1 = pop_util_expand_spikes(s1,n1,tn,w,sdf,sdtau);
      x2 = pop_util_expand_spikes(s2,n2,tn,w,sdf,sdtau);
      if (dumpname[0] != '\0'){
	append_farray_plot(dumpname,"expanded_train_1",x1,tn,1);
	append_farray_plot(dumpname,"expanded_train_2",x2,tn,1);
      }

      // Compute pair interaction, add to 'gt' and 'gnt'
      pop_util_pair_process(logf,mechp,x1,x2,tn,mask,nm,mask_nmda,nn,gt,gnt,
			    dumpname);
      
      myfree(x1);
      myfree(x2);

    }
  }

  myfree(mask);

  if (mask_nmda != NULL)
    myfree(mask_nmda);

}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_UTIL_SET_GAD_LAYER                         */
/*                                                                           */
/*  Make or zero storage for gad - adaptation conductance.                   */
/*                                                                           */
/*****************************************************************************/
void pop_util_set_gad_layer(lay,xn,yn,tnsec)
     struct pop_layer *lay;    // Layer
     int xn,yn;                // Grid dimensions
     float tnsec;              // Length of conductance
{
  int i,j,k;
  int tn;
  struct pop_cell *c;

  // WYETH - I had to set all three (four?) of these, why are they linked
  // WYETH - I had to set all three (four?) of these, why are they linked
  // WYETH - I had to set all three (four?) of these, why are they linked

  tn = my_rint(tnsec * 1000.0); // WYETH - assuming sampling is MSEC

  for(i=0;i<lay->xn;i++){
    for(j=0;j<lay->yn;j++){
      for(k=0;k<lay->zn;k++){ // For each cell
	c = &(lay->c[i][j][k]);

	if (c->gta0 != NULL){  // WYETH - How do we know it is 'tn' long????
	  //zero_farray(c->gta0,tn);  // WYETH - Trouble w/ length here???
	  zero_farray(c->gta0,c->gta_max); // NEW
	}else{
	  ; //printf("IT IS NULL for %s\n",lay->name);
	}

	if (c->gtx0 == NULL){
	  c->gtx0 = get_zero_farray(tn);    // Initialize x, i, and a,
	  if (lay->lgn_gain == 1)
	    c->ggain = get_zero_farray(tn); // Initialize x, i, and a,
	  c->gti0 = get_zero_farray(tn);    //   assuming that all 3 linked
	  if (lay->nmda_p != NULL)
	    c->gtn0 = get_zero_farray(tn);
	  c->n0 = tn;
	  // c->samp0 should already be set
	}else{
	  if (tn != c->n0)
	    exit_error("POP_UTIL_SET_GAD_LAYER","tn changed");
	  zero_farray(c->gtx0,tn);
	  if (lay->lgn_gain == 1)
	    zero_farray(c->ggain,tn);
	  zero_farray(c->gti0,tn);
	  if (lay->nmda_p != NULL)
	    zero_farray(c->gtn0,tn);
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                POPU_SUM_GAIN                              */
/*                                                                           */
/*  Add to gain signal 'gsum'.                                               */
/*                                                                           */
/*****************************************************************************/
void popu_sum_gain(logf,ln,gsum,tn,c,tscale)
     char *logf;              // Log file
     struct pop_lgn_in *ln;   // LGN input
     float *gsum;             // Add to this sum [tn]
     int tn;                  // length of time signal
     struct pop_cell *c;      // Post-syn cell
     float tscale;            // sec / time samp in 'graw' below
{
  int i,j,k;
  int xn,yn,n,xi,yi,*x,*y,ndup;
  float ***graw,*tg;
  char gtype[SLEN];
  //
  // For 'sz1_sum' below
  //
  int x1,x2,y1,y2,xc,yc,ri,ti;
  float r;    // Radius for summation field
  float sig;  // SD for weighted sum across field
  float dx2,dy,dist,w,wtot;
  float *tsum;


  strcpy(gtype,"sz1_sum");

  // Get values from LGN layer
  graw = ln->lay->rgl;
  xn   = ln->lay->xn;
  yn   = ln->lay->yn;

  ndup = my_rint(c->samp0 * tscale);
  if ((ndup < 1) || (ndup > 10))
    mylogx(logf,"POPU_SUM_GAIN","ndup out of range");
    
  if (graw == NULL)
    mylogx(logf,"POPU_SUM_GAIN","LGN 3D gain response is NULL");

  if (gsum == NULL)
    mylogx(logf,"POPU_SUM_GAIN","gsum is NULL");

  //printf("HERE inside 1\n");

  if (strcmp(gtype,"input_sum")==0){

    n = ln->cn0;  // OFF inputs
    x = ln->cx0;
    y = ln->cy0;
    for(i=0;i<n;i++){  // For each input
      xi = x[i];
      yi = y[i];
      tg = graw[xi][yi];
      for(j=0;j<tn;j++){
	k = j/2;
	gsum[j] += tg[k];
      }
    }
    
    n = ln->cn1;  // ON inputs
    x = ln->cx1;
    y = ln->cy1;
    for(i=0;i<n;i++){  // For each input
      xi = x[i];
      yi = y[i];
      tg = graw[xi][yi];
      for(j=0;j<tn;j++){
	k = j/2;
	gsum[j] += tg[k];
      }
    }
  }else if (strcmp(gtype,"sz1_sum")==0){

    //
    //  Average over a Gaussian with SD set by 'sz1', related to Gabor SD
    //

    xc = my_rint(c->rfx);
    yc = my_rint(c->rfy);

    // WYETH ATTRIB
    sig = pop_cell_attrib_get_f(c,"sz1");
    //sig = c->sz1;

    r  = 3.0 * sig;
    ri = my_rint(r);

    x1 = xc - ri;
    x2 = xc + ri;
    if (x1 < 0)
      x1 = 0;
    if (x2 >= xn)
      x2 = xn-1;

    y1 = yc - ri;
    y2 = yc + ri;
    if (y1 < 0)
      y1 = 0;
    if (y2 >= yn)
      y2 = yn-1;

    tsum = get_zero_farray(tn);
    wtot = 0.0;
    for(i=x1;i<=x2;i++){
      dx2 = (float)((i-xc)*(i-xc));
      for(j=y1;j<=y2;j++){
	dy = (float)(j-yc);
	tg = graw[i][j];
	dist = sqrt(dx2 + dy*dy);
	w = (float)func_gaussian_one(dist,0.0,sig);
	wtot += w;
	
	for(ti=0;ti<tn;ti++){
	  k = ti/2;
	  tsum[ti] += w * tg[k];
	}
      }
    }

    for(ti=0;ti<tn;ti++){
      tsum[ti] /= wtot;      // Convert Sum to Averge
      gsum[ti] += tsum[ti];  // Add average to 'gsum'
    }

  }else if (strcmp(gtype,"sz1_sum")==0){
    mylogx(logf,"POPU_SUM_GAIN","Unknown gain type.");
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                      POP_UTIL_SET_GT_ON_OFF_INPUT_LAYER                   */
/*                                                                           */
/*  For the layer 'lay', set the pre-computed excitatory conductance that    */
/*  results from the LGN ON and OFF inputs.                                  */
/*                                                                           */
/*  *** NOTE:  use '...set_gad_layer' above for layer w/o LGN input.         */
/*                                                                           */
/*****************************************************************************/
void pop_util_set_gt_on_off_input_layer(logf,lay,m,
					xn,yn,tnsec,seed,outfile,binoc_flag)
     char *logf;               // Log file
     struct pop_layer *lay;    // Post-syn layer
     struct model_struct *m;   // Model parameter pair list.
     int xn,yn;                // Grid dimensions
     float tnsec;              // Length of conductance (seconds)
     int seed;                 //
     char outfile[];           //
     int binoc_flag;           // 0-monocular left eye, 1-binocular
{
  int i,j,k,gi,ti;
  int nseed,*seedlist,si,tn;
  float *gt,*gnt,*gsum,*cw,tscale,gfrac,gadd,gmult,*sm_gain,tg,ocdom,gw;
  char dumpf[SLEN];
  struct pop_cell *c;
  struct pop_lgn_in *ln;
  struct pop_mech *mechr;
  int nm;
  float sdf,sdtau,*xsp,*mask,*g;

  //
  //  WYETH - HACK - ASSUMPTION for MESH MODEL:
  //  The post-syn layer's mech for "lgn_ex" is to be used.
  //
  if (m->ppl == NULL){
    mechr = pop_cell_get_mech_by_name(lay,"lgn_ex");
    if (mechr == NULL)
      exit_error("POP_UTIL_SET_GT_ON_OFF_INPUT_LAYER","Assumed mech missing");

    sdf = mechr->sdf;
    sdtau = mechr->sdtau;
    mask = pop_util_get_mask_diff_exp_mech(mechr,&nm);
  }else{
    exit_error("POP_UTIL_SET_GT_ON_OFF_INPUT_LAYER","Old way");
  }

  cw = (float *)NULL;  // Currently no variation in weights for connections

  // Pick randomization seed for each cell
  nseed = lay->xn * lay->yn * lay->zn;
  seedlist = get_seeds(seed,100000,nseed);

  tn = my_rint(tnsec * lay->c[0][0][0].samp0);
  tscale = onode_getpar_flt_exit(m->o,"tscale");

  if (lay->lgn_gain == 1){
    //printf(" tau = %f   tscale %f\n",lay->lgn_gain_tau,tscale);
    gfrac = exp(- tscale/lay->lgn_gain_tau);
    gadd  = lay->lgn_gain_ca;
    gmult = lay->lgn_gain_cb;
    //printf("GFRAC = %f   add,mult %f %f \n",gfrac,gadd,gmult);
  }

  si = 0; // Seedlist index
  for(i=0;i<lay->xn;i++){
    for(j=0;j<lay->yn;j++){
      for(k=0;k<lay->zn;k++){ // For each cell in post-syn layer
	c = &(lay->c[i][j][k]);

	//  For this cell, zero:  gtx0, gti0, gtn0, ggain

	popc_precomp_reset(c,lay,tn);

	gt = c->gtx0;     // AMPA conductance total
	gnt = c->gtn0;    // NMDA conductance total
	gsum = c->ggain;  // gain total
	if ((i==lay->xn/2) && (j==lay->yn/2) && (outfile[0] != '\0')){
	  strcpy(dumpf,outfile);  // Dump file name
	  printf("  Dumping %s for cell %d %d\n",dumpf,i,j);
	}else
	  dumpf[0] = '\0';  // Dump file name

	if (binoc_flag == 0){

	  if ((c->lgn_n > 0) && (c->lgnin == NULL)){
	    //
	    //  IRREGULAR MESH
	    //

	    //  Get expanded spikes, convolve w/ "lgn_ex" mask, add to total
	    xsp = popc_precomp_get_expanded_spikes(c,tn,sdf,sdtau);
	    g = convolve_with_mask_causal(xsp,tn,mask,nm); // Conv w/ mask
	    myfree(xsp); // WYETH - added 2012 May 1 to reduce mem leak
	    add_to_farray(gt,g,tn);  // Add into c->gtx0
	    myfree(g);

	    // WYETH - this shouldn't be too hard to add, if all it needs is
	    //  a pointer to the layer so it can find the gain data.
	    if (lay->lgn_gain == 1){
	      exit_error("POP_UTIL_SET_GT_ON_OFF_INPUT_LAYER",
			 "Gain not implemented for retina mesh");
	    }
	    // WYETH - This should be easy to add, just look at the
	    // block of lines for 'gt', and see '...sum_gt_spk_grid' above.
	    if (gnt != NULL){  // WYETH - THIS IS E
	      exit_error("POP_UTIL_SET_GT_ON_OFF_INPUT_LAYER",
			 "NMDA not implemented for retina mesh");
	    }

    //exit_error("WYETH HERE","Devel....");

	  }else{

	    //
	    // WYETH - the lines below could be re-done using the strategy
	    //  above, where the mask does not get rebuilt for every cell.
	    //

	    for(gi=0;gi<c->lgn_n;gi++){ // For each set of LGN inputs
	      ln = c->lgnin[gi];
	      gw = ln->gw;  // WYETH LGNW

	      // Make cw[cn

	      if (ln->cpair == NULL){
		// OFF cells
		pop_util_sum_gt_spk_grid(lay,ln->cx0,ln->cy0,ln->cn0,cw,gw,xn,
					 yn,
					 ln->lay->s[0],ln->lay->cnt[0],m,NULL,
					 ln->mech,gt,gnt,tn,dumpf,c->subclass);
		// ON cells
		pop_util_sum_gt_spk_grid(lay,ln->cx1,ln->cy1,ln->cn1,cw,gw,xn,
					 yn,
					 ln->lay->s[1],ln->lay->cnt[1],m,NULL,
					 ln->mech,gt,gnt,tn,dumpf,c->subclass);
	      }else{
		pop_util_gt_sum_grid_pair(logf,lay,ln,ln->cpair,cw,xn,yn,
					  gt,gnt,tn,dumpf);
	      }

	      if (lay->lgn_gain == 1){
		// WYETH - this should probably only run once, thus, only
		//         designed for one LGN input
		popu_sum_gain(logf,ln,gsum,tn,c,tscale);
	      }
	    }
	  }

	  //
	  //   Smooth and apply the gain control signal, if flag is set
	  //
	  if (lay->lgn_gain == 1){
	    //float gfrac; //  = 0.9;

	    sm_gain = get_farray(tn); // Storage will replace 'gsum' below

	    sm_gain[0] = gadd + gmult*gsum[0];
	    gt[0] /= sm_gain[0];   // Divide g(t) by gain expression
	    if (lay->nmda_p != NULL)
	      gnt[0] /= sm_gain[0];  // same for g_NMDA
	    for(ti=1;ti<tn;ti++){
	      tg = gadd + gmult*gsum[ti];
	      sm_gain[ti] = gfrac * sm_gain[ti-1] + (1.0-gfrac)*tg;
	      gt[ti] /= sm_gain[ti];   // Divide g(t) by gain expression
	      if (lay->nmda_p != NULL)
		gnt[ti] /= sm_gain[ti];  // same for g_NMDA
	    }

	    myfree(gsum);  // Points to c->ggain
	    c->ggain = sm_gain;

	  }

	}else{  // Binocular

	  ocdom = pop_cell_attrib_get_f(c,"ocdom");

	  for(gi=0;gi<c->lgn_n;gi++){ // For each set of LGN inputs
	    ln = c->lgnin[gi];
	    gw = ln->gw;  // WYETH LGNW

	    //
	    //  FOR LARS???
	    //  We are considering the following:
	    //    GEt the SF attrib for cell 'c'
	    //    compute gw (weight) as a function of 'sf' to compensate for
	    //    variations in LGN tuning vs. SF.
	    //

	    if (ln->cpair != NULL){
	      exit_error("POP_UTIL_SET_GT_ON_OFF_INPUT_LAYER",
			 "LGN pairing not coded for binoc yet.");
	    }
	    
	    if (ocdom == -1.0){ // Left eye
	      pop_util_sum_gt_spk_grid(lay,ln->cx1,ln->cy1,ln->cn1,cw,gw,xn,yn,
				       ln->lay->s[1],ln->lay->cnt[1],m,NULL,
				       ln->mech,gt,gnt,tn,dumpf,
				       c->subclass); // ON
	      pop_util_sum_gt_spk_grid(lay,ln->cx0,ln->cy0,ln->cn0,cw,gw,xn,yn,
				       ln->lay->s[0],ln->lay->cnt[0],m,NULL,
				       ln->mech,gt,gnt,tn,dumpf,
				       c->subclass); // OFF
	    }else if (ocdom == 1.0){ // Right eye
	      pop_util_sum_gt_spk_grid(lay,ln->cx1,ln->cy1,ln->cn1,cw,gw,xn,yn,
				       ln->lay->s[3],ln->lay->cnt[3],m,NULL,
				       ln->mech,gt,gnt,tn,dumpf,
				       c->subclass); // ON
	      pop_util_sum_gt_spk_grid(lay,ln->cx0,ln->cy0,ln->cn0,cw,gw,xn,yn,
				       ln->lay->s[2],ln->lay->cnt[2],m,NULL,
				       ln->mech,gt,gnt,tn,dumpf,
				       c->subclass); // OFF
	    }else{
	      printf("  Cell ocdom = %f\n",ocdom);
	      exit_error("POP_UTIL_SET_GT_ON_OFF_INPUT_LAYER",
			 "Ambiguous ocdom");
	    }
	  }

	  if (lay->lgn_gain == 1){
	    exit_error("POP_UTIL_SET_GT_ON_OFF_INPUT_LAYER",
		       "LGN gain not imp'd for binoc yet.");
	  }
	}

	// Use scale and bias from IFC ex params
	multiply_farray(gt,tn,c->cifcp->gx_scale);
	add_const_farray(gt,tn,c->cifcp->gx_bias);

	pop_util_add_noise(gt,tn,seedlist[si],
			   c->cifcp->gx_noise,
			   c->cifcp->gx_noise_mu,
			   c->cifcp->gx_noise_sd,
			   c->cifcp->gx_noise_tsd,tscale);
	si += 1; // Update seedlist counter

	if (c->cifcp->grect)
	  half_wave_rectify_farray(gt,tn);

	if (gnt != NULL){
	  multiply_farray(gnt,tn,c->cifcp->gx_scale);
	  // Don't add bias or noise, because they will be added,
	  //   unscaled, from AMPA input.
	}
      }
    }
  }

  myfree(mask);
  myfree(seedlist);
}
/**************************************-**************************************/
/*                                                                           */
/*                         POP_UTIL_CORR_NEARBY_SPIKES                       */
/*                                                                           */
/*****************************************************************************/
void pop_util_corr_nearby_spikes(xn,yn,s,cnt,periodms,corr_d,corr_p,tsd,seed)
     int xn,yn;       /* Size of grid of spikes */
     float ***s;      /* s[xn][yn][cnt[][]] spike times (ms) */
     int **cnt;       /* cnt[xn][yn] */
     float periodms;  /* max spike time (ms) */
     int corr_d;      /* Distance in cells:  2,3,... */
     float corr_p;    /* Probability of including neighbor spike */
     float tsd;       /* SD of Gaussian for time offset */
     int seed;        /* Randomization seed */
{
  int i,j,k,ii,jj;
  int nn,**skn,tn,d0x,d0y,maxn,n;
  float pkeep,***sk,*ts,*ss,tmax;

  if (seed > 0)
    seed = -seed;

  tmax = periodms - 1.0; /* Cut out spikes jittered beyond this time */

  nn = corr_d * corr_d - 1;                /* Number of neighbor cells */
  pkeep = 1.0/(corr_p * (float)nn + 1.0);  /* Fraction of spikes to keep */

  /*
    printf("    correl distance = %d cells, ppick = %.2f\n",corr_d,corr_p);
    printf("    Keeping %.1f%% of spikes\n",100.0*pkeep);
    printf(" CELL 30-32 has %d spikes\n",cnt[30][32]);*/
  

  /*** 1. Pick spikes to keep ***/
  maxn = max_of_2d_iarray(cnt,xn,yn);
  ts = (float *)myalloc(5*maxn*sizeof(float)); /* Temp storage, 5x for below */
  
  sk = (float ***)myalloc(xn*sizeof(float **));
  skn = (int **)myalloc(xn*sizeof(int *));
  for(i=0;i<xn;i++){
    sk[i] = (float **)myalloc(yn*sizeof(float *));
    skn[i] = (int *)myalloc(xn*sizeof(int));
    for(j=0;j<yn;j++){
      ss = s[i][j];
      n = cnt[i][j];
      tn = 0;  /* Number of spikes kept */
      for(k=0;k<n;k++){
	if (pop_util_ran2(&seed) <= pkeep){  /* Keep this spike */
	  ts[tn] = ss[k];
	  tn += 1;
	}
      }
      sk[i][j] = copy_farray(ts,tn);
      skn[i][j] = tn;
    }
  }

  /*printf(" CELL 30-32 has %d keep spikes\n",skn[30][32]);*/

  /*** 2.  Add spikes from neighbors ***/
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){ /* For each cell */
      tn = 0;  /* Total spikes so far for this cell */

      d0x = i - (int)((corr_d - 1)/2);
      d0y = j - (int)((corr_d - 1)/2);
      for(ii=d0x;ii<(d0x+corr_d);ii++){
	for(jj=d0y;jj<(d0y+corr_d);jj++){  /* For each cell in neigborhood */
	  if ((ii > 0) && (jj > 0) && (ii < (xn-1)) && (jj < (yn-1))){
	    ss = sk[ii][jj];
	    n = skn[ii][jj];
	    if ((ii == i) && (jj == j)){
	      for(k=0;k<n;k++){  /* Keep all self-spikes */
		ts[tn] = ss[k];
		tn += 1;
	      }
	    }else{
	      if (tsd > 0.0){
		for(k=0;k<n;k++){
		  if (pop_util_ran2(&seed) <= corr_p){  /* Keep this spike */
		    ts[tn] = ss[k] + tsd*pop_util_gasdev(&seed); /* jitter */
		    if ((ts[tn] > 0.0) && (ts[tn] < tmax))
		      tn += 1;
		  }
		}
	      }else{
		for(k=0;k<n;k++){
		  if (pop_util_ran2(&seed) <= corr_p){  /* no jitter */
		    ts[tn] = ss[k];
		    tn += 1;
		  }
		}
	      }
	    }
	  }
	}
      }

      /*
	if ((i==30)&&(j==32))
	printf(" CELL 30-32 had %d spikes,  Now has %d\n",cnt[30][32],tn);*/

      if (tn > 1)
	sort_farray(ts,tn);
      
      myfree(s[i][j]);
      s[i][j] = copy_farray(ts,tn);
      cnt[i][j] = tn;
    }
  }
  myfree(ts);

  free_3d_farray(sk,xn,yn,0); /* last value is not used */
  free_2d_iarray(skn,xn);

  /*
  printf("\n");
  printf("SPIKES FOR 30,32:\n");
  for(i=0;i<cnt[30][32];i++)
    printf(" %.1f",s[30][32][i]);
  printf("\n");
  printf("\n");*/

}
/**************************************-**************************************/
/*                                                                           */
/*                           POP_UTIL_INIT_BG_SPIKES                         */
/*                                                                           */
/*****************************************************************************/
void pop_util_init_bg_spikes(lay,pseed,tn,tscale)
     struct pop_layer *lay;    // Layer
     int *pseed;
     int tn;
     float tscale;
{
  int i,j,k;
  int n,s1n,s2n,s3n,done,scx_n,sci_n,pp_n;
  float duration,rmu,rsigma,sampling,f,sr1,*s1,*s2,*s3,cfrac_x,cfrac_i;
  float rate_x,rate_i,rate_sx,rate_si,rate_px,rate_pi,*scx,*sci,*pp_s,*pp_i;
  struct pop_cell *c;
  struct onode *bgo;

  //
  //  Find 'bg' <input> in 'inlist'
  //  WYETH - Here we assume there is only one input of type 'bg'
  //
  done = 0;
  i = 0;
  bgo = NULL;
  while(done == 0){
    if (onode_test_chr(lay->inlist[i],"type","bg")==1){
      bgo = lay->inlist[i];
      done = 1;
    }else{
      i += 1;
      if (i >= lay->ninlist)
	done = 1;
    }
  }
  //  WYETH - if 'bgo' is NULL, should we just save time and exit here?
  //  Or, is it possible that there could be bg params in the cells?



  rmu = 0.0;
  rsigma = 0.0;
  n = 1;
  sampling = 1.0; // Time in seconds

  duration = (float)tn * tscale * sampling;

  if (*pseed > 0)
    *pseed *= -1;


  //
  //  Check for correlated bg spikes
  //
  scx = sci = NULL;         // Shared spike trains
  scx_n = sci_n = 0;        // Number of spikes in each shared train
  rate_sx = rate_si = 0.0;  // rates of shared trains
  rate_px = rate_pi = 0.0;  // rates of private trains
  rate_x = rate_i = 0.0;

  if (bgo != NULL){


    cfrac_x = onode_getpar_flt_dflt(bgo,"ex_corr_frac",0.0);
    if (cfrac_x > 0.0){
      rate_x = onode_getpar_flt_exit(bgo,"ex_rate");

      rate_sx = rate_x * cfrac_x; // Rate for shared spikes
      rate_px = rate_x - rate_sx; // Rate for private spikes
      make_poisson_float_spikes_new(&scx,&scx_n,0.0,duration,sampling,rate_sx,
				    rmu,rsigma,n,pseed);
    }

    cfrac_i = onode_getpar_flt_dflt(bgo,"in_corr_frac",0.0);
    if (cfrac_i > 0.0){
      rate_i = onode_getpar_flt_exit(bgo,"in_rate");

      rate_si = rate_i * cfrac_i; // Rate for shared spikes
      rate_pi = rate_i - rate_si; // Rate for private spikes
      make_poisson_float_spikes_new(&sci,&sci_n,0.0,duration,sampling,rate_si,
				    rmu,rsigma,n,pseed);
    }

    /*
    printf("  rate_sx = %f\n",rate_sx);
    printf("  rate_px = %f\n",rate_px);
    printf("  rate_x  = %f\n",rate_x);
    printf("  scx_n   = %d\n",scx_n);

    printf("  rate_si = %f\n",rate_si);
    printf("  rate_pi = %f\n",rate_pi);
    printf("  rate_i  = %f\n",rate_i);
    printf("  sci_n   = %d\n",sci_n);
    */
  }


  // Nullify all spike trains
  for(i=0;i<lay->xn;i++){
    for(j=0;j<lay->yn;j++){
      for(k=0;k<lay->zn;k++){ // For each cell
	c = &(lay->c[i][j][k]);
	if (c->bg_x_s != NULL){
	  myfree(c->bg_x_s);
	  c->bg_x_s = NULL;
	}
	if (c->bg_i_s != NULL){
	  myfree(c->bg_i_s);
	  c->bg_i_s = NULL;
	}
	c->bg_x_n = -1;
	c->bg_i_n = -1;
      }
    }
  }

  for(i=0;i<lay->xn;i++){
    for(j=0;j<lay->yn;j++){
      for(k=0;k<lay->zn;k++){ // For each cell
	c = &(lay->c[i][j][k]);
	//duration = (float)tn * tscale * sampling;

	// EX
	if (c->bg_x_rate > 0.0){ // Make bg spikes
	  if (c->bg_x_corr_f > 0.0){ // We have to corr w/ another cell
	    //printf("POP_UTIL_INIT_BG_SPIKES: f=  %f\n",c->bg_x_corr_f);

	    if (c->bg_x_corr_c == NULL){ // No particular cell to corr with
	      //
	      //  Make noise shared across entire population
	      //

	      // Get spike train with private rate 'rate_px'
	      make_poisson_float_spikes_new(&pp_s,&pp_n,0.0,
					    duration,sampling,rate_px,
					    rmu,rsigma,n,pseed);

	      c->bg_x_s = merge_spike_arrays_float(pp_s,pp_n,scx,scx_n);
	      c->bg_x_n = pp_n + scx_n;
	      if (pp_s != NULL)
		myfree(pp_s);

	      //printf(" n-private = %d\n",pp_n);

	    }else{
	      //
	      //  This was written to allow correlated spikes with one
	      //  other cell in particular, and I think this was used
	      //  via old 'customize' code that is now commented out?
	      //
	      if (c->bg_x_corr_c->bg_x_n >= 0){ // Other has spikes
		
		f = c->bg_x_corr_f;
		s1 = c->bg_x_corr_c->bg_x_s;
		s1n = c->bg_x_corr_c->bg_x_n;
		spike_util_random_select_spikes(s1,s1n,f,*pseed,&s2,&s2n);

		sr1 = (1.0 - f) * c->bg_x_rate;
		if (sr1 > 0.0){
		  make_poisson_float_spikes_new(&s3,&s3n,0.0,duration,sampling,
						sr1,rmu,rsigma,n,pseed);
		}else{
		  s3 = NULL;
		  s3n = 0;
		}
		c->bg_x_s = merge_spike_arrays_float(s2,s2n,s3,s3n);
		c->bg_x_n = s2n + s3n;

		myfree(s2);
		if (s3 != NULL)
		  myfree(s3);
	      }else{
		make_poisson_float_spikes_new(&(c->bg_x_s),&(c->bg_x_n),0.0,
					      duration,sampling,c->bg_x_rate,
					      rmu,rsigma,n,pseed);
	      }
	    }
	  }else{
	    make_poisson_float_spikes_new(&(c->bg_x_s),&(c->bg_x_n),0.0,
					  duration,sampling,c->bg_x_rate,
					  rmu,rsigma,n,pseed);
	  }
	  c->bg_x_k = 0; // Index of next spike to use
	}

	// IN
	if (c->bg_i_rate > 0.0){
	  if (c->bg_i_corr_f > 0.0){ // We have to corr w/ another cell

	    if (c->bg_i_corr_c == NULL){ // No particular cell to corr with
	      //
	      //  Make noise shared across entire population
	      //

	      // Get spike train with private rate 'rate_pi'
	      make_poisson_float_spikes_new(&pp_s,&pp_n,0.0,
					    duration,sampling,rate_pi,
					    rmu,rsigma,n,pseed);

	      c->bg_i_s = merge_spike_arrays_float(pp_s,pp_n,sci,sci_n);
	      c->bg_i_n = pp_n + sci_n;
	      if (pp_s != NULL)
		myfree(pp_s);

	      //printf(" n-private = %d\n",pp_n);
	      //exit(0);

	    }else{
	      if (c->bg_i_corr_c->bg_i_n >= 0){ // Other has spikes

		f = c->bg_i_corr_f;
		s1 = c->bg_i_corr_c->bg_i_s;
		s1n = c->bg_i_corr_c->bg_i_n;
		spike_util_random_select_spikes(s1,s1n,f,*pseed,&s2,&s2n);
	      
		sr1 = (1.0 - f) * c->bg_i_rate;
		if (sr1 > 0.0){
		  make_poisson_float_spikes_new(&s3,&s3n,0.0,duration,sampling,
						sr1,rmu,rsigma,n,pseed);
		}else{
		  s3 = NULL;
		  s3n = 0;
		}
	      
		c->bg_i_s = merge_spike_arrays_float(s2,s2n,s3,s3n);
		c->bg_i_n = s2n + s3n;
		myfree(s2);
		if (s3 != NULL)
		  myfree(s3);
	      }else{
		make_poisson_float_spikes_new(&(c->bg_i_s),&(c->bg_i_n),0.0,
					      duration,sampling,c->bg_i_rate,
					      rmu,rsigma,n,pseed);
	      }
	    }
	  }else{
	    make_poisson_float_spikes_new(&(c->bg_i_s),&(c->bg_i_n),0.0,
					  duration,sampling,c->bg_i_rate,
					  rmu,rsigma,n,pseed);
	  }
	  c->bg_i_k = 0; // Index of next spike to use
	}
      }
    }
  }

  // Free any shared spike trains
  if (scx != NULL)  myfree(scx);
  if (sci != NULL)  myfree(sci);
}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_UTIL_INIT_BG_RATE                          */
/*                                                                           */
/*  Convert BG rate to a prob. value.                                        */
/*                                                                           */
/*****************************************************************************/
void pop_util_init_bg_rate(laylist,n)
     struct pop_layer **laylist;
     int n;
{
  int i,j,k,l;
  struct pop_cell *c;
  struct pop_layer *lay;

  printf("******************* POP_UTIL_INIT_BG_RATE\n");

  // Initialize simulation state and storage for all cells, all layers
  for(l=0;l<n;l++){ // For each layer
    lay = laylist[l];
    if (lay->runflag == 1){
      for(i=0;i<lay->xn;i++){
	for(j=0;j<lay->yn;j++){
	  for(k=0;k<lay->zn;k++){ // For each cell
	    c = &(lay->c[i][j][k]);

	    if (c->bg_x_rate > 0.0){
	      c->bg_x_rate /= 1000.0;
	    }
	    if (c->bg_i_rate > 0.0){
	      c->bg_i_rate /= 1000.0;
	    }
	  }
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                POP_INPUT_01                               */
/*                                                                           */
/*  Advance the state of the inputs for cell "c" based on incoming spikes.   */
/*  State is advanced in regular steps of size "DT_INPUT".                   */
/*                                                                           */
/*  NOTES                                                                    */
/*  - Spikes in the interval [T , T+DT) get used at time T+DT.               */
/*                                                                           */
/*****************************************************************************/
void pop_input_01(c,t1)
     struct pop_cell *c;   // Post-syn cell
     float t1;             // Advance state up to this time (sec)
{
  int k;
  int gi,gn,done,divide_flag,stype;
  float x_tau,x_aetau,x_tau2inv,x_totau;
  float i_tau,i_aetau,i_tau2inv,i_totau;
  float x_g_new,x_g_old,x_gd_new,x_gd_old,*x_gdata;
  float i_g_new,i_g_old,i_gd_new,i_gd_old,*i_gdata;
  float stor_x,stor_i,tval,dg_ex,dg_in,old_tsav;
  double t,dt,t0,sdt;
  struct pop_circbuf *x_sbuf,*i_sbuf,*x_gbuf,*i_gbuf;

  // mult pair processing
  struct pop_syn *syn_ptr;
  float tot,a,b,c_a,c_b;
  struct pop_circbuf *cb;

  old_tsav = -1.0;

  x_tau = c->ex1tau;  // PSC time const
  x_aetau = c->ex1amp * M_E / x_tau;
  x_tau2inv = 1.0/(x_tau*x_tau);
  x_totau = 2.0/x_tau;

  i_tau = c->in1tau;  // PSC time const
  i_aetau = c->in1amp * M_E / i_tau;
  i_tau2inv = 1.0/(i_tau*i_tau);
  i_totau = 2.0/i_tau;

  // The following 4 variables hold the ongoing, accumulated 'g' and its
  // derivative, which are used to allow g to peak and decay, thus avoiding
  // convolution with a mask.
  //
  x_g_old = c->gex1;     // Old value of g
  x_gd_old = c->gdex1;   // Old value of g deriv

  i_g_old = c->gin1;     // Old value of g
  i_gd_old = c->gdin1;   // Old value of g deriv

  x_sbuf = c->ex1s;      // Pointer to spike buffer
  i_sbuf = c->in1s;      // Pointer to spike buffer
  x_gbuf = c->ex1g;      // Pointer to g buffer
  i_gbuf = c->in1g;      // Pointer to g buffer
  x_gdata = x_gbuf->d;   // This receives the current value
  i_gdata = i_gbuf->d;   // This receives the current value

  dt = x_gbuf->dt;       // Same for inhib
  gn = x_gbuf->n;        // Same for inhib
  sdt = x_sbuf->dt;      // On 2010 Jan, was 0.0005 = 1/2 ms

  gi = x_gbuf->i;        // may be modified, same for inhib.
  t0 = x_sbuf->t + sdt;  // may be modified, same for inhib.
  t = c->gex1t;          // may be modified, same for inhib.

  if (c->savg){ // If we're saving g values for this cell
    if (c->gtex1n > 0){  // If we have some values saved already
      stor_x = c->gtex1[c->gtex1n - 1];  // take last saved ex as reference
      stor_i = c->gtin1[c->gtex1n - 1];  // take last saved in as reference
    }else
      // nothing stored yet, so use -1 to force storage of next value
      stor_x = stor_i = -1.0;
  }

  while(t < t1){
    t += dt;
    x_g_new = x_g_old + x_gd_old * dt;
    x_gd_new = x_gd_old - (x_totau * x_gd_old + x_tau2inv * x_g_old) * dt;

    i_g_new = i_g_old + i_gd_old * dt;
    i_gd_new = i_gd_old - (i_totau * i_gd_old + i_tau2inv * i_g_old) * dt;

    if (t >= t0){ // Add in new spikes at spike-store time origin
      k = x_sbuf->i; // Pointer to current time in spike buffer

      // Add Excitatory Background spikes at current time
      done = 0;
      while(!done){
	if (c->bg_x_k >= c->bg_x_n){
	  done = 1; // No bg spikes, or none left to be added
	}else{
	  //printf("HERE spike %f\n",c->bg_x_s[c->bg_x_k]);
	  if (t >= c->bg_x_s[c->bg_x_k]){
	    x_sbuf->d[k] += c->bg_x_amp;  // Add one spike here
	    c->bg_x_k += 1;
	  }else{
	    done = 1;
	  }
	}
      }

      // Add Inhibitory Background spikes at current time
      done = 0;
      while(!done){
	if (c->bg_i_k >= c->bg_i_n){
	  done = 1; // No bg spikes, or none left to be added
	}else{
	  if (t >= c->bg_i_s[c->bg_i_k]){
	    i_sbuf->d[k] += c->bg_i_amp;  // Add one spike here
	    c->bg_i_k += 1;
	  }else{
	    done = 1;
	  }
	}
      }

      // Add total spike input at this time
      x_gd_new += x_sbuf->d[k] * x_aetau;
      x_sbuf->d[k] = 0.0;   // Set spike count to zero after use

      i_gd_new += i_sbuf->d[k] * i_aetau;
      i_sbuf->d[k] = 0.0;   // Set spike count to zero after use

      k += 1;               // Advance index 1 unit around the circle
      if (k >= x_sbuf->n)
	k = 0;
      x_sbuf->i = k;

      t0 += sdt;            // t0 is sdt ahead of time origin
    }

    // Circular storage of g
    gi += 1;
    if (gi == gn)
      gi = 0;
    x_gdata[gi] = x_g_new;
    i_gdata[gi] = i_g_new;


    //
    //  Special mult pair syn processing
    //  Look at all dendritic subunits and add the contribution to g.
    //
    //  Note,  dt = 0.00005  =  1/20 msec
    //
    if (c->syn_proc_code == 1){

      tot = 0.0;
      syn_ptr = c->in;
      stype = -1;
      while (syn_ptr != NULL){
	if (syn_ptr->si != NULL){
	  if ((syn_ptr->si->syn_code == 100002) && (syn_ptr->si->ci == 0)){

	    if (stype == -1)
	      stype = syn_ptr->stype;
	    else{
	      if (syn_ptr->stype != stype){
		printf("*** old stype:  %d   now found %d\n",stype,
		       syn_ptr->stype);
		exit_error("POP_INPUT_01","syn_ptr type mismatch");
	      }
	    }

	    // Get constants, assumed to be identical for all synapses
	    c_a = syn_ptr->si->siu->s002->msi->a;
	    c_b = syn_ptr->si->siu->s002->msi->b;

	    if (syn_ptr->si->siu->s002->msi->si_comp == 1)  // Multiply
	      divide_flag = 0;
	    else
	      divide_flag = 1;  // Divide


	    // WYETH - THIS COULD BE MADE A LOT FASTER, by keeping track here
	    //  of all the values in 'pop_util_circbuf_val_at_t'.

	    cb = syn_ptr->si->siu->s002->wt1;
	    a = pop_util_circbuf_val_at_t(cb,t);
	    cb = syn_ptr->si->siu->s002->wt2;
	    b = pop_util_circbuf_val_at_t(cb,t);

	    /*** DEBUG
	    if ((t > 99.288) && (strcmp(syn_ptr->pre->name,"exs_5_5_1")==0)){
	      sleep(1);
	      printf(" ******* HERE WYETH found unit 5 5 1\n");
	      {
		int i;
		float d1[100],d2[100];

		printf("cb->n = %d\n",cb->n);
		printf("cb->dt = %f\n",(float)cb->dt);
	      
		for(i=0;i<100;i++){
		  cb = syn_ptr->si->siu->s002->wt1;
		  d1[i] = pop_util_circbuf_val_at_t(cb,t+(double)i/1000.0);
		  cb = syn_ptr->si->siu->s002->wt2;
		  d2[i] = pop_util_circbuf_val_at_t(cb,t+(double)i/1000.0);
		}
		append_farray_plot("cdump.pl","cb1",d1,100,1);
		append_farray_plot("cdump.pl","cb2",d2,100,1);
	      }
	      exit(0);
	    }
	    ***/

	    //printf("syn_ptr->w = %f\n",syn_ptr->w);

	    if (divide_flag == 1)
	      tval = a/(c_a + c_b * b);    // Divide
	    else
	      tval = a*(c_a + c_b * b);    // Multiply
	    if (tval < 0.0)
	      tval = 0.0;

	    tot += syn_ptr->w * tval; // Use synaptic weight 'w'
	  }
	}
	syn_ptr = syn_ptr->post_next;
      }
      if (stype == 1)
	x_gdata[gi] += tot;  // Add to total stored g (ex)
      else
	i_gdata[gi] += tot;  // Add to total stored g (in)
    }

    // Full time storage of dynamic input g's
    if (c->savg){
      // Only store the new g values if they have changed significantly
      // so that we don't waste memory if nothing is changing.

      dg_ex = (x_gdata[gi] - stor_x)/(x_gdata[gi] + 0.001);
      dg_in = (i_gdata[gi] - stor_i)/(i_gdata[gi] + 0.001);

      if ((dg_ex > 0.02) || (dg_ex < -0.02) ||
	  (dg_in > 0.02) || (dg_in < -0.02) || (t - old_tsav > 0.001)){
	/*if (((x_gdata[gi] - stor_x) > stor_eps) ||
	  ((x_gdata[gi] - stor_x) < -stor_eps)|| *** OLD WAY ****
	  ((i_gdata[gi] - stor_i) > stor_eps) || *** OLD WAY ****
	  ((i_gdata[gi] - stor_i) < -stor_eps)){*/
	c->gtex1[c->gtex1n] = x_gdata[gi];   // save excitatory g
	c->gtin1[c->gtex1n] = i_gdata[gi];   // save inhibitory g
	c->gtex1t[c->gtex1n] = t*1000.0; // Store time (ms) for this sample
	c->gtex1n += 1;
	if (c->gtex1n  > c->gtex1max){  // Storage over-run
	  printf("%d %d\n",c->gtex1n,c->gtex1max);
	  exit_error("POP_INPUT_01","Too many stored values");
	}
	stor_x = x_gdata[gi];  // Keep track of the last stored values
	stor_i = i_gdata[gi];
	old_tsav = t;
      }
    }

    x_g_old = x_g_new;
    x_gd_old = x_gd_new;
    i_g_old = i_g_new;
    i_gd_old = i_gd_new;
  }

  x_gbuf->i = gi;
  x_gbuf->t = t;  // Time at gbuf->i

  c->gex1t = t;
  x_sbuf->t = t0 - sdt;

  c->gex1 = x_g_new;
  c->gdex1 = x_gd_new;
  c->gin1 = i_g_new;
  c->gdin1 = i_gd_new;

  i_sbuf->i = x_sbuf->i;  // These value are copies of excitatory values
  i_sbuf->t = x_sbuf->t;
}
/**************************************-**************************************/
/*                                                                           */
/*                                POP_INPUT_01A                              */
/*                                                                           */
/*  Like "pop_input_01" but generates BG spikes on the fly.                  */
/*                                                                           */
/*****************************************************************************/
void pop_input_01a(c,t1)
     struct pop_cell *c;   // Post-syn cell
     float t1;             // Advance state up to this time (sec)
{
  int k;
  int gi,gn,done,divide_flag,stype;
  float x_tau,x_aetau,x_tau2inv,x_totau;
  float i_tau,i_aetau,i_tau2inv,i_totau;
  float x_g_new,x_g_old,x_gd_new,x_gd_old,*x_gdata;
  float i_g_new,i_g_old,i_gd_new,i_gd_old,*i_gdata;
  float stor_x,stor_i,tval,dg_ex,dg_in,old_tsav;
  double t,dt,t0,sdt;
  struct pop_circbuf *x_sbuf,*i_sbuf,*x_gbuf,*i_gbuf;

  // mult pair processing
  struct pop_syn *syn_ptr;
  float tot,a,b,c_a,c_b;
  struct pop_circbuf *cb;

  old_tsav = -1.0;

  x_tau = c->ex1tau;  // PSC time const
  x_aetau = c->ex1amp * M_E / x_tau;
  x_tau2inv = 1.0/(x_tau*x_tau);
  x_totau = 2.0/x_tau;

  i_tau = c->in1tau;  // PSC time const
  i_aetau = c->in1amp * M_E / i_tau;
  i_tau2inv = 1.0/(i_tau*i_tau);
  i_totau = 2.0/i_tau;

  // The following 4 variables hold the ongoing, accumulated 'g' and its
  // derivative, which are used to allow g to peak and decay, thus avoiding
  // convolution with a mask.
  //
  x_g_old = c->gex1;     // Old value of g
  x_gd_old = c->gdex1;   // Old value of g deriv

  i_g_old = c->gin1;     // Old value of g
  i_gd_old = c->gdin1;   // Old value of g deriv

  x_sbuf = c->ex1s;      // Pointer to spike buffer
  i_sbuf = c->in1s;      // Pointer to spike buffer
  x_gbuf = c->ex1g;      // Pointer to g buffer
  i_gbuf = c->in1g;      // Pointer to g buffer
  x_gdata = x_gbuf->d;   // This receives the current value
  i_gdata = i_gbuf->d;   // This receives the current value

  dt = x_gbuf->dt;       // Same for inhib
  gn = x_gbuf->n;        // Same for inhib
  sdt = x_sbuf->dt;      // On 2010 Jan, was 0.0005 = 1/2 ms

  gi = x_gbuf->i;        // may be modified, same for inhib.
  t0 = x_sbuf->t + sdt;  // may be modified, same for inhib.
  t = c->gex1t;          // may be modified, same for inhib.

  if (c->savg){ // If we're saving g values for this cell
    if (c->gtex1n > 0){  // If we have some values saved already
      stor_x = c->gtex1[c->gtex1n - 1];  // take last saved ex as reference
      stor_i = c->gtin1[c->gtex1n - 1];  // take last saved in as reference
    }else
      // nothing stored yet, so use -1 to force storage of next value
      stor_x = stor_i = -1.0;
  }


  /*
  if (strcmp(c->name,"exs_10_10_0")==0){
    printf("t = %f  t1 = %f   dt = %f\n",(float)t,t1,dt);
    printf("c->bg_x_rate = %f\n",c->bg_x_rate);
    printf("c->bg_i_rate = %f\n",c->bg_i_rate);
    }*/



  //
  //  WYETH - Use global seed;
  //
  k = x_sbuf->i; // Pointer to current time in spike buffer
  if (pop_util_ran2(&globseed) < c->bg_x_rate){
    x_sbuf->d[k] += c->bg_x_amp;  // Add one spike here
  }
  if (pop_util_ran2(&globseed) < c->bg_i_rate){
    i_sbuf->d[k] += c->bg_i_amp;  // Add one spike here
  }


  while(t < t1){
    t += dt;
    x_g_new = x_g_old + x_gd_old * dt;
    x_gd_new = x_gd_old - (x_totau * x_gd_old + x_tau2inv * x_g_old) * dt;

    i_g_new = i_g_old + i_gd_old * dt;
    i_gd_new = i_gd_old - (i_totau * i_gd_old + i_tau2inv * i_g_old) * dt;

    if (t >= t0){ // Add in new spikes at spike-store time origin
      k = x_sbuf->i; // Pointer to current time in spike buffer

      // Add total spike input at this time
      x_gd_new += x_sbuf->d[k] * x_aetau;
      x_sbuf->d[k] = 0.0;   // Set spike count to zero after use

      i_gd_new += i_sbuf->d[k] * i_aetau;
      i_sbuf->d[k] = 0.0;   // Set spike count to zero after use

      k += 1;               // Advance index 1 unit around the circle
      if (k >= x_sbuf->n)
	k = 0;
      x_sbuf->i = k;

      t0 += sdt;            // t0 is sdt ahead of time origin
    }

    // Circular storage of g
    gi += 1;
    if (gi == gn)
      gi = 0;
    x_gdata[gi] = x_g_new;
    i_gdata[gi] = i_g_new;


    //
    //  Special mult pair syn processing
    //  Look at all dendritic subunits and add the contribution to g.
    //
    //  Note,  dt = 0.00005  =  1/20 msec
    //
    if (c->syn_proc_code == 1){

      tot = 0.0;
      syn_ptr = c->in;
      stype = -1;
      while (syn_ptr != NULL){
	if (syn_ptr->si != NULL){
	  if ((syn_ptr->si->syn_code == 100002) && (syn_ptr->si->ci == 0)){

	    if (stype == -1)
	      stype = syn_ptr->stype;
	    else{
	      if (syn_ptr->stype != stype){
		printf("*** old stype:  %d   now found %d\n",stype,
		       syn_ptr->stype);
		exit_error("POP_INPUT_01","syn_ptr type mismatch");
	      }
	    }

	    // Get constants, assumed to be identical for all synapses
	    c_a = syn_ptr->si->siu->s002->msi->a;
	    c_b = syn_ptr->si->siu->s002->msi->b;

	    if (syn_ptr->si->siu->s002->msi->si_comp == 1)  // Multiply
	      divide_flag = 0;
	    else
	      divide_flag = 1;  // Divide


	    // WYETH - THIS COULD BE MADE A LOT FASTER, by keeping track here
	    //  of all the values in 'pop_util_circbuf_val_at_t'.

	    cb = syn_ptr->si->siu->s002->wt1;
	    a = pop_util_circbuf_val_at_t(cb,t);
	    cb = syn_ptr->si->siu->s002->wt2;
	    b = pop_util_circbuf_val_at_t(cb,t);

	    /*** DEBUG
	    if ((t > 99.288) && (strcmp(syn_ptr->pre->name,"exs_5_5_1")==0)){
	      sleep(1);
	      printf(" ******* HERE WYETH found unit 5 5 1\n");
	      {
		int i;
		float d1[100],d2[100];

		printf("cb->n = %d\n",cb->n);
		printf("cb->dt = %f\n",(float)cb->dt);
	      
		for(i=0;i<100;i++){
		  cb = syn_ptr->si->siu->s002->wt1;
		  d1[i] = pop_util_circbuf_val_at_t(cb,t+(double)i/1000.0);
		  cb = syn_ptr->si->siu->s002->wt2;
		  d2[i] = pop_util_circbuf_val_at_t(cb,t+(double)i/1000.0);
		}
		append_farray_plot("cdump.pl","cb1",d1,100,1);
		append_farray_plot("cdump.pl","cb2",d2,100,1);
	      }
	      exit(0);
	    }
	    ***/

	    //printf("syn_ptr->w = %f\n",syn_ptr->w);

	    if (divide_flag == 1)
	      tval = a/(c_a + c_b * b);    // Divide
	    else
	      tval = a*(c_a + c_b * b);    // Multiply
	    if (tval < 0.0)
	      tval = 0.0;

	    tot += syn_ptr->w * tval; // Use synaptic weight 'w'
	  }
	}
	syn_ptr = syn_ptr->post_next;
      }
      if (stype == 1)
	x_gdata[gi] += tot;  // Add to total stored g (ex)
      else
	i_gdata[gi] += tot;  // Add to total stored g (in)
    }

    // Full time storage of dynamic input g's
    if (c->savg){
      // Only store the new g values if they have changed significantly
      // so that we don't waste memory if nothing is changing.

      dg_ex = (x_gdata[gi] - stor_x)/(x_gdata[gi] + 0.001);
      dg_in = (i_gdata[gi] - stor_i)/(i_gdata[gi] + 0.001);

      if ((dg_ex > 0.02) || (dg_ex < -0.02) ||
	  (dg_in > 0.02) || (dg_in < -0.02) || (t - old_tsav > 0.001)){
	/*if (((x_gdata[gi] - stor_x) > stor_eps) ||
	  ((x_gdata[gi] - stor_x) < -stor_eps)|| *** OLD WAY ****
	  ((i_gdata[gi] - stor_i) > stor_eps) || *** OLD WAY ****
	  ((i_gdata[gi] - stor_i) < -stor_eps)){*/
	c->gtex1[c->gtex1n] = x_gdata[gi];   // save excitatory g
	c->gtin1[c->gtex1n] = i_gdata[gi];   // save inhibitory g
	c->gtex1t[c->gtex1n] = t*1000.0; // Store time (ms) for this sample
	c->gtex1n += 1;
	if (c->gtex1n  > c->gtex1max){  // Storage over-run
	  printf("%d %d\n",c->gtex1n,c->gtex1max);
	  exit_error("POP_INPUT_01","Too many stored values");
	}
	stor_x = x_gdata[gi];  // Keep track of the last stored values
	stor_i = i_gdata[gi];
	old_tsav = t;
      }
    }

    x_g_old = x_g_new;
    x_gd_old = x_gd_new;
    i_g_old = i_g_new;
    i_gd_old = i_gd_new;
  }

  x_gbuf->i = gi;
  x_gbuf->t = t;  // Time at gbuf->i

  c->gex1t = t;
  x_sbuf->t = t0 - sdt;

  c->gex1 = x_g_new;
  c->gdex1 = x_gd_new;
  c->gin1 = i_g_new;
  c->gdin1 = i_gd_new;

  i_sbuf->i = x_sbuf->i;  // These value are copies of excitatory values
  i_sbuf->t = x_sbuf->t;
}
/**************************************-**************************************/
/*                                                                           */
/*                                POP_INPUT_02                               */
/*                                                                           */
/*  Advance the state of the inputs for cell "c" based on incoming spikes.   */
/*  State is advanced in regular steps of size "DT_INPUT".                   */
/*                                                                           */
/*  NOTES                                                                    */
/*  - Spikes in the interval [T , T+DT) get used at time T+DT.               */
/*                                                                           */
/*****************************************************************************/
void pop_input_02(c,t1,lay)
     struct pop_cell *c;
     float t1;             /* Advance state up to this time (sec) */
     struct pop_layer *lay;
{
  int k;
  int gi,gn,done;

  float ex_tau,ex_aetau,ex_tau2inv,ex_totau;
  float in_tau,in_aetau,in_tau2inv,in_totau;

  float n1_aetau,n1_tau2inv,n1_totau;
  float n2_aetau,n2_tau2inv,n2_totau;
  float n3_aetau,n3_tau2inv,n3_totau;

  float ex_g_new,ex_g_old,ex_gd_new,ex_gd_old,*x_gdata;
  float in_g_new,in_g_old,in_gd_new,in_gd_old,*i_gdata;

  float n1_g_new,n1_g_old,n1_gd_new,n1_gd_old,*n_gdata;
  float n2_g_new,n2_g_old,n2_gd_new,n2_gd_old;
  float n3_g_new,n3_g_old,n3_gd_new,n3_gd_old;
  float nm_g_new;

  float stor_x,stor_i,stor_n,stor_eps;

  float xsc;
  double t,dt,t0,sdt;
  struct pop_circbuf *x_sbuf,*i_sbuf;

  stor_eps = 0.1;       // WYETH the value 0.05 is arbitrary

  // NMDA - sum of three alpha functions
  n1_aetau   = lay->nmda_p->nmda1_aetau;
  n1_tau2inv = lay->nmda_p->nmda1_tau2inv;
  n1_totau   = lay->nmda_p->nmda1_totau;
  n2_aetau   = lay->nmda_p->nmda2_aetau;
  n2_tau2inv = lay->nmda_p->nmda2_tau2inv;
  n2_totau   = lay->nmda_p->nmda2_totau;
  n3_aetau   = lay->nmda_p->nmda3_aetau;
  n3_tau2inv = lay->nmda_p->nmda3_tau2inv;
  n3_totau   = lay->nmda_p->nmda3_totau;

  ex_tau = c->ex1tau;  // PSC time const
  ex_aetau = c->ex1amp * M_E / ex_tau;
  ex_tau2inv = 1.0/(ex_tau*ex_tau);
  ex_totau = 2.0/ex_tau;

  in_tau = c->in1tau;  // PSC time const
  in_aetau = c->in1amp * M_E / in_tau;
  in_tau2inv = 1.0/(in_tau*in_tau);
  in_totau = 2.0/in_tau;

  n1_g_old  = c->gex1_nmda1;  // Old value of gNMDA
  n2_g_old  = c->gex1_nmda2;  // Old value of gNMDA
  n3_g_old  = c->gex1_nmda3;  // Old value of gNMDA
  n1_gd_old = c->gdex1_nmda1; // Old value of gNMDA deriv
  n2_gd_old = c->gdex1_nmda2; // Old value of gNMDA deriv
  n3_gd_old = c->gdex1_nmda3; // Old value of gNMDA deriv

  ex_g_old  = c->gex1;       // Old value of g
  ex_gd_old = c->gdex1;      // Old value of g deriv
  in_g_old  = c->gin1;       // Old value of g
  in_gd_old = c->gdin1;      // Old value of g deriv

  x_sbuf = c->ex1s; // Pointer to spike buffer
  i_sbuf = c->in1s; // Pointer to spike buffer

  n_gdata = c->ex1g_nmda->d; // Only need one, to store summed g NMDA
  x_gdata = c->ex1g->d;
  i_gdata = c->in1g->d;   // Pointer to g buffer data

  dt = c->ex1g->dt;       // Same for inhib
  gn = c->ex1g->n;        // Same for inhib
  sdt = x_sbuf->dt;       // Same for inhib

  gi = c->ex1g->i;        // may be modified, same for inhib.
  t0 = x_sbuf->t + sdt;   // may be modified, same for inhib.
  t = c->gex1t;           // may be modified, same for inhib.

  if (c->savg){
    if (c->gtex1n > 0){
      stor_x = c->gtex1[c->gtex1n - 1];
      stor_i = c->gtin1[c->gtex1n - 1];
      stor_n = c->gtex1_nmda[c->gtex1n - 1];
    }else
      stor_x = stor_i = stor_n = -1.0; // Force storage of first values
  }

  while(t < t1){
    t += dt;

    n1_g_new = n1_g_old + n1_gd_old * dt;
    n2_g_new = n2_g_old + n2_gd_old * dt;
    n3_g_new = n3_g_old + n3_gd_old * dt;
    ex_g_new = ex_g_old + ex_gd_old * dt;
    in_g_new = in_g_old + in_gd_old * dt;
    n1_gd_new = n1_gd_old - (n1_totau * n1_gd_old + n1_tau2inv * n1_g_old) *dt;
    n2_gd_new = n2_gd_old - (n2_totau * n2_gd_old + n2_tau2inv * n2_g_old) *dt;
    n3_gd_new = n3_gd_old - (n3_totau * n3_gd_old + n3_tau2inv * n3_g_old) *dt;
    ex_gd_new = ex_gd_old - (ex_totau * ex_gd_old + ex_tau2inv * ex_g_old) *dt;
    in_gd_new = in_gd_old - (in_totau * in_gd_old + in_tau2inv * in_g_old) *dt;

    if (t >= t0){ // Add in new spikes at spike-store time origin
      k = x_sbuf->i; // Pointer to current time in spike buffer

      // New way to do background spikes
      done = 0;
      while(!done){
	if (c->bg_x_k >= c->bg_x_n){
	  done = 1; // No bg spikes, or none left to be added
	}else{
	  // printf("HERE spike %f\n",c->bg_x_s[c->bg_x_k]);
	  if (t >= c->bg_x_s[c->bg_x_k]){
	    x_sbuf->d[k] += c->bg_x_amp;  // Add one spike here
	    c->bg_x_k += 1;
	  }else{
	    done = 1;
	  }
	}
      }
      done = 0;  // Now do inhib
      while(!done){
	if (c->bg_i_k >= c->bg_i_n){
	  done = 1; // No bg spikes, or none left to be added
	}else{
	  if (t >= c->bg_i_s[c->bg_i_k]){
	    i_sbuf->d[k] += c->bg_i_amp;  // Add one spike here
	    c->bg_i_k += 1;
	  }else{
	    done = 1;
	  }
	}
      }

      xsc = x_sbuf->d[k];
      n1_gd_new += xsc * n1_aetau;  // NMDA
      n2_gd_new += xsc * n2_aetau;
      n3_gd_new += xsc * n3_aetau;
      ex_gd_new += xsc * ex_aetau;  // AMPA
      x_sbuf->d[k] = 0.0;   // Set spike count to zero after use

      in_gd_new += i_sbuf->d[k] * in_aetau;
      i_sbuf->d[k] = 0.0;   // Set spike count to zero after use

      k += 1;               // Advance index 1 unit around the circle
      if (k >= x_sbuf->n)
	k = 0;
      x_sbuf->i = k;

      t0 += sdt;                // t0 is sdt ahead of time origin
    }

    // Full time storage of dynamic input g's
    if (c->savg){
      nm_g_new = n1_g_new + n2_g_new + n3_g_new;
      if (((ex_g_new - stor_x) > stor_eps)||((ex_g_new - stor_x) < -stor_eps)||
	  ((in_g_new - stor_i) > stor_eps)||((in_g_new - stor_i) < -stor_eps)||
	  ((nm_g_new - stor_n) > stor_eps)||((nm_g_new - stor_n) < -stor_eps)){
	
	c->gtex1[c->gtex1n]      = ex_g_new;
	c->gtex1_nmda[c->gtex1n] = nm_g_new;
	c->gtin1[c->gtex1n]      = in_g_new;
	c->gtex1t[c->gtex1n] = t*1000.0; // Time in ms
	c->gtex1n += 1;
	if (c->gtex1n  > c->gtex1max)
	  exit_error("POP_INPUT_02","Too many stored values");
	
	stor_x = ex_g_new;
	stor_i = in_g_new;
	stor_n = nm_g_new;
      }
    }

    // Circular storage of g
    gi += 1;
    if (gi == gn)
      gi = 0;
    n_gdata[gi] = n1_g_new + n2_g_new + n3_g_new;
    x_gdata[gi] = ex_g_new;
    i_gdata[gi] = in_g_new;

    n1_g_old = n1_g_new;
    n2_g_old = n2_g_new;
    n3_g_old = n3_g_new;
    n1_gd_old = n1_gd_new;
    n2_gd_old = n2_gd_new;
    n3_gd_old = n3_gd_new;
    ex_g_old = ex_g_new;
    ex_gd_old = ex_gd_new;
    in_g_old = in_g_new;
    in_gd_old = in_gd_new;
  }

  c->ex1g->i = gi;
  c->ex1g->t = t;  // Time at ex1g->i

  c->gex1t = t;
  x_sbuf->t = t0 - sdt;

  c->gex1 = ex_g_new;
  c->gdex1 = ex_gd_new;
  c->gin1 = in_g_new;
  c->gdin1 = in_gd_new;

  i_sbuf->i = x_sbuf->i;  // These value are copies of excitatory values
  i_sbuf->t = x_sbuf->t;
}
/**************************************-**************************************/
/*                                                                           */
/*                              POP_UTIL_DERIV_01                            */
/*                                                                           */
/*  Fill 'dydx'.  No change to 'x' or 'v'.                                   */
/*                                                                           */
/*****************************************************************************/
void pop_util_deriv_01(x,v,dydx,c)
     float x;      // Current time
     float v;      // Current voltage
     float *dydx;  // Derivative
     struct pop_cell *c;  // cell
{
  int k;
  int n01,i0,i1;
  float dt,t,g_ex,g_in,g_ad,tt,gexd,gind;
  struct ifc_param *p;
  struct pop_circbuf *gbuf,*i_gbuf;

  //
  //  Get value for dynamic synaptic inputs at time x, using interpolation
  //
  gbuf   = c->ex1g;
  i_gbuf = c->in1g;
  if ((x > gbuf->t) || (x < (gbuf->t - (gbuf->n-1)*gbuf->dt))){
    printf("DERIV:  x=%f gbuf->t=%f  gbuf->i=%d\n",x,gbuf->t,gbuf->i);
    exit_error("POP_UTIL_DERIV_01","x not covered by gbuf");
  }
  // Compute indices before and after desired t-value
  tt = (gbuf->t - x)/gbuf->dt;
  k = (int)tt;
  dt = (float)(k+1) - tt;  // from 0 to 1; fractional time to interpolate
  i1 = gbuf->i - k;
  if (i1 > 0)
    i0 = i1 - 1;
  else if (i1 == 0)
    i0 = gbuf->n - 1;
  else{
    i1 += gbuf->n;
    i0 = i1 - 1;
  }
  // Interpolate 'dt' of the way from 'i0' to 'i1'
  gexd =   gbuf->d[i0] + dt * (  gbuf->d[i1] -   gbuf->d[i0]);
  gind = i_gbuf->d[i0] + dt * (i_gbuf->d[i1] - i_gbuf->d[i0]);


  //
  // WYETH NEW - this part added for dynamic 'gad'; 19 June 2010
  //
  if (c->ad1g != NULL){
    g_ad = pop_util_circbuf_val_at_t_interp(c->ad1g,x);
  }else{
    g_ad = 0.0;
  }


  p = c->cifcp;
  n01 = c->n0 - 1;  // n0 is length of pre-computed conductance inputs

  if (n01 <= 0){ // If the cell has no such input ...
    g_ex = 0.0;
    g_in = 0.0;
  }else{
    // Interpolate to get g's
    t = x * c->samp0;   // Time in ms (float)
    k = (int)t;     // Time in ms (int)

    if (k < n01){ // Note, k does go up to 'c->n0'
      dt = t - (float)k;
      g_ex = c->gtx0[k] +  dt * (c->gtx0[k+1]-c->gtx0[k]);
      g_in = c->gti0[k] +  dt * (c->gti0[k+1]-c->gti0[k]);
    }else{
      g_ex = c->gtx0[n01];
      g_in = c->gti0[n01];
    }
  }
  g_ex += gexd;
  g_in += gind;

  *dydx = 1000.0/p->c_x * (g_ex        * (p->v_ex     - v) +
			   g_in        * (p->v_in     - v) +
			   p->g_leak_x * (p->v_leak_x - v) +
			   g_ad        * (p->v_ad     - v));

  //
  //  WYETH - stop if value goes out of bounds
  //
  /***
  if (!((*dydx > -100000.0) && (*dydx < 100000.0))){
    printf(" g_ex   = %f\n",g_ex);
    printf("  gexd  = %f\n",gexd);
    printf(" g_in   = %f\n",g_in);
    printf("  gind  = %f\n",gind);
    printf(" g_ad   = %f\n",g_ad);
    printf(" v      = %f\n",v);

    printf(" i0,i1 = %d,%d\n",i0,i1);
    printf("   gbuf->d[i0],[i1] = %f %f \n",gbuf->d[i0],gbuf->d[i1]);
    printf(" i_gbuf->d[i0],[i1] = %f %f \n",i_gbuf->d[i0],i_gbuf->d[i1]);

    printf(" dydx = %f\n",*dydx);
    printf(" c->name %s\n",c->name);
    exit_error("POP_UTIL_DERIV_01","dydx out of bounds");
  }
  ***/

}
/**************************************-**************************************/
/*                                                                           */
/*                              POP_UTIL_DERIV_02                            */
/*                                                                           */
/*  Two-compartment model.                                                   */
/*  Inputs 'y' and 'dydx' are arrays.                                        */
/*                                                                           */
/*****************************************************************************/
void pop_util_deriv_02(float x, float *y, float *dydx, struct pop_cell *c,
		       struct pop_layer *lay)
{
  int k;
  int n01,i0,i1;
  float dt,t,g_ex,g_nm,g_in,g_ad,tt,gexd,gnmd,gind,vs,vd,gt,rfact,ff;
  struct ifc_param *p;
  struct nmda_param *tp;
  struct pop_circbuf *gbuf,*i_gbuf,*n_gbuf;

  if (global_dump_flag)
    printf("GLOBAL DUMP--------------------------------\n");

  vs = y[1]; // V soma
  vd = y[2]; // V dendr

  /*printf("x=%f\n",x);*/

  // Get value for dynamic synaptic inputs at time x
  gbuf = c->ex1g;
  n_gbuf = c->ex1g_nmda;
  i_gbuf = c->in1g;
  if ((x > gbuf->t) || (x < (gbuf->t - (gbuf->n-1)*gbuf->dt))){
    printf("DERIV:  x=%f gbuf->t=%f  gbuf->i=%d\n",x,gbuf->t,gbuf->i);
    exit_error("POP_UTIL_DERIV_01","x not covered by gbuf");
  }
  /*printf("DERIV:  x=%f gbuf->t=%f  gbuf->i=%d\n",x,gbuf->t,gbuf->i);*/
  
  /*** Compute indices before and after desired t-value ***/
  tt = (gbuf->t - x)/gbuf->dt;
  /*printf("  tt= %f\n",tt);*/
  k = (int)tt;
  /*printf("  k= %d\n",k);*/
  dt = (float)(k+1) - tt;  /* from 0 to 1 */
  /*printf("  dt= %f\n",dt);*/
  i1 = gbuf->i - k;
  if (i1 > 0)
    i0 = i1 - 1;
  else if (i1 == 0)
    i0 = gbuf->n - 1;
  else{
    i1 += gbuf->n;
    i0 = i1 - 1;
  }
  gexd =   gbuf->d[i0] + dt * (  gbuf->d[i1] -   gbuf->d[i0]);
  gnmd = n_gbuf->d[i0] + dt * (n_gbuf->d[i1] - n_gbuf->d[i0]);
  gind = i_gbuf->d[i0] + dt * (i_gbuf->d[i1] - i_gbuf->d[i0]);

  //
  // WYETH NEW - this part added for dynamic 'gad'; 19 June 2010
  //
  if (c->ad1g != NULL){
    g_ad = pop_util_circbuf_val_at_t_interp(c->ad1g,x);
  }else{
    g_ad = 0.0;
  }


  /*printf("  i0 i1=  %d %d   ",i0,i1);*/
  /*printf("  (%f) %f (%f)\n",gbuf->d[i0],gexd,gbuf->d[i1]);*/

  /*
  if (gexd != 0){
    append_farray_plot("zz.gbuf.pl","gbuf",gbuf->d,gbuf->n,1);
    exit(0);
    }*/

  /*** Pre-computed conductances ***/
  n01 = c->n0 - 1;  /* n0 is length of pre-computed conductance inputs */
  if (n01 == 0){ /*** If the cell has no such input ... ***/
    g_ex = 0.0;
    g_nm = 0.0;
    g_in = 0.0;
  }else{
    // Interpolate to get g's
    t = x * c->samp0;   // Time in ms (float)
    k = (int)t;     // Time in ms (int)
    if (k < n01){ // Note, k does go up to 'c->n0'
      dt = t - (float)k;
      g_ex = c->gtx0[k] +  dt * (c->gtx0[k+1]-c->gtx0[k]);
      g_nm = c->gtn0[k] +  dt * (c->gtn0[k+1]-c->gtn0[k]);
      g_in = c->gti0[k] +  dt * (c->gti0[k+1]-c->gti0[k]);
    }else{
      g_ex = c->gtx0[n01];
      g_nm = c->gtn0[n01];
      g_in = c->gti0[n01];
    }
  }

  /*** Find NMDA scaling factor as a function of Vm ***/
  tp = lay->nmda_p;
  k = (vd - tp->nmda_v0f) * tp->nmda_ppmv;
  if ((k < 0) || (k >= tp->nmda_n)){
    printf("vd = %f  vs = %f\n",vd,vs);
    printf("k = %d   tp->nmda_n = %d\n",k,tp->nmda_n);
    exit_error("POP_UTIL_DERIV_O2","NMDA table index");
  }
  ff = tp->nmda_gv[k];
  /*printf("    NMDA voltage-dependent factor at %.2f mV is %f\n",v,ff);*/

  g_ex += gexd + ff*(gnmd + g_nm); /* Add in AMPA and NMDA, they share v_ex */
  g_in += gind;

  p = c->cifcp;
  gt = p->g_tran;

  /* SOMA */
  dydx[1] = 1000.0/p->c_x * (g_ex        * (p->v_ex     - vs) +
			     g_in        * (p->v_in     - vs) +
			     p->g_leak_x * (p->v_leak_x - vs) +
			     g_ad        * (p->v_ad     - vs) +
			     gt          * (vd - vs));

  // Refractory period
  if (x < c->reft)
    if (dydx[1] > 0.0){ // Only suppress *increases* in Vm
      rfact = 1.0 - (c->reft - x)/p->trefr_x_s;
      dydx[1] *= rfact;  // If in refractory period, do not change Vm
      if (global_dump_flag){
	printf("  Refraction rfact = %f\n",rfact);
	printf("    x < reft  -->   %f  < %f\n",x,c->reft);
      }
    }
  
  /* DENDR */
  dydx[2] = 1000.0/p->c_x * (g_ex        * (p->v_ex     - vd) +
			     g_in        * (p->v_in     - vd) +
			     p->g_leak_x * (p->v_leak_x - vd) +
			     gt          * (vs - vd));

  if (global_dump_flag){
    printf("  gex=%f gin=%f gad=%f glk=%f gt=%f\n",g_ex,g_in,g_ad,p->g_leak_x,
	   gt);
    printf("  1000.0/c_x  =  %f\n",1000.0/p->c_x);

    printf("  v_ex - vs  =  %f\n",p->v_ex     - vs);
    printf("  v_in - vs  =  %f\n",p->v_in     - vs);
    printf("  v_lk - vs  =  %f\n",p->v_leak_x - vs);
    printf("  v_ad - vs  =  %f\n",p->v_ad     - vs);
    printf("  vd   - vs  =  %f\n",vd          - vs);

    printf("  v_ex - vd  =  %f\n",p->v_ex     - vd);
    printf("  v_in - vd  =  %f\n",p->v_in     - vd);
    printf("  v_lk - vd  =  %f\n",p->v_leak_x - vd);
    printf("  vs   - vd  =  %f\n",vs          - vd);

    printf("  old vs = %f\n",vs);
    printf("  old vd = %f\n",vd);
    printf("  soma dydx = %f\n",dydx[1]);
    printf("  dend dydx = %f\n",dydx[2]);

    printf("  x (time) = %f\n",x);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                POP_RKCK_01                                */
/*                                                                           */
/*  Modified "rkck" (NumRecInC) to call a particular "derivs" routine, and   */
/*  to have "nvar" = 1.                                                      */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void pop_rkck_01(float y, float dydx, float x, float h, float *yout,
		 float *yerr, struct pop_cell *c)
{
  static float a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.0/14336.0;
  float dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  float ak2,ak3,ak4,ak5,ak6,ytemp;
  
  ytemp=y+b21*h*dydx;
  pop_util_deriv_01(x+a2*h,ytemp,&ak2,c);

  ytemp=y+h*(b31*dydx+b32*ak2);
  pop_util_deriv_01(x+a3*h,ytemp,&ak3,c);

  ytemp=y+h*(b41*dydx+b42*ak2+b43*ak3);
  pop_util_deriv_01(x+a4*h,ytemp,&ak4,c);

  ytemp=y+h*(b51*dydx+b52*ak2+b53*ak3+b54*ak4);
  pop_util_deriv_01(x+a5*h,ytemp,&ak5,c);

  ytemp=y+h*(b61*dydx+b62*ak2+b63*ak3+b64*ak4+b65*ak5);
  pop_util_deriv_01(x+a6*h,ytemp,&ak6,c);

  *yout=y+h*(c1*dydx+c3*ak3+c4*ak4+c6*ak6);

  *yerr=h*(dc1*dydx+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6);
}
/**************************************-**************************************/
/*                                                                           */
/*                                POP_RKCK_02                                */
/*                                                                           */
/*  Like NumRec RKCK, but added 'struct pop_cell *cp' as param and into      */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void pop_rkck_02(float y[], float dydx[], int n,float x, float h, float yout[],
		 float yerr[], struct pop_cell *cp, struct pop_layer *lay,
		 void (*derivs)(float, float [], float [], struct pop_cell *,
				struct pop_layer *))
{
  int i;
  static float a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.0/14336.0;
  float dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  float *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;

  ak2=vector(1,n);
  ak3=vector(1,n);
  ak4=vector(1,n);
  ak5=vector(1,n);
  ak6=vector(1,n);
  ytemp=vector(1,n);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+b21*h*dydx[i];
  (*derivs)(x+a2*h,ytemp,ak2,cp,lay);  /* Wyeth - added 'cp' */
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
  (*derivs)(x+a3*h,ytemp,ak3,cp,lay);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
  (*derivs)(x+a4*h,ytemp,ak4,cp,lay);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  (*derivs)(x+a5*h,ytemp,ak5,cp,lay);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  (*derivs)(x+a6*h,ytemp,ak6,cp,lay);
  for (i=1;i<=n;i++)
    yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
  for (i=1;i<=n;i++)
    yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
  free_vector(ytemp,1,n);
  free_vector(ak6,1,n);
  free_vector(ak5,1,n);
  free_vector(ak4,1,n);
  free_vector(ak3,1,n);
  free_vector(ak2,1,n);
}
/*************************************---*************************************/
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
/**************************************-**************************************/
/*                                                                           */
/*                                POP_RKQS_01                                */
/*                                                                           */
/*  Modified "rkqs" (NumRecInC) to call a particular "rkck" routine, and     */
/*  to have "nvar" = 1.                                                      */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void pop_rkqs_01(float *y,       // Vm at x (gets updated)
		 float dydx,     // deriv. of Vm at x
		 float *x,       // time (gets updated)
		 float htry,     // Step size to attempt
		 float eps,      // required accuracy
		 float yscal,    // scale factor for error
		 float *hdid,    // The step size actually accomplished
		 float *hnext,   // Estimated next step size
		 struct pop_cell *c)  // Pointer to cell
{
  float errmax,h,xnew,yerr,ytemp;

  /*
    printf("  POP_RKQS_01\n");
    printf("    y = %f\n",*y);
    printf("    dydx = %f\n",dydx);
    printf("    x = %f\n",*x);
    printf("    htry = %f\n",htry);
    printf("    eps = %f\n",eps);
    printf("    yscal = %f\n",yscal);
  */

  h=htry;
  for (;;){
    pop_rkck_01(*y,dydx,*x,h,&ytemp,&yerr,c);
    /*printf("  %f %f\n",ytemp,yerr);*/
    errmax = fabs(yerr/yscal);
    errmax /= eps;
    if (errmax > 1.0){
      h = SAFETY*h*pow(errmax,PSHRNK);
      if (h < 0.1*h){  // Only if 'h' is negative and integrating in reverse
	h *= 0.1;
	printf("POP_UTIL *** HERE h = %f\n",h);
	printf("THIS SHOULD NEVER HAPPEN\n");
	exit(0);
      }
      xnew = (*x)+h;
      if (xnew == *x)
	nrerror("stepsize underflow in pop_rkqs_01");
      continue;
    }else{
      if (errmax > ERRCON)
	*hnext = SAFETY*h*pow(errmax,PGROW);
      else
	*hnext=5.0*h;
      *x += (*hdid=h);
      *y=ytemp;
      /*printf("      ---> x,y = %f %f\n",*x,*y);*/
      break;
    }
  }
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
/*************************************---*************************************/
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
/**************************************-**************************************/
/*                                                                           */
/*                                POP_RKQS_02                                */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void pop_rkqs_02(float y[], float dydx[], int n, float *x, float htry,
		 float eps, float yscal[], float *hdid, float *hnext,
		 struct pop_cell *cp, struct pop_layer *lay,
		 void (*derivs)(float, float [], float [], struct pop_cell *,
				struct pop_layer *))
{
  void pop_rkck_02(float y[], float dydx[], int n, float x, float h,
		   float yout[],float yerr[], struct pop_cell *cp,
		   struct pop_layer *lay,
		   void (*derivs)(float, float [], float [], struct pop_cell *,
				  struct pop_layer *));
  int i;
  float errmax,h,xnew,*yerr,*ytemp,wyy;

  /*printf("  IFC_RKQS\n");
    printf("    y = %f\n",y[1]);
    printf("    dydx = %f\n",dydx[1]);
    printf("    x = %f\n",*x);
    printf("    htry = %f\n",htry);
    printf("    eps = %f\n",eps);
    printf("    yscal = %f\n",yscal[1]);*/

  yerr = vector(1,n);
  ytemp = vector(1,n);
  h=htry;
  for (;;){
    pop_rkck_02(y,dydx,n,*x,h,ytemp,yerr,cp,lay,derivs);
    /*printf("  %f %f\n",ytemp[1],yerr[1]);*/
    errmax=0.0;
    for (i=1;i<=n;i++){
      //errmax = FMAX(errmax,fabs(yerr[i]/yscal[i]));
      wyy = fabs(yerr[i]/yscal[i]);
      if (wyy > errmax)
	errmax = wyy;
    }
    errmax /= eps;
    if (errmax > 1.0){
      h = SAFETY*h*pow(errmax,PSHRNK);
      if (h < 0.1*h){  // Only if 'h' is negative and integrating in reverse
	h *= 0.1;
	printf("POP_UTIL *** HERE h = %f\n",h);
	printf("THIS SHOULD NEVER HAPPEN\n");
	exit(0);
      }
      xnew = (*x)+h;
      if (xnew == *x)
	nrerror("stepsize underflow in rkqs");
      continue;
    }else{
      if (errmax > ERRCON)
	*hnext = SAFETY*h*pow(errmax,PGROW);
      else
	*hnext=5.0*h;
      *x += (*hdid=h);
      for (i=1;i<=n;i++)
	y[i]=ytemp[i];
      /*printf("      ---> x,y = %f %f\n",*x,y[1]);*/
      break;
    }
  }
  free_vector(ytemp,1,n);
  free_vector(yerr,1,n);
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
/**************************************-**************************************/
/*                                                                           */
/*                            POP_PROCESS_SPIKE_01                           */
/*                                                                           */
/*****************************************************************************/
void pop_process_spike_01(c,x,y,yd)
     struct pop_cell *c;
     float x;             // Time of spike (in seconds)
     float y;             // Voltage at time "x"
     float yd;            // Vm Dendr  (2-comp)
{
  int j;
  int toffset,savti;
  struct pop_syn *tsyn;

  // Store spike time in msec
  if (c->savs){
    if (c->ns >= c->maxscnt)      // Assure sufficient storage
      pop_util_cell_spike_inc(c,2);
    c->s[c->ns] = x*1000.0;  // Store time in msec
    c->ns += 1;
  }

  // Store Vm, use 3 points for AP:  before, peak, after
  if (c->savv){
    if ((c->vn+2) >= c->vnmax) // Need 3 storage locations
      exit_error("POP_PROCESS_SPIKE_01","Vm storage exceeded");
    c->vmt[c->vn] = x; // Save current time
    c->vm[c->vn] = y;  // Save Vm

    //
    //  WYETH - for Poisson units, should skip the rest of this, YES?
    //  WYETH - for Poisson units, should skip the rest of this, YES?
    //  WYETH - for Poisson units, should skip the rest of this, YES?
    //  WYETH - for Poisson units, should skip the rest of this, YES?
    //
    if (c->savd){
      //c->vmt[c->vn] = x;   // Save current time
      c->vmd[c->vn] = yd;    // Save Vm  WYETH - FIX
      c->vmd[c->vn+1] = yd;  // Save Vm
      c->vmd[c->vn+2] = yd;  // Save Vm
    }
    c->vn += 1;
    c->vmt[c->vn] = x;
    c->vm[c->vn] = c->cifcp->v_spike;   // Save voltage for peak of spike
    c->vn += 1;
    c->vmt[c->vn] = x;
    c->vm[c->vn] = c->cifcp->v_reset_x; // Save voltage after spike
    c->vn += 1;
    c->xsav = x;
    c->reft = x + c->cifcp->trefr_x_s;  // Set time when refraction ends
                                        // Used somehow for NMDA processing?
  }
    
  // Adaptation after spike

  // WYETH HERE NEW ad1g
  //  'c->gta0' is sent for long-term storage; 'savti' is returned,
  //  but perhaps 'n0' should be set instead?  I.e., will the monitor
  //  window be watching 'n0' if there are no inputs?  We probably want
  //  to advance 'n0' slowly, only to the time of current spike???

  // WYETH Bug Fix, 2011:  gta_pulse was too long, writing past 'c->nad'

  if (c->nad > 0){
    pop_util_circbuf_mask_update_t(c->ad1g,x,c->gta_pulse,c->nad,c->gta0,
				   &savti);
    /***   // wygad
    if (c->gta0 != NULL){
      printf(">>>_._ savti %d  (T_spk %f, bufdt %f)\n",savti,x,c->ad1g->dt);
      printf("    c->nad = %d\n",c->nad);
    }
    ***/
  }

  // Put spike in inputs for all post-synaptic cells
  tsyn = c->out;
  while(tsyn != NULL){  // For each output synapse
    pop_util_process_spike_postsyn(tsyn,x,c->input_ex_sdf,c->input_ex_sdtau);
    tsyn = tsyn->pre_next; // Follow synapses for which "c" is pre-syn
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             POP_STEP_POISS_01                             */
/*                                                                           */
/*  NOTES                                                                    */
/*  - This advances between [x,x+h].                                         */
/*                                                                           */
/*****************************************************************************/
void pop_step_poiss_01(c,x1,x2)
     struct pop_cell *c;
     float x1,x2;         // Time limits of full simulation (s)
{
  int k;
  int n01,i1;
  float x,dt,t,g_ex,g_in,g_ad,tt,gexd,gind,prob,rate;
  struct poiss_param *p;
  struct pop_circbuf *gbuf,*i_gbuf;

  //printf("c->name = %s\n",c->name);
  //printf("x1,x2 = %f  %f    c->h = %f\n",x1,x2,c->h);

  x = c->x; // For convenience

  if (c->savv && (c->vn < c->vnmax) && ((x - c->xsav) > c->dxsav)){
    c->vmt[c->vn] = x;  // Store results no more often than "dxsav"
    c->vm[c->vn] = c->y;
    c->vn += 1;
    c->xsav = x;
  }

  //
  //  Get value for dynamic synaptic inputs at time index close to 'x'
  //  Note, we aren't bothering to interpolate, unlike in 'pop_util_deriv_01'.
  //
  gbuf   = c->ex1g;
  i_gbuf = c->in1g;
  if ((x > gbuf->t) || (x < (gbuf->t - (gbuf->n-1)*gbuf->dt))){
    printf("x=%f gbuf->t=%f  gbuf->i=%d\n",x,gbuf->t,gbuf->i);
    exit_error("POP_STEP_POISS_01","x not covered by gbuf");
  }
  // Compute indices before and after desired t-value
  tt = (gbuf->t - x)/gbuf->dt;
  k = (int)(0.5 + tt);  // Round, given that 'tt' is >= 0.0
  i1 = gbuf->i - k;
  if (i1 < 0)
    i1 += gbuf->n;
  gexd =   gbuf->d[i1];
  gind = i_gbuf->d[i1];

  //
  // Note, this is using interpolation, but it seemed easiery to let it be
  // this way for now.  We should change it for consistency w/ the lines above.
  //
  if (c->ad1g != NULL){
    g_ad = pop_util_circbuf_val_at_t_interp(c->ad1g,x);
  }else{
    g_ad = 0.0;
  }

  n01 = c->n0 - 1;  // n0 is length of pre-computed conductance inputs
  if (n01 <= 0){  // If the cell has no such input ...
    g_ex = 0.0;
    g_in = 0.0;
  }else{ // Interpolate to get g's
    t = x * c->samp0;   // Time in ms (float)
    k = (int)t;     // Time in ms (int)

    if (k < n01){ // Note, k does go up to 'c->n0'
      dt = t - (float)k;
      g_ex = c->gtx0[k] +  dt * (c->gtx0[k+1]-c->gtx0[k]);
      g_in = c->gti0[k] +  dt * (c->gti0[k+1]-c->gti0[k]);
    }else{
      g_ex = c->gtx0[n01];
      g_in = c->gti0[n01];
    }
  }
  g_ex += gexd;
  g_in += gind;

  //printf("___ g_ex TOT = %f\n",g_ex);
  //printf("___ g_in TOT = %f\n",g_in);

  p = c->cpoissp;

  //
  //  Compute firing probability using Poisson <spike_gen> params
  //
  //  This is "default" method
  //
  rate = p->offset1 + p->scale*(g_ex - p->scale_in*g_in - p->scale_ad*g_ad
				+ p->offset0);
  // The above 'rate' is interpreted as spikes/s
  // Now, conver it to prob. of a spike w/i the simulation time step, HMAX

  prob = rate * HMAX;
  if (prob > 1.0){
    prob = 1.0;
  }else if (prob < 0.0){
    prob = 0.0;
  }

  x += HMAX;  // WYETH - do entire step
  c->y = prob / HMAX;  // Store prob as a rate, instead of Vm
  c->nok += 1;  // Is this needed ?

  if (prob > 0.0){
    if (pop_util_ran2(&globseed) < prob){

      // **** WYETH DEBUG
      //if (prob > 0.500)
      //printf("  HERE Prob > 0.5,  prob = %f\n",prob);

      pop_process_spike_01(c,x,c->y,0.0); // Save spikes, Vm, and do AHP
      x += c->cifcp->trefr_x/1000.0;    // Advance x over refractory period
    }
  }

  c->x = x;     // Restore values to "c"
  c->h = HMAX;  // To be used as the next step size
}
/*************************************---*************************************/
#define TINY 1.0e-30
/**************************************-**************************************/
/*                                                                           */
/*                                 POP_STEP_01                               */
/*                                                                           */
/*  NOTES                                                                    */
/*  - This will only evaluate the deriv routine between [x,x+h], inclusive.  */
/*    Therefore, all values required by the deriv routine must be known in   */
/*    that interval.                                                         */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void pop_step_01(c,x1,x2)
     struct pop_cell *c;
     float x1,x2;         // Time limits (sec)
{
  float x,y,dydx,yscal,h,hdid,hnext;

  x = c->x; // For convenience
  y = c->y;
  h = c->h;

  pop_util_deriv_01(x,y,&dydx,c); // Fill dydx, x and y are unchanged
  yscal = fabs(y)+fabs(dydx*h)+TINY; /* Scaling used to monitor accuracy, this
					general purpose choice can be modified
					if need be */

  //
  //  WYETH, move these 6 lines above the previous 2 ???  (Feb 2012)
  //
  if (c->savv && (c->vn < c->vnmax) && ((x - c->xsav) > c->dxsav)){
    c->vmt[c->vn] = x;  // Store results no more often than "dxsav"
    c->vm[c->vn] = y;
    c->vn += 1;
    c->xsav = x;
  }

  if ((x+h) > x2) // If stepsize overshoots endpoint, decrease
    h = x2-x;
  pop_rkqs_01(&y,dydx,&x,h,c->eps,yscal,&hdid,&hnext,c); // Set new y

  // Force small steps, to avoid jumping past threshold or far into future
  if (hnext > HMAX)
    hnext = HMAX;
  
  if (hdid == h)
    c->nok += 1;
  else
    c->nbad += 1;

  if (y >= c->cifcp->v_th_x){

    //
    //  WYETH HERE - ADDED 16 May 2009 - to avoid spikes occurring after
    //    the end of the simulation.  COMMENTED OUT SAME DAY
    //
    // if (x < x2)

    pop_process_spike_01(c,x,y,0.0); // Save spikes, Vm, and do AHP


    y = c->cifcp->v_reset_x;          // Voltage becomes V_reset
    x += c->cifcp->trefr_x/1000.0;    // Advance x over refractory period
    hnext = c->h1;                // Start w/ default step
  }

  if (fabs(hnext) <= c->hmin){
    printf("hnext = %e\n",hnext);
    nrerror("Step size too small in odeint");
  }
  h = hnext;

  c->x = x; // Restore values to "c"
  c->y = y;
  c->h = h;
}
#undef TINY
/*************************************---*************************************/
#define TINY 1.0e-30
/**************************************-**************************************/
/*                                                                           */
/*                                 POP_STEP_02                               */
/*                                                                           */
/*  NOTES                                                                    */
/*  - This will only evaluate the deriv routine between [x,x+h], inclusive.  */
/*    Therefore, all values required by the deriv routine must be known in   */
/*    that interval.                                                         */
/*                                                                           */
/*****************************************************************************/
void pop_step_02(lay,c,x1,x2)
     struct pop_layer *lay; /* WYETH - is this needed here? ***/
     struct pop_cell *c;
     float x1,x2;  /* Time limits (sec) */
{
  int nvar;
  float x,yscal[3],h,hdid,hnext;
  float y[3];
  float dydx[3];

  nvar = 2;
  x = c->x; /* For convenience */
  y[1] = c->y;  /* Soma Vm */
  y[2] = c->yd; /* Dend Vm */
  h = c->h;

  pop_util_deriv_02(x,y,dydx,c,lay); /* Fill dydx, x and y are unchanged */

  /* Scaling to monitor accuracy, this general purpose choice can be
     changed if need be */
  yscal[1] = fabs(y[1])+fabs(dydx[1]*h)+TINY;
  yscal[2] = fabs(y[2])+fabs(dydx[2]*h)+TINY;

  if (c->savv && (c->vn < c->vnmax) && ((x - c->xsav) > c->dxsav)){
    c->vmt[c->vn] = x;  /* Store results no more often than "dxsav" */
    c->vm[c->vn] = y[1];
    if (c->savd)
      c->vmd[c->vn] = y[2];
    c->vn += 1;
    c->xsav = x;
  }

  if ((x+h) > x2) // If stepsize overshoots endpoint, decrease
    h = x2-x;

  // Set new y
  pop_rkqs_02(y,dydx,nvar,&x,h,c->eps,yscal,&hdid,&hnext,c,lay,
	      pop_util_deriv_02);

  // Force small steps, to avoid jumping past threshold or far into future
  if (hnext > HMAX)
    hnext = HMAX;
  
  if (hdid == h)
    c->nok += 1;
  else
    c->nbad += 1;

  
  // If a spike occurs
  if (y[1] >= c->cifcp->v_th_x){
    pop_process_spike_01(c,x,y[1],y[2]); // Save spikes, Vm, and do AHP
    
    y[1] = c->cifcp->v_reset_x;       // Voltage becomes V_reset

    /** NEW WAY OF REFR, IN 'derivs' **/
    /*x += c->cifcp->trefr_x/1000.0;*/    /* Advance x over refractory period */
    /*hnext = c->h1; */               /* Start w/ default step */
  }

  if (fabs(hnext) <= c->hmin){
    printf("hnext = %e\n",hnext);
    nrerror("Step size too small in odeint");
  }
  h = hnext;

  c->x = x; // Restore values to "c"
  c->y = y[1];
  c->yd = y[2];
  c->h = h;
}
#undef TINY
/**************************************-**************************************/
/*                                                                           */
/*                           POP_UTIL_INIT_CELL_TRIAL                        */
/*                                                                           */
/*   Initialize the state of the cell for the current trial.                 */
/*                                                                           */
/*   *** - Parts of this copied from:  'model_spop_config_cell'              */
/*         WYETH - should remove overlap                                     */
/*                                                                           */
/*****************************************************************************/
void pop_util_init_cell_trial(c,x1,x2)
     struct pop_cell *c;
     float x1,x2;
{
  int i;
  int nstor;

  c->y = c->ystart;
  if (c->cifcp->g_tran > 0.0)
    c->yd = c->ystart;

  c->x = x1;
  c->h = c->h1;

  c->reft = -1.0; // Make sure we do not start in refractory state

  // Copied from 'model_spop_config_cell'
  c->gex1t = 0.0; // Overall time for computing g, shared for e/i
  c->gex1  = 0.0;  // Set starting conductance and deriv to zero
  c->gdex1 = 0.0;
  c->gin1  = 0.0;  // Set starting conductance and deriv to zero
  c->gdin1 = 0.0;

  if (c->savv){
    c->xsav = x1 - 2.0*c->dxsav; // Be sure to save first value
    c->vn = 0;
    nstor = (x2-x1)/(c->dxsav/2.0); // Should dxsav be used?  HMAX?
    if ((c->vm != NULL) && (c->vnmax != nstor)){ // Old stor wrong len
      myfree(c->vm); c->vm = NULL;
      myfree(c->vmt); c->vmt = NULL;
      if (c->savd){
	myfree(c->vmd);
	c->vmd = NULL;
      }
    }
    if (c->vm == NULL){
      c->vnmax = nstor;
      c->vm  = (float *)myalloc(c->vnmax * sizeof(float));
      c->vmt = (float *)myalloc(c->vnmax * sizeof(float));
      if (c->savd)
	c->vmd  = (float *)myalloc(c->vnmax * sizeof(float));
    }else{
      if (c->vnmax != nstor){
      }
    }
  }

  if (c->sava){  // Save g_ad
    c->gta_n = 0;  // WYETH - It appears this is never used

    c->gta_max = my_rint((x2 - x1)*c->samp0);
    c->gta_max += c->nad + 1;  // Must accommodate duration + mask length

    if (c->gta0 == NULL){  // First time
      c->gta0 = (float *)myalloc(c->gta_max*sizeof(float));
    }

    // Clear this storage
    for(i=0;i<c->gta_max;i++)
      c->gta0[i] = 0.0;
  }


  if (c->savs){
    c->ns = 0;
    if (c->s == NULL) // First time through
      c->s = (float *)myalloc(c->maxscnt*sizeof(float));
  }

  if (c->savg){
    c->gtex1[0] = c->gex1;   // Store initial value
    c->gtex1t[0] = c->gex1t;
    c->gtin1[0] = c->gin1;   // Store initial value
    c->gtex1n = 1;           // Reset counter of values currently stored
  }
  


  // Buffers to hold upcoming input spikes, and recent g
  pop_cell_zero_circbuf(c->ex1s);  // WYETH - fixed bug, see the called func.
  pop_cell_zero_circbuf(c->in1s);

  // Hold recent dynamic conductance input
  pop_cell_zero_circbuf(c->ex1g);
  pop_cell_zero_circbuf(c->in1g);
  if (c->ad1g != NULL)
    pop_cell_zero_circbuf(c->ad1g);  // WYETH - NEW

  // Clear all output synapses
  pop_cell_reset_all_out_synapses(c);

}
/**************************************-**************************************/
/*                                                                           */
/*                          POP_UTIL_RUN_LAYER_LIST                          */
/*                                                                           */
/*  Run simulation for a list of layers.                                     */
/*                                                                           */
/*****************************************************************************/
void pop_util_run_layer_list(myid,mylogf,m,r,s,laylist,n,x1,x2,seed,flag)
     int myid;
     char *mylogf;
     struct model_struct *m;        // Model parameter pair list
     struct response_struct *r;     // Response data
     struct stim_struct *s;         // Stimulus data and param pair list
     struct pop_layer **laylist;
     int n;
     float x1,x2;                   // from mpt->x1,x2  start, stop time
     int seed;
     int flag;                      // -1=run on the fly
{
  int i,j,k,l,mm,ti;
  int nstep,mon_flag,pflag;
  float dt,dtmon,tmon;
  double t;
  char tstr[SLEN];
  struct pop_cell *c;
  struct pop_layer *lay;

  pflag = 0;

  mylog(mylogf,"  POP_UTIL_RUN_LAYER_LIST\n");

  // Set seed for adding noise spikes to ex and in synaptic inputs.
  // And, now also for Poisson spike generation for Poisson units (Feb 2012).
  if (seed > 0)
    globseed = -seed;
  else
    globseed = seed;

  //printf("_______globseed = %d\n",globseed);

  // Initialize simulation state and storage for all cells, all layers
  for(l=0;l<n;l++){ // For each layer
    lay = laylist[l];
    if (lay->runflag == 1){
      for(i=0;i<lay->xn;i++){
	for(j=0;j<lay->yn;j++){
	  for(k=0;k<lay->zn;k++){ // For each cell
	    c = &(lay->c[i][j][k]);
	    //printf(" l = %d   i,j,k  %d %d %d\n",l,i,j,k);
	    pop_util_init_cell_trial(c,x1,x2);
	    //printf("     done\n");
	  }
	}
      }
    }
  }

  if (pflag) printf("  pop_util_run_layer_list   001\n");

  // Graphical monitor
  mon_flag = paramfile_get_int_param_default(r->ppl,"mon_flag",0);
  dtmon = paramfile_get_float_param_default(r->ppl,"mon_dt",0.010);
  tmon  = 0.0;    // Time since last update

  if (pflag) printf("  pop_util_run_layer_list   002\n");

  dt = 0.001; // Advance all by this step (s)
  nstep = (x2-x1)/dt;
  t = 0.0;

  sprintf(tstr,"    Run for %d time steps, stepsize %.5f sec\n",nstep,dt);
  mylog(mylogf,tstr);
  for (ti=0;ti<nstep;ti++){ // For each time step
    t += (double)dt;
    tmon += (double)dt;

    //
    //  Advance the synaptic inputs to each cell by one time step
    //
    for(l=0;l<n;l++){ // For each layer
      lay = laylist[l];
      if (lay->runflag == 1){
	for(i=0;i<lay->xn;i++){
	  for(j=0;j<lay->yn;j++){
	    for(k=0;k<lay->zn;k++){ // For each cell
	      c = &(lay->c[i][j][k]);
	      if (c->cifcp->g_tran == 0.0){
		if (flag == -1){  // Run on the fly (no pre-compute)
		  pop_input_01a(c,t+HMAX);   // Advance inputs 1 step
		}else{
		  pop_input_01(c,t+HMAX);    // Advance inputs 1 step
		}
	      }else
		pop_input_02(c,t+HMAX,lay); // Advance inputs w/ NMDA 1 step
	    }
	  }
	}
      }
    }

    if (pflag) printf("  pop_util_run_layer_list   003\n");

    //
    //  Advance Vm and output spikes of each cell by one time step
    //
    for(l=0;l<n;l++){ // For each layer
      lay = laylist[l];
      if (lay->runflag == 1){
	for(i=0;i<lay->xn;i++){
	  for(j=0;j<lay->yn;j++){
	    for(k=0;k<lay->zn;k++){ // For each cell
	      c = &(lay->c[i][j][k]);
	      mm = 0;
	      while(c->x < t){ // Advance this cell until its time reaches t

		// Advance the state by one step
		if (c->cifcp->g_tran == 0.0){
		  if (c->cpoissp != NULL){
		    pop_step_poiss_01(c,x1,x2);
		    //printf("c->x = %f\n",c->x);
		  }else
		    pop_step_01(c,x1,x2); // Advance the state by one step
		}else{
		  //if (c->y < -100.0)
		  if (mm > 195){
		    // WYETH - gta is too high
		    global_dump_flag = 1;
		    append_farray_xy_plot("z.dmp",c->vmt,c->vm,c->vn,"v");
		    append_farray_xy_plot("z.dmp",c->vmt,c->vmd,c->vn,"vd");
		    append_farray_plot("z.dmp.pl","gtx0",c->gtx0,c->n0,1);
		    append_farray_plot("z.dmp.pl","gtn0",c->gtn0,c->n0,1);
		    append_farray_plot("z.dmp.pl","gti0",c->gti0,c->n0,1);
		    if (c->gta0 != NULL){
		      append_farray_plot("z.dmp.pl","gta0",c->gta0,c->n0,1);
		    }
		    printf("WDBHere %s  x,y  %d %d\n",c->name,c->layx,c->layy);
		    exit(0);
		  }
		  pop_step_02(lay,c,x1,x2); // Advance state for 2-comp
		}
		mm += 1;
		if (mm > 200){ // WYETH was 50, until feedback model
		  printf("l,i,j,k = %d %d %d %d   mm = %d\n",l,i,j,k,mm);
		  mylog_exit(mylogf,
		    "POP_UTIL_RUN_LAYER_LIST  Too many steps in interval\n");
		}
	      }
	    }
	  }
	}
      }
    }
    //mylog(mylogf,"  time = %f\n",t);

    if (pflag) printf("  pop_util_run_layer_list   004  ti = %d\n",ti);


    // Graphical response monitor
    if (tmon > dtmon){
      /*** WYETH DEBUG ***/
      /*
	sprintf(tstr,"  time = %f\n",t);
	mylog(mylogf,tstr);*/

      if ((myid == -1) && (mon_flag)){
	// WYETH 2019 commented to avoid OpenGL/X11 dependence
	//mod_mon_monitor(m,r,s,laylist,n,(float)t);
      }
      tmon = 0.0;
    }

    if (pflag) printf("  pop_util_run_layer_list   005\n");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         POPU_MAR_WRITE_ONODE_ITEMS                        */
/*                                                                           */
/*  Given the onode 'o', write all item names and values as strings.         */
/*                                                                           */
/*   (int)   bcode                                                           */
/*   (int)   number of items                                                 */
/*   For each item                                                           */
/*     (str)   item name                                                     */
/*     (str)   item value                                                    */
/*   (int)  -bcode                                                           */
/*                                                                           */
/*****************************************************************************/
void popu_mar_write_onode_items(fout,o,bcode,vers)
     FILE *fout;
     struct onode *o;   // Onode
     int bcode;         // block code
     int vers;          // version of mar file format
{
  int nitem;
  struct onode *t,*s;

  fwrite((char *)&bcode,sizeof(int),1,fout);  // Start of block

  nitem = onode_count_otype(o,"item");        // Count number of items
  fwrite((char *)&nitem,sizeof(int),1,fout);

  s = NULL;
  t = onode_get_next_type(o->o,"item");
  while(t != NULL){

    carray_nchar_write(t->name,fout);        // item name
    carray_nchar_write(t->val,fout);         // item value string

    t = onode_get_next_type(t->next,"item");
  }

  bcode *= -1;
  fwrite((char *)&bcode,sizeof(int),1,fout);  // End of block
}
/**************************************-**************************************/
/*                                                                           */
/*                           POPU_MAR_WRITE_AREA_NAME                        */
/*                                                                           */
/*****************************************************************************/
void popu_mar_write_area_name(fout,mpt,lay,vers)
     FILE *fout;
     struct pop_top *mpt;
     struct pop_layer *lay;
     int vers;
{
  int bcode;

  bcode = 10010;  // Area name
  fwrite((char *)&bcode,sizeof(int),1,fout);  // Start of population

  carray_nchar_write(lay->area->name,fout);   // area name

  bcode *= -1;
  fwrite((char *)&bcode,sizeof(int),1,fout);  // End of population
}
/**************************************-**************************************/
/*                                                                           */
/*                             POPU_MAR_WRITE_AREA                           */
/*                                                                           */
/*****************************************************************************/
void popu_mar_write_area(fout,mpt,area,vers)
     FILE *fout;
     struct pop_top *mpt;
     struct pop_area *area;
     int vers;
{
  int bcode;

  bcode = 10001;
  fwrite((char *)&bcode,sizeof(int),1,fout);  // Start of population

  printf("      Area %s  (%d,%d)\n",area->name,area->xn,area->yn);

  carray_nchar_write(area->name,fout);   // area name
  fwrite((char *)&(area->x0),sizeof(float),1,fout);
  fwrite((char *)&(area->y0),sizeof(float),1,fout);
  fwrite((char *)&(area->xf),sizeof(float),1,fout);
  fwrite((char *)&(area->yf),sizeof(float),1,fout);
  fwrite((char *)&(area->xn),sizeof(int),1,fout);
  fwrite((char *)&(area->yn),sizeof(int),1,fout);
  fwrite((char *)&(area->umx),sizeof(float),1,fout);
  fwrite((char *)&(area->umy),sizeof(float),1,fout);
  fwrite((char *)&(area->sscale),sizeof(float),1,fout);

  bcode *= -1;
  fwrite((char *)&bcode,sizeof(int),1,fout);  // End of population
}
/**************************************-**************************************/
/*                                                                           */
/*                            POPU_MAR_WRITE_ATTRIB                          */
/*                                                                           */
/*****************************************************************************/
void popu_mar_write_attrib(fout,lay,vers)
     FILE *fout;
     struct pop_layer *lay;
     int vers;
{
  int i,j,k;
  int ai,*vatt,attn;
  float fval;
  struct pop_cell *c;

  //
  //  Write number of attributes, followed by names for each
  //
  attn = lay->attrib_fn;                    // Number of attributes
  fwrite((char *)&attn,sizeof(int),1,fout);
  for(i=0;i<attn;i++){
    carray_nchar_write(lay->attrib_fname[i],fout);     // attrib name
  }

  //
  //  Write  variation type of each attribute
  //
  vatt = popc_attrib_layer_variation(lay);  // Type of variation for each
  for(i=0;i<attn;i++){
    fwrite((char *)&(vatt[i]),sizeof(int),1,fout);
  }

  /*
  printf("Layer %s\n",lay->name);
  for(i=0;i<attn;i++){
    printf("  %d %3d %s\n",i,vatt[i],lay->attrib_fname[i]);
  }
  exit(0);
  */

  for(ai=0;ai<attn;ai++){
    if (vatt[ai] == 0){
      //
      //  0 - write one float value, no variation
      //
      fval = lay->c[0][0][0].attrib_f[ai];
      fwrite((char *)&fval,sizeof(float),1,fout);
    }else if (vatt[ai] == 1){
      //
      //  1 - write 2D array of values (xn by yn)
      //
      for(i=0;i<lay->xn;i++){
	for(j=0;j<lay->yn;j++){
	  c = &(lay->c[i][j][0]);
	  fwrite((char *)&(c->attrib_f[ai]),sizeof(float),1,fout);
	}
      }
    }else if (vatt[ai] == 2){
      //
      //  2 - write 1D array (zn) of values.
      //
      for(i=0;i<lay->zn;i++){
	c = &(lay->c[0][0][i]);
	fwrite((char *)&(c->attrib_f[ai]),sizeof(float),1,fout);
      }
    }else if (vatt[ai] == 3){
      //
      //  3 - write 3D array (xn,yn,zn) of values.
      //
      for(i=0;i<lay->xn;i++){
	for(j=0;j<lay->yn;j++){
	  for(k=0;k<lay->zn;k++){
	    c = &(lay->c[i][j][k]);
	    fwrite((char *)&(c->attrib_f[ai]),sizeof(float),1,fout);
	  }
	}
      }
    }else{
      exit_error("POPU_MAR_WRITE_ATTRIB","Unknown variation code");
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             POPU_MAR_WRITE_LAY                            */
/*                                                                           */
/*****************************************************************************/
void popu_mar_write_lay(fout,mpt,lay,vers)
     FILE *fout;
     struct pop_top *mpt;
     struct pop_layer *lay;
     int vers;
{
  int i;
  int bcode;
  int pflag = 0;

  bcode = 20001;
  fwrite((char *)&bcode,sizeof(int),1,fout);  // Start of population

  if (pflag == 1) printf("HERE __ 1\n");


  carray_nchar_write(lay->name,fout);        // pop name

  if (pflag == 1) printf("HERE __ 1.a\n");

  if (vers > 1)
    carray_nchar_write(lay->laytype,fout);   // pop type, "lgn", "default"...
  if (pflag == 1) printf("HERE __ 1.b\n");
  if (vers > 3)
    carray_nchar_write(lay->geomt,fout);     // "irregular", "default"...

  if (pflag == 1) printf("HERE __ 1.c\n");


  fwrite((char *)&(lay->xn),sizeof(int),1,fout);
  fwrite((char *)&(lay->yn),sizeof(int),1,fout);
  fwrite((char *)&(lay->zn),sizeof(int),1,fout);
  if (pflag == 1) printf("HERE __ 1.b\n");
  fwrite((char *)&(lay->x0),sizeof(float),1,fout);
  fwrite((char *)&(lay->y0),sizeof(float),1,fout);
  fwrite((char *)&(lay->xf),sizeof(float),1,fout);
  fwrite((char *)&(lay->yf),sizeof(float),1,fout);
  fwrite((char *)&(lay->oddxoff),sizeof(float),1,fout);

  if (pflag == 1) printf("HERE __ 2\n");

  if (lay->geomt != NULL){
    if (strcmp(lay->geomt,"irregular")==0){
      //
      //  Write config for irregular geometry
      //
      fwrite((char *)&(lay->irr_n),sizeof(int),1,fout);
      for(i=0;i<lay->irr_n;i++){
	fwrite((char *)&(lay->irr_x[i]),sizeof(float),1,fout);
	fwrite((char *)&(lay->irr_y[i]),sizeof(float),1,fout);
	fwrite((char *)&(lay->irr_id[i]),sizeof(int),1,fout);
      }
    }
  }

  if (pflag == 1) printf("HERE __ 3\n");

  fwrite((char *)&(lay->ninlist),sizeof(int),1,fout);  // Number of inputs

  if (vers > 2){
    popu_mar_write_attrib(fout,lay,vers);
  }

  if (lay->icon != NULL){
    popu_mar_write_onode_items(fout,lay->icon->o,20100,vers);
  }

  if (lay->area != NULL){
    popu_mar_write_area_name(fout,mpt,lay,vers);
  }

  bcode *= -1;
  fwrite((char *)&bcode,sizeof(int),1,fout);  // End of population
}
/**************************************-**************************************/
/*                                                                           */
/*                           POPU_MAR_WRITE_SYN_REG                          */
/*                                                                           */
/*****************************************************************************/
int popu_mar_write_syn_reg(fout,vers,laypre,laypost,c,nbyte_pre,iin)
     FILE *fout;
     int vers;
     struct pop_layer *laypre;
     struct pop_layer *laypost;
     struct pop_cell *c;
     int nbyte_pre;
     int iin;
{
  int cnt,ci_pre,xs,ys;
  short si;
  float w;
  struct pop_syn *syn;
  struct pop_cell *cpre;

  xs = laypre->yn * laypre->zn; // Multipliers for pre-syn index computation
  ys = laypre->zn;

  // Count the number of connections for this input
  cnt = 0;
  syn = c->in;
  while(syn != NULL){
    if (syn->inindx == iin)
      cnt += 1;
    syn = syn->post_next;
  }

  // Write the number of connections
  fwrite((char *)&cnt,sizeof(int),1,fout);

  // Write the cell index for each presyn cell
  syn = c->in;
  while(syn != NULL){
    if (syn->inindx == iin){

      cpre = syn->pre;  // Presyn cell
      ci_pre = xs*cpre->layx + ys*cpre->layy + cpre->layz; // Cell index

      //printf("xs = %d   ys=%d  cpre = %d %d\n",xs,ys,cpre->layx,cpre->layy);
      //printf("ci_pre = %d\n",ci_pre);

      // Write cell index
      if (nbyte_pre == 2){
	si = (short)ci_pre;
	fwrite((char *)&si,sizeof(short),1,fout);
      }else{  // assume 4
	fwrite((char *)&ci_pre,sizeof(int),1,fout);
      }

      // Write weight
      fwrite((char *)&syn->w,sizeof(float),1,fout);

    }
    syn = syn->post_next;
  }

  return cnt;
}
/**************************************-**************************************/
/*                                                                           */
/*                        POPU_MAR_WRITE_SYN_REG_PAIR                        */
/*                                                                           */
/*****************************************************************************/
void popu_mar_write_syn_reg_pair(fout,vers,laypre,laypre2,laypost,c,
				 nbyte_pre,nbyte_pre2,iin)
     FILE *fout;
     int vers;
     struct pop_layer *laypre,*laypre2;  // original and paired layers
     struct pop_layer *laypost;
     struct pop_cell *c;
     int nbyte_pre,nbyte_pre2;
     int iin;
{
  int cnt,ci_pre,xs,ys,xs2,ys2;
  short si;
  float w;
  struct pop_syn *syn,*syn1;
  struct pop_cell *cpre;

  xs = laypre->yn * laypre->zn; // Multipliers for pre-syn index computation
  ys = laypre->zn;

  xs2 = laypre2->yn * laypre2->zn; // Multipliers for pre-syn index computation
  ys2 = laypre2->zn;

  // Count the number of connections for this input
  cnt = 0;
  syn = c->in;
  while(syn != NULL){

    if (syn->inindx == iin){
      if (syn->si == NULL)
	exit_error("POPU_MAR_WRITE_SYN_REG_PAIR","No SI for synapse");
      else if (syn->si->ci == 0)
	cnt += 1;
    }
    syn = syn->post_next;
  }

  // Write the number of connections
  fwrite((char *)&cnt,sizeof(int),1,fout);

  // Write the cell index for each presyn cell
  syn = c->in;
  while(syn != NULL){
    if ((syn->inindx == iin) && (syn->si->ci == 0)){

      cpre = syn->pre;  // Presyn cell
      ci_pre = xs*cpre->layx + ys*cpre->layy + cpre->layz; // Cell index
      
      // Write cell index
      if (nbyte_pre == 2){
	si = (short)ci_pre;
	fwrite((char *)&si,sizeof(short),1,fout);
      }else{  // assume 4
	fwrite((char *)&ci_pre,sizeof(int),1,fout);
      }


      //  Find the paired synapse 'syn1' that has cell index 1
      syn1 = pop_cell_syn_si_find(NULL,c->in,syn->si->siu,1);

      cpre = syn1->pre;  // Presyn cell
      ci_pre = xs2*cpre->layx + ys2*cpre->layy + cpre->layz; // Cell index
      
      // Write cell index
      if (nbyte_pre2 == 2){
	si = (short)ci_pre;
	fwrite((char *)&si,sizeof(short),1,fout);
      }else{  // assume 4
	fwrite((char *)&ci_pre,sizeof(int),1,fout);
      }

      
      // Write weight
      fwrite((char *)&syn->w,sizeof(float),1,fout);

    }
    syn = syn->post_next;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            POPU_MAR_WRITE_SYN_LGN                         */
/*                                                                           */
/*****************************************************************************/
void popu_mar_write_syn_lgn(fout,vers,laypre,laypost,c,nbyte_pre)
     FILE *fout;
     int vers;
     struct pop_layer *laypre;
     struct pop_layer *laypost;
     struct pop_cell *c;
     int nbyte_pre;
{
  int i;
  int lxi,lyi,lzi,n0,n1,n,ci_pre,xs,ys;
  short si;
  float w;
  struct pop_lgn_in *lin;

  // Get pointer to appropriate LGN input
  lin = c->lgnin[0];

  n0 = lin->cn0;  // Number of OFF
  n1 = lin->cn1;  // Number of ON
  n  = n0 + n1;   // Total
  w = 1.0;        // All weights are 1.0 ***********

  xs = laypre->yn * laypre->zn; // Multipliers for pre-syn index computation
  ys = laypre->zn;


  // Write the number of connections
  fwrite((char *)&n,sizeof(int),1,fout);

  lzi = 0; // OFF
  for(i=0;i<n0;i++){

    lxi = lin->cx0[i];
    lyi = lin->cy0[i];
    ci_pre = xs*lxi + ys*lyi + lzi; // Cell index
    
    // Write cell index
    if (nbyte_pre == 2){
      si = (short)ci_pre;
      fwrite((char *)&si,sizeof(short),1,fout);
    }else{  // assume 4
      fwrite((char *)&ci_pre,sizeof(int),1,fout);
    }
    
    // Write weight
    fwrite((char *)&w,sizeof(float),1,fout);
  }


  lzi = 1; // ON
  for(i=0;i<n1;i++){

    lxi = lin->cx1[i];
    lyi = lin->cy1[i];
    ci_pre = xs*lxi + ys*lyi + lzi; // Cell index
    
    // Write cell index
    if (nbyte_pre == 2){
      si = (short)ci_pre;
      fwrite((char *)&si,sizeof(short),1,fout);
    }else{  // assume 4
      fwrite((char *)&ci_pre,sizeof(int),1,fout);
    }
    
    // Write weight
    fwrite((char *)&w,sizeof(float),1,fout);
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                             POPU_MAR_WRITE_INPUT                          */
/*                                                                           */
/*  (int)  30001                                                             */
/*  (str)  Pre-syn layer name                                                */
/*  (str)  Post-syn layer name                                               */
/*  (int)  nbyte pre (typically 2 or 4)                                      */
/*  (int)  nbyte post (typically 2 or 4)                                     */
/*  For each post-syn cell [Use post-syn layer name to determine count]      */
/*    (s/i)  Cell index Post                                                 */
/*    (int)  n synapses                                                      */
/*    For each synapse                                                       */
/*      (s/i)  Cell index Pre [Use pre-syn layer name to decode]             */
/*      (flt)  weight                                                        */
/*  (int)  -30001                                                            */
/*                                                                           */
/*****************************************************************************/
void popu_mar_write_input(fout,mpt,lay,iin,vers)
     FILE *fout;
     struct pop_top *mpt;
     struct pop_layer *lay;
     int iin;       // Input index
     int vers;
{
  int i,j,k;
  int xn,yn,zn,bcode,ci,nbyte_pre,nbyte_pre2,nbyte_post,lgnflag,sycnt;
  short si;
  struct pop_cell *c;
  struct pop_layer *lpre,*lpre2;
  struct onode *oi;

  printf("    POPU_MAR_WRITE_INPUT\n");

  oi = lay->inlist[iin];  // <input> onode

  // Determine pre-syn layer
  if (popu_input_is_paired(lay,iin)){
    bcode = 30002;

    lpre = popu_input_get_origin(mpt,lay,iin);
    printf("      (pair) Pre-syn layer:  %s\n",lpre->name);

    lpre2 = popu_input_get_paired_lay(mpt,lay,iin);
    printf("      (pair) Pre-syn layer 2:  %s\n",lpre2->name);

    //exit_error("POPU_MAR_WRITE_INPUT","Paired input not imp'd yet");

  }else if (onode_test_chr(oi,"type","bg")){
    //
    //  WYETH - write an empty block and return
    //
    bcode = 30003;
    fwrite((char *)&bcode,sizeof(int),1,fout);  // Start of input
    bcode *= -1;
    fwrite((char *)&bcode,sizeof(int),1,fout);  // End of input
    printf("      bg input\n");
    return;
  }else if (onode_test_chr(oi,"type","template")){ // WYETH Feb 2014
    //
    //  WYETH - template
    //
    bcode = 30004;
    lpre = popu_input_get_origin(mpt,lay,iin);
    lpre2 = NULL;
    printf("      Pre-syn layer:  %s\n",lpre->name);
  }else{
    bcode = 30001;
    lpre = popu_input_get_origin(mpt,lay,iin);
    lpre2 = NULL;
    printf("      Pre-syn layer:  %s\n",lpre->name);
  }


  fwrite((char *)&bcode,sizeof(int),1,fout);  // Start of input

  if (pop_util_lclass_is_lgn0(lpre))  // Determine whether an 'lgn0' layer
    lgnflag = 1;
  else
    lgnflag = 0;

  //  Write layer names
  carray_nchar_write(lpre->name,fout);
  if (lpre2 != NULL)
    carray_nchar_write(lpre2->name,fout);
  carray_nchar_write(lay->name,fout);

  xn = lay->xn;  // Post-syn dimensions
  yn = lay->yn;
  zn = lay->zn;

  if (xn*yn*zn <= 65536)
    nbyte_post = 2;     // Post-syn cell indices are SHORT (2-byte)
  else
    nbyte_post = 4;     // Post-syn cell indices are INT (4-byte)

  if (lpre->xn*lpre->yn*lpre->zn <= 65536)
    nbyte_pre = 2;
  else
    nbyte_pre = 4;

  if (lpre2 != NULL){
    if (lpre2->xn*lpre2->yn*lpre2->zn <= 65536)
      nbyte_pre2 = 2;
    else
      nbyte_pre2 = 4;
  }

  //  Write byte-length flags
  fwrite((char *)&nbyte_pre,sizeof(int),1,fout);
  if (lpre2 != NULL)
    fwrite((char *)&nbyte_pre2,sizeof(int),1,fout);
  fwrite((char *)&nbyte_post,sizeof(int),1,fout);


  sycnt = 0;
  ci = 0;  // Cell index, post-syn
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){    //  For each cell

	c = &(lay->c[i][j][k]);  // Post-syn cell

	// Write the index of this cell
	if (nbyte_post == 2){
	  si = (short)ci;
	  fwrite((char *)&si,sizeof(short),1,fout);
	}else{ // assume 4
	  fwrite((char *)&ci,sizeof(int),1,fout);
	}

	if (lgnflag){
	  //  LGN
	  popu_mar_write_syn_lgn(fout,vers,lpre,lay,c,nbyte_pre);
	}else if (lpre2 != NULL){
	  popu_mar_write_syn_reg_pair(fout,vers,lpre,lpre2,lay,c,
				      nbyte_pre,nbyte_pre2,iin);
	}else{
	  //  Regular
	  sycnt += popu_mar_write_syn_reg(fout,vers,lpre,lay,c,nbyte_pre,iin);
	}
	ci += 1;
      }
    }
  }

  printf("      %d synapses\n",sycnt);

  bcode *= -1;
  fwrite((char *)&bcode,sizeof(int),1,fout);  // End of input
}
/**************************************-**************************************/
/*                                                                           */
/*                               POPU_MAR_WRITE                              */
/*                                                                           */
/*  Write a file that contains a description of the constructed model.       */
/*                                                                           */
/*  Codes                                                                    */
/*   3000 - General                                                          */
/*  10000 - Area                                                             */
/*  20000 - Population                                                       */
/*  20100 -   Icon                                                           */
/*  30000 - Input                                                            */
/*  40000 - Icon                                                             */
/*                                                                           */
/*****************************************************************************/
void popu_mar_write(mpt,outfile,vers)
     struct pop_top *mpt;
     char *outfile;
     int vers;
{
  FILE *fopen(),*fout;
  int i,j,k;
  int fcode,n;
  struct pop_layer *lay;
  int pflag = 0;

  printf("  POPU_MAR_WRITE\n");

  fcode = 15551777;    // Endian code
  if ((vers <= 0) || (vers > 4))
    mylogx(mpt->logf,"POPU_MAR_WRITE","Unknown version number");

  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** File %s\n",outfile);
    mylogx(mpt->logf,"POPU_MAR_WRITE","Cannot open file");
  }

  fwrite((char *)&fcode,sizeof(int),1,fout);  // 15551777
  fwrite((char *)&vers,sizeof(int),1,fout);   // Version

  if (pflag == 1) printf("HERE 001\n");

  //  General
  k = 3000;
  fwrite((char *)&k,sizeof(int),1,fout);   // Start of general block
  fwrite((char *)&(mpt->xn),sizeof(int),1,fout);
  fwrite((char *)&(mpt->yn),sizeof(int),1,fout);
  fwrite((char *)&(mpt->sscale),sizeof(float),1,fout);
  fwrite((char *)&(mpt->tscale),sizeof(float),1,fout);
  k = -k;
  fwrite((char *)&k,sizeof(int),1,fout);   // End of general block

  if (pflag == 1) printf("HERE 002\n");

  //  Areas
  printf("    Writing Areas\n");
  k = 10000;
  n = mpt->narea;
  fwrite((char *)&k,sizeof(int),1,fout);   // Start of pop block
  fwrite((char *)&n,sizeof(int),1,fout);   // Number of pop's
  for(i=0;i<n;i++){
    popu_mar_write_area(fout,mpt,mpt->area[i],vers);
  }
  k = -k;
  fwrite((char *)&k,sizeof(int),1,fout);   // End of pop block

  if (pflag == 1) printf("HERE 003\n");

  //  Populations
  k = 20000;
  n = mpt->nlay;
  fwrite((char *)&k,sizeof(int),1,fout);   // Start of pop block
  fwrite((char *)&n,sizeof(int),1,fout);   // Number of pop's
  for(i=0;i<n;i++){
    popu_mar_write_lay(fout,mpt,mpt->lay[i],vers);
  }
  k = -k;
  fwrite((char *)&k,sizeof(int),1,fout);   // End of pop block

  if (pflag == 1) printf("HERE 004\n");

  //  Inputs
  k = 30000;
  fwrite((char *)&k,sizeof(int),1,fout);   // Start of input block
  for(i=0;i<mpt->nlay;i++){
    lay = mpt->lay[i];
    for(j=0;j<lay->ninlist;j++){
      popu_mar_write_input(fout,mpt,lay,j,vers);
    }
  }
  k = -k;
  fwrite((char *)&k,sizeof(int),1,fout);   // End of input block

  if (pflag == 1) printf("HERE 005\n");

  fclose(fout);

  printf("    Wrote %s\n",outfile);
}
/**************************************-**************************************/
/*                                                                           */
/*                            POPU_RSP_GEN_GET_MAP                           */
/*                                                                           */
/*  Get a pointer to the save map for the layer 'lay'.                       */
/*                                                                           */
/*****************************************************************************/
int ***popu_rsp_gen_get_map(mpt,smap,lay)
     struct pop_top *mpt;
     int ****smap;
     struct pop_layer *lay;
{
  int k;

  k = popu_get_layer_index_by_name(mpt,lay->name);
  if (k == -1){
    printf("  Layer  %s\n",lay->name);
    exit_error("POPU_RSP_GEN_GET_MAP","Layer not found");
  }

  if (smap[k] == NULL){  // Has not been created yet
    smap[k] = get_zero_3d_iarray(lay->xn,lay->yn,lay->zn);
  }

  return smap[k];
}
/**************************************-**************************************/
/*                                                                           */
/*                          POPU_RSP_GEN_GET_LAY_PTR                         */
/*                                                                           */
/*****************************************************************************/
struct pop_layer *popu_rsp_gen_get_lay_ptr(mpt,lname)
     struct pop_top *mpt;
     char *lname;
{
  struct pop_layer *pl;

  pl = pop_cell_get_layer_pointer(mpt->lay,mpt->nlay,lname);
  if (pl == NULL){
    printf("  name:  %s\n",lname);
    exit_error("POPU_RSP_GEN_GET_LAY_PTR","Population name not found");
  }

  return pl;
}
/**************************************-**************************************/
/*                                                                           */
/*                          POPU_RSP_GEN_CMAP_FROM_W                         */
/*                                                                           */
/*****************************************************************************/
int ***popu_rsp_gen_cmap_from_w(w,xn,yn,zn,optype,fcrit,rn)
     float ***w;
     int xn,yn,zn;
     char *optype;
     float fcrit;
     int *rn;
{
  int i,j,k;
  int ***cmap,n;

  cmap = get_zero_3d_iarray(xn,yn,zn);

  n = 0;
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	if (strcmp(optype,"w_agg_min")==0){
	  //if (w[i][j][k] >= fcrit){
	  if (w[i][j][k] > fcrit){  // WYETH - 0 criterion for weights
	    cmap[i][j][k] = 1;
	    n += 1;
	  }
	}else{
	  printf("  *** optype:  %s\n",optype);
	  exit_error("POPU_RSP_GEN_CMAP_FROM_W","Unknown optype");
	}
      }
    }
  }

  *rn = n;

  return cmap;
}
/**************************************-**************************************/
/*                                                                           */
/*                          POPU_RSP_GEN_APPEND_BY_MAP                       */
/*                                                                           */
/*****************************************************************************/
void popu_rsp_gen_append_by_map(cmap,xn,yn,zn,outfile,name_save,dtype,samp_str,
				name_pre,datid,sline,nnew,w)
     int ***cmap;
     int xn,yn,zn;
     char *outfile,*name_save,*dtype,*samp_str,*name_pre,*datid;
     char *sline;
     int nnew;
     float ***w;  // [xn][yn][zn] May be null
{
  int i,j,k;
  float wsum;
  char t1[SLEN],t2[SLEN3],t3[512];  // 't3' added for compiler warning (2019)

  sprintf(t1,"##  %s\n##\n",sline);
  append_string_to_file(outfile,t1);
  sprintf(t1,"##  %d units marked\n##\n",nnew);
  append_string_to_file(outfile,t1);

  if (w != NULL)
    wsum = (float)sum_3d_farray(w,0,xn,0,yn,0,zn);

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	if (cmap[i][j][k] >= 1){
	  sprintf(t1,"save_pop_unit_as_%s_%d_%d_%d %s %s",name_save,i,j,k,
		  dtype,samp_str);
	  sprintf(t2,"%s %s %d %d %d %s",t1,name_pre,i,j,k,datid);
	  if (w != NULL)
	    sprintf(t3,"%s   # w %f (%f)\n",t2,w[i][j][k],w[i][j][k]/wsum);
	  else
	    sprintf(t3,"%s\n",t2);
	  append_string_to_file(outfile,t3);
	  cmap[i][j][k] = -1;  // Indicates this has been written
	}
      }
    }
  }
  append_string_to_file(outfile,"\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                             POPU_RSP_GEN_EXPAND                           */
/*                                                                           */
/*****************************************************************************/
int popu_rsp_gen_expand(mpt,m,sline,smap)
     struct pop_top *mpt;
     struct model_struct *m;    // Model
     char *sline;               // One line from file
     int ****smap;              // [mpt->nlay][xn][yn][zn]
{
  int i;
  int xn,yn,zn,ns,nnew,xi,yi,zi,nlay,***c,nsyn;
  float minw,fcrit,***wf,ori_0,ori_1;
  char **slist,*name_save,*dtype,*samp_str,*name_pre,*name_post,*name_pop;
  char *datid,*outfile,*optype;
  char ts[SLEN];
  struct pop_layer *lpre,*lpost,*lay,**llist;

  nnew = 0;

  get_items_from_string(sline,&slist,&ns);
  if (ns > 0){

    outfile = m->marfile;

    if (compare_prefix_string_order("gen_save_pop_lay_ori_as_",
				    slist[0])==1){
      //
      //  WYETH Dec 2013
      //
      // gen_save_pop_lay_ori_as_XXX [data_type] [sampling] [pop] ...
      //   ...  [op_type] [ori_ctr] [ori_max_dev]
      //
      // Where
      //     'op_type'     - "ori_180" - computeoperation type
      //     'ori_ctr'     - the central ori value
      //     'ori_max_dev' - the maximum deviation from center value
      //
      // e.g. "ori_ctr 0   ori_max_dev 20"  takes from -20 to 20.
      //

      if (ns != 8){
	printf("  ==>  %s\n",sline);
	exit_error("POPU_RSP_GEN_EXPAND","Format error, expecting 7 args");
      }
      name_save = (char *)(slist[0] + 24);
      dtype     =          slist[1];
      samp_str  =          slist[2];
      name_pop  =          slist[3];
      optype    =          slist[4];
      ori_0     =     atof(slist[5]);  // Central value of ori range
      ori_1     =     atof(slist[6]);  // Max deviation around central value
      datid     =          slist[7];
      printf("    Expanding 'gen_save_pop_lay_ori_as_'\n");

      //  Get pointer to the layer
      lay  = popu_rsp_gen_get_lay_ptr(mpt,name_pop);

      //  Get pointer to map of connections (or create if needed)
      c = popu_rsp_gen_get_map(mpt,smap,lay);
      nnew = popc_flag_accum_lay_ori(lay,c,optype,ori_0,ori_1);

      printf("      Including %d units from '%s'  (ori crit:  %.4f %.4f)\n",
	     nnew,name_save,ori_0,ori_1);

      //  Append lines to .rsp file
      popu_rsp_gen_append_by_map(c,lay->xn,lay->yn,lay->zn,outfile,
				 name_save,dtype,samp_str,name_pre,datid,
				 sline,nnew,NULL);

    }else if (compare_prefix_string_order("gen_save_pop_lay_to_unit_as_",
				    slist[0])==1){
      if (ns != 10){
	printf("  ==>  %s\n",sline);
	exit_error("POPU_RSP_GEN_EXPAND","Format error, expecting 9 args");
      }
      name_save = (char *)(slist[0] + 28);
      dtype     =          slist[1];
      samp_str  =          slist[2];
      name_pre  =          slist[3];
      name_post =          slist[4];
      xi        =     atoi(slist[5]);
      yi        =     atoi(slist[6]);
      zi        =     atoi(slist[7]);
      minw      =     atof(slist[8]);
      datid     =          slist[9];
      printf("    Expanding 'gen_save_pop_lay_to_unit_as_'\n");

      //  Get pointers to the pre- and post-syn layers
      lpre  = popu_rsp_gen_get_lay_ptr(mpt,name_pre);
      lpost = popu_rsp_gen_get_lay_ptr(mpt,name_post);

      //  Get pointer to map of connections
      c = popu_rsp_gen_get_map(mpt,smap,lpre);
      nnew = popc_flag_accum_lay_to_unit(lpre,lpost,xi,yi,zi,c,minw);

      printf("      Including %d units from '%s'  (minw %.4f)\n",
	     nnew,name_save,minw);

      //  Append lines to .rsp file
      popu_rsp_gen_append_by_map(c,lpre->xn,lpre->yn,lpre->zn,outfile,
				 name_save,dtype,samp_str,name_pre,datid,
				 sline,nnew,NULL);

    }else if (compare_prefix_string_order("gen_save_pop_lay_to_lay_as_",
					  slist[0])==1){
      if (ns != 7){
	printf("  ==>  %s\n",sline);
	exit_error("POPU_RSP_GEN_EXPAND","Format error, expecting 6 args");
      }
      name_save = (char *)(slist[0] + 27);
      dtype     =          slist[1];
      samp_str  =          slist[2];
      name_pre  =          slist[3];
      name_post =          slist[4];
      minw      =     atof(slist[5]);
      datid     =          slist[6];
      printf("    Expanding 'gen_save_pop_lay_to_lay_as_'\n");

      //  Get pointers to the pre- and post-syn layers
      lpre  = popu_rsp_gen_get_lay_ptr(mpt,name_pre);
      lpost = popu_rsp_gen_get_lay_ptr(mpt,name_post);

      //  Get pointer to map of connections
      c = popu_rsp_gen_get_map(mpt,smap,lpre);
      nnew = popc_flag_accum_lay_to_lay(lpre,lpost,c,minw);

      printf("      Including %d units from '%s'  (minw %.4f)\n",
	     nnew,name_save,minw);

      //  Append lines to .rsp file
      popu_rsp_gen_append_by_map(c,lpre->xn,lpre->yn,lpre->zn,outfile,
				 name_save,dtype,samp_str,name_pre,datid,
				 sline,nnew,NULL);

    }else if (compare_prefix_string_order("gen_save_pop_multi_lay_as_",
					  slist[0])==1){
      if (ns < 8){
	printf("  ==>  %s\n",sline);
	exit_error("POPU_RSP_GEN_EXPAND","Format error, expecting 8+ args");
      }
      name_save = (char *)(slist[0] + 26);
      dtype     =          slist[1];
      samp_str  =          slist[2];
      // layer list is from 3 to ns-4
      optype    =          slist[ns-3];
      fcrit     =     atof(slist[ns-2]);
      datid     =          slist[ns-1];

      nlay = ns - 6;  // Must be at least two layers
      llist = (struct pop_layer **)myalloc(nlay*sizeof(struct pop_layer *));
      for(i=0;i<nlay;i++){
	llist[i] = popu_rsp_gen_get_lay_ptr(mpt,slist[i+3]);
	//printf("%d  %s\n",i,llist[i]->name);
      }

      lpre = llist[0];
      xn = lpre->xn;
      yn = lpre->yn;
      zn = lpre->zn;

      printf("    Expanding 'gen_save_pop_multi_lay_as_'\n");

      //  Get aggregate weight array
      wf = popc_wf_accum_multi_lay(llist,nlay,-1,-1,-1,optype,fcrit);

      //  Get map of cells to output *** WYETH DO WE WANT THIS TO BE PART
      //  OF THE STORED MAPS ???

      c = popu_rsp_gen_cmap_from_w(wf,xn,yn,zn,optype,fcrit,&nnew);

      printf("      Including %d units from '%s'  (fcrit %.4f)\n",
	     nnew,name_save,fcrit);

      //  Append lines to .rsp file
      popu_rsp_gen_append_by_map(c,lpre->xn,lpre->yn,lpre->zn,outfile,
				 name_save,dtype,samp_str,name_pre,datid,
				 sline,nnew,wf);
    }else if (compare_prefix_string_order("gen_save_pop_mlay_unit_as_",
					  slist[0])==1){
      if (ns < 11){
	printf("  ==>  %s\n",sline);
	exit_error("POPU_RSP_GEN_EXPAND","Format error, expecting 11+ args");
      }
      name_save = (char *)(slist[0] + 26);
      dtype     =          slist[1];
      samp_str  =          slist[2];
      // layer list is from 3 to ns-7
      xi        =     atoi(slist[ns-6]);
      yi        =     atoi(slist[ns-5]);
      zi        =     atoi(slist[ns-4]);
      optype    =          slist[ns-3];
      fcrit     =     atof(slist[ns-2]);
      datid     =          slist[ns-1];

      nlay = ns - 9;  // Must be at least two layers
      llist = (struct pop_layer **)myalloc(nlay*sizeof(struct pop_layer *));
      for(i=0;i<nlay;i++){
	llist[i] = popu_rsp_gen_get_lay_ptr(mpt,slist[i+3]);
	//printf("%d  %s\n",i,llist[i]->name);
      }

      lpre = llist[0];
      xn = lpre->xn;
      yn = lpre->yn;
      zn = lpre->zn;

      printf("    Expanding 'gen_save_pop_mlay_unit_as_'\n");

      //  Get aggregate weight array
      wf = popc_wf_accum_multi_lay(llist,nlay,xi,yi,zi,optype,fcrit);

      //  Get map of cells to output *** WYETH DO WE WANT THIS TO BE PART
      //  OF THE STORED MAPS ???

      c = popu_rsp_gen_cmap_from_w(wf,xn,yn,zn,optype,fcrit,&nnew);

      printf("      Including %d units from '%s'  (fcrit %.4f)\n",
	     nnew,name_save,fcrit);

      //  Append lines to .rsp file
      popu_rsp_gen_append_by_map(c,lpre->xn,lpre->yn,lpre->zn,outfile,
				 name_save,dtype,samp_str,name_pre,datid,
				 sline,nnew,wf);
    }

    free_2d_carray(slist,ns);
  }

  return nnew;
}
/**************************************-**************************************/
/*                                                                           */
/*                               POPU_RSP_GEN                                */
/*                                                                           */
/*****************************************************************************/
void popu_rsp_gen(mpt,m,r)
     struct pop_top *mpt;
     struct model_struct *m;        // Model
     struct response_struct *r;     // Response generating file
{
  int i;
  int ns;
  char **sdata,*outfile;
  int ****smap;

  printf("  POPU_RESP_GEN\n");

  //
  //  Hold "save maps" for each layer, indicating which units are recorded
  //
  //  Values in the save map:
  //    -1 - connection previously written for saving in .rsp outfile
  //     0 - not marked for save
  //    >0 - marked for save, but written yet
  //  
  smap = (int ****)myalloc(mpt->nlay * sizeof(int ***));
  for(i=0;i<mpt->nlay;i++)
    smap[i] = NULL;


  outfile = m->marfile;

  printf("    Writing response file:  %s\n",outfile);

  remove_file(outfile);  // Remove old version of outfile

  read_2d_carray(r->paramfile,1,1,&sdata,&ns);  // Read the .rsg file

  //
  //  For each line in the .rsg file, read and expand it, or just copy it.
  //
  for(i=0;i<ns;i++){
    if (popu_rsp_gen_expand(mpt,m,sdata[i],smap) == 0){
      //
      //  This line was not expanded
      //
      append_string_to_file(outfile,sdata[i]);  // Write it unchanged
      append_string_to_file(outfile,"\n");
    }
  }
}
