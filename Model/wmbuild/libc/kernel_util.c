/*****************************************************************************/
/*                                                                           */
/*  kernel_util.c                                                            */
/*  wyeth bair                                                               */
/*  caltech                                                                  */
/*  02/23/93                                                                 */
/*                                                                           */
/*  For multi-dimensional functions that reprsent components of neuronal     */
/*  receptive fields in space and time.                                      */
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
#include "misc_util.h"
#include "farray_util.h"
#include "func_util.h"
#include "myrand_util.h"
#include "data_util.h"
#include "paramfile_util.h"
#include "fft_util.h"

#define PFLAG 0

/**************************************-**************************************/
/*                                                                           */
/*                           KERNU_NORM_2D_FILTER                            */
/*                                                                           */
/*****************************************************************************/
void kernu_norm_2d_filter(f,xn,yn,targ,rsump,rsumn)
     float **f;            // [xn][yn]
     int xn,yn;            //
     float targ;           // Target value
     float *rsump,*rsumn;  // Return final positive and negative sums
{
  int i,j;
  float sump,sumn,tt;

  // Compute POS and NEG SUM
  sump = sumn = 0.0;
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      tt = f[i][j];
      if (tt > 0.0)
	sump += tt;
      else
	sumn += tt;
    }
  }
  multiply_2d_farray(f,xn,yn,targ/(sump-sumn));

  //printf("KERNEL_UTIL:  sump = %f  sumn = %f\n",sump,sumn);

  // Repeat the above, compute pos and neg sum
  sump = sumn = 0.0;
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      tt = f[i][j];
      if (tt > 0.0)
	sump += tt;
      else
	sumn += tt;
    }
  }

  *rsump = sump;
  *rsumn = sumn;
}
/**************************************-**************************************/
/*                                                                           */
/*                         KERNEL_UTIL_NORM_XYT_FILTER                       */
/*                                                                           */
/*****************************************************************************/
void kernel_util_norm_xyt_filter(f,x0,xn,y0,yn,z0,zn,rsump,rsumn)
     float ***f;
     int x0,xn,y0,yn,z0,zn;
     float *rsump,*rsumn;
{
  int i,j,k;
  float sump,sumn,tt;

  // Compute POS and NEG SUM
  sump = sumn = 0.0;
  for(i=x0;i<=xn;i++){
    for(j=y0;j<=yn;j++){
      for(k=z0;k<=zn;k++){
	tt = f[i][j][k];
	if (tt > 0.0)
	  sump += tt;
	else
	  sumn += tt;
      }
    }
  }
  multiply_3d_farray(f,x0,xn,y0,yn,z0,zn,1.0/sump);

  // pos and neg sum
  sump = sumn = 0.0;
  for(i=x0;i<=xn;i++){
    for(j=y0;j<=yn;j++){
      for(k=z0;k<=zn;k++){
	tt = f[i][j][k];
	if (tt > 0.0)
	  sump += tt;
	else
	  sumn += tt;
      }
    }
  }

  *rsump = sump;
  *rsumn = sumn;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 RECTANGLE                                 */
/*                                                                           */
/*****************************************************************************/
float *rectangle(n)
     int n;
{
  int i;
  float *data;
  
  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = 1.0/(float)n;
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 DECAY_EXP                                 */
/*                                                                           */
/*****************************************************************************/
float *decay_exp(n,l1)
     int n;
     float l1;
{
  int i;
  float *data,t;
  
  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++){
    t = (float)i/1000.0;
    data[i] = l1*exp(-l1*t)/1000.0;
  }
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                                    DOE                                    */
/*                                                                           */
/*  Difference of Exponential functions.                                     */
/*                                                                           */
/*****************************************************************************/
float *doe(n,l1,l2)
     int n;
     float l1,l2;
{
  int i;
  float *data,t;
  
  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++){
    t = (float)i/1000.0;
    data[i] = (l1*exp(-l1*t) - l2*exp(-l2*t))/1000.0;
  }
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                                   ALPHA                                   */
/*                                                                           */
/*  Peak is at 1/alpha.                                                      */
/*                                                                           */
/*****************************************************************************/
float *alpha(n,alpha)
     int n;
     float alpha;
{
  int i;
  float *data,t;
  
  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++){
    t = (float)i;
    data[i] = alpha*alpha*t*exp(-alpha*t);
  }
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                                    DOA                                    */
/*                                                                           */
/*  Difference of Alpha functions.                                           */
/*                                                                           */
/*****************************************************************************/
float *doa(n,a1,a2)
     int n;
     float a1,a2;
{
  int i;
  float *data,t;
  
  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++){
    t = (float)i/1000.0;
    data[i] = (a1*a1*t*exp(-a1*t) - a2*a2*t*exp(-a2*t))/1000.0;
  }
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            KERNU_O_FILTER_BOXCAR                          */
/*                                                                           */
/*****************************************************************************/
float *kernu_o_filter_boxcar(mylogf,fo,tscale,rn)
     char mylogf[];
     struct onode *fo;  // <filter>
     float tscale;      // sec / time unit
     int *rn;           // if -1, return filter length
{
  int i;
  int n,ft0,ftn,ftn2;
  float famp,famp2,*fdata;

  if (*rn == -1)
    n     = onode_getpar_int_exit(fo,"n");
  else
    n = *rn;
  ft0    = onode_getpar_int_exit(fo,"t0");
  ftn    = onode_getpar_int_exit(fo,"tn");
  famp   = onode_getpar_flt_exit(fo,"amp");

  ftn2 = 0;
  if (onode_item(fo,"tn2")==1){               // 2nd lobe
    ftn2  = onode_getpar_int_exit(fo,"tn2");
    famp2 = onode_getpar_flt_exit(fo,"amp2");
  }

  if (ft0 < 0)
    mylogx(mylogf,"KERNU_O_FILTER_BOXCAR","ft0 < 0 for <filter>");
  if (ft0+ftn+ftn2 > n)
    mylogx(mylogf,"KERNU_O_FILTER_BOXCAR",
	   "boxcar too long for 'n' in <filter>");

  fdata = boxcar_farray((float)ft0,(float)ftn,famp,n);
  if (ftn2 > 0){
    for(i=(ft0+ftn);i<(ft0+ftn+ftn2);i++)
      fdata[i] = famp2;
  }

  *rn = n;
  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                             KERNU_O_FILTER_EXP                            */
/*                                                                           */
/*****************************************************************************/
float *kernu_o_filter_exp(mylogf,fo,tscale,rn)
     char mylogf[];
     struct onode *fo;  // <filter>
     float tscale;      // sec / time unit
     int *rn;           // if -1, return filter length
{
  int n;
  float t0,t1,t2,amp2,a,tau1,tau2,*fdata,*d1,*d2,xoff,gsd,gsds,*sm;
  char *orig_str;

  if (*rn == -1)
    n = onode_getpar_int_exit(fo,"n");
  else
    n = *rn;
  t1     = onode_getpar_flt_exit(fo,"tau1");      // Peak time (s)
  t2     = onode_getpar_flt_dflt(fo,"tau2",-1.0); // Peak time (s)
  amp2   = onode_getpar_flt_dflt(fo,"amp2",0.0);  // relative to first lobe
  xoff   = onode_getpar_flt_dflt(fo,"xoff",0.0);  // x-offset (s)
  gsd    = onode_getpar_flt_dflt(fo,"gauss_sd",0.0);  // Optional Gauss smooth

  orig_str = onode_getpar_chr_dflt(fo,"origin","zero");
  if (strcmp(orig_str,"center")==0){
    t0 = (float)(int)((n-1)/2.0);
    //t0 = (float)(n-1)/2.0;  // WYETH CHANGED Jan 13, 2015

    /*** WYETH - I THINK THERE IS STILL AN ISSUE w/ the origin.  Why isn't
     *** WYETH - I THINK THERE IS STILL AN ISSUE w/ the origin.  Why isn't
     *** WYETH - I THINK THERE IS STILL AN ISSUE w/ the origin.  Why isn't
     the function 1.0 at time zero???? or is it?  ***/

  }else if (strcmp(orig_str,"zero")==0){
    t0 = 0.0;
  }else
    exit_error("KERNU_O_FILTER_EXP","Unkown 'origin' string");

  t0 += xoff / tscale;

  tau1 = t1/tscale;  // convert (s) -> (sampling units)
  tau2 = t2/tscale;  // convert (s) -> (sampling units)

  //d1 = one_sided_exp_farray(0.0,tau1,1.0,n);
  d1 = one_sided_exp_farray(t0,tau1,1.0,n);

  if ((amp2 != 0.0) && (t2 != -1.0)){   // If there is a second function
    d2 = one_sided_exp_farray(0.0,tau2,amp2,n);
    fdata = add_farrays(d1,d2,n);
    myfree(d1);
    myfree(d2);
  }else
    fdata = d1;

  if (gsd > 0.0){
    //  Smooth the trace with a Gaussian
    gsds = gsd/tscale;  // convert (s) -> (sampling units)
    sm = smooth_with_gaussian(fdata,n,gsds,0.01);
    myfree(fdata);
    fdata = sm;
  }

  if (onode_item(fo,"max")==1){
    a = onode_getpar_flt_exit(fo,"max");
    make_max_const_farray(fdata,n,a);
  }else if (onode_item(fo,"area")==1){
    a = onode_getpar_flt_exit(fo,"area");
    norm_area_farray(fdata,n,a);
  }else{
    norm_area_farray(fdata,n,1.0);  // Default
  }

  *rn = n;
  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                            KERNU_O_FILTER_ALPHA                           */
/*                                                                           */
/*****************************************************************************/
float *kernu_o_filter_alpha(mylogf,fo,tscale,rn)
     char mylogf[];
     struct onode *fo;  // <filter>
     float tscale;      // sec / time unit
     int *rn;           // if -1, return filter length
{
  int n;
  float t0,t1,t2,amp2,a,tau1,tau2,*fdata,*d1,*d2,xoff,gsd,gsds,*sm;
  char *orig_str;

  if (*rn == -1)
    n = onode_getpar_int_exit(fo,"n");
  else
    n = *rn;
  t1     = onode_getpar_flt_exit(fo,"tau1");      // Peak time (s)
  t2     = onode_getpar_flt_dflt(fo,"tau2",-1.0); // Peak time (s)
  amp2   = onode_getpar_flt_dflt(fo,"amp2",0.0);  // relative to first lobe
  xoff   = onode_getpar_flt_dflt(fo,"xoff",0.0);  // x-offset (s)
  gsd    = onode_getpar_flt_dflt(fo,"gauss_sd",0.0);  // Optional Gauss smooth

  orig_str = onode_getpar_chr_dflt(fo,"origin","zero");
  if (strcmp(orig_str,"center")==0){
    t0 = (float)(int)((n-1)/2.0);
  }else if (strcmp(orig_str,"zero")==0){
    t0 = 0.0;
  }else
    exit_error("KERNU_O_FILTER_ALPHA","Unkown 'origin' string");

  t0 += xoff / tscale;

  tau1 = t1/tscale;  // convert (s) -> (sampling units)
  tau2 = t2/tscale;  // convert (s) -> (sampling units)

  //d1 = one_sided_exp_farray(t0,tau1,1.0,n);
  d1 = alpha_farray(t0,1.0/tau1,1.0,n);

  if ((amp2 != 0.0) && (t2 != -1.0)){   // If there is a second function
    d2 = alpha_farray(t0,1.0/tau2,amp2,n);
    fdata = add_farrays(d1,d2,n);
    myfree(d1);
    myfree(d2);
  }else
    fdata = d1;

  if (gsd > 0.0){
    //  Smooth the trace with a Gaussian
    gsds = gsd/tscale;  // convert (s) -> (sampling units)
    sm = smooth_with_gaussian(fdata,n,gsds,0.01);
    myfree(fdata);
    fdata = sm;
  }

  if (onode_item(fo,"max")==1){
    a = onode_getpar_flt_exit(fo,"max");
    make_max_const_farray(fdata,n,a);
  }else if (onode_item(fo,"area")==1){
    a = onode_getpar_flt_exit(fo,"area");
    norm_area_farray(fdata,n,a);
  }else{
    norm_area_farray(fdata,n,1.0);  // Default
  }

  *rn = n;
  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                           KERNU_O_FILTER_GAUSSIAN                         */
/*                                                                           */
/*****************************************************************************/
float *kernu_o_filter_gaussian(mylogf,fo,tscale,rn)
     char mylogf[];
     struct onode *fo;  // <filter>
     float tscale;      // sec / time unit
     int *rn;           // return filter length
{
  int n;
  float m1,m2,s1,s2,amp2,a,mu1,mu2,sd1,sd2,*fdata,*d1,*d2;

  if (*rn == -1)
    n = onode_getpar_int_exit(fo,"n");
  else
    n = *rn;
  m1     = onode_getpar_flt_exit(fo,"mean");        // mean (s)
  s1     = onode_getpar_flt_exit(fo,"sd");          // SD (s)
  m2     = onode_getpar_flt_dflt(fo,"mean2",-1.0);  // mean (s)
  s2     = onode_getpar_flt_dflt(fo,"sd2",-1.0);    // SD (s)
  amp2   = onode_getpar_flt_dflt(fo,"amp2",0.0);    // relative to first lobe

  mu1 = m1/tscale;  // Convert (s) -> (sampling units)
  mu2 = m2/tscale;
  sd1 = s1/tscale;
  sd2 = s2/tscale;

  d1 = gaussian_farray(mu1,sd1,1.0,n);

  if ((amp2 != 0.0) && (m2 != -1.0)){   // If there is a second function
    d2 = gaussian_farray(mu2,sd2,amp2,n);
    fdata = add_farrays(d1,d2,n);
    myfree(d1);
    myfree(d2);
  }else
    fdata = d1;

  if (onode_item(fo,"max")==1){
    a = onode_getpar_flt_exit(fo,"max");
    make_max_const_farray(fdata,n,a);
  }else if (onode_item(fo,"area")==1){
    a = onode_getpar_flt_exit(fo,"area");
    norm_area_farray(fdata,n,a);
  }else{
    norm_area_farray(fdata,n,1.0);  // Default
  }

  *rn = n;
  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                           KERNU_O_FILTER_MAXWELL                          */
/*                                                                           */
/*****************************************************************************/
float *kernu_o_filter_maxwell(mylogf,fo,tscale,rn)
     char mylogf[];
     struct onode *fo;  // <filter>
     float tscale;      // sec / time unit
     int *rn;           // if -1, return filter length
{
  int n,hidnum;
  float m1,m2,amp2,a,s1,s2,*fdata,*d1,*d2,t0,xoff;
  char *orig_str;

  if (*rn == -1)
    n = onode_getpar_int_exit(fo,"n");
  else
    n = *rn;

  hidnum  = onode_getpar_int_dflt(fo,"hidden",0);   // Hidden code
  if (hidnum == 0){
    m1     = onode_getpar_flt_exit(fo,"m1");        // Peak time (s)
    m2     = onode_getpar_flt_dflt(fo,"m2",-1.0);   // Peak time (s)
    xoff   = onode_getpar_flt_dflt(fo,"xoff",0.0);  // x-offset (s)
    amp2   = onode_getpar_flt_dflt(fo,"amp2",0.0);  // relative to first lobe
  }else if (hidnum == 1){
    m1     =  0.040;   // Peak time (s)
    m2     =  0.060;   // Peak time (s)
    xoff   =  0.0;     // x-offset (s)
    amp2   = -0.80;    // relative to first lobe
  }else
    mylog_exit(mylogf,"KERNU_O_FILTER_MAXWELL  Bad 'hidden' index\n");

  orig_str = onode_getpar_chr_dflt(fo,"origin","zero");
  if (strcmp(orig_str,"center")==0){
    t0 = (float)(int)((n-1)/2.0);
    //t0 = (float)(n-1)/2.0; // WYETH CHANGED to make center be on a point
  }else if (strcmp(orig_str,"zero")==0){
    t0 = 0.0;
  }else
    exit_error("KERNU_O_FILTER_MAXWELL","Unkown 'origin' string");

  t0 += xoff / tscale;

  s1 = m1/(tscale*sqrt(2.0));  // peak = m1 = sqrt(2)*s;
  s2 = m2/(tscale*sqrt(2.0));  // peak = m2 = sqrt(2)*s;

  d1 = maxwell_farray(t0,s1,1.0,n);
  //d1 = maxwell_farray(0.0,s1,1.0,n);

  if ((amp2 != 0.0) && (m2 != -1.0)){   // If there is a second function
    d2 = maxwell_farray(t0,s2,amp2,n);
    //d2 = maxwell_farray(0.0,s2,amp2,n);
    fdata = add_farrays(d1,d2,n);
    myfree(d1);
    myfree(d2);
  }else
    fdata = d1;

  // **** WYETH HERE ********
  // **** WYETH HERE ********
  // **** WYETH 
  // **** WYETH  ALL of these, and 'origin' should be done in caller???
  // **** WYETH 
  // **** WYETH 
  if (onode_item(fo,"max")==1){
    a = onode_getpar_flt_exit(fo,"max");
    make_max_const_farray(fdata,n,a);
  }else if (onode_item(fo,"area")==1){
    a = onode_getpar_flt_exit(fo,"area");
    norm_area_farray(fdata,n,a);
  }else{
    norm_area_farray(fdata,n,1.0);  // Default
  }

  *rn = n;
  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                           KERNU_O_FILT2D_GAUSSIAN                         */
/*                                                                           */
/*****************************************************************************/
float **kernu_o_filt2d_gaussian(mylogf,fo,o1,o2,scale1,scale2,rxn,ryn)
     char mylogf[];
     struct onode *fo;     // <filter>
     float o1,o2;          // origin within xn,yn grid
     float scale1,scale2;  // units/pix for 1st and 2nd dimensions
     int *rxn;             // if -1, return filter length, else use this value
     int *ryn;             // if -1, return filter length, else use this value
{
  int i,j;
  int xn,yn;
  float m1,m2,s1,s2,x,y,a,**fdata;

  if (*rxn == -1){
    xn = onode_getpar_int_exit(fo,"xn");
  }else
    xn = *rxn;

  if (*ryn == -1){
    yn = onode_getpar_int_exit(fo,"yn");
  }else
    yn = *ryn;

  fdata = get_2d_farray(xn,yn);

  m1  = onode_getpar_flt_exit(fo,"mu_1");
  m2  = onode_getpar_flt_exit(fo,"mu_2");
  s1  = onode_getpar_flt_exit(fo,"sd_1");
  s2  = onode_getpar_flt_exit(fo,"sd_2");

  //printf("s1 s2  %f %f   m1,2 %f %f\n",s1,s2,m1,m2);
  //printf("scale1,2  %f %f\n",scale1,scale2);

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      x = ((float)i-o1)*scale1;
      y = ((float)j-o2)*scale2;
      fdata[i][j] = func_2d_gaussian(x,y,m1,m2,s1,s2,1);
    }
  }

  if (onode_item(fo,"max")==1){
    a = onode_getpar_flt_exit(fo,"max");
    make_max_const_2d_farray(fdata,xn,yn,a);
  }else if (onode_item(fo,"area")==1){
    a = onode_getpar_flt_exit(fo,"area");
    norm_area_2d_farray(fdata,xn,yn,a);
  }else{
    norm_area_2d_farray(fdata,xn,yn,1.0);  // Default
  }

  *rxn = xn;
  *ryn = yn;
  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                               KERNU_O_FILTER                              */
/*                                                                           */
/*  Create and return the filter described by the <filter> object.           */
/*                                                                           */
/*****************************************************************************/
float *kernu_o_filter(mylogf,fo,tscale,rn)
     char mylogf[];
     struct onode *fo;  // <filter>
     float tscale;      // sec / time unit
     int *rn;           // if -1, return filter length, else use this value
{
  char *ftype,*dumpfile;
  float *fdata;

  mylog(mylogf,"  KERNU_O_FILTER\n");

  //
  //  <filter> objects with 'dimension 1d' should come here.  Older models
  //  (before 'dimension') should also come here.
  //

  if (fo == NULL)
    mylogx(mylogf,"KERNU_O_FILTER","Null filter onode");

  ftype = onode_getpar_chr_exit(fo,"type");
  if (strcmp(ftype,"boxcar")==0){
    fdata = kernu_o_filter_boxcar(mylogf,fo,tscale,rn);
  }else if (strcmp(ftype,"exp")==0){
    fdata = kernu_o_filter_exp(mylogf,fo,tscale,rn);
  }else if (strcmp(ftype,"Gaussian")==0){
    fdata = kernu_o_filter_gaussian(mylogf,fo,tscale,rn);
  }else if (strcmp(ftype,"Maxwell")==0){
    fdata = kernu_o_filter_maxwell(mylogf,fo,tscale,rn);
  }else if (strcmp(ftype,"alpha")==0){
    fdata = kernu_o_filter_alpha(mylogf,fo,tscale,rn);
  }else{
    mylogx(mylogf,"KERNU_O_FILTER","Unknown 1d <filter> type");
  }

  dumpfile = onode_getpar_chr_dflt(fo,"plot_filename",NULL);
  if (dumpfile != NULL){
    if ((strcmp(dumpfile,"null")!=0)&&(strcmp(dumpfile,"NULL")!=0))
      append_farray_plot(dumpfile,"filter(raw_time_units)",fdata,*rn,1);
    myfree(dumpfile);
  }
 
  myfree(ftype);

  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                              KERNU_O_FILT2D                               */
/*                                                                           */
/*  Create and return the filter described by the <filter> object.           */
/*                                                                           */
/*****************************************************************************/
float **kernu_o_filt2d(mylogf,fo,o1,o2,scale1,scale2,rxn,ryn)
     char mylogf[];
     struct onode *fo;    // <filter>
     float o1,o2;         // origin of coords within (xn,yn)
     float scale1,scale2; // units/pix in 1st and 2nd dimensions
     int *rxn,*ryn;       // if -1, return filter size, else use these values
{
  char *ftype,*dumpfile;
  float **fdata;

  mylog(mylogf,"  KERNU_O_FILT2D\n");

  if (fo == NULL)
    mylogx(mylogf,"KERNU_O_FILT2D","Null filter onode");

  ftype = onode_getpar_chr_exit(fo,"type");
  if (strcmp(ftype,"Gaussian2d")==0){
    fdata = kernu_o_filt2d_gaussian(mylogf,fo,o1,o2,scale1,scale2,rxn,ryn);
  }else{
    mylogx(mylogf,"KERNU_O_FILT2D","Unknown 2d <filter> type");
  }

  dumpfile = onode_getpar_chr_dflt(fo,"plot_filename",NULL);
  if (dumpfile != NULL){
    if ((strcmp(dumpfile,"null")!=0)&&(strcmp(dumpfile,"NULL")!=0))
      write_2d_data(dumpfile,fdata,0,0,*rxn,*ryn,4,2,1,0);
    myfree(dumpfile);
  }
 
  myfree(ftype);

  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                                WA_TEMPORAL                                */
/*                                                                           */
/*  The Watson and Ahumada (1985) temporal filter (p327).                    */
/*                                                                           */
/*    f(t) = a*[f1(t) - b*f2(t)]                                             */
/*                                                                           */
/*                  u(t)                                                     */
/*    fi(t) = ---------------- (t/tua_i)^(n_i - 1) exp(-t/tua_i)             */
/*            tau_i (n_i - 1)!                                               */
/*                                                                           */
/*  Plausible parameters (Robson, 1966)                                      */
/*                                                                           */
/*    b = 0.9                                                                */
/*    tau_1 = 0.004                                                          */
/*    tau_2 = 0.0053                                                         */
/*    n1 = 9                                                                 */
/*    n2 = 10                                                                */
/*                                                                           */
/*****************************************************************************/
float *wa_temporal(zn,tscale,a,b,n1,n2,tau1,tau2)
     int zn;
     float tscale;
     float a,b;
     int n1,n2;
     float tau1,tau2;
{
  int i;
  float *filter;
  float n1m1fact,n2m1fact;
  float t,f1,f2,ttau,n1m1,n2m1;
  int tc; /* centers of dimensions */
  int phase;

  if (PFLAG) printf("  WA_TEMPORAL\n");
  phase = 0;

  tc = (zn-1)/2;
  n1m1 = (float)(n1-1);
  n2m1 = (float)(n2-1);

  n1m1fact = 1;
  for(i=2;i<n1;i++)
    n1m1fact *= (float)i;
  n2m1fact = 1;
  for(i=2;i<n2;i++)
    n2m1fact *= (float)i;

  filter = (float *)myalloc(zn*sizeof(float));
  if (phase==0){  /*** Even "cos" filter ***/
    for(i=0;i<zn;i++){
      t = (float)(i-tc)*tscale;
      if (t >= 0.0){
	ttau = t/tau1;
	f1 = pow(ttau,n1m1) / (tau1*n1m1fact) * exp(-ttau);
	ttau = t/tau2;
	f2 = pow(ttau,n2m1) / (tau2*n2m1fact) * exp(-ttau);
      }else{
	f1 = 0.0;
	f2 = 0.0;
      }
      filter[i] = a*(f1 - b*f2);
    }
  }else{
    printf("*** WA_TEMPORAL:  No such option.\n");
    exit(0);
  }
  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                                AB_TEMPORAL                                */
/*                                                                           */
/*  The Adelson and Bergen (1985) temporal filter (p291).                    */
/*                                                                           */
/*    f(t) = (kt)^n exp(-kt) [1/n! - (kt)^2/(n+2)!]                          */
/*                                                                           */
/*  Plausible parameters (???)                                               */
/*                                                                           */
/*    n = 3 and 5                                                            */
/*    k = 100.0                                                              */
/*                                                                           */
/*****************************************************************************/
float *ab_temporal(zn,tscale,n,k,toff)
     int zn;        // duration of filter in sampling units
     float tscale;  // seconds per sample
     int n;
     float k;
     float toff;
{
  int i;
  int tc;
  float *filter,nfact,np2fact,kt,t;

  if (PFLAG) printf("  AB_TEMPORAL\n");

  tc = (zn-1)/2;  // centers of dimensions

  nfact = 1;
  for(i=2;i<=n;i++)
    nfact *= (float)i;
  np2fact = 1;
  for(i=2;i<=(n+2);i++)
    np2fact *= (float)i;

  filter = (float *)myalloc(zn*sizeof(float));
  for(i=0;i<zn;i++){
    t = (float)(i-tc)*tscale - toff;
    if (t >= 0.0){
      kt = k * t;
      filter[i] = pow(kt,(float)n) * exp(-kt) * (1.0/nfact - kt*kt/np2fact);
    }else{
      filter[i] = 0.0;
    }
  }
  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                                CWQ_TEMPORAL                               */
/*                                                                           */
/*  The Chen, Wang, Qian (2001) temporal filter (p144 Eqn 3).                */
/*                                                                           */
/*    f(t) = 1/[Gamma(a) tau^a]  t^(a-1) e^(-t/tau) cos(wt + phi),  t>=0     */
/*    f(t) = 0, t<0                                                          */
/*                                                                           */
/*  This is a Gamma function envelop times a cosine wave.                    */
/*                                                                           */
/*****************************************************************************/
float *cwq_temporal(int zn,        // duration of filter in sampling units
		    float tscale,  // Samples/sec  (1000 = msec)
		    int a,         // Integer value for Gamma func
		    float tau,     // Time const for envelope (s)
		    float fr,      // Frequency for cosine (cyc/s)
		    float phi,     // Phase for cosine (deg)
		    float toff){   // Temporal offset (s)
  int i;
  int tc; // centers of dimensions
  float *f,nfact,ph,tf;
  double t,c;

  if (PFLAG) printf("  CWQ_TEMPORAL\n");

  tc = (zn-1)/2;

  if (a < 1)
    exit_error("CWQ_TEMPORAL","Param 'a' must be >= 1");
  //
  //  Gamma(a) = (a-1)!
  //
  nfact = 1;
  for(i=2;i<a;i++)
    nfact *= (float)i;

  c = 1.0 / ((double)nfact * pow(tau,a));

  ph = phi * M_PI/180;  // convert to radians - phase offset
  tf = fr * 2.0*M_PI;   // convert to radians - temporal freq.

  f = (float *)myalloc(zn*sizeof(float));
  for(i=0;i<zn;i++){
    t = (double)(i-tc)*tscale - toff;  // Time (s)
    if (t >= 0.0){
      f[i] = (float)(c * pow(t,a-1.0) * exp(-t/tau) * cos(tf*t + ph));
    }else{
      f[i] = 0.0;
    }
  }
  return f;
}
/**************************************-**************************************/
/*                                                                           */
/*                              DD_EXP_TEMPORAL                              */
/*                                                                           */
/*****************************************************************************/
float *dd_exp_temporal(zn,tscale,t0,tau1,tau2,amp2,sig)
     int zn;        /* duration of filter in sampling units */
     float tscale;  /* sec per time unit */
     float t0;      /* center (s) */
     float tau1;    /* (s) */
     float tau2;    /* (s) */
     float amp2;    /* */
     float sig;     /* SD for smoothing (s) */
{
  int i;
  int tc;
  float *filter,t,*sm,tsig;

  printf("    DD_EXP_TEMPORAL\n");
  printf("    %f %f %f %f %f %f\n",tau1,tau2,t0,tscale,amp2,sig);

  tc = (zn-1)/2;

  if (PFLAG) printf("  AB_TEMPORAL\n");

  filter = (float *)myalloc(zn*sizeof(float));
  for(i=0;i<zn;i++){
    t = ((float)(i-tc))*tscale - t0;  /* Time in sec, relative to t0 */
    if ((t < 0.00001) && (t > -0.00001))
      t = 0.0; /* Make sure the ampl of the first (pos) lobe is 1 bef smooth */
    if (t <= 0.0)
      filter[i] = exp(t/tau1);
    else
      filter[i] = -amp2 * exp(-t/tau2);
  }

  tsig = sig/tscale; /* Convert to time units */
  sm = smooth_with_gaussian(filter,zn,tsig,0.01);
  myfree(filter);
  filter = sm;

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                              FILT_01_TEMPORAL                             */
/*                                                                           */
/*  Custom temporal filter.                                                  */
/*                                                                           */
/*****************************************************************************/
float *filt_01_temporal(ppl,o,zn,tscale,pflag)
     struct param_pair_list *ppl;
     struct onode *o;
     int zn;        /* duration of filter in sampling units */
     float tscale;  /* sec per time unit */
     int pflag;
{
  int i;
  int tc;
  float *filter,t,*sm,tsig,t0,namp,sig,tm;

  if (pflag) printf("    FILT_01_TEMPORAL\n");

  tc = (zn-1)/2;

  if (ppl != NULL){
    t0    = paramfile_get_float_param_default(ppl,"dogt_t0"    ,0.020);
    sig   = paramfile_get_float_param_default(ppl,"dogt_sigsm" ,0.001);
    namp  = paramfile_get_float_param_default(ppl,"dogt_negamp",0.100);
    tm    = paramfile_get_float_param_default(ppl,"dogt_tmult" ,100.0);
  }else{
    t0    = onode_getpar_flt_dflt(o,"dogt_t0"    ,0.020);
    sig   = onode_getpar_flt_dflt(o,"dogt_sigsm" ,0.001);
    namp  = onode_getpar_flt_dflt(o,"dogt_negamp",0.100);
    tm    = onode_getpar_flt_dflt(o,"dogt_tmult" ,100.0);
  }

  filter = (float *)myalloc(zn*sizeof(float));
  for(i=0;i<zn;i++){
    t = ((float)(i-tc))*tscale;  // Time in sec, relative to t0
    if (t < 0.0)
      filter[i] = 0.0;
    else if (t <= t0)
      filter[i] = 1.0 - (t/t0)*(1.0 + namp);
    else if (t > t0)
      filter[i] = -namp * (1.0/(1.0 + (t-t0)*tm));
  }

  if (sig > 0.0){
    tsig = sig/tscale;
    sm = smooth_with_gaussian(filter,zn,tsig,0.01);
    myfree(filter);
    filter = sm;
  }

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                                RANDOM_FILTER                              */
/*                                                                           */
/*****************************************************************************/
float ***random_filter(xn,yn,zn,sscale,tscale,sr,st,nsd_r,nsd_t,pseed)
     int xn,yn,zn;         // must be odd, center is (xn-1)/2,...,(zn-1/2)
     float sscale,tscale;  // space and time scaling
     float sr,st;
     float nsd_r,nsd_t;    // SD space, time for noise smoothing (deg),(sec)
     int *pseed;
{
  int i,j,k;
  int xc,yc,tc;  // centers of dimensions
  float ***filter,*tp,*sm,tsig;
  float c,cs,ct,csgn,x,y,t,x2,y2,gs,gt;

  if (PFLAG) printf("  RANDOM_FILTER\n");

  c = 1.0/(pow(2.0*M_PI,1.5)*sr*sr*st);
  cs = 2.0*sr*sr;
  ct = 2.0*st*st;
  xc = (xn-1)/2;
  yc = (yn-1)/2;
  tc = (zn-1)/2;

  tsig = nsd_t / tscale;      // convert to raw time units

  //printf("tsig = %f\n",tsig);

  filter = f3tensor(1,xn,1,yn,1,zn);  // Use tensor for 3D FFT

  for(i=0;i<xn;i++){
    x = (float)(i-xc)*sscale;
    x2 = x*x/cs;
    for(j=0;j<yn;j++){
      y = (float)(j-yc)*sscale;
      y2 = y*y/cs;
      gs = exp(-x2-y2);

      // SEPARABLE in space and time
      /*  COMMENT OUT SIMILAR LINES IN for(k loop
      if (myrand_util_ran2(pseed) > 0.5)
	csgn = 1.0;
      else
	csgn = -1.0;
      */

      tp = &(filter[i+1][j+1][1]);
      for(k=0;k<zn;k++){
	t = (float)(k-tc)*tscale;
	
	gt = exp(-t*t/ct);   // Gaussian by default

	// ORIGINAL WAY (see above for separable filter)
	if (myrand_util_ran2(pseed) > 0.5)
	  csgn = 1.0;
	else
	  csgn = -1.0;
	
	tp[k] = csgn * c*gs*gt;
      }

      // WYETH - TEMPORAL SMOOTHING HERE
      sm = smooth_with_gaussian(tp,zn,tsig,0.05);
      for(k=0;k<zn;k++)
	tp[k] = sm[k];

      myfree(sm);
    }
  }


  // WYETH WE should smooth here
  printf("*** WYETH - spatial smoothing not implemented\n");
  printf("*** WYETH - spatial smoothing not implemented\n");


  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                                GAUSSIAN_2D                                */
/*                                                                           */
/*  Return a 2D Gaussian function.                                           */
/*                                                                           */
/*  G(x,y) =  ampl * exp(-((y-yc)^2 + (x-xc)^2)/(2*sigma^2)                  */
/*                                                                           */
/*****************************************************************************/
float **gaussian_2d(xn,yn,xc,yc,sigma,ampl)
     int xn,yn;
     float xc,yc,sigma;
     float ampl;              //  Set ampl to 0.0 for all zeros
{
  int i,j;
  float **filter,cs,x,y,x2,y2;

  if (ampl == 0.0){
    filter = get_zero_2d_farray(xn,yn);
  }else{
    cs = 2.0*sigma*sigma;
    filter = get_2d_farray(xn,yn);
    for(i=0;i<xn;i++){
      x = ((float)i-xc);
      x2 = x*x/cs;
      for(j=0;j<yn;j++){
	y = ((float)j-yc);
	y2 = y*y/cs;
	filter[i][j] = ampl * exp(-x2-y2);
      }
    }
  }
  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GAUSS_2D_NOISE_RAW                            */
/*                                                                           */
/*  Return a 2D Gaussian times random noise.                                 */
/*                                                                           */
/*  Returned values in [-1,1]                                                */
/*                                                                           */
/*****************************************************************************/
float **gauss_2d_noise_raw(xn,yn,xc,yc,sdo,sdp,theta,npix,seed)
     int xn,yn;       // size
     float xc,yc;     // Center (pix)
     float sdo,sdp;   // Gaussian SD orth and par to ori (pix)
     float theta;     // Orientation
     int npix;        // Pixel size of noise
     int seed;        // Randomization seed
{
  int i,j;
  float **filter,cso,csp,x,y,gs;
  float a,b,c,d,thetar,xx,yy; // for rotation

  cso = 2.0*sdo*sdo;
  csp = 2.0*sdp*sdp;

  filter = get_2d_farray(xn,yn);

  if (seed > 0)
    seed = -seed;

  //printf("seed = %d\n",seed);

  if (npix != 1){
    printf("  npix = %d\n",npix);
    exit_error("GAUSS_2D_NOISE_RAW","Only pixel size = 1 implemented.");
  }

  // Use rotation matrix
  thetar = -theta * M_PI/180.0;
  a =  cos(thetar); b = sin(thetar);
  c = -sin(thetar); d = cos(thetar);

  for(i=0;i<xn;i++){
    xx = (float)i-xc;
    for(j=0;j<yn;j++){
      yy = (float)j-yc;

      x = a*xx + c*yy;
      y = b*xx + d*yy;

      gs = exp(-x*x/cso - y*y/csp);

      if (gs > 0.01){
	if (myrand_util_ran2(&seed) > 0.5)
	  filter[i][j] = gs;
	else
	  filter[i][j] = -gs;
      }else{
	filter[i][j] = 0.0;
      }
    }
  }

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GABOR_2D_SPACE_RAW                            */
/*                                                                           */
/*  Return a 2D gabor function.  Center is specified.                        */
/*                                                                           */
/*  Returned values in [-1,1]                                                */
/*                                                                           */
/*****************************************************************************/
float **gabor_2d_space_raw(xn,yn,xc,yc,sdo,sdp,sf,theta,phase)
     int xn,yn;       // size
     float xc,yc;     // Center (pix)
     float sdo,sdp;   // Gaussian SD orth and par to ori (pix)
     float sf,theta;  // cyc/pix
     float phase;     // (degr)
{
  int i,j;
  float **filter,cso,csp,twopifr,x,y,gs,ph;
  float a,b,c,d,thetar,xx,yy; // for rotation

  cso = 2.0*sdo*sdo;
  csp = 2.0*sdp*sdp;
  twopifr = 2.0*M_PI*sf;
  ph = (phase/180.0 * M_PI);

  filter = get_2d_farray(xn,yn);

  // Use rotation matrix
  thetar = -theta * M_PI/180.0;
  a =  cos(thetar); b = sin(thetar);
  c = -sin(thetar); d = cos(thetar);

  for(i=0;i<xn;i++){
    xx = (float)i-xc;
    for(j=0;j<yn;j++){
      yy = (float)j-yc;

      x = a*xx + c*yy;
      y = b*xx + d*yy;

      gs = exp(-x*x/cso - y*y/csp);

      filter[i][j] = gs*cos(twopifr * x + ph);
    }
  }

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GABOR_2D_SPACE_RAW_IRREG                         */
/*                                                                           */
/*  Return a 1D array of values for a 2D gabor function evaluated at the     */
/*  coordinates specified.                                                   */
/*                                                                           */
/*****************************************************************************/
float *gabor_2d_space_raw_irreg(xf,yf,n,xc,yc,sdo,sdp,sf,theta,phase)
     float *xf,*yf;     // [n] x,y coords of points
     int n;           // number of points
     float xc,yc;     // Center (pix)
     float sdo,sdp;   // Gaussian SD orth and par to ori (pix)
     float sf,theta;  // cyc/pix
     float phase;     // (degr)
{
  int i;
  float *filter,cso,csp,twopifr,x,y,gs,ph;
  float a,b,c,d,thetar,xx,yy; // for rotation

  //printf("xc,yc = %d 

  cso = 2.0*sdo*sdo;
  csp = 2.0*sdp*sdp;
  twopifr = 2.0*M_PI*sf;
  ph = (phase/180.0 * M_PI);
  
  //filter = get_2d_farray(xn,yn);
  filter = (float *)myalloc(n*sizeof(float));

  // Use rotation matrix
  thetar = -theta * M_PI/180.0;
  a =  cos(thetar); b = sin(thetar);
  c = -sin(thetar); d = cos(thetar);

  for(i=0;i<n;i++){
    //xx = (float)i-xc;
    xx = xf[i] - xc;
    //yy = (float)j-yc;
    yy = yf[i] - yc;

    x = a*xx + c*yy;
    y = b*xx + d*yy;

    gs = exp(-x*x/cso - y*y/csp);
    filter[i] = gs*cos(twopifr * x + ph);
  }

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GAUSS_2D_SPACE_RAW                             */
/*                                                                           */
/*  Return a 2D gaussian function.  Center is specified.                     */
/*                                                                           */
/*****************************************************************************/
float **gauss_2d_space_raw(xn,yn,xc,yc,sdo,sdp,theta)
     int xn,yn;       // size
     float xc,yc;     // Center (pix)
     float sdo,sdp;   // Gaussian SD orth and par to ori (pix)
     float theta;
{
  int i,j;
  float **filter,cso,csp,x,y;
  float a,b,c,d,thetar,xx,yy;  // for rotation

  cso = 2.0*sdo*sdo;
  csp = 2.0*sdp*sdp;

  filter = get_2d_farray(xn,yn);

  // Use rotation matrix
  thetar = -theta * M_PI/180.0;
  a =  cos(thetar); b = sin(thetar);
  c = -sin(thetar); d = cos(thetar);

  for(i=0;i<xn;i++){
    xx = (float)i-xc;
    for(j=0;j<yn;j++){
      yy = (float)j-yc;

      x = a*xx + c*yy;
      y = b*xx + d*yy;
      filter[i][j] = exp(-x*x/cso - y*y/csp);
    }
  }

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                             DISK_2D_SPACE_RAW                             */
/*                                                                           */
/*  Return a 2D disk with the given diameter.  One within, zero outside.     */
/*                                                                           */
/*****************************************************************************/
float **disk_2d_space_raw(xn,yn,xc,yc,diam)
     int xn,yn;       // size
     float xc,yc;     // Center (pix)
     float diam;      // Diameter (pix)
{
  int i,j;
  float **filter,r2,x,y;

  r2 = (diam/2.0) * (diam/2.0);

  filter = get_2d_farray(xn,yn);

  for(i=0;i<xn;i++){
    x = (float)i-xc;
    for(j=0;j<yn;j++){
      y = (float)j-yc;

      if ((x*x + y*y) <= r2)
	filter[i][j] = 1.0;
      else
	filter[i][j] = 0.0;
    }
  }

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                             DOG_2D_SPACE_RAW                              */
/*                                                                           */
/*  Return a 2D DoG function.  Center is specified.                          */
/*                                                                           */
/*****************************************************************************/
float **dog_2d_space_raw(xn,yn,xc,yc,sdctr,sdsurr,asurr)
     int xn,yn;            // size
     float xc,yc;          // Center (pix)
     float sdctr,sdsurr;   // Gaussian SD for center and surround (pix)
     float asurr;          // Area of surround, relative to center
{
  float **dcent,**dsurr,**filter;

  dcent = gauss_2d_space_raw(xn,yn,xc,yc,sdctr,sdctr,0.0);
  dsurr = gauss_2d_space_raw(xn,yn,xc,yc,sdsurr,sdsurr,0.0);

  norm_area_2d_farray(dcent,xn,yn,1.0);
  norm_area_2d_farray(dsurr,xn,yn,asurr);

  filter = subtract_2d_farrays(dcent,dsurr,xn,yn);

  free_2d_farray(dcent,xn);
  free_2d_farray(dsurr,xn);

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                               GABOR_2D_SPACE                              */
/*                                                                           */
/*  Return a 2D gabor function.  Center is (xn-1)/2 and (yn-1)/2.            */
/*                                                                           */
/*  If 'sscale' is given in degrees per pixel, e.g. 0.05, then 'fr'          */
/*  is given in cycles per degree, and sr is given in degrees,               */
/*  where degrees means degrees of visual angle.                             */
/*                                                                           */
/*****************************************************************************/
float **gabor_2d_space(xn,yn,sscale,sr,fr,theta,phase)
     int xn,yn;
     float sscale,sr,fr,theta;
     float phase;  // (degr)
{
  int i,j;
  int xc,yc; // centers of dimensions
  float **filter,cs,nx,ny,twopifr,x,y,x2,y2,gs,frndotr,ph;



  // *** WYETH DOES ANYONE CALL THIS ANYMORE ???
  // *** WYETH DOES ANYONE CALL THIS ANYMORE ???
  //
  //   it was replaced in 'stf3d_01' on Feb 3, 2017
  //
  // *** WYETH DOES ANYONE CALL THIS ANYMORE ???
  // *** WYETH DOES ANYONE CALL THIS ANYMORE ???



  cs = 2.0*sr*sr;
  nx = cos(theta/180.0 * M_PI); // nx,ny is unit vector
  ny = sin(theta/180.0 * M_PI);
  twopifr = 2.0*M_PI*fr;
  xc = (xn-1)/2;
  yc = (yn-1)/2;
  ph = (phase/180.0 * M_PI);

  filter = get_2d_farray(xn,yn);

  for(i=0;i<xn;i++){
    x = (float)(i-xc)*sscale;
    x2 = x*x/cs;
    for(j=0;j<yn;j++){
      y = (float)(j-yc)*sscale;
      y2 = y*y/cs;
      frndotr = twopifr * (nx*x + ny*y);
      gs = exp(-x2-y2);
      filter[i][j] = gs*cos(frndotr + ph);
    }
  }
  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                              GABOR_SPACE_TIME                             */
/*                                                                           */
/*  Return a Gabor function in 2 spatial dimensions and time.  The           */
/*  following complex function has the even Gabor filter as its real         */
/*  part and the odd Gabor filter as its imaginary part.  Note that          */
/*  only one of these parts is returned, depending on the value of           */
/*  "phase".                                                                 */
/*                                                                           */
/*  F(r,t:Or,n,Ot,sr,st) = 1/((2Pi)^(3/2) sr^2 st) exp(-|r|^2/(2sr^2))       */
/*                         exp(-iOr n.r)exp(-t^2/(2st^2))                    */
/*                         exp(-iOt t)                                       */
/*                                                                           */
/*  where                                                                    */
/*    r  spatial location vector                                             */
/*    t  time                                                                */
/*    sr spatial spread                                                      */
/*    st temporal spread                                                     */
/*    Or (fr) spatial frequency                                              */
/*    Ot (ft) temporal frequency                                             */
/*    n  unit vector in the direction of motion sensitivity                  */
/*       [cos(theta),sin(theta)]                                             */
/*                                                                           */
/*  Assume:                                                                  */
/*    sscale  given in degrees/pixel                                         */
/*    tscale  given in seconds/pixel (or seconds/frame)                      */
/*  and all other units in degrees and seconds where applicable.             */
/*                                                                           */
/*****************************************************************************/
float ***gabor_space_time(xn,yn,zn,sscale,tscale,sr,st,fr,ft,theta,phase)
     int xn,yn,zn; /* must be odd numbers, center is (xn-1)/2,...,(zn-1/2) */
     float sscale,tscale; /* space and time scaling */
     float sr,st,fr,ft,theta;
     int phase;  /* 0 for even, 1 for odd */
{
  int i,j,k;
  float ***filter;
  float c,cs,ct,nx,ny,twopifr,twopift;
  float x,y,t,x2,y2,gs,gt,frndotr;
  int xc,yc,tc; /* centers of dimensions */

  if (PFLAG) printf("  GABOR_SPACE_TIME\n");

  c = 1.0/(pow(2.0*M_PI,1.5)*sr*sr*st);
  cs = 2.0*sr*sr;
  ct = 2.0*st*st;
  nx = cos(theta/180.0 * M_PI);
  ny = sin(theta/180.0 * M_PI);
  twopifr = 2.0*M_PI*fr;
  twopift = 2.0*M_PI*ft;
  xc = (xn-1)/2;
  yc = (yn-1)/2;
  tc = (zn-1)/2;

  filter = get_3d_farray(xn,yn,zn);
  if (phase==0){  /*** Even "cos" filter ***/
    for(i=0;i<xn;i++){
      x = (float)(i-xc)*sscale;
      x2 = x*x/cs;
      for(j=0;j<yn;j++){
	y = (float)(j-yc)*sscale;
	y2 = y*y/cs;
	frndotr = twopifr * (nx*x + ny*y);
	gs = exp(-x2-y2);
	for(k=0;k<zn;k++){
	  t = (float)(k-tc)*tscale;
	  gt = exp(-t*t/ct);
	  filter[i][j][k] = c*gs*gt*cos(frndotr - twopift*t);
	}
      }
    }
  }else if (phase==1){ /*** Odd "sin" filter ***/
    for(i=0;i<xn;i++){
      x = (float)(i-xc)*sscale;
      x2 = x*x/cs;
      for(j=0;j<yn;j++){
	y = (float)(j-yc)*sscale;
	y2 = y*y/cs;
	frndotr = twopifr * (nx*x + ny*y);
	gs = exp(-x2-y2);
	for(k=0;k<zn;k++){
	  t = (float)(k-tc)*tscale;
	  gt = exp(-t*t/ct);
	  filter[i][j][k] = c*gs*gt*sin(frndotr - twopift*t);
	}
      }
    }
  }else{
    printf("*** GABOR_SPACE_TIME:  No such option.\n");
    exit(0);
  }

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GABOR_SPACE_TIME_TENSOR                         */
/*                                                                           */
/*  Same as "gabor_space_time", but use NumRecInC tensor format.             */
/*                                                                           */
/*****************************************************************************/
float ***gabor_space_time_tensor(xn,yn,zn,sscale,tscale,srx,sry,st,fr,ft,theta,
				 phase,tfilt)
     int xn,yn,zn;         // Size of returned 3D array.  Values must be odd,
                           //    center is (xn-1)/2,...,(zn-1/2)
     float sscale,tscale;  // space and time scaling (deg/pix) and (sec/pix)
     float srx,sry;        // SDs of Gaussian in x and y dimensions (deg)
     float st;             // SD of Gaussian in time (s)
     float fr,ft;          // Spatial and temporal freq. (cyc/deg) and (Hz)
     float theta;          // Orientaiton (spatial) of filter
     float phase;          // 0 for even, 1 for odd
     char *tfilt;          // Type of temporal filter, Gaussian is default
{
  int i,j,k;
  float ***filter,*tdata;
  float c,csx,csy,ct,nx,ny,twopifr,twopift,ph0;
  float x,y,t,x2,y2,gs,gt,frndotr;
  int xc,yc,tc; // centers of dimensions
  float *zt;    // For dumping t-filter
  int dumpflag;

  if (PFLAG) printf("  GABOR_SPACE_TIME_TENSOR\n");

  dumpflag = 0; // WYETH - for dumping t-filter

  if (dumpflag == 1){
    zt = get_zero_farray(zn);
  }

  //printf("srx %f   sry %f\n",srx,sry);

  c = 1.0/(pow(2.0*M_PI,1.5)*srx*sry*st);
  csx = 2.0*srx*srx;
  csy = 2.0*sry*sry;
  ct = 2.0*st*st;
  nx = cos(theta/180.0 * M_PI);
  ny = sin(theta/180.0 * M_PI);
  twopifr = 2.0*M_PI*fr;
  twopift = 2.0*M_PI*ft;
  xc = (xn-1)/2;      // For xn=32, xc = 15, on a 0..32 scale
  yc = (yn-1)/2;
  tc = (zn-1)/2;

  ph0 = M_PI/180.0 * phase;  // Note, set 'phase' = 90 for sin

  // Pre-compute the temporal function
  if (strcmp(tfilt,"cwq_temp")==0){
    // "a" = 2
    // "tau" = st;
    tdata = cwq_temporal(zn,tscale,2,st,0.0,0.0,0.0);
  }else{
    tdata = NULL;
  }


  filter = f3tensor(1,xn,1,yn,1,zn);  // Use tensor for 3D FFT

  for(i=0;i<xn;i++){
    x = (float)(i-xc)*sscale;  // 'x' is spatial location in degrees
    x2 = x*x/csx;
    for(j=0;j<yn;j++){
      y = (float)(j-yc)*sscale;  // 'y' is spatial location in degrees
      y2 = y*y/csy;
      frndotr = twopifr * (nx*x + ny*y);
      gs = exp(-x2-y2);
      for(k=0;k<zn;k++){
	t = (float)(k-tc)*tscale;  // 't' is time in seconds

	// Temporal envelope filter 
	if (strcmp(tfilt,"maxwell")==0){
	  if (t > 0.0)
	    gt = t*t/(st*st) * exp(-t*t/ct);
	  else
	    gt = 0.0;
	}else if (tdata != NULL){
	  gt = tdata[k];
	}else{
	  gt = exp(-t*t/ct);   // Gaussian temporal window by default
	}
 
	// DUMP FOR time plot
	if (dumpflag == 1){
	  if ((i==0)&&(j==0)){
	    zt[k] = gt;
	  }
	}

	filter[i+1][j+1][k+1] = c*gs*gt*cos(frndotr - twopift*t - ph0);
      }
    }
  }

  //write_2d_data("zzz.gabor.xt.2d",filter[xc],0,0,yn,zn,4,2,1,0);
  //printf("WYETH HERE 777\n");
  //exit(0);

  // WYETH - for dumping temporal filter - NOTE '2' scale for tscale
  if (dumpflag == 1)
    append_farray_plot("zzz.GabXYT_temporal.pl","gt",zt,zn,2);

  if (tdata != NULL)
    myfree(tdata);

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                        GABOR_SPACE_TIME_TENSOR_CURVE                      */
/*                                                                           */
/*  Same as "gabor_space_time", but use NumRecInC tensor format.             */
/*                                                                           */
/*****************************************************************************/
float ***gabor_space_time_tensor_curve(xn,yn,zn,sscale,tscale,srx,sry,st,
				       fr,theta,phase,p,tfilt)
     int xn,yn,zn;         // Size of returned 3D array.  Values must be odd,
                           //    center is (xn-1)/2,...,(zn-1/2)
     float sscale,tscale;  // space and time scaling (deg/pix) and (sec/pix)
     float srx,sry;        // SDs of Gaussian in x and y dimensions (deg)
     float st;             // SD of Gaussian in time (s)
     float fr;             // Spatial freq. (cyc/deg)
     float theta;          // Orientaiton (spatial) of filter
     float phase;          // 0 for cos, 90 for sin
     float p;              // Power for curve
     char *tfilt;          // Type of temporal filter, Gaussian is default
{
  int i,j,k;
  float ***filter;
  float c,csx,csy,ct,nx,ny,twopifr,ph0;
  float x,y,t,x2,y2,gs,gt,frndotr;
  int xc,yc,tc; // centers of dimensions
  float *zt;    // For dumping t-filter
  float pht;

  float tt,crv,trad;
  float max_shift_deg,max_shift_cyc,max_shift_rad;

  //**************** WYETH - OLD REMOVE ***** NOT USED  ********* ?
  //**************** WYETH - OLD REMOVE ***** NOT USED  ********* ?
  //**************** WYETH - OLD REMOVE ***** NOT USED  ********* ?
  //**************** WYETH - OLD REMOVE ************************* ?
  //**************** WYETH - OLD REMOVE ************************* ?


  if (PFLAG) printf("  GABOR_SPACE_TIME_TENSOR_CURVE\n");

  pht = 0.0;

  c = 1.0/(pow(2.0*M_PI,1.5)*srx*sry*st);
  csx = 2.0*srx*srx;
  csy = 2.0*sry*sry;
  ct = 2.0*st*st;
  nx = cos(theta/180.0 * M_PI);
  ny = sin(theta/180.0 * M_PI);
  twopifr = 2.0*M_PI*fr;
  //twopift = 2.0*M_PI*ft;
  xc = (xn-1)/2;
  yc = (yn-1)/2;
  tc = (zn-1)/2;

  ph0 = M_PI/180.0 * phase;  // Note, set 'phase' = 90 for sin

  filter = f3tensor(1,xn,1,yn,1,zn);  // Use tensor for 3D FFT

  for(i=0;i<xn;i++){
    x = (float)(i-xc)*sscale;  // 'x' is spatial location in degrees
    x2 = x*x/csx;
    for(j=0;j<yn;j++){
      y = (float)(j-yc)*sscale;  // 'y' is spatial location in degrees
      y2 = y*y/csy;
      frndotr = twopifr * (nx*x + ny*y);
      gs = exp(-x2-y2);
      for(k=0;k<zn;k++){
	t = (float)(k-tc)*tscale;  // 't' is time in seconds

	//   Map the x-axis of our 0..1 box onto the time axis 't' in sec
	//   Note, 'tt' stands for 'transformed time'
	//     tt      t(s)
	//     0  <-->  0  
	//     1  <-->  3 SD of Gaussian, thus 3.0*st
	//
	//   Map the output of the power funtion (0..1) back to spatial coords:
	//    crv      x(deg)
	//     0  <-->  -3 SD in spatial Gaussian
	//     1  <-->  +2 SD in spatial Gaussian
	//
	if ((t > 0.0) && (t < 3.0*st)){

	  max_shift_deg = 15.0*srx;  // Maximum shift (deg)
	  max_shift_cyc = max_shift_deg * fr; // convert to cycles of grating
	  max_shift_rad = 2.0*M_PI * max_shift_cyc; // Max shift (rad)

	  tt = t/(3.0*st);
	  crv =  pow(tt,1.0/p);  // This curve value goes from 0 to 1
	  pht =  crv * max_shift_rad;   // Multiply by total desired radians

	  if ((j == yn/2) && (i == xn/2)){
	    printf("k=%d  pht = %f  (tt = %f  crv = %f  p=%f  maxd=%f)\n",
		   k,pht,tt,crv,p,max_shift_deg);
	  }
	}

	//pht = f(t);


	// Temporal envelope filter 
	if (strcmp(tfilt,"maxwell")==0){
	  if (t > 0.0)
	    gt = t*t/(st*st) * exp(-t*t/ct);
	  else
	    gt = 0.0;
	}else{
	  gt = exp(-t*t/ct);   // Gaussian temporal window by default
	}
 
	//filter[i+1][j+1][k+1] = c*gs*gt*cos(frndotr - twopift*t - ph0);
	filter[i+1][j+1][k+1] = c*gs*gt*cos(frndotr - pht - ph0);
      }
    }
  }

  //write_2d_data("zzz.gabor.xt.2d",filter[xc],0,0,yn,zn,4,2,1,0);
  //printf("WYETH HERE 777\n");
  //exit(0);

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                             KERNU_TSHIFT_POW                              */
/*                                                                           */
/*  Return a float array of length 'tn' that contains the phase shift in     */
/*  radians to be applied to a sinusoid.                                     */
/*                                                                           */
/*****************************************************************************/
float *kernu_tshift_pow(tso,tn,sscale,tscale,sf)
     struct onode *tso;    // Time shift parameter object
     int tn;               // Length of returned array, center is (tn-1/2)
     float sscale;         // space scaling (deg/pix)
     float tscale;         // time scaling (sec/pix)
     float sf;             // spatial frequency (cyc/deg)
{
  int i;
  int tc; // centers of dimensions
  float p,t,t0,t1,tdur,tt,s,*shift,shmax,shmax_cyc,shmax_rad;

  if (PFLAG) printf("  KERNU_TSHIFT_POW\n");

  p     = onode_getpar_flt_exit(tso,"exponent");
  t0    = onode_getpar_flt_exit(tso,"t0");
  t1    = onode_getpar_flt_exit(tso,"t1");
  shmax = onode_getpar_flt_exit(tso,"shift_max");
  // WYETH get units here, and do appropriate conversion??? ***

  if (strcmp("(deg)",onode_getpar_unit(tso,"shift_max"))==0){
    shmax_cyc = shmax * sf;  // deg  x  cyc/deg
  }else if (strcmp("(cyc)",onode_getpar_unit(tso,"shift_max"))==0){
    shmax_cyc = shmax;
  }else
    exit_error("KERNU_TSHIFT_POW","Unknown units on shift_max");

  tc = (tn-1)/2;   // Coordinate of center (origin) of time axis

  tdur = t1 - t0;  // Duration of time window to apply shift

  shmax_rad = 2.0*M_PI * shmax_cyc; // Max shift (rad)

  shift = get_zero_farray(tn);

  //   Map the x-axis of our 0..1 box onto the time axis 't' in sec
  //   Note, 'tt' stands for 'transformed time'
  //     tt      t(s)
  //     0  <-->  t0  
  //     1  <-->  t1
  //
  //   Map the output 's' of the power function (0..1) back to spatial coords:
  //     s       shift(deg)
  //     0  <-->  0
  //     1  <-->  shmax
  //

  for(i=0;i<tn;i++){
    t = (float)(i-tc)*tscale;  // 't' is time in seconds

    if ((t > t0) && (t <= t1)){  // If 't' falls within the shift window
      tt = (t-t0)/(float)tdur;   // 'tt' ranges from 0 to 1 in this window
      s = pow(tt,1.0/p);         // This value goes from 0 to 1
      shift[i] = s * shmax_rad;  // Convert the shift into radians
    }
  }

  return shift;
}
/**************************************-**************************************/
/*                                                                           */
/*                             KERNU_TSHIFT_SIG                              */
/*                                                                           */
/*  Return a float array of length 'tn' that contains the phase shift in     */
/*  radians to be applied to a sinusoid.                                     */
/*                                                                           */
/*****************************************************************************/
float *kernu_tshift_sig(tso,tn,sscale,tscale,sf)
     struct onode *tso;    // Time shift parameter object
     int tn;               // Length of returned array, center is (tn-1/2)
     float sscale;         // space scaling (deg/pix)
     float tscale;         // time scaling (sec/pix)
     float sf;             // spatial frequency (cyc/deg)
{
  int i;
  int tc; // centers of dimensions
  float p,t,t0,t1,tdur,tt,s,*shift,shmax,shmax_cyc,shmax_rad;

  if (PFLAG) printf("  KERNU_TSHIFT_SIG\n");

  p     = onode_getpar_flt_exit(tso,"exponent");
  t0    = onode_getpar_flt_exit(tso,"t0");
  t1    = onode_getpar_flt_exit(tso,"t1");
  shmax = onode_getpar_flt_exit(tso,"shift_max");
  // WYETH get units here, and do appropriate conversion??? ***

  if (strcmp("(deg)",onode_getpar_unit(tso,"shift_max"))==0){
    shmax_cyc = shmax * sf;  // deg  x  cyc/deg
  }else if (strcmp("(cyc)",onode_getpar_unit(tso,"shift_max"))==0){
    shmax_cyc = shmax;
  }else
    exit_error("KERNU_TSHIFT_SIG","Unknown units on shift_max");

  tc = (tn-1)/2;   // Coordinate of center (origin) of time axis

  tdur = t1 - t0;  // Duration of time window to apply shift

  shmax_rad = 2.0*M_PI * shmax_cyc; // Max shift (rad)

  shift = get_zero_farray(tn);

  //   Map the x-axis of our 0..1 box onto the time axis 't' in sec
  //   Note, 'tt' stands for 'transformed time'
  //     tt      t(s)
  //     0  <-->  t0  
  //     1  <-->  t1
  //
  //   Map output 's' of the sigmoid function (0..1) back to spatial coords:
  //     s       shift(deg)
  //     0  <-->  0
  //     1  <-->  shmax
  //
 for(i=0;i<tn;i++){
    t = (float)(i-tc)*tscale;  // 't' is time in seconds
    if ((t > t0) && (t <= t1)){  // If 't' falls within the shift window

      if (1){
	tt = (t-t0)/(float)tdur;   // 'tt' ranges from 0 to 1 in this window
	s = -1.0 + 2.0/(1.0 + exp(-(p * tt)));
      }else{
	tt = -1.0 + (t-t0)/(float)tdur;   // 'tt' ranges from -1 to 0
	s = 1.0 - 2.0/(1.0 + exp(-(p * tt)));
      }
      shift[i] = s * shmax_rad;  // Convert the shift into radians
    }
  }
  return shift;
}
/**************************************-**************************************/
/*                                                                           */
/*                                KERNU_TSHIFT                               */
/*                                                                           */
/*  Call the appropriate function to generate the desired type of phase      */
/*  shift vs. time.                                                          */
/*                                                                           */
/*****************************************************************************/
float *kernu_tshift(tso,tn,sscale,tscale,sf)
     struct onode *tso;    // Time shift parameter object
     int tn;               // Length of returned array, center is (tn-1/2)
     float sscale;         // space scaling (deg/pix)
     float tscale;         // time scaling (sec/pix)
     float sf;             // spatial frequency (cyc/deg)
{
  char *stype;
  float *tshift;

  stype = onode_getpar_chr_ptr_exit(tso,"type");
  if (strcmp(stype,"power")==0){
    tshift = kernu_tshift_pow(tso,tn,sscale,tscale,sf);
  }else if (strcmp(stype, "sigmoid")==0){
    tshift = kernu_tshift_sig(tso,tn,sscale,tscale,sf);
  }else{
     exit_error("KERNU_TSHIFT","Unknown shift type");
  }

  return tshift;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GABOR_SPACE_TIME_TSHIFT                          */
/*                                                                           */
/*  Create a 3D filter that is a Gabor function in space, and (1) is         */
/*  multiplied by the temporal function 'tamp' and (2) the cosine component  */
/*  has an added phase shift vs. time given by 'tshift'.                     */
/*                                                                           */
/*  *** JACOB:  This routine can be used to create a variety of curved       */
/*  filters - you just have to call it with different values stored in       */
/*  'tamp' and 'tshift'.                                                     */
/*                                                                           */
/*****************************************************************************/
float ***gabor_space_time_tshift(xn,yn,zn,sscale,tscale,srx,sry,
				 fr,theta,phase,tamp,tshift)
     int xn,yn,zn;         // Size of returned 3D array.  Values must be odd,
                           //    center is (xn-1)/2,...,(zn-1/2)
     float sscale,tscale;  // space and time scaling (deg/pix) and (sec/pix)
     float srx,sry;        // SDs of Gaussian in x and y dimensions (deg)
     float fr;             // Spatial freq. (cyc/deg)
     float theta;          // Spatial orientation of filter (deg)
     float phase;          // 0 for cos, 90 for sin
     float *tamp;          // [zn] Amplitude of temporal function
     float *tshift;        // [zn] shift at each time point (rad)
{
  int i,j,k;
  int xc,yc;
  float ***filter,c,csx,csy,nx,ny,twopifr,ph0,x,y,x2,y2,gs,frndotr;

  if (PFLAG) printf("  GABOR_SPACE_TIME_TSHIFT\n");

  c = 1.0/(pow(2.0*M_PI,1.5)*srx*sry);
  csx = 2.0*srx*srx;
  csy = 2.0*sry*sry;
  nx = cos(theta/180.0 * M_PI);
  ny = sin(theta/180.0 * M_PI);
  twopifr = 2.0*M_PI*fr;
  xc = (xn-1)/2;  // Central coordinate (origin) along x-dimension
  yc = (yn-1)/2;  // Central coordinate (origin) along y-dimension
  //tc = (zn-1)/2;  // Central coordinate (origin) along t-dimension

  ph0 = M_PI/180.0 * phase;  // Note, set 'phase' = 90 for sin

  filter = f3tensor(1,xn,1,yn,1,zn);  // Use tensor for 3D FFT

  for(i=0;i<xn;i++){
    x = (float)(i-xc)*sscale;  // 'x' is spatial location in degrees
    x2 = x*x/csx;
    for(j=0;j<yn;j++){
      y = (float)(j-yc)*sscale;  // 'y' is spatial location in degrees
      y2 = y*y/csy;
      frndotr = twopifr * (nx*x + ny*y);
      gs = exp(-x2-y2);
      for(k=0;k<zn;k++){
	//t = (float)(k-tc)*tscale;  // 't' is time in seconds
	filter[i+1][j+1][k+1] = c*gs*tamp[k]*cos(frndotr - tshift[k] - ph0);
      }
    }
  }

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                             ROTATE_SPACE_TIME                             */
/*                                                                           */
/*****************************************************************************/
void rotate_space_time(x,y,t,dir,v,rx,ry,rt)
     float x,y;            // spatial position (deg)
     float t;              // time (s)
     float dir;            // Direction in space (deg)
     float v;              // Velocity (deg/s), e.g., TF/SF
     float *rx,*ry,*rt;    // Rotated x,y,t
{
  float theta,phi;
  float u1,u2,v1,v2,a,b,c,ap,bp,cp,xp,yp,zp;
  float xabc1,xabc2,yabc1,yabc2;

  theta = M_PI * dir/180.0;
  phi = atan(1.0/v);

  u1 = cos(theta);   // Vector 'u' points in 'dir'
  u2 = sin(theta);

  v1 = -sin(theta);  // Vector 'v' is orth to 'u'
  v2 = cos(theta);

  a = x*u1 + y* u2;  // proj of (x,y) onto 'u'
  b = x*v1 + y* v2;  // proj of (x,y) onto 'v'
  c = t;             // unchanged

  ap = a * cos(phi) - c * sin(phi);  // Rotate (a,c)
  cp = a * sin(phi) + c * cos(phi);
  bp = b;

  xabc1 = cos(theta);
  xabc2 = -sin(theta);
  yabc1 = sin(theta);
  yabc2 = cos(theta);

  xp = ap * xabc1 + bp * xabc2; // proj of (a',b') on 'x_abc'
  yp = ap * yabc1 + bp * yabc2; // proj of (a',b') on 'y_abc'
  zp = cp;

  *rx = xp;   // Return the rotated coordinates
  *ry = yp;
  *rt = zp;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GABOR_TILT_TENSOR                             */
/*                                                                           */
/*  Built from "gabor_space_time_tensor", but here the Gaussian window       */
/*  tilts in time to match the tilt of the grating.                          */
/*                                                                           */
/*****************************************************************************/
float ***gabor_tilt_tensor(xn,yn,zn,sscale,tscale,sr,st,fr,ft,theta,
			   phase,tfilt,shift_frac)
     int xn,yn,zn;         // must be odd, center is (xn-1)/2,...,(zn-1/2)
     float sscale;         // (deg/pix)
     float tscale;         // (sec/pix)
     float sr;             // SD space (deg)
     float st;             // SD time (s)
     float fr;             // SF (cyc/deg)
     float ft;             // TF (Hz)
     float theta;          // direction (deg)
     int phase;            // 0 for even, 1 for odd
     char *tfilt;          // Type of temporal filter
     float shift_frac;     // Amount of full shift to use, 0.0 - 1.0 - 2.0...
{
  int i,j,k;
  float ***filter;
  float c,cs,ct,nx,ny,twopifr,twopift;
  float x,y,t,x2,y2,gs,gt,frndotr,shx,shy;
  int xc,yc,tc; // centers of dimensions

  if (PFLAG) printf("  GABOR_TILT_TENSOR\n");

  // The spatial shift is a function of 'theta' and the velocity implied
  // by V = TF/SF

  c = 1.0/(pow(2.0*M_PI,1.5)*sr*sr*st);
  cs = 2.0*sr*sr;
  ct = 2.0*st*st;
  nx = cos(theta/180.0 * M_PI);
  ny = sin(theta/180.0 * M_PI);
  twopifr = 2.0*M_PI*fr;
  twopift = 2.0*M_PI*ft;
  xc = (xn-1)/2;
  yc = (yn-1)/2;
  tc = (zn-1)/2;

  /*
  printf("  HERE WYETH   sr = %f\n",sr);
  printf("  HERE WYETH   st = %f\n",st);
  printf("  HERE WYETH   fr = %f\n",fr);
  printf("  HERE WYETH   ft = %f\n",ft);
  printf("  HERE WYETH   tscale = %f\n",tscale);
  printf("  HERE WYETH   theta = %f\n",theta);
  printf("               phase = %d\n\n",phase);
  */


  // Shift with time
  shx = cos(M_PI*theta/180.0) * ft/fr;  // Shift in x w/ time  (deg/s)
  shy = sin(M_PI*theta/180.0) * ft/fr;  // Shift in y w/ time  (deg/s)

  shx *= shift_frac;  // Allows us to control the amount of shift
  shy *= shift_frac;  // Allows us to control the amount of shift

  //printf("    shx,shy  %f %f\n",shx,shy);

  filter = f3tensor(1,xn,1,yn,1,zn);  // Use tensor for 3D FFT

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	t = (float)(k-tc)*tscale;  // (sec)
	x = (float)(i-xc)*sscale;  // (deg)
	y = (float)(j-yc)*sscale;  // (deg)
	
	frndotr = twopifr * (nx*x + ny*y);
	
	// Do tilt here, after 'frndotr', so this affects the window
	// and *NOT* the grating within.
	x -= t*shx;                // Tilt:  add shift w.r.t. time here
	y -= t*shy;                // Tilt:  add shift w.r.t. time here
	
	x2 = x*x/cs;
	y2 = y*y/cs;
	gs = exp(-x2-y2);
	
	// Temporal envelope filter 
	if (strcmp(tfilt,"maxwell")==0){
	  if (t > 0.0)
	    gt = t*t/(st*st) * exp(-t*t/ct);
	  else
	    gt = 0.0;
	}else{
	  gt = exp(-t*t/ct);   // Gaussian by default
	}
	
	if (phase == 0)       // Even 'cos' filter
	  filter[i+1][j+1][k+1] = c*gs*gt*cos(frndotr - twopift*t);
	else if (phase == 1)  // Odd 'sine' filter
	  filter[i+1][j+1][k+1] = c*gs*gt*sin(frndotr - twopift*t);
	else
	  exit_error("GABOR_TILT_TENSOR","No such 'phase'");
      }
    }
  }

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GABOR_TILT2_TENSOR                            */
/*                                                                           */
/*  Like 'gabor_tilt_tensor' but here a matrix rotation is used.             */
/*                                                                           */
/*****************************************************************************/
float ***gabor_tilt2_tensor(xn,yn,zn,sscale,tscale,sr,st,fr,ft,theta,
			   phase,tfilt)
     int xn,yn,zn;         // must be odd, center is (xn-1)/2,...,(zn-1/2)
     float sscale;         // (deg/pix)
     float tscale;         // (sec/pix)
     float sr;             // SD space (deg)
     float st;             // SD time (s)
     float fr;             // SF (cyc/deg)
     float ft;             // TF (Hz)
     float theta;          // direction (deg)
     int phase;            // 0 for even, 1 for odd
     char *tfilt;          // Type of temporal filter
{
  int i,j,k;
  float ***filter;
  float c,cs,ct,nx,ny,twopifr,twopift;
  float x,y,t,x2,y2,gs,gt,frndotr,shx,shy;
  int xc,yc,tc; // centers of dimensions
  float xraw,yraw,traw,arot,ang,vx,vy,**m,*v,*r,xf; // Rotation

  if (PFLAG) printf("  GABOR_TILT2_TENSOR\n");

  if ((theta != 0.0) && (theta != 180.0))
    exit_error("GABOR_TILT2_TENSOR","Only works for 0/180");

  // The spatial shift is a function of 'theta' and the velocity implied
  // by V = TF/SF

  c = 1.0/(pow(2.0*M_PI,1.5)*sr*sr*st);
  cs = 2.0*sr*sr;
  ct = 2.0*st*st;
  nx = cos(theta/180.0 * M_PI);
  ny = sin(theta/180.0 * M_PI);
  twopifr = 2.0*M_PI*fr;
  twopift = 2.0*M_PI*ft;
  xc = (xn-1)/2;
  yc = (yn-1)/2;
  tc = (zn-1)/2;

  /*
  printf("  HERE t2 WYETH   sr = %f\n",sr);
  printf("  HERE WYETH   st = %f\n",st);
  printf("  HERE WYETH   fr = %f\n",fr);
  printf("  HERE WYETH   ft = %f\n",ft);
  printf("  HERE WYETH   tscale = %f\n",tscale);
  printf("  HERE WYETH   theta = %f\n",theta);
  printf("               phase = %d\n\n",phase);
  */


  // Rotation matrix
  arot = 180.0/M_PI * atan2(ft,fr);     // Amount of rotation (deg)
  ang  = theta + 90.0;  // Theta in XY defining axis of rotation
  vx = cos(M_PI/180.0 * ang);
  vy = sin(M_PI/180.0 * ang);
  m = rotation_matrix_3d(-arot,vx,vy,0.0);
  v = get_farray(3);

  xf = 1.0/sin(M_PI/180.0 * (90.0-arot));  // x-factor, stretch due to tilt
  
  //printf("******vx,y %f %f ************ arot %f  ang %f\n",vx,vy,arot,ang);

  filter = f3tensor(1,xn,1,yn,1,zn);  // Use tensor for 3D FFT

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	traw = (float)(k-tc)*tscale;  // (sec)
	xraw = (float)(i-xc)*sscale;  // (deg)
	yraw = (float)(j-yc)*sscale;  // (deg)
	
	// Rotation Matrix method
	v[0] = xraw;
	v[1] = yraw;
	v[2] = traw;
	r = matrix_times_vector(m,v,3,3);
	x = r[0];
	y = r[1];
	t = r[2];
	myfree(r);

	//  WYETH - THIS IS A HACK - Doesn't work for theta not 0/180????
	//  WYETH - THIS IS A HACK - Doesn't work for theta not 0/180????
	//  WYETH - THIS IS A HACK - Doesn't work for theta not 0/180????
	x *= xf;  // Undo the stretching in x
	t /= xf;  // Undo the stretching in t

	/*
	if ((xraw >= 1.5) && (yraw == 0.0)){
	  if (traw > 0.12)
	    printf("%f %f %f  -->  %f %f %f\n",xraw,yraw,traw,x,y,t);
	}
	*/

	
	frndotr = twopifr * (nx*xraw + ny*yraw);  // No rotation for sin

	x2 = x*x/cs;  // Use rotated coords here
	y2 = y*y/cs;
	gs = exp(-x2-y2);

	// Temporal envelope filter, use rotated coords here
	if (strcmp(tfilt,"maxwell")==0){
	  if (t > 0.0)
	    gt = t*t/(st*st) * exp(-t*t/ct);
	  else
	    gt = 0.0;
	}else{
	  gt = exp(-t*t/ct);   // Gaussian by default
	}

	// Rectangular window
	gs = gt = 0.0;
	if (fabs(x) < 0.4)
	  gs = 1.0;
	if (fabs(t) < 0.080)
	  gt = 1.0;
	
	if (phase == 0)       // Even 'cos' filter
	  filter[i+1][j+1][k+1] = c*gs*gt*cos(frndotr - twopift*traw);
	else if (phase == 1)  // Odd 'sine' filter
	  filter[i+1][j+1][k+1] = c*gs*gt*sin(frndotr - twopift*traw);
	else
	  exit_error("GABOR_TILT_TENSOR","No such 'phase'");
      }
    }
  }

  free_2d_farray(m,3);
  myfree(v);

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GABOR_TILT3_TENSOR                            */
/*                                                                           */
/*  Like 'gabor_tilt_tensor' but here a matrix rotation is used.             */
/*                                                                           */
/*****************************************************************************/
float ***gabor_tilt3_tensor(xn,yn,zn,sscale,tscale,sr,st,fr,ft,theta,
			   phase,tfilt)
     int xn,yn,zn;         // must be odd, center is (xn-1)/2,...,(zn-1/2)
     float sscale;         // (deg/pix)
     float tscale;         // (sec/pix)
     float sr;             // SD space (deg)
     float st;             // SD time (s)
     float fr;             // SF (cyc/deg)
     float ft;             // TF (Hz)
     float theta;          // direction (deg)
     int phase;            // 0 for even, 1 for odd
     char *tfilt;          // Type of temporal filter
{
  int i,j,k;
  float ***filter;
  float c,cs,ct,nx,ny,twopifr,twopift;
  float x,y,t,x2,y2,gs,gt,frndotr;
  int xc,yc,tc;                          // centers of dimensions
  float xraw,yraw,traw;                  // Rotation
  
  float vel;
  float rota,rotb,rotc,rotd,xu,tu,rxu,rtu;

  if (PFLAG) printf("  GABOR_TILT3_TENSOR\n");

  if ((theta != 0.0) && (theta != 180.0)) // Because it doesn't rotate in y
    exit_error("GABOR_TILT3_TENSOR","Only works for 0 or 180");
  
  c = 1.0/(pow(2.0*M_PI,1.5)*sr*sr*st);
  cs = 2.0*sr*sr;
  ct = 2.0*st*st;
  nx = cos(theta/180.0 * M_PI);
  ny = sin(theta/180.0 * M_PI);
  twopifr = 2.0*M_PI*fr;
  twopift = 2.0*M_PI*ft;
  xc = (xn-1)/2;
  yc = (yn-1)/2;
  tc = (zn-1)/2;

  /*
  printf("  HERE WYETH   sr = %f\n",sr);
  printf("  HERE WYETH   st = %f\n",st);
  printf("  HERE WYETH   fr = %f\n",fr);
  printf("  HERE WYETH   ft = %f\n",ft);
  printf("  HERE WYETH   tscale = %f\n",tscale);
  printf("  HERE WYETH   theta = %f\n",theta);
  printf("               phase = %d\n\n",phase);
  */


  filter = f3tensor(1,xn,1,yn,1,zn);  // Use tensor for 3D FFT

  //  Use the velocity to square up the space
  //  and then the rotation is always a 45 degree rotation

  vel = ft/fr;  // velocity, deg/sec
  
  if (theta == 0.0){
    rota =  cos(M_PI/4.0); rotb = sin(M_PI/4.0); // rotation matrix
    rotc = -sin(M_PI/4.0); rotd = cos(M_PI/4.0);
  }else{
    rota =  cos(M_PI/-4.0); rotb = sin(M_PI/-4.0); // rotation matrix
    rotc = -sin(M_PI/-4.0); rotd = cos(M_PI/-4.0);
  }

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	traw = (float)(k-tc)*tscale;  // (sec)
	xraw = (float)(i-xc)*sscale;  // (deg)
	yraw = (float)(j-yc)*sscale;  // (deg)

	xu = xraw / vel;  // These are now in a square space
	tu = traw;        // so we want to rotate them by 45 deg

	rxu = rota*xu + rotc*tu;
	rtu = rotb*xu + rotd*tu;

	x = vel * rxu;   // Convert back to 'deg'
	t = rtu;         // still in 'sec'
	y = yraw;        // No change in y-value

	// At this point, the rotated values are (x,y,t)

	
	frndotr = twopifr * (nx*xraw + ny*yraw);  // No rotation for sin

	// Spatial Gaussian - Use rotated coords here
	x2 = x*x/cs;
	y2 = y*y/cs;
	gs = exp(-x2-y2);

	// Temporal envelope filter, use rotated coords here
	if (strcmp(tfilt,"maxwell")==0){
	  if (t > 0.0)
	    gt = t*t/(st*st) * exp(-t*t/ct);
	  else
	    gt = 0.0;
	}else{
	  gt = exp(-t*t/ct);   // Gaussian by default
	}

	// Replace Gaussian window w/ Rectangular window - for testing
	/*
	gs = gt = 0.0;
	if (fabs(x) < 0.4)
	  gs = 1.0;
	if (fabs(t) < 0.080)
	  gt = 1.0;
	*/
	
	if (phase == 0)       // Even 'cos' filter
	  filter[i+1][j+1][k+1] = c*gs*gt*cos(frndotr - twopift*traw);
	else if (phase == 1)  // Odd 'sine' filter
	  filter[i+1][j+1][k+1] = c*gs*gt*sin(frndotr - twopift*traw);
	else
	  exit_error("GABOR_TILT3_TENSOR","No such 'phase'");
      }
    }
  }

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                        GABOR_SPACE_TIME_SEP_TENSOR                        */
/*                                                                           */
/*  2D spatial Gabor function with a choice of separable temporal filters.   */
/*                                                                           */
/*****************************************************************************/
float ***gabor_space_time_sep_tensor(ppl,o,xn,yn,tn,sscale,tscale,phase0)
     struct param_pair_list *ppl;
     struct onode *o;
     int xn,yn,tn;         // must be odd, center is (xn-1)/2,...,(zn-1/2)
     float sscale,tscale;  // space and time scaling
     float phase0;         // spatial phase offset (deg)
{
  int i,j,k;
  float xc,yc; // centers of dimensions
  int dumpflag;
  float ***filter,*tt,sf,ssd,ssdo,ssdp,sfpix,sdo,sdp,**ss,sv,t0,ph,phase;
  float c,cs,nx,ny,twopifr,x,y,x2,y2,frndotr,gs,theta;
  char *tfilt;

  if (PFLAG) printf("  GABOR_SPACE_TIME_SEP_TENSOR\n");

  tt = (float *)myalloc(tn*sizeof(float));

  t0 = (float)((tn-1)/2);  // This is an integer value, stored as float
			
  if (ppl != NULL){
    sf       = paramfile_get_float_param_or_exit(ppl,"gabor_sf");
    ssdo     = paramfile_get_float_param_or_exit(ppl,"gabor_ssd_orth");
    ssdp     = paramfile_get_float_param_or_exit(ppl,"gabor_ssd_par");
    theta    = paramfile_get_float_param_or_exit(ppl,"gabor_theta");
    phase    = paramfile_get_float_param_or_exit(ppl,"gabor_phase");
    dumpflag = paramfile_get_int_param_default(ppl,"gabor_tfilt_dump",0);
    tfilt    = paramfile_get_char_param_or_exit(ppl,"gabor_tfilt");
  }else{
    sf       = onode_getpar_flt_exit(o,"gabor_sf");
    ssdo     = onode_getpar_flt_exit(o,"gabor_ssd_orth");
    ssdp     = onode_getpar_flt_exit(o,"gabor_ssd_par");
    theta    = onode_getpar_flt_exit(o,"gabor_theta");
    phase    = onode_getpar_flt_exit(o,"gabor_phase");
    dumpflag = onode_getpar_int_dflt(o,"gabor_tfilt_dump",0);
    tfilt    = onode_getpar_chr_exit(o,"gabor_tfilt");
  }

  if (strcmp(tfilt,"gaussian")==0){
    float tsd,sd;

    if (ppl != NULL)
      tsd   = paramfile_get_float_param_or_exit(ppl,"gabor_tsd");
    else
      tsd   = onode_getpar_flt_exit(o,"gabor_tsd");

    sd = tsd / tscale;
    tt = gaussian_one_farray(t0,sd,1.0,tn);
  }else if (strcmp(tfilt,"gabor")==0){
    float tsd,sd,tf,tph0;
    
    if (ppl != NULL){
      tsd   = paramfile_get_float_param_or_exit(ppl,"gabor_tsd");
      tf    = paramfile_get_float_param_or_exit(ppl,"gabor_tf");
      tph0  = paramfile_get_float_param_or_exit(ppl,"gabor_tphase");
    }else{
      tsd   = onode_getpar_flt_exit(o,"gabor_tsd");
      tf    = onode_getpar_flt_exit(o,"gabor_tf");
      tph0  = onode_getpar_flt_exit(o,"gabor_tphase");
    }
    
    sd = tsd / tscale;
    tt = gaussian_one_farray(t0,sd,1.0,tn);
    tf *= 2.0*M_PI * tscale;
    tph0 *= M_PI/180.0;
    
    for(i=0;i<tn;i++){
      x = (float)(i-t0);
      tt[i] *= cos(tf*x - tph0);
    }
    
  }else if (strcmp(tfilt,"maxwell_diff")==0){
    float tmax,s,s2,tmax2,tamp2,*tt2;
    
    if (ppl != NULL){
      tmax = paramfile_get_float_param_or_exit(ppl,"gabor_maxwell_tmax");
      tmax2 = paramfile_get_float_param_or_exit(ppl,"gabor_maxwell_tmax_2");
      tamp2 = paramfile_get_float_param_or_exit(ppl,"gabor_maxwell_tamp_2");
    }else{
      tmax = onode_getpar_flt_exit(o,"gabor_maxwell_tmax");
      tmax2 = onode_getpar_flt_exit(o,"gabor_maxwell_tmax_2");
      tamp2 = onode_getpar_flt_exit(o,"gabor_maxwell_tamp_2");
    }
    s  = tmax  / sqrt(2.0) / tscale;
    s2 = tmax2 / sqrt(2.0) / tscale;
    
    tt  = maxwell_farray(t0,s ,1.0,tn);
    tt2 = maxwell_farray(t0,s2,1.0,tn);
    for(i=0;i<tn;i++)
      tt[i] = tt[i] - tamp2 * tt2[i];
    myfree(tt2);
    
    //write_farray_plot("zzz.t.pl",tt,tn);
  }else if (strcmp(tfilt,"maxwell")==0){
    float tmax,s;

    if (ppl != NULL)
      tmax = paramfile_get_float_param_or_exit(ppl,"gabor_maxwell_tmax");
    else
      tmax = onode_getpar_flt_exit(o,"gabor_maxwell_tmax");

    s = tmax / sqrt(2.0) / tscale;
    tt = maxwell_farray(t0,s,1.0,tn);
  }else
    exit_error("GABOR_SPACE_TIME_SEP_TENSOR","Unknown 'gabor_tfilt'");


  // Dump the temporal filter
  if (dumpflag == 1){
    float tsum;
    char tname[SLEN];
    
    tsum = sum_farray(tt,tn,0,tn);
    
    printf("GABOR_SPACE_TIME_SEP_TENSOR  tfilter sum =  %f\n",tsum);
    
    sprintf(tname,"%s_tscale_units",tfilt);
    append_farray_plot("zzz.tfilt.pl",tname,tt,tn,1);
  }

  myfree(tfilt);

  //  NOTE - temporal filter may be biphasic, thus has arbitrary area
  //  so we will NOT normalize it

  xc = (xn-1)/2;
  yc = (yn-1)/2;
  sdo = ssdo / sscale;  // Convert to pixels
  sdp = ssdp / sscale;  // Convert to pixels
  sfpix = sf * sscale;  // Convert to cyc/pix
  ss = gabor_2d_space_raw(xn,yn,xc,yc,sdo,sdp,sfpix,theta,phase+phase0);
  //write_2d_data("zzz.gabor.2d",ss,0,0,xn,yn,4,2,1,0);

  filter = f3tensor(1,xn,1,yn,1,tn);  // Use tensor for 3D FFT
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      sv = ss[i][j];
      for(k=0;k<tn;k++){
	filter[i+1][j+1][k+1] = sv * tt[k];
      }
    }
  }

  myfree(tt);
  free_2d_farray(ss,xn);

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GAUSS_SPACE_EXP_TIME_TENSOR                       */
/*                                                                           */
/*  A 2D spatial Gaussian and a 1D decalying exponential.                    */
/*                                                                           */
/*     *** TIME ZERO IS AT BEGINNING OF ARRAY, for now                       */
/*                                                                           */
/*****************************************************************************/
float ***gauss_space_exp_time_tensor(mylogf,xn,yn,zn,sscale,tscale,sig,tau)
     char *mylogf;
     int xn,yn,zn; // must be odd numbers, center is (xn-1)/2,...,(zn-1/2)
     float sscale,tscale; // space and time scaling
     float sig;  // Gaussian SD
     float tau;  // Exponential decay
{
  int i,j,k;
  int xc,yc,tc;
  float ***filter,cs,x,y,t,x2,y2,gs,gt;
  double dsum;
  char tstr[SLEN];

  mylog(mylogf,"  GAUSS_SPACE_EXP_TIME_TENSOR\n");

  sprintf(tstr,"    xn,yn,zn  %d %d %d\n",xn,yn,zn);
  mylog(mylogf,tstr);
  sprintf(tstr,"    sscale %f   tscale %f\n",sscale,tscale);
  mylog(mylogf,tstr);
  sprintf(tstr,"    sig %f   tau %f\n",sig,tau);
  mylog(mylogf,tstr);

  cs = 2.0*sig*sig;
  xc = (xn-1)/2;
  yc = (yn-1)/2;
  //tc = (zn-1)/2;

  tc = 0; // WYETH - time zero is at beginning, for MOD_DOG_UTIL

  filter = f3tensor(1,xn,1,yn,1,zn); // useful for 3D FFT

  dsum = 0.0;
  for(i=0;i<xn;i++){
    x = (float)(i-xc)*sscale;
    x2 = x*x/cs;
    for(j=0;j<yn;j++){
      y = (float)(j-yc)*sscale;
      y2 = y*y/cs;
      gs = exp(-x2-y2);
      for(k=0;k<zn;k++){
	t = (float)(k-tc)*tscale;
	if (t >= 0.0)
	  gt = exp(-t/tau);
	else
	  gt = 0.0;

	filter[i+1][j+1][k+1] = gs*gt;
	dsum += gs*gt;
	
	/***if (isinf(dsum)){
	    sprintf(tstr,"    gs %f   gt %f\n",gs,gt);
	    mylog(mylogf,tstr);
	    }***/
      }
    }
  }

  sprintf(tstr,"    Sum of filter before norm: %lf\n",dsum);
  mylog(mylogf,tstr);

  for(i=0;i<xn;i++){ // Make area be 1.0
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	filter[i+1][j+1][k+1] /= dsum;
      }
    }
  }

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                               EXP_TIME_TENSOR                             */
/*                                                                           */
/*  1D decaying exponential.                                                 */
/*                                                                           */
/*     *** TIME ZERO IS AT BEGINNING OF ARRAY, for now                       */
/*                                                                           */
/*****************************************************************************/
float ***exp_time_tensor(mylogf,xn,yn,zn,sscale,tscale,tau)
     char *mylogf;
     int xn,yn,zn; // must be odd numbers, center is (xn-1)/2,...,(zn-1/2)
     float sscale,tscale; // space and time scaling
     float tau;  // Exponential decay
{
  int i,j,k;
  int xc,yc,tc;
  float ***filter,cs,x,y,t,gt;
  double dsum;
  char tstr[SLEN];

  mylog(mylogf,"  EXP_TIME_TENSOR\n");

  sprintf(tstr,"    xn,yn,zn  %d %d %d\n",xn,yn,zn);
  mylog(mylogf,tstr);
  sprintf(tstr,"    sscale %f   tscale %f\n",sscale,tscale);
  mylog(mylogf,tstr);
  sprintf(tstr,"    tau %f\n",tau);
  mylog(mylogf,tstr);

  xc = (xn-1)/2;
  yc = (yn-1)/2;

  tc = (zn-1)/2; // WYETH - this works for Reichardt model

  filter = f3tensor(1,xn,1,yn,1,zn); // useful for 3D FFT

  // WYETH - this could be much more efficient
  dsum = 0.0;
  for(i=0;i<xn;i++){
    x = (float)(i-xc)*sscale;       // not used
    for(j=0;j<yn;j++){
      y = (float)(j-yc)*sscale;     // not used
      for(k=0;k<zn;k++){
	t = (float)(k-tc)*tscale;
	//printf("t=%f\n",t);
	if (t >= 0.0)
	  gt = exp(-t/tau);
	else
	  gt = 0.0;

	if ((i == xc) && (j == yc)){
	  filter[i+1][j+1][k+1] = gt;
	  dsum += gt;
	  //  if (dsum > 0.0)
	  //  printf("dsum = %f\n",dsum);
	}else{
	  filter[i+1][j+1][k+1] = 0.0;
	}
	
      }
    }
  }

  sprintf(tstr,"    Sum of filter before norm: %lf\n",dsum);
  mylog(mylogf,tstr);

  for(i=0;i<xn;i++){ // Make area be 1.0
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	filter[i+1][j+1][k+1] /= dsum;
      }
    }
  }

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                                DELTA_TENSOR                               */
/*                                                                           */
/*  Delta-function.                                                          */
/*                                                                           */
/*****************************************************************************/
float ***delta_tensor(mylogf,xn,yn,zn,sscale,tscale)
     char *mylogf;
     int xn,yn,zn;
     float sscale,tscale; // space and time scaling
{
  int i,j,k;
  int xc,yc,tc;
  float ***filter;
  char tstr[SLEN];

  mylog(mylogf,"  EXP_TIME_TENSOR\n");

  sprintf(tstr,"    xn,yn,zn  %d %d %d\n",xn,yn,zn);
  mylog(mylogf,tstr);
  sprintf(tstr,"    sscale %f   tscale %f\n",sscale,tscale);
  mylog(mylogf,tstr);

  // WYETH BUGFIX 14th Oct, 2008 - changed to work w/ mod_dog_util
  //   so that the stim and response are aligned
  //xc = 1+(xn-1)/2;
  //yc = 1+(yn-1)/2;

  xc = xn/2;
  yc = yn/2;
  tc = (zn-1)/2; // WYETH - this works for Reichardt model

  filter = f3tensor(1,xn,1,yn,1,zn); // useful for 3D FFT
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	filter[i+1][j+1][k+1] = 0.0;
      }
    }
  }

  filter[xc+1][yc+1][tc+1] = 1.0;

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_OPPONENT_RANDOM                           */
/*                                                                           */
/*****************************************************************************/
void get_opponent_random(xn,yn,tn,sscale,tscale,sr,st,nsd_r,nsd_t,seed,
			 fpe,fpo,fne,fno)
     int xn,yn,tn;
     float sscale,tscale;
     float sr,st;       // Spatial and temporal SDs for Gaussian envelop
     float nsd_r,nsd_t; // smoothing of noise SDs in space and time
     int seed;          // randomization seed
     float ****fpe,****fpo,****fne,****fno;
{
  int i,j,k;
  int irev,seed2;

  if (seed > 0)
    seed = -seed;

  // WYETH - for making fpo be different
  seed2 = seed - 113;

  *fpe = random_filter(xn,yn,tn,sscale,tscale,sr,st,nsd_r,nsd_t,&seed);
  *fpo = random_filter(xn,yn,tn,sscale,tscale,sr,st,nsd_r,nsd_t,&seed2);
  //
  // Make 'fpo' be a copy of 'fpe'
  //
  /*
  *fpo = f3tensor(1,xn,1,yn,1,tn);
  for(i=1;i<=xn;i++)
    for(j=1;j<=yn;j++)
      for(k=1;k<=tn;k++)
	(*fpo)[i][j][k] = (*fpe)[i][j][k];
  */

  //
  //  Make mirror images of "pref" filters
  //
  *fne = f3tensor(1,xn,1,yn,1,tn);
  *fno = f3tensor(1,xn,1,yn,1,tn);
  for(i=1;i<=xn;i++){
    // WYETH - THIS SEEMS LIKE BAD CENTERING for EVEN xn:  (xn-1)/2
    irev = (xn-1) - i;  // Mirror around center, xc = (xn-1)/2
    if (irev < 1)
      irev = 1;
    else if (irev > xn)
      irev = xn;

      for(j=1;j<=yn;j++)
	for(k=1;k<=tn;k++){
	  (*fne)[i][j][k] = (*fpe)[irev][j][k];
	  (*fno)[i][j][k] = (*fpo)[irev][j][k];
	}
  }

  printf("=====> fpe SUM:  %f\n",sum_3d_farray(*fpe,1,xn,1,yn,1,tn));
  printf("=====> fpo SUM:  %f\n",sum_3d_farray(*fpo,1,xn,1,yn,1,tn));
  printf("=====> fne SUM:  %f\n",sum_3d_farray(*fne,1,xn,1,yn,1,tn));
  printf("=====> fno SUM:  %f\n",sum_3d_farray(*fno,1,xn,1,yn,1,tn));

}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_OPPONENT_RANDOM_X                          */
/*                                                                           */
/*  Extra filter - PE, PO, PX, NE, NO, NX                                    */
/*                                                                           */
/*****************************************************************************/
void get_opponent_random_x(xn,yn,tn,sscale,tscale,sr,st,nsd_r,nsd_t,seed,
			   fpe,fpo,fpx,fne,fno,fnx)
     int xn,yn,tn;
     float sscale,tscale;
     float sr,st;       // Spatial and temporal SDs for Gaussian envelop
     float nsd_r,nsd_t; // smoothing of noise SDs in space and time
     int seed;          // randomization seed
     float ****fpe,****fpo,****fpx,****fne,****fno,****fnx;
{
  int i,j,k;
  int irev,seed2,seed3;

  if (seed > 0)
    seed = -seed;

  // WYETH - for making fpo be different
  seed2 = seed - 2113;
  seed3 = seed - 93653;

  printf("seed = %d\n",seed);
  printf("seed2 = %d\n",seed2);
  printf("seed3 = %d\n",seed3);

  *fpe = random_filter(xn,yn,tn,sscale,tscale,sr,st,nsd_r,nsd_t,&seed);
  *fpo = random_filter(xn,yn,tn,sscale,tscale,sr,st,nsd_r,nsd_t,&seed2);
  *fpx = random_filter(xn,yn,tn,sscale,tscale,sr,st,nsd_r,nsd_t,&seed3);

  //
  //  Make mirror images of "pref" filters
  //
  *fne = f3tensor(1,xn,1,yn,1,tn);
  *fno = f3tensor(1,xn,1,yn,1,tn);
  *fnx = f3tensor(1,xn,1,yn,1,tn);
  for(i=1;i<=xn;i++){
    // WYETH - THIS SEEMS LIKE BAD CENTERING for EVEN xn:  (xn-1)/2
    irev = (xn-1) - i;  // Mirror around center, xc = (xn-1)/2
    if (irev < 1)
      irev = 1;
    else if (irev > xn)
      irev = xn;

      for(j=1;j<=yn;j++)
	for(k=1;k<=tn;k++){
	  (*fne)[i][j][k] = (*fpe)[irev][j][k];
	  (*fno)[i][j][k] = (*fpo)[irev][j][k];
	  (*fnx)[i][j][k] = (*fpx)[irev][j][k];
	}
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_QUAD_OPPONENT_GABOR                         */
/*                                                                           */
/*  Added 'ph0' parameter on Oct 22, 2016.                                   */
/*                                                                           */
/*****************************************************************************/
void get_quad_opponent_gabor(xn,yn,tn,sscale,tscale,sr,st,fr,ft,theta,
			     disp_off,disp_adj,tfilt,
			     fpe,fpo,fne,fno)
     int xn,yn,tn;
     float sscale,tscale,sr,st,fr,ft,theta;
     float disp_off;   // Disparity phase offset to add to all filter phases
     int   disp_adj;   // 1 - adjust disparity offset for opposite direction
     char *tfilt;
     float ****fpe,****fpo,****fne,****fno;
{
  float ntheta,pph,nph,twrap;

  ntheta = theta + 180.0;

  pph = nph = 0.0;

  //
  //  WYETH NOTE:  Could replace this with two calls to:
  //  MISC_UTIL.C:  disparity_phase_shift_adjust(aflag,theta,disp_off)
  //

  if (disp_adj >= 1){
    //
    //  If we want to make these filters be offset for disparity (from some
    //  other set of filters, not among these), then we have to shift 
    //  oppositely directed filters in opposite directions.
    //
    if (disp_off != 0.0){  // If we want a phase offset
      twrap = theta;
      wrap_float(&twrap,360.0);  // Make sure 'twrap' is in [ 0,360.0 )

      if (disp_adj == 2){
	if ((twrap > 90.0) && (twrap <= 270.0)){
	  pph = -disp_off;
	  nph =  disp_off;
	}else{
	  pph =  disp_off;  // This range includes direction (twrap) = 0
	  nph = -disp_off;
	}
      }else{  // if (disp_adj == 1){
	pph = disp_off * cos( theta * M_PI/180.0);
	nph = disp_off * cos(ntheta * M_PI/180.0);
	//printf("  theta:  %f   pph  %f\n",theta,pph);
	//printf("   anti:  %f   nph  %f\n",ntheta,nph);
      }
    }
  }

  //
  //  *** 'sr' is used twice below, for sd_orth and sd_par
  //
  *fpe = gabor_space_time_tensor(xn,yn,tn,sscale,tscale,sr,sr,st,fr,ft,theta,
				 pph     ,tfilt);
  *fpo = gabor_space_time_tensor(xn,yn,tn,sscale,tscale,sr,sr,st,fr,ft,theta,
				 pph+90.0,tfilt);
  *fne = gabor_space_time_tensor(xn,yn,tn,sscale,tscale,sr,sr,st,fr,ft,ntheta,
				 nph     ,tfilt);
  *fno = gabor_space_time_tensor(xn,yn,tn,sscale,tscale,sr,sr,st,fr,ft,ntheta,
				 nph+90.0,tfilt);
}
/**************************************-**************************************/
/*                                                                           */
/*                        GET_QUAD_OPPONENT_GABOR_CURVE                      */
/*                                                                           */
/*****************************************************************************/
void get_quad_opponent_gabor_curve(xn,yn,tn,sscale,tscale,sr,st,fr,theta,
				   pwr,tfilt,fpe,fpo,fne,fno)
     int xn,yn,tn;
     float sscale,tscale,sr,st,fr,theta,pwr;
     char *tfilt;
     float ****fpe,****fpo,****fne,****fno;
{
  float ntheta;

  printf("pwr = %f\n",pwr);

  ntheta = theta + 180.0;

  //
  //  *** 'sr' is used twice below, for sd_orth and sd_par
  //

  *fpe = gabor_space_time_tensor_curve(xn,yn,tn,sscale,tscale,sr,sr,st,fr,
				       theta,0.0,pwr,tfilt);
  *fpo = gabor_space_time_tensor_curve(xn,yn,tn,sscale,tscale,sr,sr,st,fr,
				       theta,90.0,pwr,tfilt);
  *fne = gabor_space_time_tensor_curve(xn,yn,tn,sscale,tscale,sr,sr,st,fr,
				       ntheta,0.0,pwr,tfilt);
  *fno = gabor_space_time_tensor_curve(xn,yn,tn,sscale,tscale,sr,sr,st,fr,
				       ntheta,90.0,pwr,tfilt);
}
/**************************************-**************************************/
/*                                                                           */
/*                       GET_QUAD_OPPONENT_GABOR_TSHIFT                      */
/*                                                                           */
/*****************************************************************************/
void get_quad_opponent_gabor_tshift(mylogf,fo,xn,yn,tn,sscale,tscale,
				    fpe,fpo,fne,fno)
     char *mylogf;
     struct onode *fo;     // filter parameter object
     int xn,yn,tn;         // dimensions of filters to return
     float sscale,tscale;  // deg/pix and sec/frame
     float ****fpe,****fpo,****fne,****fno; // [xn][yn][tn] returned filters
{
  float sdx,sdy,sf,theta,atheta,*tshift,*tamp;
  char *dumpfile;
  struct onode *tso;   // Shift vs. time parameter object
  struct onode *tfo;   // Temporal filter onode

  //printf("  GET_QUAD_OPPONENT_GABOR_TSHIFT\n");

  //
  //  Read parameters from the 3D <filter> object
  //
  theta     = onode_getpar_flt_exit(fo,"direction");
  sdx = sdy = onode_getpar_flt_exit(fo,"ssd");
  sf        = onode_getpar_flt_exit(fo,"sf");

  atheta = theta + 180.0;  // Set the anti-preferred direction

  //
  //  Compute the phase shift values:  'tshift[tn]'
  //
  tso = onode_child_get_unique(fo,"shift_v_t");
  tshift = kernu_tshift(tso,tn,sscale,tscale,sf);

  //
  //  Optionally write a plot of the shift function
  //
  dumpfile = onode_getpar_chr_dflt(tso,"plot_filename",NULL);
  if (dumpfile != NULL){
    if ((strcmp(dumpfile,"null")!=0)&&(strcmp(dumpfile,"NULL")!=0))
      append_farray_plot(dumpfile,"phase_shift_vs_time",tshift,tn,1);
    myfree(dumpfile);
  }

  //
  //  Compute the temporal window amplitude:  'tamp[tn]'
  //
  tfo = onode_child_get_unique(fo,"filter");  // The 1d temporal filter
  tamp = kernu_o_filter(mylogf,tfo,tscale,&tn);  // ('tn' is not changed)

  //
  //  Create the four curved filters
  //
  *fpe = gabor_space_time_tshift(xn,yn,tn,sscale,tscale,sdx,sdy,sf,theta,0.0,
				 tamp,tshift);
  *fpo = gabor_space_time_tshift(xn,yn,tn,sscale,tscale,sdx,sdy,sf,theta,90.0,
				 tamp,tshift);
  *fne = gabor_space_time_tshift(xn,yn,tn,sscale,tscale,sdx,sdy,sf,atheta,0.0,
				 tamp,tshift);
  *fno = gabor_space_time_tshift(xn,yn,tn,sscale,tscale,sdx,sdy,sf,atheta,90.0,
				 tamp,tshift);

  myfree(tshift);  // Free storage for these 1D arrays
  myfree(tamp);
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_ARB_OPPONENT_GABOR                          */
/*                                                                           */
/*  Arbitrary (rather than quadrature) phase shift.                          */
/*                                                                           */
/*****************************************************************************/
void get_arb_opponent_gabor(xn,yn,tn,sscale,tscale,sr,st,fr,ft,theta,tfilt,ph,
			    fpe,fpo,fne,fno)
     int xn,yn,tn;
     float sscale,tscale,sr,st,fr,ft,theta;
     char *tfilt;
     float ph;  // Phase shift (deg) 90.0 for quadrature
     float ****fpe,****fpo,****fne,****fno;
{
  float ntheta;

  //
  //  *** 'sr' is used twice below, for sd_orth and sd_par
  //

  ntheta = theta + 180.0;
  *fpe = gabor_space_time_tensor(xn,yn,tn,sscale,tscale,sr,sr,st,fr,ft,theta,
				 0.0,tfilt);
  *fpo = gabor_space_time_tensor(xn,yn,tn,sscale,tscale,sr,sr,st,fr,ft,theta,
				 ph,tfilt);
  *fne = gabor_space_time_tensor(xn,yn,tn,sscale,tscale,sr,sr,st,fr,ft,ntheta,
				 0.0,tfilt);
  *fno = gabor_space_time_tensor(xn,yn,tn,sscale,tscale,sr,sr,st,fr,ft,ntheta,
				 ph,tfilt);
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_QUAD_OPP_TILT_GABOR                         */
/*                                                                           */
/*****************************************************************************/
void get_quad_opp_tilt_gabor(xn,yn,tn,sscale,tscale,sr,st,fr,ft,a,tfilt,ss,
			     rotflag,fpe,fpo,fne,fno)
     int xn,yn,tn;
     float sscale,tscale,sr,st,fr,ft,a;  // 'a' is angle theta
     char *tfilt;
     float ss;      // Shift scale, or shift fraction, to control shift
     int rotflag;
     float ****fpe,****fpo,****fne,****fno;
{
  float na;

  na = a + 180.0;

  if (rotflag == 1){
    *fpe = gabor_tilt3_tensor(xn,yn,tn,sscale,tscale,sr,st,fr,ft,a,0,tfilt);
    *fpo = gabor_tilt3_tensor(xn,yn,tn,sscale,tscale,sr,st,fr,ft,a,1,tfilt);
    *fne = gabor_tilt3_tensor(xn,yn,tn,sscale,tscale,sr,st,fr,ft,na,0,tfilt);
    *fno = gabor_tilt3_tensor(xn,yn,tn,sscale,tscale,sr,st,fr,ft,na,1,tfilt);
    /*
    *fpe = gabor_tilt2_tensor(xn,yn,tn,sscale,tscale,sr,st,fr,ft,a,0,tfilt);
    *fpo = gabor_tilt2_tensor(xn,yn,tn,sscale,tscale,sr,st,fr,ft,a,1,tfilt);
    *fne = gabor_tilt2_tensor(xn,yn,tn,sscale,tscale,sr,st,fr,ft,na,0,tfilt);
    *fno = gabor_tilt2_tensor(xn,yn,tn,sscale,tscale,sr,st,fr,ft,na,1,tfilt);
    */
  }else{
    *fpe = gabor_tilt_tensor(xn,yn,tn,sscale,tscale,sr,st,fr,ft,a,0,tfilt,ss);
    *fpo = gabor_tilt_tensor(xn,yn,tn,sscale,tscale,sr,st,fr,ft,a,1,tfilt,ss);
    *fne = gabor_tilt_tensor(xn,yn,tn,sscale,tscale,sr,st,fr,ft,na,0,tfilt,ss);
    *fno = gabor_tilt_tensor(xn,yn,tn,sscale,tscale,sr,st,fr,ft,na,1,tfilt,ss);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_QUAD_GABOR                              */
/*                                                                           */
/*****************************************************************************/
void get_quad_gabor(ppl,o,xn,yn,tn,sscale,tscale,fe,fo)
     struct param_pair_list *ppl;
     struct onode *o;
     int xn,yn,tn;
     float sscale,tscale;
     float ****fe,****fo;
{
  float ph_off;

  ph_off = 0;
  *fe = gabor_space_time_sep_tensor(ppl,o,xn,yn,tn,sscale,tscale,ph_off);

  ph_off = 90.0;
  *fo = gabor_space_time_sep_tensor(ppl,o,xn,yn,tn,sscale,tscale,ph_off);

}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_PPL_GABOR                               */
/*                                                                           */
/*****************************************************************************/
float ***get_ppl_gabor(ppl,xn,yn,tn,sscale,tscale)
     struct param_pair_list *ppl;    // Contains filter parameters
     int xn,yn,tn;                   // Dimensions of 3D filter to create
     float sscale,tscale;            // spatial and temporal units
{
  float dir,tf,sf,sdo,sdt,ph;
  float ***f;
  char *ftype,tfilt[SLEN];

  dir   = param_getf_exit(ppl,"dir");
  tf    = param_getf_exit(ppl,"tf");
  sf    = param_getf_exit(ppl,"sf");
  sdo   = param_getf_exit(ppl,"s_sd");
  sdt   = param_getf_exit(ppl,"t_sd");
  ph    = param_getf_exit(ppl,"ph");
  ftype = param_getc_exit(ppl,"type");

  if (strcmp(ftype,"Gabor")==0)
    strcpy(tfilt,"gaussian");
  else if (strcmp(ftype,"Gabor_causal_t_1")==0)
    strcpy(tfilt,"cwq_temp");
  else
    exit_error("GET_PPL_GABOR","Unknown filter type");
    
  f = gabor_space_time_tensor(xn,yn,tn,sscale,tscale,sdo,sdo,sdt,sf,tf,dir,
			      ph,tfilt);

  myfree(ftype);

  return f;
}
/**************************************-**************************************/
/*                                                                           */
/*                            DOG_SPACE_TIME_TENSOR                          */
/*                                                                           */
/*  Return a difference of Gaussians function in 2 spatial dimensions and    */
/*  time.                                                                    */
/*                                                                           */
/*  Set 'outfile' to "" to avoid dumping diagnostic plots.                   */
/*                                                                           */
/*  gtsig - SD of temporal Gaussian to multiply temporal filter (ms)         */
/*  coffx - center x offset relative to surround                             */
/*  coffy - center y offset relative to surround                             */
/*                                                                           */
/*****************************************************************************/
float ***dog_space_time_tensor(xn,yn,zn,sscale,tscale,sig1,sig2,amp1,amp2,
			       coffx,coffy,sdelay,gtsig,gtcent,
			       outfile,ppl,o,pflag,hilb_flag,t0flag,toff,
			       tau_ad)
     int xn,yn,zn; // must be odd numbers, center is (xn-1)/2,...,(zn-1/2)
     float sscale,tscale; /* space and time scaling */
     float sig1,sig2,amp1,amp2,coffx,coffy,sdelay;
     float gtsig;
     float gtcent;                // Center of multiplying Gaussian
     char outfile[];
     struct param_pair_list *ppl; // Model parameter pair list
     struct onode *o;             // onode version
     int pflag;
     int hilb_flag;  // 1-take hilbert transform, before Gauss mult
                     // WARNING:  this spoils causality
     int t0flag;     // 0-t0 is at beginning, 1-t0 is at middle
     float toff;     // Time offset in seconds for AB temporal filter
     float tau_ad;   // Addition long time constant exp decay, for adaptation
{
  int i,j,k;
  int xc,yc,d; // centers of dimensions
  int t2,ti;
  float ***filter;
  float **g1,**g2,*t,tsig,tcent,*gt,sump,sumn,sumg1,sumg2,tt;
  // FOR OTHER TEMPORAL FILTER
  float a,b,tau1,tau2,cxc,cyc;
  int n1,n2,kflag,tst;
  char *tfilt;

  if (pflag){
    printf("  DOG_SPACE_TIME_TENSOR\n");
    printf("    xn,yn,zn   %d,%d,%d\n",xn,yn,zn);
  }

  kflag = 0;

  // WYETH BUGFIX 14th Oct, 2008 - changed to work w/ mod_dog_util
  //   so that the stim and response are aligned
  //xc = (xn-1)/2;
  //yc = (yn-1)/2;
  xc = xn/2;
  yc = yn/2;

  //  Make DELTA-Function, and RETURN
  if (sig1 == -1.0){
    if (pflag) printf("  *** sig1 is -1:  USING DELTA FUNCTION at t=0.\n");
    filter = f3tensor(1,xn,1,yn,1,zn); //  For 3D FFT
    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){
	for(k=0;k<zn;k++){
	  filter[i+1][j+1][k+1] = 0.0;
	}
      }
    }
    filter[xc+1][yc+1][1] = 11.0; // WYETH - are the +1's right?,   WHY 11??
    return filter;
  }

  cxc = (float)xc + coffx/sscale; // center pixel x-coord for CENTER
  cyc = (float)yc + coffy/sscale; // center pixel y-coord for CENTER

  d = my_rint(sdelay);

  tst = my_rint(1000.0*tscale); // msec per time unit

  //
  //  Make two 2D Gaussian filters, 'g1' and 'g2'
  //
  if (pflag) printf("    amp1= %.2f  amp2= %.2f\n",amp1,amp2);
  tsig = sig1/sscale; // In pixels
  if (pflag) printf("    sig1 = %.2f  sscale = %.4f  ==> %.2f pixels\n",
		    sig1,sscale,tsig);
  g1 = gaussian_2d(xn,yn,cxc,cyc,tsig,amp1);
  tsig = sig2/sscale; // In pixels
  if (pflag) printf("    sig2 = %.2f  sscale = %.4f  ==> %.2f pixels\n",
		    sig2,sscale,tsig);
  g2 = gaussian_2d(xn,yn,(float)xc,(float)yc,tsig,amp2);

  if (outfile != NULL){
    remove_file(outfile);
    append_farray_plot(outfile,"Gauss_1_raw",g1[my_rint(cxc)],yn,1);
    append_farray_plot(outfile,"Gauss_2_raw",g2[xc],yn,1);
  }

  //
  //  Default temporal filter is 'ab'
  //
  if (ppl != NULL)
    tfilt = paramfile_get_char_param_default(ppl,"dogt_filter","ab");
  else{
    if (onode_test_ostr(o,"dogt_filter"))
      exit_error("KERNEL_UTIL, DOG_SPACE_TIME_TENSOR",
		 "Use 'tfilter', not 'dogt_filter'");
    tfilt = onode_getpar_chr_dflt(o,"tfilter","ab");
  }

  if (t0flag == 0) // Set t-initial
    ti = zn-1;
  else
    ti = zn-1 - zn/2;

  if (strcmp(tfilt,"wa")==0){ // Watson/Ahumada
    a = 1.0;
    b = 0.8;
    tau1 = 0.002;  // 0.004;
    tau2 = 0.0026; //0.0053;
    n1 = 9;
    n2 = 10;
    t = wa_temporal(2*zn-1,tscale,a,b,n1,n2,tau1,tau2);
  }else if (strcmp(tfilt,"ab")==0){ // Adelson Bergen
    float tk;
    int tn;

    if (ppl != NULL){
      tk = paramfile_get_float_param_or_exit(ppl,"tab_k");
      tn = paramfile_get_int_param_or_exit(ppl,"tab_n");
    }else{
      tk = onode_getpar_flt_exit(o,"tab_k");
      tn = onode_getpar_int_exit(o,"tab_n");
    }

    t = ab_temporal(2*zn-1,tscale,tn,tk,toff); // centers at middle
  }else if (strcmp(tfilt,"dexp")==0){
    float t0,tau1,tau2,tamp2,sigs;

    if (ppl != NULL){
      t0   = paramfile_get_float_param_default(ppl,"dogt_t0"  ,0.020);
      tau1 = paramfile_get_float_param_default(ppl,"dogt_tau1",0.003);
      tau2 = paramfile_get_float_param_default(ppl,"dogt_tau2",0.010);
      tamp2 = paramfile_get_float_param_default(ppl,"dogt_amp2",0.500);
      sigs = paramfile_get_float_param_default(ppl,"dogt_sigs",0.004);
    }else{
      t0   = onode_getpar_flt_dflt(o,"dogt_t0"  ,0.020);
      tau1 = onode_getpar_flt_dflt(o,"dogt_tau1",0.003);
      tau2 = onode_getpar_flt_dflt(o,"dogt_tau2",0.010);
      tamp2 = onode_getpar_flt_dflt(o,"dogt_amp2",0.500);
      sigs = onode_getpar_flt_dflt(o,"dogt_sigs",0.004);
    }

    t = dd_exp_temporal(2*zn-1,tscale,t0,tau1,tau2,tamp2,sigs);
  }else if (strcmp(tfilt,"tfilt_01")==0){
    t = filt_01_temporal(ppl,o,2*zn-1,tscale,pflag);
  }else if (strcmp(tfilt,"maxwell")==0){
    float tmax,s;

    if (ppl != NULL){
      tmax = paramfile_get_float_param_or_exit(ppl,"maxwell_tmax");
    }else{
      tmax = onode_getpar_flt_exit(o,"maxwell_tmax");
    }
    s = tmax / sqrt(2.0) / tscale;
    t = maxwell_farray((float)(zn-1),s,1.0,2*zn-1);

  }else if (strcmp(tfilt,"delta")==0){
    t = get_zero_farray(2*zn-1);
    t[zn-1] = 1.0;
  }else{
    t = NULL;
    exit_error("DOG_SPACE_TIME_TENSOR","Unknown temporal filter");
  }
  myfree(tfilt); // Free the name of the temporal filter

  if (outfile != NULL){
    append_farray_plot(outfile,"raw_t_filter_(ms)",&(t[ti]),zn,tst);
  }

  
  if (hilb_flag == 1){  // WYETH - hilbert transform
    float *thil;

    thil = fft_hilbert(t,2*zn-1);
    //append_farray_plot("zz.out.pl","t-Hilbert",thil,2*zn-1,tst);
    myfree(t);
    t = thil;
  }

  //  Multiply temporal filter 't' by Gaussian window
  if (gtsig > 0.0){
    tsig = gtsig / (1000.0*tscale);
    tcent = gtcent / (1000.0*tscale);
    gt = gaussian_one_farray((float)(zn-1)+tcent,tsig,1.0,2*zn-1);
    multiply_farrays_in_place(t,gt,2*zn-1);
    if (outfile != (char *)NULL){
      append_farray_plot(outfile,"Temporal_Gaussian_(ms)",&(gt[ti]),zn,tst);
    }
    myfree(gt);
  }


  //
  //  April 29, 2011 - added to allow long-time adaptation in LGN
  //  Motivated by Douglas building a model of two-stage adaptation.
  //
  if (tau_ad > 0.0){
    int tn;
    float *texp,tau;

    tn = 2*zn-1;

    norm_area_farray(t,tn,1.0);  // Will exit if area is zero.

    tau = tau_ad / tscale;
    texp = one_sided_exp_farray((float)zn,tau,1.0,tn);
    norm_area_farray(texp,tn,-1.0);  // Will exit if area is zero.

    //
    //  Put zeros in 2nd half of filter, to avoid interference during the
    //  wrap-around assumption of the FFT convolution.
    //
    for(i=(zn + zn/2);i<tn;i++)
      texp[i] = 0.0;

    add_to_farray(t,texp,tn);

    myfree(texp);
  }


  make_max_const_farray(t,2*zn-1,1.0);

  if (outfile != (char *)NULL){
    append_farray_plot(outfile,"Temporal_filter_(ms)",&(t[ti]),zn,tst);
  }
    
  filter = f3tensor(1,xn,1,yn,1,zn); // use for 3D FFT
  sump = 0.0;
  sumn = 0.0;
  sumg1 = 0.0;
  sumg2 = 0.0;

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      sumg1 += g1[i][j];
      sumg2 += g2[i][j];
      for(k=0;k<zn;k++){
	t2 = ti + k - d;
	if ((t2 < 0) || (t2 >= (2*zn-1))){ // Check for delay overrun
	  tt = t[ti+k]*g1[i][j];  // If overrun, assume t[t2] is 0
	}else{
	  tt = t[ti+k]*g1[i][j] - t[t2]*g2[i][j];
	}
	filter[i+1][j+1][k+1] = tt;
	if (tt > 0.0)
	  sump += tt;
	else
	  sumn += tt;
      }
    }
  }

  if (pflag){
    if (amp2 != 0.0){
      printf("    positive and negative sums: %.4f/%.4f = %.4f\n",
	     sump,sumn,sump/-sumn);
      printf("    ratio of Gaussian area (ctr/sur): %.4f/%.4f = %.4f\n",
	     sumg1,sumg2,sumg1/sumg2);
    }else{
      printf("    positive and negative sums: %.4f and %.4f\n",
	     sump,sumn);
      printf("    Gaussian area, ctr and sur: %.4f and %.4f\n",
	     sumg1,sumg2);
    }
  }

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                              SEPARABLE_SPACE                              */
/*                                                                           */
/*  Return a 2-D Gabor spatial filter that is seperable in x and y.          */
/*                                                                           */
/*    sf     spatial frequency                                               */
/*    sigma  spatial spread                                                  */
/*                                                                           */
/*  *** NOTE:  On 6/7/96, a 2.0 was added in the denominator of the          */
/*             Gaussian "exp" function.                                      */
/*                                                                           */
/*****************************************************************************/
float **separable_space(xn,yn,sscale,sigma,sf,phase)
     int xn,yn; /* must be odd numbers, center is 0,0,0 */
     float sscale; /* spacial scaling */
     float sigma,sf;
     int phase;  /* 0 for even, 1 for odd */
{
  int i,j;
  float **filter;
  float twopisf;
  float x,y,gx,gy,s2;
  int xc,yc; /* centers of dimensions */

  if (PFLAG) printf("  SEPARABLE_SPACE\n");

  xc = (xn-1)/2;
  yc = (yn-1)/2;
  twopisf = 2.0*M_PI*sf;
  s2 = 2.0*sigma*sigma;

  filter = get_2d_farray(xn,yn);
  if (phase==0){  /*** Even "cos" filter ***/
    for(i=0;i<xn;i++){
      x = (float)(i-xc)*sscale;
      gx = exp(-x*x/s2);
      for(j=0;j<yn;j++){
	y = (float)(j-yc)*sscale;
	gy = exp(-y*y/s2) * cos(twopisf*y);
	filter[i][j] = gx * gy;
      }
    }
  }else{
    for(i=0;i<xn;i++){
      x = (float)(i-xc)*sscale;
      gx = exp(-x*x/s2);
      for(j=0;j<yn;j++){
	y = (float)(j-yc)*sscale;
	gy = exp(-y*y/s2) * sin(twopisf*y);
	filter[i][j] = gx * gy;
      }
    }
  }
  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                                  STF3D_01                                 */
/*                                                                           */
/*  Create a 3D spatiotemporal filter.  Use the "AB_TEMPORAL" filter         */
/*  with a Gabor spatial filter.  The phase of the temporal filter           */
/*  is determined by parameters "n" and "k", and the phase of the            */
/*  spatial Gabor filter can be even or odd.                                 */
/*                                                                           */
/*****************************************************************************/
float ***stf3d_01(xn,yn,zn,sscale,tscale,ssdo,ssdp,sf,sphase,n,k,theta,
		  writeflag)
     int   xn,yn,zn;
     float sscale,tscale;    // (deg/pix), (s/pix)
     float ssdo,ssdp;        // spatial SD orth. and parallel (deg)
     float sf;               // spatial frequency (cyc/deg)
     int sphase;             // Spatial phase flag: 0/1 ==> 0/90.0 (deg)
     int n;                  // AB temporal param 'n'
     float k;                // AB temporal param 'k'
     float theta;            // direction (deg)
     int writeflag;
{
  int i,j,t;
  float ***filter,*tfilter,**sfilter,s,*col;
  float xc,yc,sdo,sdp,sfpix;  // New for 'gabor_2d_space_raw'
  char name[SLEN];

  if (PFLAG) printf("  STF3D_01\n");

  tfilter = ab_temporal(zn,tscale,n,k,0.0);

  xc = (float)((int)((xn-1)/2));
  yc = (float)((int)((yn-1)/2));
  sdo = ssdo / sscale;  // Convert to (pix)
  sdp = ssdp / sscale;  // Convert to (pix)
  sfpix = sf * sscale;  // Convert to (cyc/pix)

  if (sphase == 0)
    sfilter = gabor_2d_space_raw(xn,yn,xc,yc,sdo,sdp,sfpix,theta, 0.0);
  else
    sfilter = gabor_2d_space_raw(xn,yn,xc,yc,sdo,sdp,sfpix,theta,90.0);

  if (writeflag){
    sprintf(name,"zzz.k%d.n%d.tfilter.pl",(int)k,n);
    write_farray_plot(name,tfilter,zn);
    sprintf(name,"zzz.sf%.2f.ph%d_tfilter.pl",sf,sphase);
    col = get_column_2d_farray(sfilter,xn,yn,yn/2);
    write_farray_plot(name,col,yn);
    myfree(col);
  }

  filter = get_3d_farray(xn,yn,zn);
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      s = sfilter[i][j];
      for(t=0;t<zn;t++)
	filter[i][j][t] = s * tfilter[t];
    }
  }
  myfree(tfilter);
  free_2d_farray(sfilter,xn);

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                                  STF3D_02                                 */
/*                                                                           */
/*  Create a 3D spatiotemporal filter.  Use the "CWQ_TEMPORAL" filter        */
/*  with a Gabor spatial filter.                                             */
/*                                                                           */
/*****************************************************************************/
float ***stf3d_02(int xn, int yn, int zn,     // Size of filter
		  float sscale, float tscale, // 
		  float theta,                // Direction [0..360]
		  float sdorth, float sdpar,  // SD orth and parallel (deg)
		  float sf, float sphase,     // SF and spatial phase (deg)
		  int a, float tau, float fr, float phi,  // Temporal pars
		  int writeflag){             // 1-Write filter plots
  int i,j,t;
  float xc,yc,sdo,sdp,sfpix;
  float ***filter,*tfilter,**sfilter,s,*col;
  char name[SLEN];

  if (PFLAG) printf("  STF3D_02\n");

  tfilter = cwq_temporal(zn,tscale,a,tau,fr,phi,0.0);

  xc = (xn-1)/2;
  yc = (yn-1)/2;
  //printf("sdorth = %f  sdpar = %f  sf = %f  theta = %f  sphase = %f\n",
  //sdorth,sdpar,sf,theta,sphase);

  sdo = sdorth / sscale;  // Convert to pixels
  sdp = sdpar / sscale;  // Convert to pixels
  sfpix = sf * sscale;  // Convert to cyc/pix
  sfilter = gabor_2d_space_raw(xn,yn,xc,yc,sdo,sdp,sfpix,theta,sphase);

  if (writeflag){
    sprintf(name,"zzz.tfilter.%d.pl",(int)phi);
    write_farray_plot(name,tfilter,zn);
    sprintf(name,"zzz.sfilter.%d.pl",(int)sphase);
    col = get_column_2d_farray(sfilter,xn,yn,yn/2);
    write_farray_plot(name,col,yn);
    myfree(col);
  }

  filter = get_3d_farray(xn,yn,zn);
  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++){
      s = sfilter[i][j];
      for(t=0;t<zn;t++)
	filter[i][j][t] = s * tfilter[t];
    }
  myfree(tfilter); free_2d_farray(sfilter,xn);

  return filter;
}
/**************************************-**************************************/
/*                                                                           */
/*                        QUAD_OPPONENT_CAUSAL_STF3D_01                      */
/*                                                                           */
/*  Return up-even, up-odd, down-even, and down-odd filters.  The            */
/*  origin of the returned filters is [(xn-1)/2,(yn-1)/2,0].                 */
/*                                                                           */
/*  *** WYETH SEE NEW VERSION BELOW (this one doesn't norm filters)          */
/*  *** WYETH SEE NEW VERSION BELOW (this one doesn't norm filters)          */
/*  *** WYETH SEE NEW VERSION BELOW (this one doesn't norm filters)          */
/*                                                                           */
/*****************************************************************************/
void quad_opponent_causal_stf3d_01(xn,yn,tn,sscale,tscale,sigma,sf,
				   n1,n2,k,fuoc,fuec,fdoc,fdec,wrflag)
     int xn,yn,tn;
     float sscale,tscale,sigma,sf;
     int n1,n2;
     float k,****fuoc,****fuec,****fdoc,****fdec;
     int wrflag;
{
  float ***filter1,***filter2,***filter3,***filter4;  /* SEPARABLE */
  float ***fuo,***fue,***fdo,***fde; /* ORIENTED */
  float sd;
  int tn2,sphase;

  if (PFLAG) printf("  QUAD_OPPONENT_CAUSAL_STF3D_01\n");

  tn2 = tn*2;

  sd = sigma;

  sphase = 1;
  filter1 = stf3d_01(xn,yn,tn2,sscale,tscale,sd,sd,sf,sphase,n1,k,90.0,wrflag);
  filter2 = stf3d_01(xn,yn,tn2,sscale,tscale,sd,sd,sf,sphase,n2,k,90.0,wrflag);
  sphase = 0;
  filter3 = stf3d_01(xn,yn,tn2,sscale,tscale,sd,sd,sf,sphase,n2,k,90.0,wrflag);
  filter4 = stf3d_01(xn,yn,tn2,sscale,tscale,sd,sd,sf,sphase,n1,k,90.0,wrflag);

  fdo = add_3d_farrays(filter1,filter3,0,xn,0,yn,0,tn2);
  fue = add_3d_farrays(filter2,filter4,0,xn,0,yn,0,tn2);
  multiply_3d_farray(filter1,0,xn,0,yn,0,tn2,-1.0);
  multiply_3d_farray(filter2,0,xn,0,yn,0,tn2,-1.0);
  fuo = add_3d_farrays(filter1,filter3,0,xn,0,yn,0,tn2);
  fde = add_3d_farrays(filter2,filter4,0,xn,0,yn,0,tn2);

  free_3d_farray(filter1,xn,yn,tn2); free_3d_farray(filter2,xn,yn,tn2);
  free_3d_farray(filter3,xn,yn,tn2); free_3d_farray(filter4,xn,yn,tn2);

  /*** Since the temporal filters are causal, keep only t>=0 ***/
  get_sub_3d_farray(fuo,xn,yn,tn2,fuoc,0,xn,0,yn,tn,tn);
  get_sub_3d_farray(fue,xn,yn,tn2,fuec,0,xn,0,yn,tn,tn);
  get_sub_3d_farray(fdo,xn,yn,tn2,fdoc,0,xn,0,yn,tn,tn);
  get_sub_3d_farray(fde,xn,yn,tn2,fdec,0,xn,0,yn,tn,tn);

  free_3d_farray(fuo,xn,yn,tn2); free_3d_farray(fue,xn,yn,tn2);
  free_3d_farray(fdo,xn,yn,tn2); free_3d_farray(fde,xn,yn,tn2);

  if (PFLAG) printf("  QUAD_OPPONENT_CAUSAL_STF3D_01  Done.\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                    QUAD_OPPONENT_CAUSAL_STF3D_01_TENSOR                   */
/*                                                                           */
/*  Like the above, but returns tensor format.                               */
/*                                                                           */
/*  *** Separable filters are normalized.                                    */
/*                                                                           */
/*****************************************************************************/
void quad_opponent_causal_stf3d_01_tensor(xn,yn,tn,sscale,tscale,sdo,sdp,sf,n1,
					  n2,k,theta,fuoc,fuec,fdoc,fdec,
					  wrflag)
     int xn,yn,tn;
     float sscale,tscale,sdo,sdp,sf;
     int n1,n2;
     float k,theta,****fuoc,****fuec,****fdoc,****fdec;
     int wrflag;
{
  float ***filter1,***filter2,***filter3,***filter4;  // Separable
  float ***fuo,***fue,***fdo,***fde;                  // Oriented
  float sum;
  int sph;

  sph = 1;  // Spatial phase flag, 1 = 90 deg
  filter1 = stf3d_01(xn,yn,tn,sscale,tscale,sdo,sdp,sf,sph,n1,k,theta,wrflag);
  sum = sum_square_3d_farray(filter1,0,xn,0,yn,0,tn);
  multiply_3d_farray(filter1,0,xn,0,yn,0,tn,1.0/sqrt(sum));
  filter2 = stf3d_01(xn,yn,tn,sscale,tscale,sdo,sdp,sf,sph,n2,k,theta,wrflag);
  sum = sum_square_3d_farray(filter2,0,xn,0,yn,0,tn);
  multiply_3d_farray(filter2,0,xn,0,yn,0,tn,1.0/sqrt(sum));

  sph = 0;  // Spatial phase flag, 0 = 0 deg
  filter3 = stf3d_01(xn,yn,tn,sscale,tscale,sdo,sdp,sf,sph,n2,k,theta,wrflag);
  sum = sum_square_3d_farray(filter3,0,xn,0,yn,0,tn);
  multiply_3d_farray(filter3,0,xn,0,yn,0,tn,1.0/sqrt(sum));
  filter4 = stf3d_01(xn,yn,tn,sscale,tscale,sdo,sdp,sf,sph,n1,k,theta,wrflag);
  sum = sum_square_3d_farray(filter4,0,xn,0,yn,0,tn);
  multiply_3d_farray(filter4,0,xn,0,yn,0,tn,1.0/sqrt(sum));

  //
  //  Now, filter1...4 contain the separable filters, and next we will
  //  add/subtract them to create the inseparable (tilted) filters.
  //
  //  Note, 'u' and 'd' probably stand for 'up' and 'down',
  //        and 'o' and 'e' stand for 'odd' and 'even'.
  //

  fuo = add_3d_farrays(filter1,filter3,0,xn,0,yn,0,tn);
  fde = add_3d_farrays(filter2,filter4,0,xn,0,yn,0,tn);
  multiply_3d_farray(filter1,0,xn,0,yn,0,tn,-1.0);
  multiply_3d_farray(filter2,0,xn,0,yn,0,tn,-1.0);
  fdo = add_3d_farrays(filter1,filter3,0,xn,0,yn,0,tn);
  fue = add_3d_farrays(filter2,filter4,0,xn,0,yn,0,tn);

  //
  //  These four filters are returned to the caller.
  //
  *fuoc = convert_3d_to_tensor(fuo,xn,yn,tn);
  *fuec = convert_3d_to_tensor(fue,xn,yn,tn);
  *fdoc = convert_3d_to_tensor(fdo,xn,yn,tn);
  *fdec = convert_3d_to_tensor(fde,xn,yn,tn);

  //
  //  Free memory in the working arrays.
  //
  free_3d_farray(filter1,xn,yn,tn); free_3d_farray(filter2,xn,yn,tn);
  free_3d_farray(filter3,xn,yn,tn); free_3d_farray(filter4,xn,yn,tn);
  free_3d_farray(fuo,xn,yn,tn); free_3d_farray(fue,xn,yn,tn);
  free_3d_farray(fdo,xn,yn,tn); free_3d_farray(fde,xn,yn,tn);
}
/**************************************-**************************************/
/*                                                                           */
/*                             BDE_CWQ_STF3D_UTIL                            */
/*                                                                           */
/*  Return filter in tensor format.                                          */
/*                                                                           */
/*****************************************************************************/
float ***bde_cwq_stf3d_util(int xn, int yn, int tn,      // filter size
			    float sscale, float tscale,  // resolution
			    float theta,                 // direction
			    float sf,                    // spatial freq
			    float sdo, float sdp,        // spatial Gauss.
			    float ph,                    // spatial phase
			    int a, float tau, float tf,  // temporal params
			    float tph,                   // temporal phase
			    float eta,                   // weight of f2
			    int write_flag){             // write output

  float sum,***f1,***f2,***f,***ftens;

  f1 = stf3d_02(xn,yn,tn,sscale,tscale,theta,sdo,sdp,sf,ph,
		a,tau,tf,tph,write_flag);
  sum = sum_square_3d_farray(f1,0,xn,0,yn,0,tn);
  multiply_3d_farray(f1,0,xn,0,yn,0,tn,1.0/sqrt(sum));
  //write_3d_data_part("zzz.f1.3d",f1,0,xn,0,yn,0,tn,4,2,1);

  ph  -= 90.0;   // Shift spatial phase by 90 deg
  tph -= 90.0;   // Shift temporal phase by 90 deg
  f2 = stf3d_02(xn,yn,tn,sscale,tscale,theta,sdo,sdp,sf,ph,
		a,tau,tf,tph,write_flag);
  sum = sum_square_3d_farray(f2,0,xn,0,yn,0,tn);
  multiply_3d_farray(f2,0,xn,0,yn,0,tn,1.0/sqrt(sum));
  //write_3d_data_part("zzz.f2.3d",f2,0,xn,0,yn,0,tn,4,2,1);

  //
  //  Add f1 + eta * f2 to get a tilted FL1 (Filter Left 1)
  //

  multiply_3d_farray(f2,0,xn,0,yn,0,tn,eta);
  f = add_3d_farrays(f1,f2,0,xn,0,yn,0,tn);
  //write_3d_data_part("zzz.f.3d",f,0,xn,0,yn,0,tn,4,2,1);

  free_3d_farray(f1,xn,yn,tn);
  free_3d_farray(f2,xn,yn,tn);

  ftens = convert_3d_to_tensor(f,xn,yn,tn);

  free_3d_farray(f,xn,yn,tn);

  return ftens;
}
/**************************************-**************************************/
/*                                                                           */
/*                            BDE_CWQ_STF3D_TENSOR                           */
/*                                                                           */
/*****************************************************************************/
void bde_cwq_stf3d_tensor(fo,mo,xn,yn,tn,sscale,tscale,fl1,fr1,fl2,fr2)
     struct onode *fo; // Filter onode
     struct onode *mo; // Model onode  WYETH 2014 March, we need some params
     int xn,yn;       // spatial size
     int tn;          // temporal size
     float sscale;    // deg / pix
     float tscale;    // samp / sec
     float ****fl1;   // [1..xn][1..yn][1..tn]
     float ****fr1;
     float ****fl2;
     float ****fr2;
{
  int a,wrflag;
  int qsopp,phopp;  // WYETH March 2014
  float sf,sdo,sdp,theta,tau,tf,phi,eta,phase,ph,phshift,tph;

  wrflag  = onode_getpar_int_dflt(fo,"write_filter",0);
  sf      = onode_getpar_flt_exit(fo,"sf");
  sdo     = onode_getpar_flt_exit(fo,"sd_orth");
  sdp     = onode_getpar_flt_exit(fo,"sd_par");
  phase   = onode_getpar_flt_exit(fo,"phase");
  theta   = onode_getpar_flt_exit(fo,"direction");
  //phshift = onode_getpar_flt_dflt(fo,"phase_shift",0.0);
  a       = onode_getpar_int_exit(fo,"alpha");
  tau     = onode_getpar_flt_exit(fo,"tau");
  tf      = onode_getpar_flt_exit(fo,"tf");
  phi     = onode_getpar_flt_exit(fo,"phi");
  eta     = onode_getpar_flt_exit(fo,"eta");

  //*************
  //*************
  //*************

  // WYETH HERE *** Pamela and I changed this code to clean up the parameters
  //  a bit  March 2014


  //binoc_shift = onode_getpar_flt_exit(mo,"phase_shift");
  //mod_me_binoc_phase1  = onode_getpar_flt_dflt(mo,"phase_1",0.0);

  phshift = onode_getpar_flt_dflt(mo,"phase_shift",0.0);
  phopp = onode_getpar_int_dflt(mo,"phase_shift_opp",0);
  qsopp = onode_getpar_int_dflt(mo,"quad_shift_opp",0);

  //*************
  //*************
  //*************

  //
  //  Compute each filter from sum of two space-time seperable filters
  //

  // GABOR EQUIV:  phas = mod_me_binoc_phase1;
  ph = -phase;
  tph = phi;
  *fl1 = bde_cwq_stf3d_util(xn,yn,tn,sscale,tscale,theta,sf,sdo,sdp,ph,
			   a,tau,tf,tph,eta,wrflag);

  // GABOR EQUIV:  phas = mod_me_binoc_phase1 + binoc_shift;
  ph = -phase - phshift;  // Add disparity offset
  tph = phi;
  *fr1 = bde_cwq_stf3d_util(xn,yn,tn,sscale,tscale,theta,sf,sdo,sdp,ph,
			   a,tau,tf,tph,eta,wrflag);

  //
  // WYETH - I have set this to +90 to be consistent with the default 
  //         Gabor BDE model.
  //
  if (qsopp == 0)
    phase += 90.0;  // Add 90 deg to spatial phase
  else
    phase -= 90.0;  // Subtract 90 deg to spatial phase

  //
  //  Repeat the above, creating a pair of L and R filters in 90 deg quadr.
  //
  // GABOR EQUIV:  phas = mod_me_binoc_phase1 + quadshift;
  ph = -phase;
  tph = phi;
  *fl2 = bde_cwq_stf3d_util(xn,yn,tn,sscale,tscale,theta,sf,sdo,sdp,ph,
			   a,tau,tf,tph,eta,wrflag);


  // GABOR EQUIV-0:  phas = mod_me_binoc_phase1 + binoc_shift + quadshift;
  // GABOR EQUIV-1:  phas = mod_me_binoc_phase1 - binoc_shift + quadshift;
  if (phopp == 0)
    ph = -phase - phshift;  // Add disparity offset
  else
    ph = -phase + phshift;  // Subtract disparity offset

  tph = phi;
  *fr2 = bde_cwq_stf3d_util(xn,yn,tn,sscale,tscale,theta,sf,sdo,sdp,ph,
			   a,tau,tf,tph,eta,wrflag);

  if (wrflag != 0){
    write_3d_data_part("zzz.fl1.3d",*fl1,1,xn,1,yn,1,tn,4,2,1);
    write_3d_data_part("zzz.fr1.3d",*fr1,1,xn,1,yn,1,tn,4,2,1);
    write_3d_data_part("zzz.fl2.3d",*fl2,1,xn,1,yn,1,tn,4,2,1);
    write_3d_data_part("zzz.fr2.3d",*fr2,1,xn,1,yn,1,tn,4,2,1);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                 GET_XT_DOG                                */
/*                                                                           */
/*  Return the 2D x-t filter that is a difference of Gaussians times a       */
/*  temporal filter determined by 'ttype'.                                   */
/*                                                                           */
/*  The temporal parameters 't1' ... are assigned according to 'ttype':      */
/*                                                                           */
/*  'ab'                                                                     */
/*     t1 = tn                                                               */
/*     t2 = 1/tk                                                             */
/*  'wa'                                                                     */
/*     t1 = a                                                                */
/*     t2 = b                                                                */
/*     t3 = n1                                                               */
/*     t4 = n2                                                               */
/*     t5 = tau1                                                             */
/*     t6 = tau2                                                             */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - amp1 and amp2 control the overall amplitude, since the temporal        */
/*    filter is not scaled here.                                             */
/*                                                                           */
/*****************************************************************************/
float **get_xt_dog(xn,tn,sscale,tscale,amp1,amp2,sig1,sig2,ttype,t1,t2,t3,t4,
		   t5,t6)
     int xn,tn;
     float sscale,tscale,amp1,amp2,sig1,sig2;
     char ttype[];
     float t1,t2,t3,t4,t5,t6;
{
  int i,j;
  float x,xc,*xdata,*tdata,**xtdog;
  
  /*** Get spatial DOG ***/
  xc = (float)xn/2.0 * sscale;
  xdata = (float *)myalloc(xn*sizeof(float));

  for(i=0;i<xn;i++){
    x = (float)i*sscale;
    xdata[i] = (float)(amp1*func_gaussian_one(x,xc,sig1) -
		       amp2*func_gaussian_one(x,xc,sig2));
  }
  append_farray_plot("yyy.xtdog","x",xdata,xn,1);
  
  /*** Get temporal filter ***/
  if (strcmp(ttype,"ab")==0){
    tdata = ab_temporal(tn,tscale,my_rint(t1),1.0/t2,0.0);
  }else if (strcmp(ttype,"wa")==0){
    tdata = wa_temporal(tn,tscale,t1,t2,my_rint(t3),my_rint(t4),t5,t6);
  }else{
    tdata = NULL;
    exit_error("GET_XT_DOG","Unknown type of temporal filter requested");
  }

  append_farray_plot("yyy.xtdog","t",tdata,tn,1);
  
  /*** Create 2D xt separable filter. ***/
  xtdog = get_2d_farray(xn,tn);
  for(i=0;i<xn;i++)
    for(j=0;j<tn;j++)
      xtdog[i][j] = xdata[i]*tdata[j];
  
  return xtdog;
}
