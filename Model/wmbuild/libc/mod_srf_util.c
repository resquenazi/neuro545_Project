/*****************************************************************************/
/*                                                                           */
/*   mod_srf_util.c                                                          */
/*   wyeth bair                                                              */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2013                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "my_util.h"
#include "farray_util.h"
#include "data_util.h"
#include "spike_util.h"
#include "paramfile_util.h"
#include "fft_util.h"
#include "kernel_util.h"
#include "mod_util.h"
#include "ifc_util.h"
#include "mod.h" // Data structures
#include "paramfile.h"

// For cluster computation  (static needed w/ initialization)
static int myid = -1;         // Process ID, -1-Single, 0-Master, >0-Slave
static char *mylogf = NULL;   // Logfile name
char ggstr[LONG_SLEN];        // Global log string


// Global variables
int     mod_srf_xn;
int     mod_srf_yn;
int     mod_srf_tn;
float   mod_srf_tscale;
float   mod_srf_sscale;     // deg/pix

int     mod_srf_mean0;        // 1 - subtract the mean (0 otherwise)
int     mod_srf_framei;       // index of single frame to process (-1 othws)
int     mod_srf_ttn;          // Response tn, may be longer than true 'tn'
float   mod_srf_norm_pow;     // 1-normalize total power to 1, 0-do not

float **mod_srf_map = NULL;   // the SRF model map [xn][yn/2]
int     mod_srf_fxn;          // x size of filter
int     mod_srf_fyn;          // y size of filter
float   mod_srf_fscale;       // Frequency scale ((cyc/deg) / pix)

/**************************************-**************************************/
/*                                                                           */
/*                              MOD_SRF_01_PREP                              */
/*                                                                           */
/*  Prepare constructs that will be needed when running trials.  Do things   */
/*  here that only need to be done once, to avoid repeating them for each    */
/*  trial.                                                                   */
/*                                                                           */
/*****************************************************************************/
void mod_srf_01_prep(m,r)
     struct model_struct *m;    // Model params
     struct response_struct *r; // Response params
{
  int xn,yn,tn,fxn,fyn;
  float sscale,tscale,o1,o2,scale1,scale2,nyqf;
  struct onode *fo;

  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  mylog(mylogf,"  MOD_SRF_01_PREP\n");

  xn     = onode_getpar_int_exit(m->o,"xn");
  yn     = onode_getpar_int_exit(m->o,"yn");
  tn     = onode_getpar_int_exit(m->o,"tn");
  tscale = onode_getpar_flt_exit(m->o,"tscale");
  sscale = onode_getpar_flt_exit(m->o,"sscale");

  mod_srf_mean0    = onode_getpar_int_dflt(m->o,"zero_mean",0);
  mod_srf_framei   = onode_getpar_int_dflt(m->o,"frame_index",-1);
  mod_srf_norm_pow = onode_getpar_flt_dflt(m->o,"norm_power",0.0);
  mod_srf_ttn      = onode_getpar_int_dflt(m->o,"tn_response",tn);

  if ((mod_srf_framei < -1) || (mod_srf_framei >= tn))
    exit_error("MOD_SRF_01_PREP","Invalid value for 'mod_srf_framei'");

  // Store in global variables for usage during run
  mod_srf_xn = xn;
  mod_srf_yn = yn;
  mod_srf_tn = tn;
  mod_srf_tscale = tscale;
  mod_srf_sscale = sscale;

  //
  //  Compute the scale of frequency to be used for the SRF map
  //
  nyqf = 1.0/(2.0*sscale); // Nyquist freq (1/2 sampling freq)
  mod_srf_fscale = nyqf / ((float)xn/2.0);
  mod_srf_fxn = xn;
  mod_srf_fyn = yn/2;

  fxn = mod_srf_fxn;
  fyn = mod_srf_fyn;

  //
  //  Create the SRF map
  //
  fo = onode_child_get_unique(m->o,"filter");
  if (fo == NULL)
    exit_error("MOD_SRF_01_PREP","Cannot find <filter>");


  o1 = (float)xn/2.0;
  o2 = 0.0;

  scale1 = scale2 = mod_srf_fscale;  // (cyc/deg) / pixel

  mod_srf_map = kernu_o_filt2d(mylogf,fo,o1,o2,scale1,scale2,&fxn,&fyn);

  // DEBUG
  write_2d_data("zz.srf.2d",mod_srf_map,0,0,fxn,fyn,4,2,1,1);
}
/**************************************-**************************************/
/*                                                                           */
/*                               MOD_SRF_01_DONE                             */
/*                                                                           */
/*  Free any storage here when done.                                         */
/*                                                                           */
/*****************************************************************************/
void mod_srf_01_done(m)
     struct model_struct *m;        // Model parameters
{
  mylog(mylogf,"  MOD_SRF_01_DONE\n");

  if (mod_srf_map != NULL)
    free_2d_farray(mod_srf_map,mod_srf_fxn);

}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_SRF_01_GET_RESPONSE                         */
/*                                                                           */
/*****************************************************************************/
void mod_srf_01_get_response(m,s,r)
     struct model_struct *m;        // Model parameters
     struct stim_struct *s;         // Stimulus parameters
     struct response_struct *r;     // Response parameters
{
  int i,j;
  int xn,yn,tn,fxn,fyn,n;
  float **d2d,**pow2d,*d,tscale,oldmean,rtot;

  xn     = mod_srf_xn;  // Short names for global variables
  yn     = mod_srf_yn;
  tn     = mod_srf_tn;
  tscale = mod_srf_tscale;
  fxn    = mod_srf_fxn;
  fyn    = mod_srf_fyn;

  if (mod_srf_framei >= 0){
    //
    //  Compute the power spectrum of a single 2d frame
    //

    // Extract the 2D frame at the specified index
    d2d = get_2d_farray(xn,yn);
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++){
	d2d[i][j] = s->d[i+1][j+1][mod_srf_framei+1];
      }

    sprintf(ggstr,"    Using only frame %d of (0..%d) frames\n",mod_srf_framei,
	    tn-1);
    mylog(mylogf,ggstr);

    if (mod_srf_mean0 == 1){
      subtract_mean_2d_farray(d2d,xn,yn,&oldmean);
      sprintf(ggstr,"    Subtracted mean (%f) from frame\n",oldmean);
      mylog(mylogf,ggstr);
    }else
      oldmean = 0.0;

    pow2d = power_2d(d2d,xn,yn,0);  // result is [xn][yn/2]

    free_2d_farray(d2d,xn);
  }else{
    //
    //  Process full 3D stimulus
    //
    exit_error("MOD_SRF_01_GET_RESPONSE","Full 3D not implemented yet");
  }

  //pow_2d_farray(pow2d,fxn,fyn,0.25);
  //write_2d_data("zz.pow.2d",pow2d,0,0,fxn,fyn,4,2,1,1);

  if (mod_srf_norm_pow > 0.0){
    sprintf(ggstr,"    Normalizing total power to %f\n",mod_srf_norm_pow);
    mylog(mylogf,ggstr);
    norm_area_2d_farray(pow2d,xn,yn/2,mod_srf_norm_pow);
  }


  //
  //  Multiply the SRF by the power spectrum
  //
  rtot = 0.0;
  for(i=0;i<mod_srf_fxn;i++){
    for(j=0;j<mod_srf_fyn;j++){
      rtot += mod_srf_map[i][j] * pow2d[i][j];
    }
  }

  printf("    rtot = %f\n",rtot);

  //
  //  The response vs. time is d[n]
  //
  n = mod_srf_ttn;  // Use the (possibly expanded) response tn
  d = get_const_farray(n,rtot);
  //append_farray_plot("zz.rawRespVT.pl","raw_stim",d,n,1);


  //
  //  Generate and return spike responses
  //
  {
    int ns,vn,seed,cpflag;
    float *fs,*v,samp;
    struct onode *sgo;

    seed = m->mseed[r->tsi];
    samp = 1.0/mod_srf_tscale;

    sgo = onode_child_get_unique(m->o,"spike_gen");
    if (sgo != NULL){
      ifc_util_poisson(mylogf,m,sgo,d,n,samp,seed,0,1,&fs,&ns,&v,&vn);
    }else{
      exit_error("MOD_SRF_01_GET_RESPONSE","Cannot find <spike_gen>");
    }

    cpflag = 0;       // Do not make copy, 'fs' will be linked to 'r'
    mod_util_resp_store_s(r,0,fs,ns,cpflag,mylogf);

    
  }
  myfree(d);
}
/**************************************-**************************************/
/*                                                                           */
/*                               MOD_SRF_RUN_01                              */
/*                                                                           */
/*****************************************************************************/
void mod_srf_run_01(m,s,r,action)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     int action;                // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  if (action == 0){
    mod_srf_01_prep(m,r);
  }else if (action == 1){
    mod_srf_01_get_response(m,s,r);
  }else if (action == -1){
    mod_srf_01_done(m);
  }
}
