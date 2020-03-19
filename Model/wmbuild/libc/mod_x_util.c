/*****************************************************************************/
/*                                                                           */
/*   mod_x_util.c                                                            */
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
#include "iarray_util.h"
#include "farray_util.h"
#include "spike_util.h"
#include "paramfile_util.h"
#include "mod_util.h"
#include "mod.h" // Data structures
#include "paramfile.h"

// For cluster computation  (static needed w/ initialization)
static int myid = -1;         // Process ID, -1-Single, 0-Master, >0-Slave
static char *mylogf = NULL;   // Logfile name
char ggstr[LONG_SLEN];        // Global log string


// Global variables for model x
int mod_x_xn;
int mod_x_yn;
int mod_x_tn;
float mod_x_tscale;
int mod_x_xi;
int mod_x_yi;
int mod_x_toff;
int mod_x_eye;    // 0-left, 1-right

/**************************************-**************************************/
/*                                                                           */
/*                                MOD_X_01_PREP                              */
/*                                                                           */
/*  Prepare constructs that will be needed when running trials.  Do things   */
/*  here that only need to be done once, to avoid repeating them for each    */
/*  trial.                                                                   */
/*                                                                           */
/*****************************************************************************/
void mod_x_01_prep(m,r)
     struct model_struct *m;    // Model params
     struct response_struct *r; // Response params
{
  int xn,yn,tn,xi,yi,toff;
  float tscale;

  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  mylog(mylogf,"  MOD_X_01_PREP\n");

  xn     = onode_getpar_int_exit(m->o,"xn");
  yn     = onode_getpar_int_exit(m->o,"yn");
  tn     = onode_getpar_int_exit(m->o,"tn");
  tscale = onode_getpar_flt_exit(m->o,"tscale");

  // Coordinates within 0..xn-1  0..yn-1
  xi = onode_getpar_int_exit(m->o,"xi");
  yi = onode_getpar_int_exit(m->o,"yi");
  toff = onode_getpar_int_dflt(m->o,"toff",0);  // msec
  mod_x_eye = onode_getpar_int_dflt(m->o,"eye_flag",0);  // 0-left, 1-right

  // Store in global variables for usage during run
  mod_x_xn = xn;
  mod_x_yn = yn;
  mod_x_tn = tn;
  mod_x_tscale = tscale;

  mod_x_xi = xi;
  mod_x_yi = yi;
  mod_x_toff = toff;
}
/**************************************-**************************************/
/*                                                                           */
/*                                MOD_X_01_DONE                              */
/*                                                                           */
/*  Free any storage here when done.                                         */
/*                                                                           */
/*****************************************************************************/
void mod_x_01_done(m)
     struct model_struct *m;        // Model parameters
{
  mylog(mylogf,"  MOD_X_01_DONE\n");

  // Nothing to do here.
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_X_01_GET_RESPONSE                          */
/*                                                                           */
/*****************************************************************************/
void mod_x_01_get_response(m,s,r)
     struct model_struct *m;        // Model parameters
     struct stim_struct *s;         // Stimulus parameters
     struct response_struct *r;     // Response parameters
{
  int i;
  int xi,yi,ss,pn;
  int tn,flag,seed,*spk,cnt,cpflag;
  float *d,*prob,mu,tscale,*fs;

  // Use convenient local variables to access global parameters
  tscale = mod_x_tscale;
  tn = mod_x_tn;
  xi = mod_x_xi;
  yi = mod_x_yi;

  // Extract the stimulus vs. time at coordinate (xi,yi)
  d = (float *)myalloc(tn*sizeof(float));

  if (mod_x_eye == 0){  // Use left eye stimulus
    for(i=0;i<tn;i++)
      d[i] = s->d[1+xi][1+yi][i+1];
  }else if (mod_x_eye == 1){  // Use right eye stimulus
    if (s->d_r == NULL)
      mylog_exit(mylogf,"MOD_X_01_GET_RESPONSE  Right eye stimulus is NULL\n");

    mylog(mylogf,"    Using Right Eye stimulus\n");
    for(i=0;i<tn;i++)
      d[i] = s->d_r[1+xi][1+yi][i+1];
  }

  // Write to a dummy file for testing and debugging
  //OLD REMOVE if (r->tsi == 0) // If this is the first trial
  if (r->gtsi == 0) // If this is the first trial
    remove_file("zz.mod_x.pl");

  append_farray_plot("zz.mod_x.pl","raw_stim",d,tn,1);

  // Rescale 'd' to have a mean value appropriate for spike generation
  mu = mean_farray(d,tn);
  multiply_farray(d,tn,0.1/mu);  // Make the mean be 0.1

  // Create array 'prob', the firing probability per time bin at the
  // temporal resolution of 1 msec bins.
  if (tscale > 0.001){
    ss = tscale/0.001;  // Turn each sample point into 'ss' points
    prob = over_sample_farray(d,tn,ss);
    pn = tn*ss;
  }else if (tscale == 0.001){
    prob = copy_farray(d,tn);
    pn = tn;
  }else{
    mylog_exit(mylogf,"MOD_X_01_GET_RESPONSE  tscale < 1ms\n");
  }

  // Append the probability array to a dummy file for testing/debugging
  append_farray_plot("zz.mod_x.pl","prob",prob,pn,1);

  // Create Poisson spikes
  seed = m->mseed[r->tsi];  // Get the randomization seed for current trial
  make_poisson_spikes_from_prob(prob,pn,1000.0,&seed,&spk,&cnt);
  // 'spk' is a list of spike times, 'cnt' is the number of spikes

  add_const_iarray(spk,cnt,mod_x_toff);

  fs = i2farray(spk,cnt);  // Convert spike times to floats
  cpflag = 0;              // Do not make copy, 'fs' will be linked to 'r'
  mod_util_resp_store_s(r,0,fs,cnt,cpflag,mylogf);

  // Free temporary storage.
  // Note, 'fs' should not be freed, it is linked from 'r'
  myfree(spk);
  myfree(prob);
  myfree(d);
}
/**************************************-**************************************/
/*                                                                           */
/*                                MOD_X_RUN_01                               */
/*                                                                           */
/*****************************************************************************/
void mod_x_run_01(m,s,r,action)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     int action;                // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  if (action == 0){
    mod_x_01_prep(m,r);
  }else if (action == 1){
    mod_x_01_get_response(m,s,r);
  }else if (action == -1){
    mod_x_01_done(m);
  }
}
