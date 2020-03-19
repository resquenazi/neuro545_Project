/*****************************************************************************/
/*                                                                           */
/*  mod_test_util.c                                                          */
/*  wyeth bair                                                               */
/*                                                                           */
/*  Operate directly on the stimulus to produce spikes for testing           */
/*  analysis routines.                                                       */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "my_util.h"
#include "farray_util.h"
#include "myrand_util.h"
#include "paramfile_util.h"
#include "mod_util.h"
#include "mod.h"
#include "paramfile.h"

// For cluster computation  (static needed w/ initialization
static int myid = -1;         // Process ID, -1-Single, 0-Master, >0-Slave
static char *mylogf = NULL;   // Logfile name
char ggstr[LONG_SLEN];        // Global log string

float mod_test_sscale;
float mod_test_tscale;
int mod_test_xn;
int mod_test_yn;
int mod_test_tn;

int mod_test_eye;         // 0-left, 1-right, 2-left_to_right

int   mod_test_winx0;
int   mod_test_winy0;
int   mod_test_winxn;
int   mod_test_winyn;
int   mod_test_dx;
int   mod_test_dy;
int   mod_test_dt;
float mod_test_thresh_hi;
float mod_test_thresh_lo;


/**************************************-**************************************/
/*                                                                           */
/*                                MOD_TEST_PREP                              */
/*                                                                           */
/*****************************************************************************/
void mod_test_prep(m)
     struct model_struct *m; // Model params
{
  int xn,yn,tn;
  float sscale,tscale;
  char *modtype;
  struct param_pair_list *mppl; // Model parameter pair list

  myid = mylog_set_id_log(m->process_id,&mylogf);  // DO THIS FIRST

  mylog(mylogf,"  MOD_TEST_PREP\n");

  mppl = m->ppl;

  modtype = paramfile_get_char_param_or_exit(mppl,"mod_type");

  sscale = paramfile_get_float_param_or_exit(mppl,"sscale");
  tscale = paramfile_get_float_param_or_exit(mppl,"tscale");
  xn = paramfile_get_int_param_or_exit(mppl,"xn");
  yn = paramfile_get_int_param_or_exit(mppl,"yn");
  tn = paramfile_get_int_param_or_exit(mppl,"tn");

  mod_test_sscale = sscale;
  mod_test_tscale = tscale;
  mod_test_xn = xn;
  mod_test_yn = yn;
  mod_test_tn = tn;

  mod_test_eye = param_geti_dflt(mppl,"eye_flag",0);
  if ((mod_test_eye < 0) || (mod_test_eye > 2))
    mylog_exit(mylogf,"MOD_TEST_PREP  bad eye_flag value, use 0, 1 or 2\n");

  mod_test_winx0      = param_geti_exit(mppl,"win_x0");
  mod_test_winy0      = param_geti_exit(mppl,"win_y0");
  mod_test_winxn      = param_geti_exit(mppl,"win_xn");
  mod_test_winyn      = param_geti_exit(mppl,"win_yn");
  mod_test_dx         = param_geti_exit(mppl,"dx");
  mod_test_dy         = param_geti_exit(mppl,"dy");
  mod_test_dt         = param_geti_exit(mppl,"dt");
  mod_test_thresh_hi  = param_getf_exit(mppl,"thresh_hi");
  mod_test_thresh_lo  = param_getf_exit(mppl,"thresh_lo");

  myfree(modtype);
}
/**************************************-**************************************/
/*                                                                           */
/*                                MOD_TEST_DONE                              */
/*                                                                           */
/*****************************************************************************/
void mod_test_done(mppl)
     struct param_pair_list *mppl; // Model parameter pair list
{
  int xn,yn,tn;

  mylog(mylogf,"  MOD_TEST_DONE\n");
  
  xn = mod_test_xn;
  yn = mod_test_yn;
  tn = mod_test_tn;
}
/**************************************-**************************************/
/*                                                                           */
/*                        MOD_TEST_XYTCORR01_GET_RESPONSE                    */
/*                                                                           */
/*****************************************************************************/
void mod_test_xytcorr01_get_response(m,s,r)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response data
{
  int i;
  int xn,yn,tn,ns,ri,xi,xi1,xi2,yi,yi1,yi2;
  int winx0,winy0,winxn,winyn,dx,dy,dt,seed,dumpflag;
  float sampling,*fs,*td,*td_prev,thresh_hi,thresh_lo,prob,*corr11,*corr00;
  float *prod,pamp,poff,pr,*tp;
  float ***sd0,***sd1;  // Pointers to stimulus data
  char *comp;
  struct param_pair_list *mppl; // Model parameter pair list

  // GET PARAMS
  mppl = m->ppl;
  xn    = mod_test_xn;
  yn    = mod_test_yn;
  tn    = mod_test_tn;
  winx0 = mod_test_winx0;
  winy0 = mod_test_winy0;
  winxn = mod_test_winxn;
  winyn = mod_test_winyn;
  dx    = mod_test_dx;
  dy    = mod_test_dy;
  dt    = mod_test_dt;
  thresh_hi = mod_test_thresh_hi;
  thresh_lo = mod_test_thresh_lo;

  //seed = paramfile_get_int_param_or_exit(mppl,"spike_seed");
  prob = paramfile_get_float_param_or_exit(mppl,"spike_prob");
  pamp = paramfile_get_float_param_default(mppl,"pr_scale",1.0);
  poff = paramfile_get_float_param_default(mppl,"pr_offset",0.0);
  /*printf(" x,y  %d %d  dx,dy,dt  %d %d %d\n",x,y,dx,dy,dt);*/

  comp = paramfile_get_char_param_or_exit(mppl,"computation");

  dumpflag = paramfile_get_int_param_default(mppl,"dump_flag",0);

  seed = m->mseed[r->tsi];  // Trial seed

  if (seed > 0)
    seed *= -1;

  // MAKE SPIKES

  // Pick the middle coord if values are -1
  if (winx0 < 0)
    xi1 = xn/2 - winxn/2;
  else
    xi1 = winx0;
  if (winy0 < 0)
    yi1 = yn/2 - winyn/2;
  else
    yi1 = winy0;

  xi2 = xi1 + winxn - 1;
  yi2 = yi1 + winyn - 1;

  if ((xi1 < 0) || (xi2 >= xn))
    mylog_exit(mylogf,"MOD_TEST_XTYCORR01_GET_RESPONSE  x out of bounds");
  if ((yi1 < 0) || (yi2 >= yn))
    mylog_exit(mylogf,"MOD_TEST_XTYCORR01_GET_RESPONSE  y out of bounds");
  if ((dt < 0) || (dt >= (tn-1)))
    mylog_exit(mylogf,"MOD_TEST_XTYCORR01_GET_RESPONSE  dt out of bounds");

  if (((xi1+dx) < 0) || ((xi2+dx) >= xn))
    mylog_exit(mylogf,"MOD_TEST_XTYCORR01_GET_RESPONSE  dx out of bounds");
  if (((yi1+dy) < 0) || ((yi2+dy) >= yn))
    mylog_exit(mylogf,"MOD_TEST_XTYCORR01_GET_RESPONSE  dy out of bounds");

  /*
    printf("    Computing correlation across:\n");
    printf("      dx %d pixels\n",dx);
    printf("      dy %d pixels\n",dy);
    printf("      dt %d time units\n",dt);
    printf("    Starting at all points within:\n");
    printf("      x:  %d to %d\n",xi1,xi2);
    printf("      y:  %d to %d\n",yi1,yi2);
    printf("    White threshold  %f\n",thresh_hi);
    printf("    Black threshold  %f\n",thresh_lo);*/

  // Compute correlation as a function of time
  corr00 = get_zero_farray(tn);
  corr11 = get_zero_farray(tn);
  prod   = get_zero_farray(tn);

  if (mod_test_eye == 0){
    sd0 = sd1 = s->d;
  }else if (mod_test_eye == 1){
    sd0 = sd1 = s->d_r;
    if (s->d_r == NULL)
      mylog_exit(mylogf,"MOD_TEST_XTYCORR01_GET_RESPONSE  right stim NULL\n");
  }else if (mod_test_eye == 2){
    sd0 = s->d;
    sd1 = s->d_r;
    if (s->d_r == NULL)
      mylog_exit(mylogf,"MOD_TEST_XTYCORR01_GET_RESPONSE  right stim NULL\n");
  }

  /*
  printf("xi1 xi2  %d  %d\n",xi1,xi2);
  printf("yi1 yi2  %d  %d\n",yi1,yi2);
  */

  /*** WYETH TESTING ****/
  /*** WYETH TESTING ****/
  /*** WYETH TESTING ****/
  /*
  for(i=1;i<=50;i+=5)
    printf("%f\n",s->d[xi1][yi1][i]);
  exit(0);
  */

  for(xi=xi1;xi<=xi2;xi++){
    for(yi=yi1;yi<=yi2;yi++){

      //printf("xi = %d   yi = %d\n",xi,yi);

      td      = sd1[xi+1+dx][yi+1+dy];
      td_prev = sd0[xi+1   ][yi+1   ];

      /*printf(" TDOFF =   %d  %d\n",x+1+dx,y+1+dy);*/

      //for(i=dt+1;i<=tn;i++){
      for(i=dt;i<tn;i++){

	// Count correlated pixels
	if ((td[i+1] >= thresh_hi) && (td_prev[i+1-dt] >= thresh_hi)){
	  corr11[i] += 1.0;
	}
	if ((td[i+1] <= thresh_lo) && (td_prev[i+1-dt] <= thresh_lo)){
	  corr00[i] += 1.0;
	}

	// Integrate the product
	prod[i] +=  (td[i+1] - 0.5) * (td_prev[i+1-dt] - 0.5);
      }
    }
  }

  if (dumpflag == 1){
    append_farray_plot("zzz.dump.pl","product",prod,tn,1);
  }

  tp = get_farray(tn);  // Keep track of prob

  // Make spikes from correlation
  fs = get_zero_farray(tn);
  ns = 0;
  for(i=0;i<tn;i++){

    // Compute spiking prob, 'pr'
    if (strcmp(comp,"simple")==0){
      if ((corr00[i] > 0.0) || (corr11[i] > 0.0))
	pr = prob;
      else
	pr = 0.0;
    }else if (strcmp(comp,"multiply")==0){
      pr = pamp * (prod[i] + poff);
    }else
      exit_error("MOD_TEST_XYTCORR01_GET_RESPONSE","Unknown computation");

    if (pr < 0.0)
      pr = 0.0;
    if (pr > 1.0)
      pr = 1.0;

    tp[i] = pr;

    //  Generate a spike w/ prob 'pr'
    if (myrand_util_ran2(&seed) < pr){
      fs[ns] = (float)i * mod_test_tscale * 1000.0;  // ms time
      ns += 1;
    }
  }

  if (dumpflag == 1){
    append_farray_plot("zzz.dump.pl","Prob",tp,tn,1);
  }

  // SAVE SPIKES
  ri = 0;
  sampling = r->samp[ri];
  multiply_farray(fs,ns,sampling/1000.0);
  mod_util_resp_store_s(r,0,fs,ns,1,mylogf); // copyflag = 1

  printf("    %d spikes\n",ns);

  myfree(fs); myfree(corr00); myfree(corr11);
  myfree(prod);
  myfree(tp);
  myfree(comp);
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_TEST_RUN_XYTCORR01                          */
/*                                                                           */
/*  Look for correlation over space and time in the stimulus.                */
/*                                                                           */
/*****************************************************************************/
void mod_test_run_xytcorr01(m,s,r,action)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     int action;                // -1-cleanup, 0-prep, 1-run
{
  char *sform;

  myid = mylog_set_id_log(m->process_id,&mylogf);  // DO THIS FIRST

  if (action == 0){
    mod_test_prep(m);
  }else if (action == 1){
    mod_test_xytcorr01_get_response(m,s,r);
  }else if (action == -1){
    mod_test_done(m->ppl);
  }
}
