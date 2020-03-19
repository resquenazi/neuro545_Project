/*****************************************************************************/
/*                                                                           */
/*   mod_lin_util.c                                                          */
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
#include "farray_util.h"
#include "carray_util.h"
#include "func_util.h"
#include "myrand_util.h"
#include "spike_util.h"
#include "paramfile_util.h"
#include "kernel_util.h"
#include "mod_util.h"
#include "ifc_util.h"
#include "mod.h" // Data structures
#include "paramfile.h"

// For cluster computation  (static needed w/ initialization
static int myid = -1;         // Process ID, -1-Single, 0-Master, >0-Slave
static char *mylogf = NULL;   // Logfile name
char ggstr[LONG_SLEN];        // Global log string

// Global for lin model
float *mod_lin_f = NULL;    // Linear filter
int mod_lin_fn;             // Length of filter

float mod_lin_tscale;
int   mod_lin_tn;

struct onode *mod_lin_sgo;   // Spike generation onode
char         *mod_lin_sgen;  // Type of spike generation, "poisson", ...

struct mod_lin_unit{
  char *name;
  float rate;         // Poisson rate for local source
  float prob;         // Prob of including spike from common input
  int delay;          // Delay of included spike (sampling units)
  int spread;         // Spread of included spike (sampling units)
  int burst_dur;      // Duration of burst (sampling units)
  float burst_prob;   // Probability of burst

  float *ps;          // Current Poisson spike train [pn], or NULL
  int    pn;          // Spike count in local Poisson train
  float *cs;          // Current spikes from common source [cn]
  int    cn;          // spike count
  float *ss;          // Current summed spike train [sn]
  int    sn;          // Spike count in summed train
};

int   mod_lin_unit_n;             // Number of units in 2nd layer
struct mod_lin_unit **mod_lin_u;  // List of second level units


/**************************************-**************************************/
/*                                                                           */
/*                              MOD_LIN_GET_UNIT                             */
/*                                                                           */
/*  Create the unit structure for this onode.                                */
/*                                                                           */
/*****************************************************************************/
struct mod_lin_unit *mod_lin_get_unit(uo)
     struct onode *uo;
{
  float d,s,bdur;
  struct mod_lin_unit *lu;

  lu = (struct mod_lin_unit *)myalloc(sizeof(struct mod_lin_unit));

  lu->name   = onode_getpar_chr_exit(uo,"name");
  lu->rate   = onode_getpar_flt_exit(uo,"rate");
  lu->prob   = onode_getpar_flt_exit(uo,"prob_c");
  d          = onode_getpar_flt_dflt(uo,"delay_c",0.0);
  s          = onode_getpar_flt_dflt(uo,"spread_c",0.0);
  bdur       = onode_getpar_flt_dflt(uo,"burst_dur",0.0);
  lu->burst_prob = onode_getpar_flt_dflt(uo,"burst_prob",0.0);

  lu->burst_dur = my_rint(bdur / mod_lin_tscale);
  lu->delay     = my_rint(   d / mod_lin_tscale);
  lu->spread    = my_rint(   s / mod_lin_tscale);

  printf("  Unit %s\n",lu->name);
  printf("    rate  %f\n",lu->rate);
  printf("    prob  %f\n",lu->prob);
  printf("    delay  %d\n",lu->delay);
  printf("    spread  %d\n",lu->spread);
  printf("    burst_dur  %d\n",lu->burst_dur);
  printf("    burst_prob  %f\n",lu->burst_prob);

  // Empty spike trains
  lu->ps = NULL;
  lu->pn = 0;
  lu->ss = NULL;
  lu->sn = 0;
  lu->cs = NULL;
  lu->cn = 0;

  return lu;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MOD_LIN_UNIT_I                              */
/*                                                                           */
/*  Return the index for the unit with 'name', or -1.                        */
/*                                                                           */
/*****************************************************************************/
int mod_lin_unit_i(name)
     char *name;
{
  int i,k;
  int n;
  struct mod_lin_unit **u;  // List of second level units
  
  n = mod_lin_unit_n;
  u = mod_lin_u;

  k = -1;
  i = 0;
  while((i < n) && (k == -1)){
    if (compare_prefix_string_order(u[i]->name,name) == 1){
      k = i;
    }else
      i += 1;
  }

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_LIN_01_PREP                              */
/*                                                                           */
/*  1.  Set parameters for the linear filter.                                */
/*  2.  Create the filter at 'mod_lin_f'.                                    */
/*                                                                           */
/*****************************************************************************/
void mod_lin_01_prep(m,r)
     struct model_struct *m;    // Model params
     struct response_struct *r; // Response params
{
  int i;
  char *sgen;
  struct onode *fo,*sgo,*t;

  mylog(mylogf,"  MOD_LIN_01_PREP\n");

  mod_lin_tn     = onode_getpar_int_exit(m->o,"tn");
  mod_lin_tscale = onode_getpar_flt_exit(m->o,"tscale");

  //
  //  Linear filter
  //
  fo = onode_child_get_unique(m->o,"filter");
  mod_lin_fn = -1;  // So that it gets set in call below
  mod_lin_f = kernu_o_filter(mylogf,fo,mod_lin_tscale,&mod_lin_fn);

  //
  //  Spike generation
  //
  mod_lin_sgo  = onode_child_get_unique(m->o,"spike_gen");
  mod_lin_sgen = onode_getpar_chr_exit(mod_lin_sgo,"type");

  //
  //  Test for 2nd level of units
  //
  mod_lin_unit_n = onode_count_otype(m->o,"unit");
  if (mod_lin_unit_n > 0){
    mod_lin_u = (struct mod_lin_unit **)myalloc(mod_lin_unit_n*
						sizeof(struct mod_lin_unit *));
    i = 0;
    t = onode_get_next_type(m->o->o,"unit");
    while(t != NULL){
      mod_lin_u[i] = mod_lin_get_unit(t);
      i += 1;
      t = onode_get_next_type(t->next,"unit");
    }
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_LIN_01_DONE                              */
/*                                                                           */
/*****************************************************************************/
void mod_lin_01_done(m)
     struct model_struct *m;        // Model parameter pair list
{
  mylog(mylogf,"  MOD_LIN_01_DONE\n");

  myfree(mod_lin_f);
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_LIN_01_GET_RESPONSE                         */
/*                                                                           */
/*****************************************************************************/
void mod_lin_01_get_response(m,s,r,k)
     struct model_struct *m;        // Model parameter pair list
     struct stim_struct *s;         // Stimulus data and param pair list
     struct response_struct *r;     // Response data
     int k;                         // Repeat number
{
  int i,j;
  int si,ti,tn,seed,*spk,ns,vn,n,*pseed,tseed,sdt,tsn,bdur;
  float tscale,samp,*conv,min,max,*fs,*v;
  float pr,fdt,*cs,*ts;
  struct mod_lin_unit *tu;

  tn     = mod_lin_tn;
  tscale = mod_lin_tscale;
  n      = mod_lin_unit_n;

  samp = 1.0/tscale;

  conv = convolve_with_mask_causal(s->d1,tn,mod_lin_f,mod_lin_fn);

  get_min_max_farray(conv,tn,&min,&max);
  sprintf(ggstr,"    min, max of raw conv:  %f %f\n",min,max);
  mylog(mylogf,ggstr);

  if (strcmp(mod_lin_sgen,"poisson")==0){
    seed = m->mseed[r->tsi];
    ifc_util_poisson(mylogf,m,mod_lin_sgo,conv,tn,samp,seed,0,1,&fs,&ns,&v,&vn);
    sprintf(ggstr,"    Made %d Poisson spikes\n",ns);
    mylog(mylogf,ggstr);
  }else
    mylogx(mylogf,"MOD_LIN_01_GET_RESPONSE","Spike generation type not imp'd");


  if (n > 0){

    // Create seeds for Poisson spike generation, and common spike picking
    pseed = get_seeds(m->mseed[r->tsi],100000,2*n);

    for(i=0;i<n;i++){

      tu = mod_lin_u[i];

      //  Free old spike trains
      if (tu->ps != NULL)
	myfree(tu->ps);
      if (tu->ss != NULL)
	myfree(tu->ss);
      if (tu->cs != NULL)
	myfree(tu->cs);

      //  Generate Poisson spikes
      tseed = pseed[i];
      make_poisson_float_spikes(&(tu->ps),&(tu->pn),0,tn,samp,tu->rate,0.0,
				0.0,1,&tseed);


      //  Pick spikes from Common input
      tseed = pseed[2*i];
      if (tseed > 0)
	tseed *= -1;

      pr  = tu->prob;
      fdt = tu->delay;

      cs = (float *)myalloc(ns*sizeof(float));
      si = 0;
      for(j=0;j<ns;j++){
	if (myrand_util_ran2(&tseed) < pr){

	  sdt = 0;  // Spread dt
	  if (tu->spread > 0) // Choose random amount for spread
	    sdt = my_rint((float)tu->spread * myrand_util_ran2(&tseed));

	  cs[si] = fs[j] + fdt + sdt;
	  si += 1;
	}
      }

      tu->cs = copy_farray(cs,si);
      tu->cn = si;
      myfree(cs);

      //
      //  Merge the independent Poisson spikes, 'ps' w/ common input 'cs'
      //
      ts = merge_spike_arrays_float(tu->ps,tu->pn,tu->cs,tu->cn);
      tsn = tu->cn + tu->pn;


      //
      //  Optional spike -> burst transformation
      //
      bdur = tu->burst_dur;
      if ((bdur > 1) && (tsn > 0)){  // If burst making and there are spikes
	pr = tu->burst_prob;
	cs = (float *)myalloc(tsn*bdur*sizeof(float));
	si = 0;
	for(j=0;j<tsn;j++){  // For each spike in merged train
	  for(ti=0;ti<bdur;ti++){  // For each spike in merged train
	    if (myrand_util_ran2(&tseed) < pr){
	      cs[si] = ts[j] + ti;
	      si += 1;
	    }
	  }
	}
	sort_farray(cs,si);
	myfree(ts);  // There must be at least 1 spike, ok to free
	ts = cs;
	tsn = si;
      }


      //
      //  Store final spikes
      //
      tu->ss = copy_farray(ts,tsn);
      tu->sn = tsn;
      if (ts != NULL)  // There could be zero spikes
	myfree(ts);
    }
    myfree(pseed);
  }


  for(i=0;i<r->n;i++){ // For each response requested
    if (strcmp(r->datid[i],"spikes")==0){
      mod_util_resp_store_s(r,i,fs,ns,1,mylogf);  // copyflag = 1
    }else if (strcmp(r->datid[i],"prob")==0){
      mod_util_resp_store_f(r,i,v,vn,1,mylogf);  // copyflag = 1
    }else if (strcmp(r->datid[i],"f")==0){
      mod_util_resp_store_f(r,i,conv,tn,1,mylogf);  // copyflag = 1
    }else if (mod_lin_unit_i(r->datid[i]) >=0){

      tu = mod_lin_u[mod_lin_unit_i(r->datid[i])];
      if (strcmp(r->datid[i],tu->name)==0){
	mod_util_resp_store_s(r,i,tu->ss,tu->sn,1,mylogf); // copyflag = 1
      }else{
	printf(" response:  %s\n",r->datid[i]);
	exit(0);
      }
    }
  }
  
  myfree(conv);
  if (fs != NULL)
    myfree(fs);
  myfree(v);
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_LIN_RUN_LIN_01                            */
/*                                                                           */
/*****************************************************************************/
void mod_lin_run_lin_01(m,s,r,k,action)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response params
     int k;                      // Trial number, when action = 1
     int action;                 // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  if (action == 0){
    mod_lin_01_prep(m,r);
  }else if (action == 1){
    mod_lin_01_get_response(m,s,r,k);
  }else if (action == -1){
    mod_lin_01_done(m);
  }
}
