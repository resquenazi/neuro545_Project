/*****************************************************************************/
/*                                                                           */
/*  ifc_util.c                                                               */
/*  wyeth bair                                                               */
/*                                                                           */
/*  Integrate and fire conductance based simulations.                        */
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
#include "plot_util.h"
#include "func_util.h"
#include "myrand_util.h"
#include "spike_util.h"
#include "paramfile_util.h"
#include "sig_util.h"
#include "mod_util.h"
#include "pop_util.h"
#include "mod.h"
#include "ifc.h"
#include "paramfile.h"

#define NRANSI
#define NR_END 1
#define FREE_ARG char*
static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define HMAX 0.0005       // Maximum step for Runge-Kutta (sec)

#define MAX_SPK 20000   // Maximum number of spikes in temporary storage

#define DUMP_FNAME  "zz.dump.spike.pl"  // File for spike gen dumping

//  For cluster computation  (static needed w/ initialization
static int myid = -1;         // Process ID, -1-Single, 0-Master, >0-Slave
static char *mylogf = NULL;   // Logfile name
char ggstr[LONG_SLEN];        // Global log string
     
//  Global variables
int   ifc_gp_dur;         // duration of gtx, gti, gta in 'samp' units
float ifc_gp_samp;        // Sampling for g_ex, g_in
float ifc_gp_v_spike;     // Voltage used to represent a spike
float ifc_gp_v_th_x;      // Spike threshold, ex

float *ifc_gp_gtx;        // g_ex(t)
float *ifc_gp_gti;        // g_in(t)
float *ifc_gp_gta;        // g_ad(t)
float *ifc_gp_gta_pulse;  // Single adaptation g pulse
int   ifc_gp_gta_npulse;  // number of points in pulse
float ifc_gp_tau_r_ad;    // tau_rise, adapt
float ifc_gp_tau_f_ad;    // tau_fall, adapt
float ifc_gp_gbar_ad;

float ifc_gp_gx_scale;    // Scale excitatory conductance (before adding bias)
float ifc_gp_gx_bias;     // Add bias to exc. cond. (after andy scaling)
float ifc_gp_gi_scale;    // Scale inhibitory conductance (before adding bias)
float ifc_gp_gi_bias;     // Add bias to inhib. cond. (after andy scaling)

int ifc_gp_grect;         // 1 - half-wave rectify gx, gi after adding noise

float ifc_gp_trefr_x;     // ex refractory time (ms)
float ifc_gp_trefr_x_sd;  // SD of Gaussian refraction time (ms)
float ifc_gp_v_reset_x;   // ex reset potential (mV)
float ifc_gp_v_ex;        // ex reversal potential (mV)
float ifc_gp_v_in;        // in reversal potential (mV)
float ifc_gp_v_leak_x;    // leak reversal potential (mV)
float ifc_gp_v_ad;        // adapt reversal potential (mV)
float ifc_gp_g_leak_x;    // leakage conductance (nS)
float ifc_gp_c_x;         // membrane cap. (pF)

/*** Global storage for ODEINT intermediate results ***/
/* Preset 'kmax' and 'dxsav' in calling program.  If kmax != 0, results are
   stored at approximate intervals 'dxsav' in the arrays 'xp' and 'yp'
   where 'kount' is output by ODEINT.  Defining declarations for these
   variables, with memory allocation for 'xp' and 'yp' should be in
   calling program. */

int kmax,kount;
float *xp,**yp;  /* xp[1..kount], yp[1..nvar][1..kount] */
float dxsav;

/*** Global storage for spike times ***/
float gspk[MAX_SPK];
int gspkn;

int ifc_odeseed;  // Seed for random variations in ODEINT

//
//  IFC_01
//
int     ifc_01_tn;     // tn
float   ifc_01_samp;   // 1.0 / tscale

char   *ifc_01_stim_use;     // How to use the 1-D stimulus input

float   ifc_01_sigsm;
float   ifc_01_tdelay;       // Time delay (s) added to spikes

float  *ifc_01_ex_mask;
float  *ifc_01_in_mask;
int     ifc_01_ex_mask_n;
int     ifc_01_in_mask_n;

float   ifc_01_bg_ex_rate;
float   ifc_01_bg_in_rate;
int    *ifc_01_bg_ex_seed;
float   ifc_01_bg_ex_amp;
float   ifc_01_bg_in_amp;
int    *ifc_01_bg_in_seed;

float   ifc_01_pre_ex_rate;
float   ifc_01_pre_ex_amp;
int    *ifc_01_pre_ex_seed;

float   ifc_01_pre_in_rate;
float   ifc_01_pre_in_amp;
int    *ifc_01_pre_in_seed;



/**************************************-**************************************/
/*                                                                           */
/*                           IFC_UTIL_RM_DUMP_FILE                           */
/*                                                                           */
/*****************************************************************************/
void ifc_util_rm_dump_file(mmpid)
     int mmpid;
{
  if (mmpid == -1){
    remove_file(DUMP_FNAME);  // Don't do this in mpi mode
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           IFC_UTIL_FREE_IFC_PARAM                         */
/*                                                                           */
/*****************************************************************************/
void ifc_util_free_ifc_param(ifp)
     struct ifc_param *ifp;
{
  if (ifp->gx_noise != (char *)NULL)
    myfree(ifp->gx_noise);

  if (ifp->gi_noise != (char *)NULL)
    myfree(ifp->gi_noise);

  myfree(ifp);
}
/**************************************-**************************************/
/*                                                                           */
/*                              IFC_UTIL_GET_PARAMS                          */
/*                                                                           */
/*****************************************************************************/
struct ifc_param *ifc_util_get_params(ppl,name)
     struct param_pair_list *ppl;
     char name[];
{
  char pname[SLEN];
  struct ifc_param *ifc;

  ifc = (struct ifc_param *)myalloc(sizeof(struct ifc_param));

  sprintf(pname,"%s_v_spike",name);
  ifc->v_spike   = paramfile_get_float_param_or_exit(ppl,pname);
  sprintf(pname,"%s_tau_r_ad",name);
  ifc->tau_r_ad  = paramfile_get_float_param_or_exit(ppl,pname);
  sprintf(pname,"%s_tau_f_ad",name);
  ifc->tau_f_ad  = paramfile_get_float_param_or_exit(ppl,pname);
  sprintf(pname,"%s_gbar_ad",name);
  ifc->gbar_ad   = paramfile_get_float_param_or_exit(ppl,pname);
  sprintf(pname,"%s_v_reset_x",name);
  ifc->v_reset_x = paramfile_get_float_param_or_exit(ppl,pname);
  sprintf(pname,"%s_v_ex",name);
  ifc->v_ex      = paramfile_get_float_param_or_exit(ppl,pname);
  sprintf(pname,"%s_v_in",name);
  ifc->v_in      = paramfile_get_float_param_or_exit(ppl,pname);
  sprintf(pname,"%s_v_ad",name);
  ifc->v_ad      = paramfile_get_float_param_or_exit(ppl,pname);
  sprintf(pname,"%s_v_th_x",name);
  ifc->v_th_x    = paramfile_get_float_param_or_exit(ppl,pname);
  sprintf(pname,"%s_v_leak_x",name);
  ifc->v_leak_x  = paramfile_get_float_param_or_exit(ppl,pname);
  sprintf(pname,"%s_g_leak_x",name);
  ifc->g_leak_x  = paramfile_get_float_param_or_exit(ppl,pname);
  sprintf(pname,"%s_c_x",name);
  ifc->c_x       = paramfile_get_float_param_or_exit(ppl,pname);
  sprintf(pname,"%s_trefr_x",name);
  ifc->trefr_x   = paramfile_get_float_param_or_exit(ppl,pname);
  ifc->trefr_x_s = ifc->trefr_x/1000.0;
  sprintf(pname,"%s_trefr_x_sd",name);
  ifc->trefr_x_sd = paramfile_get_float_param_default(ppl,pname,2.0);

  sprintf(pname,"%s_gx_scale",name);
  ifc->gx_scale  = paramfile_get_float_param_or_exit(ppl,pname);
  sprintf(pname,"%s_gx_bias",name);
  ifc->gx_bias   = paramfile_get_float_param_or_exit(ppl,pname);
  sprintf(pname,"%s_gi_scale",name);
  ifc->gi_scale  = paramfile_get_float_param_or_exit(ppl,pname);
  sprintf(pname,"%s_gi_bias",name);
  ifc->gi_bias   = paramfile_get_float_param_or_exit(ppl,pname);

  // WYETH - should change to 'spike_dump' only
  sprintf(pname,"%s_dump",name);
  ifc->dump      = paramfile_get_int_param_default(ppl,pname,0);
  if (ifc->dump == 0)
    ifc->dump      = paramfile_get_int_param_default(ppl,"spike_dump",0);

  sprintf(pname,"%s_clip_spike",name);
  ifc->clip_spike= paramfile_get_int_param_default(ppl,pname,0);

  sprintf(pname,"%s_gx_noise",name);
  if (paramfile_test_param(ppl,pname)){
    ifc->gx_noise = paramfile_get_nth_char_param_or_exit(ppl,pname,0);
    if (strcmp(ifc->gx_noise,"gfg")==0){
      ifc->gx_noise_mu = paramfile_get_nth_float_param_or_exit(ppl,pname,1);
      ifc->gx_noise_sd = paramfile_get_nth_float_param_or_exit(ppl,pname,2);
      ifc->gx_noise_tsd = paramfile_get_nth_float_param_or_exit(ppl,pname,3);
    }
  }else
    ifc->gx_noise = (char *)NULL;

  sprintf(pname,"%s_gi_noise",name);
  if (paramfile_test_param(ppl,pname)){
    ifc->gi_noise = paramfile_get_nth_char_param_or_exit(ppl,pname,0);
    if (strcmp(ifc->gi_noise,"gfg")==0){
      ifc->gi_noise_mu = paramfile_get_nth_float_param_or_exit(ppl,pname,1);
      ifc->gi_noise_sd = paramfile_get_nth_float_param_or_exit(ppl,pname,2);
      ifc->gi_noise_tsd = paramfile_get_nth_float_param_or_exit(ppl,pname,3);
    }
  }else
    ifc->gi_noise = (char *)NULL;

  sprintf(pname,"%s_grect",name);
  ifc->grect = paramfile_get_int_param_default(ppl,pname,1);
  /*  WYETH THIS IS CHANGED TO RECTIFY BY DEFAULT
      if (paramfile_test_param(ppl,pname)){
      ifc->grect = paramfile_get_int_param_or_exit(ppl,pname);
      }else
      ifc->grect = 0;*/


  /*** 2-compartment ***/

  sprintf(pname,"%s_g_tran",name);
  ifc->g_tran = paramfile_get_float_param_default(ppl,pname,0.0);
  if (ifc->g_tran != 0){
    ; /*printf("  Transfer conductance is > 0.0\n");*/
    /*** Read other 2-compartment params ***/
  }

  return ifc;
}
/**************************************-**************************************/
/*                                                                           */
/*                            IFC_UTIL_GET_PARAM_O                           */
/*                                                                           */
/*  New way - onode.                                                         */
/*                                                                           */
/*****************************************************************************/
struct ifc_param *ifc_util_get_param_o(o)
     struct onode *o;
{
  struct ifc_param *ifc;
  struct onode *no;
  /*char *noise_type;*/

  ifc = (struct ifc_param *)myalloc(sizeof(struct ifc_param));

  ifc->v_spike    = onode_getpar_flt_exit(o,"v_spike");
  //
  //  WYETH
  //  1.  Make gbar_ad a dflt w/ 0
  //  2.  if gbar_ad != 0, only then read tau's.
  //  3.  v_ad is then optional also
  //
  ifc->tau_r_ad   = onode_getpar_flt_exit(o,"tau_r_ad");
  ifc->tau_f_ad   = onode_getpar_flt_exit(o,"tau_f_ad");
  ifc->gbar_ad    = onode_getpar_flt_exit(o,"gbar_ad");
  ifc->v_reset_x  = onode_getpar_flt_exit(o,"v_reset_x");
  ifc->v_ex       = onode_getpar_flt_exit(o,"v_ex");
  ifc->v_in       = onode_getpar_flt_exit(o,"v_in");
  ifc->v_ad       = onode_getpar_flt_exit(o,"v_ad");
  ifc->v_th_x     = onode_getpar_flt_exit(o,"v_th_x");
  ifc->v_leak_x   = onode_getpar_flt_exit(o,"v_leak_x");
  ifc->g_leak_x   = onode_getpar_flt_exit(o,"g_leak_x");
  ifc->c_x        = onode_getpar_flt_exit(o,"c_x");
  ifc->trefr_x    = onode_getpar_flt_exit(o,"trefr_x");
  ifc->trefr_x_s  = ifc->trefr_x/1000.0;
  ifc->trefr_x_sd = onode_getpar_flt_dflt(o,"trefr_x_sd",2.0);

  ifc->gx_scale   = onode_getpar_flt_dflt(o,"gx_scale",1.0);
  ifc->gx_bias    = onode_getpar_flt_dflt(o,"gx_bias",0.0);
  ifc->gi_scale   = onode_getpar_flt_dflt(o,"gi_scale",1.0);
  ifc->gi_bias    = onode_getpar_flt_dflt(o,"gi_bias",0.0);

  ifc->grect      = onode_getpar_int_dflt(o,"grect",1);
  ifc->dump       = onode_getpar_int_dflt(o,"dump",0);
  ifc->clip_spike = onode_getpar_int_dflt(o,"clip_spike",0);

  /*** 2-compartment ***/
  ifc->g_tran     = onode_getpar_flt_dflt(o,"g_tran",0.0);

  if (onode_count_otype(o,"gx_noise")==1){
    no = onode_child_get_unique(o,"gx_noise");
    /*noise_type = onode_getpar_chr_exit(no,"type");*/
    ifc->gx_noise = onode_getpar_chr_exit(no,"type");
    if (strcmp(ifc->gx_noise,"gfg")==0){
      /*** WYETH - THESE are used from here for cortical pop's but
	   ifc_test for LGN will get params from onode again */
      ifc->gx_noise_mu  = onode_getpar_flt_dflt(no,"mean",0.0);
      ifc->gx_noise_sd  = onode_getpar_flt_dflt(no,"sd",0.0);
/*printf("ifc_util WYETH SD = %f\n",ifc->gx_noise_sd);*/
      ifc->gx_noise_tsd = onode_getpar_flt_dflt(no,"tsd",0.0);
    }
    /*myfree(noise_type);*/
  }else
    ifc->gx_noise = (char *)NULL;

  if (onode_count_otype(o,"gi_noise")==1){
    no = onode_child_get_unique(o,"gi_noise");
    /*noise_type = onode_getpar_chr_exit(no,"type");*/
    ifc->gi_noise = onode_getpar_chr_exit(no,"type");
    if (strcmp(ifc->gi_noise,"gfg")==0){
      ifc->gi_noise_mu  = onode_getpar_flt_dflt(no,"mean",0.0);
      ifc->gi_noise_sd  = onode_getpar_flt_dflt(no,"sd",0.0);
      ifc->gi_noise_tsd = onode_getpar_flt_dflt(no,"tsd",0.0);
    }
    /*myfree(noise_type);*/
  }else
    ifc->gi_noise = (char *)NULL;

  return ifc;
}
/**************************************-**************************************/
/*                                                                           */
/*                          IFC_UTIL_GET_PARAM_POISS_O                       */
/*                                                                           */
/*  Get parameters that might be relevant for a Poisson spiking unit.        */
/*                                                                           */
/*****************************************************************************/
struct ifc_param *ifc_util_get_param_poiss_o(o)
     struct onode *o;
{
  struct ifc_param *ifc;
  struct onode *no;

  ifc = (struct ifc_param *)myalloc(sizeof(struct ifc_param));

  ifc->v_spike    = 0.0;
  ifc->tau_r_ad   = onode_getpar_flt_dflt(o,"tau_r_ad",0.0);
  ifc->tau_f_ad   = onode_getpar_flt_dflt(o,"tau_f_ad",0.0);
  ifc->gbar_ad    = onode_getpar_flt_dflt(o,"gbar_ad",0.0);
  ifc->v_reset_x  = 0.0;
  ifc->v_ex       = 0.0;
  ifc->v_in       = 0.0;
  ifc->v_ad       = 0.0;
  ifc->v_th_x     = 0.0;
  ifc->v_leak_x   = 0.0;
  ifc->g_leak_x   = 0.0;
  ifc->c_x        = 0.0;
  ifc->trefr_x    = onode_getpar_flt_dflt(o,"trefr_x",0.0);
  ifc->trefr_x_s  = ifc->trefr_x/1000.0;
  ifc->trefr_x_sd = onode_getpar_flt_dflt(o,"trefr_x_sd",2.0);
  ifc->gx_scale   = onode_getpar_flt_dflt(o,"gx_scale",1.0);
  ifc->gx_bias    = onode_getpar_flt_dflt(o,"gx_bias",0.0);
  ifc->gi_scale   = onode_getpar_flt_dflt(o,"gi_scale",1.0);
  ifc->gi_bias    = onode_getpar_flt_dflt(o,"gi_bias",0.0);
  ifc->grect      = onode_getpar_int_dflt(o,"grect",1);
  ifc->dump       = onode_getpar_int_dflt(o,"dump",0);
  ifc->clip_spike = onode_getpar_int_dflt(o,"clip_spike",0);
  ifc->g_tran     = 0.0; // 2-compartment

  if (onode_count_otype(o,"gx_noise")==1){
    no = onode_child_get_unique(o,"gx_noise");
    ifc->gx_noise = onode_getpar_chr_exit(no,"type");
    if (strcmp(ifc->gx_noise,"gfg")==0){
      /*** WYETH - THESE are used from here for cortical pop's but
	   ifc_test for LGN will get params from onode again */
      ifc->gx_noise_mu  = onode_getpar_flt_dflt(no,"mean",0.0);
      ifc->gx_noise_sd  = onode_getpar_flt_dflt(no,"sd",0.0);
      ifc->gx_noise_tsd = onode_getpar_flt_dflt(no,"tsd",0.0);
    }
  }else
    ifc->gx_noise = (char *)NULL;

  if (onode_count_otype(o,"gi_noise")==1){
    no = onode_child_get_unique(o,"gi_noise");
    ifc->gi_noise = onode_getpar_chr_exit(no,"type");
    if (strcmp(ifc->gi_noise,"gfg")==0){
      ifc->gi_noise_mu  = onode_getpar_flt_dflt(no,"mean",0.0);
      ifc->gi_noise_sd  = onode_getpar_flt_dflt(no,"sd",0.0);
      ifc->gi_noise_tsd = onode_getpar_flt_dflt(no,"tsd",0.0);
    }
  }else
    ifc->gi_noise = (char *)NULL;

  return ifc;
}
/**************************************-**************************************/
/*                                                                           */
/*                          IFC_UTIL_GET_PARAM_POISSON                       */
/*                                                                           */
/*  For Poisson spike_gen, get params that are not used by IFC.              */
/*                                                                           */
/*****************************************************************************/
struct poiss_param *ifc_util_get_param_poisson(o)
     struct onode *o;
{
  struct poiss_param *pp;

  pp = (struct poiss_param *)myalloc(sizeof(struct poiss_param));

  pp->style    = onode_getpar_chr_dflt(o,"rate_rule","default_g");
  pp->scale_in = onode_getpar_flt_dflt(o,"scale_in",1.0);
  pp->scale_ad = onode_getpar_flt_dflt(o,"scale_ad",1.0);
  pp->offset0  = onode_getpar_flt_dflt(o,"offset0" ,0.0);
  pp->scale    = onode_getpar_flt_dflt(o,"scale"   ,1.0);
  pp->offset1  = onode_getpar_flt_dflt(o,"offset"  ,0.0);

  return pp;
}
/**************************************-**************************************/
/*                                                                           */
/*                                  IFC_RKCK                                 */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void ifc_rkck(float y[], float dydx[], int n, float x, float h, float yout[],
	      float yerr[], void (*derivs)(float, float [], float []))
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
  (*derivs)(x+a2*h,ytemp,ak2);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
  (*derivs)(x+a3*h,ytemp,ak3);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
  (*derivs)(x+a4*h,ytemp,ak4);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  (*derivs)(x+a5*h,ytemp,ak5);
  for (i=1;i<=n;i++)
    ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  (*derivs)(x+a6*h,ytemp,ak6);
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
/*                                  IFC_RKQS                                 */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void ifc_rkqs(float y[], float dydx[], int n, float *x, float htry, float eps,
	      float yscal[], float *hdid, float *hnext,
	      void (*derivs)(float, float [], float []))
{
  void ifc_rkck(float y[], float dydx[], int n, float x, float h, float yout[],
		float yerr[], void (*derivs)(float, float [], float []));
  int i;
  float errmax,h,xnew,*yerr,*ytemp;

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
    ifc_rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
    /*printf("  %f %f\n",ytemp[1],yerr[1]);*/
    errmax=0.0;
    for (i=1;i<=n;i++)
      errmax = FMAX(errmax,fabs(yerr[i]/yscal[i]));
    errmax /= eps;
    if (errmax > 1.0){
      h = SAFETY*h*pow(errmax,PSHRNK);
      if (h < 0.1*h){  // Only if 'h' is negative and integrating in reverse
	h *= 0.1;
	printf("HERE h = %f\n",h);
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
/*************************************---*************************************/
#define MAXSTP 100000 // changed from 10,000, but should be (x2-x1)*2000
#define TINY 1.0e-30
/**************************************-**************************************/
/*                                                                           */
/*                                 IFC_ODEINT                                */
/*                                                                           */
/*  ystart[nvar] - initial values, REPLACED BY VALUES AT END OF INTEGRATION  */
/*  nvar         - number of equations                                       */
/*  x1,x2        - solve in this interval (sec)                              */
/*  eps          - accuracy                                                  */
/*  h1           - guessed first stepsize (sec)                              */
/*  hmin         - minimum allowable stepsize (sec)                          */
/*  nok          - number of good steps                                      */
/*  nbad         - number of bad (but retried and fixed) steps               */
/*  derivs       - user supplied routine for right-hand side derivatives     */
/*  rkqs         - name of stepper routine to be used                        */
/*                                                                           */
/*   Modified from:                                                          */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void ifc_odeint(float ystart[], int nvar, float x1, float x2, float eps, 
		float h1, float hmin, int *nok, int *nbad,
		void (*derivs)(float, float [], float []),
		void (*rkqs)(float [], float [], int, float *, float, float, 
			     float [], float *, float *, 
			     void (*)(float, float [], float [])))
{
  int nstp,i,j,sflag,toffset;
  float xsav,x,hnext,hdid,h;
  float *yscal,*y,*dydx,dreset,d2reset,rv;
  int mymaxstep;

  yscal = vector(1,nvar);
  y     = vector(1,nvar);
  dydx  = vector(1,nvar);
  x = x1;
  h = SIGN(h1,x2-x1);
  *nok = (*nbad) = kount = 0;

  /* WYETH - this was 2000, Jun 2006, but then we added the HMAX constraint */
  /*mymaxstep = (x2-x1)*3000;*/
  mymaxstep = (x2-x1)*6000; /* Doubled for Nicolas to run Adam's stim ?? */
  /* printf("  mymaxstep = %d\n",mymaxstep);*/

  /* half of the range for a uniform reset voltage */
  dreset = (ifc_gp_v_th_x - ifc_gp_v_reset_x - 0.5);
  d2reset = 2.0 * dreset;

  for (i=1;i<=nvar;i++)
    y[i] = ystart[i];
  if (kmax > 0) /* If we're storing results */
    xsav = x-dxsav*2.0; /* Assures storage of first step */
  else
    xsav = 0.0;

  /*for (nstp=1;nstp<=MAXSTP;nstp++){*/
  for (nstp=1;nstp<=mymaxstep;nstp++){
    (*derivs)(x,y,dydx); /* Fill dydx, x and y are unchanged */

    /*printf("  %4d  x,y,dydx  %f %f %f\n",nstp,x,y[1],dydx[1]);*/

    for (i=1;i<=nvar;i++)
      /* Scaling used to monitor accuracy, this general-purpose choice can
	 be modified if need be */
      yscal[i] = fabs(y[i])+fabs(dydx[i]*h)+TINY;
    if (kmax > 0 && (kount < kmax-1) && fabs(x-xsav) > fabs(dxsav)){
      xp[++kount] = x;  /* Store intermediate result */
      for (i=1;i<=nvar;i++)
	yp[i][kount] = y[i];
      xsav = x;
    }
    if ((x+h-x2)*(x+h-x1) > 0.0) /* If stepsize can overshoot, decrease */
      h = x2-x;
    (*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs); /* Set new y */

    /*  This wants to take large step sizes, but it skips over
	threshold by many milliseconds.  Also, it jumps too far given
	that presynaptic spikes may not yet have occurred. */
    /* WYETH, use the following to enforce small steps: */

    if (hnext > HMAX)
      hnext = HMAX;

    if (hdid == h){
      ++(*nok);
      /*printf("  h= %f  (%f %f)\n",h,x1,x2);*/
    }else{
      ++(*nbad);
      /*printf("  BAD:  h= %f   did= %f\n",h,hdid);*/
    }
    
    /*** Check for spikes and reset state ************************************/
    sflag = 0;
    for (i=1;i<=nvar;i++)
      if (y[i] >= ifc_gp_v_th_x)
	sflag = 1;

    if (sflag){
      /*printf("SPIKE\n");*/
      toffset = (int)(x*ifc_gp_samp);
      if (gspkn >= MAX_SPK)
	exit_error("IFC_ODEINT","Max number of spikes exceeded");
      gspk[gspkn] = x*1000.0;  // Store time in msec
      gspkn += 1;
      
      if (kmax > 0 && (kount >= kmax-1))
	exit_error("IFC_ODEINT","kmax exceeded");
      xp[++kount] = x;
      for (i=1;i<=nvar;i++)
	yp[i][kount] = y[i]; // Save all current values
      xsav = x;
      
      if (kmax > 0 && (kount >= kmax-1))
	exit_error("IFC_ODEINT","kmax exceeded");
      xp[++kount] = x + dxsav/2.0;
      for (i=1;i<=nvar;i++){
	if (y[i] >= ifc_gp_v_th_x){
	  yp[i][kount] = ifc_gp_v_spike;     // Spike

	  /* Choose a random reset V */
	  /*
	    vre = d2reset * myrand_util_ran2(&ifc_odeseed) - dreset;
	    vre += ifc_gp_v_reset_x;
	    if (vre > ifc_gp_v_th_x - 0.5)
	    vre = ifc_gp_v_th_x - 0.5;
	    y[i] = vre;*/
	  
	  y[i] = ifc_gp_v_reset_x;          /* Reset Voltage */

	  for(j=0;j<ifc_gp_gta_npulse;j++){  /* Adaptation pulse */
	    /*printf("j=%d\n",j);*/
	    if ((j+toffset) < ifc_gp_dur)
	      ifc_gp_gta[j+toffset] += ifc_gp_gta_pulse[j];
	  }
	}else
	  yp[i][kount] = yp[i][kount-1]; // No spike - repeat last value
      }

      x += ifc_gp_trefr_x/1000.0; // Advance x over refr period
      
      // Add a rectified Gaussian value to the refractory period

      rv = ifc_gp_trefr_x_sd * nr_util_gasdev(&ifc_odeseed);

      //  Turn negative values into positive ones, thus one-sided Gaussian
      if (rv < 0.0)
	x -= rv/1000.0;
      else
	x += rv/1000.0;

      hnext = h = SIGN(h1,x2-x1);     // Start w/ default step
      // printf("   done spike\n");
    }
    /*************************************************************************/
    
    if ((x-x2)*(x2-x1) >= 0.0){ // Are we done?
      for (i=1;i<=nvar;i++)
	ystart[i] = y[i];
      if (kmax) {
	xp[++kount] = x; // Save final step
	for (i=1;i<=nvar;i++)
	  yp[i][kount] = y[i];
      }
      free_vector(dydx,1,nvar);
      free_vector(y,1,nvar);
      free_vector(yscal,1,nvar);
      return; // Normal exit
    }
    if (fabs(hnext) <= hmin){
      printf("hnext = %e\n",hnext);
      nrerror("Step size too small in odeint");
    }
    h = hnext;
    /*printf("j=%d\n",nstp);*/
  }
  printf("x = %f  nstp %d  mymaxstep %d\n",x,nstp,mymaxstep);
  nrerror("Too many steps in routine odeint");
}
#undef MAXSTP
#undef TINY
/**************************************-**************************************/
/*                                                                           */
/*                                TEST_DERIVS                                */
/*                                                                           */
/*  Fill 'dydx'.  No change to 'x' or 'y'.                                   */
/*                                                                           */
/*****************************************************************************/
void test_derivs(float x, float y[], float dydx[])
{
  int k;
  float dt,t,v,g_ex,g_in,g_lk,g_ad,v_ex,v_in,v_lk,v_ad,cm;

  v = y[1];

  // Interpolate to get g's
  t = x*ifc_gp_samp;   // Time in ms (float)
  k = (int)t;     // Time in ms (int)
  if (k < (ifc_gp_dur-1)){ // Note, k does go up to 'ifc_gp_dur'
    dt = t - (float)k;
    g_ex = ifc_gp_gtx[k] +  dt * (ifc_gp_gtx[k+1]-ifc_gp_gtx[k]);
    g_in = ifc_gp_gti[k] +  dt * (ifc_gp_gti[k+1]-ifc_gp_gti[k]);
    g_ad = ifc_gp_gta[k] +  dt * (ifc_gp_gta[k+1]-ifc_gp_gta[k]);
  }else{
    g_ex = ifc_gp_gtx[ifc_gp_dur-1];
    g_in = ifc_gp_gti[ifc_gp_dur-1];
    g_ad = ifc_gp_gta[ifc_gp_dur-1];
  }
  g_lk = ifc_gp_g_leak_x;

  v_ex = ifc_gp_v_ex;
  v_in = ifc_gp_v_in;
  v_lk = ifc_gp_v_leak_x;
  v_ad = ifc_gp_v_ad;

  cm = ifc_gp_c_x;

  dydx[1] = 1000.0/cm * (g_ex * (v_ex - v) +
			 g_in * (v_in - v) +
			 g_lk * (v_lk - v) +
			 g_ad * (v_ad - v));

  /*printf("k=%d ms    -- %f  \n",k,g_ad);*/
  /*printf("dydx=%f  x=%f   v=%f   k=%d  gex=%f  gin=%f\n",dydx[1],x,
    v,k,g_ex,g_in);*/

  /*if (k==86)
    append_farray_plot("zz.ad.pl","ad_IFC",ifc_gp_gta,ifc_gp_dur,1);*/
  
  /*printf("k=%d  gex=%f  gin=%f gad=%f  glk=%f\n",k,g_ex,g_in,g_ad,g_lk);*/

  /*printf("x,y,dydx = %f %f %f\n",x,v,dydx[1]);*/
}
/**************************************-**************************************/
/*                                                                           */
/*                               IFC_ADD_NOISE                               */
/*                                                                           */
/*  gfg - Gaussian filtered Gaussian noise                                   */
/*  gshot - Poisson shot, shot is Gaussian shaped.                           */
/*                                                                           */
/*****************************************************************************/
void ifc_add_noise(ppl,pname,data,n,seed)
     struct param_pair_list *ppl;
     char pname[];
     float *data;
     int n,seed;
{
  char *ntype;
  float *sf,*noise,*snoise,mu,sig,tsig,rate,tscale,sampling,amp;
  int tseed,ns;

  ntype = paramfile_get_nth_char_param_or_exit(ppl,pname,0);
  if (strcmp(ntype,"gfg")==0){
    mu = paramfile_get_nth_float_param_or_exit(ppl,pname,1);
    sig = paramfile_get_nth_float_param_or_exit(ppl,pname,2);
    tsig = paramfile_get_nth_float_param_or_exit(ppl,pname,3);
    noise = gaussian_corr_noise(tsig,mu,sig,n,seed);

    //append_farray_plot("zzz.noise","noise_OLD",noise,n,1);
    //printf(" ******* mu, sig, tsig = %f %f %f %d\n",mu,sig,tsig,seed);

    add_to_farray(data,noise,n);
    myfree(noise);
  }else if (strcmp(ntype,"gshot")==0){
    rate = paramfile_get_nth_float_param_or_exit(ppl,pname,1); // Hz
    sig = paramfile_get_nth_float_param_or_exit(ppl,pname,2);  // time units
    amp = paramfile_get_nth_float_param_or_exit(ppl,pname,3);  // time units
    tseed = seed;
    tscale = paramfile_get_float_param_or_exit(ppl,"tscale");
    sampling = 1.0/tscale;
    make_poisson_float_spikes(&sf,&ns,0,n,sampling,rate,0.0,0.0,1,&tseed);
    noise = expand_spike_farray(sf,ns,0,n);
    snoise = smooth_with_gaussian(noise,n,sig,0.01);
    multiply_farray(snoise,n,amp * sqrt(2.0*M_PI)*sig);
    myfree(noise); myfree(sf);
    noise = snoise;
    //append_farray_plot("zzz.noise","noise",noise,n,1);
    add_to_farray(data,noise,n);
    myfree(noise);
  }else{
    printf("ntype = %s\n",ntype);
    exit_error("IFC_ADD_NOISE","Unknown noise type");
  }
  myfree(ntype);
}
/**************************************-**************************************/
/*                                                                           */
/*                               IFC_ADD_NOISE_O                             */
/*                                                                           */
/*  gfg - Gaussian filtered Gaussian noise                                   */
/*  gshot - Poisson shot, shot is Gaussian shaped.                           */
/*                                                                           */
/*****************************************************************************/
void ifc_add_noise_o(no,data,n,seed)
     struct onode *no;
     float *data;
     int n,seed;
{
  char *ntype;
  int tseed,ns;
  float *sf,*noise,*snoise,mu,sig,tsig,rate,tscale,sampling,amp;

  ntype = onode_getpar_chr_ptr_exit(no,"type");
  if (strcmp(ntype,"gfg")==0){
    mu    = onode_getpar_flt_exit(no,"mean");
    sig   = onode_getpar_flt_exit(no,"sd");
    tsig  = onode_getpar_flt_exit(no,"tsd");
    noise = gaussian_corr_noise(tsig,mu,sig,n,seed);
    add_to_farray(data,noise,n);

    //append_farray_plot("zzz.noise","noise_NEW",noise,n,1);
    //printf(" ******* mu, sig, tsig = %f %f %f %d\n",mu,sig,tsig,seed);

    myfree(noise);
  }else if (strcmp(ntype,"gshot")==0){
    rate = onode_getpar_flt_exit(no,"rate"); // Hz
    sig  = onode_getpar_flt_exit(no,"sig");  // time units
    amp  = onode_getpar_flt_exit(no,"amp");  // time units
    tseed = seed;
    exit_error("IFC_ADD_NOISE_O","WYETH - need to get tscale passed in");
    //tscale = paramfile_get_float_param_or_exit(ppl,"tscale");
    sampling = 1.0/tscale;
    make_poisson_float_spikes(&sf,&ns,0,n,sampling,rate,0.0,0.0,1,&tseed);
    noise = expand_spike_farray(sf,ns,0,n);
    snoise = smooth_with_gaussian(noise,n,sig,0.01);
    multiply_farray(snoise,n,amp * sqrt(2.0*M_PI)*sig);
    myfree(noise); myfree(sf);
    noise = snoise;
    //append_farray_plot("zzz.noise","noise",noise,n,1);
    add_to_farray(data,noise,n);
    myfree(noise);
  }else{
    printf("ntype = %s\n",ntype);
    exit_error("IFC_ADD_NOISE_O","Unknown noise type");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                  IFC_TEST                                 */
/*                                                                           */
/*      dv                                                                   */
/*   Cm -- =  g_ex (Vex-V) + g_in (Vin-V) + g_l (Vrest-V) + g_ad (Vad-V)     */
/*      dt                                                                   */
/*                                                                           */
/*                                                                           */
/*  *rs - spike times (ms)                                                   */
/*  *rns - number of spikes                                                  */
/*                                                                           */
/*                                                                           */
/*  *** WARNING:  'gtx' and 'gti' are changed to return what was used.       */
/*                                                                           */
/*****************************************************************************/
void ifc_test(mmpid,m,pname,gtx,gti,samp,dur,seed,print_flag,dump_name,
	      rs,rns,vflag,rv,rvn)
     int mmpid; // -1-Single, 0-Master, >0-Slave
     struct model_struct *m;
     char pname[];
     float *gtx,*gti,samp;
     int dur,seed,print_flag;
     char *dump_name;      // trial name, to label plots if dumped, or NULL
     float **rs;
     int *rns;
     int vflag;   // 1-Return intracellular voltage
     float **rv;  // Return the intracellular voltage [*rvn]
     int *rvn;
{
  int i;
  int nvar,nok,nbad,dump,tseed,iseed,durms,clip_spike;
  float *ystart,x1,x2,eps,h1,hmin,*vm,tau_r_ad,tau_f_ad,*sx;
  struct ifc_param *ifp;
  char tstr[SLEN],*tlogf,*dname;
  struct onode *po,*spko,*no;

  tlogf = NULL; 
  // Set local logfile name, NULL if mmpid is -1
  (void)mylog_set_id_log(mmpid,&tlogf);
  
  ifc_gp_samp = samp;  // Samples per second
  ifc_gp_dur = dur;    // Duration in sampling units
  ifc_gp_gtx = copy_farray(gtx,dur);
  ifc_gp_gti = copy_farray(gti,dur);
  ifc_gp_gta = get_zero_farray(dur);

  if (m->ppl != NULL){
    ifp = ifc_util_get_params(m->ppl,pname);  // OLD WAY, e.g. pname='ifc'
  }else{
    po = onode_get_node_type_item_val(m->o,"pop","name",pname);
    if (po == NULL){
      printf("Population name:  %s\n",pname);
      exit_error("IFC_TEST","Cannot find population");
    }
    spko = onode_child_get_unique(po,"spike_gen");
    // WYETH - assuming that this is IFC
    ifp = ifc_util_get_param_o(spko);
  }

  ifc_gp_v_spike   = ifp->v_spike;
  ifc_gp_tau_r_ad  = ifp->tau_r_ad;
  ifc_gp_tau_f_ad  = ifp->tau_f_ad;
  ifc_gp_gbar_ad   = ifp->gbar_ad;
  ifc_gp_v_reset_x = ifp->v_reset_x;
  ifc_gp_v_ex      = ifp->v_ex;
  ifc_gp_v_in      = ifp->v_in;
  ifc_gp_v_ad      = ifp->v_ad;
  ifc_gp_v_th_x    = ifp->v_th_x;
  ifc_gp_v_leak_x  = ifp->v_leak_x;
  ifc_gp_g_leak_x  = ifp->g_leak_x;
  ifc_gp_c_x       = ifp->c_x;
  ifc_gp_trefr_x   = ifp->trefr_x;
  ifc_gp_trefr_x_sd= ifp->trefr_x_sd;
  ifc_gp_gx_scale  = ifp->gx_scale;
  ifc_gp_gx_bias   = ifp->gx_bias;
  ifc_gp_gi_scale  = ifp->gi_scale;
  ifc_gp_gi_bias   = ifp->gi_bias;
  ifc_gp_grect     = ifp->grect;
  dump             = ifp->dump;
  clip_spike       = ifp->clip_spike;

  // Scale the excit. and inhib. input, then add the bias
  multiply_farray(ifc_gp_gtx,dur,ifc_gp_gx_scale);
  add_const_farray(ifc_gp_gtx,dur,ifc_gp_gx_bias);
  multiply_farray(ifc_gp_gti,dur,ifc_gp_gi_scale);
  add_const_farray(ifc_gp_gti,dur,ifc_gp_gi_bias);

  /*
    {
    float mx,sx,mi,si;
    
    mean_sdev_farray(ifc_gp_gtx,dur,&mx,&sx);
    mean_sdev_farray(ifc_gp_gti,dur,&mi,&si);
    printf("    gtx mu %f  SD %f  (before noise)\n",mx,sx);
    printf("    gti mu %f  SD %f  (before noise)\n",mi,si);
    }
  */

  if (ifp->gx_noise != NULL){
    if (m->ppl != NULL){
      sprintf(tstr,"%s_gx_noise",pname);
      ifc_add_noise(m->ppl,tstr,ifc_gp_gtx,dur,seed);
    }else{
      no = onode_child_get_unique(spko,"gx_noise");
      if (no == NULL)
	printf("--- no = NULL -------\n");
      if (spko == NULL)
	printf("--- spko = NULL -------\n");
      ifc_add_noise_o(no,ifc_gp_gtx,dur,seed);
    }
  }
  
  // Make a seed for odeint
  tseed = seed;
  if (tseed > 0)
    tseed = -tseed;
  ifc_odeseed = 69513 * myrand_util_ran2(&tseed);
  if (ifc_odeseed > 0)
    ifc_odeseed = -ifc_odeseed;


  if (ifp->gi_noise != NULL){
    iseed = 170013 * myrand_util_ran2(&tseed); // Make 2nd seed from 1st
    if (m->ppl != NULL){
      sprintf(tstr,"%s_gi_noise",pname);
      ifc_add_noise(m->ppl,tstr,ifc_gp_gti,dur,iseed);
    }else{
      no = onode_child_get_unique(spko,"gi_noise");
      // WYETH fixed two bugs on Sept 16, 2008  _gtx -> gti, seed -> iseed
      ifc_add_noise_o(no,ifc_gp_gti,dur,iseed);
    }
  }
  if (ifc_gp_grect){
    half_wave_rectify_farray(ifc_gp_gtx,dur);
    half_wave_rectify_farray(ifc_gp_gti,dur);
  }

  // Get adaptation shapes
  tau_r_ad = ifc_gp_tau_r_ad * samp/1000.0;
  tau_f_ad = ifc_gp_tau_f_ad * samp/1000.0;
  if (ifc_gp_gbar_ad > 0.0){
    ifc_gp_gta_npulse = (int)(ifc_gp_tau_f_ad*4.0 * samp/1000.0);
    ifc_gp_gta_pulse = diff_exp_farray(0.0,ifc_gp_gbar_ad,tau_f_ad,tau_r_ad,
				       ifc_gp_gta_npulse);
  }else{
    ifc_gp_gta_npulse = 0;
    ifc_gp_gta_pulse = (float *)NULL;
  }

  if (mmpid == -1){

    //printf("HERE dump = %d\n",dump);
    //printf("dump_name = %s\n",dump_name);

    if (dump_name == NULL){

      if (dump){
	if (ifc_gp_gbar_ad > 0.0){
	  append_farray_plot(DUMP_FNAME,"gtaIMP",ifc_gp_gta_pulse,
			     ifc_gp_gta_npulse,1);
	}
	append_farray_plot(DUMP_FNAME,"gtx",gtx,dur,1);
	append_farray_plot(DUMP_FNAME,"gti",gti,dur,1);
	append_farray_plot(DUMP_FNAME,"gtx_noise",ifc_gp_gtx,dur,1);
	append_farray_plot(DUMP_FNAME,"gti_noise",ifc_gp_gti,dur,1);
      }
    }else{
      if (ifc_gp_gbar_ad > 0.0){
	sprintf(tstr,"%s_gtaIMP",dump_name);
	append_farray_plot(DUMP_FNAME,tstr,ifc_gp_gta_pulse,
			   ifc_gp_gta_npulse,1);
      }
      sprintf(tstr,"%s_gtx",dump_name);
      append_farray_plot(DUMP_FNAME,tstr,gtx,dur,1);
      sprintf(tstr,"%s_gti",dump_name);
      append_farray_plot(DUMP_FNAME,tstr,gti,dur,1);
      sprintf(tstr,"%s_gtx_noise",dump_name);
      append_farray_plot(DUMP_FNAME,tstr,ifc_gp_gtx,dur,1);
      sprintf(tstr,"%s_gti_noise",dump_name);
      append_farray_plot(DUMP_FNAME,tstr,ifc_gp_gti,dur,1);
    }
  }
  
  /************************* Setup for 'odeint' ******************************/
  
  /*** Variables, initial values. ***/
  nvar = 1;
  ystart = vector(1,nvar);
  ystart[1] = -70.0;
  x1 = 0.0;
  x2 = (float)dur/samp; /* seconds */
  dxsav = 0.001; /* Save values no more often than this (1 msec) */
  kmax = (x2-x1)/dxsav;
  /*printf("    Duration %f  dxsav %f, thus, kmax = %d\n",x2-x1,dxsav,kmax);*/

  /*** Intermediate storage. ***/
  /*     kount will hold current number of entries in vector. */  
  xp = vector(1,kmax); 
  yp = matrix(1,nvar,1,kmax);

  gspkn = 0; /* Set spike counter to zero */

  h1 = 0.001;   /* 0.001 = 1ms */
  eps = 0.0001; /* Accuracy */
  hmin = 1.0e-16; /* Minimal allowed step size. */

  /*printf("HERE WYETH  x1 x2 %f %f  h1,min %f %f  eps %f\n",
    x1,x2,h1,hmin,eps);*/
  /*printf("ystart = %f  nvar = %d\n",ystart[1],nvar);*/

  ifc_odeint(ystart,nvar,x1,x2,eps,h1,hmin,&nok,&nbad,test_derivs,ifc_rkqs);

  /*
    printf("    Done 'odeint', nok = %d  nbad = %d\n",nok,nbad);
    printf("    Number of values stored = %d (%d limit)\n",kount,kmax);
    */

  if (print_flag){
    sprintf(tstr,"    IFC_TEST  %4d spikes  (ok/bad %d/%d; stored %d of %d)\n",
	    gspkn,nok,nbad,kount,kmax);
    mylog(tlogf,tstr);
  }

  /*** Compute membrane voltage. ***/
  vm = (float *)myalloc(kount*sizeof(float));
  for(i=0;i<kount;i++){
    vm[i] = yp[1][1+i];
    xp[i+1] *= 1000.0; /* Make units be in ms */
  }

  if (mmpid == -1){
    if (dump_name == (char *)NULL){
      if (dump){
	append_farray_xy_plot(DUMP_FNAME,xp+1,vm,kount,"Vm");
	append_farray_plot(DUMP_FNAME,"gta",ifc_gp_gta,dur,1);
      }
    }else{
      sprintf(tstr,"%s_Vm",dump_name);
      append_farray_xy_plot(DUMP_FNAME,xp+1,vm,kount,tstr);
      sprintf(tstr,"%s_gta",dump_name);
      append_farray_plot(DUMP_FNAME,tstr,ifc_gp_gta,dur,1);
    }
  }

  if (vflag){ // Return Vm
    if (clip_spike)
      for(i=0;i<kount;i++){  // Set high V associated w/ spike to Vthresh
	printf("i=%d (of %d)\n",i,kount);
	if (vm[i] == ifc_gp_v_spike)
	  vm[i] = ifc_gp_v_th_x;
      }

    durms = my_rint(x2 * 1000.0); // Duration in msec
    sx = get_zero_farray(durms);
    for(i=0;i<durms;i++)
      sx[i] = (float)i;

    *rv = resample_farray(xp+1,vm,kount,sx,durms); // returned in msec
    *rvn = durms;
    /*append_farray_xy_plot("zzz.pl",sx,vf,durms,"NAME_ms");*/
    /*exit(0);*/
    myfree(sx);
  }else{
    *rv = (float *)NULL;
    *rvn = 0;
  }

  for(i=0;i<dur;i++){ // Return the conductances that were used
    gtx[i] = ifc_gp_gtx[i];
    gti[i] = ifc_gp_gti[i];
  }

  myfree(ifc_gp_gtx); myfree(ifc_gp_gti); myfree(ifc_gp_gta); myfree(vm);
  free_vector(ystart,1,nvar); free_vector(xp,1,kmax);
  free_matrix(yp,1,nvar,1,kmax);

  ifc_util_free_ifc_param(ifp);

  *rs = copy_farray(gspk,gspkn);
  *rns = gspkn;
}
/**************************************-**************************************/
/*                                                                           */
/*                              GET_GT_FROM_SARRAY                           */
/*                                                                           */
/*  Sum spike trains and convolve with EPSC shape.                           */
/*  - Returned arrays are in time units of milliseconds.                     */
/*                                                                           */
/*****************************************************************************/
float *get_gt_from_sarray(s,cnt,n,tn,t1,t2,amp)
     float **s;         // Spike times in ms
     int *cnt,n;
     int tn;            // Duration in sampling units (should be ms)
     float t1,t2,amp;   // Time constants (match tn units), ampl for EPSC
{
  int i,j;
  int nx,t;
  float *xshape,*gt,*g;

  nx = (int)(t1 * 7.0);        // Length of mask is 7 times decay time
  xshape = diff_exp_farray(0.0,amp,t1,t2,nx);

  gt = get_zero_farray(tn);
  for(i=0;i<n;i++){
    for(j=0;j<cnt[i];j++){
      t = (int)(0.5 + s[i][j]);
      gt[t] += 1.0;
    }
  }
  //append_farray_plot("zz.ifc0202.pl","summed_spikes",gt,tn,1);
  append_farray_plot("zz.ifc0202.pl","xshape",xshape,nx,1);

  g = convolve_with_mask_causal(gt,tn,xshape,nx);
  myfree(xshape);
  myfree(gt);

  /*append_farray_plot("zz.ifc0202.pl","gt",g,tn,1);*/
  
  return g;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_IFC_01_PREP                              */
/*                                                                           */
/*****************************************************************************/
void mod_ifc_01_prep(m,s)
     struct model_struct *m;   // Model params
     struct stim_struct *s;    // Stimulus data and params
{
  int seed;
  float tscale,samp;
  struct onode *mo;

  mylog(mylogf,"  MOD_IFC_01_PREP\n");

  printf("    Preparing to run %d trials\n",s->ntr);

  if (m->o == NULL)
    mylogx(mylogf,"MOD_IFC_01_PREP","Convert to .moo file");


  ifc_01_tn = onode_getpar_int_exit(m->o,"tn");

  //  Sampling
  tscale = onode_getpar_flt_exit(m->o,"tscale");
  ifc_01_samp = 1.0/tscale;

  ifc_01_tdelay = onode_getpar_flt_dflt(m->o,"tdelay",0.0);
  ifc_01_sigsm  = onode_getpar_flt_dflt(m->o,"smooth_tsd",0.0);

  ifc_01_stim_use   = onode_getpar_chr_dflt(m->o,"stim_usage","g_ex");
  if (!((strcmp(ifc_01_stim_use,"rate_ex")==0)||
	(strcmp(ifc_01_stim_use,"rate_in")==0)||
	(strcmp(ifc_01_stim_use,"g_ex")==0)||
	(strcmp(ifc_01_stim_use,"g_in")==0)))
    mylogx(mylogf,"MOD_IFC_01_PREP","Unknown 'stim_usage' value.");


  //  Check for BG inputs
  ifc_01_bg_ex_rate = onode_getpar_flt_dflt(m->o,"bg_ex_rate",-1.0);
  ifc_01_bg_in_rate = onode_getpar_flt_dflt(m->o,"bg_in_rate",-1.0);
  ifc_01_bg_ex_amp  = onode_getpar_flt_dflt(m->o,"bg_ex_amp",1.0);
  ifc_01_bg_in_amp  = onode_getpar_flt_dflt(m->o,"bg_in_amp",1.0);

  //  Check for pre-synaptic spiking input
  ifc_01_pre_ex_rate = onode_getpar_flt_dflt(m->o,"presyn_ex_rate",-1.0);
  ifc_01_pre_ex_amp  = onode_getpar_flt_dflt(m->o,"presyn_ex_amp",1.0);

  ifc_01_pre_in_rate = onode_getpar_flt_dflt(m->o,"presyn_in_rate",-1.0);
  ifc_01_pre_in_amp  = onode_getpar_flt_dflt(m->o,"presyn_in_amp",1.0);

  //
  //  If the stimulus is to be used to drive the firing rate of the presyn
  //  'ex' or 'in' inputs, then there should be no other 'rate' value for them.
  //
  if ((strcmp(ifc_01_stim_use,"rate_ex")==0) && (ifc_01_pre_ex_rate >= 0.0))
    mylogx(mylogf,"MOD_IFC_01_PREP","Stim usage conflicts with presyn_ex_rate");
  if ((strcmp(ifc_01_stim_use,"rate_in")==0) && (ifc_01_pre_in_rate >= 0.0))
    mylogx(mylogf,"MOD_IFC_01_PREP","Stim usage conflicts with presyn_in_rate");


  if ((ifc_01_bg_ex_rate > 0.0) || (ifc_01_pre_ex_rate > 0.0) ||
      (strcmp(ifc_01_stim_use,"rate_ex")==0)){
    //
    //  Get masks for PSG shapes
    //
    mo = onode_get_node_type_item_val(m->o,"mech","name","ex");
    ifc_01_ex_mask = popu_mech_onode_get_mask_alpha(mo,ifc_01_samp,
						    &ifc_01_ex_mask_n);

    //append_farray_plot("zz.mask.pl","ex",ifc_01_ex_mask,ifc_01_ex_mask_n,1);
  }else{
    ifc_01_ex_mask = NULL;
  }

  if ((ifc_01_bg_in_rate > 0.0) || (ifc_01_pre_in_rate > 0.0) ||
      (strcmp(ifc_01_stim_use,"rate_in")==0)){
    
    mo = onode_get_node_type_item_val(m->o,"mech","name","in");
    ifc_01_in_mask = popu_mech_onode_get_mask_alpha(mo,ifc_01_samp,
						    &ifc_01_in_mask_n);

    //append_farray_plot("zz.mask.pl","in",ifc_01_in_mask,ifc_01_in_mask_n,1);
  }else{
    ifc_01_in_mask = NULL;
  }

  //
  //  Create random seeds for any pre-syn and background spike input
  //
  if ((ifc_01_pre_ex_rate > 0.0)||(strcmp(ifc_01_stim_use,"rate_ex")==0)){
    seed = onode_getpar_flt_exit(m->o,"presyn_ex_seed");
    ifc_01_pre_ex_seed = get_seeds(seed,100000,s->ntr);
  }
  if ((ifc_01_pre_in_rate > 0.0)||(strcmp(ifc_01_stim_use,"rate_in")==0)){
    seed = onode_getpar_flt_exit(m->o,"presyn_in_seed");
    ifc_01_pre_in_seed = get_seeds(seed,100000,s->ntr);
  }
  if (ifc_01_bg_ex_rate > 0.0){
    seed = onode_getpar_flt_exit(m->o,"bg_ex_seed");
    ifc_01_bg_ex_seed = get_seeds(seed,100000,s->ntr);
  }
  if (ifc_01_bg_in_rate > 0.0){
    seed = onode_getpar_flt_exit(m->o,"bg_in_seed");
    ifc_01_bg_in_seed = get_seeds(seed,100000,s->ntr);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          MODEL_IFC_01_GET_RESPONSE                        */
/*                                                                           */
/*****************************************************************************/
void model_ifc_01_get_response(m,s,r,k)
     struct model_struct *m;        // Model parameter pair list
     struct stim_struct *s;         // Stimulus data and param pair list
     struct response_struct *r;     // Response data
     int k;                         // Repeat number
{
  int i;
  int tn,ns,fvn,vflag,dumpflag,s_pre_xn,s_pre_in,pflag;
  float samp,*gtx,*gti,*fs,*fv,*tstim;
  float *fs_pre_x,*fs_pre_i,*gx_pre,*gi_pre,*gx_bg,*gi_bg;
  char *dname;

  pflag = 0;

  tn = ifc_01_tn;
  samp = ifc_01_samp;

  //
  //  Smoothing the stimulus, or just copy it if no smoothing.
  //
  if (ifc_01_sigsm <= 0.0){
    tstim = copy_farray(s->d1,tn);
  }else{
    mylog(mylogf,"  MODEL_IFC_01_GET_RESPONSE  Smooth SD %f\n",ifc_01_sigsm);
    tstim = smooth_with_gaussian(s->d1,tn,ifc_01_sigsm,0.01);
  }
  if (pflag){
    printf("    Writing the stimulus trace to 'zz.tstim.pl'\n");
    append_farray_plot("zz.tstim.pl","tstim",tstim,tn,1);
  }


  //
  //  Pre-synaptic spikes an g's, and BG g's
  //
  gx_pre = gx_bg = NULL;
  gi_pre = gi_bg = NULL;
  fs_pre_x = fs_pre_i = NULL;

  //
  //  Presyn EX
  //
  if (ifc_01_pre_ex_rate > 0.0){ // Generate Poisson spikes w/ constant rate
    gx_pre = spikeu_poisson_g(tn,samp,
			      ifc_01_pre_ex_rate,
			      ifc_01_pre_ex_seed[k],
			      ifc_01_ex_mask,
			      ifc_01_ex_mask_n,&fs_pre_x,&s_pre_xn);
    multiply_farray(gx_pre,tn,ifc_01_pre_ex_amp);
    if (pflag) printf("    EX spikes made from rate:  %d spikes\n",s_pre_xn);
  }else if ((strcmp(ifc_01_stim_use,"rate_ex")==0)){  // spikes from stim prob.
    multiply_farray(tstim,tn,1.0/samp);  // Convert to "prob"
    gx_pre = spikeu_poisson_g_prob(tstim,tn,samp,
				   ifc_01_pre_ex_seed[k],
				   ifc_01_ex_mask,
				   ifc_01_ex_mask_n,&fs_pre_x,&s_pre_xn);
    if (pflag) printf("    EX spikes made from stim:  %d spikes\n",s_pre_xn);
    //append_farray_plot("zz.gx_pre.pl","gx_pre",gx_pre,tn,1);
  }else{
    gx_pre = get_zero_farray(tn);
    s_pre_xn = 0;
    if (pflag) printf("    EX spikes not made\n");
  }

  //
  //  Presyn IN
  //
  if (ifc_01_pre_in_rate > 0.0){
    gi_pre = spikeu_poisson_g(tn,samp,
			      ifc_01_pre_in_rate,
			      ifc_01_pre_in_seed[k],
			      ifc_01_in_mask,
			      ifc_01_in_mask_n,&fs_pre_i,&s_pre_in);
    multiply_farray(gi_pre,tn,ifc_01_pre_in_amp);
    if (pflag) printf("    IN spikes made from rate:  %d spikes\n",s_pre_in);
  }else if ((strcmp(ifc_01_stim_use,"rate_in")==0)){
    multiply_farray(tstim,tn,1.0/samp);  // Convert to "prob"
    gi_pre = spikeu_poisson_g_prob(tstim,tn,samp,
				   ifc_01_pre_in_seed[k],
				   ifc_01_in_mask,
				   ifc_01_in_mask_n,&fs_pre_i,&s_pre_in);
    if (pflag) printf("    IN spikes made from stim:  %d spikes\n",s_pre_in);
  }else{
    gi_pre = get_zero_farray(tn);
    s_pre_in = 0;
    if (pflag) printf("    IN spikes not made\n");
  }

  //
  //  BG
  //
  if (ifc_01_bg_ex_rate > 0.0){
    gx_bg = spikeu_poisson_g(tn,samp,
			     ifc_01_bg_ex_rate,
			     ifc_01_bg_ex_seed[k],
			     ifc_01_ex_mask,
			     ifc_01_ex_mask_n,NULL,NULL);
    multiply_farray(gx_bg,tn,ifc_01_bg_ex_amp);
    if (pflag) printf("    BG spikes made from EX rate.\n");
  }else{
    gx_bg = get_zero_farray(tn);
    if (pflag) printf("    BG spikes not made for EX.\n");
  }

  if (ifc_01_bg_in_rate > 0.0){
    gi_bg = spikeu_poisson_g(tn,samp,
			     ifc_01_bg_in_rate,
			     ifc_01_bg_in_seed[k],
			     ifc_01_in_mask,
			     ifc_01_in_mask_n,NULL,NULL);
    multiply_farray(gi_bg,tn,ifc_01_bg_in_amp);
    if (pflag) printf("    BG spikes made from IN rate.\n");
  }else{
    gi_bg = get_zero_farray(tn);
    if (pflag) printf("    BG spikes not made for IN.\n");
  }


  //
  //  Create EX conductance, 'gtx'
  //
  if (strcmp(ifc_01_stim_use,"g_ex")==0)
    gtx = copy_farray(tstim,tn);
  else
    gtx = get_zero_farray(tn);

  add_to_farray(gtx,gx_pre,tn);
  add_to_farray(gtx,gx_bg,tn);


  //
  //  Create IN conductance, 'gti'
  //
  if (strcmp(ifc_01_stim_use,"g_in")==0)
    gti = copy_farray(tstim,tn);
  else
    gti = get_zero_farray(tn);

  add_to_farray(gti,gi_pre,tn);
  add_to_farray(gti,gi_bg,tn);



  if (r->nf > 0) // Will save Vm
    vflag = 1;
  else
    vflag = 0;

  dumpflag = 0;
  //if (paramfile_test_param(m->ppl,"ifc_dump")){
  if (onode_get_item_child(m->o,"ifc_dump")!=NULL){
    exit_error("MODEL_IFC_01_GET_RESPONSE","change ifc_dump --> spike_dump");
    //}else if (paramfile_test_param(m->ppl,"spike_dump")){
  }else if (onode_get_item_child(m->o,"spike_dump")){
    dumpflag = paramfile_get_int_param_or_exit(m->ppl,"spike_dump");
  }
  if (dumpflag == 0)
    dname = (char *)NULL;
  else
    dname = r->tname;      // Current trial name, or NULL if no var params

  //printf("dumpflag = %d\n",dumpflag);

  ifc_test(-1,m,"ifc",gtx,gti,samp,tn,m->mseed[r->tsi],1,dname,
	   &fs,&ns,vflag,&fv,&fvn);
  myfree(gti);

  /* Spikes come in ms, adjust based on the 'save_spikes' sampling */
  /* WYETH - need to think about adjusting sampling ... 
     ssamp = paramfile_get_nth_float_param_or_exit(m->ppl,"save_spikes",1);
     multiply_farray(fs,ns,ssamp/1000.0);***/
  
  // Apply time delay to spikes
  if (ifc_01_tdelay != 0.0){
    add_const_farray(fs,ns,ifc_01_tdelay*1000.0);
  }


  //
  //  Store responses
  //
  for(i=0;i<r->n;i++){ // For each response requested
    if (strcmp(r->datid[i],"spikes")==0){
      mod_util_resp_store_s(r,i,fs,ns,0,mylogf); // store ptr, copyflag = 0
    }else if (strcmp(r->datid[i],"spikes_pre_ex")==0){
      mod_util_resp_store_s(r,i,fs_pre_x,s_pre_xn,0,mylogf); // copyflag = 0
    }else if (strcmp(r->datid[i],"spikes_pre_in")==0){
      mod_util_resp_store_s(r,i,fs_pre_i,s_pre_in,0,mylogf); // copyflag = 0
    }else if (strcmp(r->datid[i],"vm")==0){
      //mod_util_resp_store_f(r,i,fv,fvn,0,mylogf);
      mod_util_resp_store_f_samp(r,i,fv,fvn,samp,mylogf); // Data copied
    }else if (strcmp(r->datid[i],"bg_gx")==0){
      mod_util_resp_store_f_samp(r,i,gx_bg,tn,samp,mylogf); // Data copied
    }else if (strcmp(r->datid[i],"bg_gi")==0){
      mod_util_resp_store_f_samp(r,i,gi_bg,tn,samp,mylogf); // Data copied
    }else if (strcmp(r->datid[i],"pre_gx")==0){
      mod_util_resp_store_f_samp(r,i,gx_pre,tn,samp,mylogf); // Data copied
    }else if (strcmp(r->datid[i],"pre_gi")==0){
      mod_util_resp_store_f_samp(r,i,gi_pre,tn,samp,mylogf); // Data copied
    }else if (strcmp(r->datid[i],"tot_gx")==0){
      mod_util_resp_store_f_samp(r,i,gtx,tn,samp,mylogf); // Data copied
    }else if (strcmp(r->datid[i],"tot_gi")==0){
      mod_util_resp_store_f_samp(r,i,gti,tn,samp,mylogf); // Data copied
    }
  }

  if (gx_pre != NULL)  myfree(gx_pre);
  if (gi_pre != NULL)  myfree(gi_pre);
  if (gx_bg  != NULL)  myfree(gx_bg);
  if (gi_bg  != NULL)  myfree(gi_bg);
  // Do not free spikes, in case pointers are used above for storage
  // *** NOTE:  some spikes may be neither freed nor stored

}
/**************************************-**************************************/
/*                                                                           */
/*                          MODEL_IFC_02_GET_RESPONSE                        */
/*                                                                           */
/*****************************************************************************/
void model_ifc_02_get_response(m,s,r,k)
     struct model_struct *m;        // Model parameter pair list
     struct stim_struct *s;         // Stimulus data and param pair list
     struct response_struct *r;     // Response data
     int k;                         // Repeat number
{
  int tn,ns,fvn;
  float tscale,samp,ssamp,*gtx,*gti,*fs,*fv;
  
  tn = paramfile_get_int_param_or_exit(m->ppl,"tn");
  tscale = paramfile_get_float_param_or_exit(m->ppl,"tscale");
  samp = 1.0/tscale;

  // Determine both inhibitory and excitatory conductances based on the
  // 1D time sequence in 's->d1'
  gtx = copy_farray(s->d1,tn);
  half_wave_rectify_farray(gtx,tn);
  gti = copy_farray(s->d1,tn);
  multiply_farray(gti,tn,-1.0);
  half_wave_rectify_farray(gti,tn);
  // WYETH - could add relative shift or filtering here

  ifc_test(-1,m,"ifc",gtx,gti,samp,tn,m->mseed[r->tsi],1,(char *)NULL,
	   &fs,&ns,0,&fv,&fvn);
  // WEYTH - voltage trace, fv, is thrown away currently

  // Spikes come in ms, adjust based on the 'save_spikes' sampling
  ssamp = paramfile_get_nth_float_param_or_exit(m->ppl,"save_spikes",1);
  multiply_farray(fs,ns,ssamp/1000.0);
  
  // Store response, update counter
  if (r->ns > 0){
    //r->s[0][r->tsi] = fs;
    //r->cnt[0][r->tsi] = ns; // 'ti' is updated in 'run_gen_tuning_curve'
    r->s[0][r->gtsi] = fs;
    r->cnt[0][r->gtsi] = ns; // 'ti' is updated in 'run_gen_tuning_curve'
  }
  myfree(gtx); myfree(gti);
}
/**************************************-**************************************/
/*                                                                           */
/*                        MODEL_IFC_02_02_GET_RESPONSE                       */
/*                                                                           */
/*  Two stages of IF units.                                                  */
/*                                                                           */
/*****************************************************************************/
void model_ifc_02_02_get_response(m,s,r,k)
     struct model_struct *m;        // Model parameter pair list.
     struct stim_struct *s;         // Stimulus data and param pair list
     struct response_struct *r;     // Response data
     int k;                         // Repeat number
{
  int i;
  int tn,ns,fvn,n01,vflag,*layseed,tseed,*cnt,nxsp;
  float tscale,samp,ssamp,tdelay,*gtx,*gti,*tgtx,*tgti,*fs,**fs1,*fv;
  float t1x,t2x,ampx;

  tn = paramfile_get_int_param_or_exit(m->ppl,"tn");
  tscale = paramfile_get_float_param_or_exit(m->ppl,"tscale");
  tdelay = paramfile_get_float_param_default(m->ppl,"tdelay_sec",0.0);
  samp = 1.0/tscale;

  // Get number of units in first layer
  n01 = paramfile_get_int_param_or_exit(m->ppl,"nlayer1");
  printf("    %d units in first layer\n",n01);

  // Make a seed for each unit in layer 1
  layseed = (int *)myalloc(n01*sizeof(int));
  tseed = m->mseed[r->tsi];
  if (tseed > 0)
    tseed = -tseed;
  for(i=0;i<n01;i++){
    layseed[i] = 1723981*myrand_util_ran2(&tseed);
  }

  /*** Determine both inhibitory and excitatory conductances based on the
       1D time sequence in 's->d1' ***/
  gtx = copy_farray(s->d1,tn);
  half_wave_rectify_farray(gtx,tn);
  gti = copy_farray(s->d1,tn);
  multiply_farray(gti,tn,-1.0);
  half_wave_rectify_farray(gti,tn);
  /*** WYETH - could add relative shift or filtering here ***/

  vflag = 0;  /* 0 --> Do not return Vm in 'fv' */

  /*printf("tn = %d,  sampl = %f\n",tn,samp);*/
  fs1 = (float **)myalloc(n01*sizeof(float *));
  cnt = (int *)myalloc(n01*sizeof(int));
  for(i=0;i<n01;i++){
    printf("  Layer1 unit %d",i); fflush(stdout);
    tgtx = copy_farray(gtx,tn); /* Use copy, values changed in ifc_test */
    tgti = copy_farray(gti,tn);
    ifc_test(-1,m,"ifc",tgtx,tgti,samp,tn,layseed[i],1,(char *)NULL,
	     &fs,&ns,vflag,&fv,&fvn);
    myfree(tgtx); myfree(tgti);
    fs1[i] = fs;
    cnt[i] = ns;
  }
  myfree(gtx); myfree(gti);
  // WYETH - voltage trace, fv, is thrown away currently

  // Get time constants in msec for EPSCs
  t1x   = paramfile_get_float_param_or_exit(m->ppl,"lay2_ex_tau_f");
  t2x   = paramfile_get_float_param_or_exit(m->ppl,"lay2_ex_tau_r");
  ampx  = paramfile_get_float_param_or_exit(m->ppl,"lay2_ex_amp");

  nxsp = (int)(1000.0 * (float)tn/samp + 0.5); // Stim dur in msec
  gtx = get_gt_from_sarray(fs1,cnt,n01,nxsp,t1x,t2x,ampx);
  gti = get_zero_farray(nxsp);

  vflag = 0;  // 0 --> Do not return Vm in 'fv'
  samp = 1000.0;
  ifc_test(-1,m,"ifc_02",gtx,gti,samp,nxsp,m->mseed[r->tsi],1,(char *)NULL,
	   &fs,&ns,vflag,&fv,&fvn);
  /*** WYETH - voltage trace, fv, is thrown away currently ***/

  /* Spikes come in ms, adjust based on the 'save_spikes' sampling */
  /*
    ssamp = paramfile_get_nth_float_param_or_exit(m->ppl,"save_spikes",1);
    multiply_farray(fs,ns,ssamp/1000.0);*/
  /*** WYETH - need to think about adjusting sampling units ***/

  /*** Apply time delay to spikes ***/
  add_const_farray(fs,ns,tdelay*1000.0);
  
  /*** Store response, link to storage, copyflag = 0 ***/
  mod_util_resp_store_s(r,0,fs,ns,0,mylogf);

}
/**************************************-**************************************/
/*                                                                           */
/*                         MODEL_IFC_03_GET_RESPONSE                         */
/*                                                                           */
/*  - Use excitatory and inhibitory conductances.                            */
/*  - convolve w/ stimulus dependent filters.                                */
/*                                                                           */
/*****************************************************************************/
void model_ifc_03_get_response(m,s,r,k)
     struct model_struct *m;        /* Model parameter pair list. */
     struct stim_struct *s;         /* Stimulus data and param pair list */
     struct response_struct *r;     /* Response data */
     int k;                         /* Repeat number */
{
  int tn,ns,fn,aux1,fvn;
  float tscale,samp,ssamp,*gtx,*gti,*fs,*ff,*fv;
  float *gtxf,*gtif;  /* Filtered input */
  
  tn = paramfile_get_int_param_or_exit(m->ppl,"tn");
  tscale = paramfile_get_float_param_or_exit(m->ppl,"tscale");
  samp = 1.0/tscale;

  {
    float mu,sig,fm,fs;

    mu = paramfile_get_float_param_or_exit(s->ppl,"mu");
    sig = paramfile_get_float_param_or_exit(s->ppl,"sig");
    aux1 = paramfile_get_float_param_or_exit(s->ppl,"aux1");
    printf("  mu = %f  sig = %f  aux1 = %d\n",mu,sig,aux1);

    fs = (float)aux1;  /* 'aux1' controls kernel width */

    /*** For maxwell function s=10 -> 16 units 1/2 width ***/
    fn = 256;
    fm = 10.0; /* start time */
    ff = maxwell_farray(fm,fs,1.0,fn);
    append_farray_plot("zf.pl","ff",ff,fn,1);
  }

  /*** Determine both inhibitory and excitatory conductances based on the
    1D time sequence in 's->d1' ***/
  gtx = copy_farray(s->d1,tn);
  half_wave_rectify_farray(gtx,tn);
  gti = copy_farray(s->d1,tn);
  multiply_farray(gti,tn,-1.0);
  half_wave_rectify_farray(gti,tn);

  /*** Non-linear processing of inhibitory signal to suppress singles ***/
  {
    int i,dmn;
    float *dmask,mean,m1,tfrac,max;
    
    dmn = 32;
    dmask = get_zero_farray(dmn);
    for(i=0;i<10;i++)
      dmask[i] = 0.1;
    /*
    for(i=10;i<15;i++)
      dmask[i] = -0.6;
    for(i=0;i<32;i++)
      dmask[i] = 0.03125;*/

    /*append_farray_plot("zgi.pl","gi",gti,tn,1);*/
    gtif = convolve_with_mask_causal(gti,tn,dmask,dmn);
    myfree(gti);
    gti = gtif;
    mean = mean_farray(gti,tn);
    max = max_of_farray(gti,tn);

    tfrac = 1.2; /* Threshold fraction */
    add_const_farray(gti,tn,-tfrac*mean);
    half_wave_rectify_farray(gti,tn);
    /*m1 = mean_farray(gti,tn);*/
    m1 = max_of_farray(gti,tn);
    multiply_farray(gti,tn,max/m1);
    /*append_farray_plot("zgi.pl","Fgi",gti,tn,1);*/

    /*** Now add a transient component to the inhibition ***/
    for(i=0;i<10;i++)
      dmask[i] = 0.4;
    for(i=10;i<15;i++)
      dmask[i] = -0.6;
    for(i=15;i<dmn;i++)
      dmask[i] = 0.0;
    gtif = convolve_with_mask_causal(gti,tn,dmask,dmn);
    myfree(dmask);
    myfree(gti);
    gti = gtif;
    /*append_farray_plot("zgi.pl","Tgi",gti,tn,1);*/
  }


  gtxf = convolve_with_mask_causal(gtx,tn,ff,fn);
  myfree(gtx);
  gtx = gtxf;

  gtif = convolve_with_mask_causal(gti,tn,ff,fn);
  myfree(gti);
  gti = gtif;

  /*** Delay inhibition ***/
  /*
    {
    int nshift;
    nshift = 5;  shift_farray(gti,tn,nshift,0.0);
    }*/

  ifc_test(-1,m,"ifc",gtx,gti,samp,tn,m->mseed[r->tsi],1,(char *)NULL,&fs,&ns,
	   0,&fv,&fvn);

  // Spikes come in ms, adjust based on the 'save_spikes' sampling
  ssamp = paramfile_get_nth_float_param_or_exit(m->ppl,"save_spikes",1);
  multiply_farray(fs,ns,ssamp/1000.0);

  // Store response, update counter
  if (r->ns > 0){
    //r->s[0][r->tsi] = fs;
    //r->cnt[0][r->tsi] = ns; // 'ti' is updated in 'run_gen_tuning_curve'
    r->s[0][r->gtsi] = fs;
    r->cnt[0][r->gtsi] = ns; // 'ti' is updated in 'run_gen_tuning_curve'
  }
  myfree(gtx); myfree(gti);
}
/**************************************-**************************************/
/*                                                                           */
/*                              IFC_UTIL_POISSON                             */
/*                                                                           */
/*  WYETH - THIS SHOULD BE PUT SOMEWHERE ELSE,                               */
/*  'mod_util' PERHAPS, or make a NEW FILE:  mod_spike_util ...              */
/*                                                                           */
/*  Returns spikes in 'msec'.                                                */
/*                                                                           */
/*****************************************************************************/
void ifc_util_poisson(mylogf,m,sgo,d,n,sampling,seed,pflag,vflag,
		      rfs,rfn,rv,rvn)
     char *mylogf;              // Logfile name
     struct model_struct *m;    // Model parameters, for m->ppl
     struct onode *sgo;         // Onode in which 'spike_gen' lives
     float *d;                  // Response vs. time [n]
     int n;                     // Length of 'd'
     float sampling;            // Samples/sec of 'd'
     int seed;                  // Randomization seed
     int pflag;                 // Print flag
     int vflag;                 // (0) 1-Return poisson prob
     float **rfs;               // Returned spikes
     int    *rfn;               // Number of spikes
     float **rv;                // Returned float data
     int    *rvn;               // Length of float data
{
  int i;
  int nflag,ns,dumpflag,tn,aseed;
  float px,pa,pa0,*t,*fs,ptoff,refmu,refsig,autoscale,totspikes,area;
  float noise_mu,noise_sd,noise_tsd,*noise;
  struct onode *spko;

  if (m->ppl != NULL){  // .mod
    px = paramfile_get_float_param_default(m->ppl,"poisson_scale",1.0);
    pa = paramfile_get_float_param_default(m->ppl,"poisson_offset",0.0);
    pa0 = paramfile_get_float_param_default(m->ppl,"poisson_offset0",0.0);
    ptoff = paramfile_get_float_param_default(m->ppl,"poisson_toff_ms",0.0);
    dumpflag = paramfile_get_int_param_default(m->ppl,"spike_dump",0);
    if (dumpflag == 0) // WYETH - Ultimately get rid of 'poisson_dump'
      dumpflag = paramfile_get_int_param_default(m->ppl,"poisson_dump",0);

    autoscale = paramfile_get_float_param_default(m->ppl,"autoscale",-1.0);
    refmu  = 0.0;
    refsig = 0.0;

  }else{ // .moo

    //
    //  WYETH - get rid of 'poisson_' names eventually
    //

    if (onode_get_item_child(sgo,"offset0") != NULL)
      pa0      = onode_getpar_flt_dflt(sgo,"offset0",0.0);
    else
      pa0      = onode_getpar_flt_dflt(sgo,"poisson_offset0",0.0);

    if (onode_get_item_child(sgo,"scale") != NULL)
      px       = onode_getpar_flt_dflt(sgo,"scale"  ,1.0);
    else
      px       = onode_getpar_flt_dflt(sgo,"poisson_scale"  ,1.0);

    if (onode_get_item_child(sgo,"offset") != NULL)
      pa       = onode_getpar_flt_dflt(sgo,"offset" ,0.0);
    else
      pa       = onode_getpar_flt_dflt(sgo,"poisson_offset" ,0.0);

    if (onode_get_item_child(sgo,"toffset") != NULL){
      ptoff    = onode_getpar_flt_dflt(sgo,"toffset",0.0);
      ptoff *= 1000.0;  // Convert to msec
    }else
      ptoff    = onode_getpar_flt_dflt(sgo,"poisson_toff_ms",0.0);

    refmu  = onode_getpar_flt_dflt(sgo,"refract_mean",0.0);
    refsig = onode_getpar_flt_dflt(sgo,"refract_sd",0.0);

    dumpflag = onode_getpar_flt_dflt(sgo,"spike_dump",0);

    autoscale = onode_getpar_flt_dflt(sgo,"autoscale",-1.0);

    //  Additive noise
    if (onode_get_item_child(sgo,"noise_sd") != NULL){
      noise_sd  = onode_getpar_flt_exit(sgo,"noise_sd");
      noise_tsd = onode_getpar_flt_exit(sgo,"noise_tsd");
      noise_mu  = onode_getpar_flt_dflt(sgo,"noise_mu",0.0);
    }else{
      noise_sd = noise_mu = 0.0;
    }

  }

  if (dumpflag)
    append_farray_plot(DUMP_FNAME,"Raw_Response",d,n,1);

  // WYETH - This is an inelegant way to handle time sampling conversion
  if (my_rint(sampling) == 1000){
    t = copy_farray(d,n); // Copy the data
    tn = n;
  }else if (my_rint(sampling) == 500){
    t = over_sample_farray(d,n,2);
    tn = 2*n;
  }else{
    printf("  SAMPLING = %f\n",sampling);
    exit_error("IFC_UTIL_POISSON","Sampling not implemented for this tscale.");
  }

  add_const_farray(t,tn,pa0);  // pre-scaling offset
  multiply_farray(t,tn,px);    // scaling
  add_const_farray(t,tn,pa);   // post-scaling offset

  if (noise_sd > 0.0){
    //
    //  *** WYETH ADDED April 17, 2015 for Jacob's multi-channel model
    //  *** NOTE, The same seed is used for this noise, as for the Poisson
    //     noise, which means there might be some undesired correlation
    //   *****************
    aseed = seed / 2 + 113;  // *** WYETH - BIG HACK
    noise = gaussian_corr_noise(noise_tsd,noise_mu,noise_sd,tn,aseed);
    add_to_farray(t,noise,tn);
    myfree(noise);
  }


  half_wave_rectify_farray(t,tn);  // Set negative values to 0

  if (autoscale > 0.0){
    //
    //  WYETH - Added Apr 30, 2014.  Should we consider the case where the
    //          trace is all negative, i.e., should we set mean to zero before
    //          scaling?
    //
    area = sum_farray(t,tn,0,tn);
    if (area > 0.0){
      totspikes = tn / 1000.0 * autoscale; // (sec) * (spikes/sec) = spikes
      norm_area_farray(t,tn,totspikes);    // Scale area to desired spike count
    }
  }

  nflag = 0;
  for(i=0;i<tn;i++)  // Set values > 1 back to 1
    if (t[i] > 1.0){
      t[i] = 1.0;
      nflag += 1;
    }
  if (nflag > 0){
    sprintf(ggstr,"IFC_UTIL_POISSON  Prob > 1, %d times\n",nflag);
    mylog(mylogf,ggstr);
  }

  if (seed > 0)
    seed = -seed;  // Seed needs to be negative to initialize

  // WYETH - sampling is forced to be 1000, thus spikes in msec

  if ((refmu == 0.0) && (refsig == 0.0)){
    make_poisson_float_spikes_from_prob(t,tn,1000.0,&seed,&fs,&ns);
  }else{
    //  WYETH ADDED May 2011
    make_poisson_f_spikes_from_pr_refr(t,tn,1000.0,&seed,refmu,refsig,
				       &fs,&ns);
  }

  add_const_farray(fs,ns,ptoff);

  if (dumpflag){
    int *ss;
    //printf("SAMPLING = %f\n",sampling);  ***  500.0 when tn=0.002 
    append_farray_plot(DUMP_FNAME,"Firing_Prob",t,tn,1);
    ss = f2iarray(fs,ns);
    append_xplot_spike_sampling(DUMP_FNAME,"Spikes",ss,ns,0,tn,0,1000.0,0);
    myfree(ss);
  }

  *rfs = fs;
  *rfn = ns;

  if (vflag){  // If a float value is requested, return the spike prob
    *rv = t;
    *rvn = tn;
  }else{
    myfree(t);  // Free the copied data
    *rv = NULL;
    *rvn = -1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            IFC_UTIL_RUN_IFC_01                            */
/*                                                                           */
/*****************************************************************************/
void ifc_util_run_ifc_01(m,s,r,k,action)
     struct model_struct *m;    // Model params.
     struct stim_struct *s;     // Stimulus params.
     struct response_struct *r; // Response params.
     int k;                     // Trial number, when action = 1
     int action;                // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  mylog(mylogf,"  IFC_UTIL_RUN_IFC_01\n");

  if (action == 0){
    mod_ifc_01_prep(m,s);     // Prep
  }else if (action == 1){
    model_ifc_01_get_response(m,s,r,k);
  }else if (action == -1){
    ; // No cleanup?
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            IFC_UTIL_RUN_IFC_02                            */
/*                                                                           */
/*****************************************************************************/
void ifc_util_run_ifc_02(m,s,r,k,action)
     struct model_struct *m;    // Model params.
     struct stim_struct *s;     // Stimulus params.
     struct response_struct *r; // Response params.
     int k;                     // Trial number, when action = 1
     int action;                // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  mylog(mylogf,"  IFC_UTIL_RUN_IFC_02\n");

  if (action == 0){
    ; // No prep?
  }else if (action == 1){
    model_ifc_02_get_response(m,s,r,k);
  }else if (action == -1){
    ; // No cleanup?
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            IFC_UTIL_RUN_IFC_02_02                         */
/*                                                                           */
/*****************************************************************************/
void ifc_util_run_ifc_02_02(m,s,r,k,action)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     int k;                     // Trial number, when action = 1
     int action;                // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  mylog(mylogf,"  IFC_UTIL_RUN_IFC_02_02\n");

  if (action == 0){
    ; // No prep?
  }else if (action == 1){
    model_ifc_02_02_get_response(m,s,r,k);
  }else if (action == -1){
    ; // No cleanup?
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             IFC_UTIL_RUN_IFC_03                           */
/*                                                                           */
/*****************************************************************************/
void ifc_util_run_ifc_03(m,s,r,k,action)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     int k;                     // Trial number, when action = 1
     int action;                // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  mylog(mylogf,"  IFC_UTIL_RUN_IFC_03\n");

  if (action == 0){
    ; // No prep?
  }else if (action == 1){
    model_ifc_03_get_response(m,s,r,k);
  }else if (action == -1){
    ; // No cleanup?
  }
}
