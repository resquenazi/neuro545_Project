/*****************************************************************************/
/*                                                                           */
/*  mod_me_util.c                                                            */
/*  wyeth bair                                                               */
/*                                                                           */
/*  Motion Energy model.                                                     */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */

//  Problem:  This  'mod_me_fpo_s' is not set for ds_simp_01_get_response(
//   and it does get set in:  'model_gabor_comp_prep(m)'
//     But,     fs[i]->s = matrix(1,xn,1,2*yn);

/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "my_util.h"
#include "misc_util.h"
#include "nr_util.h"
#include "iarray_util.h"
#include "farray_util.h"
#include "carray_util.h"
#include "myrand_util.h"
#include "data_util.h"
#include "paramfile_util.h"
#include "fft_util.h"
#include "sig_util.h"
#include "kernel_util.h"
#include "stim_util.h"
#include "mod_util.h"
#include "pop_cell_util.h"  // WYETH - ADDED FOR MT POP
#include "mod_conn_util.h"  // WYETH - ADDED FOR MT POP
#include "pop_util.h"       // WYETH - ADDED FOR MT POP
#include "ifc_util.h"
#include "mod.h"            // Data structures
#include "ifc.h"            // Data structures
#include "paramfile.h"      // Data structures

int ORIG_FLAG = 0;  // Used for testing new code

// For cluster computation  (static needed w/ initialization
static int myid = -1;         // Process ID, -1-Single, 0-Master, >0-Slave
static char *mylogf = NULL;   // Logfile name
char ggstr[LONG_SLEN];        // Global log string

//
//  Filter structure
//
struct mod_me_fstruct{
  struct onode *fo;            // Onode related to this filter (or NULL)
  struct param_pair_list *ppl; // Par. list related to this filter (or NULL)
  int flag;         // 0-do not process; 1-process
  float ***f;       // [xn][yn][tn] filter
  float  **s;       // [][] for FFT
  float ***r;       // [xn][yn][tn] response
  float normc;      // Normalization constant for this filter
};

struct mod_me_fstop{
  int n;                         // Number of filters
  int flag_re;                   // 1= construct separate RE filters
  struct mod_me_fstruct  **fsl;  // [n] Filter structure list
  struct param_pair_list **fpl;  // [n] Filter param list

  int nchan;                     // Number of channels, e.g., STF
  int ndir;                      // Number of directions
  int nph;                       // Number of phases
  int ***fndx;                   // [chan][dir][phase] Index into 'fsl'

  //
  //  WYETH ADAPT - could put filter table data, if we switch the old way
  //    to use this "new way" structure.
  //
};

struct mod_me_fstop *mod_me_ftop = NULL;  // Top level filter structure


//
//  Noise
//
struct mod_me_noise{
  char *type;      // 'gfg'
  float mu;        // Mean
  float sd;        // SD
  float tsd;       // Temporal SD (raw time units)
  int seed;        // Seed for seedlist generation
  int sflag;       // Seed flag, starts at 1, gets set to 0 once seed used
  int nseed;       // Number of seed values (s->ntr)
  int *seedlist;   // [nseed] List of seed values for noise generation
  int currseed;    // Currently active seed value
};

//
//  Global for multi-filter model
//
struct mod_mem{
  int nf;             // Number of filter sets
  float *theta;       // [nf] direction
  float *ssd;         // [nf] spatial SD 
  float *sf;          // [nf] spatial frequency
  float *tsd;         // [nf] temporal SD
  float *tf;          // [nf] temporal frequency
  float *scale;       // [nf] scaling factor for filter response

  // WYETH - Can replace with arrays of ...fstruct above
  float *****f;       // Filters [nf][0..3][xn][yn][tn]
  float ****fs;       // FT of filters
  float *****r;       // Response, kept for efficient repeats

  float wf;           // Factor to control how quickly weights adapt

  int dumpflag;       // 0-no dump, 1-dump weights etc...

  int nr;             // Number of responses
  int *xi,*yi;        // [nr] coordinates of responses
  float ***mp;        // [nf][nr][tn] pref ME
  float ***ma;        // [nf][nr][tn] anti ME
};

struct mod_mem *mod_me_mem; // Global for multiple filters


// Reichardt model
float ***mod_reich_lp;      // Low-pass filter
float ***mod_reich_hp;      // High-pass filter
float  **mod_reich_lp_s;    // For FT of filters
float  **mod_reich_hp_s;    // High-pass filter
float ***mod_reich_rlp = NULL;   // response, kept for efficient repeats
float ***mod_reich_rhp = NULL;   // response, kept for efficient repeats
int      mod_reich_uxn;     // Number of units along x-dimension
int      mod_reich_uyn;     // Number of units along y-dimension
int      mod_reich_ux0;     // Offset of first unit along x-dimension
int      mod_reich_uy0;     // Offset of first unit along y-dimension
int      mod_reich_dx;      // Offset in spatial sampling units
int      mod_reich_dy;      // Offset in spatial sampling units
int      mod_reich_dump;    // Dump intermediate traces

// RD_2Gabor
int      mod_rd_2gabor_dti;  // Delay of "even" filter (so dir 0 is pref)
int      mod_rd_2gabor_rect; // 1-rectify before multiplication

//
//  Filter list (2012 Dec)  WYBINOC
//

int mod_me_fre = 0;                        // 1= construct separate RE filters
int mod_me_fn = 0;                         // Number of filters
struct mod_me_fstruct **mod_me_ff = NULL;  // [fn] List of pointers to filters
float **mod_me_ftdata;                     // [ndir]

char *mod_me_binoc_rsign = NULL;   // '++++', '++--', etc
int   mod_me_binoc_srect;          // 1-rectify monoc. simp. response, 0-do not
float mod_me_binoc_sthresh;        // threshold, subtracted before simp rect
                                   //   fraction of average filter (rect'd) area
float mod_me_binoc_athresh;        // actual threshold value
char *mod_me_binoc_nonlin = NULL;  // 'halfsq', 'square' or 'none'


int      me_xflag = 0;      // 1-extra filter (odd, even, extra), thus 3 filts.

// Global for me model
float ***mod_me_fpe = NULL; // ME filters, quadrature opponent
float ***mod_me_fpo = NULL;
float ***mod_me_fpx = NULL; // WYETH fxxx
float ***mod_me_fne = NULL;
float ***mod_me_fno = NULL;
float ***mod_me_fnx = NULL; // WYETH fxxx
float ***mod_me_fg0 = NULL; // Normalization, or gain, filters
float ***mod_me_fg1 = NULL;

float  **mod_me_fpe_s; // For FT of filters
float  **mod_me_fpo_s;
float  **mod_me_fpx_s; // WYETH fxxx
float  **mod_me_fne_s;
float  **mod_me_fno_s;
float  **mod_me_fnx_s; // WYETH fxxx
float  **mod_me_fg0_s;
float  **mod_me_fg1_s;

float ***mod_me_rpe = NULL;   // response, kept for efficient repeats
float ***mod_me_rpo = NULL;
float ***mod_me_rpx = NULL; // WYETH fxxx
float ***mod_me_rne = NULL;
float ***mod_me_rno = NULL;
float ***mod_me_rnx = NULL; // WYETH fxxx
float ***mod_me_rg0 = NULL;
float ***mod_me_rg1 = NULL;

float mod_me_sscale;
float mod_me_tscale;
int   mod_me_xn;
int   mod_me_yn;
int   mod_me_tn;
char *mod_me_modtype;           // model type
char *mod_me_filt;              // 'even' or 'odd', used by 'ds_simp_01'

float mod_me_tdelay    = 0.0;   // Additive time delay for spikes (s)
int   mod_me_oppflag   = 0;     // Opponent flag, *** 0 MEANS OPPONENT ***
float mod_me_opp_w     = 1.0;   // Opponent weight, i.e., scale anti signal

float mod_me_pow       = 1.0;   // Power following ME squaring

int   mod_me_normflag  = 0;     // Whether to normalize
int   mod_me_norm_app  = 1;     // Where to apply normalization
float mod_me_norm_pow  = 1.0;   // What power to use, after squaring
float mod_me_norm_c    = 1.0;   // Additive constant
float mod_me_norm_f    = 1.0;   // Multiplicative factor
float mod_me_norm_xn   = 1;     // Spatial integration, side length
float mod_me_norm_tsd  = 0.0;   // Temporal smoothing, Gaussian SD (sec)
float mod_me_norm_toff = 0.0;   // Temporal offset for norm filters

int mod_me_stimno_prev = -2;    // Index of previous stimulus

//
//  ME_V5 model
//
struct mod_me_v1resp{  // n=dir  nc=even/odd
  float ***raw;   // [n][nc][tn] Raw V1 resp., like Ln or Ln^2 in Rust2006
  float  **me;    // [n][tn] ME computed from raw signals
  float   *sum;   // [tn] Sum of raw V1 response
  float ***norm;  // [n][nc][tn] Normalized V1 resp., like Vn in Rust2006
  float ***opp;   // [n][nc*2][tn] Opponent V1 responses
  float ***bde;   // [n][7][tn] BDE signals, even and odd, ON and OFF

  // Adaptation state and altered response
  float  **adst;   // [nPre][nc] Adaptation state variable
                   // For "s01" - the recent mean activity
                   // For "RPHomeo" - the weights from each dir ???
  float ***adrsp;  // [nPre][nc][tn] Adapted value

  // WYETH MULT-DISP for multiple disparity, might need to use 'nc' to cover
  //   all phases, and then allow 'bde' to have another dimension, which would
  //   be multiple disparity channels.
  //

  float   *mt;    // [tn] Raw MT response, like Q(t) in RustEtAl2006
};

struct mod_me_v1_input{
  int                   flag;   // 0-no storage, 1-L, 2-L,R
  struct mod_me_v1resp *dl;     // Responses, Left stream
  struct mod_me_v1resp *dr;     // Responses, Right stream
  float                *mtb;    // [tn] MT response combined from L,R streams
  // WYETH MULT-DISP for multiple disparity, might need another dimension
  //   on 'mtb'
};

struct mod_me_v5_pop{
  float ****ssum;      // [mxn][myn][mzn][tn] Spatially summed MT response
  float ****nl;        // [mxn][myn][mzn][tn] MT nonlinear response
  int    ***seed;      // [mxn][myn][mzn] seeds for spike generation
  int       nsamp;     // Number of V1 inputs (stacks) per MT unit
  int       stackflag; // 1-stacked V1 inputs, 0-random dir vs. position
  int       zdirflag;  // 1-zn maps onto number of dir chans.

  // WYETH STF new way
  struct onode *wco;   // Weight configuration onode

  float   *w;          // [n] Weights for each channel
  float   *w_r;        // [n] Weights for each channel, RIGHT EYE
  float    r_bias;     // Weight of R.Eye in binocular combination (1 = full)
  float    negwf;      // Scale factor for negative weights

  char    *nl_type;    // Type of nonlin, "rmsm", "half-rect", "none"
  float    nl_a;       // param A
  float    nl_b;       // param B

  char    *inorm_type;   // Type of norm: "pow_norm"
  float    inorm_pow;    // power parameter
  int      inorm_pow1st; // Apply power *before* linear weight - "no subunit"
  char    *inorm_subu;   // Subunit: "stack_dir" "none"
  char    *inorm_subwi;  // Within-subunit: 'avg', 'sum', ...

  //  Create weights for inhibitory MT population
  int      mtsub_ndir;  // Number of direction channels
  float   *mtsub_w;     // [..._ndir] This should match z-axis of '..._name'
  char    *mtsub_name;  // Name of MT population, with z-axis over '..._ndir'
};

// WYETH STF
struct mod_me_v1_spatgrp{  // Spatial group

  int nx;        // horizontal positions
  int ny;        // vertical positions
  struct mod_me_v1_input ***bdgr;  // [nx][ny] Binocular direction group
                                   //    could be full (xn,yn) or (1,1)
};

//  WYETH STF
//
//  The STF channels could have been organized as a 2D array, TF by SF, but I
//  have decided to organize them as a simple 1D list where arbitrary
//  parameters can vary.  This will be useful if we want to vary other
//  paramters (e.g., SDs) or if we want a dependency between the SF and TF
//  ranges.
//
struct mod_me_v1_stf{  // STF channels
  int     n;           // Number of STF channels
  char   *config;      // Configuration:  "SFxTF"
  int     npar;        // Number of *index* parameters, e.g., 2
  int     nptot;       // Totoal number of params: 'npar' + extra linked pars
  char  **par_name;    // [nptot ar] Name of parameter, e.g., 'SF', 'TF', "veli"
  int    *nval;        // [nptot ar] Number of values for each parameter
  char ***par_val;     // [nptot ar][nval] unique parameter values
  int   **vndx;        // [n][npar] *index* par value indices for each STF chan

  int     flag_config; // 0-TF is directly available
                       // 1-TF must be looked up from dir+vel

  int   **flag_dir;    // [n][ndir] Is there a filter built for this direction?
                       // 0-No filter
                       // 1-yes
                       // 2-yes, but it is redundant w/ another

  struct mod_me_v1_spatgrp **sgr;  // [n] spatial group for each STF channel

  // *** WYETH STF  HERE
  float ***sum_norm;     // [xn][yn][tn] global normalization sum
  float ***sum_norm_r;   // [xn][yn][tn] global normalization sum, R.E.
};

struct mod_me_v1_surr{  // V1 surround
  int       n;          // Number of surround signals
  float ****sl;         // [n][xn][yn][tn] surround signals, Left eye
  float ****sr;         // [n][xn][yn][tn] surround signals, Right eye
  float     sig;        // SD for surround spatial integration
  float     sigc;       // SD inner diam of optional ANNULAR surround, or -1
  float     ampl;       // Amplitude, e.g., for subtraction
  float     omix;       // 0-DS, 1-non DS  (opponent mixing flag)
  float     delay;      // (s) delay of signal
  char     *rule;       // 'subtract'
  int       dumpi;      // [-1] Index of 3D signal to dump, -1 for none
};

struct mod_me_v5_struct{
  int     binoc;    // 0-left eye, 1-left and right eyes
  int     n;        // Number of direction channels 'ndir'
  float  *w;        // [n] Weights for each channel
  float  *w_r;      // [n] Weights for each channel, RIGHT EYE
  float   alph1;    // param 1  // V1 normalization params
  float   alph2;    // param 2
  float   alph3;    // param 3
  float   alph4;    // param 4
  float   L;        // "mean squared stimulus contrast"

  struct onode *normo;   // Normalization parameter object (pointer)

  float   v1_opp_w;    // Weight of opponent subtraction (1 = full opponent)
  float   v1_opp_tmin; // No opp. for STF chans w/ TF below this value (Hz)
  int     v1_recto;    // 1-rectify after opponency
  float   v1_bin_w;    // Weight of R.E. in binocular combination (1 = full)
  float   v1_mix_f;    // Fraction of other eye to mix (0 to 0.5) [0]
  char   *v1_mix_s;    // Stage of mixing: "linear" ...
  char   *mt_nltype;   // WYETH - ADDED FOR COMPAT with new '_pop' param
  float   mt_a;        // param a  // MT nonlinearity params
  float   mt_b;        // param b
  float   bin_hack;    // Hack for binoc model of Tailby et al. % neg MT weights
  int     v1_dflag;    // 0-no disparity, 1-disparity
  int     v1_dnorm;    // 0-abs, 1-ME
  int     v1_nc;       // Number of v1 signals per dir chan (e.g. 2: odd,even)

  // Adaptation params
  struct onode *v1ado; // Adaptation parameter object (pointer)
  int     v1_adflag;   // 0-no adaptation, 1-adaptation, 2-freeze adapt
  char   *v1_adtype;   // Adaptation type
  float   v1_adp_t1;   // Adaptation time scale param
  float   v1_adp_a1;   // Adaptation strength param
  float  *v1_adp_tv;   // [ndir/2 + 1] Target values for RPHomeo
  char   *v1_adrw;     // 'read', 'write' or 'null'
  char   *v1_adsfile;  // state file name
  int     v1_adp_wwi;  // "write_wieghts" increment (time steps)
  //
  //  Wyeth: if there is surround suppression AND adaptation, then we need
  //    to run adaptation on *all* grid points.  Thus, we cannot put the state
  //    in the 'v1resp' struct, because it only gets built for units that we
  //    need to compute later responses for.
  //
  float ***v1_adst_l;    // [n][xn][yn] adaptation state (surr ==> dense grid)
  float ***v1_adst_r;    // [n][xn][yn] adaptation state (surr ==> dense grid)
  float ****v1_ad_normw; // [xn][yn][ndir][ndir] Norm weight for RPHomeo
                         // These weights all start at some constant, and
                         //   change over time.

  struct mod_me_v1_surr *v1srnd;  // Surround suppression data

  // WYETH STF  new way
  struct mod_me_v1_stf *v1stf;  // SF-TF channels

  // WYETH OLD WAY HERE - no STF channels
  struct mod_me_v1_input ***v1resp;  // [xn][yn] or [1][1] Input signals

  float  *mt_norm;  // [tn] Norm'd MT response, like M(t) in RustEtAl2006

  int mtli0;           // Index of first MT layer in 'mod_me_mpt'
  int mt_npop;         // Number of MT pops, should relate to the value in
                       //  'mod_me_mpt->nlay', e.g. less by 2.

  struct mod_me_v5_pop **mtp;  // [mt_npop] Responses for each MT pop

  struct mod_me_noise *noise_v1_out;

  float **ftdata;      // [ndir][3] Filter table data, or NULL
};

struct mod_me_v5_struct *mod_me_v5 = NULL;
struct pop_top *mod_me_mpt;          // Global pop model structure


int mod_me_tsi;     // Trial seed index, copy of r->tsi
                    // **** BE SURE  ...get_response sets this if needed ****

/**************************************-**************************************/
/*                                                                           */
/*                            MOD_ME_UTIL_BAR_RESP                           */
/*                                                                           */
/*  Return the product of the bar with the map.                              */
/*                                                                           */
/*****************************************************************************/
float mod_me_util_bar_resp(xt,xn,tn,bx0,bxn,bt0)
     float **xt;      // RF map [xn][tn]
     int xn,tn;       // dimensions of 'xt'
     int bx0,bxn;     // location and width of bar
     int bt0;         // time location of bar (assume btn is 1)
{
  int i;
  int bx1;
  float r;

  if ((bt0 < 0) || (bt0 >= tn))
    return 0.0;

  bx1 = bx0 + bxn - 1;  // End point

  if (bx0 < 0)
    bx0 = 0;
  if (bx1 >= xn)
    bx1 = xn - 1;

  r = 0.0;
  for(i=bx0;i<=bx1;i++)
    r += xt[i][bt0];

  return r;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_ME_UTIL_BAR_SIM                           */
/*                                                                           */
/*  Determine he two-ar                                                      */
/*                                                                           */
/*****************************************************************************/
void mod_me_util_bar_sim(xte,xto,xn,tn)
     float **xte,**xto;
     int xn,tn;
{
  int i,j,k;
  int dn,ntot,dt,bwid,dx,dx2,bx0;
  float *xdata,*tww,*tbb,*twb,*tbw,b1,b2,b1o,b2o;

  write_2d_data("zzz.e.2d",xte,0,0,xn,tn,4,2,0,0);
  write_2d_data("zzz.o.2d",xto,0,0,xn,tn,4,2,0,0);

  // Set number of displacement values
  dn = 2*xn - 1;

  dt = 5;    // Time separation of bars (tscale units)
  bwid = 16;  // Bar width (sscale units)

  tww = get_zero_farray(dn);
  tbb = get_zero_farray(dn);
  twb = get_zero_farray(dn);
  tbw = get_zero_farray(dn);

  xdata = get_zero_farray(dn);
  
  ntot = 0;
  for(i=0;i<dn;i++){  // For each DX
    dx = i - (dn-1)/2;  // Offset of bars in spatial units
    dx2 = (int)(dx/2);  // Half of DX
    xdata[i] = (float)dx;
    for(j=0;j<xn;j++){  // For each point in space
      bx0 = j-dx2;  // origin of bar
      for(k=(tn/2 - dt/2-1);k<=(tn/2 - dt/2);k++){  // For each point in time
	//for(k=0;k<tn;k++){  // For each point in time

	b1 = mod_me_util_bar_resp(xte,xn,tn,bx0,bwid,k);
	b2 = mod_me_util_bar_resp(xte,xn,tn,bx0+dx,bwid,k+dt);

	b1o = mod_me_util_bar_resp(xto,xn,tn,bx0,bwid,k);
	b2o = mod_me_util_bar_resp(xto,xn,tn,bx0+dx,bwid,k+dt);

	//b1 = 0.0;  // Response to bar 1
	//b2 = 0.0;  // Response to bar 2
	
	tww[i] += ( b1 + b2) * ( b1 + b2);
	tbb[i] += (-b1 - b2) * (-b1 - b2);
	twb[i] += ( b1 - b2) * ( b1 - b2);
	tbw[i] += (-b1 + b2) * (-b1 + b2);

	/*
	tww[i] += ( b1o + b2o) * ( b1o + b2o);
	tbb[i] += (-b1o - b2o) * (-b1o - b2o);
	twb[i] += ( b1o - b2o) * ( b1o - b2o);
	tbw[i] += (-b1o + b2o) * (-b1o + b2o);
	*/

	ntot += 1;
      }
    }
  }

  append_farray_xy_plot("zz.tune.pl",xdata,tww,dn,"W_W");
  append_farray_xy_plot("zz.tune.pl",xdata,tbb,dn,"B_B");
  append_farray_xy_plot("zz.tune.pl",xdata,twb,dn,"W_B");
  append_farray_xy_plot("zz.tune.pl",xdata,tbw,dn,"B_W");
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_ME_GET_QUAD_OPP_K                           */
/*                                                                           */
/*  Construct filter from params for "me_k_..." where 'k' is 0..             */
/*                                                                           */
/*****************************************************************************/
void mod_me_get_quad_opp_k(mo,fo,rpe,rpo,rne,rno,rnc)
     struct onode *mo;
     struct onode *fo;
     float ****rpe,****rpo,****rne,****rno;
     float *rnc;  // Normalization constant
{
  int xn,yn,tn,write_filt;
  float ssd,tsd,sf,tf,theta;
  char *tfilt;

  xn = mod_me_xn;
  yn = mod_me_yn;
  tn = mod_me_tn;

  theta = onode_getpar_flt_exit(fo,"direction");
  sf    = onode_getpar_flt_exit(fo,"sf");
  ssd   = onode_getpar_flt_exit(fo,"ssd");

  tfilt = onode_getpar_chr_exit(mo,"me_tfilt");

  if ((strcmp(tfilt,"gaussian")==0)||(strcmp(tfilt,"maxwell")==0)){
    tsd   = onode_getpar_flt_exit(fo,"tsd");
    tf    = onode_getpar_flt_exit(fo,"tf");

    get_quad_opponent_gabor(xn,yn,tn,mod_me_sscale,
			    mod_me_tscale,ssd,tsd,sf,tf,theta,0.0,0,tfilt,
			    rpe,rpo,rne,rno);

    write_filt = onode_getpar_int_dflt(fo,"write_filter",-1);
    //name = onode_getpar_chr_dflt(fo,"name",NULL);
    if ((write_filt == 1)||(write_filt == 5))
      write_3d_data_part("zzz.fpe.3d",*rpe,1,xn,1,yn,1,tn,4,2,1);
    if ((write_filt == 2)||(write_filt == 5))
      write_3d_data_part("zzz.fpo.3d",*rpo,1,xn,1,yn,1,tn,4,2,1);
    if ((write_filt == 3)||(write_filt == 5))
      write_3d_data_part("zzz.fne.3d",*rne,1,xn,1,yn,1,tn,4,2,1);
    if ((write_filt == 4)||(write_filt == 5))
      write_3d_data_part("zzz.fno.3d",*rno,1,xn,1,yn,1,tn,4,2,1);


    // This normalization constant accounts for changes in the size of the
    // filter based on the temporal and spatial extent, i.e., SDs.
    *rnc = ssd*10.0 * ssd*10.0 * tsd*10.0;
  }else{
    mylog_exit(mylogf,"MOD_ME_GET_QUAD_OPP_K  unknown 'me_tfilt' value");
  }

  myfree(tfilt);
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_ME_GET_QUAD_OPP_FILTERS                       */
/*                                                                           */
/*****************************************************************************/
void mod_me_get_quad_opp_filters(mppl,rpe,rpo,rne,rno)
     struct param_pair_list *mppl;
     float ****rpe,****rpo,****rne,****rno;
{
  int n1,n2,xn,yn,tn;
  int wrflag,tilt_flag;
  float k,ssd,tsd,sf,tf,sscale,tscale,theta,shift_frac,ph;
  char *tfilt;
  
  sscale = paramfile_get_float_param_or_exit(mppl,"sscale");
  tscale = paramfile_get_float_param_or_exit(mppl,"tscale");
  xn = paramfile_get_int_param_or_exit(mppl,"xn");
  yn = paramfile_get_int_param_or_exit(mppl,"yn");
  tn = paramfile_get_int_param_or_exit(mppl,"tn");

  tilt_flag = paramfile_get_int_param_default(mppl,"me_tilt",0);
  if (tilt_flag == 1)
    shift_frac = paramfile_get_float_param_default(mppl,"me_tilt_frac",1.0);

  ph = paramfile_get_float_param_default(mppl,"me_phase_shift",90.0);

  mod_me_sscale = sscale;
  mod_me_tscale = tscale;
  mod_me_xn = xn;
  mod_me_yn = yn;
  mod_me_tn = tn;
  
  tfilt = NULL; // Default value
  if ((strcmp(mod_me_modtype,"gabor_comp")==0)){
    ;
  }else if ((strcmp(mod_me_modtype,"me_ab_01")==0)){
    theta = paramfile_get_float_param_or_exit(mppl,"me_theta");
    sf = paramfile_get_float_param_or_exit(mppl,"me_sf");
    ssd = paramfile_get_float_param_or_exit(mppl,"me_ssd");
    exit_error("MOD_ME_GET_QUAD_OPP_FILTERS",
	       "WYETH - DISCONTINUED - USE .moo file");
  }else{
    theta = paramfile_get_float_param_or_exit(mppl,"me_theta");
    sf = paramfile_get_float_param_or_exit(mppl,"me_sf");
    ssd = paramfile_get_float_param_or_exit(mppl,"me_ssd");
    //tfilt = paramfile_get_char_param_or_exit(mppl,"me_tfilt");
    tfilt = paramfile_get_char_param_default(mppl,"me_tfilt","gaussian");
  }

  if ((strcmp(mod_me_modtype,"me_gabor_01")==0) ||
      (strcmp(mod_me_modtype,"ds_liv3")==0) ||
      (strcmp(mod_me_modtype,"ds_simp_01")==0)){
    tsd = paramfile_get_float_param_or_exit(mppl,"gabor_tsd");
    tf = paramfile_get_float_param_or_exit(mppl,"gabor_tf");
    if (tilt_flag == 0){
      if (ph != 90.0){
	get_arb_opponent_gabor(xn,yn,tn,sscale,tscale,ssd,tsd,sf,tf,theta,
			       tfilt,ph,rpe,rpo,rne,rno);
      }else{
	get_quad_opponent_gabor(xn,yn,tn,sscale,tscale,ssd,tsd,sf,tf,theta,0.0,
				0,tfilt,rpe,rpo,rne,rno);
      }
    }else if (tilt_flag == 1){ // Use shifting method
      get_quad_opp_tilt_gabor(xn,yn,tn,sscale,tscale,ssd,tsd,sf,tf,theta,tfilt,
			      shift_frac,0,rpe,rpo,rne,rno);
    }else if (tilt_flag == 2){  // Use rotation method
      get_quad_opp_tilt_gabor(xn,yn,tn,sscale,tscale,ssd,tsd,sf,tf,theta,tfilt,
			      0.0,1,rpe,rpo,rne,rno);
    }
  }else if (strcmp(mod_me_modtype,"me_ab_01")==0){
    n1 = paramfile_get_int_param_or_exit(mppl,"me_ab_n1");
    n2 = paramfile_get_int_param_or_exit(mppl,"me_ab_n2");
    k = paramfile_get_float_param_or_exit(mppl,"me_ab_k");
    wrflag = paramfile_get_int_param_default(mppl,"model_write_tfilt",0);
    quad_opponent_causal_stf3d_01_tensor(xn,yn,tn,sscale,tscale,ssd,ssd,sf,
					 n1,n2,k,theta,rpo,rpe,rno,rne,wrflag);
  }else if (strcmp(mod_me_modtype,"gabor_comp")==0){
    get_quad_gabor(mppl,NULL,xn,yn,tn,sscale,tscale,rpe,rpo);
  }else
    exit_error("ME_MOD_GET_QUAD_OPP_FILTERS","Unknown model type");

  if (tfilt != NULL)
    myfree(tfilt);
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_ME_FS_CREATE                             */
/*                                                                           */
/*****************************************************************************/
struct mod_me_fstruct **mod_me_fs_create(fn)
     int fn;
{
  int i;
  struct mod_me_fstruct **fs;

  // Create filter list 'fs'
  fs = (struct mod_me_fstruct **)myalloc(fn*sizeof(struct mod_me_fstruct *));
  for(i=0;i<fn;i++){
    fs[i] = (struct mod_me_fstruct *)myalloc(sizeof(struct mod_me_fstruct));
    fs[i]->flag = 1;
    fs[i]->fo  = NULL;
    fs[i]->ppl = NULL;
    fs[i]->f = NULL;
    fs[i]->s = NULL;
    fs[i]->r = NULL;
    fs[i]->normc = -1.0;
  }

  return fs;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_ME_FS_TOP_CREATE                           */
/*                                                                           */
/*****************************************************************************/
void mod_me_fs_top_create(fn)
     int fn;
{
  int i;
  struct mod_me_fstop *fts;
  struct mod_me_fstruct **fs;

  //
  //  Top level filter structure
  //
  fts = (struct mod_me_fstop *)myalloc(sizeof(struct mod_me_fstop));
  fts->n       = fn;
  fts->flag_re = 0;
  fts->fpl     = NULL;

  fts->nchan = -1;
  fts->ndir  = -1;
  fts->nph   = -1;
  fts->fndx  = NULL;

  //
  //  Filter list
  //
  fts->fsl = mod_me_fs_create(fn);
  /*  WYETH REMOVE THIS, commented out on 2018 Apr 3
  fs = (struct mod_me_fstruct **)myalloc(fn*sizeof(struct mod_me_fstruct *));
  for(i=0;i<fn;i++){
    fs[i] = (struct mod_me_fstruct *)myalloc(sizeof(struct mod_me_fstruct));
    fs[i]->flag = 1;
    fs[i]->fo  = NULL;
    fs[i]->ppl = NULL;
    fs[i]->f = NULL;
    fs[i]->s = NULL;
    fs[i]->r = NULL;
    fs[i]->normc = -1.0;
  }
  fts->fsl = fs;
  */

  mod_me_ftop = fts;  // Save to global top
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_ME_FS_FREE_R                             */
/*                                                                           */
/*  Free response storage, typically before getting response to next stim.   */
/*                                                                           */
/*****************************************************************************/
void mod_me_fs_free_r(fs,fn)
     struct mod_me_fstruct **fs;
     int fn;
{
  int i;

  // ****** WYETH - can get rid of input argument once the old 'fs' way is
  //                      changed ????

  for(i=0;i<fn;i++)  // For each filter
    if (fs[i]->r != NULL)
      free_f3tensor(fs[i]->r,1,mod_me_xn,1,mod_me_yn,1,mod_me_tn);
}
/**************************************-**************************************/
/*                                                                           */
/*                               MOD_ME_FS_FREE                              */
/*                                                                           */
/*****************************************************************************/
void mod_me_fs_free(fs,n)
     struct mod_me_fstruct **fs;  // [n] list of filter structs
     int n;
{
  int i;
  int xn,yn,tn;

  xn = mod_me_xn;
  yn = mod_me_yn;
  tn = mod_me_tn;

  for(i=0;i<n;i++){

    if (fs[i]->fo  != NULL)
      ; // WYETH - do we need to free anything?  Or is it a pointer

    if (fs[i]->ppl != NULL)
      ; // WYETH - do we need to free anything?  Or is it a pointer

    if (fs[i]->f != NULL){
      free_f3tensor(fs[i]->f,1,xn,1,yn,1,tn);
    }

    if (fs[i]->s != NULL){
      free_matrix(fs[i]->s,1,xn,1,2*yn);
    }

    if (fs[i]->r != NULL){
      free_f3tensor(fs[i]->r,1,xn,1,yn,1,tn);
    }

    myfree(fs[i]);
  }

  myfree(fs);
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_ME_FS_FLAG_EYE                            */
/*                                                                           */
/*  Set the flags for the specified eye to 1, and the other to 0.            */
/*                                                                           */
/*****************************************************************************/
void mod_me_fs_flag_eye(fs,fn,eye)
     struct mod_me_fstruct **fs;
     int fn;
     int eye;  // 0-left, 1-right
{
  int i;
  int f_le,f_re;

  if (eye == 0){
    f_le = 1;
    f_re = 0;
  }else{
    f_le = 0;
    f_re = 1;
  }

  for(i=0;i<fn;i++){     // For each filter
    if (i < fn/2)
      fs[i]->flag = f_le;   // LE filters
    else
      fs[i]->flag = f_re;   // RE filters
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_ME_FS_GET_FI                             */
/*                                                                           */
/*  Return the filter index.                                                 */
/*                                                                           */
/*****************************************************************************/
int mod_me_fs_get_fi(eye,stfi,di,pi)
     int eye;    // 0-left 1-right
     int stfi;   // STF channel index, i.e., in 'mod_me_v1_stf' struct
     int di;     // Direction index
     int pi;     // Phase index
{
  int fi;

  if ((eye < 0) || (eye > 1))
    exit_error("MOD_ME_FS_GET_FI","Invalid 'eye'");

  if ((stfi < 0) || (stfi >= mod_me_ftop->nchan))
    exit_error("MOD_ME_FS_GET_FI","Invalid 'stfi'");

  if ((di < 0) || (di >= mod_me_ftop->ndir))
    exit_error("MOD_ME_FS_GET_FI","Invalid 'di'");

  if ((pi < 0) || (pi >= mod_me_ftop->nph))
    exit_error("MOD_ME_FS_GET_FI","Invalid 'pi'");
    

  fi = mod_me_ftop->fndx[stfi][di][pi];

  if (eye == 1){
    if (mod_me_ftop->flag_re == 1)  
      fi += mod_me_ftop->n / 2;
    // Otherwise, the RE uses the same filter as the LE
  }

  return fi;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_ME_FS_BUILD_LIST                           */
/*                                                                           */
/*****************************************************************************/
void mod_me_fs_build_list()
{
  int i;
  int xn,yn,tn,nf;
  float ss,ts;
  struct mod_me_fstruct **fsl;
  struct param_pair_list **fpl;

  mylog(mylogf,"  MOD_ME_FS_BUILD_LIST\n");

  xn = mod_me_xn;
  yn = mod_me_yn;
  tn = mod_me_tn;
  ss = mod_me_sscale;
  ts = mod_me_tscale;

  nf  = mod_me_ftop->n;
  fsl = mod_me_ftop->fsl;
  fpl = mod_me_ftop->fpl;

  mylog(mylogf,"   ");  // Space before filter number list

  for(i=0;i<nf;i++){
    fsl[i]->f     = get_ppl_gabor(fpl[i],xn,yn,tn,ss,ts);
    fsl[i]->normc = -1.0;
    fsl[i]->ppl   = fpl[i];

    /*** WRITE OUT RAW FILTERS
    {
      char fname[SLEN];
      sprintf(fname,"zzz.%d",i);
      write_3d_data_part(fname,fsl[i]->f,1,xn,1,yn,1,tn,4,2,1);
    }
    ***/

    if ((i > 0) && (i%10 == 0)){
      sprintf(ggstr," %d",i);
      mylog(mylogf,ggstr); fflush(stdout);
    }
  }
  sprintf(ggstr,"\n    Created %d filters\n",nf);
  mylog(mylogf,ggstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_ME_GET_FILT_WORK                           */
/*                                                                           */
/*****************************************************************************/
void mod_me_get_filt_work(fo,xn,yn,tn,sscale,tscale,sdo,sf,theta,tf_off,
			  disp_off,disp_adj,tfilt,ph,tilt_flag,shift_frac,
			  rpe,rpo,rne,rno)
     struct onode *fo;
     int xn,yn,tn;
     float sscale,tscale;
     float sdo,sf,theta;
     float tf_off;        // Offset for TF (for variation in right eye)
     float disp_off;      // Disparity phase offset, all 4 quadrature filters
     int   disp_adj;      // 1-correct disparity for opposite direction motion
     char *tfilt;
     float ph;            // 90-Quadrature, !=90 is other than quadrature
     char *tilt_flag;
     float shift_frac;
     float ****rpe,****rpo,****rne,****rno;
{
  char *ftype;
  float tsd,tf;

  /*
    printf("  xn,yn,tn = %d %d %d\n",xn,yn,tn);
    printf("  sscale,tscale = %f %f\n",sscale,tscale);
    printf("  sdo,sf,theta = %f %f %f\n",sdo,sf,theta);
    printf("  tf_off = %f\n",tf_off);
    printf("  disp_off = %f\n",disp_off);
    printf("  disp_adj = %d\n",disp_adj);
    printf("  tfilt = %s\n",tfilt);
    printf("  ph = %f\n",ph);
    printf("  tilt_flag = %s\n",tilt_flag);
    printf("  shift_frac = %f\n",shift_frac);
  */

  ftype = onode_getpar_chr_dflt(fo,"type","none");

  if (strcmp(ftype,"Gabor_curve_t")!=0){  // No TF for this curve
    tsd = onode_getpar_flt_exit(fo,"tsd");
    tf  = onode_getpar_flt_exit(fo,"tf");
  }

  if (strcmp(tilt_flag,"none")==0){

    if (ph != 90.0){
      get_arb_opponent_gabor(xn,yn,tn,sscale,tscale,sdo,tsd,sf,tf,theta,
			     tfilt,ph,rpe,rpo,rne,rno);
    }else{
      //
      //  Quadrature (90 deg) case
      //
      if (strcmp(ftype,"Gabor_curve_t")==0){
	get_quad_opponent_gabor_tshift(mylogf,fo,xn,yn,tn,sscale,tscale,
				       rpe,rpo,rne,rno);
      }else{
	//printf("wyeth here ________________tfilt = %s\n",tfilt);
	get_quad_opponent_gabor(xn,yn,tn,sscale,tscale,sdo,tsd,sf,tf+tf_off,
				theta,disp_off,disp_adj,tfilt,rpe,rpo,rne,rno);
      }
    }
  }else if (strcmp(tilt_flag,"shift")==0){ // Use shifting method
    get_quad_opp_tilt_gabor(xn,yn,tn,sscale,tscale,sdo,tsd,sf,tf,theta,
			    tfilt,shift_frac,0,rpe,rpo,rne,rno);
  }else if (strcmp(tilt_flag,"rotate")==0){ // Use rotation method
    get_quad_opp_tilt_gabor(xn,yn,tn,sscale,tscale,sdo,tsd,sf,tf,theta,
			    tfilt,0.0,1,rpe,rpo,rne,rno);
  }
  myfree(ftype);
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_ME_RESP_LIST_ADD                           */
/*                                                                           */
/*****************************************************************************/
void mod_me_resp_list_add()
{
  // WYETH HERE 2018 June;
  ;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_ME_PREP_BASIC                             */
/*                                                                           */
/*****************************************************************************/
void mod_me_prep_basic(mo,rss,rts,rxn,ryn,rtn)
     struct onode *mo;
     float *rss;        // Return sscale
     float *rts;        // Return tscale
     int *rxn;          // Return xn
     int *ryn;          // Return yn
     int *rtn;          // Return tn
{
  mod_me_modtype        = onode_getpar_chr_exit(mo,"mod_type");

  *rss = mod_me_sscale  = onode_getpar_flt_exit(mo,"sscale");
  *rts = mod_me_tscale  = onode_getpar_flt_exit(mo,"tscale");
  *rxn = mod_me_xn      = onode_getpar_int_exit(mo,"xn");
  *ryn = mod_me_yn      = onode_getpar_int_exit(mo,"yn");
  *rtn = mod_me_tn      = onode_getpar_int_exit(mo,"tn");

  if (power_of_two(mod_me_tn) < 0){
    exit_error("MOD_ME_PREP_BASIC","Value of 'tn' must be a power of 2");
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_ME_FILT_TEMPLATE                           */
/*                                                                           */
/*  Build synaptic connection templates for each of the 'zn' channels,       */
/*  which differ in direction and phase.                                     */
/*                                                                           */
/*****************************************************************************/
void mod_me_filt_template(stfi,rfx,rfy,rfw,rfn,rzn)
     int stfi;         // Index of STF channel
     int   ***rfx;     // [rzn][rfn] x-coords
     int   ***rfy;     // [rzn][rfn] y-coords
     float ***rfw;     // [rzn][rfn] weights
     int    **rfn;     // [rzn]      Number of coords
     int     *rzn;     // Number of templates
{
  int i,j,k;
  int di,pi,fi,zi,zn,xn,yn,tn,minx,miny,mint,maxx,maxy,maxt,ti,n,ndir,nph;
  int **tmplt,**fx,**fy,*fn;
  float ***f,vmin,vmax,z,amax,**fw,v11;
  struct mod_me_fstruct **fs;

  xn = mod_me_xn;
  yn = mod_me_yn;
  tn = mod_me_tn;

  if (mod_me_ftop == NULL){
    fs = mod_me_ff;         // OLD WAY - REMOVE when not used

    if (mod_me_fre == 1)  // There are separate filters for L and R eyes.
      zn = mod_me_fn/2;
    else
      zn = mod_me_fn;  // This is typically more than number of dir chans.

    ndir = zn;  // WYETH HACK - just make the outer loop run for all filters
    nph = 1;    // Inner loop will run once.
  }else{
    fs   = mod_me_ftop->fsl;   // New way
    ndir = mod_me_ftop->ndir;
    nph  = mod_me_ftop->nph;
    zn = ndir*nph;
  }

  fx = (int   **)myalloc(zn*sizeof(int *));
  fy = (int   **)myalloc(zn*sizeof(int *));
  fw = (float **)myalloc(zn*sizeof(float *));
  fn = (int    *)myalloc(zn*sizeof(int));

  zi = 0;
  for(di=0;di<ndir;di++){
    for(pi=0;pi<nph;pi++){

      if (mod_me_ftop == NULL) // old way 
	fi = zi;  // Use all filters in sequence
      else
	fi = mod_me_ftop->fndx[stfi][di][pi];  // Filter index

      f = fs[fi]->f;  // Pointer to filter

      //  Find 'ti' which is the time of maximum value in the 3D filter.
      //  Find 'amax' which is the maximum abs value at that time.
      //
      get_min_max_coord_3d_farray(f,1,xn,1,yn,1,tn,&minx,&miny,&mint,
				  &maxx,&maxy,&maxt);
      vmin = f[minx][miny][mint];
      vmax = f[maxx][maxy][maxt];

      if (vmax > -vmin){
	ti = maxt;
	amax = fabs(vmax);
      }else{
	ti = mint;
	amax = fabs(vmin);
      }

      //
      //  Built a template to mark the connections
      //
      tmplt = get_zero_2d_iarray(xn,yn);

      n = 0; // Number of entries in template
      for(i=0;i<xn;i++){
	for(j=0;j<yn;j++){
	  z = fabs(f[i+1][j+1][ti]) / amax;
	  if (z > 0.02){
	    tmplt[i][j] = 1;
	    n += 1;
	  }
	}
      }

      if (n == 0){
	printf("  *** Indices  dir %d   phase %d   stfi %d\n",di,pi,stfi);
	printf("  *** vmin vmax = %f %f\n",vmin,vmax);
	exit_error("MOD_ME_FILT_TEMPLATE","Zero entries in template");
      }

      fx[zi] = (int   *)myalloc(n*sizeof(int));
      fy[zi] = (int   *)myalloc(n*sizeof(int));
      fw[zi] = (float *)myalloc(n*sizeof(float));
      fn[zi] = n;

      k = 0;
      for(i=0;i<xn;i++){
	for(j=0;j<yn;j++){
	  if (tmplt[i][j] == 1){
	    fx[zi][k] = i;
	    fy[zi][k] = j;
	    v11 = f[i+1][j+1][ti] / amax;   // Values from -1 to 1
	    fw[zi][k] = 0.5*(1.0 + v11);    // Scale from 0 to 1
	    k += 1;
	  }
	}
      }
      free_2d_iarray(tmplt,xn);

      zi += 1;
    }
  }
  
  *rfx = fx;
  *rfy = fy;
  *rfw = fw;
  *rfn = fn;
  *rzn = zn;
}
/**************************************-**************************************/
/*                                                                           */
/*                        MOD_ME_GET_QUAD_OPP_FILTERS_O                      */
/*                                                                           */
/*****************************************************************************/
struct mod_me_fstruct **mod_me_get_quad_opp_filters_o(mo,rpe,rpo,rpx,rne,
						      rno,rnx,rfn)
     struct onode *mo;
     float ****rpe,****rpo,****rpx,****rne,****rno,****rnx;
     int *rfn;
{
  int i,j;
  int fn,n1,n2,xn,yn,tn;
  int wrflag,seed,phopp,qsopp,ndir,disp_adj;
  float k,sdo,sdp,tsd,sf,tf,sscale,tscale,theta,shift_frac,ph,nsd_r,nsd_t;
  float binoc_shift,phas,quadshift,pwr,phase1,tf_off,disp_off;
  char *tfilt,*tilt_flag,*ftype,*ftfile;
  struct onode *fo,*v1o;
  struct mod_me_fstruct **fs;  // [fn] List of pointers to filters

  sscale = mod_me_sscale;
  tscale = mod_me_tscale;
  xn = mod_me_xn;
  yn = mod_me_yn;
  tn = mod_me_tn;

  fs = NULL;

  phopp = onode_getpar_int_dflt(mo,"phase_shift_opp",0);
  qsopp = onode_getpar_int_dflt(mo,"quad_shift_opp",0);

  fo = onode_child_get_unique(mo,"filter");
  ftype = onode_getpar_chr_dflt(fo,"type","none");

  tilt_flag = onode_getpar_chr_dflt(fo,"tilt_type","none");
  if (strcmp(tilt_flag,"shift")==0)
    shift_frac = onode_getpar_int_dflt(fo,"tilt_frac",1.0);
  else
    shift_frac = 0.0;  // WYETH.ADDED


  // Note, this value is the phase shift for the "Quadrature" shifted
  //   filters, thus it should generally be 90, and it should only be changed
  //   to make filters that are *not* in quadrature.
  ph = onode_getpar_flt_dflt(fo,"phase_shift",90.0);


  if (qsopp == 0)
    quadshift = 90.0;
  else
    quadshift = -90.0;

  //
  //  Set parameters depending on filter type
  //  *** WYETH - THIS SHOULD BE ELIMINATED AS MUCH AS POSSIBLE,
  //      and the 'fo' should be sent deeper into the kernel utilities.
  //

  tfilt = NULL; // Default value
  if (strcmp(mod_me_modtype,"gabor_comp")==0){
    ;
  }else if (strcmp(mod_me_modtype,"rd_2gabor")==0){
    float dtsec;
      
    dtsec = onode_getpar_flt_exit(fo,"dt_even");  // (s)
    mod_rd_2gabor_dti = dtsec / mod_me_tscale;
    mod_rd_2gabor_rect = onode_getpar_int_exit(fo,"rectify");  // (s)

  }else if ((strcmp(mod_me_modtype,"me_ab_01")==0)){
    theta = onode_getpar_flt_exit(fo,"direction");
    sf = onode_getpar_flt_exit(fo,"sf");
    sdo = onode_getpar_flt_exit(fo,"sd_orth");
    sdp = onode_getpar_flt_exit(fo,"sd_par");
  }else if ((strcmp(mod_me_modtype,"random_filter")==0)){
    mylog(mylogf,"  (mod_me_util) - random filter\n");
    sdo = onode_getpar_flt_exit(fo,"ssd");
  }else if (strcmp(ftype,"Gabor_CWQ")==0){
    ;
  }else{

    if (strcmp(mod_me_modtype,"me_v5")==0)
      theta = 0.0;  // Not used for this model
    else
      theta = onode_getpar_flt_exit(fo,"direction");

    sf = onode_getpar_flt_exit(fo,"sf");

    if (onode_get_item_child(fo,"ssd") != NULL){
      sdo = onode_getpar_flt_exit(fo,"ssd");
      sdp = sdo;
    }else{
      sdo = onode_getpar_flt_exit(fo,"sd_orth");
      sdp = onode_getpar_flt_exit(fo,"sd_par");
    }
    tfilt = onode_getpar_chr_dflt(fo,"tfilt","gaussian");
  }

  //
  //  Get filters, and set 'fn'
  //
  if ((strcmp(mod_me_modtype,"binoc_filter")==0)){

    if (strcmp(ftype,"Gabor")==0){

      // WYETH Added Apr 2014, to allow different TF, thus tilt, in the
      // right eye
      // *** WYETH - PAMELA WANTS THIS TO WORK FOR MT MODEL
      tf_off = onode_getpar_flt_dflt(mo,"tf_offset",0.0);

      // Phase shift for right eye
      binoc_shift = onode_getpar_flt_exit(mo,"phase_shift");

      tsd = onode_getpar_flt_exit(fo,"tsd");
      tf  = onode_getpar_flt_exit(fo,"tf");
      phase1  = onode_getpar_flt_dflt(fo,"phase",0.0);

      // pref = left
      phas = phase1;
      *rpe = gabor_space_time_tensor(xn,yn,tn,sscale,tscale,sdo,sdp,tsd,
				     sf,tf,theta,phas,tfilt);
      phas = phase1 + quadshift; // WYETH use 'ph' (above) for 90?
      *rpo = gabor_space_time_tensor(xn,yn,tn,sscale,tscale,sdo,sdp,tsd,
				     sf,tf,theta,phas,tfilt);
      // null = right
      phas = phase1 + binoc_shift;
      *rne = gabor_space_time_tensor(xn,yn,tn,sscale,tscale,sdo,sdp,tsd,
				     sf,tf+tf_off,theta,phas,tfilt);

      if (phopp == 0){
	phas = phase1 + binoc_shift + quadshift;
      }else if (phopp == 1){
	phas = phase1 - binoc_shift + quadshift;
      }else
	exit_error("ME_MOD_GET_QUAD_OPP_FILTERS_O","Bad phase_shift_opp val");
	
      *rno = gabor_space_time_tensor(xn,yn,tn,sscale,tscale,sdo,sdp,tsd,
				     sf,tf+tf_off,theta,phas,tfilt);
    }else if (strcmp(ftype,"Gabor_CWQ")==0){

      bde_cwq_stf3d_tensor(fo,mo,xn,yn,tn,sscale,tscale,rpe,rne,rpo,rno);

      //exit_error("WYETH HERE  ** (mod_me_util)","Testing");
    }
    fn = 4;
  }else if ((strcmp(mod_me_modtype,"me_gabor_01")==0) ||
	    (strcmp(mod_me_modtype,"ds_liv3")==0) ||
	    (strcmp(mod_me_modtype,"ds_simp_01")==0)){
    mod_me_get_filt_work(fo,xn,yn,tn,sscale,tscale,sdo,sf,theta,0.0,0.0,0,
			 tfilt,ph,tilt_flag,shift_frac,rpe,rpo,rne,rno);
    fn = 4;
  }else if (strcmp(mod_me_modtype,"random_filter")==0){

    nsd_r = onode_getpar_flt_exit(fo,"noise_ssd");
    nsd_t = onode_getpar_flt_exit(fo,"noise_tsd");
    seed  = onode_getpar_int_exit(fo,"noise_seed");

    tsd = onode_getpar_flt_exit(fo,"tsd");

    // The "pref" filters are independent random, the anti-pref are
    // mirrored on the x-axis relative to the pref.
    // 'rpo' is a copy of 'rpe'

    if (me_xflag == 0){
      get_opponent_random(xn,yn,tn,sscale,tscale,sdo,tsd,nsd_r,nsd_t,seed,
			  rpe,rpo,rne,rno);
      fn = 4;
    }else{
      get_opponent_random_x(xn,yn,tn,sscale,tscale,sdo,tsd,nsd_r,nsd_t,seed,
			    rpe,rpo,rpx,rne,rno,rnx);
      fn = 6;
    }
  }else if (strcmp(mod_me_modtype,"me_ab_01")==0){
    n1 = onode_getpar_int_exit(fo,"ab_n1");
    n2 = onode_getpar_int_exit(fo,"ab_n2");
    k  = onode_getpar_flt_exit(fo,"ab_k");
    wrflag = onode_getpar_int_dflt(fo,"write_tfilt",0);
    quad_opponent_causal_stf3d_01_tensor(xn,yn,tn,sscale,tscale,sdo,sdp,sf,
					 n1,n2,k,theta,rpo,rpe,rno,rne,wrflag);
    fn = 4;
  }else if ((strcmp(mod_me_modtype,"gabor_comp")==0) ||
	    (strcmp(mod_me_modtype,"rd_2gabor")==0)){
    /*
      exit_error("ME_MOD_GET_QUAD_OPP_FILTERS_O",
      "This part not converted to ONODE yet.");
    */
    get_quad_gabor(NULL,fo,xn,yn,tn,sscale,tscale,rpe,rpo);
    fn = 2;
  }else if (strcmp(mod_me_modtype,"me_v5")==0){

    //
    // WYETH STF:  Call a routine to generate a list of filter params
    //
    //    dir sf sdo sdp tf 

    v1o = onode_child_get_unique(mo,"v1");
    ndir = onode_getpar_int_exit(v1o,"n_dir_chan");

    tf_off = onode_getpar_flt_dflt(v1o,"tf_offset",0.0);  // For right eye
    if (tf_off != 0.0)
      mod_me_fre = 1;  // Global flag to indicate separate RE filters

    //
    //  WYETH - Oct 21, 2016 - Added for Pamela to change disparity tuning.
    //
    if (onode_item(v1o,"disparity_offset") == 1){
      disp_off = onode_getpar_flt_dflt(v1o,"disparity_offset",0.0); // (deg) RE
      disp_adj = onode_getpar_int_dflt(v1o,"disparity_offset_adj",1);
      mod_me_fre = 1;  // Global flag to indicate separate RE filters
      sprintf(ggstr,"    R.E. filter phase offset = %f\n",disp_off);
      mylog(mylogf,ggstr);
    }else{
      disp_off = 0.0;  // (deg)
      disp_adj = 0;    // (deg)
    }

    if (mod_me_fre == 0)
      fn = ndir * 2;  // Make one set of filters: they are same in both eyes
    else{
      fn = ndir * 4;  // Make filters for both eyes
    }
    fs = mod_me_fs_create(fn);


    //
    //  If there is a filter table, set 'ftfile' to the file name
    //
    ftfile = onode_getpar_chr_dflt(v1o,"filter_table",NULL);
    if (ftfile != NULL){
      if (strcmp(ftfile,"NULL")==0){
	myfree(ftfile);
	ftfile = NULL;
      }
    }

    if (ftfile == NULL){
      //
      //  Build filters the usual way
      //
      for(i=0;i<ndir/2;i++){
	theta = (float)i * 360.0 / (float)ndir;
	//printf("  %2d  %f  %f\n",i,sdo,theta);
	mod_me_get_filt_work(fo,xn,yn,tn,sscale,tscale,sdo,sf,theta,0.0,
			     0.0,0,tfilt,ph,tilt_flag,shift_frac,
			     rpe,rpo,rne,rno);
	fs[i*4  ]->f = *rpe;
	fs[i*4+1]->f = *rpo;
	fs[i*4+2]->f = *rne;
	fs[i*4+3]->f = *rno;
      }
    }else{
      //
      //  WYETH ADAPT Read some filter parameters from a file
      //
      int nrec,ncol;
      char ***ftdata;

      ftfile = onode_getpar_chr_exit(v1o,"filter_table");
      sprintf(ggstr,"  Reading parameters from filter table:  %s\n",ftfile);
      mylog(mylogf,ggstr);

      read_column_data(ftfile,0,&ftdata,&nrec,&ncol);
      mod_me_v5->ftdata = get_2d_farray(ndir,ncol);

      for(i=0;i<nrec;i++){
	sprintf(ggstr,"    %2d  ",i);
	mylog(mylogf,ggstr);
	for(j=0;j<ncol;j++){
	  sprintf(ggstr," %s",ftdata[i][j]);
	  mylog(mylogf,ggstr);
	  mod_me_v5->ftdata[i][j] = atof(ftdata[i][j]);
	}
	mylog(mylogf,"\n");
      }

      //  Columns are:
      //    (1) direction (deg)
      //    (2) Spatial SD (deg)
      //    (3) Amplitude [0..1]

      for(i=0;i<ndir;i++){
	//theta = (float)i * 360.0 / (float)ndir;
	theta = mod_me_v5->ftdata[i][0];
	sdo   = mod_me_v5->ftdata[i][1];
	//printf("  %2d  %.2f  %f\n",i,theta,sdo);
	mod_me_get_filt_work(fo,xn,yn,tn,sscale,tscale,sdo,sf,theta,0.0,
			     0.0,0,tfilt,ph,tilt_flag,shift_frac,
			     rpe,rpo,rne,rno);
	if (i < ndir/2){
	  fs[i*4  ]->f = *rpe;
	  fs[i*4+1]->f = *rpo;
	  free_f3tensor(*rne,1,xn,1,yn,1,tn);
	  free_f3tensor(*rno,1,xn,1,yn,1,tn);
	}else{
	  j = i - ndir/2;
	  fs[j*4+2]->f = *rpe;
	  fs[j*4+3]->f = *rpo;
	  free_f3tensor(*rne,1,xn,1,yn,1,tn);
	  free_f3tensor(*rno,1,xn,1,yn,1,tn);
	}
      }
    }

    if (mod_me_fre == 1){  // Make filters for Right Eye

      //
      //  The LE filters will be from 0...ndir*2, then RE filters will follow.
      //
      for(i=0;i<ndir/2;i++){
	// Note, "direction" param is ignored
	theta = (float)i * 360.0 / (float)ndir;

	mod_me_get_filt_work(fo,xn,yn,tn,sscale,tscale,sdo,sf,theta,tf_off,
			     disp_off,disp_adj,tfilt,ph,tilt_flag,shift_frac,
			     rpe,rpo,rne,rno);

	fs[i*4   + ndir*2]->f = *rpe;
	fs[i*4+1 + ndir*2]->f = *rpo;
	fs[i*4+2 + ndir*2]->f = *rne;
	fs[i*4+3 + ndir*2]->f = *rno;
      }
    }

  }else
    exit_error("MOD_ME_GET_QUAD_OPP_FILTERS_O","Unknown model type");

  if (tfilt != NULL)
    myfree(tfilt);
  if (ftype != NULL)  // I think it is always true
    myfree(ftype);
  if (tilt_flag != NULL)
    myfree(tilt_flag);

  if (fs == NULL){
    //
    //  If 'fs' not created yet, create and fill it in
    //
    fs = mod_me_fs_create(fn);
    if (fn == 2){
      //printf("***************** here filling it in to 'fs'\n");
      fs[0]->f = *rpe;
      fs[1]->f = *rpo;
    }else if (fn == 4){
      fs[0]->f = *rpe;
      fs[1]->f = *rpo;
      fs[2]->f = *rne;
      fs[3]->f = *rno;
    }else if (fn == 6){
      fs[0]->f = *rpe;
      fs[1]->f = *rpo;
      fs[2]->f = *rpx;
      fs[3]->f = *rne;
      fs[4]->f = *rno;
      fs[5]->f = *rnx;
    }else{
      exit_error("MOD_ME_GET_QUAD_OPP_FILTERS_O","Unknown 'fn' value");
    }
  }

  *rfn = fn;
  return fs;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_ME_GET_NORM_FILTERS                         */
/*                                                                           */
/*****************************************************************************/
void mod_me_get_norm_filters(mppl,o,rg0,rg1)
     struct param_pair_list *mppl;
     struct onode *o;
     float ****rg0,****rg1;
{
  int xn,yn,tn,tabn,normtype,pflag,hilbert,t0flag;
  float sscale,tscale,sig1,sig2,amp1,amp2,tsig,tabk,csdt,csx,csy,tcent;
  float ***g0,***g1,sump,sumn,c,sum,toff;
  char *outfile;

  normtype = paramfile_get_int_param_or_exit(mppl,"me_norm_type");

  sscale = paramfile_get_float_param_or_exit(mppl,"sscale");
  tscale = paramfile_get_float_param_or_exit(mppl,"tscale");
  xn = paramfile_get_int_param_or_exit(mppl,"xn");
  yn = paramfile_get_int_param_or_exit(mppl,"yn");
  tn = paramfile_get_int_param_or_exit(mppl,"tn");

  sig1 = paramfile_get_float_param_or_exit(mppl,"me_dog_sig1");
  sig2 = paramfile_get_float_param_or_exit(mppl,"me_dog_sig2");
  amp1 = paramfile_get_float_param_or_exit(mppl,"me_dog_amp1");
  amp2 = paramfile_get_float_param_or_exit(mppl,"me_dog_amp2");
  tsig = paramfile_get_float_param_or_exit(mppl,"me_dog_tsig");
  tabk = paramfile_get_float_param_or_exit(mppl,"me_dog_tab_k");
  tabn =   paramfile_get_int_param_or_exit(mppl,"me_dog_tab_n");
  csdt = paramfile_get_float_param_or_exit(mppl,"me_dog_cs_delay");

  if (normtype == 1){
    csx = csy = 0.0;
    tcent = 0.0;      // Center of multiplying Gaussian
    outfile = NULL;
    pflag = 0;
    t0flag = 1;       // Time origin is in middle of array
    hilbert = 0;
    toff = mod_me_norm_toff;

    exit_error("MOD_ME_GET_NORM_FILTERS","WYETH - broken, can't send tabk/n");
    //
    //  WYETH - 'dog_space_time_tensor' now picks up params from 'mppl' or 'o'
    //

    /*
    g0 = dog_space_time_tensor(xn,yn,tn,sscale,tscale,sig1,sig2,amp1,amp2,
			       csx,csy,csdt,tsig,tcent,tabk,tabn,
			       outfile,mppl,o,pflag,hilbert,t0flag,toff,0,0);
    hilbert = 1;
    g1 = dog_space_time_tensor(xn,yn,tn,sscale,tscale,sig1,sig2,amp1,amp2,
			       csx,csy,csdt,tsig,tcent,tabk,tabn,
			       outfile,mppl,o,pflag,hilbert,t0flag,toff,0.0);
    */

    /***  SEE NORMALIZATION BELOW
	kernel_util_norm_xyt_filter(g0,1,xn,1,yn,1,tn,&sump,&sumn);
	sprintf(ggstr,"    Norm filter g0, pos. and neg. sums:  %.4f  %.4f\n",
	sump,sumn);
	mylog(mylogf,ggstr);
	kernel_util_norm_xyt_filter(g1,1,xn,1,yn,1,tn,&sump,&sumn);
	sprintf(ggstr,"    Norm filter g1, pos. and neg. sums:  %.4f  %.4f\n",
	sump,sumn);
	mylog(mylogf,ggstr);
    ***/

    mod_me_fg0 = g0;
    mod_me_fg1 = g1;
  }else
    exit_error("ME_MOD_GET_NORM_FILTERS","Unknown me_norm_type");


  // Same normalization as used for ME filters
  c = (sscale/0.1)*(sscale/0.1)*(tscale/0.002);
  sum = sum_square_3d_farray(g0,1,xn,1,yn,1,tn);
  multiply_3d_farray(g0,1,xn,1,yn,1,tn,sqrt(c/sum));
  sum = sum_square_3d_farray(g1,1,xn,1,yn,1,tn);
  multiply_3d_farray(g1,1,xn,1,yn,1,tn,sqrt(c/sum));

}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_ME_NORMALIZE_FILTER                         */
/*                                                                           */
/*****************************************************************************/
void mod_me_normalize_filter(f,x0,xn,y0,yn,z0,zn,sscale,tscale,norm_const)
     float ***f;
     int x0,xn,y0,yn,z0,zn;
     float sscale,tscale;
     float norm_const;  // Based on spatial and temporal SDs
{
  float sum,c;

  c = (sscale/0.1)*(sscale/0.1)*(tscale/0.002);
  c /= norm_const;

  sum = sum_square_3d_farray(f,x0,xn,y0,yn,z0,zn);
  multiply_3d_farray(f,x0,xn,y0,yn,z0,zn,sqrt(c/sum));
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_ME_NORMALIZE_FILTERS                         */
/*                                                                           */
/*****************************************************************************/
void mod_me_normalize_filters(x0,xn,y0,yn,z0,zn,sscale,tscale,norm_const,ntype)
     int x0,xn,y0,yn,z0,zn;
     float sscale,tscale;
     float norm_const;       // Arbitrary scale factor, or OLD
                             //    Based on spatial and temporal SDs, or '-1.0'
     int ntype;              // 1 - divide by the sum of absolute value
                             // 2 - divide by the sum of square
{
  int i,k;
  int nf;
  float c,nc,sum,fscale;
  struct mod_me_fstruct **ff;

  if (mod_me_ftop == NULL){  // OLD WAY - REMOVE when not used
    nf = mod_me_fn;
    ff = mod_me_ff;
  }else{
    nf = mod_me_ftop->n;     // New way
    ff = mod_me_ftop->fsl;
  }

  c = (sscale/0.1)*(sscale/0.1)*(tscale/0.002);

  for(i=0;i<nf;i++){

    //  WYETH 2017 Jan:  New way (ntype 1) below does NOT use 'nc'
    if (norm_const == -1.0){
      nc = c / ff[i]->normc;  // Use value stored for this filter
    }else{
      nc = c / norm_const;
    }

    if (ntype == 1){
      sum = (float)sum_abs_3d_farray(ff[i]->f,x0,xn,y0,yn,z0,zn);
      fscale = norm_const/sum;
      //printf("fscale========== %f\n",fscale);
    }else if (ntype == 2){
      sum = (float)sum_square_3d_farray(ff[i]->f,x0,xn,y0,yn,z0,zn);
      fscale = sqrt(nc/sum);
    }else
      mylog_exit(mylogf,"  MOD_ME_NORMALIZE_FILTERS  bad 'ntype'");


    if (mod_me_v5 != NULL){
      if (mod_me_v5->ftdata != NULL){
	if (nf/2 != mod_me_v5->n){
	  printf("nf = %d  ndir = %d\n",nf,mod_me_v5->n);
	  mylog_exit(mylogf,"  MOD_ME_NORMALIZE_FILTERS  nf/2 != ndir");
	}

	if ((i%4)/2 == 0)
	  k = i/4;
	else
	  k = i/4 + mod_me_v5->n/2;

	fscale *= mod_me_v5->ftdata[k][2];  // WYETH ADAPT amplitude
	//printf("__________scaling by %f\n",mod_me_v5->ftdata[k][2]);
      }
    }

    //multiply_3d_farray(ff[i]->f,x0,xn,y0,yn,z0,zn,sqrt(nc/sum));
    multiply_3d_farray(ff[i]->f,x0,xn,y0,yn,z0,zn,fscale);

    /*** WRITE OUT normalized FILTERS ***
    {
      char fname[SLEN];
      sprintf(fname,"zzz.%d",i);
      write_3d_data_part(fname,ff[i]->f,1,xn,1,yn,1,zn,4,2,1);
    }
    ***/

    {
      double dsum,dsum2;

      dsum = sum_abs_3d_farray(ff[i]->f,x0,xn,y0,yn,z0,zn);
      dsum2 = sum_square_3d_farray(ff[i]->f,x0,xn,y0,yn,z0,zn);

      sprintf(ggstr,"    Filter %d norm:  %f  sum %f  nc %f  NewSum %lf  NewSum^2 %lf\n",i,sqrt(c/sum),sum,nc,dsum,dsum2);
      mylog(mylogf,ggstr);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_ME_FILTERS_FREE                            */
/*                                                                           */
/*****************************************************************************/
void mod_me_filters_free()
{
  int nf;
  struct mod_me_fstruct **ff;

  mylog(mylogf,"  MOD_ME_FILTERS_FREE\n");
	
  if (mod_me_ftop == NULL){  // OLD WAY - REMOVE when not used
    nf = mod_me_fn;
    ff = mod_me_ff;
  }else{
    nf = mod_me_ftop->n;     // New way
    ff = mod_me_ftop->fsl;

    if (mod_me_ftop->fndx != NULL){
      free_3d_iarray(mod_me_ftop->fndx,
		     mod_me_ftop->nchan,
		     mod_me_ftop->ndir,
		     mod_me_ftop->nph);
    }
  }

  mod_me_fs_free(ff,nf);
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_ME_BINOC_SET_THRESH                          */
/*                                                                           */
/*****************************************************************************/
void mod_me_binoc_set_thresh()
{
  int i;
  double dsum;

  mylog(mylogf,"  MOD_ME_BINOC_SET_THRESH\n");

  if (mod_me_fn == 0)  // WYETH - this added with new filter structure
    exit_error("MOD_ME_BINOC_SET_THRESH","No filters (old way?)");

  if (mod_me_binoc_sthresh == 0.0){
    mod_me_binoc_athresh = 0.0;
  }else{
    dsum = 0.0;
    for(i=0;i<mod_me_fn;i++){
      dsum += sum_abs_3d_farray(mod_me_ff[i]->f,1,mod_me_xn,1,mod_me_yn,
				1,mod_me_tn);
    }
    mod_me_binoc_athresh = (float)(mod_me_binoc_sthresh * dsum / mod_me_fn);
  }

  sprintf(ggstr,"    Threshold fraction of total filter area:  %f\n",
	  mod_me_binoc_sthresh);
  mylog(mylogf,ggstr);
  sprintf(ggstr,"    Threshold actual value:  %f\n",mod_me_binoc_athresh);
  mylog(mylogf,ggstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                             MODEL_ME_01_PREP_O                            */
/*                                                                           */
/*****************************************************************************/
void model_me_01_prep_o(m)
     struct model_struct *m; // Model params
{
  int xn,yn,tn,write_filt,fn;
  float ***fpe,***fpo,***fpx,***fne,***fno,***fnx,sscale,tscale;
  char *oppstr;
  struct param_pair_list *mppl; // Model parameter pair list
  struct onode *fo,*sgen;
  struct mod_me_fstruct **fs;  // [fn] List of pointers to filters

  mylog(mylogf,"  MODEL_ME_01_PREP_O\n");
  mylog(mylogf,"    Computing even and odd, pref and null, filters.\n");

  mod_me_prep_basic(m->o,&sscale,&tscale,&xn,&yn,&tn); // Get and set globals

  if (strcmp(mod_me_modtype,"binoc_filter")==0){
    mod_me_binoc_rsign   = onode_getpar_chr_dflt(m->o,"right_sign","++++");
    mod_me_binoc_srect   = onode_getpar_int_dflt(m->o,"simp_rect",0);
    mod_me_binoc_sthresh = onode_getpar_flt_dflt(m->o,"simp_thresh",0.0);
    mod_me_binoc_nonlin  = onode_getpar_chr_dflt(m->o,"binoc_nonlin","halfsq");
  }

  fs = mod_me_get_quad_opp_filters_o(m->o,&mod_me_fpe,&mod_me_fpo,
				     &mod_me_fpx,&mod_me_fne,
				     &mod_me_fno,&mod_me_fnx,
				     &fn); // Number in '..._ff'
  mod_me_ff = fs;
  mod_me_fn = fn;

  fpe = fs[0]->f; // mod_me_fpe;
  fpo = fs[1]->f; // mod_me_fpo;
  fne = fno = fpx = fnx = NULL;
  if (fn == 4){
    fne = fs[2]->f; // mod_me_fne;
    fno = fs[3]->f; // mod_me_fno;
  }else if (fn == 6){
    fpx = fs[2]->f; // mod_me_fpx;
    fne = fs[3]->f; // mod_me_fne;
    fno = fs[4]->f; // mod_me_fno;
    fnx = fs[5]->f; // mod_me_fnx;
  }

//mod_me_normalize_filters(1,xn,1,yn,1,tn,mod_me_sscale,mod_me_tscale,1.0);
  mod_me_normalize_filters(1,xn,1,yn,1,tn,sscale,tscale,1.0,2);

  if (strcmp(mod_me_modtype,"binoc_filter")==0){
    mod_me_binoc_set_thresh();
  }

  fo = onode_child_get_unique(m->o,"filter");
  mod_me_pow = onode_getpar_flt_dflt(fo,"exponent",1.0);

  // Normalization Filters
  mod_me_normflag = onode_getpar_int_dflt(fo,"norm_type",0);
  if (mod_me_normflag == 1){

    //printf("HERE____________________________******************\n");

    mod_me_norm_pow  = onode_getpar_flt_dflt(fo,"norm_pow",1.0);
    mod_me_norm_c    = onode_getpar_flt_dflt(fo,"norm_c",1.0);
    mod_me_norm_f    = onode_getpar_flt_dflt(fo,"norm_f",1.0);
    mod_me_norm_app  = onode_getpar_int_dflt(fo,"norm_app",1);
    mod_me_norm_xn   = onode_getpar_int_dflt(fo,"norm_xn",1);
    mod_me_norm_tsd  = onode_getpar_flt_dflt(fo,"norm_tsd",0.0);
    mod_me_norm_toff = onode_getpar_flt_dflt(fo,"norm_toff",0.0);

    //mylogx(mylogf,"MODEL_ME_01_PREP_O","This part not ONODE yet");
    mod_me_get_norm_filters(mppl,m->o,&mod_me_fg0,&mod_me_fg1);
  }

  //write_filt = paramfile_get_int_param_default(mppl,"model_write_filter",-1);
  write_filt = onode_getpar_int_dflt(fo,"write_filter",-1);
  if ((write_filt == 1)||(write_filt == 5))
    write_3d_data_part("zzz.fpe.3d",fpe,1,xn,1,yn,1,tn,4,2,1);
  if ((write_filt == 2)||(write_filt == 5))
    write_3d_data_part("zzz.fpo.3d",fpo,1,xn,1,yn,1,tn,4,2,1);
  if ((write_filt == 3)||(write_filt == 5))
    write_3d_data_part("zzz.fne.3d",fne,1,xn,1,yn,1,tn,4,2,1);
  if ((write_filt == 4)||(write_filt == 5))
    write_3d_data_part("zzz.fno.3d",fno,1,xn,1,yn,1,tn,4,2,1);
  if ((write_filt == 5) && (me_xflag == 1)){
    write_3d_data_part("zzz.fpx.3d",fpx,1,xn,1,yn,1,tn,4,2,1);
    write_3d_data_part("zzz.fnx.3d",fnx,1,xn,1,yn,1,tn,4,2,1);
  }

  if (write_filt == 6){
    write_3d_data_part("zzz.fg0.3d",mod_me_fg0,1,xn,1,yn,1,tn,4,2,1);
    write_3d_data_part("zzz.fg1.3d",mod_me_fg1,1,xn,1,yn,1,tn,4,2,1);
  }

  /***
  //
  // WYETH - use X-T plane of filter to test for tuning for bars
  //

  // Extract plane
  {
    int i,j;
    float **xte,**xto;

    xte = get_2d_farray(xn,tn);
    xto = get_2d_farray(xn,tn);
    for(i=0;i<xn;i++)
      for(j=0;j<tn;j++){
	xte[i][j]  = fpe[i+1][yn/2][j+1];
	xto[i][j]  = fpo[i+1][yn/2][j+1];
	//xte[i][j]  = fne[i+1][yn/2][j+1];
	//xto[i][j]  = fno[i+1][yn/2][j+1];
      }

    mod_me_util_bar_sim(xte,xto,xn,tn);

    printf("********** WYETH - for barsim work on .../an/twobar/mod/barsim\n");
    printf("********** WYETH - for barsim work on .../an/twobar/mod/barsim\n");
    exit(0);
  }

  ***/

  // Added delay, seconds
  //mod_me_tdelay = paramfile_get_float_param_default(mppl,"tdelay",0.0);
  sgen = onode_child_get_unique(m->o,"spike_gen");

  if (onode_test_chr(sgen,"type","poisson")){
    if (onode_get_item_child(sgen,"tdelay") != NULL){
      exit_error("WYETH","don't use 'tdelay', use toffset instead");
      //mod_me_tdelay = onode_getpar_flt_dflt(sgen,"tdelay",0.0);
    }
  }else{
    // WYETH - is this needed for ifc??
    mod_me_tdelay = onode_getpar_flt_dflt(sgen,"tdelay",0.0);
  }

  if (strcmp(mod_me_modtype,"ds_simp_01")==0){
    mod_me_oppflag = onode_getpar_int_dflt(fo,"opp_flag",3);
    mod_me_filt = onode_getpar_chr_dflt(fo,"branch","odd");
  }else if (strcmp(mod_me_modtype,"ds_liv3")==0){
    mod_me_oppflag = onode_getpar_int_dflt(fo,"opp_flag",3);
    mod_me_filt = onode_getpar_chr_dflt(fo,"filt","odd");
  }else{
    //
    //  *** 0 MEANS OPPONENT, FLAG IS INITIALIZED TO ZERO
    //
    oppstr = onode_getpar_chr_ptr(fo,"opponent");
    if (oppstr == NULL){
      mod_me_oppflag = 0; // Opponent
    }else if (strcmp(oppstr,"no")==0){
      mod_me_oppflag = 3; // Non-opponent
    }else if (strcmp(oppstr,"ei")==0){
      mod_me_oppflag = 2; // E-I opponency
    }else if (strcmp(oppstr,"yes")==0){
      mod_me_oppflag = 0; // Opponent
    }
    mod_me_opp_w = onode_getpar_flt_dflt(fo,"opponent_w",1.0);
    if (mod_me_opp_w != 1.0) // WYETH - added on Aug 24, 2016
      exit_error("MODEL_ME_01_PREP_O","The param 'opponent_w' does nothing");
    // *** NOTE, the above parameter should be removed if it does nothing
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MODEL_ME_01_PREP                             */
/*                                                                           */
/*****************************************************************************/
void model_me_01_prep(m)
     struct model_struct *m; // Model params
{
  int xn,yn,tn,write_filt;
  float ***fpe,***fpo,***fne,***fno;
  struct param_pair_list *mppl; // Model parameter pair list

  mylog(mylogf,"  MODEL_ME_01_PREP\n");
  mylog(mylogf,"    Computing even and odd, pref and null, filters.\n");


  mylog(mylogf,"  *** WYETH - THIS IS OLD WAY, consider using 'moo' file\n");
  mylog(mylogf,"  *** WYETH - THIS IS OLD WAY, consider using 'moo' file\n");
  mylog(mylogf,"  *** WYETH - THIS IS OLD WAY, consider using 'moo' file\n");

  mppl = m->ppl;

  mod_me_modtype = paramfile_get_char_param_or_exit(mppl,"mod_type");


  mod_me_get_quad_opp_filters(mppl,&mod_me_fpe,&mod_me_fpo,
			      &mod_me_fne,&mod_me_fno);
  xn = mod_me_xn;
  yn = mod_me_yn;
  tn = mod_me_tn;
  fpe = mod_me_fpe;
  fpo = mod_me_fpo;
  fne = mod_me_fne;
  fno = mod_me_fno;

  exit_error("WYETH HERE","Must put filters in a list, '..._ff'");
  mod_me_normalize_filters(1,xn,1,yn,1,tn,mod_me_sscale,mod_me_tscale,1.0,2);

  mod_me_pow = paramfile_get_float_param_default(mppl,"me_power",1.0);

  // Normalization Filters
  mod_me_normflag = paramfile_get_int_param_default(mppl,"me_norm_type",0);
  if (mod_me_normflag == 1){
    mod_me_norm_pow = paramfile_get_float_param_default(mppl,"me_norm_pow",
							1.0);
    mod_me_norm_c = paramfile_get_float_param_default(mppl,"me_norm_c",1.0);
    mod_me_norm_f = paramfile_get_float_param_default(mppl,"me_norm_f",1.0);
    mod_me_norm_app = paramfile_get_int_param_default(mppl,"me_norm_app",1);
    mod_me_norm_xn = paramfile_get_int_param_default(mppl,"me_norm_xn",1);
    mod_me_norm_tsd = paramfile_get_float_param_default(mppl,"me_norm_tsd",
							0.0);
    mod_me_norm_toff = paramfile_get_float_param_default(mppl,"me_norm_toff",
							 0.0);

    exit_error("MODEL_ME_01_PREP",
	       "WYETH these filters not supported, since change to fs?");

    mod_me_get_norm_filters(mppl,m->o,&mod_me_fg0,&mod_me_fg1);
  }

  write_filt = paramfile_get_int_param_default(mppl,"model_write_filter",-1);
  if ((write_filt == 0)||(write_filt == 4))
    write_3d_data_part("zzz.fpe.3d",fpe,1,xn,1,yn,1,tn,4,2,1);
  if ((write_filt == 1)||(write_filt == 4))
    write_3d_data_part("zzz.fpo.3d",fpo,1,xn,1,yn,1,tn,4,2,1);
  if ((write_filt == 2)||(write_filt == 4))
    write_3d_data_part("zzz.fne.3d",fne,1,xn,1,yn,1,tn,4,2,1);
  if ((write_filt == 3)||(write_filt == 4))
    write_3d_data_part("zzz.fno.3d",fno,1,xn,1,yn,1,tn,4,2,1);

  if (write_filt == 5){
    write_3d_data_part("zzz.fg0.3d",mod_me_fg0,1,xn,1,yn,1,tn,4,2,1);
    write_3d_data_part("zzz.fg1.3d",mod_me_fg1,1,xn,1,yn,1,tn,4,2,1);
  }

  /***
  //
  // WYETH - use X-T plane of filter to test for tuning for bars
  //

  // Extract plane
  {
    int i,j;
    float **xte,**xto;

    xte = get_2d_farray(xn,tn);
    xto = get_2d_farray(xn,tn);
    for(i=0;i<xn;i++)
      for(j=0;j<tn;j++){
	xte[i][j]  = fpe[i+1][yn/2][j+1];
	xto[i][j]  = fpo[i+1][yn/2][j+1];
	//xte[i][j]  = fne[i+1][yn/2][j+1];
	//xto[i][j]  = fno[i+1][yn/2][j+1];
      }

    mod_me_util_bar_sim(xte,xto,xn,tn);

    printf("********** WYETH - for barsim work on .../an/twobar/mod/barsim\n");
    printf("********** WYETH - for barsim work on .../an/twobar/mod/barsim\n");
    exit(0);
  }

  ***/

  // Added delay, seconds
  mod_me_tdelay = paramfile_get_float_param_default(mppl,"tdelay",0.0);


  if ((strcmp(mod_me_modtype,"ds_simp_01")==0) ||
      (strcmp(mod_me_modtype,"ds_liv3")==0)){
    mod_me_oppflag = paramfile_get_int_param_default(mppl,"opp_flag",3);
    mod_me_filt = paramfile_get_char_param_default(mppl,"me_filt","odd");
  }else{
    mod_me_oppflag = paramfile_get_int_param_default(mppl,"opp_flag",0);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_ME_NOISE_INIT                            */
/*                                                                           */
/*  Initialize values for noise from an onode.                               */
/*                                                                           */
/*****************************************************************************/
struct mod_me_noise *mod_me_noise_init(no,nseed)
     struct onode *no;
     int nseed;           // Number of seeds to create, or 0
{
  struct mod_me_noise *noise;

  noise = (struct mod_me_noise *)myalloc(sizeof(struct mod_me_noise));

  noise->type = onode_getpar_chr_exit(no,"type");
  noise->mu   = onode_getpar_flt_dflt(no,"mean",0.0);
  noise->sd   = onode_getpar_flt_exit(no,"sd");
  noise->tsd  = onode_getpar_flt_dflt(no,"tsd",0.0);
  noise->seed = onode_getpar_int_exit(no,"seed");
  noise->sflag = 1;  // 1-must initialize seed

  noise->tsd /= (mod_me_tscale * 1000.0);  // Convert from ms to tscale units

  noise->nseed = nseed;
  if (nseed > 0){
    noise->seedlist = get_seeds(noise->seed,100000,nseed);
  }else
    noise->seedlist = NULL;

  noise->currseed = 0;

  return noise;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_ME_NOISE_ADD                             */
/*                                                                           */
/*  Add noise to the data array 'd'.                                         */
/*                                                                           */
/*****************************************************************************/
void mod_me_noise_add(d,n,np)
     float *d;  // [n] data trace to add noise to
     int n;     // length of data trace
     struct mod_me_noise *np;  // Noise parameters
{
  int i;
  float *gt;
  struct mod_me_noise *noise;

  if (np->sflag == 1){  // If seed has not been used yet
    //
    //  Set the seed for this trial; negate to re-init random num gen.
    //
    np->currseed = np->seedlist[mod_me_tsi];

    if (np->currseed > 0)
      np->currseed *= -1;
    np->sflag = 0;      // Indicate that seed has been initialized
  }

  gt = gaussian_corr_noise_pseed(np->tsd,np->mu,np->sd,n,&(np->currseed));

  for(i=0;i<n;i++)
    d[i] += gt[i];

  myfree(gt);
}
/**************************************-**************************************/
/*                                                                           */
/*                           MODEL_ME_V5_RESP_INIT                           */
/*                                                                           */
/*****************************************************************************/
struct mod_me_v1resp *model_me_v5_resp_init(n,nc,tn)
     int n;     // Number of directions
     int nc;    // Number of channels per dir channel (e.g., 2 for even, odd)
     int tn;    // Time duration
{
  struct mod_me_v1resp *r;
  
  //
  //  Create response storage
  //

  // WYETH varmoo_bug - response storage created here for V1 REMOVE THIS LINE

  r = (struct mod_me_v1resp *)myalloc(sizeof(struct mod_me_v1resp));
  
  r->raw  = get_3d_farray(n,nc,tn);
  r->me   = get_2d_farray(n,tn);
  r->sum  = (float  *)myalloc(tn*sizeof(float));
  r->norm = get_3d_farray(n,nc,tn);
  r->opp  = get_3d_farray(n,nc*2,tn); // *2 for ON and OFF
  r->bde  = get_3d_farray(n,7,0);     // 0-pointers in 3rd dimension

  r->mt   = (float  *)myalloc(tn*sizeof(float));

  //
  // Storage for V1 adaptation - ONLY IF adaptation is turned on.
  //
  if (mod_me_v5->v1_adflag > 0){
    r->adrsp  = get_3d_farray(n,nc,tn);
    r->adst   = get_zero_2d_farray(n,nc);  // State should be initialized
  }else{
    r->adrsp  = NULL;
    r->adst   = NULL;
  }
  
  return r;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_ME_STF_W_GET_QOC2016                        */
/*                                                                           */
/*  Compute MT weights as a function of SF and TF.                           */
/*                                                                           */
/*  The equations for the preferred and opposite planes are:                 */
/*                                                                           */
/*     TF = -vx SFx  -  vy SFy      (Preferred)                              */
/*     TF =  vx SFx  +  vy SFy      (Opposite)                               */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/
float mod_me_stf_w_get_qoc2016(vx,vy,sfx,sfy,tf,sig_p,sig_a,sig_stat,statk,
			       printflag)
     float vx,vy;        // Velocity for preferred plane
     float sfx,sfy,tf;   // Coordinates of point in SF-TF space
     float sig_p;        // SD for preferred weights
     float sig_a;        // SD for antipreferred weights
     float sig_stat;     // SD for reduction of low TF (static) weights
     float statk;        // Constant for low TF suppression of weights
     int printflag;  // For debugging
{
  float pp_a,pp_b,pp_c,pp_d,dist_p,dist_a,w_ex,w_in,wp_ex,wp_in,statv,w;
  struct mod_me_v1_stf *stf;

  stf = mod_me_v5->v1stf;

  // Parameters to define the prefered plane:  aX + bY + cZ + d = 0
  pp_a = vx;
  pp_b = vy;
  pp_c = 1.0;
  pp_d = 0.0;

  if (printflag == 1){
    printf("sfx  %f\n",sfx);
    printf("sfy  %f\n",sfy);
    printf("tf   %f\n",tf);
    printf("vx,y   %f %f\n",vx,vy);
  }

  // WYETH SIGN FLIP - I am thinking of changing the sign here, e.g., the
  //   sign of 'tf', so that positive TF lives on the positive VX plane.
  //   But does this also work for + and -VY? 
  //
  dist_p = distance_point_plane(sfx,sfy,tf, pp_a, pp_b,pp_c,pp_d);
  dist_a = distance_point_plane(sfx,sfy,tf,-pp_a,-pp_b,pp_c,pp_d);

  if (printflag == 1)
    printf("dist_p %f  dist_a %f\n",dist_p,dist_a);

  w_ex = exp(-dist_p * dist_p / (2.0*sig_p*sig_p));
  w_in = exp(-dist_a * dist_a / (2.0*sig_a*sig_a));

  statv = exp(-tf*tf / (2.0*sig_stat*sig_stat));

  wp_ex = w_ex * (1.0 - statk * statv);
  wp_in = w_in * (1.0 -         statv);

  w = wp_ex - wp_in;

  //printf("*******  w_ex = %8.4f    w_in = %8.4f\n",w_ex,w_in);
  //printf("******* wp_ex = %8.4f   wp_in = %8.4f\n",wp_ex,wp_in);
  //printf("*******   w = %f\n",w);

  return w;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_ME_STF_W_GET_RADCONC                        */
/*                                                                           */
/*  Compute MT weights as a function of SF and TF.                           */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/
float mod_me_stf_w_get_radconc(vx,vy,sfx,sfy,tf,rc_rad,rc_cnc,rc_cof,rc_rsd,
			       rc_tf0,printflag)
     float vx,vy;        // Velocity for preferred plane
     float sfx,sfy,tf;   // Coordinates of point in SF-TF space
     char *rc_rad;       // Radial configuration
     char *rc_cnc;       // Concentric configuration
     float rc_cof;       // Concentric offset
     float rc_rsd;       // Radial SD
     float rc_tf0;       // Value for TF 0
     int printflag;      // For debugging
{
  float v1,v2,v3,r,theta_deg,w;
  struct mod_me_v1_stf *stf;

  stf = mod_me_v5->v1stf;

  // Parameters to define the prefered plane:  aX + bY + cZ + d = 0


  //  TF = -vx SFx - vy SFy     Preferred velocity plane
  //
  //  z = ax + by   Given this plane
  //  by = -ax      is the line where it intersects z=0
  //    (b,-a)      is a vector on that line
  //    (a, b)      is a vector orthogonal to that line
  
  // Vector pointing to maximum 
  v1 = vx; // -vx;
  v2 = vy; // -vy;
  v3 = vx*vx + vy*vy;

//tf *= -1.0;

  //theta_deg = angle_3d_vectors(sfx,sfy,tf,v1,v2,v3);
  theta_deg = angle_3d_vectors(sfx,sfy,0.0,v1,v2,0.0);
  r = sqrt(sfx*sfx + sfy*sfy + tf*tf);

  if (strcmp(rc_cnc,"cosine")==0){
    w = cos(M_PI * theta_deg / 180.0);
  }

  if (rc_tf0 != 0.0){
    if ((tf < 0.01) && (tf > -0.01))
      w += rc_tf0;
  }

  if (printflag == 1){
    printf("sfx,y tf  %6.2f %6.2f %6.2f    theta %5.1f  r %5.2f   w %f\n",sfx,sfy,tf,
	   theta_deg,r,w);
    printf("   v1,2,3   %.1f %.1f %.1f\n",v1,v2,v3);
    //printf("vx,y   %f %f\n",vx,vy);
    //printf("theta_deg = %f   r = %f\n",theta_deg,r);
    //printf("  w = %f\n",w);
  }

  //exit_error("WYETH HERE","Under development");

  //printf("*******  w_ex = %8.4f    w_in = %8.4f\n",w_ex,w_in);
  //printf("******* wp_ex = %8.4f   wp_in = %8.4f\n",wp_ex,wp_in);
  //printf("*******   w = %f\n",w);

  return w;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_ME_STF_GET_PAR_CI                          */
/*                                                                           */
/*  Return a pointer to the string value of parameter 'pname' for channel    */
/*  'ci'.                                                                    */
/*                                                                           */
/*****************************************************************************/
char *mod_me_stf_get_par_ci(pname,ci)
     char *pname;   // Parameter name
     int   ci;      // Channel index
{
  int k;
  int vi;
  struct mod_me_v1_stf *stf;

  stf = mod_me_v5->v1stf;

  // Get index for 'pname' in name list
  k = search_2d_carray(stf->par_name,pname,stf->npar);
  // WYETH vlink ? k = search_2d_carray(stf->par_name,pname,stf->nptot);

  if (k < 0)
    exit_error("MOD_ME_STF_GET_PAR_CI","Unknown STF parameter name");

  vi = stf->vndx[ci][k];

  return stf->par_val[k][vi];
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_ME_STF_GET_PARLINK_CI                         */
/*                                                                           */
/*  Return a pointer to the string value of parameter 'pname' for channel    */
/*  'ci', where 'piname' is the index parameter to which 'pname' is linked.  */
/*                                                                           */
/*****************************************************************************/
char *mod_me_stf_get_parlink_ci(pname,piname,ci)
     char *pname;    // Linked parameter name
     char *piname;   // Index parameter name
     int   ci;       // Channel index
{
  int k;
  int vi;
  struct mod_me_v1_stf *stf;

  stf = mod_me_v5->v1stf;

  //
  //  Find the index for the *index* parameter (e.g., "veli").
  //
  k = search_2d_carray(stf->par_name,piname,stf->npar);

  if (k < 0)
    exit_error("MOD_ME_STF_GET_PARLINK_CI","Unknown STF parameter name");

  vi = stf->vndx[ci][k];  // This is the index in the value list

  //
  //  Now find the index for the linked parameter (e.g., "vel_x"), whose
  //    value we need to return
  //
  k = search_2d_carray(stf->par_name,pname,stf->nptot);

  if (k < 0)
    exit_error("MOD_ME_STF_GET_PARLINK_CI","Unknown STF linked par. name");

  if (k < stf->npar)
    exit_error("MOD_ME_STF_GET_PARLINK_CI","Index value too low.");

  return stf->par_val[k][vi];
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_ME_STF_GET_TF_STF_DIR                         */
/*                                                                           */
/*  Determine the 'TF' value associated with STF channel 'stfi' and          */
/*  direction channel 'diri'.                                                */
/*                                                                           */
/*****************************************************************************/
float mod_me_stf_get_tf_stf_dir(stfi,diri,flag)
     int stfi;    // STF channel index
     int diri;    // Dir channel index, or -1 to use solely 'stfi'
     int flag;    // 1-negate TF for vlink computation
{
  int ndir;
  float dir,tf,sf,sfx,sfy,velx,vely;

  ndir = mod_me_v5->n;

  if (diri <= -1){
    //
    //  Return the TF for the 'stfi' channel.
    //
    tf = atof(mod_me_stf_get_par_ci("tf",stfi));
  }else{

    // WYETH - it seems like we should check some config value to know that
    //   we are using the velocity-linked TF here, which is the case that
    //   is in this "else".  Is there a global config we can check?

    sf  = atof(mod_me_stf_get_par_ci("sf",stfi));
    dir = (float)diri / (float)ndir * 360.0;
    sfx = sf * cos(dir * M_PI/180.0);
    sfy = sf * sin(dir * M_PI/180.0);

    // WYETH SIGN FLIP - I am thinking of flipping this, so that TF is
    //   positive for positive velx

    velx = atof(mod_me_stf_get_parlink_ci("vel_x","veli",stfi));
    vely = atof(mod_me_stf_get_parlink_ci("vel_y","veli",stfi));
    tf = -velx * sfx - vely * sfy;  // Make neg TF pos?
    if (flag == 1){
      tf *= -1.0;  // WYETH - keeps +TF for dir=0 going right ?
    }

    // WYETH HERE FIX
    if (tf > 0.0) // Anything on the far side of the plane should have its
      tf *= -1.0; //   sign inverted because it represents OPPONENT motion.
  }

  return tf;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_ME_MAKE_FILTER_LIST                         */
/*                                                                           */
/*  Make a complete list of the linear filters.                              */
/*                                                                           */
/*****************************************************************************/
void mod_me_make_filter_list(m)
     struct model_struct *m;  // Model params
{
  int i,j,k,l;
  int pi,n,ndir,nstf,nph,disp_adj,neye,re_flag,nsd,tfi,vi;
  float tdir,disp_off,phd,ph,tf,tf_off,ssdf,tsdf;
  char **list_ph,**list_tsd,*ssd,*tsd,*str_sf,*str_tf,*parfile,*ftype;
  char sdir[SLEN],eyeflag[SLEN],tph[SLEN],ttf[SLEN];
  char str_ssd[SLEN],str_tsd[SLEN];
  struct onode *v1o,*fo;
  struct param_pair_list **ppa;
  struct mod_me_v1_stf *stf;

  mylog(mylogf,"  MOD_ME_MAKE_FILTER_LIST\n");

  stf = mod_me_v5->v1stf;
  nstf = stf->n;

  v1o = onode_child_get_unique(m->o,"v1");
  fo  = onode_child_get_unique(v1o,"filters");

  ftype   = onode_getpar_chr_exit(fo,"type");
  //ndir    = onode_getpar_int_exit(fo,"n_dir");
  ndir    = mod_me_v5->n;
  if (ndir % 2 != 0)
    mylog_exit(mylogf,"  MOD_ME_MAKE_FILTER_LIST  'n_dir' must be even.");


  if (onode_item(fo,"s_sd_f")){
    ssdf = onode_getpar_flt_exit(fo,"s_sd_f");
    ssd  = NULL;  // Will be set to 's_sd_f'/SF
  }else{
    ssdf = 0.0;  // flag value
    ssd  = onode_getpar_chr_ptr(fo,"s_sd");
  }

  if (onode_item(fo,"tf_list_sd")){
    list_tsd = onode_getpar_chr_list(fo,"tf_list_sd",&nsd);

    tfi = search_2d_carray(stf->par_name,"tf",stf->npar); // index for "tf"
    if (nsd != stf->nval[tfi]){
      exit_error("MOD_ME_MAKE_FILTER_LIST",
		 "Count mismatch, 'tf_list' and 'tf_list_sd'");
    }
  }else
    list_tsd = NULL;

  if (onode_item(fo,"t_sd_f")){
    tsdf = onode_getpar_flt_exit(fo,"t_sd_f");
    tsd  = NULL;  // Will be set to 't_sd_f'/TF
  }else{
    tsdf = 0.0;  // flag value
    tsd  = onode_getpar_chr_ptr(fo,"t_sd");
  }

  if (onode_item(fo,"ph_list") == 1){  // If R.E. weights are given
    list_ph = onode_getpar_chr_list(fo,"ph_list",&nph);
  }else{
    nph = 2;
    list_ph = (char **)myalloc(nph*sizeof(char *));
    list_ph[0] = strdup("0.0");
    list_ph[1] = strdup("90.0");
  }
  parfile = onode_getpar_chr_dflt(fo,"write_par_table",NULL);

  re_flag = 0;  // Assume same filters for R.E. unless changed below

  //
  //  Disparity OFFSET in R.E.
  //
  if (onode_item(v1o,"disparity_offset") == 1){
    disp_off = onode_getpar_flt_dflt(v1o,"disparity_offset",0.0); // (deg) RE
    disp_adj = onode_getpar_int_dflt(v1o,"disparity_offset_adj",1);
    re_flag = 1;  // Global flag to indicate separate RE filters
    sprintf(ggstr,"    R.E. filter phase offset = %f\n",disp_off);
    mylog(mylogf,ggstr);
    neye = 2;
  }else{
    disp_off = 0.0;  // (deg)
    disp_adj = 0;    // (deg)
    neye = 1;
  }

  tf_off = onode_getpar_flt_dflt(v1o,"tf_offset",0.0);  // For right eye
  if (tf_off != 0.0){
    re_flag = 1;  // Global flag to indicate separate RE filters
    neye = 2;
  }



  n = neye * nstf * ndir * nph;
  ppa = (struct param_pair_list **)myalloc(n*sizeof(struct param_pair_list *));


  //
  //  Create top level filter structure, and store param lists there.
  //
  mod_me_fs_top_create(n);
  mod_me_ftop->fpl = ppa;
  mod_me_ftop->flag_re = re_flag; // ******** WYETH HERE NOW


  //  Create the filter index.
  mod_me_ftop->nchan = nstf;
  mod_me_ftop->ndir  = ndir;
  mod_me_ftop->nph   = nph;
  mod_me_ftop->fndx  = get_3d_iarray(nstf,ndir,nph);


  pi = 0;
  for(i=0;i<neye;i++){

    if (re_flag == 1){
      sprintf(eyeflag,"%d",i);     // 0-Left, 1-Right
    }else{
      strcpy(eyeflag,"2");         // 2-Both
    }

    for(j=0;j<nstf;j++){       // For each STF channel
      str_sf = mod_me_stf_get_par_ci("sf",j);
      str_tf = mod_me_stf_get_par_ci("tf",j);

      //
      //  Determine the TF and the TF SD for this channel
      //
      sprintf(ttf,"%s",str_tf);  // Default value, could be changed below

      if (strcmp(eyeflag,"1")==0){  // Filter for R.E. alone
	//
	//  TF offset in R.E.
	//
	if (tf_off != 0.0){  // There is a non-zero factor for computing t_SD
	  tf = atof(str_tf) + tf_off;
	  sprintf(ttf,"%f",tf);  // WYETH - how many decimals do we want?
	}
      }
      if (list_tsd != NULL){  // Use the SD value in the 'tf_list_sd'
	// We must compute the index in 'tf_list' associated with this 'j'
	//   STF channel, and then use it for the 'list_tsd'.
	vi = stf->vndx[j][tfi];
	tsd = list_tsd[vi];
      }else if (tsdf != 0.0){
	sprintf(str_tsd,"%.5f",tsdf / atof(ttf));
	tsd = str_tsd;
      }

      for(k=0;k<ndir;k++){     // For each Direction

	stf->flag_dir[j][k] = 1;  // flag that there is a filter for this dir.

	tdir = (float)k/(float)ndir * 360.0;
	sprintf(sdir,"%d",my_rint(tdir));

	for(l=0;l<nph;l++){    // For each Phase

	  sprintf(tph,"%s",list_ph[l]); // Default value, could change below

	  if (strcmp(eyeflag,"1")==0){  // Filter for R.E. alone
	    //
	    //  Disparity OFFSET for R.E.  --- Correct for direction
	    //
	    if (disp_off != 0.0){
	      phd = disparity_phase_shift_adjust(disp_adj,tdir,disp_off);
	      ph = atof(list_ph[l]) + phd;
	      sprintf(tph,"%f",ph);  // WYETH - how many decimals do we want?
	    }
	  }

	  //
	  //  If needed, compute SD values, spatial and temporal
	  //
	  if (ssdf != 0.0){
	    sprintf(str_ssd,"%.5f",ssdf / atof(str_sf));
	    ssd = str_ssd;
	  }

	  if (parfile != NULL){
	    sprintf(ggstr,"%3d %s %5s %5s %5s %6s %s %s %s\n",pi,eyeflag,
		    ttf,str_sf,sdir,tph,ssd,tsd,ftype);
	    append_string_to_file(parfile,ggstr);
	  }
	  
	  ppa[pi] = paramfile_ppl_get_init();
	  paramfile_ppl_add_name_val(ppa[pi],"tf",ttf,1);
	  paramfile_ppl_add_name_val(ppa[pi],"sf",str_sf,1);
	  paramfile_ppl_add_name_val(ppa[pi],"dir",sdir,1);
	  paramfile_ppl_add_name_val(ppa[pi],"ph",tph,1);

	  paramfile_ppl_add_name_val(ppa[pi],"s_sd",ssd,1);
	  paramfile_ppl_add_name_val(ppa[pi],"t_sd",tsd,1);
	  paramfile_ppl_add_name_val(ppa[pi],"type",ftype,1);

	  paramfile_ppl_add_name_val(ppa[pi],"eyeflag",eyeflag,1);

	  if (i==0) // Set index for L.E. only
	    mod_me_ftop->fndx[j][k][l] = pi;

	  pi += 1;
	}
      }
    }
  }

  myfree(ftype);
  if (parfile != NULL)  myfree(parfile);
  free_2d_carray(list_ph,nph);
}
/**************************************-**************************************/
/*                                                                           */
/*                      MOD_ME_MAKE_FILTER_LIST_LINKTFVEL                    */
/*                                                                           */
/*  The STF params are 'vel' and 'SF' and the TF needs to be computed from   */
/*  those so that all filters lie on the velocity plane(s).                  */
/*                                                                           */
/*****************************************************************************/
void mod_me_make_filter_list_linktfvel(m)
     struct model_struct *m;  // Model params
{
  int i,j,k,l;
  int pi,n,ndir,nstf,nph,disp_adj,neye,re_flag,nsd,tfi;
  int nvx,nvy,zflag;
  float tdir,disp_off,phd,ph,tf,tf_off,ssdf,tsdf;
  char **list_ph,**list_tsd,*ssd,*tsd,*str_sf,*str_tf,*parfile,*ftype;
  char **list_velx,**list_vely;
  char sdir[SLEN],eyeflag[SLEN],tph[SLEN],ttf[SLEN],str_ssd[SLEN];
  char str_tsd[SLEN];
  struct onode *v1o,*fo;
  struct param_pair_list **ppa;
  struct mod_me_v1_stf *stf;

  mylog(mylogf,"  MOD_ME_MAKE_FILTER_LIST_LINKTFVEL\n");


  stf = mod_me_v5->v1stf;
  nstf = stf->n;

  v1o = onode_child_get_unique(m->o,"v1");
  fo  = onode_child_get_unique(v1o,"filters");

  ftype   = onode_getpar_chr_exit(fo,"type");
  //ndir    = onode_getpar_int_exit(fo,"n_dir");
  ndir    = mod_me_v5->n; // onode_getpar_int_exit(fo,"n_dir");
  if (ndir % 2 != 0)
    mylog_exit(mylogf,
	       "  MOD_ME_MAKE_FILTER_LIST_LINKTFVEL  'n_dir' must be even.");


  if (onode_item(fo,"s_sd_f")){
    ssdf = onode_getpar_flt_exit(fo,"s_sd_f");
    ssd  = NULL;  // Will be set to 's_sd_f'/SF
  }else{
    ssdf = 0.0;  // flag value
    ssd  = onode_getpar_chr_ptr(fo,"s_sd");
  }

  /****** WYETH NO TF_LIST_SD allowed here
  if (onode_item(fo,"tf_list_sd")){
    list_tsd = onode_getpar_chr_list(fo,"tf_list_sd",&nsd);

    tfi = search_2d_carray(stf->par_name,"tf",stf->npar); // index for "tf"
    if (nsd != stf->nval[tfi]){
      exit_error("MOD_ME_MAKE_FILTER_LIST",
		 "Count mismatch, 'tf_list' and 'tf_list_sd'");
    }
  }else
    list_tsd = NULL;
  ***/

  if (onode_item(fo,"t_sd_f")){
    tsdf = onode_getpar_flt_exit(fo,"t_sd_f");
    tsd  = NULL;  // Will be set to 't_sd_f'/TF
  }else{
    tsdf = 0.0;  // flag value
    tsd  = onode_getpar_chr_ptr(fo,"t_sd");
  }

  if (onode_item(fo,"ph_list") == 1){  // If R.E. weights are given
    list_ph = onode_getpar_chr_list(fo,"ph_list",&nph);
  }else{
    nph = 2;
    list_ph = (char **)myalloc(nph*sizeof(char *));
    list_ph[0] = strdup("0.0");
    list_ph[1] = strdup("90.0");
  }
  parfile = onode_getpar_chr_dflt(fo,"write_par_table",NULL);

  re_flag = 0;  // Assume same filters for R.E. unless changed below

  //
  //  Disparity OFFSET in R.E.
  //
  if (onode_item(v1o,"disparity_offset") == 1){
    disp_off = onode_getpar_flt_dflt(v1o,"disparity_offset",0.0); // (deg) RE
    disp_adj = onode_getpar_int_dflt(v1o,"disparity_offset_adj",1);
    re_flag = 1;  // Global flag to indicate separate RE filters
    sprintf(ggstr,"    R.E. filter phase offset = %f\n",disp_off);
    mylog(mylogf,ggstr);
    neye = 2;
  }else{
    disp_off = 0.0;  // (deg)
    disp_adj = 0;    // (deg)
    neye = 1;
  }

  //  WYETH HERE - should we allow TF offset in this model?
  //  WYETH HERE - should we allow TF offset in this model?
  //  WYETH HERE - should we allow TF offset in this model?
  tf_off = onode_getpar_flt_dflt(v1o,"tf_offset",0.0);  // For right eye
  if (tf_off != 0.0){
    exit_error("MOD_ME_MAKE_FILTER_LIST_LINKTFVEL",
	       "'tf_offset' not implemented here");
    //re_flag = 1;  // Global flag to indicate separate RE filters
    //neye = 2;
  }

  //
  //  Compute total number of filters to build.
  //
  n = neye * nstf * ndir * nph;  // may be adjusted below for shared filters

  list_velx = onode_getpar_chr_list(fo,"vel_x",&nvx);
  list_vely = onode_getpar_chr_list(fo,"vel_y",&nvy);
  if (nvx != nvy){
    mylog_exit(mylogf,
      "  MOD_ME_MAKE_FILTER_LIST_LINKTFVEL  'vel_x', '_y' counts differ.");
  }

  ppa = (struct param_pair_list **)myalloc(n*sizeof(struct param_pair_list *));

  //
  //  Create top level filter structure, and store param lists there.
  //
  mod_me_fs_top_create(n);
  mod_me_ftop->fpl = ppa;
  mod_me_ftop->flag_re = re_flag; // ******** WYETH HERE NOW


  //  Create the filter index.
  mod_me_ftop->nchan = nstf;
  mod_me_ftop->ndir  = ndir;
  mod_me_ftop->nph   = nph;
  mod_me_ftop->fndx  = get_3d_iarray(nstf,ndir,nph);

  pi = 0;
  for(i=0;i<neye;i++){

    if (re_flag == 1){
      sprintf(eyeflag,"%d",i);     // 0-Left, 1-Right
    }else{
      strcpy(eyeflag,"2");         // 2-Both
    }

    // WYETH vlink - this becomes for each "SF"
    for(j=0;j<nstf;j++){       // For each STF channel

      str_sf = mod_me_stf_get_par_ci("sf",j);
  
      zflag = 0;  // TF 0 has not yet been built on this SF-V plane

      for(k=0;k<ndir;k++){     // For each Direction

	stf->flag_dir[j][k] = 1;  // Could be changed to '2' below

	tdir = (float)k/(float)ndir * 360.0;
	sprintf(sdir,"%d",my_rint(tdir));

	tf = mod_me_stf_get_tf_stf_dir(j,k,1);
	if (tf < 0.0)
	  tf *= -1.0;  // Use positive TF
	sprintf(ttf,"%f",tf);  // WYETH - how many decimals do we want?

	if (tsdf != 0.0){
	  sprintf(str_tsd,"%.5f",tsdf / tf);
	  tsd = str_tsd;  // In this case, 'tsd' was not set above
	}

	if (fabs(tf) < 0.01){  // The TF is near zero
	  zflag += 1;

	  if (zflag > 1){
	    stf->flag_dir[j][k] = 2;  // This filter is redundant
   printf("    Redundant filter, j %d   zflag = %d\n",j,zflag);
	  }
	}

	for(l=0;l<nph;l++){    // For each Phase

	  sprintf(tph,"%s",list_ph[l]); // Default, could change below

	  if (strcmp(eyeflag,"1")==0){  // Filter for R.E. alone
	    //
	    //  Disparity OFFSET for R.E.  --- Correct for direction
	    //
	    if (disp_off != 0.0){
	      phd = disparity_phase_shift_adjust(disp_adj,tdir,disp_off);
	      ph = atof(list_ph[l]) + phd;
	      sprintf(tph,"%f",ph);  // WYETH - how many decimals here?
	    }
	  }

	  //
	  //  If needed, compute SD values, spatial and temporal
	  //
	  if (ssdf != 0.0){
	    sprintf(str_ssd,"%.5f",ssdf / atof(str_sf));
	    ssd = str_ssd;
	  }

	  if (parfile != NULL){
	    sprintf(ggstr,"%3d %s %5s %5s %5s %6s %s %s %s\n",pi,eyeflag,
		    ttf,str_sf,sdir,tph,ssd,tsd,ftype);
	    append_string_to_file(parfile,ggstr);
	  }

	  ppa[pi] = paramfile_ppl_get_init();
	  paramfile_ppl_add_name_val(ppa[pi],"tf",ttf,1);
	  paramfile_ppl_add_name_val(ppa[pi],"sf",str_sf,1);
	  paramfile_ppl_add_name_val(ppa[pi],"dir",sdir,1);
	  paramfile_ppl_add_name_val(ppa[pi],"ph",tph,1);

	  paramfile_ppl_add_name_val(ppa[pi],"s_sd",ssd,1);
	  paramfile_ppl_add_name_val(ppa[pi],"t_sd",tsd,1);
	  paramfile_ppl_add_name_val(ppa[pi],"type",ftype,1);

	  paramfile_ppl_add_name_val(ppa[pi],"eyeflag",eyeflag,1);

	  if (i==0) // Set index for L.E. only
	    mod_me_ftop->fndx[j][k][l] = pi;

	  pi += 1;
	}
      }
    }
  }

  if (pi != n)
    mylog_exit(mylogf,"  MOD_ME_MAKE_FILTER_LIST_LINKTFVEL  Bad 'pi'.");


  myfree(ftype);
  if (parfile != NULL)  myfree(parfile);
  free_2d_carray(list_ph,nph);
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_ME_STF_CONFIG                            */
/*                                                                           */
/*  Configure a bank of spatio-temporal channels.                            */
/*  Developed to extend the ME_V5 model.                                     */
/*                                                                           */
/*****************************************************************************/
void mod_me_stf_config(m)
     struct model_struct *m; // Model params
{
  int i,j;
  int ndir;
  char tstr[SLEN];
  struct onode *v1o,*fo;
  struct mod_me_v1_stf *stf;

  mylog(mylogf,"  MOD_ME_STF_CONFIG\n");

  v1o  = onode_child_get_unique(m->o,"v1");
  fo   = onode_child_get_unique(v1o,"filters");
  ndir = mod_me_v5->n;  // WYETH OLD: onode_getpar_int_exit(fo,"n_dir");

  stf = (struct mod_me_v1_stf *)myalloc(sizeof(struct mod_me_v1_stf));

  stf->config = onode_getpar_chr_exit(fo,"config");
  if (strcmp(stf->config,"SFxTF")==0){

    stf->npar = 2;
    stf->nptot = stf->npar;  // There are only index parameters
    stf->flag_config = 0;    // TF is directly available

    stf->par_name = (char **)myalloc(stf->nptot*sizeof(char *));
    stf->par_name[0] = strdup("tf");
    stf->par_name[1] = strdup("sf");

    stf->nval = (int *)myalloc(stf->nptot*sizeof(int));

    stf->par_val = (char ***)myalloc(stf->nptot*sizeof(char **));
    stf->par_val[0] = onode_getpar_chr_list(fo,"tf_list",&(stf->nval[0]));
    stf->par_val[1] = onode_getpar_chr_list(fo,"sf_list",&(stf->nval[1]));

    // Use only the *index* pars ('npar' of them) for this
    get_nd_cross_product_list_iarray(stf->npar,stf->nval,
				     &(stf->vndx),&(stf->n));
    //  Example of order for TF: 0,1 and SF 0,1,2:
    //    0 0
    //    1 0
    //    0 1
    //    1 1
    //    0 2
    //    1 2   Thus, TF changes more rapidly than SF

    /*
      for(i=0;i<stf->n;i++){
      for(j=0;j<stf->npar;j++){
      printf("  %d",stf->vndx[i][j]);
      }
      printf("\n");
      }
    */
  }else if (strcmp(stf->config,"link_TF_vel")==0){

    stf->npar = 2;
    stf->nptot = stf->npar + 2;  // For 'vel_x' and 'vel_y'
    stf->flag_config = 1;        // TF must be calculated

    stf->nval     =    (int *)myalloc(stf->nptot*sizeof(int));
    stf->par_val  = (char ***)myalloc(stf->nptot*sizeof(char **));
    stf->par_name =  (char **)myalloc(stf->nptot*sizeof(char *));
    stf->par_name[0] = strdup("veli");
    stf->par_name[1] = strdup("sf");
    stf->par_name[2] = strdup("vel_x");
    stf->par_name[3] = strdup("vel_y");

    //
    //  Set velocity values  - Use index values of 0, 1, ... for 'veli'
    //
    stf->par_val[0] = onode_getpar_chr_list(fo,"vel_x",&(stf->nval[0]));
    stf->par_val[2] = onode_getpar_chr_list(fo,"vel_x",&(stf->nval[2]));
    stf->par_val[3] = onode_getpar_chr_list(fo,"vel_y",&(stf->nval[3]));
    for(i=0;i<stf->nval[0];i++){
      myfree(stf->par_val[0][i]);  // Replace the strings with index values
      sprintf(tstr,"%d",i);
      stf->par_val[0][i] = strdup(tstr);
    }

    //
    // Set SF values (given in the param list in the .moo file)
    //
    stf->par_val[1] = onode_getpar_chr_list(fo,"sf_list",&(stf->nval[1]));

    // Use only the index parameters for this
    get_nd_cross_product_list_iarray(stf->npar,stf->nval,
				     &(stf->vndx),&(stf->n));

    // The number of STF channels is number of velocities times number of SFs
    stf->n = stf->nval[0] * stf->nval[1];

  }else
    exit_error("MOD_ME_STF_CONFIG","Unknown 'config' in <filters>");

  stf->flag_dir = get_zero_2d_iarray(stf->n,ndir);

  stf->sgr = (struct mod_me_v1_spatgrp **)myalloc(stf->n*
					  sizeof(struct mod_me_v1_spatgrp *));

  stf->sum_norm   = get_3d_farray(mod_me_xn,mod_me_yn,mod_me_tn);
  stf->sum_norm_r = get_3d_farray(mod_me_xn,mod_me_yn,mod_me_tn);

  mod_me_v5->v1stf = stf;

  //
  //  Create an array of 'param_pair_list' structs, one for each filter
  //  Create top level filter structure, and store param lists there.
  //  Create the filter index table.
  //
  if (strcmp(stf->config,"SFxTF")==0){
    mod_me_make_filter_list(m);
  }else if (strcmp(stf->config,"link_TF_vel")==0){
    mod_me_make_filter_list_linktfvel(m);
  }else
    exit_error("MOD_ME_STF_CONFIG","Unknown 'config'");
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_ME_BDGR_2D_INIT                           */
/*                                                                           */
/*  Allocate storage for V1 responses at each spatial location.              */
/*                                                                           */
/*****************************************************************************/
struct mod_me_v1_input ***mod_me_bdgr_2d_init(txn,tyn,tzn,cmap)
     int txn,tyn;       // Size of spatial array
     int tzn;           // Depth of 'cmap'
     int ***cmap;       // [txn][tyn][tzn]
{
  int i,j,k;
  int tn;
  struct mod_me_v5_struct *mv5;
  struct mod_me_v1_input ***r;  // [xn][yn] or [1][1] Input signals

  mv5 = mod_me_v5;
  tn = mod_me_tn;

  r = (struct mod_me_v1_input ***)myalloc(txn*
					  sizeof(struct mod_me_v1_input **));
  for(i=0;i<txn;i++){
    r[i] = (struct mod_me_v1_input **)myalloc(tyn*
					     sizeof(struct mod_me_v1_input *));
    for(j=0;j<tyn;j++){
      r[i][j] = (struct mod_me_v1_input *)myalloc(
					       sizeof(struct mod_me_v1_input));
      r[i][j]->flag = 0;  // For now, no storage
    }
  }

  for(i=0;i<txn;i++){
    for(j=0;j<tyn;j++){
      for(k=0;k<tzn;k++){
	//
	//  If any filter/channel at this position is needed, make storage
	//  *** WYETH this is inefficient, but a quick fix for now Jun2016
	//
	if ((cmap[i][j][k] > 0) && (r[i][j]->flag == 0)){
	  r[i][j]->flag = 2;
	  r[i][j]->dl = model_me_v5_resp_init(mv5->n,mv5->v1_nc,tn);
	  // WYETH varmoo_bug; I think we make this RE storage whether we
	  //   will need it or not.
	  r[i][j]->dr = model_me_v5_resp_init(mv5->n,mv5->v1_nc,tn);
	  r[i][j]->mtb  = (float  *)myalloc(tn*sizeof(float));
	}
      }
      if (r[i][j]->flag == 0){  // Nothing to save at this [i][j]
	r[i][j]->dl  = NULL;
	r[i][j]->dr  = NULL;
	r[i][j]->mtb = NULL;
      }
    }
  }
  return r;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_ME_V5_POP_CONN                            */
/*                                                                           */
/*  Create inputs to MT units.                                               */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_pop_conn(mtpop,pl,io,disto)
     struct mod_me_v5_pop *mtpop;  // MT pop structure
     struct pop_layer *pl;         // MT layer to configure
     struct onode *io;             // 'input' onode
     struct onode *disto;          // 'distrib' onode
{
  int i,j,k;
  int sti,mxn,myn,mzn,mx0,my0,xn,yn,nposs,nn,xi,yi,zi,ci,even_out,pflag,seed;
  int maxrep,evflag,stype,nseed,nsamp,nstf,*cx,*cy,cn,*seedlist,cnt,orii,flagd;
  int hvflag,dumpflag,evit;
  short inindex;
  float cdist,minw,normw,sdo,sdp,xc,yc,ori,tdelay,x0,y0,velx,vely,eopar;
  float t_sf,t_tf,dirv1,t_sfx,t_sfy,ww,vx,vy,dsdp,dsda,ltfsd,ltfk;
  char *shape,*dconf,*dirmap,*wconf,*wfile,tstr[SLEN];
  struct mod_me_v5_struct *mv5;
  struct pop_cell *pre,*post;
  struct pop_syn *tsyn;

  // WYETH HERE - implementing new weight distrib., alternative to QOC
  char *rc_rad,*rc_cnc;
  float rc_cof,rc_rsd,rc_tf0;


  mv5 = mod_me_v5;
  xn  = mod_me_xn;
  yn  = mod_me_yn;

  eopar = -1.0;  // Dummy value for optional parameter.  May be set below.

  if (mv5->v1stf != NULL){  // New way
    nstf = mv5->v1stf->n;
  }else{
    //exit_error("here wyeth","stopping");
    nstf = 1;
  }

  mxn = pl->xn;
  myn = pl->yn;
  mzn = pl->zn;

  if ((mxn > xn) || (myn > yn))
    mylog_exit(mylogf,"MOD_ME_V5_POP_CONN  MT popolation is too large");

  shape    = onode_getpar_chr_dflt(disto,"shape","Gaussian");
  cdist    = onode_getpar_flt_dflt(disto,"cdist",0.0);
  minw     = onode_getpar_flt_dflt(disto,"minw",0.0);
  normw    = onode_getpar_flt_dflt(disto,"normw",0.0);
  evflag   = onode_getpar_int_dflt(disto,"even_out",0);
  if (evflag == 1)
    eopar = onode_getpar_flt_dflt(disto,"even_out_sig",2.0);
  evit     = onode_getpar_int_dflt(disto,"even_out_itmax",-1);
  dumpflag = onode_getpar_int_dflt(disto,"dumpflag",0);
  seed     = onode_getpar_int_dflt(disto,"seed",1777);
  nsamp    = onode_getpar_int_exit(disto,"nsamp");
  dconf    = onode_getpar_chr_dflt(io,"dir_config","stacked");
  pflag    = onode_getpar_int_dflt(io,"print_connections",0);
  dirmap   = onode_getpar_chr_dflt(io,"dir_map","none");
  hvflag   = onode_getpar_int_dflt(io,"horiz_flag",1);

  mtpop->nsamp = nsamp;

  sprintf(ggstr,"    %s  cdist %f  minw %f  normw %f  seed %d\n",
	  shape,cdist,minw,normw,seed);
  mylog(mylogf,ggstr);

  sdo = sdp = cdist / mod_me_sscale;  // Make both SDs the same
  xc = yc = 0.0;
  ori = 0.0;
  maxrep = 1;
  stype = 1;     // 1- excitatory synapse
  tdelay = 0.0;  // No synaptic dealy

  if (strcmp(dconf,"stacked")==0){  // Use all dir's at each location
    nseed = mxn*myn*mzn * nstf;
    mtpop->stackflag = 1;
  }else if (strcmp(dconf,"random")==0){  // Vary dir by location randomly
    nseed = mxn*myn*mzn*mv5->n * nstf;
    mtpop->stackflag = 0;
  }else{
    mylog_exit(mylogf,"MOD_ME_V5_POP_CONN  Unknown 'dir_config' value");
  }
  seedlist = get_seeds(seed,99999,nseed);

  if (strcmp(dirmap,"z_dir_chan")==0){
    mtpop->zdirflag = 1;
    if (pl->zn != mv5->n)
      exit_error("MOD_ME_V5_POP_CONN",
		 "For 'z_dir_chan', 'zn' must equal 'n_dir_chan'");
    orii = pop_cell_attrib_index_f(pl,"ori");
    if (orii < 0)
      exit_error("MOD_ME_V5_POP_CONN","No index found for 'ori' attrib");
  }else
    mtpop->zdirflag = 0;

  wfile = NULL;
  if (mtpop->wco != NULL){

    wconf = onode_getpar_chr_ptr_exit(mtpop->wco,"type");

    if (strcmp(wconf,"QOC_2016")==0){
      //
      //  Read the parameters for the QOC weight scheme
      //
      vx    = onode_getpar_flt_exit(mtpop->wco,"velocity_x");
      vy    = onode_getpar_flt_exit(mtpop->wco,"velocity_y");
      dsdp  = onode_getpar_flt_exit(mtpop->wco,"dist_sd_pref");
      dsda  = onode_getpar_flt_exit(mtpop->wco,"dist_sd_anti");
      ltfsd = onode_getpar_flt_exit(mtpop->wco,"low_tf_sd");
      ltfk  = onode_getpar_flt_exit(mtpop->wco,"low_tf_k");
      wfile = onode_getpar_chr_dflt(mtpop->wco,"write_w_table",NULL);

      //printf(" velx %f  vely %f\n",vx,vy);
      //printf(" dsdp %f  dsda %f\n",dsdp,dsda);
      //printf(" lowtf_sd %f  lowtf_k %f\n",ltfsd,ltfk);

    }else if (strcmp(wconf,"Rad_Conc")==0){
      //
      //  Radial - Concentric
      //
      vx     = onode_getpar_flt_exit(mtpop->wco,"velocity_x");
      vy     = onode_getpar_flt_exit(mtpop->wco,"velocity_y");
      rc_rad = onode_getpar_chr_dflt(mtpop->wco,"radial","const");
      rc_cnc = onode_getpar_chr_dflt(mtpop->wco,"concentric","const");
      rc_cof = onode_getpar_flt_dflt(mtpop->wco,"conc_offset",0.0);
      rc_rsd = onode_getpar_flt_dflt(mtpop->wco,"radial_sd",0.0);
      rc_tf0 = onode_getpar_flt_dflt(mtpop->wco,"tf0_offset",0.0);
      wfile  = onode_getpar_chr_dflt(mtpop->wco,"write_w_table",NULL);

    }else{
      printf("  *** weight_config = %s\n",wconf);
      mylog_exit(mylogf,"MOD_ME_V5_POP_CONN  Unknown 'weight_config'.");
    }
  }


  // *** WYETH NOTE - Am adding area offset to layer geom offset.
  // *** What really is the convention?
  x0 = pl->area->x0 + pl->x0;
  y0 = pl->area->y0 + pl->y0;

  ww = 1.0;  // Default synaptic weight

  cnt = 0;
  for(i=0;i<mxn;i++){
    for(j=0;j<myn;j++){
      for(k=0;k<mzn;k++){  // For each MT unit

	//printf("--------------------=====> UUUUUNIT:   i,j,k = %d %d %d \n",
	//i,j,k);

	post = &(pl->c[i][j][k]);

	if (strcmp(dirmap,"z_dir_chan")==0){
	  // Set a value for the cell attrib 'ori' based on the z-index
	  post->attrib_f[orii] = (float)k/(float)mv5->n * 360.0;
	}

	xc = x0 + (float)post->layx * pl->xf; // Note, these are (float) vals
	yc = y0 + (float)post->layy * pl->yf;

	//
	// WYETH vlink - IS THIS OK THAT WE HAVE EACH "STI" in a different
	//   pop layer???  For VLINK, STI does not have a *SINGLE* TF value
	//   and thus the "t_tf" value below will have to be calculated at
	//   a different place.  Maybe this is not the best way to do this?
	//

	for(sti=0;sti<nstf;sti++){  // For each STF channel

	  // *** WYETH - UNCLEAR WHAT TO DO WITH SEEDS FOR SPATIAL RANDOMNESS
	  // ***   below.  How many do we need?  Do we change across STF chans?
	  // ***   This may require some thought.

	  if (mv5->v1stf != NULL){  // New way
	    t_sf = atof(mod_me_stf_get_par_ci("sf",sti));
	    if (mv5->v1stf->flag_config == 0)
	      // TF is directly available, so look it up for this channel
	      t_tf = atof(mod_me_stf_get_par_ci("tf",sti));
	  }

	  cx = cy = NULL;
	  for(zi=0;zi<mv5->n;zi++){  // For each direction channel

	    if (mv5->v1stf == NULL){
	      flagd = 1;
	    }else{
	      if (mv5->v1stf->flag_dir[sti][zi] == 1)
		flagd = 1;
	      else
		flagd = 0;
	    }

	    //if (mv5->v1stf->flag_dir[sti][zi] == 1){
	    if (flagd == 1){
	      //
	      //  This direction channel exists for this STF channel
	      //

	      if (mtpop->wco != NULL){

		//  Use SF, TF and dir value to compute the weight
		dirv1 = (float)zi / (float)mv5->n * 360.0;
		t_sfx = t_sf * cos(dirv1 * M_PI/180.0);
		t_sfy = t_sf * sin(dirv1 * M_PI/180.0);


		if (mv5->v1stf->flag_config == 1){
		  // Compute TF for this channel and dir.
		  t_tf = mod_me_stf_get_tf_stf_dir(sti,zi,1);  // or 1?
		}

		//printf("  Dir %f\n",dirv1);
		//printf("    sfx  %f\n",t_sfx);
		//printf("    sfy  %f\n",t_sfy);

		if (strcmp(wconf,"QOC_2016")==0){
		  ww = mod_me_stf_w_get_qoc2016(vx,vy,t_sfx,t_sfy,t_tf,dsdp,
						dsda,ltfsd,ltfk,0);
		}else if (strcmp(wconf,"Rad_Conc")==0){

		  if ((i==0)&&(j==0)&&(k==0))
		    // WYETH - Set last param, printflag 1 here for debug
		  ww = mod_me_stf_w_get_radconc(vx,vy,t_sfx,t_sfy,-t_tf,rc_rad,
						rc_cnc,rc_cof,rc_rsd,rc_tf0,0);
		  else
		  ww = mod_me_stf_w_get_radconc(vx,vy,t_sfx,t_sfy,-t_tf,rc_rad,
						rc_cnc,rc_cof,rc_rsd,rc_tf0,0);
		}

		if (mtpop->negwf != 1.0){
		  if (ww < 0.0){
		    ww *= mtpop->negwf;
		  }
		}

		if (wfile != NULL){
		  if (i==0 && j==0 && k==0){
		    sprintf(tstr,"%3d %6.2f %6.2f %3d %9.6f\n",sti,t_tf,t_sf,
			    my_rint(dirv1),ww);
		    append_string_to_file(wfile,tstr);
		    //printf(" ci %d  di %d  ww = %f\n",sti,zi,ww);
		  }
		}
	      }

	      if ((mtpop->stackflag == 0) || (cx == NULL)){
		//
		//  If we need to pick new coordinates for this 'zi'
		//
		if (cx != NULL){
		  myfree(cx);
		  myfree(cy);
		}

		if (nsamp == 2){

		  cn = 2;
		  cx = get_zero_iarray(cn);
		  cy = get_zero_iarray(cn);

		  if (hvflag == 1){  // Horizontal
		    cx[0] = xn/4+1;
		    cy[0] = yn/2;
		    cx[1] = xn/4 + xn/2 - 2;
		    cy[1] = yn/2;
		  }else{ // Vertical
		    cx[0] = xn/2;
		    cy[0] = yn/4+1;
		    cx[1] = xn/2;
		    cy[1] = yn/4 + yn/2 - 2;
		  }

		}else{
		  nposs = mod_conn_gauss_02(mylogf,xn,yn,shape,xc,yc,sdo,sdp,
					    ori,seedlist[cnt],maxrep,minw,
					    nsamp,evflag,eopar,evit,dumpflag,
					    &cx,&cy,&cn);
		  cnt += 1;
		}
	      }

	      for(ci=0;ci<cn;ci++){  // For each input
		xi = cx[ci];
		yi = cy[ci];

		// ****** WYETH HERE ****** 'zi' should probably be mult'd
		// ****** WYETH HERE ******  because there are multiple filters
		// ****** WYETH HERE ******  for each direction channel.
		// ****** WYETH HERE ******  - Well, as it is now, the MAR
		// ****** WYETH HERE ******    file does not show the right
		// ****** WYETH HERE ******    connectivity, only the first
		// ****** WYETH HERE ******    n-dir filters are used.
		// ****** WYETH HERE ******   June 2016

		pre = &(mod_me_mpt->lay[sti+1]->c[xi][yi][zi]);  // V1 STF layer
		if (pflag == 1){
		  sprintf(ggstr,
			  "    Connecting V1 (%d,%d,%d) to %s (%d,%d,%d)\n",
			  xi,yi,zi,pl->name,i,j,k);
		  mylog(mylogf,ggstr);
		}

		inindex = sti;  // WYETH - be careful to match this index
                                //   to the inlist prepared in ..._pop_prep
		tsyn = pop_cell_add_synapse(pre,post,stype,ww,tdelay,inindex);
	      }
	    }
	  }
	  if (cx != NULL){
	    myfree(cx);  // WYETH - I think this is the right place to free?
	    myfree(cy);
	    cx = cy = NULL;
	  }
	}
      }
    }
  }

  myfree(dirmap);
  myfree(seedlist);
}
/**************************************-**************************************/
/*                                                                           */
/*                            MODEL_ME_V5_POP_PREP                           */
/*                                                                           */
/*****************************************************************************/
void model_me_v5_pop_prep(m)
     struct model_struct *m; // Model params
{
  int i,j,k;
  int **fx,**fy,*fn,zn,popcnt,re_flag,nfilt,nv1,uflag;
  float **fw;
  char *pname,**aname,lname[SLEN];
  struct onode *ao,*po,*io,*disto,*v5o,*nlo,*sgen,*tio,*nrmo,*inso;
  struct pop_layer *pl,*subpl;
  struct mod_me_v5_struct *mv5;
  struct pop_top *mpt;
  struct mod_me_v5_pop *mtpop;

  mylog(mylogf,"  MODEL_ME_V5_POP_PREP\n");

  mv5 = mod_me_v5;

  //
  //  The existence of the <pop> object indicates there is at least one
  //  MT population with spatial summation.
  //
  popcnt = onode_count_otype(m->o,"pop");
  if (popcnt < 1){
    mod_me_mpt = NULL;
    mylog(mylogf,"    No MT populations.\n");
    return;
  }else{
    if (mv5->w != NULL)
      exit_error("MODEL_ME_V5_POP_PREP","Remove old format '<v5>' in .moo");
  }
  sprintf(ggstr,"    Found %d population(s).\n",popcnt);
  mylog(mylogf,ggstr);


  //
  //  Create a 'pop_top' structure
  //
  mod_me_mpt = popu_make_pop_top(mylogf,mod_me_sscale,mod_me_tscale,
				 mod_me_xn,mod_me_yn,mod_me_tn);
  mpt = mod_me_mpt;  // Local name

  mpt->name = strdup("ME_V5");

  //
  //  Areas
  //
  ao = onode_child_get_unique(m->o,"area");  // <area>
  if (ao == NULL)
    mylog_exit(mylogf,"MOD_ME_V5_POP_PREP  Must define <area> at top level.");

  mpt->narea = 2;  // 'v1' and 'mt'
  mpt->area = (struct pop_area **)myalloc(mpt->narea*
					  sizeof(struct pop_area *));
  mpt->area[0] = popu_make_area_default("v1",mpt->xn,mpt->yn,mpt->sscale,20.0);
  mpt->area[1] = popu_make_area(logf,ao,mod_me_sscale);


  //
  //  Layers
  //

  //
  //  WYETH STF
  //  Make the number of V1 layers be the number of STF channels.
  //
  if (mv5->v1stf != NULL){  // New way
    nv1 = mv5->v1stf->n;
  }else
    nv1 = 1;

  mv5->mtli0 = 1 + nv1;   // Index of 1st MT layer in 'mpt'

  mpt->nlay = mv5->mtli0 + popcnt;  // 'stim', 'v1', and all 'mt' pops.
  mpt->lay = (struct pop_layer **)myalloc(mpt->nlay*
					  sizeof(struct pop_layer *));

  //
  //  Layer[0]  -  Stimulus layer
  //
  pl = popu_make_layer_cart("stim",mpt->xn,mpt->yn,1);
  myfree(pl->laytype);
  pl->laytype = strdup("stim");
  pl->geomt = strdup("default");
  pl->x0 = pl->y0 = 0.0;
  pl->xf = pl->yf = 1.0;
  mpt->lay[0] = pl;
  mpt->lay[0]->ninlist = 0;
  mpt->lay[0]->inlist = NULL;

  //
  //  Layer[1..nv1]  -  V1
  //
  for(i=0;i<nv1;i++){

    if (mv5->v1stf == NULL)  // New way
      sprintf(lname,"v1");  // OLD WAY
    else
      sprintf(lname,"v1_c%d",i);

    mod_me_filt_template(i,&fx,&fy,&fw,&fn,&zn);

    //
    //  Make a layer that gets input from the stim layer.
    //    The 'zn' dimension covers directions and phases.
    //
    //  WYETH Heterog 2019 *** HERE use pl->cell->attrib...  ETC.
    //    to set up heterogenous values.
    //
    pl = popu_make_layer_template(lname,mpt->xn,mpt->yn,zn,mpt->area[0],
				  mpt->lay[0],fx,fy,fw,fn);
    mpt->lay[1+i] = pl;

    free_2d_iarray(fx,zn);
    free_2d_iarray(fy,zn);
    free_2d_farray(fw,zn);  //** WYETH - BUG HERE - some null rows ???
    myfree(fn);
  }

  //  WYETH heterog 2019 We can search in the mv5->normo object to see if
  //    heterogeneous normalization parameters are set.
  //
  if (mv5->normo != NULL){
    int nhet,a2_seed,nrv,*shi;
    float a2,a2_max,rval;
    char **aname;
    struct pop_cell *c;

    a2     = onode_getpar_flt_exit(mv5->normo,"alpha_2");
    a2_max = onode_getpar_flt_dflt(mv5->normo,"alpha_2_max",a2);
    if (a2_max > a2){
      a2_seed = onode_getpar_flt_exit(mv5->normo,"alpha_2_seed");
      printf("    Using heterogeneous a2 values:  %f to %f\n",a2,a2_max);
      printf("    Seed = %d\n",a2_seed);

      // (1) Add an attrib parameter 'a2' to each layer (or just 1st ???)
      pl = mpt->lay[1];  // Lay 1 is 1st V1 layer
      nhet = 1;          // Assume 1 attrib to add to pop
      aname = (char **)myalloc(nhet*sizeof(char *));  // Make name list
      aname[0] = strdup("a2");
      pop_cell_attrib_init_f(pl,nhet,aname,1);  // 1=create cell attrib storage

      // (2) Get shuffled integers from 0 to nrv-1 (number of x,y locations)
      nrv = mpt->xn * mpt->yn;  // Number of random values needed
      shi = gr_get_shuffle_index(nrv,3,a2_seed);  // shuffled ints

      // Use the set of evenly spaced values from 'a2' to 'a2_max' and
      //   assign them randomly across the (x,y) array.
      k = 0;
      for(i=0;i<mpt->xn;i++){
	for(j=0;j<mpt->yn;j++){
	  rval = a2 + (float)shi[k]/(float)(nrv-1) * (a2_max - a2);
	  c = &(pl->c[i][j][0]);
	  pop_cell_attrib_set(c,"a2",rval,"MODEL_ME_V5_POP_PREP");
	  k += 1;
	}
      }
      printf("*** WYETH - next step is to make use of the 'a2' attribs\n");
      printf("*** WYETH - next step is to make use of the 'a2' attribs\n");
      printf("*** WYETH - next step is to make use of the 'a2' attribs\n");
      printf("*** WYETH - next step is to make use of the 'a2' attribs\n");
      exit_error("WYETH heterog HERE","Under development");
    }
  }
  


  //
  //  Layer[2] ...  -  MT population(s)
  //

  //  Storage for MT population params and responses
  mv5->mtp = (struct mod_me_v5_pop **)myalloc(popcnt*
					      sizeof(struct mod_me_v5_pop *));
  po = onode_get_next_type(m->o->o,"pop");
  for(i=0;i<popcnt;i++){

    mtpop = (struct mod_me_v5_pop *)myalloc(sizeof(struct mod_me_v5_pop));

    io    = onode_child_get_unique(po,"input");
    disto = onode_child_get_unique(io,"distrib");

    if (io == NULL)
      mylog_exit(mylogf,"MOD_ME_V5_POP_PREP  Missing <input> in <pop>.");
    if (disto == NULL)
      mylog_exit(mylogf,"MOD_ME_V5_POP_PREP  Missing <distrib> in <input>.");

    pname = onode_getpar_chr_exit(po,"name");
    pl = popu_make_pop_layer(pname);
    pl->laytype = strdup("default");
    mpt->lay[mv5->mtli0 + i] = pl;

    popu_init_geom(mylogf,m->o,po,pl,-1);  // Initialize layer geometry

    mtpop->ssum = get_4d_farray(pl->xn,pl->yn,pl->zn,0); // 3D array of ptrs
    mtpop->nl   = get_4d_farray(pl->xn,pl->yn,pl->zn,0); // 3D array of ptrs
    mtpop->seed = NULL;
    mtpop->wco  = NULL;
    mtpop->w    = NULL;
    mtpop->w_r  = NULL;

    pl->area = mpt->area[1];

    pl->c = popu_make_init_cells(pl->xn,pl->yn,pl->zn,pl); // Init cell array

    sprintf(ggstr,"    Created layer %s with %d x %d x %d units.\n",pl->name,
	    pl->xn,pl->yn,pl->zn);
    mylog(mylogf,ggstr);

    //
    //  Set layer attrib list to have 1 value: 'ori'
    //  Create storage in each cell to hold attribute value
    //
    aname = (char **)myalloc(sizeof(char *));
    aname[0] = strdup("ori");
    pop_cell_attrib_init_f(pl,1,aname,1);
    free_2d_carray(aname,1);

    //
    //  <me_v5>  (These params were formerly in a global <v5> structure
    //
    v5o = onode_child_get_unique(po,"me_v5");
    if (v5o != NULL){
      //
      //  Old way "<me_v5>" in <pop>
      //
      mtpop->r_bias  = onode_getpar_flt_dflt(v5o,"right_eye_w",1.0);
      mtpop->negwf   = 1.0;  // Wyeth - added later

      mtpop->w   = get_farray(mv5->n);
      mtpop->w_r = get_farray(mv5->n);

      if (onode_getpar_nvals(v5o,"w") != mv5->n)
	exit_error("MODEL_ME_V5_POP_PREP","Bad length of 'w' list");

      for(j=0;j<mv5->n;j++)
	mtpop->w[j] = onode_getpar_nth_flt_exit(v5o,"w",j);

      if (onode_item(v5o,"w_r") == 1){  // If R.E. weights are given

	if (onode_getpar_nvals(v5o,"w_r") != mv5->n)
	  exit_error("MODEL_ME_V5_POP_PREP","Bad length of 'w_r' list");

	for(j=0;j<mv5->n;j++)
	  mtpop->w_r[j] = onode_getpar_nth_flt_exit(v5o,"w_r",j);
      }else{
	for(j=0;j<mv5->n;j++)
	  mtpop->w_r[j] = mtpop->w[j];
      }
    }else{
      //
      //  New way "<mt_stf>"
      //
      v5o = onode_child_get_unique(po,"mt_stf");
      if (v5o == NULL)
	mylog_exit(mylogf,"MOD_ME_V5_POP_PREP  Missing <mt_stf> in <pop>.");

      mtpop->wco = onode_child_get_unique(v5o,"weight_config");

      mtpop->r_bias = onode_getpar_flt_dflt(v5o,"right_eye_w",1.0);
      mtpop->negwf  = onode_getpar_flt_dflt(v5o,"negative_w_scale",1.0);
    }

    //
    //  'v5o' maybe old or new way, but should have '<nonlin>' in it
    //
    nlo = onode_child_get_unique(v5o,"nonlin");
    if (nlo == NULL)
      mylog_exit(mylogf,"MOD_ME_V5_POP_PREP  Missing <nonlin> in <me_v5>.");

    mtpop->nl_type = onode_getpar_chr_exit(nlo,"type");
    if (strcmp(mtpop->nl_type,"rmsm")==0){
      mtpop->nl_a = onode_getpar_flt_exit(nlo,"A");
      mtpop->nl_b = onode_getpar_flt_exit(nlo,"B");
    }else if (strcmp(mtpop->nl_type,"half-rect")==0){
      mtpop->nl_a = onode_getpar_flt_dflt(nlo,"slope",1.0);
      mtpop->nl_b = onode_getpar_flt_dflt(nlo,"threshold",0.0);
    }else if (strcmp(mtpop->nl_type,"none")==0){
      mtpop->nl_a = 0.0;
      mtpop->nl_b = 0.0;
    }else{
      exit_error("MODEL_ME_V5_POP_PREP","Unknown 'type' in <nonlin>");
    }

    //
    //  MT-POW-NORM
    //
    nrmo = onode_child_get_unique(v5o,"norm");
    if (nrmo != NULL){
      mtpop->inorm_type   = onode_getpar_chr_exit(nrmo,"type");
      mtpop->inorm_pow    = onode_getpar_flt_exit(nrmo,"power");
      mtpop->inorm_pow1st = onode_getpar_int_dflt(nrmo,"power_first",0);
      mtpop->inorm_subu   = onode_getpar_chr_exit(nrmo,"subunit");
      mtpop->inorm_subwi  = onode_getpar_chr_exit(nrmo,"subunit_within");
    }else{
      mtpop->inorm_type   = NULL;
      mtpop->inorm_pow    = 1.0;
      mtpop->inorm_pow1st = 0;
      mtpop->inorm_subu   = strdup("none");
      mtpop->inorm_subwi  = strdup("sum");
    }

    //
    //  Inhibitory subunit
    //
    inso = onode_child_get_unique(v5o,"subunit");
    if (inso != NULL){
      mtpop->mtsub_name  = onode_getpar_chr_exit(inso,"pop_name");
      subpl = pop_util_get_layer_by_name(mpt,mtpop->mtsub_name);

      if (onode_getpar_nvals(inso,"w") != subpl->zn)
	exit_error("MODEL_ME_V5_POP_PREP","<subunit> Bad length of 'w' list");

      mtpop->mtsub_ndir = subpl->zn;
      mtpop->mtsub_w    = get_farray(subpl->zn);

      for(j=0;j<subpl->zn;j++){
	mtpop->mtsub_w[j] = onode_getpar_nth_flt_exit(inso,"w",j);
      }
    }else{
      mtpop->mtsub_ndir = 0;
      mtpop->mtsub_name = NULL;
      mtpop->mtsub_w = NULL;
    }


    //
    //  Input
    //
    pl->ninlist = nv1;  // One input for each STF channel
    pl->inlist = (struct onode **)myalloc(pl->ninlist*sizeof(struct onode *));

    if (mv5->v1stf == NULL){   // Old way
      pl->inlist[0] = io;      // Store a pointer to the onode
    }else{
      //
      // WYETH STF -  Make one input for each STF channel
      //
      for(j=0;j<pl->ninlist;j++){
	tio = onode_copy(io);
	sprintf(lname,"v1_c%d",j);
	uflag = onode_update_value(tio,"pop_origin",lname);
	if (uflag != 1)
	  exit_error("WYETH","lname update failed");

	pl->inlist[j] = tio;        // Store a pointer to the onode
      }
    }

    mod_me_v5_pop_conn(mtpop,pl,io,disto);

    //
    //  Spike_gen - just verify that it is OK
    //
    sgen = onode_child_get_unique(po,"spike_gen");
    if (onode_test_chr(sgen,"type","poisson")){
      ; // OK
    }else
      exit_error("MODEL_ME_V5_POP_PREP","Unknown <spike_gen> type");

    mv5->mtp[i] = mtpop;

    myfree(pname);

    po = onode_get_next_type(po->next,"pop");
  }

  mv5->mt_npop = popcnt;  // Number of MT layers

  if ((m->action == 10) && (m->marfile != NULL)){
    popu_mar_write(mpt,m->marfile,4);  // Version 4
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                        MOD_ME_V5_ADSTATE_READ_WRITE                       */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_adstate_read_write(cmd)
     char *cmd;  // 'read' or 'write'
{
  FILE *fopen(),*fp;
  int i,j,k,l;
  int ns,xn,yn,xi,yi,eyeflag,cnt,stfi,fcode,eyecode,nstf,ndir,nsch,fxn,fyn;
  float ***tst;
  char *sfile;
  char t1[SLEN],t2[SLEN];
  struct mod_me_v5_struct *mv5;
  struct mod_me_v1resp *v1s;
  struct mod_me_v1_input ***v1rr;

  mylog(mylogf,"  MOD_ME_V5_ADSTATE_READ_WRITE\n");

  mv5 = mod_me_v5;
  xn = mod_me_xn;
  yn = mod_me_yn;

  sfile = mv5->v1_adsfile;

  stfi = 0;  // WYETH HACK
  if (mv5->v1stf == NULL){
    v1rr = mv5->v1resp; // old way
  }else{
    v1rr = mv5->v1stf->sgr[stfi]->bdgr; // new way
  }

  if (mod_me_mpt == NULL){
    //
    //  No spatial summation, responses for the central point only
    //
    //if (eyeflag == 0)
    //v1s = v1rr[0][0]->dl;
    //else
    //v1s = v1rr[0][0]->dr;

    exit_error("MODEL_ME_V5_ADSTATE_READ_WRITE","Not implemented yet.");

  }else{
    //
    //  Spatial summation - compute responses at flagged V1 positions
    //

    if (mv5->v1stf != NULL){ 
      // new way
      exit_error("MODEL_ME_V5_ADSTATE_READ_WRITE","Not implemented yet.");
      /*
      for(ci=0;ci<mv5->v1stf->n;ci++){
	v1rr = mv5->v1stf->sgr[ci]->bdgr;

	for(i=0;i<xn;i++){
	  for(j=0;j<yn;j++){
	    if (v1rr[i][j]->flag == 2){
	      if (eyeflag == 0)
		v1s = v1rr[i][j]->dl;
	      else
		v1s = v1rr[i][j]->dr;

	    }
	  }
	}
      }*/
    } // else{
      // old way
  }

  eyeflag = 0;  // WYETH - HACK

  //
  //  Process header
  //
  if (strcmp(cmd,"write")==0){
    //
    //  State file format:
    //
    //  ME_V5_State [flag] // 1-sparse, only list units w/ non-zero state.
    //  eye  [flag] // 0-left, 1-right, 2-both
    //  nSTF [n]    // Number of STF channels
    //  ndir [n]    // Number of direction channels
    //  nsch [n]    // Number of sub-channels per direction channel
    //  [xn] [yn]   // spatial grid
    //  STF [i]     // STF channel
    //  [xi] [xj]  [dir0_sc0] [dir0_sc1] ... [dir1_sc0] [dir1_sc1] ... "\n"
    //     That is, all subchannels are listed for a given dir before next dir.

    if ((fp = fopen(sfile,"w")) == NULL){
      printf("  *** Cannot open file %s.\n",sfile);
      exit(0);
    }
    mylog(mylogf,"    Writing adaptation state file:\n");
    sprintf(ggstr,"      %s\n",mv5->v1_adsfile);
    mylog(mylogf,ggstr);

    fprintf(fp,"ME_V5_State 1\n");
    fprintf(fp,"eye %d\n",0);        // 0-left, 1-right, 2-both
    fprintf(fp,"nSTF %d\n",1);       // Number of STF channels
    fprintf(fp,"ndir %d\n",mv5->n);  // direction channels
    fprintf(fp,"nsch %d\n",1);       // sub-channels per direction
    fprintf(fp,"xn yn %d %d\n",xn,yn);  // spatial grid
    fprintf(fp,"STF %d\n",0);        // Start of first STF block

    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){

	if (mv5->v1srnd == NULL){

	  if (strcmp(mv5->v1_adtype,"RPHomeo")==0){
	    exit_error("MODEL_ME_V5_ADSTATE_READ_WRITE",
		       "Not implemented yet - w/o surround.");
	    //************ WYETH FIX THIS FOR FUTURE
	    //************ WYETH FIX THIS FOR FUTURE
	    //************ WYETH FIX THIS FOR FUTURE
	    //************ WYETH FIX THIS FOR FUTURE
	  }

	  //
	  //  No surround
	  //
	  if (v1rr[i][j]->flag == 2){
	    if (eyeflag == 0)
	      v1s = v1rr[i][j]->dl;
	    else
	      v1s = v1rr[i][j]->dr;

	    if (mv5->v1_dflag == 0){

	      fprintf(fp,"%3d %3d   ",i,j);
	      for(k=0;k<mv5->n;k++){ // For each dir channel
		fprintf(fp," %f",v1s->adst[k][0]);
	      }
	      fprintf(fp,"\n");
		  
	    }else{
	      exit_error("MODEL_ME_V5_ADSTATE_READ_WRITE",
			 "Not implemented yet.");
	    }
	  }
	}else{
	  //
	  //  Surround
	  //
	  if (eyeflag == 0)
	    tst = mv5->v1_adst_l;
	  else
	    tst = mv5->v1_adst_r;

	  fprintf(fp,"%3d %3d   ",i,j);

	  if (strcmp(mv5->v1_adtype,"s01")==0){
	    for(k=0;k<mv5->n;k++){ // For each dir channel
	      fprintf(fp," %f",tst[k][i][j]);
	    }
	  }else if (strcmp(mv5->v1_adtype,"RPHomeo")==0){
	    for(k=0;k<mv5->n;k++){ // For each dir channel
	      for(l=0;l<mv5->n;l++){ // For each dir channel
		fprintf(fp," %f",mv5->v1_ad_normw[i][j][k][l]);
	      }
	    }
	  }
	  fprintf(fp,"\n");
	}
      }
    }

  }else if (strcmp(cmd,"read")==0){

    if ((fp = fopen(sfile,"r")) == NULL){
      printf("  *** Cannot open file %s.\n",sfile);
      exit(0);
    }
    sprintf(ggstr,"    Reading from file:  %s\n",sfile);
    mylog(mylogf,ggstr);

    ns = fscanf(fp,"%s %d",t1,&fcode);     // 'fcode' should be 1
    ns = fscanf(fp,"%s %d",t1,&eyecode);   // 0-left, 1-right, 2-both
    ns = fscanf(fp,"%s %d",t1,&nstf);      // Number of STF channels
    ns = fscanf(fp,"%s %d",t1,&ndir);      // direction channels
    ns = fscanf(fp,"%s %d",t1,&nsch);      // sub-channels per direction
    ns = fscanf(fp,"%s %s %d %d",t1,t2,&fxn,&fyn);  // spatial grid

    cnt = 0;  // Count number of locations

    while (fscanf(fp,"%s %s",t1,t2) != EOF){
      if (strcmp(t1,"STF")==0){
	stfi = atoi(t2);
	//printf("found stfi = %d\n",stfi);
      }else{
	xi = atoi(t1);
	yi = atoi(t2);
	//printf("found xi,yi = %d %d\n",xi,yi);

	if (mv5->v1srnd == NULL){
	  //
	  //  No surround
	  //

	  if (strcmp(mv5->v1_adtype,"RPHomeo")==0){
	    exit_error("MODEL_ME_V5_ADSTATE_READ_WRITE",
		       "Not implemented yet - w/o surround.");
	    //************ WYETH FIX THIS FOR FUTURE
	    //************ WYETH FIX THIS FOR FUTURE
	    //************ WYETH FIX THIS FOR FUTURE
	    //************ WYETH FIX THIS FOR FUTURE
	  }

	  if (v1rr[xi][yi]->flag != 2){
	    exit_error("MOD_ME_V5_ADSTATE_READ_WRITE","Expecting flag 2");
	  }
	  if (eyeflag == 0)
	    v1s = v1rr[xi][yi]->dl;
	  else
	    v1s = v1rr[xi][yi]->dr;

	  for(i=0;i<ndir;i++){ // For each dir channel
	    ns = fscanf(fp,"%f",&(v1s->adst[i][0]));
	  }
	}else{
	  //
	  //  Read into the dense storage
	  //

	  if (strcmp(mv5->v1_adtype,"s01")==0){
	    if (eyeflag == 0)
	      tst = mv5->v1_adst_l;
	    else
	      tst = mv5->v1_adst_r;
	    
	    for(i=0;i<ndir;i++){ // For each dir channel
	      ns = fscanf(fp,"%f",&(tst[i][xi][yi]));
	    }
	  }else if (strcmp(mv5->v1_adtype,"RPHomeo")==0){
	    for(k=0;k<mv5->n;k++){ // For each dir channel
	      for(l=0;l<mv5->n;l++){ // For each dir channel
		ns = fscanf(fp," %f",&(mv5->v1_ad_normw[xi][yi][k][l]));
	      }
	    }
	  }
	}
	cnt += 1;
      }
    }

    sprintf(ggstr,"    Read states for %d locations.\n",cnt);
    mylog(mylogf,ggstr);

    if (mv5->v1srnd != NULL){
      if (cnt < xn*yn)
	exit_error("MOD_ME_V5_ADSTATE_READ_WRITE","Need dense state");
    }
	
  }else{
    exit_error("MOD_ME_V5_ADSTATE_READ_WRITE","Unknown command");
  }

  fclose(fp);
}
/**************************************-**************************************/
/*                                                                           */
/*                          MODEL_ME_V5_FREE_V1RSP                           */
/*                                                                           */
/*****************************************************************************/
void model_me_v5_free_v1rsp(v)
     struct mod_me_v1resp *v;
{
  int ndir,nc,tn;

  ndir = mod_me_v5->n;
  nc   = mod_me_v5->v1_nc;
  tn   = mod_me_tn;

  if (v->raw   != NULL)  free_3d_farray(v->raw,ndir,nc,tn);
  if (v->me    != NULL)  free_2d_farray(v->me,ndir);
  if (v->sum   != NULL)  myfree(v->sum);
  if (v->norm  != NULL)  free_3d_farray(v->norm,ndir,nc,tn);
  if (v->opp   != NULL)  free_3d_farray(v->opp,ndir,nc*2,tn);
  if (v->bde   != NULL)  free_3d_farray(v->bde,ndir,7,tn);
  if (v->adst  != NULL)  free_2d_farray(v->adst,ndir);
  if (v->adrsp != NULL)  free_3d_farray(v->adrsp,ndir,nc,tn);
  if (v->mt    != NULL)  myfree(v->mt);

  myfree(v);
}
/**************************************-**************************************/
/*                                                                           */
/*                          MODEL_ME_V5_FREE_V1INP                           */
/*                                                                           */
/*****************************************************************************/
void model_me_v5_free_v1inp(v)
     struct mod_me_v1_input *v;
{
  if (v->flag == 0)
    return;  // There is no storage

  if (v->dl  != NULL)  model_me_v5_free_v1rsp(v->dl);
  if (v->dr  != NULL)  model_me_v5_free_v1rsp(v->dr);
  if (v->mtb != NULL)  myfree(v->mtb);

  myfree(v);
}
/**************************************-**************************************/
/*                                                                           */
/*                          MODEL_ME_V5_FREE_MTPOP                           */
/*                                                                           */
/*****************************************************************************/
void model_me_v5_free_mtpop(p,pl)
     struct mod_me_v5_pop *p;   // V5_ME pop
     struct pop_layer *pl;      // layer in 'ifc.h'

{
  int lxn,lyn,lzn,tn;

  tn = mod_me_tn;

  lxn = pl->xn;
  lyn = pl->yn;
  lzn = pl->zn;

  free_4d_farray(p->ssum,lxn,lyn,lzn,tn);
  free_4d_farray(p->nl,lxn,lyn,lzn,tn);
  free_3d_iarray(p->seed,lxn,lyn,lzn);

  if (p->w   != NULL)  myfree(p->w);
  if (p->w_r != NULL)  myfree(p->w_r);

  if (p->nl_type     != NULL)  myfree(p->nl_type);
  if (p->inorm_type  != NULL)  myfree(p->inorm_type);
  if (p->inorm_subu  != NULL)  myfree(p->inorm_subu);
  if (p->inorm_subwi != NULL)  myfree(p->inorm_subwi);

  if (p->mtsub_w     != NULL)  myfree(p->mtsub_w);
  if (p->mtsub_name  != NULL)  myfree(p->mtsub_name);

  //printf("lxn = %d  lyn = %d  lzn = %d\n",lxn,lyn,lzn);
}
/**************************************-**************************************/
/*                                                                           */
/*                              MODEL_ME_V5_FREE                             */
/*                                                                           */
/*****************************************************************************/
void model_me_v5_free()
{
  int i,j;
  int xn,yn,tn,ndir;
  struct mod_me_v5_struct *mv5;

  mylog(mylogf,"  MODEL_ME_V5_FREE\n");

  mv5 = mod_me_v5;  // Short name

  xn = mod_me_xn;
  yn = mod_me_xn;
  tn = mod_me_xn;
  ndir = mv5->n;

  if (mv5->w          != NULL)  myfree(mv5->w);
  if (mv5->w_r        != NULL)  myfree(mv5->w);
  if (mv5->v1_mix_s   != NULL)  myfree(mv5->v1_mix_s);
  if (mv5->mt_nltype  != NULL)  myfree(mv5->mt_nltype);
  if (mv5->v1_adtype  != NULL)  myfree(mv5->v1_adtype);
  if (mv5->v1_adp_tv  != NULL)  myfree(mv5->v1_adp_tv);
  if (mv5->v1_adrw    != NULL)  myfree(mv5->v1_adrw);
  if (mv5->v1_adsfile != NULL)  myfree(mv5->v1_adsfile);

  if (mv5->v1_adst_l != NULL)
    free_3d_farray(mv5->v1_adst_l,ndir,xn,yn);
  if (mv5->v1_adst_r != NULL)
    free_3d_farray(mv5->v1_adst_r,ndir,xn,yn);

  if (mv5->v1_ad_normw != NULL)
    free_4d_farray(mv5->v1_ad_normw,xn,yn,ndir,ndir);

  if (mv5->v1srnd != NULL){

    if (mv5->v1srnd->sl != NULL)
      free_4d_farray(mv5->v1srnd->sl,ndir,xn,yn,tn);

    if (mv5->v1srnd->sr != NULL)
      free_4d_farray(mv5->v1srnd->sr,ndir,xn,yn,tn);

    if (mv5->v1srnd->rule != NULL)
      myfree(mv5->v1srnd->rule);

    myfree(mv5->v1srnd);
  }

  if (mv5->v1resp != NULL){

    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){
	model_me_v5_free_v1inp(mv5->v1resp[i][j]);
      }
      myfree(mv5->v1resp[i]);
    }
    myfree(mv5->v1resp);
  }

  if (mv5->v1stf != NULL){
    exit_error("WYETH HERE","Under development - V1STF not implemented.");
  }

  if (mv5->mt_norm != NULL)
    myfree(mv5->mt_norm);

  if (mv5->mtp != NULL){
    //
    //  NOTE, 'mv5->mtp[i]' matches 'mod_pop_top->lay[i0 + i]'
    //
    for(i=0;i<mv5->mt_npop;i++){
      model_me_v5_free_mtpop(mv5->mtp[i],mod_me_mpt->lay[mv5->mtli0 + i]);
    }

    myfree(mv5->mtp);
  }

  if (mv5->noise_v1_out != NULL){
    if (mv5->noise_v1_out->type != NULL)
      myfree(mv5->noise_v1_out->type);
    if (mv5->noise_v1_out->seedlist != NULL)
      myfree(mv5->noise_v1_out->seedlist);

    myfree(mv5->noise_v1_out);
  }

  if (mv5->ftdata != NULL)
    free_2d_farray(mv5->ftdata,ndir);

  myfree(mv5);

  //exit_error("WYETH HERE","Under development 'v1resp'");
}
/**************************************-**************************************/
/*                                                                           */
/*                              MODEL_ME_V5_PREP                             */
/*                                                                           */
/*****************************************************************************/
void model_me_v5_prep(m,s)
     struct model_struct *m;  // Model params
     struct stim_struct *s;   // Stimulus params
{
  int i,j,k;
  int xn,yn,tn,write_filt,fn,***cmap,txn,tyn,tzn,nsyn,nv1,normtype,ndir;
  float ***fpe,***fpo,***fne,***fno,*w,fsc,sscale,tscale,tval;
  char *oppstr;
  char tstr[SLEN];
  struct onode *fo,*v1o,*mto,*sgen,*io,*disto,*fso,*v1sr,*v1ad,*no;
  struct mod_me_fstruct **fs;  // [fn] List of pointers to filters
  struct mod_me_v5_struct *mv5;
  struct pop_layer *pl,*pls;
  struct pop_top *mpt;
  struct mod_me_v5_pop *mtpop;

  mylog(mylogf,"  MODEL_ME_V5_PREP\n");

  mod_me_prep_basic(m->o,&sscale,&tscale,&xn,&yn,&tn); // Get and set globals

  mod_me_tsi = -1;

  if (m->ppl != NULL)
    exit_error("MODEL_ME_V5_PREP","m->ppl is not NULL");

  //
  //  Create structure for ME_V5 model params and responses
  //
  mv5 = (struct mod_me_v5_struct *)myalloc(sizeof(struct mod_me_v5_struct));
  mod_me_v5 = mv5;
  mv5->ftdata = NULL;

  v1o = onode_child_get_unique(m->o,"v1");
  fso = onode_child_get_unique(v1o,"filters");

  if (fso != NULL){
    //
    //  New 'STF' way
    //
    mv5->n = onode_getpar_int_exit(fso,"n_dir");

    mod_me_stf_config(m);  // 'mod_me_ftop' is created
    mod_me_fs_build_list();

    fsc = onode_getpar_flt_dflt(fso,"norm_scale",1.0); // Filter scale const.
    normtype = onode_getpar_int_exit(fso,"norm_type"); // Normalization type

    // *** This will cause the raw responses to be passed through w/o
    //   opponency in worker routines called by 'mod_me_v5_raw_top'.
    mod_me_oppflag = 3; // Non-opponent; WYETH - is this used for ME_V5 mod ??
  }else{
    //
    //  OLD WAY
    //
    mv5->n     = onode_getpar_int_exit(v1o,"n_dir_chan");
    mv5->v1stf = NULL;

    mylog(mylogf,"    Creating all filters\n");

    // (1) Set xn, yn, tn.
    // (2) Create even and odd filters for every direction.
    fs = mod_me_get_quad_opp_filters_o(m->o,&mod_me_fpe,&mod_me_fpo,
				       &mod_me_fpx,&mod_me_fne,
				       &mod_me_fno,&mod_me_fnx,
				       &fn); // Number in '..._ff'
    mod_me_ff = fs;  // [fn] Store list of filter structures and
    mod_me_fn = fn;  // number of filter structs globally


    // Needed for normalization below
    fo = onode_child_get_unique(m->o,"filter");
    fsc = onode_getpar_flt_dflt(fo,"scale_sqrt",1.0);
    if (fsc != 1.0)
      fsc = 1.0/fsc;
    // Original way, before 'norm_type' added

    normtype = onode_getpar_int_dflt(fo,"norm_type",2);
  }
  ndir = mv5->n;

  modu_mds_add_3d(m,"v1","raw" ,'f',xn,yn,ndir,NULL,"full","no",NULL,-1);
  modu_mds_add_3d(m,"v1","norm",'f',xn,yn,ndir,NULL,"full","no",NULL,-1);

  modu_mds_add_3d(m,"v1r","raw" ,'f',xn,yn,ndir,NULL,"full","binoc",NULL,-1);
  modu_mds_add_3d(m,"v1r","norm",'f',xn,yn,ndir,NULL,"full","binoc",NULL,-1);


  //
  //  Normalize the filters.
  //
  mod_me_normalize_filters(1,xn,1,yn,1,tn,sscale,tscale,fsc,normtype);

  if (fso == NULL){  // Only works for old way

    //
    //  Write out filters for inspection, if requested.
    //
    write_filt = onode_getpar_int_dflt(fo,"write_filter",-1);
    if ((write_filt > 0) && write_filt <= 5){
      fpe = fs[0]->f;
      fpo = fs[1]->f;
      fne = fs[2]->f;
      fno = fs[3]->f;
  
      if ((write_filt == 1)||(write_filt == 5))
	write_3d_data_part("zzz.fpe.3d",fpe,1,xn,1,yn,1,tn,4,2,1);
      if ((write_filt == 2)||(write_filt == 5))
	write_3d_data_part("zzz.fpo.3d",fpo,1,xn,1,yn,1,tn,4,2,1);
      if ((write_filt == 3)||(write_filt == 5))
	write_3d_data_part("zzz.fne.3d",fne,1,xn,1,yn,1,tn,4,2,1);
      if ((write_filt == 4)||(write_filt == 5))
      write_3d_data_part("zzz.fno.3d",fno,1,xn,1,yn,1,tn,4,2,1);

      printf("  Exiting after writing filter(s) to .3d file\n");
      printf("  Use 'norm3d <filename> <scale_factor>' to view.\n");
      printf("  You must have X-windows (X11) running.\n");
      printf("\n");
      exit(0);
    }
    write_filt = onode_getpar_int_dflt(fo,"write_filter_index",-1);
    if (write_filt >= 0){
      if (write_filt >= mod_me_fn)
	exit_error("MODEL_ME_V5_PREP","'write_filer_index' too large");

      fpe = fs[write_filt]->f;
      sprintf(tstr,"zzz.filt_%d.3d",write_filt);
      write_3d_data_part(tstr,fpe,1,xn,1,yn,1,tn,4,2,1);

      printf("  Exiting after writing filter to .3d file\n");
      printf("  Use 'norm3d <filename> <scale_factor>' to view.\n");
      printf("  You must have X-windows (X11) running.\n");
      printf("\n");
      exit(0);
    }

    //
    //  Set 'mod_me_oppflag' which is used, e.g., by 'mod_me_v5_raw_resp'
    //    *** Note: 0 is OPPONENT, flag is initialized to 0
    //
    oppstr = onode_getpar_chr_ptr(fo,"opponent");
    if (oppstr != NULL){
      if (strcmp(oppstr,"no")==0){
	mod_me_oppflag = 3; // Non-opponent, 
      }
    }
  }

  //  'mto' will NOT be NULL if this is the old format (before MT <pop>)
  mto = onode_child_get_unique(m->o,"v5");
  if (mto != NULL){
    //
    //  Verify 'sgen'
    //
    sgen = onode_child_get_unique(m->o,"spike_gen");
    if (onode_test_chr(sgen,"type","poisson")){
      ; // OK
    }else
      exit_error("MODEL_ME_V5_PREP","Unknown <spike_gen> type");
  }


  //
  //  Get V1 configuration
  //
  mv5->normo = onode_child_get_unique(v1o,"norm");
  if (mv5->normo == NULL){
    //
    //  Old way
    //
    mv5->alph1     = onode_getpar_flt_exit(v1o,"rmsm_alpha_1");
    mv5->alph2     = onode_getpar_flt_exit(v1o,"rmsm_alpha_2");
    mv5->alph3     = onode_getpar_flt_dflt(v1o,"rmsm_alpha_3",0.0);
    mv5->alph4     = onode_getpar_flt_dflt(v1o,"alpha_4",0.0);
    if (mv5->alph3 != 0.0)
      mv5->L         = onode_getpar_flt_exit(v1o,"rmsm_L");
    else
      mv5->L = 0.0;
  }

  mv5->v1_opp_w    = onode_getpar_flt_dflt(v1o,"opponent_w",0.0);
  mv5->v1_recto    = onode_getpar_int_dflt(v1o,"opponent_rect",0);
  mv5->v1_opp_tmin = onode_getpar_flt_dflt(v1o,"opponent_min_tf",0.0);
  if (mv5->v1_opp_tmin < 0.0)
    exit_error("MODEL_ME_V5_PREP",
	       "The value for 'opponent_min_tf' should be >= 0.0");
  mv5->v1_bin_w    = onode_getpar_flt_dflt(v1o,"right_eye_w",1.0);
  mv5->v1_mix_f    = onode_getpar_flt_dflt(v1o,"binoc_balance",1.0);

  //
  //  V1 surround suppression
  //
  v1sr = onode_child_get_unique(v1o,"surround");
  if (v1sr != NULL){
    mv5->v1srnd = (struct mod_me_v1_surr *)myalloc(
					   sizeof(struct mod_me_v1_surr));
    mv5->v1srnd->sig   = onode_getpar_flt_exit(v1sr,"sigma");
    mv5->v1srnd->sigc  = onode_getpar_flt_dflt(v1sr,"sigma_idiam",-1.0);
    mv5->v1srnd->ampl  = onode_getpar_flt_exit(v1sr,"ampl");
    mv5->v1srnd->omix  = onode_getpar_flt_exit(v1sr,"oppmix_w");
    mv5->v1srnd->delay = onode_getpar_flt_exit(v1sr,"delay");
    mv5->v1srnd->dumpi = onode_getpar_int_dflt(v1sr,"write_raw_index",-1);
    mv5->v1srnd->rule  = onode_getpar_chr_exit(v1sr,"rule");
    mv5->v1srnd->n     = ndir;
    mv5->v1srnd->sl = (float ****)myalloc(mv5->v1srnd->n*sizeof(float ***));
    for(i=0;i<mv5->v1srnd->n;i++)
      mv5->v1srnd->sl[i] = get_3d_farray(xn,yn,tn);
    mv5->v1srnd->sr = NULL; // When do we know if there are also RE signals?
  }else{
    mv5->v1srnd = NULL;
  }


  if (onode_get_item_child(v1o,"early_mix_f") != NULL){
    // This error message tells users to stop using the old way...
    // Used to be 0-0.5, now mapped to 1-0.5 (1 is now pure, 0.5 is 50/50 mix)
    exit_error("MODEL_ME_V5_PREP","Convert 'early_mix_f' to 'binoc_balance'");
  }

  mv5->v1_dflag  = onode_getpar_int_dflt(v1o,"disparity_flag",0);

  mv5->v1_mix_s  = onode_getpar_chr_dflt(v1o,"early_mix_stage","none");
  if ((mv5->v1_dflag == 1) && (strcmp(mv5->v1_mix_s,"none")!=0)){
    // Don't allow 'early_mix_stage' to work for the disparity model
    printf("*** disparity_flag = 1\n");
    exit_error("MODEL_ME_V5_PREP","'early_mix_stage' should be 'none'");
  }

  mv5->v1_dnorm  = onode_getpar_int_dflt(v1o,"disparity_norm",0);
  mv5->binoc = 0;  // Assume left eye only, until stimulus is checked

  if ((mv5->v1_mix_f < 0.5) || (mv5->v1_mix_f > 1.0))
    exit_error("MODEL_ME_V5_PREP","'binoc_balance' must be in [0.5,1.0]");

  //
  //  Noise
  //
  mv5->noise_v1_out = NULL;
  no = onode_get_node_type_item_val(v1o,"noise","name","output");
  if (no != NULL){
    mv5->noise_v1_out = mod_me_noise_init(no,s->ntr);
  }


  //
  //  Get MT configuration
  //
  if (mto != NULL){
    //
    //  OLD WAY - keep for models w/o MT pop ? ?
    //
    printf("  *** WARNING:  Old format file contains a top-level <v5> node.\n");

    mv5->w   = get_farray(ndir);  // We are now allocated *both* weight
    mv5->w_r = get_farray(ndir);  //   arrays, even if one one set of weights.

    if (onode_getpar_nvals(mto,"rmsm_w") != ndir)
      exit_error("MODEL_ME_V5_PREP","Bad length of 'rmsm_w' list");

    for(i=0;i<ndir;i++)
      mv5->w[i] = onode_getpar_nth_flt_exit(mto,"rmsm_w",i);

    if (onode_item(mto,"rmsm_w_r") == 1){  // If R.E. weights are given

      if (onode_getpar_nvals(mto,"rmsm_w_r") != ndir)
	exit_error("MODEL_ME_V5_PREP","Bad length of 'rmsm_w_r' list");

      for(i=0;i<ndir;i++)
	mv5->w_r[i] = onode_getpar_nth_flt_exit(mto,"rmsm_w_r",i);
    }else{
      for(i=0;i<ndir;i++)
	mv5->w_r[i] = mv5->w[i];
    }

    mv5->mt_a = onode_getpar_flt_exit(mto,"rmsm_A");
    mv5->mt_b = onode_getpar_flt_exit(mto,"rmsm_B");
    if (mv5->mt_a > 0.0)
      mv5->mt_nltype = strdup("rmsm");  // WYETH - for compat. w/ new pop way
    else
      mv5->mt_nltype = strdup("none");


    mv5->bin_hack = onode_getpar_flt_dflt(mto,"tailby_weight",-1.0);

  }else{
    mv5->w   = NULL;  // Flag to indicate that MT values are not here.
    mv5->w_r = NULL;  // Flag to indicate that MT values are not here.
    mv5->mt_nltype = NULL;
  }


  if (mv5->v1_dflag == 1)
    mv5->v1_nc = 2;  // To store even and odd responses separately
  else
    mv5->v1_nc = 1;  // Combine even and odd early


  //
  //  Spatial Summation
  //
  model_me_v5_pop_prep(m); // .marfile written here; why not out here at end?

  if (mv5->v1stf != NULL)  // New way
    nv1 = mv5->v1stf->n;
  else
    nv1 = 1;

  mpt = mod_me_mpt;

  if (mpt != NULL){  // If there is an MT population

    if (mv5->v1stf != NULL){  // New way
      pls = pop_util_get_layer_by_name(mpt,"v1_c0");   // V1 pop
    }else{
      pls = pop_util_get_layer_by_name(mpt,"v1");   // old way
    }
    //  Scan all inputs, and make storage where needed
    cmap = get_zero_3d_iarray(xn,yn,pls->zn);

    for(i=mv5->mtli0;i<mpt->nlay;i++){
      pl = mpt->lay[i]; // MT layer (post-syn)

      for(j=1;j<=nv1;j++){
	pls = mpt->lay[j];  // V1 layer (pre-syn)
	nsyn = popc_flag_accum_lay_to_lay(pls,pl,cmap,0.0);
	if (i == mv5->mtli0)
	  sprintf(ggstr,"    %d V1 units (of %d) are pre-syn to MT pop %s.\n",
		  nsyn,xn*yn,pl->name);
	else
	  sprintf(ggstr,"    %d additional V1 units pre-syn to MT pop %s.\n",
		  nsyn,pl->name);
	mylog(mylogf,ggstr);
      }
    }

    txn = xn;
    tyn = yn;
    tzn = pls->zn;
  }else{
    txn = tyn = tzn = 1;
    cmap = get_3d_iarray(1,1,1);  // Dummy
    cmap[0][0][0] = 1;
  }

  //
  //  Initialize adaptation parameters
  //    *** NOTE, this must set 'adflag' before response storage alloc. below
  //

  //  Set default values (important to know what to free)
  mv5->v1_adflag  = 0;
  mv5->v1_adtype  = NULL;
  mv5->v1_adp_tv  = NULL;
  mv5->v1_adrw    = NULL;
  mv5->v1_adsfile = NULL;
  mv5->v1_adst_l = mv5->v1_adst_r = NULL;
  mv5->v1_ad_normw = NULL;

  v1ad = onode_child_get_unique(v1o,"adapt");
  mv5->v1ado = v1ad;
  if (v1ad != NULL){
    //
    //  Now that we keep a pointer to the V1 adapt object, we do not have
    //  to read out all param values here.
    //

    mv5->v1_adflag  = onode_getpar_int_exit(v1ad,"flag");
    mv5->v1_adtype  = onode_getpar_chr_exit(v1ad,"type");
    if (strcmp(mv5->v1_adtype,"s01")==0){
      mv5->v1_adp_t1  = onode_getpar_flt_exit(v1ad,"tau");  // sec
      mv5->v1_adp_a1  = onode_getpar_flt_exit(v1ad,"chalf");
    }else if (strcmp(mv5->v1_adtype,"RPHomeo")==0){
      mv5->v1_adp_t1  = onode_getpar_flt_exit(v1ad,"alpha");

      if (onode_getpar_nvals(v1ad,"target") != (ndir/2 + 1)){
	printf("  *** List length should be (ndir/2 + 1) = %d\n",ndir/2+1);
	exit_error("MODEL_ME_V5_PREP","Bad length of 'target' list");
      }

      mv5->v1_adp_tv = (float *)myalloc((ndir/2+1)*sizeof(float));
      for(j=0;j<(ndir/2 + 1);j++)
	mv5->v1_adp_tv[j] = onode_getpar_nth_flt_exit(v1ad,"target",j);

      tval = onode_getpar_flt_dflt(v1ad,"target_scale",1.0);
      multiply_farray(mv5->v1_adp_tv,(ndir/2 + 1),tval);

      tval = onode_getpar_flt_dflt(v1ad,"target_offset",0.0);
      add_const_farray(mv5->v1_adp_tv,(ndir/2 + 1),tval);

      //printf("____________________________xn,yn ndir %d %d %d\n",xn,yn,ndir);
      mv5->v1_ad_normw = get_const_4d_farray(xn,yn,ndir,ndir,1.0); // RPHomeo
//exit_error("HERE WYETH","weights set to 11111111111111111111");
      
      //print_farray(mv5->v1_adp_tv,(ndir/2 + 1),"TARGET VALUES");
    }else{
      mylog_exit(mylogf,"MODEL_ME_V5_PREP  Unknown adapt type");
    }
    mv5->v1_adrw    = onode_getpar_chr_dflt(v1ad,"read_write","null");
    mv5->v1_adsfile = onode_getpar_chr_dflt(v1ad,"state_file","null");
    mv5->v1_adp_wwi = onode_getpar_int_dflt(v1ad,"write_weights",0);

    if (mv5->v1srnd != NULL){
      //
      //  If adaptation AND surround, we need full grid of adapt states.
      //
      mv5->v1_adst_l = get_zero_3d_farray(ndir,xn,yn);
      mv5->v1_adst_r = get_zero_3d_farray(ndir,xn,yn);
    }
  }


  //
  //  Allocate storage for V1 responses
  //
  if (mv5->v1stf != NULL){
    //
    //  WYETH STF - New way
    //
    for(i=0;i<mv5->v1stf->n;i++){
      mv5->v1stf->sgr[i] = (struct mod_me_v1_spatgrp *)myalloc(
			                    sizeof(struct mod_me_v1_spatgrp));
      mv5->v1stf->sgr[i]->nx = txn;
      mv5->v1stf->sgr[i]->ny = tyn;
      // This next line will make storage for both LE and RE, whether
      //   we need the RE or not - related to WYETH varmoo_bug
      mv5->v1stf->sgr[i]->bdgr = mod_me_bdgr_2d_init(txn,tyn,tzn,cmap);
    }
  }else{
    //  ORIG WAY
    mv5->v1resp = mod_me_bdgr_2d_init(txn,tyn,tzn,cmap);
  }
  free_3d_iarray(cmap,txn,tyn,tzn);


  //
  //  Read adaptation state from file.
  //    *** NOTE, storage must have already been created above.
  //
  if (v1ad != NULL){
    if (strcmp(mv5->v1_adrw,"read")==0){
      mod_me_v5_adstate_read_write("read");
    }
  }

  mv5->mt_norm = NULL;  // default
  if (mpt == NULL){
    mv5->mt_norm  = (float  *)myalloc(tn*sizeof(float));
  }
  //  Storage for MT population - was already done in '..._v5_pop_prep()'


  mylog(mylogf,"    ME_V5 params\n");
  sprintf(ggstr,"      N_dir = %d\n",ndir);             mylog(mylogf,ggstr);

  if (mv5->normo == NULL){
    sprintf(ggstr,"      alpha_1 = %f\n",mv5->alph1);     mylog(mylogf,ggstr);
    sprintf(ggstr,"      alpha_2 = %f\n",mv5->alph2);     mylog(mylogf,ggstr);
    sprintf(ggstr,"      alpha_3 = %f\n",mv5->alph3);     mylog(mylogf,ggstr);
    sprintf(ggstr,"      alpha_4 = %f\n",mv5->alph4);     mylog(mylogf,ggstr);
  }else{
    // New way
    ; // WYETH - print out norm constants here ...
  }

  sprintf(ggstr,"      opponent_w = %f\n",mv5->v1_opp_w); mylog(mylogf,ggstr);

  if (mod_me_mpt == NULL){
    sprintf(ggstr,"      r_eye_w = %f\n",mv5->v1_bin_w);  mylog(mylogf,ggstr);
    sprintf(ggstr,"      MT nonlin A = %f\n",mv5->mt_a);  mylog(mylogf,ggstr);
    sprintf(ggstr,"      MT nonlin B = %f\n",mv5->mt_b);  mylog(mylogf,ggstr);
    mylog(mylogf,"      Weights:\n       ");
    for(i=0;i<ndir;i++){
      sprintf(ggstr," %.2f",mv5->w[i]);  mylog(mylogf,ggstr);
    }
    mylog(mylogf,"\n");
  }else{
    sprintf(ggstr,"      MT populations: %d\n",mv5->mt_npop);
    mylog(mylogf,ggstr);

    for(j=0;j<mv5->mt_npop;j++){
      mtpop = mv5->mtp[j];
      sprintf(ggstr,"       Name: %s\n",mod_me_mpt->lay[2+j]->name);
      mylog(mylogf,ggstr);
      sprintf(ggstr,"         r_eye_w = %f\n",mtpop->r_bias);
      mylog(mylogf,ggstr);
      sprintf(ggstr,"         MT nonlin A = %f\n",mtpop->nl_a);
      mylog(mylogf,ggstr);
      sprintf(ggstr,"         MT nonlin B = %f\n",mtpop->nl_b);
      mylog(mylogf,ggstr);
      if (mtpop->w != NULL){
	mylog(mylogf,"         Weights:\n       ");
	for(i=0;i<ndir;i++){
	  sprintf(ggstr," %.2f",mtpop->w[i]);  mylog(mylogf,ggstr);
	}
	mylog(mylogf,"\n");
      }
    }
  }

  if (m->action == 20){
    modu_mds_write(m,NULL,0);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MODEL_ME_V5_DONE                             */
/*                                                                           */
/*****************************************************************************/
void model_me_v5_done(m,s,r)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response data
{
  struct mod_me_v5_struct *mv5;

  mylog(mylogf,"  MODEL_ME_V5_DONE\n");

  mv5 = mod_me_v5;

  //if (mv5->v1_adflag == 1){
  if ((mv5->v1_adflag == 1) || (mv5->v1_adflag == 2)){
    if (strcmp(mv5->v1_adrw,"write")==0){
      mod_me_v5_adstate_read_write("write");
    }
  }

  //
  //  To free an ME_V5 model, we have to free three main things:
  //
  //  (1) The front end filters
  //  (2) The contents of 'mod_me_v5' which is a 'mod_me_v5_struct'
  //  (3) The contents of 'mod_me_mpt' which is a 'pop_top' struct.
  //

  mod_me_filters_free();  // WYETH HERE NEW

  model_me_v5_free();

  popu_free_pop_top(mod_me_mpt);


  mylog(mylogf,"  MODEL_ME_V5_DONE_end\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                           MODEL_FILTER_XYT_PREP                           */
/*                                                                           */
/*****************************************************************************/
void model_filter_xyt_prep(m)
     struct model_struct *m; // Model params
{
  int i,j,k;
  int xn,yn,tn,write_filt;
  float **xyfilt,*tfilt,sscale,tscale;
  struct onode *mo,*fso,*fto;

  mylog(mylogf,"  MODEL_FILTER_XYT_PREP\n");

  mo = m->o;

  if (mo != NULL)
    mod_me_modtype = onode_getpar_chr_exit(m->o,"mod_type");
  else
    exit_error("MODEL_FILTER_XYT_PREP","Head onode is null");

  sscale = onode_getpar_flt_exit(mo,"sscale");
  tscale = onode_getpar_flt_exit(mo,"tscale");
  xn     = onode_getpar_int_exit(mo,"xn");
  yn     = onode_getpar_int_exit(mo,"yn");
  tn     = onode_getpar_int_exit(mo,"tn");

  mod_me_sscale = sscale;
  mod_me_tscale = tscale;
  mod_me_xn     = xn;
  mod_me_yn     = yn;
  mod_me_tn     = tn;

  //  Get the spatial filter
  fso = onode_get_node_type_item_val(mo,"filter","name","spatial");
  xyfilt = mod_util_o_get_filter_2d(mylogf,fso,sscale,&xn,&yn);

  //  Get the temporal filter
  fto = onode_get_node_type_item_val(mo,"filter","name","temporal");
  tfilt = kernu_o_filter(mylogf,fto,tscale,&tn);

  mod_me_fpe = f3tensor(1,xn,1,yn,1,tn);  // Use tensor for 3D FFT
  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++)
      for(k=0;k<tn;k++)
	mod_me_fpe[i+1][j+1][k+1] = xyfilt[i][j] * tfilt[k];


  // WYETH - SHOULD PUT FILTER IN A LIST ??? ****/
  // WYETH - SHOULD PUT FILTER IN A LIST ??? ****/
  // WYETH - SHOULD PUT FILTER IN A LIST ??? ****/
  // and then use ...._filter(... instead
  mod_me_normalize_filter(mod_me_fpe,1,xn,1,yn,1,tn,
			  mod_me_sscale,mod_me_tscale,1.0);

  write_filt = onode_getpar_int_dflt(mo,"model_write_filter",-1);
  if (write_filt == 1)
    write_3d_data_part("zzz.filter.3d",mod_me_fpe,1,xn,1,yn,1,tn,4,2,1);

  mod_me_pow = onode_getpar_flt_dflt(mo,"power",1.0);

  // Make the FT of the filters
  mod_me_fpe_s = matrix(1,xn,1,2*yn);
  contort_real_3d_farray(mod_me_fpe,1,1,1,xn,yn,tn);
  rlft3(mod_me_fpe,mod_me_fpe_s,xn,yn,tn,1);
}
/**************************************-**************************************/
/*                                                                           */
/*                           MODEL_GABOR_COMP_PREP                           */
/*                                                                           */
/*****************************************************************************/
void model_gabor_comp_prep(m)
     struct model_struct *m; // Model params
{
  int xn,yn,tn,write_filt;
  float ***fpe,***fpo,sscale,tscale;
  struct param_pair_list *mppl; // Model parameter pair list
  struct mod_me_fstruct **fs;  // [fn] List of pointers to filters

  // WYETH ERROR - this was already done, commented out June 7, 2016
  //myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  mylog(mylogf,"  MODEL_GABOR_COMP_PREP\n");
  mylog(mylogf,"    Computing even and odd, pref and null, filters.\n");

  mod_me_prep_basic(m->o,&sscale,&tscale,&xn,&yn,&tn); // Get and set globals

  mppl = m->ppl;

  if (m->o != NULL)
    mod_me_modtype = onode_getpar_chr_exit(m->o,"mod_type");
  else
    mod_me_modtype = paramfile_get_char_param_or_exit(mppl,"mod_type");


  // Note, only the first two filters are returned, the second are untouched
  if (m->o != NULL){
    fs = mod_me_get_quad_opp_filters_o(m->o,&mod_me_fpe,&mod_me_fpo,&mod_me_fpx,
				       &mod_me_fne,&mod_me_fno,&mod_me_fnx,
				       &mod_me_fn);
    mod_me_ff = fs;
  }else{
    mod_me_get_quad_opp_filters(mppl,&mod_me_fpe,&mod_me_fpo,
				&mod_me_fne,&mod_me_fno);
    exit_error("WYETH---","filters need to be put in a list, ...ff");
  }

  // exit_error("WYETH abc","fix this - use filters from list...");

  fpe = mod_me_fpe;  // fs[0]->f;
  fpo = mod_me_fpo;  // fs[1]->f;

  mod_me_normalize_filters(1,xn,1,yn,1,tn,sscale,tscale,1.0,2);

  if (m->o != NULL)
    write_filt = onode_getpar_int_dflt(m->o,"model_write_filter",-1);
  else
    write_filt = paramfile_get_int_param_default(mppl,"model_write_filter",-1);
  if ((write_filt == 0)||(write_filt == 4))
    write_3d_data_part("zzz.fe.3d",fpe,1,xn,1,yn,1,tn,4,2,1);
  if ((write_filt == 1)||(write_filt == 4))
    write_3d_data_part("zzz.fo.3d",fpo,1,xn,1,yn,1,tn,4,2,1);


  // Added delay, seconds  *** WYETH - is this needed, given spike_gen dt???
  if (m->o != NULL)
    mod_me_tdelay = onode_getpar_flt_dflt(m->o,"tdelay",0.0);
  else
    mod_me_tdelay = paramfile_get_float_param_default(mppl,"tdelay",0.0);



  if (m->o != NULL)
    mod_me_pow = onode_getpar_flt_dflt(m->o,"power",1.0);
  else
    mod_me_pow = paramfile_get_float_param_default(mppl,"power",1.0);

  // Make the FT of the filters
  mod_me_fpe_s = matrix(1,xn,1,2*yn);
  mod_me_fpo_s = matrix(1,xn,1,2*yn);
  contort_real_3d_farray(mod_me_fpe,1,1,1,xn,yn,tn);
  contort_real_3d_farray(mod_me_fpo,1,1,1,xn,yn,tn);
  rlft3(mod_me_fpe,mod_me_fpe_s,xn,yn,tn,1);
  rlft3(mod_me_fpo,mod_me_fpo_s,xn,yn,tn,1);
}
/**************************************-**************************************/
/*                                                                           */
/*                                MODEL_ME_FT                                */
/*                                                                           */
/*  Prepare the FT of the filter.                                            */
/*                                                                           */
/*****************************************************************************/
void model_me_ft(f0,f1,f2,f3,xn,yn,tn,rf0s,rf1s,rf2s,rf3s)
     float ***f0,***f1,***f2,***f3;
     int xn,yn,tn;
     float ***rf0s,***rf1s,***rf2s,***rf3s;
{
  float **f0_s,**f1_s,**f2_s,**f3_s;
  
  mylog(mylogf,"  MODEL_ME_FT\n");

  f0_s = matrix(1,xn,1,2*yn);
  f1_s = matrix(1,xn,1,2*yn);
  f2_s = matrix(1,xn,1,2*yn);
  f3_s = matrix(1,xn,1,2*yn);
  contort_real_3d_farray(f0,1,1,1,xn,yn,tn);
  contort_real_3d_farray(f1,1,1,1,xn,yn,tn);
  contort_real_3d_farray(f2,1,1,1,xn,yn,tn);
  contort_real_3d_farray(f3,1,1,1,xn,yn,tn);
  rlft3(f0,f0_s,xn,yn,tn,1);
  rlft3(f1,f1_s,xn,yn,tn,1);
  rlft3(f2,f2_s,xn,yn,tn,1);
  rlft3(f3,f3_s,xn,yn,tn,1);

  *rf0s = f0_s;
  *rf1s = f1_s;
  *rf2s = f2_s;
  *rf3s = f3_s;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MODEL_ME_01_PREP_FT                            */
/*                                                                           */
/*  NOTES                                                                    */
/*  - filters should be normalized before this is called.                    */
/*                                                                           */
/*****************************************************************************/
void model_me_01_prep_ft()
{
  int i;
  int xn,yn,tn,fn;
  struct mod_me_fstruct **fs;

  mylog(mylogf,"  MODEL_ME_01_PREP_FT\n");

  xn = mod_me_xn;
  yn = mod_me_yn;
  tn = mod_me_tn;

  if (mod_me_ftop == NULL){
    fs = mod_me_ff;         // OLD WAY - REMOVE when not used
    fn = mod_me_fn;
  }else{
    fs = mod_me_ftop->fsl;   // New way
    fn = mod_me_ftop->n;
  }

  for(i=0;i<fn;i++){
    //printf("fs___==============>  i=%d  xn,yn %d %d\n",i,xn,yn);
    fs[i]->s = matrix(1,xn,1,2*yn);
    contort_real_3d_farray(fs[i]->f,1,1,1,xn,yn,tn);
    rlft3(fs[i]->f,fs[i]->s,xn,yn,tn,1);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MODEL_MEM_PREP                               */
/*                                                                           */
/*****************************************************************************/
void model_mem_prep(m)
     struct model_struct *m; // Model params
{
  int i;
  int xn,yn,tn,write_filt,n,rtype,fn;
  float ***fpe,***fpo,***fne,***fno,sscale,tscale,norm_const;
  char ggstr[LONG_SLEN],tstr[SLEN];
  struct param_pair_list *mppl; /* Model parameter pair list. */
  struct mod_mem *mem;
  struct mod_me_fstruct **fs;  // [fn] List of pointers to filters
  struct onode *fo;

  mylog(mylogf,"  MODEL_MEM_PREP\n");
  mylog(mylogf,"    Building filters.\n");

  mod_me_prep_basic(m->o,&sscale,&tscale,&xn,&yn,&tn); // Get and set globals

  //
  //  Count the number of <filter> objects.
  //
  n =  onode_count_otype(m->o,"filter");
  sprintf(ggstr,"    Found %d filter definitions.\n",n);
  mylog(mylogf,ggstr);

  // Create struct to hold filters and filter parameters
  mem = (struct mod_mem *)myalloc(sizeof(struct mod_mem));
  mem->nf = n;
  mem->theta = (float *)myalloc(n*sizeof(float));
  mem->ssd = (float *)myalloc(n*sizeof(float));
  mem->sf = (float *)myalloc(n*sizeof(float));

  mem->tsd = (float *)myalloc(n*sizeof(float));
  mem->tf = (float *)myalloc(n*sizeof(float));
  mem->scale = (float *)myalloc(n*sizeof(float));
  mem->f = (float *****)myalloc(n*sizeof(float ****));
  mem->fs = (float ****)myalloc(n*sizeof(float ***));
  mem->r = (float *****)myalloc(n*sizeof(float ****));

  // Whether to dump the intermediate weighted responses
  mem->dumpflag = onode_getpar_int_dflt(m->o,"mem_dump_wr",0);


  //  Global filter list
  mod_me_fn = fn = n*4;
  mod_me_ff = fs = mod_me_fs_create(fn);

  fo = onode_get_next_type(m->o->o,"filter");
  for(i=0;i<n;i++){ // Build each filter
    mem->f[i] = (float ****)myalloc(4*sizeof(float ***));
    mem->fs[i] = (float ***)myalloc(4*sizeof(float **));
    mem->r[i] = (float ****)myalloc(4*sizeof(float ***));
    mem->r[i][0] = NULL;
    mem->r[i][1] = NULL;
    mem->r[i][2] = NULL;
    mem->r[i][3] = NULL;

    // *********  WYETH
    // *********  WYETH
    // *********  WYETH  Ultimately, do away with the 'mem' filter struct,
    // *********  WYETH    and use the 'fs' instead.
    // *********  WYETH
    // *********  WYETH

    mod_me_get_quad_opp_k(m->o,fo,&(mem->f[i][0]),&(mem->f[i][1]),
			  &(mem->f[i][2]),&(mem->f[i][3]),&norm_const);

    fs[i*4  ]->f = mem->f[i][0];  // Save in global list
    fs[i*4+1]->f = mem->f[i][1];
    fs[i*4+2]->f = mem->f[i][2];
    fs[i*4+3]->f = mem->f[i][3];

    fs[i*4  ]->normc = norm_const;
    fs[i*4+1]->normc = norm_const;
    fs[i*4+2]->normc = norm_const;
    fs[i*4+3]->normc = norm_const;

    mem->scale[i] = onode_getpar_flt_exit(fo,"ampl");

    fo = onode_get_next_type(fo->next,"filter");

    /***
	sprintf(tstr,"zzz.f%d_0.3d",i);
	write_3d_data_part(tstr,mem->f[i][0],1,xn,1,yn,1,tn,4,2,1);
	sprintf(tstr,"zzz.f%d_1.3d",i);
	write_3d_data_part(tstr,mem->f[i][1],1,xn,1,yn,1,tn,4,2,1);
	sprintf(tstr,"zzz.f%d_2.3d",i);
	write_3d_data_part(tstr,mem->f[i][2],1,xn,1,yn,1,tn,4,2,1);
	sprintf(tstr,"zzz.f%d_3.3d",i);
	write_3d_data_part(tstr,mem->f[i][3],1,xn,1,yn,1,tn,4,2,1);***/
  }

  mod_me_normalize_filters(1,xn,1,yn,1,tn,sscale,tscale,-1.0,2);

  mem->wf = onode_getpar_flt_exit(m->o,"mem_w_f");


  //
  // Space for ME response traces
  //
  rtype = onode_getpar_int_exit(m->o,"mem_rtype");

  if (rtype == 1){  // Single response at center
    mem->nr = 1;             // Number of responses
    mem->xi = (int *)myalloc(mem->nr * sizeof(int));
    mem->yi = (int *)myalloc(mem->nr * sizeof(int));
    mem->xi[0] = xn/2;
    mem->yi[0] = yn/2;
  }else if (rtype == 2){
    //
    //  Multiple responses in a line, should be orthogonal to filter
    //
    mem->nr = 10;             // Number of responses
    mem->xi = (int *)myalloc(mem->nr * sizeof(int));
    mem->yi = (int *)myalloc(mem->nr * sizeof(int));
    for(i=0;i<mem->nr;i++){
      mem->xi[i] = xn/2 - mem->nr/2 + i;
      mem->yi[i] = yn/2;
      //printf("    response %d at coords:   %d %d\n",i,mem->xi[i],mem->yi[i]);
    }
  }
  mem->mp = get_3d_farray(n,mem->nr,tn);
  mem->ma = get_3d_farray(n,mem->nr,tn);

  mod_me_mem = mem; // Global for multiple filters
}
/**************************************-**************************************/
/*                                                                           */
/*                              MODEL_REICH_PREP                             */
/*                                                                           */
/*****************************************************************************/
void model_reich_prep(m)
     struct model_struct *m; // Model params
{
  int xn,yn,tn,write_filt,xc,yc,tc;
  float hp_tau,lp_tau,sscale,tscale;
  struct onode *rdo;

  mylog(mylogf,"  MODEL_REICH_PREP\n");

  mod_me_modtype = onode_getpar_chr_exit(m->o,"mod_type");

  sscale = onode_getpar_flt_exit(m->o,"sscale");
  tscale = onode_getpar_flt_exit(m->o,"tscale");
  xn = onode_getpar_int_exit(m->o,"xn");
  yn = onode_getpar_int_exit(m->o,"yn");
  tn = onode_getpar_int_exit(m->o,"tn");

  mod_me_sscale = sscale;
  mod_me_tscale = tscale;
  mod_me_xn = xn;
  mod_me_yn = yn;
  mod_me_tn = tn;

  rdo = onode_child_get_unique(m->o,"filter");

  // Build filters, normalize
  lp_tau = onode_getpar_flt_exit(rdo,"lp_tau");
  hp_tau = onode_getpar_flt_exit(rdo,"hp_tau");
  mod_reich_lp = exp_time_tensor(mylogf,xn,yn,tn,sscale,tscale,lp_tau);
  mod_reich_hp = exp_time_tensor(mylogf,xn,yn,tn,sscale,tscale,hp_tau);
  //mod_reich_hp = delta_tensor(mylogf,xn,yn,tn,sscale,tscale);



  // HP = 1 - decay exp
  xc = (xn-1)/2; // Center of filter
  yc = (yn-1)/2;
  tc = (tn-1)/2;
  multiply_3d_farray(mod_reich_hp,1,xn,1,yn,1,tn,-1.0);
  mod_reich_hp[xc+1][yc+1][tc+1] += 1.0;
  //mod_reich_hp[xc][yc][tc] += 1.0; THIS IS WRONG.

  // WYETH - these are commented out because we want to keep area = 1
  //mod_me_normalize_filter(mod_reich_lp,1,xn,1,yn,1,tn,sscale,tscale,1.0);
  //mod_me_normalize_filter(mod_reich_hp,1,xn,1,yn,1,tn,sscale,tscale,1.0);


  //  TO DUMP the t-filters:
  //  append_farray_plot("zzz.tfile.pl","lp",mod_reich_lp[xc+1][yc+1],tn,1);
  //  append_farray_plot("zzz.tfile.pl","hp",mod_reich_hp[xc+1][yc+1],tn,1);


  write_filt = onode_getpar_int_dflt(rdo,"model_write_filter",-1);
  if (write_filt == 0)
    write_3d_data_part("zzz.lp.3d",mod_reich_lp,1,xn,1,yn,1,tn,4,2,1);
  if (write_filt == 1)
    write_3d_data_part("zzz.hp.3d",mod_reich_hp,1,xn,1,yn,1,tn,4,2,1);

  // FT of filters
  mod_reich_lp_s = matrix(1,xn,1,2*yn);
  mod_reich_hp_s = matrix(1,xn,1,2*yn);
  contort_real_3d_farray(mod_reich_lp,1,1,1,xn,yn,tn);
  contort_real_3d_farray(mod_reich_hp,1,1,1,xn,yn,tn);
  rlft3(mod_reich_lp,mod_reich_lp_s,xn,yn,tn,1);
  rlft3(mod_reich_hp,mod_reich_hp_s,xn,yn,tn,1);

  mod_me_tdelay = onode_getpar_flt_dflt(rdo,"tdelay",0.0);

  mod_me_oppflag = onode_getpar_int_dflt(rdo,"opp_flag",1);

  mod_reich_uxn = onode_getpar_int_exit(rdo,"n_unit_x");
  mod_reich_uyn = onode_getpar_int_dflt(rdo,"n_unit_y",1);
  mod_reich_dx  = onode_getpar_int_dflt(rdo,"dx",1);
  mod_reich_dy  = 0;

  mod_reich_ux0 = (xn - mod_reich_uxn) / 2;
  mod_reich_uy0 = (xn - mod_reich_uyn) / 2;

  // Check that the array of RD units falls within the grid
  if ((mod_reich_ux0 < 0) ||
      (mod_reich_uy0 < 0) ||
      (mod_reich_ux0 + mod_reich_uxn + mod_reich_dx > xn) ||
      (mod_reich_uy0 + mod_reich_uyn + mod_reich_dy > yn)){
    mylog_exit(mylogf,"MODEL_REICH_PREP  integration bounds error");
  }

  mod_reich_dump  = onode_getpar_int_dflt(rdo,"dump",0);
}
/**************************************-**************************************/
/*                                                                           */
/*                              MODEL_ME_01_DONE                             */
/*                                                                           */
/*****************************************************************************/
void model_me_01_done(mppl)
     struct param_pair_list *mppl;   // Model parameter pair list
{
  int xn,yn,tn;

  mylog(mylogf,"  MODEL_ME_01_DONE\n");
  
  xn = mod_me_xn;
  yn = mod_me_yn;
  tn = mod_me_tn;

  /*** WYETH - COMMENTED OUT OCT 2013 - because of crash - MUST FIX ***/
  /*** WYETH - COMMENTED OUT OCT 2013 - because of crash - MUST FIX ***/
  /*** WYETH - COMMENTED OUT OCT 2013 - because of crash - MUST FIX ***/
  /*
  free_f3tensor(mod_me_fpe,1,xn,1,yn,1,tn); mod_me_fpe = (float ***)NULL;
  free_f3tensor(mod_me_fpo,1,xn,1,yn,1,tn); mod_me_fpo = (float ***)NULL;
  free_f3tensor(mod_me_fne,1,xn,1,yn,1,tn); mod_me_fne = (float ***)NULL;
  free_f3tensor(mod_me_fno,1,xn,1,yn,1,tn); mod_me_fno = (float ***)NULL;
  */

  if (mod_me_fpe_s != NULL){ // If FT was used
    free_matrix(mod_me_fpe_s,1,xn,1,2*yn); mod_me_fpe_s = (float **)NULL;
    free_matrix(mod_me_fpo_s,1,xn,1,2*yn); mod_me_fpo_s = (float **)NULL;
    free_matrix(mod_me_fne_s,1,xn,1,2*yn); mod_me_fne_s = (float **)NULL;
    free_matrix(mod_me_fno_s,1,xn,1,2*yn); mod_me_fno_s = (float **)NULL;
  }

  // I separated these commands for 'ds_simp_01'
  if (mod_me_rpe != NULL){
    free_f3tensor(mod_me_rpe,1,xn,1,yn,1,tn);
  }
  if (mod_me_rpo != NULL)
    free_f3tensor(mod_me_rpo,1,xn,1,yn,1,tn);

  if (mod_me_rne != NULL){
    free_f3tensor(mod_me_rne,1,xn,1,yn,1,tn);
    free_f3tensor(mod_me_rno,1,xn,1,yn,1,tn);
  }

  if (mod_me_normflag > 0){
    free_f3tensor(mod_me_fg0,1,xn,1,yn,1,tn); mod_me_fg0 = (float ***)NULL;
    free_f3tensor(mod_me_fg1,1,xn,1,yn,1,tn); mod_me_fg1 = (float ***)NULL;
    if (mod_me_fg0_s != NULL){
      free_matrix(mod_me_fg0_s,1,xn,1,2*yn); mod_me_fg0_s = (float **)NULL;
      free_matrix(mod_me_fg1_s,1,xn,1,2*yn); mod_me_fg1_s = (float **)NULL;
    }
    if (mod_me_rg0 != NULL){
      free_f3tensor(mod_me_rg0,1,xn,1,yn,1,tn);
      free_f3tensor(mod_me_rg1,1,xn,1,yn,1,tn);
    }
  }

  myfree(mod_me_modtype);
}
/**************************************-**************************************/
/*                                                                           */
/*                       MOD_ME_COMPUTE_RESPONSE_FROM_FT                     */
/*                                                                           */
/*  The precomputed FT of the filters is passed in, in the filter struct.    */
/*                                                                           */
/*****************************************************************************/
void mod_me_compute_response_from_ft(data,fs,fn,xn,yn,tn)
     float ***data;               // [xn][yn][tn] Data to convolve with filters
     struct mod_me_fstruct **fs;  // [fn] List of pointers to filters
     int fn;                      // Number of filters
     int xn,yn,tn;                // Dimensions of data and filters
{
  int i;
  float fc,**speq,***r,**rs;

  mylog(mylogf,"  MOD_ME_COMPUTE_RESPONSE_FROM_FT\n");
  fc = 2.0/((float)(xn*yn*tn));

  mylog(mylogf,"    Stimulus FFT\n");

  speq = matrix(1,xn,1,2*yn);
  contort_real_3d_farray(data,1,1,1,xn,yn,tn);
  rlft3(data,speq,xn,yn,tn,1);

  if (myid == -1){
    printf("    Stim FT times ME FT.\n");
    printf("      filter:  ");
    fflush(stdout);
  }else{
    mylog(mylogf,"    Stim FT times ME FT\n");
    mylog(mylogf,"    Filter:  ");
  }

  for(i=0;i<fn;i++){  //  For each filter
    if (fs[i]->flag == 1){  // This is '1' if filter apples to current eye
      three_d_fft_prod(data,speq,fs[i]->f,fs[i]->s,xn,yn,tn,&r,&rs);
      rlft3(r,rs,xn,yn,tn,-1); free_matrix(rs,1,xn,1,2*yn);
      multiply_3d_farray(r,1,xn,1,yn,1,tn,fc);
      contort_real_3d_farray(r,1,1,1,xn,yn,tn);
      fs[i]->r = r;  // Store response in filter struct.
    
      if (myid == -1){
	printf("%d ",i);
	fflush(stdout);
      }else{
	sprintf(ggstr,"%d ",i);
	mylog(mylogf,ggstr);
      }
    }
  }

  if (myid == -1){
    printf("\n");
  }else{
    mylog(mylogf,"\n");
  }

  free_matrix(speq,1,xn,1,2*yn);
}
/**************************************-**************************************/
/*                                                                           */
/*                   MOD_ME_COMPUTE_RESPONSE_QUAD_OPP_FROM_FT                */
/*                                                                           */
/*  The precomputed FT of the filters is passed.                             */
/*                                                                           */
/*****************************************************************************/
void mod_me_compute_response_quad_opp_from_ft(data,fpe,fpo,fne,fno,
					      fpe_s,fpo_s,fne_s,fno_s,
					      xn,yn,tn,rrpe,rrpo,rrne,rrno,
					      rrg0,rrg1)
     float ***data,***fpe,***fpo,***fne,***fno;
     float **fpe_s,**fpo_s,**fne_s,**fno_s;
     int xn,yn,tn;
     float ****rrpe,****rrpo,****rrne,****rrno;
     float ****rrg0,****rrg1;
{
  int i;
  int fn;
  struct mod_me_fstruct **fs;

  mylog(mylogf,"  MOD_ME_COMPUTE_RESPONSE_QUAD_OPP_FROM_FT\n");

  //
  //  Create and load filter structure 'fs'
  //
  fn = 4;
  if (mod_me_normflag > 0)
    fn += 2;
  fs = (struct mod_me_fstruct **)myalloc(fn*sizeof(struct mod_me_fstruct *));
  for(i=0;i<fn;i++)
    fs[i] = (struct mod_me_fstruct *)myalloc(sizeof(struct mod_me_fstruct));
  fs[0]->f = fpe;  fs[0]->s = fpe_s;
  fs[1]->f = fpo;  fs[1]->s = fpo_s;
  fs[2]->f = fne;  fs[2]->s = fne_s;
  fs[3]->f = fno;  fs[3]->s = fno_s;
  if (mod_me_normflag > 0){
    fs[4]->f = mod_me_fg0;  fs[4]->s = mod_me_fg0_s;
    fs[5]->f = mod_me_fg1;  fs[5]->s = mod_me_fg1_s;
  }


  //
  //  Get response for each filter
  //
  mod_me_compute_response_from_ft(data,fs,fn,xn,yn,tn);


  //
  //  Store results, and discard temporary 'fs' structure
  //
  *rrpe = fs[0]->r;
  *rrpo = fs[1]->r;
  *rrne = fs[2]->r;
  *rrno = fs[3]->r;
  if (mod_me_normflag > 0){
    *rrg0 = fs[4]->r;
    *rrg1 = fs[5]->r;
  }else{
    *rrg0 = *rrg1 = NULL;
  }
  for(i=0;i<fn;i++)
    myfree(fs[i]);
  myfree(fs);
}
/**************************************-**************************************/
/*                                                                           */
/*                   MOD_ME_COMPUTE_RESPONSE_QUAD_OPP_FROM_FT                */
/*                                                                           */
/*  The precomputed FT of the filters is passed.                             */
/*                                                                           */
/*****************************************************************************/
void mod_me_compute_response_quad_opp_from_ft_x(data,fpe,fpo,fpx,fne,fno,fnx,
						fpe_s,fpo_s,fpx_s,
						fne_s,fno_s,fnx_s,
						xn,yn,tn,rrpe,rrpo,rrpx,
						rrne,rrno,rrnx,
						rrg0,rrg1)
     float ***data,***fpe,***fpo,***fpx,***fne,***fno,***fnx;
     float **fpe_s,**fpo_s,**fpx_s,**fne_s,**fno_s,**fnx_s;
     int xn,yn,tn;
     float ****rrpe,****rrpo,****rrpx,****rrne,****rrno,****rrnx;
     float ****rrg0,****rrg1;
{
  float **speq,***rpe,***rpo,***rpx,***rne,***rno,***rnx;
  float **pe_s,**po_s,**px_s,**ne_s,**no_s,**nx_s,fc;

  float ***rg0,***rg1;
  float **g0_s,**g1_s;


  mylog(mylogf,"  MOD_ME_COMPUTE_RESPONSE_QUAD_OPP_FROM_FT_X\n");
  fc = 2.0/((float)(xn*yn*tn));

  if (myid == -1){
    printf("      Stimulus FFT.");
    fflush(stdout);
  }else
    mylog(mylogf,"    Stimulus FFT\n");

  speq = matrix(1,xn,1,2*yn);
  contort_real_3d_farray(data,1,1,1,xn,yn,tn);
  rlft3(data,speq,xn,yn,tn,1);

  if (myid == -1){
    printf("  Stim FT times ME FT.");
    fflush(stdout);
  }else
    mylog(mylogf,"    Stim FT times ME FT\n");

  three_d_fft_prod(data,speq,fpe,fpe_s,xn,yn,tn,&rpe,&pe_s);
  three_d_fft_prod(data,speq,fpo,fpo_s,xn,yn,tn,&rpo,&po_s);
  three_d_fft_prod(data,speq,fpx,fpx_s,xn,yn,tn,&rpx,&px_s);
  three_d_fft_prod(data,speq,fne,fne_s,xn,yn,tn,&rne,&ne_s);
  three_d_fft_prod(data,speq,fno,fno_s,xn,yn,tn,&rno,&no_s);
  three_d_fft_prod(data,speq,fnx,fnx_s,xn,yn,tn,&rnx,&nx_s);

  if (mod_me_normflag > 0){
    // NOTE, the gain filters are global, instead of passed in
    three_d_fft_prod(data,speq,mod_me_fg0,mod_me_fg0_s,xn,yn,tn,&rg0,&g0_s);
    three_d_fft_prod(data,speq,mod_me_fg1,mod_me_fg1_s,xn,yn,tn,&rg1,&g1_s);
  }

  free_matrix(speq,1,xn,1,2*yn);

  if (myid == -1){
    printf("  Inverse FFT.");
    fflush(stdout);
  }else
    mylog(mylogf,"    Inverse FFT\n");

  // WYETH - DEBUG MUSTAFA
  sprintf(ggstr,"   xn,yn,tn  %d %d %d",xn,yn,tn);
  mylog(mylogf,ggstr);

  rlft3(rpe,pe_s,xn,yn,tn,-1); free_matrix(pe_s,1,xn,1,2*yn);
  rlft3(rpo,po_s,xn,yn,tn,-1); free_matrix(po_s,1,xn,1,2*yn);
  rlft3(rpx,px_s,xn,yn,tn,-1); free_matrix(px_s,1,xn,1,2*yn);
  rlft3(rne,ne_s,xn,yn,tn,-1); free_matrix(ne_s,1,xn,1,2*yn);
  rlft3(rno,no_s,xn,yn,tn,-1); free_matrix(no_s,1,xn,1,2*yn);
  rlft3(rnx,nx_s,xn,yn,tn,-1); free_matrix(nx_s,1,xn,1,2*yn);
  multiply_3d_farray(rpe,1,xn,1,yn,1,tn,fc);
  multiply_3d_farray(rpo,1,xn,1,yn,1,tn,fc);
  multiply_3d_farray(rpx,1,xn,1,yn,1,tn,fc);
  multiply_3d_farray(rne,1,xn,1,yn,1,tn,fc);
  multiply_3d_farray(rno,1,xn,1,yn,1,tn,fc);
  multiply_3d_farray(rnx,1,xn,1,yn,1,tn,fc);
  contort_real_3d_farray(rpe,1,1,1,xn,yn,tn);
  contort_real_3d_farray(rpo,1,1,1,xn,yn,tn);
  contort_real_3d_farray(rpx,1,1,1,xn,yn,tn);
  contort_real_3d_farray(rne,1,1,1,xn,yn,tn);
  contort_real_3d_farray(rno,1,1,1,xn,yn,tn);
  contort_real_3d_farray(rnx,1,1,1,xn,yn,tn);

  if (mod_me_normflag > 0){
    mylog(mylogf,"    'normflag' Computation\n");

    rlft3(rg0,g0_s,xn,yn,tn,-1); free_matrix(g0_s,1,xn,1,2*yn);
    rlft3(rg1,g1_s,xn,yn,tn,-1); free_matrix(g1_s,1,xn,1,2*yn);
    multiply_3d_farray(rg0,1,xn,1,yn,1,tn,fc);
    multiply_3d_farray(rg1,1,xn,1,yn,1,tn,fc);
    contort_real_3d_farray(rg0,1,1,1,xn,yn,tn);
    contort_real_3d_farray(rg1,1,1,1,xn,yn,tn);
  }else{
    rg0 = rg1 = NULL;
  }

  if (myid == -1)
    printf("  Done.\n");
  else
    mylog(mylogf,"    Done.\n");

  *rrpe = rpe;
  *rrpo = rpo;
  *rrpx = rpx;
  *rrne = rne;
  *rrno = rno;
  *rrnx = rnx;

  *rrg0 = rg0;
  *rrg1 = rg1;
}
/**************************************-**************************************/
/*                                                                           */
/*                 MOD_ME_COMPUTE_RESPONSE_GABOR_COMP_FROM_FT                */
/*                                                                           */
/*  The precomputed FT of the filters is passed.                             */
/*                                                                           */
/*****************************************************************************/
void mod_me_compute_response_gabor_comp_from_ft(data,fpe,fpo,fpe_s,fpo_s,
						xn,yn,tn,rrpe,rrpo)
     float ***data,***fpe,***fpo;
     float **fpe_s,**fpo_s;
     int xn,yn,tn;
     float ****rrpe,****rrpo;
{
  float **speq,***rpe,***rpo;
  float **pe_s,**po_s,fc;

  mylog(mylogf,"  MOD_ME_COMPUTE_RESPONSE_GABOR_COMP_FROM_FT\n");
  fc = 2.0/((float)(xn*yn*tn));

  if (myid == -1){
    printf("      Stimulus FFT.");
    fflush(stdout);
  }

  speq = matrix(1,xn,1,2*yn);
  contort_real_3d_farray(data,1,1,1,xn,yn,tn);
  rlft3(data,speq,xn,yn,tn,1);

  if (myid == -1){
    printf("  Stim FT times ME FT.");
    fflush(stdout);
  }
  three_d_fft_prod(data,speq,fpe,fpe_s,xn,yn,tn,&rpe,&pe_s);
  three_d_fft_prod(data,speq,fpo,fpo_s,xn,yn,tn,&rpo,&po_s);

  free_matrix(speq,1,xn,1,2*yn);

  if (myid == -1){
    printf("  Inverse FFT.");
    fflush(stdout);
  }
  rlft3(rpe,pe_s,xn,yn,tn,-1); free_matrix(pe_s,1,xn,1,2*yn);
  rlft3(rpo,po_s,xn,yn,tn,-1); free_matrix(po_s,1,xn,1,2*yn);

  multiply_3d_farray(rpe,1,xn,1,yn,1,tn,fc);
  multiply_3d_farray(rpo,1,xn,1,yn,1,tn,fc);

  contort_real_3d_farray(rpe,1,1,1,xn,yn,tn);
  contort_real_3d_farray(rpo,1,1,1,xn,yn,tn);

  if (myid == -1)
    printf("  Done.\n");

  *rrpe = rpe;
  *rrpo = rpo;
}
/**************************************-**************************************/
/*                                                                           */
/*                      MOD_REICH_COMPUTE_RESPONSE_FROM_FT                   */
/*                                                                           */
/*  The precomputed FT of the filters is passed.                             */
/*                                                                           */
/*****************************************************************************/
void mod_reich_compute_response_from_ft(data,flp,fhp,flp_s,fhp_s,
					xn,yn,tn,rrhp,rrlp)
     float ***data,***flp,***fhp;
     float **flp_s,**fhp_s;
     int xn,yn,tn;
     float ****rrhp,****rrlp;
{
  float **speq,***rlp,***rhp;
  float **lp_s,**hp_s,fc;

  mylog(mylogf,"  MOD_REICH_COMPUTE_RESPONSE_FROM_FT\n");
  fc = 2.0/((float)(xn*yn*tn));

  if (myid == -1){
    printf("      Stimulus FFT.");
    fflush(stdout);
  }

  speq = matrix(1,xn,1,2*yn);
  contort_real_3d_farray(data,1,1,1,xn,yn,tn);
  rlft3(data,speq,xn,yn,tn,1);

  if (myid == -1){
    printf("  Stim FT times ME FT.");
    fflush(stdout);
  }
  three_d_fft_prod(data,speq,flp,flp_s,xn,yn,tn,&rlp,&lp_s);
  three_d_fft_prod(data,speq,fhp,fhp_s,xn,yn,tn,&rhp,&hp_s);

  free_matrix(speq,1,xn,1,2*yn);

  if (myid == -1){
    printf("  Inverse FFT.");
    fflush(stdout);
  }
  rlft3(rlp,lp_s,xn,yn,tn,-1); free_matrix(lp_s,1,xn,1,2*yn);
  rlft3(rhp,hp_s,xn,yn,tn,-1); free_matrix(hp_s,1,xn,1,2*yn);

  multiply_3d_farray(rlp,1,xn,1,yn,1,tn,fc);
  multiply_3d_farray(rhp,1,xn,1,yn,1,tn,fc);

  contort_real_3d_farray(rlp,1,1,1,xn,yn,tn);
  contort_real_3d_farray(rhp,1,1,1,xn,yn,tn);

  if (myid == -1)
    printf("  Done.\n");

  *rrlp = rlp;
  *rrhp = rhp;
}
/**************************************-**************************************/
/*                                                                           */
/*                        MOD_ME_COMPUTE_RESPONSE_MEM                        */
/*                                                                           */
/*  The precomputed FT of the filters is passed.                             */
/*                                                                           */
/*****************************************************************************/
void mod_me_compute_response_mem(data,mem,xn,yn,tn)
     float ***data;
     struct mod_mem *mem;
     int xn,yn,tn;
{
  int i;
  int n;
  float **speq,**pe_s,**po_s,**ne_s,**no_s,fc;

  mylog(mylogf,"  MOD_ME_COMPUTE_RESPONSE_MEM\n");
  fc = 2.0/((float)(xn*yn*tn));

  if (myid == -1){
    printf("      Stimulus FFT.");
    fflush(stdout);
  }
  speq = matrix(1,xn,1,2*yn);
  contort_real_3d_farray(data,1,1,1,xn,yn,tn);
  rlft3(data,speq,xn,yn,tn,1);

  if (myid == -1){
    printf("  Stim FT times ME FT.");
    fflush(stdout);
  }

  n = mem->nf;
  for(i=0;i<n;i++){

    three_d_fft_prod(data,speq,mem->f[i][0],mem->fs[i][0],xn,yn,tn,
		     &(mem->r[i][0]),&pe_s);
    three_d_fft_prod(data,speq,mem->f[i][1],mem->fs[i][1],xn,yn,tn,
		     &(mem->r[i][1]),&po_s);
    three_d_fft_prod(data,speq,mem->f[i][2],mem->fs[i][2],xn,yn,tn,
		     &(mem->r[i][2]),&ne_s);
    three_d_fft_prod(data,speq,mem->f[i][3],mem->fs[i][3],xn,yn,tn,
		     &(mem->r[i][3]),&no_s);

    rlft3(mem->r[i][0],pe_s,xn,yn,tn,-1); free_matrix(pe_s,1,xn,1,2*yn);
    rlft3(mem->r[i][1],po_s,xn,yn,tn,-1); free_matrix(po_s,1,xn,1,2*yn);
    rlft3(mem->r[i][2],ne_s,xn,yn,tn,-1); free_matrix(ne_s,1,xn,1,2*yn);
    rlft3(mem->r[i][3],no_s,xn,yn,tn,-1); free_matrix(no_s,1,xn,1,2*yn);
    multiply_3d_farray(mem->r[i][0],1,xn,1,yn,1,tn,fc);
    multiply_3d_farray(mem->r[i][1],1,xn,1,yn,1,tn,fc);
    multiply_3d_farray(mem->r[i][2],1,xn,1,yn,1,tn,fc);
    multiply_3d_farray(mem->r[i][3],1,xn,1,yn,1,tn,fc);

    contort_real_3d_farray(mem->r[i][0],1,1,1,xn,yn,tn);
    contort_real_3d_farray(mem->r[i][1],1,1,1,xn,yn,tn);
    contort_real_3d_farray(mem->r[i][2],1,1,1,xn,yn,tn);
    contort_real_3d_farray(mem->r[i][3],1,1,1,xn,yn,tn);
    
    if (myid == -1)
      printf("    Done filter %d.\n",i);
  }

  free_matrix(speq,1,xn,1,2*yn);
}
/**************************************-**************************************/
/*                                                                           */
/*                     MOD_ME_COMPUTE_RESPONSE_SIMP_FROM_FT                  */
/*                                                                           */
/*  The precomputed FT of the filters is passed.                             */
/*                                                                           */
/*****************************************************************************/
void mod_me_compute_response_simp_from_ft(data,fpe,fpe_s,xn,yn,tn,rrpe)
     float ***data,***fpe;
     float **fpe_s;
     int xn,yn,tn;
     float ****rrpe;
{
  float **speq,***rpe;
  float **pe_s,fc;

  mylog(mylogf,"  MOD_ME_COMPUTE_RESPONSE_SIMP_FROM_FT\n");
  fc = 2.0/((float)(xn*yn*tn));

  if (myid == -1){
    printf("      Stimulus FFT.");
    fflush(stdout);
  }
  speq = matrix(1,xn,1,2*yn);
  contort_real_3d_farray(data,1,1,1,xn,yn,tn);

  rlft3(data,speq,xn,yn,tn,1);

  if (myid == -1){
    printf("  Stim FT times ME FT.");
    fflush(stdout);
  }

  three_d_fft_prod(data,speq,fpe,fpe_s,xn,yn,tn,&rpe,&pe_s);

  free_matrix(speq,1,xn,1,2*yn);

  if (myid == -1){
    printf("  Inverse FFT.");
    fflush(stdout);
  }
  rlft3(rpe,pe_s,xn,yn,tn,-1); free_matrix(pe_s,1,xn,1,2*yn);
  multiply_3d_farray(rpe,1,xn,1,yn,1,tn,fc);
  contort_real_3d_farray(rpe,1,1,1,xn,yn,tn);

  if (myid == -1)
    printf("  Done.\n");

  *rrpe = rpe;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MODEL_ME_HANDOVER_BASIC                         */
/*                                                                           */
/*****************************************************************************/
void mod_me_handover_basic(m,r,mep_raw,men_raw,mep,men,rpe,rpo,rne,rno,
			   xn,yn,tn)
     struct model_struct *m;      // Model params
     struct response_struct *r;   // Response data
     float *mep_raw,*men_raw;     // ME pref,anti, before power/gain/etc [tn]
     float *mep,*men;             // ME signals, after power/gain [tn]
     float ***rpe,***rpo,***rne,***rno; // Raw P/A even/odd [xn][yn][1..tn]
     int xn,yn,tn;                // signal duration
{
  int i;
  float *t;

  for(i=0;i<r->n;i++){ // For each response requested

    //printf("=====> %s\n",r->datid[i]);

    if (strcmp(r->datid[i],"me_pref")==0){
      mod_util_resp_store_f_samp(r,i,mep,tn,1.0/mod_me_tscale,mylogf);
    }else if (strcmp(r->datid[i],"me_anti")==0){
      mod_util_resp_store_f_samp(r,i,men,tn,1.0/mod_me_tscale,mylogf);
    }else if (strcmp(r->datid[i],"me_opp")==0){
      if ((mep != NULL) && (men != NULL)){
	t = diff_farrays(mep,men,tn); // Opponent ME
	mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
	myfree(t);
      }
    }else if (strcmp(r->datid[i],"me_pref_raw")==0){
      mod_util_resp_store_f_samp(r,i,mep_raw,tn,1.0/mod_me_tscale,mylogf);
    }else if (strcmp(r->datid[i],"me_anti_raw")==0){
      mod_util_resp_store_f_samp(r,i,men_raw,tn,1.0/mod_me_tscale,mylogf);
    }else if (strcmp(r->datid[i],"filter_pref_even")==0){
      if (rpe != NULL){
	t = &(rpe[xn/2][yn/2][1]);
	mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
      }
    }else if (strcmp(r->datid[i],"filter_pref_odd")==0){
      if (rpo != NULL){
	t = &(rpo[xn/2][yn/2][1]);
	mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
      }
    }else if (strcmp(r->datid[i],"filter_anti_even")==0){
      if (rne != NULL){
	t = &(rne[xn/2][yn/2][1]);
	mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
      }
    }else if (strcmp(r->datid[i],"filter_anti_odd")==0){
      if (rno != NULL){
	t = &(rno[xn/2][yn/2][1]);
	mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
      }
    }else if (strcmp(r->datid[i],"filter_out")==0){
      if (strcmp(mod_me_modtype,"ds_simp_01")==0){
	// This is the filter output (possibly scaled)
	mod_util_resp_store_f_samp(r,i,mep,tn,1.0/mod_me_tscale,mylogf);
      }else if (strcmp(mod_me_modtype,"ds_liv3")==0){
	mod_util_resp_store_f_samp(r,i,mep_raw,tn,1.0/mod_me_tscale,mylogf);
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_ME_RESP_BINOC_THRESH                         */
/*                                                                           */
/*  Apply any threshold, sign inversion and rectification.                   */
/*                                                                           */
/*****************************************************************************/
float *mod_me_resp_binoc_thresh(float *d,       // [n] response data
				int n,          // length of data
				float thresh,   // threshold
				int signflag){  // -1 neg, 0-full, 1-pos
  int i;
  float t,*y,fsgn;

  fsgn = 1.0;
  if (signflag == -1)
    fsgn = -1.0;

  y = (float *)myalloc(n*sizeof(float));

  if (mod_me_binoc_srect == 1){
    for(i=0;i<n;i++){
      t = fsgn*d[i] - thresh;
      if (t >= 0.0)
	y[i] = t;
      else
	y[i] = 0.0;
    }
  }else{
    for(i=0;i<n;i++)
      y[i] = fsgn*d[i] - thresh;
  }

  return y;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_ME_RESP_BINOC_SIG                          */
/*                                                                           */
/*  Return one of the 4 binocular simple signals in the BDE model.           */
/*                                                                           */
/*****************************************************************************/
float *mod_me_resp_binoc_sig(int xi, int yi,  // Coord w/i spatial field
			     char *sid,       // Signal ID:  b1p, b1n, b2p, b2n
			     float thresh,    // Threshold value
			     int sqflag){     // 0-none, 1-halfsquare, 2-square
  int i;
  int tn,signflag,si;
  float ty,sr,*tl,*tr,*ttl,*ttr,*y;

  tn = mod_me_tn;

  //
  //  Get the raw left and right eye monocular siganls, 'tl' and 'tr'
  //
  if ((strcmp(sid,"b1p")==0) || (strcmp(sid,"b1n")==0)){
    tl = &(mod_me_ff[0]->r[xi][yi][1]);   // Left even
    tr = &(mod_me_ff[2]->r[xi][yi][1]);   // Right even

    if (strcmp(sid,"b1p")==0)  // Sign bit in ..._rsign
      si = 0;
    else
      si = 1;

  }else if ((strcmp(sid,"b2p")==0) || (strcmp(sid,"b2n")==0)){
    tl = &(mod_me_ff[1]->r[xi][yi][1]);   // Left odd
    tr = &(mod_me_ff[3]->r[xi][yi][1]);   // Right odd

    if (strcmp(sid,"b2p")==0)  // Sign bit in ..._rsign
      si = 2;
    else
      si = 3;

  }else{
    exit_error("MOD_ME_RESP_BINOC_SIG","Unknown sid");
  }

  if (sid[2] == 'p')
    signflag = 1;
  else if (sid[2] == 'n')
    signflag = -1;

  if (mod_me_binoc_rsign[si] == '+')
    sr = 1.0;
  else
    sr = -1.0;

  ttl = mod_me_resp_binoc_thresh(tl,tn,thresh,signflag);
  ttr = mod_me_resp_binoc_thresh(tr,tn,thresh,signflag);

  y = (float *)myalloc(tn*sizeof(float));
  for(i=0;i<tn;i++){
    if (sqflag == 1){  // Half-square (Half-wave rectify before square)
      ty = ttl[i] + sr*ttr[i];
      if (ty > 0.0)
	y[i] = ty*ty;
      else
	y[i] = 0.0;
    }else if (sqflag == 2)  // Full square
      y[i] = (ttl[i] + sr*ttr[i]) * (ttl[i] + sr*ttr[i]);
    else
      y[i] = ttl[i] + sr*ttr[i];  // No non-lin.
  }

  return y;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_ME_GEN_POISSON                            */
/*                                                                           */
/*****************************************************************************/
void mod_me_gen_poisson(m,r,sgo,d,n,seed,rfs,rns,rv,rvn)
     struct model_struct *m;      // Model params
     struct response_struct *r;   // Response data
     struct onode *sgo;           // "<spike_gen>" onode
     float *d;                    // Data trace [n]
     int n;                       // length of data
     int seed;
     float **rfs;                 // Return float spikes [rns]
     int    *rns;                 // Return number of spikes, times in ms
     float **rv;                  // Return float prob trace [rvn]
     int    *rvn;                 // Return length of prob. trace
{
  int ns,vn;
  float *fs,*v,samp,*t;

  samp = 1.0/mod_me_tscale;

  t = copy_farray(d,n);

  // 'vflag' = 1, so v and vn are returned
  // Spikes are returned in msec

  //WYETH REMOVE sgo = onode_child_get_unique(m->o,"spike_gen");
  ifc_util_poisson(mylogf,m,sgo,t,n,samp,seed,0,1,&fs,&ns,&v,&vn);

  sprintf(ggstr,"    Made %d Poisson spikes.\n",ns);
  mylog(mylogf,ggstr);

  // Assuming spikes are returned in milliseconds
  add_const_farray(fs,ns,mod_me_tdelay*1000.0); // Added latency

  *rfs = fs;
  *rns = ns;
  *rv  = v;
  *rvn = vn;

  myfree(t);  // WYETH BUG FIXED HERE Oct 20, 2016 //
  //myfree(d);  // WYETH BUG FIXED HERE Oct 20, 2016 //
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_ME_GEN_POISSON_DUMMY                        */
/*                                                                           */
/*  FOR DEBUGGING                                                            */
/*  Can remove.                                                              */
/*                                                                           */
/*****************************************************************************/
void mod_me_gen_poisson_dummy(m,r,sgo,d,n,seed,rfs,rns,rv,rvn)
     struct model_struct *m;      // Model params
     struct response_struct *r;   // Response data
     struct onode *sgo;           // "<spike_gen>" onode
     float *d;                    // Data trace [n]
     int n;                       // length of data
     int seed;
     float **rfs;                 // Return float spikes [rns]
     int    *rns;                 // Return number of spikes, times in ms
     float **rv;                  // Return float prob trace [rvn]
     int    *rvn;                 // Return length of prob. trace
{
  int ns,vn;
  float *fs,*v;

  printf("____________MOD_ME_GEN_POISSON_DUMMY_______________\n");

  ns = 2;
  fs = (float *)myalloc(ns*sizeof(float));
  fs[0] = 100.0;
  fs[1] = 100.0;

  vn = n;
  v = get_zero_farray(vn);
  
  *rfs = fs;
  *rns = ns;
  *rv  = v;
  *rvn = vn;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_ME_HANDOVER_BINOC                          */
/*                                                                           */
/*****************************************************************************/
void mod_me_handover_binoc(m,r,xi,yi,binen,xn,yn,tn)
     struct model_struct *m;      // Model params
     struct response_struct *r;   // Response data
     int xi,yi;                   // Coordinate from which to take response
     float *binen;                // Binocular energy signal, before scaling
     int xn,yn,tn;                // signal duration
{
  int i;
  int sqflag;
  float thresh,*t,*t1,*t2,*tpos,*tneg;

  thresh = mod_me_binoc_athresh;  // One threshold for all filters

  if (strcmp(mod_me_binoc_nonlin,"halfsq")==0)
    sqflag = 1;
  else if (strcmp(mod_me_binoc_nonlin,"square")==0)
    sqflag = 2;
  else if (strcmp(mod_me_binoc_nonlin,"none")==0)
    sqflag = 0;

  for(i=0;i<r->n;i++){ // For each response requested
    //printf("=====> %s\n",r->datid[i]);
    if (strcmp(r->datid[i],"binoc_energy")==0){
      mod_util_resp_store_f_samp(r,i,binen,tn,1.0/mod_me_tscale,mylogf);
    }else if (strcmp(r->datid[i],"filter_left_1")==0){  // E.g., even
      t = &(mod_me_ff[0]->r[xi][yi][1]);
      mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
    }else if (strcmp(r->datid[i],"filter_left_2")==0){  // E.g., odd
      t = &(mod_me_ff[1]->r[xi][yi][1]);
      mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
    }else if (strcmp(r->datid[i],"filter_right_1")==0){
      t = &(mod_me_ff[2]->r[xi][yi][1]);
      mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
    }else if (strcmp(r->datid[i],"filter_right_2")==0){
      t = &(mod_me_ff[3]->r[xi][yi][1]);
      mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
    }else if (strcmp(r->datid[i],"binoc_1_pos")==0){
      t = mod_me_resp_binoc_sig(xi,yi,"b1p",thresh,0);  // 0-no square
      mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
      myfree(t);
    }else if (strcmp(r->datid[i],"binoc_1_neg")==0){
      t = mod_me_resp_binoc_sig(xi,yi,"b1n",thresh,0);  // 0-no square
      mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
      myfree(t);
    }else if (strcmp(r->datid[i],"binoc_2_pos")==0){
      t = mod_me_resp_binoc_sig(xi,yi,"b2p",thresh,0);  // 0-no square
      mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
      myfree(t);
    }else if (strcmp(r->datid[i],"binoc_2_neg")==0){
      t = mod_me_resp_binoc_sig(xi,yi,"b2n",thresh,0);  // 0-no square
      mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
      myfree(t);
    }else if (strcmp(r->datid[i],"binoc_1_pos_sq")==0){
      t = mod_me_resp_binoc_sig(xi,yi,"b1p",thresh,sqflag);
      mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
      myfree(t);
    }else if (strcmp(r->datid[i],"binoc_1_neg_sq")==0){
      t = mod_me_resp_binoc_sig(xi,yi,"b1n",thresh,sqflag);
      mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
      myfree(t);
    }else if (strcmp(r->datid[i],"binoc_2_pos_sq")==0){
      t = mod_me_resp_binoc_sig(xi,yi,"b2p",thresh,sqflag);
      mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
      myfree(t);
    }else if (strcmp(r->datid[i],"binoc_2_neg_sq")==0){
      t = mod_me_resp_binoc_sig(xi,yi,"b2n",thresh,sqflag);
      mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
      myfree(t);
    }else if (strcmp(r->datid[i],"binoc_1")==0){
      tpos = mod_me_resp_binoc_sig(xi,yi,"b1p",thresh,sqflag);
      tneg = mod_me_resp_binoc_sig(xi,yi,"b1n",thresh,sqflag);
      t = add_farrays(tpos,tneg,tn);
      mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
      myfree(tpos); myfree(tneg); myfree(t);
    }else if (strcmp(r->datid[i],"binoc_2")==0){
      tpos = mod_me_resp_binoc_sig(xi,yi,"b2p",thresh,sqflag);
      tneg = mod_me_resp_binoc_sig(xi,yi,"b2n",thresh,sqflag);
      t = add_farrays(tpos,tneg,tn);
      mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
      myfree(tpos); myfree(tneg); myfree(t);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_ME_V5_HANDOVER                            */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_handover(m,r)
     struct model_struct *m;      // Model params
     struct response_struct *r;   // Response data
{
  int i;
  int xi,yi,ci,li,tn,id,xn,yn,seed,ns,vn,stfi,chi,strn;
  float *t,*tt,*fs,*v;
  char *tid,*plname;
  struct mod_me_v5_struct *mv5;
  struct mod_me_v1_input ***v1resp;
  struct mod_me_v1resp *tv1;
  struct pop_layer *pl;
  struct mod_me_v5_pop *mtresp;
  struct onode *sgo,*popo;

  mv5 = mod_me_v5;

  if (mod_me_ftop == NULL){
    v1resp = mv5->v1resp;        // old way
  }else{
    stfi = 0; // *** WYETH FIX
    v1resp = mv5->v1stf->sgr[stfi]->bdgr; // new way
  }

  xn = mod_me_xn;
  yn = mod_me_yn;
  
  tt = NULL;  // Optional temporary storage

  tn = mod_me_tn;

  for(i=0;i<r->n;i++){ // For each response requested
    //printf("=====> %s\n",r->datid[i]);
    //printf("   ==> %d\n",r->rformat[i]);

    if (r->rformat[i] == 0){
      //
      //  'save_as_'
      //

      if (compare_prefix_string_order("v1r_raw_", r->datid[i]) ||
	  compare_prefix_string_order("v1r_norm_",r->datid[i]) ||
	  compare_prefix_string_order("v1r_opp_", r->datid[i]) ||
	  compare_prefix_string_order("v1r_bde_", r->datid[i]) ||
	  compare_prefix_string_order("v1r_me_", r->datid[i]) ||
	  compare_prefix_string_order("v1_binoc_",r->datid[i]) ||
	  (strcmp(r->datid[i],"v1r_sum")==0)){
	if (mv5->binoc == 0){
	  printf("  *** Data ID requires Right Eye:  %s\n",r->datid[i]);
	  exit_error("MOD_ME_V5_HANDOVER","Stimulus is left eye only");
	}
      }

      /*** WYETH HERE REMOVE THIS, eventually ... Aug 24 ***/
      if (compare_prefix_string_order("v1_binoc_",r->datid[i])){
	exit_error("MOD_ME_V5_HANDOVER","No 'v1_binoc_...' response");
      }


      if (compare_prefix_string_order("v1_raw_",  r->datid[i]) ||
	  compare_prefix_string_order("v1_rawodd_", r->datid[i]) ||
	  compare_prefix_string_order("v1r_raw_", r->datid[i]) ||
	  compare_prefix_string_order("v1r_rawodd_", r->datid[i]) ||
	  compare_prefix_string_order("v1_norm_", r->datid[i]) ||
	  compare_prefix_string_order("v1r_norm_",r->datid[i]) ||
	  compare_prefix_string_order("v1_normodd_", r->datid[i]) ||
	  compare_prefix_string_order("v1r_normodd_",r->datid[i]) ||

	  // NEW
	  compare_prefix_string_order("v1_opp_even_pos_",  r->datid[i]) ||
	  compare_prefix_string_order("v1_opp_even_neg_",  r->datid[i]) ||
	  compare_prefix_string_order("v1_opp_odd_pos_",  r->datid[i]) ||
	  compare_prefix_string_order("v1_opp_odd_neg_",  r->datid[i]) ||
	  compare_prefix_string_order("v1r_opp_even_pos_",  r->datid[i]) ||
	  compare_prefix_string_order("v1r_opp_even_neg_",  r->datid[i]) ||
	  compare_prefix_string_order("v1r_opp_odd_pos_",  r->datid[i]) ||
	  compare_prefix_string_order("v1r_opp_odd_neg_",  r->datid[i]) ||
	  compare_prefix_string_order("v1_bde_neg_",  r->datid[i]) ||
	  compare_prefix_string_order("v1_bde_pos_",  r->datid[i]) ||
	  compare_prefix_string_order("v1r_bde_neg_",  r->datid[i]) ||
	  compare_prefix_string_order("v1r_bde_pos_",  r->datid[i]) ||

	  compare_prefix_string_order("v1_opp_",  r->datid[i]) ||
	  compare_prefix_string_order("v1r_opp_", r->datid[i]) ||
	  compare_prefix_string_order("v1_oppodd_",  r->datid[i]) ||
	  compare_prefix_string_order("v1r_oppodd_", r->datid[i]) ||
	  compare_prefix_string_order("v1_bde_",  r->datid[i]) ||
	  compare_prefix_string_order("v1r_bde_", r->datid[i]) ||
	  compare_prefix_string_order("v1_bd_",  r->datid[i]) ||
	  compare_prefix_string_order("v1r_bd_", r->datid[i]) ||
	  compare_prefix_string_order("v1_binoc_", r->datid[i])){

	id = carray_util_underscore_int_final(r->datid[i]);
	if ((id < 1) || (id > mv5->n)){
	  printf("  Data ID:  %s\n",r->datid[i]);
	  printf("    Index:  %d\n",id);
	  exit_error("MOD_ME_V5_HANDOVER","Response index out of range");
	}
	if (compare_prefix_string_order("v1_raw_",  r->datid[i])){
	  t = v1resp[0][0]->dl->raw[id-1][0];   //WYDISP
	}else if (compare_prefix_string_order("v1_rawodd_", r->datid[i])){
	  t = v1resp[0][0]->dl->raw[id-1][1];
	}else if (compare_prefix_string_order("v1r_raw_", r->datid[i])){
	  t = v1resp[0][0]->dr->raw[id-1][0];
	}else if (compare_prefix_string_order("v1r_rawodd_", r->datid[i])){
	  t = v1resp[0][0]->dr->raw[id-1][1];
	}else if (compare_prefix_string_order("v1_norm_", r->datid[i])){
	  t = v1resp[0][0]->dl->norm[id-1][0];
	}else if (compare_prefix_string_order("v1_normodd_", r->datid[i])){
	  t = v1resp[0][0]->dl->norm[id-1][1];
	}else if (compare_prefix_string_order("v1r_norm_", r->datid[i])){
	  t = v1resp[0][0]->dr->norm[id-1][0];
	}else if (compare_prefix_string_order("v1r_normodd_", r->datid[i])){
	  t = v1resp[0][0]->dr->norm[id-1][1];

	  //
	  //  NEW WAY
	  //
	}else if (compare_prefix_string_order("v1_opp_even_pos_",r->datid[i])){
	  t = v1resp[0][0]->dl->opp[id-1][0];
	}else if (compare_prefix_string_order("v1_opp_odd_pos_",r->datid[i])){
	  t = v1resp[0][0]->dl->opp[id-1][1];
	}else if (compare_prefix_string_order("v1_opp_even_neg_",r->datid[i])){
	  t = v1resp[0][0]->dl->opp[id-1][2];
	}else if (compare_prefix_string_order("v1_opp_odd_neg_",r->datid[i])){
	  t = v1resp[0][0]->dl->opp[id-1][3];
	}else if (compare_prefix_string_order("v1r_opp_even_pos_",r->datid[i])){
	  t = v1resp[0][0]->dr->opp[id-1][0];
	}else if (compare_prefix_string_order("v1r_opp_odd_pos_",r->datid[i])){
	  t = v1resp[0][0]->dr->opp[id-1][1];
	}else if (compare_prefix_string_order("v1r_opp_even_neg_",r->datid[i])){
	  t = v1resp[0][0]->dr->opp[id-1][2];
	}else if (compare_prefix_string_order("v1r_opp_odd_neg_",r->datid[i])){
	  t = v1resp[0][0]->dr->opp[id-1][3];
	  
	  //
	  //  OLD WAY  *** WYETH - NOT UPDATED FOR new 'v1resp'
	  //
	}else if (compare_prefix_string_order("v1_opp_", r->datid[i])){
	  exit_error("MOD_ME_V5_HANDOVER",
	     "Use 'v1_opp_even_pos_..' instead of 'v1_opp_..' in .rsp\n");
	}else if (compare_prefix_string_order("v1_oppodd_", r->datid[i])){
	  exit_error("MOD_ME_V5_HANDOVER",
	     "Use 'v1_opp_odd_pos_..' instead of 'v1_oppodd_..' in .rsp\n");
	}else if (compare_prefix_string_order("v1r_opp_", r->datid[i])){
	  exit_error("MOD_ME_V5_HANDOVER",
	     "Use 'v1r_opp_even_pos_...' instead of 'v1r_opp_...' (.rsp)\n");
	}else if (compare_prefix_string_order("v1r_oppodd_", r->datid[i])){
	  exit_error("MOD_ME_V5_HANDOVER",
	     "Use 'v1r_opp_odd_pos_...' instead of 'v1r_oppodd_...' (.rsp)\n");

	  //
	  //  NEW WAY
	  //
	}else if (compare_prefix_string_order("v1_bde_", r->datid[i])){
	  t = v1resp[0][0]->dl->bde[id-1][6];
	}else if (compare_prefix_string_order("v1r_bde_", r->datid[i])){
	  t = v1resp[0][0]->dr->bde[id-1][6];
	}else if (compare_prefix_string_order("v1_bde_pos_", r->datid[i])){
	  t = v1resp[0][0]->dl->bde[id-1][2];
	}else if (compare_prefix_string_order("v1r_bde_pos_", r->datid[i])){
	  t = v1resp[0][0]->dr->bde[id-1][2];
	}else if (compare_prefix_string_order("v1_bde_neg_", r->datid[i])){
	  t = v1resp[0][0]->dl->bde[id-1][5];
	}else if (compare_prefix_string_order("v1r_bde_neg_", r->datid[i])){
	  t = v1resp[0][0]->dr->bde[id-1][5];

	  //
	  //  OLD WAY  *** WYETH - NOT UPDATED FOR new 'v1resp'
	  //
	}else if (compare_prefix_string_order("v1_bde_", r->datid[i])){
	  exit_error("MOD_ME_V5_HANDOVER",
		     "Use 'v1_bde_pos_...' instead of 'v1_bde_...' (.rsp)\n");
	  //t = mv5->v1_bde[id-1][2];
	}else if (compare_prefix_string_order("v1r_bde_", r->datid[i])){
	  exit_error("MOD_ME_V5_HANDOVER",
		     "Use 'v1r_bde_pos_...' instead of 'v1r_bde_...' (.rsp)\n");
	  //t = mv5->v1r_bde[id-1][2];
	}else if (compare_prefix_string_order("v1_bd_even_", r->datid[i])){
	  //t = mv5->v1_bde[id-1][0];
	  //t = mv5->v1resp[0][0]->dl->bde[id-1][0];
	  t = v1resp[0][0]->dl->bde[id-1][0];
	}else if (compare_prefix_string_order("v1r_bd_even_", r->datid[i])){
	  //t = mv5->v1r_bde[id-1][0];
	  //t = mv5->v1resp[0][0]->dr->bde[id-1][0];
	  t = v1resp[0][0]->dr->bde[id-1][0];
	}else if (compare_prefix_string_order("v1_bd_odd_", r->datid[i])){
	  //t = mv5->v1_bde[id-1][1];
	  //t = mv5->v1resp[0][0]->dl->bde[id-1][1];
	  t = v1resp[0][0]->dl->bde[id-1][1];
	}else if (compare_prefix_string_order("v1r_bd_odd_", r->datid[i])){
	  //t = mv5->v1r_bde[id-1][1];
	  //t = mv5->v1resp[0][0]->dr->bde[id-1][1];
	  t = v1resp[0][0]->dr->bde[id-1][1];


	}else if (compare_prefix_string_order("v1_bd_even2_", r->datid[i])){
	  t = v1resp[0][0]->dl->bde[id-1][0];
	  tt = copy_farray(t,tn);	square_farray(tt,tn); t = tt;
	}else if (compare_prefix_string_order("v1r_bd_even2_", r->datid[i])){
	  t = v1resp[0][0]->dr->bde[id-1][0];
	  tt = copy_farray(t,tn);	square_farray(tt,tn); t = tt;
	}else if (compare_prefix_string_order("v1_bd_odd2_", r->datid[i])){
	  t = v1resp[0][0]->dl->bde[id-1][1];
	  tt = copy_farray(t,tn);	square_farray(tt,tn); t = tt;
	}else if (compare_prefix_string_order("v1r_bd_odd2_", r->datid[i])){
	  t = v1resp[0][0]->dr->bde[id-1][1];
	  tt = copy_farray(t,tn);	square_farray(tt,tn); t = tt;
	}else if (compare_prefix_string_order("v1_binoc_", r->datid[i])){
	  exit_error("MOD_ME_V5_HANDOVER",
		     "Do not use 'v1_binoc_...' in .rsp file (Aug 2016).");
	  // t = v1resp[0][0]->binoc[id-1];
	}

	//    mv5 = mod_me_v5;
	// v1resp = mv5->v1resp;

      }else if (strcmp(r->datid[i],"v1_sum")==0){
	t = v1resp[0][0]->dl->sum;
      }else if (strcmp(r->datid[i],"v1r_sum")==0){
	t = v1resp[0][0]->dr->sum;
      }else if (strcmp(r->datid[i],"mt")==0){
	t = v1resp[0][0]->mtb;
      }else if (strcmp(r->datid[i],"mt_nl")==0){
	t = mv5->mt_norm;
      }else if ((strcmp(r->datid[i],"spikes")==0) ||
		(strcmp(r->datid[i],"poisson_prob")==0)){

	t = mv5->mt_norm;
	if (t == NULL){
	  printf("  *** For dataID = %s\n",r->datid[i]);
	  exit_error("MOD_ME_V5_HANDOVER","Signal 'mt_norm' is null");
	}

	seed = m->mseed[r->tsi];
	sgo = onode_child_get_unique(m->o,"spike_gen");
	
	mod_me_gen_poisson(m,r,sgo,t,tn,seed,&fs,&ns,&v,&vn);
	//mod_me_gen_poisson_dummy(m,r,sgo,t,tn,seed,&fs,&ns,&v,&vn);

	if (strcmp(r->datid[i],"spikes")==0){
	  mod_util_resp_store_s(r,i,fs,ns,1,mylogf);  // copyflag = 1
	}else if (strcmp(r->datid[i],"poisson_prob")==0){
	  mod_util_resp_store_f(r,i,v,vn,1,mylogf);   // copyflag = 1
	}
	myfree(fs);
	myfree(v);

	t = NULL;  // Prevent any storage of 't' below.
      }else
	t = NULL;

    }else if ((r->rformat[i] == 1) || (r->rformat[i] == 3)){
      //
      //  'save_pop_unit_as_'
      //  'save_pop_layer_as_'
      //

      xi = r->xi[i];
      yi = r->yi[i];
      ci = r->zi[i];

      plname = r->plname[i];

      if (mod_me_mpt == NULL){
	//
	//  There is no MT population, so only take responses at [0][0]
	//
	if ((xi != 0) || (yi != 0))
	  exit_error("MOD_ME_V5_HANDOVER",
		     "xi and yi must be 0, model has no MT summation");
      }

      if (compare_prefix_string_order("v1",plname)){
	//if ((strcmp(plname,"v1")==0)||(strcmp(plname,"v1r")==0)){
	if ((xi < 0) || (xi >= xn) || (yi < 0) || (yi >= yn))
	  exit_error("MOD_ME_V5_HANDOVER",
	      "xi and yi must be in [0,xn-1] and [0,yn-1], respectively");
	if (v1resp[xi][yi]->flag == 0){
	  printf("  ***  xi,yi = %d, %d\n",xi,yi);
	  exit_error("MOD_ME_V5_HANDOVER",
		     "No V1 response was saved for the requested location.");
	}

	//
	//  Attempt to dis-allow inappropriate response names.
	//
	tid = r->datid[i];
	if (mv5->v1_dflag == 1){
	  if ((strcmp(tid,"raw")==0)||
	      (strcmp(tid,"norm")==0)||
	      (strcmp(tid,"opp")==0)){
	    printf("  *** Response data ID:  %s  (pop name:  %s)\n",tid,
		   plname);
	    printf("  *** You must specify '_even' or '_odd' for disparity.\n");
	    exit_error("MOD_ME_V5_HANDOVER","Inappropriate 'datid'");
	  }
	}else{
	  if ((strcmp(tid,"raw_even"    )==0)||
	      (strcmp(tid,"raw_odd"     )==0)||
	      (strcmp(tid,"norm_even"   )==0)||
	      (strcmp(tid,"norm_odd"    )==0)||
	      (strcmp(tid,"opp_pos_even")==0)||
	      (strcmp(tid,"opp_pos_odd" )==0)||
	      (strcmp(tid,"opp_neg_even")==0)||
	      (strcmp(tid,"opp_neg_odd" )==0)||

	      (strcmp(tid,"bd_pos_even" )==0)||
	      (strcmp(tid,"bd_pos_odd"  )==0)||
	      (strcmp(tid,"bd_neg_even" )==0)||
	      (strcmp(tid,"bd_neg_odd"  )==0)||
	      (strcmp(tid,"bde_pos"     )==0)||
	      (strcmp(tid,"bde_neg"     )==0)||
	      (strcmp(tid,"bde"         )==0)){
	    printf("  *** Response data ID:  %s  (pop name:  %s)\n",tid,
		   plname);
	    printf("  *** This is not a disparity model.\n");
	    printf("  *** Please try one of:  'raw', 'norm' or 'opp'\n");
	    exit_error("MOD_ME_V5_HANDOVER","Inappropriate 'datid'");
	  }
	}


	if (strcmp(plname,"v1")==0){
	  tv1 = v1resp[xi][yi]->dl;
	}else if (strcmp(plname,"v1r")==0){
	  tv1 = v1resp[xi][yi]->dr;
	}else{
	  //
	  //  Decode the name of the form:  'v1r_c12'
	  //
	  strn = strlen(plname);
	  if (strn < 5){
	    printf("  *** plname:  %s\n",plname);
	    exit_error("MOD_ME_V5_HANDOVER","Bad 'plname'");
	  }
	  if (plname[2] == 'r'){
	    if (plname[4] != 'c')
	      exit_error("MOD_ME_V5_HANDOVER","Expecting 'c' at pos 4");
	    chi = atoi(&(plname[5]));
	    if (chi >= mv5->v1stf->n){
	      printf("Offending channel index:  %d\n",chi);
	      exit_error("MOD_ME_V5_HANDOVER","STF channel index too large");
	    }
	    tv1 = mv5->v1stf->sgr[chi]->bdgr[xi][yi]->dr; // new way
	  }else{
	    if (plname[3] != 'c')
	      exit_error("MOD_ME_V5_HANDOVER","Expecting 'c' at pos 3");
	    chi = atoi(&(plname[4]));
	    if (chi >= mv5->v1stf->n){
	      printf("Offending channel index:  %d\n",chi);
	      exit_error("MOD_ME_V5_HANDOVER","STF channel index too large");
	    }
	    tv1 = mv5->v1stf->sgr[chi]->bdgr[xi][yi]->dl; // new way
	  }
	}

	/***

	    The following lists the names and the data structure
	    elements that are returned for the following save
	    commands:
            
	    'save_pop_unit_as_'
	    'save_pop_layer_as_'

            Note, 'xi' 'yi' and 'ci' are the 3 integer coordinate values

	    Pop Name
	     v1            tv1 = v1resp[xi][yi]->dl
	     v1r           tv1 = v1resp[xi][yi]->dr;

            Data ID
	     raw           tv1->raw[ci][0]
	     raw_even       "
	     raw_odd       tv1->raw[ci][1]
	     norm          tv1->norm[ci][0]
	     norm_even      "
	     norm_odd      tv1->norm[ci][1]
	     opp           tv1->opp[ci][0]
	     opp_pos_even   "
	     opp_pos_odd   tv1->opp[ci][1]
	     opp_neg_even  tv1->opp[ci][2]
	     opp_neg_odd   tv1->opp[ci][3]
	     sum           tv1->sum            'ci' is ignored
	     sum_glob      v1stf->sum_norm[xi][yi]
	     sum_glob_r    v1stf->sum_norm_r[xi][yi]
	     bd_pos_even   tv1->bde[ci][0]
	     bd_pos_odd    tv1->bde[ci][1]
	     bde_pos       tv1->bde[ci][2]
	     bd_neg_even   tv1->bde[ci][3]
	     bd_neg_odd    tv1->bde[ci][4]
	     bde_neg       tv1->bde[ci][5]
	     bde           tv1->bde[ci][6]
	     mtsub_monoc   tv1->mt;
	     mtsub_binoc   v1resp[xi][yi]->mtb;
	***/

	if ((strcmp(r->datid[i],"raw"     )==0)||
	    (strcmp(r->datid[i],"raw_even")==0)){
	  t = tv1->raw[ci][0];
	}else if (strcmp(r->datid[i],"raw_odd")==0){
	  t = tv1->raw[ci][1];
	}else if (strcmp(r->datid[i],"surr")==0){
	  if (mv5->v1srnd == NULL)
	    exit_error("MOD_ME_V5_HANDOVER","No surround signal");

	  if (strcmp(plname,"v1")==0){
	    t = mv5->v1srnd->sl[ci][xi][yi];
	  }else if (strcmp(plname,"v1r")==0){
	    t = mv5->v1srnd->sr[ci][xi][yi];
	  }

	}else if (strcmp(r->datid[i],"adresp")==0){
	  if (mv5->v1_adflag == 0)
	    exit_error("MOD_ME_V5_HANDOVER","No adaptation");
	  t = tv1->adrsp[ci][0];

	}else if ((strcmp(r->datid[i],"norm"     )==0)||
		  (strcmp(r->datid[i],"norm_even")==0)){
	  t = tv1->norm[ci][0];
	}else if (strcmp(r->datid[i],"norm_odd")==0){
	  t = tv1->norm[ci][1];
	}else if (strcmp(r->datid[i],"opp_pos_even")==0){
	  t = tv1->opp[ci][0];
	}else if ((strcmp(r->datid[i],"opp")==0)||
		  (strcmp(r->datid[i],"opp_pos_even")==0)){
	  t = tv1->opp[ci][0];
	}else if (strcmp(r->datid[i],"opp_pos_odd")==0){
	  t = tv1->opp[ci][1];
	}else if (strcmp(r->datid[i],"opp_neg_even")==0){
	  t = tv1->opp[ci][2];
	}else if (strcmp(r->datid[i],"opp_neg_odd")==0){
	  t = tv1->opp[ci][3];

	}else if (strcmp(r->datid[i],"sum")==0){
	  // 'ci' is ignored.
	  t = tv1->sum;
	}else if (strcmp(r->datid[i],"sum_glob")==0){
	  t = mv5->v1stf->sum_norm[xi][yi];
	}else if (strcmp(r->datid[i],"sum_glob_r")==0){
	  t = mv5->v1stf->sum_norm_r[xi][yi];
	}else if (strcmp(r->datid[i],"bd_pos_even")==0){
	  t = tv1->bde[ci][0];
	}else if (strcmp(r->datid[i],"bd_pos_odd")==0){
	  t = tv1->bde[ci][1];
	}else if (strcmp(r->datid[i],"bde_pos")==0){
	  t = tv1->bde[ci][2];
	}else if (strcmp(r->datid[i],"bd_neg_even")==0){
	  t = tv1->bde[ci][3];
	}else if (strcmp(r->datid[i],"bd_neg_odd")==0){
	  t = tv1->bde[ci][4];
	}else if (strcmp(r->datid[i],"bde_neg")==0){
	  t = tv1->bde[ci][5];
	}else if (strcmp(r->datid[i],"bde")==0){
	  t = tv1->bde[ci][6];
	}else if (strcmp(r->datid[i],"mtsub_monoc")==0){
	  t = tv1->mt;
	}else if (strcmp(r->datid[i],"mtsub_binoc")==0){
	  t = v1resp[xi][yi]->mtb;
	}else{
	  printf("  *** datid:  %s\n",r->datid[i]);
	  exit_error("MOD_ME_V5_HANDOVER","Unknown 'datid'");
	}
	/********** WYETH - THE ABOVE LIST IS INCOMPLETE ***************/
	/********** WYETH - THE ABOVE LIST IS INCOMPLETE ***************/

      }else{

	//  For name 'mt' populations:
	//
	//    raw            mtresp->ssum[xi][yi][ci]
	//    norm           mtresp->nl[xi][yi][ci]
	//    spikes         mtresp->nl[xi][yi][ci]   ==> SPIKE GEN.
	//    poisson_prob   mtresp->nl[xi][yi][ci]   ==> SPIKE GEN.
	//

	if (mod_me_mpt != NULL){
	  li = popu_get_layer_index_by_name(mod_me_mpt,plname);
	  if (li < 0){
	    printf("  %s\n",plname);
	    mylog_exit(mylogf,"MOD_ME_V5_HANDOVER  Unknown population name.");
	  }
	  pl = mod_me_mpt->lay[li];  // + mv5->mtli0];
	  mtresp = mv5->mtp[li - mv5->mtli0];

	  if ((xi < 0) || (xi >= pl->xn) || (yi < 0) || (yi >= pl->yn))
	    exit_error("MOD_ME_V5_HANDOVER",
		       "xi or yi exceed size of 'mt' population");
	}else{
	  if ((xi != 0) || (yi != 0))
	    exit_error("MOD_ME_V5_HANDOVER",
		       "xi and yi must be 0, there is no MT population");
	}


	if (strcmp(r->datid[i],"raw")==0){
	  t = mtresp->ssum[xi][yi][ci];
	}else if (strcmp(r->datid[i],"norm")==0){
	  if (mod_me_mpt != NULL){
	    t = mtresp->nl[xi][yi][ci];
	  }else{
	    t = mv5->mt_norm;  // If there is no MT population
	  }
	}else if ((strcmp(r->datid[i],"spikes")==0) ||
		  (strcmp(r->datid[i],"poisson_prob")==0)){

	  if (mod_me_mpt != NULL){
	    t    = mtresp->nl[xi][yi][ci];
	    seed = mtresp->seed[xi][yi][ci];
	    popo = onode_get_node_type_item_val(m->o,"pop","name",plname);
	    sgo  = onode_child_get_unique(popo,"spike_gen");
	  }else{
	    //
	    //  WYETH - there is no MT pop.
	    //
	    t    = mv5->mt_norm;
	    seed = m->mseed[r->tsi];
	    sgo  = onode_child_get_unique(m->o,"spike_gen");
	  }

	  mod_me_gen_poisson(m,r,sgo,t,tn,seed,&fs,&ns,&v,&vn);

	  if (strcmp(r->datid[i],"spikes")==0){
	    mod_util_resp_store_s(r,i,fs,ns,1,mylogf);  // copyflag = 1
	  }else if (strcmp(r->datid[i],"poisson_prob")==0){
	    mod_util_resp_store_f(r,i,v,vn,1,mylogf);   // copyflag = 1
	  }
	  myfree(fs);
	  myfree(v);

	  t = NULL;  // Prevent any storage of 't' below.

	}else{
	  printf("  *** datid:  %s\n",r->datid[i]);
	  exit_error("MOD_ME_V5_HANDOVER","Unknown 'datid'");
	}
      }
    }else{
      exit_error("MOD_ME_V5_HANDOVER","Unkown response condition.");
    }



    if (t != NULL){
      //printf("i = %d (before store) \n",i);
      // Original data can be freed, it is not retained by this call
      mod_util_resp_store_f_samp(r,i,t,tn,1.0/mod_me_tscale,mylogf);
    }else{
      ; //printf("  NOTE:  t is null for dataID:  %s\n",r->datid[i]);
    }
    
    if (tt != NULL){  // This needs to be checked for each response
      myfree(tt);
      tt = NULL;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_ME_GEN_HAND_POISSON                         */
/*                                                                           */
/*****************************************************************************/
    void mod_me_gen_hand_poisson(m,r,dp,da,n,diff_flag)
     struct model_struct *m;      // Model params
     struct response_struct *r;   // Response data
     float *dp;                   // Data trace [n] preferred
     float *da;                   // Data trace [n] anti-pref
     int n;                       // length of data
     int diff_flag;               // 0-dp, 1-(dp-da)
{
  int i;
  int ns,vn,ri,seed;
  float *fs,*v,samp,*d;
  struct onode *sgo;

  samp = 1.0/mod_me_tscale;
  seed = m->mseed[r->tsi];

  if (diff_flag == 1){
    if ((dp == NULL) || (da == NULL))
      mylogx(mylogf,"MOD_ME_GEN_HAND_POISSON","Null data received");
    d = subtract_farrays(dp,da,n);   // Pref - Anti
  }else{
    d = copy_farray(dp,n);           // Pref only
  }


  // 'vflag' = 1, so v and vn are returned
  // Spikes are returned in msec

  if (m->ppl == NULL){
    sgo = onode_child_get_unique(m->o,"spike_gen");
    ifc_util_poisson(mylogf,m,sgo,d,n,samp,seed,0,1,&fs,&ns,&v,&vn);
  }else{
    ifc_util_poisson(mylogf,m,NULL,d,n,samp,seed,0,1,&fs,&ns,&v,&vn);
  }

  sprintf(ggstr,"    Made %d Poisson spikes\n",ns);
  mylog(mylogf,ggstr);

  // Assuming spikes are returned in milliseconds
  add_const_farray(fs,ns,mod_me_tdelay*1000.0); // Retinal latency


  for(i=0;i<r->n;i++){ // For each response requested
    if (strcmp(r->datid[i],"spikes")==0){
      mod_util_resp_store_s(r,i,fs,ns,1,mylogf);  // copyflag = 1
    }else if (strcmp(r->datid[i],"poisson_prob")==0){
      mod_util_resp_store_f(r,i,v,vn,1,mylogf);  // copyflag = 1
    }
  }

  myfree(d);
  myfree(fs);
  myfree(v);
}
/**************************************-**************************************/
/*                                                                           */
/*                            MODEL_ME_GEN_HAND_IFC                          */
/*                                                                           */
/*****************************************************************************/
void mod_me_gen_hand_ifc(m,r,dp,da,tn)
     struct model_struct *m;      // Model params
     struct response_struct *r;   // Response data
     float *dp;                   // Data trace [n] preferred
     float *da;                   // Data trace [n] anti-pref
     int tn;                      // length of data
{
  int i;
  int ns,fvn,seed;
  float *gti,*gtx,*fs,*fv,samp,samp_out;

  samp = 1.0/mod_me_tscale;
  seed = m->mseed[r->tsi];

  if (mod_me_oppflag == 0){
    gtx = diff_farrays(dp,da,tn);  // Opponent ME
    gti = get_zero_farray(tn);
    //append_farray_plot("xxx.pl","ME_opp",gtx,tn,1);
  }else if (mod_me_oppflag == 2){
    gtx = diff_farrays(dp,da,tn);  // Opponent ME
    gti = diff_farrays(da,dp,tn);  // -gtx
  }else if (mod_me_oppflag == 3){
    gtx = copy_farray(dp,tn);      // Use only pref.
    gti = get_zero_farray(tn);
  }else{
    gtx = copy_farray(dp,tn);
    gti = copy_farray(da,tn);
  }

  // Get spikes in milliseconds, get Vm ('fv'), and gtx, gti are updated
  ifc_test(myid,m,"ifc",gtx,gti,samp,tn,seed,1,NULL,&fs,&ns,1,&fv,&fvn);

  for(i=0;i<r->n;i++){ // For each response requested

    if (strcmp(r->datid[i],"spikes")==0){

      add_const_farray(fs,ns,mod_me_tdelay*1000.0);  // Add latency in msec
      samp_out = r->samp[i];  // Adjust spike sampling for output
      multiply_farray(fs,ns,samp_out/1000.0);
      mod_util_resp_store_s(r,i,fs,ns,1,mylogf);  // copyflag = 1

    }else if (strcmp(r->datid[i],"ifc_gi")==0){
      mod_util_resp_store_f_samp(r,i,gti,tn,samp,mylogf);
    }else if (strcmp(r->datid[i],"ifc_gx")==0){
      mod_util_resp_store_f_samp(r,i,gtx,tn,samp,mylogf);
    }else if (strcmp(r->datid[i],"ifc_vm")==0){
      mod_util_resp_store_f_samp(r,i,fv,fvn,samp,mylogf);
    }
  }

  myfree(fs);
  myfree(fv);
  myfree(gtx);
  myfree(gti);
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_ME_RESP_BINOC_GAIN_TEST                       */
/*                                                                           */
/*****************************************************************************/
void mod_me_resp_binoc_gain_test(s)
     struct stim_struct *s;      // Stimulus params
{
  int ti,xi,yi;
  int tn,xn,yn,t0;
  float ***fr0,***fr1,***fr2,***fr3;
  float *sum0,*sum1,*sum2,*sum3,*tot;
  char tname[SLEN];

  float sl,sr;

  tn = mod_me_tn;
  xn = mod_me_xn;
  yn = mod_me_yn;

  tot  = get_zero_farray(tn);
  for(ti=1;ti<=tn;ti++){
    t0 = ti-1;
    for(xi=1;xi<=xn;xi++){
      for(yi=1;yi<=yn;yi++){
	sl =  s->d[xi][yi][ti] - 0.5;
	sr =  s->d_r[xi][yi][ti] - 0.5;
	tot[t0] += (sl*sl);
	tot[t0] += (sr*sr);
      }
    }
  }

  sprintf(tname,"Tot_%d",s->stimno);
  append_farray_plot("zz.gain_binoc.pl",tname,tot,tn,1);
  

  //
  //  Integrate the filter response
  //
  /*
  fr0 = mod_me_ff[0]->r;
  fr1 = mod_me_ff[1]->r;
  fr2 = mod_me_ff[2]->r;
  fr3 = mod_me_ff[3]->r;

  sum0 = get_zero_farray(tn);
  sum1 = get_zero_farray(tn);
  sum2 = get_zero_farray(tn);
  sum3 = get_zero_farray(tn);
  tot  = get_zero_farray(tn);

  for(ti=1;ti<=tn;ti++){
    t0 = ti-1;
    for(xi=1;xi<=xn;xi++){
      for(yi=1;yi<=yn;yi++){
	sum0[t0] += fr0[xi][yi][ti] * fr0[xi][yi][ti];
	sum1[t0] += fr1[xi][yi][ti] * fr1[xi][yi][ti];
	sum2[t0] += fr2[xi][yi][ti] * fr2[xi][yi][ti];
	sum3[t0] += fr3[xi][yi][ti] * fr3[xi][yi][ti];
      }
    }
    tot[t0] = sum0[t0] + sum1[t0] + sum2[t0] + sum3[t0];
  }

  //remove_file("zz.gain_binoc.pl");
  append_farray_plot("zz.gain_binoc.pl","sum0",sum0,tn,1);
  append_farray_plot("zz.gain_binoc.pl","sum1",sum1,tn,1);
  append_farray_plot("zz.gain_binoc.pl","sum2",sum2,tn,1);
  append_farray_plot("zz.gain_binoc.pl","sum3",sum3,tn,1);
  sprintf(tname,"Tot_%d",s->stimno);
  append_farray_plot("zz.gain_binoc.pl",tname,tot,tn,1);
  */

}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_ME_RESP_BINOC_FILT                          */
/*                                                                           */
/*****************************************************************************/
void mod_me_resp_binoc_filt(xi,yi,rr)
     int xi,yi;   // Coordinates of response to process
     float **rr;  // [tn] Return the response
{
  int i;
  int tn,sqflag;
  float thresh,*pos1,*neg1,*pos2,*neg2,*r;

  tn = mod_me_tn;
  thresh = mod_me_binoc_athresh;

  if (strcmp(mod_me_binoc_nonlin,"halfsq")==0)
    sqflag = 1;
  else if (strcmp(mod_me_binoc_nonlin,"square")==0)
    sqflag = 2;
  else if (strcmp(mod_me_binoc_nonlin,"none")==0)
    sqflag = 0;

  pos1 = mod_me_resp_binoc_sig(xi,yi,"b1p",thresh,sqflag);  // 1=halfsquare
  neg1 = mod_me_resp_binoc_sig(xi,yi,"b1n",thresh,sqflag);
  pos2 = mod_me_resp_binoc_sig(xi,yi,"b2p",thresh,sqflag);
  neg2 = mod_me_resp_binoc_sig(xi,yi,"b2n",thresh,sqflag);

  r = (float *)myalloc(tn*sizeof(float));
  for(i=0;i<tn;i++)
    r[i] = pos1[i] + neg1[i] + pos2[i] + neg2[i];

  *rr = r;

  myfree(pos1); myfree(neg1); myfree(pos2); myfree(neg2);
}
/**************************************-**************************************/
/*                                                                           */
/*                         MODEL_ME_01_GET_RESPONSE                          */
/*                                                                           */
/*  Calling sequence:                                                        */
/*    1.  model_me_01_prep                                                   */
/*    2.  model_me_01_prep_ft                                                */
/*    3.  this routine                                                       */
/*        ...                                                                */
/*    N.  model_me_01_done                                                   */
/*                                                                           */
/*****************************************************************************/
void model_me_01_get_response(m,s,r,k)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response data
     int k;                      // repeat number, k=0 means new stim
{
  int i,jx,jy;
  int xn,yn,tn,fn,x,y,tnr,diff_flag;
  float ***rpe,***rpo,***rne,***rno,*mep,*men,*mep_raw,*men_raw;
  float ***rg0,***rg1,*r1pe,*r1po,*r1ne,*r1no,*gg;
  char *spikegen;
  struct onode *sgen;
  struct mod_me_fstruct **fs;  // Filter list

  xn = mod_me_xn;
  yn = mod_me_yn;
  tn = mod_me_tn;
  fs = mod_me_ff;
  fn = mod_me_fn;

  /*
  if ((strcmp(mod_me_modtype,"binoc_filter")==0)){
    // WYETH TESTING
    mod_me_resp_binoc_gain_test(s);
  }
  */

  //
  //  Compute new response for each filter, only if needed
  //
  if ((k==0) || (fs[0]->r == NULL) || (s->stimno != mod_me_stimno_prev)){

    for(i=0;i<fn;i++)  // For each filter
      if (fs[i]->r != NULL)
	free_f3tensor(fs[i]->r,1,xn,1,yn,1,tn);  // Free old response

    if ((strcmp(mod_me_modtype,"binoc_filter")==0)){
      fs[0]->flag = fs[1]->flag = 1;  // Left eye
      fs[2]->flag = fs[3]->flag = 0;
      mod_me_compute_response_from_ft(s->d,fs,fn,xn,yn,tn);   // Get new resp
      fs[0]->flag = fs[1]->flag = 0;
      fs[2]->flag = fs[3]->flag = 1;  // Right eye
      mod_me_compute_response_from_ft(s->d_r,fs,fn,xn,yn,tn); // Get new resp
    }else{
      mod_me_compute_response_from_ft(s->d,fs,fn,xn,yn,tn);   // Get new resp
    }
  }

  //
  //  Convenient pointers - Compat w/ old way
  //
  rpe = fs[0]->r;
  rpo = fs[1]->r;
  rne = fs[2]->r;
  rno = fs[3]->r;
  if (fn == 6){
    rg0 = fs[4]->r;
    rg1 = fs[5]->r;
  }

  //
  //  Compute desired quantity from filters
  //
  x = xn/2;
  y = yn/2;

  if ((strcmp(mod_me_modtype,"binoc_filter")==0)){

    // WYETH TESTING
    //mod_me_resp_binoc_gain_test(s);

    // Compute response at (x,y) from '..._ff' struct
    mod_me_resp_binoc_filt(x,y,&mep_raw);
    men_raw = get_zero_farray(tn);

  }else{
    //
    //  ME Squaring
    //
    mep_raw = (float *)myalloc(tn*sizeof(float));
    men_raw = (float *)myalloc(tn*sizeof(float));

    for(i=0;i<tn;i++){
      mep_raw[i] = rpe[x][y][1+i]*rpe[x][y][1+i]+rpo[x][y][1+i]*rpo[x][y][1+i];
      men_raw[i] = rne[x][y][1+i]*rne[x][y][1+i]+rno[x][y][1+i]*rno[x][y][1+i];
    }
  }
  tnr = tn;


  //  WYETH - will this interact w/ normalization of the filters??
  //          for example, they are normalized by the sum of their squares.
  //        - added 2006 Jun 19, allows sqrt, etc, after squaring
  if (mod_me_pow != 1.0){
    mep = (float *)myalloc(tn*sizeof(float));
    men = (float *)myalloc(tn*sizeof(float));
    mylog(mylogf,"  WYETH NEW Applying power %f after squaring.\n",mod_me_pow);
    for(i=0;i<tnr;i++){
      mep[i] = pow(mep_raw[i],mod_me_pow);
      men[i] = pow(men_raw[i],mod_me_pow);
    }
  }else{
    mep = copy_farray(mep_raw,tn);
    men = copy_farray(men_raw,tn);
  }


  //
  //   Gain control
  //
  if (mod_me_normflag > 0){
    float *smooth,tnorm,tt;
    float mean,sdev,tsig;
    int sp1,sp2;

    sp1 = x - mod_me_norm_xn/2;
    sp2 = sp1 + mod_me_norm_xn;

    gg  = (float *)myalloc(tn*sizeof(float));
    for(i=0;i<tn;i++){
      tnorm = 0.0;
      
      for(jx=sp1;jx<sp2;jx++){ // Spatial integration
	for(jy=sp1;jy<sp2;jy++){
	  tt = (rg0[jx][jy][1+i]*rg0[jx][jy][1+i] +
		rg1[jx][jy][1+i]*rg1[jx][jy][1+i]);
	  tnorm += pow(tt,mod_me_norm_pow);
	}
      }
      gg[i] = mod_me_norm_c + 
	mod_me_norm_f * tnorm / (float)(mod_me_norm_xn*mod_me_norm_xn);
      if (gg[i] < 1.0)
	gg[i] = 1.0;
    }

    if (myid == -1)
      append_farray_plot("zz.norm.pl","gg",gg,tn,1);

    
    // Temporal integration
    if (mod_me_norm_tsd > 0.0){
      tsig = mod_me_norm_tsd / mod_me_tscale;
      smooth = smooth_with_gaussian(gg,tn,tsig,0.01);
      myfree(gg);
      gg = smooth;
    }

    if (myid == -1)
      append_farray_plot("zz.norm.pl","gg_Smooth",gg,tn,1);

    if (mod_me_norm_app == 1){
      for(i=0;i<tn;i++){
	mep[i] /= gg[i];
	men[i] /= gg[i];
      }
    }else{
      exit_error("MODEL_ME_01_GET_RESPONSE","Bad me_norm_app type");
    }

    mean_sdev_farray(gg+50,tn-100,&mean,&sdev);
    sprintf(ggstr,"  GAIN SIGNAL STATS:  %f %f\n",mean,sdev);
    mylog(mylogf,ggstr);
    mean_sdev_farray(mep+50,tn-100,&mean,&sdev);
    sprintf(ggstr,"  ME-P SIGNAL STATS:  %f %f\n",mean,sdev);
    mylog(mylogf,ggstr);

    myfree(gg);
  }
  // Now, 'mep' and 'men' are PREF and NULL ME signals, length 'tnr'

  mod_me_stimno_prev = s->stimno;


  //
  //  Handover any basic float responses
  //
  if ((strcmp(mod_me_modtype,"binoc_filter")==0)){
    mod_me_handover_binoc(m,r,x,y,mep,xn,yn,tn);
  }else{
    mod_me_handover_basic(m,r,mep_raw,men_raw,mep,men,rpe,rpo,rne,rno,xn,yn,tn);
  }


  //
  //  Generate/handover spike responses
  //
  if (m->ppl == NULL){
    sgen = onode_child_get_unique(m->o,"spike_gen");
    spikegen = onode_getpar_chr_exit(sgen,"type");
  }else
    spikegen = paramfile_get_char_param_or_exit(m->ppl,"spike_gen");
  if (strcmp(spikegen,"poisson")==0){
    if (mod_me_oppflag == 0)
      diff_flag = 1;
    else
      diff_flag = 0;
    mod_me_gen_hand_poisson(m,r,mep,men,tn,diff_flag);
  }else if (strcmp(spikegen,"ifc")==0){
    mod_me_gen_hand_ifc(m,r,mep,men,tn);
  }else{
    sprintf(ggstr,"spike_gen = %s\n",spikegen);
    mylog(mylogf,ggstr);
    mylog_exit(mylogf,"MODEL_ME_01_GET_RESPONSE  Unknown spike_gen.");
  }
  myfree(spikegen);

  myfree(mep);     // WYETH - Added Jun 18, 2009, I think it's OK to free
  myfree(men);     // WYETH - Added Jun 18, 2009, I think it's OK to free
  myfree(mep_raw); // WYETH - Added Jun 18, 2009, I think it's OK to free
  myfree(men_raw); // WYETH - Added Jun 18, 2009, I think it's OK to free

  //
  //  Exit if there are responses that have not been filled
  //
  mod_util_resp_unfilled_exit(r,mylogf);
}
/**************************************-**************************************/
/*                                                                           */
/*                         MODEL_ME_01_GET_RESPONSE_X                        */
/*                                                                           */
/*  June 2010 - to test the three-way multiplication of filters.             */
/*                                                                           */
/*****************************************************************************/
void model_me_01_get_response_x(m,s,r,k)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response data
     int k;                      // repeat number, k=0 means new stim
{
  int i,jx,jy;
  int xn,yn,tn,x,y,tnr,diff_flag;
  float ***fpe,***fpo,***fpx,***fne,***fno,***fnx;
  float **fpe_s,**fpo_s,**fpx_s,**fne_s,**fno_s,**fnx_s;
  float ***rpe,***rpo,***rpx,***rne,***rno,***rnx,*mep,*men,*mep_raw,*men_raw;
  float ***rg0,***rg1,*gg;
  char *spikegen;
  struct onode *sgen;

  xn = mod_me_xn;
  yn = mod_me_yn;
  tn = mod_me_tn;
  fpe   = mod_me_fpe;    fpo   = mod_me_fpo;    fpx   = mod_me_fpx;
  fne   = mod_me_fne;    fno   = mod_me_fno;    fnx   = mod_me_fnx;
  fpe_s = mod_me_fpe_s;  fpo_s = mod_me_fpo_s;  fpx_s = mod_me_fpx_s;
  fne_s = mod_me_fne_s;  fno_s = mod_me_fno_s;  fnx_s = mod_me_fnx_s;

  if ((k==0) || (mod_me_rpe == NULL) || (s->stimno != mod_me_stimno_prev)){
    // Compute new response

    if (mod_me_rpe != NULL){
      free_f3tensor(mod_me_rpe,1,xn,1,yn,1,tn);
      free_f3tensor(mod_me_rpo,1,xn,1,yn,1,tn);
      free_f3tensor(mod_me_rpx,1,xn,1,yn,1,tn);
      free_f3tensor(mod_me_rne,1,xn,1,yn,1,tn);
      free_f3tensor(mod_me_rno,1,xn,1,yn,1,tn);
      free_f3tensor(mod_me_rnx,1,xn,1,yn,1,tn);
    }
    if (mod_me_rg0 != NULL){
      free_f3tensor(mod_me_rg0,1,xn,1,yn,1,tn);
      free_f3tensor(mod_me_rg1,1,xn,1,yn,1,tn);
    }
    // Next call modifies s->d
    mod_me_compute_response_quad_opp_from_ft_x(s->d,fpe,fpo,fpx,fne,fno,fnx,
					       fpe_s,fpo_s,fpx_s,
					       fne_s,fno_s,fnx_s,
					       xn,yn,tn,&rpe,&rpo,&rpx,
					       &rne,&rno,&rnx,
					       &rg0,&rg1);

    // WYETH - debug mustafa - IT NEVER GETS HERE
    mylog(mylogf,"  done:  mod_me_compute_response_quad_opp_from_ft\n");

    // Store for future 'repeats'
    mod_me_rpe = rpe; mod_me_rpo = rpo; mod_me_rpx = rpx;
    mod_me_rne = rne; mod_me_rno = rno; mod_me_rnx = rnx;
    mod_me_rg0 = rg0; mod_me_rg1 = rg1;
  }else{
    // Use stored response
    rpe = mod_me_rpe; rpo = mod_me_rpo; rpx = mod_me_rpx;
    rne = mod_me_rne; rno = mod_me_rno; rnx = mod_me_rnx;
    rg0 = mod_me_rg0; rg1 = mod_me_rg1;
  }

  x = xn/2;
  y = yn/2;
  mep_raw = (float *)myalloc(tn*sizeof(float));
  men_raw = (float *)myalloc(tn*sizeof(float));

  // Multiplying
  mylog(mylogf,"****  WYETH HERE multiplying three\n");
  for(i=0;i<tn;i++){
    mep_raw[i] = rpe[x][y][1+i]*rpo[x][y][1+i]*rpx[x][y][1+i];
    men_raw[i] = rne[x][y][1+i]*rno[x][y][1+i]*rnx[x][y][1+i];
  }


  tnr = tn;
  


  //  WYETH - will this interact w/ normalization of the filters??
  //          for example, they are normalized by the sum of their squares.
  //        - added 2006 Jun 19, allows sqrt, etc, after squaring
  if (mod_me_pow != 1.0){
    mep = (float *)myalloc(tn*sizeof(float));
    men = (float *)myalloc(tn*sizeof(float));
    sprintf(ggstr,"  WYETH - NEW Applying power %f after squaring.\n",
	    mod_me_pow);
    mylog(mylogf,ggstr);
    
    for(i=0;i<tnr;i++){
      mep[i] = pow(mep_raw[i],mod_me_pow);
      men[i] = pow(men_raw[i],mod_me_pow);
    }
  }else{
    mep = copy_farray(mep_raw,tn);
    men = copy_farray(men_raw,tn);
  }


  //
  //   Gain control
  //
  if (mod_me_normflag > 0){
    float *smooth,tnorm,tt;
    float mean,sdev,tsig;
    int sp1,sp2;

    sp1 = x - mod_me_norm_xn/2;
    sp2 = sp1 + mod_me_norm_xn;

    gg  = (float *)myalloc(tn*sizeof(float));
    for(i=0;i<tn;i++){
      tnorm = 0.0;
      
      for(jx=sp1;jx<sp2;jx++){ // Spatial integration
	for(jy=sp1;jy<sp2;jy++){
	  tt = (rg0[jx][jy][1+i]*rg0[jx][jy][1+i] +
		rg1[jx][jy][1+i]*rg1[jx][jy][1+i]);
	  tnorm += pow(tt,mod_me_norm_pow);
	}
      }
      gg[i] = mod_me_norm_c + 
	mod_me_norm_f * tnorm / (float)(mod_me_norm_xn*mod_me_norm_xn);
      if (gg[i] < 1.0)
	gg[i] = 1.0;
    }

    if (myid == -1)
      append_farray_plot("zz.norm.pl","gg",gg,tn,1);

    
    // Temporal integration
    if (mod_me_norm_tsd > 0.0){
      tsig = mod_me_norm_tsd / mod_me_tscale;
      smooth = smooth_with_gaussian(gg,tn,tsig,0.01);
      myfree(gg);
      gg = smooth;
    }
    
    if (myid == -1)
      append_farray_plot("zz.norm.pl","gg_Smooth",gg,tn,1);
    
    if (mod_me_norm_app == 1){
      for(i=0;i<tn;i++){
	mep[i] /= gg[i];
	men[i] /= gg[i];
      }
    }else{
      exit_error("MODEL_ME_01_GET_RESPONSE","Bad me_norm_app type");
    }
    
    mean_sdev_farray(gg+50,tn-100,&mean,&sdev);
    sprintf(ggstr,"  GAIN SIGNAL STATS:  %f %f\n",mean,sdev);
    mylog(mylogf,ggstr);
    mean_sdev_farray(mep+50,tn-100,&mean,&sdev);
    sprintf(ggstr,"  ME-P SIGNAL STATS:  %f %f\n",mean,sdev);
    mylog(mylogf,ggstr);
    
    myfree(gg);
  }
  // Now, 'mep' and 'men' are PREF and NULL ME signals, length 'tnr'

  mod_me_stimno_prev = s->stimno;


  //
  //  Handover any basic float responses
  //
  mod_me_handover_basic(m,r,mep_raw,men_raw,mep,men,rpe,rpo,rne,rno,xn,yn,tn);


  //
  //  Generate/handover spike responses
  //
  if (m->ppl == NULL){
    sgen = onode_child_get_unique(m->o,"spike_gen");
    spikegen = onode_getpar_chr_exit(sgen,"type");
  }else
    spikegen = paramfile_get_char_param_or_exit(m->ppl,"spike_gen");
  if (strcmp(spikegen,"poisson")==0){
    if (mod_me_oppflag == 0)
      diff_flag = 1;
    else
      diff_flag = 0;
    mod_me_gen_hand_poisson(m,r,mep,men,tn,diff_flag);
  }else if (strcmp(spikegen,"ifc")==0){
    mod_me_gen_hand_ifc(m,r,mep,men,tn);
  }else{
    printf("spike_gen = %s\n",spikegen);
    mylog_exit(mylogf,"MODEL_ME_01_GET_RESPONSE  Unknown spike_gen.");
  }
  myfree(spikegen);

  myfree(mep);     // WYETH - Added Jun 18, 2009, I think it's OK to free
  myfree(men);     // WYETH - Added Jun 18, 2009, I think it's OK to free
  myfree(mep_raw); // WYETH - Added Jun 18, 2009, I think it's OK to free
  myfree(men_raw); // WYETH - Added Jun 18, 2009, I think it's OK to free

  //
  //  Exit if there are responses that have not been filled
  //
  mod_util_resp_unfilled_exit(r,mylogf);
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_ME_V5_ADAPT_V1_TRACE_S01                      */
/*                                                                           */
/*  Run adaptation algorithm "s01" on the raw data 'd' to produce the        */
/*  adapted trace 'dad' using (and perhaps evolving) the stored adapt state  */
/*  'rastt'.                                                                 */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_adapt_v1_trace_s01(d,dad,tn,rastt,frac,chalf,aflag)
     float *d;      // [tn] Raw data trace to be adapted
     float *dad;    // [tn] Location to write adapted signal
     int tn;        // duration of data traces
     float *rastt;  // adapted state value, THIS GETS OVERWRITTEN
     float  frac;   // parameter for fractional decay
     float  chalf;  // half-rise of non-linearity
     int aflag;     // 1-evolved adaptation, 2-freeze adaptation
{
  int i;
  float asum,totw,adf;

  if (aflag == 1){
    //
    //  Allow adaptation state to evolve and influence response
    //
    totw  = 1.0 / (1.0 - frac);

    asum = *rastt;  // Get the stored value of adaptation state
    for(i=0;i<tn;i++){
      asum = (asum * frac) + d[i]/totw;  // Decay sum, add resp(t=i)
      adf = 1.0 - asum/(asum + chalf);  // compute adaptation factor
      dad[i] = adf * d[i];    // Compute adapted response
      //dad[i] = asum;    // Store the adapt state variable
    }
    *rastt = asum;  // Store the new adaptation state

  }else if (aflag == 2){
    //
    //  Freeze adaptation (do not evolve) but allow it to influence response
    //
    asum = *rastt;  // Get the stored value of adaptation state
    adf = 1.0 - asum/(asum + chalf);  // compute adaptation factor
    for(i=0;i<tn;i++){
      dad[i] = adf * d[i];    // Compute adapted response
      //dadpt[i] = asum;     // Store the adapt state variable
    }
  }else{
    exit_error("MOD_ME_V5_ADAPT_V1_TRACE_S01","Unknown adapt 'aflag' value");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_ME_V1SURR_GETPAR                           */
/*                                                                           */
/*  Return parameter values for V1 surround signal computation.              */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_v1surr_getpar(romix,rsig,rsigc,rtol,rdt)
     float *romix;     // Amount of opponent direction to mix in surround
     float *rsig;      // SD for smoothing (pix)
     float *rsigc;     // SD inner if Annulus (whole in center) (pix) or -1.0
     float *rtol;      // Tolerance for smoothing
     int   *rdt;       // Delay of surround in raw time units
{
  int dt;
  float omix,sigma,sigc,tolerance;
  struct mod_me_v5_struct *mv5;  // V5 params and responses

  mv5 = mod_me_v5;

  omix = mv5->v1srnd->omix;  // Amount of opponent dir. to mix in surr. signal

  sigma = mv5->v1srnd->sig / mod_me_sscale;
  tolerance = 0.04;
  dt = my_rint(mv5->v1srnd->delay / mod_me_tscale);

  sprintf(ggstr,"    Surround SD: %.2f deg (%.1f pix)\n",mv5->v1srnd->sig,
	  sigma);
  mylog(mylogf,ggstr);

  if (mv5->v1srnd->sigc > 0.0){
    sigc  = mv5->v1srnd->sigc / mod_me_sscale;
    sprintf(ggstr,"    Annulus inner SD: %.2f deg (%.1f pix)\n",
	    mv5->v1srnd->sigc,sigc);
    mylog(mylogf,ggstr);

    if (sigc > sigma){
      mylog_exit(mylogf,"MOD_ME_V1SURR_GETPAR  Center SD > Surr SD.");
    }

  }else
    sigc = -1.0;

  sprintf(ggstr,"    Surround delay:  %.3f s (%d time units)\n",
	  mv5->v1srnd->delay,dt);
  mylog(mylogf,ggstr);

  *romix = omix;
  *rsig  = sigma;
  *rsigc = sigc;
  *rtol  = tolerance;
  *rdt   = dt;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_ME_V5_ME_RAW_VECT                           */
/*                                                                           */
/*  WYETH building this for the RPHomeo method                               */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_me_raw_vect(eyeflag,x,y,ti,mepow,mevect)
     int eyeflag;       // 0-left, 1-right
     int x;             // x-coord
     int y;             // y-coord
     int ti;            // Time index
     float mepow;       // Wyeth 2018Apr18
     float *mevect;     // Over-write this vector to return
  
{
  int i;
  int fn,re_flag,foff,ndir;
  float ve,vo;
  struct mod_me_fstruct **fs;  // [fn] List of pointers to filters

  ndir = mod_me_v5->n;         // Number of direction channels

  if (mod_me_ftop == NULL){
    fs  = mod_me_ff;
    fn  = mod_me_fn;
    re_flag = mod_me_fre;
  }else{
    fs = mod_me_ftop->fsl;   // New way
    fn = mod_me_ftop->n;
    re_flag = mod_me_ftop->flag_re;
  }

  foff = 0;  // Offset index in 'fs' array.
  if (eyeflag == 1){
    if (re_flag == 1){
      foff = fn/2;  // R. Eye uses second half of filter set
    }
  }

  if (mod_me_v5->v1stf != NULL) // *** Need different 'fs' below for STF
    exit_error("MOD_ME_V5_ME_RAW_VECT","Not impl'd for Multi STF yet");

  for(i=0;i<ndir/2;i++){
    // Preferred direction
    ve = fs[i*4 + 0 + foff]->r[x+1][y+1][ti+1];
    vo = fs[i*4 + 1 + foff]->r[x+1][y+1][ti+1];
    if (mepow == 1.0)
      mevect[i] = ve*ve + vo*vo;
    else if (mepow == 0.5)
      mevect[i] = sqrt(ve*ve + vo*vo);
    else
      mevect[i] = pow(ve*ve + vo*vo,mepow);

    // Anti-preferred direction
    ve = fs[i*4 + 2 + foff]->r[x+1][y+1][ti+1];
    vo = fs[i*4 + 3 + foff]->r[x+1][y+1][ti+1];
    if (mepow == 1.0)  // Wyeth 2018Apr18
      mevect[i + ndir/2] = ve*ve + vo*vo;
    else if (mepow == 0.5)
      mevect[i + ndir/2] = sqrt(ve*ve + vo*vo);
    else
      mevect[i + ndir/2] = pow(ve*ve + vo*vo,mepow);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_ME_V5_V1_SIMUL                            */
/*                                                                           */
/*  Simultaneous computation of raw, norm and surround response.             */
/*                                                                           */
/*  WYETH HERE RPHomeo                                                       */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_v1_simul(eyeflag)
     int eyeflag;                  // 0-left, 1-right
{
  int i,j;
  int ndir,fi,fj,fio,ti,xn,yn,tn,ddi,save_flag,surrc_i,sdelay;
  float sum,sigma,sigc,tolerance,a,a2,eps,delta,alpha,srt,omix,merpow;
  float *me,**nw,***dnorm,***smt,****surrc,*targ,**tp,****sv_surr;
  char tstr[SLEN],fname[SLEN],*zstr;
  struct mod_me_v5_struct *mv5;  // V5 params and responses
  struct mod_me_v1resp *v1s;     // Where to store responses

  mylog(mylogf,"  MOD_ME_V5_V1_SIMUL\n");

  xn  = mod_me_xn;
  yn  = mod_me_yn;
  tn  = mod_me_tn;
  mv5 = mod_me_v5;
  ndir = mv5->n;

  if (mv5->normo != NULL){  // Wyeth 2018Apr18
    merpow = onode_getpar_flt_dflt(mv5->normo,"me_raw_power",1.0);
    a2     = onode_getpar_flt_exit(mv5->normo,"alpha_2");
    eps    = onode_getpar_flt_exit(mv5->normo,"alpha_4");

    if (mv5->v1srnd != NULL){
      a = mv5->v1srnd->ampl;
    }else{
      a = 0.0;
    }
      
    sprintf(ggstr,"  a    = %f\n",a);
    mylog(mylogf,ggstr);
    sprintf(ggstr,"  a2   = %f\n",a2);
    mylog(mylogf,ggstr);
    sprintf(ggstr,"  eps  = %f\n",eps);
    mylog(mylogf,ggstr);
    sprintf(ggstr,"  me raw pow = %f\n",merpow);
    mylog(mylogf,ggstr);
  }else{
    exit_error("MOD_ME_V5_V1_SIMUL",
	       "Please use <norm> in <v1> to set 'alpha_2' and 'alpha_4'");
  }

  alpha = mv5->v1_adp_t1; // Learning rate for updating normalization weights.


  mod_me_v5_v1surr_getpar(&omix,&sigma,&sigc,&tolerance,&sdelay);
  if (sigc > 0.0){
    mylog_exit(mylogf,"MOD_ME_V5_V1_SIMUL  Annular surround not yet impl'd.");
  }


  me = (float *)myalloc(ndir*sizeof(float));  // Raw ME signal

  dnorm = get_3d_farray(ndir,xn,yn);

  //
  //  Circular storage for surround
  //
  surrc = get_4d_farray(sdelay,ndir,xn,yn);
  surrc_i = 0;   // time index for next write

  if (eyeflag == 0)
    sv_surr = mv5->v1srnd->sl;
  else
    sv_surr = mv5->v1srnd->sr;

  for(ti=0;ti<tn;ti++){

    //
    //  (1) Compute 2D arrays of Norm'd response for each dir channel
    //
    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){

	if (mv5->v1resp[i][j]->flag == 2) // WYETH 2 means storage for 2 eyes
	  save_flag = 1;
	else
	  save_flag = 0;

	if (eyeflag == 0)
	  v1s = mv5->v1resp[i][j]->dl;
	else
	  v1s = mv5->v1resp[i][j]->dr;

	// Get raw ME at all directions
	mod_me_v5_me_raw_vect(eyeflag,i,j,ti,merpow,me);  // Fills 'me[ndir]'

	if (save_flag){  // Save raw ME for this location
	  for(fi=0;fi<ndir;fi++){
	    //printf("fi = %d\n",fi);
	    v1s->raw[fi][0][ti] = me[fi];
	  }
	}

	//  v1_ad_normw  [xn][yn][ndir][ndir] Norm weight for RPHomeo
	nw = mv5->v1_ad_normw[i][j];  // 

	// Get Norm response in each dir-chan
	for(fi=0;fi<ndir;fi++){
	  sum = 0.0;
	  for(fj=0;fj<ndir;fj++){
	    //printf("nw[i][j] = %f\n",nw[fi][fj]);
	    sum += nw[fi][fj] * me[fj];
	  }
	  sum /= (float)ndir;

	  if (ti >= sdelay){
	    srt = surrc[surrc_i][fi][i][j];
	  }else{
	    srt = 0.0;  // Too early, no surround signal yet
	  }

	  if (save_flag){
	    sv_surr[fi][i][j][ti] = srt;  // Save surround signal for .rsp
	  }


	  //
	  //  NORMALIZATION
	  //
	  dnorm[fi][i][j] = me[fi] / (a*srt + a2*sum + eps);

	  if (save_flag){
	    v1s->norm[fi][0][ti] = dnorm[fi][i][j];
	  }
	}
      }
    }

    targ = mv5->v1_adp_tv;  // Target response product [ndir/2 + 1]


    //
    //  (2) Weight update
    //
    if (mv5->v1_adflag == 1){  // If we are adapting... (and not frozen)
      //
      //  Update Norm Weights
      //
      for(i=0;i<xn;i++){
	for(j=0;j<yn;j++){  // For every location
	  
	  nw = mv5->v1_ad_normw[i][j];

	  for(fi=0;fi<ndir;fi++){  // For each direction

	    if (mv5->v1_adp_wwi > 0){
	      //
	      //  Dump out the weights for the central unit.  This is
	      //    controlled by "write_weights" - set to 0 to turn off.
	      //
	      if ((i==xn/2) && (j==yn/2) && ((ti%mv5->v1_adp_wwi)==0)){

		sprintf(tstr,"%d,",ti);
		zstr = get_leading_zero_string_from_int(fi,2);
		sprintf(fname,"zzz.rph.%s.w.pl",zstr);
		append_farray_plot(fname,tstr,nw[fi],ndir,1);
		myfree(zstr);
	      }
	    }

	    for(fj=0;fj<ndir;fj++){  // For each weight to other directions
	      ddi = abs(fi - fj); // Difference in direction index
	      if (ddi > ndir/2)
		ddi = ndir - ddi;
	      //printf("fi,j  %d %d  ddi = %d\n\n",fi,fj,ddi);

	      //if ((dnorm[fi][i][j] > 0.01) && (dnorm[fj][i][j] > 0.01)){
	      delta = dnorm[fi][i][j]*dnorm[fj][i][j] - targ[ddi];
	      //}else{
	      //delta = 0.0;
	      //}

	      //
	      // WYETH - must add a decay term so that the weights decay toward
	      //   their initial uniform distribution.
	      //
	      nw[fi][fj] += alpha * delta;

	      //
	      //  Wyeth - don't allow negative Norm Weights
	      //
	      if (nw[fi][fj] < 0.0)  // THIS IS A RANDOM GUESS????? 
		nw[fi][fj] = 0.0;
	    }
	  }
	}
      }
    }else if (mv5->v1_adflag == 2){
      ; // Don't change the weights - frozen weight condition.
    }


    //
    //  (3) Surround signal computation
    //

    //  Create temporary array of spatially smooth signals for each dir.
    //  Now smooth Norm arrays to get SURR{t} response.
    //  Might not need all (just pref-anti axis), but do all for simplicity.
    smt = (float ***)myalloc(ndir*sizeof(float **));
    for(fi=0;fi<ndir;fi++)
      smt[fi] = smooth_2d_with_gaussian(dnorm[fi],xn,yn,sigma,tolerance);

    for(fi=0;fi<ndir;fi++){  // Optional mixing of opponent directions

      fio = fi + ndir/2;  // 'fio' is the index opponent to 'fi'
      if (fio >= ndir)
	fio -= ndir;

      tp = surrc[surrc_i][fi];
      for(i=0;i<xn;i++){
	for(j=0;j<yn;j++){  // For every location
	  tp[i][j] = smt[fi][i][j] + omix*smt[fio][i][j];
	}
      }
    }
    free_3d_farray(smt,ndir,xn,yn);


    // Increment surround pointer, wrapping around if needed.
    surrc_i += 1;
    if (surrc_i >= sdelay)
      surrc_i = 0;
  }

  myfree(me);
  free_3d_farray(dnorm,ndir,xn,yn);
  free_4d_farray(surrc,sdelay,ndir,xn,yn); // Free circular buffer for Surr
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_ME_V5_V1SURR_SUM                           */
/*                                                                           */
/*  Compute spatial sum for V1 surround signals.                             */
/*                                                                           */
/*  Fill values in:                                                          */
/*    mod_me_v5->v1srnd->sr                                                  */
/*                     ->sl                                                  */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_v1surr_sum(eyeflag)
     int eyeflag;                  // 0-left, 1-right
{
  int i,j,k;
  int d,fi,fio,fn,xn,yn,tn,dt,foff,re_flag,aflag;
  float sigma,sigc,tolerance,omix,tvp,tvn,merpow,merthr,meroff;
  float ***tpe,***tpo,***tne,***tno,**mep,**men,**smep,**smen,****ss;
  float **ctrp,**ctrn,tval,sigs2,sigc2,sigd2;
  float ***admep,***admen,*ttp,*ttn,*tta,frac,chalf,*padst,***tst;
  struct mod_me_fstruct **fs;    // Filter list
  struct mod_me_v5_struct *mv5;  // V5 params and responses
  struct mod_me_v1_surr *surr;   // Pointer

  xn  = mod_me_xn;
  yn  = mod_me_yn;
  tn  = mod_me_tn;
  mv5 = mod_me_v5;
  surr = mv5->v1srnd;

  if (surr == NULL)       // There is no V1 surround suppression
    return;

  if (surr->ampl == 0.0)  // Surround suppression is deactivated
    return;

  mylog(mylogf,"  MOD_ME_V5_V1SURR_SUM\n");

  if (mv5->normo != NULL){  // Wyeth 2018Apr18
    merpow = onode_getpar_flt_dflt(mv5->normo,"me_raw_power",1.0);
    merthr = onode_getpar_flt_dflt(mv5->normo,"me_raw_thresh",0.0);
    meroff = onode_getpar_flt_dflt(mv5->normo,"me_raw_offset",0.0);
  }else{
    merpow = 1.0;
    merthr = 0.0;
    meroff = 0.0;
  }
  //printf("merpow %f   surAmp %f\n",merpow,surr->ampl);

  mod_me_v5_v1surr_getpar(&omix,&sigma,&sigc,&tolerance,&dt);
  //printf("    omix %f sigm %f tol %f  dt %d\n",omix,sigma,tolerance,dt);

  if (mod_me_ftop != NULL){
    exit_error("MOD_ME_V5_V1SURR_SUM","Too new, can't handle STF");
    fs = mod_me_ftop->fsl;
    fn = mod_me_ftop->n;
    re_flag = mod_me_ftop->flag_re;
  }else{
    fs  = mod_me_ff;
    fn  = mod_me_fn;
    re_flag = mod_me_fre;
  }

  foff = 0;  // Offset index in 'fs' array.
  if (eyeflag == 0){
    ss = surr->sl; // Store result in left eye data
  }else if (eyeflag == 1){
    //
    //  Allocate storage the first time we need it.
    //
    if (surr->sr == NULL){
      surr->sr = (float ****)myalloc(surr->n*sizeof(float ***));
      for(i=0;i<surr->n;i++)
	surr->sr[i] = get_3d_farray(xn,yn,tn);
    }
    ss = surr->sr; // Right eye

    if (re_flag == 1){
      foff = fn/2;  // R. Eye uses second half of filter set
    }
  }

  if (mod_me_mpt == NULL)
    exit_error("MOD_ME_V5_V1SURR_SUM","Old way - not implemented");

  if (mv5->v1stf != NULL){
    // Wyeth - the 'fi' index values will have to be reworked for this,
    //   see examples in 'mod_me_v5_raw_resp' below.
    exit_error("MOD_ME_V5_V1SURR_SUM","Not implemented yet for STF models.");
  }

  mep = get_2d_farray(xn,yn);  // Local arrays to hold one spatial frame
  men = get_2d_farray(xn,yn);

  //
  //  New way - with adaptation
  //
  aflag = mv5->v1_adflag;

  if (aflag > 0){
    admep = get_3d_farray(xn,yn,tn);
    admen = get_3d_farray(xn,yn,tn);

    ttp = (float *)myalloc(tn*sizeof(float));
    ttn = (float *)myalloc(tn*sizeof(float));

    if (strcmp(mv5->v1_adtype,"s01")==0){
      chalf = mv5->v1_adp_a1;
      frac  = exp(-mod_me_tscale/mv5->v1_adp_t1);
    }else{
      exit_error("MOD_ME_V5_V1SURR_SUM","Unknown 'aflag'");
    }
  }

  for(fi=0;fi<mv5->n/2;fi++){  //  For HALF of the direction channels

    fio = fi + mv5->n/2;

    tpe = fs[fi*4 + 0 + foff]->r;
    tpo = fs[fi*4 + 1 + foff]->r;
    tne = fs[fi*4 + 2 + foff]->r;
    tno = fs[fi*4 + 3 + foff]->r;

    //
    //  For each time point, create a 2D (spatial) ME frame
    //    Write into 'admen' and 'admep'
    //
    if (aflag > 0){
      for(i=1;i<=xn;i++){
	for(j=1;j<=yn;j++){

	  for(k=1;k<=tn;k++){  // For each frame in the destination array
	    ttp[k-1] = tpe[i][j][k]*tpe[i][j][k] + tpo[i][j][k]*tpo[i][j][k];
	    ttn[k-1] = tne[i][j][k]*tne[i][j][k] + tno[i][j][k]*tno[i][j][k];

	    // This offset defaults to zero  2018Oct29
	    ttp[k-1] += meroff;  // Add offset
	    ttn[k-1] += meroff;  // Add offset

	    //  Threshold to suppress small values, 2018Oct19
	    if (ttp[k-1] < merthr)
	      ttp[k-1] = 0.0;
	    if (ttn[k-1] < merthr)
	      ttn[k-1] = 0.0;

	    if (merpow == 0.5){ // Wyeth 2018Apr18
	      ttp[k-1] = sqrt(ttp[k-1]);
	      ttn[k-1] = sqrt(ttn[k-1]);
	    }else if (merpow == 1.0){
	      ;
	    }else{
	      ttp[k-1] = pow(ttp[k-1],merpow);
	      ttn[k-1] = pow(ttn[k-1],merpow);
	    }

	  }

	  if (eyeflag == 0)
	    tst = mv5->v1_adst_l;
	  else
	    tst = mv5->v1_adst_r;

	  tta = admep[i-1][j-1];
	  padst = &(tst[fi][i-1][j-1]);
	  mod_me_v5_adapt_v1_trace_s01(ttp,tta,tn,padst,frac,chalf,aflag);

	  tta = admen[i-1][j-1];
	  padst = &(tst[fio][i-1][j-1]);
	  mod_me_v5_adapt_v1_trace_s01(ttn,tta,tn,padst,frac,chalf,aflag);
	}
      }
    }

    //  NOW DO SURROUND STUFF FOR THIS CHANNEL (and it's opponent channel)

    for(k=0;k<tn;k++){  // For each frame in the destination array

      d = k+1 - dt;  // Compute index in the source data, factoring in delay

      if (d >= 1){

	//
	//  For each time point, create a 2D (spatial) ME frame
	//
	if (aflag == 0){
	  //
	  //  No adaptation
	  //
	  for(i=1;i<=xn;i++){
	    for(j=1;j<=yn;j++){
	      mep[i-1][j-1] = tpe[i][j][d]*tpe[i][j][d] + 
		tpo[i][j][d]*tpo[i][j][d];
	      men[i-1][j-1] = tne[i][j][d]*tne[i][j][d] +
		tno[i][j][d]*tno[i][j][d];

	      // This offset defaults to zero  2018Oct29
	      mep[i-1][j-1] += meroff;  // Add offset
	      men[i-1][j-1] += meroff;  // Add offset

	      //  Threshold to suppress small values, 2018Oct19
	      //    Reason: when using low 'merpow' values, small input values
	      //    were being amplified too much.
	      if (mep[i-1][j-1] < merthr)
		mep[i-1][j-1] = 0.0;
	      if (men[i-1][j-1] < merthr)
		men[i-1][j-1] = 0.0;
	      
	      if (merpow == 0.5){ // Wyeth 2018Apr18
		mep[i-1][j-1] = sqrt(mep[i-1][j-1]);
		men[i-1][j-1] = sqrt(men[i-1][j-1]);
	      }else if (merpow == 1.0){
		;
	      }else{
		mep[i-1][j-1] = pow(mep[i-1][j-1],merpow);
		men[i-1][j-1] = pow(men[i-1][j-1],merpow);
	      }
	    }
	  }
	}else{
	  //
	  //  Adaptation
	  //
	  for(i=0;i<xn;i++){
	    for(j=0;j<yn;j++){
	      mep[i][j] = admep[i][j][d];
	      men[i][j] = admen[i][j][d];
	    }
	  }
	}

	//  Spatially smooth the ME frame
	//
	//  Note, this does *not* do wrap-around smoothing, the mask is
	//  truncated and renormalized.
	//
	smep = smooth_2d_with_gaussian(mep,xn,yn,sigma,tolerance);
	smen = smooth_2d_with_gaussian(men,xn,yn,sigma,tolerance);

	if (sigc > 0.0){
	  //
	  //  Smooth w/ smaller center Gaussian, subtract for ANNULUS
	  //
	  //  Given that the original surround signal, 'S', is computed
	  //  with a Gaussian of area 1, as is the new center signal, 'C',
	  //  we need to adjust the weights when doing the subtraction, 
	  //  to simulate the case where the equivalent annulus would have
	  //  area = 1.
	  //
	  //               sig_s^2 * S  -  sig_c^2 * C
	  //   Surr_adj =  ---------------------------
	  //                    sig_s^2 - sig_c^2
	  //
	  sigs2 = sigma * sigma;
	  sigc2 = sigc  * sigc;
	  sigd2 = sigs2 - sigc2;
	  //printf("sigd2 = %f  (sig_s  sig_c %f  %f\n",sigd2,sigma,sigc);

	  ctrp = smooth_2d_with_gaussian(mep,xn,yn,sigc,tolerance);
	  ctrn = smooth_2d_with_gaussian(men,xn,yn,sigc,tolerance);

	  for(i=0;i<xn;i++){
	    for(j=0;j<yn;j++){

	      tval = (sigs2 * smep[i][j] - sigc2 * ctrp[i][j]) / sigd2;
	      if (tval < 0.0)
		tval = 0.0;
	      smep[i][j] = tval;

	      tval = (sigs2 * smen[i][j] - sigc2 * ctrn[i][j]) / sigd2;
	      if (tval < 0.0)
		tval = 0.0;
	      smen[i][j] = tval;

	    }
	  }
	  free_2d_farray(ctrp,xn);
	  free_2d_farray(ctrn,xn);
	}

	//  Store the smoothed frame in a 3D x,y,t array

	if (omix == 0.0){
	  //
	  //  Purely DS surround signal
	  //
	  for(i=0;i<xn;i++){
	    for(j=0;j<yn;j++){
	      ss[fi ][i][j][k] = smep[i][j];
	      ss[fio][i][j][k] = smen[i][j];
	    }
	  }
	}else{
	  //
	  //  Some mixing of surround with opposite direction
	  //
	  for(i=0;i<xn;i++){
	    for(j=0;j<yn;j++){
	      tvp = smep[i][j] + omix * smen[i][j];
	      tvn = smen[i][j] + omix * smep[i][j];
	      ss[fi ][i][j][k] = tvp;
	      ss[fio][i][j][k] = tvn;
	    }
	  }
	}

	free_2d_farray(smep,xn);
	free_2d_farray(smen,xn);
      }else{
	//
	//  Fill surround signal with '0' up until the delay time
	//
	for(i=0;i<xn;i++){
	  for(j=0;j<yn;j++){
	    ss[fi ][i][j][k] = 0.0;
	    ss[fio][i][j][k] = 0.0;
	  }
	}
      }
    }
  }

  if (aflag > 0){
    free_3d_farray(admep,xn,yn,tn);
    free_3d_farray(admen,xn,yn,tn);
    myfree(ttp);
    myfree(ttn);
  }

  free_2d_farray(mep,xn);
  free_2d_farray(men,xn);

  //
  //  Dump surround
  //
  if (surr->dumpi > -1){
    i = surr->dumpi;
    if (eyeflag == 0)
      write_3d_data_part("zzz.surr.3d",ss[i],0,xn,0,yn,0,tn,4,2,1);
    else
      // WYETH - this probably cannot happen if 'exit' below after left eye.
      write_3d_data_part("zzz.surr.re.3d",ss[i],0,xn,0,yn,0,tn,4,2,1);
    //write_3d_data_part("zzz.fs0.3d",fs[0]->r,1,xn,1,yn,1,tn,4,2,1);
    exit_error("MOD_ME_V5_V1SURR_SUM","Exiting after writing raw data");
  }

  mylog(mylogf,"    finished surround signal computation.\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_ME_V5_RAW_RESP                            */
/*                                                                           */
/*  Get the raw response traces from the filters, so that the filters can    */
/*  be re-run on another (the right eye) stimulus.                           */
/*                                                                           */
/*  Fill  mod_me_v5 -> v1_raw[]                                              */
/*                  -> v1_sum[]                                              */
/*  Or,                                                                      */
/*                  -> v1r_raw[]                                             */
/*                  -> v1r_sum[]                                             */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_raw_resp(eyeflag,stfi,v1s,x,y)
     int eyeflag;                  // 0-left, 1-right
     int stfi;                     // STF channel index
     struct mod_me_v1resp *v1s;    // Where to store responses
     int x,y;                      // [1..xn] , [1..yn]
{
  int i,k;
  int fi,fj,xn,yn,tn,fn,foff,re_flag;
  float t_mep,t_men,*tpe,*tpo,*tne,*tno,merpow,merthr,meroff;
  float *tr,*tr_opp,*lsum,***traw,**tme,*pstate;
  struct mod_me_fstruct **fs;    // Filter list
  struct mod_me_v5_struct *mv5;  // V5 params and responses
  //  Adaptation-related:
  float chalf,asum,*tra,adf;
  float frac,totw;

  xn  = mod_me_xn;
  yn  = mod_me_yn;
  tn  = mod_me_tn;
  mv5 = mod_me_v5;

  if (mv5->normo != NULL){  // Wyeth 2018Apr18
    merpow = onode_getpar_flt_dflt(mv5->normo,"me_raw_power",1.0);
    merthr = onode_getpar_flt_dflt(mv5->normo,"me_raw_thresh",0.0);
    meroff = onode_getpar_flt_dflt(mv5->normo,"me_raw_offset",0.0);
  }else{
    merpow = 1.0;
    merthr = 0.0;
    meroff = 0.0;
  }

  if (mod_me_ftop == NULL){
    fs  = mod_me_ff;
    fn  = mod_me_fn;
    re_flag = mod_me_fre;
  }else{
    fs = mod_me_ftop->fsl;   // New way
    fn = mod_me_ftop->n;
    re_flag = mod_me_ftop->flag_re;
  }

  foff = 0;  // Offset index in 'fs' array.
  if (eyeflag == 1){
    if (re_flag == 1){
      foff = fn/2;  // R. Eye uses second half of filter set
    }
  }

  traw = v1s->raw;
  tme  = v1s->me;

  if (mv5->v1_dflag == 0){

    //  Compute 'v1_raw' signals
    //    * Consider SQRT(.) ?
    //    * Consider non-opponent version ?
    for(fi=0;fi<mv5->n/2;fi++){

      fj = fi + mv5->n/2;  // index for opponent channel

      if (mv5->v1stf == NULL){
	tpe = fs[fi*4 + 0 + foff]->r[x][y];
	tpo = fs[fi*4 + 1 + foff]->r[x][y];
	tne = fs[fi*4 + 2 + foff]->r[x][y];
	tno = fs[fi*4 + 3 + foff]->r[x][y];
      }else{

	k = mod_me_fs_get_fi(eyeflag,stfi,fi,0); // 0 - EVEN phase
	//printf("stfi %d dir %d   k = %d",stfi,fi,k);
	tpe = fs[k]->r[x][y];
	k = mod_me_fs_get_fi(eyeflag,stfi,fi,1); // 1 - ODD phase
	tpo = fs[k]->r[x][y];

	// use 'fj' index for opponent channel

	k = mod_me_fs_get_fi(eyeflag,stfi,fj,0); // 0 - EVEN phase
	//printf("  OPPk = %d\n",k);
	tne = fs[k]->r[x][y];
	k = mod_me_fs_get_fi(eyeflag,stfi,fj,1); // 1 - ODD phase
	tno = fs[k]->r[x][y];
      }

      // Pointers to storage for responses
      tr     = traw[fi][0];  // Response for direction 'fi'
      tr_opp = traw[fj][0];  // Response for opposite direction

      for(i=0;i<tn;i++){
	// Compute Pref and Null ME signals
	t_mep = tpe[1+i]*tpe[1+i] + tpo[1+i]*tpo[1+i];
	t_men = tne[1+i]*tne[1+i] + tno[1+i]*tno[1+i];

	// This offset defaults to zero  2018Oct29
	t_mep += meroff;  // Add offset
	t_men += meroff;  // Add offset

	if (t_mep < merthr)
	  t_mep = 0.0;
	if (t_men < merthr)
	  t_men = 0.0;

	if (merpow == 0.5){ // Wyeth 2018Apr18
	  t_mep = sqrt(t_mep);
	  t_men = sqrt(t_men);
	}else if (merpow == 1.0){
	  ;
	}else{
	  t_mep = pow(t_mep,merpow);
	  t_men = pow(t_men,merpow);
	}

	if (mod_me_oppflag == 3){
	  //printf("____________MOD_ME_OPPFLAG == 3\n");
	  //
	  //  Let the V1 raw signal be *non* opponent
	  //
	  tr[i]     = t_mep;
	  tr_opp[i] = t_men;
	}else if (mod_me_oppflag == 0){
	  //
	  //  Let V1 raw signal be the half-wave rectified opponent sum
	  //
	  if (t_mep > t_men){
	    tr[i]     = t_mep - t_men;
	    tr_opp[i] = 0.0;
	  }else{
	    tr[i] = 0.0;
	    tr_opp[i] = t_men - t_mep;  // Opposite direction
	  }
	}else{
	  exit_error("MOD_ME_V5_RAW_RESP","Bad value for opp_flag");
	}
      }

      /*** WYETH DEBUG ****/
      /*  REMOVE
      if (fj == 6){
	char tstr[SLEN];

	if (max_of_farray(traw[fj],tn) < 0.01){
	  printf("x,y = %d %d  fi=%d  prep_count = %d  myID %d\n",x,y,fj,
		 prep_count,myid);

	  sprintf(tstr,"ID_%d__pcnt_%d",myid,prep_count);
	  append_farray_plot("zzz.dump.pl",tstr,tpe,tn,1);

	  printf("stim_min, max = %f %f\n",stim_min,stim_max);

	  mylog(mylogf,"    *** PROBLEM HERE ***.\n");


	  exit_error("HERE WYETH","max is < 0.01");
	}
      }*/
      /*** WYETH DEBUG ****/
    }

    // ************************************************************ BEGIN ADAPT
    //
    //  Adaptation
    //
    if (mv5->v1_adflag > 0){

      if (strcmp(mv5->v1_adtype,"s01")==0){

	chalf = mv5->v1_adp_a1;
	frac  = exp(-mod_me_tscale/mv5->v1_adp_t1);

	for(fi=0;fi<mv5->n;fi++){

	  //
	  //  WYETH BUG - this logic is WRONG
	  //    The surround has already adapted and change 'v1_adst_l' but
	  //    we don't want to start with that value here !!!!
	  //
	  if (mv5->v1srnd == NULL){
	    pstate = &(v1s->adst[fi][0]);
	  }else{
	    pstate = &(mv5->v1_adst_l[fi][x-1][y-1]);
	  }

	  mod_me_v5_adapt_v1_trace_s01(traw[fi][0],v1s->adrsp[fi][0],tn,
				       pstate,frac,chalf,
				       mv5->v1_adflag);
	}
      }else{
	exit_error("MOD_ME_V5_RAW_RESP","Unknown adaptation type");
      }
    }
    // ************************************************************** END ADAPT


    //
    //  Compute sum over all directions, for UNTUNED normalization signal
    //
    lsum = v1s->sum;

    for(i=0;i<tn;i++)  // Zero the summed signal
      lsum[i] = 0.0;

    for(fi=0;fi<mv5->n;fi++){
      if (mv5->v1_adflag > 0){
	tr = v1s->adrsp[fi][0];   // Adapted signal
      }else{
	tr = traw[fi][0];         // Raw signal
      }
      for(i=0;i<tn;i++){
	lsum[i] += tr[i];  // Sum all channels for norm below
      }
    }


  }else{
    //
    //  DISPARITY OPTION - Keep even and odd raw responses
    //
    //WYDISP

    // WYETH STF
    //   we would want to do the below for each FILTER[sf][tf] ... [dir][ph]
    //   
    // Filter order OLD
    //   for each dir:  pe po ne no  [R.E., if different, is 2nd half of list]
    //   

    for(fi=0;fi<mv5->n/2;fi++){  // For each direction around HALF circle

      fj = fi+mv5->n/2;  // index for opponent channel

      if (mv5->v1stf == NULL){
	tpe = fs[fi*4 + 0 + foff]->r[x][y];  // Get 4 quad opp filter response
	tpo = fs[fi*4 + 1 + foff]->r[x][y];  //   at the (x,y) spatial location
	tne = fs[fi*4 + 2 + foff]->r[x][y];
	tno = fs[fi*4 + 3 + foff]->r[x][y];
      }else{

	k = mod_me_fs_get_fi(eyeflag,stfi,fi,0); // 0 - EVEN phase
	tpe = fs[k]->r[x][y];
	k = mod_me_fs_get_fi(eyeflag,stfi,fi,1); // 1 - ODD phase
	tpo = fs[k]->r[x][y];

	// use 'fj' index for opponent channel

	k = mod_me_fs_get_fi(eyeflag,stfi,fj,0); // 0 - EVEN phase
	tne = fs[k]->r[x][y];
	k = mod_me_fs_get_fi(eyeflag,stfi,fj,1); // 1 - ODD phase
	tno = fs[k]->r[x][y];
      }

      for(i=0;i<tn;i++){
	traw[fi][0][i] = tpe[1+i];  // Save raw filter response (... ->raw)
	traw[fi][1][i] = tpo[1+i];
	traw[fj][0][i] = tne[1+i];
	traw[fj][1][i] = tno[1+i];

	t_mep = tpe[1+i]*tpe[1+i] + tpo[1+i]*tpo[1+i];
	t_men = tne[1+i]*tne[1+i] + tno[1+i]*tno[1+i];

	tme[fi][i] = sqrt(t_mep);  // Save ME signal (... ->me)
	tme[fj][i] = sqrt(t_men);

	// Sum of square (or abs) of all four channels (even/odd, pref/null)
	if (mv5->v1_dnorm == 0){
	  // Sum of Abs(...)
	  lsum[i] += fabs(tpe[1+i]) + fabs(tpo[1+i]) +
	    fabs(tne[1+i]) + fabs(tno[1+i]);
	}else if (mv5->v1_dnorm == 1){
	  // Sum of squares
	  lsum[i] += tpe[1+i]*tpe[1+i] + tpo[1+i]*tpo[1+i] +
	    tne[1+i]*tne[1+i] + tno[1+i]*tno[1+i];
	}else{
	  exit_error("MOD_ME_V5_RAW_RESP",
		     "Unknown disparity norm flag value\n");
	}
      }
    }
  }  // End disparity option
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_ME_V5_RAW_TOP                             */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_raw_top(eyeflag)
     int eyeflag;                  // 0-left, 1-right
{
  int i,j;
  int ci,xn,yn,tn,stfi;
  float ***nsum;
  struct mod_me_v5_struct *mv5;  // V5 params and responses
  struct mod_me_v1resp *v1s;     // Where to store responses
  struct mod_me_v1_input ***v1rr;

  xn  = mod_me_xn;
  yn  = mod_me_yn;
  tn  = mod_me_tn;
  mv5 = mod_me_v5;

  // ************************************
  stfi = 0;  // WYETH STF fix this  ***********************
  // *************
  if (mv5->v1stf == NULL){
    v1rr = mv5->v1resp; // old way
  }else{
    v1rr = mv5->v1stf->sgr[stfi]->bdgr; // new way
  }


  if (mod_me_mpt == NULL){
    //
    //  No spatial summation, compute responses for the central point only
    //

    if (eyeflag == 0)
      v1s = v1rr[0][0]->dl;
    else
      v1s = v1rr[0][0]->dr;

    mod_me_v5_raw_resp(eyeflag,stfi,v1s,xn/2,yn/2);
  }else{
    //
    //  Spatial summation - compute responses at flagged V1 positions
    //

    if (mv5->v1stf != NULL){
      //
      // new way
      //

      //  Zero the normalization sum for the appropriate eye.
      if (eyeflag == 0)
	nsum = mv5->v1stf->sum_norm;
      else
	nsum = mv5->v1stf->sum_norm_r;
      zero_3d_farray(nsum,xn,yn,tn);

      for(ci=0;ci<mv5->v1stf->n;ci++){
	v1rr = mv5->v1stf->sgr[ci]->bdgr;

	for(i=0;i<xn;i++){
	  for(j=0;j<yn;j++){
	    if (v1rr[i][j]->flag == 2){
	      if (eyeflag == 0)
		v1s = v1rr[i][j]->dl;
	      else
		v1s = v1rr[i][j]->dr;

	      mod_me_v5_raw_resp(eyeflag,ci,v1s,i+1,j+1);
	      add_to_farray(nsum[i][j],v1s->sum,tn);  // Add to norm. sum
	    }
	  }
	}
      }
    }else{
      //
      // old way
      //
      for(i=0;i<xn;i++){
	for(j=0;j<yn;j++){
	  if (v1rr[i][j]->flag == 2){
	    if (eyeflag == 0)
	      v1s = v1rr[i][j]->dl;
	    else
	      v1s = v1rr[i][j]->dr;
	    mod_me_v5_raw_resp(eyeflag,stfi,v1s,i+1,j+1);
	  }
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_ME_V5_NORMV1_C01                           */
/*                                                                           */
/*  Low level normalization.                                                 */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_normv1_c01(dr,dn,n,a1,a2,a4,c3,c3r,nsum,ds,surr_type,s1)
     float *dr;       // [n] raw data
     float *dn;       // [n] normalized data (returned values)
     int n;           // length of array
     float a1,a2,a4;  // Normalization constants
     float c3,c3r;    // Normalization constants
     float *nsum;     // [n] Normalization sum
     float *ds;       // [n] Surround signal
     int surr_type;   // 0-none, 1-subtr, 2-divide, ...
     float s1;        // Surround amplitude
{
  int i;
  float tt;

  if (surr_type == 0){ //  No surround suppression

    for(i=0;i<n;i++)
      dn[i] = dr[i] / (a1 * dr[i] + a2 * nsum[i] + c3 + a4);

  }else if (surr_type == 1){  //  Subtract from raw signal

    for(i=0;i<n;i++){
      tt = (dr[i] - s1*ds[i]);
      if (tt < 0.0)
	tt = 0.0;
      dn[i] = tt / (a1*dr[i] + a2*nsum[i] + c3 + a4);
    }

  }else if (surr_type == 2){  //  Use the surround in the denominator
    for(i=0;i<n;i++){
      dn[i] = dr[i] / (s1*ds[i] + a1*dr[i] + a2*nsum[i] + c3 + a4);
    }

  }else{
    exit_error("MOD_ME_V5_NORMV1_C01","Unknown value for 'surr_type'");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_ME_V5_NORM_V1_WORK                          */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_norm_v1_work(v1s,xi,yi,a1,a2,a4,c3,c3r,nos,norml,normr,
			    surr_type,s1)
     struct mod_me_v1_input *v1s;
     int xi,yi;            // (x,y) location (added for surround)
     float a1,a2,a4;       // Normalization constants
     float c3,c3r;         // Normalization constants
     float nos;            // Normalization output scale
     float *norml,*normr;  // Normalization sums
     int surr_type;        // 0-none, 1-subtr, 2-divide, ...
     float s1;             // Surround amplitude
{
  int i,j;
  int ndir,tn;
  float *tr1,*tr2,*trs,*nsum,tt;

  ndir = mod_me_v5->n;         // Number of direction channels
  tn   = mod_me_tn;

  trs  = NULL;  // Default value

  if (norml == NULL)
    nsum = v1s->dl->sum;  // old way
  else
    nsum = norml;         // new way

  for(i=0;i<ndir;i++){

    if (mod_me_v5->v1_adflag > 0){
      tr1 = v1s->dl->adrsp[i][0];  // adaptation
    }else{
      tr1 = v1s->dl->raw[i][0];  // Pointer to raw response
    }

    tr2 = v1s->dl->norm[i][0]; // Pointer to norm'd resp

    if (surr_type > 0)
      trs = mod_me_v5->v1srnd->sl[i][xi][yi]; // Pointer to surround signal

    //printf("_new__new__new__new__new__new__new__new__new_\n");
    mod_me_v5_normv1_c01(tr1,tr2,tn,a1,a2,a4,c3,c3r,nsum,trs,surr_type,s1);
    multiply_farray(tr2,tn,nos);  // Apply <norm> "output_scale"
  }

  if (mod_me_v5->binoc == 1){  // Compute a separate right eye signal

    if (normr == NULL)
      nsum = v1s->dr->sum;
    else
      nsum = normr;

    for(i=0;i<ndir;i++){

      if (mod_me_v5->v1_adflag > 0){
	tr1 = v1s->dr->adrsp[i][0];  // adaptation
      }else{
	tr1 = v1s->dr->raw[i][0];   // Ptr to raw resp
      }

      tr2 = v1s->dr->norm[i][0];  // Ptr to norm'd resp

      if (surr_type > 0)
	trs = mod_me_v5->v1srnd->sr[i][xi][yi]; // Pointer to surround signal

      mod_me_v5_normv1_c01(tr1,tr2,tn,a1,a2,a4,c3,c3r,nsum,trs,surr_type,s1);
      multiply_farray(tr2,tn,nos);  // Apply <norm> "output_scale"
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_ME_V5_NORM_V1_DISP                          */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_norm_v1_disp(v1s,a1,a2,a4,c3,c3r,nos,norml,normr)
     struct mod_me_v1_input *v1s;
     float a1,a2,a4;       // Normalization constants
     float c3,c3r;         // Normalization constants
     float nos;            // Normalization output scale
     float *norml,*normr;  // Normalization sums
{
  int i,j;
  int fi,fio,ndir,tn;
  float *tr1,*tr2,*tr3,*tr3off,*tr1r,*tr2r,*tr3r,*tr3roff,*tte,*tto;
  float *lsum,*lsumr,*mel,*mer,oppw,f,f1;
  struct mod_me_v5_struct *mv5;  // V5 params and responses

  mv5 = mod_me_v5;

  ndir = mv5->n;         // Number of direction channels
  tn   = mod_me_tn;

  //
  //
  //  (1) NORM: Fill in 'norm' values within 'mod_me_v1_input' struct
  //
  if (norml == NULL)
    lsum  = v1s->dl->sum;
  else
    lsum  = norml;

  if (normr == NULL)
    lsumr = v1s->dr->sum;
  else
    lsumr = normr;

  for(fi=0;fi<ndir;fi++){

    // These contain Sqrt(ME), for normalization
    mel = v1s->dl->me[fi];
    mer = v1s->dr->me[fi];

    for(j=0;j<mv5->v1_nc;j++){
      tr1  = v1s->dl->raw[fi][j];   // Ptr to raw resp
      tr2  = v1s->dl->norm[fi][j];  // Ptr to norm'd resp
      tr1r = v1s->dr->raw[fi][j];   // Ptr to raw resp
      tr2r = v1s->dr->norm[fi][j];  // Ptr to norm'd resp

      if (mv5->v1_dnorm == 0){
	for(i=0;i<tn;i++){  // does <norm> "output_scale"
	  tr2[i]  = nos*tr1[i]  / (a1*fabs(tr1[i])  + a2*lsum[i]  + c3  + a4);
	  tr2r[i] = nos*tr1r[i] / (a1*fabs(tr1r[i]) + a2*lsumr[i] + c3r + a4);
	}
      }else if (mv5->v1_dnorm == 1){
	//  Use energy instead of abs val.
	for(i=0;i<tn;i++){  // does <norm> "output_scale"
	  tr2[i]  = nos*tr1[i]  / (a1*mel[i] + a2*lsum[i]  + c3  + a4);
	  tr2r[i] = nos*tr1r[i] / (a1*mer[i] + a2*lsumr[i] + c3r + a4);
	}
      }
    }
  }


  //
  //  (2) OPP: Fill in 'opp' values within 'mod_me_v1_input' struct
  //  WYDISP - Opponency
  //
  oppw = mv5->v1_opp_w;
  for(fi=0;fi<ndir;fi++){

    fio = fi + ndir/2;         // Index of opponent channels
    if (fio >= ndir)
      fio -= ndir;

    // These contain Sqrt(ME), for opponent signal
    mel = v1s->dl->me[fio];  // ME signals (sqrt), for norm.
    mer = v1s->dr->me[fio];

    for(j=0;j<mv5->v1_nc;j++){
      tr1     = v1s->dl->norm[fi][j];   // Left eye norm'd
      tr2     = v1s->dl->norm[fio][j];  // Opponent for left
      tr3     = v1s->dl->opp[fi][j];    // opponent response
      tr3off  = v1s->dl->opp[fi][j+2];  // opponent response

      tr1r    = v1s->dr->norm[fi][j];  // Right eye norm'd
      tr2r    = v1s->dr->norm[fio][j]; // Opponent for right
      tr3r    = v1s->dr->opp[fi][j];   // compute opp. resp
      tr3roff = v1s->dr->opp[fi][j+2]; // compute opp. resp

      for(i=0;i<tn;i++){
	tr3[i]     =  tr1[i]  - oppw * mel[i];  // ON channel
	tr3r[i]    =  tr1r[i] - oppw * mer[i];
	tr3off[i]  = -tr1[i]  - oppw * mel[i];  // OFF channel
	tr3roff[i] = -tr1r[i] - oppw * mer[i];
      }

      if (mv5->v1_recto == 1){

	// ******* WYETH HERE:  Keep the other side after the opponent
	// ******* WYETH HERE:    Thus double the number of channels

	half_wave_rectify_farray(tr3,tn);
	half_wave_rectify_farray(tr3r,tn);
	half_wave_rectify_farray(tr3off,tn);
	half_wave_rectify_farray(tr3roff,tn);
      }
    }
  }


  //
  //  (3) BDE: Fill in 'bde' values within 'mod_me_v1_input' struct
  //  WYDISP - Binocular combination
  //    This leaves an (im)balanced binocular combination in the even [0]
  //    indices of mv5->v1_bde[fi][2] and mv5->v1r_bde[fi][2].
  //
  for(fi=0;fi<ndir;fi++){
    //
    //  WYETH *** Adding the squared signals is the same as taking the
    //            1/2-square of the signal and adding it to
    //            the 1/2-square of the negative signal.
    //
    f = mv5->v1_mix_f;
    f1 = 1.0 - f;

    //                           f * L-Even   +   f1 * R-Even
    tte = add_scale_farrays(v1s->dl->opp[fi][0],
			    v1s->dr->opp[fi][0],f,f1,tn);
    //                           f * L-Odd    +   f1 * R-Odd
    tto = add_scale_farrays(v1s->dl->opp[fi][1],
			    v1s->dr->opp[fi][1],f,f1,tn);

    if (v1s->dl->bde[fi][0] != NULL){
      myfree(v1s->dl->bde[fi][0]);  // Free the old storage
      myfree(v1s->dl->bde[fi][1]);
      myfree(v1s->dl->bde[fi][2]);
    }
    v1s->dl->bde[fi][0] = tte;
    v1s->dl->bde[fi][1] = tto;
    v1s->dl->bde[fi][2] = add_squared_farrays(tto,tte,tn);


    //
    //  Now do the same BDE computation, but with the "OFF" channel,
    //    which would otherwise have been lost because of the
    //    rectification following the Opponency above
    //   
    //
    //                           f * L-Even   +   f1 * R-Even
    tte = add_scale_farrays(v1s->dl->opp[fi][0+2],
			    v1s->dr->opp[fi][0+2],f,f1,tn);
    //                           f * L-Odd    +   f1 * R-Odd
    tto = add_scale_farrays(v1s->dl->opp[fi][1+2],
			    v1s->dr->opp[fi][1+2],f,f1,tn);

    if (v1s->dl->bde[fi][0] != NULL){
      myfree(v1s->dl->bde[fi][3]);  // Free the old storage
      myfree(v1s->dl->bde[fi][4]);
      myfree(v1s->dl->bde[fi][5]);
      myfree(v1s->dl->bde[fi][6]);
    }
    v1s->dl->bde[fi][3] = tte;
    v1s->dl->bde[fi][4] = tto;
    v1s->dl->bde[fi][5] = add_squared_farrays(tto,tte,tn);
    //  Scale by 1/2 to account for having twice as many channels (ON + OFF)
    v1s->dl->bde[fi][6] =
      add_scale_farrays(v1s->dl->bde[fi][2],
			v1s->dl->bde[fi][5],0.5,0.5,tn);

    //
    //  Now repeat above for 'v1s->dr'
    //

    //                           f1 * L-Even  +   f * R-Even
    tte = add_scale_farrays(v1s->dl->opp[fi][0],
			    v1s->dr->opp[fi][0],f1,f,tn);
    //                           f1 * L-Odd   +   f * R-Odd
    tto = add_scale_farrays(v1s->dl->opp[fi][1],
			    v1s->dr->opp[fi][1],f1,f,tn);

    if (v1s->dr->bde[fi][0] != NULL){
      myfree(v1s->dr->bde[fi][0]);  // Free the old storage
      myfree(v1s->dr->bde[fi][1]);
      myfree(v1s->dr->bde[fi][2]);
    }
    v1s->dr->bde[fi][0] = tte;
    v1s->dr->bde[fi][1] = tto;
    v1s->dr->bde[fi][2] = add_squared_farrays(tto,tte,tn);

    //
    //  Opponent "OFF" channel for Right eye
    //   
    //                           f1 * L-Even  +   f * R-Even
    tte = add_scale_farrays(v1s->dl->opp[fi][0+2],
			    v1s->dr->opp[fi][0+2],f1,f,tn);
    //                           f1 * L-Odd   +   f * R-Odd
    tto = add_scale_farrays(v1s->dl->opp[fi][1+2],
			    v1s->dr->opp[fi][1+2],f1,f,tn);

    if (v1s->dr->bde[fi][0] != NULL){
      myfree(v1s->dr->bde[fi][3]);  // Free the old storage
      myfree(v1s->dr->bde[fi][4]);
      myfree(v1s->dr->bde[fi][5]);
      myfree(v1s->dr->bde[fi][6]);
    }
    v1s->dr->bde[fi][3] = tte;
    v1s->dr->bde[fi][4] = tto;
    v1s->dr->bde[fi][5] = add_squared_farrays(tto,tte,tn);
    //  Scale by 1/2 to account for having twice as many channels (ON + OFF)
    v1s->dr->bde[fi][6] =
      add_scale_farrays(v1s->dr->bde[fi][2],
			v1s->dr->bde[fi][5],0.5,0.5,tn);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_ME_V5_NORM_V1                             */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_norm_v1(a1,a2,a4,c3,c3r,nos)
     float a1,a2,a4;       // Normalization constants
     float c3,c3r;         // Normalization constants
     float nos;            // Normalization output scale
{
  int i,j;
  int ci,stfi,surt;
  float *lsum,*rsum,sura;
  struct mod_me_v5_struct *mv5;  // V5 params and responses
  struct mod_me_v1_input *v1i;
  struct mod_me_v1_input ***v1rr;

  mv5 = mod_me_v5;

  // ************************************
  // ************************************
  stfi = 0;  // WYETH STF fix this  ***********************
  // *************
  // *************

  if (mv5->v1stf == NULL){
    v1rr = mv5->v1resp; // old way
  }else{
    v1rr = mv5->v1stf->sgr[stfi]->bdgr; // new way
  }

  if (mv5->v1srnd == NULL){
    surt = 0;
    sura = 0.0; // Wyeth 2018 Jun 28 (previously, no value assigned here.
  }else{
    sura = mv5->v1srnd->ampl;
    if (sura == 0.0)
      surt = 0;
    else if (strcmp(mv5->v1srnd->rule,"subtract")==0)
      surt = 1;
    else if (strcmp(mv5->v1srnd->rule,"divide")==0)
      surt = 2;
    else
      exit_error("MOD_ME_V5_NORM_V1","Unknown surround 'rule'.");
  }

  if (mod_me_mpt == NULL){  // No spatial summation
    if (mv5->v1_dflag == 0)
      mod_me_v5_norm_v1_work(v1rr[0][0],0,0,a1,a2,a4,c3,c3r,nos,NULL,NULL,
			     surt,sura);
    else
      mod_me_v5_norm_v1_disp(v1rr[0][0],a1,a2,a4,c3,c3r,nos,NULL,NULL);

  }else{  // Spatial summation - compute responses at flagged V1 positions
    if (mv5->v1stf != NULL){
      //
      // new way
      //
      for(ci=0;ci<mv5->v1stf->n;ci++){
	v1rr = mv5->v1stf->sgr[ci]->bdgr;

	for(i=0;i<mod_me_xn;i++){
	  for(j=0;j<mod_me_yn;j++){

	    lsum = mv5->v1stf->sum_norm[i][j];
	    rsum = mv5->v1stf->sum_norm_r[i][j];

	    if (v1rr[i][j]->flag == 2){
	      if (mv5->v1_dflag == 0)
		mod_me_v5_norm_v1_work(v1rr[i][j],i,j,a1,a2,a4,c3,c3r,nos,
				       lsum,rsum,surt,sura);
	      else
		// WYETH - does this work for new STF ???
		mod_me_v5_norm_v1_disp(v1rr[i][j],a1,a2,a4,c3,c3r,nos,
				       lsum,rsum);
	    }
	    // else - do not process this (i,j) location - data not needed
	  }
	}
      }
    }else{
      // Old way
      for(i=0;i<mod_me_xn;i++){
	for(j=0;j<mod_me_yn;j++){
	  if (v1rr[i][j]->flag == 2){
	    if (mv5->v1_dflag == 0)

	      // WYETH Heterog 2019 *** HERE is where we could vary
	      //   this 'a2' parameter.

	      mod_me_v5_norm_v1_work(v1rr[i][j],i,j,a1,a2,a4,c3,c3r,nos,
				     NULL,NULL,surt,sura);
	    else
	      mod_me_v5_norm_v1_disp(v1rr[i][j],a1,a2,a4,c3,c3r,nos,
				     NULL,NULL);
	  }
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_ME_V5_OPP_V1_WORK                           */
/*                                                                           */
/*   Fill:  v1s->dl->opp[0..ndir-1][0]                                       */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_opp_v1_work(stfi,v1s)
     int stfi;                        // STF index
     struct mod_me_v1_input *v1s;
{
  int i;
  int fi,fio,ndir,tn;
  float *tr1,*tr2,*tr3,oppw,tf;
  struct mod_me_v5_struct *mv5;  // V5 params and responses

  mv5  = mod_me_v5;
  ndir = mv5->n;         // Number of direction channels
  tn   = mod_me_tn;

  oppw = mv5->v1_opp_w;  // Regular opponency, may change below

  //
  //  Opponent sum for non-disparity case
  //
  for(fi=0;fi<ndir;fi++){

    if (stfi >= 0){
      if (mod_me_v5->v1stf->flag_config == 0){
	tf = mod_me_stf_get_tf_stf_dir(stfi,-1,1);
      }else if (mod_me_v5->v1stf->flag_config == 1){
	tf = mod_me_stf_get_tf_stf_dir(stfi,fi,1);
      }

      if (fabs(tf) < mv5->v1_opp_tmin){
	oppw = 0.0;              // set opponent weight to zero
      }else{
	oppw = mv5->v1_opp_w;    // Regular opponency calculation
      }
    }

    fio = fi + ndir/2;         // Index of opponent channels
    if (fio >= ndir)
      fio -= ndir;

    tr1 = v1s->dl->norm[fi][0];   // Left eye norm'd
    tr2 = v1s->dl->norm[fio][0];  // Opponent for left
    tr3 = v1s->dl->opp[fi][0];    // compute opp. resp
    for(i=0;i<tn;i++)
      tr3[i] = tr1[i] - oppw * tr2[i];

    if (mv5->v1_recto == 1)
      half_wave_rectify_farray(tr3,tn);

    if (mv5->binoc == 1){
      tr1 = v1s->dr->norm[fi][0];   // Right eye norm'd
      tr2 = v1s->dr->norm[fio][0];  // Opponent for right
      tr3 = v1s->dr->opp[fi][0];    // compute opp. resp
      for(i=0;i<tn;i++)
	tr3[i] = tr1[i] - oppw * tr2[i];

      if (mv5->v1_recto == 1)
	half_wave_rectify_farray(tr3,tn);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_ME_V5_OPP_V1                             */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_opp_v1()
{
  int i,j;
  int ci;
  struct mod_me_v5_struct *mv5;  // V5 params and responses
  struct mod_me_v1_input ***v1rr;

  mv5 = mod_me_v5;

  if (mv5->v1stf == NULL){
    v1rr = mv5->v1resp; // old way
  }else{
    // WYETH STF fix this  ***********************
    // WYETH STF fix this  ***********************
    v1rr = mv5->v1stf->sgr[0]->bdgr; // new way
  }

  if (mod_me_mpt == NULL){  // No spatial summation

    mod_me_v5_opp_v1_work(-1,v1rr[0][0]);

  }else{  // Spatial summation - compute responses at flagged V1 positions

    if (mv5->v1stf != NULL){
      //
      // new way
      //

      for(ci=0;ci<mv5->v1stf->n;ci++){

	//printf("_________%d_____  OPP_FLAG %d\n",ci,opp_flag);

	v1rr = mv5->v1stf->sgr[ci]->bdgr;

	for(i=0;i<mod_me_xn;i++){
	  for(j=0;j<mod_me_yn;j++){
	    if (v1rr[i][j]->flag == 2){
	      // WYETH vlink - send in 'ci' to compute TF in there
	      mod_me_v5_opp_v1_work(ci,v1rr[i][j]);
	    }
	  }
	}
      }
    }else{
      // old way
      for(i=0;i<mod_me_xn;i++){
	for(j=0;j<mod_me_yn;j++){
	  if (v1rr[i][j]->flag == 2){
	    mod_me_v5_opp_v1_work(-1,v1rr[i][j]);
	  }
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_ME_V5_MT_NL                              */
/*                                                                           */
/*  Apply a nonlinearity to the raw data 'd'.                                */
/*                                                                           */
/*****************************************************************************/
float *mod_me_v5_mt_nl(d,tn,nltype,a,b)
     float *d;
     int tn;
     char *nltype;   // "half_rect", "RMSM"
     float a,b;
{
  int i;
  float *t,x;

  t = (float *)myalloc(tn*sizeof(float));

  if (strcmp(nltype,"rmsm")==0){
    for(i=0;i<tn;i++)
      t[i] = a * exp(b * d[i]);  // Non-linear transformation
  }else if (strcmp(nltype,"half-rect")==0){

    for(i=0;i<tn;i++){
      x = d[i] - b;       // 'b' is threshold
      if (x > 0.0)
	t[i] = a * d[i];  // 'a' is slope
      else
	t[i] = 0.0;
    }
  }else if (strcmp(nltype,"none")==0){
    for(i=0;i<tn;i++)
      t[i] = d[i];  // no transform.
  }else{
    printf("  *** nltype:  %s\n",nltype);
    exit_error("MOD_ME_VT_MT_NL","Unknown MT nonlinearity type.");
  }

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_ME_V5_MT_WORK                             */
/*                                                                           */
/*  For a given v1 location, 'v1s', use the MT weights to integrate signals  */
/*  across direction channels within each stream.  These are like local      */
/*  'monocular' mt subunits.  Then, combine across both streams.             */
/*                                                                           */
/*  Fill:  v1s->dl->mt                                                       */
/*            ->dr->mt                                                       */
/*            ->mtb       the sum of the previous two, with RE bias weight   */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_mt_work(v1s,whack_flag,whack,w,wr,binw)
     struct mod_me_v1_input *v1s;
     int whack_flag;   // Tailby hack
     float whack;      // Tailby weight
     float *w;         // weights, Left eye, or Both
     float *wr;        // weights, Right eye
     float binw;       // Bias for right stream in L+R at MT level
{
  int i;
  int fi,ndir,tn;
  float *tr1,*tr2,*tr3,wt;
  struct mod_me_v5_struct *mv5;  // V5 params and responses

  mv5 = mod_me_v5;
  ndir = mv5->n;         // Number of direction channels
  tn   = mod_me_tn;

  //  Compute:       ...->dl->mt
  //    from either:  ->dl->bde (if 'dflag' is 1) or ->dl->opp
  //
  //  Compute:       ...->dr->mt
  //    from either:  ->dr->bde (if 'dflag' is 1) or ->dr->opp
  //
  //  Compute:       ...->mtb  from  ->dr->mt and ->dl->mt


  tr2 = v1s->dl->mt;  // Left MT stream
  for(i=0;i<tn;i++)   // reset MT signal to zero
    tr2[i] = 0.0;

  //
  //  Apply 'w' over direction channels in LEFT stream
  //
  for(fi=0;fi<ndir;fi++){
    // WYETH DIRMAP 2016 - could rotate the 'fi' index here as a function of
    //   position - we would need to send in another param to specify the
    //   rotation, or send in rotated weights.
    wt = w[fi];

    if ((whack_flag == 1) && (wt < 0.0)) // Tailby Hack - eventually remove?
      wt *= whack;   // Rescale negative weights

    //
    //  Pick signal, and add noise to that signal if needed
    //
    if (mv5->v1_dflag == 1){
      tr1 = v1s->dl->bde[fi][6];  // 'bde'
    }else{
      tr1 = v1s->dl->opp[fi][0];  // 'opp'
    }

    // Add noise to 'tr1'
    if (mv5->noise_v1_out != NULL)
      mod_me_noise_add(tr1,tn,mv5->noise_v1_out);

    for(i=0;i<tn;i++)
      tr2[i] += wt * tr1[i];   // linear sum in MT
  }

  tr2 = v1s->dr->mt;  // Right MT stream
  for(i=0;i<tn;i++)   // reset MT signal to zero
    tr2[i] = 0.0;

  //
  //  Apply 'wr' over direction channels in RIGHT stream
  //
  for(fi=0;fi<ndir;fi++){
    wt = wr[fi];

    if ((whack_flag == 1) && (wt < 0.0)) // Tailby Hack - eventually remove?
	wt *= whack;   // Rescale negative weights

    if (mv5->v1_dflag == 1)
      tr1 = v1s->dr->bde[fi][6];  // 'bde'
    else
      tr1 = v1s->dr->opp[fi][0];  // 'opp'

    if (mv5->noise_v1_out != NULL) // Add noise to 'tr1'
      mod_me_noise_add(tr1,tn,mv5->noise_v1_out);

    for(i=0;i<tn;i++)
      tr2[i] += wt * tr1[i];   // linear sum in MT
  }

  tr1 = v1s->mtb;     // compute the MT response combined across streams
  tr2 = v1s->dr->mt;  // Previously computed right eye signal
  tr3 = v1s->dl->mt;  // Previously computed left eye signal
  for(i=0;i<tn;i++){
    tr1[i] = tr3[i] + binw * tr2[i];  // Sum of left + right
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                MOD_ME_V5_MT                               */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_mt(whack_flag,whack,w,wr,r_bias)
     int whack_flag;   // Tailby hack
     float whack;      // Tailby weight
     float *w;         // Weights
     float *wr;        // Weights, R.Eye (optional)
     float r_bias;     // Bias for RE stream
{
  int i,j;
  int stfi;
  struct mod_me_v5_struct *mv5;  // V5 params and responses
  struct mod_me_v1_input ***v1rr;

  mv5 = mod_me_v5;


  // ************************************
  // ************************************
  stfi = 0;  // WYETH STF fix this  ***********************
  // *************
  // *************
  if (mv5->v1stf == NULL){
    v1rr = mv5->v1resp; // old way
  }else{
    v1rr = mv5->v1stf->sgr[stfi]->bdgr; // new way
  }


  if (mod_me_mpt == NULL){  // No spatial summation

    mod_me_v5_mt_work(v1rr[0][0],whack_flag,whack,w,wr,r_bias);

  }else{
    //
    // Spatial summation - compute responses at flagged V1 positions
    //

    // ***********************
    // ***********************
    // ***********************
    exit_error("MOD_ME_V5_MT","WYETH - THIS SHOULD NOT HAPPEN\n");
    // This code cannot run now, because 'mod_me_mpt' == NULL if this
    //    routine is called.

    for(i=0;i<mod_me_xn;i++){
      for(j=0;j<mod_me_yn;j++){
	if (v1rr[i][j]->flag == 2){  // 2: use both 'dl' and 'dr'
	  // WYETH DIRMAP 2016 - could rotate the 'w' here according to a
	  //   direction map (or could rotate it inside the routine)
	  mod_me_v5_mt_work(v1rr[i][j],whack_flag,whack,w,wr,r_bias);
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_MEU_GET_GROUP_INP                          */
/*                                                                           */
/*  Examine the synaptic inputs to MT cell 'c' and organize the inputs       */
/*  into appropriate groups for normalization and integration.               */
/*                                                                           */
/*****************************************************************************/
void mod_meu_get_group_inp(mtpop,c,k,mv5,gtype,rgl,rgr,rwl,rwr,rgcnt,rgn)
     struct mod_me_v5_pop *mtpop;   // where to write MT pop response
     struct pop_cell *c;
     int k;                         // z-index of MT cell
     struct mod_me_v5_struct *mv5;  // V5 params and responses
     char *gtype;   // 'none', 'stack_dir'
     float ****rgl; // [rng][rgcnt[]][tn] ptr to resp for each group, Left
     float ****rgr; // [rng][rgcnt[]][tn] ptr to resp for each group, Right
     float ***rwl;  // [rng][rgcnt[]] weight for each response, Left
     float ***rwr;  // [rng][rgcnt[]] weight for each response, Right
     int **rgcnt;   // [rng] Number of reponses per group
     int *rgn;      // Number of groups
{
  int i,j;
  int ii,stfi,li,xn,yn,wi,xi,yi,zi,nin,gi,ng,**sloc,nloc,*gcnt,*tgcnt;
  int wsyn_flag;
  float ***gpl,***gpr,**gwl,**gwr;
  char *wconf;
  struct pop_syn *s,*t;
  struct mod_me_v1_input *v1s,***v1rr;

  xn = mod_me_xn;
  yn = mod_me_yn;

  if (mv5->v1stf == NULL) // else 'v1rr' will be set below
    v1rr = mv5->v1resp; // old way

  //  Decide whether to use the stored synaptic weight
  wsyn_flag = 0;
  if (mtpop->wco != NULL){
    wconf = onode_getpar_chr_ptr_exit(mtpop->wco,"type");
    if ((strcmp(wconf,"QOC_2016")==0) ||
	(strcmp(wconf,"Rad_Conc")==0)){
      wsyn_flag = 1;
    }
  }


  if (strcmp(gtype,"stack_dir")==0){

// PAMELA SAYS - if this is on, tere is no effect.

    //
    //  Initialize a stack location array.  Value will indicate number of
    //  synaptic inputs at each location.
    //
    sloc = get_zero_2d_iarray(xn,yn);
    nloc = 0;  // Number of stacks
  }

  //
  //  Count the total number of inputs.  For stacks, make a location array.
  //
  nin = 0;
  s = c->in;
  while(s != NULL){  // For each input to this MT unit

    if (strcmp(gtype,"stack_dir")==0){
      //
      //  Build an array of stack location
      //
      xi = s->pre->layx;  // (x,y) Coords within V1 array
      yi = s->pre->layy;
      if ((xi < 0)||(xi >= xn)||(yi < 0)||(yi >= yn))
	exit_error("MOD_MEU_GET_GROUP_INP","Coordinate value error");

      if (sloc[xi][yi] == 0)
	nloc += 1;  // This is a new stack location
      sloc[xi][yi] += 1;  // Count number of inputs at this location
    }
    nin += 1;

    s = s->post_next;
  }
  //printf("found %d inputs\n",nin);
  //if (strcmp(gtype,"stack_dir")==0)
  //printf("  There are %d stack locations\n",nloc);

  //
  //  Set up the structure of the groups, and allocate storage
  //
  if (strcmp(gtype,"none")==0){
    ng = nin;
    gcnt = get_const_iarray(ng,1);
  }else if (strcmp(gtype,"stack_dir")==0){
    ng = nloc;    //  Set 'ng' to the number of stacks
    gcnt = (int *)myalloc(ng*sizeof(int));

    gi = 0;
    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){
	if (sloc[i][j] > 0){
	  gcnt[gi] = sloc[i][j];  // Input count for this group
	  sloc[i][j] = gi;        // Group index
	  gi += 1;
	}else{
	  sloc[i][j] = -1;        // Flag indicates no group here
	}
      }
    }
    if (gi != ng)
      exit_error("MOD_ME_UTIL"," gi != ng, this should not happen");

    tgcnt = get_zero_iarray(ng);
  }

  //
  //  Create storage for response pointers and weights for all groups
  //
  gpl = (float ***)myalloc(ng*sizeof(float **));
  gpr = (float ***)myalloc(ng*sizeof(float **));
  gwl = (float **)myalloc(ng*sizeof(float *));
  gwr = (float **)myalloc(ng*sizeof(float *));
  for(i=0;i<ng;i++){
    gpl[i] = (float **)myalloc(gcnt[i]*sizeof(float *));
    gpr[i] = (float **)myalloc(gcnt[i]*sizeof(float *));
    gwl[i] = (float *)myalloc(gcnt[i]*sizeof(float));
    gwr[i] = (float *)myalloc(gcnt[i]*sizeof(float));
  }

  //
  //  Fill in response pointers and weights for each group
  //
  gi = 0;
  ii = 0;
  s = c->in;
  while(s != NULL){  // For each input to this MT unit

    if (mv5->v1stf != NULL){  // new way
      li = popu_get_layer_index_by_name(mod_me_mpt,s->pre->pl->name);
      stfi = li - 1;  // Assuming layer index mapping to STF index
      v1rr = mv5->v1stf->sgr[stfi]->bdgr;
    }

    xi = s->pre->layx;  // (x,y) Coords within V1 array
    yi = s->pre->layy;
    zi = s->pre->layz;  // Direction index

    if (strcmp(gtype,"stack_dir")==0){
      gi = sloc[xi][yi];   // Look up group index for this spatial location
      ii = tgcnt[gi];  // Get number of responses at this location so far
      tgcnt[gi] += 1;
    }

    //printf("gi = %d  ii = %d\n",gi,ii);

    v1s = v1rr[xi][yi];  // All responses at this location
    if (v1s->flag == 0)
      mylog_exit(mylogf,"MOD_ME_V5_POP  input not saved.");

    if (mv5->v1_dflag == 1){
      gpl[gi][ii] = v1s->dl->bde[zi][6];  // 'bde'
      gpr[gi][ii] = v1s->dr->bde[zi][6];
    }else{
      gpl[gi][ii] = v1s->dl->opp[zi][0];  // 'opp'
      gpr[gi][ii] = v1s->dr->opp[zi][0];

      // WYETH varmoo_bug -- NOTE, the 'dr->' above may be UNINITIALIZED
      //   and not used.  A condition on 'mv5->binoc == 1' can be used here
      //   to prevent computing 'gpr' if desired.

/***
      // WYETH varmoo_bug ???
      if (v1s->dl->opp[zi][0][0] > 1000.0){
	printf("v1s->dl->opp[zi][0][0] = %f\n",v1s->dl->opp[zi][0][0]);
	exit_error("WWWWWWWWW","Left Too big\n");
      }
      if (v1s->dr->opp[zi][0][0] > 1000.0){
	printf("v1s->dr->opp[zi][0][0] = %f\n",v1s->dr->opp[zi][0][0]);
	exit_error("WWWWWWWWW","Right Too big\n");
      }
***/

    }


    if (mtpop->zdirflag == 1){ // WYETH DIRMAP 2016 - rotate 'zi' for a dirmap
      // **** WYETH should move out of here and set syn w in prep routines ****
      if (wsyn_flag == 1) // *** WYETH we should move this out of here
	mylog_exit(mylogf,"MOD_ME_V5_POP  dir-map not imp'd w/ syn w.");

      wi = zi - k;
      if (wi < 0)
	wi += mv5->n;
    }else
      wi = zi;

    if ((wi < 0) || (wi >= mv5->n)) // WYETH - SHOULD NOT HAPPEN - REMOVE?
      mylog_exit(mylogf,"MOD_ME_V5_POP  'wi' out of range.");

    if (wsyn_flag == 1){
      gwl[gi][ii] = s->w;    // Use the stored synaptic weight
      gwr[gi][ii] = s->w;    // **** WYETH NOTE, same w for R.E. as L.E.
    }else{
      gwl[gi][ii] = mtpop->w[wi];    // Choose weight for this dir channel
      gwr[gi][ii] = mtpop->w_r[wi];  // Choose weight for this dir channel
    }

    if (strcmp(gtype,"none")==0)
      gi += 1;

    s = s->post_next;
  } // End of loop over inputs to this MT unit


  //
  //  Clean up
  //
  if (strcmp(gtype,"stack_dir")==0){
    free_2d_iarray(sloc,xn);
    myfree(tgcnt);
  }

  *rgl = gpl;
  *rgr = gpr;
  *rwl = gwl;
  *rwr = gwr;
  *rgcnt = gcnt;
  *rgn = ng;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_MEU_POP_MT_INTEGRATE                        */
/*                                                                           */
/*****************************************************************************/
void mod_meu_pop_mt_integrate(gtype,within,gl,gr,wl,wr,gcnt,gn,tn,normflag,
			      mtpow,sl,sr)
     char *gtype;  // 'none', 'stack_dir'
     char *within; // 'sum', 'avg', ...
     float ***gl;  // [gn][gcnt[]][tn] ptr to resp for each group, Left
     float ***gr;  // [gn][gcnt[]][tn] ptr to resp for each group, Right
     float **wl;   // [gn][gcnt[]] weight for each response, Left
     float **wr;   // [gn][gcnt[]] weight for each response, Right
     int *gcnt;    // [gn] Number of reponses per group
     int gn;       // Number of groups
     int tn;       // Length of response in time units
     int normflag; // 0-no norm; 1-power norm; 2-No_Subunit_power_1st
     float mtpow;  // Power for normalization
     float *sl;    // Response, left stream; Should be all zeros coming in
     float *sr;    // Response, right stream; Should be all zeros coming in
{
  int i,j;
  int ti,tsign,nn;
  float posval,*tsuml,*tsumr,*tl,*tr,twl,twr,const_n;

  //
  //  *** WYETH varmoo_bug - NOTE, RE signals are not initialized, and should
  //  ***                 not be calculated unless 'mv5->binoc == 1'
  //

  tsuml = get_farray(tn);
  tsumr = get_farray(tn);

  //
  //  Apply the linear weights within the group, then apply the power
  //
  if (normflag < 2){
    for(i=0;i<gn;i++){  // For each group

      for(ti=0;ti<tn;ti++){  // Zero the sum for this group
	tsuml[ti] = 0.0;
	tsumr[ti] = 0.0;
      }

      nn = gcnt[i];  // Number of response traces in this group
      for(j=0;j<nn;j++){  // Combine linearly within the group
	tl  = gl[i][j];
	tr  = gr[i][j];
	twl = wl[i][j];
	twr = wr[i][j];
	for(ti=0;ti<tn;ti++){
	  tsuml[ti] += twl * tl[ti];   // linear sum
	  tsumr[ti] += twr * tr[ti];   // linear sum
	}
      }
      if ((strcmp(within,"avg")==0) && (nn > 1)){
	//
	//  Compute the average within the subunit
	//
	for(ti=0;ti<tn;ti++){
	  tsuml[ti] /= (float)nn;  // Change sum to average
	  tsumr[ti] /= (float)nn;
	}
      }

      //
      //  Add to non-linear total across groups
      //
      if (normflag == 0){  // No normalization
	for(ti=0;ti<tn;ti++){
	  sl[ti] += tsuml[ti];
	  sr[ti] += tsumr[ti];
	}
      }else if (normflag == 1){  // power normalization
	for(ti=0;ti<tn;ti++){
	  posval = tsuml[ti];
	  tsign = 1;
	  if (posval < 0.0){
	    tsign = -1;
	    posval = -posval;
	  }
	  sl[ti] += tsign * pow(posval,mtpow);

	  posval = tsumr[ti];
	  tsign = 1;
	  if (posval < 0.0){
	    tsign = -1;
	    posval = -posval;
	  }
	  sr[ti] += tsign * pow(posval,mtpow);
	}
      }
    }
  }else if (normflag == 2){
    //
    //  When there are no subunits, there is 1 chan per group.
    //
    //  Here we apply the power first, then multiply by the weights.
    //
    for(i=0;i<gn;i++){  // For each channel (group of 1)

      if (gcnt[i] != 1) // Must be one channel per group
	mylog_exit(mylogf,"MOD_MEU_POP_MT_INTEGRATE  'gcnt' != 1");

      tl = gl[i][0];
      tr = gr[i][0];

      twl = wl[i][0];
      twr = wr[i][0];

      for(ti=0;ti<tn;ti++){
	posval = tl[ti];
	tsign = 1;
	if (posval < 0.0){
	  tsign = -1;
	  posval = -posval;
	}
	sl[ti] += twl * tsign * pow(posval,mtpow);  // Power, then weight

	posval = tr[ti];
	tsign = 1;
	if (posval < 0.0){
	  tsign = -1;
	  posval = -posval;
	}
	sr[ti] += twr * tsign * pow(posval,mtpow);
      }
    }
  }else{
    mylog_exit(mylogf,"MOD_MEU_POP_MT_INTEGRATE  Unknown 'normflag' value");
  }

  //
  //  Normalize the sum
  //
  if (normflag == 0){
    //
    // Simple linear average
    //
    for(ti=0;ti<tn;ti++){
      sl[ti] /= (float)gn;
      sr[ti] /= (float)gn;
    }
  }else if (normflag == 1){  // power normalization
    //
    //  Now raise the sum to the 1/pow
    //

    const_n = pow((float)gn,1.0/mtpow);

    for(ti=0;ti<tn;ti++){
      posval = sl[ti];
      tsign = 1;
      if (posval < 0.0){
	tsign = -1;
	posval = -posval;
      }
      sl[ti] = tsign * pow(posval,1.0/mtpow) / const_n;

      posval = sr[ti];
      tsign = 1;
      if (posval < 0.0){
	tsign = -1;
	posval = -posval;
      }
      sr[ti] = tsign * pow(posval,1.0/mtpow) / const_n;
    }
  }

  myfree(tsuml);
  myfree(tsumr);
}
/**************************************-**************************************/
/*                                                                           */
/*                               MOD_ME_V5_POP                               */
/*                                                                           */
/*  Compute responses for all MT units.  New way - does not assume the       */
/*  V1 inputs are 'stacked' at each x,y position.                            */
/*                                                                           */
/*  Fills:                                                                   */
/*    mtpop->ssum                                                            */
/*         ->nl                                                              */
/*                                                                           */
/*  Alternative for subunits:     WYETH HERE                                 */
/*    The only storage of the synaptic weights is in                         */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_pop(pl,mtpop)
     struct pop_layer *pl;          // MT layer
     struct mod_me_v5_pop *mtpop;   // where to write MT pop response
{
  int i,j,k;
  int ti,tn,mt_pow_flag,*gcnt,gn,li,mxn,myn,mzn;
  float *trr,*trl,mta,mtb,mtpow,***gl,***gr,**wl,**wr,*subsig,w;
  char *gtype,*within;
  struct mod_me_v5_struct *mv5;  // V5 params and responses
  struct pop_cell *c;
  struct mod_me_v5_pop *mtsub;   // MT pop acting as subunits

  mylog(mylogf,"  MOD_ME_V5_POP\n");

  mv5 = mod_me_v5;
  tn  = mod_me_tn;

  //
  //  MT - Inhib-Subunit
  //
  if (mtpop->mtsub_name != NULL){

    li = popu_get_layer_index_by_name(mod_me_mpt,mtpop->mtsub_name);
    if (li < 0){
      printf("  %s\n",mtpop->mtsub_name);
      mylog_exit(mylogf,"MOD_ME_V5_POP  Unknown population name.");
    }
    mtsub = mv5->mtp[li - mv5->mtli0];

    subsig = get_zero_farray(tn);

    mxn = mod_me_mpt->lay[li]->xn;
    myn = mod_me_mpt->lay[li]->yn;
    mzn = mod_me_mpt->lay[li]->zn;

    for(i=0;i<mxn;i++){
      for(j=0;j<myn;j++){
	for(k=0;k<mzn;k++){
	  w = mtpop->mtsub_w[k];
	  for(ti=0;ti<tn;ti++){
	    subsig[ti] += w * mtsub->nl[i][j][k][ti];
	  }
	}
      }
    }
    multiply_farray(subsig,tn,1.0/(float)(mxn*myn));  // Average over space

    //exit_error("WYETH HERE mod_me_v5_pop","Under development");
  }else{
    subsig = NULL;
  }


  //
  //  MT-POW-NORM  Prepare for possible normalization
  //
  mt_pow_flag = 0;
  gtype  = mtpop->inorm_subu;
  within = mtpop->inorm_subwi;

  if (mtpop->inorm_pow != 1.0){
    if (mtpop->inorm_pow1st == 1)
      mt_pow_flag = 2;
    else
      mt_pow_flag = 1;
    mtpow = mtpop->inorm_pow;
    sprintf(ggstr,"   Using MT normalization with power = %f\n",mtpow);
    mylog(mylogf,ggstr);
    sprintf(ggstr,"   Subunit type:  %s\n",gtype);
    mylog(mylogf,ggstr);
  }

  mta = mtpop->nl_a;
  mtb = mtpop->nl_b;

  trr = (float *)myalloc(tn*sizeof(float *));
  trl = (float *)myalloc(tn*sizeof(float *));

  for(i=0;i<pl->xn;i++){     // For each MT cell
    for(j=0;j<pl->yn;j++){
      for(k=0;k<pl->zn;k++){

	c = &(pl->c[i][j][k]);  // Pointer to the MT unit

	for(ti=0;ti<tn;ti++){
	  trr[ti] = 0.0;     // Zero the total left MT sum
	  trl[ti] = 0.0;     // Zero the total right MT sum
	}

	//
	//  Get pointers to all input signals to be integrated, and weights
	//  for each input.  Inputs are organized in 'ng' groups.
	//
	mod_meu_get_group_inp(mtpop,c,k,mv5,gtype,&gl,&gr,&wl,&wr,&gcnt,&gn);

	//
	//  Integrate input signals within each group, then across groups.
	//
	mod_meu_pop_mt_integrate(gtype,within,gl,gr,wl,wr,gcnt,gn,tn,
				 mt_pow_flag,mtpow,trl,trr);

/***
// WYETH varmoo_bug ???
if (trl[0] > 1000.0){
  printf("trl[0] = %f\n",trl[0]);
  exit_error("HERE WYETH","trl TOO BIG");
}
if (trr[0] > 1000.0){
  printf("trr[0] = %f\n",trr[0]);
  exit_error("HERE WYETH","trr TOO BIG");
}
***/

	//
	//  Free storage
	//
	myfree(gcnt);
	free_2d_farray(wl,gn);
	free_2d_farray(wr,gn);
	free_2d_farray(gl,gn);  // Don't free last dimension
	free_2d_farray(gr,gn);  // Don't free last dimension


	//
	//  Combine left and right stream signals
	//

	// WYETH varmoo_bug June 10, 2018 DOES THIS FIX IT ???
	if (mod_me_v5->binoc == 1){  // Use the right eye signal
	  multiply_farray(trr,tn,mtpop->r_bias);  // Imbalance for right stream
	  mtpop->ssum[i][j][k] = add_farrays(trl,trr,tn);  // Add L+R streams
	}else{
	  mtpop->ssum[i][j][k] = copy_farray(trl,tn);  // Use only L stream
	}

	if (subsig != NULL){
	  //
	  //  Add the subunit signal
	  //  WYETH - this should allow an MT Component cell to apply
	  //          inhibition to an MT Pattern cell
	  //
	  for(ti=0;ti<tn;ti++){
	    mtpop->ssum[i][j][k][ti] += subsig[ti];
	  }
	}

	mtpop->nl[i][j][k] = mod_me_v5_mt_nl(mtpop->ssum[i][j][k],tn,
					     mtpop->nl_type,mta,mtb);
      }
    }
  }

  myfree(trr);
  myfree(trl);

  if (subsig != NULL)
    myfree(subsig);
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_ME_V5_MIX_WORKER                           */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_mix_worker(s1,s2,n,tn,f)
     float ***s1,***s2;  // [n][nc][tn] Signals to mix  //WYDISP
     int n;              // Number of signals
     int tn;             // Length of signals
     float f;            // Fraction of main stream, 0.5 to 1.0
{
  int i,j;
  float *t1,*t2,f1;

  f1 = 1.0 - f;

  for(i=0;i<n;i++){

    t1 = copy_farray(s1[i][0],tn);
    t2 = copy_farray(s2[i][0],tn);

    for(j=0;j<tn;j++){
      s1[i][0][j]  = f*t1[j] + f1*t2[j];
      s2[i][0][j]  = f*t2[j] + f1*t1[j];
    }

    myfree(t1);
    myfree(t2);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_ME_V5_MIX_V1                             */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_mix_v1(mylogf,stage_str)
     char *mylogf;
     char *stage_str;  // This string indicates where this was called from
{
  int i,j;
  int tn,ndir,xn,yn,stfi,ci;
  struct mod_me_v5_struct *mv5;  // V5 params and responses
  struct mod_me_v1_input ***v1rr;

  mv5 = mod_me_v5;
  xn  = mod_me_xn;
  yn  = mod_me_yn;

  if (mv5->v1stf == NULL){
    v1rr = mv5->v1resp; // old way
  }else{
    if (mv5->v1stf->n > 1){
      exit_error("MOD_ME_V5_MIX_V1",
		 "Mixing not implemented for multiple STF channels");
    }
    stfi = 0;  // *** WYETH mutliple STF channels not implemented yet
    v1rr = mv5->v1stf->sgr[stfi]->bdgr; // new way
  }

  //
  //  Return right away if we are not going to do mixing
  //
  if ((strcmp(stage_str,mv5->v1_mix_s)!=0) ||  // Not the mixing stage
      (mv5->binoc != 1) ||                     // Not a binoc model
      (mv5->v1_mix_f >= 1.0)){                 // No mixing requested
    return;
  }
  if (mod_me_ftop != NULL)
    exit_error("MOD_ME_V5_MIX_V1","New system needed here.");

  ndir = mv5->n;   // Number of direction channels
  tn = mod_me_tn;

  if (mod_me_mpt == NULL){
    //
    //  No spatial summation - Center response only
    //

    if (mv5->v1stf->n > 1)
      exit_error("MOD_ME_V5_MIX_V1","Not implemented for multi STF");


    if (strcmp(stage_str,"linear")==0){  // Mix the raw filter outputs
      mylog(mylogf,"    *** Mixing linear V1 signals\n");
      mod_me_v5_mix_worker(v1rr[0][0]->dl->raw,
			   v1rr[0][0]->dr->raw,ndir,tn,mv5->v1_mix_f);

      /*** WYETH FIX ****/
      /*** WYETH FIX ****/
      printf("*** WARNING:  the normalization sum has not been mixed\n");

    }else if (strcmp(stage_str,"normalized")==0){
      mylog(mylogf,"    *** Mixing normalized V1 signals\n");
      mod_me_v5_mix_worker(v1rr[0][0]->dl->norm,
			   v1rr[0][0]->dr->norm,ndir,tn,mv5->v1_mix_f);
    }else if (strcmp(stage_str,"opponent")==0){
      mylog(mylogf,"    *** Mixing opponent V1 signals\n");
      mod_me_v5_mix_worker(v1rr[0][0]->dl->opp,
			   v1rr[0][0]->dr->opp,ndir,tn,mv5->v1_mix_f);
    }else{
      mylog_exit(mylogf,"MOD_ME_V5_MIX_V1  Unknown mixing stage");
    }
  }else{
    //
    //  Spatial summation - mix all inputs
    //

    for(ci=0;ci<mv5->v1stf->n;ci++){
      v1rr = mv5->v1stf->sgr[ci]->bdgr;

      for(i=0;i<xn;i++){
	for(j=0;j<yn;j++){
	  if (v1rr[i][j]->flag == 2){

	    if (strcmp(stage_str,"linear")==0){  // Mix the raw filter outputs
	      if (i+j == 0)
		mylog(mylogf,"    *** Mixing linear V1 signals\n");
	      mod_me_v5_mix_worker(v1rr[i][j]->dl->raw,
				   v1rr[i][j]->dr->raw,
				   ndir,tn,mv5->v1_mix_f);

	      /*** WYETH FIX ****/
	      /*** WYETH FIX ****/
	      if (i+j == 0)
		printf("*** WARNING:  normalization sum has not been mixed\n");

	    }else if (strcmp(stage_str,"normalized")==0){
	      if (i+j == 0)
		mylog(mylogf,"    *** Mixing normalized V1 signals\n");
	      mod_me_v5_mix_worker(v1rr[i][j]->dl->norm,
				   v1rr[i][j]->dr->norm,
				   ndir,tn,mv5->v1_mix_f);
	    }else if (strcmp(stage_str,"opponent")==0){
	      if (i+j == 0)
		mylog(mylogf,"    *** Mixing opponent V1 signals\n");
	      mod_me_v5_mix_worker(v1rr[i][j]->dl->opp,
				   v1rr[i][j]->dr->opp,
				   ndir,tn,mv5->v1_mix_f);
	    }else{
	      mylog_exit(mylogf,"MOD_ME_V5_MIX_V1  Unknown mixing stage");
	    }
	  }
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          MODEL_ME_V5_GET_RESPONSE                         */
/*                                                                           */
/*  Calling sequence:                                                        */
/*    1.  model_me_01_prep                                                   */
/*    2.  model_me_01_prep_ft                                                */
/*    3.  this routine                                                       */
/*        ...                                                                */
/*    N.  model_me_01_done                                                   */
/*                                                                           */
/*****************************************************************************/
void model_me_v5_get_response(m,s,r,k)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response data
     int k;                      // repeat number, k=0 means new stim
{
  int i;
  int xn,yn,tn,ndir,fn,n,whack_flag,monoc_gate,*lseedlist;
  int re_flag,stfi,nc,rawflag;
  float stimpar_l,stimpar_lr,whack,rmsml;
  float a1,a2,c3,c3r,a4,nos,mta,mtb;
  char tstr[SLEN];
  struct mod_me_fstruct **fs;    // Filter list
  struct mod_me_v5_struct *mv5;  // V5 params and responses
  struct pop_layer *pl;
  struct mod_me_v5_pop *mtpop;

  mod_me_tsi = r->tsi;

  xn   = mod_me_xn;
  yn   = mod_me_yn;
  tn   = mod_me_tn;
  mv5  = mod_me_v5;
  ndir = mv5->n;           // Number of direction channels

  if (mod_me_ftop == NULL){
    fs  = mod_me_ff;         // OLD WAY - REMOVE when not used
    fn  = mod_me_fn;
    re_flag = mod_me_fre;
  }else{
    fs = mod_me_ftop->fsl;   // New way
    fn = mod_me_ftop->n;
    re_flag = mod_me_ftop->flag_re;
  }

  //
  //  Initialize noise
  //
  if (mv5->noise_v1_out != NULL)
    mv5->noise_v1_out->sflag = 1;  // To indicate that seed must be initialized

  rawflag = 0;  // Use default flow for raw response / early adapt computation
  if (mv5->v1_adflag > 0){
    if (strcmp(mv5->v1_adtype,"RPHomeo")==0)
      rawflag = 1;  // Use alternative flow for raw response /early adapt
  }

  //
  //  Compute new response for each filter, only if needed
  //
  if ((k==0) || (fs[0]->r == NULL) || (s->stimno != mod_me_stimno_prev)){

    mod_me_fs_free_r(fs,fn);  // Free old response, if it exists

    if (re_flag == 1)  // There are separate filters for L and R eyes.
      mod_me_fs_flag_eye(fs,fn,0);  // 0-Left eye; Flag filters to use

    // Fill in the 3D responses in fs->r, and then
    //   store relevant response traces, so 3D response arrays can be freed.

    mod_me_compute_response_from_ft(s->d,fs,fn,xn,yn,tn);   // Get new resp
    // ADAPT HERE ??
    //...v1_adapt();
    if (rawflag == 1){
      mod_me_v5_v1_simul(0);  // RPHomeo
    }else{
      mod_me_v5_v1surr_sum(0);  // 0-left eye;  Surround
      mod_me_v5_raw_top(0);     // 0-left eye;
    }


    if (s->d_r != NULL){

      mv5->binoc = 1;
      mylog(mylogf,"    right eye\n");

      if (re_flag == 0){           // One set of filters, shared by both eyes
	mod_me_fs_free_r(fs,fn);     // Free old response, if it exists
      }else{                       // Separate filter sets for the two eyes
	mod_me_fs_flag_eye(fs,fn,1);  // 1-Right eye; Flag filters to use
      }

      mod_me_compute_response_from_ft(s->d_r,fs,fn,xn,yn,tn); // Right eye
      if (rawflag == 1){
	mod_me_v5_v1_simul(1);  // RPHomeo
      }else{
	mod_me_v5_v1surr_sum(1);  // 1-right eye
	mod_me_v5_raw_top(1);     // 1-right eye
      }
    }else{
      mv5->binoc = 0;
      if (re_flag == 1)
	exit_error("MODEL_ME_V5_GET_RESPONSE",
		   "Stimulus is monocular, but there are RE filters.");
    }
  }


  mod_me_v5_mix_v1(mylogf,"linear");

  //
  //  Set constant parameters for 'v1_norm' signals
  //
  // *** NOTE:  'stimmod' parameter is not used in stimulus generation ***
  //
  if (mv5->normo == NULL){
    //
    // old way  (no <norm> object)
    //
    stimpar_l  = param_getf_dflt(s->ppl,"stimmod_msc",mv5->L);
    if (stimpar_l < 0.0)
      stimpar_l = mv5->L;
    sprintf(ggstr,"   STIMPAR: rmsm_L = %f\n",stimpar_l);
    mylog(mylogf,ggstr);
    if (stimpar_l < 0.0)
      mylog_exit(mylogf,
	 "MODEL_ME_V5_GET_RESPONSE  Either rmsm_L or stimpar_msc < 0");

    a1 = mv5->alph1;
    a2 = mv5->alph2 / (float)ndir;
    c3 = mv5->alph3 * stimpar_l;  // L is "mean squared contast of stimulus"
    a4 = mv5->alph4;
    nos = 1.0;  // <norm> "output_scale"

    if (mv5->binoc == 1){  // Compute a separate right eye signal
      stimpar_lr = param_getf_dflt(s->ppl,"stimmod_msc_r",stimpar_l);
      if (stimpar_lr < 0.0)
	stimpar_lr = stimpar_l;
      c3r = mv5->alph3 * stimpar_lr;  // L is "mean squared contast of stimulus"
    }

  }else{
    //
    // new way, using <norm> object [Pamela uses this for the STF models]
    //
    // Note, L is "mean squared contast of stimulus"

    if (onode_test_chr(mv5->normo,"type","rmsm")==1){

      rmsml = 0.0;  // Default value
      c3 = onode_getpar_flt_dflt(mv5->normo,"alpha_3",0.0);
      if (c3 != 0.0){
	rmsml = onode_getpar_flt_exit(mv5->normo,"rmsm_L");
	c3 *= rmsml;
      }

      //rmsml = onode_getpar_flt_exit(mv5->normo,"rmsm_L");

      stimpar_l  = param_getf_dflt(s->ppl,"stimmod_msc",rmsml);
      if (stimpar_l < 0.0)
	stimpar_l = rmsml;
      sprintf(ggstr,"   STIMPAR: rmsm_L = %f\n",stimpar_l);
      mylog(mylogf,ggstr);
      if (stimpar_l < 0.0)
	mylog_exit(mylogf,
	  "MODEL_ME_V5_GET_RESPONSE  Either rmsm_L or stimpar_msc < 0");

      a1  = onode_getpar_flt_exit(mv5->normo,"alpha_1");
      a2  = onode_getpar_flt_exit(mv5->normo,"alpha_2") / (float)ndir;
      //c3 = onode_getpar_flt_exit(mv5->normo,"alpha_3") * stimpar_l;
      a4  = onode_getpar_flt_dflt(mv5->normo,"alpha_4",0.0);
      nos = onode_getpar_flt_dflt(mv5->normo,"output_scale",1.0);

      if (mv5->binoc == 1){  // Compute a separate right eye signal
	stimpar_lr = param_getf_dflt(s->ppl,"stimmod_msc_r",stimpar_l);
	if (stimpar_lr < 0.0)
	  stimpar_lr = stimpar_l;
	c3r = onode_getpar_flt_dflt(mv5->normo,"alpha_3",0.0) * stimpar_lr;
      }

    }else if (onode_test_chr(mv5->normo,"type","local_stf")==1){

      //  Compute the number of channels to divide the sum
      nc = ndir;
      if (mv5->v1stf != NULL)
	nc *= mv5->v1stf->n;

      a1  = onode_getpar_flt_exit(mv5->normo,"alpha_1");
      a2  = onode_getpar_flt_exit(mv5->normo,"alpha_2") / (float)nc;
      c3  = c3r = 0.0;
      a4  = onode_getpar_flt_dflt(mv5->normo,"alpha_3",0.0);
      nos = onode_getpar_flt_dflt(mv5->normo,"output_scale",1.0);
    }else{
      mylog_exit(mylogf,"MODEL_ME_V5_GET_RESPONSE  Unknown <norm> type");
    }
  }

  //
  //  V1 normalize.   *** This handles several steps for DISPARITY case ***
  //
  //  V1 surround signal (if present) is applied here.
  //
  if (rawflag == 0)  // Default flow
    mod_me_v5_norm_v1(a1,a2,a4,c3,c3r,nos);


  //
  //  (C) Optional binoc integration [AS SHOWN IN OUR FIGURE]
  //
  mod_me_v5_mix_v1(mylogf,"normalized");


  if (mv5->v1_dflag == 0){
    //
    //  (D) Opponent sum
    //
    mod_me_v5_opp_v1();
  }


  mod_me_v5_mix_v1(mylogf,"opponent");  // Optional mixing


  //
  //  Binocular sum
  //
  whack_flag = 0;
  monoc_gate = 0;
  if (mod_me_mpt == NULL){
    //
    //  Tailby weight is only considered for old way, and even this perhaps
    //  should be eliminated.
    //
    whack = mv5->bin_hack;  // weight factor for Tailby hack
    if (whack >= 0.0){
      whack_flag = 1;
      sprintf(tstr,"  *** Tailby Hack: scale negative weights by %f\n",whack);
      mylog(mylogf,tstr);
    }
  }


  //
  //  MT
  //
  if (mod_me_mpt == NULL){
    //
    //  No spatial summation.
    //    In this case, note that the parameter 'right_eye_w' from <v1> is
    //    being used here for the eye bias in the MT summation.  This param
    //    has a default value of 1.0 if it is not present.
    //

    if (mv5->w != NULL){ // OLD WAY
      mta = mv5->mt_a;  // Constants for MT nonlinearity
      mtb = mv5->mt_b;
    }else{
      // *** WYETH - Or is this always to be done the old way ???
      exit_error("WYETH","New way not implemented yet");
    }

    //  Compute MT signals
    mod_me_v5_mt(whack_flag,whack,mv5->w,mv5->w_r,mv5->v1_bin_w);

    if (mod_me_ftop == NULL){ // old way
      mv5->mt_norm = mod_me_v5_mt_nl(mv5->v1resp[0][0]->mtb,tn,
				     mv5->mt_nltype,mta,mtb);
    }else{
      stfi = 0; // *** WYETH FIX
      mv5->mt_norm = mod_me_v5_mt_nl(mv5->v1stf->sgr[stfi]->bdgr[0][0]->mtb,tn,
				     mv5->mt_nltype,mta,mtb);
    }
  }else{
    //
    //  Spatial summation
    //
    for(i=0;i<mv5->mt_npop;i++){  // Compute responses for each 

      pl = mod_me_mpt->lay[mv5->mtli0 + i];  // Layer from 'pop_top'
      mtpop = mv5->mtp[i];           // MT response struct

      mod_me_v5_pop(pl,mtpop);  // NEW WAY Jun 2016, no stack
    }
  }


  //
  //  Create seeds for spike generation, used in handover (below)
  //
  if (mod_me_mpt != NULL){

    //  Create a seed for each layer, from the trial seed
    lseedlist = get_seeds(m->mseed[r->tsi],1000000,mv5->mt_npop);

    for(i=0;i<mv5->mt_npop;i++){  // For each layer

      pl = mod_me_mpt->lay[mv5->mtli0 + i];  // Layer from 'pop_top'
      mtpop = mv5->mtp[i];                   // MT pop struct

      if (mtpop->seed != NULL)
	free_3d_iarray(mtpop->seed,pl->xn,pl->yn,pl->zn);

      //  **** WYETH HERE - MUST FIX SEED PROBLEM ******
      mtpop->seed = get_seeds_3d(lseedlist[i],100000,pl->xn,pl->yn,pl->zn);
    }
    myfree(lseedlist);
  }

  //
  //  Response Handover
  //
  mod_me_v5_handover(m,r);


  mod_me_stimno_prev = s->stimno;

  //
  //  Exit if there are responses that have not been filled
  //
  mod_util_resp_unfilled_exit(r,mylogf);
}
/**************************************-**************************************/
/*                                                                           */
/*                       MODEL_FILTER_XYT_GET_RESPONSE                       */
/*                                                                           */
/*****************************************************************************/
void model_filter_xyt_get_response(m,s,r,k)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response data
     int k;                      // repeat number, k=0 means new stim
{
  int i;
  int xn,yn,tn,x,y;
  float ***fpe,**fpe_s,***rpe,*mep,*resp_e;
  char *spikegen;
  struct onode *sgo;

  xn    = mod_me_xn;
  yn    = mod_me_yn;
  tn    = mod_me_tn;
  fpe   = mod_me_fpe;
  fpe_s = mod_me_fpe_s;

  if ((k==0) || (mod_me_rpe == NULL) || (s->stimno != mod_me_stimno_prev)){
    // Must recompute
    if (mod_me_rpe != NULL){
      free_f3tensor(mod_me_rpe,1,xn,1,yn,1,tn);
    }
    // Next call modifies s->d

    mod_me_compute_response_simp_from_ft(s->d,fpe,fpe_s,xn,yn,tn,&rpe);

    mod_me_rpe = rpe;  // Store for future 'repeats'
  }else{
    rpe = mod_me_rpe;  // Use stored response
  }

  x = xn/2; // Use the unit at the center of the model
  y = yn/2;
  resp_e = &(rpe[x][y][1]);


  // Apply a non-linearity, perhaps a square root, if desired
  if (mod_me_pow != 1.0){
    mep = (float *)myalloc(tn*sizeof(float));
    for(i=0;i<tn;i++){
      mep[i] = pow(resp_e[i],mod_me_pow);
    }
  }else
    mep = copy_farray(resp_e,tn);

  // Now, 'mep' is 'tn' = 'tnr' long, and contains the response trace

  mod_me_stimno_prev = s->stimno;

  //
  //  Response Handover
  //
  for(i=0;i<r->n;i++){ // For each response requested
    if (strcmp(r->datid[i],"spikes")==0){
      //
      //  Generate spikes
      //
      if (m->ppl == NULL){
	sgo = onode_child_get_unique(m->o,"spike_gen");
	spikegen = onode_getpar_chr_exit(sgo,"type");
      }else
	spikegen = paramfile_get_char_param_or_exit(m->ppl,"spike_gen");
      
      if (strcmp(spikegen,"poisson")==0){

	//
	//  *** WYETH - PROBLEM HERE - THIS ROUTINE LOOKS TO SAVE OTHER
	//  ***   responses
	//
	mod_me_gen_hand_poisson(m,r,mep,NULL,tn,0);
      }else{
	printf("spikegen = %s\n",spikegen);
	mylog_exit(mylogf,
		   "MODEL_FILTER_XYT_GET_RESPONSE  unknown spike_gen");
      }
      myfree(spikegen);

      //
      //  WYETH - CHOICES BELOW ARE FOR gabor_comp
      //  *** SHOULD specify CHOICES FOR rd_2gabor, which can use 
      //      'even', 'odd', and 'response', but not the squared.
      //
    }else if (strcmp(r->datid[i],"response")==0){
      mod_util_resp_store_f_samp(r,i,resp_e,tn,1.0/mod_me_tscale,mylogf);
    }else if (strcmp(r->datid[i],"response_power")==0){
      mod_util_resp_store_f_samp(r,i,mep,tn,1.0/mod_me_tscale,mylogf);
    }else{
      mylog_exit(mylogf,"MOD_FILTER_XYT_GET_RESPONSE unknown response item");
    }
  }

  myfree(mep);
}
/**************************************-**************************************/
/*                                                                           */
/*                       MODEL_GABOR_COMP_GET_RESPONSE                       */
/*                                                                           */
/*****************************************************************************/
void model_gabor_comp_get_response(m,s,r,k)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response data
     int k;                      // repeat number, k=0 means new stim
{
  int i;
  int xn,yn,tn,x,y,ri,dti;
  float ***fpe,***fpo,**fpe_s,**fpo_s,***rpe,***rpo,*mep;
  float *resp_e,*resp_o,*resp_e2,*resp_o2,sampling;
  char *spikegen;

  xn = mod_me_xn;
  yn = mod_me_yn;
  tn = mod_me_tn;
  fpe = mod_me_fpe;
  fpo = mod_me_fpo;

  fpe_s = mod_me_fpe_s;  fpo_s = mod_me_fpo_s;

  if ((k==0) || (mod_me_rpe == NULL) || (s->stimno != mod_me_stimno_prev)){
    // Must recompute
    if (mod_me_rpe != NULL){
      free_f3tensor(mod_me_rpe,1,xn,1,yn,1,tn);
      free_f3tensor(mod_me_rpo,1,xn,1,yn,1,tn);
    }
    // Next call modifies s->d
    mod_me_compute_response_gabor_comp_from_ft(s->d,fpe,fpo,fpe_s,fpo_s,
					       xn,yn,tn,&rpe,&rpo);

    mod_me_rpe = rpe; mod_me_rpo = rpo;  // Store for future 'repeats'
  }else{
    rpe = mod_me_rpe; rpo = mod_me_rpo;  // Use stored response
  }

  x = xn/2; // Use the unit at the center of the model
  y = yn/2;
  mep = (float *)myalloc(tn*sizeof(float));

  resp_e = &(rpe[x][y][1]);
  resp_o = &(rpo[x][y][1]);
  resp_e2 = (float *)myalloc(tn*sizeof(float));
  resp_o2 = (float *)myalloc(tn*sizeof(float));

  if (strcmp(mod_me_modtype,"rd_2gabor")==0){
    //
    //  rd_2gabor
    //

    if (mod_rd_2gabor_rect == 1){
      half_wave_rectify_farray(resp_e,tn);
      half_wave_rectify_farray(resp_o,tn);
    }

    dti = mod_rd_2gabor_dti;
    for(i=0;i<tn;i++){
      if (i >= dti)
	mep[i] = resp_e[i-dti] * resp_o[i];  // RD multiply
      else
	mep[i] = resp_o[i];  // No data to multiply against
    }
  }else{
    //
    //  "gabor_comp"
    //
    for(i=0;i<tn;i++){
      resp_e2[i] = resp_e[i]*resp_e[i];  // Squaring
      resp_o2[i] = resp_o[i]*resp_o[i];
      mep[i] = resp_e2[i] + resp_o2[i];
    }
  }

  // Apply a non-linearity, perhaps a square root, if desired
  if (mod_me_pow != 1.0){
    for(i=0;i<tn;i++){
      mep[i] = pow(mep[i],mod_me_pow);
    }
  }

  // Now, 'mep' is 'tn' = 'tnr' long, and contains the response trace

  mod_me_stimno_prev = s->stimno;


  //
  //  Response Handover
  //
  for(i=0;i<r->n;i++){ // For each response requested
    if (strcmp(r->datid[i],"spikes")==0){
      {// Generate spikes
	char *spikegen;
	struct onode *sgo;

	if (m->ppl == NULL){
	  sgo = onode_child_get_unique(m->o,"spike_gen");
	  spikegen = onode_getpar_chr_exit(sgo,"type");
	}else
	  spikegen = paramfile_get_char_param_or_exit(m->ppl,"spike_gen");

	if (strcmp(spikegen,"poisson")==0){
	  mod_me_gen_hand_poisson(m,r,mep,NULL,tn,0);
	}else{
	  printf("spikegen = %s\n",spikegen);
	  mylog_exit(mylogf,
		     "MODEL_GABOR_COMP_GET_RESPONSE  unknown spike_gen");
	}
	myfree(spikegen);
      }
      //
      //  WYETH - CHOICES BELOW ARE FOR gabor_comp
      //  *** SHOULD specify CHOICES FOR rd_2gabor, which can use 
      //      'even', 'odd', and 'response', but not the squared.
      //
    }else if (strcmp(r->datid[i],"response")==0){
      mod_util_resp_store_f_samp(r,i,mep,tn,1.0/mod_me_tscale,mylogf);
    }else if (strcmp(r->datid[i],"even")==0){
      mod_util_resp_store_f_samp(r,i,resp_e,tn,1.0/mod_me_tscale,mylogf);
    }else if (strcmp(r->datid[i],"even_squared")==0){
      mod_util_resp_store_f_samp(r,i,resp_e2,tn,1.0/mod_me_tscale,mylogf);
    }else if (strcmp(r->datid[i],"odd")==0){
      mod_util_resp_store_f_samp(r,i,resp_o,tn,1.0/mod_me_tscale,mylogf);
    }else if (strcmp(r->datid[i],"odd_squared")==0){
      mod_util_resp_store_f_samp(r,i,resp_o2,tn,1.0/mod_me_tscale,mylogf);
    }else{
      mylog_exit(mylogf,"MOD_GABOR_COMP_GET_RESPONSE unknown response item\n");
    }
  }

  myfree(mep);
  myfree(resp_e2);
  myfree(resp_o2);
}
/**************************************-**************************************/
/*                                                                           */
/*                        MODEL_DS_SIMP_01_GET_RESPONSE                      */
/*                                                                           */
/*  Calling sequence:                                                        */
/*    1.  model_me_01_prep                                                   */
/*    2.  model_me_01_prep_ft                                                */
/*    3.  this routine                                                       */
/*        ...                                                                */
/*    N.  model_me_01_done                                                   */
/*                                                                           */
/*****************************************************************************/
void model_ds_simp_01_get_response(m,s,r,k)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response data
     int k;                      // repeat number, k=0 means new stim
{
  int i;
  int xn,yn,tn,x,y,qflag;
  float ***fpe,**fpe_s,***rpe,*mep,tr,***tfil;
  float ***fpo,**fpo_s,***rpo;
  char *spikegen;
  struct onode *sgen;
  struct mod_me_fstruct **fs;

  xn    = mod_me_xn;
  yn    = mod_me_yn;
  tn    = mod_me_tn;

  fs = mod_me_ff;

  if (strcmp(mod_me_filt,"two")==0)
    qflag = 1;
  else
    qflag = 0;


  if ((strcmp(mod_me_filt,"even")==0)|| qflag){
    fpe   = mod_me_fpe;
    fpe_s = fs[0]->s; // OLD WAY:  mod_me_fpe_s;
  }else if (strcmp(mod_me_filt,"odd")==0){
    fpe   = mod_me_fpo;
    fpe_s = fs[1]->s; // OLD WAY:  mod_me_fpo_s;
  }else{
    exit_error("MODEL_DS_SIMP_01_GET_RESPONSE","Unknown 'me_filt'");
  }
  if (qflag){
    fpo   = mod_me_fpo;
    fpo_s = fs[1]->s; // OLD WAY mod_me_fpo_s;
  }

  if ((k==0) || (mod_me_rpe == NULL) || (s->stimno != mod_me_stimno_prev)){
    // Compute new response

    if (mod_me_rpe != NULL){
      free_f3tensor(mod_me_rpe,1,xn,1,yn,1,tn);
      if (qflag)
	free_f3tensor(mod_me_rpo,1,xn,1,yn,1,tn);
    }

    // Next call modifies s->d

    mod_me_compute_response_simp_from_ft(s->d,fpe,fpe_s,xn,yn,tn,&rpe);

    if (qflag)
      mod_me_compute_response_simp_from_ft(s->d,fpo,fpo_s,xn,yn,tn,&rpo);
    
    mod_me_rpe = rpe; // Store for future 'repeats'
    if (qflag)
      mod_me_rpo = rpo; // Store for future 'repeats'
  }else{
    rpe = mod_me_rpe; // Use stored response
    if (qflag)
      rpo = mod_me_rpo;
  }


  x = xn/2;
  y = yn/2;
  mep = (float *)myalloc(tn*sizeof(float));
  if (qflag == 1){
    for(i=0;i<tn;i++){
      tr  = rpe[x][y][1+i];           // Even
      if (tr <= 0.0)                  // Rectify
	tr = 0.0;
      if (mod_me_pow != 1.0){
	mep[i] = pow(tr,mod_me_pow);  // Raise to power
      }else
	mep[i] = tr;

      tr  = rpo[x][y][1+i];           // Odd
      if (tr <= 0.0)
	tr = 0.0;
      if (mod_me_pow != 1.0){
	mep[i] += pow(tr,mod_me_pow);
      }else
	mep[i] += tr;
    }
  }else{
    //
    //  Regular 1-branch model
    //
    for(i=0;i<tn;i++){
      tr = rpe[x][y][1+i];
      if (mod_me_pow != 1.0){
	mep[i] = pow(tr,mod_me_pow);
      }else
	mep[i] = tr;
    }
  }

  mod_me_stimno_prev = s->stimno;

  //
  //  Handover any basic float responses
  //
  //  Uses "filter_out" for raw filter data
  mod_me_handover_basic(m,r,NULL,NULL,mep,NULL,rpe,NULL,NULL,NULL,xn,yn,tn);
  if (qflag)
    mod_me_handover_basic(m,r,NULL,NULL,mep,NULL,rpe,rpo,NULL,NULL,xn,yn,tn);


  //
  //  Generate/handover spike responses
  //
  if (m->ppl == NULL){
    sgen = onode_child_get_unique(m->o,"spike_gen");
    spikegen = onode_getpar_chr_exit(sgen,"type");
  }else
    spikegen = paramfile_get_char_param_or_exit(m->ppl,"spike_gen");
  if (strcmp(spikegen,"poisson")==0){
    mod_me_gen_hand_poisson(m,r,mep,NULL,tn,0);
  }else  if (strcmp(spikegen,"ifcn")==0){
    mod_me_gen_hand_ifc(m,r,mep,NULL,tn);  // gti will be zero
  }else{
    printf("spike_gen = %s\n",spikegen);
    mylog_exit(mylogf,"MODEL_DS_SIMP_01_GET_RESPONSE  Unknown spike_gen.");
  }
  myfree(spikegen);
  myfree(mep);


  //
  //  Exit if there are responses that have not been filled
  //
  mod_util_resp_unfilled_exit(r,mylogf);
}
/**************************************-**************************************/
/*                                                                           */
/*                         MODEL_DS_LIV3_GET_RESPONSE                        */
/*                                                                           */
/*  Built from '...ds_simp_01..." to emulate the 3-subunit model of          */
/*  Livingstone and Conway (2003).                                           */
/*                                                                           */
/*****************************************************************************/
void model_ds_liv3_get_response(m,s,r,k)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response data
     int k;                      // repeat number, k=0 means new stim
{
  int i;
  int xn,yn,tn,x,y,dy,locflag;
  float ***fpe,**fpe_s,***rpe,*mep,tr,***tfil,yoff,ssd;
  float *s1,*s2,*s3;  // Three signals from filter
  char *spikegen;

  xn    = mod_me_xn;
  yn    = mod_me_yn;
  tn    = mod_me_tn;

  if (strcmp(mod_me_filt,"even")==0){
    fpe   = mod_me_fpe;
    fpe_s = mod_me_fpe_s;
  }else if (strcmp(mod_me_filt,"odd")==0){
    fpe   = mod_me_fpo;
    fpe_s = mod_me_fpo_s;
  }else{
    exit_error("MODEL_DS_SIMP_01_GET_RESPONSE","Unknown 'me_filt'");
  }


  if ((k==0) || (mod_me_rpe == NULL) || (s->stimno != mod_me_stimno_prev)){
    // Compute new response

    if (mod_me_rpe != NULL){
      free_f3tensor(mod_me_rpe,1,xn,1,yn,1,tn);
    }

    // Next call modifies s->d
    mod_me_compute_response_simp_from_ft(s->d,fpe,fpe_s,xn,yn,tn,&rpe);

    mod_me_rpe = rpe; // Store for future 'repeats'
  }else{
    rpe = mod_me_rpe; // Use stored response
  }

  //
  // WYETH - Should be done in a 'prep' routine
  //
  yoff = paramfile_get_float_param_or_exit(m->ppl,"liv3_offset_sd");
  ssd  = paramfile_get_float_param_or_exit(m->ppl,"me_ssd");
  locflag = paramfile_get_float_param_or_exit(m->ppl,"liv3_loc_flag");

  x = xn/2;
  y = yn/2;
  dy = my_rint(yoff * ssd / mod_me_sscale);

  if ((y+dy) >= yn)
    mylogx(mylogf,"MODEL_DS_LIV3_GET_RESPONSE","liv3_offset_sd too large");
  
  //printf(" YOFFSET = %f SD,  %d pixels\n",yoff,dy);

  s1 = (float *)myalloc(tn*sizeof(float));
  s2 = (float *)myalloc(tn*sizeof(float));
  s3 = (float *)myalloc(tn*sizeof(float));
  mep = (float *)myalloc(tn*sizeof(float));
  for(i=0;i<tn;i++){
    s1[i] = rpe[x][y-dy][1+i];
    s2[i] = rpe[x][y   ][1+i];
    s3[i] = rpe[x][y+dy][1+i];
  }

  if (locflag == 0){  // Sum first
    for(i=0;i<tn;i++){
      mep[i] = s1[i] + s2[i] + s3[i];
    }
    half_wave_rectify_farray(mep,tn);
  }else if (locflag == 1){
    half_wave_rectify_farray(s1,tn); // Rectify first
    half_wave_rectify_farray(s2,tn);
    half_wave_rectify_farray(s3,tn);
    for(i=0;i<tn;i++){
      mep[i] = s1[i] + s2[i] + s3[i];
    }
  }else if (locflag == 2){            // Use only middle subunit
    half_wave_rectify_farray(s2,tn);
    for(i=0;i<tn;i++){
      mep[i] = s2[i];
    }
  }else{
    mylogx(mylogf,"MODEL_DS_LIV3_GET_RESPONSE","Bad 'liv3_loc_flag' value");
  }
  square_farray(mep,tn);

  myfree(s1); myfree(s2); myfree(s3);
  
  mod_me_stimno_prev = s->stimno;

  //
  //  Handover any basic float responses
  //
  //  Uses "filter_out" for raw filter data
  mod_me_handover_basic(m,r,NULL,NULL,mep,NULL,rpe,NULL,NULL,NULL,xn,yn,tn);


  //
  //  Generate/handover spike responses
  //
  spikegen = paramfile_get_char_param_or_exit(m->ppl,"spike_gen");
  if (strcmp(spikegen,"poisson")==0){
    mod_me_gen_hand_poisson(m,r,mep,NULL,tn,0);
  }else  if (strcmp(spikegen,"ifcn")==0){
    mod_me_gen_hand_ifc(m,r,mep,NULL,tn);  // gti will be zero
  }else{
    printf("spike_gen = %s\n",spikegen);
    mylog_exit(mylogf,"MODEL_DS_SIMP_01_GET_RESPONSE  Unknown spike_gen.");
  }
  myfree(spikegen);
  myfree(mep);

  //
  //  Exit if there are responses that have not been filled
  //
  mod_util_resp_unfilled_exit(r,mylogf);
}
/**************************************-**************************************/
/*                                                                           */
/*                         MODEL_MEM_REGULATE_RESPONSE                       */
/*                                                                           */
/*****************************************************************************/
float *model_mem_regulate_response(m,topp)
     struct model_struct *m;     // Model params
     float **topp;               // [nf][tn]
{
  int i,j;
  int nf,tn,imax;
  float **wf,*maxa,f,f1,*r;
  char tstr[SLEN];
  struct mod_mem *mem;

  mem = mod_me_mem; // Global structure
  tn = mod_me_tn;

  nf = mem->nf;

  // 'f' is the factor that controls how quickly weights change (was 0.9)
  f = mem->wf;
  f1 = 1.0 - f;

  //
  // Compute weights 'wf' for each filter over all time
  //
  //    The current algorithm does not maintain a constant sum of
  //  filter weights, but it tries to push one to 1 and the rest to 0.
  //
  wf = get_2d_farray(nf,tn);
  for(j=0;j<nf;j++)
    wf[j][0] = 1.0/(float)nf;  // Start all weights at same value
  
  maxa = get_farray(nf);
  for(i=1;i<tn;i++){
    // Find channel 'imax' has max absolute value
    for(j=0;j<nf;j++)
      maxa[j] = fabs(topp[j][i]);
    imax = max_index_farray(maxa,nf);

    for(j=0;j<nf;j++){ // For each ME filter
      if (j == imax)
	wf[j][i] = f * wf[j][i-1] + f1 * 1.0;  // Decay towards 1
      else
	wf[j][i] = f * wf[j][i-1] + f1 * 0.0;  // Decay towards 0
    }
  }
  myfree(maxa);

  //
  // Dump weights
  //
  if ((myid == -1) && (mem->dumpflag == 1)){
    for(j=0;j<nf;j++){
      sprintf(tstr,"W_%d",j);
      append_farray_plot("zzz.resp.pl",tstr,wf[j],tn,1);
    }
  }

  //
  // Create final weighted ressponse 'r'
  //
  r = get_zero_farray(tn);
  for(i=0;i<tn;i++){
    for(j=0;j<nf;j++){
      r[i] += wf[j][i] * topp[j][i];
    }
  }

  free_2d_farray(wf,nf);

  return r;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MODEL_MEM_GET_RESPONSE                          */
/*                                                                           */
/*****************************************************************************/
void model_mem_get_response(m,s,r,k)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response data
     int k;                      // repeat number, k=0 means new stim
{
  int i,j;
  int xn,yn,tn,x,y,ns;
  float *mp,*ma,*mr0,*mr1,*mr2,*mr3,*wr;
  float *totp,*tota,**topp,*fs;
  char tstr[SLEN],*spikegen;
  struct mod_mem *mem;
  struct onode *sgo;

  mem = mod_me_mem; // Global structure
  xn = mod_me_xn;
  yn = mod_me_yn;
  tn = mod_me_tn;

  //  If 0, use opponent signal into excitatory pathway,
  //  if 1, use Pref signal excitatory, and Null signal inhibitory
  //  oppflag = paramfile_get_int_param_default(m->ppl,"opp_flag",0);
  
  if ((k==0) || (mem->r[0][0] == NULL) || (s->stimno != mod_me_stimno_prev)){
    if (mem->r[0][0] != NULL){  // Free previous responses
      for(i=0;i<mem->nf;i++){
	free_f3tensor(mem->r[i][0],1,xn,1,yn,1,tn);
	free_f3tensor(mem->r[i][1],1,xn,1,yn,1,tn);
	free_f3tensor(mem->r[i][2],1,xn,1,yn,1,tn);
	free_f3tensor(mem->r[i][3],1,xn,1,yn,1,tn);
      }
    }
    // Next call modifies s->d, computes new responses:  mem->r[nf][4]
    mod_me_compute_response_mem(s->d,mem,xn,yn,tn);
  }

  // Compute Pref and Anti response traces for each coordinate
  for(i=0;i<mem->nf;i++){  // For each filter
    for(j=0;j<mem->nr;j++){  // For each response coordinate
      x = mem->xi[j];
      y = mem->yi[j];
      mr0 = mem->r[i][0][x][y];
      mr1 = mem->r[i][1][x][y];
      mr2 = mem->r[i][2][x][y];
      mr3 = mem->r[i][3][x][y];

      mp = mem->mp[i][j];
      ma = mem->ma[i][j];

      for(k=0;k<tn;k++){
	mp[k] = mr0[1+k]*mr0[1+k] + mr1[1+k]*mr1[1+k];
	ma[k] = mr2[1+k]*mr2[1+k] + mr3[1+k]*mr3[1+k];
      }

      if ((myid == -1) && (mem->dumpflag == 1)){
	if (mem->nr == 1){
	  sprintf(tstr,"rf_%d__%d_%d_pref",i,x,y);
	  append_farray_plot("zzz.resp.pl",tstr,mp,tn,1);
	  sprintf(tstr,"rf_%d__%d_%d_anti",i,x,y);
	  append_farray_plot("zzz.resp.pl",tstr,ma,tn,1);
	}
      }

    }
  }

  // Get total oppenent signals in each filter channel
  topp = (float **)myalloc(mem->nf*sizeof(float *));
  
  if (mem->nr > 1){
    
    totp = get_zero_farray(tn);
    tota = get_zero_farray(tn);
    
    for(i=0;i<mem->nf;i++){  // For each filter channel
      for(j=0;j<mem->nr;j++){  // For each response coordinate
	for(k=0;k<tn;k++){
	  totp[k] += mem->mp[i][j][k];
	  tota[k] += mem->ma[i][j][k];
	}
      }
      multiply_farray(totp,tn,1.0/(float)mem->nr);
      multiply_farray(tota,tn,1.0/(float)mem->nr);

      topp[i] = diff_farrays(totp,tota,tn);

      // Apply scaling factor from .mod file to each channel
      if (mem->scale[i] != 1.0)
	multiply_farray(topp[i],tn,mem->scale[i]);
	
      /*
      sprintf(tstr,"Avg_%d_pref",i);
      append_farray_plot("zzz.resp.pl",tstr,totp,tn,1);
      sprintf(tstr,"Avg_%d_anti",i);
      append_farray_plot("zzz.resp.pl",tstr,tota,tn,1);*/

      if ((myid == -1) && (mem->dumpflag == 1)){
	sprintf(tstr,"T_opp_%d",i);
	append_farray_plot("zzz.resp.pl",tstr,topp[i],tn,1);
      }

      zero_farray(totp,tn);
      zero_farray(tota,tn);
    }
  }

  // Get the weighted response
  wr = model_mem_regulate_response(m,topp);

  // WYETH - testing
  if ((myid == -1) && (mem->dumpflag == 1)){
    append_farray_plot("zzz.resp.pl","WR",wr,tn,1);
  }

  mod_me_stimno_prev = s->stimno;

  // GENERATE SPIKES
  sgo = onode_child_get_unique(m->o,"spike_gen");
  spikegen = onode_getpar_chr_exit(sgo,"type");
  if (strcmp(spikegen,"poisson")==0){
    mod_me_gen_hand_poisson(m,r,wr,NULL,tn,0);
  }else{
    printf("spikegen = %s\n",spikegen);
    mylog_exit(mylogf,"MODEL_MEM_GET_RESPONSE  unknown spike_gen");
  }
  myfree(spikegen);

}
/**************************************-**************************************/
/*                                                                           */
/*                          MODEL_REICH_GET_RESPONSE                         */
/*                                                                           */
/*****************************************************************************/
void model_reich_get_response(m,s,r,k)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response data
     int k;                      // repeat number, k=0 means new stim
{
  int i,j;
  int xn,yn,tn,ri,x0,x1,y0,y1,ti,dx,n_unit;
  float ***rlp,***rhp,*tlp_a,*tlp_b,*thp_a,*thp_b,*tt;
  char pname[SLEN];

  xn = mod_me_xn;
  yn = mod_me_yn;
  tn = mod_me_tn;

  if ((k==0) || (mod_reich_rlp == NULL) ||(s->stimno != mod_me_stimno_prev)){
    // Must recompute
    if (mod_reich_rlp != NULL){
      free_f3tensor(mod_reich_rlp,1,xn,1,yn,1,tn);
      free_f3tensor(mod_reich_rhp,1,xn,1,yn,1,tn);
    }
    // Next call modifies s->d
    mod_reich_compute_response_from_ft(s->d,mod_reich_lp,mod_reich_hp,
				       mod_reich_lp_s,mod_reich_hp_s,
				       xn,yn,tn,&rhp,&rlp);
    mod_reich_rlp = rlp; // Store for future 'repeats'
    mod_reich_rhp = rhp;

    mod_me_stimno_prev = s->stimno;
  }else{
    rlp = mod_reich_rlp;  // Use stored response
    rhp = mod_reich_rhp;
  }
  
  //  RD opponent subunit:
  //
  //     left  right      left  right    |
  //     A  B  A  B       B  A  B  A     |
  //     |  |  |  |       |  |  |  |     |
  //    LP HP  HP LP     LP HP  HP LP    |
  //     |/      \|       |/      \|     |
  //     X        X       X        X     |
  //      \      /         \      /      |
  //        subtr            subtr       |
  //
  tt = get_zero_farray(tn);    // Accumulate response over RD units

  x0 = mod_reich_ux0;          // Limits of integration over RD units
  x1 = x0 + mod_reich_uxn - 1;
  y0 = mod_reich_uy0;
  y1 = y0 + mod_reich_uyn - 1;
  dx = mod_reich_dx;

  for(i=x0;i<=x1;i++){
    for(j=y0;j<=y1;j++){
      tlp_a = rlp[i   ][j];  // Position 'A'
      thp_a = rhp[i   ][j];  // 
      tlp_b = rlp[i+dx][j];  // Position 'B'
      thp_b = rhp[i+dx][j];  // 

      if (mod_me_oppflag == 1){
	for(ti=0;ti<tn;ti++){
	  tt[ti] += tlp_a[ti]*thp_b[ti] - thp_a[ti]*tlp_b[ti];
	}
      }else{
	for(ti=0;ti<tn;ti++){
	  tt[ti] += tlp_a[ti]*thp_b[ti];
	}
      }
    }
  }

  // Change total into an average
  n_unit = mod_reich_uxn * mod_reich_uyn;
  multiply_farray(tt,tn,1.0/(float)n_unit);


  if (mod_reich_dump == 1){ // Dump LP and HP outputs at central position
    int x,y;
    float *tla,*tra,*tdiff,*tdump;
    
    x = xn/2;
    y = yn/2;

    tla = (float *)myalloc(tn*sizeof(float));
    tra = (float *)myalloc(tn*sizeof(float));
    tdiff = (float *)myalloc(tn*sizeof(float));
    for(ti=0;ti<tn;ti++){
      tla[ti] = rlp[x][y][ti] * rhp[x+dx][y][ti];
      tra[ti] = rhp[x][y][ti] * rlp[x+dx][y][ti];
      tdiff[ti] = tla[ti] - tra[ti];
    }

    sprintf(pname,"%s_HP_ctr",r->tname);
    append_farray_plot("zzz.rd.pl",pname,rhp[x][y],tn,1);
    sprintf(pname,"%s_LP_ctr+%d",r->tname,dx);
    append_farray_plot("zzz.rd.pl",pname,rlp[x+dx][y],tn,1);
    sprintf(pname,"%s_right_arm",r->tname);
    append_farray_plot("zzz.rd.pl",pname,tra,tn,1);
    sprintf(pname,"%s_HP_ctr+%d",r->tname,dx);
    append_farray_plot("zzz.rd.pl",pname,rhp[x+dx][y],tn,1);
    sprintf(pname,"%s_LP_ctr",r->tname);
    append_farray_plot("zzz.rd.pl",pname,rlp[x][y],tn,1);
    sprintf(pname,"%s_left_arm",r->tname);
    append_farray_plot("zzz.rd.pl",pname,tla,tn,1);
    sprintf(pname,"%s_la-ra",r->tname);
    append_farray_plot("zzz.rd.pl",pname,tdiff,tn,1);
    sprintf(pname,"%s_Tot_Response",r->tname);
    append_farray_plot("zzz.rd.pl",pname,tt,tn,1);

    myfree(tla);
    myfree(tra);
    myfree(tdiff);
  }

  //
  //  Response Handover
  //
  for(i=0;i<r->n;i++){ // For each response requested
    if (strcmp(r->datid[i],"spikes")==0){
      {// Generate spikes
	char *spikegen;
	struct onode *sgo;
	
	sgo = onode_child_get_unique(m->o,"spike_gen");
	
	spikegen = onode_getpar_chr_exit(sgo,"type");
	if (strcmp(spikegen,"poisson")==0){
	  mod_me_gen_hand_poisson(m,r,tt,NULL,tn,0);
	}else{
	  printf("spikegen = %s\n",spikegen);
	  mylog_exit(mylogf,"MODEL_REICH_GET_RESPONSE  unknown spike_gen");
	}
	myfree(spikegen);
      }
    }else if (strcmp(r->datid[i],"response")==0){
      mod_util_resp_store_f_samp(r,i,tt,tn,1.0/mod_me_tscale,mylogf);
    }else{
      mylog_exit(mylogf,"MOD_REICH_GET_RESPONSE  unknown response item");
    }
  }

  myfree(tt);
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_ME_RUN_ME_GABOR_01                          */
/*                                                                           */
/*  Motion energy model with Gabor (and perhaps AB temporal) filters.        */
/*                                                                           */
/*****************************************************************************/
void mod_me_run_me_gabor_01(m,s,r,k,action)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     int k;                     // Trial number, when action = 1
     int action;                // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  if (action == 0){
    if (m->ppl == NULL)
      model_me_01_prep_o(m);
    else
      model_me_01_prep(m);
    model_me_01_prep_ft();
  }else if (action == 1){

    if (me_xflag == 0)
      model_me_01_get_response(m,s,r,k);
    else
      model_me_01_get_response_x(m,s,r,k);

  }else if (action == -1){
    model_me_01_done(m->ppl);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_ME_V5_RSPINFO                            */
/*                                                                           */
/*****************************************************************************/
void mod_me_v5_rspinfo()
{
  //         xn  yn  dir
  //  v1     32  32   12   raw, surr, norm, sum
  //  v1r    32  32   12   raw, surr, norm, sum
  //  
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_ME_RUN_ME_V5                              */
/*                                                                           */
/*  V5 model build on motion energy filters.                                 */
/*                                                                           */
/*****************************************************************************/
void mod_me_run_me_v5(m,s,r,k,action)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     int k;                     // Repeat number, when action = 1
     int action;                // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  //
  // WYETH - 2018 rspinf
  //
  //   The model should be able to tell information about responses it
  //   can return.
  //
  //   20 - rspinf
  //    (1) call prep
  //    (2) call resp_info
  //
  //
  //

  if ((action ==  0) ||
      (action == 10) ||  // 10 - write .mar file
      (action == 20)){   // 20 - write MDS list
    model_me_v5_prep(m,s);
    if (action == 0)
      model_me_01_prep_ft();
  }else if (action == 1){
    model_me_v5_get_response(m,s,r,k);
  }else if (action == -1){
    model_me_v5_done(m,s,r);  // May write adaptation state
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_ME_RUN_FILTER_XYT                          */
/*                                                                           */
/*  Use for 'GaborSimp'                                                      */
/*                                                                           */
/*****************************************************************************/
void mod_me_run_filter_xyt(m,s,r,k,action)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     int k;                     // Trial number, when action = 1
     int action;                // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  if (action == 0){
    model_filter_xyt_prep(m);
  }else if (action == 1){
    model_filter_xyt_get_response(m,s,r,k);
  }else if (action == -1){
    ; //model_filter_xyt_done(m->ppl);  // WYETH - clean up
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_ME_RUN_GABOR_COMP                           */
/*                                                                           */
/*  Motion energy model with Gabor (and perhaps AB temporal) filters.        */
/*                                                                           */
/*****************************************************************************/
void mod_me_run_gabor_comp(m,s,r,k,action)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     int k;                     // Trial number, when action = 1
     int action;                // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  if (action == 0){
    model_gabor_comp_prep(m);
  }else if (action == 1){
    model_gabor_comp_get_response(m,s,r,k);
  }else if (action == -1){
    ; //model_gabor_comp_done(m->ppl);  // WYETH - clean up
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_ME_RUN_DS_SIMP_01                           */
/*                                                                           */
/*  A single oriented filter, to model a DS simple cell.                     */
/*                                                                           */
/*****************************************************************************/
void mod_me_run_ds_simp_01(m,s,r,k,action)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     int k;                     // Trial number, when action = 1
     int action;                // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  if (action == 0){
    if (m->ppl == NULL)
      model_me_01_prep_o(m);
    else
      model_me_01_prep(m);
    model_me_01_prep_ft();
  }else if (action == 1){
    model_ds_simp_01_get_response(m,s,r,k);
  }else if (action == -1){
    model_me_01_done(m->ppl);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_ME_RUN_DS_LIV3                             */
/*                                                                           */
/*  Livingstone & Conway (2003).                                             */
/*                                                                           */
/*****************************************************************************/
void mod_me_run_ds_liv3(m,s,r,k,action)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     int k;                     // Trial number, when action = 1
     int action;                // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  if (action == 0){
    if (m->ppl == NULL)
      model_me_01_prep_o(m);
    else
      model_me_01_prep(m);
    model_me_01_prep_ft();
  }else if (action == 1){
    model_ds_liv3_get_response(m,s,r,k);
  }else if (action == -1){
    model_me_01_done(m->ppl);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               MOD_ME_RUN_MEM                              */
/*                                                                           */
/*  Motion energy model with Gabor (and perhaps AB temporal) filters.        */
/*                                                                           */
/*****************************************************************************/
void mod_me_run_mem(m,s,r,k,action)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     int k;                     // Trial number, when action = 1
     int action;                // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  if (action == 0){
    model_mem_prep(m);
    model_me_01_prep_ft();  // WYETH 2014 - added this line

    { // WYETH 2014 -- TEMP, replace ...
      int i;
      int fn;
      struct mod_me_fstruct **fs;

      fn = mod_me_fn;
      fs = mod_me_ff;
      for(i=0;i<mod_me_mem->nf;i++){
	mod_me_mem->fs[i][0] = fs[i*4  ]->s;
	mod_me_mem->fs[i][1] = fs[i*4+1]->s;
	mod_me_mem->fs[i][2] = fs[i*4+2]->s;
	mod_me_mem->fs[i][3] = fs[i*4+3]->s;
      }
    }

  }else if (action == 1){
    model_mem_get_response(m,s,r,k);
  }else if (action == -1){
    mylog(mylogf," *** WYETH HERE NEED TO CLEAN UP\n");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_ME_RUN_REICH                             */
/*                                                                           */
/*  Reichardt model.                                                         */
/*                                                                           */
/*****************************************************************************/
void mod_me_run_reich(m,s,r,k,action)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     int k;                     // Trial number, when action = 1
     int action;                // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  if (action == 0){
    model_reich_prep(m);
  }else if (action == 1){
    model_reich_get_response(m,s,r,k);
  }else if (action == -1){
    ; // mylog(mylogf," *** WYETH HERE NEED TO CLEAN UP\n");
  }
}
