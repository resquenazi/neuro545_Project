/*****************************************************************************/
/*                                                                           */
/*   mod_mesh_util.c                                                         */
/*   wyeth bair                                                              */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// WYETH 2019 - to remove OpenGL, X11 dependency
//#include <GL/glx.h>
//#include <GL/gl.h>  // Added for glplot_util.h GLubyte error

#include "my_util.h"
#include "nr_util.h"
#include "carray_util.h"
#include "iarray_util.h"
#include "farray_util.h"
#include "func_util.h"
#include "myrand_util.h"
#include "data_util.h"
#include "paramfile_util.h"
//#include "glplot_util.h"  // WYETH 2019 - to remove OpenGL, X11 dependency
#include "retina_util.h"
#include "stim_util.h"
#include "mod_util.h"
#include "ifc_util.h"
//#include "mod_gui_util.h" // WYETH 2019 to remove OpenGL, X11 dependency
#include "pop_util.h"
#include "mod.h"
#include "paramfile.h"
#include "ifc.h"

// For cluster computation  (static needed w/ initialization
static int myid = -1;       // Process ID, -1-Single, 0-Master, >0-Slave 
static char *mylogf = NULL; // Logfile name 
char ggstr[LONG_SLEN];      // Global log string 

// Global for mesh model
int mod_mesh_prep_flag = 0;  // 1-model has been prepped

char *mod_mesh_out_prefix;  // Prefix for output files 
int mod_mesh_reset;         // 1-reset between stimuli, default

char *mod_mesh_replay_file; // .2d filename to replay data on mosaic

int mod_mesh_type;          // Type of mesh, 0-cartesian, 1-hex/irreg

int mod_mesh_w;             // Width of mesh
int mod_mesh_h;             // Height of mesh
int mod_mesh_tmax;          // Maximum time duration of model (time units)
int mod_mesh_tn;            // Time duration for current trial
float mod_mesh_tscale;      // seconds per time unit
float mod_mesh_sscale;      // degr/spatial unit
float mod_mesh_dt;          // Time increment (time_units)
float mod_mesh_h2_dt;       // Time increment (time_units) H2 mesh
int   mod_mesh_nstep;       // Number of sub-steps per time step ()
int   mod_mesh_h2_nstep;       // Number of sub-steps per time step ()
float mod_mesh_dx;          // Space increment (1.0)
float mod_mesh_csig;        // SD Gaussian summation in center (pix)
float mod_mesh_csig_normw;  // Total BP weights should be this value
float mod_mesh_csig_min;    // Minimum value of weight to include in summation
int   mod_mesh_gcs_code;    // GC summation code, 111-SML, 010-M only, ...
int   mod_mesh_gcs_distr;   // GC summation distrib, 0-Gauss, 1-unif
float mod_mesh_gcs_prob;    // GC summation, probability 0..1
int   mod_mesh_gcs_seed;    // GC summation, seed for prob < 1.0
int   mod_mesh_gcs_n;       // GC summation, # of inputs (-1 if unset)

float mod_mesh_whc;         // Weight of H-cell signal in difference
float mod_mesh_w_h1_s;      // Weight of H1-cell signal in S-bipolar diff
float mod_mesh_w_h2_s;      // Weight of H2-cell signal in S-bipolar diff
float mod_mesh_w_h1_lm;     // Weight of H1-cell signal in L,M-bipolar diff
float mod_mesh_w_h2_lm;     // Weight of H2-cell signal in L,M-bipolar diff

float mod_mesh_off_flag;    // [0] 1-OFF cell, i.e., negate 'diff' before IFC

float mod_mesh_m_csig;      // SD Gaussian summation in center
float mod_mesh_m_csig_min;  // Minimum value of weight to include in summation
char *mod_mesh_sgen;        // Spike gen method for GC

int   mm_stim_override;     // 0, 1-override the stimulus
int   mm_stim_overdel;      // 0-no delay, 1-override at tn/2, zero before
int   mm_stim_overzci;      // -1-ignore, cone index to hold at zero
char *mm_stim_overbin;      // '0' (other CID), 'L','M','S','all','none'

int mod_mesh_gain_hpflag;   // Flag for gain-controlled HP-filtering
float *mod_mesh_lp_mask;    // Mask for low-pass filter 
int mod_mesh_lpmn;          // number of points in mask 
float mod_mesh_lptau;       // time const (time units) 
float *mod_mesh_lp2_mask;   // Mask for low-pass filter for gain2
int mod_mesh_lp2mn;         // number of points in mask 
float mod_mesh_lp2tau;      // time const (time units) 
float mod_mesh_gsig_slope;  // slope of sigmoid for lp gain 
float mod_mesh_gsig_mean;   // mean of sigmoid for lp gain 
float mod_mesh_gain_g2in;   // (0) 1-apply 'gain2' as 'gti' in IFC spike-gen
float *mod_mesh_gain1;      // Global gain, average of abs.val. over time [tn]
float *mod_mesh_gain2;      // Global gain, average of power over time [tn]

int mm_mosaic_flag;         // 0-default: gray, cartesian, 1-color random/hex
int mm_mosaic_gui;          // 1-default: show GUI for cone mosaic
int mm_mosaic_out_id;       // Cone ID for writing output
int mm_mosaic_cn;           // Number of cones in mosaic
int *mm_mosaic_cid;         // [cn]Cone ID, 0-S, 1-M, 2-L
float *mm_mosaic_cx;        // [cn] x-coord in [0..w)
float *mm_mosaic_cy;        // [cn] y-coord in [0..h)
int *mm_mosaic_cwi;         // [cn] truncated x-coord
int *mm_mosaic_cwj;         // [cn] truncated y-coord
float **mm_mosaic_cww;      // [cn][4] weights for four nearest grid points
float **mm_mosaic_cex;      // [cn][tn] cone excitation
float **mm_mosaic_gin;      // [w][h] input conductance
float mm_fsr,mm_fsg,mm_fsb; // S sensitivity to R,G,B
float mm_fmr,mm_fmg,mm_fmb; // M sensitivity to R,G,B
float mm_flr,mm_flg,mm_flb; // L sensitivity to R,G,B
float mm_mosaic_csizepix;   // Default cone size for GUI (pix)
float mm_mosaic_tsizepix;   // True cone size for GUI (pix)
float mm_mosaic_dens_mm;    // Cones per mm^2
float mm_mosaic_dens_deg;   // Cones per deg^2
float mm_mosaic_degpmm;     // Deg per mm on retina
float mm_mosaic_noise;      // Noise in cone positions, 0.0 = None
char *mm_mosaic_arrgn;      // "Spiral", "RegJit"
int   mm_mosaic_hex_w;      // Width of hex grid, for RegJit, wide central row
int   mm_mosaic_hex_h;      // Height of hex grid, for RegJit, number of rows

int    *mm_mosaic_bps_n;    // [cn] Number of BP inputs to ith RGC
int   **mm_mosaic_bps_i;    // [cn][...bps_n[]] Index of each input
float **mm_mosaic_bps_w;    // [cn][...bps_n[]] Weight of each input


//  The goal is to establish a unique list of edges and their weights:
//    mm_edge_0, ..._1, ..._w
//
//
int   **mm_edge6_cid;       // [cn][6] 6 closest edges, CID
double **mm_edge6_dist;     // [cn][6] 6 closest edges, distance [0..w),[0..h)
int    *mm_edge6_n;         // [cn] number of edges stored, 0..6
int     mm_edge_n;          // Number of unique neighbor pairs
int    *mm_edge_0;          // [mm_edge_n] CID of 1st cone in pair
int    *mm_edge_1;          // [mm_edge_n] CID of 2nd cone in pair
double *mm_edge_w;          // [mm_edge_n] weight of edge

int mod_mesh_ph_tmask_n;    // Length of mask 
float *mod_mesh_ph_tmask;   // Mask for photoreceptor time smoothing 
float mod_mesh_qno_f;       // Quantal photorecptor noise, scale factor 

int   mod_mesh_h2_flag;     // 0-just H1 mesh, 1-H1 and H2 mesh
float mod_mesh_gh; // Horizontal conductance 
float mod_mesh_gp; // Stimulus input conductance 
float mod_mesh_h2_gh; // Horizontal conductance 
float mod_mesh_h2_gp; // Stimulus input conductance 
float mod_mesh_cc; // C, Capacitance 
float *mod_mesh_hcw;        //  0-S, 1-M, 2-L, Cone-dependent H-mesh weights
float mod_mesh_h1cw[3];     // Cone weights for H1 mesh
float mod_mesh_h2cw[3];     // Cone weights for H2 mesh

float ***mod_mesh_stim; // Stimulus [t][h][w] 
float **mod_mesh_ts1;
float **mod_mesh_ts2;

float *mod_mesh_a; // Coef's for tri-diag 
float *mod_mesh_b;
float *mod_mesh_c;
float *mod_mesh_r;

float ***mod_mesh_u;      // Storage for solution of Horiz cell mesh
float **mod_mesh_uo;
float **mod_mesh_un;
float *mod_mesh_tcol;

float **mod_mesh_hirr;    // [cn][tn] Irregular (hex) horizontal mesh solution
float **mod_mesh_h2rr;    // [cn][tn] Same thing, but for the H2 Mesh

float ***mod_mesh_diff = NULL;   // response, kept for efficient repeats 
float  **mod_mesh_diffm = NULL;  // same, for mosaic

float **mod_mesh_ad1;     // Adaptation state for Stim-Horiz diff.
float  *mod_mesh_ad1m;    // [cn] Adaptation state for Stim-Horiz diff. MOSAIC
float mod_mesh_ad1_tau;   // Adaptation time constant (s)
float mod_mesh_ad1_c;     // Multiplicative constant based on ..._ad1_tau
float mod_mesh_ad1_a;     // Adapt ampl. [1.0] for incomplete, set < 1.0
int   mod_mesh_ad1_t0;    // Time index to set ad1 reference value
float mod_mesh_ad1_init_l;  // Adaptation initial value for L-BP
float mod_mesh_ad1_init_m;  // Adaptation initial value for M-BP
float mod_mesh_ad1_init_s;  // Adaptation initial value for S-BP

float mod_mesh_m_ad1_tau; // Adaptation time constant (s)
float mod_mesh_m_ad1_c;   // Multiplicative constant based on ..._ad1_tau
float mod_mesh_m_ad1_a;   // Adapt ampl. [1.0] for incomplete, set < 1.0

int   mod_mesh_disp;      // 0-no display, 1-display real-time retinal state
float mod_mesh_zoom;      // Pixel zoom for displaying 2D image in OpenGL
int mod_mesh_dispn;       // Show 2D image every n frames 
float mod_mesh_min;       // Min value in image - for normalizing display only 
float mod_mesh_max;       // Max value in image - for normalizing display only 
char *mod_mesh_outtype;   // 'diff', 'hc', what to show 

char *mod_mesh_dump_type; // Type of output to dump
char *mod_mesh_dump_file; // File name for dump
int   mod_mesh_dump_cid;  // Cone ID reference for dump
float mod_mesh_dump_tsec; // Time reference for dump (s), -1.0 for max time
int   mod_mesh_dumpp;     // Dump traces computed during 'prep' (not responses)
int   mod_mesh_dump;      // Catch-all dump flag, dump to 'mod_mesh_dumpf'
int   mod_mesh_dump_hex;  // WYETH - unify this w/ other dump commands

int   mod_mesh_slice_ti;  // Time index for slicing the H-mesh
int   mod_mesh_slice_xi;  // x-coord for slicing the H-mesh

char *mod_mesh_dumpf;     // File name for dump
char *mod_mesh_dumps;     // Stimulus name for dump file plots

/**************************************-**************************************/
/*                                                                           */
/*                    MODEL_MESH_MOSAIC_GET_IJW_FOR_CID                      */
/*                                                                           */
/*****************************************************************************/
void model_mesh_mosaic_get_ijw_for_cid(cid,rti,rtj,rw0,rw1,rw2,rw3)
     int cid;
     int *rti,*rtj;
     float *rw0,*rw1,*rw2,*rw3;
{
  *rti = mm_mosaic_cwi[cid];
  *rtj = mm_mosaic_cwj[cid];
  *rw0 = mm_mosaic_cww[cid][0];
  *rw1 = mm_mosaic_cww[cid][1];
  *rw2 = mm_mosaic_cww[cid][2];
  *rw3 = mm_mosaic_cww[cid][3];
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_MESH_RESP_GET_H                           */
/*                                                                           */
/*  Get the horizontal cell signal.                                          */
/*                                                                           */
/*****************************************************************************/
float *mod_mesh_resp_get_h(mi,xi,yi,hi)
     int mi;  // Mosiac index
     int xi;  // cartesian coordinate
     int yi;  // cartesian coordinate
     int hi;  // index of the H-mesh, e.g. 0=H1, 1=H2
{
  int i;
  int ti,tj,tn;
  float w0,w1,w2,w3,*t,***u;

  tn = mod_mesh_tn;
  u  = mod_mesh_u;

  t = get_farray(tn);

  if (mm_mosaic_flag == 1){
    if (mod_mesh_type == 0){
      if (hi != 0) // WYETH not bothering to implement H2 for old mesh type
	exit_error("MOD_MESH_RESP_GET_H","Not yet implemented for non H1");
      model_mesh_mosaic_get_ijw_for_cid(mi,&ti,&tj,&w0,&w1,&w2,&w3);
      for(i=0;i<tn;i++){
	t[i] = (w0*u[i][ti  ][tj  ] + w1*u[i][ti+1][tj  ] + // Horiz signal
		w2*u[i][ti  ][tj+1] + w3*u[i][ti+1][tj+1]);
      }
    }else{
      if (hi == 1){
	for(i=0;i<tn;i++){
	  t[i] = mod_mesh_h2rr[mi][i];
	}
      }else{
	for(i=0;i<tn;i++){
	  t[i] = mod_mesh_hirr[mi][i];
	}
      }
    }
  }else{
    exit_error("MOD_MESH_RESP_GET_H","Not imp'd yet");
  }

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_MESH_GET_MM_MOSAIC_CN                         */
/*                                                                           */
/*****************************************************************************/
int mod_mesh_get_mm_mosaic_cn(){
  return mm_mosaic_cn;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_MESH_PTR_MM_MOSAIC_CID                        */
/*                                                                           */
/*****************************************************************************/
int *mod_mesh_ptr_mm_mosaic_cid(){
  return mm_mosaic_cid;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_MESH_PTR_MM_MOSAIC_CX                         */
/*                                                                           */
/*****************************************************************************/
float *mod_mesh_ptr_mm_mosaic_cx(){
  return mm_mosaic_cx;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_MESH_PTR_MM_MOSAIC_CY                         */
/*                                                                           */
/*****************************************************************************/
float *mod_mesh_ptr_mm_mosaic_cy(){
  return mm_mosaic_cy;
}
/**************************************-**************************************/
/*                                                                           */
/*                       MOD_MESH_GET_MM_MOSAIC_DENS_MM                      */
/*                                                                           */
/*****************************************************************************/
float mod_mesh_get_mm_mosaic_dens_mm(){
  return mm_mosaic_dens_mm;
}
/**************************************-**************************************/
/*                                                                           */
/*                       MOD_MESH_GET_MM_MOSAIC_DENS_DEG                     */
/*                                                                           */
/*****************************************************************************/
float mod_mesh_get_mm_mosaic_dens_deg(){
  return mm_mosaic_dens_deg;
}
/**************************************-**************************************/
/*                                                                           */
/*                       MOD_MESH_GET_MM_MOSAIC_DEGPMM                       */
/*                                                                           */
/*****************************************************************************/
float mod_mesh_get_mm_mosaic_degpmm(){
  return mm_mosaic_degpmm;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_MESH_GET_MM_RESP                           */
/*                                                                           */
/*  Let outside routines (e.g., mod_pop_util) request response traces.       */
/*                                                                           */
/*  NOTES                                                                    */
/*  1.  New storage is created and returned, thus caller should free.        */
/*                                                                           */
/*****************************************************************************/
float *mod_mesh_get_mm_resp(xi,yi,zi,dtype,rn,rsamp)
     int xi;        // x-index, or Cone index
     int yi;        // y-index, or 0 (if mosaic)
     int zi;        // 0-OFF, 1-ON
     char *dtype;   // "photo", "horiz", "h2", "diff"
     int *rn;       // length of returned data in sampling units
     float *rsamp;  // samples/s
{
  int tn;
  float *d;

  tn = mod_mesh_tn;

  if (mm_mosaic_flag == 1){
    if ((xi >= 0) && (xi < mm_mosaic_cn) && (yi == 0)){
      if (strcmp(dtype,"photo")==0){
	d = copy_farray(mm_mosaic_cex[xi],tn);
      }else if (strcmp(dtype,"horiz")==0){
	d = mod_mesh_resp_get_h(xi,-1,-1,0);
      }else if (strcmp(dtype,"h2")==0){
	d = mod_mesh_resp_get_h(xi,-1,-1,1);
      }else if (strcmp(dtype,"diff")==0){
	d = copy_farray(mod_mesh_diffm[xi],tn);
      }else
	exit_error("MOD_MESH_GET_MM_RESP","Unknown data type");
    }else{
      exit_error("MOD_MESH_GET_MM_RESP","Bad mosaic index");
    }
  }else{
    if ((xi >= 0) && (xi < mod_mesh_w) && (yi >= 0) && (yi < mod_mesh_h)){
      if (strcmp(dtype,"photo")==0){
	d = copy_farray(mod_mesh_stim[xi][yi],tn);
      }else if (strcmp(dtype,"horiz")==0){
	d = mod_mesh_resp_get_h(-1,xi,yi,0);
      }else if (strcmp(dtype,"diff")==0){
	d = copy_farray(mod_mesh_diff[xi][yi],tn);
      }else
	exit_error("MOD_MESH_GET_MM_RESP","Unknown data type");
    }else{
      exit_error("MOD_MESH_GET_MM_RESP","Bad grid index");
    }
  }

  *rn = tn;
  *rsamp = 1.0 / mod_mesh_tscale;

  return d;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_MESH_GET_RGC_CONN                          */
/*                                                                           */
/*  Return a list of the connections and weights coming from BPs to the      */
/*  RGC 'k'.                                                                 */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_get_rgc_conn(k,zi,rn,rci,rw)
     int k;       // RGC cell index
     int zi;      // z index, 0-OFF, 1-ON  *** UNUSED CURRENTLY
     int *rn;     // Number of connections returned
     int **rci;   // [*rn]
     float **rw;  // [*rn]
{
  int nsyn;

  if (mm_mosaic_flag != 1)
    exit_error("MOD_MESH_GET_RCG_CONN","Mosaic flag not 1");

  nsyn = mm_mosaic_bps_n[k];

  *rn = nsyn;
  *rci = copy_iarray(mm_mosaic_bps_i[k],nsyn);
  *rw  = copy_farray(mm_mosaic_bps_w[k],nsyn);
}
/**************************************-**************************************/
/*                                                                           */
/*                              BPS_OVERLAP_COUNT                            */
/*                                                                           */
/*****************************************************************************/
int bps_overlap_count(g1,g2)
     int g1,g2;  // CID for two RGCs
{
  int i,j,k;
  int cnt,n1,n2,*ilist1,*ilist2;

  n1 = mm_mosaic_bps_n[g1];  // Number of BP inputs to g1
  n2 = mm_mosaic_bps_n[g2];  // Number of BP inputs to g2
  ilist1 = mm_mosaic_bps_i[g1];
  ilist2 = mm_mosaic_bps_i[g2];

  cnt = 0;
  for(i=0;i<n1;i++){
    k = ilist1[i];
    for(j=0;j<n2;j++){
      if (k == ilist2[j])
	cnt += 1;
    }
  }

  return cnt;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_MESH_MOSAIC_GET_CIRC_W                        */
/*                                                                           */
/*                                                                           */
/*          |.....                                                           */
/*          |    |**.        C = area of *'s                                 */
/*         y|    |*** .                                                      */
/*          |    ----- .                                                     */
/*          |           .                                                    */
/*          -------------                                                    */
/*               x                                                           */
/*                                                                           */
/*  If quadrant I of the unit circle is broken by a vertical line at x and   */
/*  a horizontal line at y, then compute the area of the circle that lies    */
/*  to the upper right of the point (x,y).                                   */
/*                                                                           */
/*  WYETH - THIS MIGHT BE USEFUL FOR COMPUTING pix-to-cone WEIGHTs..         */
/*                                                                           */
/*****************************************************************************/
float mod_mesh_mosaic_get_circ_w(x,y)
     float x,y;  //  0 < x,y < 1
{
  float a,b,c,r;
  float th1,th2,th3,alpha,beta;


  if ((x <= 0) || (x >= 1.0) || (y <= 0) || (y >= 1.0)){
    printf("  *** x,y = %f %f\n",x,y);
    exit_error("MOD_MESH_MOSAIC_GET_CIRC_W","Bad x,y values");
  }

  th1 = asin(y);
  th2 = atan2(y,x);
  th3 = acos(x);

  printf("th1,2,3  %f %f %f\n",th1,th2,th3);
  printf("th1,2,3  %f %f %f (deg)\n",180.0*th1/M_PI,180.0*th2/M_PI,180.0*th3/
	 M_PI);

  alpha = th3 - th2;
  beta  = th2 - th1;

  printf("alpha, beta   %f %f (deg)\n",180.0*alpha/M_PI,180.0*beta/M_PI);

  r = sqrt(x*x + y*y);

  printf("r  %f\n",r);

  a = r/2.0 * sin(alpha);
  b = r/2.0 * sin(beta);

  printf("a,b  %f %f\n",a,b);

  c = (alpha + beta)/2.0 - (a+b);

  printf("  c = %f\n",c);

  return c;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_MESH_OVERRIDE_TEST                          */
/*                                                                           */
/*****************************************************************************/
int mod_mesh_override_test(ci)
     int ci;  // Cone index [0..cn-1]
{
  int flag,cid;
  char *s;

  s = mm_stim_overbin;

  cid = mm_mosaic_cid[ci];

  flag = 0;

  if (ci == mm_stim_overzci){  // mm_stim_overzci is '-1' to ignore
    flag = 0;
  }else if (strcmp(s,"all")==0){
    flag = 1;
  }else if (strcmp(s,"none")==0){
    flag = 0;
  }else if (strcmp(s,"S")==0){
    if (cid == 0)
      flag = 1;
  }else if (strcmp(s,"M")==0){
    if (cid == 1)
      flag = 1;
  }else if (strcmp(s,"L")==0){
    if (cid == 2)
      flag = 1;
  }else if (is_int_string(s)){
    if (ci == atoi(s))
      flag = 1;
  }else{
    mylog_ival(mylogf,"*** cone index",ci);
    mylog_ival(mylogf,"*** cone ID",cid);
    mylog_cval(mylogf,"*** stim_override_binary",s);
    mylogx(mylogf,"MOD_MESH_OVERRIDE_TEST",
	   "Unknown value for 'stim_override_binary'");
  }

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_MESH_MPT_MAR                             */
/*                                                                           */
/*  Create an 'mpt' structure for the mesh model, then write it as a .mar    */
/*  file.                                                                    */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_mpt_mar(m)
     struct model_struct *m;    // Model params
{
  int i;
  int n;
  float *cx,*cy;
  struct pop_top *mpt;
  struct pop_layer *pl;

  mylog(mylogf,"  MOD_MESH_MPT_MAR\n");

  // Get 'pop_top' structure
  mpt = popu_make_pop_top(mylogf,mod_mesh_sscale,mod_mesh_tscale,
			  mod_mesh_w,mod_mesh_h,mod_mesh_tn);

  mpt->nlay = 4;  // P, H, BP, G
  mpt->lay = (struct pop_layer **)myalloc(mpt->nlay*
					  sizeof(struct pop_layer *));

  n = mm_mosaic_cn;
  cx = mm_mosaic_cx;
  cy = mm_mosaic_cy;
  mpt->lay[0] = popu_make_layer_irreg("cone",n,1,cx,cy,cx,cy);
  mpt->lay[1] = popu_make_layer_cart("horiz",mod_mesh_w,mod_mesh_h,1);
  mpt->lay[2] = popu_make_layer_irreg("bp",n,1,cx,cy,cx,cy);
  mpt->lay[3] = popu_make_layer_irreg("rgc",n,1,cx,cy,cx,cy);


  popu_mar_write(mpt,m->marfile,1);
}
/**************************************-**************************************/
/*                                                                           */
/*                          MODEL_MESH_MOSAIC_GET_W4                         */
/*                                                                           */
/*  Get weights and index values for sampling points on the cartesian        */
/*  grid, 'xn' by 'yn', for each cone.                                       */
/*                                                                           */
/*     2  3                                                                  */
/*                                                                           */
/*     0  1                                                                  */
/*                                                                           */
/*****************************************************************************/
void model_mesh_mosaic_get_w4(xn,yn,cx,cy,cn,rwi,rwj,rww)
     int xn,yn;
     float *cx,*cy;
     int cn;
     int **rwi,**rwj;
     float ***rww;
{
  int i;
  int *cwi,*cwj,xi,yi;
  float **cww,dx,dy;

  cwi = (int *)myalloc(cn*sizeof(int));
  cwj = (int *)myalloc(cn*sizeof(int));
  cww = get_2d_farray(cn,4);
  
  for(i=0;i<cn;i++){
    xi = (int)cx[i];
    yi = (int)cy[i];

    if ((xi < 0) || (xi >= xn) || (yi < 0) || (yi >= yn))
      exit_error("MODEL_MESH_MOSAIC_GET_W4","Should not happen");

    cwi[i] = xi; // Only store the coords for point '0'
    cwj[i] = yi;

    dx = cx[i] - (float)xi;
    dy = cy[i] - (float)yi;

    cww[i][0] = (1.0 - dx) * (1.0 - dy);
    cww[i][1] =        dx  * (1.0 - dy);
    cww[i][2] = (1.0 - dx) *        dy;
    cww[i][3] =        dx  *        dy;
  }
  *rwi = cwi; *rwj = cwj; *rww = cww;
}
/**************************************-**************************************/
/*                                                                           */
/*                        MODEL_MESH_MOSAIC_GET_CONEX4                       */
/*                                                                           */
/*  Get the cone excitation levels given the stimulus.                       */
/*                                                                           */
/*****************************************************************************/
void model_mesh_mosaic_get_conex4(s,tr_seed)
     struct stim_struct *s;    // Stimulus data and param pair list 
     int tr_seed;              // Trial seed
{
  int i,j,k;
  int tn,cn,ti,tj,tid,*cid,***dc;
  float w0,w1,w2,w3,*cxt,ar,ag,ab,colorw,fws,fwm,fwl,*smooth,ts,tnoise;
  float *s0,*s1,*s2,*s3,**cex,***sm,***d;
  int *si0,*si1,*si2,*si3;
  float wr,wg,wb,tr,tg,tb,oval;
  int qseed;
  float qnof;

  mylog(mylogf,"  MODEL_MESH_MOSAIC_GET_CONEX4\n");

  qseed = 29 + tr_seed + tr_seed/2;
  if (qseed > 0)
    qseed = -qseed;
  qnof = mod_mesh_qno_f;

  tn  = mod_mesh_tn;
  cn  = mm_mosaic_cn;
  cid = mm_mosaic_cid;

  if (mm_stim_override == 1){

    mylog(mylogf,"  WARNING: overriding stimulus values.\n");
    
    for(i=0;i<cn;i++){  // for each cone 
      if (mod_mesh_override_test(i)){
	oval = 1.0;
      }else
	oval = 0.0;
      cxt = mm_mosaic_cex[i];
      if (mm_stim_overdel == 0){
	for(j=0;j<tn;j++)   // for all time
	  cxt[j] = oval;
      }else{
	for(j=0;j<tn;j++){  // zeros before tn/2, override after
	  if (j >= tn/2)
	    cxt[j] = oval;
	  else
	    cxt[j] = 0.0;
	}
      }
    }

  }else if ((s->d != NULL) && (s->dc == NULL)){  // OLD WAY

    ar = paramfile_get_float_param_default(s->ppl,"amp_r",1.0);
    ag = paramfile_get_float_param_default(s->ppl,"amp_g",1.0);
    ab = paramfile_get_float_param_default(s->ppl,"amp_b",1.0);
    
    fws = mm_fsr * ar  +  mm_fsg * ag  +  mm_fsb * ab;
    fwm = mm_fmr * ar  +  mm_fmg * ag  +  mm_fmb * ab;
    fwl = mm_flr * ar  +  mm_flg * ag  +  mm_flb * ab;
    
    sprintf(ggstr,"    Weights for stimulus (s,m,l)  %f %f %f\n",fws,fwm,fwl);
    mylog(mylogf,ggstr);
    
    for(i=0;i<cn;i++){  // for each cone 
      
      tid = cid[i];
      if (tid == 0)
	colorw = fws;        // S-cone Blue 
      else if (tid == 1)
	colorw = fwm;        // M-cone Green 
      else
	colorw = fwl;        // L-cone Red 
      
      w0 = colorw * mm_mosaic_cww[i][0];
      w1 = colorw * mm_mosaic_cww[i][1];
      w2 = colorw * mm_mosaic_cww[i][2];
      w3 = colorw * mm_mosaic_cww[i][3];
      
      ti = mm_mosaic_cwi[i] + 1;  // Add 1 because stimulus data starts at 1
      tj = mm_mosaic_cwj[i] + 1;
      
      s0 = &(s->d[ti  ][tj  ][1]);  // temp ptr to stim data d[1..]...
      s1 = &(s->d[ti+1][tj  ][1]);
      s2 = &(s->d[ti  ][tj+1][1]);
      s3 = &(s->d[ti+1][tj+1][1]);
      
      cxt = mm_mosaic_cex[i];
      for(j=0;j<tn;j++){  // for all time
	
	ts = w0*s0[j] + w1*s1[j] + w2*s2[j] + w3*s3[j]; // Weighted raw stim
	
	// Noise amplitude is scaled by SQRT of signal
	ts += sqrt(ts) * qnof * nr_util_gasdev(&qseed);
	
	if (ts < 0.0)
	  ts = 0.0;
	
	cxt[j] = ts;  // stimulus signal + quantal noise
      }
    }
  }else{  // NEW WAY

    //  WYETH - this sacrifices performance to avoid representing the
    //  full 3D stimulus as separate R,G, and B arrays.  Thus, it could
    //  be made to run faster


    if (s->dc == NULL)
      exit_error("NOT COLOR","Stimulus not 3c");
    
    for(i=0;i<cn;i++){  // for each cone 

      tid = cid[i];
      if (tid == 0){         // S-cone Blue 
	wr = mm_fsr;
	wg = mm_fsg;
	wb = mm_fsb;
      }else if (tid == 1){   // M-cone Blue 
	wr = mm_fmr;
	wg = mm_fmg;
	wb = mm_fmb;
      }else{                 // L-cone Red 
	wr = mm_flr;
	wg = mm_flg;
	wb = mm_flb;
      }

      // Weights for each of the four nearest stimulus pixels, and origin
      model_mesh_mosaic_get_ijw_for_cid(i,&ti,&tj,&w0,&w1,&w2,&w3);

      w0 /= 255.0;  // Account for data scaling
      w1 /= 255.0;
      w2 /= 255.0;
      w3 /= 255.0;

      si0 = s->dc[ti  ][tj  ];  // temp ptr to stimulus data
      si1 = s->dc[ti+1][tj  ];  // For 3c data, values range from 0..255
      si2 = s->dc[ti  ][tj+1];
      si3 = s->dc[ti+1][tj+1];

      cxt = mm_mosaic_cex[i];
      for(j=0;j<tn;j++){  // for all time 

	data_util_t11_get_rgb_float(si0[j],&tr,&tg,&tb); // color components
	ts  = w0 * (wr*tr + wg*tg + wb*tb);
	data_util_t11_get_rgb_float(si1[j],&tr,&tg,&tb); // color components
	ts += w1 * (wr*tr + wg*tg + wb*tb);
	data_util_t11_get_rgb_float(si2[j],&tr,&tg,&tb); // color components
	ts += w2 * (wr*tr + wg*tg + wb*tb);
	data_util_t11_get_rgb_float(si3[j],&tr,&tg,&tb); // color components
	ts += w3 * (wr*tr + wg*tg + wb*tb);

	// Noise amplitude is scaled by SQRT of signal
	if (qnof > 0.0){
	  ts += sqrt(ts) * qnof * nr_util_gasdev(&qseed);
	}

	if (ts < 0.0)
	  ts = 0.0;

	cxt[j] = ts;  // stimulus signal + quantal noise
      }
    }
  }
  //append_farray_plot("zdump.pl","cex[0]",mm_mosaic_cex[0],tn,1);

  cex = mm_mosaic_cex;

  // Apply photoreceptor temporal filter
  if (mm_stim_override == 0){
    for(i=0;i<cn;i++){

      // Assumed that stimulus zero value holds before t=0
      smooth = convolve_with_mask_causal_x0(cex[i],tn,mod_mesh_ph_tmask,
					    mod_mesh_ph_tmask_n);
      myfree(cex[i]);
      cex[i] = smooth;
    }
  }

  if (mod_mesh_type == 0){  // Cartesian h-mesh

    // Fill in 'mod_mesh_stim[tn][xn][yn]' for resistive grid
    sm = mod_mesh_stim;
    for(i=0;i<tn;i++)
      for(j=0;j<mod_mesh_w;j++)
	for(k=0;k<mod_mesh_h;k++)
	  sm[i][j][k] = 0.0;

    //
    //  WYETH - the problem is that if sscale becomes fine w.r.t. cones,
    //  then there will be gaps in this cartesian 'sm' array.
    //
    for(i=0;i<cn;i++){  // for each cone
      model_mesh_mosaic_get_ijw_for_cid(i,&ti,&tj,&w0,&w1,&w2,&w3);
      
      cxt = cex[i];
      for(j=0;j<tn;j++){  // for all time 
	
	//  'mm_mosaic_gin' is weight

	sm[j][ti  ][tj  ] += w0 * cxt[j] / mm_mosaic_gin[ti  ][tj  ];
	sm[j][ti+1][tj  ] += w1 * cxt[j] / mm_mosaic_gin[ti+1][tj  ];
	sm[j][ti  ][tj+1] += w2 * cxt[j] / mm_mosaic_gin[ti  ][tj+1];
	sm[j][ti+1][tj+1] += w3 * cxt[j] / mm_mosaic_gin[ti+1][tj+1];
	// WYETH - March 2011 - SHOULD WE KEEP WEIGHT HERE?  Not all
	// H-cells will get sames number/weight of inputs????
      }
    }

    /*** WYETH - DEBUG  THIS dumps out the cone to H-cell weight matrix
    {
      float ***t;
      int n;
      n = mod_mesh_w;
      t = get_xyt_from_txy_3d_farray(sm,tn,n,n);
      write_3d_data("zzz.sm.3d",t,n,n,tn,4,2,1);
      exit(0);
    }
    ***/
  }else{
    ; // 'mod_mesh_stim' is not needed, 'mm_mosaic_cex' is used instead
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_MESH_MOSAIC_EDGE_ADD                         */
/*                                                                           */
/*  Add the edge to the global list, keeping the list sorted from farthest   */
/*  to nearest.                                                              */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_mosaic_edge_add(cid1,cid2,w)
     int cid1;     // CID first index
     int cid2;     // CID second index
     double w;     // distance
{
  int i,k;
  int n,posi,*cidlist,nmax;
  double *distlist;

  //  mm_edge6_dist[] has farthest at [0], and IS SORTED

  nmax = 6;

  n        = mm_edge6_n[cid1];
  distlist = mm_edge6_dist[cid1];
  cidlist  = mm_edge6_cid[cid1];

  posi = -1;
  k = 0;
  while(posi == -1){
    if (w < distlist[k]){  // This is worse
      k += 1;
      if (k == n)
	posi = k;  // This might be 'nmax', which requires moving
    }else{
      posi = k;
    }
  }

  if (n == nmax){         // Full list, put it above this spot and pop top
    if (posi == 0)           // Value was larger than all in full list
      return;                // no change
      
    for(i=1;i<posi;i++){
      distlist[i-1] = distlist[i];  // Shift values up
      cidlist[i-1]  = cidlist[i];
    }
    posi = posi-1;
    
  }else{                  // Put at this spot and move down
    for(i=n-1;i>=posi;i--){
      distlist[i+1] = distlist[i];
      cidlist[i+1] = cidlist[i];
    }
    n += 1;
  }

  distlist[posi] = w;
  cidlist[posi] = cid2;
  mm_edge6_n[cid1] = n;

  //mm_edge6_cid[i][mm_edge6_n[i]] = j;
  //mm_edge6_dist[i][mm_edge6_n[i]] = w;
  //mm_edge6_n[i] += 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_MESH_CONE_MOSAIC_EDGE                        */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_cone_mosaic_edge()
{
  int i,j,k;
  int cn,ntot,cid;
  float *cx,*cy;
  double crit_dist_sq,dx,dy,x,y,dsq;

  mylog(mylogf,"  MOD_MESH_CONE_MOSAIC_EDGE\n");

  cx = mm_mosaic_cx;
  cy = mm_mosaic_cy;

  // Use 1.6 times the distance from cone 0 to cone 1 as the critical
  // distance for finding neighboring cones.  Keep this as a squared
  // value for ease of computation.
  dx = (cx[0] - cx[1]);
  dy = (cy[0] - cy[1]);
  crit_dist_sq = 1.2*1.2*(dx*dx + dy*dy);

  cn = mm_mosaic_cn;

  mm_edge6_cid = get_2d_iarray(cn,6);   // CIDs for six closest edges
  mm_edge6_dist = get_2d_darray(cn,6);  // Distances for six closest edges
  mm_edge6_n = get_zero_iarray(cn);     // Number of edges stored

  ntot = 0;  // total number of neighbor edges (will include duplicates)
  for(i=0;i<cn;i++){
    x = cx[i];
    y = cy[i];
    for(j=0;j<cn;j++){
      if (j != i){
	dx = x - cx[j];
	dy = y - cy[j];
	dsq = dx*dx + dy*dy;
	if (dsq < crit_dist_sq){

	  //  Add it to the list, which is sorted from farthest to closest
	  mod_mesh_mosaic_edge_add(i,j,sqrt(dsq));

	}
      }
    }
    ntot += mm_edge6_n[i];
  }

  sprintf(ggstr,"    Found %d edges, including duplicates\n",ntot);
  mylog(mylogf,ggstr);

  mm_edge_0 = get_zero_iarray(ntot/2);  // CID of first cone in pair
  mm_edge_1 = get_zero_iarray(ntot/2);  // CID of first cone in pair
  mm_edge_w = get_darray(ntot/2);       // Weight of edge

  //
  //  From the edge data, find the set of unique edges
  //
  k = 0;
  for(i=0;i<cn;i++){
    for(j=0;j<mm_edge6_n[i];j++){
      cid = mm_edge6_cid[i][j];
      if (i < cid){
	mm_edge_0[k] = i;
	mm_edge_1[k] = cid;
	mm_edge_w[k] = mm_edge6_dist[i][j];
	k += 1;
      }
    }
  }
  mm_edge_n = k;

  sprintf(ggstr,"    Found %d unique neighbor pairs\n",mm_edge_n);
  mylog(mylogf,ggstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                       MOD_MESH_CONE_MOSAIC_EDGE_REGJIT                    */
/*                                                                           */
/*  Fill in the "edge6" arrays to know neighbors.                            */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_cone_mosaic_edge_regjit()
{
  int i;
  int cn,i0,i1,n0,n1,nue;
  double dist,dx,dy;

  mylog(mylogf,"  MOD_MESH_CONE_MOSAIC_EDGE_REGJIT\n");

  cn = mm_mosaic_cn;

  mm_edge6_cid  = get_2d_iarray(cn,6);  // CIDs for six closest edges
  mm_edge6_dist = get_2d_darray(cn,6);  // Distances for six closest edges
  mm_edge6_n    = get_zero_iarray(cn);  // Number of edges stored

  nue = mm_edge_n;  // Number of unique edges between neighbors
  for(i=0;i<nue;i++){
    i0 = mm_edge_0[i];  // CID for endpoints of edge
    i1 = mm_edge_1[i];

    dx = mm_mosaic_cx[i0] - mm_mosaic_cx[i1];  // Distance along edge
    dy = mm_mosaic_cy[i0] - mm_mosaic_cy[i1];
    dist = sqrt(dx*dx + dy*dy);

    n0 = mm_edge6_n[i0];  // Current number of neighbors for both endpoints
    n1 = mm_edge6_n[i1];

    mm_edge6_cid[i0][n0] = i1;
    mm_edge6_dist[i0][n0] = dist;
    mm_edge6_n[i0] += 1;

    mm_edge6_cid[i1][n1] = i0;
    mm_edge6_dist[i1][n1] = dist;
    mm_edge6_n[i1] += 1;
  }

  sprintf(ggstr,"    Found %d unique neighbor pairs\n",mm_edge_n);
  mylog(mylogf,ggstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                        MODEL_MESH_CONE_MOSAIC_PREP                        */
/*                                                                           */
/*****************************************************************************/
void model_mesh_cone_mosaic_prep(ro)
     struct onode *ro;  // retina
{
  int i,j;
  int xn,yn,tn,cn,cnt,ti,tj,ncust,seedc,seedx,*cid,out_id,regs;
  float *cx,*cy,dnoise,mupix,dens,lmr,prs,*tww,ecc;
  float distmm,degpermm,conerad,frob;
  char *mofile,*arrgn;
  struct onode *mo;  // mosaic node

  mylog(mylogf,"  MODEL_MESH_CONE_MOSAIC_PREP\n");

  mo = onode_child_get_unique(ro,"cone_mosaic");

  out_id = onode_getpar_int_exit(mo,"out_id");
  mm_mosaic_out_id = out_id;

  if (myid == -1)
    mm_mosaic_gui = onode_getpar_int_dflt(mo,"gui",1);
  else
    mm_mosaic_gui = 0;

  xn = mod_mesh_w;
  yn = mod_mesh_h;
  tn = mod_mesh_tn;

  // Given average sphere diam of 21.45 mm, there are 5.34 degr/mm. 
  degpermm = onode_getpar_flt_exit(mo,"degpermm");
  cnt = onode_getpar_int_exit(mo,"ncone");
  ecc = onode_getpar_flt_exit(mo,"ecc");
  if (ecc == -1.0){
    dens = onode_getpar_flt_exit(mo,"dens");
    sprintf(ggstr,"    Density specified as %.2f cones/mm^2\n",dens);
    mylog(mylogf,ggstr);
  }else{
    distmm = ecc/degpermm;
    dens = retutil_get_cone_density(0.95,distmm);
    sprintf(ggstr,"    Eccentricity %.2f degr implies %.2f mm\n",ecc,distmm);
    mylog(mylogf,ggstr);
    sprintf(ggstr,"    Density is %.2f cones/mm^2\n",dens);
    mylog(mylogf,ggstr);
  }
  mupix = 1000.0 * mod_mesh_sscale / degpermm;  // deg/pix / deg/mm 
  sprintf(ggstr,"    Derived microns per spatial unit = %f\n",mupix);
  mylog(mylogf,ggstr);



  mofile = onode_getpar_chr_dflt(mo,"file",(char *)NULL);
  if (mofile != NULL){
    retutil_read_mosaic(mofile,&cx,&cy,&cid,&cn);
    myfree(mofile);
  }else{

    seedx  = onode_getpar_int_exit(mo,"seed_loc");
    dnoise = onode_getpar_flt_exit(mo,"noise");
    lmr    = onode_getpar_flt_exit(mo,"lm_ratio");
    prs    = onode_getpar_flt_exit(mo,"prob_s");
    regs   = onode_getpar_int_dflt(mo,"regular_s",0);
    seedc  = onode_getpar_int_exit(mo,"seed_col");
    arrgn  = onode_getpar_chr_dflt(mo,"arrangement","Spiral");

    if (strcmp(arrgn,"Spiral")==0){
      retutil_get_mosaic_color(xn,yn,cnt,dens,mupix,dnoise,seedx,seedc,prs,lmr,
			       mylogf,&cx,&cy,&cid,&cn);
    }else if (strcmp(arrgn,"RegJit")==0){
      retutil_get_mosaic_regjit(xn,yn,cnt,dens,mupix,dnoise,seedx,seedc,prs,lmr,
				regs,mylogf,&cx,&cy,&cid,&cn,
				&mm_edge_0,&mm_edge_1,&mm_edge_w,&mm_edge_n,
				&mm_mosaic_hex_w,&mm_mosaic_hex_h);
    }
    mm_mosaic_noise = dnoise;
    mm_mosaic_arrgn = arrgn;
  }

  mm_mosaic_dens_mm  = dens;
  mm_mosaic_dens_deg = dens / (degpermm * degpermm);
  mm_mosaic_degpmm   = degpermm;


  //
  //  Allow customised L,M,S assignment by cone ID
  //
  mofile = onode_getpar_chr_dflt(mo,"custom_cid_file",(char *)NULL);
  if (mofile != NULL){
    if ((strcmp(mofile,"null")!=0) && (strcmp(mofile,"NULL")!=0)){
      ncust = retutil_custom_mosaic_cid_file(mofile,cid,cn);
      sprintf(ggstr,"    %d cones customized from file %s\n",ncust,mofile);
      mylog(mylogf,ggstr);
    }
    myfree(mofile);
  }

  //
  //  Swap L and M cones
  //
  /***
  for(i=0;i<cn;i++){
    if (cid[i] == 1)
      cid[i] = 2;
    else if (cid[i] == 2)
      cid[i] = 1;
  }
  ***/


  // Keep for run time
  mm_mosaic_cn  = cn;
  mm_mosaic_cid = cid;
  mm_mosaic_cx  = cx;
  mm_mosaic_cy  = cy;
  mm_mosaic_cex = get_2d_farray(cn,tn);

  mod_mesh_ad1m = get_zero_farray(cn); // Adaptation state for mosaic

  // Get weights to map cones onto cartesian grid
  model_mesh_mosaic_get_w4(xn,yn,cx,cy,cn,&mm_mosaic_cwi,&mm_mosaic_cwj,
			   &mm_mosaic_cww);

  if (mod_mesh_type == 0){
    // Not all grid points will have the same conductance
    mm_mosaic_gin = get_zero_2d_farray(xn,yn);
    for(i=0;i<cn;i++){  // for each cone
      ti  = mm_mosaic_cwi[i];
      tj  = mm_mosaic_cwj[i];
      tww = mm_mosaic_cww[i];

      mm_mosaic_gin[ti  ][tj  ] += tww[0];
      mm_mosaic_gin[ti+1][tj  ] += tww[1];
      mm_mosaic_gin[ti  ][tj+1] += tww[2];
      mm_mosaic_gin[ti+1][tj+1] += tww[3];
    }
  }else{
    mm_mosaic_gin = NULL;
  }

  // Color sensitivity matrix for S,M,L to R,G,B
  mm_fsr = onode_getpar_nth_flt_exit(mo,"S_rgb",0);
  mm_fsg = onode_getpar_nth_flt_exit(mo,"S_rgb",1);
  mm_fsb = onode_getpar_nth_flt_exit(mo,"S_rgb",2);
  mm_fmr = onode_getpar_nth_flt_exit(mo,"M_rgb",0);
  mm_fmg = onode_getpar_nth_flt_exit(mo,"M_rgb",1);
  mm_fmb = onode_getpar_nth_flt_exit(mo,"M_rgb",2);
  mm_flr = onode_getpar_nth_flt_exit(mo,"L_rgb",0);
  mm_flg = onode_getpar_nth_flt_exit(mo,"L_rgb",1);
  mm_flb = onode_getpar_nth_flt_exit(mo,"L_rgb",2);

  // Compute Frobenious Matrix Norm, and normalize to 1.0
  frob  = mm_fsr*mm_fsr + mm_fsg*mm_fsg + mm_fsb*mm_fsb;
  frob += mm_fmr*mm_fmr + mm_fmg*mm_fmg + mm_fmb*mm_fmb;
  frob += mm_flr*mm_flr + mm_flg*mm_flg + mm_flb*mm_flb;
  frob  = sqrt(frob);

  sprintf(ggstr,"    Original Frobenious norm:  %f\n",frob);
  mylog(mylogf,ggstr);
  sprintf(ggstr,"     %8.4f %8.4f %8.4f\n",mm_flr,mm_flg,mm_flb);
  mylog(mylogf,ggstr);
  sprintf(ggstr,"     %8.4f %8.4f %8.4f\n",mm_fmr,mm_fmg,mm_fmb);
  mylog(mylogf,ggstr);
  sprintf(ggstr,"     %8.4f %8.4f %8.4f\n",mm_fsr,mm_fsg,mm_fsb);
  mylog(mylogf,ggstr);
  mm_fsr /= frob;
  mm_fsg /= frob;
  mm_fsb /= frob;
  mm_fmr /= frob;
  mm_fmg /= frob;
  mm_fmb /= frob;
  mm_flr /= frob;
  mm_flg /= frob;
  mm_flb /= frob;
  frob  = mm_fsr*mm_fsr + mm_fsg*mm_fsg + mm_fsb*mm_fsb; // Check Frob. Norm
  frob += mm_fmr*mm_fmr + mm_fmg*mm_fmg + mm_fmb*mm_fmb;
  frob += mm_flr*mm_flr + mm_flg*mm_flg + mm_flb*mm_flb;
  frob  = sqrt(frob);
  sprintf(ggstr,"    Final Frobenious norm:  %f\n",frob);
  mylog(mylogf,ggstr);
  sprintf(ggstr,"     %8.4f %8.4f %8.4f\n",mm_flr,mm_flg,mm_flb);
  mylog(mylogf,ggstr);
  sprintf(ggstr,"     %8.4f %8.4f %8.4f\n",mm_fmr,mm_fmg,mm_fmb);
  mylog(mylogf,ggstr);
  sprintf(ggstr,"     %8.4f %8.4f %8.4f\n",mm_fsr,mm_fsg,mm_fsb);
  mylog(mylogf,ggstr);


  //
  //  For use in GUI: estimated distance in pixels
  //
  mm_mosaic_csizepix = 0.4 * sqrt(1.0/(sqrt(3.0)/2.0*dens)) * 1000.0/mupix;
  conerad = 0.5 * retutil_get_cone_diam_um(ecc,"Tyler85-OD");
  sprintf(ggstr,"    Cone diameter at eccentricity %.2f degr is %.2f um\n",ecc,
	 conerad*2.0);
  mylog(mylogf,ggstr);
  mm_mosaic_tsizepix = conerad / mupix; // Estimated true cone size


  //
  //  Compute near edges for each cone
  //
  if (mod_mesh_type == 1){

    // For the original spiral algorithm, we need to discover edges
    if (strcmp(arrgn,"Spiral")==0){
      mod_mesh_cone_mosaic_edge();
    }else if (strcmp(arrgn,"RegJit")==0){
      mod_mesh_cone_mosaic_edge_regjit();
    }
    mod_mesh_hirr = NULL;
    mod_mesh_h2rr = NULL;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_MESH_BP_SUM_PICK_PROB                        */
/*                                                                           */
/*  Picking strategy:  each unit within the spatial region is chosen with    */
/*  an independent probability.  The overall number of inputs chosen is      */
/*  random.                                                                  */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_bp_sum_pick_prob(cx,cy,pseed,ti,tw,rn,rtotw)
     float cx,cy;   // Center location of post-syn unit
     int   *pseed;  // Pointer to seed;
     int   *ti;     // [*rn] index of inputs
     float *tw;     // [*rn] weight of inputs
     int   *rn;     // Number of inputs chosen
     float *rtotw;  // Total weight
{
  int i;
  int cn,n,flag;
  float x,y,w,totw;

  cn = mm_mosaic_cn;  // Number of units

  n = 0;
  totw = 0.0;
  for(i=0;i<cn;i++){ // For each possible BP input

    x = mm_mosaic_cx[i];     // Coordinates for 'i'th mosaic index
    y = mm_mosaic_cy[i];
    w = mod_util_cone_weight(mylogf,x,y,cx,cy,mod_mesh_csig,mm_mosaic_cid[i],
			     mod_mesh_gcs_code,mod_mesh_gcs_distr);

    flag = 1;
    if (w < mod_mesh_csig_min){
      flag = 0;
    }else if (mod_mesh_gcs_prob < 1.0){
      if (myrand_util_ran2(pseed) >= mod_mesh_gcs_prob)
	flag = 0;
    }

    if (flag == 1){
      ti[n] = i;
      tw[n] = w;
      n += 1;
      totw += w;
    }
  }

  *rn = n;
  *rtotw = totw;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_MESH_BP_SUM_PICK_N                          */
/*                                                                           */
/*  Picking strategy:  of all units within the spatial regions, chose 'n'    */
/*  at random.                                                               */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_bp_sum_pick_n(cid,cx,cy,n,seed,ti,tw,rtotw,rncan)
     int    cid;    // Cone index
     float  cx,cy;  // Center location of post-syn unit
     int    n;      // Number of inputs to pick
     int    seed;   // Seed for shuffling
     int   *ti;     // [*rn] index of inputs
     float *tw;     // [*rn] weight of inputs
     float *rtotw;  // Total weight
     int   *rncan;  // Number of candiates
{
  int i;
  int cn,tn,*tlist;
  float x,y,w,totw;

  cn = mm_mosaic_cn;  // Number of units

  //
  //  Create the list of candidate units, from which to pick 'n'
  //
  tlist = (int *)myalloc(cn*sizeof(int));

  tn = 0;
  totw = 0.0;
  for(i=0;i<cn;i++){ // For each possible BP input

    x = mm_mosaic_cx[i];     // Coordinates for 'j'th mosaic index
    y = mm_mosaic_cy[i];
    w = mod_util_cone_weight(mylogf,x,y,cx,cy,mod_mesh_csig,mm_mosaic_cid[i],
			     mod_mesh_gcs_code,mod_mesh_gcs_distr);

    if (w > 0.0){  // Add to list if within region
      tlist[tn] = i;
      tn += 1;
    }
  }

  //  Exit if there are too few units in the candidate poola
  if (tn < n){
    printf("  %d candidates, %d requested (cid = %d)\n",tn,n,cid);
    exit_error("MOD_MESH_BP_SUM_PICK_N","Too few candidates");
  }

  shuffle_iarray(tlist,tn,2,seed);

  for(i=0;i<n;i++){
    ti[i] = tlist[i];
    tw[i] = 1.0;
  }
  
  *rtotw = (float)n;
  *rncan = tn;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_MESH_BP_SUM_INIT                           */
/*                                                                           */
/*  Initialise the summation of BP outputs.                                  */
/*                                                                           */
/*  *** WYETH - CURRENTLY THERE IS NO DIFFERENCE FOR 'ON' vs. 'OFF' ???      */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_bp_sum_init(){

  int i;
  int cn,n,*ti,totmembyte,seed,*shuf_seed,ncan,ncan_max;
  float cx,cy,totw,*tw;

  mylog(mylogf,"  MOD_MESH_BP_SUM_INIT\n");

  cn = mm_mosaic_cn;  // Number of units

  ti = (int *)myalloc(cn*sizeof(int));      // Temporary input array
  tw = (float *)myalloc(cn*sizeof(float));  // Temporary weight array

  // probabilistic sampling
  if ((mod_mesh_gcs_prob < 1.0) || (mod_mesh_gcs_n > 0)){
    if (mod_mesh_gcs_seed > 0.0)
      seed = -mod_mesh_gcs_seed;  // Set seed to initialize
    else
      seed = mod_mesh_gcs_seed;  // Set seed to initialize
  }

  if (mod_mesh_gcs_n > 0){
    shuf_seed = get_seeds(seed,100000,cn);
    n = mod_mesh_gcs_n;
    ncan_max = 0;
  }else
    shuf_seed = NULL;


  //
  //  Storage for inputs
  //
  mm_mosaic_bps_n = get_zero_iarray(cn);
  mm_mosaic_bps_i = (int **)myalloc(cn*sizeof(int *));
  mm_mosaic_bps_w = (float **)myalloc(cn*sizeof(float *));

  totmembyte = 0;  // Keep track of total memory used

  for(i=0;i<cn;i++){ // For each GC that sums BP inputs

    cx = mm_mosaic_cx[i];  // Set center of Gaussian at position 'mi'
    cy = mm_mosaic_cy[i];

    //  Get indices 'ti', weights 'tw' and number of inputs 'n'
    if (mod_mesh_gcs_n > 0){
      mod_mesh_bp_sum_pick_n(i,cx,cy,n,shuf_seed[i],ti,tw,&totw,&ncan);
      if (ncan > ncan_max)
	ncan_max = ncan;
    }else
      mod_mesh_bp_sum_pick_prob(cx,cy,&seed,ti,tw,&n,&totw);

    mm_mosaic_bps_n[i] = n;
    mm_mosaic_bps_i[i] = copy_iarray(ti,n);
    mm_mosaic_bps_w[i] = copy_farray(tw,n);

    totmembyte += n*8;
  }

  if (mod_mesh_gcs_n > 0){
    sprintf(ggstr,"    Max candidates:  %d\n",ncan_max);
    mylog(mylogf,ggstr);
  }

  sprintf(ggstr,"    Memory used:  %d bytes\n",totmembyte);
  mylog(mylogf,ggstr);

  myfree(ti);
  myfree(tw);
  if (shuf_seed != NULL)
    myfree(shuf_seed);
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_MESH_AD1_RESET                            */
/*                                                                           */
/*  Initialise the reference values for slow adaptation.  This can be        */
/*  done in a cone-ID dependent manner.                                      */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_ad1_reset(){

  int i,j;
  int w,h,t0;
  float *photo,*horiz,*h2,whc,aval,ws_h1,ws_h2,wlm_h1,wlm_h2;

  t0 = mod_mesh_ad1_t0;

  if (mm_mosaic_flag == 0){

    w = mod_mesh_w;
    h = mod_mesh_h;

    if ((mod_mesh_ad1_init_l == 0.0) &&
	(mod_mesh_ad1_init_m == 0.0) &&
	(mod_mesh_ad1_init_s == 0.0)){

      for(i=0;i<h;i++){
	for(j=0;j<w;j++){
	  mod_mesh_ad1[i][j] = 0.0;
	}
      }

    }else{
      mylogx(mylogf,"MOD_MESH_AD1_RESET","No CIDs for cartesian retina?");
    }
  }else if (mm_mosaic_flag == 1){

    if (t0 >= 0){
      //
      //  set to P-H diff at index for specified time
      //

      if (mod_mesh_h2_flag == 0){
	whc = mod_mesh_whc;
      }else{
	ws_h1  = mod_mesh_w_h1_s;     //   otherwise, these 4 are used.
	ws_h2  = mod_mesh_w_h2_s;
	wlm_h1 = mod_mesh_w_h1_lm;
	wlm_h2 = mod_mesh_w_h2_lm;
      }

      for(i=0;i<mm_mosaic_cn;i++){
	photo = mm_mosaic_cex[i];   // Pointer to cone signal
	horiz = mod_mesh_hirr[i];   // Pointer to H-signal

	if (mod_mesh_h2_flag == 0){
	  aval = photo[t0] - whc*horiz[t0];
	}else{
	  h2 = mod_mesh_h2rr[i];   // Pointer to H-signal
	  if (mm_mosaic_cid[i] == 0){ // This is an s-cone bipolar
	    aval = photo[t0] -  ws_h1*horiz[t0] -  ws_h2*h2[t0];    // Raw diff
	  }else{
	    aval = photo[t0] - wlm_h1*horiz[t0] - wlm_h2*h2[t0];    // Raw diff
	  }
	}
	mod_mesh_ad1m[i] = aval;
      }

    }else if ((mod_mesh_ad1_init_l == 0.0) &&
	      (mod_mesh_ad1_init_m == 0.0) &&
	      (mod_mesh_ad1_init_s == 0.0)){

      //
      //  set to 0
      //
      for(i=0;i<mm_mosaic_cn;i++)
	mod_mesh_ad1m[i] = 0.0;

    }else{

      mylog(mylogf,"    Setting cone-specific 'ad1' adaptation values.");

      //
      //  Set to user specified values
      //
      for(i=0;i<mm_mosaic_cn;i++){
	if (mm_mosaic_cid[i] == 0){
	  mod_mesh_ad1m[i] = mod_mesh_ad1_init_s;
	}else if (mm_mosaic_cid[i] == 1){
	  mod_mesh_ad1m[i] = mod_mesh_ad1_init_m;
	}else if (mm_mosaic_cid[i] == 2){
	  mod_mesh_ad1m[i] = mod_mesh_ad1_init_l;
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MODEL_MESH_RGC_01_PREP                          */
/*                                                                           */
/*****************************************************************************/
void model_mesh_rgc_01_prep(m,s,r)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
{
  int i;
  int w,h,tmax,n,nm1,auto_cent,upflag,win_w,win_h,xi;
  float *mask0,*mask1,*mask2,m1,m2,m1s,m2s,a2,alpha,cx,cy,ad1_t0;
  char fname[SLEN],xval[SLEN],yval[SLEN],tstr[SLEN],*win_title;
  double avg_dist;
  struct onode *ro,*rgo,*rmgo,*spko,*h1o,*h2o;

  myid = mylog_set_id_log(m->process_id,&mylogf);  // DO THIS FIRST

  mylog(mylogf,"  MODEL_MESH_RGC_01_PREP\n");


  mod_mesh_prep_flag = 1;

  // Real-time display params, from .rsp file
  if (m->marfile == NULL){
    mod_mesh_disp  = param_geti_dflt(r->ppl,"mesh_display_flag",0);
    mod_mesh_zoom  = param_getf_exit(r->ppl,"mesh_display_zoom");
    mod_mesh_dispn = param_getf_exit(r->ppl,"mesh_display_n");
    mod_mesh_min   = param_getf_exit(r->ppl,"mesh_display_min");
    mod_mesh_max   = param_getf_exit(r->ppl,"mesh_display_max");
    mod_mesh_replay_file = param_getc_dflt(r->ppl,"mesh_replay_file","null");
  }else{
    // Set safe params for MAR option
    mod_mesh_disp  = 0;
    mod_mesh_zoom  = 8.0;
    mod_mesh_dispn = -1;
    mod_mesh_min   = 0.0;
    mod_mesh_max   = 1.0;
    mod_mesh_replay_file = strdup("null");
  }

  // WYETH - Added Jan 2010
  if (r != NULL){
    win_w   = param_geti_dflt(r->ppl,"mesh_gui_winw",480);
    win_h   = param_geti_dflt(r->ppl,"mesh_gui_winh",480);
  }else
    win_w = win_h = 480;  // Safe values for MAR option

  // Relevant onode pointers
  ro   = onode_child_get_unique(m->o,"retina0");
  rgo  = onode_get_node_type_item_val(m->o,"pop","name","rgc");
  rmgo = onode_get_node_type_item_val(m->o,"pop","name","rgc_m");

  // Get outfile prefix
  mod_mesh_out_prefix = onode_getpar_chr_dflt(ro,"outfile","zz");
  //paramfile_get_char_param_default(r->ppl,"outfile","zz");

  // Override stimulus
  mm_stim_override = onode_getpar_int_dflt(ro,"stim_override",0);
  if (mm_stim_override == 1){
    mm_stim_overbin = onode_getpar_chr_dflt(ro,"stim_override_binary","0");
    mm_stim_overdel = onode_getpar_int_dflt(ro,"stim_override_delay",0);
    mm_stim_overzci = onode_getpar_int_dflt(ro,"stim_override_zero_ci",-1);
  }

  mod_mesh_w      = onode_getpar_int_exit(m->o,"xn");
  mod_mesh_h      = onode_getpar_int_exit(m->o,"yn");
  mod_mesh_tmax   = onode_getpar_int_exit(m->o,"tn");
  mod_mesh_tn     = mod_mesh_tmax;
  mod_mesh_tscale = onode_getpar_flt_exit(m->o,"tscale");
  mod_mesh_sscale = onode_getpar_flt_exit(m->o,"sscale");
  mod_mesh_reset  = onode_getpar_int_dflt(m->o,"mod_trial_reset",1);

  {
    float tsl,hsl;

    tsl = onode_getpar_flt_dflt(ro,"mesh_dump_slice_t",-1.0);

    if (tsl == -1.0){
      mod_mesh_slice_ti = -1;  // Do not dump a slice
      mod_mesh_slice_xi = -1;  // Do not dump a slice
    }else{
      mod_mesh_slice_ti = my_rint(tsl/mod_mesh_tscale);

      mylog_xci(mylogf,"MODEL_MESH_RGC_01_PREP","slice_t coord",
		mod_mesh_slice_ti,0,mod_mesh_tn);  // Range check w/ exit

      hsl = onode_getpar_flt_dflt(ro,"mesh_dump_slice_x",0.0);
      mod_mesh_slice_xi = mod_mesh_w/2 + my_rint(hsl/mod_mesh_sscale);

      mylog_xci(mylogf,"MODEL_MESH_RGC_01_PREP","slice_x coord",
		mod_mesh_slice_xi,0,mod_mesh_w);  // Range check w/ exit
    }
  }

  mod_mesh_dump_type = onode_getpar_chr_dflt(ro,"mesh_dump_type","null");
  mod_mesh_dump_cid  = onode_getpar_int_dflt(ro,"mesh_dump_cid",-1);
  mod_mesh_dump_tsec = onode_getpar_flt_dflt(ro,"mesh_dump_tsec",-1.0);


  mod_mesh_dumpp  = onode_getpar_int_dflt(ro,"mesh_dump_prep",0);

  mod_mesh_dump  = onode_getpar_int_dflt(ro,"mesh_dump_all",0);
  if (myid != -1)
    mod_mesh_dump = 0;  // No dumping if not single processor mode

  if (mod_mesh_dump){  // Make main dump file name
    sprintf(fname,"%s.dump.pl",mod_mesh_out_prefix);
    mod_mesh_dumpf = strdup(fname);
  }else{
    mod_mesh_dumpf = NULL;
    mod_mesh_dumps = NULL;
  }

  mod_mesh_dump_hex  = onode_getpar_int_dflt(ro,"mesh_dump_hex",0);

  mm_mosaic_flag = onode_getpar_int_dflt(ro,"cone_mosaic_flag",0);
  mod_mesh_dump_file = onode_getpar_chr_dflt(ro,"mesh_dump_file","null");

  mod_mesh_outtype = onode_getpar_chr_dflt(ro,"mesh_display_type","hc");


  // Type of horizontal mesh
  mod_mesh_type = onode_getpar_int_dflt(ro,"mesh_type",0);

  h1o  = onode_get_node_type_item_val(ro,"h_mesh","name","h1");
  if (h1o == NULL)
    exit_error("MODEL_MESH_RGC_01_PREP","Cannot find <h_mesh> object");
  mod_mesh_gh = onode_getpar_flt_exit(h1o,"gh");
  mod_mesh_gp = onode_getpar_flt_exit(h1o,"gp");
  mod_mesh_dt = onode_getpar_flt_exit(h1o,"dt");

  h2o  = onode_get_node_type_item_val(ro,"h_mesh","name","h2");
  if (h2o != NULL){
    mod_mesh_h2_flag = 1;
    mod_mesh_h2_gh = onode_getpar_flt_exit(h2o,"gh");
    mod_mesh_h2_gp = onode_getpar_flt_exit(h2o,"gp");
    mod_mesh_h2_dt = onode_getpar_flt_exit(h2o,"dt");
  }else
    mod_mesh_h2_flag = 0;

  /*
  mod_mesh_gh = onode_getpar_flt_exit(ro,"mesh_gh");
  mod_mesh_gp = onode_getpar_flt_exit(ro,"mesh_gp");
  mod_mesh_dt = onode_getpar_flt_exit(ro,"mesh_dt");
  */
  if (mod_mesh_type == 0){
    //mod_mesh_cc = onode_getpar_flt_exit(ro,"mesh_c");
    mod_mesh_cc = onode_getpar_flt_exit(h1o,"c");
    mod_mesh_nstep = -1;
    mod_mesh_h2_nstep = -1;
  }else{
    mod_mesh_cc = -1.0; // Not used by hex-mesh

    // Set 'mod_mesh_dt' to be the number of substeps
    if ((mod_mesh_dt > 1.0) || (mod_mesh_dt <= 0.0)){
      //mod_mesh_dt = 1.0;
      mod_mesh_nstep = 1;
      mod_mesh_h2_nstep = 1;
    }else{
      //mod_mesh_dt = (float)my_rint(1.0/mod_mesh_dt);

      mod_mesh_nstep = my_rint(1.0/mod_mesh_dt);
      if (mod_mesh_h2_flag == 1)
	mod_mesh_h2_nstep = my_rint(1.0/mod_mesh_h2_dt);
    }
  }

  //
  //  Read in cone weights to drive inputs onto two H-meshes
  //
  mod_mesh_h1cw[0] = onode_getpar_flt_dflt(h1o,"w_s",1.0);
  mod_mesh_h1cw[1] = onode_getpar_flt_dflt(h1o,"w_m",1.0);
  mod_mesh_h1cw[2] = onode_getpar_flt_dflt(h1o,"w_l",1.0);
  if (mod_mesh_h2_flag == 1){
    mod_mesh_h2cw[0] = onode_getpar_flt_dflt(h2o,"w_s",0.0);
    mod_mesh_h2cw[1] = onode_getpar_flt_dflt(h2o,"w_m",0.0);
    mod_mesh_h2cw[2] = onode_getpar_flt_dflt(h2o,"w_l",0.0);
  }

  mod_mesh_dx = 1.0;

  if (mod_mesh_h2_flag == 0){
    mod_mesh_whc = onode_getpar_flt_exit(ro,"mesh_whc");
  }else{
    mod_mesh_w_h1_s  = onode_getpar_flt_exit(ro,"bipolar_s_wh1");
    mod_mesh_w_h2_s  = onode_getpar_flt_exit(ro,"bipolar_s_wh2");
    mod_mesh_w_h1_lm = onode_getpar_flt_exit(ro,"bipolar_lm_wh1");
    mod_mesh_w_h2_lm = onode_getpar_flt_exit(ro,"bipolar_lm_wh2");

    // WYETH - GET RID OF THIS???
    mod_mesh_whc = 1.0; // Needed for testing...
  }


  mod_mesh_off_flag = onode_getpar_int_dflt(ro,"mesh_off_flag",0);

  tmax = mod_mesh_tmax;
  w = mod_mesh_w;
  h = mod_mesh_h;
  if (w != h) exit_error("MOD_MESH_RGC_01_PREP","width must equal height");
  n = w;

  mod_mesh_a = get_zero_farray(n); // tri-diag coef's 
  mod_mesh_b = get_zero_farray(n);
  mod_mesh_c = get_zero_farray(n);
  mod_mesh_r = get_zero_farray(n);

  if (mod_mesh_type == 0){
    mod_mesh_u = get_3d_farray(tmax,h,w);     // Solution 
    mod_mesh_uo = get_zero_2d_farray(h,w); // Old values 
    mod_mesh_un = get_zero_2d_farray(h,w); // New values 
    mod_mesh_tcol = get_zero_farray(h);    // Temporary column 
    mod_mesh_hirr = NULL;
    mod_mesh_h2rr = NULL;
  }else{
    mod_mesh_u = NULL;
    mod_mesh_uo = NULL;
    mod_mesh_un = NULL;
    mod_mesh_tcol = NULL;
    // 'mod_mesh_hirr' will be set later
  }

  mod_mesh_ad1 = get_zero_2d_farray(h,w); // Adaptation state
  mod_mesh_ad1_tau    = onode_getpar_flt_dflt(ro,"ad1_tau",0.0); // (s)
  mod_mesh_ad1_a      = onode_getpar_flt_dflt(ro,"ad1_amp",0.0); // Amp.
  mod_mesh_ad1_init_l = onode_getpar_flt_dflt(ro,"ad1_init_l",0.0);
  mod_mesh_ad1_init_m = onode_getpar_flt_dflt(ro,"ad1_init_m",0.0);
  mod_mesh_ad1_init_s = onode_getpar_flt_dflt(ro,"ad1_init_s",0.0);

  ad1_t0     = onode_getpar_flt_dflt(ro,"ad1_t0",-1.0); // (s)
  if (ad1_t0 < 0.0)
    mod_mesh_ad1_t0 = -1;
  else{
    mod_mesh_ad1_t0 = my_rint(ad1_t0 / mod_mesh_tscale);
  }


  //   Given:  x_new = c*x_old, set c=e^(-1/n) to get 1/e decline on nth step
  mod_mesh_ad1_c = exp(-1.0/(mod_mesh_ad1_tau / mod_mesh_tscale));
  sprintf(ggstr,"    Adaptation constant (ad1_c) = %f\n",mod_mesh_ad1_c);
  mylog(mylogf,ggstr);


  mod_mesh_stim = get_3d_farray(tmax,h,w);
  mod_mesh_ts1 = get_2d_farray(h,w); // Interpolated stim for 1st half-step 
  mod_mesh_ts2 = get_2d_farray(h,w); // Interpolated stim for 2nd half-step 

  // Temporal filter for cones
  //m1 = onode_getpar_flt_exit(ro,"mod_mesh_photo_m1");
  //m2 = onode_getpar_flt_exit(ro,"mod_mesh_photo_m2");
  
  if (onode_test_int(ro,"photo_delta",1)){
    nm1 = 4;  // Could be 1
    mask0 = get_zero_farray(nm1);
    mask0[0] = 1.0;
  }else{
    m1s = onode_getpar_flt_exit(ro,"mod_mesh_photo_m1"); // (sec)
    m2s = onode_getpar_flt_exit(ro,"mod_mesh_photo_m2"); // (sec)
    m1 = m1s / mod_mesh_tscale;
    m2 = m2s / mod_mesh_tscale;
    a2 = onode_getpar_flt_exit(ro,"mod_mesh_photo_a2");
    if (m2 > m1)
      nm1 = 5*(int)(0.5+m2); // Mask length 
    else
      nm1 = 5*(int)(0.5+m1); // Mask length 
    mask1 = maxwell_farray(0.0,m1,1.0,nm1);
    mask2 = maxwell_farray(0.0,m2,1.0,nm1);
    multiply_farray(mask2,nm1,a2);
    mask0 = subtract_farrays(mask1,mask2,nm1);
    norm_area_farray(mask0,nm1,1.0);
    
    if (mod_mesh_dumpp){
      sprintf(fname,"%s.prep.dump.pl",mod_mesh_out_prefix);
      append_farray_plot(fname,"maxwell_1",mask1,nm1,1);
      append_farray_plot(fname,"maxwell_2",mask2,nm1,1);
      append_farray_plot(fname,"Filter_M1-M2_v_tunits",mask0,nm1,1);
    }
    myfree(mask1); myfree(mask2);
  }
      
  mod_mesh_ph_tmask_n = nm1;
  mod_mesh_ph_tmask   = mask0;

  // Quantal noise for cones

  if (onode_test_ostr(ro,"mesh_qno_f") == 1)
    // OLD NAME - WYETH REMOVE EVENTUALLY
    mod_mesh_qno_f = onode_getpar_flt_dflt(ro,"mesh_qno_f",1.0);
  else
    mod_mesh_qno_f = onode_getpar_flt_dflt(ro,"noise_factor",1.0);

  mod_mesh_gcs_distr= onode_getpar_int_dflt(rgo,"sum_bp_distrib",0);

  // WYETH - allow deg or stimPix units (Feb 2013)
  if (onode_item(rgo,"sum_bp_dist_deg") == 1){
    mod_mesh_csig  = onode_getpar_flt_exit(rgo,"sum_bp_dist_deg");
    mod_mesh_csig /= mod_mesh_sscale;
  }else{
    mod_mesh_csig     = onode_getpar_flt_exit(rgo,"sum_bp_dist");
  }

  mod_mesh_csig_normw = onode_getpar_flt_exit(rgo,"sum_bp_normw");
  mod_mesh_csig_min   = onode_getpar_flt_dflt(rgo,"sum_bp_minw",0.02);
  mod_mesh_gcs_code   = onode_getpar_int_dflt(rgo,"sum_bp_ccode",111);
  mod_mesh_gcs_prob   = onode_getpar_flt_dflt(rgo,"sum_bp_prob",1.0);
  mod_mesh_gcs_seed   = onode_getpar_int_dflt(rgo,"sum_bp_seed",17473);
  mod_mesh_gcs_n      = onode_getpar_int_dflt(rgo,"sum_bp_n",-1);

  spko = onode_child_get_unique(rgo,"spike_gen");
  mod_mesh_sgen = onode_getpar_chr_exit(spko,"type");

  if (rmgo != NULL){

    if (onode_item(rmgo,"sum_bp_dist_deg") == 1){
      mod_mesh_m_csig  = onode_getpar_flt_exit(rmgo,"sum_bp_dist_deg");
      mod_mesh_m_csig /= mod_mesh_sscale;
    }else{
      mod_mesh_m_csig     = onode_getpar_flt_exit(rmgo,"sum_bp_dist");
    }

    mod_mesh_m_csig_min = onode_getpar_flt_dflt(rmgo,"sum_bp_minw",0.02);

    mod_mesh_m_ad1_tau = onode_getpar_flt_dflt(ro,"m_ad1_tau",0.0); // (s)
    mod_mesh_m_ad1_a   = onode_getpar_flt_dflt(ro,"m_ad1_amp",0.0); // Amp.
    mod_mesh_m_ad1_c = exp(-1.0/(mod_mesh_m_ad1_tau / mod_mesh_tscale));
  }


  //  If we are going to be appending to the ifc_util DUMP_FNAME, then remove
  //  it before we start.
  if (onode_test_ostr(rgo,"spike_dump_xi") == 1){
    xi = onode_getpar_int_exit(rgo,"spike_dump_xi");
    if (xi >= 0){
      ifc_util_rm_dump_file();
    }
  }


  // Storage for global gain signals
  mod_mesh_gain1 = get_zero_farray(mod_mesh_tmax);
  mod_mesh_gain2 = get_zero_farray(mod_mesh_tmax);

  // Temporal filter for low-pass for gain-control
  mod_mesh_gain_hpflag = onode_getpar_int_dflt(ro,"mesh_gain_hp_flag",0);
  if (mod_mesh_gain_hpflag == 1){
    mod_mesh_lptau = onode_getpar_flt_exit(ro,"mesh_lp_tau");
    mod_mesh_lptau /= mod_mesh_tscale;  // Convert from s to time units

    mod_mesh_lpmn = mod_mesh_lptau * 7.0;
    alpha = 1.0/mod_mesh_lptau;
    mod_mesh_lp_mask = alpha_farray(0.0,alpha,1.0,mod_mesh_lpmn);
    norm_area_farray(mod_mesh_lp_mask,mod_mesh_lpmn,1.0);

    if (mod_mesh_dumpp){
      sprintf(fname,"%s.prep.dump.pl",mod_mesh_out_prefix);
      append_farray_plot(fname,"exp_lp",mod_mesh_lp_mask,mod_mesh_lpmn,1);
    }

    mod_mesh_gsig_slope = onode_getpar_flt_exit(ro,"mesh_gsig_slope");
    mod_mesh_gsig_mean  = onode_getpar_flt_exit(ro,"mesh_gsig_mean");
    if (mod_mesh_dumpp){
      int cfn,i;
      float *cf,*cfx,xxn;
      
      cfn = 100;
      cf = (float *)myalloc(cfn*sizeof(float));
      cfx = (float *)myalloc(cfn*sizeof(float));
      for(i=0;i<cfn;i++){
	cfx[i] = mod_mesh_gsig_mean + 
	  5.0 * (float)(i-cfn/2)/(float)(cfn/2) / mod_mesh_gsig_slope;
	cf[i] = 0.5 + 0.5*func_sigmoid_tanh(cfx[i],mod_mesh_gsig_slope,
					    mod_mesh_gsig_mean);
      }
      sprintf(fname,"%s.prep.dump.pl",mod_mesh_out_prefix);
      //append_farray_plot(fname,"sigmoid",cf,cfn,1);
      append_farray_xy_plot(fname,cfx,cf,cfn,"sigmoid");

      myfree(cf);
    }

  }else{
    mod_mesh_lp_mask = NULL;
  }

  mod_mesh_gain_g2in = onode_getpar_int_dflt(ro,"mesh_gain_inhib_g2",0);
  if (mod_mesh_gain_g2in == 1){
    mod_mesh_lp2tau = onode_getpar_flt_exit(ro,"mesh_lp2_tau");
    mod_mesh_lp2tau /= mod_mesh_tscale;  // Convert from s to time units

    mod_mesh_lp2mn = mod_mesh_lp2tau * 7.0;
    alpha = 1.0/mod_mesh_lp2tau;
    mod_mesh_lp2_mask = alpha_farray(0.0,alpha,1.0,mod_mesh_lp2mn);
    norm_area_farray(mod_mesh_lp2_mask,mod_mesh_lp2mn,1.0);

    if (mod_mesh_dumpp){
      sprintf(fname,"%s.prep.dump.pl",mod_mesh_out_prefix);
      append_farray_plot(fname,"exp_lp2",mod_mesh_lp2_mask,mod_mesh_lp2mn,1);
    }
  }else{
    mod_mesh_lp2_mask = NULL;
  }

  if (mm_mosaic_flag == 1){
    //float wmin,wmax;  WYETH DEBUG REMOVE
    float ch;

    model_mesh_cone_mosaic_prep(ro);
    if (mod_mesh_type == 1){
      mod_mesh_hirr = get_2d_farray(mm_mosaic_cn,mod_mesh_tmax);
      if (mod_mesh_h2_flag == 1)
	mod_mesh_h2rr = get_2d_farray(mm_mosaic_cn,mod_mesh_tmax);

      ch = mod_mesh_gh * 1000.0 * (mod_mesh_tscale / mod_mesh_nstep);
      //printf("_________ch = %f\n",ch);

      //
      //  Use distance values 'mm_edge_w' to set irregular mesh weights
      //
      avg_dist = 0.0;
      for(i=0;i<mm_edge_n;i++){
	avg_dist += mm_edge_w[i];
      }
      avg_dist /= (double)mm_edge_n;  // Average weight

      /* WYETH DEBUG REMOVE
      wmin = 10000.0;
      wmax = -1.0;*/
      
      for(i=0;i<mm_edge_n;i++){
	if (mm_edge_w[i] <= 0.0)
	  mylogx(mylogf,"MODEL_MESH_RGC_01_PREP","Edge distance <= zero");

	//if (i < 100) printf("w[%d]i = %f   ",i,mm_edge_w[i]);
	//mm_edge_w[i] = mod_mesh_gh * avg_dist / mm_edge_w[i];
	mm_edge_w[i] = ch * avg_dist / mm_edge_w[i];
	//if (i < 100) printf("    w[%d]o = %f\n",i,mm_edge_w[i]);

	/* WYETH DEBUG remove
	if (mm_edge_w[i] > wmax)
	  wmax = mm_edge_w[i];
	else if (mm_edge_w[i] < wmin)
	  wmin = mm_edge_w[i];
	*/
      }
      /* WYETH DEBUG REMOVE
      printf("mod_mesh_gh = %f\n",mod_mesh_gh);
      printf("w AVG: %f\n",avg_dist);
      printf("w min: %f\n",wmin);
      printf("w max: %f\n",wmax);
      */
    }

    //
    //  Write the mosaic coordinates and CID to a file
    //
    if (strcmp(mod_mesh_dump_type,"mosaic_coord")==0){
      if (strcmp(mod_mesh_dump_file,"null")==0)  // Name of output file
	sprintf(fname,"zz.mosaic.txt");
      else
	sprintf(fname,"%s",mod_mesh_dump_file);
      
      for(i=0;i<mm_mosaic_cn;i++){
	sprintf(tstr,"%5d %d %8.2f %8.2f\n",i,mm_mosaic_cid[i],mm_mosaic_cx[i],
		mm_mosaic_cy[i]);
	append_string_to_file(fname,tstr);
      }
      printf("  Appended mosaic data to:  %s\n",fname);
      printf("\n");
      exit(0);
    }
  }else{
     mm_mosaic_out_id = -1;
  }

  //  Spatial summation of BP signals
  mod_mesh_bp_sum_init();


/*
  //
  //  To Dump a list of pairs of units for use in .rsp files or analysis
  //  
  if ((myid == -1) && (paramfile_test_param(r->ppl,"region_outfile"))){
    int i,j,k;
    int ne6,nov,ucnt,pcnt;
    float x0,x1,y0,y1,xt,yt;
    char tstr[SLEN],*regofile;
    FILE *fout;

    regofile = param_getc_exit(r->ppl,"region_outfile");
    x0 = param_getf_dflt(r->ppl,"region_stim_x0",0.0);
    x1 = param_getf_dflt(r->ppl,"region_stim_x1",0.0);
    y0 = param_getf_dflt(r->ppl,"region_stim_y0",0.0);
    y1 = param_getf_dflt(r->ppl,"region_stim_y1",0.0);

    remove_file(regofile);

    fout = my_fopen(regofile,"a","MODEL_MESH_RGC_01_PREP");

    ucnt = 0;  // Units within region
    pcnt = 0;  // Pairs within region
    for(i=0;i<mm_mosaic_cn;i++){
      xt = mm_mosaic_cx[i];
      yt = mm_mosaic_cy[i];
      if ((xt >= x0) && (xt <= x1) && (yt >= y0) && (yt <= y1)){
	ne6 = mm_edge6_n[i];  // Number of neighbors (up to 6)
	for(j=0;j<ne6;j++){
	  k = mm_edge6_cid[i][j];  // CID of neighbor
	  if (k > i){  // Only for higher CID neighbors
	    nov = bps_overlap_count(i,k);
	    sprintf(tstr,"%6d %6d %6d %3d\n",pcnt,i,k,nov);
	    fprintf(fout,"%s",tstr);
	    pcnt += 1;
	  }
	}
	ucnt += 1;
      }
    }
    fclose(fout);
    sprintf(ggstr,"    %d units within region, %d pairs\n",ucnt,pcnt);
    mylog(mylogf,ggstr);
  }

*/

  // This must occur after the cone mosaic has been defined
  mod_mesh_ad1_reset();  // Reset the slow adaptation reference value

  // Automatic Stimulus Centering (mosaic must be built first)
  auto_cent = onode_getpar_int_dflt(ro,"stim_auto_center_mi",-1);
  if ((auto_cent >= 0) && (m->marfile == NULL)){

    if (paramfile_test_param(s->ppl,"cx") &&
	paramfile_test_param(s->ppl,"cy")){

      // Set 'cx' and 'cy' to be the center value in deg
      cx = mm_mosaic_cx[auto_cent] - (float)(mod_mesh_w-1)/2.0;  // [0..w]
      cy = mm_mosaic_cy[auto_cent] - (float)(mod_mesh_h-1)/2.0;  // [0..h]
      cx *= mod_mesh_sscale;
      cy *= mod_mesh_sscale;

      /*
      printf(" mm_cx = %8.4f (pix)\n",mm_mosaic_cx[auto_cent]);
      printf(" mm_cy = %8.4f (pix)\n",mm_mosaic_cy[auto_cent]);
      printf("new cx = %8.4f (deg)\n",cx);
      printf("new cy = %8.4f (deg)\n",cy);
      */

      sprintf(xval,"%.4f",cx);
      sprintf(yval,"%.4f",cy);

      upflag = paramfile_list_update(s->ppl,"cx",xval);
      upflag = paramfile_list_update(s->ppl,"cy",yval);

      update_const_param(mylogf,"cx",xval,'\0',s->cname,s->cval,s->ctype,
			  s->ncon);
      update_const_param(mylogf,"cy",yval,'\0',s->cname,s->cval,s->ctype,
			 s->ncon);

      //exit(0);
    }else{
      mylogx(mylogf,"MODEL_MESH_RGC_01_PREP",
	     "Auto-center requires .stim cy and cy");
    }
  }

  if (m->marfile != NULL){
    mod_mesh_mpt_mar(m);
    return;
  }


  if (mm_mosaic_gui == 1){
    // *** WYETH 2019 disabled to remove OpenGL, X11 dependency
    exit_error("MODEL_MESH_RGC_01_PREP",
	       "Parameter 'gui' disabled");

    /***
    int tcn;
    int ttn;
    int nbytes;
    int tcode;
    float **rpydata;

    //
    //  WYETH - for replay of responses on mosaic GUI - new Jan 2010
    //
    if (strcmp(mod_mesh_replay_file,"null")!=0){

      if (r != NULL)
	win_title = param_getc_dflt(r->ppl,"mesh_replay_title","Replay");
      else
	win_title = strdup("Replay");
      
      read_2d_data(mod_mesh_replay_file,&rpydata,&tcn,&ttn,&nbytes,&tcode);
      if (tcn != mm_mosaic_cn){
	printf("  ncones = %d    tn = %d\n",tcn,ttn);
	exit_error("MOD_MESH_RGC_01_PREP","Replay file differs in n cones");
      }else
	printf(" **** Read replay file\n");
    }else{
      ttn = 0;
      rpydata = NULL;
      win_title = strdup("Cone Mosiac");
    }

    mod_gui_view_cones(s,mod_mesh_w,mod_mesh_h,mod_mesh_tmax,mod_mesh_sscale,
		       mod_mesh_tscale,mm_mosaic_cx,mm_mosaic_cy,
		       mm_mosaic_cid,mm_mosaic_cn,
		       mm_mosaic_csizepix,mm_mosaic_tsizepix,mm_mosaic_out_id,
		       mod_mesh_csig,mod_mesh_csig_min,mod_mesh_gcs_code,
		       mod_mesh_gcs_distr,rpydata,ttn,win_w,win_h,win_title);
    ***/
  }


  // 1=Init window
  if (mod_mesh_disp){
    exit_error("MOD_MESH_RGC_01_PREP",
	       "Parameter 'mesh_display_flag' disabled");
    // WYETH 2019 Sept 22, this is being disabled to remove OpenGL
    //   and X11 dependencies from 'wm'
    //glplot_open_show_2d(mod_mesh_uo,n,n,mod_mesh_zoom,1,-1.0,-1.0);
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                           MODEL_MESH_RGC_01_DONE                          */
/*                                                                           */
/*****************************************************************************/
void model_mesh_rgc_01_done()
{
  int w,h,tn;

  if (mod_mesh_prep_flag == 0){  // Model never prep'd, no 'mylogf' is set.
    return;  // There is no clean up to do.
  }

  mylog(mylogf,"  MODEL_MESH_RGC_01_DONE\n");


  tn = mod_mesh_tn;
  w = mod_mesh_w;
  h = mod_mesh_h;

  myfree(mod_mesh_a);
  myfree(mod_mesh_b);
  myfree(mod_mesh_c);
  myfree(mod_mesh_r);

  if (mod_mesh_type == 0){
    free_3d_farray(mod_mesh_u,tn,h,w);
    free_2d_farray(mod_mesh_uo,h);
    free_2d_farray(mod_mesh_un,h);
  }else if (mod_mesh_type == 1){
    free_2d_farray(mod_mesh_hirr,mm_mosaic_cn);
    if (mod_mesh_h2rr != NULL)
      free_2d_farray(mod_mesh_h2rr,mm_mosaic_cn);
  }

  free_3d_farray(mod_mesh_stim,mod_mesh_tmax,h,w);
  free_2d_farray(mod_mesh_ts1,h);
  free_2d_farray(mod_mesh_ts2,h);

  if (mod_mesh_dump_file != (char *)NULL)
    myfree(mod_mesh_dump_file);

  if (mod_mesh_dump_type != (char *)NULL)
    myfree(mod_mesh_dump_type);

  // 2=Close window 
  if (mod_mesh_disp){
    exit_error("MOD_MESH_RGC_01_DONE",
	       "Parameter 'mesh_display_flag' disabled");
    // WYETH 2019 disabled to remove OpenGL and X11 dependencies
    //glplot_open_show_2d((float **)NULL,0,0,0.0,2,-1.0,-1.0);
}

}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_MESH_SPIKE_GEN                            */
/*                                                                           */
/*  Ideally, all spike generation would be handled here.                     */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_spike_gen(m,dex,din,tn,seed,samp,sgen,pop_name,dumptr,vflag,
			rs,rcnt,rv,rvn)
     struct model_struct *m;    // Model params
     float *dex;                // [tn] excitatory data
     float *din;                // [tn] inhibitory data
     int tn;                    // data length
     int seed;                  // spike gen randomization seed
     float samp;                // samples/sec for dex and din
     char *sgen;                // 'ifc', 'poisson'
     char *pop_name;            // population name, to retrieve spike params

     char *dumptr;              // Name to ID dump record, or NULL for none

     int vflag;                 // 1-return 'v', 0-do not

     // int **rs;                // [*rcnt] returned spike times
     // WYETH BUG FIX Apr, 10 2013

     float **rs;                // [*rcnt] returned spike times
     int *rcnt;                 // Return number of spike generated
     float **rv;                // [*rvn] return voltage (or prob)
     int *rvn;                  // Length of returned 'rv'
{
  int i;
  struct onode *po,*spko;
  char tname[SLEN];

  // Convenient onode pointers
  po = onode_get_node_type_item_val(m->o,"pop","name",pop_name);
  if (po == NULL)
    mylog_exit(mylogf,"MOD_MESH_SPIKE_GEN  Bad pop_name.\n");
  spko = onode_child_get_unique(po,"spike_gen");


  /*
  // Get copy of the 'diff' signal for unit (i,j)
  td = mod_mesh_diff[i][j];

  if (onflag == 0)
    for(ti=0;ti<tn;ti++)
      diff[ti] = -td[ti];  // For OFF response, negative 'diff'
  else
    for(ti=0;ti<tn;ti++)
      diff[ti] = td[ti];
  */

  if (strcmp(sgen,"poisson")==0){

    //
    //  WYETH *** THE ABILITY TO DUMP using 'spike_dump_xi' was added here
    //  in Sept 2011.  HOWEVER *** This is a mess because there is no
    //  consistent way by which these responses are being dumped, or can
    //  be sent back as responses (meaning things seem to operated
    //  differently for IFC, and differently depending who is calling
    //  this routine?
    //
    if (dumptr != NULL)
      vflag = 1;
    //printf("__........._____________NOIT NUILLLLLLLLLLLLLLLLLLL\n");

    ifc_util_poisson(mylogf,m,spko,dex,tn,samp,seed,0,vflag,rs,rcnt,rv,rvn);

    if (dumptr != NULL){
      sprintf(tname,"Prob_%s",dumptr);
      printf("  *** Appending Poisson rate to 'zz.poisson.dump.pl'\n");
      append_farray_plot("zz.poisson.dump.pl",tname,*rv,*rvn,1);
    }

  }else if (strcmp(sgen,"ifc")==0){
    //dumpname_ptr = (char *)NULL;

    // This changes 'gtx' and 'gti' to reflect what was used.
    ifc_test(myid,m,pop_name,dex,din,samp,tn,seed,0,dumptr,
	     rs,rcnt,vflag,rv,rvn);
  }else{
    mylog_exit(mylogf,"MOD_MESH_SPIKE_GEN  Unknown 'sgen'\n");
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                         MODEL_MESH_PHOTO_FILTER                           */
/*                                                                           */
/*  Apply filtering to the visual stimulus to mimic photoreceptors.          */
/*  Data 'd' is [1..w][1..h][1..t].                                          */
/*                                                                           */
/*****************************************************************************/
void model_mesh_photo_filter(d,tr_seed)
     float ***d;                   // Stimulus data, all indices [1..]
     int tr_seed;                  // Trial seed
{
  int i,j,k;
  int h,w,tn,nm1,qseed;
  float *mask1,*sm,qnof,*td,*sp;
  char tstr[SLEN];

  w = mod_mesh_w;
  h = mod_mesh_h;
  tn = mod_mesh_tn;

  qseed = 29 + tr_seed + tr_seed/2;
  if (qseed > 0)
    qseed = -qseed;
  qnof = mod_mesh_qno_f;

  nm1 = mod_mesh_ph_tmask_n;
  mask1 = mod_mesh_ph_tmask;

  mylog(mylogf,"    Convolving visual stimulus w/ photoreceptor filter...\n");

  td = get_farray(tn);  // Temporary array for adding noise, smoothing

  for(i=0;i<h;i++){
    for(j=0;j<w;j++){
      
      /* if ((i==h/2) && (j==w/2)){
	 float min,max;
	 get_min_max_farray(d[i+1][j+1]+1,t,&min,&max);
	 printf("Stimulus min max at center is = %f %f\n",min,max);
	 }*/
      
      sp = &(d[i+1][j+1][1]);  // Stimulus pointer (stim data starts at 1)

      for(k=0;k<tn;k++){
	// Noise amplitude is scaled by SQRT of signal
	td[k] = sp[k] + sqrt(sp[k]) * qnof * nr_util_gasdev(&qseed);
	if (td[k] < 0.0)
	  td[k] = 0.0;
      }

      // Assumed that stimulus zero value holds before t=0
      sm = convolve_with_mask_causal_x0(td,tn,mask1,nm1);

      for(k=0;k<tn;k++){
	mod_mesh_stim[k][i][j] = sm[k];
      }
      myfree(sm);
    }
  }
  
  myfree(td);
  mylog(mylogf,"      Done.\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_MESH_SPATIAL_SUM                           */
/*                                                                           */
/*                                                                           */
/*  NOTES                                                                    */
/*  - 'mod_mesh_diffm' and 'mod_mesh_diff' must be ready BEFORE calling.     */
/*                                                                           */
/*****************************************************************************/
float *mod_mesh_spatial_sum(mi,xi,yi,pflag)
     int mi;                  // mosaic index
     int xi,yi;               // grid index (instead of 'mi')
     int pflag;               // 0-no print, 1-print
{
  int i,j,k;
  int tn,nin;
  float csig,totw,weight,cx,cy,x,y,*d,*s;
  char tname[SLEN];

  tn = mod_mesh_tn;
  csig = mod_mesh_csig;  // SD for spatial summation

  s = get_zero_farray(tn);  // Summed signal to be returned
  totw = 0.0;

  if (mm_mosaic_flag == 1){

    if (1){

      // WYETH - 2013 May 20
      // Note - 'totw' can/should be factored in during prep

      nin = mm_mosaic_bps_n[mi]; // Number of inputs to unit 'mi'

      for(i=0;i<nin;i++){         // For each input to this unit
	j      = mm_mosaic_bps_i[mi][i];  // Index of this input
	weight = mm_mosaic_bps_w[mi][i];  // Weight of this input

	d = mod_mesh_diffm[j];  // Pointer to raw diff array
	totw += weight;
	for(k=0;k<tn;k++){  // The 'k'th time point of this input to 'mi'
	  s[k] += weight * d[k];
	}
      }
    }else{
      //
      //  *** WYETH - OLD WAY - - - REMOVE
      //  *** WYETH - OLD WAY - - - REMOVE
      //
      cx = mm_mosaic_cx[mi];  // Set center of Gaussian at position 'mi'
      cy = mm_mosaic_cy[mi];
      for(j=0;j<mm_mosaic_cn;j++){ // For each mosaic index

	x = mm_mosaic_cx[j];     // Coordinates for 'j'th mosaic index
	y = mm_mosaic_cy[j];
	d = mod_mesh_diffm[j];  // Pointer to raw diff array

	//weight = mod_util_cone_weight(mylogf,x,y,cx,cy,csig,mm_mosaic_cid[mi],
	weight = mod_util_cone_weight(mylogf,x,y,cx,cy,csig,mm_mosaic_cid[j],
				      mod_mesh_gcs_code,mod_mesh_gcs_distr);

	if (weight >= mod_mesh_csig_min){

	  //printf("mi %d  j %d  w %f\n",mi,j,weight);

	  totw += weight;
	  for(i=0;i<tn;i++){
	    s[i] += weight * d[i];
	  }
	}
      }
    }
    if (pflag == 1){
      sprintf(ggstr,"  Total weight (mosiac index %d) = %f\n",mi,totw);
      mylog(mylogf,ggstr);
    }
  }else{
    for(j=0;j<mod_mesh_w;j++){
      for(k=0;k<mod_mesh_h;k++){
	//d = mod_mesh_diff[xi][yi];  // Pointer to raw diff array
	d = mod_mesh_diff[j][k];  // Pointer to raw diff array
	weight = func_2d_gaussian((float)j,(float)k,(float)xi,(float)yi,
				  csig,csig,0);
	if (weight >= mod_mesh_csig_min){
	  totw += weight;
	  for(i=0;i<tn;i++){
	    s[i] += weight * d[i];
	  }
	}
      }
    }
    if (pflag == 1){
      sprintf(ggstr,"    Total weight at (%d,%d) = %f\n",xi,yi,totw);
      mylog(mylogf,ggstr);
    }
  }

  //
  //  WYETH DEBUG
  //
  /**
     if ((mi > 3220) && (mi < 3260)){
     sprintf(tname,"sum_%d",mi);
     append_farray_plot("zdump_new.pl",tname,s,tn,1);
     }
  **/

  //multiply_farray(s,tn,1.0/totw);
  multiply_farray(s,tn,mod_mesh_csig_normw/totw);

  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_MESH_GET_RESP                            */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_get_resp(m,mi,xi,yi,seed,pop_name,onflag,rs,rsn,rv,rvn,rgx,rgx0,
		       rgn)
     struct model_struct *m;  // Model parameter pair list
     int mi;                  // mosaic index
     int xi,yi;               // grid index (instead of 'mi')
     int seed;                // Random seed
     char *pop_name;          // Population name used for IFC params
     int onflag;              // 1-on response, 0-off response
     float **rs;              // returned spikes [rsn]
     int    *rsn;             // number of spikes
     float **rv;              // returned Vm [rvn]
     int    *rvn;             // length of Vm
     float **rgx;             // returned G_ex [rgn]
     float **rgx0;            // returned G_ex [rgn] RAW version
     int    *rgn;             // length of gx
{
  int i;
  int tcnt,tvn,tn;
  float *ts,*tv,*gtx,*gti,*gtx0,samp;
  char *dmpf,tname[SLEN];

  // WYETH - THIS ASSUMES Diff has been computed already.
  // WYETH - THIS ASSUMES Diff has been computed already.

  if (mm_mosaic_flag == 1){
    if ((mi < 0) || (mi >= mm_mosaic_cn))
      mylog_exit(mylogf,"MOD_MESH_GET_RESP  Mosaic index out of bounds.\n");
    sprintf(tname,"%s_pre_GTX_mi_%d",mod_mesh_dumps,mi);
  }else{
    if ((xi < 0) || (xi >= mod_mesh_w) || (yi < 0) || (yi >= mod_mesh_h))
      mylog_exit(mylogf,"MOD_MESH_GET_RESP  Grid indices out of bounds.\n");
    sprintf(tname,"%s_pre_GTX_xy_%d_%d",mod_mesh_dumps,xi,yi);
  }

  tn   = mod_mesh_tn;
  samp = 1.0/mod_mesh_tscale; // Sampling for 'diff'

  // Determine whether to dump output for an LGN unit
  dmpf = NULL;  // Dumpfile name - WYETH - for now, no dump file for IFC call

  //  The inhibitory input will be zero, unless gain control is used below
  gti = get_zero_farray(tn);

  //
  //  Get 'gtx' from the pre-computed diff data
  //
  if (mod_mesh_csig > 0.0){
    gtx = mod_mesh_spatial_sum(mi,xi,yi,1);   // Spatial summation, 1-print
  }else{
    if (mm_mosaic_flag == 1)
      gtx = copy_farray(mod_mesh_diffm[mi],tn);    // Copy data into 'gtx'
    else
      gtx = copy_farray(mod_mesh_diff[xi][yi],tn);
  }
  if (onflag == 0)
    multiply_farray(gtx,tn,-1.0); // Negate for OFF response

  gtx0 = copy_farray(gtx,tn);  // Make a copy to return for response handover

  if (mod_mesh_dump){  // Plot the signal before scale/offset for GTX
    append_farray_plot(mod_mesh_dumpf,tname,gtx,tn,1);
  }

  //
  //  Spike Generation
  //

  // Must refresh 'gti' since it gets changed
  if (mod_mesh_gain_g2in){  // Inhibotory gain control
    for(i=0;i<tn;i++)
      gti[i] = mod_mesh_gain2[i];
  }else{
    for(i=0;i<tn;i++)
      gti[i] = 0.0; // Zero inhibitory input
  }

  //  Call either 'ifc_test' or 'ifc_util_poisson'; gtx and gti could change
  mod_mesh_spike_gen(m,gtx,gti,tn,seed,samp,mod_mesh_sgen,pop_name,NULL,0,
		     &ts,&tcnt,&tv,&tvn);

  if (dmpf != NULL)
    myfree(dmpf);

  myfree(gti);

  *rs = ts;
  *rsn = tcnt;
  *rv = tv;
  *rvn = tvn;
  *rgx = gtx;
  *rgx0 = gtx0;
  *rgn = tn;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_MESH_RESP_GET_PHAD                         */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_resp_get_phad(m,r,i)
     struct model_struct *m;        // Model data
     struct response_struct *r;     // Response data
     int i;                         // Response index
{
  int mi,xi,yi;

  //
  //
  //  WYETH - NOT IMPLEMENTED YET (started this, but found 'gtx' could be
  //                               used instead)
  //

  if (mm_mosaic_flag == 1)
    mi = r->xi[i];  // Mosaic index
  else{
    xi = r->xi[i];
    yi = r->yi[i];
    mi = mod_mesh_w * yi + xi; // derived index
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_MESH_RESP_COUNT                           */
/*                                                                           */
/*****************************************************************************/
int mod_mesh_resp_count(r)
     struct response_struct *r;     // Response data
{
  int i;
  int cnt;

  cnt = 0;
  for(i=0;i<r->n;i++){ // For each response requested

    if ((r->rformat[i] >= 0) && (r->rformat[i] <= 3)){

      if ((strcmp(r->datid[i],"gc_spikes")==0) ||
	  (strcmp(r->datid[i],"gc_vm")==0) ||
	  (strcmp(r->datid[i],"gc_gx_raw")==0) ||
	  (strcmp(r->datid[i],"horiz")==0) ||
	  (strcmp(r->datid[i],"photo")==0) ||
	  (strcmp(r->datid[i],"gc_gx")==0)){

	cnt += 1;
	
	// WYETH - THIS MAY OVER-COUNT RESPONSES ?
	//   1. same cell may have multiple responses, but needs only 1 seed
	//   2. cells from other pop's might have same 'datid' names?
	
      }
    }
  }
  return cnt;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_MESH_RESP_HANDOVER                         */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_resp_handover(m,r)
     struct model_struct *m;        // Model data
     struct response_struct *r;     // Response data
{
  int i;
  int n,cnt,mi,tsn,tvn,tgn,**sn,**vn,**gn,xi,yi,zi,zn,si,*seedlist,nseed;
  float ***s,***v,***g,***g0,*ts,*tv,*tg,*tg0,*t;
  char *pop_name;

  zn = 2;  // For ON and OFF at each location

  if (mm_mosaic_flag == 1) // Set 'n' to be the number of cones
    n = mm_mosaic_cn;
  else
    n = mod_mesh_w * mod_mesh_h;

  // Each element will be NULL (or 0 for int) until its response is computed
  s = get_3d_farray(n,zn,0);     // pointers to NULL
  v = get_3d_farray(n,zn,0);     // pointers to NULL
  g = get_3d_farray(n,zn,0);     // pointers to NULL
  g0 = get_3d_farray(n,zn,0);     // pointers to NULL
  sn = get_zero_2d_iarray(n,zn);
  vn = get_zero_2d_iarray(n,zn);
  gn = get_zero_2d_iarray(n,zn);

  // Get enough randomization seeds for all responses
  nseed = mod_mesh_resp_count(r);
  seedlist = get_seeds(m->mseed[r->tsi],1000000,nseed);

  //
  //  Check each entry in the response structure to see if it is retinal.
  //  If so, then store the appropriate data for the request.
  //
  cnt = 0;
  si = 0;   // Seed index
  for(i=0;i<r->n;i++){ // For each response requested

    //  Get a pointer to the relevant cell
    if ((r->rformat[i] >= 0) && (r->rformat[i] <= 3)){

      if ((strcmp(r->datid[i],"gc_spikes")==0) ||
	  (strcmp(r->datid[i],"gc_vm")==0) ||
	  (strcmp(r->datid[i],"gc_gx_raw")==0) ||
	  (strcmp(r->datid[i],"horiz")==0) ||
	  (strcmp(r->datid[i],"photo")==0) ||
	  (strcmp(r->datid[i],"gc_gx")==0)){
	// WYETH - WHEN ADDING HERE, seed "mod_mesh_resp_count" ABOVE ****

	pop_name = r->plname[i];

	if (mm_mosaic_flag == 1)
	  mi = r->xi[i];  // Mosaic index
	else{
	  xi = r->xi[i];
	  yi = r->yi[i];
	  mi = mod_mesh_w * yi + xi; // derived index
	}

	zi = r->zi[i];  // 'zi' acts as 'onflag'

	// Check if we need to compute responses for location 'mi'
	if ((g[mi][zi] == NULL) && 
	    ((strcmp(r->datid[i],"gc_spikes")==0) ||
	     (strcmp(r->datid[i],"gc_vm")==0) ||
	     (strcmp(r->datid[i],"gc_gx_raw")==0) ||
	     (strcmp(r->datid[i],"gc_gx")==0))){

	  // Note, don't let 'horiz' or 'photo' in here, because they
	  // have population names that don't exist as 'pop' objects

	  mod_mesh_get_resp(m,mi,xi,yi,seedlist[si],pop_name,zi,
			    &ts,&tsn,&tv,&tvn,&tg,&tg0,&tgn);
	  s[mi][zi] = ts;
	  v[mi][zi] = tv;
	  g[mi][zi] = tg;
	  g0[mi][zi] = tg0;
	  sn[mi][zi] = tsn;
	  vn[mi][zi] = tvn;
	  gn[mi][zi] = tgn;
	  cnt += 1;      // How many responses have been computed
	  si += 1;
	}

	if (strcmp(r->datid[i],"gc_spikes")==0){
	  mod_util_resp_store_s(r,i,s[mi][zi],sn[mi][zi],1,mylogf); // 1-cpflag
	}else if (strcmp(r->datid[i],"gc_vm")==0){
	  mod_util_resp_store_f(r,i,v[mi][zi],vn[mi][zi],1,mylogf);
	}else if (strcmp(r->datid[i],"gc_gx")==0){
	  mod_util_resp_store_f(r,i,g[mi][zi],gn[mi][zi],1,mylogf);
	}else if (strcmp(r->datid[i],"gc_gx_raw")==0){
	  mod_util_resp_store_f(r,i,g0[mi][zi],gn[mi][zi],1,mylogf);
	}else if (strcmp(r->datid[i],"horiz")==0){
	  t = mod_mesh_resp_get_h(mi,xi,yi,0);
	  mod_util_resp_store_f(r,i,t,mod_mesh_tn,0,mylogf); // 0-do not copy
	}else if (strcmp(r->datid[i],"photo")==0){
	  if (mm_mosaic_flag == 1)
	    t = mm_mosaic_cex[mi];
	  else
	    t = mod_mesh_stim[xi][yi];
	  mod_util_resp_store_f(r,i,t,mod_mesh_tn,1,mylogf); // 0-do not copy
	}

	/*  WYETH - haven't implemented these yet
      }else if (strcmp(r->datid[i],"phdiff_ad")==0){
	mod_mesh_resp_get_phad(m,r,i);
	*/
      }else{
	sprintf(ggstr,"  *** DataID:  %s\n",r->datid[i]);
	mylog(mylogf,ggstr);
	mylog_exit(mylogf,"MOD_MESH_RESP_HANDOVER  Unknown data ID.\n");
      }
    }
  }

  sprintf(ggstr,"    Responses computed for %d locations\n",cnt);
  mylog(mylogf,ggstr);

  free_3d_farray_null(s,n,zn,0);
  free_3d_farray_null(v,n,zn,0);
  free_3d_farray_null(g,n,zn,0);
  free_3d_farray_null(g0,n,zn,0);

  free_2d_farray(sn,n);
  free_2d_farray(vn,n);
  free_2d_farray(gn,n);
}
/**************************************-**************************************/
/*                                                                           */
/*                                MOD_MESH_RUN                               */
/*                                                                           */
/*  CONVENTIONS:                                                             */
/*  - first spatial subscript selects row                                    */
/*                                                                           */
/*  Notation                                                                 */
/*    u - old value at j,k                                                   */
/*    uN - old value at j,k+1                                                */
/*    uS - old value at j,k-1                                                */
/*    uE - old value at j+1,k                                                */
/*    uW - old value at j-1,k                                                */
/*    U  - new value at j,k                                                  */
/*    UN - ...                                                               */
/*    p  - visual stimulus on old step                                       */
/*    P  - visual stimulus on new step                                       */
/*                                                                           */
/*    gh - horiz. conductance                                                */
/*    gp - stimu. conductance                                                */
/*    dt - time step                                                         */
/*    dx - spatial step                                                      */
/*                                                                           */
/*  Use ADI Crank-Nicholson method (NumRecInC V2 p855-856).                  */
/*  ADI = Alternating Direction Implicit Method                              */
/*        is a Finite Difference Method                                      */
/*                                                                           */
/*  First half-step:                                                         */
/*                                                                           */
/*    U-u     1   [ a*gh (UE-2U+UW + uN-2u+uS)                ]              */
/*  C ---  =  - * [ ---- --------------------- + gp*(p+P-u-U) ]              */
/*     dt     2   [  2           dx^2                         ]              */
/*                                                                           */
/*  let  A= dt/(2C)                                                          */
/*       B= a*gh/(2*dx*dx)                                                   */
/*                                                                           */
/*  -AB*UE + (1+2AB+Agp)*U - AB*UW = AB*(uN+uS) + (1-2AB-Agp)*u + Agp*(p+P)  */
/*                                                                           */
/*  Second half-step:                                                        */
/*                   ... (uE-2u+uW + UN-2U+US) ...                           */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_run()
{
  int i,j;
  int tn,ti,n,nm1,nm2,cnt;
  float alpha,alph2,beta,beta1m,gamma,cap,gh_prime;
  float *a,*b,*c,*r,dx,gh,gp,**uo,**un,**ut,*tpf1,*tpf2;
  double dt,tf,tfrac1,tfrac2;
  float *tu,*tun,*tus,*tcol,***u,***s,**ts1,**ts2,*tts1,*tts2;

  mylog(mylogf,"  MOD_MESH_RUN\n");
  //system("date");
  //printf("xn = %d\n",mod_mesh_w);

  n = mod_mesh_w;       // simulation size, should = mod_mesh_h
  tn = mod_mesh_tn;

  u = mod_mesh_u;       // solution
  uo = mod_mesh_uo;     // old values
  un = mod_mesh_un;     // new values
  tcol = mod_mesh_tcol; // temporary column

  s = mod_mesh_stim;    // stimulus 
  ts1 = mod_mesh_ts1;   // Interpolated stim. for 1st half-step 
  ts2 = mod_mesh_ts2;   // Interpolated stim. for 2nd half-step 

  gh = mod_mesh_gh;     // Conductance in horizontal mesh (?) 
  gp = mod_mesh_gp;     // Conductance to photoreceptor (?) 
  cap = mod_mesh_cc;    // Conductance to photoreceptor (?) 
  dt = (double)mod_mesh_dt; // Time step (time units) 
  dx = mod_mesh_dx;

  for(i=0;i<n;i++)      // Start w/ the old value of the mesh at t=0
    for(j=0;j<n;j++)    // WYETH - NEW Mar 25, 2009, remove advance in H sigal
      u[0][i][j] = uo[i][j];
  
  
  gh_prime = gh / (2.0*cap);
  
  gamma = dt * gp/cap;
  alpha = gh_prime * dt/(dx*dx);
  alph2 = alpha/2.0;
  beta = alpha + gp*dt/(2.0*cap);
  beta1m = 1.0 - beta;

  nm1 = n-1;
  nm2 = n-2;

  // Make tri-diag coef arrays - these never change.
  a = mod_mesh_a; // Arrays of length "n" = "mod_mesh_w"
  b = mod_mesh_b;
  c = mod_mesh_c;
  r = mod_mesh_r;
  for(i=0;i<n;i++){ // This would only have to go to nm2, I think
    a[i] = c[i] = -alpha/2.0;
    b[i] = 1.0 + beta;
  }
  c[0] *= 2.0; //  Account for edge
  a[nm2-1] *= 2.0; // Account for edge, only nm2 are used

  ti = 0;   // Simulation time elapsed, integer 
  tf = 0.0; // Simulation time elapsed, float 
  cnt = 0;  // Count each step 
  while(ti < tn){
    // Compute stimuli for 1st and 2nd half steps.
    // ts1 - avg of old and new stim for 1st half step 
    // ts2 - avg of old and new stim for 2nd half step 
    tfrac1 = tf - (double)ti + dt/4.0;   // Between 0 and 0.5 
    tfrac2 = tfrac1 + dt/2.0; // Between 0.5 and 1 
    for(i=0;i<n;i++){
      tpf1 = s[ti][i]; // temp pointers for speed 
      if (ti >= (tn-1)) // So stim only has to be 't' long 
	tpf2 = s[ti][i];
      else
	tpf2 = s[ti+1][i];
      tts1 = ts1[i];
      tts2 = ts2[i];
      for(j=0;j<n;j++){
	tts1[j] = tpf1[j] + tfrac1 * (tpf2[j] - tpf1[j]);
	tts2[j] = tpf1[j] + tfrac2 * (tpf2[j] - tpf1[j]);
      }
    }
    
    /*******************/
    /* FIRST HALF STEP */
    /*******************/
    for(i=1;i<nm1;i++){ // For each row 
      tu = uo[i];
      tun = uo[i-1];
      tus = uo[i+1];
      tts1 = ts1[i];
      for(j=1;j<nm1;j++) // Go along each row 
	r[j] = beta1m*tu[j] + alph2*(tun[j] + tus[j]) + gamma*tts1[j];
      tridag(a-1,b-1,c-1,r,un[i],nm2); // Solve tri-diag; r, un begin 1
      // Reflect into left and right edges
      un[i][0] = un[i][2];
      un[i][n-1] = un[i][n-3];
    }
    for(j=0;j<n;j++){ // Reflect into top and bottom edges
      un[0][j] = un[2][j];
      un[n-1][j] = un[n-3][j];
    }
    
    ut = uo; uo = un; un = ut; // Swap old and new 
    /********************/
    /* SECOND HALF STEP */
    /********************/
    for(j=1;j<nm1;j++){ // For each column
      for(i=1;i<nm1;i++){ // For each element down the column
	r[i] = beta1m*uo[i][j] + alph2*(uo[i][j-1] + uo[i][j+1]) +
               gamma*ts2[i][j];
      }
      tridag(a-1,b-1,c-1,r,tcol,nm2); // Solve tri-diag
      for(i=1;i<=nm2;i++) // Put column back in 2D array
	un[i][j] = tcol[i];
      un[0][j] = un[2][j]; // Reflect into top of column
      un[n-1][j] = un[n-3][j]; // Reflect into bottom of column
    }
    for(i=0;i<n;i++){ // Reflect into left and right edges
      un[i][0] = un[i][2];
      un[i][n-1] = un[i][n-3];
    }

    tf += dt;
    if (tf >= (float)(ti+1)){

      ti += 1; // Integeter time units elapsed

      if (ti%1000==0){
	sprintf(ggstr,"%8d time units\n",ti);
	mylog(mylogf,ggstr);
      }

      // WYETH - CHANGED Mar 25th, 2009 SO THAT 'u' is NOT Advanced !!!
      // WYETH - Mar 25th, 2009, note that u had EARLIER storage than stim
      if (ti < tn){
	for(i=0;i<n;i++) // Save the state, data for t stored at t-1
	  for(j=0;j<n;j++)
	    u[ti][i][j] = un[i][j];
      }

    }

    if (mod_mesh_disp){
      exit_error("MOD_MESH_RUN",
		 "Parameter 'mesh_display_flag' disabled");
      // WYETH 2019 disabled to remove OpenGL and X11 dependencies
      /*
      if ((mod_mesh_dispn > 0)&&(cnt % mod_mesh_dispn)==0){
	if (strcmp("diff",mod_mesh_outtype)==0){  // Show difference
	  float **diff;
	  diff = subtract_2d_farrays(ts1,un,n,n);
	  glplot_open_show_2d(diff,n,n,4.0,0,-1.0,-1.0);
	  free_2d_farray(diff,n);
	}else if (strcmp("hc",mod_mesh_outtype)==0){ // Show Horiz Cell
	  glplot_open_show_2d(un,n,n,mod_mesh_zoom,0,
			      mod_mesh_min,mod_mesh_max);
	}
      }
      */
    }
      
    ut = uo; uo = un; un = ut; // Swap old and new
    cnt += 1;
  }

  sprintf(ggstr,"%8d time units simulated.\n",ti);
  mylog(mylogf,ggstr);

  //system("date");

  // Dump 'u' to 3d data
  if ((strcmp(mod_mesh_dump_type,"3d")!=0) &&
      (strcmp(mod_mesh_dump_file,"null")!=0)){
    float ***dxyt;

    // Dump 'u' to 3d data
    //  u == mod_mesh_u     [t][h][w]  H-cell solution
    //  s == mod_mesh_stim  [t][h][w]  stimulus

    if (strcmp("diff",mod_mesh_outtype)==0){  // stim - hcell difference
      float ***tdiff;

      //
      //  WYETH - DOES THIS MAKE ANY SENSE FOR COLOR STIMULI?
      //  what is being subtracted if the stimulus is in color?
      //
      printf(" WYETH - This was developed, but never tested, does it work?\n");
      exit(0);


      tdiff = subtract_3d_farrays(s,u,0,tn,0,n,0,n);
      dxyt = get_xyt_from_txy_3d_farray(tdiff,tn,n,n);
      free_3d_farray(tdiff,tn,n,n);
      printf("  DUMPING 'diff':  stim - h-cell\n");
    }else{
      dxyt = get_xyt_from_txy_3d_farray(u,tn,n,n);
      printf("  DUMPING 'hc':  h-cell response\n");
    }
    write_3d_data(mod_mesh_dump_file,dxyt,n,n,tn,4,2,1);
    free_3d_farray(dxyt,n,n,tn);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_MESH_RUN_IRREG                            */
/*                                                                           */
/*  Compute the H-cell signal using a simple difference algorithm across     */
/*  edges.  The values from 'mm_mosaic_cex' are used, rather than the        */
/*  'mod_mesh_stim' values used by the cartesian mesh.                       */
/*                                                                           */
/*  'mod_mesh_hirr' will be set at the end.                                  */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_run_irreg(hv,gp,gh,dt,nstep)
     float **hv;
     float gp,gh,dt;   // Conductances and time step
     int nstep;        // Number of time steps
{
  int i;
  int ti,tt,cn,tn,tim1,ne,ci0,ci1;
  float ch,cp,diff,cval,cw_cp;
  float **cex;       // [cn][tn]  Cone excitation
  float **thv;       // [cn][nstep]  Horizontal cell value

  mylog(mylogf,"  MOD_MESH_RUN_IRREG\n");

  cn  = mm_mosaic_cn;    // Number of cones
  tn  = mod_mesh_tn;     // Duration
  ne  = mm_edge_n;       // Number of edges
  cex = mm_mosaic_cex;   // [cn][tn] Cone exitation

  // Values are given in reference to ms, so multiply by 1000 for sec.
  // Think of this as:  (change/sec) * (sec)

  ch = gh * 1000.0 * (mod_mesh_tscale / nstep);
  cp = gp * 1000.0 * (mod_mesh_tscale / nstep);

  sprintf(ggstr,"    ch,cp  %f %f   dt %f  nstep %d\n",ch,cp,dt,
	  mod_mesh_nstep);
  mylog(mylogf,ggstr);

  thv = get_2d_farray(cn,nstep);  // Storage in between time steps

  //
  //  For all time
  //
  for(ti=1;ti<tn;ti++){

    tim1 = ti-1;

    for(tt=0;tt<nstep;tt++){  // Use multiple steps to get next value

      //
      //  Add contribution from cone that drives this mesh point.
      //  Note, cone at 'ti' is compared to H at used to drive 'hv' at 'ti'.
      //
      if (tt == 0){  // FIRST small STEP
	for(i=0;i<cn;i++){  // Each cone drives an H-cell value
	  cw_cp = cp * mod_mesh_hcw[mm_mosaic_cid[i]];  // Cone-dependent factor
	  diff = cw_cp * (cex[i][ti] - hv[i][tim1]);
	  thv[i][tt] = hv[i][tim1] + diff;
	}
      }else{
	// Note, interpolating 'cex' changes little here, thus not done.
	for(i=0;i<cn;i++){  // Each cone drives an H-cell value
	  cw_cp = cp * mod_mesh_hcw[mm_mosaic_cid[i]];  // Cone-dependent factor
	  diff = cw_cp * (cex[i][ti] - thv[i][tt-1]);
	  thv[i][tt] = thv[i][tt-1] + diff;
	}
      }


      //
      //  Add contributions from all edges
      //
      if (mm_mosaic_noise != 0.0){
	if (tt == 0){  // FIRST small STEP
	  for(i=0;i<ne;i++){   // For each edge
	    ci0 = mm_edge_0[i];
	    ci1 = mm_edge_1[i];
	    diff = mm_edge_w[i] * (hv[ci0][tim1] - hv[ci1][tim1]);
	    hv[ci0][ti] -= diff;
	    hv[ci1][ti] += diff;
	  }
	}else{
	  for(i=0;i<ne;i++){   // For each edge
	    ci0 = mm_edge_0[i];
	    ci1 = mm_edge_1[i];
	    diff = mm_edge_w[i] * (thv[ci0][tt-1] - thv[ci1][tt-1]);
	    thv[ci0][tt] -= diff;
	    thv[ci1][tt] += diff;
	  }
	}
      }else{
	//
	//  Same as above, but for CONSTANT 'ch' = 'mod_mesh_gh'
	//
	if (tt == 0){  // FIRST small STEP
	  for(i=0;i<ne;i++){   // For each edge
	    ci0 = mm_edge_0[i];
	    ci1 = mm_edge_1[i];
	    diff = ch * (hv[ci0][tim1] - hv[ci1][tim1]);
	    hv[ci0][ti] -= diff;
	    hv[ci1][ti] += diff;
	  }
	}else{
	  for(i=0;i<ne;i++){   // For each edge
	    ci0 = mm_edge_0[i];
	    ci1 = mm_edge_1[i];
	    diff = ch * (thv[ci0][tt-1] - thv[ci1][tt-1]);
	    thv[ci0][tt] -= diff;
	    thv[ci1][tt] += diff;
	  }
	}
      }
    }

    for(i=0;i<cn;i++){
      hv[i][ti] = thv[i][nstep-1];   // Store the final value for this 'ti'
    }
  }

  free_2d_farray(thv,cn);

  // Convert to cartesian 3D


  if ((strcmp(mod_mesh_dump_type,"h_v_dist")==0) ||
      (strcmp(mod_mesh_dump_type,"h2_v_dist")==0)){ // WYETH NOT RIGHT??
    int cid,t0;
    char fname[SLEN];

    cid = mod_mesh_dump_cid;       // Cone ID for reference point

    // Set time index 't0' for analysis
    if (mod_mesh_dump_tsec == -1.0)  // -1 indicates Last time point
      t0 = tn - 1;
    else{
      t0 = my_rint(mod_mesh_dump_tsec/mod_mesh_tscale);  
      if ((t0 < 0) || (t0 >= tn)){
	printf("mod_mesh_dump_tsec = %f\n",mod_mesh_dump_tsec);
	printf("tn = %d   tscale = %f\n",tn,mod_mesh_tscale);
	exit_error("MOD_MESH_RUN_IRREG","  mod_mesh_dump_tsec out of bounds");
      }
    }
    if (strcmp(mod_mesh_dump_file,"null")==0)  // Name of output file
      sprintf(fname,"zz.dist.pl");
    else
      sprintf(fname,"%s",mod_mesh_dump_file);

    retutil_val_dist_plot(hv,cn,mm_mosaic_cx,mm_mosaic_cy,cid,t0,fname);

    printf("  Exiting after writing dump file, %s\n\n",fname);
    exit(0);
  }else if (strcmp(mod_mesh_dump_type,"time")==0){
    int cid;
    char fname[SLEN],pname[SLEN];

    cid = mod_mesh_dump_cid;

    if (strcmp(mod_mesh_dump_file,"null")==0)  // Name of output file
      sprintf(fname,"zz.cid_%d_v_time.pl",cid);
    else
      sprintf(fname,"%s",mod_mesh_dump_file);

    sprintf(pname,"cid_%d",cid);

    append_farray_plot(fname,pname,hv[cid],tn,1);

    printf("  Exiting after writing dump file, %s\n\n",fname);
    exit(0);
  }

  if ((strcmp(mod_mesh_dump_type,"3d_h")==0) ||  // new
      (mod_mesh_dump_hex == 1)){                 // old condition
    int j;
    float ***d3;
    int xn,yn,cid;
    float *xcut;

    append_farray_plot("zz.hex.cut.pl","hex_0",hv[0],tn,1);
    if (cn == 1613){
      append_farray_plot("zz.hex.cut.pl","hex_37",hv[37],tn,1);
      append_farray_plot("zz.hex.cut.pl","hex_1575",hv[1575],tn,1);
      append_farray_plot("zz.hex.cut.pl","hex_1612",hv[1612],tn,1);
      append_farray_plot("zz.hex.cut.pl","hex_806_CTR",hv[806],tn,1);
    }else if (cn == 6481){
      // 3203 ... 3277  along middle row
      append_farray_plot("zz.hex.cut.pl","hex_73",hv[73],tn,1);
      append_farray_plot("zz.hex.cut.pl","hex_6407",hv[6407],tn,1);
      append_farray_plot("zz.hex.cut.pl","hex_6480",hv[6480],tn,1);
      append_farray_plot("zz.hex.cut.pl","hex_3240_CTR",hv[3240],tn,1);
    }

    if (strcmp(mm_mosaic_arrgn,"RegJit")==0){
      //exit_error("MOD_MESH_RUN_IRREG","  Not yet imp'd Cartesian here");

      xn = mod_mesh_w;
      yn = mod_mesh_h;
      d3 = retutil_hex_to_3d_new(hv,cn,tn,mm_mosaic_cx,mm_mosaic_cy,xn,yn);

    }else{
      d3 = retutil_hex_to_3d(hv,cn,tn,mm_mosaic_cx,mm_mosaic_cy,&xn,&yn);
    }
    write_3d_data("zzz.hexh.3d",d3,xn,yn,tn,4,2,1);


    //
    //  Plot a horizontal trace going through cone 0
    //
    {
      int subn;
      float *subdata;
      char tname[SLEN];

      j = yn/2;

      xcut = get_farray(xn);
      for(i=0;i<xn;i++)
	xcut[i] = d3[i][j][tn-1];
      subsample_farray(xcut,xn,&subdata,&subn,2,0);
      sprintf(tname,"y=%d_tn-1",j);
      append_farray_plot("zz.hex.cut.pl",tname,subdata,subn,1);

      append_farray_plot("zz.hex.cut2.pl","hex_t=tn-1",xcut,xn,1);

      /*
      for(i=0;i<xn;i++)
	xcut[i] = d3[i][j][tn-2];
      subsample_farray(xcut,xn,&subdata,&subn,2,0);
      append_farray_plot("zz.h0.xcut.pl","sub_hex_t=tn-2",subdata,subn,1);
      //append_farray_plot("zz.h0.xcut.pl","hex_t=tn-2",xcut,xn,1);

      for(i=0;i<xn;i++)
	xcut[i] = d3[i][j][tn-3];
      subsample_farray(xcut,xn,&subdata,&subn,2,0);
      append_farray_plot("zz.h0.xcut.pl","sub_hex_t=tn-3",subdata,subn,1);
      //append_farray_plot("zz.h0.xcut.pl","hex_t=tn-3",xcut,xn,1);
      */

      printf("  Exiting after writing dump files\n");
      exit(0);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MODEL_MESH_RESP_0                             */
/*                                                                           */
/*  First steps for processing the stimulus to get the responses.            */
/*    0.  stimulus duration setting                                          */
/*    1.  photoreceptor noise and filtering                                  */
/*    2.  optional state resetting                                           */
/*    3.  horizontal cell mesh                                               */
/*                                                                           */
/*****************************************************************************/
void model_mesh_resp_0(m,s,r)
     struct model_struct *m;        // Model parameter pair list
     struct stim_struct *s;         // Stimulus data and param pair list
     struct response_struct *r;     // Response data
{
  int i,j;
  int tn,h,w;

  // 0. Set effective 'tn' for this trial, checking for variable duration
  mod_mesh_tn = mod_util_set_var_tn(mylogf,s,mod_mesh_tmax,mod_mesh_tscale);
  tn = mod_mesh_tn;
  w = mod_mesh_w;
  h = mod_mesh_h;


  // 1. Photoreceptor noise and filtering, 'mod_mesh_stim' gets filled
  if (mm_mosaic_flag == 0){
    model_mesh_photo_filter(s->d,m->mseed[r->tsi]); // Sep08 s->d is UNCHANGED
  }else if (mm_mosaic_flag == 1){
    model_mesh_mosaic_get_conex4(s,m->mseed[r->tsi]);
  }else
    exit_error("MODEL_MESH_RESP_0","Unknown mosaic_flag value");


  // 2. Initialize start values of H-mesh for first or all stimuli
  if ((r->tsi == 0) || (mod_mesh_reset == 1)){
    if (mod_mesh_type == 0){
      for(i=0;i<h;i++){
	for(j=0;j<w;j++){
	  mod_mesh_uo[i][j] = mod_mesh_stim[0][i][j];  // OR COULD USE 0.5?
	}
      }
    }else if (mod_mesh_type == 1){
      for(i=0;i<mm_mosaic_cn;i++){
	mod_mesh_hirr[i][0] = mm_mosaic_cex[i][0];
      }
      if (mod_mesh_h2_flag == 1){
	for(i=0;i<mm_mosaic_cn;i++)
	  mod_mesh_h2rr[i][0] = mm_mosaic_cex[i][0];
      }
    }else{
      mylogx(mylogf,"MODEL_MESH_RESP_0","Unknown mesh type");
    }
  }


  // 3. Run the horizontal cell mesh
  if (mod_mesh_type == 1){
    mod_mesh_hcw = mod_mesh_h1cw;
    mod_mesh_run_irreg(mod_mesh_hirr,mod_mesh_gp,mod_mesh_gh,
		       mod_mesh_dt,mod_mesh_nstep);

    if (mod_mesh_h2_flag == 1){
      mod_mesh_hcw = mod_mesh_h2cw;
      mod_mesh_run_irreg(mod_mesh_h2rr,mod_mesh_h2_gp,mod_mesh_h2_gh,
		       mod_mesh_h2_dt,mod_mesh_h2_nstep);
    }
  }else if (mod_mesh_type == 0){
    mod_mesh_run();       // Photo conductance constant
  }


  // 4. Optional resetting of 'ad1' between trials
  //    Note, ad1 depends on P and H, thus computed after H mesh.
  if ((r->tsi == 0) || (mod_mesh_reset == 1)){
    mod_mesh_ad1_reset();  // Reset the slow adaptation reference value
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_MESH_RESP_HP_GAIN                           */
/*                                                                           */
/*  Apply a gain-controlled high-pass filter to the data 'd'.  The filter    */
/*  is applied by subtracting a scaled, smoothed version of 'd' from 'd'     */
/*  itself.  The scaling is controlled by the 'gain' signal, after the       */
/*  'gain' signal has been squashed by a sigmoid.                            */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_resp_hp_gain(m,s,ci,cj,d,gain,dumpflag)
     struct model_struct *m;   // Model parameters
     struct stim_struct *s;    // Stimulus parameters
     int ci,cj;                // Coordinates for this response
     float *d;                 // response data [tn]
     float *gain;              // gain signal [tn]
     int dumpflag;             // 1-dump plots
{
  int i;
  int tn;
  float *dhp,*cf,*lp;
  char name[SLEN];

  tn = mod_mesh_tn;

  // Get low-pass of d, for gain-control
  if (mod_mesh_gain_hpflag == 1){

    //mylog(mylogf,"  MOD_MESH_RESP_HP_GAIN\n");

    // Get a smoothed version of 'd' 
    // This convolution assumes d[0] applies for d[i<0]
    lp = convolve_with_mask_causal_x0(d,tn,mod_mesh_lp_mask,mod_mesh_lpmn);

    dhp = get_zero_farray(tn); // Hold the final high-passed signal
    cf = (float *)myalloc(tn*sizeof(float)); // Temporary for plotting
    for(i=0;i<tn;i++){

      // Use a fixed sigmoid to compress the gain signal in [0..1]
      cf[i] = 0.5 + 0.5*func_sigmoid_tanh(gain[i],mod_mesh_gsig_slope,
					  mod_mesh_gsig_mean);

      // Subtract the gain-scaled smoothed signal from the original
      dhp[i] = d[i] - cf[i] * lp[i];
    }

    if (dumpflag){
      sprintf(name,"%s_LP(t)",mod_mesh_dumps);
      append_farray_plot(mod_mesh_dumpf,name,lp,tn,1);

      sprintf(name,"%s_PHADG(t)",mod_mesh_dumps);
      append_farray_plot(mod_mesh_dumpf,name,dhp,tn,1);

      sprintf(name,"%s_Gain1S(t)",mod_mesh_dumps);
      append_farray_plot(mod_mesh_dumpf,name,cf,tn,1);
    }
    myfree(lp); myfree(cf);
    
    for(i=0;i<tn;i++)
      d[i] = dhp[i];
    
    myfree(dhp);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_MESH_GAIN_GLOBAL                           */
/*                                                                           */
/*  Compute the average of the absolute value and the average of the         */
/*  squared value of the adapted PH difference signal.                       */
/*                                                                           */
/*  - Currently, the average is taken across the ENTIRE grid.                */
/*  - The gain state is advanced here and should NOT BE CHANGED ELSEWHERE.   */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_gain_global(m,rsig,rsig2)
     struct model_struct *m;   // Model parameters
     float **rsig;    // Avg of all absolute values (stim minus HC)
     float **rsig2;   // Avg of all squared values (stim minus HC)
{
  int i,j,k;
  int tn,h,w,cn,ti,tj;
  float *sig,*sig2,diff,tdiff,***u,whc,fc,avg,*lp;
  float ad1c,tad,aa,w0,w1,w2,w3,ht;
  char fname[SLEN];

  mylog(mylogf,"  MOD_MESH_GAIN_GLOBAL  Computing global power\n");

  /***
  tn   = mod_mesh_tn;
  h    = mod_mesh_h;
  w    = mod_mesh_w;
  whc  = mod_mesh_whc;
  ad1c = mod_mesh_ad1_c;
  aa   = mod_mesh_ad1_a;
  u    = mod_mesh_u;

  sig  = get_zero_farray(tn);
  sig2 = get_zero_farray(tn);
  sig[0] = sig2[0] = 0.0;       // WYETH - this is redundant??


  if (mm_mosaic_flag == 1){

    if (mod_mesh_type == 0){
      for(j=0;j<mm_mosaic_cn;j++){  // For each photoreceptor
      
	// Get cartesian index and weights for four corners on grid
	model_mesh_mosaic_get_ijw_for_cid(j,&ti,&tj,&w0,&w1,&w2,&w3);
      
	tad = mod_mesh_ad1m[j];  // For speed
      
	for(i=0;i<tn;i++){  // For all time
	
	  // Get Horiz signal at Cone j's location
	  ht = (w0*u[i][ti  ][tj  ] + w1*u[i][ti+1][tj  ] +
		w2*u[i][ti  ][tj+1] + w3*u[i][ti+1][tj+1]);
	  tdiff = mm_mosaic_cex[j][i] - whc*ht;      // Raw diff
	  diff = tdiff - tad;                        // Subtract adaptive offset
	  tad = ad1c * (tad - aa*tdiff) + aa*tdiff;  // Adapt offset to diff
	
	  if (diff > 0.0)
	    sig[i] += diff;
	  else
	    sig[i] -= diff;
	  sig2[i] += diff*diff;
	}
	mod_mesh_ad1m[j] = tad; // Restore
      }
    }else{
      exit_error("MOD_MESH_GAIN_GLOBAL","WYETH not yet imp'd for hex mesh");
    }
    fc = (float)mm_mosaic_cn;
  }else{
    for(j=0;j<h;j++)
      for(k=0;k<w;k++){  // For each photoreceptor
	
	tad = mod_mesh_ad1[j][k]; // For speed
	for(i=0;i<tn;i++){
	  tdiff = mod_mesh_stim[i][j][k] - whc*u[i][j][k]; // Raw diff
	  diff = tdiff - tad;                      // Subtract adaptive offset
	  tad = ad1c * (tad - aa*tdiff) + aa*tdiff;  // Adapt offset to diff
	  
	  if (diff > 0.0)
	    sig[i] += diff;
	  else
	    sig[i] -= diff;
	  sig2[i] += diff*diff;
	}
	mod_mesh_ad1[j][k] = tad; // Restore
      }
    fc = (float)(h*w);
  }

  for(i=0;i<tn;i++){  // Divide by the number of signals added
    sig[i] /= fc;
    sig2[i] /= fc;
  }

  //  Apply an Alpha-function smoothing to gain2
  if (mod_mesh_lp2tau > 0.0){
    lp = convolve_with_mask_causal_x0(sig2,tn,
				      mod_mesh_lp2_mask,mod_mesh_lp2mn);
    myfree(sig2);
    sig2 = lp;
  }

  if (mod_mesh_dump){

    // Print out time-averaged values of the signals
    avg = mean_farray(sig,tn);
    sprintf(ggstr,"    Global gain time average abs.val.  %f\n",avg);
    mylog(mylogf,ggstr);

    avg = mean_farray(sig2,tn);
    sprintf(ggstr,"    Global gain time average power     %f\n",avg);
    mylog(mylogf,ggstr);

    // Plot gain signals
    sprintf(fname,"%s.mesh.pl",mod_mesh_out_prefix);
    append_farray_plot(fname,"Gain1(t)",sig,tn,1);

    sprintf(fname,"%s.mesh.pl",mod_mesh_out_prefix);
    append_farray_plot(fname,"Gain2(t)",sig2,tn,1);
  }

  ***/

  *rsig = sig; *rsig2 = sig2;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_MESH_COMPUTE_DIFF_ALL                        */
/*                                                                           */
/*  1.  Compute responses for all spatial elements, 'mod_mesh_diff(m)'       */
/*  2.  Advance the adaptive state, 'mod_mesh_ad1'                           */
/*  3.  Apply gain-controlled hp-filter                                      */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_compute_diff_all(m,s,r)
     struct model_struct *m;      // Model parameter pair list
     struct stim_struct *s;       // Stimulus data and param pair list
     struct response_struct *r;   // Response data
{
  int i,j,ti,tj;
  int tn,xn,yn,save_flag,out_id;
  float ht,*diff,***u,whc,ad1c,dd,tad,fc,*lp,aa,w0,w1,w2,w3;
  float *save_ht,*save_tad,*save_d,*save_p;
  float *photo,*horiz,*h2,ws_h1,ws_h2,wlm_h1,wlm_h2;
  char tstr[SLEN];

  mylog(mylogf,"  MOD_MESH_COMPUTE_DIFF_ALL\n");

  //  1.  Set stim duration, 'mod_mesh_tn'
  //  2.  Photoreceptor noise and filtering, 'mod_mesh_stim' is set
  //  3.  Horizontal cell mesh, 'mod_mesh_u' OR 'mod_mesh_hirr' is set
  //                                            and possibly 'mod_mesh_h2rr'
  model_mesh_resp_0(m,s,r);

  tn   = mod_mesh_tn;   // Do this after 'model_mesh_resp_0'

  whc  = mod_mesh_whc;         // This is used if no H2 mesh,
  ws_h1 = mod_mesh_w_h1_s;     //   otherwise, these 4 are used.
  ws_h2 = mod_mesh_w_h2_s;
  wlm_h1 = mod_mesh_w_h1_lm;
  wlm_h2 = mod_mesh_w_h2_lm;

  xn   = mod_mesh_w;
  yn   = mod_mesh_h;
  u    = mod_mesh_u;
  ad1c = mod_mesh_ad1_c;
  aa   = mod_mesh_ad1_a;

  //printf("________________WYETH mod_mesh_util.c  --  ad1c = %f\n",ad1c);
  //printf("________________WYETH mod_mesh_util.c  --  aa   = %f\n",aa);
  // exit(0);

  if (mod_mesh_dump == 0){
    out_id = -1;                    // Don't dump output
  }else{
    out_id = mm_mosaic_out_id;      // Dump for this cone position
    save_ht  = get_farray(tn);
    save_tad = get_farray(tn);

    if (mm_mosaic_flag == 0)
      mylog_exit(mylogf,
		 "MOD_MESH_COMPUTE_DIFF_ALL  Dump not imp'd for rect. grid.");
  }

  for(i=0;i<tn;i++){ // Zero the gain signal
    mod_mesh_gain1[i] = 0.0;
    mod_mesh_gain2[i] = 0.0;
  }

  if (mm_mosaic_flag == 1){

    //
    //  MOSAIC
    //
    if (mod_mesh_diffm == NULL) // Storage for global gain
      mod_mesh_diffm = get_2d_farray(mm_mosaic_cn,tn);

    if (mod_mesh_type == 0){  // Cartesian H-mesh
      //
      //  This is the original way, w/ H-mesh matched to the stimulus
      //
      for(j=0;j<mm_mosaic_cn;j++){
	if (j == out_id)  // Are we saving traces for this cone?
	  save_flag = 1;
	else
	  save_flag = 0;

	diff = mod_mesh_diffm[j]; // Pointer to storage
	model_mesh_mosaic_get_ijw_for_cid(j,&ti,&tj,&w0,&w1,&w2,&w3);
	tad = mod_mesh_ad1m[j];  // pointer to adaptation signal 

	for(i=0;i<tn;i++){
	  ht = (w0*u[i][ti  ][tj  ] + w1*u[i][ti+1][tj  ] + // Horiz signal
		w2*u[i][ti  ][tj+1] + w3*u[i][ti+1][tj+1]);
	  dd = mm_mosaic_cex[j][i] - whc*ht;     // Raw diff
	  diff[i] = dd - tad;                    // Subtr adaptive offset
	  tad = ad1c * (tad - aa*dd) + aa*dd;    // Adapt tad to aa*dd
	  
	  if (save_flag){
	    save_ht[i]  = ht;
	    save_tad[i] = tad;
	  }
	
	  if (diff[i] > 0.0) // GAIN computation
	    mod_mesh_gain1[i] += diff[i];
	  else
	    mod_mesh_gain1[i] -= diff[i];
	  mod_mesh_gain2[i] += diff[i]*diff[i];
	}
	mod_mesh_ad1m[j] = tad;  // Restore
      }
    }else if (mod_mesh_type == 1){  // Hex/Irreg H-mesh

      //
      //  Hexagonal / Irregular H-mesh
      //
      for(j=0;j<mm_mosaic_cn;j++){
	if (j == out_id)  // Are we saving traces for this cone?
	  save_flag = 1;
	else
	  save_flag = 0;

	diff = mod_mesh_diffm[j];   // Pointer to storage
	photo = mm_mosaic_cex[j];   // Pointer to cone signal
	horiz = mod_mesh_hirr[j];   // Pointer to H-signal
	if (mod_mesh_h2_flag == 1)
	  h2 = mod_mesh_h2rr[j];   // Pointer to H2-signal
	tad   = mod_mesh_ad1m[j];   // Adaptation initial value



	for(i=0;i<tn;i++){

	  if (mod_mesh_h2_flag == 0){
	    // ORIGINAL LINE for H1 alone
	    dd = photo[i] - whc*horiz[i];        // Raw diff
	  }else{
	    //
	    //  NEW for H1 and H2 meshes
	    //
	    if (mm_mosaic_cid[j] == 0){ // This is an s-cone bipolar
	      dd = photo[i] -  ws_h1*horiz[i] -  ws_h2*h2[i];    // Raw diff
	    }else{
	      dd = photo[i] - wlm_h1*horiz[i] - wlm_h2*h2[i];    // Raw diff
	    }
	  }


	  diff[i] = dd - tad;                  // Subtr adaptive offset
	  tad = ad1c * (tad - aa*dd) + aa*dd;  // Adapt tad to aa*dd
	  // Previous line is essentially:  ad1c*OLD + (1-ad1c)*NEW
	  
	  if (save_flag){
	    //save_ht[i]  = ht;  // *** NOT NEEDED, ALREADY SAVED
	    save_tad[i] = tad;
	  }

	  //
	  //  Update two gain signals
	  //
	  if (diff[i] > 0.0)
	    mod_mesh_gain1[i] += diff[i];
	  else
	    mod_mesh_gain1[i] -= diff[i];
	  mod_mesh_gain2[i] += diff[i]*diff[i];
	}
	mod_mesh_ad1m[j] = tad;  // Save updated adaptive offset


	// ****** WYETH DEBUG ************
	// ****** WYETH DEBUG ************
	/***
	if (j == 0){
	  append_farray_plot("zzz.pl","P",mm_mosaic_cex[j],tn,1);
	  append_farray_plot("zzz.pl","H",mod_mesh_hirr[j],tn,1);
	  append_farray_plot("zzz.pl","Diff",mod_mesh_diffm[j],tn,1);
	  exit(0);
	}
	***/
	// ****** WYETH DEBUG ************
	// ****** WYETH DEBUG ************

      }
    }


    fc = (float)mm_mosaic_cn;

  }else{
    //
    //  CARTESIAN GRID
    //
    if (mod_mesh_diff == NULL) // Storage for global gain
      mod_mesh_diff = get_3d_farray(xn,yn,tn);
    
    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){

	diff = mod_mesh_diff[i][j]; // Pointer to storage
	
	//  - Get 'diff' response for (i,j)th unit
	//  - This should match the functionality in 'mod_mesh_resp_ph_diff'

	tad = mod_mesh_ad1[i][j]; // For speed, and don't change original
	for(ti=0;ti<tn;ti++){
	  ht = u[ti][i][j];
	  dd = mod_mesh_stim[ti][i][j] - whc*ht; // Raw difference
	  diff[ti] = dd - tad;                   // Subtract adaptive offset
	  tad = ad1c * (tad - aa*dd) + aa*dd;    // Adapt offset to diff 'dd'
	  
	  if (diff[ti] > 0.0) // GAIN computation
	    mod_mesh_gain1[ti] += diff[ti];
	  else
	    mod_mesh_gain1[ti] -= diff[ti];
	  mod_mesh_gain2[ti] += diff[ti]*diff[ti];
	}
	mod_mesh_ad1[i][j] = tad; // Restore
      }
    }
    fc = (float)(xn*yn); // To compute average gain signal
  }


  // DUMP OUTPUT
  if (out_id >= 0){

    save_p = mm_mosaic_cex[out_id];     // Photoreceptor signal
    save_d = mod_mesh_diffm[out_id];    // P-H diff

    sprintf(tstr,"%s_P",mod_mesh_dumps);
    append_farray_plot(mod_mesh_dumpf,tstr,save_p,tn,1);

    sprintf(tstr,"%s_H",mod_mesh_dumps);
    if (mod_mesh_type == 1){
      append_farray_plot(mod_mesh_dumpf,tstr,mod_mesh_hirr[out_id],tn,1);

      if (mod_mesh_h2_flag == 1){
	// WYETH - IS THIS USEFUL ??
	sprintf(tstr,"%s_H2",mod_mesh_dumps);
	append_farray_plot(mod_mesh_dumpf,tstr,mod_mesh_h2rr[out_id],tn,1);
      }

    }else
      append_farray_plot(mod_mesh_dumpf,tstr,save_ht,tn,1);

    sprintf(tstr,"%s_PHAD",mod_mesh_dumps);
    append_farray_plot(mod_mesh_dumpf,tstr,save_d,tn,1);

    sprintf(tstr,"%s_ad1",mod_mesh_dumps);
    append_farray_plot(mod_mesh_dumpf,tstr,save_tad,tn,1);

    myfree(save_ht);
    myfree(save_tad);

    // Dump 'diff' for entire mosaic
    // ,0,1 ==> no transpose, printflag
    write_2d_data("zzz.phad.2d",mod_mesh_diffm,0,0,mm_mosaic_cn,tn,4,2,0,1);

  }

  // Dump a slice of the horizontal mesh:  u[ti][xi][yi]
  if (mod_mesh_slice_ti != -1){
    if (mod_mesh_type == 0){
      sprintf(tstr,"H_t%d_x%d",mod_mesh_slice_ti,mod_mesh_slice_xi);
      append_farray_plot("zzz.hsect.pl",tstr,
			 u[mod_mesh_slice_ti][mod_mesh_slice_xi],yn,1);
    }else{
      mylogx(mylogf,"MOD_MESH_COMPUTE_DIFF_ALL","H-slice not imp'd for irr");
    }
  }

  for(i=0;i<tn;i++){
    mod_mesh_gain1[i] /= fc;  // Divide by number of signals
    mod_mesh_gain2[i] /= fc;
  }

  //  Apply an Alpha-function smoothing to gain2
  if (mod_mesh_lp2tau > 0.0){
    printf("********************* YWETHEHT -----\n");
    if (mod_mesh_gain2 == NULL)
      printf(". . . is NULL\n");
    else
      printf(". . . is NULL\n");
    if (mod_mesh_lp2_mask == NULL)
      printf("lp2 MASK. . . is NULL\n");
    else
      printf("lp2 MASK. . . is nnooooootttt  NULL\n");
    
    lp = convolve_with_mask_causal_x0(mod_mesh_gain2,tn,
				      mod_mesh_lp2_mask,mod_mesh_lp2mn);
    myfree(mod_mesh_gain2);
    mod_mesh_gain2 = lp;
  }

  // Apply gain-controlled high-pass filter to the 'diff' signal
  if (mod_mesh_gain_hpflag == 1){
    if (mm_mosaic_flag == 1){
      printf("*-----------------****************  HERE WYETH GAIN\n");
      for(i=0;i<mm_mosaic_cn;i++){
	//  WYETH - consider making a fast version of this function
	//          that uses a lookup table for the sigmoid, etc.
	if (i == out_id)
	  mod_mesh_resp_hp_gain(m,s,i,-1,mod_mesh_diffm[i],mod_mesh_gain1,1);
	else
	  mod_mesh_resp_hp_gain(m,s,i,-1,mod_mesh_diffm[i],mod_mesh_gain1,0);
      }
    }else{
      if (mod_mesh_gain_hpflag == 1){
	for(i=0;i<xn;i++){
	  for(j=0;j<yn;j++){
	    //  WYETH - consider making a fast version of this function
	    //          that uses a lookup table for the sigmoid, etc.
	    mod_mesh_resp_hp_gain(m,s,i,j,mod_mesh_diff[i][j],mod_mesh_gain1,
				  0);
	  }
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_MESH_GET_SPIKES_2D_FLAG                       */
/*                                                                           */
/*  Get spike trains for any spatial element where flag = 1.                 */
/*                                                                           */
/*  NOTES                                                                    */
/*  - Before calling this, must do 'mod_mesh_compute_diff_all'               */
/*  - Spikes are returned in 'msec'                                          */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_get_spikes_2d_flag(m,s,r,pop_name,onflag,xn,yn,flag,
				 rs2,rcnt2)
     struct model_struct *m;      // Model parameter pair list
     struct stim_struct *s;       // Stimulus data and param pair list
     struct response_struct *r;   // Response data
     char *pop_name;              // Population name, e.g., "rgc" or "rgc_p"
     int onflag;                  // 1=ON, 0=OFF
     int xn,yn;                   // size of grid
     int **flag;                  // 1-make spikes for this unit [xn][yn]
     float ****rs2;  // [xn][yn][rcnt2] return spikes at flagged coords
     int ***rcnt2;   // [xn][yn] number of spikes at flagged coords
{
  int i,j,ti;
  int tn,**cnt2,unit_count,spike_count,vn,seed,*seedlist;
  float ***s2,*gti,*gtx,*v,samp,*td;
  char *sgen;
  struct onode *po,*spko;

  mylog(mylogf,"  MOD_MESH_GET_SPIKES_2D_FLAG\n");
  sprintf(ggstr,"      Layer %s\n",pop_name);
  mylog(mylogf,ggstr);


  //printf("    Pop %s  xn,yn  %d %d\n",pop_name,xn,yn);

  if (mm_mosaic_flag == 1)
    mylog_exit(mylogf,"MOD_MESH_GET_SPIKES_2D_FLAG  Mosaic not imp'd.\n");

  if ((xn != mod_mesh_w) || (yn != mod_mesh_h))
    mylog_exit(mylogf,"MOD_MESH_GET_SPIKES_2D_FLAG  Dimension mismatch.\n");

  // Convenient onode pointers
  po = onode_get_node_type_item_val(m->o,"pop","name",pop_name);
  if (po == NULL)
    mylog_exit(mylogf,"MOD_MESH_GET_SPIKES_2D_FLAG  Bad pop_name.\n");
  spko = onode_child_get_unique(po,"spike_gen");
  sgen = onode_getpar_chr_exit(spko,"type");


  if (myid == -1){ // WYETH - REMOVE; check if dumping output for an LGN unit
    if (onode_test_ostr(po,"ifc_dump_coord")){
      printf("***** WYETH - 'ifc_dump_coord' is to be phased out - remove?\n");
    }
  }

  // Storage for spike trains
  cnt2 = get_zero_2d_iarray(xn,yn);
  s2 = (float ***)myalloc(xn*sizeof(float **));
  for(i=0;i<xn;i++)
    s2[i] = (float **)myalloc(yn*sizeof(float *));

  tn   = mod_mesh_tn;
  samp = 1.0/mod_mesh_tscale; // Sampling for 'diff'

  gtx = get_zero_farray(tn);
  gti = get_zero_farray(tn);  // Set inhibitory input to zero
  seedlist = get_seeds(m->mseed[r->tsi],1000000,xn*yn);


  unit_count = 0;   // Units for which spike trains are generated
  spike_count = 0;  // Total of all spike counts
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      if (flag[i][j]){

	td = mod_mesh_diff[i][j];  //   'diff' signal for unit (i,j)

	if (onflag == 0)
	  for(ti=0;ti<tn;ti++)
	    gtx[ti] = -td[ti];  // For OFF response, negative 'diff'
	else
	  for(ti=0;ti<tn;ti++)
	    gtx[ti] = td[ti];

	if (mod_mesh_gain_g2in){  // Inhibotory gain control (Not for Poisson)
	  for(ti=0;ti<tn;ti++)
	    gti[ti] = mod_mesh_gain2[ti];
	}else{
	  for(ti=0;ti<tn;ti++)
	    gti[ti] = 0.0; // Zero in input; 'gti' can change in 'ifc_test'
	}

	seed = seedlist[i*yn + j];

	//  Call either 'ifc_test' or 'ifc_util_poisson'
	mod_mesh_spike_gen(m,gtx,gti,tn,seed,samp,sgen,pop_name,NULL,0,
			   &(s2[i][j]),&(cnt2[i][j]),&v,&vn);

	unit_count += 1;
	spike_count += cnt2[i][j];
      }else{ // Flag was not set
	s2[i][j] = NULL;  // Important for freeing s2 later
	cnt2[i][j] = -1;  // -1 indicates spikes never computed
      }
    }
  }
  myfree(gtx);
  myfree(gti);
  myfree(sgen);
  myfree(seedlist);

  sprintf(ggstr,"      Spikes generated at %d of %d units, %.4f spk/cell\n",
	  unit_count,xn*yn,(float)spike_count/(float)unit_count);
  mylog(mylogf,ggstr);

  *rs2 = s2; *rcnt2 = cnt2;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_MESH_GET_SPIKES_3D_FLAG                       */
/*                                                                           */
/*  Get spike trains for any spatial element where flag = 1.                 */
/*                                                                           */
/*  NOTES                                                                    */
/*  - Before calling this, must do 'mod_mesh_compute_diff_all'               */
/*  - Spikes are returned in msec                                            */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_get_spikes_3d_flag(m,s,r,pop_name,cmap,xn,yn,zn,rs,rcnt)
     struct model_struct *m;      // Model parameter pair list
     struct stim_struct *s;       // Stimulus data and param pair list
     struct response_struct *r;   // Response data
     char *pop_name;              // Population name, e.g., "rgc" or "rgc_p"
     int ***cmap;                 // 1-make spikes for this unit [xn][yn][zn]
     int xn,yn,zn;                // size of grid
     float *****rs;     // [zn][xn][yn][rcnt2] return spikes at flagged coords
     int ****rcnt;      // [zn][xn][yn] number of spikes at flagged coords
{
  int i,j,k,ti;
  int tn,***cnt,unit_count,spike_count,vn,seed,*seedlist,dumpi,dumpk;
  float ****s1,*gti,*gtx,*v,samp,*td;
  char *sgen,*dump_ptr,dump_name[SLEN];
  struct onode *po,*spko;

  mylog(mylogf,"  MOD_MESH_GET_SPIKES_3D_FLAG\n");
  sprintf(ggstr,"      Layer %s\n",pop_name);
  mylog(mylogf,ggstr);

  //printf("    Pop %s  xn,yn  %d %d\n",pop_name,xn,yn);

  if (mm_mosaic_flag != 1)
    mylogx(mylogf,"MOD_MESH_GET_SPIKES_3D_FLAG","Intended for mosaic.\n");

  if ((xn != mm_mosaic_cn) || (yn != 1) || (zn != 2))
    mylogx(mylogf,"MOD_MESH_GET_SPIKES_3D_FLAG","Bad cmap size.\n");


  // Convenient onode pointers
  po = onode_get_node_type_item_val(m->o,"pop","name",pop_name);
  if (po == NULL)
    mylog_exit(mylogf,"MOD_MESH_GET_SPIKES_3D_FLAG  Bad pop_name.\n");
  spko = onode_child_get_unique(po,"spike_gen");
  sgen = onode_getpar_chr_exit(spko,"type");

  dumpi = dumpk = -1;
  if (myid == -1){ // WYETH - We're still using this in sc1.moo
    if (onode_test_ostr(po,"spike_dump_xi")){
      dumpi = onode_getpar_int_exit(po,"spike_dump_xi");
      dumpk = onode_getpar_int_exit(po,"spike_dump_zi");
      sprintf(dump_name,"%s_%d_0_%d.%s",pop_name,dumpi,dumpk,s->name);
      //printf("***** WYETH 'ifc_dump_coord' is to be phased out - remove?\n");
    }
  }

  // Storage for spike trains
  cnt = get_zero_3d_iarray(zn,xn,yn);
  s1 = (float ****)myalloc(zn*sizeof(float ***));
  for(i=0;i<zn;i++){
    s1[i] = (float ***)myalloc(xn*sizeof(float **));
    for(j=0;j<xn;j++)
      s1[i][j] = (float **)myalloc(yn*sizeof(float *));
  }

  tn   = mod_mesh_tn;
  samp = 1.0/mod_mesh_tscale; // Sampling for 'diff'

  gtx = get_zero_farray(tn);
  gti = get_zero_farray(tn);  // Set inhibitory input to zero
  seedlist = get_seeds(m->mseed[r->tsi],1000000,xn*zn);

  unit_count = 0;   // Units for which spike trains are generated
  spike_count = 0;  // Total of all spike counts
  for(k=0;k<zn;k++){
    for(i=0;i<xn;i++){
      if (cmap[i][0][k]){  // y = 0

	//
	//  *** WYETH - check this
	//
	if (mod_mesh_csig > 0.0){
	  td = mod_mesh_spatial_sum(i,-1,-1,0);  // Spatial summation, 0-quiet
	}else{
	  td = mod_mesh_diffm[i];  //   'diff' signal for cone 'i'
	}

	if (k == 0)
	  for(ti=0;ti<tn;ti++)
	    gtx[ti] = -td[ti];  // For OFF response, negative 'diff'
	else
	  for(ti=0;ti<tn;ti++)
	    gtx[ti] = td[ti];

	// WYETH added 2012 Apr to reduce memory leak?
	if (mod_mesh_csig > 0.0){
	  myfree(td);
	}

	if (mod_mesh_gain_g2in){  // Inhib gain control (Not for Poisson)
	  for(ti=0;ti<tn;ti++)
	    gti[ti] = mod_mesh_gain2[ti];
	}else{
	  for(ti=0;ti<tn;ti++)
	    gti[ti] = 0.0; // Zero in input; 'gti' can change in 'ifc_test'
	}

	seed = seedlist[k*xn + i];

	//printf("%s  %4d %4d seed = %d\n",pop_name,i,k,seed);

	if ((i == dumpi) && (k == dumpk)){
	  dump_ptr = dump_name;
	  //printf("*** dump_ptr = %s\n",dump_ptr);
	}else
	  dump_ptr = NULL;

	//  Call either 'ifc_test' or 'ifc_util_poisson'
	mod_mesh_spike_gen(m,gtx,gti,tn,seed,samp,sgen,pop_name,dump_ptr,0,
			   &(s1[k][i][0]),&(cnt[k][i][0]),&v,&vn);

	unit_count += 1;
	spike_count += cnt[k][i][0];
      }else{ // Flag was not set
	s1[k][i][0] = NULL;  // Important for freeing s2 later
	cnt[k][i][0] = -1;  // -1 indicates spikes never computed
      }
    }
  }
  myfree(gtx);
  myfree(gti);
  myfree(sgen);
  myfree(seedlist);

  sprintf(ggstr,"      Spikes generated at %d of %d units, %.4f spk/cell\n",
	  unit_count,xn*zn,(float)spike_count/(float)unit_count);
  mylog(mylogf,ggstr);

  *rs = s1;
  *rcnt = cnt;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_MESH_RESP_PH_DIFF                          */
/*                                                                           */
/*  Return the difference between the photoreceptor signal 'mod_mesh_stim'   */
/*  and the horizontal cell signal 'mod_mesh_u'.                             */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_resp_ph_diff(m,rht,rdiff,name,xi,yi)
     struct model_struct *m;    // Model parameter pair list
     float **rht;               // Horiz signal
     float **rdiff;             // Photo - Horiz
     char *name;
     int xi,yi;
{
  int i,j,k;
  int tn,cid,ti,tj,tti,ttj,h,w,cn;
  float *ht,*diff,***u,whc,csig,weight,totw,*p,ad1c,*att,dd,tad;
  float w0,w1,w2,w3,tw0,tw1,tw2,tw3,aa;
  char fname[SLEN],tstr[SLEN];

  // *********** WYETH --- this is no longer used **********/
  // *********** WYETH --- this is no longer used **********/
  // *********** WYETH --- this is no longer used **********/
  // *********** WYETH --- this is no longer used **********/
  // Why are we keeping it (2011 Oct)?

  //  2008 Sep. WYETH - it doesn't make sense that the returned arrays 
  //    start with 0.0 as their first value.  The stim data starts at i=1,
  //    but these responses should probably start at t=0.

  /***

  sprintf(fname,"%s.mesh.pl",mod_mesh_out_prefix);

  whc  = mod_mesh_whc;
  csig = mod_mesh_csig;

  h    = mod_mesh_h;      // Height
  w    = mod_mesh_w;      // Width

  tn   = mod_mesh_tn;
  u    = mod_mesh_u;      // Horiz cell signal
  ad1c = mod_mesh_ad1_c;  // exp(-1.0/(mod_mesh_ad1_tau / mod_mesh_tscale));
  aa   = mod_mesh_ad1_a;  // Amplitude of adaptation
  
  p = get_zero_farray(tn);   // Temporary storage for plotting
  att = get_zero_farray(tn); // Temporary storage for plotting

  ht   = get_zero_farray(tn);  // To be returned
  diff = get_zero_farray(tn);
  ht[0]   = 0.0;
  diff[0] = 0.0;

  if (mm_mosaic_flag == 1){
    //
    //  MOSAIC
    //

    if (mod_mesh_type == 1) // WYETH - UPDATE for hex mesh??
      exit_error("MOD_MESH_RESP_PH_DIFF","Not imp'd for hex-mesh yet");

    cn = mm_mosaic_cn;
    cid = mm_mosaic_out_id;  // Index of cone for output
    model_mesh_mosaic_get_ijw_for_cid(cid,&ti,&tj,&w0,&w1,&w2,&w3);
    if (csig <= 0.0){

      tad = mod_mesh_ad1m[cid]; // pointer to adaptation signal
      for(i=0;i<tn;i++){  // WYETH changed from i=1... Sept 2008
	// Weighted H-cell signal
	ht[i] = (w0*u[i][ti  ][tj  ] + w1*u[i][ti+1][tj  ] +
		 w2*u[i][ti  ][tj+1] + w3*u[i][ti+1][tj+1]);
	dd = mm_mosaic_cex[cid][i] - whc*ht[i]; // Raw diff
	diff[i] = dd - tad; // Subtr adaptive offset
	tad = ad1c * (tad - aa*dd) + aa*dd;     // Adapt tad to aa * diff 'dd'
	att[i] = tad;                           // For plotting 'ad1'
	p[i]   = mm_mosaic_cex[cid][i];
      }
    }else{
      exit_error("MOD_MESH_RESP_PH_DIFF","Not imp'd for mosaic yet");
      totw = 0.0;
      for(j=0;j<cn;j++){
	model_mesh_mosaic_get_ijw_for_cid(j,&tti,&ttj,&tw0,&tw1,&tw2,&tw3);
	// WYETH - why not use float cone positions instead of
	//   tti,ttj,ti,tj here ??? See email to/from Jimmy Jia
	weight = func_2d_gaussian((float)tti,(float)ttj,(float)ti,
				  (float)tj,csig,csig,0);
	if (weight > mod_mesh_csig_min){
	  totw += weight;
	  for(i=0;i<tn;i++){  // WYETH changed from i=1... Sept 2008
	    ht[i] = (tw0*u[i][tti  ][ttj  ] + tw1*u[i][tti+1][ttj  ] +
		     tw2*u[i][tti  ][ttj+1] + tw3*u[i][tti+1][ttj+1]);
	    diff[i] += weight * (mm_mosaic_cex[j][i] - whc*ht[i]);
	  }
	}
      }
      sprintf(ggstr,"  Total weight = %f\n",totw);
      mylog(mylogf,ggstr);
    }
  }else{
    //
    //  CARTESIAN GRID
    //
    if (csig <= 0.0){ // Center is single element, no spatial summation

      tad = mod_mesh_ad1[xi][yi]; // For speed, and we don't change original
      for(i=0;i<tn;i++){
	ht[i] = u[i][xi][yi];
	dd = mod_mesh_stim[i][xi][yi] - whc*ht[i]; // Raw difference
	diff[i] = dd - tad;                      // Subtract adaptive offset
	tad = ad1c * (tad - aa*dd) + aa*dd;      // Adapt tad to aa * diff 'dd'
	att[i] = tad;                            // For plotting 'ad1'
	p[i] = mod_mesh_stim[i][xi][yi];         // For plotting photo signal
      }
      // Do not restore 'tad' value,  '...gain_global' will set later
    }else{  // Sum over Gaussian for center signal
      totw = 0.0;
      for(j=0;j<h;j++)
	for(k=0;k<w;k++){
	  weight = func_2d_gaussian((float)j,(float)k,(float)xi,(float)yi,
				    csig,csig,1);
	  if (weight > mod_mesh_csig_min){
	    totw += weight;
	    tad = mod_mesh_ad1[j][k]; // For speed
	    for(i=0;i<tn;i++){
	      ht[i] = u[i][j][k];
	      //ORIG: diff[i] += weight * (mod_mesh_stim[i][j][k] - whc*ht[i]);
	      dd = mod_mesh_stim[i][j][k] - whc*ht[i]; // Raw difference
	      diff[i] += weight * (dd - tad);          // Subtract offset
	      tad = ad1c * (tad - aa*dd) + aa*dd;      // Adapt 'ad1' to 'dd'
	    }
	    // Do not restore 'tad' value, it will be done elsewhere
	  }
	}
      sprintf(ggstr,"  Total weight = %f\n",totw);
      mylog(mylogf,ggstr);
    }
  }

  if (mod_mesh_dump){
    sprintf(tstr,"%s_P",name);
    append_farray_plot(fname,tstr,p,tn,1);

    sprintf(tstr,"%s_H",name);
    append_farray_plot(fname,tstr,ht,tn,1);

    sprintf(tstr,"%s_PHAD",name);
    append_farray_plot(fname,tstr,diff,tn,1);

    sprintf(tstr,"%s_ad1",name);
    append_farray_plot(fname,tstr,att,tn,1);
  }

  ***/

  myfree(p);
  *rht = ht; *rdiff = diff;
}
/**************************************-**************************************/
/*                                                                           */
/*                      MODEL_MESH_RGC_01_GET_RESPONSE                       */
/*                                                                           */
/*  Calling sequence:                                                        */
/*    1.  model_mesh_rgc_01_prep                                             */
/*    2.  this routine                                                       */
/*        ...                                                                */
/*    N.  model_mesh_rgc_01_done                                             */
/*                                                                           */
/*****************************************************************************/
void model_mesh_rgc_01_get_response(m,s,r,k)
     struct model_struct *m;        // Model parameter pair list
     struct stim_struct *s;         // Stimulus data and param pair list
     struct response_struct *r;     // Response data
     int k;                         // Repeat number
{
  mod_mesh_dumps = s->name;  // Pointer

  mod_mesh_compute_diff_all(m,s,r);
  // WYETH, this used to call 'mod_mesh_resp_ph_diff'.  Why are we keeping
  // that routine around, if it is not being used?
  // Is it needed for the cartesian version of the model?

  // WYETH, this also used to call 'mod_mesh_gain_global'...

  mod_mesh_resp_handover(m,r);
}
/**************************************-**************************************/
/*                                                                           */
/*                       MODEL_MESH_MAG_01_GET_RESPONSE                      */
/*                                                                           */
/*  Calling sequence:                                                        */
/*    1.  model_mesh_mag_01_prep                                             */
/*    2.  this routine                                                       */
/*        ...                                                                */
/*    N.  model_mesh_rgc_01_done                                             */
/*                                                                           */
/*****************************************************************************/
void model_mesh_mag_01_get_response(m,s,r,k)
     struct model_struct *m;        // Model parameter pair list. 
     struct stim_struct *s;         // Stimulus data and param pair list 
     struct response_struct *r;     // Response data 
     int k;                         // Repeat number 
{
  int i,j,jj;
  int w,h,tn,rvn;
  float *rt,*diff,csig,whc,weight,totw;
  float *rv;

  csig = mod_mesh_csig;
  whc = mod_mesh_whc;

  w = mod_mesh_w;
  h = mod_mesh_h;
  tn = mod_mesh_tn;

  // Photoreceptor processing, 'model_mesh_stim' gets filled
  model_mesh_photo_filter(s->d,m->mseed[r->tsi]);

  // Clear the 'old' state
  for(i=0;i<h;i++)
    for(j=0;j<w;j++)
      mod_mesh_uo[i][j] = 0.0;

  mod_mesh_run();

  rt = get_zero_farray(tn);
  diff = get_zero_farray(tn);
  rt[0] = 0.0;
  diff[0] = 0.0;

  if (csig <= 0.0){
    for(i=1;i<tn;i++){
      //rt[i] = mod_mesh_u[i-1][h/2][w/2];
      rt[i] = mod_mesh_u[i][h/2][w/2]; // CHANGED Mar 25, 2009
      diff[i] = mod_mesh_stim[i][h/2][w/2] - whc*rt[i];
    }
  }else{
    totw = 0.0;
    for(j=0;j<h;j++)
      for(jj=0;jj<w;jj++){
	weight = func_2d_gaussian((float)j,(float)jj,(float)h/2.0,
				  (float)w/2.0,csig,csig,1);
	totw += weight;
	for(i=1;i<tn;i++){
	  //rt[i] = mod_mesh_u[i-1][j][jj];
	  rt[i] = mod_mesh_u[i][j][jj];  // CHANGED Mar 25, 2009
	  diff[i] += weight * (mod_mesh_stim[i][j][jj] - whc*rt[i]);
	}
      }
    sprintf(ggstr,"  Total weight = %f\n",totw);
    mylog(mylogf,ggstr);
  }

  if (myid == -1){
    append_farray_plot("zz.mesh.pl","rt",rt,tn,1);
    append_farray_plot("zz.mesh.pl","diff",diff,tn,1);
  }

  // IFC Model
  {
    int ns,nmask;
    float *gti,*gtx,samp,*fs,sampling,*mask;
    float *pow,powtau,*spow;

    // Compute the power in the diff signal, integrate in window 
    pow = copy_farray(diff,tn);
    square_farray(pow,tn);
    powtau = 20.0;             // WYETH - make this a parameter 
    nmask = powtau * 4;
    mask = one_sided_exp_farray(0.0,powtau,1.0,nmask);
    norm_area_farray(mask,nmask,1.0);
    spow = convolve_with_mask_causal(pow,tn,mask,nmask);

    if (myid == -1){
      append_farray_plot("zz.mesh.gtx.pl","pow",pow,tn,1);
      append_farray_plot("zz.mesh.gtx.pl","spow",spow,tn,1);
    }
    
    // WYETH - 'spow' ranges from 0 to 0.071

    gti = copy_farray(spow,tn);
    multiply_farray(gti,tn,1.0/0.071);
    //gti = get_zero_farray(tn);

    // This signal should be TF-tuned

    gtx = diff; // Or rt? 
    samp = 1.0/mod_mesh_tscale;

    if (myid == -1)
      append_farray_plot("zz.mesh.gtx.pl","gtx",gtx,tn,1);

    //printf("  samp = %f  dur = %d  tscale = %f\n",samp,tn,mod_mesh_tscale);

    ifc_test(myid,m,"ifc",gtx,gti,samp,tn,m->mseed[r->tsi],1,(char *)NULL,
	     &fs,&ns,0,&rv,&rvn);
    /**** To replace call above for debuggin
	  fs = get_zero_farray(3);
	  ns = 3;
	  fs[0] = 2.0;
	  fs[1] = 10.0;
	  fs[2] = 50.0;
    */
    
    myfree(gti); myfree(gtx);
    
    // Spikes come in ms, adjust based on the 'save_spikes' sampling 
    exit_error("MODEL_MESH_MAG_01_GET_RESPONSE",
	       "Old way, m->ppl, no longer supported");
    sampling = 1000.0; // WYETH 
    //sampling = paramfile_get_nth_float_param_or_exit(m->ppl,"save_spikes",1);
    multiply_farray(fs,ns,sampling/1000.0);
    
    // Store response, update counter
    if (r->ns > 0){
      //r->s[0][r->tsi] = fs;
      //r->cnt[0][r->tsi] = ns; // 'tsi' is updated in 'run_gen_tuning_curve'
      r->s[0][r->gtsi] = fs;
      r->cnt[0][r->gtsi] = ns; // 'gtsi' is updated in 'run_gen_tuning_curve'
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_MESH_RUN_MESH_RGC_01                        */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_run_mesh_rgc_01(m,s,r,k,action)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     int k;                     // Trial number, when action = 1 
     int action;                // -1-cleanup, 0-prep, 1-run 
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // DO THIS FIRST

  mylog(mylogf,"  MOD_MESH_RUN_MESH_RGC_01\n");

  if (action == 0){
    if (m->ppl == NULL)
      model_mesh_rgc_01_prep(m,s,r);
    else{ // WYETH OLD WAY - remove
      mylog(mylogf,"  *** See .../libc/archive/mod_mesh_util for old code\n");
      mylog_exit(mylogf,"MOD_MESH_RUN_MESH_RGC_01  Old .mod not supported\n");
    }
  }else if (action == 1){
    model_mesh_rgc_01_get_response(m,s,r,k);
  }else if (action == -1){
    model_mesh_rgc_01_done();
  }else if (action == 10){
    model_mesh_rgc_01_prep(m,s,r); // r,s are NULL, write MAR file
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_MESH_RUN_MESH_MAG_01                        */
/*                                                                           */
/*****************************************************************************/
void mod_mesh_run_mesh_mag_01(m,s,r,k,action)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     int k;                     // Trial number, when action = 1
     int action;                // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  mylog(mylogf,"  MOD_MESH_RUN_MESH_MAG_01\n");
  
  if (action == 0){
    model_mesh_rgc_01_prep(m,s,r);
  }else if (action == 1){
    model_mesh_mag_01_get_response(m,s,r,k);
  }else if (action == -1){
    model_mesh_rgc_01_done();
  }
}
