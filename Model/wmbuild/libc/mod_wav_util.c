/*****************************************************************************/
/*                                                                           */
/*   mod_wav_util.c                                                         */
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
#include "nr_util.h"
#include "carray_util.h"
#include "iarray_util.h"
#include "farray_util.h"
#include "func_util.h"
#include "myrand_util.h"
#include "data_util.h"
#include "ndata_util.h"
#include "paramfile_util.h"
#include "retina_util.h"
#include "stim_util.h"
#include "mod_util.h"
#include "ifc_util.h"
#include "pop_cell_util.h"
#include "pop_util.h"
#include "ifc.h"
#include "mod.h"
#include "paramfile.h"
#include "ndata.h"

// For cluster computation  (static needed w/ initialization
static int myid = -1;       // Process ID, -1-Single, 0-Master, >0-Slave 
static char *mylogf = NULL; // Logfile name 
char ggstr[LONG_SLEN];      // Global log string 

// Global for wave model

int   mod_wav_xn;          // Width of mesh
int   mod_wav_yn;          // Height of mesh
int   mod_wav_zn;          // Number of RGC populations JBCODE
int   mod_wav_tn;          // Time duration for current trial
double mod_wav_tscale;     // seconds per time unit
float mod_wav_sscale;      // degr/spatial unit

int   mod_wav_reset_time;  // 1 - reset time index after 500000 steps JBCODE
                           // for waves to stabilise. 0 - don't.  JBCODE
int   mod_wav_tt;          // Time to (re)start the simulation
int   mod_wav_seed;        // Track seed across multiple epochs
int   mod_wav_spike_full;  // Flag to indicate it is time to end trial
                           //   and return responses in indefinite run mode

float mod_wav_umpix;       // Microns per pixel
float mod_wav_am_tau;      // Decay time (s)
float mod_wav_am_fc;       // Precomputed decay constant (tscale units)
float mod_wav_am_rad;      // 
float mod_wav_am_thr;      // 
float mod_wav_gc_tau;      // 
float mod_wav_gc_fc;       // Precomputed decay constant (tscale units)
float mod_wav_gc_rad;      // Radius (um)
float mod_wav_gc_thr;      // 
int   mod_wav_gc_sgen_n;   // Pre-compute for speed
float mod_wav_gc_rate;     // GC spike rate (spk/s)
int   mod_wav_gc_seed;     // Seed for generating GC spikes
float mod_wav_gc_prob;     // Firing prob per sampling unit
// JBCODE next two lines had +1 dimension, first dim is the population index
int ****mod_wav_gc_spk;    // Spike times for GC output [pop][xn][yn][nspk]
int  ***mod_wav_gc_cnt;    // Spike counts for GC output [pop][xn][yn][nspk]
int    mod_wav_gc_smax;    // Size of spike storage
float  mod_wav_gc_samp;    // Samples per second
float mod_wav_refr_mu;     // 
float mod_wav_refr_sd;     // 
int   mod_wav_refr_seed;   // 
int **mod_wav_refr_am;     // Refractory periods for Amacrine cells
float mod_wav_spont_prob;  // Probability per time step of firing
float mod_wav_act_dur;     // Firing duration (s)
int   mod_wav_act_durn;    // Firing duration (time units)
float mod_wav_pop_delay;   // Time delay between RGC populations (s) JBCODE


float **mod_wav_x_am;      // Excitation level, Amacrine
float **mod_wav_x_gc;      // Excitation level, Ganglion
float **mod_wav_synin;     // Synaptic input
int   **mod_wav_am_state;  // 0-not firing, >0-Firing

int  **mod_wav_con_aa_n;   // [xn][yn] Number of connections
int ***mod_wav_con_aa_x;   // [xn][yn][...con_aa_n] x-index
int ***mod_wav_con_aa_y;   // [xn][yn][...con_aa_n] y-index

float mod_wav_init_blank;  // (s) initial duration to blank spike output

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
/*                              MOD_WAV_UTIL_RAN2                            */
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
float mod_wav_util_ran2(idum)
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
/*                         MOD_WAV_MAKE_CONN_LIST_RAD                        */
/*                                                                           */
/*  Given a radius 'rad' and a grid of cells on [xn,yn], make the lists,     */
/*  'rx' and 'ry' of x- and y-coordinates of units that are within 'rad'     */
/*  of each unit.  Also return 'rn', the number of units in each list.       */
/*                                                                           */
/*****************************************************************************/
void mod_wav_make_conn_list_rad(xn,yn,rad,rx,ry,rn)
     int xn,yn;           // size of grid
     float rad;           // maximum connection radius
     int ****rx,****ry;   // [xn][yn][rn[][]]  x,y-index lists
     int ***rn;           // [xn][yn] number of connections in list
{
  int i,j,ii,jj;
  int ***x,***y,**n,*tx,*ty,cn;
  int ir,r2,d2,x0,x1,y0,y1;

  mylog(mylogf,"  MOD_WAV_MAKE_CONN_LIST_RAD\n");

  printf("    xn %d\n",xn);
  printf("    yn %d\n",yn);
  printf("    rad %f\n",rad);

  //
  //  Amacrine connection list
  //
  n = get_zero_2d_iarray(xn,yn);
  x = get_2d_pointer_iarray(xn,yn);
  y = get_2d_pointer_iarray(xn,yn);

  ir = my_rint(rad + 0.999);    // Radius as an integer, upper bound
  r2 = my_rint(rad*rad);        // Squared radius

  tx = (int *)myalloc(xn*yn*sizeof(int));
  ty = (int *)myalloc(xn*yn*sizeof(int));

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){

      //printf("ir %d   x,y,  %d %d\n",ir,i,j);
      
      // Search within a rectangular region [x0,x1] X [y0,y1]
      x0 = i - ir;
      x1 = i + ir;
      if (x0 < 0)    x0 = 0;
      if (x1 >= xn)  x1 = xn-1;

      y0 = j - ir;
      y1 = j + ir;
      if (y0 < 0)    y0 = 0;
      if (y1 >= yn)  y1 = yn-1;

      //
      //  SHOULD WE (DO WE?) BLock inputs from SELF???
      //
      cn = 0;
      for(ii=x0;ii<=x1;ii++){
	for(jj=y0;jj<=y1;jj++){
	  d2 = (ii-i)*(ii-i) + (jj-j)*(jj-j);
	  if (d2 <= r2){
	    tx[cn] = ii;
	    ty[cn] = jj;
	    cn += 1;
	  }
	}
      }

      x[i][j] = copy_iarray(tx,cn);
      y[i][j] = copy_iarray(ty,cn);
      n[i][j] = cn;
    }
  }

  myfree(tx);
  myfree(ty);

  *rx = x;
  *ry = y;
  *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                        MOD_WAV_DEBUG_WRITE_CONN_IMAGE                     */
/*                                                                           */
/*****************************************************************************/
void mod_wav_debug_write_conn_image(xi,yi)
     int xi,yi;  // Write the weight map for this unit
{
  int i;
  int xn,yn,nbytes,tcode,transpose,pflag;
  int **data,*tcx,*tcy,tcn;

  xn = mod_wav_xn;
  yn = mod_wav_yn;

  data = get_zero_2d_iarray(xn,yn);  // Connection matrix

  tcx = mod_wav_con_aa_x[xi][yi];
  tcy = mod_wav_con_aa_y[xi][yi];
  tcn = mod_wav_con_aa_n[xi][yi];

  for(i=0;i<tcn;i++){
    data[tcx[i]][tcy[i]] += 1;
  }

  nbytes = 4;
  tcode = 1;  // 1-int
  transpose = 1;  // For correct ordering in 'show2d'
  pflag = 0;
  write_2d_data("zz.conn.2d",data,0,0,xn,yn,nbytes,tcode,transpose,pflag);

  free_2d_iarray(data,xn);
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_WAV_FBARS_MAKE_MPT                          */
/*                                                                           */
/*****************************************************************************/
void mod_wav_fbars_make_mpt(m,s,r)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
{
  int i,j,k;
  int stype,xn,yn,*ctx,*cty,ctn;
  short inindex;
  float w,tdelay;
  struct pop_top *mpt;
  struct pop_layer *pl;
  struct pop_syn *tsyn;
  struct pop_cell *pre,*post;
  struct onode *ot,*o;

  mylog(mylogf,"  MOD_WAV_FBARS_MAKE_MPT\n");

  // Basic global params, determine duration of sim in seconds
  mpt = (struct pop_top *)myalloc(sizeof(struct pop_top));

  if  (mylogf == NULL){
    mpt->logf = NULL;
  }else
    mpt->logf = strdup(mylogf);

  mpt->tscale = onode_getpar_flt_exit(m->o,"tscale");
  mpt->sscale = onode_getpar_flt_exit(m->o,"sscale");
  mpt->xn = onode_getpar_int_exit(m->o,"xn");
  mpt->yn = onode_getpar_int_exit(m->o,"yn");
  mpt->tn = onode_getpar_int_exit(m->o,"tn");
  mpt->x1 = 0.0;
  mpt->x2 = mpt->tscale * (float)mpt->tn;
  mpt->out_prefix = NULL;
  if (r != NULL){
    if (r->ppl != NULL)
      mpt->out_prefix = param_getc_dflt(r->ppl,"outfile","zz");
  }

  mpt->binoc = 0;   // Assume monocular left eye
  mpt->retflag = 0; // Default?

  mpt->conn_rw = NULL;
  mpt->conn_read_dir = NULL;

  mpt->narea = 0;
  mpt->area = NULL;
  mpt->nmap = 0;
  mpt->map = NULL;

  //
  //  Layers
  //
  mpt->nlay = 2;
  mpt->lay = (struct pop_layer **)myalloc(mpt->nlay*
					  sizeof(struct pop_layer *));

  pl = popu_make_layer_cart("amacrine",mpt->xn,mpt->yn,1);
  pl->geomt = strdup("default");
  pl->x0 = pl->y0 = 0.0;
  pl->xf = pl->yf = 1.0;
  mpt->lay[0] = pl;
  mpt->lay[0]->ninlist = 1;
  mpt->lay[0]->inlist = (struct onode **)myalloc(sizeof(struct onode *));

  // Create an onode for the input
  ot = paramfile_onode_create_onode();
  ot->otype = strdup("input");
  o = paramfile_onode_get_init_onode("item",NULL,"pop_origin","amacrine",
				     1,NULL,NULL,1);
  onode_insert_child_at_end(ot,o);
  mpt->lay[0]->inlist[0] = ot;

  pl = popu_make_layer_cart("rgc",mpt->xn,mpt->yn,1);
  pl->geomt = strdup("default");
  pl->x0 = pl->y0 = 0.0;
  pl->xf = pl->yf = 1.0;
  mpt->lay[1] = pl;
  mpt->lay[1]->ninlist = 1;
  mpt->lay[1]->inlist = (struct onode **)myalloc(sizeof(struct onode *));

  // Create an onode for the input
  ot = paramfile_onode_create_onode();
  ot->otype = strdup("input");
  o = paramfile_onode_get_init_onode("item",NULL,"pop_origin","amacrine",
				     1,NULL,NULL,1);
  onode_insert_child_at_end(ot,o);
  mpt->lay[1]->inlist[0] = ot;

  stype = 1;
  tdelay = 0.0;
  inindex = 0;
  w = 1.0;

  xn = mod_wav_xn;
  yn = mod_wav_yn;

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){  // For the (i,j)th amacrine cell

      //printf("LayName = %s\n",mpt->lay[0]->name);
      //printf("Cell = %s\n",mpt->lay[0]->c[0][0][0].name);

      //if (mpt->lay[0]->c == NULL)
      //printf("ITS NULLLLLLLLLLLLLLLLLL\n");

      pre = &(mpt->lay[0]->c[i][j][0]);

      //printf("HERE WYET __________  00 2\n");

      ctx = mod_wav_con_aa_x[i][j];   // List of x-coords
      cty = mod_wav_con_aa_y[i][j];   // List of y-coords
      ctn = mod_wav_con_aa_n[i][j];   // Number of connections
      //printf("HERE WYET __________  00 3\n");

      for(k=0;k<ctn;k++){
	// i,j sends input to ctx[k], cty[k]

	//printf("HERE WYET __________  00 4\n");

	post = &(mpt->lay[0]->c[ctx[k]][cty[k]][0]);

	//printf("%d %d %d  making synapse...\n",i,j,k);
	tsyn = pop_cell_add_synapse(pre,post,stype,w,tdelay,inindex);
	//printf("  done\n");

	//  Same input to RGC
	post = &(mpt->lay[1]->c[ctx[k]][cty[k]][0]);
	tsyn = pop_cell_add_synapse(pre,post,stype,w,tdelay,inindex);

	//syn[ctx[k]][cty[k]] += 1.0;   // Active unit contributes to inputs
      }
    }
  }

  popu_mar_write(mpt,m->marfile,4);  // Version 4
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_WAV_FBARS_PREP                            */
/*                                                                           */
/*****************************************************************************/
void mod_wav_fbars_prep(m,s,r)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
{
  int i,j;
  int xn,yn,zn,seed,tau_ti,nunits;
  float rad_pix,sr,rmu,rsd;

  mylog(mylogf,"  MOD_WAV_FBARS_PREP\n");

  m->run_mode = 8; // *** Activate 'indefinite' run mode

  mod_wav_xn         = onode_getpar_int_exit(m->o,"xn");
  mod_wav_yn         = onode_getpar_int_exit(m->o,"yn");
  mod_wav_zn         = onode_getpar_int_exit(m->o,"zn");  // JBCODE
  mod_wav_tn         = onode_getpar_int_exit(m->o,"tn");
  mod_wav_tscale     = onode_getpar_dbl_exit(m->o,"tscale");
  mod_wav_sscale     = onode_getpar_flt_exit(m->o,"sscale");
  mod_wav_umpix      = onode_getpar_flt_exit(m->o,"um_per_pix");
  mod_wav_am_tau     = onode_getpar_flt_exit(m->o,"am_tau");
  mod_wav_am_rad     = onode_getpar_flt_exit(m->o,"am_input_rad");
  mod_wav_am_thr     = onode_getpar_flt_exit(m->o,"am_thresh");
  mod_wav_gc_tau     = onode_getpar_flt_exit(m->o,"gc_tau");
  mod_wav_gc_rad     = onode_getpar_flt_exit(m->o,"gc_input_rad");
  mod_wav_gc_thr     = onode_getpar_flt_exit(m->o,"gc_thresh");
  mod_wav_gc_rate    = onode_getpar_flt_exit(m->o,"gc_rate_active");
  mod_wav_gc_seed    = onode_getpar_int_exit(m->o,"gc_spike_seed");
  mod_wav_refr_mu    = onode_getpar_flt_exit(m->o,"refract_mean");
  mod_wav_refr_sd    = onode_getpar_flt_exit(m->o,"refract_sd");
  mod_wav_refr_seed  = onode_getpar_int_exit(m->o,"refract_seed");
  sr                 = onode_getpar_flt_exit(m->o,"spontaneous_rate");
  mod_wav_act_dur    = onode_getpar_flt_exit(m->o,"activation_dur");
  mod_wav_pop_delay  = onode_getpar_flt_dflt(m->o,"population_delay",0.0);
  mod_wav_init_blank = onode_getpar_flt_dflt(m->o,"initial_blank",0.0); // (s)

  // JBCODE (prev line)
  mod_wav_act_durn   = my_rint(mod_wav_act_dur / mod_wav_tscale);
  mod_wav_spont_prob = sr * mod_wav_tscale;  // prob/s --> prob/tscale

  printf("    act_durn = %d\n",mod_wav_act_durn);
  printf("    spont_prob = %f\n",mod_wav_spont_prob);

  mod_wav_tt = 0;  // Start at time zero
  mod_wav_reset_time = 1; // UNUSED

  // Compute the fractional constant, which if multiplied on every time
  // step will cause exponential decay according to 'mod_wav_am_tau'
  tau_ti = mod_wav_am_tau / mod_wav_tscale;  // Tau in raw time units
  mod_wav_am_fc  = exp(-1.0/tau_ti);
  printf("    _am_fc = %f\n",mod_wav_am_fc);

  tau_ti = mod_wav_gc_tau / mod_wav_tscale;  // Tau in raw time units
  mod_wav_gc_fc  = exp(-1.0/tau_ti);
  printf("    _gc_fc = %f\n",mod_wav_gc_fc);

  xn = mod_wav_xn;
  yn = mod_wav_yn;
  zn = mod_wav_zn;  // JBCODE

  //
  //  Storage for excitation levels
  //
  mod_wav_x_am     = get_zero_2d_farray(xn,yn);
  mod_wav_x_gc     = get_zero_2d_farray(xn,yn);
  mod_wav_synin    = get_zero_2d_farray(xn,yn);
  mod_wav_am_state = get_zero_2d_iarray(xn,yn);

  //
  //  Amacrine refractory periods
  //
  seed = mod_wav_refr_seed;
  if (seed > 0)
    seed *= -1;

  rmu = mod_wav_refr_mu / mod_wav_tscale;
  rsd = mod_wav_refr_sd / mod_wav_tscale;
  printf("    rmu = %f\n",rmu);
  printf("    rsd = %f\n",rsd);


  mod_wav_refr_am = get_zero_2d_iarray(xn,yn);
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      // Note, 'nr_util_gasdev' calls 'nr_util_ran2'
      mod_wav_refr_am[i][j] = my_rint(rmu + rsd * nr_util_gasdev(&seed));
      if (mod_wav_refr_am[i][j] < 0)
	mod_wav_refr_am[i][j] = 0;
    }
  }

  //
  //  Amacrine connection list
  //
  rad_pix = mod_wav_am_rad / mod_wav_umpix;
  mod_wav_make_conn_list_rad(xn,yn,rad_pix,
			     &mod_wav_con_aa_x,
			     &mod_wav_con_aa_y,
			     &mod_wav_con_aa_n);

  //
  //  Write out example(s) of amacrine connection maps, for verification
  //
  //mod_wav_debug_write_conn_image(xn/2,yn/2);


  //
  //  Storage for Ganglion cell spikes
  //
  //  These are somewhat arbitrary, hard-coded guidelines:
  //  (1) Never request more than 200 MByte (thus 50 M ints)
  //  (2) Do not allow more than 10,000 spikes per unit per trial.
  //  (3) Spikes are stored at msec resolution (samp = 1000.0)
  //
  nunits = mod_wav_xn * mod_wav_yn * mod_wav_zn;
  mod_wav_gc_smax = 50000000 / nunits; // 50,000,000 / nunits
  if (mod_wav_gc_smax > 10000){
    mod_wav_gc_smax = 10000;
  }

  // mod_wav_gc_smax = 200; WYETH - DEBUG - to force shorter trials

  printf("    Maximum number of spikes/unit per trial:  %d\n",mod_wav_gc_smax);

  //  Initialize stored spike counts to zero
  mod_wav_gc_cnt = get_zero_3d_iarray(mod_wav_zn,mod_wav_xn,mod_wav_yn);
  mod_wav_gc_spk = (int ****)myalloc(zn*sizeof(int ***));
  for (i=0;i<zn;i++)
    mod_wav_gc_spk[i] = get_3d_iarray(mod_wav_xn,mod_wav_yn,mod_wav_gc_smax);

  mod_wav_gc_samp = 1000.0;

  // Number of GC sampling units per wave simulation time step
  mod_wav_gc_sgen_n = my_rint(mod_wav_gc_samp * mod_wav_tscale);

  printf("    Number of spike sampling units per time step:  %d\n",
	 mod_wav_gc_sgen_n);

  // Correct the 'smax' value so that we never run out of room to store
  // spikes during one 'mod_wav_tscale' epoch.
  mod_wav_gc_smax -= mod_wav_gc_sgen_n;

  mod_wav_gc_prob = mod_wav_gc_rate / mod_wav_gc_samp;

  if (mod_wav_gc_seed > 0){
    mod_wav_gc_seed = -mod_wav_gc_seed;
  }

  //
  //   Write .mar file and return
  //
  if ((m->action == 10) && (m->marfile != NULL)){
    mod_wav_fbars_make_mpt(m,s,r);
    return;
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_WAV_FBARS_DONE                            */
/*                                                                           */
/*****************************************************************************/
void mod_wav_fbars_done()
{
  mylog(mylogf,"  MOD_WAV_FBARS_DONE\n");

  // This currently does nothing.
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_WAV_WRITE_NDATA                            */
/*                                                                           */
/*****************************************************************************/
void mod_wav_write_ndata()
{
  int i,j,k,pi;
  int n,nchan,ti;
  struct ndata_struct *nd;
  char tstr[SLEN];

  printf("  MOD_WAV_WRITE_NDATA\n");

  // Create and fill 'nd'
  nd = get_ndata();
  nd->class = strdup("MODEL");

  // One spike channel for each RGC
  nchan = mod_wav_xn * mod_wav_yn * mod_wav_zn;

  nd->nconst = 0; // Const params
  /*
  nd->cname = (char **)myalloc(nd->nconst*sizeof(char *));
  nd->ctype = (char *)myalloc(nd->nconst*sizeof(char));
  nd->cval = (char **)myalloc(nd->nconst*sizeof(char *));
  i = 0;
  nd->ctype[i] = ;
  nd->cname[i] = ;
  nd->cval[i] = ;
  */


  nd->nvar = 0; // Var params
  /*
  nd->vname = (char **)myalloc(nd->nvar*sizeof(char *));
  nd->vtype = (char *)myalloc(nd->nvar*sizeof(char));
  for(i=0;i<nd->nvar;i++){
    nd->vtype[i] = s->vtype[i];
    nd->vname[i] = strdup(s->vname[i]);
  }
  */

  nd->ntable = 0;

  n = 1;
  nd->ntrial = n;
  nd->t = (struct ndtrial_struct *)myalloc(n*sizeof(struct ndtrial_struct));
  ti = 0; // Trial save index will differ from 'i' for partial save
  nd->t[ti].tcode = i; // might not be 'ti' for partial save
  nd->t[ti].tref = 0;
  nd->t[ti].nparam = nd->nvar;  // 0 here
  nd->t[ti].pname = NULL;
  nd->t[ti].pval = NULL;

  nd->t[ti].nrec = nchan;
  nd->t[ti].r = (struct ndrec_struct *)myalloc(nd->t[ti].nrec*
						 sizeof(struct ndrec_struct));
  k = 0;
  for(pi=0;pi<mod_wav_zn;pi++){  // Population index
    for(i=0;i<mod_wav_xn;i++){
      for(j=0;j<mod_wav_yn;j++){

	// Next line creates storage, otherwise all trials point to one name,
	//   and this creates problems when freeing.

	sprintf(tstr,"rgc_%d_%d",i,j);
	
	nd->t[ti].r[k].name = strdup(tstr);
	nd->t[ti].r[k].rcode = 0;
	nd->t[ti].r[k].sampling = mod_wav_gc_samp;
	nd->t[ti].r[k].t0 = 0;
	nd->t[ti].r[k].tn = 200000; // WYETH HERE FIX
    
	nd->t[ti].r[k].rtype = 0;
	nd->t[ti].r[k].n = mod_wav_gc_cnt[pi][i][j];
	nd->t[ti].r[k].p = mod_wav_gc_spk[pi][i][j];
	k += 1;
      }
    }
  }

  write_ndata("zz.wav.nd",nd);
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_WAV_MAKE_SPIKES                            */
/*                                                                           */
/*****************************************************************************/
void mod_wav_make_spikes(ti,xi,yi,zi)
     int ti;  // Time index for spikes
     int xi;  // Unit x-coord
     int yi;  // Unit y-coord
     int zi;  // Population index
{
  int i;
  int n,tn,*ts,t,delay;

  //t = ti * mod_wav_gc_sgen_n;  // Time offset for spike time
  t = (ti - mod_wav_tt) * mod_wav_gc_sgen_n;  // Time offset for spike time

  ts = mod_wav_gc_spk[zi][xi][yi];
  tn = mod_wav_gc_cnt[zi][xi][yi];

  delay = my_rint(mod_wav_pop_delay * (float)mod_wav_gc_sgen_n /
		  mod_wav_tscale);

  for(i=0;i<mod_wav_gc_sgen_n;i++){
    if (mod_wav_gc_prob > mod_wav_util_ran2(&mod_wav_gc_seed)){
      ts[tn] = t + i + (zi * delay);
      tn += 1;
      if (tn >= mod_wav_gc_smax){

	mod_wav_spike_full = 1;  // This will cause the trial to break
	//
	//  WYETH - WHEN THIS HAPPENS, WE need to start a new trial in the
	//   ndata file, which is a sequential trial - NEW FILE CLASS
	//
	//exit_error("MOD_WAV_MAKE_SPIKES","WYETH - MAKE NEW TRIAL *** FIX");

	//ndata_append_trial(phsav->nd_fname,phsav->t);
	//ndata_overwrite_ntr(phsav->nd_fname,phsav->nd_ntr_addr,phsav->ntr);

      }
    }
  }
  mod_wav_gc_cnt[zi][xi][yi] = tn;  // Update number of spikes

}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_WAV_FBARS_GET_RESPONSE                        */
/*                                                                           */
/*****************************************************************************/
void mod_wav_fbars_get_response(m,s,r)
     struct model_struct *m;        // Model parameter pair list
     struct stim_struct *s;         // Stimulus data and param pair list
     struct response_struct *r;     // Response data
{
  int i,j,k;
  int ti,xn,yn,zn,tn,*ctx,*cty,ctn,actn,**am_state,xi,yi,zi,ns,rmu,rsd,seed;
  int cnt_re,cnt_actnew,cnt_spont,done;
  float amfc,amthr,**xam,**xgc,**syn,***amr;
  float gcfc,gcthr;

  mylog(mylogf,"  MOD_WAV_FBARS_GET_RESPONSE\n");

  xn = mod_wav_xn;
  yn = mod_wav_yn;
  zn = mod_wav_zn;
  tn = mod_wav_tn;

  xam = mod_wav_x_am;
  xgc = mod_wav_x_gc;
  syn = mod_wav_synin;
  am_state = mod_wav_am_state;

  amthr = mod_wav_am_thr;
  amfc  = mod_wav_am_fc;
  gcthr = mod_wav_gc_thr;
  gcfc  = mod_wav_gc_fc;
  actn  = mod_wav_act_durn;  // Number of times steps cell stays active

  rmu = mod_wav_refr_mu / mod_wav_tscale;
  rsd = mod_wav_refr_sd / mod_wav_tscale;
  seed = mod_wav_refr_seed;
  if (seed>0)
    seed *= -1;

  mod_wav_spike_full = 0;
  //
  //  First time processing - only do this at the beginning of the simulation
  //
  if (mod_wav_tt == 0){
    //  WYETH - THIS CAN BE DANGEROUS, because we are assuming that nothing
    //  else will use the 'myrand_util_ran2' if/when this returns a partial
    //  trial in 'indefinite' mode.
    mod_wav_seed = m->mseed[r->tsi];  // Model trial seed
    if (mod_wav_seed > 0)
      mod_wav_seed = -mod_wav_seed;

    //  *** WYETH - CANNOT MAKE A LONG ARRAY, if doing multiple trials
    //
    // WYETH - commenting out AMR
    //amr = (float ***)myalloc(tn*sizeof(float **));

    //  Reset state values
    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){
	xam[i][j] = 0.0;
	xgc[i][j] = 0.0;
	am_state[i][j] = 0;
      }
    }
  }

  //  Reset number of stored spikes
  for(k=0;k<zn;k++){
    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){
	mod_wav_gc_cnt[k][i][j] = 0;
      }
    }
  }

  ti = mod_wav_tt;
  if (ti < tn)
    done = 0;
  else
    done = 1;

  while(done == 0){
    
    for(i=0;i<xn;i++){  // Clear the synaptic input
      for(j=0;j<yn;j++){
	syn[i][j] = 0.0;
      }
    }

    //  These counters are not needed, but useful for debugging
    cnt_re = cnt_actnew = cnt_spont = 0;

    // Compute inputs for this time step
    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){

	//
	//  If inactive units exceed threshold, mark them as active
	//  Or, if they activate spontaneously, mark them as active
	//
	if (am_state[i][j] == 0){  // Not active, not refractory
	  if (xam[i][j] > amthr){
	    am_state[i][j] = actn;  // Number of time steps to stay active
	    cnt_actnew += 1;
	  }else{
	    //
	    //  An inactive, non-refractory unit can activate spontaneously
	    //
	    if (mod_wav_spont_prob > myrand_util_ran2(&mod_wav_seed)){
	      am_state[i][j] = actn;
	      //printf("(ti %d) Spont %d %d  Pr %f\n",ti,i,j,
	      //mod_wav_spont_prob);
	      cnt_spont += 1;
	    }
	  }
	}else if (am_state[i][j] < 0){  // This unit is refractory
	  am_state[i][j] += 1;
	}

	//
	//  Process active units
	//
	if (am_state[i][j] > 0){
	  ctx = mod_wav_con_aa_x[i][j];   // List of x-coords
	  cty = mod_wav_con_aa_y[i][j];   // List of y-coords
	  ctn = mod_wav_con_aa_n[i][j];   // Number of connections
	  for(k=0;k<ctn;k++){
	    syn[ctx[k]][cty[k]] += 1.0;   // Active unit contributes to inputs
	  }
	  am_state[i][j] -= 1;            // Count down to zero
	  if (am_state[i][j] == 0){
	    // Set this unit refractory, period varies across units
	    // am_state[i][j] = -mod_wav_refr_am[i][j];
	    am_state[i][j] = -my_rint(rmu + rsd * nr_util_gasdev(&seed));
	    if (am_state[i][j] > 0)
	      am_state[i][j] = 0;
	    //printf("r %d\n",mod_wav_refr_am[i][j]);
	    cnt_re += 1;
	  }
	}
      }
    }

    // WYETH DEBUG
    //printf("%4d  act %6d  refr %6d  spont %6d\n",ti,
    //cnt_actnew,cnt_re,cnt_spont);

    //
    //  Update unit states based on input and exponential decay
    //
    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){
	xam[i][j] = xam[i][j] * amfc  + syn[i][j];
	xgc[i][j] = xgc[i][j] * gcfc  + syn[i][j];

	// let waves stabilise befor recording spikes
	if (ti >= my_rint(mod_wav_init_blank/mod_wav_tscale)){
	  if (xgc[i][j] > gcthr){
	    for (k=0;k<zn;k++){  // For each population
	      mod_wav_make_spikes(ti,i,j,k); // Generate spikes for this GC
	    }
	  }
	}
      }
    }

    ti += 1;
    if (ti >= tn)
      done = 1;    // End because the full simulation time is finished

    if (mod_wav_spike_full == 1)
      done = 1;          // Return to caller to write out responses so far
    
    // WYETH - commenting out AMR
    //amr[ti] = copy_2d_farray(xam,xn,yn);
  }


  /*  // WYETH - commenting out AMR
  //  Dump Amacrine response to .3d file 
  {
    float ***xyt;
    xyt = get_xyt_from_txy_3d_farray(amr,tn,xn,yn);
    write_3d_data_part("zz.r.3d",xyt,0,xn,0,yn,0,tn,4,2,1);
    free_3d_farray(xyt,xn,yn,tn);
  }
  */

  //
  //  HAND OVER RESPONSE
  //
  for(i=0;i<r->n;i++){
    if (strcmp(r->plname[i],"gc")==0){
      xi = r->xi[i];
      yi = r->yi[i];
      zi = r->zi[i];

      ns = mod_wav_gc_cnt[zi][xi][yi];
      r->s[i][0] = i2farray(mod_wav_gc_spk[zi][xi][yi],ns);
      r->cnt[i][0] = ns;
    }else{
      printf("  r->plname[%d] = %s\n",i,r->plname[i]);
      exit_error("MOD_WAV_...","Unexpected response request");
    }
  }

  //printf("mod_wav_tt = %d \n",mod_wav_tt);
  //printf("mod_wav_tscale = %.20lf \n",mod_wav_tscale);

  //                   (int)         (double)
  r->ttref = (double)mod_wav_tt * mod_wav_tscale;  // Start time in *seconds*
  r->ts0 = 0.0;
  r->tsn = (double)(ti - mod_wav_tt) * mod_wav_tscale;  // Duration (s)

  //printf("r->ttref = %.20lf \n",r->ttref);
  //printf("r->tsn = %.20lf \n",r->tsn);

  mod_wav_tt = ti;  // Start for next trial

  //myfree(amr);  // WYETH - commenting out AMR

  if (ti >= tn)  // The simulation is complete
    m->run_mode = 9; // Stop 'indefinite' run mode
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_WAV_RUN_FBARS                             */
/*                                                                           */
/*****************************************************************************/
void mod_wav_run_fbars(m,s,r,action)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     int action;                // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  if (action == 0){
    mod_wav_fbars_prep(m,s,r);
  }else if (action == 1){
    mod_wav_fbars_get_response(m,s,r);
  }else if (action == -1){
    mod_wav_fbars_done();
  }else if (action == 10){
    mod_wav_fbars_prep(m,s,r);  // r,s are NULL, write MAR file
  }
}
