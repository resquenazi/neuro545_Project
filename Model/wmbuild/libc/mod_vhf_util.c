/*****************************************************************************/
/*                                                                           */
/*   mod_vhf_util.c                                                          */
/*   wyeth bair                                                              */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 2014                                          */
/*    Incorporates C2 layer                                                  */
/*                                                                           */
/*  NOTES                                                                    */
/*   - 'mod_vhf_resp_s2' - see notes w/i to make more efficient.             */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>  // For sleep
#include <math.h>
#include <string.h>

// user defined functions
#include "my_util.h"
#include "nr_util.h"
#include "iarray_util.h"
#include "carray_util.h"
#include "farray_util.h"
#include "myrand_util.h"
#include "data_util.h"
#include "min_util.h"
#include "spike_util.h"
#include "paramfile_util.h"
#include "fft_util.h"
#include "stim_util.h"
#include "stm_util.h"
#include "kernel_util.h"
#include "mod_util.h"
#include "mm_util.h"
#include "pop_cell_util.h"
#include "pop_util.h"      // For .mar file
#include "ifc.h"           // For .mar file
#include "mod.h"           // Data structures
#include "paramfile.h"

#define KSMALL 0.0001   // Cadieu et al (2007) p1734: k=0.0001
#define MAX_PAR  100    // Max number of parameters in model fit (e.g. > 25+3)

// For cluster computation  (static needed w/ initialization)
static int myid = -1;         // Process ID, -1-Single, 0-Master, >0-Slave
static int numproc = -1;      // Number of processors (MM)
static char *mylogf = NULL;   // Logfile name
char ggstr[LONG_SLEN];        // Global log string

int mod_vhf_action;           // Action value from caller

struct vhf_grid{
  char *name;       // Name of this grid
  int   xn;         // Number of units along x-axis
  int   yn;         // Number of units along y-axis
  int   dx;         // x-offset between consecutive units (pix)
  int   dy;         // y-offset between consecutive units (pix)
  int   w;          // Extent of integrated RF, width (pix)
  int   h;          // Extent of integrated RF, height (pix)
  int   ntile;      // Number of tile sizes to integrate
  int  *tile_list;  // [ntile] List of tile size indices to integrate

  //
  //  Derived values
  //
  int    totpixw;   // Number of pixels on x-axis covered by all RFs in grid
  int    totpixh;   // Number of pixels on y-axis covered by all RFs in grid

  int   *s1xn;      // [ntile] Width of S1 unit field (pix)
  int   *s1yn;      // [ntile] Height of S1 unit field (pix)

  int ***s1x0;      // [ntile][xn][yn] lower left S1 unit in RF, x-coord
  int ***s1y0;      // [ntile][xn][yn] lower left S1 unit in RF, y-coord

  float  avgtw;     // Average tile width (pix)

  int    ndx0;      // Start index of this grid within mod_vhf_c1_ndx
};

struct vhf_c1w{
  int gi;            // C1 Grid index
  int xi;            // x-coord
  int yi;            // y-coord
  int ci;            // ori index  /
  float w;           // weight
};

struct vhf_s2_par{
  //
  //  These parameters are constant accross all S2 units
  //
  char  *nltype;        // Type of nonlin, e.g.,
                        //     'sigmoid0'  - Original Cadieu sigmoid
                        //     'none'
  float  s;             // Sigmoid param
  float  a;             // Sigmoid param
  float  b;             // Sigmoid param
  float  eps2;          // Normalization epsilon for S2
                        // *** WYETH see KSMALL, which may be used instead
  int    nw;            // Number of C1 weights
  struct vhf_c1w **c1w; // [nw] pointers to C1 coord and weight
};

struct vhf_s2{
  //
  //  These responses vary across S2 units
  //
  float *****c1rs;  // [nstim][c1grid][x][y][nch] C1 resp, all stim, for fit
};

struct vhf_c2{
  int s2_xn;            // number of S2 units along x-direction
  int s2_dx;            // (pix) lateral offset of neighbor S2 units
  struct vhf_s2 ***s2;  // [s2_xn][s2_xn] pointers to struct w/ c1 responses
};

struct vhf_chans{
  char *name;
  char *type;
  int nori;
  int nph;
  float ph0;
  int maxph_flag;      // 1-max over phase, 0-do not, 2-checkerboard phase max
  int nch;             // number of c1 chans created, either nori, or nph*nori
  int nch0;            // index w/i ...._rs...[*]

  float       f_o;     // Divisor to convert tilew to Gabor SD orth
  float       f_p;     // Divisor to convert tilew to Gabor SD para
  float       f_f;     // Factor to convert tilew to Gabor SF
  float  *****f;       // [size][ori][phase][x][y] Filters
  float ******ft;      // [size][ori][phase][1][x][y] FT of filters
  float  *****fts;     // [size][ori][phase][1][x] FT of filters, extra
  float  *****r;       // [size][ori][phase][x][y] Response

  //
  //  S1 normalization
  //
  int     norm_flag;   // 0-no S1 normalization, 1-local image squared
  float   norm_eps;    // Small number added to norm denominator
  int   **norm_mask_x; // [nsz][...mask_n[]] relative x-coord, norm RF
  int   **norm_mask_y; // [nsz][...mask_n[]] relative x-coord, norm RF
  int    *norm_mask_n; // [nsz] number of points in norm RF
};

//
//  Global variables - model specific
//
int         mod_vhf_xn;      // Width of stimulus (pix)
int         mod_vhf_yn;      // Height of stimulus (pix)
int         mod_vhf_tn;      // Duration of stimulus (frames)
float       mod_vhf_tscale;  // (s/frame)
float       mod_vhf_sscale;  // (deg/pix)
char       *mod_vhf_outfile_prefix;  // Prefix for output files
int         mod_vhf_nsz;     // Number of sizes in S1
int        *mod_vhf_tile_w;  // [nsz] Width of S1 tiles (pix)
int        *mod_vhf_tile_id; // [nsz] "name" of tile size
int        *mod_vhf_tile_x0; // [nsz] Region of S1 units used in computations
int        *mod_vhf_tile_y0; //   defined by the lower left corner (x0,y0)
int        *mod_vhf_tile_xn; //   and the width 'xn' and height 'yn'
int        *mod_vhf_tile_yn; //
int         mod_vhf_ft_resp_flag;  // 1-FFT to get filter responses

//
//  Each of these 'channels' refers to a set of channels.
//    For example, if there is only one set of 4 orientation channels, then 
//    ...s1_nchans is 1 here.
//
int                 mod_vhf_s1_nchans;  // Number of <channels>
struct vhf_chans  **mod_vhf_chans;      // [nchans] list of pointers

//
//  S1 normalization
//
int     mod_vhf_s1_norm_flag;   // 1-compute ..._s1_norm_t2, 0-do not compute
float **mod_vhf_s1_norm_t2;     // [xn][yn] square of current stimulus frame


// These are computed when C1 grids are configured

//
//  C1
//
int         mod_vhf_c1_nch;   // Total number of channels at each C1 grid point
int         mod_vhf_c1_ntot;  // Total number of C1 units across all grids
int         mod_vhf_c1_gridn; // Number of C1 grids
float ******mod_vhf_c1_r;     // [i][j][c1grid][x][y][nch] C1 resp, single stim
                              //    i,j are indices into the S2 grid, eg 3 x 3
struct vhf_grid  **mod_vhf_c1_grid;   // [gridn] pointers

int       **mod_vhf_c1_ndx;   // [_ntot][4] index values for all c1 units
                              //        [0..3]  x, y, ori, grid
int        *mod_vhf_c1_flag;  // [_ntot] flag which c1 units are picked

int         mod_vhf_c1_s1_x0;  // Coords of lower left S1 unit of all
int         mod_vhf_c1_s1_y0;  //   those that are used by any C1 units.
int         mod_vhf_c1_s1_xn;  // Width of S1 region used by C1
int         mod_vhf_c1_s1_yn;  // Height of S1 region used by C1

//
//  S2
//
struct vhf_s2_par  *mod_vhf_s2;      // Current S2 configuration
float         **mod_vhf_s2_r;    // [i][j] S2 resp,single stim

//  C2
struct vhf_c2  *mod_vhf_c2;      // Current S2 configuration
float           mod_vhf_c2_r;    // for storing C2 response  

//
//  Fitting
//
int         mod_vhf_fit_c1_max_n; // Maximum no. of C1 inputs to S2
int         mod_vhf_fit_nc1;      // Number of C1 inputs for current fitting
int         mod_vhf_fit_nrsp;     // Number of responses to fit
float      *mod_vhf_fit_r;        // [..._nrsp] Target response vector
float   ****mod_vhf_fit_c1r = NULL;  // [.._nrsp][s2n][s2n][.._nc1] C1 outputs
int         mod_vhf_fit_xn;     // Number of params being fit
float      *mod_vhf_fit_x;      // [..._xn] Current vector of param values
float      *mod_vhf_fit_dw;     // [nc1] Pre-alloc'd for accumulating derivs.
float    ***mod_vhf_fit_dj;     // [nrsp][s2n][s2n] Norm. const. from C1 resps.

float      *mod_vhf_fit_opt = NULL;  // [MAX_PAR] Best fit params so far
int        *mod_vhf_fit_ndx;    // [nc1] Index values fit so far
float      *mod_vhf_fit_mse;    // [nc1] MSE value after this fit stage

float      *mod_vhf_fit_mse_ts; // [nc1] MSE value on Test data
float      *mod_vhf_fit_rvl;    // [nc1] r-value after this fit stage
float      *mod_vhf_fit_rvl_ts; // [nc1] r-value on Test data

int         mod_vhf_boot_split_i; // index of current split
int        *mod_vhf_boot_shuff;   // [..._nrsp] shuffle index for responses
float      *mod_vhf_fit_train_r;  // [..._train_n] Target response vector
int         mod_vhf_fit_train_n;  // [..._ntrain] Target response vector

//
//  MM Fitting
//
struct vhf_mmfit{
  int    n;        // Number of jobs on fit list
  int    ndone;    // Number of jobs complete (end when this equals 'n')
  int    nci;      // Number of 'ci' values for each fit (e.g., 2)
  int  **ci;       // [n][nci]  Indices of C1 units for this config
  int   *stat;     // [n] Job status:  0-unassigned, 1-assigned, 2-done
  float *mse;      // [n] MSE for done jobs
  float  mse_min;  // Lowest MSE so far
  int    nxti;     // Index of next job to assign
  int   *pstat;    // [numproc] Processor state, 0..n-1 Job# assigned, or -1
};

struct vhf_mmfit *mod_vhf_mmfit;


/**************************************-**************************************/
/*                                                                           */
/*                               MOD_VHF_FNAME                               */
/*                                                                           */
/*  Return an output filename, based on a parameter name 'pname' in an       */
/*  onode 'o'.                                                               */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_fname(o,pname,suffix,rname)
     struct onode *o;  // Onode that contains parameter 'pname'
     char *pname;      // Parameter name
     char *suffix;     // Suffix to add to filename, e.g., ".best.fit"
     char  rname[];    // Write the final name into the string
{
  char *outfile;

  outfile = onode_getpar_chr_dflt(o,pname,"*");

  if (strcmp(outfile,"*")==0)
    sprintf(rname,"%s%s",mod_vhf_outfile_prefix,suffix);
  else if (strcmp(outfile,"NULL")==0){
    rname[0] = '\0';  // NULL
  }else
    sprintf(rname,"%s",outfile);

  myfree(outfile);
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_VHF_CHAN_INDEX                             */
/*                                                                           */
/*  Return two indices corresponding to the global channel index.            */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_chan_index(gci,rci,ri)
     int  gci;  // Global channel index [0..mod_vhf_c1_nch-1]
     int *rci;  // Index in 'mod_vhf_chans', i.e., [0..mod_vhf_s1_nchans]
     int *ri;   // Index within <channel> set
{
  int i;
  int done,max_ci;

  *rci = -1;
  *ri = -1;

  i = 0;
  done = 0;
  while(done == 0){
    //
    //  'i' runs across the sets of channels.
    //

    // For this, the ith set, of channels, the maximum 'ci' value would be
    //   the starting channel number (nch0) plus the number (nch) - 1
    max_ci = mod_vhf_chans[i]->nch0 + mod_vhf_chans[i]->nch - 1;
    if (gci <= max_ci){
      *rci = i;

      // *** WYETH BUG - May 30 2016
      // *** WYETH BUG - I think this reverses the indices of the channels
      // *** WYETH BUG - I think this reverses the indices of the channels
      // *** WYETH BUG - I think this reverses the indices of the channels

      // WYETH OLD WRONG - BUG
      //*ri = max_ci - gci; // *** WYETH IS THIS RIGHT ???? THe 'k' values
                        //  in the PP plot maker seems to get wrong ori indices

      *ri = gci - mod_vhf_chans[i]->nch0;  // NEW FIXED May 20 2016

      done = 1;
    }else
      i += 1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_VHF_NDX_FOR_GRID                           */
/*                                                                           */
/*  Given a C1 grid index 'gi' return the starting value for this grid in    */
/*  the global index list 'mod_vhf_c1_ndx' and the number of points in the   */
/*  list for this grid.                                                      */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_ndx_for_grid(gi,r0,rn)
     int gi;   // Index of grid
     int *r0;  // Return first index in 'mod_vhf_c1_ndx' for grid 'gi'
     int *rn;  // Return number of points in 'mod_vhf_c1_ndx' for grid 'gi'
{
  int i;
  int gxn,gyn,c0,cn;

  if ((gi < 0) || (gi >= mod_vhf_c1_gridn))
    exit_error("MOD_VHF_NDX_FOR_GRID","Grid index out of range");

  c0 = 0;
  for(i=0;i<gi;i++){
    gxn = mod_vhf_c1_grid[i]->xn;      // width of grid
    gyn = mod_vhf_c1_grid[i]->yn;      // height of grid
    c0 += gxn * gyn * mod_vhf_c1_nch;  // Total number of C1 units
  }

  gxn = mod_vhf_c1_grid[gi]->xn;       // width of grid
  gyn = mod_vhf_c1_grid[gi]->yn;       // height of grid
  cn  = gxn * gyn * mod_vhf_c1_nch;    // Total number of C1 units

  *rn = cn;
  *r0 = c0;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_VHF_C1_NDX_GET                            */
/*                                                                           */
/*  Return the index into 'mod_vhf_c1_ndx' given the grid, position and      */
/*  channel indices.                                                         */
/*                                                                           */
/*****************************************************************************/
int mod_vhf_c1_ndx_get(gi,xi,yi,ci)
     int gi;  // grid index
     int xi;  // x-index w/i grid 'gi'
     int yi;  // y-index w/i grid 'gi'
     int ci;  // channel index at this grid point
{
  int i;
  int i0,gxn,gyn,nch;

  i0  = mod_vhf_c1_grid[gi]->ndx0;
  gxn = mod_vhf_c1_grid[gi]->xn;
  gyn = mod_vhf_c1_grid[gi]->yn;
  nch = mod_vhf_c1_nch;

  i = i0 + gxn*(gyn*nch) + gyn*nch + ci;

  //
  //  WYETH - Test this, can remove after a while
  //
  if ((mod_vhf_c1_ndx[i][0] != xi)||
      (mod_vhf_c1_ndx[i][1] != yi)||
      (mod_vhf_c1_ndx[i][2] != ci)||
      (mod_vhf_c1_ndx[i][3] != gi))
    exit_error("MOD_VHF_C1_NDX_GET","Mismatch - should not happen");

  return i;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_VHF_C1_NDX_PRINT                           */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_c1_ndx_print(k)
     int k;  // Global channel index
{
  printf("      GCI %d =  Grid %d  (%d, %d)  Chan %d\n",k,mod_vhf_c1_ndx[k][3],
	 mod_vhf_c1_ndx[k][0],mod_vhf_c1_ndx[k][1],mod_vhf_c1_ndx[k][2]);
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_VHF_C1_RESP_GCI                            */
/*                                                                           */
/*  Return the array of responses for channel 'gci' across all stimuli.      */
/*                                                                           */
/*****************************************************************************/
float *mod_vhf_c1_resp_gci(s2i,s2j,gci)
     int s2i,s2j;  // Index into S2 grid
     int gci;      // Global channel index
{
  int xi,yi,ci,gi,i;
  float *tdata,*****trs;

  trs = mod_vhf_c2->s2[s2i][s2j]->c1rs;

  tdata = get_farray(mod_vhf_fit_nrsp);

  xi = mod_vhf_c1_ndx[gci][0];
  yi = mod_vhf_c1_ndx[gci][1];
  ci = mod_vhf_c1_ndx[gci][2];
  gi = mod_vhf_c1_ndx[gci][3];

  for(i=0;i<mod_vhf_fit_nrsp;i++)  // For each stimulus
    tdata[i] = trs[i][gi][xi][yi][ci]; // [nstim][c1grid][x][y][chn]

  return tdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_VHF_RELATE_C1_TARG                          */
/*                                                                           */
/*  Examine statistical relationships between the c1 responses and the       */
/*  target neuronal responses.                                               */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_relate_c1_targ()
{
  int i;
  int nst,nc1,vmax_i,vmin_i,s2i,s2j;
  float *targ,*tc1,*c1vtarg,x,vmax,vmin;
  char tstr[SLEN];

  printf("  MOD_VHF_RELATE_C1_TARG\n");

  s2i = s2j = mod_vhf_c2->s2_xn / 2;
  printf("    *** WARNING:  Using s2 unit at location %d %d\n",s2i,s2j);

  nst   = mod_vhf_fit_nrsp;  // # of stimuli = # of target responses
  nc1   = mod_vhf_c1_ntot;   // Total number of C1 units (across all grids)

  targ = mod_vhf_fit_r;       // [nst] Pointer to target responsesa

  c1vtarg = get_farray(nc1);  // [nc1] Store r-value here


  remove_file("zzz.c1.table.txt");
  for(i=0;i<nc1;i++){
    sprintf(tstr,"%4d   Grid %d  (%d, %d)  Chan %d\n",i,mod_vhf_c1_ndx[i][3],
	    mod_vhf_c1_ndx[i][0],mod_vhf_c1_ndx[i][1],mod_vhf_c1_ndx[i][2]);
    append_string_to_file("zzz.c1.table.txt",tstr);
  }

  //
  //  Compute r-values between the C1 responses and the target responses
  //
  vmax = -2.0;  // Smaller than anything
  vmin =  2.0;  // Larger than anything
  for(i=0;i<nc1;i++){
    tc1 = mod_vhf_c1_resp_gci(s2i,s2j,i);

    if (sum_farray(tc1,nst,0,nst) == 0){
      printf("  *** Channel %d has all zero responses.\n",i);
      x = -2.0;  // Flag value
    }else{
      simple_pearsn(targ-1,tc1-1,nst,&x);  // Arrays must start at [1]

      if (x > vmax){
	vmax = x;
	vmax_i = i;
      }else if (x < vmin){
	vmin = x;
	vmin_i = i;
      }
    }

    c1vtarg[i] = x;


    myfree(tc1);
  }
  printf("    Largest r-value:   %5.2f  (gci = %d)\n",vmax,vmax_i);
  mod_vhf_c1_ndx_print(vmax_i);
  printf("    Smallest r-value:  %5.2f  (gci = %d)\n",vmin,vmin_i);
  mod_vhf_c1_ndx_print(vmin_i);

  remove_file("zzz.c1vtarg.pl");
  append_farray_plot("zzz.c1vtarg.pl","r__c1_vs_targ",c1vtarg,nc1,1);


  //
  //   Append the scatter plot for the largest positive correlation
  //
  tc1 = mod_vhf_c1_resp_gci(s2i,s2j,vmax_i);
  sprintf(tstr,"c1_%d__v__targ",vmax_i);
  append_farray_xy_plot("zzz.c1vtarg.pl",targ,tc1,nst,tstr);
  myfree(tc1);

  exit(0);
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_VHF_NORM_MASK                             */
/*                                                                           */
/*  Return coordinates for all points that lie within the disk of radius     */
/*  'maskr'.                                                                 */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_norm_mask(maskr,rxi,ryi,rn)
     float maskr;  // Mask radius (pix)
     int **rxi;    // [*rn] x-coords for mask
     int **ryi;    // [*rn] y-coords for mask
     int  *rn;     // Number of points in mask
{
  int i,j,k;
  int *xi,*yi,d2,dmax,n,r;

  dmax = (int)(maskr*maskr);
  r = my_rint(maskr);

  //
  //  Count the number of points within the radius 'r'
  //
  n = 0;
  for(i=-r;i<=r;i++){
    for(j=-r;j<=r;j++){
      d2 = i*i + j*j;
      if (d2 <= dmax)
	n += 1;
    }
  }

  xi = (int *)myalloc(n*sizeof(int));
  yi = (int *)myalloc(n*sizeof(int));

  //
  //  Save coordinates for all points in the disk
  //
  k = 0;
  for(i=-r;i<=r;i++){
    for(j=-r;j<=r;j++){
      d2 = i*i + j*j;
      if (d2 <= dmax){
	xi[k] = i;
	yi[k] = j;
	k += 1;
      }
    }
  }

  *rxi = xi;
  *ryi = yi;
  *rn  = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_VHF_RESP_OP                              */
/*                                                                           */
/*  Initialize or reset the S1 and C1 response storage.                      */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_resp_op(optype,nn)
     char *optype;      // "s1_init", "c1s_init"
     int nn;            // integer parameter
{
  int ii,jj,ci,i,j,k;
  int gxn,gyn,nsz,nori,nph,s2n;
  float *****r;

  if (strcmp(optype,"s1_init")==0){
    //
    //  Allocate memory for S1 response storage
    //
    for(ci=0;ci<mod_vhf_s1_nchans;ci++){

      nsz  = mod_vhf_nsz;                // *** 'nsz' is still constant
      nori = mod_vhf_chans[ci]->nori;
      nph  = mod_vhf_chans[ci]->nph;

      r = (float *****)myalloc(nsz*sizeof(float ****));
      for(i=0;i<nsz;i++){
	r[i] = (float ****)myalloc(nori*sizeof(float ***));
	for(j=0;j<nori;j++){
	  r[i][j] = (float ***)myalloc(nph*sizeof(float **));
	  for(k=0;k<nph;k++){
	    r[i][j][k] = get_2d_farray(mod_vhf_xn,mod_vhf_yn);
	  }
	}
      }

      mod_vhf_chans[ci]->r = r;

      //
      //  Storage for S1 normalization, to hold squared stimulus values
      //
      mod_vhf_s1_norm_t2 = get_2d_farray(mod_vhf_xn,mod_vhf_yn);
    }

  }else if (strcmp(optype,"c1s_init")==0){
    //
    //  Allocate memory for C1 response storage
    //
    s2n = mod_vhf_c2->s2_xn;

    for(ii=0;ii<s2n;ii++){
      for(jj=0;jj<s2n;jj++){

	r = (float *****)myalloc(nn*sizeof(float ****));  // Number of stim.
	for(i=0;i<nn;i++){
	  r[i] = (float ****)myalloc(mod_vhf_c1_gridn*sizeof(float ***));
	  for(j=0;j<mod_vhf_c1_gridn;j++){
	    gxn = mod_vhf_c1_grid[j]->xn;
	    gyn = mod_vhf_c1_grid[j]->yn;
	    
	    r[i][j] = get_3d_farray(gxn,gyn,mod_vhf_c1_nch);
	  }
	}

	mod_vhf_c2->s2[ii][jj]->c1rs = r;   // [nstim][c1grid][gx][gy][nch]
      }
    }
  }else{
    printf("  optype:  %s\n",optype);
    mylog_exit(mylogf,"MOD_VHF_RESP_OP  Unknoown optype.");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_VHF_RESP_S1_NORM                           */
/*                                                                           */
/*  Normalize the responses in 'mod_vhf_chans[]->r'                          */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_resp_s1_norm(ci)
     int ci;  // Index to a set of channels
{
  int i,j,k;
  int oi,pi,si,xi,yi,xn,yn,nsz,nori,nph,x0,x1,y0,y1,*maskx,*masky,maskn;
  int s2n,s2dx,adj0,adj1;
  float ****tr1,tot,eps;
  struct vhf_chans *chns;

  chns = mod_vhf_chans[ci];

  if (chns->norm_flag != 1)
    exit_error("MOD_VHF_RESP_S1_NORM","Mask norm flag - bad value");

  xn   = mod_vhf_xn;
  yn   = mod_vhf_yn;
  nsz  = mod_vhf_nsz;

  nori = chns->nori;
  nph  = chns->nph;
  eps  = chns->norm_eps;

  s2n  = mod_vhf_c2->s2_xn;     // Size of S2 grid, e.g., 3 (x 3)
  s2dx = mod_vhf_c2->s2_dx;     // E.g. 30 pix
  // Adjust for S2 grid
  adj0 = -s2dx * (int)(s2n/2);  // Add this to start index
  adj1 =  s2dx * (int)(s2n-1);  // Add this to length of region
  //printf(" adjust________adj0,1 =  %d %d\n",adj0,adj1);


  for(k=0;k<nsz;k++){  // For each size

    // Get the norm mask points for this size
    maskx = chns->norm_mask_x[k]; // x-coords of mask points
    masky = chns->norm_mask_y[k]; // y-coords of mask points
    maskn = chns->norm_mask_n[k]; // Number of mask points

    tr1 = chns->r[k];  // Ptr to responses for this size [ori][phase][x][y]

    //x0  = mod_vhf_tile_x0[k];  // NOTE, This _x0[] is used here only
    //x1  = x0 + mod_vhf_tile_xn[k] - 1;
    //y0  = mod_vhf_tile_y0[k];
    //y1  = y0 + mod_vhf_tile_yn[k] - 1;
    x0  = adj0 + mod_vhf_tile_x0[k];  // NOTE, This _x0[] is used here only
    x1  = x0 + adj1 + mod_vhf_tile_xn[k] - 1;
    y0  = adj0 + mod_vhf_tile_y0[k];
    y1  = y0 + adj1 + mod_vhf_tile_yn[k] - 1;

    // ***** WYETH BUG HERE  **************
    // ***** WYETH BUG HERE  **************
    // ***** WYETH BUG HERE  **************  *** These ranges only cover the
    // ***** WYETH BUG HERE  **************      units for the ***middle***
    // ***** WYETH BUG HERE  **************  S2 UNIT ***********

    //printf("TILE REGION:  (%d)  x0,1 %d %d   y0,1 %d %d\n",k,x0,x1,y0,y1);

    for(i=x0;i<=x1;i++){
      for(j=y0;j<=y1;j++){  // For each S1 location used for this size

	if ((i < 0) || (j < 0)){
	  printf("  k = %d\n",k);
	  printf("  mod_vhf_tile_x0[k] = %d\n",mod_vhf_tile_x0[k]);
	  printf("  adj0 = %d\n",adj0);
	  printf("  i = %d   j = %d\n",i,j);
	  exit_error("MOD_VHF_RESP_S1_NORM","index error");
	}

	// 1D loop over the 2D normalization region
	tot = 0.0;
	for(si=0;si<maskn;si++){  // For each point in the mask
	  xi = i + maskx[si];
	  yi = j + masky[si];
	  tot += mod_vhf_s1_norm_t2[xi][yi];  // coords relative to (i,j)
	}
	//tot /= (float)maskn;  //*** THIS NOT NEEDED if S1 filters not norm'd
	tot = sqrt(tot);
	tot += eps; // To avoid divide by zero

	/*
	if ((i==(xn/2))&&(j==(yn/2)))
	  printf("tot = %f  (n=%d)  tr1[oi][pi][%d][%d] = %f\n",tot,maskn,
		 i,j,tr1[0][0][i][j]);
	*/
	//printf("tot = %f tr1[oi][pi][64][64] = %f\n",tot,tr1[0][0][64][64]);

	// Apply this norm value to all oris and phases at this location
	for(oi=0;oi<nori;oi++){
	  for(pi=0;pi<nph;pi++){
	    tr1[oi][pi][i][j] /= tot;
	    //tr1[oi][pi][i][j] = 100.0;  // DEBUG
	  }
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_VHF_FFT_RESP                             */
/*                                                                           */
/*  Fill '...chans->r' with the response of all filters to a single 2D stim  */
/*  frame that is held within a 3D array 'stim1'.                            */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_fft_resp(stim1)
     float ***stim1;  // One frame from stimulus: [1][1..xn][1..yn]
                      //   Thus, this is a 2D frame w/i a 3D 'f3tensor'
{
  int i,j,k;
  int ci,xi,yi,xn,yn,tn,nsz,nori,nph;
  float **speq,***ft,**fts,fc,***r,**rs,**tr;
  struct vhf_chans *chns;

  xn   = mod_vhf_xn;
  yn   = mod_vhf_yn;
  tn   = mod_vhf_tn;
  nsz  = mod_vhf_nsz;


  fc = 2.0/((float)(xn*yn));

  //  Compute the FFT of the stimulus 'stim1'
  speq = matrix(1,1,1,2*xn);
  contort_real_3d_farray(stim1,1,1,1,1,xn,yn);
  rlft3(stim1,speq,1,xn,yn,1);


  for(ci=0;ci<mod_vhf_s1_nchans;ci++){  // For each set of S1 channels

    chns = mod_vhf_chans[ci];
    nori = chns->nori;
    nph  = chns->nph;

    //
    //  For each filter, compute product of FT with stim FT, inverse transform
    //
    for(i=0;i<nsz;i++){
      for(j=0;j<nori;j++){
	for(k=0;k<nph;k++){

	  ft  = chns->ft[i][j][k];
	  fts = chns->fts[i][j][k];

	  three_d_fft_prod(stim1,speq,ft,fts,1,xn,yn,&r,&rs);
	  // Note, r[1][1..xn][1..yn]
	  rlft3(r,rs,1,xn,yn,-1);
	  free_matrix(rs,1,1,1,2*xn);
	  multiply_3d_farray(r,1,1,1,xn,1,yn,fc);
	  contort_real_3d_farray(r,1,1,1,1,xn,yn);

	  //
	  //  Copy the response into global array
	  //
	  tr = chns->r[i][j][k];  // Storage already exists here
	  for(xi=0;xi<xn;xi++){
	    for(yi=0;yi<yn;yi++){
	      tr[xi][yi] = r[1][xi+1][yi+1];
	    }
	  }

	  free_f3tensor(r,1,1,1,xn,1,yn);
	}
      }
    }
  }
  free_matrix(speq,1,1,1,2*xn);
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_VHF_RESP_S1                              */
/*                                                                           */
/*  Compute the S1 response, setting values in mod_vhf_s1_r                  */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_resp_s1(stim1)
     float ***stim1;   // NOTE This stimulus may be modified here
{
  int i,j;

  //
  //  (1) Compute stimulus squared for normalization in step (3), if needed.
  //
  if (mod_vhf_s1_norm_flag == 1){
    for(i=0;i<mod_vhf_xn;i++){
      for(j=0;j<mod_vhf_yn;j++){
	mod_vhf_s1_norm_t2[i][j] = stim1[1][1+i][1+j] * stim1[1][1+i][1+j];
      }
    }
  }

  //
  //  (2) Compute the filter responses.
  //
  if (mod_vhf_ft_resp_flag == 1){
    mod_vhf_fft_resp(stim1);  // 'stim1' is modified here
  }else
    mylogx(ggstr,"MOD_VHF_01_GET_RESPONSE",
	   "non-FFT method not implemented yet");

  //
  //  (3) Normalize the S1 responses, if indicated.
  //
  for(i=0;i<mod_vhf_s1_nchans;i++)  // For each set of S1 channels
    if (mod_vhf_chans[i]->norm_flag != 0)
      mod_vhf_resp_s1_norm(i);  // Norm the S1 responses

}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_VHF_RESP_C1_ALL                           */
/*                                                                           */
/*  Compute responses for all C1 units.                                      */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_resp_c1_all(resp,ii,jj)
     float ****resp;   // [c1grid][x][y][chn] Store responses here
     int ii;           // index in S2 grid, for computing offset
     int jj;           // index in S2 grid, for computing offset
{
  int i,j,k,l;
  int ti,xi,yi,si,pi,ci,gn,gxn,gyn,x0,x1,y0,y1,tin,nori,nph,s2xn,x_offset,ch0;
  int odd_point;
  float smax,**rxy;
  struct vhf_grid *gt;
  struct vhf_chans *chns;

  //printf("  MOD_VHF_RESP_C1_ALL\n");

  s2xn     = mod_vhf_c2->s2_xn;  // E.g. 3 x 3 grid of S2 units
  x_offset = mod_vhf_c2->s2_dx;  // E.g. 30 pix
  gn       = mod_vhf_c1_gridn;   // Number of C1 grids

  // For each set of channels
  for(ci=0;ci<mod_vhf_s1_nchans;ci++){

    chns = mod_vhf_chans[ci];
    nori = chns->nori;
    nph  = chns->nph;
    ch0 = chns->nch0;

    for(i=0;i<gn;i++){  // For each C1 grid
      gt = mod_vhf_c1_grid[i];
      gxn = gt->xn;
      gyn = gt->yn;
      tin = gt->ntile;

      for(j=0;j<gxn;j++){
	for(k=0;k<gyn;k++){

	  // WYETH - Mar 15, 2018 - This 'odd_point' flag is being added to
	  // allow for not maxing over phase by randomly choosing one phase
	  // or the other in each channel
	  //
	  odd_point = (j+k)%2;  // 0-even 1-odd; Checkerboard pattern
	  //printf("%d %d   odd: %d\n",j,k,odd_point);

	  for(l=0;l<nori;l++){  // For each orientation at each grid location
	    //
	    //  Find max response over all tiles, phases, and (x,y) positions
	    //
	      
	    // *** WYETH HERE - add a flag array condition for only those
	    //                  C1 units that have non-zero weights to S2
	      
	    // Note, If maxing over opposite phases, then the max can never
	    //   be negative.

	    // Can there not be negative responses???
	    smax = 0.0;  // In case we are maxing over phase, 0 it now.

	    odd_point = 1 - odd_point;  // Alternate 0 and 1 w/ each ori

	    for(pi=0;pi<nph;pi++){  // For each phase

	      if ((chns->maxph_flag != 2) ||
		  ((chns->maxph_flag == 2) && ((pi%2) == odd_point))){
		//  if maxph_flag is 2, then one do the even or odd index phase


		// ******* WYETH TESTING -------- REMOVE
		/*
		if (chns->maxph_flag == 2){
		  printf("CHECKERBOARD:  %d grid  (x,y,ori): %d %d %d  max over phase %d\n",i,j,k,l,pi);
		}*/


		if (chns->maxph_flag == 0)
		  smax = 0.0;

		for(ti=0;ti<tin;ti++){  // For each entry in the tile list
		  // Get the index of this tile in the list of all S1 tiles
		  si = get_index_search_iarray(mod_vhf_tile_id,mod_vhf_nsz,
					       gt->tile_list[ti]);

		  // Coords defining rectangle of S1 inputs for this C1 location
		  if ((s2xn != 3)&&(s2xn != 0)&&(s2xn != 1)&&(s2xn != 5))
		    //
		    //  *** WYETH BUG HERE, I think, this shouldn't have
		    //      half values (rounded to int) when s2xn = 3 ????
		    //      Works for even but not odd??
		    exit_error("WYETH EXIT","Must fix this ???");
		  x0 = gt->s1x0[ti][j][k] + (ii-(s2xn/2))*x_offset;
		  x1 = x0 + gt->s1xn[ti] - 1;
		  //if ((ti == 0)&&(pi==0)&&(l==0)&&(k==0)&&(j==(gxn-1)))
		  //printf("ii,jj %d %d  s2xn %d  x_offs %d   x0,x1 %d %d\n",
		  //ii,jj,s2xn,x_offset,x0,x1);

		  y0 = gt->s1y0[ti][j][k] + (jj-(s2xn/2))*x_offset;;
		  y1 = y0 + gt->s1yn[ti] - 1;

		  if ((x0 < 0)||(y0 < 0))
		    exit_error("MOD_VHF_RESP_C1_ALL",
			       "Index into response array is negative");

		  // Convenient pointer
		  rxy = chns->r[si][l][pi];  // [size][ori][phase][x][y] Resp.

		  for(xi=x0;xi<=x1;xi++){  // For each column in S1 inputs
		    for(yi=y0;yi<=y1;yi++){  // For each row in S1 inputs
		      if (rxy[xi][yi] > smax){
			smax = rxy[xi][yi];
		      }
		    }
		  }
		}

		if (chns->maxph_flag == 0){
		  resp[i][j][k][ch0+l*nph+pi] = smax;  // [c1grid][x][y][chn]
		}
	      }
	    } // end for 'pi' each phase


	    //if (chns->maxph_flag == 1){
	    if ((chns->maxph_flag == 1)||(chns->maxph_flag == 2)){
	      resp[i][j][k][ch0+l] = smax;  // [c1grid][x][y][chn]

	      //printf("  smax[%d][%d][%d][%d] = %f\n",i,j,k,ch0+l,smax);
	      if (ch0 + l >= mod_vhf_c1_nch)
		exit_error("WYETH","This should not happen");
	    }
	  }
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_VHF_FILE_C1                              */
/*                                                                           */
/*  Read / write C1 responses from / to a file.                              */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_file_c1(fito,fmode)
     struct onode *fito;
     char *fmode;
{
  FILE *fopen(),*fin,*fout;
  int i,j,k,ii,jj;
  int si,gi,gn,gxn,gyn,nch,nstm,rv,s2n,ti;
  float ***r;
  char *fname,ts[SLEN],vers[SLEN];

  mylog(mylogf,"  MOD_VHF_FILE_C1\n");

  fname = onode_getpar_chr_exit(fito,"file_c1");

  nch  = mod_vhf_c1_nch;
  nstm = mod_vhf_fit_nrsp;
  gn   = mod_vhf_c1_gridn;
  s2n  = mod_vhf_c2->s2_xn;

  if (strcmp(fmode,"write")==0){

    if ((fout = fopen(fname,"w")) == NULL){
      sprintf(ggstr,"  *** %s\n",fname);
      mylog(mylogf,ggstr);
      mylogx(mylogf,"MOD_VHF_FILE_C1","Could not open file to write");
    }
    sprintf(ggstr,"    Write C1 responses to file:  %s\n",fname);
    mylog(mylogf,ggstr);

    fprintf(fout,"Version 2.0\n");
    fprintf(fout,"Stimuli %d\n",nstm);
    fprintf(fout,"Grids %d\n",gn);
    //fprintf(fout,"C1_tot %d\n",-1);  // WYETH ADD THIS ???
    fprintf(fout,"S2_size %d\n",s2n);

    for(si=0;si<nstm;si++){  // For each stimulus
      fprintf(fout,"Stimulus %d\n",si);

      for(ii=0;ii<s2n;ii++){
	for(jj=0;jj<s2n;jj++){  // For each S2 unit
	  for(gi=0;gi<gn;gi++){  // For each grid
	    gxn = mod_vhf_c1_grid[gi]->xn;   // width of grid
	    gyn = mod_vhf_c1_grid[gi]->yn;   // height of grid
	    r = mod_vhf_c2->s2[ii][jj]->c1rs[si][gi]; // [nstim][c1grid][...]
	    for(i=0;i<gxn;i++){
	      for(j=0;j<gyn;j++){
		for(k=0;k<nch;k++){
		  fprintf(fout," %f",r[i][j][k]);
		}
	      }
	    }
	  }
	}
      }

      fprintf(fout,"\n");
    }
    fclose(fout);
  }else if (strcmp(fmode,"read")==0){

    if ((fin = fopen(fname,"r")) == NULL){
      sprintf(ggstr,"  *** %s\n",fname);
      mylog(mylogf,ggstr);
      mylogx(ggstr,"MOD_VHF_FILE_C1","Could not open file to read");
    }
    sprintf(ggstr,"    Reading C1 responses from file:  %s\n",fname);
    mylog(mylogf,ggstr);

    rv = fscanf(fin,"%s %s",ts,vers);
    check_string_exit_carray(ts,"Version");
    rv = fscanf(fin,"%s %d",ts,&nstm);
    check_string_exit_carray(ts,"Stimuli");
    rv = fscanf(fin,"%s %d",ts,&gn);
    check_string_exit_carray(ts,"Grids");

    if (strcmp(vers,"1.0")==0){
      s2n = 1;
    }else{
      rv = fscanf(fin,"%s %d",ts,&ti);
      check_string_exit_carray(ts,"S2_size");
      if (ti != s2n){
	mylogx(mylogf,"MOD_VHF_FILE_C1","Mismatch of S2 grid size");
      }
    }

    for(si=0;si<nstm;si++){  // For each stimulus

      rv = fscanf(fin,"%s %d",ts,&si);  // Stimulus index
      check_string_exit_carray(ts,"Stimulus");

      for(ii=0;ii<s2n;ii++){
	for(jj=0;jj<s2n;jj++){  // For each S2 unit

	  for(gi=0;gi<gn;gi++){  // For each grid
	    gxn = mod_vhf_c1_grid[gi]->xn;   // width of grid
	    gyn = mod_vhf_c1_grid[gi]->yn;   // height of grid

	    r = mod_vhf_c2->s2[ii][jj]->c1rs[si][gi]; // [nstim][c1grid][...]

	    for(i=0;i<gxn;i++){
	      for(j=0;j<gyn;j++){
		for(k=0;k<nch;k++){
		  rv = fscanf(fin,"%f",&(r[i][j][k]));
		}
	      }
	    }
	  }
	}
      }
    }
    fclose(fin);
  }

  myfree(fname);
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_VHF_RESP_S2                              */
/*                                                                           */
/*  Response of a "Selectivity" unit.                                        */
/*                                                                           */
/*****************************************************************************/
float mod_vhf_resp_s2(nltype,s,a,b,eps,w,x,n,ru)
     char *nltype; // 'sigmoid0', 'none', etc...
     float s,a,b;  // Three sigmoid params
     float eps;    // Normalization constant (to avoid division by zero)
     float *w;     // [n] Weights for C1 inputs 
     float *x;     // [n] Drive of C1 inputs
     int n;        // Number of C1 inputs
     float *ru;    // Return the normalized value, before the sigmoid
{
  int i;
  float sum,sumx2,u,gu;

  // *** WYETH - WHY NOT USE precomputed normalization values, ...dj here?
  // *** WYETH - WHY NOT USE precomputed normalization values, ...dj here?
  // *** WYETH - WHY NOT USE precomputed normalization values, ...dj here?

  //
  //  Eqn 1 (Cadieu et al., 2007)
  //
  sum = 0.0;
  sumx2 = 0.0;
  for(i=0;i<n;i++){  // For each weighted input
    sum += w[i] * x[i];    // SUM w*x
    sumx2 += x[i] * x[i];  // SUM x^2
  }
  u = sum / (sqrt(sumx2) + eps);

  // "sigmoid0"
  if (strcmp(nltype,"sigmoid0")==0){
    gu = s / (1.0 + exp(-a*(u-b)));
  }else{
    gu = u;
  }
    
  *ru = u;  // Also, return the value before it goes through the sigmoid

  return gu;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_VHF_RESP_C2_MAXIJ                          */
/*                                                                           */
/*  Determine which S2 unit has the maximum response.                        */
/*                                                                           */
/*****************************************************************************/
float mod_vhf_resp_c2_maxij(s,a,b,eps,w,x,s2n,n,ri,rj)
     float s,a,b;  // Three sigmoid params
     float eps;    // Normalization constant (to avoid division by zero)
     float *w;     // [n] Weights for C1 inputs 
     float ***x;   // [s2n][s2n][n] Drive of C1 inputs
     int s2n;      // dimension of S2 grid
     int n;        // Number of C1 inputs
     int *ri,*rj;  // Return the S2 coords of the max response
{
  int i,ii,jj;
  int maxi,maxj,midi;
  float sum,sumx2,u,gu,s2max,*tx;

  // *** WYETH NEW - WHY NOT USE precomputed normalization values, ...dj here?
  // *** WYETH NEW - WHY NOT USE precomputed normalization values, ...dj here?
  // *** WYETH NEW - WHY NOT USE precomputed normalization values, ...dj here?

  s2max = -1000.0;  // *** WYETH - MUST BE LOWER THAN ANY POSSIBLE VALUE

  midi = s2n/2;

  maxi = maxj = -1;

  for(ii=0;ii<s2n;ii++){
    for(jj=0;jj<s2n;jj++){

      tx = x[ii][jj];  // Pointer to the C1 inputs for this S2 unit

      //
      //  Eqn 1 (Cadieu et al., 2007)
      //
      sum = 0.0;
      sumx2 = 0.0;
      for(i=0;i<n;i++){  // For each weighted input
	sum += w[i] * tx[i];     // SUM w*x
	sumx2 += tx[i] * tx[i];  // SUM x^2
      }
      u = sum / (sqrt(sumx2) + eps);
      gu = s / (1.0 + exp(-a*(u-b)));

      if (gu > s2max){
	s2max = gu;
	maxi = ii;
	maxj = jj;
      }else if ((ii == midi) && (jj == midi)){
	if (gu == s2max){  // Use the middle unit in case of tie for max
	  //printf("HERE ****************** ****************\n");
	  maxi = ii;
	  maxj = jj;
	}
      }
    }
  }

  if ((maxi == -1)||(maxj == -1))
    exit_error("MOD_VHF_RESP_C2_MAXIJ","Negative index value");

  *ri = maxi;
  *rj = maxj;
  return s2max;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_VHF_RESP_MAX_IJ_QUICK                        */
/*                                                                           */
/*  Determine which S2 unit has the maximum response.                        */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_resp_max_ij_quick(s,a,b,eps,w,x,s2n,n,ri,rj)
     float s,a,b;  // Three sigmoid params
     float eps;    // Normalization constant (to avoid division by zero)
     float *w;     // [n] Weights for C1 inputs 
     float ***x;   // [s2n][s2n][n] Drive of C1 inputs
     int s2n;      // dimension of S2 grid
     int n;        // Number of C1 inputs
     int *ri,*rj;  // Return the S2 coords of the max response
{
  int i,ii,jj;
  int maxi,maxj,midi;
  float sum,sumx2,u,gu,s2max,*tx;

  // *** WYETH NEW - WHY NOT USE precomputed normalization values, ...dj here?
  // *** WYETH NEW - WHY NOT USE precomputed normalization values, ...dj here?
  // *** WYETH NEW - WHY NOT USE precomputed normalization values, ...dj here?

  s2max = -1000.0;  // *** WYETH - MUST BE LOWER THAN ANY POSSIBLE VALUE

  midi = s2n/2;

  maxi = maxj = -1;

  for(ii=0;ii<s2n;ii++){
    for(jj=0;jj<s2n;jj++){

      tx = x[ii][jj];  // Pointer to the C1 inputs for this S2 unit

      //
      //  Eqn 1 (Cadieu et al., 2007)
      //
      sum = 0.0;
      sumx2 = 0.0;
      for(i=0;i<n;i++){  // For each weighted input
	sum += w[i] * tx[i];     // SUM w*x
	sumx2 += tx[i] * tx[i];  // SUM x^2
      }
      if (1){
	u = sum / (sqrt(sumx2) + eps);
      }else{
	if (sum >= 0.0)   // Changes result in 5th decimal place
	  u = sum*sum / sumx2;
	else
	  u = -sum*sum / sumx2;
      }


      //**** WYETH - SHOULD THIS BE a*s ???
      //**** WYETH - SHOULD THIS BE a*s ???
      //**** WYETH - SHOULD THIS BE a*s ???
      //**** WYETH - SHOULD THIS BE a*s ???
      if (a >= 0.0)
	gu = u;
      else
	gu = -u;

      //gu = s / (1.0 + exp(-a*(u-b)));

      //printf("gu = %f   (s,a,b) %f %f %f\n",gu,s,a,b);

      if (gu > s2max){
	s2max = gu;
	maxi = ii;
	maxj = jj;
      }else if ((ii == midi) && (jj == midi)){
	if (gu == s2max){  // Use the middle unit in case of tie for max
	  //printf("HERE ****************** ****************\n");
	  maxi = ii;
	  maxj = jj;
	}
      }
    }
  }

  if ((maxi == -1)||(maxj == -1)){
    mylog(mylogf,"******************  MOD_VHF_RESP_C2_MAXIJ - Negative**\n");

    printf(" (s,a,b) %f %f %f   w[0] = %f\n",s,a,b,w[0]);
    exit_error("MOD_VHF_RESP_C2_MAX_IJ_QUICK","Negative index value");
  }

  *ri = maxi;
  *rj = maxj;
  //return s2max;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_VHF_RESP_C2                              */
/*                                                                           */
/*  Response of C2 unit.                                                     */
/*                                                                           */
/*****************************************************************************/
float mod_vhf_resp_c2(s,a,b,eps,w,x,dj,s2n,n)
     float s,a,b;  // Three sigmoid params
     float eps;    // Normalization constant (to avoid division by zero)
     float *w;     // [n] Weights for C1 inputs 
     float ***x;   // [s2n][s2n][n] Drive of C1 inputs
     float **dj;   // [s2n][s2n] norm value
     int s2n;      // dimension of S2 grid
     int n;        // Number of C1 inputs
{
  int i,ii,jj;
  float sum,sumx2,u,gu,s2max,*tx;

  s2max = -100000.0;  // *** WYETH - MUST BE LOWER THAN ANY POSSIBLE VALUE

  for(ii=0;ii<s2n;ii++){
    for(jj=0;jj<s2n;jj++){

      tx = x[ii][jj];  // Pointer to the C1 inputs for this S2 unit

      //
      //  Eqn 1 (Cadieu et al., 2007)
      //
      sum = 0.0;
      sumx2 = 0.0;
      for(i=0;i<n;i++){  // For each weighted input
	sum += w[i] * tx[i];  // SUM w*x
	//REPLACE?  sumx2 += tx[i] * tx[i];  // SUM x^2
      }
      //REPLACE?  u = sum / (sqrt(sumx2) + eps);
      u = sum / dj[ii][jj];

      gu = s / (1.0 + exp(-a*(u-b)));

      if (gu > s2max)
	s2max = gu;
    }
  }

  return s2max;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_VHF_RESP_C2_FAST                           */
/*                                                                           */
/*  Faster version of the above.                                             */
/*  (1)  Use precomputed 'dj' values for normalization.                      */
/*  (2)  Track the maximal value *before* the nonlinearity, and apply        */
/*       the nonlinearity only once at the end.                              */
/*                                                                           */
/*****************************************************************************/
float mod_vhf_resp_c2_fast(s,a,b,eps,w,x,dj,s2n,n)
     float s,a,b;  // Three sigmoid params
     float eps;    // Normalization constant (to avoid division by zero)
     float *w;     // [n] Weights for C1 inputs 
     float ***x;   // [s2n][s2n][n] Drive of C1 inputs
     float **dj;   // [s2n][s2n] norm value
     int s2n;      // dimension of S2 grid
     int n;        // Number of C1 inputs
{
  int i,ii,jj;
  float sum,u,gu,s2max,*tx;
  //float sumx2;  // *** NOTE, only slight diff. whether we use dj or sumx2

  s2max = -100000.0;  // *** WYETH - MUST BE LOWER THAN ANY POSSIBLE VALUE

  for(ii=0;ii<s2n;ii++){
    for(jj=0;jj<s2n;jj++){
      tx = x[ii][jj];  // Pointer to the C1 inputs for this S2 unit

      //
      //  Eqn 1 (Cadieu et al., 2007)
      //
      sum = 0.0;
      //sumx2 = 0.0;
      for(i=0;i<n;i++){  // For each weighted input
	sum += w[i] * tx[i];  // SUM w*x
	//sumx2 += tx[i] * tx[i];  // SUM x^2
      }
      //u = sum / (sqrt(sumx2) + eps);
      u = sum / dj[ii][jj];

      if (a*s < 0.0)  // If the sigmoid curve is reversed ('a' xor 's' is < 0)
	u = -u;

      if (u > s2max)
	s2max = u;
    }
  }

  if (a*s < 0.0)  // If the sigmoid curve is reversed ('a' xor 's' is < 0)
    s2max = -s2max;

  gu = s / (1.0 + exp(-a*(s2max-b)));

  return gu;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_VHF_FUNC_MIN_01                            */
/*                                                                           */
/*****************************************************************************/
float mod_vhf_func_min_01(x)
     float *x;
{
  int i;
  int nrsp,nc1,s2n;
  float mse,d,s,a,b,eps,*w,*targ,****c1r;
  float **dj;

  nc1  = mod_vhf_fit_nc1;      // Number of C1 inputs for current fitting
  nrsp = mod_vhf_fit_train_n;  // Number of responses to fit
  targ = mod_vhf_fit_train_r;  // [nrsp] Target response vector
  c1r  = mod_vhf_fit_c1r;      // [nrsp][s2i][s2j][nc1] C1 outputs
  s2n = mod_vhf_c2->s2_xn;     // Size of S2 grid

  s = x[0];     // First 3 values in vector are sigmoid params
  a = x[1];
  b = x[2];
  w = &(x[3]);  // remaining (nc1) values are weights for C1 inputs

  eps = mod_vhf_s2->eps2;

  mse =  0.0;
  for(i=0;i<nrsp;i++){  // For each response being fit

    dj = mod_vhf_fit_dj[i];

    /*
    if (     mod_vhf_resp_c2(s,a,b,eps,w,c1r[i],dj,s2n,nc1) != 
	mod_vhf_resp_c2_fast(s,a,b,eps,w,c1r[i],dj,s2n,nc1)){
      printf("a = %f\n",a);
      exit_error("Not equal","not equal");
    }
    */

    //d = targ[i] - mod_vhf_resp_c2(s,a,b,eps,w,c1r[i],dj,s2n,nc1);  // diff
    d = targ[i] - mod_vhf_resp_c2_fast(s,a,b,eps,w,c1r[i],dj,s2n,nc1);

    mse += d*d;
  }
  mse /= nrsp;
  //printf("mse = %f\n",mse);

  return mse;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_VHF_FUNC_MIN_01_D                           */
/*                                                                           */
/*  Return the derivatives of func 01 for each param in 'rd'.                */
/*                                                                           */
/*****************************************************************************/
int mod_vhf_func_min_01_d(x,rd)
     float *x;   // [3 + nc1] Evaluated at this point
     float *rd;  // [3 + nc1] Returned derivatives
{
  int i,j;
  int nrsp,nc1,ii,jj,s2n;
  float s,a,b,*w,*targ,****c1r,***d,*tr,pod,n2,n2s,n2sa,eps,maxv;
  float uj,fj,hj,vj,pj,sum,d_s,d_a,d_b,*d_w;

  //printf("STARTING_____MIN_01_D\n");

  nc1  = mod_vhf_fit_nc1;      // Number of C1 inputs for current fitting
  nrsp = mod_vhf_fit_train_n;  // Number of responses to fit
  targ = mod_vhf_fit_train_r;  // [nrsp] Target response vector
  c1r  = mod_vhf_fit_c1r;      // [nrsp][s2n][s2n][nc1] C1 outputs
  d_w  = mod_vhf_fit_dw;       // Pre-alloc'd array for accumulating derivs.
  d    = mod_vhf_fit_dj;       // [nrsp][s2n][s2n] Norm. const. from C1 resps.

  s2n  = mod_vhf_c2->s2_xn;
  eps  = mod_vhf_s2->eps2;

  s = x[0];     // First 3 values in vector are sigmoid params
  a = x[1];
  b = x[2];
  w = &(x[3]);  // remaining (nc1) values are weights for C1 inputs

  // Zero the accumulation variables
  d_s = 0.0;
  d_a = 0.0;
  d_b = 0.0;
  for(i=0;i<nc1;i++)
    d_w[i] = 0.0;

  // WYETH DEBUG - This was added to find cases where the params are out of
  //   range.  It would be good to specify ranges in the .moo <fit> so that
  //   this could be done systematically
  //
  if ((s > 100.0) || (s < -100.0)){
    printf("  *** MOD_VHF_FUNC_MIN_01_D  Note, 's' is too big, stopping.");
    printf("    s= %f  a= %f  b= %f\n",s,a,b);
    printf("    w[0] = %f\n",w[0]);
    //exit_error("MOD_VHF_FUNC_MIN_01_D","s is too big");
    return -10;  // Indicate that a parameter value was out of range
  }

  for(j=0;j<nrsp;j++){  // For each stimulus

    //
    //  WYETH DERIV MAX - We would have to compute which S2 unit has the
    //    max response for this stimulus, and then use the appropriate C1
    //    responses.  This means keeping 9 times more C1 responses around.
    //
    //    WYETH - we could drop the gradient descent altogether ???
    //
    //  WYETH -- I THINK THE NORM CONSTANTS, E.G. the 'd[..]' array, would
    //     have to be computed for all 9 units as well.
    //

    //  To compute the Max responses, we don't want to have to go through the
    // non-lin, so if we simply establish whether the nonlin is
    // rising or falling (is it allowed to fall?) then we can just use the
    // magnitude of the weights times the C1 inputs.  Do we have to include
    // the normalization?? - Perhaps, because that value can vary across
    // S2 units.

    //maxv = mod_vhf_resp_c2_maxij(s,a,b,eps,w,c1r[j],s2n,nc1,&ii,&jj);

    mod_vhf_resp_max_ij_quick(s,a,b,eps,w,c1r[j],s2n,nc1,&ii,&jj);

    /*
    if ((ii == 1) && (jj == 1))
      printf("**..** Stim %3d  max i,j   %d  %d    max = %f\n",j,ii,jj,maxv);
    else // if ((ii != 1) || (jj != 1))
      printf("****** Stim %3d  max i,j   %d  %d    max = %f\n",j,ii,jj,maxv);
    */

    //ii = jj = 0;

    tr = c1r[j][ii][jj];

    sum = 0.0;  // Compute dot product of weights times C1 responses
    for(i=0;i<nc1;i++)
      sum += w[i] * tr[i];  // WYETH DERIV MAX - 'tr[i]' would vary w/ stimulus
                            //    to reflect which S2 is max?

    //uj = sum/d[j];  // d[j] are constants that do not depend on params
    uj = sum/d[j][ii][jj];  // d[j] are constants that do not depend on params
    fj = exp(-a*(uj - b));
    hj = 1.0 + fj;
    vj = s / hj - targ[j];  // *** INSERT SHUFFLE INDEX ***
    pj = fj * vj/(hj*hj);

    d_s += vj/hj;
    d_a += pj*(uj - b);
    d_b += pj;

    //pod = pj/d[j];
    pod = pj/d[j][ii][jj];
    for(i=0;i<nc1;i++)
      d_w[i] += pod * tr[i];
  }

  n2   = 2.0/(float)nrsp;  // Constants
  n2s  = n2*s;
  n2sa = n2s*a;

  rd[0] =  n2   * d_s;  //  Sigmoid param 's'
  rd[1] =  n2s  * d_a;  //  Sigmoid param 'a'
  rd[2] = -n2sa * d_b;  //  Sigmoid param 'b'

  for(i=0;i<nc1;i++){
    rd[3+i] = n2sa * d_w[i];  // Weights
  }

  //
  //  The returned value should be 1 (success) or >= -10 to indicate some
  //  error.  This value is checked by the caller, e.g., minu_pr().
  //
  return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_VHF_FIT_PREP                             */
/*                                                                           */
/*  Prepare to run the minimization routines for 'nw' C1 weights.            */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_fit_prep(nw)
     int nw;
{
  mylog(mylogf,"  MOD_VHF_FIT_PREP\n");

  sprintf(ggstr,"    nw = %d\n",nw);
  mylog(mylogf,ggstr);

  mod_vhf_fit_nc1 = nw;      // Number of C1 weights
  mod_vhf_fit_xn  = nw + 3;  // Total number of parameters, w's + sigmoid

  if (mod_vhf_fit_dw != NULL)
    myfree(mod_vhf_fit_dw);
  mod_vhf_fit_dw = (float *)myalloc(mod_vhf_fit_nc1*sizeof(float));

  if (mod_vhf_fit_x != NULL)
    myfree(mod_vhf_fit_x);
  mod_vhf_fit_x = (float *)myalloc(mod_vhf_fit_xn*sizeof(float));
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_VHF_FIT_PREP_C1RS                          */
/*                                                                           */
/*  Prepare 'mod_vhf_fit_c1r'.                                               */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_fit_prep_c1rs(i1,i2,wi,sqs)
     int i1;     // C1 index (into cilist), or -1 to indicate 'next prep'
     int i2;     // C1 index (into cilist), or -1
     int wi;     // C1 weight list index, used when 'i2' is -1
     float ***sqs;  // [nstm][s2n][s2n] Accumulated sum square responses
{
  int i,j,k,ii,fi;
  int gi1,gi2,xi1,xi2,yi1,yi2,ci1,ci2,s2n,ks,nstm,n0;
  float *****tc1rs,****t,tr0,tr1;

  s2n  = mod_vhf_c2->s2_xn;
  nstm = mod_vhf_fit_train_n;
  t    = mod_vhf_fit_c1r;

  // -----------------------------------------------------------------------
  // *** WYETH
  //  Why not restucture the C1 response array as follows:
  //  xxx[cilist_i][stim][s2i][s2j]  *** WYETH Figure this out ***
  //        OR
  //  xxx[cilist_i][s2i][s2j][stim]  *** WYETH Figure this out ***
  // -----------------------------------------------------------------------


  if (i1 >= 0){
    xi1 = mod_vhf_c1_ndx[i1][0];
    yi1 = mod_vhf_c1_ndx[i1][1];
    ci1 = mod_vhf_c1_ndx[i1][2];
    gi1 = mod_vhf_c1_ndx[i1][3];
  }
  if (i2 >= 0){
    xi2 = mod_vhf_c1_ndx[i2][0];
    yi2 = mod_vhf_c1_ndx[i2][1];
    ci2 = mod_vhf_c1_ndx[i2][2];
    gi2 = mod_vhf_c1_ndx[i2][3];
  }

  for(i=0;i<s2n;i++){
    for(j=0;j<s2n;j++){
      tc1rs = mod_vhf_c2->s2[i][j]->c1rs;

      if ((i1 >= 0) && (i2 >= 0)){
	//
	//  Prep for 'fit_min2'
	//

	for(k=0;k<nstm;k++){
	  ks = mod_vhf_boot_shuff[k];  // Index to shuffle responses arrays

	  tr0 = tc1rs[ks][gi1][xi1][yi1][ci1];
	  tr1 = tc1rs[ks][gi2][xi2][yi2][ci2];
	  t[k][i][j][0] = tr0;
	  t[k][i][j][1] = tr1;

  //if ((i==2) && (j==0) && (k < 10))
	  //printf("WXO (%d,%d) si %d stim %d    %f %f\n",i,j,k,ks,tr0,tr1);

	  // Compute denominator "d" of "u", the normalized input
	  mod_vhf_fit_dj[k][i][j] = sqrt(tr0*tr0 + tr1*tr1) + KSMALL;
	}

      }else if ((i1 >= 0) && (i2 == -1)){
	//
	//  Prep for 'fit_next'
	//

	for(k=0;k<nstm;k++){
	  ks = mod_vhf_boot_shuff[k];  // Index to shuffle responses arrays
	  tr0 = tc1rs[ks][gi1][xi1][yi1][ci1];  // n0 = n1-1
	  t[k][i][j][wi] = tr0;
	  
	  // Compute denominator "d" of "u", the normalized input
	  mod_vhf_fit_dj[k][i][j] = sqrt(sqs[k][i][j] + tr0*tr0) + KSMALL;
	}
      }else{

	n0 = wi;  // The number of pre-existing C1 units, on which to build

	//
	//  Set responses to all current C1 input units to all stimuli
	//
	for(ii=0;ii<n0;ii++){  // For each C1 input unit already included
	  fi = mod_vhf_fit_ndx[ii];  // Get 'ndx' val for already chosen C1 unit
	  xi1 = mod_vhf_c1_ndx[fi][0];
	  yi1 = mod_vhf_c1_ndx[fi][1];
	  ci1 = mod_vhf_c1_ndx[fi][2];
	  gi1 = mod_vhf_c1_ndx[fi][3];

	  //
	  //  This could be moved outside of this routine, since it is being
	  //  redone each time.
	  //

	  for(k=0;k<nstm;k++){
	    ks = mod_vhf_boot_shuff[k];  // Index to shuffle responses arrays
	    t[k][i][j][ii] = tc1rs[ks][gi1][xi1][yi1][ci1];
	  }
	}

	//
	//  Fill 'sqs' to hold constants to be used in denominator "d" of "u"
	//
	for(k=0;k<nstm;k++){
	  sqs[k][i][j] = 0.0;
	  for(ii=0;ii<n0;ii++){  // For each C1 unit
	    sqs[k][i][j] += t[k][i][j][ii] * t[k][i][j][ii];
	  }
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_VHF_TEST_DERIVS                           */
/*                                                                           */
/*  For testing, compare the deriv function values to the actual function.   */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_test_derivs()
{
  int i;
  int i1,i2,oi1,oi2,xi1,xi2,yi1,yi2,gi,nx,nstm,s2n,dflag;
  float fv,*y,*xx,*deriv,*dd;
  float *x;   // [3 + nc1] Evaluated at this point

  printf("  MOD_VHF_TEST_DERIVS\n");

  mod_vhf_fit_prep(2);  // Prepare for fitting with this many weights

  nstm = mod_vhf_fit_nrsp;
  mod_vhf_fit_train_n = nstm;
  printf("  *** Over-riding ..._train_n for ..._prep_c1rs()\n");

  s2n  = mod_vhf_c2->s2_xn;

  x = mod_vhf_fit_x;

  mod_vhf_fit_c1r = get_4d_farray(nstm,s2n,s2n,2);  // C1 responses


  oi1 = 0;
  oi2 = 0;
  xi1 = 0;
  xi2 = 1;  // Neighbors along x-axis
  yi1 = 0;
  yi2 = 0;
  gi = 1;  // Middle group

  i1 = mod_vhf_c1_ndx_get(gi,xi1,yi1,oi1);
  i2 = mod_vhf_c1_ndx_get(gi,xi2,yi2,oi2);

  mod_vhf_fit_prep_c1rs(i1,i2,-1,NULL);

  // Parameter vector
  x[0] = 1.0;  // 's' of sigmoid
  x[1] = 1.0;  // 'a' of sigmoid
  x[2] = 1.0;  // 'b' of sigmoid
  x[3] = 1.0;  // w0
  x[4] = 1.0;  // w1

  nx = 100;
  xx = (float *)myalloc(nx*sizeof(float));
  y  = (float *)myalloc(nx*sizeof(float));
  dd = (float *)myalloc(nx*sizeof(float));

  deriv = (float *)myalloc(mod_vhf_fit_xn*sizeof(float));

  //
  //  The MSE vs. s is parabolic
  //
  for(i=0;i<nx;i++){
    fv = -4.0 + 10.0*(float)i/(float)(nx-1);  // From -4.0 to 10.0
    x[0] = fv;  // s
    xx[i] = fv;
    y[i] = mod_vhf_func_min_01(x);

    dflag = mod_vhf_func_min_01_d(x,deriv);
    dd[i] = deriv[0];
  }
  append_farray_xy_plot("zzz.test_deriv.pl",xx,y,nx,"MSE_v_s");
  append_farray_xy_plot("zzz.test_deriv.pl",xx,dd,nx,"dds_v_s");


  x[0] = x[1] = x[2] = x[3] = x[4] = 1.0;  // Reset all param values

  //
  //  The MSE vs. a
  //
  for(i=0;i<nx;i++){
    fv = -30.0 + 60.0*(float)i/(float)(nx-1);  // From -10 to 20
    x[1] = fv; // a
    xx[i] = fv;
    y[i] = mod_vhf_func_min_01(x);

    dflag = mod_vhf_func_min_01_d(x,deriv);
    dd[i] = deriv[1];
  }
  append_farray_xy_plot("zzz.test_deriv.pl",xx,y,nx,"MSE_v_a");
  append_farray_xy_plot("zzz.test_deriv.pl",xx,dd,nx,"dds_v_a");


  x[0] = x[1] = x[2] = x[3] = x[4] = 1.0;  // Reset all param values

  //
  //  The MSE vs. b
  //
  for(i=0;i<nx;i++){
    fv = -10.0 + 20.0*(float)i/(float)(nx-1);  // From 0 to 3
    x[2] = fv; // b
    xx[i] = fv;
    y[i] = mod_vhf_func_min_01(x);

    dflag = mod_vhf_func_min_01_d(x,deriv);
    dd[i] = deriv[2];
  }
  append_farray_xy_plot("zzz.test_deriv.pl",xx,y,nx,"MSE_v_b");
  append_farray_xy_plot("zzz.test_deriv.pl",xx,dd,nx,"dds_v_b");


  x[0] = x[1] = x[2] = x[3] = x[4] = 1.0;  // Reset all param values

  //
  //  The MSE vs. w0
  //
  for(i=0;i<nx;i++){
    fv = -20.0 + 40.0*(float)i/(float)(nx-1);  // From -10 to 20
    x[3] = fv;  // w0
    xx[i] = fv;
    y[i] = mod_vhf_func_min_01(x);

    dflag = mod_vhf_func_min_01_d(x,deriv);
    dd[i] = deriv[3];
  }
  append_farray_xy_plot("zzz.test_deriv.pl",xx,y,nx,"MSE_v_w0");
  append_farray_xy_plot("zzz.test_deriv.pl",xx,dd,nx,"dds_v_w0");


  x[0] = x[1] = x[2] = x[3] = x[4] = 1.0;  // Reset all param values

  //
  //  The MSE vs. w1
  //
  for(i=0;i<nx;i++){
    fv = -20.0 + 40.0*(float)i/(float)(nx-1);  // From -10 to 20
    x[4] = fv;  // w1
    xx[i] = fv;
    y[i] = mod_vhf_func_min_01(x);

    dflag = mod_vhf_func_min_01_d(x,deriv);
    dd[i] = deriv[4];
  }
  append_farray_xy_plot("zzz.test_deriv.pl",xx,y,nx,"MSE_v_w1");
  append_farray_xy_plot("zzz.test_deriv.pl",xx,dd,nx,"dds_v_w1");

  myfree(xx);
  myfree(y);
  myfree(dd);
  myfree(deriv);
  free_4d_farray(mod_vhf_fit_c1r,nstm,s2n,s2n,2);
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_VHF_DERIV_TEST_S2                          */
/*                                                                           */
/*  Test the derivative function to see if this S2 config is a minimum.      */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_deriv_test_s2()
{
  int i,j,ii,jj;
  int nw,gi,xi,yi,oi,nx,nstm,s2n,dflag;
  float fv,*y,*xx,****t,*deriv,*dd,ssq,*x,vold;
  float *****tc1rs;
  char tstr[SLEN];

  printf("  MOD_VHF_DERIV_TEST_S2\n");

  s2n = mod_vhf_c2->s2_xn;

  nw = mod_vhf_s2->nw;

  mod_vhf_fit_prep(nw);  // Prepare for fitting with this many weights

  nstm = mod_vhf_fit_nrsp;

  x = mod_vhf_fit_x;

  //
  //  Create response array for C1 inputs to S2 unit
  //
  //t = get_2d_farray(nstm,nw);    // [nstim][nc1] Model responses to stimuli
  t = get_4d_farray(nstm,s2n,s2n,nw); // [nstim]..[nc1] C1 responses
  mod_vhf_fit_c1r = t;

  for(ii=0;ii<s2n;ii++){
    for(jj=0;jj<s2n;jj++){

      tc1rs = mod_vhf_c2->s2[ii][jj]->c1rs;  // [nstim][c1grid][gx][gy][nch]

      for(j=0;j<nw;j++){
	gi = mod_vhf_s2->c1w[j]->gi;
	xi = mod_vhf_s2->c1w[j]->xi;
	yi = mod_vhf_s2->c1w[j]->yi;
	oi = mod_vhf_s2->c1w[j]->ci;
	for(i=0;i<nstm;i++){
	  t[i][ii][jj][j] = tc1rs[i][gi][xi][yi][oi];
	}
      }

      //
      //  Compute denominator "d" of "u", the normalized input
      //
      for(i=0;i<nstm;i++){
	ssq = 0.0;
	for(j=0;j<nw;j++){
	  ssq += t[i][ii][jj][j] * t[i][ii][jj][j];
	}
	mod_vhf_fit_dj[i][ii][jj] = sqrt(ssq) + KSMALL;
      }
    }
  }

  // Parameter vector
  if (strcmp(mod_vhf_s2->nltype,"sigmoid0")!=0)
    exit_error("MOD_VHF_DERIV_TEST_S2","Only works for 'sigmoid0' s2 nonlin");
  x[0] = mod_vhf_s2->s;  // 's' of sigmoid
  x[1] = mod_vhf_s2->a;  // 'a' of sigmoid
  x[2] = mod_vhf_s2->b;  // 'b' of sigmoid
  for(i=0;i<nw;i++)
    x[3+i] = mod_vhf_s2->c1w[i]->w;  // weight


  //
  //  Storage for deriv plots
  //
  nx = 100;
  xx = (float *)myalloc(nx*sizeof(float));
  y  = (float *)myalloc(nx*sizeof(float));
  dd = (float *)myalloc(nx*sizeof(float));

  deriv = (float *)myalloc(mod_vhf_fit_xn*sizeof(float));

  //
  //  The MSE vs. s
  //
  vold = x[0];
  for(i=0;i<nx;i++){
    fv = -5.0 + 10.0*(float)i/(float)(nx-1);  // From -5 to +5
    x[0] = vold + fv;  // s
    xx[i] = fv;
    y[i] = mod_vhf_func_min_01(x);

    //***WYETH - NaN here since boostrap added ?????  Do we need to fix?
    //***WYETH - NaN here since boostrap added ?????  on non-bootstrap .moo's?
    //printf("______y[i] = %f\n",y[i]);

    dflag = mod_vhf_func_min_01_d(x,deriv);
    dd[i] = deriv[0];
  }
  x[0] = vold;  // RESTORE THE VARIED PARAM
  append_farray_xy_plot("zzz.test_deriv.pl",xx,y,nx,"MSE_v_s");
  append_farray_xy_plot("zzz.test_deriv.pl",xx,dd,nx,"dds_v_s");


  //
  //  The MSE vs. a
  //
  vold = x[1];
  for(i=0;i<nx;i++){
    fv = -20.0 + 40.0*(float)i/(float)(nx-1);  // From -20 to 20
    x[1] = vold + fv; // a
    xx[i] = fv;
    y[i] = mod_vhf_func_min_01(x);

    dflag = mod_vhf_func_min_01_d(x,deriv);
    dd[i] = deriv[1];
  }
  x[1] = vold;  // RESTORE THE VARIED PARAM
  append_farray_xy_plot("zzz.test_deriv.pl",xx,y,nx,"MSE_v_a");
  append_farray_xy_plot("zzz.test_deriv.pl",xx,dd,nx,"dds_v_a");


  //
  //  The MSE vs. b
  //
  vold = x[2];
  for(i=0;i<nx;i++){
    fv = -10.0 + 20.0*(float)i/(float)(nx-1);  // From -10 to 10
    x[2] = vold + fv; // b
    xx[i] = fv;
    y[i] = mod_vhf_func_min_01(x);

    dflag = mod_vhf_func_min_01_d(x,deriv);
    dd[i] = deriv[2];
  }
  x[2] = vold;  // RESTORE THE VARIED PARAM
  append_farray_xy_plot("zzz.test_deriv.pl",xx,y,nx,"MSE_v_b");
  append_farray_xy_plot("zzz.test_deriv.pl",xx,dd,nx,"dds_v_b");


  //
  //  The MSE vs. w[j]
  //
  for(j=0;j<nw;j++){
    vold = x[3+j];
    for(i=0;i<nx;i++){
      fv = -10.0 + 20.0*(float)i/(float)(nx-1);  // From -10 to 10
      x[3+j] = vold + fv;  // w0
      xx[i] = fv;
      y[i] = mod_vhf_func_min_01(x);

      dflag = mod_vhf_func_min_01_d(x,deriv);
      dd[i] = deriv[3+j];
    }
    x[3+j] = vold;  // RESTORE THE VARIED PARAM
    sprintf(tstr,"MSE_v_w[%d]",j);
    append_farray_xy_plot("zzz.test_deriv.pl",xx,y,nx,tstr);
    sprintf(tstr,"D_v_w[%d]",j);
    append_farray_xy_plot("zzz.test_deriv.pl",xx,dd,nx,tstr);
  }

  myfree(xx);
  myfree(y);
  myfree(dd);
  myfree(deriv);
  free_4d_farray(t,nstm,s2n,s2n,nw);
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_VHF_S2_RESP_ALL                            */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_s2_resp_all()
{
  int i,j;
  int nw,gi,xi,yi,oi,nx,nstm,s2i,s2j;
  float *x,*w,*r,u;
  float *****tc1rs;
  struct vhf_s2_par *s2;

  // *** CURRENTLY ONLY USED TO TEST RESULTS AFTER FITTING
  // *** MAY NOT BE NEEDED GOING FORWARD

  printf("  MOD_VHF_S2_RESP_ALL\n");

  printf("  *** WARNING:  Using only central S2 unit\n");
  s2i = s2j = mod_vhf_c2->s2_xn / 2;
  tc1rs = mod_vhf_c2->s2[s2i][s2j]->c1rs;   // [nstim][c1grid][gx][gy][nch]


  s2 = mod_vhf_s2;     // Current S2 configuration

  nw = s2->nw;
  x = (float *)myalloc(nw*sizeof(float));
  w = (float *)myalloc(nw*sizeof(float));

  nstm = mod_vhf_fit_nrsp;
  r = (float *)myalloc(nstm*sizeof(float));  // Storage for all responses

  //
  //  Compute S2 response for each stimulus
  //
  for(i=0;i<nstm;i++){  // For each stimulus
    for(j=0;j<nw;j++){
      gi = s2->c1w[j]->gi;
      xi = s2->c1w[j]->xi;
      yi = s2->c1w[j]->yi;
      oi = s2->c1w[j]->ci;
      x[j] = tc1rs[i][gi][xi][yi][oi];  // C1 response for jth input
      w[j] = s2->c1w[j]->w;                     // Weight for jth input
    }
    r[i] = mod_vhf_resp_s2(s2->nltype,s2->s,s2->a,s2->b,s2->eps2,w,x,nw,&u);
  }

  append_farray_xy_plot("zzz.s2_v_targ.pl",mod_vhf_fit_r,r,nstm,"S2_v_targ");

  printf("    Mean neuronal response:  %f\n",mean_farray(mod_vhf_fit_r,nstm));
  printf("    Mean    S2    response:  %f\n",mean_farray(r,nstm));
  printf("\n");

  myfree(r);
  myfree(x);
  myfree(w);
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_VHF_PP_WRITE                             */
/*                                                                           */
/*  Write out a .pp file for the model.                                      */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_pp_write(outfile,fitflag)
     char *outfile;
     int fitflag;
{
  int i,k;
  int n,xi,yi,ci,gi,chi;
  float w,pp_x,pp_y,pp_ori,pp_w,pp_h,pp_lw,plw_in,ppin,xdeg,icon_w,icon_l;
  float faw,wmin,wmax,wrange,lw_min,lw_d;
  char tstr[SLEN3],cstr[SLEN];
  struct vhf_grid *gt;

  lw_min = 0.5;  // Line width for min weigh value (pts)
  lw_d   = 5.0;  // Line width delta to add to get to max weight (pts)

  plw_in = 3.0;  // Plot width (inches) on paper
  ppin = (float)mod_vhf_xn / plw_in;  // Model stim pixels per paper inch
  xdeg = (float)mod_vhf_xn * mod_vhf_sscale; // deg width of full x-axis

  if (fitflag == 1)
    n = mod_vhf_fit_nc1;  // Number of C1 inputs
  else
    n = mod_vhf_s2->nw;   // Number of C1 inputs

  //
  // (0) Determine the min and max of |w| across weights
  //
  for(i=0;i<n;i++){
    if (fitflag == 1)
      w = mod_vhf_fit_opt[3+i];
    else
      w =  mod_vhf_s2->c1w[i]->w;

    faw = fabs(w);
    if ((i==0)||(faw < wmin))
      wmin = faw;
    if ((i==0)||(faw > wmax))
      wmax = faw;
  }
  wrange = wmax - wmin;

  remove_file(outfile);

  //
  // (1) Write header for .pp file
  //
  append_string_to_file(outfile,"beginheader\n");
  append_string_to_file(outfile,"plots 1\n");
  append_string_to_file(outfile,"x 3.0\n");
  append_string_to_file(outfile,"y 6.0\n");
  sprintf(tstr,"width %.2f\n",plw_in);
  append_string_to_file(outfile,tstr);
  sprintf(tstr,"height %.2f\n",plw_in);
  append_string_to_file(outfile,tstr);
  append_string_to_file(outfile,"axis_fontsize 12\n");
  append_string_to_file(outfile,"axis_fontsize_num 10\n");
  sprintf(tstr,"htext 1.0 1.0 14 %s\n",outfile);
  append_string_to_file(outfile,tstr);
  append_string_to_file(outfile,"endheader\n\n");


  //
  // (2) Plot type and ellipses for each C1 input
  //
  append_string_to_file(outfile,"type normal\n");
  for(i=0;i<n;i++){
    if (fitflag == 1){
      k  = mod_vhf_fit_ndx[i];  // Index in list of all candidates
      xi = mod_vhf_c1_ndx[k][0];
      yi = mod_vhf_c1_ndx[k][1];
      ci = mod_vhf_c1_ndx[k][2];  // Channel index
      gi = mod_vhf_c1_ndx[k][3];
      w = mod_vhf_fit_opt[3+i];
    }else{
      xi = mod_vhf_s2->c1w[i]->xi;
      yi = mod_vhf_s2->c1w[i]->yi;
      ci = mod_vhf_s2->c1w[i]->ci;
      gi = mod_vhf_s2->c1w[i]->gi;
      w = mod_vhf_s2->c1w[i]->w;
    }

    gt = mod_vhf_c1_grid[gi];


    //
    //  Write data for postplot
    //
    //   #        x   y  ori   w   h  lw_pts   r   g   b   rgb_fill
    //   ellipse 1.0 1.0  45  0.4 0.2  1.0    0.0 0.0 0.0  -1 -1 -1
    //

    // *** WYETH HERE *** Need to compute these values for postplot
    // *** WYETH HERE *** Need to compute these values for postplot
    // *** WYETH HERE *** Need to compute these values for postplot

    //  mod_vhf_sscale  (deg/pix)
    

    //     lower left s1 unit    width in s1 units
    pp_x = gt->s1x0[0][xi][yi] + gt->s1xn[0]/2;   // (pix)
    pp_y = gt->s1y0[0][xi][yi] + gt->s1yn[0]/2;   // (pix)
    //printf("grid %d  pp_x %f  ppy_y %f\n",gi,pp_x,pp_y);
    // Grid 0: pp_x =  28 52 76 100  (of 128)

    pp_x /= ppin;  // Convert to inches on plot
    pp_y /= ppin;  // Convert to inches on plot


    mod_vhf_chan_index(ci,&chi,&k);

    // WYETH - HACK - this assumes we're averaging over phase ****
    // WYETH - HACK - this assumes we're averaging over phase ****
    // WYETH - HACK - this assumes we're averaging over phase ****
    // WYETH - HACK - this assumes we're averaging over phase ****
    pp_ori = k * 180.0 / mod_vhf_chans[chi]->nori;  // Orientation (deg)

    icon_w = 4.0 * gt->avgtw / mod_vhf_chans[chi]->f_o;  // +-2 * SD_Orth (pix)
    icon_l = 4.0 * gt->avgtw / mod_vhf_chans[chi]->f_p;  // +-2 * SD_Para (pix)

    pp_w = icon_w / ppin;    // Ellipse long axis (inch)
    pp_h = icon_l / ppin;    // Ellipse short axis (inch)


    pp_lw = lw_min + lw_d * (fabs(w) - wmin)/wrange; // Line thickness (pts)

    if (w >= 0.0){
      sprintf(cstr,"0.8 0 0 -1 -1 -1");  // Red for positive weights
    }else{
      sprintf(cstr,"0 0 1.0 -1 -1 -1");  // Green for negative weights
    }

    sprintf(tstr,"ellipse  %.3f %.3f  %3.0f %.3f %.3f  %.1f  %s\n",
	    pp_x,pp_y,pp_ori,pp_w,pp_h,pp_lw,cstr);
    append_string_to_file(outfile,tstr);
  }

  //
  // (3) Write body of plot
  //
  append_string_to_file(outfile,"begin_axis\n");
  sprintf(tstr,"  range 0 %.2f\n",xdeg);
  append_string_to_file(outfile,tstr);
  append_string_to_file(outfile,"  label Visual angle (deg)\n");
  append_string_to_file(outfile,"  div_size 1\n");  // 1 deg division size
  append_string_to_file(outfile,"  div_minor 1\n");
  append_string_to_file(outfile,"end_axis\n");
  append_string_to_file(outfile,"begin_axis vertical\n");
  sprintf(tstr,"  range 0 %.2f\n",xdeg);
  append_string_to_file(outfile,tstr);
  append_string_to_file(outfile,"  label Visual angle (deg)\n");
  append_string_to_file(outfile,"  div_size 1\n");
  append_string_to_file(outfile,"  div_minor 1\n");
  append_string_to_file(outfile,"end_axis\n");
  append_string_to_file(outfile,"data_pair 0\n\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_VHF_FIT_RUN_VAR                            */
/*                                                                           */
/*  Vary the initial parameters over a specified range.                      */
/*                                                                           */
/*****************************************************************************/
float mod_vhf_fit_run_var(vps,xinit,plotflag)
     struct var_par_struct *vps;
     float *xinit;
     int plotflag;
{
  int i,j;
  int niter,n,nv,nst,i_min,xi;
  float tlrnc,v_min,minv,*x,*p;
  char tstr[SLEN];

  tlrnc = 0.00001;     // allowed error on min
  n = mod_vhf_fit_xn;  // Number of parameter dimensions
  x = mod_vhf_fit_x;   // [n] Starting value for all params

  nst = vps->m_n;      // Number of different starting points for minimization
  nv = vps->m_nvar;    // Number of params to vary, <= n

  i_min = -1;          // No minimum saved yet

  for(i=0;i<nst;i++){  // For each start

    // *** WYETH REMOVE
    //printf("HERE WYETH i= 500\n");
    //if (i==0)
    //i=499;

    //
    //  Reset all params to same starting point, then set new vals for varpars
    //
    for(j=0;j<n;j++) // For each parameter dimension
      x[j] = xinit[j];  // Set initial value


    for(j=0;j<nv;j++){ // For each variable param

      xi = atoi(vps->m_name[j]);  // The name must be the index

      // ****** WYETH - BIG HACK - THIS USES ONLY THE LAST PARAMETER
      //   Need to solve ths - so that it doesn't use hard-coded indices
      // ****** WYETH - BIG HACK - THIS USES ONLY THE LAST PARAMETER

      if (xi == 5){
	//printf(" **********changing xi from %d to %d\n",xi,n-1);
	xi = n-1;
      }

      if (xi < n)  // Only use params that are relevant
	x[xi] = atof(vps->m_vval[i][j]); // Set variable value
    }

    //*** WYETH HACK - DON'T START AT SLOPE ZERO ***//
    //*** WYETH HACK - DON'T START AT SLOPE ZERO ***//
    //*** WYETH HACK - DON'T START AT SLOPE ZERO ***//
    if (x[1] == 0.0)
      printf("*** STARTING WITH SLOPE ZERO\n");

    /***
	for(j=0;j<n;j++){ // For each parameter dimension
	printf("...x[%d] = %f\n",j,x[j]);
	}***/

    minu_pr(x,n,tlrnc,mod_vhf_func_min_01,mod_vhf_func_min_01_d,&niter,&minv);


    // *** WYETH HERE changed on Dec 30, 2015, because I think that niter = 0
    //      is OK.
    //if (niter > 0){
    if (niter >= 0){
      if ((i_min == -1) || (minv < v_min)){  // Keep track of the best minimum
	v_min = minv;
	i_min = i;
	//printf("      start %d  new min:  %f\n",i,minv);
      }
      if (plotflag == 1){
	sprintf(tstr,"%d %f\n",i,minv);
	append_string_to_file("zzz.iter.pl",tstr);
      }
    }else if (niter == -2){
      mylog(mylogf,"*** NaN was encountered.  This start is ignored.\n");
      //printf("niter = %d\n",niter);
    }else if (niter == -1){
      ;  // there were too many iterations
      //printf("niter = %d\n",niter);
    }else if (niter == -10){
      mylog(mylogf,"*** MOD_VHF_FIT_RUN_VAR  Param out of bounds, skipped.\n");
      //printf("niter = %d\n",niter);
    }else{ // else there were too many iterations ???
      ; // Often there are 0 iterations, so we won't report this.  I don't
      //     recall why this is at the moment...
      //mylog(mylogf,"*** MOD_VHF_FIT_RUN_VAR  Other condition.\n");
      //printf("OTHER  niter = %d\n",niter);
    } 
  }

  if (i_min == -1)
    exit_error("MOD_VHF_FIT_RUN_VAR","No valid starts");

  //
  //  Re-run at the minimum, so vector 'mod_vhf_fit_x' has param values
  //
  for(j=0;j<n;j++)
    x[j] = xinit[j];

  for(j=0;j<nv;j++){
    xi    = atoi(vps->m_name[j]);  // The name must be the index
    if (xi == 5){
      //********** BIG HACK - see above
      //********** BIG HACK - see above
      //printf(" **********changing xi from %d to %d\n",xi,n-1);
      xi = n-1;
    }
    x[xi] = atof(vps->m_vval[i_min][j]);
  }

  minu_pr(x,n,tlrnc,mod_vhf_func_min_01,mod_vhf_func_min_01_d,&niter,&minv);

  if (niter > 1)
    printf("      Number of iterations:  %d\n",niter);

  return minv;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_VHF_FIT_PERF                             */
/*                                                                           */
/*  Compute performance metrics for the current model fit, for either        */
/*  train or test data.                                                      */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_fit_perf(ntrain,tflag,outfile,si,rr,rmse)
     int ntrain;     // Number of responses in training set
     int tflag;      // 0-training data, 1-testing data
     char *outfile;  // Name of output file, or NULL
     int si;         // Split index
     float *rr;      // Return the r-value
     float *rmse;    // Return the mse value
{
  int i,j,k,ii,jj;
  int i0,xi,yi,oi,gi,ish,nc1,nresp,s2n;
  float sig_s,sig_a,sig_b,eps,mse,totsq;
  float ***x,*w,*r,*rtarg,rval,*****tc1rs,**tdj;
  char plname[SLEN];

  s2n = mod_vhf_c2->s2_xn;
  nc1 = mod_vhf_fit_nc1;

  x   = get_3d_farray(s2n,s2n,nc1);
  tdj = get_2d_farray(s2n,s2n);  // Hold normalization constants
  w   = (float *)myalloc(nc1*sizeof(float));

  if (tflag == 0){
    i0 = 0;
    nresp = ntrain;
  }else{
    i0 = ntrain;
    nresp = mod_vhf_fit_nrsp - ntrain;
  }

  r     = (float *)myalloc(nresp*sizeof(float)); // Model responses
  rtarg = (float *)myalloc(nresp*sizeof(float)); // Neuronal responses

  //
  //  Compute S2 response for each stimulus
  //
  mse = 0.0;
  for(i=0;i<nresp;i++){  // For each stimulus
    ish = mod_vhf_boot_shuff[i+i0];

    //
    //  *** WYETH NOTE - This is similar to, but does not fit neatly into
    //                   'mod_vhf_fit_prep_c1rs'
    //
    for(j=0;j<nc1;j++){
      k = mod_vhf_fit_ndx[j];  // Index in list of all candidates
      xi = mod_vhf_c1_ndx[k][0];
      yi = mod_vhf_c1_ndx[k][1];
      oi = mod_vhf_c1_ndx[k][2];
      gi = mod_vhf_c1_ndx[k][3];

      for(ii=0;ii<s2n;ii++){
	for(jj=0;jj<s2n;jj++){
	  tc1rs = mod_vhf_c2->s2[ii][jj]->c1rs;
	  x[ii][jj][j] = tc1rs[ish][gi][xi][yi][oi]; // C1 resp for jth input
	}
      }

      w[j] = mod_vhf_fit_opt[3+j];             // Weight for jth input
    }

    for(ii=0;ii<s2n;ii++){
      for(jj=0;jj<s2n;jj++){
	totsq = 0.0;
	for(j=0;j<nc1;j++)
	  totsq += x[ii][jj][j] * x[ii][jj][j];
	tdj[ii][jj] = sqrt(totsq) + KSMALL;
      }
    }

    sig_s = mod_vhf_fit_opt[0];
    sig_a = mod_vhf_fit_opt[1];
    sig_b = mod_vhf_fit_opt[2];
    eps = mod_vhf_s2->eps2;

    //r[i] = mod_vhf_resp_c2(sig_s,sig_a,sig_b,eps,w,x,dj,s2n,nc1);
    r[i] = mod_vhf_resp_c2_fast(sig_s,sig_a,sig_b,eps,w,x,tdj,s2n,nc1);

    rtarg[i] = mod_vhf_fit_r[ish];
    mse += (r[i] - rtarg[i]) * (r[i] - rtarg[i]);
  }

  mse /= (float)nresp;

  if (outfile != NULL){
    if (tflag == 0)
      sprintf(plname,"Train_s2_v_targ_%d",si);
    else
      sprintf(plname,"Test_s2_v_targ_%d",si);

    append_farray_xy_plot(outfile,rtarg,r,nresp,plname);
  }

  simple_pearsn(rtarg-1,r-1,nresp,&rval);  // Arrays must start at [1]

  free_3d_farray(x,s2n,s2n,nc1);
  myfree(w);
  myfree(r);
  myfree(rtarg);
  free_2d_farray(tdj,s2n);

  *rr = rval;
  *rmse = mse;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_VHF_FIT_PERF_SAVE                          */
/*                                                                           */
/*  Save fit performance metrics for each current fit.                       */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_fit_perf_save(k)
     int k;  // Index
{
  int i;
  float tr,tmse;

  i = mod_vhf_boot_split_i;

  mod_vhf_fit_perf(mod_vhf_fit_train_n,0,NULL,i,&tr,&tmse);  // 0-train
  mod_vhf_fit_rvl[k] = tr;
  mod_vhf_fit_mse[k] = tmse;

  mod_vhf_fit_perf(mod_vhf_fit_train_n,1,NULL,i,&tr,&tmse);  // 1-test
  mod_vhf_fit_rvl_ts[k] = tr;
  mod_vhf_fit_mse_ts[k] = tmse;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_VHF_VAR_CONFIG                            */
/*                                                                           */
/*  NOTE:  Only fill 'rvps' if '*rvps' is NULL.                              */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_var_config(fito,stage,npar,rxinit,rvps)
     struct onode *fito;              // <fit> object
     int stage;                       // 1,2
     int npar;                        // Number of parameters to be fit
     float **rxinit;                  // Return initial values
     struct var_par_struct **rvps;    // Return var param structure
{
  int i;
  float *xinit;
  struct onode *stgo,*inito,*varo;

  if (stage == 1)
    stgo = onode_child_get_unique(fito,"stage_1");  // Object for stage 1 pars
  else if (stage == 2)
    stgo = onode_child_get_unique(fito,"stage_2");  // Object for stage 2 pars

  inito = onode_child_get_unique(stgo,"initial_values");
  xinit = (float *)myalloc(npar*sizeof(float));

  if (stage == 1){
    if (npar != 5)
      mylogx(mylogf,"MOD_VHF_VAR_CONFIG","Hack - hardcoded to 2 weights");

    xinit[0] = onode_getpar_flt_exit(inito,"s");
    xinit[1] = onode_getpar_flt_exit(inito,"a");
    xinit[2] = onode_getpar_flt_exit(inito,"b");
    xinit[3] = onode_getpar_flt_exit(inito,"w0");
    xinit[4] = onode_getpar_flt_exit(inito,"w1");
  }else if (stage == 2){
    //
    //   *** 'mod_vhf_fit_opt' must be already set
    //
    for(i=0;i<npar-1;i++)  // Three sigmoid plus 'n1' input weights
      xinit[i] = mod_vhf_fit_opt[i];

    xinit[npar-1] = onode_getpar_flt_exit(inito,"wn");  // Val. for new weights

  }else{
    mylogx(mylogf,"MOD_VHF_VAR_CONFIG","Bad 'stage' number");
  }

  if (*rvps == NULL){
    varo = onode_child_get_unique(inito,"var_param");
    if (varo != NULL){
      // Extract the model var params from onodes, build var par records
      *rvps = mod_util_var_par(mylogf,NULL,varo);
    }else
      mylogx(mylogf,"MOD_VHF_VAR_CONFIG","No <var_param> found in <fit>");
  }

  *rxinit = xinit;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_VHF_FIT_SLAVE_MIN2                          */
/*                                                                           */
/*****************************************************************************/
float mod_vhf_fit_slave_min2(i1,i2,vps,xinit)
     int i1,i2;
     struct var_par_struct *vps;  // WYETH - should be global?
     float *xinit;                // WYETH - should be global?
{
  float mse;

  mod_vhf_fit_prep_c1rs(i1,i2,-1,NULL);

  mse = mod_vhf_fit_run_var(vps,xinit,0);

  return mse;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_VHF_FIT_SLAVE_NEXT                          */
/*                                                                           */
/*****************************************************************************/
float mod_vhf_fit_slave_next(i1,vps,xinit,sqs)
     int i1;
     struct var_par_struct *vps;  // WYETH - should be global?
     float *xinit;                // WYETH - should be global?
     float ***sqs;
{
  int n0;
  float mse;

  n0 = mod_vhf_fit_nc1 - 1;    // Number of C1 inputs in base

  mod_vhf_fit_prep_c1rs(i1,-1,n0,sqs);

  mse = mod_vhf_fit_run_var(vps,xinit,0);

  return mse;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_VHF_FIT_SLAVE_TOP                          */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_fit_slave_top(fito)
     struct onode *fito;
{
  int i,j;
  int tn,nprep,*clist,nci,s2n;
  float *xinit1,*xinit2,***sqs,mse;
  char *cstr;
  struct var_par_struct *vps1,*vps2;

  /* Variables to be set, and the (#) showing where they get are set below:

     mod_vhf_fit_train_n   Already done by:  'mod_vhf_fit_boot()'
     mod_vhf_boot_shuff    (1)
     mod_vhf_fit_train_r   (1.1)
     mod_vhf_fit_xn        (2)
     mod_vhf_fit_nc1       (2)
     mod_vhf_fit_dw        (2) storage allocated
     mod_vhf_fit_x         (2) storage allocated - Has final param vals.
     mod_vhf_fit_c1r       (0) storage allocated
     vps                   (3)
     xinit                 (3)
     mod_vhf_fit_dj        (4.1)

     NOTE:  the above list includes what is needed by 'mod_vhf_func_min_01'
            and 'mod_vhf_func_min_01_d' 
  */

  xinit1 = NULL;
  xinit2 = NULL;
  vps1 = NULL;
  vps2 = NULL;

  s2n = mod_vhf_c2->s2_xn;
  
  //
  //  (0) Allocate storage
  //
  clist = (int *)myalloc(mod_vhf_c1_ntot*sizeof(int));
  sqs = get_3d_farray(mod_vhf_fit_train_n,s2n,s2n);

  mod_vhf_fit_c1r = get_4d_farray(mod_vhf_fit_train_n,s2n,s2n,
				  mod_vhf_fit_c1_max_n);


  while(1){
    mylog(mylogf,"  Waiting to receive a command ...\n");

    cstr = mm_cmd_recv(0,mylogf);
    //printf("SLAVE  Got command  %s\n",cstr);

    if (strcmp(cstr,"get_shuffle_index")==0){
      //
      //  (1) Prepare for the current data partition
      //
      mm_data_recv_iarray(0,mylogf,mod_vhf_boot_shuff,mod_vhf_fit_nrsp,&tn);

      if (tn != mod_vhf_fit_nrsp)
	mylogx(mylogf,"MOD_VHF_FIT_SLAVE_TOP","Bad shuffle index length");

      for(i=0;i<mod_vhf_fit_train_n;i++){  // (1.1)
	j = mod_vhf_boot_shuff[i];
	mod_vhf_fit_train_r[i] = mod_vhf_fit_r[j];
      }

    }else if (strcmp(cstr,"prep_n")==0){
      //
      //  (2) Prepare for fitting with 'n' C1 weights
      //
      mm_data_recv_iarray(0,mylogf,&nprep,1,&tn);
      mod_vhf_fit_prep(nprep);
      
    }else if (strcmp(cstr,"prep_stage_1")==0){
      //
      //  (3) Prepare var params and initial values
      //
      if (xinit1 != NULL)
	myfree(xinit1);

      if (vps1 != NULL)
	mylogx(mylogf,"MOD_VHF_FIT_SLAVE_TOP","Shouldn't prep vps1 twice");

      mod_vhf_var_config(fito,1,mod_vhf_fit_xn,&xinit1,&vps1);

    }else if (strcmp(cstr,"prep_stage_2")==0){
      //
      //  (*) Prepare var params and initial values
      //
      if (xinit2 != NULL)
	myfree(xinit2);

      mod_vhf_var_config(fito,2,mod_vhf_fit_xn,&xinit2,&vps2);

    }else if (strcmp(cstr,"get_c1_list")==0){
      //
      //  (4) Get the list of indices for the C1 units to fit
      //
      mm_data_recv_iarray(0,mylogf,clist,mod_vhf_c1_ntot,&nci);
      //printf("SLAVE  Got C1 index lists\n");

      //
      //  Run minimization
      //
      if (nci == 2){

	if (nci != nprep){
	  printf(" nci = %d  nprep = %d\n",nci,nprep);
	  mylogx(mylogf,"MOD_VHF_FIT_SLAVE_TOP","Array length mismatch");
	}

	mse = mod_vhf_fit_slave_min2(clist[0],clist[1],vps1,xinit1);  // (4.1)
      }else if (nci == 1){
	mse = mod_vhf_fit_slave_next(clist[0],vps2,xinit2,sqs);  // (4.2)
      }else
	mylogx(mylogf,"MOD_VHF_FIT_SLAVE_TOP","Only works for 2 weights");
	
      //printf("(SLAVE) mse = %f\n",mse);

      //
      //  Send the mse back
      //
      mm_data_send_farray(0,1,&mse,1,mylogf);
      //printf("(SLAVE) sent mse  = %f\n",mse);

    }else if (strcmp(cstr,"prep_c1_base")==0){
      //
      //  Receive list of C1 indices for inputs added so far
      //
      mm_data_recv_iarray(0,mylogf,mod_vhf_fit_ndx,mod_vhf_c1_ntot,&nci);
      //printf("SLAVE  Got C1 base list, n= %d\n",nci);

      if ((nci + 1) != mod_vhf_fit_nc1)
	mylogx(mylogf,"MOD_VHF_FIT_SLAVE_TOP","Base n+1 not prepared");

      mm_data_recv_farray(0,mylogf,mod_vhf_fit_opt,MAX_PAR,&tn,NULL);

      // *** WYETH HERE *** print out the optimal params, match to 'wm'
      /*
      {
	int i;
	for(i=0;i<tn;i++)
	  printf(" opt_par[%d] = %f\n",i,mod_vhf_fit_opt[i]);
      }
      */
      //printf("tn = %d   ..._xn = %d\n",tn,mod_vhf_fit_xn);

      if ((tn+1) != mod_vhf_fit_xn)
	mylogx(mylogf,"MOD_VHF_FIT_SLAVE_TOP","Total param mismatch");


      //  WYETH - Sequence of calls should be:
      //             'prep_n', this, 'prep_stage_2', then 'get_c1_list'
      //
      //  WYETH - Now go prepare the 'sqs' array for these
      //          and most of the 'mod_vhf_fit_c1r' array

      mod_vhf_fit_prep_c1rs(-1,-1,nci,sqs);


    }else if (strcmp(cstr,"get_par_val")==0){
      //
      //  (5) Master requests that we send back the parameter values.
      //
      mm_data_send_farray(0,1,mod_vhf_fit_x,mod_vhf_fit_xn,mylogf);
      //printf("(SLAVE) sent par vals\n");

    }else if (strcmp(cstr,"exit")==0){
      exit(0);  // Normal exit
    }else{
      sprintf(ggstr,"  *** command string:  %s\n",cstr);
      mylog(mylogf,ggstr);
      mylogx(mylogf,"MOD_VHF_FIT_SLAVE_TOP","Unknown command");
    }

    myfree(cstr);
  }

  //myfree(sqs);  // This line is never reached, because of exit() above
  free_3d_farray(sqs,mod_vhf_fit_train_n,s2n,s2n);
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_VHF_FIT_MASTER_ASSIGN                         */
/*                                                                           */
/*  Assign the next job to process 'pid'.                                    */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_fit_master_assign(mmfit,pid)
     struct vhf_mmfit *mmfit;
     int pid;
{
  int k;

  k = mmfit->nxti; // Next job number to be assigned
  if (k < mmfit->n){
    mm_cmd_send(pid,1,"get_c1_list",mylogf);
    mm_data_send_iarray(pid,1,mmfit->ci[k],mmfit->nci,mylogf);
    mmfit->stat[k] = 1;     // 1 indicates that job is now assigned
    mmfit->pstat[pid] = k;  // Record job index assigned to proc 'i'
    mmfit->nxti += 1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_VHF_FIT_MASTER_RUN                          */
/*                                                                           */
/*  Run all the jobs loaded into the 'mmfit' list, and return the index      */
/*  of the best, and the MSE value.                                          */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_fit_master_run(mmfit,ri,rmse)
     struct vhf_mmfit *mmfit;     // Information about jobs to be run
     int *ri;                     // Index of best configuration
     float *rmse;                 // MSE minimum value
{
  int i;
  int min_i,done,tn,pid;
  float mse,mse_min;

  //  Make initial assignments to processors
  for(i=1;i<numproc;i++)  // 0 is master, so start at 1
    mod_vhf_fit_master_assign(mmfit,i);  // Does nothing if no more jobs

  min_i = -1;
  done = 0;
  while(done == 0){

    mm_data_recv_farray(-1,mylogf,&mse,1,&tn,&pid);  // Wait for slave result

    if ((min_i == -1) || (mse < mse_min)){ // Is this best MSE so far
      mse_min = mse;  
      min_i = mmfit->pstat[pid];  // Config index that achieved the min

      if (mmfit->nci == 1){
	printf("    new min at (%d)  %f\n",mmfit->ci[min_i][0],mse);
      }else if (mmfit->nci == 2)
	printf("    new min at (%d, %d)  %f\n",mmfit->ci[min_i][0],
	       mmfit->ci[min_i][1],mse);

      mm_cmd_send(pid,1,"get_par_val",mylogf); // Ask slave to send par vals

      // Receive the param values
      mm_data_recv_farray(pid,mylogf,mod_vhf_fit_x,mod_vhf_fit_xn,&tn,NULL);

      if (tn != mod_vhf_fit_xn)
	mylogx(mylogf,"MOD_VHF_FIT_MASTER_RUN","Param count error");

      //printf("MASTER BEST MSE:  %f\n",mse_min);
    }

    mmfit->stat[mmfit->pstat[pid]] = 2;   // Job status: completed

    mmfit->ndone += 1;  // Track total number of jobs finsed so far

    mod_vhf_fit_master_assign(mmfit,pid);  // Assign next job (if more to do)

    if (mmfit->ndone == mmfit->n)
      done = 1;
  }

  //printf("MASTER__________returning mse:  %f\n",mse_min);

  *ri = min_i;
  *rmse = mse_min;
}
/**************************************-**************************************/
/*                                                                           */
/*                        MOD_VHF_FIT_MASTER_MMFIT_FREE                      */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_fit_master_mmfit_free(mmfit)
     struct vhf_mmfit *mmfit;
{
  free_2d_iarray(mmfit->ci,mmfit->n);
  myfree(mmfit->stat);
  myfree(mmfit->mse);
  myfree(mmfit->pstat);

  myfree(mmfit);
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_VHF_FIT_MASTER_MIN2                          */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_fit_master_min2(cil0,ciln,ri1,ri2,rmse)
     int cil0;   // Start index in list of C1 units
     int ciln;   // Number of C1 units
     int *ri1;   // Index of 1st C1 unit
     int *ri2;   // Index of 2nd C1 unit
     float *rmse;   // MSE value
{
  int i,j,k;
  int min_i;
  float mse_min;
  struct vhf_mmfit *mmfit;

  //
  //  (1) Tell slaves to prepare for fitting 2 weights
  //
  mm_cmd_send(1,numproc-1,"prep_n",mylogf);
  i = 2;
  mm_data_send_iarray(1,numproc-1,&i,1,mylogf);

  //
  //  (2) Tell slaves to prepare for stage_1 parameter variation
  //
  if (mod_vhf_boot_split_i == 0)
    mm_cmd_send(1,numproc-1,"prep_stage_1",mylogf);

  //
  //  (3) Build a runlist with all pairs to run
  //
  mmfit = (struct vhf_mmfit *)myalloc(sizeof(struct vhf_mmfit));
  mmfit->n       = (ciln * ciln - ciln)/2;
  mmfit->ndone   = 0;
  mmfit->nci     = 2;
  mmfit->ci      = get_2d_iarray(mmfit->n,mmfit->nci);
  mmfit->stat    = get_zero_iarray(mmfit->n);
  mmfit->mse     = (float *)myalloc(mmfit->n*sizeof(float));
  mmfit->mse_min = -1.0;
  mmfit->nxti    = 0;
  mmfit->pstat   = get_const_iarray(numproc,-1);

  k = 0;
  for(i=0;i<(ciln-1);i++){
    for(j=i+1;j<ciln;j++){
      mmfit->ci[k][0] = i + cil0;
      mmfit->ci[k][1] = j + cil0;
      k += 1;
    }
  }
  if (k != mmfit->n)
    exit_error("MOD_VHF_FIT_MASTER_MIN2","Counting error");

  //
  //  (4) Run all the jobs and find the configuration w/ minimum MSE
  //
  mod_vhf_fit_master_run(mmfit,&min_i,&mse_min);

  *ri1 = mmfit->ci[min_i][0];
  *ri2 = mmfit->ci[min_i][1];
  *rmse = mse_min;

  mod_vhf_fit_master_mmfit_free(mmfit);
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_VHF_FIT_MIN2                             */
/*                                                                           */
/*  Find the two C1 filters that give the minimum MSE.                       */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_fit_min2(fito)
     struct onode *fito;
{
  int i,j,k;
  int nstm,gi1,gi2,gi,min_i,min_j,s2n,**cilist,cil0,ciln;
  float **mse,minval,*minpar,*xinit;
  struct var_par_struct *vps;
  struct onode *stgo;

  printf("  MOD_VHF_FIT_MIN2\n");

  stgo = onode_child_get_unique(fito,"stage_1");  // Object for stage 1 pars
  gi = onode_getpar_int_dflt(stgo,"grid_index",1); // WYETH *** change to exit

  //
  //  Determine index of C1 grid to use, and the range of global C1 indices
  //
  if (gi >= 0)
    mod_vhf_ndx_for_grid(gi,&cil0,&ciln);  // start and length in 'cilist'
  else{
    // use all units
    cil0 = 0;
    ciln = mod_vhf_c1_ntot;
  }

  //  Allocate storage and set:
  //     mod_vhf_fit_nc1 
  //     mod_vhf_fit_xn
  mod_vhf_fit_prep(2);  // Prepare to fit w/ two weights

  minpar = (float *)myalloc(mod_vhf_fit_xn*sizeof(float));


  //
  //  --- --- Preparation done above here may be needed for MM work --- ---
  //


  if (myid == 0){
    //
    //  Mutli-processor (mm) mode
    //
    mod_vhf_fit_master_min2(cil0,ciln,&min_i,&min_j,&minval);

    // fill minpar
    for(k=0;k<5;k++)
      minpar[k] = mod_vhf_fit_x[k];

  }else{
    //
    //  Single-processor (wm) mode
    //
    nstm = mod_vhf_fit_train_n;
    s2n = mod_vhf_c2->s2_xn;
    mod_vhf_fit_c1r = get_4d_farray(nstm,s2n,s2n,2);  // C1 responses

    // WYETH - We don't actually use this array once it is filled
    mse = get_const_2d_farray(ciln,ciln,-1.0);  // Init all MSE values to -1

    vps = NULL;  // Next line will fill 'vps' only if it is NULL
    mod_vhf_var_config(fito,1,mod_vhf_fit_xn,&xinit,&vps);

    //  Compute the MSE for each unique pair
    for(i=0;i<(ciln-1);i++){
      for(j=i+1;j<ciln;j++){

	mod_vhf_fit_prep_c1rs(i+cil0,j+cil0,-1,NULL);

	mse[i][j] = mod_vhf_fit_run_var(vps,xinit,0); // plotflag = 0

	//  Track the minimum value
	if (((i==0)&&(j==1)) || (mse[i][j] < minval)){
	  min_i = i + cil0;
	  min_j = j + cil0;
	  minval = mse[i][j];
	  for(k=0;k<5;k++)
	    minpar[k] = mod_vhf_fit_x[k];
	  printf("    new min at (%d, %d)  %f\n",min_i,min_j,minval);
	}
      }
    }
    free_2d_farray(mse,ciln);
  }

  cilist = mod_vhf_c1_ndx;  // [cn][xi,yi,ori,grid] Indices for all C1 units
  
  gi1 = cilist[min_i][3];  // Grid for 1st C1 unit
  gi2 = cilist[min_j][3];  // Grid for 2nd C1 unit

  sprintf(ggstr,"    Minimum MSE found for pair:  %d %d  in grid %d %d\n",
	  min_i,min_j,gi1,gi2);
  mylog(mylogf,ggstr);
  sprintf(ggstr,"     %3d:  xi, yi, oi =  %d %d %d\n",min_i,cilist[min_i][0],
	  cilist[min_i][1],cilist[min_i][2]);
  mylog(mylogf,ggstr);
  sprintf(ggstr,"     %3d:  xi, yi, oi =  %d %d %d\n",min_j,cilist[min_j][0],
	  cilist[min_j][1],cilist[min_j][2]);
  mylog(mylogf,ggstr);
  for(k=0;k<5;k++){
    sprintf(ggstr,"     param[%d] = %f\n",k,minpar[k]);
    mylog(mylogf,ggstr);
  }


  //
  //  Store the optimal params to be used as initial values in the next step
  //
  for(i=0;i<5;i++)  // Three sigmoid plus two input weights
    mod_vhf_fit_opt[i] = minpar[i];


  //
  //  Store the units that were picked.
  //
  mod_vhf_c1_flag[min_i] = 1;
  mod_vhf_c1_flag[min_j] = 1;
  mod_vhf_fit_ndx[0] = min_i;
  mod_vhf_fit_ndx[1] = min_j;

  // *** WYETH REMOVE THESE
  //mod_vhf_fit_mse[0] = minval;  // Store best MSE after this stage
  //mod_vhf_fit_mse[1] = minval;


  //
  //  Save performance metrics
  //
  mod_vhf_fit_perf_save(0);  // This sets [0], which is used to fill [1] below
  mod_vhf_fit_mse[1]    = mod_vhf_fit_mse[0];
  mod_vhf_fit_rvl[1]    = mod_vhf_fit_rvl[0];
  mod_vhf_fit_mse_ts[1] = mod_vhf_fit_mse_ts[0];
  mod_vhf_fit_rvl_ts[1] = mod_vhf_fit_rvl_ts[0];
  if (mod_vhf_fit_mse[0] != minval)
    exit_error("MOD_VHF_FIT_MIN2","Inconsistent 'minval' should not happen.");


  if (mod_vhf_fit_c1r != NULL)
    free_4d_farray(mod_vhf_fit_c1r,nstm,s2n,s2n,2);
  myfree(minpar);
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_VHF_FIT_MASTER_NEXT                          */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_fit_master_next(nw,ri,rmse)
     int nw;        // Total number of weights, including one to be fit
     int *ri;       // Index of C1 unit to add
     float *rmse;   // MSE value
{
  int i,j,k;
  int min_i;
  float mse_min;
  struct vhf_mmfit *mmfit;

  //
  //  (1) Tell slaves to prepare for fitting 'nw' weights
  //
  mm_cmd_send(1,numproc-1,"prep_n",mylogf);
  mm_data_send_iarray(1,numproc-1,&nw,1,mylogf);

  //
  //  (2) Send base configuration
  //
  mm_cmd_send(1,numproc-1,"prep_c1_base",mylogf);
  // Send list of C1 indices
  mm_data_send_iarray(1,numproc-1,mod_vhf_fit_ndx,nw-1,mylogf);
  // Send optimal params
  mm_data_send_farray(1,numproc-1,mod_vhf_fit_opt,mod_vhf_fit_xn-1,mylogf);


  //
  //  (2) Tell slaves to prepare for stage_2 parameter variation
  //
  mm_cmd_send(1,numproc-1,"prep_stage_2",mylogf);

  //
  //  (*) Count number of C1 units to be tested
  //
  k = 0;
  for(i=0;i<mod_vhf_c1_ntot;i++){
    if (mod_vhf_c1_flag[i] == 0){  // Hasn't been picked yet
      k += 1;
    }
  }

  //
  //  (3) Build a runlist with units to test
  //
  mmfit = (struct vhf_mmfit *)myalloc(sizeof(struct vhf_mmfit));
  mmfit->n       = k;
  mmfit->ndone   = 0;
  mmfit->nci     = 1;
  mmfit->ci      = get_2d_iarray(mmfit->n,mmfit->nci);
  mmfit->stat    = get_zero_iarray(mmfit->n);
  mmfit->mse     = (float *)myalloc(mmfit->n*sizeof(float));
  mmfit->mse_min = -1.0;
  mmfit->nxti    = 0;
  mmfit->pstat   = get_const_iarray(numproc,-1);

  k = 0;
  for(i=0;i<mod_vhf_c1_ntot;i++){
    if (mod_vhf_c1_flag[i] == 0){  // Hasn't been picked yet
      mmfit->ci[k][0] = i;
      k += 1;
    }
  }

  //
  //  (4) Run all the jobs and find the configuration w/ minimum MSE
  //
  mod_vhf_fit_master_run(mmfit,&min_i,&mse_min);

  *ri = mmfit->ci[min_i][0];
  *rmse = mse_min;

  mod_vhf_fit_master_mmfit_free(mmfit);
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_VHF_FIT_NEXT                             */
/*                                                                           */
/*  Find the two C1 filters that give the minimum MSE.                       */
/*                                                                           */
/*****************************************************************************/
int mod_vhf_fit_next(fito)
     struct onode *fito;
{
  int i,k;
  int nstm,cn,min_i,n0,n1,npar,s2n;
  float *mse,minval,*minpar,*xinit,***sqs;
  struct var_par_struct *vps;

  printf("  MOD_VHF_FIT_NEXT\n");

  s2n = mod_vhf_c2->s2_xn;

  if (mod_vhf_fit_nc1 >= mod_vhf_fit_c1_max_n)  // Stop adding new C1 inputs
    return 1;  // Signal to stop fitting new inputs

  cn = mod_vhf_c1_ntot;
  // WYETH - Is NC1  meant to keep this state?
  nstm = mod_vhf_fit_train_n;
  n0 = mod_vhf_fit_nc1;    // Number of C1 input units so far
  n1 = n0 + 1;             // New total number being fit

  //  Allocate storage and set:
  //     mod_vhf_fit_nc1
  //     mod_vhf_fit_xn
  mod_vhf_fit_prep(n1);  // Prepare to fit w/ 'n1' weights

  npar = mod_vhf_fit_xn;  // The total number of params being fit


  //  To store param values at current minimum
  minpar = (float *)myalloc(npar*sizeof(float));

  if (myid == 0){
    //
    //  Mutli-processor (mm) mode
    //
    mod_vhf_fit_master_next(n1,&min_i,&minval);

    // fill minpar
    for(k=0;k<npar;k++)
      minpar[k] = mod_vhf_fit_x[k];


  }else{
    //
    //  Single-processor (wm) mode
    //
    vps = NULL; // Next line will fill 'vps' only if it is NULL
    mod_vhf_var_config(fito,2,npar,&xinit,&vps);

    sqs = get_3d_farray(nstm,s2n,s2n);
    mod_vhf_fit_c1r = get_4d_farray(nstm,s2n,s2n,n1);

    //  Set 'mod_vhf_fit_c1r' and 'sqs'
    mod_vhf_fit_prep_c1rs(-1,-1,n0,sqs);

    //
    //  Loop over all C1 candidates for the 'n1'th input
    //
    mse = get_const_farray(cn,-1.0);  // Initialize all MSE values to -1
    minval = -1.0;
    for(i=0;i<cn;i++){
      if (mod_vhf_c1_flag[i] == 0){  // Hasn't been picked yet

	mod_vhf_fit_prep_c1rs(i,-1,n0,sqs);

	mse[i] = mod_vhf_fit_run_var(vps,xinit,0); // plotflag = 0

	//  Track the minimum value
	if ((minval < 0.0) || (mse[i] < minval)){
	  min_i = i;
	  minval = mse[i];
	  for(k=0;k<npar;k++)
	    minpar[k] = mod_vhf_fit_x[k];
	  printf("    new min for unit %d,  %f\n",i,minval);
	}
      }
    }
    myfree(mse); // This is essentially unused

    mod_util_var_par_free(vps);
    free_4d_farray(mod_vhf_fit_c1r,nstm,s2n,s2n,n1);
    free_3d_farray(sqs,nstm,s2n,s2n);
    myfree(xinit);
  }

  //
  //  Store the optimal params to be used as initial values in the next step
  //
  for(i=0;i<npar;i++)  // Three sigmoid plus two input weights
    mod_vhf_fit_opt[i] = minpar[i];


  mod_vhf_c1_flag[min_i] = 1;   // This one has just been picked
  mod_vhf_fit_ndx[n0] = min_i;  // Add it to the picked 'ndx' list

  //mod_vhf_fit_mse[n0] = minval; // The best MSE value achieved at this point
  mod_vhf_fit_perf_save(n0);


  printf("    Minimum MSE found for unit %d (of %d)\n",min_i,cn);

  printf("     %3d:  xi, yi, oi =  %d %d %d\n",min_i,mod_vhf_c1_ndx[min_i][0],
	 mod_vhf_c1_ndx[min_i][1],mod_vhf_c1_ndx[min_i][2]);
  for(k=0;k<npar;k++)
    printf("     param[%d] = %f\n",k,minpar[k]);

  myfree(minpar);

  return 0;  // Signal to keep fitting new inputs
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_VHF_FIT_WRITE_PARAM                         */
/*                                                                           */
/*  Write a text file containing the parameters for the best fit.            */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_fit_write_param(outfile,si)
     char *outfile;   // output file name
     int si;          // Split index
{
  int i,k;
  int nc1,xi,yi,oi,gi;
  float w,mse;
  char tstr[SLEN];
  struct vhf_grid *gt;

  nc1 = mod_vhf_fit_nc1;

  printf("    Done fitting split %d\n",si);
  printf("      %d C1 units have been chosen.\n",nc1);
  printf("      Writing weights and MSE values to:  %s\n",outfile);
  printf("  -------------------------------------------------------------\n");

  sprintf(tstr,"==================== SPLIT %d ====================\n",si);
  append_string_to_file(outfile,tstr);
  sprintf(tstr,"  sigmoid_s      %7.4f\n",mod_vhf_fit_opt[0]);
  append_string_to_file(outfile,tstr);
  sprintf(tstr,"  sigmoid_alpha  %7.4f\n",mod_vhf_fit_opt[1]);
  append_string_to_file(outfile,tstr);
  sprintf(tstr,"  sigmoid_beta   %7.4f\n",mod_vhf_fit_opt[2]);
  append_string_to_file(outfile,tstr);

  for(i=0;i<nc1;i++){
    k = mod_vhf_fit_ndx[i];  // Index in list of all cadidates
    xi = mod_vhf_c1_ndx[k][0];
    yi = mod_vhf_c1_ndx[k][1];
    oi = mod_vhf_c1_ndx[k][2];
    gi = mod_vhf_c1_ndx[k][3];

    gt = mod_vhf_c1_grid[gi];

    w = mod_vhf_fit_opt[3+i];

    mse = mod_vhf_fit_mse[i];  // Index in list of all cadidates

    sprintf(tstr,"  %s_%d_%d_%d    %7.4f   # MSE_%d: %6.4f\n",gt->name,xi,yi,oi,
	    w,i,mse);
    append_string_to_file(outfile,tstr);
  }
  append_string_to_file(outfile,"\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_VHF_FIT_SPLIT_OUT                           */
/*                                                                           */
/*  Write output following the fit to a particular split.                    */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_fit_split_out(fito,ntrain,ntest,spli,rtr,mtr,rts,mts)
     struct onode *fito;  // <fitp> onode
     int ntrain;          // Number of training points
     int ntest;           // Number of test points
     int spli;            // Split index
     float *rtr;          // Return TRAIN r-value
     float *mtr;          // Return TRAIN MSE
     float *rts;          // Return TEST  r-value
     float *mts;          // Return TEST  MSE
{
  int nc1;
  char tname[SLEN],t1[SLEN],t2[SLEN],tstr[SLEN];

  //
  //  Compare performance
  //
  mod_vhf_fname(fito,"outfile_boot",".boot.rval",t1); // Build outfile name
  mod_vhf_fname(fito,"outfile_boot",".boot.mse" ,t2); // Build outfile name
  mod_vhf_fname(fito,"outfile_s2vtarg",".s2vtarg.pl",tname);

  mod_vhf_fit_perf(ntrain,0,tname,spli,rtr,mtr);  // 0-train
  printf("    TRAIN  r %f  mse %f (n=%d)\n",*rtr,*mtr,ntrain);
  if (ntest > 0){
    mod_vhf_fit_perf(ntrain,1,tname,spli,rts,mts);  // 1-test
    printf("    TEST   r %f  mse %f (n=%d)\n",*rts,*mts,ntest);
  }else{
    *rts = 0.0;  // There was no test set, so simply set these to zero.
    *mts = 0.0;
  }
  //printf(" ==>  ptrain (mse) = %f\n",mod_vhf_fit_mse[mod_vhf_fit_nc1-1]);

  if (spli == 0){
    append_string_to_file(t1,"  TRAIN    TEST\n");
    append_string_to_file(t2,"  TRAIN    TEST\n");
  }
  sprintf(tstr,"%f %f\n",*rtr,*rts);
  append_string_to_file(t1,tstr);
  sprintf(tstr,"%f %f\n",*mtr,*mts);
  append_string_to_file(t2,tstr);


  //
  //  Append fit params to outfile
  //
  mod_vhf_fname(fito,"outfile_best",".best.fit",tname); // Build outfile name
  if (tname[0] == '\0')
    exit_error("MOD_VHF_FIT_SPLIT_OUT","NULL filename for 'outfile_best'");
  if (spli == 0)
    remove_file(tname);
  mod_vhf_fit_write_param(tname,spli);


  //
  //  Write out a .pp file to make a PostPlot map of the C1 RFs
  //
  sprintf(tstr,".%d.pp",spli);
  mod_vhf_fname(fito,"outfile_pp",tstr,tname); // Build outfile name
  if (tname[0] != '\0')
    mod_vhf_pp_write(tname,1);


  //
  //  Plot the MSE value vs. C1 count
  //
  mod_vhf_fname(fito,"outfile_mse",".mse",tname); // Build outfile name
  if (tname[0] != '\0'){
    nc1 = mod_vhf_fit_nc1;
    if (spli == 0)
      remove_file(tname);
    sprintf(tstr,"MSE_vs_C1cnt_Split_%d.pp",spli);
    append_farray_plot(tname,tstr,mod_vhf_fit_mse,nc1,1);

    sprintf(tstr,"rVal_vs_C1cnt_Split_%d.pp",spli);
    append_farray_plot(tname,tstr,mod_vhf_fit_rvl,nc1,1);

    if (ntest > 0){
      sprintf(tstr,"Test__MSE_vs_C1cnt_Split_%d.pp",spli);
      append_farray_plot(tname,tstr,mod_vhf_fit_mse_ts,nc1,1);

      sprintf(tstr,"Test__rVal_vs_C1cnt_Split_%d.pp",spli);
      append_farray_plot(tname,tstr,mod_vhf_fit_rvl_ts,nc1,1);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_VHF_FIT_BOOT                             */
/*                                                                           */
/*  Loop over multiple train / test splits for fitting and cross-validation. */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_fit_boot(fito)
     struct onode *fito;
{
  int i,j;
  int js,flag,nfb,fracn,ntrain,ntest,seed,nshuffle;
  float trfrac,*mse_trn,*mse_tst,*rvl_trn,*rvl_tst,mu_tr,mu_ts,sd_tr,sd_ts;
  char tstr[SLEN],tname[SLEN],*s1;
  struct onode *boo,*stgo;

  mylog(mylogf,"  MOD_VHF_FIT_BOOT\n");

  mod_vhf_boot_split_i = -1;
  nshuffle = 5;

  //
  //  There must be a <bootstrap> object with 'fracn' > 1, or else
  //  we will fit to *all* of the data.
  //
  fracn = 0;
  boo = onode_child_get_unique(fito,"bootstrap");
  if (boo != NULL){
    fracn  = onode_getpar_int_exit(boo,"split_frac");
    if (fracn > 1){
      nfb    = onode_getpar_int_exit(boo,"split_n");
      seed   = onode_getpar_int_exit(boo,"split_seed");
    }
    trfrac = (float)(fracn-1)/(float)fracn;  // Train fraction, e.g., 5/6
  }
  if (fracn == 0){
    mylog(mylogf,"    Fitting to all data, no bootstrap.\n");
    nfb    = 1;     // Only run one fit
    trfrac = 1.0;   // Use all data
    seed   = -1;
  }

  ntrain = my_rint(trfrac * (float)mod_vhf_fit_nrsp);
  ntest  = mod_vhf_fit_nrsp - ntrain;
  sprintf(ggstr,"    ntrain = %d   ntest = %d\n",ntrain,ntest);
  mylog(mylogf,ggstr);

  mse_trn = get_farray(nfb);
  rvl_trn = get_farray(nfb);
  mse_tst = get_farray(nfb);
  rvl_tst = get_farray(nfb);

  mod_vhf_fit_train_n = ntrain;  // Set the global value

  // Create list of indices to be shuffled
  mod_vhf_boot_shuff = get_permutation_iarray(mod_vhf_fit_nrsp);

  // Set max number of C1 inputs (needed before ...fit_slave_top)
  stgo = onode_child_get_unique(fito,"stage_2");
  mod_vhf_fit_c1_max_n = onode_getpar_int_exit(stgo,"max_c1_n");


  if (myid > 0){ // MM fit slave
    // ***** WYETH FIX THIS FOR THE '..._train_n' = all case
    // ***** WYETH FIX THIS FOR THE '..._train_n' = all case
    // ***** WYETH FIX THIS FOR THE '..._train_n' = all case
    mod_vhf_fit_slave_top(fito);  //  Go and wait for commands from master
    // THIS WILL NEVER RETURN HERE.
  }


  for(i=0;i<nfb;i++){  // For each split
    mod_vhf_boot_split_i = i;

    //
    //  Initialize for new fit
    //
    zero_iarray(mod_vhf_c1_flag,mod_vhf_c1_ntot);  // No C1's picked yet

    //
    //  (1) Set the (global) target response vector 'mod_vhf_fit_train_r[]'
    //
    if (ntrain < mod_vhf_fit_nrsp){
      // shuffle if we are not using the entire data set
      shuffle_iarray_return_seed(mod_vhf_boot_shuff,mod_vhf_fit_nrsp,nshuffle,
				 &seed);
    }
    for(j=0;j<ntrain;j++){
      js = mod_vhf_boot_shuff[j];
      mod_vhf_fit_train_r[j] = mod_vhf_fit_r[js];
    }


    if (myid == 0){ // Master MM fit
      mm_cmd_send(1,numproc-1,"get_shuffle_index",mylogf);  // Tell all
      // Send shuffle data
      mm_data_send_iarray(1,numproc-1,mod_vhf_boot_shuff,mod_vhf_fit_nrsp,
			  mylogf);
    }

    mod_vhf_fit_min2(fito); // find best two filters in the middle C1 grid

    flag = 0;
    while(flag == 0)
      flag = mod_vhf_fit_next(fito);  //  Fit additional units


    //
    // (3) Write .best.fit and .pp files for this split
    //
    mod_vhf_fit_split_out(fito,ntrain,ntest,i,&(rvl_trn[i]),&(mse_trn[i]),
			  &(rvl_tst[i]),&(mse_tst[i]));
  }

  //
  //  Statistics across all splits
  //
  mean_sdev_farray(rvl_trn,nfb,&mu_tr,&sd_tr);
  mean_sdev_farray(rvl_tst,nfb,&mu_ts,&sd_ts);

  printf("    Average performance:  %f (train) %f (test)\n",mu_tr,mu_ts);

  mod_vhf_fname(fito,"outfile_avgperf",".avg",tname); // Build outfile name

  sprintf(tstr,"%s  %f %f\n",mod_vhf_outfile_prefix,mu_tr,mu_ts);
  append_string_to_file(tname,tstr);


  myfree(mod_vhf_boot_shuff);
}
/**************************************-**************************************/
/*                                                                           */
/*                               MOD_VHF_01_FIT                              */
/*                                                                           */
/*  Fit the model to a set of target stimulus-response pairs.                */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_01_fit(m,s,r,fito)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     struct onode *fito;
{
  int i,j,k,ii,jj;
  int xn,yn,tn,n,tdn,xi,yi,fit_flag,tflag,nstim,tsamp,nrpt,s2n;
  float w,tscale,sscale,stim_samp,mse,*tdata,***stm,***stim1,**tframe;
  char *infile,tstr[SLEN],*fmode,*filemode,*outfile;
  struct vhf_grid *gt;

  mylog(mylogf,"  MOD_VHF_01_FIT\n");

  fit_flag = onode_getpar_int_dflt(fito,"fit",0);
  if (fit_flag == 0){
    mylog(mylogf,"    Fitting is turned off\n");
    return;
  }

  fmode = onode_getpar_chr_exit(fito,"mode");

  xn     = mod_vhf_xn;
  yn     = mod_vhf_yn;
  tn     = mod_vhf_tn;
  tscale = mod_vhf_tscale;
  sscale = mod_vhf_sscale;

  s2n    = mod_vhf_c2->s2_xn;

  // Initialize pointers to memory that may be re-sized during fitting
  mod_vhf_fit_x  = NULL;     // [..._xn] Current vector of param values
  mod_vhf_fit_dw = NULL;

  //
  //  Target responses - read from a file
  //
  infile = onode_getpar_chr_exit(fito,"data_file");
  tflag = read_farray(infile,&tdata,&tdn);
  if (tflag == 0){
    sprintf(ggstr,"    'data_file' = %s\n",infile);
    mylog(mylogf,ggstr);
    mylog_exit(mylogf,"MOD_VHF_01_FIT  Could not open data file");
  }
  sprintf(ggstr,"    %d responses in file:  %s\n",tdn,infile);
  mylog(mylogf,ggstr);
  mod_vhf_fit_r = tdata;  // [nstim] Target response data


  //
  //  C1 responses - read from a file, or pre-compute
  //
  filemode = onode_getpar_chr_dflt(fito,"file_c1_mode","none");
  if (strcmp(filemode,"read")==0){

    mod_vhf_fit_nrsp = tdn;  // Number of responses to fit

    //  Compute responses for all C1 units to all stimuli
    mod_vhf_resp_op("c1s_init",tdn); // Make all storage for S2 c1rs

    mod_vhf_file_c1(fito,"read");
    nstim = tdn;
  }else{

    nrpt = param_geti_exit(s->ppl,"stim_nrpt");
    nstim = s->ntr / nrpt;
    sprintf(ggstr,"    NSTIM = %d  (nrpt = %d)\n",nstim,nrpt);
    mylog(mylogf,ggstr);
    if (nstim != tdn)
      mylogx(mylogf,"MOD_VHF_01_FIT","Response count not equal to stimuli");


    //
    //  Generate stimuli for target responses, stm[nstim][xn][yn]
    //
    if (nstim != tdn)
      mylogx(mylogf,"MOD_VHF_01_FIT","Number of stimuli and responses differ");

    stim_samp = param_getf_exit(s->ppl,"stim_samp"); // Max samples/s
    tsamp = stm_get_tsamp(tscale,stim_samp);  // Pattern frames per real frame

    stm = (float ***)myalloc(nstim*sizeof(float **));

    for(i=0;i<nstim;i++){

      //printf("%2d %3d\n",shid[i],shrotdeg[i]);a

      //printf(" stim %d\n",i);
      mod_util_get_3d_stim(s,i,xn,yn,tn,sscale,tscale,tsamp,NULL);

      //  Copy the first frame into the stimulus array
      tframe = get_2d_farray(xn,yn);
      for(xi=0;xi<xn;xi++)
	for(yi=0;yi<yn;yi++)
	  tframe[xi][yi] = s->d[xi+1][yi+1][1];
      stm[i] = tframe;

      free_f3tensor(s->d,1,xn,1,yn,1,tn);

      // Write examples of the stimuli
      /**
	 if (i < 10){
	 sprintf(tstr,"zzz.Stim_%d.2d",i);
	 write_2d_data(tstr,stm[i],0,0,xn,yn,4,2,1,0);
	 }
      **/
    }

    sprintf(ggstr,"    Number of stimuli:  %d\n",nstim);
    mylog(mylogf,ggstr);
    sprintf(ggstr,"    Size of each image:  %d x %d\n",xn,yn);
    mylog(mylogf,ggstr);

    //
    //  Compute responses for all C1 units to all stimuli
    //
    mod_vhf_resp_op("c1s_init",nstim); // Make all storage for S2 c1rs

    for(i=0;i<nstim;i++){

      //  Extract the 'i'th stimulus into a tensor array, compute response
      stim1 = f3tensor(1,1,1,xn,1,yn);
      for(j=0;j<xn;j++){
	for(k=0;k<yn;k++){
	  stim1[1][j+1][k+1] = stm[i][j][k];
	}
      }

      // Write into the global s1 response array - no storage created
      mod_vhf_resp_s1(stim1);  // Compute response for S1 units

      // Use the global s1 response array to compute the c1 response array
      //printf("mod_vhf_c2->s2_xn = %d\n",mod_vhf_c2->s2_xn);

      for(ii=0;ii<s2n;ii++){
	for(jj=0;jj<s2n;jj++){
	  mod_vhf_resp_c1_all(mod_vhf_c2->s2[ii][jj]->c1rs[i],ii,jj);
	}
      }

      free_f3tensor(stim1,1,1,1,xn,1,yn);

      sprintf(ggstr,"  Done stim %d\n",i);
      mylog(mylogf,ggstr);
    }
    mod_vhf_fit_nrsp = nstim;  // Number of responses to fit
  }

  mod_vhf_fit_dj = get_3d_farray(mod_vhf_fit_nrsp,s2n,s2n);
  mod_vhf_fit_train_r = (float *)myalloc(mod_vhf_fit_nrsp*sizeof(float));

  //
  //  Check for writing C1 responses to a file
  //
  filemode = onode_getpar_chr_dflt(fito,"file_c1_mode","none");
  if (strcmp(filemode,"write")==0){
    mod_vhf_file_c1(fito,"write");
    printf("  Exiting after writing response file.\n\n");
    exit(0);
  }
  myfree(filemode);

  //
  //  Examine any relationships between the raw C1 responses and the
  //  response data.
  //
  if (strcmp(fmode,"c1_stat")==0)
    mod_vhf_relate_c1_targ();  // Exit when done


  if (0){
    mod_vhf_test_derivs();
    exit(0);
  }

  if (strcmp(fmode,"test_s2")==0){
    mod_vhf_deriv_test_s2();
    mod_vhf_s2_resp_all();
    exit(0);
  }


  //
  //  Do bootstrap fitting
  //
  mod_vhf_fit_boot(fito);


  myfree(infile);

  if (fit_flag == 1){
    printf("  Exiting after fit.\n");
    if (myid == 0){ // Master MM fit
      mm_cmd_send(1,numproc-1,"exit",mylogf);
      sleep(2); // Wait for sub-processes to exit
    }
    printf("\n");
    exit(0);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_VHF_MAKE_MPT                             */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_make_mpt(m,s,r)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
{
  int i,j,k;
  int xi,yi,stype,xn,yn,*ctx,*cty,ctn,txn,tyn,tzn,nori,nph;
  short inindex;
  float w,tdelay,**f;
  char tstr[SLEN];
  struct pop_top *mpt;
  struct pop_layer *pl;
  struct pop_syn *tsyn;
  struct pop_cell *pre,*post;
  struct onode *ot,*o;
  struct vhf_grid *gt;

  int nsz,gn;

  mylog(mylogf,"  MOD_VHF_MAKE_MPT\n");

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

  mpt->narea = 1;
  mpt->area = (struct pop_area **)myalloc(sizeof(struct pop_area *));

  printf(" WYETH *********** TESTING NEW CODE HERE\n");
  mpt->area[0] = popu_make_area_default("v1",mpt->xn,mpt->yn,mpt->sscale,20.0);
  /***
  mpt->area[0] = (struct pop_area *)myalloc(sizeof(struct pop_area));
  mpt->area[0]->name = strdup("v1");
  mpt->area[0]->x0   = 0.0;
  mpt->area[0]->y0   = 0.0;
  mpt->area[0]->xf   = 1.0;
  mpt->area[0]->yf   = 1.0;
  mpt->area[0]->xn   = mpt->xn;
  mpt->area[0]->yn   = mpt->yn;
  mpt->area[0]->umx  = 20.0; // WYETH ?
  mpt->area[0]->umy  = 20.0;
  mpt->area[0]->sscale = mpt->sscale;
  ***/

  mpt->nmap = 0;
  mpt->map = NULL;

  nsz  = mod_vhf_nsz;

  exit_error("MOD_VHF_MAKE_MPT","WYETH - must fix for <channels>");
  //nori = mod_vhf_nori;
  //nph  = mod_vhf_nph;
  gn   = mod_vhf_c1_gridn;

  //
  //  Layers
  //
  mpt->nlay = 1 + nsz + gn;  // 1 for stim pix layer
  //mpt->nlay = 1 + nsz; //+ gn;  // 1 for stim pix layer
  mpt->lay = (struct pop_layer **)myalloc(mpt->nlay*
					  sizeof(struct pop_layer *));

  pl = popu_make_layer_cart("stim",mpt->xn,mpt->yn,1);
  myfree(pl->laytype);
  pl->laytype = strdup("stim");
  pl->geomt = strdup("default");
  pl->x0 = pl->y0 = 0.0;
  pl->xf = pl->yf = 1.0;
  mpt->lay[0] = pl;
  mpt->lay[0]->ninlist = 0;
  mpt->lay[0]->inlist = NULL;

  for(i=0;i<nsz;i++){
    sprintf(tstr,"s1_%d",mod_vhf_tile_w[i]);
    txn = mpt->xn - mod_vhf_tile_w[i] + 1;
    tyn = mpt->yn - mod_vhf_tile_w[i] + 1;
    tzn = nori * nph;

    // *** WYETH - USE THIS ROUTINE TO SIMPLIFY THE CODE below
    // *** WYETH - USE THIS ROUTINE TO SIMPLIFY
    // *** WYETH - USE THIS ROUTINE TO SIMPLIFY
//struct pop_layer *popu_make_layer_template(name,xn,yn,zn,pa,lpre,fx,fy,fw,fn)

    
    pl = popu_make_layer_cart(tstr,txn,tyn,tzn);
    pl->geomt = strdup("default");
    pl->area = mpt->area[0];
    pl->x0 = pl->y0 = (int)(mod_vhf_tile_w[i]/2.0);
    pl->xf = pl->yf = 1.0;
    mpt->lay[i+1] = pl;
    mpt->lay[i+1]->ninlist = 1;
    mpt->lay[i+1]->inlist = (struct onode **)myalloc(sizeof(struct onode *));

    // Create an onode for the input
    ot = paramfile_onode_create_onode();
    ot->otype = strdup("input");
    o = paramfile_onode_get_init_onode("item",NULL,"pop_origin","stim",
				       1,NULL,NULL,1);
    onode_insert_child_at_end(ot,o);

    // *** WYETH - I ASSUME this 'template' type means that only the connections
    //   for the first unit need to be stored.
    o = paramfile_onode_get_init_onode("item",NULL,"type","template",
				       1,NULL,NULL,1);
    onode_insert_child_at_end(ot,o);
    mpt->lay[i+1]->inlist[0] = ot;

    //
    //  Write connections for first unit
    //
    stype = 1;
    tdelay = 0.0;
    inindex = 0;
    w = 1.0;

    xn = mod_vhf_xn;
    yn = mod_vhf_yn;

    post = &(mpt->lay[i+1]->c[0][0][0]);

    exit_error("MOD_VHF_MAKE_MPT","WYETH - 2 must 2 fix 2 for 2 <channels>");
    //f = mod_vhf_f[i][0][0];    // [size][ori][phase][x][y] Filters
    for(xi=0;xi<xn;xi++){
      for(yi=0;yi<yn;yi++){
	if (fabs(f[xi][yi]) > 0.05){
	  //w = fabs(f[xi][yi]);
	  w = 0.5*(1.0 + f[xi][yi]);
	  //printf(" f[][] = %f\n",f[xi][yi]);

	  pre = &(mpt->lay[0]->c[xi][yi][0]);

	  tsyn = pop_cell_add_synapse(pre,post,stype,w,tdelay,inindex);
	}
      }
    }
  }

  k = nsz + 1;  // number of layers defined so far

  for(i=0;i<gn;i++){
    gt = mod_vhf_c1_grid[i];
    sprintf(tstr,"%s",gt->name);
    txn = gt->xn;
    tyn = gt->yn;
    tzn = nori;
    pl = popu_make_layer_cart(tstr,txn,tyn,tzn);
    pl->geomt = strdup("default");
    pl->area = mpt->area[0];
    pl->x0 = 10; // WYETH HACK
    pl->y0 = 20; // WYETH HACK -- NOTE error in iModel Browser - y is inverted?
    pl->xf = gt->dx;
    pl->yf = gt->dy;

    mpt->lay[i+k] = pl;
    mpt->lay[i+k]->ninlist = 1;
    mpt->lay[i+k]->inlist = (struct onode **)myalloc(sizeof(struct onode *));

    // Create an onode for the input
    ot = paramfile_onode_create_onode();
    ot->otype = strdup("input");

    /******************************WYETH HERE ********/
    /******************************WYETH HERE ********/
    sprintf(tstr,"s1_%d",32);  // WYETH HACK - NEED TO COMPUTE TRUE pre-SYN

    o = paramfile_onode_get_init_onode("item",NULL,"pop_origin",tstr,
				       1,NULL,NULL,1);
    onode_insert_child_at_end(ot,o);
    o = paramfile_onode_get_init_onode("item",NULL,"type","template",
				       1,NULL,NULL,1);
    onode_insert_child_at_end(ot,o);
    mpt->lay[i+k]->inlist[0] = ot;

    //
    //  Write connections for first unit
    //
    stype = 1;
    tdelay = 0.0;
    inindex = 0;
    w = 1.0;

    //xn = mod_vhf_xn;
    //yn = mod_vhf_yn;

    post = &(mpt->lay[i+k]->c[0][0][0]);

    /******************************WYETH HACK ********/
    for(xi=0;xi<10;xi++){
      for(yi=0;yi<10;yi++){
	pre = &(mpt->lay[1]->c[xi+45][yi+45][0]); // HACK
	tsyn = pop_cell_add_synapse(pre,post,stype,w,tdelay,inindex);
      }
    }
  }

  popu_mar_write(mpt,m->marfile,4);  // Version 4
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_VHF_CONFIG_S1_TILE                          */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_config_s1_tile(s1o)
     struct onode *s1o;  // Contains the <s1> object from .moo file
{
  int i,k;
  int xn,yn,nsz,*ntile;
  struct onode *to;

  mylog(mylogf,"  MOD_VHF_CONFIG_S1_TILE\n");

  xn   = mod_vhf_xn;  // Convenient values
  yn   = mod_vhf_yn;

  //
  //  Read the S1 tile description
  //
  nsz  = mod_vhf_nsz  = onode_count_otype(s1o,"tile");
  sprintf(ggstr,"    %d S1 tile sizes\n",nsz);
  mylog(mylogf,ggstr);

  mod_vhf_tile_w  = (int *)myalloc(nsz*sizeof(int));
  mod_vhf_tile_id = (int *)myalloc(nsz*sizeof(int));
  mod_vhf_tile_x0 = (int *)myalloc(nsz*sizeof(int));
  mod_vhf_tile_y0 = (int *)myalloc(nsz*sizeof(int));
  mod_vhf_tile_xn = (int *)myalloc(nsz*sizeof(int));
  mod_vhf_tile_yn = (int *)myalloc(nsz*sizeof(int));

  k = 0;
  to = onode_get_next_type(s1o->o,"tile");
  while(to != NULL){
    mod_vhf_tile_w[k]  = onode_getpar_int_exit(to,"size");
    mod_vhf_tile_id[k] = onode_getpar_int_exit(to,"name");
    k += 1;
    to = onode_get_next_type(to->next,"tile");
  }


  //
  //  Print out some statistics about the S1 configation
  //

  //  Compute number of S1 units along x-dimension for each tile size
  ntile = (int *)myalloc(nsz*sizeof(int));
  for(i=0;i<nsz;i++)
    ntile[i] = 1 + xn - mod_vhf_tile_w[i];

  //  Compute total number of RFs
  k = 0;
  for(i=0;i<nsz;i++){
    k += ntile[i] * ntile[i];  // Number for size band 'i'
    sprintf(ggstr,"    tile %d  %3d pix:    %d x %d\n",mod_vhf_tile_id[i],
	   mod_vhf_tile_w[i],ntile[i],ntile[i]);
    mylog(mylogf,ggstr);
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_VHF_CONFIG_S1_CHANS                         */
/*                                                                           */
/*****************************************************************************/
struct vhf_chans *mod_vhf_config_s1_chans(co)
     struct onode *co;  // Contains the <channels> object
{
  int i,j,k,ii,jj;
  int xn,yn,nori,nph,nsz,*ntile,ntot,dumpflag;
  float theta,phase,sdo,sdp,sf,xc,yc,fw,**t,maskr,norm_mask_sd_f,norm_targ;
  float *****f,******ft,*****fts,***tft,**tfts,*****rr,sump,sumn,*tcol;
  char tname[SLEN];
  struct vhf_chans *ch;

  mylog(mylogf,"  MOD_VHF_CONFIG_S1_CHANS\n");

  ch = (struct vhf_chans *)myalloc(sizeof(struct vhf_chans));

  ch->name   = onode_getpar_chr_exit(co,"name");
  ch->type   = onode_getpar_chr_dflt(co,"type","gabor");
  ch->nori   = onode_getpar_int_dflt(co,"nori",1);
  ch->nph    = onode_getpar_int_exit(co,"nphase");
  ch->ph0    = onode_getpar_flt_dflt(co,"phase0",0.0);
  ch->maxph_flag  = onode_getpar_int_dflt(co,"max_over_phase",1);
  dumpflag   = onode_getpar_int_dflt(co,"write_filters_exit",0);

  if ((ch->maxph_flag == 1)||(ch->maxph_flag == 2))
    ch->nch = ch->nori;
  else
    ch->nch = ch->nph * ch->nori;

  if (strcmp(ch->type,"dog")==0){
    // *** WYETH HACK - these struct names (for Gabor) are  re-used for DoG
    // *** WYETH HACK - these struct names (for Gabor) are  re-used for DoG
    ch->f_o       = onode_getpar_flt_exit(co,"dog_sd_surr_f");
    ch->f_p       = onode_getpar_flt_exit(co,"dog_sd_center_f");
    ch->f_f       = onode_getpar_flt_dflt(co,"dog_amp_surr",1.0);
  }else{
    ch->f_o       = onode_getpar_flt_exit(co,"gabor_sd_orth_f");
    ch->f_p       = onode_getpar_flt_exit(co,"gabor_sd_para_f");
    ch->f_f       = onode_getpar_flt_exit(co,"gabor_sf_f");
  }
  norm_targ       = onode_getpar_flt_exit(co,"filter_norm_targ");
  norm_mask_sd_f  = onode_getpar_flt_dflt(co,"norm_mask_sd_f",2.5);
  ch->norm_flag   = onode_getpar_int_dflt(co,"norm_type",0);
  ch->norm_eps    = onode_getpar_flt_exit(co,"norm_eps");

  xn   = mod_vhf_xn;
  yn   = mod_vhf_yn;
  nori = ch->nori;
  nph  = ch->nph;
  nsz  = mod_vhf_nsz;

  if (ch->norm_flag != 0){
    ch->norm_mask_x = (int **)myalloc(nsz*sizeof(int *));
    ch->norm_mask_y = (int **)myalloc(nsz*sizeof(int *));
    ch->norm_mask_n = (int  *)myalloc(nsz*sizeof(int));
    mod_vhf_s1_norm_flag = 1;  // Set to 1 if *any* <channels> needs this
  }

  /***
  //
  //  (2) Print out some statistics about the S1 configation
  //

  //  Compute number of S1 units along x-dimension for each tile size
  ntile = (int *)myalloc(nsz*sizeof(int));
  for(i=0;i<nsz;i++)
    ntile[i] = 1 + xn - mod_vhf_tile_w[i];

  //  Compute total number of RFs
  k = 0;
  for(i=0;i<nsz;i++){
    k += ntile[i] * ntile[i];  // Number for size band 'i'
    printf("    tile %d  %3d pix:    %d x %d\n",mod_vhf_tile_id[i],
	   mod_vhf_tile_w[i],ntile[i],ntile[i]);
  }
  ntot = k * nori * nph;
  printf("    n_ori   = %d\n",nori);
  printf("    n_phase = %d\n",nph);
  printf("    Total number of RFs to apply:  %d\n",ntot);
  ***/


  //
  //  (3) Create the S1 filters
  //
  f = (float *****)myalloc(nsz*sizeof(float ****));
  ft = (float ******)myalloc(nsz*sizeof(float *****));
  fts = (float *****)myalloc(nsz*sizeof(float ****));
  for(i=0;i<nsz;i++){
    fw = (float)mod_vhf_tile_w[i];

    sdo = fw / ch->f_o;  // (pix)
    sdp = fw / ch->f_p;  // (pix)
    sf  = ch->f_f / fw;  // (cyc/pix)

    // This can go outside the loop
    xc = xn / 2.0;
    yc = yn / 2.0;

    // Get the 2D filters
    f[i]   = (float  ****)myalloc(nori*sizeof(float ***));
    ft[i]  = (float *****)myalloc(nori*sizeof(float ****));
    fts[i] = (float  ****)myalloc(nori*sizeof(float ***));
    for(j=0;j<nori;j++){
      f[i][j]   = (float  ***)myalloc(nph*sizeof(float **));
      ft[i][j]  = (float ****)myalloc(nph*sizeof(float ***));
      fts[i][j] = (float  ***)myalloc(nph*sizeof(float **));
      theta = j * 180.0/nori;

      for(k=0;k<nph;k++){
	phase = ch->ph0 + k * (360.0/(float)nph);

	if (strcmp(ch->type,"dog")==0){
	  f[i][j][k] = dog_2d_space_raw(xn,yn,xc,yc,sdp,sdo,ch->f_f);
	  if (my_rint(phase) == 180)
	    multiply_2d_farray(f[i][j][k],xn,yn,-1.0);
	  else if (my_rint(phase) != 0)
	    exit_error("MOD_VHF_CONFIG_S1_TILE","DoG phase must be 0 or 180");
	    
	}else{
	  // "gabor"
	  f[i][j][k] = gabor_2d_space_raw(xn,yn,xc,yc,sdo,sdp,sf,theta,phase);
	}

	// Scale filter so that the sum of its absolute value is 'norm_targ'
	if (norm_targ != 0.0)
	  kernu_norm_2d_filter(f[i][j][k],xn,yn,norm_targ,&sump,&sumn);
	ft[i][j][k] = f3tensor(1,1,1,xn,1,yn);
	fts[i][j][k] = matrix(1,1,1,2*xn);

	// Convenient pointers
	t    =   f[i][j][k];
	tft  =  ft[i][j][k];
	tfts = fts[i][j][k];

	//
	//  Compute the FT of the filters
	//
	for(ii=1;ii<=xn;ii++){
	  for(jj=1;jj<=yn;jj++){
	    tft[1][ii][jj] = t[ii-1][jj-1];
	  }
	}

	contort_real_3d_farray(tft,1,1,1,1,xn,yn);
	rlft3(tft,tfts,1,xn,yn,1);
      }
    }


    //
    //  Get the normalization masks
    //
    if (ch->norm_flag != 0){
      if (sdp >= sdo)  // Mask radius is, e.g., 2.5 * max SD
	maskr = sdp * norm_mask_sd_f;
      else
	maskr = sdo * norm_mask_sd_f;
      mod_vhf_norm_mask(maskr,&(ch->norm_mask_x[i]),
                              &(ch->norm_mask_y[i]),
                              &(ch->norm_mask_n[i]));
      sprintf(ggstr,"    maskr = %f   tile_w = %f\n",maskr,fw);
      mylog(mylogf,ggstr);
    }
  }

  ch->f = f;
  ch->ft = ft;
  ch->fts = fts;

  //
  //  Write out filters
  //
  if (dumpflag == 1){
    remove_file("zzz.f_xcut.pl");

    for(i=0;i<nsz;i++){
      for(j=0;j<nori;j++){
	for(k=0;k<nph;k++){
	  t = f[i][j][k];  // Pointer to 2D filter
	  sprintf(tname,"zzz.f_%d_%d_%d.2d",i,j,k);
	  write_2d_data(tname,t,0,0,xn,yn,4,2,1,0);

	  printf("    %s  Area = %f\n",tname,sum_2d_farray(t,0,xn,0,yn));

	  tcol = get_column_2d_farray(t,xn,yn,yn/2);
	  append_farray_plot("zzz.f_xcut.pl",tname,tcol,xn,1);
	  myfree(tcol);
	}
      }
    }
    printf("  Wrote filters to 'zzz.f_...2d'.  Exiting.\n\n");
    exit(0);
  }

  return ch;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_VHF_CONFIG_S1_TOP                           */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_config_s1_top(s1o)
     struct onode *s1o;  // Contains the <s1> object from .moo file
{
  int i;
  int nchans;
  struct onode *to;

  mylog(mylogf,"  MOD_VHF_CONFIG_S1_TOP\n");

  //
  //  Get the basic tile information, this will be needed to build the filters
  //  in each <channels> below.
  //
  mod_vhf_config_s1_tile(s1o);

  nchans = onode_count_otype(s1o,"channels");
  mod_vhf_chans = (struct vhf_chans **)myalloc(nchans*
					       sizeof(struct vhf_chans *));
  mod_vhf_s1_nchans = nchans;  // Save this globally
  sprintf(ggstr,"    Number of <channels> objects:  %d\n",nchans);
  mylog(mylogf,ggstr);

  //
  //  Build each set of channels
  //
  to = onode_get_next_type(s1o->o,"channels");
  i = 0;
  mod_vhf_c1_nch = 0;
  while(to != NULL){
    mod_vhf_chans[i] = mod_vhf_config_s1_chans(to); // Build a 'vhf_chans'
    mod_vhf_chans[i]->nch0 = mod_vhf_c1_nch;

    mod_vhf_c1_nch += mod_vhf_chans[i]->nch;
    to = onode_get_next_type(to->next,"channels");
    i += 1;
  }
  sprintf(ggstr,"    %d channels at each C1 grid point\n",mod_vhf_c1_nch);
  mylog(mylogf,ggstr);

  if (mod_vhf_c1_nch == 0)
    mylog_exit(mylogf,"MOD_VHF_CONFIG_S1_TOP  No <channels> object found.");

  //
  //  Create storage for responses
  //
  mod_vhf_resp_op("s1_init",0);

}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_VHF_CONFIG_C1                             */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_config_c1(c1o)
     struct onode *c1o;  // Contains the <c1> object from .moo file
{
  int i,j,k,ii,jj;
  int xi,yi,ti,gi,ci,gn,gxn,gyn,gw,gov,ntot,nch,ntile,s1w,xx0,yy0,s2xn;
  float ******rr,s1w_avg;
  struct onode *to;
  struct vhf_grid **g,*gt;

  mylog(mylogf,"  MOD_VHF_CONFIG_C1\n");

  nch = mod_vhf_c1_nch;  // Number of channels at each C1 grid point

  gn = onode_count_otype(c1o,"grid");  // Number of C1 grids, e.g., there are
                                       // typically 3:  2x2, 3x3, and 4x4

  sprintf(ggstr,"    %d C1 grids\n",gn);  // Print out the number of grids
  mylog(mylogf,ggstr);

  //  Create a structure 'g' to hold all the C1 grid information.
  g = (struct vhf_grid **)myalloc(gn*sizeof(struct vhf_grid *));

  k = 0;
  ntot = 0;
  to = onode_get_next_type(c1o->o,"grid");  // 'to' contains grid parameters
  while(to != NULL){  // For each grid (the current grid is defined by 'to')
    gt = g[k] = (struct vhf_grid *)myalloc(sizeof(struct vhf_grid));
    gt->xn = onode_getpar_int_exit(to,"xn");  // Width of grid in C1 units
    gt->yn = onode_getpar_int_dflt(to,"yn",gt->xn);  // Height of grid

    gt->name = onode_getpar_chr_exit(to,"name");

    gxn = gt->xn;
    gyn = gt->yn;

    ntot += gxn * gyn;  // Accumulate total number of C1 positions in all grids

    gt->dx = onode_getpar_int_exit(to,"dx");  // Separation of C1 units in grid
    gt->dy = onode_getpar_int_dflt(to,"dy",gt->dx);

    gt->w = onode_getpar_int_exit(to,"width");
    gt->h = onode_getpar_int_dflt(to,"height",gt->w);

    // Total number of stim pixels along x-axis covered by RFs in this grid
    gt->totpixw = gt->w + (gxn - 1)*gt->dx;
    gt->totpixh = gt->h + (gyn - 1)*gt->dy;

    sprintf(ggstr,"    %d x %d  size %d  (%d pix total width)\n",gxn,gyn,gt->w,
	    gt->totpixw);
    mylog(mylogf,ggstr);

    gt->tile_list = onode_getpar_int_list(to,"tiles",&ntile);
    gt->ntile = ntile;

    gt->s1xn = (int   *)myalloc(ntile*sizeof(int));
    gt->s1yn = (int   *)myalloc(ntile*sizeof(int));
    gt->s1x0 = get_3d_iarray(ntile,gxn,gyn);
    gt->s1y0 = get_3d_iarray(ntile,gxn,gyn);


    s1w_avg = 0.0;
    for(i=0;i<ntile;i++){
      // Get the index of this tile in the list of all S1 tiles
      ti = get_index_search_iarray(mod_vhf_tile_id,mod_vhf_nsz,
				   gt->tile_list[i]);
      s1w = mod_vhf_tile_w[ti];  // Width of S1 tile
      gt->s1xn[i] = gt->w - s1w + 1;
      gt->s1yn[i] = gt->h - s1w + 1;  // Note, 's1w' because tiles are square
      //
      // Ex:    3 =   6   -  4  + 1,  Thus, 3 tiles of length 4 gives RF of 6
      //  1  ----
      //  2   ----
      //  3    ----
      //     123456
      //
      if ((gt->s1xn[i] < 1) || (gt->s1yn[i] < 1))
	mylog_exit(mylogf,"MOD_VHF_CONFIG_C1  Impossible grid geometry.");

      s1w_avg += (float)s1w;


      //
      //  Q:  Where does the left-most C1 unit take it's leftmost S1 tile?
      //
      //    . . s s s s    Example, 4x4 S1 tile, with 2 pix margin 
      //    . . s s * s    Diagram shows lower left corner
      //    . . s s s s 
      //    . . s s s s    ASSUMPTION:  that the resulting convolution 
      //    . . . . . .    value for this tile ends up at the pixel below
      //    . . . . . .    element (size/2,size/2) of this tile, marked * here,
      //    0 1 2 3 4 5    which for odd size would be the central element
      //
      //  A:  margin + size/2  (here:  2 + 4/2 = 4)
      //
      xx0 = (mod_vhf_xn - gt->totpixw)/2 + (int)(s1w/2);
      yy0 = (mod_vhf_yn - gt->totpixh)/2 + (int)(s1w/2);

      //
      //  Define the rectangular region of these tiles that are used by C1
      //  *** NOTE, if there is a grid of S2 units, these coordinates will
      //      need to be updated to account for all S2 units.
      //
      //  THESE are the values that will need to be used in the update:
      // mod_vhf_c2->s2_xn    e.g. 3x3 grid of S2 units
      // mod_vhf_c2->s2_dx    e.g. 30 pix offset
      //

      mod_vhf_tile_x0[ti] = xx0;  // Left most S1 position
      mod_vhf_tile_y0[ti] = yy0;

      mod_vhf_tile_xn[ti] = gt->s1xn[i] + (gxn-1)*gt->dx;  // Width in S1 units
      mod_vhf_tile_yn[ti] = gt->s1yn[i] + (gyn-1)*gt->dy;  // Height in S1 units

      if ((xx0 < 0) || (yy0 < 0))
	mylog_exit(mylogf,"MOD_VHF_CONFIG_C1  grid exceeds stimulus area.");

      sprintf(ggstr,
	      "      tile %d   %d x %d S1 units (lower left S1: %d, %d)\n",
	      gt->tile_list[i],gt->s1xn[i],gt->s1yn[i],xx0,yy0);
      mylog(mylogf,ggstr);
      
      for(xi=0;xi<gxn;xi++){
	for(yi=0;yi<gyn;yi++){
	  gt->s1x0[i][xi][yi] = xx0 + xi*gt->dx;
	  gt->s1y0[i][xi][yi] = yy0 + yi*gt->dy;
	}
      }
    }

    gt->avgtw = s1w_avg / (float)ntile;  // Avg tile width (pix)
    //printf("gt->avgtw = %f\n",gt->avgtw);

    //printf("  xn,yn %d %d  wid,ovr %d %d\n",gxn,gyn,gw,gov);

    k += 1;
    to = onode_get_next_type(to->next,"grid");
  }

  mod_vhf_c1_ntot = ntot * nch;
  sprintf(ggstr,"    C1 has %d total units\n",mod_vhf_c1_ntot);
  mylog(mylogf,ggstr);

  //
  //  Create 6-D array to store responses  r[s2xn][s2xn][gn][gxn][gyn][nch]
  //
  s2xn = mod_vhf_c2->s2_xn;    // the S2 RF shifts

  rr = (float ******)myalloc(s2xn*sizeof(float *****));
  for(ii=0;ii<s2xn;ii++){
    rr[ii] = (float *****)myalloc(s2xn*sizeof(float ****));
    for(jj=0;jj<s2xn;jj++){
      rr[ii][jj] = (float ****)myalloc(gn*sizeof(float ***));
      for(i=0;i<gn;i++){
	gxn = g[i]->xn;
	rr[ii][jj][i] = (float ***)myalloc(gxn*sizeof(float **));
	for(j=0;j<gxn;j++){
	  gyn = g[i]->yn;
	  rr[ii][jj][i][j] = (float **)myalloc(gyn*sizeof(float *));
	  for(k=0;k<gyn;k++){
	    rr[ii][jj][i][j][k] = (float *)myalloc(nch*sizeof(float));
	  }
	}
      }
    }
  }

  //
  //  The next chunk of code conceptually defines a single 1-D list to organize
  //  all C1 units.  For each C1 unit, say it is index 'k' within a 1-D list, 
  //  we need to know where it falls within the higher dimensional
  //  representation, which includes knowing the following 4 things:
  //
  //     mod_vhf_c1_ndx[k][0]  - The x-index within the local C1 grid
  //                    " [1]  - The y-index within the local C1 grid
  //                    " [2]  - The 'channel' index [0,...,nch-1]
  //                    " [3]  - The index of the grid,
  //                               e.g., 0= 2x2, 1= 3x3, 2= 4x4

  //
  //  Set up a list of indices for all C1 units
  //
  mod_vhf_c1_ndx = get_2d_iarray(mod_vhf_c1_ntot,4);
  k = 0;
  for(gi=0;gi<gn;gi++){  // For each C1 grid
    gxn = g[gi]->xn;     // width of grid
    gyn = g[gi]->yn;     // height of grid

    g[gi]->ndx0 = k;     // Save start index for this grid w/i the ..._c1_ndx

    for(i=0;i<gxn;i++){
      for(j=0;j<gyn;j++){
	for(ci=0;ci<nch;ci++){
	  mod_vhf_c1_ndx[k][0] = i;   // x-axis index
	  mod_vhf_c1_ndx[k][1] = j;   // y-axis index
	  mod_vhf_c1_ndx[k][2] = ci;  // channel index
	  mod_vhf_c1_ndx[k][3] = gi;  // grid ID
	  k += 1;
	}
      }
    }
  }

  // Set up a flag array to tell which C1 units have been picked
  mod_vhf_c1_flag = get_const_iarray(mod_vhf_c1_ntot,0);

  // Array to hold indices w.r.t. 'mod_vhf_c1_ndx' for chosen C1 inputs
  mod_vhf_fit_ndx = (int *)myalloc(mod_vhf_c1_ntot*sizeof(int));

  // Array to hold best MSE values after each C1 unit is chosen
  mod_vhf_fit_mse = (float *)myalloc(mod_vhf_c1_ntot*sizeof(float));

  // WYETH - THREE MORE ARRAYS TO HOLD PERFORMANCE vs. C1 unit count
  mod_vhf_fit_mse_ts = (float *)myalloc(mod_vhf_c1_ntot*sizeof(float));
  mod_vhf_fit_rvl    = (float *)myalloc(mod_vhf_c1_ntot*sizeof(float));
  mod_vhf_fit_rvl_ts = (float *)myalloc(mod_vhf_c1_ntot*sizeof(float));


  mod_vhf_c1_grid  = g;   // [gn] pointers
  mod_vhf_c1_gridn = gn;
  mod_vhf_c1_r     = rr;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_VHF_CONFIG_S2                             */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_config_s2(s2o)
     struct onode *s2o;  // Contains the <s2> object from .moo file
{
  int i,k;
  int ng,n,sn,*ival;
  float sig_s,sig_a,sig_b;
  char *outfile;
  struct onode *to;
  struct vhf_grid  **g;
  struct vhf_s2_par  *s2;

  mylog(mylogf,"  MOD_VHF_CONFIG_S2\n");

  s2 = (struct vhf_s2_par *)myalloc(sizeof(struct vhf_s2_par));

  //  Let's read in an overall 'type' of non-linearity

  s2->nltype = onode_getpar_chr_dflt(s2o,"nonlin","sigmoid0"); //  'sigmoid0'
  if (strcmp(s2->nltype,"sigmoid0")==0){
    s2->s    = onode_getpar_flt_exit(s2o,"sigmoid_s");
    s2->a    = onode_getpar_flt_exit(s2o,"sigmoid_alpha");
    s2->b    = onode_getpar_flt_exit(s2o,"sigmoid_beta");
  }else if (strcmp(s2->nltype,"none")==0){
    ; // Do nothing
  }else
    exit_error("MOD_VHF_CONFIG_S2","Unknown 'nonlin' in <s2>");
    
  
  // *** WYETH WARNING see 'KSMALL' which might be used instead
  s2->eps2 = onode_getpar_flt_dflt(s2o,"norm_eps",0.0001);
    

  outfile = onode_getpar_chr_dflt(s2o,"postplot_file",NULL);
  if (outfile != NULL){
    if (strcmp(outfile,"NULL")==0){  // If it is "NULL" set to NULL
      myfree(outfile);
      outfile = NULL;
    }
  }

  //
  //  OPTIONAL:  Plot the sigmoid
  //
  if (0){
    float x,s;
    char tstr[SLEN];

    if (strcmp(s2->nltype,"sigmoid0")!=0)
      exit_error("MOD_VHF_CONFIG_S2","Only works for 'sigmoid0' s2 nonlin");

    n = 200;
    for(i=0;i<n;i++){
      x = (float)(i-n/2)/(float)(n/10);
      s = s2->s / (1.0 + exp(-s2->a*(x-s2->b)));
      sprintf(tstr,"%f %f\n",x,s);
      append_string_to_file("zzz.sigmoid.pl",tstr);
    }
    exit(0);
  }

  ng = mod_vhf_c1_gridn;
  g  = mod_vhf_c1_grid;

  //
  //  Count number of C1 inputs to this S2 unit
  //
  s2->nw = 0;
  to = onode_get_next_type(s2o->o,"item");
  while(to != NULL){
    for(i=0;i<ng;i++){
      if (compare_prefix_string_order(g[i]->name,to->name) == 1){
	s2->nw += 1;
      }
    }
    to = onode_get_next_type(to->next,"item");
  }

  //
  //  Allocate storage for the weights and coordinates from C1 to S2
  //
  s2->c1w = (struct vhf_c1w **)myalloc(s2->nw*sizeof(struct vhf_c1w *));
  for(i=0;i<s2->nw;i++)
    s2->c1w[i] = (struct vhf_c1w *)myalloc(sizeof(struct vhf_c1w));


  //
  //  Read and save the coordinates and weights for inputs to S2
  //
  k = 0;
  to = onode_get_next_type(s2o->o,"item");
  while(to != NULL){
    for(i=0;i<ng;i++){
      if (compare_prefix_string_order(g[i]->name,to->name) == 1){
	sn = strlen(g[i]->name);
	ival = carray_util_parse_underscore(&(to->name[sn+1]),&n);
	if (n != 3)
	  mylog_exit(mylogf,"MOD_VHF_CONFIG_S2  Expecting 3 indices.");

	s2->c1w[k]->gi = i;
	s2->c1w[k]->xi = ival[0];
	s2->c1w[k]->yi = ival[1];
	s2->c1w[k]->ci = ival[2];
	s2->c1w[k]->w  = atof(to->val);
	k += 1;

	//printf(" %s %d %d %d  :  %s\n",g[i]->name,ival[0],ival[1],ival[2],
	//to->val);
      }
    }
    to = onode_get_next_type(to->next,"item");
  }

  mylog(mylogf,"    S2 configuration:\n");

  if (strcmp(s2->nltype,"sigmoid0")==0){
    mylog(mylogf,"      Sigmoid\n");
    sprintf(ggstr,"        s  %f\n",s2->s);
    mylog(mylogf,ggstr);
    sprintf(ggstr,"        a  %f\n",s2->a);
    mylog(mylogf,ggstr);
    sprintf(ggstr,"        b  %f\n",s2->b);
    mylog(mylogf,ggstr);
  }else if (strcmp(s2->nltype,"none")==0){
    mylog(mylogf,"      No nonlinearity applied\n");
  }

  sprintf(ggstr,"      C1 inputs: %d\n",s2->nw);
  mylog(mylogf,ggstr);
  mylog(mylogf,"        # grid xi yi chan  weight\n");
  for(i=0;i<s2->nw;i++){
    sprintf(ggstr,"       %2d  %2d  %2d %2d %2d  %9.5f\n",i,s2->c1w[i]->gi,
	    s2->c1w[i]->xi,s2->c1w[i]->yi,s2->c1w[i]->ci,s2->c1w[i]->w);
    mylog(mylogf,ggstr);
  }

  mod_vhf_s2 = s2;  // Store global

  if (outfile != NULL)
    mod_vhf_pp_write(outfile,0);

  // WYETH - CALL TEST DERIV HERE
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_VHF_CONFIG_C2                             */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_config_c2(c2o)
     struct onode *c2o;  // Contains the <c2> object from .moo file
{
  int i,j;
  int gxn,dx; // ,gwidth;
  struct onode *go;

  mylog(mylogf,"  MOD_VHF_CONFIG_C2\n");

  mod_vhf_c2 = (struct vhf_c2 *)myalloc(sizeof(struct vhf_c2));
  if (c2o != NULL){
    go = onode_child_get_unique(c2o,"grid");  // Pointer to the grid node
    gxn = onode_getpar_int_exit(go,"xn"); // Grid width in S2 units, e.g. 3 x 3
    dx  = onode_getpar_int_exit(go,"dx"); // Offset of neighboring S2's (30 pix)
  }else{
    //
    //  No 'C2' definition
    //
    gxn = 1;
    dx = 0;
  }

  //
  //  Create S2 grid
  //   - to hold C1 responses for pre-computing and fitting
  //
  mod_vhf_c2->s2 = (struct vhf_s2 ***)myalloc(gxn*sizeof(struct vhf_s2 **));
  for(i=0;i<gxn;i++){
    mod_vhf_c2->s2[i] = (struct vhf_s2 **)myalloc(gxn*sizeof(struct vhf_s2 *));
    for(j=0;j<gxn;j++){
      mod_vhf_c2->s2[i][j] = (struct vhf_s2 *)myalloc(sizeof(struct vhf_s2));
    }
  }

  //gwidth = (gxn-1)*dx + onode_getpar_int_exit(go,"width"); 

  mod_vhf_c2->s2_xn = gxn;  // Store these values in global params
  mod_vhf_c2->s2_dx = dx;

  // Make response storage for the grid of S2 units
  mod_vhf_s2_r = get_2d_farray(gxn,gxn);

  sprintf(ggstr,"    Grid is  %d x %d\n",gxn,gxn);
  mylog(mylogf,ggstr);
  sprintf(ggstr,"    dx = %d pixels\n",dx);
  mylog(mylogf,ggstr);
  //printf("    gwidth = %d pixels\n", gwidth);
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_VHF_01_PREP                              */
/*                                                                           */
/*  Prepare constructs that will be needed when running trials.  Do things   */
/*  here that only need to be done once, to avoid repeating them for each    */
/*  trial.                                                                   */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_01_prep(m,s,r)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
{
  int xn,yn,tn;
  float tscale,sscale;
  struct onode *fito,*s1o,*c1o,*s2o,*c2o;

  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  mylog(mylogf,"  MOD_VHF_01_PREP\n");

  xn     = mod_vhf_xn     = onode_getpar_int_exit(m->o,"xn");
  yn     = mod_vhf_yn     = onode_getpar_int_exit(m->o,"yn");
  tn     = mod_vhf_tn     = onode_getpar_int_exit(m->o,"tn");
  tscale = mod_vhf_tscale = onode_getpar_flt_exit(m->o,"tscale");
  sscale = mod_vhf_sscale = onode_getpar_flt_exit(m->o,"sscale");

  mod_vhf_ft_resp_flag = 1;  // 1-Use FFT - FOR NOW, THIS IS NO OTHER OPTION

  if (mod_vhf_fit_opt == NULL)
    mod_vhf_fit_opt = (float *)myalloc(MAX_PAR*sizeof(float));

  mod_vhf_outfile_prefix = param_getc_dflt(r->ppl,"outfile","zz");

  //
  //  S1 layer
  //
  s1o = onode_child_get_unique(m->o,"s1");
  mod_vhf_config_s1_top(s1o);


  //
  //  C2 layer
  //
  c2o = onode_child_get_unique(m->o,"c2");
  mod_vhf_config_c2(c2o);  // Will handle NULL case


  //
  //  C1 layer
  //
  c1o = onode_child_get_unique(m->o,"c1");
  if (c1o == NULL)
    mylog_exit(mylogf,"MOD_VHF_01_PREP  <c1> not found in .moo file.");
  else
    mod_vhf_config_c1(c1o);


  //
  //  S2 layer
  //
  s2o = onode_child_get_unique(m->o,"s2");
  if (s2o != NULL){
    mod_vhf_config_s2(s2o);
  }else{
    //
    // Note, even if the code is running in fitting mode, there is a
    // parameter, norm_eps, that is stored in 'mod_vhf_s2' which is needed.
    // Thus, we will insist (for now) that an <s2> object is present.
    // *** HOWEVER - it appears that 'norm_eps' is not used, and that KSMALL
    // may be used instead.  THIS SHOULD BE CLEANED UP.
    //
    exit_error("MOD_VHF_01_PREP","There is no <s2> object");
  }


  //
  //   Write .mar file and return
  //
  if ((m->action == 10) && (m->marfile != NULL)){
    mod_vhf_make_mpt(m,s,r);
    return;
  }


  //
  //  Fit to data
  //
  fito = onode_child_get_unique(m->o,"fit");
  if (fito != NULL)
    mod_vhf_01_fit(m,s,r,fito);

}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_VHF_01_DONE                              */
/*                                                                           */
/*  Free any storage here when done.                                         */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_01_done(m)
     struct model_struct *m;        // Model parameters
{
  mylog(mylogf,"  MOD_VHF_01_DONE\n");

  // Nothing to do here.
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_VHF_HANDOVER                             */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_handover(r)
     struct response_struct *r; // Response params
{
  int i;
  int ii,jj,xi,yi,ci,gi,tn,s2n,*s2i,rn;
  float tscale,*t,fval;

  tscale = mod_vhf_tscale;
  tn     = mod_vhf_tn;
  s2n    = mod_vhf_c2->s2_xn;

  //
  //  Compute C2 response from the S2 responses
  //
  mod_vhf_c2_r = max_of_2d_farray(mod_vhf_s2_r,s2n,s2n);
  //sprintf(ggstr,"  RespS2_Mid  %f\n",mod_vhf_s2_r[s2n/2][s2n/2]);
  //mylog(mylogf,ggstr);
  //sprintf(ggstr,"  RespC2   %f\n",mod_vhf_c2_r);
  //mylog(mylogf,ggstr);

  // Use convenient local variables to access global parameters


  for(i=0;i<r->n;i++){ // For each response requested

    t = NULL;

    if (r->rformat[i] == 0){  // 'save_as_'
      //
      //  This could be S2, if there is no C2, or the final C2 response
      //

      fval = mod_vhf_c2_r;
    }else if (r->rformat[i] == 1){  // 'save_pop_unit_as_'

      //
      //###                                   s2x_y gi      x,y  ci
      //save_pop_unit_as_c1_011_110  f 1000  c1_0_1_1       1 1   0  response
      //save_pop_unit_as_c1_000_120  f 1000  c1_0_0_0       2 1   0  response
      //save_pop_unit_as_c1_000_130  f 1000  c1_0_0_0       3 1   0  response
      //save_pop_unit_as_c1_000_140  f 1000  c1_0_0_0       4 1   0  response
      //

      s2i = carray_util_parse_underscore(&(r->plname[i][3]),&rn);
      if (rn != 3)
	mylog_exit(mylogf,"MOD_VHF_HANDOVER  c1_x_y_z needs 3 indices");
	
      //ii = jj = 0;
      
      ii = s2i[0];
      jj = s2i[1];
      gi = s2i[2];
      myfree(s2i);

      //printf("ii,jj = %d %d   gi = %d\n",ii,jj,gi);
      //exit(0);

      //printf("PLNAME = %s\n",r->plname[i]);

      //gi = carray_util_underscore_int_final(r->plname[i]);
      xi = r->xi[i];
      yi = r->yi[i];
      ci = r->zi[i];


      //printf("g %d   x,y, %d %d  ci  %d\n",gi,xi,yi,ci);

      fval = mod_vhf_c1_r[ii][jj][gi][xi][yi][ci];
    }

    t = get_const_farray(tn,fval);
    mod_util_resp_store_f_samp(r,i,t,tn,1.0/tscale,mylogf);
    myfree(t);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_VHF_01_GET_RESPONSE                          */
/*                                                                           */
/*  This routine computes the response of the VHF model to the stimlus       */
/*  stored in 's'.  The response is stored in 'r'.  The model has already    */
/*  been configured by a call to 'mod_vhf_01_prep()'.                        */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_01_get_response(m,s,r)
     struct model_struct *m;        // Model parameters
     struct stim_struct *s;         // Stimulus parameters
     struct response_struct *r;     // Response parameters
{
  int i,j,ii,jj;
  int xn,yn,xi,yi,gi,oi,nw,s2n;
  float *w,*x,***stim,u;
  char tstr[SLEN];
  struct vhf_s2_par *s2;

  xn = mod_vhf_xn;
  yn = mod_vhf_yn;

  //
  //  (1) Extract the first frame of the visual stimulus movie.  The first
  //      frame is the only one that is used here.
  //
  stim = f3tensor(1,1,1,xn,1,yn);  // Allocate a 2D array
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      stim[1][i+1][j+1] = s->d[i+1][j+1][1];  // Fill with 1st frame of stim.
    }
  }

  //  Debug:  write out a copy of the visual stimulus.
  //write_2d_data("zz.s0.9.2d",stim[1],1,1,xn,yn,4,2,1,0);


  //
  //  (2) Compute responses for all S1 filters.
  //
  mod_vhf_resp_s1(stim);  // NOTE: 'stim' is modified within this call.

  // Debug.
  if (0){
    float **t2;

    t2 = mod_vhf_chans[0]->r[0][0][0];  // [size][ori][phase][x][y] Response
    write_2d_data("zz.rs.2d",t2,0,0,xn,yn,4,2,1,0);
    exit(0);
  }

  // Debug:  Examine some responses
  // The filter
  //write_2d_data("zz.fi.9.2d",mod_vhf_f[0][0][0],0,0,xn,yn,4,2,1,0);
  // The FT of filter
  //write_2d_data("zz.ft.9.2d",mod_vhf_ft[0][0][0][1],1,1,xn,yn,4,2,1,0);
  // The response
  //write_2d_data("zz.rs.9.2d",mod_vhf_r[0][0][0],0,0,xn,yn,4,2,1,0);


  //
  //  (3) Use S1 responses to compute C1 responses beneath each S2 unit.
  //
  s2n = mod_vhf_c2->s2_xn;
  for(ii=0;ii<s2n;ii++){
    for(jj=0;jj<s2n;jj++){
      mod_vhf_resp_c1_all(mod_vhf_c1_r[ii][jj],ii,jj);
    }
  }


  //
  //  (4) Compute S2 responses from the C1 responses.
  //
  s2 = mod_vhf_s2;
  nw = s2->nw;

  w = get_farray(nw);
  x = get_farray(nw);

  for(j=0;j<s2n;j++){
    jj = (s2n-1) - j;
    for(ii=0;ii<s2n;ii++){
      for(i=0;i<nw;i++){
	w[i] = s2->c1w[i]->w;  // Get the 'i'th C1 weight
	gi   = s2->c1w[i]->gi;
	xi   = s2->c1w[i]->xi;
	yi   = s2->c1w[i]->yi;
	oi   = s2->c1w[i]->ci;
	x[i] = mod_vhf_c1_r[ii][jj][gi][xi][yi][oi];
      }

      mod_vhf_s2_r[ii][jj] = mod_vhf_resp_s2(s2->nltype,s2->s,s2->a,s2->b,
					     s2->eps2,w,x,nw,&u);

      if (0){
	sprintf(tstr,"%f\n",u);
	append_string_to_file("zzz.u.data",tstr);
      }

      // *** WYETH DEBUG REMOVE
      //printf("nw = %d   w[0] = %f\n",nw,w[0]);
      //printf("(ii,jj %d,%d)  x[0] = %f  ==>   %f\n",
      //  ii,jj,x[0],mod_vhf_s2_r[ii][jj]);

      //printf("  %f",mod_vhf_s2_r[ii][jj]);
    }
    //printf("\n");
  }
  myfree(w);
  myfree(x);
  free_f3tensor(stim,1,1,1,xn,1,yn);


  //
  //  (5) Return responses to the caller by storing values within the
  //      response structure 'r'.  For details on 'r', see the definition of
  //      'response_struct' in the file '.../libc/mod.h'
  //
  mod_vhf_handover(r);
}
/**************************************-**************************************/
/*                                                                           */
/*                               MOD_VHF_RUN_01                              */
/*                                                                           */
/*****************************************************************************/
void mod_vhf_run_01(m,s,r,action)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     int action;                // -1-cleanup, 0-prep, 1-run, 2-fit, 10-mar
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // DO THIS FIRST
  numproc = m->process_n;

  mod_vhf_action = action;

  if ((action == 0) || (action == 2) || (action == 10)){
    mod_vhf_01_prep(m,s,r);
  }else if (action == 1){
    mod_vhf_01_get_response(m,s,r);
  }else if (action == -1){
    mod_vhf_01_done(m);
  }
}
