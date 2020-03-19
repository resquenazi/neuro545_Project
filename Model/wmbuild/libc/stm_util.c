/*****************************************************************************/
/*                                                                           */
/*  stm_util.c                                                               */
/*  wyeth bair                                                               */
/*                                                                           */
/*  Fundamental low-level support for stimulus common to                     */
/*   - labr                                                                  */
/*   - wm model                                                              */
/*   - nda analysis                                                          */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2011                                     */
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
#include "stmh.h"

#include "stm_pasushape.h"   // Coordinates for Anitha's shape stimuli

int glob_tot = 0;

/**************************************-**************************************/
/*                                                                           */
/*                               STM_GET_TSAMP                               */
/*                                                                           */
/*  Return the number of raw frames per stimulus pattern frame.              */
/*  For example, if there is a stimulus pattern only every 5th frame, then   */
/*  'tsamp' = 5.                                                             */
/*                                                                           */
/*    For example, tscale 0.002 and stim_samp 1000.0   -->  1                */
/*    For example, tscale 0.001 and stim_samp 1000.0   -->  1                */
/*    For example, tscale 0.002 and stim_samp  100.0   -->  5                */
/*                                                                           */
/*****************************************************************************/
int stm_get_tsamp(tscale,stim_samp)     
     double tscale;     // Temporal sampling (s/frame)    [from .moo file]
     double stim_samp;  // Pattern frames per second      [from .stm file]
{
  int tsamp;

  tsamp = my_rint(1.0/(tscale*stim_samp)); // Time units per stim frame
  if (tsamp <= 0)
    tsamp = 1; // Minimum value for tsamp

  //
  //  WYETH - we should probably complain if 1/tscale is not a multiple of
  //    stim_samp ???
  //

  return tsamp;
}
/**************************************-**************************************/
/*                                                                           */
/*                                STM_GET_TI0                                */
/*                                                                           */
/*  Determine the index of the first frame on which the stimulus pattern     */
/*  will be shown.  This is the closest pattern frame to 'st0'.              */
/*                                                                           */
/*     Raw Frames     - underlying frames                                    */
/*                        occur every 'tscale' seconds                       */
/*     Pattern Frames - frames on which light is displayed                   */
/*                        occur every 'tscale' * 'tsamp' seconds             */
/*                                                                           */
/*****************************************************************************/
int stm_get_ti0(tn,tscale,st0,tsamp)
     int tn;         // Total simulation duration (frames) [.moo]
     double tscale;  // Temporal sampling (s/frame)        [.moo]
     double st0;     // Stimulus pattern start (s)         [.stm]
     int tsamp;      // raw frames per stim pattern frame, see 'stm_get_tsamp'
{
  int t0i,t0_patti;
  double sim_dur_sec;
  double dt_patt_sec;

  sim_dur_sec = tn * tscale;  // Total simulation duration (s)

  if (st0 < 0.0)
    exit_error("STM_GET_TI0","Stimulus start time, 'st0', is negative");

  if (st0 >= sim_dur_sec)
    exit_error("STM_GET_TI0","Stimulus start time exceeds simulation time");

  //  Compute the time between pattern frames (s)
  dt_patt_sec = tsamp * tscale;

  //  Compute the index of the pattern frame closest to 'st0'
  t0_patti = my_rint(st0 / dt_patt_sec);

  //  Compute the raw frame index of the first pattern frame
  t0i = t0_patti * tsamp;

  return t0i;
}
/**************************************-**************************************/
/*                                                                           */
/*                                STM_GET_T1I                                */
/*                                                                           */
/*  Determine the index of the last frame on which the stimulus pattern      */
/*  will be shown.                                                           */
/*                                                                           */
/*     Raw Frames     - underlying frames                                    */
/*                        occur every 'tscale' seconds                       */
/*     Pattern Frames - frames on which light is displayed                   */
/*                        occur every 'tscale' * 'tsamp' seconds             */
/*     Unique Frames  - frames on which the pattern can be altered           */
/*                        occur every 'tscale' * 'tsamp' * 'dwell' seconds   */
/*                                                                           */
/*****************************************************************************/
int stm_get_t1i(tn,tscale,st0,stn,tsamp)
     int tn;         // Total simulation duration (frames) [.moo]
     double tscale;  // Temporal sampling (s/frame)        [.moo]
     double st0;     // Stimulus pattern start (s)         [.stm]
     double stn;     // Stimulus pattern duration (s)      [.stm]
     int tsamp;      // raw frames per stim pattern frame, see 'stm_get_tsamp'
{
  int t0i,t1i,dur_raw_ideal;

  //  Compute the raw frame index of the first pattern frame
  t0i = stm_get_ti0(tn,tscale,st0,tsamp);

  //  Compute the ideal duration of the stimulus in raw frames
  dur_raw_ideal = my_rint(stn / tscale);

  //  Compute the raw frame index where the stimulus ends
  t1i = t0i + dur_raw_ideal;
  if (t1i >= tn)
    t1i = tn - 1;

  return t1i;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 STM_GET_N                                 */
/*                                                                           */
/*  Determine the number of UNIQUE frames that will be required, so that     */
/*  each is shown for 'dwell' frames.  This allows the possibility that      */
/*  the pattern stimulus ends before the last frame during the stimulus      */
/*  epoch or ends before the end of the simulation time.                     */
/*                                                                           */
/*     Raw Frames     - underlying frames                                    */
/*                        occur every 'tscale' seconds                       */
/*     Pattern Frames - frames on which light is displayed                   */
/*                        occur every 'tscale' * 'tsamp' seconds             */
/*     Unique Frames  - frames on which the pattern can be altered           */
/*                        occur every 'tscale' * 'tsamp' * 'dwell' seconds   */
/*                                                                           */
/*****************************************************************************/
int stm_get_n(tn,tscale,st0,stn,dwell,tsamp)
     int tn;         // Total simulation duration (frames) [.moo]
     double tscale;  // Temporal sampling (s/frame)        [.moo]
     double st0;     // Stimulus pattern start (s)         [.stm]
     double stn;     // Stimulus pattern duration (s)      [.stm]
     int dwell;      // Repeat each pattern frame this many times [.stm]
     int tsamp;      // raw frames per stim pattern frame, see 'stm_get_tsamp'
{
  int n;
  int t0i,t1i,dur_raw;

  //  Compute the raw frame index of the first pattern frame
  t0i = stm_get_ti0(tn,tscale,st0,tsamp);

  //  Compute the raw frame index where the stimulus ends
  t1i = stm_get_t1i(tn,tscale,st0,stn,tsamp);

  //  Compute the actual duration in raw frames
  dur_raw = t1i - t0i + 1;

  //  Compute the total number of unique pattern frames
  n = (int)(dur_raw / (tsamp * dwell));

  if (n <= 0)
    exit_error("STIM_GET_N","Number of unique stimulus frames is <= 0");


  //  Note, all 'n' unique pattern frames will appear, and they will all be
  //  shown for 'dwell' frames.  There may be some remaining frames that do
  //  not have a stimulus pattern, because there was not enough time for the
  //  completion of the 'dwell'

  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                          STM_UTIL_COUNT_NEW_FRAMES                        */
/*                                                                           */
/*  Determine how many new frames (with novel patterns) will be created,     */
/*  to estimate number of random draws needed.                               */
/*                                                                           */
/*****************************************************************************/
int stm_util_count_new_frames(zn,tscale,t0,tn,tsamp,dwell)
     int zn;
     float tscale,t0,tn;
     int tsamp;             // Time units per stim frame
     int dwell;
{
  int i;
  int cnt,dn;
  float t;

  /*
    printf("  STM_UTIL_COUNT_NEW_FRAMES\n");
    printf("    zn = %d\n",zn);
    printf("    tscale = %f\n",tscale);
    printf("    tn = %f\n",tn);
    printf("    t0 = %f\n",t0);
    printf("    tsamp = %d\n",tsamp);
    printf("    dwell = %d\n",dwell);
  */

  cnt = 0;
  dn = dwell;
  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;
    if (((i-1)%tsamp) == 0){
      if ((t >= t0) && (t < t0+tn)){
	if (dn == dwell){
	  cnt += 1;
	  dn = 0;
	}
	dn += 1;
      }
    }
  }
  return cnt;
}
/**************************************-**************************************/
/*                                                                           */
/*                              STM_GET_RAN_SEQ                              */
/*                                                                           */
/*  Get several typical types of random sequences.                           */
/*                                                                           */
/*****************************************************************************/
void stm_get_ran_seq(mylogf,distrib,ord,seed,fps,stn,dt,v1,v2,v3,rs,rsn)
     char *mylogf;  // Log file, or NULL
     int distrib;   // 2-Binary, 3-ternary, 7-Mseq
     int ord;       // Order for mseq
     int seed;      // Randomization seed
     float fps;     // Frames per second
     float stn;     // Stimulus duration (sec)
     int dt;        // Dwell time (frames)
     float v1,v2;   // Values for binary distribution
     float v3;      // Middle value for ternary distrib
     float **rs;    // Returned sequence
     int *rsn;      // Returned sequence length (frames)
{
  int i;
  int *seqi,seqn;
  float *seqf;

  /*
    printf("  STM_GET_RAN_SEQ\n");
    printf("    distrib  %d\n",distrib);
    printf("    ord      %d\n",ord);
    printf("    seed     %d\n",seed);
    printf("    fps      %f\n",fps);
    printf("    stn      %f\n",stn);
    printf("    dt       %d\n",dt);
    printf("    v1,2     %f %f\n",v1,v2);
    printf("    v3       %f\n",v3);
  */


  if (distrib == 2){
    seqn = (int)(stn * fps / (float)dt);
    printf("    seqn     %d\n",seqn);
    seqi = myrand_get_std_bin_seq(seqn,seed);
  }else if (distrib == 3){
    seqn = (int)(stn * fps / (float)dt);
    seqi = myrand_get_std_tern_seq(seqn,seed);
  }else if (distrib == 7){
    seqn = 1;
    for(i=0;i<ord;i++)
      seqn *= 2;
    seqi = myrand_get_std_mseq(seqn,seed,ord);
  }else
    mylog_exit(mylogf,"STM_GET_RAN_SEQ  Unknown distrib value");

  seqf = (float *)myalloc(seqn*sizeof(float));
  if (distrib == 3){
    for(i=0;i<seqn;i++){
      if (seqi[i] == -1)   // Store luminance values
	seqf[i] = v1;
      else if (seqi[i] == 1)   // Store luminance values
	seqf[i] = v2;
      else
	seqf[i] = v3;   // seqi[i] is zero
    }
  }else{
    for(i=0;i<seqn;i++){
      if (seqi[i] == 1)   // Store luminance values
	seqf[i] = v1;
      else
	seqf[i] = v2;
    }
  }
  myfree(seqi);

  *rs  = seqf;  // Float sequence
  *rsn = seqn;  // Length of sequence
}
/**************************************-**************************************/
/*                                                                           */
/*                              STM_GET_RAN_SEQ_2D                           */
/*                                                                           */
/*  Returned array is [*rsn][n].                                             */
/*  Thus, there are 'n' random values at each time point.                    */
/*                                                                           */
/*****************************************************************************/
void stm_get_ran_seq_2d(mylogf,n,distrib,seed,fps,stn,dt,v1,v2,v3,rs,rsn)
     char *mylogf;  // Log file, or NULL
     int n;         // Length of second dimension
     int distrib;   // 2-Binary, 3-ternary
     int seed;      // Randomization seed
     float fps;     // Frames per second
     float stn;     // Stimulus duration (sec)
     int dt;        // Dwell time (frames)
     float v1,v2;   // Values for binary distribution
     float v3;      // Middle value for ternary distribution
     float ***rs;   // Returned sequence [rsn][n]
     int *rsn;      // Returned sequence length (frames)
{
  int i,j,k;
  int *seqi,seqn;
  float **seqf;


  // Adding 1 here because value can be rounded too low.
  seqn = 1 + (int)(stn * fps / (float)dt);

  //printf("stn = %f  fps = %f  dt = %d\n",stn,fps,dt);
  //printf("seqn= %d\n",seqn);

  if (distrib == 2){
    seqi = myrand_get_std_bin_seq(n*seqn,seed);
  }else if (distrib == 3){
    seqi = myrand_get_std_tern_seq(n*seqn,seed);
  }else{
    mylog_exit(mylogf,"STM_GET_RAN_SEQ_2D  Unknown distrib value");
  }

  seqf = get_2d_farray(seqn,n);
  if (distrib == 3){
    k = 0;
    for(i=0;i<seqn;i++){
      for(j=0;j<n;j++){
	if (seqi[k] == -1)   // Store luminance values
	  seqf[i][j] = v1;
	else if (seqi[k] == 1)   // Store luminance values
	  seqf[i][j] = v2;
	else
	  seqf[i][j] = v3;   // seqi[i] is zero
	  //seqf[i][j] = 0.0;   // seqi[i] is zero
	k += 1;
      }
    }
  }else{
    k = 0;
    for(i=0;i<seqn;i++){
      for(j=0;j<n;j++){
	if (seqi[k] == 1)   // Store luminance values
	  seqf[i][j] = v1;
	else
	  seqf[i][j] = v2;
	k += 1;
      }
    }
  }
  myfree(seqi);

  *rs  = seqf;  // Float sequence
  *rsn = seqn;  // Length of sequence
}
/**************************************-**************************************/
/*                                                                           */
/*                           STM_SHAPE_FILL_RECURSIVE                        */
/*                                                                           */
/*  If the current point is not filled, fill it and then call the routine    */
/*  on the neighboring 4 pixels.                                             */
/*                                                                           */
/*  ***                                                                      */
/*  2012 Sep 26 - appears to be a limit to recursion of about 130,000 calls  */
/*                                                                           */
/*****************************************************************************/
void stm_shape_fill_recursive(data,xn,yn,i,j,sval)
     float **data;       // Writing shape into this 2D array [xn][yn]
     int xn,yn;          // Size of the 2D drawing frame (pix)
     int i,j;            // Fill from this point
     float sval;         // Luminance value for drawing the shape
{
  // 
  //  If the coordinate (i,j) lies within the 'data' array
  // 

  // WYETH _ DEBUG
  //printf("___LGOB_TOTL = %d  i,j %d %d\n",glob_tot,i,j);
  glob_tot += 1;
  if (glob_tot > 180000){
    //return;
    printf("  *** The shape is too large in number of pixels\n");
    exit_error("STM_SHAPE_FILL_RECURSIVE","Too many recursive calls");
  }
  //printf("___xn,yn  %d %d  i,j  %d %d  sval %f\n",xn,yn,i,j,sval);

  // Fill it, and call the same routine on its 4 neighbors
  data[i][j] = sval;
  if (i+1 < xn)
    if (data[i+1][j] != sval)
      stm_shape_fill_recursive(data,xn,yn,i+1,j,sval);
  if (i > 0)
    if (data[i-1][j] != sval)
      stm_shape_fill_recursive(data,xn,yn,i-1,j,sval);
  if (j+1 < yn)
    if (data[i][j+1] != sval)
      stm_shape_fill_recursive(data,xn,yn,i,j+1,sval);
  if (j > 0)
    if (data[i][j-1] != sval)
      stm_shape_fill_recursive(data,xn,yn,i,j-1,sval);
}
/**************************************-**************************************/
/*                                                                           */
/*                        STM_SHAPE_FILL_RING_RECUR                          */
/*                                                                           */
/*****************************************************************************/
void stm_shape_fill_ring_recur(data,xn,yn,i,j,sval)
     float **data;       // Writing shape into this 2D array [xn][yn]
     int xn,yn;          // Size of the 2D drawing frame (pix)
     int i,j;            // Fill from this point
     float sval;         // Luminance value for drawing the shape
{
  // 
  //  If the coordinate (i,j) lies within the 'data' array
  // 

  // WYETH _ DEBUG
  //printf("___LGOB_TOTL = %d  i,j %d %d\n",glob_tot,i,j);
  glob_tot += 1;
  if (glob_tot > 180000)
    return;
  /*
  if (glob_tot > 600000){
    printf("  *** The shape is too large in number of pixels\n");
    exit_error("STM_SHAPE_FILL_RECURSIVE","Too many recursive calls");
  }
  */
  //printf("___xn,yn  %d %d  i,j  %d %d  sval %f\n",xn,yn,i,j,sval);

  //*** WYETH - UNDER CONSTRUCTION  ***/
  //*** WYETH - UNDER CONSTRUCTION  ***/
  //*** WYETH - UNDER CONSTRUCTION  ***/
  //*** WYETH - UNDER CONSTRUCTION  ***/

  /*

  while(!done){

    // Check whether ring 'ri' is empty
    for(ii=i-ri;ii<=i+ri;ii++)

    if (next_ring_empty)
      fill this ring
      else{
	done = 1;
	call_recurs on this righ
      }
  }

  */

  /*
  if ((i >= 0) && (i < xn)){
    if ((j >= 0) && (j < yn)){
      //
      //  If it has not already been filled
      //
      if (data[i][j] != sval){
  */
	//
	//  Fill it, and call the same routine on its 4 neighbors
	//
	data[i][j] = sval;
	if (i+1 < xn)
	  if (data[i+1][j] != sval)
	    stm_shape_fill_recursive(data,xn,yn,i+1,j,sval);
	if (i > 0)
	  if (data[i-1][j] != sval)
	    stm_shape_fill_recursive(data,xn,yn,i-1,j,sval);
	if (j+1 < yn)
	  if (data[i][j+1] != sval)
	    stm_shape_fill_recursive(data,xn,yn,i,j+1,sval);
	if (j > 0)
	  if (data[i][j-1] != sval)
	    stm_shape_fill_recursive(data,xn,yn,i,j-1,sval);
/*
      }
    }
  }
*/
}
/**************************************-**************************************/
/*                                                                           */
/*                           STM_SHAPE_FILL_NONRECUR                         */
/*                                                                           */
/*  Use a non-recursive algorithm, because recursion seems to be limited in  */
/*  how many times the routine can be called.                                */
/*                                                                           */
/*****************************************************************************/
void stm_shape_fill_nonrecur(data,xn,yn,i,j,sval)
     float **data;       // Writing shape into this 2D array [xn][yn]
     int xn,yn;          // Size of the 2D drawing frame (pix)
     int i,j;            // Fill from this point
     float sval;         // Luminance value for drawing the shape
{
  int done;
  // 
  //  If the coordinate (i,j) lies within the 'data' array
  // 

  // ************** WYETH ************
  // ************** WYETH ************
  // ************** WYETH ************  THIS DOES NOT APPEAR TO BE USED?
  // ************** WYETH ************    Jan 1, 2016
  // ************** WYETH ************

  // WYETH _ DEBUG
  //printf("___LGOB_TOTL = %d  i,j %d %d\n",glob_tot,i,j);
  //glob_tot += 1;
  // printf("___xn,yn  %d %d  i,j  %d %d  sval %f\n",xn,yn,i,j,sval);

  done = 0;
  while(done == 0){

    if ((i >= 0) && (i < xn)){
      if ((j >= 0) && (j < yn)){
	//
	//  If it has not already been filled
	//
	if (data[i][j] != sval){
	  //
	  //  Fill it, and call the same routine on its 4 neighbors
	  //
	  data[i][j] = sval;
	  stm_shape_fill_recursive(data,xn,yn,i+1,j,sval);
	  stm_shape_fill_recursive(data,xn,yn,i-1,j,sval);
	  stm_shape_fill_recursive(data,xn,yn,i,j+1,sval);
	  stm_shape_fill_recursive(data,xn,yn,i,j-1,sval);
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           STM_SHAPE_DRAW_CONTOUR                          */
/*                                                                           */
/*  Given a set of coordinates that define points on a closed shape, draw    */
/*  the pixels that fall on (closest to) the boundary.                       */
/*                                                                           */
/*****************************************************************************/
void stm_shape_draw_contour(data,xn,yn,sx,sy,sn,size,xc,yc,sval)
     float **data;       // Write shape into 2D array [xn][yn]
     int xn,yn;          // Size of frame (pix)
     float *sx,*sy;      // Coordinates of points on shape contour [sn]
     int sn;             // Number of shape points
     float size;         // Size of shape (scaling applied to pix)
     float xc,yc;        // Center of shape in [0..xn-1 , 0..yn-1]
     float sval;         // Luminance value to write for shape
{
  int i;
  int xi,yi;
  float fx,fy;

  for(i=0;i<sn;i++){  // For each point in the contour

    fx = sx[i]*size + xc;  // Float x-coord, after scaling and centering
    fy = sy[i]*size + yc;  // Float y-coord, after scaling and centering

    xi =  my_rint(fx);     // Round to nearest pixel coordinates
    yi =  my_rint(fy);

    // If the rounded coordinates lie in the frame, set this pixel
    if ((xi >= 0) && (xi < xn) && (yi >= 0) && (yi < yn))
      data[xi][yi] = sval;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                      STM_SHAPE_DRAW_CONTOUR_ANTIALIAS                     */
/*                                                                           */
/*  Given a set of coordinates that define points on a closed shape, draw    */
/*  a Gaussian at each location, thereby creating a smooth outline.          */
/*                                                                           */
/*****************************************************************************/
void stm_shape_draw_contour_antialias(data,xn,yn,sx,sy,sn,size,xc,yc,sval,
				      aa_sd)
     float **data;       // Write shape into 2D array [xn][yn]
     int xn,yn;          // Size of frame (pix)
     float *sx,*sy;      // Coordinates of points on shape contour [sn]
     int sn;             // Number of shape points
     float size;         // Size of shape (deg)
     float xc,yc;        // Center of shape in [0..xn-1 , 0..yn-1]
     float sval;         // Value to write for shape
     float aa_sd;        // Standard Deviation of Gaussian for anti-aliasing
{
  int i,j,k;
  int aa_rad,xi,yi,x0,x1,y0,y1;
  float fx,fy,gc,dx,dy,gx,gy,g;

  if ((sval > 0.0) && (sval <= 0.5)){
    exit_error("STM_SHAPE_DRAW_CONTOUR_ANTIALIAS",
	       "Wyeth - this might not work for such 'sval'");
  }

  //
  //  Determine a fixed radius in pixels that will contain any drawing for
  //  anti-aliasing around any given point.  One pixel is the minimum.
  //
  aa_rad = my_rint(aa_sd * 2.0);    // Anti-alias radius (pix)
  if (aa_rad < 1)
    aa_rad = 1;

  gc = -1.0/(2.0 * aa_sd * aa_sd);  // Constant for Gaussian equation

  for(i=0;i<sn;i++){  // For each point on the shape boundary

    fx = sx[i]*size + xc;  // Float x-coord
    fy = sy[i]*size + yc;  // Float y-coord

    xi =  my_rint(fx);  // Round to nearest pixel
    yi =  my_rint(fy);

    //
    //  Determine the loop coordinates:  [x0,x1] and [y0,y1] for drawing
    //
    x0 = xi - aa_rad;
    x1 = xi + aa_rad;
    if (x0 < 0)
      x0 = 0;
    if (x1 >= xn)
      x1 = xn-1;

    y0 = yi - aa_rad;
    y1 = yi + aa_rad;
    if (y0 < 0)
      y0 = 0;
    if (y1 >= yn)
      y1 = yn-1;

    for(j=x0;j<=x1;j++){
      dx = (float)j - fx;  // x-Distance from true point to current pixel
      //printf("dx = %f  gc= %f\n",dx,gc);
      gx = exp(gc * dx*dx);  // Gaussian in x-axis
      for(k=y0;k<=y1;k++){
	dy = (float)k - fy;  // y-Dist. from true point to current pixel
	gy = exp(gc * dy*dy);   // Gaussian in y-axis

	if (sval == 0.0)
	  g = 1.0 - gx * gy;    // xy Gaussian times stimulus value
	else
	  g = gx * gy * sval;    // xy Gaussian times stimulus value
	
	//
	//  WYETH The Following MIGHT BE A PROBLEM w/ Dark on light background
	//
	if (sval > 0.0){
	  if (g > data[j][k]){  // Only if we would be increasing brightness,
	    //printf("data[j]k] = %f  g= %f\n",data[j][k],g);
	    data[j][k] = g;     // thus, use largest value cast from all points
	  }
	}else{
	  if (g < data[j][k]){  // Only if we would be decreasing brightness,
	    data[j][k] = g;     // thus, use largest value cast from all points
	  }
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             STM_SHAPE_DRAW_FILL                           */
/*                                                                           */
/*  Given a set of coordinates that define points on a closed shape, fill    */
/*  the pixels enclosed by the shape with the given value.  Also, if the     */
/*  Standard Deviation for anti-aliasing, 'aa_sd' is > 0, then do anti-      */
/*  aliasing.                                                                */
/*                                                                           */
/*****************************************************************************/
void stm_shape_draw_fill(mylogf,data,xn,yn,sx,sy,sn,size,xc,yc,sval,bgval,
			 x0,y0,aa_sd)
     char *mylogf;       // Log file for printing helpful output, or NULL
     float **data;       // Write shape into this 2D array [xn][yn]
     int xn,yn;          // Size of the frame to be drawn into (pix)
     float *sx,*sy;      // Coordinates of points on shape contour [sn]
     int sn;             // Number of shape points
     float size;         // Size scaling of raw shape coordinates
     float xc,yc;        // Center of shape in [0..xn-1 , 0..yn-1]
     float sval;         // Luminance value to write for shape
     float bgval;        // Background value
     float x0,y0;        // One point known to be on inside of raw shape
     float aa_sd;        // SD for anti-aliasing (0.0 for none)
{
  int i,j;
  int xin,yin,cnt,si,sj,xci,yci;
  float txc,tyc,**td;

  // Determine the inside point after scaling and translating
  xin = my_rint(x0*size + xc);
  yin = my_rint(y0*size + yc);

  if ((yin < 0) || (xin < 0)  || (xin >= xn) ||  (yin >= yn)){
    // ****
    // ****
    // ****  WYETH - This is a hack - we need to figure out an interior
    // ****    point that is within bounds.
    // ****
    //printf("  *** Interior point:  (x,y) = %d %d\n",xin,yin);
    //exit_error("STM_SHAPE_DRAW_FILL","Interior point is out of bounds");

    //
    //  *** WYETH HERE new strategy:   draw into a temporary array the shape
    //       centered at 0,0, then copy this to the desired array
    //      Limitations:  be careful about pixel rounding issues?
    //

    // *** WYETH - note, these are relative to about xn/2 and yn/2
    //xc = (float)(xn-1)/2.0 + cx/sscale;  // Center w.r.t. 0..xn-1 (pix)
    //yc = (float)(yn-1)/2.0 + cy/sscale;  // Center w.r.t. 0..yn-1 (pix)

    xci = (int)(xc - (xn-1)/2.0);          //  E.g.,  5 = (int)5.5;
    yci = (int)(yc - (yn-1)/2.0);

    txc = xc - (float)xci;  //  E.g.,  0.5 = 5.5 - 5.0
    tyc = yc - (float)yci;

    //printf("  *** Interior point:  (x0,y0) = %f %f\n",x0,y0);

    xin = my_rint(x0*size + txc);
    yin = my_rint(y0*size + tyc);

    td = get_const_2d_farray(xn,yn,bgval);

    // Draw hard edge
    stm_shape_draw_contour(td,xn,yn,sx,sy,sn,size,txc,tyc,sval);
    // Fill it

    //printf("  *** Interior point:  (x,y) = %d %d\n",xin,yin);
    glob_tot = 0;
    stm_shape_fill_recursive(td,xn,yn,xin,yin,sval);

    // Copy 'td' into 'data', applying appropriate pixel shift
    for(i=0;i<xn;i++){
      si = i + xci;
      if ((si >= 0) && (si < xn)){
	for(j=0;j<yn;j++){
	  sj = j + yci;
	  if ((sj >= 0) && (sj < yn)){
	    data[si][sj] = td[i][j];
	  }
	}
      }
    }

    free_2d_farray(td,xn);

  }else{

    // Draw in the hard edge
    stm_shape_draw_contour(data,xn,yn,sx,sy,sn,size,xc,yc,sval);
    //printf("sn = %d   %f %f   %f %f\n",sn,sx[0],sy[0],sx[sn-1],sy[sn-1]);

    //printf("HERE  xn,yn %d %d  xin,yin %d %d   sval %f\n",xn,yn,xin,yin,sval);

    // Call recursive fill routine on one point known to be inside
    glob_tot = 0;
    stm_shape_fill_recursive(data,xn,yn,xin,yin,sval);
    //printf("_________GLOB_TOT = %d\n",glob_tot);

  }

  //
  // Fill any isolated pixels.  I believe that there could be some points
  // that end up being isolated and unfilled.  This tries to correct any
  // such points by filling them in.
  //
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      if (data[i][j] == bgval){
	cnt = 0;
	if (i > 0)
	  if (data[i-1][j] == sval)
	    cnt += 1;
	if (i < (xn-1))
	  if (data[i+1][j] == sval)
	    cnt += 1;
	if (j > 0)
	  if (data[i][j-1] == sval)
	    cnt += 1;
	if (j < (yn-1))
	  if (data[i][j+1] == sval)
	    cnt += 1;
	if (cnt == 4)
	  data[i][j] = sval;  // Fill, if all neighbors are on
      }
    }
  }

  if (aa_sd > 0.0){  // If we are using anti-aliasing
    //
    //  WYETH - it would be more sophisticated to remove the contour and
    //  then redraw an anti-aliased contour, but this is a bit harder, because
    //  some pixels on the inside edge of the contour are lower.  I would need
    //  special checking for this.  Thus, I have decided to do the simpler
    //  thing below.  This causes the shape to be a bit larger than expected.
    //
    // Remove the contour
    //stm_shape_draw_contour(data,xn,yn,sx,sy,sn,size,xc,yc,bgval);

    // Redraw the antialiased contour
    stm_shape_draw_contour_antialias(data,xn,yn,sx,sy,sn,size,xc,yc,sval,
				     aa_sd);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             STM_NOISE_GET_GRID                            */
/*                                                                           */
/*  Fill the frame with a grid of boxes having 'gdata' values.               */
/*                                                                           */
/*****************************************************************************/
float **stm_noise_get_grid(xn,yn,gw,gh,pixw,pixh,x0,y0,gdata,bgval)
     int xn,yn;     // Size of frame (pix)
     int gw,gh;     // Grid dimensions (boxes)
     int pixw,pixh; // Box size (pix)
     int x0,y0;     // lower left corner
     float **gdata; // values for each box [gw][gh]
     float bgval;
{
  int i,j;
  int i0,i1,j0,j1,gdi,gdj,gdi0,gdj0,bx0,by0,iwrap,jwrap,jwrap0;
  float **g,gval;

  /*
  printf("xn,yn  %d, %d\n",xn,yn);
  printf("gw,gh  %d, %d\n",gw,gh);
  printf("pixw,y  %d, %d\n",pixw,pixh);
  printf("x0,y0  %d, %d\n",x0,y0);
  printf("bgval  %f\n",bgval);
  */

  g = get_const_2d_farray(xn,yn,bgval);

  i0 = x0;
  i1 = x0 + gw*pixw - 1;
  if (i0 < 0){
    i0 = 0;
    iwrap = (i0-x0) % pixw;  // Remainder
    gdi0  = (i0-x0) / pixw;
  }else{
    gdi0 = 0;
    iwrap = 0;
  }

  if (i1 >= xn)
    i1 = xn-1;

  j0 = y0;
  j1 = y0 + gh*pixh - 1;
  if (j0 < 0){
    j0 = 0;
    jwrap0 = (j0-y0) % pixh;
    gdj0   = (j0-y0) / pixh;
  }else{
    gdj0 = 0;
    jwrap0 = 0;
  }

  if (j1 >= yn)
    j1 = yn-1;

  /*
  printf("i:  %d - %d\n",i0,i1);
  printf("j:  %d - %d\n",j0,j1);
  printf("gdi0  %d\n",gdi0);
  printf("gdj0  %d\n",gdj0);
  printf("iwrap  %d\n",iwrap);
  printf("jwrap0  %d\n",jwrap0);
  */

  gdi = gdi0;
  for(i=i0;i<=i1;i++){
    gdj = gdj0;
    gval = gdata[gdi][gdj];
    jwrap = jwrap0;
    for(j=j0;j<=j1;j++){
      g[i][j] = gval;

      jwrap += 1;
      //if (jwrap == pixh){  OLD, WRONG
      if ((jwrap == pixh) && (j < j1)){  // Don't enter on last rep
	jwrap = 0;
	gdj += 1;  // Move to next index in 'gdata'
	gval = gdata[gdi][gdj];
      }

    }
    iwrap += 1;
    if (iwrap == pixw){
      iwrap = 0;
      gdi += 1;  // Move to next index in 'gdata'
    }
  }

  return g;
}
/**************************************-**************************************/
/*                                                                           */
/*                             STM_VANHAT_GET_PATCH                          */
/*                                                                           */
/*  Fill the frame with a grid of boxes having 'gdata' values.               */
/*                                                                           */
/*****************************************************************************/
float **stm_vanhat_get_patch(seed,xn,yn,imagedir,normtype)
     int seed;        // Random seed to pick image
     int xn,yn;       // Size of image patch
     char *imagedir;  // Directory where images are located
     char *normtype;  // "-1 to 1"
{
  int i,j;
  int image_n,image_w,image_h,xi,yi,ii,**idata,*idt,done;
  float **fdata,dmin,dmax,fc,*fdt,fval;

  image_n = 100;
  image_w = 1536;
  image_h = 1024;

  //
  //  Select coordinates and image randomly
  //
  if (seed > 0)
    seed = -seed;

  fdata = get_2d_farray(xn,yn);

  done = 0;
  while(done == 0){  // Until we find a non-constant image patch
    xi = (int)((image_w - xn + 1) * myrand_util_ran2(&seed));  // x-index
    yi = (int)((image_h - yn + 1) * myrand_util_ran2(&seed));  // y-index
    ii = 1 + (int)(image_n * myrand_util_ran2(&seed));         // image index

    idata = read_vanhateren(imagedir,ii,xi,yi,xn,yn);

    //
    //  Convert to float and find min, max at same time
    //
    dmax = dmin = (float)idata[0][0];
    for(i=0;i<xn;i++){
      idt = idata[i];
      fdt = fdata[i];
      for(j=0;j<yn;j++){
	fdt[j] = fval = (float)idt[j];
	if (fval > dmax)
	  dmax = fval;
	else if (fval < dmin)
	  dmin = fval;
      }
    }
    free_2d_iarray(idata,xn);

    if (dmax > dmin)
      done = 1;
    else{
      ; //printf("===> CONSTANT IMAGE, val = %f\n, try again...",dmin);
    }
  }

  if (strcmp(normtype,"-1 to 1")==0){
    fc = 2.0 / (dmax - dmin);
    for(i=0;i<xn;i++){
      fdt = fdata[i];
      for(j=0;j<yn;j++){
	fdt[j] = -1.0 + fc * (fdt[j] - dmin);
      }
    }
  }else{
    exit_error("STM_VANHAT_GET_PATCH","Unknown 'normtype'");
  }

  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                            STM_RANSTEP_GET_STIM                           */
/*                                                                           */
/*  Return jump sequence in terms of fraction of cycle to jump.              */
/*                                                                           */
/*****************************************************************************/
float *stm_ranstep_get_stim(mylogf,sf,jsz,jzero,jn,seed,stn,fps,dt,iflag,rseqn)
     char *mylogf;  // Log file, or NULL
     float sf;      // spatial frequency of sine wave
     float jsz;     // 'step_size' from param file
     int jzero;     // 'step_zero' from param file
     int jn;        // 'step_n' from param file
     int seed;      // randomization seed
     float stn;     // stimulus duration (sec)
     float fps;     // video frames per second
     int dt;        // 'dt' from param file, number of frames to hold step
     int iflag;     // 0-return step size, 1-return stimulus INDEX VALUE
     int *rseqn;    // return length of sequence
{
  int i;
  int nstep,seqn;
  int *seqi;
  float jfrac,*seqf;

  //printf("sf = %f   jsz = %f\n",sf,jsz);

  // Get number of unique jump sizes 'nstep' and unit jump size 'jfrac'
  if (jsz > 0.0){
    jfrac = jsz * sf;      // Fraction of cycle to jump
    if (jzero == 1){
      nstep = jn*2 + 1;    // Use zero
    }else{
      nstep = jn*2;        // Don't use zero
    }
  }else{  // If 'jsz' = 0, use evenly spaced positions in 0...360
    jsz = 0.0;  // Just to be sure it isn't negative
    jfrac = 1.0/(float)jn;  // Not needed, not used, below
    nstep = jn;
  }

  seqn = (int)(stn * fps / (float)dt);

  //printf("  stn = %f\n",stn);
  //printf("  fps = %f\n",fps);
  //printf("  dt  = %d\n",dt);
  

  seqi = myrand_get_std_unif_int_seq(seqn,seed,nstep); // 0..nstep -1
 
  if (iflag == 1)
    seqf = i2farray(seqi,seqn);
  else{
    seqf = (float *)myalloc(seqn*sizeof(float));
    if (jsz == 0.0){
      for(i=0;i<seqn;i++)
	seqf[i] = (float)seqi[i]/(float)jn;        // 0, 1/jn, ... (jn-1)/jn
    }else if (jzero == 1){
      for(i=0;i<seqn;i++)
	seqf[i] = jfrac * (float)(seqi[i] - jn);   // jf*[-jn ...0... +jn]
    }else{
      for(i=0;i<seqn;i++){
	if (seqi[i] < jn)
	  seqf[i] = jfrac * (float)(seqi[i]-jn);   // jf*[-jn ... -1]
	else
	  seqf[i] = jfrac * (float)(seqi[i]+1-jn); // jf*[   1 ... +jn]
      }
    }
  }

  *rseqn = seqn;
  return seqf;
}
/**************************************-**************************************/
/*                                                                           */
/*                          STM_BARNOISE_GET_TEMPLATE                        */
/*                                                                           */
/*  Return a template for a set of contiguous bars.  The value at each       */
/*  pixel ranges from 0 (background) to 1 (first bar) ... to n (last bar)    */
/*  to n+1 (background, again).                                              */
/*                                                                           */
/*  Values that are not exact integers indicate that the color is on the     */
/*  smoothed boundary between bars, and should be a mixture of the two       */
/*  bar values.                                                              */
/*                                                                           */
/*  No "smoothing" of values is applied at the long ends of the bars.        */
/*                                                                           */
/*****************************************************************************/
float **stm_barnoise_get_template(mylogf,xn,yn,nbar,w,h,theta,bord,xc,yc,hoff)
     char *mylogf;  // Log file, or NULL
     int xn,yn;     // Size of template (pix)
     int nbar;      // Number of bars
     float w,h;     // Width and height (length) of bars (pix)
     float theta;   // Orientation of bars
     float bord;    // Width of border to smooth between bars
     float xc,yc;   // Desired center of barfield (pix)
     float hoff;    // Horizontal offset
{
  int i,j;
  int txi;
  float **t,a,b,c,d,thetar,xctrf,yctrf,xx,yy,x,y,ymin,ymax;
  double x0,dx,tx,marg,xtra,bord_frac;

  //mylog(mylogf,"  STM_BARNOISE_GET_TEMPLATE\n");

  while(theta >= 360.0)
    theta -= 360.0;
  while(theta < 0.0)
    theta += 360.0;

  thetar = -theta * M_PI/180.0;
  a =  cos(thetar); b = sin(thetar); // rotation matrix
  c = -sin(thetar); d = cos(thetar);

  t = get_2d_farray(xn,yn);

  xctrf = ((float)xn-1.0)/2.0;  // Exact float center of field
  yctrf = ((float)yn-1.0)/2.0;

  ymin = -(float)h/2.0;  // Upper extent of bar
  ymax =  (float)h/2.0;  // Lower extent of bar

  dx = (double)w;
  x0 = -(double)(nbar+2)*dx/2.0;

  /*printf("    ymin = %f\n",ymin);
    printf("    ymax = %f\n",ymax);
    printf("    xctrf = %f\n",xctrf);
    printf("    yctrf = %f\n",yctrf);
    printf("    nbar = %d\n",nbar);
    printf("    bar_wid = %f\n",w);
    printf("    bar_len = %f\n",h);
    printf("    border = %f\n",bord);  */

  bord_frac = bord/w;

  for(i=0;i<xn;i++){
    xx = (float)i - xctrf;   // relative to center of field
    xx -= xc;
    for(j=0;j<yn;j++){
      yy = (float)j - yctrf; // relative to center of field
      yy -= yc;

      //
      //  (x,y) position lies on a canonical set of vertical bars.
      //
      x = a*xx + c*yy;  // rotate back
      y = b*xx + d*yy;

      x -= hoff;  // WYETH - to allow movement of bars orthogonal to ori

      if ((y > ymax) || (y < ymin)){  // Above or below bars
	t[i][j] = 0.0;  // Use background color
      }else{
	tx = (x - x0)/dx;  // Distance in 'bars' from 0th bar start
	txi = (int)tx;     // integer bar index
	marg = tx - (double)txi;  // marginal amount
	if (marg < bord_frac)
	  // Slightly into new bar
	  xtra = -0.5 * (1.0 - marg/bord_frac);  // Amount to add to index
	else if (marg > (1.0-bord_frac))
	  // Almost out of current bar
	  xtra = 0.5 * (marg - (1.0-bord_frac))/bord_frac;
	else
	  xtra = 0.0;  // We're in middle of bar, no extra to add.

	/*
	if (j == yn/2)
	  printf("i = %d  xx %.2f  x %.2f  tx %f  txi %d  xtra %f\n",i,
	  xx,x,tx,txi,xtra);*/

	if (tx < 0.5){  // Don't allow smoothing beyond half-bar beyond end
	  t[i][j] = 0.0;  // background
	}else if (tx >= ((double)nbar+1.5)){
	  t[i][j] = 0.0;  // background
	}else{
	  t[i][j] = (float)txi + (float)xtra;
	}
      }
    }
  }
  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                            STM_UTIL_GRID_COORDS                           */
/*                                                                           */
/*  Return the coordinates of a rectangular grid tilted by theta.            */
/*                                                                           */
/*****************************************************************************/
void stm_util_grid_coords(mylogf,nx,ny,theta,dx,dy,rcx,rcy)
     char *mylogf;    // Log file, or NULL
     int   nx;        // Number of elements horizontally
     int   ny;        // Number of elements vertically
     float theta;     // Direction
     float dx;        // Horizontal offset
     float dy;        // Vertical offset
     float ***rcx,***rcy; // [nbar][npos] grid coords w.r.t. (0,0)
{
  int i,j;
  float a,b,c,d,x0,y0,x,y,**cx,**cy;

  // Allow theta to be 0..360 to show mirror image using 0,180
  while (theta >= 360.0)
    theta -= 360.0;
  while (theta < 0.0)
    theta += 360.0;

  // rotation matrix
  a =  cos(M_PI/180.0*theta); b = sin(M_PI/180.0*theta);
  c = -sin(M_PI/180.0*theta); d = cos(M_PI/180.0*theta);

  cx = get_2d_farray(nx,ny);
  cy = get_2d_farray(nx,ny);

  // Coords for bar at lower left of array, centered on 0,0 in pixels
  x0 = -((float)(nx-1) * dx)/2.0;
  y0 = -((float)(ny-1) * dy)/2.0;

  // Get bar coords before rotation, centered on 0,0, in pixels
  x = x0;
  y = y0;
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){

      // Rotate
      cx[i][j] = a*x + c*y;
      cy[i][j] = b*x + d*y;

      y += dy;
    }
    x += dx;
    y = y0;
  }

  *rcx = cx;
  *rcy = cy;
}
/**************************************-**************************************/
/*                                                                           */
/*                           STM_UTIL_GRID_COORDS_1D                         */
/*                                                                           */
/*  Return the coordinates of a rectangular grid tilted by theta.            */
/*                                                                           */
/*****************************************************************************/
void stm_util_grid_coords_1d(mylogf,nx,ny,theta,dx,dy,rcx,rcy)
     char *mylogf;      // Log file, or NULL
     int   nx;          // Number of elements horizontally
     int   ny;          // Number of elements vertically
     float theta;       // Direction
     float dx;          // Horizontal offset
     float dy;          // Vertical offset
     float **rcx,**rcy; // [nx*ny] grid coords w.r.t. (0,0)
{
  int i,j,k;
  float **cx,**cy,*cx1,*cy1;

  stm_util_grid_coords(mylogf,nx,ny,theta,dx,dy,&cx,&cy);

  cx1 = (float *)myalloc(nx*ny*sizeof(float));
  cy1 = (float *)myalloc(nx*ny*sizeof(float));

  k = 0;
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      cx1[k] = cx[i][j];
      cy1[k] = cy[i][j];
      k += 1;
    }
  }

  free_2d_farray(cx,nx);
  free_2d_farray(cy,nx);

  *rcx = cx1;
  *rcy = cy1;
}
/**************************************-**************************************/
/*                                                                           */
/*                          STM_BARRAY_GET_BAR_CENTERS                       */
/*                                                                           */
/*  Return the coordinates of the bar centers for all bars in the grid.      */
/*                                                                           */
/*****************************************************************************/
void stm_barray_get_bar_centers(mylogf,nbar,npos,theta,barw,barl,bar_gap,
				 pos_gap,rcx,rcy)
     char *mylogf;  // Log file, or NULL
     int nbar;      // Number of bars
     int npos;      // Number of positions
     float theta;   // Direction
     float barw;    // Bar width (deg)
     float barl;    // Bar length (deg)
     float bar_gap; // Bar gap (deg)
     float pos_gap; // Pos gap (deg)
     float ***rcx,***rcy; // [nbar][npos] coords (deg) for bar centers wrt 0,0
{

  theta += 90.0;  // To match 'direction' convention
  stm_util_grid_coords(mylogf,nbar,npos,theta,barl+bar_gap,barw+pos_gap,
		       rcx,rcy);
}
/**************************************-**************************************/
/*                                                                           */
/*                             STM_BARRAY_GET_SEQ                            */
/*                                                                           */
/*  The returned array is [n][nbar] and contains an integer to indicate the  */
/*  state of each bar.                                                       */
/*                                                                           */
/*  Bar positions are numbered from 1..npos from bottom to top of array      */
/*  and 0 indicates no bar.                                                  */
/*                                                                           */
/*****************************************************************************/
int **stm_barray_get_seq(mylogf,n,nbar,npos,seqflag,seed)
     char *mylogf;  // Log file, or NULL
     int n;         // Number of bar patterns in sequence
     int nbar;      // Number of bars
     int npos;      // Number of positions
     int seqflag;   // 0-regular, 1-one bar at a time, 99-draw all bars
     int seed;      // for bar sequence
{
  int i,j,k;
  int **seqi;
  float fran;

  seqi = get_zero_2d_iarray(n,nbar);

  if (seed > 0)
    seed = -seed;

  // Get random sequence
  if (seqflag == 0){  // Pick each bar independently
    fran = (float)(npos + 1);

    for(i=0;i<n;i++){
      for(j=0;j<nbar;j++){
	seqi[i][j] = (int)(fran * myrand_util_ran2(&seed));
      }
    }
  }else if (seqflag == 1){  // Pick one bar, or none
    fran = (float)(npos * nbar + 1);
    for(i=0;i<n;i++){
      k = (int)(fran * myrand_util_ran2(&seed));
      j = k/npos;
      if (j < nbar){  // Change this bar to a non-zero position value
	seqi[i][j] = (k % npos) + 1;
      }
    }
  }else if (seqflag == 2){  // Pick ends to be the same

    fran = (float)(npos + 1);

    for(i=0;i<n;i++){
      for(j=0;j<(nbar-1);j++){
	seqi[i][j] = (int)(fran * myrand_util_ran2(&seed));
      }
      seqi[i][nbar-1] = seqi[i][0];
    }

  }else if (seqflag == 99){  // Put '1' everywhere
    for(i=0;i<n;i++){
      for(j=0;j<nbar;j++){
	seqi[i][j] = 1;
      }
    }
  }else{
    mylog_exit(mylogf,"STM_BARRAY_GET_SEQ  Unknown seqflag");
  }

  return seqi;
}
/**************************************-**************************************/
/*                                                                           */
/*                             STM_STAT1D_GET_SEQ                            */
/*                                                                           */
/*  Return the list of ON and OFF durations in 'tscale' time units.          */
/*  Return the index list for blockwise randomization.                       */
/*                                                                           */
/*****************************************************************************/
void stm_stat1d_get_seq(mylogf,nblock,seqflag,seed,tscale,rndur,rton,rtoff,rx)
     char *mylogf;  // Log file, or NULL
     int nblock;    // Number of blocks
     int seqflag;   // 0: 1s...32s, 
     int seed;      // for bar sequence
     float tscale;  // time unit (sec)
     int *rndur;    // Number of stimuli (number of durations)
     int **rton;    // ON durations, tscale units [ndur]
     int **rtoff;   // OFF durations, tscale units [ndur]
     int ***rx;     // Indices [nblock][ndur]
{
  int i;
  int ndur,*ton,*toff,*seedlist,**indx;

  mylog(mylogf,"  STM_STAT1D_GET_SEQ\n");

  if (seqflag == 0){
    ndur = 6;
    ton = (int *)myalloc(ndur*sizeof(int));
    toff = (int *)myalloc(ndur*sizeof(int));
    ton[0] = my_rint( 1.0/tscale);  toff[0] = my_rint(16.0/tscale);
    ton[1] = my_rint( 2.0/tscale);  toff[1] = my_rint( 2.0/tscale);
    ton[2] = my_rint( 4.0/tscale);  toff[2] = my_rint( 4.0/tscale);
    ton[3] = my_rint( 8.0/tscale);  toff[3] = my_rint( 8.0/tscale);
    ton[4] = my_rint(16.0/tscale);  toff[4] = my_rint( 8.0/tscale);
    ton[5] = my_rint(32.0/tscale);  toff[5] = my_rint(16.0/tscale); // ndur-1
  }else{
    mylog_exit(mylogf,"STM_STAT1D_GET_SEQ  Unknown seqflag");
  }

  seedlist = get_seeds(seed,1000000,nblock);

  indx = (int **)myalloc(nblock*sizeof(int *));
  for(i=0;i<nblock;i++)
    indx[i] = get_shuffle_index(ndur,3,seedlist[i]);

  myfree(seedlist);

  *rndur = ndur;
  *rton  = ton;
  *rtoff = toff;
  *rx    = indx;
}
/**************************************-**************************************/
/*                                                                           */
/*                          STM_DIRMOD_GET_TEST_PAT                          */
/*                                                                           */
/*  Return the set of patterns [n][pdur] to be added to the sequence.        */
/*                                                                           */
/*  PATTERN TYPES "test_pat"                                                 */
/*    0 - alternate between all P and all A                                  */
/*    1 - alternate between all P and 3A,P...                                */
/*                                                                           */
/*****************************************************************************/
float **stm_dirmod_get_test_pat(mylogf,n,pdur,t_etf,t_pat,n_anti,fps,sflag,
				degflag)
     char *mylogf;
     int n;          // Number of patterns
     int pdur;       // pattern duration (frames)
     float t_etf;    // TF for jumps
     int t_pat;      // Pattern type
     int n_anti;     // Number of anti frames
     float fps;      // Frames per second
     int sflag;      // 0-sign for labr, 1-sign for model
     int degflag;    // 0-use fraction of a cycle,  1-use 0..360 degrees
{
  int i,j;
  float **pat,x,frac;
  char ggstr[SLEN];

  //sprintf(ggstr,"    pdur =   %d\n",pdur); mylog(mylogf,ggstr);
  //sprintf(ggstr,"    n    =   %d\n",n); mylog(mylogf,ggstr);
  //sprintf(ggstr,"    t_etf = %.2f\n",t_etf); mylog(mylogf,ggstr);
  //sprintf(ggstr,"    t_pat =   %d\n",t_pat); mylog(mylogf,ggstr);

  frac = t_etf/fps;  // Fraction of cycle to jump
  if (degflag)
    frac *= 360.0;   // Degrees

  // NOTE:  -frac IS PREF JUMP for LABR STIMULI
  if (sflag == 1)
    frac *= -1.0;
  
  pat = get_const_2d_farray(n,pdur,-frac);  // Start w/ all pref jumps
  
 if (t_pat == 0){
    for(i=0;i<n;i++){
      if (i%2 == 0){
	;                       // Let it be all 'P' jumps
      }else{
	for(j=0;j<pdur;j++)
	  pat[i][j] = frac;     // All 'A' jumps
      }
    }
  }else if (t_pat == 1){
    for(i=0;i<n;i++){
      if (i%2 == 0){
	;                       // Let it be all 'P' jumps
      }else{
	for(j=0;j<n_anti;j++){
	  if (j < pdur){
	    pat[i][j] = frac;   // Start w/ n_anti 'A' jumps
	  }
	}
      }
    }
  }else{
    mylog_exit(mylogf,"STM_DIRMOD_GET_TEST_PAT  Unknown pattern type.\n");
  }

  return pat;
}
/**************************************-**************************************/
/*                                                                           */
/*                             STM_DIRMOD_TEST_SEQ                           */
/*                                                                           */
/*  For 'dirmod_test' stimulus.  Add test patterns into the sequence.        */
/*                                                                           */
/*****************************************************************************/
void stm_dirmod_test_seq(mylogf,t_dur,t_int,t_etf,t_pat,n_anti,fps,sflag,
			 degflag,rseqf,rseqn)
     char *mylogf;
     int t_dur;      // Test duration (frames)
     int t_int;      // interval between test patterns (frames)
     float t_etf;    // TF for jumps
     int t_pat;      // Pattern type
     int n_anti;     // Number of anti frames
     float fps;      // Frames per second
     int sflag;      // 0-sign for labr, 1-sign for model
     int degflag;    // 0-frac of cycle, 1-degrees
     float **rseqf;
     int *rseqn;
{
  int i,j,k;
  int n,nt,nn,ti;
  float *seqf,**pat;
  char ggstr[SLEN];

  //sprintf(ggstr,"    t_dur =   %d\n",t_dur); mylog(mylogf,ggstr);
  //sprintf(ggstr,"    t_int =   %d\n",t_int); mylog(mylogf,ggstr);

  n = *rseqn;        /* Original sequence length */
  nt = n / t_int;    /* Number of test sequences to add */
  nn = n + nt*t_dur; /* New number of total frames */

  // Get patterns to insert
  pat = stm_dirmod_get_test_pat(mylogf,nt,t_dur,t_etf,t_pat,n_anti,fps,sflag,
				degflag);

  seqf = (float *)myalloc(nn*sizeof(float));

  //sprintf(ggstr," nt=%d\n",nt);  mylog(mylogf,ggstr);
  //sprintf(ggstr," n=%d\n",n);  mylog(mylogf,ggstr);
  
  k = 0;
  for(i=0;i<nt;i++){
    for(j=0;j<t_int;j++){
      ti = i*t_int + j;
      if (ti >= n)
	mylog_exit(mylogf,"STM_DIRMOD_TEST_SEQ  Should not happen\n");
      seqf[k] = (*rseqf)[ti];   // Copy original sequence
      k += 1;
    }
    for(j=0;j<t_dur;j++){       // Insert pattern
      seqf[k] = pat[i][j];
      k += 1;
    }
  }

  ti += 1;
  while (ti < n){ // Copy the end of the original sequence
    seqf[k] = (*rseqf)[ti];
    ti += 1;
    k += 1;
    if (k > nn)
      mylog_exit(mylogf,"STM_DIRMOD_TEST_SEQ  Should not happen for k\n");
  }
  free_2d_farray(pat,nt);

  myfree(*rseqf); // Free old pattern
  *rseqf = seqf;
  *rseqn = nn;
}
/**************************************-**************************************/
/*                                                                           */
/*                             STM_WY4ST_GET_SEQ                             */
/*                                                                           */
/*  The intended meaning of the values in the returned state sequence is     */
/*                                                                           */
/*        cs                                                                 */
/*    0 - 00                                                                 */
/*    1 - 01                                                                 */
/*    2 - 10                                                                 */
/*    3 - 11  where '1' indicates ori parallel to CRF preferred              */
/*                                                                           */
/*****************************************************************************/
void stm_wy4st_get_seq(mylogf,n,seed,ststr,rst)
     char *mylogf;
     int n;         // Number of states in sequence
     int seed;      // Randomization seed
     char *ststr;   // state string, a sequence of 0,1,2,3
     int **rst;     // Returned list ints, 0..3 [n]
{
  int i,k;
  int *st,bit_s,bit_c;
  char ggstr[SLEN];

  //sprintf(ggstr,"  STM_WY4ST_GET_SEQ\n");
  //mylog(mylogf,ggstr);

  //printf("    seed = %d\n",seed);
  //printf("    n = %d\n",n);
  //printf("    state string = %s\n",ststr);

  st = get_zero_iarray(n);

  if (ststr != NULL){ // Use the state string
    
    if (strlen(ststr) != n)
      mylog_exit(mylogf,"STM_WY4ST_GET_SEQ - state string length is not n");

    for(i=0;i<n;i++){
      st[i] = (int)(ststr[i]) - (int)'0';
    }
  }else if ((seed >= 0) && (seed < 256)){  // 'seed' used to encode sequence
    //
    //  Bits encode C+S state, for example with n=4:
    //
    //            cs cs cs cs
    //    152  =  10 01 10 00
    //
    for(i=0;i<n;i++){
      bit_s = (seed >> (i*2  )) % 2;
      bit_c = (seed >> (i*2+1)) % 2;
      
      //printf("bits C S = %d %d\n",bit_c,bit_s);

      k = n-1-i;

      if (bit_c == 0){
	if (bit_s == 0)
	  st[k] = 0;      //  Orth center, Orth surr
	else
	  st[k] = 1;      //  Orth center, Iso surr
      }else{
	if (bit_s == 0)
	  st[k] = 2;      //  Center w/ Orth surround
	else
	  st[k] = 3;      //  Center w/ Iso surround
      }
    }
  }else{
    st = myrand_get_std_unif_int_seq(n,seed,4);
  }

  // WYETH HERE DEBUG - print out state list
  //for(i=0;i<n;i++)
  //printf("  %d state  %d\n",i,st[i]);

  *rst = st;
}
/**************************************-**************************************/
/*                                                                           */
/*                               STM_GPATCH_GET                              */
/*                                                                           */
/*****************************************************************************/
struct stmh_gpatch *stm_gpatch_get(x0,y0,dir,sf,tf,sdp,sdo,ph0,amp)
     float x0,y0;
     float dir;
     float sf;
     float tf;
     float sdp;
     float sdo;
     float ph0;
     float amp;
{
  struct stmh_gpatch *gp;

  gp = (struct stmh_gpatch *)myalloc(sizeof(struct stmh_gpatch));

  gp->ptype   = 1;
  gp->x0      = x0;
  gp->y0      = y0;
  gp->dir     = dir;
  gp->sf      = sf;
  gp->tf      = tf;
  gp->sd_par  = sdp;
  gp->sd_orth = sdo;
  gp->ph0     = ph0;
  gp->amp     = amp;
  gp->seed    = -1;
  gp->pixw    = -1;

  return gp;
}
/**************************************-**************************************/
/*                                                                           */
/*                           STM_GPATCH_CONFIG_SIMP                          */
/*                                                                           */
/*  Return lists of parameters that will configure a gabor patch display.    */
/*                                                                           */
/*  This routine handles "simple" configurations, which are those where      */
/*  the patches are at the same grid points in the LE and RE, and have       */
/*  the opposite direction of movement.                                      */
/*                                                                           */
/*****************************************************************************/
void stm_gpatch_config_simp(gxn,gyn,boxw,boxh,seed,n,dirbal,regflag,
			    rcx,rcy,rdiri,rph1,rph2)
     int gxn,gyn;       // grid dimensions
     float boxw,boxh;   // grid spacing (deg)
     int seed;          // randomization of positions and/or starting phases
     int n;             // Number of patches
     float dirbal;      // Direction balance: fraction assigned to dir1
     int regflag;       // 0-random, 1-checkerboard
     float **rcx;       // [n] x-coord
     float **rcy;       // [n] y-coord
     int   **rdiri;     // [n] direction *index*, e.g., 1 or 2
     float **rph1;      // [n] initial phase
     float **rph2;      // [n] initial phase
{
  int i,j,k;
  int **flag,cnt,nless,*diri;
  float *cx,*cy,*ph1,*ph2,gx0,gy0,x,y,dir_more,dir_less;

  if (regflag == 1)  // Regular grid has patches in all locations.
    n = gxn*gyn;

  cx   = (float *)myalloc(n*sizeof(float));
  cy   = (float *)myalloc(n*sizeof(float));
  diri = (int   *)myalloc(n*sizeof(int));
  ph1  = (float *)myalloc(n*sizeof(float));
  ph2  = (float *)myalloc(n*sizeof(float));

  if (seed > 0)
    seed *= -1;

  //
  //  Set the 'flag' array to signal which grid points will have patches.
  //
  if (regflag == 1){
    //
    //  Regular - patches at all grid locations
    //
    flag = get_const_2d_iarray(gxn,gyn,1);
  }else if (n > (gxn * gyn)/2){
    flag = get_const_2d_iarray(gxn,gyn,1);
    cnt = gxn*gyn;  // Number of flags set
    while(cnt > n){  // While we have too many flags set, unset one randomly
      i = (int)((float)gxn * myrand_util_ran2(&seed));
      j = (int)((float)gyn * myrand_util_ran2(&seed));
      if (flag[i][j] == 1){
	flag[i][j] = 0;
	cnt -= 1;
      }
    }
  }else{
    flag = get_const_2d_iarray(gxn,gyn,0);
    cnt = 0;
    while(cnt < n){  // While we have too few flags set, set one randomly
      i = (int)((float)gxn * myrand_util_ran2(&seed));
      j = (int)((float)gyn * myrand_util_ran2(&seed));
      if (flag[i][j] == 0){
	flag[i][j] = 1;
	cnt += 1;
      }
    }
  }

  //
  //  Set the grid coordinates for the flagged patches, and
  //    pick the random starting phases.
  //
  gx0 = -((float)(gxn - 1) * boxw) / 2.0; // x-coord of left grid point (deg)
  gy0 = -((float)(gyn - 1) * boxh) / 2.0; // y-coord of bottom point (deg)
  x = gx0;
  y = gy0;

  k = 0;
  for(i=0;i<gxn;i++){
    for(j=0;j<gyn;j++){
      if (flag[i][j] == 1){

	cx[k] = x;
	cy[k] = y;
	ph1[k] = 360.0 * myrand_util_ran2(&seed);
	ph2[k] = 360.0 * myrand_util_ran2(&seed);

	k += 1;
      }
      y += boxh;
    }
    x += boxw;
    y = gy0;
  }


  if (regflag == 1){
    //
    //  Checkerboard
    //
    k = 0;
    for(i=0;i<gxn;i++){
      for(j=0;j<gyn;j++){
	if (i%2 == 0){
	  if (j%2 == 0)
	    diri[k] = 1;       // flag indicated which direction
	  else
	    diri[k] = 2;
	}else{
	  if (j%2 == 0)
	    diri[k] = 2;
	  else
	    diri[k] = 1;
	}
	k += 1;
      }
    }
  }else{
    //
    //  Distribute 'dir1' and 'dir2' according to 'dirbal'
    //
    if ((dirbal < 0.0) || (dirbal > 1.0))
      exit_error("STM_GPATCH_CONFIG_SIMP","Bad 'dir_balance', must be [0..1]");

    if (dirbal <= 0.5){
      dir_more = 2;
      dir_less = 1;
      nless = my_rint(dirbal * (float)n);
    }else{
      dir_more = 1;
      dir_less = 2;
      nless = n - my_rint(dirbal * (float)n);
    }

    for(i=0;i<n;i++)
      diri[i] = dir_more;  // Start all with the more numerous value

    cnt = 0;   // Number of less numerous values so far
    while(cnt < nless){  // need more 'dir1' patches
      i = (int)((float)n * myrand_util_ran2(&seed));
      if (diri[i] != dir_less){
	diri[i] = dir_less;
	cnt += 1;
      }
    }
  }

  free_2d_iarray(flag,gxn);

  *rcx   = cx;
  *rcy   = cy;
  *rdiri = diri;
  *rph1  = ph1;
  *rph2  = ph2;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_COHERENCE_DOT_MOVIE                          */
/*                                                                           */
/*  PARAMETERS:                                                              */
/*                                                                           */
/*  sscale       deg/pix                                                     */
/*  framerate    Video rate (frames/s)                                       */
/*  coher        Motion coherence, from 0 to 1  (Britten et al 1992).        */
/*  speed        deg/s                                                       */
/*  theta        Angle in degrees.                                           */
/*  sizex        Size of rectangular frame for making/wrapping dots (pix)    */
/*  sizey        Size of rectangular frame for making/wrapping dots (pix)    */
/*  nframes      Time duration of stimulus (video frames).                   */
/*  dpf          Dots per frame.                                             */
/*  dt           Frames between refresh of signal dots (AKA ref_frames).     */
/*  seed         Address of randomization seed.                              */
/*                                                                           */
/*  Dot coordinates are returned in pixels (float) within a rectangle        */
/*  of 'sizex' by 'sizey' pixels.                                            */
/*                                                                           */
/*****************************************************************************/
void stm_get_coherence_dot_movie(sscale,framerate,coher,speed,theta,sizex,
				 sizey,nframes,dpf,dt,seed,rdotx,rdoty,rdotn)
     float sscale,framerate,coher,speed,theta,sizex,sizey;
     int nframes,dpf,dt,seed;
     float ***rdotx,***rdoty;
     int **rdotn;
{
  int i,j,k;
  int *dotn,n;
  float **dotx,**doty,x,y,dx,dy;

  if (seed > 0)
    seed = -seed;

  dotn = (int *)myalloc(nframes*sizeof(int));
  dotx = (float **)myalloc(nframes*sizeof(float *));
  doty = (float **)myalloc(nframes*sizeof(float *));
  for(i=0;i<nframes;i++){
    dotn[i] = dpf;
    dotx[i] = (float *)myalloc(dpf*sizeof(float));
    doty[i] = (float *)myalloc(dpf*sizeof(float));
  }

  dy = speed * (float)dt/framerate / sscale * sin(theta/180.0*M_PI);
  dx = speed * (float)dt/framerate / sscale * cos(theta/180.0*M_PI);

  for (i=0;i<nframes;i++){
    n = 0; // number of dots created so far for the ith frame
    k = i - dt; // k is index to dots in a previous frame
    if (k >= 0){ // START WITH DOTS FROM EARLIER FRAME, if that frame exists
      for (j=0;j<dpf;j++) // for each dot
	if (myrand_util_ran2(&seed) < coher){ // REPLOT AS SIGNAL DOT
	  x = dotx[k][j] + dx;
	  y = doty[k][j] + dy;
	  // Take care of wrap-around
	  if (x < 0.0) x += sizex; else if (x >= sizex) x -= sizex;
	  if (y < 0.0) y += sizey; else if (y >= sizey) y -= sizey;
	  dotx[i][n] = x; // store the coordinates
	  doty[i][n] = y;
	  n += 1;
	}
    }

    while (n < dpf){ // Choose additional, randomly placed, dots
      dotx[i][n] = myrand_util_ran2(&seed)*sizex;
      doty[i][n] = myrand_util_ran2(&seed)*sizey;
      n += 1;
    }
  }

  *rdotn = dotn; *rdotx = dotx; *rdoty = doty;
}
/**************************************-**************************************/
/*                                                                           */
/*                             STM_MASK_DOT_MOVIE                            */
/*                                                                           */
/*****************************************************************************/
void stm_mask_dot_movie(dotx,doty,dotn,nframes,cx,cy,w,h,sscale,aptype)
     float **dotx,**doty;    // Dot coords (pixels)
     int *dotn,nframes;
     float cx,cy,w,h;        // Aperture center, width, and height (deg)
     float sscale;           // Spatial scale (deg/pix)
     int aptype;
{
  int i,j;
  int n,dflag;
  float x,y,r2pix,tx,ty,apcx,apcy,wo2,ho2;

  // For rectangular aperture
  wo2 = w/(2.0*sscale);   // w/2 (pixels)
  ho2 = h/(2.0*sscale);   // h/2 (pixels)

  // For circular aperture
  r2pix = wo2*wo2;        // pixel radius squared

  // Aperture center in pixel coords
  apcx = cx/sscale;
  apcy = cy/sscale;

  for (i=0;i<nframes;i++){
    n = 0; // number of dots kept so far in the ith frame
    for (j=0;j<dotn[i];j++){
      tx = dotx[i][j];
      ty = doty[i][j];
      x = tx - apcx;      // Coords relative to ap center (pixels)
      y = ty - apcy;
      //y *= -1;    // INVERT y-coordinate

      dflag = 0;
      if (aptype == 0){ // none
	dflag = 1;
      }else if (aptype == 1){ // circular
	if ((x*x + y*y) <= r2pix)
	  dflag = 1;
      }else if (aptype == 2){ // rectangular
	if ((x >= -wo2) && (x <= wo2) && (y >= -ho2) && (y <= ho2))
	  dflag = 1;
      }else{
	printf("*** Aperture type %d not implemented yet\n",aptype);
	dflag = 1;
      }

      if (dflag){
	dotx[i][n] = tx; // Over-write old coordinates
	doty[i][n] = ty;
	n += 1;
      }
    }
    dotn[i] = n;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             STM_DOT_MOVIE_SHIFT                           */
/*                                                                           */
/*  Add the offsets to the dot coordinates.                                  */
/*                                                                           */
/*****************************************************************************/
void stm_dot_movie_shift(dotx,doty,dotn,nframes,dx,dy)
     float **dotx,**doty;   // [nframes][dotn[...]] Dot coords (pixels)
     int *dotn;             // [nframes] number of dots per frame
     int nframes;           // Number of frames
     float dx,dy;           // offset values to add to dot coordinates
{
  int i,j;
  int nd;

  for(i=0;i<nframes;i++){
    nd = dotn[i];  // Number of dots in this frame
    for(j=0;j<nd;j++){
      dotx[i][j] += dx;
      doty[i][j] += dy;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             STM_DOT_MOVIE_COPY                            */
/*                                                                           */
/*  Return a copy of the dot movie, including dot amplitude 'dota' unless    */
/*  'dota' is NULL, in which case do not return the amplitude.               */
/*                                                                           */
/*****************************************************************************/
void stm_dot_movie_copy(dotx,doty,dota,dotn,nframes,rdotx,rdoty,rdota,rdotn)
     float **dotx,**doty;      // [nframes][dotn[...]] Dot coords (pixels)
     float **dota;             // [nframes][dotn[...]] Dot amplitude or NULL
     int *dotn;                // [nframes] number of dots per frame
     int nframes;              // Number of frames
     float ***rdotx,***rdoty;  // Return copy of dot coords
     float ***rdota;           // Return copy of dot amplitudes
     int **rdotn;              // Return copy of number of dots per frame
{
  int i;
  int nd,*cdn;
  float **cdx,**cdy,**cda;

  cdx = (float **)myalloc(nframes*sizeof(float *));
  cdy = (float **)myalloc(nframes*sizeof(float *));
  if (dota != NULL)
    cda = (float **)myalloc(nframes*sizeof(float *));
  cdn = (int *)myalloc(nframes*sizeof(int));

  for(i=0;i<nframes;i++){
    nd = dotn[i];  // Number of dots in this frame
    cdx[i] = copy_farray(dotx[i],nd);
    cdy[i] = copy_farray(doty[i],nd);
    if (dota != NULL)
      cda[i] = copy_farray(dota[i],nd);
    cdn[i] = nd;
  }

  *rdotx = cdx;
  *rdoty = cdy;
  if (dota != NULL)
    *rdota = cda;
  *rdotn = cdn;
}
/**************************************-**************************************/
/*                                                                           */
/*                              STM_DOTMOV_INIT                              */
/*                                                                           */
/*  Create, initialize and return a 'dotmov' structure.                      */
/*  If a 'flag' array is created here, all flag values will be set to 0.     */
/*                                                                           */
/*****************************************************************************/
struct stmh_dotmov *stm_dotmov_init(nframes,dpf,aflag,zflag,fflag)
     int nframes;  // Number of frames in the dot movie
     int dpf;      // Number of dots per frame (DPF)
     int aflag;    // 1-create amplitude storage; 0-don't
     int zflag;    // 1-create a z-dimension; 0-don't
     int fflag;    // 1-create a flag array; 0-don't
{
  int i;
  struct stmh_dotmov *dm;

  //printf("nframes = %d\n",nframes);
  //printf("dpf = %d\n",dpf);
  //printf("a,z,f flags = %d %d %d\n",aflag,zflag,fflag);

  dm = (struct stmh_dotmov *)myalloc(sizeof(struct stmh_dotmov));

  dm->nframes = nframes;
  dm->cnt =    (int *)myalloc(nframes*sizeof(int));
  dm->x   = (float **)myalloc(nframes*sizeof(float *));
  dm->y   = (float **)myalloc(nframes*sizeof(float *));

  if (aflag == 1)
    dm->a   = (float **)myalloc(nframes*sizeof(float *));
  else
    dm->a = NULL;

  if (zflag == 1)
    dm->z   = (float **)myalloc(nframes*sizeof(float *));
  else
    dm->z = NULL;

  if (fflag == 1)
    dm->flag = (int **)myalloc(nframes*sizeof(int *));
  else
    dm->flag = NULL;


  for(i=0;i<nframes;i++){
    dm->cnt[i] = dpf;
    dm->x[i] = (float *)myalloc(dpf*sizeof(float));
    dm->y[i] = (float *)myalloc(dpf*sizeof(float));
    if (aflag == 1)
      dm->a[i] = (float *)myalloc(dpf*sizeof(float));
    if (zflag == 1)
      dm->z[i] = (float *)myalloc(dpf*sizeof(float));
    if (fflag == 1)
      dm->flag[i] = get_zero_iarray(dpf);
  }

  return dm;
}
/**************************************-**************************************/
/*                                                                           */
/*                            STM_DOTMOV_GET_PARTS                           */
/*                                                                           */
/*  Return separate coordinate arrays from the dotmov struct.                */
/*                                                                           */
/*****************************************************************************/
void stm_dotmov_get_parts(dm,fflag,aflag,rx,ry,ra,rcnt)
     struct stmh_dotmov *dm;
     int fflag;    // 1-return only dots with flag = 1;  0-return all dots
     int aflag;    // 1-return 'ra' amplitude values; 0-don't
     float ***rx;  // [dm->nframes][rcnt[]] x-coords 
     float ***ry;  // [dm->nframes][rcnt[]] y-coords 
     float ***ra;  // [dm->nframes][rcnt[]] amplitude; or NULL if aflag = 0
     int  **rcnt;  // [dm->nframes] number of dots per frame
{
  int i,j,k;
  int n,*dotn,ntot,nval;
  float **dotx,**doty,**dota;

  n = dm->nframes;
  dotn = get_zero_iarray(n);
  dotx = (float **)myalloc(n*sizeof(float *));
  doty = (float **)myalloc(n*sizeof(float *));
  if (aflag == 1)
    dota = (float **)myalloc(n*sizeof(float *));
  else
    dota = NULL;

  if (fflag == 1){
    //
    //  Return data for only those dots with 'flag' = 1
    //

    for(i=0;i<n;i++){  // For each frame

      ntot = dm->cnt[i];  // Total dots
      nval = sum_iarray(dm->flag[i],ntot,0,ntot);  // Number valid for copying

      dotx[i] = (float *)myalloc(nval*sizeof(float));
      doty[i] = (float *)myalloc(nval*sizeof(float));
      if (aflag == 1)
	dota[i] = (float *)myalloc(nval*sizeof(float));

      k = 0;
      for(j=0;j<ntot;j++){
	if (dm->flag[i][j] == 1){  // For each valid dot

	  dotx[i][k] = dm->x[i][j];   // Save it, to be returned
	  doty[i][k] = dm->y[i][j];
	  if (aflag == 1)
	    dota[i][k] = dm->a[i][j];

	  k += 1;
	}
      }
      if (k != nval)
	exit_error("STM_DOTMOV_GET_PARTS","WYETH - this should not happen");

      dotn[i] = nval;
    }
  }else{
    //
    //  Return all data
    //

    exit_error("STM_DOTMOV_GET_PARTS","WYETH - not implemented yet");
  }

  *rx   = dotx;
  *ry   = doty;
  *ra   = dota;
  *rcnt = dotn;
}
/**************************************-**************************************/
/*                                                                           */
/*                          STM_DOTMOV_GET_RAND_XYA                          */
/*                                                                           */
/*  Return a movie with random x,y coordinates and random a, z = 0.          */
/*                                                                           */
/*****************************************************************************/
struct stmh_dotmov *stm_dotmov_get_rand_xya(szx,szy,nframes,ndots,seed,z,bipol)
     float szx,szy;  // Range of x,y coordinates
     int nframes;    // Number of frames
     int ndots;      // Number of dots to create in one frame
     int seed;       // Controls random x,y positions of dots
     float z;        // Depth value
     int bipol;      // 1-create bipolar dots; 0-create bright dots
{
  int i,j;
  struct stmh_dotmov *dm;

  if (seed > 0)
    seed = -seed;

  dm = stm_dotmov_init(nframes,ndots,1,1,1);  // Create dotmov struct, w/ flags

  for(i=0;i<nframes;i++){
    for(j=0;j<ndots;j++){
      dm->x[i][j] = szx * myrand_util_ran2(&seed);
      dm->y[i][j] = szy * myrand_util_ran2(&seed);
      dm->z[i][j] = z;
      dm->a[i][j] = 1.0;   // Set all amplitude values to 1.0
      if (bipol == 1){
	if (myrand_util_ran2(&seed) < 0.5)
	  dm->a[i][j] = -1.0;  // dark
      }

      dm->flag[i][j] = 1;  // Set all flag values to 1
    }
  }

  return dm;
}
/**************************************-**************************************/
/*                                                                           */
/*                           STM_DOTMOV_GET_DEEP_F0                          */
/*                                                                           */
/*  Return a single frame of dots evenly spread in depth.                    */
/*                                                                           */
/*****************************************************************************/
struct stmh_dotmov *stm_dotmov_get_deep_f0(szx,szy,ndots,pseed,z0,z1,bipol)
     float szx,szy;  // Range of x,y coordinates
     int ndots;      // Number of dots to create in one frame
     int *pseed;     // Pointer to seed controlling random x,y coords. of dots
     float z0,z1;    // Start and end of depth range
     int bipol;      // 1-create bipolar dots; 0-create bright dots
{
  int i;
  float dz,z;
  struct stmh_dotmov *dm;

  dm = stm_dotmov_init(1,ndots,1,1,0);  // Create dotmov struct, no flags

  dz = (z1 - z0) / (float)ndots;

  z = z0 + dz;  // This should give 'z' as:  [-1.975, ..., 2.000]
  for(i=0;i<ndots;i++){
    dm->x[0][i] = szx * myrand_util_ran2(pseed);
    dm->y[0][i] = szy * myrand_util_ran2(pseed);
    dm->z[0][i] = z;
    dm->a[0][i] = 1.0;  // Set all amplitude values to 1.0
    if (bipol == 1){
      if (myrand_util_ran2(pseed) < 0.5)
	dm->a[0][i] = -1.0;  // dark
    }

    z += dz;
  }

  return dm;
}
/**************************************-**************************************/
/*                                                                           */
/*                             STM_DOTMOV_GET_DEEP                           */
/*                                                                           */
/*  Return a dot movie where dots move in depth according to the temporal    */
/*  sequency 'tseq'.                                                         */
/*                                                                           */
/*  Sign conventions:                                                        */
/*  (1) positive 'tseq' value increases dot 'z' toward observer              */
/*  (2) positive 'tseq' value moves dots to right (positively) in left eye   */
/*                                                                           */
/*****************************************************************************/
void stm_dotmov_get_deep(sscale,sizex,sizey,framerate,speed,maxshift,denv,
			 nframes,dpf,dt,seed,tseq,tsn,bipol,dotcorr,rdml,rdmr)
     float sscale;        // (deg/pix)
     float sizex,sizey;   // (pix) extent of x,y dot coordinates
     float framerate;     // (frames/s)
     float speed;         // (deg/s) dot speed in each eye
     float maxshift;      // (deg) maximum disparity before dot is hidden
     int denv;            // Depth envelope: 0-sqr 1-tent 2-sqrt_tent 3-cos
                          //   4-circ
     int nframes;         // number of movie frames
     int dpf;             // visible dots per movie frame
     int dt;              // dwell: number of frames to repeat motion step
     int seed;            // controls x,y dot placement
     float *tseq;         // [tsn] time sequence of movements in 'speed' units
     int tsn;             // length of 'tseq'
     int bipol;           // 1-create amplitude data; 0-don't
     int dotcorr;         // 1-normal, 0-spatially uncorr, -1-opp.sign, Note,
                          //   1 vs. 0 matters here, but -1 not processed here
     struct stmh_dotmov **rdml,**rdmr;  // Return left and right eye movies
{
  int i,j,k;
  int dpf2;
  float z,z0,z1,ppf,xpz,dx,*xref,*yref,xoff,aa,fz;
  struct stmh_dotmov *dm0,*dml,*dmr,*dm0r;

  //
  //  Check the relationship between the motion sequence, dwell and nframes
  //
  if ((tsn * dt) < nframes)
    exit_error("STM_DOTMOV_GET_DEEP","Motion sequence too short");

  dpf2 = 2*dpf;  // Generate twice as many dots, because 1/2 will be hidden

  //
  //  Generate a single frame of dots, evenly in depth, but random in x,y.
  //  The depth values will range from -2 to 2, with the intent that
  //  in the final movie, only dots in the range from -1 to 1 will be shown.
  //
  if (seed > 0)
    seed = -seed;

  z0 = -2.0;
  z1 =  2.0;
  dm0 = stm_dotmov_get_deep_f0(sizex,sizey,dpf2,&seed,z0,z1,bipol);
  if (dotcorr == 0){
    //
    //  Make a spatially independent reference for the RE.
    //  The same seed chain is continued here.
    //  This will have exactly the same 'z' values as 'dm0'.
    //  It will also use the amplitude values of 'dm0'.
    //
    dm0r = stm_dotmov_get_deep_f0(sizex,sizey,dpf2,&seed,z0,z1,bipol);
  }

  //
  //  Compute conversion factor to scale 'tseq' motion to dot coordinates
  //
  ppf = speed / (framerate * sscale) ;  // pix/frame

  //
  //  Compute conversion factor relating z (depth) to x
  //  Divide by 2 here to make 'maxshift' be the max disparity between
  //    the eyes, rather than the maximum shift from the reference.
  //
  xpz = (maxshift/2.0) / sscale;  // x per z (pix per depth)

  //
  //  Generate two movies, one each for the L and R eyes, by applying the
  //  motion sequence 'tseq' to the original dot frame 'dm0'.
  //
  dml = stm_dotmov_init(nframes,dpf2,1,1,1);  // dot movie for left eye
  dmr = stm_dotmov_init(nframes,dpf2,1,1,1);  // dot movie for right eye


  xref = dm0->x[0];  // Reference x-coords
  yref = dm0->y[0];  // Reference y-coords

  k = 0;
  for(i=0;i<nframes;i++){
    dx = ppf*tseq[k];  // Amount to move dots from previous frame

    for(j=0;j<dpf2;j++){  // For each dot
      if (i==0)
	z = dm0->z[0][j];    // Get original z-coord (depth)
      else
	z = dml->z[i-1][j];  // Get previous frame z-coord of dot (depth)
                             //   which is same for left and right eye

      // update z, wrapping around
      z += dx / xpz;
      while (z > 2.0)
	z -= 4.0;
      while (z <= -2.0)
	z += 4.0;

      xoff = z * xpz;  // The x-offset is a function of depth, 'z'

      dml->x[i][j] = xref[j] + xoff;
      dml->y[i][j] = yref[j];
      dml->z[i][j] = z;

      if (dotcorr == 0){
	// Use x,y coords from a separate reference movie
	dmr->x[i][j] = dm0r->x[0][j] - xoff;
	dmr->y[i][j] = dm0r->y[0][j];
      }else{
	//  dotcorr may be 1 or -1 (same or opposite sign), both use same x,y
	dmr->x[i][j] = xref[j] - xoff;
	dmr->y[i][j] = yref[j];
      }
      dmr->z[i][j] = z;

      aa = dm0->a[0][j];
      fz = fabs(z);

      if (fz > 1.0)
	aa = 0.0;
      else{
	if (denv == 1)
	  aa *= 1.0 - fz;         // Tent function: max at z=0
	else if (denv == 2)
	  aa *= sqrt(1.0 - fz);   // Sqrt(tent)
	else if (denv == 3)
	  aa *= cos(z*M_PI*0.5);  // Cosine
	else if (denv == 4)
	  aa *= sqrt(1.0 - z*z);  // Semi Circle
      }
      dml->a[i][j] = aa;
      dmr->a[i][j] = aa;

      if ((z > -1.0) && (z <= 1.0)){
	dml->flag[i][j] = 1;  // Dot is in visible range
	dmr->flag[i][j] = 1;
      }else{
	dml->flag[i][j] = 0;  // Dot is out of visible range
	dmr->flag[i][j] = 0;
      }

    }
    if (((i+1) % dt) == 0)  // E.g., k should be 0,0,1,1,2,2, ... etc if dt=2
      k += 1;
  }

// ********** WYETH FREE dm0, and anything else ????????
// ********** WYETH FREE dm0, and anything else ????????
// ********** WYETH FREE dm0, and anything else ????????
// ********** WYETH FREE dm0, and anything else ????????

  *rdml = dml;
  *rdmr = dmr;
}
/**************************************-**************************************/
/*                                                                           */
/*                           STM_DOTMOV_GET_PLANES                           */
/*                                                                           */
/*  Return a dot movie where several planes defined by incoherent RDSs       */
/*  move in depth according to the sequence 'tseq'.                          */
/*                                                                           */
/*****************************************************************************/
void stm_dotmov_get_planes(sscale,sizex,sizey,framerate,speed,maxshift,denv,
			   nframes,dpf,nplanes,dt,seed,tseq,tsn,bipol,rdml,
			   rdmr)
     float sscale;        // (deg/pix)
     float sizex,sizey;   // (pix) extent of x,y dot coordinates
     float framerate;     // (frames/s)
     float speed;         // (deg/s) dot speed in each eye
     float maxshift;      // (deg) maximum shift before dot is hidden
     int denv;            // Depth envelope: 0-sqr 1-tent 2-sqrt_tent 3-cos
                          //   4-circ
     int nframes;         // number of movie frames
     int dpf;             // visible dots per movie frame
     int nplanes;         // number of planes of dots
     int dt;              // dwell: number of frames to repeat motion step
     int seed;            // controls x,y dot placement
     float *tseq;         // [tsn] time sequence of movements in 'speed' units
     int tsn;             // length of 'tseq'
     int bipol;           // 1-bipolar, 0-unipolar
     struct stmh_dotmov **rdml,**rdmr;  // Return left and right eye movies
{
  int i,j,k;
  float z,ppf,xpz,dx,*xref,*yref,xoff,aa,fz;
  struct stmh_dotmov *dml,*dmr;

  //
  //  Check the relationship between the motion sequence, dwell and nframes
  //
  if ((tsn * dt) < nframes)
    exit_error("STM_DOTMOV_GET_PLANES","Motion sequence too short");

  //
  //  Generate a single frame of dots, evenly in depth, but random in x,y.
  //  The depth values will range from -2 to 2, with the intent that
  //  in the final movie, only dots in the range from -1 to 1 will be shown.
  //
  z = 0.0;
  dml = stm_dotmov_get_rand_xya(sizex,sizey,nframes,dpf,seed,z,bipol);

  if (nplanes == 2){
    //
    //  WYETH - could extend to arbitrary number of planes?
    //
    for(i=0;i<dpf/2;i++)
      dml->z[0][i] = 1.0;  // Move depth of half of dots to a second plane
  }

  //
  //  Compute conversion factor to scale 'tseq' motion to dot coordinates
  //
  ppf = speed / (framerate * sscale) ;  // pix/frame

  //
  //  Compute conversion factor relating z (depth) to x
  //
  xpz = maxshift / sscale;  // x per z (pix per depth)

  //
  //  Generate storage for movie for right eye.
  //
  dmr = stm_dotmov_init(nframes,dpf,1,1,1);


  k = 0;
  for(i=0;i<nframes;i++){

    xref = dml->x[i];  // right eye x-coords
    yref = dml->y[i];  // right eye y-coords

    dx = ppf*tseq[k];  // Amount to move dots from previous frame

    for(j=0;j<dpf;j++){  // For each dot
      if (i == 0)
	z = dml->z[i][j];    // z-coord (depth) from left eye
      else
	z = dml->z[i-1][j];  // Get previous frame z-coord of dot (depth)


      // update z, wrapping around
      z += dx / xpz;
      while (z > 1.0)
	z -= 2.0;
      while (z <= -1.0)
	z += 2.0;

      xoff = z * xpz;  // The x-offset is a function of depth, 'z'

      dmr->x[i][j] = xref[j] - xoff;
      dmr->y[i][j] = yref[j];
      dmr->z[i][j] = z;
      dml->z[i][j] = z;  // WYETH - strange - we're changing the initial z??

      aa = dml->a[i][j];
      fz = fabs(z);

      if (fz > 1.0)
	aa = 0.0;    // This should not happen, z is from -1 to 1
      else{
	if (denv == 1)
	  aa *= 1.0 - fz;         // Tent function: max at z=0
	else if (denv == 2)
	  aa *= sqrt(1.0 - fz);   // Sqrt(tent)
	else if (denv == 3)
	  aa *= cos(z*M_PI*0.5);  // Cosine
	else if (denv == 4)
	  aa *= sqrt(1.0 - z*z);  // Semi Circle
      }
      dml->a[i][j] = aa;
      dmr->a[i][j] = aa;

      dml->flag[i][j] = 1;  // Dot is in visible range
      dmr->flag[i][j] = 1;

    }
    if (((i+1) % dt) == 0)  // E.g., k should be 0,0,1,1,2,2, ... etc if dt=2
      k += 1;
  }

  *rdml = dml;
  *rdmr = dmr;
}
/**************************************-**************************************/
/*                                                                           */
/*                              STM_STAR_GET_TRIG                            */
/*                                                                           */
/*  Trigger times are returned in 'sampling' units (typically msec) and are  */
/*  specified relative to time 0 being the start of the first adapt period.  */
/*  The 'tflag' is used to select triggers for specific conditions:          */
/*                                                                           */
/*    ad   - the beginning of each adapt period                              */
/*    ad0  - the beginning of adapt for state0 = 0                           */
/*    ad1  - the beginning of adapt for state0 = 1                           */
/*                                                                           */
/*****************************************************************************/
void stm_star_get_trig(stn,fps,sampling,dt_ad,dt_t0,dt_t1,state0,ttype,
		       rtrig,rntrig)
     float stn;       // stimulus duration (s)
     float fps;       // frames per second
     float sampling;  // samples/sec for returned triggers, 1000.0 for msec
     int dt_ad;       // adapt duration (frames)
     int dt_t1;       // test0 duration (frames)
     int dt_t0;       // test1 duration (frames)
     int state0;      // starting state
     char *ttype;     // type of trigger to return, see comments above
     int **rtrig;     // trigger times
     int *rntrig;     // number of trigggers
{
  int i,k;
  int frpcyc,ncyc,*t,first_flag;
  double dtsec,t0s;

  frpcyc = dt_ad + dt_t0 + dt_t1;  // Frames/cycle
  ncyc = (int)(stn*fps) / frpcyc;  // Number of cycles

  dtsec = (double)frpcyc / (double)fps;  // seconds for one full cycle

  t = (int *)myalloc(ncyc*sizeof(int));  // Maximum, all may not be used

  k = 0;
  first_flag = state0;
  for(i=0;i<ncyc;i++){

    t0s = (double)i*dtsec;  // Start time of adapt on ith cycle (sec)
    
    if (strcmp(ttype,"ad")==0){        // All adapt
      t[k] = my_rint(sampling * t0s);
      k += 1;
    }else if (strcmp(ttype,"ad0")==0){
      if (first_flag == 0){
	t[k] = my_rint(sampling * t0s);
	k += 1;
      }
    }else if (strcmp(ttype,"ad1")==0){
      if (first_flag == 1){
	t[k] = my_rint(sampling * t0s);
	k += 1;
      }
    }else{
      printf("  trig_type  %s\n",ttype);
      exit_error("STM_STAR_GET_TRIG","Unknown trig_type");
    }

    first_flag = 1 - first_flag;  // toggles between 0,1
  }
  
  *rntrig = k;   // Number of triggers
  *rtrig = t;    // trigger values (msec)
}
/**************************************-**************************************/
/*                                                                           */
/*                       STM_PASUPATHY_SHAPE_IDROT_LIST                      */
/*                                                                           */
/*  Return two lists, 'rid' contains the shape IDs, e.g., [1..51], and       */
/*  'rrot' contains the degrees of rotation [0..360).                        */
/*                                                                           */
/*****************************************************************************/
void stm_pasupathy_shape_idrot_list(mylogf,stype,rn,rid,rrot)
     char *mylogf;  // Log file, or NULL
     char *stype;   // Type of shape list: "pc370"
     int  *rn;      // Length of returned lists
     int **rid;     // [*rn] Shape ID
     int **rrot;    // [*rn] Shape rotation, deg
{
  int i,j,k;
  int n,nr,ddeg,*shid,*rot;

  if (strcmp(stype,"pc370")==0){
    //
    //  pc370:  Pasupathy and Connor (2001), 370 shapes
    //
    n = 370;               // Number of unique stimuli
    shid = (int *)myalloc(n*sizeof(int));  // [n] Shape ID
    rot  = (int *)myalloc(n*sizeof(int));  // [n] Sape rotation (deg)

    k = 0;
    for(i=0;i<51;i++){

      nr = pasu_shape_nrot[i];  // Number of rotations per shape
      //
      //  ************* WYETH BUG ***************  Apr 1
      //  For 1 rotations, should be:  0
      //  For 2 rotations, should be:  0, 45
      //  For 4 rotations, should be:  0, 45, 90, 135
      //  For 8 rotations, should be:  0, 45, 90, 135, 180, 225, 270, 315
      //
      //ddeg = 360 / nr;    // ERROR
      ddeg = 45;            // Corrected on 2017 Oct 27
      for(j=0;j<nr;j++){
	shid[k] = i+1;
	rot[k] = j*ddeg;
	k += 1;
      }
    }
    if (k != 370)
      exit_error("STM_PASUPATHY_SHAPE_IDROT_LIST","Should not happen");
  }else{
    printf("  stype:  %s\n",stype);
    exit_error("STM_PASUPATHY_SHAPE_IDROT_LIST","Bad 'stype'");
  }

  *rn = n;
  *rid = shid;
  *rrot = rot;
}
/**************************************-**************************************/
/*                                                                           */
/*                        STM_PASUPATHY_SHAPE_GET_CONTROL                    */
/*                                                                           */
/*  Return the values and number of control points for each shape.           */
/*  NOTE:  Storage is created - called should free 'rx' and 'ry'.            */
/*                                                                           */
/*****************************************************************************/
void stm_pasupathy_shape_get_control(mylogf,shid,rx,ry,rn)
     char *mylogf;  // Log file, or NULL
     int shid;      // Shape ID
     float **rx;    // Return x-coords [*rn]
     float **ry;    // Return y-coords [*rn]
     int    *rn;    // Return number of points
{
  int i;
  int n;
  float *t,*x,*y;

  if ((shid < 1) || (shid > 51)){
    mylogx(mylogf,"STM_PASUPATHY_SHAPE_GET_CONTROL",
	   "Shape ID must be [1..51]");
  }

  if      (shid ==  1)  t = pasu_shape_01;
  else if (shid ==  2)  t = pasu_shape_02;
  else if (shid ==  3)  t = pasu_shape_03;
  else if (shid ==  4)  t = pasu_shape_04;
  else if (shid ==  5)  t = pasu_shape_05;
  else if (shid ==  6)  t = pasu_shape_06;
  else if (shid ==  7)  t = pasu_shape_07;
  else if (shid ==  8)  t = pasu_shape_08;
  else if (shid ==  9)  t = pasu_shape_09;
  else if (shid == 10)  t = pasu_shape_10;
  else if (shid == 11)  t = pasu_shape_11;
  else if (shid == 12)  t = pasu_shape_12;
  else if (shid == 13)  t = pasu_shape_13;
  else if (shid == 14)  t = pasu_shape_14;
  else if (shid == 15)  t = pasu_shape_15;
  else if (shid == 16)  t = pasu_shape_16;
  else if (shid == 17)  t = pasu_shape_17;
  else if (shid == 18)  t = pasu_shape_18;
  else if (shid == 19)  t = pasu_shape_19;
  else if (shid == 20)  t = pasu_shape_20;
  else if (shid == 21)  t = pasu_shape_21;
  else if (shid == 22)  t = pasu_shape_22;
  else if (shid == 23)  t = pasu_shape_23;
  else if (shid == 24)  t = pasu_shape_24;
  else if (shid == 25)  t = pasu_shape_25;
  else if (shid == 26)  t = pasu_shape_26;
  else if (shid == 27)  t = pasu_shape_27;
  else if (shid == 28)  t = pasu_shape_28;
  else if (shid == 29)  t = pasu_shape_29;
  else if (shid == 30)  t = pasu_shape_30;
  else if (shid == 31)  t = pasu_shape_31;
  else if (shid == 32)  t = pasu_shape_32;
  else if (shid == 33)  t = pasu_shape_33;
  else if (shid == 34)  t = pasu_shape_34;
  else if (shid == 35)  t = pasu_shape_35;
  else if (shid == 36)  t = pasu_shape_36;
  else if (shid == 37)  t = pasu_shape_37;
  else if (shid == 38)  t = pasu_shape_38;
  else if (shid == 39)  t = pasu_shape_39;
  else if (shid == 40)  t = pasu_shape_40;
  else if (shid == 41)  t = pasu_shape_41;
  else if (shid == 42)  t = pasu_shape_42;
  else if (shid == 43)  t = pasu_shape_43;
  else if (shid == 44)  t = pasu_shape_44;
  else if (shid == 45)  t = pasu_shape_45;
  else if (shid == 46)  t = pasu_shape_46;
  else if (shid == 47)  t = pasu_shape_47;
  else if (shid == 48)  t = pasu_shape_48;
  else if (shid == 49)  t = pasu_shape_49;
  else if (shid == 50)  t = pasu_shape_50;
  else if (shid == 51)  t = pasu_shape_51;

  n = pasu_shape_n[shid-1] / 2;  // Number of coords

  x = get_farray(n);
  y = get_farray(n);

  for(i=0;i<n;i++){
    x[i] = t[i*2];
    y[i] = t[i*2+1];
  }

  *rn = n;
  *rx = x;
  *ry = y;
}
/**************************************-**************************************/
/*                                                                           */
/*                       STM_PASUPATHY_YOSH_GET_CONTROL                      */
/*                                                                           */
/*  Return the values and number of control points for each shape.           */
/*  NOTE:  Storage is created - called should free 'rx' and 'ry'.            */
/*                                                                           */
/*****************************************************************************/
void stm_pasupathy_yosh_get_control(mylogf,shid,rx,ry,rn)
     char *mylogf;  // Log file, or NULL
     int shid;      // Shape ID
     float **rx;    // Return x-coords [*rn]
     float **ry;    // Return y-coords [*rn]
     int    *rn;    // Return number of points
{
  int i;
  int n;
  float *t,*x,*y;

  if ((shid < 1) || (shid > 43)){
    mylogx(mylogf,"STM_PASUPATHY_YOSH_GET_CONTROL",
	   "Shape ID must be [1..43]");
  }

  if      (shid ==  1)  t =                pasu_shape_32; //yosh_shape_00; 
  else if (shid ==  2)  t =                pasu_shape_33; //yosh_shape_01;
  else if (shid ==  3)  t =                pasu_shape_34; //yosh_shape_02;
  else if (shid ==  4)  t =                pasu_shape_35; //yosh_shape_03;
  else if (shid ==  5)  t =                pasu_shape_36; //yosh_shape_04;
  else if (shid ==  6)  t = yosh_shape_05;
  else if (shid ==  7)  t = yosh_shape_06;
  else if (shid ==  8)  t =                pasu_shape_24; //yosh_shape_07;
  else if (shid ==  9)  t =                pasu_shape_25; //yosh_shape_08;
  else if (shid == 10)  t =                pasu_shape_26; //yosh_shape_09;
  else if (shid == 11)  t =                pasu_shape_27; //yosh_shape_10;
  else if (shid == 12)  t =                pasu_shape_28; //yosh_shape_11;
  else if (shid == 13)  t =                pasu_shape_29; //yosh_shape_12;
  else if (shid == 14)  t =                pasu_shape_30; //yosh_shape_13;
  else if (shid == 15)  t =                pasu_shape_31; //yosh_shape_14;
  else if (shid == 16)  t = yosh_shape_15; // Similar to 40, but not exact
  else if (shid == 17)  t = yosh_shape_16;
  else if (shid == 18)  t = yosh_shape_17;
  else if (shid == 19)  t = yosh_shape_18;
  else if (shid == 20)  t = yosh_shape_19;
  else if (shid == 21)  t = yosh_shape_20; // Similar to 43, but not exact
  else if (shid == 22)  t =                pasu_shape_19; //yosh_shape_21;
  else if (shid == 23)  t =                pasu_shape_16; //yosh_shape_22;
  else if (shid == 24)  t =                pasu_shape_17; //yosh_shape_23;
  else if (shid == 25)  t =                pasu_shape_18; //yosh_shape_24;
  else if (shid == 26)  t =                pasu_shape_21; //yosh_shape_25;
  else if (shid == 27)  t =                pasu_shape_22; //yosh_shape_26;
  else if (shid == 28)  t =                pasu_shape_23; //yosh_shape_27;
  else if (shid == 29)  t =                pasu_shape_50; //yosh_shape_28;
  else if (shid == 30)  t =                pasu_shape_38; //yosh_shape_29;
  else if (shid == 31)  t = yosh_shape_30;
  else if (shid == 32)  t =                pasu_shape_41; //yosh_shape_31;
  else if (shid == 33)  t = yosh_shape_32;
  else if (shid == 34)  t = yosh_shape_33;
  else if (shid == 35)  t = yosh_shape_34;
  else if (shid == 36)  t =                pasu_shape_05; //yosh_shape_35;
  else if (shid == 37)  t =                pasu_shape_06; //yosh_shape_36;
  else if (shid == 38)  t =                pasu_shape_07; //yosh_shape_37;
  else if (shid == 39)  t =                pasu_shape_10; //yosh_shape_38;
  else if (shid == 40)  t =                pasu_shape_11; //yosh_shape_39;
  else if (shid == 41)  t =                pasu_shape_14; //yosh_shape_40;
  else if (shid == 42)  t =                pasu_shape_15; //yosh_shape_41;
  else if (shid == 43)  t = yosh_shape_42;

  n = yosh_shape_n[shid-1] / 2;  // Number of coords

  x = get_farray(n);
  y = get_farray(n);

  for(i=0;i<n;i++){
    x[i] = t[i*2];
    y[i] = t[i*2+1];
  }

  *rn = n;
  *rx = x;
  *ry = y;
}
/**************************************-**************************************/
/*                                                                           */
/*                       STM_PASUPATHY_YASM_GET_CONTROL                      */
/*                                                                           */
/*  Return the values and number of control points for each shape.           */
/*  NOTE:  Storage is created - called should free 'rx' and 'ry'.            */
/*                                                                           */
/*****************************************************************************/
void stm_pasupathy_yasm_get_control(mylogf,shid,rx,ry,rn)
     char *mylogf;  // Log file, or NULL
     int shid;      // Shape ID
     float **rx;    // Return x-coords [*rn]
     float **ry;    // Return y-coords [*rn]
     int    *rn;    // Return number of points
{
  int i;
  int n;
  float *t,*x,*y;

  if ((shid < 1) || (shid > 13)){
    mylogx(mylogf,"STM_PASUPATHY_YASM_GET_CONTROL",
	   "Shape ID must be [1..13]");
  }

  if      (shid ==  1)  t = yasm_shape_01;
  else if (shid ==  2)  t = yasm_shape_02;
  else if (shid ==  3)  t = yasm_shape_03;
  else if (shid ==  4)  t = yasm_shape_04;
  else if (shid ==  5)  t = yasm_shape_05;
  else if (shid ==  6)  t = yasm_shape_06;
  else if (shid ==  7)  t = yasm_shape_07;
  else if (shid ==  8)  t = yasm_shape_08;
  else if (shid ==  9)  t = yasm_shape_09;
  else if (shid == 10)  t = yasm_shape_10;
  else if (shid == 11)  t = yasm_shape_11;
  else if (shid == 12)  t = yasm_shape_12;
  else if (shid == 13)  t = yasm_shape_13;

  n = yasm_shape_n[shid-1] / 2;  // Number of coords

  x = get_farray(n);
  y = get_farray(n);

  for(i=0;i<n;i++){
    x[i] = t[i*2];
    y[i] = t[i*2+1];
  }

  *rn = n;
  *rx = x;
  *ry = y;
}
/**************************************-**************************************/
/*                                                                           */
/*                         STM_PASUPATHY_SHAPE_COMPARE                       */
/*                                                                           */
/*  Compare the control points for two shapes.                               */
/*                                                                           */
/*****************************************************************************/
void stm_pasupathy_shape_compare(sset1,shid1,sset2,shid2)
     char *sset1;   // "pc2001", "yasm"
     int shid1;     // Shape ID
     char *sset2;   // "pc2001", "yasm"
     int shid2;     // Shape ID
{
  int i;
  int sn,sn2,flag;
  float *sx,*sy,*sx2,*sy2;

  if (strcmp(sset1,"yasm")==0)
    stm_pasupathy_yasm_get_control(NULL,shid1,&sx,&sy,&sn);
  else if (strcmp(sset1,"yosh")==0)
    stm_pasupathy_yosh_get_control(NULL,shid1,&sx,&sy,&sn);
  else
    stm_pasupathy_shape_get_control(NULL,shid1,&sx,&sy,&sn);

  if (strcmp(sset2,"yasm")==0)
    stm_pasupathy_yasm_get_control(NULL,shid2,&sx2,&sy2,&sn2);
  else if (strcmp(sset2,"yosh")==0)
    stm_pasupathy_yosh_get_control(NULL,shid2,&sx2,&sy2,&sn2);
  else
    stm_pasupathy_shape_get_control(NULL,shid2,&sx2,&sy2,&sn2);

  flag = 0;

  if (sn != sn2){
    printf("    ==>  Point counts differ:  %d %d\n",sn,sn2);
    flag = 1;
  }else{
    for(i=0;i<sn;i++){
      if (sx[i] != sx2[i]){
	printf("    ==>  x-coord %d differs:  %f %f\n",i,sx[i],sx2[i]);
	flag = 1;
      }
      if (sy[i] != sy2[i]){
	printf("    ==>  y-coord %d differs:  %f %f\n",i,sy[i],sy2[i]);
	flag = 1;
      }
    }
  }

  if (flag == 1)
    printf("    _____________ %s %d   %s %d ___________\n",sset1,shid1,
	   sset2,shid2);
  else
    printf("    ==>  SHAPES ARE THE SAME\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                     STM_PASUPATHY_SHAPE_GET_INSIDE_POINT                  */
/*                                                                           */
/*  Return the value of a point that lies inside of the shape.               */
/*  The coordinates are given in the control point coordinate system.        */
/*                                                                           */
/*****************************************************************************/
void stm_pasupathy_shape_get_inside_point(shape_set,shid,theta,rx0,ry0)
     char *shape_set; // "pc2001", "yasm"
     int shid;        // Shape ID
     float theta;     // rotation
     float *rx0;      // Return x-coords [*rn]
     float *ry0;      // Return y-coords [*rn]
{
  float xx[1],yy[1];

  if (strcmp(shape_set,"pc2001")==0){
    if ((shid >= 8) && (shid <= 10)){
      *rx0 = 0.4;
      *ry0 = 0.4;
    }else if (shid == 11){  // OK
      *rx0 = -0.4;
      *ry0 =  0.4;
    }else if (shid == 12){  // OK
      *rx0 =  0.0;
      *ry0 =  0.2;
    }else if (shid == 14){  // OK
      *rx0 =  0.0;
      *ry0 =  0.2;
    }else if (shid == 15){  // OK
      *rx0 =  0.0;
      *ry0 =  0.2;
    }else{
      *rx0 = 0.0;
      *ry0 = 0.0;
    }
  }else if (strcmp(shape_set,"yosh")==0){
    if (shid == 39){
      *rx0 = 0.4;
      *ry0 = 0.4;
    }else if (shid == 40){  // OK
      *rx0 = -0.4;
      *ry0 =  0.4;
    }else if (shid == 41){  // OK
      *rx0 =  0.0;
      *ry0 =  0.2;
    }else if (shid == 42){  // OK
      *rx0 =  0.0;
      *ry0 =  0.2;
    }else{
      *rx0 = 0.0;
      *ry0 = 0.0;
    }
  }else if (strcmp(shape_set,"yasm")==0){
    *rx0 = 0.0;
    *ry0 = -0.2;
  }else{
    *rx0 = 0.0;
    *ry0 = 0.0;
  }

  if (theta != 0.0){
    xx[0] = *rx0;
    yy[0] = *ry0;
    farray_rotate_coords(theta,xx,yy,1);
    *rx0 = xx[0];
    *ry0 = yy[0];
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                           STM_PASUPATHY_SPLINE_EQN                        */
/*                                                                           */
/*  Given the 4 control coordinates, and the contour parameter 't', return   */
/*  the values of the shape coords at position 't' between 'p2' and 'p3'.    */
/*                                                                           */
/*****************************************************************************/
double stm_pasupathy_spline_eqn(p1,p2,p3,p4,t)
     double p1,p2,p3,p4;  // Coordinates (x or y) of four successive points
     double t;            // value of contour parameter, from 0..1
{
  double t3,t2,v;
  
  t2 = t*t;
  t3 = t2*t;

  v = (p1*(    -t3 + 3.0*t2 - 3.0*t + 1.0) +
       p2*( 3.0*t3 - 6.0*t2         + 4.0) +
       p3*(-3.0*t3 + 3.0*t2 + 3.0*t + 1.0) +
       p4*      t3)/6.0;

  return v;
}
/**************************************-**************************************/
/*                                                                           */
/*                         STM_PASUPATHY_SHAPE_GET_COORDS                    */
/*                                                                           */
/*  Get a set of finely spaced points for the shape.                         */
/*                                                                           */
/*****************************************************************************/
void stm_pasupathy_shape_get_coords(mylogf,shape_set,shid,theta,normflag,
				    px,rx,ry,rn)
     char *mylogf;  // Log file, or NULL
     char *shape_set;  // "pc2001", "yasm"
     int shid;      // Shape ID
     float theta;   // Rotation (deg)
     int normflag;  // 0-raw coords; 1-recenter
     int px;        // Point expansion: # of fine points / control point
     float **rx;    // Return x-coords [*rn]
     float **ry;    // Return y-coords [*rn]
     int    *rn;    // Return number of points
{
  int i,j,k;
  int sn,fn;
  float *sx,*sy,*tsx,*tsy,*fx,*fy,tmin,tmax;
  double x1,x2,x3,x4,y1,y2,y3,y4,t;

  // Get the control coords
  if (strcmp(shape_set,"yasm")==0)
    stm_pasupathy_yasm_get_control(mylogf,shid,&sx,&sy,&sn);
  else if (strcmp(shape_set,"yosh")==0)
    stm_pasupathy_yosh_get_control(mylogf,shid,&sx,&sy,&sn);
  else
    stm_pasupathy_shape_get_control(mylogf,shid,&sx,&sy,&sn);

  if (normflag == 1){
    //
    //  Recenter the shape
    //
    get_min_max_farray(sx,sn,&tmin,&tmax);
    add_const_farray(sx,sn,-(tmax+tmin)/2.0);

    get_min_max_farray(sy,sn,&tmin,&tmax);
    add_const_farray(sy,sn,-(tmax+tmin)/2.0);
  }

  /*
    mylog_ival(mylogf," (stm_util)  sn =",sn);
    mylog_ival(mylogf," (stm_util)  shid =",shid);
    mylog_ival(mylogf," (stm_util)  px =",px);
    mylog_ival(mylogf," (stm_util)  normflag =",normflag);
    mylog_fval(mylogf," (stm_util)  theta =",theta);
    mylog_cval(mylogf," (stm_util)  set =",shape_set);
  */

  // Rotate by the angle theta
  farray_rotate_coords(theta,sx,sy,sn);

  // Make arrays with first three points appended to end, for easy looping
  tsx = get_farray(sn+2);
  tsy = get_farray(sn+2);
  for(i=0;i<sn;i++){
    tsx[i] = sx[i];
    tsy[i] = sy[i];
  }

  tsx[sn  ] = sx[1];  // The first point, [0], is already repeated at end
  tsx[sn+1] = sx[2];

  tsy[sn  ] = sy[1];
  tsy[sn+1] = sy[2];


  // Create finely sampled points
  fn = px * (sn-1);  // First control point is repeated at end
  fx = get_farray(fn);
  fy = get_farray(fn);

  k = 0;
  for(i=0;i<(sn-1);i++){

    x1 = tsx[i];
    x2 = tsx[i+1];
    x3 = tsx[i+2];
    x4 = tsx[i+3];

    y1 = tsy[i];
    y2 = tsy[i+1];
    y3 = tsy[i+2];
    y4 = tsy[i+3];

    for(j=0;j<px;j++){
      t = (double)j/(double)px;
      fx[k] = stm_pasupathy_spline_eqn(x1,x2,x3,x4,t);
      fy[k] = stm_pasupathy_spline_eqn(y1,y2,y3,y4,t);
      k += 1;
    }
  }

  *rx = fx;
  *ry = fy;
  *rn = fn;

  myfree(tsx);
  myfree(tsy);
  myfree(sx);
  myfree(sy);
}
/**************************************-**************************************/
/*                                                                           */
/*                             STM_PASUPATHY_SHAPE                           */
/*                                                                           */
/*  Write the specified shape into the 'data' frame.                         */
/*                                                                           */
/*****************************************************************************/
void stm_pasupathy_shape(mylogf,data,xn,yn,xc,yc,shape_set,shid,theta,size,
			 fgval,bgval,aa_sd,fill_flag,norm_flag)
     char *mylogf;       // Log file, or NULL
     float **data;       // Write shape into 2D array [xn][yn]
                         //   NOTE: this is already filled w/ 'bgval'
     int xn,yn;          // Size of frame (pix)
     float xc,yc;        // Center of shape in [0..xn-1 , 0..yn-1]
     char *shape_set;    // "pc2001", "yasm"
     int shid;           // Shape ID
     float theta;        // Rotation of shape (deg) [0..360)
     float size;         // No. of pixels corresponding to shape 2 diam (pix)
     float fgval;        // Value to write inside shape
     float bgval;        // Background value (outside of shape)
     float aa_sd;        // SD for anti-aliasing (0.0 for none)
     int fill_flag;      // 1-fill, 0-outline
     int norm_flag;      // 0-none, 1-make area 1.0
{
  int sn,normflag,px;
  float *sx,*sy,sz_scale,x0,y0;
  char ggstr[SLEN];

  if (0){
    //
    //  Wyeth - this code was used to make comparisons, before .h file
    //  was finalized.  This can be removed.
    //

    printf("WYETH DEBUG ==== CALLING COMPARE\n");
    stm_pasupathy_shape_compare("yosh",1,"pc2001",32);
    stm_pasupathy_shape_compare("yosh",2,"pc2001",33);
    stm_pasupathy_shape_compare("yosh",3,"pc2001",34);
    stm_pasupathy_shape_compare("yosh",4,"pc2001",35);
    stm_pasupathy_shape_compare("yosh",5,"pc2001",36);

    // These are not the same, control points are rotated, and different
    //printf("6/51: ");
    //stm_pasupathy_shape_compare("yosh",6,"pc2001",51);
    //printf("7/48: ");
    //stm_pasupathy_shape_compare("yosh",7,"pc2001",48);

    stm_pasupathy_shape_compare("yosh",8,"pc2001",24);
    stm_pasupathy_shape_compare("yosh",9,"pc2001",25);
    stm_pasupathy_shape_compare("yosh",10,"pc2001",26);
    stm_pasupathy_shape_compare("yosh",11,"pc2001",27);
    stm_pasupathy_shape_compare("yosh",12,"pc2001",28);
    stm_pasupathy_shape_compare("yosh",13,"pc2001",29);
    stm_pasupathy_shape_compare("yosh",14,"pc2001",30);
    stm_pasupathy_shape_compare("yosh",15,"pc2001",31);

    //  These are not the same, the "yosh" shape is less rounded.
    stm_pasupathy_shape_compare("yosh",16,"pc2001",40);

    //  These are not the same, the "yosh" shape is less rounded.
    stm_pasupathy_shape_compare("yosh",21,"pc2001",43);

    stm_pasupathy_shape_compare("yosh",22,"pc2001",19);
    stm_pasupathy_shape_compare("yosh",23,"pc2001",16);
    stm_pasupathy_shape_compare("yosh",24,"pc2001",17);
    stm_pasupathy_shape_compare("yosh",25,"pc2001",18);
    stm_pasupathy_shape_compare("yosh",26,"pc2001",21);
    stm_pasupathy_shape_compare("yosh",27,"pc2001",22);
    stm_pasupathy_shape_compare("yosh",28,"pc2001",23);
    stm_pasupathy_shape_compare("yosh",29,"pc2001",50);
    stm_pasupathy_shape_compare("yosh",30,"pc2001",38);

    stm_pasupathy_shape_compare("yosh",32,"pc2001",41);

    stm_pasupathy_shape_compare("yosh",36,"pc2001",5);
    stm_pasupathy_shape_compare("yosh",37,"pc2001",6);
    stm_pasupathy_shape_compare("yosh",38,"pc2001",7);
    stm_pasupathy_shape_compare("yosh",39,"pc2001",10);
    stm_pasupathy_shape_compare("yosh",40,"pc2001",11);
    stm_pasupathy_shape_compare("yosh",41,"pc2001",14);
    stm_pasupathy_shape_compare("yosh",42,"pc2001",15);

    exit(0);
  }


  if (strcmp(shape_set,"yasm")==0){
    sz_scale = size / (1.7 * 1.85);  // Adjust for diam of BigCirc
    px = (int)(size / 4);  // How many fine points to make for each shape point
  }else{
    sz_scale = size / (pasu_shape_02[5] * 1.85);  // Adjust for diam of BigCirc
    px = (int)(size / 2);  // How many fine points to make for each shape point
  }

  //sprintf(ggstr,"    %d fine points per control point\n",px);
  //mylog(mylogf,ggstr);

  normflag = 0;  // 0-no norm.  WYETH unclear what we want to do here.

  stm_pasupathy_shape_get_coords(mylogf,shape_set,shid,theta,normflag,px,
				 &sx,&sy,&sn);

  // Coordinates of a point known to be on the inside of the shape,
  // w.r.t. 'sx' and 'sy', i.e., before scale and translation
  stm_pasupathy_shape_get_inside_point(shape_set,shid,theta,&x0,&y0);

  if (fill_flag == 1){
    stm_shape_draw_fill(mylogf,data,xn,yn,sx,sy,sn,sz_scale,xc,yc,fgval,bgval,
			x0,y0,aa_sd);
  }else{
    stm_shape_draw_contour_antialias(data,xn,yn,sx,sy,sn,sz_scale,xc,yc,fgval,
				     aa_sd);
  }

  if (norm_flag == 0){
    ; // Do nothing
  }else if (norm_flag == 1){
    if (bgval != 0.0)
      mylog(mylogf,"*** STM_PASUPATHY_SHAPE: norm interacts with bgval");
    norm_area_2d_farray(data,xn,yn,1.0);  // Make the area 1.0
  }else{
    exit_error("STM_PASUPATHY_SHAPE","Bad value for 'norm_flag'");
  }

  myfree(sx);
  myfree(sy);
}
/**************************************-**************************************/
/*                                                                           */
/*                               STM_GABOR_BRITT                             */
/*                                                                           */
/*  Provide random positions and times for gabor patches.                    */
/*                                                                           */
/*****************************************************************************/
void stm_gabor_britt(mylogf,npos,npatch,seed,soff,nstim,rposlist)
     char *mylogf;      // Log file, or NULL
     int npos;          // Number of spatial positions (e.g., 25 for 5x5 grid)
     int npatch;        // Number of Gabor patches shown at one time
     int seed;          // Randomization seed
     int soff;          // Sequence offset
     int nstim;         // Number of stimuli to return from random sequence
     int ***rposlist;   // [nstim][2] Positions for each epoch
{
  int i,j,k;
  int si,nunique,**ulist,**poslist,*shuf;

  //mylog(mylogf,"  STM_GABOR_BRITT\n");

  if (npatch == 2){
    nunique = (npos * npos - npos) / 2;
    ulist = get_2d_iarray(nunique,2);
    k = 0;
    for(i=0;i<npos;i++){
      for(j=0;j<npos;j++){
	if (i < j){
	  ulist[k][0] = i;
	  ulist[k][1] = j;
	  k += 1;
	}
      }
    }
    if (k != nunique)
      exit_error("STM_GABOR_BRITT","This should not happen");


    if ((soff < 0) || (soff >= nunique))
      exit_error("STM_GABOR_BRIT","Bad value for 'ranseq_offset'");

    poslist = get_2d_iarray(nstim,2);
    shuf = get_shuffle_index(nunique,5,seed);

    k = soff;
    for(i=0;i<nstim;i++){
      si = shuf[k];
      poslist[i][0] = ulist[si][0];
      poslist[i][1] = ulist[si][1];
      k += 1;
      if (k >= nunique)  // Wrap around to the beginning if needed
	k = 0;
    }

    myfree(shuf);
    free_2d_iarray(ulist,nunique);

  }else if (npatch == 1){

    nunique = npos;  // I.e., the number of positions

    if ((soff < 0) || (soff >= nunique)){
      printf("  *** soff = %d, but should be less than 'nunique'.\n",soff);
      printf("  *** nunique = %d\n",nunique);
      exit_error("STM_GABOR_BRIT","Bad value for 'ranseq_offset'");
    }

    poslist = get_2d_iarray(nstim,2);

    shuf = get_shuffle_index(nunique,5,seed);

    k = soff;
    for(i=0;i<nstim;i++){
      si = shuf[k];
      poslist[i][0] = si;
      poslist[i][1] = -1;
      k += 1;
      if (k >= nunique)  // Wrap around to the beginning if needed
	k = 0;
    }

    myfree(shuf);
  }
  *rposlist = poslist;
}
/**************************************-**************************************/
/*                                                                           */
/*                           STM_GABOR_BRITT_TIMING                          */
/*                                                                           */
/*  Return the timing values for presentation of one stimulus epoch, and     */
/*  the total number of stimulus epochs per trial.                           */
/*                                                                           */
/*****************************************************************************/
void stm_gabor_britt_timing(durfc,durmp,dursp,tscale,zn,t0,tsamp,n_epoch,
			    rn_fc,rn_rmp,rn_sep,rn_tot,rnstim)
     float durfc;       // (s) duration at full contrast
     float durmp;       // (s) duration of ramp
     float dursp;       // (s) additional separation
     float tscale;      // (s/frame)
     int   zn;          // (frames) duration of full trial
     float t0;          // (s) start time of stimulus
     int   tsamp;       // Time units per stimulus frame
     int   n_epoch;     // Desired number of stimuli, or '0' for max possible
     int *rn_fc;        // frames at full contrast
     int *rn_rmp;       // frames of ramp (contrast between 0 and full)
     int *rn_sep;       // frames of additional separation
     int *rn_tot;       // frames total assigned to one stimulus epoch
     int *rnstim;       // Total number of stimuli to show
{
  int t0i;

  *rn_fc  = durfc / tscale;      // Frames at full contrast
  *rn_rmp = durmp / tscale - 1;  // Frames to ramp up
  *rn_sep = dursp / tscale;      // Additional frames to separate stimuli
  if (*rn_rmp < 0)
    *rn_rmp = 0;

  *rn_tot = *rn_fc + 2 * *rn_rmp + *rn_sep + 2;

  t0i = stm_get_ti0(zn,tscale,t0,tsamp);  // Time index at which stim starts

  *rnstim = (zn - t0i) / *rn_tot;  // Maximum possible epochs per trial

  if (n_epoch > *rnstim)
    exit_error("STM_GABOR_BRITT_TIMING","'n_epoch' is too large");
  else if (n_epoch > 0){
    *rnstim = n_epoch;
  }
}
