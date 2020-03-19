/*****************************************************************************/
/*                                                                           */
/*  noise_util.c                                                             */
/*  Wyeth Bair                                                               */
/*                                                                           */
/*  This file contains routines related to the analysis of data from the     */
/*  velocity noise stimuli used in experiments in the Movshon lab on V1 and  */
/*  MT neurons.                                                              */
/*                                                                           */
/*  Much of this code was copied from "noise_anal.c".                        */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "my_util.h"
#include "misc_util.h"
#include "carray_util.h"
#include "iarray_util.h"
#include "farray_util.h"
#include "myrand_util.h"
#include "spike_util.h"
#include "ndata_util.h"
#include "fft_util.h"
#include "stm_util.h"
#include "ndata.h"

#include "rand.h"        // WYETH - "tap" arrays moved to .h file, 17/nov/10


/**************************************-**************************************/
/*                                                                           */
/*                              GET_MSEQ_TAP_ARRAY                           */
/*                                                                           */
/*  Return a pointer to the appropriate tap array, and the length.           */
/*                                                                           */
/*****************************************************************************/
void get_mseq_tap_array(morder,ra,rn)
     int morder,**ra,*rn;
{
  *rn = 40; // Default
  if (morder==8){
    *ra = sstap08;
    *rn = 16;
  }else if (morder==9)
    *ra = sstap09;
  else if (morder==10)
    *ra = sstap10;
  else if (morder==11)
    *ra = sstap11;
  else if (morder==12)
    *ra = sstap12;
  else if (morder==13)
    *ra = sstap13;
  else if (morder==14)
    *ra = sstap14;
  else if (morder==15)
    *ra = sstap15;
  else if (morder==16)
    *ra = sstap16;
}
/**************************************-**************************************/
/*                                                                           */
/*                           NOISE_UTIL_GET_SEED_ARRAY                       */
/*                                                                           */
/*  Return a pointer to the seed list.                                       */
/*                                                                           */
/*****************************************************************************/
void noise_util_get_seed_array(ra,rn)
     int **ra,*rn;
{
  *ra = ssseed40;
  *rn = 40;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_VALUE_OPPOSITE_TAP_VALUE                     */
/*                                                                           */
/*  Return the tap value that is "opposite" (assuming the n values           */
/*  were arranged in a circle) to the specified value.                       */
/*                                                                           */
/*****************************************************************************/
int get_value_opposite_tap_value(tapval,morder)
     int tapval,morder;
{
  int k,kopp,*taparray,tan;

  get_mseq_tap_array(morder,&taparray,&tan); /* Get array of correct order */
  k = get_index_search_iarray(taparray,tan,tapval); /* Get index of tapval. */
  kopp = (k + tan/2) % tan; /* Compute index of value opposite k. */
  return taparray[kopp];
}
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
/*                             NOISE_UTIL_RAN2                               */
/*                                                                           */
/*  I changed the "long" to "int" so this works the same on DEC Alpha and    */
/*  Sun.  1997.10.17.                                                        */
/*                                                                           */
/*  Long period (>2*10^18) random number generator of L'Ecuyer with Bays-    */
/*  Durham shuffle and added safeguards.  Returns a uniform deviate between  */
/*  0.0 and 1.0 (exclusive of the endpoint values).  Call with idum a        */
/*  negative integer to initialize; thereafter, do not alter idum between    */
/*  successive deviates in a sequence.  RNMX should approximate the largest  */
/*  floating value that is less than 1.                                      */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float noise_util_ran2(idum)
     int *idum;
     /*long *idum;*/
{
  int j;
  /*long k;*/
  int k;
  /***  static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];***/
  static int idum2=123456789;
  static int iy=0;
  static int iv[NTAB];
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
/*                            NOISE_UTIL_02_RAN2                             */
/*                                                                           */
/*  I changed the "long" to "int" so this works the same on DEC Alpha and    */
/*  Sun.  1997.10.17.                                                        */
/*                                                                           */
/*  Long period (>2*10^18) random number generator of L'Ecuyer with Bays-    */
/*  Durham shuffle and added safeguards.  Returns a uniform deviate between  */
/*  0.0 and 1.0 (exclusive of the endpoint values).  Call with idum a        */
/*  negative integer to initialize; thereafter, do not alter idum between    */
/*  successive deviates in a sequence.  RNMX should approximate the largest  */
/*  floating value that is less than 1.                                      */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float noise_util_02_ran2(idum)
     int *idum;
     /*long *idum;*/
{
  int j;
  /*long k;*/
  int k;
  /***  static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];***/
  static int idum2=123456789;
  static int iy=0;
  static int iv[NTAB];
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
/*                              GASDEV  (p288-290)                           */
/*                                                                           */
/*  Returns a normally distributed deviate with zero mean and unit variance  */
/*  using ran1(idum) as the source of uniform deviates.                      */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float noise_util_gasdev(idum)
     int *idum;
     /*long *idum;*/ /* Changed to make it the same on the DEC Alpha */
{
  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;
  
  /*if (iset == 0){*/
  if ((iset == 0) || (*idum <= 0)){  /* Wyeth, reset state if reseeded */
    do {
      v1=2.0*noise_util_ran2(idum)-1.0;
      v2=2.0*noise_util_ran2(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  }else{
    iset=0;
    return gset;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                        GET_STIMULUS_WORD_INDEX_FORMAT                     */
/*                                                                           */
/*  Recode the stimulus into work indices.                                   */
/*  Eg, a binary stimulus would be recoded into two letter words using       */
/*  0,1,2,3 at each location in the array.  The last location would be -1,   */
/*  since the last location has a word without a second letter.              */
/*                                                                           */
/*  Binary stimuli must be + and - values.                                   */
/*  Ternary stimuli must be +, 0, - values.                                  */
/*                                                                           */
/*****************************************************************************/
void get_stimulus_word_index_format(data,ntr,n,sord,swn,rwif,rwifn)
     float **data;  /* Stimulus data */
     int ntr;       /* Number of trials */
     int *n;        /* Duration of trial (allowed to vary across trials) */
     int sord;      /* Stimulus order */
     int swn;       /* Letters per stimulus word */
     int ***rwif,**rwifn; /* Returned stimulus word IDs */
{
  int i,j,k;
  int **wif,*wifn,wpt,wi,base,base0;

  if ((sord < 2) || (sord > 3))
    exit_error("GET_STIMULUS_WORD_INDEX_FORMAT","Order out of range");
  if (swn < 1)
    exit_error("GET_STIMULUS_WORD_INDEX_FORMAT","Word length < 1");

  wif = (int **)myalloc(ntr*sizeof(int *));
  wifn = (int *)myalloc(ntr*sizeof(int));

  base0 = 1;
  for(i=0;i<swn;i++)
    base0 *= sord;
  /*printf("base0 = %d\n",base0);*/

  for(i=0;i<ntr;i++){ /* For each trial */
    if (swn > n[i])
      exit_error("GET_STIMULUS_WORD_INDEX_FORMAT","Word length > n");
    wpt = (int)(n[i]-swn+1);
    wif[i] = (int *)myalloc(wpt*sizeof(int));
    wifn[i] = wpt;
    
    for(j=0;j<wpt;j++){ /* For each start location of a stimulus word */
      base = base0/sord;
      wi = 0;
      for(k=0;k<swn;k++){
	/*printf(" %2d",(int)data[i][j+k]);*/
	if (sord == 2){
	  if (data[i][j+k] > 0.0)
	    wi += base;
	}else if (sord == 3){
	  if (data[i][j+k] > 0.0)
	    wi += 2*base;
	  else if (data[i][j+k] == 0.0)
	    wi += 1*base;
	}
	base /= sord;
      }
      /*printf(" %3d\n",wi);*/
      wif[i][j] = wi;
    }
  }
  /*
    append_iarray_plot("wif.pl","wif",wif[0],wifn[0],1);
    append_farray_plot("wif.pl","stim",data[0],n[0],1);*/

  *rwif = wif; *rwifn = wifn;
}
/**************************************-**************************************/
/*                                                                           */
/*                            CONVERT_JUMP_TO_POSITION                       */
/*                                                                           */
/*  x0 - starting phase                                                      */
/*  xjump - maxjump                                                          */
/*  xtot - total, typically 1024                                             */
/*                                                                           */
/*****************************************************************************/
void convert_jump_to_position(data,n,x0,jumpsize,period)
     float *data; /* Delta, or jump */
     int n,x0,jumpsize,period;
{
  int i,k;
  int t;

  t = x0;
  for(i=0;i<n;i++){
    t += my_rint(data[i]);
    if (t >= period)
      t -= period;
    else if (t < 0)
      t += period;
    k = t/jumpsize;
    /***
      if (jumpsize == 256)
      printf("t,k = %d %d\n",t,k);***/
    data[i] = (float)k;
  }

  /**printf("*** WYETH - changing jump seq.\n");
     for(i=0;i<n;i++)
     data[i] = ((i/8)%255)%period;**/
}
/**************************************-**************************************/
/*                                                                           */
/*                       CONVERT_POSITION_TO_TOWARD_AWAY                     */
/*                                                                           */
/*  Convert the stimulus from position values to ternary (-1,0,+1) values    */
/*  depending on whether each jump is toward, neutral, or away from the      */
/*  preferred stimulus position.                                             */
/*                                                                           */
/*  ppref - position number of preferred position                            */
/*  pnull - position number of null position                                 */
/*  npos - total number of positions                                         */
/*                                                                           */
/*****************************************************************************/
void convert_position_to_toward_away(data,n,ppref,pnull,npos)
     float *data; /* Delta, or jump */
     int n;
     float ppref,pnull;
     int npos;
{
  int i,j,k;
  float *table,dp,dn,dpn,fnpos,*temp;

  fnpos = (float)npos;

  table = (float *)myalloc(npos*sizeof(float));
  for(i=0;i<npos;i++){
    dp = get_circular_diff((float)i,ppref,fnpos);
    dn = get_circular_diff((float)i,pnull,fnpos);
    dpn = get_circular_interval_length_containing_point(ppref,pnull,i,fnpos);
    if (dp < dn)
      table[i] = (dpn - dp)/dpn;
    else
      table[i] = dn/dpn;
  }
  printf("npos = %d  pref=%f  null=%f\n",npos,ppref,pnull);
  append_farray_plot("zzz.table.pl","table",table,npos,1);
  /*append_farray_plot("zzz.table.pl","data",data,n,1);*/

  /*** WYETH - this is inefficient, can do it w/o temp ***/
  temp = (float *)myalloc(n*sizeof(float));
  temp[0] = 0.0;
  for(i=1;i<n;i++){
    j = data[i-1];
    k = data[i];
    if ((k < 0) || (k >= npos) || (j < 0) || (j >= npos)){
      printf("i=%d  j=%d  k=%d  npos = %d\n",i,j,k,npos);
      exit_error("CONVERT_POSITION_TO_TOWARD_AWAY","Illegal position value");
    }
    if (table[j] < table[k])
      temp[i] = 1.0;
    else if (table[j] > table[k])
      temp[i] = -1.0;
    else
      temp[i] = 0.0;
  }
  for(i=0;i<n;i++)
    data[i] = temp[i];
  myfree(temp);
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_STIM_SCREEN_Y_CONST                        */
/*                                                                           */
/*  Return the vertical position of the stimulus (center) on the screen.     */
/*                                                                           */
/*****************************************************************************/
int get_stim_screen_y_const(nd)
     struct ndata_struct *nd;
{
  int y,tot,flag_y,flag_apy,flag_pat1y,flag_scry;
  char *y_str,*apy_str,*pat1y_str,*scry_str;
  
  flag_y     = ndata_get_const_param_int(nd,"y",        &y_str);
  flag_apy   = ndata_get_const_param_int(nd,"apy",    &apy_str);
  flag_pat1y = ndata_get_const_param_int(nd,"pat1y",&pat1y_str);
  flag_scry  = ndata_get_const_param_int(nd,"scry",  &scry_str);

  tot = flag_y + flag_apy + flag_pat1y + flag_scry;

  y = -1; /* to avoid -Wall compiler warning */
  if (tot == 0)
    exit_error("GET_STIM_SCREEN_Y_CONST","Cannot find screen y position");
  else if (tot > 1)
    exit_error("GET_STIM_SCREEN_Y_CONST","Found multiple screen y positions");
  else{
    if (flag_y == 1){
      y = atoi(y_str); myfree(y_str);
    }else if (flag_apy == 1){
      y = atoi(apy_str); myfree(apy_str);
    }else if (flag_pat1y == 1){
      y = atoi(pat1y_str); myfree(pat1y_str);
    }else if (flag_scry == 1){
      y = atoi(scry_str); myfree(scry_str);
    }
  }

  return y;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_WN_STIM_PARAMS                            */
/*                                                                           */
/*  Return the stimulus sequence for the "k"the trial.  This routine covers  */
/*  the following p-files:                                                   */
/*    wnone100.p                                                             */
/*    wnmat100.p                                                             */
/*    wnall100.p                                                             */
/*                                                                           */
/*****************************************************************************/
void get_wn_stim_params(nd,k,awake_flag,rframerate,rjumpseed,rdtframe,rperiod,
			rdistrib,rmu,rsigma,rminjump,rmaxjump)
     struct ndata_struct *nd;
     int k,awake_flag;
     int *rframerate,*rjumpseed,*rdtframe,*rperiod,*rdistrib,*rmu,*rsigma;
     int *rminjump,*rmaxjump;
{
  int flag;

  flag = ndata_get_vc_param_int(nd,k,"framrate",rframerate);
  if (!flag) exit_error("GET_WN_STIM_PARAMS","Param framerate not found");
  flag = ndata_get_vc_param_int(nd,k,"jumpseed",rjumpseed);
  if (!flag) exit_error("GET_WN_STIM_PARAMS","Param jumpseed not found");
  flag = ndata_get_vc_param_int(nd,k,"dt",rdtframe);
  if (!flag) exit_error("GET_WN_STIM_PARAMS","Param dt not found");
  flag = ndata_get_vc_param_int(nd,k,"period",rperiod);
  if (!flag) exit_error("GET_WN_STIM_PARAMS","Param period not found");

  /*** This catches the "distrib" value for vbin2md.p ***/
  flag = ndata_get_vc_param_int(nd,k,"vdpflag",rdistrib);
  if (!flag)
    flag = ndata_get_vc_param_int(nd,k,"distrib",rdistrib);

  if (!flag) exit_error("GET_WN_STIM_PARAMS","Param distrib not found");
  flag = ndata_get_vc_param_int(nd,k,"mu",rmu);
  if (!flag) exit_error("GET_WN_STIM_PARAMS","Param mu not found");
  flag = ndata_get_vc_param_int(nd,k,"sigma",rsigma);
  if (!flag) exit_error("GET_WN_STIM_PARAMS","Param sigma not found");
  if (awake_flag){
    flag = ndata_get_vc_param_int(nd,k,"minjump",rminjump);
    if (!flag) exit_error("GET_WN_STIM_PARAMS","Param minjump not found");
  }else{
    flag = ndata_get_vc_param_int(nd,k,"min_jump",rminjump);
    if (!flag) exit_error("GET_WN_STIM_PARAMS","Param min_jump not found");
  }
  if (awake_flag){
    flag = ndata_get_vc_param_int(nd,k,"maxjump",rmaxjump);
    if (!flag) exit_error("GET_WN_STIM_PARAMS","Param maxjump not found");
  }else{
    flag = ndata_get_vc_param_int(nd,k,"max_jump",rmaxjump);
    if (!flag) exit_error("GET_WN_STIM_PARAMS","Param max_jump not found");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_WN_TWO_STIM_PARAMS                         */
/*                                                                           */
/*  Return the stimulus sequence for the "k"the trial.  This routine covers  */
/*  the following p-files:                                                   */
/*    vbin2md.p                                                              */
/*    vbin2.p                                                                */
/*    wntwo100.p ?                                                           */
/*                                                                           */
/*****************************************************************************/
void get_wn_two_stim_params(nd,k,expflag,rframerate,rjumpseed,rdtframe,rperiod,
			    rdistrib,rmu,rsigma,rmnj0,rmnj1,rmxj0,rmxj1)
     struct ndata_struct *nd;
     int k,expflag;
     int *rframerate,*rjumpseed,*rdtframe,*rperiod,*rdistrib,*rmu,*rsigma;
     int *rmnj0,*rmnj1,*rmxj0,*rmxj1;
{
  int flag;

  flag = ndata_get_vc_param_int(nd,k,"framrate",rframerate);
  if (!flag) exit_error("GET_WN_STIM_PARAMS","Param framerate not found");
  flag = ndata_get_vc_param_int(nd,k,"jumpseed",rjumpseed);
  if (!flag) exit_error("GET_WN_STIM_PARAMS","Param jumpseed not found");
  flag = ndata_get_vc_param_int(nd,k,"dt",rdtframe);
  if (!flag) exit_error("GET_WN_STIM_PARAMS","Param dt not found");
  flag = ndata_get_vc_param_int(nd,k,"period",rperiod);
  if (!flag) exit_error("GET_WN_STIM_PARAMS","Param period not found");

  /*** This catches the "distrib" value for vbin2md.p ***/
  flag = ndata_get_vc_param_int(nd,k,"vdpflag",rdistrib);
  if (!flag)
    flag = ndata_get_vc_param_int(nd,k,"distrib",rdistrib);
  if (!flag) exit_error("GET_WN_STIM_PARAMS","Param distrib not found");

  flag = ndata_get_vc_param_int(nd,k,"mu",rmu);
  if (!flag) exit_error("GET_WN_STIM_PARAMS","Param mu not found");
  flag = ndata_get_vc_param_int(nd,k,"sigma",rsigma);
  if (!flag) exit_error("GET_WN_STIM_PARAMS","Param sigma not found");
  if (expflag){
    flag = ndata_get_vc_param_int(nd,k,"minjump0",rmnj0);
    if (!flag) exit_error("GET_WN_STIM_PARAMS","Param minjump0 not found");
    flag = ndata_get_vc_param_int(nd,k,"minjump1",rmnj1);
    if (!flag) exit_error("GET_WN_STIM_PARAMS","Param minjump1 not found");
    flag = ndata_get_vc_param_int(nd,k,"maxjump0",rmxj0);
    if (!flag) exit_error("GET_WN_STIM_PARAMS","Param maxjump0 not found");
    flag = ndata_get_vc_param_int(nd,k,"maxjump1",rmxj1);
    if (!flag) exit_error("GET_WN_STIM_PARAMS","Param maxjump1 not found");
  }else{
    flag = ndata_get_vc_param_int(nd,k,"min_jmp0",rmnj0);
    if (!flag) exit_error("GET_WN_STIM_PARAMS","Param minjump0 not found");
    flag = ndata_get_vc_param_int(nd,k,"min_jmp1",rmnj1);
    if (!flag) exit_error("GET_WN_STIM_PARAMS","Param minjump1 not found");
    flag = ndata_get_vc_param_int(nd,k,"max_jmp0",rmxj0);
    if (!flag) exit_error("GET_WN_STIM_PARAMS","Param maxjump0 not found");
    flag = ndata_get_vc_param_int(nd,k,"max_jmp1",rmxj1);
    if (!flag) exit_error("GET_WN_STIM_PARAMS","Param maxjump1 not found");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_WY2PAT_STIM_PARAMS                         */
/*                                                                           */
/*  Feb 21, 2001 - adding 'rrpt' as return value for 'wydt.p'.               */
/*                                                                           */
/*****************************************************************************/
void get_wy2pat_stim_params(nd,k,awake_flag,rframerate,rseed,rdtframe,
			    rperiod,rdistrib,rmu,rsigma,rminv,rmaxv,rrpt)
     struct ndata_struct *nd;
     int k,awake_flag;
     int *rframerate,*rseed,*rdtframe,*rperiod,*rdistrib,*rmu,*rsigma;
     int *rminv,*rmaxv,*rrpt;
{
  int flag;
  char *pfile;

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0)
    exit_error("GET_WY2PAT_STIM_PARAMS","Const param pfile not found");

  flag = ndata_get_vc_param_int(nd,k,"framrate",rframerate);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param framerate not found");
  flag = ndata_get_vc_param_int(nd,k,"seed",rseed);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param seed not found");
  flag = ndata_get_vc_param_int(nd,k,"dt",rdtframe);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param dt not found");
  if ((strcmp(pfile,"wydt.p")==0)){
    flag = ndata_get_vc_param_int(nd,k,"rpt",rrpt);
    if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","rpt not found");
  }else
    *rrpt = 1;

  flag = ndata_get_vc_param_int(nd,k,"duration",rperiod);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param duration not found");
  flag = ndata_get_vc_param_int(nd,k,"distrib",rdistrib);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param distrib not found");

  /*** WYETH - 'wypnphal.p' is built on a wypnori parameter name framework */
  if ((strcmp(pfile,"wypnori.p")==0)||(strcmp(pfile,"wypnori3.p")==0)||
      (strcmp(pfile,"wypnoril.p")==0)||(strcmp(pfile,"wypnphal.p")==0)){
    flag = ndata_get_vc_param_int(nd,k,"order",rmu);
    if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param order not found");
    flag = ndata_get_vc_param_int(nd,k,"bgflag",rsigma);
    if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param bgflag not found");
    flag = ndata_get_vc_param_int(nd,k,"orinull",rminv);
    if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param orinull not found");
    flag = ndata_get_vc_param_int(nd,k,"oripref",rmaxv);
    if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param oripref not found");
    *rminv = -1; /* Don't want the ori's */
    *rmaxv = 1;
  }else if ((strcmp(pfile,"wypnpha.p")==0)||(strcmp(pfile,"wyphas8.p")==0)||
	    (strcmp(pfile,"wypnpha3.p")==0)||(strcmp(pfile,"wyph3.p")==0)||
	    (strcmp(pfile,"wyphasms.p")==0)||
	    (strcmp(pfile,"wyphasm.p")==0)||(strcmp(pfile,"wyphasm8.p")==0)){
    flag = ndata_get_vc_param_int(nd,k,"order",rmu);
    if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param order not found");
    flag = ndata_get_vc_param_int(nd,k,"bgflag",rsigma);
    if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param bgflag not found");
    flag = ndata_get_vc_param_int(nd,k,"nullph",rminv);
    if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param nullph not found");
    flag = ndata_get_vc_param_int(nd,k,"prefph",rmaxv);
    if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param prefph not found");
    *rminv = -1; /* Don't want the ori's */
    *rmaxv = 1;
  }else if ((strcmp(pfile,"wy4st.p")==0)||(strcmp(pfile,"x4st.p")==0)||
	    (strcmp(pfile,"wy4stg.p")==0)||
	    (strcmp(pfile,"wy9st.p")==0)||
	    (strcmp(pfile,"wy4stp.p")==0)||
	    (strcmp(pfile,"wy4std.p")==0)||
	    (strcmp(pfile,"wy2flat.p")==0)|| /*** No mu or sigma in xfile ***/
	    (strcmp(pfile,"wy4sum.p")==0)){ /*** No mu or sigma in xfile ***/
    *rmu = *rsigma = -1;
    flag = ndata_get_vc_param_int(nd,k,"minj1",rminv);
    if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param minj1 not found");
    flag = ndata_get_vc_param_int(nd,k,"maxj1",rmaxv);
    if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param maxj1 not found");
  }else if ((strcmp(pfile,"wmcon.p")==0)||(strcmp(pfile,"wmco1.p")==0)||
	    (strcmp(pfile,"wmlum.p")==0)||(strcmp(pfile,"wmlu1.p")==0)||
	    (strcmp(pfile,"wmlu4.p")==0)||(strcmp(pfile,"wmco4.p")==0)||
	    (strcmp(pfile,"wmcon160.p")==0)||(strcmp(pfile,"wmco1160.p")==0)||
	    (strcmp(pfile,"wmlum160.p")==0)||(strcmp(pfile,"wmlu1160.p")==0)){
    /*** "sigma" was replaced by "modval" in the xfile. ***/
    flag = ndata_get_vc_param_int(nd,k,"mu",rmu);
    if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param mu not found");
    flag = ndata_get_vc_param_int(nd,k,"modval",rsigma);
    if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param modval not found");
    flag = ndata_get_vc_param_int(nd,k,"minj",rminv);
    if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param minj not found");
    flag = ndata_get_vc_param_int(nd,k,"maxj",rmaxv);
    if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param maxj not found");
  }else{
    flag = ndata_get_vc_param_int(nd,k,"mu",rmu);
    if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param mu not found");
    flag = ndata_get_vc_param_int(nd,k,"sigma",rsigma);
    if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param sigma not found");
    flag = ndata_get_vc_param_int(nd,k,"minj",rminv);
    if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param minj not found");
    flag = ndata_get_vc_param_int(nd,k,"maxj",rmaxv);
    if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param maxj not found");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_WY2PAT_STIM_PARAMS_02                        */
/*                                                                           */
/*****************************************************************************/
void get_wy2pat_stim_params_02(nd,k,awake_flag,rframerate,rseed,rperiod,
			       rdt1,rrpt1,rdist1,rmu1,rsig1,rmin1,rmax1,
			       rdt2,rrpt2,rdist2,rmu2,rsig2,rmin2,rmax2)
     struct ndata_struct *nd;
     int k,awake_flag;
     int *rframerate,*rseed,*rperiod;
     int *rdt1,*rrpt1,*rdist1,*rmu1,*rsig1,*rmin1,*rmax1;
     int *rdt2,*rrpt2,*rdist2,*rmu2,*rsig2,*rmin2,*rmax2;
{
  int flag;
  char *pfile;

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0)
    exit_error("GET_WY2PAT_STIM_PARAMS_02","Const param pfile not found");

  flag = ndata_get_vc_param_int(nd,k,"framrate",rframerate);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_02","framerate not found");
  flag = ndata_get_vc_param_int(nd,k,"seed",rseed);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_02","seed not found");
  flag = ndata_get_vc_param_int(nd,k,"duration",rperiod);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_02","period not found");

  flag = ndata_get_vc_param_int(nd,k,"dt1",rdt1);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_02","dt1 not found");
  flag = ndata_get_vc_param_int(nd,k,"rpt1",rrpt1);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_02","rpt1 not found");
  flag = ndata_get_vc_param_int(nd,k,"distrib1",rdist1);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_02","distrib1 not found");
  flag = ndata_get_vc_param_int(nd,k,"mu1",rmu1);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_02","mu1 not found");
  flag = ndata_get_vc_param_int(nd,k,"sigma1",rsig1);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_02","sigma1 not found");
  flag = ndata_get_vc_param_int(nd,k,"minj1",rmin1);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_02","minj1 not found");
  flag = ndata_get_vc_param_int(nd,k,"maxj1",rmax1);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_02","maxj1 not found");

  flag = ndata_get_vc_param_int(nd,k,"dt2",rdt2);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_02","dt2 not found");
  flag = ndata_get_vc_param_int(nd,k,"rpt2",rrpt2);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_02","rpt2 not found");
  flag = ndata_get_vc_param_int(nd,k,"distrib2",rdist2);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_02","distrib2 not found");
  flag = ndata_get_vc_param_int(nd,k,"mu2",rmu2);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_02","mu2 not found");
  flag = ndata_get_vc_param_int(nd,k,"sigma2",rsig2);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_02","sigma2 not found");
  flag = ndata_get_vc_param_int(nd,k,"minj2",rmin2);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_02","minj2 not found");
  flag = ndata_get_vc_param_int(nd,k,"maxj2",rmax2);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_02","maxj2 not found");
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_WY2PAT_STIM_PARAMS_03                        */
/*                                                                           */
/*****************************************************************************/
void get_wy2pat_stim_params_03(nd,k,rframerate,rseed,rdtframe,
			       rperiod,rdistrib,rmu,rsigma,rminv,rmaxv,
			       rwalkflag,rsteptype,ramp0)
     struct ndata_struct *nd;
     int k;
     int *rframerate,*rseed,*rdtframe,*rperiod,*rdistrib,*rmu,*rsigma;
     int *rminv,*rmaxv,*rwalkflag,*rsteptype,*ramp0;
{
  int flag,stimtype;
  char *pfile;
  
  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0)
    exit_error("GET_WY2PAT_STIM_PARAMS","Const param pfile not found");

  flag = ndata_get_vc_param_int(nd,k,"stimtype",&stimtype);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param stimtype not found");
  *rwalkflag = 23-stimtype;
  
  flag = ndata_get_vc_param_int(nd,k,"framrate",rframerate);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param framerate not found");
  flag = ndata_get_vc_param_int(nd,k,"seed",rseed);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param seed not found");
  flag = ndata_get_vc_param_int(nd,k,"dt",rdtframe);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param dt not found");
  flag = ndata_get_vc_param_int(nd,k,"duration",rperiod);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param period not found");
  flag = ndata_get_vc_param_int(nd,k,"distrib",rdistrib);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param distrib not found");

  flag = ndata_get_vc_param_int(nd,k,"mu",rmu);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param mu not found");
  flag = ndata_get_vc_param_int(nd,k,"sigma",rsigma);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param sigma not found");
  flag = ndata_get_vc_param_int(nd,k,"minj",rminv);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param minj not found");
  flag = ndata_get_vc_param_int(nd,k,"maxj",rmaxv);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param maxj not found");

  flag = ndata_get_vc_param_int(nd,k,"steptype",rsteptype);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param steptype not found");
  flag = ndata_get_vc_param_int(nd,k,"amp0",ramp0);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS","Param amp0 not found");
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_WY2PAT_STIM_PARAMS_04                        */
/*                                                                           */
/*  For `wy2flat2.p' - used for ternary center or surround in LGN expt.      */
/*                                                                           */
/*****************************************************************************/
void get_wy2pat_stim_params_04(nd,k,awake_flag,rframerate,rseed,rperiod,
			       rdt1,rrpt1,rdist1,rpat1amp,
			       rdt2,rrpt2,rdist2,rpat2amp)
     struct ndata_struct *nd;
     int k,awake_flag;
     int *rframerate,*rseed,*rperiod;
     int *rdt1,*rrpt1,*rdist1,*rpat1amp;
     int *rdt2,*rrpt2,*rdist2,*rpat2amp;
{
  int flag;
  char *pfile;

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0)
    exit_error("GET_WY2PAT_STIM_PARAMS_04","Const param pfile not found");

  flag = ndata_get_vc_param_int(nd,k,"framrate",rframerate);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_04","framerate not found");
  flag = ndata_get_vc_param_int(nd,k,"seed",rseed);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_04","seed not found");
  flag = ndata_get_vc_param_int(nd,k,"duration",rperiod);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_04","period not found");

  flag = ndata_get_vc_param_int(nd,k,"dt1",rdt1);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_04","dt1 not found");
  flag = ndata_get_vc_param_int(nd,k,"rpt1",rrpt1);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_04","rpt1 not found");
  flag = ndata_get_vc_param_int(nd,k,"distrib1",rdist1);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_04","distrib1 not found");
  flag = ndata_get_vc_param_int(nd,k,"pat1amp",rpat1amp);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_04","pat1amp not found");

  flag = ndata_get_vc_param_int(nd,k,"dt2",rdt2);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_04","dt2 not found");
  flag = ndata_get_vc_param_int(nd,k,"rpt2",rrpt2);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_04","rpt2 not found");
  flag = ndata_get_vc_param_int(nd,k,"distrib2",rdist2);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_04","distrib2 not found");
  flag = ndata_get_vc_param_int(nd,k,"pat2amp",rpat2amp);
  if (!flag) exit_error("GET_WY2PAT_STIM_PARAMS_04","pat2amp not found");
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_WY2PAT_REL_STIM_PARAMS                       */
/*                                                                           */
/*****************************************************************************/
void get_wy2pat_rel_stim_params(nd,k,rframerate,rseed,rperiod,
				rtf1,rtf2,rdeltamin,rdeltamax,rncycles)
     struct ndata_struct *nd;
     int k,*rframerate,*rseed,*rperiod;
     int *rtf1,*rtf2,*rdeltamin,*rdeltamax,*rncycles;
{
  int flag;
  char *pfile;

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0)
    exit_error("GET_WY2PAT_STIM_PARAMS_02","Const param pfile not found");

  flag = ndata_get_vc_param_int(nd,k,"tf1",rtf1);
  if (!flag) exit_error("GET_WY2PAT_REL_STIM_PARAMS","tf1 not found");
  if (*rtf1 < 0) /*** 'tf1' is set to -tf for stationary (non-drifting) ***/
    *rtf1 *= -1;
  flag = ndata_get_vc_param_int(nd,k,"tf2",rtf2);
  if (!flag) exit_error("GET_WY2PAT_REL_STIM_PARAMS","tf2 not found");
  flag = ndata_get_vc_param_int(nd,k,"deltamin",rdeltamin);
  if (!flag) exit_error("GET_WY2PAT_REL_STIM_PARAMS","deltamin not found");
  flag = ndata_get_vc_param_int(nd,k,"deltamax",rdeltamax);
  if (!flag) exit_error("GET_WY2PAT_REL_STIM_PARAMS","deltamax not found");
  flag = ndata_get_vc_param_int(nd,k,"ncycles",rncycles);
  if (!flag) exit_error("GET_WY2PAT_REL_STIM_PARAMS","ncycles not found");

  flag = ndata_get_vc_param_int(nd,k,"framrate",rframerate);
  if (!flag) exit_error("GET_WY2PAT_REL_STIM_PARAMS","framrate not found");
  flag = ndata_get_vc_param_int(nd,k,"seed",rseed);
  if (!flag) exit_error("GET_WY2PAT_REL_STIM_PARAMS","seed not found");
  flag = ndata_get_vc_param_int(nd,k,"duration",rperiod);
  if (!flag) exit_error("GET_WY2PAT_REL_STIM_PARAMS","duration not found");
}
/**************************************-**************************************/
/*                                                                           */
/*                        GET_WY2PAT_SIMPLE_STIM_PARAMS                      */
/*                                                                           */
/*****************************************************************************/
void get_wy2pat_simple_stim_params(nd,k,rframerate,rseed,rperiod,rtf1,rtf2,
				   rnphase,rncycles,rpat1amp,rpat2amp)
     struct ndata_struct *nd;
     int k,*rframerate,*rseed,*rperiod;
     int *rtf1,*rtf2,*rnphase,*rncycles,*rpat1amp,*rpat2amp;
{
  int flag,stimtype;
  char *pfile;

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0)
    exit_error("GET_WY2PAT_SIMPLE_STIM_PARAMS","Const param pfile not found");

  flag = ndata_get_vc_param_int(nd,k,"tf1",rtf1);
  if (!flag) exit_error("GET_WY2PAT_SIMPLE_STIM_PARAMS","tf1 not found");
  if (*rtf1 < 0) /*** 'tf1' is set to -tf for stationary (non-drifting) ***/
    *rtf1 *= -1;
  flag = ndata_get_vc_param_int(nd,k,"tf2",rtf2);
  if (!flag) exit_error("GET_WY2PAT_SIMPLE_STIM_PARAMS","tf2 not found");
  flag = ndata_get_vc_param_int(nd,k,"nphase",rnphase);
  if (!flag) exit_error("GET_WY2PAT_SIMPLE_STIM_PARAMS","nphase not found");
  flag = ndata_get_vc_param_int(nd,k,"ncycles",rncycles);
  if (!flag) exit_error("GET_WY2PAT_SIMPLE_STIM_PARAMS","ncycles not found");

  flag = ndata_get_vc_param_int(nd,k,"pat1amp",rpat1amp);
  if (!flag) exit_error("GET_WY2PAT_SIMPLE_STIM_PARAMS","pat1amp not found");
  flag = ndata_get_vc_param_int(nd,k,"pat2amp",rpat2amp);
  if (!flag) exit_error("GET_WY2PAT_SIMPLE_STIM_PARAMS","pat2amp not found");
  if (strcmp(pfile,"x4simp.p")==0){
    flag = ndata_get_vc_param_int(nd,k,"stimtype",&stimtype);
    if (!flag) exit_error("GET_WY2PAT_SIMPLE_STIM_PARAM","stimtype not found");
    if (stimtype == 1)
      *rpat2amp = 0;
    else if (stimtype == 2)
      *rpat1amp = 0;
  }

  flag = ndata_get_vc_param_int(nd,k,"framrate",rframerate);
  if (!flag) exit_error("GET_WY2PAT_SIMPLE_STIM_PARAMS","framrate not found");
  flag = ndata_get_vc_param_int(nd,k,"seed",rseed);
  if (!flag) exit_error("GET_WY2PAT_SIMPLE_STIM_PARAMS","seed not found");
  flag = ndata_get_vc_param_int(nd,k,"duration",rperiod);
  if (!flag) exit_error("GET_WY2PAT_SIMPLE_STIM_PARAMS","duration not found");
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_WY2FLAT_REL_STIM_PARAMS                      */
/*                                                                           */
/*****************************************************************************/
void get_wy2flat_rel_stim_params(nd,k,rframerate,rseed,rperiod,
				 rdeltamin,rdeltamax,rdwell)
     struct ndata_struct *nd;
     int k,*rframerate,*rseed,*rperiod;
     int *rdeltamin,*rdeltamax,*rdwell;
{
  int flag;
  char *pfile;

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0)
    exit_error("GET_WY2PAT_STIM_PARAMS_02","Const param pfile not found");

  flag = ndata_get_vc_param_int(nd,k,"deltamin",rdeltamin);
  if (!flag) exit_error("GET_WY2PAT_REL_STIM_PARAMS","deltamin not found");
  flag = ndata_get_vc_param_int(nd,k,"deltamax",rdeltamax);
  if (!flag) exit_error("GET_WY2PAT_REL_STIM_PARAMS","deltamax not found");
  flag = ndata_get_vc_param_int(nd,k,"dwell",rdwell);
  if (!flag) exit_error("GET_WY2PAT_REL_STIM_PARAMS","dwell not found");

  flag = ndata_get_vc_param_int(nd,k,"framrate",rframerate);
  if (!flag) exit_error("GET_WY2PAT_REL_STIM_PARAMS","framrate not found");
  flag = ndata_get_vc_param_int(nd,k,"seed",rseed);
  if (!flag) exit_error("GET_WY2PAT_REL_STIM_PARAMS","seed not found");
  flag = ndata_get_vc_param_int(nd,k,"duration",rperiod);
  if (!flag) exit_error("GET_WY2PAT_REL_STIM_PARAMS","duration not found");
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_WY2PA_STIM_PARAMS                         */
/*                                                                           */
/*****************************************************************************/
void get_wy2pa_stim_params(nd,k,rframerate,rseed,rdur,rantidur,rprefdur,
			   rtoverlap)
     struct ndata_struct *nd;
     int k,*rframerate,*rseed,*rdur,*rantidur,*rprefdur,*rtoverlap;
{
  int flag;
  char *pfile;

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0)
    exit_error("GET_WY2PA_STIM_PARAMS","Const param pfile not found");

  flag = ndata_get_vc_param_int(nd,k,"framrate",rframerate);
  if (!flag) exit_error("GET_WY2PA_STIM_PARAMS","framrate not found");
  flag = ndata_get_vc_param_int(nd,k,"seed",rseed);
  if (!flag) exit_error("GET_WY2PA_STIM_PARAMS","seed not found");
  flag = ndata_get_vc_param_int(nd,k,"duration",rdur);
  if (!flag) exit_error("GET_WY2PA_STIM_PARAMS","duration not found");
  flag = ndata_get_vc_param_int(nd,k,"antidur",rantidur);
  if (!flag) exit_error("GET_WY2PA_STIM_PARAMS","antidur not found");
  flag = ndata_get_vc_param_int(nd,k,"prefdur",rprefdur);
  if (!flag) exit_error("GET_WY2PA_STIM_PARAMS","prefdur not found");
  flag = ndata_get_vc_param_int(nd,k,"toverlap",rtoverlap);
  if (!flag) exit_error("GET_WY2PA_STIM_PARAMS","toverlap not found");
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_WYTFV_STIM_PARAMS                         */
/*                                                                           */
/*****************************************************************************/
void get_wytfv_stim_params(nd,k,rframerate,rseed,rdur,rzflag,rt1,rt2,rt3,rntf,
			   rmaxtf)
     struct ndata_struct *nd;
     int k,*rframerate,*rseed,*rdur,*rzflag,*rt1,*rt2,*rt3,*rntf,*rmaxtf;
{
  int flag;
  char *pfile;

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0)
    exit_error("GET_WYTFV_STIM_PARAMS","Const param pfile not found");

  flag = ndata_get_vc_param_int(nd,k,"framrate",rframerate);
  if (!flag) exit_error("GET_WYTFV_STIM_PARAMS","framrate not found");
  flag = ndata_get_vc_param_int(nd,k,"seed",rseed);
  if (!flag) exit_error("GET_WYTFV_STIM_PARAMS","seed not found");
  flag = ndata_get_vc_param_int(nd,k,"duration",rdur);
  if (!flag) exit_error("GET_WYTFV_STIM_PARAMS","duration not found");

  flag = ndata_get_vc_param_int(nd,k,"zflag",rzflag);
  if (!flag) exit_error("GET_WYTFV_STIM_PARAMS","zflag not found");
  flag = ndata_get_vc_param_int(nd,k,"t1",rt1);
  if (!flag) exit_error("GET_WYTFV_STIM_PARAMS","t1 not found");
  flag = ndata_get_vc_param_int(nd,k,"t2",rt2);
  if (!flag) exit_error("GET_WYTFV_STIM_PARAMS","t2 not found");
  flag = ndata_get_vc_param_int(nd,k,"t3",rt3);
  if (!flag) exit_error("GET_WYTFV_STIM_PARAMS","t3 not found");
  flag = ndata_get_vc_param_int(nd,k,"ntf",rntf);
  if (!flag) exit_error("GET_WYTFV_STIM_PARAMS","ntf not found");
  flag = ndata_get_vc_param_int(nd,k,"maxtf",rmaxtf);
  if (!flag) exit_error("GET_WYTFV_STIM_PARAMS","maxtf not found");
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_ADDOT_STIM_PARAMS                         */
/*                                                                           */
/*****************************************************************************/
void get_addot_stim_params(nd,k,rseed,rt1nstim,rt2nstim,rt1ori,rt2ori)
     struct ndata_struct *nd;
     int k,*rseed,*rt1nstim,*rt2nstim,*rt1ori,*rt2ori;
{
  int flag;
  char *pfile;

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0)
    exit_error("GET_WY2PAT_STIM_PARAMS_02","Const param pfile not found");

  flag = ndata_get_vc_param_int(nd,k,"seed",rseed);
  if (!flag) exit_error("GET_ADDOT_STIM_PARAMS","seed not found");
  flag = ndata_get_vc_param_int(nd,k,"t1nstim",rt1nstim);
  if (!flag) exit_error("GET_ADDOT_STIM_PARAMS","rt1nstim not found");
  flag = ndata_get_vc_param_int(nd,k,"t2nstim",rt2nstim);
  if (!flag) exit_error("GET_ADDOT_STIM_PARAMS","rt2nstim not found");
  flag = ndata_get_vc_param_int(nd,k,"t1ori",rt1ori);
  if (!flag) exit_error("GET_ADDOT_STIM_PARAMS","rt1ori not found");
  flag = ndata_get_vc_param_int(nd,k,"t2ori",rt2ori);
  if (!flag) exit_error("GET_ADDOT_STIM_PARAMS","rt2ori not found");
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_WN_ACTUAL_DT                            */
/*                                                                           */
/*  Return the time between frames in milliseconds for the given frame rate  */
/*  parameter.  These values were estimated from sync pulses in data files.  */
/*                                                                           */
/*****************************************************************************/
float get_wn_actual_dt(framerate,awake_flag)
     int framerate,awake_flag;
{
  float dt;

  /*** WYETH - new feature as of Nov 19, 2001 ***/
  /**** WYETH *****/ /* NOTE: check for 'actual_dt' in .nd file */
  /**** WYETH *****/ /* NOTE: check for 'actual_dt' in .nd file */
  /**** WYETH *****/ /* NOTE: check for 'actual_dt' in .nd file */
  /**** WYETH *****/ /* NOTE: check for 'actual_dt' in .nd file */

  if (awake_flag){
    if (framerate == 100)
      dt = 10.0105; /*** Awake lab. ***/
    else{
      printf("*** GET_WN_ACTUAL_DT: WARNING, Estimating dt.\n");
      dt = 1000.0/(double)framerate;
    }
  }else if (framerate == 50)
    dt = 19.902; /*** OK, based on 500 frame epoch. ***/
  else if (framerate == 53)
    dt = 18.6053; /*** This is good. ***/
  else if (framerate == 80)
    dt = 12.3522; /*** This is good. ***/
  else if (framerate == 100)
    dt = 10.0375; /*** This is good, ACUTE LAB - pep.  maybe 10.03733 ***/
  else if (framerate == 160)
    dt = 6.255625; /*** This is OK, ACUTE LAB - pep, (from 800 frames) ***/
  else{
    printf("*** GET_WN_ACTUAL_DT: WARNING, Estimating time between frames.\n");
    dt = 1000.0/(double)framerate;
  }
  return dt;
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_SGI_ACTUAL_DT                           */
/*                                                                           */
/*  Takes the framerate and returns the # of milliseconds per frame.         */
/*                                                                           */
/* NOTES: This is based on empirical testing and would change if the screen  */
/* resolution was changed, so these numbers should be updated if the system  */
/* is altered from its standard state.                                       */
/*                                                                           */
/*  hflag (hardware flag)                                                    */
/*   0 - Octane (circa Dec 2000)                                             */
/*                                                                           */
/*****************************************************************************/
float get_sgi_actual_dt(framerate,hflag)
     int framerate,hflag;
{
  float msperframe;

  /*** WYETH - new feature as of Nov 19, 2001 ***/
  /**** WYETH *****/ /* NOTE: check for 'actual_dt' in .nd file */
  /**** WYETH *****/ /* NOTE: check for 'actual_dt' in .nd file */
  /**** WYETH *****/ /* NOTE: check for 'actual_dt' in .nd file */

  /*** Gottfried Dec 2000, at 100Hz = 10.000358 ***/

  msperframe = -1.0; /* to avoid -Wall compiler warning */
  if (hflag == 0){
    
    if (framerate == 100) { /* for gottfried */
      msperframe = 10.000358;
    }else
      exit_error("GET_SGI_ACTUAL_DT","Unknown framerate");
    
    /*** WYETH - see read_pep.c ***/
    /*** The following value has been current for a while, Feb 2002 ***/
    /* dtstr = strdup("10.000275");*/ /* SGI Octane2 V10 (senta) */
    
  }else{ /*** Old numbers for O2 ***/
    if (framerate == 100) { /* for ortrud - 1024x731_100 */
      msperframe = 10.0157;
    }else if (framerate == 72) { /* for 1280x1024_72 - standard SGI */
      msperframe = 13.8485;
    }else if (framerate == 75) { /* for 1280x1024_75 - standard SGI */
      msperframe = 13.3448;
    }else {
      exit_error("GET_SGI_ACTUAL_DT",
		 "The ms per frame is unknown for that framerate");
    }
  }

  return msperframe;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_ACTUAL_SCREEN_DIMENSIONS                      */
/*                                                                           */
/*  Return the width and height of the screen in pixels.                     */
/*                                                                           */
/*****************************************************************************/
void get_actual_screen_dimensions(pf,framerate,flag,rw,rh)
     char pf[]; /* pfile */
     int framerate,flag,*rw,*rh;
{
  int w,h;

  if (framerate == 50){ /* wnmat050 */
    w = 880;
    h = 643;
  }else if (framerate == 160){ /* wymat160 */
    w = 784;
    h = 461;
  }else{
    w = 1024;
    h = 731;
  }

  *rw = w; *rh = h;
}
/**************************************-**************************************/
/*                                                                           */
/*                              GET_ACTUAL_LATENCY                           */
/*                                                                           */
/*  Return the time between the sync pulse indicating stimulus onset and     */
/*  the actual change in luminance on the video screen (as measured by a     */
/*  photometer).                                                             */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - See /u/u9/wyeth/analysis/latency_pep/ for latency estimates.           */
/*  - This routine does not apply to 'ori' like stimuli, for which the       */
/*    time to the top of the screen is 0~ms, and to middle is 4ms@100Hz.     */
/*  - This latency does *not* include the probably 2--6~ms delay associated  */
/*    with spike discrimination (ideally if delay on Bak is set to 0, then   */
/*    latency is time between initial slope/level trigger and window B       */
/*    (assuming dual window A+B criterion) which would be 1--2 ms.           */
/*                                                                           */
/*****************************************************************************/
int get_actual_latency(pf,dt)
     char pf[];
     float dt;
{
  int lat;

  if ((strcmp(pf,"MODEL.p")==0)||
      (strcmp(pf,"MODEL_wyamc.p")==0)||(strcmp(pf,"MODEL_new.p")==0)){
    lat = 0;
  }else if ((strcmp(pf, "wypnori.p")==0)||(strcmp(pf, "wypnpha.p")==0)||
	    (strcmp(pf,"wypnori3.p")==0)||(strcmp(pf,"wypnpha3.p")==0)||
	    (strcmp(pf,   "wy4st.p")==0)||(strcmp(pf,   "wyone.p")==0)||
	    (strcmp(pf,  "wyonem.p")==0)||
	    (strcmp(pf, "wyone23.p")==0)||(strcmp(pf,   "wy9st.p")==0)||
	    (strcmp(pf,  "wy4stg.p")==0)||(strcmp(pf,"wypnoril.p")==0)||
	    (strcmp(pf, "wy2flat.p")==0)||(strcmp(pf,"wypnphal.p")==0)||
	    (strcmp(pf,  "wy4sum.p")==0)||(strcmp(pf,"wyphstat.p")==0)||
	    (strcmp(pf,   "wycon.p")==0)||(strcmp(pf,  "wymtf5.p")==0)||
	    (strcmp(pf, "wymtf5p.p")==0)||(strcmp(pf,"wymtf4ph.p")==0)||
	    (strcmp(pf,   "wyall.p")==0)||
	    (strcmp(pf,  "wymtf9.p")==0)||(strcmp(pf,  "wymsf9.p")==0)||
	    (strcmp(pf, "wymtf7t.p")==0)||
	    (strcmp(pf,  "wymsf6.p")==0)||(strcmp(pf, "wymtf9p.p")==0)||
	    (strcmp(pf,    "wydt.p")==0)||
	    (strcmp(pf, "wymsf11.p")==0)||(strcmp(pf,  "wyamc5.p")==0)||
	    (strcmp(pf,  "wyamc7.p")==0)||(strcmp(pf,"wyamcsf5.p")==0)||
	    (strcmp(pf,   "wyamc.p")==0)||(strcmp(pf,"wyamcpha.p")==0)||
	    (strcmp(pf, "wyphasm.p")==0)||(strcmp(pf,"wyphasm8.p")==0)||
	    (strcmp(pf,"wyphasms.p")==0)||
	    (strcmp(pf,   "wyph3.p")==0)||
	    (strcmp(pf,  "wymatm.p")==0)||(strcmp(pf, "wyphas8.p")==0)){
    lat = 16;
  }else if ((strcmp(pf,    "x4st.p")==0)||(strcmp(pf,   "gpori.p")==0)||
	    (strcmp(pf,  "x4stp.p")==0)||(strcmp(pf,  "x4stp2.p")==0)||
	    (strcmp(pf,  "x4std.p")==0)){
    lat = 14;
  }else if ((strcmp(pf, "wncon100.p")==0)){
    lat = 16;
  }else if (strcmp(pf,"wnmat050.p")==0){
    lat = 30;
  }else if (strcmp(pf,"wymat160.p")==0){
    lat = 10;
  }else{
    lat = my_rint(1.5*dt + 1.0);
    printf("  *** Warning GET_ACTUAL_LATENCY: unknown pfile, using %d\n",lat);
    printf("  *** according to the rule 1.5*dt + 1.0 (dt=%.2f)\n",dt);
  }
  
  return lat;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_LATENCY_SCREEN_POS                         */
/*                                                                           */
/*  Take into account the position of the stimulus on the screen to          */
/*  compute the time between the sync pulse indicating stimulus onset and    */
/*  the actual change in luminance on the video screen.                      */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Value returned is in units of milliseconds.                            */
/*                                                                           */
/*****************************************************************************/
int get_latency_screen_pos(pf,framerate,ypos)
     char pf[];
     int framerate,ypos;
{
  int lat,delta,w,h;
  float fh;

  /*** WYETH - NEVER USED ***/
  /*** WYETH - NEVER USED ***/

  lat = get_actual_latency(pf,1000.0/(float)framerate);
  get_actual_screen_dimensions(pf,framerate,0,&w,&h);

  fh = (float)ypos/(float)h;
  delta = my_rint(fh * (1000.0/(float)framerate - 2.0));

  printf("DELTA for scrpos = %d  (ypos = %d)\n",delta,ypos);

  return (lat + delta);
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_SYNC_CODE_INT                             */
/*                                                                           */
/*  Given a list of syncs, extract an integer value from the end of the      */
/*  list.                                                                    */
/*                                                                           */
/*****************************************************************************/
int get_sync_code_int(data,n,nbit)
     int *data,n,nbit;
{
  int i,j,k;
  int tmin,flag,diff,nzer,scode;
  float dt,fdiff;

  if (n < 2)
    exit_error("GET_SYNC_CODE_INT","Too few syncs");

  tmin = data[1]-data[0];
  for(i=2;i<n;i++){
    diff = data[i]-data[i-1];
    if (diff < tmin)
      tmin = diff;
  }
  /*printf("tmin = %d\n",tmin);*/
  if ((tmin >= 9)&&(tmin <= 12)){
    printf("    Assuming 100 Hz framerate.\n");
    dt = get_wn_actual_dt(100,0);
  }else if ((tmin >= 5)&&(tmin <= 8)){
    printf("    Assuming 160 Hz framerate.\n");
    dt = get_wn_actual_dt(160,0);
  }else{
    printf("    tmin = %d\n",tmin);
    exit_error("GET_SYNC_CODE_INT","tmin not close to 10");
    dt = 0.0; /* to avoid -Wall compiler warning */
  }

  flag = 0;
  scode = 0;
  k = 0;
  printf("      ");
  for(i=1;i<n;i++){
    fdiff = (float)(data[i]-data[i-1]);
    if (flag == 0){
      /*printf("fdiff = %.2f  dt=%.2f  fabs=%.2f\n",
	fdiff,dt,(float)fabs((double)(fdiff-dt)));*/
      if (fabs(fdiff - dt) < (dt*1.3)) /* This is two syncs 1 frame apart */
	flag = 1;
    }else{
      /*printf("fdiff = %.2f\n",fdiff);*/
      nzer = my_rint(fdiff/dt) - 1;
      for(j=0;j<nzer;j++){
	printf("0");
	scode = scode << 1;
	k += 1;
      }
      scode = scode << 1;
      scode += 1;
      k += 1;
      printf("1");
    }
  }
  if (k < nbit)
    scode = scode << (nbit-k);
  for(i=0;i<nbit-k;i++)
    printf("0");
  printf("  %d\n",scode);

  return scode;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_WY2PAT_STIM_FOR_PARAMS                        */
/*                                                                           */
/*  This will get two stimulus sequences.                                    */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - "fdt" is the exact time between video frames in milliseconds.          */
/*  - This fills in with 0's between 'dt' but not when 'rpt' is used for     */
/*    stimulus dwell (thus, different from _wn_ stimuli.                     */
/*                                                                           */
/*****************************************************************************/
void get_wy2pat_stim_for_params(framerate,jumpseed,period,
				dist1,mu1,sig1,minj1,maxj1,dt1,rpt1,
				dist2,mu2,sig2,minj2,maxj2,dt2,rpt2,
				fdt,rxdata,ry1data,ry2data,rn)
     int framerate,jumpseed,period,dist1,mu1,sig1,minj1,maxj1,dt1,rpt1;
     int dist2,mu2,sig2,minj2,maxj2,dt2,rpt2;
     float fdt,**rxdata,**ry1data,**ry2data;
     int *rn;
{
  int i;
  int n,dtcnt1,dtcnt2,rptcnt1,rptcnt2;
  int ms_r1,ms_r2,ms_b1,ms_b2,ms_q,ms_t1,ms_t2,ms_ord;
  float fjump,jump1,jump2,*xdata,*y1data,*y2data,jrange1,jrange2;

  jrange1 = (float)(maxj1 - minj1 + 1);
  jrange2 = (float)(maxj2 - minj2 + 1);

  if ((dist1 == 7) || (dist2 == 7))
    dist2 = dist1;


  ms_r1 = ms_r2 = ms_q = ms_t1 = ms_t2 = 0; /* to avoid -Wall warning */
  if (dist1 == 7){ /*** m-sequence. ***/
    ms_ord = mu1; /* Order of the m-sequence. */
    n = 1;
    for(i=0;i<ms_ord;i++)
      n *= 2;
    ms_r1 = ms_r2 = 1;
    ms_q = 1 << (ms_ord-1); /* Put a 1 in the "ms_ord"th position */
    ms_t1 = jumpseed;
    ms_t2 = get_value_opposite_tap_value(ms_t1,ms_ord);

    /*** WYETH - REMOVE THIS ONCE TESTED. ***/
    /**printf(" TAPVAL:  %d  should be opposite of %d (ord %d)\n",ms_t1,ms_t2,
	   ms_ord);**/
  }else{
    n = framerate * period; /* Number of stimulus frames. */
    if (jumpseed > 0) /*** Initialize random number generator. ***/
      jumpseed = -jumpseed;
    noise_util_ran2(&jumpseed);
  }

  if (dt1 < 0) /* values > 0 signify between-frame blanks. */
    dt1 = -dt1;
  if (dt2 < 0)
    dt2 = -dt2;

  xdata = (float *)myalloc(n*sizeof(float));
  y1data = (float *)myalloc(n*sizeof(float));
  y2data = (float *)myalloc(n*sizeof(float));

  jump1 = jump2 = -1.0; /* To avoid -Wall compiler warning */
  dtcnt1 = dt1;
  dtcnt2 = dt2;
  rptcnt1 = rpt1;
  rptcnt2 = rpt2;
  for(i=0;i<n;i++){
    xdata[i] = (float)((double)i*fdt);
    if (dtcnt1 == dt1){
      dtcnt1 = 0;
      if (rptcnt1 == rpt1){
	rptcnt1 = 0;
	if (dist1 == 0) /*** Uniform ***/
	  jump1 = (float)((int)(noise_util_ran2(&jumpseed)*jrange1)+minj1);
	else if (dist1 == 1){ /*** Gaussian ***/
	  fjump = sig1 * noise_util_gasdev(&jumpseed) + mu1;
	  if (fjump > 0.0)
	    jump1 = (float)((int)(0.5+fjump));
	  else
	    jump1 = (float)((int)(-0.5+fjump));
	}else if (dist1 == 2){ /*** Binary ***/
	  if (noise_util_ran2(&jumpseed) > 0.5)
	    jump1 = (float)minj1;
	  else
	    jump1 = (float)maxj1;
	}else if (dist1 == 3){
	  fjump = noise_util_ran2(&jumpseed);
	  if (fjump < 0.33333)
	    jump1 = (float)minj1;
	  else if (fjump < 0.66667)
	    jump1 = 0.0;
	  else
	    jump1 = (float)maxj1;
	}else if (dist1 == 4){ /* Modified Ternary distribution. */
	  fjump = noise_util_ran2(&jumpseed);
	  if (fjump < 0.25)
	    jump1 = (float)minj1;
	  else if (fjump < 0.75)
	    jump1 = 0.0;
	  else
	    jump1 = (float)maxj1;
	}else if (dist1 == 5){ /* 5-valued distribution. */
	  fjump = noise_util_ran2(&jumpseed);
	  if (fjump < 0.20)
	    jump1 = (float)minj1;
	  else if (fjump < 0.40)
	    jump1 = (float)(minj1/2);
	  else if (fjump < 0.60)
	    jump1 = 0;
	  else if (fjump < 0.80)
	    jump1 = (float)(maxj1/2);
	  else
	    jump1 = (float)maxj1;
	}else if (dist1 == 6){ /* 5-valued distribution. */
	  fjump = noise_util_ran2(&jumpseed);
	  if (fjump < 0.20)
	    jump1 = (float)minj1;
	  else if (fjump < 0.40)
	    jump1 = (float)(minj1/4);
	  else if (fjump < 0.60)
	    jump1 = 0;
	  else if (fjump < 0.80)
	    jump1 = (float)(maxj1/4);
	  else
	    jump1 = (float)maxj1;
	}else if (dist1 == 7){ /* Binary m-sequence */
	  ms_b1 = ms_r1 & 1;
	  if ((ms_r1 & ms_q) != ms_q) /* If n-1 bit is not set */
	    ms_r1 = ms_r1 << 1;
	  else
	    ms_r1 = (ms_r1 << 1) ^ ms_t1;
	  if (ms_b1 == 0)
	    jump1 = (float)minj1;
	  else
	    jump1 = (float)maxj1;
	}else{
	  exit_error("GET_WY2PAT_STIM_FOR_PARAMS","Unknown distribution");
	}
      }
      y1data[i] = jump1;
      rptcnt1 += 1;
    }else
      y1data[i] = 0.0; /* No movement */

    if (dtcnt2 == dt2){
      dtcnt2 = 0;
      if (rptcnt2 == rpt2){
	rptcnt2 = 0;
	if (dist2 == 0) /*** Uniform ***/
	  jump2 = (float)((int)(noise_util_ran2(&jumpseed)*jrange2)+minj2);
	else if (dist2 == 1){ /*** Gaussian ***/
	  fjump = sig2 * noise_util_gasdev(&jumpseed) + mu2;
	  if (fjump > 0.0)
	    jump2 = (float)((int)(0.5+fjump));
	  else
	    jump2 = (float)((int)(-0.5+fjump));
	}else if (dist2 == 2){ /*** Binary ***/
	  if (noise_util_ran2(&jumpseed) > 0.5)
	    jump2 = (float)minj2;
	  else
	    jump2 = (float)maxj2;
	}else if (dist2 == 3){
	  fjump = noise_util_ran2(&jumpseed);
	  if (fjump < 0.33333)
	    jump2 = (float)minj2;
	  else if (fjump < 0.66667)
	    jump2 = 0.0;
	  else
	    jump2 = (float)maxj2;
	}else if (dist2 == 4){ /* Modified Ternary distribution. */
	  fjump = noise_util_ran2(&jumpseed);
	  if (fjump < 0.25)
	    jump2 = (float)minj2;
	  else if (fjump < 0.75)
	    jump2 = 0.0;
	  else
	    jump2 = (float)maxj2;
	}else if (dist2 == 5){ /* 5-valued distribution. */
	  fjump = noise_util_ran2(&jumpseed);
	  if (fjump < 0.20)
	    jump2 = (float)minj2;
	  else if (fjump < 0.40)
	    jump2 = (float)(minj2/2);
	  else if (fjump < 0.60)
	    jump2 = 0;
	  else if (fjump < 0.80)
	    jump2 = (float)(maxj2/2);
	  else
	    jump2 = (float)maxj2;
	}else if (dist2 == 6){ /* 5-valued distribution. */
	  fjump = noise_util_ran2(&jumpseed);
	  if (fjump < 0.20)
	    jump2 = (float)minj2;
	  else if (fjump < 0.40)
	    jump2 = (float)(minj2/4);
	  else if (fjump < 0.60)
	    jump2 = 0;
	  else if (fjump < 0.80)
	    jump2 = (float)(maxj2/4);
	  else
	    jump2 = (float)maxj2;
	}else if (dist2 == 7){ /* Binary m-sequence */
	  ms_b2 = ms_r2 & 1;
	  if ((ms_r2 & ms_q) != ms_q) /* If n-1 bit is not set */
	    ms_r2 = ms_r2 << 1;
	  else
	    ms_r2 = (ms_r2 << 1) ^ ms_t2;
	  if (ms_b2 == 0)
	    jump2 = (float)minj2;
	  else
	    jump2 = (float)maxj2;
	}else{
	  jump2 = 0.0;
	  exit_error("GET_WY2PAT_STIM_FOR_PARAMS","Unknown distribution");
	}
      }
      y2data[i] = jump2;
      rptcnt2 += 1;
    }else
      y2data[i] = 0.0; /* No movement */

    dtcnt1 += 1;
    dtcnt2 += 1;
  }
  *rxdata = xdata; *ry1data = y1data; *ry2data = y2data; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_WN_STIM_FOR_PARAMS                          */
/*                                                                           */
/*  NOTES                                                                    */
/*  - Fills stim with 0's in between 'dt'.                                   */
/*  **** WYETH - there is a difficulty of doing PTA on stimuli with rpt>1.   */
/*               A pattern like 1111 will trigger many times at 1 frame      */
/*               intervals.                                                  */
/*  - 'dt' is used for 'xdata' values.                                       */
/*                                                                           */
/*****************************************************************************/
void get_wn_stim_for_params(framerate,jumpseed,dtframe,period,distrib,mu,sigma,
			    minjump,maxjump,rpt,dt,rxdata,rydata,rn)
     int framerate,jumpseed,dtframe,period,distrib,mu,sigma,minjump,maxjump;
     int rpt;
     float dt,**rxdata,**rydata;
     int *rn;
{
  int i;
  int n,posdtframe,count,rptcnt;
  float fjump,*xdata,*ydata,jumprange,jump;
  int ms_r,ms_b,ms_q,ms_t,ms_ord;

  /*printf("  *** GET_WN_STIM_FOR_PARAMS\n");*/

  jumprange = -1.0; /* To avoid -Wall compiler warning */
  ms_r = ms_q = ms_t = 0;

  if (distrib == 7){ /*** m-sequence. ***/
    ms_ord = mu; /* Order of the m-sequence. */
    n = 1;
    for(i=0;i<ms_ord;i++)
      n *= 2;
    ms_r = 1;
    ms_q = 1 << (ms_ord-1); /* Put a 1 in the "ms_ord"th position */
    ms_t = jumpseed;
  }else{
    jumprange = (float)(maxjump - minjump + 1);
    n = framerate * period; /* Number of stimulus frames. */
    if (jumpseed > 0) /*** Initialize random number generator. ***/
      jumpseed = -jumpseed;
    noise_util_ran2(&jumpseed);
  }

  posdtframe = dtframe;
  if (posdtframe < 0) /* dtframe values < 0 signify between-frame blanks. */
    posdtframe *= -1;

  jump = -1.0; /* Avoid -Wall warning */
  xdata = (float *)myalloc(n*sizeof(float));
  ydata = (float *)myalloc(n*sizeof(float));
  count = posdtframe;
  rptcnt = rpt;
  for(i=0;i<n;i++){
    xdata[i] = (float)((double)i*dt);
    if (count == posdtframe){
      if (rptcnt == rpt){
	rptcnt = 0;
	if (distrib == 0) /*** Uniform ***/
	  jump = (float)((int)(noise_util_ran2(&jumpseed)*jumprange) +
			     minjump);
	else if (distrib == 1){ /*** Gaussian ***/
	  fjump = sigma * noise_util_gasdev(&jumpseed) + mu;
	  if (fjump > 0.0)
	    jump = (float)((int)(0.5+fjump));
	  else
	    jump = (float)((int)(-0.5+fjump));
	}else if (distrib == 2){ /*** Binary ***/
	  if (noise_util_ran2(&jumpseed) > 0.5)
	    jump = (float)minjump;
	  else
	    jump = (float)maxjump;
	}else if (distrib == 3){
	  fjump = noise_util_ran2(&jumpseed);
	  if (fjump < 0.33333)
	    jump = (float)minjump;
	  else if (fjump < 0.66667)
	    jump = 0.0;
	  else
	    jump = (float)maxjump;
	}else if (distrib == 4){ /* Modified Ternary distribution. */
	  fjump = noise_util_ran2(&jumpseed);
	  if (fjump < 0.25)
	    jump = (float)minjump;
	  else if (fjump < 0.75)
	    jump = 0.0;
	  else
	    jump = (float)maxjump;
	}else if (distrib == 5){ /* 5-valued distribution. */
	  fjump = noise_util_ran2(&jumpseed);
	  if (fjump < 0.20)
	    jump = (float)minjump;
	  else if (fjump < 0.40)
	    jump = (float)(minjump/2);
	  else if (fjump < 0.60)
	    jump = 0;
	  else if (fjump < 0.80)
	    jump = (float)(maxjump/2);
	  else
	    jump = (float)maxjump;
	}else if (distrib == 6){ /* 5-valued distribution. */
	  fjump = noise_util_ran2(&jumpseed);
	  if (fjump < 0.20)
	    jump = (float)minjump;
	  else if (fjump < 0.40)
	    jump = (float)(minjump/4);
	  else if (fjump < 0.60)
	    jump = 0;
	  else if (fjump < 0.80)
	    jump = (float)(maxjump/4);
	  else
	    jump = (float)maxjump;
	}else if (distrib == 7){ /* Binary m-sequence */
	  ms_b = ms_r & 1;
	  if ((ms_r & ms_q) != ms_q) /* If n-1 bit is not set */
	    ms_r = ms_r << 1;
	  else
	    ms_r = (ms_r << 1) ^ ms_t;
	  if (ms_b == 0)
	    jump = (float)minjump;
	  else
	    jump = (float)maxjump;
	}else
	  exit_error("GET_WN_STIM_FOR_PARAMS","Unknown distribution");
      } /* wyeth - new */
      ydata[i] = jump; /* wyeth - new */
      rptcnt += 1; /* wyeth - new */
      
      count = 0;
    }else
      ydata[i] = 0.0; /* No movement */

    /*printf("ydata[i] = %.2f\n",ydata[i]);*/

    count += 1;
  }
  *rxdata = xdata; *rydata = ydata; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                        GET_AMP_MOD_STIM_FOR_PARAMS                        */
/*                                                                           */
/*  Get stimulus for amplitude modulation of contrast (e.g. wyamc7, wyamc5)  */
/*  where stimulus may be reflected at endpoints of range.                   */
/*                                                                           */
/*  walkflag - determined from stimulus type.                                */
/*  amp0 - pb->pat2amp                                                       */
/*  agflag - pb->sigma2                                                      */
/*                                                                           */
/*****************************************************************************/
void get_amp_mod_stim_for_params(framerate,jumpseed,dtframe,period,distrib,mu,
				 sigma,minjump,maxjump,walkflag,steptype,amp0,
				 dt,rxdata,rydata,rn,pos_flag)
     int framerate,jumpseed,dtframe,period,distrib,mu,sigma,minjump,maxjump;
     int walkflag,steptype,amp0;
     float dt,**rxdata,**rydata;
     int *rn,pos_flag;
{
  int i;
  int n,posdtframe,count,curramp,maxampv,oldcurramp;
  float fjump,*xdata,*ydata,jumprange;
  int ms_r,ms_b,ms_q,ms_t,ms_ord,NTTAB;

  /*printf("    pos_flag = %d\n",pos_flag);*/

  jumprange = -1.0; /* Avoid -Wall warning */
  ms_r = ms_q = ms_t = 0;

  if (distrib == 7){ /*** m-sequence. ***/
    ms_ord = mu; /* Order of the m-sequence. */
    n = 1;
    for(i=0;i<ms_ord;i++)
      n *= 2;
    ms_r = 1;
    ms_q = 1 << (ms_ord-1); /* Put a 1 in the "ms_ord"th position */
    ms_t = jumpseed;
  }else{
    jumprange = (float)(maxjump - minjump + 1);
    n = framerate * period; /* Number of stimulus frames. */
    if (jumpseed > 0) /*** Initialize random number generator. ***/
      jumpseed = -jumpseed;
    noise_util_ran2(&jumpseed);
  }

  NTTAB = 129;
  if (steptype == 0) /*** From wy2patx.c:  (pb->sigma2 == 0) ***/
    maxampv = NTTAB-1; /* Arithmetic progression */
  else
    maxampv = 32; /* Geometric progression */

  curramp = amp0;  /*** In wy2patx.c:  pb->pat2amp */
  if (curramp < -maxampv)
    curramp = -maxampv;
  else if (curramp > maxampv)
    curramp = maxampv;

  posdtframe = dtframe;
  if (posdtframe < 0) /* dtframe values < 0 signify between-frame blanks. */
    posdtframe *= -1;

  xdata = (float *)myalloc(n*sizeof(float));
  ydata = (float *)myalloc(n*sizeof(float));
  count = posdtframe;
  for(i=0;i<n;i++){
    xdata[i] = (float)((double)i*dt);
    if (count == posdtframe){
      if (distrib == 0) /*** Uniform ***/
	ydata[i] = (float)((int)(noise_util_ran2(&jumpseed)*jumprange) +
			   minjump);
      else if (distrib == 1){ /*** Gaussian ***/
	fjump = sigma * noise_util_gasdev(&jumpseed) + mu;
	if (fjump > 0.0)
	  ydata[i] = (float)((int)(0.5+fjump));
	else
	  ydata[i] = (float)((int)(-0.5+fjump));
      }else if (distrib == 2){ /*** Binary ***/
	if (noise_util_ran2(&jumpseed) > 0.5)
	  ydata[i] = (float)minjump;
	else
	  ydata[i] = (float)maxjump;
      }else if (distrib == 3){
	fjump = noise_util_ran2(&jumpseed);
	if (fjump < 0.33333)
	  ydata[i] = (float)minjump;
	else if (fjump < 0.66667)
	  ydata[i] = 0.0;
	else
	  ydata[i] = (float)maxjump;
      }else if (distrib == 4){ /* Modified Ternary distribution. */
	fjump = noise_util_ran2(&jumpseed);
	if (fjump < 0.25)
	  ydata[i] = (float)minjump;
	else if (fjump < 0.75)
	  ydata[i] = 0.0;
	else
	  ydata[i] = (float)maxjump;
      }else if (distrib == 5){ /* 5-valued distribution. */
	fjump = noise_util_ran2(&jumpseed);
	if (fjump < 0.20)
	  ydata[i] = (float)minjump;
	else if (fjump < 0.40)
	  ydata[i] = (float)(minjump/2);
	else if (fjump < 0.60)
	  ydata[i] = 0;
	else if (fjump < 0.80)
	  ydata[i] = (float)(maxjump/2);
	else
	  ydata[i] = (float)maxjump;
      }else if (distrib == 6){ /* 5-valued distribution. */
	fjump = noise_util_ran2(&jumpseed);
	if (fjump < 0.20)
	  ydata[i] = (float)minjump;
	else if (fjump < 0.40)
	  ydata[i] = (float)(minjump/4);
	else if (fjump < 0.60)
	  ydata[i] = 0;
	else if (fjump < 0.80)
	  ydata[i] = (float)(maxjump/4);
	else
	  ydata[i] = (float)maxjump;
      }else if (distrib == 7){ /* Binary m-sequence */
	ms_b = ms_r & 1;
	if ((ms_r & ms_q) != ms_q) /* If n-1 bit is not set */
	  ms_r = ms_r << 1;
	else
	  ms_r = (ms_r << 1) ^ ms_t;
	if (ms_b == 0)
	  ydata[i] = (float)minjump;
	else
	  ydata[i] = (float)maxjump;
      }else
	exit_error("GET_AMP_MOD_STIM_FOR_PARAMS","Unknown distribution");
      count = 0;
    }else
      ydata[i] = 0.0;

    oldcurramp = curramp;

    if (walkflag==1)
      curramp += ydata[i];
    else
      curramp = ydata[i];

    /*printf("curramp %d ",curramp);*/
    while((curramp < -maxampv)||(curramp > maxampv)){
      if (curramp > maxampv) /*** Keep curramp in range 0..NTTAB-1 ***/
	curramp = maxampv-(curramp-maxampv); /* Reflect the jump. */
      if (curramp < -maxampv) /*** Keep curramp in range 0..NTTAB-1 ***/
	curramp = -maxampv-(curramp+maxampv); /* Reflect the jump. */
      /*printf(" wrapped to %d   ydata: %.1f to",curramp,ydata[i]);*/

      ydata[i] = (float)(curramp - oldcurramp);
      
      /*printf(" %.1f\n",ydata[i]);*/
    }
    /*printf("\n");*/
    
    /*printf("ydata[i] = %.2f %d\n",ydata[i],curramp);*/
    if (pos_flag == 1)
      ydata[i] = (float)curramp;
    
    count += 1;
  }

  /*printf("\nNO CORRECTIONS\n");*/

  *rxdata = xdata; *rydata = ydata; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                       GET_4STATE_REL_STIM_FOR_PARAMS                      */
/*                                                                           */
/*****************************************************************************/
void get_4state_rel_stim_for_params(framerate,jumpseed,period,tf1,tf2,deltamin,
				    deltamax,ncycles,rsdt,rstim,rn)
     int framerate,jumpseed,period,tf1,tf2,deltamin,deltamax,ncycles;
     float *rsdt;
     int **rstim,*rn;
{
  int i,j;
  int n,ns,count,dwell,dwell2,nfmin,range,*stim;
  float fdt;

  if (jumpseed > 0) /*** Initialize random number generator. ***/
    jumpseed = -jumpseed;
  noise_util_ran2(&jumpseed);

  if (tf1 == 0) /* wy2flatr.p */
    dwell = tf2; /* For wy2flatr.p, get dwell from tf2. */
  else /* wy4strel.p, wy4strev.p */
    dwell = 1024/tf1 * ncycles; /* Number of frames per step. */
  dwell2 = 2*dwell;
  nfmin = period * framerate; /* Minimum duration requested. */
  n = dwell*4; /* Nframes to go once through all 4 steps. */
  while (n < nfmin)
    n += dwell*4; /* Nframes to go once through all 4 steps. */
  range = deltamax-deltamin+1;

  ns = n/dwell2; /* Number of time-varied transitions. */
  stim = (int *)myalloc(ns*sizeof(int));

  count = dwell2;
  j = 0;
  for(i=0;i<n;i++){
    if (count == dwell2){ /* Pick new dt, set appropriate phase */
      stim[j] = (int)(noise_util_ran2(&jumpseed)*range) + deltamin;
      j += 1;
      count = 0;
    }
    count += 1;
  }
  if (j != ns)
    exit_error("GET_4STATE_REL_STIM_FOR_PARAMS","counting error");

  fdt = get_wn_actual_dt(framerate,0);
  *rsdt = (float)dwell2*fdt;
  *rstim = stim; *rn = ns;
}
/**************************************-**************************************/
/*                                                                           */
/*                      GET_4STATE_SIMPLE_STIM_FOR_PARAMS                    */
/*                                                                           */
/*****************************************************************************/
void get_4state_simple_stim_for_params(pfile,framerate,seed,period,tf1,tf2,
				       nphase,ncycles,rt0,rsdt,rstim,rn)
     char pfile[];
     int framerate,seed,period,tf1,tf2,nphase,ncycles;
     float *rt0,*rsdt;
     int **rstim,*rn;
{
  int i;
  int n,ns,dwell,*stim,t0,tpara,torth,dframe,ttot,ph0,ranph;
  float fdt;

  if (seed > 0) /*** Initialize random number generator. ***/
    seed = -seed;
  noise_util_ran2(&seed);

  /*** WYETH WYETH - what happens for static???? ***/
  if (tf1 == 0) /* WAS for wy2flatr.p */
    dwell = tf2; /* WAS For wy2flatr.p, get dwell from tf2. */
  else /* wy4strel.p wy4strev.p */
    dwell = 1024/tf1; /* Number of frames per step. */

  dframe = dwell/nphase;
  if (dframe < 1)
    dframe = 1;
  tpara = dwell*ncycles;
  torth = dwell*ncycles + dframe;

  /*** Compute 'ns', the number of "stimulus" transitions during the trial ***/
  n = period * framerate;
  ttot = tpara+torth; /* nframes between 'stimulus' transition */
  ns = (n-dwell)/ttot; /* Number of 'stimuli' during the trial */
  stim = (int *)myalloc(ns*sizeof(int));

  printf("dwell=%d\n",dwell);
  printf("ncycles=%d\n",ncycles);
  printf("n=%d\n",n);
  printf("ttot=%d\n",ttot);
  printf("framerate=%d\n",framerate);
  printf("NUMBER OF TRANSITIONS PER TRIAL = %d\n",ns);

  /*** Find the starting phase ***/
  ranph = (int)(noise_util_ran2(&seed)*nphase);
  t0 = dframe * ranph; /* Number of frames before '0' */
  /*count = -t0;*/
  ph0 = 1 + ranph; /* '1' accounts for dframe during torth */

  printf("INITIAL PHASE = %d\n",ph0);

  /*** On each transition, phase moves forward by one quanta. ***/
  for(i=0;i<ns;i++)
    stim[i] = (ph0+i)%nphase;

  if ((strcmp(pfile,"x4st.p")==0)||(strcmp(pfile,"x4stp.p")==0)||
      (strcmp(pfile,"x4stp2.p")==0)||(strcmp(pfile,"x4std.p")==0)){
    printf("  *** WARNING:  Est. SGI dt OLD WAY, PROBABLY wrong.\n");
    fdt = get_sgi_actual_dt(framerate,0);
  }else
    fdt = get_wn_actual_dt(framerate,0);

  *rt0 = (float)(t0+torth)*fdt; /* Actual time before first transition */
  *rsdt = (float)ttot*fdt;      /* Actual time between stimuli */
  *rstim = stim; *rn = ns;      /* Stimulus phase numbers */
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_WY2PA_ANTIDUR_LIST                        */
/*                                                                           */
/*  Return the list of antipreferred durations.                              */
/*                                                                           */
/*****************************************************************************/
void get_wy2pa_antidur_list(antidur,prefdur,ralist,rn)
     int antidur,prefdur,**ralist,*rn;
{
  int n,*alist;

  if (antidur == 0){
    n = 8;
    alist = (int *)myalloc(n*sizeof(int));
    alist[0] = 0;
    alist[1] = 1;
    alist[2] = 2;
    alist[3] = 3;
    alist[4] = 4;
    alist[5] = 6;
    alist[6] = 8;
    alist[7] = 10;
  }else{
    if (antidur < 0){
      n = 1;
      alist = (int *)myalloc(n*sizeof(int));
      alist[0] = prefdur;
    }else{
      n = 2;
      alist = (int *)myalloc(n*sizeof(int));
      alist[0] = 0;
      alist[1] = antidur;
    }
  }
  *ralist = alist; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_WY2PA_STIM_FOR_PARAMS                       */
/*                                                                           */
/*  fdt - actual float time between frames                                   */
/*  rdt - returned float time between pref movements.                        */
/*  rtanti - anti duration times                                             */
/*  rpstate - which stim does which movement                                 */
/*    00 - patch 0 anti, 0 pref                                              */
/*    01 - patch 0 anti, 1 pref                                              */
/*    10 - patch 1 anti, 0 pref                                              */
/*    11 - patch 1 anti, 1 pref                                              */
/*                                                                           */
/*****************************************************************************/
void get_wy2pa_stim_for_params(pfile,seed,framerate,dur,antidur,prefdur,
			       toverlap,fdt,rdt,rpstate,rtanti,rn)
     char pfile[];
     int seed,framerate,dur,antidur,prefdur,toverlap;
     float fdt,*rdt;
     int **rpstate,**rtanti,*rn;
{
  int i,k;
  int n,ns,*stim,p0,p1,a0,a1,count,full,antimax,tanti,apat,ppat,jump1,jump2;
  int *pstate,printflag,tstatic,pflag;

  printflag = 0;

  a0 = tanti = apat = ppat = -1; /* Avoid -Wall warning */

  tstatic = 4;
  pflag = 0;
  if (antidur < 0){ /* Use this for simultaneous pref stimuli */
    pflag = 1;
    tstatic = -antidur;
    antidur = prefdur;
    antimax = 0;
  }else if (antidur == 0) /* Use this to implement other anti sequences */
    antimax = 10;
  else
    antimax = antidur;

  if (pflag){
    p0 = tstatic;
    p1 = p0+prefdur-1;
    a1 = p1;
    full = tstatic + prefdur;
    toverlap = 0;
  }else{
    p0 = antimax + 4 - toverlap;    /* First frame of pref */
    p1 = p0+prefdur-1;              /* Last frame of pref */
    full = p0+prefdur + 2;
    a1 = p0-1 + toverlap;           /* Last frame of antipref */
  }

  if (seed > 0) /*** Initialize random number generator. ***/
    seed = -seed;

  if (printflag)
    printf("SEED = %d\n",seed);

  n = dur * framerate;

  ns = n/full;

  stim = (int *)myalloc(ns*sizeof(int));
  pstate = (int *)myalloc(ns*sizeof(int));

  k = -1; /* Count number of stimuli */
  count = full;
  for(i=0;i<n;i++){
    if (count == full){
      if (printflag)
	printf("------------\n");
      count = 0;
      if (antidur == 0){
	tanti = (int)(noise_util_ran2(&seed)*8); /* dur of antipulse */
	if (tanti == 5)
	  tanti = 6;
	else if (tanti == 6)
	  tanti = 8;
	else if (tanti == 7)
	  tanti = 10;
      }else{
	if (noise_util_ran2(&seed) < 0.5)
	  tanti = 0;
	else
	  tanti = antidur;
      }
      a0 = a1 - tanti + 1;
      if (noise_util_ran2(&seed) < 0.5) /* 'apat' patch does anti motion */
	apat = 0;
      else
	apat = 1;
      if (noise_util_ran2(&seed) < 0.5) /* 'ppat' patch does pref motion */
	ppat = 0;
      else
	ppat = 1;

      k += 1;
      if (pflag == 1){
	if ((apat == 0)&&(ppat == 0))
	  pstate[k] = 0;
	else if ((apat == 1)&&(ppat == 1))
	  pstate[k] = 1;
	else if ((apat == 0)&&(ppat == 1)){
	  if (tanti == 0)
	    pstate[k] = 1;
	  else
	    pstate[k] = 2;
	}else{ /* 1 0 */
	  if (tanti == 0)
	    pstate[k] = 0;
	  else
	    pstate[k] = 2;
	}
	stim[k] = antidur;  /* Always antidur */
      }else{
	stim[k] = tanti;
	pstate[k] = apat*2 + ppat;
      }
      /*printf("%d %d\n",stim[k],pstate[k]);*/
    }

    jump1 = jump2 = 0;
    if ((count >= a0)&&(count <= a1)){  /* Show anti */
      if (apat)	jump1 = -1; /* replaced minj/maxj with -1/+1 */
      else	jump2 = -1;
    }
    if ((count >= p0)&&(count <= p1)){  /* Show pref (may overwrite anti) */
      if (ppat)	jump1 = 1;
      else	jump2 = 1;
    }

    if (printflag==1){
      printf("%d %d  %d\n",apat,ppat,tanti);
    }else if (printflag==2){
      printf("%2d  ",count);
      if (jump1 == 1)  printf("p");
      else if (jump1 == -1) printf("a");
      else printf(".");
      
      if (jump2 == 1)  printf("p\n");
      else if (jump2 == -1) printf("a\n");
      else printf(".\n");
    }

    count += 1;
  }
  if (count == full)
    k += 1;

  printf("  k=%d  ns = %d  Both are no. of stimuli.\n",k,ns);

  *rdt = (float)full*fdt;    /* Actual time between stimuli */
  *rpstate = pstate; *rtanti = stim; *rn = ns;   /* Stimulus values */
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_WYTFV_STIM_FOR_PARAMS                       */
/*                                                                           */
/*  fdt - actual float time between frames                                   */
/*  rdt - returned float time between pref movements.                        */
/*                                                                           */
/*****************************************************************************/
void get_wytfv_stim_for_params(pfile,seed,framerate,dur,zflag,t1,t2,t3,ntf,
			       tfmax,fdt,rdt,rstima,rstimp,rn)
     char pfile[];
     int seed,framerate,dur,zflag,t1,t2,t3,ntf,tfmax;
     float fdt,*rdt;
     int **rstimp,**rstima,*rn;
{
  int i,k;
  int ns,*stima,*stimp;
  int cnt,nper,ncyc,nframes,tfmin,tf1,tf2;
  int *ranlist1,*ranlist2,*rilist;
    
  if (seed > 0)
    seed = -seed;

  nper = t1+t2+t3; /* Total period in frames */
  ncyc = (dur * framerate)/nper; /* No. of cycles in duration */
  nframes = ncyc * nper; /* Total number of frames to show */

  ranlist1 = (int *)myalloc(ncyc*sizeof(int));
  ranlist2 = (int *)myalloc(ncyc*sizeof(int));
  rilist = (int *)myalloc(ncyc*sizeof(int));

  ns = ncyc;
  stima = (int *)myalloc(ns*sizeof(int));
  stimp = (int *)myalloc(ns*sizeof(int));

  /*** Put the TF values in 'tflist' to be chosen from later. ***/
  if (zflag == 1)
    ntf += 1;

  tfmin = tfmax;
  for(i=1;i<ntf;i++)
    tfmin = tfmin/2;

  tf1 = tfmax;
  tf2 = tfmax;
  if (ntf*ntf > ncyc) /*** ERROR, stimulus duration too short ***/
    tf1 = tf2 = tfmax = tfmin = 0; /* Show all static stimulus as error */
  for(i=0;i<ncyc;i++){
    rilist[i] = i;
    if ((zflag == 1)&&(tf1 == tfmin))
      ranlist1[i] = 0;
    else
      ranlist1[i] = tf1;
    
    if ((zflag == 1)&&(tf2 == tfmin))
      ranlist2[i] = 0;
    else
      ranlist2[i] = tf2;
    
    if (tf2 == tfmin){
      tf2 = tfmax;
      if (tf1 == tfmin)
	tf1 = tfmax;
      else
	tf1 = tf1/2;
    }else
      tf2 = tf2/2;
  }
    
  shuffle_iarray(rilist,ncyc,1,seed);

  cnt = nper;
  k = 0; /* Number of cycles done */
  for(i=0;i<nframes;i++){
    if (cnt == nper){ /* Pick new TFs (jumps) */
      stima[k] = -ranlist1[rilist[k]]; /* Anti */
      stimp[k] =  ranlist2[rilist[k]]; /* Pref */
      k += 1;
      cnt = 0;
    }
    cnt += 1;
  }
  myfree(ranlist1); myfree(ranlist2); myfree(rilist);

  *rdt = (float)nper*fdt;    /* Actual time between stimuli */
  *rstima = stima; *rstimp = stimp; *rn = ns;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_WN_STIM_FOR_TRIAL                          */
/*                                                                           */
/*  Return the stimulus sequence for the "k"the trial.  This routine covers  */
/*  the following p-files:                                                   */
/*    wnall100.p                                                             */
/*    wncon100.p                                                             */
/*    wnmat050.p wnmat053.p wnmat080.p wnmat100.p                            */
/*    wnmdt100.p                                                             */
/*    wnone053.p wnone080.p wnone100.p                                       */
/*    wnposdt.p                                                              */
/*    wnsiz100.p                                                             */
/*                                                                           */
/*  Notes:                                                                   */
/*  - xdata values are returned in milliseconds.                             */
/*                                                                           */
/*****************************************************************************/
void get_wn_stim_for_trial(nd,k,awake_flag,pfile,rxdata,rydata,rn)
     struct ndata_struct *nd;
     int k,awake_flag;
     char *pfile;
     float **rxdata,**rydata;
     int *rn;
{
  int i,j;
  int framerate,jumpseed,dtframe,period,distrib,mu,sigma,minjump,maxjump;
  int wavelen,flag,wvo2,count,rpt;
  float dt;

  get_wn_stim_params(nd,k,awake_flag,&framerate,&jumpseed,&dtframe,&period,
		     &distrib,&mu,&sigma,&minjump,&maxjump);
  if (strcmp(pfile,"MODEL.p")==0)
    dt = 1000.0/framerate;
  else
    dt = get_wn_actual_dt(framerate,awake_flag);
  rpt = 1; /* WYETH */
  get_wn_stim_for_params(framerate,jumpseed,dtframe,period,distrib,mu,sigma,
			 minjump,maxjump,rpt,dt,rxdata,rydata,rn);

  if (awake_flag){
    flag = ndata_get_vc_param_int(nd,k,"wavelen",&wavelen);
    if (!flag) exit_error("GET_WN_STIM_FOR_TRIAL","Param wavelen not found");
  }else{
    flag = ndata_get_vc_param_int(nd,k,"min_wvln",&wavelen);
    if (!flag){
      flag = ndata_get_vc_param_int(nd,k,"sper",&wavelen);
      if (!flag)
	exit_error("GET_WN_STIM_FOR_TRIAL","Neither min_wvln nor sper found");
    }
  }

  count = 0;
  /*** Watch out for spatial GWN. ***/
  if ((strcmp(pfile,"wynoise.p")!=0)&&(strcmp(pfile,"wngnv.p")!=0)){
    wvo2 = wavelen/2;
    for(i=0;i<*rn;i++){
      j = my_rint((*rydata)[i]);
      if ((j > wvo2)||(j <= -wvo2)){
	while (j > wvo2)
	  j -= wavelen;
	while (j <= -wvo2)
	  j += wavelen;
	printf("Change %d to %d (wvln=%d)\n",(int)((*rydata)[i]),j,wavelen);
	(*rydata)[i] = (float)j;
	count += 1;
      }
    }
  }
  if (count > 0)
    printf("  GET_WN_STIM_FOR_TRIAL:  %d jump values wrapped around.\n",count);
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_WY2PAT_STIM_FOR_TRIAL_01                      */
/*                                                                           */
/*  Made from:  get_wn_stim_for_trial.                                       */
/*  Return the stimulus sequence for the "k"the trial.  This routine covers  */
/*  the following p-files:                                                   */
/*    wypnori.p                                                              */
/*    wypnpha.p                                                              */
/*    wy4st.p                                                                */
/*    wy4stg.p                                                               */
/*    wy2flat.p                                                              */
/*     ...                                                                   */
/*                                                                           */
/*  Notes:                                                                   */
/*  - xdata values are returned in milliseconds.                             */
/*  - 'wyamc7' comes here, but apparently there is no adjustment for         */
/*    reflection which occurs for stimtype '22'.                             */
/*                                                                           */
/*****************************************************************************/
void get_wy2pat_stim_for_trial_01(nd,k,awake_flag,pfile,rxdata,rydata,rn)
     struct ndata_struct *nd;
     int k,awake_flag;
     char *pfile;
     float **rxdata,**rydata;
     int *rn;
{
  int i;
  int framerate,jumpseed,dtframe,period,distrib,mu,sigma,minjump,maxjump;
  int wavelen,flag,count,rpt;
  float dt,x;

  /*printf("GET_WY2PAT_STIM_FOR_TRIAL_01\n");*/

  get_wy2pat_stim_params(nd,k,awake_flag,&framerate,&jumpseed,&dtframe,
			 &period,&distrib,&mu,&sigma,&minjump,&maxjump,&rpt);

  if ((strcmp(pfile,"MODEL.p")==0)||(strcmp(pfile,"MODEL_new.p")==0))
    dt = 1000.0/framerate;
  else
    dt = get_wn_actual_dt(framerate,awake_flag);
  get_wn_stim_for_params(framerate,jumpseed,dtframe,period,distrib,mu,sigma,
			 minjump,maxjump,rpt,dt,rxdata,rydata,rn);
  flag = ndata_get_vc_param_int(nd,k,"sper",&wavelen);
  if (!flag) exit_error("GET_WY2PAT_STIM_FOR_TRIAL_01","Param sper not found");

  /***for(i=0;i<*rn;i++)
    printf("%.2f %.2f\n",(*rxdata)[i],(*rydata)[i]);***/

  /*** Adjust large jumps that wrap around, e.g. for Gaussian distrib. ***/
  /*** Jumps are given in SPER units where 1024 = 1 full cycle of jump. ***/
  count = 0;
  for(i=0;i<*rn;i++){
    x = (*rydata)[i];
    if ((x > 512.0) || (x <= -512.0)){
      count += 1;
      /*printf("x = %f\n",x);*/
      while (x > 512.0){
	printf(" ****** ADJUSTING %f\n",x);
	x -= 1024.0;
      }
      while (x <= -512.0){
	printf(" ****** ADJUSTING %f\n",x);
	x += 1024.0;
      }
      (*rydata)[i] = x;
    }
  }
  if (count > 0){
    printf("  GET_WY2PAT_STIM_FOR_TRIAL_01:  %d stimulus jump values",count);
    printf(" adjusted for wrap-around\n");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_WY2PAT_STIM_FOR_TRIAL_02                      */
/*                                                                           */
/*  Return the stimulus sequences for the "k"the trial.  This routine        */
/*  covers the following p-files:                                            */
/*    wy4stdir.p                                                             */
/*     ...                                                                   */
/*                                                                           */
/*  Notes:                                                                   */
/*  - xdata values are returned in milliseconds.                             */
/*                                                                           */
/*****************************************************************************/
void get_wy2pat_stim_for_trial_02(nd,k,awake_flag,pfile,rxdata,
				  ry1data,ry2data,rn,rfdt)
     struct ndata_struct *nd;
     int k,awake_flag;
     char *pfile;
     float **rxdata,**ry1data,**ry2data;
     int *rn;
     float *rfdt;
{
  int i;
  int count,framerate,jumpseed,period,dt1,rpt1,dist1,mu1,sig1,minj1,maxj1;
  int dt2,rpt2,dist2,mu2,sig2,minj2,maxj2;
  float x;

  get_wy2pat_stim_params_02(nd,k,awake_flag,&framerate,&jumpseed,&period,
			    &dt1,&rpt1,&dist1,&mu1,&sig1,&minj1,&maxj1,
			    &dt2,&rpt2,&dist2,&mu2,&sig2,&minj2,&maxj2);

  if ((strcmp(pfile,"MODEL.p")==0)||(strcmp(pfile,"MODEL_new.p")==0))
    *rfdt = 1000.0/framerate;
  else
    *rfdt = get_wn_actual_dt(framerate,awake_flag);

  get_wy2pat_stim_for_params(framerate,jumpseed,period,
			     dist1,mu1,sig1,minj1,maxj1,dt1,rpt1,
			     dist2,mu2,sig2,minj2,maxj2,dt2,rpt2,
			     *rfdt,rxdata,ry1data,ry2data,rn);

  /*** Adjust large jumps that wrap around, e.g. for Gaussian distrib. ***/
  /*** Jumps are given in SPER units where 1024 = 1 full cycle of jump. ***/
  count = 0;
  for(i=0;i<*rn;i++){
    x = (*ry1data)[i];
    if ((x > 512.0) || (x <= -512.0)){
      count += 1;
      while (x > 512.0)
	x -= 1024.0;
      while (x <= -512.0)
	x += 1024.0;
      (*ry1data)[i] = x;
    }
    x = (*ry2data)[i];
    if ((x > 512.0) || (x <= -512.0)){
      count += 1;
      while (x > 512.0)
	x -= 1024.0;
      while (x <= -512.0)
	x += 1024.0;
      (*ry2data)[i] = x;
    }
  }

  if (count > 0){
    printf("  GET_WY2PAT_STIM_FOR_TRIAL_02:  %d stimulus jump values",count);
    printf(" adjusted for wrap-around\n");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_WY2PAT_STIM_FOR_TRIAL_03                      */
/*                                                                           */
/*  Return the stimulus sequence for the "k"the trial.  This routine covers  */
/*  the following p-files:                                                   */
/*    wyamc5.p                                                               */
/*    wyamc7.p                                                               */
/*     ...                                                                   */
/*                                                                           */
/*  Notes:                                                                   */
/*  - xdata values are returned in milliseconds.                             */
/*                                                                           */
/*****************************************************************************/
void get_wy2pat_stim_for_trial_03(nd,k,awake_flag,pos_flag,pfile,rxdata,
				  rydata,rn)
     struct ndata_struct *nd;
     int k,awake_flag,pos_flag;
     char *pfile;
     float **rxdata,**rydata;
     int *rn;
{
  int framerate,jumpseed,dtframe,period,distrib,mu,sigma,minjump,maxjump;
  int wavelen,flag,walkflag,steptype,amp0;
  float dt;

  get_wy2pat_stim_params_03(nd,k,&framerate,&jumpseed,&dtframe,&period,
			    &distrib,&mu,&sigma,&minjump,&maxjump,&walkflag,
			    &steptype,&amp0);

  if ((strcmp(pfile,"MODEL.p")==0)||(strcmp(pfile,"MODEL_new.p")==0))
    dt = 1000.0/framerate;
  else
    dt = get_wn_actual_dt(framerate,awake_flag);

  get_amp_mod_stim_for_params(framerate,jumpseed,dtframe,period,distrib,mu,
			      sigma,minjump,maxjump,walkflag,steptype,amp0,dt,
			      rxdata,rydata,rn,pos_flag);

  flag = ndata_get_vc_param_int(nd,k,"sper",&wavelen);
  if (!flag) exit_error("GET_WY2PAT_STIM_FOR_TRIAL_03","Param sper not found");

}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_WY2PAT_STIM_FOR_TRIAL_04                      */
/*                                                                           */
/*  Return the stimulus sequence for the "k"the trial.  This routine covers  */
/*  the following p-files:                                                   */
/*    wy2flat2.p                                                             */
/*                                                                           */
/*  Notes:                                                                   */
/*  - xdata values are returned in milliseconds.                             */
/*                                                                           */
/*****************************************************************************/
void get_wy2pat_stim_for_trial_04(nd,k,awake_flag,pfile,rxdata,rydata,rn)
     struct ndata_struct *nd;
     int k,awake_flag;
     char *pfile;
     float **rxdata,**rydata;
     int *rn;
{
  int framerate,seed,period;
  int dist1,dist2,mu1,mu2,sig1,sig2,minj1,minj2,maxj1,maxj2,dt1,dt2,rpt1,rpt2;
  int pat1amp,pat2amp;
  float fdt,*y1,*y2;

  /*printf("GET_WY2PAT_STIM_FOR_TRIAL_04\n");*/

  get_wy2pat_stim_params_04(nd,k,awake_flag,&framerate,&seed,&period,
			    &dt1,&rpt1,&dist1,&pat1amp,
			    &dt2,&rpt2,&dist2,&pat2amp);

  if ((strcmp(pfile,"MODEL.p")==0)||(strcmp(pfile,"MODEL_new.p")==0))
    fdt = 1000.0/framerate;
  else
    fdt = get_wn_actual_dt(framerate,awake_flag);

  mu1 = mu2 = sig1 = sig2 = 0;
  minj1 = minj2 = -1;
  maxj1 = maxj2 = 1;
  get_wy2pat_stim_for_params(framerate,seed,period,
			     dist1,mu1,sig1,minj1,maxj1,dt1,rpt1,
			     dist2,mu2,sig2,minj2,maxj2,dt2,rpt2,
			     fdt,rxdata,&y1,&y2,rn);
  if (pat1amp == 0){
    *rydata = y2;
    myfree(y1);
  }else{
    *rydata = y1;
    myfree(y2);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_WN_MD_STIM_FOR_TRIAL_PAST                     */
/*                                                                           */
/*  Return the stimulus sequence for the "k"the trial.  This routine covers  */
/*  the following p-files:                                                   */
/*    vbin2md.p                                                              */
/*                                                                           */
/*  Notes:                                                                   */
/*  - xdata values are returned in milliseconds.                             */
/*  - Page 0 shows first, thus, ydata[i] is for pattern0 when i%2=0.         */
/*                                                                           */
/*****************************************************************************/
void get_wn_md_stim_for_trial_PAST(nd,k,awake_flag,rxdata,rydata,rn)
     struct ndata_struct *nd;
     int k,awake_flag;
     float **rxdata,**rydata;
     int *rn;
{
  int i;
  int framerate,jumpseed,dtframe,period,vdpflag,mu,sigma,mnj0,mnj1,mxj0,mxj1;
  int n,count,minjump,maxjump;
  int moveflag,pflag,page,tjump;
  float dt,*xdata,*ydata;

  get_wn_two_stim_params(nd,k,awake_flag,&framerate,&jumpseed,&dtframe,&period,
			 &vdpflag,&mu,&sigma,&mnj0,&mnj1,&mxj0,&mxj1);
  n = framerate * period; /* Number of stimulus frames. */
  dt = get_wn_actual_dt(framerate,awake_flag);
  if (jumpseed > 0) /*** Initialize random number generator. ***/
    jumpseed = -jumpseed;
  noise_util_ran2(&jumpseed);

  tjump = -1; /* avoid -Wall warning */

  xdata = (float *)myalloc(n*sizeof(float));
  ydata = (float *)myalloc(n*sizeof(float));
  page = 0;
  pflag = 1;
  count = moveflag = 0;
  for(i=0;i<n;i++){
    xdata[i] = (float)((double)i*dt);

    if (page){
      minjump = mnj1;
      maxjump = mxj1;
    }else{
      minjump = mnj0;
      maxjump = mxj0;
    }

    if ((count == 0)&&(moveflag)){ /*** Pick a delay. ***/
      count = (int)(noise_util_ran2(&jumpseed)*4);
      if (count == 3) /*** Convert 0,1,2,3 --> 0,1,2,4 ***/
	count = 4;
      moveflag = 0;
    }
    if ((count == 0)&&(!moveflag)){ /*** Start a new movement. ***/
      if (vdpflag == 2)
	pflag = 1-pflag;
      else
	pflag = vdpflag;
      count = (int)(noise_util_ran2(&jumpseed)*9);
      if (count == 0){
	tjump = 0;
	count = 8;
      }else{
	if (count < 5)
	  tjump = minjump;
	else{
	  tjump = maxjump;
	  count -= 4;
	}
	if (count == 3)  /*** Convert 1,2,3,4 --> 1,2,6,4 ***/
	  count = 6;
      }
      moveflag = 1;
    }
    if (count > 0){ /*** Finish the current period. ***/
      if (pflag==page){
	ydata[i] = tjump;
	count -= 1;
      }else{
	ydata[i] = 0;
      }
    }else{
      printf("WYETH count = 0\n");  /*** THIS NEVER HAPPENS. ***/
      ydata[i] = ydata[i] - 1;
    }

    if (page)
      page = 0;
    else
      page = 1;
  }
  *rxdata = xdata; *rydata = ydata; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*   **** FUTURE              GET_WN_MD_STIM_FOR_TRIAL                       */
/*                                                                           */
/*  Return the stimulus sequence for the "k"the trial.  This routine covers  */
/*  the following p-files:                                                   */
/*    vbin2md.p                                                              */
/*                                                                           */
/*  Notes:                                                                   */
/*  - xdata values are returned in milliseconds.                             */
/*  - Page 0 shows first, thus, ydata[i] is for pattern0 when i%2=0.         */
/*                                                                           */
/*****************************************************************************/
void get_wn_md_stim_for_trial(nd,k,awake_flag,rxdata,rydata,rn)
     struct ndata_struct *nd;
     int k,awake_flag;
     float **rxdata,**rydata;
     int *rn;
{
  int i;
  int framerate,jumpseed,dtframe,period,vdpflag,mu,sigma,mnj0,mnj1,mxj0,mxj1;
  int n,count,pcount,minjump,maxjump;
  int moveflag,pflag,page,tjump;
  float dt,*xdata,*ydata;

  get_wn_two_stim_params(nd,k,awake_flag,&framerate,&jumpseed,&dtframe,&period,
			 &vdpflag,&mu,&sigma,&mnj0,&mnj1,&mxj0,&mxj1);
  n = framerate * period; /* Number of stimulus frames. */
  dt = get_wn_actual_dt(framerate,awake_flag);
  if (jumpseed > 0) /*** Initialize random number generator. ***/
    jumpseed = -jumpseed;
  noise_util_ran2(&jumpseed);

  tjump = -1; /* avoid -Wall warning */

  xdata = (float *)myalloc(n*sizeof(float));
  ydata = (float *)myalloc(n*sizeof(float));
  page = 0;
  pflag = 1;
  count = pcount = moveflag = 0;
  for(i=0;i<n;i++){
    xdata[i] = (float)((double)i*dt);

    if (page){
      minjump = mnj1;
      maxjump = mxj1;
    }else{
      minjump = mnj0;
      maxjump = mxj0;
    }

    if ((count == 0)&&(moveflag)){ /*** Pick a delay. ***/
      count = (int)(noise_util_ran2(&jumpseed)*4);
      /*** 0123->0124, 0123->0123 if same pat ***/
      if ((count == 3)&&((pcount == 2)||(vdpflag!=2)))
	count = 4;
      moveflag = 0;
      tjump = 0; /*** WYETH, THIS WAS LEFT OUT ON THE VERY FIRST p-FILE. ***/
    }
    if ((count == 0)&&(!moveflag)){ /*** Start a new movement. ***/
      if (vdpflag == 2){
	if (pcount == 2){
	  pflag = 1-pflag;
	  pcount = 0;
	}
      }else if (vdpflag == 3)
	pflag = 1-pflag;
      else
	pflag = vdpflag;
      count = (int)(noise_util_ran2(&jumpseed)*7) + 1;
      if (count == 0){
	tjump = 0;
	count = 4;
      }else{
	if (count < 4)
	  tjump = minjump;
	else{
	  tjump = maxjump;
	  count -= 3;
	}
	if (count == 3)  /*** Convert 1,2,3 --> 1,2,4 ***/
	  count = 4;
      }
      moveflag = 1;
      pcount += 1;
    }
    if (count > 0){ /*** Finish the current period. ***/
      if (pflag==page){
	if (moveflag) /*** WYETH - THIS CONDITION LEFT OUT ***/
	  ydata[i] = tjump;
	else
	  ydata[i] = 0;
	count -= 1;
      }else{
	ydata[i] = 0;
      }
    }else{
      printf("WYETH count = 0\n");  /*** THIS NEVER HAPPENS. ***/
      ydata[i] = ydata[i] - 1;
    }

    if (page)
      page = 0;
    else
      page = 1;
  }
  *rxdata = xdata; *rydata = ydata; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_WY2PA_STIM_FOR_TRIAL                       */
/*                                                                           */
/*****************************************************************************/
void get_wy2pa_stim_for_trial(nd,k,pfile,rpstate,rtanti,rns,radur,rpdur,rdt)
     struct ndata_struct *nd;
     int k;
     char *pfile;
     int **rpstate,**rtanti,*rns,*radur,*rpdur;
     float *rdt;
{
  int framerate,seed,dur,adur,pdur,tov;
  float fdt;

  get_wy2pa_stim_params(nd,k,&framerate,&seed,&dur,&adur,&pdur,&tov);

  if (strcmp(pfile,"MODEL.p")==0)
    fdt = 1000.0/framerate;
  else
    fdt = get_wn_actual_dt(framerate,0);

  get_wy2pa_stim_for_params(pfile,seed,framerate,dur,adur,pdur,tov,
				       fdt,rdt,rpstate,rtanti,rns);
  *radur = adur; *rpdur = pdur;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_WYTFV_STIM_FOR_TRIAL                       */
/*                                                                           */
/*****************************************************************************/
void get_wytfv_stim_for_trial(nd,k,pfile,rstima,rstimp,rns,rdt)
     struct ndata_struct *nd;
     int k;
     char *pfile;
     int **rstima,**rstimp,*rns;
     float *rdt;
{
  int framerate,seed,dur,zflag,t1,t2,t3,ntf,maxtf;
  float fdt;
  
  get_wytfv_stim_params(nd,k,&framerate,&seed,&dur,&zflag,&t1,&t2,&t3,&ntf,
			&maxtf);

  if (strcmp(pfile,"MODEL.p")==0)
    fdt = 1000.0/framerate;
  else
    fdt = get_wn_actual_dt(framerate,0);

  get_wytfv_stim_for_params(pfile,seed,framerate,dur,zflag,t1,t2,t3,ntf,maxtf,
			    fdt,rdt,rstima,rstimp,rns);
}
/**************************************-**************************************/
/*                                                                           */
/*                          NU_RANSTEP_GET_STIM_GROUP                        */
/*                                                                           */
/*****************************************************************************/
void nu_ranstep_get_stim_group(nd,ndg,k,iflag,norm_flag,rstim,rns,rdt)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k;           // Group number
     int iflag;       // 1-use integers for stimuli
     int norm_flag;   // 1-zscore
     float ***rstim;
     int **rns;
     float *rdt;
{
  int i,j;
  int n,jz,jn,seed,dti,tns,*ns,ti;
  float **stim,sf,jsz,stn,fps,*newstim,*tstim,dt,ssamp;

  printf("  NU_RANSTEP_GET_STIM_GROUP\n");

  n = ndg->cnt[k];
  stim = (float **)myalloc(n*sizeof(float *));
  ns = (int *)myalloc(n*sizeof(int));

  if (strcmp(nd->class,"MODEL")==0){
    ssamp = ndata_get_const_param_float(nd,"stim_samp");
    fps = ssamp; // Frames per second

    //printf(" ****** WYETH WARNING - stim_samp =  %f\n",ssamp);
    //fps = 500.0;
    //printf(" ****** WYETH WARNING - hardcoded to %f fps\n",fps);
  }else{
    fps = ndata_get_const_param_float(nd,"labrcon_fps");
  }
  dt = 1.0/fps * 1000.0;  // Send 'dt' back in milliseconds

  for(i=0;i<n;i++){
    ti = ndg->tnum[k][i];  // Trial index
    sf  = ndata_get_vc_param_float_or_exit(nd,ti,"sf",
					   "NU_RANSTEP_GET_STIM_GROUP");
    jsz = ndata_get_vc_param_float_or_exit(nd,ti,"step_size",
					   "NU_RANSTEP_GET_STIM_GROUP");
    jz  = ndata_get_vc_param_int_or_exit(nd,ti,"step_zero",
					 "NU_RANSTEP_GET_STIM_GROUP");
    jn  = ndata_get_vc_param_int_or_exit(nd,ti,"step_n",
					 "NU_RANSTEP_GET_STIM_GROUP");
    seed= ndata_get_vc_param_int_or_exit(nd,ti,"seed",
					 "NU_RANSTEP_GET_STIM_GROUP");
    stn = ndata_get_vc_param_float_or_exit(nd,ti,"stn",
					   "NU_RANSTEP_GET_STIM_GROUP");
    dti = ndata_get_vc_param_int_or_exit(nd,ti,"dt",
					 "NU_RANSTEP_GET_STIM_GROUP");

    //iflag = 0;
    tstim = stm_ranstep_get_stim(NULL,sf,jsz,jz,jn,seed,stn,fps,dti,iflag,
				 &tns);

    // append_farray_plot("zzz.stim","stim",tstim,tns,1);


    //if (jsz != 0)
    //exit_error("NU_RANSTEP_GET_STIM_GROUP","WYETH - jsz != 0, not imp'd");

    if (iflag == 0){
      newstim = get_farray(tns);
      // 0,1,2,3,4,5,6,7  -->  0,-1,-2,-3,-4=4, 3, 2, 1
      for(j=1;j<tns;j++){
	if (tstim[j] < 0.5)
	  newstim[j] = -tstim[j];
	else if (tstim[j] == 0.5) // Contrast reversal, ambiguous motion
	  newstim[j] = 0.0;
	else
	  newstim[j] = 1.0 - tstim[j];
      }
      //multiply_farray(newstim,tns,-1.0);
      myfree(tstim);
      tstim = newstim;

      if (norm_flag > 0){
	z_score_farray(tstim,tns);  // WYETH - is this OK?
      }
    }
    stim[i] = tstim;
    ns[i] = tns;
  }

  // WYETH - assuming that all 'dti' values are the same
  dt = dt * (float)dti;


  
  printf("    dti = %d\n",dti);
  printf("     dt = %f\n",dt);



  *rstim = stim;
  *rns = ns;
  *rdt = dt;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_STIM_AT_RESOLUTION_SIMPLE                     */
/*                                                                           */
/*  Return the stimulus at the specified sampling resolution.                */
/*  "res_in" and "res_out" are in samples per second.                        */
/*                                                                           */
/*  For example, res_in = 53.75 for 18.605 msec/frame stimulus.              */
/*               res_out = 1000.0 to return msec resolution.                 */
/*                                                                           */
/*****************************************************************************/
void get_stim_at_resolution_simple(data,n,res_in,res_out,rdata,rn)
     float *data;
     int n;
     float res_in,res_out;   // (samples/s)
     float **rdata;
     int *rn;
{
  int i,k;
  int nn;
  float *ydata;

  nn = (int)(((float)n-0.5) * (res_out/res_in)); // Length of new array
  ydata = (float *)myalloc(nn*sizeof(float));
  for(i=0;i<nn;i++){
    k = (int)(0.5 + (float)i/(res_out/res_in)); // Index in old array
    ydata[i] = data[k];
    if (k >= n)
      exit_error("GET_STIM_AT_RESOLUTION_SIMPLE","This should never happen");
  }

  *rn = nn;
  *rdata = ydata;
}
/**************************************-**************************************/
/*                                                                           */
/*                      GET_STIM_AT_RESOLUTION_PIECEWISE                     */
/*                                                                           */
/*****************************************************************************/
void get_stim_at_resolution_piecewise(data,n,res_in,res_out,rdata,rn)
     float *data;
     int n;
     float res_in,res_out,**rdata;
     int *rn;
{
  int i,k;
  int nn;
  float x,y,*ydata;

  nn = (int)(((float)n-0.5) * (res_out/res_in)); // Length of new array
  ydata = (float *)myalloc(nn*sizeof(float));
  for(i=0;i<nn;i++){
    x = (float)i/(res_out/res_in); // Time in old array
    k = (int)x;
    if (k<(n-1))
      y = data[k] + (x-(float)k)*(data[k+1]-data[k]);
    else
      y = data[n-1];
    ydata[i] = y;
    if (k >= n)
      exit_error("GET_STIM_AT_RESOLUTION_PIECEWISE","Should never happen");
  }
  *rn = nn;
  *rdata = ydata;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_STIM_AT_RESOLUTION_IMPULSE                    */
/*                                                                           */
/*  Like "get_stim_at_resolution_simple" but set only one value to stimulus  */
/*  value, don't use boxcars.                                                */
/*                                                                           */
/*****************************************************************************/
void get_stim_at_resolution_impulse(data,n,res_in,res_out,rdata,rn)
     float *data;
     int n;
     float res_in,res_out,**rdata;
     int *rn;
{
  int i,k;
  int nn;
  float *ydata;

  /*** EXPAND THE DATA ARRAY TO msec TIME SCALE ***/
  nn = (int)((float)n * (res_out/res_in));
  ydata = get_zero_farray(nn);
  for(i=0;i<n;i++){
    k = (int)((float)i*(res_out/res_in));
    if (k<nn)
      ydata[k] = data[i];
  }
  *rn = nn;
  *rdata = ydata;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_STIM_AT_RESOLUTION                          */
/*                                                                           */
/*  Return the stimulus at the specified sampling resolution.  "res_in" and  */
/*  "res_out" are in samples per second.                                     */
/*                                                                           */
/*  For example, res_in = 53.75 for 18.605 msec/frame stimulus.              */
/*               res_out = 1000.0 to return msec resolution.                 */
/*                                                                           */
/*****************************************************************************/
void get_stim_at_resolution(data,n,res_in,res_out,rdata,rn)
     float *data;
     int n;
     float res_in,res_out,**rdata;
     int *rn;
{
  int i;
  int nn,expf,nexp;
  float *expdata,*xexpdata,*ydata,*xdata;

  //
  // Expand the data to beyond the desired resolution
  //
  expf = (int)(1.0 + res_out/res_in);
  //printf("    Expanding input data by %d before subsampling. \n",expf);
  expdata = fft_interpolate(data,n,expf);
  nexp = expf*n;
  xexpdata = (float *)myalloc(nexp*sizeof(float));
  for(i=0;i<nexp;i++)
    xexpdata[i] = (float)i * res_out/(res_in*(float)expf);

  //
  // Subsample the data at the desired resolution
  //
  nn = (int)((float)n * res_out/res_in);
  xdata = (float *)myalloc(nn*sizeof(float));
  for(i=0;i<nn;i++)
    xdata[i] = (float)i;
  ydata = resample_farray(xexpdata,expdata,nexp,xdata,nn);
  myfree(expdata); myfree(xexpdata); myfree(xdata);

  *rn = nn; *rdata = ydata;
}
/**************************************-**************************************/
/*                                                                           */
/*                              ADJUST_TERNARY_STIM                          */
/*                                                                           */
/*****************************************************************************/
void adjust_ternary_stim(nd,k,awake_flag,data,n,pfile)
     struct ndata_struct *nd;
     int k,awake_flag;
     float *data;
     int n;
     char pfile[];
{
  int i;
  int maxjump,flag;

  if ((strcmp(pfile,"wypnori.p")==0)||(strcmp(pfile,"wypnori3.p")==0)||
      (strcmp(pfile,"wypnoril.p")==0)||(strcmp(pfile,"wypnphal.p")==0)||
      (strcmp(pfile,"wyph3.p")==0)||(strcmp(pfile,"wy2flat2.p")==0)||
      (strcmp(pfile,"wypnpha.p")==0)||(strcmp(pfile,"wypnpha3.p")==0)){
    maxjump = 1;
  }else if (awake_flag){
    flag = ndata_get_vc_param_int(nd,k,"maxjump",&maxjump);
    if (!flag) exit_error("ADJUST_TERNARY_STIM","Param maxjump not found");
  }else{
    flag = ndata_get_vc_param_int(nd,k,"max_jump",&maxjump);
    if (!flag){
      /*printf("  NOTE:  max_jump not found, looking for maxj\n");*/
      flag = ndata_get_vc_param_int(nd,k,"maxj",&maxjump);
      if (!flag) exit_error("ADJUST_TERNARY_STIM","Param maxj not found");
    }
  }
  for(i=0;i<n;i++)
    data[i] = (float)my_rint(data[i]/maxjump);

  /*append_farray_plot("zzz.stim","stim",data,n,1);*/
}
/**************************************-**************************************/
/*                                                                           */
/*                          ADJUST_DOUBLE_TERNARY_STIM                       */
/*                                                                           */
/*  This handles the case where two ternary stimuli with possibly different  */
/*  jumpsizes are interleaved, as in "vbin2.p" or "vbin2md.p".               */
/*                                                                           */
/*****************************************************************************/
void adjust_double_ternary_stim(nd,k,awake_flag,data,n)
     struct ndata_struct *nd;
     int k,awake_flag;
     float *data;
     int n;
{
  int i,j;
  int mxj0,mxj1,mnj0,mnj1,flag;

  if (awake_flag){
    flag = ndata_get_vc_param_int(nd,k,"maxjump0",&mxj0);
    if (!flag) exit_error("ADJUST_TERNARY_STIM","Param maxjump0 not found");
    flag = ndata_get_vc_param_int(nd,k,"maxjump1",&mxj1);
    if (!flag) exit_error("ADJUST_TERNARY_STIM","Param maxjump1 not found");
    flag = ndata_get_vc_param_int(nd,k,"minjump0",&mnj0);
    if (!flag) exit_error("ADJUST_TERNARY_STIM","Param minjump0 not found");
    flag = ndata_get_vc_param_int(nd,k,"minjump1",&mnj1);
    if (!flag) exit_error("ADJUST_TERNARY_STIM","Param minjump1 not found");
  }else{
    flag = ndata_get_vc_param_int(nd,k,"max_jmp0",&mxj0);
    if (!flag) exit_error("ADJUST_TERNARY_STIM","Param max_jmp0 not found");
    flag = ndata_get_vc_param_int(nd,k,"max_jmp1",&mxj1);
    if (!flag) exit_error("ADJUST_TERNARY_STIM","Param max_jmp1 not found");
    flag = ndata_get_vc_param_int(nd,k,"min_jmp0",&mnj0);
    if (!flag) exit_error("ADJUST_TERNARY_STIM","Param min_jmp0 not found");
    flag = ndata_get_vc_param_int(nd,k,"min_jmp1",&mnj1);
    if (!flag) exit_error("ADJUST_TERNARY_STIM","Param min_jmp1 not found");
  }
  for(i=0;i<n;i++){
    j = my_rint(data[i]);
    if ((j == mxj0)||(j == mxj1))
      data[i] = 1;
    else if ((j == mnj0)||(j == mnj1))
      data[i] = -1;
    else
      data[i] = 0;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                       GET_NEW_NOISE_SIMPLE_STIM_GROUP                     */
/*                                                                           */
/*  Notes:                                                                   */
/*  - Assuming same dt for all stimuli.                                      */
/*  - binflag = 1-Binary, 3-Ternary, 4-Gaussian, 9-Position (not jump)       */
/*  - Dangerous to have 'binflag' 1 if stimulus is ternary, because mean     */
/*    value is used to pick 1's and 0's for binary.                          */
/*                                                                           */
/*****************************************************************************/
void get_new_noise_simple_stim_group(nd,ndg,k,awake_flag,binflag,rstim,rns,rdt)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k,awake_flag,binflag;
     float ***rstim;
     int **rns;
     float *rdt;
{
  int i,j;
  int flag,*ns,n,invflag,maxj,phase;
  float *tx,**stim,dt,mean,sdev;
  char *pfile,*invert,*filename;

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0)
    exit_error("GET_NEW_NOISE_SIMPLE_STIM_GROUP","Param pfile not found");
  if ((strcmp(pfile, "wypnori.p")==0)||(strcmp(pfile, "wypnpha.p")==0)||
      (strcmp(pfile,"wypnori3.p")==0)||(strcmp(pfile,"wypnpha3.p")==0)||
      (strcmp(pfile,   "wy4st.p")==0)||(strcmp(pfile,   "wyone.p")==0)||
      (strcmp(pfile,  "wyonem.p")==0)||
      (strcmp(pfile, "wyone23.p")==0)||(strcmp(pfile,"wypnoril.p")==0)||
      (strcmp(pfile,  "wy4stg.p")==0)||(strcmp(pfile,"wypnphal.p")==0)||
      (strcmp(pfile, "wy2flat.p")==0)||(strcmp(pfile,    "x4st.p")==0)||
      (strcmp(pfile,   "x4stp.p")==0)||(strcmp(pfile,   "x4std.p")==0)||
      (strcmp(pfile,  "x4stp2.p")==0)||
      (strcmp(pfile,  "wy4sum.p")==0)||(strcmp(pfile,"wyphstat.p")==0)||
      (strcmp(pfile,   "wycon.p")==0)||(strcmp(pfile,  "wymtf5.p")==0)||
      (strcmp(pfile, "wymtf5p.p")==0)||(strcmp(pfile,"wymtf4ph.p")==0)||
      (strcmp(pfile,"wymat160.p")==0)||(strcmp(pfile,   "wyall.p")==0)||
      (strcmp(pfile,  "wymtf9.p")==0)||(strcmp(pfile,  "wymsf9.p")==0)||
      (strcmp(pfile, "wymtf7t.p")==0)||
      (strcmp(pfile,  "wymsf6.p")==0)||(strcmp(pfile, "wymtf9p.p")==0)||
      (strcmp(pfile,    "wydt.p")==0)||
      (strcmp(pfile, "wymsf11.p")==0)||(strcmp(pfile,  "wyamc5.p")==0)||
      (strcmp(pfile,  "wyamc7.p")==0)||(strcmp(pfile,"wyamcsf5.p")==0)||
      (strcmp(pfile,   "wyamc.p")==0)||(strcmp(pfile,"wyamcpha.p")==0)||
      (strcmp(pfile, "wyphasm.p")==0)||(strcmp(pfile,"wyphasm8.p")==0)||
      (strcmp(pfile,"wyphasms.p")==0)||
      (strcmp(pfile,   "wyph3.p")==0)||(strcmp(pfile,"wy2flat2.p")==0)||
      (strcmp(pfile,"MODEL_wyamc.p")==0)||(strcmp(pfile,"MODEL_new.p")==0)||
      (strcmp(pfile,  "wymatm.p")==0)||(strcmp(pfile, "wyphas8.p")==0)){
    flag = ndata_get_const_param(nd,"filename",&filename);
    if (flag > 0){
      if (strcmp(filename,"452l032.p07")==0)
	exit_error("GET_NEW_NOISE_SIMPLE_STIM_GROUP","Special case");
    }

    n = ndg->cnt[k];
    stim = (float **)myalloc(n*sizeof(float *));
    ns = (int *)myalloc(n*sizeof(int));

    invflag = ndata_get_const_param(nd,"invert",&invert);
    if (invflag == 1){
      myfree(invert);
      printf("    Invert flag == 1\n");
    }

    dt = -1.0; /* Avoid -Wall warning */
    for(i=0;i<n;i++){
      j = ndg->tnum[k][i];
      if ((strcmp(pfile,"wyamc5.p")==0)||(strcmp(pfile,"wyamc7.p")==0)||
	  (strcmp(pfile,"MODEL_wyamc.p")==0)){
	get_wy2pat_stim_for_trial_03(nd,j,awake_flag,1,pfile,&tx,
				     &stim[i],&ns[i]);
      }else if (strcmp(pfile,"wy2flat2.p")==0){
	get_wy2pat_stim_for_trial_04(nd,j,awake_flag,pfile,&tx,&stim[i],
				     &ns[i]);
      }else{
	get_wy2pat_stim_for_trial_01(nd,j,awake_flag,pfile,&tx,&stim[i],
				     &ns[i]);
      }
      dt = (tx[ns[0]-1] - tx[0])/(float)(ns[0]-1); /* Do this before free. */
      myfree(tx);

      /*append_farray_plot("zz.stim","before",stim[i],ns[i],1);*/

      if (invflag)
	multiply_farray(stim[i],ns[i],-1.0);

      if (binflag == 1){
	mean_sdev_farray(stim[i],ns[i],&mean,&sdev);
	binary_threshold_farray(stim[i],ns[i],mean);
      }else if (binflag == 3){
	adjust_ternary_stim(nd,j,awake_flag,stim[i],ns[i],pfile);
      }else if (binflag == 9){ /* Convert stimulus to position, from jump */
	/*** WYETH ***/
	flag = ndata_get_vc_param_int(nd,j,"phase",&phase);
	if (!flag) exit_error("***","Param phase not found");
	flag = ndata_get_vc_param_int(nd,j,"maxj",&maxj);
	if (!flag) exit_error("***","Param maxj not found");
	convert_jump_to_position(stim[i],ns[i],phase,maxj,1024);
	printf(" WYETH - converting to position \n");
      }
    }
    *rstim = stim; *rns = ns; *rdt = dt;
    myfree(pfile);
  }else{
    printf("*** pfile = %s\n",pfile);
    exit_error("GET_NEW_NOISE_SIMPLE_STIM_GROUP","Unrecognized pfile");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                        GET_NOISE_SIMPLE_STIM_GROUP                        */
/*                                                                           */
/*  Notes:                                                                   */
/*  - Assuming same dt for all stimuli.                                      */
/*  - binflag = 1-Binary, 3-Ternary, 4-Gaussian                              */
/*                                                                           */
/*****************************************************************************/
void get_noise_simple_stim_group(nd,ndg,k,awake_flag,binflag,rstim,rns,rdt)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k,awake_flag,binflag;
     float ***rstim;
     int **rns;
     float *rdt;
{
  int i,j;
  int flag,*ns,n,negflag,invflag,nsl;
  float *tx,**stim,dt,mean,sdev;
  char *pfile,*invert,*stimtype,**slist;

  if ((strcmp(nd->class,"labradoodle_01")==0) || 
      (strcmp(nd->class,"MODEL")==0)){
    flag = ndata_get_const_param(nd,"stim_type",&stimtype);
    if (flag == 0)
      exit_error("GET_NOISE_SIMPLE_STIM_GROUP",
		 "Const param 'stim_type' not found");
    get_items_from_string(stimtype,&slist,&nsl);
    myfree(stimtype);
    stimtype = strdup(slist[0]);
    free_2d_carray(slist,nsl);

    if (strcmp(stimtype,"ranstep")==0){
      int iflag;

      iflag = 1; // Use indices for stimulus, not step sizes
      nu_ranstep_get_stim_group(nd,ndg,k,iflag,0,rstim,rns,rdt);

      printf("NOISE UTIL get_noise_simple_stim_group    %f\n",*rdt);

    }else{
      printf("stimtype = %s\n",stimtype);
      exit_error("GET_NOISE_STIM_GROUP","Unknown stimtype");
    }
  }else{
    flag = ndata_get_const_param(nd,"pfile",&pfile);
    if (flag == 0)
      exit_error("GET_NOISE_SIMPLE_STIM_GROUP","Const param pfile not found");
    if ((strcmp(pfile,"wnall100.p")==0)||(strcmp(pfile,"wncon100.p")==0)||
	(strcmp(pfile,"wnmat050.p")==0)||(strcmp(pfile,"wnmat053.p")==0)||
	(strcmp(pfile,"wnmat080.p")==0)||(strcmp(pfile,"wnmat100.p")==0)||
	(strcmp(pfile,"wnmdt100.p")==0)||(strcmp(pfile,"wnone053.p")==0)||
	(strcmp(pfile,"wnone080.p")==0)||(strcmp(pfile,"wnone100.p")==0)||
	(strcmp(pfile, "wnposdt.p")==0)||(strcmp(pfile,"wnsiz100.p")==0)||
	(strcmp(pfile,   "wngnv.p")==0)||(strcmp(pfile, "wynoise.p")==0)||
	(strcmp(pfile,  "MODEL.p")==0) ||(strcmp(pfile, "wnnovar.p")==0)||
	(strcmp(pfile, "wnjump.p")==0)||
	(strcmp(pfile,   "vbin1.p")==0)||(strcmp(pfile,"vbinmat.p")==0)){
      n = ndg->cnt[k];
      stim = (float **)myalloc(n*sizeof(float *));
      ns = (int *)myalloc(n*sizeof(int));
      
      invflag = ndata_get_const_param(nd,"invert",&invert);
      if (invflag == 1){
	myfree(invert);
	printf("    Invert flag == 1\n");
      }
      
      dt = -1.0; // Avoid -Wall warning
      for(i=0;i<n;i++){
	j = ndg->tnum[k][i];
	get_wn_stim_for_trial(nd,j,awake_flag,pfile,&tx,&stim[i],&ns[i]);
	dt = (tx[ns[0]-1] - tx[0])/(float)(ns[0]-1); /* Do this before free. */
	myfree(tx);
	
	negflag = 0;
	if (strcmp(pfile,"vbin1.p")==0)
	  negflag =  nd->t[j].tcode % 2; // Null direction
	else if ((strcmp(pfile,"wngnv.p")==0)||(strcmp(pfile,"wynoise.p")==0)||
		 (strcmp(pfile,"wnnovar.p")==0)||
		 (strcmp(pfile,"wnone080.p")==0))
	  negflag = ndata_trial_get_rcode_chan(&(nd->t[j]),"unit0"); // Null
	if ((negflag&&(!invflag))||((!negflag)&&invflag))
	  multiply_farray(stim[i],ns[i],-1.0);
	
	if (binflag == 1){
	  mean_sdev_farray(stim[i],ns[i],&mean,&sdev);
	  binary_threshold_farray(stim[i],ns[i],mean);
	}else if (binflag == 3){
	  adjust_ternary_stim(nd,j,awake_flag,stim[i],ns[i],pfile);
	}
      }
      *rstim = stim; *rns = ns; *rdt = dt;
      myfree(pfile);
    }else
      get_new_noise_simple_stim_group(nd,ndg,k,awake_flag,binflag,rstim,rns,
				      rdt);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                        GET_NOISE_SIMPLE_MD_STIM_GROUP                     */
/*                                                                           */
/*  Notes:                                                                   */
/*  - Assuming same dt for all stimuli.                                      */
/*                                                                           */
/*****************************************************************************/
void get_noise_simple_md_stim_group(nd,ndg,k,awake_flag,rstim,rns,rdt)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k,awake_flag;
     float ***rstim;
     int **rns;
     float *rdt;
{
  int i,j;
  int flag,*ns,n;
  float *tx,**stim,dt;
  char *pfile;

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0)
    exit_error("GET_NOISE_SIMPLE_MD_STIM_GROUP","Const param pfile not found");
  if (strcmp(pfile,"vbin2md.p")==0){
    n = ndg->cnt[k];
    stim = (float **)myalloc(n*sizeof(float *));
    ns = (int *)myalloc(n*sizeof(int));
    dt = -1.0; /* Avoid -Wall warning */
    for(i=0;i<n;i++){
      j = ndg->tnum[k][i];
      get_wn_md_stim_for_trial(nd,j,awake_flag,&tx,&stim[i],&ns[i]);
      dt = (tx[ns[0]-1] - tx[0])/(float)(ns[0]-1); /* Do this before free. */
      myfree(tx);
      adjust_double_ternary_stim(nd,j,awake_flag,stim[i],ns[i]);
    }
    *rstim = stim; *rns = ns; *rdt = dt;
    myfree(pfile);
  }else
    exit_error("GET_NOISE_SIMPLE_MD_STIM_GROUP","Unrecognized pfile");
}
/**************************************-**************************************/
/*                                                                           */
/*                      GET_NOISE_SIMPLE_STIM_VEL_GROUP                      */
/*                                                                           */
/*  vflag:                                                                   */
/*  0 - simple 2-point derivative, first value set to zero.                  */
/*                                                                           */
/*****************************************************************************/
void get_noise_simple_stim_vel_group(nd,ndg,k,awake_flag,binflag,vflag,
				     rvel,rns,rdt)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k,awake_flag,binflag,vflag;
     float ***rvel;
     int **rns;
     float *rdt;
{
  int i;
  int n,*ns;
  float **stim,**vel;
  
  get_noise_simple_stim_group(nd,ndg,k,awake_flag,binflag,&stim,&ns,rdt);

  n = ndg->cnt[k];
  vel = (float **)myalloc(n*sizeof(float *));
  if (vflag == 0){
    for(i=0;i<n;i++){
      vel[i] = derivative_simple_farray(stim[i],ns[i]);
      round_to_int_farray(vel[i],ns[i]);
    }
    free_2d_farray(stim,n);
  }else
    exit_error("GET_NOISE_SIMPLE_STIM_VEL_GROUP","Option not implemented");
  
  *rvel = vel; *rns = ns;
}
/**************************************-**************************************/
/*                                                                           */
/*                      GET_SPECIAL_STIM_GROUP_452L032P07                    */
/*                                                                           */
/*****************************************************************************/
void get_special_stim_group_452l032p07(nd,ndg,k,resolution,fft_flag,
					  norm_flag,rstim,rns,rdt)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k;
     float resolution;
     int norm_flag,fft_flag;
     float ***rstim;
     int **rns;
     float *rdt;
{
  int i,j,l;
  int seed,ng,n,nstim,*allns,*ns;
  float dt,tf,**allstim,**stim,res_in;

  printf("  GET_SPECIAL_STIM_GROUP_452l032p07\n");

  dt = get_wn_actual_dt(100,0);

  seed = 0; /*** Use 300 calls to gasdev ***/
  for(i=0;i<300;i++)
    tf = noise_util_gasdev(&seed);

  /*** Create stimulus for each trial distrib=1, dt=1, mu=0.0, sigma=64.0 ***/
  n = nd->ntrial;
  nstim = 30 * 100;
  allstim = get_2d_farray(n,nstim);
  allns = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++){ /* for each trial in the *file* */
    allstim[i] = (float *)myalloc(nstim*sizeof(float));
    for(j=0;j<nstim;j++){
      tf = 64.0 * noise_util_gasdev(&seed);
      if (tf > 0.0)
	allstim[i][j] = (float)((int)(0.5+tf));
      else
	allstim[i][j] = (float)((int)(-0.5+tf));
    }
    allns[i] = nstim;
  }

  ng = ndg->cnt[k];
  stim = get_2d_farray(ng,nstim);
  ns = (int *)myalloc(ng*sizeof(int));
  for(i=0;i<ng;i++){
    res_in = 1000.0/dt;
    l = ndg->tnum[k][i];
    if (fft_flag == 0)
      get_stim_at_resolution_simple(allstim[l],allns[l],res_in,resolution,
				    &stim[i],&ns[i]);
    else if (fft_flag == 1)
      get_stim_at_resolution(allstim[l],allns[l],res_in,resolution,
			     &stim[i],&ns[i]);
    else if (fft_flag == 2)
      get_stim_at_resolution_impulse(allstim[l],allns[l],res_in,resolution,
				     &stim[i],&ns[i]);
    else if (fft_flag == 3)
      get_stim_at_resolution_piecewise(allstim[l],allns[l],res_in,resolution,
				     &stim[i],&ns[i]);
    else
      exit_error("GET_SPECIAL_STIM_GROUP_452L032P07","Unknown fft_flag value");

    if (norm_flag > 0)
      norm_variance_farray(stim[i],ns[i]);
  }

  free_2d_farray(allstim,n);
  myfree(allns);

  *rstim = stim; *rns = ns; *rdt = dt;
}
/**************************************-**************************************/
/*                                                                           */
/*                     GET_WY4ST_CRACK_PARAMS_452L004P17                     */
/*                                                                           */
/*  Not yet.                                                                 */
/*                                                                           */
/*****************************************************************************/
void get_wy4st_crack_params_452l004p17(nstim,stimlist,nlist,escflag,dstate,
				       nss,rlostlist,rnlost,rdwell,rdt)
     int *nstim,**stimlist,**nlist,*escflag,*dstate,*nss,**rlostlist,*rnlost;
     int *rdwell;
     float *rdt;
{
  int *lostlist;

  *nstim = 39;
  lostlist = get_zero_iarray(*nstim);
  *stimlist = (int *)myalloc(*nstim*sizeof(int));
  *nlist = (int *)myalloc(*nstim*sizeof(int));
  /* For dwell 8, 10, or 16, n = 38, 30, or 19 */
  (*stimlist)[0] = 0  ; (*nlist)[0] = 19; /* (3*100 - 1)/16 + 1 */
  (*stimlist)[1] =   2; (*nlist)[1] = 19;
  (*stimlist)[2] =   2; (*nlist)[2] = 188;
  (*stimlist)[3] = 0  ; (*nlist)[3] = 188; /*188 = (30*100 - 1)/16 + 1*/
  (*stimlist)[4] =  1 ; (*nlist)[4] = 188;
  (*stimlist)[5] =   2; (*nlist)[5] = 188;
  (*stimlist)[6] = 0  ; (*nlist)[6] = 188;
  (*stimlist)[7] =  1 ; (*nlist)[7] = 188; /* Special reset seed to 0 here. */
  (*stimlist)[8] =   2; (*nlist)[8] = 188;

  (*stimlist)[9] = 0  ; (*nlist)[9] = 188;
  (*stimlist)[10]=   2; (*nlist)[10]= 188;
  (*stimlist)[11]=  1 ; (*nlist)[11]= 188;
  (*stimlist)[12]=  1 ; (*nlist)[12]= 188;
  (*stimlist)[13]=   2; (*nlist)[13]= 188;
  (*stimlist)[14]= 0  ; (*nlist)[14]= 188;
  (*stimlist)[15]=  1 ; (*nlist)[15]= 188;
  (*stimlist)[16]=   2; (*nlist)[16]= 188;
  (*stimlist)[17]= 0  ; (*nlist)[17]= 188;
  (*stimlist)[18]= 0  ; (*nlist)[18]= 188;
  (*stimlist)[19]=  1 ; (*nlist)[19]= 188;
  (*stimlist)[20]=   2; (*nlist)[20]= 188;
  (*stimlist)[21]= 0  ; (*nlist)[21]= 188;
  (*stimlist)[22]=  1 ; (*nlist)[22]= 188;
  (*stimlist)[23]=   2; (*nlist)[23]= 188;
  (*stimlist)[24]= 0  ; (*nlist)[24]= 188;
  (*stimlist)[25]=   2; (*nlist)[25]= 188;
  (*stimlist)[26]=  1 ; (*nlist)[26]= 188;
  (*stimlist)[27]=   2; (*nlist)[27]= 188;
  (*stimlist)[28]=  1 ; (*nlist)[28]= 188;
  (*stimlist)[29]= 0  ; (*nlist)[29]= 188;
  (*stimlist)[30]=  1 ; (*nlist)[30]= 188;
  (*stimlist)[31]=   2; (*nlist)[31]= 188;
  (*stimlist)[32]= 0  ; (*nlist)[32]= 188;
  (*stimlist)[33]=   2; (*nlist)[33]= 188;
  (*stimlist)[34]= 0  ; (*nlist)[34]= 188;
  (*stimlist)[35]=  1 ; (*nlist)[35]= 188;
  (*stimlist)[36]=  1 ; (*nlist)[36]= 188;
  (*stimlist)[37]= 0  ; (*nlist)[37]= 188;
  (*stimlist)[38]=   2; (*nlist)[38]= 188;

  *escflag = 1; /* Assume "esc" key pressed before this trial, 511 ==> -511 */
  *dstate = 30;
  *nss = 188;
  *rdwell = 16;
  *rdt = get_wn_actual_dt(100,0);
  *rlostlist = lostlist; *rnlost = 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                     GET_WY4ST_CRACK_PARAMS_452L013P07                     */
/*                                                                           */
/*  Cracked - all trials.                                                    */
/*                                                                           */
/*  Tried:                                                                   */
/*  - esc 1  30,30,30                                                        */
/*  - esc 1  30,30,300                                                       */
/*  - esc 1  30,30,375                                                       */
/*  - esc 1  30,300,300                                                      */
/*  - esc 1  30,375,375                                                      */
/*  - esc 1  38,38,38                                                        */
/*  - esc 1  38,38,375 - something here.                                     */
/*                                                                           */
/*****************************************************************************/
void get_wy4st_crack_params_452l013p07(nstim,stimlist,nlist,escflag,dstate,
				       nss,rlostlist,rnlost,rdwell,rdt)
     int *nstim,**stimlist,**nlist,*escflag,*dstate,*nss,**rlostlist,*rnlost;
     int *rdwell;
     float *rdt;
{
  int *lostlist;

  *nstim = 27;
  lostlist = get_zero_iarray(*nstim);
  *stimlist = (int *)myalloc(*nstim*sizeof(int));
  *nlist = (int *)myalloc(*nstim*sizeof(int));
  /* For dwell 8, 10, or 16, n = 38, 30, or 19 */
  (*stimlist)[0] = 0; (*nlist)[0] = 38; /* (3*100 - 1)/8 + 1 */
  (*stimlist)[1] = 0; (*nlist)[1] = 38;
  (*stimlist)[2] = 0; (*nlist)[2] = 375;

  (*stimlist)[3] =  1 ; (*nlist)[3] = 375; /*375 = (30*100 - 1)/8 + 1*/
  (*stimlist)[4] = 0  ; (*nlist)[4] = 375;
  (*stimlist)[5] =   2; (*nlist)[5] = 375;
  (*stimlist)[6] =   2; (*nlist)[6] = 375;
  (*stimlist)[7] =  1 ; (*nlist)[7] = 375; /* Special reset seed to 0 here. */
  (*stimlist)[8] = 0  ; (*nlist)[8] = 375;
  (*stimlist)[9] = 0  ; (*nlist)[9] = 375;
  (*stimlist)[10]=  1 ; (*nlist)[10]= 375;
  (*stimlist)[11]=   2; (*nlist)[11]= 375;
  (*stimlist)[12]=  1 ; (*nlist)[12]= 375;
  (*stimlist)[13]=   2; (*nlist)[13]= 375;
  (*stimlist)[14]= 0  ; (*nlist)[14]= 375;
  (*stimlist)[15]=  1 ; (*nlist)[15]= 375;
  (*stimlist)[16]=   2; (*nlist)[16]= 375;
  (*stimlist)[17]= 0  ; (*nlist)[17]= 375;
  (*stimlist)[18]=   2; (*nlist)[18]= 375;
  (*stimlist)[19]= 0  ; (*nlist)[19]= 375;
  (*stimlist)[20]=  1 ; (*nlist)[20]= 375;
  (*stimlist)[21]= 0  ; (*nlist)[21]= 375;
  (*stimlist)[22]=  1 ; (*nlist)[22]= 375;
  (*stimlist)[23]=   2; (*nlist)[23]= 375;
  (*stimlist)[24]= 0  ; (*nlist)[24]= 375;
  (*stimlist)[25]=  1 ; (*nlist)[25]= 375;
  (*stimlist)[26]=   2; (*nlist)[26]= 375;

  *escflag = 1; /* Assume "esc" key pressed before this trial, 511 ==> -511 */
  *dstate = 24;
  *nss = 375;
  *rdwell = 8;
  *rdt = get_wn_actual_dt(100,0);
  *rlostlist = lostlist; *rnlost = 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                     GET_WY4ST_CRACK_PARAMS_452L013P08                     */
/*                                                                           */
/*  Cracked - all trials.                                                    */
/*                                                                           */
/*  (This starts with the entire p07 experiment.)                            */
/*                                                                           */
/*****************************************************************************/
void get_wy4st_crack_params_452l013p08(nstim,stimlist,nlist,escflag,dstate,
				       nss,rlostlist,rnlost,rdwell,rdt)
     int *nstim,**stimlist,**nlist,*escflag,*dstate,*nss,**rlostlist,*rnlost;
     int *rdwell;
     float *rdt;
{
  int *lostlist;

  *nstim = 39;
  lostlist = get_zero_iarray(*nstim);
  *stimlist = (int *)myalloc(*nstim*sizeof(int));
  *nlist = (int *)myalloc(*nstim*sizeof(int));
  /* For dwell 8, 10, or 16, n = 38, 30, or 19 */
  (*stimlist)[0] = 0; (*nlist)[0] = 38; /* (3*100 - 1)/8 + 1 */
  (*stimlist)[1] = 0; (*nlist)[1] = 38;
  (*stimlist)[2] = 0; (*nlist)[2] = 375;

  (*stimlist)[3] =  1 ; (*nlist)[3] = 375; /*375 = (30*100 - 1)/8 + 1*/
  (*stimlist)[4] = 0  ; (*nlist)[4] = 375;
  (*stimlist)[5] =   2; (*nlist)[5] = 375;
  (*stimlist)[6] =   2; (*nlist)[6] = 375;
  (*stimlist)[7] =  1 ; (*nlist)[7] = 375; /* Special reset seed to 0 here. */
  (*stimlist)[8] = 0  ; (*nlist)[8] = 375;
  (*stimlist)[9] = 0  ; (*nlist)[9] = 375;
  (*stimlist)[10]=  1 ; (*nlist)[10]= 375;
  (*stimlist)[11]=   2; (*nlist)[11]= 375;
  (*stimlist)[12]=  1 ; (*nlist)[12]= 375;
  (*stimlist)[13]=   2; (*nlist)[13]= 375;
  (*stimlist)[14]= 0  ; (*nlist)[14]= 375;
  (*stimlist)[15]=  1 ; (*nlist)[15]= 375;
  (*stimlist)[16]=   2; (*nlist)[16]= 375;
  (*stimlist)[17]= 0  ; (*nlist)[17]= 375;
  (*stimlist)[18]=   2; (*nlist)[18]= 375;
  (*stimlist)[19]= 0  ; (*nlist)[19]= 375;
  (*stimlist)[20]=  1 ; (*nlist)[20]= 375;
  (*stimlist)[21]= 0  ; (*nlist)[21]= 375;
  (*stimlist)[22]=  1 ; (*nlist)[22]= 375;
  (*stimlist)[23]=   2; (*nlist)[23]= 375;
  (*stimlist)[24]= 0  ; (*nlist)[24]= 375;
  (*stimlist)[25]=  1 ; (*nlist)[25]= 375;
  (*stimlist)[26]=   2; (*nlist)[26]= 375;

  (*stimlist)[27]= 0  ; (*nlist)[27]= 188; /*** 452l013.p08 starts here. ***/
  (*stimlist)[28]=   2; (*nlist)[28]= 188;
  (*stimlist)[29]=  1 ; (*nlist)[29]= 188;
  (*stimlist)[30]=   2; (*nlist)[30]= 188;
  (*stimlist)[31]=  1 ; (*nlist)[31]= 188;
  (*stimlist)[32]= 0  ; (*nlist)[32]= 188;
  (*stimlist)[33]=  1 ; (*nlist)[33]= 188;
  (*stimlist)[34]= 0  ; (*nlist)[34]= 188;
  (*stimlist)[35]=   2; (*nlist)[35]= 188;
  (*stimlist)[36]=   2; (*nlist)[36]= 188;
  (*stimlist)[37]=  1 ; (*nlist)[37]= 188;
  (*stimlist)[38]= 0  ; (*nlist)[38]= 188;

  *escflag = 1; /* Assume "esc" key pressed before this trial, 511 ==> -511 */
  *dstate = 12;
  *nss = 188;
  *rdwell = 16;
  *rdt = get_wn_actual_dt(100,0);
  *rlostlist = lostlist; *rnlost = 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                     GET_WY4ST_CRACK_PARAMS_452L019P10                     */
/*                                                                           */
/*  Cracked - all trials.                                                    */
/*                                                                           */
/*  Tried:                                                                   */
/*  - esc 1  19,19,19     - trial c/s 0 not broken                           */
/*  - esc 1  30,19,19     - all cracked.                                     */
/*                                                                           */
/*****************************************************************************/
void get_wy4st_crack_params_452l019p10(nstim,stimlist,nlist,escflag,dstate,
				       nss,rlostlist,rnlost,rdwell,rdt)
     int *nstim,**stimlist,**nlist,*escflag,*dstate,*nss,**rlostlist,*rnlost;
     int *rdwell;
     float *rdt;
{
  int *lostlist;

  *nstim = 15;
  lostlist = get_zero_iarray(*nstim);
  *stimlist = (int *)myalloc(*nstim*sizeof(int));
  *nlist = (int *)myalloc(*nstim*sizeof(int));
  /* Could be dwell was 10, or 16 for first 2, so 30 or 19 */
  (*stimlist)[0] = 0; (*nlist)[0] = 30; /* (3*100 - 1)/8  + 1 */
  (*stimlist)[1] = 0; (*nlist)[1] = 19;
  (*stimlist)[2] = 0; (*nlist)[2] = 19;

  (*stimlist)[3] = 0; (*nlist)[3] = 375;  /* (30*100 - 1)/8 + 1*/
  (*stimlist)[4] = 2; (*nlist)[4] = 375;
  (*stimlist)[5] = 1; (*nlist)[5] = 375;
  (*stimlist)[6] = 1; (*nlist)[6] = 375;
  (*stimlist)[7] = 0; (*nlist)[7] = 375;
  (*stimlist)[8] = 2; (*nlist)[8] = 375;
  (*stimlist)[9] = 0; (*nlist)[9] = 375;
  (*stimlist)[10]= 1; (*nlist)[10]= 375;
  (*stimlist)[11]= 2; (*nlist)[11]= 375;
  (*stimlist)[12]= 1; (*nlist)[12]= 375;
  (*stimlist)[13]= 2; (*nlist)[13]= 375;
  (*stimlist)[14]= 0; (*nlist)[14]= 375;

  *escflag = 1; /* Assume "esc" key pressed before this trial, 511 ==> -511 */
  *dstate = 12;
  *nss = 375;
  *rdwell = 8;
  *rdt = get_wn_actual_dt(100,0);
  *rlostlist = lostlist; *rnlost = 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                     GET_WY4ST_CRACK_PARAMS_452L021P08                     */
/*                                                                           */
/*  Cracked - all trials.                                                    */
/*                                                                           */
/*****************************************************************************/
void get_wy4st_crack_params_452l021p08(nstim,stimlist,nlist,escflag,dstate,
				       nss,rlostlist,rnlost,rdwell,rdt)
     int *nstim,**stimlist,**nlist,*escflag,*dstate,*nss,**rlostlist,*rnlost;
     int *rdwell;
     float *rdt;
{
  int *lostlist;

  *nstim = 14;
  lostlist = get_zero_iarray(*nstim);
  *stimlist = (int *)myalloc(*nstim*sizeof(int));
  *nlist = (int *)myalloc(*nstim*sizeof(int));
  /* Could be dwell was 10, or 16 for first 2, so 30 or 19 */
  (*stimlist)[0] = 0; (*nlist)[0] = 30; /* (3*100 - 1)/10 + 1 */
  (*stimlist)[1] = 0; (*nlist)[1] = 19; /* (3*100 - 1)/16 + 1 */
  (*stimlist)[2] = 2; (*nlist)[2] = 188; /*188 = (30*100 - 1)/16 + 1*/
  (*stimlist)[3] = 1; (*nlist)[3] = 188;
  (*stimlist)[4] = 0; (*nlist)[4] = 188;
  (*stimlist)[5] = 0; (*nlist)[5] = 188;
  (*stimlist)[6] = 2; (*nlist)[6] = 188;
  (*stimlist)[7] = 1; (*nlist)[7] = 188;
  (*stimlist)[8] = 0; (*nlist)[8] = 188;
  (*stimlist)[9] = 2; (*nlist)[9] = 188;
  (*stimlist)[10]= 1; (*nlist)[10]= 188;
  (*stimlist)[11]= 2; (*nlist)[11]= 188;
  (*stimlist)[12]= 1; (*nlist)[12]= 188;
  (*stimlist)[13]= 0; (*nlist)[13]= 188;

  *escflag = 1; /* Assume "esc" key pressed before this trial, 511 ==> -511 */
  *dstate = 12;
  *nss = 188;
  *rdwell = 16;
  *rdt = get_wn_actual_dt(100,0);
  *rlostlist = lostlist; *rnlost = 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                     GET_WY4ST_CRACK_PARAMS_452L022P09                     */
/*                                                                           */
/*  *** Probably no reboot, so must use state from 452l021.p08.              */
/*  NOT  YET.                                                                */
/*                                                                           */
/*****************************************************************************/
void get_wy4st_crack_params_452l022p09(nstim,stimlist,nlist,escflag,dstate,
				       nss,rlostlist,rnlost,rdwell,rdt)
     int *nstim,**stimlist,**nlist,*escflag,*dstate,*nss,**rlostlist,*rnlost;
     int *rdwell;
     float *rdt;
{
  int *lostlist;

  *nstim = 16;
  lostlist = get_zero_iarray(*nstim);
  *stimlist = (int *)myalloc(*nstim*sizeof(int));
  *nlist = (int *)myalloc(*nstim*sizeof(int));
  /* Could be dwell was 10, or 16 for first 2, so 30 or 19 */
  (*stimlist)[0] = 0; (*nlist)[0] = 19; /* (3*100 - 1)/10 + 1 */
  (*stimlist)[1] = 0; (*nlist)[1] = 19; /* (3*100 - 1)/16 + 1 */
  (*stimlist)[2] = 0; (*nlist)[2] = 19; /*188 = (30*100 - 1)/16 + 1*/
  (*stimlist)[3] = 0; (*nlist)[3] = 19;
  (*stimlist)[4] = 2; (*nlist)[4] = 188;
  (*stimlist)[5] = 0; (*nlist)[5] = 188;
  (*stimlist)[6] = 1; (*nlist)[6] = 188;
  (*stimlist)[7] = 0; (*nlist)[7] = 188;
  (*stimlist)[8] = 1; (*nlist)[8] = 188;
  (*stimlist)[9] = 2; (*nlist)[9] = 188;
  (*stimlist)[10]= 2; (*nlist)[10]= 188;
  (*stimlist)[11]= 1; (*nlist)[11]= 188;
  (*stimlist)[12]= 0; (*nlist)[12]= 188;
  (*stimlist)[13]= 1; (*nlist)[13]= 188;
  (*stimlist)[14]= 2; (*nlist)[14]= 188;
  (*stimlist)[15]= 0; (*nlist)[15]= 188;

  *escflag = 1; /* Assume "esc" key pressed before this trial, 511 ==> -511 */
  *dstate = 12;
  *nss = 188;
  *rdwell = 16;
  *rdt = get_wn_actual_dt(100,0);
  *rlostlist = lostlist; *rnlost = 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                     GET_WY4ST_CRACK_PARAMS_452L024P08                     */
/*                                                                           */
/*  Cracked - all trials.                                                    */
/*                                                                           */
/*****************************************************************************/
void get_wy4st_crack_params_452l024p08(nstim,stimlist,nlist,escflag,dstate,
				       nss,rlostlist,rnlost,rdwell,rdt)
     int *nstim,**stimlist,**nlist,*escflag,*dstate,*nss,**rlostlist,*rnlost;
     int *rdwell;
     float *rdt;
{
  int *lostlist;

  *nstim = 15;
  lostlist = get_zero_iarray(*nstim);
  *stimlist = (int *)myalloc(*nstim*sizeof(int));
  *nlist = (int *)myalloc(*nstim*sizeof(int));
  /* Could be dwell was 10, or 16 for first 3, so 30 or 19 */
  (*stimlist)[0] = 0; (*nlist)[0] = 30; /* (3*100 - 1)/10 + 1 */
  (*stimlist)[1] = 0; (*nlist)[1] = 19; /* (3*100 - 1)/16 + 1 */
  (*stimlist)[2] = 0; (*nlist)[2] = 19;
  (*stimlist)[3] = 0; (*nlist)[3] = 188; /*188 = (30*100 - 1)/16 + 1*/
  (*stimlist)[4] = 2; (*nlist)[4] = 188;
  (*stimlist)[5] = 1; (*nlist)[5] = 188;
  (*stimlist)[6] = 1; (*nlist)[6] = 188;
  (*stimlist)[7] = 2; (*nlist)[7] = 188;
  (*stimlist)[8] = 0; (*nlist)[8] = 188;
  (*stimlist)[9] = 2; (*nlist)[9] = 188;
  (*stimlist)[10]= 1; (*nlist)[10]= 188;
  (*stimlist)[11]= 0; (*nlist)[11]= 188;
  (*stimlist)[12]= 2; (*nlist)[12]= 188;
  (*stimlist)[13]= 1; (*nlist)[13]= 188;
  (*stimlist)[14]= 0; (*nlist)[14]= 188;

  *escflag = 1; /* Assume "esc" key pressed before this trial, 511 ==> -511 */
  *dstate = 12;
  *nss = 188;
  *rdwell = 16;
  *rdt = get_wn_actual_dt(100,0);
  *rlostlist = lostlist; *rnlost = 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                     GET_WY4ST_CRACK_PARAMS_452L025P09                     */
/*                                                                           */
/*  *** Tricky, unknown ordering because of "rerun".                         */
/*  NOT YET                                                                  */
/*                                                                           */
/*  Attempted:                                                               */
/*  - 30,30,30,300,300,300  (all permutations of 0,1,2 for the 3 300's)      */
/*  - 30,30,300,300,300,19  (all permutations of 0,1,2 for the 3 300's)      */
/*                                                                           */
/*****************************************************************************/
void get_wy4st_crack_params_452l025p09(nstim,stimlist,nlist,escflag,dstate,
				       nss,rlostlist,rnlost,rdwell,rdt)
     int *nstim,**stimlist,**nlist,*escflag,*dstate,*nss,**rlostlist,*rnlost;
     int *rdwell;
     float *rdt;
{
  int *lostlist;

  *nstim = 21;
  lostlist = get_zero_iarray(*nstim);
  lostlist[13] = lostlist[18] = lostlist[19] = 1;
  *stimlist = (int *)myalloc(*nstim*sizeof(int));
  *nlist = (int *)myalloc(*nstim*sizeof(int));
  /* Could be dwell was 10, or 16 for first 2, so 30 or 19 */
  (*stimlist)[0] = 0; (*nlist)[0] = 30; /* (3*100 - 1)/10 + 1 */
  (*stimlist)[1] = 0; (*nlist)[1] = 30; /* (3*100 - 1)/16 + 1 */

  /*** Some order of 0,1,2 in the next 3 spots.  I think dwell was 10. */
  (*stimlist)[2] = 2; (*nlist)[2] = 300; /*188 = (30*100 - 1)/16 + 1*/
  (*stimlist)[3] = 1; (*nlist)[3] = 300; /*  rerun loss - guess order  */
  (*stimlist)[4] = 0; (*nlist)[4] = 300; /*  rerun loss - guess order  */

  (*stimlist)[5] = 0; (*nlist)[5] = 19; /*  rerun loss - guess order  */

  (*stimlist)[6] = 0; (*nlist)[6] = 188; /*  0  */ /* Dwell now 16. */
  (*stimlist)[7] = 1; (*nlist)[7] = 188; /*  1  */
  (*stimlist)[8] = 2; (*nlist)[8] = 188; /*  2  */
  (*stimlist)[9] = 2; (*nlist)[9] = 188; /*  3  */
  (*stimlist)[10]= 1; (*nlist)[10]= 188; /*  4  */
  (*stimlist)[11]= 0; (*nlist)[11]= 188; /*  5  */
  (*stimlist)[12]= 0; (*nlist)[12]= 188; /*  6  */
  (*stimlist)[13]= 2; (*nlist)[13]= 188; /*    Lost on interrupt */
  (*stimlist)[14]= 2; (*nlist)[14]= 188; /*  7  */
  (*stimlist)[15]= 1; (*nlist)[15]= 188; /*  8  */
  (*stimlist)[16]= 2; (*nlist)[16]= 188; /*  9  */
  (*stimlist)[17]= 0; (*nlist)[17]= 188; /* 10  */
  (*stimlist)[18]= 0; (*nlist)[18]= 188; /*    Lost on interrupt */
  (*stimlist)[19]= 0; (*nlist)[19]= 188; /*    Lost on interrupt */
  (*stimlist)[20]= 1; (*nlist)[20]= 188; /* 11  */

  *escflag = 1; /* Assume "esc" key pressed before this trial, 511 ==> -511 */
  *dstate = 12;
  *nss = 188;
  *rdwell = 16;
  *rdt = get_wn_actual_dt(100,0);
  *rlostlist = lostlist; *rnlost = 3;
}
/**************************************-**************************************/
/*                                                                           */
/*                     GET_WY4ST_CRACK_PARAMS_452L026P09                     */
/*                                                                           */
/*  Cracked - all trials.                                                    */
/*                                                                           */
/*  Tried:                                                                   */
/*  - esc 1  19 -  Both of these worked, doesn't matter                      */
/*  - esc 1  30 -                                                            */
/*                                                                           */
/*****************************************************************************/
void get_wy4st_crack_params_452l026p09(nstim,stimlist,nlist,escflag,dstate,
				       nss,rlostlist,rnlost,rdwell,rdt)
     int *nstim,**stimlist,**nlist,*escflag,*dstate,*nss,**rlostlist,*rnlost;
     int *rdwell;
     float *rdt;
{
  int *lostlist;

  *nstim = 13;
  lostlist = get_zero_iarray(*nstim);
  *stimlist = (int *)myalloc(*nstim*sizeof(int));
  *nlist = (int *)myalloc(*nstim*sizeof(int));
  /* Could be dwell was 10, or 16 for first 3, so 30 or 19 */
  (*stimlist)[0] = 0; (*nlist)[0] = 19; /* (3*100 - 1)/10 + 1 */

  (*stimlist)[1] = 1; (*nlist)[1] = 188; /* (3*100 - 1)/16 + 1 */
  (*stimlist)[2] = 2; (*nlist)[2] = 188;
  (*stimlist)[3] = 0; (*nlist)[3] = 188; /*188 = (30*100 - 1)/16 + 1*/
  (*stimlist)[4] = 0; (*nlist)[4] = 188;
  (*stimlist)[5] = 2; (*nlist)[5] = 188;
  (*stimlist)[6] = 1; (*nlist)[6] = 188;
  (*stimlist)[7] = 0; (*nlist)[7] = 188;
  (*stimlist)[8] = 2; (*nlist)[8] = 188;
  (*stimlist)[9] = 1; (*nlist)[9] = 188;
  (*stimlist)[10]= 0; (*nlist)[10]= 188;
  (*stimlist)[11]= 2; (*nlist)[11]= 188;
  (*stimlist)[12]= 1; (*nlist)[12]= 188;

  *escflag = 1; /* Assume "esc" key pressed before this trial, 511 ==> -511 */
  *dstate = 12;
  *nss = 188;
  *rdwell = 16;
  *rdt = get_wn_actual_dt(100,0);
  *rlostlist = lostlist; *rnlost = 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                     GET_WY4ST_CRACK_PARAMS_452L034P10                     */
/*                                                                           */
/*  Cracked - all trials.                                                    */
/*                                                                           */
/*  Tried:                                                                   */
/*  - esc 1,2  19  19                                                        */
/*  - esc 1,2  30  19                                                        */
/*  - esc 1,2  30  30                                                        */
/*  - esc 1    30 188                                                        */
/*  - esc 1    10  19                                                        */
/*  - esc 1    10  10                                                        */
/*  - esc 1    30  10                                                        */
/*  - esc 1    19 188                                                        */
/*  - esc 1   188 188                                                        */
/*  - esc 1    30 300                                                        */
/*  - esc 1    38  19                                                        */
/*  - esc 1   188  19                                                        */
/*                                                                           */
/*****************************************************************************/
void get_wy4st_crack_params_452l034p10(nstim,stimlist,nlist,escflag,dstate,
				       nss,rlostlist,rnlost,rdwell,rdt)
     int *nstim,**stimlist,**nlist,*escflag,*dstate,*nss,**rlostlist,*rnlost;
     int *rdwell;
     float *rdt;
{
  int *lostlist;

  *nstim = 20;
  lostlist = get_zero_iarray(*nstim);
  *stimlist = (int *)myalloc(*nstim*sizeof(int));
  *nlist = (int *)myalloc(*nstim*sizeof(int));
  /* Could be dwell was 10, or 16 for first 3, so 30 or 19 */
  (*stimlist)[0] = 0; (*nlist)[0] = 300; /* (3*100 - 1)/10 + 1 */
  (*stimlist)[1] = 0; (*nlist)[1] = 19; /* (3*100 - 1)/16 + 1 */

  (*stimlist)[2] = 0  ; (*nlist)[2] = 188;
  (*stimlist)[3] =   2; (*nlist)[3] = 188; /*188 = (30*100 - 1)/16 + 1*/
  (*stimlist)[4] =  1 ; (*nlist)[4] = 188;
  (*stimlist)[5] = 0  ; (*nlist)[5] = 188;
  (*stimlist)[6] =   2; (*nlist)[6] = 188;
  (*stimlist)[7] =  1 ; (*nlist)[7] = 188;
  (*stimlist)[8] =  1 ; (*nlist)[8] = 188;
  (*stimlist)[9] = 0  ; (*nlist)[9] = 188;
  (*stimlist)[10]=   2; (*nlist)[10]= 188;
  (*stimlist)[11]= 0  ; (*nlist)[11]= 188;
  (*stimlist)[12]=   2; (*nlist)[12]= 188;
  (*stimlist)[13]=  1 ; (*nlist)[13]= 188;
  (*stimlist)[14]=   2; (*nlist)[14]= 188;
  (*stimlist)[15]= 0  ; (*nlist)[15]= 188;
  (*stimlist)[16]=  1 ; (*nlist)[16]= 188;
  (*stimlist)[17]= 0  ; (*nlist)[17]= 188;
  (*stimlist)[18]=  1 ; (*nlist)[18]= 188;
  (*stimlist)[19]=   2; (*nlist)[19]= 188;

  *escflag = 1; /* Assume "esc" key pressed before this trial, 511 ==> -511 */
  *dstate = 18;
  *nss = 188;
  *rdwell = 16;
  *rdt = get_wn_actual_dt(100,0);
  *rlostlist = lostlist; *rnlost = 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_SPECIAL_WY4ST_GROUP                          */
/*                                                                           */
/*****************************************************************************/
void get_special_wy4st_group(filename,nd,ndg,k,rstim,rns,rdwell,rdt)
     char filename[];
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k,***rstim,**rns,*rdwell;
     float *rdt;
{
  int i,j,m;
  int *pseed,s511,s0,sn511,*stimlist,*nlist,nstim,firstflag;
  int dstate,**ss,*sstype,nss,escflag,**gstim,*gns,*lostlist,nlost,nskip;
  float frv;

  printf("  GET_SPECIAL_WY4ST_GROUP\n");

  if (strcmp(filename,"452l004.p17")==0)
    get_wy4st_crack_params_452l004p17(&nstim,&stimlist,&nlist,&escflag,&dstate,
				      &nss,&lostlist,&nlost,rdwell,rdt);
  else if (strcmp(filename,"452l019.p10")==0)
    get_wy4st_crack_params_452l019p10(&nstim,&stimlist,&nlist,&escflag,&dstate,
				      &nss,&lostlist,&nlost,rdwell,rdt);
  else if (strcmp(filename,"452l021.p08")==0)
    get_wy4st_crack_params_452l021p08(&nstim,&stimlist,&nlist,&escflag,&dstate,
				      &nss,&lostlist,&nlost,rdwell,rdt);
  else if (strcmp(filename,"452l013.p07")==0)
    get_wy4st_crack_params_452l013p07(&nstim,&stimlist,&nlist,&escflag,&dstate,
				      &nss,&lostlist,&nlost,rdwell,rdt);
  else if (strcmp(filename,"452l013.p08")==0)
    get_wy4st_crack_params_452l013p08(&nstim,&stimlist,&nlist,&escflag,&dstate,
				      &nss,&lostlist,&nlost,rdwell,rdt);
  else if (strcmp(filename,"452l024.p08")==0)
    get_wy4st_crack_params_452l024p08(&nstim,&stimlist,&nlist,&escflag,&dstate,
				      &nss,&lostlist,&nlost,rdwell,rdt);
  else if (strcmp(filename,"452l025.p09")==0)
    get_wy4st_crack_params_452l025p09(&nstim,&stimlist,&nlist,&escflag,&dstate,
				      &nss,&lostlist,&nlost,rdwell,rdt);
  else if (strcmp(filename,"452l026.p09")==0)
    get_wy4st_crack_params_452l026p09(&nstim,&stimlist,&nlist,&escflag,&dstate,
				      &nss,&lostlist,&nlost,rdwell,rdt);
  else if (strcmp(filename,"452l034.p10")==0)
    get_wy4st_crack_params_452l034p10(&nstim,&stimlist,&nlist,&escflag,&dstate,
				      &nss,&lostlist,&nlost,rdwell,rdt);
  else
    exit_error("GET_SPECIAL_WY4ST_GROUP","Unknown filename");

  ss = get_2d_iarray(dstate,nss);
  sstype = (int *)myalloc(dstate*sizeof(int));

  sn511 = 8; /* God-given value that resides at location -511 after reboot. */
  s511 = 0;  /* Value at location 511 after ldp. */
  s0 = 0;    /* Value at location 0 after ldp. */
  nskip = nstim - (dstate + nlost); /* No. of stimuli before final "run" */
  m = 0; /* Count the number of saved stimuli "ss" */
  firstflag = 0;
  for(i=0;i<nstim;i++){
    if (i==0){ /* Set the seed. */
      firstflag = 1;  /* Stays 1 until the state ID changes */
    }else if (firstflag == 1){
      if (stimlist[i] != stimlist[i-1])
	firstflag = 0;
      else
	if (i >= escflag)
	  firstflag = 0;
    }
    
    if (firstflag == 1)
      pseed = &s511;  /* First trial always points to location 511. */
    else
      if ((stimlist[i] == 0)||(stimlist[i] == 2))
	pseed = &sn511;
      else
	pseed = &s0;

    if (((strcmp(filename,"452l013.p07")==0)||
	 (strcmp(filename,"452l013.p08")==0)) && (i==7)) /*** Special ***/
      s0 = 0;
    
    /* Run the random generator. */
    for(j=0;j<nlist[i];j++){
      frv = noise_util_ran2(pseed);
      if ((i >= nskip)&&(lostlist[i] == 0)){
	if (frv < 0.25)
	  ss[m][j] = 3;           /* state 3 = PP (center,surround) */
	else if (frv < 0.50)
	  ss[m][j] = 2;           /* state 2 = PA */
	else if (frv < 0.75)
	  ss[m][j] = 1;           /* state 1 = AP */
	else
	  ss[m][j] = 0;           /* state 0 = AA */
      }
    }
    if ((i >= nskip)&&(lostlist[i] == 0)){ /* Count number of kept stimuli */
      sstype[m] = stimlist[i];
      m += 1;
    }
  }

  /*** Adjust state values in cases where there is no center or no surr. ***/
  /*** Report center state as 0,1 ***/
  for(i=0;i<dstate;i++){
    /*printf("sstype[%d] = %d\n",i,sstype[i]);*/
    if (sstype[i] == 1)
      for(j=0;j<nss;j++)
	if (ss[i][j] >= 2)
	  ss[i][j] = 1;
	else
	  ss[i][j] = 0;
    else if (sstype[i] == 2)
      for(j=0;j<nss;j++){
	if (ss[i][j] == 3)
	  ss[i][j] = 1;
	if (ss[i][j] == 2)
	  ss[i][j] = 0;
      }
  }

  gstim = (int **)myalloc(ndg->cnt[k]*sizeof(int *));
  gns = (int *)myalloc(ndg->cnt[k]*sizeof(int));
  for(i=0;i<ndg->cnt[k];i++){
    j = ndg->tnum[k][i];
    gstim[i] = copy_iarray(ss[j],nss);
    gns[i] = nss;
  }

  myfree(stimlist); myfree(nlist); free_2d_iarray(ss,dstate);
  myfree(lostlist); myfree(sstype);
  *rstim = gstim; *rns = gns;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_TWO_NOISE_STIM_GROUP                         */
/*                                                                           */
/*  This handles pfiles for "wy2pat.out" that use two stimuli.               */
/*                                                                           */
/*****************************************************************************/
void get_two_noise_stim_group(nd,ndg,k,awake_flag,resolution,fft_flag,
			      norm_flag,rstim1,rstim2,rns,rdt)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k,awake_flag;
     float resolution;
     int norm_flag,fft_flag;
     float ***rstim1,***rstim2;
     int **rns;
     float *rdt;
{
  int i;
  int flag,*ns,tn,n;
  float *tx,*ty1,*ty2,**stim1,**stim2,res_in,dt,oldmean;
  char *pfile;

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0)
    exit_error("GET_TWO_NOISE_STIM_GROUP","Const param pfile not found");
  if ((strcmp(pfile, "wy2pcs9.p")==0)){
    n = ndg->cnt[k];
    stim1 = (float **)myalloc(n*sizeof(float *));
    stim2 = (float **)myalloc(n*sizeof(float *));
    ns = (int *)myalloc(n*sizeof(int));
    for(i=0;i<n;i++){
      get_wy2pat_stim_for_trial_02(nd,ndg->tnum[k][i],awake_flag,pfile,&tx,
				   &ty1,&ty2,&tn,&dt);
      /*dt = (tx[tn-1] - tx[0])/(float)(tn-1);*/
      res_in = 1000.0/dt;
      
      if (fft_flag == 0){
	get_stim_at_resolution_simple(ty1,tn,res_in,resolution,&stim1[i],
				      &ns[i]);
	get_stim_at_resolution_simple(ty2,tn,res_in,resolution,&stim2[i],
				      &ns[i]);
      }else if (fft_flag == 1){
	get_stim_at_resolution(ty1,tn,res_in,resolution,&stim1[i],&ns[i]);
	get_stim_at_resolution(ty2,tn,res_in,resolution,&stim2[i],&ns[i]);
      }else if (fft_flag == 2){
	get_stim_at_resolution_impulse(ty1,tn,res_in,resolution,&stim1[i],
				       &ns[i]);
	get_stim_at_resolution_impulse(ty2,tn,res_in,resolution,&stim2[i],
				       &ns[i]);
      }else if (fft_flag == 3){
	get_stim_at_resolution_piecewise(ty1,tn,res_in,resolution,&stim1[i],
				       &ns[i]);
	get_stim_at_resolution_piecewise(ty2,tn,res_in,resolution,&stim2[i],
				       &ns[i]);
      }else
	exit_error("GET_TWO_NOISE_STIM_GROUP","Unknown fft_flag value");

      if (norm_flag > 0){
	subtract_mean_farray(stim1[i],ns[i],&oldmean);
	norm_variance_quiet_farray(stim1[i],ns[i]);
	subtract_mean_farray(stim2[i],ns[i],&oldmean);
	norm_variance_quiet_farray(stim2[i],ns[i]);
      }

      myfree(tx); myfree(ty1); myfree(ty2);
    }
    *rstim1 = stim1; *rstim2 = stim2; *rns = ns; *rdt = dt;
    myfree(pfile);
  }else{
    printf("pfile = %s\n",pfile);
    exit_error("GET_TWO_NOISE_STIM_GROUP","Unrecognized pfile");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_NEW_NOISE_STIM_GROUP                         */
/*                                                                           */
/*  This handles a new set of pfiles, those that use "wy2pat.out".           */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Set 'fft_flag' >= 10 to convert jumps to position.                     */
/*                                                                           */
/*****************************************************************************/
void get_new_noise_stim_group(nd,ndg,k,awake_flag,resolution,fft_flag,
			      norm_flag,rstim,rns,rdt)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k,awake_flag;
     float resolution;
     int norm_flag,fft_flag;
     float ***rstim;
     int **rns;
     float *rdt;
{
  int i,j;
  int flag,*ns,tn,n,maxj,phase,pos_flag;
  float *tx,*ty,**stim,res_in,dt,oldmean;
  char *pfile,*filename;

  if (fft_flag >= 10){ /*** Position values. ***/
    fft_flag -= 10;
    pos_flag = 1;
  }else
    pos_flag = 0;

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0)
    exit_error("GET_NEW_NOISE_STIM_GROUP","Const param pfile not found");
  if ((strcmp(pfile, "wypnori.p")==0)||(strcmp(pfile, "wypnpha.p")==0)||
      (strcmp(pfile,"wypnori3.p")==0)||(strcmp(pfile,"wypnpha3.p")==0)||
      (strcmp(pfile,   "wy4st.p")==0)||(strcmp(pfile,   "wyone.p")==0)||
      (strcmp(pfile,  "wyonem.p")==0)||
      (strcmp(pfile, "wyone23.p")==0)||(strcmp(pfile,"wypnoril.p")==0)||
      (strcmp(pfile,  "wy4stg.p")==0)||(strcmp(pfile,"wypnphal.p")==0)||
      (strcmp(pfile, "wy2flat.p")==0)||(strcmp(pfile,    "x4st.p")==0)||
      (strcmp(pfile,   "x4stp.p")==0)||(strcmp(pfile,   "x4std.p")==0)||
      (strcmp(pfile,  "x4stp2.p")==0)||
      (strcmp(pfile,  "wy4sum.p")==0)||(strcmp(pfile,"wyphstat.p")==0)||
      (strcmp(pfile,   "wycon.p")==0)||(strcmp(pfile, "wyphasm.p")==0)||
      (strcmp(pfile,"wyphasms.p")==0)||
      (strcmp(pfile,"wyphasm8.p")==0)||(strcmp(pfile, "wyph3.p")==0)||
      (strcmp(pfile,  "wymatm.p")==0)||(strcmp(pfile, "wyphas8.p")==0)||
      (strcmp(pfile,  "wymat160.p")==0)||(strcmp(pfile, "wyall.p")==0)||
      (strcmp(pfile,  "wymtf9.p")==0)||(strcmp(pfile,  "wymsf9.p")==0)||
      (strcmp(pfile, "wymtf7t.p")==0)||
      (strcmp(pfile,  "wymsf6.p")==0)||(strcmp(pfile, "wymtf9p.p")==0)||
      (strcmp(pfile,    "wydt.p")==0)||
      (strcmp(pfile, "wymsf11.p")==0)||(strcmp(pfile,  "wymtf5.p")==0)||
      (strcmp(pfile, "wymtf5p.p")==0)||(strcmp(pfile,"wymtf4ph.p")==0)||
      (strcmp(pfile,   "wmcon.p")==0)||(strcmp(pfile,   "wmco1.p")==0)||
      (strcmp(pfile,   "wmco4.p")==0)||(strcmp(pfile,   "wmlu4.p")==0)||
      (strcmp(pfile,   "wmlum.p")==0)||(strcmp(pfile,   "wmlu1.p")==0)||
      (strcmp(pfile,  "wyamc7.p")==0)||(strcmp(pfile,"wyamcsf5.p")==0)||
      (strcmp(pfile,   "wyamc.p")==0)||(strcmp(pfile,"wyamcpha.p")==0)||
      (strcmp(pfile,  "wyamc5.p")==0)||
      (strcmp(pfile,"MODEL_wyamc.p")==0)||(strcmp(pfile,"MODEL_new.p")==0)||
      (strcmp(pfile,"wmcon160.p")==0)||(strcmp(pfile,"wmco1160.p")==0)||
      (strcmp(pfile,"wmlum160.p")==0)||(strcmp(pfile,"wmlu1160.p")==0)){
    /*** Handle any special case processing. ***/
    flag = ndata_get_const_param(nd,"filename",&filename);
    if (flag > 0){
      if (strcmp(filename,"452l032.p07")==0){
	get_special_stim_group_452l032p07(nd,ndg,k,resolution,fft_flag,
					  norm_flag,rstim,rns,rdt);
	return;
      }
    }

    dt = -1.0; /* Avoid -Wall warning */
    
    n = ndg->cnt[k];
    stim = (float **)myalloc(n*sizeof(float *));
    ns = (int *)myalloc(n*sizeof(int));
    for(i=0;i<n;i++){
      j = ndg->tnum[k][i];
      if ((strcmp(pfile,"wyamc5.p")==0)||(strcmp(pfile,"wyamc7.p")==0)||
	  (strcmp(pfile,"MODEL_wyamc.p")==0)){
	/*printf(" **** USING  pos_flag = %d\n",pos_flag);*/
	get_wy2pat_stim_for_trial_03(nd,j,awake_flag,pos_flag,pfile,&tx,&ty,
				     &tn);
      }else if (strcmp(pfile,"wy2flat2.p")==0){
	get_wy2pat_stim_for_trial_04(nd,j,awake_flag,pfile,&tx,&ty,&tn);
      }else{
	get_wy2pat_stim_for_trial_01(nd,j,awake_flag,pfile,&tx,&ty,&tn);

	if (pos_flag){ /*** Position values. ***/
	  flag = ndata_get_vc_param_int(nd,j,"phase",&phase); /* Def. zero? */
	  if (!flag) exit_error("***","Param phase not found");
	  flag = ndata_get_vc_param_int(nd,j,"maxj",&maxj);
	  if (!flag) exit_error("***","Param maxj not found");
	  convert_jump_to_position(ty,tn,phase,maxj,1024);
	}
      }

      dt = (tx[tn-1] - tx[0])/(float)(tn-1);
      res_in = 1000.0/dt;

      if (fft_flag == 0)
	get_stim_at_resolution_simple(ty,tn,res_in,resolution,&stim[i],&ns[i]);
      else if (fft_flag == 1)
	get_stim_at_resolution(ty,tn,res_in,resolution,&stim[i],&ns[i]);
      else if (fft_flag == 2)
	get_stim_at_resolution_impulse(ty,tn,res_in,resolution,&stim[i],
				       &ns[i]);
      else if (fft_flag == 3)
	get_stim_at_resolution_piecewise(ty,tn,res_in,resolution,&stim[i],
					 &ns[i]);
      else
	exit_error("GET_NEW_NOISE_STIM_GROUP","Unknown fft_flag value");

      if (norm_flag > 0){
	subtract_mean_farray(stim[i],ns[i],&oldmean);
	norm_variance_farray(stim[i],ns[i]);
      }

      myfree(tx); myfree(ty);
    }
    *rstim = stim; *rns = ns; *rdt = dt;
    myfree(pfile);
  }else{
    printf("pfile = %s\n",pfile);
    exit_error("GET_NEW_NOISE_STIM_GROUP","Unrecognized pfile");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_NOISE_STIM_GROUP                            */
/*                                                                           */
/*  This handles pfiles that use "noise.out"                                 */
/*  Called by 'STA_NDA'.                                                     */
/*                                                                           */
/*****************************************************************************/
void get_noise_stim_group(nd,ndg,k,awake_flag,resolution,fft_flag,norm_flag,
			  rstim,rns,rdt)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k,awake_flag;
     float resolution;
     int norm_flag,fft_flag;
     float ***rstim;
     int **rns;
     float *rdt;
{
  int i,j;
  int flag,*ns,tn,n,angle,nsl,ti;
  /*int *ti,j;*/
  float *tx,*ty,**stim,res_in,dt;
  char *pfile,*stimtype,**slist;

  dt = -1.0; /* Avoid -Wall warning */

  if (strcmp(nd->class,"labradoodle_01")==0){
    // Get the 'stim_type'
    flag = ndata_get_const_param(nd,"stim_type",&stimtype);
    if (flag == 0)
      exit_error("GET_NOISE_STIM_GROUP","Const param 'stim_type' not found");
    get_items_from_string(stimtype,&slist,&nsl);
    myfree(stimtype);
    stimtype = strdup(slist[0]);
    free_2d_carray(slist,nsl);

    if (strcmp(stimtype,"ranstep")==0){
      //int jz,jn,seed,dti;
      //float sf,jsz,stn,fps,diff,*newstim;

      float **tstim;
      int *tns,iflag;

      iflag = 0; // Use values, not indices, for stimulus
      nu_ranstep_get_stim_group(nd,ndg,k,iflag,norm_flag,&tstim,&tns,&dt);

      n = ndg->cnt[k];
      stim = (float **)myalloc(n*sizeof(float *));
      ns = (int *)myalloc(n*sizeof(int));
      res_in = 1000.0/dt;

      for(i=0;i<n;i++)
	get_stim_at_resolution_simple(tstim[i],tns[i],res_in,resolution,
				      &(stim[i]),&(ns[i]));

      *rstim = stim; *rns = ns; *rdt = dt;
    }else{
      printf("stimtype = %s\n",stimtype);
      exit_error("GET_NOISE_STIM_GROUP","Unknown stimtype");
    }
  }else{
    flag = ndata_get_const_param(nd,"pfile",&pfile);
    if (flag == 0)
      exit_error("GET_NOISE_STIM_GROUP","Const param pfile not found");
    
    if ((strcmp(pfile,"wnall100.p")==0)||(strcmp(pfile,"wncon100.p")==0)||
	(strcmp(pfile,"wnmat050.p")==0)||(strcmp(pfile,"wnmat053.p")==0)||
	(strcmp(pfile,"wnmat080.p")==0)||(strcmp(pfile,"wnmat100.p")==0)||
	(strcmp(pfile,"wnmdt100.p")==0)||(strcmp(pfile,"wnone053.p")==0)||
	(strcmp(pfile,"wnone080.p")==0)||(strcmp(pfile,"wnone100.p")==0)||
	(strcmp(pfile, "wnposdt.p")==0)||(strcmp(pfile,"wnsiz100.p")==0)||
	(strcmp(pfile,   "MODEL.p")==0)||(strcmp(pfile,  "wnjump.p")==0)||
	(strcmp(pfile,   "vbin1.p")==0)||(strcmp(pfile, "vbinmat.p")==0)){
      n = ndg->cnt[k];
      stim = (float **)myalloc(n*sizeof(float *));
      ns = (int *)myalloc(n*sizeof(int));
      for(i=0;i<n;i++){
	get_wn_stim_for_trial(nd,ndg->tnum[k][i],awake_flag,pfile,&tx,&ty,&tn);
	/*** Next line is for printing out stimulus. ***/
	/*append_farray_plot("zzz.stim",ndg->name[k],ty,tn,1);*/
	
	if (strcmp(pfile,"vbin1.p")==0){
	  flag = ndata_get_var_value_trial_int(nd,ndg->tnum[k][i],"angle"
					       ,&angle);
	  if ((nd->t[ndg->tnum[k][i]].tcode % 2) == 1) /*** NULL angle ***/
	    multiply_farray(ty,tn,-1.0);
	}
	
	/***  USED TO PRINT OUT STIMULUS DATA FOR WEB PAGE
	      ti = f2iarray(ty,tn);
	      printf("%d\n",tn);
	      for(j=0;j<tn;j++)
	      printf("%d ",ti[j]);
	      printf("\n");
	      myfree(ti);***/
	
	/*write_farray_xy_plot("zz.xy.stim.pl",tx,ty,tn);*/
	/*exit(0);*/
	
	/*** write_farray_xy_plot("zz.xy.stim.pl",tx,ty,tn);
	     exit(0);*/
	
	dt = (tx[tn-1] - tx[0])/(float)(tn-1);
	res_in = 1000.0/dt;
	
	if (fft_flag == 0)
	  get_stim_at_resolution_simple(ty,tn,res_in,resolution,&stim[i],
					&ns[i]);
	else if (fft_flag == 1)
	  get_stim_at_resolution(ty,tn,res_in,resolution,&stim[i],&ns[i]);
	else if (fft_flag == 2)
	  get_stim_at_resolution_impulse(ty,tn,res_in,resolution,&stim[i],
					 &ns[i]);
	else if (fft_flag == 3)
	  get_stim_at_resolution_piecewise(ty,tn,res_in,resolution,&stim[i],
					   &ns[i]);
	else
	  exit_error("GET_NOISE_STIM_GROUP","Unknown fft_flag value");
	
	if (norm_flag > 0){
	  norm_variance_farray(stim[i],ns[i]);
	}
	
	myfree(tx); myfree(ty);
      }
      *rstim = stim; *rns = ns; *rdt = dt;
      myfree(pfile);
    }else // Maybe this is a new pfile
      get_new_noise_stim_group(nd,ndg,k,awake_flag,resolution,fft_flag,
			       norm_flag,rstim,rns,rdt);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_MODEL_GRID_STIM_GROUP                        */
/*                                                                           */
/*  Return the stimulus for 'wm' stimulus type 'noise'.                      */
/*                                                                           */
/*  Return:                                                                  */
/*    stim[trial][time][gridbox]  - stimulus values, at spike resolution     */
/*     nst[trial]                 - time duration, at spike res (e.g. msec)  */
/*      ns                        - number of stimuli (i.e. grid boxes)      */
/*      dt                        - dwell time                               */
/*                                                                           */
/*****************************************************************************/
void get_model_grid_stim_group(nd,ndg,k,sampling,rgnx,rgny,rstim,rnst,rns,rdt)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k;
     float sampling;    // samples/s for the expanded stimulus
     int *rgnx,*rgny;
     float ****rstim;
     int **rnst,*rns;
     float *rdt;
{
  int i,j;
  int fi,ti,xi,xj,sk,*nst,ns,n,gnx,gny,dwell,ri,seed,dwellms;
  float ***stim,fpow,***d3ran;
  char prname[SLEN],*stim_type,*set_name,*image_dir;

  int moo_tn,tsamp,*seedlist,nframes,si;
  float moo_tscale,stim_samp,st0,stn;

  strcpy(prname,"GET_MODEL_GRID_STIM_GROUP");

  moo_tn     = ndata_get_const_param_int(nd,"MOO_tn");
  moo_tscale = ndata_get_const_param_float(nd,"MOO_tscale");
  stim_samp  = ndata_get_const_param_float(nd,"stim_samp");
  stim_type  = ndata_get_const_param_char(nd,"stim_type");
  tsamp = stm_get_tsamp(moo_tscale,stim_samp);

  gnx = gny = ns = dwell = 0; // Avoid -Wall warning

  n = ndg->cnt[k];
  stim = (float ***)myalloc(n*sizeof(float **));
  nst = (int *)myalloc(n*sizeof(int));

  for(i=0;i<n;i++){      // For each trial in this group
    j = ndg->tnum[k][i]; //   get trial number

    st0   = ndata_get_vc_param_float_or_exit(nd,j,"st0",prname);
    stn   = ndata_get_vc_param_float_or_exit(nd,j,"stn",prname);
    gnx   = ndata_get_vc_param_int_or_exit(nd,j,"grid_w",prname);
    gny   = ndata_get_vc_param_int_or_exit(nd,j,"grid_h",prname);
    dwell = ndata_get_vc_param_int_or_exit(nd,j,"dwell",prname);
    seed  = ndata_get_vc_param_int_or_exit(nd,j,"seed",prname);
    fpow  =          ndata_get_vc_flt_dflt(nd,j,"power_law_2d",0.0);

    // Unique frame count
    nframes = stm_get_n(moo_tn,moo_tscale,st0,stn,dwell,tsamp);
    seedlist = get_seeds(seed,100000,nframes);  // One seed per uniqe frame

    // Create entire 3D random value array now, so that it can be
    // Transformed, if needed.

    d3ran = (float ***)myalloc(nframes*sizeof(float **));

    if (strcmp(stim_type,"noise")==0){
      for(fi=0;fi<nframes;fi++){
	d3ran[fi] = myrand_get_std_bin_float_2d(gnx,gny,seedlist[fi]);
	if (fpow != 0.0){
	  power_law_2d_transform(d3ran[fi],gnx,gny,fpow,"1 max 0.5 mid");
	}
      }
    }else if (strcmp(stim_type,"image_set")==0){
      set_name = ndata_get_const_param_char(nd,"set_name");
      image_dir = ndata_get_const_param_char(nd,"set_dir");
      if (strcmp(set_name,"vanhateren")==0){
	for(fi=0;fi<nframes;fi++){
	  //if (fi == 0)
	  //printf("WY Seed = %d\n",seedlist[fi]);
	  d3ran[fi] = stm_vanhat_get_patch(seedlist[fi],gnx,gny,image_dir,
					   "-1 to 1");
	  //printf("fi = %d\n",fi);
	}
      }else{
	exit_error("GET_MODEL_GRID_STIM_GROUP","Unknown 'set_name'");
      }
      myfree(image_dir);
      myfree(set_name);
    }else{
      exit_error("GET_MODEL_GRID_STIM_GROUP","Unknown 'stim_type'");
    }

    // Duration of one unique pattern presentation (spike time units)
    dwellms = my_rint(tsamp * dwell * moo_tscale *    // seconds *
		      sampling);       //   samp/sec
    //nd->t[j].r[ri].sampling);       //   samp/sec

    // Stimulus duration (in spike sampling units)
    nst[i] = nframes * dwellms;  // Unique frames * sampling units per frame

    //printf("    %dx%d dwell %d ms  period %d seed %d\n",gnx,gny,dwellms,
    //nst[i],seed);


    if (i==0)
      ns = gnx * gny;
    else
      if (ns != (gnx*gny))
	exit_error(prname,"Number of grid boxes varies within group");

    // Fill 'stim' for this trial
    si = 0;  // Seed index
    stim[i] = get_2d_farray(nst[i],ns);
    for(ti=0;ti<nst[i];ti++){ // For each time unit

      sk = 0; // grid box count

      if (ti%dwellms == 0){ // Pick new random values

	if (ti == -1){
	  printf("%f %f %f %f\n",d3ran[0][0][0],d3ran[1][0][0],d3ran[2][0][0],
		 d3ran[3][0][0]);
	}
	
	for(xi=0;xi<gnx;xi++){
	  for(xj=0;xj<gny;xj++){
	    stim[i][ti][sk] = d3ran[si][xi][xj];
	    sk += 1;

	  }
	}
	si += 1;
      }else // Fill in with the value at the previous time step
	for(xi=0;xi<gnx;xi++)
	  for(xj=0;xj<gny;xj++){
	    stim[i][ti][sk] = stim[i][ti-1][sk];
	    sk += 1;
	  }
    }
    myfree(seedlist);
    free_3d_farray(d3ran,nframes,gnx,gny);
  }

  myfree(stim_type);

  *rgnx = gnx; *rgny = gny;
  *rstim = stim; *rnst = nst; *rns = ns; *rdt = (float)dwell;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_LABR_GRID_STIM_GROUP                        */
/*                                                                           */
/*  Return the stimulus for labradoodle 'grid' stimulus.                     */
/*                                                                           */
/*  Return:                                                                  */
/*    stim[trial][time][gridbox]  - stimulus values                          */
/*     nst[trial]                 - time duration                            */
/*      ns                        - number of stimuli (i.e. grid boxes)      */
/*      dt                        - dwell time                               */
/*                                                                           */
/*****************************************************************************/
void get_labr_grid_stim_group(nd,ndg,k,recname,rgnx,rgny,rstim,rnst,rns,rdt)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k;
     char recname[]; // This is needed to get period, this is a hack
     int *rgnx,*rgny;
     float ****rstim;
     int **rnst,*rns;
     float *rdt;
{
  int i,j;
  int ti,tk,xi,xj,sk,*nst,ns,n,gnx,gny,dwell,seed,tsi,*tseq,tn;
  float ***stim,fps,stn,***newstim;
  char prname[SLEN];

  strcpy(prname,"GET_LABR_GRID_STIM_GROUP");

  gnx = gny = ns = dwell = 0; /* Avoid -Wall warning */

  n = ndg->cnt[k];
  stim = (float ***)myalloc(n*sizeof(float **));
  nst = (int *)myalloc(n*sizeof(int));

  for(i=0;i<n;i++){      // For each trial in this group
    j = ndg->tnum[k][i]; //   get trial number

    gnx   = ndata_get_vc_param_int_or_exit(nd,j,"grid_xn",prname);
    gny   =          ndata_get_vc_int_dflt(nd,j,"grid_yn",gnx);
    dwell = ndata_get_vc_param_int_or_exit(nd,j,"dt",prname);
    seed  = ndata_get_vc_param_int_or_exit(nd,j,"seed",prname);
    fps   =          ndata_get_vc_flt_exit(nd,j,"labrcon_fps",prname);
    stn   =          ndata_get_vc_flt_exit(nd,j,"stn",prname);
    printf("    Frames/s = %f\n",fps);
    printf("    Duration = %f sec\n",stn);

    // Length of random sequence
    tn = (int)(fps * stn) / dwell;

    // Get number of grid boxes 'ns'
    if (i==0)
      ns = gnx * gny;
    else
      if (ns != (gnx*gny))
	exit_error(prname,"Number of grid boxes varies within group");

    // Get random sequence (same as in ss_grid_test/prep)
    // In 'tseq' run through each grid box, then move to next time step
    tseq = myrand_get_std_bin_seq(tn*ns,seed);


    // Compute stimulus duration in msec (at resolution of spike times)
    nst[i] = (int)(1000.0 * (float)(tn*dwell) / fps);

    printf("    nst[i] = %d\n",nst[i]);
    printf("    %dx%d dwell %d period %d seed %d\n",gnx,gny,dwell,nst[i],seed);

    // Fill 'stim' for this trial
    stim[i] = get_2d_farray(nst[i],ns);

    for(ti=0;ti<nst[i];ti++){ // For each time unit

      // Determine time index 'tk' for tseq from msec index 'ti'
      tk = (int)((double)ti * (double)fps / 1000.0 / (double)dwell);
      if (tk >= tn){
	printf(" ti=%d (of %d) tk = %d  >=  tn %d\n",ti,nst[i],tk,tn);
	exit_error("HERE WYETH","tk is too large");
      }
      
      sk = 0; // grid box count
      for(xi=0;xi<gnx;xi++){
	for(xj=0;xj<gny;xj++){
	  
	  tsi = tk*ns + sk;  // Index into 'tseq', 'tk' time, 'sk' gridbox

	  if (tseq[tsi] == 0)
	    stim[i][ti][sk] = -1.0;
	  else
	    stim[i][ti][sk] =  1.0;
	  
	  sk += 1;
	}
      }
    }
  }

  /*** NEVER WORKED, in that the STA was all -1,
  //   - I tried making this in lab on the fly...
  //
  //   WYETH - big hack
  //
  printf("****** WYETH WARNING:  Hack to change to change\n");

  newstim = get_3d_farray(n,nst[0],gnx*gny);
  for(i=0;i<n;i++){      // For each trial in this group
    
    sk = 0; // Grid box count
    for(xi=0;xi<gnx;xi++){
      for(xj=0;xj<gny;xj++){
	
	
	newstim[i][0][sk] = 1.0;  // First box is a change
	for(ti=1;ti<nst[i];ti++){ // For each time unit
	  if (stim[i][ti][sk] == stim[i][ti-1][sk])
	    newstim[i][ti][sk] = -1.0;  // No change
	  else
	    newstim[i][ti][sk] =  1.0;  // Change
	}
	
	sk += 1;
      }
    }
  }
  stim = newstim;
  ***/

  *rgnx = gnx; *rgny = gny;
  *rstim = stim; *rnst = nst; *rns = ns; *rdt = 1000.0/fps * (float)dwell;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_NDATA_STIM_GROUP                          */
/*                                                                           */
/*  Get the stimulus stored on 'stimchan' for the 'k'th group.  The 'rstim'  */
/*  should be freed by the caller.                                           */
/*                                                                           */
/*  Notes:                                                                   */
/*  - binflag = 0-default, 1-Binary, 3-Ternary, 4-Gaussian                   */
/*                                                                           */
/*****************************************************************************/
void get_ndata_stim_group(nd,ndg,k,stimchan,binflag,rstim,rns,rdt,rt0)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k;
     char stimchan[];
     int binflag;
     float ***rstim;   /* stim[i][ns[i]] stimulus values for ith trial */
     int **rns;        /* ns[i] number of stimulus values for ith trial */
     float *rdt;       /* time between stimulus values (msec) */
     int *rt0;         /* stimulus onset (msec) */
{
  int i,j;
  int *ns,n,t0;
  float **stim,dt,sampling,*px; /*mean,sdev*/
  struct ndtrial_struct *t;

  t0 = 0;     /* Avoid -Wall warning */
  dt = -1.0;

  n = ndg->cnt[k];
  stim = (float **)myalloc(n*sizeof(float *));
  ns = (int *)myalloc(n*sizeof(int));

  for(i=0;i<n;i++){      /* For each trial in this group */
    j = ndg->tnum[k][i]; /*   get trial number */
    t = &(nd->t[j]);     /*   get pointer to trial */
    px = ndata_trial_getp_adata_chan(t,stimchan,&ns[i]); /* get stim. */
    stim[i] = copy_farray(px,ns[i]);

    sampling = ndata_trial_get_sampling_chan(t,stimchan);
    dt = 1000.0/sampling; /* return 'dt' in msec */

    t0 = my_rint(1000.0 * ndata_trial_get_start_sec_chan(t,stimchan));

    if (binflag == 1){
      /*mean_sdev_farray(stim[i],ns[i],&mean,&sdev);*/
      /*binary_threshold_farray(stim[i],ns[i],mean);*/
      ;
    }else if (binflag == 3){
      /*adjust_ternary_stim(nd,j,awake_flag,stim[i],ns[i],pfile);*/
      ;
    }
  }
  *rstim = stim; *rns = ns; *rdt = dt; *rt0 = t0;
}
/**************************************-**************************************/
/*                                                                           */
/*                        GET_NDATA_EXPANDED_STIM_GROUP                      */
/*                                                                           */
/*  Called by 'STA_NDA' to get a stimulus stored in the ndata file.          */
/*                                                                           */
/*****************************************************************************/
void get_ndata_expanded_stim_group(nd,ndg,k,stimchan,resolution,fft_flag,
				   norm_flag,rstim,rns,rdt,rt0)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k;
     char stimchan[];
     float resolution;
     int norm_flag,fft_flag;
     float ***rstim;
     int **rns;
     float *rdt;
     int *rt0;
{
  int i;
  int n,tn,tt0,*ns,*tns;
  float res_in,tdt,*t,**stim,**tstim;

  get_ndata_stim_group(nd,ndg,k,stimchan,0,&tstim,&tns,&tdt,&tt0); //binflg=0

  n = ndg->cnt[k];
  stim = (float **)myalloc(n*sizeof(float *));
  ns = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++){
    res_in = 1000.0/tdt;
    t = tstim[i];
    tn = tns[i];
    if (fft_flag == 0)
      get_stim_at_resolution_simple(t,tn,res_in,resolution,&stim[i],&ns[i]);
    else if (fft_flag == 1)
      get_stim_at_resolution(t,tn,res_in,resolution,&stim[i],&ns[i]);
    else if (fft_flag == 2)
      get_stim_at_resolution_impulse(t,tn,res_in,resolution,&stim[i],&ns[i]);
    else if (fft_flag == 3)
      get_stim_at_resolution_piecewise(t,tn,res_in,resolution,&stim[i],&ns[i]);
    else
      exit_error("GET_NDATA_EXPANDED_STIM_GROUP","Unknown fft_flag value");

    if (norm_flag > 0){ /* Wyeth - changed to z-score to get rid of mean */
      z_score_farray(stim[i],ns[i]);
      /*norm_variance_farray(stim[i],ns[i]);*/
    }
  }
  free_2d_farray(tstim,n); myfree(tns);
  *rstim = stim; *rns = ns; *rdt = tdt; *rt0 = tt0;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_STD_STIM_GROUP                            */
/*                                                                           */
/*  Get a standard random stimulus sequence.                                 */
/*                                                                           */
/*****************************************************************************/
void get_std_stim_group(nd,ndg,k,stimstr,seedname,rstim,rns,rdt,rt0)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k;
     char *stimstr;
     char *seedname;   // e.g., 'seed', 'seed1', 'seed2' ...
     float ***rstim;   /* stim[i][ns[i]] stimulus values for ith trial */
     int **rns;        /* ns[i] number of stimulus values for ith trial */
     float *rdt;       /* time between stimulus values (msec) */
     int *rt0;         /* stimulus onset (msec) */
{
  int i,j,l;
  int *ns,n,t0,tn,nsl,nn,*t,seed,order,flag,mseq_ord_flag,mseq_mu_flag;
  float **stim,dt;
  char **slist,stype[SLEN],*tstr;

  n = ndg->cnt[k];
  stim = (float **)myalloc(n*sizeof(float *));
  ns = (int *)myalloc(n*sizeof(int));

  get_items_from_string(stimstr,&slist,&nsl);
  strcpy(stype,slist[0]);
  dt = atof(slist[1]);  /* In sampling units OR ms??? */
  t0 = atoi(slist[2]);  /* Start time in sampling units OR ms??? */
  tn = atoi(slist[3]);  /* Duration in sampling units OR ms??? */
  /*printf("TYPE = %s  DT= %f  t0= %d  tn= %d\n",stype,dt,t0,tn);*/
  free_2d_carray(slist,nsl);

  mseq_ord_flag = mseq_mu_flag = 0;
  nn = (float)tn/dt;
  for(i=0;i<n;i++){      /* For each trial in this group */
    j = ndg->tnum[k][i]; /*   get trial number */

    //flag = ndata_get_var_or_const_param(nd,j,"seed",&tstr);

    seed = ndata_get_vc_param_int_or_exit(nd,j,seedname,"GET_STD_STIM_GROUP");

    /*    
    flag = ndata_get_var_or_const_param(nd,j,seedname,&tstr);
    if (flag == 1){
      seed = ndata_get_vc_param_int_or_exit(nd,j,"seed","GET_STD_STIM_GROUP");
      //printf("trial %d  SEED = %d\n",j,seed);
    }else{
      // Added for 'dmmask'
      printf("  Reading seed as 'seed1'\n");
      seed = ndata_get_vc_param_int_or_exit(nd,j,"seed1","GET_STD_STIM_GROUP");
    }
    */

    ns[i] = nn;

    if ((strcmp(stype,"binary")==0) || (strcmp(stype,"mseq")==0)){
      if (strcmp(stype,"binary")==0)
	t = myrand_get_std_bin_seq(ns[i],seed);
      else if (strcmp(stype,"mseq")==0){

	flag = ndata_get_var_or_const_param(nd,j,"mseq_ord",&tstr);
	if (flag == 1){
	  mseq_ord_flag += 1;
	  //printf("  Reading mseq order as 'mseq_ord'\n");
	  order = ndata_get_vc_param_int_or_exit(nd,j,"mseq_ord",
						 "GET_STD_STIM_GROUP");
	}else{
	  mseq_mu_flag += 1;
	  //printf("  Reading mseq order as 'mu'\n");
	  order = ndata_get_vc_param_int_or_exit(nd,j,"mu",
						 "GET_STD_STIM_GROUP");
	}

	t = myrand_get_std_mseq(ns[i],seed,order);
      }
      stim[i] = (float *)myalloc(ns[i]*sizeof(float));
      for(l=0;l<ns[i];l++) // Convert from 0/1 to -1/1
	if (t[l] == 0)
	  stim[i][l] = -1;
	else
	  stim[i][l] = 1;
      myfree(t);
    }else if (strcmp(stype,"gaussian")==0){
      stim[i] = myrand_get_std_gauss_seq(ns[i],seed);
    }else
      exit_error("GET_STD_STIM_GROUP","Unknown standard stim type");
  }

  // If there were m-sequences, report how the order was set
  if ((mseq_ord_flag > 0) || (mseq_mu_flag > 0)){
    if (mseq_ord_flag == 0)
      printf("      Order of m-seq was determined from param 'mu'\n");
    else if (mseq_mu_flag == 0)
      printf("      Order of m-seq was determined from param 'mseq_order'\n");
    else
      exit_error("GET_STD_STIM_GROUP","Inconsistent m-seq order encoding.");
  }


  *rstim = stim; *rns = ns; *rdt = dt; *rt0 = t0;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_STD_STIM_GROUP_EXPANDED                      */
/*                                                                           */
/*****************************************************************************/
void get_std_stim_group_expanded(nd,ndg,k,stimstr,seedname,resolution,fft_flag,
				 norm_flag,rstim,rns,rdt,rt0)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k;
     char stimstr[];
     char *seedname;
     float resolution;
     int norm_flag,fft_flag;
     float ***rstim;
     int **rns;
     float *rdt;
     int *rt0;
{
  int i;
  int n,tn,tt0,*ns,*tns;
  float res_in,tdt,*t,**stim,**tstim;

  //printf("(noise_util)----------IN HERE seedname ==>%s<==\n",seedname);

  get_std_stim_group(nd,ndg,k,stimstr,seedname,&tstim,&tns,&tdt,&tt0);

  /*append_farray_plot("zz.stim","stim",tstim[0],tns[0],1);*/

  n = ndg->cnt[k];
  stim = (float **)myalloc(n*sizeof(float *));
  ns = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++){
    res_in = 1000.0/tdt;
    t = tstim[i];
    tn = tns[i];
    if (fft_flag == 0)
      get_stim_at_resolution_simple(t,tn,res_in,resolution,&stim[i],&ns[i]);
    else if (fft_flag == 1)
      get_stim_at_resolution(t,tn,res_in,resolution,&stim[i],&ns[i]);
    else if (fft_flag == 2)
      get_stim_at_resolution_impulse(t,tn,res_in,resolution,&stim[i],&ns[i]);
    else if (fft_flag == 3)
      get_stim_at_resolution_piecewise(t,tn,res_in,resolution,&stim[i],&ns[i]);
    else
      exit_error("GET_NDATA_EXPANDED_STIM_GROUP","Unknown fft_flag value");

    if (norm_flag > 0){
      z_score_farray(stim[i],ns[i]);
    }
  }
  free_2d_farray(tstim,n); myfree(tns);
  *rstim = stim; *rns = ns; *rdt = tdt; *rt0 = tt0;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_4STATE_STIM_TRIAL                          */
/*                                                                           */
/*  "csflag" - center-surround flag.                                         */
/*  0 - return full 4 states, indicating center and surround values.         */
/*  1 - return 2 states, indicating center value only.                       */
/*  2 - return 2 states, indicating surround value only.                     */
/*                                                                           */
/*  'dt' is estimated and returned, but not *used* here.                     */
/*                                                                           */
/*****************************************************************************/
void get_4state_stim_trial(nd,pfile,k,rstate,rdwell,rn,rdt)
     struct ndata_struct *nd;
     char pfile[];
     int k,**rstate,*rdwell,*rn;
     float *rdt;
{
  int i,j;
  int *s,ns,n,framerate,seed,dwell,period,count1,flag,pat1amp,pat2amp,csflag;
  int rpt2,patflag,ternflag,tf;
  float frv,dt,tf1,stn,tscale,con1,con2;

  ternflag = 0;

  flag = ndata_get_vc_param_int(nd,k,"seed",&seed);
  if (!flag) exit_error("GET_4STATE_STIM_TRIAL","Param seed not found");

  if (strcmp(nd->class,"MODEL")==0){
    framerate = 1000; /* Make up some value */
    dt = 1.0;       /* Actual dt in msec */
    tscale = 0.001;  /* Make up some value */
    flag = ndata_get_vc_param_float(nd,k,"tf1",&tf1);
    if (!flag) exit_error("GET_4STATE_STIM_TRIAL","tf1 not found");
    dwell = my_rint(1.0/(tf1*tscale));
    printf("    dwell = %d\n",dwell);

    flag = ndata_get_vc_param_float(nd,k,"stn",&stn);
    if (!flag) exit_error("GET_4STATE_STIM_TRIAL","Prm stn not found");
    period = (int)stn;  /* May truncate */

    flag = ndata_get_vc_param_float(nd,k,"contrast1",&con1);
    if (!flag) exit_error("GET_4STATE_STIM_TRIAL","Prm contrast1 not found");
    flag = ndata_get_vc_param_float(nd,k,"contrast2",&con2);
    if (!flag) exit_error("GET_4STATE_STIM_TRIAL","Prm contrast2 not found");
    pat1amp = (int)(con1*100.0);
    pat2amp = (int)(con2*100.0);
    if ((pat1amp > 0)&&(pat2amp > 0))
      csflag = 0; /* Center and surround */
    else if (pat2amp == 0)
      csflag = 1; /* Center only */
    else
      csflag = 2;  /* Surround only */

  }else{

    flag = ndata_get_vc_param_int(nd,k,"framrate",&framerate);
    if (!flag) exit_error("GET_4STATE_STIM_TRIAL","framerate not found");

    flag = ndata_get_vc_param_int(nd,k,"dwell",&dwell);
    if (!flag){
      flag = ndata_get_vc_param_int(nd,k,"rpt1",&dwell);
      if (!flag)
	exit_error("GET_4STATE_STIM_TRIAL","Param dwell (rpt1) not found");
      flag = ndata_get_vc_param_int(nd,k,"rpt2",&rpt2);
      if (flag)
	if (rpt2 != dwell)
	  exit_error("GET_4STATE_STIM_TRIAL","rpt1 != rpt2");
    }
    if (dwell == -1){
      flag = ndata_get_vc_param_int(nd,k,"maxj1",&tf);
      if (!flag) exit_error("GET_4STATE_STIM_TRIAL","maxj1 not found");
      dwell = 1024/tf;
    }
    
    flag = ndata_get_vc_param_int(nd,k,"duration",&period);
    if (!flag) exit_error("GET_4STATE_STIM_TRIAL","Prm duration not found");
    
    if (strcmp(pfile,"wy9st.p")==0)
      ternflag = 1;
    
    flag = ndata_get_vc_param_float(nd,k,"actual_dt",&dt);
    if (!flag){
      if ((strcmp(pfile,"x4st.p")==0)||(strcmp(pfile,"x4stp.p")==0)||
	  (strcmp(pfile,"x4stp2.p")==0)||(strcmp(pfile,"x4std.p")==0)){
	printf("  *** WARNING:  Estimating SGI dt old way, PROBABLY wrong.\n");
	dt = get_sgi_actual_dt(framerate,0); /* 0 means OCTANE (not O2) */
      }else{
	dt = get_wn_actual_dt(framerate,0);
      }
    }
    
    if (strcmp(pfile,"wy4sum.p")==0){
      flag = ndata_get_var_value_trial_int(nd,k,"patflag",&patflag);
      if (!flag) exit_error("GET_4STATE_STIM_TRIAL","Prm patflag not found");
      csflag = patflag;
    }else if ((strcmp(pfile,"x4st.p")==0)||(strcmp(pfile,"x4stp.p")==0)||
	      (strcmp(pfile,"x4stp2.p")==0)||
	      (strcmp(pfile,"x4std.p")==0)){
      flag = ndata_get_vc_param_int(nd,k,"pat1amp",&pat1amp);
      if (!flag) exit_error("GET_4STATE_STIM_TRIAL","Prm pat1amp not found");
      flag = ndata_get_vc_param_int(nd,k,"pat2amp",&pat2amp);
      if (!flag) exit_error("GET_4STATE_STIM_TRIAL","Prm pat2amp not found");
      if ((pat1amp > 0)&&(pat2amp > 0))
	csflag = 0; /* Center and surround */
      else if (pat2amp == 0)
	csflag = 1; /* Center only */
      else
	csflag = 2;  /* Surround only */
    }else{ /* For wy4st, etc */
      /*flag = ndata_get_var_value_trial_int(nd,k,"pat1amp",&pat1amp);*/
      flag = ndata_get_vc_param_int(nd,k,"pat1amp",&pat1amp);
      if (!flag) exit_error("GET_4STATE_STIM_TRIAL","Prm pat1amp not found");
      /*flag = ndata_get_var_value_trial_int(nd,k,"pat2amp",&pat2amp);*/
      flag = ndata_get_vc_param_int(nd,k,"pat2amp",&pat2amp);
      if (!flag) exit_error("GET_4STATE_STIM_TRIAL","Prm pat2amp not found");
      if ((pat1amp > 0)&&(pat2amp > 0))
	csflag = 0; /* Center and surround */
      else if (pat2amp == 0)
	csflag = 1; /* Center only */
      else
	csflag = 2;  /* Surround only */
    }
  }

  n = framerate * period; /* Number of stimulus frames. */
  if (seed > 0) /*** Initialize random number generator. ***/
    seed = -seed;
  noise_util_ran2(&seed);

  ns = n/dwell + 1;
  s = (int *)myalloc(ns*sizeof(int));
  count1 = dwell;
  j = 0;
  for(i=0;i<n;i++){
    if (count1 == dwell){ /*** Choose new state. ***/
      frv = noise_util_ran2(&seed);
      if (ternflag == 0){
	if (frv < 0.25)
	  s[j] = 3;               /* state 3 = PP (center,surround) */
	else if (frv < 0.50)
	  s[j] = 2;               /* state 2 = PA */
	else if (frv < 0.75)
	  s[j] = 1;               /* state 1 = AP */
	else
	  s[j] = 0;               /* state 0 = AA */
	j += 1;
	count1 = 0;
      }else{
	/*gflag = 0;*/
	if (frv < 0.33333){
	  frv = noise_util_ran2(&seed);
	  if (frv < 0.33333){
	    s[j] = 3; /* vdpk = 0;*/ /* PP */
	  }else if (frv < 0.66666){
	    s[j] = 2; /* vdpk = 1;*/ /* PA */
	  }else{
	    s[j] = 6; /* vdpk = 1;*/ /* P0 */
	    /*gflag = 2;*/
	  }
	}else if (frv < 0.66666){
	  frv = noise_util_ran2(&seed);
	  if (frv < 0.33333){
	    s[j] = 1; /* vdpk = 2;*/ /* AP */
	  }else if (frv < 0.66666){
	    s[j] = 0; /* vdpk = 3;*/ /* AA */
	  }else{
	    s[j] = 7; /* vdpk = 3;*/ /* A0 */
	    /*gflag = 1;*/
	  }
	}else{
	  frv = noise_util_ran2(&seed);
	  if (frv < 0.33333){
	    s[j] = 5; /* vdpk = 2;*/ /* 0P */
	    /*gflag = 2;*/
	  }else if (frv < 0.66666){
	    s[j] = 8; /* vdpk = 3;*/ /* 0A */
	    /*gflag = 3;*/
	  }else{
	    s[j] = 4; /* vdpk = 3;*/ /* 00 */
	    /*gflag = 2;*/
	  }
	}

	j += 1;
	count1 = 0;
      }
    }
    count1 += 1;
  }
  if (j > ns)
    exit_error("GET_4STATE_STIM_TRIAL","j too large");

  if (ternflag == 0){
    if (csflag == 1) /*** Report center state as 0,1 ***/
      for(i=0;i<j;i++)
	if (s[i] >= 2)
	  s[i] = 1;
	else
	  s[i] = 0;
    else if (csflag == 2)
      for(i=0;i<j;i++){
	if (s[i] == 3)
	  s[i] = 1;
	if (s[i] == 2)
	  s[i] = 0;
      }
  }else{
    if (csflag == 1) /*** Report center state as 0,1,2 ***/
      for(i=0;i<j;i++)
	if (s[i]==0 || s[i]==1 || s[i]==7)
	  s[i] = 0;
	else if (s[i]==2 || s[i]==3 || s[i]==6)
	  s[i] = 1;
	else
	  s[i] = 2;
    else if (csflag == 2) /* Surr */
      for(i=0;i<j;i++){
	if (s[i]==0 || s[i]==2 || s[i]==8)
	  s[i] = 0;
	else if (s[i]==1 || s[i]==3 || s[i]==5)
	  s[i] = 1;
	else
	  s[i] = 2;
      }
  }
  
  /***
printf("seed = %d\n",kseed);
for(i=0;i<j;i++)
  printf("%d ",s[i]);
printf("\n");***/

  printf("    n=%d  dwell = %d  j=%d  dt= %f\n",n,dwell,j,dt);
  
  *rstate = s; *rdwell = dwell; *rn = j; *rdt = dt;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_4STATE_STIM_LIM_TRIAL                         */
/*                                                                           */
/*  "csflag" - center-surround flag.                                         */
/*  0 - return full for states, indicating center and surround values.       */
/*  1 - return 2 states, indicating center value only.                       */
/*  2 - return 2 states, indicating surround value only.                     */
/*                                                                           */
/*****************************************************************************/
void get_4state_stim_lim_trial(nd,k,rstate,rdwell,rn,rdt)
     struct ndata_struct *nd;
     int k,**rstate,*rdwell,*rn;
     float *rdt;
{
  int i,j;
  int *s,ns,n,framerate,dwell,period,count1,flag,pat1amp,pat2amp,csflag;
  int tf1,tf2,ncycles,nfmin,step;
  float dt;

  flag = ndata_get_vc_param_int(nd,k,"framrate",&framerate);
  if (!flag) exit_error("GET_4STATE_STIM_LIM_TRIAL","framerate not found");
  flag = ndata_get_vc_param_int(nd,k,"tf1",&tf1);
  if (!flag) exit_error("GET_4STATE_STIM_LIM_TRIAL","tf1 not found");
  flag = ndata_get_vc_param_int(nd,k,"tf2",&tf2);
  if (!flag) exit_error("GET_4STATE_STIM_LIM_TRIAL","tf2 not found");
  if (tf1 != tf2)
    exit_error("GET_4STATE_STIM_LIM_TRIAL","tf1 != tf2");
  flag = ndata_get_vc_param_int(nd,k,"ncycles",&ncycles);
  if (!flag) exit_error("GET_4STATE_STIM_LIM_TRIAL","ncycles not found");
  
  flag = ndata_get_vc_param_int(nd,k,"duration",&period);
  if (!flag) exit_error("GET_4STATE_STIM_LIM_TRIAL","duration not found");
  dt = get_wn_actual_dt(framerate,0);

  flag = ndata_get_vc_param_int(nd,k,"pat1amp",&pat1amp);
  if (!flag) exit_error("GET_4STATE_STIM_LIM_TRIAL","Param pat1amp not found");
  flag = ndata_get_vc_param_int(nd,k,"pat2amp",&pat2amp);
  if (!flag) exit_error("GET_4STATE_STIM_LIM_TRIAL","Param pat2amp not found");

  if ((pat1amp > 0)&&(pat2amp > 0))
    csflag = 0; /* Center and surround */
  else if (pat2amp == 0)
    csflag = 1; /* Center only */
  else
    csflag = 2;  /* Surround only */

  dwell = 1024/tf1 * ncycles; /* Number of frames per step. */
  nfmin = period * framerate; /* Minimum duration requested. */
  n = dwell*4; /* Nframes to go once through all 4 steps. */
  while (n < nfmin)
    n += dwell*4; /* Nframes to go once through all 4 steps. */

  ns = n/dwell;
  s = (int *)myalloc(ns*sizeof(int));
  count1 = dwell;
  step = -1;

  j = 0;
  for(i=0;i<n;i++){
    if (count1 == dwell){ /*** Choose new state. ***/
      count1 = 0;
      step += 1; /* Update step */
      if (step >=4)
	step = 0;

      if (step == 0)
	s[j] = 0; /* NN null center null surround */
      else if (step == 1)
	s[j] = 2; /* PN */
      else if (step == 2)
	s[j] = 3; /* PP */
      else if (step == 3)
	s[j] = 2; /* PN */
      j += 1;

    }
    count1 += 1;
  }
  if (j != ns)
    exit_error("GET_4STATE_STIM_LIM_TRIAL","j too large");

  if (csflag == 1) /*** Report center state as 0,1 ***/
    for(i=0;i<j;i++)
      if (s[i] >= 2)
	s[i] = 1;
      else
	s[i] = 0;
  else if (csflag == 2)
    for(i=0;i<j;i++){
      if (s[i] == 3)
	s[i] = 1;
      if (s[i] == 2)
	s[i] = 0;
    }
  
  *rstate = s; *rdwell = dwell; *rn = ns; *rdt = dt;
}
/**************************************-**************************************/
/*                                                                           */
/*                       GET_4STATE_STIM_WY2PAT_TRIAL                        */
/*                                                                           */
/*  "csflag" - center-surround flag.                                         */
/*  0 - return full for states, indicating center and surround values.       */
/*  1 - return 2 states, indicating center value only.                       */
/*  2 - return 2 states, indicating surround value only.                     */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - "dwell" is set to "rpt1".                                              */
/*                                                                           */
/*****************************************************************************/
void get_4state_stim_wy2pat_trial(nd,k,rstate,rdwell,rn,rfdt)
     struct ndata_struct *nd;
     int k,**rstate,*rdwell,*rn;
     float *rfdt;
{
  int i,j;
  int *s,ns,n,flag,csflag;
  int dwell,rpt2,awake_flag,maxj1,maxj2,pat1amp,pat2amp;
  float *xdata,*y1data,*y2data;
  char *pfile;

  /*printf("  GET_4STATE_STIM_WY2PAT_TRIAL\n");*/

  flag = ndata_get_vc_param_int(nd,k,"rpt1",&dwell);
  if (!flag) exit_error("GET_4STATE_STIM_WY2PAT_TRIAL","Param rpt1 not found");
  flag = ndata_get_vc_param_int(nd,k,"rpt2",&rpt2);
  if (!flag) exit_error("GET_4STATE_STIM_WY2PAT_TRIAL","Param rpt2 not found");
  if (rpt2 != dwell)
    exit_error("GET_4STATE_STIM_WY2PAT_TRIAL","rpt1 != rpt2");

  /*** Var params. ***/
  /*flag = ndata_get_var_value_trial_int(nd,k,"pat1amp",&pat1amp);*/
  flag = ndata_get_vc_param_int(nd,k,"pat1amp",&pat1amp);
  if (!flag) exit_error("GET_4STATE_STIM_WY2PAT_TRIAL","pat1amp not found");
  /*flag = ndata_get_var_value_trial_int(nd,k,"pat2amp",&pat2amp);*/
  flag = ndata_get_vc_param_int(nd,k,"pat2amp",&pat2amp);
  if (!flag) exit_error("GET_4STATE_STIM_WY2PAT_TRIAL","pat2amp not found");

  /*** Const. (or var?) params. ***/
  flag = ndata_get_vc_param_int(nd,k,"maxj1",&maxj1);
  if (!flag) exit_error("GET_4STATE_STIM_WY2PAT_TRIAL","maxj1 not found");
  flag = ndata_get_vc_param_int(nd,k,"maxj2",&maxj2);
  if (!flag) exit_error("GET_4STATE_STIM_WY2PAT_TRIAL","maxj2 not found");

  if ((pat1amp > 0)&&(pat2amp > 0))
    csflag = 0; /* Center and surround */
  else if (pat2amp == 0)
    csflag = 1; /* Center only */
  else
    csflag = 2;  /* Surround only */

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0)
    exit_error("GET_NEW_NOISE_SIMPLE_STIM_GROUP","Param pfile not found");

  awake_flag = 0;
  get_wy2pat_stim_for_trial_02(nd,k,awake_flag,pfile,&xdata,&y1data,&y2data,
			       &n,rfdt);

  /*** Now deduce sequence of 4 states based on numerical jump values. ***/
  ns = n/dwell; /* A partial last state is ignored. */

  s = (int *)myalloc(ns*sizeof(int));
  for(i=0;i<ns;i++){
    j = i*dwell;
    if (y1data[j] == (float)maxj1){
      if (y2data[j] == (float)maxj2)
	s[i] = 3;               /* state 3 = pref/pref (center/surround) */
      else
	s[i] = 2;               /* state 2 = pref/null */
    }else{
      if (y2data[j] == (float)maxj2)
	s[i] = 1;               /* state 1 = null/pref */
      else
	s[i] = 0;               /* state 0 = null/null */
    }
  }

  if (csflag == 1) /*** Report center state as 0,1 ***/
    for(i=0;i<ns;i++)
      if (s[i] >= 2)
	s[i] = 1;
      else
	s[i] = 0;
  else if (csflag == 2)
    for(i=0;i<ns;i++){
      if (s[i] == 3)
	s[i] = 1;
      if (s[i] == 2)
	s[i] = 0;
    }

  myfree(xdata); myfree(y1data); myfree(y2data); myfree(pfile);
  
  *rstate = s; *rdwell = dwell; *rn = ns;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_FOUR_STATE_STIM_GROUP                         */
/*                                                                           */
/*  This handles the pfile "wy4st.p".                                        */
/*                                                                           */
/*  The value of 'dt' is determined and returned, but not used for           */
/*  calculations.                                                            */
/*                                                                           */
/*****************************************************************************/
void get_four_state_stim_group(nd,ndg,k,rstim,rns,rdwell,rdt)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k,***rstim,**rns,*rdwell;
     float *rdt;
{
  int i;
  int n,flag,**st,*ns,dwell;
  float dt;
  char *pfile,*filename;

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0)
    exit_error("GET_FOUR_STATE_STIM_GROUP","Const param pfile not found");
  if ((strcmp(pfile,"wy4st.p")==0) ||
      (strcmp(pfile,"wy4stg.p")==0) ||
      (strcmp(pfile,"wy4stio.p")==0) ||
      (strcmp(pfile,"wy9st.p")==0) ||
      (strcmp(pfile,"wy4sum.p")==0) ||
      (strcmp(pfile,"x4st.p")==0) ||
      (strcmp(pfile,"x4stp.p")==0) ||
      (strcmp(pfile,"x4stp2.p")==0) ||
      (strcmp(pfile,"x4std.p")==0) ||
      (strcmp(pfile,"wy2flat.p")==0) ||
      (strcmp(pfile,"wy2cspha.p")==0) ||
      (strcmp(pfile,"wy4stlim.p")==0) ||
      (strcmp(pfile,"wy4stfar.p")==0) ||
      (strcmp(pfile,"wy4stdir.p")==0)){
    flag = ndata_get_const_param(nd,"filename",&filename);
    if (flag > 0){
      if ((strcmp(filename,"452l004.p17")==0)||
          (strcmp(filename,"452l013.p07")==0)||
	  (strcmp(filename,"452l013.p08")==0)||
	  (strcmp(filename,"452l019.p10")==0)||
	  (strcmp(filename,"452l021.p08")==0)||
	  (strcmp(filename,"452l024.p08")==0)||
	  (strcmp(filename,"452l026.p09")==0)||
	  (strcmp(filename,"452l025.p09")==0)||
	  (strcmp(filename,"452l034.p10")==0)){
	get_special_wy4st_group(filename,nd,ndg,k,rstim,rns,rdwell,rdt);
	myfree(filename);
	return;
      }
    }

    n = ndg->cnt[k];
    st = (int **)myalloc(n*sizeof(int *));
    ns = (int *)myalloc(n*sizeof(int));
    for(i=0;i<n;i++){
      if ((strcmp(pfile,"wy4st.p")==0)||(strcmp(pfile,"wy2flat.p")==0)||
	  (strcmp(pfile,"wy4stg.p")==0)||(strcmp(pfile,"wy4stio.p")==0)||
	  (strcmp(pfile,"x4stp.p")==0)||(strcmp(pfile,"x4std.p")==0)||
	  (strcmp(pfile,"x4stp2.p")==0)||
	  (strcmp(pfile,"x4st.p")==0)||(strcmp(pfile,"wy4sum.p")==0)||
	  (strcmp(pfile,"wy9st.p")==0)){
	get_4state_stim_trial(nd,pfile,ndg->tnum[k][i],&st[i],&dwell,
				  &ns[i],&dt);
      }else if ((strcmp(pfile,"wy4stlim.p")==0)||
	       (strcmp(pfile,"wy4stfar.p")==0))
	get_4state_stim_lim_trial(nd,ndg->tnum[k][i],&st[i],&dwell,
				  &ns[i],&dt);
      else /*if (strcmp(pfile,"wy4stdir.p")==0)*/
	get_4state_stim_wy2pat_trial(nd,ndg->tnum[k][i],&st[i],&dwell,
				     &ns[i],&dt);
      /*** Next line is for printing out stimulus. ***/
      /*append_iarray_plot("zzz.stim",ndg->name[k],st[i],ns[i],1);*/
    }
    *rstim = st; *rns = ns; *rdwell = dwell; *rdt = dt;
    myfree(pfile);
  }else
    exit_error("GET_FOUR_STATE_STIM_GROUP","Inappropriate pfile");
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_TRIG_FOR_4STATE_STIM                         */
/*                                                                           */
/*  Compute triggers for the specified pattern.  "dt" is the time between    */
/*  states in "st".  The trigger time is given relative to the last state    */
/*  in the pattern "ipat".                                                   */
/*                                                                           */
/*****************************************************************************/
void get_trig_for_4state_stim(st,ns,n,dt,ipat,npat,rs,rcnt)
     int **st,*ns,n;
     float dt;
     int *ipat,npat,***rs,**rcnt;
{
  int i,j,k;
  int **s,*cnt,*t,nt,flag;

  s = (int **)myalloc(n*sizeof(int *));
  cnt = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++){ /* For each stimulus. */
    t = (int *)myalloc(ns[i]*sizeof(int));
    nt = 0;
    for(j=0;j<(ns[i]-npat);j++){ /* For each starting state. */
      flag = 1; /*** Assume a match ***/
      for(k=0;k<npat;k++) /* For each pattern state. */
	if (st[i][j+k] != ipat[k])
	  flag = 0;
      if (flag == 1){
	t[nt] = my_rint(dt*(float)(j+npat-1));
	nt += 1;
      }
    }
    s[i] = copy_iarray(t,nt);
    cnt[i] = nt;
    myfree(t);
  }
  *rs = s; *rcnt = cnt;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_4STATE_REL_TRIG_TRIAL                        */
/*                                                                           */
/*  Fill in the "trig" and "ntrig" arrays with the trigger times for the     */
/*  variously timed transitions.                                             */
/*                                                                           */
/*****************************************************************************/
void get_4state_rel_trig_trial(pfile,nd,k,trig,ntrig,ntt)
     char pfile[];
     struct ndata_struct *nd;
     int k,***trig,**ntrig,ntt;
{
  int i,j;
  int framerate,seed,period,tf1,tf2,deltamin,deltamax,ncycles;
  int *tdata,tn,tcount,n1,n2,*t1,*t2,dwell;
  float sdt;


  if ((strcmp(pfile,"wy4strel.p")==0)||(strcmp(pfile,"wy4strev.p")==0)||
      (strcmp(pfile,"wy4stref.p")==0)){
    get_wy2pat_rel_stim_params(nd,k,&framerate,&seed,&period,
			       &tf1,&tf2,&deltamin,&deltamax,&ncycles);
  }else if (strcmp(pfile,"wy2flatr.p")==0){
    get_wy2flat_rel_stim_params(nd,k,&framerate,&seed,&period,
				&deltamin,&deltamax,&dwell);
    tf1 = 0;
    tf2 = dwell;
    ncycles = 0;
  }

  get_4state_rel_stim_for_params(framerate,seed,period,tf1,tf2,deltamin,
				 deltamax,ncycles,&sdt,&tdata,&tn);

  if (strcmp(pfile,"wy4stref.p")==0){
    tcount = deltamax - deltamin + 1;
    if (tcount != ntt)
      exit_error("GET_4STATE_REL_TRIG_TRIAL","Bad ntt value");
    
    t1 = (int *)myalloc(tn*sizeof(int));
    for(i=0;i<tcount;i++){
      n1 = 0;
      for(j=0;j<tn;j++)
	if (tdata[j] == (i+deltamin)){
	  t1[n1] = (int)(((float)j+0.5)*sdt);
	  n1 += 1;
	}
      trig[i][k] = copy_iarray(t1,n1);
      ntrig[i][k] = n1;
    }
    myfree(t1);
  }else{ /*** For wy4strel, etc, where there are 0-3 and 1-2 transitions ***/
    tcount = deltamax - deltamin + 1;
    if (tcount != (ntt/2))
      exit_error("GET_4STATE_REL_TRIG_TRIAL","Bad ntt value");
    
    t1 = (int *)myalloc(tn*sizeof(int));
    t2 = (int *)myalloc(tn*sizeof(int));
    for(i=0;i<tcount;i++){
      n1 = n2 = 0;
      for(j=0;j<tn;j++)
	if (tdata[j] == (i+deltamin)){
	  if ((j%2) == 0){
	    t1[n1] = (int)(((float)j+0.5)*sdt);
	    n1 += 1;
	  }else{
	    t2[n2] = (int)(((float)j+0.5)*sdt);
	    n2 += 1;
	  }
	}
      trig[i][k] = copy_iarray(t1,n1);
      ntrig[i][k] = n1;
      trig[i+tcount][k] = copy_iarray(t2,n2);
      ntrig[i+tcount][k] = n2;
    }
    myfree(t1); myfree(t2);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_TRIG_4STATE_REL_GROUP                       */
/*                                                                           */
/*  Compute triggers for each of the different timing values for both of     */
/*  the transitions for stimuli specified by "wy4strel.p".                   */
/*                                                                           */
/*****************************************************************************/
void get_trig_4state_rel_group(nd,ndg,k,rtrig,rntrig,rntt,rttname)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k,****rtrig,***rntrig,*rntt;
     char ***rttname;
{
  int i;
  int ***trig,**ntrig,ntt,n,stimtype;
  int flag,framerate,seed,period,ncycles,deltamin,deltamax,tf1,tf2,dwell;
  char **ttname,*pfile,tname[SLEN];

  stimtype = ndata_get_const_param_int(nd,"stimtype");
  printf("    stimtype %d\n",stimtype); /* 8 is wy4strel, 15 is wy4stref */

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0) exit_error("GET_TRIG_4STATE_REL_GROUP","'pfile' not found");
  if ((strcmp(pfile,"wy4strel.p")==0)||(strcmp(pfile,"wy4strev.p")==0)||
      (strcmp(pfile,"wy4stref.p")==0)){
    /*** Look at first trial in group to determine number of timing types. ***/
    get_wy2pat_rel_stim_params(nd,ndg->tnum[k][0],&framerate,&seed,&period,
			       &tf1,&tf2,&deltamin,&deltamax,&ncycles);
  }else if (strcmp(pfile,"wy2flatr.p")==0){
    get_wy2flat_rel_stim_params(nd,ndg->tnum[k][0],&framerate,&seed,&period,
				&deltamin,&deltamax,&dwell);
  }else{
    printf("pfile ==>%s<==\n",pfile);
    exit_error("GET_TRIG_4STATE_REL_GROUP","Wrong pfile");
  }

  n = ndg->cnt[k];
  if (stimtype == 8)
    ntt = 2 * (deltamax - deltamin + 1);
  else
    ntt = deltamax - deltamin + 1;

  trig = get_2d_pointer_iarray(ntt,n);
  ntrig = get_2d_iarray(ntt,n);
  ttname = (char **)myalloc(ntt*sizeof(char *));

  if (stimtype == 8){
    for(i=0;i<ntt/2;i++){
      sprintf(tname,"0-3_t=%d",deltamin+i);
      ttname[i] = strdup(tname);
      sprintf(tname,"1-2_t=%d",deltamin+i);
      ttname[i+ntt/2] = strdup(tname);
    }
  }else{ /* stimtype = 15, wy4stref */
    for(i=0;i<ntt;i++){
      sprintf(tname,"0-3_t=%d",deltamin+i);
      ttname[i] = strdup(tname);
    }
  }

  for(i=0;i<n;i++)
    get_4state_rel_trig_trial(pfile,nd,ndg->tnum[k][i],trig,ntrig,ntt);

  *rtrig = trig; *rntrig = ntrig; *rntt = ntt; *rttname = ttname;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_4STATE_SIMPLE_TRIG_TRIAL                      */
/*                                                                           */
/*  Fill in the "trig" and "ntrig" arrays with the trigger times for the     */
/*  variously timed transitions.                                             */
/*                                                                           */
/*****************************************************************************/
void get_4state_simple_trig_trial(pfile,nd,k,tindex,trig,ntrig,nph)
     char pfile[];
     struct ndata_struct *nd;
     int k,tindex,***trig,**ntrig,nph;
{
  int i,j;
  int framerate,seed,period,tf1,tf2,ncycles,nphase,pat1amp,pat2amp;
  int *tdata,tn,n1,*t1;
  float sdt,st0;

  get_wy2pat_simple_stim_params(nd,k,&framerate,&seed,&period,
				&tf1,&tf2,&nphase,&ncycles,&pat1amp,&pat2amp);
  get_4state_simple_stim_for_params(pfile,framerate,seed,period,tf1,tf2,nphase,
				    ncycles,&st0,&sdt,&tdata,&tn);

  t1 = (int *)myalloc(tn*sizeof(int));
  for(i=0;i<nphase;i++){ /*** WYETH - changed tcount to nphase?????????????*/
    n1 = 0;
    for(j=0;j<tn;j++) /* For each trigger */
      if (tdata[j] == i){
	/*t1[n1] = (int)(((float)j+0.5)*sdt);*/
	t1[n1] = (int)(st0 + (float)j*sdt);
	n1 += 1;
      }
    trig[i][tindex] = copy_iarray(t1,n1);
    ntrig[i][tindex] = n1;
  }
  myfree(t1);
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_TRIG_4STATE_SIMPLE_GROUP                      */
/*                                                                           */
/*  Compute triggers for each of the different phase types for 'wy4simp.p'.  */
/*                                                                           */
/*****************************************************************************/
void get_trig_4state_simple_group(nd,ndg,k,rtrig,rntrig,rnph,rphname)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k,****rtrig,***rntrig,*rnph;
     char ***rphname;
{
  int i,j,l;
  int ***trig,**ntrig,n;
  int flag,framerate,seed,period,ncycles,nphase,tf1,tf2,pat1amp,pat2amp;
  char **phname,*pfile,tname[SLEN3],tstr[SLEN];

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0) exit_error("GET_TRIG_4STATE_SIMPLE_GROUP","pfile not found");
  if ((strcmp(pfile,"wy4simp.p")!=0)&&(strcmp(pfile,"x4simp.p")!=0)){
    printf("*** pfile = %s\n",pfile);
    exit_error("GET_TRIG_4STATE_SIMPLE_GROUP","inappropriate pfile");
  }

  /*** Look at first trial in group to determine number of timing types. ***/
  get_wy2pat_simple_stim_params(nd,ndg->tnum[k][0],&framerate,&seed,&period,
				&tf1,&tf2,&nphase,&ncycles,&pat1amp,&pat2amp);
  n = ndg->cnt[k];

  trig = get_2d_pointer_iarray(nphase,n);
  ntrig = get_2d_iarray(nphase,n);
  phname = (char **)myalloc(nphase*sizeof(char *));

  printf("AMP = %d %d\n",pat1amp,pat2amp);
  if ((pat1amp > 0)&&(pat2amp > 0))
    strcpy(tstr,"CS");
  else if (pat1amp == 0)
    strcpy(tstr,"S");
  else
    strcpy(tstr,"C");

  for(i=0;i<nphase;i++){
    sprintf(tname,"%d/%d_%s",i,nphase,tstr);
    phname[i] = strdup(tname);
  }
  for(i=0;i<n;i++)
    get_4state_simple_trig_trial(pfile,nd,ndg->tnum[k][i],i,trig,ntrig,nphase);

  for(i=0;i<nphase;i++){ /* phase */
    for(j=0;j<n;j++){ /* trial */
      printf("%4d %4d  %4d  ",i,j,ntrig[i][j]);
      for(l=0;l<ntrig[i][j];l++)
	printf("%d ",trig[i][j][l]);
      printf("\n");
    }
  }

  *rtrig = trig; *rntrig = ntrig; *rnph = nphase; *rphname = phname;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_WY2PA_STIM_GROUP                           */
/*                                                                           */
/*  This handles the pfile "wy2pa.p".                                        */
/*                                                                           */
/*****************************************************************************/
void get_wy2pa_stim_group(nd,ndg,k,rpstate,rtanti,rns,radur,rpdur,rdt)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k,***rpstate,***rtanti,**rns,*radur,*rpdur;
     float *rdt;
{
  int i;
  int n,flag,**pstate,**tanti,*ns,adur,pdur;
  float dt;
  char *pfile;

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0)
    exit_error("GET_WY2PA_STIM_GROUP","Const param pfile not found");
  if ((strcmp(pfile,"wy2pa.p")==0)||(strcmp(pfile,"wy2pp.p")==0)||
      (strcmp(pfile,"wy2ppc3.p")==0)||(strcmp(pfile,"wy2ppsep.p")==0)){
    n = ndg->cnt[k];
    pstate = (int **)myalloc(n*sizeof(int *));
    tanti = (int **)myalloc(n*sizeof(int *));
    ns = (int *)myalloc(n*sizeof(int));
    for(i=0;i<n;i++)
      get_wy2pa_stim_for_trial(nd,ndg->tnum[k][i],pfile,&pstate[i],&tanti[i],
			       &ns[i],&adur,&pdur,&dt);
    *rpstate = pstate; *rtanti = tanti; *rns = ns; *rdt = dt; *radur = adur;
    *rpdur = pdur;
  }else
    exit_error("GET_WY2PA_STIM_GROUP","Inappropriate pfile");

  myfree(pfile);
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_WYTFV_STIM_GROUP                           */
/*                                                                           */
/*  This handles the pfile "wytfv.p".                                        */
/*                                                                           */
/*****************************************************************************/
void get_wytfv_stim_group(nd,ndg,k,rstima,rstimp,rns,rdt,rtflist,rntf)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k,***rstima,***rstimp,**rns;
     float *rdt;
     int **rtflist,*rntf;
{
  int i,j;
  int n,flag,**stima,**stimp,*ns,ntf,zflag,maxtf,*tflist;
  float dt;
  char *pfile;

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0)
    exit_error("GET_WYTFV_STIM_GROUP","Const param pfile not found");
  if (strcmp(pfile,"wytfv.p")==0){
    n = ndg->cnt[k];

    j = ndg->tnum[k][0];
    flag = ndata_get_vc_param_int(nd,j,"zflag",&zflag);
    if (!flag) exit_error("GET_WYTFV_STIM_GROUP","zflag not found");
    flag = ndata_get_vc_param_int(nd,j,"ntf",&ntf);
    if (!flag) exit_error("GET_WYTFV_STIM_GROUP","ntf not found");
    flag = ndata_get_vc_param_int(nd,j,"maxtf",&maxtf);
    if (!flag) exit_error("GET_WYTFV_STIM_GROUP","maxtf not found");

    if (zflag == 1)
      ntf += 1;
    tflist = (int *)myalloc(ntf*sizeof(int));
    tflist[0] = maxtf;
    for(i=1;i<ntf;i++)
      tflist[i] = tflist[i-1]/2;
    if (zflag == 1)
      tflist[ntf-1] = 0;

    stima = (int **)myalloc(n*sizeof(int *));
    stimp = (int **)myalloc(n*sizeof(int *));
    ns = (int *)myalloc(n*sizeof(int));
    for(i=0;i<n;i++)
      get_wytfv_stim_for_trial(nd,ndg->tnum[k][i],pfile,&stima[i],&stimp[i],
			       &ns[i],&dt);
    *rstima = stima; *rstimp = stimp; *rns = ns; *rdt = dt; *rtflist = tflist;
    *rntf = ntf;
  }else
    exit_error("GET_WYTFV_STIM_GROUP","Inappropriate pfile");

  myfree(pfile);
}
/**************************************-**************************************/
/*                                                                           */
/*                              WRITE_WYTFV_STIM                             */
/*                                                                           */
/*  For dumping out stimulus sequence.  E.G. for matlab analysis.            */
/*                                                                           */
/*****************************************************************************/
void write_wytfv_stim(outfile,stima,stimp,ns,n,dt,tflist,ntf)
     char outfile[];
     int **stima,**stimp,*ns,n;
     float dt;
     int *tflist,ntf;
{
  FILE *fopen(),*fout;
  int i,j;

  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** File %s\n",outfile);
    exit_error("WRITE_WYTFV_STIM","Cannot open file");
  }

  if (!is_const_iarray(ns,n))
    exit_error("WRITE_WYTFV_STIM","Number of stimuli vary across trials");

  fprintf(fout,"dt %.6f\n",dt);
  fprintf(fout,"ntf %d\n",ntf);
  for(i=0;i<ntf;i++){
    fprintf(fout,"%d\n",tflist[i]);
  }
  fprintf(fout,"ntrial %d\n",n);
  fprintf(fout,"nstim %d\n",ns[0]);
  
  for(i=0;i<n;i++){
    for(j=0;j<ns[0];j++){
      fprintf(fout,"%d %d\n",stima[i][j],stimp[i][j]);
    }
  }

  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_ADDOT_STIM_SEQ_GROUP                        */
/*                                                                           */
/*  Compute the directions of motion associated with the 'addot' stimuli.    */
/*  These are adaptation dot stimuli developed for Adam Kohn.                */
/*                                                                           */
/*  Return 2d array rstim[trial][stimulus] which has the orientations of     */
/*  the stimuli listed in the order in which they occurred for the first     */
/*  and second test periods.                                                 */
/*                                                                           */
/*****************************************************************************/
void get_addot_stim_seq_group(nd,ndg,k,rstim,rnstim)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k,***rstim,*rnstim;
{
  int i,j,l;
  int n,**stim,nstim,flag,tseed,seed,nstim1,nstim2,ori1,ori2,dotseed;
  int *t1shuf,*t2shuf;
  char *pfile;

  n = ndg->cnt[k]; /* Number of trials in group */

  flag = ndata_get_const_param(nd,"pfile",&pfile);
  if (flag == 0) exit_error("GET_ADDOT_STIM_SEQ_GROUP","'pfile' not found");

  /*** Look at first trial in group to determine number of stimuli ***/
  get_addot_stim_params(nd,ndg->tnum[k][0],&seed,&nstim1,&nstim2,&ori1,&ori2);
  printf("  nstim1,2= %d %d  ori1,2= %d %d\n",nstim1,nstim2,ori1,ori2);

  nstim = nstim1 + nstim2;
  stim = get_2d_iarray(n,nstim);

  for(i=0;i<n;i++){ /*** For each trial ***/
    get_addot_stim_params(nd,ndg->tnum[k][i],&seed,&nstim1,&nstim2,&ori1,
			  &ori2);
    if (seed > 0)
      tseed = -seed;
    else
      tseed = seed;
    
    t1shuf = get_shuffle_index_return_seed(nstim1,3,&tseed);
    dotseed = -tseed; /* used for a different random generator for dots */
    printf("  dotseed = %d\n",dotseed);
    l = 0;
    for(j=0;j<nstim1;j++){
      stim[i][l] = ori1 + 360.0 * (float)t1shuf[j]/(float)nstim1;
      l += 1;
    }
    myfree(t1shuf);
    
    t2shuf = get_shuffle_index_return_seed(nstim2,3,&tseed);
    for(j=0;j<nstim2;j++){
      stim[i][l] = ori2 + 360.0 * (float)t2shuf[j]/(float)nstim2;
      l += 1;
    }
    myfree(t2shuf);
  }

  myfree(pfile);

  *rstim = stim; *rnstim = nstim;
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_SEED40_INDEX                            */
/*                                                                           */
/*  Return the index of the "seed" in the seed40 list, -1 if not found.      */
/*                                                                           */
/*****************************************************************************/
int get_seed40_index(seed)
     int seed;
{
  int i,k;

  k = -1;
  for(i=0;i<40;i++)
    if (ssseed40[i] == seed)
      k = i;

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_MAP_STIM_PARAMS                           */
/*                                                                           */
/*  Return parameters for the "map" stimulus for the "k"th trial.            */
/*                                                                           */
/*****************************************************************************/
void get_map_stim_params(nd,k,rframerate,rseed,rndir,rtheta,rduration,rdtstim,
			 rdtisi,rntf,rtfmax,rtfratio,rnsper,rsper,rspratio,
			 rnhpat,rnvpat,rpatsep)
     struct ndata_struct *nd;
     int k;
     int *rframerate,*rseed,*rndir,*rtheta,*rduration,*rdtstim,*rdtisi;
     int *rntf,*rtfmax,*rtfratio,*rnsper,*rsper,*rspratio,*rnhpat,*rnvpat;
     int *rpatsep;
{
  int flag;

  flag = ndata_get_vc_param_int(nd,k,"framrate",rframerate);
  if (!flag) exit_error("GET_MAP_STIM_PARAMS","Param framerate not found");

  flag = ndata_get_vc_param_int(nd,k,"seed",rseed);
  if (!flag) exit_error("GET_MAP_STIM_PARAMS","Param seed not found");
  flag = ndata_get_vc_param_int(nd,k,"ndir",rndir);
  if (!flag) exit_error("GET_MAP_STIM_PARAMS","Param ndir not found");
  flag = ndata_get_vc_param_int(nd,k,"theta",rtheta);
  if (!flag) exit_error("GET_MAP_STIM_PARAMS","Param theta not found");
  flag = ndata_get_vc_param_int(nd,k,"duration",rduration);
  if (!flag) exit_error("GET_MAP_STIM_PARAMS","Param duration not found");
  flag = ndata_get_vc_param_int(nd,k,"dt",rdtstim);
  if (!flag) exit_error("GET_MAP_STIM_PARAMS","Param dt not found");
  flag = ndata_get_vc_param_int(nd,k,"dtisi",rdtisi);
  if (!flag) exit_error("GET_MAP_STIM_PARAMS","Param dtisi not found");

  flag = ndata_get_vc_param_int(nd,k,"ntf",rntf);
  if (!flag) exit_error("GET_MAP_STIM_PARAMS","Param ntf not found");
  flag = ndata_get_vc_param_int(nd,k,"tfmax",rtfmax);
  if (!flag) exit_error("GET_MAP_STIM_PARAMS","Param tfmax not found");
  flag = ndata_get_vc_param_int(nd,k,"tfratio",rtfratio);
  if (!flag) exit_error("GET_MAP_STIM_PARAMS","Param tfratio not found");

  flag = ndata_get_vc_param_int(nd,k,"nsper",rnsper);
  if (!flag) exit_error("GET_MAP_STIM_PARAMS","Param nsper not found");
  flag = ndata_get_vc_param_int(nd,k,"sper",rsper);
  if (!flag) exit_error("GET_MAP_STIM_PARAMS","Param sper not found");
  flag = ndata_get_vc_param_int(nd,k,"spratio",rspratio);
  if (!flag) exit_error("GET_MAP_STIM_PARAMS","Param spratio not found");

  flag = ndata_get_vc_param_int(nd,k,"nhpat",rnhpat);
  if (!flag) exit_error("GET_MAP_STIM_PARAMS","Param nhpat not found");
  flag = ndata_get_vc_param_int(nd,k,"nvpat",rnvpat);
  if (!flag) exit_error("GET_MAP_STIM_PARAMS","Param nvpat not found");
  flag = ndata_get_vc_param_int(nd,k,"patsep",rpatsep);
  if (!flag) exit_error("GET_MAP_STIM_PARAMS","Param patsep not found");
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_MAP_STIM_FOR_TRIAL                         */
/*                                                                           */
/*  Return the stimulus sequence for the "k"the trial.  This routine covers  */
/*  the following p-files:                                                   */
/*    wymap10.p                                                              */
/*                                                                           */
/*  "ddir"  direction                                                        */
/*  "dtf"   temporal frequency                                               */
/*  "dsper" spatial period                                                   */
/*  "dhloc" horizontal location                                              */
/*  "dvloc" vertical location                                                */
/*                                                                           */
/*  *** WARNING:  "nshuffle" is hard-coded to 3.                             */
/*                                                                           */
/*****************************************************************************/
void get_map_stim_for_trial(nd,k,rhloc,rvloc,rdir,rtf,rsper,rn,
			    rvalh,rvalv,rvaldir,rvaltf,rvalsper,
			    rnh,rnv,rndir,rntf,rnsper,rdtstim,rdtisi)
     struct ndata_struct *nd;
     int k;
     int **rhloc,**rvloc,**rdir,**rtf,**rsper,*rn;
     int **rvalh,**rvalv,**rvaldir,**rvaltf,**rvalsper;
     int *rnh,*rnv,*rndir,*rntf,*rnsper,*rdtstim,*rdtisi;
{
  int i;
  int framerate,seed,ndir,theta,duration,dtstim,dtisi;
  int ntf,tfmax,tfratio,nsper,sper,spratio,nhpat,nvpat,patsep;
  int *ddir,*dtf,*dsper,*dhloc,*dvloc;
  int *valdir,*valtf,*valsper,*valx,*valy;
  int n,*shuffle,nshuffle,pos,npos,nother;

  get_map_stim_params(nd,k,&framerate,&seed,&ndir,&theta,&duration,&dtstim,
		      &dtisi,&ntf,&tfmax,&tfratio,&nsper,&sper,&spratio,&nhpat,
		      &nvpat,&patsep);

  npos = nhpat * nvpat;
  nother = ndir * ntf * nsper;
  n = npos * nother;
  nshuffle = 3;
  shuffle = get_shuffle_index(n,nshuffle,seed);

  ddir = (int *)myalloc(n*sizeof(int));
  dtf = (int *)myalloc(n*sizeof(int));
  dsper = (int *)myalloc(n*sizeof(int));
  dhloc = (int *)myalloc(n*sizeof(int));
  dvloc = (int *)myalloc(n*sizeof(int));

  for(i=0;i<n;i++){
    k = shuffle[i];
    dsper[i] = k % nsper;
    pos = (int)(k/nother);
    dhloc[i] = pos % nhpat;
    dvloc[i] = (int)(pos/nvpat);

    ddir[i] = (int)(k % nother/(ntf*2)); /*** VDP number, ndir/2 ***/
    dtf[i] = (int)((k % (nsper * ntf * 2))/(nsper)); /*** TF and direction ***/
    if (dtf[i] >= ntf){
      ddir[i] = ddir[i] + ndir/2;
      dtf[i] = dtf[i] - ntf;
    }
  }

  valdir = (int *)myalloc(ndir*sizeof(int));
  valtf = (int *)myalloc(ntf*sizeof(int));
  valsper = (int *)myalloc(nsper*sizeof(int));
  valx = (int *)myalloc(nhpat*sizeof(int));
  valy = (int *)myalloc(nvpat*sizeof(int));
  for(i=0;i<ndir;i++)
    valdir[i] = theta + (float)i/(float)ndir*360.0;
  valtf[0] = tfmax;
  for(i=1;i<ntf;i++)
    valtf[i] = valtf[i-1]/tfratio;
  valsper[0] = sper;
  for(i=1;i<nsper;i++)
    valsper[i] = valsper[i-1]/spratio;
  for(i=0;i<nhpat;i++)
    valx[i] = i*patsep - (int)(((float)nhpat - 1.0)/2.0 * patsep);
  for(i=0;i<nvpat;i++)
    valy[i] = i*patsep - (int)(((float)nvpat - 1.0)/2.0 * patsep);

  /***
    for(i=0;i<n;i++){
    dhloc[i] = valx[dhloc[i]];
    dvloc[i] = valy[dvloc[i]];
    ddir[i] = valdir[ddir[i]];
    dtf[i] = valtf[dtf[i]];
    dsper[i] = valsper[dsper[i]];
    }***/

  *rhloc=dhloc; *rvloc=dvloc; *rdir=ddir; *rtf=dtf; *rsper=dsper; *rn=n;
  *rvalh=valx; *rvalv=valy; *rvaldir=valdir; *rvaltf=valtf; *rvalsper=valsper;
  *rnh=nhpat; *rnv=nvpat; *rndir=ndir; *rntf=ntf; *rnsper=nsper;
  *rdtstim=dtstim; *rdtisi=dtisi;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_MAP_STIM_GROUP                            */
/*                                                                           */
/*****************************************************************************/
void get_map_stim_group(nd,ndg,k,rhloc,rvloc,rdir,rtf,rsper,rn,
			rvalh,rvalv,rvaldir,rvaltf,rvalsper,rnh,rnv,rndir,rntf,
			rnsper,rdtstim,rdtisi)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k;
     int ***rhloc,***rvloc,***rdir,***rtf,***rsper,**rn;
     int ***rvalh,***rvalv,***rvaldir,***rvaltf,***rvalsper;
     int **rnh,**rnv,**rndir,**rntf,**rnsper,**rdtstim,**rdtisi;

{
  int i;
  int flag,n;
  int **hloc,**vloc,**dir,**tf,**sper,*nstim;
  int **valh,**valv,**valdir,**valtf,**valsper;
  int *nh,*nv,*ndir,*ntf,*nsper,*dtstim,*dtisi;
  char *pfile;

  flag = ndata_get_const_param(nd,"pfile",&pfile);

  if (flag == 0)
    exit_error("GET_MAP_STIM_GROUP","Const param pfile not found");
  if ((strcmp(pfile,"wymap10.p")==0)||(strcmp(pfile,"wymap20.p")==0)||
      (strcmp(pfile,"wymap.p")==0)||(strcmp(pfile,"wyspat20.p")==0)){
    n = ndg->cnt[k];

    hloc = (int **)myalloc(n*sizeof(int *));
    vloc = (int **)myalloc(n*sizeof(int *));
    dir = (int **)myalloc(n*sizeof(int *));
    tf = (int **)myalloc(n*sizeof(int *));
    sper = (int **)myalloc(n*sizeof(int *));
    nstim = (int *)myalloc(n*sizeof(int));

    valh = (int **)myalloc(n*sizeof(int *));
    valv = (int **)myalloc(n*sizeof(int *));
    valdir = (int **)myalloc(n*sizeof(int *));
    valtf = (int **)myalloc(n*sizeof(int *));
    valsper = (int **)myalloc(n*sizeof(int *));

    nh = (int *)myalloc(n*sizeof(int));
    nv = (int *)myalloc(n*sizeof(int));
    ndir = (int *)myalloc(n*sizeof(int));
    ntf = (int *)myalloc(n*sizeof(int));
    nsper = (int *)myalloc(n*sizeof(int));
    dtstim = (int *)myalloc(n*sizeof(int));
    dtisi = (int *)myalloc(n*sizeof(int));

    for(i=0;i<n;i++)
      get_map_stim_for_trial(nd,ndg->tnum[k][i],&hloc[i],&vloc[i],&dir[i],
			     &tf[i],&sper[i],&nstim[i],&valh[i],&valv[i],
			     &valdir[i],&valtf[i],&valsper[i],&nh[i],&nv[i],
			     &ndir[i],&ntf[i],&nsper[i],&dtstim[i],&dtisi[i]);

    *rhloc=hloc; *rvloc=vloc; *rdir=dir; *rtf=tf; *rsper=sper; *rn=nstim;
    *rvalh=valh; *rvalv=valv; *rvaldir=valdir; *rvaltf=valtf;
    *rvalsper=valsper; *rnh=nh; *rnv=nv; *rndir=ndir; *rntf=ntf; *rnsper=nsper;
    *rdtstim=dtstim; *rdtisi=dtisi;
    myfree(pfile);
  }else
    exit_error("GET_MAP_STIM_GROUP","Unrecognized pfile");
}
/**************************************-**************************************/
/*                                                                           */
/*                                GET_MAP_PSTH                               */
/*                                                                           */
/*  Return the PSTHs for the specified stimulus conditions.                  */
/*                                                                           */
/*    toffset  - time stimulus begins relative to spike time 0.              */
/*    period   - period of PSTH in sampling units.                           */
/*    binsize  - for PSTH, in sampling units.                                */
/*    frate    - frame rate, with high precision                             */
/*                                                                           */
/*****************************************************************************/
void get_map_psth(s,cnt,n,hloc,vloc,dir,tf,sper,nstim,valh,valv,valdir,valtf,
		  valsper,nh,nv,ndir,ntf,nsper,dtstim,dtisi,toffset,period,
		  binsize,hpos,vpos,tdir,ttf,tsper,tframe,rpsth,rnbin)
     int **s,*cnt,n;
     int **hloc,**vloc,**dir,**tf,**sper,*nstim;
     int **valh,**valv,**valdir,**valtf,**valsper,*nh,*nv,*ndir,*ntf,*nsper;
     int *dtstim,*dtisi,toffset,period,binsize,hpos,vpos,tdir,ttf,tsper;
     float tframe,*******rpsth;
     int *rnbin;
{
  int i,j,k,l,m;
  int npsth,nbin;
  float start,******fpsth; /* [hloc][vloc][dir][tf][sper] */

  printf("  GET_MAP_PSTH\n");
  printf("    nh=%d nv=%d ndir=%d ntf=%d nsper=%d\n",nh[0],nv[0],ndir[0],
	 ntf[0],nsper[0]);

  npsth = nh[0]*nv[0]*ndir[0]*ntf[0]*nsper[0];
  nbin = period/binsize;
  printf("    Computing %d PSTHs with %d bins.\n",npsth,nbin);

  fpsth = (float ******)myalloc(nh[0]*sizeof(float *****));
  for(i=0;i<nh[0];i++){
    fpsth[i] = (float *****)myalloc(nv[0]*sizeof(float ****));
    for(j=0;j<nv[0];j++){
      fpsth[i][j] = (float ****)myalloc(ndir[0]*sizeof(float ***));
      for(k=0;k<ndir[0];k++){
	fpsth[i][j][k] = (float ***)myalloc(ntf[0]*sizeof(float **));
	for(l=0;l<ntf[0];l++){
	  fpsth[i][j][k][l] = (float **)myalloc(nsper[0]*sizeof(float *));
	  for(m=0;m<nsper[0];m++)
	    fpsth[i][j][k][l][m] = get_zero_farray(nbin);
	}
      }
    }
  }

  printf("    Trial:");
  for(i=0;i<n;i++){ /*** For each trial. ***/
    for(j=0;j<nstim[i];j++){ /*** For each stimulus. ***/
      /*printf("%d %d %d %d\n",hloc[i][j],vloc[i][j],dir[i][j],tf[i][j]);*/
      start = (float)(j*(dtstim[0]+dtisi[0])*tframe);
      increment_psth_float_bin_trial(s[i],cnt[i],start,(float)period,
				     (float)binsize,
       fpsth[hloc[i][j]][vloc[i][j]][dir[i][j]][tf[i][j]][sper[i][j]],nbin);
    }
    printf(" %d",i);
    fflush(stdout);
  }
  printf("\n");

  *rpsth = fpsth; *rnbin = nbin;
}
/**************************************-**************************************/
/*                                                                           */
/*                          NOISE_UTIL_STAR_GET_TRIG                         */
/*                                                                           */
/*  Return trigger times for 'star' stimulus.                                */
/*  Values returned are INTEGER in 'sampling' units.                         */
/*                                                                           */
/*****************************************************************************/
void noise_util_star_get_trig(nd,ndg,gi,sampling,ttype,rtrig,rntrig)
     struct ndata_struct *nd;      //
     struct ndgroup_struct *ndg;   //
     int gi;                       // Group index
     float sampling;               // for returned triggers, typ. 1000.0
     char *ttype;                  // trigger type, see STM_STAR_GET_TRIG
     int ***rtrig,**rntrig;        // triggers, counts [trial][
{
  int i;
  int n,state0,dta,dt0,dt1,*t,tn,**trig,*ntrig;
  float fps,stn;
  
  //printf("    NOISE_UTIL_STAR_GET_TRIG\n");

  fps = ndata_get_const_param_float(nd,"labrcon_fps");

  n = ndg->cnt[gi];  // Number of trials in group

  trig = (int **)myalloc(n*sizeof(int *));
  ntrig = (int *)myalloc(n*sizeof(int));

  for(i=0;i<n;i++){ // For each trial

    dta = ndata_get_vc_param_float_or_exit(nd,i,"dt_adapt",
					   "NOISE_UTIL_STAR_GET_TRIG");
    dt0 = ndata_get_vc_param_float_or_exit(nd,i,"dt_test0",
					   "NOISE_UTIL_STAR_GET_TRIG");
    dt1 = ndata_get_vc_param_float_or_exit(nd,i,"dt_test1",
					   "NOISE_UTIL_STAR_GET_TRIG");
    state0= ndata_get_vc_param_int_or_exit(nd,i,"start_state",
					   "NOISE_UTIL_STAR_GET_TRIG");
    stn = ndata_get_vc_param_int_or_exit(nd,i,"stn",
					 "NOISE_UTIL_STAR_GET_TRIG");

    //printf("  Trial %3d  state0 %d\n",i,state0);
    //printf("      dt_adpat  %d frames\n",dta);
    //printf("      dt_test0  %d frames\n",dt0);
    //printf("      dt_test1  %d frames\n",dt1);
    
    stm_star_get_trig(stn,fps,sampling,dta,dt0,dt1,state0,ttype,&t,&tn);

    printf("      Got %d triggers\n",tn);
    
    trig[i] = t;
    ntrig[i] = tn;
  }

  *rtrig = trig;
  *rntrig = ntrig;
}
