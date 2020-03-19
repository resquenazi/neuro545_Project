/*****************************************************************************/
/*                                                                           */
/*  spike_util.c                                                             */
/*  wyeth bair                                                               */
/*  caltech                                                                  */
/*  02/23/93                                                                 */
/*                                                                           */
/*  This file contains utilities to operate on spike trains.                 */
/*  Ideally, these routines should allow for negative spike times.           */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "my_util.h"
#include "nr_util.h"
#include "misc_util.h"
#include "carray_util.h"
#include "iarray_util.h"
#include "farray_util.h"
#include "plot_util.h"
#include "myrand_util.h"
#include "spike_util.h"
#include "spike_data.h"    // Only used for struct event_trial

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
/*                              SPIKE_UTIL_RAN2                              */
/*                                                                           */
/*  Long period (>2*10^18) random number generator of L'Ecuyer               */
/*  with Bays-Durham shuffle and added safeguards.  Returns a                */
/*  uniform deviate between 0.0 and 1.0 (exclusive of the                    */
/*  endpoint values).  Call with idum a negative integer to                  */
/*  initialize; thereafter, do not alter idum between                        */
/*  successive deviates in a sequence.  RNMX should approximate              */
/*  the largest floating value that is less than 1.                          */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float spike_util_ran2(idum)
     int *idum;
{
  int j;
  int k;
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
/*                                 SU_GASDEV                                 */
/*                                                                           */
/*    ********** Initialize "spike_util_ran2" before calling ***********     */
/*                                                                           */
/*  Returns a normally distributed deviate with zero mean and unit           */
/*  variance, using ran1(idum) as the source of uniform deviates.            */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float su_gasdev(idum)
     int *idum;
{
  float ran1();
  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;
  
  /*if (iset==0){*/
  if ((iset == 0) || (*idum <= 0)){  /* Wyeth, reset state if reseeded */
    do {
      v1=2.0*spike_util_ran2(idum)-1.0;
      v2=2.0*spike_util_ran2(idum)-1.0;
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
/*                               GAUSSIAN_INT                                */
/*                                                                           */
/*   Return a Gaussian deviate, N(mu,sigma), rounded to an integer.          */
/*                                                                           */
/*  Copied from "make_spikes.c" 9/13/96.                                     */
/*                                                                           */
/*****************************************************************************/
int gaussian_int(mu,sigma,seed)
     float mu,sigma;
     int *seed;
{
  return(my_rint(su_gasdev(seed)*sigma + mu));
}
/**************************************-**************************************/
/*                                                                           */
/*                                 SU_EXPDEV                                 */
/*                                                                           */
/*****************************************************************************/
float su_expdev(lambda,idum)
     float lambda;
     int *idum;
{
  if (lambda == 0.0)
    return 0.0;
  else
    return -log(spike_util_ran2(idum))/lambda;
}
/**************************************-**************************************/
/*                                                                           */
/*                             SU_GAUSSIAN_FLOAT                             */
/*                                                                           */
/*   Return a Gaussian deviate, N(mu,sigma).                                 */
/*                                                                           */
/*****************************************************************************/
float su_gaussian_float(mu,sigma,idum)
     float mu,sigma;
     int *idum;
{
  return(su_gasdev(idum)*sigma + mu);
}
/**************************************-**************************************/
/*                                                                           */
/*                              FREE_SARRAY                                  */
/*                                                                           */
/*****************************************************************************/
void free_sarray(s,cnt,n)
     int **s,*cnt,n;
{
  int i;

  for(i=0;i<n;i++)
    myfree(s[i]);
  myfree(s); myfree(cnt);
}
/**************************************-**************************************/
/*                                                                           */
/*                            WRITE_SPIKE_PLOT                               */
/*                                                                           */
/*  Output plot has x-axis units msec.                                       */
/*                                                                           */
/*****************************************************************************/
void write_spike_plot(outfile,s,n,start,period,sampling)
     char outfile[];
     int *s,n,start,period;
     float sampling;
{
  FILE *fopen(),*fout;
  int i;

  fout = fopen(outfile,"w");
  fprintf(fout,"%d 0\n",(int)((float)start*1000.0/sampling));
  for (i=0;i<n;i++)
    if ((s[i] >= start) && (s[i] < start+period)){
      fprintf(fout,"%d 0\n",(int)((float)s[i]*1000.0/sampling));
      fprintf(fout,"%d 1\n",(int)((float)s[i]*1000.0/sampling));
      fprintf(fout,"%d 0\n",(int)((float)s[i]*1000.0/sampling));
    }
  fprintf(fout,"%d 0\n",(int)((float)(start+period-1)*1000.0/sampling));
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                         WRITE_FLOAT_SPIKE_PLOT                            */
/*                                                                           */
/*  Output plot has x-axis units msec.                                       */
/*                                                                           */
/*****************************************************************************/
void write_float_spike_plot(outfile,s,n,start,period,sampling)
     char outfile[];
     float *s;
     int n,start,period;
     float sampling;
{
  FILE *fopen(),*fout;
  int i;

  fout = fopen(outfile,"w");
  fprintf(fout,"%d 0\n",(int)((float)start*1000.0/sampling));
  for (i=0;i<n;i++)
    if ((s[i] >= (float)start) && (s[i] < (float)(start+period))){
      fprintf(fout,"%.3e 0\n",s[i]*1000.0/sampling);
      fprintf(fout,"%.3e 1\n",s[i]*1000.0/sampling);
      fprintf(fout,"%.3e 0\n",s[i]*1000.0/sampling);
    }
  fprintf(fout,"%d 0\n",(int)((float)(start+period-1)*1000.0/sampling));
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                                COUNT_SPIKES                               */
/*                                                                           */
/*  Count and return the number of spikes in [start,start+period-1].         */
/*                                                                           */
/*****************************************************************************/
int count_spikes(data,n,start,period)
     int *data,n;
     int start,period;
{
  int i;
  int scount;

  scount = 0;
  for (i=0;i<n;i++)
    if ((data[i] >= start) && (data[i] < start+period))
      scount += 1;

  return scount;
}
/**************************************-**************************************/
/*                                                                           */
/*                             COUNT_SPIKES_FLOAT                            */
/*                                                                           */
/*  Count and return the number of spikes in [start,start+period-1].         */
/*                                                                           */
/*****************************************************************************/
int count_spikes_float(data,n,start,period)
     float *data;
     int n;
     float start,period;
{
  int i;
  int scount;

  scount = 0;
  for (i=0;i<n;i++)
    if ((data[i] >= start) && (data[i] < start+period))
      scount += 1;

  return scount;
}
/**************************************-**************************************/
/*                                                                           */
/*                              COUNT_SPIKES_DBL                             */
/*                                                                           */
/*  Count and return the number of spikes in [start,start+period-1].         */
/*                                                                           */
/*****************************************************************************/
int count_spikes_dbl(data,n,start,period)
     double *data;
     int n;
     double start,period;
{
  int i;
  int scount;

  scount = 0;
  for (i=0;i<n;i++)
    if ((data[i] >= start) && (data[i] < start+period))
      scount += 1;

  return scount;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MEAN_SPIKE_RATE                              */
/*                                                                           */
/*  Return the mean spike rate, in spikes/sec, for the spike train segment.  */
/*                                                                           */
/*****************************************************************************/
float mean_spike_rate(data,n,start,period,sampling)
     int *data,n;
     int start,period;
     float sampling;
{
  float sr;

  sr = (float)count_spikes(data,n,start,period);
  sr /= (float)period/sampling;

  return sr;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MEAN_SPIKE_RATE_FLOAT                          */
/*                                                                           */
/*  Return the mean spike rate, in spikes/sec, for the spike train segment.  */
/*                                                                           */
/*****************************************************************************/
float mean_spike_rate_float(data,n,start,period,sampling)
     float *data;
     int n;
     float start,period;
     float sampling;
{
  float sr;

  sr = (float)count_spikes_float(data,n,start,period);
  sr /= period/sampling;

  return sr;
}
/**************************************-**************************************/
/*                                                                           */
/*                         TOTAL_SPIKES_SARRAY                               */
/*                                                                           */
/*  Count and return the number of spikes in [start,start+period-1].         */
/*                                                                           */
/*****************************************************************************/
int total_spikes_sarray(s,cnt,n,start,period)
     int **s,*cnt,n,start,period;
{
  int i;
  int total;

  total = 0;
  for (i=0;i<n;i++)
    total += count_spikes(s[i],cnt[i],start,period);
  return total;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_MAX_TIME_SARRAY                           */
/*                                                                           */
/*  Return the largest spike time value in sarray.                           */
/*                                                                           */
/*****************************************************************************/
int get_max_time_sarray(s,cnt,n)
     int **s,*cnt,n;
{
  int i;
  int max,tmax;

  max = max_of_iarray(s[0],cnt[0]);
  for(i=1;i<n;i++){
    tmax = max_of_iarray(s[i],cnt[i]);
    if (tmax > max)
      max = tmax;
  }
  return max;
}
/**************************************-**************************************/
/*                                                                           */
/*                        CONVERT_SAMPLING_UNITS_TRIAL                       */
/*                                                                           */
/*  Return a spike train that has the sampling units converted by            */
/*  multiplying by a factor "f" and rounding to the nearest integer.         */
/*                                                                           */
/*****************************************************************************/
int *convert_sampling_units_trial(s,n,f)
     int *s,n;
     float f;
{
  int i;
  int *r;
  float sf;

  r = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++){
    //*** WYETH - fixed a possible problem here w/ negative vs. pos. times.
    //*** WYETH    on Jun 7, 2016  made this match  my_rint(...)
    //r[i] = (int)(0.5 + ((float)s[i]*f)); // *** WYETH: may mess up neg. times

    sf = (float)s[i]*f;
    if (sf >= 0.0)
      r[i] = (int)(0.5 + sf);
    else
      r[i] = (int)(-0.5 + sf);
  }

  return r;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_SARRAY_CONVERT_SAMPLING                       */
/*                                                                           */
/*  Return a spike array that has the sampling units converted by            */
/*  multiplying by a factor "f" and rounding to the nearest integer.         */
/*                                                                           */
/*****************************************************************************/
int **get_sarray_convert_sampling(s,cnt,n,f)
     int **s,*cnt,n;
     float f;
{
  int i;
  int **r;

  r = (int **)myalloc(n*sizeof(int *));
  for(i=0;i<n;i++)
    r[i] = convert_sampling_units_trial(s[i],cnt[i],f);
  return r;
}
/**************************************-**************************************/
/*                                                                           */
/*                           SARRAY_CONVERT_SAMPLING                         */
/*                                                                           */
/*  Convert the sampling units of the given sarray by multiplying by a       */
/*  factor "f" and rounding to the nearest integer.                          */
/*                                                                           */
/*****************************************************************************/
void sarray_convert_sampling(s,cnt,n,f)
     int **s,*cnt,n;
     float f;
{
  int i,j;
  float t;

  for(i=0;i<n;i++)
    for(j=0;j<cnt[i];j++){
      t = (float)s[i][j];
      if (t >= 0.0)
	s[i][j] = (int)(0.5 + t*f);
      else
	s[i][j] = (int)(-0.5 + t*f);
    }
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_FIRST_SPIKE                             */
/*                                                                           */
/*  Return the first spike at or after time "t", return a value less than t  */
/*  if there is no such spike.                                               */
/*                                                                           */
/*****************************************************************************/
int get_first_spike(data,n,t)
     int *data,n,t;
{
  int i;
  int first;

  first = t-1;
  i = 0;
  while ((first < t)&&(i<n)){
    if (data[i] >= t)
      first = data[i];
    else
      i += 1;
  }
  return first;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_FIRST_SPIKE_IN_WINDOW                      */
/*                                                                           */
/*  Return the time of first spike in the window.  Return a value less than  */
/*  t if there is no such spike.                                             */
/*                                                                           */
/*****************************************************************************/
int get_first_spike_in_window(data,n,t,tn)
     int *data,n,t,tn;
{
  int i,k;
  int first;

  first = t-1;
  i = 0;
  while ((first < t)&&(i<n)){
    k = data[i] - t;
    if ((k >= 0)&&(k < tn))
      first = data[i];
    else
      i += 1;
  }
  return first;
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_LAST_SPIKE                              */
/*                                                                           */
/*  Return the last spike at or before time "t", return a value              */
/*  greater than t if there is no such spike.                                */
/*                                                                           */
/*****************************************************************************/
int get_last_spike(data,n,t)
     int *data,n,t;
{
  int i;
  int last;

  last = t+1;
  i = n-1;
  while ((last > t)&&(i>=0)){
    if (data[i] <= t)
      last = data[i];
    else
      i -= 1;
  }
  return last;
}
/**************************************-**************************************/
/*                                                                           */
/*                               EXTRACT_SPIKES                              */
/*                                                                           */
/*   Return a spike time array that contains only those spikes within the    */
/*   specified window.                                                       */
/*                                                                           */
/*****************************************************************************/
int *extract_spikes(spike_data,n,start,period,count,zero_flag)
     int *spike_data,n,start,period;
     int *count;
     int zero_flag;  // if 1, subtract "start" from spike times
{
  int i,k,t;
  int *data;
  
  *count = count_spikes(spike_data,n,start,period);
  data = (int *)myalloc(*count*sizeof(int));
  k = 0;
  for (i=0;i<n;i++){
    t = spike_data[i] - start;
    if ((t >= 0) && (t < period)){
      if (zero_flag)
	data[k] = t;
      else
	data[k] = spike_data[i];
      k += 1;
    }
  }
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                        EXTRACT_INT_SPIKES_FROM_FLOAT                      */
/*                                                                           */
/*   Return a spike time array that contains only those spikes within the    */
/*   specified window.                                                       */
/*                                                                           */
/*****************************************************************************/
int *extract_int_spikes_from_float(s,n,start,period,count,zero_flag,tscale)
     float *s;
     int n;
     float start,period;
     int *count;
     int zero_flag;  // if 1, subtract "start" from spike times
     float tscale;   // Multiply by this factor before rounding to int
{
  int i,k;
  int *data;
  
  *count = count_spikes_float(s,n,start,period);
  data = (int *)myalloc(*count*sizeof(int));

  k = 0;
  for (i=0;i<n;i++){
    // Use same condition here as in 'count_spikes_float'
    if ((s[i] >= start) && (s[i] < start+period)){
      if (zero_flag)
	data[k] = my_rint(tscale*(s[i] - start));
      else
	data[k] = my_rint(tscale*s[i]);
      k += 1;
    }
  }
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                         EXTRACT_INT_SPIKES_FROM_DBL                       */
/*                                                                           */
/*   Return a spike time array that contains only those spikes within the    */
/*   specified window.                                                       */
/*                                                                           */
/*****************************************************************************/
int *extract_int_spikes_from_dbl(s,n,start,period,count,zero_flag,tscale)
     double *s;
     int n;
     double start,period;
     int *count;
     int zero_flag;  // if 1, subtract "start" from spike times
     double tscale;   // Multiply by this factor before rounding to int
{
  int i,k;
  int *data;

  *count = count_spikes_dbl(s,n,start,period);
  data = (int *)myalloc(*count*sizeof(int));

  k = 0;
  for (i=0;i<n;i++){
    // Use same condition here as in 'count_spikes_float'
    if ((s[i] >= start) && (s[i] < start+period)){
      if (zero_flag)
	data[k] = my_rint(tscale*(s[i] - start));
      else
	data[k] = my_rint(tscale*s[i]);
      k += 1;
    }
  }

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                               ADD_CONST_SARRAY                            */
/*                                                                           */
/*  Add a constant value to all times in the sarray.                         */
/*                                                                           */
/*****************************************************************************/
void add_const_sarray(s,cnt,n,c)
     int **s,*cnt,n,c;
{
  int i,j;

  for(i=0;i<n;i++)
    for(j=0;j<cnt[i];j++)
      s[i][j] += c;
}
/**************************************-**************************************/
/*                                                                           */
/*                            COUNT_UNIQUE_SPIKES                            */
/*                                                                           */
/*  Count and return the number of unique spike times in                     */
/*  [start,start+duration-1].                                                */
/*                                                                           */
/*  *** NOTE:  This routine assumes spike times are non-decreasing.          */
/*                                                                           */
/*****************************************************************************/
int count_unique_spikes(data,n,start,duration)
     int *data,n;
     int  start,duration;
{
  int i;
  int scount;
  int t,t_old;
  
  scount = 0;
  t_old = -1;
  for (i=0;i<n;i++){
    t = data[i] - start;
    if ((t >= 0)&&(t < duration))
      if (t != t_old)
	scount += 1;
    t_old = t;
  }
  return scount;
}
/**************************************-**************************************/
/*                                                                           */
/*                            EXTRACT_UNIQUE_SPIKES                          */
/*                                                                           */
/*  Return the spike train with only the unique spike times.                 */
/*                                                                           */
/*  *** NOTE:  This routine assumes spike times are non-decreasing.          */
/*                                                                           */
/*****************************************************************************/
int *extract_unique_spikes(data,n,start,duration,n_unique,zero_flag)
     int *data,n;
     int start,duration;
     int *n_unique;
     int zero_flag;  /* if 1, subtract "start" from spike times */
{
  int i;
  int t,t_old,pos,*udata;

  *n_unique = count_unique_spikes(data,n,start,duration);
  udata = (int *)myalloc(*n_unique*sizeof(int));
  
  pos = 0;
  t_old = -1;
  for(i=0;i<n;i++){
    t = data[i] - start;
    if ((t >= 0)&&(t < duration))
      if (t != t_old){
	if (zero_flag)
	  udata[pos] = t;
	else
	  udata[pos] = data[i];
	pos += 1;
      }
    t_old = t;
  }
  return udata;
}	
/**************************************-**************************************/
/*                                                                           */
/*                              GET_UNIQUE_SARRAY                            */
/*                                                                           */
/*  Make an sarray which does not have duplicate spike times.                */
/*                                                                           */
/*****************************************************************************/
void get_unique_sarray(s,cnt,n,period,rs,rcnt)
     int **s,*cnt,n,period,***rs,**rcnt;
{
  int i;
  int **ts,*tcnt;

  ts = (int **)myalloc(n*sizeof(int *));
  tcnt = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++)
    ts[i] = extract_unique_spikes(s[i],cnt[i],0,period,&tcnt[i],0);

  *rs = ts; *rcnt = tcnt;
}	
/**************************************-**************************************/
/*                                                                           */
/*                        SPIKEU_EXPAND_SPIKES_WSD_ACCUM                     */
/*                                                                           */
/*  Expand the float spike train using the weight and synaptic depression    */
/*  parameters given.                                                        */
/*                                                                           */
/*****************************************************************************/
void spikeu_expand_spikes_wsd_accum(s,n,xacc,tn,w,sdf,sdtau)
     float *s;      // Spike train [n]
     int n;         // spike count
     float *xacc;   // Accumulation array [tn]
     int tn;        // length of expanded array
     float w;       // spike weight
     float sdf;     // Synaptic depression, frac reduction after spike
     float sdtau;   // Synaptic depression, recovery time const
{
  int i;
  int t,tbad,tmax;
  float a,anew,dt;

  tbad = tmax = 0;
  
  a = 1.0;      // Amplitude is full
  dt = 10000.0; // Last spike was a long time ago (10 sec)
  for(i=0;i<n;i++){
    t = my_rint(s[i]);
    if ((t >= 0)&&(t < tn)){

      // synaptic depression
      if (i > 0)
	dt = s[i] - s[i-1];
      anew = 1.0 - (1.0 - a*sdf) * exp(-dt/sdtau);
      a = anew;

      xacc[t] += anew * w;        // Use 'w' to weight each spike

    }else{
      if (t == tn)
	tmax += 1;
      else{
	tbad += 1;
	// WYETH FIX THE CAUSE OF THIS
	// WYETH FIX THE CAUSE OF THIS
	//printf("  Spike out of range:  %d  (%f)\n",t,s[i]);
      }
    }
  }
  
  // This happens frequently, don't report it
  // if (tmax > 0)
  //  printf("    %d spike times were too large by 1 time unit.\n",tmax);

  // WYETH FIX THE CAUSE OF THIS
  // WYETH FIX THE CAUSE OF THIS

  //if (tbad > 0)
  //printf("*** %d spike times were out of range.\n",tbad);
}
/**************************************-**************************************/
/*                                                                           */
/*                             EXPAND_SPIKE_ARRAY                            */
/*                                                                           */
/*   Change part or all of an array of spike times into an array             */
/*   over time indicating the number of spikes which ocurred at each         */
/*   time.  NOTE: if the spike array has non-unique entries, the             */
/*   expanded array will reflect this by having counts equal to the          */
/*   number of occurrences of a particular spike time.                       */
/*                                                                           */
/*****************************************************************************/
float *expand_spike_array(spike_data,n,start,period)
     int *spike_data,n;
     int start,period;
{
  int i,t;
  float *data;

  data = get_zero_farray(period);
  for (i=0;i<n;i++){
    t = spike_data[i] - start;
    if ((t >= 0) && (t < period))
      data[t] += 1.0;
  }
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            EXPAND_SPIKE_INT_ARRAY                         */
/*                                                                           */
/*   Change part or all of an array of spike times into an array             */
/*   over time indicating the number of spikes which ocurred at each         */
/*   time.  NOTE: if the spike array has non-unique entries, the             */
/*   expanded array will reflect this by having counts equal to the          */
/*   number of occurrences of a particular spike time.                       */
/*                                                                           */
/*****************************************************************************/
int *expand_spike_int_array(spike_data,n,start,period)
     int *spike_data,n;
     int start,period;
{
  int i,t;
  int *data;

  data = get_zero_iarray(period);
  for (i=0;i<n;i++){
    t = spike_data[i] - start;
    if ((t >= 0) && (t < period))
      data[t] += 1;
  }
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             EXPAND_SPIKE_FARRAY                           */
/*                                                                           */
/*   Change part or all of an array of spike times into an array             */
/*   over time indicating the number of spikes which ocurred at each         */
/*   time.  NOTE: if the spike array has non-unique entries, the             */
/*   expanded array will reflect this by having counts equal to the          */
/*   number of occurrences of a particular spike time.                       */
/*                                                                           */
/*****************************************************************************/
float *expand_spike_farray(spike_data,n,start,period)
     float *spike_data;
     int n,start,period;
{
  int     i,t;
  float  *data;

  data = get_zero_farray(period);
  for (i=0;i<n;i++){
    t = my_rint(spike_data[i]) - start;
    if ((t >= 0) && (t < period))
      data[t] += 1.0;
  }
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                     GET_TIMES_FROM_EXPANDED_INT_SPIKES                    */
/*                                                                           */
/*  Return an array with the times of "1"s in the data array.                */
/*                                                                           */
/*****************************************************************************/
void get_times_from_expanded_int_spikes(data,n,start,rs,rcnt)
     int *data,n,start,**rs,*rcnt;
{
  int i,j,k;
  int *t;

  if (min_of_iarray(data,n) < 0)
    exit_error("GET_TIMES_FROM_EXPANDED_INT_SPIKES","NEGATIVE VALUE FOUND");

  *rcnt = sum_iarray(data,n,0,n);
  t = (int *)myalloc(*rcnt*sizeof(int));
  k = 0;
  for (i=0;i<n;i++)
    for(j=0;j<data[i];j++){
      t[k] = i + start;
      k += 1;
    }
  *rs = t;
}
/**************************************-**************************************/
/*                                                                           */
/*                    GET_TIMES_FROM_EXPANDED_FLOAT_SPIKES                   */
/*                                                                           */
/*  Return an array with the times of "1"s in the data array.                */
/*                                                                           */
/*****************************************************************************/
void get_times_from_expanded_float_spikes(data,n,start,s,cnt)
     float *data;
     int n,start,**s,*cnt;
{
  int i,j,k;
  int *t;

  if (min_of_farray(data,n) < 0)
    exit_error("GET_TIMES_FROM_EXPANDED_INT_SPIKES","NEGATIVE VALUE FOUND");

  *cnt = sum_farray(data,n,0,n);
  t = (int *)myalloc(*cnt*sizeof(int));
  k = 0;
  for (i=0;i<n;i++)
    for(j=0;j<(int)(data[i]);j++){
      t[k] = i + start;
      k += 1;
    }
  *s = t;
}
/**************************************-**************************************/
/*                                                                           */
/*                          HISTOGRAM_EXPAND_SARRAY                          */
/*                                                                           */
/*  Recode the spike train by combining 'binsize' time units into one time   */
/*  unit.                                                                    */
/*         rn - returns the number of bins in the recoded train              */
/*     rtrunc - returns the number of bins truncated                         */
/*     rshift - returns the number of spikes shifted to avoid truncation     */
/*                                                                           */
/*****************************************************************************/
int *histogram_expand_sarray(s,n,start,period,binsize,binflag,shiftflag,
			     rn,rtrunc,rshift)
     int *s,n;
     int start,period;
     int binsize;       /* Number of time units to include in each bin */
     int binflag;       /* 1-result is forced to be binary */
     int shiftflag;     /* Shift spikes to neighboring bins to ensure */
     int *rn,*rtrunc,*rshift;
{
  int i,k;
  int hn,t,ntrunc,nshift,*x;

  hn = period / binsize; /* Number of bins in expanded array */
  x = get_zero_iarray(hn);

  for (i=0;i<n;i++){
    t = s[i] - start;
    if ((t >= 0) && (t < period)){
      k = (int)(t/binsize);
      if ((k >= 0)&&(k < hn)){
	x[k] += 1;
      }
    }
  }
  
  ntrunc = 0;
  nshift = 0;
  if (binflag == 1){
    if (shiftflag == 0){       /*** Binarize, do not try to shift spieks ***/
      for(i=0;i<hn;i++)
	if (x[i] > 1){
	  x[i] = 1;
	  ntrunc += 1;
	}
    }else{
      for(i=0;i<hn;i++)
	if (x[i] > 1){
	  if (i > 0)           /*** Try to move 1 to the left ***/
	    if (x[i-1] == 0){
	      x[i-1] += 1;
	      x[i] -= 1;
	      nshift += 1;
	    }
	  if (x[i] > 1)        /*** If still too big ***/
	    if (i < (hn-1))         /*** Try to move it to the right ***/
	      if (x[i+1] == 0){
		x[i+1] += 1;
		x[i] -= 1;
		nshift += 1;
	      }
	  if (x[i] > 1){        /*** If still > 1, set to 1 ***/
	    x[i] = 1;
	    ntrunc += 1;
	  }
	}
    }

    for(i=0;i<hn;i++)
      if (x[i] > 1){
	exit_error("HISTOGRAM_EXPAND_SARRAY","Truncation failed");
      }
    
    /*
      printf("  %d bins truncated\n",ntrunc);
      if (shiftflag == 1)
      printf("  %d spikes shifted\n",nshift);
      */
  }

  *rn = hn; *rtrunc = ntrunc; *rshift = nshift;
  return x;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MERGE_SPIKE_ARRAYS_FLOAT                        */
/*                                                                           */
/*  Merge all spikes and sort the merged array.                              */
/*                                                                           */
/*****************************************************************************/
float *merge_spike_arrays_float(s1,n1,s2,n2)
     float *s1;
     int n1;
     float *s2;
     int n2;
{
  int i;
  int n;
  float *s;

  n = n1 + n2;
  s = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n1;i++)
    s[i] = s1[i];
  for(i=0;i<n2;i++)
    s[n1+i] = s2[i];

  if (n > 1)
    sort_farray(s,n);

  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MERGE_SPIKE_ARRAYS                            */
/*                                                                           */
/*****************************************************************************/
int *merge_spike_arrays(s1,n1,s2,n2,start,period,rn)
     int *s1,n1,*s2,n2;
     int start,period;
     int *rn;
{
  int i,j;
  int t,pos,n;
  int *merged,*expanded;

  expanded = get_zero_iarray(period);

  n = 0;
  for (i=0;i<n1;i++){ /*** MARK ALL SPIKES IN EXPANDED ARRAYS ***/
    t = s1[i]-start;
    if ((t>=0) && (t<period)){
      expanded[t] += 1;
      n += 1;
    }
  }
  for (i=0;i<n2;i++){
    t = s2[i]-start;
    if ((t>=0) && (t<period)){
      expanded[t] += 1;
      n += 1;
    }
  }

  merged = get_zero_iarray(n); /*** Make the merged point array ***/
  pos = 0;
  for (i=0;i<period;i++)
    for (j=0;j<expanded[i];j++){ /* Maybe multiple spikes at one point */
      merged[pos] = i+start;
      pos += 1;
    }
  if (pos != n)
    exit_error("MERGE_SPIKE_ARRAY","This never happens");
  
  myfree(expanded);
  *rn = n;
  return merged;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MERGE_SARRAY                                 */
/*                                                                           */
/*  Make one spike train from all spike trains in sarray.                    */
/*                                                                           */
/*  *** Apparently, spike times are not changed, but only spikes in          */
/*      the window [start,start+period-1] are included.                      */
/*                                                                           */
/*****************************************************************************/
int *merge_sarray(s,cnt,n,start,period,rn)
     int **s,*cnt,n;
     int start,period;
     int *rn;
{
  int i,j;
  int t,pos,count;
  int *merged,*expanded;

  expanded = get_zero_iarray(period);

  count = 0;
  for(i=0;i<n;i++)
    for (j=0;j<cnt[i];j++){ /*** MARK ALL SPIKES IN EXPANDED ARRAYS ***/
      t = s[i][j]-start;
      if ((t>=0) && (t<period)){
	expanded[t] += 1;
	count += 1;
      }
    }

  merged = get_zero_iarray(count); /*** Make the merged point array ***/
  pos = 0;
  for (i=0;i<period;i++)
    for (j=0;j<expanded[i];j++){ /* Maybe multiple spikes at one point */
      merged[pos] = i+start;
      pos += 1;
    }
  if (pos != count)
    exit_error("MERGE_SPIKE_ARRAY","This never happens");

  myfree(expanded);
  *rn = count;
  return merged;
}
/**************************************-**************************************/
/*                                                                           */
/*                          FLOAT_TO_INT_SPIKE                               */
/*                                                                           */
/*   Convert an array of floating point spike times to integer               */
/*   spike times using the conversion factor c.                              */
/*                                                                           */
/*****************************************************************************/
int *float_to_int_spike(fdata,fn,n,start,duration,c)
     float *fdata;
     int fn,*n;
     int start,duration; /* in sampling units */
     float c; /* to convert to sampling units */
{
  int i;
  int *s,t;

  s = (int *)myalloc(fn*sizeof(int));

  *n = 0;
  for (i=0;i<fn;i++){
    t = (int)(0.5 + c*fdata[i]);
    if ((t>=start)&&(t<(start+duration))){
      s[*n] = t;
      *n += 1;
    }
  }
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                                PRINT_SARRAY                               */
/*                                                                           */
/*  Print out the contents of the sarray in a simple format.                 */
/*                                                                           */
/*****************************************************************************/
void print_sarray(s,cnt,n,t,tn)
     int **s,*cnt,n,t,tn;
{
  int i,j;
  int *ss,sn;

  printf("\n%d\n",n);
  for(i=0;i<n;i++){
    ss = extract_spikes(s[i],cnt[i],t,tn,&sn,0);
    printf("%d\n",sn);
    for(j=0;j<sn;j++)
      printf(" %d",ss[j]);
    printf("\n");
    myfree(ss);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               COMPARE_TRAINS                              */
/*                                                                           */
/*  Compare the spike trains, return a string that describes them.           */
/*                                                                           */
/*****************************************************************************/
char *compare_trains(s1,n1,s2,n2)
     int *s1,n1;   // s1[n1]  spike times
     int *s2,n2;   // s2[n2]  spike times
{
  int i,j;
  char *t;

  t = NULL;

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                             PRINT_SARRAY_APPEND                           */
/*                                                                           */
/*  Append spike counts and times to a file in a simple text format.         */
/*                                                                           */
/*****************************************************************************/
void print_sarray_append(outfile,s,cnt,n,t,tn)
     char outfile[];
     int **s,*cnt,n,t,tn;
{
  int i,j;
  int *ss,sn;
  FILE *fout,*fopen();

  if ((fout = fopen(outfile,"a"))==NULL){
    exit_error("PRINT_SARRAY_APPEND","Cannot open file");
  }
  
  fprintf(fout,"\n%d\n",n);
  for(i=0;i<n;i++){
    ss = extract_spikes(s[i],cnt[i],t,tn,&sn,0);
    fprintf(fout,"%d\n",sn);
    for(j=0;j<sn;j++)
      fprintf(fout," %d",ss[j]);
    fprintf(fout,"\n");
    myfree(ss);
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                    ALIGN_TO_FIRST_SPIKE_IN_WINDOW_SARRAY                  */
/*                                                                           */
/*  This is used to align the responses in "discrete events" to measure      */
/*  the within-event firing modulations, under the assumption that there     */
/*  is potentially an overall jitter moving the event around.                */
/*                                                                           */
/*  Make the first spike in the window be at "tnew", and move all other      */
/*  spikes accordingly.  If there is no spike in the window, do not          */
/*  modify the spike train.                                                  */
/*                                                                           */
/*****************************************************************************/
void align_to_first_spike_in_window_sarray(s,cnt,n,t,tn,tnew)
     int **s,*cnt,n,t,tn,tnew;
{
  int i,k;

  for(i=0;i<n;i++){
    k = get_first_spike_in_window(s[i],cnt[i],t,tn);
    if (k > t)
      add_const_iarray(s[i],cnt[i],tnew-k);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             SPIKE_COUNT_SARRAYS                           */
/*                                                                           */
/*  Return the farrays having the spike counts for the two sarrays.          */
/*                                                                           */
/*****************************************************************************/
void spike_count_sarrays(s1,cnt1,s2,cnt2,n,start,period,rsc1,rsc2)
     int **s1,*cnt1,**s2,*cnt2,n;
     int start,period;
     float **rsc1,**rsc2;
{
  int i;
  float *c1,*c2,z,mean,sd1,sd2;

  c1 = (float *)myalloc(n*sizeof(float));
  c2 = (float *)myalloc(n*sizeof(float));

  for (i=0;i<n;i++){
    c1[i] = (float)count_spikes(s1[i],cnt1[i],start,period);
    c2[i] = (float)count_spikes(s2[i],cnt2[i],start,period);
  }

  *rsc1 = c1;
  *rsc2 = c2;
}
/**************************************-**************************************/
/*                                                                           */
/*                           SPIKE_COUNT_R_SARRAYS                           */
/*                                                                           */
/*  Compute Pearson's r value for the spike counts.   Return the             */
/*  value -2.0 if there is no variance in either list of counts.             */
/*                                                                           */
/*****************************************************************************/
void spike_count_r_sarrays(s1,cnt1,s2,cnt2,n,rr,prob,start,period,sflag)
     int **s1,*cnt1,**s2,*cnt2,n;
     float *rr,*prob;
     int start,period,sflag;
{
  int i;
  float *c1,*c2,z,mean,sd1,sd2;

  c1 = (float *)myalloc(n*sizeof(float));
  c2 = (float *)myalloc(n*sizeof(float));

  for (i=0;i<n;i++){
    c1[i] = (float)count_spikes(s1[i],cnt1[i],start,period);
    c2[i] = (float)count_spikes(s2[i],cnt2[i],start,period);
  }
  mean_sdev_farray(c1,n,&mean,&sd1);
  mean_sdev_farray(c2,n,&mean,&sd2);

  if ((sd1 == 0.0)||(sd2 == 0.0)||(n < 3))
    *rr = *prob = -2.0;
  else
    if (sflag){
      simple_pearsn(c1-1,c2-1,(long)n,rr);
      *prob = 1.0;
    }else
      pearsn(c1-1,c2-1,(long)n,rr,prob,&z);

  myfree(c1); myfree(c2);
}
/**************************************-**************************************/
/*                                                                           */
/*                           SPIKE_TRIG_STIM_SARRAY                          */
/*                                                                           */
/*  Return all instances of     the float "data" within the region           */
/*  "t0" and "tn".                 "doffset" specifies the start time of     */
/*  the stimulus "data" relative to the spike times.  "dperiod" is the       */
/*  length of the "data" array.                                              */
/*                                                                           */
/*  "sstart" and "speriod" can be used to select the segment of spike data   */
/*  from which spikes will be considered.  "doffset" is in time units        */
/*  stimon==0.                                                               */
/*                                                                           */
/*****************************************************************************/
void spike_trig_stim_sarray(s,cnt,n,data,doffset,dperiod,sstart,speriod,
			    t0,tn,rtdata,rnspikes)
     int **s,*cnt,n;
     float **data;
     int doffset,dperiod,sstart,speriod,t0,tn;
     float ***rtdata;
     int *rnspikes;
{
  int i,j;
  int t,tspikes,nspikes;
  float **tdata,**ttdata;

  tspikes = 0;
  for(i=0;i<n;i++)
    tspikes += cnt[i];
  tdata = (float **)myalloc(tspikes*sizeof(float *));

  nspikes = 0;
  for(i=0;i<n;i++)
    for(j=0;j<cnt[i];j++){
      t = s[i][j];
      if ((t >= sstart)&&(t < (sstart+speriod))){
	t += t0 - doffset; // index for start of avg in data
	if ((t>=0) && ((t+tn)<=dperiod)){
	  tdata[nspikes] = copy_farray(data[i]+t,tn);
	  nspikes += 1;
	}
      }
    }
  
  if (nspikes > 0){
    ttdata = (float **)myalloc(nspikes*sizeof(float *));
    for(i=0;i<nspikes;i++)
      ttdata[i] = tdata[i];
  }else{
    ttdata = NULL;
  }

  myfree(tdata); // Only free first dimension
  *rtdata = ttdata;
  *rnspikes = nspikes;
}
/**************************************-**************************************/
/*                                                                           */
/*                           SPIKE_TRIG_AVG_SARRAY                           */
/*                                                                           */
/*  Return the average value of the float "data" within the region           */
/*  "avg_start" and "avg_period".  "doffset" specifies the start time of     */
/*  the stimulus "data" relative to the spike times.  "dperiod" is the       */
/*  length of the "data" array.                                              */
/*                                                                           */
/*  "sstart" and "speriod" can be used to select the segment of spike data   */
/*  from which spikes will be considered.  "doffset" is in time units        */
/*  stimon==0.                                                               */
/*                                                                           */
/*****************************************************************************/
void spike_trig_avg_sarray(s,cnt,n,data,doffset,dperiod,sstart,speriod,
			   avg_start,avg_period,rmean,rsdev,rnspikes)
     int **s,*cnt,n;
     float **data;
     int doffset,dperiod,sstart,speriod,avg_start,avg_period;
     float **rmean,**rsdev;
     int *rnspikes;
{
  int i,j;
  int t,tspikes,nspikes;
  float *mean,*sdev,**tdata;

  /*printf(" STIM MAX = %f\n",max_of_2d_farray(data,n,dperiod));*/

  tspikes = 0;
  for(i=0;i<n;i++)
    tspikes += cnt[i];
  tdata = (float **)myalloc(tspikes*sizeof(float *));

  nspikes = 0;
  for(i=0;i<n;i++)
    for(j=0;j<cnt[i];j++){
      t = s[i][j];
      if ((t >= sstart)&&(t < (sstart+speriod))){
	t += avg_start - doffset; // index for start of avg in data
	if ((t>=0) && ((t+avg_period)<=dperiod)){
	  tdata[nspikes] = data[i] + t;
	  nspikes += 1;
	}
      }
    }
  if (nspikes > 0){
    mean_sdev_2d_farray(tdata,nspikes,avg_period,&mean,&sdev);
    /*
    if (max_of_farray(mean,avg_period) > 1.0){
      printf("max = %f\n",max_of_farray(mean,avg_period));
      exit(0);
    }
    */
  }else{
    mean = get_zero_farray(avg_period);
    sdev = get_zero_farray(avg_period);
  }
    
  myfree(tdata); // Only free first dimension
  *rmean = mean; *rsdev = sdev; *rnspikes = nspikes;
}
/**************************************-**************************************/
/*                                                                           */
/*                          VARCOND_CHECK_CONDITION                          */
/*                                                                           */
/*  This is a utility routine for 'spike_trig_avg_sarray_varcond'            */
/*                                                                           */
/*****************************************************************************/
int varcond_check_condition(cond,t0,data,n,outfile)
     struct spike_util_cond_struct *cond;    /* Conditional RVS */
     int t0;          /* t0 is relative to start of 'data' */
     float *data;
     int n;
     char *outfile;
{
  int k;
  int flag,minmaxflag,tn;
  float x,mu,sd,vcrit,min,max;
  char ct1[SLEN],*ctype;

  /*
    printf("  VARCOND_CHECK_CONDITION\n");
    printf("    ctype = %s\n",cond->ctype);
    printf("    win0  = %d\n",cond->win0);
    printf("    winn  = %d\n",cond->winn);
    printf("    vcrit = %f\n",cond->vcrit);*/

  tn = cond->winn;

  if ((t0 < 0) || (t0+tn > n)){
    printf(" t0 = %d   cond->winn = %d\n",t0,tn);
    exit_error("VARCOND_CHECK_CONDITION","Condition window error");
  }

  ctype = cond->ctype;
  k = strlen(ctype);
  if ((ctype[k-3]=='m')&&(ctype[k-2]=='i')&&(ctype[k-1]=='n'))
    minmaxflag = -1;
  else if ((ctype[k-3]=='m')&&(ctype[k-2]=='a')&&(ctype[k-1]=='x'))
    minmaxflag = 1;
  else
    minmaxflag = 0;
  
  flag = 0;
  if (minmaxflag != 0){
    strcpy(ct1,ctype);
    ct1[k-4] = '\0';
    /*printf(" ct1 = %s   minmaxflag = %d\n",ct1,minmaxflag);*/
    
    /*** Compute the value to test. ***/
    if (strcmp(ct1,"abs_val_avg")==0){
      x = (float)sum_abs_farray(data+t0,tn)/(float)tn;
    }else if (strcmp(ct1,"sqr_avg")==0){
      x = (float)sum_square_farray(data+t0,tn)/(float)tn;
    }else if (strcmp(ct1,"var")==0){
      mean_sdev_farray(data+t0,tn,&mu,&sd);
      x = sd*sd;
    }else if (strcmp(ct1,"sd")==0){
      mean_sdev_farray(data+t0,tn,&mu,&sd);
      x = sd;
    }else if (strcmp(ct1,"avg")==0){
      mean_sdev_farray(data+t0,tn,&mu,&sd);
      x = mu;
    }else if (strcmp(ct1,"val")==0){
      if (minmaxflag == -1)
	x = min_of_farray(data+t0,tn);
      else
	x = max_of_farray(data+t0,tn);
    }else{
      printf("ctype = %s\n",ctype);
      exit_error("VARCOND_CHECK_CONDITION","Unknown condition type");
      x = 0.0;
    }
    
    if (((minmaxflag == -1) && (x >= cond->vcrit)) ||
	((minmaxflag ==  1) && (x <= cond->vcrit))){
      flag = 1;
      
      /*** Write out the scores for accepted spikes ***/
      /***{
	  char tstr[SLEN];
	  
	  sprintf(tstr,"%f\n",x);
	  append_string_to_file("zzz.VARCOND",tstr);
	  }***/
    }
  }else if ((strcmp(ctype,"range")==0)||(strcmp(ctype,"range_x")==0)){
    /*printf("  lo,hi =   %f  %f\n",cond->vcrit,cond->vcrit1);*/
    min = min_of_farray(data+t0,tn);
    max = max_of_farray(data+t0,tn);
    /*printf("    min,max = %f  %f\n",min,max);*/
    
    if (strcmp(ctype,"range")==0){
      if ((min >= cond->vcrit) && (max <= cond->vcrit1)){
	flag = 1;

	/* WYETH ***/
	if (outfile != NULL){
	  char ttstr[SLEN],tname[SLEN];
	  float tmean,tsd;
	  
	  mean_sdev_farray(data+t0,tn,&tmean,&tsd);
	  sprintf(ttstr,"%f %f\n",tmean,tsd);
	  sprintf(tname,"%s.dump",outfile);
	  append_string_to_file(tname,ttstr);
	}
      }
    }else{ /*** RANGE EXCLUDE 'range_x' ***/
      if ((min < cond->vcrit) && (max > cond->vcrit1)){
	flag = 1;
      }
      /*printf(" SATISFIED ---- lo,hi = %f  %f\n",cond->vcrit,cond->vcrit1);*/
    }
  }else{
    printf("ctype = %s\n",ctype);
    exit_error("VARCOND_CHECK_CONDITION","Unknown condition type");
  }
	    
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                      SPIKE_TRIG_AVG_SARRAY_VARCOND                        */
/*                                                                           */
/*  Only use those traces which have variance in the stimulus window that    */
/*  matches the criterion.                                                   */
/*                                                                           */
/*  The idea is to create an STA for low variance conditions, and one for    */
/*  high variance conditions and compare these.                              */
/*                                                                           */
/*****************************************************************************/
void spike_trig_avg_sarray_varcond(s,cnt,n,data,doffset,dperiod,sstart,speriod,
				   avgt0,avgtn,ctype,c0,cn,vcrit,
				   rmean,rsdev,rnspikes)
     int **s,*cnt,n;
     float **data;
     int doffset,dperiod,sstart,speriod,avgt0,avgtn;
     char ctype[];
     int c0,cn;
     float vcrit,**rmean,**rsdev;
     int *rnspikes;
{
  int i,j;
  int t,tspikes,nspikes,nreject,x0;
  float *mean,*sdev,**tdata,r;
  struct spike_util_cond_struct cond;    /* Conditional RVS */

  cond.ctype = strdup(ctype); /*** WYETH cond should come in as param ***/
  cond.win0 = c0;
  cond.winn = cn;
  cond.vcrit = vcrit;

  printf("  SPIKE_TRIG_AVG_SARRAY_VARCOND\n");
  printf("    %s c0=%d cn=%d vcrit=%f\n",ctype,c0,cn,vcrit);

  printf("  ************* WYETH CHECK THIS CODE since struct added \n");
  printf("  ************* WYETH CHECK THIS CODE since struct added \n");
  printf("  ************* WYETH CHECK THIS CODE since struct added \n");

  tspikes = 0;
  for(i=0;i<n;i++)
    tspikes += cnt[i];
  tdata = (float **)myalloc(tspikes*sizeof(float *));

  x0 = c0 - avgt0;    /* Start of criterion window rel. to STA window */

  nspikes = nreject = 0;
  for(i=0;i<n;i++) /* For each trial */


    for(j=0;j<cnt[i];j++){ /* For each spike */
      t = s[i][j];
      if ((t >= sstart)&&(t < (sstart+speriod))){
	t += avgt0 - doffset; /* index for start of avg in data */
	if ((t>=0) && ((t+avgtn)<=dperiod)){
	  tdata[nspikes] = data[i] + t;
	  if (varcond_check_condition(&cond,x0,tdata[nspikes],avgtn,NULL))
	    nspikes += 1;
	  else
	    nreject += 1;
	}
      }
    }
  if (nspikes > 0)
    mean_sdev_2d_farray(tdata,nspikes,avgtn,&mean,&sdev);
  else{
    mean = get_zero_farray(avgtn);
    sdev = get_zero_farray(avgtn);
  }

  r = 100.0 * (float)nspikes/(float)(nspikes+nreject);
  printf("    %.2f %% spikes accepted (%d of %d)\n",r,nspikes,nreject+nspikes);

  myfree(cond.ctype);
    
  myfree(tdata); /* Only free first dimension */
  *rmean = mean; *rsdev = sdev; *rnspikes = nspikes;
}
/**************************************-**************************************/
/*                                                                           */
/*                       SPIKE_TRIG_AVG_CONDITION_SARRAY                     */
/*                                                                           */
/*  This is like "spike_trig_avg_sarray", but conditions are set on the      */
/*  data traces to be averaved.                                              */
/*                                                                           */
/*****************************************************************************/
void spike_trig_avg_condition_sarray(s,cnt,n,data,doffset,dperiod,sstart,
				     speriod,avgt0,avgtn,min,max,
				     rmean,rsdev,rnspikes)
     int **s,*cnt,n;
     float **data;
     int doffset,dperiod,sstart,speriod,avgt0,avgtn;
     float min,max,**rmean,**rsdev;
     int *rnspikes;
{
  int i,j;
  int t,tspikes,nspikes;
  float *mean,*sdev,**tdata,*tf,tmin,tmax;

  tspikes = 0;
  for(i=0;i<n;i++)
    tspikes += cnt[i];
  tdata = (float **)myalloc(tspikes*sizeof(float *));

  nspikes = 0;
  for(i=0;i<n;i++)
    for(j=0;j<cnt[i];j++){
      t = s[i][j];
      if ((t >= sstart)&&(t < (sstart+speriod))){
	t += avgt0 - doffset; /* index for start of avg in data */
	if ((t>=0) && ((t+avgtn)<=dperiod)){
	  tf = data[i] + t;
	  get_min_max_farray(tf,avgtn,&tmin,&tmax);
	  if ((tmin >= min)&&(tmax <= max)){
	    tdata[nspikes] = tf;
	    nspikes += 1;
	  }
	}
      }
    }
  mean_sdev_2d_farray(tdata,nspikes,avgtn,&mean,&sdev);
  myfree(tdata); /* Only free first dimension */
  *rmean = mean; *rsdev = sdev; *rnspikes = nspikes;
}
/**************************************-**************************************/
/*                                                                           */
/*                      SPIKE_TRIG_DISTRIB_SARRAY                            */
/*                                                                           */
/*  Return the distribution of int "data" values within the region           */
/*  "avg_start" and "avg_period" relative to a spike.  The STA               */
/*  would simply be the mean.  The values "start" and "period"               */
/*  are specifications for "data" relative to the spike times.               */
/*                                                                           */
/*  "vstart" and "period" can be used to select the segment of data          */
/*  from which spikes will be considered.  "vstart" is in time units         */
/*  relative to stimon==0.                                                   */
/*                                                                           */
/*  *** WARNING:  Data must be int.                                          */
/*                                                                           */
/*****************************************************************************/
void spike_trig_distrib_sarray(s,cnt,n,data,doffset,dperiod,sstart,speriod,
			       avg_start,avg_period,vmin,vmax,g80_flag,
			       rdistrib)
     int **s,*cnt,n;
     int **data;
     int doffset,dperiod,sstart,speriod,avg_start,avg_period,vmin,vmax;
     int g80_flag;
     float ***rdistrib;
{
  int i,j,k;
  int t,dk,nspikes,nv;
  float **distrib;

  printf("  SPIKE_TRIG_DISTRIB_SARRAY\n");

  nv = vmax - vmin + 1;
  distrib = get_zero_2d_farray(avg_period,nv);

  nspikes = 0;
  if (g80_flag == 0)
    for(i=0;i<n;i++)
      for(j=0;j<cnt[i];j++){
	t = s[i][j];
	if ((t >= sstart)&&(t < (sstart+speriod))){
	  t = s[i][j] - doffset + avg_start; /* start index of avg in data */
	  if ((t>=0) && ((t+avg_period)<=dperiod)){
	    nspikes += 1;
	    for(k=t;k<(t+avg_period);k++){
	      dk = data[i][k] - vmin;
	      if ((dk >= 0)&&(dk < nv)){
		distrib[k-t][dk] += 1.0;
	      }else
		exit_error("SPIKE_TRIG_DISTRIB_SARRAY",
			   "Stimulus value out of bounds");
	    }
	  }
	}
      }
  else
    for(i=0;i<n;i++)
      for(j=0;j<cnt[i];j++){
	t = s[i][j];
	if ((t >= sstart)&&(t < (sstart+speriod))){
	  t = s[i][j] - doffset + avg_start; /* start index of avg in data */
	  if ((t>=0) && ((t+avg_period)<=dperiod) &&
	      (((t>=1248)&&((t+avg_period)<13600))||(t>=14847))){
	    
	    if ((t>=1248)&&((t+avg_period)<13600))
	      t -= 12;  /*** Correct for one frame ***/
	    else
	      t -= 25;  /*** Correct for two frames ***/

	    nspikes += 1;
	    for(k=t;k<(t+avg_period);k++){
	      dk = data[i][k] - vmin;
	      if ((dk >= 0)&&(dk < nv)){
		distrib[k-t][dk] += 1.0;
	      }else
		exit_error("SPIKE_TRIG_DISTRIB_SARRAY",
			   "Stimulus value out of bounds");
	    }
	  }
	}
      }


  printf("    %d spikes used as triggers\n",nspikes);
  multiply_2d_farray(distrib,avg_period,nv,1.0/(float)nspikes);

  for(i=0;i<nv;i++) /*** Fill distrib[0] with the unconditional distrib. ***/
    distrib[0][i] = 0.0;
  nspikes = 0;
  for(i=0;i<n;i++)
    for(j=0;j<dperiod;j++){
      dk = data[i][j] - vmin;
      if ((dk >= 0)&&(dk < nv)){
	distrib[0][dk] += 1.0;
	nspikes += 1;
      }else
	exit_error("SPIKE_TRIG_DISTRIB_SARRAY",
		   "Stimulus value out of bounds");
    }
  multiply_farray(distrib[0],nv,1.0/(float)nspikes);
  
  *rdistrib = distrib;
}
/**************************************-**************************************/
/*                                                                           */
/*                          SPIKE_TRIG_2D_AVG_TRIAL                          */
/*                                                                           */
/*  Return the average value of the float "data[t][dn]" within the region    */
/*  "avg_start" and "avg_period".  "doffset" specifies the start time of     */
/*  the stimulus "data" relative to the spike times.  "dperiod" is the       */
/*  length of the "data" array.                                              */
/*                                                                           */
/*  "sstart" and "speriod" can be used to select the segment of spike data   */
/*  from which spikes will be considered.  "doffset" is in time units        */
/*  stimon==0.                                                               */
/*                                                                           */
/*****************************************************************************/
void spike_trig_2d_avg_trial(st,n,data,doffset,dperiod,dn,sstart,speriod,
			     avg_start,avg_period,rmean,rnspikes)
     int *st,n;
     float **data;
     int doffset,dperiod,dn,sstart,speriod,avg_start,avg_period;
     float ***rmean;
     int *rnspikes;
{
  int j;
  int t,nspikes;
  float **mean,***tdata;

  tdata = (float ***)myalloc(n*sizeof(float **));
  nspikes = 0;
  for(j=0;j<n;j++){
    t = st[j];
    if ((t >= sstart)&&(t < (sstart+speriod))){
      t += avg_start - doffset; // index for start of avg in data
      if ((t>=0) && ((t+avg_period)<=dperiod)){
        tdata[nspikes] = data + t; // Or:  tdata[nspikes] = &(data[t]);
	nspikes += 1;
      }
    }
  }
  
  if (nspikes > 0)
    mean_3d_farray(tdata,nspikes,avg_period,dn,&mean);
  else
    mean = get_zero_2d_farray(avg_period,dn);

  myfree(tdata); // Free only first dimension
  *rmean = mean; *rnspikes = nspikes;
}
/**************************************-**************************************/
/*                                                                           */
/*                       SECOND_ORDER_WIENER_KERNEL_SARRAY                   */
/*                                                                           */
/*  **** WYETH this is not a complete implementation **********************  */
/*                                                                           */
/*  Compute the second order Wiener kernel for the spike data and the given  */
/*  stimulus array.                                                          */
/*                                                                           */
/*  "sstart" and "speriod" can be used to select the segment of spike data   */
/*  from which spikes will be considered.  "doffset" is the start time of    */
/*  the stimulus relative to time zero for the spikes.                       */
/*                                                                           */
/*****************************************************************************/
void second_order_wiener_kernel_sarray(s,cnt,n,data,doffset,dperiod,sstart,
				       speriod,k0,kn,rw2)
     int **s,*cnt,n;
     float **data;
     int doffset,dperiod,sstart,speriod,k0,kn;
     float ***rw2;
{
  int i,j,k,l;
  int t1,t2,count;
  float **w2;

  w2 = get_zero_2d_farray(kn,kn);

  for(i=0;i<kn;i++){ /* Tau_1 */
    printf("  tau_1 = %d\n",i);
    for(j=i;j<kn;j++){ /* Tau_2 */
      count = 0;
      for(k=0;k<n;k++) /* for each trial */
	for(l=0;l<cnt[k];l++){ /* for each spike */
	  t1 = s[k][l] + k0 + i - doffset;
	  t2 = s[k][l] + k0 + j - doffset;
	  if ((t1>0)&&(t1<dperiod)&&(t2>0)&&(t2<dperiod)){
	    w2[i][j] += data[k][t1]*data[k][t2];
	    count += 1;
	  }
	}
      w2[i][j] /= (float)count;
    }
  }
  *rw2 = w2;
}
/**************************************-**************************************/
/*                                                                           */
/*                     GET_PATTERN_TRIG_SARRAY_FROM_SARRAY                   */
/*                                                                           */
/*  Return the sarray containing all spike trains triggered on a particular  */
/*  stimulus pattern.                                                        */
/*                                                                           */
/*  NOTES:                                                                   */
/*   - Returned spike times are given relative to start of window = 0.       */
/*                                                                           */
/*****************************************************************************/
void get_pattern_trig_sarray_from_sarray(s,cnt,n,stim,ns,start,expf,pat,npat,
					 win0,winn,rs,rcnt,rn)
     int **s,*cnt,n;
     float **stim;
     int *ns;
     int start; /* stim start time, relative to spikes */
     float expf; /* sampling units per stim value */
     int *pat,npat; /* stimulus trigger pattern */
     int win0,winn; /* start and length of window for averaging */
     int ***rs,**rcnt,*rn;
{
  int i,j,k,l;
  int **ss,*scnt,sn,period,match,t;

  sn = 0; /*** Count the number of matches. ***/
  for(l=0;l<n;l++){ /*** For each stimulus/spike-train pair. ***/
    period = (int)((float)ns[l]*expf);
    for(i=0;i<(ns[l]-npat+1);i++){
      k = 0;
      match = 1;
      while((k<npat)&&(match))
	if ((int)stim[l][i+k] != pat[k])
	  match = 0;
	else
	  k += 1;
      /**** WYETH - 'start' was unused!!!! added 3/15/01 ***/
      /*t = win0 + (int)((float)(i+npat-1)*expf);*/
      t = win0 + start + (int)((float)(i+npat-1)*expf);
      if ((k==npat)&&(t+winn < period)&&(t>=0)){ /*** Found the pattern. ***/
	sn += 1;
      }
    }
  }
  ss = (int **)myalloc(sn*sizeof(int *));
  scnt = (int *)myalloc(sn*sizeof(int));
  j = 0;
  for(l=0;l<n;l++){ /*** For each stimulus/spike-train pair. ***/
    period = (int)((float)ns[l]*expf);
    for(i=0;i<(ns[l]-npat+1);i++){
      k = 0;
      match = 1;
      while((k<npat)&&(match))
	if ((int)stim[l][i+k] != pat[k])
	  match = 0;
	else
	  k += 1;
      /**** WYETH - 'start' was unused!!!! added 3/15/01 ***/
      /*t = win0 + (int)((float)(i+npat-1)*expf);*/
      t = win0 + start + (int)((float)(i+npat-1)*expf);
      if ((k==npat)&&(t+winn < period)&&(t>=0)){ /*** Found the pattern. ***/
	ss[j] = extract_spikes(s[l],cnt[l],t,winn,&scnt[j],1);
	j += 1;
      }
    }
  }
  *rs = ss; *rcnt = scnt; *rn = sn;
}
/**************************************-**************************************/
/*                                                                           */
/*                           PATTERN_TRIG_AVG_SARRAY                         */
/*                                                                           */
/*  Return the average response during the period relative to the end of     */
/*  the specified stimulus pattern.                                          */
/*                                                                           */
/*****************************************************************************/
int pattern_trig_avg_sarray(s,cnt,n,stim,ns,start,expf,pat,npat,win0,winn,rav)
     int **s,*cnt,n; /* n - number of spike trains (number of "trials") */
                     /* cnt[0..n-1] - number of spikes on each trial */
                     /* s[0..n-1][0..cnt[0..n-1]-1] - spike times */
     float **stim;   /* stim[0..n-1][0..ns[n]]*/
     int *ns;        /* ns[0..n-1] - number of stim values for each trial */
     int start;      /* stim start time, relative to spikes */
     float expf;     /* sampling units (relative to spikes) per stim value */
     int *pat;       /* pat[0..npat-1] - stimulus trigger pattern */
     int npat;       /* npat - length of stimulus pattern */
     int win0,winn;  /* start and length of window for returned response */
     float **rav;    /* pointer to average response to be returned */
{
  int i,j,k,l;
  int period,pcount,match,t;
  float *sx,*avg;

  /*printf("start=%d expf=%.4f winn,0= %d %d\n",start,expf,win0,winn);*/

  avg = get_zero_farray(winn);  // Allocate and fill response array with 0

  pcount = 0; // Number of times stimulus pattern 'pat' is found in 'stim'

  for(l=0;l<n;l++){ // For each stimulus/spike-train pair (each trial)
    period = (int)((float)ns[l]*expf);

    // printf("cnt[%d] = %d  period = %d\n",l,cnt[l],period);

    // Get array 'period' long of 0's, but 1's where spikes occurred
    sx = expand_spike_array(s[l],cnt[l],0,period);

    for(i=0;i<(ns[l]-npat+1);i++){  // Check every offset within stimulus

      k = 0;  // position in 'pat'
      match = 1; // Assume it matches, until we find it doesn't
      while((k<npat)&&(match)){  // Does 'pat' match 'stim' at offset 'i'?
	if ((int)stim[l][i+k] != pat[k])
	  match = 0; // 'pat' no longer matches 'stim' here
	else
	  k += 1;    // 'pat' still matches stim, keep going
      }

      // Determine time 't' in response 'sx' corresponding to this pattern
      t = win0 + start + (int)((float)(i+npat-1)*expf);

      if ((k==npat)&&(t+winn < period)&&(t>=0)){ // The pattern was found
	pcount += 1;
	for(j=0;j<winn;j++)  // Add this response into the average
	  avg[j] += sx[t+j];
      }
    }
    myfree(sx);
  }
  if (pcount > 0)
    multiply_farray(avg,winn,1.0/(float)pcount); // Divide total to get avg

  *rav = avg;     // Return the average resposne in the pointer
  return pcount;  // Return the total number of times the pattern was found
}
/**************************************-**************************************/
/*                                                                           */
/*                             PATTERN_TRIG_AVG_SELF                         */
/*                                                                           */
/*  Return the average STIMULUS during the period relative to the end of     */
/*  the specified stimulus pattern.                                          */
/*                                                                           */
/*****************************************************************************/
int pattern_trig_avg_self(n,stim,ns,expf,pat,npat,win0,winn,rav)
     int n;
     float **stim;
     int *ns;
     float expf; /* sampling units per stim value */
     int *pat,npat; /* stimulus trigger pattern */
     int win0,winn; /* start and length of window for averaging */
     float **rav;
{
  int i,j,k,l;
  int period,pcount,match,t,ixf;
  float *sx,*avg;

  /*printf("expf=%.4f winn,0= %d %d\n",expf,win0,winn);*/

  ixf = my_rint(expf);

  avg = get_zero_farray(winn);
  pcount = 0;
  for(l=0;l<n;l++){ /*** For each stimulus/spike-train pair. ***/
    period = ns[l]*ixf;
    sx = get_zero_farray(period); /*** make expanded stimulus ***/
    k = 0;
    for(i=0;i<ns[l];i++){
      for(j=0;j<ixf;j++){
	sx[k] = stim[l][i];
	k += 1;
      }
    }

    for(i=0;i<(ns[l]-npat+1);i++){
      k = 0;
      match = 1;
      while((k<npat)&&(match))
	if ((int)stim[l][i+k] != pat[k])
	  match = 0;
	else
	  k += 1;

      t = win0 + (i+npat-1)*ixf;
      
      if ((k==npat)&&(t+winn < period)&&(t>=0)){ /*** Found the pattern. ***/
	pcount += 1;
	for(j=0;j<winn;j++)
	  avg[j] += sx[t+j];
      }
    }
    myfree(sx);
  }
  if (pcount > 0)
    multiply_farray(avg,winn,1.0/(float)pcount);

  *rav = avg;
  return pcount;
}
/**************************************-**************************************/
/*                                                                           */
/*                        MULTI_PATTERN_TRIG_AVG_SARRAY                      */
/*                                                                           */
/*  Return the average response during the period relative to the end of     */
/*  the specified stimulus patterns.                                         */
/*                                                                           */
/*****************************************************************************/
void multi_pattern_trig_avg_sarray(s,cnt,n,stim,ns,start,expf,pat,pcnt,np,
				   win0,winn,evenflag,rav,rnav)
     int **s,*cnt,n;
     float **stim;
     int *ns;
     int start;          // stim start time, relative to spikes
     float expf;         // sampling units per stim value
     int **pat,*pcnt,np; // trigger patterns
     int win0,winn;      // start and length of window for averaging
     int evenflag;
     float ***rav;
     int **rnav;
{
  int i,j,k,l,p;
  int i0,iinc,period,*pcount,match,t;
  float *sx,**avg;

  //printf("start=%d expf=%.4f winn,0= %d %d\n",start,expf,win0,winn);

  // Can limit search to patterns starting on even or odd elements
  i0 = 0;
  iinc = 1;
  if (evenflag >= 0){
    iinc = 2;
    if (evenflag == 0)
      i0 = 1;
    else if (evenflag == 1)
      i0 = 0;
  }

  avg = get_zero_2d_farray(np,winn);
  pcount = get_zero_iarray(np);
  for(l=0;l<n;l++){ // For each stimulus/spike-train pair
    period = (int)((float)ns[l]*expf);
    sx = expand_spike_array(s[l],cnt[l],0,period);

    //append_farray_plot("zzz.pl","sx",sx,period,1);
    //append_farray_plot("zzz.pl","stim",stim[l],ns[l],1);
    //exit(0);

    for(p=0;p<np;p++){ // For each pattern
      for(i=i0;i<(ns[l]-pcnt[p]+1);i+=iinc){
	k = 0;
	match = 1;
	while((k < pcnt[p]) && match)
	  if ((int)stim[l][i+k] != pat[p][k]) // WYETH - ints match floats?
	    match = 0;
	  else
	    k += 1;

	t = win0 + start + (int)((float)(i+pcnt[p]-1)*expf);
	if ((k==pcnt[p])&&(t+winn < period)&&(t>=0)){ // Found pattern
	  pcount[p] += 1;
	  for(j=0;j<winn;j++)
	    avg[p][j] += sx[t+j];
	}
      }
    }
    myfree(sx);
  }
  for(i=0;i<np;i++){
    if (pcount[i] > 0)
      multiply_farray(avg[i],winn,1.0/(float)pcount[i]);
    /***
    printf("pcount[%d] = %d\n",i,pcount[i]);
    for(j=0;j<pcnt[i];j++)
      printf("%d ",pat[i][j]);
    printf("\n");***/
  }

  *rav = avg;
  *rnav = pcount;
}
/**************************************-**************************************/
/*                                                                           */
/*                        MULTI_PATTERN_TRIG_AVG_FARRAY                      */
/*                                                                           */
/*  Return the average response during the period relative to the end of     */
/*  the specified stimulus patterns.                                         */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Returns mean *and* SD.                                                 */
/*                                                                           */
/*****************************************************************************/
void multi_pattern_trig_avg_farray(fdata,n,fdn,stim,ns,start,expf,pat,pcnt,np,
				   win0,winn,evenflag,rav,rsd,rnav)
     float **fdata;
     int n,fdn;
     float **stim;
     int *ns;
     int start; /* stim start time, relative to spikes */
     float expf; /* sampling units per stim value */
     int **pat,*pcnt,np; /* trigger patterns */
     int win0,winn; /* start and length of window for averaging */
     int evenflag;
     float ***rav,***rsd;
     int **rnav;
{
  int i,k,l,p;
  int i0,iinc,period,*pcount,match,t;
  float **avg,**sd,***fp;

  /*printf("start=%d expf=%.4f winn,0= %d %d\n",start,expf,win0,winn);*/

  /*** Can limit search to patterns starting on even or odd elements. ***/
  i0 = 0;
  iinc = 1;
  if (evenflag >= 0){
    iinc = 2;
    if (evenflag == 0)
      i0 = 1;
    else if (evenflag == 1)
      i0 = 0;
  }

  fp = (float ***)myalloc(np*sizeof(float **));
  for(i=0;i<np;i++)
    fp[i] = (float **)myalloc(fdn*sizeof(float *)); /* May not be max */

  /***avg = get_zero_2d_farray(np,winn); BEFORE SD ***/
  avg = (float **)myalloc(np*sizeof(float *));
  sd = (float **)myalloc(np*sizeof(float *));
  pcount = get_zero_iarray(np);
  for(l=0;l<n;l++){ /*** For each stimulus/spike-train pair. ***/
    period = (int)((float)ns[l]*expf);
    for(p=0;p<np;p++){ /* For each pattern. */
      for(i=i0;i<(ns[l]-pcnt[p]+1);i+=iinc){
	k = 0;
	match = 1;
	while((k<pcnt[p])&&(match))
	  if ((int)stim[l][i+k] != pat[p][k])
	    match = 0;
	  else
	    k += 1;
	t = win0 + start + (int)((float)(i+pcnt[p]-1)*expf);
	if ((k==pcnt[p])&&(t+winn < period)&&(t>=0)){ /*** Found pattern. ***/
	  fp[p][pcount[p]] = &(fdata[l][t]); /*&(sx[t]);*/
	  pcount[p] += 1;
	}
      }
    }
  }

  /*** FOR SD ***/
  for(i=0;i<np;i++){
    if (pcount[i] > 0)
      mean_sdev_2d_farray(fp[i],pcount[i],winn,&(avg[i]),&(sd[i]));
  }
  for(i=0;i<np;i++)
    myfree(fp[i]);
  myfree(fp);

  *rav = avg; *rsd = sd; *rnav = pcount;
}
/**************************************-**************************************/
/*                                                                           */
/*                         AVG_SPIKES_AT_TRIG_SARRAY                         */
/*                                                                           */
/*  Effectively performs a cross-correlation between the trigger array "ts"  */
/*  and the spike array "s".  This is meant to be used when the triggers     */
/*  are small in number relative to the spikes.                              */
/*                                                                           */
/*  The result is smoothed if "sigma" > 0.0.                                 */
/*                                                                           */
/*****************************************************************************/
void avg_spikes_at_trig_sarray(s,cnt,n,ts,tcnt,toffset,start,period,avg_start,
			   avg_period,sigma,ravg,rnavg)
     int **s,*cnt,n,**ts,*tcnt;
     int toffset,start,period,avg_start,avg_period;
     float sigma,**ravg;
     int *rnavg;               // Number of triggers included in average
{
  int i,j,k;
  int ia,ib,ss,t,navg;
  float *avg,*smooth;

  avg = get_zero_farray(avg_period);
  navg = 0;

  /*printf(" WYETH - is start used correctly here??? *******\n");*/
  /*printf(" WYETH - it was left out of several routines. *******\n");*/

  for(i=0;i<n;i++){
    for(j=0;j<tcnt[i];j++){ // For each trigger
      t = ts[i][j] + toffset;
      ia = t + avg_start;
      ib = ia + avg_period - 1;
      if ((ia >= start)&&(ib < start+period)){
	navg += 1;
	for(k=0;k<cnt[i];k++){ // For each spike
	  ss = s[i][k] - ia;
	  if ((ss >= 0)&&(ss < avg_period))
	    avg[ss] += 1.0;
	}
      }
    }
  }

  /*append_farray_plot("zzz.zzz.pl","avg",avg,avg_period,1);*/

  if (navg > 0)
    multiply_farray(avg,avg_period,1.0/(float)navg);

  if (sigma > 0.0){
    smooth = smooth_with_gaussian(avg,avg_period,sigma,0.01);
    myfree(avg);
    avg = smooth;
  }

  *ravg = avg; *rnavg = navg;
}
/**************************************-**************************************/
/*                                                                           */
/*                     SPIKEU_PATT_TRIG_AVG_EVCODE_SARRAY                    */
/*                                                                           */
/*  Given an event code stimulus, where the stimulus values 'stim' and       */
/*  times 'stimt' are integers, return the PTA of the spikes for the         */
/*  integer pattern 'patt'.                                                  */
/*                                                                           */
/*  The result is smoothed if "sigma" > 0.0.                                 */
/*                                                                           */
/*****************************************************************************/
void spikeu_patt_trig_avg_evcode_sarray(s,cnt,n,stim,stimt,nstim,patt,npat,
					toffset,start,period,
					avg_start,avg_period,sigma,
					ravg,rnavg)
     int **s,*cnt,n;    // Spike times; number of spikes; trials
     int **stim;        // [n][nstim] stimulus code values
     int **stimt;       // [n][nstim] stimulus times
     int  *nstim;       // [n] Number of stimulus codes
     int  *patt;        // [npat] Pattern values
     int   npat;        // length of pattern
     int   toffset;     //
     int   start;       //
     int   period;      //
     int   avg_period;  // 
     int   avg_start;   // 
     float sigma;       // SD for Gaussian smoothing
     float **ravg;      // [avg_period] Average spike response
     int    *rnavg;     // Number of instances in average
{
  int i,j;
  int *trig,*trign,trn;
  int **trigt;

  trigt = (int **)myalloc(n*sizeof(int *));
  trig  = (int  *)myalloc(n*sizeof(int));
  trign = (int  *)myalloc(n*sizeof(int));

  for(i=0;i<n;i++){

    // Get the indices at which pattern occurred
    iarray_pattern_search(stim[i],nstim[i],patt,npat,&trig,&trn);

    trigt[i] = (int *)myalloc(trn*sizeof(int));

    for(j=0;j<trn;j++)
      trigt[i][j] = stimt[i][trig[j]];  // Save times associated with patterns

    trign[i] = trn;

    /*
    printf(" stimulus  nstim[%d] = %d\n",i,nstim[i]);
    for(j=0;j<nstim[i];j++)
      printf("    %d\n",stim[i][j]);
    */

    /*
    printf(" TRIGGERS (pat[0] %d) = %d\n",patt[0],trn);
    for(j=0;j<trn;j++)
      printf("  %d\n",trigt[i][j]);
    */

  }

  avg_spikes_at_trig_sarray(s,cnt,n,trigt,trign,toffset,start,period,avg_start,
			    avg_period,sigma,ravg,rnavg);

  free_2d_iarray(trigt,n);
  myfree(trig);
  myfree(trign);
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_COINCIDENCE_SPIKES                           */
/*                                                                           */
/*  Compute spike trains that have only those spikes which fall              */
/*  within "winsize" units of spikes in the other trial.                     */
/*                                                                           */
/*****************************************************************************/
void get_coincidence_spikes(d1,n1,d2,n2,rd1,rn1,rd2,rn2,start,period,winsize)
     int *d1,n1,*d2,n2,**rd1,*rn1,**rd2,*rn2;
     int start,period,winsize;
{
  int i,j;
  int *t1,*t2,tn1,tn2,done,dt,match;

  t1 = (int *)myalloc(period*sizeof(int));
  t2 = (int *)myalloc(period*sizeof(int));

  i = j = 0; /*** Advance pointers to first valid spikes ***/
  done = 0;
  while(!done && (i<n1))
    if (d1[i] < start)
      i += 1;
    else
      done = 1;
  done = 0;
  while(!done && (j<n2))
    if (d2[j] < start)
      j += 1;
    else
      done = 1;

  match = 0; /* Avoid -Wall warning */
  tn1 = tn2 = 0;
  if ((i<n1)&&(j<n2)){
    match = done = 0;
    dt = d1[i] - d2[j];
    if ((dt >= -winsize)&&(dt <= winsize)){
      match = 1;
      t1[tn1] = d1[i];
      tn1 += 1;
      t2[tn2] = d2[j];
      tn2 += 1;
    }
  }else
    done = 1;

  while(!done){
    while(match && !done){
      if ((i < n1-1)&&(j < n2-1)){
	if ((d2[j+1] - d1[i]) <= (d1[i+1] - d2[j])){
	  j += 1;
	  dt = d1[i] - d2[j];
	  if ((dt >= -winsize)&&(dt <= winsize)){
	    t2[tn2] = d2[j];
	    tn2 += 1;
	  }else
	    match = 0;
	}else{
	  i += 1;
	  dt = d1[i] - d2[j];
	  if ((dt >= -winsize)&&(dt <= winsize)){
	    t1[tn1] = d1[i];
	    tn1 += 1;
	  }else
	    match = 0;
	}
      }else if (i < n1-1){
	i += 1;
	dt = d1[i] - d2[j];
	if ((dt >= -winsize)&&(dt <= winsize)){
	  t1[tn1] = d1[i];
	  tn1 += 1;
	}else
	  match = 0;
      }else if (j < n2-1){
	j += 1;
	dt = d1[i] - d2[j];
	if ((dt >= -winsize)&&(dt <= winsize)){
	  t2[tn2] = d2[j];
	  tn2 += 1;
	}else
	  match = 0;
      }else
	done = 1;
    }
    while(!match && !done){
      if (d1[i] <= d2[j])
	i += 1;
      else
	j += 1;
      if ((i >= n1)||(j >= n2))
	done = 1;
      else if ((d1[i] >= start+period)||(d2[j] >= start+period))
	done = 1;
      else{
	dt = d1[i] - d2[j];
	if ((dt >= -winsize)&&(dt <= winsize)){
	  match = 1;
	  t1[tn1] = d1[i];
	  tn1 += 1;
	  t2[tn2] = d2[j];
	  tn2 += 1;
	}
      }
    }
  }
/***
  printf("T1:  ");
  for(i=0;i<tn1;i++)
    printf("%d ",t1[i]);
  printf("\n");
  printf("T2:  ");
  for(i=0;i<tn2;i++)
    printf("%d ",t2[i]);
  printf("\n");
***/

  *rd1 = copy_iarray(t1,tn1);
  *rd2 = copy_iarray(t2,tn2);
  myfree(t1); myfree(t2);
  *rn1 = tn1; *rn2 = tn2;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_COINCIDENCE_SARRAYS                         */
/*                                                                           */
/*  Compute sarrays that have only those spikes which fall within "winsize"  */
/*  units of spikes on the other sarray.                                     */
/*                                                                           */
/*****************************************************************************/
void get_coincidence_sarrays(s1,cnt1,s2,cnt2,n,cs1,ccnt1,cs2,ccnt2,
			     start,period,winsize)
     int **s1,*cnt1,**s2,*cnt2,n,***cs1,**ccnt1,***cs2,**ccnt2;
     int start,period,winsize;
{
  int i;
  int **c1,**c2,*cn1,*cn2;

  c1 = (int **)myalloc(n*sizeof(int *));
  c2 = (int **)myalloc(n*sizeof(int *));
  cn1 = (int *)myalloc(n*sizeof(int));
  cn2 = (int *)myalloc(n*sizeof(int));

  for(i=0;i<n;i++)
    get_coincidence_spikes(s1[i],cnt1[i],s2[i],cnt2[i],&c1[i],&cn1[i],
			   &c2[i],&cn2[i],start,period,winsize);
  
  *cs1 = c1; *cs2 = c2; *ccnt1 = cn1; *ccnt2 = cn2;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_SUPPRESS_SPIKES_SARRAYS                      */
/*                                                                           */
/*  Remove the spikes in 's2' occuring in the window (win0,winn) around      */
/*  spikes in 's1'.                                                          */
/*                                                                           */
/*****************************************************************************/
void get_suppress_spikes_sarrays(s1,cnt1,s2,cnt2,n,win0,winn,rs,rcnt)
     int **s1,*cnt1,**s2,*cnt2,n,win0,winn,***rs,**rcnt;
{
  int i,j,k;
  int **s,*cnt,flag,t,t1;

  s = (int **)myalloc(n*sizeof(int *));
  cnt = (int *)myalloc(n*sizeof(int));

  for(i=0;i<n;i++){
    s[i] = (int *)myalloc(cnt2[i]*sizeof(int));
    cnt[i] = 0;
    for(j=0;j<cnt2[i];j++){
      flag = 1;
      t = s2[i][j];
      for(k=0;k<cnt1[i];k++){
	t1 = s1[i][k] + win0;
	if ((t >= t1)&&(t < t1+winn))
	  flag = 0;
      }
      if (flag){
	s[i][cnt[i]] = t;
	cnt[i] += 1;
      }
    }
  }
  *rs = s; *rcnt = cnt;
}
/**************************************-**************************************/
/*                                                                           */
/*                          TRANSMISSION_STAT_SARRAYS                        */
/*                                                                           */
/*  Compute statistics related to transmission between two spike trains.     */
/*                                                                           */
/*  'tr' - prob that spike occurs in 's2' (in win0...win0+winn) relative to  */
/*         a spike in 's1'.                                                  */
/*                                                                           */
/*  'sflag'                                                                  */
/*   0 - simple calculation, no special processing for duplicate triggers    */
/*                                                                           */
/*  Examples:                                                                */
/*    simultaneous spikes - winn=0 win0=1                                    */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/
void transmission_stat_sarrays(s1,cnt1,s2,cnt2,n,start,period,win0,winn,
			       sflag,rntrig,rntrans)
     int **s1,*cnt1,**s2,*cnt2,n,start,period,win0,winn,sflag;
     int *rntrig,*rntrans;
{
  int i,j,k;
  int flag,t,t2,a,b,spm1,ntrig,ntrans;

  spm1 = start + period - 1;

  if (winn < 0)
    exit_error("TRANSMISSION_STAT_SARRAYS","winn must be greater than 0");

  ntrig = ntrans = 0;
  for(i=0;i<n;i++){
    for(j=0;j<cnt1[i];j++){
      t = s1[i][j];
      a = t + win0;
      b = a + winn - 1;
      if ((t >= start)&&(a >= start)&&(t<=spm1)&&(b<=spm1)){
	ntrig += 1;
	flag = 0;
	for(k=0;k<cnt2[i];k++){
	  t2 = s2[i][k];
	  if ((t2 >= a)&&(t2 <= b)){
	    flag += 1;
	    /*printf("a,b,t2  %d %d   %d\n",a,b,t2);*/
	  }
	}
	if (flag > 0)
	  ntrans += 1;
      }
    }
  }
  /*printf("  ntrig,ntrans = %d %d\n",ntrig,ntrans);*/
  
  *rntrig = ntrig; *rntrans = ntrans;
}
/**************************************-**************************************/
/*                                                                           */
/*                                SHIFT_SARRAY                               */
/*                                                                           */
/*  Shift the records in the sarray by "k" positions, where k < n.           */
/*  Trial 0 becomes trial k, trial 1 becomes trial k+1, etc.                 */
/*                                                                           */
/*****************************************************************************/
void shift_sarray(s,cnt,n,k)
     int **s,*cnt,n,k;
{
  int i,j;
  int **ts,*tcnt;

  if (k >= n)
    exit_error("SHIFT_SARRAY","Shift k > n");

  ts = (int **)myalloc(n*sizeof(int *));
  tcnt = (int *)myalloc(n*sizeof(int));

  for(i=0;i<n;i++){
    j = i+k; /* take from jth, put at ith */
    if (j >= n)
      j -= n;
    ts[i] = s[j];
    tcnt[i] = cnt[j];
  }
  for(i=0;i<n;i++){ 
    s[i] = ts[i]; /*** This should switch the pointers in s. ***/
    cnt[i] = tcnt[i];
  }
  myfree(ts); myfree(tcnt);
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_FILTERED_SARRAY                               */
/*                                                                           */
/*  Return an "sarray" that has only those trials that satisfy the           */
/*  specified conditions.  Apply conditions with values >= 0.0.              */
/*                                                                           */
/*    nsd   accept spike trains that are within this number of               */
/*          standard deviations from the mean spike count.                   */
/*    tmin  accept spike trains that have a spike at or before this          */
/*          time.                                                            */
/*    tmax  accept spike trains that have a spike at or after this           */
/*          time.                                                            */
/*    minspikes  accept spike trains with this many or more spikes.          */
/*                                                                           */
/*  *** NOTE:  Storage is allocated for spike arrays!  This is               */
/*  different than the usual "sarray" which is simply a list of              */
/*  pointers to storage existing in a multi_channel structure.               */
/*                                                                           */
/*****************************************************************************/
void get_filtered_sarray(s,cnt,n,rs,rcnt,rn,nsd,tmin,tmax,minspikes)
     int **s,*cnt,n,***rs,**rcnt,*rn; /* spike data, and return data */
     float nsd;
     int tmin,tmax,minspikes; 
{
  int i,k;
  int **ts,*tcnt,tn,*flag;
  int pflag = 0;
  float mean,sdev;

  flag = (int *)myalloc(n*sizeof(int)); /*** Set all flags ***/
  for(i=0;i<n;i++)
    flag[i] = 1;

  if (nsd >= 0.0){ /*** Flag by number of standard deviations from mean ***/
    mean_sdev_iarray(cnt,n,&mean,&sdev);
    for(i=0;i<n;i++)
      if (nsd < fabs(((float)cnt[i] - mean)/sdev)){
	flag[i] = 0;
	if (pflag) printf("  (NSD) i = %d\n",i);
      }
  }
  if (tmin >= 0){
    for(i=0;i<n;i++)
      if (cnt[i] > 0){
	if (s[i][0] > tmin){
	  flag[i] = 0;
	  if (pflag) printf("  (TMIN) i = %d\n",i);
	}
      }else
	flag[i] = 0;
  }
  if (tmax >= 0){
    for(i=0;i<n;i++)
      if (cnt[i] > 0){
	if (s[i][cnt[i]-1] < tmax){
	  flag[i] = 0;
	  if (pflag) printf("  (TMAX) i = %d\n",i);
	}
      }else
	flag[i] = 0;
  }
  if (minspikes >= 0){
    for(i=0;i<n;i++)
      if (cnt[i] < minspikes){
	flag[i] = 0;
	if (pflag) printf("  (MINSPIKES) i = %d\n",i);
      }
  }
  
  tn = sum_iarray(flag,n,0,n);
  ts = (int **)myalloc(tn*sizeof(int *));
  tcnt = (int *)myalloc(tn*sizeof(int));
  k = 0;
  for(i=0;i<n;i++)
    if (flag[i]){
      ts[k] = copy_iarray(s[i],cnt[i]);
      tcnt[k] = cnt[i];
      k += 1;
    }
  myfree(flag);

  *rs = ts; *rcnt = tcnt; *rn = tn;
}
/**************************************-**************************************/
/*                                                                           */
/*                     GET_SC_FILTERED_DOUBLE_SARRAY                         */
/*                                                                           */
/*  Return corresponding sarrays with only those trials which have           */
/*  the minimum spike count for both arrays.                                 */
/*                                                                           */
/*  *** NOTE:  Storage is allocated for spike arrays!  This is               */
/*  different than the usual "sarray" which is simply a list of              */
/*  pointers to storage existing in a multi_channel structure.               */
/*                                                                           */
/*****************************************************************************/
void get_sc_filtered_double_sarray(s1,cnt1,s2,cnt2,n,rs1,rcnt1,rs2,rcnt2,rn,
				   minspikes)
     int **s1,*cnt1,**s2,*cnt2,n,***rs1,**rcnt1,***rs2,**rcnt2,*rn;
     int minspikes;
{
  int i,k;
  int **ts1,*tcnt1,**ts2,*tcnt2,tn,*flag;

  flag = (int *)myalloc(n*sizeof(int)); /*** Set all flags ***/
  for(i=0;i<n;i++)
    flag[i] = 1;

  for(i=0;i<n;i++)
    if ((cnt1[i] < minspikes)||(cnt2[i] < minspikes))
      flag[i] = 0;
  
  tn = sum_iarray(flag,n,0,n);
  ts1 = (int **)myalloc(tn*sizeof(int *));
  ts2 = (int **)myalloc(tn*sizeof(int *));
  tcnt1 = (int *)myalloc(tn*sizeof(int));
  tcnt2 = (int *)myalloc(tn*sizeof(int));
  k = 0;
  for(i=0;i<n;i++)
    if (flag[i]){
      ts1[k] = copy_iarray(s1[i],cnt1[i]);
      ts2[k] = copy_iarray(s2[i],cnt2[i]);
      tcnt1[k] = cnt1[i];
      tcnt2[k] = cnt2[i];
      k += 1;
    }
  myfree(flag);

  *rs1 = ts1; *rcnt1 = tcnt1; *rs2 = ts2; *rcnt2 = tcnt2; *rn = tn;
}
/**************************************-**************************************/
/*                                                                           */
/*                         CHECK_TRIAL_CONDITION_SU                          */
/*                                                                           */
/*  This routine renamed (old name: check_trial_condition)!                  */
/*  ******** USE "check_trial_condition" in "multi_channel.c"                */
/*  11/10/95.                                                                */
/*                                                                           */
/*****************************************************************************/
int check_trial_condition_su(coherency,direction,series_type,value)
     int coherency,direction;
     char series_type[];
     int value;
{
  int cflag,dflag,allflag;
  int match;
  
  cflag = dflag = allflag = 0;
  if (strcmp(series_type,"cseries")==0)
    cflag = 1;
  else if (strcmp(series_type,"dseries")==0)
    dflag = 1;
  else if (strcmp(series_type,"single")==0)
    allflag = 1;

  if ((cflag && (value == coherency))||
      (dflag && (value == direction))|| allflag)
    match = 1;
  else
    match = 0;
  return match;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_SARRAY_FROM_EVENTS                            */
/*                                                                           */
/*  Convert an array of "struct event_trial" to an "sarray".                 */
/*                                                                           */
/*****************************************************************************/
void get_sarray_from_events(etrial,n,spikes,count)
     struct event_trial *etrial;
     int n; /* number of event trials */
     int ***spikes,**count;
{
  int i,j;
  int **s,*cnt;
  
  cnt = (int *)myalloc(n*sizeof(int));
  s = (int **)myalloc(n*sizeof(int *));
  for (i=0;i<n;i++){
    cnt[i] = etrial[i].nev;
    s[i] = (int *)myalloc(etrial[i].nev*sizeof(int));
    for(j=0;j<cnt[i];j++)
      s[i][j] = (etrial[i].ev[j].t[0] +
		 etrial[i].ev[j].t[etrial[i].ev[j].n-1])/2;
  }
  *spikes = s;
  *count = cnt;
}
/**************************************-**************************************/
/*                                                                           */
/*                            SU_VARDUR_RATE_SARRAY                          */
/*                                                                           */
/*****************************************************************************/
float *su_vardur_rate_sarray(s,cnt,n,start,vdur,sampling)
     int **s,*cnt,n;
     int start;
     int *vdur;         // [n] duration of each trial
     float sampling;
{
  int i;
  float *data;

  data = (float *)myalloc(n*sizeof(float));

  for(i=0;i<n;i++){
    data[i] = (float)count_spikes(s[i],cnt[i],start,vdur[i]) /
              (float)(vdur[i]) * sampling;
  }

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MEAN_SPIKE_COUNT_SARRAY                         */
/*                                                                           */
/*****************************************************************************/
float mean_spike_count_sarray(s,cnt,n,start,period,sampling,rsdev)
     int **s,*cnt,n;
     int start,period;
     float sampling,*rsdev;
{
  int i;
  float avg_sc,*data;

  data = (float *)myalloc(n*sizeof(float));

  avg_sc = 0.0;
  for (i=0;i<n;i++)
    data[i] = (float)count_spikes(s[i],cnt[i],start,period);
  mean_sdev_farray(data,n,&avg_sc,rsdev);
  myfree(data);

  return avg_sc;
}
/**************************************-**************************************/
/*                                                                           */
/*                    MEAN_SPIKE_COUNT_SARRAY_WINDOW                         */
/*                                                                           */
/*****************************************************************************/
float mean_spike_count_sarray_window(s,cnt,n,start,period,winwidth,rn,rsdev)
     int **s,*cnt,n,start,period,winwidth,*rn;
     float *rsdev;
{
  int i,j,k;
  int nwin,tn;
  float *data,mean;

  tn = period/winwidth;
  nwin = n * tn;
  data = (float *)myalloc(nwin*sizeof(float));

  k = 0;
  for (i=0;i<n;i++)
    for (j=0;j<tn;j++){
      data[k] = (float)count_spikes(s[i],cnt[i],start+j*winwidth,winwidth);
      k += 1;
    }
  mean_sdev_farray(data,nwin,&mean,rsdev);
  myfree(data);
  *rn = nwin;

  return mean;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MEAN_SPIKE_RATE_SARRAY                          */
/*                                                                           */
/*****************************************************************************/
float mean_spike_rate_sarray(s,cnt,n,start,period,sampling,sdev)
     int **s,*cnt,n; /* spike data, counts */
     int start,period;
     float sampling,*sdev;
{
  int i;
  float avg_sr;
  float *data;

  data = (float *)myalloc(n*sizeof(float));
  avg_sr = 0.0;
  for (i=0;i<n;i++){
    data[i] = ((float)count_spikes(s[i],cnt[i],start,period) /
	       (float)period * sampling);
  }
  mean_sdev_farray(data,n,&avg_sr,sdev);
  myfree(data);

  return avg_sr;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_NORM_SPIKE_COUNT_SARRAY                      */
/*                                                                           */
/*  Return the farray containing z-scored spike counts for the sarray.       */
/*                                                                           */
/*****************************************************************************/
float *get_norm_spike_count_sarray(s,cnt,n,start,period,sampling)
     int **s,*cnt,n;
     int start,period;
     float sampling;
{
  int i;
  float *data;

  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = (float)count_spikes(s[i],cnt[i],start,period);
  z_score_farray(data,n);
  
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_CYCLE_TRIG_AVG_TRIAL                         */
/*                                                                           */
/*  Return the cycle triggered average firing rate in spikes/sec.            */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Because we are rounding the length of the cycle triggerred average     */
/*    *down*, there may be spikes that are ignored.                          */
/*                                                                           */
/*****************************************************************************/
void get_cycle_trig_avg_trial(data,n,start,ncyc,period,sampling,ravg,rn)
     int *data,n,start,ncyc;
     float period,sampling,**ravg;
     int *rn;
{
  int i,k;
  int m,navg,nignore,nuse;
  float t,duration,*avg;

  navg = (int)period;
  avg = get_zero_farray(navg);
  duration = (float)ncyc * period;
  nignore = nuse = 0;
  for(i=0;i<n;i++){
    t = (float)(data[i] - start);
    if ((t>=0.0)&&(t<duration)){
      m = (int)(t/period); /* Added parens. */
      k = (int)(t-(float)m*period);
      if ((k < 0)||(k >= navg))
	nignore += 1;
      else{
	avg[k] += 1.0;
	nuse += 1;
      }
    }
  }
  if (nignore > 0)
    printf("  GET_CYCLE_TRIG_AVG_TRIAL: ignored %d spikes, used %d.\n",nignore,
	   nuse);
  multiply_farray(avg,navg,sampling/(float)ncyc);
  *ravg = avg; *rn = navg;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_FOURIER_HARMONIC_SARRAY                       */
/*                                                                           */
/*  Compute the mean and std. dev. for the Fourier harmonic for the          */
/*  fundamental "period" specified.  Amplitude is in spikes/sec, phase is    */
/*  in degrees.                                                              */
/*                                                                           */
/*  start - start time of epoch to analyze, typically 0.                     */
/*  fdur - (float) duration of epoch to analyze.                             */
/*  period - length of cycle (fundamental period, typically time betw. syncs)*/
/*  sampling - number samples per second, typically 1000.0                   */
/*  order - index of fourier component to compute, 1=f1, 2=f2, etc, no 0.    */
/*                                                                           */
/*****************************************************************************/
void get_fourier_harmonic_sarray(s,cnt,n,start,fdur,period,sampling,order,
                                 ramplmean,ramplsd,rphmean,rphsd,rampl)
     int **s,*cnt,n,start;
     float fdur,period,sampling;
     int order;
     float *ramplmean,*ramplsd,*rphmean,*rphsd;
     float **rampl;
{
  int i,j;
  float *ampl,*phase,t,todd,teven,fc,r,theta,*d;

  ampl = (float *)myalloc(n*sizeof(float));
  phase = (float *)myalloc(n*sizeof(float));

  fc = 2.0*M_PI / (period/(float)order);

  //printf("n=%d  cnt[0]=%d  fdur=%f period=%f sampling=%f\n",n,cnt[0],fdur,
  //period,sampling);

  for(i=0;i<n;i++){ // For each trial of spikes
    todd = teven = 0.0;
    for(j=0;j<cnt[i];j++){ // For each spike on this trial
      t = (float)(s[i][j]-start);
      if ((t > 0.0) && (t < fdur)){
        t *= fc;
        todd += sin(t);
        teven += cos(t);
      }
    }
    ampl[i] = 2.0*sqrt(todd*todd + teven*teven) / (fdur/sampling);
    phase[i] = (float)atan2(todd,teven) * 180.0/M_PI;
  }
  mean_sdev_farray(ampl,n,ramplmean,ramplsd);

  // Use vector average to find rough direction of mean
  vector_average_farray(ampl,phase,n,2,2,&r,&theta);
  // Compute mean around the direction of the vector average
  d = get_farray(n);
  for(i=0;i<n;i++){
    d[i] = get_signed_circular_diff(phase[i]+180.0,theta+180.0,360.0);
    /*printf("d[i] = %f  (%f %f)\n",d[i],phase[i]+180.0,theta+180.0);*/
  }
  mean_sdev_farray(d,n,rphmean,rphsd);
  myfree(d);

  /*** WYETH REMOVE printf("MEAN= %f\n",*rphmean); ***/
  *rphmean += theta;

  /***printf("--------  r = %f  theta = %f  ====>  %f\n",r,theta,*rphmean);***/

  *rampl = ampl; // Return individual ampl data points for each trial
  myfree(phase);
}
/**************************************-**************************************/
/*                                                                           */
/*                        GET_FOURIER_HARMONIC_SARRAY_OLD                    */
/*                                                                           */
/*  Compute the mean and std. dev. for the Fourier harmonic for the          */
/*  fundamental "period" specified.  Amplitude is in spikes/sec, phase is    */
/*  in degrees.                                                              */
/*                                                                           */
/*  *** WYETH - replaced by new routine above.  Eventually, remove this.     */
/*                                                                           */
/*****************************************************************************/
void get_fourier_harmonic_sarray_old(s,cnt,n,start,ncyc,period,sampling,order,
				     ramplmean,ramplsd,rphmean,rphsd)
     int **s,*cnt,n,start,ncyc;
     float period,sampling;
     int order;
     float *ramplmean,*ramplsd,*rphmean,*rphsd;
{
  int i;
  int nav;
  float *ampl,*phase,*avg;

  ampl = (float *)myalloc(n*sizeof(float));
  phase = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++){
    get_cycle_trig_avg_trial(s[i],cnt[i],start,ncyc,period,sampling,&avg,&nav);
    get_fourier_harmonic_farray(avg,nav,order,period,&ampl[i],&phase[i]);
    myfree(avg);
  }
  mean_sdev_farray(ampl,n,ramplmean,ramplsd);
  mean_sdev_farray(phase,n,rphmean,rphsd);
  myfree(ampl); myfree(phase);
}
/**************************************-**************************************/
/*                                                                           */
/*                     INCREMENT_PSTH_FLOAT_BIN_TRIAL                        */
/*                                                                           */
/*  Increment the PSTH using the given trial.                                */
/*                                                                           */
/*****************************************************************************/
void increment_psth_float_bin_trial(data,n,start,period,binsize,psth,nbin)
     int *data,n;
     float start,period,binsize;
     float *psth;
     int nbin;
{
  int i,k;
  float t;

  for(i=0;i<n;i++){
    t = (float)data[i] - start;
    k = (int)(t/binsize);
    if ((t>=0.0)&&(k >= 0)&&(k<nbin))
      psth[k] += 1.0;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             PSTH_FLOAT_BIN_TRIAL                          */
/*                                                                           */
/*  Compute the number of spikes per bin (with floating point bin            */
/*  size) for one trial.                                                     */
/*                                                                           */
/*  This routine will be used for the analysis of the BINARY noise           */
/*  data.  Note the stimulus period is not a multiple of millisec.           */
/*                                                                           */
/*****************************************************************************/
void psth_float_bin_trial(data,n,start,period,binsize,rpsth,rn)
     int *data,n;
     float start,period,binsize;
     float **rpsth;
     int *rn;
{
  int i,k;
  int nbin;
  float t,*psth;

  nbin = (int)(period/binsize);
  psth = get_zero_farray(nbin);

  for(i=0;i<n;i++){
    t = (float)data[i] - start;
    k = (int)(t/binsize);
    if ((t>=0.0)&&(k >= 0)&&(k<nbin))
      psth[k] += 1.0;
  }
  *rpsth = psth;
  *rn = nbin;
}
/**************************************-**************************************/
/*                                                                           */
/*                            PSTH_FLOAT_BIN_SARRAY                          */
/*                                                                           */
/*  Compute the number of spikes per bin (with floating point bin size) for  */
/*  all spike trains.                                                        */
/*                                                                           */
/*****************************************************************************/
void psth_float_bin_sarray(s,cnt,n,start,period,binsize,rpsth,rn)
     int **s,*cnt,n;
     float start,period,binsize;
     float **rpsth;
     int *rn;
{
  int i,j,k;
  int nbin;
  float t,*psth;

  nbin = (int)(period/binsize);
  psth = get_zero_farray(nbin);

  for(i=0;i<n;i++) /* For each spike train. */
    for(j=0;j<cnt[i];j++){ /* For each spike time. */
      t = (float)s[i][j] - start;
      k = (int)(t/binsize);
      if ((t>=0.0)&&(k >= 0)&&(k<nbin))
	psth[k] += 1.0;
    }
  *rpsth = psth;
  *rn = nbin;
}
/**************************************-**************************************/
/*                                                                           */
/*                         SU_PSTH_FLOAT_BIN_SARRAY_VDUR                     */
/*                                                                           */
/*  Like "PSTH_FLOAT_BIN_SARRAY" but for variable duration trials.           */
/*                                                                           */
/*****************************************************************************/
void su_psth_float_bin_sarray_vdur(s,cnt,n,start,vdur,binsize,rpsth,rn,rpsthn)
     int **s,*cnt,n;
     float start;
     int *vdur;           // trial duration [n]
     float binsize;
     float **rpsth;       // Number of counts per bin [rn]
     int *rn;             // Number of bins
     int **rpsthn;        // Number of trials covering each bin [rn]
{
  int i,j,k,ii;
  int nbin,period,*psthn,max_bin_index;
  float t,*psth;

  period = max_of_iarray(vdur,n);

  nbin = (int)(period/binsize);
  psth  = get_zero_farray(nbin);
  psthn = get_zero_iarray(nbin);  // Count number of trials contributing

  for(i=0;i<n;i++){ // For each spike train

    // Do not allow contributions to bins not fully within trial duration
    max_bin_index = (int)(vdur[i]/binsize);
    
    for(j=0;j<cnt[i];j++){ // For each spike time
      t = (float)s[i][j] - start;
      k = (int)(t/binsize);
      if ((t >= 0.0) && (k >= 0) && (k < max_bin_index)){
	psth[k] += 1.0;
      }
    }
    for(ii=0;ii<max_bin_index;ii++)  // Number of trials that contributed
      psthn[ii] += 1;
  }

  *rpsth = psth;
  *rn = nbin;
  *rpsthn = psthn;
}
/**************************************-**************************************/
/*                                                                           */
/*                             SIMPLE_PST_SARRAY                             */
/*                                                                           */
/*****************************************************************************/
int *simple_pst_sarray(s,cnt,n,start,period)
     int **s,*cnt,n,start,period;
{
  int i,j;
  int *pst,t;

  pst = get_zero_iarray(period);
  for (i=0;i<n;i++)
    for (j=0;j<cnt[i];j++){
      t = s[i][j] - start;
      if ((t>=0)&&(t<period))
	pst[t] += 1;
    }
  return pst;
}
/**************************************-**************************************/
/*                                                                           */
/*                               PSTH_PROB_SARRAY                            */
/*                                                                           */
/*****************************************************************************/
float *psth_prob_sarray(s,cnt,n,start,period)
     int **s,*cnt,n,start,period;
{
  int *ipsth;
  float *psth;

  ipsth = simple_pst_sarray(s,cnt,n,start,period); // Get PSTH
  psth = i2farray(ipsth,period); // Convert to floats
  myfree(ipsth);
  multiply_farray(psth,period,1.0/(float)n); // Represent as probability
  return psth;
}
/**************************************-**************************************/
/*                                                                           */
/*                            SIMPLE_PSTH_BIN_SARRAY                         */
/*                                                                           */
/*****************************************************************************/
int *simple_psth_bin_sarray(s,cnt,n,start,period,binsize,npsth)
     int **s,*cnt,n,start,period,binsize,*npsth;
{
  int i,j;
  int *psth,t,b;

  *npsth = (period+binsize-1)/binsize;
  psth = get_zero_iarray(*npsth);

  for (i=0;i<n;i++)
    for (j=0;j<cnt[i];j++){
      t = s[i][j] - start;
      if ((t>=0)&&(t<period)){
	b = t / binsize;
	if (b > *npsth) exit_error("SIMPLE_PSTH_BIN_SARRAY","b too big");
	psth[b] += 1.0;
      }
    }
  return psth;
}
/**************************************-**************************************/
/*                                                                           */
/*                            ADAPTIVE_PST_SARRAY                            */
/*                                                                           */
/*  Compute the PST for the data using the adaptive window method            */
/*  described in Eyal Seidemann's master's thesis.                           */
/*                                                                           */
/*****************************************************************************/
float *adaptive_pst_sarray(spikes,count,n,start,period,sampling,ncrit)
     int **spikes,*count,n; /* spike data, counts */
     int start,period;      /* in sampling units */
     float sampling;        /* samples per second */
     int ncrit;             /* minimum number of spikes to make estimate */
{
  int i,j;
  int t,total,*sum;
  int done,a,b;
  float *pst;

  pst = (float *)myalloc(period*sizeof(float));
  sum = get_zero_iarray(period);

  for (i=0;i<n;i++)
    for (j=0;j<count[i];j++){
      t = spikes[i][j] - start;
      if ((t>=0)&&(t<period))
	sum[t] += 1;
    }
  for (i=0;i<period;i++){
    a = b = i;
    done = 0;
    total = sum[a];
    while (!done)
      if (total >= ncrit){
	done = 1;
	pst[i] = (float)total/((float)((b-a+1)*n)/sampling);
      }else{
	a -= 1;
	b += 1;
	if ((a<0) || (b>=period)){
	  done = 1;
	  pst[i] = -1.0;
	}else
	  total += (sum[a] + sum[b]);
      }
  }
  myfree(sum);
  return pst;
}
/**************************************-**************************************/
/*                                                                           */
/*                       SPIKEU_SARRAY_CONDITION_SARRAY                      */
/*                                                                           */
/*  Use the spike times of 's1' to conditionally accept spikes of 's2'.      */
/*  Return 'rs' and 'rcnt' containing the accepted sarray.                   */
/*                                                                           */
/*  NOTES                                                                    */
/*  1. Used for conditional cross-correlation in nda_corr.c                  */
/*                                                                           */
/*****************************************************************************/
void spikeu_sarray_condition_sarray(s1,cnt1,s2,cnt2,n,tn,toff,trad,rs,rcnt)
     int **s1;       // [cnt1[]][n]
     int *cnt1;      // [n]
     int **s2;       // [cnt2[]][n]
     int *cnt2;      // [n]
     int n;          // Number of trials
     int tn;         // Trial duration
     int toff;       // Temporal offset
     int trad;       // Temporal radius
     int ***rs;      // [rcnt[]][n]
     int **rcnt;     // [n]
{
  int i,j,k;
  int k0,k1,t;
  int **x,**s3,*cnt3,*ts;

  printf("  SPIKEU_SARRAY_CONDITION_SARRAY\n");
  printf("    toff = %d\n",toff);
  printf("    trad = %d\n",trad);

  if (trad < 0)
    exit_error("SPIKEU_SARRAY_CONDITION_SARRAY","trad < 0");

  //
  //  (1) Expand the conditioning spikes, s1
  //
  x = get_zero_2d_iarray(n,tn);

  for(i=0;i<n;i++){  // For each trial
    for (j=0;j<cnt1[i];j++){  // For each spike
      t = s1[i][j] + toff;
      
      k0 = t - trad;
      if (k0 < 0)
	k0 = 0;
      k1 = t + trad;
      if (k1 >= tn)
	k1 = tn-1;
      
      for(k=k0;k<=k1;k++){
	x[i][k] = 1;
      }
    }
  }

  s3 = (int **)myalloc(n*sizeof(int *));
  cnt3 = (int *)myalloc(n*sizeof(int));
  ts = (int *)myalloc(tn*sizeof(int));

  for(i=0;i<n;i++){  // For each trial

    k = 0;
    for (j=0;j<cnt2[i];j++){  // For each spike
      t = s2[i][j];

      if ((t >= 0) && (t < tn)){
	if (x[i][t] == 1){
	  ts[k] = t;
	  k += 1;
	}
      }
    }
    s3[i] = copy_iarray(ts,k);
    cnt3[i] = k;
  }

  free_2d_iarray(x,n);
  myfree(ts);

  *rs = s3;
  *rcnt = cnt3;
}
/**************************************-**************************************/
/*                                                                           */
/*                        GET_FIRST_SPIKE_DISTRIB_SARRAY                     */
/*                                                                           */
/*  Return the farray containing the number of first spikes at each sample   */
/*  time.  Also, compute the precision (SD) and reliabilty (fraction of      */
/*  trials with at least one spike in the window.                            */
/*                                                                           */
/*****************************************************************************/
float *get_first_spike_distrib_sarray(s,cnt,n,start,period,rmu,rsd,rrel)
     int **s,*cnt,n;
     int start,period;
     float *rmu,*rsd,*rrel;
{
  int i,k;
  int count,t;
  float *psth,*tdata;

  psth = get_zero_farray(period);
  tdata = (float *)myalloc(n*sizeof(float));
  count = 0;
  for(i=0;i<n;i++){ /* For each spike train. */
    k = get_first_spike_in_window(s[i],cnt[i],start,period);
    t = k-start;
    if (t >= 0){  /*** WYETH - bug fixed here, used to be: t > 0 ***/
      psth[t] += 1.0;
      tdata[count] = (float)k;
      count += 1;
    }
  }
  mean_sdev_farray(tdata,count,rmu,rsd);
  myfree(tdata);
  *rrel = (float)count/(float)n;
  return psth;
}
/**************************************-**************************************/
/*                                                                           */
/*                           SIMPLE_ISI_SARRAY                               */
/*                                                                           */
/*  Return the Inter-Spike Interval Histogram for the spikes trains.         */
/*  The 0th bin contains intervals of length 0.                              */
/*                                                                           */
/*****************************************************************************/
int *simple_isi_sarray(spikes,count,n,start,period,max_isi,nn,ntot)
     int **spikes,*count,n; /* spike data, counts */
     int start,period;      /* in sampling units */
     int max_isi;           /* given in sampling units */
     int *nn,*ntot; /* return length of ISI and total intervals */
{
  int i,j;
  int interval;
  int *isi,nisi;

  nisi = max_isi+1;
  isi = get_zero_iarray(nisi);

  *ntot = 0;
  for (i=0;i<n;i++)
    for (j=0;j<(count[i]-1);j++){
      if ((spikes[i][j] >= start)&&(spikes[i][j+1] < (start+period))){
	interval = spikes[i][j+1] - spikes[i][j];
	if (interval < nisi){
	  isi[interval]  += 1;
	  *ntot += 1;
	}
      }
    }
  *nn = nisi;
  return isi;
}
/**************************************-**************************************/
/*                                                                           */
/*                        GET_LONGEST_ISI_SARRAY                             */
/*                                                                           */
/*  Return a list of the start times and lengths of the longest "k"          */
/*  ISIs.                                                                    */
/*                                                                           */
/*****************************************************************************/
void get_longest_isi_sarray(s,cnt,n,start,period,k,tisi,pisi)
     int **s,*cnt,n,start,period,k;
     int **tisi,**pisi; /* Start times and periods of longest ISIs */
{
  int i,j;
  int *isi,*t,*p,nisi,ntot,count,min_count,min_isi,tt,ti;
  
  isi = simple_isi_sarray(s,cnt,n,start,period,period,&nisi,&ntot);

  min_isi = nisi-1; /* Find the smallest length among the largest "k" ISIs */
  count = isi[min_isi];
  while((min_isi>0)&&(count<k)){
    min_isi -= 1;
    count += isi[min_isi];
  }
  if (count<k)
    exit_error("GET_LONGEST_ISI_SARRAY","Not enough ISIs");

  min_count = k-(count-isi[min_isi]);

  t = (int *)myalloc(k*sizeof(int));
  p = (int *)myalloc(k*sizeof(int));

  ti = 0;
  for (i=0;i<n;i++){ /* for each trial */
    for (j=0;j<(cnt[i]-1);j++){ /* for each ISI */
      if ((s[i][j] >= start) && (s[i][j+1] <= start+period)){
	tt = s[i][j+1] - s[i][j];
	if (tt < 0)
	  exit_error("GET_LONGEST_ISI_SARRAY","Spikes out of order");
	if (tt > min_isi){
	  t[ti] = s[i][j];
	  p[ti] = tt;
	  ti += 1;
	}else if ((tt == min_isi) && (min_count > 0)){
	  t[ti] = s[i][j];
	  p[ti] = tt;
	  ti += 1;
	  min_count -= 1;
	}
      }
    }
  }
  myfree(isi);
  *tisi = t; *pisi = p;
} 
/**************************************-**************************************/
/*                                                                           */
/*                        GET_FRACTION_ISI_SARRAY                            */
/*                                                                           */
/*  Return the length "k" of the ISI such that the integral from 0           */
/*  to "k" of the ISI distribution accounts for a fraction "f" of the        */
/*  distribution.  "k" is the minimum such value.                            */
/*                                                                           */
/*****************************************************************************/
int get_fraction_isi_sarray(s,cnt,n,start,period,f)
     int **s,*cnt,n,start,period;
     float f;
{
  int k;
  int *isi,nisi,ntot,total;
  
  isi = simple_isi_sarray(s,cnt,n,start,period,period,&nisi,&ntot);

  k = 0;
  total = isi[k];
  while((float)total/(float)ntot < f){
    k += 1;
    total += isi[k];
  }
  myfree(isi);
  return k;
} 
/**************************************-**************************************/
/*                                                                           */
/*                     GET_TIME_FRACTION_ISI_SARRAY                          */
/*                                                                           */
/*  Return the length "k" of the ISI such that the integral from 0           */
/*  to "k" of the ISI *time* distribution accounts for a fraction "f"        */
/*  of the distribution.  The time distribution is the fraction of           */
/*  time in the spike train accounted for by each ISI value.                 */
/*                                                                           */
/*****************************************************************************/
int get_time_fraction_isi_sarray(s,cnt,n,start,period,f)
     int **s,*cnt,n,start,period;
     float f;
{
  int i,k;
  int *isi,nisi,ntot;
  float *fisi,total;
  
  isi = simple_isi_sarray(s,cnt,n,start,period,period,&nisi,&ntot);
  fisi = i2farray(isi,nisi);
  for(i=1;i<nisi;i++)
    fisi[i] *= (float)i;
  norm_area_farray(fisi,nisi,1.0);

  k = 0;
  total = fisi[k];
  while(total < f){
    k += 1;
    total += fisi[k];
  }
  myfree(isi); myfree(fisi);
  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                              POWER_BURST_METRIC                           */
/*                                                                           */
/*  This routine should compute the spectrum shape measure 'P' for a set of  */
/*  power spectra 'pwr' as described in Bair et al. 1994.                    */
/*                                                                           */
/*****************************************************************************/
void power_burst_metric(pwr,m,n,dx,rp,rcat)
     float **pwr;
     int m,n;
     float dx,*rp;
     char **rcat;
{
  int i,j;
  float *mean,*t,fb,fnb,*rb,*rnb,p;
  int bstart,bstop,width; /* PASS IN THESE PARAMS */
  int imax,imin,imin0;
  int nbstart;
  
  bstart = 5; /* ORIGINAL: 3;*/ /* burst peak parameters */
  nbstart = 3;
  bstop = 13;
  width = 7;

  mean_2d_farray(pwr,m,n,&mean);
  t = get_zero_farray(n); /*** Compute running integral in window 'width' ***/
  for(i=0;i<(n-width+1);i++)
    for(j=0;j<width;j++)
      t[i] += mean[i+j];
  
  imax = bstart + max_coord_farray(t+bstart,bstop-bstart+1);
  imin = (imax-1) + min_coord_farray(t+(imax-1),(n-width)-(imax-1)+1);
  imin0 = nbstart + min_coord_farray(t+nbstart,bstop-nbstart+1);
  printf("  imin0= %d  imin= %d  imax= %d\n",imin0,imin,imax);

  /* Compute ratios for each stimulus condition. */
  rb = (float *)myalloc(m*sizeof(float));
  rnb = (float *)myalloc(m*sizeof(float));
  fb = 1.0;

  printf("********** WYETH - what shoud fnb start at? (Using 1)\n");
  fnb = 1.0;

  for(i=0;i<m;i++){
    rb[i] = sum_farray(pwr[i],n,imax,width) / sum_farray(pwr[i],n,imin,width);
    rnb[i] = sum_farray(pwr[i],n,imin0,width) / (float)width;
    printf("  %f %f\n",rb[i],rnb[i]);
    if (rb[i] > 1.0)
      fb += 1.0;
    if (rnb[i] < 1.0)
      fnb += 1.0;
  }
  fb /= (float)m;
  fnb /= (float)m;

  if (fb >= 0.9){
    printf("  Burst\n");
    *rcat = strdup("b");
    p = mean_farray(rb,m);
  }else if (fnb >= 0.9){
    printf("  Non-burst\n");
    *rcat = strdup("nb");
    p = mean_farray(rnb,m);
  }else if (fb >= 0.5){
    printf("  Mixed, burst\n");
    *rcat = strdup("mb");
    p = mean_farray(rb,m);
  }else{
    printf("  Mixed, non-burst\n");
    *rcat = strdup("mnb");
    p = mean_farray(rnb,m);
  }
  
  printf("  p = %.4f\n",p);
  *rp = p; 
  
  myfree(mean); myfree(t); myfree(rb); myfree(rnb);
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_TUNING_TIME                             */
/*                                                                           */
/*  Compute the first time (in sampling units) when there is a               */
/*  significant difference in the responses in the two sarrays.              */
/*                                                                           */
/*****************************************************************************/
int get_tuning_time(s1,cnt1,n1,s2,cnt2,n2,start,period,k,sigma,signif)
     int **s1,*cnt1,n1,**s2,*cnt2,n2; /* spike data, counts */
     int start,period,k;
     float sigma,signif; /* smoothing for signif. array; signif. level */
{
  int i;
  int *psth1,*psth2,nsig,ttime;
  float *f1,*f2,t,*sig,*smsig;
  
  printf("  GET_TUNING_TIME\n");

  psth1 = simple_pst_sarray(s1,cnt1,n1,start,period);
  psth2 = simple_pst_sarray(s2,cnt2,n2,start,period);
  f1 = i2farray(psth1,period);
  f2 = i2farray(psth2,period);
  for(i=0;i<period;i+=k){     /*** Add in extra spikes in both ***/
    f1[i] += 10.0/(float)n1;  /*** To be sure we have non-zero variance ***/
    f2[i] += 10.0/(float)n1;  /*** for the NumRec routine. ***/
  }

  /*** Numerical Recipes test for different mean ***/
  nsig = period - k + 1;
  sig = (float *)myalloc(nsig*sizeof(float));
  for(i=0;i<nsig;i++)
    tutest(f1-1+i,k,f2-1+i,k,&t,&sig[i]);
  /*write_farray_plot("zz.sig",sig,nsig);*/
  smsig = smooth_with_gaussian(sig,nsig,sigma,0.01);
  /*write_farray_plot("zz.smsig",smsig,nsig);*/

  ttime = -1;
  for(i=0;i<nsig;i++)
    if ((ttime == -1)&&(smsig[i] < signif))
      ttime = i;

  return ttime;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_RESPONSE_ONSET                             */
/*                                                                           */
/*  Compute the first time (in sampling units) when there is a               */
/*  response above (or below) the background level.  The background          */
/*  level is established from the period starting at "bstart" last-          */
/*  ing "bperiod".                                                           */
/*                                                                           */
/*****************************************************************************/
int get_response_onset(s,cnt,n,start,period,bstart,bperiod,k,sigma,signif)
     int **s,*cnt,n,start,period,bstart,bperiod,k;
     float sigma,signif; /* smoothing for signif. array; signif. level */
{
  int i,j;
  int *psth,nsig,ttime,r;
  float *fpsth,*sig,*smsig,bg;
  
  printf("  GET_RESPONSE_ONSET\n");

  psth = simple_pst_sarray(s,cnt,n,start,period);
  fpsth = i2farray(psth,period);
  bg = 0.0;
  for(i=0;i<bperiod;i++)
    bg += fpsth[bstart+i];
  bg += 10.0/(float)n; /*** ADD IN A FEW SPIKES ***/
  bg /= (float)(bperiod*n);

  /*** Numerical Recipes test for different mean ***/
  nsig = period - k + 1;
  sig = (float *)myalloc(nsig*sizeof(float));
  for(i=0;i<nsig;i++){
    r = 0;
    for(j=0;j<k;j++)
      r += psth[i+j];
    sig[i] = binomial_signif(r,n*k,bg);
  }
  /*write_farray_plot("zon.sig",sig,nsig);*/
  smsig = smooth_with_gaussian(sig,nsig,sigma,0.01);
  /*write_farray_plot("zon.smsig",smsig,nsig);*/

  ttime = -1;
  for(i=0;i<nsig;i++)
    if ((ttime == -1)&&(smsig[i] < signif))
      ttime = i;

  return ttime;
}
/**************************************-**************************************/
/*                                                                           */
/*                         SHIFTED_DIFFERENCE_FLOATS                         */
/*                                                                           */
/*  Compute something like the Hamming distance between two float            */
/*  arrays for all shifts.  This is like cross correlation, but we           */
/*  compute the absolute value of the difference rather than the             */
/*  product.                                                                 */
/*                                                                           */
/*****************************************************************************/
float *shifted_difference_floats(data1,data2,n)
     float *data1,*data2;
     int n;
{
  int i,j;
  int start,stop;
  float *dist;
  int dn;
  
  dn = 2*n-1;
  dist = get_zero_farray(dn);

  for(i=-(n-1);i<=(n-1);i++){
    start = 0;
    stop = n-1;
    if (i>0)
      start = i;
    if (i<0)
      stop = n+i-1;
    for(j=start;j<=stop;j++)
      dist[i+n-1] += fabs(data1[j]-data2[j-i]);
    dist[i+n-1] /= (float)(stop-start+1);
  }
  return dist;
}
/**************************************-**************************************/
/*                                                                           */
/*                              SPIKES_DISTANCE                              */
/*                                                                           */
/*  Compute something like the Hamming distance between spike trains         */
/*  after smoothing with a Gaussian of standard deviation "sigma".           */
/*                                                                           */
/*****************************************************************************/
float *distance_trial(spikes1,n1,spikes2,n2,start,period,sigma)
     int *spikes1,n1,*spikes2,n2;
     int start,period;
     float sigma;
{
  float *s1,*s2;   /* expanded spike arrays */
  float *sm1,*sm2; /* smoothed spike arrays */
  float *dist;

  printf("  SPIKES_DISTANCE\n");

  s1 = expand_spike_array(spikes1,n1,start,period);
  s2 = expand_spike_array(spikes2,n2,start,period);

  sm1 = smooth_with_gaussian(s1,period,sigma,0.01);
  sm2 = smooth_with_gaussian(s2,period,sigma,0.01);

  myfree(s1);
  myfree(s2);

  dist = shifted_difference_floats(sm1,sm2,period);

  myfree(sm1);
  myfree(sm2);

  return dist;
}
/**************************************-**************************************/
/*                                                                           */
/*                              LOG_PROB_SPIKES                              */
/*                                                                           */
/*    *** THIS WORKS, BUT ISN'T REALLY WHAT WE WANT ***                      */
/*                                                                           */
/*  Compute the probability of the spike train based on the array            */
/*  "prob" which gives the probability of firing a spike at a given          */
/*  time.                                                                    */
/*                                                                           */
/*****************************************************************************/
float log_prob_spikes(data,n,prob,duration)
     int *data,n;
     float *prob;
     int duration;
{
  int i;
  float p;

  p = 0.0;
  for (i=0;i<n;i++)
    p += log(prob[data[i]]);

  return p;
}
/**************************************-**************************************/
/*                                                                           */
/*                              LOG_PROB_SARRAY                              */
/*                                                                           */
/*  Compute the probability of the kth spike train among the set of          */
/*  spike trains.                                                            */
/*                                                                           */
/*  *** WYETH - kill this routine at some time in future.                    */
/*                                                                           */
/*****************************************************************************/
float *log_prob_sarray(spikes,count,n,start,period)
     int **spikes,*count,n; /* spike data, counts */
     int start,period;      /* in sampling units */
{
  int i,j;
  int t;
  /*struct pst_struct *pst;*/
  float *prob,*logp;

  exit_error("LOG_PROB_SARRAY","No longer maintained");

  /***pst = pst_sarray(spikes,count,n,start,period,0.0,1);***/

  prob = (float *)myalloc(period*sizeof(float));
  logp = (float *)myalloc(n*sizeof(float));

  for (i=0;i<n;i++){
    for (j=0;j<period;j++)      /* copy counts from PST */
      prob[j] = 0.0; /*** pst->bin[j]; ***/
    for (j=0;j<count[i];j++){   /* subtract spikes of current trial */
      t = spikes[i][j] - start;
      if ((t>=0)&&(t<period))
	prob[t] -= 1.0;
    }
    for (j=0;j<period;j++){      /* normalize */
      prob[j] /= (float)(n-1);
      if (prob[j] == 0.0)         /* do not allow p=0 or p=1 */
	prob[j] = 1.0/(2.0*(float)(n-1));
      if (prob[j] == 1.0)
	prob[j] = 1.0 - 1.0/(2.0*(float)(n-1));
    }
    logp[i] = log_prob_spikes(spikes[i],count[i],prob,period);
  }
  return logp;
}
/**************************************-**************************************/
/*                                                                           */
/*                            COMPARE_SPIKE_TO_PST                           */
/*                                                                           */
/*  Compute the integral of the product of a spike train with the            */
/*  given PST (any float array, actually) as a function of the width         */
/*  of the spikes.  The PST is assumed to have its average value             */
/*  at points beyond the range of the PST array.                             */
/*                                                                           */
/*  Spikes are now represented as boxcar functions of varying width.         */
/*  Note:  widths for boxcar functions should be odd.                        */
/*                                                                           */
/*****************************************************************************/
float *compare_spike_to_pst(s,cnt,pst,n,width_list,nwidth)
     int *s,cnt;
     float *pst;
     int n,*width_list,nwidth;
{
  int i,j,k;
  int width,mid,pos;
  float *sum;
  float mean,sdev;
  double dsum,height;
  
  mean_sdev_farray(pst,n,&mean,&sdev);
  
  sum = get_zero_farray(nwidth);
  for(i=0;i<nwidth;i++){ /*** For each width. ***/
    width = width_list[i];
    mid = (width - 1)/2;
    height = 1.0/(float)width;
    dsum = 0.0;
    for(j=0;j<cnt;j++) /*** For each spike. ***/
      for(k=0;k<width;k++){
	pos = s[j] - mid + k; /*** Spike time - middle + index ***/
	if ((pos >= 0) && (pos < n))
	  dsum += height * pst[pos];
	else
	  dsum += height * mean;
      }
    sum[i] = (float)dsum;
  }
  return sum;
}
/**************************************-**************************************/
/*                                                                           */
/*                               AVG_ISI_TRIAL                               */
/*                                                                           */
/*  Compute average ISI in this trial (in sampling units).                   */
/*                                                                           */
/*****************************************************************************/
float avg_isi_trial(data,n,start,period)
     int *data,n;
     int start,period;
{
  int i;
  int count;
  float mean;
  
  count = 0;
  mean = 0.0;
  for(i=0;i<(n-1);i++)
    if ((data[i] >= start) && (data[i+1] < (start+period))){
      mean += (float)(data[i+1] - data[i]);
      count += 1;
    }
  mean /= (float)count;
  
  return mean;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_TRIAL_INTERVALS                            */
/*                                                                           */
/*  Return the float array containing the interval length at each            */
/*  sampling unit.  The interval length is the number of sampling            */
/*  units between the two spikes that surround the given sampling            */
/*  point.  We assign the length of the interval between "t1" and            */
/*  "t2" to the sampling units [t1,t2).  Before the first spike and          */
/*  at and after the last spike, -1 is used to indicate NO DATA.             */
/*                                                                           */
/*  Return "min" and "max", the first and last positions that are            */
/*  not -0.0001.                                                             */
/*                                                                           */
/*****************************************************************************/
float *get_trial_intervals(data,n,start,period,pmin,pmax)
     int *data,n;
     int start,period;
     int *pmin,*pmax;
{
  int i,j;
  float *interval,t;
  int t1,t2,min,max;

  interval = (float *)myalloc(period*sizeof(float));
  for(i=0;i<period;i++)
    interval[i] = -0.0001;

  min = period;
  max = -1;
  for (i=0;i<n-1;i++){
    t1 = data[i] - start;
    t2 = data[i+1] - start;
    t = (float)(t2-t1);
    if (t == 0.0) exit_error("GET_TRIAL_INTERVALS","Zero length interval");

    for(j=t1;j<t2;j++)
      if ((j>=0)&&(j<period)){
	interval[j] = t;
	if (j < min)
	  min = j;
	if (j > max)
	  max = j;
      }
  }
  *pmin = min;
  *pmax = max;
  return interval;
}
/**************************************-**************************************/
/*                                                                           */
/*                          INTERVAL_THRESHOLD                               */
/*                                                                           */
/*  This routine is a special purpose thresholding algorithm that            */
/*  operates on an farray that represents a spike train in its               */
/*  "interval" format (see "GET_TRIAL_INTERVALS").                           */
/*                                                                           */
/*  The mean of values below the threshold is computed.  This mean           */
/*  is substituted for those values.  If "flag" = 1, the mean is             */
/*  substituted for the above threshold values.                              */
/*                                                                           */
/*  "flag" = 0   Remove fluctuation in subthreshold values.                  */
/*  "flag" = 1   Set all above-threshold values to the subthreshold          */
/*               mean.                                                       */
/*                                                                           */
/*  *** This routine modifies "data".                                        */
/*                                                                           */
/*****************************************************************************/
void interval_threshold(data,n,thresh,flag)
     float *data;
     int n;
     float thresh;
     int flag;
{
  int i,k;
  float mean;

  k = 0;          /*** Compute the subthreshold mean. ***/
  mean = 0.0;
  for(i=0;i<n;i++)
    if (data[i] < thresh){
      mean += data[i];
      k += 1;
    }
  mean /= (float)k;
  
  for(i=0;i<n;i++)
    if (flag == 1){
      if (data[i] >= thresh)
	data[i] = mean;
    }else{
      if (data[i] < thresh)
	data[i] = mean;
    }
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_INTERVALS_XY_PLOT                            */
/*                                                                           */
/*  Return a list of (x,y) coordinates that make a plot of the inter-        */
/*  val representation of the spike train.                                   */
/*                                                                           */
/*****************************************************************************/
void get_intervals_xy_plot(data,n,start,period,rx,ry,rn)
     int *data,n;
     int start,period;
     int **rx,**ry,*rn;
{
  int i,k;
  int *x,*y,maxn,pmin,pmax;
  float *interv;

  maxn = 3*n; /*** There should only be about 2*n points ***/
  x = (int *)myalloc(maxn*sizeof(int));
  y = (int *)myalloc(maxn*sizeof(int));

  interv = get_trial_intervals(data,n,start,period,&pmin,&pmax);

  k = 0;
  x[k] = start;
  y[k] = interv[k];
  k = 1;
  for(i=1;i<period;i++){ /*** NOTE, This can give redundant point ***/
    if (interv[i] != interv[i-1]){
      x[k] = i-1+start;
      y[k] = (int)(interv[i-1] + 0.5);
      k += 1;
      x[k] = i+start;
      y[k] = (int)(interv[i] + 0.5);
      k += 1;
    }
  }
  x[k] = start+period-1;
  y[k] = (int)(interv[period-1] + 0.5);
  k += 1;

  myfree(interv);
  *rx = x; *ry = y; *rn = k;
}
/**************************************-**************************************/
/*                                                                           */
/*                           ISI_MEAN_SDEV_SKEW_TRIAL                        */
/*                                                                           */
/*  Compute the mean, standard deviation, and skewness of the ISI distrib.   */
/*                                                                           */
/*****************************************************************************/
void isi_mean_sdev_skew_trial(data,n,start,period,rm1,rm2,rm3,rnused)
     int *data,n;
     int start,period;
     float *rm1,*rm2,*rm3;
     int *rnused;
{
  int i,k;
  int spp,t0,t1;
  float *idata;

  if (!is_nondecreasing_iarray(data,n))
    exit_error("ISI_MOMENTS_TRIAL","Interval list not non-decreasing");
  
  idata = (float *)myalloc(n*sizeof(float));
  spp = start + period;
  k = 0;
  for(i=0;i<(n-1);i++){ /*** Use intervals with both spikes in window ***/
    t0 = data[i];
    t1 = data[i+1];
    if ((t0 >= start) && (t1 < spp)){
      idata[k] = (float)(t1-t0);
      k += 1;
    }
  }

  *rnused = k;
  if (k >= 2)
    mean_sdev_skew_farray(idata,k,rm1,rm2,rm3);
  else
    *rm1 = *rm2 = *rm3 = -1.0;
  
  myfree(idata);
}
/**************************************-**************************************/
/*                                                                           */
/*                       AVG_RATE_INTERVALS_TRIAL                            */
/*                                                                           */
/*  Compute average spike rate based on intervals rather than spike          */
/*  occurrence times.                                                        */
/*                                                                           */
/*  It appears that we assign the length of the interval between             */
/*  "t1" and "t2" to the time interval [t1,t2).                              */
/*                                                                           */
/*****************************************************************************/
float *avg_rate_intervals_trial(data,n,start,period,sampling)
     int *data,n;
     int start,period;
     float sampling;
{
  int i,j;
  float *avg,rate;
  int pos;
  int t1,t2;

  avg = get_zero_farray(period);  /* if no spikes, all zeros are returned */

  pos = 0;
  rate = 0.0;
  for (i=0;i<n-1;i++){
    t1 = data[i] - start;
    t2 = data[i+1] - start;
    if ((t2 >= 0)&&(t2 <= period)){
      if ((t2-t1) == 0)
	exit_error("AVG_RATE_INTERVALS_SARRAY","Zero interval");
      rate = sampling/(float)(t2-t1);
      for (j=pos;j<t2;j++)
	avg[j] += rate; /* average the spike rates, not intervals */
      pos = t2;
    }
  }
  for (j=pos;j<period;j++)
    avg[j] += rate; /* fill in using the last value of rate */
  pos = period;

  return avg;
}
/**************************************-**************************************/
/*                                                                           */
/*                      AVG_RATE_INTERVALS_SARRAY                            */
/*                                                                           */
/*  Compute average spike rate based on intervals rather than spike          */
/*  occurrence times.                                                        */
/*                                                                           */
/*****************************************************************************/
float *avg_rate_intervals_sarray(spikes,count,n,start,period,sampling)
     int **spikes,*count,n; /* spike data, counts */
     int start,period;      /* in sampling units */
     float sampling;
{
  int i,j,k;
  float *avg,rate;
  int pos,t1,t2;

  avg = get_zero_farray(period);  /* if no spikes, all zeros are returned */

  for (i=0;i<n;i++){
    pos = 0;
    rate = 0.0;
    for (j=0;j<count[i]-1;j++){
      t1 = spikes[i][j] - start;
      t2 = spikes[i][j+1] - start;
      if ((t2 >= 0)&&(t2 <= period)){
	if ((t2-t1) == 0)
	  exit_error("AVG_RATE_INTERVALS_SARRAY","Zero interval");
	rate = sampling/(float)(t2-t1);
	for (k=pos;k<t2;k++)
	  avg[k] += rate; /* average the spike rates, not intervals */
	pos = t2;
      }
    }
    for (k=pos;k<period;k++)
      avg[k] += rate; /* fill in using the last value of rate */
    pos = period;
  }
  for (i=0;i<period;i++)
    avg[i] /= (float)n;
  return avg;
}
/**************************************-**************************************/
/*                                                                           */
/*                             INTERVAL_CV_SARRAY                            */
/*                                                                           */
/*  Compute the coefficient of variation (CV) for the distribution           */
/*  of ISIs occuring at a particular sampling time over all trials.          */
/*  In other words, if time during a trial runs horizontally and             */
/*  trials are arranged vertically, compute the CV by averaging              */
/*  vertically over ISIs for each discrete time point.                       */
/*                                                                           */
/*  Intervals which are "open" at the beginning and end of the trial         */
/*  are not included in the average.                                         */
/*                                                                           */
/*****************************************************************************/
float *interval_cv_sarray(s,cnt,n,period)
     int **s,*cnt,n; /* spike data, counts */
     int period; /* start is assumed to be zero */
{
  int i,j;
  int pos;
  float *cv;
  float *data,mean,sdev;
  int *ptr; /* point to current spike in train */

  printf("INTERVAL_CV_SARRAY\n");
  printf("  period = %d\n",period);

  data = (float *)myalloc(n*sizeof(float));
  cv = (float *)myalloc(period*sizeof(float));
  ptr = (int *)myalloc(n*sizeof(int));
  for (i=0;i<n;i++)
    ptr[i] = -1;

  for(i=0;i<period;i++){
    pos = 0;
    for(j=0;j<n;j++)
      if ((ptr[j]+1) < cnt[j]){ /* if not done with jth trial */
	if (i >= s[j][ptr[j]+1]) /* adjust jth pointer for time i */
	  ptr[j] = ptr[j] + 1;
	if ((ptr[j]>=0) && ((ptr[j]+1) < cnt[j])){
	  data[pos] = 1.0/(float)(s[j][ptr[j]+1] - s[j][ptr[j]]);
	  pos += 1;
	}
      }
    mean_sdev_farray(data,pos,&mean,&sdev);
    cv[i] = sdev/mean;
  }
  myfree(data);
  myfree(ptr);

  return cv;
}
/**************************************-**************************************/
/*                                                                           */
/*                             AVG_TRIAL_CV_SARRAY                           */
/*                                                                           */
/*  Compute the coefficient of variation (CV) for the distribution of ISIs   */
/*  on a trial-by-trial basis, and return the average of these values.       */
/*                                                                           */
/*  I have arbitrarily chosen 3 ISIs to be the minimum per trial for a       */
/*  valid CV computation, i.e. trials with less than 3 ISIs are ignored.     */
/*                                                                           */
/*****************************************************************************/
void avg_trial_cv_sarray(s,cnt,n,start,period,cv_mean,cv_sdev)
     int **s,*cnt,n;
     int start,period;
     float *cv_mean,*cv_sdev;
{
  int i,j;
  int max,ipos,cpos,count,t,spp;
  float *isi_data,*cv_data;
  float mean,sdev;

  max = max_of_iarray(cnt,n);

  isi_data = (float *)myalloc(max*sizeof(float));
  cv_data = (float *)myalloc(n*sizeof(float));

  spp = start + period;
  cpos = 0;
  count = 0;
  for(i=0;i<n;i++){
    ipos = 0;
    for(j=0;j<(cnt[i]-1);j++){
      t = s[i][j];
      if ((t >= start) && (t < spp)){
	isi_data[ipos] = (float)(s[i][j+1]-t);
	ipos += 1;
      }
    }
    if (ipos > 2){
      mean_sdev_farray(isi_data,ipos,&mean,&sdev);
      cv_data[cpos] = sdev/mean;
      cpos += 1;
    }else
      count += 1;
  }
  if (count > 0)
    printf("  *** %d trials not included:  less than 3 ISIs.\n",count);
  mean_sdev_farray(cv_data,cpos,cv_mean,cv_sdev);

  myfree(isi_data);
  myfree(cv_data);
}
/**************************************-**************************************/
/*                                                                           */
/*                             AVG_TRIAL_CV2_SARRAY                          */
/*                                                                           */
/*  Compute the average value of CV2 for each trial, and return the mean     */
/*  and std. dev. of this value.  For consecutive ISIs a and b,              */
/*                                                                           */
/*                   2 * |a-b|                                               */
/*             CV2 = ---------.                                              */
/*                     a+b                                                   */
/*                                                                           */
/*  This measure was proposed by Gary Holt in Klab.  It should be less       */
/*  sensitive to slow (relative to the length of an ISI) rate variations     */
/*  during the trial.                                                        */
/*                                                                           */
/*****************************************************************************/
void avg_trial_cv2_sarray(s,cnt,n,start,period,cv_mean,cv_sdev)
     int **s,*cnt,n;
     int start,period;
     float *cv_mean,*cv_sdev;
{
  int i,j,k;
  int max,nisi,count,t,spp;
  float *isi_data,*cv_data,*cv;
  float mean,sdev;

  max = max_of_iarray(cnt,n);

  isi_data = (float *)myalloc(max*sizeof(float));
  cv_data = (float *)myalloc(max*sizeof(float));
  cv = (float *)myalloc(n*sizeof(float));

  spp = start + period;
  k = 0;
  count = 0;
  for(i=0;i<n;i++){
    nisi = 0;
    for(j=0;j<(cnt[i]-1);j++){
      t = s[i][j];
      if ((t >= start) && (t < spp)){
	isi_data[nisi] = (float)(s[i][j+1]-t);
	nisi += 1;
      }
    }
    if (nisi > 2){
      for(j=0;j<(nisi-1);j++)
	cv_data[j] = 2.0*fabs(isi_data[j]-isi_data[j+1])/
	  (float)(isi_data[j]+isi_data[j+1]);
      mean_sdev_farray(cv_data,(nisi-1),&mean,&sdev);
      cv[k] = mean;
      k += 1;
    }else
      count += 1;
  }
  if (count > 0)
    printf("  *** %d trials not included:  less than 3 ISIs.\n",count);
  mean_sdev_farray(cv,k,cv_mean,cv_sdev);

  myfree(isi_data); myfree(cv_data); myfree(cv);
}
/**************************************-**************************************/
/*                                                                           */
/*                              GET_GAP_SARRAY                               */
/*                                                                           */
/*  Return an sarray containing the start times of ISI's that are            */
/*  between (including) the specified sizes.                                 */
/*                                                                           */
/*****************************************************************************/
void get_gap_sarray(s,cnt,n,min,max,rs,rcnt)
     int **s,*cnt,n,min,max;
     int ***rs,**rcnt;
{
  int i,j,k;
  int **sg,*cntg,*stemp,tn,t;

  sg = (int **)myalloc(n*sizeof(int *));
  cntg = (int *)myalloc(n*sizeof(int));

  tn = max_of_iarray(cnt,n);
  stemp = (int *)myalloc(tn*sizeof(int));
  for (i=0;i<n;i++){ // for each trial
    k = 0;
    for (j=0;j<(cnt[i]-1);j++){
      t = s[i][j+1] - s[i][j];
      if ((t>=min)&&(t<=max)){
	stemp[k] = s[i][j];
	k += 1;
      }
    }
    sg[i] = copy_iarray(stemp,k);
    cntg[i] = k;
  }
  myfree(stemp);
  *rs = sg; *rcnt = cntg;
} 
/**************************************-**************************************/
/*                                                                           */
/*                       FIRST_LAST_COUNT_EVENTS_TRIAL                       */
/*                                                                           */
/*   Get the index in the spike array of the first spike of the first        */
/*   "event" and the last spike of the last "event" to consider.             */
/*   All events should be totally contained within the specified             */
/*   duration, i.e., they should be separated from the ends by               */
/*   at least "d" blank sampling units, or be separated from other           */
/*   spikes beyond the ends by at least "d" blank sampling units.            */
/*                                                                           */
/*****************************************************************************/
void first_last_count_events_trial(s,n,d,start,duration,first,last,count)
     int *s,n,d;
     int start,duration;
     int *first,*last,*count;
{
  int i;
  int previous,t;

  *count = 0;
  *first = -1; /* These are indexes in the spike array. */
  *last = -1;
  previous = start-1;
  for(i=0;i<n;i++){ /*** FIND FIRST INDEX ***/
    t = s[i];
    if ((t >= start) && (t < (start+duration)) && (t > (d+previous)) &&
	(*first < 0))
      *first = i;
    else
      previous = t;
  }
  if (*first >= 0){
    previous = start+duration;
    for(i=(n-1);i>=0;i--){ /*** FIND LAST INDEX ***/
      t = s[i];
      if ((t >= start) && (t < (start+duration)) && (t < (previous-d)) &&
	  (*last < 0))
	*last = i;
      else
	previous = t;
    }
    if (*first <= *last){
      *count = 1;
      for(i=*first;i<*last;i++)
	if ((s[i+1] - s[i]) > d)
	  *count += 1;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_EVENTS_SARRAY                             */
/*                                                                           */
/*  Return an "event_trial" structure for the given spike trains.            */
/*  Notice that only events totally isolated within the given spike          */
/*  trains are used (see "first_last_count_events").  Events are             */
/*  defined as the maximum length trains of consecutive spikes with          */
/*  no intervals longer than "d" sampling units.                             */
/*                                                                           */
/*  Modified by WYETH on Thu Mar 24 03:13:20 PST 1994.                       */
/*                                                                           */
/*****************************************************************************/
struct event_trial *get_events_sarray(s,cnt,n,start,duration,d)
     int **s,*cnt,n,start,duration;
     int d; // max within-burst interval in sampling units
{
  int i,j,k;
  int ecount;     /* number of events processed in this trial */
  int scount;     /* number of spikes processed in this event */
  int end_burst;  /* flag indicating end of burst */
  int first,last,count;
  struct event_trial *etrial;

  etrial = (struct event_trial *)myalloc(n*sizeof(struct event_trial));

  for (i=0;i<n;i++){ /* for each trial */
    first_last_count_events_trial(s[i],cnt[i],d,start,duration,
				  &first,&last,&count);
    etrial[i].nev = count;
    etrial[i].ev = (struct event_struct *)
      myalloc(count*sizeof(struct event_struct));
    if (count > 0){
      k = first;
      ecount = 0;
      scount = 0;
      while (k<=last){ /* process the kth spike in this (ith) trial */
	scount += 1;
	k+=1;
	end_burst = 0;
	while (!end_burst) /* until the end of this event */
	  if (k>last) /* at end of valid spikes */
	    end_burst = 1;
	  else if (s[i][k] > (d+s[i][k-1])) /* interval too long */
	    end_burst = 1;
	  else{ /* count this spike in the event */
	    k += 1;
	    scount += 1;
	  }
	etrial[i].ev[ecount].n = scount;
	etrial[i].ev[ecount].t = (int *)myalloc(scount*sizeof(int));
	for(j=0;j<scount;j++)
	  etrial[i].ev[ecount].t[j] = s[i][k-scount+j];
	ecount += 1;
	scount = 0;
      }
    }
  }
  return etrial;
}
/**************************************-**************************************/
/*                                                                           */
/*                                COUNT_EVENTS                               */
/*                                                                           */
/*****************************************************************************/
int count_events(data,n,start,period,d)
     int *data,n,start,period,d;
{
  int i,k;
  int count,done;
  
  k = 0;
  done = 0;
  while (!done){
    if (k >= n)
      done = 1;
    else if ((data[k] >= start) && (data[k] < start+period))
      done = 1;
    else
      k += 1;
  }
  
  count = 0;
  if (k < n){
    for (i=k;i<(n-1);i++) /* count inter-event pauses */
      if (data[i+1] < (start+period))
	if (data[i+1] - data[i] > d)
	  count += 1;
    count += 1;
  }
  
  return count;
}
/**************************************-**************************************/
/*                                                                           */
/*                         COUNT_WEIGHTED_EVENTS                             */
/*                                                                           */
/*  List of "option" values:                                                 */
/*   -1  count multi-spike events only                                       */
/*    0  count number of spikes in event (equiv to spike count)              */
/*    1  count each event as 1                                               */
/*    2  count square of number of spikes in event                           */
/*    3  count square root of spikes in event                                */
/*                                                                           */
/*****************************************************************************/
float count_weighted_events(data,n,start,period,d,option)
     int *data,n,start,period,d,option;
{
  int i,k;
  int done,spikes;
  float weight;
  
  k = 0;
  done = 0;
  while (!done){
    if (k >= n)
      done = 1;
    else if ((data[k] >= start) && (data[k] < start+period))
      done = 1;
    else
      k += 1;
  }
  
  weight = 0.0;
  if (k < n){
    spikes = 1;
    for (i=k;i<(n-1);i++) /* count inter-event pauses */
      if (data[i+1] < (start+period)){
	if (data[i+1] - data[i] > d){
	  if ((option==-1)&&(spikes > 1))
	    weight += 1.0;
	  else if (option==0)
	    weight += (float)spikes;
	  else if (option==1)
	    weight += 1.0;
	  else if (option==2)
	    weight += (float)(spikes*spikes);
	  else if (option==3)
	    weight += sqrt((float)spikes);

	  spikes = 1;
	}else
	  spikes += 1; /* another spike in this burst */
      }
    
    if ((option==-1)&&(spikes > 1)) /* PROCESS FINAL EVENT */
      weight += 1.0;
    else if (option==0)
      weight += (float)spikes;
    else if (option==1)
      weight += 1.0;
    else if (option==2)
      weight += (float)(spikes*spikes);
    else if (option==3)
      weight += sqrt((float)spikes);
  }
  return weight;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MAKE_POISSON_FLOAT_SPIKES                         */
/*                                                                           */
/*  This routine may produce spike trains that have multiple copies of a     */
/*  given occurrence time (after rounding to sampling units).                */
/*                                                                           */
/*  lambda - firing rate (spikes/sec)                                        */
/*  rmu - Mean of Gaussian refractory period                                 */
/*  rsigma - SD of Gaussian refractory period                                */
/*  k - place a spike only on every "k" poisson point.                       */
/*                                                                           */
/*  Notes:                                                                   */
/*  - Refractory time is chosen until a non-negative value is returned.      */
/*  - "seed" should point to a negative number to initialize random          */
/*    sequence.  If seed is passed subsequently, it can be the positive      */
/*    value unchanged from the previous call.  BE CAREFUL that no other      */
/*    routines are executed between calls that might make use of the         */
/*    "spike_util_ran2" procedure.                                           */
/*                                                                           */
/*****************************************************************************/
void make_poisson_float_spikes(data,n,start,duration,sampling,lambda,rmu,
			       rsigma,k,seed)
     float **data;
     int *n;
     int start,duration;    // (sampling units)
     float sampling;        // samples per second
     float lambda;          // firing rate (spikes/s)
     float rmu,rsigma;      // refractory mean and SD (s)
     int k;                 // store kth spike only
     int *seed;             // random seed
{  
  int i;
  float *parray,ftime,grefract,rmus,rsigs;
  int spike_count,max_spikes;

  /*
    printf("  MAKE_POISSON_FLOAT_SPIKES\n");
    printf("    start %d\n",start);
    printf("    duration %d\n",duration);
    printf("    sampling %f\n",sampling);
    printf("    lambda %f\n",lambda);
    printf("    refr_mu %f\n",rmu);
    printf("    refr_sigma %f\n",rsigma);
    printf("    k %d\n",k);
    printf("    seed %d\n",*seed);
  */

  if (lambda == 0.0){
    *n = 0;
    *data = NULL;
    return;
  }else if (lambda < 0.0){
    exit_error("MAKE_POISSON_FLOAT_SPIKES","Lambda is negative");
  }

  rmus = rmu*sampling;      // Convert to sampling units
  rsigs = rsigma*sampling;  // Convert to sampling units

  max_spikes = 10*duration; // an arbitrary limit
  parray = (float *)myalloc(max_spikes*sizeof(float));
  spike_count = 0;

  // generate first interval (sum of k intervals)
  ftime = 0.0;
  for (i=0;i<k;i++)
    ftime += su_expdev(lambda/sampling,seed);

  //while ((int)rint(ftime) < duration){
  //
  //  WYETH ******   Something changed here, Sept 2007, low lambda created
  //                 long spike times that were negative when rounded to int
  //                 using line above.  Don't know why this happened?
  //
  //  WYETH ******   Something changed here, Sept 2007, low lambda created
  //                 long spike times that were negative when rounded to int
  //                 using line above.  Don't know why this happened?
  //
  while (ftime < (float)duration){
    //printf("rint(ftime) = %d  \n",(int)rint(ftime));
    if (spike_count >= max_spikes){
      printf("  lambda = %f\n",lambda);
      printf("  duration= %d,  max_spikes = %d  spike_count = %d\n",
	     duration,max_spikes,spike_count);
      exit_error("MAKE_POISSON_FLOAT_SPIKES","too many spikes");
    }
    parray[spike_count] = (float)start + ftime;
    spike_count += 1;
    // get next interval (sum of k intervals)
    for (i=0;i<k;i++){
      ftime += su_expdev(lambda/sampling,seed);
      if ((rsigma > 0.0) && (rmu > 0.0)){
	//grefract = su_gaussian_float(rmu*sampling,rsigma*sampling,seed);
	grefract = su_gaussian_float(rmus,rsigs,seed);
	while (grefract < 0.0)
	  grefract = su_gaussian_float(rmus,rsigs,seed);
	//grefract = su_gaussian_float(rmu*sampling,rsigma*sampling,seed);
	ftime += grefract;
      }
    }
  }
  *n = spike_count;
  *data = parray;
}
/**************************************-**************************************/
/*                                                                           */
/*                        MAKE_POISSON_FLOAT_SPIKES_NEW                      */
/*                                                                           */
/*  *** 'duration' is FLOAT here.                                            */
/*                                                                           */
/*  This routine may produce spike trains that have multiple copies of a     */
/*  given occurrence time (after rounding to sampling units).                */
/*                                                                           */
/*  lambda - firing rate (spikes/sec)                                        */
/*  rmu - Mean of Gaussian refractory period                                 */
/*  rsigma - SD of Gaussian refractory period                                */
/*  k - place a spike only on every "k" poisson point.                       */
/*                                                                           */
/*  Notes:                                                                   */
/*  - Refractory time is chosen until a non-negative value is returned.      */
/*  - "seed" should point to a negative number to initialize random          */
/*    sequence.  If seed is passed subsequently, it can be the positive      */
/*    value unchanged from the previous call.  BE CAREFUL that no other      */
/*    routines are executed between calls that might make use of the         */
/*    "spike_util_ran2" procedure.                                           */
/*                                                                           */
/*****************************************************************************/
void make_poisson_float_spikes_new(data,n,start,duration,sampling,lambda,rmu,
				   rsigma,k,seed)
     float **data;            // [*n] Returned spike times
     int *n;                  // Returned number of spikes
     float start,duration;    // sampling units
     float sampling;          // samples per second
     float lambda,rmu,rsigma; // time unit is SECONDS here
     int k;                   // store kth spike only
     int *seed;               // random seed
{  
  int i;
  float *parray,*ts,ftime,grefract;
  int spike_count,max_spikes;

  if (lambda == 0.0){
    *n = 0;
    *data = NULL;
    return;
  }else if (lambda < 0.0)
    exit_error("MAKE_POISSON_FLOAT_SPIKES_NEW","Lambda is negative");

  max_spikes = (int)(10.0 * lambda * duration/sampling); // 10x mean
  parray = (float *)myalloc(max_spikes*sizeof(float));
  spike_count = 0;

  // generate first interval (sum of k intervals)
  ftime = 0.0;
  for (i=0;i<k;i++)
    ftime += su_expdev(lambda/sampling,seed);

  while (ftime < duration){
    if (spike_count >= max_spikes){
      printf("  *** MAKE_POISSON_FLOAT_SPIKES:  too many spikes.  Exiting.\n");
      exit(0);
    }
    parray[spike_count] = start + ftime;
    spike_count += 1;
    // get next interval (sum of k intervals)
    for (i=0;i<k;i++){
      ftime += su_expdev(lambda/sampling,seed);
      if ((rsigma > 0.0) && (rmu > 0.0)){
	grefract = su_gaussian_float(rmu*sampling,rsigma*sampling,seed);
	while (grefract < 0.0)
	  grefract = su_gaussian_float(rmu*sampling,rsigma*sampling,seed);
	ftime += grefract;
      }
    }
  }
  *n = spike_count;
  *data = copy_farray(parray,spike_count);
  myfree(parray);
}
/**************************************-**************************************/
/*                                                                           */
/*                     MAKE_POISSON_SPIKES_INT_REFRACT                       */
/*                                                                           */
/*   INSTEAD, USE make_poisson_float_spikes ??                               */
/*                                                                           */
/*   Fill "data" with poisson distributed spikes having parameter            */
/*   "lambda" with a refractory period "refract".  The parameter             */
/*   "k" allows spikes to be placed at every "k"th poisson spike             */
/*   only, thus achieving gamma distributed spikes.  "refract"               */
/*   must be 1 or greater.                                                   */
/*                                                                           */
/*  *** The firing rate is correct only if refract == 0.                     */
/*                                                                           */
/*****************************************************************************/
int make_poisson_spikes_int_refract(data,period,lambda,sampling,refract,k,
				    seed)
     int *data,period;  // expanded spike array
     float lambda;      // Firing rate in spikes/sec
     float sampling;    // Samples per second
     int refract,k;
     int *seed;
{  
  int i;
  int time,n;
  float ftime;

  if (*seed < 0){
    *seed = -*seed;
    spike_util_ran2(seed);
  }

  // generate first interval (sum of k intervals)
  ftime = 0.0;
  for (i=0;i<k;i++)
    ftime += sampling*su_expdev(lambda,seed);
  time = my_rint(ftime);

  n = 0;
  while (time < period) {
    data[time] = 1;
    n += 1;
    while (my_rint(ftime) <= time+refract)
      for (i=0;i<k;i++)
	ftime += sampling*su_expdev(lambda,seed);
    time = my_rint(ftime);
  }
  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MAKE_GAUSSIAN_SPIKES                           */
/*                                                                           */
/*   Fill data array with spike train having Gaussian inter-spike-           */
/*   interval density function with mean "mu" and standard deviation         */
/*   sigma.  The first spike is uniformly distributed from 0 to "mu".        */
/*   "mu" must be less than "n".                                             */
/*                                                                           */
/*   *** Copied from "make_spikes.c" 9/13/96.                                */
/*                                                                           */
/*****************************************************************************/
int make_gaussian_spikes(data,n,mu,sigma,seed)
     int *data,n;
     float mu,sigma;
     int *seed;
{
  int i;
  int count,time,rnum;

  if (*seed < 0){
    *seed = -*seed;
    spike_util_ran2(seed);
  }

  /* clear all spikes */
  count = 0;
  for (i=0;i<n;i++) {
    if (data[i] == 1)
      count -= 1;
    data[i] = 0;
  }
  
  /* pick random starting time between 0 and mu, uniform */
  /* create spikes until time runs out */
  time = my_rint(spike_util_ran2(seed)*mu);
  while (time < n) {
    data[time] = 1;
    count += 1;
    rnum = gaussian_int(mu,sigma,seed);
    while (rnum <= 0)
      rnum = gaussian_int(mu,sigma,seed);
    time += rnum;
  }
  return count;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MAKE_UNIFORM_SPIKES                           */
/*                                                                           */
/*   Fill data array with spike train having Poisson characteristics.        */
/*   Use uniform density function to generate spikes.                        */
/*                                                                           */
/*   *** Copied from "make_spikes.c" 9/13/96.                                */
/*                                                                           */
/*****************************************************************************/
int make_uniform_spikes(data,n,dens,seed)
     int *data,n;
     float dens; /* probability of a spike in a msec bin */
     int *seed;
{
  int i;
  int count;

  if (*seed < 0){
    *seed = -*seed;
    spike_util_ran2(seed);
  }
  
  /* clear all spikes */
  count = 0;
  for (i=0;i<n;i++) {
    if (data[i] == 1)
      count -= 1;
    data[i] = 0;
  }
  
  for (i=0;i<n;i++)
    if (spike_util_ran2(seed) < dens){
      data[i] = 1;
      count += 1;
    }
  return count;
}
/**************************************-**************************************/
/*                                                                           */
/*                           SPIKE_TO_GAUSSIAN_BURST                         */
/*                                                                           */
/*   Replace each occurrence of a spike with a burst of spikes               */
/*   centered at the original spike.                                         */
/*                                                                           */
/*   *** Copied from "make_spikes.c" 9/13/96.                                */
/*                                                                           */
/*****************************************************************************/
int spike_to_gaussian_burst(data,n,width_mu,width_sigma,mu,sigma,seed)
     int *data,n;
     float width_mu,width_sigma; /* length of burst */
     float mu,sigma;             /* spike density within burst */
     int *seed;
{
  int i;
  int start,*old_data,count,burst_width;

  if (*seed < 0){
    *seed = -*seed;
    spike_util_ran2(seed);
  }
  
  old_data = (int *)myalloc(n*sizeof(int));
  for (i=0;i<n;i++) {
    old_data[i] = data[i];
    data[i] = 0;
  }

  count = 0;
  for (i=0;i<n;i++)
    if (old_data[i] == 1) {
      burst_width = gaussian_int(width_mu,width_sigma,seed);
      if (burst_width > 0) {
	start = i - burst_width / 2;
	if (start < 0) {burst_width += start; start = 0;}
	if ((start+burst_width) > n) burst_width = n-start;
	count += make_gaussian_spikes(data+start,burst_width,mu,sigma,
				      seed);
      }
    }
  myfree(old_data);
  return count;
}
/**************************************-**************************************/
/*                                                                           */
/*                           SPIKE_TO_UNIFORM_BURST                          */
/*                                                                           */
/*   Replace each occurrence of a spike with a burst of spikes               */
/*   centered at the original spike.  The burst width is Gaussian.           */
/*   The interburst spikes are poisson.                                      */
/*                                                                           */
/*   *** Copied from "make_spikes.c" 9/13/96.                                */
/*                                                                           */
/*****************************************************************************/
int spike_to_uniform_burst(data,n,dens,width_mu,width_sigma,seed)
     int *data,n;
     float dens; /* density of spikes within burst */
     float width_mu,width_sigma;  /* burst width mean and std.dev. */
     int *seed;
{
  int i;
  int start,*old_data,count,burst_width;
  
  if (*seed < 0){
    *seed = -*seed;
    spike_util_ran2(seed);
  }

  old_data = (int *)myalloc(n*sizeof(int));
  for (i=0;i<n;i++) {
    old_data[i] = data[i];
    data[i] = 0;
  }
  
  count = 0;
  for (i=0;i<n;i++)
    if (old_data[i] == 1) {
      burst_width = gaussian_int(width_mu,width_sigma,seed);
      if (burst_width > 0) {
	start = i - burst_width / 2;
	if (start < 0) {burst_width += start; start = 0;}
	if ((start+burst_width) > n) burst_width = n-start;
	count += make_uniform_spikes(data+start,burst_width,dens,seed);
      }
    }
  return count;
}
/**************************************-**************************************/
/*                                                                           */
/*                        SPIKE_UTIL_RANDOM_SELECT_SPIKES                    */
/*                                                                           */
/*  Create a new spike train by randomly selecting spikes in 's' with        */
/*  probability 'p'.                                                         */
/*                                                                           */
/*****************************************************************************/
void spike_util_random_select_spikes(s,n,p,seed,rs,rn)
     float *s;      /* Pick spikes from this train [n] */
     int n;
     float p;       /* Pick spikes with this probability */
     int seed;
     float **rs;    /* Return new storage with selected spikes [rn] */
     int *rn;
{
  int i,k;
  float *t;

  if (seed > 0)
    seed *= -1;

  t = (float *)myalloc(n*sizeof(float));
  k = 0;
  for(i=0;i<n;i++){
    if (spike_util_ran2(&seed) < p){
      t[k] = s[i];
      k += 1;
    }
  }

  *rn = k;
  *rs = copy_farray(t,k);
  myfree(t);
}
/**************************************-**************************************/
/*                                                                           */
/*                        MAKE_RANDOM_SARRAY_FROM_ISI                        */
/*                                                                           */
/*  Return an "sarray" in which ISIs are chosen randomly from the given isi  */
/*  probabilities "pisi".                                                    */
/*                                                                           */
/*****************************************************************************/
void make_random_sarray_from_isi(pisi,nisi,rs,rcnt,n,period,seed1,seed2)
     float *pisi;
     int nisi,***rs,**rcnt,n,period;
     int seed1,seed2;
{
  int i,j,k;
  int *isilist,imax,**s,*cnt,kisi,t,*tdata;

  isilist = (int *)myalloc(10000*sizeof(int));
  imax = 0;
  for(i=0;i<nisi;i++){
    k = (int)(10000.0 * pisi[i]);
    for(j=0;j<k;j++){
      isilist[imax] = i;
      imax += 1;
      if (imax > 10000)
	exit_error("MAKE_RANDOM_SARRAY_FROM_ISI","imax TOO LARGE");
    }
  }
  shuffle_iarray(isilist,imax,3,seed1);

  if (seed2 > 0) seed2 = -seed2;
  spike_util_ran2(&seed2);

  cnt = (int *)myalloc(n*sizeof(int));
  s = (int **)myalloc(n*sizeof(int *));
  tdata = (int *)myalloc(period*sizeof(int));
  for(i=0;i<n;i++){
    k = 0;
    t = 0;
    while(t<period){
      kisi = (int)(spike_util_ran2(&seed2) * (float)imax);
      if ((kisi < 0) || (kisi >= imax))
	exit_error("MAKE_RANDOM_SARRAY_FROM_ISI","kisi OUT OF RANGE");
      t += isilist[kisi];
      if (t < period){
	tdata[k] = t;
	k += 1;
      }
    }
    s[i] = copy_iarray(tdata,k);
    cnt[i] = k;
  }
  myfree(tdata);
  *rs = s;
  *rcnt = cnt;
}
/**************************************-**************************************/
/*                                                                           */
/*                        MAKE_RANDOM_SPIKES_FROM_PROB                       */
/*                                                                           */
/*  Return a spike train having spikes chosen according to the probability   */
/*  of firing sent in "prob".                                                */
/*                                                                           */
/*  *** NOTE:  "seed" should point to a negative value on the first call to  */
/*             initialize the random number generator.                       */
/*                                                                           */
/*****************************************************************************/
void make_random_spikes_from_prob(prob,period,spikes,n,seed)
     float *prob;
     int period,**spikes,*n;
     int *seed;
{
  int i,k;
  int *tdata;

  tdata = (int *)myalloc(period*sizeof(int));

  /************************************************************** TESTING
    fdata = (float *)myalloc(period*sizeof(float));
    k = 0;
    for(i=0;i<period;i++){
    fdata[i] = spike_util_ran2(seed);
    if (prob[i] > fdata[i]){
    tdata[k] = i;
    k += 1;
    }
    }
    mean_sdev_farray(fdata,period,&mean,&sdev);
    mean_sdev_farray(prob,period,&mp,&sdev);
    printf("muRan= %.4f  muProb= %.4f   n= %d\n",mean,mp,k);
    myfree(fdata);
    **********************************************************************/

  k = 0;
  for(i=0;i<period;i++){
    if (prob[i] > spike_util_ran2(seed)){
      tdata[k] = i;
      k += 1;
    }
  }
  *spikes = copy_iarray(tdata,k);
  *n = k;

  myfree(tdata);
}
/**************************************-**************************************/
/*                                                                           */
/*                    MAKE_RANDOM_SPIKES_FROM_PROB_REFRACT                   */
/*                                                                           */
/*  Return a spike train having spikes chosen according to the probability   */
/*  of firing sent in "prob".                                                */
/*                                                                           */
/*  *** NOTE:  "seed" should point to a negative value on the first call to  */
/*             initialize the random number generator.                       */
/*                                                                           */
/*  The refractory period is implemented by multiplying a 1-r(t) times       */
/*  the probability of firing, where r(t) is a Gaussian with SD "rsig" and   */
/*  maximum "rmax" <= 1, and r(t) is shifted relative to the last spike.     */
/*                                                                           */
/*****************************************************************************/
void make_random_spikes_from_prob_refract(prob,rsig,rmax,period,spikes,n,seed)
     float *prob,rsig,rmax;
     int period,**spikes,*n;
     int *seed;
{
  int i,k;
  int *tdata,nr,c,tr;
  float *r,p;

  discrete_gaussian(rsig,0.01,&r,&nr); // Get the refractory function
  c = (nr-1)/2; // Center index of Gaussian in "r" array
  if (rmax > 1.0)
    rmax = 1.0;
  make_max_const_farray(r,nr,rmax);
  //append_farray_plot("zzz.pl","mask",r,nr,1); exit(0);
  tdata = (int *)myalloc(period*sizeof(int));
  k = 0;
  tr = nr; // Time since last spike, set to maximum so no refract at first
  for(i=0;i<period;i++){
    p = prob[i];
    if ((c+tr) < nr)
      p *= 1.0 - r[c+tr];
    if (p > spike_util_ran2(seed)){
      tdata[k] = i;
      k += 1;
      tr = 0;
    }else
      tr += 1;
  }
  *spikes = copy_iarray(tdata,k);
  *n = k;

  myfree(tdata); myfree(r);
}
/**************************************-**************************************/
/*                                                                           */
/*                        MAKE_RANDOM_SARRAY_FROM_PROB                       */
/*                                                                           */
/*  Return an "sarray" that has spike chosen according to the probability    */
/*  of firing sent in "prob".  "n" is set before calling.                    */
/*                                                                           */
/*****************************************************************************/
void make_random_sarray_from_prob(prob,period,spikes,count,n,seed)
     float *prob;
     int period,***spikes,**count,n;
     int seed;
{
  int i,j,k;
  int **s,*cnt,*tdata;

  if (seed > 0)
    seed = -seed;
  spike_util_ran2(&seed);

  cnt = (int *)myalloc(n*sizeof(int));
  s = (int **)myalloc(n*sizeof(int *));
  tdata = (int *)myalloc(period*sizeof(int));
  for(i=0;i<n;i++){
    k = 0;
    for(j=0;j<period;j++)
      if (prob[j] > spike_util_ran2(&seed)){
	tdata[k] = j;
	k += 1;
      }
    s[i] = copy_iarray(tdata,k);
    cnt[i] = k;
  }
  myfree(tdata);
  *spikes = s;
  *count = cnt;
}
/**************************************-**************************************/
/*                                                                           */
/*                       MAKE_RANDOM_SARRAY_FROM_SARRAY                      */
/*                                                                           */
/*  Return an "sarray" that has spike chosen according to the prob-          */
/*  ability of firing computed from the PST of the given sarray.             */
/*  Make "nrec" records.                                                     */
/*                                                                           */
/*****************************************************************************/
void make_random_sarray_from_sarray(s,cnt,n,period,spikes,count,nrec,seed)
     int **s,*cnt,n,period,***spikes,**count,nrec;
     int seed;
{
  int i;
  int *pst;
  float *prob;

  pst = simple_pst_sarray(s,cnt,n,0,period);
  prob = (float *)myalloc(period*sizeof(float));
  for(i=0;i<period;i++)
    prob[i] = (float)pst[i]/(float)n;

  make_random_sarray_from_prob(prob,period,spikes,count,nrec,seed);
  myfree(pst); myfree(prob);
}
/**************************************-**************************************/
/*                                                                           */
/*                     MAKE_POISSON_FLOAT_SPIKES_FROM_PROB                   */
/*                                                                           */
/*  Use a time warping algorithm to make Poisson spikes.                     */
/*                                                                           */
/*  *** NOTE: This routine can return duplicate spike times.                 */
/*  *** NOTE: 'pseed' must point to NEGATIVE value to INITIALIZE.            */
/*                                                                           */
/*****************************************************************************/
void make_poisson_float_spikes_from_prob(prob,n,sampling,pseed,rs,rn)
     float *prob;    // Probability of a spike in each time bin
     int n;          // Number of time bins
     float sampling; // Time bins per second
     int *pseed;     // Pointer to randomization seed (value gets modified)
     float **rs;     // Returned array of spike times, floating point
     int *rn;        // Length of returned spike time array
{
  int i,j,k;
  int fn,*warp,t;
  float mu,lambda,*fs,*fint;

  // Compute the integral of the firing rate, divide by total.
  fint = (float *)myalloc(n*sizeof(float));
  fint[0] = prob[0];
  for(i=1;i<n;i++)
    fint[i] = fint[i-1] + prob[i];
  multiply_farray(fint,n,1.0/fint[n-1]);
  //write_farray_plot("zzz.fint.pl",fint,n);
  //append_farray_plot("zzz.fint.pl","fint",fint,n,1);

  j = 0;
  warp = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++){
    // WYETH line below was changed 2011 May, to avoid spikes at t=0
    //   when the prob was 0 for an epoch at the beginning.  This got
    //   rid of warp[0] = 0, in such a case.
    //while (fint[j] < ((float)i/(float)n))
    while (fint[j] <= ((float)i/(float)n))
      j += 1;
    warp[i] = j;
  }
  myfree(fint);
  //write_iarray_plot("zzz.warp.pl",warp,n);
  //append_iarray_plot("zzz.warp.pl","warp",warp,n,1);


  // Get homogeneous Poisson spikes with mean rate matched to 'prob'
  mu = mean_farray(prob,n); // Mean firing probability.
  lambda = mu * sampling;
  k = 1;
  make_poisson_float_spikes(&fs,&fn,0,n,sampling,lambda,0.0,0.0,k,pseed);

  // Warp to inhomogenous
  for(j=0;j<fn;j++){
    t = my_rint(fs[j]);
    if (t == n)  // Wyeth - Round can give values = 'n', check added Nov2009
      t = n-1;
    else if (t > n)
      exit_error("MAKE_POISSON_FLOAT_SPIKES_FROM_PROB","t too large");
    fs[j] = (float)warp[t];
  }
  myfree(warp);
  
  *rs = fs;
  *rn = fn;
}
/**************************************-**************************************/
/*                                                                           */
/*                      MAKE_POISSON_F_SPIKES_FROM_PR_REFR                   */
/*                                                                           */
/*  Like 'make_poisson_float_spikes_from_prob', but w/ refractory period.    */
/*                                                                           */
/*  *** NOTE: Using the refractory period with the time warping algorithm    */
/*            will have complex effects.  The time varying mean rate will    */
/*            not be preserved.  The refractory period will be shorter       */
/*            when the rate is higher.                                       */
/*                                                                           */
/*  Created for use with the 'LinearCommon.moo' model.                       */
/*                                                                           */
/*****************************************************************************/
void make_poisson_f_spikes_from_pr_refr(prob,n,sampling,pseed,refmu,refsig,
					rs,rn)
     float *prob;    // Probability of a spike in each time bin
     int n;          // Number of time bins
     float sampling; // Time bins per second
     int *pseed;     // Pointer to randomization seed (value gets modified)
     float refmu;    // Refractory period mu (s)
     float refsig;   // Refractory period sigma (s)
     float **rs;     // Returned array of spike times, floating point
     int *rn;        // Length of returned spike time array
{
  int i,j,k;
  int fn,*warp,t;
  float mu,lambda,*fs,*fint;

  // Compute the integral of the firing rate, divide by total.
  fint = (float *)myalloc(n*sizeof(float));
  fint[0] = prob[0];
  for(i=1;i<n;i++)
    fint[i] = fint[i-1] + prob[i];
  multiply_farray(fint,n,1.0/fint[n-1]);
  //write_farray_plot("zzz.fint.pl",fint,n);
  //append_farray_plot("zzz.fint.pl","fint",fint,n,1);

  j = 0;
  warp = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++){
    while (fint[j] < ((float)i/(float)n))
      j += 1;
    warp[i] = j;
  }
  myfree(fint);
  //write_iarray_plot("zzz.warp.pl",warp,n);
  //append_iarray_plot("zzz.warp.pl","warp",warp,n,1);


  // Get homogeneous Poisson spikes with mean rate matched to 'prob'
  mu = mean_farray(prob,n); // Mean firing probability.
  lambda = mu * sampling;
  k = 1;
  make_poisson_float_spikes(&fs,&fn,0,n,sampling,lambda,refmu,refsig,k,pseed);

  // Warp to inhomogenous
  for(j=0;j<fn;j++){
    t = my_rint(fs[j]);
    if (t == n)  // Wyeth - Round can give values = 'n', check added Nov2009
      t = n-1;
    else if (t > n)
      exit_error("MAKE_POISSON_FLOAT_SPIKES_FROM_PROB","t too large");
    fs[j] = (float)warp[t];
  }
  myfree(warp);
  
  *rs = fs;
  *rn = fn;
}
/**************************************-**************************************/
/*                                                                           */
/*                       MAKE_POISSON_SPIKES_FROM_PROB                       */
/*                                                                           */
/*  Use a time warping algorithm to make Poisson spikes.                     */
/*                                                                           */
/*  *** NOTE: This routine can return duplicate spike times.                 */
/*                                                                           */
/*****************************************************************************/
void make_poisson_spikes_from_prob(prob,n,sampling,pseed,rs,rn)
     float *prob;    // Probability of a spike in each time bin
     int n;          // Number of time bins
     float sampling; // Time bins per second
     int *pseed;     // Pointer to randomization seed (value gets modified)
     int **rs,*rn;   // Returned array of spike times and array length
{
  int *s,fn;
  float *fs;

  make_poisson_float_spikes_from_prob(prob,n,sampling,pseed,&fs,&fn);
  s = f2iarray(fs,fn);
  myfree(fs);
  
  *rs = s;
  *rn = fn;
}
/**************************************-**************************************/
/*                                                                           */
/*                      MAKE_POISSON_SARRAY_FROM_PROB                        */
/*                                                                           */
/*  Use a time warping algorithm to make Poisson spikes.                     */
/*                                                                           */
/*  *** NOTE: This routine can return duplicate spike times.                 */
/*                                                                           */
/*****************************************************************************/
void make_poisson_sarray_from_prob(prob,period,spikes,count,n,seed)
     float *prob;
     int period,***spikes,**count,n;
     int seed;
{
  int i,j,k;
  int **s,*cnt,*warp,t;
  float mu,lambda,sampling,*fst,*fint;

  /*** Compute the integral of the firing rate, divide by total. ***/
  fint = (float *)myalloc(period*sizeof(float));
  fint[0] = prob[0];
  for(i=1;i<period;i++)
    fint[i] = fint[i-1] + prob[i];
  multiply_farray(fint,period,1.0/fint[period-1]);
  /*write_farray_plot("zzz.fint.pl",fint,period);*/

  j = 0;
  warp = (int *)myalloc(period*sizeof(int));
  for(i=0;i<period;i++){
    while (fint[j] < ((float)i/(float)period))
      j += 1;
    warp[i] = j;
  }
  myfree(fint);
  /*write_iarray_plot("zzz.warp.pl",warp,period);*/

  if (seed > 0)
    seed = -seed;
  spike_util_ran2(&seed);

  cnt = (int *)myalloc(n*sizeof(int));
  s = (int **)myalloc(n*sizeof(int *));

  k = 1;
  sampling = 1000.0;
  mu = mean_farray(prob,period); // Mean firing probability
  lambda = mu * sampling;
  for(i=0;i<n;i++){
    /*** Get homogeneous Poisson spikes. ***/
    make_poisson_float_spikes(&fst,&cnt[i],0,period,sampling,lambda,0.0,
			      0.0,k,&seed);
    /*** Warp to inhomogenous. ***/
    for(j=0;j<cnt[i];j++){
      t = my_rint(fst[j]);
      fst[j] = (float)warp[t];
    }
    s[i] = f2iarray(fst,cnt[i]);
    myfree(fst);
  }
  myfree(warp);

  *spikes = s;
  *count = cnt;
}
/**************************************-**************************************/
/*                                                                           */
/*                    MAKE_INTEGRATE_FIRE_SPIKES_FROM_FARRAY                 */
/*                                                                           */
/*  Use integrate and fire model to generate spike times.                    */
/*                                                                           */
/*  vrest - zero.                                                            */
/*  tau - number of steps to decay to 1/e = 0.36788                          */
/*  vthresh - threshold value for spike generation                           */
/*  vthrel - threshold elevation value                                       */
/*  vreset - value to                                                        */
/*  nmu - noise mean (Gaussian noise)                                        */
/*  nsd - noise standard deviation                                           */
/*  pseed - pointer to randomization seed                                    */
/*  rs - returned spike times                                                */
/*  rn - returned number of spikes                                           */
/*  saveflag - 0-nosave, 1-overwrite, 2-append                               */
/*  outfile - output file name                                               */
/*                                                                           */
/*****************************************************************************/
void make_integrate_fire_spikes_from_farray(d,period,tau,vthresh,vthrel,eltau,
					    vreset,nmu,nsd,pseed,saveflag,
					    outfile,rs,rn)
     float *d; /* Input drive */
     int period;
     float tau,vthresh,vthrel,eltau,vreset,nmu,nsd;
     int saveflag;
     char outfile[];
     int *pseed,**rs,*rn;
{
  int i;
  int *t,*s,n;
  float v,vrest,dtau,deltau,thr;
  float *sav_v,*sav_thr;

  if (vthresh <= 0.0)
    exit_error("MAKE_INTEGRATE_FIRE_SPIKES_FROM_FARRAY","Vthresh <= 0.0");
  if (vreset >= vthresh)
    exit_error("MAKE_INTEGRATE_FIRE_SPIKES_FROM_FARRAY","Vreset >= Vthresh");
  if (vthrel < vthresh)
    exit_error("MAKE_INTEGRATE_FIRE_SPIKES_FROM_FARRAY","Vthrel < Vthresh");

  if (saveflag > 0){
    sav_v = (float *)myalloc(period*sizeof(float));
    sav_thr = (float *)myalloc(period*sizeof(float));
  }else
    sav_v = sav_thr = NULL; /* Avoid -Wall warning */

  if (tau <= 0.0)
    dtau = 1.0;
  else
    dtau = 1.0 - exp(-1.0/tau);
  
  if (eltau <= 0.0)
    deltau = 1.0;
  else
    deltau = 1.0 - exp(-1.0/eltau);

  /*printf("deltau = %f   dtau = %f\n",deltau,dtau);*/

  if (*pseed > 0)
    *pseed *= -1;

  t = (int *)myalloc(period*sizeof(int));

  vrest = 0.0;
  v = vrest;
  n = 0;
  thr = vthresh;
  for(i=0;i<period;i++){
    v += d[i]; /* Add input */
    if (nsd > 0.0) /* Add noise */
      v += su_gasdev(pseed)*nsd + nmu;

    if (saveflag > 0){
      sav_v[i] = v;
      sav_thr[i] = thr;
    }

    if (v > thr){
      t[n] = i;
      n += 1;
      v = vreset;
      thr = vthrel;
    }else{
      v -= dtau*(v-vrest);
      thr -= deltau*(thr-vthresh);
    }
  }
  /*printf("    %d spikes\n",n);*/
  
  s = copy_iarray(t,n);
  myfree(t);

  if (saveflag > 0){
    if (saveflag == 1)
      write_farray_plot(outfile,sav_v,period);
    else
      append_farray_plot(outfile,"v",sav_v,period,1);
    append_farray_plot(outfile,"vthr",sav_thr,period,1);
    myfree(sav_v); myfree(sav_thr);
  }
  *rs = s; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                              SPIKEU_POISSON_G                             */
/*                                                                           */
/*  Return a conductance trace that is generted from the given Poisson rate  */
/*  and PSG mask.                                                            */
/*                                                                           */
/*  The spike train is returned if 'rs' is not null.                         */
/*                                                                           */
/*****************************************************************************/
float *spikeu_poisson_g(tn,sampling,rate,seed,mask,maskn,rs,rsn)
     int tn;          // Duration (sampling units)
     float sampling;  // Sampling (samples/s), e.g., 1000 for ms
     float rate;      // Spike rate (spikes/s)
     int seed;        // Randomization seed
     float *mask;     // post-syn conductance waveform [maskn]
     int maskn;       // length of mask
     float **rs;      // Spike train [rsn]
     int *rsn;        // number of spikes
{
  int n,*s;
  float *fs,*xs,*gt;

  if (rate > 0.0){

    if (seed > 0)
      seed = -seed;
    make_poisson_float_spikes(&fs,&n,0,tn,sampling,rate,0.0,0.0,1,&seed);
    s = f2iarray(fs,n);  // Convert to int
    xs = expand_spike_array(s,n,0,tn);  // Expand
    gt = convolve_with_mask_causal(xs,tn,mask,maskn);  // Convolve

    myfree(s);

    if (rs != NULL){  // Return the spikes
      *rs = fs;
      *rsn = n;
    }else{            // Discard the spikes
      myfree(fs);
      //*rs = NULL;
      //*rsn = 0;
    }
    myfree(xs);
  }else{
    gt = get_zero_farray(tn);
    if (rs != NULL){
      *rs = NULL;
      *rsn = 0;
    }
  }

  return gt;
}
/**************************************-**************************************/
/*                                                                           */
/*                            SPIKEU_POISSON_G_PROB                          */
/*                                                                           */
/*  Return a conductance trace that is generted from the given time-varying  */
/*  rate for an inhomogeneous Poisson process, and give the PSG mask.        */
/*                                                                           */
/*  The spike train is returned if 'rs' is not null.                         */
/*                                                                           */
/*****************************************************************************/
float *spikeu_poisson_g_prob(prob,tn,sampling,seed,mask,maskn,rs,rsn)
     float *prob;     // Probability vs. time [tn]
     int tn;          // Duration (sampling units)
     float sampling;  // Sampling (samples/s), e.g., 1000 for ms
     int seed;        // Randomization seed
     float *mask;     // post-syn conductance waveform [maskn]
     int maskn;       // length of mask
     float **rs;      // Spike train [rsn]
     int *rsn;        // number of spikes
{
  int n,*s;
  float *fs,*xs,*gt;

  if (tn <= 0)
    exit_error("SPIKEU_POISSON_G_PROB","Data length error");
  if (maskn <= 0)
    exit_error("SPIKEU_POISSON_G_PROB","Mask length error");

  if (seed > 0)
    seed = -seed;

  make_poisson_float_spikes_from_prob(prob,tn,sampling,&seed,&fs,&n);

  s = f2iarray(fs,n);  // Convert to int
  xs = expand_spike_array(s,n,0,tn);  // Expand
  gt = convolve_with_mask_causal(xs,tn,mask,maskn);  // Convolve

  //append_farray_plot("zz.gx_pre.pl","su_mask",mask,maskn,1);
  //append_farray_plot("zz.gx_pre.pl","su_expanded_spikes",xs,tn,1);
  //append_farray_plot("zz.gx_pre.pl","su_gex",gt,tn,1);

  if (rs != NULL){  // Return the spikes
    *rs = fs;
    *rsn = n;
  }else{            // Discard the spikes
    myfree(fs);
  }
  myfree(xs);
  myfree(s);

  return gt;
}
/**************************************-**************************************/
/*                                                                           */
/*                                ROC_SARRAYS                                */
/*                                                                           */
/*  Typically, if first sarray has more spikes, the ROC value will           */
/*  be > 0.5.                                                                */
/*                                                                           */
/*****************************************************************************/
float roc_sarrays(s1,cnt1,n1,s2,cnt2,n2,start,period)
     int **s1,*cnt1,n1,**s2,*cnt2,n2,start,period;
{
  int i;
  float *c1,*c2,roc;

  c1 = (float *)myalloc(n1*sizeof(float));
  c2 = (float *)myalloc(n2*sizeof(float));

  for(i=0;i<n1;i++)
    c1[i] = (float)count_spikes(s1[i],cnt1[i],start,period);
  for(i=0;i<n2;i++)
    c2[i] = (float)count_spikes(s2[i],cnt2[i],start,period);

  roc = roc_farrays(c1,c2,n1,n2);
  myfree(c1); myfree(c2);
  return roc;
}
/**************************************-**************************************/
/*                                                                           */
/*                            PICK_CROSS_SPIKES                              */
/*                                                                           */
/*  Select spikes from "chan1" which match a condition on the spikes         */
/*  from "chan2".                                                            */
/*                                                                           */
/*****************************************************************************/
void pick_cross_spikes(s1,cnt1,s2,cnt2,n,start,period,lag,width,rs,rcnt)
     int **s1,*cnt1,**s2,*cnt2,n,start,period;
     int lag,width; /* must be a spike within "width" starting at "lag" */
     int ***rs,**rcnt; /* returned sarray */
{
  int i,j;
  int **s,*cnt,*ts,tcnt,t;

  s = (int **)myalloc(n*sizeof(int *));
  cnt = (int *)myalloc(n*sizeof(int));
  ts = (int *)myalloc(period*sizeof(int));

  for(i=0;i<n;i++){
    cnt[i] = 0;
    for(j=0;j<cnt1[i];j++){ /* for each spike in ith train */
      t = s1[i][j];
      tcnt = count_spikes(s2[i],cnt2[i],t+lag,width);
      if (tcnt > 0){
	ts[cnt[i]] = t;
	cnt[i] += 1;
      }
    }
    s[i] = extract_spikes(ts,cnt[i],start,period,&t,0);
    if (t!=cnt[i])
      exit_error("PICK_CROSS_SPIKES","Counts not equal");
  }
      
  myfree(ts);
  *rs = s;
  *rcnt = cnt;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_STATIONARITY_SARRAY                          */
/*                                                                           */
/*  Return the minimum and maximum values of the smoothed spike              */
/*  count as a function of index during the sarray (thus, sequency           */
/*  during the experiment.                                                   */
/*                                                                           */
/*****************************************************************************/
void get_stationarity_sarray(s,cnt,n,min,max)
     int **s,*cnt,n;
     float *min,*max;
{
  float *fcnt,*smooth,sigma;

  if (n < 2)
    exit_error("GET_STATIONARITY_SARRAY","Too few trials");
  sigma = 10.0; /* standard deviation of smoothing (in trials) */
  fcnt = i2farray(cnt,n);
  smooth = smooth_with_gaussian(fcnt,n,sigma,0.01);
  *min = min_of_farray(smooth,n);
  *max = max_of_farray(smooth,n);
  myfree(fcnt); myfree(smooth);
}
/**************************************-**************************************/
/*                                                                           */
/*                           SPIKE_WORD_RECODE_SARRAY                        */
/*                                                                           */
/*  Recode the sarray as an expanded data array that gives the spike word    */
/*  index at each time.                                                      */
/*                                                                           */
/*****************************************************************************/
void spike_word_recode_sarray(s,cnt,n,start,period,nlet,dt,binflag,shiftflag,
			      rswi,rwpt)
     int **s,*cnt,n,start,period;
     int nlet;  /* Number of letters in the word */
     int dt;    /* Number of sampling units per letter */
     int binflag; /* Limit to binary words */
     int shiftflag; /* Move spikes to neighboring bin to maintain binary */
     int ***rswi;   /* Word indices at each time */
     int *rwpt;     /* Words per trial */
{
  int i,j,k,l;
  int *x,wpt,wlen,*wd,base0,base,wi;
  int **swi;

  if (binflag != 1)
    exit_error("SPIKE_WORD_RECODE_SARRAY","Must set binflag to 1");

  /*** Make storage for the probability distribution ***/
  base0 = 1;
  for(i=0;i<nlet;i++)
    base0 *= 2;

  wlen = nlet * dt;         /* Word length in spike time units */
  wpt = period - wlen + 1;  /* Words per trial */
  swi = get_zero_2d_iarray(n,wpt);
  /*printf("  base0= %d  wlen= %d  wpt= %d  nlet=%d  dt=%d\n",base0,wlen,wpt,
    nlet,dt);*/
  
  wd = (int *)myalloc(nlet*sizeof(int));
  for(i=0;i<n;i++){ /* For each trial */
    x = expand_spike_int_array(s[i],cnt[i],start,period);
    for(j=0;j<wpt;j++){ /* For each starting time */
      for(k=0;k<nlet;k++){ /* For each letter in the word */
	wd[k] = 0;
	for(l=0;l<dt;l++){
	  wd[k] += x[j+k*dt+l];
	}
      }
      /* Now the word 'wd[nlet]' tells the # spikes for each letter */
      if (binflag == 1){
	if (shiftflag == 0){
	  for(k=0;k<nlet;k++)
	    if (wd[k] > 1)
	      wd[k] = 1;
	}else{
	  exit_error("SPIKE_WORD_RECODE_SARRAY","Must set shiftflag to 0");
	}
      }
      wi = 0; /* word index */
      base = base0/2;
      for(k=0;k<nlet;k++){
	/*printf(" %d",wd[k]);*/
	wi += wd[k] * base;
	base /= 2;
      }
      /*printf("  %d\n",wi);*/
      if (wi >= base0){
	exit_error(" **** HERE HERE","wi > base0");
      }
      swi[i][j] = wi;
    }
    myfree(x);
  }
  myfree(wd);
  *rswi = swi; *rwpt = wpt;
}
/**************************************-**************************************/
/*                                                                           */
/*                          SPIKE_WORD_DISTRIB_SARRAY                        */
/*                                                                           */
/*  Return the distribution of spike words.                                  */
/*                                                                           */
/*****************************************************************************/
void spike_word_distrib_sarray(s,cnt,n,start,period,nlet,dt,binflag,shiftflag,
			       rp,rpn)
     int **s,*cnt,n,start,period;
     int nlet;  /* Number of letters in the word */
     int dt;    /* Number of sampling units per letter */
     int binflag; /* Limit to binary words */
     int shiftflag; /* Move spikes to neighboring bin to maintain binary */
     float **rp;
     int *rpn;
{
  int i,j,k;
  int *x,xn,ntrunc,nshift,tot_trunc,tot_shift,np,base,wi,nw;
  float *p;

  if (binflag != 1)
    exit_error("SPIKE_WORD_DISTRIB_SARRAY","Must set binflag to 1");

  /*** Make storage for the probability distribution ***/
  np = 1;
  for(i=0;i<nlet;i++)
    np *= 2;
  p = get_zero_farray(np);

  tot_trunc = tot_shift = 0;
  for(i=0;i<n;i++){
     x = histogram_expand_sarray(s[i],cnt[i],start,period,dt,binflag,
       shiftflag,&xn,&ntrunc,&nshift);

    tot_trunc += ntrunc;
    tot_shift += nshift;
    /***  FOR CHECKING
      append_iarray_plot("zzz.pl","x",x,xn,1);
      {
      float *xf;
      xf = expand_spike_array(s[i],cnt[i],start,period);
      append_farray_plot("zzz.pl","fx",xf,period,1);
      }***/

    /*** Catalog each word in the string ***/
    nw = xn - (nlet - 1);
    for(j=0;j<nw;j++){ /* For each start index in the code array */
      wi = 0; /* word index */
      base = np/2;
      for(k=0;k<nlet;k++){
	/*printf("j+k=%d %d",j+k,x[j+k]);*/
	wi += x[j+k] * base;
	base /= 2;
      }
      /*printf("  %d\n",wi);*/
      if (wi >= np){
	exit_error(" **** HERE HERE","wi > np");
      }
      p[wi] += 1.0;
    }

    myfree(x);
  }

  if (binflag == 1){
    printf("      Total bins truncated:  %d\n",tot_trunc);
    if (shiftflag == 1)
      printf("      Total spikes shifted:  %d\n",tot_shift);
  }
  *rp = p; *rpn = np;
}
/**************************************-**************************************/
/*                                                                           */
/*                     STIM_CONDITION_WORD_DISTRIB_SARRAY                    */
/*                                                                           */
/*  Return the frequency of occurrence of spike words conditioned on the     */
/*  occurrence of stimulus patterns.                                         */
/*                                                                           */
/*  pr[ti][si][ri]                                                           */
/*    ti - time offset index [0..tn] for offsets from t0 to t0+tn-1          */
/*    si - stimulus pattern index  sord^swn                                  */
/*    ri - response word index                                               */
/*                                                                           */
/*  "start" and "period" can be used to select the segment of spike data     */
/*  from which spikes will be considered.  "doffset" is the start time of    */
/*  the stimulus relative to time zero for the spikes.                       */
/*                                                                           */
/*  Returns:                                                                 */
/*    pr[tn][sn][rn] - array of response counts                              */
/*    sn - number of unique stimulus words                                   */
/*    rn - number of unique response words                                   */
/*    (tn is not returned, because it is sent in)                            */
/*                                                                           */
/*****************************************************************************/
void stim_condition_word_distrib_sarray(s,cnt,n,start,period,shiftflag,binflag,
					stim,stimn,stimt0,stimdt,
					t0,tn,sord,swn,rord,rwn,rdt,
					r_pr,r_sn,r_rn)
     int **s,*cnt,n,start,period;  /* Spike response */
     int shiftflag,binflag;
     int **stim,*stimn;  /* Stimulus, encoded as word ID */
     int stimt0; /* Time of stimulus onset relative to spike time 0 */
     float stimdt; /* Exact time of one stimulus frame (s sampling units) */
     int t0,tn;  /* Start and duration of temporal offsets to consider */
     int sord;   /* Stimulus order (2-binary, 3-ternary, etc) */
     int swn;    /* Stimulus word length */
     int rord;   /* Stimulus order (2-binary) */
     int rwn;    /* Response word length */
     int rdt;    /* Response letter dt (temporal discretization) */
     float ****r_pr;  /* Returned array of counts */
     int *r_sn,*r_rn; /* 2nd and 3rd dimension lengths for 'pr' */
{
  int i,j,k;
  int si,ri,sn,rn,wpt,tw0,**rwi,rwin;
  float ***pr;

  printf("  STIM_CONDITION_WORD_DISTRIB_SARRAY\n");
  printf("    %d trials, start %d, period %d\n",n,start,period);
  printf("    Stimulus onset relative to spike t0:  %d\n",stimt0);
  printf("    Stimulus frame exact dt:  %f\n",stimdt);

  spike_word_recode_sarray(s,cnt,n,start,period,rwn,rdt,binflag,shiftflag,
			   &rwi,&rwin);
  /*append_iarray_plot("zz.rwi.pl","rwi_0",rwi[0],rwin,1);*/

  if (rord != 2)
    exit_error("STIM_CONDITION_WORD_DISTRIB_SARRAY","Response must be binary");

  sn = 1; /*** Determine number of stimulus words ***/
  for(i=0;i<swn;i++)
    sn *= sord;

  rn = 1; /*** Determine number of response words ***/
  for(i=0;i<rwn;i++)
    rn *= rord;

  pr = get_zero_3d_farray(tn,sn,rn);
  printf("    Probability matrix is  %d t  x  %d stim  x  %d resp\n",tn,sn,rn);

  for(i=0;i<n;i++){ /*** For each trial ***/
    wpt = stimn[i];
    for(j=0;j<wpt;j++){ /*** For each stimulus word in the trial ***/
      si = stim[i][j]; /* Get index of this stimulus word */
      if ((si < 0)|| (si >= sn))
	exit_error("STIM_CONDITION_WORD_DISTRIB_SARRAY","si out of range");
      tw0 = stimt0 + my_rint((float)j*stimdt); /* Word start wrt spikes */
      tw0 += t0; /* Add information window lag time */
      for(k=0;k<tn;k++){ /*** For each temporal offset ***/
	if (tw0+k < rwin){ /* If full response word exists at this offset */
	  ri = rwi[i][tw0+k];
	  if ((ri < 0)|| (ri >= rn))
	    exit_error("STIM_CONDITION_WORD_DISTRIB_SARRAY","ri out of range");
	  pr[k][si][ri] += 1.0; /* Count this response word */
	}
      }
    }
  }
  free_2d_iarray(rwi,n);
  *r_pr = pr; *r_sn = sn; *r_rn = rn;
}
/**************************************-**************************************/
/*                                                                           */
/*                       SPIKE_UTIL_STA_R_V_STIM_SARRAY                      */
/*                                                                           */
/*  This method can be thought of as running along s(t) = STA(t) * STIM(t)   */
/*  and associating a mean response (in a window) with each value of s(t)    */
/*  at the time resolution of the signal (typically 1 ms).                   */
/*                                                                           */
/*****************************************************************************/
void spike_util_sta_r_v_stim_sarray(s,cnt,n,start,period,sampling,
				    stim,stimn,stimt0,stscale,
				    sta,stan,toff,stn,normstim,normr,
				    outfile,nbin,
				    mindata,dump_nbin,dump_stscale,write_n,
				    var_flag,sigsm,cond,
				    rxx,ryy,rnn)
     int **s,*cnt,n;     /* Spike array */
     int start,period;   /* Analysis epoch, in terms of spike times */
     float sampling;  /* */
     float **stim;  /* Stimulus [n][stimn[n]], must match spike sampling */
     int *stimn;    /* Length of stimulus */
     int stimt0;    /* Time of stimulus onset relative to spike time 0 */
     float stscale; /* Scale stimulus */
     float *sta;    /* STA [stan] */
     int stan;      /* length of STA */
     int toff;      /* Time offset of response averaging window */
     int stn;       /* Time duration of response averaging window */
     int normstim;  /* Normalize stimulus */
     int normr;     /* Normalize response (divide by mean) */
     char outfile[]; /* Output file */
     int nbin;       /* Number of bins for x-axis */
     int mindata;    /* Minimum number of points per bin */
     int dump_nbin;  /* Dump stim/resp for this many upper bins */
     float dump_stscale;  /* Dump stim/resp for this 'stscale' value */
     int write_n;    /* If 1, write the number of samples per bin */
     int var_flag;   /* Compute variance, or dump all values */
                     /*   1 - write SD and number of values */
                     /*   4 - dump all values */
     float sigsm;    /* Sigma for Gaussian response window */
     struct spike_util_cond_struct *cond;    /* Conditional RVS */
     float **rxx,**ryy;  /* Return x and y coords */
     int *rnn;           /* Return number of values */
{
  int i,j,k,t;
  int sn,*rvsn,t1,st0,maxi,nmax,cflag,tc0,ccn0,ccn1,t0;
  float min,max,x,y,sr,srsd,mean,sd,*sx,*smooth;
  float *conv,*star,*rvs,**rvs2d,**rvs2dx,*rx,*xx,*yy,*yn,*zstim;
  char name[SLEN],tstr[SLEN3],tstr1[SLEN];
  int run_flag,ndump;

  ndump = 200;

  printf("    SPIKE_UTIL_STA_R_V_STIM_SARRAY\n");
  printf("      Stimulus scaling:  %f\n",stscale);
  printf("      STA is %d long\n",stan);
  printf("      Stim is %d long\n",stimn[0]);
  printf("      %d trials\n",n);
  printf("      stimt0 = %d\n",stimt0);

  /*** Response vs. stimulus, binned ***/
  rvs = get_zero_farray(nbin);
  rvsn = get_zero_iarray(nbin);

  if (var_flag > 0){
    nmax = stimn[0] * n;  /* Maximum possible */
    rvs2d = get_2d_farray(nbin,nmax);
    rvs2dx = get_2d_farray(nbin,nmax);
  }

  star = copy_farray(sta,stan);
  reverse_farray(star,stan);
  if (normstim == 3){  /* Don't bother normalizing STA if norm'lzing later */
    norm_area_farray(star,stan,1.0);
  }

  /* Analysis start, converted to stim time (may be clipped) */
  t0 = start - stimt0;
  t1 = t0 + period;

  if (t0 < 0){
    t0 = 0;  /* Cannot start in stimulus before t=0, thus we clip here */
  }
  
  ccn0 = ccn1 = 0;
  cond->nsuc = cond->nfail = 0;
  for(i=0;i<n;i++){ /* For each trial */

    /*** Convolve stimulus with STA ***/
    sn = stimn[i];

    if (t1 > (sn-40)){ /* Stay away from end of stimulus */
      printf("  t1 = %d  sn-40 = %d\n",t1,sn-40);
      printf("  *** Decrease 'rvs_period', currently %d\n",period);
      exit_error("SPIKE_UTIL_STA_R_V_STIM_SARRAY","Too close to stim end");
    }

    conv = convolve_with_mask_causal(stim[i],sn,star,stan);

    /***  FOR DUMPING EXAMPLES 
    {
      char ttstr[SLEN];

      sprintf(ttstr,"conv_%d",i);
      append_farray_plot("zzz.dump.pl",ttstr,conv,sn,1);
      sprintf(ttstr,"stim_%d",i);
      append_farray_plot("zzz.dump.pl",ttstr,stim[i],sn,1);
      sprintf(ttstr,"sta_%d",i);
      append_farray_plot("zzz.dump.pl",ttstr,star,stan,1);

      sprintf(ttstr,"spk_%d",i);
      append_xplot_spike_sampling("zzz.dump.pl",ttstr,s[i],cnt[i],0,period,
				  0,sampling,0);
				  }***/

    if (normstim == 1){
      printf("    z-scoring stimulus %d\n",i);
      z_score_farray(conv,sn);
      if (cond->ctype != NULL){
	if (strcmp(cond->data,"stim")==0){
	  zstim = copy_farray(stim[i],sn);
	  printf(" WYETH HERE STIM STIM STIM STIM\n");
	  z_score_farray(zstim,sn);
	}
      }
    }else if (normstim == 2){
      z_score_farray(conv,sn);
      multiply_farray(conv,sn,stscale);
      mean_sdev_farray(conv,sn,&mean,&sd);
      printf("WYETH HERE Normstim = 2 ======>  mean = %f  SD = %f\n",mean,sd);
      if (cond->ctype != NULL){
	if (strcmp(cond->data,"stim")==0){
	  zstim = copy_farray(stim[i],sn);
	  z_score_farray(zstim,sn);
	}
      }
    }else
      multiply_farray(conv,sn,stscale);

    if (sigsm > 0.0){
      sx = expand_spike_array(s[i],cnt[i],start,period);
      smooth = smooth_with_gaussian(sx,period,sigsm,0.01);
      myfree(sx);
      sx = smooth;
      multiply_farray(sx,period,sampling);
    }

    if (i==0){
      get_min_max_farray(conv,sn,&min,&max);

      /*** Write out the sdev (for Nicolas) ***/
      /*
	mean_sdev_farray(conv,sn,&mean,&sd);
	sprintf(tstr,"%s.sd",outfile);
	sprintf(tstr1,"%.1f %.4f\n",stscale,sd);
	append_string_to_file(tstr,tstr1);*/

      printf("      min max = %f %f\n",min,max);
      /* Expand min,max to account for other trials */
      min *= 1.1;
      max *= 1.1;
    }

    /*** Extract stim values and firing rates ***/
    run_flag = 0; /* Starting a run of high stimulus values */
    for(t=t0;t<t1;t++){  /* Time w/i stimulus */
      x = conv[t];
      st0 = t + stimt0 + toff;

      if (sigsm > 0.0)
	y = sx[st0];    /* Use smoothed (Gaussian) estimate */
      else
	y = mean_spike_rate(s[i],cnt[i],st0,stn,sampling); /* boxcar */

      k = (int)((x-min)/(max-min) * (float)nbin);

      /*** Conditional RVS ***/
      if (cond->ctype != NULL){
	tc0 = t + cond->win0;
	/* Condition time window must be valid */
	if ((tc0 >= 0) && (tc0+cond->winn <= sn)){
	  if (strcmp(cond->data,"stim")==0){
	    cflag = varcond_check_condition(cond,tc0,zstim,sn,outfile);
	  }else{
	    cflag = varcond_check_condition(cond,tc0,conv,sn,outfile);
	  }
	  if (cflag){
	    ccn1 += 1;
	    cond->nsuc += 1;

	    /*
	    printf("t=%d  (i=%d)\n",t,i);
	    if ((ccn1 > 100) && (i > 0)){
	      write_farray_plot("zz.conv.pl",conv,sn);
	      exit(0);
	      }*/

	  }else{
	    ccn0 += 1;
	    cond->nfail += 1;
	  }
	}else
	  cflag = 0;
      }

      /*** Handle cases that are off by 1 bin ***/
      if (k == -1)
	k = 0;
      else if (k == nbin)
	k = nbin-1;

      if (cflag == 0){
	;  /* Condition not satisfied.  Do nothing. */
      }else if ((k < 0) || (k >= nbin)){
	printf("*** k = %d  (outside of bin range)\n",k);
	/*exit_error("SPIKE_UTIL_STA_R_V_STIM_SARRAY","Shouldn't happen");*/
      }else{
	rvs[k] += y;
	if (var_flag > 0){
	  rvs2d[k][rvsn[k]] = y;
	  rvs2dx[k][rvsn[k]] = x;
	}
	rvsn[k] += 1;
	
	/*** Dump out examples ***/
	if ((dump_nbin > 0) && (stscale == dump_stscale) &&
	    (k >= (nbin-dump_nbin))){
	  if (run_flag == 0){
	    printf("  Start of run, trial %d, k %d, time %d\n",i,k,st0);
	    run_flag = 1;     /* The start of a 'run' of high values */
	    
	    /* 'stim[i]' is 'sn' long */
	    append_farray_plot_offset("zzz.rvs.dump.pl","stim",
				      &(stim[i][t-160]),ndump,-160);
	    
	    append_xplot_spike_sampling("zzz.rvs.dump.pl","spikes",s[i],cnt[i],
					t-160,ndump+160,2,sampling,t);
	  }
	}else{
	  run_flag = 0;
	}
      }
    }
    if (sigsm > 0.0)
      myfree(sx);

    myfree(conv);
  }

  if (cond->ctype != NULL){
    printf("      Condition satisfied: %7d\n",ccn1);
    printf("      Condition failed:    %7d\n",ccn0);
    printf("        Success rate =  %.2f %%\n",100.0*
	   (float)ccn1/((float)(ccn1+ccn0)));
    cond->fsuc = (float)cond->nsuc/((float)(cond->nsuc+cond->nfail));
  }

  sr = mean_spike_rate_sarray(s,cnt,n,start,period,sampling,&srsd);

  /*** Compute the Resp vs. Stim curve ***/
  rx = (float *)myalloc(nbin*sizeof(float));
  for(i=0;i<nbin;i++){
    rx[i] = min + (max-min)/(float)nbin * (float)(i+0.5);
    if (rvsn[i] > 0){
      rvs[i] /= (float)rvsn[i];
      if (normr == 1)
	rvs[i] /= sr;
    }else
      rvs[i] = -1.0;

    /*** if (stscale == 32.0)
	 printf("bin[%d] = %d\n",i,rvsn[i]);***/
  }

  xx = (float *)myalloc(nbin*sizeof(float));
  yy = (float *)myalloc(nbin*sizeof(float));
  yn = (float *)myalloc(nbin*sizeof(float));
  k = 0;
  for(i=0;i<nbin;i++){
    if (rvsn[i] >= mindata){
      xx[k] = rx[i];
      if (rvs[i] > 0.0001)  /* Avoid zero values, for log plot */
	yy[k] = rvs[i];
      else
	yy[k] = 0.0001;
      yn[k] = (float)rvsn[i];
      k += 1;
      maxi = i;
    }
  }
  printf("      Bin %d is highest with > %d entries.\n",maxi,mindata);

  sprintf(name,"scale_%.2f",stscale);
  if (var_flag == 1){        /* Write SD and N */
    sprintf(tstr,"%s %d\n",name,nbin);
    append_string_to_file(outfile,tstr);
    for(i=0;i<nbin;i++){
      mean_sdev_farray(rvs2d[i],rvsn[i],&mean,&sd);
      sprintf(tstr,"%3d %5d %6.2f %8.4f\n",i,rvsn[i],mean,sd);
      append_string_to_file(outfile,tstr);
    }
  }else if ((var_flag == 4)||(var_flag == 5)){  /* Dump all values */
    FILE *fout,*fopen();
    
    if ((fout = fopen(outfile,"a"))==NULL)
      exit_error("SPIKE_UTIL_STA_R_V_STIM_SARRAY","Cannot open file");
    fprintf(fout,"/newplot\n");
    fprintf(fout,"/plotname %s\n",name);

    if (var_flag == 4){
      for(i=0;i<nbin;i++){
	for(j=0;j<rvsn[i];j++){
	  fprintf(fout,"%d %.2f\n",i,rvs2d[i][j]);
	}
      }
    }else{
      for(i=0;i<nbin;i++){
	for(j=0;j<rvsn[i];j++){
	  fprintf(fout,"%.2e %.2f\n",rvs2dx[i][j],rvs2d[i][j]);
	}
      }
    }
    fclose(fout);
  }else{
    /*** Default output, for xplot ***/
    append_farray_xy_plot(outfile,xx,yy,k,name);
    if (write_n == 1){
      sprintf(tstr,"N_%.2f",stscale);
      append_farray_xy_plot(outfile,xx,yn,k,tstr);
    }
  }

  if (var_flag > 0){
    free_2d_farray(rvs2d,nbin);
    free_2d_farray(rvs2dx,nbin);
  }

  myfree(rx); myfree(yn);

  /*myfree(xx); myfree(yy);*/

  *rxx = xx; *ryy = yy; *rnn = k;
}
/**************************************-**************************************/
/*                                                                           */
/*                   SPIKE_UTIL_STA_R_V_STIM_BAYES_SARRAY                    */
/*                                                                           */
/*  This method is based on Bayes rule and is designed to implement the      */
/*  procedure described by equation 7 of Brenner et al. 2000.                */
/*                                                                           */
/*****************************************************************************/
void spike_util_sta_r_v_stim_bayes_sarray(s,cnt,n,start,period,sampling,
					  stim,stimn,stimt0,stscale,
					  sta,stan,toff,stn,t0,normstim,normr,
					  outfile,nbin,mindata,
					  rxx,ryy,rnn)
     int **s,*cnt,n,start,period;  /* Spike array */
     float sampling;  /* */
     float **stim;  /* Stimulus [n][stimn[n]], must match spike sampling */
     int *stimn;    /* Length of stimulus */
     int stimt0;    /* Time of stimulus onset relative to spike time 0 */
     float stscale; /* Scale stimulus */
     float *sta;    /* STA [stan] */
     int stan;      /* length of STA */
     int toff;      /* Time offset of response averaging window */
     int stn;       /* Time duration of response averaging window */
     int t0;        /* Time in stim to begin processing */
     int normstim;  /* Normalize stimulus */
     int normr;     /* Normalize response (divide by mean) */
     char outfile[]; /* Output file */
     int nbin;       /* Number of bins for x-axis */
     int mindata;    /* Minimum number of points per bin */
     float **rxx,**ryy;  /* Return x and y coords */
     int *rnn;           /* Return number of values */
{
  int i,j,k,t;
  int sn,*rvsn,t1,psn,psgsn;
  float min,max,x,y,sr,srsd,mean,sdev;
  float *conv,*star,*ps,*psgs,*rvs,*rx,*xx,*yy;
  char name[SLEN];

  printf("    SPIKE_UTIL_STA_R_V_STIM_BAYES_SARRAY\n");
  printf("      Stimulus scaling:  %f\n",stscale);
  printf("      STA is %d long\n",stan);
  printf("      Stim is %d long\n",stimn[0]);
  printf("      %d trials\n",n);
  printf("      stimt0 = %d\n",stimt0);

  /*** Response vs. stimulus, binned ***/
  ps   = get_zero_farray(nbin);  /* Prob of stim value */
  psn  = 0;
  psgs = get_zero_farray(nbin);  /* Prob of stim value given a spike */
  psgsn = 0;
  rvs  = get_zero_farray(nbin);
  rvsn = get_zero_iarray(nbin);

  /*** Prepare the STA for convolution ***/
  star = copy_farray(sta,stan);
  reverse_farray(star,stan);
  norm_area_farray(star,stan,1.0);

  for(i=0;i<n;i++){ /* For each trial */
    /*** Convolve stimulus with STA ***/
    sn = stimn[i];
    t1 = sn - 40;  /* Stay away from end of stimulus */

    conv = convolve_with_mask_causal(stim[i],sn,star,stan);

    /*** Determine how to scale the stimulus values (horizontal axis) ***/
    if (normstim == 1)
      z_score_farray(conv,sn);
    else
      multiply_farray(conv,sn,stscale);

    /*** Find the range of stimulus values */
    if (i==0){
      get_min_max_farray(conv,sn,&min,&max);
      printf("  min max = %f %f\n",min,max);
      /* Expand min,max to account for other trials */
      min *= 1.1;
      max *= 1.1;
    }

    /*** Histogram stim values to form Pr(stim) ***/
    for(t=t0;t<t1;t++){
      x = conv[t];
      k = (int)((x-min)/(max-min) * (float)nbin);
      if ((k < 0) || (k >= nbin)){
	printf("k = %d\n",k);
	/*exit_error("SPIKE_UTIL_STA_R_V_STIM_BAYES_SARRAY","Unexpected");*/
      }else{
	ps[k] += 1.0;
	psn += 1;
      }
    }

    /*** Histogram stim values conditional on a spike:  Pr(stim|spike) ***/
    for(j=0;j<cnt[i];j++){
      t = s[i][j] - stimt0 - toff;
      if ((t >= t0) && (t < t1)){
	x = conv[t];
	k = (int)((x-min)/(max-min) * (float)nbin);
	if ((k < 0) || (k >= nbin)){
	  printf("k = %d\n",k);
	  /*exit_error("SPIKE_UTIL_STA_R_V_STIM_BAYES_SARRAY","Unexpected");*/
	}else{
	  psgs[k] += 1.0;
	  psgsn += 1;
	}
      }
    }
    myfree(conv);
  }

  /*** for normalizing plot in figure ***/
  mean_sdev_freq_hist_farray(ps,nbin,&mean,&sdev);
  printf("  mean = %f  sdev = %f\n",mean,sdev);

  multiply_farray(ps,nbin,1.0/(float)psn); /* Get probability */
  append_farray_plot(outfile,"ps",ps,nbin,1);

  multiply_farray(psgs,nbin,1.0/(float)psgsn); /* Get probability */
  append_farray_plot(outfile,"psgs",psgs,nbin,1);

  sr = mean_spike_rate_sarray(s,cnt,n,start,period,sampling,&srsd);

  rx = (float *)myalloc(nbin*sizeof(float));
  for(i=0;i<nbin;i++){
    if (ps[i] > 0.0){
      rvs[i] = sr * psgs[i] / ps[i];

      if (normr == 1)  /* Divide by mean firing rate */
	rvs[i] /= sr;

    }else
      rvs[i] = -1.0;

    rx[i] = min + (max-min)/(float)nbin * (float)(i+0.5);
  }

  xx = (float *)myalloc(nbin*sizeof(float));
  yy = (float *)myalloc(nbin*sizeof(float));
  k = 0;
  for(i=0;i<nbin;i++){
    if (ps[i] > 0.0){
      xx[k] = rx[i];
      yy[k] = rvs[i];
      k += 1;
    }
  }

  sprintf(name,"scale_%.2f",stscale);
  append_farray_xy_plot(outfile,xx,yy,k,name);

  myfree(rx); myfree(rvs); myfree(ps); myfree(psgs);

  /*myfree(xx); myfree(yy); */

  *rxx = xx; *ryy = yy; *rnn = k;
}
