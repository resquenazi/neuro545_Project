/*****************************************************************************/
/*                                                                           */
/*   sig_util.c                                                              */
/*   wyeth bair                                                              */
/*   caltech                                                                 */
/*   02/03/94                                                                */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "my_util.h"
#include "farray_util.h"
#include "fft_util.h"

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
/*                               SIG_UTIL_RAN2                               */
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
float sig_util_ran2(idum)
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
/*                              SIG_GASDEV  (p288-290)                       */
/*                                                                           */
/*    ********** Initialize "sig_util_ran2" before calling ***********       */
/*                                                                           */
/*  Returns a normally distributed deviate with zero mean and unit           */
/*  variance, using ran1(idum) as the source of uniform deviates.            */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float sig_gasdev(idum)
     int *idum;
{
  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;
  
  /*if (iset == 0) {*/
  if ((iset == 0) || (*idum <= 0)){  /* Wyeth, reset state if reseeded */
    do {
      v1=2.0*sig_util_ran2(idum)-1.0;
      v2=2.0*sig_util_ran2(idum)-1.0;
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
/*                             SIG_EXPDEV                                    */
/*                                                                           */
/*****************************************************************************/
float sig_expdev(lambda,idum)
     float lambda;
     int *idum;
{
  if (lambda == 0.0)
    return 0.0;
  else
    return(-log(sig_util_ran2(idum))/lambda);
}
/**************************************-**************************************/
/*                                                                           */
/*                         SIG_GAUSSIAN_FLOAT                                */
/*                                                                           */
/*   Return a Gaussian deviate, N(mu,sigma).                                 */
/*                                                                           */
/*****************************************************************************/
float sig_gaussian_float(mu,sigma,idum)
     float mu,sigma;
     int *idum;
{
  return(sig_gasdev(idum)*sigma + mu);
}
/**************************************-**************************************/
/*                                                                           */
/*                         MEAN_AUTOCORR_NOISE                               */
/*                                                                           */
/*  Compute the mean over all trials of the autocorrelation of the           */
/*  noise.  Noise is defined as the difference between a given               */
/*  trial and the mean.                                                      */
/*                                                                           */
/*****************************************************************************/
void mean_autocorr_noise(data,m,n,acmean,acsdev)
     float **data;
     int m,n;  /* m trials, n data points per trial */
     float **acmean,**acsdev;
{
  int i,j;
  float **ac,*mean,*sdev,*noise;

  printf("  MEAN_AUTOCORR_NOISE\n");

  mean_sdev_2d_farray(data,m,n,&mean,&sdev); /*** get signal mean and sdev ***/

  ac = (float **)myalloc(m*sizeof(float *));
  noise = (float *)myalloc(n*sizeof(float));
  for(i=0;i<m;i++){  /*** For each trial ***/
    printf("%d ",i); fflush(stdout);
    for(j=0;j<n;j++)   /*** compute noise ***/
      noise[j] = data[i][j] - mean[j];
    z_score_farray(noise,n);
    ac[i] = autocorr_farray(noise,n);
  }
  mean_sdev_2d_farray(ac,m,n,acmean,acsdev); /*** get noise mean and sdev ***/

  for(i=0;i<m;i++)
    myfree(ac[i]);
  myfree(ac);
  
  myfree(mean); myfree(sdev); myfree(noise);
}
/**************************************-**************************************/
/*                                                                           */
/*                            CROSSCORR_MEAN                                 */
/*                                                                           */
/*  Compute the cross-correlation of the mean of the signals.                */
/*   **************** NEVER USED!!! NEVER TESTED !!! **************          */
/*   **************** NEVER USED!!! NEVER TESTED !!! **************          */
/*   **************** NO ENTRY IN .h FILE          ! **************          */
/*                                                                           */
/*****************************************************************************/
void crosscorr_mean(data1,data2,m,n,ccorr)
     float **data1,**data2;
     int m,n;  /* m trials, n data points per trial */
     float **ccorr;
{
  int ncc;
  float *cc,*mean1,*mean2,*sdev1,*sdev2;

  printf("  CROSSCORR_MEAN\n");

  mean_sdev_2d_farray(data1,m,n,&mean1,&sdev1); /* get signal mean and sdev */
  mean_sdev_2d_farray(data2,m,n,&mean2,&sdev2);

  z_score_farray(mean1,n);
  z_score_farray(mean2,n);
  cc = correlate(mean1,mean2,n,n,&ncc); /*** in "fft_util.c" ***/

  *ccorr = cc;
  myfree(mean1); myfree(mean2); myfree(sdev1); myfree(sdev2);
}
/**************************************-**************************************/
/*                                                                           */
/*                         MEAN_CROSSCORR_NOISE                              */
/*                                                                           */
/*  Compute the mean over all trials of the cross-correlation of the         */
/*  noise.  Noise is defined as the difference between a given               */
/*  trial and the mean.                                                      */
/*                                                                           */
/*****************************************************************************/
void mean_crosscorr_noise(data1,data2,m,n,ccmean,ccsdev)
     float **data1,**data2;
     int m,n;  /* m trials, n data points per trial */
     float **ccmean,**ccsdev;
{
  int i,j,si;
  int ncc,shift_flag;
  float **cc,*mean1,*mean2,*sdev1,*sdev2,*noise1,*noise2;
  float mm1,mm2,ss;

  printf("  MEAN_CROSSCORR_NOISE\n");
  shift_flag = 1;

  mean_sdev_2d_farray(data1,m,n,&mean1,&sdev1); /* get signal mean and sdev */
  mean_sdev_2d_farray(data2,m,n,&mean2,&sdev2);

  mean_sdev_farray(mean1,n,&mm1,&ss);
  mean_sdev_farray(mean2,n,&mm2,&ss);

  cc = (float **)myalloc(m*sizeof(float *));
  noise1 = (float *)myalloc(n*sizeof(float));
  noise2 = (float *)myalloc(n*sizeof(float));
  for(i=0;i<m;i++){  /*** For each trial ***/
    si = (i+m/3) % m; /*** For shift predictor ***/
    if (shift_flag)
      printf("%d(%d) ",i,si);
    else
      printf("%d ",i);
    fflush(stdout);
    for(j=0;j<n;j++){   /*** compute noise ***/
      noise1[j] = data1[i][j] - mean1[j];
      /*noise1[j] = data1[i][j] - mm1;*/
      if (shift_flag)
	noise2[j] = data2[si][j] - mean2[j];
      /*noise2[j] = data2[si][j] - mm2;*/
      else
	noise2[j] = data2[i][j] - mean2[j];
      /*noise2[j] = data2[i][j] - mm2;*/
    }
    /**z_score_farray(noise1,n);
      z_score_farray(noise2,n);**/
    cc[i] = correlate(noise1,noise2,n,n,&ncc); /*** in "fft_util.c" ***/
  }
  printf("\n");
  mean_sdev_2d_farray(cc,m,ncc,ccmean,ccsdev); /* get c-corr mean and sdev */

  for(i=0;i<m;i++)
    myfree(cc[i]);
  myfree(cc);
  
  myfree(mean1); myfree(mean2); myfree(sdev1); myfree(sdev2);
  myfree(noise1); myfree(noise2);
}
/**************************************-**************************************/
/*                                                                           */
/*                            GAUSSIAN_WHITE_NOISE                           */
/*                                                                           */
/*  Return an array of Gaussian white noise with the specified               */
/*  mean "mu" and standard deviation "sigma".                                */
/*                                                                           */
/*****************************************************************************/
float *gaussian_white_noise(mu,sigma,n,idum)
     float mu,sigma;
     int n;
     int idum;
{
  int i;
  float *noise;

  if (idum > 0)
    idum = -idum;
  
  noise = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    noise[i] = sig_gaussian_float(mu,sigma,&idum);
  return noise;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GAUSSIAN_WHITE_NOISE_PSEED                        */
/*                                                                           */
/*  Return an array of Gaussian white noise with the specified               */
/*  mean "mu" and standard deviation "sigma".                                */
/*                                                                           */
/*****************************************************************************/
float *gaussian_white_noise_pseed(mu,sigma,n,pseed)
     float mu,sigma;
     int n;
     int *pseed;
{
  int i;
  float *noise;

  noise = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    noise[i] = sig_gaussian_float(mu,sigma,pseed);
  return noise;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GAUSSIAN_CORR_NOISE                            */
/*                                                                           */
/*  (1) Create Gaussian white noise.                                         */
/*  (2) Smooth with Gaussian SD="corr_sigma".                                */
/*  (3) Set mean and SD.                                                     */
/*                                                                           */
/*****************************************************************************/
float *gaussian_corr_noise(corr_sigma,mu,sigma,n,seed)
     float corr_sigma,mu,sigma;
     int n;
     int seed;
{
  int i;
  float *gwn,*data,mean,sdev,cf,x;

  if (sigma > 0.0){
    gwn = gaussian_white_noise(0.0,1.0,n,seed);
    data = smooth_with_gaussian(gwn,n,corr_sigma,0.01);
    myfree(gwn);
    mean_sdev_farray(data,n,&mean,&sdev);

    cf = sigma/sdev; // Subtract mean, scale SD to 'sigma', make mean 'mu'
    for(i=0;i<n;i++){
      x = data[i] - mean;
      data[i] = mu + x*cf;
    }
  }else
    data = get_const_farray(n,mu);
  
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GAUSSIAN_CORR_NOISE_PSEED                         */
/*                                                                           */
/*  (1) Create Gaussian white noise.                                         */
/*  (2) Smooth with Gaussian SD="corr_sigma".                                */
/*  (3) Set mean and SD.                                                     */
/*                                                                           */
/*****************************************************************************/
float *gaussian_corr_noise_pseed(corr_sigma,mu,sigma,n,pseed)
     float corr_sigma,mu,sigma;
     int n;
     int *pseed;  // Pointer to seed
{
  int i;
  float *gwn,*data,mean,sdev,cf,x;

  if (sigma > 0.0){
    gwn = gaussian_white_noise_pseed(0.0,1.0,n,pseed);
    data = smooth_with_gaussian(gwn,n,corr_sigma,0.01);
    myfree(gwn);
    mean_sdev_farray(data,n,&mean,&sdev);

    cf = sigma/sdev; // Subtract mean, scale SD to 'sigma', make mean 'mu'
    for(i=0;i<n;i++){
      x = data[i] - mean;
      data[i] = mu + x*cf;
    }
  }else
    data = get_const_farray(n,mu);
  
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                              LEAKY_INTEGRATE                              */
/*                                                                           */
/*  Return the array of data convolved with e^(-t/tau), where tau            */
/*  is in sampling units.                                                    */
/*                                                                           */
/*****************************************************************************/
float *leaky_integrate(data,n,tau,tolerance,nm)
     float *data;
     int n;
     float tau,tolerance;
     int *nm; /* return the length of the mask */
{
  int i;
  int nmask,nfull;
  float *mask,*result,mval;

  nmask = 0;
  mval = exp(-((float)nmask/tau));
  while(mval > tolerance){
    nmask += 1;
    mval = exp(-((float)nmask/tau));
  }
  nmask += 2;
  nfull = 2*nmask - 1;

  mask = get_zero_farray(nfull);
  for(i=0;i<nmask;i++)
    mask[nmask-1+i] = exp(-((float)i/tau));
  result = convolve_with_mask(data,n,mask,nfull);

  *nm = nmask;
  return result;
}
