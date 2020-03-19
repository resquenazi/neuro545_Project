/*****************************************************************************/
/*                                                                           */
/*  myrand_util.c                                                            */
/*  wyeth bair                                                               */
/*  caltech                                                                  */
/*  01/19/95                                                                 */
/*                                                                           */
/*  This file contains utilities involving random numbers.                   */
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

/*****************************************************************************/
/*                                                                           */
/*  From "util.c" in /usr/pkg/vnl/yarl/src/.                                 */
/*                                                                           */
/*****************************************************************************/
#define	LENG	17
static long ary[LENG] = {
  0x4B14EA50L,
  0x53C4A8E0L,
  0x67B1FA98L,
  0x9CFB3BB5L,
  0x18761AF1L,
  0x7970CD66L,
  0xDBAFE136L,
  0x3C31FC3EL,
  0x697B37DEL,
  0x07BC568BL,
  0xCAFD3967L,
  0xA8F48722L,
  0x4AB26824L,
  0xA479EE47L,
  0x5C7246E2L,
  0x954BF297L,
  0x20A713ADL,
};
static int i1 = 0;
static int i2 = 12;
/**************************************-**************************************/
/*                                                                           */
/*                                  YARLRAND                                 */
/*                                                                           */
/*  From "util.c" in /usr/pkg/vnl/yarl/src/.                                 */
/*                                                                           */
/*  Returns numbers between 0..32767 inclusive.                              */
/*                                                                           */
/*****************************************************************************/
int yarlrand()
{
  if(++i1 >= LENG)
    i1 = 0;
  if(++i2 >= LENG)
    i2 = 0;
  ary[i1] += ary[i2];
  return ((ary[i1] >> 15) & 0x7FFF);
}
/**************************************-**************************************/
/*                                                                           */
/*                                  YARLSRAND                                */
/*                                                                           */
/*  From "util.c" in /usr/pkg/vnl/yarl/src/.                                 */
/*                                                                           */
/*****************************************************************************/
void yarlsrand(seed)
     int seed;
{
  register int i;

  for(i=0;i<LENG;i++){
    ary[i] = 0x55555555L;
    if(seed & 1)
      ary[i] = 0xCCCCCCCCL;
    if(seed & 0x8000)
      ary[i] ^= 0xF0F0F0F0L;
    seed >>= 1;
  }
  i1 = 0;
  i2 = 12;
  for(i = 1; i < 32; i += i)
    do {
      (void)yarlrand();
      ary[i1] += ary[i1] >> i;
      ary[i2] += ary[i2] << i;
    } while (i1);
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
/*                              MYRAND_UTIL_RAN2                             */
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
float myrand_util_ran2(idum)
     int *idum; // Changed from long for consistency on SUNS, DEC ALPHA
{
  int j;
  int k; // Changed from long to int
  static int idum2=123456789; // Changed from long to int
  static int iy=0; // Changed from long to int
  static int iv[NTAB]; // Changed from long to int
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
/*                           MYRAND_GET_POISSON_COUNT                        */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Caller must be sure that 'myrand_util_ran2' was seeded. *************  */
/*                                                                           */
/*****************************************************************************/
int myrand_get_poisson_count(mu,rseed)
     float mu;
     int *rseed;
{  
  int n;
  float ftime;

  if (mu == 0.0){
    return 0;
  }else if (mu < 0.0)
    exit_error("MYRAND_GET_POISSON_COUNT","Mu is negative");

  ftime = 0.0;
  n = 0;
  while (ftime < mu){
    ftime += -log(myrand_util_ran2(rseed)); // Rate is one event per unit
    n += 1;
  }
  
  return n-1;
}
/**************************************-**************************************/
/*                                                                           */
/*                        MYRAND_GET_POISSON_COUNT_MEAN                      */
/*                                                                           */
/*  Return the mean of 'n' trials with rate 'mu'.                            */
/*                                                                           */
/*****************************************************************************/
float myrand_get_poisson_count_mean(float mu,      // Mean count
				    int n,         // Number of trials
				    int *rseed){   // Pointer to rand seed
  int i;
  float rv,ftime;

  rv = 0.0;

  if ((n > 0) && (mu > 0)){
    for(i=0;i<n;i++){
      rv += (float)myrand_get_poisson_count(mu,rseed);
    }
    rv /= (float)n;
  }

  return rv;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MYRAND_POISSON_COUNT                           */
/*                                                                           */
/*****************************************************************************/
int myrand_poisson_count(mu,seed)
     float mu;
     int seed;
{  
  int n;
  float ftime;

  if (seed > 0)
    seed = -seed;

  if (mu == 0.0){
    return 0;
  }else if (mu < 0.0)
    exit_error("MYRAND_POISSON_COUNT","Mu is negative");

  // generate first interval (sum of k intervals)
  ftime = 0.0;
  n = 0;
  while (ftime < mu){
    ftime += -log(myrand_util_ran2(&seed)); // Rate is one event per unit
    n += 1;
  }

  return n-1;
}
/**************************************-**************************************/
/*                                                                           */
/*                                   HPSORT                                  */
/*                                                                           */
/*  Sorts an array ra[1..n] into ascending numerical order using             */
/*  the Heapsort algorithm.  n is input; ra is replaced on output            */
/*  by its sorted rearrangement.                                             */
/*                                                                           */
/*  *** Modified by Wyeth to rearrange the array rx in parallel.             */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void hpsort(n,ra,rx)
     unsigned long n;
     float ra[];
     int *rx; /*wyeth*/
{
  unsigned long i,ir,j,l;
  float rra;
  int rrx; /*wyeth*/
  
  if (n < 2) return;
  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1) {
      rra=ra[--l];
      rrx=rx[l]; /*wyeth*/
    } else {
      rra=ra[ir];
      rrx=rx[ir]; /*wyeth*/
      ra[ir]=ra[1];
      rx[ir]=rx[1]; /*wyeth*/
      if (--ir == 1) {
	ra[1]=rra;
	rx[1]=rrx; /*wyeth*/
	break;
      }
    }
    i=l;
    j=l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
	ra[i]=ra[j];
	rx[i]=rx[j]; /*wyeth*/
	i=j;
	j <<= 1;
      } else j=ir+1;
    }
    ra[i]=rra;
    rx[i]=rrx; /*wyeth*/
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               SHUFFLE_IARRAY                              */
/*                                                                           */
/*  Shuffle an array of ints.                                                */
/*                                                                           */
/*****************************************************************************/
void shuffle_iarray(data,n,nshuffle,seed)
     int *data,n,nshuffle;
     int seed;
{
  int i,j;
  float *r; // random numbers to sort

  if (seed > 0) // Reseed random number generator
    seed = -seed;
  //myrand_util_ran2(&seed); // REMOVED, see old NYU version below

  r = (float *)myalloc(n*sizeof(float));
  for (i=0;i<nshuffle;i++){
    for (j=0;j<n;j++)
      r[j] = myrand_util_ran2(&seed);
    hpsort((unsigned long)n,r-1,data-1);
  }
  myfree(r);
}
/**************************************-**************************************/
/*                                                                           */
/*                             SHUFFLE_IARRAY_NYU                            */
/*                                                                           */
/*  Shuffle an array of ints.                                                */
/*                                                                           */
/*****************************************************************************/
void shuffle_iarray_nyu(data,n,nshuffle,seed)
     int *data,n,nshuffle;
     int seed;
{
  int i,j;
  float *r; /* random numbers to sort */

  /***  printf("  SHUFFLE_IARRAY\n");
    printf("    Shuffling %d times.\n",nshuffle);***/

  if (seed > 0) /*** Reseed random number generator. ***/
    seed = -seed;
  myrand_util_ran2(&seed); // THIS LINE IS IMPORTANT FOR 'tmon_grid'

  r = (float *)myalloc(n*sizeof(float));
  for (i=0;i<nshuffle;i++){
    for (j=0;j<n;j++)
      r[j] = myrand_util_ran2(&seed);
    hpsort((unsigned long)n,r-1,data-1);
  }
  myfree(r);
}
/**************************************-**************************************/
/*                                                                           */
/*                         SHUFFLE_IARRAY_RETURN_SEED                        */
/*                                                                           */
/*  Shuffle an array of ints.                                                */
/*                                                                           */
/*****************************************************************************/
void shuffle_iarray_return_seed(data,n,nshuffle,seed)
     int *data,n,nshuffle;
     int *seed;
{
  int i,j;
  float *r; // random numbers to sort
  
  if (*seed > 0) // Reseed random number generator
    *seed = -*seed;
  //myrand_util_ran2(seed);  WYETH REMOVED Apr 2007 [Keep this note]
  
  r = (float *)myalloc(n*sizeof(float));
  for (i=0;i<nshuffle;i++){
    for (j=0;j<n;j++)
      r[j] = myrand_util_ran2(seed);
    hpsort((unsigned long)n,r-1,data-1);
  }
  myfree(r);
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_SHUFFLE_INDEX                             */
/*                                                                           */
/*  Return a shuffled list of the numbers from 0 to n-1.                     */
/*                                                                           */
/*****************************************************************************/
int *get_shuffle_index(n,nshuffle,seed)
     int n,nshuffle;
     int seed;
{
  int i;
  int *s;

  s = (int *)myalloc(n*sizeof(int));
  for (i=0;i<n;i++)
    s[i] = i;
  shuffle_iarray(s,n,nshuffle,seed);
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_SHUFFLE_INDEX_RETURN_SEED                     */
/*                                                                           */
/*  Return a shuffled list of the numbers from 0 to n-1.                     */
/*                                                                           */
/*****************************************************************************/
int *get_shuffle_index_return_seed(n,nshuffle,rseed)
     int n,nshuffle;
     int *rseed;
{
  int i;
  int *s;

  s = (int *)myalloc(n*sizeof(int));
  for (i=0;i<n;i++)
    s[i] = i;
  shuffle_iarray_return_seed(s,n,nshuffle,rseed);
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 GET_SEEDS                                 */
/*                                                                           */
/*  Return an array of seeds from a single seed.                             */
/*                                                                           */
/*****************************************************************************/
int *get_seeds(seed,scale,n)
     int seed,scale,n;
{
  int i;
  int *s;

  if (seed == 0)
    exit_error("(myrand_util) GET_SEEDS","seed = 0");

  if (seed > 0) // Reseed random number generator
    seed = -seed;
  
  s = (int *)myalloc(n*sizeof(int));
  for (i=0;i<n;i++){

    s[i] = (int)(myrand_util_ran2(&seed) * (double)scale);  // WYETH NEW
    //s[i] = (int)(myrand_util_ran2(&seed) * (float)scale);
    // The line above gave different results on MacBook Pro and koch.
    while(s[i] == 0)
      s[i] = (int)(myrand_util_ran2(&seed) * (double)scale);
  }
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_SEEDS_2D                                */
/*                                                                           */
/*  Return a 2D array of seeds from a single seed.                           */
/*                                                                           */
/*****************************************************************************/
int **get_seeds_2d(seed,scale,xn,yn)
     int seed,scale;
     int xn,yn;
{
  int i,j;
  int **s;

  if (seed == 0)
    exit_error("(myrand_util) GET_SEEDS_2D","seed = 0");

  if (seed > 0) // Reseed random number generator
    seed = -seed;

  s = get_2d_iarray(xn,yn);

  for (i=0;i<xn;i++){
    for (j=0;j<yn;j++){

      s[i][j] = (int)(myrand_util_ran2(&seed) * (double)scale);
      while(s[i][j] == 0)
	s[i][j] = (int)(myrand_util_ran2(&seed) * (double)scale);
    }
  }
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_SEEDS_3D                                */
/*                                                                           */
/*  Return a 3D array of seeds from a single seed.                           */
/*                                                                           */
/*****************************************************************************/
int ***get_seeds_3d(seed,scale,xn,yn,zn)
     int seed,scale;
     int xn,yn,zn;
{
  int i,j,k;
  int ***s;

  if (seed == 0)
    exit_error("(myrand_util) GET_SEEDS_3D","seed = 0");

  if (seed > 0) // Reseed random number generator
    seed = -seed;

  s = get_3d_iarray(xn,yn,zn);

  for (i=0;i<xn;i++){
    for (j=0;j<yn;j++){
      for (k=0;k<zn;k++){

	s[i][j][k] = (int)(myrand_util_ran2(&seed) * (double)scale);
	while(s[i][j][k] == 0)
	  s[i][j][k] = (int)(myrand_util_ran2(&seed) * (double)scale);
      }
    }
  }
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                              GET_RANDOM_FLOATS                            */
/*                                                                           */
/*****************************************************************************/
float *get_random_floats(seed,n)
     int seed,n;
{
  int i;
  float *d;
  
  if (seed > 0) // Reseed random number generator
    seed = -seed;
  
  d = (float *)myalloc(n*sizeof(float));
  for (i=0;i<n;i++)
    d[i] = myrand_util_ran2(&seed);
  return d;
}
/**************************************-**************************************/
/*                                                                           */
/*                        MYRAND_GET_STD_UNIF_INT_SEQ                        */
/*                                                                           */
/*****************************************************************************/
int *myrand_get_std_unif_int_seq(n,seed,m)
     int n,seed;
     int m;  // Uniform values between 0 and m-1
{
  int i;
  int *s;
  float fm;
  
  fm = (float)m;
  
  if (seed == 0)
    exit_error("MYRAND_GET_STD_UNIF_INT_SEQ","seed = 0");
  else if (seed > 0)
    seed *= -1;
  
  s = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++)
    s[i] = (int)(fm * myrand_util_ran2(&seed));
  
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                        MYRAND_GET_STD_UNIF_FLOAT_SEQ                      */
/*                                                                           */
/*****************************************************************************/
float *myrand_get_std_unif_float_seq(n,seed,m)
     int n,seed;
     float m;  // Uniform values between 0 and m
{
  int i;
  float *s;

  if (seed == 0)
    exit_error("MYRAND_GET_STD_UNIF_FLOAT_SEQ","seed = 0");
  else if (seed > 0)
    seed *= -1;

  s = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    s[i] = m * myrand_util_ran2(&seed);

  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MYRAND_GET_STD_BIN_SEQ                          */
/*                                                                           */
/*****************************************************************************/
int *myrand_get_std_bin_seq(n,seed)
     int n,seed;
{
  int i;
  int *s;

  if (seed == 0)
    exit_error("MYRAND_GET_STD_BIN_SEQ","seed = 0");
  else if (seed > 0)
    seed *= -1;

  s = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++)
    if (myrand_util_ran2(&seed) > 0.5)
      s[i] = 1;
    else
      s[i] = 0;

  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MYRAND_GET_STD_TERN_SEQ                         */
/*                                                                           */
/*****************************************************************************/
int *myrand_get_std_tern_seq(n,seed)
     int n,seed;
{
  int i;
  int *s;
  float rv;

  if (seed == 0)
    exit_error("MYRAND_GET_STD_TERN_SEQ","seed = 0");
  else if (seed > 0)
    seed *= -1;

  s = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++){
    rv = myrand_util_ran2(&seed);
    if (rv < 0.33333)
      s[i] = -1;
    else if (rv < 0.66667)
      s[i] = 0;
    else
      s[i] = 1;
  }

  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MYRAND_GET_STD_QUAD_SEQ                         */
/*                                                                           */
/*  Returned values are:  0,1,2,3                                            */
/*                                                                           */
/*****************************************************************************/
int *myrand_get_std_quad_seq(n,seed)
     int n,seed;
{
  int i;
  int *s;
  float rv;

  if (seed == 0)
    exit_error("MYRAND_GET_STD_QUAD_SEQ","seed = 0");
  else if (seed > 0)
    seed *= -1;

  s = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++){
    rv = myrand_util_ran2(&seed);
    if (rv < 0.25)
      s[i] = 0;
    else if (rv < 0.5)
      s[i] = 1;
    else if (rv < 0.75)
      s[i] = 2;
    else
      s[i] = 3;
  }

  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MYRAND_GET_STD_MSEQ                           */
/*                                                                           */
/*****************************************************************************/
int *myrand_get_std_mseq(n,tap_seed,ord)
     int n;
     int tap_seed;  // Tap register seed value
     int ord;       // Order of m-sequence
{
  int i;
  int *s,nm;
  int ms_r,ms_b,ms_q,ms_t;

  nm = 1;
  for(i=0;i<ord;i++)
    nm *= 2;
  if (n > nm){
    printf("  *** m-seq order %d had %d values\n",ord,nm);
    printf("  *** %d values were requested\n",n);
    exit_error("MYRAND_GET_STD_MSEQ","requested too many values for order");
  }

  /*
    printf("myrand_util WYETH HERE  tap %d  ord %d   nm %d\n",tap_seed,ord,nm);
  */

  ms_r = 1;
  ms_q = 1 << (ord-1); /* Put a 1 in the "ord"th position */
  ms_t = tap_seed;

  s = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++){
    ms_b = ms_r & 1;
    if ((ms_r & ms_q) != ms_q) /* If n-1 bit is not set */
      ms_r = ms_r << 1;
    else
      ms_r = (ms_r << 1) ^ ms_t;
    
    if (ms_b == 0)
      s[i] = 0;
    else
      s[i] = 1;
  }
  
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MYRAND_GET_STD_GAUSS_SEQ                         */
/*                                                                           */
/*****************************************************************************/
float *myrand_get_std_gauss_seq(n,seed)
     int n,seed;
{
  int i;
  float *s;
  
  if (seed == 0)
    exit_error("MYRAND_GET_STD_GAUSS_SEQ","seed = 0");
  else if (seed > 0)
    seed *= -1;
  
  s = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    s[i] = nr_util_gasdev(&seed);   // ...gasdev calls nr_util_ran2
  
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MYRAND_GET_STD_BIN_FLOAT_2D                       */
/*                                                                           */
/*****************************************************************************/
float **myrand_get_std_bin_float_2d(xn,yn,seed)
     int xn,yn,seed;
{
  int i,j,k;
  int *seqi;
  float **f2d;

  seqi = myrand_get_std_bin_seq(xn*yn,seed);

  f2d = get_2d_farray(xn,yn);

  k = 0;
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      if (seqi[k] == 0)
	f2d[i][j] = -1.0;
      else
	f2d[i][j] = 1.0;
      k += 1;
    }
  }
  myfree(seqi);

  return f2d;
}
/**************************************-**************************************/
/*                                                                           */
/*                        MYRAND_GET_STD_QUAD_FLOAT_2D                       */
/*                                                                           */
/*  Returned values are:  -1, -1/3, 1/3, 1                                   */
/*                                                                           */
/*****************************************************************************/
float **myrand_get_std_quad_float_2d(xn,yn,seed)
     int xn,yn,seed;
{
  int i,j,k;
  int *seqi;
  float **f2d;

  seqi = myrand_get_std_quad_seq(xn*yn,seed);

  f2d = get_2d_farray(xn,yn);

  k = 0;
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      if (seqi[k] == 0)
	f2d[i][j] = -1.0;
      else if (seqi[k] == 1)
	f2d[i][j] = -1.0/3.0;
      else if (seqi[k] == 2)
	f2d[i][j] = 1.0/3.0;
      else
	f2d[i][j] = 1.0;

      k += 1;
    }
  }
  myfree(seqi);

  return f2d;
}
/**************************************-**************************************/
/*                                                                           */
/*                       GET_DEVIATE_LIST_FROM_DESCRIPTION                   */
/*                                                                           */
/*  Return a random number based on the given description.                   */
/*                                                                           */
/*  Distribution discription formats:                                        */
/*                                                                           */
/*  uniform_int <seed> <min> <max>                                           */
/*  - Generate a random integer between the <min> and <max> values,          */
/*    inclusive.                                                             */
/*  uniform_float <seed> <min> <max>                                         */
/*  - Generate a random float between the <min> and <max> values,            */
/*    including <min> but not <max>.                                         */
/*  uniform_list <seed> <item1> <item2> ...                                  */
/*  - Select randomly between the items in the list.                         */
/*                                                                           */
/*****************************************************************************/
char **get_deviate_list_from_description(desc,n)
     char *desc;
     int n;
{
  int i,k;
  int ns,seed;
  float min,range;
  char **slist,**ranlist,tstr[SLEN];

  ranlist = (char **)myalloc(n*sizeof(char *));
  get_items_from_string(desc,&slist,&ns);
  if (strcmp(slist[0],"uniform_int")==0){
    seed = atoi(slist[1]);
    min = atof(slist[2]);
    range = atof(slist[3]) - min;
    if (seed > 0)
      seed *= -1;
    for(i=0;i<n;i++){
      sprintf(tstr,"%d",my_rint(min + range*myrand_util_ran2(&seed)));
      ranlist[i] = strdup(tstr);
    }
  }else if (strcmp(slist[0],"uniform_float")==0){
    seed = atoi(slist[1]);
    min = atof(slist[2]);
    range = atof(slist[3]) - min;
    if (seed > 0)
      seed *= -1;
    for(i=0;i<n;i++){
      sprintf(tstr,"%.6f",min + range*myrand_util_ran2(&seed));
      ranlist[i] = strdup(tstr);
    }
  }else if (strcmp(slist[0],"uniform_list")==0){
    seed = atoi(slist[1]);
    min = 0.0;
    range = (float)(ns-2); /* Subtract for name and seed. */
    if (seed > 0)
      seed *= -1;
    for(i=0;i<n;i++){
      k = my_rint(min + range*myrand_util_ran2(&seed));
      ranlist[i] = strdup(slist[2+k]);
    }
  }else
    exit_error("GET_DEVIATE_LIST_FROM_DESCRIPTION","Unknown distrib. type");

  free_2d_carray(slist,ns);

  return ranlist;
}
/**************************************-**************************************/
/*                                                                           */
/*                            EJRAND_StimRandShort                           */
/*                                                                           */
/*  EJ's code, via Greg H's lab.                                             */
/*                                                                           */
/*  Greg says:  My attempt to implement EJ's random number generation        */
/*    routine as a MEX file.  Hopefully this will help in speeding           */
/*    up getWhtnsStim.m.  GDLH 7/23/01                                       */
/*                                                                           */
/*  This replicates the original Macintosh Toolbox Random() routine.         */
/*  The funciton is defined here explicitly for portability                  */
/*  and independence from changes in the MacOS.  EJC 1999-12-22              */
/*  return value formerly 'short'                                            */
/*                                                                           */
/*  Returned values are (almost) evenly distributed between                  */
/*  0 and 65535 except that 32768 is skipped.                                */
/*                                                                           */
/*****************************************************************************/
int StimRandShort(int *seed)
{
  int temp1,temp2,temp3,result;

  temp1 = (*seed & 0xFFFF) * 0x41A7;
  temp2 = (*seed >> 16) * 0x41A7 + (temp1 >> 16);
  temp3 = temp2 * 2 >> 16;
  temp1 = (temp1 & 0xFFFF) - 0x7FFFFFFF;
  temp2 &= 0x7FFF;
  temp1 += (temp2 << 16) | (temp2 >> 16) + temp3;
  if (temp1 < 0)
    temp1 += 0x7FFFFFFF;
  *seed = temp1;
  result = temp1 & 0xFFFF;  // W: Limit result to at most 65,535 (2^16 - 1)
  if (result == 0x8000)     // W: Why is this number chosen to be zero?
    result = 0;

  return result;
}
/**************************************-**************************************/
/*                                                                           */
/*                               EJRAND_GETNUMS                              */
/*                                                                           */
/*  EJ's code, via Greg H's lab.                                             */
/*                                                                           */
/*****************************************************************************/
int ejrand_getnums(int *op,    // [iter] Fill this array with random numbers
		   int iter,   // The number of random nums to generate
		   int seed)   // The initial seed to generate the sequence
{
  int i;

  for (i=0;i<iter;i++)
    op[i] = StimRandShort(&seed);

  return seed;
}
/**************************************-**************************************/
/*                                                                           */
/*                        MYRAND_GET_STD_QUAD_EJRAND_2D                      */
/*                                                                           */
/*  Returned values are:  -1, -1/3, 1/3, 1                                   */
/*                                                                           */
/*****************************************************************************/
float **myrand_get_std_quad_ejrand_2d(xn,yn,seed)
     int xn,yn,seed;
{
  int i,j,jj,k;
  int *seqi,tseed;
  float **f2d;

  seqi = (int *)myalloc(xn*yn*sizeof(int));  // Storage for 1D sequence

  tseed = ejrand_getnums(seqi,xn*yn,seed);  // Fill 'seqi' (ignore 'tseed')

  f2d = get_2d_farray(xn,yn);

  //  *** I flipped the y-axis order here to match Patrick's Matlab analysis.

  k = 0;
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      jj = yn-1 - j;   // Flip the y-axis order
      if (seqi[k] % 4 == 0)
	f2d[i][jj] = -1.0;
      else if (seqi[k] % 4 == 1)
	f2d[i][jj] = -1.0/3.0;
      else if (seqi[k] % 4 == 2)
	f2d[i][jj] = 1.0/3.0;
      else
	f2d[i][jj] = 1.0;

      //printf("seqi[k] = %d   f2d[i][j] = %f \n",seqi[k],f2d[i][j]);

      k += 1;
    }
  }
  myfree(seqi);

  return f2d;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MYRAND_EJRAND_TEST                            */
/*                                                                           */
/*****************************************************************************/
void myrand_ejrand_test(n,seed)
     int n;      // Compute this many random numbers
     int seed;   // Starting seed value
{
  int i,k;
  int maxval,mod4[4];
  float *hist;

  printf("  MYRAND_EJRAND_TEST\n");
  printf("    Histogram plot of %d random values for seed %d\n",n,seed);

  maxval = 65536;
  hist = get_zero_farray(maxval);
  mod4[0] = mod4[1] = mod4[2] = mod4[3] = 0;

  for(i=0;i<n;i++){
    k = StimRandShort(&seed);
    if ((k < 0) || (k >= maxval)){
      exit_error("MYRAND_EJRAND_TEST","Value out of bounds");
    }
    hist[k] += 1.0;
    mod4[k%4] += 1;
  }

  printf("   0:  %d\n",mod4[0]);
  printf("   1:  %d\n",mod4[1]);
  printf("   2:  %d\n",mod4[2]);
  printf("   3:  %d\n",mod4[3]);

  append_farray_plot("zzz.ejrand.hist.pl","hist",hist,maxval,1);

  printf("   Done.  Exiting.\n\n");

  exit(0);
}
/**************************************-**************************************/
/*                                                                           */
/*                              EJRAND_GET_SEEDS                             */
/*                                                                           */
/*  Written by Wyeth:  this will extract frame seeds for 'stimx' style       */
/*  stimulus generation.                                                     */
/*                                                                           */
/*****************************************************************************/
int *ejrand_get_seeds(int seed,  // Overall initial seed
		      int vpf,   // Values per frame, e.g. 121 for 11 x 11
		      int n)     // Number of frames
{
  int i;
  int *slist,*t;

  slist = (int *)myalloc(n*sizeof(int));    // List to keep 1st seed values
  t = (int *)myalloc(vpf*sizeof(int));      // Temp list of values

  for (i=0;i<n;i++){
    slist[i] = seed;
    seed = ejrand_getnums(t,vpf,slist[i]);
  }

  myfree(t);

  return slist;
}


/*************************************---*************************************/
// Period parameters
#define GENRAND_N 624
#define GENRAND_M 397
#define MATRIX_A 0x9908b0dfUL   // constant vector a
#define UPPER_MASK 0x80000000UL // most significant w-r bits
#define LOWER_MASK 0x7fffffffUL // least significant r bits

static unsigned long genrand_mt[GENRAND_N];   // array for the state vector
static int           genrand_mti=GENRAND_N+1; // mti==N+1 means mt[N] not init
/**************************************-**************************************/
/*                                                                           */
/*                                 INIT_GENRAND                              */
/*                                                                           */
/*  WYETH:  On Oct 2, 2019, I downloaded and reformmated this code from:     */
/*  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html   */
/*  This notices below apply to all functions below with 'genrand' in name.  */
/*                                                                           */
/*   A C-program for MT19937, with initialization improved 2002/1/26.        */
/*   Coded by Takuji Nishimura and Makoto Matsumoto.                         */
/*                                                                           */
/*   Before using, initialize the state by using init_genrand(seed)          */
/*   or genrand_init_by_array(init_key, key_length).                         */
/*                                                                           */
/*   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,       */
/*   All rights reserved.                                                    */
/*                                                                           */
/*  Redistribution and use in source and binary forms, with or without       */
/*  modification, are permitted provided that the following conditions       */
/*  are met:                                                                 */
/*                                                                           */
/*   1. Redistributions of source code must retain the above copyright       */
/*      notice, this list of conditions and the following disclaimer.        */
/*   2. Redistributions in binary form must reproduce the above copyright    */
/*      notice, this list of conditions and the following disclaimer in the  */
/*      documentation and/or other materials provided with the distribution. */
/*   3. The names of its contributors may not be used to endorse or promote  */
/*      products derived from this software without specific prior written   */
/*      permission.                                                          */
/*                                                                           */
/*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS  */
/*  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED    */
/*  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTIC-  */
/*  ULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR   */
/*  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,    */
/*  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,      */
/*  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR       */
/*  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF   */
/*  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING     */
/*  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS       */
/*  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             */
/*                                                                           */
/*  Any feedback is very welcome.                                            */
/*  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html                 */
/*  email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)                 */
/*                                                                           */
/*  Note, the real... versions below are due to Isaku Wada, 2002/01/09 added */
/*                                                                           */
/*****************************************************************************/
void genrand_init(unsigned long s){  // initializes genrand_mt[N] with a seed
  genrand_mt[0]= s & 0xffffffffUL;
  for (genrand_mti=1; genrand_mti<GENRAND_N; genrand_mti++) {
    genrand_mt[genrand_mti] = 
      (1812433253UL * (genrand_mt[genrand_mti-1] ^ (genrand_mt[genrand_mti-1] >> 30)) + genrand_mti); 
    // See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
    // In the previous versions, MSBs of the seed affect
    // only MSBs of the array genrand_mt[].
    // 2002/01/09 modified by Makoto Matsumoto
    genrand_mt[genrand_mti] &= 0xffffffffUL;
    // for >32 bit machines
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            GENRAND_INIT_BY_ARRAY                          */
/*                                                                           */
/*  ORIGINAL NOTES FROM DISTRIBUTION:                                        */
/*    Initialize by an array with array-length                               */
/*    init_key is the array for initializing keys                            */ 
/*    key_length is its length                                               */
/*    slight change for C++, 2004/2/26                                       */
/*                                                                           */
/*****************************************************************************/
void genrand_init_by_array(unsigned long init_key[],
			   int key_length){
  int i, j, k;

  genrand_init(19650218UL);
  i=1; j=0;
  k = (GENRAND_N>key_length ? GENRAND_N : key_length);

  for (; k; k--) {
    genrand_mt[i] = (genrand_mt[i] ^
		     ((genrand_mt[i-1] ^ (genrand_mt[i-1] >> 30)) *
		      1664525UL)) + init_key[j] + j; // non linear
    genrand_mt[i] &= 0xffffffffUL; // for WORDSIZE > 32 machines
    i++; j++;
    if (i>=GENRAND_N) { genrand_mt[0] = genrand_mt[GENRAND_N-1]; i=1; }
    if (j>=key_length) j=0;
  }

  for (k=GENRAND_N-1; k; k--) {
    genrand_mt[i] = (genrand_mt[i] ^
		     ((genrand_mt[i-1] ^ (genrand_mt[i-1] >> 30)) *
		      1566083941UL))- i; // non linear
    genrand_mt[i] &= 0xffffffffUL; // for WORDSIZE > 32 machines
    i++;
    if (i>=GENRAND_N) { genrand_mt[0] = genrand_mt[GENRAND_N-1]; i=1; }
  }

  genrand_mt[0] = 0x80000000UL; // MSB is 1; assuring non-zero initial array
}
/**************************************-**************************************/
/*                                                                           */
/*                               GENRAND_INT32                               */
/*                                                                           */
/*  Generates a random number on [0,0xffffffff]-interval                     */
/*                                                                           */
/*****************************************************************************/
unsigned long genrand_int32(void){
  unsigned long y;
  static unsigned long mag01[2]={0x0UL, MATRIX_A};
  // mag01[x] = x * MATRIX_A  for x=0,1

  if (genrand_mti >= GENRAND_N) { // generate N words at one time
    int kk;

    if (genrand_mti == GENRAND_N+1) // if genrand_init() has not been called,
      genrand_init(5489UL);         //   a default initial seed is used

    for (kk=0;kk<GENRAND_N-GENRAND_M;kk++) {
      y = (genrand_mt[kk]&UPPER_MASK)|(genrand_mt[kk+1]&LOWER_MASK);
      genrand_mt[kk] = genrand_mt[kk+GENRAND_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for (;kk<GENRAND_N-1;kk++) {
      y = (genrand_mt[kk]&UPPER_MASK)|(genrand_mt[kk+1]&LOWER_MASK);
      genrand_mt[kk] = genrand_mt[kk+(GENRAND_M-GENRAND_N)] ^ (y >> 1)
	^ mag01[y & 0x1UL];
    }
    y = (genrand_mt[GENRAND_N-1]&UPPER_MASK)|(genrand_mt[0]&LOWER_MASK);
    genrand_mt[GENRAND_N-1] = genrand_mt[GENRAND_M-1] ^ (y >> 1)
      ^ mag01[y & 0x1UL];

    genrand_mti = 0;
  }
  
  y = genrand_mt[genrand_mti++];

  // Tempering
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  return y;
}
/**************************************-**************************************/
/*                                                                           */
/*                               GENRAND_INT31                               */
/*                                                                           */
/*  Generates a random number on [0,0x7fffffff]-interval                     */
/*                                                                           */
/*****************************************************************************/
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}
/**************************************-**************************************/
/*                                                                           */
/*                               GENRAND_REAL1                               */
/*                                                                           */
/*  Generates a random number on [0,1]-real-interval                         */
/*                                                                           */
/*****************************************************************************/
double genrand_real1(void)
{
  return genrand_int32()*(1.0/4294967295.0);  // divided by 2^32-1
}
/**************************************-**************************************/
/*                                                                           */
/*                               GENRAND_REAL2                               */
/*                                                                           */
/*  Generates a random number on [0,1)-real-interval                         */
/*                                                                           */
/*****************************************************************************/
double genrand_real2(void){
  return genrand_int32()*(1.0/4294967296.0);   // divided by 2^32
}
/**************************************-**************************************/
/*                                                                           */
/*                               GENRAND_REAL3                               */
/*                                                                           */
/*  Generates a random number on (0,1)-real-interval                         */
/*                                                                           */
/*****************************************************************************/
double genrand_real3(void){
  return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); // div by 2^32
}
/**************************************-**************************************/
/*                                                                           */
/*                               GENRAND_RES53                               */
/*                                                                           */
/*  Generates a random number on [0,1) with 53-bit resolution                */
/*                                                                           */
/*****************************************************************************/
double genrand_res53(void){ 
  unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
  return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/*
int main(void)
{
    int i;
    unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
    genrand_init_by_array(init, length);
    printf("1000 outputs of genrand_int32()\n");
    for (i=0; i<1000; i++) {
      printf("%10lu ", genrand_int32());
      if (i%5==4) printf("\n");
    }
    printf("\n1000 outputs of genrand_real2()\n");
    for (i=0; i<1000; i++) {
      printf("%10.8f ", genrand_real2());
      if (i%5==4) printf("\n");
    }
    return 0;
}
*/

/**************************************-**************************************/
/*                                                                           */
/*                             GR_SHUFFLE_IARRAY                             */
/*                                                                           */
/*  Shuffle an array of ints.                                                */
/*                                                                           */
/*****************************************************************************/
void gr_shuffle_iarray(data,n,nshuffle,seed)
     int *data,n,nshuffle;
     int seed;
{
  int i,j;
  float *r; // random numbers to sort

  genrand_init((unsigned long)seed);

  r = (float *)myalloc(n*sizeof(float));
  for (i=0;i<nshuffle;i++){
    for (j=0;j<n;j++)
      r[j] = genrand_real1();
    hpsort((unsigned long)n,r-1,data-1);
  }
  myfree(r);
}
/**************************************-**************************************/
/*                                                                           */
/*                            GR_GET_SHUFFLE_INDEX                           */
/*                                                                           */
/*  Return a shuffled list of the numbers from 0 to n-1.                     */
/*                                                                           */
/*****************************************************************************/
int *gr_get_shuffle_index(n,nshuffle,seed)
     int n,nshuffle;
     int seed;
{
  int i;
  int *s;

  s = (int *)myalloc(n*sizeof(int));
  for (i=0;i<n;i++)
    s[i] = i;
  gr_shuffle_iarray(s,n,nshuffle,seed);
  return s;
}
