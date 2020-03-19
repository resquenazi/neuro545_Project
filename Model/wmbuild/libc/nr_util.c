/*****************************************************************************/
/*                                                                           */
/*  nr_util.c                                                                */
/*  wyeth bair                                                               */
/*  caltech                                                                  */
/*  03/12/93                                                                 */
/*                                                                           */
/*  This file contains utilities from Numerical Recipes in C.  Some have     */
/*  been modified.                                                           */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#include "my_util.h"

/**************************************-**************************************/
/*                                                                           */
/*                                  NRERROR                                  */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void nrerror(error_text)
     char error_text[];
{
  printf("Numerical Recipes run-time error...\n");
  printf("%s\n",error_text);
  printf("...now exiting to system...\n");
  exit(1);
}
/**************************************-**************************************/
/*                                                                           */
/*                               FREE_F3TENSOR                               */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void free_f3tensor(t,nrl,nrh,ncl,nch,ndl,ndh)
     float ***t;
     long nrl,nrh,ncl,nch,ndl,ndh;
{
  free((char*) (t[nrl][ncl]+ndl-1));
  free((char*) (t[nrl]+ncl-1));
  free((char*) (t+nrl-1));
}
/**************************************-**************************************/
/*                                                                           */
/*                                  F3TENSOR                                 */
/*                                                                           */
/*  allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh]      */
/*                                                                           */
/*  The float storage is a single 1D array.                                  */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float ***f3tensor(nrl,nrh,ncl,nch,ndl,ndh)
     long nrl,nrh;  // row low, hi
     long ncl,nch;  // col low, hi
     long ndl,ndh;  // depth lo, hi
{
  //long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  long i,j,nrow,ncol,ndep;
  float ***t;

  nrow = nrh-nrl+1;
  ncol = nch-ncl+1;
  ndep = ndh-ndl+1;

  // allocate pointers to pointers to rows
  t = (float ***)malloc((size_t)((nrow+1)*sizeof(float **)));
  if (!t) nrerror("allocation failure 1 in f3tensor()");
  t += 1;
  t -= nrl;

  // allocate pointers to rows and set pointers to them
  t[nrl]=(float **) malloc((size_t)((nrow*ncol+1)*sizeof(float*)));
  if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
  t[nrl] += 1;
  t[nrl] -= ncl;

  // allocate rows and set pointers to them
  t[nrl][ncl]=(float *)malloc((size_t)((nrow*ncol*ndep+1)*sizeof(float)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
  t[nrl][ncl] += 1;
  t[nrl][ncl] -= ndl;

  for(j=ncl+1;j<=nch;j++)
    t[nrl][j] = t[nrl][j-1] + ndep;

  for(i=nrl+1;i<=nrh;i++){
    t[i] = t[i-1] + ncol;
    t[i][ncl] = t[i-1][ncl] + ncol*ndep;
    for(j=ncl+1;j<=nch;j++)
      t[i][j] = t[i][j-1] + ndep;
  }

  // return pointer to array of pointers to rows
  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_F3TENSOR_FROM_DATA                          */
/*                                                                           */
/*  Take data from a 3d array and make a NumRec f3tensor.                    */
/*                                                                           */
/*****************************************************************************/
void get_f3tensor_from_data(data,xn,yn,tn,rtdata,x0,xd,y0,yd,t0,td)
     float ***data;
     int xn,yn,tn;
     float ****rtdata;
     int x0,xd,y0,yd,t0,td; // start and duration in each dimension
{
  int i,j,k;
  float ***tdata;

  tdata = f3tensor(1,xd,1,yd,1,td);
  for(i=1;i<=xd;i++)
    for(j=1;j<=yd;j++)
      for(k=1;k<=td;k++)
	tdata[i][j][k] = data[x0+i-1][y0+j-1][t0+k-1];

  *rtdata = tdata;
}
/**************************************-**************************************/
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void free_matrix(m,nrl,nrh,ncl,nch)
     float **m;
     long nrl,nrh,ncl,nch;
{
  free((char*) (m[nrl]+ncl-1));
  free((char*) (m+nrl-1));
}
/**************************************-**************************************/
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float **matrix(nrl,nrh,ncl,nch)
     long nrl,nrh,ncl,nch;
     /* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float **m;
  
  /* allocate pointers to rows */
  m=(float **) malloc((size_t)((nrow+1)*sizeof(float*)));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += 1;
  m -= nrl;
  
  /* allocate rows and set pointers to them */
  m[nrl]=(float *) malloc((size_t)((nrow*ncol+1)*sizeof(float)));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += 1;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}
/**************************************-**************************************/
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float *vector(nl,nh)
     int nl,nh;
{
  float *v;
  
  v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
  if (!v) nrerror("allocation failure in vector()");
  return v-nl;
}
/**************************************-**************************************/
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void free_vector(v,nl,nh)
     float *v;
     int nl,nh;
{
  free((char*) (v+nl));
}
/**************************************-**************************************/
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
int *ivector(nl,nh)
     long nl,nh;
     /* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;
  
  v=(int *)malloc((size_t) ((nh-nl+1+1)*sizeof(int)));
  if (!v) nrerror("allocation failure in ivector()");
  return v-nl+1;
}
/**************************************-**************************************/
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void free_ivector(v,nl,nh)
     int *v;
     long nl,nh;
     /* free an int vector allocated with ivector() */
{
  free((char*) (v+nl-1));
}
/**************************************-**************************************/
/*                                                                           */
/*                                   AVEVAR                                  */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void avevar(data,n,ave,svar)
     float data[],*ave,*svar;
     int n;
{
  int j;
  float s;
  
  *ave=(*svar)=0.0;
  for (j=1;j<=n;j++) *ave += data[j];
  *ave /= n;
  for (j=1;j<=n;j++) {
    s=data[j]-(*ave);
    *svar += s*s;
  }
  *svar /= (n-1);
}
/**************************************-**************************************/
/*                                                                           */
/*                                   ERFCC                                   */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float erfcc(x)
     float x;
{
  float t,z,ans;
  
  z=fabs(x);
  t=1.0/(1.0+0.5*z);
  ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
	    t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
            t*(-0.82215223+t*0.17087277)))))))));
  return x >= 0.0 ? ans : 2.0-ans;
}
/**************************************-**************************************/
/*                                                                           */
/*                                  GAMMLN                                   */
/*                                                                           */
/*C     (C) Copr. 1986-92 Numerical Recipes Software #.3.                    */
/*                                                                           */
/*****************************************************************************/
float gammln(xx)
     float xx;
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			  24.01409824083091,-1.231739572450155,
			  0.1208650973866179e-2,-0.5395239384953e-5};
  int j;

  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}
/*************************************---*************************************/
#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
/**************************************-**************************************/
/*                                                                           */
/*                                   BETACF                                  */
/*                                                                           */
/*  Continued fraction used by "betai".                                      */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float betacf(a,b,x)
     float a,b,x;
{
  int m,m2;
  float aa,c,d,del,h,qab,qam,qap;
  
  qab=a+b;
  qap=a+1.0;
  qam=a-1.0;
  c=1.0;
  d=1.0-qab*x/qap;
  if (fabs(d) < FPMIN) d=FPMIN;
  d=1.0/d;
  h=d;
  for (m=1;m<=MAXIT;m++) {
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  if (m > MAXIT) nrerror("a or b too big, or MAXIT too small in betacf");
  return h;
}
#undef MAXIT
#undef EPS
#undef FPMIN
/**************************************-**************************************/
/*                                                                           */
/*                                   BETAI                                   */
/*                                                                           */
/*  Incomplete beta function.                                                */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float betai(a,b,x)
     float a,b,x;
{
  float betacf();
  float gammln();
  float bt;
  
  if (x < 0.0 || x > 1.0) nrerror("Bad x in routine betai");
  if (x == 0.0 || x == 1.0) bt=0.0;
  else
    bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
  if (x < (a+1.0)/(a+b+2.0))
    return bt*betacf(a,b,x)/a;
  else
    return 1.0-bt*betacf(b,a,1.0-x)/b;
}
/*************************************---*************************************/
#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
/**************************************-**************************************/
/*                                                                           */
/*                                  MYBETCF                                  */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float mybetacf(a,b,x,flag)
     float a,b,x;
     int *flag;
{
  int m,m2;
  float aa,c,d,del,h,qab,qam,qap;
  
  qab=a+b;
  qap=a+1.0;
  qam=a-1.0;
  c=1.0;
  d=1.0-qab*x/qap;
  if (fabs(d) < FPMIN) d=FPMIN;
  d=1.0/d;
  h=d;
  for (m=1;m<=MAXIT;m++) {
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  if (m > MAXIT)
    *flag = 1;
  else
    *flag = 0;

  return h;
}
#undef MAXIT
#undef EPS
#undef FPMIN
/**************************************-**************************************/
/*                                                                           */
/*                                  MYBETAI                                  */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float mybetai(a,b,x)
     float a,b,x;
{
  float betacf();
  float gammln();
  float bt,bcf;
  int flag;
  
  if (x < 0.0 || x > 1.0) nrerror("Bad x in routine betai");
  if (x == 0.0 || x == 1.0) bt=0.0;
  else
    bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
  if (x < (a+1.0)/(a+b+2.0)){
    bcf = mybetacf(a,b,x,&flag);
    if (flag)
      return -1.0;
    else
      return bt*bcf/a;
  }else{
    bcf = mybetacf(b,a,1.0-x,&flag);
    if (flag)
      return -1.0;
    else
      return 1.0-bt*bcf/b;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                   CRANK                                   */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void crank(n,w,s)
     float w[],*s;
     int n;
{
  int j=1,ji,jt;
  float t,rank;
  
  *s=0.0;
  while (j < n) {
    if (w[j+1] != w[j]) {
      w[j]=j;
      ++j;
    } else {
      for (jt=j+1;jt<=n;jt++)
	if (w[jt] != w[j]) break;
      rank=0.5*(j+jt-1);
      for (ji=j;ji<=(jt-1);ji++) w[ji]=rank;
      t=jt-j;
      *s += t*t*t-t;
      j=jt;
    }
  }
  if (j == n) w[n]=n;
}
/*************************************---*************************************/
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
/**************************************-**************************************/
/*                                                                           */
/*                                   RAN1                                    */
/*                                                                           */
/*  Returns a uniform random deviate between 0.0 and 1.0.  Set "idum" to     */
/*  any negative value to initialize or reinitialize the sequence.           */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float ran1(idum)
     long *idum;
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  
  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
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
/*                                NR_UTIL_RAN2                               */
/*                                                                           */
/*  To be used within one routine, i.e., when no other routine could inter-  */
/*  fer by trying to reseed.                                                 */
/*                                                                           */
/*  COPIED FROM NumRecInC, Second Edition, page 282.                         */
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
float nr_util_ran2(idum)
     int *idum; /*** Changed from long for consistency on SUNS, DEC ALPHA ***/
{
  int j;
  int k; /*** Changed from long to int. ***/
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
/*                            NR_UTIL_GASDEV  (p288-290)                     */
/*                                                                           */
/*  Returns a normally distributed deviate with zero mean and unit variance  */
/*  using ran2(idum) as the source of uniform deviates.                      */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3                        */
/*                                                                           */
/*****************************************************************************/
float nr_util_gasdev(idum)
     int *idum;
{
  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;
  
  /*if (iset == 0){*/
  if ((iset == 0) || (*idum <= 0)){  /* Wyeth, reset state if reseeded */
    do {
      v1=2.0*nr_util_ran2(idum)-1.0;
      v2=2.0*nr_util_ran2(idum)-1.0;
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
/*************************************---*************************************/
#define NRANSI
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50
/**************************************-**************************************/
/*                                                                           */
/*                                   SORT                                    */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void sort(n,arr)
     unsigned long n;
     float arr[];
{
  unsigned long i,ir=n,j,k,l=1;
  int jstack=0,*istack;
  float a,temp;
  
  istack=ivector(1,NSTACK);
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	a=arr[j];
	for (i=j-1;i>=1;i--) {
	  if (arr[i] <= a) break;
	  arr[i+1]=arr[i];
	}
	arr[i+1]=a;
      }
      if (jstack == 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      SWAP(arr[k],arr[l+1])
	if (arr[l+1] > arr[ir]) {
	  SWAP(arr[l+1],arr[ir])
	    }
      if (arr[l] > arr[ir]) {
	SWAP(arr[l],arr[ir])
	  }
      if (arr[l+1] > arr[l]) {
	SWAP(arr[l+1],arr[l])
	  }
      i=l+1;
      j=ir;
      a=arr[l];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i],arr[j]);
      }
      arr[l]=arr[j];
      arr[j]=a;
      jstack += 2;
      if (jstack > NSTACK) nrerror("NSTACK too small in sort.");
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
  free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI
/**************************************-**************************************/
/*                                                                           */
/*                                   SORT2                                   */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void sort2(n,ra,rb)
     int n;
     float ra[],rb[];
{
  int l,j,ir,i;
  float rrb,rra;
  
  l=(n >> 1)+1;
  ir=n;
  for (;;){
    if (l > 1){
      rra=ra[--l];
      rrb=rb[l];
    }else{
      rra=ra[ir];
      rrb=rb[ir];
      ra[ir]=ra[1];
      rb[ir]=rb[1];
      if (--ir == 1){
	ra[1]=rra;
	rb[1]=rrb;
	return;
      }
    }
    i=l;
    j=l << 1;
    while (j <= ir){
      if (j < ir && ra[j] < ra[j+1]) ++j;
      if (rra < ra[j]){
	ra[i]=ra[j];
	rb[i]=rb[j];
	j += (i=j);
      }
      else j=ir+1;
    }
    ra[i]=rra;
    rb[i]=rrb;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                 HPSORT_INT                                */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void hpsort_int(n,ra)
     unsigned long n;
     int ra[];
{
  unsigned long i,ir,j,l;
  int rra;
  
  if (n < 2) return;
  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1) {
      rra=ra[--l];
    } else {
      rra=ra[ir];
      ra[ir]=ra[1];
      if (--ir == 1) {
	ra[1]=rra;
	break;
      }
    }
    i=l;
    j=l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
	ra[i]=ra[j];
	i=j;
	j <<= 1;
      } else j=ir+1;
    }
    ra[i]=rra;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                  GASDEV                                   */
/*                                                                           */
/*  Returns a normally distributed deviate with zero mean and unit           */
/*  variance, using ran1(idum) as the source of uniform deviates.            */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float gasdev(idum)
     long *idum;
{
  float ran1();
  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;
  
  /*if  (iset == 0) {*/
  if ((iset == 0) || (*idum <= 0)){  /* Wyeth, reset state if reseeded */
    do {
      v1=2.0*ran1(idum)-1.0;
      v2=2.0*ran1(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}
/*************************************---*************************************/
#define TINY 1.0e-20
/**************************************-**************************************/
/*                                                                           */
/*                                   PEARSN                                  */
/*                                                                           */
/*  Given two arrays x[1..n] and y[1..n], this routine computes their        */
/*  correlation coefficient "r", the significance level at which the null    */
/*  hypothesis of zero correlation is disproved ("prob" whose small value    */
/*  indicates a significant correlation), and Fisher's "z", whose value      */
/*  can be used in further statistical tests as described in NumRecInC.      */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void pearsn(x,y,n,r,prob,z)
     float x[],y[];
     unsigned long n;
     float *r,*prob,*z;
{
  float betai(); //(float a, float b, float x);
  float erfcc(); //float x);
  unsigned long j;
  float yt,xt,t,df;
  float syy=0.0,sxy=0.0,sxx=0.0,ay=0.0,ax=0.0;
  
  for (j=1;j<=n;j++) {
    ax += x[j];
    ay += y[j];
  }
  ax /= n;
  ay /= n;
  for (j=1;j<=n;j++) {
    xt=x[j]-ax;
    yt=y[j]-ay;
    sxx += xt*xt;
    syy += yt*yt;
    sxy += xt*yt;
  }
  *r=sxy/sqrt(sxx*syy);
  *z=0.5*log((1.0+(*r)+TINY)/(1.0-(*r)+TINY));
  df=n-2;
  t=(*r)*sqrt(df/((1.0-(*r)+TINY)*(1.0+(*r)+TINY)));
  *prob=betai(0.5*df,0.5,df/(df+t*t));
}
#undef TINY
/**************************************-**************************************/
/*                                                                           */
/*                               SIMPLE_PEARSN                               */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void simple_pearsn(x,y,n,r) // Modification by Wyeth to ignore stat.
     float x[],y[];
     unsigned long n;
     float *r;
{
  unsigned long j;
  float yt,xt;
  float syy=0.0,sxy=0.0,sxx=0.0,ay=0.0,ax=0.0;

  for (j=1;j<=n;j++){
    ax += x[j];
    ay += y[j];
  }
  ax /= n;
  ay /= n;
  for (j=1;j<=n;j++){
    xt=x[j]-ax;
    yt=y[j]-ay;
    sxx += xt*xt;
    syy += yt*yt;
    sxy += xt*yt;
  }
  *r=sxy/sqrt(sxx*syy);
}
/**************************************-**************************************/
/*                                                                           */
/*                               SIMPLE_PEARSN0                              */
/*                                                                           */
/*   THIS VERSION USES ARRAYS FROM 0..n-1                                    */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void simple_pearsn0(x,y,n,r) // Modification by Wyeth to ignore stat
     float x[],y[];
     unsigned long n;
     float *r;
{
  unsigned long j;
  float yt,xt;
  float syy=0.0,sxy=0.0,sxx=0.0,ay=0.0,ax=0.0;
  
  for (j=0;j<n;j++) {
    ax += x[j];
    ay += y[j];
  }
  ax /= n;
  ay /= n;
  for (j=0;j<n;j++) {
    xt=x[j]-ax;
    yt=y[j]-ay;
    sxx += xt*xt;
    syy += yt*yt;
    sxy += xt*yt;
  }
  *r=sxy/sqrt(sxx*syy);
}
/*************************************---*************************************/
#define NRANSI
/**************************************-**************************************/
/*                                                                           */
/*                                   SPEAR                                   */
/*                                                                           */
/*  Given two data arrays, data1[1..n] and data2[1..n], this routine         */
/*  returns their sum-squared difference of ranks as "d", the number of      */
/*  standard deviations by which "d" deviates from its null hypothesis       */
/*  expected value as "zd", the two-sided significance level of this         */
/*  deviation as "probd", Spearmans's rank correlation "rs", and the         */
/*  two-sided significance level of its devitaion from zero as "probrs".     */
/*  The external routines "crank" and "sort2" are used.  A small value of    */
/*  either "probd" or "probrs" indicates a significant correlation ("rs"     */
/*  positive) or anticorrelation ("rs" negative).                            */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void spear(data1,data2,n,d,zd,probd,rs,probrs)
     float data1[],data2[];
     unsigned long n;
     float *d,*zd,*probd,*rs,*probrs;
{
  float betai();
  void crank();
  float erfcc();
  void sort2();
  unsigned long j;
  float vard,t,sg,sf,fac,en3n,en,df,aved,*wksp1,*wksp2;
  
  wksp1=vector(1,n);
  wksp2=vector(1,n);
  for (j=1;j<=n;j++) {
    wksp1[j]=data1[j];
    wksp2[j]=data2[j];
  }
  sort2(n,wksp1,wksp2);
  crank(n,wksp1,&sf);
  sort2(n,wksp2,wksp1);
  crank(n,wksp2,&sg);
  *d=0.0;
  for (j=1;j<=n;j++)
    *d += (wksp1[j]-wksp2[j]) * (wksp1[j]-wksp2[j]);
  //*d += SQR(wksp1[j]-wksp2[j]);
  en=n;
  en3n=en*en*en-en;
  aved=en3n/6.0-(sf+sg)/12.0;
  fac=(1.0-sf/en3n)*(1.0-sg/en3n);
  vard=((en-1.0)*en*en*(en+1.0)*(en+1.0)/36.0)*fac;
  //vard=((en-1.0)*en*en*SQR(en+1.0)/36.0)*fac;
  *zd=(*d-aved)/sqrt(vard);
  *probd=erfcc(fabs(*zd)/1.4142136);
  *rs=(1.0-(6.0/en3n)*(*d+(sf+sg)/12.0))/sqrt(fac);
  fac=(*rs+1.0)*(1.0-(*rs));
  if (fac > 0.0) {
    t=(*rs)*sqrt((en-2.0)/fac);
    df=en-2.0;
    *probrs=betai(0.5*df,0.5,df/(df+t*t));
  } else
    *probrs=0.0;
  free_vector(wksp2,1,n);
  free_vector(wksp1,1,n);
}
#undef NRANSI
/*************************************---*************************************/
#define ITMAX 100
#define EPS 3.0e-7
/**************************************-**************************************/
/*                                                                           */
/*                                   GSER                                    */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void gser(gamser,a,x,gln)
     float *gamser,a,x,*gln;
{
  float gammln();
  int n;
  float sum,del,ap;
  
  *gln=gammln(a);
  if (x <= 0.0) {
    if (x < 0.0) nrerror("x less than 0 in routine gser");
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    nrerror("a too large, ITMAX too small in routine gser");
    return;
  }
}
#undef ITMAX
#undef EPS
/*************************************---*************************************/
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
/**************************************-**************************************/
/*                                                                           */
/*                                    GCF                                    */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void gcf(gammcf,a,x,gln)
     float *gammcf,a,x,*gln;
{
  float gammln();
  int i;
  float an,b,c,d,del,h;
  
  *gln=gammln(a);
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for (i=1;i<=ITMAX;i++) {
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=b+an/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
  *gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN
/**************************************-**************************************/
/*                                                                           */
/*                                   GAMMQ                                   */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float gammq(a,x)
     float a,x;
{
  void gcf();
  void gser();
  float gamser,gammcf,gln;
  
  if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
  if (x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return 1.0-gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return gammcf;
  }
}
/*************************************---*************************************/
#define NRANSI
#define TINY 1.0e-30
/**************************************-**************************************/
/*                                                                           */
/*                                  CNTAB1                                   */
/*                                                                           */
/*  Given a two-dimensional contingency table in the form of an integer      */
/*  array nn[1..ni][1..nj], this routine returns the chi-square "chisq",     */
/*  the number of degrees of freedom "df", the significance level "prob"     */
/*  (small values indicating a significant association), and two measures    */
/*  of asssociation, Cramer's V "cramv" and the contingency coefficient C    */
/*  "ccc".                                                                   */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void cntab1(nn,ni,nj,chisq,df,prob,cramrv,ccc)
     int **nn,ni,nj;
     float *chisq,*df,*prob,*cramrv,*ccc;
{
  float gammq();
  int nnj,nni,j,i,minij;
  float sum=0.0,expctd,*sumi,*sumj,temp;
  
  sumi=vector(1,ni);
  sumj=vector(1,nj);
  nni=ni;
  nnj=nj;
  for (i=1;i<=ni;i++) {
    sumi[i]=0.0;
    for (j=1;j<=nj;j++) {
      sumi[i] += nn[i][j];
      sum += nn[i][j];
    }
    if (sumi[i] == 0.0) --nni;
  }
  for (j=1;j<=nj;j++) {
    sumj[j]=0.0;
    for (i=1;i<=ni;i++) sumj[j] += nn[i][j];
    if (sumj[j] == 0.0) --nnj;
  }
  *df=nni*nnj-nni-nnj+1;
  *chisq=0.0;
  for (i=1;i<=ni;i++) {
    for (j=1;j<=nj;j++) {
      expctd=sumj[j]*sumi[i]/sum;
      temp=nn[i][j]-expctd;
      *chisq += temp*temp/(expctd+TINY);
    }
  }
  *prob=gammq(0.5*(*df),0.5*(*chisq));
  minij = nni < nnj ? nni-1 : nnj-1;
  *cramrv=sqrt(*chisq/(sum*minij));
  *ccc=sqrt(*chisq/(*chisq+sum));
  free_vector(sumj,1,nj);
  free_vector(sumi,1,ni);
}
#undef TINY
#undef NRANSI
/*************************************---*************************************/
#define NRANSI
#define TINY 1.0e-30
/**************************************-**************************************/
/*                                                                           */
/*                                  CNTAB2                                   */
/*                                                                           */
/*  Given a two-dimensional contingency table in the form of an integer      */
/*  array nn[i][j], where i labels the x variable and ranges from 1 to ni,   */
/*  j labels the y variable and ranges from 1 to nj, this routine returns    */
/*  the entropy "h" of the whole table, the entropy "hx" of the x            */
/*  distribution, the entropy "hy" of the y distribution, the entropy        */
/*  "hygx" of y given x, the entropy "hxgy" of x given y, the dependency     */
/*  "uygx" of y on x (eq. 14.4.15), the dependency "uxgy" of x on y (eq.     */
/*  14.4.16), and the symmetrical dependency "uxy" (eq. 14.4.17).            */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void cntab2(nn,ni,nj,h,hx,hy,hygx,hxgy,uygx,uxgy,uxy)
     int **nn,ni,nj;
     float *h,*hx,*hy,*hygx,*hxgy,*uygx,*uxgy,*uxy;
{
  int i,j;
  float sum=0.0,p,*sumi,*sumj;
  
  sumi=vector(1,ni);
  sumj=vector(1,nj);
  for (i=1;i<=ni;i++) {
    sumi[i]=0.0;
    for (j=1;j<=nj;j++) {
      sumi[i] += nn[i][j];
      sum += nn[i][j];
    }
  }
  for (j=1;j<=nj;j++) {
    sumj[j]=0.0;
    for (i=1;i<=ni;i++)
      sumj[j] += nn[i][j];
  }
  *hx=0.0;
  for (i=1;i<=ni;i++)
    if (sumi[i]) {
      p=sumi[i]/sum;
      *hx -= p*log(p);
    }
  *hy=0.0;
  for (j=1;j<=nj;j++)
    if (sumj[j]) {
      p=sumj[j]/sum;
      *hy -= p*log(p);
    }
  *h=0.0;
  for (i=1;i<=ni;i++)
    for (j=1;j<=nj;j++)
      if (nn[i][j]) {
	p=nn[i][j]/sum;
	*h -= p*log(p);
      }
  *hygx=(*h)-(*hx);
  *hxgy=(*h)-(*hy);
  *uygx=(*hy-*hygx)/(*hy+TINY);
  *uxgy=(*hx-*hxgy)/(*hx+TINY);
  *uxy=2.0*(*hx+*hy-*h)/(*hx+*hy+TINY);
  free_vector(sumj,1,nj);
  free_vector(sumi,1,ni);
}
#undef TINY
#undef NRANSI
/*************************************---*************************************/
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
/**************************************-**************************************/
/*                                                                           */
/*                                   FOURN                                   */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void fourn(data,nn,ndim,isign)
     float data[];
     unsigned long nn[];
     int ndim,isign;
{
  int idim;
  unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
  unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
  float tempi,tempr;
  double theta,wi,wpi,wpr,wr,wtemp;
  
  for (ntot=1,idim=1;idim<=ndim;idim++)
    ntot *= nn[idim];
  nprev=1;
  for (idim=ndim;idim>=1;idim--) {
    n=nn[idim];
    nrem=ntot/(n*nprev);
    ip1=nprev << 1;
    ip2=ip1*n;
    ip3=ip2*nrem;
    i2rev=1;
    for (i2=1;i2<=ip2;i2+=ip1) {
      if (i2 < i2rev) {
	for (i1=i2;i1<=i2+ip1-2;i1+=2) {
	  for (i3=i1;i3<=ip3;i3+=ip2) {
	    i3rev=i2rev+i3-i2;
	    SWAP(data[i3],data[i3rev]);
	    SWAP(data[i3+1],data[i3rev+1]);
	  }
	}
      }
      ibit=ip2 >> 1;
      while (ibit >= ip1 && i2rev > ibit) {
	i2rev -= ibit;
	ibit >>= 1;
      }
      i2rev += ibit;
    }
    ifp1=ip1;
    while (ifp1 < ip2) {
      ifp2=ifp1 << 1;
      theta=isign*6.28318530717959/(ifp2/ip1);
      wtemp=sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (i3=1;i3<=ifp1;i3+=ip1) {
	for (i1=i3;i1<=i3+ip1-2;i1+=2) {
	  for (i2=i1;i2<=ip3;i2+=ifp2) {
	    k1=i2;
	    k2=k1+ifp1;
	    tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
	    tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
	    data[k2]=data[k1]-tempr;
	    data[k2+1]=data[k1+1]-tempi;
	    data[k1] += tempr;
	    data[k1+1] += tempi;
	  }
	}
	wr=(wtemp=wr)*wpr-wi*wpi+wr;
	wi=wi*wpr+wtemp*wpi+wi;
      }
      ifp1=ifp2;
    }
    nprev *= n;
  }
}
#undef SWAP
/*************************************---*************************************/
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
/**************************************-**************************************/
/*                                                                           */
/*                                  MY_FOURN                                 */
/*                                                                           */
/*  Array index starts from zero.                                            */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void my_fourn(data,nn,ndim,isign)
     float data[];   // Index from 0..
     unsigned long nn[];
     int ndim,isign;
{
  int idim;
  unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
  unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
  float tempi,tempr;
  double theta,wi,wpi,wpr,wr,wtemp;
  
  for (ntot=1,idim=1;idim<=ndim;idim++)
    ntot *= nn[idim];
  nprev=1;
  for (idim=ndim;idim>=1;idim--) {
    n=nn[idim];
    nrem=ntot/(n*nprev);
    ip1=nprev << 1;
    ip2=ip1*n;
    ip3=ip2*nrem;
    i2rev = 1;
    for (i2=1;i2<=ip2;i2+=ip1) {
      if (i2 < i2rev) {
	for (i1=i2;i1<=i2+ip1-2;i1+=2){
	  for (i3=i1;i3<=ip3;i3+=ip2){
	    i3rev = i2rev + i3 - i2;

	    //SWAP(data[i3],data[i3rev]);
	    //SWAP(data[i3+1],data[i3rev+1]);

	    SWAP(data[i3-1],data[i3rev-1]);
	    SWAP(data[i3],data[i3rev]);
	  }
	}
      }
      ibit = ip2 >> 1;
      while (ibit >= ip1 && i2rev > ibit) {
	i2rev -= ibit;
	ibit >>= 1;
      }
      i2rev += ibit;
    }
    ifp1=ip1;
    while (ifp1 < ip2) {
      ifp2 = ifp1 << 1;
      theta = isign*6.28318530717959/(ifp2/ip1);
      wtemp = sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi = sin(theta);
      wr = 1.0;
      wi = 0.0;
      for (i3=1;i3<=ifp1;i3+=ip1){
	for (i1=i3;i1<=i3+ip1-2;i1+=2){
	  for (i2=i1;i2<=ip3;i2+=ifp2){
	    k1 = i2;
	    k2 = k1 + ifp1;

	    /*
	    tempr = (float)wr * data[k2]  - (float)wi * data[k2+1];
	    tempi = (float)wr * data[k2+1]+ (float)wi * data[k2];
	    data[k2]   = data[k1]   - tempr;
	    data[k2+1] = data[k1+1] - tempi;
	    data[k1] += tempr;
	    data[k1+1] += tempi;
	    */

	    tempr = (float)wr * data[k2-1] - (float)wi * data[k2];
	    tempi = (float)wr * data[k2]   + (float)wi * data[k2-1];
	    data[k2-1]  = data[k1-1] - tempr;
	    data[k2]    = data[k1]   - tempi;
	    data[k1-1] += tempr;
	    data[k1]   += tempi;
	  }
	}
	wtemp = wr;
	wr = wtemp*wpr -    wi*wpi + wr;
	wi =    wi*wpr + wtemp*wpi + wi;
      }
      ifp1 = ifp2;
    }
    nprev *= n;
  }
}
#undef SWAP
/**************************************-**************************************/
/*                                                                           */
/*                                   RLFT3                                   */
/*                                                                           */
/*  Given a three-dimensional real array data[1..nn1][1..nn2][1..nn3]        */
/*  (where nn1=1 for the case of a logically two-dimensional array), this    */
/*  routine returns (for isign=1) the complex fast Fourier transform as two  */
/*  complex arrays: On output, "data" contains the zero and positive         */
/*  frequency values of the third frequency component, while                 */
/*  speq[1..nn1][1..2*nn2] contains the Nyquist critical frequency values    */
/*  of the third frequency component.  First (and second) frequency          */
/*  components are stored for zero, positive, and negative frequencies, in   */
/*  standard wrap-around order.  See text for description of how complex     */
/*  values are arranged.  For isign=-1, the inverse transform (times         */
/*  nn1*nn2*nn3/2 as a constant multiplicative factor) is performed, with    */
/*  output "data" (viewed as a real array) deriving from input "data"        */
/*  (viewed as complex) and "speq".  The dimensions nn1, nn2, nn3 must       */
/*  always be integer powers of 2.                                           */
/*                                                                           */
/*  What is returned may be thought of as a complex array 'SPEC' with the    */
/*  following layout:                                                        */
/*                                                                           */
/*    SPEC[1..nn1][1..nn2][1..nn3/2]                                         */
/*                                                                           */
/*  In the 1st two dimensions, the spectrum is packed in wrap-around order:  */
/*                                                                           */
/*   1     f = 0                                                             */
/*   2     f = 1/(N*D) D= 1 pix    [smallest positive freq]                  */
/*   ...	                                                             */
/*  N/2    f = (N/2 - 1) / (N*D)                                             */
/*  N/2+1  f = +- 1/(2*D)          [combintation cutoff frequency]           */
/*  N/2+2  f = -(N/2 - 1) / (N*D)  [most negative freq, other than cutoff]   */
/*  ...                                                                      */
/*   N     f = -1/(N*D)            [smallest negative freq]                  */
/*                                                                           */
/*  Re(SPEC[i1][i2][i3]) = data[i1][i2][2*i3-1]                              */
/*  Im(SPEC[i1][i2][i3]) = data[i1][i2][2*i3  ]                              */
/*                                                                           */
/*  For values in the "plane" SPEC[1..nn1][1..nn2][nn3/2+1], which are all   */
/*  values at the cutoff freq in the 3rd dimension, these are stored in      */
/*  'speq' as:                                                               */
/*                                                                           */
/*    Re(SPEC[i1][i2][nn3/2+1]) = speq[i1][2*i2-1]                           */
/*    Im(SPEC[i1][i2][nn3/2+1]) = speq[i1][2*i2  ]                           */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void rlft3(data,speq,nn1,nn2,nn3,isign)
     float ***data,**speq;
     unsigned long nn1,nn2,nn3;
     int isign;
{
  void fourn();
  unsigned long i1,i2,i3,j1,j2,j3,nn[4],ii3;
  double theta,wi,wpi,wpr,wr,wtemp;
  float c1,c2,h1r,h1i,h2r,h2i;

  if (1+&data[nn1][nn2][nn3]-&data[1][1][1] != nn1*nn2*nn3)
    nrerror("rlft3: problem with dimensions or contiguity of data array\n");
  c1=0.5;
  c2 = -0.5*isign;
  theta=isign*(6.28318530717959/nn3);
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  nn[1]=nn1;
  nn[2]=nn2;
  nn[3]=nn3 >> 1;
  if (isign == 1) {
    fourn(&data[1][1][1]-1,nn,3,isign);
    for (i1=1;i1<=nn1;i1++)
      for (i2=1,j2=0;i2<=nn2;i2++) {
	speq[i1][++j2]=data[i1][i2][1];
	speq[i1][++j2]=data[i1][i2][2];
      }
  }
  for (i1=1;i1<=nn1;i1++) {
    j1=(i1 != 1 ? nn1-i1+2 : 1);
    wr=1.0;
    wi=0.0;
    for (ii3=1,i3=1;i3<=(nn3>>2)+1;i3++,ii3+=2) {
      for (i2=1;i2<=nn2;i2++) {
	if (i3 == 1) {
	  j2=(i2 != 1 ? ((nn2-i2)<<1)+3 : 1);
	  h1r=c1*(data[i1][i2][1]+speq[j1][j2]);
	  h1i=c1*(data[i1][i2][2]-speq[j1][j2+1]);
	  h2i=c2*(data[i1][i2][1]-speq[j1][j2]);
	  h2r= -c2*(data[i1][i2][2]+speq[j1][j2+1]);
	  data[i1][i2][1]=h1r+h2r;
	  data[i1][i2][2]=h1i+h2i;
	  speq[j1][j2]=h1r-h2r;
	  speq[j1][j2+1]=h2i-h1i;
	} else {
	  j2=(i2 != 1 ? nn2-i2+2 : 1);
	  j3=nn3+3-(i3<<1);
	  h1r=c1*(data[i1][i2][ii3]+data[j1][j2][j3]);
	  h1i=c1*(data[i1][i2][ii3+1]-data[j1][j2][j3+1]);
	  h2i=c2*(data[i1][i2][ii3]-data[j1][j2][j3]);
	  h2r= -c2*(data[i1][i2][ii3+1]+data[j1][j2][j3+1]);
	  data[i1][i2][ii3]=h1r+wr*h2r-wi*h2i;
	  data[i1][i2][ii3+1]=h1i+wr*h2i+wi*h2r;
	  data[j1][j2][j3]=h1r-wr*h2r+wi*h2i;
	  data[j1][j2][j3+1]= -h1i+wr*h2i+wi*h2r;
	}
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
  }

  if (isign == -1)
    fourn(&data[1][1][1]-1,nn,3,isign);

}
/**************************************-**************************************/
/*                                                                           */
/*                            NRTEMP_GET_3D_FARRAY                           */
/*                                                                           */
/*  Can be used as 'get_2d_pointer_farray' by setting d3 to 0.               */
/*                                                                           */
/*****************************************************************************/
float ***nrtemp_get_3d_farray(d1,d2,d3)
     int d1,d2,d3;
{
  int i1,i2;
  float ***data;

  //
  // ********* WYETH REMOVE THIS ROUTINE
  // ********* WYETH REMOVE THIS ROUTINE
  // ********* WYETH REMOVE THIS ROUTINE
  //
  
  data = (float ***)myalloc(d1*sizeof(float **));
  if (d2 >= 0)
    for(i1=0;i1<d1;i1++){
      data[i1] = (float **)myalloc(d2*sizeof(float *));
      if (d3 >= 0){
	for(i2=0;i2<d2;i2++)
	  data[i1][i2] = (float *)myalloc(d3*sizeof(float)); // NULL if d3=0
      }
    }
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                                  MY_RLFT3                                 */
/*                                                                           */
/*  Here I modified NumRec's 'rlft3' routine to eliminate teh need for       */
/*  array pointer arithmetic, so this can be more easily converted to Java.  */
/*                                                                           */
/*****************************************************************************/
void my_rlft3(tdata,tspeq,nn1,nn2,nn3,isign)
     float ***tdata,**tspeq;
     long nn1,nn2,nn3;  // Changed from 'unsigned long'
     int isign;
{
  int i,j,k,l;
  float *t1d;

  void my_fourn();
  //void fourn();
  long i1,i2,i3,j1,j2,j3,ii3;
  //long nn[4];
  unsigned long nn[4];
  double theta,wi,wpi,wpr,wr,wtemp;
  float c1,c2,h1r,h1i,h2r,h2i;
  float ***data,**speq;

  //printf("  MY_RLFT3\n");

  // COPY TO ZERO-INDEX ARRAY
  //data = f3tensor(1,nn1,1,nn2,1,nn3);
  data = nrtemp_get_3d_farray(nn1,nn2,nn3);
  for(i=0;i<nn1;i++){
    for(j=0;j<nn2;j++){
      for(k=0;k<nn3;k++){
	data[i][j][k] = tdata[i+1][j+1][k+1];
      }
    }
  }

  // COPY TO ZERO-INDEX ARRAY
  speq = (float **)myalloc(nn1*sizeof(float *));
  for(i=0;i<nn1;i++)
    speq[i] = (float *)myalloc((2*nn2)*sizeof(float));
  for(i=0;i<nn1;i++){
    for(j=0;j<(2*nn2);j++){
      speq[i][j] = tspeq[i+1][j+1];
    }
  }

  
  c1 = 0.5;
  c2 = -0.5*isign;
  theta = isign*(6.28318530717959/nn3);
  wtemp = sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi = sin(theta);
  nn[1] = nn1;
  nn[2] = nn2;
  //nn[3] = nn3 >> 1;
  nn[3] = nn3 / 2; //>> 1;

  t1d = get_farray(nn1*nn2*nn3);  // WYETH added

  if (isign == 1){
    l = 0; // Copy to a 1D array
    for(i=0;i<nn1;i++){
      for(j=0;j<nn2;j++){
	for(k=0;k<nn3;k++){
	  t1d[l] = data[i][j][k]; // NEW
	  l += 1;
	}
      }
    }
    my_fourn(t1d,nn,3,isign);
    l = 0; // Copy back to 3D array
    for(i=0;i<nn1;i++){
      for(j=0;j<nn2;j++){
	for(k=0;k<nn3;k++){
	  data[i][j][k] = t1d[l]; // NEW
	  l += 1;
	}
      }
    }

    for (i1=0;i1<nn1;i1++)
      for (i2=0,j2=0;i2<nn2;i2++){
	speq[i1][j2] = data[i1][i2][0]; // NEW
	j2 += 1;
	speq[i1][j2] = data[i1][i2][1]; // NEW
	j2 += 1;
      }
  }

  for (i1=0;i1<nn1;i1++){  // NEW
    if (i1 != 0) // NEW
      j1 = nn1 - i1;
    else
      j1 = 0;

    wr = 1.0;
    wi = 0.0;
    for (ii3=0,i3=1;i3<=(nn3>>2)+1;i3++,ii3+=2){ // NEW
      for (i2=0;i2<nn2;i2++){  // NEW
	if (i3 == 1) {
	  if (i2 != 0)  // NEW
	    j2 = ((nn2-(i2+1))<<1) + 2;
	  else
	    j2 = 0;

	  // NEW	  
	  h1r =  c1*(data[i1][i2][0] + speq[j1][j2]); 
	  h1i =  c1*(data[i1][i2][1] - speq[j1][j2+1]);
	  h2i =  c2*(data[i1][i2][0] - speq[j1][j2]);
	  h2r = -c2*(data[i1][i2][1] + speq[j1][j2+1]);
	  data[i1][i2][0] = h1r + h2r;
	  data[i1][i2][1] = h1i + h2i;
	  speq[j1][j2]    = h1r - h2r;
	  speq[j1][j2+1]  = h2i - h1i;

	}else{

	  if (i2 != 0) // NEW
	    j2 = nn2 - i2;
	  else
	    j2 = 0;

	  j3 = nn3 + 3 - (i3<<1) - 1;  // NEW
	  h1r =  c1*(data[i1][i2][ii3]   + data[j1][j2][j3]);
	  h1i =  c1*(data[i1][i2][ii3+1] - data[j1][j2][j3+1]);
	  h2i =  c2*(data[i1][i2][ii3]   - data[j1][j2][j3]);
	  h2r = -c2*(data[i1][i2][ii3+1] + data[j1][j2][j3+1]);
	  data[i1][i2][ii3]   =  h1r + wr*h2r - wi*h2i;
	  data[i1][i2][ii3+1] =  h1i + wr*h2i + wi*h2r;
	  data[j1][j2][j3]    =  h1r - wr*h2r + wi*h2i;
	  data[j1][j2][j3+1]  = -h1i + wr*h2i + wi*h2r;
	}
      }
      wtemp = wr;
      wr = wr*wpr -    wi*wpi + wr;
      wi = wi*wpr + wtemp*wpi + wi;
    }
  }
  if (isign == -1){
    l = 0; // Copy into a 1D array
    for(i=0;i<nn1;i++){
      for(j=0;j<nn2;j++){
	for(k=0;k<nn3;k++){
	  t1d[l] = data[i][j][k];
	  l += 1;
	}
      }
    }
    my_fourn(t1d,nn,3,isign);
    l = 0; // Copy back to 3D array
    for(i=0;i<nn1;i++){
      for(j=0;j<nn2;j++){
	for(k=0;k<nn3;k++){
	  data[i][j][k] = t1d[l];
	  l += 1;
	}
      }
    }
  }
  myfree(t1d);  // WYETH Added, free the temporary 1D array

  // COPY BACK FROM ZERO-INDEX ARRAY
  for(i=0;i<nn1;i++){
    for(j=0;j<nn2;j++){
      for(k=0;k<nn3;k++){
	tdata[i+1][j+1][k+1] = data[i][j][k]; // NEW
      }
    }
  }
  for(i=0;i<nn1;i++){
    for(j=0;j<(2*nn2);j++){
      tspeq[i+1][j+1] = speq[i][j];
    }
  }

}
/*************************************---*************************************/
#define NRANSI
/**************************************-**************************************/
/*                                                                           */
/*                                    FIT                                    */
/*                                                                           */
/*  Given a set of data points x[1..ndata], y[1..ndata] with individual      */
/*  standard deviations sig[1..ndata], fit them to a straight line y=a+bx    */
/*  by minimizing Chi^2.  Returned are a,b and their respective probably     */
/*  uncertainties siga and sigb, the chi-square chi2, and the                */
/*  goodness-of-fit probability q (that the fit would hvae chi^2 this        */
/*  large or larger). If mwt=0 on input, then the standard deviations are    */
/*  assumed to be unavailable: q is returned as 1.0 and the normalization    */
/*  of chi^2 is to unit standard deviation on all points.                    */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void fit(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q)
     float x[],y[];
     int ndata;
     float sig[];
     int mwt;
     float *a,*b,*siga,*sigb,*chi2,*q;
{
  float gammq();
  int i;
  float wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;
  
  *b=0.0;
  if (mwt) {
    ss=0.0;
    for (i=1;i<=ndata;i++) {
      //wt=1.0/SQR(sig[i]);
      wt=1.0/(sig[i]*sig[i]);
      ss += wt;
      sx += x[i]*wt;
      sy += y[i]*wt;
    }
  } else {
    for (i=1;i<=ndata;i++) {
      sx += x[i];
      sy += y[i];
    }
    ss=ndata;
  }
  sxoss=sx/ss;
  if (mwt) {
    for (i=1;i<=ndata;i++) {
      t=(x[i]-sxoss)/sig[i];
      st2 += t*t;
      *b += t*y[i]/sig[i];
    }
  } else {
    for (i=1;i<=ndata;i++) {
      t=x[i]-sxoss;
      st2 += t*t;
      *b += t*y[i];
    }
  }
  *b /= st2;
  *a=(sy-sx*(*b))/ss;
  *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
  *sigb=sqrt(1.0/st2);
  *chi2=0.0;
  if (mwt == 0) {
    for (i=1;i<=ndata;i++)
      *chi2 += (y[i]-(*a)-(*b)*x[i]) * (y[i]-(*a)-(*b)*x[i]);
    //*chi2 += SQR(y[i]-(*a)-(*b)*x[i]);
    *q=1.0;
    sigdat=sqrt((*chi2)/(ndata-2));
    *siga *= sigdat;
    *sigb *= sigdat;
  } else {
    for (i=1;i<=ndata;i++)
      *chi2 += ((y[i]-(*a)-(*b)*x[i])/sig[i]) * ((y[i]-(*a)-(*b)*x[i])/sig[i]);
    //*chi2 += SQR((y[i]-(*a)-(*b)*x[i])/sig[i]);
    *q=gammq(0.5*(ndata-2),0.5*(*chi2));
  }
}
#undef NRANSI
/*************************************---*************************************/
#define NRANSI
/**************************************-**************************************/
/*                                                                           */
/*                                  PYTHAG                                   */
/*                                                                           */
/*  Computes sqrt(a^2 + b^2) without destructive underflow or overflow.      */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
float pythag(a,b)
     float a,b;
{
  float absa,absb;

  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+(absb/absa)*(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
  //if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  //else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}
#undef NRANSI
/*************************************---*************************************/
#define NRANSI
/**************************************-**************************************/
/*                                                                           */
/*                                  TRIDAG                                   */
/*                                                                           */
/*  Solves for a vector u[1..n] the tridiagonal linear set given by          */
/*  equation (2.4.1).  a[1..n], b[1..n], c[1..n], and r[1..n] are input      */
/*  vectors and are not modified.                                            */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void tridag(float a[], float b[], float c[], float r[], float u[],
	unsigned long n)
{
  unsigned long j;
  float bet,*gam;

  gam=vector(1,n);
  if (b[1] == 0.0) nrerror("Error 1 in tridag");
  u[1]=r[1]/(bet=b[1]);
  for (j=2;j<=n;j++) {
    gam[j]=c[j-1]/bet;
    bet=b[j]-a[j]*gam[j];
    if (bet == 0.0) nrerror("Error 2 in tridag");
    u[j]=(r[j]-a[j]*u[j-1])/bet;
  }
  for (j=(n-1);j>=1;j--)
    u[j] -= gam[j+1]*u[j+1];
  free_vector(gam,1,n);
}
#undef NRANSI
/**************************************-**************************************/
/*                                                                           */
/*                                  EIGSRT                                   */
/*                                                                           */
/*  Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n] as output   */
/*  from "jacobi" or "tqli", this routine sorts the eigenvalues into         */
/*  descending order, and rearranges the columns of "v" correspondingly.     */
/*  The method is straight insertion.                                        */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void eigsrt(d,v,n)
     float d[],**v;
     int n;
{
  int k,j,i;
  float p;
  
  for (i=1;i<n;i++) {
    p=d[k=i];
    for (j=i+1;j<=n;j++)
      if (d[j] >= p) p=d[k=j];
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
      for (j=1;j<=n;j++) {
	p=v[j][i];
	v[j][i]=v[j][k];
	v[j][k]=p;
      }
    }
  }
}
/*************************************---*************************************/
#define NRANSI
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);
/**************************************-**************************************/
/*                                                                           */
/*                                  JACOBI                                   */
/*                                                                           */
/*  Computes all eigenvalues and eigenvectors of a real symmetric matrix     */
/*  a[1..n][1..n].  On output, elements of "a" above the diagonal are        */
/*  destroyed.  d[1..n] returns the eigenvalues of "a".  v[1..n][1..n] is    */
/*  a matrix whose columns contain, on output, the normalized eigenvectors   */
/*  of "a".  "nrot" returns the number of Jacobi rotations that were         */
/*  required.                                                                */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void jacobi(a,n,d,v,nrot)
     float **a;
     int n;
     float d[],**v;
     int *nrot;
{
  int j,iq,ip,i;
  float tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
  
  b=vector(1,n);
  z=vector(1,n);
  for (ip=1;ip<=n;ip++) {
    for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
  for (ip=1;ip<=n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++)
	sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0) {
      free_vector(z,1,n);
      free_vector(b,1,n);
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++) {
	g=100.0*fabs(a[ip][iq]);
	if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
	    && (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
	  a[ip][iq]=0.0;
	else if (fabs(a[ip][iq]) > tresh) {
	  h=d[iq]-d[ip];
	  if ((float)(fabs(h)+g) == (float)fabs(h))
	    t=(a[ip][iq])/h;
	  else {
	    theta=0.5*h/(a[ip][iq]);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq]=0.0;
	  for (j=1;j<=ip-1;j++) {
	    ROTATE(a,j,ip,j,iq)
	    }
	  for (j=ip+1;j<=iq-1;j++) {
	    ROTATE(a,ip,j,j,iq)
	    }
	  for (j=iq+1;j<=n;j++) {
	    ROTATE(a,ip,j,iq,j)
	    }
	  for (j=1;j<=n;j++) {
	    ROTATE(v,j,ip,j,iq)
	    }
	  ++(*nrot);
	}
      }
    }
    for (ip=1;ip<=n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  nrerror("Too many iterations in routine jacobi");
}
#undef ROTATE
#undef NRANSI
/**************************************-**************************************/
/*                                                                           */
/*                                   TRED2                                   */
/*                                                                           */
/*  Householder reduction of a real, symmetric matrix a[1..n][1..n].  On     */
/*  output, "a" is replaced by the orthogonal matrix Q effecting the         */
/*  transformation.  d[1..n] returns the diagonal elements of the            */
/*  tridiagonal matrix, and e[1..n] the off-diagonal elements, with          */
/*  e[1]=0.  Several statments, as noted in comments, can be omitted if      */
/*  only eigenvalues are to be found, in which case "a" contains no useful   */
/*  information on output.  Otherwise they are to be included.               */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void tred2(a,n,d,e)
     float **a;
     int n;
     float d[],e[];
{
  int l,k,j,i;
  float scale,hh,h,g,f;
  
  for (i=n;i>=2;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 1) {
      for (k=1;k<=l;k++)
	scale += fabs(a[i][k]);
      if (scale == 0.0)
	e[i]=a[i][l];
      else {
	for (k=1;k<=l;k++) {
	  a[i][k] /= scale;
	  h += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
	e[i]=scale*g;
	h -= f*g;
	a[i][l]=f-g;
	f=0.0;
	for (j=1;j<=l;j++) {
	  a[j][i]=a[i][j]/h;
	  g=0.0;
	  for (k=1;k<=j;k++)
	    g += a[j][k]*a[i][k];
	  for (k=j+1;k<=l;k++)
	    g += a[k][j]*a[i][k];
	  e[j]=g/h;
	  f += e[j]*a[i][j];
	}
	hh=f/(h+h);
	for (j=1;j<=l;j++) {
	  f=a[i][j];
	  e[j]=g=e[j]-hh*f;
	  for (k=1;k<=j;k++)
	    a[j][k] -= (f*e[k]+g*a[i][k]);
	}
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  d[1]=0.0;
  e[1]=0.0;
  /* Contents of this loop can be omitted if eigenvectors not
     wanted except for statement d[i]=a[i][i]; */
  for (i=1;i<=n;i++) {
    l=i-1;
    if (d[i]) {
      for (j=1;j<=l;j++) {
	g=0.0;
	for (k=1;k<=l;k++)
	  g += a[i][k]*a[k][j];
	for (k=1;k<=l;k++)
	  a[k][j] -= g*a[k][i];
      }
    }
    d[i]=a[i][i];
    a[i][i]=1.0;
    for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
  }
}
/*************************************---*************************************/
#define NRANSI
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
/**************************************-**************************************/
/*                                                                           */
/*                                   TQLI                                    */
/*                                                                           */
/*  QL algorithm with implicit shifts, to determine the eigenvalues and      */
/*  eigenvectors of a real, symmetric, tridiagonal matrix, or of a real,     */
/*  symmetric matrix previously reduced by "tred2".  On input, d[1..n]       */
/*  contains the diagonal elements of the tridiagonal matrix.  On output,    */
/*  it returns the eigenvalues.  The vector e[1..n] inputs the subdiagonal   */
/*  elements of the tridiagonal matrix, with e[1] arbitrary.  On output      */
/*  "e" is destroyed.  When finding only the eigenvalues, several lines may  */
/*  be omitted, as noted in the comments.  If the eigenvectors of a          */
/*  tridiagonal matrix are desired, the matrix z[1..n][1..n] is input as     */
/*  the identity matrix.  If the eigenvectors of a matrix that has been      */
/*  reduced by "tred2" are required, then "z" is input as the matrix output  */
/*  by "tred2".  In either case, the kth column of "z" returns the           */
/*  normalized eigenvector corresponding to d[k].                            */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void tqli(d,e,n,z)
     float d[],e[];
     int n;
     float **z;
{
  float pythag();
  int m,l,iter,i,k;
  float s,r,p,g,f,dd,c,b;
  
  for (i=2;i<=n;i++) e[i-1]=e[i];
  e[n]=0.0;
  for (l=1;l<=n;l++) {
    iter=0;
    do {
      for (m=l;m<=n-1;m++) {
	dd=fabs(d[m])+fabs(d[m+1]);
	if ((float)(fabs(e[m])+dd) == dd) break;
      }
      if (m != l) {
	if (iter++ == 30) nrerror("Too many iterations in tqli");
	g=(d[l+1]-d[l])/(2.0*e[l]);
	r=pythag(g,1.0);
	g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--) {
	  f=s*e[i];
	  b=c*e[i];
	  e[i+1]=(r=pythag(f,g));
	  if (r == 0.0) {
	    d[i+1] -= p;
	    e[m]=0.0;
	    break;
	  }
	  s=f/r;
	  c=g/r;
	  g=d[i+1]-p;
	  r=(d[i]-g)*s+2.0*c*b;
	  d[i+1]=g+(p=s*r);
	  g=c*r-b;
	  for (k=1;k<=n;k++) {
	    f=z[k][i+1];
	    z[k][i+1]=s*z[k][i]+c*f;
	    z[k][i]=c*z[k][i]-s*f;
	  }
	}
	if (r == 0.0 && i >= l) continue;
	d[l] -= p;
	e[l]=g;
	e[m]=0.0;
      }
    } while (m != l);
  }
}
#undef SIGN
#undef NRANSI
/*************************************---*************************************/
#define NRANSI
/**************************************-**************************************/
/*                                                                           */
/*                                  MRQCOF                                   */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs)
     float x[],y[],sig[];
     int ndata;
     float a[];
     int ia[],ma;
     float **alpha,beta[],*chisq;
     void (*funcs)();  /*** float, float [], float *, float [], int)) ***/
{
  int i,j,k,l,m,mfit=0;
  float ymod,wt,sig2i,dy,*dyda;
  
  dyda=vector(1,ma);
  for (j=1;j<=ma;j++)
    if (ia[j]) mfit++;
  for (j=1;j<=mfit;j++){
    for (k=1;k<=j;k++) alpha[j][k]=0.0;
    beta[j]=0.0;
  }
  *chisq=0.0;
  for (i=1;i<=ndata;i++){
    (*funcs)(x[i],a,&ymod,dyda,ma);
    sig2i = 1.0/(sig[i]*sig[i]);
    dy = y[i]-ymod;
    for (j=0,l=1;l<=ma;l++){
      if (ia[l]){
	wt = dyda[l]*sig2i;
	for (j++,k=0,m=1;m<=l;m++)
	  if (ia[m]) alpha[j][++k] += wt*dyda[m];
	beta[j] += dy*wt;
      }
    }
    *chisq += dy*dy*sig2i;
  }

  for (j=2;j<=mfit;j++)
    for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
  free_vector(dyda,1,ma);
}
#undef NRANSI
/*************************************---*************************************/
#define NRANSI
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
/**************************************-**************************************/
/*                                                                           */
/*                                  GAUSSJ                                   */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void gaussj(a,n,b,m)
     float **a;
     int n;
     float **b;
     int m;
{
  int *indxc,*indxr,*ipiv;
  int i,icol,irow,j,k,l,ll;
  float big,dum,pivinv,temp;
  
  indxc=ivector(1,n);
  indxr=ivector(1,n);
  ipiv=ivector(1,n);
  for (j=1;j<=n;j++) ipiv[j]=0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if (ipiv[j] != 1)
	for (k=1;k<=n;k++) {
	  if (ipiv[k] == 0) {
	    if (fabs(a[j][k]) >= big) {
	      big=fabs(a[j][k]);
	      irow=j;
	      icol=k;
	    }
	  } else if (ipiv[k] > 1) nrerror("gaussj: Singular Matrix-1");
	}
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
	for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
	}
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix-2");
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=1;l<=n;l++) a[icol][l] *= pivinv;
    for (l=1;l<=m;l++) b[icol][l] *= pivinv;
    for (ll=1;ll<=n;ll++)
      if (ll != icol) {
	dum=a[ll][icol];
	a[ll][icol]=0.0;
	for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
	for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
      }
  }
  for (l=n;l>=1;l--) {
    if (indxr[l] != indxc[l])
      for (k=1;k<=n;k++)
	SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  }
  free_ivector(ipiv,1,n);
  free_ivector(indxr,1,n);
  free_ivector(indxc,1,n);
}
#undef NRANSI
#undef SWAP
/*************************************---*************************************/
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
/**************************************-**************************************/
/*                                                                           */
/*                                  COVSRT                                   */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void covsrt(covar,ma,ia,mfit)
     float **covar;
     int ma;
     int ia[];
     int mfit;
{
  int i,j,k;
  float swap;
  
  for (i=mfit+1;i<=ma;i++)
    for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
  k=mfit;
  for (j=ma;j>=1;j--) {
    if (ia[j]) {
      for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
	for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
	  k--;
    }
  }
}
#undef SWAP
/*************************************---*************************************/
#define NRANSI
/**************************************-**************************************/
/*                                                                           */
/*                                   MRQMIN                                  */
/*                                                                           */
/*  Levenberg-Marquardt method, attempting to reduce the value chisq of a    */
/*  fit between a set of data points x[1..ndata],y[1..ndata] with            */
/*  individual standard deviations sig[1..ndata], and a nonlinear function   */
/*  dependent on ma coefficients a[1..ma].  The input array ia[1..ma]        */
/*  indicates by nonzero entries those components of a that should be        */
/*  fitted for, and by zero entries those components that should be held     */
/*  fixed at their input values.  The program returns current best-fit       */
/*  values for the parameters a[1..ma], and chisq.  The arrays               */
/*  covar[1..ma][1..ma], alpha[1..ma][1..ma] are used as working space       */
/*  during most iterations.  Supply a routine funcs(x,a,yfit,dyda,ma) that   */
/*  evaluates the fitting function yfit, and its derivatives dyda[1..ma]     */
/*  with respoect to the fitting parameters a at x.  On the first call,      */
/*  provide an initial guess for the parameters a, and set alamda < 0 for    */
/*  initialization (which then sets alamda=0.001).  If a step succeeds       */
/*  chisq becomes smaller and alamda decreases by a factor of 10.  If a      */
/*  step fails alamda grows by a factor of 10.  You must call this routine   */
/*  repeatedly until convergence is achieved.  Then, make one final call     */
/*  with alamda=0, so that covar[1..ma][1..ma] returns the covariance        */
/*  matrix, and alpha the curvature matrix.  (Parameters held fixed will     */
/*  return zero covariances.)                                                */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,chisq,funcs,alamda)
     float x[],y[],sig[];
     int ndata;
     float a[];
     int ia[],ma;
     float **covar,**alpha,*chisq;
     void (*funcs)(); /*** float, float [], float *, float [], int),  ***/
     float *alamda;
{
  void covsrt(); /***float **covar, int ma, int ia[], int mfit);***/
  void gaussj(); /***float **a, int n, float **b, int m);***/
  void mrqcof(); /***float x[], float y[], float sig[], int ndata, float a[],
		   int ia[], int ma, float **alpha, float beta[], float *chisq,
		   void (*funcs)(float, float [], float *, float [], int));***/
  int j,k,l,m;
  static int mfit;
  static float ochisq,*atry,*beta,*da,**oneda;

  if (*alamda < 0.0) {
    atry=vector(1,ma);
    beta=vector(1,ma);
    da=vector(1,ma);
    for (mfit=0,j=1;j<=ma;j++)
      if (ia[j]) mfit++;
    oneda=matrix(1,mfit,1,1);
    *alamda=0.001;
    mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs);
    ochisq=(*chisq);
    for (j=1;j<=ma;j++) atry[j]=a[j];
  }
  for (j=0,l=1;l<=ma;l++) {
    if (ia[l]) {
      for (j++,k=0,m=1;m<=ma;m++) {
	if (ia[m]) {
	  k++;
	  covar[j][k]=alpha[j][k];
	}
      }
      covar[j][j]=alpha[j][j]*(1.0+(*alamda));
      oneda[j][1]=beta[j];
    }
  }
  gaussj(covar,mfit,oneda,1);
  for (j=1;j<=mfit;j++) da[j]=oneda[j][1];
  if (*alamda == 0.0) {
    covsrt(covar,ma,ia,mfit);
    free_matrix(oneda,1,mfit,1,1);
    free_vector(da,1,ma);
    free_vector(beta,1,ma);
    free_vector(atry,1,ma);
    return;
  }
  for (j=0,l=1;l<=ma;l++)
    if (ia[l]) atry[l]=a[l]+da[++j];
  mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs);
  if (*chisq < ochisq) {
    *alamda *= 0.1;
    ochisq=(*chisq);
    for (j=0,l=1;l<=ma;l++) {
      if (ia[l]) {
	for (j++,k=0,m=1;m<=ma;m++) {
	  if (ia[m]) {
	    k++;
	    alpha[j][k]=covar[j][k];
	  }
	}
	beta[j]=da[j];
	a[l]=atry[l];
      }
    }
  }else{
    *alamda *= 10.0;
    *chisq=ochisq;
  }
}
#undef NRANSI
/**************************************-**************************************/
/*                                                                           */
/*                                  FGAUSS                                   */
/*                                                                           */
/*  y(x;a) is the sum of na/3 gaussians.  The amplitude, center, and width   */
/*  of the gaussians are stored in consecutive locations of a: a[i] = Bk,    */
/*  a[i+1] = Ek, a[i+2] = Gk, k=1,...,na/3.  The dimensions of the arrays    */
/*  are a[1..na], dyda[1..na].                                               */
/*                                                                           */
/*          K         (    (x-Ek)^2 )                                        */
/*  y(x) = SUM Bk exp ( -  (----)   )                                        */
/*         k=1        (    ( Gk )   )                                        */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void fgauss(x,a,y,dyda,na)
     float x,a[],*y,dyda[];
     int na;
{
  int i;
  float fac,ex,arg;
  
  *y=0.0;
  for (i=1;i<=na-1;i+=3) {
    arg=(x-a[i+1])/a[i+2];
    ex=exp(-arg*arg);
    fac=a[i]*ex*2.0*arg;
    *y += a[i]*ex;
    dyda[i]=ex;
    dyda[i+1]=fac/a[i+2];
    dyda[i+2]=fac*arg/a[i+2];
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                  FGAUSS4                                  */
/*                                                                           */
/*****************************************************************************/
void fgauss4(x,a,y,dyda,na) // By wyeth, to add a yoffset param.
     float x,a[],*y,dyda[];
     int na;
{
  int i;
  float fac,ex,arg;
  
  *y=0.0;
  for (i=1;i<=na-1;i+=4){
    arg = (x-a[i+1])/(sqrt(2.0)*a[i+2]);
    ex = exp(-arg*arg);
    fac = a[i]*ex*2.0*arg;
    *y += a[i+3] + a[i]*ex;
    dyda[i] = ex;
    dyda[i+1] = fac/(sqrt(2.0)*a[i+2]);
    dyda[i+2] = fac*arg/(sqrt(2.0)*a[i+2]);
    dyda[i+3] = 1.0;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                   FSIN4                                   */
/*                                                                           */
/*  Fitting to a sine wave with up to 4 parameters that vary.  Evaluate      */
/*  the function at x:                                                       */
/*                                                                           */
/*    y = b + a*sin(fx-p)                                                    */
/*                                                                           */
/*  and return 'y' and the derivatives 'dyda' with respect to each of the    */
/*  four parameters.                                                         */
/*                                                                           */
/*    a[1] = a = amplitude                                                   */
/*    a[2] = b = vertical offset                                             */
/*    a[3] = p = phase                                                       */
/*    a[4] = f = frequency                                                   */
/*                                                                           */
/*****************************************************************************/
void fsin4(x,a,y,dyda)
     float x,a[],*y,dyda[];
{
  float arg;

  arg = a[4]*x - a[3];  /* fx+p */
  *y = a[2] + a[1]*sin(arg);

  printf(" ********* WYETH --------- DERIV WAS WRONG HERE.\n");
  exit_error("FEXP3","Check this");
  
  //dyda[1] = arg;
  dyda[1] = sin(arg);

  dyda[2] = 1.0;
  dyda[3] = -a[1]*cos(arg);
  dyda[4] = x*a[1]*cos(arg);
}
/**************************************-**************************************/
/*                                                                           */
/*                                   FEXP3                                   */
/*                                                                           */
/*  Fitting to an exponential with up to 3 params that vary.  Evaluate       */
/*  the function at x:                                                       */
/*                                                                           */
/*    y = b + a*exp(-x/t)                                                    */
/*                                                                           */
/*  and return 'y' and the derivatives 'dyda' with respect to each of the    */
/*  four parameters.                                                         */
/*                                                                           */
/*    a[1] = a = amplitude                                                   */
/*    a[2] = b = vertical offset                                             */
/*    a[3] = t = tau, time constant                                          */
/*                                                                           */
/*****************************************************************************/
void fexp3(x,a,y,dyda)
     float x,a[],*y,dyda[];
{
  float arg;

  arg = -x / a[3];  /* (-x/t) */
  *y = a[2] + a[1]*exp(arg);

  dyda[1] = exp(arg);
  dyda[2] = 1.0;
  dyda[3] = a[1] * exp(arg) * x / (a[3]*a[3]);
}
/**************************************-**************************************/
/*                                                                           */
/*                                  SIGPSYCH                                 */
/*                                                                           */
/*  y(x) = d - (d-0.5) e^(-(x/a)^b)                                          */
/*                                                                           */
/*****************************************************************************/
void sigpsych(x,p,y,dyda,na) /* Sigmoidal function for psychophysics */
     float x,p[],*y,dyda[];
     int na;
{
  float a,b,d,ex;

  /*printf("SIGPSYCH\n");*/

  a = p[1];
  b = p[2];
  d = p[3];

  /*printf("a=%f b=%f d=%f x=%f   x/a = %f\n",a,b,d,x,x/a);*/

if (a < 0.0){
  *y = -1.0;
  dyda[1] = -1.0/a;
  dyda[2] = 1.0;
  dyda[3] = 1.0;
  return;
}
  
  ex = exp(-pow(x/a,b));
  *y = d - (d - 0.5) * ex;
  dyda[1] = -(d-0.5)*(b/a)*pow(x/a,b) * ex;
  if (x==0.0)
    dyda[2] = 0.0;
  else
    dyda[2] = (d-0.5) * pow(x/a,b) * log(x/a) * ex;
  dyda[3] = 1.0 - ex;

  /*printf("   dy1 = %f  dy2 = %f  dy3 = %f\n",dyda[1],dyda[2],dyda[3]);*/
}
/**************************************-**************************************/
/*                                                                           */
/*****************************************************************************/
void flat_line(x,a,y,dyda,na) // By wyeth, to add a yoffset param
     float x,a[],*y,dyda[];
     int na;
{
  *y = a[1];
  dyda[1] = 1.0;
}
/**************************************-**************************************/
/*                                                                           */
/*                                  MY_TTEST                                 */
/*                                                                           */
/*  Independent 1-sample t-test, to test if mean is equal to 'mu'.           */
/*  Here, I have used betai by analogy to the NumRec "ttest" and "tutest".   */
/*                                                                           */
/*****************************************************************************/
void my_ttest(data,n,mu,rt,rprob)
     float *data;  // ******  WYETH [1...n] ********************
     int n;
     float mu;    // Mean to test against
     float *rt;    // t-value
     float *rprob; // prob value 
{
  void avevar();
  float betai();
  float ave1,var1,df,t;

  avevar(data,n,&ave1,&var1);

  df = n-1;
  t = (ave1 - mu) / sqrt(var1/(float)n);  // Wikipedia

  *rt = t;
  *rprob = betai(0.5*df,0.5,df/(df + t*t));
}
/**************************************-**************************************/
/*                                                                           */
/*                                   TTEST                                   */
/*                                                                           */
/*  Added Dec 30, 2008.                                                      */
/*                                                                           */
/*  T-test for difference of mean, when variance is assumed to be the same.  */
/*                                                                           */
/*C   (C) Copr. 1986-92 Numerical Recipes Software #.3.                      */
/*                                                                           */
/*****************************************************************************/
void ttest(data1,n1,data2,n2,t,prob)
     float data1[];
     unsigned long n1;
     float data2[];
     unsigned long n2;
     float *t;
     float *prob;
{
  void avevar();
  float betai();
  float var1,var2,svar,df,ave1,ave2;
  
  avevar(data1,n1,&ave1,&var1);
  avevar(data2,n2,&ave2,&var2);

  df=n1+n2-2;

  svar=((n1-1)*var1+(n2-1)*var2)/df;

  *t=(ave1-ave2)/sqrt(svar*(1.0/n1+1.0/n2));
  *prob=betai(0.5*df,0.5,df/(df+(*t)*(*t)));
}
/**************************************-**************************************/
/*                                                                           */
/*                                   TUTEST                                  */
/*                                                                           */
/*  Given the arrays data1[1..n1] and data2[1..n2], this routine returns     */
/*  Student's t as t, and its significance as prob, small values of prob     */
/*  indicating that the arrays have significantly different means.  The      */
/*  data arrays are allowed to be drawn from populations with unequal        */
/*  variances.                                                               */
/*                                                                           */
/*  This routine and its subroutines are from Numerical Recipes in C.        */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void tutest(data1,n1,data2,n2,t,prob)
     float data1[],data2[],*t,*prob;
     int n1,n2;
{
  float var1,var2,df,ave1,ave2;
  void avevar();
  float betai();
  
  avevar(data1,n1,&ave1,&var1);
  avevar(data2,n2,&ave2,&var2);
  *t=(ave1-ave2)/sqrt(var1/n1+var2/n2);
  df=(var1/n1+var2/n2)*(var1/n1+var2/n2)/
    ((var1/n1)*(var1/n1)/(n1-1)+(var2/n2)*(var2/n2)/(n2-1));
  *prob=betai(0.5*df,0.5,df/(df+(*t)*(*t)));
  //df=SQR(var1/n1+var2/n2)/(SQR(var1/n1)/(n1-1)+SQR(var2/n2)/(n2-1));
  //*prob=betai(0.5*df,0.5,df/(df+SQR(*t)));
}
/*************************************---*************************************/
#define NRANSI
/**************************************-**************************************/
/*                                                                           */
/*                                  SPLINE                                   */
/*                                                                           */
/*  Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e.,  */
/*  y_i = f(x_i), with x_1 < x_2 < ... < x_N, and given values yp1 and ypn   */
/*  for the first derivative of the interpolating function at points 1 and   */
/*  n, respectively, this routine returns an array y2[1..n] that contains    */
/*  the second derivatives of the interpolating function at the tabulated    */
/*  points x_i.  If yp1 and/or ypn are equal to 1x10^30 or larger, the       */
/*  routine is signaled to set the corresponding boundary condition for a    */
/*  natural spline, with zero second derivative on that boundary.            */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void spline(x,y,n,yp1,ypn,y2)
     float x[],y[];
     int n;
     float yp1,ypn,y2[];
{
  int i,k;
  float p,qn,sig,un,*u;
  
  u=vector(1,n-1);
  if (yp1 > 0.99e30)
    y2[1]=u[1]=0.0;
  else {
    y2[1] = -0.5;
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i=2;i<=n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  free_vector(u,1,n-1);
}
#undef NRANSI
/**************************************-**************************************/
/*                                                                           */
/*                                  SPLINT                                   */
/*                                                                           */
/*  Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with  */
/*  the xa_i's in order), and given the array y2a[1..n], which is the        */
/*  output from "spline" above, and given a value of x, this routine         */
/*  returns a cubic-spline interpolated value y.                             */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void splint(xa,ya,y2a,n,x,y)
     float xa[],ya[],y2a[];
     int n;
     float x,*y;
{
  int klo,khi,k;
  float h,b,a;
  
  klo=1;
  khi=n;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h == 0.0){
    printf("*** SPLINT nrerror\n");
    nrerror("Bad xa input to routine splint");
  }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}
