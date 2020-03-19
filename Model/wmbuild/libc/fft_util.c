/*****************************************************************************/
/*                                                                           */
/*   fft_util.c                                                              */
/*   wyeth bair                                                              */
/*   caltech                                                                 */
/*   02/07/93                                                                */
/*                                                                           */
/*   This file contains routines for computing the fft, power                */
/*   spectrum, autocorrelation and cross-correlation of trains               */
/*   of action potentials in the multi_channel (.mc) format.                 */
/*                                                                           */
/*   This code was taken from the "fft.c" program, which only                */
/*   works for the old single unit data format.                              */
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
#include "farray_util.h"
#include "myrand_util.h"
#include "spike_util.h"

/**************************************-**************************************/
/*                                                                           */
/*                                AB_TO_RTHETA                               */
/*                                                                           */
/*  Convert to phase and modulus.                                            */
/*                                                                           */
/*****************************************************************************/
void ab_to_rtheta(a,b,rr,rtheta)
     float a,b,*rr,*rtheta;
{
  float r,theta;

  r = sqrt(a*a + b*b);
  theta = atan2(b,a);

  *rr = r; *rtheta = theta;
}
/**************************************-**************************************/
/*                                                                           */
/*                                RTHETA_TO_AB                               */
/*                                                                           */
/*  Convert to phase and modulus.                                            */
/*                                                                           */
/*****************************************************************************/
void rtheta_to_ab(r,theta,ra,rb)
     float r,theta,*ra,*rb;
{
  float a,b;

  a = r*cos(theta);
  b = r*sin(theta);

  *ra = a; *rb = b;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 ROTATE_AB                                 */
/*                                                                           */
/*  Add phi to the phase angle of the complex number.  Phi in radians.       */
/*                                                                           */
/*****************************************************************************/
void rotate_ab(a,b,phi)
     float *a,*b,phi;
{
  float r,theta;

  ab_to_rtheta(*a,*b,&r,&theta);
  theta += phi;
  rtheta_to_ab(r,theta,a,b);
}
/**************************************-**************************************/
/*                                                                           */
/*                               THREE_D_POWER                               */
/*                                                                           */
/*  Compute the power spectrum from the arrays returned by NumRec's          */
/*  "rlft3" routine.                                                         */
/*                                                                           */
/*  The output spectrum comes back "logically" as a complex 3-d              */
/*  array which we will call SPEC[1..nn1][1..nn2][1..nn3/2+1]                */
/*  (cf. eq. 12.5.3).  In the first 2 of 3 dimensions, the respective        */
/*  frequency values f1 or f2 are stored in standard wrap-around             */
/*  order: 0 freq first, smallest positive f next, then the smallest         */
/*  negative f in the last index.  The 3rd dimension returns only            */
/*  the positive half of the frequency spectrum.  See Fig. p527.             */
/*                                                                           */
/*  Remember, H(-n) = H(n)*, where n is a vector of wave numbers,            */
/*  "*" denotes complex conjugate, and H is the 3-d transform.               */
/*                                                                           */
/*  The subscript range SPEC[1..nn1][1..nn2][1.nn3/2] is returned in         */
/*  the input array data[1..nn1][1..nn2][1..nn3], with the following         */
/*  correspondence:                                                          */
/*                                                                           */
/*    Re(SPEC[i1][i2][i3]) = data[i1][i2][2*i3-1]                            */
/*    Im(SPEC[i1][i2][i3]) = data[i1][i2][2*i3]                              */
/*                                                                           */
/*  The remaining "plane" of values, SPEC[1..nn1][1..nn2][nn3/2+1],          */
/*  is returned in the 2-d float array speq[1..nn1][1..2*nn2], with          */
/*  the correspondence:                                                      */
/*                                                                           */
/*    Re(SPEC[i1][i2][nn3/2+1]) = speq[i1][2*i2-1]                           */
/*    Im(SPEC[i1][i2][nn3/2+1]) = speq[i1][2*i2]                             */
/*                                                                           */
/*  The following is an example of the ordering of frequency data            */
/*  in the returned array along the first two dimensions in the case         */
/*  where n1=n2=8:                                                           */
/*                                                                           */
/*    INDEX  FREQUENCY                                                       */
/*      0    largest negative = -(n1/2 - 1)/(n1*delta)                       */
/*      1                                                                    */
/*      2    smallest negative                                               */
/*      3    0                                                               */
/*      4    smallest positive                                               */
/*      5                                                                    */
/*      6                                                                    */
/*      7    Nyquist = +- 1/(2*delta)                                        */
/*                                                                           */
/*****************************************************************************/
float ***three_d_power(data,speq,n1,n2,n3,zero_flag)
     float ***data,**speq;
     int n1,n2,n3;
     int zero_flag;  // set DC power to zero if 1
{
  int i,j,k;
  int pi,pj;
  int c1,c2;    // the output power[c1][c2][0] has data for freq=0
  float ***power;
  float a,b;

  //printf("  THREE_D_POWER\n");  // Used by 'mod_srf_util', thus no printing

  if (n1 == 1) // n1=1 for 2D fft's
    c1 = 0;
  else
    c1 = n1/2 - 1;
  c2 = n2/2 - 1;
  power = get_3d_farray(n1,n2,n3/2+1);

  for(i=1;i<=n1;i++)
    for(j=1;j<=n2;j++)
      for(k=1;k<=n3/2;k++){
	a = data[i][j][2*k-1];
	b = data[i][j][2*k];
	if (i <= c1+2)  pi = c1+i-1;
	else            pi = i - (c1+3);
	if (j <= c2+2)  pj = c2+j-1;
	else            pj = j - (c2+3);
	power[pi][pj][k-1] = a*a + b*b;
      }

  for(i=1;i<=n1;i++)
    for(j=1;j<=n2;j++){
      a = speq[i][2*j-1];
      b = speq[i][2*j];
      if (i <= c1+2) pi = c1+i-1;
      else           pi = i - (c1+3);
      if (j <= c2+2) pj = c2+j-1;
      else           pj = j - (c2+3);
      power[pi][pj][n3/2] = a*a + b*b;
    }
  
  if (zero_flag)
    power[c1][c2][0] = 0.0;

  return power;
}
/**************************************-**************************************/
/*                                                                           */
/*                              THREE_D_FFT_PROD                             */
/*                                                                           */
/*  Compute the product of two FFTs returned from NumRec's "rlft3".          */
/*  The data "data" and "speq" are complex arrays.  Array indices            */
/*  begin at 1.                                                              */
/*                                                                           */
/*  See comments in "three_d_power" for data format.                         */
/*                                                                           */
/*****************************************************************************/
void three_d_fft_prod(data1,speq1,data2,speq2,n1,n2,n3,rdata,rspeq)
     float ***data1,**speq1,***data2,**speq2;
     int n1,n2,n3;
     float ****rdata,***rspeq;
{
  int i,j,k;
  int k2,j2;
  float *t1,*t2,*t,t1r,t1i,t2r,t2i,*s1,*s2,*s;
  float ***data,**speq;

  data = f3tensor(1,n1,1,n2,1,n3); // USE THIS FOR 3D FFT
  speq = matrix(1,n1,1,2*n2);

  for(i=1;i<=n1;i++){
    s1 = speq1[i];
    s2 = speq2[i];
    s  =  speq[i];
    
    for(j=1;j<=n2;j++){

      j2 = 2*j;
      t1r = s1[j2-1];
      t1i = s1[j2];
      t2r = s2[j2-1];
      t2i = s2[j2];

      s[j2-1] = t1r*t2r - t1i*t2i;  // REAL
      s[j2]   = t1r*t2i + t1i*t2r;  // IMG

      /*

      speq[i][2*j-1] = speq1[i][2*j-1]*speq2[i][2*j-1] -
                       speq1[i][2*j]  *speq2[i][2*j]; // REAL
      speq[i][2*j]   = speq1[i][2*j-1]*speq2[i][2*j] +
                       speq1[i][2*j]  *speq2[i][2*j-1]; // IMAG
      */

      t1 = data1[i][j];
      t2 = data2[i][j];
      t  =  data[i][j];

      for(k=1;k<=n3/2;k++){

	k2 = 2*k;
	t1r = t1[k2-1];
	t1i = t1[k2];
	t2r = t2[k2-1];
	t2i = t2[k2];

	// (a+bi)(c+di) = (ac-bd) + (ad+bc)i
	t[k2-1] = t1r*t2r - t1i*t2i;  // REAL
	t[k2]   = t1r*t2i + t1i*t2r;  // IMG
      }
    }
  }

  *rdata = data;
  *rspeq = speq;
}
/**************************************-**************************************/
/*                                                                           */
/*                              GET_POWER_LAW_FFT                            */
/*                                                                           */
/*  Return a 3D function in Fourier space, formatted like the output of      */
/*  NumRec 'rlft3', that will be used to multiply another FFT.               */
/*                                                                           */
/*  This will be used to reshape a power spectrum,                           */
/*  for example, this can be used to create 1/f power spectrum noise.        */
/*                                                                           */
/*  See comments in "three_d_power" for data format.                         */
/*                                                                           */
/*  *** NOTE at f=0, the value is 1.0  ************                          */
/*                                                                           */
/*****************************************************************************/
void get_power_law_fft(fpow,n1,n2,n3,rdata,rspeq)
     float fpow;                // Raise 'f' to this power, e.g., -1.0
     int n1,n2,n3;              // n1 = 1 for 2D data
     float ****rdata,***rspeq;
{
  int i,j,k;
  float ***data,**speq,f1,f2,f3,fdist;

  data = f3tensor(1,n1,1,n2,1,n3); /*** USE THIS FOR 3D FFT ***/
  speq = matrix(1,n1,1,2*n2);  // f3 component is at Nyquist

  for(i=1;i<=n1;i++){

    // Compute f1 given index 'i',  Note, the smallest non-zero f is 1
    if (i <= n1/2)
      f1 = (float)(i-1);     // 0,1, ... (N/2-1)   [Positive]
    else
      f1 = (float)(n1-i+1);  // N/2, N/2-1 ... 1   [Negative]

    for(j=1;j<=n2;j++){

      // Compute f2 given index 'j'
      if (j <= n2/2)
	f2 = (float)(j-1);
      else
	f2 = (float)(n2-j+1);

      f3 = n3/2;  // N/2 = Cutoff;  All of 'speq' is at f3 Cutoff

      fdist = sqrt(f1*f1 + f2*f2 + f3*f3);
      if (fdist > 0.0)
	speq[i][2*j-1] = pow(fdist,fpow);    // Real
      else
	speq[i][2*j-1] = 1.0;
      
      speq[i][2*j]   = 0.0;                // Imag

      for(k=1;k<=n3/2;k++){ /*** (a+bi)(c+di) = (ac-bd) + (ad+bc)i ***/

	// Compute f3 given index 'k'
	f3 = (float)(k-1);

	fdist = sqrt(f1*f1 + f2*f2 + f3*f3);

	// WYETH HERE DEBUG
	// WYETH HERE DEBUG
	//if (j == 3)
	//printf("fdist = %f fpow= %f  PRW: %f\n",fdist,fpow,pow(fdist,fpow));


	if (fdist > 0.0)
	  data[i][j][2*k-1] = pow(fdist,fpow);   // Real
	else
	  data[i][j][2*k-1] = 1.0;

	data[i][j][2*k]   = 0.0;               // Imag
      }
    }
  }
  *rdata = data;
  *rspeq = speq;
}
/**************************************-**************************************/
/*                                                                           */
/*                             PRINT_COMPLEX_ARRAY                           */
/*                                                                           */
/*****************************************************************************/
void print_complex_array(complex_array,n)
     float *complex_array;
     int n;
{
  int i;

  printf("          Real        Imag\n");
  printf("        ========    ========\n");
  for (i=0;i<n;i++)
    printf("%6d  %8.4f    %8.4f\n",i,complex_array[2*i],complex_array[2*i+1]);
  printf("\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                             MAKE_COMPLEX_ARRAY                            */
/*                                                                           */
/*****************************************************************************/
float *make_complex_array(real_data,n)
     float *real_data;
     int n;
{
  int i;
  float *complex_data;

  complex_data = (float *)myalloc((2*n)*sizeof(float));

  for(i=0;i<n;i++){
    complex_data[2*i  ] = real_data[i];
    complex_data[2*i+1] = 0.0;
  }
  return complex_data;
}
/**************************************-**************************************/
/*                                                                           */
/*                     MAKE_COMPLEX_ARRAY_FROM_REAL_IMAG                     */
/*                                                                           */
/*****************************************************************************/
float *make_complex_array_from_real_imag(real_data,imag_data,n)
     float *real_data,*imag_data;
     int n;
{
  int i;
  float *complex_data;

  complex_data = (float *)myalloc((2*n)*sizeof(float));

  for(i=0;i<n;i++){
    complex_data[2*i  ] = real_data[i];
    complex_data[2*i+1] = imag_data[i];
  }
  return complex_data;
}
/**************************************-**************************************/
/*                                                                           */
/*                                    REAL                                   */
/*                                                                           */
/*   Return array containing real part of complex array.                     */
/*                                                                           */
/*****************************************************************************/
float *real(complex_signal,n)
     float *complex_signal;
     int n;
{
  int i;
  float *real_part;

  real_part = (float *)myalloc(n*sizeof(float));
  for (i=0;i<n;i++)
    real_part[i] = complex_signal[2*i];
  return real_part;
}
/**************************************-**************************************/
/*                                                                           */
/*                                    IMAG                                   */
/*                                                                           */
/*   Return array containing imaginary part of complex array.                */
/*                                                                           */
/*****************************************************************************/
float *imag(complex_signal,n)
     float *complex_signal;
     int n;
{
  int i;
  float *imag_part;

  imag_part = (float *)myalloc(n*sizeof(float));
  for (i=0;i<n;i++)
    imag_part[i] = complex_signal[2*i+1];
  return imag_part;
}
/**************************************-**************************************/
/*                                                                           */
/*                             PLOT_COMPLEX_ARRAY                            */
/*                                                                           */
/*****************************************************************************/
void plot_complex_array(outfile,name,append_flag,cdata,n)
     char outfile[],name[];
     int append_flag;
     float *cdata;
     int n;
{
  float *rdata,*idata;

  rdata = real(cdata,n);
  idata = imag(cdata,n);
  
  if (append_flag)
    append_farray_xy_plot(outfile,rdata,idata,n,name);
  else
    write_farray_xy_plot(outfile,name,rdata,idata,n);

  myfree(rdata); myfree(idata);
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_MIN_POWER_OF_TWO                           */
/*                                                                           */
/*   Returns the minimum number that is a power of two and greater           */
/*   than or equal to the given number.                                      */
/*                                                                           */
/*****************************************************************************/
int get_min_power_of_two(x)
     int  x;
{
  int total;
  
  total = 1;
  while (total < x)
    total = total * 2;
  return total;
}
/**************************************-**************************************/
/*                                                                           */
/*                            CONTORT_REAL_FARRAY                            */
/*                                                                           */
/*   Swap the first and second halves of a float array of even length.       */
/*                                                                           */
/*****************************************************************************/
void contort_real_farray(data,n)
     float *data;
     int n;
{
  int i;
  int half;
  float temp;

  half = n/2;
  if (2*half != n){
    printf(" n = %d\n",n);
    exit_error("CONTORT_REAL_FARRAY","Array length not even");
  }
  for (i=0;i<half;i++) {
    temp         = data[i];
    data[i]      = data[i+half];
    data[i+half] = temp;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           CONTORT_REAL_2D_FARRAY                          */
/*                                                                           */
/*   Swap quadrants of a float array of even lengths.                        */
/*                                                                           */
/*     0.......8......15                                                     */
/*     ----------------               ----------------                       */
/*  0 |                |             |----        ++++|                      */
/*  . |       o*       |             |---          +++|                      */
/*  . |     ooo***     |             |-              +|                      */
/*  . |    oooo****    |  Becomes -> |                |                      */
/*  4 |    ++++----    |             |                |                      */
/*  . |     +++---     |             |*              o|                      */
/*  . |       +-       |             |***          ooo|                      */
/*  7 |                |             |****        oooo|                      */
/*     ----------------               ----------------                       */
/*                                                                           */
/*   The coordinate [x0,y0] contains the origin.                             */
/*                                                                           */
/*****************************************************************************/
void contort_real_2d_farray(data,x0,y0,xn,yn)
     float **data;
     int x0,y0,xn,yn;
{
  int i,j;
  int halfx,halfy;
  float temp;

  halfx = xn/2;
  halfy = yn/2;
  if ((2*halfx != xn)||(2*halfy != yn)){
    printf("  xn = %d   yn = %d\n",xn,yn);
    exit_error("CONTORT_REAL_2D_FARRAY","Array length not even");
  }

  for (i=0;i<halfx;i++){
    for (j=0;j<halfy;j++){
      temp = data[x0+i][y0+j]; // SWAP FIRST AND THIRD QUADRANTS
      data[x0+i][y0+j] = data[x0+i+halfx][y0+j+halfy];
      data[x0+i+halfx][y0+j+halfy] = temp;
      temp = data[x0+i+halfx][y0+j]; // SWAP SECOND AND FOURTH QUADRANTS
      data[x0+i+halfx][y0+j] = data[x0+i][y0+j+halfy];
      data[x0+i][y0+j+halfy] = temp;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           CONTORT_REAL_3D_FARRAY                          */
/*                                                                           */
/*  Translate the contents of the farray so that [s1+n1/2,s2+n2/2,s3+n3/2]   */
/*  becomes [s1,s2,s3], and data wraps around as needed.                     */
/*                                                                           */
/*  Modified to handle a 2D array, i.e. when n1 == 1.                        */
/*                                                                           */
/*****************************************************************************/
void contort_real_3d_farray(data,s1,s2,s3,n1,n2,n3)
     float ***data;
     int s1,s2,s3,n1,n2,n3;
{
  int i,j,k;
  int half1;
  float temp;

  if (n1 == 1){
    half1 = 0;
  }else{
    half1 = n1/2;
    if (2*half1 != n1){
      printf(" n1 = %d\n",n1);
      exit_error("CONTORT_REAL_3D_FARRAY","Array length not even");
    }
  }

  for(i=0;i<n1;i++)
    contort_real_2d_farray(data[s1+i],s2,s3,n2,n3);
  for(i=0;i<half1;i++)
    for(j=0;j<n2;j++)
      for(k=0;k<n3;k++){
	temp = data[s1+i][s2+j][s3+k];
	data[s1+i][s2+j][s3+k] = data[s1+i+half1][s2+j][s3+k];
	data[s1+i+half1][s2+j][s3+k] = temp;
      }
}
/**************************************-**************************************/
/*                                                                           */
/*                        CONTORT_REAL_3D_FARRAY_SPACE                       */
/*                                                                           */
/*  This performs an operation similar to 'contort_real_2d_farray' on the    */
/*  2D slice for each value of t on data[x][y][t].                           */
/*                                                                           */
/*****************************************************************************/
void contort_real_3d_farray_space(data,s1,s2,s3,n1,n2,n3)
     float ***data;
     int s1,s2,s3,n1,n2,n3;
{
  int i,j,k;
  int halfx,halfy;
  float temp;

  halfx = n1/2;
  halfy = n2/2;
  if ((2*halfx != n1)||(2*halfy != n2))
    exit_error("CONTORT_REAL_3D_FARRAY_SPACE","Array length not even");

  for (i=0;i<halfx;i++){
    for (j=0;j<halfy;j++){
      for (k=s3;k<(s3+n3);k++){
	temp = data[s1+i][s2+j][k]; /*** SWAP 1st AND 3rd QUADRANTS ***/
	data[s1+i][s2+j][k] = data[s1+i+halfx][s2+j+halfy][k];
	data[s1+i+halfx][s2+j+halfy][k] = temp;
	temp = data[s1+i+halfx][s2+j][k]; /*** SWAP 2nd AND 4TH QUADRANTS ***/
	data[s1+i+halfx][s2+j][k] = data[s1+i][s2+j+halfy][k];
	data[s1+i][s2+j+halfy][k] = temp;
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                  POWER_2D                                 */
/*                                                                           */
/*****************************************************************************/
float **power_2d(data,xn,yn,zero_flag)
     float **data;
     int xn,yn;
     int zero_flag;
{
  int i,j;
  float ***tdata,***tpow;
  float **speq,**pow;

  tdata = convert_2d_to_tensor(data,xn,yn);
  contort_real_3d_farray(tdata,1,1,1,1,xn,yn);

  speq = matrix(1,1,1,2*xn);
  rlft3(tdata,speq,1,xn,yn,1); // compute FFT of 2D array (n1=1)

  tpow = three_d_power(tdata,speq,1,xn,yn,zero_flag);

  pow = get_2d_farray(xn,yn/2);

  for(i=0;i<xn;i++){
    for(j=0;j<yn/2;j++){
      pow[i][j] = tpow[0][i][j];
    }
  }

  free_matrix(speq,1,1,1,2*xn);
  free_f3tensor(tdata,1,1,1,xn,1,yn);
  free_3d_farray(tpow,1,xn,yn/2+1);

  return pow;
}
/**************************************-**************************************/
/*                                                                           */
/*                             APPLY_POWER_LAW_2D                            */
/*                                                                           */
/*  Take the FFT of the 2D array, multiply it by a power function, and       */
/*  invert the transform.                                                    */
/*                                                                           */
/*****************************************************************************/
float **apply_power_law_2d(data,xn,yn,fpow)
     float **data;
     int xn,yn;
     float fpow;
{
  int i,j;
  float ***tdata,***pl_data,***r;
  float **speq,**pl_speq,**r_s,fc,**r2d;

  tdata = convert_2d_to_tensor(data,xn,yn);
  contort_real_3d_farray(tdata,1,1,1,1,xn,yn);

  speq = matrix(1,1,1,2*xn);
  rlft3(tdata,speq,1,xn,yn,1); // compute FFT of 2D array (n1=1)
  //my_rlft3(tdata,speq,1,xn,yn,1); // test version for Java conversion

  // Get the power-law mask to multiply in Freq. Domain
  get_power_law_fft(fpow,1,xn,yn,&pl_data,&pl_speq);

  // Multiply in Frequency domain
  three_d_fft_prod(tdata,speq,pl_data,pl_speq,1,xn,yn,&r,&r_s);

  free_matrix(speq   ,1,1,1,2*xn);
  free_matrix(pl_speq,1,1,1,2*xn);
  free_f3tensor(tdata  ,1,1,1,xn,1,yn);
  free_f3tensor(pl_data,1,1,1,xn,1,yn);

  //  Inverse FFT
  rlft3(r,r_s,1,xn,yn,-1);
  //my_rlft3(r,r_s,1,xn,yn,-1);  // test version for java conversion
  free_matrix(r_s,1,1,1,2*xn);
  fc = 2.0/((float)(1*xn*yn));
  multiply_3d_farray(r,1,1,1,xn,1,yn,fc);

  contort_real_3d_farray(r,1,1,1,1,xn,yn);

  r2d = get_2d_farray(xn,yn);
  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++)
      r2d[i][j] = r[1][i+1][j+1];

  free_f3tensor(r,1,1,1,xn,1,yn);

  return r2d;
}
/**************************************-**************************************/
/*                                                                           */
/*                             APPLY_POWER_LAW_3D                            */
/*                                                                           */
/*  Take the FFT of the 3D array, multiply it by a power function, and       */
/*  invert the transform.                                                    */
/*                                                                           */
/*****************************************************************************/
float ***apply_power_law_3d(data,xn,yn,tn,fpow,tensor_flag)
     float ***data;
     int xn,yn,tn;
     float fpow;
     int tensor_flag;
{
  int i,j,k;
  float ***tdata,***pl_data,***r;
  float **speq,**pl_speq,**r_s,fc,***r3d;

  if (tensor_flag == 0){
    tdata = convert_3d_to_tensor(data,xn,yn,tn);
  }else
    tdata = data;
  
  contort_real_3d_farray(tdata,1,1,1,xn,yn,tn);

  speq = matrix(1,xn,1,2*yn);
  rlft3(tdata,speq,xn,yn,tn,1); // compute FFT of 2D array (n1=1)

  // Get the power-law mask to multiply in Freq. Domain
  get_power_law_fft(fpow,xn,yn,tn,&pl_data,&pl_speq);

  // Multiply in Frequency domain
  three_d_fft_prod(tdata,speq,pl_data,pl_speq,xn,yn,tn,&r,&r_s);

  free_matrix(speq   ,1,xn,1,2*yn);
  free_matrix(pl_speq,1,xn,1,2*yn);
  if (tensor_flag == 0)
    free_f3tensor(tdata  ,1,xn,1,yn,1,tn);
  free_f3tensor(pl_data,1,xn,1,yn,1,tn);

  //  Inverse FFT
  rlft3(r,r_s,xn,yn,tn,-1);
  free_matrix(r_s,1,xn,1,2*yn);
  fc = 2.0/((float)(xn*yn*tn));
  multiply_3d_farray(r,1,xn,1,yn,1,tn,fc);

  contort_real_3d_farray(r,1,1,1,xn,yn,tn);

  if (tensor_flag == 0){
    r3d = get_3d_farray(xn,yn,tn);
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	for(k=0;k<tn;k++)
	  r3d[i][j][k] = r[i+1][j+1][k+1];
    free_f3tensor(r,1,xn,1,yn,1,tn);
  }else{
    r3d = r;
  }

  return r3d;
}
/**************************************-**************************************/
/*                                                                           */
/*                           POWER_LAW_2D_TRANSFORM                          */
/*                                                                           */
/*  Replace the data with with a frame-by-frame power law transform.         */
/*                                                                           */
/*****************************************************************************/
void power_law_2d_transform(data,xn,yn,fpow,rescale)
     float **data;
     int xn,yn;
     float fpow;
     char *rescale;   // "none", "1 max 0.5 mid"
{
  int i,j;
  float **pow_data,min,max,amp;

  pow_data = apply_power_law_2d(data,xn,yn,fpow);

  get_min_max_2d_farray(pow_data,xn,yn,&min,&max);

  //  Set 'amp' to the largest of the two extreme values
  if (-min > max)
    amp = -min;
  else
    amp = max;

  if (strcmp(rescale,"1 max 0.5 mid")==0){
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	data[i][j] = pow_data[i][j] / amp;
  }else{
    exit_error("POWER_LAW_2D_TRANSFORM","Unknown 'rescale' option");
  }

  free_2d_farray(pow_data,xn);
}
/*************************************---*************************************/
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
/**************************************-**************************************/
/*                                                                           */
/*                                    FFT                                    */
/*                                                                           */
/*  This procedure copied from "Numercial Recipes in C, The Art of           */
/*  Scientific Computing", William H. Press, Saul A. Teukolsky,              */
/*  Brian P. Flannery, and William T. Vetterling, pp 411-412,                */
/*  Cambridge University Press, 1988.                                        */
/*                                                                           */
/*  Replaces "data" by its discrete Fourier transform, if "isign"            */
/*  is input as 1; or replaces data by nn times its inverse                  */
/*  discrete Fourier transform, if "isign" is input as -1.  "data"           */
/*  is a complex array of length "nn", input as a real array                 */
/*  "data[1..2*nn].  "nn" MUST be an integer power of 2 (this is             */
/*  not checked for!).                                                       */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void fft(data,nn,isign)
     float data[];
     int   nn;
     int   isign;
{
  int     n,mmax,m,j,istep,i;
  double  wtemp,wr,wpr,wpi,wi,theta;
  float   tempr,tempi;

  n = nn << 1;
  j = 1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m = n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax = 2;

  while (n > mmax) {
    istep = mmax << 1;
    theta = isign*(6.28318530717959/mmax);
    wtemp = sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m=1;m<mmax;m+=2){
      for (i=m;i<=n;i+=istep){
	j = i + mmax;
	tempr = wr * data[j] - wi*data[j+1];
	tempi = wr * data[j+1] + wi*data[j];
	data[j] = data[i] - tempr;
	data[j+1] = data[i+1] - tempi;
	data[i] += tempr;
	data[i+1] += tempi;
      }
      wr = (wtemp=wr)*wpr-wi*wpi+wr;
      wi = wi*wpr + wtemp*wpi + wi;
    }
    mmax = istep;
  }
}
#undef SWAP
/*************************************---*************************************/
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
/**************************************-**************************************/
/*                                                                           */
/*                                   FOUR1                                   */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void four1(data,nn,isign)
     float data[];
     unsigned long nn;
     int isign;
{
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  float tempr,tempi;

  //
  //
  //  *** WYETH THIS IS SAME as fft above except for TYPE unsigned long
  //  *** WYETH THIS IS SAME as fft above except for TYPE unsigned long
  //  *** WYETH THIS IS SAME as fft above except for TYPE unsigned long
  //   should replace fft() above with this???
  //
  
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
	j=i+mmax;
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i] += tempr;
	data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}
#undef SWAP
/*************************************---*************************************/
#define  WINDOW(j,a,b) (1.0-fabs((((j)-1)-(a))*(b))) // Parzen
/**************************************-**************************************/
/*                                                                           */
/*                                  SPCTRM                                   */
/*                                                                           */
/*   This procedure copied from "Numercial Recipes in C, The Art of          */
/*   Scientific Computing", William H. Press, Saul A. Teukolsky,             */
/*   Brian P. Flannery, and William T. Vetterling, pp 445-446,               */
/*   Cambridge University Press, 1988.                                       */
/*                                                                           */
/*   Reads data from input stream specified by file pointer "fp" and         */
/*   (modified to read from the array "data")                                */
/*   returns as p[j] the data's power (mean square amplitude) at             */
/*   frequency (j-1)/(2*m) cycles per gridpoint, for j=1,2,...,m,            */
/*   based on (2*k+1)*m data points (if "ovrlap is set TRUE (1)) or          */
/*   4*k*m data points (if "ovrlap" is set FALSE (0)).  The number           */
/*   of segments of the data is 2*k in both cases: the routine calls         */
/*   "fft" k-times, each call with 2 partitions each of 2*m real             */
/*   data points.                                                            */
/*                                                                           */
/*   wyeth's comments:                                                       */
/*   NOTE:  apparently this routine throws away the value at the             */
/*          cutoff frequency.  The values at locations 2m+1 and              */
/*          2m+2 are not put in the p[] array and thus not returned.         */
/*          I have modified the routine to return this value.                */
/*                                                                           */
/*   NOTE:  the final normalization "den" = sumw * k * 4 * m.                */
/*          This is correct since 2*k segments of length 2*m are             */
/*          averaged.  The "sumw" normalizes for the sum of the              */
/*          window squared.  (See page 442 in NumRecInC.)                    */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void spctrm(data,p,m,k,ovrlap) /* data used to be fp */
     float  *data; /* used to be:  FILE *fp; STARTS AT 0 */
     float   p[];
     int     m,k,ovrlap;
{
  int mm,m44,m43,m4,kk,joffn,joff,j2,j;
  int data_pos = 0; // current position in data array
  float w,facp,facm,*w1,*w2,sumw=0.0,den=0.0,wyy;
  float *vector();
  void fft(),free_vector();

  mm=m+m;             /* useful factors */
  m43=(m4=mm+mm)+3;
  m44=m43+1;
  w1=vector(1,m4);
  w2=vector(1,m);
  facm=m-0.5;
  facp=1.0/(m+0.5);

  for (j=1;j<=mm;j++){
    wyy = WINDOW(j,facm,facp);
    sumw += wyy * wyy;
    //sumw += SQR(WINDOW(j,facm,facp));
  }
  for (j=1;j<=m+1;j++) p[j]=0.0; /* init spectrum to zero */  /*wyeth*/
  if (ovrlap)
    for (j=1;j<=m;j++) {
      w2[j] = data[data_pos];
      data_pos += 1;
      /*fscanf(fp,"%f",&w2[j]);*/
    }
  for (kk=1;kk<=k;kk++) {
    for (joff= -1;joff<=0;joff++) {
      if (ovrlap) {
	for (j=1;j<=m;j++) w1[joff+j+j]=w2[j];
	for (j=1;j<=m;j++) {
	  w2[j] = data[data_pos];
	  data_pos += 1;
	  /*fscanf(fp,"%f",&w2[j]);*/
	}
	joffn=joff+mm;
	for (j=1;j<=m;j++) w1[joffn+j+j]=w2[j];
      } else {
	for (j=joff+2;j<=m4;j+=2) {
	  w1[j] = data[data_pos];
	  data_pos += 1;
	  /*fscanf(fp,"%f",&w1[j]);*/
	}
      }
    }
    for (j=1;j<=mm;j++) {
      j2=j+j;
      w=WINDOW(j,facm,facp);
      w1[j2] *= w;
      w1[j2-1] *= w;
    }
    fft(w1,mm,1);                     // fft the windowed data
    //p[1] += (SQR(w1[1])+SQR(w1[2]));  // sum results into prev. segs.
    p[1] += (w1[1]*w1[1] + w1[2]*w1[2]);  // sum results into prev. segs.
    //p[m+1] += (SQR(w1[2*m+1])+SQR(w1[2*m+2])); // cutoff freq.  WYETH
    p[m+1] += (w1[2*m+1]*w1[2*m+1] + w1[2*m+2]*w1[2*m+2]); // cutoff freq.
    for (j=2;j<=m;j++) {
      j2 = j+j;
      p[j] += (w1[j2]    *w1[j2]     + w1[j2-1]  *w1[j2-1] + 
	       w1[m44-j2]*w1[m44-j2] + w1[m43-j2]*w1[m43-j2]);
      //p[j] += (SQR(w1[j2]) + SQR(w1[j2-1])
      //+ SQR(w1[m44-j2]) + SQR(w1[m43-j2]));
    }
    den += sumw;
  }
  den *= m4;
  for (j=1;j<=m+1;j++) p[j] /= den;  /*wyeth*/
  free_vector(w2,1,m);
  free_vector(w1,1,m4);
}
#undef WINDOW
/**************************************-**************************************/
/*                                                                           */
/*                                  TWOFFT                                   */
/*                                                                           */
/*   Given two real input arrays "data1[1..n]" and "data2[1..n]",            */
/*   this routine calls "fft" and returns two complex output                 */
/*   arrays, "fft1" and "fft2", each of complex length n (i.e. real          */
/*   dimensions [1..2n]), which contain the discrete Fourier                 */
/*   transforms of the respective datas.  "n" must be an integer             */
/*   power of two.                                                           */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void twofft(data1,data2,fft1,fft2,n)
     float data1[],data2[],fft1[],fft2[];
     int n;
{
  int nn3,nn2,jj,j;
  float rep,rem,aip,aim;
  void fft();
  
  nn3=1+(nn2=2+n+n);
  for (j=1,jj=2;j<=n;j++,jj+=2) {
    fft1[jj-1]=data1[j];
    fft1[jj]=data2[j];
  }
  fft(fft1,n,1);
  fft2[1]=fft1[2];
  fft1[2]=fft2[2]=0.0;
  for (j=3;j<=n+1;j+=2) {
    rep=0.5*(fft1[j]+fft1[nn2-j]);
    rem=0.5*(fft1[j]-fft1[nn2-j]);
    aip=0.5*(fft1[j+1]+fft1[nn3-j]);
    aim=0.5*(fft1[j+1]-fft1[nn3-j]);
    fft1[j]=rep;
    fft1[j+1]=aim;
    fft1[nn2-j]=rep;
    fft1[nn3-j] = -aim;
    fft2[j]=aip;
    fft2[j+1] = -rem;
    fft2[nn2-j]=aip;
    fft2[nn3-j]=rem;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                   REALFT                                  */
/*                                                                           */
/*   Calculates the FT of a set of 2n real-valued data points.               */
/*   Replaces this data (which is stored in array "data[1..2n]" by           */
/*   the positive frequency half of its complex FT.  The real-valued         */
/*   first and last components of the complex transform are returned         */
/*   as elements "data[1]" and "data[2]" respectively.  "n" must be          */
/*   a power of 2.  This routine also calculates the inverse                 */
/*   transform of a complex data array if it is the transform of             */
/*   real data.  (Result in this case must be multiplied by 1/n.)            */
/*                                                                           */
/*C   (C) Copr. 1986-92 Numerical Recipes Software #.3.                      */
/*                                                                           */
/*****************************************************************************/
void realft(data,n,isign)
     float data[];
     unsigned long n;
     int isign;
{
  void four1();
  unsigned long i,i1,i2,i3,i4,np3;
  float c1=0.5,c2,h1r,h1i,h2r,h2i;
  double wr,wi,wpr,wpi,wtemp,theta;
  
  theta=3.141592653589793/(double) (n>>1);
  if (isign == 1) {
    c2 = -0.5;
    four1(data,n>>1,1);
  } else {
    c2=0.5;
    theta = -theta;
  }
  wtemp=sin(0.5*theta);
  wpr = -2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  np3=n+3;
  for (i=2;i<=(n>>2);i++) {
    i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
    h1r=c1*(data[i1]+data[i3]);
    h1i=c1*(data[i2]-data[i4]);
    h2r = -c2*(data[i2]+data[i4]);
    h2i=c2*(data[i1]-data[i3]);
    data[i1]=h1r+wr*h2r-wi*h2i;
    data[i2]=h1i+wr*h2i+wi*h2r;
    data[i3]=h1r-wr*h2r+wi*h2i;
    data[i4] = -h1i+wr*h2i+wi*h2r;
    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wi*wpr+wtemp*wpi+wi;
  }
  if (isign == 1) {
    data[1] = (h1r=data[1])+data[2];
    data[2] = h1r-data[2];
  } else {
    data[1]=c1*((h1r=data[1])+data[2]);
    data[2]=c1*(h1r-data[2]);
    four1(data,n>>1,-1);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                   CORREL                                  */
/*                                                                           */
/*   Computes the correlation of two real data sets "data1[1..n]"            */
/*   and "data2[1..n]", each of length "n" including any user-               */
/*   supplied zero-padding.  "n" must be a power of two.  The                */
/*   result is returned as the first "n" points in "ans[1..2*n]"             */
/*   stored in wraparound order, i.e. correlations at increasingly           */
/*   negative lags are in "ans[n]" on down to "ans[n/2+1]", while            */
/*   correlations at increasingly positive lags are in "ans[1]"              */
/*   (zero lag) on up to "ans[n/2]".  Note that "ans" must be                */
/*   supplied in the calling program with length at least 2*n                */
/*   since it is also used as working space.  Sign convention of             */
/*   this routine: if "data1" lags "data2", i.e. is shifted to the           */
/*   right of it, then "ans" will show a peak at positive lags.              */
/*                                                                           */
/*C   (C) Copr. 1986-92 Numerical Recipes Software #.3.                      */
/*                                                                           */
/*****************************************************************************/
void correl(data1,data2,n,ans)
     float data1[],data2[];
     int n;
     float ans[];
{
  int no2,i;
  float dum,*fft,*vector();
  void twofft(),realft(),free_vector();

  fft = vector(1,2*n);
  twofft(data1,data2,fft,ans,n);
  no2=n/2;
  for(i=2;i<=n+2;i+=2){
    ans[i-1] = (fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/no2;
    ans[i]   = (fft[i]*dum-fft[i-1]*ans[i])/no2;
  }
  ans[2]=ans[n+1];
  /*realft(ans,no2,-1);*/
  realft(ans,n,-1); // NumRecInC 2nd Ed.
  free_vector(fft,1,2*n);
}
/**************************************-**************************************/
/*                                                                           */
/*                                 CORRELATE                                 */
/*                                                                           */
/*  Compute the cross-correlation between two floating arrays.  Pad          */
/*  with zeros to avoid overlap.  Arrays may be different lengths,           */
/*  but are assumed to have zeroth elements aligned at zero lag.             */
/*                                                                           */
/*****************************************************************************/
float *correlate(data1,data2,n1,n2,ncc)
     float *data1,*data2;
     int n1,n2,*ncc;
{
  int i;
  int np2,nn;
  float *pad1,*pad2,*ans,*cc;

  if (n1 > n2)
    nn = n1;
  else
    nn = n2;
  cc = (float *)myalloc(2*nn*sizeof(float)); // only need 2*nn-1
  np2 = 2 * get_min_power_of_two(nn);

  pad1 = pad_end_farray(data1,n1,np2,0.0);
  pad2 = pad_end_farray(data2,n2,np2,0.0);
  ans = (float *)myalloc(2*np2*sizeof(float));
  correl(pad1-1,pad2-1,np2,ans-1);

  for (i=1;i<nn;i++){
    cc[i-1] = ans[i+np2-nn]; // most negative lag at 0
    cc[i+nn-1] = ans[i];
  }
  cc[nn-1] = ans[0]; // zero lag at period-1
  *ncc = 2*nn - 1;
  return cc;
}
/**************************************-**************************************/
/*                                                                           */
/*                                   CONVLV                                  */
/*                                                                           */
/*  Convolves or deconvolves a real data set data[1..n] (including           */
/*  any user-supplied zero padding) with a response function                 */
/*  respns[1..n].  The response function must be stored in wrap              */
/*  around order in the first m elements of respns, where m is an            */
/*  odd integer <= n.  Wrap around order means that the first half           */
/*  of the array contains the impulse response function at positive          */
/*  times, while the second half of the array contains the impulse           */
/*  response function at negative times, counting down from the              */
/*  highest element respns[m].  On input isign is +1 for convolution,        */
/*  -1 for deconvolution.  The answer is returned in the first n             */
/*  components of ans.  However, ans must be supplied in the calling         */
/*  program with dimensions [1..2*n], for consistency with twofft.           */
/*  n MUST be an integer power of two.                                       */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void convlv(data,n,respns,m,isign,ans)
     float data[];
     unsigned long n;
     float respns[];
     unsigned long m;
     int isign;
     float ans[];
{
  void realft(),twofft(),free_vector();
  unsigned long i,no2;
  float dum,mag2,*fft,*vector();

  //printf("n,m  = %d %d\n",n,m);

  fft=vector(1,n<<1);
  for (i=1;i<=(m-1)/2;i++)          /* put respns in array of length n */
    respns[n+1-i]=respns[m+1-i];
  for (i=(m+3)/2;i<=n-(m-1)/2;i++)  /* pad with zeros */
    respns[i]=0.0;
  twofft(data,respns,fft,ans,n);    /* FFT both at once */
  no2=n>>1;
  for (i=2;i<=n+2;i+=2){
    if (isign == 1){
      ans[i-1]=(fft[i-1]*(dum=ans[i-1])-fft[i]*ans[i])/no2;
      ans[i]=(fft[i]*dum+fft[i-1]*ans[i])/no2;   /* multiply to convolve */
    }else if (isign == -1){
      //if ((mag2=SQR(ans[i-1])+SQR(ans[i])) == 0.0)  /* divide to deconvolve */
      // divide to deconvolve
      if ((mag2 = ans[i-1]*ans[i-1] + ans[i]*ans[i]) == 0.0)
	exit_error("CONVLV","Deconvolving at response zero");
      ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/mag2/no2;
      ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/mag2/no2;
    }else
      exit_error("CONVLV","No meaning for ISIGN");
  }
  ans[2]=ans[n+1];           /* Pack last element with first for realft */
  realft(ans,n,-1);          /* Inverse transform back to time domain */
  free_vector(fft,1,n<<1);
}
/**************************************-**************************************/
/*                                                                           */
/*                                 DECONVOLVE                                */
/*                                                                           */
/*  Modification of 'convlv' from "Numerical Recipes in C".                  */
/*                                                                           */
/*****************************************************************************/
void deconvolve(outfile,data,n,respns,m,ans)
     char outfile[];
     float data[];
     unsigned long n;
     float respns[];
     unsigned long m;
     float ans[];
{
  int k;
  unsigned long i,no2;
  void realft(),twofft(),free_vector();
  float mag2,*fft,*vector();
  float m1,m2,*powd,*powr,*powa;

  printf("  DECONVOLVE\n");

  append_farray_plot(outfile,"data",data+1,n,1);
  append_farray_plot(outfile,"respns",respns+1,m,1);
  
  fft=vector(1,2*n);
  for (i=1;i<=(m-1)/2;i++)          /* put respns in array of length n */
    respns[n+1-i]=respns[m+1-i];
  for (i=(m+3)/2;i<=n-(m-1)/2;i++)  /* pad with zeros */
    respns[i]=0.0;
  twofft(data,respns,fft,ans,n);    /* FFT both at once */

  powd = (float *)myalloc(n*sizeof(float)); /* data */
  powr = (float *)myalloc(n*sizeof(float)); /* respnse */
  powa = (float *)myalloc(n*sizeof(float)); /* ans */

  no2=n/2;
  for(i=2;i<=n+2;i+=2){
    k = i/2 - 1;

    powd[k] = fft[i-1]*fft[i-1] + fft[i]*fft[i];
    powr[k] = ans[i-1]*ans[i-1] + ans[i]*ans[i];
    mag2 = powr[k];
    if (mag2 < 0.00001){
      printf("***DECONVOLVE k=%d\n",k);
      mag2 = 1.0;
      /*exit_error("DECONVOLVE","Cannot divide by 0.0");*/
    }
    
    m1 = fft[i-1]*ans[i-1] + fft[i]  *ans[i];
    m2 = fft[i]  *ans[i-1] - fft[i-1]*ans[i];
    ans[i-1] = m1/mag2/no2;
    ans[i] = m2/mag2/no2;

    powa[k] = ans[i-1]*ans[i-1] + ans[i]*ans[i];
  }
  ans[2]=ans[n+1];           /* Pack last element with first for realft */
  realft(ans,n,-1);          /* Inverse transform back to time domain */
  free_vector(fft,1,2*n);

  append_farray_plot(outfile,"pow_data",powd,n/2,1);
  append_farray_plot(outfile,"pow_respns",powr,n/2,1);
  append_farray_plot(outfile,"pow_ans",powa,n/2,1);
  append_farray_plot(outfile,"decon",ans+1,n,1);
}
/**************************************-**************************************/
/*                                                                           */
/*                                FFT_CONVOLVE                               */
/*                                                                           */
/*  Convolve two float arrays if isign = 1, deconvolve if -1.                */
/*                                                                           */
/*  NOTE: "response" is assumed to have data stored with its origin          */
/*        at the center of the array, and m must be ODD.                     */
/*                                                                           */
/*****************************************************************************/
float *fft_convolve(data,n,response,m,isign)
     float *data,*response;
     int n,m;
     int isign;
{
  int i;
  int np2;
  float *p2data;     /* data sent to "convlv", padded, power of two */
  float *wrap_resp;  /* "response" in wrap-around order */
  float *ans;        /* answer from "convlv" */
  float *n_ans;      /* short answer to return to caller */
  int center;        /* origin of response array */
  int wrap_padding;  /* no. 0's appended to data to avoid wrap overlap */

  if (m%2 == 0){
    printf("  *** m = %d\n",m);
    exit_error("FFT_CONVOLVE","m is not odd");
  }
  center = (m-1)/2;          /* index of the origin of the response */
  wrap_padding = 0;
  for (i=1;i<=center;i++)
    if ((response[center+i] != 0.0)||(response[center-i] != 0.0))
      wrap_padding = i;

  if (n+wrap_padding>=m) /* padded data must be at least as long as response */
    np2 = get_min_power_of_two(n+wrap_padding); /* and a power of two */
  else
    np2 = get_min_power_of_two(m);

  wrap_resp = (float *)myalloc(np2*sizeof(float)); /* long response array */
  wrap_resp[0] = response[center];  /* make wrap-around response */
  for (i=1;i<=center;i++){
    wrap_resp[i] = response[center+i];
    wrap_resp[i+center] = response[i-1];
  }
  p2data = get_zero_farray(np2); /* make padded power of two data */
  for (i=0;i<n;i++)
    p2data[i] = data[i];
  ans = (float *)myalloc(2*np2*sizeof(float)); /* NumRec needs extra space */

  convlv(p2data-1,np2,wrap_resp-1,m,isign,ans-1); /* use the NumRec routine */

  n_ans = (float *)myalloc(n*sizeof(float)); /* get the short answer */
  for(i=0;i<n;i++)
    n_ans[i] = ans[i];

  myfree(wrap_resp);
  myfree(p2data);
  myfree(ans);

  return n_ans;
}
/**************************************-**************************************/
/*                                                                           */
/*                                FFT_HILBERT                                */
/*                                                                           */
/*  Compute the Hilbert transform of the data.  Assume that "n" is           */
/*  odd and the origin is in the center of the array.                        */
/*                                                                           */
/*  I'm not sure what to do with the DC component, so I set it to 0.         */
/*  For the f_c component, I rotate it with the negative frequencies,        */
/*  although it contains the values for both negative and positive           */
/*  f_c.  Maybe it too should be set to 0.                                   */
/*                                                                           */
/*  The Hilbert transform of a real function should be real, of an           */
/*  odd should be even, and of an even should be odd.                        */
/*                                                                           */
/*  Here I take the Fourier transform of the Hilbert kernel to be            */
/*  the transfer function:                                                   */
/*                                                                           */
/*    -i sgn(w)                                                              */
/*                                                                           */
/*****************************************************************************/
float *fft_hilbert(data,n)
     float *data;
     int n;
{
  int i;
  int np2,center;
  float re,im,*wrap_data,*imarray,*rearray,*cdata,*unwrap;

  if (n%2 == 0){
    printf("  *** n = %d\n",n);
    exit_error("FFT_HILBERT","n is not odd");
  }
  center = (n-1)/2;         /* index of the origin of the data */
  np2 = get_min_power_of_two(n);

  wrap_data = get_zero_farray(np2);
  wrap_data[0] = data[center];  /* make wrap-around response */
  for (i=1;i<=center;i++){
    wrap_data[i] = data[center+i];
    wrap_data[i+center] = data[i-1];
  }

  cdata = make_complex_array(wrap_data,np2);
  fft(cdata-1,np2,1);
  /* Hilbert transform destroys DC component? */
  cdata[0] = 0.0;
  cdata[1] = 0.0;
  for(i=1;i<np2/2;i++){
    re = cdata[2*i];
    im = cdata[2*i+1];
    cdata[2*i] = im;  /* multiply by -i */
    cdata[2*i+1] = -re;
  }
  /*** NOT CLEAR WHAT TO DO WITH np2/2, which has plus and minus freq */
  for(i=np2/2;i<np2;i++){
    re = cdata[2*i];
    im = cdata[2*i+1];
    cdata[2*i] = -im;  /* multiply by i */
    cdata[2*i+1] = re;
  }
  fft(cdata-1,np2,-1);

  imarray = imag(cdata,np2);
  rearray = real(cdata,np2);
  multiply_farray(rearray,np2,1/(float)np2);
  multiply_farray(imarray,np2,1/(float)np2);

  /**write_farray_plot("zz.im",imarray,np2);
    write_farray_plot("zz.re",rearray,np2);**/

  unwrap = (float *)myalloc(n*sizeof(float));
  for(i=0;i<=center;i++)
    unwrap[center+i] = rearray[i];
  for(i=0;i<center;i++)
    unwrap[i] = rearray[np2-center+i];

  myfree(wrap_data); myfree(cdata); myfree(imarray); myfree(rearray);
  return unwrap;
}
/**************************************-**************************************/
/*                                                                           */
/*                              FFT_INTERPOLATE                              */
/*                                                                           */
/*  Expand the data by a factor "f" by multiplying the Fourier               */
/*  transform by a boxcar function.  Ends are padded with zeros to           */
/*  a power of two.                                                          */
/*                                                                           */
/*****************************************************************************/
float *fft_interpolate(data,n,f)
     float *data;
     int n,f;
{
  int i;
  int fn,p2n;
  float *fdata,*p2data,*cdata;

  // Create a sparse, expanded array containing the original samples
  fn = n*f;
  fdata = get_zero_farray(fn);
  for(i=0;i<n;i++)
    fdata[i*f] = data[i];

  // Compute the fft
  p2n = get_min_power_of_two(fn);
  p2data = pad_end_farray(fdata,fn,p2n,0.0); // no padding occurs if n==p2n
  cdata = make_complex_array(p2data,p2n);
  fft(cdata-1,p2n,1);

  // "Multiply" the spectrum by a boxcar function, scaled by "f"
  for(i=0;i<2*p2n;i++) // set high frequencies to zero
    if ((i >= p2n/f)&&(i<(2*p2n-p2n/f))) // set high frequencies to zero
      cdata[i] = 0.0;
    else
      cdata[i] *= f; //(float)f; Scale spectrum
  
  fft(cdata-1,p2n,-1);
  for (i=0;i<fn;i++)
    fdata[i] = cdata[2*i]/(float)p2n;

  myfree(cdata); myfree(p2data);
  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_FARRAY_AT_RESOLUTION                         */
/*                                                                           */
/*  Return the data at the specified sampling resolution.  "res_in" and      */
/*  "res_out" are in samples per second.                                     */
/*                                                                           */
/*  For example, res_in = 53.75 for 18.605 msec/frame stimulus.              */
/*               res_out = 1000.0 to return msec resolution.                 */
/*                                                                           */
/*****************************************************************************/
void get_farray_at_resolution(data,n,res_in,res_out,rdata,rn)
     float *data;
     int n;
     float res_in,res_out,**rdata;
     int *rn;
{
  int i;
  int nn,expf,nexp;
  float *expdata,*xexpdata,*ydata,*xdata;

  /*** Expand the data to beyond the desired resolution. ***/
  expf = (int)(1.0 + res_out/res_in);
  expdata = fft_interpolate(data,n,expf);
  nexp = expf*n;
  xexpdata = (float *)myalloc(nexp*sizeof(float));
  for(i=0;i<nexp;i++)
    xexpdata[i] = (float)i * res_out/(res_in*(float)expf);

  /*** Subsample the data at the desired resolution. ***/
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
/*                         GET_2D_FARRAY_AT_RESOLUTION                       */
/*                                                                           */
/*****************************************************************************/
void get_2d_farray_at_resolution(data,m,n,res_in,res_out,rdata,rn)
     float **data;
     int m,n;
     float res_in,res_out,***rdata;
     int *rn;
{
  int i;
  int tn;
  float **t;

  t = (float **)myalloc(m*sizeof(float *));
  for(i=0;i<m;i++)
    get_farray_at_resolution(data[i],n,res_in,res_out,&(t[i]),&tn);

  *rn = tn; *rdata = t;
}
/**************************************-**************************************/
/*                                                                           */
/*                          FFT_SMOOTH_WITH_GAUSSIAN                         */
/*                                                                           */
/*  Convolve the data with a Gaussian of standard deviation "sigma".         */
/*  There is no tolerance specified, since we will use a mask as             */
/*  large as the data array.  The result has "n" entries, and edges          */
/*  are handled by mirroring the data about the 0th and n-1st                */
/*  elements in the array.                                                   */
/*                                                                           */
/*****************************************************************************/
float *fft_smooth_with_gaussian(data,n,sigma)
     float *data;
     int n;
     float sigma;
{
  int i;
  float *mask,*long_data,*smooth,*result;
  int nmask,cmask,nlong;

  nmask = n;
  if (nmask%2 == 0) /* the "response" array must have odd length */
    nmask = n-1;

  if (sigma > 0.0){
    cmask = (nmask-1)/2;
    mask = (float *)myalloc(nmask*sizeof(float));
    for (i=0;i<nmask;i++)
      mask[i] = gaussian(sigma,(float)(i-cmask));
  }else{
    cmask = 0;
    nmask = 1;
    mask = (float *)myalloc(nmask*sizeof(float));
    mask[0] = 1.0;
  }

  nlong = n + nmask-1;
  long_data = (float *)myalloc(nlong*sizeof(float));
  for (i=0;i<cmask;i++)
    long_data[cmask-1-i] = data[i+1]; /* mirror data about 0th element */
  for (i=0;i<n;i++)
    long_data[cmask+i] = data[i]; /* fill in the main data */
  for (i=0;i<cmask;i++)
    long_data[cmask+n+i] = data[n-2-i]; /* mirror data about n-1st element */
    
  smooth = fft_convolve(long_data,nlong,mask,nmask,1);
  result = (float *)myalloc(n*sizeof(float));
  for (i=0;i<n;i++)
    result[i] = smooth[cmask+i];

  myfree(mask);
  myfree(long_data);
  myfree(smooth);

  return result;
}
/**************************************-**************************************/
/*                                                                           */
/*                         FFT_SMOOTH_2D_WITH_GAUSSIAN                       */
/*                                                                           */
/*  Smoothing by a Gaussian in 2D is done by smoothing in horizontal         */
/*  then vertical with 1D Gaussians, both with standard deviation            */
/*  sigma.                                                                   */
/*                                                                           */
/*****************************************************************************/
float **fft_smooth_2d_with_gaussian(data,m,n,sigma,tolerance)
     float **data;
     int m,n;
     float sigma,tolerance;
{
  int i,j;
  float *column,*smooth_column;
  float **smooth;

  printf("  FFT_SMOOTH_2D_WITH_GAUSSIAN\n");

  smooth = (float **)myalloc(m*sizeof(float *));
  /*printf("    Horizontal:");*/
  for (i=0;i<m;i++){
    if (i%100 == 0){
      printf(" %d",i);
      fflush(stdout);
    }
    smooth[i] = fft_smooth_with_gaussian(data[i],n,sigma);
  }
  
  /*printf("\n    Vertical:");*/
  column = (float *)myalloc(m*sizeof(float));
  for (j=0;j<n;j++){
    if (j%100 == 0){
      printf(" %d",j);
      fflush(stdout);
    }
    for (i=0;i<m;i++)
      column[i] = smooth[i][j];
    smooth_column = fft_smooth_with_gaussian(column,m,sigma);
    for (i=0;i<m;i++)
      smooth[i][j] = smooth_column[i];
    myfree(smooth_column);
  }
  printf("\n");
  myfree(column);
  return smooth;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 AUTO_CORR                                 */
/*                                                                           */
/*  The result is "autoc[0..2n-1]".                                          */
/*                                                                           */
/*****************************************************************************/
float *auto_corr(data,n)
     float *data;
     int    n;
{
  float *autoc;

  autoc = (float *)myalloc(2*n*sizeof(float));
  correl(data-1,data-1,n,autoc-1);
  return autoc;
}
/**************************************-**************************************/
/*                                                                           */
/*                              AUTOCORR_FARRAY                              */
/*                                                                           */
/*  Compute the autocorrelation of a floating array.  Pad with               */
/*  zeros.  Return the values from lag 0 to lag n-1.                         */
/*                                                                           */
/*****************************************************************************/
float *autocorr_farray(data,n)
     float *data;
     int n;
{
  int i;
  int p2n;
  float *ac,*acraw,*pdata;

  p2n = 2*get_min_power_of_two(n); // Make padded data
  pdata = get_zero_farray(p2n);
  for(i=0;i<n;i++)
    pdata[i] = data[i];

  acraw = auto_corr(pdata,p2n); // Compute raw autocorrelation
  myfree(pdata);

  ac = (float *)myalloc(n*sizeof(float)); // return positive lags
  for(i=0;i<n;i++)
    ac[i] = acraw[i];
  myfree(acraw);
  return ac;
}
/**************************************-**************************************/
/*                                                                           */
/*                              BANDPASS_FARRAY                              */
/*                                                                           */
/*  Set the spectrum of the farray to zero above the specified freq.         */
/*                                                                           */
/*****************************************************************************/
float *bandpass_farray(data,n,sampling,freq)
     float *data;
     int n;
     float sampling; /* samples per second */
     float freq;
{
  int i;
  int p2n;
  float *p2data,*cdata,f,*fdata;

  p2n = get_min_power_of_two(n);
  p2data = pad_end_farray(data,n,p2n,0.0); /* no padding occurs if n==p2n */
  cdata = make_complex_array(p2data,p2n);
  fft(cdata-1,p2n,1);

  for (i=1;i<p2n/2;i++){
    f = (float)i*sampling/(float)p2n;
    if (f > freq){
      cdata[2*i] = 0.0;
      cdata[2*i+1] = 0.0;
      cdata[2*(p2n-i)] = 0.0;
      cdata[2*(p2n-i)+1] = 0.0;
    }
  }
  f = sampling/2.0;
  if (f > freq){
    cdata[p2n] = 0.0; /*** Cutoff ***/
    cdata[p2n+1] = 0.0; /*** Cutoff ***/
  }
  
  fft(cdata-1,p2n,-1);
  fdata = real(cdata,p2n);
  multiply_farray(fdata,p2n,1.0/(float)p2n);
  myfree(p2data); myfree(cdata);

  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                                HIPASS_FARRAY                              */
/*                                                                           */
/*  Set the spectrum of the farray to zero below the specified freq.         */
/*                                                                           */
/*****************************************************************************/
float *hipass_farray(data,n,sampling,freq)
     float *data;
     int n;
     float sampling; // samples per second
     float freq;
{
  int i;
  int p2n;
  float *p2data,*cdata,f,*fdata;

  p2n = get_min_power_of_two(n);
  p2data = pad_end_farray(data,n,p2n,0.0); // no padding occurs if n==p2n
  cdata = make_complex_array(p2data,p2n);
  fft(cdata-1,p2n,1);

  for (i=1;i<p2n/2;i++){
    f = (float)i*sampling/(float)p2n;
    if (f < freq){
      cdata[2*i] = 0.0;
      cdata[2*i+1] = 0.0;
      cdata[2*(p2n-i)] = 0.0;
      cdata[2*(p2n-i)+1] = 0.0;
    }
  }
  f = sampling/2.0;
  if (f < freq){
    cdata[p2n] = 0.0;   // Cutoff
    cdata[p2n+1] = 0.0; // Cutoff
  }
  
  fft(cdata-1,p2n,-1);
  fdata = real(cdata,p2n);
  multiply_farray(fdata,p2n,1.0/(float)p2n);
  myfree(p2data); myfree(cdata);

  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GAUSSIAN_BANDPASS_FARRAY                         */
/*                                                                           */
/*  Multiply the spectra of the farray by a Gaussian, then invert            */
/*  and return the result.                                                   */
/*                                                                           */
/*****************************************************************************/
float *gaussian_bandpass_farray(data,n,sampling,fmu,fsigma)
     float *data;
     int n;
     float sampling; /* samples per second */
     float fmu,fsigma;
{
  int i;
  int p2n;
  float *p2data,*cdata,freq,g,*fdata;

  p2n = get_min_power_of_two(n);
  p2data = pad_end_farray(data,n,p2n,0.0); /* no padding occurs if n==p2n */
  cdata = make_complex_array(p2data,p2n);
  fft(cdata-1,p2n,1);

  g = exp(-0.5*(0.0-fmu)*(0.0-fmu)/(fsigma*fsigma));
  cdata[0] *= g; /*** DC Component ***/
  cdata[1] *= g; /*** DC Component ***/
  for (i=1;i<p2n/2;i++){
    freq = (float)i*sampling/(float)p2n;
    g = exp(-0.5*(freq-fmu)*(freq-fmu)/(fsigma*fsigma));
    cdata[2*i] *= g;
    cdata[2*i+1] *= g;
    cdata[2*(p2n-i)] *= g;
    cdata[2*(p2n-i)+1] *= g;
  }
  freq = sampling/2.0;
  g = exp(-0.5*(freq-fmu)*(freq-fmu)/(fsigma*fsigma));
  cdata[p2n] *= g; /*** Cutoff ***/
  cdata[p2n+1] *= g; /*** Cutoff ***/

  fft(cdata-1,p2n,-1);
  fdata = real(cdata,p2n);
  multiply_farray(fdata,p2n,1.0/(float)p2n);
  myfree(p2data); myfree(cdata);

  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                                    DFT                                    */
/*                                                                           */
/*****************************************************************************/
float *dft(real_data,n)
     float *real_data;
     int n;
{
  float *complex_data;

  if (power_of_two(n) == -1)
    exit_error("DFT","Data array length not a power of 2");
  
  contort_real_farray(real_data,n);
  complex_data = make_complex_array(real_data,n);
  contort_real_farray(real_data,n);

  /* fft wants array subscripts from 1 to 2n */
  fft(complex_data-1,n,1);

  return complex_data;
}
/**************************************-**************************************/
/*                                                                           */
/*                              POWER_OF_COMPLEX                             */
/*                                                                           */
/*   Return real array containing power spectrum of complex array.           */
/*   Here we simply take the square of the modulus of the DFT.               */
/*                                                                           */
/*****************************************************************************/
float *power_of_complex(complex_data,n)
     float  *complex_data;
     int     n;
{
  int     i;
  float  *power_spectrum;
  
  power_spectrum = (float *)myalloc(n*sizeof(float));
  for (i=0;i<n;i++)
    power_spectrum[i] = complex_data[2*i  ]*complex_data[2*i  ] +
                        complex_data[2*i+1]*complex_data[2*i+1];
  return power_spectrum;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MODULUS_OF_COMPLEX                            */
/*                                                                           */
/*   Return real array containing the modulus of the complex array.          */
/*                                                                           */
/*****************************************************************************/
float *modulus_of_complex(complex_data,n)
     float *complex_data;
     int n;
{
  int i;
  float *modu;
  
  modu = (float *)myalloc(n*sizeof(float));
  for (i=0;i<n;i++)
    modu[i] = sqrt(complex_data[2*i  ]*complex_data[2*i  ] +
		   complex_data[2*i+1]*complex_data[2*i+1]);
  return modu;
}
/**************************************-**************************************/
/*                                                                           */
/*                              PHASE_OF_COMPLEX                             */
/*                                                                           */
/*   Return real array containing phase spectrum of complex array.           */
/*                                                                           */
/*****************************************************************************/
float *phase_of_complex(complex_data,n,degflag)
     float *complex_data;
     int n,degflag;
{
  int i;
  float *phase,a,b;
  
  phase = (float *)myalloc(n*sizeof(float));
  for (i=0;i<n;i++){
    a = complex_data[2*i];
    b = complex_data[2*i+1];
    phase[i] = atan2(b,a);
    if (degflag)
      phase[i] *= 180.0/M_PI;
  }
  return phase;
}
/**************************************-**************************************/
/*                                                                           */
/*                                PHASE_UNWRAP                               */
/*                                                                           */
/*  Unwrap the phase, under the assumption that changes in phase are small   */
/*  between neighboring frequencies.                                         */
/*                                                                           */
/*  **** WARNING:  This is a hack, we should probably look at the complex    */
/*  data in the complex plane and find locations where the trajectory        */
/*  appears to cross the negative real axis (use 'plot_complex_array').      */
/*                                                                           */
/*****************************************************************************/
float *phase_unwrap(ph,n)
     float *ph;
     int n;
{
  int i;
  float totph,delta,*unph,dplus,dminus;

  unph = (float *)myalloc(n*sizeof(float));
  totph = 0.0;
  unph[0] = ph[0];
  for(i=1;i<n;i++){
    delta = ph[i] - ph[i-1];
    dplus = delta - 360.0;  /* *add* 360 to ph[i-1] */
    dminus = delta + 360.0; /* *subract* 360 from ph[i-1] */
    if (fabs(dplus) < fabs(delta))
      totph -= 360.0;
    else if (fabs(dminus) < fabs(delta))
      totph += 360.0;
    unph[i] = totph+ph[i];
  }
  return unph;
}
/**************************************-**************************************/
/*                                                                           */
/*                              POWER_OF_SPIKES                              */
/*                                                                           */
/*   Returns the modulus squared of the fft of the spike train.              */
/*                                                                           */
/*****************************************************************************/
float *power_of_spikes(data,n,start,period)
     int *data,n; // spike data
     int  start,period;
{
  float *fdata;
  float *cdata;
  float *pdata;
  
  // make array of floats from spike time data
  fdata = expand_spike_array(data,n,start,period);
  cdata = dft(fdata,period);
  //print_complex_array(cdata,10);
  pdata = power_of_complex(cdata,period);
  
  myfree(fdata);
  myfree(cdata);
  return pdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                          PERIODOGRAM_NORMALIZATION                        */
/*                                                                           */
/*   Takes fft results and does periodogram normalization as                 */
/*   described in Numerical Recipes In C, page 439.                          */
/*   The array returned has n/2 + 1  entries.                                */
/*                                                                           */
/*****************************************************************************/
float *periodogram_normalization(complex_data,n)
     float *complex_data;  /* complex array of fft data */
     int    n;             /* number of complex data points */
{
  int    i;
  float *pgram;
  float *power;
  float  n2;

  n2 = n*n;

  pgram = (float *)myalloc((n/2 + 1)*sizeof(float));
  power = power_of_complex(complex_data,n);
  pgram[0] = power[0];
  pgram[n/2] = power[n/2];
  for (i=1;i<=(n/2-1);i++)
    pgram[i] = power[i] + power[n-i];
  for (i=0;i<=(n/2);i++)
    pgram[i] = pgram[i] / n2;

  myfree(power);

  return pgram;
}
/**************************************-**************************************/
/*                                                                           */
/*                           PERIODOGRAM_OF_SPIKES                           */
/*                                                                           */
/*   Returns periodogram (NumRecInC) of the spike train.                     */
/*   Note: the array returned has 2/n + 1 entries.                           */
/*                                                                           */
/*****************************************************************************/
float *periodogram_of_spikes(data,n,start,period)
     int *data,n;        /* spike data */
     int  start,period;  /* window to process */
{
  float *fdata;
  float *cdata;
  float *pgram;
  
  /* make array of floats from spike time data */
  fdata = expand_spike_array(data,n,start,period);
  cdata = dft(fdata,period);
  /**print_complex_array(cdata,10);**/
  /* periodogram normalization */
  pgram = periodogram_normalization(cdata,period);

  myfree(fdata);
  myfree(cdata);
  return pgram;
}
/**************************************-**************************************/
/*                                                                           */
/*                            PWRWIN_TRIAL_SEGMENT                           */
/*                                                                           */
/*  Resulting array starts at 0, goes to m.  "spctrm" uses an input          */
/*  array starting at 0.                                                     */
/*                                                                           */
/*****************************************************************************/
float *pwrwin_trial_segment(data,n,start,m,k)
     int *data,n;  /* spike data */
     int start;    /* msec start time of segment */
     int m,k;      /* data windowing */
{
  int duration;
  float *pdata;
  float *sdata;
  
  duration = (2*k+1)*m;
  pdata = (float *)myalloc((m+1)*sizeof(float)); /* start array at 1 */
  sdata = expand_spike_array(data,n,start,duration);
  spctrm(sdata,pdata-1,m,k,1); /* send spike data, get power data */
  myfree(sdata);
  return pdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                             WINDOW_SPIKE_RATE                             */
/*                                                                           */
/*  Compute the average spike rate, normalized by taking account of          */
/*  the overlapping windows used by SPCTRM.                                  */
/*                                                                           */
/*  *** Further correction should be made for the shape of the               */
/*      window!                                                              */
/*                                                                           */
/*  ****** ALSO NOTE: a better algorithm can take account of the             */
/*         fact that ALL spikes should get sent twice except for             */
/*         the ones in the first half of the first window and the            */
/*         last half of the last window.  The triangular Parzen              */
/*         window should weight each spike the same when you take            */
/*         into account both appearance of the spike, except for             */
/*         those end spikes mentioned above.                                 */
/*                                                                           */
/*****************************************************************************/
float window_spike_rate(data,n,m,k,start)
     int *data;    /* spike data */
     int n;        /* trial number */
     int m,k;      /* SPCTRM params */
     int start;    /* spike train time window */
{
  int   i;
  float spike_rate;
  int   spike_count;
  int   windows,win_period;
  
  windows = 2*k;
  win_period = 2*m;
  
  spike_count = 0;
  for(i=0;i<windows;i++)
    spike_count += count_spikes(data,n,start+i*m,win_period);
  spike_rate = (float)spike_count / ((float)(win_period*windows)/1000.0);

  return spike_rate;
}
/**************************************-**************************************/
/*                                                                           */
/*                       NORMALIZE_PWRWIN_BY_SPIKE_RATE                      */
/*                                                                           */
/*  On 3/4/93 this was changed to be consistent with                         */
/*  "pwrwin_trial_segment" and "avg_pwr_sarray" in "pcorr_util.c"            */
/*  in that they assume the power spectrum array is indexed from             */
/*  zero to m, having m+1 entries.                                           */
/*                                                                           */
/*****************************************************************************/
void normalize_pwrwin_by_spike_rate(pdata,sdata,n,m,k,start,sampling)
     float *pdata;     /* pwrwin spectrum for "sdata" */ 
     int   *sdata,n;   /* spike data */
     int    m,k,start; /* windowing parameters */
     float  sampling;  /* samples per second */
{
  int    i;
  float  c,rate;
  
  rate = window_spike_rate(sdata,n,m,k,start);
  if (rate > 0.0){
    c = (float)m*sampling/rate;
    for(i=0;i<=m;i++)
      pdata[i] *= c;
  }else
    for(i=0;i<=m;i++)
      pdata[i] = 0.0;
}
/**************************************-**************************************/
/*                                                                           */
/*                             NORMALIZE_AUTOCORR                            */
/*                                                                           */
/*  Normalize the autocorrelation using the total number of spikes           */
/*  in the period.  The total number of spikes is stored in                  */
/*  "acdata[0]".                                                             */
/*                                                                           */
/*****************************************************************************/
void normalize_autocorr(acdata,period)
     float *acdata;  /* autocorrelation data */
     int    period;  /* windowing duration */
{
  int    i;
  float  spikes_per_bin,spb2;
  float  overlap;
  
  if (acdata[0] > 0.0){
    spikes_per_bin = acdata[0]/(float)period;
    spb2 = spikes_per_bin * spikes_per_bin;
    for(i=1;i<period;i++){ /* let zero bin remain as number of spikes */
      overlap = (float)(period - i);
      acdata[i] -= overlap * spb2;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            WRITE_PWRWIN_RECORD                            */
/*                                                                           */
/*****************************************************************************/
void write_pwrwin_record(fout,data,name,corr,dir,resp,start,m,k)
     FILE   *fout;
     float  *data;
     char    name[];
     int     corr,dir,resp,start;
     int     m,k;
{
  int i;

  fprintf(fout,"%s\n %d %d %d %d %d %d\n",name,corr,dir,resp,start,k,m);
  for (i=0;i<m;i++)  /* don't write cutoff datum, index (m) */
    fprintf(fout,"%.3e ",data[i]);
  fprintf(fout,"\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_PADDED_SPIKES                             */
/*                                                                           */
/*   Return padded and expanded version of the spike train for use           */
/*   in correlation.  Only process the window defined by "start"             */
/*   and "duration".                                                         */
/*                                                                           */
/*****************************************************************************/
float *get_padded_spikes(data,n,start,period)
     int *data,n; /* spike data */
     int  start,period;
{
  int     j;
  float  *sdata;
  float  *padded;
  
  sdata = expand_spike_array(data,n,start,period);
  padded = (float *)myalloc((2*period)*sizeof(float));

  for (j=0;j<period;j++)
    padded[j] = sdata[j];
  for (j=period;j<2*period;j++)
    padded[j] = 0.0;
  
  myfree(sdata);
  return padded;
}
/**************************************-**************************************/
/*                                                                           */
/*                                AUTO_SPIKES                                */
/*                                                                           */
/*   Use spike data within the window specified by "start" and               */
/*   "period".  If "period" is not a power of 2, then use the next           */
/*   higher power of 2.  Pad with zero's on each side if "pad_flag".         */
/*   The result is "acdata[0..period-1]" representing lags from              */
/*   0 to period-1.                                                          */
/*                                                                           */
/*****************************************************************************/
float *auto_spikes(data,n,start,period,pad_flag)
     int *data,n;        /* spike time array */
     int  start,period;  /* process spikes within this window */
     int  pad_flag;      /* if 1, pad with zeros */
{
  int    j;
  float *acdata;
  float *ac;
  float *sdata;
  int    scount;
  float *padded;
  int   *p2data,p2n,p2period;

  acdata = (float *)myalloc(period*sizeof(float));
  scount = count_spikes(data,n,start,period);
  if (scount > 0){ /* if there are spikes */
    if (power_of_two(period) == -1){
      p2data = extract_spikes(data,n,start,period,&p2n,0);
      p2period = get_min_power_of_two(period);
    }else{
      p2data = data;
      p2n = n;
      p2period = period;
    }
    if (pad_flag){ /* pad beginning and end with zeros */
      padded = get_padded_spikes(p2data,p2n,start,p2period);
      ac = auto_corr(padded,2*p2period);
      myfree(padded);
    }else{ /* no padding, except power of 2, data is WRAPPED AROUND */
      sdata = expand_spike_array(p2data,p2n,start,p2period);
      ac = auto_corr(sdata,p2period);
      myfree(sdata);
    }
    /* The first "period" elements represent all positive lags. */
    /* Negative lags are redundant */
    for (j=0;j<period;j++)
      acdata[j] = ac[j];
    myfree(ac);
  }else{ /* set the spectrum to zero */
    for (j=0;j<period;j++)
      acdata[j] = 0.0;
  }
  return acdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                            AUTO_CONVOLVE_SPIKES                           */
/*                                                                           */
/*   Use spike data within the window specified by "start" and               */
/*   "period".  If "period" is not a power of 2, then use the next           */
/*   higher power of 2.  Pad with zero's on each side if "pad_flag".         */
/*   The result is "acdata[0..period-1]" representing lags from              */
/*   0 to period-1.                                                          */
/*                                                                           */
/*****************************************************************************/
float *auto_convolve_spikes(data,n,start,period)
     int *data,n;        /* spike time array */
     int  start,period;  /* process spikes within this window */
{
  int    i;
  float *result;
  int    scount;
  float *sdata,*cntr_data;
  int odd_per;
  
  scount = count_spikes(data,n,start,period);
  if (scount > 0){ /* if there are spikes */
    odd_per = period;
    if (period%2 == 0)  /* if the period is even, add one spot */
      odd_per = period + 1;
    sdata = expand_spike_array(data,n,start,odd_per);
    if (period < odd_per)
      sdata[period] = 0.0;
    
    cntr_data = get_zero_farray(2*odd_per-1);
    for (i=0;i<period;i++) /* leave out possible extra spot */
      cntr_data[(odd_per-1)/2+i] = sdata[i];
      
    result = fft_convolve(cntr_data,2*odd_per-1,sdata,odd_per,1);
    myfree(sdata);
    myfree(cntr_data);
  }else /* set the spectrum to zero */
    result = get_zero_farray(2*period-1);

  return result;
}
/**************************************-**************************************/
/*                                                                           */
/*                                CROSS_SPIKES                               */
/*                                                                           */
/*   Use spike data from 0 to m-1 where m is a power of two.  Pad            */
/*   both sides with m/2 zero's if "pad_flag" is set to 1.  The              */
/*   result is an "m" element array with lags from -m/2 to m/2.              */
/*                                                                           */
/*   The result has 2*period-1 entries, starting with -(period-1),           */
/*                Entry    Lag                                               */
/*                    0    -(period-1)                                       */
/*             period-1    0                                                 */
/*           2*period-2    period-1                                          */
/*                                                                           */
/*****************************************************************************/
float *cross_spikes(data1,n1,data2,n2,start1,start2,period,pad_flag)
     int *data1,n1,*data2,n2;
     int  start1,start2;
     int  period;
     int  pad_flag;
{
  int    k;
  float *ccdata,*cc;
  float *sdata1,*sdata2;
  float *padded1,*padded2;
  int   *p2data1,p2n1,*p2data2,p2n2,p2period;

  ccdata = (float *)myalloc(2*period*sizeof(float)); // only need 2*per-1

  //  *** WYETH ***
  //  *** WYETH ***
  //  *** WYETH - THE LOGIC OF THE 'if else' SEEMS WRONG TO ME ***********
  //    E.g., why are we extracting a new spike time array in the 'if' and
  //      using hte original in the 'else'?
  //    Why isn't the 'if' just setting 'p2period' ???
  //
  if (power_of_two(period) == -1){
    p2data1 = extract_spikes(data1,n1,start1,period,&p2n1,0);
    p2data2 = extract_spikes(data2,n2,start2,period,&p2n2,0);
    p2period = get_min_power_of_two(period);
  }else{
    p2data1 = data1;
    p2data2 = data2;
    p2n1 = n1;
    p2n2 = n2;
    p2period = period;
  }

  if (pad_flag){ // pad beginning and end with zeros
    padded1 = get_padded_spikes(p2data1,p2n1,start1,p2period);
    padded2 = get_padded_spikes(p2data2,p2n2,start2,p2period);

    cc = (float *)myalloc(2*2*p2period*sizeof(float));
    correl(padded1-1,padded2-1,2*p2period,cc-1);

    myfree(padded1); myfree(padded2);
    // put data in order from most negative lag to most positive
    for (k=1;k<period;k++){
      ccdata[k-1] = cc[k+2*p2period-period]; // most negative lag at 0
      ccdata[k+period-1] = cc[k];
    }
    ccdata[period-1] = cc[0]; // zero lag at period-1
  }else{ // no padding
    sdata1 = expand_spike_array(p2data1,p2n1,start1,p2period);
    sdata2 = expand_spike_array(p2data2,p2n2,start2,p2period);
    cc = (float *)myalloc(2*p2period*sizeof(float));
    correl(sdata1-1,sdata2-1,p2period,cc-1);
    myfree(sdata1); myfree(sdata2);
    for (k=1;k<period;k++){
      ccdata[k-1] = cc[k+p2period-period];
      ccdata[k+period-1] = cc[k];
    }
    ccdata[period-1] = cc[0];
  }
  myfree(cc);

  return ccdata; // returned array will start at index 0
}
/**************************************-**************************************/
/*                                                                           */
/*                            WRITE_AUTOCORR_RECORD                          */
/*                                                                           */
/*****************************************************************************/
void write_autocorr_record(fout,data,name,corr,dir,resp,start,n)
     FILE *fout;
     float *data;
     char name[];
     int corr,dir,resp,start;
     int n;
{
  int i;
  
  fprintf(fout,"%s\n %d %d %d %d %d\n",name,corr,dir,resp,start,n);
  for (i=0;i<n;i++)
    fprintf(fout,"%.3e ",data[i]);
  fprintf(fout,"\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                               DATA_TRANS_INV                              */
/*                                                                           */
/*   Takes same input file format as "plain_power".                          */
/*                                                                           */
/*   Compute the Fourier transform of a data file, then compute              */
/*   the inverse transform of that transform.  Write the data,               */
/*   the transform, and the inverse transform of the transform               */
/*   to the output file.                                                     */
/*                                                                           */
/*****************************************************************************/
void data_trans_inv(infile,outfile,duration,mean_flag)
     char infile[];
     char outfile[];
     int duration;
     int mean_flag;
{
  FILE *fopen(),*fin;
  int i;
  int n,ns;
  float sampling,*data,*cdata,*pdata,*tdata,mean,sdev;

  printf("  DATA_TRANS_INV\n");

  // get the input data
  fin = fopen(infile,"r");
  ns = fscanf(fin,"%d %f",&n,&sampling);

  if (duration > 0)
    n = duration;
  
  if (power_of_two(n) == -1){
    printf("n = %d\n",n);
    exit_error("DATA_TRANS_INV","Data length must be integer power of two");
  }

  data = (float *)myalloc(n*sizeof(float));
  for (i=0;i<n;i++)
    ns = fscanf(fin,"%f",&data[i]);
  fclose(fin);
  
  if (mean_flag){
    mean_sdev_farray(data,n,&mean,&sdev);
    /*printf("    Subtracting mean\n");*/
    for (i=0;i<n;i++)
      data[i] -= mean;
  }
  
  /* contort_real_farray(data,n); */

  /* INPUT */
  append_farray_plot(outfile,"Input",data,n,1);

  /*
    fout = fopen(outfile,"w");
    fprintf(fout,"%d\n",6);
    
    fprintf(fout,"Input  %d\n",n);
    for (i=0;i<n;i++)
    fprintf(fout,"%f ",data[i]);
    fprintf(fout,"\n");*/


  cdata = make_complex_array(data,n);
  fft(cdata-1,n,1);
  pdata = power_of_complex(cdata,n);

  /* TRANSFORM (real and imaginary */

  tdata = get_zero_farray(n);
  for (i=0;i<n;i++)
    tdata[i] = cdata[2*i];
  append_farray_plot(outfile,"TransformReal",tdata,n,1);

  for (i=0;i<n;i++)
    tdata[i] = cdata[2*i+1];
  append_farray_plot(outfile,"TransformImag",tdata,n,1);

  for (i=0;i<n;i++)
    tdata[i] = 180.0/M_PI * atan2(cdata[2*i+1],cdata[2*i]);
  append_farray_plot(outfile,"Transform_Phase(deg)",tdata,n,1);

  for (i=0;i<n;i++)
    tdata[i] = sqrt(cdata[2*i+1]*cdata[2*i+1] + cdata[2*i]*cdata[2*i]);
  append_farray_plot(outfile,"Transform_Amp",tdata,n,1);

  /* POWER (real) */

  append_farray_plot(outfile,"Power",pdata,n,1);

  /*fprintf(fout,"Power  %d\n",n);
    for (i=0;i<n;i++)
    fprintf(fout,"%f ",pdata[i]);
    fprintf(fout,"\n");*/

  fft(cdata-1,n,-1);

  /* INVERSE TRANSFORM (real and imaginary) */

  for (i=0;i<n;i++)
    tdata[i] = cdata[2*i]/(float)n;
  append_farray_plot(outfile,"InvReal",tdata,n,1);

  for (i=0;i<n;i++)
    tdata[i] = cdata[2*i+1]/(float)n;
  append_farray_plot(outfile,"InvImag",tdata,n,1);


  /*
    fprintf(fout,"InvReal  %d\n",n);
    for (i=0;i<n;i++)
    fprintf(fout,"%f ",cdata[2*i]/(float)n);
    fprintf(fout,"\n");
    fprintf(fout,"InvImag  %d\n",n);
    for (i=0;i<n;i++)
    fprintf(fout,"%f ",cdata[2*i+1]/(float)n);
    fprintf(fout,"\n");
    fclose(fout);*/

  myfree(tdata);

  /*printf("  ***Use PLAIN style in MYPLOT.\n");*/
}
/**************************************-**************************************/
/*                                                                           */
/*                                POWER_FARRAY                               */
/*                                                                           */
/*  If the input array is not a power of two, pad it with zeros to           */
/*  the next greater power of two.  Compute the modulus squared of           */
/*  the fft.  Return the one-sided spectrum and the x-coordinates            */
/*  in Hz.                                                                   */
/*                                                                           */
/*  WYETH:  Is the center of the array at data[0]?                           */
/*                                                                           */
/*****************************************************************************/
void power_farray(data,n,sampling,power,xfreq,pn)
     float *data;
     int n;
     float sampling; /* samples per second */
     float **power,**xfreq; /* returned power and x coordinates */
     int *pn; /* number of points in returned spectrum */
{
  int i;
  int p2n;
  float *p2data,*cdata,*pdata;

  p2n = get_min_power_of_two(n);
  p2data = pad_end_farray(data,n,p2n,0.0); /* no padding occurs if n==p2n */
  cdata = make_complex_array(p2data,p2n);
  fft(cdata-1,p2n,1);
  pdata = power_of_complex(cdata,p2n);

  *pn = p2n/2 + 1;
  *power = (float *)myalloc(*pn*sizeof(float));
  *xfreq = (float *)myalloc(*pn*sizeof(float));
  (*power)[0] = pdata[0];
  (*xfreq)[0] = 0.0;
  for (i=1;i<p2n/2;i++){
    (*power)[i] = pdata[i] + pdata[p2n-i];
    (*xfreq)[i] = sampling/p2n*(float)i;
  }
  (*power)[p2n/2] = 2.0*pdata[p2n/2];
  (*xfreq)[p2n/2] = sampling/2.0; /*** WYETH ??? ***/

  myfree(p2data); myfree(cdata); myfree(pdata);
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_MOD_PHASE_FARRAY                          */
/*                                                                           */
/*  If the input array is not a power of two, pad it with zeros to the next  */
/*  greater power of two.  Compute the modulus of the fft and the phase in   */
/*  radians.                                                                 */
/*                                                                           */
/*  NOTES                                                                    */
/*  - Phase may wrap-around, is there a simple fix?                          */
/*                                                                           */
/*****************************************************************************/
void get_mod_phase_farray(data,n,sampling,degflag,plotflag,rmod,rpha,rxfreq,rn)
     float *data;
     int n;
     float sampling; /* samples per second */
     int degflag,plotflag;
     float **rmod,**rpha,**rxfreq;
     int *rn;
{
  int i;
  int p2n,nn;
  float *p2data,*cdata,*pmod,*ppha,*mod,*pha,*xfreq;

  p2n = get_min_power_of_two(n);
  p2data = pad_end_farray(data,n,p2n,0.0); /* no padding occurs if n==p2n */
  cdata = make_complex_array(p2data,p2n);
  fft(cdata-1,p2n,1);
  pmod = modulus_of_complex(cdata,p2n);
  ppha = phase_of_complex(cdata,p2n,degflag);

  if (plotflag)
    plot_complex_array("zzz.complex.pl","cx",1,cdata,100);
  
  /*** Make one-sided estimate. ***/
  nn = p2n/2 + 1;
  mod = (float *)myalloc(nn*sizeof(float));
  pha = (float *)myalloc(nn*sizeof(float));
  xfreq = (float *)myalloc(nn*sizeof(float));
  mod[0] = 2.0*pmod[0]; /*** WYETH - apparently this gets doubled ***/
  pha[0] = ppha[0];
  xfreq[0] = 0.0;
  for (i=1;i<p2n/2;i++){
    mod[i] = pmod[i] + pmod[p2n-i];
    pha[i] = ppha[i];  /**** Adding them = 0, they are negatives */
    xfreq[i] = sampling/p2n*(float)i;
    /*printf("%f %f %f\n",ppha[i],ppha[p2n-i],ppha[i] + ppha[p2n-i]);*/
  }
  mod[p2n/2] = 2.0*pmod[p2n/2];
  pha[p2n/2] = 2.0*ppha[p2n/2]; /*** Why double this? ***/
  xfreq[p2n/2] = sampling/2.0; /*** WYETH ??? ***/

  *rmod = mod; *rpha = pha; *rxfreq = xfreq; *rn = nn;

  myfree(p2data); myfree(cdata); myfree(pmod); myfree(ppha);
}
/**************************************-**************************************/
/*                                                                           */
/*                                   POWER                                   */
/*                                                                           */
/*   Takes same input file format as "plain_power".                          */
/*                                                                           */
/*   Compute the power spectrum from the fft in the obvious way.             */
/*                                                                           */
/*****************************************************************************/
void power(infile,outfile,duration,mean_flag,one_sided)
     char infile[];
     char outfile[];
     int duration;
     int mean_flag;
     int one_sided;
{
  FILE *fopen(),*fin,*fout;
  int i;
  int n,ns;
  float sampling,*data,*cdata,*pdata,mean,sdev;

  printf("  POWER\n");

  /* get the input data */
  fin = fopen(infile,"r");
  ns = fscanf(fin,"%d %f",&n,&sampling);

  if (duration>0)
    n = duration;
  
  if (power_of_two(n) == -1)
    exit_error("POWER","Data length must be integer power of two");

  data = (float *)myalloc(n*sizeof(float));
  for (i=0;i<n;i++)
    ns = fscanf(fin,"%f",&data[i]);
  fclose(fin);

  if (mean_flag){
    mean_sdev_farray(data,n,&mean,&sdev);
    printf("    Subtracting mean\n");
    for (i=0;i<n;i++)
      data[i] -= mean;
  }
  /* contort_real_farray(data,n); */
  
  cdata = make_complex_array(data,n);
  fft(cdata-1,n,1);
  pdata = power_of_complex(cdata,n);

  /* POWER (real) */
  fout = fopen(outfile,"w");
  if (one_sided){
    fprintf(fout,"%.4f %.3e\n",0.0,pdata[0]);
    for (i=1;i<n/2;i++)
      fprintf(fout,"%.4f %.3e\n",sampling/n*(float)i,
	      pdata[i]+pdata[n-i]);
    fprintf(fout,"%.4f %.3e\n",sampling/2.0,pdata[n/2]*2.0); /* ??? */
  }else{
    for (i=n/2;i<n;i++)
      fprintf(fout,"%.4f %.3e\n",sampling/n*(float)(i-n),pdata[i]);
    for (i=0;i<=n/2;i++)
      fprintf(fout,"%.4f %.3e\n",sampling/n*(float)i,pdata[i]);
  }
  fclose(fout);
  
  printf("  *** Output is ready for XPLOT.\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                                PLAIN_PWRWIN                               */
/*                                                                           */
/*   Compute the power spectrum (periodogram, see NumRecInC) using           */
/*   data windowing of a "plain" ascii data file.  The file format           */
/*   for input is as follows:                                                */
/*                                                                           */
/*      <n> <samples per second>                                             */
/*      <d1>                                                                 */
/*      <d2>                                                                 */
/*      ...                                                                  */
/*      <dn>                                                                 */
/*                                                                           */
/*   where "n" is the number of data points and the sampling rate            */
/*   is given in Hz, and may be a floating point number.                     */
/*                                                                           */
/*   The data points need not be on separate lines, just separated           */
/*   by white space.                                                         */
/*                                                                           */
/*****************************************************************************/
void plain_pwrwin(infile,outfile,duration,m,mean_flag)
     char infile[];
     char outfile[];
     int duration,m,mean_flag;
{
  FILE *fopen(),*fin,*fout;
  int i,k;
  int n,period,ns;
  float *data,*pdata,mean,sdev,sampling;

  printf("  PLAIN_PWRWIN\n");

  /* get the input data */
  fin = fopen(infile,"r");
  ns = fscanf(fin,"%d %f",&n,&sampling);
  data = (float *)myalloc(n*sizeof(float));
  for (i=0;i<n;i++)
    ns = fscanf(fin,"%f",&data[i]);
  fclose(fin);

  if (duration > 0)
    n = duration;
  
  if (mean_flag){
    mean_sdev_farray(data,n,&mean,&sdev);
    printf("    Subtracting mean\n");
    for (i=0;i<n;i++)
      data[i] -= mean;
  }
  k = (int)(((float)n/(float)m - 1.0)/2.0);
  period = (2*k+1)*m;
  
  printf("    Data windowing: m=%d k=%d n=%d\n",m,k,period);
  printf("    Read %d values.\n",n);

  pdata = (float *)myalloc((m+1)*sizeof(float));
  spctrm(data,pdata-1,m,k,1); /* send spike data, get power data */

  printf("    Normalizing by m*s, where s is sampling rate.\n");
  for (i=0;i<=m;i++)
    pdata[i] *= (float)m*sampling;

  printf("    Writing spectrum to %s.  Use \"pwrwin\" option in myplot.\n",
	 outfile);
  fout = fopen(outfile,"w");
  fprintf(fout,"%d %.4f\n",1,sampling);
  write_pwrwin_record(fout,pdata,infile,0,0,0,0,m,k);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                                PWRWIN_DATA                                */
/*                                                                           */
/*   Like the function above, but take an array of floats as input           */
/*   and return the spectrum.  Returned spectrum has m+1 entries.            */
/*                                                                           */
/*****************************************************************************/
float *pwrwin_data(data,sampling,n,m,mean_flag)
     float *data,sampling;
     int n,m,mean_flag;
{
  int i,k;
  float *tdata,*pdata;
  float mean,sdev;

  //printf("  PWRWIN_DATA\n");

  tdata = copy_farray(data,n);
  if (mean_flag){
    mean_sdev_farray(data,n,&mean,&sdev);
    //printf("    Subtracting mean\n");
    for (i=0;i<n;i++)
      tdata[i] -= mean;
  }
  k = (int)(((float)n/(float)m - 1.0)/2.0);

  /***
    int period;
    period = (2*k+1)*m;
    printf("    Data windowing: m=%d k=%d   Using %d of %d values.\n",
    m,k,period,n);
    ***/

  pdata = (float *)myalloc((m+1)*sizeof(float));
  spctrm(tdata,pdata-1,m,k,1); // send spike data, get power data
  myfree(tdata);

  //printf("    Normalizing by m*s, where s is sampling rate.\n");
  for (i=0;i<=m;i++)
    pdata[i] *= (float)m*sampling;
  //printf("    Norm --- done\n");

  return pdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                          AVG_PWRWIN_ROWS_2D_FARRAY                        */
/*                                                                           */
/*  Compute the average of the "pwrwin" spectrum for the rows in the         */
/*  2D farray.  The windowing parameter is "mm".                             */
/*                                                                           */
/*****************************************************************************/
void avg_pwrwin_rows_2d_farray(data,m,n,mm,sampling,pmean,psdev,pfreq)
     float **data;
     int m,n,mm;
     float sampling,**pmean,**psdev,**pfreq;
{
  int i;
  int k,period;
  float **pwr;

  printf("  AVG_PWRWIN_ROWS_2D_FARRAY\n");
  k = (int)(((float)n/(float)mm - 1.0)/2.0);
  period = (2*k+1)*mm;
  printf("    Data windowing: m=%d k=%d   Using %d of %d values.\n",
	 mm,k,period,n);

  pwr = (float **)myalloc(m*sizeof(float *));
  for(i=0;i<m;i++)
    pwr[i] = pwrwin_data(data[i],sampling,n,mm,0); /* mean_flag = 0 */

  mean_sdev_2d_farray(pwr,m,mm+1,pmean,psdev);
  for(i=0;i<m;i++)
    myfree(pwr[i]);
  myfree(pwr);

  (*pfreq) = (float *)myalloc((mm+1)*sizeof(float));
  for(i=0;i<=mm;i++)
    (*pfreq)[i] = (float)i/(float)mm * sampling/2.0;
}
/**************************************-**************************************/
/*                                                                           */
/*                       COMPUTE_RESPONSE_SINGLE_FROM_FT                     */
/*                                                                           */
/*  The precomputed FT of the filters is passed.                             */
/*  Apparently, "data" is modified here.                                     */
/*                                                                           */
/*****************************************************************************/
void compute_response_single_from_ft(data,f,f_s,xn,yn,tn,pflag,rr)
     float ***data,***f,**f_s;
     int xn,yn,tn;
     int pflag;
     float ****rr;
{
  float **speq,***r,**r_s,fc;

  if (pflag) printf("    COMPUTE_RESPONSE_SINGLE_FROM_FT\n");
  fc = 2.0/((float)(xn*yn*tn));

  if (pflag) printf("      Stimulus FFT."); fflush(stdout);
  speq = matrix(1,xn,1,2*yn);
  contort_real_3d_farray(data,1,1,1,xn,yn,tn);
  rlft3(data,speq,xn,yn,tn,1);

  if (pflag) printf("  Stim FT times filter FT."); fflush(stdout);
  three_d_fft_prod(data,speq,f,f_s,xn,yn,tn,&r,&r_s);

  free_matrix(speq,1,xn,1,2*yn);

  if (pflag) printf("  Inverse FFT.\n");
  rlft3(r,r_s,xn,yn,tn,-1); free_matrix(r_s,1,xn,1,2*yn);
  multiply_3d_farray(r,1,xn,1,yn,1,tn,fc);

  // WYETH - this is necessary to get correct answer from mfilter.c
  // which calls mod_dog_util's model_dog_01_get_response
  // REPLACE:  contort_real_3d_farray(r,1,1,1,xn,yn,tn); w/ BELOW:
  contort_real_3d_farray_space(r,1,1,1,xn,yn,tn);

  *rr = r;
}
/**************************************-**************************************/
/*                                                                           */
/*                      COMPUTE_RESPONSE_SINGLE_FROM_FT_2                    */
/*                                                                           */
/*  The precomputed FT of the filters is passed.                             */
/*  Apparently, "data" is modified here.                                     */
/*                                                                           */
/*  The 'data' and 'd_s' are already in FT form.                             */
/*                                                                           */
/*****************************************************************************/
void compute_response_single_from_ft_2(data,d_s,f,f_s,xn,yn,tn,pflag,rr)
     float ***data,**d_s,***f,**f_s;
     int xn,yn,tn;
     int pflag;
     float ****rr;
{
  float **speq,***r,**r_s,fc;

  if (pflag) printf("    COMPUTE_RESPONSE_SINGLE_FROM_FT_2\n");
  fc = 2.0/((float)(xn*yn*tn));

  if (pflag) printf("      Stim FT times filter FT."); fflush(stdout);
  three_d_fft_prod(data,d_s,f,f_s,xn,yn,tn,&r,&r_s);

  if (pflag) printf("  Inverse FFT.\n");
  rlft3(r,r_s,xn,yn,tn,-1); free_matrix(r_s,1,xn,1,2*yn);
  multiply_3d_farray(r,1,xn,1,yn,1,tn,fc);

  contort_real_3d_farray_space(r,1,1,1,xn,yn,tn);

  *rr = r;
}
/**************************************-**************************************/
/*                                                                           */
/*                     COMPUTE_RESPONSE_SINGLE_FROM_FT_RET_S                 */
/*                                                                           */
/*  The precomputed FT of the filters is passed.                             */
/*  Apparently, "data" is modified here.                                     */
/*                                                                           */
/*****************************************************************************/
float **compute_response_single_from_ft_ret_s(data,f,f_s,xn,yn,tn,pflag,rr)
     float ***data,***f,**f_s;
     int xn,yn,tn;
     int pflag;
     float ****rr;
{
  float **speq,***r,**r_s,fc;

  if (pflag) printf("    COMPUTE_RESPONSE_SINGLE_FROM_FT_RET_S\n");
  fc = 2.0/((float)(xn*yn*tn));

  if (pflag) printf("      Stimulus FFT."); fflush(stdout);
  speq = matrix(1,xn,1,2*yn);
  contort_real_3d_farray(data,1,1,1,xn,yn,tn);
  rlft3(data,speq,xn,yn,tn,1);

  if (pflag) printf("  Stim FT times filter FT."); fflush(stdout);
  three_d_fft_prod(data,speq,f,f_s,xn,yn,tn,&r,&r_s);

  if (pflag) printf("  Inverse FFT.\n");
  rlft3(r,r_s,xn,yn,tn,-1); free_matrix(r_s,1,xn,1,2*yn);
  multiply_3d_farray(r,1,xn,1,yn,1,tn,fc);

  contort_real_3d_farray_space(r,1,1,1,xn,yn,tn);

  *rr = r;

  return speq;
}
/**************************************-**************************************/
/*                                                                           */
/*                    OPTIMAL_LINEAR_FILTER_SPIKE_TRAIN                      */
/*                                                                           */
/*  Implement the equation:                                                  */
/*                                                                           */
/*            Sxs(w)                                                         */
/*    H(w) =  ------                                                         */
/*            Sxx(w)                                                         */
/*                                                                           */
/*  and ultimately determine h(t), where                                     */
/*                                                                           */
/*    Ssx(w) is the Fourier Transform (FT) of the cross-correlation          */
/*           of the stimulus and spike train                                 */
/*    Sxx(w) is the FT of the autocorrelation of the spike train,            */
/*           i.e., the power spectrum of the spike train.                    */
/*    H(w) is the FT of h(t), which is the optimal linear filter             */
/*         to convolve with the spike train to estimate the stim.            */
/*                                                                           */
/*****************************************************************************/
void optimal_linear_filter_spike_train(s,ns,start,period,nwin,stim,rht,rn)
     int *s,ns,start,period;
     int nwin;     /* Length of cross-corr to use */
     float *stim;
     float **rht;  /* Returned filter */
     int *rn;      /* Returned length */
{
  int ncc,ncc_work,n;
  float *sdata,*cc,*ccwork,*cdata,scount,oldmean;
  float *cc_real,*cc_imag,*ac_real,*ac_imag,*h_real,*h_imag,*hw,*ht;

  printf("  OPTIMAL_LINEAR_FILTER_SPIKE_TRAIN\n");
  n = period;
  ncc_work = 2*nwin;

  /***************************************************************************/
  /*** Compute cross-correlation of spikes and stimulus ***/
  scount = (float)count_spikes(s,ns,start,n);
  printf("    %d spikes\n",(int)scount);
  sdata = expand_spike_array(s,ns,start,n);
  subtract_mean_farray(sdata,n,&oldmean);
  /*append_farray_plot("zz.input","sdata",sdata,n,1);*/
  subtract_mean_farray(stim,n,&oldmean);
  cc = correlate(stim,sdata,n,n,&ncc); /*** CORRELATE ***/

  /***************************************************************************/
  /*** Compute the FT of a segment of this cross-correlation ***/
  ccwork = cc + (ncc-1)/2 - nwin;
  /*write_farray_plot("zz.cross",cc,ncc);
    write_farray_plot("zz.cross",ccwork,ncc_work);*/
  cdata = dft(ccwork,ncc_work);
  cc_real = real(cdata,ncc_work); /* Even function */
  cc_imag = imag(cdata,ncc_work); /* Even function */
  /*append_farray_plot("zz.dft","real",cc_real,ncc_work,1);
    append_farray_plot("zz.dft","imag",cc_imag,ncc_work,1);*/
  myfree(cc);

  /***************************************************************************/
  /*** Compute auto-correlation of spikes ***/
  cc = correlate(sdata,sdata,n,n,&ncc);
  ccwork = cc + (ncc-1)/2 - nwin;
  cdata = dft(ccwork,ncc_work);
  ac_real = real(cdata,ncc_work);
  ac_imag = imag(cdata,ncc_work);
  /*append_farray_plot("zz.spk.dft","real",ac_real,ncc_work,1);
    append_farray_plot("zz.spk.dft","imag",ac_imag,ncc_work,1);*/

  /***************************************************************************/
  /*** Divide in the Fourier domain, FT-1 to get h(t) ***/
  h_real = divide_farrays(cc_real,ac_real,ncc_work);
  h_imag = divide_farrays(cc_imag,ac_real,ncc_work);
  hw = make_complex_array_from_real_imag(h_real,h_imag,ncc_work);
  fft(hw-1,ncc_work,-1);

  ht = real(hw,ncc_work);
  contort_real_farray(ht,ncc_work);
  multiply_farray(ht,ncc_work,1.0/((float)ncc_work));
  /*write_farray_plot("zz.HT",ht,ncc_work);*/

  myfree(sdata);
  myfree(cc); myfree(cdata); myfree(cc_real); myfree(cc_imag);
  myfree(ac_real); myfree(ac_imag); myfree(hw); myfree(h_real); myfree(h_imag);

  *rn = ncc_work;
  *rht = ht;
}
/**************************************-**************************************/
/*                                                                           */
/*                                   AKSTIM                                  */
/*                                                                           */
/*   This is the function that needs to be called. allocates memory to       */
/*   store the stimulus in. tau, tseed, std1, std2, period, offset,          */
/*   floorScale, are parameters that define the stimlus. nfram is the        */
/*   length of the stimlus.                                                  */
/*                                                                           */
/*****************************************************************************/
float *akstim(tau,tseed,std1,std2,period,offset,floorp,scalep,nfram)
     int tau;
     int tseed;
     int std1;
     int std2;
     int period;
     int offset;
     int floorp;
     float scalep;
     int nfram;
{
  int i,j,k;
  int seed,nshuffle;
  int akstd1,akstd2,nrep,state;
  double tempd;
  float tempf;
  float q[8192];
  float r[8192];
  float ans[16384];
  float* stim = NULL;


  tseed = -tseed;
  if (tseed < 0)
    seed = tseed;
  else {
    printf("\n\nAK stimulus: Seed must be positive\n\n");
    return (float *)NULL;
  }

  akstd1 = std1;
  akstd2 = std2;

  /* This is the data series for steps in phase advance */
  nshuffle=3;
  /*  q = (float *)myalloc(2*nfram*sizeof(float)); */
  for (i=0;i<nshuffle;i++){  
    for (j=0;j<8192;j++)  
      q[j] = 20.0*(myrand_util_ran2(&seed))-10.0; /* Uniform from -10 to 10 */
  }


  /* This is the filter to smooth the data series */
  /* Data are stored in wrap around order  */
  for (j=0;j<2501;j++){              /* Indices 0-2500 are the filter */
    tempd=(double)j;
    r[j]=(float)exp(-(tempd)/tau);
  }
  for (j=5000;j>2500;j--)    /* Indices 2501-5000 are the wraparound filter */
    r[j]=r[-(j-5000)];
  for (j=5001;j<8192;j++)    /* Zero padding */
    r[j]=0.0;

  tempf=0.0;                 /* Normalize sum of coefficients to 1 */
  for (j=0;j<8192;j++)
    tempf=tempf+r[j];
  for (j=0;j<8192;j++) {
    r[j]=r[j]/tempf;
  }

    
  /*convlv(q,8192 ,r,5001,1,ans);*/
  convlv(q-1,8192 ,r-1,5001,1,ans-1);


  /*  for(j = 0; j < 100; j++) {
    if(j<100) mexPrintf("%3.2f", ans[j]);
  }
  */
  
  if (period==0){    /* For Brenner like stimulus, only one variance  */
    /*  Make a mean of 0 */
    tempf=0.0;

    /* Shift array by 1 because 1st entry contaminated */
    /**for (j=1;j<nfram;j++)
       tempf=tempf+ans[j];**/

    for (j=0;j<nfram;j++)
      tempf += ans[j];
    for (j=0;j<nfram;j++){
      ans[j] -= tempf/(float)nfram;
    }
    
    /*  Scale by std  */
    tempd=0.0;
    tempf=0.0;
    for (j=0;j<nfram;j++){
      tempd=(double)ans[j];
      tempf=(float)(tempf+pow(tempd,2.0));
    }
    tempf=tempf/nfram;
    tempd=(double)tempf;
    tempf=(float)pow(tempd,0.5);          /* std */
    for (j=0;j<nfram;j++)
      ans[j]=(akstd1*(ans[j]/tempf))+offset;
  }
  
  if (period!=0){    /* For Fairhall like stimulus where variance changes  */
    printf("Shouldn't see this if period is zero\n");
    nrep=nfram/period;
    state=0;

    /* Shift array by 1 because 1st entry contaminated */
    /*for (j=0;j<nfram;j++)
      ans[j]=ans[j+1];*/
    /* Subtract mean and normalize */
    for (j=0;j<nrep;j++){
      tempf=0.0;
      for (k=0;k<period;k++)
	tempf=tempf+ans[(j*period)+k];
      for (k=0;k<period;k++)
	ans[(j*period)+k]=ans[(j*period)+k]-(tempf/(period-1.0));
      tempd=0.0; 	/*  Scale by std  */
      tempf=0.0;
      for (k=0;k<period;k++){
	tempd=(double)ans[(j*period)+k];
	tempf=(float)(tempf+pow(tempd,2.0));
      }
      tempf=tempf/period;
      tempd=(double)tempf;
      tempf=(float)pow(tempd,0.5);          /* std */
      if (state==0){
	for (k=0;k<period;k++)
	  ans[(j*period)+k]=(akstd1*(ans[(j*period)+k]/tempf))+offset;
      }
      if (state==1){
	for (k=0;k<period;k++)
	  ans[(j*period)+k]=(akstd2*(ans[(j*period)+k]/tempf))+offset;
      }
      if (state==0)
	state=1;
      else if (state==1)
	state=0;
    }
  }

  stim = (float *)malloc(nfram*sizeof(float));

  if (floorp == 0) {
    for(i=0;i<nfram;i++){
      stim[i] = ans[i] * scalep;
    }
  }else{
    for(i=0;i<nfram;i++) {
      int temp;
      temp = floorf(ans[i]*scalep);
      stim[i] = temp;
    }
  }

  /*
    if(floorscale == 0)
    memcpy(stim,ans,sizeof(float)*nfram);
    else{
    for(i = 0; i < nfram; i++) {
    int temp;
    temp = ans[i];
    stim[i] = temp * floorscale;
    }
    }
  */

  return stim;
}
