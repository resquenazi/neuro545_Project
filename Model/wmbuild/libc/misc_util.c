/*****************************************************************************/
/*                                                                           */
/*  misc_util.c                                                              */
/*  wyeth bair                                                               */
/*  caltech                                                                  */
/*  11/20/93                                                                 */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "my_util.h"
#include "nr_util.h"

/**************************************-**************************************/
/*                                                                           */
/*                             GET_PATH_WITHOUT_NAME                         */
/*                                                                           */
/*  Return the part of "name" up to and including the last "/" in a new      */
/*  string.                                                                  */
/*                                                                           */
/*****************************************************************************/
char *get_path_without_name(name)
     char *name;
{
  int i;
  int nch;
  char *ptr,*p;

  ptr = strrchr(name,'/');
  if (ptr!=NULL){
    nch = (int)(ptr - name) + 1;
    p = (char *)myalloc((nch+1)*sizeof(char));  // +1 for '\0'
    for(i=0;i<nch;i++)
      p[i] = name[i];
    p[nch] = '\0';
  }else
    p = NULL;
  
  return p;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_NAME_WITHOUT_PATH                         */
/*                                                                           */
/*  Copy the part of "name" following the last "/" to "short_name".  If      */
/*  there is no "/", copy the entire name.                                   */
/*                                                                           */
/*****************************************************************************/
void get_name_without_path(name,short_name)
     char name[],short_name[];
{
  char *ptr;

  ptr = strrchr(name,'/');
  if (ptr!=NULL)
    strcpy(short_name,(ptr+1));
  else
    strcpy(short_name,name);
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_FILENAME_EXTENSION                        */
/*                                                                           */
/*  Copy the part of "name" following the last "." to "extension".  If       */
/*  there is no ".", set first char to NULL.                                 */
/*                                                                           */
/*****************************************************************************/
void get_filename_extension(name,extension)
     char name[],extension[];
{
  char *ptr;

  ptr = strrchr(name,'.');
  if (ptr!=NULL)
    strcpy(extension,(ptr+1));
  else
    extension[0] = '\0';
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_FILENAME_WITHOUT_EXTENSION                    */
/*                                                                           */
/*****************************************************************************/
void get_filename_without_extension(name,shortname)
     char name[],shortname[];
{
  char *ptr;

  strcpy(shortname,name);
  ptr = strrchr(shortname,'.');
  if (ptr!=NULL)
    *ptr = '\0'; // replace '.' with null.
}
/**************************************-**************************************/
/*                                                                           */
/*                            CUMULATIVE_BINOMIAL                            */
/*                                                                           */
/*  Return the probability of "k" or more successes in "n" trials            */
/*  with "p" probability of success per trial.                               */
/*                                                                           */
/*****************************************************************************/
float cumulative_binomial(k,n,p)
     int k,n;
     float p;
{
  float a,b,x;
  float t;

  a = (float)k;   // NumRecInC 2nd Ed p229
  b = (float)(n-k+1);
  x = p;
  t = mybetai(a,b,x);
  if (t==-1.0){
    printf("  *** k=%d n=%d p=%f\n",k,n,p);
    exit_error("CUMULATIVE_BINOMIAL","Error from mybetai");
  }
  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                                  BINOMIAL                                 */
/*                                                                           */
/*  Return the probabilility of 'k' successes in 'n' trials given prob p.    */
/*                                                                           */
/*****************************************************************************/
float binomial(k,n,p)
     int k,n;
     float p;
{
  float p1,p2;

  //printf("BINOMIAL\n");
  //printf("  k = %d   n= %d  p= %f\n",k,n,p);

  if (k <= 0)
    p1 = 1.0;
  else
    p1 = cumulative_binomial(k,n,p);

  if (k >= n)
    p2 = 0.0;
  else
    p2 = cumulative_binomial(k+1,n,p);

  //printf("  p1, p2  =  %f  %f      diff = %f \n",p1,p2,p1-p2);

  return p1 - p2;
}
/**************************************-**************************************/
/*                                                                           */
/*                             BINOMIAL_SIGNIF                               */
/*                                                                           */
/*  Return the probability of "k" or more successes in "n" trials            */
/*  with "p" probability of success per trial IF np <= k, otherwise          */
/*  return the probability of "k" or fewer successes.                        */
/*                                                                           */
/*****************************************************************************/
float binomial_signif(k,n,p)
     int k,n;
     float p;
{
  int nn;

  if (n>5000){
    /*printf("*** BINOMIAL_SIGNIF:  Changing %d/%d to ",k,n);*/
    nn = 5000 + my_rint(pow((float)(n-5000),0.8));
    k = my_rint((float)k * (float)nn/(float)n);
    n = nn;
    /*printf("%d/%d\n",k,n);*/
  }
  
  if ((float)k >= (float)n*p)
    return cumulative_binomial(k,n,p); /* k or more successes */
  else
    return cumulative_binomial(n-k,n,1.0-p); /* n-k or more failures */
}
/**************************************-**************************************/
/*                                                                           */
/*                            DISTANCE_POINT_PLANE                           */
/*                                                                           */
/*  Compute the distance of point (x,y,z) from the plane defined by          */
/*                                                                           */
/*     aX + bY + cZ + D = 0                                                  */
/*                                                                           */
/*****************************************************************************/
float distance_point_plane(x,y,z,a,b,c,d)
     float x,y,z,a,b,c,d;
{
  float dist;

  dist = fabs(a*x + b*y + c*z + d) / sqrt(a*a + b*b + c*c);

  return dist;
}
/**************************************-**************************************/
/*                                                                           */
/*                                PROJECTION_2D                              */
/*                                                                           */
/*  Return the projection of u on v, where u and v are 2D vectors.           */
/*                                                                           */
/*****************************************************************************/
void projection_2d(u1,u2,v1,v2,r1,r2)
     float u1,u2,v1,v2,*r1,*r2;
{
  float m,d;
  
  d = (v1*v1 + v2*v2);
  if (d==0.0){
    m = 0.0;
    exit_error("PROJECTION_2D","Vector has length 0");
  }else
    m = (u1*v1 + u2*v2)/d;

  *r1 = m*v1;
  *r2 = m*v2;
}
/**************************************-**************************************/
/*                                                                           */
/*                                DIFFERENCE_LOG                             */
/*                                                                           */
/*  Return the difference between the log of two numbers.  The base is       */
/*  specified by "base":                                                     */
/*    0 - natural log                                                        */
/*    2 - base 2                                                             */
/*   10 - base 10                                                            */
/*                                                                           */
/*****************************************************************************/
double difference_log(a,b,base)
     double a,b;
     int base;
{
  double d;

  if (base == 0)
    d = log(a) - log(b);
  else if (base == 2)
    d = my_log2(a) - my_log2(b);
  else if (base == 10)
    d = log10(a) - log10(b);
  else{
    d = 0.0;
    exit_error("LOG_DIFF","Unknown base");
  }

  return d;
}
/**************************************-**************************************/
/*                                                                           */
/*                              ANGLE_3D_VECTORS                             */
/*                                                                           */
/*  Return the angle between two 3D vectors.                                 */
/*                                                                           */
/*****************************************************************************/
float angle_3d_vectors(a1,a2,a3,b1,b2,b3)
     float a1,a2,a3;
     float b1,b2,b3;
{
  float dp,ms,theta;

  dp = a1*b1 + a2*b2 + a3*b3;

  ms = sqrt(a1*a1 + a2*a2 + a3*a3) * sqrt(b1*b1 + b2*b2 + b3*b3);

  if (ms == 0.0)
    exit_error("ANGLE_3D_VECTORS","Vector magnitude is zero"); 

  theta = 180.0 / M_PI * acos(dp / ms);

  return theta;
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_CIRCULAR_DIFF                           */
/*                                                                           */
/*  Return the difference between two numbers given on a circular scale      */
/*  from [0,d).                                                              */
/*                                                                           */
/*****************************************************************************/
float get_circular_diff(a,b,d)
     float a,b,d;
{
  float x;

  x = (a-b);
  if (x < 0.0)
    x *= -1;
  if (x > d/2.0)
    x = d-x;

  return x;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_SIGNED_CIRCULAR_DIFF                       */
/*                                                                           */
/*  Return the difference between two numbers given on a circular scale      */
/*  from [0,d).  Report the difference as positive if getting from b to      */
/*  a by the shortest route can be CCW, negative otherwise.                  */
/*                                                                           */
/*****************************************************************************/
float get_signed_circular_diff(a,b,d)
     float a,b,d;
{
  float x;

  x = (a-b);

  //  printf("a,b  %f %f      x  %f\n",a,b,x);


  if (x > d/2.0)
    x = -(d-x);
  else if (x < -d/2)
    x = d+x;

  //  printf("a,b  %f %f      x  %f\n",a,b,x);

  return x;
}
/**************************************-**************************************/
/*                                                                           */
/*                 GET_CIRCULAR_INTERVAL_LENGTH_CONTAINING_POINT             */
/*                                                                           */
/*  The points 'a' and 'b' divide the circle (of diameter 'd') into two      */
/*  intervals.  Return the length of the interval containing 'c'.            */
/*                                                                           */
/*****************************************************************************/
float get_circular_interval_length_containing_point(a,b,c,d)
     float a,b,c,d;
{
  float dab,dac,dbc,t;

  dab = get_circular_diff(a,b,d); // Distance between a and b
  dac = get_circular_diff(a,c,d); // a and c
  dbc = get_circular_diff(b,c,d); // b and c

  t = fabs((dac + dbc) - dab);

  if (t < 0.01*d)
    return dab;
  else
    return (d - dab);
}
/**************************************-**************************************/
/*                                                                           */
/*                                  SWAP_INT                                 */
/*                                                                           */
/*****************************************************************************/
void swap_int(a,b)
     int *a,*b;
{
  int t;
  
  t = *a;
  *a = *b;
  *b = t;
}  
/**************************************-**************************************/
/*                                                                           */
/*                                 SWAP_FLOAT                                */
/*                                                                           */
/*****************************************************************************/
void swap_float(a,b)
     float *a,*b;
{
  float t;

  t = *a;
  *a = *b;
  *b = t;
}  
/**************************************-**************************************/
/*                                                                           */
/*                                 WRAP_FLOAT                                */
/*                                                                           */
/*  Modify the value 'f' so that it is now wrapped around in the interval    */
/*  [0,size).                                                                */
/*                                                                           */
/*****************************************************************************/
void wrap_float(f,size)
     float *f,size;
{
  while(*f < 0.0)
    *f += size;

  while(*f >= size)
    *f -= size;
}  
/**************************************-**************************************/
/*                                                                           */
/*                                WRAP_DOUBLE                                */
/*                                                                           */
/*  Modify the value 'f' so that it is now wrapped around in the interval    */
/*  [0,size).                                                                */
/*                                                                           */
/*****************************************************************************/
void wrap_double(f,size)
     double *f,size;
{
  while(*f < 0.0)
    *f += size;

  while(*f >= size)
    *f -= size;
}  
/**************************************-**************************************/
/*                                                                           */
/*                                 DIFF_SDEV                                 */
/*                                                                           */
/*  Return the distance between a and b in units of "sdev".                  */
/*                                                                           */
/*****************************************************************************/
float diff_sdev(a,b,sdev)
     float a,b,sdev;
{
  if (a > b)
    return (a-b)/sdev;
  else
    return (b-a)/sdev;
}  
/**************************************-**************************************/
/*                                                                           */
/*                                 INTERP_TWO                                */
/*                                                                           */
/*  Return the y value for x0 on the line between (x1,y1) and (x2,y2).       */
/*                                                                           */
/*****************************************************************************/
float interp_two(x0,x1,x2,y1,y2)
     float x0,x1,x2,y1,y2;
{
  float f,y0;

  f = (x0-x1)/(x2-x1);
  y0 = y1 + f*(y2-y1);

  return y0;
}  
/**************************************-**************************************/
/*                                                                           */
/*                          MISC_UTIL_GET_ROUND_NUMBER                       */
/*                                                                           */
/*  Get a number that is close to the number but is a "round" number,        */
/*  meaning numbers like 1 2 5 10 20 50 etc, as opposed to 2.5, 7, 13 etc.   */
/*                                                                           */
/*****************************************************************************/
int misc_util_get_round_number(x)
     float x;
{
  int ten,y;
  float xt;

  xt = x;

  ten = 1;

  if (xt < 0.0){
    xt = -xt;
    ten = -1;
  }

  while(xt > 10.0){
    xt /= 10.0;
    ten *= 10;
  }

  if (xt < 1.5)
    y = 1;
  else if (xt < 3.5)
    y = 2;
  else if (xt < 7.5)
    y = 5;
  else
    y = 10;
  
  return y * ten;
}  
/**************************************-**************************************/
/*                                                                           */
/*                        DISPARITY_PHASE_SHIFT_ADJUST                       */
/*                                                                           */
/*  Used by 'mod_me_util'                                                    */
/*                                                                           */
/*****************************************************************************/
float disparity_phase_shift_adjust(aflag,theta,disp_off)
     int   aflag;         // 0 - No adjustment
                          // 1 - correct for direction
                          // 2 - partial correction - Should remove ?
     float theta;         // direction of motion
     float disp_off;      // phase offset for 0 deg (horiz. right) motion
{
  float phd,twrap;

  phd = 0.0;

  if (aflag >= 1){
    //
    //  If we want to make these filters be offset for disparity (from some
    //  other set of filters, not among these), then we have to shift 
    //  oppositely directed filters in opposite directions.
    //
    if (disp_off != 0.0){  // If we want a phase offset
      twrap = theta;
      wrap_float(&twrap,360.0);  // Make sure 'twrap' is in [ 0,360.0 )

      if (aflag == 2){
	if ((twrap > 90.0) && (twrap <= 270.0)){
	  phd = -disp_off;
	}else{
	  phd =  disp_off;  // This range includes direction (twrap) = 0
	}
      }else{  // if (aflag == 1){
	phd = disp_off * cos(theta * M_PI/180.0);
      }
    }
  }

  return phd;
}
