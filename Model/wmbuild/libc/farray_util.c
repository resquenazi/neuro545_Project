/*****************************************************************************/
/*                                                                           */
/*  farray_util.c                                                            */
/*  wyeth bair                                                               */
/*  caltech                                                                  */
/*  02/23/93                                                                 */
/*                                                                           */
/*  These routines do not depend on "iarray_util.c" or "carray_util.c".      */
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
#include "misc_util.h"

/**************************************-**************************************/
/*                                                                           */
/*                               PRINT_FARRAY                                */
/*                                                                           */
/*****************************************************************************/
void print_farray(d,n,name)
     float *d;    // [n] data
     int n;       // length of array
     char *name;  // name of array
{
  int i;
  float **fdata;

  printf("  PRINT_FARRAY\n");
  printf("    %s %d\n",name,n);
  for(i=0;i<n;i++)
    printf(" %f",d[i]);
  printf("\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_2D_FARRAY                               */
/*                                                                           */
/*  GET_POINTER_FARRAY if n = 0, second dim points to NULL.                  */
/*                                                                           */
/*****************************************************************************/
float **get_2d_farray(m,n)
     int m,n;
{
  int i;
  float **fdata;
  
  fdata = (float **)myalloc(m*sizeof(float *));
  for(i=0;i<m;i++)
    fdata[i] = (float *)myalloc(n*sizeof(float));  // NULL ptr if n=0
  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_2D_DARRAY                               */
/*                                                                           */
/*  GET_POINTER_FARRAY if n = 0, second dim points to NULL.                  */
/*                                                                           */
/*****************************************************************************/
double **get_2d_darray(m,n)
     int m,n;
{
  int i;
  double **ddata;
  
  ddata = (double **)myalloc(m*sizeof(double *));
  for(i=0;i<m;i++)
    ddata[i] = (double *)myalloc(n*sizeof(double));  // NULL ptr if n=0
  return ddata;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_CONST_2D_FARRAY                           */
/*                                                                           */
/*****************************************************************************/
float **get_const_2d_farray(m,n,c)
     int m,n;
     float c;
{
  int i,j;
  float **fdata;
  
  fdata = (float **)myalloc(m*sizeof(float *));
  for(i=0;i<m;i++){
    fdata[i] = (float *)myalloc(n*sizeof(float));
    for(j=0;j<n;j++)
      fdata[i][j] = c;
  }
  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_ZERO_2D_FARRAY                            */
/*                                                                           */
/*****************************************************************************/
float **get_zero_2d_farray(m,n)
     int m,n;
{
  return get_const_2d_farray(m,n,0.0);
}
/**************************************-**************************************/
/*                                                                           */
/*                               FREE_2D_FARRAY                              */
/*                                                                           */
/*****************************************************************************/
void free_2d_farray(data,n)
     float **data;
     int n;
{
  int i;

  for(i=0;i<n;i++){
    if (data[i] != NULL)
      myfree(data[i]);
    else{
      exit_error("FREE_2D_FARRAY","Found null row");
    }
  }
  if (data != NULL)
    myfree(data);
  else
    exit_error("FREE_2D_FARRAY","Found null data");
}
/**************************************-**************************************/
/*                                                                           */
/*                          FREE_2D_FARRAY_IGNORE_NULL                       */
/*                                                                           */
/*****************************************************************************/
void free_2d_farray_ignore_null(data,n)
     float **data;
     int n;
{
  int i;

  for(i=0;i<n;i++){
    if (data[i] != NULL)
      myfree(data[i]);
  }
  if (data != NULL)
    myfree(data);
  else
    exit_error("FREE_2D_FARRAY","Found null data");
}
/**************************************-**************************************/
/*                                                                           */
/*                                COPY_FARRAY                                */
/*                                                                           */
/*  Return a copy of a float array.                                          */
/*                                                                           */
/*****************************************************************************/
float *copy_farray(data,n)
     float *data;
     int n;
{
  int i;
  float *cdata;
  
  if (n < 0)
    exit_error("COPY_FARRY","n < 0");
  
  cdata = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    cdata[i] = data[i];
  return cdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                                COPY_DARRAY                                */
/*                                                                           */
/*  Return a copy of a float array.                                          */
/*                                                                           */
/*****************************************************************************/
double *copy_darray(data,n)
     double *data;
     int n;
{
  int i;
  double *cdata;
  
  if (n < 0)
    exit_error("COPY_FARRY","n < 0");
  
  cdata = (double *)myalloc(n*sizeof(double));
  for(i=0;i<n;i++)
    cdata[i] = data[i];
  return cdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                           UNSIGNED_CHAR_TO_FARRAY                         */
/*                                                                           */
/*****************************************************************************/
float *unsigned_char_to_farray(data,n)
     unsigned char *data;
     int n;
{
  int i;
  float *fdata;
  
  fdata = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    fdata[i] = (float)data[i];
  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                                  I2FARRAY                                 */
/*                                                                           */
/*  Return a floating point array containing the iarray data.                */
/*                                                                           */
/*****************************************************************************/
float *i2farray(data,n)
     int *data,n;
{
  int i;
  float *fdata;

  fdata = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    fdata[i] = (float)data[i];
  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                                  D2FARRAY                                 */
/*                                                                           */
/*  Return a floating point array containing the darray data.                */
/*                                                                           */
/*****************************************************************************/
float *d2farray(data,n)
     double *data;
     int n;
{
  int i;
  float *fdata;

  fdata = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    fdata[i] = (float)data[i];
  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 I2FARRAY_2D                               */
/*                                                                           */
/*  Return a floating point array containing the iarray data.                */
/*                                                                           */
/*****************************************************************************/
float **i2farray_2d(data,n,m)
     int **data,n,m;
{
  int i,j;
  float **fdata;

  fdata = (float **)myalloc(n*sizeof(float *));
  for(i=0;i<n;i++){
    fdata[i] = (float *)myalloc(m*sizeof(float));
    for(j=0;j<m;j++)
      fdata[i][j] = (float)data[i][j];
  }
  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                               REVERSE_FARRAY                              */
/*                                                                           */
/*  Reverse the values in the farray, i.e., make the last one first.         */
/*                                                                           */
/*****************************************************************************/
void reverse_farray(data,n)
     float *data;
     int n;
{
  int i;
  float t;
  
  for(i=0;i<n/2;i++){
    t = data[i];
    data[i] = data[n-1-i];
    data[n-1-i] = t;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_REGULAR_SAMPLE_FARRAY                        */
/*                                                                           */
/*  Return an array with n equally spaced values from a to b.                */
/*                                                                           */
/*****************************************************************************/
float *get_regular_sample_farray(a,b,n)
     float a,b;
     int n;
{
  int i;
  float *data,d;

  data = (float *)myalloc(n*sizeof(float));
  d = b-a;
  for(i=0;i<n;i++)
    data[i] = (float)i/(float)(n-1)*d + a;

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                               COPY_2D_FARRAY                              */
/*                                                                           */
/*****************************************************************************/
float **copy_2d_farray(data,m,n)
     float **data;
     int n,m;
{
  int i;
  float **cdata;
  
  cdata = (float **)myalloc(m*sizeof(float *));
  for(i=0;i<m;i++)
    cdata[i] = copy_farray(data[i],n);
  return cdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_TRANSPOSE_2D_FARRAY                        */
/*                                                                           */
/*****************************************************************************/
float **get_transpose_2d_farray(data,m,n)
     float **data;
     int n,m;
{
  int i,j;
  float **t;

  t = get_2d_farray(n,m);
  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      t[i][j] = data[j][i];

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                             TRANSPOSE_2D_FARRAY                           */
/*                                                                           */
/*  This free the original storage, mallocs new storage, and may change the  */
/*  dimensions of the array.                                                 */
/*                                                                           */
/*****************************************************************************/
void transpose_2d_farray(rdata,m,n)
     float ***rdata;
     int n,m;
{
  float **t;

  t = get_transpose_2d_farray(*rdata,m,n);
  free_2d_farray(*rdata,m);
  *rdata = t;
}
/**************************************-**************************************/
/*                                                                           */
/*                               COPY_3D_FARRAY                              */
/*                                                                           */
/*****************************************************************************/
float ***copy_3d_farray(data,x0,xn,y0,yn,z0,zn)
     float ***data;
     int x0,xn,y0,yn,z0,zn;
{
  int i,j;
  float ***cdata;
  
  cdata = (float ***)myalloc(xn*sizeof(float **));
  for(i=0;i<xn;i++){
    cdata[i] = (float **)myalloc(yn*sizeof(float *));
    for(j=0;j<yn;j++)
      cdata[i][j] = copy_farray(data[x0+i][y0+j]+z0,zn);
  }
  return cdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_COLUMN_2D_FARRAY                           */
/*                                                                           */
/*  Return a sub-array of the 2d array.                                      */
/*                                                                           */
/*****************************************************************************/
float *get_column_2d_farray(data,xn,yn,col)
     float **data;
     int xn,yn,col;
{
  int i;
  float *tdata;

  tdata = (float *)myalloc(xn*sizeof(float));
  for(i=0;i<xn;i++)
    tdata[i] = data[i][col];

  return tdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_COLUMN_2D_DOUBLE_FARRAY                       */
/*                                                                           */
/*  Return a sub-array of the 2d double array.                               */
/*                                                                           */
/*****************************************************************************/
float *get_column_2d_double_farray(data,xn,yn,col)
     double **data;
     int xn,yn,col;
{
  int i;
  float *tdata;

  tdata = (float *)myalloc(xn*sizeof(float));
  for(i=0;i<xn;i++)
    tdata[i] = (float)data[i][col];
  return tdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_SUB_2D_FARRAY                             */
/*                                                                           */
/*  Return a sub-array of the 2d array.                                      */
/*                                                                           */
/*****************************************************************************/
void get_sub_2d_farray(data,xn,yn,rdata,x0,xd,y0,yd)
     float **data;
     int xn,yn;
     float ***rdata;
     int x0,xd,y0,yd; // start and duration in each dimension
{
  int i,j;
  float **tdata;

  tdata = get_2d_farray(xd,yd);
  for(i=0;i<xd;i++)
    for(j=0;j<yd;j++)
      tdata[i][j] = data[x0+i][y0+j];

  *rdata = tdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                                GET_3D_FARRAY                              */
/*                                                                           */
/*  Can be used as 'get_2d_pointer_farray' by setting d3 to 0.               */
/*                                                                           */
/*****************************************************************************/
float ***get_3d_farray(d1,d2,d3)
     int d1,d2,d3;
{
  int i1,i2;
  float ***data;
  
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
/*                             GET_ZERO_3D_FARRAY                            */
/*                                                                           */
/*****************************************************************************/
float ***get_zero_3d_farray(m,n,o)
     int m,n,o;
{
  int i,j;
  float ***fdata;
  
  fdata = (float ***)myalloc(m*sizeof(float **));
  for(i=0;i<m;i++){
    fdata[i] = (float **)myalloc(n*sizeof(float *));
    for(j=0;j<n;j++)
      fdata[i][j] = get_zero_farray(o);
  }
  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                               ZERO_3D_FARRAY                              */
/*                                                                           */
/*****************************************************************************/
void zero_3d_farray(fdata,m,n,o)
     float ***fdata;
     int m,n,o;
{
  int i,j,k;
  
  for(i=0;i<m;i++)
    for(j=0;j<n;j++)
      for(k=0;k<o;k++)
	fdata[i][j][k] = 0.0;
}
/**************************************-**************************************/
/*                                                                           */
/*                              FREE_3D_FARRAY                               */
/*                                                                           */
/*****************************************************************************/
void free_3d_farray(fdata,m,n,o)
     float ***fdata;
     int m,n,o;
{
  int i,j;
  
  for(i=0;i<m;i++){
    for(j=0;j<n;j++)
      myfree(fdata[i][j]);
    myfree(fdata[i]);
  }
  myfree(fdata);
}
/**************************************-**************************************/
/*                                                                           */
/*                             FREE_3D_FARRAY_NULL                           */
/*                                                                           */
/*  Allows the third dimension to contain NULL entries.                      */
/*  AKA:   free_2d_pointer_farray                                            */
/*                                                                           */
/*****************************************************************************/
void free_3d_farray_null(fdata,m,n,o)
     float ***fdata;
     int m,n,o;           // 'o' is not used currently
{
  int i,j;
  
  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      if (fdata[i][j] != NULL)
	myfree(fdata[i][j]);
    }
    myfree(fdata[i]);
  }
  myfree(fdata);
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_XYT_FROM_TXY_3D_FARRAY                        */
/*                                                                           */
/*****************************************************************************/
float ***get_xyt_from_txy_3d_farray(data,n1,n2,n3)
     float ***data;
     int n1,n2,n3;
{
  int i,j,k;
  float ***t;

  t = get_3d_farray(n2,n3,n1);
  for(i=0;i<n1;i++)
    for(j=0;j<n2;j++)
      for(k=0;k<n3;k++)
	t[j][k][i] = data[i][j][k];

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_TXY_FROM_XYT_3D_FARRAY                        */
/*                                                                           */
/*****************************************************************************/
float ***get_txy_from_xyt_3d_farray(data,n1,n2,n3)
     float ***data;
     int n1,n2,n3;
{
  int i,j,k;
  float ***t;

  t = get_3d_farray(n3,n1,n2);
  for(i=0;i<n1;i++)
    for(j=0;j<n2;j++)
      for(k=0;k<n3;k++)
	t[k][i][j] = data[i][j][k];

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                              GET_SUB_3D_FARRAY                            */
/*                                                                           */
/*  WYETH - SAME AS COPY_3D_FARRAY above ????                                */
/*  WYETH - SAME AS COPY_3D_FARRAY above ????                                */
/*  WYETH - SAME AS COPY_3D_FARRAY above ????                                */
/*                                                                           */
/*  Return a sub-array of the 3d array.                                      */
/*                                                                           */
/*****************************************************************************/
void get_sub_3d_farray(data,xn,yn,tn,rdata,x0,xd,y0,yd,t0,td)
     float ***data;
     int xn,yn,tn;
     float ****rdata;
     int x0,xd,y0,yd,t0,td; // start and duration in each dimension
{
  int i,j,k;
  float ***tdata;

  tdata = get_3d_farray(xd,yd,td);
  for(i=0;i<xd;i++)
    for(j=0;j<yd;j++)
      for(k=0;k<td;k++)
	tdata[i][j][k] = data[x0+i][y0+j][t0+k];

  *rdata = tdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_1D_FROM_3D_FARRAY                          */
/*                                                                           */
/*  Return an array spanning one dimension of the 3D array.                  */
/*                                                                           */
/*  Set 'sx' and 'nx' (x=1,2,3) to the start index and number of points to   */
/*  extract for the appropriate dimension, x.  Other nx values should be 1   */
/*  and other sx values should specify coordinates in the plane orthogonal   */
/*  to the extracted 1d array.                                               */
/*                                                                           */
/*****************************************************************************/
float *get_1d_from_3d_farray(data,s1,n1,s2,n2,s3,n3)
     float ***data;
     int s1,n1,s2,n2,s3,n3;
{
  int i;
  float *t;

  if (n1==1 && n2==1){
    t = copy_farray(&(data[s1][s2][s3]),n3);
  }else if (n1==1 && n3==1){
    t = (float *)myalloc(n2*sizeof(float));
    for(i=0;i<n2;i++)
      t[i] = data[s1][s2+i][s3];
  }else if (n2==1 && n3==1){
    t = (float *)myalloc(n1*sizeof(float));
    for(i=0;i<n1;i++)
      t[i] = data[s1+i][s2][s3];
  }else{
    t = NULL;
    exit_error("GET_1D_FROM_3D_FARRAY","Invalid 1d array specification");
  }

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_2D_FROM_3D_FARRAY                           */
/*                                                                           */
/*  Return a sub-array of the 3d array.                                      */
/*                                                                           */
/*****************************************************************************/
float **get_2d_from_3d_farray(data,x0,xn,y0,yn,k)
     float ***data;
     int x0,xn,y0,yn,k;
{
  int i,j;
  float **tdata;

  tdata = get_2d_farray(xn,yn);
  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++)
      tdata[i][j] = data[i+x0][j+y0][k];

  return tdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_4D_FARRAY                               */
/*                                                                           */
/*  Set the nth index to -1 to avoid malloc'ing storage.                     */
/*                                                                           */
/*****************************************************************************/
float ****get_4d_farray(d1,d2,d3,d4)
     int d1,d2,d3,d4;
{
  int i1,i2,i3;
  float ****fdata;

// *** WYETH - consider rewriting so that unallocated pointers are set to NULL
// *** WYETH - consider rewriting so that unallocated pointers are set to NULL
// *** WYETH - consider rewriting so that unallocated pointers are set to NULL

  fdata = (float ****)myalloc(d1*sizeof(float ***));
  if (d2 >= 0)
    for(i1=0;i1<d1;i1++){
      fdata[i1] = (float ***)myalloc(d2*sizeof(float **));
      if (d3 >= 0)
	for(i2=0;i2<d2;i2++){
	  fdata[i1][i2] = (float **)myalloc(d3*sizeof(float *));
	  if (d4 >= 0){ // *** WYETH myalloc gives NULL if d4==0 anyway
	                // *** thus there is redundancy/confusion here?
	    for(i3=0;i3<d3;i3++){
	      fdata[i1][i2][i3] = (float *)myalloc(d4*sizeof(float));
	    }
	  }else{
	    for(i3=0;i3<d3;i3++){
	      fdata[i1][i2][i3] = NULL;
	    }
	  }
	}
    }
  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_CONST_4D_FARRAY                           */
/*                                                                           */
/*****************************************************************************/
float ****get_const_4d_farray(d1,d2,d3,d4,c)
     int d1,d2,d3,d4;
     float c;
{
  int i1,i2;
  float ****fdata;

  fdata = (float ****)myalloc(d1*sizeof(float ***));
  for(i1=0;i1<d1;i1++){
    fdata[i1] = (float ***)myalloc(d2*sizeof(float **));
    for(i2=0;i2<d2;i2++){
      fdata[i1][i2] = get_const_2d_farray(d3,d4,c);
    }
  }
  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                              READ_4D_FARRAY                               */
/*                                                                           */
/*  I wrote this as an example for Timothy Ma.  Have not yet used it         */
/*  on Nov 2, 2017.                                                          */
/*                                                                           */
/*****************************************************************************/
float ****read_4d_farray(infile,rn1,rn2,rn3,rn4)
     char *infile;
     int *rn1,*rn2,*rn3,*rn4;
{
  FILE *fopen(),*fin;
  int i,j,k,l;
  int n1,n2,n3,n4,nn;
  float ****fdata;

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** Cannot open file %s.\n",infile);
    exit(0);
  }
    
  nn = fscanf(fin,"%d %d %d %d",&n1,&n2,&n3,&n4);

  fdata = (float ****)myalloc(n1*sizeof(float ***));
  for(i=0;i<n1;i++){
    fdata[i] = (float ***)myalloc(n2*sizeof(float **));
    for(j=0;j<n2;j++){
      fdata[i][j] = (float **)myalloc(n3*sizeof(float *));
      for(k=0;k<n3;k++){
	fdata[i][j][k] = (float *)myalloc(n4*sizeof(float));
      }
    }
  }

  for(i=0;i<n1;i++)
    for(j=0;j<n2;j++)
      for(k=0;k<n3;k++)
	for(l=0;l<n4;l++)
	  nn = fscanf(fin,"%f",&(fdata[i][j][k][l]));
    
  fclose(fin);

  *rn1 = n1;
  *rn2 = n2;
  *rn3 = n3;
  *rn4 = n4;

  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                           READ_3D_VARIABLE_FARRAY                         */
/*                                                                           */
/*  I wrote this as an example for Timothy Ma.  Have not yet used it         */
/*  on Nov 2, 2017.                                                          */
/*                                                                           */
/*****************************************************************************/
float ***read_3d_variable_farray(infile,rn1,rn2,rcnt)
     char *infile;   // Input file name
     int *rn1,*rn2;  // (pointers to) Lengths of first two dimensions
     int ***rcnt;    // (pointer to) 2D array of lengths [*rn1][*rn2]
{
  FILE *fopen(),*fin;
  int i,j,k;
  int n1,n2,n3,**cnt,nn;
  float ***fdata;

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** Cannot open file %s.\n",infile);
    exit(0);
  }
    
  nn = fscanf(fin,"%d %d",&n1,&n2);  // Header has first two dimensions

  //
  //  Allocate storage for the 2D int array called 'cnt'
  //  Think of this as:  cnt[n1][n2]
  //
  cnt = (int **)myalloc(n1*sizeof(int *));
  for(i=0;i<n1;i++)
    cnt[i] = (int *)myalloc(n2*sizeof(int));

  //
  //  Now create the 2D array of float POINTERS
  //
  fdata = (float ***)myalloc(n1*sizeof(float **));
  for(i=0;i<n1;i++)
    fdata[i] = (float **)myalloc(n2*sizeof(float *));


  //
  //  Now read the data from the data file
  //
  for(i=0;i<n1;i++){
    for(j=0;j<n2;j++){  // For each position in the 2D matrix

      nn = fscanf(fin,"%d",&n3);  // Header has first two dimensions
      fdata[i][j] = (float *)myalloc(n3*sizeof(float));
      for(k=0;k<n3;k++){
	nn = fscanf(fin,"%f",&(fdata[i][j][k]));
      }
      cnt[i][j] = n3;
    }
  }

  fclose(fin);

  *rn1 = n1;
  *rn2 = n2;
  *rcnt = cnt;

  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                              FREE_4D_FARRAY                               */
/*                                                                           */
/*****************************************************************************/
void free_4d_farray(fdata,d1,d2,d3,d4)
     float ****fdata;
     int d1,d2,d3,d4;
{
  int i1,i2,i3;

  for(i1=0;i1<d1;i1++){
    for(i2=0;i2<d2;i2++){
      for(i3=0;i3<d3;i3++)
	myfree(fdata[i1][i2][i3]);
      myfree(fdata[i1][i2]);
    }
    myfree(fdata[i1]);
  }
  myfree(fdata);
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_5D_FARRAY                               */
/*                                                                           */
/*  Set the nth index to -1 to avoid malloc'ing storage.                     */
/*                                                                           */
/*****************************************************************************/
float *****get_5d_farray(d1,d2,d3,d4,d5)
     int d1,d2,d3,d4,d5;
{
  int i1,i2,i3,i4;
  float *****fdata;

  fdata = (float *****)myalloc(d1*sizeof(float ****));
  if (d2 >= 0)
    for(i1=0;i1<d1;i1++){
      fdata[i1] = (float ****)myalloc(d2*sizeof(float ***));
      if (d3 >= 0)
	for(i2=0;i2<d2;i2++){
	  fdata[i1][i2] = (float ***)myalloc(d3*sizeof(float **));
	  if (d4 >= 0)
	    for(i3=0;i3<d3;i3++){
	      fdata[i1][i2][i3] = (float **)myalloc(d4*sizeof(float *));
	      if (d5 >= 0)
		for(i4=0;i4<d4;i4++)
		  fdata[i1][i2][i3][i4] = (float *)myalloc(d5*sizeof(float));
	    }
	}
    }
  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                              FREE_5D_FARRAY                               */
/*                                                                           */
/*****************************************************************************/
void free_5d_farray(fdata,d1,d2,d3,d4,d5)
     float *****fdata;
     int d1,d2,d3,d4,d5;
{
  int i1,i2,i3,i4;

  for(i1=0;i1<d1;i1++){
    for(i2=0;i2<d2;i2++){
      for(i3=0;i3<d3;i3++){
	for(i4=0;i4<d4;i4++)
	  myfree(fdata[i1][i2][i3][i4]);
	myfree(fdata[i1][i2][i3]);
      }
      myfree(fdata[i1][i2]);
    }
    myfree(fdata[i1]);
  }
  myfree(fdata);
}
/**************************************-**************************************/
/*                                                                           */
/*                               EQUAL_FARRAYS                               */
/*                                                                           */
/*****************************************************************************/
int equal_farrays(data1,data2,n)
     float *data1,*data2;
     int n;
{
  int k;
  int flag;

  flag = 1;
  k = 0;
  while((flag == 1)&&(k < n)){
    if (data1[k] != data2[k])
      flag = 0;
    else
      k += 1;
  }
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                                ZERO_FARRAY                                */
/*                                                                           */
/*****************************************************************************/
void zero_farray(data,n)
     float *data;
     int n;
{
  int i;
  
  for(i=0;i<n;i++)
    data[i] = 0.0;
}
/**************************************-**************************************/
/*                                                                           */
/*                            ROUND_TO_INT_FARRAY                            */
/*                                                                           */
/*****************************************************************************/
void round_to_int_farray(data,n)
     float *data;
     int n;
{
  int i;
  
  for(i=0;i<n;i++)
    data[i] = (float)my_rint(data[i]);
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_ROUND_TO_INT_FARRAY                        */
/*                                                                           */
/*  Return an iarray containing the float data, rounding to the nearest      */
/*  integer.                                                                 */
/*                                                                           */
/*****************************************************************************/
int *get_round_to_int_farray(data,n)
     float *data;
     int n;
{
  int i;
  int *idata;
  
  idata = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++)
    idata[i] = my_rint(data[i]);
  return idata;
}
/**************************************-**************************************/
/*                                                                           */
/*                             ROTATION_MATRIX_3D                            */
/*                                                                           */
/*  Return the 3x3 rotation matrix to rotate by 'theta' around vector xyz.   */
/*                                                                           */
/*****************************************************************************/
float **rotation_matrix_3d(theta,vx,vy,vz)
     float theta;         // Angle (deg)
     float vx,vy,vz;      // vector
{
  float **m,x,y,z,d,a;

  d = sqrt(vx*vx + vy*vy + vz*vz);  // Length of vector

  x = vx/d;  // Unit vector
  y = vy/d;
  z = vz/d;

  a = M_PI/180.0 * theta;  // Angle (radians)

  m = get_2d_farray(3,3);

  m[0][0] = 1.0 + (1.0-cos(a)) * (x*x - 1.0);
  m[0][1] = -z*sin(a) + (1.0-cos(a))*x*y;
  m[0][2] =  y*sin(a) + (1.0-cos(a))*x*z;

  m[1][0] =  z*sin(a) + (1.0-cos(a))*x*y;
  m[1][1] = 1.0 + (1.0-cos(a)) * (y*y - 1.0);
  m[1][2] = -x*sin(a) + (1.0-cos(a))*y*z;

  m[2][0] = -y*sin(a) + (1.0-cos(a))*x*z;
  m[2][1] =  x*sin(a) + (1.0-cos(a))*y*z;
  m[2][2] = 1.0 + (1.0-cos(a)) * (z*z - 1.0);

  return m;
}
/**************************************-**************************************/
/*                                                                           */
/*                            FARRAY_ROTATE_COORDS                           */
/*                                                                           */
/*  Rotation matrix is applied to all coordinates.                           */
/*                                                                           */
/*****************************************************************************/
void farray_rotate_coords(theta,x,y,n)
     float theta;
     float *x,*y;
     int n;
{
  int i;
  float tx,ty;
  float a,b,c,d;

  a =  cos(M_PI/180.0*theta); b = sin(M_PI/180.0*theta); // rotation matrix
  c = -sin(M_PI/180.0*theta); d = cos(M_PI/180.0*theta);

  for(i=0;i<n;i++){
    tx = x[i];
    ty = y[i];
    x[i] = a*tx + c*ty;
    y[i] = b*tx + d*ty;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              ROTATE_FOUR_PAIRS                            */
/*                                                                           */
/*  Rotation matrix is applied to four pairs of x,y coordinates.             */
/*                                                                           */
/*****************************************************************************/
void rotate_four_pairs(theta,px0,py0,px1,py1,px2,py2,px3,py3)
     float theta;
     float *px0,*py0,*px1,*py1,*px2,*py2,*px3,*py3;
{
  float x,y;
  float a,b,c,d;

  a =  cos(M_PI/180.0*theta); b = sin(M_PI/180.0*theta); // rotation matrix
  c = -sin(M_PI/180.0*theta); d = cos(M_PI/180.0*theta);

  x = *px0;
  y = *py0;
  *px0 = a*x + c*y;
  *py0 = b*x + d*y;

  x = *px1;
  y = *py1;
  *px1 = a*x + c*y;
  *py1 = b*x + d*y;

  x = *px2;
  y = *py2;
  *px2 = a*x + c*y;
  *py2 = b*x + d*y;

  x = *px3;
  y = *py3;
  *px3 = a*x + c*y;
  *py3 = b*x + d*y;
}
/**************************************-**************************************/
/*                                                                           */
/*                                ROTATE_FARRAY                              */
/*                                                                           */
/*  Shift all elements in the array to the right (to higher indices) by      */
/*  k indices;                                                               */
/*                                                                           */
/*****************************************************************************/
void rotate_farray(data,n,k)
     float *data;
     int n,k;
{
  int i,j;
  float *t;

  while(k<0) // Change negative rotations into equivalent positive ones.
    k += n;
  while(k>=n)
    k -= n;

  if ((n > 0)&&(k != 0)){
    t = (float *)myalloc(n*sizeof(float));
    for(i=0;i<n;i++){
      j = (i+k)%n;
      t[j] = data[i];
    }
    for(i=0;i<n;i++)
      data[i] = t[i];
    myfree(t);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                 SHIFT_FARRAY                              */
/*                                                                           */
/*  Shift all elements in the array to the right (to higher indices) by      */
/*  k indices.  k < 0 indicates leftward shifts.                             */
/*                                                                           */
/*****************************************************************************/
void shift_farray(data,n,k,fillvalue)
     float *data;
     int n,k;
     float fillvalue;
{
  int i;

  if (k > n)
    k = n;
  if (k < -n)
    k = -n;

  if (n > 0){
    if (k > 0){
      for(i=n-1;i>=k;i--)
	data[i] = data[i-k];
      for(i=0;i<k;i++)
	data[i] = fillvalue;
    }else if (k < 0){
      k = -k;
      for(i=0;i<(n-k);i++)
	data[i] = data[i+k];
      for(i=n-k;i<n;i++)
	data[i] = fillvalue;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_Y_AT_X_FARRAY                           */
/*                                                                           */
/*  Return the (x,y) pair closest to the given x value.                      */
/*                                                                           */
/*****************************************************************************/
void get_y_at_x_farray(xdata,ydata,n,x,rx,ry)
     float *xdata,*ydata;
     int n;
     float x,*rx,*ry;
{
  int i;
  int imin;
  float diff,min;

  if (n <= 0)
    exit_error("GET_Y_AT_X_FARRAY","No data in array");

  imin = 0;
  min = fabs(xdata[0] - x);
  for(i=1;i<n;i++){
    diff = fabs(xdata[i] - x);
    if (diff < min){
      imin = i;
      min = diff;
    }
  }

  *rx = xdata[imin];
  *ry = ydata[imin];
}
/**************************************-**************************************/
/*                                                                           */
/*                           INTERP_MASK_CIRCULAR_W                          */
/*                                                                           */
/*  Get the weight for the x-position of the given mask.                     */
/*                                                                           */
/*****************************************************************************/
float interp_mask_circular_w(x,x0,period,masktype,maskpar1,maskpar2)
     float x;                  // x-coord for weight  [0..period)
     float x0;                 // mask is centered here  [0..period)
     float period;             // Distance around circle
     char *masktype;           // Interpolation mask [mn]
     float maskpar1,maskpar2;  // Parameters for mask
{
  float dx,w;

  dx = x - x0;

  if (dx > period/2.0)
    dx -= period;
  if (dx <= -period/2.0)
    dx += period;

  //DEBUG printf("x-x0 = %f - %f      dx %f   period %f\n",x,x0,dx,period);

  if (strcmp(masktype,"hanning")==0){
    //  maskpar1 - wavelength of cosine wave
    if ((dx > -maskpar1/2.0) && (dx < maskpar1/2.0))
      w = 0.5 * (1.0 + cos(2.0*M_PI * dx/maskpar1));
    else
      w = 0.0;  // Weights are zero outside of one period
  }else
    exit_error("INTERP_MASK_CIRCULAR_W","Unknown mask type");

  return w;
}
/**************************************-**************************************/
/*                                                                           */
/*                         INTERP_MASK_CIRCULAR_FARRAY                       */
/*                                                                           */
/*  Return a set of y-values that are interpolated from the x,y data,        */
/*  using a mask of the type specified.                                      */
/*                                                                           */
/*****************************************************************************/
float *interp_mask_circular_farray(xdata,ydata,n,period,xnew,nnew,
				   masktype,maskpar1,maskpar2)
     float *xdata,*ydata;      // Original (x,y) data [n]
     int n;                    // Number of points in original data
     float period;             // Distance around circle
     float *xnew;              // X-coords for interpolated data [nnew]
     int nnew;                 // Number of points in interpolated array
     char *masktype;           // Interpolation mask [mn]
     float maskpar1,maskpar2;  // Parameters for mask
{
  int i,j;
  float *d,totw,maskw;

  d = (float *)myalloc(nnew*sizeof(float));

  for(i=0;i<nnew;i++){  // For each new point

    totw = 0.0;  // Total weight
    for(j=0;j<n;j++){  // Add up weights of each original point

      // Determine mask weight
      maskw = interp_mask_circular_w(xdata[j],xnew[i],period,masktype,
				     maskpar1,maskpar2);

      // DEBUG if (i==0)
      //printf("%2d (%8.2f) x %f   w %f\n",j,xnew[i],xdata[j],maskw);
      
      d[i] += maskw * ydata[j];
      totw += maskw;
    }
    d[i] /= totw;

    //if (i==0)
    //exit(0);
  }

  return d;
}
/**************************************-**************************************/
/*                                                                           */
/*                              LIN_INTERP_FARRAY                            */
/*                                                                           */
/*  Use linear interpolation to return a value for the specified             */
/*  coordinate.  Return the closest value if "inx" is outside of domain of   */
/*  array.                                                                   */
/*                                                                           */
/*****************************************************************************/
float lin_interp_farray(x,y,n,inx)
     float *x,*y;
     int n;
     float inx;
{
  int i;
  
  if (inx < x[0])
    return y[0];
  else if (inx >= x[n-1])
    return y[n-1];
  else{
    i=0;
    while(i<(n-1)){
      if (x[i] == inx)
	return y[i];
      else if ((x[i] < inx)&&(x[i+1] > inx)){
	return y[i] + (inx-x[i])/(x[i+1]-x[i]) * (y[i+1]-y[i]);
      }else
	i += 1;
    }
  }
  printf("n = %d  inx = %f\n",n,inx);
  printf("x[0]= %f   x[n-1]= %f\n",x[0],x[n-1]);
  exit_error("LIN_INTERP_FARRAY","Error");
  return 0.0; // This makes lint happy
}
/**************************************-**************************************/
/*                                                                           */
/*                            INTERPOLATE_2D_FARRAY                          */
/*                                                                           */
/*  Use linear interpolation to return a value for the specified             */
/*  coordinate.  Return the closest value if "inx" is outside of domain of   */
/*  array.                                                                   */
/*                                                                           */
/*****************************************************************************/
float interpolate_2d_farray(data,m,n,x,y)
     float **data;
     int m,n;
     float x,y;
{
  int x1,x2,y1,y2;
  float dx,dy,val,val1,val2;

  if ((x<0) || (x>(float)(m-1)))
    exit_error("INTERPOLATE_2D_FARRAY","x is out of range");
  if ((y<0) || (y>(float)(n-1)))
    exit_error("INTERPOLATE_2D_FARRAY","y is out of range");

  x1 = (int)x;
  x2 = x1 + 1;
  y1 = (int)y;
  y2 = y1 + 1;

  dx = x - (float)x1;
  dy = y - (float)y1;

  /*printf("%d %d %d %d\n",x1,y1,x2,y2);*/

  if ((x == (float)(m-1)) && (y == (float)(n-1)))
    val = data[x1][y1];
  else if (x == (float)(m-1)){
    val1 = data[x1][y1];
    val2 = data[x1][y2];
    val  = (1-dy)*val1         + dy*val2;
  }else if (y == (float)(n-1)){
    val1 = (1-dx)*data[x1][y1] + dx*data[x2][y1];
    val  = val1;
  }else{
    val1 = (1-dx)*data[x1][y1] + dx*data[x2][y1];
    val2 = (1-dx)*data[x1][y2] + dx*data[x2][y2];
    val  = (1-dy)*val1         + dy*val2;
  }

  /*printf("val = %f\n",val);*/

  return val;
}
/**************************************-**************************************/
/*                                                                           */
/*                                RESAMPLE_FARRAY                            */
/*                                                                           */
/*  Return the array of y values at the specified "sx" x values.             */
/*  Use linear interpolation.                                                */
/*                                                                           */
/*****************************************************************************/
float *resample_farray(x,y,n,sx,sn)
     float *x,*y;
     int n;
     float *sx;
     int sn;
{
  int i,k;
  float *sy;

  sy = (float *)myalloc(sn*sizeof(float));

  if ((sx[0] < x[0])||(sx[sn-1] > x[n-1])){
    printf("%f %f    %f %f  n=%d\n",sx[0],x[0],sx[sn-1],x[n-1],n);
    exit_error("RESAMPLE_FARRAY","Sampling point out of range");
  }

  k = 0;
  for(i=0;i<sn;i++){
    while(x[k+1] < sx[i])
      k += 1;
    sy[i] = y[k] + (sx[i]-x[k])/(x[k+1]-x[k]) * (y[k+1]-y[k]);
  }
  return sy;
}
/**************************************-**************************************/
/*                                                                           */
/*                              OVER_SAMPLE_FARRAY                           */
/*                                                                           */
/*  Return an array 'n*s' long with each value in 'x' repeated 's' times.    */
/*                                                                           */
/*****************************************************************************/
float *over_sample_farray(x,n,s)
     float *x;
     int n,s;
{
  int i,j,k;
  float *y;

  if (s < 1)
    exit_error("OVER_SAMPLE_FARRAY","k < 1");
  if (n < 1)
    exit_error("OVER_SAMPLE_FARRAY","n < 1");

  y = (float *)myalloc(n*s*sizeof(float));
  k = 0;
  for(i=0;i<n;i++)
    for(j=0;j<s;j++){
      y[k] = x[i];
      k++;
    }
  return y;
}
/**************************************-**************************************/
/*                                                                           */
/*                          OVER_SAMPLE_FARRAY_INTERP                        */
/*                                                                           */
/*  Return an array (n-1)*s+1 long, where each value, except the last, is    */
/*  replaced with 'n' values using linear interpolation.                     */
/*                                                                           */
/*****************************************************************************/
float *over_sample_farray_interp(x,n,s,rn)
     float *x;
     int n,s;
     int *rn;    // (n-1)*s + 1
{
  int i,j,k;
  int nn;
  float *y,dy;

  if (s < 1)
    exit_error("OVER_SAMPLE_FARRAY","k < 1");
  if (n < 1)
    exit_error("OVER_SAMPLE_FARRAY","n < 1");

  nn = (n-1)*s + 1;
  
  y = (float *)myalloc(nn*sizeof(float));
  k = 0;
  for(i=0;i<(n-1);i++){
    dy = (x[i+1] - x[i])/(float)s;
    for(j=0;j<s;j++){
      y[k] = x[i] + (float)j*dy;
      k++;
    }
  }
  y[nn-1] = x[n-1];

  *rn = nn;

  return y;
}
/**************************************-**************************************/
/*                                                                           */
/*                              SUBSAMPLE_FARRAY                             */
/*                                                                           */
/*****************************************************************************/
void subsample_farray(data,n,rdata,rn,ss,origin)
     float *data;    // original data
     int n;
     float **rdata;  // returned data
     int *rn;
     int ss;         // original points per new sampling unit
     int origin;     // a subsample must be on this point
{
  int i,pos;
  int start;

  check_int_range(0,n-1,origin,"subsample_farray: origin");
  check_int_range(0,n-1,ss,"subsample_farray: ss");

  if (n > 0){
    start = origin - (origin/ss)*ss;
    *rn = origin/ss + 1 + (n-origin-1)/ss;
  }else{
    start = 0;
    *rn = 0;
  }
  *rdata = (float *)myalloc(*rn*sizeof(float));
  pos = 0;
  for (i=start;i<n;i+=ss){
    (*rdata)[pos] = data[i];
    pos += 1;
  }
  if (pos != *rn)
    exit_error("SUBSAMPLE_FARRAY","pos != rn");
}
/**************************************-**************************************/
/*                                                                           */
/*                            SUBSAMPLE_2D_FARRAY                            */
/*                                                                           */
/*****************************************************************************/
void subsample_2d_farray(data,m,n,rdata,rm,rn,ss,originx,originy)
     float **data;         // original data
     int m,n;
     float ***rdata;       // returned data
     int *rm,*rn;
     int ss;               // original points per new sampling unit
     int originx,originy;  // a subsample must be on this point
{
  int i,pos;
  int startx;

  startx = originx - (originx/ss)*ss;
  *rm = originx/ss + 1 + (m-originx-1)/ss;

  *rdata = (float **)myalloc(*rm*sizeof(float *));
  pos = 0;
  for (i=startx;i<m;i+=ss){
    subsample_farray(data[i],n,&(*rdata)[pos],rn,ss,originy);
    pos += 1;
  }
  if (pos != *rm)
    exit_error("SUBSAMPLE_2D_FARRAY","pos != rm");
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_BINARY_THRESHOLD_FARRAY                       */
/*                                                                           */
/*  Set all values greater than or equal to the critical value to 1.0, and   */
/*  set all other values to 0.0.  Return the new array.                      */
/*                                                                           */
/*****************************************************************************/
float *get_binary_threshold_farray(data,n,crit)
     float *data;
     int n;
     float crit;
{
  int i;
  float *t;

  t = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    if (data[i] >= crit)
      t[i] = 1.0;
    else
      t[i] = 0.0;

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                           BINARY_THRESHOLD_FARRAY                         */
/*                                                                           */
/*  Set all values greater than or equal to the critical value to 1.0, and   */
/*  set all other values to 0.0.  Return the new array.                      */
/*                                                                           */
/*****************************************************************************/
void binary_threshold_farray(data,n,crit)
     float *data;
     int n;
     float crit;
{
  int i;

  for(i=0;i<n;i++)
    if (data[i] >= crit)
      data[i] = 1.0;
    else
      data[i] = 0.0;
}
/**************************************-**************************************/
/*                                                                           */
/*                    GET_THRESHOLD_CROSSING_INDICES_FARRAY                  */
/*                                                                           */
/*  Return the indices in the farray where the threshold 'crit' is first     */
/*  crossed.                                                                 */
/*                                                                           */
/*****************************************************************************/
void get_threshold_crossing_indices_farray(data,n,crit,rndx,rn)
     float *data;
     int n;
     float crit;
     int **rndx,*rn;
{
  int i;
  int *t,count,flag;
  float *b;

  /*** Get binarized data ***/
  b = get_binary_threshold_farray(data,n,crit);

  count = 0;
  flag = 0; /* 0-past value is below threshold, 1-past value above thresh */
  for(i=0;i<n;i++){
    if ((flag==0) && (b[i] > 0.0))
      count += 1;
    if (b[i] > 0.0)
      flag = 1;
    else
      flag = 0;
  }

  t = (int *)myalloc(count*sizeof(int));
  count = 0;
  flag = 0;
  for(i=0;i<n;i++){
    if ((flag==0) && (b[i] > 0.0)){
      t[count] = i;
      count += 1;
    }
    if (b[i] > 0.0)
      flag = 1;
    else
      flag = 0;
  }

  myfree(b);

  *rndx = t; *rn = count;
}
/**************************************-**************************************/
/*                                                                           */
/*                              DISCRETIZE_FARRAY                            */
/*                                                                           */
/*  Use the "ncrit" values in "crit" to discretize the farray, assigning     */
/*  "dval[i]" to the ith range (0 < i <= crit).  The array "crit" is         */
/*  assumed to be in ascending order.                                        */
/*                                                                           */
/*  For example, for x < crit[0], x is mapped to dval[0].                    */
/*  For crit[0] <= x < crit[1], x is mapped to dval[1].                      */
/*     . . .                                                                 */
/*  For x >= crit[ncrit-1], x is mapped to dval[ncrit].                      */
/*                                                                           */
/*****************************************************************************/
void discretize_farray(data,n,crit,ncrit,dval)
     float *data;
     int n;
     float *crit;
     int ncrit;
     float *dval; /* Has ncrit+1 entries. */
{
  int i;
  int done,k;
  
  for(i=0;i<n;i++){
    k = 0;
    done = 0;
    while(!done){
      if (k == ncrit)
	done = 1;
      else if (data[i] < crit[k])
	done = 1;
      else
	k += 1;
    }
    data[i] = dval[k];
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              ADD_CONST_FARRAY                             */
/*                                                                           */
/*****************************************************************************/
void add_const_farray(data,n,c)
     float *data;
     int n;
     float c;
{
  int i;

  for(i=0;i<n;i++)
    data[i] += c;
}
/**************************************-**************************************/
/*                                                                           */
/*                             ADD_CONST_2D_FARRAY                           */
/*                                                                           */
/*****************************************************************************/
void add_const_2d_farray(data,x0,xn,y0,yn,fconst)
     float **data;
     int x0,xn,y0,yn;
     float fconst;
{
  int i,j;

  if (fconst == 0.0)
    return;

  for(i=x0;i<(x0+xn);i++)
    for(j=y0;j<(y0+yn);j++)
      data[i][j] += fconst;
}
/**************************************-**************************************/
/*                                                                           */
/*                              ADD_SCALE_FARRAYS                            */
/*                                                                           */
/*****************************************************************************/
float *add_scale_farrays(data1,data2,f1,f2,n)
     float *data1,*data2;
     float f1,f2;
     int n;
{
  int i;
  float *data;

  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = f1*data1[i] + f2*data2[i];
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             ADD_SQUARED_FARRAYS                           */
/*                                                                           */
/*****************************************************************************/
float *add_squared_farrays(d1,d2,n)
     float *d1,*d2;
     int n;
{
  int i;
  float *d;

  d = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    d[i] = d1[i]*d1[i] + d2[i]*d2[i];

  return d;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 ADD_FARRAYS                               */
/*                                                                           */
/*****************************************************************************/
float *add_farrays(data1,data2,n)
     float *data1,*data2;
     int n;
{
  int i;
  float *data;

  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = data1[i] + data2[i];
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 ADD_DARRAYS                               */
/*                                                                           */
/*****************************************************************************/
double *add_darrays(data1,data2,n)
     double *data1,*data2;
     int n;
{
  int i;
  double *data;

  data = (double *)myalloc(n*sizeof(double));
  for(i=0;i<n;i++)
    data[i] = data1[i] + data2[i];
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                                ADD_TO_FARRAY                              */
/*                                                                           */
/*  Add 'data2' into 'data1'.  'data2' is unchanged.                         */
/*                                                                           */
/*****************************************************************************/
void add_to_farray(data1,data2,n)
     float *data1,*data2;
     int n;
{
  int i;

  for(i=0;i<n;i++)
    data1[i] += data2[i];
}
/**************************************-**************************************/
/*                                                                           */
/*                            ACCUMULATE_2D_FARRAY                           */
/*                                                                           */
/*****************************************************************************/
void accumulate_2d_farray(total,data,n,m)
     float **total,**data;
     int n,m;
{
  int i,j;

  // ******* WYETH - this def is not consistent w/ other accumulate...
  // ******* WYETH - this def is not consistent w/ other accumulate...
  // ******* WYETH - this def is not consistent w/ other accumulate...

  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      total[i][j] += data[i][j];
}
/**************************************-**************************************/
/*                                                                           */
/*                              SUBTRACT_FARRAYS                             */
/*                                                                           */
/*****************************************************************************/
float *subtract_farrays(data1,data2,n)
     float *data1,*data2;
     int n;
{
  int i;
  float *data;

  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = data1[i] - data2[i];
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                              SUBTRACT_DARRAYS                             */
/*                                                                           */
/*****************************************************************************/
double *subtract_darrays(data1,data2,n)
     double *data1,*data2;
     int n;
{
  int i;
  double *data;

  data = (double *)myalloc(n*sizeof(double));
  for(i=0;i<n;i++)
    data[i] = data1[i] - data2[i];
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                         SUBTRACT_INTEGRATE_FARRAYS                        */
/*                                                                           */
/*****************************************************************************/
float subtract_integrate_farrays(data1,data2,n)
     float *data1,*data2;
     int n;
{
  int i;
  float x;

  x = 0.0;
  for(i=0;i<n;i++)
    x += data1[i] - data2[i];
  return x;
}
/**************************************-**************************************/
/*                                                                           */
/*                               ADD_2D_FARRAYS                              */
/*                                                                           */
/*****************************************************************************/
float **add_2d_farrays(data1,data2,m,n)
     float **data1,**data2;
     int m,n;
{
  int i,j;
  float **data;

  data = (float **)myalloc(m*sizeof(float *));
  for(i=0;i<m;i++){
    data[i] = (float *)myalloc(n*sizeof(float));
    for(j=0;j<n;j++)
      data[i][j] = data1[i][j] + data2[i][j];
  }

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             SUBTRACT_2D_FARRAYS                           */
/*                                                                           */
/*****************************************************************************/
float **subtract_2d_farrays(data1,data2,m,n)
     float **data1,**data2;
     int m,n;
{
  int i;
  float **data;

  data = (float **)myalloc(m*sizeof(float *));
  for(i=0;i<m;i++)
    data[i] = subtract_farrays(data1[i],data2[i],n);
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            SUBTRACT_3D_FARRAYS                            */
/*                                                                           */
/*****************************************************************************/
float ***subtract_3d_farrays(data1,data2,x0,xn,y0,yn,z0,zn)
     float ***data1,***data2;
     int x0,xn,y0,yn,z0,zn;
{
  int i,j,k;
  float ***data;

  data = get_3d_farray(xn,yn,zn);
  for(i=x0;i<(x0+xn);i++)
    for(j=y0;j<(y0+yn);j++)
      for(k=z0;k<(z0+zn);k++)
	data[i-x0][j-y0][k-z0] = data1[i][j][k] - data2[i][j][k];

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MEAN_SQUARE_ERROR_FARRAYS                         */
/*                                                                           */
/*****************************************************************************/
float mean_square_error_farrays(data1,data2,n)
     float *data1,*data2;
     int n;
{
  int i;
  float e;

  e = 0.0;
  for(i=0;i<n;i++)
    e += (data1[i]-data2[i])*(data1[i]-data2[i]);
  e /= (float)n;

  return e;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 ADD_3D_FARRAYS                            */
/*                                                                           */
/*****************************************************************************/
float ***add_3d_farrays(data1,data2,x0,xn,y0,yn,z0,zn)
     float ***data1,***data2;
     int x0,xn,y0,yn,z0,zn;
{
  int i,j,k;
  float ***data;

  data = get_3d_farray(xn,yn,zn);
  for(i=x0;i<(x0+xn);i++)
    for(j=y0;j<(y0+yn);j++)
      for(k=z0;k<(z0+zn);k++)
	data[i-x0][j-y0][k-z0] = data1[i][j][k] + data2[i][j][k];

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 ADD_3D_TENSORS                            */
/*                                                                           */
/*****************************************************************************/
float ***add_3d_tensors(data1,data2,xn,yn,zn)
     float ***data1,***data2;
     int xn,yn,zn;
{
  int i,j,k;
  float ***data;

  data = f3tensor(1,xn,1,yn,1,zn);
  for(i=1;i<=xn;i++)
    for(j=1;j<=yn;j++)
      for(k=0;k<=zn;k++)
	data[i][j][k] = data1[i][j][k] + data2[i][j][k];
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                               AND_3D_TENSORS                              */
/*                                                                           */
/*  If both values are greater than zero, use 't' else use 'f'.              */
/*                                                                           */
/*****************************************************************************/
float ***and_3d_tensors(data1,data2,xn,yn,zn,t,f)
     float ***data1,***data2;
     int xn,yn,zn;
     float t,f;
{
  int i,j,k;
  float ***data;

  data = f3tensor(1,xn,1,yn,1,zn);
  for(i=1;i<=xn;i++)
    for(j=1;j<=yn;j++)
      for(k=0;k<=zn;k++)
	if ((data1[i][j][k] > 0.0)&&(data2[i][j][k] > 0.0))
	  data[i][j][k] = t;
	else
	  data[i][j][k] = f;
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                               AND_2D_FARRAYS                              */
/*                                                                           */
/*  If both values are greater than zero, use 't' else use 'f'.              */
/*                                                                           */
/*****************************************************************************/
float **and_2d_farrays(data1,data2,xn,yn,t,f)
     float **data1,**data2;
     int xn,yn;
     float t,f;  // Values for True and False
{
  int i,j;
  float **data;

  data = get_2d_farray(xn,yn);
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      if ((data1[i][j] > 0.0) && (data2[i][j] > 0.0))
	data[i][j] = t;
      else
	data[i][j] = f;
    }
  }

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            ADD_CONST_3D_FARRAY                            */
/*                                                                           */
/*****************************************************************************/
void add_const_3d_farray(data,x0,xn,y0,yn,z0,zn,fconst)
     float ***data;
     int x0,xn,y0,yn,z0,zn;
     float fconst;
{
  int i,j,k;

  for(i=x0;i<(x0+xn);i++)
    for(j=y0;j<(y0+yn);j++)
      for(k=z0;k<(z0+zn);k++)
	data[i][j][k] += fconst;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MULTIPLY_FARRAY                              */
/*                                                                           */
/*****************************************************************************/
void multiply_farray(data,n,x)
     float *data;
     int n;
     float x;
{
  int i;

  for(i=0;i<n;i++)
    data[i] *= x;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MULTIPLY_DARRAY                              */
/*                                                                           */
/*****************************************************************************/
void multiply_darray(data,n,x)
     double *data;
     int n;
     double x;
{
  int i;

  for(i=0;i<n;i++)
    data[i] *= x;
}
/**************************************-**************************************/
/*                                                                           */
/*                              GET_SCALE_FARRAY                             */
/*                                                                           */
/*****************************************************************************/
float *get_scale_farray(data,n,x)
     float *data;
     int n;
     float x;
{
  int i;
  float *cdata;

  if (n < 0)
    exit_error("GET_SCALE_FARRAY","n < 0");

  cdata = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    cdata[i] = x * data[i];

  return cdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MULTIPLY_2D_FARRAY                            */
/*                                                                           */
/*****************************************************************************/
void multiply_2d_farray(data,xn,yn,x)
     float **data;
     int xn,yn;
     float x;
{
  int i,j;

  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++)
      data[i][j] *= x;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MULTIPLY_3D_FARRAY                             */
/*                                                                           */
/*****************************************************************************/
void multiply_3d_farray(data,x0,xn,y0,yn,z0,zn,x)
     float ***data;
     int x0,xn,y0,yn,z0,zn;
     float x;
{
  int i,j,k;

  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++)
      for(k=0;k<zn;k++)
	data[x0+i][y0+j][z0+k] *= x;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MULTIPLY_FARRAYS                             */
/*                                                                           */
/*****************************************************************************/
float *multiply_farrays(data1,data2,n)
     float *data1,*data2;
     int n;
{
  int i;
  float *prod;

  prod = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    prod[i] = data1[i]*data2[i];
  return prod;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MULTIPLY_DARRAYS                             */
/*                                                                           */
/*****************************************************************************/
double *multiply_darrays(data1,data2,n)
     double *data1,*data2;
     int n;
{
  int i;
  double *prod;

  prod = (double *)myalloc(n*sizeof(double));
  for(i=0;i<n;i++)
    prod[i] = data1[i]*data2[i];
  return prod;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MULTIPLY_FARRAYS_IN_PLACE                       */
/*                                                                           */
/*  The first farray is changed.                                             */
/*                                                                           */
/*****************************************************************************/
void multiply_farrays_in_place(data1,data2,n)
     float *data1,*data2;
     int n;
{
  int i;

  for(i=0;i<n;i++)
    data1[i] *= data2[i];
}
/**************************************-**************************************/
/*                                                                           */
/*                          MULTIPLY_INTEGRATE_FARRAYS                       */
/*                                                                           */
/*****************************************************************************/
float multiply_integrate_farrays(data1,data2,n)
     float *data1,*data2;
     int n;
{
  int i;
  float x;

  x = 0.0;
  for(i=0;i<n;i++)
    x += data1[i]*data2[i];
  return x;
}
/**************************************-**************************************/
/*                                                                           */
/*                             DOT_PROD_2D_FARRAYS                           */
/*                                                                           */
/*****************************************************************************/
double dot_prod_2d_farrays(data1,data2,x0,y0,xn,yn)
     float **data1,**data2;
     int x0,y0;
     int xn,yn;
{
  int i,j;
  int x1,y1;
  double x;

  x1 = x0 + xn - 1;
  y1 = y0 + yn - 1;

  x = 0.0;
  for(i=x0;i<=x1;i++)
    for(j=y0;j<=y1;j++)
      x += data1[i][j] * data2[i][j];

  return x;
}
/**************************************-**************************************/
/*                                                                           */
/*                               DIVIDE_FARRAYS                              */
/*                                                                           */
/*****************************************************************************/
float *divide_farrays(data1,data2,n)
     float *data1,*data2;
     int n;
{
  int i;
  float *div;

  div = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    div[i] = data1[i]/data2[i];
  return div;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MULTIPLY_2D_FARRAYS                           */
/*                                                                           */
/*****************************************************************************/
float **multiply_2d_farrays(data1,data2,xn,yn)
     float **data1,**data2;
     int xn,yn;
{
  int i,j;
  float **prod;

  prod = get_2d_farray(xn,yn);
  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++)
      prod[i][j] = data1[i][j] * data2[i][j];

  return prod;
}
/**************************************-**************************************/
/*                                                                           */
/*                          HALF_WAVE_RECTIFY_FARRAY                         */
/*                                                                           */
/*****************************************************************************/
void half_wave_rectify_farray(data,n)
     float *data;
     int n;
{
  int i;

  for(i=0;i<n;i++)
    if (data[i] < 0.0)
      data[i] = 0.0;
}
/**************************************-**************************************/
/*                                                                           */
/*                            ABSOLUTE_VALUE_FARRAY                          */
/*                                                                           */
/*****************************************************************************/
void absolute_value_farray(data,n)
     float *data;
     int n;
{
  int i;

  for(i=0;i<n;i++)
    if (data[i] < 0.0)
      data[i] *= -1.0;
}
/**************************************-**************************************/
/*                                                                           */
/*                               SQUARE_FARRAY                               */
/*                                                                           */
/*****************************************************************************/
void square_farray(data,n)
     float *data;
     int n;
{
  int i;

  for(i=0;i<n;i++)
    data[i] *= data[i];
}
/**************************************-**************************************/
/*                                                                           */
/*                               SQUARE_DARRAY                               */
/*                                                                           */
/*****************************************************************************/
void square_darray(data,n)
     double *data;
     int n;
{
  int i;

  for(i=0;i<n;i++)
    data[i] *= data[i];
}
/**************************************-**************************************/
/*                                                                           */
/*                                 POW_FARRAY                                */
/*                                                                           */
/*****************************************************************************/
void pow_farray(data,n,p)
     float *data;
     int n;
     float p;
{
  int i;

  for(i=0;i<n;i++)
    data[i] = pow(data[i],p);
}
/**************************************-**************************************/
/*                                                                           */
/*                                POW_2D_FARRAY                              */
/*                                                                           */
/*****************************************************************************/
void pow_2d_farray(data,xn,yn,p)
     float **data;
     int xn,yn;
     float p;
{
  int i,j;
  float d;

  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++){
      d = data[i][j];
      data[i][j] = pow(d,p);
    }
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_EXP_FARRAY                              */
/*                                                                           */
/*  Base should be 0 (natural), 2 or 10.                                     */
/*                                                                           */
/*****************************************************************************/
float *get_exp_farray(data,n,base)
     float *data;
     int n,base;
{
  int i;
  float *t;

  t = (float *)myalloc(n*sizeof(float));
  if (base==0)
    for(i=0;i<n;i++)
      t[i] = exp(data[i]);
  else if (base==2)
    for(i=0;i<n;i++)
      t[i] = pow(2.0,data[i]); /* exp2(data[i]); */
  else if (base==10)
    for(i=0;i<n;i++)
      t[i] = pow(10.0,data[i]); /* exp10(data[i]); */
  else
    exit_error("GET_EXP_FARRAY","Unknown base value");

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_LOG_FARRAY                              */
/*                                                                           */
/*  Base should be 0 (natural), 2 or 10.                                     */
/*                                                                           */
/*****************************************************************************/
float *get_log_farray(data,n,base)
     float *data;
     int n,base;
{
  int i;
  float *t;

  t = (float *)myalloc(n*sizeof(float));
  if (base==0)
    for(i=0;i<n;i++)
      t[i] = log(data[i]);
  else if (base==2)
    for(i=0;i<n;i++)
      t[i] = my_log2(data[i]);
  else if (base==10)
    for(i=0;i<n;i++)
      t[i] = log10(data[i]);
  else
    exit_error("GET_LOG_FARRAY","Unknown base value");

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                                SUM_FARRAY                                 */
/*                                                                           */
/*****************************************************************************/
float sum_farray(data,n,start,duration)
     float *data;
     int n,start,duration;
{
  int i;
  int stop;
  float sum;

  stop = start+duration-1;
  if ((start < 0)||(stop >= n)||(duration < 0))
    exit_error("SUM_FARRAY","Parameter error");
  
  sum = 0.0;
  for(i=start;i<=stop;i++)
    sum += data[i];
  return sum;
}
/**************************************-**************************************/
/*                                                                           */
/*                              SUM_SQUARE_FARRAY                            */
/*                                                                           */
/*****************************************************************************/
double sum_square_farray(data,n)
     float *data;
     int n;
{
  int i;
  double sum;

  sum = 0.0;
  for(i=0;i<n;i++)
    sum += data[i]*data[i];
  return sum;
}
/**************************************-**************************************/
/*                                                                           */
/*                            SUM_SQUARE_2D_FARRAY                           */
/*                                                                           */
/*****************************************************************************/
double sum_square_2d_farray(data,n,m)
     float **data;
     int n,m;
{
  int i,j;
  double x,sum;
  
  sum = 0.0;
  for(i=0;i<n;i++)
    for(j=0;j<m;j++){
      x = data[i][j];
      sum += x*x;
    }
  return sum;
}
/**************************************-**************************************/
/*                                                                           */
/*                               SUM_2D_FARRAY                               */
/*                                                                           */
/*****************************************************************************/
double sum_2d_farray(data,x0,xn,y0,yn)
     float **data;
     int x0,xn,y0,yn;
{
  int i1,i2;
  double sum;

  sum = 0.0;
  for(i1=x0;i1<(x0+xn);i1++)
    for(i2=y0;i2<(y0+yn);i2++)
      sum += data[i1][i2];
  return sum;
}
/**************************************-**************************************/
/*                                                                           */
/*                               SUM_3D_FARRAY                               */
/*                                                                           */
/*****************************************************************************/
double sum_3d_farray(data,x0,xn,y0,yn,z0,zn)
     float ***data;
     int x0,xn,y0,yn,z0,zn;
{
  int i1,i2,i3;
  double sum;

  sum = 0.0;
  for(i1=x0;i1<(x0+xn);i1++)
    for(i2=y0;i2<(y0+yn);i2++)
      for(i3=z0;i3<(z0+zn);i3++)
	sum += data[i1][i2][i3];
  return sum;
}
/**************************************-**************************************/
/*                                                                           */
/*                             SUM_ABS_3D_FARRAY                             */
/*                                                                           */
/*****************************************************************************/
double sum_abs_3d_farray(data,x0,xn,y0,yn,z0,zn)
     float ***data;
     int x0,xn,y0,yn,z0,zn;
{
  int i1,i2,i3;
  double sum;

  sum = 0.0;
  for(i1=x0;i1<(x0+xn);i1++)
    for(i2=y0;i2<(y0+yn);i2++)
      for(i3=z0;i3<(z0+zn);i3++){
	if (data[i1][i2][i3] > 0.0)
	  sum += data[i1][i2][i3];
	else
	  sum -= data[i1][i2][i3];
      }

  return sum;
}
/**************************************-**************************************/
/*                                                                           */
/*                            SUM_SQUARE_3D_FARRAY                           */
/*                                                                           */
/*****************************************************************************/
double sum_square_3d_farray(data,x0,xn,y0,yn,z0,zn)
     float ***data;
     int x0,xn,y0,yn,z0,zn;
{
  int i1,i2,i3;
  double sum;

  sum = 0.0;
  for(i1=x0;i1<(x0+xn);i1++)
    for(i2=y0;i2<(y0+yn);i2++)
      for(i3=z0;i3<(z0+zn);i3++)
	sum += data[i1][i2][i3]*data[i1][i2][i3];
  return sum;
}
/**************************************-**************************************/
/*                                                                           */
/*                               SUM_ABS_FARRAY                              */
/*                                                                           */
/*****************************************************************************/
double sum_abs_farray(data,n)
     float *data;
     int n;
{
  int i;
  double sum;

  sum = 0.0;
  for(i=0;i<n;i++)
    if (data[i] > 0.0)
      sum += data[i];
    else
      sum -= data[i];
  return sum;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_SEQUENCE_FARRAY                           */
/*                                                                           */
/*****************************************************************************/
float *get_sequence_farray(x0,step,n)
     double x0,step;
     int n;
{
  int i;
  float *data;

  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = x0 + (double)i * step;
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MAX_OF_FARRAY                               */
/*                                                                           */
/*****************************************************************************/
float max_of_farray(data,n)
     float *data;
     int n;
{
  int i;
  float max;
  
  if (n < 0) {
    printf("(max_of_farray):  Array length < 0.\n");
    return(0.0);
  }
  max = data[0];
  for (i=1;i<n;i++)
    if (data[i] > max)
      max = data[i];
  return max;
}
/**************************************-**************************************/
/*                                                                           */
/*                                MIN_OF_FARRAY                              */
/*                                                                           */
/*****************************************************************************/
float min_of_farray(data,n)
     float *data;
     int n;
{
  int i;
  float min;
  
  if (n < 0) {
    printf("(min_of_farray):  Array length < 0.\n");
    return(0.0);
  }
  min = data[0];
  for (i=1;i<n;i++){
    if (data[i] < min)
      min = data[i];
  }
  return min;
}
/**************************************-**************************************/
/*                                                                           */
/*                              GET_MIN_MAX_FARRAY                           */
/*                                                                           */
/*****************************************************************************/
void get_min_max_farray(data,n,rmin,rmax)
     float *data;
     int n;
     float *rmin,*rmax;
{
  *rmin = min_of_farray(data,n);
  *rmax = max_of_farray(data,n);
}
/**************************************-**************************************/
/*                                                                           */
/*                              MAX_INDEX_FARRAY                             */
/*                                                                           */
/*****************************************************************************/
int max_index_farray(data,n)
     float *data;
     int n;
{
  int i;
  int imax;
  float max;
  
  if (n < 0)
    exit_error("MAX_INDEX_FARRAY","Array length < 0");

  max = data[0];
  imax = 0;
  for (i=1;i<n;i++)
    if (data[i] > max){
      max = data[i];
      imax = i;
    }
  return imax;
}
/**************************************-**************************************/
/*                                                                           */
/*                        GET_INDEX_MIN_LEVEL_RIGHT_FARRAY                   */
/*                                                                           */
/*  Return the index of the first value >= "t" at or to the right of "k".    */
/*                                                                           */
/*****************************************************************************/
int get_index_min_level_right_farray(data,t,n,k)
     float *data,t;
     int n,k;
{
  int i,j;

  j = -1;
  i = k;
  while((i<n)&&(j<0)){
    if (data[i] >= t)
      j = i;
    else
      i += 1;
  }
  return j;
}
/**************************************-**************************************/
/*                                                                           */
/*                        GET_INDEX_MAX_LEVEL_RIGHT_FARRAY                   */
/*                                                                           */
/*  Return the index of the first value <= "t" at or to the right of "k".    */
/*                                                                           */
/*****************************************************************************/
int get_index_max_level_right_farray(data,t,n,k)
     float *data,t;
     int n,k;
{
  int i,j;

  j = -1;
  i = k;
  while((i<n)&&(j<0)){
    if (data[i] <= t)
      j = i;
    else
      i += 1;
  }
  return j;
}
/**************************************-**************************************/
/*                                                                           */
/*                       GET_INDEX_LOCAL_MIN_RIGHT_FARRAY                    */
/*                                                                           */
/*  Return the index of the local minimum to the right of k.                 */
/*                                                                           */
/*****************************************************************************/
int get_index_local_min_right_farray(data,n,k)
     float *data;
     int n,k;
{
  int i,j;
  int done;

  j = k; /* Index of current minimum */
  i = k+1; /* Current position */
  done = 0;
  while(!done){
    if (i >= n)
      done = 1;
    else{
      if (data[i] > data[i-1])
	done = 1;
      else if (data[i] < data[i-1])
	j = i;
      i += 1;
    }
  }
  return j;
}
/**************************************-**************************************/
/*                                                                           */
/*                       GET_INDEX_LOCAL_MIN_LEFT_FARRAY                     */
/*                                                                           */
/*  Return the index of the local minimum to the left of k.                  */
/*                                                                           */
/*****************************************************************************/
int get_index_local_min_left_farray(data,n,k)
     float *data;
     int n,k;
{
  int i,j;
  int done;

  j = k; /* Index of current minimum */
  i = k-1; /* Current position */
  done = 0;
  while(!done){
    if (i < 0)
      done = 1;
    else{
      if (data[i] > data[i+1])
	done = 1;
      else if (data[i] < data[i+1])
	j = i;
      i -= 1;
    }
  }
  return j;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_MAGNITUDE_RANGE_FARRAY                      */
/*                                                                           */
/*  Return (max-min) value in the segment of the farray.                     */
/*                                                                           */
/*****************************************************************************/
float get_magnitude_range_farray(data,n,t0,tn)
     float *data;
     int n,t0,tn;
{
  float min,max;

  if (!((t0>=0)&&(t0<n)&&((t0+tn)<=n)))
    exit_error("GET_MAGNITUDE_RANGE_FARRAY","Invalid subarray specification");
  
  min = min_of_farray(data+t0,tn);
  max = max_of_farray(data+t0,tn);
  return (max - min);
}
/**************************************-**************************************/
/*                                                                           */
/*                        GET_INDEX_NEAREST_VALUE_FARRAY                     */
/*                                                                           */
/*  Return the index of the (first) value in the array which is closest to   */
/*  the key value.                                                           */
/*                                                                           */
/*****************************************************************************/
int get_index_nearest_value_farray(data,n,key)
     float *data;
     int n;
     float key;
{
  int i,k;
  float diff;
  
  k = -1;
  if (n >= 1){
    k = 0;
    diff = fabs(data[0]-key);
    for (i=1;i<n;i++)
      if (fabs(data[i]-key) < diff){
	diff = fabs(data[i]-key);
	k = i;
      }
  }
  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MAX_COORD_FARRAY                             */
/*                                                                           */
/*  Return the index of the *first* occurrence of the maximum value in the   */
/*  floating point array.                                                    */
/*                                                                           */
/*  WYETH - This could be renamed "get_index_max_farray".                    */
/*                                                                           */
/*****************************************************************************/
int max_coord_farray(data,n)
     float *data;
     int n;
{
  int i,k;
  float max;
  
  if (n < 0) {
    printf("(max_coord_farray):  Array length < 0.\n");
    return(-1);
  }
  k = 0;
  max = data[0];
  for (i=1;i<n;i++)
    if (data[i] > max){
      max = data[i];
      k = i;
    }
  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MIN_COORD_FARRAY                             */
/*                                                                           */
/*****************************************************************************/
int min_coord_farray(data,n)
     float *data;
     int n;
{
  int i;
  int k;
  float min;
  
  if (n < 0) {
    printf("(min_coord_farray):  Array length < 0.\n");
    return(-1);
  }
  k = 0;
  min = data[0];
  for (i=1;i<n;i++)
    if (data[i] < min){
      min = data[i];
      k = i;
    }
  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MAX_OF_2D_FARRAY                             */
/*                                                                           */
/*****************************************************************************/
float max_of_2d_farray(data,m,n)
     float **data;
     int m,n;
{
  int i;
  float *mdata,max;

  if ((n<0)||(m<0)){
    printf("(max_of_2d_farray):  Array length < 0.\n");
    return(0.0);
  }
  mdata = (float *)myalloc(m*sizeof(float));
  for (i=0;i<m;i++)
    mdata[i] = max_of_farray(data[i],n);
  max = max_of_farray(mdata,m);
  myfree(mdata);
  return max;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MIN_OF_2D_FARRAY                            */
/*                                                                           */
/*****************************************************************************/
float min_of_2d_farray(data,m,n)
     float **data;
     int m,n;
{
  int i;
  float *mdata,min;

  if ((n<0)||(m<0)){
    printf("(min_of_2d_farray):  Array length < 0.\n");
    return(0.0);
  }
  mdata = (float *)myalloc(m*sizeof(float));
  for (i=0;i<m;i++)
    mdata[i] = min_of_farray(data[i],n);
  min = min_of_farray(mdata,m);
  myfree(mdata);
  return min;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_MIN_MAX_2D_FARRAY                          */
/*                                                                           */
/*****************************************************************************/
void get_min_max_2d_farray(data,m,n,rmin,rmax)
     float **data;
     int m,n;
     float *rmin,*rmax;
{
  int i,j;
  float min,max;

  if ((n<=0)||(m<=0))
    exit_error("GET_MIN_MAX_2D_FARRAY","Array length <= 0");

  min = max = data[0][0];
  for (i=0;i<m;i++)
    for (j=0;j<n;j++)
      if (data[i][j] < min)
	min = data[i][j];
      else if (data[i][j] > max)
	max = data[i][j];

  *rmin = min; *rmax = max;
}
/**************************************-**************************************/
/*                                                                           */
/*                          ABSOLUTE_VALUE_2D_FARRAY                         */
/*                                                                           */
/*****************************************************************************/
void absolute_value_2d_farray(data,m,n)
     float **data;
     int m,n;
{
  int i,j;

  if ((n<=0)||(m<=0))
    exit_error("ABSOLUTE_VALUE_2D_FARRAY","Array length <= 0");

  for (i=0;i<m;i++)
    for (j=0;j<n;j++)
      if (data[i][j] < 0.0)
	data[i][j] *= -1.0;
}
/**************************************-**************************************/
/*                                                                           */
/*                         HALF_WAVE_RECTIFY_2D_FARRAY                       */
/*                                                                           */
/*****************************************************************************/
void half_wave_rectify_2d_farray(data,n,m)
     float **data;
     int n,m;
{
  int i,j;

  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      if (data[i][j] < 0.0)
	data[i][j] = 0.0;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MAX_COORD_2D_FARRAY                           */
/*                                                                           */
/*****************************************************************************/
void max_coord_2d_farray(data,m,n,im,in)
     float **data;
     int m,n,*im,*in;
{
  int i,j;
  float max;
  int x,y;

  if ((n<0)||(m<0)){
    printf("(max_coord_2d_farray):  Array length < 0.\n");
  }
  max = data[0][0];
  x = y = 0;
  for (i=0;i<m;i++)
    for (j=0;j<n;j++)
      if (data[i][j] > max){
	x = i;
	y = j;
	max = data[i][j];
      }
  *im = x;
  *in = y;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MIN_COORD_2D_FARRAY                           */
/*                                                                           */
/*****************************************************************************/
void min_coord_2d_farray(data,m,n,im,in)
     float  **data;
     int m,n,*im,*in;
{
  int i,j;
  float min;
  int x,y;

  if ((n<0)||(m<0)){
    printf("(min_coord_2d_farray):  Array length < 0.\n");
  }
  min = data[0][0];
  x = y = 0;
  for (i=0;i<m;i++)
    for (j=0;j<n;j++)
      if (data[i][j] < min){
	x = i;
	y = j;
	min = data[i][j];
      }
  *im = x;
  *in = y;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MIN_MAX_COORD_2D_FARRAY                         */
/*                                                                           */
/*****************************************************************************/
void min_max_coord_2d_farray(data,m,n,rminx,rminy,rmaxx,rmaxy)
     float  **data;
     int m,n;
     int *rminx,*rminy,*rmaxx,*rmaxy;
{
  int i,j;
  float min,max;

  if ((n<0)||(m<0))
    exit_error("MIN_MAX_COORD_2D_FARRAY","Array length < 0");

  min = max = data[0][0];
  *rminx = *rminy = *rmaxx = *rmaxy = 0;
  for (i=0;i<m;i++){
    for (j=0;j<n;j++)
      if (data[i][j] < min){
	*rminx = i;
	*rminy = j;
	min = data[i][j];
      }else if (data[i][j] > max){
	*rmaxx = i;
	*rmaxy = j;
	max = data[i][j];
      }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_MIN_MAX_3D_FARRAY                          */
/*                                                                           */
/*****************************************************************************/
void get_min_max_3d_farray(data,x0,xn,y0,yn,z0,zn,rmin,rmax)
     float ***data;
     int x0,xn,y0,yn,z0,zn;
     float *rmin,*rmax;
{
  int i,j,k;
  float min,max;

  min = max = data[x0][y0][z0];
  for (i=x0;i<(x0+xn);i++)
    for (j=y0;j<(y0+yn);j++)
      for (k=z0;k<(z0+zn);k++)
	if (data[i][j][k] < min)
	  min = data[i][j][k];
	else if (data[i][j][k] > max)
	  max = data[i][j][k];
  *rmin = min;
  *rmax = max;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_MAX_COORD_3D_FARRAY                         */
/*                                                                           */
/*****************************************************************************/
void get_max_coord_3d_farray(data,x0,xn,y0,yn,z0,zn,rx,ry,rz)
     float ***data;
     int x0,xn,y0,yn,z0,zn;
     int *rx,*ry,*rz;
{
  int i,j,k;
  int x,y,z;
  float max;

  max = data[x0][y0][z0];
  x = x0;
  y = y0;
  z = z0;
  for (i=x0;i<(x0+xn);i++)
    for (j=y0;j<(y0+yn);j++)
      for (k=z0;k<(z0+zn);k++)
	if (data[i][j][k] > max){
	  max = data[i][j][k];
	  x = i;
	  y = j;
	  z = k;
	}

  *rx = x; *ry = y; *rz = z;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_MIN_MAX_COORD_3D_FARRAY                       */
/*                                                                           */
/*****************************************************************************/
void get_min_max_coord_3d_farray(data,x0,xn,y0,yn,z0,zn,ra,rb,rc,rx,ry,rz)
     float ***data;
     int x0,xn,y0,yn,z0,zn;
     int *ra,*rb,*rc;        // Coord of min
     int *rx,*ry,*rz;        // Coord of max
{
  int i,j,k;
  int a,b,c,x,y,z;
  float max,min;

  min = data[x0][y0][z0];
  a = x0;
  b = y0;
  c = z0;
  max = data[x0][y0][z0];
  x = x0;
  y = y0;
  z = z0;
  for (i=x0;i<(x0+xn);i++)
    for (j=y0;j<(y0+yn);j++)
      for (k=z0;k<(z0+zn);k++)
	if (data[i][j][k] > max){
	  max = data[i][j][k];
	  x = i;
	  y = j;
	  z = k;
	}else if (data[i][j][k] < min){
	  min = data[i][j][k];
	  a = i;
	  b = j;
	  c = k;
	}

  *ra = a; *rb = b; *rc = c;
  *rx = x; *ry = y; *rz = z;
}
/**************************************-**************************************/
/*                                                                           */
/*                          FIND_EXTREMA_COORD_FARRAY                        */
/*                                                                           */
/*  Find local minima, maxima, and horizontal SADDLE points, but do          */
/*  not include endpoints.  Return the array containing the coordin-         */
/*  ates of the extrema.                                                     */
/*                                                                           */
/*****************************************************************************/
int *find_extrema_coord_farray(data,n,nex)
     float *data;
     int n,*nex;
{
  int i,k;
  int *extrema;

  extrema = (int *)myalloc(n*sizeof(int));

  k = 0;
  for(i=1;i<(n-1);i++)
    if (((data[i-1] >= data[i]) && (data[i+1] >= data[i])) ||
	((data[i-1] <= data[i]) && (data[i+1] <= data[i]))){
      extrema[k] = i;
      k += 1;
    }
  *nex = k;
  return extrema;
}
/**************************************-**************************************/
/*                                                                           */
/*                        STRICTLY_INCREASING_FARRAY                         */
/*                                                                           */
/*  Return "1" if the farray is stricly increasing on the interval           */
/*  specified, "0" otherwise.                                                */
/*                                                                           */
/*****************************************************************************/
int strictly_increasing_farray(data,n,start,period)
     float *data;
     int n,start,period;
{
  int i;
  int stop;

  stop = start + period - 1;
  for(i=start;i<stop;i++)
    if (data[i+1] <= data[i])
      return 0;
  return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                        STRICTLY_DECREASING_FARRAY                         */
/*                                                                           */
/*  Return "1" if the farray is stricly decreasing on the interval           */
/*  specified, "0" otherwise.                                                */
/*                                                                           */
/*****************************************************************************/
int strictly_decreasing_farray(data,n,start,period)
     float *data;
     int n,start,period;
{
  int i;
  int stop;

  stop = start + period - 1;
  for(i=start;i<stop;i++)
    if (data[i+1] >= data[i])
      return 0;
  return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                             SPLINE_INTERP_FARRAY                          */
/*                                                                           */
/*  Use spline interpolation to return y-values at the specified x-values.   */
/*                                                                           */
/*****************************************************************************/
float *spline_interp_farray(xdata,ydata,ndata,yp1,ypn,x,n)
     float *xdata,*ydata;
     int ndata;
     float yp1,ypn;    // 1st deriv at endpoints, or >=10^30 for 2nd deriv 0
     float *x;         // x-values for new points
     int n;
{
  int i;
  float *y2,*y;

  y2 = (float *)myalloc(ndata*sizeof(float));
  spline(xdata-1,ydata-1,ndata,yp1,ypn,y2-1); // NumRecInC function

  y = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    splint(xdata-1,ydata-1,y2-1,ndata,x[i],&y[i]);

  myfree(y2);
  return y;
}
/**************************************-**************************************/
/*                                                                           */
/*                             SPLINE_MAX_FARRAY                             */
/*                                                                           */
/*  Return the maximum value of the interpolated array.                      */
/*                                                                           */
/*****************************************************************************/
void spline_max_farray(xdata,ydata,n,nsamp,rx,ry)
     float *xdata,*ydata;   // [n] data
     int n;                 // length of data
     int nsamp;             // Number of sub-samples for interpolation
     float *rx,*ry;         // Return (x,y) values at peak
{
  int i,j,k;
  int ns,imax;
  float *xspline,*yspline,aderiv,bderiv;
  //void append_farray_xy_plot();

  ns = nsamp * (n-1) + 1;  // Number of values in spline interp
  xspline = get_farray(ns);

  // x-coords for spline interp
  k = 0;
  for(i=1;i<n;i++){
    for(j=0;j<nsamp;j++){
      xspline[k] = xdata[i-1]+(xdata[i]-xdata[i-1])*(float)j/(float)nsamp;
      k += 1;
    }
  }
  xspline[k] = xdata[n-1];

  aderiv = bderiv = 2E30; // Use "natural" spline, see nr_util.c
  yspline = spline_interp_farray(xdata,ydata,n,aderiv,bderiv,xspline,ns);

  imax = max_coord_farray(yspline,ns);

  *rx = xspline[imax];
  *ry = yspline[imax];
  //append_farray_xy_plot("zz.spline.pl",xspline,yspline,ns,"spl");

  myfree(xspline);
  myfree(yspline);
}
/**************************************-**************************************/
/*                                                                           */
/*                               PAD_END_FARRAY                              */
/*                                                                           */
/*****************************************************************************/
float *pad_end_farray(data,n,npad,padval)
     float *data;
     int n,npad;
     float padval;
{
  int i;
  float *pdata;

  pdata = (float *)myalloc(npad*sizeof(float));
  for(i=0;i<n;i++)
    pdata[i] = data[i];
  for(i=n;i<npad;i++)
    pdata[i] = padval;
  return pdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                                DIFF_FARRAYS                               */
/*                                                                           */
/*****************************************************************************/
float *diff_farrays(data1,data2,n)
     float *data1,*data2;
     int n;
{
  int i;
  float *diff;

  diff = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    diff[i] = data1[i] - data2[i];
  return diff;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 DIFF_DARRAYS                              */
/*                                                                           */
/*****************************************************************************/
double *diff_darrays(data1,data2,n)
     double *data1,*data2;
     int n;
{
  int i;
  double *diff;

  diff = (double *)myalloc(n*sizeof(double));
  for(i=0;i<n;i++)
    diff[i] = data1[i] - data2[i];
  return diff;
}
/**************************************-**************************************/
/*                                                                           */
/*                                ROC_FARRAYS                                */
/*                                                                           */
/*****************************************************************************/
float roc_farrays(data1,data2,n1,n2)
     float *data1,*data2;
     int n1,n2;
{
  int i,j;
  float area;
  
  area = 0.0;
  for(i=0;i<n1;i++)
    for(j=0;j<n2;j++)
      if (data1[i] > data2[j])
	area += 2.0;
      else if (data1[i] == data2[j])
	area += 1.0;
  area /= (float)(n1*n2*2);
  
  return area;
}
/**************************************-**************************************/
/*                                                                           */
/*                                WRITE_XVIEW                                */
/*                                                                           */
/*****************************************************************************/
void write_xview(data,m,n,outfile)
     float **data;
     int m,n;
     char outfile[];
{
  FILE *fopen(),*fout;
  int i,j;
  int xmax,ymax;
  int datatype; /* 3 for int, 4 for float ? */
  float delta;
  
  printf("  WRITE_XVIEW\n");
  printf("    Writing %s\n",outfile);

  xmax = m-1;
  ymax = n-1;
  delta = 1.0;
  datatype = 4;

  fout = fopen(outfile,"w");

  fwrite((char *)&xmax,sizeof(int),1,fout);
  fwrite((char *)&ymax,sizeof(int),1,fout);
  fwrite((char *)&delta,sizeof(float),1,fout);
  fwrite((char *)&datatype,sizeof(int),1,fout);
  for (i=0;i<n;i++)
    for (j=0;j<m;j++)
      fwrite((char *)&data[j][i],sizeof(float),1,fout);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                           WRITE_XVIEW_TRANSPOSE                           */
/*                                                                           */
/*****************************************************************************/
void write_xview_transpose(data,m,n,outfile)
     float **data;
     int m,n;
     char outfile[];
{
  FILE *fopen(),*fout;
  int i,j;
  int xmax,ymax;
  int datatype; /* 3 for int, 4 for float ? */
  float delta;
  
  printf("  WRITE_XVIEW_TRANSPOSE\n");
  printf("    Writing %s\n",outfile);

  xmax = m-1;
  ymax = n-1;
  delta = 1.0;
  datatype = 4;

  fout = fopen(outfile,"w");

  fwrite((char *)&ymax,sizeof(int),1,fout);
  fwrite((char *)&xmax,sizeof(int),1,fout);
  fwrite((char *)&delta,sizeof(float),1,fout);
  fwrite((char *)&datatype,sizeof(int),1,fout);
  for (j=0;j<m;j++)
    for (i=0;i<n;i++)
      fwrite((char *)&data[j][i],sizeof(float),1,fout);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                                 READ_XVIEW                                */
/*                                                                           */
/*****************************************************************************/
void read_xview(data,m,n,infile)
     float ***data;
     int *m,*n;
     char infile[];
{
  FILE *fopen(),*fin;
  int i,j;
  int datatype; /* 3 for int, 4 for float ? */
  int mm,nn,nr;
  float delta,**d;
  
  printf("  READ_XVIEW\n");
  printf("    Reading %s\n",infile);

  fin = fopen(infile,"r");
  
  nr = fread((char *)&mm,sizeof(int),1,fin);
  nr = fread((char *)&nn,sizeof(int),1,fin);
  nr = fread((char *)&delta,sizeof(float),1,fin);
  nr = fread((char *)&datatype,sizeof(int),1,fin);
  
  mm += 1; /* convert from max subscript to total count */
  nn += 1;
  d = (float **)myalloc(mm*sizeof(float *));
  for (j=0;j<mm;j++)
    d[j] = (float *)myalloc(nn*sizeof(float));
  for (i=0;i<nn;i++)
    for (j=0;j<mm;j++)
      nr = fread((char *)&d[j][i],sizeof(float),1,fin);
  fclose(fin);
  *data = d; *m = mm; *n = nn;
}
/**************************************-**************************************/
/*                                                                           */
/*                             WRITE_XVIEW_MOVIE                             */
/*                                                                           */
/*  Assume that the data has first two dimensions x and y, and then          */
/*  use the third dimension, z, as the movie "frame" dimension.              */
/*                                                                           */
/*****************************************************************************/
void write_xview_movie(data,m,n,frames,outfile,transpose)
     float ***data;
     int m,n,frames;
     char outfile[];
     int transpose; /* if 1, transpose the array before writing */
{
  FILE *fopen(),*fout;
  int i,j,k;
  int xmax,ymax;
  int datatype; /* 3 for int, 4 for float ? */
  float delta;
  
  printf("  WRITE_XVIEW_MOVIE\n");
  printf("    Writing %s\n",outfile);

  xmax = m-1;
  ymax = n-1;
  delta = 1.0;
  datatype = 4;

  fout = fopen(outfile,"w");

  fwrite((char *)&xmax,sizeof(int),1,fout);
  fwrite((char *)&ymax,sizeof(int),1,fout);
  fwrite((char *)&delta,sizeof(float),1,fout);
  fwrite((char *)&datatype,sizeof(int),1,fout);
  if (transpose == 1){
    for (k=0;k<frames;k++)
      for (j=0;j<n;j++)
	for (i=0;i<m;i++)
	  fwrite((char *)&data[i][j][k],sizeof(float),1,fout);
  }else{
    for (k=0;k<frames;k++)
      for (i=0;i<m;i++)
	for (j=0;j<n;j++)
	  fwrite((char *)&data[i][j][k],sizeof(float),1,fout);
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                          WRITE_XVIEW_MOVIE_TENSOR                         */
/*                                                                           */
/*  Assume that the data has first two dimensions x and y, and then          */
/*  use the third dimension, z, as the movie "frame" dimension.              */
/*                                                                           */
/*****************************************************************************/
void write_xview_movie_tensor(data,m,n,frames,outfile,transpose)
     float ***data;
     int m,n,frames;
     char outfile[];
     int transpose; /* if 1, transpose the array before writing */
{
  FILE *fopen(),*fout;
  int i,j,k;
  int xmax,ymax;
  int datatype; /* 3 for int, 4 for float ? */
  float delta;
  
  printf("  WRITE_XVIEW_MOVIE_TENSOR\n");
  printf("    Writing %s\n",outfile);

  xmax = m-1;
  ymax = n-1;
  delta = 1.0;
  datatype = 4;

  fout = fopen(outfile,"w");

  fwrite((char *)&xmax,sizeof(int),1,fout);
  fwrite((char *)&ymax,sizeof(int),1,fout);
  fwrite((char *)&delta,sizeof(float),1,fout);
  fwrite((char *)&datatype,sizeof(int),1,fout);
  if (transpose == 1){
    for (k=1;k<=frames;k++)
      for (j=1;j<=n;j++)
	for (i=1;i<=m;i++)
	  fwrite((char *)&data[i][j][k],sizeof(float),1,fout);
  }else{
    for (k=1;k<=frames;k++)
      for (i=1;i<=m;i++)
	for (j=1;j<=n;j++)
	  fwrite((char *)&data[i][j][k],sizeof(float),1,fout);
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                              WRITE_HIPS_MOVIE                             */
/*                                                                           */
/*  Assume that the data has first two dimensions x and y, and then          */
/*  use the third dimension, z, as the movie "frame" dimension.              */
/*                                                                           */
/*****************************************************************************/
void write_hips_movie(data,x0,xn,y0,yn,z0,zn,outfile,transpose,yrevflag,
		      scaleflag)
     float ***data;
     int x0,xn,y0,yn,z0,zn;
     char outfile[];
     int transpose; /* if 1, transpose the array before writing */
     int yrevflag; /* "xvanim" thinks "y" increases down the screen */
     int scaleflag; /* if 1, scale image to 0..255 */
{
  FILE *fopen(),*fout;
  int i,j,k;
  unsigned char t;
  float min,max,***t3d;

  printf("  WRITE_HIPS_MOVIE\n");
  printf("    Writing %s\n",outfile);

  /*** THIS ROUTINE CAN'T HANDLE NON-ZERO START INDEXES ***/
  if ((x0 != 0)||(y0 != 0)||(z0 != 0))
    exit_error("WRITE_HIPS_MOVIE","Start index not 0");

  if (scaleflag){
    get_min_max_3d_farray(data,x0,xn,y0,yn,z0,zn,&min,&max);
    t3d = copy_3d_farray(data,x0,xn,y0,yn,z0,zn);
    x0 = y0 = z0 = 0;
    add_const_3d_farray(t3d,x0,xn,y0,yn,z0,zn,-min);
    multiply_3d_farray(t3d,x0,xn,y0,yn,z0,zn,255.0/(max-min));
  }else
    t3d = data;

  fout = fopen(outfile,"w");
  fprintf(fout,"HIPS\n");
  fprintf(fout,"WRITE_HIPS_MOVIE\n");
  fprintf(fout,"WYETH\n");
  fprintf(fout,"%d\n",zn);
  fprintf(fout,"00-00-00\n");
  fprintf(fout,"%d\n%d\n",xn,yn); /*** ROWS, COLS ***/
  fprintf(fout,"%d\n%d\n",xn,yn); /*** Rows, cols in region of interest. ***/
  fprintf(fout,"%d\n%d\n",0,0); /*** First row, first col in region of int ***/
  fprintf(fout,"%d\n",0); /*** Pixel format, 0==bytes ***/
  fprintf(fout,"%d\n",1); /*** numcolor ***/
  fprintf(fout,"%d\n\n",1); /*** sizehist ***/
  fprintf(fout,"%d\n\n",1); /*** sizedesc ***/
  fprintf(fout,"%d\n%d\n",0,0); /*** ??? ***/

  if (transpose){
    for (k=z0;k<(z0+zn);k++)
      if (yrevflag)
	for (j=(y0+yn-1);j>=y0;j--)
	  for (i=x0;i<(x0+xn);i++){
	    t = (unsigned char)t3d[i][j][k];
	    fwrite(&t,sizeof(unsigned char),1,fout);
	  }
      else
	for (j=y0;j<(y0+yn);j++)
	  for (i=x0;i<(x0+xn);i++){
	    t = (unsigned char)t3d[i][j][k];
	    fwrite(&t,sizeof(unsigned char),1,fout);
	  }
  }else{
    for (k=z0;k<(z0+zn);k++)
      for (i=x0;i<(x0+xn);i++)
	if (yrevflag)
	  for (j=(y0+yn-1);j>=y0;j--){
	    t = (unsigned char)t3d[i][j][k];
	    fwrite(&t,sizeof(unsigned char),1,fout);
	  }
	else
	  for (j=y0;j<(y0+yn);j++){
	    t = (unsigned char)t3d[i][j][k];
	    fwrite(&t,sizeof(unsigned char),1,fout);
	  }
  }
  fclose(fout);

  if (scaleflag)
    free_3d_farray(t3d,xn,yn,zn);
}
/**************************************-**************************************/
/*                                                                           */
/*                           CONVERT_2D_TO_TENSOR                            */
/*                                                                           */
/*  Put the 2D data into a 3D NumRecInC tensor.  The first dimension has     */
/*  length 1.                                                                */
/*                                                                           */
/*****************************************************************************/
float ***convert_2d_to_tensor(data,d2,d3)
     float **data;
     int d2,d3;
{
  int j,k;
  int d1;
  float ***tdata;

  d1 = 1;
  tdata = f3tensor(1,d1,1,d2,1,d3);
  for (j=1;j<=d2;j++)
    for (k=1;k<=d3;k++)
      tdata[1][j][k] = data[j-1][k-1];
  return tdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                           CONVERT_3D_TO_TENSOR                            */
/*                                                                           */
/*  The NumRecInC "tensor" data type is used for the "rlft3" routine.        */
/*                                                                           */
/*****************************************************************************/
float ***convert_3d_to_tensor(data,d1,d2,d3)
     float ***data;
     int d1,d2,d3;
{
  int i,j,k;
  float ***tdata;
  
  tdata = f3tensor(1,d1,1,d2,1,d3);
  for (i=1;i<=d1;i++)
    for (j=1;j<=d2;j++)
      for (k=1;k<=d3;k++)
	tdata[i][j][k] = data[i-1][j-1][k-1];
  return tdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                       GET_FRACTIONAL_CHANGE_XY_FARRAY                     */
/*                                                                           */
/*  Return the first (x,y) pair, and every pair that represents a            */
/*  change of greater than "f" fraction of the total range of the            */
/*  "ydata" relative to the last pair included.                              */
/*                                                                           */
/*  This routine aids in the creation of xy plots that have a sparse         */
/*  set of points.                                                           */
/*                                                                           */
/*****************************************************************************/
void get_fractional_change_xy_farray(xdata,ydata,n,f,rx,ry,rn)
     float *xdata,*ydata;
     int n;
     float f,**rx,**ry;
     int *rn;
{
  int i,k;
  float *x,*y,min,max,c;

  x = (float *)myalloc(n*sizeof(float));
  y = (float *)myalloc(n*sizeof(float));

  min = min_of_farray(ydata,n);
  max = max_of_farray(ydata,n);
  c = f*(max-min);

  k = 0;
  x[k] = xdata[0];
  y[k] = ydata[0];
  for(i=1;i<(n-1);i++){
    if (((ydata[i] - y[k]) > c)||((ydata[i] - y[k]) < -c)){
      k += 1;
      x[k] = xdata[i];
      y[k] = ydata[i];
    }
  }
  k += 1;
  x[k] = xdata[n-1];
  y[k] = ydata[n-1];
  k += 1;

  *rx = x; *ry = y; *rn = k;
}
/**************************************-**************************************/
/*                                                                           */
/*                        WRITE_FARRAY_SDEV_PLOT                             */
/*                                                                           */
/*****************************************************************************/
void write_farray_sdev_plot(outfile,data,sdev,n,xoffset)
     char outfile[];
     float *data,*sdev;
     int n,xoffset;
{
  FILE *fopen(),*fout;
  int i;

  fout = fopen(outfile,"w");
  for (i=0;i<n;i++)
    fprintf(fout,"%d %.4e %.4e\n",i+xoffset,data[i],sdev[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                       APPEND_FARRAY_SDEV_PLOT                             */
/*                                                                           */
/*****************************************************************************/
void append_farray_sdev_plot(outfile,name,data,sdev,n)
     char outfile[],name[];
     float *data,*sdev;
     int n;
{
  FILE *fopen(),*fout;
  int i;

  if ((fout = fopen(outfile,"a")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("APPEND_FARRAY_SDEV_PLOT","Cannot open file");
  }

  fprintf(fout,"/newplot\n");
  fprintf(fout,"/plotname %s\n",name);
  for (i=0;i<n;i++)
    fprintf(fout,"%d %.4e %.4e\n",i,data[i],sdev[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                    APPEND_FARRAY_SDEV_PLOT_OFFSET                         */
/*                                                                           */
/*****************************************************************************/
void append_farray_sdev_plot_offset(outfile,name,data,sdev,n,xoffset)
     char outfile[],name[];
     float *data,*sdev;
     int n,xoffset;
{
  FILE *fopen(),*fout;
  int i;

  /*fout = fopen(outfile,"a");*/
  if ((fout = fopen(outfile,"a")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("APPEND_FARRAY_SDEV_PLOT_OFFSET","Cannot open file");
  }


  fprintf(fout,"/newplot\n");
  fprintf(fout,"/plotname %s\n",name);
  for (i=0;i<n;i++)
    fprintf(fout,"%d %.4e %.4e\n",i+xoffset,data[i],sdev[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                             WRITE_FARRAY_XY_PLOT                          */
/*                                                                           */
/*****************************************************************************/
void write_farray_xy_plot(outfile,name,xdata,ydata,n)
     char outfile[],name[];
     float *xdata,*ydata;
     int n;
{
  FILE *fopen(),*fout;
  int i;

  /*fout = fopen(outfile,"w");*/
  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("WRITE_FARRAY_XY_PLOT","Cannot open file");
  }

  fprintf(fout,"/newplot\n");
  fprintf(fout,"/plotname %s\n",name);
  for (i=0;i<n;i++)
    fprintf(fout,"%.4e %.4e\n",xdata[i],ydata[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                         WRITE_FARRAY_XY_PLOT_D7                           */
/*                                                                           */
/*  Has 7 decimal places instead of the usual 4.                             */
/*                                                                           */
/*****************************************************************************/
void write_farray_xy_plot_d7(outfile,name,xdata,ydata,n)
     char outfile[],name[];
     float *xdata,*ydata;
     int n;
{
  FILE *fopen(),*fout;
  int i;

  /*fout = fopen(outfile,"w");*/
  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("WRITE_FARRAY_XY_PLOT","Cannot open file");
  }

  fprintf(fout,"/newplot\n");
  fprintf(fout,"/plotname %s\n",name);
  for (i=0;i<n;i++)
    fprintf(fout,"%.7e %.7e\n",xdata[i],ydata[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                           WRITE_FARRAY_XY_SEQ_PLOT                        */
/*                                                                           */
/*****************************************************************************/
void write_farray_xy_seq_plot(outfile,name,data,n,x0,dx)
     char outfile[],name[];
     float *data;
     int n;
     float x0,dx;
{
  FILE *fopen(),*fout;
  int i;

  /*fout = fopen(outfile,"w");*/
  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("WRITE_FARRAY_XY_SEQ_PLOT","Cannot open file");
  }

  fprintf(fout,"/newplot\n");
  fprintf(fout,"/plotname %s\n",name);
  for (i=0;i<n;i++)
    fprintf(fout,"%.4e %.4e\n",x0*(float)i*dx,data[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                       WRITE_FARRAY_XY_SDEV_PLOT                           */
/*                                                                           */
/*****************************************************************************/
void write_farray_xy_sdev_plot(outfile,xdata,ydata,sdev,n)
     char outfile[];
     float *xdata,*ydata,*sdev;
     int n;
{
  FILE *fopen(),*fout;
  int i;

  /*fout = fopen(outfile,"w");*/
  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("WRITE_FARRAY_XY_SDEV_PLOT","Cannot open file");
  }

  for (i=0;i<n;i++){
    fprintf(fout,"%.6e %.4e %.4e\n",xdata[i],ydata[i],sdev[i]);
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                         APPEND_FARRAY_XY_SDEV_PLOT                        */
/*                                                                           */
/*****************************************************************************/
void append_farray_xy_sdev_plot(outfile,name,xdata,ydata,sdev,n)
     char outfile[],name[];
     float *xdata,*ydata,*sdev;
     int n;
{
  FILE *fopen(),*fout;
  int i;

  /*fout = fopen(outfile,"a");*/
  if ((fout = fopen(outfile,"a")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("APPEND_FARRAY_XY_SDEV_PLOT","Cannot open file");
  }

  fprintf(fout,"/newplot\n");
  fprintf(fout,"/plotname %s\n",name);
  for (i=0;i<n;i++)
    fprintf(fout,"%.4e %.4e %.4e\n",xdata[i],ydata[i],sdev[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                            APPEND_FARRAY_XY_PLOT                          */
/*                                                                           */
/*****************************************************************************/
void append_farray_xy_plot(outfile,xdata,ydata,n,name)
     char outfile[];
     float *xdata,*ydata;
     int n;
     char name[];
{
  FILE *fopen(),*fout;
  int i;

  if ((fout = fopen(outfile,"a")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("APPEND_FARRAY_XY_PLOT","Cannot open file");
  }

  fprintf(fout,"/newplot\n");
  fprintf(fout,"/plotname %s\n",name);
  for (i=0;i<n;i++){
    fprintf(fout,"%.4f %.4e\n",xdata[i],ydata[i]);
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                        APPEND_FARRAY_XY_PLOT_NONEW                        */
/*                                                                           */
/*  If there is no 'newplot' in a file with only one plot, xplot will not    */
/*  list its name twice.                                                     */
/*                                                                           */
/*****************************************************************************/
void append_farray_xy_plot_nonew(outfile,xdata,ydata,n,name)
     char outfile[];
     float *xdata,*ydata;
     int n;
     char name[];
{
  FILE *fopen(),*fout;
  int i;

  /*fout = fopen(outfile,"a");*/
  if ((fout = fopen(outfile,"a")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("APPEND_FARRAY_XY_PLOT_NONEW","Cannot open file");
  }

  fprintf(fout,"/plotname %s\n",name);
  for (i=0;i<n;i++)
    fprintf(fout,"%.4f %.4e\n",xdata[i],ydata[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                         APPEND_FARRAY_XY_PLOT_OFFSET                      */
/*                                                                           */
/*****************************************************************************/
void append_farray_xy_plot_offset(outfile,name,xdata,ydata,n,xoffset)
     char outfile[],name[];
     float *xdata,*ydata;
     int n;
     float xoffset;
{
  FILE *fopen(),*fout;
  int i;

  /*fout = fopen(outfile,"a");*/
  if ((fout = fopen(outfile,"a")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("APPEND_FARRAY_XY_PLOT_OFFSET","Cannot open file");
  }

  fprintf(fout,"/newplot\n");
  fprintf(fout,"/plotname %s\n",name);
  for (i=0;i<n;i++)
    fprintf(fout,"%.4f %.4e\n",xdata[i]+xoffset,ydata[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                          APPEND_FARRAY_XY_SEQ_PLOT                        */
/*                                                                           */
/*****************************************************************************/
void append_farray_xy_seq_plot(outfile,name,data,n,x0,dx)
     char outfile[],name[];
     float *data;
     int n;
     float x0,dx;
{
  FILE *fopen(),*fout;
  int i;

  /*fout = fopen(outfile,"a");*/
  if ((fout = fopen(outfile,"a")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("APPEND_FARRAY_XY_SEQ_PLOT","Cannot open file");
  }

  fprintf(fout,"/newplot\n");
  fprintf(fout,"/plotname %s\n",name);
  for (i=0;i<n;i++)
    fprintf(fout,"%.4f %.4e\n",x0+(float)i*dx,data[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                            WRITE_FARRAY_PLOT                              */
/*                                                                           */
/*****************************************************************************/
void write_farray_plot(outfile,data,n)
     char outfile[];
     float *data;
     int n;
{
  FILE *fopen(),*fout;
  int i;

  /*fout = fopen(outfile,"w");*/
  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("WRITE_FARRAY_PLOT","Cannot open file");
  }

  for (i=0;i<n;i++)
    fprintf(fout,"%d %.4e\n",i,data[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                         WRITE_FARRAY_PLOT_NAME                            */
/*                                                                           */
/*****************************************************************************/
void write_farray_plot_name(outfile,name,data,n)
     char outfile[],name[];
     float *data;
     int n;
{
  FILE *fopen(),*fout;
  int i;

  /*fout = fopen(outfile,"w");*/
  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("WRITE_FARRAY_PLOT_NAME","Cannot open file");
  }

  fprintf(fout,"/newplot\n");
  fprintf(fout,"/plotname %s\n",name);
  for (i=0;i<n;i++)
    fprintf(fout,"%d %.4e\n",i,data[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                         WRITE_FARRAY_SCALE_PLOT                           */
/*                                                                           */
/*****************************************************************************/
void write_farray_scale_plot(outfile,data,n,xscale,yscale)
     char outfile[];
     float *data;
     int n;
     float xscale,yscale;
{
  FILE *fopen(),*fout;
  int i;

  /*fout = fopen(outfile,"w");*/
  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("WRITE_FARRAY_SCALE_PLOT","Cannot open file");
  }

  for (i=0;i<n;i++)
    fprintf(fout,"%.4e %.2e\n",(float)i*xscale,data[i]*yscale);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                        WRITE_FARRAY_PLOT_OFFSET                           */
/*                                                                           */
/*****************************************************************************/
void write_farray_plot_offset(outfile,data,n,xoffset)
     char outfile[];
     float *data;
     int n,xoffset;
{
  FILE *fopen(),*fout;
  int i;

  /*fout = fopen(outfile,"w");*/
  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("WRITE_FARRAY_PLOT_OFFSET","Cannot open file");
  }

  for (i=0;i<n;i++)
    fprintf(fout,"%d %.4e\n",i+xoffset,data[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                       APPEND_FARRAY_PLOT_OFFSET                           */
/*                                                                           */
/*****************************************************************************/
void append_farray_plot_offset(outfile,name,data,n,xoffset)
     char outfile[],name[];
     float *data;
     int n,xoffset;
{
  FILE *fopen(),*fout;
  int i;

  /*fout = fopen(outfile,"a");*/
  if ((fout = fopen(outfile,"a")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("APPEND_FARRAY_PLOT_OFFSET","Cannot open file");
  }

  fprintf(fout,"/newplot\n");
  fprintf(fout,"/plotname %s\n",name);
  for (i=0;i<n;i++)
    fprintf(fout,"%d %.4e\n",i+xoffset,data[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_HISTOGRAM_FARRAY                          */
/*                                                                           */
/*  Make a histogram of the data, return it and the number of points out of  */
/*  bounds.                                                                  */
/*                                                                           */
/*****************************************************************************/
float *get_histogram_farray(data,n,min,max,nbin,rbinsize,rnpob)
     float *data;
     int n;
     float min,max;
     int nbin;
     float *rbinsize;
     int *rnpob;  // Number of Points Out of Bounds
{
  int i,k;
  int npob;
  float *hist,binsize;

  binsize = (max-min)/(float)nbin;
  hist = get_zero_farray(nbin);
  npob = 0;
  for (i=0;i<n;i++){
    k = (data[i]-min)/binsize;
    if ((k >= 0)&&(k < nbin))
      hist[k] += 1.0;
    else{
      //printf("(farray_util.c) GET_HISTOGRAM_FARRAY  data[%d] = %f\n",i,
      // data[i]);
      npob += 1;
    }
  }
  *rnpob = npob; *rbinsize = binsize;
  return hist;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_HISTOGRAM_AUTO_FARRAY                        */
/*                                                                           */
/*  Make a histogram of the data, return it and the number of points out of  */
/*  bounds.                                                                  */
/*                                                                           */
/*****************************************************************************/
float *get_histogram_auto_farray(data,n,nbin,rmin,rmax,rbinsize)
     float *data;         // Values to be histogrammed [n]
     int n;               // Number of values
     int nbin;            // Number of bins to use in histogram
     float *rmin,*rmax;   // Return min and max bounds of histogram
     float *rbinsize;     // Return binsize
{
  int npob;
  float dmin,dmax,drange,*hh,binsize;

  get_min_max_farray(data,n,&dmin,&dmax);
  if (dmin == dmax)
    exit_error("GET_HISTOGRAM_AUTO_FARRAY","dmin = dmax");

  drange = dmax - dmin;
  dmax = dmax + 0.1*drange;
  dmin = dmin - 0.1*drange;
  hh = get_histogram_farray(data,n,dmin,dmax,nbin,&binsize,&npob);
  if (npob != 0)
    exit_error("GET_HISTOGRAM_AUTO_FARRAY","Points are out of bounds");

  *rmin = dmin;
  *rmax = dmax;
  *rbinsize = binsize;

  return hh;
}
/**************************************-**************************************/
/*                                                                           */
/*                            WRITE_FARRAY_HIST_PLOT                         */
/*                                                                           */
/*  WYETH --- THIS SHOULD GO INTO "plot_util".                               */
/*  WYETH --- THIS SHOULD USE GET_HISTOGRAM_FARRAY.                          */
/*                                                                           */
/*****************************************************************************/
void write_farray_hist_plot(outfile,data,n,min,max,nbin,npob)
     char outfile[];
     float *data;
     int n;
     float min,max;
     int nbin;
     int *npob;  /* Number of Points Out of Bounds */
{
  FILE *fopen(),*fout;
  int i;
  int k;
  float *hist,binsize;

  binsize = (max-min)/(float)nbin;
  hist = get_zero_farray(nbin);
  *npob = 0;
  for (i=0;i<n;i++){
    k = (data[i]-min)/binsize;
    if ((k >= 0)&&(k < nbin))
      hist[k] += 1.0;
    else
      *npob += 1;
  }

  /*fout = fopen(outfile,"w");*/
  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("WRITE_FARRAY_HIST_PLOT","Cannot open file");
  }

  fprintf(fout,"%.4f %.4f\n",min,0.0);
  for (i=0;i<nbin;i++){
    fprintf(fout,"%.4f %.4f\n",(float)i*binsize+min,hist[i]);
    fprintf(fout,"%.4f %.4f\n",(float)(i+1)*binsize+min,hist[i]);
    fprintf(fout,"%.4f %.4f\n",(float)(i+1)*binsize+min,0.0);
  }
  fclose(fout);
  myfree(hist);
}
/**************************************-**************************************/
/*                                                                           */
/*                             APPEND_FARRAY_PLOT                            */
/*                                                                           */
/*****************************************************************************/
void append_farray_plot(outfile,name,data,n,xscale)
     char outfile[],name[];
     float *data;
     int n,xscale;
{
  FILE *fopen(),*fout;
  int i;

  if ((fout = fopen(outfile,"a")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("APPEND_FARRAY_PLOT","Cannot open file");
  }

  if (name != NULL){
    fprintf(fout,"/newplot\n");
    fprintf(fout,"/plotname %s\n",name);
  }
  for (i=0;i<n;i++)
    fprintf(fout,"%d %.3e\n",i*xscale,data[i]);

  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                         APPEND_FARRAY_PLOT_XSCALE                         */
/*                                                                           */
/*****************************************************************************/
void append_farray_plot_xscale(outfile,name,data,n,xscale)
     char outfile[],name[];
     float *data;
     int n;
     float xscale;
{
  FILE *fopen(),*fout;
  int i;

  if ((fout = fopen(outfile,"a")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("APPEND_FARRAY_PLOT_XSCALE","Cannot open file");
  }

  if (name != NULL){
    fprintf(fout,"/newplot\n");
    fprintf(fout,"/plotname %s\n",name);
  }
  for (i=0;i<n;i++)
    fprintf(fout,"%.3e %.3e\n",(float)i*xscale,data[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                          APPEND_FARRAY_PLOT_COLOR                         */
/*                                                                           */
/*  Colors                                                                   */
/*    0 black                                                                */
/*    1 red                                                                  */
/*    2 bright green                                                         */
/*    3 yellow                                                               */
/*    4 blue (hard to focus on)                                              */
/*    5 purple                                                               */
/*    6 light blue                                                           */
/*    7 white                                                                */
/*    8 dark gray (hard to see)                                              */
/*    9 brownish/purple                                                      */
/*   10 pale green                                                           */
/*   11 olive green                                                          */
/*   12 pale blue                                                            */
/*  ...                                                                      */
/*  255                                                                      */
/*                                                                           */
/*  Pointstyles                                                              */
/*   0 open square                                                           */
/*   1 fill square                                                           */
/*   2 open circle                                                           */
/*   3 open circle large                                                     */
/*   4 open triangle                                                         */
/*   5 fill triangle                                                         */
/*   6 x                                                                     */
/*   7 +                                                                     */
/*   8 .                                                                     */
/*   9 nothing                                                               */
/*                                                                           */
/*****************************************************************************/
void append_farray_plot_color(outfile,name,data,n,xoffset,xscale,color,pstyle)
     char outfile[],name[];
     float *data;
     int n,xoffset,xscale,color,pstyle;
{
  FILE *fopen(),*fout;
  int i;

  if ((fout = fopen(outfile,"a")) == NULL){
    printf("  *** Filename %s\n",outfile);
    exit_error("APPEND_FARRAY_PLOT_COLOR","Cannot open file");
  }

  if (name != NULL){
    fprintf(fout,"/newplot\n");
    fprintf(fout,"/plotname %s\n",name);
  }
  if (color >= 0)
    fprintf(fout,"/color %d\n",color);
  if (pstyle >= 0)
    fprintf(fout,"/pointstyle %d\n",pstyle);
  for (i=0;i<n;i++)
    fprintf(fout,"%d %.3e\n",xoffset + i*xscale,data[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                         FREE_FARRAY_PLOT_FILE_OLD                         */
/*                                                                           */
/*  Free storage from 'read_farray_plot_file_old'.                           */
/*                                                                           */
/*****************************************************************************/
void free_farray_plot_file_old(xdata,ydata,zdata,name,m,n)
     float **xdata,**ydata,**zdata;
     char **name;
     int m,n;
{
  free_2d_farray(xdata,m);
  free_2d_farray(ydata,m);
  free_2d_farray(zdata,m);
  free_2d_carray(name,m);
}
/**************************************-**************************************/
/*                                                                           */
/*                       READ_FARRAY_PLOT_FILE_OLD                           */
/*                                                                           */
/*  **** Doesn't handle variable length plots in one file.                   */
/*                                                                           */
/*  Read a sequence of float pairs (or triples) separated by lines           */
/*  beginning with "/".                                                      */
/*                                                                           */
/*  FOR EXAMPLE:                                                             */
/*                                                                           */
/*    /newplot                                                               */
/*    /plotname zz.plot.1                                                    */
/*    1.0 2.34e+01 23.42                                                     */
/*    ...                                                                    */
/*    100.0 8.09e-03 203.32                                                  */
/*    /newplot                                                               */
/*    /plotname zz.plot.2                                                    */
/*    1.0 0.88e+01 230.1                                                     */
/*    ...                                                                    */
/*                                                                           */
/*****************************************************************************/
void read_farray_plot_file_old(infile,rxdata,rydata,rzdata,rname,m,n)
     char infile[];
     float ***rxdata,***rydata,***rzdata;
     char ***rname;
     int *m,*n;
{
  int j,k;
  int flag,plot_count,len;
  FILE *fopen(),*fin;
  char line[SLEN],t1[SLEN],t2[SLEN],**name;
  float **xdata,**ydata,**zdata;

  printf("  READ_FARRAY_PLOT_FILE_OLD\n");
  printf("    Reading %s.\n",infile);

  if ((fin = fopen(infile,"r")) != NULL){
    len = plot_count = 0;
    flag = 1;
    while(fgets(line,SLEN,fin) != NULL){
      if ((line[0] == '/')&&(line[0] != '\n')){
	flag = 1;
      }else
	if (flag == 1){
	  plot_count += 1;
	  flag = 0;
	  len = 1;
	}else
	  len += 1;
    }
    fclose(fin);
    printf("    %d plots found.  Plot length = %d.\n",plot_count,len);
    
    xdata = get_2d_farray(plot_count,len);
    ydata = get_2d_farray(plot_count,len);
    zdata = get_2d_farray(plot_count,len);
    name = get_2d_carray(plot_count,SLEN); /* Def'd. in "my_util.c" */
    fin = fopen(infile,"r");
    flag = 1;
    j = -1;
    k = 0;
    while(fgets(line,SLEN,fin) != NULL){
      if ((line[0] == '/')&&(line[0] != '\n')){
	flag = 1;
	sscanf(line,"%1s %8s",t1,t2);
	if (strcmp(t2,"plotname")==0)
	  sscanf(line,"%1s %8s %s",t1,t2,name[j+1]);
      }else{
	if (flag == 1){
	  j += 1;
	  k = 0;
	  flag = 0;
	}
	sscanf(line,"%f %f %f\n",&xdata[j][k],&ydata[j][k],&zdata[j][k]);
	k += 1;
      }
    }
    fclose(fin);

    *rxdata = xdata;
    *rydata = ydata;
    *rzdata = zdata;
    *rname = name;
    *m = plot_count;
    *n = len;
  }else{
    *m = -1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                READ_FARRAY                                */
/*                                                                           */
/*****************************************************************************/
int read_farray(infile,data,n)
     char infile[];
     float **data;
     int *n;
{
  FILE *fopen(),*fin;
  int i;
  int ns;
  float temp;

  if ((fin = fopen(infile,"r"))==NULL)
    return 0;

  *n = 0;
  while (fscanf(fin,"%f",&temp)!= EOF)
    *n += 1;
  fclose(fin);

  *data = (float *)myalloc(*n*sizeof(float));

  fin = fopen(infile,"r");
  for(i=0;i<*n;i++)
    ns = fscanf(fin,"%f",&(*data)[i]);
  fclose(fin);

  return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                        READ_FARRAY_PTR_UNTIL_STRING                       */
/*                                                                           */
/*****************************************************************************/
void read_farray_ptr_until_string(fin,estr,nmax,rdata,rn)
     FILE *fin;       // Read from this file pointer
     char *estr;      // Until this string is found; NULL for EOF
     int nmax;        // Maximum number of data values
     float **rdata;   // [*rn] Returned data
     int *rn;         // Number of values returned
{
  int n,done;
  float *tdata;
  char tstr[SLEN];

  tdata = (float *)myalloc(nmax*sizeof(float));

  n = 0;

  if (estr == NULL){
    //
    //  Read until EOF
    //
    while (fscanf(fin,"%s",tstr) != EOF){
      tdata[n] = atof(tstr);
      n += 1;
    }
  }else{
    //
    //  Read until the string 'estr' is found
    //
    done = 0;
    while (done == 0){
      if (fscanf(fin,"%s",tstr) == EOF){
	done = 1;
	printf("  ***  End string:  %s\n",estr);
	exit_error("READ_FARRAY_PTR_UNTIL_STRING","End string not found");
      }else{
	if (strcmp(tstr,estr)==0){
	  done = 1;
	}else{
	  tdata[n] = atof(tstr);
	  n += 1;
	}
      }
    }
  }

  if (n == 0){
    *rdata = NULL;
  }else{
    *rdata = copy_farray(tdata,n);
  }

  *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                               READ_2D_FARRAY                              */
/*                                                                           */
/*****************************************************************************/
int read_2d_farray(infile,n2,data,rn1)
     char infile[];  // File name
     int n2;         // Each farray will have this many elements
     float ***data;  // [*n1][n2]
     int *rn1;       // The number of farrays will be determined
{
  FILE *fopen(),*fin;
  int i,j;
  int n,n1,ns;
  float temp;

  if ((fin = fopen(infile,"r"))==NULL)
    return 0;

  n = 0;
  while (fscanf(fin,"%f",&temp)!= EOF)
    n += 1;
  fclose(fin);

  n1 = n/n2;
  if (n1*n2 != n){
    printf("  *** %d floats found, must be a multiple of %d\n",n,n2);
    exit_error("READ_2D_FARRAY","File size is not a multiple of 'n2'\n");
  }

  *data = get_2d_farray(n1,n2);

  fin = fopen(infile,"r");
  for(i=0;i<n1;i++)
    for(j=0;j<n2;j++)
      ns = fscanf(fin,"%f",&(*data)[i][j]);
  fclose(fin);

  *rn1 = n1;

  return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                              READ_2COL_FARRAY                             */
/*                                                                           */
/*****************************************************************************/
int read_2col_farray(infile,xdata,ydata,n)
     char infile[];
     float **xdata,**ydata;
     int *n;
{
  FILE *fopen(),*fin;
  int i;
  int ns;
  float temp;

  if ((fin = fopen(infile,"r"))==NULL)
    return 0;

  *n = 0;
  while (fscanf(fin,"%f %f",&temp,&temp)!= EOF)
    *n += 1;
  fclose(fin);

  *xdata = (float *)myalloc(*n*sizeof(float));
  *ydata = (float *)myalloc(*n*sizeof(float));

  fin = fopen(infile,"r");
  for(i=0;i<*n;i++)
    ns = fscanf(fin,"%f %f",&(*xdata)[i],&(*ydata)[i]);
  fclose(fin);

  return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                                WRITE_FARRAY                               */
/*                                                                           */
/*****************************************************************************/
void write_farray(outfile,data,n)
     char outfile[];
     float *data;
     int n;
{
  FILE *fopen(),*fout;
  int i;
  
  if ((fout = fopen(outfile,"w"))==NULL){
    printf("  *** File %s\n",outfile);
    exit_error("WRITE_FARRAY","Can't open file");
  }
  for(i=0;i<n;i++)
    fprintf(fout,"%f\n",data[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                               APPEND_FARRAY                               */
/*                                                                           */
/*****************************************************************************/
void append_farray(outfile,data,n)
     char outfile[];
     float *data;
     int n;
{
  FILE *fopen(),*fout;
  int i;
  
  if ((fout = fopen(outfile,"a"))==NULL){
    printf("  *** File %s\n",outfile);
    exit_error("APPEND_FARRAY","Can't open file");
  }
  for(i=0;i<n;i++)
    fprintf(fout,"%f\n",data[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                             APPEND_LINE_FARRAY                            */
/*                                                                           */
/*****************************************************************************/
void append_line_farray(outfile,data,n)
     char outfile[];
     float *data;
     int n;
{
  FILE *fopen(),*fout;
  int i;
  
  if ((fout = fopen(outfile,"a"))==NULL){
    printf("  *** File %s\n",outfile);
    exit_error("APPEND_LINE_FARRAY","Can't open file");
  }

  if (n > 0){
    for(i=0;i<(n-1);i++)
      fprintf(fout,"%f ",data[i]);
    fprintf(fout,"%f\n",data[n-1]);
  }else{
    fprintf(fout,"\n");
  }


  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                             WRITE_FARRAY_FORMAT                           */
/*                                                                           */
/*  User specifies the printf format string.                                 */
/*                                                                           */
/*****************************************************************************/
void write_farray_format(outfile,data,n,pstr)
     char outfile[];
     float *data;
     int n;
     char pstr[];
{
  FILE *fopen(),*fout;
  int i;
  
  if ((fout = fopen(outfile,"w"))==NULL){
    printf("  *** File %s\n",outfile);
    exit_error("WRITE_FARRAY_FORMAT","Can't open file");
  }
  for(i=0;i<n;i++)
    fprintf(fout,pstr,data[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                            APPEND_FARRAY_FORMAT                           */
/*                                                                           */
/*  User specifies the printf format string.                                 */
/*                                                                           */
/*****************************************************************************/
void append_farray_format(outfile,data,n,pstr)
     char outfile[];
     float *data;
     int n;
     char pstr[];
{
  FILE *fopen(),*fout;
  int i;
  
  if ((fout = fopen(outfile,"a"))==NULL){
    printf("  *** File %s\n",outfile);
    exit_error("APPEND_FARRAY_FORMAT","Can't open file");
  }
  for(i=0;i<n;i++)
    fprintf(fout,pstr,data[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                               WRITE_2COL_FARRAY                           */
/*                                                                           */
/*****************************************************************************/
void write_2col_farray(outfile,data1,data2,n)
     char outfile[];
     float *data1,*data2;
     int n;
{
  FILE *fopen(),*fout;
  int i;
  
  if ((fout = fopen(outfile,"w"))==NULL){
    printf("  *** File %s\n",outfile);
    exit_error("WRITE_2COL_FARRAY","Can't open file");
  }
  for(i=0;i<n;i++)
    fprintf(fout,"%f %f\n",data1[i],data2[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                              APPEND_2COL_FARRAY                           */
/*                                                                           */
/*****************************************************************************/
void append_2col_farray(outfile,data1,data2,n)
     char outfile[];
     float *data1,*data2;
     int n;
{
  FILE *fopen(),*fout;
  int i;
  
  if ((fout = fopen(outfile,"a"))==NULL){
    printf("  *** File %s\n",outfile);
    exit_error("WRITE_2COL_FARRAY","Can't open file");
  }
  for(i=0;i<n;i++)
    fprintf(fout,"%f %f\n",data1[i],data2[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                               WRITE_3COL_FARRAY                           */
/*                                                                           */
/*****************************************************************************/
void write_3col_farray(outfile,data1,data2,data3,n)
     char outfile[];
     float *data1,*data2,*data3;
     int n;
{
  FILE *fopen(),*fout;
  int i;
  
  if ((fout = fopen(outfile,"w"))==NULL){
    printf("  *** File %s\n",outfile);
    exit_error("WRITE_3COL_FARRAY","Can't open file");
  }
  for(i=0;i<n;i++)
    fprintf(fout,"%f %f %f\n",data1[i],data2[i],data3[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                              APPEND_2D_FARRAY                             */
/*                                                                           */
/*****************************************************************************/
void append_2d_farray(outfile,data,n1,n2,transpose)
     char outfile[];
     float **data;
     int n1,n2;
     int transpose;
{
  FILE *fopen(),*fout;
  int i,j;
  
  if ((fout = fopen(outfile,"a"))==NULL){
    printf("  *** File %s\n",outfile);
    exit_error("WRITE_2COL_FARRAY","Cannot open file to append");
  }
  if (transpose == 0){
    for(i=0;i<n1;i++){
      fprintf(fout,"%f",data[i][0]);
      for(j=1;j<n2;j++){
	fprintf(fout," %f",data[i][j]);
      }
      fprintf(fout,"\n");
    }
  }else{
    for(j=0;j<n2;j++){
      fprintf(fout,"%f",data[0][j]);
      for(i=1;i<n1;i++){
	fprintf(fout," %f",data[i][j]);
      }
      fprintf(fout,"\n");
    }
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                               WRITE_2D_FARRAY                             */
/*                                                                           */
/*****************************************************************************/
void write_2d_farray(outfile,data,n1,n2,transpose)
     char outfile[];
     float **data;
     int n1,n2;
     int transpose;
{
  remove_file(outfile);
  append_2d_farray(outfile,data,n1,n2,transpose);
}
/**************************************-**************************************/
/*                                                                           */
/*                             ADD_FARRAY_TO_FILE                            */
/*                                                                           */
/*****************************************************************************/
void add_farray_to_file(outfile,data,n)
     char outfile[];
     float *data;
     int n;
{
  int i;
  int ndata;
  float *fdata;

  if (read_farray(outfile,&fdata,&ndata) == 1){
    if (ndata != n)
      exit_error("ADD_FARRAY_TO_FILE","Stored array differs");
    for(i=0;i<n;i++)
      fdata[i] += data[i];
  }else
    fdata = data;
  
  write_farray(outfile,fdata,n);
}
/**************************************-**************************************/
/*                                                                           */
/*                               MEAN_ABS_FARRAY                             */
/*                                                                           */
/*  Mean of the absolute value of the farray.                                */
/*                                                                           */
/*****************************************************************************/
float mean_abs_farray(data,n)
     float *data;
     int n;
{
  int i;
  float m;
  
  m = 0.0;
  for(i=0;i<n;i++)
    m += fabs(data[i]);
  if (n>0)
    m /= (float)n;
  
  return m;
}
/**************************************-**************************************/
/*                                                                           */
/*                                  MEAN_FARRAY                              */
/*                                                                           */
/*  Compute the mean of the farray.                                          */
/*                                                                           */
/*****************************************************************************/
float mean_farray(data,n)
     float *data;
     int n;
{
  int i;
  float m;
  
  m = 0.0;
  for(i=0;i<n;i++)
    m += data[i];
  if (n>0)
    m /= (float)n;
  
  return m;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MEAN_SDEV_FARRAY                            */
/*                                                                           */
/*  Take an array of floats and compute the mean and standard deviation.     */
/*                                                                           */
/*****************************************************************************/
void mean_sdev_farray(data,n,mean,sdev)
     float *data;
     int n;
     float *mean,*sdev;
{
  int i;
  float m,s;

  m = 0.0;
  for(i=0;i<n;i++)
    m += data[i];
  if (n>0) m /= (float)n;

  s = 0.0;
  for(i=0;i<n;i++)
    s += ((data[i]-m)*(data[i]-m));
  if (n>1) s = sqrt(s/(float)(n-1));

  *mean = m; *sdev = s;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MEAN_ADEV_FARRAY                            */
/*                                                                           */
/*  Take an array of floats and compute the mean and the mean of the         */
/*  absolute deviation from the mean.                                        */
/*                                                                           */
/*****************************************************************************/
void mean_adev_farray(data,n,mean,adev)
     float *data;
     int n;
     float *mean,*adev;
{
  int i;
  float m,s;

  if (n < 1)
    exit_error("MEAN_ADEV_FARRAY","n < 1");

  m = 0.0;
  for(i=0;i<n;i++)
    m += data[i];
  if (n>0) m /= (float)n;

  s = 0.0;
  for(i=0;i<n;i++)
    s += fabs(data[i]-m);

  *mean = m; *adev = s;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MEAN_SDEV_DARRAY                            */
/*                                                                           */
/*  Take an array of double and compute the mean and standard deviation.     */
/*                                                                           */
/*****************************************************************************/
void mean_sdev_darray(data,n,mean,sdev)
     double *data;
     int n;
     double *mean,*sdev;
{
  int i;
  double m,s;

  m = 0.0;
  for(i=0;i<n;i++)
    m += data[i];
  if (n>0) m /= (double)n;

  s = 0.0;
  for(i=0;i<n;i++)
    s += ((data[i]-m)*(data[i]-m));
  if (n>1) s = sqrt(s/(double)(n-1));

  *mean = m; *sdev = s;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MEAN_SDEV_SKEW_FARRAY                         */
/*                                                                           */
/*  Compute mean, std. dev., and skewness.                                   */
/*                                                                           */
/*****************************************************************************/
void mean_sdev_skew_farray(data,n,mean,sdev,skew)
     float *data;
     int n;
     float *mean,*sdev,*skew;
{
  int i;
  float m,s,t,sk;
  
  m = 0.0;
  for(i=0;i<n;i++)
    m += data[i];
  if (n>0)
    m /= (float)n;
  
  s = 0.0;
  for(i=0;i<n;i++)
    s += ((data[i]-m)*(data[i]-m));
  if (n>1)
    s = sqrt(s/(float)(n-1));

  t = 0.0; /*** WYETH, CHECK THIS FORMULA ***/
  for(i=0;i<n;i++)
    t += ((data[i]-m)*(data[i]-m)*(data[i]-m));
  if (n>1)
    sk = t/(float)(n-1)/(s*s*s);
  else
    sk = 0.0;
  
  *mean = m; *sdev = s; *skew = sk;
}
/**************************************-**************************************/
/*                                                                           */
/*                                SKEW_FARRAY                                */
/*                                                                           */
/*  Compute exact skewness of farray.                                        */
/*                                                                           */
/*****************************************************************************/
float skew_farray(data,n)
     float *data;
     int n;
{
  int i;
  float m,d,m2,m3,sk;
  
  m = 0.0;
  for(i=0;i<n;i++)
    m += data[i];
  if (n>0)
    m /= (float)n;
  
  m2 = 0.0;
  m3 = 0.0;
  for(i=0;i<n;i++){
    d = data[i]-m;
    m2 += d*d;
    m3 += d*d*d;
  }
  if (n>1){
    m2 /= (float)n;
    m3 /= (float)n;
    sk = m3 / sqrt(m2*m2*m2);
  }

  printf("  mean = %f\n",m);

  return sk;
}
/**************************************-**************************************/
/*                                                                           */
/*                        WEIGHTED_MEAN_SDEV_FARRAY                          */
/*                                                                           */
/*  Take an array of floats and compute the mean and standard                */
/*  deviation.                                                               */
/*                                                                           */
/*  *** Warning:  multiplying by 1/(n-1) has been replaced with              */
/*                multiplying by n/(n-1) * 1/total_weight.                   */
/*                                                                           */
/*****************************************************************************/
void weighted_mean_sdev_farray(data,n,weight,mean,sdev)
     float *data;
     int n;
     float *weight,*mean,*sdev;
{
  int i;
  float m,s,total_weight;
  
  m = 0.0;
  total_weight = 0.0;
  for(i=0;i<n;i++){
    m += weight[i] * data[i];
    total_weight += weight[i];
  }
  if (n>0)
    m /= (float)total_weight;
  
  s = 0.0;
  for(i=0;i<n;i++)
    s += weight[i] * ((data[i]-m)*(data[i]-m));
  if (n>1)
    s = sqrt(s*(float)n/(float)(n-1) / total_weight);
  
  *mean = m; *sdev = s;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MEAN_TWO_SIDED_SDEV_FARRAY                      */
/*                                                                           */
/*  Compute a "standard deviation" separately for the left and right sides   */
/*  of the mean.  Count values at the mean as 1/2 weight to each side.       */
/*                                                                           */
/*****************************************************************************/
void mean_two_sided_sdev_farray(data,n,mean,rsdl,rsdr)
     float *data;
     int n;
     float *mean,*rsdl,*rsdr;
{
  int i;
  float m,sl,sr,dev,nr,nl;
  
  m = 0.0;
  for(i=0;i<n;i++)
    m += data[i];
  if (n>0) m /= (float)n;
  
  nr = nl = 0.0;
  sr = sl = 0.0;
  for(i=0;i<n;i++){
    dev = data[i] - m;
    if (dev > 0.0){
      sr += dev*dev;
      nr += 1.0;
    }else if (dev < 0.0){
      sl += dev*dev;
      nl += 1.0;
    }else{
      /* dev == 0.0, so no need to increment sr, sl. */
      nr += 0.5;
      nl += 0.5;
    }
  }
  if (nl>1)
    sl = sqrt(sl/(nl-1.0));
  if (nr>1)
    sr = sqrt(sr/(nr-1.0));
  
  *mean = m; *rsdl = sl; *rsdr = sr;
}
/**************************************-**************************************/
/*                                                                           */
/*                            SUBTRACT_MEAN_FARRAY                           */
/*                                                                           */
/*  Take an array of floats, subtract the mean.  This used to be called      */
/*  "make_mean_zero_farray".                                                 */
/*                                                                           */
/*****************************************************************************/
void subtract_mean_farray(data,n,oldmean)
     float *data;
     int n;
     float *oldmean;
{
  int i;
  float mean,sdev;
  
  mean_sdev_farray(data,n,&mean,&sdev);
  for(i=0;i<n;i++)
    data[i] -= mean;
  *oldmean = mean;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_ZERO_MEAN_FARRAY                          */
/*                                                                           */
/*****************************************************************************/
float *get_zero_mean_farray(data,n)
     float *data;
     int n;
{
  float *zd,oldmean;
  
  zd = copy_farray(data,n);
  subtract_mean_farray(zd,n,&oldmean);
  return zd;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MAKE_MAX_CONST_FARRAY                         */
/*                                                                           */
/*  Rescale the farray so that its maximum value is "c".                     */
/*                                                                           */
/*****************************************************************************/
void make_max_const_farray(data,n,c)
     float *data;
     int n;
     float c;
{
  float max;
  
  max = max_of_farray(data,n);
  if (max == 0.0)
    exit_error("MAKE_MAX_CONST_FARRAY","Max is zero");
  multiply_farray(data,n,c/max);
}
/**************************************-**************************************/
/*                                                                           */
/*                          SCALE_MEAN_TO_VALUE_FARRAY                       */
/*                                                                           */
/*  Rescale the farray so that its mean value is "c".                        */
/*                                                                           */
/*****************************************************************************/
void scale_mean_to_value_farray(data,n,c)
     float *data;
     int n;
     float c;
{
  float mean;

  mean = mean_farray(data,n);
  if (mean == 0.0)
    exit_error("SCALE_MEAN_TO_VALUE_FARRAY","Mean is zero");
  multiply_farray(data,n,c/mean);
}
/**************************************-**************************************/
/*                                                                           */
/*                               Z_SCORE_FARRAY                              */
/*                                                                           */
/*  Take an array of floats, subtract the mean, divide by std.dev.           */
/*                                                                           */
/*****************************************************************************/
void z_score_farray(data,n)
     float *data;
     int n;
{
  int i;
  float mean,sdev;
  
  mean_sdev_farray(data,n,&mean,&sdev);
  for(i=0;i<n;i++)
    data[i] = (data[i] - mean) / sdev;
}
/**************************************-**************************************/
/*                                                                           */
/*                           Z_SCORE_EXACT_FARRAY                            */
/*                                                                           */
/*  This accounts for the "n-1" division in the "sdev" routine.              */
/*                                                                           */
/*****************************************************************************/
void z_score_exact_farray(data,n)
     float *data;
     int n;
{
  int i;
  float mean,sdev,sdev_star;

  mean_sdev_farray(data,n,&mean,&sdev);
  sdev_star = sdev * sqrt((float)(n-1)/(float)n);
  for(i=0;i<n;i++)
    data[i] = (data[i] - mean) / sdev_star;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MEAN_SDEV_FREQ_HIST_FARRAY                      */
/*                                                                           */
/*  Compute the mean and standard deviation of the frequency histogram.      */
/*                                                                           */
/*****************************************************************************/
void mean_sdev_freq_hist_farray(data,n,mean,sdev)
     float *data;
     int n;
     float *mean,*sdev;
{
  int i;
  float m,s,weight,totw;
  
  m = 0.0;
  totw = 0.0;
  for(i=0;i<n;i++){
    weight = data[i];
    m += weight * (float)i;
    totw += weight;
  }
  if (totw > 0)
    m /= totw;
  
  s = 0.0;
  for(i=0;i<n;i++){
    weight = data[i];
    s += weight*(((float)i-m)*((float)i-m));
  }
  if (totw > 1)
    s = sqrt(s/(float)(totw-1.0));
  else
    s = -1.0;
  
  *mean = m;
  *sdev = s;
}
/**************************************-**************************************/
/*                                                                           */
/*                        PARTIAL_CORRELATIONS_FARRAYS                       */
/*                                                                           */
/*  See ~/math/prob for definition of partial correlation.                   */
/*                                                                           */
/*****************************************************************************/
void partial_correlations_farrays(data,data1,data2,n,rp1,rp2,rp1z,rp2z)
     float *data,*data1,*data2;
     int n;
     float *rp1,*rp2;              // r-values
     float *rp1z,*rp2z;            // Fisher-Z transformation values
{
  float r1,r2,r12,p1,p2,p1z,p2z,tiny;

  simple_pearsn(data-1,data1-1,n,&r1);
  simple_pearsn(data-1,data2-1,n,&r2);
  simple_pearsn(data1-1,data2-1,n,&r12);

  // Partial corr. coeff. betw. data and data1, controlling for data2
  p1 = (r1 - r2*r12)/sqrt((1.0-r2*r2)*(1.0-r12*r12));

  // Partial corr. coeff. betw. data and data2, controlling for data1
  p2 = (r2 - r1*r12)/sqrt((1.0-r1*r1)*(1.0-r12*r12));

  // Fisher-Z transformations
  tiny = 1.0e-20;
  //p1z = 0.5 * log((1.0 + p1 + tiny) / (1.0 - p1 + tiny));
  //p2z = 0.5 * log((1.0 + p2 + tiny) / (1.0 - p2 + tiny));

  //
  //  Equation from Smith, Majaj, Movshon (2005), see their page 8
  //
  //  *** WE THINK THEIR EQUATION IS WRONG
  //p1z = 0.5 * log((1.0 + p1 + tiny) / (1.0 - p2 + tiny));
  //p2z = 0.5 * log((1.0 + p2 + tiny) / (1.0 - p1 + tiny));
  //  *** AND HAVE REPLACED IT WITH THIS:
  p1z = 0.5 * log((1.0 + p1 + tiny) / (1.0 - p1 + tiny));
  p2z = 0.5 * log((1.0 + p2 + tiny) / (1.0 - p2 + tiny));
  p1z /= sqrt(1.0/(float)(n-3));
  p2z /= sqrt(1.0/(float)(n-3));


  *rp1 = p1;
  *rp2 = p2;
  *rp1z = p1z;
  *rp2z = p2z;
}
/**************************************-**************************************/
/*                                                                           */
/*                             VECTOR_AVERAGE_FARRAY                         */
/*                                                                           */
/*****************************************************************************/
void vector_average_farray(a,b,n,inflag,outflag,ra,rb)
     float *a,*b;
     int n;
     int inflag;    // 0-cartesian 1-polar_rad 2-polar-deg
     int outflag;   // 0-cartesian 1-polar_rad 2-polar-deg
     float *ra,*rb;
{
  int i;
  float *x,*y,r,theta,xm,xsd,ym,ysd;

  if (inflag == 2){
    x = (float *)myalloc(n*sizeof(float));
    y = (float *)myalloc(n*sizeof(float));
    for(i=0;i<n;i++){
      x[i] = a[i] * cos(b[i]*M_PI/180.0);
      y[i] = a[i] * sin(b[i]*M_PI/180.0);
    }
    mean_sdev_farray(x,n,&xm,&xsd);
    mean_sdev_farray(y,n,&ym,&ysd);
    theta = atan2(ym,xm);
    r = sqrt(xm*xm + ym*ym);

    myfree(x); myfree(y);
  }else{
    r = theta = 0.0;
    exit_error("VECTOR_AVERAGE_FARRAY","inflag not implemented");
  }

  if (outflag == 2){
    *ra = r;
    *rb = theta * 180.0/M_PI;
  }else
    exit_error("VECTOR_AVERAGE_FARRAY","outflag not implemented");

}
/**************************************-**************************************/
/*                                                                           */
/*                   CENTROID_CIRCULAR_VARIANCE_POLAR_FARRAY                 */
/*                                                                           */
/*  Compute the centroid and circular variance.                              */
/*                                                                           */
/*  Note, set 'theta_max' to 180 to match Ringach, Shapley, Hawken (2002)    */
/*  who map 360 deg directions onto a 180 deg orientation circle.   This     */
/*  is presumably necessary so that sharp ori-tuning has a max value.        */
/*                                                                           */
/*  Circ_var is 0 if all points except one point are zero.                   */
/*  Circ_var is 1 is all points have the same non-zero value.                */
/*  Circ_var is "nan" if all values are 0.0.                                 */
/*  Circ_var increases (toward 1.0) if a constant baseline is added.         */
/*  Circ_var does not change if data are multiplied by non-zero value.       */
/*                                                                           */
/*****************************************************************************/
void centroid_circular_variance_polar_farray(amp,theta,n,theta_max,rcentr,
					     rcentt,rcv)
     float *amp;       // amplitude [n]
     float *theta;     // angle 0..theta_max (deg) [n]
     int n;            // number of data points
     float theta_max;  // e.g., 360.0 for dir, 180.0 for ori
     float *rcentr;    // Centroid 'r'
     float *rcentt;    // Centroid 'theta'
     float *rcv;       // Circular variance
{
  int i;
  float tot,x,y,xm,ym,cv,tm2,cr,ctheta,tm2pi;

  tm2pi = 2.0 * M_PI / theta_max;

  tot = x = y = 0.0;
  for(i=0;i<n;i++){
    tot += amp[i];
    x += amp[i] * cos(tm2pi * theta[i]);
    y += amp[i] * sin(tm2pi * theta[i]);
  }
  cv = 1.0 - sqrt(x*x + y*y) / tot;

  xm = x/(float)n;
  ym = y/(float)n;

  cr     = sqrt(xm*xm + ym*ym);    // Centroid length
  ctheta = atan2(ym,xm) / tm2pi;   // direction (-theta_max/2..theta_max/2)
  if (ctheta < 0.0)
    ctheta += theta_max;           // Centroid direction 0..theta_max

  *rcv = cv;
  *rcentr = cr;
  *rcentt = ctheta;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 MEAN_2D_FARRAY                            */
/*                                                                           */
/*  Compute the row that is the mean of all rows (i.e., the ith entry in     */
/*  the new row is the mean of the entries in the ith column.)               */
/*  compute the standard deviations for that row.                            */
/*                                                                           */
/*****************************************************************************/
void mean_2d_farray(data,m,n,rmean)
     float **data;
     int m,n;
     float **rmean;
{
  int i,j;
  float *mu;
  
  if (m > 1){
    mu = get_zero_farray(n);
    for(i=0;i<n;i++){
      for(j=0;j<m;j++)
	mu[i] += data[j][i];
      mu[i] /= (float)m;
    }
  }else if (m == 1)
    mu = copy_farray(data[0],n);
  else{
    mu = NULL;
    exit_error("MEAN_2D_FARRAY","Array has zero rows");
  }

  *rmean = mu;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MEAN_3D_FARRAY                              */
/*                                                                           */
/*  Compute the 2d array that is the mean of all 2d arrays                   */
/*                                                                           */
/*****************************************************************************/
void mean_3d_farray(data,n1,n2,n3,rmean)
     float ***data;
     int n1,n2,n3;
     float ***rmean;
{
  int i,j,k;
  float **mu;
  
  if (n1 > 1){
    mu = get_zero_2d_farray(n2,n3);
    for(i=0;i<n2;i++)
      for(j=0;j<n3;j++){
	for(k=0;k<n1;k++)
	  mu[i][j] += data[k][i][j];
	mu[i][j] /= (float)n1;
      }
  }else if (n1 == 1)
    mu = copy_2d_farray(data[0],n2,n3);
  else{
    mu = NULL;
    exit_error("MEAN_3D_FARRAY","Array has zero rows");
  }

  *rmean = mu;
}
/**************************************-**************************************/
/*                                                                           */
/*                          SINGLE_MEAN_SDEV_2D_FARRAY                       */
/*                                                                           */
/*  Compute the mean and sdev of all elements contained in the 2d array.     */
/*                                                                           */
/*  *** WARNING:  divides by 'n*m', not 'n*m-1'                              */
/*                                                                           */
/*****************************************************************************/
void single_mean_sdev_2d_farray(data,m,n,rm,rsd)
     float **data;
     int m,n;
     float *rm,*rsd;
{
  int i,j;
  float mu,sd;

  mu = 0.0;
  for(i=0;i<m;i++)
    for(j=0;j<n;j++)
      mu += data[i][j];
  mu /= (float)(m*n);

  sd = 0.0;
  for(i=0;i<m;i++)
    for(j=0;j<n;j++)
      sd += ((data[i][j] - mu)*(data[i][j] - mu));

  if ((m*n) > 1)
    sd = sqrt(sd/(float)(m*n));
  else
    sd = 0.0;

  *rm = mu; *rsd = sd;
}
/**************************************-**************************************/
/*                                                                           */
/*                           SINGLE_MEAN_3D_FARRAY                           */
/*                                                                           */
/*  Compute the mean of the 3D farray.                                       */
/*                                                                           */
/*****************************************************************************/
void single_mean_3d_farray(data,x0,xn,y0,yn,z0,zn,rmu)
     float ***data;
     int x0,xn,y0,yn,z0,zn;
     float *rmu;
{
  int i,j,k;
  int cnt;
  float mu;

  //printf("%d %d %d %d %d %d\n",x0,xn,y0,yn,z0,zn);

  cnt = 0;
  mu = 0.0;
  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++)
      for(k=0;k<zn;k++){
	mu += data[x0+i][y0+j][z0+k];
	cnt += 1;
      }
  if (cnt > 0)
    mu /= (float)cnt;
  else
    exit_error("SINGLE_MEAN_3D_FARRAY","No values");

  *rmu = mu;
}
/**************************************-**************************************/
/*                                                                           */
/*                         SINGLE_MEAN_SDEV_3D_FARRAY                        */
/*                                                                           */
/*  Compute the mean and SD of the 3D farray.                                */
/*                                                                           */
/*****************************************************************************/
void single_mean_sdev_3d_farray(data,x0,xn,y0,yn,z0,zn,rmu,rsd)
     float ***data;
     int x0,xn,y0,yn,z0,zn;
     float *rmu,*rsd;
{
  int i,j,k;
  int cnt;
  float mu,sd;

  //printf("%d %d %d %d %d %d\n",x0,xn,y0,yn,z0,zn);

  cnt = 0;
  mu = 0.0;
  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++)
      for(k=0;k<zn;k++){
	mu += data[x0+i][y0+j][z0+k];
	cnt += 1;
      }
  if (cnt > 0)
    mu /= (float)cnt;
  else
    exit_error("SINGLE_MEAN_SDEV_3D_FARRAY","No values");


  sd = 0.0;
  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++)
      for(k=0;k<zn;k++){
	sd += ((data[i][j][k] - mu)*(data[i][j][k] - mu));
      }

  if (cnt > 1)
    sd = sqrt(sd/(float)cnt);
  else
    sd = 0.0;

  *rmu = mu; *rsd = sd;
}
/**************************************-**************************************/
/*                                                                           */
/*                           SUBTRACT_MEAN_2D_FARRAY                         */
/*                                                                           */
/*  Take an array of floats, subtract the mean.  This used to be called      */
/*  "make_mean_zero_farray".                                                 */
/*                                                                           */
/*****************************************************************************/
void subtract_mean_2d_farray(data,n,m,oldmean)
     float **data;
     int n,m;
     float *oldmean;
{
  int i,j;
  float mean,sdev;

  single_mean_sdev_2d_farray(data,n,m,&mean,&sdev);

  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      data[i][j] -= mean;
  *oldmean = mean;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MEAN_SDEV_2D_FARRAY                           */
/*                                                                           */
/*  Compute the row that is the mean of all rows (i.e., the ith entry in     */
/*  the new row is the mean of the entries in the ith column.)  Also,        */
/*  compute the standard deviations for that row.                            */
/*                                                                           */
/*  Returned arrays have 'n' entries.                                        */
/*                                                                           */
/*****************************************************************************/
void mean_sdev_2d_farray(data,m,n,mean,sdev)
     float **data;
     int m,n;
     float **mean,**sdev;
{
  int i,j;
  float *mu,*sigma,*t;
  
  sigma = get_zero_farray(n);
  if (m > 1){
    mu = get_zero_farray(n);
    t = (float *)myalloc(m*sizeof(float));
    for(i=0;i<n;i++){
      for(j=0;j<m;j++)
	t[j] = data[j][i];
      mean_sdev_farray(t,m,&(mu[i]),&(sigma[i]));
    }
    myfree(t);
  }else if (m == 1)
    mu = copy_farray(data[0],n);
  else{
    mu = NULL;
    exit_error("MEAN_SDEV_2D_FARRAY","Array has zero rows");
  }

  *mean = mu; *sdev = sigma;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MEAN_SD_TRIAL_2D_FARRAY                         */
/*                                                                           */
/*  Think of each 'data[i]' as a trial which has a mean value.  Compute      */
/*  the MEAN and SD of these mean values across trials.                      */
/*                                                                           */
/*****************************************************************************/
void mean_sd_trial_2d_farray(data,n,m,t0,tn,rmu,rsd)
     float **data;        // [n][m]
     int n,m;             // 
     int t0,tn;           // Start averaging at 't0', go for 'tn'
     float *rmu,*rsd;
{
  int i;
  float *tm,tsd;

  if ((t0 < 0) || ((t0+tn) > m)){
    printf("  Data length is %d units\n",m);
    printf("  t0 %d   tn %d\n",t0,tn);
    exit_error("MEAN_SD_TRIAL_2D_FARRAY","Indices out of bounds");
  }

  tm = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++){
    mean_sdev_farray(&(data[i][t0]),tn,&(tm[i]),&tsd);
  }
  mean_sdev_farray(tm,n,rmu,rsd);
  myfree(tm);
}
/**************************************-**************************************/
/*                                                                           */
/*                            AVG_1ST_DIM_2D_FARRAY                          */
/*                                                                           */
/*  Average over the 1st array dimension, and return the array that is the   */
/*  average.                                                                 */
/*                                                                           */
/*  Returned array has 'n' entries.                                          */
/*                                                                           */
/*****************************************************************************/
float *avg_1st_dim_2d_farray(data,m,n)
     float **data;  // [m][n]
     int m,n;
{
  int i,j;
  float sum,*mu;

  mu = get_farray(n);

  for(j=0;j<n;j++){
    sum = 0.0;
    for(i=0;i<m;i++)
      sum += data[i][j];
    mu[j] = sum/(float)m;
  }

  return mu;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MEAN_SDEV_COL_2D_FARRAY                         */
/*                                                                           */
/*  Compute the col that is the mean of all cols (i.e., the ith entry in     */
/*  the new col is the mean of the entries in the ith row.)  Also, compute   */
/*  the standard deviations.                                                 */
/*                                                                           */
/*  Returned arrays have 'm' entries.                                        */
/*                                                                           */
/*****************************************************************************/
void mean_sdev_col_2d_farray(data,m,n,mean,sdev)
     float **data;
     int m,n;
     float **mean,**sdev;
{
  int i;
  float *mu,*sigma;
  
  sigma = get_zero_farray(m);
  mu = get_zero_farray(m);
  if (n > 1){
    for(i=0;i<m;i++)
      mean_sdev_farray(data[i],n,&(mu[i]),&(sigma[i]));
  }else if (n == 1)
    for(i=0;i<m;i++){
      mu[i] = data[i][0];
      sigma[i] = 0.0;
    }
  else{
    mu = NULL;
    exit_error("MEAN_SDEV_COL_2D_FARRAY","Array has zero cols");
  }

  *mean = mu; *sdev = sigma;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MEAN_SDEV_2D_FLAG_FARRAY                         */
/*                                                                           */
/*  Compute the mean and sdev for only those rows with flag[i] = 1.          */
/*                                                                           */
/*****************************************************************************/
void mean_sdev_2d_flag_farray(data,m,n,flag,mean,sdev)
     float **data;
     int m,n,*flag;
     float **mean,**sdev;
{
  int i,j;
  int count;
  float **t;

  count = 0;
  for(i=0;i<n;i++)
    if (flag[i] == 1)
      count += 1;

  t = (float **)myalloc(count*sizeof(float *));
  j = 0;
  for(i=0;i<n;i++)
    if (flag[i] == 1){
      t[j] = data[i];
      j += 1;
    }
  mean_sdev_2d_farray(t,count,n,mean,sdev);
  myfree(t);
  
}  
/**************************************-**************************************/
/*                                                                           */
/*                         MEAN_SDEV_2D_SEGMENT_FARRAY                       */
/*                                                                           */
/*  Compute the row that is the mean of all rows (i.e., the ith entry in     */
/*  the new row is the mean of the entries in the ith column.)  Also,        */
/*  compute the standard deviations for that row.                            */
/*                                                                           */
/*****************************************************************************/
void mean_sdev_2d_segment_farray(data,m,t0,tn,mean,sdev)
     float **data;
     int m,t0,tn;
     float **mean,**sdev;
{
  int i,j;
  int tfin;
  float *mu,*sigma,*t;

  tfin = t0+tn;
  sigma = get_zero_farray(tn);
  if (m > 1){
    mu = get_zero_farray(tn);
    t = (float *)myalloc(m*sizeof(float));
    for(i=t0;i<tfin;i++){
      for(j=0;j<m;j++)
	t[j] = data[j][i];
      mean_sdev_farray(t,m,&(mu[i-t0]),&(sigma[i-t0]));
    }
    myfree(t);
  }else if (m == 1)
    mu = copy_farray(data[0]+t0,tn);
  else{
    mu = NULL;
    exit_error("MEAN_SDEV_2D_SEGMENT_FARRAY","Array has zero rows");
  }

  *mean = mu; *sdev = sigma;
}
/**************************************-**************************************/
/*                                                                           */
/*                        WEIGHTED_MEAN_SDEV_2D_FARRAY                       */
/*                                                                           */
/*  Like "mean_sdev_2d_farray", but each row has an associated "weight".     */
/*                                                                           */
/*****************************************************************************/
void weighted_mean_sdev_2d_farray(data,m,n,weight,mean,sdev)
     float **data;
     int m,n;
     float *weight,**mean,**sdev;
{
  int i,j;
  float *mu,*sigma,*t;
  
  mu = get_zero_farray(n);
  sigma = get_zero_farray(n);
  t = (float *)myalloc(m*sizeof(float));
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      t[j] = data[j][i];
    weighted_mean_sdev_farray(t,m,weight,&(mu[i]),&(sigma[i]));
  }
  myfree(t);
  *mean = mu; *sdev = sigma;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MEAN_SDEV_DIFF_2D_FARRAY                        */
/*                                                                           */
/*  Compute the mean of the difference between corresponding rows in the     */
/*  2D farrays.                                                              */
/*                                                                           */
/*****************************************************************************/
void mean_sdev_diff_2d_farray(data1,data2,m,n,rmean,rsdev)
     float **data1,**data2;
     int m,n;
     float **rmean,**rsdev;
{
  int i;
  float **diff;

  diff = (float **)myalloc(m*sizeof(float *));
  for(i=0;i<m;i++)
    diff[i] = diff_farrays(data1[i],data2[i],n);
  mean_sdev_2d_farray(diff,m,n,rmean,rsdev);

  free_2d_farray(diff,m);
}
/**************************************-**************************************/
/*                                                                           */
/*                      WEIGHTED_MEAN_SDEV_DIFF_2D_FARRAY                    */
/*                                                                           */
/*  Compute the weighted mean of the difference between corresponding rows   */
/*  in the 2D farrays.                                                       */
/*                                                                           */
/*****************************************************************************/
void weighted_mean_sdev_diff_2d_farray(data1,data2,m,n,weight,rmean,rsdev)
     float **data1,**data2;
     int m,n;
     float *weight,**rmean,**rsdev;
{
  int i;
  float **diff;

  diff = (float **)myalloc(m*sizeof(float *));
  for(i=0;i<m;i++)
    diff[i] = diff_farrays(data1[i],data2[i],n);
  weighted_mean_sdev_2d_farray(diff,m,n,weight,rmean,rsdev);

  free_2d_farray(diff,m);
}
/**************************************-**************************************/
/*                                                                           */
/*                           MAKE_MAX_CONST_2D_FARRAY                        */
/*                                                                           */
/*  Rescale the farray so that its maximum value is "c".                     */
/*                                                                           */
/*****************************************************************************/
void make_max_const_2d_farray(data,m,n,c)
     float **data;
     int m,n;
     float c;
{
  float max;
  
  max = max_of_2d_farray(data,m,n);
  if (max == 0.0)
    exit_error("MAKE_MAX_CONST_2D_FARRAY","Max is zero");
  multiply_2d_farray(data,m,n,c/max);
}
/**************************************-**************************************/
/*                                                                           */
/*                             NORM_01_2D_FARRAY                             */
/*                                                                           */
/*  Rescale the farray to map into [0..1].                                   */
/*                                                                           */
/*****************************************************************************/
void norm_01_2d_farray(data,m,n,flag)
     float **data;
     int m,n;
     int flag;  // 1-force the value '0' to map to '0.5', thus 0 is gray
{
  int i,j;
  float dmin,dmax,drange,d,ss;
  
  get_min_max_2d_farray(data,m,n,&dmin,&dmax);

  if (flag == 0){
    drange = dmax - dmin;
    
    for(i=0;i<m;i++){
      for(j=0;j<n;j++){
	d = data[i][j];
	data[i][j] = (d - dmin)/drange;
      }
    }
  }else if (flag == 1){

    if (fabs(dmax) > fabs(dmin))
      ss = 2.0 * fabs(dmax);
    else
      ss = 2.0 * fabs(dmin);

    for(i=0;i<m;i++){
      for(j=0;j<n;j++){
	d = data[i][j] / ss;   // Now in range -0.5 to 0.5
	data[i][j] = 0.5 + d;  // Now 0..1, with 0 --> 0.5
      }
    }
  }else if (flag == -1){
    //
    //  Truncate, but do not scale
    //
    for(i=0;i<m;i++){
      for(j=0;j<n;j++){
	if (data[i][j] < 0.0){
	  printf("  data[i][j] = %f (< 0.0)\n",data[i][j]);
	  data[i][j] = 0.0;
	}else if (data[i][j] > 1.0){
	  printf("  data[i][j] = %f (> 1.0)\n",data[i][j]);
	  data[i][j] = 1.0;
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             NORM_01_3D_FARRAY                             */
/*                                                                           */
/*  Rescale the farray to map into [0..1].                                   */
/*                                                                           */
/*****************************************************************************/
void norm_01_3d_farray(data,xn,yn,zn,flag)
     float ***data;
     //int m,n;
     int xn,yn,zn;
     int flag;  // 1-force the value '0' to map to '0.5', thus 0 is gray
{
  int i,j,k;
  float dmin,dmax,drange,d,ss;
  
  get_min_max_3d_farray(data,0,xn,0,yn,0,zn,&dmin,&dmax);

  if (flag == 0){
    drange = dmax - dmin;
    
    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){
	for(k=0;k<zn;k++){
	  d = data[i][j][k];
	  data[i][j][k] = (d - dmin)/drange;
	}
      }
    }
  }else if (flag == 1){

    if (fabs(dmax) > fabs(dmin))
      ss = 2.0 * fabs(dmax);
    else
      ss = 2.0 * fabs(dmin);

    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){
	for(k=0;k<zn;k++){
	  d = data[i][j][k] / ss;   // Now in range -0.5 to 0.5
	  data[i][j][k] = 0.5 + d;  // Now 0..1, with 0 --> 0.5
	}
      }
    }
  }else if (flag == -1){
    //
    //  Truncate, but do not scale
    //
    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){
	for(k=0;k<zn;k++){
	  if (data[i][j][k] < 0.0){
	    printf("  data[i][j][k] = %f (< 0.0)\n",data[i][j][k]);
	    data[i][j][k] = 0.0;
	  }else if (data[i][j][k] > 1.0){
	    printf("  data[i][j][k] = %f (> 1.0)\n",data[i][j][k]);
	    data[i][j][k] = 1.0;
	  }
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           RECT_TO_POLAR_2D_FARRAY                         */
/*                                                                           */
/*  nt - number of samples around circle                                     */
/*  nr - number of samples along raduis                                      */
/*                                                                           */
/*  Assume that the circle circumscribes the square.                         */
/*                                                                           */
/*****************************************************************************/
float **rect_to_polar_2d_farray(data,m,n,nt,nr)
     float **data;
     int m,n,nt,nr;
{
  int i,j;
  float **pol,theta,r,x,y,rad,xc,yc;

  if (m!=n)
    exit_error("RECT_TO_POLAR_2D_FARRAY","Array must be square");

  rad = sqrt(n*n + n*n)/2.0;
  xc = yc = (float)n/2.0;
  pol = get_2d_farray(nt,nr);

  /* Assume that 0,0 is the upper left corner of the rectangular array. */
  /* The first radius is taken from the center of the array downward.   */
  /* Subsequent radii are taken CCW from the center outward. */
  for(i=0;i<nt;i++){
    theta = (float)i/(float)nt * 2.0*M_PI;
    for(j=0;j<nr;j++){
      r = (float)j/(float)(nr-1) * rad;
      x = xc + r*cos(theta);
      y = yc + r*sin(theta);
      /*printf("x,y = %.2f %.2f\n",x,y);*/
      if ((x < 0.0)||(y < 0.0)||(x > (float)(m-1))||(y > (float)(n-1))){
	pol[i][j] = pol[i][j-1]; /* Use closer in value on radius */
      }else{
	pol[i][j] = interpolate_2d_farray(data,m,n,x,y);
      }
    }
  }
  return pol;
}
/**************************************-**************************************/
/*                                                                           */
/*                           POLAR_TO_RECT_2D_FARRAY                         */
/*                                                                           */
/*  nt - number of samples around circle                                     */
/*  nr - number of samples along raduis                                      */
/*                                                                           */
/*  Assume that the square just fits within the circle.                      */
/*                                                                           */
/*****************************************************************************/
float **polar_to_rect_2d_farray(polar,nt,nr,m,n)
     float **polar;
     int nt,nr,m,n;
{
  int i,j;
  float **data,t,tt,r,rr,x,y,rad,xc,yc;

  if (m!=n)
    exit_error("POLAR_TO_RECT_2D_FARRAY","Array must be square");

  rad = sqrt(n*n + n*n)/2.0;
  xc = yc = (float)n/2.0;
  data = get_2d_farray(m,n);

  for(i=0;i<m;i++){ /* Vertical on image */
    x = (float)(i-xc);
    for(j=0;j<n;j++){ /* Horizontal on image */
      y = (float)(j-yc);
      r = sqrt(y*y + x*x);
      rr = r/rad * (float)(nr-1);
      t = atan2(y,x);
      if (t < 0.0)
	t += 2.0*M_PI;
      /*printf("xyt = %.2f %.2f %.4f  %.4f\n",x,y,t,atan2(0.0,0.0));*/
      tt = t/(M_PI*2.0) * (float)nt;

      /*printf("%3d %3d  %d %d\n",i,j,(int)tt,(int)rr);*/

      if ((tt < 0.0)||(rr < 0.0)||(tt > (float)(nt-1))||(rr > (float)(nr-1))){
	printf("%3d %3d  %.2f %.2f\n",i,j,tt,rr);
	if (tt > (float)(nt-1))
	  tt = (float)(nt-1);
	else
	  exit_error("POLAR_TO_RECT_2D_FARRAY","Out of bounds");
      }
      
      data[i][j] = interpolate_2d_farray(polar,nt,nr,tt,rr);
    }
  }
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                              NORM_AREA_FARRAY                             */
/*                                                                           */
/*  Normalize the area of the farray.                                        */
/*                                                                           */
/*****************************************************************************/
void norm_area_farray(data,n,norm)
     float *data;
     int n;
     float norm;
{
  int i;
  float area,c;

  area = sum_farray(data,n,0,n);
  if (area == 0.0)
    exit_error("NORM_AREA_FARRAY","Area = 0.0");
  c = norm/area;
  for(i=0;i<n;i++)
    data[i] *= c;
}
/**************************************-**************************************/
/*                                                                           */
/*                            NORM_AREA_2D_FARRAY                            */
/*                                                                           */
/*  Normalize the area of the farray.                                        */
/*                                                                           */
/*****************************************************************************/
void norm_area_2d_farray(data,m,n,norm)
     float **data;
     int m,n;
     float norm;
{
  int i,j;
  double sum,c;

  sum = 0.0;
  for(i=0;i<m;i++)
    for(j=0;j<n;j++)
      sum += data[i][j];

  if (sum == 0.0)
    exit_error("NORM_AREA_2D_FARRAY","Area = 0.0");

  c = norm/sum;
  for(i=0;i<m;i++)
    for(j=0;j<n;j++)
      data[i][j] *= c;
}
/**************************************-**************************************/
/*                                                                           */
/*                          NORM_POSITIVE_AREA_FARRAY                        */
/*                                                                           */
/*  Normalize the positive area of the farray.                               */
/*                                                                           */
/*****************************************************************************/
void norm_positive_area_farray(data,n,norm)
     float *data;
     int n;
     float norm;
{
  int i;
  float area,c;

  area = 0.0; /*** Compute positive area ***/
  for(i=0;i<n;i++)
    if (data[i] > 0.0)
      area += data[i];

  if (area == 0.0)
    exit_error("NORM_AREA_FARRAY","Area = 0.0");
  c = norm/area;
  for(i=0;i<n;i++)
    data[i] *= c;
}
/**************************************-**************************************/
/*                                                                           */
/*                             NORM_MAX_FARRAY                               */
/*                                                                           */
/*  Normalize the farray so that the maximum value is "norm" by              */
/*  multiplying by norm/max where "max" is the maximum value in              */
/*  the array.  If the maximum value is 0.0, exit.                           */
/*                                                                           */
/*****************************************************************************/
void norm_max_farray(data,n,norm)
     float *data;
     int n;
     float norm;
{
  int i;
  float max,c;

  max = max_of_farray(data,n);
  if (max == 0.0)
    exit_error("NORM_MAX_FARRAY","Max = 0.0");
  c = norm/max;
  for(i=0;i<n;i++)
    data[i] *= c;
}
/**************************************-**************************************/
/*                                                                           */
/*                            NORM_VARIANCE_FARRAY                           */
/*                                                                           */
/*  Modify "data" to have variance 1.                                        */
/*                                                                           */
/*****************************************************************************/
void norm_variance_farray(data,n)
     float *data;
     int n;
{
  float mean,sdev;

  mean_sdev_farray(data,n,&mean,&sdev);
  if (sdev == 0.0)
    printf("*** NORM_VARIANCE_FARRAY  Warning: sdev=0, cannot normalize.\n");
  else
    multiply_farray(data,n,1.0/sdev);
}
/**************************************-**************************************/
/*                                                                           */
/*                         NORM_VARIANCE_QUIET_FARRAY                        */
/*                                                                           */
/*  Modify "data" to have variance 1.                                        */
/*                                                                           */
/*****************************************************************************/
void norm_variance_quiet_farray(data,n)
     float *data;
     int n;
{
  float mean,sdev;

  mean_sdev_farray(data,n,&mean,&sdev);
  if (sdev > 0.0)
    multiply_farray(data,n,1.0/sdev);
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_CUMULATIVE_FARRAY                          */
/*                                                                           */
/*****************************************************************************/
float *get_cumulative_farray(data,n)
     float *data;
     int n;
{
  int i;
  float *c;

  c = (float *)myalloc(n*sizeof(float *));
  c[0] = data[0];
  for(i=1;i<n;i++)
    c[i] = c[i-1] + data[i];
  
  return c;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_FOURIER_HARMONIC_FARRAY                       */
/*                                                                           */
/*  Ampl and theta are returned for the "order"th component.  Theta is in    */
/*  degrees from -180 to 180.                                                */
/*                                                                           */
/*****************************************************************************/
void get_fourier_harmonic_farray(data,n,order,period,rampl,rtheta)
     float *data;
     int n,order;
     float period,*rampl,*rtheta;
{
  int i;
  float x,todd,teven;

  todd = teven = 0.0;
  for(i=0;i<n;i++){
    x = (float)i/(period/(float)order) * 2.0*M_PI;
    todd += 2.0 * sin(x) * data[i];
    teven += 2.0 * cos(x) * data[i];
  }
  teven /= (float)n;
  todd /= (float)n;
  *rampl = sqrt(todd*todd + teven*teven);

  *rtheta = (float)atan2(todd,teven) * 180.0/M_PI;
}
/**************************************-**************************************/
/*                                                                           */
/*                     GET_FOURIER_HARMONIC_TRIAL_2D_FARRAY                  */
/*                                                                           */
/*  Think of each 'data[i]' as a trial which has an AMP and PHASE.  Compute  */
/*  the MEAN and SD of these values across trials.                           */
/*                                                                           */
/*****************************************************************************/
void get_fourier_harmonic_trial_2d_farray(data,n,m,t0,tn,period,sampling,
					  order,ramplmean,ramplsd,
					  rphmean,rphsd,rampl)
     float **data;       // [n][m]
     int n,m;            //
     int t0,tn;          // start index, length of data segment (array units)
     float period;       // 1/frequency, (array units)
     float sampling;     // samples per second
     int order;          // 1-fundamental
     float *ramplmean,*ramplsd;  // amplitude mean, SD
     float *rphmean,*rphsd;      // phase mean, SD
     float **rampl;              // amplitude for each trial [n]
{
  int i;
  float *ampl,*phase,*td,aperiod,r,theta,*d;

  if (period <= 0.0){  // Return all zeros
    printf("  *** (get_fourier_harmonic_trial_2d_farray) %s\n",
	   "WARNING:  period is <= 0.0");
    *ramplmean = 0.0;
    *ramplsd = 0.0;
    *rphmean = 0.0;
    *rphsd = 0.0;
    *rampl = get_zero_farray(n);
    return;
  }

  /*
  printf("__________________\n");
  printf("  n,m  %d,%d\n",n,m);
  printf("t0,tn %d,%d\n",t0,tn);
  printf("period %f\n",period);
  printf("sampling %f\n",sampling);
  printf("------------------\n");
  */

  //aperiod = period*sampling;  // cycle period (array units)
  aperiod = period;  // cycle period (array units)

  //printf("  aperiod = %f\n",aperiod);

  ampl = (float *)myalloc(n*sizeof(float));
  phase = (float *)myalloc(n*sizeof(float));

  for(i=0;i<n;i++){ // For each data segment
    td = &(data[i][t0]);
    get_fourier_harmonic_farray(td,tn,order,aperiod,&(ampl[i]),&(phase[i]));
    //printf("ampl,phase   %f   %f\n",ampl[i],phase[i]);
  }
  mean_sdev_farray(ampl,n,ramplmean,ramplsd);

  // Use vector average to find rough direction of mean
  vector_average_farray(ampl,phase,n,2,2,&r,&theta);

  // Compute mean around the direction of the vector average
  d = get_farray(n);
  for(i=0;i<n;i++){
    d[i] = get_signed_circular_diff(phase[i]+180.0,theta+180.0,360.0);
  }
  mean_sdev_farray(d,n,rphmean,rphsd);
  myfree(d);

  *rphmean += theta;

  *rampl = ampl; // Return individual ampl data points for each trial
  myfree(phase);
}
/**************************************-**************************************/
/*                                                                           */
/*                                SORT_FARRAY                                */
/*                                                                           */
/*  Call the NumRecInC sort (quicksort) program.                             */
/*                                                                           */
/*****************************************************************************/
void sort_farray(data,n)
     float *data;
     int n;
{
  sort(n,data-1);
}
/**************************************-**************************************/
/*                                                                           */
/*                           INDEX_BUBBLE_SORT_FARRAY                        */
/*                                                                           */
/*  For longer arrays, use "INDEX_SORT_FARRAY" below, it might be quicker.   */
/*                                                                           */
/*****************************************************************************/
int *index_bubble_sort_farray(data,n)
     float *data;
     int    n;
{
  int i;
  int *index,t;
  int done,bottom;
  
  index = (int *)myalloc(n*sizeof(int));
  for (i=0;i<n;i++)
    index[i] = i;

  bottom = n-1;
  done = 0;
  while (!done){
    done = 1;
    for (i=0;i<bottom;i++)
      if (data[index[i]] > data[index[i+1]]){
	done = 0;
	t = index[i];
	index[i] = index[i+1];
	index[i+1] = t;
      }
    bottom -= 1;
  }
  return index;
}
/**************************************-**************************************/
/*                                                                           */
/*                             FARRAY_SORT_CARRAY                            */
/*                                                                           */
/*****************************************************************************/
void farray_sort_carray(data,n)
     char **data;  // Sort this based on its 'atof' values
     int    n;
{
  int i;
  int *index;
  float *fd;
  char **fc;

  fd = get_farray(n);
  fc = (char **)myalloc(n*sizeof(char *));
  for(i=0;i<n;i++){
    fd[i] = atof(data[i]);  // float version of 'data'
    fc[i] = data[i];        // copy of 'data'
  }

  index = index_bubble_sort_farray(fd,n);

  for(i=0;i<n;i++){
    data[i] = fc[index[i]];
  }
  
  myfree(fd);
  myfree(fc);
  myfree(index);
}
/*************************************---*************************************/
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp
#define M 7
#define NSTACK 50
/**************************************-**************************************/
/*                                                                           */
/*                              SORT_THREE_FARRAY                            */
/*                                                                           */
/*  Quicksort based on the first farray, rearrange the other two.            */
/*                                                                           */
/*  This is a modified version of "sort2" from:                              */
/*                                                                           */
/*C  (C) Copr. 1986-92 Numerical Recipes Software #.3.                       */
/*                                                                           */
/*****************************************************************************/
void sort_three_farray(wd1,wd2,wd3,n)
     float *wd1,*wd2,*wd3;
     int n;
{
  unsigned long i,ir=n,j,k,l=1;
  int *istack,jstack=0;
  float a,b,c,temp;
  float *d1,*d2,*d3;

  d1 = wd1-1;
  d2 = wd2-1;
  d3 = wd3-1;
  
  istack = (int *)myalloc((NSTACK+1)*sizeof(int));
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	a=d1[j];
	b=d2[j];
	c=d3[j];
	for (i=j-1;i>=1;i--) {
	  if (d1[i] <= a) break;
	  d1[i+1]=d1[i];
	  d2[i+1]=d2[i];
	  d3[i+1]=d3[i];
	}
	d1[i+1]=a;
	d2[i+1]=b;
	d3[i+1]=c;
      }
      if (!jstack) {
	myfree(istack);
	return;
      }
      ir=istack[jstack];
      l=istack[jstack-1];
      jstack -= 2;
    }else{
      k=(l+ir) >> 1;
      SWAP(d1[k],d1[l+1]);
      SWAP(d2[k],d2[l+1]);
      SWAP(d3[k],d3[l+1]);
      if (d1[l+1] > d1[ir]) {
	SWAP(d1[l+1],d1[ir]);
	SWAP(d2[l+1],d2[ir]);
	SWAP(d3[l+1],d3[ir]);
      }
      if (d1[l] > d1[ir]) {
	SWAP(d1[l],d1[ir]);
	SWAP(d2[l],d2[ir]);
	SWAP(d3[l],d3[ir]);
      }
      if (d1[l+1] > d1[l]) {
	SWAP(d1[l+1],d1[l]);
	SWAP(d2[l+1],d2[l]);
	SWAP(d3[l+1],d3[l]);
      }
      i=l+1;
      j=ir;
      a=d1[l];
      b=d2[l];
      c=d3[l];
      for (;;) {
	do i++; while (d1[i] < a);
	do j--; while (d1[j] > a);
	if (j < i) break;
	SWAP(d1[i],d1[j]);
	SWAP(d2[i],d2[j]);
	SWAP(d3[i],d3[j]);
      }
      d1[l]=d1[j];
      d1[j]=a;
      d2[l]=d2[j];
      d2[j]=b;
      d3[l]=d3[j];
      d3[j]=c;
      jstack += 2;
      if (jstack > NSTACK)
	exit_error("SORT_THREE_FARRAY","NSTACK too small");
      if (ir-i+1 >= j-l){
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      }else{
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
}
#undef M
#undef NSTACK
#undef SWAP
/**************************************-**************************************/
/*                                                                           */
/*                             INDEX_SORT_FARRAY                             */
/*                                                                           */
/*  Do not alter the input array.  Return a list of indices for the sort.    */
/*                                                                           */
/*****************************************************************************/
int *index_sort_farray(data,n)
     float *data;
     int    n;
{
  int i;
  int *index;
  float *findex,*dcopy,*dummy;

  dcopy = copy_farray(data,n);
  dummy = get_zero_farray(n);

  findex = (float *)myalloc(n*sizeof(float));
  for (i=0;i<n;i++)
    findex[i] = (float)i;

  sort_three_farray(dcopy,findex,dummy,n);

  // Wyeth - chose not to use "f2iarray_round()" here, to avoid dependency
  //         on 'iarray_util'
  index = (int *)myalloc(n*sizeof(int));
  for (i=0;i<n;i++)
    index[i] = my_rint(findex[i]);

  myfree(dcopy);
  myfree(dummy);
  myfree(findex);

  return index;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 GAUSSIAN                                  */
/*                                                                           */
/*****************************************************************************/
float gaussian(sigma,x)
     double sigma,x;
{
  return(1.0/(sqrt(2.0*M_PI)*sigma)*exp(-0.5*(x/sigma)*(x/sigma)));
}
/**************************************-**************************************/
/*                                                                           */
/*                              DISCRETE_GAUSSIAN                            */
/*                                                                           */
/*  This procedure generates an array of values describing a discrete        */
/*  gaussian of standard deviation "sigma" such that the height of the       */
/*  sample furthest from the mean, zero, is a smaller fraction than "trunc"  */
/*  of the height at the mean.                                               */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - "width" returned width will be an odd number.                          */
/*                                                                           */
/*****************************************************************************/
void discrete_gaussian(sigma,trunc,dest,width)
     float sigma,trunc,**dest;
     int   *width;
{
  int i,j;
  float *buffer,peak,area;

  i = 0;
  if (sigma > 0.0){
    peak = gaussian(sigma,(float)0);
    while ((gaussian(sigma,(float)i)/peak) > trunc)
      i += 1;
  }
  *width = 2*i+1;
  
  buffer = (float *)myalloc(*width*sizeof(float));

  if (sigma > 0.0)
    for (j=0;j<*width;j++)
      buffer[j] = gaussian(sigma,(float)(j-i));
  else
    buffer[0] = 1.0;

  // normalize mask area to 1
  area = 0.0;
  for(i=0;i<*width;i++)
    area += buffer[i];
  for(i=0;i<*width;i++)
    buffer[i] /= area;
  
  *dest = buffer;
}
/**************************************-**************************************/
/*                                                                           */
/*                              FIT_TO_GAUSSIAN                              */
/*                                                                           */
/*  Use the NumRec Levenberg-Marquardt routine to fit the data to a          */
/*  guassian using the given initial parameters for A, M, S.                 */
/*                                                                           */
/*                      (x-M)^2                                              */
/*       f(x) = A exp - -------                                              */
/*                        S^2                                                */
/*                                                                           */
/*  The flags indicate by non-zero entries which of a,m,s are to be fitted.  */
/*  This routine sends "fgauss", the NumRec sum of Gaussians routine, to     */
/*  the fitter.                                                              */
/*                                                                           */
/*****************************************************************************/
void fit_to_gaussian(x,y,sd,n,a,m,s,aflag,mflag,sflag)
     float *x,*y,*sd; /* Data with standard deviations. */
     int n;                      /* Number of data points. */
     float *a,*m,*s;             /* Initial guesses for parameters */
     int aflag,mflag,sflag;
{
  int i;
  int done,itcount,ma,*ia,pflag;
  float diff,old_chisq,old_alamda;
  float alamda,chisq,*coeff,**covar,**alpha;

  pflag = 0;

  if (pflag)
    printf("  FIT_TO_GAUSSIAN\n");
  
  for(i=0;i<n;i++)
    if (sd[i] <= 0.0)
      exit_error("FIT_TO_GAUSSIAN","Standard deviation <= 0.0");

  ma = 3; // PREPARATION FOR "mrqmin"
  covar = matrix(1,ma,1,ma); alpha = matrix(1,ma,1,ma);
  coeff = vector(1,ma);
  coeff[1] = *a; coeff[2] = *m; coeff[3] = *s;
  ia = ivector(1,ma);
  ia[1] = aflag; ia[2] = mflag; ia[3] = sflag;

  alamda = -1.0; // For initialization (first) call to mrqmin
  mrqmin(x-1,y-1,sd-1,n,coeff,ia,ma,covar,alpha,&chisq,fgauss,&alamda);
  old_chisq = chisq; old_alamda = alamda;
  itcount = 0;
  done = 0;
  while(!done){
    mrqmin(x-1,y-1,sd-1,n,coeff,ia,ma,covar,alpha,&chisq,fgauss,&alamda);
    itcount += 1;
    diff = chisq - old_chisq;
    if (pflag)
      printf("    %4d %9.5f %9.5f\n",itcount,chisq,diff);
    if ((diff <= 0.0) && (diff > -0.0001) && (alamda <= old_alamda))
      done = 1;
    old_chisq = chisq; old_alamda = alamda;
  }
  alamda = 0.0;   // final call to get alpha, covar matrices, alamda=0.0
  mrqmin(x-1,y-1,sd-1,n,coeff,ia,ma,covar,alpha,&chisq,fgauss,&alamda);

  *a = coeff[1]; *m = coeff[2]; *s = coeff[3];
  free_matrix(covar,1,ma,1,ma); free_matrix(alpha,1,ma,1,ma);
  free_vector(coeff,1,ma);
  free_ivector(ia,1,ma);
}
/**************************************-**************************************/
/*                                                                           */
/*                              FIT_TO_GAUSSIAN_4                            */
/*                                                                           */
/*  Like 'fit_to_gaussian' above, but 1/2 factor is in the exponent, and     */
/*  a y-offset value is included (thus, 4 parameters can vary).              */
/*                                                                           */
/*                          (x-M)^2                                          */
/*       f(x) = C + A exp - -------                                          */
/*                          2 * S^2                                          */
/*                                                                           */
/*  NOTES:                                                                   */
/*  1.  chisq is returned.                                                   */
/*                                                                           */
/*****************************************************************************/
void fit_to_gaussian_4(x,y,sd,n,a,m,s,c,aflag,mflag,sflag,cflag,rchisq,rprob)
     float *x,*y,*sd;             /* Data with standard deviations. */
     int n;                       /* Number of data points. */
     float *a,*m,*s,*c;           /* Initial guesses for parameters */
     int aflag,mflag,sflag,cflag; /* Whether to fit or hold param. */
     float *rchisq,*rprob;
{
  int i;
  int done,itcount,ma,*ia,pflag;
  float diff,old_chisq,old_alamda;
  float alamda,chisq,*coeff,**covar,**alpha;

  pflag = 0;
  if (pflag)
    printf("  FIT_TO_GAUSSIAN_4\n");
  
  for(i=0;i<n;i++)
    if (sd[i] <= 0.0)
      exit_error("FIT_TO_GAUSSIAN_4","Standard deviation <= 0.0");

  ma = 4; /*** PREPARATION FOR "mrqmin" ***/
  covar = matrix(1,ma,1,ma); alpha = matrix(1,ma,1,ma);
  coeff = vector(1,ma);
  coeff[1] = *a; coeff[2] = *m; coeff[3] = *s; coeff[4] = *c;
  ia = ivector(1,ma);
  ia[1] = aflag; ia[2] = mflag; ia[3] = sflag; ia[4] = cflag;

  alamda = -1.0; /*** For initialization (first) call to mrqmin ***/
  mrqmin(x-1,y-1,sd-1,n,coeff,ia,ma,covar,alpha,&chisq,fgauss4,&alamda);
  old_chisq = chisq; old_alamda = alamda;
  itcount = 0;
  done = 0;
  while(!done){
    mrqmin(x-1,y-1,sd-1,n,coeff,ia,ma,covar,alpha,&chisq,fgauss4,&alamda);
    itcount += 1;
    diff = chisq - old_chisq;
    if (pflag)
      printf("    %4d %9.5f %9.5f\n",itcount,chisq,diff);
    if ((diff <= 0.0) && (diff > -0.0001) && (alamda <= old_alamda))
      done = 1;
    old_chisq = chisq; old_alamda = alamda;
  }
  alamda = 0.0;   /* final call to get alpha, covar matrices, alamda=0.0 */
  mrqmin(x-1,y-1,sd-1,n,coeff,ia,ma,covar,alpha,&chisq,fgauss4,&alamda);

  *rprob = gammq(0.5*(float)(n-ma),0.5*chisq); /* page 660 NumRecInC 2nd Ed. */
  *a = coeff[1]; *m = coeff[2]; *s = coeff[3]; *c = coeff[4]; *rchisq = chisq;
  free_matrix(covar,1,ma,1,ma); free_matrix(alpha,1,ma,1,ma);
  free_vector(coeff,1,ma);
  free_ivector(ia,1,ma);
}
/**************************************-**************************************/
/*                                                                           */
/*                              FIT_TO_FLAT_LINE                             */
/*                                                                           */
/*  Fit to a flat (horizontal) line.                                         */
/*                                                                           */
/*****************************************************************************/
void fit_to_flat_line(x,y,sd,n,a,rchisq,rprob)
     float *x,*y,*sd;             /* Data with standard deviations. */
     int n;                       /* Number of data points. */
     float *a;                    /* Initial guesses for constant value */
     float *rchisq,*rprob;
{
  int i;
  int done,itcount,ma,*ia,pflag;
  float diff,old_chisq,old_alamda;
  float alamda,chisq,*coeff,**covar,**alpha;

  pflag = 0;
  if (pflag)
    printf("  FIT_TO_FLAT_LINE\n");
  
  for(i=0;i<n;i++)
    if (sd[i] <= 0.0)
      exit_error("FIT_TO_FLAT_LINE","Standard deviation <= 0.0");

  ma = 1; /*** PREPARATION FOR "mrqmin" ***/
  covar = matrix(1,ma,1,ma); alpha = matrix(1,ma,1,ma);
  coeff = vector(1,ma);
  coeff[1] = *a;
  ia = ivector(1,ma);
  ia[1] = 1;

  alamda = -1.0; /*** For initialization (first) call to mrqmin ***/
  mrqmin(x-1,y-1,sd-1,n,coeff,ia,ma,covar,alpha,&chisq,flat_line,&alamda);
  old_chisq = chisq; old_alamda = alamda;
  itcount = 0;
  done = 0;
  while(!done){
    mrqmin(x-1,y-1,sd-1,n,coeff,ia,ma,covar,alpha,&chisq,flat_line,&alamda);
    itcount += 1;
    diff = chisq - old_chisq;
    if (pflag)
      printf("    %4d %9.5f %9.5f  %.8f\n",itcount,chisq,diff,coeff[1]);
    if ((diff <= 0.0) && (diff > -0.0001) && (alamda <= old_alamda))
      done = 1;
    old_chisq = chisq; old_alamda = alamda;
  }
  alamda = 0.0;   /* final call to get alpha, covar matrices, alamda=0.0 */
  mrqmin(x-1,y-1,sd-1,n,coeff,ia,ma,covar,alpha,&chisq,flat_line,&alamda);

  *rprob = gammq(0.5*(float)(n-ma),0.5*chisq); /* page 660 NumRecInC 2nd Ed. */
  *a = coeff[1]; *rchisq = chisq;
  free_matrix(covar,1,ma,1,ma); free_matrix(alpha,1,ma,1,ma);
  free_vector(coeff,1,ma);
  free_ivector(ia,1,ma);
}
/**************************************-**************************************/
/*                                                                           */
/*                               FIT_TO_SIGMOID                              */
/*                                                                           */
/*****************************************************************************/
void fit_to_sigmoid(x,y,sd,n,a,b,d,aflag,bflag,dflag,rdf,rchisq,rprob)
     float *x,*y,*sd;        /* Data with standard deviations. */
     int n;                  /* Number of data points. */
     float *a,*b,*d;         /* Initial guesses for parameters */
     int aflag,bflag,dflag;  /* Whether to fit or hold param. */
     int *rdf;
     float *rchisq,*rprob;
{
  int i;
  int done,itcount,ma,*ia,pflag;
  float diff,old_chisq,old_alamda;
  float alamda,chisq,*coeff,**covar,**alpha;

  pflag = 0;
  if (pflag)
    printf("  FIT_TO_SIGMOID\n");
  
  for(i=0;i<n;i++)
    if (sd[i] <= 0.0)
      exit_error("FIT_TO_SIGMOID","Standard deviation <= 0.0");

  ma = 3; /*** PREPARATION FOR "mrqmin" ***/
  covar = matrix(1,ma,1,ma); alpha = matrix(1,ma,1,ma);
  coeff = vector(1,ma);
  coeff[1] = *a; coeff[2] = *b; coeff[3] = *d;
  ia = ivector(1,ma);
  ia[1] = aflag; ia[2] = bflag; ia[3] = dflag;

  alamda = -1.0; /*** For initialization (first) call to mrqmin ***/
  mrqmin(x-1,y-1,sd-1,n,coeff,ia,ma,covar,alpha,&chisq,sigpsych,&alamda);
  old_chisq = chisq;
  old_alamda = alamda;
  itcount = 0;
  done = 0;
  while(!done){
    mrqmin(x-1,y-1,sd-1,n,coeff,ia,ma,covar,alpha,&chisq,sigpsych,&alamda);
    itcount += 1;
    diff = chisq - old_chisq;
    if (pflag)
      printf("    %4d %9.5f %9.5f %f\n",itcount,chisq,diff,alamda);
    if ((diff <= 0.0) && (diff > -0.0001) && (alamda <= old_alamda))
      done = 1;
    old_chisq = chisq;
    old_alamda = alamda;
  }
  alamda = 0.0;   /* final call to get alpha, covar matrices, alamda=0.0 */
  mrqmin(x-1,y-1,sd-1,n,coeff,ia,ma,covar,alpha,&chisq,sigpsych,&alamda);

  *rdf = n-ma; /* Degrees of freedom. */
  *rprob = gammq(0.5*(float)(n-ma),0.5*chisq); /* page 660 NumRecInC 2nd Ed. */
  *a = coeff[1]; *b = coeff[2]; *d = coeff[3]; *rchisq = chisq;
  free_matrix(covar,1,ma,1,ma); free_matrix(alpha,1,ma,1,ma);
  free_vector(coeff,1,ma);
  free_ivector(ia,1,ma);
}
/**************************************-**************************************/
/*                                                                           */
/*                                 FIT_TO_SINE                               */
/*                                                                           */
/*  Fit to:                                                                  */
/*                                                                           */
/*    y = b + a*sin(fx-p)                                                    */
/*                                                                           */
/*  NumRec params:                                                           */
/*    a[0] = a = amplitude                                                   */
/*    a[1] = b = vertical offset                                             */
/*    a[2] = p = phase                                                       */
/*    a[3] = f = frequency                                                   */
/*                                                                           */
/*****************************************************************************/
void fit_to_sine(x,y,sd,n,a,b,p,f,aflag,bflag,pflag,fflag,rdf,rchisq,rprob)
     float *x,*y,*sd;        /* Data with standard deviations. */
     int n;                  /* Number of data points. */
     float *a,*b,*p,*f;           /* Initial guesses for parameters */
     int aflag,bflag,pflag,fflag; /* Whether to fit or hold param. */
     int *rdf;
     float *rchisq,*rprob;
{
  int i;
  int done,itcount,ma,*ia,prflag;
  float diff,old_chisq,old_alamda;
  float alamda,chisq,*coeff,**covar,**alpha;

  prflag = 0;
  if (prflag)
    printf("  FIT_TO_SINE\n");
  
  for(i=0;i<n;i++)
    if (sd[i] <= 0.0)
      exit_error("FIT_TO_SINE","Standard deviation <= 0.0");

  ma = 4; /*** PREPARATION FOR "mrqmin" ***/
  covar = matrix(1,ma,1,ma); alpha = matrix(1,ma,1,ma);
  coeff = vector(1,ma);
  coeff[1] = *a; coeff[2] = *b; coeff[3] = *p; coeff[4] = *f;
  ia = ivector(1,ma);
  ia[1] = aflag; ia[2] = bflag; ia[3] = pflag; ia[4] = fflag;

  /**printf("n=%d\n",n);
    for(i=0;i<n;i++)
    printf("%f %f %f\n",x[i],y[i],sd[i]);**/

  alamda = -1.0; /*** For initialization (first) call to mrqmin ***/
  mrqmin(x-1,y-1,sd-1,n,coeff,ia,ma,covar,alpha,&chisq,fsin4,&alamda);
  old_chisq = chisq;
  old_alamda = alamda;
  itcount = 0;
  done = 0;
  while(!done){
    mrqmin(x-1,y-1,sd-1,n,coeff,ia,ma,covar,alpha,&chisq,fsin4,&alamda);
    itcount += 1;
    diff = chisq - old_chisq;
    if (prflag)
      printf("    %4d %9.5f %9.5f %f\n",itcount,chisq,diff,alamda);
    if ((diff <= 0.0) && (diff > -0.0001) && (alamda <= old_alamda))
      done = 1;
    old_chisq = chisq;
    old_alamda = alamda;
  }
  alamda = 0.0;   /* final call to get alpha, covar matrices, alamda=0.0 */
  mrqmin(x-1,y-1,sd-1,n,coeff,ia,ma,covar,alpha,&chisq,fsin4,&alamda);

  *rdf = n-ma; /* Degrees of freedom. */
  *rprob = gammq(0.5*(float)(n-ma),0.5*chisq); /* page 660 NumRecInC 2nd Ed. */
  *a = coeff[1]; *b = coeff[2]; *p = coeff[3]; *f = coeff[4]; *rchisq = chisq;
  free_matrix(covar,1,ma,1,ma); free_matrix(alpha,1,ma,1,ma);
  free_vector(coeff,1,ma);
  free_ivector(ia,1,ma);
}
/**************************************-**************************************/
/*                                                                           */
/*                              FIT_TO_DECAY_EXP                             */
/*                                                                           */
/*  Fit to:                                                                  */
/*                                                                           */
/*    y = b + a * e^(- x/t)                                                  */
/*                                                                           */
/*  NumRec params:                                                           */
/*    a[0] = a = amplitude                                                   */
/*    a[1] = b = vertical offset                                             */
/*    a[2] = t = time constant (tau)                                         */
/*                                                                           */
/*****************************************************************************/
void fit_to_decay_exp(x,y,sd,n,a,b,t,aflag,bflag,tflag,rdf,rchisq,rprob)
     float *x,*y,*sd;        // Data with standard deviations
     int n;                  // Number of data points
     float *a,*b,*t;         // Initial guesses for parameters
     int aflag,bflag,tflag;  // Whether to fit or hold param
     int *rdf;
     float *rchisq,*rprob;
{
  int i;
  int done,itcount,ma,*ia,prflag;
  float diff,old_chisq,old_alamda;
  float alamda,chisq,*coeff,**covar,**alpha;

  prflag = 0;
  if (prflag)
    printf("  FIT_TO_DECAY_EXP\n");
  
  for(i=0;i<n;i++)
    if (sd[i] <= 0.0)
      exit_error("FIT_TO_DECAY_EXP","Standard deviation <= 0.0");

  ma = 3; // PREPARATION FOR "mrqmin"
  covar = matrix(1,ma,1,ma); alpha = matrix(1,ma,1,ma);
  coeff = vector(1,ma);
  coeff[1] = *a; coeff[2] = *b; coeff[3] = *t;
  ia = ivector(1,ma);
  ia[1] = aflag; ia[2] = bflag; ia[3] = tflag;

  /**printf("n=%d\n",n);
    for(i=0;i<n;i++)
    printf("%f %f %f\n",x[i],y[i],sd[i]);**/

  alamda = -1.0; // For initialization (first) call to mrqmin
  mrqmin(x-1,y-1,sd-1,n,coeff,ia,ma,covar,alpha,&chisq,fexp3,&alamda);
  old_chisq = chisq;
  old_alamda = alamda;
  itcount = 0;
  done = 0;
  while(!done){
    mrqmin(x-1,y-1,sd-1,n,coeff,ia,ma,covar,alpha,&chisq,fexp3,&alamda);
    itcount += 1;
    diff = chisq - old_chisq;
    if (prflag)
      printf("    %4d %9.5f %9.5f %f\n",itcount,chisq,diff,alamda);
    if ((diff <= 0.0) && (diff > -0.0001) && (alamda <= old_alamda))
      done = 1;
    old_chisq = chisq;
    old_alamda = alamda;
  }
  alamda = 0.0;   // final call to get alpha, covar matrices, alamda=0.0
  mrqmin(x-1,y-1,sd-1,n,coeff,ia,ma,covar,alpha,&chisq,fexp3,&alamda);

  *rdf = n-ma; // Degrees of freedom.
  *rprob = gammq(0.5*(float)(n-ma),0.5*chisq); // page 660 NumRecInC 2nd Ed.
  *a = coeff[1]; *b = coeff[2]; *t = coeff[3]; *rchisq = chisq;
  free_matrix(covar,1,ma,1,ma); free_matrix(alpha,1,ma,1,ma);
  free_vector(coeff,1,ma);
  free_ivector(ia,1,ma);
}
/**************************************-**************************************/
/*                                                                           */
/*                          CONVOLVE_WITH_MASK_RENORM                        */
/*                                                                           */
/*  NOTE: This convolution routine is used to handle edges with Gaussian     */
/*  smoothing---NOT FOR GENERAL USE.                                         */
/*                                                                           */
/*  NOTE:  "mask" is reversed during convolution.                            */
/*                                                                           */
/*  WARNING:  "mask" is assumed to have odd "width" since the center value   */
/*  is taken as the origin of the mask.                                      */
/*                                                                           */
/*****************************************************************************/
float *convolve_with_mask_renorm(data,n,mask,width)
     float *data;
     int n;
     float *mask;
     int width;
{
  int i,j;
  int p,mid;
  float *s,wtot,wpart;

  s = get_zero_farray(n);

  // Compute the total area (weight) of the mask
  wtot = 0.0;
  for (i=0;i<width;i++)
    wtot += mask[i];

  mid = (width-1)/2;  // Define the origin of the mask (the mid point)
  for(i=0;i<n;i++){
    if ((i<mid) || (i>=(n-mid))){ // if the mask over-extends data
      wpart = 0.0;  // Compute the partial weight that is used
      for(j=0;j<width;j++){
	p = i + (mid-j); // "+" reverses mask during convolution!
	if ((p >= 0) && (p < n)){  // If we are within the data
	  s[i] += (mask[j] * data[p]);  // sum up the product
	  wpart += mask[j];  // Update the weight used
	}
      }
      if (wpart != 0.0)
	s[i] *= wtot/wpart;  // Divide by the weight used
    }else{
      //  The mask is completely contained within the data, so this is easy
      for(j=0;j<width;j++)
	s[i] += (mask[j] * data[i+(mid-j)]);
    }
  }
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                            CONVOLVE_WITH_MASK                             */
/*                                                                           */
/*  Only use the portion of the mask that overlaps valid data.               */
/*                                                                           */
/*  NOTE:  "mask" is reversed during convolution.                            */
/*                                                                           */
/*  WARNING:  "mask" is assumed to have odd "width" since the center value   */
/*  is taken as the origin of the mask.                                      */
/*                                                                           */
/*****************************************************************************/
float *convolve_with_mask(data,n,mask,width)
     float *data;
     int n;
     float *mask;
     int width;
{
  int i,j;
  int p,mid;
  float *s;

  s = get_zero_farray(n);

  mid = (width-1)/2;
  for(i=0;i<n;i++){
    if ((i<mid) || (i>=(n-mid))){ // mask over-extends data
      for(j=0;j<width;j++){
	p = i + (mid-j); // "+" reverses mask during convolution!
	if ((p >= 0) && (p < n))
	  s[i] += (mask[j] * data[p]);
      }
    }else
      for(j=0;j<width;j++)
	s[i] += (mask[j] * data[i+(mid-j)]);
  }
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                         CONVOLVE_WITH_MASK_CAUSAL                         */
/*                                                                           */
/*  Only use the portion of the mask that overlaps valid data.               */
/*  NOTE:  "mask" is reversed during convolution.                            */
/*  NOTE:  value of data is implicitly zero before 'data' begins.            */
/*  WARNING:  "mask" is assumed to have the time origin at the 0th index.    */
/*                                                                           */
/*****************************************************************************/
float *convolve_with_mask_causal(data,n,mask,width)
     float *data;
     int n;
     float *mask;
     int width;
{
  int i,j;
  int wm1;
  float *s;

  if (n <= width)
    exit_error("CONVOLVE_WITH_MASK_CAUSAL","data is not longer than mask");

  s = get_zero_farray(n);

  wm1 = width - 1;
  for(i=0;i<wm1;i++){
    for(j=0;j<=i;j++){
      s[i] += (mask[j] * data[i-j]);
    }
  }
  for(i=wm1;i<n;i++)
    for(j=0;j<width;j++)
      s[i] += (mask[j] * data[i-j]);
  
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                       CONVOLVE_WITH_MASK_CAUSAL_X0                        */
/*                                                                           */
/*  Assume the value of data is 'data[0]' before the data begins.            */
/*  See notes above for 'convolve_with_mask_causal'.                         */
/*                                                                           */
/*****************************************************************************/
float *convolve_with_mask_causal_x0(data,n,mask,width)
     float *data;
     int n;
     float *mask;
     int width;
{
  int i,j;
  int wm1;
  float *s;

  if (n <= width)
    exit_error("CONVOLVE_WITH_MASK_CAUSAL_X0","data is not longer than mask");

  s = get_zero_farray(n);

  wm1 = width - 1;
  for(i=0;i<wm1;i++){
    for(j=0;j<=i;j++){
      s[i] += (mask[j] * data[i-j]);
    }
    for(j=i+1;j<width;j++){
      s[i] += (mask[j] * data[0]);  // Assume data[t<0] = data[0]
    }
  }
  for(i=wm1;i<n;i++)
    for(j=0;j<width;j++)
      s[i] += (mask[j] * data[i-j]);
  
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                         CONVOLVE_CIRCULAR_WITH_MASK                       */
/*                                                                           */
/*  Like "convolve_with_mask" but assumes the data is circular.              */
/*                                                                           */
/*  NOTE:  "mask" is reversed during convolution.                            */
/*                                                                           */
/*  WARNING:  "mask" is assumed to have odd "width" since the center value   */
/*  is taken as the origin of the mask.                                      */
/*                                                                           */
/*****************************************************************************/
float *convolve_circular_with_mask(data,n,mask,width)
     float *data;
     int n;
     float *mask;
     int width;
{
  int i,j,k;
  int mid;
  float *s;

  mid = (width-1)/2;
  s = get_zero_farray(n);
  for(i=0;i<n;i++){
    for(j=0;j<width;j++){
      k = i+(mid-j);
      while(k < 0)
	k += n;
      while(k >= n)
	k -= n;
      s[i] += (mask[j] * data[k]);
    }
  }
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                          CONVOLVE_SPARSE_WITH_MASK                        */
/*                                                                           */
/*  Only use the portion of the mask that overlaps valid data, which         */
/*  is indicated by a "1" in the "valid" array.                              */
/*                                                                           */
/*****************************************************************************/
float *convolve_sparse_with_mask(data,valid,n,mask,width)
     float *data;
     unsigned char *valid;
     int    n;
     float *mask;
     int    width;
{
  int     i,j;
  int     p;
  float  *s;
  int     mid;
  float   wtot,wpart; /* total and partial weight of mask */
  unsigned char *tvalid;

  tvalid = (unsigned char *)myalloc(n*sizeof(unsigned char));
  for (i=0;i<n;i++)
    tvalid[i] = 0;

  s = get_zero_farray(n);

  wtot = 0.0;
  for (i=0;i<width;i++)
    wtot += mask[i];

  mid = (width-1)/2;
  for(i=0;i<n;i++){
    wpart = 0.0;
    for(j=0;j<width;j++){     /* position in the mask */
      p = i+(mid-j);        /* position in the data */
      if ((p>=0) && (p<n) && (valid[p]==1)){
	tvalid[i] = 1;
	s[i] += (mask[j] * data[p]);
	wpart += mask[j];
      }
    }
    if (wpart != 0.0)
      s[i] *= wtot/wpart;
  }
  for (i=0;i<n;i++)
    valid[i] = tvalid[i];
  myfree(tvalid);
  return(s);
}
/**************************************-**************************************/
/*                                                                           */
/*                            SMOOTH_WITH_GAUSSIAN                           */
/*                                                                           */
/*  Smooth the array of floats with a discrete Gaussian mask of the          */
/*  specified standard deviation and tolerance.                              */
/*                                                                           */
/*****************************************************************************/
float *smooth_with_gaussian(data,n,sigma,tolerance)
     float *data;
     int n;
     float sigma,tolerance;
{
  int width;
  float *mask,*smooth;

  if (sigma > 0.0){
    discrete_gaussian(sigma,tolerance,&mask,&width);
    smooth = convolve_with_mask_renorm(data,n,mask,width);
    myfree(mask);
  }else
    smooth = copy_farray(data,n);

  return smooth;
}
/**************************************-**************************************/
/*                                                                           */
/*                         SMOOTH_WITH_GAUSSIAN_CIRCULAR                     */
/*                                                                           */
/*  Smooth the array of floats with a discrete Gaussian mask of the          */
/*  specified standard deviation and tolerance.                              */
/*                                                                           */
/*****************************************************************************/
float *smooth_with_gaussian_circular(data,n,sigma,tolerance)
     float *data;
     int n;
     float sigma,tolerance;
{
  int width;
  float *mask,*smooth;

  if (sigma > 0.0){
    discrete_gaussian(sigma,tolerance,&mask,&width);
    smooth = convolve_circular_with_mask(data,n,mask,width);
    myfree(mask);
  }else
    smooth = copy_farray(data,n);

  return smooth;
}
/**************************************-**************************************/
/*                                                                           */
/*                      SMOOTH_WITH_GAUSSIAN_OMIT_CENTER                     */
/*                                                                           */
/*  Smooth the array of floats with a discrete Gaussian mask of the          */
/*  specified standard deviation and tolerance, but set center value in the  */
/*  Gaussian mask to zero, and renormalize.                                  */
/*                                                                           */
/*****************************************************************************/
float *smooth_with_gaussian_omit_center(data,n,sigma,tolerance)
     float *data;
     int n;
     float sigma,tolerance;
{
  int width;
  float *mask,*smooth;

  if (sigma > 0.0){
    discrete_gaussian(sigma,tolerance,&mask,&width);
    mask[(width-1)/2] = 0.0;
    norm_area_farray(mask,width,1.0);
    smooth = convolve_with_mask_renorm(data,n,mask,width);
    myfree(mask);
  }else{
    smooth = NULL;
    exit_error("SMOOTH_WITH_GAUSSIAN_OMIT_CENTER","Sigma <= 0.0");
  }

  return smooth;
}
/**************************************-**************************************/
/*                                                                           */
/*                         SMOOTH_CIRCULAR_WITH_GAUSSIAN                     */
/*                                                                           */
/*  Smooth the array of floats with a discrete Gaussian mask of the          */
/*  specified standard deviation and tolerance.  Assume data lies on a       */
/*  circle.                                                                  */
/*                                                                           */
/*****************************************************************************/
float *smooth_circular_with_gaussian(data,n,sigma,tolerance)
     float *data;
     int n;
     float sigma,tolerance;
{
  int width;
  float *mask,*smooth;

  if (sigma > 0.0){
    discrete_gaussian(sigma,tolerance,&mask,&width);
    smooth = convolve_circular_with_mask(data,n,mask,width);
    myfree(mask);
  }else
    smooth = copy_farray(data,n);

  return smooth;
}
/**************************************-**************************************/
/*                                                                           */
/*                         SMOOTH_SPARSE_WITH_GAUSSIAN                       */
/*                                                                           */
/*  Smooth the array of floats with a discrete Gaussian mask of the          */
/*  specified standard deviation and tolerance.                              */
/*                                                                           */
/*****************************************************************************/
float *smooth_sparse_with_gaussian(data,valid,n,sigma,tolerance)
     float *data;
     unsigned char *valid;
     int n;
     float sigma,tolerance;
{
  float *mask;
  int width;
  float *smooth;
  
  discrete_gaussian(sigma,tolerance,&mask,&width);
  smooth = convolve_sparse_with_mask(data,valid,n,mask,width);

  myfree(mask);
  return(smooth);
}
/**************************************-**************************************/
/*                                                                           */
/*                           SMOOTH_2D_WITH_GAUSSIAN                         */
/*                                                                           */
/*  Smoothing by a Gaussian in 2D is done by smoothing in horizontal then    */
/*  vertical with 1D Gaussians, both with standard deviation sigma.          */
/*                                                                           */
/*****************************************************************************/
float **smooth_2d_with_gaussian(data,m,n,sigma,tolerance)
     float **data;
     int m,n;
     float sigma,tolerance;
{
  int i,j;
  float *column,*smooth_column;
  float **smooth;

  if (m*n > 100000) printf("  SMOOTH_2D_WITH_GAUSSIAN\n");
  
  smooth = (float **)myalloc(m*sizeof(float *));
  if (m*n > 100000) printf("    Horizontal:");
  for (i=0;i<m;i++){
    if ((i%100 == 0) && (m*n > 100000)){
      printf(" %d",i);
      fflush(stdout);
    }
    smooth[i] = smooth_with_gaussian(data[i],n,sigma,tolerance);
  }
  if (m*n > 100000) printf("\n    Vertical:");
  column = (float *)myalloc(m*sizeof(float));
  for (j=0;j<n;j++){
    if ((j%100 == 0) && (m*n > 100000)){
      printf(" %d",j);
      fflush(stdout);
    }
    for (i=0;i<m;i++)
      column[i] = smooth[i][j];
    smooth_column = smooth_with_gaussian(column,m,sigma,tolerance);
    for (i=0;i<m;i++)
      smooth[i][j] = smooth_column[i];
    myfree(smooth_column);
  }
  if (m*n > 100000) printf("\n");
  myfree(column);
  return(smooth);
}
/**************************************-**************************************/
/*                                                                           */
/*                       SMOOTH_2D_WITH_GAUSSIAN_CIRCULAR                    */
/*                                                                           */
/*  Smoothing by a Gaussian in 2D is done by smoothing in horizontal then    */
/*  vertical with 1D Gaussians, both with standard deviation sigma.          */
/*                                                                           */
/*  Wrap around.                                                             */
/*                                                                           */
/*****************************************************************************/
float **smooth_2d_with_gaussian_circular(data,m,n,sigma,tolerance)
     float **data;
     int m,n;
     float sigma,tolerance;
{
  int i,j;
  float *column,*smooth_column;
  float **smooth;

  //if (m*n > 10000) printf("  SMOOTH_2D_WITH_GAUSSIAN\n");
  
  smooth = (float **)myalloc(m*sizeof(float *));
  //if (m*n > 10000) printf("    Horizontal:");
  for (i=0;i<m;i++){
    /*
    if ((i%100 == 0) && (m*n > 10000)){
      printf(" %d",i);
      fflush(stdout);
    }*/
    smooth[i] = smooth_with_gaussian_circular(data[i],n,sigma,tolerance);
  }
  //if (m*n > 10000) printf("\n    Vertical:");
  column = (float *)myalloc(m*sizeof(float));
  for (j=0;j<n;j++){
    /*
    if ((j%100 == 0) && (m*n > 10000)){
      printf(" %d",j);
      fflush(stdout);
    }
    */
    for (i=0;i<m;i++)
      column[i] = smooth[i][j];
    smooth_column = smooth_with_gaussian_circular(column,m,sigma,tolerance);
    for (i=0;i<m;i++)
      smooth[i][j] = smooth_column[i];
    myfree(smooth_column);
  }
  //if (m*n > 10000) printf("\n");
  myfree(column);
  return(smooth);
}
/**************************************-**************************************/
/*                                                                           */
/*                       SMOOTH_SPARSE_2D_WITH_GAUSSIAN                      */
/*                                                                           */
/*  Smoothing by a Gaussian in 2D is done by smoothing in horizontal         */
/*  then vertical with 1D Gaussians, both with standard deviation            */
/*  sigma.                                                                   */
/*                                                                           */
/*****************************************************************************/
float **smooth_sparse_2d_with_gaussian(data,valid,m,n,sigma,tolerance)
     float **data;
     unsigned char **valid;
     int m,n;
     float sigma,tolerance;
{
  int i,j;
  float *column,*smooth_column;
  unsigned char *vcolumn;
  float **smooth;

  printf("  SMOOTH_SPARSE_2D_WITH_GAUSSIAN\n");

  smooth = (float **)myalloc(m*sizeof(float *));
  printf("    Horizontal:");
  for (i=0;i<m;i++){
    if (i%100 == 0){
      printf(" %d",i);
      fflush(stdout);
    }
    smooth[i] = smooth_sparse_with_gaussian(data[i],valid[i],n,
					    sigma,tolerance);
  }
  printf("\n    Vertical:");
  column = (float *)myalloc(m*sizeof(float));
  vcolumn = (unsigned char *)myalloc(m*sizeof(unsigned char));
  for (j=0;j<n;j++){
    if (j%100 == 0){
      printf(" %d",j);
      fflush(stdout);
    }
    for (i=0;i<m;i++){
      column[i] = smooth[i][j];
      vcolumn[i] = valid[i][j];
    }
    smooth_column = smooth_sparse_with_gaussian(column,vcolumn,m,sigma,
						tolerance);
    for (i=0;i<m;i++)
      smooth[i][j] = smooth_column[i];
    myfree(smooth_column);
  }
  printf("\n");
  myfree(column);
  myfree(vcolumn);
  return(smooth);
}
/**************************************-**************************************/
/*                                                                           */
/*                          SMOOTH_3D_2D_WITH_GAUSSIAN                       */
/*                                                                           */
/*  Assuming data[xn][yn][tn], smooth each x,y frame in 2D.                  */
/*                                                                           */
/*****************************************************************************/
float ***smooth_3d_2d_with_gaussian(data,xn,yn,tn,sigma,tolerance)
     float ***data;
     int xn,yn,tn;
     float sigma,tolerance;
{
  int i,j,k;
  float ***dsm,**d2d,**d2s;

  dsm = get_3d_farray(xn,yn,tn);

  for(k=0;k<tn;k++){
    d2d = get_2d_from_3d_farray(data,0,xn,0,yn,k);

    d2s = smooth_2d_with_gaussian(d2d,xn,yn,sigma,tolerance);
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	dsm[i][j][k] = d2s[i][j];

    free_2d_farray(d2d,xn);
    free_2d_farray(d2s,xn);
  }

  return dsm;
}
/**************************************-**************************************/
/*                                                                           */
/*                          DERIVATIVE_SIMPLE_FARRAY                         */
/*                                                                           */
/*  Compute the first derivative of the farray using:  1/(2 dx) [-1 1].      */
/*                                                                           */
/*****************************************************************************/
float *derivative_simple_farray(data,n)
     float *data;
     int n;
{
  int i;
  float *deriv;

  if (n > 1){
    deriv = (float *)myalloc(n*sizeof(float));
    deriv[0] = 0;
    for(i=1;i<n;i++)
      deriv[i] = data[i] - data[i-1];
  }else{
    deriv = NULL;
    exit_error("DERIVATIVE_SIMPLE_FARRAY","Array too short\n");
  }

  return deriv;
}
/**************************************-**************************************/
/*                                                                           */
/*                              DERIVATIVE_FARRAY                            */
/*                                                                           */
/*  Compute the first derivative of the farray using:  1/(2 dx) [-1 0 1],    */
/*  and 1/dx [-1 1] at the the endpoints (dx = 1).                           */
/*                                                                           */
/*****************************************************************************/
float *derivative_farray(data,n)
     float *data;
     int n;
{
  int i;
  float *deriv;

  if (n > 1){
    deriv = (float *)myalloc(n*sizeof(float));
    deriv[0] = data[1] - data[0];
    deriv[n-1] = data[n-1] - data[n-2];
    for(i=1;i<(n-1);i++)
      deriv[i] = 0.5 * (data[i+1] - data[i-1]);
  }else{
    deriv = NULL;
    exit_error("DERIVATIVE_FARRAY","Array too short\n");
  }

  return deriv;
}
/**************************************-**************************************/
/*                                                                           */
/*                          DERIVATIVE_SMOOTH_FARRAY                         */
/*                                                                           */
/*  Compute the first derivative of the farray using:  1/(2 dx) [-1 0 1],    */
/*  and 1/dx [-1 1] at the the endpoints (dx = 1).                           */
/*                                                                           */
/*****************************************************************************/
float *derivative_smooth_farray(data,n,sigma,tolerance)
     float *data;
     int n;
     float sigma,tolerance;
{
  float *deriv,*t;

  t = smooth_with_gaussian(data,n,sigma,tolerance);
  deriv = derivative_farray(t,n);
  myfree(t);

  return deriv;
}
/**************************************-**************************************/
/*                                                                           */
/*                        DERIVATIVE_SMOOTH_2D_FARRAY                        */
/*                                                                           */
/*****************************************************************************/
float **derivative_smooth_2d_farray(data,m,n,sigma,tolerance)
     float **data;
     int m,n;
     float sigma,tolerance;
{
  int i;
  float **deriv;

  deriv = (float **)myalloc(m*sizeof(float *));
  for(i=0;i<m;i++)
    deriv[i] = derivative_smooth_farray(data[i],n,sigma,tolerance);

  return deriv;
}
/**************************************-**************************************/
/*                                                                           */
/*                            SECOND_DERIVATIVE_FARRAY                       */
/*                                                                           */
/*  Compute the second derivative of the farray using                        */
/*    1/dx^2 [1 -2 1]                                                        */
/*  and use 0 at the endpoints (dx = 1).                                     */
/*                                                                           */
/*****************************************************************************/
float *second_derivative_farray(data,n)
     float *data;
     int n;
{
  int i;
  float *deriv;

  if (n > 1){
    deriv = (float *)myalloc(n*sizeof(float));
    deriv[0] = deriv[n-1] = 0.0;
    for(i=1;i<(n-1);i++)
      deriv[i] = data[i+1] -2.0*data[i] + data[i-1];
  }else{
    deriv = NULL;
    exit_error("SECOND_DERIVATIVE_FARRAY","Array too short\n");
  }

  return deriv;
}
/**************************************-**************************************/
/*                                                                           */
/*                       GET_SMOOTH_ABS_DERIV_XY_DATA                        */
/*                                                                           */
/*  From two farrays which contain horizontal and vertical positions,        */
/*  compute the absolute value of the velocity.                              */
/*                                                                           */
/*****************************************************************************/
float *get_smooth_abs_deriv_xy_data(xdata,ydata,n,sigma)
     float *xdata,*ydata;
     int n;
     float sigma;
{
  int i;
  float *dx,*dy,*dabs,*txdata,*tydata;

  if (sigma > 0.0){ /*** Smooth the data. ***/
    txdata = smooth_with_gaussian(xdata,n,sigma,0.01);
    tydata = smooth_with_gaussian(ydata,n,sigma,0.01);
    dx = derivative_farray(txdata,n);
    dy = derivative_farray(tydata,n);
    myfree(txdata); myfree(tydata);
  }else{
    dx = derivative_farray(xdata,n);
    dy = derivative_farray(ydata,n);
  }
  dabs = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    dabs[i] = sqrt(dx[i]*dx[i] + dy[i]*dy[i]);
  myfree(dx); myfree(dy);

  return dabs;
}
/**************************************-**************************************/
/*                                                                           */
/*                      GET_PARADIAGONAL_2D_FLOAT                            */
/*                                                                           */
/*  We imagine the 2D array to have its first subscript along the            */
/*  horizontal axis and its second subscript along the vertical              */
/*  axis of the usual Euclidean plane.  The array of points within           */
/*  the bounds of the square array along the line y=x-lag is                 */
/*  returned.  For example, in this array                                    */
/*                                                                           */
/*                  4  a b c d e                                             */
/*                  3  b c d e f                                             */
/*                  2  c d e f g                                             */
/*                  1  d e f g h                                             */
/*                  0  e f g h i                                             */
/*                                                                           */
/*                     0 1 2 3 4                                             */
/*                                                                           */
/*  The f's represent the paradiagonal at lag 1, while the d's               */
/*  represent the paradiagonal at lag -1.                                    */
/*                                                                           */
/*****************************************************************************/
float *get_paradiagonal_2d_float(data,n,lag,npd)
     float **data;  // [n][n]
     int n,lag;
     int *npd; // the length of the paradiagonal array returned
{
  int i;
  float *pd;

  if (abs(lag) >= n)
    exit_error("GET_PARADIAGONAL_2D_FLOAT","Lag too big");

  *npd = n-abs(lag);
  pd = (float *)myalloc(*npd*sizeof(float));
  for (i=0;i<*npd;i++)
    if (lag >= 0)
      pd[i] = data[i+lag][i];
    else
      pd[i] = data[i][-lag+i];

  return pd;
}
/**************************************-**************************************/
/*                                                                           */
/*                      GET_BAND_DIAGONAL_2D_FLOAT                           */
/*                                                                           */
/*  See comments in "get_paradiagonal_2d_float".  Here, we average           */
/*  perpendicular to the diagonal in the band [left,right].                  */
/*                                                                           */
/*    . . . * o .   The '*' and 'o' pattern represents the band              */
/*    . . * o * o   from [-2,1] and results in a length 9 array              */
/*    . * o * o .   of values averaged along a line of slope '\',            */
/*    * o * o . .   i.e., perpendicular to the main diagonal.                */
/*    o * o . . .                                                            */
/*    . o . . . .                                                            */
/*                                                                           */
/*****************************************************************************/
float *get_band_diagonal_2d_float(data,n,left,right,npd)
     float **data;
     int n,left,right; /* average the diagonal over [left,right] */
     int *npd; /* the length of the paradiagonal array returned */
{
  int i;
  int lag,tn,start,max_abs;
  float *pd,*temp;


  if (abs(left) > abs(right))
    max_abs = abs(left);
  else
    max_abs = abs(right);

  if (max_abs >= n)
    exit_error("GET_BAND_DIAGONAL_2D_FLOAT","Lag too big");

  *npd = n - max_abs;
  pd = get_zero_farray(*npd);

  for (lag=left;lag<=right;lag++){
    temp = get_paradiagonal_2d_float(data,n,lag,&tn);
    start = (max_abs - abs(lag)) / 2;
    for(i=0;i<*npd;i++)
      pd[i] += temp[i+start];
    myfree(temp);
  }
  for (i=0;i<*npd;i++)
    pd[i] /= (float)(right-left+1); // normalize for average

  return pd;
}
/**************************************-**************************************/
/*                                                                           */
/*                    GET_PARADIAG_PROJECTION_2D_FLOAT                       */
/*                                                                           */
/*  This is used for obtaining the Cross-correlogram from a 2d float         */
/*  array based on the JPSTH.  The projection is a length normalized         */
/*  average.                                                                 */
/*                                                                           */
/*****************************************************************************/
float *get_paradiag_projection_2d_float(data,n,max_lag)
     float **data;
     int n,max_lag;
{
  int i;
  float *pd,*cc;
  int npd;
  float sdev;

  if (max_lag >= n)
    exit_error("GET_PARADIAG_PROJECTION_2D_FLOAT","Lag too big");

  cc = (float *)myalloc((2*max_lag+1)*sizeof(float));
  for(i=-max_lag;i<=max_lag;i++){
    pd = get_paradiagonal_2d_float(data,n,i,&npd);
    mean_sdev_farray(pd,npd,&cc[i+max_lag],&sdev);
    myfree(pd);
  }
  return cc;
}
/**************************************-**************************************/
/*                                                                           */
/*                               PROJECT32                                   */
/*                                                                           */
/*  Project, i.e. average, a 3-d matrix onto a 2-d matrix w.r.t. a           */
/*  particular dimension.                                                    */
/*                                                                           */
/*****************************************************************************/
float **project32(data,n1,n2,n3,dim,rd1,rd2)
     float ***data;
     int n1,n2,n3,dim;
     int *rd1,*rd2; /* return dimensions of returned array */
{
  int i,j,k;
  float **p,sum;
  int pn1,pn2;

  /*printf("  PROJECT32\n");*/

  if (dim==1){
    pn1=n2; pn2=n3;
  }else if (dim==2){
    pn1=n1; pn2=n3;
  }else if (dim==3){
    pn1=n1; pn2=n2;
  }else{
    pn1 = pn2 = 0;
    exit_error("PROJECT32","dimension error");
  }
  
  p = (float **)myalloc(pn1*sizeof(float *));
  for (i=0;i<pn1;i++)
    p[i] = (float *)myalloc(pn2*sizeof(float));

  if (dim==1)
    for(i=0;i<n2;i++)
      for(j=0;j<n3;j++){
	sum = 0.0;
	for(k=0;k<n1;k++)
	  sum += data[k][i][j];
	p[i][j] = sum/(float)n1;
      }
  else if (dim==2)
    for(i=0;i<n1;i++)
      for(j=0;j<n3;j++){
	sum = 0.0;
	for(k=0;k<n2;k++)
	  sum += data[i][k][j];
	p[i][j] = sum/(float)n2;
      }
  else if (dim==3)
    for(i=0;i<n1;i++)
      for(j=0;j<n2;j++){
	sum = 0.0;
	for(k=0;k<n3;k++)
	  sum += data[i][j][k];
	p[i][j] = sum/(float)n3;
      }
  *rd1 = pn1;
  *rd2 = pn2;
  return p;
}
/**************************************-**************************************/
/*                                                                           */
/*                               PROJECT43                                   */
/*                                                                           */
/*  Project, i.e. average, a 4-d matrix onto a 3-d matrix w.r.t. a           */
/*  particular dimension.                                                    */
/*                                                                           */
/*****************************************************************************/
float ***project43(data,n1,n2,n3,n4,dim,rd1,rd2,rd3)
     float ****data;
     int n1,n2,n3,n4,dim;
     int *rd1,*rd2,*rd3; /* return dimensions of returned array */
{
  int i,j,k,l;
  float ***p,sum;
  int pn1,pn2,pn3,na;

  if (dim==1){
    na=n1; pn1=n2; pn2=n3; pn3=n4;
  }else if (dim==2){
    na=n2; pn1=n1; pn2=n3; pn3=n4;
  }else if (dim==3){
    na=n3; pn1=n1; pn2=n2; pn3=n4;
  }else if (dim==4){
    na=n4; pn1=n1; pn2=n2; pn3=n3;
  }else{
    pn1 = pn2 = pn3 = 0;
    exit_error("PROJECT43","dimension error");
  }
  
  p = (float ***)myalloc(pn1*sizeof(float **));
  for (i=0;i<pn1;i++){
    p[i] = (float **)myalloc(pn2*sizeof(float *));
    for (j=0;j<pn2;j++){
      p[i][j] = (float *)myalloc(pn3*sizeof(float));
    }
  }

  for(i=0;i<pn1;i++){
    for(j=0;j<pn2;j++){
      for(k=0;k<pn3;k++){
	sum = 0.0;
	for(l=0;l<na;l++){
	  if (dim==1)
	    sum += data[l][i][j][k];
	  else if (dim==2)
	    sum += data[i][l][j][k];
	  else if (dim==3)
	    sum += data[i][j][l][k];
	  else if (dim==4)
	    sum += data[i][j][k][l];
	}
	p[i][j][k] = sum/(float)na;
      }
    }
  }

  *rd1 = pn1;
  *rd2 = pn2;
  *rd3 = pn3;
  return p;
}
/**************************************-**************************************/
/*                                                                           */
/*                          WIDTH_AT_FRAC_HEIGHT                             */
/*                                                                           */
/*  More general than width at half height.                                  */
/*                                                                           */
/*  Find the maximum value between p1 and p2, inclusive.  Search from the    */
/*  maximum to the left for the first value less than or equal to frac*max.  */
/*  Do the same to the right.                                                */
/*                                                                           */
/*  Values returned:                                                         */
/*    "max" the maximum value.                                               */
/*    "imax" the position of the maximum.                                    */
/*    "left" the position of the left fractional height.                     */
/*    "right" the position of the right fractional height.                   */
/*                                                                           */
/*****************************************************************************/
void width_at_frac_height(data,n,p1,p2,frac,max,imax,left,right)
     float *data;
     int n,p1,p2;
     float frac,*max;
     int *imax,*left,*right;
{
  int i;
  int peakcnt;

  *imax = max_coord_farray(data+p1,p2-p1+1);
  *imax += p1;
  *max = data[*imax];

  /*** Adjust *imax in case peak is flat ***/
  i = *imax + 1;
  peakcnt = 0;
  while (i<n){
    if (data[i] == *max){
      peakcnt += 1;
      i += 1;
    }else
      i = n;
  }
  *imax += (int)(peakcnt/2);

  *left = -1; *right = -1;
  for(i=*imax;i>=0;i--)
    if ((data[i] <= frac*(*max))&&(*left == -1))
      *left = i;
  for(i=*imax;i<n;i++)
    if ((data[i] <= frac*(*max))&&(*right == -1))
      *right = i;
}
/**************************************-**************************************/
/*                                                                           */
/*                       WIDTH_AT_FRAC_HEIGHT_INTERP                         */
/*                                                                           */
/*  Extended 'width_at_frac_height' to interpolate.                          */
/*                                                                           */
/*  Find the maximum value between p1 and p2, inclusive.  Search from the    */
/*  maximum to the left for the first value less than or equal to frac*max.  */
/*  Do the same to the right.                                                */
/*                                                                           */
/*  Values returned:                                                         */
/*    "max" the maximum value.                                               */
/*    "imax" the position of the maximum.                                    */
/*    "left" the position of the left fractional height.                     */
/*    "right" the position of the right fractional height.                   */
/*                                                                           */
/*****************************************************************************/
void width_at_frac_height_interp(data,n,p1,p2,frac,max,imax,fleft,fright)
     float *data;
     int n,p1,p2;
     float frac,*max;
     int *imax;
     float *fleft,*fright;
{
  int i;
  int peakcnt,left,right;

  *imax = max_coord_farray(data+p1,p2-p1+1);
  *imax += p1;
  *max = data[*imax];

  /*** Adjust *imax in case peak is flat ***/
  i = *imax + 1;
  peakcnt = 0;
  while (i<n){
    if (data[i] == *max){
      peakcnt += 1;
      i += 1;
    }else
      i = n;
  }
  *imax += (int)(peakcnt/2);

  left = -1; right = -1;
  for(i=*imax;i>=0;i--)
    if ((data[i] <= frac*(*max))&&(left == -1))
      left = i;
  for(i=*imax;i<n;i++)
    if ((data[i] <= frac*(*max))&&(right == -1))
      right = i;

  if (left == -1)
    *fleft = -1.0;
  else{
    *fleft = (float)left;
    *fleft += (frac*(*max) - data[left])/ (data[left+1] - data[left]);
  }

  if (right == -1)
    *fright = -1.0;
  else{
    *fright = (float)(right - 1);
    *fright += (data[right-1] - frac*(*max))/ (data[right-1] - data[right]);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          FIND_TALLEST_MIDDLE_PEAK                         */
/*                                                                           */
/*  Find the tallest local peak in the range.  Neither end-point can be      */
/*  a peak, the peak must have lower values to the left and right.           */
/*  Return the index for the peak, or -1 if no such peak exists.             */
/*                                                                           */
/*****************************************************************************/
int find_tallest_middle_peak(data,n,p1,p2)
     float *data;
     int n,p1,p2;
{
  int i;
  int imax;

  imax = -1;
  for(i=(p1+1);i<p2;i++)
    if ((data[i] > data[i-1]) && (data[i] > data[i+1])){
      if (imax == -1)
	imax = i;
      else if (data[imax] < data[i])
	imax = i;
    }

  return imax;
}
/**************************************-**************************************/
/*                              Blake                                        */
/*                          FIND_NL_BIN_INDEX                                */
/*                                                                           */
/* Find the indices of the closest filter level across a time varying signal.*/
/* For use in the non-linear filtering step.                                 */
/*                                                                           */
/*****************************************************************************/
void find_nl_bin_index(indices,numfilt,lengthx,x,tres,tau,sums)
     int *indices;
     int numfilt,lengthx;
     float *x;
     float tres,tau;
     float *sums;
{
  int i,j,level;
  float min,sum;	
  
  // get the minimum of x
  min = min_of_farray(x,lengthx);
  
  // get the indices
  for(i=0; i<lengthx; i++)
    {
      sum = 0;
      for(j=0; j<=i; j++) {
	sum += (1/pow(tau,j))*x[i-j]*x[i-j];
      }
      if(sum < 0.0) {sum = 0;}
      sums[i] = sum;
      level = my_rint(sum/tres);
      if(level > numfilt-1) indices[i] = numfilt-1;
      else if(level < 0) indices[i] = 0;
      else indices[i] = level;
      //EDITif(level > numfilt-1) indices[i] = 0;
      //else if(level < 0) indices[i] = numfilt-1;
      //else indices[i] = abs(level - (num_filt - 1));
    }
}
/**************************************-**************************************/
/*                                                                           */
/*                       WIDTH_AT_FRAC_HEIGHT_CHECK_DOUBLE                   */
/*                                                                           */
/*  A modified version of "width_at_frac_height".                            */
/*  This checks for a "double" peak (for STA analysis) and adjusts the       */
/*  'left' and 'right' values accordingly.                                   */
/*                                                                           */
/*  A secondary double peak is one which lies within 'dt' of the original    */
/*  peak and achieves 'pf' fraction of the original peaks height.            */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - If '*rimin' is not -1, there was a double peak.                        */
/*                                                                           */
/*****************************************************************************/
void width_at_frac_height_check_double(data,n,p1,p2,dt,frac,pf,rmax,rimax,
				       rleft,rright,rimax2,rimin)
     float *data;
     int n,p1,p2,dt;
     float frac,pf,*rmax;
     int *rimax,*rleft,*rright,*rimax2,*rimin;
{
  int i;
  int imax,imax2,left,right,imin,iminl,iminr,ipr,ipl,tl,tr,ttl,ttr;
  float max;

  if ((p1 < 0)||(p2 >= n))
    exit_error("WIDTH_AT_FRAC_HEIGHT_CHECK_DOUBLE","Invalid range");

  imax = max_coord_farray(data+p1,p2-p1+1);
  imax += p1;
  max = data[imax];

  /*** Get left and right endpoints for double peak search. ***/
  tl = imax - dt;
  if (tl < 0)
    tl = 0;
  tr = imax + dt;
  if (tr > p2)
    tr = p2;

  /*** Search for secondary peaks to the left and right. ***/
  ipl = find_tallest_middle_peak(data,n,tl,imax-1);
  ipr = find_tallest_middle_peak(data,n,imax+1,tr);

  /*** If there are secondary peaks that are tall enough, adjust. ***/
  ttl = ttr = imax;
  imin = imax2 = -1;
  if (ipl != -1){
    if (data[ipl] >= (pf * max)){
      ttl = ipl;
      /*printf("DOUBLE PEAK TO LEFT\n");*/
      iminl = min_coord_farray(data+ipl,imax-ipl+1);
      if (iminl != -1)
	iminl += ipl;
      imin = iminl;
      imax2 = ipl;
    }
  }
  if (ipr != -1){
    if (data[ipr] >= (pf * max)){
      ttr = ipr;
      /*printf("DOUBLE PEAK TO RIGHT\n");*/
      iminr = min_coord_farray(data+imax,ipr-imax+1);
      if (iminr != -1)
	iminr += imax; /* Express index relative to data[0]. */

      if (imax2 == -1){
	imin = iminr;
	imax2 = ipr;
      }else if (data[ipr] > data[imax2]){
	imin = iminr;
	imax2 = ipr;
      }
    }
  }

  left = right = -1;
  for(i=ttl;i>=0;i--)
    if ((data[i] <= frac*max)&&(left == -1))
      left = i;
  for(i=ttr;i<n;i++)
    if ((data[i] <= frac*max)&&(right == -1))
      right = i;
  
  *rmax = max; *rimax = imax; *rleft = left; *rright = right; *rimin = imin;
  *rimax2 = imax2;
}
/**************************************-**************************************/
/*                                                                           */
/*                            FIRST_LEFT_FRAC_RISE                           */
/*                                                                           */
/*  Find the first point, searching left to right, that equals or exceeds    */
/*  "frac" of the peak value.                                                */
/*                                                                           */
/*  Values returned:                                                         */
/*    "max" the maximum value.                                               */
/*    "imax" the position of the maximum.                                    */
/*    "left" the position of the left fractional height.                     */
/*                                                                           */
/*****************************************************************************/
void first_left_frac_rise(data,n,p1,p2,frac,max,imax,left)
     float *data;
     int n,p1,p2;
     float frac,*max;
     int *imax,*left;
{
  int i;

  *imax = max_coord_farray(data+p1,p2-p1+1);
  *imax += p1;
  *max = data[*imax];
  *left = -1;
  for(i=0;i<=*imax;i++)
    if ((data[i] >= frac*(*max))&&(*left == -1))
      *left = i;
}
/**************************************-**************************************/
/*                                                                           */
/*                          FIRST_LEFT_FRAC_RISE_INTERP                      */
/*                                                                           */
/*  Find the first point, searching left to right, that equals or exceeds    */
/*  "frac" of the peak value.                                                */
/*                                                                           */
/*  Values returned:                                                         */
/*    "max" the maximum value.                                               */
/*    "imax" the position of the maximum.                                    */
/*    "left" the position of the left fractional height.                     */
/*                                                                           */
/*  *** WARNING: WYETH - starts search for left rise at 'p1', UNLIKE         */
/*      'first_left_frac_rise'.                                              */
/*                                                                           */
/*****************************************************************************/
void first_left_frac_rise_interp(data,n,p1,p2,frac,max,imax,left,fleft)
     float *data;
     int n,p1,p2;
     float frac,*max;
     int *imax,*left;
     float *fleft;
{
  int i;
  
  *imax = max_coord_farray(data+p1,p2-p1+1);
  *imax += p1;
  *max = data[*imax];
  *left = -1;
  for(i=p1;i<=*imax;i++)
    if ((data[i] >= frac*(*max))&&(*left == -1)){
      *left = i;
    }

  if (*left == -1)
    *fleft = -1.0;
  else if (*left == 0){
    if (data[0] == frac*(*max))
      *fleft = 0.0;
    else{
      *left = -1;
      *fleft = -1.0;
    }
  }else{
    *fleft = (float)(*left - 1);
    *fleft += (frac*(*max) - data[*left-1])/ (data[*left] - data[*left-1]);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           LAST_LEFT_FRAC_RISE_INTERP                      */
/*                                                                           */
/*  Find the first point, searching right to left, that equals 'frac' of     */
/*  the peak value.                                                          */
/*                                                                           */
/*  Values returned:                                                         */
/*    "max" the maximum value.                                               */
/*    "imax" the position of the maximum.                                    */
/*    "left" the position of the left fractional height.                     */
/*                                                                           */
/*  *** WARNING: WYETH - starts search for left rise at 'p1', UNLIKE         */
/*      'first_left_frac_rise'.                                              */
/*                                                                           */
/*****************************************************************************/
void last_left_frac_rise_interp(data,n,p1,p2,frac,max,imax,left,fleft)
     float *data;
     int n,p1,p2;
     float frac,*max;
     int *imax,*left;
     float *fleft;
{
  int i,k;
  
  k = p1 + max_coord_farray(data+p1,p2-p1+1);
  *imax = k;
  *max = data[k];
  *left = -1;

  /*** Search from max to the left ***/
  for(i=k;i>=0;i--)
    if ((data[i] <= frac*(*max))&&(*left == -1)){
      /*printf("data[%d] = %f\n",i,data[i]);*/
      *left = i;
    }

  if (*left == -1)
    *fleft = -1.0;
  else{
    *fleft = (float)*left;
    *fleft += (frac*(*max) - data[*left])/ (data[*left+1] - data[*left]);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                       WIDTH_AT_FRAC_HEIGHT_CIRC_FARRAY                    */
/*                                                                           */
/*  Like "width_at_frac_height" but assumes the farray wraps around, so the  */
/*  search for left or right rise points might have to wrap around.          */
/*                                                                           */
/*  Find the maximum value between p1 and p2, inclusive.  Search from the    */
/*  maximum to the left for the first value less than or equal to frac*max.  */
/*  Do the same to the right.                                                */
/*                                                                           */
/*  Values returned:                                                         */
/*    "max" the maximum value.                                               */
/*    "imax" the position of the maximum.                                    */
/*    "left" the position of the left fractional height.                     */
/*    "right" the position of the right fractional height.                   */
/*                                                                           */
/*****************************************************************************/
void width_at_frac_height_circ_farray(data,n,p1,p2,frac,max,imax,left,right)
     float *data;
     int n,p1,p2;
     float frac,*max;
     int *imax,*left,*right;
{
  int i;

  *imax = max_coord_farray(data+p1,p2-p1+1);
  *imax += p1;
  *max = data[*imax];

  *left = -1;
  for(i=*imax;i>=0;i--)
    if ((data[i] <= frac*(*max))&&(*left == -1))
      *left = i;
  if (*left == -1) // Wrap around to upper part of array
    for(i=n-1;i>*imax;i--)
      if ((data[i] <= frac*(*max))&&(*left == -1))
	*left = i;

  *right = -1;
  for(i=*imax;i<n;i++)
    if ((data[i] <= frac*(*max))&&(*right == -1))
      *right = i;
  if (*right == -1) // Wrap around to lower part of array
    for(i=0;i<*imax;i++)
      if ((data[i] <= frac*(*max))&&(*right == -1))
	*right = i;
}
/**************************************-**************************************/
/*                                                                           */
/*                            WIDTH_AT_FRAC_AREA                             */
/*                                                                           */
/*  Determine the width (and locations) for which "frac" fraction of the     */
/*  area of the peak is included.  Maximum area is determined by searching   */
/*  from the specified index "p" until the area "crit" away from the         */
/*  current point becomes less than the area up to the current point.  This  */
/*  is done in both directions, and the maximum area in that range is taken. */
/*                                                                           */
/*  Values returned:                                                         */
/*    "amax" the maximum area.                                               */
/*    "left" the position of the left fractional area.                       */
/*    "right" the position of the right fractional area.                     */
/*                                                                           */
/*  *** WYETH - I think this requires the data to drop below zero.           */
/*                                                                           */
/*****************************************************************************/
void width_at_frac_area(data,n,p,crit,frac,amax,left,right)
     float *data;
     int n,p,crit;
     float frac,*amax;
     int *left,*right;
{
  int i;
  float *total,max,maxleft,maxright,eps;

  total = (float *)myalloc(n*sizeof(float));
  total[p] = data[p]/2.0;
  max = maxright = -1.0;
  for(i=p+1;i<n;i++){ /*** TOTAL RIGHT ***/
    total[i] = total[i-1] + data[i];
    if (total[i] > max)
      max = total[i];
    if ((i >= p+crit)&&(maxright < 0.0))
      if (total[i] < total[i-crit])
	maxright = max;
  }
  max = maxleft = -1.0;
  for(i=p-1;i>=0;i--){ /*** TOTAL LEFT ***/
    total[i] = total[i+1] + data[i];
    if (total[i] > max)
      max = total[i];
    if ((i <= p-crit)&&(maxleft < 0.0))
      if (total[i] < total[i+crit])
	maxleft = max;
  }
  *amax = maxleft + maxright;

  eps = 0.5 * (*amax * (1.0-frac));

  *left = p;
  while(total[*left] < (maxleft-eps))
    *left -= 1;
  *right = p;
  while(total[*right] < (maxright-eps))
    *right += 1;

  myfree(total);
}
/**************************************-**************************************/
/*                                                                           */
/*                            ESTIMATE_AREA_PEAK                             */
/*                                                                           */
/*  Determine the area of the peak by searching away from index "p" until    */
/*  the area increase in the last "crit" steps is less than "frac" of the    */
/*  area measured so far.                                                    */
/*                                                                           */
/*  Values returned:                                                         */
/*    "rarea" the maximum area.                                              */
/*    "rleft" starting position for measuring area.                          */
/*    "rright" end position for measuring area.                              */
/*                                                                           */
/*****************************************************************************/
void estimate_area_peak(data,n,p,crit,frac,rarea,rleft,rright)
     float *data;
     int n,p,crit;
     float frac,*rarea;
     int *rleft,*rright;
{
  int i;
  int left,right;
  float *total;

  total = (float *)myalloc(n*sizeof(float));
  total[p] = data[p]/2.0;
  left = right = -1;
  for(i=p+1;i<n;i++){ // TOTAL RIGHT
    total[i] = total[i-1] + data[i];
    if ((i >= p+crit)&&(right < 0.0))
      if ((total[i] - total[i-crit])/total[i-crit] < frac)
	right = i;
  }
  for(i=p-1;i>=0;i--){ // TOTAL LEFT
    total[i] = total[i+1] + data[i];
    if ((i <= p-crit)&&(left < 0))
      if ((total[i] - total[i+crit])/total[i+crit] < frac)
	left = i;
  }
  *rarea = total[left] + total[right];
  *rleft = left; *rright = right;

  myfree(total);
}
/**************************************-**************************************/
/*                                                                           */
/*                       FIRST_LEFT_MIN_BELOW_CRIT_FARRAY                    */
/*                                                                           */
/*  Find the first minimum point that is less than "crit" and is to the      */
/*  left (lower index) side of "p1".                                         */
/*                                                                           */
/*****************************************************************************/
int first_left_min_below_crit_farray(data,n,p1,crit)
     float *data;
     int n,p1;
     float crit;
{
  int i,k;
  int done;

  i = -1;
  k = p1-1;
  done = 0;
  while(!done && (i==-1)){
    if (k < 1)
      done = 1;
    else if ((data[k] < crit)&&(data[k] <= data[k+1])&&(data[k] <= data[k-1]))
      i = k;
    else
      k -= 1;
  }
  return i;
}
/**************************************-**************************************/
/*                                                                           */
/*                               FIND_PEAK_START                             */
/*                                                                           */
/*  Attempt to estimate the position to the left of the peak where the peak  */
/*  begins.  This is used for computing time (or latency) of beginning of    */
/*  response.                                                                */
/*                                                                           */
/*****************************************************************************/
int find_peak_start(data,n,p1,p2,frac)
     float *data;
     int n,p1,p2;
     float frac;
{
  int k;
  int imax,left,right,width,d2max;
  float max,*d2;

  width_at_frac_height(data,n,p1,p2,frac,&max,&imax,&left,&right);
  width = right-left;
  d2 = second_derivative_farray(data,n);
  k = left-width;

  if ((k > 0)&&(k < n-width))
    d2max = k + max_coord_farray(d2+k,width);
  else
    d2max = -1;
  myfree(d2);

  return d2max;
}
/**************************************-**************************************/
/*                                                                           */
/*                            FIRST_ZERO_BEFORE_PEAK                         */
/*                                                                           */
/*  Find the first zero value searching leftward from the peak.              */
/*                                                                           */
/*****************************************************************************/
int first_zero_before_peak(data,n,p1,p2)
     float *data;
     int n,p1,p2;
{
  int i,k;

  i = p1 + max_coord_farray(data+p1,p2-p1+1);
  k = -1;
  while((i>=p1)&&(k==-1)){
    if (data[i] <= 0.0)
      k = i;
    i -= 1;
  }
  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                         FIRST_FRAC_BEFORE_PEAK_INTERP                     */
/*                                                                           */
/*  Find the first value searching leftward from the peak that is 'frac'     */
/*  fraction of the peak value, interpolate.                                 */
/*                                                                           */
/*****************************************************************************/
float first_frac_before_peak_interp(data,n,frac,p1,p2)
     float *data;
     int n;
     float frac;
     int p1,p2;
{
  int i,k;
  float cf,ff;

  i = p1 + max_coord_farray(data+p1,p2-p1+1);
  cf = frac * data[i];
  k = -1;
  while((i>=p1)&&(k==-1)){
    if (data[i] <= cf)
      k = i;
    i -= 1;
  }
  ff = (float)k + (cf - data[k])/(data[k+1] - data[k]);
  return ff;
}
/**************************************-**************************************/
/*                                                                           */
/*                              GET_ENTROPY_FARRAY                           */
/*                                                                           */
/*  Return the entropy in bits, interpreting the farray as a set of          */
/*  probabilities.                                                           */
/*                                                                           */
/*****************************************************************************/
float get_entropy_farray(data,n)
     float *data;
     int n;
{
  int i;
  float h;

  h = 0.0;
  for(i=0;i<n;i++){
    if (data[i] > 0.0)
      h += data[i] * (float)my_log2(1.0/(double)data[i]);
  }
  return h;
}
/**************************************-**************************************/
/*                                                                           */
/*                            COVAR_MATRIX_2D_FARRAY                         */
/*                                                                           */
/*  The covariance matrix is returned as c[m][m].  The input                 */
/*  matrix "r" is modified by having its mean subtracted.  The mean          */
/*  is returned in "rmean".                                                  */
/*                                                                           */
/*****************************************************************************/
float **covar_matrix_2d_farray(r,n,m,rmean)
     float **r;
     int n,m;
     float **rmean;
{
  int i,j,k;
  float **c;
  float *mean; // mean response
  float *t,sdev;

  printf("  COVAR_MATRIX_2D_FARRAY\n");
  
  mean = (float *)myalloc(m*sizeof(float));
  t = (float *)myalloc(n*sizeof(float));
  for (i=0;i<m;i++){
    for (j=0;j<n;j++)
      t[j] = r[j][i];
    mean_sdev_farray(t,n,&mean[i],&sdev);
  }
  /***write_farray_plot("zz.mr",mean,m);***/

  c = (float **)myalloc(m*sizeof(float *));
  for (i=0;i<m;i++)
    c[i] = get_zero_farray(m);

  for (i=0;i<n;i++) /* SUBTRACT MEAN FROM ALL RESPONSES */
    for (j=0;j<m;j++)
      r[i][j] -= mean[j];
  for (i=0;i<n;i++) /* COMPUTE HALF OF COVARIANCE MATRIX */
    for (j=0;j<m;j++)
      for (k=0;k<=j;k++)
	c[j][k] += r[i][k] * r[i][j];
  for (j=0;j<m;j++)
    for (k=0;k<=j;k++){
      c[j][k] /= (float)n;
      c[k][j] = c[j][k]; /* MATRIX IS SYMMETRICAL */
    }
  myfree(t);
  *rmean = mean;
  return c;
}
/**************************************-**************************************/
/*                                                                           */
/*                      GET_EIGENVECTORS_COVAR_FARRAY                        */
/*                                                                           */
/*  Return the eigenvectors as the columns of the matrix "z" and the         */
/*  eigenvalues in "d", both are sorted.  The input matrix "c" is            */
/*  "m" by "m".                                                              */
/*                                                                           */
/*****************************************************************************/
void get_eigenvectors_covar_farray(c,m,z,d,jacobi_flag)
     float **c;
     int m;
     float **z,*d;
     int jacobi_flag;
{
  int i,j;
  int nrot;
  float **a,**v,*e;

  printf("  GET_EIGENVECTORS_COVAR_FARRAY\n");

  a = matrix(1,m,1,m); /* NumRec Indexing [1..] */
  for (i=0;i<m;i++)
    for (j=0;j<m;j++)
      a[i+1][j+1] = c[i][j];

  if (jacobi_flag){ /*** JACOBI METHOD ***/
    printf("    Jacobi method\n");
    v = matrix(1,m,1,m); /* NumRec Indexing [1..] */
    jacobi(a,m,d-1,v,&nrot);
    printf("    %d Jacobi rotations performed.\n",nrot);
    eigsrt(d-1,v,m);
    free_matrix(a,1,m,1,m);
  }else{ /*** TRIDIAGONAL QL IMPLICIT METHOD ***/
    printf("    Tridiagonal QL Implicit method\n");
    e = (float *)myalloc(m*sizeof(float));
    tred2(a,m,d-1,e-1);
    tqli(d-1,e-1,m,a);
    myfree(e);
    eigsrt(d-1,a,m);
    v = a;
  }
  for(i=0;i<m;i++) /* Switch back from NumRec indexing */
    for(j=0;j<m;j++)
      z[i][j] = v[i+1][j+1];
  free_matrix(v,1,m,1,m);
}
/**************************************-**************************************/
/*                                                                           */
/*                                MATRIX_PRINT                               */
/*                                                                           */
/*****************************************************************************/
void matrix_print(a,m,n,dp)
     float **a;
     int m,n;
     int dp;      // Decimal precision
{
  int i,j;

  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      if (dp == 0)
	printf(" %d",(int)a[i][j]);
      else if (dp == 1)
	printf(" %.1f",a[i][j]);
      else if (dp == 2)
	printf(" %.2f",a[i][j]);
      else if (dp == 3)
	printf(" %.3f",a[i][j]);
      else if (dp == 4)
	printf(" %.4f",a[i][j]);
      else
	printf(" %f",a[i][j]);
    }
    printf("\n");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MATRIX_TRANSPOSE                             */
/*                                                                           */
/*****************************************************************************/
float **matrix_transpose(a,m,n)
     float **a;
     int m,n;
{
  int i,j;
  float **mt;

  mt = get_2d_farray(n,m);

  for(i=0;i<m;i++)
    for(j=0;j<n;j++)
      mt[j][i] = a[i][j];

  return mt;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MATRIX_TIMES_VECTOR                           */
/*                                                                           */
/*  Multiply the matrix 'a' times column vector 'x'                          */
/*                                                                           */
/*****************************************************************************/
float *matrix_times_vector(a,x,m,n)
     float **a,*x;   //  a[m][n], x[x]
     int m,n;        // 'a' is m x n, 'x' is n
{
  int i;
  float *y;

  y = (float *)myalloc(n*sizeof(float));

  for(i=0;i<m;i++)
    y[i] = multiply_integrate_farrays(a[i],x,n);

  return y;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MATRIX_MULTIPLY                             */
/*                                                                           */
/*  Multiply the matrix 'a' times 'b'.                                       */
/*                                                                           */
/*****************************************************************************/
float **matrix_multiply(a,b,m,s,n)
     float **a,**b;  // Multiply 'a' x 'b'
     int m,s,n;      // 'a' is m x s, 'b' is s by n
{
  int i,j,k;
  float **mn;

  mn = get_zero_2d_farray(m,n);

  for(i=0;i<m;i++)
    for(j=0;j<n;j++){
      for(k=0;k<s;k++)
	mn[i][j] += a[i][k] * b[k][j];
    }

  return mn;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MATRIX_INVERT_3x3                            */
/*                                                                           */
/*  Return the inverse of the 3x3 matrix.  The first index is ROW index.     */
/*                                                                           */
/*****************************************************************************/
float **matrix_invert_3x3(m,rd)
     float **m;  // Matrix, first index is ROW index.
     float *rd;  // return the determinant; if 0, singlular (no determinant)
{
  float a,b,c,d,e,f,g,h,i;
  float deta,**mi,**mt;

  //       a b c
  //  A =  d e f   ,   det A = aei - afh - bdi + bfg + cdh - ceg
  //       g h i

  a = m[0][0];  b = m[0][1];  c = m[0][2];
  d = m[1][0];  e = m[1][1];  f = m[1][2];
  g = m[2][0];  h = m[2][1];  i = m[2][2];

  // Compute the determinant
  deta = a*e*i - a*f*h - b*d*i + b*f*g + c*d*h - c*e*g;

  mi = get_2d_farray(3,3);

  mi[0][0] =  (e*i - f*h)/deta;
  mi[0][1] = -(d*i - f*g)/deta;
  mi[0][2] =  (d*h - e*g)/deta;

  mi[1][0] = -(b*i - c*h)/deta;
  mi[1][1] =  (a*i - c*g)/deta;
  mi[1][2] = -(a*h - b*g)/deta;

  mi[2][0] =  (b*f - c*e)/deta;
  mi[2][1] = -(a*f - c*d)/deta;
  mi[2][2] =  (a*e - b*d)/deta;

  mt = matrix_transpose(mi,3,3);
  
  free_2d_farray(mi,3);

  *rd = deta;
  return mt;
}
/**************************************-**************************************/
/*                                                                           */
/*                            PRINT_MIN_MAX_FARRAY                           */
/*                                                                           */
/*  For debugging - print out stats about the farray.                        */
/*                                                                           */
/*****************************************************************************/
void print_min_max_farray(fa,n,name)
     float *fa;
     int n;
     char *name;
{
  int i;
  int n_nan;
  float min,max,mu,sd;

  printf("===> ARRAY:  %s\n",name);

  //
  //  Check for NaN values
  //
  n_nan = 0;
  for(i=0;i<n;i++){
    if (isnan(fa[i])){
      printf("  Element %d is NaN\n",i);
      n_nan += 1;
      if (n_nan > 10){
	printf(" ...  exiting\n");
	exit(0);
      }
    }
  }

  get_min_max_farray(fa,n,&min,&max);

  printf("  Length %d  min,max  %f %f\n",n,min,max);

  mean_sdev_farray(fa,n,&mu,&sd);
  printf("  Mean (SD)  %f (%f)\n",mu,sd);

}
/**************************************-**************************************/
/*                                                                           */
/*                           PRINT_MIN_MAX_2D_FARRAY                         */
/*                                                                           */
/*****************************************************************************/
void print_min_max_2d_farray(fa,xn,yn,name)
     float **fa;
     int xn,yn;
     char *name;
{
  int i,j;
  int n_nan;
  float min,max,mu,sd;

  if (name != NULL){
    printf("  PRINT_MIN_MAX_2D_FARRAY\n");
    printf("    Array:  %s\n",name);
  }

  //
  //  Check for NaN values
  //
  n_nan = 0;
  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++)
      if (isnan(fa[i][j])){
	printf("  Element %d %d  is NaN\n",i,j);
	n_nan += 1;
	if (n_nan > 10){
	  printf(" ...  exiting\n");
	  exit(0);
	}
      }
  
  get_min_max_2d_farray(fa,xn,yn,&min,&max);
  single_mean_sdev_2d_farray(fa,xn,yn,&mu,&sd);

  printf("  %d x %d  min,max =  %f %f   mu,sd = %f %f\n",xn,yn,min,max,mu,sd);
}
/**************************************-**************************************/
/*                                                                           */
/*                           PRINT_MIN_MAX_3D_FARRAY                         */
/*                                                                           */
/*****************************************************************************/
void print_min_max_3d_farray(fa,x0,xn,y0,yn,z0,zn,name)
     float ***fa;
     int x0,xn,y0,yn,z0,zn;
     char *name;
{
  int i,j,k;
  int n_nan;
  float min,max;

  printf("===> ARRAY:  %s\n",name);

  //
  //  Check for NaN values
  //
  n_nan = 0;
  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++)
      for(k=0;k<zn;k++)
	if (isnan(fa[i+x0][j+y0][k+z0])){
	  printf("  Element %d %d %d  is NaN\n",i,j,k);
	  n_nan += 1;
	  if (n_nan > 10){
	    printf(" ...  exiting\n");
	    exit(0);
	  }
	}

  get_min_max_3d_farray(fa,x0,xn,y0,yn,z0,zn,&min,&max);

  printf("  [%d,%d,%d]  min,max  %f %f\n",xn,yn,zn,min,max);
}
/**************************************-**************************************/
/*                                                                           */
/*                             DISTRIB_BOX_MATCH                             */
/*                                                                           */
/*  Given a set of box integrals 'x', find the upper limits to use for       */
/*  each box to make a piece-wise constant approximation to a piece-wise     */
/*  constant distribution defined by 'd'.                                    */
/*                                                                           */
/*                                                                           */
/*        222                                                                */
/*    1111222333333      <-- The area of this distrib must be 1.0            */
/*    111122233333344                                                        */
/*    1111222333333445                                                       */
/*       ^  ^     ^ ^^   <-- These show the chosen upper limits for the      */
/*    0..............1       boxes, 1, 2, 3, 4 and 5, which had weights      */
/*                           of 12, 12, 18, 4 and 1, respectively, to        */
/*                           approximate a distribution that looks like      */
/*                           the resulting histogram shown.                  */
/*                                                                           */
/*****************************************************************************/
float *distrib_box_match(d,dn,d0,dx,w,wn,dumpfile)
     float *d;  // [dn] Points that define a distrib.
     int dn;    // Number of points in 'd'
     float d0;  // First 'd' value applies from d0 to d0+dd
     float dx;  // Width on x-axis for each 'd' value
     float *w;  // [wn] Weights for 'wn' boxes to use to approximate 'd'
     int wn;    // Number of boxes
     char *dumpfile; // File to write plots
{
  int i,k;
  int done,npl;
  float *rlim,*dnorm,*d_cdf,*wnorm,targ,extra,xx,wid,prev,*px,*py;

  dnorm = copy_farray(d,dn);
  norm_area_farray(dnorm,dn,1.0);

  d_cdf = get_cumulative_farray(dnorm,dn);  // d_cdf[0] will have dnorm[0]

  wnorm = copy_farray(w,wn);
  norm_area_farray(wnorm,wn,1.0);

  rlim = (float *)myalloc(wn*sizeof(float));

  k = 0;  // Highest index into 'd' with less area than accumulated so far
  targ = 0.0;
  for(i=0;i<(wn-1);i++){
    targ += wnorm[i];

    done = 0;
    while(done == 0){
      if (d_cdf[k] < targ){  // d_cdf[0] has "sum" of first value
	k += 1;
	if (k >= dn)
	  done = 1;
      }else
	done = 1;
    }

    if (k == 0)  // Because d_cdf[0] contains "sum" of first item
      extra = targ;  // Cumulative amount is 0.
    else
      extra = targ - d_cdf[k-1];
    xx = dx * extra / dnorm[k];

    rlim[i] = d0 + k*dx + xx;

    //printf("i %d  targ %f  k %d  extra %f  xval %f  xx %f\n",i,
    //targ,k,extra,rlim[i],xx);
  }
  rlim[wn-1] = d0 + dn*dx;

  if (dumpfile != NULL){
    if (strcmp(dumpfile,"NULL")!=0){
      //
      //  Plot the distrib data to be matched.
      //
      npl = 1 + dn*3;
      px = get_farray(npl);
      py = get_farray(npl);
      px[0] = d0;
      py[0] = 0.0;
      for(i=0;i<dn;i++){
	px[i*3+1] = d0 + dx * i;
	px[i*3+2] = d0 + dx * (i+1);
	px[i*3+3] = d0 + dx * (i+1);
	py[i*3+1] = d[i];
	py[i*3+2] = d[i];
	py[i*3+3] = 0.0;
      }
      append_farray_xy_plot(dumpfile,px,py,npl,"raw_distrib");
      myfree(px);
      myfree(py);

      //
      //  Same thing, normalized
      //
      npl = 1 + dn*3;
      px = get_farray(npl);
      py = get_farray(npl);
      px[0] = d0;
      py[0] = 0.0;
      for(i=0;i<dn;i++){
	px[i*3+1] = d0 + dx * i;
	px[i*3+2] = d0 + dx * (i+1);
	px[i*3+3] = d0 + dx * (i+1);
	py[i*3+1] = dnorm[i];
	py[i*3+2] = dnorm[i];
	py[i*3+3] = 0.0;
      }
      append_farray_xy_plot(dumpfile,px,py,npl,"target");
      myfree(px);
      myfree(py);

      //
      //  Plot the resulting match
      //
      npl = 1 + wn*3;
      px = get_farray(npl);
      py = get_farray(npl);
      px[0] = d0;
      py[0] = 0.0;
      prev = d0;
      for(i=0;i<wn;i++){
	px[i*3+1] = px[i*3];  // Same as previous
	px[i*3+2] = rlim[i];
	px[i*3+3] = rlim[i];

	wid = (rlim[i] - prev)/dx;
	py[i*3+1] = wnorm[i]/wid;
	py[i*3+2] = wnorm[i]/wid;
	py[i*3+3] = 0.0;

	prev = rlim[i];
      }
      append_farray_xy_plot(dumpfile,px,py,npl,"match");
      myfree(px);
      myfree(py);
    }
  }

  myfree(dnorm);
  myfree(wnorm);
  myfree(d_cdf);

  return rlim;
}
