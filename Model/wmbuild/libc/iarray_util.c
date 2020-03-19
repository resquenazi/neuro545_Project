/*****************************************************************************/
/*                                                                           */
/*  iarray_util.c                                                            */
/*  wyeth bair                                                               */
/*  caltech                                                                  */
/*  03/19/93                                                                 */
/*                                                                           */
/*  All routine names must have "iarray" in the title.                       */
/*                                                                           */
/*  Types of iarray data structures:                                         */
/*                                                                           */
/*      iarray - 1D                                                          */
/*   2d_iarray - 2D rectangular                                              */
/*  2dv_iarray - 2D variable length in 2nd dimension, specified by *cnt      */
/*   3d_iarray - 3D rectangular                                              */
/*  3dv_iarray - 3D variable length in 3rd dimension, specified by **cnt     */
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

/**************************************-**************************************/
/*                                                                           */
/*                                ZERO_IARRAY                                */
/*                                                                           */
/*****************************************************************************/
void zero_iarray(data,n)
     int *data,n;
{
  int i;

  for(i=0;i<n;i++)
    data[i] = 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_CONST_2D_IARRAY                           */
/*                                                                           */
/*****************************************************************************/
int **get_const_2d_iarray(m,n,c)
     int m,n;
     int c;
{
  int i,j;
  int **data;
  
  data = (int **)myalloc(m*sizeof(int *));
  for(i=0;i<m;i++){
    data[i] = (int *)myalloc(n*sizeof(int));
    for(j=0;j<n;j++)
      data[i][j] = c;
  }
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_ZERO_2D_IARRAY                            */
/*                                                                           */
/*****************************************************************************/
int **get_zero_2d_iarray(m,n)
     int m,n;
{
  return get_const_2d_iarray(m,n,0);
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_2D_IARRAY                               */
/*                                                                           */
/*****************************************************************************/
int **get_2d_iarray(m,n)
     int m,n;
{
  int i;
  int **data;
  
  data = (int **)myalloc(m*sizeof(int *));
  for(i=0;i<m;i++)
    data[i] = (int *)myalloc(n*sizeof(int));
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_2D_POINTER_IARRAY                          */
/*                                                                           */
/*****************************************************************************/
int ***get_2d_pointer_iarray(m,n)
     int m,n;
{
  int i;
  int ***a;

  a = (int ***)myalloc(m*sizeof(int **));
  for(i=0;i<m;i++)
    a[i] = (int **)myalloc(n*sizeof(int *));
  
  return a;
}
/**************************************-**************************************/
/*                                                                           */
/*                               FREE_2D_IARRAY                              */
/*                                                                           */
/*****************************************************************************/
void free_2d_iarray(data,n)
     int **data;
     int n;
{
  int i;
  
  for(i=0;i<n;i++)
    myfree(data[i]);
  myfree(data);
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_3D_IARRAY                               */
/*                                                                           */
/*****************************************************************************/
int ***get_3d_iarray(d1,d2,d3)
     int d1,d2,d3;
{
  int i1,i2;
  int ***data;
  
  data = (int ***)myalloc(d1*sizeof(int **));
  if (d2 >= 0)
    for(i1=0;i1<d1;i1++){
      data[i1] = (int **)myalloc(d2*sizeof(int *));
      if (d3 >= 0)
	for(i2=0;i2<d2;i2++)
	  data[i1][i2] = (int *)myalloc(d3*sizeof(int));
    }
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_ZERO_3D_IARRAY                            */
/*                                                                           */
/*****************************************************************************/
int ***get_zero_3d_iarray(d1,d2,d3)
     int d1,d2,d3;
{
  int i1;
  int ***data;
  
  data = (int ***)myalloc(d1*sizeof(int **));
  if (d2 >= 0)
    for(i1=0;i1<d1;i1++)
      data[i1] = get_zero_2d_iarray(d2,d3);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                               FREE_3D_IARRAY                              */
/*                                                                           */
/*****************************************************************************/
void free_3d_iarray(data,d1,d2,d3)
     int ***data;
     int d1,d2,d3;
{
  int i1,i2;
  
  for(i1=0;i1<d1;i1++){
    for(i2=0;i2<d2;i2++)
      myfree(data[i1][i2]);
    myfree(data[i1]);
  }
  myfree(data);
}
/**************************************-**************************************/
/*                                                                           */
/*                                COPY_IARRAY                                */
/*                                                                           */
/*  Return a copy of an integer array.                                       */
/*                                                                           */
/*****************************************************************************/
int *copy_iarray(data,n)
     int *data,n;
{
  int i;
  int *cdata;
  
  cdata = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++)
    cdata[i] = data[i];
  return cdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                               COPY_2D_IARRAY                              */
/*                                                                           */
/*****************************************************************************/
int **copy_2d_iarray(data,xn,yn)
     int **data;  // [xn][yn]
     int xn,yn;
{
  int i;
  int **cdata;
  
  cdata = (int **)myalloc(xn*sizeof(int *));
  for(i=0;i<xn;i++)
    cdata[i] = copy_iarray(data[i],yn);

  return cdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                               APPEND_IARRAY                               */
/*                                                                           */
/*  Append the second iarray to the first.                                   */
/*                                                                           */
/*****************************************************************************/
void append_iarray(pd1,n1,d2,n2)
     int **pd1;     // Pointer to first array, [n1], to be modified
     int n1;
     int *d2;       // Append this data, [n2], to first array
     int n2;
{
  int i;
  int *t;

  t = (int *)myalloc((n1+n2)*sizeof(int));
  for(i=0;i<n1;i++)
    t[i] = (*pd1)[i];
  for(i=0;i<n2;i++)
    t[n1+i] = d2[i];

  myfree(*pd1);
  *pd1 = t;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 SUM_IARRAY                                */
/*                                                                           */
/*****************************************************************************/
int sum_iarray(data,n,start,duration)
     int *data,n,start,duration;
{
  int i;
  int sum,stop;

  stop = start+duration-1;
  if ((start < 0)||(stop >= n)||(duration < 0))
    exit_error("SUM_IARRAY","Parameter error");
  
  sum = 0;
  for(i=start;i<=stop;i++)
    sum += data[i];

  return sum;
}
/**************************************-**************************************/
/*                                                                           */
/*                              IS_CONST_IARRAY                              */
/*                                                                           */
/*  Return 1 if the iarray contains the same value everywhere.               */
/*                                                                           */
/*****************************************************************************/
int is_const_iarray(data,n)
     int *data,n;
{
  int i;
  int flag;
  
  flag = 1;
  for(i=1;i<n;i++)
    if (data[i] != data[0])
      flag = 0;
  
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                           IS_CONST_COL_2D_IARRAY                          */
/*                                                                           */
/*  Return 1 if data[i][k] contains the save value for every i.              */
/*                                                                           */
/*****************************************************************************/
int is_const_col_2d_iarray(data,n,k)
     int **data,n,k;
{
  int i;
  int flag;
  
  flag = 1;
  for(i=1;i<n;i++)
    if (data[i][k] != data[0][k])
      flag = 0;
  
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                           IS_NONDECREASING_IARRAY                         */
/*                                                                           */
/*  Return 1 if the values are in non-decreasing order.                      */
/*                                                                           */
/*****************************************************************************/
int is_nondecreasing_iarray(data,n)
     int *data,n;
{
  int i;
  int flag;
  
  flag = 1;
  i = 0;
  while((i<(n-1))&&(flag==1))
    if (data[i] > data[i+1])
      flag = 0;
    else
      i += 1;
  
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                            IS_NONNEGATIVE_IARRAY                          */
/*                                                                           */
/*  Return 1 if all values are >= 0.                                         */
/*                                                                           */
/*****************************************************************************/
int is_nonnegative_iarray(data,n)
     int *data,n;
{
  int i;
  
  for(i=0;i<n;i++){
    if (data[i] < 0)
      return 0;
  }
  return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                            CONVERT_2DV_TO_IARRAY                          */
/*                                                                           */
/*  Convert a 2dv iarray to a 1D iarray.                                     */
/*                                                                           */
/*****************************************************************************/
int *convert_2dv_to_iarray(data,n,cnt,rn)
     int **data,n,*cnt;
     int *rn; /* length of returned array */
{
  int i,j,k;
  int *d,nn;

  if (!is_nonnegative_iarray(cnt,n))
    exit_error("CONVERT_2DV_TO_IARRAY","Negative count found");

  nn = sum_iarray(cnt,n,0,n);
  
  d = (int *)myalloc(nn*sizeof(int));
  k = 0;
  for(i=0;i<n;i++)
    for(j=0;j<cnt[i];j++){
      d[k] = data[i][j];
      k += 1;
    }

  *rn = nn;
  return d;
}
/**************************************-**************************************/
/*                                                                           */
/*                               REVERSE_IARRAY                              */
/*                                                                           */
/*  Reverse the order of the elements of the array.                          */
/*                                                                           */
/*****************************************************************************/
void reverse_iarray(data,n)
     int *data,n;
{
  int i;
  int t,mid;

  mid = n/2;
  for(i=0;i<mid;i++){
    t = data[i];
    data[i] = data[n-1-i];
    data[n-1-i] = t;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               ROTATE_IARRAY                               */
/*                                                                           */
/*  Shift all elements in the array to the right (to higher indices) by      */
/*  k indices;                                                               */
/*                                                                           */
/*****************************************************************************/
void rotate_iarray(data,n,k)
     int *data;
     int n,k;
{
  int i,j;
  int *t;

  while(k<0) /*** Change negative rotations into equivalent positive ones. ***/
    k += n;
  while(k>=n)
    k -= n;

  if ((n > 0)&&(k != 0)){
    t = (int *)myalloc(n*sizeof(int));
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
/*                           GET_STRING_FROM_IARRAY                          */
/*                                                                           */
/*  Make a string that contains the elements of the iarray separated by      */
/*  'sepstr'.                                                                */
/*                                                                           */
/*****************************************************************************/
char *get_string_from_iarray(data,n,sepstr)
     int *data,n;
     char sepstr[];
{
  int i,k;
  int slen,nsep;
  char tstr[1024],*t,c0;

  nsep = strlen(sepstr);
  slen = (n-1)*nsep + 1; // Count n-1 separators pluse the \0 at end
  for(i=0;i<n;i++){
    sprintf(tstr,"%d",data[i]);
    slen += strlen(tstr);
  }

  t = (char *)myalloc(slen*sizeof(char));

  k = 0;
  for(i=0;i<(n-1);i++){
    sprintf(tstr,"%d",data[i]);
    sprintf((char *)(t+k),"%s%s",tstr,sepstr);
    k += strlen(tstr) + strlen(sepstr);
  }
  // Nov 2009 - 'cc' on Duke Mac Mini didn't like next line:
  //sprintf((char *)(t+k),"%d\0",data[i]);

  c0 = '\0';
  sprintf((char *)(t+k),"%d%c",data[i],c0);

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                          COUNT_OCCURRENCES_IARRAY                         */
/*                                                                           */
/*  Return the number of times "k" occurs in the array.                      */
/*                                                                           */
/*****************************************************************************/
int count_occurrences_iarray(data,n,k)
     int *data,n,k;
{
  int i;
  int count;

  count = 0;
  for(i=0;i<n;i++)
    if (data[i] == k)
      count += 1;
  
  return count;
}
/**************************************-**************************************/
/*                                                                           */
/*                           IARRAY_PATTERN_SEARCH                           */
/*                                                                           */
/*  Return the indices where the given pattern begins in the iarray.         */
/*                                                                           */
/*****************************************************************************/
void iarray_pattern_search(data,n,pat,np,rlist,rn)
     int *data,n;   // Search in this array
     int *pat,np;   // for this pattern
     int **rlist;   // [*rn] list of indices where pattern begins
     int *rn;       // Number of times pattern was matched
{
  int i,k;
  int ntry,match,cnt,*tlist;

  ntry = n - np + 1;  // Number of positions to check

  tlist = (int *)myalloc(ntry*sizeof(int));

  cnt = 0;
  for(i=0;i<ntry;i++){
    
    k = 0;  // position in 'pat'
    match = 1; // Assume it matches, until we find it doesn't
    while((k < np) && match){  // Does 'pat' match 'stim' at offset 'i'?
      if (data[i+k] != pat[k])
	match = 0; // 'pat' no longer matches 'data' here
      else
	k += 1;    // 'pat' still matches 'data', keep going
    }

    if (match == 1){
      tlist[cnt] = i;
      cnt += 1;
    }
  }

  *rlist = copy_iarray(tlist,cnt);
  *rn = cnt;

  myfree(tlist);
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_INDEX_MIN_SEARCH_IARRAY                       */
/*                                                                           */
/*  Return the index of the first occurrence of the maximum value in iarray  */
/*  which is *less than* or equal to "k".  Return -1 if no such value.       */
/*                                                                           */
/*****************************************************************************/
int get_index_min_search_iarray(data,n,k)
     int *data,n,k;
{
  int i;
  int t,max,tmin;

  i = 0;
  t = tmin = -1;
  max = k+1; /* Added to suppress compiler warning */
  while((i<n)&&(t==-1)){
    if (data[i]==k)
      t = i;
    else{
      if (data[i] < k){
	if (tmin == -1){
	  tmin = i;
	  max = data[i];
	}else if (data[i] > max){
	  tmin = i;
	  max = data[i];
	}
      }
      i+=1;
    }
  }
  
  if (t >= 0)
    return t;
  else
    return tmin;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_INDEX_SEARCH_IARRAY                         */
/*                                                                           */
/*  Return the index of the first occurrence of the element in the iarray.   */
/*  Return -1 if not found.                                                  */
/*                                                                           */
/*****************************************************************************/
int get_index_search_iarray(data,n,k)
     int *data,n,k;
{
  int i;
  int t;

  i = 0;
  t = -1;
  while((i<n)&&(t==-1))
    if (data[i]==k)
      t = i;
    else
      i+=1;
  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_INDEX_SEARCH_NTH_IARRAY                      */
/*                                                                           */
/*  Return the index of the "nn"th occurrence of the element in the iarray.  */
/*  Return -1 if not found.  Note, nn=1 indicates the first occurrence.      */
/*                                                                           */
/*****************************************************************************/
int get_index_search_nth_iarray(data,n,k,nn)
     int *data,n,k,nn;
{
  int i;
  int t,count;

  i = count = 0;
  t = -1;
  while((i<n)&&(t==-1))
    if (data[i]==k){
      count += 1;
      if (count == nn)
	t = i;
      else
	i+=1;
    }else
      i+=1;
  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                        GET_INDEX_SEARCH_RANGE_IARRAY                      */
/*                                                                           */
/*  Return the index of the first occurrence of an element in the iarray     */
/*  that is between min and max, inclusive.  Return -1 if not found.         */
/*                                                                           */
/*****************************************************************************/
int get_index_search_range_iarray(data,n,min,max)
     int *data,n,min,max;
{
  int i;
  int t;

  i = 0;
  t = -1;
  while((i<n)&&(t==-1))
    if ((data[i]>=min)&&(data[i]<=max))
      t = i;
    else
      i+=1;
  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                                DIFF_IARRAYS                               */
/*                                                                           */
/*  Return 1 if iarrays are different in length or content.                  */
/*  Other search words:  EQUAL_IARRAYS  SAME                                 */
/*                                                                           */
/*****************************************************************************/
int diff_iarrays(data1,n1,data2,n2)
     int *data1,n1,*data2,n2;
{
  int i;
  int flag;

  flag = 0;
  if (n1 != n2)
    flag = 1;

  for(i=0;i<n1;i++)
    if (data1[i] != data2[i])
      flag = 1;

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_INDEX_MAX_DIFF_IARRAYS                        */
/*                                                                           */
/*  Return the index at which the iarrays have the maximum difference, or    */
/*  -1 if there is no difference.                                            */
/*                                                                           */
/*****************************************************************************/
int get_index_max_diff_iarrays(data1,data2,n)
     int *data1,*data2,n;
{
  int i,k;
  int max;

  k = -1;
  max = 0;
  for(i=0;i<n;i++)
    if (abs(data1[i]-data2[i]) > max){
      k = i;
      max = abs(data1[i]-data2[i]);
    }
  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_MAX_DIFF_IARRAYS                           */
/*                                                                           */
/*  Return the maximum value ifferent in length or content.                  */
/*                                                                           */
/*****************************************************************************/
int get_max_diff_iarrays(data1,data2,n)
     int *data1,*data2,n;
{
  int k;

  k = get_index_max_diff_iarrays(data1,data2,n);
  if (k>=0)
    return(data1[k] - data2[k]);
  else
    return 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                                ADD_IARRAYS                                */
/*                                                                           */
/*****************************************************************************/
int *add_iarrays(data1,data2,n)
     int *data1,*data2;
     int n;
{
  int i;
  int *data;

  data = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++)
    data[i] = data1[i] + data2[i];
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                              ADD_CONST_IARRAY                             */
/*                                                                           */
/*  Add a constant value to the integer array.                               */
/*                                                                           */
/*****************************************************************************/
void add_const_iarray(data,n,c)
     int *data,n,c;
{
  int i;
  
  for(i=0;i<n;i++)
    data[i] += c;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MULTIPLY_IARRAY                             */
/*                                                                           */
/*  Multiply integer array by a float, round to nearest integer.             */
/*                                                                           */
/*****************************************************************************/
int *multiply_iarray(data,n,f)
     int *data,n;
     float f;
{
  int i;
  int *mdata;
  
  mdata = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++)
    mdata[i] = (int)(f * (float)data[i] + 0.5);
  return mdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                              SUBSAMPLE_IARRAY                             */
/*                                                                           */
/*****************************************************************************/
void subsample_iarray(data,n,rdata,rn,ss,origin)
     int *data,n;     /* original data */
     int **rdata,*rn; /* returned data */
     int ss;          /* original points per new sampling unit */
     int origin;      /* a subsample must be on this point */
{
  int i;
  int pos,start;

  check_int_range(0,n-1,origin,"subsample_iarray: origin");
  check_int_range(0,n-1,ss,"subsample_iarray: ss");

  if (n > 0){
    start = origin - (origin/ss)*ss;
    *rn = origin/ss + 1 + (n-origin-1)/ss;
  }else{
    start = 0;
    *rn = 0;
  }
  *rdata = (int *)myalloc(*rn*sizeof(int));
  pos = 0;
  for (i=start;i<n;i+=ss){
    (*rdata)[pos] = data[i];
    pos += 1;
  }
  if (pos != *rn)
    exit_error("SUBSAMPLE_IARRAY","pos != rn");
}
/**************************************-**************************************/
/*                                                                           */
/*                         COUNT_MIN_RUNS_MIDDLE_IARRAY                      */
/*                                                                           */
/*  Count all runs of "x" of length at least "min_run".  Do not count runs   */
/*  that include the first or last element.                                  */
/*                                                                           */
/*****************************************************************************/
int count_min_runs_middle_iarray(data,n,x,min_run)
     int *data,n,x,min_run;
{
  int k;
  int flag,count,rlen;

  if (min_run < 1)
    exit_error("COUNT_MIN_RUNS_MIDDLE_IARRAY","min_run is less than 1");

  k = 0; /*** Count the number of displacements ***/
  flag = 0;
  while((k<(n-1))&&(!flag)){ /*** Can't start counting until we see a 0 ***/
    if (data[k] != x)
      flag = 1;
    k += 1;
  }
  count = 0; /*** Count every run of length min_run or longer. ***/
  rlen = 0;
  while(k<n){
    if (data[k]!=x){
      if (rlen >= min_run) /* End of a long run. */
	count += 1;
      rlen = 0;
    }else
      rlen += 1;
    k += 1;
  }
  return count;
}
/**************************************-**************************************/
/*                                                                           */
/*                        REPLACE_MIN_RUNS_MIDDLE_IARRAY                     */
/*                                                                           */
/*  Replace all runs of "x" less than "min_run" long with "y".               */
/*  For example:                                                             */
/*    replace_min_runs_iarray(data,n,0,1,3)                                  */
/*  fills in all runs of 2 or less 0's with 1's.  Thus,                      */
/*    01111011110001111001111101111001111110000                              */
/*  becomes                                                                  */
/*    01111111110001111111111111111111111110000                              */
/*                                                                           */
/*  Note:  Runs including the first or last element of the array are not     */
/*         filled, with the idea that we don't know how long this run        */
/*         might be if the data extended beyond the array.  Thus, the        */
/*         word "middle" in the name.                                        */
/*                                                                           */
/*****************************************************************************/
void replace_min_runs_middle_iarray(data,n,x,y,min_run)
     int *data,n,x,y,min_run;
{
  int i,j,k;

  k = 0; /* Not currently at a run.  k stores the current run length. */
  for(i=1;i<n;i++)
    if ((data[i] != x)&&(k > 0)&&(k < min_run)){ /*** Fill in short run. ***/
      for(j=i-k;j<i;j++)
	data[j] = y;
      k = 0;
    }else if ((data[i] == x)&&(data[i-1] != x)){ /*** Start of a run. ***/
      k = 1;
    }else if ((data[i] == x)&&(k >= 1)){ /*** In a run. ***/
      k += 1;
    }
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_COLUMN_2D_IARRAY                           */
/*                                                                           */
/*  Return an iarray which is the "k"th column of the 2d iarray.             */
/*                                                                           */
/*****************************************************************************/
int *get_column_2d_iarray(data,m,n,k)
     int **data,m,n,k;
{
  int i;
  int *col;
  
  col = (int *)myalloc(m*sizeof(int));
  for(i=0;i<m;i++)
    col[i] = data[i][k];
  return col;
}
/**************************************-**************************************/
/*                                                                           */
/*                       GET_DECIMAL_NUMBER_FROM_IARRAY                      */
/*                                                                           */
/*  For example, convert array:  0 0 3 0 1  to the number 301.               */
/*                                                                           */
/*****************************************************************************/
int get_decimal_number_from_iarray(data,n)
     int *data,n;
{
  int i,k;
  int total;

  if (n<=0)
    return 0;

  k = 1;
  total = 0;
  for(i=(n-1);i>=0;i--){
    total += k*data[i];
    k *= 10;
  }
  return total;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_BIT_IARRAY_FROM_INT                         */
/*                                                                           */
/*  Return the array of "n" low order bits for the integer "k".              */
/*                                                                           */
/*****************************************************************************/
void get_bit_iarray_from_int(k,n,rdata)
     int k,n,**rdata;
{
  int i;
  int *data;

  data = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++)
    data[n-1-i] = (k>>i) & 1;

  *rdata = data;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 COUNT_DIGITS                              */
/*                                                                           */
/*  Return the number of digits in the number.  A negative sign counts as    */
/*  one digit.                                                               */
/*                                                                           */
/*****************************************************************************/
int count_digits(k)
     int k;
{
  int count;

  count = 1;
  if (k < 0){
    count += 1;
    k = -k;
  }
  while ((k/10) >= 1){
    k = k/10;
    count += 1;
  }
  return count;
}
/**************************************-**************************************/
/*                                                                           */
/*                           PRINT_IARRAY_80_COLUMN                          */
/*                                                                           */
/*****************************************************************************/
void print_iarray_80_column(fout,data,n)
     FILE *fout;
     int *data,n;
{
  int i;
  int char_count,digits;
  
  if (n > 0){
    char_count = 0;
    for (i=0;i<n;i++){
      digits = count_digits(data[i]);
      if ((char_count + digits + 1) >= 80){
	fprintf(fout,"\n");
	char_count = 0;
      }
      if (char_count == 0){
	fprintf(fout,"%d",data[i]);
	char_count += digits;
      }else{
	fprintf(fout," %d",data[i]);
	char_count += digits + 1;
      }
    }
    fprintf(fout,"\n");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               HEAP_SORT_IARRAY                            */
/*                                                                           */
/*****************************************************************************/
void heap_sort_iarray(data,n)
     int *data,n;
{
  hpsort_int((unsigned long)n,data-1);
}
/**************************************-**************************************/
/*                                                                           */
/*                              BUBBLE_SORT_IARRAY                           */
/*                                                                           */
/*****************************************************************************/
void bubble_sort_iarray(data,n)
     int *data,n;
{
  int i;
  int done,bottom,t;
  
  bottom = n-1;
  done = 0;
  while (!done){
    done = 1;
    for (i=0;i<bottom;i++)
      if (data[i] > data[i+1]){
	done = 0;
	t = data[i];
	data[i] = data[i+1];
	data[i+1] = t;
      }
    bottom -= 1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                        INDEX_BUBBLE_SORT_IARRAY_DOUBLE                    */
/*                                                                           */
/*  Sort based on the most signif 'dhi' and least signif 'dlo'.              */
/*                                                                           */
/*****************************************************************************/
int *index_bubble_sort_iarray_double(dhi,dlo,n)
     int *dhi,*dlo,n;
{
  int i;
  int *index,done,bottom,t;

  index = (int *)myalloc(n*sizeof(int));
  for (i=0;i<n;i++)
    index[i] = i;

  bottom = n-1;
  done = 0;
  while (!done){
    done = 1;
    for (i=0;i<bottom;i++)
      if ((dhi[index[i]] > dhi[index[i+1]]) ||
	  ((dhi[index[i]] == dhi[index[i+1]]) &&
	   (dlo[index[i]]  > dlo[index[i+1]]))){
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
/*                           INDEX_BUBBLE_SORT_IARRAY                        */
/*                                                                           */
/*****************************************************************************/
int *index_bubble_sort_iarray(data,n)
     int *data,n;
{
  int i;
  int *index,done,bottom,t;
  
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
/*                           GET_PERMUTATION_IARRAY                          */
/*                                                                           */
/*  Return the array A with A[i] = i, for i=0..n-1.                          */
/*                                                                           */
/*****************************************************************************/
int *get_permutation_iarray(n)
     int n;
{
  int i;
  int *p;

  p = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++)
    p[i] = i;

  return p;
}
/**************************************-**************************************/
/*                                                                           */
/*                              PERMUTE_IARRAY                               */
/*                                                                           */
/*  Rearrange the elements of "data" according to the permutation            */
/*  array "perm".   "perm[i]" contains the index of the element in           */
/*  "data" that should be put at the ith position.                           */
/*                                                                           */
/*  "perm" must contain the numbers from 0..n-1 with no repeats.             */
/*                                                                           */
/*****************************************************************************/
void permute_iarray(data,n,perm)
     int *data,n,*perm;
{
  int i;
  int *t;

  t = copy_iarray(data,n);
  for(i=0;i<n;i++)
    data[i] = t[perm[i]];

  myfree(t);
}
/**************************************-**************************************/
/*                                                                           */
/*                        GET_NEXT_PERMUTATION_PATTERN                       */
/*                                                                           */
/*  Given:    1 1 0 0 0    1 0 0 0 1                                         */
/*  Returen:  1 0 1 0 0    0 1 1 0 0                                         */
/*                                                                           */
/*****************************************************************************/
int get_next_permutation_pattern(data,n)
     int *data,n;
{
  int i,k;
  int done,p;

  k = n-2;
  p = -1;
  while((p < 0)&&(k >= 0)){
    if ((data[k] == 1)&&(data[k+1] == 0))
      p = k;
    k -= 1;
  }
  if (p >= 0){
    data[p+1] = 1; data[p] = 0; /* Move this "1" to the right. */
    /* Slide all 1's right of this back to the left. */
    done = 0;
    while(!done){
      done = 1;
      for(i=p+3;i<n;i++)
	if ((data[i]==1)&&(data[i-1]==0)){
	  data[i] = 0;
	  data[i-1] = 1;
	  done = 0;
	}
    }
  }

  if (p == -1)
    return 0;
  else
    return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                       INCLUDE_COMPONENT_PERMUTATION                       */
/*                                                                           */
/*  For all -1 in "data", replace with "k" if "perm" at that location is 1.  */
/*                                                                           */
/*****************************************************************************/
void include_component_permutation(data,n,perm,np,k)
     int *data,n,*perm,np,k;
{
  int i,j;

  j = 0; /* Position in perm. */
  for(i=0;i<n;i++){
    if (data[i] == -1){
      if (perm[j] == 1)
	data[i] = k;
      j += 1;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_ALL_PERMUTATIONS_IARRAY                       */
/*                                                                           */
/*  Return an array that contains every possible permutation of the number   */
/*  of each integer 0,1,2,... specified.                                     */
/*                                                                           */
/*  "k" - The number of different integers involved.                         */
/*  "cnt" - the multiplicity of each of the "k" integers.                    */
/*                                                                           */
/*****************************************************************************/
void get_all_permutations_iarray(k,cnt,rperm,rn)
     int k,*cnt,***rperm,*rn;
{
  int i,j,j0,j1,l,m,m1;
  int n,**perm,nperm,*t,flag,***pc,*nc;
  int nr0,nr1;

  n = sum_iarray(cnt,k,0,k); /* Length of any permutation. */

  /*** Compute the total number of permutations, and the number of each ***/
  /*** of the "k" component permutations. ***/
  nc = (int *)myalloc(k*sizeof(int));
  j = n;
  nperm = my_rint(n_choose_m(j,cnt[0]));
  nc[0] = nperm;
  for(i=1;i<k;i++){
    j -= cnt[i-1];
    nc[i] = my_rint(n_choose_m(j,cnt[i]));
    nperm *= nc[i];
  }

  /*** Get each of the component permutations in 0,1 form ***/
  l = n;
  pc = (int ***)myalloc((k-1)*sizeof(int **));
  for(i=0;i<(k-1);i++){
    pc[i] = (int **)myalloc(nc[i]*sizeof(int *));
    t = get_zero_iarray(l);
    for(j=0;j<cnt[i];j++)
      t[j] = 1;

    pc[i][0] = copy_iarray(t,l);
    m = 1;
    flag =  get_next_permutation_pattern(t,l);
    while(flag == 1){
      pc[i][m] = copy_iarray(t,l);
      flag =  get_next_permutation_pattern(t,l);
      m += 1;
    }
    myfree(t);
    l -= cnt[i];
  }

  /*** Create and initialize the permutation array. ***/
  perm = get_2d_iarray(nperm,n);
  for(i=0;i<nperm;i++)
    for(j=0;j<n;j++)
      perm[i][j] = -1;

  /*** Combine the components into the final permutation array. ***/
  l = n;
  nr1 = nperm;
  for(i=0;i<(k-1);i++){
    m1 = 0;
    nr0 = nperm/nr1;
    nr1 = nperm/(nc[i]*nr0);
    /*** Note, nr0 * nc[i] * nr1 == nperm. ***/
    for(j0=0;j0<nr0;j0++) /* Outer "prev" repeats. */
      for(j=0;j<nc[i];j++)
	for(j1=0;j1<nr1;j1++){ /* Inner "next" repeats. */
	  include_component_permutation(perm[m1],n,pc[i][j],l,i);
	  m1 += 1;
	}
    l -= cnt[i];
  }

  for(i=0;i<nperm;i++)
    for(j=0;j<n;j++)
      if (perm[i][j] == -1)
	perm[i][j] = k-1;

  for(i=0;i<(k-1);i++)
    free_2d_iarray(pc[i],nc[i]);
  myfree(nc);
  
  *rperm = perm; *rn = nperm;
}
/**************************************-**************************************/
/*                                                                           */
/*                         CONVERT_MULTI_INDEX_TO_1D                         */
/*                                                                           */
/*  Given an set of indices 'ndx' specifying a point in an n-dimensional     */
/*  space, return an index into a 1D array for which the first dimension     */
/*  is taken as the "most significant bit".                                  */
/*                                                                           */
/*****************************************************************************/
int convert_multi_index_to_1d(n,cnt,ndx)
     int n;      // The number of dimensions
     int *cnt;   // [n] the length of each dimension
     int *ndx;   // [n] the index along each dimension
{
  int i,k;
  int nn;

  nn = 1;
  k = 0;
  for(i=(n-1);i>=0;i--){   // For each dimension
    k += ndx[i] * nn;
    nn *= cnt[i];
  }

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                       GET_ND_CROSS_PRODUCT_LIST_IARRAY                    */
/*                                                                           */
/*  Get the n-dimensional cross-product list for the given numbers of        */
/*  items.                                                                   */
/*                                                                           */
/*  n - the number of component sets.                                        */
/*  cnt - the number of items in each component set [n].                     */
/*  list - the array of all cross-product sets [ntot][n].                    */
/*                                                                           */
/*****************************************************************************/
void get_nd_cross_product_list_iarray(n,cnt,rlist,rntot)
     int n,*cnt,***rlist,*rntot;
{
  int i,j,k;
  int done,**list,ntot,*counter;

  /*
    printf("GET_ND_CROSS_PRODUCT_LIST_IARRAY\n");
    printf("    n= %d\n",n);
    printf("    Counts: ");
    for(i=0;i<n;i++)
    printf(" %d",cnt[i]);
    printf("\n");*/

  ntot = 1; // Compute total number of cross-product sets.
  for(i=0;i<n;i++)
    ntot *= cnt[i];

  list = get_2d_iarray(ntot,n);
  counter = get_zero_iarray(n);

  k = 0;
  for(i=0;i<ntot;i++){
    for(j=0;j<n;j++){
      list[i][j] = counter[j];
    }
    counter[k] += 1;
    done = 0; // Handle overflow
    while (!done){
      if (counter[k] == cnt[k]){ // Overflow
	counter[k] = 0;
	k += 1;
	if (k < n)
	  counter[k] += 1;
	else
	  done = 1;
      }else
	done = 1;
    }
    k = 0; // Reset counter to increment low order index
  }
  myfree(counter);

  *rlist = list; *rntot = ntot;
}
/**************************************-**************************************/
/*                                                                           */
/*                     GET_ND_CROSS_PRODUCT_LIST_IARRAY_REV                  */
/*                                                                           */
/*  Make the last in list change most frequently, rather than first, as      */
/*  done in the non-rev version.                                             */
/*                                                                           */
/*  n - the number of component sets.                                        */
/*  cnt - the number of items in each component set [n].                     */
/*  list - the array of all cross-product sets [ntot][n].                    */
/*                                                                           */
/*  *** WYETH - Nov 8, 2019                                                  */
/*  *** This is now the routine called to vary parameters in 'wm'            */
/*                                                                           */
/*****************************************************************************/
void get_nd_cross_product_list_iarray_rev(n,cnt,rlist,rntot)
     int n,*cnt,***rlist,*rntot;
{
  int i,j,k;
  int done,**list,ntot,*counter;

  ntot = 1; // Compute total number of cross-product sets.
  for(i=0;i<n;i++)
    ntot *= cnt[i];

  list = get_2d_iarray(ntot,n);
  counter = get_zero_iarray(n);

  //k = 0;
  k = n-1;
  for(i=0;i<ntot;i++){
    for(j=0;j<n;j++){
      list[i][j] = counter[j];
    }
    counter[k] += 1;
    done = 0; // Handle overflow
    while (!done){
      if (counter[k] == cnt[k]){ // Overflow
	counter[k] = 0;
	//k += 1;
	//if (k < n)
	k -= 1;
	if (k >= 0)
	  counter[k] += 1;
	else
	  done = 1;
      }else
	done = 1;
    }
    //k = 0; // Reset counter to increment low order index
    k = n-1; // Reset counter to increment high order index
  }
  myfree(counter);

  *rlist = list; *rntot = ntot;
}
/**************************************-**************************************/
/*                                                                           */
/*                              GET_UNIQUE_IARRAY                            */
/*                                                                           */
/*  *** WARNING: n^2 algorithm, could be done faster with sorting.           */
/*                                                                           */
/*  Return an iarray with only one occurrence of each entry in the input     */
/*  iarray.                                                                  */
/*                                                                           */
/*****************************************************************************/
int *get_unique_iarray(data,n,rn)
     int *data;
     int n;
     int *rn;
{
  int i,j,k;
  int *udata,*rdata,temp,done;
  
  udata = (int *)myalloc(n*sizeof(int));
  udata[0] = data[0]; // first one is unique
  k = 1;
  for(i=1;i<n;i++){
    temp = data[i];
    j = 0;
    done = 0;
    while(!done){
      if (j==k)
	done = 1;
      else if (temp == udata[j])
	done = 1;
      else
	j += 1;
    }
    if (j==k){ // Found a new value
      udata[k] = temp;
      k += 1;
    }
  }
  rdata = copy_iarray(udata,k);
  myfree(udata);

  *rn = k;
  return rdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_UNIQUE_SORT_IARRAY                         */
/*                                                                           */
/*  Return a sorted iarray with only one occurrence of each entry in the     */
/*  input array.                                                             */
/*                                                                           */
/*****************************************************************************/
int *get_unique_sort_iarray(data,n,rn)
     int *data,n,*rn;
{
  int i,k;
  int *sdata,*udata;

  sdata = copy_iarray(data,n);
  heap_sort_iarray(sdata,n);

  k = 1; /* First one is unique */
  for(i=1;i<n;i++)
    if (sdata[i] != sdata[i-1])
      k += 1;

  /*printf("GET_UNIQUE_SORT_IARRAY:  %d unique values\n",k);*/
  
  udata = (int *)myalloc(k*sizeof(int));
  udata[0] = sdata[0];
  k = 1; /* First one is unique */
  for(i=1;i<n;i++)
    if (sdata[i] != sdata[i-1]){
      udata[k] = sdata[i];
      k += 1;
    }
  
  myfree(sdata);

  *rn = k;
  return udata;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_UNIQUE_COUNT_IARRAY                        */
/*                                                                           */
/*  Count the number of times each unique value occurs in the iarray.        */
/*  Returns:                                                                 */
/*    runiq[rn] - unique values, sorted                                      */
/*    rcnt[rn]  - number of times each value occurred.                       */
/*    rn        - number of unique values.                                   */
/*                                                                           */
/*****************************************************************************/
void get_unique_count_iarray(data,n,runiq,rcnt,rn)
     int *data,n;
     int **runiq,**rcnt,*rn;
{
  int i,j,k;
  int *uniq,*cnt,nu;

  uniq = get_unique_sort_iarray(data,n,&nu);
  cnt = get_zero_iarray(nu);
  for(i=0;i<nu;i++){
    k = uniq[i];
    for(j=0;j<n;j++){
      if (data[j] == k)
	cnt[i] += 1;
    }
  }

  *runiq = uniq; *rcnt = cnt; *rn = nu;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MEAN_SDEV_FREQ_HIST                           */
/*                                                                           */
/*  Compute the mean and standard deviation of the frequency histogram.      */
/*                                                                           */
/*****************************************************************************/
void mean_sdev_freq_hist(data,n,mean,sdev)
     int *data,n;
     float *mean,*sdev;
{
  int i;
  float m,s,weight,totw;
  
  m = 0.0;
  totw = 0.0;
  for(i=0;i<n;i++){
    weight = (float)data[i];
    m += weight * (float)i;
    totw += weight;
  }
  if (totw>0) m /= totw;
  
  s = 0.0;
  for(i=0;i<n;i++){
    weight = (float)data[i];
    s += weight*(((float)i-m)*((float)i-m));
  }
  if (totw>1)
    s = sqrt(s/(float)(totw-1.0));
  else
    s = -1.0;
  
  *mean = m;
  *sdev = s;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MEAN_SDEV_IARRAY                             */
/*                                                                           */
/*  Take an array of ints and compute the mean and standard                  */
/*  deviation.                                                               */
/*                                                                           */
/*****************************************************************************/
void mean_sdev_iarray(data,n,mean,sdev)
     int *data,n;
     float *mean,*sdev;
{
  int i;
  float m,s;
  
  m = 0.0;
  for(i=0;i<n;i++)
    m += (float)data[i];
  if (n>0) m /= (float)n;
  
  s = 0.0;
  for(i=0;i<n;i++)
    s += (((float)data[i]-m)*((float)data[i]-m));
  if (n>1) s = sqrt(s/(float)(n-1));
  
  *mean = m;
  *sdev = s;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MAKE_MEAN_ZERO_IARRAY                          */
/*                                                                           */
/*  Subtract the mean from each element of the iarray.                       */
/*                                                                           */
/*****************************************************************************/
void make_mean_zero_iarray(data,n,oldmean)
     int *data,n;
     float *oldmean;
{
  int i;
  float mean,sdev,fdata;

  mean_sdev_iarray(data,n,&mean,&sdev);
  for(i=0;i<n;i++){
    fdata = (float)data[i];
    data[i] = (int)(fdata - mean + 0.5);
  }
  *oldmean = mean;
}
/**************************************-**************************************/
/*                                                                           */
/*                                  F2IARRAY                                 */
/*                                                                           */
/*  Return an iarray containing the float data.                              */
/*                                                                           */
/*****************************************************************************/
int *f2iarray(data,n)
     float *data;
     int n;
{
  int i;
  int *idata;

  idata = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++)
    idata[i] = (int)data[i];
  return idata;
}
/**************************************-**************************************/
/*                                                                           */
/*                               F2IARRAY_ROUND                              */
/*                                                                           */
/*  Round the float array to integer data.                                   */
/*                                                                           */
/*****************************************************************************/
int *f2iarray_round(data,n)
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
/*                           COUNT_NON_ZERO_IARRAY                           */
/*                                                                           */
/*  Return the number of array elements that are non-zero.  Used for         */
/*  counting flagged elements.                                               */
/*                                                                           */
/*****************************************************************************/
int count_non_zero_iarray(data,n)
     int *data,n;
{
  int i;
  int count;

  count = 0;
  for(i=0;i<n;i++)
    if (data[i] != 0)
      count += 1;
  return count;
}
/**************************************-**************************************/
/*                                                                           */
/*                            SUM_COLUMN_2D_IARRAY                           */
/*                                                                           */
/*  Return the sum of the "k"th column of the 2d iarray.                     */
/*                                                                           */
/*****************************************************************************/
int sum_column_2d_iarray(data,m,n,k)
     int **data,m,n,k;
{
  int sum,*col;
  
  col = get_column_2d_iarray(data,m,n,k);
  sum = sum_iarray(col,m,0,m);
  myfree(col);

  return sum;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MIN_OF_IARRAY                               */
/*                                                                           */
/*****************************************************************************/
int min_of_iarray(data,n)
     int *data,n;
{
  int i,min;

  if (n > 0){
    min = data[0];
    for(i=0;i<n;i++)
      if (data[i] < min)
	min = data[i];
  }else
    min = 0;
  return min;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MAX_OF_IARRAY                               */
/*                                                                           */
/*****************************************************************************/
int max_of_iarray(data,n)
     int *data,n;
{
  int i,max;

  if (n > 0){
    max = data[0];
    for(i=0;i<n;i++)
      if (data[i] > max)
	max = data[i];
  }else
    max = 0;
  return max;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MAX_COORD_IARRAY                             */
/*                                                                           */
/*****************************************************************************/
int max_coord_iarray(data,n)
     int *data;
     int n;
{
  int i,k,max;
  
  if (n < 0){
    printf("(max_coord_iarray):  Array length < 0.\n");
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
/*                              MIN_COORD_IARRAY                             */
/*                                                                           */
/*****************************************************************************/
int min_coord_iarray(data,n)
     int *data;
     int n;
{
  int i,k,min;
  
  if (n < 0) {
    printf("(min_coord_iarray):  Array length < 0.\n");
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
/*                             MAX_OF_2D_IARRAY                              */
/*                                                                           */
/*****************************************************************************/
int max_of_2d_iarray(data,n,m)
     int **data,n,m;
{
  int i,j;
  int max;

  if ((n <= 0) || (m <= 0))
    exit_error("MAX_OF_2D_IARRAY","Empty array");
  
  max = data[0][0];
  for(i=0;i<n;i++){
    for(j=0;j<m;j++){
      if (data[i][j] > max)
	max = data[i][j];
    }
  }
  
  return max;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MIN_MAX_2D_IARRAY                             */
/*                                                                           */
/*****************************************************************************/
void min_max_2d_iarray(data,xn,yn,rmin,rmax)
     int **data;
     int xn,yn;
     int *rmin,*rmax;
{
  int i,j;
  int min,max;

  min = max = data[0][0];
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      if (data[i][j] > max)
	max = data[i][j];
      else if (data[i][j] < min)
	min = data[i][j];
    }
  }

  *rmin = min;
  *rmax = max;
}
/**************************************-**************************************/
/*                                                                           */
/*                           PRINT_MIN_MAX_2D_IARRAY                         */
/*                                                                           */
/*****************************************************************************/
void print_min_max_2d_iarray(data,xn,yn,name)
     int **data;
     int xn,yn;
     char *name;
{
  int i,j;
  int n_nan;
  int min,max;

  printf("===> ARRAY:  %s\n",name);

  //
  //  Check for NaN values
  //
  n_nan = 0;
  min = max = data[0][0];
  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++){
      if (isnan((float)(data[i][j]))){
	printf("  Element %d %d  is NaN\n",i,j);
	n_nan += 1;
	if (n_nan > 10){
	  printf(" ...  exiting\n");
	  exit(0);
	}
      }else if (data[i][j] > max){
	max = data[i][j];
      }else if (data[i][j] < min){
	min = data[i][j];
      }
    }

  printf("  [%d,%d]  min,max  %d %d\n",xn,yn,min,max);
}
/**************************************-**************************************/
/*                                                                           */
/*                              PRINT_2D_IARRAY                              */
/*                                                                           */
/*  Print the contents of the iarray in matrix format.                       */
/*                                                                           */
/*****************************************************************************/
void print_2d_iarray(data,m,n)
     int **data;
     int m,n;
{
  int i,j;

  for(i=0;i<m;i++){
    for(j=0;j<n;j++)
      printf("%d ",data[i][j]);
    printf("\n");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                        READ_IARRAY_PTR_UNTIL_STRING                       */
/*                                                                           */
/*****************************************************************************/
void read_iarray_ptr_until_string(fin,estr,nmax,rdata,rn)
     FILE *fin;       // Read from this file pointer
     char *estr;      // Until this string is found; NULL for EOF
     int nmax;        // Maximum number of data values
     int **rdata;     // [*rn] Returned data
     int *rn;         // Number of values returned
{
  int n,done;
  int *tdata;
  char tstr[SLEN];

  tdata = (int *)myalloc(nmax*sizeof(int));

  n = 0;

  if (estr == NULL){
    //
    //  Read until EOF
    //
    while (fscanf(fin,"%s",tstr) != EOF){
      tdata[n] = atoi(tstr);
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
	exit_error("READ_IARRAY_PTR_UNTIL_STRING","End string not found");
      }else{
	if (strcmp(tstr,estr)==0){
	  done = 1;
	}else{
	  tdata[n] = atoi(tstr);
	  n += 1;
	}
      }
    }
  }

  if (n == 0){
    *rdata = NULL;
  }else{
    *rdata = copy_iarray(tdata,n);
  }

  *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                                READ_IARRAY                                */
/*                                                                           */
/*****************************************************************************/
void read_iarray(infile,data,n)
     char infile[];
     int **data,*n;
{
  FILE *fopen(),*fin;
  int i;
  int temp,ns;

  if ((fin = fopen(infile,"r"))==NULL){
    printf("  *** File %s\n",infile);
    exit_error("READ_IARRAY","Can't open file");
  }
  *n = 0;
  while (fscanf(fin,"%d",&temp)!= EOF)
    *n += 1;
  fclose(fin);

  *data = (int *)myalloc(*n*sizeof(int));

  fin = fopen(infile,"r");
  for(i=0;i<*n;i++)
    ns = fscanf(fin,"%d",&(*data)[i]);
  fclose(fin);
}
/**************************************-**************************************/
/*                                                                           */
/*                                WRITE_IARRAY                               */
/*                                                                           */
/*****************************************************************************/
void write_iarray(outfile,data,n)
     char outfile[];
     int *data,n;
{
  FILE *fopen(),*fout;
  int i;

  if ((fout = fopen(outfile,"w"))==NULL){
    printf("  *** File %s\n",outfile);
    exit_error("WRITE_IARRAY","Can't open file");
  }
  for (i=0;i<n;i++)
    fprintf(fout,"%d\n",data[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                             WRITE_IARRAY_PLOT                             */
/*                                                                           */
/*****************************************************************************/
void write_iarray_plot(outfile,data,n)
     char outfile[];
     int *data,n;
{
  FILE *fopen(),*fout;
  int i;

  fout = fopen(outfile,"w");
  for (i=0;i<n;i++)
    fprintf(fout,"%d %d\n",i,data[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                            WRITE_IARRAY_XY_PLOT                           */
/*                                                                           */
/*****************************************************************************/
void write_iarray_xy_plot(outfile,xdata,ydata,n)
     char outfile[];
     int *xdata,*ydata;
     int n;
{
  FILE *fopen(),*fout;
  int i;

  fout = fopen(outfile,"w");
  for (i=0;i<n;i++)
    fprintf(fout,"%d %d\n",xdata[i],ydata[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                           APPEND_IARRAY_XY_PLOT                           */
/*                                                                           */
/*****************************************************************************/
void append_iarray_xy_plot(outfile,name,xdata,ydata,n)
     char outfile[],name[];
     int *xdata,*ydata;
     int n;
{
  FILE *fopen(),*fout;
  int i;

  fout = fopen(outfile,"a");
  fprintf(fout,"/newplot\n");
  fprintf(fout,"/plotname %s\n",name);
  for (i=0;i<n;i++)
    fprintf(fout,"%d %d\n",xdata[i],ydata[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                            APPEND_IARRAY_PLOT                             */
/*                                                                           */
/*****************************************************************************/
void append_iarray_plot(outfile,name,data,n,xscale)
     char outfile[],name[];
     int *data,n,xscale;
{
  FILE *fopen(),*fout;
  int i;

  fout = fopen(outfile,"a");
  fprintf(fout,"/newplot\n");
  fprintf(fout,"/plotname %s\n",name);
  for (i=0;i<n;i++)
    fprintf(fout,"%d %d\n",i*xscale,data[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                         APPEND_IARRAY_PLOT_OFFSET                         */
/*                                                                           */
/*****************************************************************************/
void append_iarray_plot_offset(outfile,name,data,n,xscale,xoffset)
     char outfile[],name[];
     int *data,n,xscale,xoffset;
{
  FILE *fopen(),*fout;
  int i;

  fout = fopen(outfile,"a");
  fprintf(fout,"/newplot\n");
  fprintf(fout,"/plotname %s\n",name);
  for (i=0;i<n;i++)
    fprintf(fout,"%d %d\n",xoffset + i*xscale,data[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                             LIN_INTERP_IARRAY                             */
/*                                                                           */
/*  Use linear interpolation to return a value for the specified             */
/*  coordinate.  Return the closest value if "inx" is outside of             */
/*  domain of array.                                                         */
/*                                                                           */
/*****************************************************************************/
float lin_interp_iarray(x,y,n,inx)
     int *x,*y;
     int n;
     float inx;
{
  int i;

  printf("*** NOT TESTED*** Copied from similar routine in farray_util.c\n");

  if (inx < (float)x[0])
    return (float)y[0];
  else if (inx > (float)x[n-1])
    return (float)y[n-1];
  else{
    i=0;
    while(i<(n-1)){
  printf("%d %d %d %f\n",i,x[i],y[i],inx);
      if ((float)x[i] == inx)
	return (float)y[i];
      else if (((float)x[i] < inx)&&((float)x[i+1] > inx)){
	return (float)y[i] + (inx-(float)x[i])/(float)(x[i+1]-x[i]) *
	  (float)(y[i+1]-y[i]);
      }else
	i += 1;
    }
  }
  exit_error("LIN_INTERP_IARRAY","Error");
  return 0.0; /*** This makes lint happy. ***/
}
/**************************************-**************************************/
/*                                                                           */
/*                              PASCAL_TRIANGLE                              */
/*                                                                           */
/*  Bernoulli triangle.  Binomial coefficents.                               */
/*                                                                           */
/*  n=1             1                                                        */
/*  n=2            1 1                                                       */
/*  n=3           1 2 1                                                      */
/*  n=4          1 3 3 1                                                     */
/*  n=5         1 4 6 4 1                                                    */
/*  ...            ...                                                       */
/*                                                                           */
/*****************************************************************************/
int *pascal_triangle(n)
     int n;
{
  int i,j;
  int *p1,*p2,*t;

  if (n < 1)
    return NULL;

  p1 = get_zero_iarray(n);
  p2 = get_zero_iarray(n);

  p1[0] = 1;
  for(i=1;i<n;i++){
    p2[0] = 1;
    for(j=1;j<=i;j++){
      p2[j] = p1[j-1]+p1[j];
    }
    t = p1;
    p1 = p2;
    p2 = t;
  }

  myfree(p2);

  return p1;
}
