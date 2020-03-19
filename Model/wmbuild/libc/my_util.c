/*****************************************************************************/
/*                                                                           */
/*  my_util.c                                                                */
/*  wyeth bair                                                               */
/*  caltech                                                                  */
/*  02/23/93                                                                 */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>  // for 'scandir'

#include "my_util.h"

/**************************************-**************************************/
/*                                                                           */
/*                                   MYMINI                                  */
/*                                                                           */
/*****************************************************************************/
int mymini(a,b)
     int a,b;
{
  if (a <= b)
    return a;
  return b;
}
/**************************************-**************************************/
/*                                                                           */
/*                                   MYMAXI                                  */
/*                                                                           */
/*****************************************************************************/
int mymaxi(a,b)
     int a,b;
{
  if (a >= b)
    return a;
  return b;
}
/**************************************-**************************************/
/*                                                                           */
/*                                   MYMINF                                  */
/*                                                                           */
/*****************************************************************************/
float myminf(a,b)
     float a,b;
{
  if (a <= b)
    return a;
  return b;
}
/**************************************-**************************************/
/*                                                                           */
/*                                   MYMAXF                                  */
/*                                                                           */
/*****************************************************************************/
float mymaxf(a,b)
     float a,b;
{
  if (a >= b)
    return a;
  return b;
}
/**************************************-**************************************/
/*                                                                           */
/*                                MIN_OF_THRE                                */
/*                                                                           */
/*****************************************************************************/
float min_of_three(a,b,c)
     float a,b,c;
{
  float min;

  if (a <= b)
    min = a;
  else
    min = b;

  if (c < min)
    min = c;

  return min;
}
/**************************************-**************************************/
/*                                                                           */
/*                                MAX_OF_THRE                                */
/*                                                                           */
/*****************************************************************************/
float max_of_three(a,b,c)
     float a,b,c;
{
  float max;

  if (a >= b)
    max = a;
  else
    max = b;

  if (c > max)
    max = c;

  return max;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 EXIT_ERROR                                */
/*                                                                           */
/*****************************************************************************/
void exit_error(proc_name,statement)
     char proc_name[],statement[];
{
  printf("  *** %s:  %s.  Exiting.\n",proc_name,statement);
  exit(1);
}
/**************************************-**************************************/
/*                                                                           */
/*                                     EE                                    */
/*                                                                           */
/*****************************************************************************/
void ee(proc_name,statement)
     char proc_name[],statement[];
{
  printf("  *** %s:  %s.  Exiting.\n",proc_name,statement);
  exit(1);
}
/**************************************-**************************************/
/*                                                                           */
/*                                  MYALLOC                                  */
/*                                                                           */
/*  On 4/23/95, this routine returns a NULL pointer if n <= 0.  This         */
/*  change was made because on the DEC Alpha, malloc returns NULL for        */
/*  n=0, but apparently the SUNs did not return a NULL value (thus,          */
/*  things worked on the SUNs, but stopped on the Alpha.)                    */
/*                                                                           */
/*****************************************************************************/
void *myalloc(n)
     int n;
{
  void *temp;

  if (n > 0){
    temp = (void *)malloc((unsigned)n);

    if (temp == NULL)
      exit_error("MYALLOC","Null pointer from malloc");
  }else
    temp = NULL;
    
  return temp;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 GET_FARRAY                                */
/*                                                                           */
/*****************************************************************************/
float *get_farray(n)
     int n;
{
  return (float *)myalloc(n*sizeof(float));
}
/**************************************-**************************************/
/*                                                                           */
/*                                 GET_DARRAY                                */
/*                                                                           */
/*****************************************************************************/
double *get_darray(n)
     int n;
{
  return (double *)myalloc(n*sizeof(double));
}
/**************************************-**************************************/
/*                                                                           */
/*                                  MY_FOPEN                                 */
/*                                                                           */
/*****************************************************************************/
FILE *my_fopen(fname,mode,caller)
     char *fname;     // File to open
     char *mode;      // "a", "r", ...
     char *caller;    // Name of calling routine
{
  FILE *fp,*fopen();
  int flag;

  if ((fp = fopen(fname,mode)) == NULL){
    printf("*** (my_util.c, my_fopen) caller:  %s\n",caller);
    printf("*** Cannot open file:  %s\n",fname);
    exit(0);
  }

  return fp;
}
/**************************************-**************************************/
/*                                                                           */
/*                                   MYLOG                                   */
/*                                                                           */
/*  Append a string to a log file, or print the string if no log file name   */
/*  is given.                                                                */
/*                                                                           */
/*  Return 0 if the specified log file name could not be opened.             */
/*                                                                           */
/*****************************************************************************/
int mylog(fname,s)
     char fname[];    // Output logfile to append
     char s[];        // String to log
{
  FILE *fout,*fopen();
  int flag;

  flag = -1;
  if (fname == NULL){
    printf("%s",s);
    flag = 1;
  }else if (fname[0] == '\0'){
    printf("%s",s);
    flag = 1;
  }else{
    if ((fout = fopen(fname,"a")) != NULL){
      fprintf(fout,"%s",s);
      fclose(fout);
      flag = 1;
    }else
      flag = 0;
  }
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 MYLOG_EXIT                                */
/*                                                                           */
/*  Append a string to a log file, or print the string if no log file name   */
/*  is given, and exit.                                                      */
/*                                                                           */
/*****************************************************************************/
void mylog_exit(fname,s)
     char fname[];    // Output logfile to append
     char s[];        // String to log
{
  (void)mylog(fname,s);
  (void)mylog(fname,"\n");  // Add a carriage return
  exit(0);
}
/**************************************-**************************************/
/*                                                                           */
/*                                   MYLOGX                                  */
/*                                                                           */
/*  Like MYLOG_EXIT, but with multiple parameters.                           */
/*                                                                           */
/*****************************************************************************/
void mylogx(fname,procname,str)
     char fname[];      // Output logfile to append
     char procname[];   // Procedure that has exited
     char str[];        // String to log
{
  char tstr[LONG_SLEN];

  sprintf(tstr,"*** %s  %s  Exiting.\n",procname,str);

  (void)mylog(fname,tstr);
  exit(0);
}
/**************************************-**************************************/
/*                                                                           */
/*                                 MYLOG_XCI                                 */
/*                                                                           */
/*  Exit conditional on the value of an integer.                             */
/*                                                                           */
/*****************************************************************************/
void mylog_xci(fname,procname,parname,a,lo,hi)
     char *fname;       // Output logfile to append
     char *procname;    // Procedure that has exited
     char *parname;     // Procedure that has exited
     int  a;            // Value to check
     int  lo,hi;        // Range, inclusive
{
  char tstr[LONG_SLEN];

  if ((a < lo) || (a > hi)){
    
    sprintf(tstr,"*** The value %d of %s is not in the range %d, %d\n",
	    a,parname,lo,hi);
    (void)mylog(fname,tstr);

    sprintf(tstr,"*** %s  Exiting.\n",parname);
    (void)mylog(fname,tstr);
    exit(0);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                 MYLOG_IVAL                                */
/*                                                                           */
/*****************************************************************************/
void mylog_ival(fname,msg,ival)
     char fname[];      // Output logfile to append
     char *msg;         // Message to print before the value
     int ival;          // value to print
{
  char tstr[LONG_SLEN];

  sprintf(tstr,"%s %d\n",msg,ival);

  (void)mylog(fname,tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                                 MYLOG_FVAL                                */
/*                                                                           */
/*****************************************************************************/
void mylog_fval(fname,msg,fval)
     char fname[];      // Output logfile to append
     char *msg;         // Message to print before the value
     float fval;        // value to print
{
  char tstr[LONG_SLEN];

  sprintf(tstr,"%s %f\n",msg,fval);

  (void)mylog(fname,tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                                 MYLOG_CVAL                                */
/*                                                                           */
/*****************************************************************************/
void mylog_cval(fname,msg,cval)
     char fname[];      // Output logfile to append
     char *msg;         // Message to print before the value
     char *cval;        // value to print
{
  char tstr[LONG_SLEN];

  sprintf(tstr,"%s %s\n",msg,cval);

  (void)mylog(fname,tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                                  MY_RINT                                  */
/*                                                                           */
/*  Round to nearest integer.                                                */
/*                                                                           */
/*****************************************************************************/
int my_rint(x)
     float x;
{
  int n;

  if (x >= 0.0)
    n = (int)(0.5 + x);
  else
    n = (int)(-0.5 + x);

  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MY_RINT_DOUBLE                              */
/*                                                                           */
/*  Round to nearest integer.                                                */
/*                                                                           */
/*****************************************************************************/
int my_rint_double(x)
     double x;
{
  int n;

  if (x >= 0.0)
    n = (int)(0.5 + x);
  else
    n = (int)(-0.5 + x);

  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                                  MY_LOG2                                  */
/*                                                                           */
/*****************************************************************************/
double my_log2(x)
     double x;
{
  return log(x)/log(2.0);
}
/**************************************-**************************************/
/*                                                                           */
/*                                  MY_SIGN                                  */
/*                                                                           */
/*****************************************************************************/
int my_sign(x)
     double x;
{
  if (x >= 0.0)
    return 1;
  else
    return -1;
}
/**************************************-**************************************/
/*                                                                           */
/*                                POWER_OF_TWO                               */
/*                                                                           */
/*   Returns the log base 2 of an integer if it is an integral power of two, */
/*   -1 otherwise.                                                           */
/*                                                                           */
/*****************************************************************************/
int power_of_two(x)
     int x;
{
  int r,p,total;
  
  r = -1;
  if (x > 0){
    p = 0;
    total = 1;
    while (total < x){
      total = total * 2;
      p += 1;
    }
    if (total==x)
      r = p;
  }
  return r;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 N_CHOOSE_M                                */
/*                                                                           */
/*  *** NO GOOD FOR LARGE "n" ***                                            */
/*  Return the number of ways to choose "m" of "n" items.                    */
/*           n!                                                              */
/*       ---------                                                           */
/*       m! (n-m)!                                                           */
/*                                                                           */
/*****************************************************************************/
float n_choose_m(n,m)
     int m,n;
{
  float prob;
  int i,nn;

  if ((m == n) || (m == 0))
    prob = 1.0;
  else if (m > n){
    prob = 0.0;
    exit_error("N_CHOOSE_M","m > n");
  }else{
    if ((n-m) < m)
      nn = n-m;
    else
      nn = m;
    prob = 1;
    for (i=0;i<nn;i++)
      prob *= (float)(n-i)/(float)(nn-i);
  }
  return prob;
}
/**************************************-**************************************/
/*                                                                           */
/*                            COUNT_ONES_BINARY_INT                          */
/*                                                                           */
/*  Return the number of ones in a binary integer.                           */
/*                                                                           */
/*****************************************************************************/
int count_ones_binary_int(x)
     int x;
{
  int i;
  int n,nx,t;

  if (x < 0)
    exit_error("COUNT_ONES_BINARY_INT","VALUE IS NEGATIVE");

  nx = 8*sizeof(int);

  n = 0;
  t = x;
  for(i=0;i<nx;i++){
    t = x >> i;
    if (t & 1)
      n += 1;
  }

  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_BINARY_ENTROPY                            */
/*                                                                           */
/*  Return the entropy in bits, given the probability 'p'.                   */
/*                                                                           */
/*****************************************************************************/
float get_binary_entropy(p)
     float p;
{
  float e;

  if ((p < 0.0)||(p > 1.0))
    exit_error("GET_BINARY_ENTROPY","p out of range");

  e = 0.0;
  if ((p > 0.0) && (p < 1.0)){
    e +=      p  * (float)my_log2(1.0/(double)p);
    e += (1.0-p) * (float)my_log2(1.0/(double)(1.0-p));
  }
  return e;
}
/**************************************-**************************************/
/*                                                                           */
/*                               REVERSE_BYTES                               */
/*                                                                           */
/*  For changing endian.  Rerverse the ordering of the "n" bytes beginning   */
/*  at address "p".                                                          */
/*                                                                           */
/*****************************************************************************/
void reverse_bytes(p,n)
     unsigned char *p;
     int n;
{
  int i;
  unsigned char t;

  for(i=0;i<n/2;i++){
    t = *(p+i);
    *(p+i) = *(p+n-1-i);
    *(p+n-1-i) = t;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                 REV_FREAD                                 */
/*                                                                           */
/*  This version of "fread" changes endian.  "fread" returns 0 on end of     */
/*  file.                                                                    */
/*                                                                           */
/*****************************************************************************/
int rev_fread(ptr,size,nitems,stream)
     char *ptr;
     int size,nitems;
     FILE *stream;
{
  int i;
  int nread;

  nread = fread(ptr,size,nitems,stream);
  for(i=0;i<nitems;i++)
    reverse_bytes((unsigned char *)(ptr + i*size),size);

  return nread;
}
/**************************************-**************************************/
/*                                                                           */
/*                               REVFLAG_FREAD                               */
/*                                                                           */
/*  This version of "fread" changes endian depeding on 'revflag'.            */
/*                                                                           */
/*****************************************************************************/
int revflag_fread(ptr,size,nitems,stream,revflag)
     char *ptr;
     int size,nitems;
     FILE *stream;
     int revflag;
{
  int i;
  int nread;

  nread = fread(ptr,size,nitems,stream);
  if (revflag == 1)
    for(i=0;i<nitems;i++)
      reverse_bytes((unsigned char *)(ptr + i*size),size);
  return nread;
}
/**************************************-**************************************/
/*                                                                           */
/*                               REVFLAG_FWRITE                              */
/*                                                                           */
/*  Reverse the byte order if 'revflag' is 1.                                */
/*  Return the number of ITEMS (not bytes) written.                          */
/*                                                                           */
/*****************************************************************************/
int revflag_fwrite(ptr,size,nitems,stream,revflag)
     char *ptr;
     int size,nitems;
     FILE *stream;
     int revflag;
{
  int i;
  int nwrite;

  if (revflag == 0){
    nwrite = fwrite((char *)ptr,size,nitems,stream);
  }else{
    nwrite = 0;
    for(i=0;i<nitems;i++){
      reverse_bytes((unsigned char *)(ptr + i*size),size);
      nwrite += fwrite((char *)(ptr + i*size),size,1,stream);
      reverse_bytes((unsigned char *)(ptr + i*size),size);
    }
  }
  return nwrite;
}
/**************************************-**************************************/
/*                                                                           */
/*                                GET_2D_CARRAY                              */
/*                                                                           */
/*****************************************************************************/
char **get_2d_carray(m,n)
     int m,n;
{
  int i;
  char **c;

  c = (char **)myalloc(m*sizeof(char *));
  for(i=0;i<m;i++)
    c[i] = (char *)myalloc(n*sizeof(char));

  return c;
}
/**************************************-**************************************/
/*                                                                           */
/*                               FREE_2D_CARRAY                              */
/*                                                                           */
/*****************************************************************************/
void free_2d_carray(c,n)
     char **c;
     int n;
{
  int i;

  if (c != NULL){
    for(i=0;i<n;i++){
      if (c[i] != NULL)
	myfree(c[i]);
    }
    myfree(c);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              GET_ZERO_FARRAY                              */
/*                                                                           */
/*****************************************************************************/
float *get_zero_farray(n)
     int n;
{
  int i;
  float *f;

  f = (float *)myalloc(n*sizeof(float));
  for (i=0;i<n;i++)
    f[i] = 0.0;
  return f;
}
/**************************************-**************************************/
/*                                                                           */
/*                              GET_CONST_FARRAY                             */
/*                                                                           */
/*****************************************************************************/
float *get_const_farray(n,c)
     int n;
     float c;
{
  int i;
  float *f;

  f = (float *)myalloc(n*sizeof(float));
  for (i=0;i<n;i++)
    f[i] = c;
  return f;
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_ZERO_IARRAY                             */
/*                                                                           */
/*****************************************************************************/
int *get_zero_iarray(n)
     int n;
{
  int i,*ia;

  ia = (int *)myalloc(n*sizeof(int));
  for (i=0;i<n;i++)
    ia[i] = 0;
  return ia;
}
/**************************************-**************************************/
/*                                                                           */
/*                              GET_CONST_IARRAY                             */
/*                                                                           */
/*****************************************************************************/
int *get_const_iarray(n,c)
     int n,c;
{
  int i,*ia;

  ia = (int *)myalloc(n*sizeof(int));
  for (i=0;i<n;i++)
    ia[i] = c;
  return ia;
}
/**************************************-**************************************/
/*                                                                           */
/*                               CHECK_INT_RANGE                             */
/*                                                                           */
/*****************************************************************************/
void check_int_range(a,b,x,name)
     int a,b,x;
     char *name;
{
  if ((x<a)||(x>b)){
    printf("  *** CHECK_INT_RANGE:  (%s) %d not in [%d,%d]\n",name,x,a,b);
    exit_error("CHECK_INT_RANGE","Range error");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                REMOVE_FILE                                */
/*                                                                           */
/*****************************************************************************/
void remove_file(filename)
     char filename[];
{
  int rval;
  char command[SLEN];

  sprintf(command,"rm -f %s",filename);
  rval = system(command);
  if (rval == -1)
    exit_error("REMOVE_FILE","System command returned -1");
}
/**************************************-**************************************/
/*                                                                           */
/*                                   IS_DIR                                  */
/*                                                                           */
/*  Return 1 if it is a directory.                                           */
/*                                                                           */
/*****************************************************************************/
int is_dir(fname)
     char *fname;
{
  struct stat stat_p;

  // WYETH THIS DOESN"T SEEM TO WORK????
  // WYETH THIS DOESN"T SEEM TO WORK????   Nov 2009
  // WYETH THIS DOESN"T SEEM TO WORK????

  stat(fname,&stat_p);
  /*
    printf("fname %s  st_mode  %o\n",fname,(int)(stat_p.st_mode));
    printf("ISDIR:  %d\n",S_ISDIR(stat_p.st_mode));
    printf("ISREG:  %d\n",S_ISREG(stat_p.st_mode));*/
  
  return S_ISDIR(stat_p.st_mode);
}
/**************************************-**************************************/
/*                                                                           */
/*                                  IS_FILE                                  */
/*                                                                           */
/*  Return 1 if it is a regular file.                                        */
/*                                                                           */
/*****************************************************************************/
int is_file(fname)
     char *fname;
{
  struct stat stat_p;

  stat(fname,&stat_p);
  
  return S_ISREG(stat_p.st_mode);
}
/**************************************-**************************************/
/*                                                                           */
/*                                  FILELIST                                 */
/*                                                                           */
/*****************************************************************************/
char **filelist(dir,rn)
     char  *dir;
     int   *rn;
{
  int i;
  int n;
  char **flist;
  struct dirent **namelist;
  
  n = scandir(dir,&namelist,0,alphasort);  // See man scandir
  if (n < 0){
    n = -1;
    flist = NULL;
  }else if (n == 0){
    flist = NULL;
  }else{
    flist = (char **)myalloc(n*sizeof(char *));
    for(i=0;i<n;i++){
      /*printf("%s\n", namelist[i]->d_name);*/
      flist[i] = strdup(namelist[i]->d_name);
      free(namelist[i]);
    }
    free(namelist);
  }

  *rn = n;
  return flist;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MAKE_FILELIST                               */
/*                                                                           */
/*****************************************************************************/
char **make_filelist(dir,rn)
     char  *dir;
     int   *rn;
{
  int i;
  int rval;
  FILE *fopen(),*fin;
  char command[SLEN],**flist,name[SLEN];

  sprintf(command,"ls -1 %s > .zz.wytemp.filelist",dir);
  rval = system(command);
  if (rval == -1)
    exit_error("MAKE_FILELIST","System command (ls) returned -1");

  // count the files in the data directory
  fin = fopen(".zz.wytemp.filelist","r");
  *rn = 0;
  while (fscanf(fin,"%*s") != EOF)
    *rn += 1;
  fclose(fin);

  // store all filenames
  flist = (char **)myalloc(*rn*sizeof(char *));
  fin = fopen(".zz.wytemp.filelist","r");
  for (i=0;i<*rn;i++) {
    rval = fscanf(fin,"%s",name);
    flist[i] = strdup(name);
  }
  fclose(fin);

  rval = system("rm .zz.wytemp.filelist");
  if (rval == -1)
    exit_error("MAKE_FILELIST","System command (rm) returned -1");
  
  return flist;
}
/**************************************-**************************************/
/*                                                                           */
/*                                FILE_EXISTS                                */
/*                                                                           */
/*  This may not be technically correct, but it should suffice for           */
/*  determining the existence of a file.                                     */
/*                                                                           */
/*  Return 1 if the file exists (i.e., can be opened to read), 0 otherwise.  */
/*                                                                           */
/*****************************************************************************/
int file_exists(fname)
     char *fname;
{
  FILE *fopen(),*fin;

  if ((fin = fopen(fname,"r"))==NULL){
    return 0;
  }else{
    fclose(fin);
    return 1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               FILE_NOT_EMPTY                              */
/*                                                                           */
/*  Return 1 if the file exists and is not empty.                            */
/*                                                                           */
/*****************************************************************************/
int file_not_empty(fname)
     char *fname;
{
  FILE *fopen(),*fin;
  char c;

  if ((fin = fopen(fname,"r"))==NULL)
    return 0;
  else{
    if (fread(&c,sizeof(char),1,fin) == 1){
      fclose(fin);
      return 1;
    }else{
      fclose(fin);
      return 0;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                CREATE_FILE                                */
/*                                                                           */
/*****************************************************************************/
void create_file(fname)
     char *fname;
{
  FILE *fopen(),*fout;

  if ((fout = fopen(fname,"w"))==NULL)
    exit_error("CREATE_FILE","Could not open file");
  else
    fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                               PROCESS_EXISTS                              */
/*                                                                           */
/*  Return 1 if the process exists.                                          */
/*                                                                           */
/*****************************************************************************/
int process_exists(pid)
     int pid;
{
  int  priority;
  extern int errno;

  errno = 0;
  priority = getpriority(PRIO_PROCESS,pid);
  if (errno > 0){
    return 0;
  }else
    return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                                MY_TRUE_TEST                               */
/*                                                                           */
/*****************************************************************************/
int my_true_test()
{
  return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                               SMART_DIV_SIZE                              */
/*                                                                           */
/*****************************************************************************/
double smart_div_size(ymax)
     double ymax;
{
  int i;
  double ty,logyi,y1,ydiv,pf;
  int smart_divisors[6];
  
  ty = ymax;
  if (ty < 0.0)
    ty = -ty;
  
  logyi = floor(log10(ty));
  
  pf = pow(10.0,logyi);
  
  // 'y1' should now be in the range from 1 to 10 ??
  y1  = ty / pf;
  
  //System.out.println("logyi = " + logyi + "  ty = " + ty);
  //System.out.println("y1 = " + y1);
  
  for(i=0;i<=5;i++)
    smart_divisors[i] = 0;
  
  if (y1 > 8.0){
    ydiv = 8.0;
    smart_divisors[4] = 1;
    smart_divisors[2] = 1;
  }else if (y1 > 6.0){
    ydiv = 6.0;
    smart_divisors[3] = 1;
    smart_divisors[2] = 1;
  }else if (y1 > 5.0){
    ydiv = 5.0;
    smart_divisors[5] = 1;
  }else if (y1 > 4.0){
    ydiv = 4.0;
    smart_divisors[4] = 1;
    smart_divisors[2] = 1;
  }else if (y1 > 3.0){
    ydiv = 3.0;
    smart_divisors[3] = 1;
  }else if (y1 > 2.0){
    ydiv = 2.0;
    smart_divisors[2] = 1;
  }else{
    ydiv = 1.0;
  }
  
  ydiv *= pf;
  
  return ydiv;
}
