/*****************************************************************************/
/*                                                                           */
/*  carray_util.c                                                            */
/*  wyeth bair                                                               */
/*  caltech                                                                  */
/*  03/19/93                                                                 */
/*                                                                           */
/*  "2d_carray" refers to char **a.                                          */
/*  "3d_carray" refers to char ***a.                                         */
/*  also, 'slist' can refer to a list of strings (char **)                   */
/*  "string_list_carray" refers to an array of variable-length arrays of     */
/*                      characters, thus ***a, where a[n][cnt[n]].           */
/*                                                                           */
/*  All names should have "carray" or "string" in them.                      */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "my_util.h"
#include "carray_util.h" // Defines 'dll_slist' and 'strel'

/**************************************-**************************************/
/*                                                                           */
/*                             DLL_STRING_GET_INIT                           */
/*                                                                           */
/*****************************************************************************/
struct dll_slist *dll_string_get_init()
{
  struct dll_slist *d;

  d = (struct dll_slist *)myalloc(sizeof(struct dll_slist));

  d->n = 0;
  d->first = d->last = NULL;

  return d;
}
/**************************************-**************************************/
/*                                                                           */
/*                               DLL_STRING_FREE                             */
/*                                                                           */
/*****************************************************************************/
void dll_string_free(d)
     struct dll_slist *d;
{
  struct strel *t,*t1;

  t = d->first;
  while(t!=NULL){
    t1 = t;
    t = t->next;
    myfree(t1->s);
    myfree(t1);
  }
  myfree(d);
}
/**************************************-**************************************/
/*                                                                           */
/*                               DLL_STRING_PRINT                            */
/*                                                                           */
/*****************************************************************************/
void dll_string_print(d)
     struct dll_slist *d;
{
  struct strel *t;

  t = d->first;
  while(t!=NULL){
    printf("%s\n",t->s);
    t = t->next;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             DLL_STRING_WRITE_FILE                         */
/*                                                                           */
/*****************************************************************************/
void dll_string_write_file(d,outfile)
     struct dll_slist *d;
     char outfile[];
{
  struct strel *t;
  FILE *fopen(),*fout;

  if ((fout = fopen(outfile,"w"))==NULL)
    exit_error("DLL_STRING_WRITE_FILE","Cannot open file");

  t = d->first;
  while(t!=NULL){
    fprintf(fout,"%s\n",t->s);
    t = t->next;
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                             DLL_STRING_INSERT                             */
/*                                                                           */
/*  Insert 'str' into the DLL 'd' if it is not already in the list.          */
/*                                                                           */
/*****************************************************************************/
void dll_string_insert(d,str)
     struct dll_slist *d;
     char *str;
{
  struct strel *t,*tt;

  if (d->n == 0){
    d->first = (struct strel *)myalloc(sizeof(struct strel));
    d->first->s = strdup(str);
    d->first->prev = d->first->next = NULL;
    d->n += 1;
    d->last = d->first;
  }else{
    t = d->first;
    while((strcmp(t->s,str)<0)&&(t->next != NULL))
      t = t->next;
    if (strcmp(t->s,str)==0){ /*** Found it, done ***/
      ;
    }else{ /*** Add new string to list. ***/
      tt = (struct strel *)myalloc(sizeof(struct strel));
      tt->s = strdup(str);
      d->n += 1;
      if (strcmp(t->s,str)<0){ /*** Append. ***/
	tt->prev = t;
	tt->next = NULL;
	t->next = tt;
	d->last = tt;
      }else{
	tt->prev = t->prev;
	tt->next = t;
	t->prev = tt;
	if (tt->prev != NULL){ /*** Insert. ***/
	  tt->prev->next = tt;
	}else{
	  d->first = tt; /*** Prepend. ***/
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              DLL_STRING_SEARCH                            */
/*                                                                           */
/*  Return a POINTER to the string element, or NULL if not found.            */
/*                                                                           */
/*****************************************************************************/
struct strel *dll_string_search(d,str)
     struct dll_slist *d;
     char *str;
{
  struct strel *t,*tt;

  tt = NULL;
  if (d->n > 0){
    t = d->first;
    while((strcmp(t->s,str)!=0)&&(t->next != NULL))
      t = t->next;
    if (strcmp(t->s,str)==0){ // Found it, done
      tt = t;
    }
  }
  return tt;
}
/**************************************-**************************************/
/*                                                                           */
/*                             DLL_STRING_APPEND                             */
/*                                                                           */
/*  Append 'str' to the DLL.                                                 */
/*                                                                           */
/*****************************************************************************/
void dll_string_append(d,str,repflag)
     struct dll_slist *d;
     char *str;
     int repflag;  // 1-allow repetitions
{
  int fflag;
  struct strel *t,*tt;

  if (d->n == 0){
    d->first = (struct strel *)myalloc(sizeof(struct strel));
    d->first->s = strdup(str);
    d->first->prev = d->first->next = NULL;
    d->n += 1;
    d->last = d->first;
  }else{
    if (dll_string_search(d,str) == NULL)
      fflag = 0;
    else
      fflag = 1;

    if ((fflag == 0) || ((fflag == 1) && repflag)){
      tt = (struct strel *)myalloc(sizeof(struct strel));
      tt->s = strdup(str);
      tt->prev = d->last;
      tt->next = NULL;
      d->last->next = tt;
      d->last = tt;
      d->n += 1;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            DLL_STRING_GET_SLIST                           */
/*                                                                           */
/*****************************************************************************/
char **dll_string_get_slist(d)
     struct dll_slist *d;
{
  int i;
  int n;
  char **slist;
  struct strel *t;

  n = d->n;
  slist = (char **)myalloc(n*sizeof(char *));

  i = 0;
  t = d->first;
  while(t!=NULL){
    slist[i] = strdup(t->s);
    t = t->next;
    i += 1;
  }

  return slist;
}
/**************************************-**************************************/
/*                                                                           */
/*                        DLL_STRING_GET_NTH_STRING_PTR                      */
/*                                                                           */
/*  Return a pointer to the string for the nth element, or NULL if there is  */
/*  no nth element.                                                          */
/*                                                                           */
/*****************************************************************************/
char *dll_string_get_nth_string_ptr(d,n)
     struct dll_slist *d;
     int n;
{
  int i;
  char *s;
  struct strel *t;

  if (d == NULL)
    exit_error("DLL_STRING_GET_NTH_STRING_PTR","dll is null");

  if ((d->first == NULL) && (d->n > 0))
    exit_error("DLL_STRING_GET_NTH_STRING_PTR","First is null");


  s = NULL;

  if (n < d->n){
    t = d->first;
    for(i=0;i<n;i++){
      t = t->next;
      if (t == NULL)
	exit_error("DLL_STRING_GET_NTH_STRING_PTR","Unexpected null");
    }
    s = t->s;
  }

  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                              STRING_COUNT_CHAR                            */
/*                                                                           */
/*****************************************************************************/
int string_count_char(s,c)
     char *s;   // In this string
     char c;    // Count the number of times this char occurs
{
  int i;
  int n,cnt;

  cnt = 0;
  n = strlen(s);
  for(i=0;i<n;i++)
    if (s[i] == c)
      cnt += 1;

  return cnt;
}
/**************************************-**************************************/
/*                                                                           */
/*                            STRING_MAKE_LOWER_CASE                         */
/*                                                                           */
/*****************************************************************************/
void string_make_lower_case(s)
     char *s;
{
  int i,k;
  int n;

  /* A-Z = 65 - 90 */
  /* a-z = 97 - 122 */

  n = strlen(s);
  for(i=0;i<n;i++){
    k = (int)s[i];
    if ((k >= 65) && (k <= 90))
      s[i] = (char)(k+32);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                IS_DIGIT_CHAR                              */
/*                                                                           */
/*  Return 1 if the character is in 0..9                                     */
/*                                                                           */
/*****************************************************************************/
int is_digit_char(c)
     char c;
{
  if (((int)c >= (int)'0') && ((int)c <= (int)'9'))
    return 1;
  else
    return 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                               IS_NUMERIC_CHAR                             */
/*                                                                           */
/*  Return 1 if the character is the first character of a number.            */
/*                                                                           */
/*****************************************************************************/
int is_numeric_char(c)
     char c;
{
  int i;
  int n;
  char numchar[13];

  strcpy(numchar,"0123456789.-");
  n = strlen(numchar);

  for(i=0;i<n;i++)
    if (c==numchar[i])
      return 1;

  return 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 IS_NAME_CHAR                              */
/*                                                                           */
/*  Return 1 if the character is alpha, digit or '_'.                        */
/*  Used in names of C routines.                                             */
/*                                                                           */
/*****************************************************************************/
int is_name_char(c)
     char c;
{
  int i;
  int n;
  char namechar[SLEN];

  strcpy(namechar,
	 "1234567890_ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
  n = strlen(namechar);
  
  for(i=0;i<n;i++)
    if (c == namechar[i])
      return 1;

  return 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                               IS_NUMBER_STRING                            */
/*                                                                           */
/*  Return 1 if the string is a number.                                      */
/*                                                                           */
/*****************************************************************************/
int is_number_string(s)
     char *s;
{
  int i;
  int n,pflag,flag,okend;
  char c,numchar[13];

  flag = 1; // Assume it is a number
  n = strlen(s);
  if (n < 1)
    flag = 0;

  /* pflag  0  -  First char, must find '-', '.' or digit
            1  -  After '-', must find '.' or digit
            2  -  Past first digit before any decimal
            3  -  After '.', must find digit
            4  -  Past first digit after decimal
            5  -  After 'e'
            6  -  after '+' or '-', following 'e', must find digit
            7  -  after first number in exponent
  */

  pflag = 0;
  okend = 0; // Is it OK if it ends here
  i = 0;
  while ((i < n) && (flag == 1)){
    c = s[i];
    if (pflag == 0){ // Check first char
      if (is_digit_char(c)){
	pflag = 2;
	okend = 1;
      }else if (c == '.')
	pflag = 3;
      else if (c == '-')
	pflag = 1;
      else
	flag = 0;
      /*printf("pf 0  flag=%d\n",flag);*/
    }else if (pflag == 1){
      okend = 0;
      if (is_digit_char(c)){
	pflag = 2;
	okend = 1;
      }else if (c == '.')
	pflag = 3;
      else
	flag = 0;
      /*printf("pf 1  flag=%d\n",flag);*/
    }else if (pflag == 2){
      if (c == '.'){
	pflag = 3;
	okend = 1;
      }else if (c == 'e'){
	pflag = 5;
	okend = 0;
      }else if (!is_digit_char(c))
	flag = 0;
      /*printf("pf 2  flag=%d\n",flag);*/
    }else if (pflag == 3){
      if (is_digit_char(c)){
	pflag = 4;
	okend = 1;
      }else
	flag = 0;
      /*printf("pf 3  flag=%d\n",flag);*/
    }else if (pflag == 4){
      if (c == 'e'){
	pflag = 5;
	okend = 0;
      }else if (!is_digit_char(c))
	flag = 0;
      /*printf("pf 4  flag=%d\n",flag);*/
    }else if (pflag == 5){
      if ((c == '-')||(c == '+')){
	pflag = 6;
	okend = 0;
      }else
	flag = 0;
      /*printf("pf 5  flag=%d\n",flag);*/
    }else if (pflag == 6){
      if (is_digit_char(c)){
	pflag = 7;
	okend = 1;
      }else
	flag = 0;
      /*printf("pf 6  flag=%d\n",flag);*/
    }else if (pflag == 7){
      if (!is_digit_char(c))
	flag = 0;
      /*printf("pf 7  flag=%d\n",flag);*/
    }
    i += 1;
  }

  if (!okend){ // If number did not end OK, set flag to 0
    flag = 0;
  }

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_2D_POINTER_CARRAY                          */
/*                                                                           */
/*****************************************************************************/
char ***get_2d_pointer_carray(m,n)
     int m,n;
{
  int i;
  char ***c;

  c = (char ***)myalloc(m*sizeof(char **));
  for(i=0;i<m;i++)
    c[i] = (char **)myalloc(n*sizeof(char *));

  return c;
}
/**************************************-**************************************/
/*                                                                           */
/*                           FREE_2D_POINTER_CARRAY                          */
/*                                                                           */
/*  (This is really free_3d_carray.)                                         */
/*                                                                           */
/*****************************************************************************/
void free_2d_pointer_carray(data,m,n)
     char ***data;
     int m,n;
{
  int i,j;

  if (data != NULL){
    for(i=0;i<m;i++){
      for(j=0;j<n;j++)
	myfree(data[i][j]);
      myfree(data[i]);
    }
    myfree(data);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                      FREE_2D_POINTER_VAR_COUNT_CARRAY                     */
/*                                                                           */
/*  WYETH - CHECK THIS !!!   THIS ROUTINE SEEMS TO HAVE DISAPPEARED!!!       */
/*  And then  I rewrote it.                                                  */
/*                                                                           */
/*  SEEMS TO BE THE SAME AS 'free_string_list_carray' below.                 */
/*                                                                           */
/*****************************************************************************/
void free_2d_pointer_var_count_carray(data,n,cnt)
     char ***data;
     int n,*cnt;
{
  int i,j;

  if (data != NULL){
    for(i=0;i<n;i++){
      for(j=0;j<cnt[i];j++)
	myfree(data[i][j]);
      myfree(data[i]);
    }
    myfree(data);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          FREE_STRING_LIST_CARRAY                          */
/*                                                                           */
/*****************************************************************************/
void free_string_list_carray(data,n,cnt)
     char ***data;
     int n,*cnt;
{
  int i,j;

  for(i=0;i<n;i++){
    for(j=0;j<cnt[i];j++)
      myfree(data[i][j]);
    myfree(data[i]);
  }
  myfree(data);
}
/**************************************-**************************************/
/*                                                                           */
/*                            APPEND_NULL_TO_CARRAY                          */
/*                                                                           */
/*  Append the null char, '\0', to the carray, so it is a string.            */
/*                                                                           */
/*****************************************************************************/
void append_null_to_carray(pdata,n)
     char **pdata;
     int n;
{
  int i;
  char *t;

  t = (char *)myalloc((n+1)*sizeof(char));
  for(i=0;i<n;i++)
    t[i] = (*pdata)[i];
  t[n] = '\0';
  myfree(*pdata);
  *pdata = t;
}
/**************************************-**************************************/
/*                                                                           */
/*                               COPY_2D_CARRAY                              */
/*                                                                           */
/*****************************************************************************/
char **copy_2d_carray(c,n)
     char **c;
     int n;
{
  int i;
  char **t;

  t = (char **)myalloc(n*sizeof(char *));
  for(i=0;i<n;i++)
    t[i] = strdup(c[i]);

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                              CONCAT_2D_CARRAYS                            */
/*                                                                           */
/*****************************************************************************/
char **concat_2d_carrays(c1,n1,c2,n2,rn)
     char **c1;
     int n1;
     char **c2;
     int n2,*rn;
{
  int i;
  int n;
  char **t;

  n = n1 + n2;
  t = (char **)myalloc(n*sizeof(char *));
  for(i=0;i<n1;i++)
    t[i] = strdup(c1[i]);
  for(i=0;i<n2;i++)
    t[n1+i] = strdup(c2[i]);

  *rn = n;
  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                            IS_CONST_2D_CARRAY                             */
/*                                                                           */
/*  Return 1 if all strings in the list are the same, 0 otherwise.           */
/*                                                                           */
/*  WYETH - CAN'T handle NULL strings.                                       */
/*                                                                           */
/*****************************************************************************/
int is_const_2d_carray(c,n)
     char **c;
     int n;
{
  int i;
  int done;

  done = 0;
  i = 1;
  while(!done){
    if (i>=n)
      done = 1;
    else{
      if (strcmp(c[i],c[0])==0)
	i += 1;
      else{
	done = 1;
      }
    }
  }
  if (i==n){
    return 1;
  }else{
    return 0;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               IS_SUFFIX_STRING                            */
/*                                                                           */
/*****************************************************************************/
int is_suffix_string(s,suf)
     char s[],suf[];
{
  int i,k;
  int n1,n2,min;

  k = 0;

  n1 = strlen(s);
  n2 = strlen(suf);
  if (n1 >= n2){
    k = 1;
    for(i=0;i<n2;i++)
      if (s[n1-n2 + i] != suf[i])
	k = 0;
  }

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                            COMPARE_PREFIX_STRING                          */
/*                                                                           */
/*****************************************************************************/
int compare_prefix_string(s1,s2)
     char s1[],s2[];
{
  int i,k;
  int n1,n2,min;

  n1 = strlen(s1);
  n2 = strlen(s2);
  if (n1 <= n2)
    min = n1;
  else
    min = n2;

  k = 1;
  for(i=0;i<min;i++)
    if (s1[i] != s2[i])
      k = 0;
    
  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                          COMPARE_PREFIX_STRING_ORDER                      */
/*                                                                           */
/*  Return 1 if 's1' is a prefix of 's2', 0 otherwise.                       */
/*                                                                           */
/*****************************************************************************/
int compare_prefix_string_order(s1,s2)
     char s1[],s2[];
{
  int i,k;
  int n1,n2,min;

  n1 = strlen(s1);
  n2 = strlen(s2);
  if (n1 <= n2)
    min = n1;
  else
    return 0; // A longer string is not a prefix

  k = 1;
  for(i=0;i<min;i++)
    if (s1[i] != s2[i])
      k = 0;
    
  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                             CARRAY_REPLACE_CHAR                           */
/*                                                                           */
/*  Replace all occurrences of character 'co' with 'cn' in string 's'.       */
/*                                                                           */
/*****************************************************************************/
void carray_replace_char(s,co,cn)
     char *s;   // String in which to replace character
     char co;   // Old char (to be replaced)
     char cn;   // New char to write over 'co'
{
  char *t;

  t = strchr(s,co);
  while(t != NULL){
    *t = cn;
    t = strchr(s,co);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              REPLACE_HTML_GTLT                            */
/*                                                                           */
/*  Replace '>' and '<' with '& g t ;' and '& l t ;'                         */
/*  Return the new string, or NULL if no replacements were made.             */
/*                                                                           */
/*****************************************************************************/
char *replace_html_gtlt(s)
     char *s;
{
  int i,k;
  int n,nn,cnt;
  char *newstr;

  //
  //  Count number of replacements to make
  //
  n = strlen(s);
  cnt = 0;
  for(i=0;i<n;i++){
    if ((s[i] == '<') || (s[i] == '>'))
      cnt += 1;
  }

  newstr = NULL;

  if (cnt > 0){
    nn = n + cnt*3 + 1;   // Length of new string needs additional chars
    newstr = (char *)myalloc(nn*sizeof(char));

    k = 0;
    for(i=0;i<n;i++){
      if (s[i] == '<'){
	newstr[k]   = '&';
	newstr[k+1] = 'l';
	newstr[k+2] = 't';
	newstr[k+3] = ';';
	k += 4;
      }else if (s[i] == '>'){
	newstr[k]   = '&';
	newstr[k+1] = 'g';
	newstr[k+2] = 't';
	newstr[k+3] = ';';
	k += 4;
      }else{
	newstr[k] = s[i];
	k += 1;
      }
    }
    newstr[k] = '\0';
    k += 1;
    if (k != nn)
      exit_error("WYETH","This should not happen");

  }
  return newstr;
}
/**************************************-**************************************/
/*                                                                           */
/*                          REPLACE_HTML_GTLT_2D_CARRAY                      */
/*                                                                           */
/*  Replace '>' and '<' with '& g t ;' and '& l t ;'                         */
/*                                                                           */
/*****************************************************************************/
char **replace_html_gtlt_2d_carray(c,n)
     char **c;
     int n;
{
  int i;
  char **cnew;

  cnew = (char **)myalloc(n*sizeof(char *));

  for(i=0;i<n;i++){
    if (c[i] != NULL){
      cnew[i] = replace_html_gtlt(c[i]);
      if (cnew[i] == NULL)  // No replacements were made
	cnew[i] = strdup(c[i]);
    }else
      cnew[i] = NULL;
  }

  return cnew;
}
/**************************************-**************************************/
/*                                                                           */
/*                            REPLACE_STRING_SEGMENT                         */
/*                                                                           */
/*  Replace a specified segment of a string with a new string, of possibly   */
/*  different length.                                                        */
/*                                                                           */
/*****************************************************************************/
char *replace_string_segment(s,newstr,k,n)
     char *s;        // Within this string
     char *newstr;   // Insert this in place of
     int k;          // The string starting here
     int n;          //   and having this length
{
  int i;
  int nn,tn,sn,en;
  char *repstr;

  tn = strlen(s);
  if (tn < k+n)
    exit_error("REPLACE_STRING_SEGMENT","String too short");

  nn = strlen(newstr);

  sn = tn - n + nn + 1;  // +1 for null char
  repstr = (char *)myalloc(sn*sizeof(char));

  en = tn - (k+n);

  //printf("sn= %d en= %d tn= %d nn= %d n= %d  k=%d\n",sn,en,tn,nn,n,k);

  for(i=0;i<k;i++)
    repstr[i] = s[i];  // Copy everything before the segment

  for(i=0;i<nn;i++)
    repstr[k+i] = newstr[i];  // Write the new string

  for(i=0;i<en;i++)
    repstr[k+nn+i] = s[k+n+i];  // Copy everything after the segment

  repstr[sn-1] = '\0'; // Null to terminate string

  return repstr;
}
/**************************************-**************************************/
/*                                                                           */
/*                               REPLACE_STRING                              */
/*                                                                           */
/*  Replace all occurrences of 'oldstr' with 'newstr' in 's'.                */
/*  The original is unchanged.                                               */
/*  The number of replacements is returned in 'rn'.                          */
/*  If there are no replacements, NULL is returned.                          */
/*                                                                           */
/*****************************************************************************/
char *replace_string(s,oldstr,newstr,rn)
     char *s;
     char *oldstr;
     char *newstr;
     int *rn;
{
  int i,k;
  int n,replen,nnew,nold;
  char *tp,*cp,*repstr,*tt;

  if (oldstr == NULL)
    exit_error("REPLACE_STRING","oldstr is NULL");
  if (newstr == NULL)
    exit_error("REPLACE_STRING","newstr is NULL");
  if (s == NULL)
    exit_error("REPLACE_STRING","string is NULL");

  nnew = strlen(newstr);
  nold = strlen(oldstr);

  //printf("     s:  %s\n",s);
  //printf("   old:  %s\n",oldstr);
  //printf("   new:  %s\n",newstr);

  n = 0;  // Count number of replacements
  tp = strstr(s,oldstr);
  while(tp != NULL){
    n += 1;
    tp += nold;
    tt = strstr(tp,oldstr);
    tp = tt;
  }

  if (n == 0){
    repstr = NULL;
  }else{

    // New length is original length + n * diff, plus 1 for NULL
    replen = 1 + strlen(s) + n * (nnew - nold);
    repstr = (char *)myalloc(replen*sizeof(char));

    cp = s;  // current pointer in 's' to copy from
    k = 0;   // current index in 'repstr' to copy to
    tp = strstr(s,oldstr);
    while(tp != NULL){

      while(cp != tp){  // Copy everything up to the occurrence
	repstr[k] = *cp;
	k += 1;
	cp += 1;
      }
      cp += nold;  // Jump to end of occurrence

      for(i=0;i<nnew;i++){  // Put in new string
	repstr[k] = newstr[i];
	k += 1;
      }

      tp += nold;  // Jump past the 'old' string
      tt = strstr(tp,oldstr);  // Find next occurrence
      tp = tt;
    }
    while(*cp != '\0'){
      repstr[k] = *cp;
      k += 1;
      cp += 1;
    }

    repstr[k] = '\0';  // Termination
    k += 1;
  }

  *rn = n;
  return repstr;
}
/**************************************-**************************************/
/*                                                                           */
/*                           REPLACE_STRING_2D_CARRAY                        */
/*                                                                           */
/*  Replace all occurrences of 'oldstr' with 'newstr' in the list of         */
/*  strings 'c'.  New storage is returned, the original 'c' is unchanged.    */
/*                                                                           */
/*****************************************************************************/
char **replace_string_2d_carray(c,n,oldstr,newstr,rn)
     char **c;       // [n] List of strings
     int n;          // Number of strings in list
     char *oldstr;   // String to be replaced
     char *newstr;   // String to insert
     int *rn;        // Number of replacements made
{
  int i;
  int ntot,nrep;
  char **cnew;

  cnew = (char **)myalloc(n*sizeof(char *));

  ntot = 0;
  for(i=0;i<n;i++){
    if (c[i] != NULL){
      cnew[i] = replace_string(c[i],oldstr,newstr,&nrep);
      if (nrep == 0)
	cnew[i] = strdup(c[i]);
    }else
      cnew[i] = NULL;
    ntot += nrep;
  }

  *rn = ntot;
  return cnew;
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_STRING_FROM_INT                         */
/*                                                                           */
/*  Return the string for the integer 'k'.                                   */
/*                                                                           */
/*****************************************************************************/
char *get_string_from_int(k)
     int k;
{
  char tstr[SLEN];

  sprintf(tstr,"%d",k);
  return strdup(tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                      GET_LEADING_ZERO_STRING_FROM_INT                     */
/*                                                                           */
/*  E.g., 7 -->  007                                                         */
/*                                                                           */
/*****************************************************************************/
char *get_leading_zero_string_from_int(k,n)
     int k;  // The number
     int n;  // How many decimal places
{
  int i;
  int klen;
  char tstr[SLEN],*t;

  t = (char *)myalloc((n+1)*sizeof(char));
  for(i=0;i<n;i++)
    t[i] = '0';
  t[n] = '\0';

  sprintf(tstr,"%d",k);
  klen = strlen(tstr);
  if (klen > n)
    exit_error("GET_LEADING_ZERO_STRING_FROM_INT","n is too small");
  for(i=0;i<klen;i++)
    t[n-1-i] = tstr[klen-1-i];

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_BIT_STRING_FROM_INT                         */
/*                                                                           */
/*  Return the string of '0' and '1' chars corresponding to the last n bits  */
/*  of the number k.                                                         */
/*                                                                           */
/*****************************************************************************/
char *get_bit_string_from_int(k,n)
     int k,n;
{
  int i;
  int *data;
  char *t;

  data = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++)
    data[n-1-i] = (k>>i) & 1;

  t = (char *)myalloc((n+1)*sizeof(char));
  for(i=0;i<n;i++)
    if (data[i])
      t[i] = '1';
    else
      t[i] = '0';
  t[n] = '\0';

  myfree(data);
  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                            PRINT_BIT_STRING_FOR_INT                       */
/*                                                                           */
/*  Print the low order "n" bits of the number "k" from highest to lowest    */
/*  left to right.                                                           */
/*                                                                           */
/*****************************************************************************/
void print_bit_string_for_int(k,n)
     int k,n;
{
  char *s;

  s = get_bit_string_from_int(k,n);
  printf("%s",s);
  myfree(s);
}
/**************************************-**************************************/
/*                                                                           */
/*                         CONVERT_01_STRING_TO_IARRAY                       */
/*                                                                           */
/*  For example, change "01101" to the integer array:  0 1 1 0 1.            */
/*                                                                           */
/*****************************************************************************/
void convert_01_string_to_iarray(cstr,invert_flag,rdata,rn)
     char cstr[];
     int invert_flag,**rdata,*rn;
{
  int i;
  int *data,n;

  n = strlen(cstr);
  data = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++){
    data[i] = (cstr[i]=='1');
    if (invert_flag == 1)
      data[i] = 1-data[i];
  }

  *rdata = data; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                        CONVERT_n0p_STRING_TO_IARRAY                       */
/*                                                                           */
/*  For example, change "n0ppn" to the integer array:  -1 0 1 1 -1.          */
/*                                                                           */
/*****************************************************************************/
void convert_n0p_string_to_iarray(cstr,invert_flag,rdata,rn)
     char cstr[];
     int invert_flag,**rdata,*rn;
{
  int i;
  int *data,n;

  n = strlen(cstr);
  data = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++){
    if (cstr[i]=='n')
      data[i] = -1;
    else if (cstr[i]=='0')
      data[i] = 0;
    else if (cstr[i]=='p')
      data[i] = 1;
    else
      exit_error("CONVERT_n0p_STRING_TO_IARRAY","Bad ternary value");

    if (invert_flag == 1)
      data[i] = -data[i];
  }
  
  *rdata = data; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                       CONVERT_DIGIT_STRING_TO_IARRAY                      */
/*                                                                           */
/*  For example, change "012..8910..." to:  0 1 2 ... 8 9 1 0 ...            */
/*                                                                           */
/*****************************************************************************/
void convert_digit_string_to_iarray(cstr,invert_flag,rdata,rn)
     char cstr[];
     int invert_flag,**rdata,*rn;
{
  int i;
  int *data,n;

  n = strlen(cstr);
  data = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++){
    data[i] = (int)(cstr[i]) - (int)'0';
  }
  
  *rdata = data; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                       CONVERT_MULTI_STRING_TO_IARRAY                      */
/*                                                                           */
/*  For example, change "o1p2-b-b-o1p1" to an integer index array.           */
/*                                                                           */
/*****************************************************************************/
void convert_multi_string_to_iarray(cstr,ndim,dn,nbl,blank,rdata,rn)
     char cstr[];
     int ndim;     // Number of dimensions
     int *dn;      // [ndim] Length of each dimension
     int nbl;      // [ndim] Length of each dimension
     char blank;   // Character used to indicate blank
     int **rdata,*rn;
{
  int i,j,k; 
  int *data,n,slen,si,nstim,nstimreg,*cdn;
  char *t,dummy;

  slen = strlen(cstr);
  if (slen <= 0)
    exit_error("CONVERT_MULTI_STRING_TO_IARRAY","Empty string");

  n = 1;
  for(i=0;i<slen;i++)
    if (cstr[i] == '-')
      n += 1;
  /*printf("  string %s has %d words\n",cstr,n);*/

  /*** Compute total number of stimuli ***/
  nstim = 1;
  for(i=0;i<ndim;i++){
    nstim *= dn[i];
    /*printf("dn[%d] = %d\n",i,dn[i]);*/
  }
  nstimreg = nstim;
  if (nbl > 0)
    nstim += 1;  /* All blanks are the same stimulus */

  /*** Cumulative array */
  cdn = get_zero_iarray(ndim);
  for(i=0;i<ndim;i++){
    cdn[i] = 1;
    for(j=i+1;j<ndim;j++)
      cdn[i] *= dn[j];
  }

  /*printf("  %d stimuli in total.\n",nstim);*/

  data = (int *)myalloc(n*sizeof(int));
  t = cstr;
  printf("    %s  =  ",cstr);
  for(i=0;i<n;i++){
    if (t[0] == blank){
      if (nbl == 0)
	exit_error("CONVERT_MULTI_STRING_TO_IARRAY","No blanks");
      data[i] = nstim - 1; /* Last stimulus category */
      t += 2; /* Move past separator */
    }else{
      /*printf("IN HERE t[0]=%c  blk=%c\n",t[0],blank);*/
      si = 0;
      for(j=0;j<ndim;j++){
	sscanf(t,"%c%d",&dummy,&k);
	t += 1; /* Skip over the character */
	while(is_digit_char(*t)){
	  t += 1; /* Skip over all digits */
	}
	/*printf(" dim %d:  %d\n",j,k);*/
	si += (k-1)*cdn[j];
      }
      if (si >= nstimreg){
	printf("    si = %d  k = %d\n",si,k);
	exit_error("CONVERT_MULTI_STRING_TO_IARRAY","Stim index too high");
      }
	
      data[i] = si;
      t += 1; // Move past separator
    }
    printf(" %d",data[i]);
  }
  printf("\n");
  
  *rdata = data; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                            APPEND_STRING_TO_FILE                          */
/*                                                                           */
/*****************************************************************************/
void append_string_to_file(outfile,data)
     char outfile[];
     char *data;
{
  FILE *fopen(),*fout;

  if ((fout = fopen(outfile,"a"))==NULL){
    printf("  FILENAME:  %s\n",outfile);
    exit_error("APPEND_STRING_TO_FILE","Cannot open file");
  }
  fprintf(fout,"%s",data);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                               IS_INT_STRING                               */
/*                                                                           */
/*  Return 1 if the string represents an integer.                            */
/*                                                                           */
/*****************************************************************************/
int is_int_string(sdata)
     char *sdata;
{
  int i;
  int slen;
  char sdigit[11];

  strcpy(sdigit,"0123456789");

  slen = strlen(sdata);
  for(i=0;i<slen;i++)
    if (strchr(sdigit,sdata[i]) == NULL)
      return 0;
  return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                               IS_FLOAT_STRING                             */
/*                                                                           */
/*  Return 1 if the string represents a floating point value.                */
/*                                                                           */
/*  *** WYETH - more checking should be done here e.g.                       */
/*    - 2 decimal points not allowed                                         */
/*                                                                           */
/*****************************************************************************/
int is_float_string(sdata)
     char *sdata;
{
  int i;
  int slen;
  char sdigit[11];

  sdigit[0] = '0';
  sdigit[1] = '1';
  sdigit[2] = '2';
  sdigit[3] = '3';
  sdigit[4] = '4';
  sdigit[5] = '5';
  sdigit[6] = '6';
  sdigit[7] = '7';
  sdigit[8] = '8';
  sdigit[9] = '9';
  sdigit[10] = '.';

  // Wyeth, 'cc' on Duke Mac Mini Nov 2009 didn't like this line:
  //strcpy(sdigit,"0123456789.");

  slen = strlen(sdata);
  for(i=0;i<slen;i++){
    if ((i==0) && (sdata[i] == '-'))  // WYETH - Added Apr 7, 2014
      ; // Skip over a leading '-' sign
    else if (strchr(sdigit,sdata[i]) == NULL)
      return 0;
  }

  return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                               COUNT_DECIMALS                              */
/*                                                                           */
/*  Return the number of places after the LAST decimal:                      */
/*                                                                           */
/*    "23"   -->  0                                                          */
/*    "1.0"  -->  1                                                          */
/*    "9.15" -->  2                                                          */
/*           ...                                                             */
/*                                                                           */
/*****************************************************************************/
int count_decimals(sdata)
     char *sdata;
{
  int n;
  char *tdata;

  tdata = strrchr(sdata,'.');
  if (tdata == NULL){
    n = 0;
  }else{
    n = strlen(sdata) - ((long)tdata - (long)sdata + 1);
  }

  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                            COUNT_ITEMS_IN_STRING                          */
/*                                                                           */
/*****************************************************************************/
int count_items_in_string(sdata)
     char *sdata;
{
  int n;
  char white[4],*tdata;

  white[0] = (char)9;   // Horiz. Tab
  white[1] = (char)32;  // Space
  white[2] = '\n';      // New line
  white[3] = '\0';

  tdata = strdup(sdata);
  if (strtok(tdata,white) != NULL){
    n = 1;
    while(strtok(NULL,white)!=NULL)
      n += 1;
  }else
    n = 0;
  myfree(tdata);

  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_ITEMS_FROM_STRING                          */
/*                                                                           */
/*****************************************************************************/
void get_items_from_string(sdata,rlist,rn)
     char *sdata;
     char ***rlist;
     int *rn;
{
  int i;
  int n;
  char **list,white[5],*tdata;

  white[0] = (char)9;   // Horiz. Tab
  white[1] = (char)32;  // Space
  white[2] = '\n';      // New line
  white[3] = '\0';
  white[4] = (char)13;  // Carriage return

  tdata = strdup(sdata);
  if (strtok(tdata,white) != NULL){
    n = 1;
    while(strtok(NULL,white)!=NULL)
      n += 1;
  }else
    n = 0;
  myfree(tdata);

  /*
  printf("GET_ITEMS_FROM_STRING:  %s\n",sdata);
  printf("    n= %d\n",n);
  for(i=0;i<strlen(sdata);i++){
    printf("    s[%d] = %d\n",i,(int)sdata[i]);
  }
  */

  tdata = strdup(sdata);
  if (n > 0){
    list = (char **)myalloc(n*sizeof(char *));
    list[0] = strdup(strtok(tdata,white));
    for(i=1;i<n;i++){
      list[i] = strdup(strtok(NULL,white));
    }
  }else
    list = NULL;
  myfree(tdata);

  *rlist = list; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_FIRST_ITEM_FROM_STRING                       */
/*                                                                           */
/*  Return NULL if no items.                                                 */
/*                                                                           */
/*  *** WYETH - this can be very inefficient.                                */
/*                                                                           */
/*****************************************************************************/
char *get_first_item_from_string(s)
     char *s;
{
  int n;
  char **slist,*t;

  get_items_from_string(s,&slist,&n);

  if (n < 1)
    t = (char *)NULL;
  else
    t = strdup(slist[0]);

  free_2d_carray(slist,n);

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_IARRAY_FROM_STRING                         */
/*                                                                           */
/*****************************************************************************/
void get_iarray_from_string(sdata,rdata,rn)
     char *sdata;
     int **rdata,*rn;
{
  int i;
  int n,*data;
  char **list;

  get_items_from_string(sdata,&list,&n);
  data = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++)
    data[i] = atoi(list[i]);

  *rdata = data; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_FARRAY_FROM_STRING                         */
/*                                                                           */
/*****************************************************************************/
void get_farray_from_string(sdata,rdata,rn)
     char *sdata;
     float **rdata;
     int *rn;
{
  int i;
  int n;
  float *data;
  char **list;
  
  get_items_from_string(sdata,&list,&n);
  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = atof(list[i]);
  
  *rdata = data; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_FARRAY_FROM_CARRAY                         */
/*                                                                           */
/*****************************************************************************/
float *get_farray_from_carray(sdata,n)
     char **sdata;
     int n;
{
  int i;
  float *data;
  
  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++)
    data[i] = atof(sdata[i]);
  
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_IARRAY_FROM_CARRAY                         */
/*                                                                           */
/*****************************************************************************/
int *get_iarray_from_carray(sdata,n)
     char **sdata;
     int n;
{
  int i;
  int *data;
  
  data = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++)
    data[i] = atoi(sdata[i]);
  
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MAKE_STRING_FROM_ITEMS                         */
/*                                                                           */
/*  Return the string consisting of the items separated by spaces.           */
/*                                                                           */
/*****************************************************************************/
char *make_string_from_items(cdata,n)
     char **cdata;
     int n;
{
  int i,j,k;
  int sn;
  char *s;

  sn = n; // n-1 spaces plus 1 null character
  for(i=0;i<n;i++)
    sn += strlen(cdata[i]);
  s = (char *)myalloc(sn*sizeof(char));

  k = 0;
  for(i=0;i<n;i++){
    for(j=0;j<strlen(cdata[i]);j++){
      s[k] = cdata[i][j];
      k += 1;
    }
    s[k] = ' ';
    k += 1;
  }
  s[k-1] = '\0'; // Overwrite the final ' '

  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                           CHECK_STRING_EXIT_CARRAY                        */
/*                                                                           */
/*****************************************************************************/
void check_string_exit_carray(s1,s2)
     char s1[],s2[];
{
  if (strcmp(s1,s2)!=0){
    printf("*** %s and %s not equal\n",s1,s2);
    exit_error("CHECK_STRING_EXIT_CARRAY","Strings not equal");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            BUBBLE_SORT_2D_CARRAY                          */
/*                                                                           */
/*****************************************************************************/
void bubble_sort_2d_carray(data,n)
     char **data;
     int n;
{
  int i;
  int done,bottom;
  char *t;
  
  bottom = n-1;
  done = 0;
  while (!done){
    done = 1;
    for (i=0;i<bottom;i++)
      if (strcmp(data[i],data[i+1]) > 0){
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
/*                          IS_UNIQUE_SORTED_2D_CARRAY                       */
/*                                                                           */
/*  Assuming the array is sorted, return 1 if the array has no duplicates.   */
/*                                                                           */
/*****************************************************************************/
int is_unique_sorted_2d_carray(c,n)
     char **c;
     int n;
{
  int i;

  for(i=1;i<n;i++)
    if (strcmp(c[i-1],c[i])==0)
      return 0;

  return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                            IS_UNIQUE_2D_CARRAY                            */
/*                                                                           */
/*  Return 1 if the array has no duplicates.                                 */
/*                                                                           */
/*****************************************************************************/
int is_unique_2d_carray(c,n)
     char **c;
     int n;
{
  int flag;
  char **cs;

  cs = copy_2d_carray(c,n);
  bubble_sort_2d_carray(cs,n);  // Wyeth this had been left out until Dec 2010
  flag = is_unique_sorted_2d_carray(cs,n);
  free_2d_carray(cs,n);
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                           CARRAY_GET_DUPLICATE_PTR                        */
/*                                                                           */
/*  Return a pointer to the first duplicate string, or NULL if none.         */
/*                                                                           */
/*****************************************************************************/
char *carray_get_duplicate_ptr(c,n)
     char **c;
     int n;
{
  int i,j;
  
  i = 0;
  while(i<(n-1)){
    j = i+1;
    while(j<n){
      if (strcmp(c[i],c[j])==0)
	return c[i];
      j += 1;
    }
    i += 1;
  }
  
  return NULL;
}
/**************************************-**************************************/
/*                                                                           */
/*                              SEARCH_2D_CARRAY                             */
/*                                                                           */
/*  Search for the string "s" in the list "c".  If found, return index,      */
/*  else return -1.                                                          */
/*                                                                           */
/*****************************************************************************/
int search_2d_carray(c,s,n)
     char **c,*s;
     int n;
{
  int k;
  int done;

  done = 0;
  k = 0;
  while((!done)&&(k<n)){
    if (strcmp(c[k],s)==0)
      done = 1;
    else
      k += 1;
  }
  if (done)
    return k;
  else
    return -1;
}
/**************************************-**************************************/
/*                                                                           */
/*                            SEARCH_2D_CARRAY_EXIT                          */
/*                                                                           */
/*  Search for the string "s" in the list "c".  If found, return index,      */
/*  else exit.                                                               */
/*                                                                           */
/*****************************************************************************/
int search_2d_carray_exit(c,s,n)
     char **c,*s;
     int n;
{
  int k;

  k = search_2d_carray(c,s,n);
  if (k < 0){
    printf("  string = %s\n",s);
    exit_error("SEARCH_2D_CARRAY_EXIT","String not found in list");
  }

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                             SEARCH_2D_CARRAY_FLAG                         */
/*                                                                           */
/*  Like "SEARCH_2D_CARRAY" but only look at elements for which              */
/*  "flag[i]" != 1.                                                          */
/*                                                                           */
/*****************************************************************************/
int search_2d_carray_flag(c,s,n,flag)
     char **c,*s;
     int n,*flag;
{
  int k;
  int done;

  done = 0;
  k = 0;
  while((!done)&&(k<n)){
    if (flag[k] != 0)
      if (strcmp(c[k],s)==0)
	done = 1;
      else
	k += 1;
    else
      k += 1;
  }
  if (done)
    return k;
  else
    return -1;
}
/**************************************-**************************************/
/*                                                                           */
/*                         UPDATE_KEY_VALUE_PAIR_CARRAY                      */
/*                                                                           */
/*  Add new values "newval" to a list of old values "oldval" based on the    */
/*  key "newname" and "oldname".   If any "newname" exists in the "oldname"  */
/*  list, then use the "newval" for that name.                               */
/*                                                                           */
/*  Return a new list formed by combining the old and new lists.             */
/*                                                                           */
/*****************************************************************************/
void update_key_value_pair_carray(oldname,oldval,nold,newname,newval,nnew,
				  rname,rval,rn)
     char **oldname,**oldval;
     int nold;
     char **newname,**newval;
     int nnew;
     char ***rname,***rval;
     int *rn;
{
  int i,k;
  int *flag,count,n;
  char **name,**val;

  count = 0;
  flag = get_zero_iarray(nold);
  for(i=0;i<nold;i++){
    flag[i] = search_2d_carray(newname,oldname[i],nnew);
    if (flag[i] < 0)
      count += 1;
  }
  /*printf("adding %d new names\n",count);*/
  n = count + nnew;
  name = (char **)myalloc(n*sizeof(char *));
  val = (char **)myalloc(n*sizeof(char *));
  for(i=0;i<nnew;i++){
    name[i] = strdup(newname[i]);
    val[i] = strdup(newval[i]);
  }

  k = nnew;
  for(i=0;i<nold;i++){
    if (flag[i] < 0){
      name[k] = strdup(oldname[i]);
      val[k] = strdup(oldval[i]);
      k += 1;
    }
  }
  myfree(flag);

  *rname = name; *rval = val; *rn = n;
}     
/**************************************-**************************************/
/*                                                                           */
/*                                PRINT_CARRAY                               */
/*                                                                           */
/*****************************************************************************/
void print_carray(cdata,n)
     char **cdata;
     int n;
{
  int i;

  for(i=0;i<n;i++)
    printf("%s\n",cdata[i]);
}
/**************************************-**************************************/
/*                                                                           */
/*                           PRINT_2D_POINTER_CARRAY                         */
/*                                                                           */
/*  Print the contents of the carray in matrix format.                       */
/*                                                                           */
/*****************************************************************************/
void print_2d_pointer_carray(cdata,m,n)
     char ***cdata;
     int m,n;
{
  int i,j;

  if ((m>0) && (n>0)){
    for(i=0;i<m;i++){
      for(j=0;j<n;j++)
	printf("%s ",cdata[i][j]);
      printf("\n");
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          PRINT_STRING_LIST_CARRAY                         */
/*                                                                           */
/*  Print the contents of the carray in matrix format.                       */
/*                                                                           */
/*****************************************************************************/
void print_string_list_carray(data,n,cnt)
     char ***data;
     int n,*cnt;
{
  int i,j;
  
  for(i=0;i<n;i++){
    for(j=0;j<cnt[i];j++)
      printf("%s ",data[i][j]);
    printf("\n");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          CARRAY_STRING_READ_COUNT                         */
/*                                                                           */
/*  Count how many strings can be read from this file using 'fscanf'.        */
/*                                                                           */
/*****************************************************************************/
int carray_string_read_count(infile)
     char *infile;
{
  FILE *fopen(),*fin;
  int n;
  char tstr[SLEN_MAX];

  if ((fin = fopen(infile,"r"))==NULL){
    printf("  filename = %s\n",infile);
    exit_error("CARRAY_STRING_READ_COUNT","Cannot open file");
  }

  n = 0;
  while (fscanf(fin,"%s",tstr) != EOF)
    n += 1;

  fclose(fin);

  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 READ_CARRAY                               */
/*                                                                           */
/*****************************************************************************/
void read_carray(infile,data,n)
     char infile[];
     char **data;
     int *n;
{
  FILE *fopen(),*fin;
  int i;

  if ((fin = fopen(infile,"r"))==NULL){
    printf("  filename = %s\n",infile);
    exit_error("READ_CARRAY","Cannot open file");
  }

  *n = 0;
  while (getc(fin) != EOF) // Must read into int for EOF condition
    *n += 1;
  fclose(fin);

  fin = fopen(infile,"r");
  *data = (char *)myalloc(*n*sizeof(char));
  for(i=0;i<*n;i++)
    (*data)[i] = (char)getc(fin);
  fclose(fin);
}
/**************************************-**************************************/
/*                                                                           */
/*                               READ_2D_CARRAY                              */
/*                                                                           */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - If the file cannot be opened, return -1 in *rn.                        */
/*  - If 'remove_nl' is 1, remove the '\n' from the end of strings.          */
/*  - WARNING:  this has a limited line length, SLINE_MAX, and will exit     */
/*    if that is exceeded.                                                   */
/*                                                                           */
/*****************************************************************************/
void read_2d_carray(infile,remove_nl,blank_flag,rdata,rn)
     char infile[];
     int remove_nl;
     int blank_flag;  // 0-drop blank lines, 1-keep blank lines
     char ***rdata;
     int *rn;
{
  FILE *fopen(),*fin;
  int k;
  int n,ns;
  char **data,**slist,line[SLINE_MAX];

  if ((fin = fopen(infile,"r"))==NULL){
    *rdata = NULL;
    *rn = -1;
    return;
  }

  //printf("HERE WYETH reading %s\n",infile);
  n = 0;
  while((fgets(line,SLINE_MAX,fin) != NULL)){
    if ((strlen(line) >= SLINE_MAX-2))
      exit_error("READ_2D_CARRAY","LINE TOO LONG");
    get_items_from_string(line,&slist,&ns);
    if (ns > 0){
      n += 1;
      free_2d_carray(slist,ns);
    }else if (blank_flag == 1){
      n += 1;                    // Count this blank line
    }
  }
  fclose(fin);
  //printf("HERE WYETH done\n");

  if ((fin = fopen(infile,"r"))==NULL){
    *rdata = NULL;
    *rn = -1;
    return;
  }
  data = (char **)myalloc(n*sizeof(char *));
  k = 0;
  //printf("HERE WYETH  001\n");
  while((fgets(line,SLINE_MAX,fin) != NULL)){
    get_items_from_string(line,&slist,&ns);
    if ((ns > 0) || (blank_flag == 1)){
      if (remove_nl){

	//
	//  WYETH - this modified to handle the value '13' that appears
	//  in some ASCII files near the end of line.  This shows up as
	//  ^M in some editors.
	//

	if (strlen(line) > 0){  // If there is at least one character
	  if (line[strlen(line)-1]=='\n')  // Newline
	    line[strlen(line)-1] = '\0';   // remove it
	  if (strlen(line) > 0){  // If there is at least one character
	    if (line[strlen(line)-1]==(char)13)  // Carriage return
	      line[strlen(line)-1] = '\0';   // remove it
	  }
	}
      }
      data[k] = strdup(line);
      k += 1;
      if (ns > 0)
	free_2d_carray(slist,ns);
    }
  }
  //printf("HERE WYETH  002\n");
  fclose(fin);
  *rdata = data; *rn = k;
}
/**************************************-**************************************/
/*                                                                           */
/*                            SEARCH_FILE_2D_CARRAY                          */
/*                                                                           */
/*  Return 1 if 's' is one of the lines of 'infile', 0 otherwise.            */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - File is a list of strings separated by carriage returns.               */
/*  - Strings may contain spaces.                                            */
/*  - Carriage return is ignored in matching process.                        */
/*                                                                           */
/*****************************************************************************/
int search_file_2d_carray(infile,s)
     char infile[];
     char *s;
{
  FILE *fopen(),*fin;
  int ns;
  char line[SLINE_MAX];
  
  if ((fin = fopen(infile,"r"))==NULL)
    exit_error("SEARCH_FILE_2D_CARRAY","Cannot open file");

  while((fgets(line,SLINE_MAX,fin) != NULL)){
    if ((strlen(line) >= SLINE_MAX-2))
      exit_error("SEARCH_FILE_2D_CARRAY","LINE TOO LONG");
    ns = strlen(line);
    if (ns > 0){
      if (line[ns-1]=='\n')
	line[ns-1] = '\0';
    }
    ns = strlen(line);
    if (strcmp(s,line)==0)
      return 1;
  }
  fclose(fin);
  return 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                                WRITE_CARRAY                               */
/*                                                                           */
/*****************************************************************************/
void write_carray(outfile,data,n)
     char outfile[];
     char *data;
     int n;
{
  FILE *fopen(),*fout;
  int i;

  fout = fopen(outfile,"w");
  for(i=0;i<n;i++)
    fprintf(fout,"%c",data[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                             WRITE_CARRAY_SEGMENT                          */
/*                                                                           */
/*****************************************************************************/
void write_carray_segment(outfile,data,n,start,period)
     char outfile[];
     char *data;
     int n,start,period;
{
  FILE *fopen(),*fout;
  int i;

  fout = fopen(outfile,"w");
  for(i=start;i<(start+period);i++)
    fprintf(fout,"%c",data[i]);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                           ENCODE_SLIST_AS_CARRAY                          */
/*                                                                           */
/*****************************************************************************/
void encode_slist_as_carray(slist,ns,rc,rn)
     char **slist;
     int ns;
     char **rc;
     int *rn;
{
  int i,j,k;
  int n,tn;
  char *t,*s,tstr[SLEN];

  exit_error("ENCODE_SLIST_AS_CARRAY","WYETH - NEVER TESTED");

  /*** Determine length of carray ***/
  n = 0;
  for(i=0;i<ns;i++){
    n += strlen(slist[i]);
    n += 1; /* For terminator chars after each string */
  }
  sprintf(tstr,"%d",ns);
  n += strlen(tstr) + 1;  /* length of header number and '\0' */

  s = (char *)myalloc(n*sizeof(char));
  sprintf(s,"%d",ns); /* Number of strings followed by '\0' */

  k = strlen(s) + 1;
  for(i=0;i<ns;i++){ /*** Store each string followed by '\0' ***/
    t = slist[i];
    tn = strlen(t);
    for(j=0;j<tn;j++){
      s[k] = t[j];
      k += 1;
    }
    s[k] = '\0';
    k += 1;
  }

  *rc = s;
  *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                              CARRAY_NCHAR_READ                            */
/*                                                                           */
/*****************************************************************************/
void carray_nchar_read(p,fin,revflag)
     char **p;
     FILE *fin;
     int revflag;
{
  int n;
  int nread;

  nread = revflag_fread((char *)&n,sizeof(int),1,fin,revflag);
  if (nread == 0)
    exit_error("CARRAY_NCHAR_READ","Error");
  *p = (char *)myalloc((n+1)*sizeof(char));
  nread = fread((char *)(*p),sizeof(char),n,fin);
  (*p)[n] = '\0';
}
/**************************************-**************************************/
/*                                                                           */
/*                              CARRAY_NCHAR_WRITE                           */
/*                                                                           */
/*****************************************************************************/
void carray_nchar_write(p,fout)
     char *p;
     FILE *fout;
{
  int n;

  n = strlen(p);
  fwrite((char *)&n,sizeof(int),1,fout);
  fwrite((char *)p,sizeof(char),n,fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                   GET_EXCLUSIVE_INDEX_COMPARE_2D_CARRAYS                  */
/*                                                                           */
/*  Get an index into "cdata1" for each item that is not in "cdata2".        */
/*                                                                           */
/*****************************************************************************/
void get_exclusive_index_compare_2d_carrays(cdata1,n1,cdata2,n2,rindex,rn)
     char **cdata1;
     int n1;
     char **cdata2;
     int n2,**rindex,*rn;
{
  int i,j,k;
  int *ndx,*t,done;

  k = 0; /* Count the number of "new" items. */
  ndx = (int *)myalloc(n1*sizeof(int));
  for(i=0;i<n1;i++){ /*** For each string in the first array. ***/
    done = 0;
    j = 0;
    while(!done){ /*** Search for that string in the second array. ***/
      if (j >= n2)
	done = 1;
      else if (strcmp(cdata1[i],cdata2[j])==0)
	done = 1;
      else
	j += 1;
    }
    if (j >= n2){ /*** This string was not found, so it is new. ***/
      ndx[k] = i;
      k += 1;
    }
  }

  t = (int *)myalloc(k*sizeof(int));
  for(i=0;i<k;i++)
    t[i] = ndx[i];

  *rindex = t;
  *rn = k;
}
/**************************************-**************************************/
/*                                                                           */
/*                     GET_SHARED_INDEX_COMPARE_2D_CARRAYS                   */
/*                                                                           */
/*  Create an index containing "n1" entries which points to the position of  */
/*  each value of "cdata1" in "cdata2", or -1 if the value is not in         */
/*  "cdata2".                                                                */
/*                                                                           */
/*****************************************************************************/
int *get_shared_index_compare_2d_carrays(cdata1,n1,cdata2,n2)
     char **cdata1;
     int n1;
     char **cdata2;
     int n2;
{
  int i;
  int *ndx;

  ndx = (int *)myalloc(n1*sizeof(int));
  for(i=0;i<n1;i++) // For each string in the first array.
    ndx[i] = search_2d_carray(cdata2,cdata1[i],n2);
  return ndx;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_COLUMN_3D_CARRAY                           */
/*                                                                           */
/*  Get the "k" column from the 2D array of strings.                         */
/*                                                                           */
/*****************************************************************************/
char **get_column_3d_carray(cdata,k,nrow,ncol)
     char ***cdata;
     int k,nrow,ncol;
{
  int i;
  char **data;

  data = (char **)myalloc(nrow*sizeof(char *));
  for(i=0;i<nrow;i++)
    data[i] = strdup(cdata[i][k]);
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                       GET_UNIQUE_SLIST_FROM_3D_CARRAY                     */
/*                                                                           */
/*  Get the "k" column from the 2D array of strings.                         */
/*                                                                           */
/*****************************************************************************/
char **get_unique_slist_from_3d_carray(cdata,xn,yn,rn)
     char ***cdata;
     int xn,yn;
     int *rn;          /* Number of unique elements */
{
  int i,j,k;
  int un,uflag;
  char **u,**t,*tt;

  if ((xn*yn) < 1)
    exit_error("GET_UNIQUE_STRINGS_FROM_3D_CARRAY","No elements");

  un = 1;
  t = (char **)myalloc(xn*yn*sizeof(char *));
  t[0] = cdata[0][0];

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      tt = cdata[i][j];  /* Current string */
      uflag = 1;         /* Assume unique ... */
      for(k=0;k<un;k++){
	if (strcmp(tt,t[k])==0)  /* until it is found in the list */
	  uflag = 0;
      }
      if (uflag){     /* Add this new element to the unique list */
	/*printf("unique:  %s\n",tt);*/
	t[un] = tt;
	un += 1;
      }
    }
  }

  u = (char **)myalloc(un*sizeof(char *));
  for(i=0;i<un;i++)
    u[i] = strdup(t[i]);

  myfree(t);

  *rn = un;
  return u;
}
/**************************************-**************************************/
/*                                                                           */
/*                         CARRAY_UTIL_PARSE_UNDERSCORE                      */
/*                                                                           */
/*  Return an array of ints for a string as follows:  <int>_<int>_..._<int>  */
/*                                                                           */
/*****************************************************************************/
int *carray_util_parse_underscore(s,rn)
     char *s;   // String to parse
     int *rn;   // Number of values in returned arary
{
  int i,k;
  int n,cnt,*ival;
  char *a,*tstr;

  tstr = strdup(s);
  n = strlen(s);

  cnt = 0;
  for(i=0;i<n;i++){
    if (tstr[i] == '_')
      cnt += 1;
  }

  ival = (int *)myalloc((cnt+1)*sizeof(int));

  k = 0;
  a = tstr;
  for(i=0;i<n;i++){
    if (tstr[i] == '_'){
      tstr[i] = ' ';
      sscanf(a,"%d",&(ival[k]));
      //printf("____a = %s\n",a);
      k += 1;
      a = &(tstr[i+1]);
    }
  }
  sscanf(a,"%d",&(ival[k]));
      
  *rn = cnt + 1;
  return ival;
}
/**************************************-**************************************/
/*                                                                           */
/*                      CARRAY_UTIL_UNDERSCORE_INT_FINAL                     */
/*                                                                           */
/*  Return the integer at the end of the string "..._<int>"                  */
/*                                                                           */
/*****************************************************************************/
int carray_util_underscore_int_final(s)
     char *s;   // String to parse
{
  int i;
  int n,ui,ival;
  char *a;

  n = strlen(s);

  ui = -1; // Get the index of the final underscore
  for(i=0;i<n;i++){
    if (s[i] == '_')
      ui = i;
  }

  if (ui == -1){
    printf("   string:  %s\n",s);
    exit_error("CARRAY_UTIL_UNDERSCORE_INT_FINAL",
	       "No underscore found in string");
  }else if (ui > (n-2)){
    printf("   string:  %s\n",s);
    exit_error("CARRAY_UTIL_UNDERSCORE_INT_FINAL",
	       "Underscore index too large");
  }

  a = &(s[ui+1]);
  sscanf(a,"%d",&ival);

  return ival;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_UNDERSCORE_EQUAL_FLOAT                        */
/*                                                                           */
/*  'Underscore equal' is a type of string that I will use for names which   */
/*  carry values.  The format is 'name1=val1_name2=val2_..._namen=valn'      */
/*                                                                           */
/*****************************************************************************/
int carray_get_float_val_underscore_eq(s,name,rval)
     char s[],name[];
     float *rval;
{
  int i,j,k;
  int slen,done;
  char tname[SLEN],vstr[SLEN];

  /*printf("s= %s\n",s);*/

  slen = strlen(s);
  j = 0;
  k = 0;  /* Pointer to current letter in 's' */
  done = 0;
  while(!done){
    if (k >=slen){
      return 0;
    }else{
      if (s[k] == '='){
	tname[j] = '\0'; /* Put a null at the end of the tname */
	/*printf("TNAME = %s\n",tname);*/
	if (strcmp(name,tname)==0){ /* We've found it, now get value */
	  i = 0;
	  k += 1; /* Skip over the '=' */
	  while(!done){
	    if (k >= slen){ /* End of the input string */
	      done = 1;
	      vstr[i] = '\0';
	      /*printf("Vstr= %s\n",vstr);*/
	    }else if (s[k] == '_'){ /* End of the value */
	      done = 1;
	      vstr[i] = '\0';
	      /*printf("Vstr= %s\n",vstr);*/
	    }else{
	      vstr[i] = s[k];
	      i += 1;
	      k += 1;
	    }
	  }
	  done = 1;
	}else{ /* We haven't found it, advance to next '_', start again */
	  while(s[k] != '_'){
	    k += 1;
	    if (k >= slen)
	      return 0; /* Not found */
	  }
	  k += 1; /* Skip over the '_' */
	  j = 0; /* Reset 'tname' letter counter */
	}
      }else{
	tname[j] = s[k];
	k += 1;
	j += 1;
      }
    }
  }
  *rval = atof(vstr);
  return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                       CARRAY_ITEM_LIST_REMOVE_COMMENT                     */
/*                                                                           */
/*  Remove anything at or after the first item beginning with '#', and       */
/*  reformat as an item list separated by single spaces.                     */
/*                                                                           */
/*****************************************************************************/
char *carray_item_list_remove_comment(s)
     char s[];
{
  int k;
  int sn,done;
  char **slist,*ss;

  get_items_from_string(s,&slist,&sn);
  
  if (sn > 1){  // There is more than one string
    k = 0;
    done = 0;
    while(!done){
      if (k >= sn)
	done = 1;
      else if (slist[k][0] == '#')
	done = 1;
      else
	k++;
    }
    if (k < sn){  // If there are some comments
      ss = make_string_from_items(slist,k);
    }else
      ss = strdup(s);
  }else
    ss = strdup(s);

  free_2d_carray(slist,sn);
  
  return ss;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MSTREL_GET_INIT                              */
/*                                                                           */
/*****************************************************************************/
struct math_strel *mstrel_get_init()
{
  struct math_strel *ms;

  ms = (struct math_strel *)myalloc(sizeof(struct math_strel));
  ms->etype = NULL;
  ms->fval = 0.0;
  ms->ival = 0;
  ms->sval = NULL;
  ms->prev = ms->next = NULL;

  return ms;
}
/**************************************-**************************************/
/*                                                                           */
/*                                MSTREL_FREE                                */
/*                                                                           */
/*****************************************************************************/
void mstrel_free(ms)
     struct math_strel *ms;
{
  if (ms->etype != NULL)
    myfree(ms->etype);
  if (ms->sval != NULL)
    myfree(ms->sval);

  myfree(ms);
}
/**************************************-**************************************/
/*                                                                           */
/*                                MSTREL_PUSH                                */
/*                                                                           */
/*****************************************************************************/
void mstrel_push(ps,e)
     struct math_strel **ps;  /* Pointer to top element */
     struct math_strel *e;    /* Element to push */
{
  if (*ps == NULL){  /* Stack is empty, start w/ this element */
    *ps = e;
  }else{
    (*ps)->prev = e;
    e->next = *ps;
    *ps = e;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                 MSTREL_POP                                */
/*                                                                           */
/*****************************************************************************/
struct math_strel *mstrel_pop(ps)
     struct math_strel **ps;  /* Pointer to top element */
{
  struct math_strel *e;

  if (*ps == NULL)  /* Stack is empty, start w/ this element */
    e = NULL;
  else{
    e = (*ps);             /* return the top element */

    /*printf("  e  %s %s\n",e->etype,e->sval);*/

    if (e->next != NULL) 
      e->next->prev = NULL;  /* new top element will have no 'prev' */
    *ps = e->next;        /* second element is top */

    e->next = e->prev = NULL;  /* A popped element is isolated */
  }

  return e;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MSTREL_FREE_STACK                             */
/*                                                                           */
/*****************************************************************************/
void mstrel_free_stack(ps)
     struct math_strel **ps;  /* pointer to stack */
{
  struct math_strel *e;

  e = mstrel_pop(ps);
  while(e != NULL){
    mstrel_free(e);
    e = mstrel_pop(ps);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MSTREL_STACK_SIZE                             */
/*                                                                           */
/*****************************************************************************/
int mstrel_stack_size(s)
     struct math_strel *s;  /* top element */
{
  int n;
  struct math_strel *t;

  t = s;
  n = 0;
  while(t != NULL){
    t = t->next;
    n += 1;
  }

  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MSTREL_ETYPE_N                               */
/*                                                                           */
/*  Return '1' if the 'n'th element is of type 'etype'.                      */
/*                                                                           */
/*****************************************************************************/
int mstrel_etype_n(s,n,etype)
     struct math_strel *s;  /* top element */
     int n;
     char *etype;
{
  int k;
  struct math_strel *t;

  t = s;
  k = 1;
  while(t != NULL){
    if (k == n){
      if (strcmp(t->etype,etype)==0)
	return 1;
      else
	return 0;
    }
    t = t->next;
    k += 1;
  }
  
  return 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MSTREL_STACK_PRINT                            */
/*                                                                           */
/*****************************************************************************/
void mstrel_stack_print(s)
     struct math_strel *s;  /* top element */
{
  struct math_strel *t;

  printf("--------------STACK------------\n");
  t = s;
  while(t != NULL){
    printf("  %s %s %f\n",t->etype,t->sval,t->fval);
    t = t->next;
  }
  printf("-------------------------------\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                          MSTREL_PAREN_MATCH_REMOVE                        */
/*                                                                           */
/*  Look at second element back, if '(' then remove it.                      */
/*                                                                           */
/*****************************************************************************/
int mstrel_paren_match_remove(ps)
     struct math_strel **ps;  /* Pointer to pointer to first stack element */
{
  int flag;
  struct math_strel *e1,*e2;

  flag = 0;
  if (mstrel_etype_n(*ps,2,"group")){
    if (strcmp((*ps)->next->sval,"(")==0){
      e1 = mstrel_pop(ps);
      e2 = mstrel_pop(ps);
      myfree(e2);
      mstrel_push(ps,e1);
      flag = 1;
    }
  }
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                                MSTREL_CALC                                */
/*                                                                           */
/*  Do the calculation using the top 3 elements of the stack, and put the    */
/*  result on the top of the stack.                                          */
/*                                                                           */
/*****************************************************************************/
void mstrel_calc(ps,errstr)
     struct math_strel **ps;  /* Pointer to pointer to first stack element */
     char *errstr;
{
  float v1,v2,a;
  struct math_strel *e1,*e2,*e3,*e;

  e1 = mstrel_pop(ps);
  if (e1 == NULL)
    strcpy(errstr,"MSTREL_CALC - No elements on stack");
  if (strcmp(e1->etype,"fval")!=0)
    strcpy(errstr,"MSTREL_CALC - Top element should be value");
  v2 = e1->fval;

  e2 = mstrel_pop(ps);
  if (e2 == NULL)
    strcpy(errstr,"MSTREL_CALC - No operator element");
  if ((strcmp(e2->etype,"op1")!=0) && (strcmp(e2->etype,"op2")!=0))
    strcpy(errstr,"MSTREL_CALC - Second element should be operator");
  
  e3 = mstrel_pop(ps);
  if (e3 == NULL)
    strcpy(errstr,"MSTREL_CALC - No third element");
  if (strcmp(e3->etype,"fval")!=0)
    strcpy(errstr,"MSTREL_CALC - Third element should be value");
  v1 = e3->fval;
  
  /*printf("CALC:  v1 v2 op =  %f %f %s\n",v1,v2,e2->sval);*/

  if (strcmp(e2->sval,"+")==0)
    a = v1 + v2;
  else if (strcmp(e2->sval,"-")==0)
    a = v1 - v2;
  else if (strcmp(e2->sval,"*")==0)
    a = v1 * v2;
  else if (strcmp(e2->sval,"/")==0)
    a = v1 / v2;

  e = mstrel_get_init();
  e->etype = strdup("fval");
  e->fval = a;
  mstrel_push(ps,e);

  mstrel_free(e1);
  mstrel_free(e2);
  mstrel_free(e3);
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_NEXT_MATH_STRING_ELEMENT                      */
/*                                                                           */
/*****************************************************************************/
struct math_strel *get_next_math_string_element(pmstr)
     char **pmstr;
{
  int k;
  char c,*t,tstr[SLEN];
  struct math_strel *ms;

  ms = mstrel_get_init();

  t = *pmstr;  /* t points to beginning of string */
  while(t[0] == ' '){
    t += 1;
  }
  c = t[0];

  if ((c == '(') || (c == ')')){
    ms->etype = strdup("group");
    sprintf(tstr,"%c",t[0]);
    ms->sval = strdup(tstr);
    *pmstr = t+1;
  }else if ((c == '+') || (c == '-')){
    ms->etype = strdup("op1");
    sprintf(tstr,"%c",t[0]);
    ms->sval = strdup(tstr);
    *pmstr = t+1;
  }else if ((c == '*') || (c == '/')){
    ms->etype = strdup("op2");
    sprintf(tstr,"%c",t[0]);
    ms->sval = strdup(tstr);
    *pmstr = t+1;
  }else if (is_digit_char(c)){
    ms->etype = strdup("fval");

    k = 0;
    while (is_digit_char(*t) || (*t == '.')){
      tstr[k] = *t;
      t += 1;
      k += 1;
    }
    tstr[k] = '\0';
    ms->sval = strdup(tstr);
    ms->fval = atof(ms->sval);
    *pmstr = t;
  }else if (c == '\0'){
    return NULL;
  }else{
    ms->etype = strdup("var");
    
    k = 0;
    while (is_name_char(*t)){
      tstr[k] = *t;
      t += 1;
      k += 1;
    }
    tstr[k] = '\0';
    ms->sval = strdup(tstr);
    *pmstr = t;
  }
  return ms;
}
/**************************************-**************************************/
/*                                                                           */
/*                                MSTR_EVALUATE                              */
/*                                                                           */
/*  Evaluate a math string.                                                  */
/*                                                                           */
/*****************************************************************************/
float mstr_evaluate(mstr,nv,vname,vval,rer)
     char *mstr;
     int nv;         /* Number of variables */
     char **vname;   /* Names of variables [nv] */
     char **vval;    /* Values of variables [nv] */
     char **rer;     /* Return error string (if not NULL, caller must free */
{
  int k;
  int done,pflag,ssize;
  float x;
  char errstr[SLEN];
  struct math_strel *ms,*ms1,*mstack;

  pflag = 0;
  mstack = NULL;
  errstr[0] = '\0';   /* Empty the error string */
  *rer = NULL;        /* Set return error string to NULL */

  done = 0;
  while(!done){
    ms = get_next_math_string_element(&mstr);
    if (ms == NULL){
      done = 1;
      while (mstrel_stack_size(mstack) >= 3){
	if (pflag)
	  printf("  calling calc:  end of input reached\n");
	mstrel_calc(&mstack,errstr);
	if (errstr[0] != '\0'){
	  *rer = strdup(errstr);
	  return 0.0;
	}
      }
      ssize = mstrel_stack_size(mstack);
      if (ssize != 1){
	sprintf(errstr,"MSTR_EVALUATE - Final stack size %d",ssize);
	*rer = strdup(errstr);
	return 0.0;
      }
      x = mstack->fval;
      mstrel_free_stack(&mstack);
    }else{
      if (strcmp(ms->etype,"op1")==0){
	if (!mstrel_etype_n(mstack,1,"fval")){
	  if (strcmp(ms->sval,"+")==0)
	    mstrel_free(ms); /* Throw away a unitary plus sign */
	  else{  /* Assume unitary minus, put "-1.0 *" on stack */
	    ms1 = mstrel_get_init();
	    ms1->etype = strdup("fval");
	    ms1->fval = -1.0;
	    ms1->sval = strdup("-1.0");
	    mstrel_push(&mstack,ms1);

	    myfree(ms->etype);
	    myfree(ms->sval);
	    ms->etype = strdup("op2");
	    ms->sval = strdup("*");
	    mstrel_push(&mstack,ms);
	  }
	}else{
	  if (mstrel_stack_size(mstack) >= 3){
	    if (!mstrel_etype_n(mstack,2,"group")){
	      if (pflag)
		printf("  calling calc - op1 found, and stack >= 3\n");
	      mstrel_calc(&mstack,errstr);
	      if (errstr[0] != '\0'){
		*rer = strdup(errstr);
		return 0.0;
	      }
	    }
	  }
	  mstrel_push(&mstack,ms);
	}
      }else if (strcmp(ms->etype,"op2")==0){
	mstrel_push(&mstack,ms);
      }else if (strcmp(ms->etype,"group")==0){
	if (strcmp(ms->sval,"(")==0)
	  mstrel_push(&mstack,ms);
	else if (strcmp(ms->sval,")")==0){
	  if (!mstrel_paren_match_remove(&mstack)){
	    if (mstrel_stack_size(mstack) >= 3){
	      if (pflag)
		printf("  calling calc - parenthetical end\n");
	      mstrel_calc(&mstack,errstr);
	      if (errstr[0] != '\0'){
		*rer = strdup(errstr);
		return 0.0;
	      }
	      mstrel_paren_match_remove(&mstack);
	    }
	  }
	  if (mstrel_etype_n(mstack,2,"op2")){
	    if (pflag)
	      printf("  calling calc - op2 before parenthetical\n");
	    mstrel_calc(&mstack,errstr);
	    if (errstr[0] != '\0'){
	      *rer = strdup(errstr);
	      return 0.0;
	    }
	  }
	}
      }else{
	/**** FOR NOW IGNORING VARs ***/
	/*** Must be a numeric value ***/
	mstrel_push(&mstack,ms);

	if (!is_float_string(ms->sval)){
	  k = search_2d_carray(vname,ms->sval,nv);
	  if (k < 0){
	    sprintf(errstr,"MSTR_EVALUATE - Variable name '%s' not found",
		    ms->sval);
	    *rer = strdup(errstr);
	    return 0.0;
	  }
	  ms->fval = atof(vval[k]);
	  myfree(ms->etype);
	  ms->etype = strdup("fval");
	  myfree(ms->sval);
	  ms->sval = strdup(vval[k]);
	}

	
	if (mstrel_stack_size(mstack) >= 3){
	  if (mstrel_etype_n(mstack,2,"op2")){
	    if (pflag)
	      printf("  calling calc - value pushed and previous is op2\n");
	    mstrel_calc(&mstack,errstr);
	    if (errstr[0] != '\0'){
	      *rer = strdup(errstr);
	      return 0.0;
	    }
	  }
	}
      }
    }
    if (pflag)
      mstrel_stack_print(mstack);
  }

  return x;
}
