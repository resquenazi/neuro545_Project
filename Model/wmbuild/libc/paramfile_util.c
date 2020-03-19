/*****************************************************************************/
/*                                                                           */
/*  paramfile_util.c                                                         */
/*  Wyeth Bair                                                               */
/*                                                                           */
/*  Routines for reading a list of parameter names and values.               */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "my_util.h"
#include "carray_util.h"
#include "data_util.h"
#include "paramfile.h"

struct onode_parse_state {
  int j;         // Line number
  int i;         // Item number, 0..ns-1
                 // XML - current char index
  char **slist;  // Current list of items [ns]
                 // XML - not used
  int ns;        // Number of items
};

static int totocnt = 0;

/**************************************-**************************************/
/*                                                                           */
/*                           PARAMFILE_PPL_GET_INIT                          */
/*                                                                           */
/*****************************************************************************/
struct param_pair_list *paramfile_ppl_get_init()
{
  struct param_pair_list *p;

  p = (struct param_pair_list *)myalloc(sizeof(struct param_pair_list));

  p->n = 0;
  p->p = NULL;

  p->nc = 0;
  p->c = NULL;

  p->ns = 0;
  p->s = NULL;

  return p;
}
/**************************************-**************************************/
/*                                                                           */
/*                             PARAMFILE_FREE_PAIR                           */
/*                                                                           */
/*****************************************************************************/
void paramfile_free_pair(pp)
     struct param_pair *pp;
{
  myfree(pp->name);
  myfree(pp->value);
  myfree(pp);
}
/**************************************-**************************************/
/*                                                                           */
/*                             PARAMFILE_FREE_PPL                            */
/*                                                                           */
/*****************************************************************************/
void paramfile_free_ppl(ppl)
     struct param_pair_list *ppl;
{
  int i;
  struct param_pair *pp,*next;

  pp = ppl->p;
  while(pp != NULL){
    next = pp->next;
    paramfile_free_pair(pp);
    pp = next;
  }

  pp = ppl->c;
  while(pp != NULL){
    next = pp->next;
    paramfile_free_pair(pp);
    pp = next;
  }

  if (ppl->s != NULL){
    for(i=0;i<ppl->ns;i++)
      myfree(ppl->s[i]);
    myfree(ppl->s);
  }

  myfree(ppl);
}
/**************************************-**************************************/
/*                                                                           */
/*                            PARAMFILE_PPL_TO_CARRAY                        */
/*                                                                           */
/*  Convert a ppl to a carray for transmission.                              */
/*                                                                           */
/*****************************************************************************/
void paramfile_ppl_to_carray(ppl,rs,rn)
     struct param_pair_list *ppl;
     char **rs;
     int *rn;
{
  int i,j,k;
  int n,np,nc,ns,nn;
  char *s,*t,tstr[SLEN];
  struct param_pair *p;

  //  ppl->p
  np = ppl->n;
  n = 0;
  p = ppl->p;
  for(i=0;i<np;i++){
    n += strlen(p->name);
    n += strlen(p->value);
    n += 2; // For terminator chars after each string
    p = p->next;
  }

  //  ppl->c
  nc = ppl->nc;
  p = ppl->c;
  for(i=0;i<nc;i++){
    n += strlen(p->name);
    n += strlen(p->value);
    n += 2; // For terminator chars after each string
    p = p->next;
  }

  //  ppl->s
  ns = ppl->ns;
  for(i=0;i<ns;i++){
    n += strlen(ppl->s[i]);
    n += 1; // For terminator chars after each string
  }

  sprintf(tstr,"%d %d %d",np,nc,ns);
  n += strlen(tstr) + 1;  // length of header number and '\0'

  s = (char *)myalloc(n*sizeof(char));
  sprintf(s,"%d %d %d",np,nc,ns); // Start with # of pairs followed by '\0'

  p = ppl->p;
  k = strlen(s) + 1;
  for(i=0;i<np;i++){
    // Store the 'name' followed by '\0'
    t = p->name;
    nn = strlen(t);
    for(j=0;j<nn;j++){
      s[k] = t[j];
      k += 1;
    }
    s[k] = '\0';
    k += 1;

    // Store the 'value' followed by '\0'
    t = p->value;
    nn = strlen(t);
    for(j=0;j<nn;j++){
      s[k] = t[j];
      k += 1;
    }
    s[k] = '\0';
    k += 1;

    p = p->next;
  }

  p = ppl->c;
  for(i=0;i<nc;i++){
    // Store the 'name' followed by '\0'
    t = p->name;
    nn = strlen(t);
    for(j=0;j<nn;j++){
      s[k] = t[j];
      k += 1;
    }
    s[k] = '\0';
    k += 1;

    // Store the 'value' followed by '\0'
    t = p->value;
    nn = strlen(t);
    for(j=0;j<nn;j++){
      s[k] = t[j];
      k += 1;
    }
    s[k] = '\0';
    k += 1;

    p = p->next;
  }

  for(i=0;i<ns;i++){
    // Store the 'name' followed by '\0'
    t = ppl->s[i];
    nn = strlen(t);
    for(j=0;j<nn;j++){
      s[k] = t[j];
      k += 1;
    }
    s[k] = '\0';
    k += 1;
  }

  *rs = s; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                            PARAMFILE_CARRAY_TO_PPL                        */
/*                                                                           */
/*  Convert a ppl to a carray for transmission.                              */
/*                                                                           */
/*****************************************************************************/
struct param_pair_list *paramfile_carray_to_ppl(c,cn)
     char *c;   // [cn] array of characters
     int cn;    // WYETH - Bug found June 2, 2016, this had been '*cn'
{
  int i;
  int np,nc,ns;
  char *t;
  struct param_pair *p;
  struct param_pair_list *ppl;

  ppl = paramfile_ppl_get_init();

  sscanf(c,"%d %d %d",&np,&nc,&ns);
  ppl->n  = np;
  ppl->nc = nc;
  ppl->ns = ns;

  p = NULL; // Avoid -Wall warning

  t = strchr(c,'\0') + 1;
  for(i=0;i<np;i++){
    if (i==0){
      ppl->p = (struct param_pair *)myalloc(sizeof(struct param_pair));
      p = ppl->p;
    }else{
      p->next = (struct param_pair *)myalloc(sizeof(struct param_pair));
      p = p->next;
    }
    p->name = strdup(t);

    t = strchr(t,'\0') + 1;
    p->value = strdup(t);

    t = strchr(t,'\0') + 1;
    p->next = NULL;
  }

  for(i=0;i<nc;i++){
    if (i==0){
      ppl->c = (struct param_pair *)myalloc(sizeof(struct param_pair));
      p = ppl->c;
    }else{
      p->next = (struct param_pair *)myalloc(sizeof(struct param_pair));
      p = p->next;
    }
    p->name = strdup(t);

    t = strchr(t,'\0') + 1;
    p->value = strdup(t);

    t = strchr(t,'\0') + 1;
    p->next = NULL;
  }

  if (ns > 0){
    ppl->s = (char **)myalloc(ns*sizeof(char *));
    for(i=0;i<ns;i++){
      ppl->s[i] = strdup(t);
      t = strchr(t,'\0') + 1;
    }
  }

  //
  //  Verify that the number of characters read from 'c' is 'cn'
  //
  if (cn != (int)((long)t - (long)c))
    exit_error("PARAMFILE_CARRAY_TO_PPL","Character count error");
  
  return ppl;
}
/**************************************-**************************************/
/*                                                                           */
/*                              PARAMFILE_PRINT_PPL                          */
/*                                                                           */
/*****************************************************************************/
void paramfile_print_ppl(ppl)
     struct param_pair_list *ppl;
{
  int i;
  struct param_pair *p;

  p = ppl->p;
  for(i=0;i<ppl->n;i++){
    printf("%s %s\n",p->name,p->value);
    p = p->next;
  }
  p = ppl->c;
  for(i=0;i<ppl->nc;i++){
    printf("%s %s\n",p->name,p->value);
    p = p->next;
  }

  if (ppl->ns > 0){
    printf("INLINE DATA:  %d lines\n",ppl->ns);
    for(i=0;i<ppl->ns;i++)
      printf("%s\n",ppl->s[i]);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             PARAMFILE_TEST_PARAM                          */
/*                                                                           */
/*  Return 1 if the parameter is in the 'ppl', 0 otherwise.                  */
/*                                                                           */
/*****************************************************************************/
int paramfile_test_param(ppl,name)
     struct param_pair_list *ppl;
     char *name;
{
  int k;
  char *tname,*tc;
  struct param_pair *tp;

  if (ppl == NULL){ // WYETH this added in case onode model calls
    printf("*** PARAMFILE_TEST_PARAM:  Warning, ppl is null.  Test %s\n",
	   name);
    return 0;
  }

  //
  // Check for '.' notation, such as "VARGEN_seed.6"
  //
  tname = strdup(name);    // Make a copy so we can modify
  tc = strchr(tname,'.');  // Search for '.' in the name
  if (tc != NULL){
    tc[0] = '\0';          // Change the '.' into a string terminator
  }

  k = 0;
  if (ppl->p != NULL){
    tp = ppl->p;
    while((strcmp(tp->name,tname)!=0) && (tp->next != NULL))
      tp = tp->next;
    if (strcmp(tp->name,tname)==0)
      k = 1;
  }

  // Free the copy of the name
  if (tc != NULL)
    tc[0] = '.';  // Do we really need to overwrite the '/0' before freeing?
  myfree(tname);

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 PARAM_TEST                                */
/*                                                                           */
/*  *** SHORT NAME for "PARAMFILE_TEST_PARAM"                                */
/*                                                                           */
/*  Return 1 if the parameter is in the 'ppl', 0 otherwise.                  */
/*                                                                           */
/*****************************************************************************/
int param_test(ppl,name)
     struct param_pair_list *ppl;
     char *name;
{
  int k;
  char *tname,*tc;
  struct param_pair *tp;

  if (ppl == NULL){ // WYETH this added in case onode model calls
    printf("*** PARAM_TEST:  Warning, ppl is null.  Test %s\n",
	   name);
    return 0;
  }

  //
  // Check for '.' notation, such as "VARGEN_seed.6"
  //
  tname = strdup(name);    // Make a copy so we can modify
  tc = strchr(tname,'.');  // Search for '.' in the name
  if (tc != NULL){
    tc[0] = '\0';          // Change the '.' into a string terminator
  }

  k = 0;
  if (ppl->p != NULL){
    tp = ppl->p;
    while((strcmp(tp->name,tname)!=0) && (tp->next != NULL))
      tp = tp->next;
    if (strcmp(tp->name,tname)==0)
      k = 1;
  }

  // Free the copy of the name
  if (tc != NULL)
    tc[0] = '.';  // Do we really need to overwrite the '/0' before freeing?
  myfree(tname);

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                    PARAMFILE_GET_NEXT_PREFIX_POINTER                      */
/*                                                                           */
/*  Return the pointer to the next node which has a parameter name that      */
/*  begins with the specified prefix.  Return NULL if not found.             */
/*                                                                           */
/*****************************************************************************/
struct param_pair *paramfile_get_next_prefix_pointer(cp,prefix)
     struct param_pair *cp;        // Current pointer
     char *prefix;
{
  int done;
  struct param_pair *tp;

  done = 0;
  tp = cp->next;
  while(!done){
    if (tp == NULL){
      done = 1;
    }else{
      if (compare_prefix_string_order(prefix,tp->name))
	done = 1;
      else
	tp = tp->next;
    }
  }
  return tp;
}
/**************************************-**************************************/
/*                                                                           */
/*                    PARAMFILE_GET_FIRST_PREFIX_POINTER                     */
/*                                                                           */
/*  Return the pointer to the first node which has a parameter name that     */
/*  begins with the specified prefix.  Return NULL if not found.             */
/*                                                                           */
/*****************************************************************************/
struct param_pair *paramfile_get_first_prefix_pointer(ppl,prefix)
     struct param_pair_list *ppl;
     char *prefix;
{
  struct param_pair *tp;

  tp = ppl->p;
  if (tp != NULL){
    if (compare_prefix_string_order(prefix,tp->name) != 1)
      tp = paramfile_get_next_prefix_pointer(tp,prefix);
  }
  return tp;
}
/**************************************-**************************************/
/*                                                                           */
/*                        PARAMFILE_COUNT_PREFIX_PARAM                       */
/*                                                                           */
/*  Count the number of parameters that begin with the prefix.               */
/*                                                                           */
/*****************************************************************************/
int paramfile_count_prefix_param(ppl,prefix)
     struct param_pair_list *ppl;
     char *prefix;
{
  int n;
  struct param_pair *pp;

  n = 0;
  pp = paramfile_get_first_prefix_pointer(ppl,prefix);
  while(pp != NULL){
    n += 1;
    pp = paramfile_get_next_prefix_pointer(pp,prefix);
  }
  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                          PARAMFILE_GET_PARAM_POINTER                      */
/*                                                                           */
/*  Return the pointer to the node with the specified parameter name.        */
/*                                                                           */
/*****************************************************************************/
struct param_pair *paramfile_get_param_pointer(ppl,pname)
     struct param_pair_list *ppl;
     char *pname;
{
  struct param_pair *tp;

  tp = NULL;
  if (ppl->p != NULL){
    tp = ppl->p;
    while((strcmp(tp->name,pname)!=0) && (tp->next != NULL))
      tp = tp->next;
    if (strcmp(tp->name,pname)!=0)
      tp = NULL;
  }
  return tp;
}
/**************************************-**************************************/
/*                                                                           */
/*                      PARAMFILE_GET_PARAM_POINTER_OR_EXIT                  */
/*                                                                           */
/*  Return the pointer to the node with the specified parameter name.        */
/*                                                                           */
/*****************************************************************************/
struct param_pair *paramfile_get_param_pointer_or_exit(ppl,pname)
     struct param_pair_list *ppl;
     char *pname;
{
  struct param_pair *tp;

  if (ppl->p == NULL){
    printf("  pname:  %s\n",pname);
    exit_error("PARAMFILE_GET_PARAM_POINTER_OR_EXIT","Empty parameter list");
  }

  tp = ppl->p;
  while((strcmp(tp->name,pname)!=0) && (tp->next != NULL))
    tp = tp->next;
  if (strcmp(tp->name,pname)!=0){
    printf("name = %s\n",pname);
    exit_error("PARAMFILE_GET_PARAM_POINTER_OR_EXIT","Param name not found");
  }

  return tp;
}
/**************************************-**************************************/
/*                                                                           */
/*                      PARAMFILE_GET_PARAM_POINTER_PREV                     */
/*                                                                           */
/*  Return the pointer to the node *before* the specified parameter name.    */
/*                                                                           */
/*****************************************************************************/
struct param_pair *paramfile_get_param_pointer_prev(ppl,pname,rndx)
     struct param_pair_list *ppl;
     char *pname;
     int *rndx;  // -1 - the desired node was not found
                 //  0 - the index of the desired node
{
  int k;
  struct param_pair *tp,*tprev;

  k = -1;

  tprev = NULL;

  tp = NULL;
  if (ppl->p != NULL){
    tp = ppl->p;
    k = 0;
    while((strcmp(tp->name,pname)!=0) && (tp->next != NULL)){
      tprev = tp;
      tp = tp->next;
      k += 1;
    }
    if (strcmp(tp->name,pname)!=0){
      tp = NULL;
      k = -1;       // Not found, we ran off the end of the list
    }
  }

  *rndx = k;

  return tprev;
}
/**************************************-**************************************/
/*                                                                           */
/*                            PARAMFILE_LIST_DELETE                          */
/*                                                                           */
/*  Remove the parameter with 'name' from the list.                          */
/*  Returned value:                                                          */
/*    -2  - list was empty                                                   */
/*    -1  - 'name' not found in list                                         */
/*     0,1,2...  - index of item before it was removed                       */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/
int paramfile_list_delete(ppl,name)
     struct param_pair_list *ppl;
     char *name;
{
  int k;
  struct param_pair *tp,*tt;

  if (ppl->n < 1){
    return -2;  // List is empty
  }

  tp = paramfile_get_param_pointer_prev(ppl,name,&k);
  if (k == 0){
    tt = ppl->p;        // Item to be deleted
    ppl->p = tt->next;  // Link to next
  }else{
    tt = tp->next;      // Item to be deleted
    tp->next = tt->next;      // Link to next
  }
  ppl->n -= 1;               // Reduce the count

  //printf("DELETING:  %s  %s   TARGNAME: %s   tpName %s\n",
  //tt->name,tt->value,name,tp->name);

  paramfile_free_pair(tt);   // Free the data for deleted element

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                            PARAMFILE_LIST_UPDATE                          */
/*                                                                           */
/*  Update the parameter value *if* it exists in the list.  Return 1 if so,  */
/*  0 otherwise.                                                             */
/*                                                                           */
/*****************************************************************************/
int paramfile_list_update(ppl,name,val)
     struct param_pair_list *ppl;
     char *name,*val;
{
  int k;
  struct param_pair *tp;

  tp = paramfile_get_param_pointer(ppl,name);
  if (tp == NULL)
    k = 0;
  else{
    myfree(tp->value);
    tp->value = strdup(val);
    k = 1;
  }
  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                         PARAMFILE_LIST_UPDATE_INDEX                       */
/*                                                                           */
/*  Update the 'k'th element of the value list for the parameter *if* it     */
/*  exists in the list.  Return 1 if so, 0 otherwise.                        */
/*                                                                           */
/*****************************************************************************/
int paramfile_list_update_index(ppl,name,val,k)
     struct param_pair_list *ppl;
     char *name,*val;
     int k;              // Index for multiple valued parameter, 1..n
{
  int flag,ns;
  char **slist;
  struct param_pair *tp;

  tp = paramfile_get_param_pointer(ppl,name);
  if (tp == NULL){
    flag = 0;
    //printf("HERE WYETH - not found  - NAME = %s\n",name);
  }else{
    get_items_from_string(tp->value,&slist,&ns);

    // Check that there are multiple values, and the right number
    if (ns <= 1){
      printf("  *** Item name:  %s\n",name);
      exit_error("PARAMFILE_LIST_UPDATE_INDEX","Expecting multiple values");
    }else if (ns < k){
      printf("  *** %s has %d values, but %d was requested\n",name,ns,k);
      exit_error("PARAMFILE_LIST_UPDATE_INDEX","Too few values");
    }

    myfree(slist[k-1]);        // Free storage for value to be replaced
    slist[k-1] = strdup(val);  // Replace with new value string
    tp->value = make_string_from_items(slist,ns);  // Rebuild value string

    free_2d_carray(slist,ns);
    flag = 1;
  }
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                         PARAMFILE_UPDATE_PAIR_IN_LIST                     */
/*                                                                           */
/*  Update the parameter value *if* it exists in the list, otherwise, do     */
/*  nothing.                                                                 */
/*                                                                           */
/*****************************************************************************/
void paramfile_update_pair_in_list(ppl,pp)
     struct param_pair_list *ppl;
     struct param_pair *pp;
{
  struct param_pair *tp;

  if (ppl->p != NULL){
    tp = ppl->p;
    while((strcmp(tp->name,pp->name)!=0) && (tp->next != NULL))
      tp = tp->next;
    if (strcmp(tp->name,pp->name)==0){  // Update the value
      myfree(tp->value);
      tp->value = strdup(pp->value);
      // paramfile_free_pair(pp);  WYETH - removed 2008 July 27
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                       PARAMFILE_UPDATE_ADD_PAIR_TO_LIST                   */
/*                                                                           */
/*  Add the parameter pair to the list if a parameter of the same name       */
/*  does not already exist.                                                  */
/*                                                                           */
/*****************************************************************************/
void paramfile_update_add_pair_to_list(ppl,pp)
     struct param_pair_list *ppl;
     struct param_pair *pp;
{
  struct param_pair *tp;

  if (ppl->p == NULL){
    ppl->p = pp;
    ppl->n = 1;
  }else{
    tp = ppl->p;
    while((strcmp(tp->name,pp->name)!=0) && (tp->next != NULL)){
      tp = tp->next;
    }
    if (strcmp(tp->name,pp->name)!=0){
      tp->next = pp;
      ppl->n += 1;
    }else{ /* Update the value */
      myfree(tp->value);
      tp->value = strdup(pp->value);
      paramfile_free_pair(pp);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         PARAMFILE_ADD_PAIR_TO_LIST                        */
/*                                                                           */
/*  Add the parameter pair to the list if a parameter of the same name       */
/*  does not already exist.                                                  */
/*                                                                           */
/*****************************************************************************/
void paramfile_add_pair_to_list(ppl,pp)
     struct param_pair_list *ppl;
     struct param_pair *pp;
{
  struct param_pair *tp;

  if (ppl->p == NULL){
    ppl->p = pp;
    ppl->n = 1;
  }else{
    tp = ppl->p;
    while((strcmp(tp->name,pp->name)!=0) && (tp->next != NULL))
      tp = tp->next;
    if (strcmp(tp->name,pp->name)!=0){
      tp->next = pp;
      ppl->n += 1;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                        PARAMFILE_ADD_PAIR_TO_LIST_C                       */
/*                                                                           */
/*  Add the parameter pair to the list 'c'.                                  */
/*                                                                           */
/*****************************************************************************/
void paramfile_add_pair_to_list_c(ppl,pp)
     struct param_pair_list *ppl;
     struct param_pair *pp;
{
  struct param_pair *tp;

  if (ppl->c == NULL){
    ppl->c = pp;
    ppl->nc = 1;
  }else{
    tp = ppl->c;
    while(tp->next != NULL)
      tp = tp->next;
    
    tp->next = pp;
    ppl->nc += 1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         PARAMFILE_PPL_ADD_NAME_VAL                        */
/*                                                                           */
/*****************************************************************************/
void paramfile_ppl_add_name_val(p,name,val,flag)
     struct param_pair_list *p;
     char *name,*val;
     int flag;  // 0 - add or update if it exists
		// 1 - add or exit if it exists
{
  int i;
  struct param_pair *pp;

  pp = (struct param_pair *)myalloc(sizeof(struct param_pair));
  pp->name = strdup(name);
  pp->value = strdup(val);
  pp->next = NULL;

  if (flag == 1){ /* Exit if it exists */
    if (paramfile_test_param(p,name) == 1)
      exit_error("PARAMFILE_PPL_ADD_NAME_VAL","Name already in list");
  }

  // For flag 0 or 1, it is now safe to add
  paramfile_update_add_pair_to_list(p,pp);
}
/**************************************-**************************************/
/*                                                                           */
/*                              PARAMFILE_INCLUDE                            */
/*                                                                           */
/*****************************************************************************/
void paramfile_include(ppl,infile,modname)
     struct param_pair_list *ppl; // add to this param list
     char infile[];               // file to read from
     char modname[];              // module to include
{
  int i;
  int n,ns,nold,module_flag,nestdepth;
  char **data,**slist,*tag,*incfile,*imodname,*tname;
  struct param_pair *pp;
  
  printf("  PARAMFILE_INCLUDE\n");
  printf("    Including %s from %s\n",modname,infile);

  read_2d_carray(infile,1,0,&data,&n);
  if (n == -1){
    printf("%s\n",infile);
    exit_error("PARAMFILE_INCLUDE","Cannot open parameter file");
  }

  nestdepth = 0;   // depth of nesting inside target module
  module_flag = 0; // becomes 1 when inside module, 2 when done
  for(i=0;i<n;i++){
    get_items_from_string(data[i],&slist,&ns);
    if (ns > 0){
      // Skip comments (#), HTML commands (<), or other tags
      if (slist[0][0] == '<'){
	tag = myxml_get_tag_from_string(data[i]);
	if (strcmp(tag,"include")==0){
	  incfile = myxml_get_value_from_string(data[i],"file");
	  imodname = myxml_get_value_from_string(data[i],"name");
	  printf("FOUND INCLUDE,  file==>%s<==\n",incfile);
	  printf("FOUND INCLUDE,  name==>%s<==\n",imodname);
	  paramfile_include(ppl,incfile,imodname);
	  myfree(incfile);
	  myfree(imodname);
	}else if (strcmp(tag,"module")==0){
	  tname = myxml_get_value_from_string(data[i],"name");
	  if (strcmp(tname,modname)==0)
	    module_flag = 1; // The module to include begins
	  else
	    if (module_flag == 1)
	      nestdepth += 1; // Keep track of nested modules
	  myfree(tname);
	}else if (strcmp(tag,"/module")==0){
	  if (nestdepth == 0)
	    module_flag = 2;     // Exited the module
	  else if (module_flag == 1)
	    nestdepth -= 1;      // Just exited a module
	}
	myfree(tag);
      }else if ((slist[0][0] != '#') && (module_flag == 1)){
	pp = (struct param_pair *)myalloc(sizeof(struct param_pair));
	pp->name = strdup(slist[0]);
	if (ns == 1) // No value given with this name
	  pp->value = NULL;
	else if (ns == 2)
	  pp->value = strdup(slist[1]);
	else
	  pp->value = make_string_from_items(&(slist[1]),ns-1);
	pp->next = NULL;
	if (strcmp(pp->name,"CUSTOMIZE")==0) // Should generalize keyword
	  paramfile_add_pair_to_list_c(ppl,pp);
	else{
	  nold = ppl->n;
	  paramfile_add_pair_to_list(ppl,pp);
	  if (ppl->n == nold){
	    printf("*** WARNING PARAMFILE_INCLUDE %s appears ",pp->name);
	    printf("more than once in %s.\n    Later occurrences ignored.\n",
		   infile);
	  }
	}
      }

      free_2d_carray(slist,ns);
    }
  }
  free_2d_carray(data,n);
}
/**************************************-**************************************/
/*                                                                           */
/*                              PARAMFILE_READ                               */
/*                                                                           */
/*****************************************************************************/
struct param_pair_list *paramfile_read(infile,pflag,xflag)
     char infile[];
     int pflag;
     int xflag;  // 1-exit if duplicate name
{
  int i,j;
  int n,ns,nold,ninl,i0;
  char **data,**slist,*tag,*incfile,*modname,tstr[SLEN];
  struct param_pair_list *ppl;
  struct param_pair *pp;

  if (pflag){
    printf("  PARAMFILE_READ\n");
    printf("    Reading parameters from %s\n",infile);
  }

  read_2d_carray(infile,1,0,&data,&n);
  if (n == -1){
    printf("%s\n",infile);
    exit_error("PARAMFILE_READ","Cannot open parameter file");
  }

  ppl = paramfile_ppl_get_init();

  for(i=0;i<n;i++){  // For each of 'n' lines in 'data'

    get_items_from_string(data[i],&slist,&ns);
    //printf("**** ns=%d   %s\n",ns,data[i]);
    if (ns > 0){
      // Skip comments (#), HTML commands (<), or other tags
      if (slist[0][0] == '<'){
	tag = myxml_get_tag_from_string(data[i]);
	if (strcmp(tag,"include")==0){
	  incfile = myxml_get_value_from_string(data[i],"file");
	  modname = myxml_get_value_from_string(data[i],"name");
	  paramfile_include(ppl,incfile,modname);
	  myfree(incfile);
	  myfree(modname);
	}
	myfree(tag);
      }else if ((int)(slist[0][0]) == 13){  // Vertical tab
	; // Ignore this, which occurs in some files downloaded from web
      }else if (strcmp(slist[0],"INLINE")==0){
	//
	//  Everything following this line is data.
	//  (1) Create one param pair to hold this INLINE instruction
	//  (2) Store the data strings in 'ppl->s'
	//

	pp = (struct param_pair *)myalloc(sizeof(struct param_pair));
	pp->name = strdup("INLINE");
	pp->value = strdup(slist[1]);  // This is INLINE <type>, e.g. 'VARFILE'
	pp->next = NULL;
	paramfile_add_pair_to_list(ppl,pp);

	i0   = i+1;   // Inline data starts at 'i0'
	ninl = n-i0;  // The number of lines of INLINE data
	if (pflag)
	  printf("    INLINE:  %d lines\n",ninl);

	ppl->ns = ninl;
	ppl->s = (char **)myalloc(ninl*sizeof(char *));

	for(j=0;j<ninl;j++)
	  ppl->s[j] = strdup(data[i0+j]);

	i = n;  // Force loop to end

      }else if (slist[0][0] != '#'){
	pp = (struct param_pair *)myalloc(sizeof(struct param_pair));
	pp->name = strdup(slist[0]);

	if (ns == 1) // No value given with this name
	  pp->value = NULL;
	else if (ns == 2)
	  pp->value = strdup(slist[1]);
	else
	  pp->value = make_string_from_items(&(slist[1]),ns-1);
	pp->next = NULL;
	if (strcmp(pp->name,"CUSTOMIZE")==0) // Should generalize keyword
	  paramfile_add_pair_to_list_c(ppl,pp);
	else{
	  nold = ppl->n;
	  paramfile_add_pair_to_list(ppl,pp);
	  if (ppl->n == nold){
	    //printf("==>%s<==\n",pp->name);
	    //printf("%d\n",(int)(pp->name[0]));
	    if (xflag == 1){
	      printf("  *** param name:  %s\n",pp->name);
	      exit_error("PARAMFILE_READ","Duplicate parameter name");
	    }else{
	      printf("*** WARNING PARAMFILE_READ %s appears ",pp->name);
	      printf("more than once in %s.\n    Later occurrences ignored.\n",
		     infile);
	    }
	  }
	}
      }
      free_2d_carray(slist,ns);
    }
  }
  free_2d_carray(data,n);

  return ppl;
}
/**************************************-**************************************/
/*                                                                           */
/*                          PARAMFILE_UPDATE_FROM_LIST                       */
/*                                                                           */
/*****************************************************************************/
void paramfile_update_from_list(p,slist,ns)
     struct param_pair_list *p;
     char **slist;
     int ns;
{
  int i;
  struct param_pair *pp;
  
  /*printf("  PARAMFILE_UPDATE_FROM_LIST\n");*/

  /*** Add or update parameter definitions from command line ***/
  for(i=0;i<ns;i+=2){
    if ((i+1) < ns){
      pp = (struct param_pair *)myalloc(sizeof(struct param_pair));
      pp->name = strdup(slist[i]);
      pp->value = strdup(slist[i+1]);
      pp->next = NULL;
      paramfile_update_add_pair_to_list(p,pp);
    }else{
      printf("*** WARNING last parameter name ignored---no value.\n");
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           PARAMFILE_INIT_FROM_LIST                        */
/*                                                                           */
/*****************************************************************************/
struct param_pair_list *paramfile_init_from_list(slist,ns)
     char **slist;
     int ns;
{
  struct param_pair_list *p;

  p = paramfile_ppl_get_init();
  
  paramfile_update_from_list(p,slist,ns);
  
  return p;
}
/**************************************-**************************************/
/*                                                                           */
/*                         PARAMFILE_GET_INITIAL_PARAMS                      */
/*                                                                           */
/*****************************************************************************/
struct param_pair_list *paramfile_get_initial_params(paramfile,slist,ns,xflag)
     char paramfile[];
     char **slist;
     int ns;
     int xflag;  // 1-exit if duplicate param names found, 0-ignore
{
  int i;
  struct param_pair_list *p;
  struct param_pair *pp;

  printf("  PARAMFILE_GET_INITIAL_PARAMS\n");

  // Get parameter definitions from file
  p = paramfile_read(paramfile,1,xflag);

  /*** WYETH BELOW SHOULD USE 'update_from_list' above !!! ****/

  /*** Add or update parameter definitions from command line ***/
  for(i=0;i<ns;i+=2){
    if ((i+1) < ns){
      pp = (struct param_pair *)myalloc(sizeof(struct param_pair));
      pp->name = strdup(slist[i]);
      pp->value = strdup(slist[i+1]);
      pp->next = NULL; // Wyeth - added
      paramfile_update_add_pair_to_list(p,pp);
    }else{
      printf("*** WARNING last parameter name ignored---no value.\n");
    }
  }

  return p;
}
/**************************************-**************************************/
/*                                                                           */
/*                  PARAMFILE_GET_INITIAL_PARAMS_UPDATE_ONLY                 */
/*                                                                           */
/*****************************************************************************/
struct param_pair_list *paramfile_get_initial_params_update_only(paramfile,
								 slist,ns,
								 xflag)
     char paramfile[];
     char **slist;
     int ns;
     int xflag;
{
  int i,k;
  int flag,uflag;
  char *tc,*tname;
  struct param_pair_list *p;
  struct param_pair *pp;

  printf("  PARAMFILE_GET_INITIAL_PARAMS_UPDATE_ONLY\n");

  // Get parameter definitions from file
  p = paramfile_read(paramfile,1,xflag);

  // Add or update parameter definitions from command line
  for(i=0;i<ns;i+=2){
    if ((i+1) < ns){
      tname = strdup(slist[i]);
      tc = strchr(tname,'.');  // Search for '.' in the name
      if (tc == NULL){ // Simple name
	pp = (struct param_pair *)myalloc(sizeof(struct param_pair));
	pp->name = tname;
	pp->value = strdup(slist[i+1]);
	pp->next = NULL;
	paramfile_update_pair_in_list(p,pp);
	paramfile_free_pair(pp); //  tname freed here - added 2008 July 27
      }else{
	//
	// There is a '.' in the name, so this could be the format like:
	// "VARGEN_seed.6  1777".  Note, the .6 is integer in [1..n]
	//
	flag = is_int_string(&(tc[1]));
	if (flag == 1){
	  k = atoi(&(tc[1])); // Item index
	  //printf("====>  Looking for item number %d\n",k);

	  tc[0] = '\0';  // Change the '.' into a string terminator

	  uflag = paramfile_list_update_index(p,tname,slist[i+1],k);
	  // returns 1 if updated, 0 otherwise
	  if (uflag == 1)
	    ; // printf("  Successful update.\n");
	}
      }

    }else{
      printf("*** WARNING last parameter name ignored---no value.\n");
    }
  }

  return p;
}
/**************************************-**************************************/
/*                                                                           */
/*                           PARAMFILE_PREFIX_DELETE                         */
/*                                                                           */
/*  Remove all elements with a name beginning with 'pname'.                  */
/*  Return the number of elements removed.                                   */
/*                                                                           */
/*****************************************************************************/
int paramfile_prefix_delete(ppl,pname)
     struct param_pair_list *ppl;
     char *pname;  // prefix to remove
{
  int k;
  int ndx;
  struct param_pair *tp,*tt;

  k = 0;
  tp = paramfile_get_first_prefix_pointer(ppl,pname);
  while(tp != NULL){
    //printf("k__ = %d   tp->name = %s    val= %s\n",k,tp->name,tp->value);
    ndx = paramfile_list_delete(ppl,tp->name);
    //printf(" k = %d  ndx = %d\n",k,ndx);
    if (ndx < 0)
      exit_error("PARAMFILE_PREFIX_DELETE","This should not happen");
    
    tp = paramfile_get_first_prefix_pointer(ppl,pname);
    k += 1;
  }

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                   PARAMFILE_GET_NTH_INT_PARAM_PTR_OR_EXIT                 */
/*                                                                           */
/*  Return the 'nth' value for this pointer.                                 */
/*                                                                           */
/*****************************************************************************/
int paramfile_get_nth_int_param_ptr_or_exit(tp,n)
     struct param_pair *tp;
     int n;
{
  int k;
  int ns;
  char **slist;

  if (tp == NULL)
    exit_error("PARAMFILE_GET_INT_PARAM_PTR_OR_EXIT","Null pointer");

  get_items_from_string(tp->value,&slist,&ns);

  k = -1; // Avoid -Wall warning
  if (ns <= n)
    exit_error("PARAMFILE_GET_NTH_INT_PARAM_PTR_OR_EXIT","Not enough values");
  else if (slist[n][0] == '#')
    exit_error("PARAMFILE_GET_NTH_INT_PARAM_PTR_OR_EXIT","nth is comment");
  else
    k = atoi(slist[n]);
  free_2d_carray(slist,ns);

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                      PARAMFILE_GET_NTH_INT_PARAM_OR_EXIT                  */
/*                                                                           */
/*  Return the 'nth' value for the param 'name'.                             */
/*                                                                           */
/*****************************************************************************/
int paramfile_get_nth_int_param_or_exit(ppl,name,n)
     struct param_pair_list *ppl;
     char *name;
     int n;
{
  int k;
  struct param_pair *tp;

  tp = paramfile_get_param_pointer(ppl,name);
  k = paramfile_get_nth_int_param_ptr_or_exit(tp,n);

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                       PARAMFILE_GET_INT_PARAM_OR_EXIT                     */
/*                                                                           */
/*****************************************************************************/
int paramfile_get_int_param_or_exit(ppl,name)
     struct param_pair_list *ppl;
     char *name;
{
  int k;
  struct param_pair *tp;
  int ns;
  char **slist;

  tp = paramfile_get_param_pointer_or_exit(ppl,name);
  get_items_from_string(tp->value,&slist,&ns);
  if (ns > 1){
    if (slist[1][0] != '#'){ /* Assume this is not a comment. */
      printf("*** WARNING PARAMFILE_GET_INT_PARAM_OR_EXIT ");
      printf("parameter %s has multiple values, using first.\n",name);
    }
  }
  k = atoi(slist[0]);
  free_2d_carray(slist,ns);

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                   PARAMFILE_GET_NTH_FLOAT_PARAM_PTR_OR_EXIT               */
/*                                                                           */
/*  Return the 'nth' value for this pointer.                                 */
/*                                                                           */
/*****************************************************************************/
float paramfile_get_nth_float_param_ptr_or_exit(tp,n)
     struct param_pair *tp;
     int n;
{
  int ns;
  char **slist;
  float x;

  x = 0.0; // Avoid -Wall warning

  if (tp == NULL)
    exit_error("PARAMFILE_GET_NTH_FLOAT_PARAM_PTR_OR_EXIT","Null pointer");

  get_items_from_string(tp->value,&slist,&ns);
  if (ns <= n)
    exit_error("PARAMFILE_GET_NTH_FLOAT_PARAM_PTR_OR_EXIT","Too few values");
  else if (slist[n][0] == '#')
    exit_error("PARAMFILE_GET_NTH_FLOAT_PARAM_PTR_OR_EXIT","nth is comment");
  else
    x = atof(slist[n]);
  free_2d_carray(slist,ns);

  return x;
}
/**************************************-**************************************/
/*                                                                           */
/*                    PARAMFILE_GET_NTH_FLOAT_PARAM_OR_EXIT                  */
/*                                                                           */
/*  Item numbers start at 0.                                                 */
/*                                                                           */
/*****************************************************************************/
float paramfile_get_nth_float_param_or_exit(ppl,name,n)
     struct param_pair_list *ppl;
     char *name;
     int n;
{
  float x;
  struct param_pair *tp;

  tp = paramfile_get_param_pointer(ppl,name);
  x = paramfile_get_nth_float_param_ptr_or_exit(tp,n);

  return x;
}
/**************************************-**************************************/
/*                                                                           */
/*                      PARAMFILE_GET_FLOAT_PARAM_OR_EXIT                    */
/*                                                                           */
/*****************************************************************************/
float paramfile_get_float_param_or_exit(ppl,name)
     struct param_pair_list *ppl;
     char *name;
{
  struct param_pair *tp;
  float x;
  int ns;
  char **slist;

  tp = paramfile_get_param_pointer_or_exit(ppl,name);
  get_items_from_string(tp->value,&slist,&ns);
  if (ns > 1){
    if (slist[1][0] != '#'){ // Assume this is not a comment
      printf("*** WARNING PARAMFILE_GET_FLOAT_PARAM_OR_EXIT ");
      printf("parameter %s has multiple values, using first.\n",name);
    }
  }
  x = atof(slist[0]);
  free_2d_carray(slist,ns);

  return x;
}
/**************************************-**************************************/
/*                                                                           */
/*                  PARAMFILE_GET_NTH_CHAR_PARAM_PP_OR_EXIT                  */
/*                                                                           */
/*  Return the 'nth' value for this pointer.                                 */
/*  Returned storage should be freed.                                        */
/*                                                                           */
/*****************************************************************************/
char *paramfile_get_nth_char_param_pp_or_exit(tp,n)
     struct param_pair *tp;
     int n;
{
  int ns;
  char **slist,*t;

  t = NULL; // Avoid -Wall warning

  if (tp == NULL)
    exit_error("PARAMFILE_GET_NTH_CHAR_PARAM_PTR_OR_EXIT","Null pointer");

  get_items_from_string(tp->value,&slist,&ns);
  if (ns <= n)
    exit_error("PARAMFILE_GET_NTH_CHAR_PARAM_PTR_OR_EXIT","Not enough values");
  else if (slist[n][0] == '#')
    exit_error("PARAMFILE_GET_NTH_CHAR_PARAM_PTR_OR_EXIT","nth is comment");
  else
    t = strdup(slist[n]);
  free_2d_carray(slist,ns);

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                     PARAMFILE_GET_NTH_CHAR_PARAM_OR_EXIT                  */
/*                                                                           */
/*****************************************************************************/
char *paramfile_get_nth_char_param_or_exit(ppl,name,n)
     struct param_pair_list *ppl;
     char *name;
     int n;
{
  struct param_pair *tp;
  char *t;

  tp = paramfile_get_param_pointer(ppl,name);
  t = paramfile_get_nth_char_param_pp_or_exit(tp,n);

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                      PARAMFILE_GET_CHAR_PARAM_OR_EXIT                     */
/*                                                                           */
/*****************************************************************************/
char *paramfile_get_char_param_or_exit(ppl,name)
     struct param_pair_list *ppl;
     char *name;
{
  struct param_pair *tp;
  int ns;
  char **slist,*t;

  tp = paramfile_get_param_pointer_or_exit(ppl,name);
  get_items_from_string(tp->value,&slist,&ns);
  if (ns > 1){
    if (slist[1][0] != '#'){ // Assume this is not a comment.
      printf("*** WARNING PARAMFILE_GET_CHAR_PARAM_OR_EXIT ");
      printf("parameter %s has multiple values, using first.\n",name);
    }
  }
  t = strdup(slist[0]);
  free_2d_carray(slist,ns);

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                       PARAMFILE_GET_INT_PARAM_DEFAULT                     */
/*                                                                           */
/*****************************************************************************/
int paramfile_get_int_param_default(ppl,name,defval)
     struct param_pair_list *ppl;
     char *name;
     int defval;
{
  int k;
  struct param_pair *tp;
  int ns;
  char **slist;

  k = 0; // Avoid -Wall warning

  if (ppl->p == NULL){
    //exit_error("PARAMFILE_GET_INT_PARAM_DEFAULT","Empty parameter list");
    k = defval;
  }else{
    tp = ppl->p;
    while((strcmp(tp->name,name)!=0) && (tp->next != NULL))
      tp = tp->next;
    if (strcmp(tp->name,name)!=0)
      k = defval;
    else{
      get_items_from_string(tp->value,&slist,&ns);
      if (ns > 1){
	if (slist[1][0] != '#'){ /* Assume this is not a comment. */
	  printf("*** WARNING PARAMFILE_GET_INT_PARAM_DEFAULT ");
	  printf("parameter %s has multiple values, using first.\n",name);
	}
      }
      k = atoi(slist[0]);
      free_2d_carray(slist,ns);
    }
  }
  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                      PARAMFILE_GET_FLOAT_PARAM_DEFAULT                    */
/*                                                                           */
/*****************************************************************************/
float paramfile_get_float_param_default(ppl,name,defval)
     struct param_pair_list *ppl;
     char *name;
     float defval;
{
  struct param_pair *tp;
  float x;
  int ns;
  char **slist;

  x = 0.0; /* Avoid -Wall warning */

  if (ppl->p == NULL){
    //exit_error("PARAMFILE_GET_FLOAT_PARAM_DEFAULT","Empty parameter list");
    x = defval;
  }else{
    tp = ppl->p;
    while((strcmp(tp->name,name)!=0) && (tp->next != NULL))
      tp = tp->next;
    if (strcmp(tp->name,name)!=0)
      x = defval;
    else{
      get_items_from_string(tp->value,&slist,&ns);
      if (ns > 1){
	if (slist[1][0] != '#'){ /* Assume this is not a comment. */
	  printf("*** WARNING PARAMFILE_GET_FLOAT_PARAM_DEFAULT ");
	  printf("parameter %s has multiple values, using first.\n",name);
	}
      }
      x = atof(slist[0]);
      free_2d_carray(slist,ns);
    }
  }
  return x;
}
/**************************************-**************************************/
/*                                                                           */
/*                      PARAMFILE_GET_CHAR_PARAM_DEFAULT                     */
/*                                                                           */
/*****************************************************************************/
char *paramfile_get_char_param_default(ppl,name,defval)
     struct param_pair_list *ppl;
     char *name;
     char *defval;
{
  struct param_pair *tp;
  int ns;
  char **slist,*t;

  t = NULL; // Avoid -Wall warning

  if (ppl->p == NULL){
    //exit_error("PARAMFILE_GET_CHAR_PARAM_DEFAULT","Empty parameter list");
    t = NULL;
  }else{
    tp = ppl->p;
    while((strcmp(tp->name,name)!=0) && (tp->next != NULL))
      tp = tp->next;
    if (strcmp(tp->name,name)!=0){
      if (defval == NULL)
	t = NULL;
      else
	t = strdup(defval);
    }else{
      get_items_from_string(tp->value,&slist,&ns);
      if (ns > 1){
	if (slist[1][0] != '#'){ // Assume this is not a comment
	  printf("*** WARNING PARAMFILE_GET_CHAR_PARAM_DEFAULT ");
	  printf("parameter %s has multiple values, using first.\n",name);
	}
      }
      t = strdup(slist[0]);
      free_2d_carray(slist,ns);
    }
  }
  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                              PARAM_GETI_DFLT                              */
/*                                                                           */
/*****************************************************************************/
int param_geti_dflt(ppl,name,defval)
     struct param_pair_list *ppl;
     char *name;
     int defval;
{
  if (ppl == NULL)  // WYETH 2014 Jan
    return defval;
  else
    return paramfile_get_int_param_default(ppl,name,defval);
}
/**************************************-**************************************/
/*                                                                           */
/*                              PARAM_GETF_DFLT                              */
/*                                                                           */
/*****************************************************************************/
float param_getf_dflt(ppl,name,defval)
     struct param_pair_list *ppl;
     char *name;
     float defval;
{
  if (ppl == NULL)
    return defval;
  else
    return paramfile_get_float_param_default(ppl,name,defval);
}
/**************************************-**************************************/
/*                                                                           */
/*                              PARAM_GETC_DFLT                              */
/*                                                                           */
/*  This allocates storage - should free returned value.                     */
/*                                                                           */
/*****************************************************************************/
char *param_getc_dflt(ppl,name,defval)
     struct param_pair_list *ppl;
     char *name;
     char *defval;
{
  if (ppl == NULL){
    if (defval == NULL)  //
      return NULL;       //  Wyeth - added this part of condition Sept 2016
    else                 //
      return strdup(defval);
  }else
    return paramfile_get_char_param_default(ppl,name,defval);
}
/**************************************-**************************************/
/*                                                                           */
/*                              PARAM_GETI_EXIT                              */
/*                                                                           */
/*****************************************************************************/
int param_geti_exit(ppl,name)
     struct param_pair_list *ppl;
     char *name;
{
  return paramfile_get_int_param_or_exit(ppl,name);
}
/**************************************-**************************************/
/*                                                                           */
/*                              PARAM_GETF_EXIT                              */
/*                                                                           */
/*****************************************************************************/
float param_getf_exit(ppl,name)
     struct param_pair_list *ppl;
     char *name;
{
  return paramfile_get_float_param_or_exit(ppl,name);
}
/**************************************-**************************************/
/*                                                                           */
/*                              PARAM_GETC_EXIT                              */
/*                                                                           */
/*****************************************************************************/
char *param_getc_exit(ppl,name)
     struct param_pair_list *ppl;
     char *name;
{
  return paramfile_get_char_param_or_exit(ppl,name);
}
/**************************************-**************************************/
/*                                                                           */
/*                             PARAM_GETI_S2_EXIT                            */
/*                                                                           */
/*****************************************************************************/
void param_geti_s2_exit(ppl,name,r1,r2,rflag)
     struct param_pair_list *ppl;
     char *name;
     int *r1,*r2;  // Values for two params, "s1_[name]", "s2_[name]"
     int *rflag;     // set to '2' if dual names, otherwise do not change
{
  char tname[LONG_SLEN];

  if (param_test(ppl,name)){
    *r1 = paramfile_get_int_param_or_exit(ppl,name);
    *r2 = *r1;
  }else{
    sprintf(tname,"s1_%s",name);
    *r1 = paramfile_get_int_param_or_exit(ppl,tname);
  }

  sprintf(tname,"s2_%s",name);
  if (param_test(ppl,tname)){
    *r2 = paramfile_get_int_param_or_exit(ppl,tname);
    *rflag = 2;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             PARAM_GETF_S2_EXIT                            */
/*                                                                           */
/*****************************************************************************/
void param_getf_s2_exit(ppl,name,r1,r2,rflag)
     struct param_pair_list *ppl;
     char  *name;
     float *r1,*r2;  // Values for two params, "s1_[name]", "s2_[name]"
     int   *rflag;   // set to '2' if dual names, otherwise do not change
{
  char tname[LONG_SLEN];

  if (param_test(ppl,name)){
    *r1 = paramfile_get_float_param_or_exit(ppl,name);
    *r2 = *r1;
  }else{
    sprintf(tname,"s1_%s",name);
    *r1 = paramfile_get_float_param_or_exit(ppl,tname);
  }

  sprintf(tname,"s2_%s",name);
  if (param_test(ppl,tname)){
    *r2 = paramfile_get_float_param_or_exit(ppl,tname);
    *rflag = 2;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                PARAM_GETF_2                               */
/*                                                                           */
/*****************************************************************************/
void param_getf_2(ppl,name,name1,name2,r1,r2)
     struct param_pair_list *ppl;
     char *name,*name1,*name2; // e.g., "con" "con1" "con2"
     float *r1,*r2;  // Values for two params
{
  char tname[LONG_SLEN];

  if (param_test(ppl,name)){  // If the first name exists, use it
    *r1 = paramfile_get_float_param_or_exit(ppl,name);
    *r2 = *r1;
  }else{                      // Use separate names
    *r1 = paramfile_get_float_param_or_exit(ppl,name1);
    *r2 = paramfile_get_float_param_or_exit(ppl,name2);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           PARAMFILE_GET_FLT_2OF3                          */
/*                                                                           */
/*****************************************************************************/
int paramfile_get_flt_2of3(cname,ppl,p1,p2,p3,r1,r2,r3)
     char *cname;                     // Calling routine
     struct param_pair_list *ppl;
     char *p1,*p2,*p3;                // three parameter names
     float *r1,*r2,*r3;               // returned values
{
  int k,n;

  k = -1;  // Index of parameter that is missing
  n = 0;   // Number of parameters found

  if (paramfile_test_param(ppl,p1)){
    *r1 = paramfile_get_float_param_or_exit(ppl,p1);
    n += 1;
  }else
    k = 1;

  if (paramfile_test_param(ppl,p2)){
    *r2 = paramfile_get_float_param_or_exit(ppl,p2);
    n += 1;
  }else
    k = 2;

  if (paramfile_test_param(ppl,p3)){
    *r3 = paramfile_get_float_param_or_exit(ppl,p3);
    n += 1;
  }else
    k = 3;

  if (n != 2){
    printf("  *** Called from %s\n",cname);
    printf("  *** Please set exactly two of:  %s %s %s\n",p1,p2,p3);
    exit_error("PARAMFILE_GET_FLT_2OF3","Wrong parameter combination");
  }

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                              PARAM_GETF_LRD                               */
/*                                                                           */
/*  One-step setting of two variables with a base name 'foo', and actual     */
/*  names:  foo_L  foo_R  foo_disp                                           */
/*                                                                           */
/*  If 'foo' is set, then set foo_L and foo_R to be this value.              */
/*  Otherwise, two of the three must be set, and the third will be           */
/*  computed by solving:                                                     */
/*                                                                           */
/*     foo_R = foo_L + foo_disp                                              */
/*                                                                           */
/*  This is used for binocular stimuli, where a parameter may be set by      */
/*  specifying values for the Left Eye, Right Eye or one eye and a           */
/*  disparity value.                                                         */
/*                                                                           */
/*****************************************************************************/
void param_getf_lrd(cname,ppl,pname,rl,rr)
     char *cname;                   // Calling routine
     struct param_pair_list *ppl;
     char *pname;                   // base name of parameter
     float *rl,*rr;                 // returned values, left, right, disp
{
  int res;
  float vl,vr,vd;
  char lstr[SLEN],rstr[SLEN],dstr[SLEN];

  sprintf(lstr,"%s_L",pname);
  sprintf(rstr,"%s_R",pname);
  sprintf(dstr,"%s_disp",pname);

  if (paramfile_test_param(ppl,pname) == 1){
    vl = vr = param_getf_exit(ppl,pname);

    // Tell user not to set redundant paramter names
    if ((paramfile_test_param(ppl,lstr) == 1)||
	(paramfile_test_param(ppl,rstr) == 1)||
	(paramfile_test_param(ppl,dstr) == 1)){
      printf("  *** Setting the stimulus parameter: %s\n",pname);
      printf("      precludes the use of:  %s %s %s\n",lstr,rstr,dstr);
      exit_error("PARAM_GETF_LRD","Please remove extra stimulus parameters");
    }
  }else{
    res = paramfile_get_flt_2of3(cname,ppl,lstr,rstr,dstr,&vl,&vr,&vd);

    if (res == 1)
      vl = vr - vd;
    else if (res == 2)
      vr = vl + vd;
  }

  *rl = vl;
  *rr = vr;
}
/**************************************-**************************************/
/*                                                                           */
/*                       PARAMFILE_GET_COMMENT_FOR_PARAM                     */
/*                                                                           */
/*  Return the comment string for 'name', or "" if no comment.               */
/*  Exits if 'name' not found.                                               */
/*  First '#' is not included in returned string.                            */
/*                                                                           */
/*****************************************************************************/
char *paramfile_get_comment_for_param(ppl,name)
     struct param_pair_list *ppl;
     char *name;
{
  int i;
  int ns,cflag,fflag;
  char **slist,*t,tstr[2048];
  struct param_pair *tp;

  tp = paramfile_get_param_pointer_or_exit(ppl,name);

  get_items_from_string(tp->value,&slist,&ns);
  cflag = 0;  /* Found "#" */
  fflag = 0;  /* First item has been appended */
  tstr[0] = '\0';
  for(i=0;i<ns;i++){
    if (cflag){
      if (fflag)
	strcat(tstr," ");
      else
	fflag = 1;
      strcat(tstr,slist[i]);
    }else{
      if (slist[i][0] == '#'){
	cflag = 1;
	if (strlen(slist[i]) > 1){
	  fflag = 1;
	  strcat(tstr,slist[i]+1);
	}
      }
    }
  }
  t = strdup(tstr);
  free_2d_carray(slist,ns);

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                        PARAMFILE_GET_TYPE_FROM_COMMENT                    */
/*                                                                           */
/*  Look for the character c in '[c]' following the '#' comment sign.        */
/*                                                                           */
/*****************************************************************************/
char paramfile_get_type_from_comment(ppl,name,dflt)
     struct param_pair_list *ppl;
     char *name;
     char dflt;   /* Default value to use if not found */
{
  int i;
  int n,state;
  char *comstr,c;

  comstr = paramfile_get_comment_for_param(ppl,name);
  n = strlen(comstr);

  c = dflt;  /* Use default, unless [x] found below */

  i = 0;
  state = 0;
  while(state >= 0){
    if (i >= n){  /* We have run off of the string, no type specifier */
      state = -1;
      c = dflt;  /* Revert to default */
    }else if (state == 0){  /* Haven't found the '[' yet */
      if (comstr[i] == '[')
	state = 1;
      else if ((comstr[i] == ' ')||(comstr[i] == '#'))
	; /* Stay in the state, keep looking for '[' */
      else
	state = -1; /* Done, there is no type specifier */
    }else if (state == 1){
      c = comstr[i];
      state = 2;
    }else if (state == 2){
      if (comstr[i] == ']')
	state = -1; /* Done, we found it */
      else{
	c = dflt;  /* Revert to default */
	state = -1; /* Done, we found it */
      }
    }
    /*printf("%d CHAR %c\n",state,comstr[i]);*/
    i += 1;
  }

  if (!((c=='i')||(c=='f')||(c=='c'))){
    printf("In:  %s\n",comstr);
    exit_error("PARAMFILE_GET_TYPE_FROM_COMMENT","Invalid type string");
  }

  myfree(comstr);

  return c;
}
/**************************************-**************************************/
/*                                                                           */
/*                        PARAMFILE_COUNT_VALUES_POINTER                     */
/*                                                                           */
/*  Count the number of values for this pointer.  Return -1 if the pointer   */
/*  is NULL.  Do not count values at or following any value beginning        */
/*  with '#'.                                                                */
/*                                                                           */
/*****************************************************************************/
int paramfile_count_values_pointer(tp)
     struct param_pair *tp;
{
  int n,sn,done;
  char **slist;

  if (tp != NULL){
    get_items_from_string(tp->value,&slist,&sn);

    n = 0;
    done = 0;
    while(!done){
      if (n >= sn)
	done = 1;
      else if (slist[n][0] == '#')
	done = 1;
      else
	n++;
    }
    free_2d_carray(slist,sn);
  }else
    n = -1;
  
  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                         PARAMFILE_COUNT_VALUES_PARAM                      */
/*                                                                           */
/*  Count the number of values for this parameter.  Return -1 if the param   */
/*  is not found.  Do not count values at or following any value beginning   */
/*  with '#'.                                                                */
/*                                                                           */
/*****************************************************************************/
int paramfile_count_values_param(ppl,name)
     struct param_pair_list *ppl;
     char *name;
{
  int n;
  struct param_pair *tp;

  tp = paramfile_get_param_pointer(ppl,name);
  n = paramfile_count_values_pointer(tp);
  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                         PARAMFILE_GET_SLIST_POINTER                       */
/*                                                                           */
/*****************************************************************************/
int paramfile_get_slist_pointer(ppl,tp,rslist)
     struct param_pair_list *ppl;
     struct param_pair *tp;
     char ***rslist;
{
  int i;
  int n,done;
  char **slist;

  *rslist = NULL;
  n = 0;

  if (tp != NULL){
    get_items_from_string(tp->value,&slist,&n);

    i = 0;
    done = 0;
    while(!done){
      if (i >= n)
	done = 1;
      else if (slist[i][0] == '#')
	done = 1;
      else
	i++;
    }
    if (i<n){ // There is a comment on this line
      *rslist = copy_2d_carray(slist,i);
      free_2d_carray(slist,n);
      n = i;
    }else
      *rslist = slist;
  }
  
  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                         PARAMFILE_GET_SLIST_PARAM                         */
/*                                                                           */
/*****************************************************************************/
int paramfile_get_slist_param(ppl,name,rslist)
     struct param_pair_list *ppl;
     char *name,***rslist;
{
  int n;
  struct param_pair *tp;

  *rslist = NULL;
  n = 0;

  tp = paramfile_get_param_pointer(ppl,name);
  n = paramfile_get_slist_pointer(ppl,tp,rslist);

  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                            PARAMFILE_MAKE_FNAME                           */
/*                                                                           */
/*  Build the name of the output file.                                       */
/*                                                                           */
/*  Return (char *)NULL if "name" is not found.                              */
/*  If "name" is found and its value is                                      */
/*    == "*"  -->  return "prefix"+"extension"                               */
/*    != "*"  -->  return "name"                                             */
/*                                                                           */
/*****************************************************************************/
char *paramfile_make_fname(ppl,name,prefix,extension)
     struct param_pair_list *ppl;
     char *name;       // Look up this parameter name
     char *prefix;     // If "name" is "*", then use the prefix and extension
     char *extension;  // file extension
{
  char *tname,tstr[SLEN];

  if (paramfile_test_param(ppl,name)){
    tname = paramfile_get_nth_char_param_or_exit(ppl,name,0);
    if (strcmp(tname,"*")==0)
      sprintf(tstr,"%s%s",prefix,extension);
    else
      sprintf(tstr,"%s",tname);
    myfree(tname);
  }else
    return (char *)NULL;

  return strdup(tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                              PARAMFILE_READ_DIFF                          */
/*                                                                           */
/*  List any differences in the parameter files.                             */
/*                                                                           */
/*****************************************************************************/
void paramfile_read_diff(infile1,infile2)
     char infile1[],infile2[];
{
  int i,j;
  int ndiff,nomit,n1,n2,match,nc1,nc2;
  char *s1,*s2,*tname;
  struct param_pair *p,*q;
  struct param_pair_list *p1,*p2;

  p1 = paramfile_read(infile1,1,0);
  p2 = paramfile_read(infile2,1,0);

  printf("\n  EXPLICIT CONFLICTS\n\n");
  p = p1->p;
  ndiff = 0;
  for(i=0;i<p1->n;i++){
    tname = p->name;
    if (paramfile_test_param(p2,tname) == 1){
      n1 = paramfile_count_values_param(p1,tname);
      n2 = paramfile_count_values_param(p2,tname);
      if ((n1==1)&&(n2==1)){
	s1 = paramfile_get_char_param_or_exit(p1,tname);
	s2 = paramfile_get_char_param_or_exit(p2,tname);
	if (strcmp(s1,s2)!=0){
	  printf("      %30s  %16s %16s\n",tname,s1,s2);
	  ndiff += 1;
	}
	myfree(s1); myfree(s2);
      }else if (n1==n2){
	j = 0;
	while(j<n1){
	  s1 = paramfile_get_nth_char_param_or_exit(p1,tname,j);
	  s2 = paramfile_get_nth_char_param_or_exit(p2,tname,j);
	  if (strcmp(s1,s2)!=0){
	    q = paramfile_get_param_pointer(p2,tname);
	    printf("      %30s  MULTI-VALUED PARAMETER\n",tname);
	    printf("      %30s  %s\n","(file1)",p->value);
	    printf("      %30s  %s\n","(file2)",q->value);
	    ndiff += 1;
	    j = n1; /*** DON'T CHECK THE REST ***/
	  }
	  myfree(s1); myfree(s2);
	  j += 1;
	}
      }else{
	q = paramfile_get_param_pointer(p2,tname);
	printf("      %30s  NUMBER OF VALUES DIFFER\n",tname);
	printf("      %30s  %s\n","(file1)",p->value);
	printf("      %30s  %s\n","(file2)",q->value);
	ndiff += 1;
      }
    }
    p = p->next;
  }
  if (ndiff == 0)
    printf("      None.\n\n");


  printf("\n  OMISSIONS\n\n");
  p = p1->p;
  nomit = 0;
  for(i=0;i<p1->n;i++){
    tname = p->name;
    if (paramfile_test_param(p2,tname) != 1){
      n1 = paramfile_count_values_param(p1,tname);
      if (n1 == 1){
	s1 = paramfile_get_char_param_or_exit(p1,tname);
	printf("      %30s  %16s\n",tname,s1);
	myfree(s1);
      }else
	printf("      %30s  (multiple values)\n",tname);
      nomit += 1;
    }
    p = p->next;
  }

  p = p2->p;
  for(i=0;i<p2->n;i++){
    tname = p->name;
    if (paramfile_test_param(p1,tname) != 1){
      n2 = paramfile_count_values_param(p2,tname);
      if (n2 == 1){
	s2 = paramfile_get_char_param_or_exit(p2,tname);
	printf("      %30s                   %16s\n",tname,s2);
	myfree(s2);
      }else
	printf("      %30s                   (multiple values)\n",tname);
      nomit += 1;
    }
    p = p->next;
  }

  if (ndiff == 0)
    printf("      None.\n\n");


  printf("\n  CUSTOMIZE differences\n\n");

  nc1 = 0;
  p = p1->c;
  for(i=0;i<p1->nc;i++){
    j = 0;
    match = 0;
    q = p2->c;
    for(j=0;j<p2->nc;j++){
      if ((strcmp(p->name,q->name)==0) && (strcmp(p->value,q->value)==0))
	match = 1;
      q = q->next;
    }
    if (match == 0){
      printf("    (file1)  %s %s\n",p->name,p->value);
      nc1 += 1;
    }
    p = p->next;
  }

  nc2 = 0;
  p = p2->c;
  for(i=0;i<p2->nc;i++){
    j = 0;
    match = 0;
    q = p1->c;
    for(j=0;j<p1->nc;j++){
      if ((strcmp(p->name,q->name)==0) && (strcmp(p->value,q->value)==0))
	match = 1;
      q = q->next;
    }
    if (match == 0){
      printf("    (file2)  %s %s\n",p->name,p->value);
      nc2 += 1;
    }
    p = p->next;
  }
  if ((nc1 == 0) && (nc2 == 0))
    printf("      None.\n\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                            PARAMFILE_ONODE_FREE                           */
/*                                                                           */
/*  Free this onode and recursively free the entire tree structure below.    */
/*                                                                           */
/*****************************************************************************/
void paramfile_onode_free(o)
     struct onode *o;
{
  int i;
  struct onode *t,*next;

  if (o == NULL)
    return;

  //
  //  Free all tag values
  //
  if (o->tv_n > 0){
    for(i=0;i<o->tv_n;i++){
      if (o->tv_name[i] != NULL)
	myfree(o->tv_name[i]);
      if (o->tv[i] != NULL)
	myfree(o->tv[i]);
    }
    // WYETH - Added 11 April 2010
    if (o->tv_name != NULL)
      myfree(o->tv_name);
    if (o->tv != NULL)
      myfree(o->tv);
  }

  //
  //  Recursively free each child
  //
  t = o->o;
  while(t != NULL){
    next = t->next;
    paramfile_onode_free(t);
    t = next;
  }

  if (o->otype != NULL)   myfree(o->otype);
  if (o->otag != NULL)    myfree(o->otag);
  if (o->name != NULL)    myfree(o->name);
  if (o->val  != NULL)    myfree(o->val);
  if (o->unit != NULL)    myfree(o->unit);
  if (o->comm != NULL)    myfree(o->comm);


  myfree(o); // WYETH - ADDED 11 Apr 2010

  totocnt -= 1;

  if (totocnt < 0)
    printf("totocnt = %d\n",totocnt);

}
/**************************************-**************************************/
/*                                                                           */
/*                            PARAMFILE_ONODE_WRITE                          */
/*                                                                           */
/*****************************************************************************/
void paramfile_onode_write(o,k,fp)
     struct onode *o;
     int k;             // nesting level
     FILE *fp;
{
  int i;
  struct onode *t;

  for(i=0;i<k;i++)
    fprintf(fp,"  ");
  if (o->tv_n == 0)
    fprintf(fp,"<%s>\n",o->otype);
  else{
    fprintf(fp,"<%s",o->otype);
    for(i=0;i<o->tv_n;i++){
      fprintf(fp," %s=\"%s\"",o->tv_name[i],o->tv[i]);
    }
    fprintf(fp,">\n");
  }
  
  t = o->o;
  while(t != NULL){
    if (strcmp(t->otype,"item")==0){
      for(i=0;i<k+1;i++)
	fprintf(fp,"  ");
      fprintf(fp,"%s %s %s %s\n",t->name,t->val,t->unit,t->comm);
    }else{
      paramfile_onode_write(t,k+1,fp);
    }
    t = t->next;
  }

  for(i=0;i<k;i++)
    fprintf(fp,"  ");
  fprintf(fp,"</%s>\n",o->otype);

}
/**************************************-**************************************/
/*                                                                           */
/*                         PARAMFILE_ONODE_WRITE_FILE                        */
/*                                                                           */
/*****************************************************************************/
void paramfile_onode_write_file(o,outfile)
     struct onode *o;
     char *outfile;
{
  FILE *open,*fout;

  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** File %s\n",outfile);
    exit_error("PARAMFILE_ONODE_WRITE","Cannot open file");
  }

  paramfile_onode_write(o,0,fout);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                            PARAMFILE_ONODE_PRINT                          */
/*                                                                           */
/*****************************************************************************/
void paramfile_onode_print(o,k,iflag)
     struct onode *o;
     int k;             // nesting level
     int iflag;         // 0-summarise item counts,  1-print all items
{
  int i;
  int nitem;
  struct onode *t;

  for(i=0;i<k;i++)
    printf("  ");

  if (o->tv_n == 0)
    printf("<%s>\n",o->otype);
  else{
    printf("<%s",o->otype);
    for(i=0;i<o->tv_n;i++){
      printf(" %s=\"%s\"",o->tv_name[i],o->tv[i]);
    }
    printf(">\n");
  }

  nitem = 0;  // Item count
  t = o->o;
  while(t != NULL){
    if (strcmp(t->otype,"item")==0){
      if (iflag){
	for(i=0;i<k+1;i++)
	  printf("  ");
	printf("%s %s %s %s\n",t->name,t->val,t->unit,t->comm);
      }else{
	nitem += 1;
      }
    }else{
      paramfile_onode_print(t,k+1,iflag);
    }

    t = t->next;
  }
  if (iflag == 0){ // Print item count
    for(i=0;i<k+1;i++)
      printf("  ");
    printf("%d items.\n",nitem);
  }else{
    for(i=0;i<k;i++)
      printf("  ");
    printf("</%s>\n",o->otype);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                        PARAMFILE_ONODE_CREATE_ONODE                       */
/*                                                                           */
/*****************************************************************************/
struct onode *paramfile_onode_create_onode()
{
  struct onode *o;
  
  o = (struct onode *)myalloc(sizeof(struct onode));

  o->otype = NULL;

  o->resolve = 0;     // Nothing to resolve

  o->tv_n    = 0;     // Tag value assignments
  o->tv_name = NULL;
  o->tv      = NULL;

  o->otag = NULL;
  o->n = 0;
  o->o    = NULL;

  o->prev = NULL;
  o->next = NULL;

  o->name = NULL;
  o->val  = NULL;
  o->unit = NULL;
  o->comm = NULL;

  o->nval = 0;

  //
  //  Keep track of total onodes created and freed.
  //
  totocnt += 1;
  if (totocnt > 100000)
    printf("totocnt = %d\n",totocnt);

  return o;
}
/**************************************-**************************************/
/*                                                                           */
/*                           PARAMFILE_ONODE_CARRAY                          */
/*                                                                           */
/*  Return a carray that describes the onode.                                */
/*                                                                           */
/*****************************************************************************/
char *paramfile_onode_carray(o,rn)
     struct onode *o;
     int *rn;             // Return the length of the carray
{
  int i,j,k;
  int n,tn,*cn;
  char tstr_n[SLEN],tstr_ntv[SLEN],tstr_nres[SLEN],*s,*so,**cs,*t;
  struct onode *op;

  // Get carray 's' for this onode

  // First, compute the number of characters, 'n'
  n = 0;
  if (o->otype != NULL)    n += 1 + strlen(o->otype); else    n += 1;
  if (o->otag  != NULL)    n += 1 + strlen(o->otag);  else    n += 1;
  if (o->name  != NULL)    n += 1 + strlen(o->name);  else    n += 1;
  if (o->val   != NULL)    n += 1 + strlen(o->val);   else    n += 1;
  if (o->unit  != NULL)    n += 1 + strlen(o->unit);  else    n += 1;
  if (o->comm  != NULL)    n += 1 + strlen(o->comm);  else    n += 1;

  sprintf(tstr_nres,"%d",o->resolve); // Resolve flag
  n += 1 + strlen(tstr_nres);
  
  sprintf(tstr_n,"%d",o->n); // Number of children
  n += 1 + strlen(tstr_n);

  //
  //  WYETH - this is a new block of code, for TAG VALUES
  //  April 26, 2008
  //
  sprintf(tstr_ntv,"%d",o->tv_n); // Tag value list
  n += 1 + strlen(tstr_ntv);
  for(i=0;i<o->tv_n;i++){
    n += 1 + strlen(o->tv_name[i]);
    n += 1 + strlen(o->tv[i]);
  }
  //  WYETH - END NEW BLOCK


  // Second, create string 's' to hold this onode
  s = (char *)myalloc(n*sizeof(char));  // Create string to hold onode
  t = s;

  if (o->otype != NULL){
    strcpy(t,o->otype);   t += 1 + strlen(o->otype);
  }else{
    *t = '\0'; t += 1;
  }

  strcpy(t,tstr_nres);    t += 1 + strlen(tstr_nres);  // Resolve flag

  //
  //
  //  WYETH - this is a new block of code, for TAG VALUES
  //  April 26, 2008
  //
  strcpy(t,tstr_ntv);         t += 1 + strlen(tstr_ntv);  // No. of tag vals
  for(i=0;i<o->tv_n;i++){
    strcpy(t,o->tv_name[i]);  t += 1 + strlen(o->tv_name[i]);
    strcpy(t,o->tv[i]);       t += 1 + strlen(o->tv[i]);
  }
  //  WYETH - END NEW BLOCK

  if (o->otag != NULL){
    strcpy(t,o->otag);    t += 1 + strlen(o->otag);
  }else{
    *t = '\0'; t += 1;
  }

  strcpy(t,tstr_n);       t += 1 + strlen(tstr_n);   // Number of children

  if (o->name != NULL){
    strcpy(t,o->name);    t += 1 + strlen(o->name);
  }else{
    *t = '\0'; t += 1;
  }
  if (o->val != NULL){
    strcpy(t,o->val);     t += 1 + strlen(o->val);
  }else{
    *t = '\0'; t += 1;
  }
  if (o->unit != NULL){
    strcpy(t,o->unit);    t += 1 + strlen(o->unit);
  }else{
    *t = '\0'; t += 1;
  }
  if (o->comm != NULL){
    strcpy(t,o->comm);    t += 1 + strlen(o->comm);
  }else{
    *t = '\0'; t += 1;
  }

  // NOTE, 'nval' will be derived from the value string, thus is not stored


  // Third, combine 's' with the strings for any children
  if (strcmp(o->otype,"item")==0){  // 's' will be returned
    tn = n;
    so = s;
  }else{ // Concatenate 's' with all children strings
    cs = (char **)myalloc(o->n*sizeof(char *));
    cn = (int *)myalloc(o->n*sizeof(int));

    tn = n;
    i = 0;
    op = o->o;
    while(op != NULL){
      cs[i] = paramfile_onode_carray(op,&(cn[i]));
      tn += cn[i];
      i += 1;
      op = op->next;
    }

    so = (char *)myalloc(tn*sizeof(char));
    k = 0;
    for(j=0;j<n;j++){
      so[k] = s[j];
      k += 1;
    }
    for(i=0;i<o->n;i++){
      for(j=0;j<cn[i];j++){
	so[k] = cs[i][j];
	k += 1;
      }
      myfree(cs[i]);
    }
    myfree(cs);
    myfree(cn);
  }

  *rn = tn;
  return so;
} 
/**************************************-**************************************/
/*                                                                           */
/*                       PARAMFILE_GET_ONODE_FROM_CARRAY                     */
/*                                                                           */
/*  Return an onode structure built from the carray.                         */
/*                                                                           */
/*****************************************************************************/
struct onode *paramfile_get_onode_from_carray(c)
     char **c;  // Address of pointer
{
  int i;
  char *t,*tstr;
  struct onode *op,*ot;

  op = paramfile_onode_create_onode();

  t = *c;
  if (*t != '\0'){
    op->otype = strdup(t);
    t += strlen(op->otype);
  }
  t += 1;

  tstr = strdup(t);
  op->resolve = atoi(tstr);  // Resolve flag
  t += 1 + strlen(tstr);
  myfree(tstr);


  //
  //
  //  WYETH - this is a new block of code, for TAG VALUES
  //  April 26, 2008
  //
  tstr = strdup(t); // Number of tag assignments
  op->tv_n = atoi(tstr);
  t += 1 + strlen(tstr);
  myfree(tstr);

  op->tv_name = (char **)myalloc(op->tv_n * sizeof(char *));
  op->tv      = (char **)myalloc(op->tv_n * sizeof(char *));
  for(i=0;i<op->tv_n;i++){
    op->tv_name[i] = strdup(t);
    t += strlen(op->tv_name[i]);
    t += 1;

    op->tv[i] = strdup(t);
    t += strlen(op->tv[i]);
    t += 1;
  }
  //
  //  WYETH - END NEW BLOCK
  //



  if (*t != '\0'){
    op->otag = strdup(t);
    t += strlen(op->otag);
  }
  t += 1;
  
  tstr = strdup(t);
  op->n = atoi(tstr);
  t += 1 + strlen(tstr);
  myfree(tstr);

  if (*t != '\0'){
    op->name = strdup(t);
    t += strlen(op->name);
  }
  t += 1;
  if (*t != '\0'){
    op->val = strdup(t);
    t += strlen(op->val);
  }
  t += 1;
  if (*t != '\0'){
    op->unit = strdup(t);
    t += strlen(op->unit);
  }
  t += 1;
  if (*t != '\0'){
    op->comm = strdup(t);
    t += strlen(op->comm);
  }
  t += 1;

  // WYETH NEW nval is computed from the value string.
  if (op->val != NULL)
    op->nval =  count_items_in_string(op->val);


  *c = t;

  if (op->n > 0){
    op->o = paramfile_get_onode_from_carray(c);
    ot = op->o;
    ot->prev = NULL;
    for(i=1;i<op->n;i++){
      ot->next = paramfile_get_onode_from_carray(c);
      ot->next->prev = ot;
      ot = ot->next;
    }
  }

  return op;
} 
/**************************************-**************************************/
/*                                                                           */
/*                      PARAMFILE_ONODE_INSERT_NEXT_ONODE                    */
/*                                                                           */
/*  Insert 'onew' right after 'o' in a doubly linked list.                   */
/*                                                                           */
/*****************************************************************************/
void paramfile_onode_insert_next_onode(o,onew)
     struct onode *o;
     struct onode *onew;
{
  if (onew->next != NULL)
    exit_error("PARAMFILE_ONODE_INSERT_NEXT_ONODE","Next is not null");
  if (onew->prev != NULL)
    exit_error("PARAMFILE_ONODE_INSERT_NEXT_ONODE","Prev is not null");

  onew->prev = o;
  onew->next = o->next;
  if (o->next != NULL)
    o->next->prev = onew;

  o->next = onew;
}
/**************************************-**************************************/
/*                                                                           */
/*                           ONODE_INSERT_CHILD_AT_END                       */
/*                                                                           */
/*  Insert 'onew' at the end of the list of children of 'o'.                 */
/*                                                                           */
/*****************************************************************************/
void onode_insert_child_at_end(o,onew)
     struct onode *o;
     struct onode *onew;
{
  struct onode *t;

  if (o->n == 0){
    o->o = onew;
    o->n = 1;
  }else{
    // Advance 't' to the last child node
    t = o->o;
    while(t->next != NULL)
      t = t->next;

    paramfile_onode_insert_next_onode(t,onew);
    o->n += 1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                       PARAMFILE_ONODE_GET_INIT_ONODE                      */
/*                                                                           */
/*****************************************************************************/
struct onode *paramfile_onode_get_init_onode(otype,otag,name,val,nval,unit,
					     comm,cpflag)
     char *otype;
     char *otag;
     char *name;
     char *val;
     int nval;
     char *unit;
     char *comm;
     int cpflag;
{
  struct onode *o;

  o = paramfile_onode_create_onode();

  if (cpflag){
    o->otype = strdup(otype);
    if (otag != NULL)
      o->otag  = strdup(otag);
    o->name  = strdup(name);
    o->val   = strdup(val);
    if (unit != NULL)
      o->unit  = strdup(unit);
    if (comm != NULL)
      o->comm  = strdup(comm);
  }else{
    o->otype = otype;
    o->otag  = otag;
    o->name  = name;
    o->val   = val;
    o->unit  = unit;
    o->comm  = comm;
  }
  o->nval  = nval;

  return o;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 ONODE_COPY                                */
/*                                                                           */
/*  Return a copy of the onode (all children are copied recursively).        */
/*                                                                           */
/*****************************************************************************/
struct onode *onode_copy(o)
     struct onode *o;
{
  int tn;
  struct onode *t,*tt,*onew;

  t = paramfile_onode_create_onode();

  if (o->otype != NULL)  t->otype = strdup(o->otype);
  if (o->otag  != NULL)  t->otag  = strdup(o->otag);
  if (o->name  != NULL)  t->name  = strdup(o->name);
  if (o->val   != NULL)  t->val   = strdup(o->val);
  if (o->unit  != NULL)  t->unit  = strdup(o->unit);
  if (o->comm  != NULL)  t->comm  = strdup(o->comm);

  t->nval = o->nval;

  tn = o->tv_n;
  t->tv_n = tn;
  if (tn > 0){
    t->tv_name = copy_2d_carray(o->tv_name,tn);
    t->tv      = copy_2d_carray(o->tv     ,tn);
  }

  tt = o->o;
  while(tt != NULL){
    onew = onode_copy(tt);
    onode_insert_child_at_end(t,onew);
    tt = tt->next;
  }

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                              ONODE_TREE_COUNT                             */
/*                                                                           */
/*  Return the number of onodes found that satisfy the condition, including  */
/*  'o' and all onodes below it.                                             */
/*                                                                           */
/*****************************************************************************/
int onode_tree_count(o,cflag,cstr)
     struct onode *o;
     int cflag;         // What to count:
                        //   0 - all nodes
                        //   1 - item nodes
                        //   2 - non-item nodes
                        //   3 - node type specified in 'cstr'
                        //   4 - node tag assignment, lhs specified in 'cstr'
                        //   5 - resolve > 0
     char *cstr;        // string value associated w/ cflag
{
  int n;
  struct onode *t;

  if (cflag == 0)
    n = 1;
  else if ((cflag == 1) && (strcmp(o->otype,"item")==0))
    n = 1;
  else if ((cflag == 2) && (strcmp(o->otype,"item")!=0))
    n = 1;
  else if ((cflag == 3) && (strcmp(o->otype,cstr)==0))
    n = 1;
  else if ((cflag == 4) && (search_2d_carray(o->tv_name,cstr,o->tv_n) >= 0))
    n = 1;
  else if ((cflag == 5) && (o->resolve > 0))
    n = 1;
  else
    n = 0;  // Do not count this node

  t = o->o;
  while(t != NULL){
    n += onode_tree_count(t,cflag,cstr);
    t = t->next;
  }

  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                            ONODE_TREE_FILL_LIST                           */
/*                                                                           */
/*  Fill the list 'olist', which must already be the appropriate length,     */
/*  with pointers to all onodes that satisfy the condition.                  */
/*                                                                           */
/*****************************************************************************/
void onode_tree_fill_list(o,cflag,cstr,olist,rk)
     struct onode *o;
     int cflag;         // What to include:
                        //   0 - all nodes
                        //   1 - item nodes
                        //   2 - non-item nodes
                        //   3 - node type specified in 'cstr'
                        //   4 - node tag assignment, lhs specified in 'cstr'
                        //   5 - resolve > 0
     char *cstr;        // string value associated w/ cflag
     struct onode **olist;
     int *rk;           // Next pointer in list to fill
{
  int flag;
  struct onode *t;

  if (cflag == 0)
    flag = 1;
  else if ((cflag == 1) && (strcmp(o->otype,"item")==0))
    flag = 1;
  else if ((cflag == 2) && (strcmp(o->otype,"item")!=0))
    flag = 1;
  else if ((cflag == 3) && (strcmp(o->otype,cstr)==0))
    flag = 1;
  else if ((cflag == 4) && (search_2d_carray(o->tv_name,cstr,o->tv_n) >= 0))
    flag = 1;
  else if ((cflag == 5) && (o->resolve > 0))
    flag = 1;
  else
    flag = 0;  // Do not count this node

  if (flag){
    olist[*rk] = o;
    *rk += 1;
  }

  t = o->o;
  while(t != NULL){
    onode_tree_fill_list(t,cflag,cstr,olist,rk);
    t = t->next;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            ONODE_TREE_GET_LIST                            */
/*                                                                           */
/*  Fill the list 'olist', which must already be the appropriate length,     */
/*  with pointers to all onodes that satisfy the condition.                  */
/*                                                                           */
/*****************************************************************************/
struct onode **onode_tree_get_list(o,cflag,cstr,rn)
     struct onode *o;
     int cflag;         // What to include:
                        //   0 - all nodes
                        //   1 - item nodes
                        //   2 - non-item nodes
                        //   3 - node type specified in 'cstr'
                        //   4 - node tag assignment, lhs specified in 'cstr'
     char *cstr;        // string value associated w/ cflag
     int *rn;           // Number of onodes in list
{
  int i;
  int nid,tn;
  struct onode **olist;

  nid = onode_tree_count(o,cflag,cstr);
  olist = (struct onode **)myalloc(nid*sizeof(struct onode *));
  tn = 0;
  onode_tree_fill_list(o,cflag,cstr,olist,&tn);

  *rn = nid;
  return olist;
}
/**************************************-**************************************/
/*                                                                           */
/*                         PARAMFILE_ONODE_EXTRACT_ITEM                      */
/*                                                                           */
/*****************************************************************************/
struct onode *paramfile_onode_extract_item(data,n,ops)
     char **data;  // List of lines [n] from the infile
     int n;        // Number of lines in file
     struct onode_parse_state *ops;
{
  int k;
  int nval;
  struct onode *o;
  char *otype,*otag,*tname,*tval,*tunit,*tcomm;

  //  WYETH HERE  - major changes 20 Jan 2009 to allow multiple valued items
  //  WYETH HERE  - major changes 20 Jan 2009 to allow multiple valued items
  
  if (ops->ns < 2){
    printf("  %s\n",ops->slist[0]);
    exit_error("PARAMFILE_ONODE_EXTRACT_ITEM","No value for item");
  }else{
    tname = strdup(ops->slist[ops->i]);
    //tval  = strdup(ops->slist[ops->i+1]); // Wait, incase multiple values
    tunit = NULL;
    tcomm = NULL;
    
    k = ops->i+2;
    nval = 1;
    while(k < ops->ns){
      
      if (ops->slist[k][0] == '('){  // Units
	tunit  = strdup(ops->slist[k]);
	
	if (ops->ns > (k+1)){
	  if (ops->slist[k+1][0] == '#'){  // Comment
	    tcomm = make_string_from_items(&(ops->slist[k+1]),ops->ns-(k+1));
	    k = ops->ns; // End the loop
	  }else{
	    exit_error("PARAMFILE_ONODE_EXTRACT_ITEM","Error after units");
	  }
	}
	k = ops->ns; // End the loop;   THIS LINE added 29 JAN 2009; 
      }else if (ops->slist[k][0] == '#'){  // Comment
	tcomm = make_string_from_items(&(ops->slist[k]),ops->ns-k);
	k = ops->ns; // End the loop
      }else{
	k += 1;
	nval += 1;  // Multiple values for this parameter name
      }
    }
    
    if (nval == 1)
      tval = strdup(ops->slist[ops->i+1]);
    else{
      tval = make_string_from_items(&(ops->slist[1]),nval);
    }
    
    otype = strdup("item");
    otag = NULL;
    o = paramfile_onode_get_init_onode(otype,otag,tname,tval,nval,tunit,
				       tcomm,0);
  }
  return o;
}
/**************************************-**************************************/
/*                                                                           */
/*                           PARAMFILE_ONODE_READ                            */
/*                                                                           */
/*  Recursively read and build a list of OBJECTs from 'data' strings.        */
/*  The pointer is returned for the object that contains the list.           */
/*                                                                           */
/*  OBJECT is any one of the following or a list of OBJECTS                  */
/*                                                                           */
/*    1.  name value                                                         */
/*    2.  name value (units)                                                 */
/*    3.  name value (units)  # Comments                                     */
/*    4.  name value          # Comments                                     */
/*    5.  [tag]                                                              */
/*    6.  <otype [otag]>                                                     */
/*         OBJECT                                                            */
/*         ...                                                               */
/*        </otype>  OR  <>                                                   */
/*                                                                           */
/*****************************************************************************/
struct onode *paramfile_onode_read(data,n,ops)
     char **data;  // List of lines [n] from the infile
     int n;        // Number of lines in file
     struct onode_parse_state *ops;  // State includes:  i,j,ns,slist
{
  int i;
  int fflag = 1;
  int done,pflag,nass,resflag;
  char *tag,**aname,**aval;
  struct onode *o,*ot,*olast;
  
  pflag = 0;
  
  o = paramfile_onode_create_onode();  // Get empty node
  olast = NULL;                        // Last node added to linked list
  
  done = 0;
  while(!done){

    if (ops->i >= ops->ns){
      if (ops->ns > 0)
	free_2d_carray(ops->slist,ops->ns);
      ops->j += 1;

      if (ops->j < n){
	get_items_from_string(data[ops->j],&(ops->slist),&(ops->ns));
	ops->i = 0; // Start at beginning of list
      }else{
	done = 1;
      }
    }else if (ops->ns > 0){
      // Skip comments (#), HTML commands (<), or other tags

      if (fflag){
	if (pflag)
	  printf("PARAMFILE_ONODE_READ ---- starting w/ slist[0] = %s\n",
		 ops->slist[0]);
	fflag = 0;
      }

      if (ops->slist[ops->i][0] == '#'){
	ops->i = ops->ns; // Skip rest of line
      }else if (ops->slist[ops->i][0] == '<'){
	if (ops->slist[ops->i][1] == '>'){
	  done = 1;
	  if (pflag) printf("  END TAG:  <>\n");
	  ops->i = ops->ns; // Skip rest of line
	}else if (ops->slist[ops->i][1] == '/'){
	  done = 1;
	  tag = myxml_get_tag_from_string(ops->slist[ops->i]);
	  if (pflag) printf("  END TAG:  %s\n",tag);
	  /*** WYETH - should check for appropriate name ***/
	  ops->i = ops->ns; // Skip rest of line
	}else{

	  // Make nodelist for this new tag
	  // Note, 'tag' is set to lower case, new storage is returned
	  tag = myxml_get_tag_from_string(ops->slist[ops->i]);
	  if (pflag) printf("    tag=  %s\n",tag);

	  if (strcmp(tag,"include")==0){
	    if (pflag) printf("    ignoring INCLUDE\n");
	    ops->i = ops->ns; /* Skip rest of line */
	    myfree(tag); // WYETH - added 26 April, 2008
	  }else{

	    // Get the tag assignments, and check for issues to resolve
	    nass = myxml_get_tag_assignments(data[ops->j],&aname,&aval);
	    resflag = 0;
	    for(i=0;i<nass;i++){
	      if (strcmp(aname[i],"COPY")==0)
		resflag = 1;
	    }

	    ops->i = ops->ns; // Skip rest of line

	    ot = paramfile_onode_read(data,n,ops);
	    /*ot->name = tag;*/
	    ot->otype = tag;

	    ot->tv_n    = nass;
	    ot->tv_name = aname;
	    ot->tv      = aval;
	    ot->resolve = resflag;  // There are issues to resolve
	    
	    
	    if (pflag) printf("    inserting list node  %s\n",ot->otype);
	    if (o->n == 0){
	      o->o = ot;
	      o->n = 1;
	      olast = ot;
	    }else{
	      paramfile_onode_insert_next_onode(olast,ot);
	      olast = ot;
	      o->n += 1;
	    }
	    if (pflag) printf("      done insert\n");
	    if (pflag) printf("        i= %d  j= %d\n",ops->i,ops->j);
	  }
	}
      }else if (ops->slist[ops->i][0] == '['){
	ops->i += 1;
      }else{ // name of an item

	ot = paramfile_onode_extract_item(data,n,ops);

	if (pflag) printf("    inserting item node  %s\n",ot->name);

	if (o->n == 0){
	  o->o = ot;
	  o->n = 1;
	  olast = ot;
	}else{
	  paramfile_onode_insert_next_onode(olast,ot);
	  olast = ot;
	  o->n += 1;
	}
	ops->i = ops->ns; // Ignore the rest of the line
      }
    }
  }

  return o;
}
/**************************************-**************************************/
/*                                                                           */
/*                         PARAMFILE_ONODE_FILE_READ                         */
/*                                                                           */
/*****************************************************************************/
struct onode *paramfile_onode_file_read(infile,pflag)
     char infile[];
     int pflag;
{
  int i,j;
  int n;
  char **data;
  struct onode *o;
  struct onode_parse_state *ops;

  if (pflag){
    printf("  PARAMFILE_ONODE_FILE_READ\n");
    printf("    Reading parameters from %s\n",infile);
  }

  read_2d_carray(infile,1,0,&data,&n);
  if (n == -1){
    printf("%s\n",infile);
    exit_error("PARAMFILE_ONODE_FILE_READ","Cannot open file");
  }

  ops = (struct onode_parse_state *)myalloc(sizeof(struct onode_parse_state));
  ops->j = -1;  // WYETH - July 2009, to read first line of file
  ops->i = 0;
  ops->ns = 0;
  ops->slist = NULL;

  i = j = 0;

  o = paramfile_onode_read(data,n,ops);

  o->otype = strdup("TOP");

  //paramfile_onode_print(o,0,0);

  return o;
}
/**************************************-**************************************/
/*                                                                           */
/*                         PARAMFILE_XML_SKIP_COMMENT                        */
/*                                                                           */
/*  Skip over any sequence of comments with optional white space before and  */
/*  in between.                                                              */
/*                                                                           */
/*****************************************************************************/
void paramfile_xml_skip_comment(data,n,ops)
     char **data;  // List of lines [n] from infile, end in '/0' (not '/n')
     int n;        // Number of lines in file
     struct onode_parse_state *ops;  // State includes:  i,j,ns,slist
{
  int i,j;
  int done,alldone,com_start,com_done,pflag;
  char *tstr;

  pflag = 0;

  //printf("HERE 0\n");

  j = ops->j;  // Line number
  i = ops->i;  // Char number
  tstr = data[j];

  com_start = com_done = 0;

  alldone = 0;
  while(!alldone){

    // Skip all white space until '<'
    done = 0;
    while(!done){
      if (i >= strlen(tstr)){
	i = 0;
	j += 1;
	if (j >= n){
	  done = 1;
	  alldone = 1;
	}else
	  tstr = data[j];

      }else if ((tstr[i] == (char)9) ||   // tab
		(tstr[i] == (char)32) ||  // space
		(tstr[i] == '\n')){       // newline
	i += 1;
      }else{
	done = 1;
	//printf("HERE done 1\n");
      }
    }

    //printf("HERE 2\n");

    if (!alldone){
      
      if ((i+3) >= strlen(tstr))
	alldone = 1;  // This cannot be a comment because the line ends
      else if ((tstr[i] != '<') || (tstr[i+1] != '!') || (tstr[i+2] != '-') ||
	       (tstr[i+3] != '-')){
	alldone = 1;  // This is not a comment start symbol

	//printf("HERE 3 - alldone, not a comment start\n");

      }else{
	// This is a comment start, skip to the end

	//printf("HERE 4 - this is a comment start\n");

	com_start += 1;  // A comment has started
	com_done = 0;

	i += 3;
	done = 0;
	while(!done){
	  //printf("i %d   ==>%s<==\n",i,tstr);
	  if (i >= strlen(tstr)){
	    i = 0;
	    j += 1;
	    if (j >= n){
	      done = 1;
	      alldone = 1;
	    }else
	      tstr = data[j];

	  }else if (tstr[i] == '>'){
	    if (i >= 2){
	      if ((tstr[i-2] == '-') && (tstr[i-1] == '-')){
		done = 1;
		com_done = 1;  // A comment has ended
		//printf("COMM -------- END DONE!!!\n");
	      }
	    }
	    i += 1;
	  }else{
	    i += 1;
	  }
	}
      }
    }
  }

  if (com_start >= 1){
    if (com_done == 1){
      ops->i = i;
      ops->j = j;
      if (pflag) printf("Skipped %d comments\n",com_start);
    }else{
      exit_error("PARAMFILE_XML_SKIP_COMMENT","Comment never ended");
    }
  }
  // If there was no comment, do not advance the pointers.
    
}
/**************************************-**************************************/
/*                                                                           */
/*                          PARAMFILE_ONODE_XML_READ                         */
/*                                                                           */
/*  Recursively read and build a list of OBJECTs from 'data' strings         */
/*  which represents an XML file.                                            */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/
struct onode *paramfile_onode_xml_read(data,n,ops)
     char **data;  // List of lines [n] from infile, end in '/0' (not '/n')
     int n;        // Number of lines in file
     struct onode_parse_state *ops;  // State includes:  i,j,ns,slist
{
  int i,k;
  int done,pflag,dk,nonwhite,closedflag;
  char ttag[LONG_SLEN];
  char *tstr,*tend,dstr[SLEN_MAX];
  struct onode *o,*ot,*olast;
  
  pflag = 0;

  //
  //  We are at an '<' that opens a new node, create new empty node:
  //
  o = paramfile_onode_create_onode();  // Get empty node
  olast = NULL;                        // Last node added to linked list

  //
  // Store our tag in 'otype', advance to the tag close '>'
  //
  paramfile_xml_skip_comment(data,n,ops);  // Skip over one or more comments
  o->otype = myxml_get_tag_from_string_adv_index(data[ops->j],&(ops->i));
  ops->i += 1;  // Skip over closing '>'

  if (pflag) printf("otype ==>%s<===\n",o->otype);

  nonwhite = 0;    // Haven't found any non-white data yet
  closedflag = 0;  // We just passed an opening tag
  dk = 0;          // data counter is 0

  i = ops->i;
  tstr = data[ops->j];

  done = 0;
  while(!done){
    
    if (pflag) printf("HERE 0.  i= %d  j= %d\n",i,ops->j);
    
    //  ops->i   - index of character on line
    //  ops->j   - index of line in file
    
    if (i >= strlen(tstr)){  // if at last character
      ops->j += 1;  // go to next line
      ops->i = 0;   // back to first char
      
      if (ops->j >= n){
	done = 1;
      }else{
	tstr = data[ops->j];  // Current line
	i = ops->i;           // index of current char
      }

      if (pflag) printf("HERE 1.   i= %d   j= %d   n= %d\n",i,ops->j,n);

    }else{

      if (pflag) printf("HERE 2.\n");
      
      if (tstr[i] == '<'){  // Tag or end tag
	if (tstr[i+1] == '/'){  // End tag
	  if (pflag) printf("END TAG\n");

	  i += 1;  // Point to '/'

	  tend = myxml_get_tagend_from_string_adv_index(tstr,&i);
	  if (strcmp(tend,o->otype)!=0){
	    printf("  *** ==>%s<== does not match ==>%s<==\n",tend,o->otype);
	    exit_error("wyeth","Tag end does not match");
	  }

	  if (pflag) printf("Found matching endtag:  %s\n",tend);
	  if (nonwhite == 1){  // There was data content
	    dstr[dk] = '\0';
	    o->val = strdup(dstr);  // Store the content
	  }

	  if (tend != NULL)
	    myfree(tend); // WYETH - Added 11 April 2010, try to fix mem leak
	  done = 1;  // Return current node

	}else{ // New tag
	  
	  // New start - make recursive call

	  if (pflag) printf("NEW start - recurse\n");

	  ops->i = i;
	  ot = paramfile_onode_xml_read(data,n,ops);
	  i = ops->i;

	  ot->tv_n    = 0;
	  ot->tv_name = NULL;
	  ot->tv      = NULL;
	  ot->resolve = 0;    // There are issues to resolve
	    
	  if (pflag) printf("    inserting list node  %s\n",ot->otype);
	  if (o->n == 0){
	    o->o = ot;
	    o->n = 1;
	    olast = ot;
	  }else{
	    paramfile_onode_insert_next_onode(olast,ot);
	    olast = ot;
	    o->n += 1;
	  }
	  if (pflag) printf("      done insert\n");
	  if (pflag) printf("        i= %d  j= %d\n",ops->i,ops->j);

	  closedflag = 0; // Tag just closed
	
	  //exit(0);
	}
      }else if (closedflag == 1){
	i += 1;
      }else if (closedflag == 0){
	//
	//  Treat this as potential content within a tag
	//
	if ((tstr[i] == (char)9) ||   // tab
	    (tstr[i] == (char)32) ||  // space
	    (tstr[i] == '\n')){       // newline
	  ; // White space
	}else{
	  nonwhite = 1;
	}

	if (pflag) printf("  w-here i=%d  %s\n",i,tstr);

	dstr[dk] = tstr[i];
	dk += 1;
	if (dk >= SLEN_MAX)
	  exit_error("wyeth","Data exceeds arbitrary limit: 16384");
	i += 1;  // Move to next char
	
      }
    }
  }

  ops->i = i;

  return o;
}
/**************************************-**************************************/
/*                                                                           */
/*                        PARAMFILE_ONODE_XML_FILE_READ                      */
/*                                                                           */
/*  Read a simple XML file.                                                  */
/*                                                                           */
/*****************************************************************************/
struct onode *paramfile_onode_xml_file_read(infile,pflag)
     char infile[];
     int pflag;
{
  int n;
  char **data;
  struct onode *o;
  struct onode_parse_state *ops;

  if (pflag){
    printf("  PARAMFILE_ONODE_XML_FILE_READ\n");
    printf("    Reading parameters from %s\n",infile);
  }

  read_2d_carray(infile,1,0,&data,&n);
  if (n == -1){
    printf("%s\n",infile);
    exit_error("PARAMFILE_ONODE_XML_FILE_READ","Cannot open file");
  }
  
  ops = (struct onode_parse_state *)myalloc(sizeof(struct onode_parse_state));
  ops->j = 0;  // WYETH - July 2009, to read first line of file
  ops->i = 0;
  ops->ns = 0;
  ops->slist = NULL;

  o = paramfile_onode_xml_read(data,n,ops);

  o->otype = strdup("TOP");

  if (pflag)
    paramfile_onode_print(o,0,0);

  if (n > 0){ // WYETH added 11 April 2010, trying to solve mcontrol mem leak.
    free_2d_carray(data,n);
  }

  return o;
}
/**************************************-**************************************/
/*                                                                           */
/*                            ONODE_GET_TAG_VAL_PTR                          */
/*                                                                           */
/*****************************************************************************/
char *onode_get_tag_val_ptr(o,vname)
     struct onode *o;
     char *vname;
{
  int k;

  k = search_2d_carray(o->tv_name,vname,o->tv_n);

  if (k < 0)
    return NULL;
  else
    return o->tv[k];

}
/**************************************-**************************************/
/*                                                                           */
/*                            ONODE_GET_NEXT_TYPE                            */
/*                                                                           */
/*****************************************************************************/
struct onode *onode_get_next_type(o,otype)
     struct onode *o;
     char *otype;
{
  int done;
  struct onode *t,*tr;

  tr = NULL;
  t = o;
  done = 0;
  while(!done){
    if (t == NULL)
      done = 1;
    else if (strcmp(t->otype,otype)==0){
      tr = t;
      done = 1;
    }else{
      t = t->next;
    }
  }

  return tr;
}
/**************************************-**************************************/
/*                                                                           */
/*                            ONODE_GET_ITEM_CHILD                           */
/*                                                                           */
/*  Of all children of otype 'item', return the one having 'name',           */
/*  or NULL if not found.                                                    */
/*                                                                           */
/*****************************************************************************/
struct onode *onode_get_item_child(o,name)
     struct onode *o;
     char *name;
{
  struct onode *t,*s;

  s = NULL;
  t = onode_get_next_type(o->o,"item");
  while(t != NULL){
    if (strcmp(t->name,name)==0){
      s = t;
      t = NULL;
    }else
      t = onode_get_next_type(t->next,"item");
  }

  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                                ONODE_ITEM                                 */
/*                                                                           */
/*  Return 1 if item child is found, 0 otherwise.                            */
/*                                                                           */
/*****************************************************************************/
int onode_item(o,name)
     struct onode *o;
     char *name;
{
  if (onode_get_item_child(o,name) == NULL)
    return 0;
  else
    return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                           ONODE_CHILD_GET_UNIQUE                          */
/*                                                                           */
/*  Return the unique child having 'otype', NULL if none, exit if multiple.  */
/*                                                                           */
/*****************************************************************************/
struct onode *onode_child_get_unique(o,otype)
     struct onode *o;
     char *otype;
{
  struct onode *t,*s;

  s = onode_get_next_type(o->o,otype);
  if (s != NULL){
    t = onode_get_next_type(s->next,otype);
    if (t != NULL){
      printf("  otype %s\n",otype);
      exit_error("ONODE_CHILD_GET_UNIQUE","Multiple children of otype");
    }
  }

  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                             ONODE_GET_NTH_VAL                             */
/*                                                                           */
/*  Storage is created and returned, must be freed.                          */
/*                                                                           */
/*****************************************************************************/
char *onode_get_nth_val(o,n)
     struct onode *o;
     int n;             // Return the 'n' th value, 0,1,...
{
  int ns;
  char **slist,*tval;

  if (n >= o->nval)
    exit_error("ONODE_GET_NTH_VAL","Not enough values");

  get_items_from_string(o->val,&slist,&ns);
  if (ns != o->nval)
    exit_error("ONODE_GET_NTH_VAL","Stored 'n' error");

  tval = strdup(slist[n]);

  free_2d_carray(slist,ns);

  return tval;
}
/**************************************-**************************************/
/*                                                                           */
/*                             ONODE_GETPAR_NVALS                            */
/*                                                                           */
/*  Return the number of values for the parameter 'name' in 'o'.             */
/*                                                                           */
/*****************************************************************************/
int onode_getpar_nvals(o,name)
     struct onode *o;
     char *name;
{
  struct onode *t;
  
  t = onode_get_item_child(o,name);
  if (t == NULL){
    printf("*** param name %s\n",name);
    exit_error("ONODE_GETPAR_NVALS","Parameter not found");
  }

  return t->nval;
}
/**************************************-**************************************/
/*                                                                           */
/*                           ONODE_GETPAR_INT_DFLT                           */
/*                                                                           */
/*****************************************************************************/
int onode_getpar_int_dflt(o,name,defval)
     struct onode *o;
     char *name;
     int defval;
{
  int k;
  struct onode *t;
  
  t = onode_get_item_child(o,name);
  if (t == NULL)
    k = defval;
  else{
    k = atoi(t->val);
  }

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                           ONODE_GETPAR_INT_EXIT                           */
/*                                                                           */
/*****************************************************************************/
int onode_getpar_int_exit(o,name)
     struct onode *o;
     char *name;
{
  int k;
  struct onode *t;
  
  t = onode_get_item_child(o,name);
  if (t == NULL){
    printf("*** param name %s\n",name);
    exit_error("ONODE_GETPAR_INT_EXIT","Parameter not found");
  }else{
    k = atoi(t->val);
  }

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                           ONODE_GETPAR_INT_LIST                           */
/*                                                                           */
/*****************************************************************************/
int *onode_getpar_int_list(o,name,rn)
     struct onode *o;
     char *name;
     int *rn;
{
  int i;
  int *ival,ns;
  char **slist;
  struct onode *t;
  char *tval;
  
  t = onode_get_item_child(o,name);
  if (t == NULL){
    printf("*** param name %s\n",name);
    exit_error("ONODE_GETPAR_INT_LIST","Parameter not found");
  }else{

    get_items_from_string(t->val,&slist,&ns);
    if (ns != t->nval){
      printf("name = %s,  t->val = %s\n",name,t->val);
      printf("ns = %d  t->nval = %d\n",ns,t->nval);
      exit_error("ONODE_GETPAR_INT_LIST","Stored 'n' error");
    }

    ival = (int *)myalloc(ns*sizeof(int));
    for(i=0;i<ns;i++)
      ival[i] = atoi(slist[i]);

    free_2d_carray(slist,ns);
  }

  *rn = ns;
  return ival;
}
/**************************************-**************************************/
/*                                                                           */
/*                           ONODE_GETPAR_FLT_DFLT                           */
/*                                                                           */
/*****************************************************************************/
float onode_getpar_flt_dflt(o,name,defval)
     struct onode *o;
     char *name;
     float defval;
{
  float k;
  struct onode *t;
  
  t = onode_get_item_child(o,name);
  if (t == NULL)
    k = defval;
  else{
    k = atof(t->val);
  }

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                           ONODE_GETPAR_FLT_EXIT                           */
/*                                                                           */
/*****************************************************************************/
float onode_getpar_flt_exit(o,name)
     struct onode *o;
     char *name;
{
  float x;
  struct onode *t;
  
  t = onode_get_item_child(o,name);
  if (t == NULL){
    printf("*** param name %s\n",name);
    exit_error("ONODE_GETPAR_FLT_EXIT","Parameter not found");
  }else{
    x = atof(t->val);
  }
  return x;
}
/**************************************-**************************************/
/*                                                                           */
/*                         ONODE_GETPAR_NTH_FLT_EXIT                         */
/*                                                                           */
/*****************************************************************************/
float onode_getpar_nth_flt_exit(o,name,n)
     struct onode *o;
     char *name;
     int n;
{
  float x;
  struct onode *t;
  char *tval;
  
  t = onode_get_item_child(o,name);
  if (t == NULL){
    printf("*** param name %s\n",name);
    exit_error("ONODE_GETPAR_FLT_EXIT","Parameter not found");
  }else{
    tval = onode_get_nth_val(t,n);
    x = atof(tval);
    myfree(tval);
  }
  return x;
}
/**************************************-**************************************/
/*                                                                           */
/*                           ONODE_GETPAR_FLT_LIST                           */
/*                                                                           */
/*****************************************************************************/
float *onode_getpar_flt_list(o,name,rn)
     struct onode *o;
     char *name;
     int *rn;
{
  int i;
  int ns;
  float *fval;
  char **slist,*tval;
  struct onode *t;
  
  t = onode_get_item_child(o,name);
  if (t == NULL){
    printf("*** param name %s\n",name);
    exit_error("ONODE_GETPAR_FLT_LIST","Parameter not found");
  }else{

    get_items_from_string(t->val,&slist,&ns);
    if (ns != t->nval){
      printf("name = %s,  t->val = %s\n",name,t->val);
      printf("ns = %d  t->nval = %d\n",ns,t->nval);
      exit_error("ONODE_GETPAR_FLT_LIST","Stored 'n' error");
    }

    fval = (float *)myalloc(ns*sizeof(float));
    for(i=0;i<ns;i++)
      fval[i] = atof(slist[i]);

    free_2d_carray(slist,ns);
  }

  *rn = ns;
  return fval;
}
/**************************************-**************************************/
/*                                                                           */
/*                           ONODE_GETPAR_DBL_EXIT                           */
/*                                                                           */
/*****************************************************************************/
double onode_getpar_dbl_exit(o,name)
     struct onode *o;
     char *name;
{
  double x;
  struct onode *t;
  
  t = onode_get_item_child(o,name);
  if (t == NULL){
    printf("*** param name %s\n",name);
    exit_error("ONODE_GETPAR_DBL_EXIT","Parameter not found");
  }else{
    x = atof(t->val);  // 'atof' returns double
  }
  return x;
}
/**************************************-**************************************/
/*                                                                           */
/*                           ONODE_GETPAR_CHR_DFLT                           */
/*                                                                           */
/*****************************************************************************/
char *onode_getpar_chr_dflt(o,name,defval)
     struct onode *o;
     char *name;
     char *defval;
{
  struct onode *t;
  char *s;
  
  t = onode_get_item_child(o,name);
  if (t == NULL)
    if (defval == NULL)
      s = NULL;
    else
      s = strdup(defval);
  else{
    s = strdup(t->val);
  }
  
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                           ONODE_GETPAR_CHR_NULL                           */
/*                                                                           */
/*  Get the char param value, but convert "NULL" to (char *)NULL.  Return    */
/*  NULL if param not found                                                  */
/*                                                                           */
/*****************************************************************************/
char *onode_getpar_chr_null(o,name)
     struct onode *o;
     char *name;
{
  struct onode *t;
  char *s;
 
  t = onode_get_item_child(o,name);
  if (t == NULL)
    s = NULL;
  else{
    if (strcmp(t->val,"NULL")==0)
      s = NULL;
    else
      s = strdup(t->val);
  }
  
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                           ONODE_GETPAR_CHR_EXIT                           */
/*                                                                           */
/*****************************************************************************/
char *onode_getpar_chr_exit(o,name)
     struct onode *o;
     char *name;
{
  struct onode *t;
  char *s;
  
  t = onode_get_item_child(o,name);
  if (t == NULL){
    printf("*** %s\n",name);
    exit_error("ONODE_GETPAR_CHR_EXIT","Parameter not found");
  }else{
    s = strdup(t->val);
  }
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                            ONODE_GETPAR_CHR_PTR                           */
/*                                                                           */
/*****************************************************************************/
char *onode_getpar_chr_ptr(o,name)
     struct onode *o;
     char *name;
{
  struct onode *t;
  char *s;
  
  s = NULL;
  t = onode_get_item_child(o,name);
  if (t != NULL)
    s = t->val;
  
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                         ONODE_GETPAR_CHR_PTR_EXIT                         */
/*                                                                           */
/*****************************************************************************/
char *onode_getpar_chr_ptr_exit(o,name)
     struct onode *o;
     char *name;
{
  char *s;
  
  s = onode_getpar_chr_ptr(o,name);
  if (s == NULL){
    printf("*** %s\n",name);
    exit_error("ONODE_GETPAR_CHR_PTR_EXIT","Parameter not found");
  }

  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                      ONODE_GETPAR_CHILD_CHR_PTR_EXIT                      */
/*                                                                           */
/*  Return a pointer to the value of the item "parname" in the unique        */
/*  child "childname" of "o".                                                */
/*                                                                           */
/*****************************************************************************/
char *onode_getpar_child_chr_ptr_exit(o,childname,parname)
     struct onode *o;
     char *childname;   // Find this child in "o", or exit
     char *parname;     // Find this parameter in the child, or exit
{
  char *s;
  struct onode *oc;

  oc = onode_child_get_unique(o,childname);
  if (oc == NULL){
    printf("*** %s\n",childname);
    exit_error("ONODE_GETPAR_CHILD_CHR_PTR_EXIT","Child node not found");
  }
  
  s = onode_getpar_chr_ptr(oc,parname);
  if (s == NULL){
    printf("*** %s\n",parname);
    exit_error("ONODE_GETPAR_CHILD_CHR_PTR_EXIT","Parameter not in child");
  }

  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                           ONODE_GETPAR_CHR_LIST                           */
/*                                                                           */
/*****************************************************************************/
char **onode_getpar_chr_list(o,name,rn)
     struct onode *o;
     char *name;
     int *rn;
{
  int ns;
  char **slist,*tval;
  struct onode *t;
  
  t = onode_get_item_child(o,name);
  if (t == NULL){
    printf("*** param name %s\n",name);
    exit_error("ONODE_GETPAR_FLT_LIST","Parameter not found");
  }else{

    get_items_from_string(t->val,&slist,&ns);
    if (ns != t->nval){
      printf("name = %s,  t->val = %s\n",name,t->val);
      printf("ns = %d  t->nval = %d\n",ns,t->nval);
      exit_error("ONODE_GETPAR_CHR_LIST","Stored 'n' error");
    }
  }

  *rn = ns;
  return slist;
}
/**************************************-**************************************/
/*                                                                           */
/*                             ONODE_GETPAR_UNIT                             */
/*                                                                           */
/*  This returns a POINTER, no storage is allocated.                         */
/*                                                                           */
/*****************************************************************************/
char *onode_getpar_unit(o,name)
     struct onode *o;
     char *name;
{
  struct onode *t;
  char *s;
  
  s = NULL;
  t = onode_get_item_child(o,name);
  if (t != NULL)
    s = t->unit;
  
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                              ONODE_COUNT_OTYPE                            */
/*                                                                           */
/*****************************************************************************/
int onode_count_otype(o,otype)
     struct onode *o;
     char *otype;
{
  int n;
  struct onode *t;

  n = 0;
  t = o->o;
  while(t!=NULL){
    if (strcmp(t->otype,otype)==0)
      n += 1;
    t = t->next;
  }
  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                         ONODE_GET_NODE_TYPE_ITEM_VAL                      */
/*                                                                           */
/*  Return the onode of 'otype' that has a child item w/ "name" and "val".   */
/*                                                                           */
/*****************************************************************************/
struct onode *onode_get_node_type_item_val(o,otype,item_name,item_val)
     struct onode *o;
     char *otype;
     char *item_name;
     char *item_val;
{
  struct onode *t,*oi,*of;

  of = NULL;
  t = onode_get_next_type(o->o,otype);
  while(t != NULL){
    oi = onode_get_item_child(t,item_name);
    if (oi != NULL){                      // has item w/ name
      if (strcmp(oi->val,item_val)==0){   // value is correct
	if (of == NULL){
	  of = t;
	}else{
	  exit_error("ONODE_GET_NODE_TYPE_ITEM_VAL","Multiple matches");
	}
      }
    }
    t = onode_get_next_type(t->next,otype);
  }

  return of;
}
/**************************************-**************************************/
/*                                                                           */
/*                           ONODE_PARSE_OSTR_NEXT                           */
/*                                                                           */
/*****************************************************************************/
char *onode_parse_ostr_next(s,ptype,rstr,rtype)
     char *s;
     int ptype;     // type of item previously parsed
     char **rstr;   // pointer to string
     int *rtype;    // 1-item, 2-otype, 3-otype_dot, 4-nameval
{
  int i;
  int done,stype;
  char *t,*tr;

  tr = NULL;

  i = 0;
  done = 0;
  while(!done){
    if (s[i] == '\0'){
      if (i==0)	exit_error("ONODE_PARSE_OSTR","Ostr format error, EOS");
      done = 1;
      t = strdup(s);
      stype = 1;
      tr = &(s[i]);
    }else if (s[i] == '.'){
      if (i==0)	exit_error("ONODE_PARSE_OSTR","Ostr format error, period");
      s[i] = '\0';
      t = strdup(s);
      s[i] = '.';
      stype = 3;
      done = 1;
      tr = &(s[i+1]);
    }else if (s[i] == '/'){
      if (i==0)	exit_error("ONODE_PARSE_OSTR","Ostr format error, (");
      s[i] = '\0';
      t = strdup(s);
      s[i] = '/';
      if (ptype == 3)
	stype = 4;
      else
	stype = 2;
      done = 1;
      tr = &(s[i+1]);
    }
    i += 1;
  }

  *rstr = t;
  *rtype = stype;

  return tr;
}
/**************************************-**************************************/
/*                                                                           */
/*                          ONODE_GET_ITEM_FOR_OSTR                          */
/*                                                                           */
/*  ostr examples                                                            */
/*                                                                           */
/*    tn                                                                     */
/*    pop.in/geometry/xn                                                     */
/*                                                                           */
/*****************************************************************************/
struct onode *onode_get_item_for_ostr(o,oname)
     struct onode *o;
     char *oname;
{
  int k;
  int stype,pflag;
  char *s,*t,*otype;
  struct onode *ot;

  pflag = 0;

  ot = o;
  stype = 0;
  t = oname;
  while(*t != '\0'){
    t = onode_parse_ostr_next(t,stype,&s,&stype);
    if (pflag) printf("0:  %s  %d\n",s,stype);

    if (stype == 1){
      ot = onode_get_item_child(ot,s);
      if (ot != NULL){
	if (pflag) printf("  1:  ot->otype = %s\n",ot->otype);
      }
    }else if (stype == 2){
      ot = onode_child_get_unique(ot,s);
      if (ot == NULL){
	return NULL;
      }
      if (pflag) printf("  2:  ot->otype = %s\n",ot->otype);
    }else if (stype == 3){
      otype = strdup(s);  // Save until name+val is parsed
    }else if (stype == 4){
      ot = onode_get_node_type_item_val(ot,otype,"name",s);
      if (ot != NULL){
	if (pflag) printf("  4:  ot->otype = %s\n",ot->otype);
      }else{
	return NULL;  // WYETH BUG FIX  2010 Dec 11
      }
      myfree(otype);
    }else{
      exit_error("ONODE_GET_ITEM_FOR_OSTR","Unknown stype");
    }
    if (pflag) printf("  5:  before free\n");
    myfree(s);
    if (pflag) printf("  6:  after free\n");
  }

  if (ot == NULL){
    if (pflag) printf("  %s not found.\n",oname);
  }else{
    if (pflag) printf("  Found, ot->name = %s\n",ot->name);
  }

  return ot;
}
/**************************************-**************************************/
/*                                                                           */
/*                              ONODE_TEST_OSTR                              */
/*                                                                           */
/*  Return 1 if 'ostr' describes a valid onode item.                         */
/*                                                                           */
/*****************************************************************************/
int onode_test_ostr(o,ostr)  // FOR COMMAND LINE ARGUMENT STRINGS
     struct onode *o;
     char *ostr;
{
  struct onode *ot;

  ot = onode_get_item_for_ostr(o,ostr);
  if (ot != NULL)
    return 1;
  else
    return 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                               ONODE_TEST_INT                              */
/*                                                                           */
/*  Return 1 if 'item' exists AND has the value 'ival'.                      */
/*                                                                           */
/*****************************************************************************/
int onode_test_int(o,item,ival)
     struct onode *o;
     char *item;         // item name
     int ival;           // item value
{
  int flag,k;
  struct onode *t;

  flag = 0;

  t = onode_get_item_child(o,item);
  if (t != NULL){
    k = atoi(t->val);
    if (k == ival)
      flag = 1;
  }

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                               ONODE_TEST_CHR                              */
/*                                                                           */
/*  Return 1 if 'item' exists AND has the value 'cval'.                      */
/*                                                                           */
/*****************************************************************************/
int onode_test_chr(o,item,cval)
     struct onode *o;
     char *item;         // item name
     char *cval;         // item value
{
  int flag;
  struct onode *t;

  flag = 0;

  t = onode_get_item_child(o,item);
  if (t != NULL){
    if (strcmp(t->val,cval)==0)
      flag = 1;
  }

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                              ONODE_MAKE_FNAME                             */
/*                                                                           */
/*  Designed after 'paramfile_make_fname'                                    */
/*                                                                           */
/*  Return (char *)NULL if "name" is not found.                              */
/*  If "name" is found and its value is                                      */
/*    == "*"  -->  return "prefix"+"extension"                               */
/*    != "*"  -->  return "name"                                             */
/*                                                                           */
/*****************************************************************************/
char *onode_make_fname(o,name,prefix,extension)
     struct onode *o;  // Onode
     char *name;       // Look up this parameter name
     char *prefix;     // If "name" is "*", then use the prefix and extension
     char *extension;  // file extension
{
  char *tname,tstr[SLEN];

  if (onode_test_ostr(o,name)){
    tname = onode_getpar_chr_exit(o,name);
    if (strcmp(tname,"*")==0)
      sprintf(tstr,"%s%s",prefix,extension);
    else
      sprintf(tstr,"%s",tname);
    myfree(tname);
  }else
    return (char *)NULL;

  return strdup(tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                             ONODE_UPDATE_VALUE                            */
/*                                                                           */
/*****************************************************************************/
int onode_update_value(o,ostr,new_val)
     struct onode *o;  // E.g., this may be the top model onode
     char *ostr;
     char *new_val;
{
  int flag;
  struct onode *ot;

  ot = onode_get_item_for_ostr(o,ostr);
  if (ot != NULL){
    if (ot->val != NULL)
      myfree(ot->val);
    ot->val = strdup(new_val);
    flag = 1;
  }else
    flag = 0;

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                           ONODE_UPDATE_FROM_SLIST                         */
/*                                                                           */
/*****************************************************************************/
int *onode_update_from_slist(o,mylogf,slist,ns)
     struct onode *o;  // E.g., this may be the top model onode
     char *mylogf;
     char **slist;
     int ns;
{
  int i;
  int *ulist,flag;
  char *tname,*tval;

  mylog(mylogf,"  ONODE_UPDATE_FROM_SLIST\n");

  if (ns > 0)
    ulist = get_zero_iarray(ns);  // Flag array to indicate which were used
  else
    ulist = NULL;

  for(i=0;i<ns;i+=2){
    if ((i+1) < ns){
      
      tname = slist[i];
      tval  = slist[i+1];

      flag = onode_update_value(o,tname,tval);
      if (flag == 1)
	ulist[i] = 1; // Set flag to indicate this param was updated
      // Else, the param may not have been for the moo file...

    }else{
      printf("*** WARNING last parameter name ignored---no value.\n");
    }
  }

  return ulist;
}
/**************************************-**************************************/
/*                                                                           */
/*                        ONODE_CONFLICT_WITH_CHILDREN                       */
/*                                                                           */
/*  A conflict exists if 'o' has the same 'otype' and name item as one of    */
/*  the children of 'po'.  If neither have a name item, this is the same     */
/*  has if they have the same name item.                                     */
/*                                                                           */
/*****************************************************************************/
int onode_conflict_with_children(po,o)
     struct onode *po;  // Check children of this parent node
     struct onode *o;   // to see if this node is a conflisource
{
  int conflict;
  char *name1,*name2;
  struct onode *t;

  conflict = 0;  // Assume no conflict
  t = po->o;
  while(t!=NULL){
    if (strcmp(t->otype,o->otype)==0){
      if (strcmp(t->otype,"item")==0){
	if (strcmp(t->name,o->name)==0)  // Items of same name
	  conflict = 1;
      }else{
	name1 = onode_getpar_chr_ptr(t,"name");
	name2 = onode_getpar_chr_ptr(o,"name");
	if ((name1 == NULL) && (name2 == NULL))  // same onodes w/o names
	  conflict = 1;
	else if (strcmp(name1,name2)==0)         // same onodes w/ same name
	  conflict = 1;
      }
    }
    t = t->next;
  }
  return conflict;
}
/**************************************-**************************************/
/*                                                                           */
/*                          ONODE_COPY_CONTENTS_NOVEL                        */
/*                                                                           */
/*  Copy any novel 'onodes' in 'os' to 'od'.                                 */
/*                                                                           */
/*****************************************************************************/
void onode_copy_contents_novel(od,os)
     struct onode *od;  // destination
     struct onode *os;  // source
{
  struct onode *t,*oc;

  t = os->o;
  while(t!=NULL){
    // Check if this onode exists in the destination
    if (onode_conflict_with_children(od,t) == 0){
      oc = onode_copy(t);
      onode_insert_child_at_end(od,oc);
      //printf("  NO CONFLICT FOR  %s\n",t->otype);
    }else
      ; //printf("  *** ---> CONFLICT FOR  %s\n",t->otype);
    t = t->next;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           PARAMFILE_ONODE_RESOLVE                         */
/*                                                                           */
/*****************************************************************************/
void paramfile_onode_resolve(o,mylogf)
     struct onode *o;
     char *mylogf;
{
  int i,k;
  int nid,ncop,tn,todo,nfail;
  char *cid,**id,*pname,tstr[SLEN];
  struct onode **olist,**cop_list;

  mylog(mylogf,"  PARAMFILE_ONODE_RESOLVE\n");

  // Get a list of pointers to all onodes with a tag value for 'ID'
  olist = onode_tree_get_list(o,4,"ID",&nid);

  id = (char **)myalloc(nid*sizeof(char *)); // Make a list of ID values
  for(i=0;i<nid;i++){
    id[i] = onode_get_tag_val_ptr(olist[i],"ID");
    sprintf(tstr,"    ID %2d  %s = %s\n",i,olist[i]->otype,id[i]);
    mylog(mylogf,tstr);
  }

  // Get a list of pointers to all onodes with a tag value for 'COPY'
  cop_list = onode_tree_get_list(o,4,"COPY",&ncop);

  todo = ncop;  // How many nodes left to resolve
  k = 0;        // index in copy list
  nfail = 0;    // Number of consecutive nodes that could not be resolved
  while(todo > 0){
    cid = onode_get_tag_val_ptr(cop_list[k],"COPY"); // ID of node to copy
    i = search_2d_carray(id,cid,nid);                // Find this ID in list

    tn = onode_tree_count(cop_list[k],5,NULL); // Check for issues to resolve
    if (tn == 1){ // Only the top node counted
      onode_copy_contents_novel(cop_list[k],olist[i]);

      pname = onode_getpar_chr_ptr(cop_list[k],"name");
      if (pname != NULL)
	sprintf(tstr,"    Resolving 'COPY' in <%s> %s\n",olist[i]->otype,
		pname);
      else
	sprintf(tstr,"    Resolving 'COPY' in <%s>\n",olist[i]->otype);
      mylog(mylogf,tstr);

      nfail = 0;
      todo -= 1;   // One less to do
      cop_list[k]->resolve = 0;
    }else if (tn > 1){
      nfail += 1;  // Count the number of failures
      if (nfail >= ncop)
	exit_error("PARAMFILE_ONODE_RESOLVE","Nested COPY impossible");
    }

    
    k += 1;      // Move on to the next onode, wrapping around if needed
    if (k >= ncop)
      k = 0;
  }

  if (olist != NULL)
    myfree(olist);
  if (id != NULL)
    myfree(id);
  if (cop_list != NULL)
    myfree(cop_list);
}
/**************************************-**************************************/
/*                                                                           */
/*                       MYXML_GET_NODE_TAG_WITH_CHILD_VAL                   */
/*                                                                           */
/*  Return the first node having the 'tagname' and having a child of type    */
/*  'childtag' with value 'childval'.                                        */
/*                                                                           */
/*****************************************************************************/
struct onode *myxml_get_node_tag_with_child_val(o,tagname,childtag,childval)
     struct onode *o;
     char *tagname;    // Find node that has type 'otype'
     char *childtag;   // and has a child with this 'otype'
     char *childval;   // and that child has this val
{
  struct onode *t,*oc,*of;

  of = NULL;
  t = onode_get_next_type(o->o,tagname);
  while(t != NULL){
    //printf("  t->otype = %s\n",t->otype);
    oc = onode_child_get_unique(t,childtag);
    if (oc != NULL){                      // has child of correct type
      //printf("    oc->val = %s\n",oc->val);
      if (strcmp(oc->val,childval)==0){      // child value is correct
	if (of == NULL){
	  of = t;
	}else{
	  exit_error("MYXML_GET_NODE_TAG_WITH_CHILD_VAL","Multiple matches");
	}
      }
    }
    t = onode_get_next_type(t->next,tagname);
  }

  return of;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MYXML_GET_TAG_VAL                            */
/*                                                                           */
/*  Return the POINTER to the value of the child having 'tagname'            */
/*                                                                           */
/*****************************************************************************/
char *myxml_get_tag_val(o,tagname)
     struct onode *o;
     char *tagname;    // Find node that has type 'otype'
{
  struct onode *oc;
  char *t;

  t = NULL;

  oc = onode_child_get_unique(o,tagname);
  if (oc != NULL){
    t = oc->val;
  }

  return t;
}
