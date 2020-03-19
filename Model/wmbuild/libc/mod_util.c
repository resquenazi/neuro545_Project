/*****************************************************************************/
/*                                                                           */
/*  mod_util.c                                                               */
/*  wyeth bair                                                               */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h> // For getlogin
#include <unistd.h>

#include "my_util.h"
#include "nr_util.h"
#include "carray_util.h"
#include "iarray_util.h"
#include "farray_util.h"
#include "func_util.h"
#include "data_util.h"
#include "paramfile_util.h"
#include "fft_util.h"
#include "kernel_util.h"
#include "stim_util.h"
#include "paramfile.h"
#include "mod.h"
#include "ndata.h"

/**************************************-**************************************/
/*                                                                           */
/*                           MYLOG_GET_GLOBAL_FNAME                          */
/*                                                                           */
/*****************************************************************************/
char *mylog_get_global_fname(id)
     int id;
{
  char tstr[SLEN];

  sprintf(tstr,"/tmp/zz.%d.mm.P_%d.log",(int)getuid(),id);

  // WYETH - WHY DOESN"T THIS WORK???
  //sprintf(tstr,"/tmp/zz.%s.mmm.P_%d.log",getlogin(),id);

  return strdup(tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                          MYLOG_GET_GLOBAL_JA_FNAME                        */
/*                                                                           */
/*****************************************************************************/
char *mylog_get_global_ja_fname()
{
  char tstr[SLEN];

  sprintf(tstr,"/tmp/zzz.%d.mm.ja",(int)getuid());

  return strdup(tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                               MYLOG_SET_ID_LOG                            */
/*                                                                           */
/*  Used by models called by 'mm' and 'wm'.                                  */
/*                                                                           */
/*****************************************************************************/
int mylog_set_id_log(myid,rmylogf)
     int myid;
     char **rmylogf;
{
  char *mylogf;

  // if (myid == -1){  // *** WYETH - old way - master gets a log file
  if ((myid == -1) || (myid == 0)){  // Give master a NULL log file
    mylogf = (char *)NULL;
  }else{
    mylogf = mylog_get_global_fname(myid);
  }

  if (*rmylogf == NULL){
    *rmylogf = mylogf; // Use the new name
  }else{
    if (strcmp(*rmylogf,mylogf)==0){
      myfree(mylogf);    // They are the same, keep the old name
    }else{
      myfree(*rmylogf);
      *rmylogf = mylogf; // Use the new name
    }
  }

  return myid;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_UTIL_SAMPLING_EQUAL                         */
/*                                                                           */
/*  Return 1 if the two sampling values are within ~1/1000th of each other.  */
/*                                                                           */
/*****************************************************************************/
int mod_util_sampling_equal(s1,s2)
     float s1,s2;
{
  int eqflag;

  if (my_rint(1000.0*s1) == my_rint(1000.0*s2))
    eqflag = 1;  // Close enough to equal
  else
    eqflag = 0;

  return eqflag;
}
/**************************************-**************************************/
/*                                                                           */
/*                        MOD_UTIL_SUBSTITUTE_VAR_PARAMS                     */
/*                                                                           */
/*  Fill in values in 's->ppl' for all var params for the 'k'th trial.       */
/*  If 'k' is <0, do nothing.                                                */
/*                                                                           */
/*****************************************************************************/
void mod_util_substitute_var_params(s,k,mylogf,nv,vname,vval)
     struct stim_struct *s;     // Stimulus params
     int k;                     // kth trial
     char mylogf[];
     int nv;         // Number of variables
     char **vname;   // Names of variables [nv]
     char **vval;    // Values of variables [nv]
{
  int i;
  int numflag;
  float new_val_f;
  char tstr[SLEN],*errstr,*new_val;
  struct param_pair *tp;

  /*** WYETH - WHY NOT USE CODE LIKE THIS, copied from 'lab_runtag_check':
       uflag = paramfile_list_update(ppl,pname,lm->setval[i]);
       if (uflag == 0)
       exit_error("LAB_RUNTAG_CHECK","This should not happen");
  ***/


  // WYETH 2019 - This comment is irrelevant now that 'mod_mon_util'
  //   has been eliminated.
  /*** WYETH - should I set 's->stimno' here so 'mod_mon_util' knows? ***/


  if (k >= 0){
    for(i=0;i<s->nvar;i++){ // For each var param

      // Check that the var param exists in the 'PPL' and free the old value
      tp = paramfile_get_param_pointer(s->ppl,s->vname[i]);
      if (tp == NULL){
	sprintf(tstr,"  %s\n",s->vname[i]);
	mylog(mylogf,tstr);
	mylog_exit(mylogf,"MOD_UTIL_SUBSTITUTE_VAR_PARAMS  Null points\n");
      }
      myfree(tp->value);
      /*** WARNING: the old value will include all on the line that follows the
	   name in the original file.  This includes comments.  Such
	   information will be replaced by a single value here. ***/


      // Handle expressions with variable names
      new_val = s->vval[k][i];

      // A numeric var param w/ non-numeric value must be an expression
      if ((s->vtype[i] != 'c') && (!is_float_string(new_val))){
	new_val_f = mstr_evaluate(new_val,nv,vname,vval,&errstr);
	if (errstr != NULL){
	  mylog(mylogf,errstr);
	  mylog(mylogf,"\n");
	  mylog_exit(mylogf,"MOD_UTIL_SUBSTITUTE_VAR_PARAMS mstr_evaluate\n");
	}
	sprintf(tstr,"%f",new_val_f);
	tp->value = strdup(tstr);

	sprintf(tstr,"    Evaluating expression, %s = %s\n",s->vname[i],
		tp->value);
	mylog(mylogf,tstr);
      }else{
	tp->value = strdup(new_val);
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_UTIL_SET_VAR_TN                           */
/*                                                                           */
/*  When 'tn' varies, do some checking.                                      */
/*                                                                           */
/*****************************************************************************/
int mod_util_set_var_tn(mylogf,s,tn_max,tscale)
     char mylogf[];             // Log file
     struct stim_struct *s;     // Stimulus params
     int tn_max;                // Maximum allowed 'tn'
     float tscale;              // Model tscale, (sec/sample)
{
  int vtn;
  float stim_tn;  // Nov 2009, Bug fix: this had been declared 'int'
  char ggstr[SLEN];

  stim_tn = param_getf_dflt(s->ppl,"stim_tn",0.0);
  if (stim_tn == 0.0){
    vtn = tn_max;    // Does not vary, use global max
  }else{
    vtn = (int)(stim_tn / tscale);
    sprintf(ggstr,"  Variable duration, %f sec, tn %d\n",stim_tn,vtn);
    mylog(mylogf,ggstr);

    if (vtn > tn_max){
      sprintf(ggstr,"  *** stim_tn requires %d steps, but tn is %d.\n",
	      vtn,tn_max);
      mylog(mylogf,ggstr);
      mylog_exit(mylogf,"MOD_UTIL_SET_VAR_TN  tn max exceeded.");
    }
  }
  return vtn;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_UTIL_VAR_PAR_REC_FREE                         */
/*                                                                           */
/*****************************************************************************/
void mod_util_var_par_rec_free(vpr)
     struct var_par_rec *vpr;
{
  free_2d_carray(vpr->name,vpr->n);
  myfree(vpr->ptype);
  free_2d_pointer_carray(vpr->val,vpr->n,vpr->nval);
  myfree(vpr);
}
/**************************************-**************************************/
/*                                                                           */
/*                        MOD_UTIL_VAR_PAR_LOOKUP_READ                       */
/*                                                                           */
/*  Read all <varpar_lookup> objects.                                        */
/*                                                                           */
/*****************************************************************************/
void mod_util_var_par_lookup_read(v,vpo)
     struct var_par_struct *v;   // Variable parameter configuration
     struct onode *vpo;          // Main <var_param> onode
{
  int j,k;
  int n,nv,nval;
  struct onode *t,*tco;

  n = onode_count_otype(vpo,"varpar_lookup");
  if (n == 0)
    return;  // There are no such objects.

  v->vplk_n = n;
  v->vplk_name = (char **)myalloc(n*sizeof(char *));
  v->vplk_nval = (int *)myalloc(n*sizeof(int));
  v->vplk_key  = (char ***)myalloc(n*sizeof(char **));
  v->vplk_val  = (char ***)myalloc(n*sizeof(char **));

  t = onode_get_next_type(vpo->o,"varpar_lookup");
  k = 0;
  while(t != NULL){
    nv = t->n - 1;  // Count all items except "name"
    v->vplk_name[k] = onode_getpar_chr_exit(t,"name");
    v->vplk_nval[k] = nv;
    v->vplk_key[k] = (char **)myalloc(nv*sizeof(char *));
    v->vplk_val[k] = (char **)myalloc(nv*sizeof(char *));

    //  Now, scan the children within this object:
    //
    //  <varpar_lookup>
    //    name w
    //    1   0 0 0
    //    2   0 1 2
    //    3   2 1 0
    //  </varpar_lookup>
    //
    j = 0;
    tco = t->o; // head node in child list
    nval = -1;
    while(tco != NULL){
      if (strcmp(tco->otype,"item")!=0)
	exit_error("MOD_UTIL_VAR_PAR_LOOKUP_READ","Expecting 'item' node");
      if (strcmp(tco->name,"name")!=0){
	v->vplk_key[k][j] = strdup(tco->name);
	v->vplk_val[k][j] = strdup(tco->val);
	if (nval == -1){
	  nval = tco->nval;
	}else{
	  if (tco->nval != nval)
	    exit_error("MOD_UTIL_VAR_PAR_LOOKUP_READ",
		       "Number of values varies in <varpar_lookup>");
	}

	//printf("(mod_util) WYETH HERE   key %s  val %s  (nval %d)\n",
	//v->vplk_key[k][j],v->vplk_val[k][j],tco->nval);

	j += 1;
      }
      tco = tco->next; // head node in child list
    }

    t = onode_get_next_type(t->next,"varpar_lookup");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                       MOD_UTIL_VAR_PAR_LOOKUP_REPLACE                     */
/*                                                                           */
/*  Replace the values of a var par with its lookup values.                  */
/*                                                                           */
/*****************************************************************************/
void mod_util_var_par_lookup_replace(name,vals,n,v)
     char *name;   // Parameter name to look up
     char **vals;  // [n] parameter values, to be replaced
     int n;        // Number of parameter values
     struct var_par_struct *v;   // Variable parameter configuration
{
  int i,k,vi;
  int nkey;
  char **vkey,**vval;

  k = search_2d_carray(v->vplk_name,name,v->vplk_n);
  if (k < 0)
    exit_error("MOD_UTIL_VAR_PAR_LOOKUP_REPLACE","Name not found");

  nkey = v->vplk_nval[k];  // Number of keys for this param name
  vkey = v->vplk_key[k];   // [nkey] List of key values
  vval = v->vplk_val[k];   // [nkey] List of values associated with key

  //printf("-----------REPLACING Values for:  %s\n",name);
  //printf("  there are %d values.  There are %d keys\n",n,nkey);

  for(i=0;i<n;i++){
    //printf("  oldval[%d]:  %s\n",i,vals[i]);

    vi = search_2d_carray(vkey,vals[i],nkey);
    if (vi < 0){
      printf("  *** Param name:  %s   Invalid key:  %s\n",name,vals[i]);
      exit_error("MOD_UTIL_VAR_PAR_LOOKUP_REPLACE","Key not found");
    }

    myfree(vals[i]);
    vals[i] = strdup(vval[vi]);
    //printf("     new[%d]:  %s\n",i,vals[i]);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                       MOD_UTIL_VAR_PAR_LOOKUP_REVERSE                     */
/*                                                                           */
/*  Reverse Lookup:                                                          */
/*  Return a pointer to the key value associated with 'val' if 'pname' is    */
/*  one of the varpar_lookup parameters, or NULL if it is not.               */
/*                                                                           */
/*****************************************************************************/
char *mod_util_var_par_lookup_reverse(pname,val,v)
     char *pname;   // Parameter name to look up
     char *val;     // Value of the parameter, for which key is sought
     struct var_par_struct *v;   // Variable parameter configuration
{
  int k,vi;

  k = search_2d_carray(v->vplk_name,pname,v->vplk_n);
  if (k < 0)
    return NULL;  // 'pname' is not associated with a <varpar_lookup>

  vi = search_2d_carray(v->vplk_val[k],val,v->vplk_nval[k]);
  if (vi < 0)
    exit_error("MOD_UTIL_VAR_PAR_LOOKUP_REVERSE","Unknown value for lookup");

  return v->vplk_key[k][vi];
}
/**************************************-**************************************/
/*                                                                           */
/*                        MOD_UTIL_VAR_PAR_REC_GET_ONE                       */
/*                                                                           */
/*****************************************************************************/
struct var_par_rec *mod_util_var_par_rec_get_one(name,n,ptype,val,cpflag)
     char *name;   // Param name
     int n;        // Number of values
     char ptype;   // 'f', 'i', 'c'
     char **val;   // [n] Value list
     int cpflag;   // 0-no copy, 1-copy val, 2-copy name and val
{
  int i;
  struct var_par_rec *v;

  exit_error("(mod_util)  WYETH HERE","Unused, Untested");

  v = (struct var_par_rec *)myalloc(sizeof(struct var_par_rec));

  v->n = 1;
  v->name = (char **)myalloc(sizeof(char *));
  if (cpflag < 2)
    v->name[0] = strdup(name);
  else
    v->name[0] = name;

  v->ptype[0] = ptype;
  v->nval = n;

  v->val = (char ***)myalloc(sizeof(char **));
  if (cpflag > 0){
    v->val[0] = val;
  }else{
    v->val[0] = (char **)myalloc(n*sizeof(char *));
    for(i=0;i<n;i++)
      v->val[0][i] = strdup(val[i]);
  }

  return v;
}  
/**************************************-**************************************/
/*                                                                           */
/*                      MOD_UTIL_VAR_PAR_REC_FROM_ONODE                      */
/*                                                                           */
/*****************************************************************************/
struct var_par_rec *mod_util_var_par_rec_from_onode(o,v)
     struct onode *o;
     struct var_par_struct *v;
{
  int i,j;
  int nval,iflag;
  char *sval,*s0,*s1,tstr[SLEN];
  float v0,v1,vf;
  struct onode *t;
  struct var_par_rec *vpr;

  vpr = (struct var_par_rec *)myalloc(sizeof(struct var_par_rec));

  vpr->n = o->n;  // Number of linked var param
  vpr->val = (char ***)myalloc(vpr->n*sizeof(char **));
  vpr->name = (char **)myalloc(vpr->n*sizeof(char *));
  vpr->ptype = (char *)myalloc(vpr->n*sizeof(char));

  t = o->o;  // First child
  for(i=0;i<vpr->n;i++){

    //printf("NAME = %s\n",t->name);
    //printf("  type = %s\n",t->otype);

    vpr->name[i]  = strdup(t->name);
    vpr->ptype[i] = 'c';               // Assuming char for now

    //
    //  Determine the number of values associated with this line
    //
    if ((strcmp(o->otype,"var_moo")==0)||(strcmp(o->otype,"var_list")==0)){
      //
      //  A simple list of values:  [name]  <v0> ... <vn-1>
      //
      nval = t->nval;
    }else if (strcmp(o->otype,"var_range_n")==0){
      //
      //  A range and count:  [name]  <v0> <v1> <n>
      //
      sval = onode_get_nth_val(t,2);
      nval = atoi(sval);
      myfree(sval);
    }

    //
    //  Do not allow variation in 'nval' across linked params
    //
    if (i==0)
      vpr->nval = nval;
    else{
      if (vpr->nval != nval)
	exit_error("MOD_UTIL_VAR_PAR_REC_FROM_ONODE","nval != nvar");
    }

    //
    //  Extract, or create, the list of values for this param
    //
    vpr->val[i] = (char **)myalloc(nval*sizeof(char *));
    if ((strcmp(o->otype,"var_moo")==0)||(strcmp(o->otype,"var_list")==0)){
      for(j=0;j<nval;j++)
	vpr->val[i][j] = onode_get_nth_val(t,j);

      if (vpr->name[i][0] == '$'){
	//
	//  The values for this variable need to be replaced with those
	//  from the <varpar_lookup> record.
	//

	//  Remove the '$' from the name
	myfree(vpr->name[i]);
	vpr->name[i]  = strdup(&(t->name[1]));

	mod_util_var_par_lookup_replace(vpr->name[i],vpr->val[i],nval,v);

	//exit_error("WYETH HERE","Under development");
      }

    }else if (strcmp(o->otype,"var_range_n")==0){

      s0 = onode_get_nth_val(t,0);
      s1 = onode_get_nth_val(t,1);
      if (is_int_string(s0) && is_int_string(s0))
	iflag = 1;
      else
	iflag = 0;
      v0 = atof(s0);
      v1 = atof(s1);

      for(j=0;j<nval;j++){
	vf = v0 + (v1-v0)*(float)j/(float)(nval-1);
	if (iflag == 1)
	  sprintf(tstr,"%d",my_rint(vf));
	else
	  sprintf(tstr,"%f",vf);
	vpr->val[i][j] = strdup(tstr);
      }
      myfree(s0);
      myfree(s1);
    }

    t = t->next;  // Next child
  }

  return vpr;
}
/**************************************-**************************************/
/*                                                                           */
/*                       MOD_UTIL_VAR_PAR_REC_CHECK_MOO                      */
/*                                                                           */
/*  Return 1 if the param names in the record are all valid, else exit.      */
/*                                                                           */
/*****************************************************************************/
int mod_util_var_par_rec_check_moo(mylogf,m,vpr)
     char *mylogf;             // Log file
     struct model_struct *m;   // Model params
     struct var_par_rec *vpr;  // Var par record
{
  int i;
  char ggstr[SLEN];

  for(i=0;i<vpr->n;i++){
    if (onode_test_ostr(m->o,vpr->name[i]) == 0){
      sprintf(ggstr,
	      "*** MOD_UTIL_VAR_PAR_REC_CHECK_MOO   Var name not found:  %s",
	      vpr->name[i]);
      mylog_exit(mylogf,ggstr);
    }
  }

  return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_UTIL_VAR_PAR_FREE                           */
/*                                                                           */
/*****************************************************************************/
void mod_util_var_par_free(vps)
     struct var_par_struct *vps;
{
  int i;

  for(i=0;i<vps->m_var_n;i++)
    mod_util_var_par_rec_free(vps->m_var[i]);

  for(i=0;i<vps->m_one_n;i++)
    mod_util_var_par_rec_free(vps->m_one[i]);

  // The contents of these name pointers was freed by '...rec_free' above
  //   thus we do not use 'free_2d_carray' here
  myfree(vps->m_name);

  free_2d_pointer_carray(vps->m_vval,vps->m_n,vps->m_nvar);

  myfree(vps);
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_UTIL_VAR_PAR_GET_INIT                         */
/*                                                                           */
/*****************************************************************************/
struct var_par_struct *mod_util_var_par_get_init()
{
  struct var_par_struct *v;

  v = (struct var_par_struct *)myalloc(sizeof(struct var_par_struct));
  v->m_var_n = 0;
  v->m_one_n = 0;
  v->m_var = NULL;
  v->m_one = NULL;

  v->vplk_n = 0;
  v->vplk_name = NULL;
  v->vplk_nval = NULL;
  v->vplk_key  = NULL;
  v->vplk_val  = NULL;

  return v;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_UTIL_VAR_PAR_CROSS_TABLE                      */
/*                                                                           */
/*  Make the cross-product table for the given moo var param configuration.  */
/*                                                                           */
/*  *** NOTE SEE ALSO:  (stim_util.c) GET_CROSS_PRODUCT_PARAMS               */
/*  *** THESE TWO ROUTINES CAN / SHOULD be merged.                           */
/*        (would have to add ability to handle singletons here)              */
/*  *** Ultimate - put this as a 'carray_util' ?                             */
/*                                                                           */
/*****************************************************************************/
void mod_util_var_par_cross_table(mylogf,ndim,vlist,rn,rlist)
     char *mylogf;
     int ndim;                    // Number of variable dimensions
     struct var_par_rec **vlist;  // [ndim] params and values for each dim.
     int *rn;                     // total number of model configs
     char ****rlist;              // [*rn][npar] var param values
{
  int i,j,k,l;
  int **list,nlist,*cnt,npar,pi;
  char ***vv;

  cnt = (int *)myalloc(ndim*sizeof(int)); // # of vals along each dimension
  npar = 0;
  for(i=0;i<ndim;i++){
    cnt[i] = vlist[i]->nval;
    npar  += vlist[i]->n;
  }

  get_nd_cross_product_list_iarray(ndim,cnt,&list,&nlist);

  vv = get_2d_pointer_carray(nlist,npar);

  for(i=0;i<nlist;i++){

    pi = 0;  // Param index (over all var pars in all dimensions)
    for(j=0;j<ndim;j++){  // For each dimension
      for(k=0;k<vlist[j]->n;k++){  // For each linked par in this dimension
	l = list[i][j];
	vv[i][pi] = strdup(vlist[j]->val[k][l]);
	pi += 1;
      }
    }
  }

  free_2d_iarray(list,nlist);
  myfree(cnt);

  *rn = nlist;
  *rlist = vv;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_UTIL_VAR_PAR                             */
/*                                                                           */
/*****************************************************************************/
struct var_par_struct *mod_util_var_par(mylogf,m,vpo)
     char *mylogf;            // Log file
     struct model_struct *m;  // Model params, or NULL
     struct onode *vpo;       // Main <var_param> onode
{
  int i,j,k;
  int tn;
  char ggstr[SLEN],***vval;
  struct var_par_struct *v;
  struct onode *t;

  mylog(mylogf,"  MOD_UTIL_VAR_PAR\n");

  v = mod_util_var_par_get_init();

  //
  //  Read data for all "varpar_lookup" objects (returns if none).
  //
  mod_util_var_par_lookup_read(v,vpo);


  //
  //  Count <var_moo> and other similar objects
  //
  v->m_var_n = onode_count_otype(vpo,"var_moo");
  sprintf(ggstr,"    %d <var_moo> records found\n",v->m_var_n);
  mylog(mylogf,ggstr);

  tn = onode_count_otype(vpo,"var_list");  // Alias for "var_moo"
  sprintf(ggstr,"    %d <var_list> records found\n",tn);
  mylog(mylogf,ggstr);
  v->m_var_n += tn;

  tn = onode_count_otype(vpo,"var_range_n");
  sprintf(ggstr,"    %d <var_range_n> records found\n",tn);
  mylog(mylogf,ggstr);
  v->m_var_n += tn;
  
  //
  //  Create/store <var_moo> records
  //
  v->m_var = (struct var_par_rec **)myalloc(v->m_var_n*
					    sizeof(struct var_par_rec *));
  k = 0;
  t = onode_get_next_type(vpo->o,"var_moo");
  while(t != NULL){
    v->m_var[k] = mod_util_var_par_rec_from_onode(t,v);
    if (m != NULL) // If any of the names do not exist in moo file, exit
      (void)mod_util_var_par_rec_check_moo(mylogf,m,v->m_var[k]);
    k += 1;
    t = onode_get_next_type(t->next,"var_moo");
  }

  t = onode_get_next_type(vpo->o,"var_list");
  while(t != NULL){
    v->m_var[k] = mod_util_var_par_rec_from_onode(t,v);
    if (m != NULL) // If any of the names do not exist in moo file, exit
      (void)mod_util_var_par_rec_check_moo(mylogf,m,v->m_var[k]);
    k += 1;
    t = onode_get_next_type(t->next,"var_list");
  }

  t = onode_get_next_type(vpo->o,"var_range_n");
  while(t != NULL){
    v->m_var[k] = mod_util_var_par_rec_from_onode(t,v);
    if (m != NULL) // If any of the names do not exist in moo file, exit
      (void)mod_util_var_par_rec_check_moo(mylogf,m,v->m_var[k]);
    k += 1;
    t = onode_get_next_type(t->next,"var_range_n");
  }

  //
  //  Build the moo var param table
  //
  v->m_n = 1;      // Number of total model configs
  v->m_nvar = 0;   // Number of total var params
  for(i=0;i<v->m_var_n;i++){
    v->m_n *= v->m_var[i]->nval;
    v->m_nvar += v->m_var[i]->n;
  }
  sprintf(ggstr,"    %d total model configurations\n",v->m_n);
  mylog(mylogf,ggstr);
  sprintf(ggstr,"    %d total model variable parameters\n",v->m_nvar);
  mylog(mylogf,ggstr);

  //
  //  Create/fill the list of all moo var param names
  //
  v->m_name = (char **)myalloc(v->m_nvar*sizeof(char *));
  k = 0;
  for(i=0;i<v->m_var_n;i++){
    for(j=0;j<v->m_var[i]->n;j++){
      v->m_name[k] = v->m_var[i]->name[j];  // Set pointer to name string
      k += 1;
    }
  }
  for(i=0;i<v->m_nvar;i++){
    sprintf(ggstr,"      %d %s\n",i,v->m_name[i]);
    mylog(mylogf,ggstr);
  }

  mod_util_var_par_cross_table(mylogf,v->m_var_n,v->m_var,&tn,&vval);
  if (tn != v->m_n)
    exit_error("MOD_UTIL","This should not happen");

  // Print the table for debugging
  /***
  for(i=0;i<v->m_n;i++){
    for(j=0;j<v->m_nvar;j++){
      printf(" %s",vval[i][j]);
    }
    printf("\n");
    }***/

  v->m_vval = vval;

  return v;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MODU_MDS_WRITE                              */
/*                                                                           */
/*  Write out a list of the model data structure entries.                    */
/*                                                                           */
/*****************************************************************************/
void modu_mds_write(m,outfile,vflag)
     struct model_struct *m;    // Model
     char *outfile;             // filename or NULL to print
     int vflag;                 // Version format flag
{
  FILE *fopen(),*fout;
  struct model_data_struct *md;

  printf("  MODU_MDS_WRITE\n");

  if (outfile != NULL){
    if ((fout = fopen(outfile,"w")) == NULL){
      printf("  *** File %s\n",outfile);
      exit_error("MODU_MDS_WRITE","Cannot open file");
    }
  }

  if (vflag == 0){
    md = m->mds;
    while(md != NULL){
      if (outfile != NULL){
	fprintf(fout,"    %s %s\n",md->pname,md->dataid);
      }else{
	printf("    %s %s\n",md->pname,md->dataid);
      }
      md = md->next;
    }
  }else{
    exit_error("MODU_MDS_WRITE","Unknown 'vflag'");
  }

  if (outfile != NULL)
    fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                              MODU_MDS_CREATE                              */
/*                                                                           */
/*  Write out a list of the model data structure entries.                    */
/*                                                                           */
/*****************************************************************************/
struct model_data_struct *modu_mds_create(pname,dataid,dtype,ndim,dlist,
					  dptr,davail,stimdep,ospec)
     char *pname;   // population name
     char *dataid;  // Data ID
     char  dtype;   // 'f' (float), 's' (spike)
     int   ndim;    // dimension of population, e.g., 3
     int  *dlist;   // [dim] length of each dimension
     void *dptr;    // Pointer to multi-dimensional data, or NULL
     char *davail;  // 'full', 'part_nullsig', 'part_unknow', 'other'
     char *stimdep; // 'no', 'binoc' (do not use NULL)
     char *ospec;   // Other specification, or NULL
{
  struct model_data_struct *md;

  /*** WYETH DEBUG
  printf("pname = %s\n",pname);
  printf("  dataid = %s\n",dataid);
  printf("  dtype = %c\n",dtype);
  printf("  ndim = %d\n",ndim);
  {
    int i;
    for(i=0;i<ndim;i++)
      printf("    %d",dlist[i]);
    printf("\n");
  }
  printf("  davail = %s\n",davail);
  printf("  stimdep = %s\n",stimdep);
  printf("  ospec = %s\n",ospec);
  ***/

  md = (struct model_data_struct *)myalloc(sizeof(struct model_data_struct));

  md->pname   = strdup(pname);
  md->dataid  = strdup(dataid);
  md->dtype   = dtype;
  md->ndim    = ndim;
  md->dlist   = copy_iarray(dlist,ndim);
  md->dptr    = dptr;
  md->davail  = strdup(davail);
  md->stimdep = strdup(stimdep);
  if (ospec != NULL)   md->ospec  = strdup(ospec);
  else                 md->ospec  = NULL;

  md->prev = NULL;
  md->next = NULL;

  return md;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MODU_MDS_GET_PTR_N                            */
/*                                                                           */
/*  Return a pointer to the nth element in the list.                         */
/*                                                                           */
/*****************************************************************************/
struct model_data_struct *modu_mds_get_ptr_n(m,n)
     struct model_struct *m;    // Model
     int n;
{
  int i;
  struct model_data_struct *t;

  if ((n < 0) || (n >= m->mds_n)){
    printf("  n = %d\n",n);
    printf("  m->mds_n = %d\n",m->mds_n);
    exit_error("MODU_MDS_GET_PTR_N","Index out of range");
  }

  t = m->mds;

  for(i=0;i<n;i++){
    t = t->next;
  }

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MODU_MDS_INSERT                              */
/*                                                                           */
/*  Add a MDS record to the list.                                            */
/*                                                                           */
/*****************************************************************************/
void modu_mds_insert(m,mds,posi)
     struct model_struct      *m;    // Model
     struct model_data_struct *mds;  // record to insert
     int                       posi; // Position index, 0=head, -1=tail
{
  struct model_data_struct *t;

  if (m->mds == NULL){

    if ((posi < -1) || (posi > 0))
      exit_error("MODU_MDS_INSERT","Index position is invalid");
    m->mds = mds;
    m->mds->prev = m->mds->next = NULL;
  }else if (posi == -1){
    t = modu_mds_get_ptr_n(m,m->mds_n-1);  // Pointer to last element

    mds->next = NULL;
    mds->prev = t;
    t->next = mds;
  }else{

    t = modu_mds_get_ptr_n(m,posi);
    mds->next = t;
    mds->prev = t->prev;

    if (t->prev != NULL)
      t->prev->next = mds;
    else
      m->mds = mds;  // New head of list
    t->prev = mds;
  }

  m->mds_n += 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MODU_MDS_ADD_3D                              */
/*                                                                           */
/*  Add a MDS record to the list.                                            */
/*                                                                           */
/*****************************************************************************/
void modu_mds_add_3d(m,pname,dataid,dtype,xn,yn,zn,dptr,davail,stimdep,ospec,
		     posi)
     struct model_struct *m;    // Model
     char *pname;   // population name
     char *dataid;  // Data ID
     char  dtype;   // 'f' (float), 's' (spike)
     int   xn,yn,zn;  // dimensions of 'dptr'
     void *dptr;    // Pointer to multi-dimensional data, or NULL
     char *davail;  // 'full', 'part_nullsig', 'part_unknow', 'other'
     char *stimdep; // 'no', 'binoc' (do not use NULL)
     char *ospec;   // Other specification, or NULL
     int   posi;    // Position index, 0=head, -1=tail
{
  int *dlist;
  struct model_data_struct *mds;

  //printf("stimdep = %s\n",stimdep);
  //printf("ospec = %s\n",ospec);
  //printf("posi = %d\n",posi);

  dlist = (int *)myalloc(3*sizeof(int));
  dlist[0] = xn;
  dlist[1] = yn;
  dlist[2] = zn;

  mds = modu_mds_create(pname,dataid,dtype,3,dlist,dptr,davail,stimdep,ospec);

  modu_mds_insert(m,mds,posi);

  myfree(dlist);
}
/**************************************-**************************************/
/*                                                                           */
/*                        MOD_UTIL_RESP_UNFILLED_COUNT                       */
/*                                                                           */
/*  Count the number of unfilled responses.                                  */
/*                                                                           */
/*****************************************************************************/
int mod_util_resp_unfilled_count(r)
     struct response_struct *r;
{
  int i;
  int n,ri,ti;

  //ti = r->tsi;   // Trial index
  ti = r->gtsi;   // Trial index

  n = 0;
  for(i=0;i<r->n;i++){

    ri = r->ri[i];  // Response storage index

    if (strcmp(r->rtype[i],"s")==0){
      if (r->cnt[ri][ti] == -1)
	n += 1;
    }else if (strcmp(r->rtype[i],"f")==0){
      if (r->fcnt[ri][ti] == -1)
	n += 1;
    }else{
      exit_error("MOD_UTIL_RESP_UNFILLED_COUNT","Unknown response type.");
    }
  }
  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_UTIL_RESP_UNFILLED_EXIT                       */
/*                                                                           */
/*  If there are unfilled responses, report them and exit.                   */
/*                                                                           */
/*****************************************************************************/
void mod_util_resp_unfilled_exit(r,mylogf)
     struct response_struct *r;
     char *mylogf;
{
  int i;
  int n,ri,ti;

  //ti = r->tsi;   // Trial index
  ti = r->gtsi;   // Trial index

  n = 0;
  for(i=0;i<r->n;i++){

    ri = r->ri[i];  // Response storage index

    if (strcmp(r->rtype[i],"s")==0){
      if (r->cnt[ri][ti] == -1){
	printf("  UNFILLED RESPONSE  %s %s\n",r->nd_name[i],r->datid[i]);
	n += 1;
      }
    }else if (strcmp(r->rtype[i],"f")==0){
      if (r->fcnt[ri][ti] == -1){
	printf("  UNFILLED RESPONSE  %s %s\n",r->nd_name[i],r->datid[i]);
	n += 1;
      }
    }else{
      exit_error("MOD_UTIL_RESP_UNFILLED_EXIT","Unknown response type.");
    }
  }

  if (n > 0)
    mylog_exit(mylogf,"MOD_UTIL_RESP_UNFILLED_EXIT  Unfilled responses.");

}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_UTIL_RESP_GET_RTYPE_PTR                       */
/*                                                                           */
/*  Return a pointer to the 'rtype' string.                                  */
/*                                                                           */
/*****************************************************************************/
char *mod_util_resp_get_rtype_ptr(r,name)
     struct response_struct *r;
     char name[];
{
  int k,ri;
  char *rtype;

  k = search_2d_carray(r->plname,name,r->n);
  if (k >= 0)
    rtype = r->rtype[k];
  else
    rtype = (char *)NULL;

  return rtype;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_UTIL_RESP_GET_RI                           */
/*                                                                           */
/*  Return -1 if the 'name' is not found.                                    */
/*                                                                           */
/*****************************************************************************/
int mod_util_resp_get_ri(r,name)
     struct response_struct *r;
     char name[];
{
  int k,ri;

  k = search_2d_carray(r->plname,name,r->n);
  if (k >= 0)
    ri = r->ri[k];
  else
    ri = -1;

  return ri;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_UTIL_RESP_GET_SAMP                         */
/*                                                                           */
/*  Return sampling, or 0.0 if 'name' not found.                             */
/*                                                                           */
/*****************************************************************************/
float mod_util_resp_get_samp(r,name)
     struct response_struct *r;
     char name[];
{
  int k;
  float samp;

  k = search_2d_carray(r->plname,name,r->n);
  if (k >= 0){
    samp = r->samp[k];
  }else
    samp = 0.0;

  return samp;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_UTIL_RESP_SET_SAMP                         */
/*                                                                           */
/*  Overwrite the sampling value.                                            */
/*                                                                           */
/*****************************************************************************/
void mod_util_resp_set_samp(r,name,samp)
     struct response_struct *r;
     char name[];
     float samp;
{
  int k;

  k = search_2d_carray(r->plname,name,r->n);
  if (k >= 0)
    r->samp[k] = samp;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_UTIL_RESP_CHECK_NAME                         */
/*                                                                           */
/*  Return 1 if 'name' is found and matches 'rtype', 0 if not found.         */
/*  Exit if rtype conflict, logging to mylogf.                               */
/*                                                                           */
/*****************************************************************************/
int mod_util_resp_check_name(r,name,rtype,mylogf)
     struct response_struct *r;
     char name[];
     char rtype[];
     char mylogf[];
{
  int flag;
  char tstr[SLEN],*rp;

  rp = mod_util_resp_get_rtype_ptr(r,name);
  if (rp != NULL){
    flag = 1;
    if (strcmp(rp,rtype)!=0){
      sprintf(tstr,"*** rsp file:  %s   should be:  %s\n",rp,rtype);
      mylog(mylogf,tstr);
      sprintf(tstr,"MOD_UTIL_RESP_CHECK_NAME Wrong rtype for %s\n",name);
      mylog_exit(mylogf,tstr);
    }
  }else
    flag = 0;
  
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_UTIL_RESP_GET_CHECK_NAME                      */
/*                                                                           */
/*  Return 1 if 'name' is found and matches 'rtype', 0 if not found.         */
/*  Exit if rtype conflict, logging to mylogf.                               */
/*                                                                           */
/*****************************************************************************/
int mod_util_resp_get_check_name(r,name,rtype,mylogf,rri)
     struct response_struct *r;
     char name[];
     char rtype[];
     char mylogf[];
     int *rri;
{
  int ri,flag;
  char tstr[SLEN],*rp;
  
  ri = mod_util_resp_get_ri(r,name);
  if (ri >= 0){
    flag = 1;
    rp = mod_util_resp_get_rtype_ptr(r,name);
    if (strcmp(rp,rtype)!=0){
      sprintf(tstr,"*** rsp file:  %s   should be:  %s\n",rp,rtype);
      mylog(mylogf,tstr);
      sprintf(tstr,"MOD_UTIL_RESP_GET_CHECK_NAME Wrong rtype for %s\n",name);
      mylog_exit(mylogf,tstr);
    }
  }else
    flag = 0;

  *rri = ri;
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_UTIL_RESP_CHECK_STORE_F                       */
/*                                                                           */
/*  Return 1 if 'name' is in the response list.                              */
/*  Exit if rtype conflict.                                                  */
/*  Exit if storage pointer not NULL.                                        */
/*  Copy the data if 'cpflag' = 1.                                           */
/*                                                                           */
/*  WYETH - OLD WAY - phase this out                                         */
/*  WYETH - OLD WAY - phase this out                                         */
/*  WYETH - OLD WAY - phase this out                                         */
/*  WYETH - OLD WAY - phase this out                                         */
/*****************************************************************************/
int mod_util_resp_check_store_f(r,name,mylogf,data,n,cpflag)
     struct response_struct *r;
     char name[];
     char mylogf[];
     float *data;
     int n;
     int cpflag;
{
  int ri,ti,flag;
  char tstr[SLEN];

  exit_error("MOD_UTIL.C","THIS NOT USED - WYETH - REMOVE");

  flag = mod_util_resp_get_check_name(r,name,"f",mylogf,&ri);
  if (flag){
    //ti = r->tsi;  // Trial index
    ti = r->gtsi;  // Trial index
    if (r->f[ri][ti] != NULL){
      sprintf(tstr,"MOD_UTIL_RESP_CHECK_STORE_F Already stored %s\n",name);
      mylog_exit(mylogf,tstr);
    }
    if (cpflag == 1)
      r->f[ri][ti] = copy_farray(data,n);
    else
      r->f[ri][ti] = data;
    r->fcnt[ri][ti] = n;
  }
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_UTIL_RESP_CHECK_STORE_S                       */
/*                                                                           */
/*  Return 1 if 'name' is in the response list.                              */
/*  Exit if rtype conflict.                                                  */
/*  Exit if storage pointer not NULL.                                        */
/*  Copy the data if 'cpflag' = 1.                                           */
/*                                                                           */
/*  WYETH - OLD WAY - phase this out                                         */
/*  WYETH - OLD WAY - phase this out                                         */
/*  WYETH - OLD WAY - phase this out                                         */
/*  WYETH - OLD WAY - phase this out                                         */
/*                                                                           */
/*****************************************************************************/
int mod_util_resp_check_store_s(r,name,mylogf,s,cnt,cpflag)
     struct response_struct *r;
     char name[];
     char mylogf[];
     float *s;
     int cnt;
     int cpflag;
{
  int ri,ti,flag;
  char tstr[SLEN];

  flag = mod_util_resp_get_check_name(r,name,"s",mylogf,&ri);

  if (flag == 0){
    printf("  Looking for name:  %s\n",name);
    exit_error("MOD_UTIL_RESP_CHECK_STORE_S","flag is 0");
  }if (ri < 0)
    exit_error("MOD_UTIL_RESP_CHECK_STORE_S","ri < 0");

  if (flag){
    //ti = r->tsi;  // Trial index
    ti = r->gtsi;  // Trial index
    if (r->s[ri][ti] != NULL){
      sprintf(tstr,"MOD_UTIL_RESP_CHECK_STORE_S Already stored %s\n",name);
      mylog_exit(mylogf,tstr);
    }
    if (cpflag == 1)
      r->s[ri][ti] = copy_farray(s,cnt);
    else
      r->s[ri][ti] = s;
    r->cnt[ri][ti] = cnt;
  }

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                        MOD_UTIL_RESP_DERIVE_CURRENT                       */
/*                                                                           */
/*  Compute the current associated with the conductance, voltage, and        */
/*  V_rev.                                                                   */
/*                                                                           */
/*****************************************************************************/
float *mod_util_resp_derive_current(g,gsamp,gx,gn,tscale,v,vx,vn,vmax,vrev,
				    mylogf,rn)
     float *g;        // Conductance nS [gn]
     float gsamp;     // samples/s, if regular, 0.0 otherwise
     float *gx;       // times for conductance, NULL if regular
     int  gn;         // length of 'g'
     float tscale;    // Multiply time by this to get sec
     float *v;        // voltage mV [vn]
     float *vx;       // times for voltage (sec)
     int  vn;         // length of voltage trace
     float vmax;      // Set any 'v' higher than this to this value
     float  vrev;     // reversal potential for g
     char mylogf[];
     int *rn;         // Number of values returned
{
  int i,k;
  int done;
  float *di,*vt,vm,t;

  // Make a copy of 'v' and remove any spikes
  vt = copy_farray(v,vn);
  for(i=0;i<vn;i++)
    if (vt[i] > vmax)
      vt[i] = vmax;

  di = (float *)myalloc(gn*sizeof(float));

  k = 0;
  i = 0;
  done = 0;
  if (gsamp > 0.0){     // Regularly sampled g
    while(!done && (k<gn)){
      t = tscale * (float)i / gsamp;  // Time
      if (t <= (vx[vn-1] + 0.00001)){ // Sometimes t is sligtly greater
	vm = lin_interp_farray(vx,vt,vn,t);
	di[k] = g[i] * (vrev - vm) / 1000.0; // nA = nS * mV / 1000.0
	k += 1;
      }else{
	//printf("t = %f   vx[[vn-1] = %f\n",t,vx[vn-1]);
	done = 1;  // There are no more Vm values
      }
      i += 1;
    }
  }else{                // Sample times in 'gx'
    while(!done && (i<gn)){
      t = tscale * gx[i];  // Time
      if (t <= vx[vn-1]){
	vm = lin_interp_farray(vx,vt,vn,t);
	di[i] = g[i] * (vrev - vm) / 1000.0; // pA = nS * mV
	k += 1;
      }else{
	done = 1;
      }
      i += 1;
    }
  }
  myfree(vt);

  *rn = k;
  return di;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_UTIL_RESP_STORE_S                          */
/*                                                                           */
/*  Copy the data if 'cpflag' = 1.                                           */
/*                                                                           */
/*  NEW WAY (phase out 'mod_util_resp_check_store_s' )                       */
/*                                                                           */
/*****************************************************************************/
void mod_util_resp_store_s(r,k,s,cnt,cpflag,mylogf)
     struct response_struct *r;
     int k;          // Response index
     float *s;
     int cnt;
     int cpflag;
     char mylogf[];
{
  int ri,ti;
  char tstr[SLEN];

  if (strcmp(r->rtype[k],"s")!=0)  // WYETH Added Sep 8, 2015
    exit_error("MOD_UTIL_RESP_STORE_S",
	       "In .rsp file, 'spikes' must have type 's'");

  //ti = r->tsi;    // Trial index
  ti = r->gtsi;    // Trial index
  ri = r->ri[k];  // pointer within "s" response list
  //printf("  _store_s   Trial index:  %d   k %d  ri %d\n",ti,k,ri);
  if (r->s[ri][ti] != NULL){
    //printf("  Trial index:  %d\n",ti);
    sprintf(tstr,"*** WARNING:  MOD_UTIL_RESP_STORE_S Already stored %s\n",
	    r->nd_name[k]);
    mylog(mylogf,tstr);
    // Allow overwrite in case of mm following a suspension.
    sprintf(tstr,"    Trial %d results are being over-written\n",ti);
    mylog(mylogf,tstr);
    myfree(r->s[ri][ti]);
  }

  if (cpflag == 1){
    r->s[ri][ti] = copy_farray(s,cnt);
  }else{
    r->s[ri][ti] = s;
  }

  r->cnt[ri][ti] = cnt;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_UTIL_RESP_STORE_F                          */
/*                                                                           */
/*  Copy the data if 'cpflag' = 1.                                           */
/*                                                                           */
/*****************************************************************************/
void mod_util_resp_store_f(r,k,data,n,cpflag,mylogf)
     struct response_struct *r;
     int k;          // Response index
     float *data;
     int n;
     int cpflag;
     char mylogf[];
{
  int ri,ti;
  char tstr[SLEN];

  if (data == NULL)  // Return without storing
    return;

  //ti = r->tsi;   // Trial index
  ti = r->gtsi;   // Trial index
  ri = r->ri[k]; // pointer within "s" response list
  //printf("  _store_f   Trial index:  %d   k %d  ri %d\n",ti,k,ri);
  if (r->f[ri][ti] != NULL){
    sprintf(tstr,"*** WARNING:  MOD_UTIL_RESP_STORE_F Already stored %s\n",
	    r->nd_name[k]);
    mylog(mylogf,tstr);
    // Allow overwrite in case of mm following a suspension
    sprintf(tstr,"    Trial %d results are being over-written\n",ti);
    mylog(mylogf,tstr);
    myfree(r->f[ri][ti]);
  }
  if (cpflag == 1)
    r->f[ri][ti] = copy_farray(data,n);
  else
    r->f[ri][ti] = data;

  r->fcnt[ri][ti] = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_UTIL_RESP_STORE_F_SAMP                       */
/*                                                                           */
/*  Original 'data' can be freed by caller, it is not retained here.         */
/*                                                                           */
/*****************************************************************************/
void mod_util_resp_store_f_samp(r,k,data,n,samp,mylogf)
     struct response_struct *r;
     int k;          // Response index
     float *data;
     int n;
     float samp;     // Sampling units in which 'data' is supplied
     char mylogf[];
{
  int xfact,nn;
  float ssamp,fxfact,diff,*odata;

  if (data == NULL)  // Return without storing
    return;

  ssamp = r->samp[k]; // Sampling to be used when storing in 'r'

  if (samp != ssamp){  // Need to convert
    fxfact = ssamp/samp;     // Expansion factor (float)
    xfact = my_rint(fxfact); // Expansion factor (int)
    if (fxfact < 1.0){
      // WYETH - no undersamping has been implemented yet
      mylog_exit(mylogf,"MOD_UTIL_RESP_STORE_F_SAMP  response subsampling\n");
    }
    diff = fxfact - (float)xfact;
    if ((diff > -0.0001) && (diff < 0.0001)){
      odata = over_sample_farray_interp(data,n,xfact,&nn);
    }else{
      mylog_exit(mylogf,"MOD_UTIL_RESP_STORE_F_SAMP  non-integer factor\n");
    }
  }else{
    odata = copy_farray(data,n);
    nn = n;
  }

  // WYETH - for original debugging.
  //printf("n = %d   nn = %d\n",n,nn);
  //append_farray_plot("zzz.pl","orig",data,n,1);
  //append_farray_plot("zzz.pl","OVER",odata,nn,1);
  //exit(0);

  mod_util_resp_store_f(r,k,odata,nn,0,mylogf); // 0-use original
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_UTIL_RESP_STORE_F_REG                        */
/*                                                                           */
/*  Make the data have regular sampling.                                     */
/*                                                                           */
/*  Do we know the appropriate length, so all traces will have same n???     */
/*                                                                           */
/*  This creates new storage.                                                */
/*                                                                           */
/*****************************************************************************/
void mod_util_resp_store_f_reg(r,k,xdata,ydata,n,samp_in,mylogf)
     struct response_struct *r;
     int k;           /* Response index */
     float *xdata;    /* x-coords, may be irregularly spaced */
     float *ydata;    /* y-coords */
     int n;
     float samp_in;  /* samples/s for stored trace */
     char mylogf[];
{
  int i;
  int nsamp;
  float *tx,*regx,*regy,samp_out,xf,*xtx,*xty;

  samp_out = r->samp[k];     // Samples/s for stored data

  // Adjust sampling units of x-axis if necessary
  if (samp_in != samp_out){
    tx = (float *)myalloc(n*sizeof(float));
    xf = samp_out/samp_in;
    for(i=0;i<n;i++)
      tx[i] = xdata[i] * xf;
  }else
    tx = xdata;

  nsamp = r->tsn * samp_out; // No. of regularly-spaced samples to store

  regx = (float *)myalloc(nsamp * sizeof(float));
  for(i=0;i<nsamp;i++)
    regx[i] = (float)i;

  // Sometimes, n=1 because the value didn't change during the simulation
  if (tx[n-1] < regx[nsamp-1]){ // If needed, add a sample beyond end
    xtx = (float *)myalloc((n+1)*sizeof(float));
    xty = (float *)myalloc((n+1)*sizeof(float));
    for(i=0;i<n;i++){
      xtx[i] = tx[i];
      xty[i] = ydata[i];
    }
    xtx[n] = 1.0 + regx[nsamp-1];
    xty[n] = xty[n-1];
    regy = resample_farray(xtx,xty,n+1,regx,nsamp);
    myfree(xtx); myfree(xty);
  }else
    regy = resample_farray(tx,ydata,n,regx,nsamp);

  mod_util_resp_store_f(r,k,regy,nsamp,0,mylogf);  // regy is saved here

  if (samp_in != samp_out)
    myfree(tx);

  myfree(regx);
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_UTIL_RESP_POP_UNIT                          */
/*                                                                           */
/*  Search for unit xi,yi,zi of population 'pop_name' in the response list.  */
/*  Return the index in the response lists, or -1 if not found.              */
/*                                                                           */
/*****************************************************************************/
int mod_util_resp_pop_unit(r,mylogf,pop_name,xi,yi,zi)
     struct response_struct *r;
     char mylogf[];
     char pop_name[];
     int xi,yi,zi;
{
  int i;
  int ti;

  ti = -1;
  for(i=0;i<r->n;i++){
    // 'save_pop_unit_as_'  or 'save_pop_grid_as_' ...
    if ((r->rformat[i] == 1)||(r->rformat[i] == 2)||(r->rformat[i] == 3)){
      if ((strcmp(r->plname[i],pop_name)==0) && 
	  (r->xi[i] == xi) && (r->yi[i] == yi) && (r->zi[i] == zi)){
	ti = i;
      }
    }
  }
  return ti;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_UTIL_CONE_WEIGHT                           */
/*                                                                           */
/*  Return the weight for cone at x,y given its type and the type code.      */
/*                                                                           */
/*****************************************************************************/
float mod_util_cone_weight(mylogf,x,y,cx,cy,sd,cid,ccode,distrib)
     char mylogf[];
     float x,y;        // Location of cone
     float cx,cy;      // Center of Gaussian weight function
     float sd;         // SD of Gaussian
     int cid;          // Cone type:  0-S, 1-M, 2-L
     int ccode;        // Code:  111-SML, 001-L only, ...
     int distrib;      // 0-Gauss, 1-unif
     //float prob;       // Prob. of making a connection  WYETH REMOVE
{
  int tc,sflag,mflag,lflag;
  float w,dx,dy;

  sflag = mflag = lflag = 0;
  tc = ccode;
  if (tc >= 100){
    sflag = 1;
    tc -= 100;
  }
  if (tc >= 10){
    mflag = 1;
    tc -= 10;
  }
  if (tc >= 1)
    lflag = 1;

  if (((cid == 0) && sflag) ||
      ((cid == 1) && mflag) ||
      ((cid == 2) && lflag)){
    if (distrib == 0){
      w = func_2d_gaussian(x,y,cx,cy,sd,sd,0);   // 1.0 at center
    }else if (distrib == 1){
      dx = cx - x;
      dy = cy - y;

      if ((dx*dx + dy*dy) < (sd*sd))
	w = 1.0;
      else
	w = 0.0;
    }
  }else{
    w = 0.0;
  }

  return w;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_UTIL_O_GET_FILT2D_GABOR                       */
/*                                                                           */
/*****************************************************************************/
float **mod_util_o_get_filt2d_gabor(mylogf,fo,sscale,rxn,ryn)
     char mylogf[];
     struct onode *fo;  // <filter>
     float sscale;      // (deg/pix)
     int *rxn;          // if -1, return filter length, else use this value
     int *ryn;          // if -1, return filter length, else use this value
{
  int xn,yn,hidnum;
  float theta,sdp_deg,sdo_deg,sf_deg,phase,xc_deg,yc_deg;
  float xc,yc,sdp,sdo,sf,a,**fdata;

  //******* WYETH - OTHER ..._O_GET_FILT... stuff was moved to "KERNEL_UTIL"
  //******* WYETH - OTHER ..._O_GET_FILT... stuff was moved to "KERNEL_UTIL"
  //******* WYETH - OTHER ..._O_GET_FILT... stuff was moved to "KERNEL_UTIL"
  //******* WYETH - OTHER ..._O_GET_FILT... stuff was moved to "KERNEL_UTIL"

  if (*rxn == -1){
    xn = onode_getpar_int_exit(fo,"xn");
  }else
    xn = *rxn;

  if (*ryn == -1){
    yn = onode_getpar_int_exit(fo,"yn");
  }else
    yn = *ryn;


  hidnum  = onode_getpar_int_dflt(fo,"hidden",0);   // Hidden code
  if (hidnum == 0){
    theta   = onode_getpar_flt_exit(fo,"theta");      // orientation (0..360)
    sdp_deg = onode_getpar_flt_exit(fo,"sd_par");     // SD parallel (deg)
    sdo_deg = onode_getpar_flt_exit(fo,"sd_orth");    // SD orthogonal (deg)
    sf_deg  = onode_getpar_flt_exit(fo,"sf");         // SF (cyc/deg)
    phase   = onode_getpar_flt_exit(fo,"phase");      // SF (0..360)
    xc_deg  = onode_getpar_flt_dflt(fo,"xc",0.0);     // centering (deg)
    yc_deg  = onode_getpar_flt_dflt(fo,"yc",0.0);     // centering (deg)
  }else if (hidnum == 1){
    theta   = 260.0;    // orientation (0..360)
    sdp_deg = 0.110;    // SD parallel (deg)
    sdo_deg = 0.075;    // SD orthogonal (deg)
    sf_deg  = 6.0;      // SF (cyc/deg)
    phase   = 0.0;      // Spatial phase (0..360)
    xc_deg  = (float)(*rxn)/2.0 * sscale - sdp_deg*3.0;  // centering (deg)
    yc_deg  = (float)(*ryn)/2.0 * sscale - sdp_deg*3.0;  // centering (deg)
  }else
    mylog_exit(mylogf,"MOD_UTIL_O_GET_FILT2D_GABOR  Bad 'hidden' index\n");


  xc  =  xc_deg/sscale;   // Convert (deg) -> (sampling units)
  yc  =  yc_deg/sscale;
  sdp = sdp_deg/sscale;
  sdo = sdo_deg/sscale;
  sf  =  sf_deg*sscale;

  xc += (float)(xn-1)/2.0;  // Center of grid
  yc += (float)(yn-1)/2.0;  // Center of grid

  fdata = gabor_2d_space_raw(xn,yn,xc,yc,sdo,sdp,sf,theta,phase);

  if (onode_item(fo,"max")==1){
    a = onode_getpar_flt_exit(fo,"max");
    make_max_const_2d_farray(fdata,xn,yn,a);
  }else if (onode_item(fo,"area")==1){
    a = onode_getpar_flt_exit(fo,"area");
    norm_area_2d_farray(fdata,xn,yn,a);
  }else{
    norm_area_2d_farray(fdata,xn,yn,1.0);  // Default
  }

  *rxn = xn;
  *ryn = yn;
  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_UTIL_O_GET_FILTER_2D                         */
/*                                                                           */
/*  Create and return the filter described by the <filter> object.           */
/*                                                                           */
/*****************************************************************************/
float **mod_util_o_get_filter_2d(mylogf,fo,sscale,rxn,ryn)
     char mylogf[];
     struct onode *fo;  // <filter>
     float sscale;      // (deg/pix)
     int *rxn,*ryn;     // if -1, return filter size, else use these values
{
  char *ftype,*dumpfile;
  float **fdata;

  //******* WYETH - OTHER ..._O_GET_FILT... stuff was moved to "KERNEL_UTIL"
  //******* WYETH - OTHER ..._O_GET_FILT... stuff was moved to "KERNEL_UTIL"
  //******* WYETH - OTHER ..._O_GET_FILT... stuff was moved to "KERNEL_UTIL"
  //******* WYETH - OTHER ..._O_GET_FILT... stuff was moved to "KERNEL_UTIL"
  // ************ WYETH OLD WAY ***************
  // ************ WYETH OLD WAY ***************
  // ************ WYETH OLD WAY ***************
  // ************ WYETH OLD WAY ***************
  // ************ WYETH OLD WAY ***************    NEW is BELOW
  // ************ WYETH OLD WAY ***************
  // ************ WYETH OLD WAY ***************

  mylog(mylogf,"  MOD_UTIL_O_GET_FILTER_2D\n");

  if (fo == NULL)
    mylogx(mylogf,"MOD_UTIL_O_GET_FILTER_2D","Null filter onode");

  ftype = onode_getpar_chr_exit(fo,"type");
  if (strcmp(ftype,"Gabor")==0){
    fdata = mod_util_o_get_filt2d_gabor(mylogf,fo,sscale,rxn,ryn);
  }else{
    mylogx(mylogf,"MOD_UTIL_O_GET_FILTER_2D","Unknown 2d <filter> type");
  }

  dumpfile = onode_getpar_chr_dflt(fo,"plot_filename",NULL);
  if (dumpfile != NULL){
    if ((strcmp(dumpfile,"null")!=0)&&(strcmp(dumpfile,"NULL")!=0))
      write_2d_data(dumpfile,fdata,0,0,*rxn,*ryn,4,2,1,0);
    myfree(dumpfile);
  }
 
  myfree(ftype);

  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_UTIL_STIM_SMOOTH                           */
/*                                                                           */
/*  Convolve the stimulus in space and/or time.                              */
/*                                                                           */
/*****************************************************************************/
void mod_util_stim_smooth(s,xn,yn,tn,sscale,tscale,pflag,mylogf)
     struct stim_struct *s;  // Stimulus params
     int xn,yn,tn;
     float sscale,tscale;
     int pflag;
     char *mylogf;
{
  float ***smooth;
  char ggstr[SLEN],*sform;

  //
  //  WYETH - THIS ONLY WORKS FOR '3d' SFORM  (3c is done in 'stim_util')
  //
  //  *** In the future perhaps this should also happen lower down, in
  //      stim_util ??  If different stim-types require special methods???
  //
  sform = paramfile_get_char_param_default(s->ppl,"stim_form","3d");
  if (strcmp(sform,"3d")!=0){
    myfree(sform);
    return;
  }

  smooth = stim_util_blur(mylogf,s->ppl,s->d,xn,yn,tn,sscale,tscale);
  if (smooth != NULL){
    free_f3tensor(s->d,1,xn,1,yn,1,tn);  // Free original stimulus
    s->d = smooth;                       // Replace with smoothed stimulus
  }
      
  myfree(sform);
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_UTIL_GET_3D_STIM                          */
/*                                                                           */
/*  Stimulus data is [x][y][t].                                              */
/*                                                                           */
/*  Get stimulus data for the 'k'th trial.  If 'k' is -1, var params are     */
/*  ignored.                                                                 */
/*                                                                           */
/*****************************************************************************/
void mod_util_get_3d_stim(s,k,xn,yn,tn,sscale,tscale,tsamp,mylogf)
     struct stim_struct *s;  // Stimulus params
     int k,xn,yn,tn;
     float sscale,tscale;
     int tsamp;
     char *mylogf;
{
  int pflag;
  float ***d;
  char *stimtype,*stimadapt,ggstr[LONG_SLEN];

  pflag = 0;
  if (mylogf == NULL)
    pflag = 1;

  // Next line does nothing if k < 0
  mod_util_substitute_var_params(s,k,mylogf,0,NULL,NULL);

  d = get_stim_3d_ppl(mylogf,s,xn,yn,tn,sscale,tscale,tsamp);

  //
  //  Check if NULL stim is OK.
  //
  stimtype = paramfile_get_char_param_or_exit(s->ppl,"stim_type");
  if (d == NULL){
    if (strcmp(stimtype,"null")==0){
      s->d = d;
      return;
    }else{
      sprintf(ggstr,"stimtype:  %s\n",stimtype);
      mylog(mylogf,ggstr);
      mylog_exit(mylogf,"MOD_UTIL_GET_3D_STIM  Stim is NULL\n");
    }
  }
  myfree(stimtype);


  //
  //  Spatial Blurring:   use 3D conv for blurring here
  //
  // printf("*** WYETH HERE DELTA_FUNC STIM!!!! ********** \n");
  // s->d = get_stim_impulse(xn,yn,tn,5,5,2,10.0);
  s->d = d;
  mod_util_stim_smooth(s,xn,yn,tn,sscale,tscale,pflag,mylogf);
  d = s->d;  // In case old one was freed

  //
  //  Stimulus adaptation
  //
  stimadapt = paramfile_get_char_param_default(s->ppl,"stim_adapt","none");
  if (strcmp(stimadapt,"meangray")==0){
    int i,j,k;
    int nmask,cn;
    float mu,adtau,*emask,*td,*ad;

    mylog(mylogf,"  *** Subjecting stimulus to adaptation.\n");

    adtau = paramfile_get_float_param_default(s->ppl,"stim_adapt_tau",-1.0);

    if (adtau < 0.0){
      for(i=1;i<=xn;i++){
	for(j=1;j<=yn;j++){
	  mu = mean_farray(&(d[i][j][1]),tn);
	  for(k=1;k<=tn;k++)
	    d[i][j][k] += 0.5 - mu;
	}
      }
    }else{
      nmask = (int)(adtau / tscale) * 10 + 1; /* Must be odd */
      cn = (nmask-1)/2;

      /*** WYETH - something strangehere - the line below doesn't work on 2nd
	time around.  I replaced it with the get_zero_farray call. ???? ***/
      /*emask = (float *)myalloc(nmask*sizeof(float));*/

      emask = get_zero_farray(nmask);
      for(i=cn;i<nmask;i++)
	emask[i] = func_one_sided_exp((double)(i-cn),0.0,
				      (double)(adtau/tscale));
      norm_area_farray(emask,nmask,1.0);
      append_farray_plot("zzz.emask","emask",emask,nmask,1);

      td = (float *)myalloc((cn+tn)*sizeof(float));
      for(i=1;i<=xn;i++){
	for(j=1;j<=yn;j++){

	  for(k=0;k<cn;k++){  /*** Pad beginning with initial value ***/
	    td[k] = d[i][j][1];
	  }
	  for(k=1;k<=tn;k++)  /*** Pad beginning with initial value ***/
	    td[k+cn-1] = d[i][j][k];
	  if ((i==10)&&(j==10))
	    append_farray_plot("zzz.emask","before",td+cn,tn,1);
	  ad = fft_convolve(td,cn+tn,emask,nmask,1);
	  if ((i==10)&&(j==10))
	    append_farray_plot("zzz.emask","adstim",ad+cn,tn,1);
	  
	  for(k=1;k<=tn;k++)
	    d[i][j][k] += 0.5 - ad[cn+k-1];

	  if ((i==10)&&(j==10))
	    append_farray_plot("zzz.emask","FINAL",&(d[i][j][1]),tn,1);

	  myfree(ad);
	}
      }
      myfree(emask);
      myfree(td);
    }
  }
  myfree(stimadapt);

  s->d = d;

  s->d_r  = NULL; // No stimulus for right eye
  s->dc   = NULL;
  s->dc_r = NULL;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_UTIL_GET_3D_STIM_BINOC                        */
/*                                                                           */
/*  Create  s->d[xn][yn][tn]                                                 */
/*          s->d_r[xn][yn][tn]                                               */
/*                                                                           */
/*  Get stimulus data for the 'k'th trial.  If 'k' is -1, var params are     */
/*  ignored.                                                                 */
/*                                                                           */
/*****************************************************************************/
void mod_util_get_3d_stim_binoc(s,k,xn,yn,tn,sscale,tscale,tsamp,mylogf)
     struct stim_struct *s;  // Stimulus params
     int k,xn,yn,tn;
     float sscale,tscale;
     int tsamp;
     char *mylogf;
{
  int pflag;
  char *stimtype,ggstr[LONG_SLEN];

  pflag = 0;
  if (mylogf == NULL)
    pflag = 1;

  // Next line does nothing if k < 0
  mod_util_substitute_var_params(s,k,mylogf,0,NULL,NULL);

  stimtype = paramfile_get_char_param_or_exit(s->ppl,"stim_type");

  get_stim_3d_b_ppl(mylogf,stimtype,s->ppl,xn,yn,tn,sscale,tscale,tsamp,
		    &(s->d),&(s->d_r));

  if ((s->d == NULL)||(s->d_r == NULL)){
    sprintf(ggstr,"stimtype %s\n",stimtype);
    mylog(mylogf,ggstr);
    mylog_exit(mylogf,"MOD_UTIL_GET_3D_STIM_BINOC  Stim is NULL\n");
  }
  myfree(stimtype);


  // use 3D conv for blurring here - WYETH - see MOD_UTIL_GET_3D_STIM above

  // stim_adapt here - WYETH - see MOD_UTIL_GET_3D_STIM above


  s->dc   = NULL;
  s->dc_r = NULL;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_UTIL_GET_3C_STIM                           */
/*                                                                           */
/*  Stimulus data is [x][y][t].                                              */
/*                                                                           */
/*  Get stimulus data for the 'k'th trial.  If 'k' is -1, var params are     */
/*  ignored.                                                                 */
/*                                                                           */
/*****************************************************************************/
void mod_util_get_3c_stim(s,k,xn,yn,tn,sscale,tscale,tsamp,mylogf)
     struct stim_struct *s;  // Stimulus params
     int k,xn,yn,tn;
     float sscale,tscale;
     int tsamp;
     char *mylogf;
{
  int pflag;
  int ***d;
  char *stimtype,*stimadapt,ggstr[LONG_SLEN];

  pflag = 0;
  if (mylogf == NULL)
    pflag = 1;

  // Next line does nothing if k < 0
  mod_util_substitute_var_params(s,k,mylogf,0,NULL,NULL);

  stimtype = paramfile_get_char_param_or_exit(s->ppl,"stim_type");

  d = get_stim_3c_ppl(mylogf,stimtype,s->ppl,xn,yn,tn,sscale,tscale,tsamp);

  if (d == NULL){
    sprintf(ggstr,"stimtype %s\n",stimtype);
    mylog(mylogf,ggstr);
    mylog_exit(mylogf,"MOD_UTIL_GET_3C_STIM  Stim is NULL\n");
  }

  s->dc = d;

  // If smoothing has been requested, this will replace 's->dc'
  mod_util_stim_smooth(s,xn,yn,tn,sscale,tscale,pflag,mylogf);  // ONLY FOR 3D


  s->dc_r = NULL;
  s->d    = NULL;
  s->d_r  = NULL; // No stimulus for right eye
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_UTIL_GET_3RGB_STIM                          */
/*                                                                           */
/*  Stimulus data is [x][y][t][rgb].                                         */
/*                                                                           */
/*  Get stimulus data for the 'k'th trial.  If 'k' is -1, var params are     */
/*  ignored.                                                                 */
/*                                                                           */
/*****************************************************************************/
void mod_util_get_3rgb_stim(s,k,xn,yn,tn,sscale,tscale,tsamp,mylogf)
     struct stim_struct *s;  // Stimulus params
     int k,xn,yn,tn;
     float sscale,tscale;
     int tsamp;
     char *mylogf;
{
  int pflag;
  float ****d;
  char *stimtype,*stimadapt,ggstr[LONG_SLEN];

  pflag = 0;
  if (mylogf == NULL)
    pflag = 1;

  // Next line does nothing if k < 0
  mod_util_substitute_var_params(s,k,mylogf,0,NULL,NULL);

  stimtype = paramfile_get_char_param_or_exit(s->ppl,"stim_type");

  d = get_stim_3rgb_ppl(mylogf,s,stimtype,s->ppl,xn,yn,tn,sscale,tscale,tsamp);

  if (d == NULL){
    sprintf(ggstr,"stimtype %s\n",stimtype);
    mylog(mylogf,ggstr);
    mylog_exit(mylogf,"MOD_UTIL_GET_3D_RGB_STIM  Stim is NULL\n");
  }

  s->drgb_l = d;

  // If smoothing has been requested, this will replace 's->dc'
  // *** Note - this smoothing currently (2017 Aug) does nothing unless
  //            the stim_form is "3d"
  mod_util_stim_smooth(s,xn,yn,tn,sscale,tscale,pflag,mylogf);  // ONLY FOR 3D

  s->drgb_r  = NULL;
  s->dc      = NULL;
  s->dc_r    = NULL;
  s->d       = NULL;
  s->d_r     = NULL; // No stimulus for right eye
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_UTIL_GET_1D_STIM                           */
/*                                                                           */
/*  Stimulus data is [t].                                                    */
/*                                                                           */
/*  Get stimulus data for the 'k'th trial.  If 'k' is -1, var params are     */
/*  ignored.                                                                 */
/*                                                                           */
/*****************************************************************************/
float *mod_util_get_1d_stim(s,k,tn,tscale,tsamp,mylogf)
     struct stim_struct *s;  // Stimulus params
     int k,tn;
     float tscale;
     int tsamp;
     char *mylogf;
{
  float *d;
  char *stimtype,ggstr[LONG_SLEN];

  // Next line does nothing if k<0
  mod_util_substitute_var_params(s,k,mylogf,0,NULL,NULL); 

  stimtype = paramfile_get_char_param_or_exit(s->ppl,"stim_type");
  if (strcmp(stimtype,"pulse1d")==0){
    d = get_stim_1d_ppl_pulse(s->ppl,tn,tscale,tsamp);
  }else if (strcmp(stimtype,"sine1d")==0){
    d = get_stim_1d_ppl_sine(s->ppl,tn,tscale,tsamp);
  }else if (strcmp(stimtype,"vsksos1d")==0){
    d = get_stim_1d_ppl_vsksos(s->ppl,tn,tscale,tsamp);
  }else if (strcmp(stimtype,"noise1d")==0){
    d = get_stim_1d_ppl_noise(s->ppl,tn,tscale,tsamp);
  }else if (strcmp(stimtype,"stat1")==0){
    d = get_stim_1d_ppl_stat1(s->ppl,tn,tscale,tsamp);
  }else{
    d = NULL;
    sprintf(ggstr,"stimtype %s\n",stimtype);
    mylog(mylogf,ggstr);
    mylog_exit(mylogf,"MOD_UTIL_GET_1D_STIM_DATA  Unknown stim_type\n");
  }
  myfree(stimtype);

  return d;
}
