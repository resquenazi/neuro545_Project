/*****************************************************************************/
/*                                                                           */
/*  stim_util.c                                                              */
/*  wyeth bair                                                               */
/*  NYU                                                                      */
/*  05/29/96                                                                 */
/*                                                                           */
/*  Routines here are concerned with creating a 3D stimulus that will be     */
/*  used in a model of motion processing.  Some routines are concerned with  */
/*  defining the model.                                                      */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "my_util.h"
#include "misc_util.h"
#include "nr_util.h"
#include "carray_util.h"
#include "iarray_util.h"
#include "farray_util.h"
#include "func_util.h"
#include "myrand_util.h"
#include "data_util.h"
#include "paramfile_util.h"
#include "fft_util.h"
#include "sig_util.h"
#include "noise_util.h"
#include "stm_util.h"
#include "kernel_util.h"
#include "stim_util.h"

#include "paramfile.h"
#include "mod.h"
#include "stmh.h"


#define VSK_sos_unit 0.015258789 // 1/65.536 Hz (Victor/Shapley/Knight 1977)

#define NPAT 256 /* For making spatial pattern for noise stimuli. */
/*** This table will have positive and negative ampls, unlike wy2patx.c
  in which ttab holds only positive ampls, and negative ampls are created
  by shifting. ***/
#define NAMP 257 /* +/- 128 levels and 0;  Like 2*NTTAB-1 in wy2patx.c */

/*** For sum of sinusoid stimuli of  Victor, Shapley, Knight (1977) ***/
/*** In units of 1/65.536 Hz ***/
int VSK_sos_n[3] = {6,8,8}; /* Number of entries in stimuli 0,1,2 */
int VSK_sos_0[6] = {42,72,162,352,802,1402};
int VSK_sos_1[8] = {14,30,62,126,254,510,1022,2046};
int VSK_sos_2[8] = {15,31,63,127,255,511,1023,2047};

/**************************************-**************************************/
/*                                                                           */
/*                           GET_CROSS_PRODUCT_PARAMS                        */
/*                                                                           */
/*  Get the n-dimensional cross-product list for the given numbers of        */
/*  items.  Also, if 'nrpt' > 1, then repeat each unique combination that    */
/*  many times.                                                              */
/*                                                                           */
/*  n - the number of component sets.                                        */
/*  cnt - the number of items in each component set [n].                     */
/*  list - the array of all cross-product sets [ntot][n].                    */
/*                                                                           */
/*****************************************************************************/
void get_cross_product_params(mylogf,n,cnt,val,nrpt,nvp,nvl,nlinkval,nvs,
			      rn,rlist)
     char *mylogf;
     int n;          // Number of parameters
     int *cnt;       // [n] number of values for each parameter
     int nrpt;       // number of times to repeat each unique selection
     int nvp;        // number of linked pairs
     int nvl;        // number of linked variable parameters
     int nlinkval;   // number of values for linked variable params.
     int nvs;        // number of VARSINGLE values, 0 or 1
     char ***val;    // Values, [n][cnt[i]]
     int *rn;        // length of list returned
     char ****rlist; // [*rn][n] var param values (except VARSINGLE is empty)
{
  int i,j,k,l;
  int **list,ntot,nlist,t,tn,nsing,*cnt_new,li;
  char ***vv;

  if (n==0)
    mylog_exit(mylogf,"GET_CROSS_PRODUCT_PARAMS  No variable params\n");

  nsing = n - nvp*2 - nvl - nvs; // Number of single variable dimensions

  // Calculate the number of independent dimensions, 'tn'
  if (nvl > 1)
    tn = nsing + nvp + 1; // Count all linked vars as one
  else
    tn = nsing + nvp;

  //printf("tn = %d\n",tn);

  //
  //  Must make new 'cnt' list if there are pairs.
  //  The 'cnt' list tells the number of values for each independent axis
  //
  cnt_new = (int *)myalloc(tn*sizeof(int));
  for(i=0;i<nsing;i++)
    cnt_new[i] = cnt[i];
  for(i=0;i<nvp;i++)
    cnt_new[nsing+i] = cnt[nsing + 2*i];
  if (nvl > 0)
    cnt_new[tn-1] = cnt[n-nvs-1];


  // Data for linked variables is stored at end of lists

  // ********** WYETH - I HAVE REVERSED TRIAL ORDER for consistency with
  //   ******** WMPI PYTHON interface 2019 Nov 8th
  //get_nd_cross_product_list_iarray(tn,cnt_new,&list,&nlist);
  get_nd_cross_product_list_iarray_rev(tn,cnt_new,&list,&nlist);
  //ntot = nlist * nrpt;
  ntot = (nlist + nvs) * nrpt;  // Add one extra stim per repeat for 'nvs'
  vv = get_2d_pointer_carray(ntot,n);

  t = 0; // Trial number
  for(i=0;i<nlist;i++){
    for(j=0;j<nrpt;j++){

      li = 0;
      for(k=0;k<nsing;k++){  // Singles
	l = list[i][li];
	vv[t][k] = strdup(val[k][l]);
	li += 1;
      }
      for(k=0;k<nvp;k++){  // Pairs
	l = list[i][li];
	vv[t][nsing+2*k]   = strdup(val[nsing+2*k][l]);
	vv[t][nsing+2*k+1] = strdup(val[nsing+2*k+1][l]);
	li += 1;
      }

      if ((li != tn) && (nvl == 0)){
	printf("li = %d    tn == %d  nvl =%d\n",li,tn,nvl);
	mylog_exit(mylogf,"GET_CROSS_PRODUCT_PARAMS  bad values, nvl = 0\n");
      }else if ((li != (tn-1)) && (nvl > 0)){
	printf("li = %d    tn == %d\n",li,tn);
	mylog_exit(mylogf,"GET_CROSS_PRODUCT_PARAMS  bad values, nvl > 0\n");
      }

      l = list[i][li];
      for(k=0;k<nvl;k++){    // 'varlinks'
	vv[t][nsing+2*nvp+k] = strdup(val[nsing+2*nvp+k][l]);
      }

      // If there is a final VARSINGLE_ param, then use the default value
      if (nvs == 1)
	vv[t][n-1] = strdup(val[n-1][0]); // Use 0th (thus default) value

      t += 1;
    }
  }
  free_2d_iarray(list,nlist);

  //
  //  *** NOTE, all 'VARSINGLE_' stimuli are at end of list, and are 
  //  left empty here, to be filled in by the caller.
  //

  myfree(cnt_new);

  *rn = ntot;
  *rlist = vv;
}
/**************************************-**************************************/
/*                                                                           */
/*                          STIM_UTIL_GET_VARGEN_LIST                        */
/*                                                                           */
/*****************************************************************************/
int stim_util_get_vargen_list(mylogf,pp,rval)
     char *mylogf;
     struct param_pair *pp;
     char ***rval;
{
  int i;
  int n,dec,seed;
  float m,a,x;
  char *gtype,**val,pstr[SLEN],temp[SLEN];

  gtype = paramfile_get_nth_char_param_pp_or_exit(pp,0);
  dec   = paramfile_get_nth_int_param_ptr_or_exit(pp,1);
  n     = paramfile_get_nth_int_param_ptr_or_exit(pp,2);
  m     = paramfile_get_nth_float_param_ptr_or_exit(pp,3);
  a     = paramfile_get_nth_float_param_ptr_or_exit(pp,4);
  seed  = paramfile_get_nth_int_param_ptr_or_exit(pp,5);

  if (seed > 0)
    seed = -seed;

  val = (char **)myalloc(n*sizeof(char *));

  if (strcmp(gtype,"uniform")==0){
    if (dec > 0)
      sprintf(pstr,"%%.%df",dec);
    else
      sprintf(pstr,"%%d");

    for(i=0;i<n;i++){
      x = a + m*myrand_util_ran2(&seed);
      if (dec > 0)
	sprintf(temp,pstr,x);
      else
	sprintf(temp,pstr,(int)x);
      val[i] = strdup(temp);
    }
  }else
    mylog_exit(mylogf,"STIM_UTIL_GET_VARGEN_LIST  Unknown type\n");

  myfree(gtype);
  
  *rval = val;
  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                        STIM_UTIL_GET_VARGENPAIR_LIST                      */
/*                                                                           */
/*****************************************************************************/
int stim_util_get_vargenpair_list(mylogf,pp,rval,rval2)
     char *mylogf;
     struct param_pair *pp;
     char ***rval;          // List of values for 1st parameter
     char ***rval2;         // List of values for 2nd parameter
{
  int i,j,k;
  int n,ti,tio,*tap,tapn,seed,*seedlist;
  char *gtype,**val1,**val2,temp[SLEN];

  n     = paramfile_get_nth_int_param_ptr_or_exit(pp,1);
  gtype = paramfile_get_nth_char_param_pp_or_exit(pp,2);

  if (strcmp(gtype,"list_opp_mseq_tap_11")==0){

    // starting index for first paramname
    k = paramfile_get_nth_int_param_ptr_or_exit(pp,3);

    // Get pointer to the order 11 tap array
    get_mseq_tap_array(11,&tap,&tapn);

    if ((k >= tapn) || (k < 0))
      mylog_exit(mylogf,"STIM_UTIL_GET_VARGENPAIR_LIST  Index error\n");

    val1 = (char **)myalloc(n*sizeof(char *));
    val2 = (char **)myalloc(n*sizeof(char *));

    ti = k;
    for(i=0;i<n;i++){

      sprintf(temp,"%d",tap[ti]);
      val1[i] = strdup(temp);

      tio = ti+ tapn/2;  // Tap index opposite (halfway round the list)
      if (tio >= tapn)
	tio -= tapn;

      sprintf(temp,"%d",tap[tio]);
      val2[i] = strdup(temp);

      ti += 1;
      if (ti == tapn)
	ti = 0;
    }
  }else if (strcmp(gtype,"unif_100000")==0){

    // Get random seed to choose seeds
    seed = paramfile_get_nth_int_param_ptr_or_exit(pp,3);

    // Get enough seeds to fill two lists of length 'n'
    seedlist = get_seeds(seed,100000,2*n);

    val1 = (char **)myalloc(n*sizeof(char *));  // Create the lists
    val2 = (char **)myalloc(n*sizeof(char *));

    for(i=0;i<n;i++){
      sprintf(temp,"%d",seedlist[2*i]);
      val1[i] = strdup(temp);

      sprintf(temp,"%d",seedlist[2*i+1]);
      val2[i] = strdup(temp);
    }
    myfree(seedlist);
  }else{
    printf(" *** gtype = %s\n",gtype);
    mylog_exit(mylogf,"STIM_UTIL_GET_VARGENPAIR_LIST  Unknown type\n");
  }

  myfree(gtype);

  *rval  = val1;
  *rval2 = val2;
  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                            STIMU_VARSEQ_GET_LIST                          */
/*                                                                           */
/*****************************************************************************/
int stimu_varseq_get_list(mylogf,pp,rval)
     char *mylogf;
     struct param_pair *pp;
     char ***rval;
{
  int i;
  int n,ndec,fv2si,ival,done;
  char *cv1,*cv2,*seqtype,**vlist;
  char tstr[SLEN],fstr[SLEN];
  float finc,fv1,fv2,fscale,x;

  n = 0;
  vlist = NULL;

  seqtype = paramfile_get_nth_char_param_pp_or_exit(pp,0);

  if (strcmp(seqtype,"range_inc")==0){
    cv1 = paramfile_get_nth_char_param_pp_or_exit(pp,1);
    cv2 = paramfile_get_nth_char_param_pp_or_exit(pp,2);
    finc = paramfile_get_nth_float_param_ptr_or_exit(pp,3);
    //printf("cv1 = %s  cv2 = %s   finc = %f\n",cv1,cv2,finc);

    ndec = count_decimals(cv1);

    fscale = pow((double)10.0,(double)ndec);  // Scale factor for rounding

    fv1 = atof(cv1);                // start value, float
    fv2 = atof(cv2);                // end value, float
    fv2si = my_rint(fv2 * fscale);  // end value, scaled, rounded to int

    //
    //  Count the number of values
    //
    x = fv1;
    done = 0;
    while(done == 0){
      ival = my_rint(x * fscale);  // Scaled, integer value
      if (ival > fv2si){
	done = 1;
      }else if (ival == fv2si){
	n += 1;
	done = 1;
      }else{
	n += 1;
	x += finc;
      }
    }
    if (n <= 0){
      exit_error("STIMU_VARSEQ_GET_N","No values defined by sequence");
    }
    myfree(cv1);
    myfree(cv2);

    //
    //  Store the values as character strings in 'vlist[n]'
    //
    if (rval != NULL){
      vlist = (char **)myalloc(n*sizeof(char *));
      x = fv1;
      i = 0;
      done = 0;
      while(done == 0){
	ival = my_rint(x * fscale);  // Scaled, integer value
	if (ival > fv2si){
	  done = 1;
	}else{
	  sprintf(fstr,"%%.%df",ndec);
	  sprintf(tstr,fstr,x);
	  //printf("  tstr ==>%s<==\n",tstr);
	  vlist[i] = strdup(tstr);
	  i += 1;
	  x += finc;
	  if (ival == fv2si){
	    done = 1;
	  }
	}
      }
    }
  }else{
    printf("  *** sequence type:  %s\n",seqtype);
    exit_error("STIMU_VARSEQ_GET_LIST","Unknown sequence type.");
  }

  myfree(seqtype);

  if (rval != NULL)
    *rval = vlist;

  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                             STIMU_VARSEQ_GET_N                            */
/*                                                                           */
/*****************************************************************************/
int stimu_varseq_get_n(mylogf,pp)
     char *mylogf;
     struct param_pair *pp;
{
  int n;

  n = stimu_varseq_get_list(mylogf,pp,NULL);

  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                               FREE_VAR_PARAMS                             */
/*                                                                           */
/*****************************************************************************/
void free_var_params(mylogf,vval,vname,vn,n,val,vcnt)
     char *mylogf;
     char **vname;          // var param names [vn]
     char ***vval;          // var param values for all trials [n][vn]
     int vn;                // number of var params
     int n;                 // total number of trials
     char ***val;           // Unique values [vn][vcnt[]]
     int *vcnt;             // Number of unique values [vn]
{
  mylog(mylogf,"  FREE_VAR_PARAMS\n");

  free_2d_pointer_carray(vval,n,vn);
  free_2d_carray(vname,vn);
  free_2d_pointer_var_count_carray(val,vn,vcnt);
  myfree(vcnt);
}
/**************************************-**************************************/
/*                                                                           */
/*                             PREP_VPAR_VFILE_VL                            */
/*                                                                           */
/*  Get var params and values from a VARFILE, or INLINE, and format them     */
/*  to be put into the VARLINK spot within PREP_VAR_PARAMS.                  */
/*                                                                           */
/*****************************************************************************/
void prep_vpar_vfile_vl(mylogf,sppl,iflag,rnvar,rname,rnval,rval)
     char *mylogf;
     struct param_pair_list *sppl;   // Stimulus parameter list
     int iflag;               // 0-VARFILE, 1-INLINE VAR_TABLE
     int *rnvar;              // number of var params
     char ***rname;           // var param names [rnvar]
     int *rnval;              // Number of linked values
     char ****rval;           // values for all params [rnvar][rnval]
{
  FILE *fopen(),*fin;
  int i,j;
  int nvar,nstim,ns;
  char *infile,ts[SLEN],**vname,***val,**slist;

  if (iflag == 0){
    //
    //  Open the VARFILE
    //
    infile = paramfile_get_char_param_or_exit(sppl,"VARFILE");
    if ((fin = fopen(infile,"r")) == NULL){
      printf("  *** Cannot open file %s.\n",infile);
      mylog_exit(mylogf,"Cannot open VARFILE filename'\n");
    }
    sprintf(ts,"    Reading %s\n",infile);
    mylog(mylogf,ts);


    //
    //  Read param number and names
    //
    ns = fscanf(fin,"%s %d",ts,&nvar);
    if (strcmp(ts,"npar")!=0)
      mylog_exit(mylogf,"Expecting constant string 'npar'\n");

    vname = (char **)myalloc(nvar*sizeof(char *));
    for(i=0;i<nvar;i++){
      ns = fscanf(fin,"%s",ts);
      vname[i] = strdup(ts);
    }

    //
    //  Read param values for each stimulus
    //
    ns = fscanf(fin,"%s %d",ts,&nstim);
    if (strcmp(ts,"nstim")!=0)
      mylog_exit(mylogf,"Expecting constant string 'nstim'\n");

    val = get_2d_pointer_carray(nvar,nstim);
    for(i=0;i<nstim;i++){
      for(j=0;j<nvar;j++){
	ns = fscanf(fin,"%s",ts);
	val[j][i] = strdup(ts);
      }
    }

    fclose(fin);
    myfree(infile);

  }else{
    //  The var table data came from 'INLINE VAR_TABLE' and is stored
    //  in 'ppl->s' with the following format (example for 3 pars):
    //
    //   i   Data string
    //  ---  ----------------------
    //  [0]  npar 3
    //  [1]  si shape_id rotation
    //  [2]  nstim 370
    //  [3]  1 1 0
    //  [4]  2 2 0
    //  [5]  3 3 0
    //  ...  ...

    //
    //  Extract param number and names from ppl string data
    //
    get_items_from_string(sppl->s[1],&slist,&nvar);
    vname = (char **)myalloc(nvar*sizeof(char *));
    for(i=0;i<nvar;i++)
      vname[i] = strdup(slist[i]);
    free_2d_carray(slist,nvar);

    //
    //  Read param values for each stimulus
    //
    sscanf(sppl->s[2],"%*s %d",&nstim);

    val = get_2d_pointer_carray(nvar,nstim);
    for(i=0;i<nstim;i++){

      get_items_from_string(sppl->s[3+i],&slist,&ns);

      if (ns != nvar){
	printf("ns = %d   nvar = %d\n",ns,nvar);
	exit_error("PREP_VPAR_VFILE","Mismatch in INLINE data count");
      }
      
      for(j=0;j<nvar;j++)
	val[j][i] = strdup(slist[j]);

      free_2d_carray(slist,ns);
    }
  }

  *rnvar = nvar;
  *rname = vname;
  *rnval = nstim;
  *rval = val;       // var param values for all trials [nvar][nstim]
}
/**************************************-**************************************/
/*                                                                           */
/*                              PREP_VPAR_VFILE                              */
/*                                                                           */
/*  Get var params and values from a VARFILE.                                */
/*                                                                           */
/*  This was written in Sep 2014 to allow a long list of linked var params   */
/*  for Pasupathy Lab shape stimuli.                                         */
/*                                                                           */
/*****************************************************************************/
void prep_vpar_vfile(mylogf,sppl,rptflag,iflag,rvval,rvname,rvtype,rvn,rn,rval,
		     rvcnt,rnvl,rnvp)
     char *mylogf;
     struct param_pair_list *sppl;   // Stimulus parameter list
     int rptflag;             // 0-ignore 'stim_nrpt', 1-use 'stim_nrpt
     int iflag;               // 0-VARFILE, 1-INLINE VAR_TABLE
     char ****rvval;          // var param values for all trials [rn][rvn]
     char ***rvname;          // var param names
     char **rvtype;           // var param types
     int *rvn;                // number of var params
     int *rn;                 // total number of trials
     char ****rval;           // Unique values [rvn][rvcnt[]]
     int **rvcnt;             // Number of unique values [rvn]
     int *rnvl;               // Number of linked values
     int *rnvp;               // Number of linked pairs
{
  FILE *fopen(),*fin;
  int i,j;
  int nvar,nstim,*vcnt,nvl,nvp,tvn,nrpt,nlinkval,nvs,ns;
  char *infile,ts[SLEN],**vname,***tvval,***val,*vtype,**slist;

  // *************** WYETH OLD WAY REMOVE ********** 2019 June 4
  // *************** WYETH OLD WAY REMOVE ********** 2019 June 4
  // *************** WYETH OLD WAY REMOVE ********** 2019 June 4

  if (rptflag == 1)
    nrpt = paramfile_get_int_param_or_exit(sppl,"stim_nrpt");
  else
    nrpt = 1; // Added for labradoodle

  prep_vpar_vfile_vl(mylogf,sppl,iflag,&nvar,&vname,&nstim,&val);

  nvp = 0;     // No linked pairs
  nvl = nvar;  // Number of linked vars is 'nvar'
  vcnt = get_const_iarray(nvar,nstim);
  nlinkval = nstim;
  nvs = 0;     // No varsingle values


  // Set types for 'var' params
  vtype = (char *)myalloc(nvar*sizeof(char));
  for(i=0;i<nvar;i++){
    vtype[i] = paramfile_get_type_from_comment(sppl,vname[i],'f');
  }

  get_cross_product_params(mylogf,nvar,vcnt,val,nrpt,nvp,nvl,nlinkval,nvs,
			   &tvn,&tvval);

  *rvval = tvval; *rvname = vname; *rvn = nvar; *rn = tvn; *rvtype = vtype;
  *rval = val; *rvcnt = vcnt; *rnvl = nvl;  // Added for labradoodle
  *rnvp = nvp;
}
/**************************************-**************************************/
/*                                                                           */
/*                               PREP_VAR_PARAMS                             */
/*                                                                           */
/*  Determine the number of variable parameters and the total number of      */
/*  trials, including repeats.  Return values for var params for all trials  */
/*  and names of var params and their types.                                 */
/*                                                                           */
/*  VAR_<pname> <val> <val> ... <val>                                        */
/*  VARGEN_<pname> <type> <decimals> <n> <mult> <add> <seed>                 */
/*  VARLINK_<pname> <val> <val> ... <val>                                    */
/*                                                                           */
/*  NOTES:                                                                   */
/*    - Put 'varlink' params last in list.                                   */
/*                                                                           */
/*****************************************************************************/
void prep_var_params(mylogf,sppl,rptflag,rvval,rvname,rvtype,rvn,rn,rval,
		     rvcnt,rnvl,rnvp)
     char *mylogf;
     struct param_pair_list *sppl;   // Stimulus parameter list
     int rptflag;             // 0-ignore 'stim_nrpt', 1-use 'stim_nrpt
     char ****rvval;          // var param values for all trials [rn][rvn]
     char ***rvname;          // var param names
     char **rvtype;           // var param types
     int *rvn;                // number of var params
     int *rn;                 // total number of trials
     char ****rval;           // Unique values [rvn][rvcnt[]]
     int **rvcnt;             // Number of unique values [rvn]
     int *rnvl;               // Number of linked values
     int *rnvp;               // Number of linked pairs
{
  int i,j,k;
  int n,nprod,nvar,nv,nvg,nvp,nvl,nlinkval,nrpt,tvn,*vcnt,nvs,nvq,vlf_flag;
  char **vname,***tvval,***val,*vtype,tstr[SLEN],*par2,*istr,*stimtype;
  char **vlf_name,***vlf_val;
  struct param_pair *pp;
  char ggstr[LONG_SLEN];

  mylog(mylogf,"  PREP_VAR_PARAMS\n");

  //
  //  ********** WYETH - this 'VARFILE' (and INLINE) should be made to load
  //  **********   into the varlinks (below), and thus it could be used with
  //  **********   other VAR_... lines
  //

  vlf_flag = 0; // Assume no VARFILE
  if (paramfile_test_param(sppl,"VARFILE")){
    //
    //  If this .stm file uses a VARFILE, handle that separately and return
    //
    vlf_flag = 1;
    if (1){
      prep_vpar_vfile_vl(mylogf,sppl,0,&nvl,&vlf_name,&nlinkval,&vlf_val);
    }else{
      /// OLD WAY - REMOVE ***
      prep_vpar_vfile(mylogf,sppl,rptflag,0,rvval,rvname,rvtype,rvn,rn,rval,
		      rvcnt,rnvl,rnvp);
      return;
    }
  }else if (paramfile_test_param(sppl,"INLINE")){
    //
    //  If this is an "INLINE VAR_TABLE", then call the VARFILE code with
    //    a flag value of '1'.
    //
    istr = param_getc_exit(sppl,"INLINE");
    if (strcmp(istr,"VAR_TABLE")==0){
      vlf_flag = 1;
      if (1){
	prep_vpar_vfile_vl(mylogf,sppl,0,&nvl,&vlf_name,&nlinkval,&vlf_val);
      }else{
	/// OLD WAY - REMOVE ***
	prep_vpar_vfile(mylogf,sppl,rptflag,1,rvval,rvname,rvtype,rvn,rn,rval,
			rvcnt,rnvl,rnvp);
	return;
      }
    }
  }


  nprod = 1;
  nv = 0;
  pp = paramfile_get_first_prefix_pointer(sppl,"VAR_");
  while(pp != NULL){
    strcpy(tstr,&(pp->name[4])); // parname starts at 4th char
    n = paramfile_count_values_pointer(pp);
    sprintf(ggstr,"    `%s' has %d values\n",tstr,n);
    mylog(mylogf,ggstr);
    nprod *= n;
    pp = paramfile_get_next_prefix_pointer(pp,"VAR_");
    nv += 1;
    if (!paramfile_test_param(sppl,tstr))
      mylog_exit(mylogf,"PREP_VAR_PARAMS  var param has no default value\n");
  }

  // *** WYETH NEW 2016 Jan
  // *** WYETH NEW     VARSEQ_[pname]  range_inc  <v1> <v2> <incr>
  // *** WYETH NEW     VARSEQ_[pname]  range_n    <v1> <v2> <n>
  // *** WYETH NEW     VARSEQ_[pname]  init_inc_n <v1> <v2> <incr>
  nvq = 0;
  pp = paramfile_get_first_prefix_pointer(sppl,"VARSEQ_");
  while(pp != NULL){
    strcpy(tstr,&(pp->name[7])); // parname starts at 7th char
    n = paramfile_count_values_pointer(pp);
    if (n != 4)
      mylog_exit(mylogf,"PREP_VAR_PARAMS  VARSEQ needs 4 values\n");
    n = stimu_varseq_get_n(mylogf,pp);  // Get number of values in sequence
    sprintf(ggstr,"    `%s' has %d values\n",tstr,n);
    mylog(mylogf,ggstr);
    nprod *= n;
    pp = paramfile_get_next_prefix_pointer(pp,"VARSEQ_");
    nvq += 1;
    if (!paramfile_test_param(sppl,tstr))
      mylog_exit(mylogf,"PREP_VAR_PARAMS  varseq param has no default val\n");
  }

  nvg = 0;
  pp = paramfile_get_first_prefix_pointer(sppl,"VARGEN_");
  while(pp != NULL){
    strcpy(tstr,&(pp->name[7])); // parname starts at 7th char
    n = paramfile_count_values_pointer(pp);
    if (n!=6)
      mylog_exit(mylogf,"PREP_VAR_PARAMS  VARGEN needs 6 values\n");
    n = paramfile_get_nth_int_param_ptr_or_exit(pp,2);
    nprod *= n;
    pp = paramfile_get_next_prefix_pointer(pp,"VARGEN_");
    nvg += 1;
    if (!paramfile_test_param(sppl,tstr))
      mylog_exit(mylogf,"PREP_VAR_PARAMS  vargen param has no default val\n");
  }

  nvp = 0;
  pp = paramfile_get_first_prefix_pointer(sppl,"VARGENPAIR_");
  while(pp != NULL){
    strcpy(tstr,&(pp->name[11]));  // parname1 starts at 11th character
    n = paramfile_count_values_pointer(pp);
    if (n < 3)
      mylog_exit(mylogf,"PREP_VAR_PARAMS  VARGENPAIR needs >= 3 values\n");
    n = paramfile_get_nth_int_param_ptr_or_exit(pp,1);
    nprod *= n;  // Each pair counts as 1 dimension
    par2 = paramfile_get_nth_char_param_pp_or_exit(pp,0);
    pp = paramfile_get_next_prefix_pointer(pp,"VARGENPAIR_");
    nvp += 1;  // Each will later count as two
    if (!paramfile_test_param(sppl,tstr))
      mylog_exit(mylogf,"PREP_VAR_PARAMS  vargenpair param1 has no default\n");
    if (!paramfile_test_param(sppl,par2))
      mylog_exit(mylogf,"PREP_VAR_PARAMS  vargenpair param2 has no default\n");
  }


  //  If a VARFILE was found above, then counts are already set, otherwise set
  //    nvl - Number of var linked variables
  //    nlinkval - Number of values for link variables (one value for all)
  //
  if (vlf_flag == 0){
    nvl = 0;
    nlinkval = -1;
    pp = paramfile_get_first_prefix_pointer(sppl,"VARLINK_");
    while(pp != NULL){
      strcpy(tstr,&(pp->name[8]));
      if (nlinkval == -1){
	nlinkval = paramfile_count_values_pointer(pp);
	sprintf(ggstr,"    Linked variables have %d values\n",nlinkval);
	mylog(mylogf,ggstr);
      }else{
	n = paramfile_count_values_pointer(pp);
	if (n != nlinkval){
	  printf("  *** %s has %d values\n",tstr,n);
	  mylog_exit(mylogf,
	     "PREP_VAR_PARAMS  All varlink params need same no. of values\n");
	}
      }
      pp = paramfile_get_next_prefix_pointer(pp,"VARLINK_");
      nvl += 1;
      if (!paramfile_test_param(sppl,tstr))
	mylog_exit(mylogf,
		   "PREP_VAR_PARAMS  varlink param has no default val\n");
    }
  }

  if (nlinkval > -1)
    nprod *= nlinkval;

  if (nvl == 1){
    mylogx(mylogf,"PREP_VAR_PARAMS",
	   "The .stm file has a single VARLINK_, please use VAR_\n");
  }

  //
  //  WYETH - NEW - VARSINGLE_
  //
  //  *** Aug 1, 2014 ********** PROBLEMS IF the varsingle param is
  //       also a VARLINK or VAR param - it gets stored twice.
  //  *** MUST FIX
  //  *** MUST FIX
  //  *** MUST FIX
  //
  nvs = 0;
  pp = paramfile_get_first_prefix_pointer(sppl,"VARSINGLE_");
  if (pp != NULL){
    nvs = 1;
    pp = paramfile_get_next_prefix_pointer(pp,"VARSINGLE_");
    if (pp != NULL){
      mylogx(mylogf,"PREP_VAR_PARAMS",
	     "The .stm file can have only one VARSINGLE_\n");
    }
  }

  //nvar = nv + nvg + nvp*2 + nvl + nvs;
  nvar = nv + nvg + nvp*2 + nvq + nvl + nvs;

  if (rptflag == 1)
    nrpt = paramfile_get_int_param_or_exit(sppl,"stim_nrpt");
  else
    nrpt = 1; // Added for labradoodle

  if (nvar == 0){
    mylog(mylogf,"    No variable parameters.\n");
    tvval = NULL;
    vname = NULL;
    vcnt = NULL;
    val = NULL;
    nvar = 0;
    tvn = nrpt; // Total number of trials
  }else{
    sprintf(ggstr,"    %d var params, %d unique trials\n",nvar,nprod);
    mylog(mylogf,ggstr);
    vcnt = get_zero_iarray(nvar);
    vname = (char **)myalloc(nvar*sizeof(char *));
    val = (char ***)myalloc(nvar*sizeof(char **));

    for(i=0;i<nv;i++){
      if (i==0)
	pp = paramfile_get_first_prefix_pointer(sppl,"VAR_");
      else
	pp = paramfile_get_next_prefix_pointer(pp,"VAR_");
      vname[i] = strdup(&(pp->name[4]));
      vcnt[i] = paramfile_get_slist_pointer(sppl,pp,&(val[i]));
    }

    k = nv;
    for(i=0;i<nvq;i++){  // VARSEQ_
      if (i==0)
	pp = paramfile_get_first_prefix_pointer(sppl,"VARSEQ_");
      else
	pp = paramfile_get_next_prefix_pointer(pp,"VARSEQ_");
      vname[k] = strdup(&(pp->name[7]));
      vcnt[k] = stimu_varseq_get_list(mylogf,pp,&(val[k]));
      k += 1;
    }

    for(i=0;i<nvg;i++){  // VARGEN_
      if (i==0)
	pp = paramfile_get_first_prefix_pointer(sppl,"VARGEN_");
      else
	pp = paramfile_get_next_prefix_pointer(pp,"VARGEN_");
      vname[k] = strdup(&(pp->name[7]));
      vcnt[k] = stim_util_get_vargen_list(mylogf,pp,&(val[k]));
      k += 1;
    }

    // Put 'vargenpair' params 2nd to last in list
    for(i=0;i<nvp;i++){
      if (i==0)
	pp = paramfile_get_first_prefix_pointer(sppl,"VARGENPAIR_");
      else
	pp = paramfile_get_next_prefix_pointer(pp,"VARGENPAIR_");

      par2 = paramfile_get_nth_char_param_pp_or_exit(pp,0);

      vname[k] = strdup(&(pp->name[11]));
      vname[k+1] = strdup(par2);
      vcnt[k] = stim_util_get_vargenpair_list(mylogf,pp,&(val[k]),&(val[k+1]));
      vcnt[k+1] = vcnt[k];
      k += 2;
    }

    // Put 'varlink' params (almost) last in list
    if (vlf_flag == 1){
      //  Values came from a VARFILE above.
      for(i=0;i<nvl;i++){
	vname[k] = strdup(vlf_name[i]);
	vcnt[k]  = nlinkval;
	val[k]   = vlf_val[i]; // WYETH - Not sure about freeing this later??
	k += 1;
      }
    }else{
      for(i=0;i<nvl;i++){
	if (i==0)
	  pp = paramfile_get_first_prefix_pointer(sppl,"VARLINK_");
	else
	  pp = paramfile_get_next_prefix_pointer(pp,"VARLINK_");
	vname[k] = strdup(&(pp->name[8]));
	vcnt[k] = paramfile_get_slist_pointer(sppl,pp,&(val[k]));

	k += 1;
      }
    }

    // Put 'varsingle' param VERY last in list
    if (nvs == 1){
      pp = paramfile_get_first_prefix_pointer(sppl,"VARSINGLE_");
      vname[k] = strdup(&(pp->name[10]));
      vcnt[k] = 2;  // Takes only two values, the default and VARSINGLE_ one
      val[k] = (char **)myalloc(2*sizeof(char *));

      // Put the VARSINGLE_ value second, at index 1
      val[k][1] = carray_item_list_remove_comment(pp->value);

      pp = paramfile_get_param_pointer(sppl,vname[k]);
      if (pp == NULL){
	mylogx(mylogf,"PREP_VAR_PARAMS",
	       "Cannot find default value for VARSINGLE_\n");
      }
      // Put the default value first, at index 0
      val[k][0] = carray_item_list_remove_comment(pp->value);

      //printf("HERE WYETH___________  val[%d][0] = %s\n",k,val[k][0]);
      //printf("HERE WYETH___________  val[%d][1] = %s\n",k,val[k][1]);
      //exit(0);

      k += 1;
    }

    if (nrpt == 1)
      sprintf(ggstr,"    %d repeat of each trial\n",nrpt);
    else
      sprintf(ggstr,"    %d repeats of each trial\n",nrpt);
    mylog(mylogf,ggstr);

    // This fills in all var param values EXCEPT those for 'VARSINGLE_',
    // which are not part of the 'product' matrix.
    get_cross_product_params(mylogf,nvar,vcnt,val,nrpt,nvp,nvl,nlinkval,nvs,
			     &tvn,&tvval);
    //
    //  Fill in any 'VARSINGLE_' stimuli using default values for all
    //  but the last (varsingle) value.
    //
    if (nvs == 1){ // There are 'nrpt' blanks at the end to be filled
      for(j=0;j<(nvar-1);j++){
	pp = paramfile_get_param_pointer(sppl,vname[j]);
	for(i=0;i<nrpt;i++){
	  tvval[tvn-nrpt + i][j] = carray_item_list_remove_comment(pp->value);
	}
      }
      j = nvar - 1;  // The VARSINGLE_ name
      for(i=0;i<nrpt;i++){
	tvval[tvn-nrpt + i][j] = strdup(val[nvar-1][1]);
      }
    }

    sprintf(ggstr,"    %d total stimuli\n",tvn);
    mylog(mylogf,ggstr);
  }

  /*
  printf("----------------------------------------------------------------\n");
  for(i=0;i<tvn;i++){
    printf("%d",i);
    for(j=0;j<nvar;j++){
      printf(" %s",tvval[i][j]);
    }
    printf("\n");
  }
  printf("----------------------------------------------------------------\n");
  exit(0);
*/


  // Set types for 'var' params
  vtype = (char *)myalloc(nvar*sizeof(char));
  for(i=0;i<nvar;i++){
    vtype[i] = paramfile_get_type_from_comment(sppl,vname[i],'f');
  }

  //
  //  For 'frameset' stimuli, set filename params to be 'c'
  //
  stimtype = param_getc_exit(sppl,"stim_type");
  if (strcmp(stimtype,"frameset")==0){
    k = search_2d_carray(vname,"fst_file_1",nvar);
    if (k >= 0)
      vtype[k] = 'c';
    k = search_2d_carray(vname,"fst_file_2",nvar);
    if (k >= 0)
      vtype[k] = 'c';
  }

  *rvval = tvval; *rvname = vname; *rvn = nvar; *rn = tvn; *rvtype = vtype;
  *rval = val; *rvcnt = vcnt; *rnvl = nvl;  // Added for labradoodle
  *rnvp = nvp;
}
/**************************************-**************************************/
/*                                                                           */
/*                            STIMU_VAR_PARAM_CHECK                          */
/*                                                                           */
/*  If 'pname' is a var param, return the number of values and the list of   */
/*  values that it takes.                                                    */
/*                                                                           */
/*****************************************************************************/
void stimu_var_param_check(sppl,pname,rn,rval)
     struct param_pair_list *sppl;   // Stimulus parameter list
     char *pname;                    // param name to check
     int *rn;                        // number of unique values
     char ***rval;                   // Unique values [*rn]
{
  int k;
  int n,vn,*vcnt,nvl,nvp;
  char **vname,***vval,***val,*vtype;

  //printf("  STIMU_VAR_PARAM_CHECK\n");

  prep_var_params(NULL,sppl,0,&vval,&vname,&vtype,&vn,&n,&val,
		  &vcnt,&nvl,&nvp);

  k = search_2d_carray(vname,pname,vn);
  if (k < 0){
    *rn = 0;  // No variable params
    *rval = NULL;
  }else{
    *rn = vcnt[k];
    *rval = copy_2d_carray(val[k],vcnt[k]);
  }

  free_var_params(NULL,vval,vname,vn,n,val,vcnt);
  myfree(vtype);

}
/**************************************-**************************************/
/*                                                                           */
/*                              PREP_CONST_PARAMS                            */
/*                                                                           */
/*  Determine the names, numbers, types, and values for const params in a    */
/*  format that is useful for writing ndata output.                          */
/*                                                                           */
/*****************************************************************************/
void prep_const_params(mylogf,s,moo_par,moo_val,nmoo,rcval,rcname,rctype,rncon)
     char *mylogf;
     struct stim_struct *s;  // Stimulus params
     char **moo_par;         // Name [nmoo] model params from command line
     char **moo_val;         // Value [nmoo] model params from command line
     int nmoo;               // Number of model  params from command line
     char ***rcval;          // const param values [rncon]
     char ***rcname;         // const param names [rncon]
     char **rctype;          // const param types [rncon]
     int *rncon;             // number of const params
{
  int i;
  int ncon,nvar;
  char **cname,**cval,*ctype,*vtype,tstr[SLEN];
  struct param_pair *pp;

  mylog(mylogf,"  PREP_CONST_PARAMS\n");

  // Everything in 's->ppl' that isn't in 's->vname' is const
  ncon = nmoo;  // Count includes command line model pars
  pp = s->ppl->p;
  while(pp != NULL){
    if (search_2d_carray(s->vname,pp->name,s->nvar) == -1){
      ncon += 1;
    }
    pp = pp->next;
  }
  //printf("ncon = %d\n",ncon);

  cname = (char **)myalloc(ncon*sizeof(char *));
  ctype = (char  *)myalloc(ncon*sizeof(char));
  cval  = (char **)myalloc(ncon*sizeof(char *));

  for(i=0;i<nmoo;i++){
    sprintf(tstr,"MOO_%s",moo_par[i]);
    cname[i] = strdup(tstr);
    ctype[i] = 'c';
    cval[i]  = strdup(moo_val[i]);
  }

  i = nmoo;
  pp = s->ppl->p;
  while(pp != NULL){
    if (search_2d_carray(s->vname,pp->name,s->nvar) == -1){
      cname[i] = strdup(pp->name);
      //cval[i] = strdup(pp->value);  WYETH, CHANGED Apr 2007
      cval[i] = carray_item_list_remove_comment(pp->value);

      // Use 'f' as the default type for all parameters
      ctype[i] = paramfile_get_type_from_comment(s->ppl,pp->name,'f');

      i += 1;
    }
    pp = pp->next;
  }

  *rcval = cval; *rcname = cname; *rctype = ctype; *rncon = ncon;
}
/**************************************-**************************************/
/*                                                                           */
/*                             UPDATE_CONST_PARAM                            */
/*                                                                           */
/*  Update the value and/or type of the const param having 'tname'.          */
/*                                                                           */
/*****************************************************************************/
void update_const_param(mylogf,tname,tval,ttype,cname,cval,ctype,n)
     char *mylogf;      // log file
     char *tname;       // Look up this name
     char *tval;        // new value, or NULL
     char ttype;        // new type, or '\0'
     char **cname;      // const param names [n]
     char **cval;       // const param values [n]
     char *ctype;       // const param types [n]
     int n;             // number of const params
{
  int k;
  char *ggstr;

  mylog(mylogf,"  UPDATE_CONST_PARAM\n");

  // Everything in 's->ppl' that isn't in 's->vname' is const
  k = search_2d_carray(cname,tname,n);
  if (k == -1){
    sprintf(ggstr,"Name:  %s\n",tname);
    mylog(mylogf,ggstr);
    mylogx(mylogf,"UPDATE_CONST_PARAM","Cannot find parameter name\n");
  }

  if (tval != NULL){
    if (cval[k] != NULL)
      myfree(cval[k]);
    cval[k] = strdup(tval);
    //printf("UPDATE_CONST_PARAM  cval = %s\n",tval);
  }
  if (ttype != '\0'){
    ctype[k] = ttype;
    //printf("UPDATE_CONST_PARAM  ctype = %c\n",ttype);
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_DOT_OPP_MO                              */
/*                                                                           */
/*  Get a pair of movies of dots moving in opposite directions, crossing in  */
/*  the middle frame of the movie.                                           */
/*                                                                           */
/*  For Sanada and DeAngelis (2014) stimulus.                                */
/*                                                                           */
/*****************************************************************************/
void get_dot_opp_mo(sscale,framerate,speed,disp0,theta,sizex,sizey,nframes,dpf,
		    seed,corr,rdotx1,rdoty1,rdotx2,rdoty2,rdotn,rampl)
     float sscale;      // (deg/pix)
     float framerate;   // (frames/s)
     float speed;       // (deg/s)
     float disp0;       // Disparity at middle frame (deg)
     float theta;       // direction (deg)
     float sizex,sizey; // Size of dot field (deg)
     int nframes;       // Length of movie (frames)
     int dpf;           // Dots per frame
     int seed;          // Randomization seed for dot positions.
     int corr;          // 1-correlated, 0-uncorrelated dots across eyes
     float ***rdotx1,***rdoty1,***rdotx2,***rdoty2;  // *[nframes][dpf]
     int **rdotn;                                    // *[nframes]
     float ***rampl;                                 // *[nframes][dpf]
{
  int i,j;
  int *dotn;
  float **dotx1,**doty1,**dotx2,**doty2,**dota,x,y,dx0,dy0,tc,a;
  double dx,dy,x1,y1,x2,y2;

  if (seed > 0)
    seed = -seed;

  dotn  = (int *)myalloc(nframes*sizeof(int));
  dotx1 = (float **)myalloc(nframes*sizeof(float *));
  doty1 = (float **)myalloc(nframes*sizeof(float *));
  dotx2 = (float **)myalloc(nframes*sizeof(float *));
  doty2 = (float **)myalloc(nframes*sizeof(float *));
  dota  = (float **)myalloc(nframes*sizeof(float *));
  for(i=0;i<nframes;i++){
    dotn[i]  = dpf;
    dotx1[i] = (float *)myalloc(dpf*sizeof(float));
    doty1[i] = (float *)myalloc(dpf*sizeof(float));
    dotx2[i] = (float *)myalloc(dpf*sizeof(float));
    doty2[i] = (float *)myalloc(dpf*sizeof(float));
    dota[i]  = (float *)myalloc(dpf*sizeof(float));
  }

  // Compute time origin (time at which dots cross paths)
  tc = (float)(nframes-1)/2.0;  // (4-1)/2 = 1.5

  // Change in x and y position (pixels/frame)
  dx = (double)speed / (framerate * sscale) * cos(theta/180.0*M_PI);
  dy = (double)speed / (framerate * sscale) * sin(theta/180.0*M_PI);

  dx0 = disp0 * cos(theta/180.0*M_PI);
  dy0 = disp0 * sin(theta/180.0*M_PI);

  //
  // Create first dot frame
  //
  for(j=0;j<dpf;j++){
    x = dx0 + myrand_util_ran2(&seed)*sizex;  // Select mid-point value
    y = dy0 + myrand_util_ran2(&seed)*sizey;
    if (myrand_util_ran2(&seed) < 0.5)
      a = -1.0;
    else
      a = 1.0;

    x1 = x - dx*tc;  // Compute start of trajectory
    y1 = y - dy*tc;

    if (corr != 1){
      // Pick new positions for dots in other eye
      x = dx0 + myrand_util_ran2(&seed)*sizex;  // Select mid-point value
      y = dy0 + myrand_util_ran2(&seed)*sizey;
    }
    x2 = x + dx*tc;
    y2 = y + dy*tc;
      
    wrap_double(&x1,(double)sizex); // Wrap values within [0,size*)
    wrap_double(&x2,(double)sizex);
    wrap_double(&y1,(double)sizey);
    wrap_double(&y2,(double)sizey);

    for (i=0;i<nframes;i++){  // For all subsequent frames

      dotx1[i][j] = (float)x1;  // Store current position
      doty1[i][j] = (float)y1;
      dotx2[i][j] = (float)x2;
      doty2[i][j] = (float)y2;
      dota[i][j]  = a;          // Store dot amplitude

      x1 += dx;  // Update the dot positions
      y1 += dy;
      x2 -= dx;
      y2 -= dy;

      if (x1 < 0.0) x1 += sizex; else if (x1 >= sizex) x1 -= sizex; // Wrap
      if (x2 < 0.0) x2 += sizex; else if (x2 >= sizex) x2 -= sizex;
      if (y1 < 0.0) y1 += sizey; else if (y1 >= sizey) y1 -= sizey;
      if (y2 < 0.0) y2 += sizey; else if (y2 >= sizey) y2 -= sizey;
    }
  }
  *rdotx1 = dotx1;
  *rdoty1 = doty1;
  *rdotx2 = dotx2;
  *rdoty2 = doty2;
  *rampl  = dota;
  *rdotn  = dotn;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_COH_DOT_MOVIE_AMPL                          */
/*                                                                           */
/*  Like 'GET_COHERENCE_DOT_MOVIE' but also returns amplitude values, which  */
/*  are randomly assigned to +/- 1.0.                                        */
/*                                                                           */
/*****************************************************************************/
void get_coh_dot_movie_ampl(sscale,framerate,coher,speed,theta,sizex,sizey,
			    nframes,dpf,dt,seed,rdotx,rdoty,rdotn,rampl)
     float sscale,framerate,coher,speed,theta,sizex,sizey;
     int nframes,dpf,dt,seed;
     float ***rdotx,***rdoty;
     int **rdotn;
     float ***rampl;
{
  int i,j,k;
  int *dotn,n;
  float **dotx,**doty,x,y,dx,dy,**dota,a;

  if (seed > 0)
    seed = -seed;

  dotn = (int *)myalloc(nframes*sizeof(int));
  dotx = (float **)myalloc(nframes*sizeof(float *));
  doty = (float **)myalloc(nframes*sizeof(float *));
  dota = (float **)myalloc(nframes*sizeof(float *));
  for(i=0;i<nframes;i++){
    dotn[i] = dpf;
    dotx[i] = (float *)myalloc(dpf*sizeof(float));
    doty[i] = (float *)myalloc(dpf*sizeof(float));
    dota[i] = (float *)myalloc(dpf*sizeof(float));
  }

  dy = speed * (float)dt/framerate / sscale * sin(theta/180.0*M_PI);
  dx = speed * (float)dt/framerate / sscale * cos(theta/180.0*M_PI);

  for (i=0;i<nframes;i++){
    n = 0; // number of dots created so far for the ith frame
    k = i - dt; // k is index to dots in a previous frame
    if (k >= 0){ // START WITH DOTS FROM EARLIER FRAME, if that frame exists
      for (j=0;j<dpf;j++) // for each dot
	if (myrand_util_ran2(&seed) < coher){ // REPLOT AS SIGNAL DOT
	  x = dotx[k][j] + dx;
	  y = doty[k][j] + dy;
	  // Take care of wrap-around
	  if (x < 0.0) x += sizex; else if (x >= sizex) x -= sizex;
	  if (y < 0.0) y += sizey; else if (y >= sizey) y -= sizey;
	  dotx[i][n] = x; // store the coordinates
	  doty[i][n] = y;
	  dota[i][n] = dota[k][j];  // Carry amplitude of dot
	  n += 1;
	}
    }
    while (n < dpf){ // Choose additional, randomly placed, dots
      dotx[i][n] = myrand_util_ran2(&seed)*sizex;
      doty[i][n] = myrand_util_ran2(&seed)*sizey;
      if (myrand_util_ran2(&seed) < 0.5)
	dota[i][n] = -1.0;
      else
	dota[i][n] = 1.0;
      n += 1;
    }
  }
  *rdotn = dotn; *rdotx = dotx; *rdoty = doty; *rampl = dota;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_COHERENCE_DOT_MOVIE                          */
/*                                                                           */
/*  PARAMETERS:                                                              */
/*                                                                           */
/*  sscale       deg/pix                                                     */
/*  framerate    Video rate (frames/s)                                       */
/*  coher        Motion coherence, from 0 to 1  (Britten et al 1992).        */
/*  speed        deg/s                                                       */
/*  theta        Angle in degrees.                                           */
/*  sizex        Size of rectangular frame for making/wrapping dots (pix)    */
/*  sizey        Size of rectangular frame for making/wrapping dots (pix)    */
/*  nframes      Time duration of stimulus (video frames).                   */
/*  dpf          Dots per frame.                                             */
/*  dt           Frames between refresh of signal dots (AKA ref_frames).     */
/*  seed         Address of randomization seed.                              */
/*                                                                           */
/*  Dot coordinates are returned in pixels (float) within a rectangle        */
/*  of 'sizex' by 'sizey' pixels.                                            */
/*                                                                           */
/*****************************************************************************/
void get_coherence_dot_movie(sscale,framerate,coher,speed,theta,sizex,sizey,
			     nframes,dpf,dt,seed,rdotx,rdoty,rdotn)
     float sscale,framerate,coher,speed,theta,sizex,sizey;
     int nframes,dpf,dt,seed;
     float ***rdotx,***rdoty;
     int **rdotn;
{
  int i,j,k;
  int *dotn,n;
  float **dotx,**doty,x,y,dx,dy;

  if (seed > 0)
    seed = -seed;

  dotn = (int *)myalloc(nframes*sizeof(int));
  dotx = (float **)myalloc(nframes*sizeof(float *));
  doty = (float **)myalloc(nframes*sizeof(float *));
  for(i=0;i<nframes;i++){
    dotn[i] = dpf;
    dotx[i] = (float *)myalloc(dpf*sizeof(float));
    doty[i] = (float *)myalloc(dpf*sizeof(float));
  }

  dy = speed * (float)dt/framerate / sscale * sin(theta/180.0*M_PI);
  dx = speed * (float)dt/framerate / sscale * cos(theta/180.0*M_PI);

  for (i=0;i<nframes;i++){
    n = 0; // number of dots created so far for the ith frame
    k = i - dt; // k is index to dots in a previous frame
    if (k >= 0){ // START WITH DOTS FROM EARLIER FRAME, if that frame exists
      for (j=0;j<dpf;j++) // for each dot
	if (myrand_util_ran2(&seed) < coher){ // REPLOT AS SIGNAL DOT
	  x = dotx[k][j] + dx;
	  y = doty[k][j] + dy;
	  // Take care of wrap-around
	  if (x < 0.0) x += sizex; else if (x >= sizex) x -= sizex;
	  if (y < 0.0) y += sizey; else if (y >= sizey) y -= sizey;
	  dotx[i][n] = x; // store the coordinates
	  doty[i][n] = y;
	  n += 1;
	}
    }
    while (n < dpf){ // Choose additional, randomly placed, dots
      dotx[i][n] = myrand_util_ran2(&seed)*sizex;
      doty[i][n] = myrand_util_ran2(&seed)*sizey;
      n += 1;
    }
  }
  *rdotn = dotn; *rdotx = dotx; *rdoty = doty;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_PAIRABLE_DOT_MOVIE                          */
/*                                                                           */
/*  Implement the paired dots of Qian et al. (1994).                         */
/*                                                                           */
/*  PARAMETERS:                                                              */
/*                                                                           */
/*  sscale       deg/pix                                                     */
/*  framerate    Video rate (frames/s)                                       */
/*  dist         Total distance during lifetime (deg)                        */
/*  offx         Spatial x-offset (deg)                                      */
/*  offy         Spatial y-offset (deg)                                      */
/*  theta        Angle in degrees.                                           */
/*  sizex        Size of rectangular frame for making/wrapping dots (pix)    */
/*  sizey        Size of rectangular frame for making/wrapping dots (pix)    */
/*  lifet        Lifetime of dots (frames)                                   */
/*  nframes      Time duration of stimulus (video frames).                   */
/*  dpf          Dots per frame.                                             */
/*  dt           Frames between refresh of signal dots (AKA ref_frames).     */
/*  seed         Address of randomization seed.                              */
/*                                                                           */
/*  Dot coordinates are returned in pixels (float) within a rectangle        */
/*  of 'sizex' by 'sizey' pixels.                                            */
/*                                                                           */
/*****************************************************************************/
void get_pairable_dot_movie(sscale,framerate,dist,offx,offy,theta,sizex,sizey,
			    lifet,nframes,dpf,seed,rdotx,rdoty,rdotn)
     float sscale,framerate,dist,offx,offy,theta,sizex,sizey;
     int lifet,nframes,dpf,seed;
     float ***rdotx,***rdoty;
     int **rdotn;
{
  int i,j;
  int *dotn,*age;
  float **dotx,**doty,x,y,dx,dy,th,pixpfr,tdist,af,xoff,yoff,soffx,soffy;

  if (seed > 0)
    seed = -seed;

  dotn = (int    *)myalloc(nframes*sizeof(int));
  dotx = (float **)myalloc(nframes*sizeof(float *));
  doty = (float **)myalloc(nframes*sizeof(float *));
  for(i=0;i<nframes;i++){
    dotn[i] = dpf;
    dotx[i] = (float *)myalloc(dpf*sizeof(float));
    doty[i] = (float *)myalloc(dpf*sizeof(float));
  }

  age = (int *)myalloc(dpf*sizeof(int));

  th = theta /180.0*M_PI;   // convert theta to radians

  tdist = dist / sscale;  // Total distance in pix
  pixpfr = tdist / (float)(lifet - 1);  // pixel movement per frame

  dx = pixpfr * cos(th);  // step size in pixels
  dy = pixpfr * sin(th);

  xoff = -tdist/2.0 * cos(th);
  yoff = -tdist/2.0 * sin(th);

  soffx = offx / sscale;  // Offset (pix)
  soffy = offy / sscale;

  //
  //  Generate 'dpf' dots at random positions and with random ages.
  //
  for(j=0;j<dpf;j++){ // for each dot

    age[j] = (int)((float)lifet * myrand_util_ran2(&seed));
    af = (float)age[j];

    //  AgeAdj   Random_Reference_Location     Motion  Offset
    x = af*dx + myrand_util_ran2(&seed)*sizex + xoff + soffx;
    y = af*dy + myrand_util_ran2(&seed)*sizey + yoff + soffy;

    //  Handle wrap-around
    if (x < 0.0) x += sizex; else if (x >= sizex) x -= sizex;
    if (y < 0.0) y += sizey; else if (y >= sizey) y -= sizey;
    dotx[0][j] = x;
    doty[0][j] = y;
  }

  for(i=1;i<nframes;i++){
    for(j=0;j<dpf;j++){ // for each dot
      age[j] += 1;
      if (age[j] < lifet){
	x = dotx[i-1][j] + dx;
	y = doty[i-1][j] + dy;
	if (x < 0.0) x += sizex; else if (x >= sizex) x -= sizex;
	if (y < 0.0) y += sizey; else if (y >= sizey) y -= sizey;
	dotx[i][j] = x;
	doty[i][j] = y;

      }else{
	//  Create a new dot somewhere else
	age[j] = 0;

	//   Random_Reference_Location     Motion  Offset
	x = myrand_util_ran2(&seed)*sizex + xoff + soffx;
	y = myrand_util_ran2(&seed)*sizey + yoff + soffy;

	if (x < 0.0) x += sizex; else if (x >= sizex) x -= sizex;
	if (y < 0.0) y += sizey; else if (y >= sizey) y -= sizey;
	dotx[i][j] = x;
	doty[i][j] = y;
      }
    }
  }
  myfree(age);

  *rdotn = dotn; *rdotx = dotx; *rdoty = doty;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MASK_DOT_MOVIE                              */
/*                                                                           */
/*****************************************************************************/
void mask_dot_movie(dotx,doty,dota,dotn,nframes,x0,y0,cx,cy,w,h,sscale,aptype)
     float **dotx,**doty;    // Dot coords (pixels)
     float **dota;           // Dot amplitude, or NULL
     int *dotn,nframes;
     float x0,y0;            // Origin w/i dot coords (pixels)
     float cx,cy,w,h;        // Aperture center, width, and height (deg)
     float sscale;           // Spatial scale (deg/pix)
     int aptype;
{
  int i,j;
  int n,dflag;
  float x,y,r2pix,tx,ty,apcx,apcy,wo2,ho2;

  // For rectangular aperture
  wo2 = w/(2.0*sscale);   // w/2 (pixels)
  ho2 = h/(2.0*sscale);   // h/2 (pixels)

  // For circular aperture
  r2pix = wo2*wo2;        // pixel radius squared

  // Aperture center in pixel coords
  apcx = x0 + cx/sscale;
  apcy = y0 + cy/sscale;

  for (i=0;i<nframes;i++){
    n = 0; // number of dots kept so far in the ith frame
    for (j=0;j<dotn[i];j++){
      tx = dotx[i][j];
      ty = doty[i][j];
      x = tx - apcx;      // Coords relative to ap center (pixels)
      y = ty - apcy;
      /*y *= -1;*/ /*** INVERT y-coordinate ***/

      dflag = 0;
      if (aptype == 1){ // circular
	if ((x*x + y*y) <= r2pix)
	  dflag = 1;
      }else if (aptype == 2){ // rectangular
	if ((x >= -wo2) && (x <= wo2) && (y >= -ho2) && (y <= ho2))
	  dflag = 1;
      }else{
	printf("*** Aperture type %d not implemented yet\n",aptype);
	dflag = 1;
      }

      if (dflag){
	dotx[i][n] = tx; // Over-write old coordinates
	doty[i][n] = ty;
	if (dota != NULL)
	  dota[i][n] = dota[i][j];
	n += 1;
      }
    }
    dotn[i] = n;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            TRANSLATE_DOT_MOVIE                            */
/*                                                                           */
/*****************************************************************************/
void translate_dot_movie(dotx,doty,dotn,nframes,xoff,yoff)
     float **dotx,**doty;    // Dot coords (pixels)
     int *dotn,nframes;
     float xoff,yoff;        // Add these offsets (pix)
{
  int i,j;

  for (i=0;i<nframes;i++){
    for (j=0;j<dotn[i];j++){
      dotx[i][j] += xoff;
      doty[i][j] += yoff;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_DOT_STEREOGRAM                            */
/*                                                                           */
/*  PARAMETERS:                                                              */
/*                                                                           */
/*  sscale       deg/pix                                                     */
/*  framerate    Video rate (frames/s)                                       */
/*  coher        Motion coherence, from 0 to 1  (Britten et al 1992).        */
/*  speed        deg/s                                                       */
/*  theta        Angle in degrees.                                           */
/*  sizex        Size of rectangular frame for making/wrapping dots (pix)    */
/*  sizey        Size of rectangular frame for making/wrapping dots (pix)    */
/*  nframes      Time duration of stimulus (video frames).                   */
/*  dpf          Dots per frame.                                             */
/*  dt           Frames between refresh of signal dots (AKA ref_frames).     */
/*  seed         Address of randomization seed.                              */
/*                                                                           */
/*  Dot coordinates are returned in pixels (float) within a rectangle        */
/*  of 'sizex' by 'sizey' pixels.                                            */
/*                                                                           */
/*****************************************************************************/
void get_dot_stereogram(ptype,sscale,framerate,coher,speed,theta,tf,waveform,
			sizex,sizey,ox1,ox2,oy1,oy2,dx,dy,dx0,dy0,nframes,
			ddens,dt,seed,dotcorr,rdxl,rdyl,rdnl,rdxr,rdyr,rdnr,
			rdampl,rdampr)
     int ptype;         // Patch type;
     float sscale,framerate,coher,speed,theta,tf;
     char *waveform;    // "const", "sine", "ramp"
     float sizex,sizey; // patch width, height (pix)
     float ox1,ox2;     // Patch boundaries along x-axis
     float oy1,oy2;     // Patch boundaries along y-axis
     float dx,dy;       // (pix) to shift patch (end of range, if 'tf' !=0)
     float dx0,dy0;     // Start of range for shifting (if 'tf' !=0)
     int nframes;
     float ddens;       // Dots per deg squared
     int dt,seed;
     int dotcorr;       // 0: uncorr. spatial positions across eyes
                        // -1: opposite sign - NOT HANDLED HERE.
     float ***rdxl,***rdyl;
     int **rdnl;
     float ***rdxr,***rdyr;
     int **rdnr;
     float ***rdampl,***rdampr;  // *[nframes][rdn_]
{
  int i;
  int fi,n,tn,seed_f,dpf,*dnl,*dnr,*dnf;
  float **dxl,**dyl,**dxr,**dyr,**dxf,**dyf,*tx,*ty,*ta;
  float x,y,a,vdx,vdy,tsec,tcyc,nx1,nx2,ny1,ny2;
  float **dampl,**dampr,**dampf,xa,xb,ya,yb,oszx,oszy;

  //printf("  Patch (x0,y0) = %f,%f\n",px0,py0);
  //printf("        w,h  = %f %f\n",pw,ph);
  //printf("  Size x,y   = %f %f\n",sizex,sizey);
  //printf("    dx0,dx   = %f %f\n",dx0,dx);
  //printf("    dy0,dy   = %f %f\n",dy0,dy);

  //
  //  If disparity is not constant, set min and range to vary over
  //
  if (tf > 0.0){
    // Used for sinusoidal waveform
    xa = dx0 + (dx - dx0)/2.0;  // Mean disparity
    xb = (dx - dx0)/2.0;        // Amplitude of variation
    ya = dy0 + (dy - dy0)/2.0;
    yb = (dy - dy0)/2.0;
    //printf("xa = %f  xb = %f  (dx = %f  dx0 = %f)\n",xa,xb,dx,dx0);
    //printf("ya = %f  yb = %f\n",ya,yb);
  }
  //
  //  WYETH - THESE PROBABLY NEED TO DEPEND ON TF, because dx and dy are
  //          then changing with time ***  2013 Sep 16
  //
  /*
  nx1 = ox1 + dx;
  nx2 = ox2 + dx;
  ny1 = oy1 + dy;
  ny2 = oy2 + dy;
  */

  seed_f = 7 * seed + 1013;  // For picking fill dots

  if (ptype == 0){  // Full field patch
    //
    // WYETH - New way May 2016: create one movie that is over-sized,
    //         copy it to the R.E. and shift the dots in the R.E. and/or
    //         the L.E.
    //

    oszx = my_rint(sizex + fabs(dx) + 1.0);  // Larger size to account for shift
    oszy = my_rint(sizey + fabs(dy) + 1.0);  // Larger size to account for shift
    
    // Dots per frame, within oversized frame
    dpf = (int)(oszx*oszy * sscale*sscale * ddens); // Area X density

    get_coh_dot_movie_ampl(sscale,framerate,coher,speed,theta,oszx,oszy,
			   nframes,dpf,dt,seed,&dxl,&dyl,&dnl,&dampl);

    stm_dot_movie_copy(dxl,dyl,dampl,dnl,nframes,&dxr,&dyr,&dampr,&dnr);

    if (dx < 0.0)       // Shift dots leftward in right movie
      stm_dot_movie_shift(dxr,dyr,dnr,nframes,dx,0.0);
    else if (dx > 0.0)  // Shift dots lefward in left movie
      stm_dot_movie_shift(dxl,dyl,dnl,nframes,-dx,0.0);

    if (dy < 0.0)       // Shift dots down in right movie
      stm_dot_movie_shift(dxr,dyr,dnr,nframes,0.0,dy);
    else if (dx > 0.0)  // Shift dots down in left movie
      stm_dot_movie_shift(dxl,dyl,dnl,nframes,0.0,-dy);


    *rdxl = dxl;  *rdyl = dyl;  *rdnl = dnl;
    *rdxr = dxr;  *rdyr = dyr;  *rdnr = dnr;
    *rdampl = dampl;
    *rdampr = dampr;
    return;
  }

  dpf = (int)(sizex*sizey * sscale*sscale * ddens); // Area X density

  // Get left dot movie - used for main stimulus (ampl +/- 1)
  get_coh_dot_movie_ampl(sscale,framerate,coher,speed,theta,sizex,sizey,
			 nframes,dpf,dt,seed,&dxl,&dyl,&dnl,&dampl);

  // Get movie for filling in
  get_coh_dot_movie_ampl(sscale,framerate,coher,speed,theta,sizex,sizey,
			 nframes,dpf,dt,seed_f,&dxf,&dyf,&dnf,&dampf);

  if (dotcorr == 0){  // Make L and R be completely uncorrelated
    dxr = dxf;
    dyr = dyf;
    dnr = dnf;
    dampr = dampf;
  }else{

    // Space for right movie
    dxr   = (float **)myalloc(nframes*sizeof(float *));
    dyr   = (float **)myalloc(nframes*sizeof(float *));
    dnr   =    (int *)myalloc(nframes*sizeof(int));
    dampr = (float **)myalloc(nframes*sizeof(float *));

    // Temporary dot storage
    tx = (float *)myalloc(dpf*2 * sizeof(float));
    ty = (float *)myalloc(dpf*2 * sizeof(float));
    ta = (float *)myalloc(dpf*2 * sizeof(float));

    for(fi=0;fi<nframes;fi++){  // frame index

      //
      //  Compute the disparity shift for this frame
      //
      if (tf > 0.0){
	tsec = fi / framerate;
	// WYETH HERE 2013 Sept 16
	//vdx = dx * sin(tf*tsec*2.0*M_PI);
	//vdy = dy * sin(tf*tsec*2.0*M_PI);

	if (strcmp(waveform,"ramp")==0){
	  tcyc = (tsec * tf) - (int)(tsec * tf);

	  vdx = dx0 + (dx - dx0)*tcyc;
	  vdy = dy0 + (dy - dy0)*tcyc;

	  //printf("%d  %f  %f   disp: %f\n",fi,tsec,tcyc,vdx);

	}else if (strcmp(waveform,"sine")==0){
	  vdx = xa + xb*sin(tf*tsec*2.0*M_PI);
	  vdy = ya + yb*sin(tf*tsec*2.0*M_PI);
	}else{
	  printf("  *** patch_waveform:  %s\n",waveform);
	  exit_error("GET_DOT_STEREOGRAM","Unknown 'patch_waveform'");
	}

      }else{
	vdx = dx;
	vdy = dy;
      }

      nx1 = ox1 + vdx;
      nx2 = ox2 + vdx;
      ny1 = oy1 + vdy;
      ny2 = oy2 + vdy;

      //printf("%d  %f\n",fi,vdx);

      // ******
      // ******  WYETH - This method is fine for *NEGATIVE* vdx, for things
      // ******    that appear behind the aperture, but this is not correct
      // ******    for positive vdx, where the patch floats in front.  In the
      // ******    latter case, there are essentially two "clipping regions"
      // ******    which are offset equally and in opposite directions in the
      // ******    two images.
      // ******

      // Copy dots from L to R, shifting some, keeping some, removing some
      n = dnl[fi];  // number of dots in left frame
      tn = 0;       // number of dots in right frame

      for(i=0;i<n;i++){
	x =   dxl[fi][i];
	y =   dyl[fi][i];
	a = dampl[fi][i];
	if ((x >= ox1) && (x <= ox2) && (y >= oy1) && (y <= oy2)){
	  tx[tn] = x + vdx;  // Dot is in the original rectangle, so shift it
	  ty[tn] = y + vdy;
	  ta[tn] = a;
	  if ((tx[tn] >= ox1) && (tx[tn] <= ox2) &&
	      (ty[tn] >= oy1) && (ty[tn] <= oy2))
	    tn += 1;
	}else{
	  tx[tn] = x;  // Dot is background, keep it as is
	  ty[tn] = y;
	  ta[tn] = a;
	  tn += 1;
	}
      }

      // Move filler dots here
      n = dnf[fi];  // number of dots in filler frame
      for(i=0;i<n;i++){
	x =   dxf[fi][i];
	y =   dyf[fi][i];
	a = dampf[fi][i];
	if ((x >= ox1) && (x <= ox2) && (y >= oy1) && (y <= oy2)){
	  // Dot is in the origin rectangle
	  if (!((x >= nx1) && (x <= nx2) && (y >= ny1) && (y <= ny2))){
	    // and not in target, move this in as filler
	    tx[tn] = x;  // Dot is background, keep it as is
	    ty[tn] = y;
	    ta[tn] = a;
	    tn += 1;
	  }
	}
      }

      // Store this frame for right
      dxr[fi]   = copy_farray(tx,tn);
      dyr[fi]   = copy_farray(ty,tn);
      dampr[fi] = copy_farray(ta,tn);
      dnr[fi] = tn;
    }

    myfree(tx); myfree(ty); myfree(ta);
    free_2d_farray(  dxf,nframes);
    free_2d_farray(  dyf,nframes);
    free_2d_farray(dampf,nframes);
    myfree(dnf);
  }

  *rdxl = dxl;  *rdyl = dyl;  *rdnl = dnl;
  *rdxr = dxr;  *rdyr = dyr;  *rdnr = dnr;
  *rdampl = dampl;
  *rdampr = dampr;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_RDS_DISP_MOD                              */
/*                                                                           */
/*  PARAMETERS:                                                              */
/*                                                                           */
/*  sscale       deg/pix                                                     */
/*  framerate    Video rate (frames/s)                                       */
/*  coher        Motion coherence, from 0 to 1  (Britten et al 1992).        */
/*  speed        deg/s                                                       */
/*  theta        Angle in degrees.                                           */
/*  sizex        Size of rectangular frame for making/wrapping dots (pix)    */
/*  sizey        Size of rectangular frame for making/wrapping dots (pix)    */
/*  nframes      Time duration of stimulus (video frames).                   */
/*  dpf          Dots per frame.                                             */
/*  dt           Frames between refresh of signal dots (AKA ref_frames).     */
/*  seed         Address of randomization seed.                              */
/*                                                                           */
/*  Dot coordinates are returned in pixels (float) within a rectangle        */
/*  of 'sizex' by 'sizey' pixels.                                            */
/*                                                                           */
/*****************************************************************************/
void get_rds_disp_mod(sscale,framerate,disp0,disp_range,disp_sign,theta,tf,
		      sizex,sizey,dx,dy,dx0,dy0,nframes,dpf,seed,corr,
		      rdxl,rdyl,rdnl,rdxr,rdyr,rdnr,rdampl,rdampr)
     float sscale,framerate;
     float disp0;               // Disparity at center of range
     float disp_range;          // E.g., 1.0 ==> range from -1.0 to 1.0
     int disp_sign;             // 1 - approach, -1 - recede
     float theta,tf,sizex,sizey;
     float dx,dy;   // amount to shift rectangle
     float dx0,dy0;  // Start of range for shifting, if 'tf' not 0
     int nframes,dpf,seed;
     int corr;                   // 1-correlation, 0-uncorrelated across eyes
     float ***rdxl,***rdyl;      // *[nframes][dpf]
     int **rdnl;                 // *[nframes]
     float ***rdxr,***rdyr;      // *[nframes][dpf]
     int **rdnr;                 // *[nframes]
     float ***rdampl,***rdampr;  // *[nframes][dpf]
{
  int i,k;
  int fi,n,tn,seed_f,*dotn,nfpc,*dn;
  float **dotx1,**doty1,**dotx2,**doty2,**dota,speed;
  float **dlx,**dly,**drx,**dry,**da;

  //
  //  Given disp0, disp1, tf, compute the speed
  //
  speed = (disp_range * 2.0) * tf;   //  Units:  (deg/cyc)*(cyc/sec) = deg/sec

  if (disp_sign != 1)
    speed *= -1.0;  // Change the direction of motion

  nfpc = my_rint(framerate / tf);  // Units:  (frames/s) * (s/cyc) = frame/cyc


  //
  //  Get one cycle of the movie
  //
  get_dot_opp_mo(sscale,framerate,speed,disp0,theta,sizex,sizey,
		 nfpc,dpf,seed,corr,&dotx1,&doty1,&dotx2,&doty2,&dotn,&dota);

  //
  //  The final stimulus may contain multiple repeats of the movie
  //
  dlx = (float **)myalloc(nframes*sizeof(float *));
  dly = (float **)myalloc(nframes*sizeof(float *));
  drx = (float **)myalloc(nframes*sizeof(float *));
  dry = (float **)myalloc(nframes*sizeof(float *));
  da  = (float **)myalloc(nframes*sizeof(float *));
  dn  = (int    *)myalloc(nframes*sizeof(int));

  k = 0;
  for(i=0;i<nframes;i++){
    dlx[i] = copy_farray(dotx1[k],dotn[k]);
    dly[i] = copy_farray(doty1[k],dotn[k]);
    drx[i] = copy_farray(dotx2[k],dotn[k]);
    dry[i] = copy_farray(doty2[k],dotn[k]);
    da[i]  = copy_farray(dota[k], dotn[k]);
    dn[i]  = dotn[k];
    k += 1;
    if (k == nfpc){
      k = 0; // Start again at beginning of movie
    }
  }

  free_2d_farray(dotx1,nfpc);
  free_2d_farray(doty1,nfpc);
  free_2d_farray(dotx2,nfpc);
  free_2d_farray(doty2,nfpc);
  free_2d_farray(dota,nfpc);
  myfree(dotn);

  *rdxl = dlx;
  *rdyl = dly;
  *rdnl = dn;

  *rdxr = drx;
  *rdyr = dry;
  *rdnr = copy_iarray(dn,nframes);

  *rdampl = da;
  *rdampr = copy_2d_farray(da,nframes,dpf);
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_VSK_SUM_OF_SINUSOIDS                         */
/*                                                                           */
/*  Return the value of the Victor, Shapley, Knight (1977) sum of sinusoid   */
/*  signal at time 't' in seconds.                                           */
/*                                                                           */
/*  NOTES                                                                    */
/*  - All components have the sample amplitude.                              */
/*  - The returned signal varies from -1 to 1.                               */
/*  - 'p' specifies the pertubation component from 1..n, -1=no perturbing    */
/*                                                                           */
/*****************************************************************************/
float get_vsk_sum_of_sinusoids(snum,p,pratio,t)
     int snum,p;
     float pratio,t;
{
  int i;
  int n;
  float f,y,a,amp,pamp,npamp;

  if ((snum <0)||(snum > 2))
    exit_error("GET_VSK_SUM_OF_SINUSOIDS","Unknown sum index");

  n = VSK_sos_n[snum];

  amp = 1.0/(float)n; // Regular weight for each component
  npamp = 1.0/(pratio + (float)(n-1));
  pamp = pratio*npamp;

  y = 0.0;
  for(i=0;i<n;i++){
    if (p < 1)
      a = amp; // No perturbation
    else if ((i+1) == p)
      a = pamp; // Perturbation amplitude
    else
      a = npamp; // Non-perturbation amplitude

    //printf("a=%.4f amp= %.4f  npamp= %.4f pamp= %.4f\n",a,amp,npamp,pamp);

    if (snum == 0)      f = (float)VSK_sos_0[i];
    else if (snum == 1) f = (float)VSK_sos_1[i];
    else if (snum == 2) f = (float)VSK_sos_2[i];
    else                f = 0.0; // Won't happen
    f *= VSK_sos_unit;
    
    y += a * sin(f*t*M_PI*2.0);
  }
  return y;
}
/**************************************-**************************************/
/*                                                                           */
/*                               STIMU_RANGRID                               */
/*                                                                           */
/*  Return x,y coordinates for 'n' points chosen from a grid that is         */
/*  'gw' by 'gw' points, without replacement.  Add uniform noise of radius   */
/*  'jitr' to the coordinates.                                               */
/*                                                                           */
/*****************************************************************************/
void stimu_rangrid(sscale,x0,y0,gw,gn,n,jitr,seed,rx,ry)
     float sscale;  // (deg/pix)
     float x0,y0;   // Lower left corner of grid (deg)
     float gw;      // grid width (deg)
     int gn;        // Number of grid boxes
     int n;         // Number of coords to pick
     float jitr;    // Amount of uniform noise to add to coordinates (radius)
     int seed;      // Randomization seed
     float **rx;    // [n] Return x-coords (pix)
     float **ry;    // [n] Return y-coords (pix)
{
  int i;
  int *shuff,xi,yi,jseed;
  float *x,*y,jr2,dg,jx,jy,xdeg,ydeg;

  /*
    printf("ssacle = %f\n",sscale);
    printf("gw = %f\n",gw);
    printf("gn = %d\n",gn);
    printf("n = %d\n",n);
    printf("seed = %d\n",seed);
    printf("jitr = %f\n",jitr);
  */

  if (n > gn*gn)
    exit_error("STIMU_RANGRID","Too many grid points requested");

  x = (float *)myalloc(n*sizeof(float));
  y = (float *)myalloc(n*sizeof(float));

  //
  //  Pick a seed for adding jitter noise
  //
  if (seed > 0)
    jseed = -(seed + 12307);
  else
    jseed = seed - 12307;

  shuff = get_shuffle_index(gn*gn,7,seed);

  jr2 = jitr*2.0;          // Jitter diameter (deg)
  dg = gw/(float)(gn-1);   // Grid box size (deg)

  for(i=0;i<n;i++){
    xi = shuff[i] / gn;
    yi = shuff[i] - (xi * gn);

    jx = -jitr + jr2*myrand_util_ran2(&jseed);  // (deg)
    jy = -jitr + jr2*myrand_util_ran2(&jseed);  // (deg)

    xdeg = (float)xi*dg + jx;  // Coord rel. to lower left grid box = 0 (deg)
    ydeg = (float)yi*dg + jy;  // Coord rel. to lower left grid box = 0 (deg)

    //printf("%d (x,y) = %f %f\n",i,xdeg,ydeg);

    x[i] = (x0 + xdeg) / sscale;  // final coord (pix)
    y[i] = (y0 + ydeg) / sscale;  // final coord (pix)
  }

  myfree(shuff);
  *rx = x;
  *ry = y;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_STIM_1D_PPL_PULSE                          */
/*                                                                           */
/*  Return a 1D pulse.                                                       */
/*                                                                           */
/*****************************************************************************/
float *get_stim_1d_ppl_pulse(ppl,n,tscale,tsamp)
     struct param_pair_list *ppl;
     int n;
     float tscale;
     int tsamp;
{
  int i;
  int i0;
  float *data;
  float t0,tn,amp0,amp1,offset0,t;

  t0       = paramfile_get_float_param_or_exit(ppl,"st0");
  tn       = paramfile_get_float_param_or_exit(ppl,"stn");
  amp0     = paramfile_get_float_param_or_exit(ppl,"amp0");

  if (paramfile_test_param(ppl,"amp1"))
    amp1     = paramfile_get_float_param_or_exit(ppl,"amp1");
  else{
    offset0  = paramfile_get_float_param_or_exit(ppl,"offset0");
    amp1 = amp0 + offset0;
  }

  data = (float *)myalloc(n*sizeof(float));

  i0 = my_rint(t0 / tscale);
  for(i=0;i<n;i++){
    t = (float)(i-i0) * tscale; // Time in sec
    if (((i%tsamp) == 0) && (t >= 0.0) && (t < tn)){
      data[i] = amp1;
    }else
      data[i] = amp0;
  }
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_STIM_1D_PPL_SINE                          */
/*                                                                           */
/*  Return a 1D sinusoidal grating.                                          */
/*                                                                           */
/*    tscale - temporal scale (sec/frame)                                    */
/*    tf - temporal frequency (cyc/sec)                                      */
/*    phase - phase offset (deg)                                             */
/*    theta - orientation of grating (deg)                                   */
/*                                                                           */
/*  Returned array dimensions start at 0.                                    */
/*                                                                           */
/*****************************************************************************/
float *get_stim_1d_ppl_sine(ppl,n,tscale,tsamp)
     struct param_pair_list *ppl;
     int n;
     float tscale;
     int tsamp;
{
  int i;
  int i0;
  float *data,ph,f;
  float t0,tn,tf,phase,ampl,offset,bgval,t;

  t0       = paramfile_get_float_param_or_exit(ppl,"st0");
  tn       = paramfile_get_float_param_or_exit(ppl,"stn");
  tf       = paramfile_get_float_param_or_exit(ppl,"tf");
  phase    = paramfile_get_float_param_or_exit(ppl,"phase");
  ampl     = paramfile_get_float_param_or_exit(ppl,"ampl");
  offset   = paramfile_get_float_param_or_exit(ppl,"offset");
  bgval    = paramfile_get_float_param_or_exit(ppl,"bgval");

  data = (float *)myalloc(n*sizeof(float));

  f = tf * 2.0*M_PI;      // Freq in radians/sec
  ph = phase*M_PI/180.0;  // Phase in radians

  i0 = my_rint(t0 / tscale);
  for(i=0;i<n;i++){
    t = (float)(i-i0) * tscale; // Time in sec
    if (((i%tsamp) == 0) && (t >= 0.0) && (t < tn)){
      data[i] = offset + ampl * sin(f*t-ph);
    }else
      data[i] = bgval; // Wyeth - this confounds background and 'no stim'
  }
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_STIM_1D_PPL_VSK_SOS                         */
/*                                                                           */
/*  Sum of sinusoids.                                                        */
/*                                                                           */
/*    tscale - temporal scale (sec/frame)                                    */
/*    tf - temporal frequency (cyc/sec)                                      */
/*    phase - phase offset (deg)                                             */
/*    theta - orientation of grating (deg)                                   */
/*                                                                           */
/*  Returned array dimensions start at 0.                                    */
/*                                                                           */
/*****************************************************************************/
float *get_stim_1d_ppl_vsksos(ppl,n,tscale,tsamp)
     struct param_pair_list *ppl;
     int n;
     float tscale;
     int tsamp;
{
  int i;
  int i0,vskn,pcomp;
  float *data,pratio;
  float t0,tn,ampl,offset,bgval,t;

  t0     = paramfile_get_float_param_or_exit(ppl,"st0");
  tn     = paramfile_get_float_param_or_exit(ppl,"stn");
  ampl   = paramfile_get_float_param_or_exit(ppl,"ampl");
  offset = paramfile_get_float_param_or_exit(ppl,"offset");
  bgval = paramfile_get_float_param_or_exit(ppl,"bgval");
  vskn   =   paramfile_get_int_param_or_exit(ppl,"vskn");
  pcomp  =   paramfile_get_int_param_or_exit(ppl,"pcomp",-1);
  pratio = paramfile_get_float_param_or_exit(ppl,"pratio",1.0);

  data = (float *)myalloc(n*sizeof(float));

  i0 = my_rint(t0 / tscale);
  for(i=0;i<n;i++){
    t = (float)(i-i0) * tscale; /* Time in sec */
    if (((i%tsamp) == 0) && (t >= 0) && (t < tn)){
      data[i] = offset + ampl * get_vsk_sum_of_sinusoids(vskn,pcomp,pratio,t);
    }else
      data[i] = bgval; /* Wyeth - this confounds background and 'no stim' */
  }
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_STIM_1D_PPL_NOISE                         */
/*                                                                           */
/*  Return a 1D sinusoidal grating.                                          */
/*                                                                           */
/*    tscale - temporal scale (sec/frame)                                    */
/*    tf - temporal frequency (cyc/sec)                                      */
/*    phase - phase offset (deg)                                             */
/*    theta - orientation of grating (deg)                                   */
/*                                                                           */
/*  Returned array dimensions start at 0.                                    */
/*                                                                           */
/*****************************************************************************/
float *get_stim_1d_ppl_noise(ppl,n,tscale,tsamp)
     struct param_pair_list *ppl;
     int n;
     float tscale;
     int tsamp;
{
  int i;
  int i0,dt,seed,distrib,*rseq,nrseq,cnt,ri;
  float *data,rv,*gseq;
  float t0,tn,mu,sig,bgval,t;

  t0       = paramfile_get_float_param_or_exit(ppl,"st0");
  tn       = paramfile_get_float_param_or_exit(ppl,"stn");
  mu       = paramfile_get_float_param_or_exit(ppl,"mu");
  sig      = paramfile_get_float_param_or_exit(ppl,"sig");
  bgval   = paramfile_get_float_param_or_exit(ppl,"bgval");
  dt       =   paramfile_get_int_param_or_exit(ppl,"dt");
  seed     =   paramfile_get_int_param_or_exit(ppl,"seed");
  distrib  =   paramfile_get_int_param_or_exit(ppl,"distrib");

  if ((distrib < 1) || (distrib > 2))
    exit_error("GET_STIM_1D_PPL_NOISE","Error - must use distrib 1 or 2");

  data = (float *)myalloc(n*sizeof(float));

  //  Test run to see how many RV's are needed
  cnt = dt;
  ri = 0;
  i0 = my_rint(t0 / tscale);
  for(i=0;i<n;i++){
    t = (float)(i-i0) * tscale; // Time in sec
    if (((i%tsamp) == 0) && (t >= 0) && (t < tn)){
      if (cnt==dt){
	ri += 1;
	cnt = 0;
      }
      cnt += 1;
    }
  }
  nrseq = ri;
  //printf("  nrseq = %d\n",nrseq);

  rseq = NULL;
  gseq = NULL;
  if (distrib == 1)
    gseq = myrand_get_std_gauss_seq(nrseq,seed);  // Guassian - floats
  else if (distrib == 2)
    rseq = myrand_get_std_bin_seq(nrseq,seed);    // Binary - ints

  rv = -1.0; // Avoid -Wall warning

  cnt = dt;
  ri = 0;
  i0 = my_rint(t0 / tscale);
  for(i=0;i<n;i++){
    t = (float)(i-i0) * tscale; // Time in sec
    if (((i%tsamp) == 0) && (t >= 0) && (t < tn)){
      if (cnt==dt){
	if (ri >= nrseq)
	  exit_error("GET_STIM_1D_PPL_NOISE","Exceeded random sequence");

	if (distrib == 1){ // Gaussian
	  rv = gseq[ri] * sig;
	}else if (distrib == 2){
	  if (rseq[ri])
	    rv = sig;
	  else
	    rv = -sig;
	}

	ri += 1;
	cnt = 0;
      }
      data[i] = mu + rv;

      cnt += 1;
    }else{
      data[i] = bgval;
    }
  }

  if (ri != nrseq){
    printf(" ri = %d   nrseq = %d\n",ri,nrseq);
    exit_error("GET_STIM_1D_PPL_NOISE","Excess random sequence");
  }

  if (distrib == 1)
    myfree(gseq);
  else if (distrib == 2)
    myfree(rseq);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_STIM_1D_PPL_STAT1                           */
/*                                                                           */
/*  Emulate the 'stat1_anti' experiments.                                    */
/*                                                                           */
/*****************************************************************************/
float *get_stim_1d_ppl_stat1(ppl,n,tscale,tsamp)
     struct param_pair_list *ppl;
     int n;
     float tscale;
     int tsamp;
{
  int i,j,k;
  int seed,nblk,ndur,bi,*ton,*toff,ts,t0,t1,done,**ndx,seqflag;
  int nmask,cn;
  float *data,tau,tsep,ampl,*emask,*ad;

  tau   = paramfile_get_float_param_or_exit(ppl,"tau");  // decay time (s)
  seed  = paramfile_get_int_param_or_exit(ppl,"seed");
  nblk  = paramfile_get_int_param_or_exit(ppl,"nblock"); // Num repeats
  ampl  = paramfile_get_float_param_or_exit(ppl,"ampl"); // -1..1
  tsep  = paramfile_get_float_param_or_exit(ppl,"tsep"); // t between stim (s)
  seqflag  = paramfile_get_int_param_or_exit(ppl,"seq_type");

  // Get ON/OFF durations for each stimulus, and blockwise random indices
  stm_stat1d_get_seq(NULL,nblk,seqflag,seed,tscale,&ndur,&ton,&toff,&ndx);

  ts = my_rint(tsep/tscale);

  data = (float *)myalloc(n*sizeof(float));

  bi = 0; // Block index
  i = 0;  // Time index
  done = 0;
  while(!done){
    j = 0;  // Stimulus index, w/i block
    while(done == 0){
      t1 = ton[ndx[bi][j]];
      t0 = toff[ndx[bi][j]];

      if ((i+t1+t0+ts) >= n)  // Not enough time for this
	done = 1;
      else{
	for(k=0;k<ts;k++){   // Separation
	  data[i] = 0.0;
	  i += 1;
	}
	for(k=0;k<t1;k++){   // Stimulus
	  data[i] = ampl;
	  i += 1;
	}
	for(k=0;k<t0;k++){   // After-period
	  data[i] = 0.0;
	  i += 1;
	}
      }
      
      j += 1;
      if (j >= ndur)
	done = 2;
    }

    if (done == 2)  // Ended after using all stimuli in this block
      done = 0;
    else{
      printf("*** GET_STIM_1D_PPL_STAT1  Warning, Cannot fit\n");
      printf("*** %d blocks in duration %f\n",nblk,(float)n*tscale);
    }

    bi += 1;
    if (bi >= nblk)
      done = 1;
  }
  free_2d_iarray(ndx,nblk);  myfree(ton);  myfree(toff);
  //append_farray_plot("zzz.stat1.pl","data_raw",data,n,1);


  //
  // Adaptation
  //

  // Make decaying exponential mask
  nmask = (int)(tau/tscale) * 10 + 1; /* Must be odd */
  cn = (nmask-1)/2;
  emask = get_zero_farray(nmask);
  for(i=cn;i<nmask;i++)
    emask[i] = func_one_sided_exp((double)(i-cn),0.0,(double)(tau/tscale));
  norm_area_farray(emask,nmask,1.0);
  //append_farray_plot("zzz.emask","emask",emask,nmask,1);


  // Apply mask
  ad = fft_convolve(data,n,emask,nmask,1);
  //append_farray_plot("zzz.stat1.pl","convolved",ad,n,1);
  for(i=0;i<n;i++)
    data[i] +=  -ad[i];
  //append_farray_plot("zzz.stat1.pl","data_ad",data,n,1);
  myfree(ad); myfree(emask);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MAKE_FRAME_FLAT                             */
/*                                                                           */
/*  Fill frame "frame_num" with "value".                                     */
/*                                                                           */
/*****************************************************************************/
void make_frame_flat(data,x0,xn,y0,yn,z0,zn,frame_num,value)
     float ***data;
     int x0,xn,y0,yn,z0,zn,frame_num;
     float value;
{
  int i,j;

  for(i=x0;i<(x0+xn);i++)
    for(j=y0;j<(y0+yn);j++)
      data[i][j][frame_num] = value;
}
/**************************************-**************************************/
/*                                                                           */
/*                        MAKE_FRAME_SCALE_FRAME_ZMID                        */
/*                                                                           */
/*  Fill frame "frame_num" with a scaled version of frame 'k', but before    */
/*  scaling, add 1, then subtract 1 after.                                   */
/*                                                                           */
/*  ****** THIS IS FOR DECAY TO BLACK of zero-mid frames.                    */
/*                                                                           */
/*****************************************************************************/
void make_frame_scale_frame_zmid(data,x0,xn,y0,yn,z0,zn,frame_num,k,scale,
				 blankval)
     float ***data;
     int x0,xn,y0,yn,z0,zn,frame_num,k;
     float scale;
     float blankval;  // typically -1.0, for black, could be 0 for gray
{
  int i,j;
  float yoff;

  // blankval -1  ==>  1.0
  // blankval  0  ==>  0.0

  yoff = -blankval;

  for(i=x0;i<(x0+xn);i++)
    for(j=y0;j<(y0+yn);j++)
      data[i][j][frame_num] = scale * (data[i][j][k] + yoff) - yoff;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MAKE_FRAME_FROM_3D                            */
/*                                                                           */
/*  Fill 'data' frame 'fi' with values from frame 'k' in 'd'.                */
/*                                                                           */
/*****************************************************************************/
void make_frame_from_3d(data,x0,xn,y0,yn,z0,zn,fi,d,dxn,dyn,dtn,k,winx0,winy0,
			winxn,winyn,destx0,desty0,bgval,maxlum)
     float ***data;                 /* Destination of stimulus data */
     int x0,xn,y0,yn,z0,zn,fi;      /* */
     float ***d;                    /* Source of data */
     int dxn,dyn,dtn,k;             /* */
     int winx0,winy0,winxn,winyn;   /* Window to take from 'd' */ 
     int destx0,desty0;             /* Where to put the window in 'data' */
     float bgval,maxlum;           /* Value for filling; maximum value */
{     
  int i,j;
  int winx1,winy1;

  /*
    printf("dat====>  %d %d %d %d %d %d %d\n",x0,xn,y0,yn,z0,zn,fi);
    printf("d  ====>  %d %d %d %d\n",dxn,dyn,dtn,k);
    printf("win====>  %d %d %d %d\n",winx0,winy0,winxn,winyn);
    printf("dst====>  %d %d\n",destx0,desty0);*/

  winx1 = winx0 + winxn;
  winy1 = winy0 + winyn;

  if ((winx0 < 0) || (winy0 < 0) || (winx1 > dxn) || (winy1 > dyn))
    exit_error("MAKE_FRAME_FROM_3D","Source window coordinate error");

  if ((destx0 < x0) || (desty0 < y0) || ((destx0+winxn) > (x0+xn)) ||
      ((desty0+winyn) > (y0+yn)))
    exit_error("MAKE_FRAME_FROM_3D","Destination coordinate error");

  // If filling-in around the window is needed
  if ((destx0 > x0) || (desty0 > y0) || (winxn < xn)  || (winyn < yn))
    make_frame_flat(data,x0,xn,y0,yn,z0,zn,fi,bgval);

  for(i=winx0;i<winx1;i++)
    for(j=winy0;j<winy1;j++)
      data[destx0+i][desty0+j][fi] = d[i][j][k];
}
/**************************************-**************************************/
/*                                                                           */
/*                           STIMU_APPLY_MASK1_VAL                           */
/*                                                                           */
/*  For a single-frame mask.                                                 */
/*                                                                           */
/*  Mix the image data 'd' with the mask 'm' using the value 'v' as such:    */
/*                                                                           */
/*    final_value(i,j,k) = m(i,j) * d(i,j,k)  +  [1 - m(i,j)] v              */
/*                                                                           */
/*****************************************************************************/
void stimu_apply_mask1_val(d,m,xn,yn,tn,v)
     float ***d;   // [1..xn][1..yn][1..tn] image data to be masked
     float **m;    // [0..xn-1][0..yn-1] mask to apply, values should be [0..1]
     int xn,yn,tn; // Stimulus data size
     float v;      // Masking value
{
  int i,j,k;
  float mt,dnew;

  // **** WYETH - NOTE  TENSOR RANGE   [1...xn]  for STIMULUS but Not MASK
  // **** WYETH - NOTE  TENSOR RANGE   [1...yn]  for STIMULUS but Not MASK
  // **** WYETH - NOTE  TENSOR RANGE   [1...tn]  for STIMULUS but Not MASK

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      //mt = m[1+i][1+j];
      mt = m[i][j];
      if (mt < 1.0){  // Otherwise, don't change the image
	for(k=0;k<tn;k++){
	  if (mt <= 0.0)
	    dnew = v;
	  else
	    dnew = mt * d[1+i][1+j][1+k] + (1.0 - mt) * v;
	  d[1+i][1+j][1+k] = dnew;
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            STIMU_APPLY_MASK_VAL                           */
/*                                                                           */
/*  For a multi-frame mask.                                                  */
/*                                                                           */
/*  Mix the image data 'd' with the mask 'm' using the value 'v' as such:    */
/*                                                                           */
/*    final_value(i,j,k) = m(i,j,k) * d(i,j,k)  +  [1 - m(i,j,k)] v          */
/*                                                                           */
/*****************************************************************************/
void stimu_apply_mask_val(d,m,xn,yn,tn,v)
     float ***d;   // [xn][yn][tn] image data to be masked
     float ***m;   // [xn][yn][tn] mask ot apply, values should be [0..1]  
     int xn,yn,tn; // Image and mask size
     float v;      // Masking value
{
  int i,j,k;
  float mt,dnew;

  // ******* WYETH - NOTE  TENSOR RANGE   [1...xn]
  // ******* WYETH - NOTE  TENSOR RANGE   [1...yn]  for BOTH stim and mask
  // ******* WYETH - NOTE  TENSOR RANGE   [1...tn]

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<tn;k++){
	//printf("i,j,k %d %d %d\n",i,j,k);
	mt = m[1+i][1+j][1+k];
	//printf("mt = %f\n",mt);
	if (mt < 1.0){  // Otherwise, don't change the image
	  if (mt <= 0.0)
	    dnew = v;
	  else
	    dnew = mt * d[1+i][1+j][1+k] + (1.0 - mt) * v;
	  d[1+i][1+j][1+k] = dnew;
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               STIMU_FST_MASK                              */
/*                                                                           */
/*  Mix the image data 'd' with the mask 'm' using the value 'v' as such:    */
/*                                                                           */
/*    final_value(i,j) = m(i,j) * d(i,j)  +  [1 - m(i,j)] v                  */
/*                                                                           */
/*****************************************************************************/
void stimu_fst_mask(sdata,xn,yn,tn,mask,xm,ym,tm,val)
     float ***sdata;   // [xn][yn][tn] stim data to be masked.  INDEXING 1...
     int xn,yn,tn;     // Stimulus size
     float ***mask;    // [xm][ym][tm] mask to apply.  INDEXING 1...
     int xm,ym,tm;     // Mask size
     float val;        // Masking value
{
  int i,j;
  float **m2d;

  //printf("  STIMU_FST_MASK\n");
  //printf("    mask:  %d %d %d (xn yn tn)\n",xm,ym,tm);
  //printf("   sdata:  %d %d %d (xn yn tn)\n",xn,yn,tn);
  //printf("   making values:  %f\n",val);

  if ((xn != xm) || (yn != ym))
    exit_error("STIMU_FST_MASK","Stimulus and mask differ in 'xn' or 'yn'");

  if (tm == 1){
    //
    //  The mask is only one frame long
    //
    m2d = get_2d_from_3d_farray(mask,1,xn,1,yn,1);  // Get the 2D mask

    stimu_apply_mask1_val(sdata,m2d,xn,yn,tn,val);  // Apply 2D mask

    free_2d_farray(m2d,xn);  // Free the 2D mask
  }else{
    if (tn != tm)
      exit_error("STIMU_FST_MASK",
		 "Multi-frame mask does not match stimulus 'tn'");

    stimu_apply_mask_val(sdata,mask,xn,yn,tn,val);  // Apply 3D mask
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            FRAME_APPLY_APERTURE                           */
/*                                                                           */
/*  Draw an aperture.                                                        */
/*                                                                           */
/*  aptype                                                                   */
/*   -1 - hole (draw gray disk)                                              */
/*    0 - no aperture                                                        */
/*    1 - circular (w is diameter)                                           */
/*    2 - rectangular                                                        */
/*    3 - tiltrect                                                           */
/*  xc,yc - center of aperture                                               */
/*  w,h   - width and height of aperture                                     */
/*  r     - angle of aperture                                                */
/*                                                                           */
/*****************************************************************************/
void frame_apply_aperture(data,x0,xn,y0,yn,z0,zn,frame_num,aptype,xc,yc,w,h,r,
			  value)
     float ***data;
     int x0,xn,y0,yn,z0,zn,frame_num,aptype;
     float xc,yc;   // Center of aperture (pix) wrt 0..xn, 0..yn
     float w,h;     // width and height of aperture (pix)
     float r;       // angle of aperture (deg)
     float value;   // value to assign outside of aperture
{
  int i,j;
  int flag;
  float xf,yf,r2,t2,x1,x2,y1,y2;
  float a,b,c,d,thetar,xo,yo;

  if (aptype == 0)
    return;
  else if (aptype == 3){
    printf("   *** Aptype 3 is Gaussian, implemented elsewhere.\n");
    return;
  }

  r2 = w*w/4.0;
  x1 = -w/2.0;  x2 = w/2.0;
  y1 = -h/2.0;  y2 = h/2.0;

  if (aptype == 4){
    thetar = -r * M_PI/180.0;
    a =  cos(thetar); b = sin(thetar); // rotation matrix
    c = -sin(thetar); d = cos(thetar);
  }

  for(i=0;i<xn;i++){
    xf = (float)i-xc;
    for(j=0;j<yn;j++){
      yf = (float)j-yc;  // WYETH BUG FIXED
      //yf = (float)j-xc;
      
      flag = 0;
      if (aptype == 1){
	t2 = xf*xf + yf*yf;
	if (t2 > r2)
	  flag = 1;
      }else if (aptype == -1){
	t2 = xf*xf + yf*yf;
	if (t2 <= r2)
	  flag = 1;
      }else if (aptype == 2){
	if ((xf < x1)||(xf >= x2)||(yf < y1)||(yf >= y2))
	  flag = 1;
      }else if (aptype == -2){
	if ((xf >= x1)&&(xf < x2)&&(yf >= y1)&&(yf < y2))
	  flag = 1;
      }else if (aptype == 3){
	;  // Aptype 3 is Gaussian, implemented elsewhere
      }else if (aptype == 4){

	xo = a*xf + c*yf;  // rotate back
	yo = b*xf + d*yf;

	flag = 1;  // Assume background
	if ((xo > -w/2.0) && (xo < w/2.0) &&
	    (yo > -h/2.0) && (yo < h/2.0))
	  flag = 0;  // foreground

      }else
	exit_error("FRAME_APPLY_APERTURE","Unknown aptype");
      
      if (flag)
	data[1+i][1+j][frame_num] = value;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            APPLY_APERTURE_SCALE                           */
/*                                                                           */
/*  Draw an aperture.                                                        */
/*                                                                           */
/*  aptype                                                                   */
/*   -1 - hole (draw gray disk)                                              */
/*    0 - no aperture                                                        */
/*    1 - circular (w is diameter)                                           */
/*    2 - rectangular                                                        */
/*    3 - tiltrect                                                           */
/*  xc,yc - center of aperture                                               */
/*  w,h   - width and height of aperture                                     */
/*  r     - angle of aperture                                                */
/*                                                                           */
/*****************************************************************************/
void apply_aperture_scale(data,xn,yn,aptype,xc,yc,w,h,r,value,sval)
     float **data;
     int xn,yn;
     int aptype;
     float xc,yc;   // Center of aperture (pix) wrt 0..xn, 0..yn
     float w,h;     // width and height of aperture (pix)
     float r;       // angle of aperture (deg)
     float value;   // value to assign outside of aperture
     float sval;    // scale result by this number
{
  int i,j;
  int flag;
  float xf,yf,r2,t2,x1,x2,y1,y2;
  float a,b,c,d,thetar,xo,yo;

  /*
  printf("  xn,yn = %d %d\n",xn,yn);
  printf("  aptype = %d\n",aptype);
  printf("  xc,yc = %f %f\n",xc,yc);
  printf("  w,h = %f %f\n",w,h);
  printf("  r = %f\n",r);
  printf("  value = %f\n",value);
  printf("  sval = %f\n",sval);
  */

  if (aptype == 0){
    if (sval != 1.0){
      for(i=0;i<xn;i++){
	for(j=0;j<yn;j++){
	  data[i][j] *= sval;    // Scale image everywhere
	}
      }
    }
    return;
  }else if (aptype == 3){
    printf("   *** Aptype 3 is Gaussian, implemented elsewhere.\n");
    return;
  }

  r2 = w*w/4.0;
  x1 = -w/2.0;  x2 = w/2.0;
  y1 = -h/2.0;  y2 = h/2.0;

  if (aptype == 4){
    thetar = -r * M_PI/180.0;
    a =  cos(thetar); b = sin(thetar); // rotation matrix
    c = -sin(thetar); d = cos(thetar);
  }

  for(i=0;i<xn;i++){
    xf = (float)i-xc;
    for(j=0;j<yn;j++){
      yf = (float)j-yc;
      //yf = (float)j-xc;
      
      flag = 0;
      if (aptype == 1){
	t2 = xf*xf + yf*yf;
	if (t2 > r2)
	  flag = 1;
      }else if (aptype == -1){
	t2 = xf*xf + yf*yf;
	if (t2 <= r2)
	  flag = 1;
      }else if (aptype == 2){
	if ((xf < x1)||(xf >= x2)||(yf < y1)||(yf >= y2))
	  flag = 1;
      }else if (aptype == -2){
	if ((xf >= x1)&&(xf < x2)&&(yf >= y1)&&(yf < y2))
	  flag = 1;
      }else if (aptype == 3){
	; // Aptype 3 is Gaussian, implemented elsewhere
      }else if (aptype == 4){

	xo = a*xf + c*yf;  // rotate back
	yo = b*xf + d*yf;

	flag = 1;  // Assume background
	if ((xo > -w/2.0) && (xo < w/2.0) &&
	    (yo > -h/2.0) && (yo < h/2.0))
	  flag = 0;  // foreground

      }else
	exit_error("APPLY_APERTURE","Unknown aptype");
      
      if (flag)
	//data[i][j] = sval * value; // Apply aperture WYETH REMOVED 2014 Nov
	data[i][j] = value;   // Apply value outside aperture
      else
	data[i][j] *= sval;   // Scale image inside aperture

    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               MULTIPLY_FRAME                              */
/*                                                                           */
/*****************************************************************************/
void multiply_frame(data,x0,xn,y0,yn,z0,zn,frame_num,value)
     float ***data;
     int x0,xn,y0,yn,z0,zn,frame_num;
     float value;
{
  int i,j;

  for(i=x0;i<(x0+xn);i++)
    for(j=y0;j<(y0+yn);j++)
      data[i][j][frame_num] *= value;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MULTIPLY_FRAME_GAUSSIAN                        */
/*                                                                           */
/*  Multiply frame "frame_num" with a Gaussian.                              */
/*                                                                           */
/*    x0 - x center of Gaussian (degr)                                       */
/*    y0 - y center of Gaussian (degr)                                       */
/*    sigma - SD of Gaussian (degr)                                          */
/*    sscale - degrees (of vis. angle) / pixel                               */
/*                                                                           */
/*****************************************************************************/
void multiply_frame_gaussian(data,i1,i2,i3,n1,n2,n3,frame_num,xc,yc,sigma,
			     sscale)
     float ***data;
     int i1,i2,i3,n1,n2,n3,frame_num;
     float xc,yc;                     // center (pix)  [Mar 16, 2009]
     float sigma;                     // deg
     float sscale;                    // deg/pix
{
  int i,j;
  float x,y;

  for(i=0;i<n1;i++){
    x = ((float)i-xc)*sscale;
    for(j=0;j<n2;j++){
      y = ((float)j-yc)*sscale;
      data[i1+i][i2+j][frame_num] *= 
	(float)func_2d_gaussian(x,y,0.0,0.0,sigma,sigma,0);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          MULTIPLY_FRAME_GAUSSIAN_0                        */
/*                                                                           */
/*  Multiply frame "frame_num" with a Gaussian.                              */
/*                                                                           */
/*    x0 - x center of Gaussian (degr)                                       */
/*    y0 - y center of Gaussian (degr)                                       */
/*    sigma - SD of Gaussian (degr)                                          */
/*    sscale - degrees (of vis. angle) / pixel                               */
/*                                                                           */
/*****************************************************************************/
void multiply_frame_gaussian_0(data,n1,n2,xc,yc,sigma,sscale)
     float **data;
     int n1,n2;
     float xc,yc;                     // center (pix)  [Mar 16, 2009]
     float sigma;                     // deg
     float sscale;                    // deg/pix
{
  int i,j;
  float x,y;

  for(i=0;i<n1;i++){
    x = ((float)i-xc)*sscale;
    for(j=0;j<n2;j++){
      y = ((float)j-yc)*sscale;
      data[i][j] *= (float)func_2d_gaussian(x,y,0.0,0.0,sigma,sigma,0);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MULTIPLY_FRAME_SINUSOID                        */
/*                                                                           */
/*  Multiply frame "frame_num" with a sinusoid.                              */
/*                                                                           */
/*    freq - cycles/degree of vis. angle                                     */
/*    theta - degrees (in the circle)                                        */
/*    sscale - degrees (of vis. angle) / pixel                               */
/*    phase - in distance (degrees of *visual angle*)                        */
/*                                                                           */
/*****************************************************************************/
void multiply_frame_sinusoid(data,i1,i2,i3,n1,n2,n3,frame_num,freq,theta,
			     sscale,phase)
     float ***data;
     int i1,i2,i3,n1,n2,n3,frame_num;
     float freq,theta,sscale,phase;
{
  int i,j;
  int xc,yc;
  float x,y,nx,ny,twopif,ph;

  twopif = 2.0*M_PI*freq;
  ph = freq * phase * 2.0*M_PI;
  nx = cos(theta/180.0 * M_PI);
  ny = sin(theta/180.0 * M_PI);
  xc = (n1-1)/2;
  yc = (n2-1)/2;

  for(i=0;i<n1;i++){
    x = (float)(i-xc)*sscale;
    for(j=0;j<n2;j++){
      y = (float)(j-yc)*sscale;
      data[i1+i][i2+j][frame_num] *= cos(ph + twopif*(nx*x+ny*y));
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          DRAW_CIRC_SINE_IN_FRAME                          */
/*                                                                           */
/*  Draw a circular aperture sinewave in "frame_num".                        */
/*  Do not alter the background.                                             */
/*                                                                           */
/*    freq_c - cycles/degree of vis. angle                                   */
/*    theta_c - degrees (in the circle)                                      */
/*    sscale - degrees (of vis. angle) / pixel                               */
/*    phase_c - in degrees of a cycle [0,360)                                */
/*    con1 - contrast from -1 to 1.                                          */
/*    cx - center of patch along horizontal axis (degr)                      */
/*    cy - center of patch along vertical axis (degr)                        */
/*                                                                           */
/*  *** NOTE:  All values written into the frame are in [-1,1].              */
/*                                                                           */
/*****************************************************************************/
void draw_circ_sine_in_frame(data,x0,xn,y0,yn,z0,zn,frame_num,cx,cy,size,freq,
			     theta,phase,sscale,con)
     float ***data;
     int x0,xn,y0,yn,z0,zn,frame_num;
     float cx,cy,size,freq,theta,phase,sscale,con;
{
  int i,j;
  int xc,yc,i1,i2,j1,j2;
  float x,y,nx,ny,ph,twopif,r1s;

  r1s = size*size/4.0;

  twopif = 2.0*M_PI*freq;
  ph = phase/180.0 * M_PI;
  nx = cos(theta/180.0 * M_PI);
  ny = sin(theta/180.0 * M_PI);
  xc = (xn-1)/2 + (int)(cx/sscale);
  yc = (yn-1)/2 + (int)(cy/sscale);

  // *** WYETH - VERY INEFFICIENT - WE SHOULD USE THE PATCH DIAM TO
  // *** CONSTRAIN THE SIZE HERE, not just the full window size
  // *** SEE "   MAKE_FRAME_ADD_GPATCH_SINE   " below.
  // *** SEE "   MAKE_FRAME_ADD_GPATCH_SINE   " below.

  // Define the rectangular region which covers the circle and is within
  //   the frame boundary
  i1 = xc - xn/2 - 1; // Left x value
  if (i1 < 0) i1 = 0;
  i2 = xc + xn/2 + 1; // Right x value
  if (i2 >= xn) i2 = xn-1;

  j1 = yc - yn/2 - 1; // Min y value
  if (j1 < 0) j1 = 0;
  j2 = yc + yn/2 + 1; // Max y value
  if (j2 >= yn) j2 = yn-1;

  // Draw the sine wave in the circle
  for(i=i1;i<=i2;i++){
    x = (float)(i-xc)*sscale;
    for(j=j1;j<=j2;j++){
      y = (float)(j-yc)*sscale;
      if ((x*x + y*y) < r1s) // Inner patch
	data[x0+i][y0+j][frame_num] = con*cos(twopif*(nx*x+ny*y) - ph);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          MAKE_FRAME_CIRC_CS_SINE                          */
/*                                                                           */
/*  Fill frame "frame_num" with a circular patch of grating surrounded by a  */
/*  circular annulus.                                                        */
/*                                                                           */
/*    freq_c - cycles/degree of vis. angle                                   */
/*    freq_s - cycles/degree of vis. angle                                   */
/*    theta_c - degrees (in the circle)                                      */
/*    theta_s - degrees (in the circle)                                      */
/*    sscale - degrees (of vis. angle) / pixel                               */
/*    phase_c - in degrees of a cycle [0,360)                                */
/*    phase_s - in degrees of a cycle [0,360)                                */
/*    con1 - contrast from -1 to 1.                                          */
/*    con2 - contrast from -1 to 1.                                          */
/*                                                                           */
/*  *** NOTE:  All values written into the frame are in [-1,1].              */
/*                                                                           */
/*****************************************************************************/
void make_frame_circ_cs_sine(data,x0,xn,y0,yn,z0,zn,frame_num,size,idiam,odiam,
			     freq_c,freq_s,theta_c,theta_s,phase_c,phase_s,
			     sscale,con1,con2,bgampl,cx,cy,cs_flag)
     float ***data;
     int x0,xn,y0,yn,z0,zn,frame_num;
     float size,idiam,odiam,freq_c,freq_s,theta_c,theta_s,phase_c,phase_s;
     float sscale,con1,con2,bgampl;
     float cx,cy; // Center coords (deg)
     int cs_flag; // 0-neither, 1-ctr, 2-surr, 3-both
{
  int i,j;
  float xc,yc;
  float x,y,nx1,ny1,ph1,nx2,ny2,ph2,twopif1,twopif2,r1s,r2s,r3s;

  r1s = size*size/4.0;
  r2s = idiam*idiam/4.0;
  r3s = odiam*odiam/4.0;


  twopif1 = 2.0*M_PI*freq_c;
  twopif2 = 2.0*M_PI*freq_s;
  ph1 = phase_c/180.0 * M_PI;  // ph = freq * phase * 2.0*M_PI;
  nx1 = cos(theta_c/180.0 * M_PI);
  ny1 = sin(theta_c/180.0 * M_PI);
  ph2 = phase_s/180.0 * M_PI;  // ph = freq * phase * 2.0*M_PI;
  nx2 = cos(theta_s/180.0 * M_PI);
  ny2 = sin(theta_s/180.0 * M_PI);

  xc = (float)(xn-1)/2.0 *sscale + cx;  // Center (deg)
  yc = (float)(yn-1)/2.0 *sscale + cy;  // Center (deg)

  for(i=0;i<xn;i++){
    x = (float)i*sscale - xc;
    for(j=0;j<yn;j++){
      y = (float)j*sscale - yc;
      if ((x*x + y*y) < r1s){ // Inner patch
	if ((cs_flag == 1) || (cs_flag == 3))
	  data[x0+i][y0+j][frame_num] = con1*cos(twopif1*(nx1*x+ny1*y) - ph1);
	else
	  data[x0+i][y0+j][frame_num] = bgampl; // Background
      }else if ((x*x + y*y) < r2s) // between center and surround
	data[x0+i][y0+j][frame_num] = bgampl; // Background
      else if ((x*x + y*y) < r3s){ // surround
	if ((cs_flag == 2) || (cs_flag == 3))
	  data[x0+i][y0+j][frame_num] = con2*cos(twopif2*(nx2*x+ny2*y) - ph2);
	else
	  data[x0+i][y0+j][frame_num] = bgampl; // Background
      }else
	data[x0+i][y0+j][frame_num] = bgampl; // Background
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MAKE_FRAME_DISK_FLAT                           */
/*                                                                           */
/*  Fill frame "frame_num" with a circular disk surrounded by background     */
/*  color.                                                                   */
/*                                                                           */
/*    sscale - degrees (of vis. angle) / pixel                               */
/*                                                                           */
/*  *** NOTE:  All values written into the frame are in [-1,1].              */
/*                                                                           */
/*****************************************************************************/
void make_frame_disk_flat(data,x0,xn,y0,yn,z0,zn,frame_num,size,
			  sscale,ampl,bgval)
     float ***data;
     int x0,xn,y0,yn,z0,zn,frame_num;
     float size,sscale,ampl,bgval;
{
  int i,j;
  int xc,yc;
  float x,y,r1s;

  r1s = size*size/4.0;

  xc = (xn-1)/2;
  yc = (yn-1)/2;

  for(i=0;i<xn;i++){
    x = (float)(i-xc)*sscale;
    for(j=0;j<yn;j++){
      y = (float)(j-yc)*sscale;
      if ((x*x + y*y) < r1s) // Inner patch
	data[x0+i][y0+j][frame_num] = ampl;
      else 
	data[x0+i][y0+j][frame_num] = bgval; // Background
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          MAKE_FRAME_CIRC_CS_FLAT                          */
/*                                                                           */
/*  Fill frame "frame_num" with a circular disk surrounded by a circular     */
/*  annulus.                                                                 */
/*                                                                           */
/*    sscale - degrees (of vis. angle) / pixel                               */
/*                                                                           */
/*  *** NOTE:  All values written into the frame are in [-1,1].              */
/*                                                                           */
/*****************************************************************************/
void make_frame_circ_cs_flat(data,x0,xn,y0,yn,z0,zn,frame_num,size,idiam,odiam,
			     sscale,campl,sampl,bgval)
     float ***data;
     int x0,xn,y0,yn,z0,zn,frame_num;
     float size,idiam,odiam,sscale,campl,sampl,bgval;
{
  int i,j;
  int xc,yc;
  float x,y,r1s,r2s,r3s;

  r1s = size*size/4.0;
  r2s = idiam*idiam/4.0;
  r3s = odiam*odiam/4.0;

  xc = (xn-1)/2;
  yc = (yn-1)/2;

  for(i=0;i<xn;i++){
    x = (float)(i-xc)*sscale;
    for(j=0;j<yn;j++){
      y = (float)(j-yc)*sscale;
      if ((x*x + y*y) < r1s) // Inner patch
	data[x0+i][y0+j][frame_num] = campl;
      else if ((x*x + y*y) < r2s) // between center and surround
	data[x0+i][y0+j][frame_num] = bgval; // Background
      else if ((x*x + y*y) < r3s) // surround
	data[x0+i][y0+j][frame_num] = sampl;
      else 
	data[x0+i][y0+j][frame_num] = bgval; // Background
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MAKE_FRAME_TRIANGLE                           */
/*                                                                           */
/*  Fill frame "frame_num" with a triangle surrounded by background          */
/*  color.                                                                   */
/*                                                                           */
/*    sscale - degrees (of vis. angle) / pixel                               */
/*                                                                           */
/*  *** NOTE:  All values written into the frame are in [-1,1].              */
/*                                                                           */
/*****************************************************************************/
void make_frame_add_triangle(data,x0,xn,y0,yn,z0,zn,frame_num,sbot,sleft,sright,
			     xctr,yctr,sscale,ampl,bgval)
     float ***data;
     int x0,xn,y0,yn,z0,zn,frame_num;
     float sbot,sleft,sright;    // Length of sides (deg)
     float xctr,yctr;            // position of center
     float sscale,ampl,bgval;
{
  int i,j;
  int xc,yc;
  float x,y;
  float a,b,c,ta,theta_a,theta_b,ml,mr,h,wl;
  float ax,ay,bx,by,cx,cy,yleft,yright;

  //
  //  *** WYETH - TRIANGES THAT OVERHANG TO left or right are not handled ***
  //  *** WYETH - TRIANGES THAT OVERHANG TO left or right are not handled ***
  //  *** WYETH - TRIANGES THAT OVERHANG TO left or right are not handled ***
  //

  c = sbot;
  a = sleft;
  b = sright;

  //  Triangle
  //            cx,cy
  //             /C\
  //            /   \
  //         a /     \ b
  //          /       \
  //         /         \
  //  bx,by /B_________A\ ax,ay
  //              c
  //

  // Determine equations of lines
  theta_a = acos((b*b + c*c - a*a) / (2.0*b*c));     // Angle A, opposite to a
  theta_b = acos((a*a + c*c - b*b) / (2.0*a*c));     // Angle B, opposite to b
  //printf("    theta_a = %f\n",theta_a * 180.0/M_PI);
  //printf("    theta_b = %f\n",theta_b * 180.0/M_PI);
  ml = tan(theta_b);                                 // slope of left side
  mr = -tan(theta_a);                                // slope of right side
  //printf("    slope_left = %f\n",ml);
  //printf("    slope_right = %f\n",mr);
  h = a * sin(theta_b);                              // height of triangle
  //printf("    height = %f\n",h);
  wl = a * cos(theta_b);                             // width to left of peak

  ax =  sbot / 2.0;                                  // Point a, x-coord
  ay = -h    / 2.0;                                  // Point a, y-coord
  bx = -sbot / 2.0;                                  // Point b, x-coord
  by = -h    / 2.0;                                  // Point b, y-coord
  cx =  ax + wl;                                     // Point c, x-coord
  cy =  h    / 2.0;                                  // Point c, y-coord

  xc = (xn-1)/2;
  yc = (yn-1)/2;
  
  for(i=0;i<xn;i++){
    x = (float)(i-xc)*sscale;

    yleft  = by + ml * (x - bx);  // Y-coord on left side
    yright = ay + mr * (x - ax);  // Y-coord on right side

    for(j=0;j<yn;j++){
      y = (float)(j-yc)*sscale;
      
      if ((y >= ay) && (y <= yleft) && (y <= yright)){
	data[x0+i][y0+j][frame_num] = ampl;
      }
    }
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                              MAKE_FRAME_EDGE_01                           */
/*                                                                           */
/*  Wyeth - make stimuli w.r.t. collaboration w/ Bartlett.                   */
/*                                                                           */
/*****************************************************************************/
void make_frame_edge_01(data,x0,xn,y0,yn,z0,zn,frame_num,amp,sd,pseed)
     float ***data;
     int x0,xn,y0,yn,z0,zn,frame_num;
     float amp,sd;
     int *pseed;
{
  int i,j,j1,j2;
  int k,tmp_xn,tmp_yn,ymid,yoff,yoffr,cnt;
  float **gt,x,y,dx,s,ss,yran;

  k = (int)(3*sd);
  tmp_xn = 2*k - 1;
  tmp_yn = k;

  // Build Gaussian template
  gt = get_2d_farray(tmp_xn,tmp_yn);
  for(i=0;i<tmp_xn;i++){
    for(j=0;j<tmp_yn;j++){
      x = i;
      y = j;
      gt[i][j] = func_2d_gaussian(x,y,(float)k,0.0,sd,sd,0);  // max = 1
      if (gt[i][j] > 1.0)
	printf("%f\n",gt[i][j]);
    }
  }
  
  for(i=x0;i<(x0+xn);i++)
    for(j=y0;j<(y0+yn);j++)
      data[i][j][frame_num] = 0.0;
  
  //data[1][1][1] = -1.0;
  //data[1][1][2] =  1.0;

  ymid = (y0+yn)/2;

  //ss = -1.0 + 2.0*myrand_util_ran2(pseed);
  ss = 0.3 * myrand_util_ran2(pseed);

  dx = tmp_xn;
  cnt = 0;
  //for(i=x0+k;i<(x0+xn-k);i+=dx){
  for(i=x0+k;i<(x0+xn-k);i++){

    if (cnt == 0){
      yran = myrand_util_ran2(pseed);
      if (yran < 0.333)
	yoffr = ymid - 3;
      else if (yran > 0.667)
	yoffr = ymid + 3;
      else
	yoffr = ymid;
    }

    yoffr = ymid;  // WYETH - fix the yoff

    if (cnt == 1){
      yoff = ymid;
      s = ss;
    }else{
      yoff = yoffr;
      if (ss > 0)
	s = 1.0;
      else
	s = -1.0;
    }

    s = -0.3 + 0.6 * myrand_util_ran2(pseed);

    
    for(j1=0;j1<tmp_xn;j1++)
      for(j2=0;j2<tmp_yn;j2++)
	data[i+j1-k][yoff-j2][frame_num] += -s * gt[j1][j2];

    for(j1=0;j1<tmp_xn;j1++)
      for(j2=0;j2<tmp_yn;j2++)
	data[i+j1-k][yoff+j2+1][frame_num] += s * gt[j1][j2];

    cnt += 1;
  }
  
  free_2d_farray(gt,tmp_xn);

}
/**************************************-**************************************/
/*                                                                           */
/*                             MAKE_FRAME_SINUSOID                           */
/*                                                                           */
/*  Fill frame "frame_num" with a cosine centered in the middle of the       */
/*  frame.  For example, for 1..32, the center is 16.                        */
/*                                                                           */
/*    freq - cycles/degree of vis. angle                                     */
/*    theta - degrees (in the circle)                                        */
/*    sscale - degrees (of vis. angle) / pixel                               */
/*    phase - in degrees of a cycle [0,360)                                  */
/*                                                                           */
/*  *** NOTE:  All values written into the frame are in [-1,1].              */
/*                                                                           */
/*****************************************************************************/
void make_frame_sinusoid(data,x0,xn,y0,yn,z0,zn,frame_num,freq,theta,sscale,
			 phase,xc,yc)
     float ***data;
     int x0,xn,y0,yn,z0,zn,frame_num;
     float freq,theta,sscale,phase;
     float xc,yc;                       // Center (deg)
{
  int i,j;
  float x,y,nx,ny,twopif,ph;

  twopif = 2.0*M_PI*freq;
  ph = phase/180.0 * M_PI;
  nx = cos(theta/180.0 * M_PI);
  ny = sin(theta/180.0 * M_PI);

  //printf("twopi %f   freq %f\n",twopif,freq);

  // WYETH - if you want 'gobal phase' use this, but convert to deg:
  //xc = (xn-1)/2;
  //yc = (yn-1)/2;

  for(i=0;i<xn;i++){
    x = (float)i*sscale - xc;
    for(j=0;j<yn;j++){
      y = (float)j*sscale - yc;
      data[x0+i][y0+j][frame_num] = cos(twopif*(nx*x+ny*y) - ph);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MAKE_FRAME_SINUSOID_0                           */
/*                                                                           */
/*  *** Like MAKE_FRAME_SINUSOID above, but not for tensor.                  */
/*                                                                           */
/*  *** NOTE:  All values written into the frame are in [-1,1].              */
/*                                                                           */
/*****************************************************************************/
void make_frame_sinusoid_0(data,xn,yn,freq,theta,sscale,phase,xc,yc,ampl)
     float **data;  // [xn][yn] Draw into this 2D image
     int xn,yn;
     float freq;    // cycles/degree of vis. angle
     float theta;   // degrees (in the circle)
     float sscale;  // degrees (of vis. angle) / pixel
     float phase;   // in degrees of a cycle [0,360)
     float xc,yc;   // Center (deg)
     float ampl;    // Amplitude
{
  int i,j;
  float x,y,nx,ny,twopif,ph;

  twopif = 2.0*M_PI*freq;
  ph = phase/180.0 * M_PI;
  nx = cos(theta/180.0 * M_PI);
  ny = sin(theta/180.0 * M_PI);

  // WYETH - if you want 'gobal phase' use this, but convert to deg:
  //xc = (xn-1)/2;
  //yc = (yn-1)/2;

  for(i=0;i<xn;i++){
    x = (float)i*sscale - xc;
    for(j=0;j<yn;j++){
      y = (float)j*sscale - yc;
      data[i][j] = ampl * cos(twopif*(nx*x+ny*y) - ph);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          MAKE_FRAME_SINUSOID_SCALE                        */
/*                                                                           */
/*                                                                           */
/*  *** WYETH NOTE:  On Dec 21, 2015, I found that the comment below was     */
/*      incorrect - that the sine was not centered in the middle, but was    */
/*      in fact centered at pixel 0,0.  I updated the code so that this      */
/*      comment should now be true.                                          */
/*  Fill frame "frame_num" with a cosine centered in the middle of the       */
/*  frame.  For example, for 1..32, the center is 17.                        */
/*                                                                           */
/*    freq - cycles/degree of vis. angle                                     */
/*    theta - degrees (in the circle)                                        */
/*    sscale - degrees (of vis. angle) / pixel                               */
/*    phase - in degrees of a cycle [0,360)                                  */
/*                                                                           */
/*  NOTE:  All values written into the frame are in [-c,c].                  */
/*                                                                           */
/*****************************************************************************/
void make_frame_sinusoid_scale(data,x0,xn,y0,yn,z0,zn,frame_num,freq,theta,
			       sscale,phase,xc,yc,c,phorg)
     float ***data;
     int x0,xn,y0,yn,z0,zn,frame_num;
     float freq,theta,sscale,phase;
     float xc,yc;                       // Center (deg)
     float c;                           // Scale factor
     int   phorg;                       // Phase origin 0-LL corner, 1-center
{
  int i,j;
  float x,y,nx,ny,twopif,ph;

  int cxi,cyi;

  if (phorg == 1){   // Set phase origin of sinusoid to (xn/2, yn/2)
    cxi = xn/2;
    cyi = yn/2;
  }else{
    cxi = cyi = 0;   // This was the only way before Dec 21, 2015
  }

  twopif = 2.0*M_PI*freq;
  ph = phase/180.0 * M_PI;
  nx = cos(theta/180.0 * M_PI);
  ny = sin(theta/180.0 * M_PI);

  // I believe subtracting 'xc' and 'yc' is important so that the grating will
  // maintain the same phase as it is translated in the image.

  for(i=0;i<xn;i++){
    //x = (float)i*sscale - xc;    // OLD
    x = (float)(i-cxi)*sscale - xc;  // phase constant at center (xn/2,yn/2)
    for(j=0;j<yn;j++){
      //y = (float)j*sscale - yc;
      y = (float)(j-cyi)*sscale - yc; // phase constant at center (xn/2,yn/2)
      data[x0+i][y0+j][frame_num] = c * cos(twopif*(nx*x+ny*y) - ph);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             VDP_TILT_RECTANGLE                            */
/*                                                                           */
/*  Draw a solid filled rectangle.                                           */
/*                                                                           */
/*                                                                           */
/*                                  |                                        */
/*                                  |                                        */
/*                                  |         .(x1,y1)                       */
/*                                  |     .  / .                             */
/*      *             Eqn. for side | .     /   .                            */
/*       *              y=m1*x+b1 . |      /     .                           */
/*        *       d           .     |     /       .                          */
/*         *--------------.         |   r/         . y=m2*x+b2               */
/*          *         .             |   /           .                        */
/*    here,  *    .                 |  /  phi   .    .                       */
/*    ytop=y4 .    .                | /     .         .                      */
/*     (x4,y4) .    .               |/  .   theta      .                     */
/*       _______.____.______________|___________________.___________         */
/*               .    .             |                    .                   */
/*                .    .            |                     .(x2,y2)           */
/*                 .hwid.th         |                 .  here, ybot=y2       */
/*                  .____.          |             .                          */
/*                   .    .         |         .                              */
/*          y=m2*x+b4 .    .        |     .                                  */
/*                     .    .       | .    y=m1*x+b3                         */
/*                      .    .    . |                                        */
/*                       .    .     |                                        */
/*                        .         |                                        */
/*                  (x3,y3)         |                                        */
/*                                  |                                        */
/*                                  |                                        */
/*                                                                           */
/*****************************************************************************/
void vdp_tilt_rectangle(xcenter,ycenter,width,height,angle,k)
     int xcenter,ycenter,width,height,angle,k;
{
  int ip; // Pixel ip;
  short ncols,row;
  int jp; // register Pixel jp;
  register short col;

  int i;
  int x1,y1,x2,y2,x3,y3,x4,y4,ytop,ybot,left,right,dx1,dy1,dx2,dy2;
  float fh,fw,theta,phi,r,m1,m2,b1,b2,b3,b4,ml,bl,mr,br;
  int ycenter2,xmax,nn;
  int ROWINC;

  /// WYETH - THIS COPIED FROM OLD CODE
  /// WYETH - THIS COPIED FROM OLD CODE
  /// WYETH - THIS COPIED FROM OLD CODE
  /// WYETH - THIS COPIED FROM OLD CODE

  ROWINC = 128; // WYETH THIS IS ARBITRAY, MUST FIX
  ROWINC = 128; // WYETH THIS IS ARBITRAY, MUST FIX

  /***
  if (angle >= 180)
    angle -= 180;
  if (angle >= 90){ // Swap the width and height.
    angle -= 90;
    i = width;
    width = height;
    height = i;
  }

  if (angle == 0){
    vdpfrect(ycenter-height/2,xcenter-width/2,height,width,k);
    return;
  }

  xmax = getncols();
  ycenter2 = 2*ycenter;
  theta = (float)angle*M_PI/180.0;
  fh = (float)height; fw = (float)width;

  phi = atan(fh/fw);
  r = sqrt(fw*fw/4.0 + fh*fh/4.0);
  x1 = (int)(0.5 + r*cos(theta + phi)) + xcenter;
  y1 = (int)(0.5 + r*sin(theta + phi)) + ycenter;
  dx1 = (int)(0.5 + fw*cos(theta)); dy1 = (int)(0.5 + fw*sin(theta));
  dx2 = (int)(0.5 + fh*cos(theta)); dy2 = (int)(0.5 + fh*sin(theta));
  x2 = x1+dy2; y2 = y1-dx2; x3 = x2-dx1; y3 = y2-dy1; x4 = x1-dx1; y4 = y1-dy1;
  m1 = (float)dy1/(float)dx1;  m2 = -(float)dx2/(float)dy2;
  b1 = (float)y1 - m1*(float)x1;  b2 = (float)y2 - m2*(float)x2;
  b3 = (float)y3 - m1*(float)x3;  b4 = (float)y3 - m2*(float)x3;

  if (y2 > y4){
    ytop = y2; ybot = y4; ml = m1; bl = b1; mr = m1; br = b3;
  }else{
    ytop = y4; ybot = y2; ml = m2; bl = b4; mr = m2; br = b2;
  }

  ip = vdppaddr(ycenter2-y1,0);

  nn = y1-ytop; //for(i=y1;i>ytop;i--){
  for(row=0;row<nn;row++){
    i = y1 - row;
    left = (int)(0.5 + ((float)i - b1)/m1);
    right = (int)(0.5 + ((float)i - b2)/m2);
    if (left < 0) left = 0;
    if (right >= xmax) right = xmax - 1;

    col = left; // Starting screen pixel column for this row.
    jp = ip + col; // Starting address in VDP memory for this row.
    ncols = right-left;
    while(ncols--)
      psetainc(jp,k); // set value at and increment "jp"
    ip += ROWINC;
  }

  nn = ytop-ybot; //for(i=ytop;i>ybot;i--){
  for(row=0;row<nn;row++){
    i = ytop - row;
    left = (int)(0.5 + ((float)i - bl)/ml);
    right = (int)(0.5 + ((float)i - br)/mr);
    if (left < 0) left = 0;
    if (right >= xmax) right = xmax - 1;

    col = left; // Starting screen pixel column for this row.
    jp = ip + col; // Starting address in VDP memory for this row.
    ncols = right-left;
    while(ncols--)
      psetainc(jp,k);
    ip += ROWINC;
  }

  nn = ybot-y3+1; // for(i=ybot;i>=y3;i--){
  for(row=0;row<nn;row++){
    i = ybot - row;
    left = (int)(0.5 + ((float)i - b4)/m2);
    right = (int)(0.5 + ((float)i - b3)/m1);
    if (left < 0) left = 0;
    if (right >= xmax) right = xmax - 1;

    col = left; // Starting screen pixel column for this row.
    jp = ip + col; // Starting address in VDP memory for this row.
    ncols = right-left;
    while(ncols--)
      psetainc(jp,k);
    ip += ROWINC;
  }
  ***/
}
/**************************************-**************************************/
/*                                                                           */
/*                            MAKE_FRAME_BAR_ADD_0                           */
/*                                                                           */
/*  Add a bar into a zero-background frame                                   */
/*                                                                           */
/*****************************************************************************/
void make_frame_bar_add_0(data,xn,yn,cbx,cby,theta,len,wid,sdaa,sscale,ampl,
			  wflag)
     float **data;       // [xn][yn] image data
     int xn,yn;          // size of image
     float cbx,cby;      // Bar center (deg)
     float theta;        // Orientation (deg)
     float len,wid;      // Length and width of bar (deg)
     float sdaa;         // SD for anti-aliasing (pix)
     float sscale;       // (deg/pix)
     float ampl;         // can be positive or negative
     int   wflag;        // 0-add, 1-max, -1-min
{
  int i,j;
  int xc,yc;
  float x,y,xx,yy,f,a,b,c,d,thetar,lo2,wo2,db,dbx,dby,sd,sd4,newval;

  /*
  printf("cbx,cby = %f %f\n",cbx,cby);
  printf("len,wid = %f %f\n",len,wid);
  printf("sscale = %f\n",sscale);
  printf("theta = %f\n",theta);
  printf("ampl = %f\n",ampl);
  exit(0);
  */

  lo2 = len/2.0;
  wo2 = wid/2.0;
  sd = sdaa * sscale;  // (deg)
  sd4 = sd * 4.0;      // (deg)
  
  while(theta >= 180.0)
    theta -= 180.0;
  while(theta < 0.0)
    theta += 180.0;

  thetar = -theta * M_PI/180.0;
  a =  cos(thetar); b = sin(thetar); // rotation matrix
  c = -sin(thetar); d = cos(thetar);

  xc = (xn-1)/2; // Center of field, integer (pixels)
  yc = (yn-1)/2; // Note, integer arithmetic !

  for(i=0;i<xn;i++){
    xx = (float)(i-xc)*sscale - cbx;  // relative to bar center (deg)
    for(j=0;j<yn;j++){
      yy = (float)(j-yc)*sscale - cby; // relative to bar center (deg)

      x = a*xx + c*yy;  // rotate back
      y = b*xx + d*yy;

      // Compute distance from bar edge, 'db'
      if (x > 0.0)
	dbx = wo2 - x;  // +/- indicates inside/outside
      else
	dbx = x + wo2;

      if (y > 0.0)
	dby = lo2 - y;  // +/- indicates inside/outside
      else
	dby = y + lo2;

      if (dbx < 0.0){   // x outside
	if (dby < 0.0){
	  db = -sqrt(dbx*dbx + dby*dby);  // Both outside
	}else{
	  db = dbx;
	}
      }else{ // x inside
	if (dby < 0.0){  // y outside
	  db = dby;
	}else{ // Both inside - take the smallest distance
	  if (dby < dbx)
	    db = dby;
	  else
	    db = dbx;
	}
      }

      if (sdaa > 0.0){
	if (db > sd)
	  f = 1.0; // 1.0;  // far inside
	else if (db < -sd)
	  f = 0.0;  // far outside
	else // Use sinusoidal sigmoid
	  f = 0.5 + 0.5*sin(db/sd * 0.25*M_PI);
      }else{
	if (db >= 0.0)
	  f = 1.0; // inside
	else
	  f = 0.0; // outside
      }

      newval = f * ampl;

      if (wflag == 1){
	//
	//  Write if the new value is larger
	//
	if (newval > data[i][j])
	  data[i][j] = newval;
      }else if (wflag == -1){
	//
	//  Write if the new value is smaller
	//
	if (newval < data[i][j])
	  data[i][j] = newval;
      }else{
	//
	//  Add the bar
	//
	data[i][j] += newval;
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               MAKE_FRAME_BAR                              */
/*                                                                           */
/*  Fill frame "frame_num" with a bar.                                       */
/*  Note, the spatial center is taken to be n/2, e.g., for 1..32, ctr = 16.  */
/*                                                                           */
/*****************************************************************************/
void make_frame_bar(data,x0,xn,y0,yn,z0,zn,frame_num,cbx,cby,theta,len,wid,
		    sdaa,sscale,ampl,bgval,flag)
     float ***data;
     int x0,xn,y0,yn,z0,zn,frame_num;
     float cbx,cby,theta,len,wid,sdaa,sscale,ampl,bgval;
     int flag;  // 0-Bar and background, 1-write bar only
{
  int i,j;
  int xc,yc;
  float x,y,xx,yy,f,a,b,c,d,thetar,lo2,wo2,db,dbx,dby,sd,sd4,newval;

  //printf("flag = %d\n",flag);

  /*
    printf("cbx, cby = %f %f\n",cbx,cby);
    printf("theta len wid = %f %f %f\n",theta,len,wid);
    printf("sdaa sscale ampl bgval = %f %f %f %f\n",sdaa,sscale,ampl,bgval);
  */

  lo2 = len/2.0;
  wo2 = wid/2.0;
  sd = sdaa * sscale;  // (deg)
  sd4 = sd * 4.0;      // (deg)
  
  while(theta >= 180.0)
    theta -= 180.0;
  while(theta < 0.0)
    theta += 180.0;

  thetar = -theta * M_PI/180.0;
  a =  cos(thetar); b = sin(thetar); // rotation matrix
  c = -sin(thetar); d = cos(thetar);

  xc = (xn-1)/2; // Center of field, integer (pixels)
  yc = (yn-1)/2; // Note, integer arithmetic !

  for(i=0;i<xn;i++){
    xx = (float)(i-xc)*sscale - cbx;  // relative to bar center (deg)
    for(j=0;j<yn;j++){
      yy = (float)(j-yc)*sscale - cby; // relative to bar center (deg)

      x = a*xx + c*yy;  // rotate back
      y = b*xx + d*yy;

      // Compute distance from bar, 'db'
      if (x > wo2)
	dbx = x - wo2;
      else if (x < -wo2)
	dbx = -wo2 - x;
      else
	dbx = 0.0;

      if (y > lo2)
	dby = y - lo2;
      else if (y < -lo2)
	dby = -lo2 - y;
      else
	dby = 0.0;

      db = sqrt(dbx*dbx + dby*dby);

      // printf("db sd4 %f %f   x,y = %f %f\n",db,sd4,x,y);

      if (db > sd4)
	f = 0.0;     // Greater than 4 SDs away
      else if (db > 0.0)
	f = func_gaussian_one(db,0.0,sd);
      else
	f = 1.0;

      newval = bgval + f * (ampl-bgval);

      if (flag == 0) // Bar and background
	data[x0+i][y0+j][frame_num] = newval;
      else if (f > 0.0){ // Add bar, do not change background

	// WYETH - this used to be a BIG HACK, only for light on dark
	// but now works for dark on light as well.
	if (ampl >= bgval){
	  if (data[x0+i][y0+j][frame_num] < newval)
	    data[x0+i][y0+j][frame_num] = newval;
	}else{
	  if (data[x0+i][y0+j][frame_num] > newval)
	    data[x0+i][y0+j][frame_num] = newval;
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MAKE_FRAME_BAR_0                             */
/*                                                                           */
/*  Like 'MAKE_FRAME_BAR' above, except this works in a regular 2D array,    */
/*  as opposed to a frame w/i a tensor.                                      */
/*                                                                           */
/*****************************************************************************/
void make_frame_bar_0(data,xn,yn,cbx,cby,theta,len,wid,
		      sdaa,sscale,ampl,bgval,flag)
     float **data;
     int xn,yn;
     float cbx,cby,theta,len,wid,sdaa,sscale,ampl,bgval;
     int flag;  // 0-Bar and background, 1-write bar only
{
  int i,j;
  int xc,yc;
  float x,y,xx,yy,f,a,b,c,d,thetar,lo2,wo2,db,dbx,dby,sd,sd4,newval;

  //printf("flag = %d\n",flag);

  /*
    printf("cbx, cby = %f %f\n",cbx,cby);
    printf("theta len wid = %f %f %f\n",theta,len,wid);
    printf("sdaa sscale ampl bgval = %f %f %f %f\n",sdaa,sscale,ampl,bgval);
  */

  lo2 = len/2.0;
  wo2 = wid/2.0;
  sd = sdaa * sscale;  // (deg)
  sd4 = sd * 4.0;      // (deg)
  
  while(theta >= 180.0)
    theta -= 180.0;
  while(theta < 0.0)
    theta += 180.0;

  thetar = -theta * M_PI/180.0;
  a =  cos(thetar); b = sin(thetar); // rotation matrix
  c = -sin(thetar); d = cos(thetar);

  xc = (xn-1)/2; // Center of field, integer (pixels)
  yc = (yn-1)/2; // Note, integer arithmetic !

  for(i=0;i<xn;i++){
    xx = (float)(i-xc)*sscale - cbx;  // relative to bar center (deg)
    for(j=0;j<yn;j++){
      yy = (float)(j-yc)*sscale - cby; // relative to bar center (deg)

      x = a*xx + c*yy;  // rotate back
      y = b*xx + d*yy;

      // Compute distance from bar, 'db'
      if (x > wo2)
	dbx = x - wo2;
      else if (x < -wo2)
	dbx = -wo2 - x;
      else
	dbx = 0.0;

      if (y > lo2)
	dby = y - lo2;
      else if (y < -lo2)
	dby = -lo2 - y;
      else
	dby = 0.0;

      db = sqrt(dbx*dbx + dby*dby);

      // printf("db sd4 %f %f   x,y = %f %f\n",db,sd4,x,y);

      if (db > sd4)
	f = 0.0;     // Greater than 4 SDs away
      else if (db > 0.0)
	f = func_gaussian_one(db,0.0,sd);
      else
	f = 1.0;

      newval = bgval + f * (ampl-bgval);

      if (flag == 0) // Bar and background
	data[i][j] = newval;
      else if (f > 0.0){ // Add bar, do not change background

	if (ampl >= bgval){
	  if (data[i][j] < newval)
	    data[i][j] = newval;
	}else{
	  if (data[i][j] > newval)
	    data[i][j] = newval;
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               MAKE_FRAME_BARS                             */
/*                                                                           */
/*  Fill frame "frame_num" with a bar.                                       */
/*  Note, the spatial center is taken to be n/2, e.g., for 1..32, ctr = 16.  */
/*                                                                           */
/*****************************************************************************/
void make_frame_bars(data,x0,xn,y0,yn,z0,zn,frame_num,nbar,bw,bl,theta,
		     bar_val,bgval)
     float ***data;
     int x0,xn,y0,yn,z0,zn;
     int frame_num;    // frame index
     int nbar;         // Number of bars
     float bw,bl;      // Bar width, bar length (pix)
     float theta;      // orientation, 0...360 deg
     float *bar_val;   // bar luminance [nbar]
     float bgval;      // background value
{
  int i,j;
  int w,l,bi,xx0,xx1,yy0,yy1,x1,y1;

  /*
    printf("bw, bl = %f %f\n",bw,bl);
    printf("theta = %f\n",theta);
    printf("nbar = %d\n",nbar);
  */

  if (theta != 0){
    exit_error("MAKE_FRAME_BARS","Direction must be zero, others not imp'd");
  }

  x1 = x0 + xn - 1;  // Highest valid x-index
  y1 = y0 + yn - 1;  // Highest valid y-index

  w = my_rint(bw);  // Bar width pix
  l = my_rint(bl);  // Bar length pix

  // Draw background
  for(i=x0;i<(x0+xn);i++){
    for(j=y0;j<(y0+yn);j++){
      data[i][j][frame_num] = bgval;
    }
  }

  yy0 = y0 + yn/2 - l/2;
  yy1 = y0 + yy0 + l - 1;

  if (yy0 < y0)
    yy0 = y0;
  if (yy1 > y1)
    yy1 = y1;

  //bar_val[0] = 0.7;

  xx0 = x0 + xn/2 - nbar*w/2;
  for(bi=0;bi<nbar;bi++){
    xx1 = xx0 + w - 1;
    for(i=xx0;i<=xx1;i++){
      if ((i >= x0) && (i <= x1)){
	for(j=yy0;j<=yy1;j++){
	  data[i][j][frame_num] = bar_val[bi];
	}
      }
    }
    xx0 += w;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MAKE_FRAME_BARS_TMPLT                           */
/*                                                                           */
/*  Fill frame "frame_num" with a bar.                                       */
/*                                                                           */
/*****************************************************************************/
float **make_frame_bars_tmplt(bartemp,xn,yn,nbar,bar_val,bgval)
     float **bartemp;  // Template for bar pattern [xn][yn]
     int xn,yn;        // Size of frame and template
     int nbar;         // Number of bars
     float *bar_val;   // bar luminance [nbar]
     float bgval;      // background value
{
  int i,j;
  int bi;
  float **d,tv,v0,v1,rem;

  d = get_2d_farray(xn,yn);

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      tv = bartemp[i][j];    // Template value
      bi = (int)tv;          // bar index
      rem = tv - (float)bi;  // Remainder, beyond integer value
      if (bi < 1){
	v0 = bgval;
	v1 = bar_val[0];  // first bar value
      }else if (bi > nbar){  // Fully beyond last bar, rem=0, does this occur?
	v0 = bgval;
	v1 = bgval;
      }else if (bi == nbar){
	v0 = bar_val[nbar-1];  // last bar value
	v1 = bgval;
      }else{
	v0 = bar_val[bi-1];   // Convert from [1..nbar] to [0..nbar-1]
	v1 = bar_val[bi];
      }

      if (rem == 0.0){
	d[i][j] = v0;  // WYETH - does this save time?
      }else{
	d[i][j] = (1.0 - rem) * v0 + rem*v1;  // Weighted average
      }
    }
  }

  return d;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MAKE_FRAME_BARPATMO                            */
/*                                                                           */
/*  Draw a bar pattern for a random bar motion stimulus.                     */
/*                                                                           */
/*****************************************************************************/
float **make_frame_barpatmo(xn,yn,x0,y0,nbar,bwpix,blpix,seq,seqn,si,pi,
			    altsign,val1,val2,bgval,idir)
     int xn,yn;        // Size of frame and template
     int x0,y0;        // Lower left corner of bar pattern
     int nbar;         // Number of bars
     int bwpix;        // bar width (pix)
     int blpix;        // bar length (pix)
     int *seq;         // sequence of bar values [seqn]
     int seqn;         // length of sequence
     int si;           // sequence index
     int pi;           // pixel index
     int altsign;      // 1 for contrast, -1 for inverse contrast
     float val1,val2;  // light and dark values [0..1]
     float bgval;      // background value [0..1]
     int idir;         // integer value for direction, 0 or 180 deg
{
  int i,j;
  int irev,x1,y1,seqval;
  float **d,bval;

  d = get_const_2d_farray(xn,yn,bgval);  // Values will range from 0..1

  x1 = x0 + nbar*bwpix;
  y1 = y0 + blpix;

  //  Do not draw above or below the stimulus grid
  if (y0 < 0)
    y0 = 0;
  if (y1 > yn)  // NOTE, only values <y1 are used in the loop below
    y1 = yn;

  for(i=x0;i<x1;i++){
    // Get bar color at this x-pos
    if (si >= seqn)
      exit_error("MAKE_FRAME_BARPATMO","Sequence index too large");

    seqval = seq[si];  // Sequence value for this bar

    //  Determine bar luminance, 'bval'
    if (altsign == 1){  // Invert contrast
      if (seqval == -1)
	bval = val1;
      else if (seqval == 1)
	bval = val2;
      else
	bval = bgval;
    }else{
      if (seqval == 1)
	bval = val1;
      else if (seqval == -1)
	bval = val2;
      else
	bval = bgval;
    }

    if ((i >= 0) && (i < xn)){  // Only draw w/i stimulus grid
      if (idir == 180){
	for(j=y0;j<y1;j++)
	  d[i][j] = bval;
      }else if (idir == 0){
	irev = x1-1 - (i-x0);
	for(j=y0;j<y1;j++)
	  d[irev][j] = bval;
      }
    }

    pi += 1;
    if (pi >= bwpix){
      pi = 0;
      si += 1;  // Move on to next bar in sequence
    }
  }

  return d;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MAKE_FRAME_ADD_GPATCH_SINE                        */
/*                                                                           */
/*  Add a sine patch function to the frame.                                  */
/*  Originally devloped for Kumbhani et al. (2015) stimulus.                 */ 
/*                                                                           */
/*****************************************************************************/
void make_frame_add_gpatch_sine(data,xn,yn,gp,sscale,phoff)
     float **data;
     int xn,yn;               // (pix)
     struct stmh_gpatch *gp;  // Gabor params
     float sscale;            // spatial scale (deg/pix)
     float phoff;             // phase offset (deg);
{
  int i,j;
  int xc,yc,i1,i2,j1,j2,xx,yy;
  float xci,yci,sf,phase,con,diam,m_crit,m_fact,m_len,radmarg;
  float x,y,nx,ny,ph,twopif,r1s;

  diam  = gp->sd_orth;  // Hack: 'sd_orth' has diameter value
  m_len = gp->sd_par;   // Hack: 'sd_par'  has taper margin size

  //
  //  Values for computing the raised cosine taper of outer edge
  //
  m_crit = diam/2.0 - m_len;
  if (m_crit < 0.0)
    exit_error("MAKE_FRAME_ADD_GPATCH_SINE","Taper margin (sd_par) too big");
  m_fact = M_PI / m_len;  // Scale margin by this value
  

  xci = (xn-1)/2;  // center of frame (pix)
  yci = (yn-1)/2;

  xci += gp->x0 / sscale;  // center of Gabor (pix)
  yci += gp->y0 / sscale;

  phase  = gp->ph0 + phoff;       // (deg)

  con = gp->amp;

  r1s = diam*diam/4.0;  // Pre-compute radius squared

  twopif = 2.0*M_PI*gp->sf;
  ph = phase/180.0 * M_PI;
  nx = cos(gp->dir/180.0 * M_PI);
  ny = sin(gp->dir/180.0 * M_PI);
  xc = xci;
  yc = yci;

  xx = my_rint(diam/sscale) + 1; // Max width of patch (pix)
  yy = xx;                       // Assume circular patch

  // Define the rectangular region which covers the circle and is within
  //   the frame boundary
  i1 = xc - xx/2 - 1; // Left x value
  if (i1 < 0) i1 = 0;
  i2 = xc + xx/2 + 1; // Right x value
  if (i2 >= xn) i2 = xn-1;

  j1 = yc - yy/2 - 1; // Min y value
  if (j1 < 0) j1 = 0;
  j2 = yc + yy/2 + 1; // Max y value
  if (j2 >= yn) j2 = yn-1;

  //printf(" X: i1,i2 = %d .. %d\n",i1,i2);
  //printf(" Y: j1,j2 = %d .. %d\n",j1,j2);


  // Draw the sine wave in the circle
  for(i=i1;i<=i2;i++){
    x = (float)(i-xc)*sscale;
    for(j=j1;j<=j2;j++){
      y = (float)(j-yc)*sscale;
      if ((x*x + y*y) < r1s){ // Inside patch
	radmarg = sqrt(x*x + y*y) - m_crit;
	if (radmarg > 0.0){ // taper the outer edge
	  data[i][j] += con*cos(twopif*(nx*x+ny*y) - ph) *
	    (0.5 + 0.5*cos(radmarg * m_fact));
	  //data[i][j] = 1.0;
	}else{
	  data[i][j] += con*cos(twopif*(nx*x+ny*y) - ph);
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         MAKE_FRAME_ADD_GPATCH_GABOR                       */
/*                                                                           */
/*  Add a Gabor function to the frame.                                       */
/*                                                                           */
/*****************************************************************************/
void make_frame_add_gpatch_gabor(data,xn,yn,gp,sscale,phoff)
     float **data;
     int xn,yn;               // (pix)
     struct stmh_gpatch *gp;  // Gabor params
     float sscale;            // spatial scale (deg/pix)
     float phoff;             // phase offset (deg);
{
  int i,j;
  float xci,yci,sdp,sdo,sf,ph;
  float x,y;
  float **f;

  xci = (xn-1)/2;  // center of frame (pix)
  yci = (yn-1)/2;

  xci += gp->x0 / sscale;  // center of Gabor (pix)
  yci += gp->y0 / sscale;

  sdp = gp->sd_par / sscale;   // (pix)
  sdo = gp->sd_orth / sscale;  // (pix)
  sf  = gp->sf * sscale;       // (cyc/pix)
  ph  = gp->ph0 + phoff;       // (deg)

  // 'ph' is negated here because the phase param is *added*
  f = gabor_2d_space_raw(xn,yn,xci,yci,sdo,sdp,sf,gp->dir,-ph);

  //printf("_____(stim_util.c)____GP_AMP = %f\n",gp->amp);

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      data[i][j] += gp->amp * f[i][j];
      //if (data[i][j] > 1.00001)
      //printf(" data[i][j] = %f\n",data[i][j]);
    }
  }

  free_2d_farray(f,xn);
}
/**************************************-**************************************/
/*                                                                           */
/*                         MAKE_FRAME_ADD_GPATCH_NOISE                       */
/*                                                                           */
/*  Add a noisy Gaussian shaped patch to the frame.                          */
/*                                                                           */
/*****************************************************************************/
void make_frame_add_gpatch_noise(data,xn,yn,gp,sscale,phoff)
     float **data;
     int xn,yn;               // (pix)
     struct stmh_gpatch *gp;  // Gabor params
     float sscale;            // spatial scale (deg/pix)
     float phoff;             // phase offset (deg);
{
  int i,j;
  float xci,yci,sdp,sdo,sf,ph;
  float x,y;
  float **f;

  xci = (xn-1)/2;  // center of frame (pix)
  yci = (yn-1)/2;

  xci += gp->x0 / sscale;  // center of Gabor (pix)
  yci += gp->y0 / sscale;

  sdp = gp->sd_par / sscale;   // (pix)
  sdo = gp->sd_orth / sscale;  // (pix)
  sf  = gp->sf * sscale;       // (cyc/pix)
  ph  = gp->ph0 + phoff;       // (deg)

  // 'ph' is negated here because the phase param is *added*
  //f = gabor_2d_space_raw(xn,yn,xci,yci,sdo,sdp,sf,gp->dir,-ph);
  f = gauss_2d_noise_raw(xn,yn,xci,yci,sdo,sdp,gp->dir,gp->pixw,gp->seed);


  //printf("_____(stim_util.c)____GP_AMP = %f\n",gp->amp);

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      data[i][j] += gp->amp * f[i][j];
      //if (data[i][j] > 1.00001)
      //printf(" data[i][j] = %f\n",data[i][j]);
    }
  }

  free_2d_farray(f,xn);
}
/**************************************-**************************************/
/*                                                                           */
/*                           MAKE_FRAME_ADD_2D_GAUSS                         */
/*                                                                           */
/*  Add a 2D Gaussian to the frame.                                          */
/*                                                                           */
/*****************************************************************************/
void make_frame_add_2d_gauss(data,xn,yn,amp,sd,x0,y0,sscale)
     float **data;
     int xn,yn;       // (pix)
     float amp;       // amplituded of Gaussian
     float sd,x0,y0;  // (deg)
     float sscale;    // spatial scale (deg/pix)

{
  int i,j;
  float xci,yci;
  float x,y;

  xci = (xn-1)/2;
  yci = (yn-1)/2;

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      x = (xci - (float)i)*sscale - x0;
      y = (yci - (float)j)*sscale - y0;
      data[i][j] += exp(-(x*x + y*y)/(0.5*sd*sd));
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MAKE_FRAME_ADD_RECT_PIX                         */
/*                                                                           */
/*  Add a value to the frame in a rectangular region (specified in pixels).  */
/*                                                                           */
/*****************************************************************************/
void make_frame_add_rect_pix(data,x0,xn,y0,yn,z0,zn,frame_num,x,y,w,h,val)
     float ***data;
     int x0,xn,y0,yn,z0,zn,frame_num;
     int x,y;    // Lower left corner of box
     int w,h;    // width and height of box
     float val;  // value to add in box
{
  int i,j;
  int x1,x2,y1,y2;

  x1 = x;
  if (x1<x0)
    x1 = x0;

  x2 = x+w-1;
  if (x2 >= (x0+xn))
    x2 = (x0+xn-1);

  y1 = y;
  if (y1<y0)
    y1 = y0;

  y2 = y+h-1;
  if (y2 >= (y0+yn))
    y2 = (y0+yn-1);

  for(i=x1;i<=x2;i++){
    for(j=y1;j<=y2;j++){
      data[i][j][frame_num] += val;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          MAKE_FRAME_ADD_RECT_PIX_0                        */
/*                                                                           */
/*  Add a value to the frame in a rectangular region (specified in pixels).  */
/*                                                                           */
/*****************************************************************************/
void make_frame_add_rect_pix_0(data,xn,yn,x,y,w,h,val)
     float **data;
     int xn,yn;
     int x,y;    // Lower left corner of box
     int w,h;    // width and height of box
     float val;  // value to add in box
{
  int i,j;
  int x1,x2,y1,y2;

  x1 = x;
  if (x1 < 0)
    x1 = 0;

  x2 = x+w-1;
  if (x2 >= xn)
    x2 = xn-1;

  y1 = y;
  if (y1 < 0)
    y1 = 0;

  y2 = y+h-1;
  if (y2 >= yn)
    y2 = yn-1;

  for(i=x1;i<=x2;i++){
    for(j=y1;j<=y2;j++){
      data[i][j] += val;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MAKE_FRAME_ADD_SINE_BUMP                        */
/*                                                                           */
/*  Add a sinusoidal bump, which goes from 0 to 'val' to 0.                  */
/*                                                                           */
/*  NOTE:  coordinates are for CENTER and are FLOAT here.                    */
/*                                                                           */
/*****************************************************************************/
void make_frame_add_sine_bump(data,x0,xn,y0,yn,z0,zn,frame_num,x,y,diam,val)
     float ***data;
     int x0,xn,y0,yn,z0,zn,frame_num;
     float x,y;   // Center of bump
     float diam;  // diam of bump (1 cycle of sine wave)
     float val;   // Amplitude of bump (bump goes from 0 to val back to 0)
{
  int i,j;
  int x1,x2,y1,y2;
  float c,do2,w,dx,dy,dist;

  // Determine the full extent of non-zero pixels for this bump

  x1 = (int)(x - diam/2.0);      // Round down, thus, to the left
  x2 = (int)(x + diam/2.0 + 1);  // Round up, approximately, thus to right
  
  y1 = (int)(y - diam/2.0);      // Round down, thus, to the bottom
  y2 = (int)(y + diam/2.0 + 1);  // Round up, approximately, thus to top

  if (x1 < x0)
    x1 = x0;

  if (x2 >= (x0+xn))
    x2 = x0 + xn - 1;

  if (y1 < y0)
    y1 = y0;
  
  if (y2 >= (y0+yn))
    y2 = y0 + yn - 1;

  //printf("%4d   x,y  %f %f  diam %f\n",frame_num,x,y,diam);

  do2 = diam/2.0;
  c = M_PI/do2;

  for(i=x1;i<=x2;i++){
    for(j=y1;j<=y2;j++){
      dx = (float)i - x;
      dy = (float)j - y;
      dist = sqrt(dx*dx + dy*dy);
      if (dist > do2)
	w = 0.0;
      else
	w = 0.5 + 0.5*cos(dist * c);

      data[i][j][frame_num] += val * w;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          MAKE_FRAME_ADD_SINE_BUMP_0                       */
/*                                                                           */
/*  Add a sinusoidal bump, which goes from 0 to 'val' to 0.                  */
/*                                                                           */
/*  NOTE:  coordinates are for CENTER and are FLOAT here.                    */
/*                                                                           */
/*****************************************************************************/
void make_frame_add_sine_bump_0(data,xn,yn,x,y,diam,val)
     float **data;
     int xn,yn;
     float x,y;   // Center of bump
     float diam;  // diam of bump (1 cycle of sine wave)
     float val;   // Amplitude of bump (bump goes from 0 to val back to 0)
{
  int i,j;
  int x1,x2,y1,y2;
  float c,do2,w,dx,dy,dist;

  // Determine the full extent of non-zero pixels for this bump

  x1 = (int)(x - diam/2.0);      // Round down, thus, to the left
  x2 = (int)(x + diam/2.0 + 1);  // Round up, approximately, thus to right
  
  y1 = (int)(y - diam/2.0);      // Round down, thus, to the bottom
  y2 = (int)(y + diam/2.0 + 1);  // Round up, approximately, thus to top

  if (x1 < 0)
    x1 = 0;

  if (x2 >= xn)
    x2 = xn-1;

  if (y1 < 0)
    y1 = 0;
  
  if (y2 >= yn)
    y2 = yn-1;

  do2 = diam/2.0;
  c = M_PI/do2;

  for(i=x1;i<=x2;i++){
    for(j=y1;j<=y2;j++){
      dx = (float)i - x;
      dy = (float)j - y;
      dist = sqrt(dx*dx + dy*dy);
      if (dist > do2)
	w = 0.0;
      else
	w = 0.5 + 0.5*cos(dist * c);

      data[i][j] += val * w;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MAKE_FRAME_BINARIZE                           */
/*                                                                           */
/*  Set to a if >= threshold, b otherwise.                                   */
/*                                                                           */
/*****************************************************************************/
void make_frame_binarize(data,x0,xn,y0,yn,z0,zn,frame_num,thresh,a,b)
     float ***data;
     int x0,xn,y0,yn,z0,zn,frame_num;
     float thresh,a,b;
{
  int i,j;

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      if (data[x0+i][y0+j][frame_num] >= thresh)
	data[x0+i][y0+j][frame_num] = a;
      else
	data[x0+i][y0+j][frame_num] = b;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MAKE_FRAME_BINARIZE_0                          */
/*                                                                           */
/*  Set to a if >= threshold, b otherwise.                                   */
/*                                                                           */
/*****************************************************************************/
void make_frame_binarize_0(data,xn,yn,thresh,a,b)
     float **data;
     int xn,yn;
     float thresh,a,b;
{
  int i,j;

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      if (data[i][j] >= thresh)
	data[i][j] = a;
      else
	data[i][j] = b;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MAKE_FRAME_ADD_DOTS                            */
/*                                                                           */
/*  Draw dots into the frame.                                                */
/*                                                                           */
/*****************************************************************************/
void make_frame_add_dots(data,x0,xn,y0,yn,z0,zn,frame_num,dotx,doty,dotn,
			 ampl,ampmin,ampmax,dsize,varamp)
     float ***data;
     int x0,xn,y0,yn,z0,zn,frame_num;
     float *dotx,*doty;               // [dotn]
     int dotn;
     float ampl,ampmin,ampmax,dsize;  // Set either 'ampl' or 'varamp'
     float *varamp;                   // [dotn] Variable dot ampl, or NULL
{
  int i,j,k;
  int nio2,ni,xi,yi,xfin,yfin;
  float x1,y1,d,damp;

  // Size of square region w/i which to draw each Gaussian dot
  nio2 = (int)(0.5 + dsize*3.0);
  if (nio2 < 1)
    nio2 = 1;
  ni = 2*nio2 + 1;

  for(i=0;i<dotn;i++){
    x1 = dotx[i];  // Coordinates for this dot
    y1 = doty[i];

    if (varamp == NULL)
      damp = ampl;
    else
      damp = varamp[i];

    xi = my_rint(x1) - nio2;
    yi = my_rint(y1) - nio2;
    xfin = xi+ni;
    yfin = yi+ni;
    for(j=xi;j<xfin;j++){
      for(k=yi;k<yfin;k++){  // Over the square region around the dot
	if ((j>=0) && (j<xn) && (k>=0) && (k<yn)){
	  d = data[x0+j][y0+k][frame_num];

	  //
	  //  Could use look-up table to speed this up, as dot size gets
	  //  larger, this gets quite slow.
	  //
	  d += damp * func_2d_gaussian((float)j,(float)k,x1,y1,dsize,dsize,0);
	  if (d > ampmax)
	    d = ampmax;
	  else if (d < ampmin)
	    d = ampmin;

	  data[x0+j][y0+k][frame_num] = d;
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MAKE_FRAME_ADD_DOTS_0                           */
/*                                                                           */
/*  Draw dots into the frame.  If 'clearflag' is 1, set the bgamp value      */
/*  before drawing dots.                                                     */
/*                                                                           */
/*****************************************************************************/
void make_frame_add_dots_0(data,xn,yn,dotx,doty,dotn,ampl,ampmin,ampmax,dsize,
			   varamp)
     float **data;
     int xn,yn;
     float *dotx,*doty;
     int dotn;
     float ampl,ampmin,ampmax,dsize;
     float *varamp;                   // [dotn] Variable dot ampl, or NULL
{
  int i,j,k;
  int nio2,ni,xi,yi,xfin,yfin;
  float x1,y1,d,damp;

  // Size of square region w/i which to draw each Gaussian dot
  nio2 = (int)(0.5 + dsize*3.0);
  if (nio2 < 1)
    nio2 = 1;
  ni = 2*nio2 + 1;

  //printf("ampmin = %f   ampmax = %f\n",ampmin,ampmax);
  //printf("data[10][10] = %f\n",data[10][10]);
  //printf("varamp[0] = %f  [1] = %f\n",varamp[0],varamp[1]);
  // ampmin = 0.0;

  for(i=0;i<dotn;i++){
    x1 = dotx[i];  // Coordinates for this dot
    y1 = doty[i];

    if (varamp == NULL)
      damp = ampl;
    else
      damp = varamp[i];

    xi = my_rint(x1) - nio2;
    yi = my_rint(y1) - nio2;
    xfin = xi+ni;
    yfin = yi+ni;
    for(j=xi;j<xfin;j++){
      for(k=yi;k<yfin;k++){  // Over the square region around the dot
	if ((j>=0) && (j<xn) && (k>=0) && (k<yn)){
	  d = data[j][k];

	  d += damp * func_2d_gaussian((float)j,(float)k,x1,y1,dsize,dsize,0);
	  if (d > ampmax)
	    d = ampmax;
	  else if (d < ampmin)
	    d = ampmin;

	  data[j][k] = d;
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         MAKE_FRAME_ADD_DOTS_CLIPR                         */
/*                                                                           */
/*  Draw dots into the frame.  If 'clearflag' is 1, set the bgamp value      */
/*  before drawing dots.                                                     */
/*                                                                           */
/*****************************************************************************/
void make_frame_add_dots_clipr(data,xn,yn,dotx,doty,dotn,ampl,ampmin,ampmax,
			       dsize,varamp,crx1,crx2,cry1,cry2)
     float **data;        // [xn][yn] draw into this 2D frame
     int xn,yn;           // 
     float *dotx,*doty;   // [dotn] dot coordinates (pix)
     int dotn;            // Number of dots
     float ampl;          // Amplitude of dot
     float ampmin,ampmax; // Truncate values outside this range
     float dsize;         // SD of Gaussian for dot (pix)
     float *varamp;       // [dotn] Variable dot ampl, or NULL
     float crx1,crx2;     // left and right edges of clip region
     float cry1,cry2;     // top and bottom edges of clip region
{
  int i,j,k;
  int nio2,ni,xi,yi,xfin,yfin,flag,inflag,ix1,ix2,iy1,iy2;
  float x1,y1,d,damp;

  // Size of square region w/i which to draw each Gaussian dot
  nio2 = (int)(0.5 + dsize*3.0);
  if (nio2 < 1)
    nio2 = 1;
  ni = 2*nio2 + 1;

  //printf("###crx1,2  y1,y2 = %f %f %f %f\n",crx1,crx2,cry1,cry2);

  ix1 = my_rint(crx1);  // Define integer pixel region corresponding to the
  ix2 = my_rint(crx2);  //   interior of the patch.
  iy1 = my_rint(cry1);
  iy2 = my_rint(cry2);

  //printf("ampmin = %f   ampmax = %f\n",ampmin,ampmax);
  //printf("data[10][10] = %f\n",data[10][10]);
  //printf("varamp[0] = %f  [1] = %f\n",varamp[0],varamp[1]);
  // ampmin = 0.0;

  for(i=0;i<dotn;i++){
    x1 = dotx[i];  // Coordinates for this dot
    y1 = doty[i];

    if (varamp == NULL)
      damp = ampl;
    else
      damp = varamp[i];

    xi = my_rint(x1) - nio2;
    yi = my_rint(y1) - nio2;
    xfin = xi+ni - 1;
    yfin = yi+ni - 1;

    if ((x1 >= crx1) && (x1 <= crx2) && (y1 >= cry1) && (y1 <= cry2)){

      inflag = 1; // Dot is inside clip region
      //
      //  Thus, we simply limit the rectangle in which the Gaussian is drawn,
      //  so that none of this dot is drawn outside the clip region.
      //
      if (xi   < ix1)  xi   = ix1;
      if (xfin > ix2)  xfin = ix2;
      if (yi   < iy1)  yi   = iy1;
      if (yfin > iy2)  yfin = iy2;
    }else{
      inflag = 0;
    }

    for(j=xi;j<=xfin;j++){
      for(k=yi;k<=yfin;k++){  // Over the square region around the dot

	flag = 1;
	if (inflag == 0){
	  //
	  //  If the dot is outside the rectangle, we need to block any pixel
	  //  that falls within the rectangle
	  //
	  if ((j >= ix1) && (j <= ix2) &&
	      (k >= iy1) && (k <= iy2))
	    flag = 0;

	  //flag = 0;  // WYETH TESTING - block dots outside clip region
	}

	if ((j>=0) && (j<xn) && (k>=0) && (k<yn) && (flag==1)){
	  d = data[j][k];

	  d += damp * func_2d_gaussian((float)j,(float)k,x1,y1,dsize,dsize,0);
	  if (d > ampmax)
	    d = ampmax;
	  else if (d < ampmin)
	    d = ampmin;

	  data[j][k] = d;
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MAKE_FRAME_GP_TRANS                          */
/*                                                                           */
/*  Fill frame "frame_num" with a translational Glass Pattern.  Dots are     */
/*  modeled as Gaussians.                                                    */
/*                                                                           */
/*    sscale - degrees (of vis. angle) / pixel                               */
/*    r - distance between dots (degr)                                       */
/*    theta - dot pair orientation                                           */
/*    amp1 - luminance of dot1                                               */
/*    amp2 - luminance of dot2                                               */
/*    ampbg - background luminance                                           */
/*    dsize - SD of Gaussian dot                                             */
/*    xf - expansion factor:  dots are picked in a field xn*xf by yn*xf,     */
/*         and only those falling in the corner xn by yn are used.           */
/*    dn - number of dot pairs per frame in expanded field.                  */
/*                                                                           */
/*  NOTES                                                                    */
/*  - The caller must handle the seed, do not reseed here.                   */
/*                                                                           */
/*****************************************************************************/
void make_frame_gp_trans(data,x0,xn,y0,yn,z0,zn,frame_num,sscale,r,theta,
			 amp1,amp2,ampbg,ampmin,ampmax,dsize,xf,dn,rseed)
     float ***data;
     int x0,xn,y0,yn,z0,zn,frame_num;
     float sscale,r,theta,amp1,amp2,ampbg,ampmin,ampmax,dsize,xf;
     int dn,*rseed;
{
  int i,j,k;
  int xi,yi,ni,nio2,xfin,yfin;
  float x,y,x1,y1,x2,y2,fxn,fyn;
  float dx,dy;

  // Dot displacments from center of pair
  dx = (r/2.0)/sscale * cos(theta/180.0 * M_PI);
  dy = (r/2.0)/sscale * sin(theta/180.0 * M_PI);

  // Size of field for dot generation
  fxn = xf * (float)xn;
  fyn = xf * (float)yn;
  
  // Set background value
  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++)
      data[x0+i][y0+j][frame_num] = ampbg;

  // Size of square region w/i which to draw each Gaussian dot
  nio2 = (int)(0.5 + dsize*3.0);
  if (nio2 < 1)
    nio2 = 1;
  ni = 2*nio2 + 1;

  //printf("DOTS DRAWN ON %d x %d pixel field\n",ni,ni);

  // For each dot pair
  for(i=0;i<dn;i++){
    // Get dot coordinates (x1,y1) (x2,y2)
    x = fxn * myrand_util_ran2(rseed);
    y = fyn * myrand_util_ran2(rseed);
    x1 = x + dx;
    y1 = y + dy;
    x2 = x - dx;
    y2 = y - dy;
    
    // Wrap dots in expansion field
    while (x1 <    0) x1 += fxn;    while (x1 >= fxn) x1 -= fxn;
    while (y1 <    0) y1 += fyn;    while (y1 >= fyn) y1 -= fyn;
    while (x2 <    0) x2 += fxn;    while (x2 >= fxn) x2 -= fxn;
    while (y2 <    0) y2 += fyn;    while (y2 >= fyn) y2 -= fyn;

    // Draw dot 1
    xi = my_rint(x1) - nio2;
    yi = my_rint(y1) - nio2;
    xfin = xi+ni;
    yfin = yi+ni;
    for(j=xi;j<xfin;j++){
      for(k=yi;k<yfin;k++){
	if ((j>=0) && (j<xn) && (k>=0) && (k<yn)){
	  data[x0+j][y0+k][frame_num] += amp1*
	    (float)func_2d_gaussian((float)j,(float)k,x1,y1,dsize,dsize,0);
	  if (data[x0+j][y0+k][frame_num] > ampmax)
	    data[x0+j][y0+k][frame_num] = ampmax;
	  else if (data[x0+j][y0+k][frame_num] < ampmin)
	    data[x0+j][y0+k][frame_num] = ampmin;
	}
      }
    }
    // Draw dot 2
    xi = my_rint(x2) - nio2;
    yi = my_rint(y2) - nio2;
    xfin = xi+ni;
    yfin = yi+ni;
    for(j=xi;j<xfin;j++){
      for(k=yi;k<yfin;k++){
	if ((j>=0) && (j<xn) && (k>=0) && (k<yn)){
	  data[x0+j][y0+k][frame_num] += amp2*
	    (float)func_2d_gaussian((float)j,(float)k,x2,y2,dsize,dsize,0);
	  if (data[x0+j][y0+k][frame_num] > ampmax)
	    data[x0+j][y0+k][frame_num] = ampmax;
	  if (data[x0+j][y0+k][frame_num] < ampmin)
	    data[x0+j][y0+k][frame_num] = ampmin;
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MAKE_FRAME_GP_TRANS_0                          */
/*                                                                           */
/*  Fill frame "frame_num" with a translational Glass Pattern.  Dots are     */
/*  modeled as Gaussians.                                                    */
/*                                                                           */
/*    sscale - degrees (of vis. angle) / pixel                               */
/*    r - distance between dots (degr)                                       */
/*    theta - dot pair orientation                                           */
/*    amp1 - luminance of dot1                                               */
/*    amp2 - luminance of dot2                                               */
/*    ampbg - background luminance                                           */
/*    dsize - SD of Gaussian dot                                             */
/*    xf - expansion factor:  dots are picked in a field xn*xf by yn*xf,     */
/*         and only those falling in the corner xn by yn are used.           */
/*    dn - number of dot pairs per frame in expanded field.                  */
/*                                                                           */
/*****************************************************************************/
void make_frame_gp_trans_0(data,xn,yn,sscale,r,theta,
			   amp1,amp2,ampbg,ampmin,ampmax,dsize,xf,dn,seed)
     float **data;
     int xn,yn;
     float sscale,r,theta,amp1,amp2,ampbg,ampmin,ampmax,dsize,xf;
     int dn,seed;
{
  int i,j,k;
  int xi,yi,ni,nio2,xfin,yfin;
  float x,y,x1,y1,x2,y2,fxn,fyn;
  float dx,dy;

  // Dot displacments from center of pair
  dx = (r/2.0)/sscale * cos(theta/180.0 * M_PI);
  dy = (r/2.0)/sscale * sin(theta/180.0 * M_PI);

  // Size of field for dot generation
  fxn = xf * (float)xn;
  fyn = xf * (float)yn;

  // Set background value
  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++)
      data[i][j] = ampbg;

  // Size of square region w/i which to draw each Gaussian dot
  nio2 = (int)(0.5 + dsize*3.0);
  if (nio2 < 1)
    nio2 = 1;
  ni = 2*nio2 + 1;

  //printf("DOTS DRAWN ON %d x %d pixel field\n",ni,ni);

  if (seed > 0)
    seed = -seed;

  // For each dot pair
  for(i=0;i<dn;i++){
    // Get dot coordinates (x1,y1) (x2,y2)
    x = fxn * myrand_util_ran2(&seed);
    y = fyn * myrand_util_ran2(&seed);
    x1 = x + dx;
    y1 = y + dy;
    x2 = x - dx;
    y2 = y - dy;
    
    // Wrap dots in expansion field
    while (x1 <    0) x1 += fxn;    while (x1 >= fxn) x1 -= fxn;
    while (y1 <    0) y1 += fyn;    while (y1 >= fyn) y1 -= fyn;
    while (x2 <    0) x2 += fxn;    while (x2 >= fxn) x2 -= fxn;
    while (y2 <    0) y2 += fyn;    while (y2 >= fyn) y2 -= fyn;

    // Draw dot 1
    xi = my_rint(x1) - nio2;
    yi = my_rint(y1) - nio2;
    xfin = xi+ni;
    yfin = yi+ni;
    for(j=xi;j<xfin;j++){
      for(k=yi;k<yfin;k++){
	if ((j>=0) && (j<xn) && (k>=0) && (k<yn)){
	  data[j][k] += amp1*
	    (float)func_2d_gaussian((float)j,(float)k,x1,y1,dsize,dsize,0);
	  if (data[j][k] > ampmax)
	    data[j][k] = ampmax;
	  else if (data[j][k] < ampmin)
	    data[j][k] = ampmin;
	}
      }
    }
    // Draw dot 2
    xi = my_rint(x2) - nio2;
    yi = my_rint(y2) - nio2;
    xfin = xi+ni;
    yfin = yi+ni;
    for(j=xi;j<xfin;j++){
      for(k=yi;k<yfin;k++){
	if ((j>=0) && (j<xn) && (k>=0) && (k<yn)){
	  data[j][k] += amp2*
	    (float)func_2d_gaussian((float)j,(float)k,x2,y2,dsize,dsize,0);
	  if (data[j][k] > ampmax)
	    data[j][k] = ampmax;
	  if (data[j][k] < ampmin)
	    data[j][k] = ampmin;
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MAKE_FRAME_BINARY_NOISE_0                       */
/*                                                                           */
/*  Fill data array with binary noise, either 'v0' or 'v1' using the given   */
/*  'pix_size'.                                                              */
/*                                                                           */
/*****************************************************************************/
void make_frame_binary_noise_0(float **data,
			       int xn,
			       int yn,
			       int pix_size,
			       float v0,
			       float v1,
			       int seed)  // Set > 0 to reseed
{
  int i,j,k,l;
  float rvalue;

  if (seed > 0)
    genrand_init((unsigned long)seed);

  for(i=0;i<xn;i+=pix_size)
    for(j=0;j<yn;j+=pix_size){
      if (genrand_real1() > 0.5)
	rvalue = v1;
      else
	rvalue = v0;
      for(k=0;k<pix_size;k++)
	for(l=0;l<pix_size;l++)
	  if (((i+k) < xn) && ((j+l) < yn))
	    data[i+k][j+l] = rvalue;
    }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MAKE_FRAME_BINARY_NOISE                       */
/*                                                                           */
/*  Fill frame "frame_num" with binary noise ("value" or -"value")           */
/*  of the given pixel size.                                                 */
/*                                                                           */
/*  *** MUST SET *pseed NEGATIVE ON FIRST CALL ***                           */
/*                                                                           */
/*****************************************************************************/
void make_frame_binary_noise(data,x0,xn,y0,yn,frame_num,pix_size,value,pseed)
     float ***data;
     int x0,xn,y0,yn,frame_num,pix_size;
     float value;
     int *pseed;
{
  int i,j,k,l;
  float rvalue;

  // WYETH - could eventually replace with
  //    'make_frame_binary_noise_0' from above to avoid 'ran2'

  for(i=0;i<xn;i+=pix_size)
    for(j=0;j<yn;j+=pix_size){
      if (myrand_util_ran2(pseed) > 0.5)
	rvalue = value;
      else
	rvalue = -value;
      for(k=0;k<pix_size;k++)
	for(l=0;l<pix_size;l++)
	  if (((i+k) < xn) && ((j+l) < yn))
	    data[x0+i+k][y0+j+l][frame_num] = rvalue;
    }
}
/**************************************-**************************************/
/*                                                                           */
/*                           STIMU_FRAME_DISK_ARRAY                          */
/*                                                                           */
/*  Draw 'n' disks into the frame at the specified coordinates.              */
/*                                                                           */
/*****************************************************************************/
void stimu_frame_disk_array(data,xn,yn,x,y,r,n,amp)
     float **data;  // [xn][yn]
     int xn,yn;     // size of image frame
     float *x,*y;   // [n] x and y coords (pix)
     float r;       // radii of disks
     int n;         // number of disks to draw
     float amp;     // color to draw
{
  int i,j;
  int di,x0,y0,x1,y1;
  float xx,yy,d2,r2;

  r2 = r*r;  // Radius squared

  for(di=0;di<n;di++){  // For each disk
    //
    //  Determine the start and width of box for this disk
    //
    xx = x[di];
    yy = y[di];

    x0 = (int)(xx - r);
    y0 = (int)(yy - r);
    x1 = my_rint(xx + r + 0.5);
    y1 = my_rint(yy + r + 0.5);

    if (x0 <  0)   x0 = 0;
    if (y0 <  0)   y0 = 0;
    if (x1 >= xn)  x1 = xn-1;
    if (y1 >= yn)  y1 = yn-1;

    for(i=x0;i<=x1;i++){
      for(j=y0;j<=y1;j++){
	d2 = ((float)i-xx)*((float)i-xx) + ((float)j-yy)*((float)j-yy);
	if (d2 <= r2)
	  data[i][j] = amp;
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             STIM_UTIL_LIMIT_3D                            */
/*                                                                           */
/*  Impose maximum and minimum values on a 3D stimulus.                      */
/*                                                                           */
/*****************************************************************************/
void stim_util_limit_3d(data,x0,xn,y0,yn,z0,zn,min,max)
     float ***data;
     int x0,xn,y0,yn,z0,zn;
     float min,max;
{
  int i,j,k;

  for(i=x0;i<(x0+xn);i++)
    for(j=y0;j<(y0+yn);j++)
      for(k=z0;k<(z0+zn);k++){
	if (data[i][j][k] < min)
	  data[i][j][k] = min;
	else if (data[i][j][k] > max)
	  data[i][j][k] = max;
      }
}
/**************************************-**************************************/
/*                                                                           */
/*                             STIM_UTIL_LIMIT_2D                            */
/*                                                                           */
/*  Impose maximum and minimum values on a 2D stimulus.                      */
/*                                                                           */
/*****************************************************************************/
void stim_util_limit_2d(data,xn,yn,min,max)
     float **data;
     int xn,yn;
     float min,max;
{
  int i,j;

  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++){
      if (data[i][j] < min)
	data[i][j] = min;
      else if (data[i][j] > max)
	data[i][j] = max;
    }
}
/**************************************-**************************************/
/*                                                                           */
/*                                STIMX_PIXEL                                */
/*                                                                           */
/*****************************************************************************/
void stimx_pixel(ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,rframe)
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // Time units per stim frame
     int ti;                       // Time index (tscale units)
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframe;              // Return 2D frame into this address
{
  int i,j;
  int pw,ph,pcx,pcy;
  float t,**tframe;

  static int px,py,px2,py2;
  static float bgval,ampl,t0,tn;

  if (task == -1){  // Clean up
    ;
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //

    t0        = param_getf_exit(ppl,"st0");
    tn        = param_getf_exit(ppl,"stn");
    pcx       = param_geti_dflt(ppl,"pix_cx",xn/2);  // (pix)
    pcy       = param_geti_dflt(ppl,"pix_cy",yn/2);  // (pix)
    px        = param_geti_exit(ppl,"pix_x");        // (pix)
    py        = param_geti_exit(ppl,"pix_y");        // (pix)
    pw        = param_geti_dflt(ppl,"pix_w",1);      // (pix)
    ph        = param_geti_dflt(ppl,"pix_h",1);      // (pix)
    ampl      = param_getf_exit(ppl,"ampl");
    bgval     = param_getf_exit(ppl,"bgval");

    px += pcx;
    py += pcy;

    px2 = px + pw - 1;
    py2 = py + ph - 1;

    if (px  < 0)    exit_error("STIMX_PIXEL","pix_x is less than zero");
    if (py  < 0)    exit_error("STIMX_PIXEL","pix_y is less than zero");
    if (px2 >= xn)  exit_error("STIMX_PIXEL","pix_x + pix_w is too large");
    if (py2 >= yn)  exit_error("STIMX_PIXEL","pix_y + pix_h is too large");

  }else if (task == 1){
    //
    //  Return stimulus frame
    //

    // Make new flat frame
    tframe = get_const_2d_farray(xn,yn,bgval);

    t = ti*tscale; // Time (s)
    if ((t >= t0) && (t < t0+tn)){
      for(i=px;i<=px2;i++)
	for(j=py;j<=py2;j++)
	  tframe[i][j] = ampl;
    }

    *rframe = tframe;
  }else{
    exit_error("STIMX_PIXEL","unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             STIMU_STD_MAKE_DRAW                           */
/*                                                                           */
/*  Standardized routine to determine whether to make and return both        */
/*  left and right stimulus frames, and whether to draw (non-blank)          */
/*  stimuli into both frames.                                                */
/*                                                                           */
/*  The 'make...' flags indicate whether a stimulus array is returned.       */
/*  The 'draw...' flags indicate whether the stimulus array is non-blank.    */
/*                                                                           */
/*****************************************************************************/
void stimu_std_make_draw(ppl,rmake_l,rmake_r,rdraw_l,rdraw_r)
     struct param_pair_list *ppl;  // Stimulus parameter list
     int *rmake_l;
     int *rmake_r;
     int *rdraw_l;
     int *rdraw_r;
{
  char *stimform;
  char *monoflag;
  int make_left,make_right,draw_left,draw_right;

  stimform = param_getc_dflt(ppl,"stim_form","3d");
  monoflag = param_getc_dflt(ppl,"stim_monocular","both"); // "left","right"

  if (strcmp(stimform,"3d_b")==0){
    make_left  = 1;
    make_right = 1;
  }else{
    make_left  = 1;
    make_right = 0;
  }

  draw_left  = make_left;
  draw_right = make_right;

  if (strcmp(monoflag,"left")==0)
    draw_right = 0;
  else if (strcmp(monoflag,"right")==0)
    draw_left = 0;

  myfree(stimform);
  myfree(monoflag);

  *rmake_l = make_left;
  *rmake_r = make_right;
  *rdraw_l = draw_left;
  *rdraw_r = draw_right;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 STIMX_BAR                                 */
/*                                                                           */
/*****************************************************************************/
void stimx_bar(ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,rframl,rframr)
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // Time units per stim frame
     int ti;                       // Time index (tscale units)
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframl;              // Return 2D frame into this address
     float ***rframr;              // Return 2D frame into this address
{
  int i,j,k;
  int ctrmot;
  float **tframl,**tframr,t,r,x,y,maxlum,bf_dw,bf_dh,x0,y0;

  static float costh,sinth,sdaa,bgval,speed_l,speed_r;
  static float t0  ,tn  ,cx ,cy, len ,wid,theta,ori,ampl;
  static float t0_2,tn_2,cx2,cy2,len2,wid2,speed2,theta2,ampl2;
  static int draw_left,draw_right,nbar,bf_nw,bf_nh;
  static float x0_l,x0_r,y0_l,y0_r,disp_x,disp_y,*bx0,*by0,*bx0r,*by0r;
  static char *stimform = NULL;
  static char *monoflag;  // *** WYETH - Does not need to be 'static' here

  if (task == -1){  // Clean up
    ;  // Nothing to do
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //
    t0       = param_getf_exit(ppl,"st0");
    tn       = param_getf_exit(ppl,"stn");
    cx       = param_getf_exit(ppl,"cx");
    cy       = param_getf_exit(ppl,"cy");
    theta    = param_getf_exit(ppl,"direction");
    ctrmot   = param_geti_dflt(ppl,"center_motion",0);  // 1-center the motion
    len      = param_getf_exit(ppl,"length");
    wid      = param_getf_exit(ppl,"width");
    ampl     = param_getf_exit(ppl,"ampl");
    bgval    = param_getf_exit(ppl,"bgval");
    maxlum   = param_getf_exit(ppl,"maxlum");
    sdaa     = param_getf_exit(ppl,"sdaa");
    stimform = param_getc_dflt(ppl,"stim_form","3d");
    monoflag = param_getc_dflt(ppl,"stim_monocular","both"); // "left","right"

    ori = param_getf_dflt(ppl,"bar_ori",theta);

    costh = cos(theta*M_PI/180.0);
    sinth = sin(theta*M_PI/180.0);

    // Default is to draw only the left frame; may change below
    draw_left  = 1;
    draw_right = 0;
    nbar = 1;        // Default is 1 bar, may change below

    if (strcmp(stimform,"3d_b")==0){
      //
      //  Which eye images need to be created?
      //
      if (strcmp(monoflag,"right")==0){
	draw_right = 1;
	draw_left = 0;
      }else if (strcmp(monoflag,"both")==0){
	draw_right = 1;
      } // Otherwise, make just the left by default

      //
      //  To create 3D motion, set velocities independently in two eyes
      //
      disp_x   = param_getf_dflt(ppl,"disparity_x",0.0);  // (deg)
      disp_y   = param_getf_dflt(ppl,"disparity_y",0.0);  // (deg)
      speed_l  = param_getf_exit(ppl,"speed_l");          // (deg/s)
      speed_r  = param_getf_exit(ppl,"speed_r");          // (deg/s)

      //
      //  Compute speeds and start points for Left and Right images
      //
      x0_l = cx          - 0.5 * tn * speed_l * costh;
      x0_r = cx + disp_x - 0.5 * tn * speed_r * costh;
      y0_l = cy          - 0.5 * tn * speed_l * sinth;
      y0_r = cy + disp_y - 0.5 * tn * speed_r * sinth;

      nbar = 1;
      bx0 = (float *)myalloc(nbar*sizeof(float));
      by0 = (float *)myalloc(nbar*sizeof(float));
      bx0[0] = x0_l;
      by0[0] = y0_l;
      bx0r = (float *)myalloc(nbar*sizeof(float));
      by0r = (float *)myalloc(nbar*sizeof(float));
      bx0r[0] = x0_r;
      by0r[0] = y0_r;

    }else{
      speed_l  = param_getf_exit(ppl,"speed");

      if (ctrmot == 1){
	x0_l = cx - 0.5 * tn * speed_l * costh;
	y0_l = cy - 0.5 * tn * speed_l * sinth;
      }else{
	x0_l = cx;
	y0_l = cy;
      }

      len2 = -1.0;   // No second bar
      bf_nw = bf_nh = 1;

      if (paramfile_test_param(ppl,"field_wid_n")){
	//
	//  Bar field
	//
	bf_nw     = param_geti_exit(ppl,"field_wid_n");
	bf_nh     = param_geti_exit(ppl,"field_len_n");
	bf_dw     = param_getf_exit(ppl,"field_wid_off");  // (deg)
	bf_dh     = param_getf_exit(ppl,"field_len_off");  // (deg)


	nbar = bf_nw * bf_nh;
	stm_util_grid_coords_1d(NULL,bf_nw,bf_nh,ori,bf_dw,bf_dh,&bx0,&by0);
	add_const_farray(bx0,nbar,x0_l);
	add_const_farray(by0,nbar,y0_l);

      }else if (paramfile_test_param(ppl,"cx2")){
	//
	//  A second, independent bar
	//
	t0_2      = param_getf_exit(ppl,"st0_2");
	tn_2      = param_getf_exit(ppl,"stn_2");
	cx2       = param_getf_exit(ppl,"cx2");
	cy2       = param_getf_exit(ppl,"cy2");
	theta2    = param_getf_exit(ppl,"direction2");
	len2      = param_getf_exit(ppl,"length2");
	wid2      = param_getf_exit(ppl,"width2");
	speed2    = param_getf_exit(ppl,"speed2");
	ampl2     = param_getf_exit(ppl,"ampl2");

	if (ctrmot == 1){
	  cx2 -=  0.5 * tn * speed2 * cos(theta2*M_PI/180.0);
	  cy2 -=  0.5 * tn * speed2 * sin(theta2*M_PI/180.0);
	}

	nbar = 1;
	bx0 = (float *)myalloc(nbar*sizeof(float));
	by0 = (float *)myalloc(nbar*sizeof(float));
	bx0[0] = x0_l;
	by0[0] = y0_l;

	// *** WYETH - BIG HACK - 2nd Bar is treated by itself, because
	//     we do not currently store all the 2nd bar params in the bar
	//     lists, bx0 and by0 (just the x,y coords) ***
	//
      }else{
	nbar = 1;
	bx0 = (float *)myalloc(nbar*sizeof(float));
	by0 = (float *)myalloc(nbar*sizeof(float));
	bx0[0] = x0_l;
	by0[0] = y0_l;
      }
    }
    myfree(monoflag);
  }else if (task == 1){
    //
    //  Return stimulus frame
    //
    t = ti*tscale; // Time (s)

    // Create new flat frame
    tframl = get_zero_2d_farray(xn,yn);
    if (strcmp(stimform,"3d_b")==0)
      tframr = get_zero_2d_farray(xn,yn);

    //
    //  Add 1st bar into frame
    //
    if ((t >= t0) && (t < t0+tn)){

      if (draw_left == 1){
	for(i=0;i<nbar;i++){
	  x = bx0[i] + speed_l * (t - t0) * costh;
	  y = by0[i] + speed_l * (t - t0) * sinth;
	  make_frame_bar_add_0(tframl,xn,yn,x,y,ori,len,wid,sdaa,
			       sscale,(ampl-bgval),0);
	}
      }

      if (draw_right == 1){
	for(i=0;i<nbar;i++){
	  x = bx0r[i] + speed_r * (t - t0) * costh;
	  y = by0r[i] + speed_r * (t - t0) * sinth;
	  make_frame_bar_add_0(tframr,xn,yn,x,y,ori,len,wid,sdaa,
			       sscale,(ampl-bgval),0);
	}
      }
    }

    if (len2 >= 0.0){  // If there is a 2nd bar
      //
      //  Add 2nd bar into frame
      //
      if ((t >= t0_2) && (t < t0_2+tn_2)){
	x = cx2 + speed2*(t-t0) * cos(theta2*M_PI/180.0);
	y = cy2 + speed2*(t-t0) * sin(theta2*M_PI/180.0);
	make_frame_bar_add_0(tframl,xn,yn,x,y,theta2,
			     len2,wid2,sdaa,sscale,(ampl2-bgval),0);
      }
    }

    *rframl = tframl;
    if (strcmp(stimform,"3d_b")==0)
      *rframr = tframr;

  }else{
    exit_error("STIMX_BAR","unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            STIMX_GABOR_BRITTEN                            */
/*                                                                           */
/*  Replicate the Britten and Heuer (1999) stimulus.                         */
/*                                                                           */
/*****************************************************************************/
void stimx_gabor_britten(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,rframl)
     char *mylogf;                 // Log file
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // Time units per stim frame
     int ti;                       // Time index (tscale units)
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframl;              // Return 2D frame into this address
{
  int i,k;
  int si,kloc,gxi,gyi,seed,n_epoch,*scnt,npos,rs_off;
  float durfc,durmp,dursp,mdd,sf,sdo,sdp,etf,kmidf,**tframl,t,maxlum,phase;
  char ggstr[LONG_SLEN];
  static int nstim,n_patch,t0i,gxn,nf_fc,nf_rmp,nf_sep,nf_tot;
  static int **poslist = NULL;
  static float bgval,t0,tn,cx,cy,theta,mdx,mdy,gx0,gy0,gdx,contrast,gamp;
  static struct stmh_gpatch *gp1 = NULL;

  if (task == -1){  // Clean up
    ;  // Nothing to do
  }else if (task == 0){

    //  Clean up from last prep.
    if (poslist != NULL)
      free_2d_iarray(poslist,nstim);  // 'nstim' should still have old value

    //
    //  Prepare for stimulus creation
    //
    t0       = param_getf_exit(ppl,"st0");
    tn       = param_getf_exit(ppl,"stn");
    cx       = param_getf_exit(ppl,"cx");
    cy       = param_getf_exit(ppl,"cy");
    theta    = param_getf_exit(ppl,"direction");
    sdo      = param_getf_exit(ppl,"sd_orth");
    sdp      = param_getf_exit(ppl,"sd_par");
    bgval    = param_getf_exit(ppl,"bgval");
    contrast = param_getf_exit(ppl,"contrast");
    maxlum   = param_getf_exit(ppl,"maxlum");
    sf       = param_getf_exit(ppl,"sf");
    n_epoch  = param_geti_exit(ppl,"n_epoch");
    n_patch  = param_geti_exit(ppl,"n_patch");
    gxn      = param_geti_exit(ppl,"grid_xn");    // Number of points
    gdx      = param_getf_exit(ppl,"grid_dx");    // (deg) separation
    durfc    = param_getf_exit(ppl,"dur_full");   // (s) at full contrast
    durmp    = param_getf_exit(ppl,"dur_ramp");   // (s) to ramp from 0 to full
    dursp    = param_getf_exit(ppl,"dur_sep");    // (s) of extra blank bewteen
    etf      = param_getf_exit(ppl,"etf");        // (Hz) jump size wrt SF
    seed     = param_geti_exit(ppl,"seed");
    rs_off   = param_geti_exit(ppl,"ranseq_offset");
    phase    = param_getf_dflt(ppl,"phase",90.0);  // spatial phase
    gamp     = param_getf_exit(ppl,"g_amp");       // amplitude of Gabor

    //
    //  Compute the number of frames to draw each Gabor, and 'nstim'
    //
    //         - - - - - -           Example:  nf_fc  6
    //       -             -                   nf_rmp 2
    //     -                 -                 nf_sep 3
    //   -                     - - - -         nf_tot = 6 + 2*2 + 3 + 2 = 15
    //   1 2 3 4 5 6 7 8 9 0 1 2 3 4 5
    //
    stm_gabor_britt_timing(durfc,durmp,dursp,tscale,zn,t0,tsamp,n_epoch,
			   &nf_fc,   // Frames at full contrast
			   &nf_rmp,  // Frames to rame up (or down)
			   &nf_sep,  // Frames of additional separation
			   &nf_tot,  // Frames assigned to one full epoch
			   &nstim);  // Number of stimulus epochs in trial

    //printf("Frames:  %d + 2*%d + %d + 2 = %d\n",nf_fc,nf_rmp,nf_sep,nf_tot);
    sprintf(ggstr,"    Will show %d stimuli per trial\n",nstim);
    mylog(mylogf,ggstr);

    t0i = stm_get_ti0(zn,tscale,t0,tsamp);  // Time index at which stim starts


    mdd = etf / sf * tscale;
    mdx = mdd * cos(theta/180.0 * M_PI);
    mdy = mdd * sin(theta/180.0 * M_PI);
    //printf("  Motion step: %f (deg)  dx,dy = %f %f \n",mdd,mdx,mdy);

    //
    //  Set up a Gabor patch with basic parameters
    //
    //  Note, if phase is 90.0, then the Gabor amplitude will be less than 1.0.
    //  This 'gamp' is used to correct this amplitude, so that contrast 1.0
    //  would give maximum (0..1 or -1..1) amplitude range.
    //
    gp1 = stm_gpatch_get(0.0,0.0,theta,sf,0.0,sdp,sdo,phase,contrast*gamp);

    //
    //  Get grid spatial constants
    //
    gx0 = cx - (float)(gxn-1) * gdx / 2.0;   // left side of grid (deg)
    gy0 = cy - (float)(gxn-1) * gdx / 2.0;   // bottom of grid (deg)

    //printf("gx0 = %f  gy0 = %f\n",gx0,gy0);


    //
    //  Generate 'poslist[epoch][npatch]' to tell positions.  Let us say that
    //  positions, p, will be assigned as follows:
    //
    //     20 21 22 23 24        So that   x = p % 5
    //     15 16 17 18 19            and   y = (int)(p / 5);
    //     10 11 12 13 14
    //      5  6  7  8  9
    //      0  1  2  3  4
    //
    npos = gxn * gxn;
    stm_gabor_britt(mylogf,npos,n_patch,seed,rs_off,nstim,&poslist);


  }else if (task == 1){
    //
    //  Return stimulus frame
    //
    t = ti*tscale; // Time (s)

    // Compute the index into the random bar sequence to use for this frame
    k = (int)((ti-t0i)/nf_tot);

    // Compute the frame index w/i this stimulus patch
    kloc  = ti - (k*nf_tot + t0i);
    kmidf = (float)kloc - (float)(nf_tot - 1)/2.0;  // 'k' w.r.t. midpoint

    //printf("kmidf = %f\n",kmidf);

    // Make new flat frame
    tframl = get_zero_2d_farray(xn,yn);

    //if ((t >= t0) && (t < t0+tn)){
    if ((k >= 0) && (k < nstim) &&
	(kloc > 0) && (kloc < (nf_tot-1-nf_sep))){
	  
      //printf("kloc = %d\n",kloc);

      //
      //  Patch 1
      //
      gxi = poslist[k][0] % gxn;  // Grid x-index [0..gxn-1]
      gyi = poslist[k][0] / gxn;  // Grid y-index [0..gxn-1]

      gp1->x0 = gx0 + gxi * gdx;  // Grid point position in image frame (deg)
      gp1->y0 = gy0 + gyi * gdx;  // Grid point position in image frame (deg)

      gp1->x0 += kmidf*mdx;  // Adjust position for movement
      gp1->y0 += kmidf*mdy;

      if (kloc <= nf_rmp){
	gp1->amp = (float)kloc/(float)(nf_rmp+1) * contrast*gamp;
	//printf("%3d ramping UP  %f\n",kloc,gp1->amp);
      }else if (kloc > (nf_rmp + nf_fc)){
	gp1->amp = (float)(nf_rmp*2 + nf_fc + 1 - kloc)/(float)(nf_rmp+1)
	  * contrast*gamp;
	if (gp1->amp < 0.0)
	  gp1->amp = 0.0;
	//printf("%d ramping DOWN %f\n",kloc,gp1->amp);
      }else
	gp1->amp = contrast*gamp;

      //printf("xi,yi %d %d   x0,y0 = %f %f\n",gxi,gyi,gp1->x0,gp1->y0);
      make_frame_add_gpatch_gabor(tframl,xn,yn,gp1,sscale,0.0);  // ph = 0.0

      if (n_patch == 2){
	gxi = poslist[k][1] % gxn;
	gyi = poslist[k][1] / gxn;

	gp1->x0 = gx0 + gxi * gdx;
	gp1->y0 = gy0 + gyi * gdx;

	gp1->x0 += kmidf*mdx;
	gp1->y0 += kmidf*mdy;

	//printf("xi,yi %d %d   x0,y0 = %f %f\n",gxi,gyi,gp1->x0,gp1->y0);
	make_frame_add_gpatch_gabor(tframl,xn,yn,gp1,sscale,0.0);  // ph = 0.0
      }
    }

    *rframl = tframl;
  }else{
    exit_error("STIMX_GABOR_BRITTEN","unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               STIMX_3BARLINK                              */
/*                                                                           */
/*****************************************************************************/
void stimx_3barlink(ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,rframl)
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // Time units per stim frame
     int ti;                       // Time index (tscale units)
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframl;              // Return 2D frame into this address
{
  float **tframl,t,r,x,y,maxlum,th1,th2,lenx;
  static float sdaa,bgval,xlen;
  static float t0,tn,cx,cy,len,wid,theta,thetafl,ampl;
  static float x0_l,x0_r,y0_l,y0_r,disp_x,disp_y;

  if (task == -1){  // Clean up
    ;  // Nothing to do
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //
    t0       = param_getf_exit(ppl,"st0");
    tn       = param_getf_exit(ppl,"stn");
    cx       = param_getf_exit(ppl,"cx");
    cy       = param_getf_exit(ppl,"cy");
    theta    = param_getf_exit(ppl,"direction");
    thetafl  = param_getf_exit(ppl,"theta_flank");
    len      = param_getf_exit(ppl,"length");
    lenx     = param_getf_dflt(ppl,"length_x",0.0);
    wid      = param_getf_exit(ppl,"width");
    ampl     = param_getf_exit(ppl,"ampl");
    bgval    = param_getf_exit(ppl,"bgval");
    maxlum   = param_getf_exit(ppl,"maxlum");
    sdaa     = param_getf_exit(ppl,"sdaa");

    xlen = len + lenx*sscale;  // extra pixel length to connect to main bar

  }else if (task == 1){
    //
    //  Return stimulus frame
    //
    t = ti*tscale; // Time (s)

    // Make new flat frame
    tframl = get_zero_2d_farray(xn,yn);

    //
    //  Add 1st bar into frame
    //
    if ((t >= t0) && (t < t0+tn)){

      //
      //  Central bar
      //
      x = cx;
      y = cy;
      make_frame_bar_add_0(tframl,xn,yn,x,y,theta,len,wid,sdaa,
			   sscale,(ampl-bgval),0);

      //
      //  2nd bar
      //
      th1 = (90.0 + theta) * M_PI/180.0;
      th2 = thetafl        * M_PI/180.0;

      x = cx + len/2.0 * (cos(th1) + cos(th1 + th2));
      y = cy + len/2.0 * (sin(th1) + sin(th1 + th2));

      make_frame_bar_add_0(tframl,xn,yn,x,y,theta+thetafl,
			   xlen,wid,sdaa,sscale,(ampl-bgval),1);

      //
      //  3rd bar
      //
      th1 = (270.0 + theta) * M_PI/180.0;

      x = cx + len/2.0 * (cos(th1) + cos(th1 - th2));
      y = cy + len/2.0 * (sin(th1) + sin(th1 - th2));

      make_frame_bar_add_0(tframl,xn,yn,x,y,theta-thetafl,
			   xlen,wid,sdaa,sscale,(ampl-bgval),1);
    }

    *rframl = tframl;

  }else{
    exit_error("STIMX_3BARLINK","unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               STIMX_BARNOISE                              */
/*                                                                           */
/*****************************************************************************/
void stimx_barnoise(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,
		    rframl,rframr)
     char *mylogf;
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // Time units per stim frame
     int ti;                       // Time index (tscale units)
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframl;              // Return 2D frame, left eye
     float ***rframr;              // Return 2D frame, right eye
{
  int i,j,k;
  int seed,seed_r,distrib,nn,mask_seed,mask_distr,ii,dx_bi,mask_dir_flag;
  float v1,v2,v3,t0,tn,bw,bl,bgval,cx,cy,**mt,**mask,v,dx_pix,dx;
  float maxlum,contrast,mask_speed,mask_con,mask_diroff;
  char *monoflag;
  // ______ Static ______
  static float **bartemp = NULL;
  static float **seq = NULL;
  static float *mask_seq = NULL;
  static float **seq_r = NULL;
  static int mask_flag = 0;
  static float bgamp,mask_speed_pix,mask_dir;
  static float dir,bwpix,blpix,bord,cxpix,cypix;
  static int nbar,nseq,t0i,dwell,make_left,make_right,mask_seqn;
  static char *stimform = NULL;

  if (task == -1){  // Clean up
    if (bartemp != NULL)
      free_2d_farray(bartemp,xn);
    if (seq != NULL)
      free_2d_farray(seq,nseq);
    if (seq_r != NULL)
      free_2d_farray(seq_r,nseq);
    if (stimform != NULL)
      myfree(stimform);
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //

    nbar     = param_geti_exit(ppl,"nbar");
    if (paramfile_test_param(ppl,"bar_width_pix") == 1)
      bwpix    = (float)param_geti_exit(ppl,"bar_width_pix");
    else{
      bw       = param_getf_exit(ppl,"bar_width");  // (deg)
      bwpix    = bw / sscale;
    }
    bl        = param_getf_exit(ppl,"bar_length");
    blpix     = bl / sscale;

    cx        = param_getf_dflt(ppl,"cx",0.0);
    cy        = param_getf_dflt(ppl,"cy",0.0);
    t0        = param_getf_exit(ppl,"st0");
    tn        = param_getf_exit(ppl,"stn");
    dir       = param_getf_exit(ppl,"direction");
    dwell     = param_geti_exit(ppl,"dwell");  // Raw tscale units
    distrib   = param_geti_exit(ppl,"distrib");
    seed      = param_geti_exit(ppl,"seed");
    seed_r    = param_geti_dflt(ppl,"seed_r",seed);
    contrast  = param_getf_exit(ppl,"contrast");
    maxlum    = param_getf_exit(ppl,"maxlum");
    bgval     = param_getf_exit(ppl,"bgval");
    bord      = param_getf_dflt(ppl,"border",0.0);  // Number of pixels smooth
    stimform  = param_getc_dflt(ppl,"stim_form","3d");
    monoflag  = param_getc_dflt(ppl,"stim_monocular","both"); // "left","right"

    if (param_test(ppl,"mask_speed") == 1){
      mask_flag = 1;  // Need to make a drifting mask
      contrast *= 2.0;   // Image will be divided by two after mask is added

      if (param_test(ppl,"mask_direction") == 1){
	//
	//  If this parameter is set, then use this as the absolute direction
	//  of the moving mask, and make the noise pattern direction be
	//  relative to this.
	//
	mask_dir = param_getf_dflt(ppl,"mask_direction");
	mask_dir_flag = 1;
	dir += mask_dir;  // Make 'dir' be relative to 'mask_dir'
      }else
	mask_dir_flag = 0;
    }

    /*
      bgamp = 2.0*bgval - 1.0; // Convert 0..1 --> -1..1
      v1 = -1.0;  // Dark bar color
      v2 = 1.0;   // Light bar color
    */
    bgamp = bgval;
    v1 = bgval + (maxlum - bgval)*contrast;
    v2 = bgval - (maxlum - bgval)*contrast;
    v3 = bgval;  // This is only used if 'distrib' is ternary
    //printf("******************* v1,2,3  %f  %f  %f\n",v1,v2,v3);

    t0i = my_rint(t0 / tscale);  // Time index at which stim starts

    cxpix = cx / sscale;  // Center x (pix)
    cypix = cy / sscale;  // Center y (pix)

    make_left = make_right = 1;
    if (strcmp(monoflag,"left")==0)
      make_right = 0;
    else if (strcmp(monoflag,"right")==0)
      make_left = 0;
    myfree(monoflag);

    bartemp = stm_barnoise_get_template(mylogf,xn,yn,nbar,bwpix,blpix,dir,bord,
					cxpix,cypix,0.0);
    //write_2d_data("zz.template.2d",bartemp,0,0,xn,yn,4,2,1,1);
    //exit(0);

    if ((distrib >= 2) || (distrib <= 3)){
      // Returned array is [nseq][nbar]
      stm_get_ran_seq_2d(NULL,nbar,distrib,seed,1.0/tscale,tn,dwell,v1,v2,v3,
			 &seq,&nseq);
      if (strcmp(stimform,"3d_b")==0){
	stm_get_ran_seq_2d(NULL,nbar,distrib,seed_r,1.0/tscale,tn,dwell,v1,v2,
			   v3,&seq_r,&nseq);
      }
    }else{
      exit_error("STIMX_BARNOISE","bad 'distrib', must be 2 or 3");
    }

    //
    //  Add a drifting mask - proposed by Quaia and Cumming (Oct 2015)
    //
    if (mask_flag == 1){
      mask_speed   = param_getf_dflt(ppl,"mask_speed",0.0);  // deg/s
      if (mask_dir_flag == 0){
	mask_diroff  = param_getf_dflt(ppl,"mask_dir_offset",90.0);
	mask_dir = dir+mask_diroff;
      }
      mask_seed    = param_geti_exit(ppl,"mask_seed");
      mask_distr   = param_geti_dflt(ppl,"mask_distrib",distrib);
      mask_con     = param_getf_dflt(ppl,"mask_contrast",contrast);

      if (mask_con != contrast){
	v1 = bgval + (maxlum - bgval)*(mask_con*2);  // For divide/2 later
	v2 = bgval - (maxlum - bgval)*(mask_con*2);
      }

      mask_speed_pix = mask_speed / sscale;

      //
      //  Compute a sufficiently long random sequence of bar values
      //
      mask_seqn = 2.0 * xn / bwpix;  // Twice the number of bars to fill frame
      mask_seqn += tn * mask_speed / sscale  / bwpix;  // plus # of new bars

      stm_get_ran_seq(mylogf,mask_distr,0,mask_seed,1.0,(float)mask_seqn,1,
		      v1,v2,v3,&mask_seq,&nn);

      if (nn != mask_seqn)
	mylog_exit(mylogf,"STIMX_BARNOISE  'mask_seqn' error\n");
    }

  }else if (task == 1){
    //
    //  Return stimulus frame
    //

    if ((bartemp == NULL) || (seq == NULL))
      exit_error("STIMX_BARNOISE","Not prepared to generate stim");

    // Compute the index into the random bar sequence to use for this frame
    k = (int)((ti-t0i)/dwell);
    if ((k < 0) || (k >= nseq))
      exit_error("STIMX_BARNOISE","index exceeds random sequence length");

    if (make_left == 1)
      *rframl = make_frame_bars_tmplt(bartemp,xn,yn,nbar,seq[k],bgamp);
    else
      *rframl = get_const_2d_farray(xn,yn,bgamp);

    if (mask_flag == 1){

      // Compute current distance moved, first in pix, then in bars
      dx_pix = (float)(ti-t0i)*tscale  *  mask_speed_pix;  // s * pix/s
      dx_bi  = (int)(dx_pix / bwpix);  // Integer number of bars moved by
      dx = dx_pix - ((float)dx_bi * bwpix);
      //printf("ti %d  dx_pix %f  dx_bar %f  bi %d\n",ti,dx_pix,dx_bar,dx_bi);

      ii = (mask_seqn - nbar) - dx_bi;  // Reverse from end of bar value array
      if (ii < 0)
	mylog_exit(mylogf,"STIMX_BARNOISE  bar index < 0\n");

      mt = stm_barnoise_get_template(mylogf,xn,yn,nbar,bwpix,blpix,
				     mask_dir,bord,cxpix,cypix,dx);
                                     //dir+mask_diroff,bord,cxpix,cypix,dx);
      mask = make_frame_bars_tmplt(mt,xn,yn,nbar,&(mask_seq[ii]),bgamp);

      for(i=0;i<xn;i++){
	for(j=0;j<yn;j++){
	  v = ((*rframl)[i][j] + mask[i][j]) / 2.0;
	  (*rframl)[i][j] = v;
	}
      }
      free_2d_farray(mask,xn);
      free_2d_farray(mt,xn);
    }

    if (strcmp(stimform,"3d_b")==0){
      if (make_right == 1)
	*rframr = make_frame_bars_tmplt(bartemp,xn,yn,nbar,seq_r[k],bgamp);
      else
	*rframr = get_const_2d_farray(xn,yn,bgamp);
    }
  }else{
    exit_error("STIMX_BARNOISE","unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               STIMX_BARPATMO                              */
/*                                                                           */
/*****************************************************************************/
void stimx_barpatmo(ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,rframe)
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // Time units per stim frame
     int ti;                       // Time index (tscale units)
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframe;              // Return 2D frame into this address
{
  int i,k;
  int seed,distrib,npatframe,npixnew,si,pi;
  float dir,t0,tn,bw,bl,cx,cy;

  static int nbar,seqn,t0i,dwell,bwpix,blpix,speedpix,x0,y0,altcon,altsign;
  static int idir;
  static int *seqi = NULL;
  static float bgval,val1,val2;

  if (task == -1){  // Clean up
    if (seqi != NULL)
      myfree(seqi);
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //

    nbar      = param_geti_exit(ppl,"nbar");
    if (param_test(ppl,"bar_width_pix") == 1)
      bwpix   = param_geti_exit(ppl,"bar_width_pix");
    else{
      bw      = param_getf_exit(ppl,"bar_width");  // (deg)
      bwpix    = bw / sscale;
    }
    bl        = param_getf_exit(ppl,"bar_length");
    cx        = param_getf_dflt(ppl,"cx",0.0);
    cy        = param_getf_dflt(ppl,"cy",0.0);
    t0        = param_getf_exit(ppl,"st0");
    tn        = param_getf_exit(ppl,"stn");
    dir       = param_getf_exit(ppl,"direction");
    dwell     = param_geti_exit(ppl,"dwell");  // Raw tscale units
    distrib   = param_geti_exit(ppl,"distrib");
    seed      = param_geti_exit(ppl,"seed");
    speedpix  = param_geti_exit(ppl,"speed_pix");  // pix/pat frame
    val1      = param_getf_dflt(ppl,"bar_val_1",1.0);
    val2      = param_getf_dflt(ppl,"bar_val_2",0.0);
    bgval     = param_getf_exit(ppl,"bgval");
    altcon    = param_geti_dflt(ppl,"alternate",0);

    idir = my_rint(dir);
    if ((idir != 0) && (idir != 180)){
      exit_error("STIMX_BARPATMO","direction must be 0 or 180 deg");
    }

    blpix    = bl / sscale;

    t0i       = my_rint(t0 / tscale);  // Time index at which stim starts

    //  (x0,y0) is lower left pixel of bar pattern
    x0 = my_rint(xn/2 + cx / sscale - (float)(nbar * bwpix)/2.0);
    y0 = my_rint(yn/2 + cy / sscale - (float)blpix/2.0);

    if (seed > 0)
      seed = -seed;

    //
    //  Number of pattern frames, with unique pattern
    //
    npatframe = (int)(tn / (dwell*tsamp*tscale));

    /*
 printf("tn = %f\n",tn);
 printf("tsamp = %d\n",tsamp);
 printf("tscale = %f\n",tscale);
 printf("npatframe = %d\n",npatframe);
 printf("nbar = %d\n",nbar);
 printf("speedpix = %d\n",speedpix);
 printf("bwpix = %d\n",bwpix);
    */

    //
    //  Compute 'seqn', the number of random values needed to set the
    //  lum/color of all bars that will ultimately appear as the pattern moves.
    //
    seqn = nbar; // Number of random values for bar lum/color
    if (speedpix > 0){
      npixnew = speedpix * npatframe;   // Total new exposed pixels
      seqn += npixnew / bwpix + 1;      // Add implied number of bars
    }

    // printf("seqn = %d\n",seqn);
    // printf("dwell = %d\n",dwell);

    if (distrib == 0){
      seqi = (int *)myalloc(seqn*sizeof(int));
      for(i=0;i<seqn;i++){ // Convert 0/1 --> -1/1
	if (i%2 == 0)
	  seqi[i] = 1;
	else
	  seqi[i] = -1;
      }
    }else if (distrib == 2){
      seqi = myrand_get_std_bin_seq(seqn,seed);
      for(i=0;i<seqn;i++) // Convert 0/1 --> -1/1
	if (seqi[i] < 1)
	  seqi[i] = -1;
    }else if (distrib == 3){
      seqi = myrand_get_std_tern_seq(seqn,seed);
    }else{
      exit_error("STIMX_BARPATMO","bad 'distrib', must be 2 or 3");
    }

    altsign = 1;  // Start with positive contrast pattern

  }else if (task == 1){
    //
    //  Return stimulus frame
    //

    if (seqi == NULL)
      exit_error("STIMX_BARPATMO","Not prepared to generate stim");

    //
    // Compute the index 'si' into the random sequence for this frame
    //
    k =  (int)((ti-t0i)/dwell);  // Number of bar movements
    si = k * speedpix / bwpix;   // Number of pixels passed / bwpix
    pi = k * speedpix % bwpix;   // pixel index

    //if ((si < 0) || (si >= seqn)){
    if ((k < 0) || (k >= npatframe)){  // This frame not w/i stimulus
      //
      //  Fill frame w/ background
      //
      *rframe = get_const_2d_farray(xn,yn,bgval);
    }else{
      if (altcon == 1){
	if (k%2 == 1)  // Do not alternate during dwell
	  altsign = -1;
	else
	  altsign = 1;
      }
    }

    *rframe = make_frame_barpatmo(xn,yn,x0,y0,nbar,bwpix,blpix,seqi,seqn,
				  si,pi,altsign,val1,val2,bgval,idir);
  }else{
    exit_error("STIMX_BARPATMO","unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                STIMX_BARRAY                               */
/*                                                                           */
/*****************************************************************************/
void stimx_barray(ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,rframe)
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // Time units per stim frame
     int ti;                       // Time index (tscale units)
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframe;              // Return 2D frame into this address
{
  int i,j,k;
  int si,bi,seed;
  int dt_off,bar_off;
  float t0,tn,contrast,maxlum,**tframe,t,x,y;
  float bar_gap,pos_gap,size;

  // Static
  static float **bar_cx,**bar_cy;
  static float theta,barw,barl,barval,bgval,sdaa,cx,cy;
  static int n,nbar,npos,t0i,dtt,seqtype,dt;
  static int **seqi;

  if (task == -1){  // Clean up
    if (bar_cx != NULL)
      free_2d_farray(bar_cx,nbar);
    if (bar_cy != NULL)
      free_2d_farray(bar_cy,nbar);
    if (seqi != NULL)
      free_2d_iarray(seqi,n);
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //
    t0       = param_getf_exit(ppl,"st0");
    tn       = param_getf_exit(ppl,"stn");
    cx       = param_getf_exit(ppl,"cx");
    cy       = param_getf_exit(ppl,"cy");
    theta    = param_getf_exit(ppl,"direction");
    barw     = param_getf_exit(ppl,"bar_w");
    bar_gap  = param_getf_exit(ppl,"bar_gap");
    pos_gap  = param_getf_exit(ppl,"pos_gap");
    nbar     = param_geti_exit(ppl,"bar_n");
    npos     = param_geti_exit(ppl,"pos_n");
    dt       = param_geti_exit(ppl,"dt");
    dt_off   = param_geti_exit(ppl,"dt_off");
    seed     = param_geti_exit(ppl,"seed");
    contrast = param_getf_exit(ppl,"contrast");
    barval   = param_getf_exit(ppl,"bar_amp");
    bgval    = param_getf_exit(ppl,"bgval");
    maxlum   = param_getf_exit(ppl,"maxlum");
    seqtype  = param_geti_exit(ppl,"seq_type");
    size     = param_getf_exit(ppl,"size");
    sdaa     = param_getf_dflt(ppl,"sdaa",0.0);

    dtt = dt + dt_off;

    barl = (size - (float)(nbar-1)*bar_gap)/(float)nbar; // Bar length (degr)
    n = (int)(tn/tscale / (float)(dt+dt_off)); // Bar patterns in random sequ.

    // Get the random sequence [n][nbar]
    seqi = stm_barray_get_seq(NULL,n,nbar,npos,seqtype,seed);

    // Get center coords for each possible bar in the grid [nbar][npos]
    // in degrees, w.r.t. 0,0  (Thus x,y, offsets must be added)
    stm_barray_get_bar_centers(NULL,nbar,npos,theta,barw,barl,bar_gap,
			       pos_gap,&bar_cx,&bar_cy);

    t0i = my_rint(t0 / tscale);  // Time index at which stim starts

  }else if (task == 1){
    //
    //  Get stimulus frame
    //

    if (seqi == NULL)
      exit_error("STIMX_BARRAY","Not prepared to generate stim");

    // Compute the index into the random bar sequence to use for this frame
    si =  (int)((ti-t0i)/dtt);

    // Determine whether a bar is on the screen now, or not
    //printf("ti-t0i = %d   dtt = %d   si = %d\n",ti-t0i,dtt,si);
    if (((ti-t0i) - dtt*si) >= dt)
      bar_off = 1;
    else
      bar_off = 0;

    if ((si < 0) || (si >= n)){
      printf("  si = %d\n",si);
      exit_error("STIMX_BARRAY","index exceeds random sequence length");
    }

    // Make new flat frame
    tframe = get_2d_farray(xn,yn);
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	tframe[i][j] = bgval;

    // Draw each bar (or gray if no bar)

    if (bar_off == 0){
      for(j=0;j<nbar;j++){
      
	bi = seqi[si][j] - 1;
	//printf("  %d bar %d    %d\n",i,j,bi);
	if (bi >= 0){ // This bar is drawn

	  if (seqtype == 99){ // Draw all bars - for demo
	    for(k=0;k<npos;k++){
	      x = bar_cx[j][k] + cx;
	      y = bar_cy[j][k] + cy;
	      make_frame_bar_0(tframe,xn,yn,x,y,theta,barl,barw,
			       sdaa,sscale,barval,bgval,1);
	    }
	  }else{
	    x = bar_cx[j][bi] + cx;
	    y = bar_cy[j][bi] + cy;
	    //printf("i=%d bar=%d  pos=%d  x,y  %f %f\n",i,j,bi,x,y);
	    make_frame_bar_0(tframe,xn,yn,x,y,theta,barl,barw,
			     sdaa,sscale,barval,bgval,1);
	  }
	}
      }
    }

    *rframe = tframe;

  }else{
    exit_error("STIMX_BARRAY","unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               STIMX_BAR2RAN                               */
/*                                                                           */
/*****************************************************************************/
void stimx_bar2ran(ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,rframe)
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // Time units per stim frame
     int ti;                       // Time index (tscale units)
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframe;              // Return 2D frame into this address
{
  int i,j;
  int seed,ri,*seedlist,bwflag;
  int xr0,yr0,xr1,yr1,hpix,wpix;
  float t0,tn,cx,cy,theta,w,h,barw,barh,contrast,maxlum,wf,hf;
  float **tframe,t,x,y;
  float xv1,yv1,xv0,yv0;

  static int nn,x0,y0,bw,bh,t0i,dwell;
  static int *xpix0,*xpix1,*ypix0,*ypix1;
  static float bgval,val1,val2,diamf;
  static float *xf0,*xf1,*yf0,*yf1;  // Float coords for bumps
  static char *shape;


  if (task == -1){  // Clean up
    if (xpix0 != NULL)  myfree(xpix0);
    if (xpix1 != NULL)  myfree(xpix1);
    if (ypix0 != NULL)  myfree(ypix0);
    if (ypix1 != NULL)  myfree(ypix1);
    if (xf0 != NULL)  myfree(xf0);
    if (xf1 != NULL)  myfree(xf1);
    if (yf0 != NULL)  myfree(yf0);
    if (yf1 != NULL)  myfree(yf1);
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //

    t0       = param_getf_exit(ppl,"st0");
    tn       = param_getf_exit(ppl,"stn");
    cx       = param_getf_exit(ppl,"cx");
    cy       = param_getf_exit(ppl,"cy");
    theta    = param_getf_exit(ppl,"direction");
    w        = param_getf_exit(ppl,"ranw");
    h        = param_getf_exit(ppl,"ranh");
    barw     = param_getf_exit(ppl,"barw");
    barh     = param_getf_exit(ppl,"barh");
    seed     = param_geti_exit(ppl,"seed");
    contrast = param_getf_exit(ppl,"contrast");
    bgval    = param_getf_exit(ppl,"bgval");
    maxlum   = param_getf_exit(ppl,"maxlum");
    dwell    = param_geti_exit(ppl,"dwell");
    //voidval  = param_getf_dflt(ppl,"stim_blank_val",0.0);
    bwflag   = param_geti_dflt(ppl,"black_white",0);
    shape    = paramfile_get_char_param_default(ppl,"shape","rect");

    if (paramfile_test_param(ppl,"stim_void") == 1){
      exit_error("STIMX_BAR2RAN","Change stim_void to stim_blank_val");
    }
    
    if (theta != 0.0)
      exit_error("STIMX_BAR2RAN","Direction must be zero");

    bw = my_rint(barw/sscale);  // Bar width in pixels
    bh = my_rint(barh/sscale);  // Bar height in pixels
    wpix = my_rint(w/sscale);   // Range in pixels
    hpix = my_rint(h/sscale);   // Range in pixels
    
    diamf = barw/sscale;    // Bump size (pix)
    wf = w/sscale;              // Range (pix)
    hf = h/sscale;              // Range (pix)
    
    
    if (bwflag == 0){
      val1 = contrast/2.0;      // white bar
      val2 = -contrast/2.0;     // dark bar
    }else if (bwflag == -1){
      // bgval should be 1.0, so these are both subtracted
      val1 = -contrast/2.0;     // dark bar
      val2 = -contrast/2.0;     // dark bar
    }else if (bwflag == 1){
      // bgval should be 0.0, so these are both added
      val1 = contrast/2.0;     // white bar
      val2 = contrast/2.0;     // white bar
    }

    // WYETH - possible violation of centering convention: (xn-1)/2  ??
    x0 = xn/2 + (cx - w/2.0 - barw/2.0)/sscale; // Pix left edge of xrange
    y0 = yn/2 + (cy - h/2.0 - barh/2.0)/sscale; // Pix bottom of y-range

    //printf("bw,bh  %d %d   x0,y0  %d %d\n",bw,bh,x0,y0);

    // Compute number of new random frames to create
    //
    //  WYETH - this differs from stimgen_bar2ran way.  That way is probably
    //          better.
    //
    nn = stm_util_count_new_frames(zn,tscale,t0,tn,tsamp,dwell);

    seedlist = get_seeds(seed,100000,4);

    if (strcmp(shape,"rect")==0){
      // Get uniform random integer values from 0..wpix-1
      xpix0 = myrand_get_std_unif_int_seq(nn,seedlist[0],wpix);
      xpix1 = myrand_get_std_unif_int_seq(nn,seedlist[1],wpix);
      if (h > 0.0){ // There is also random placement in y dimension
	ypix0 = myrand_get_std_unif_int_seq(nn,seedlist[2],hpix);
	ypix1 = myrand_get_std_unif_int_seq(nn,seedlist[3],hpix);
      }else{
	ypix0 = get_const_iarray(nn,0);
	ypix1 = get_const_iarray(nn,0);
      }
    }else if (strcmp(shape,"sine_bump")==0){
      xf0 = myrand_get_std_unif_float_seq(nn,seedlist[0],wf);
      xf1 = myrand_get_std_unif_float_seq(nn,seedlist[1],wf);
      yf0 = myrand_get_std_unif_float_seq(nn,seedlist[2],hf);
      yf1 = myrand_get_std_unif_float_seq(nn,seedlist[3],hf);
    }else{
      exit_error("STIMX_BAR2RAN","Unknown shape");
    }

    t0i = my_rint(t0 / tscale);  // Time index at which stim starts

  }else if (task == 1){
    //
    //  Get stimulus frame
    //

    if ((xpix1 == NULL) && (xf1 == NULL))
      exit_error("STIMX_BAR2RAN","Not prepared to generate stim");

    // Index into random bar position arrays
    ri = (int)((ti-t0i)/(dwell*tsamp));
    if ((ri < 0) || (ri >= nn)){
      printf("  ti = %d\n",ti);
      printf("  dwell = %d\n",dwell);
      printf("  tsamp = %d\n",tsamp);
      printf("  nn = %d\n",nn);
      printf("  ri = %d\n",ri);
      exit_error("STIMX_BAR2RAN","index exceeds random sequence length");
    }


    if (strcmp(shape,"rect")==0){
      xr1 = x0 + xpix1[ri]; // light bar
      yr1 = y0 + ypix1[ri];
      xr0 = x0 + xpix0[ri]; // dark bar
      yr0 = y0 + ypix0[ri];
    }else if (strcmp(shape,"sine_bump")==0){
      xv1 = (float)x0 + xf1[ri]; // light bar
      yv1 = (float)y0 + yf1[ri];
      xv0 = (float)x0 + xf0[ri]; // dark bar
      yv0 = (float)y0 + yf0[ri];
    }
    
    //make_frame_flat(data,1,xn,1,yn,1,zn,i,bgval);
    // Make new flat frame
    tframe = get_2d_farray(xn,yn);
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	tframe[i][j] = bgval;


    if (strcmp(shape,"rect")==0){
      make_frame_add_rect_pix_0(tframe,xn,yn,xr1,yr1,bw,bh,val1);
      make_frame_add_rect_pix_0(tframe,xn,yn,xr0,yr0,bw,bh,val2);
    }else if (strcmp(shape,"sine_bump")==0){
      make_frame_add_sine_bump_0(tframe,xn,yn,xv1,yv1,diamf,val1);
      make_frame_add_sine_bump_0(tframe,xn,yn,xv0,yv0,diamf,val2);
    }

    *rframe = tframe;

  }else{
    exit_error("STIMX_BAR2RAN","unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               STIMX_BAR1RAN                               */
/*                                                                           */
/*****************************************************************************/
void stimx_bar1ran(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,
		   rframl,rframr)
     char *mylogf;                 // Log file
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // Time units per stim frame
     int ti;                       // Time index (tscale units)
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframl;              // Return 2D frame, left eye
     float ***rframr;              // Return 2D frame, right eye
{
  int i,j,k;
  int seed_l,seed_r,ri,bwflag,*tvall,*tvalr;
  float t0,tn,cx,cy,contrast,maxlum,fi,tdx,tdy,amp;
  float **tframl,**tframr;

  static int nbar,nn,t0i,dwell,*bpos_l,*bpos_r;
  static float barw,barh,theta;
  static float bgval,val0,val1,sdaa,*bx,*by,*bval_l,*bval_r;


  if (task == -1){  // Clean up
    if (bpos_l != NULL)  myfree(bpos_l);
    if (bval_l != NULL)  myfree(bval_l);
    if (bpos_r != NULL)  myfree(bpos_r);
    if (bval_r != NULL)  myfree(bval_r);
    if (bx     != NULL)  myfree(bx);
    if (by     != NULL)  myfree(by);
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //

    t0       = param_getf_exit(ppl,"st0");
    tn       = param_getf_exit(ppl,"stn");
    cx       = param_getf_exit(ppl,"cx");
    cy       = param_getf_exit(ppl,"cy");
    theta    = param_getf_exit(ppl,"direction");
    nbar     = param_geti_exit(ppl,"nbar");
    barw     = param_getf_exit(ppl,"barw");  // (deg)
    barh     = param_getf_exit(ppl,"barh");  // (deg)
    seed_l   = param_geti_exit(ppl,"seed");
    seed_r   = param_geti_exit(ppl,"seed_r");
    contrast = param_getf_exit(ppl,"contrast");
    bgval    = param_getf_exit(ppl,"bgval");
    maxlum   = param_getf_exit(ppl,"maxlum");
    dwell    = param_geti_exit(ppl,"dwell");
    bwflag   = param_geti_dflt(ppl,"black_white",0);
    sdaa     = param_getf_dflt(ppl,"sdaa",0.0);

    //
    //  Choose color values for the 'white' and 'black' bars
    //
    amp = (maxlum - bgval) * contrast;
    if (bwflag == 0){
      val0 = -amp;     // dark bar
      val1 = amp;      // white bar
    }else if (bwflag == -1){
      val0 = val1 = -amp;     // both dark
    }else if (bwflag == 1){
      val0 = val1 =  amp;     // both light
    }

    //
    //  Compute x,y center coordinates for each of 'nbar'
    //
    bx = (float *)myalloc(nbar*sizeof(float));
    by = (float *)myalloc(nbar*sizeof(float));
    tdx = barw * cos(M_PI*theta/180.0);
    tdy = barw * sin(M_PI*theta/180.0);
    for(i=0;i<nbar;i++){
      fi = (float)i - (float)(nbar - 1)/2.0;
      bx[i] = tdx * fi + cx;
      by[i] = tdy * fi + cy;
    }

    // Compute number of new random frames to create
    // WYETH ** this differs from stimgen_bar2ran, that way is probably better
    nn = stm_util_count_new_frames(zn,tscale,t0,tn,tsamp,dwell);

    //
    //  Get random bar positions, integer values from 0..nbar-1
    //
    bpos_l = myrand_get_std_unif_int_seq(nn,seed_l,nbar);
    bpos_r = myrand_get_std_unif_int_seq(nn,seed_r,nbar);

    //
    //  Get random bar values: 0-black, 1-white
    //
    tvall = myrand_get_std_bin_seq(nn,seed_l*3+1019);
    tvalr = myrand_get_std_bin_seq(nn,seed_r*3+1019);
    bval_l = (float *)myalloc(nn*sizeof(float));
    bval_r = (float *)myalloc(nn*sizeof(float));
    for(i=0;i<nn;i++){
      if (tvall[i] == 0)
	bval_l[i] = val0;
      else
	bval_l[i] = val1;
      if (tvalr[i] == 0)
	bval_r[i] = val0;
      else
	bval_r[i] = val1;
    }
    //printf("var0 = %f  val1= %f\n",val0,val1);
    myfree(tvall);
    myfree(tvalr);

    t0i = my_rint(t0 / tscale);  // Time index at which stim starts

  }else if (task == 1){
    //
    //  Get stimulus frame
    //

    if (bpos_l == NULL)
      exit_error("STIMX_BAR1RAN","Not prepared to generate stim");

    // Index into random bar position arrays
    ri = (int)((ti-t0i)/(dwell*tsamp));
    if ((ri < 0) || (ri >= nn)){
      printf("  ti = %d\n",ti);
      printf("  dwell = %d\n",dwell);
      printf("  tsamp = %d\n",tsamp);
      printf("  nn = %d\n",nn);
      printf("  ri = %d\n",ri);
      exit_error("STIMX_BAR1RAN","index exceeds random sequence length");
    }

    // Make new flat frame
    tframl = get_2d_farray(xn,yn);
    tframr = get_2d_farray(xn,yn);
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++){
	tframr[i][j] = bgval;
	tframl[i][j] = bgval;
      }

    k = bpos_l[ri];  // Index of bar position
    make_frame_bar_add_0(tframl,xn,yn,bx[k],by[k],theta,barh,barw,sdaa,
			 sscale,bval_l[ri],0);

    k = bpos_r[ri];  // Index of bar position
    make_frame_bar_add_0(tframr,xn,yn,bx[k],by[k],theta,barh,barw,sdaa,
			 sscale,bval_r[ri],0);
    *rframl = tframl;
    *rframr = tframr;
  }else{
    exit_error("STIMX_BAR1RAN","unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                STIMX_DIRMOD                               */
/*                                                                           */
/*****************************************************************************/
void stimx_dirmod(ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,rframl,rframr)
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // Time units per stim frame
     int ti;                       // Time index (tscale units), frame number
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframl;              // Return 2D frame into this address
     float ***rframr;              // Return 2D frame for Right Eye
{
  int i,j,k;
  int distrib,seed1,seed2,ord,si,seqn;
  float t0,tn,contrast,maxlum,**tframe,**tframe1,**tframe2,t,x,y,fps;
  float **tframe1r,**tframr;
  float etf1p,etf1a,etf2p,etf2a,phase1,phase2,size1,size2;
  float v1p,v1a,v2p,v2a,bgval,*seq1,*seq2,ph,tt;

  // Static
  static int aptype,nstim,seqcn,dt,t0i,dicho,dichopp;
  static float sf1,theta1,sz1,con1;
  static float sf2,theta2,sz2,con2;
  static float xc,yc,cx,cy,bgamp,*seq1cum,*seq2cum,mid,*seq1cumr = NULL;
  static char *stimform = NULL;

  /*
    printf("xn,yn,zn = %d %d %d\n",xn,yn,zn);
    printf("sscale = %f\n",sscale);
    printf("tscale = %f\n",tscale);
    printf("tsamp  = %d\n",tsamp);
    printf("ti     = %d\n",ti);
    printf("task   = %d\n",task);
  */

  stimform = param_getc_dflt(ppl,"stim_form","3d");
  if (strcmp(stimform,"3d_b")==0){
    // 'dicho' =  -1 - Not dichoptic (0-L, 1-R, 2-split, 3-Both)    
    dicho   = param_geti_exit(ppl,"dichoptic"); 
    dichopp = param_geti_dflt(ppl,"dicho_opp",0);
  }else
    dicho = -1;  // monocular, return left eye only

  if (task == -1){  // Clean up
    if (seq1cum  != NULL)  myfree(seq1cum);
    if (seq1cumr != NULL)  myfree(seq1cumr);
    if (seq2cum  != NULL)  myfree(seq2cum);
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //

    nstim = 1;  // Assume one stimulus, unless we find params for two

    cx = param_getf_dflt(ppl,"cx",0.0);
    cy = param_getf_dflt(ppl,"cy",0.0);

    if (param_test(ppl,"sf2")){  // WYETH - OLD WAY - discourage
      printf("  *** OLD WAY - Please convert to new param names:\n");
      printf("  sf2    -->  s2_sf\n");
      printf("  etf1a  -->  s1_etfa\n");
      printf("         ...\n");
    }else{
      if (param_test(ppl,"etf")){
	etf1p  = param_getf_exit(ppl,"etf");
	etf1a  = etf1p;
	etf2p  = etf1p;
	etf2a  = etf1a;
      }else if (param_test(ppl,"s1_etf")){
	param_getf_s2_exit(ppl,"etf",&etf1p,&etf2p,&nstim);
	etf1a = etf1p;
	etf2a = etf2p;
      }else{
	param_getf_s2_exit(ppl,"etfp",&etf1p,&etf2p,&nstim);
	param_getf_s2_exit(ppl,"etfa",&etf1a,&etf2a,&nstim);
      }
      param_getf_s2_exit(ppl,"sf",&sf1,&sf2,&nstim);
      param_getf_s2_exit(ppl,"phase",&phase1,&phase2,&nstim);
      theta1   = param_getf_exit(ppl,"direction");
      param_getf_s2_exit(ppl,"contrast",&con1,&con2,&nstim);
      param_getf_s2_exit(ppl,"size",&size1,&size2,&nstim);
      param_geti_s2_exit(ppl,"seed",&seed1,&seed2,&nstim);
    }

    /*
    printf("sf1,2    = %f  %f\n",sf1,sf2);
    printf("etf1a,p  = %f  %f\n",etf1a,etf1p);
    printf("etf2a,p  = %f  %f\n",etf2a,etf2p);
    printf("con1,2   = %f  %f\n",con1,con2);
    printf("seed1,2  = %d  %d\n",seed1,seed2);
    printf("size1,2  = %f  %f\n",size1,size2);
    printf("phase1,2  = %f  %f\n",phase1,phase2);
    */


    maxlum   = param_getf_dflt(ppl,"maxlum",1.0);
    aptype   = param_geti_exit(ppl,"aptype");
    bgval    = param_getf_exit(ppl,"bgval");
    bgamp = bgval * 2.0 - 1.0; // Convert 0..1 --> -1..1

    theta2 = theta1 + 90.0;
    mid = (float)maxlum/2.0;

    sz1 = size1/sscale;
    sz2 = size2/sscale;

    //
    // Random Sequences
    //
    t0       = param_getf_exit(ppl,"st0");
    tn       = param_getf_exit(ppl,"stn");
    distrib  = param_geti_exit(ppl,"distrib");
    dt       = param_geti_exit(ppl,"dt");

    if (distrib == 7){
      ord = param_geti_exit(ppl,"mseq_ord");
    }else
      ord = -1;

    // Set float values for random sequences
    fps = 1.0/(tsamp*tscale);  // Frames/s - Pattern frames (not blanks)
    v1p =  etf1p/fps * 360.0;  // step size in degrees of spatial phase
    v1a = -etf1a/fps * 360.0;
    v2p =  etf2p/fps * 360.0;
    v2a = -etf2a/fps * 360.0;

    // Return a sequence of random values (no repeats for dt, minimal seq)
    stm_get_ran_seq(NULL,distrib,ord,seed1,fps,tn,dt,v1p,v1a,0.0,&seq1,&seqn);

    // Create cumulative phase sequence
    //
    // *** NOTE:  in the old accumulation algorithm, the first step occurred
    //     before the stimulus is shown (old 'dmmask' stim).  This is 
    //     ********* NOT TRUE ANY MORE, thus this NOW DIFFERS FROM OLD VERSION
    //
    seqcn = seqn*dt+1;    // Use each jump 'dt' times, +1 for ph0
    seq1cum  = get_farray(seqcn);
    seq1cumr = get_farray(seqcn);
    seq1cum[0]  = phase1;
    seq1cumr[0] = phase1;
    k = 1;
    for(i=0;i<seqn;i++){
      for(j=0;j<dt;j++){
	seq1cum[k]  = seq1cum[k-1]  + seq1[i];
	if (dichopp == 1){ // Make RE move opposite
	  seq1cumr[k] = seq1cumr[k-1] - seq1[i];
	}else
	  seq1cumr[k] = seq1cumr[k-1] + seq1[i];
	k += 1;
      }
    }
    myfree(seq1);

    if (nstim == 2){
      stm_get_ran_seq(NULL,distrib,ord,seed2,fps,tn,dt,v2p,v2a,0.0,&seq2,&seqn);

      // Create cumulative phase sequence
      seq2cum = get_farray(seqcn);
      seq2cum[0] = phase2;
      k = 1;
      for(i=0;i<seqn;i++){
	for(j=0;j<dt;j++){
	  seq2cum[k] = seq2cum[k-1] + seq2[i];
	  k += 1;
	}
      }
      myfree(seq2);
    }

    xc = (float)(xn-1)/2.0 + cx/sscale;  // Center w.r.t. 0..xn-1 (pix)
    yc = (float)(yn-1)/2.0 + cy/sscale;  // Center w.r.t. 0..yn-1 (pix)

    // DO WE NEED THIS?
    t0i = my_rint(t0 / tscale);  // Time index at which stim starts

  }else if (task == 1){
    //
    //  Get stimulus frame
    //

    if (seq1cum == NULL)
      exit_error("STIMX_DIRMOD","Not prepared to generate stim");

    //
    // Compute the index into the random bar sequence to use for this frame
    //
    si = (ti-t0i)/tsamp;  // Number of pattern frames from stim start


    // Make new flat frame
    tframe1 = get_2d_farray(xn,yn);
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	tframe1[i][j] = bgamp;

    if (dicho >= 0){
      //
      //  Make a second frame
      //
      tframe1r = get_2d_farray(xn,yn);
      for(i=0;i<xn;i++)
	for(j=0;j<yn;j++)
	  tframe1r[i][j] = bgamp;
    }

    // Make new flat frame
    if (nstim == 2){
      tframe2 = get_2d_farray(xn,yn);
      for(i=0;i<xn;i++)
	for(j=0;j<yn;j++)
	  tframe2[i][j] = bgamp;
      //
      //  *** WYETH - I have not yet implemented dichoptic for the 2nd patch
      //
      if (dicho >= 0)
	exit_error("STIMX_DIRMOD","'dichoptic' not implemented for 2nd patch");
    }

    if ((si >= 0) && (si < seqcn)){
      ph = seq1cum[si];
      if ((dicho == -1)||(dicho == 0)||(dicho == 3))  // LE only, or both eyes
	make_frame_sinusoid_0(tframe1,xn,yn,sf1,theta1,sscale,ph,cx,cy,1.0);
      if ((dicho == 1)||(dicho == 3)){
	// Same thing in R.Eye
	ph = seq1cumr[si];
	make_frame_sinusoid_0(tframe1r,xn,yn,sf1,theta1,sscale,ph,cx,cy,1.0);
      }
      
      if (nstim == 2){
	ph = seq2cum[si];
	make_frame_sinusoid_0(tframe2,xn,yn,sf2,theta2,sscale,ph,cx,cy,1.0);
      }
    }

    // Stimulus values are -1...1 up to now
    apply_aperture_scale(tframe1,xn,yn,aptype,xc,yc,sz1,sz1,0.0,bgamp,
			 con1*mid);
    if (dicho >= 0){
      apply_aperture_scale(tframe1r,xn,yn,aptype,xc,yc,sz1,sz1,0.0,bgamp,
			   con1*mid);
    }
    
    if (nstim == 2){
      apply_aperture_scale(tframe2,xn,yn,aptype,xc,yc,sz2,sz2,0.0,bgamp,
			   con2*mid);
    }

    if (nstim == 1){
      tframe = tframe1;
      tframr = tframe1r;
    }else if (nstim == 2){
      tframe = add_2d_farrays(tframe1,tframe2,xn,yn);
      // Values still centered on 0, but could exceed -mid...mid
      // Restrict stimulus values to -mid...mid
      for(i=0;i<xn;i++){
	for(j=0;j<yn;j++){
	  tt = tframe[i][j];
	  if (tt < -mid)
	    tframe[i][j] = -mid;
	  else if (tt > mid)
	    tframe[i][j] = mid;
	}
      }
    }

    // Final values should be in 0..maxlum
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	tframe[i][j] += mid;

    if (dicho >= 0){
      for(i=0;i<xn;i++)
	for(j=0;j<yn;j++)
	  tframr[i][j] += mid;
    }

    if (nstim == 2){
      free_2d_farray(tframe1,xn);
      free_2d_farray(tframe2,xn);
    }

    *rframl = tframe;
    if (dicho >= 0){
      *rframr = tframr;
    }
    
  }else{
    exit_error("STIMX_DIRMOD","unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                STIMX_PLAID                                */
/*                                                                           */
/*  Sum of two gratings.                                                     */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
void stimx_plaid(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,
		 rframl,rframr)
     char *mylogf;
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // Time units per stim frame
     int ti;                       // Time index (tscale units)
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframl;              // Return 2D frame, left eye
     float ***rframr;              // Return 2D frame, right eye
{
  int i,j;
  int tpseudi;
  float ph1,ph2,dir,theta1,theta2,t,size,th1t,th2t;
  float **tframe,**tframe_r,**tframe1,**tframe2,**tframe3,**tframe4;
  char *stimform;
  static int aptype,pldflag,dicho,meanzero,t0i,tpseudo,tpseud_lval,tpseud_rval;
  static int phdispfl;
  static float cx1,cx2,cy1,cy2,sf1,sf2,th1,th2,tf1,tf2,phase1,phase2;
  static float xc1,yc1,xc2,yc2,phdisp1,phdisp2;
  static float t0,tn,t0_2,tn_2,con1,con2,mid,sz,v0v,v1v,bgval,bgamp,maxlum;

  if (task == -1){  // Clean up
    ;
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //
    param_getf_2(ppl,"cx","cx1","cx2",&cx1,&cx2);
    param_getf_2(ppl,"cy","cy1","cy2",&cy1,&cy2);
    t0       = param_getf_exit(ppl,"st0");
    tn       = param_getf_exit(ppl,"stn");
    t0_2     = param_getf_dflt(ppl,"st0_2",t0);
    tn_2     = param_getf_dflt(ppl,"stn_2",tn);
    dir      = param_getf_exit(ppl,"direction");
    size     = param_getf_exit(ppl,"size");
    maxlum   = param_getf_exit(ppl,"maxlum");
    aptype   = param_geti_exit(ppl,"aptype");
    pldflag  = param_geti_dflt(ppl,"plaid_flag",0);
    bgval    = param_getf_exit(ppl,"bgval");
    stimform = param_getc_dflt(ppl,"stim_form","3d");
    meanzero = param_geti_dflt(ppl,"stim_mean_zero",0);
    phdisp1  = param_getf_dflt(ppl,"ph_disp_1",0.0);
    phdisp2  = param_getf_dflt(ppl,"ph_disp_2",0.0);
    phdispfl = param_geti_dflt(ppl,"phase_disp_flag",1);
    tpseudo  = param_geti_dflt(ppl,"t_pseudo",0);
    tpseudi  = param_geti_dflt(ppl,"t_pseudo_dich",0);

    param_getf_2(ppl,"sf","sf1","sf2",&sf1,&sf2);
    param_getf_2(ppl,"tf","tf1","tf2",&tf1,&tf2);
    param_getf_2(ppl,"phase","phase1","phase2",&phase1,&phase2);
    param_getf_2(ppl,"theta","theta1","theta2",&theta1,&theta2);
    param_getf_2(ppl,"contrast","contrast1","contrast2",&con1,&con2);

    if (strcmp(stimform,"3d_b")==0)
      dicho = param_geti_exit(ppl,"dichoptic");
    else
      dicho = -1;  // monocular, return left eye only

    th1 = dir + theta1;
    th2 = dir + theta2;

    if (meanzero == 0){
      mid = maxlum/2.0;
      bgamp = bgval * 2.0 - 1.0; // Convert 0..1 --> -1..1

      v0v = maxlum/2.0;  // The mean of the grating
      v1v = maxlum/2.0;  // The amplitude of a 100% contrast grating
    }else{
      // WYETH 2014 Nov - I'm trying these new def's for mean zero
      v0v = 0.0;      // The mean of the grating
      v1v = 1.0;      // The amplitude of a 100% contrast grating
      bgamp = bgval;  // 'bgval' should already be in -1..1
      //printf("WYETH HERE ___________ bgval = %f\n",bgval);
    }

    // WYETHWYETH - T-alt
    t0i = stm_get_ti0(zn,tscale,t0,tsamp);  // Time index at which stim starts

    if (tpseudi == 0){
      tpseud_lval = 0;
      tpseud_rval = 1;
    }else{
      tpseud_lval = 1;
      tpseud_rval = 0;
    }
      

    sz = size/sscale;

    myfree(stimform);

  }else if (task == 1){
    //
    //  Get stimulus frame
    //

    // Make new flat frame
    tframe1 = get_2d_farray(xn,yn);
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	tframe1[i][j] = bgamp;

    // Make new flat frame
    tframe2 = get_2d_farray(xn,yn);
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	tframe2[i][j] = bgamp;

    xc1 = (float)(xn-1)/2.0 + cx1/sscale;  // Center w.r.t. 0..xn-1 (pix)
    yc1 = (float)(yn-1)/2.0 + cy1/sscale;  // Center w.r.t. 0..yn-1 (pix)
    xc2 = (float)(xn-1)/2.0 + cx2/sscale;  // Center w.r.t. 0..xn-1 (pix)
    yc2 = (float)(yn-1)/2.0 + cy2/sscale;  // Center w.r.t. 0..yn-1 (pix)

    t = (float)ti*tscale;


    // Was 'i' when i=1 on first frame, so now off by 1?
    ph1 = phase1 + (float)ti*tscale*tf1 * 360.0;
    ph2 = phase2 + (float)ti*tscale*tf2 * 360.0;

    if ((t >= t0) && (t < t0+tn)){
      if ((tpseudo < 1) ||
	  (((ti - t0i) / tpseudo) % 2 == 0)){
	make_frame_sinusoid_0(tframe1,xn,yn,sf1,th1,sscale,ph1,cx1,cy1,
			      con1*v1v);
      }
    }

    if ((t >= t0_2) && (t < t0_2+tn_2)){
      if ((tpseudo < 1) ||
	  (((ti - t0i) / tpseudo) % 2 == 1)){
	make_frame_sinusoid_0(tframe2,xn,yn,sf2,th2,sscale,ph2,cx2,cy2,
			      con2*v1v);
      }
    }


    if ((dicho == 4) || (dicho == 5)){
      //
      //  Allow 4 gratings (plaids in both eyes)
      //    4 - with phase disparity shifts for both
      //    5 - opposite direction in both eyes
      //
      //tframe3 = get_2d_farray(xn,yn);
      //tframe4 = get_2d_farray(xn,yn);

      // Make new flat frame
      tframe3 = get_2d_farray(xn,yn);
      for(i=0;i<xn;i++)
	for(j=0;j<yn;j++)
	  tframe3[i][j] = bgamp;

      // Make new flat frame
      tframe4 = get_2d_farray(xn,yn);
      for(i=0;i<xn;i++)
	for(j=0;j<yn;j++)
	  tframe4[i][j] = bgamp;


      //printf(" %d %d\n",tpseud_lval,tpseud_rval);

      th1t = th1;
      th2t = th2;
      if (dicho == 5){
	th1t += 180.0;
	th2t += 180.0;
      }

      if ((t >= t0) && (t < t0+tn)){
	if (phdispfl == 1){
	  // Value is interpreted as max horiz phase shift (deg)
	  //   and is adjusted according to the ori/direction of the grating
	  ph1 += phdisp1 * cos(th1t * M_PI/180.0);
	  //printf("  Theta:  %f   Ph_disp  %f\n",
	  //   th1t,phdisp1*cos(th1t*M_PI/180.0));
	}else
	  // Value is interpreted as raw, orthogonal phase shift, regardless
	  //   of direction.  ***NOTE horiz. disparity will change w/ dir!
	  ph1 += phdisp1;

	if ((tpseudo < 1) ||
	    (((ti - t0i) / tpseudo) % 2 == tpseud_lval)){
	  make_frame_sinusoid_0(tframe3,xn,yn,sf1,th1t,sscale,ph1,cx1,cy1,
				con1*v1v);
	}
      }

      if ((t >= t0_2) && (t < t0_2+tn_2)){
	if (phdispfl == 1)
	  // Value is interpreted as max horiz phase shift (deg)
	  //   and is adjusted according to the ori/direction of the grating
	  ph2 += phdisp2 * cos(th2t * M_PI/180.0);
	else
	  // Value is interpreted as raw, orthogonal phase shift, regardless
	  //   of direction.  ***NOTE horiz. disparity will change w/ dir!
	  ph2 += phdisp2;

	if ((tpseudo < 1) ||
	    (((ti - t0i) / tpseudo) % 2 == tpseud_rval)){
	  make_frame_sinusoid_0(tframe4,xn,yn,sf2,th2t,sscale,ph2,cx2,cy2,
				con2*v1v);
	}
      }
      if (pldflag == 2){  // Binarize the existing frame
	make_frame_binarize_0(tframe3,xn,yn,0.0,con2*v1v,-con1*v1v);
	make_frame_binarize_0(tframe4,xn,yn,0.0,con2*v1v,-con2*v1v);
      }
    }


    // Stim values are now centered around zero

    if (pldflag == 2){  // Binarize the existing frame
      make_frame_binarize_0(tframe1,xn,yn,0.0,con1*v1v,-con1*v1v);
      make_frame_binarize_0(tframe2,xn,yn,0.0,con2*v1v,-con2*v1v);
    }

    if (dicho == 2){
      //
      //  Keep the plaid split between two eyes
      //
      add_const_2d_farray(tframe1,0,xn,0,yn,v0v);
      add_const_2d_farray(tframe2,0,xn,0,yn,v0v);

      apply_aperture_scale(tframe1,xn,yn,aptype,xc1,yc1,sz,sz,0.0,bgval,1.0);
      apply_aperture_scale(tframe2,xn,yn,aptype,xc2,yc2,sz,sz,0.0,bgval,1.0);

      *rframl = tframe1;
      *rframr = tframe2;
    }else if ((dicho == 4) || (dicho == 5)){
      //
      //  WYETH - NEW Jun 3, 2016 (dicho 5 added Oct 17)
      //

      tframe   = add_2d_farrays(tframe1,tframe2,xn,yn);
      tframe_r = add_2d_farrays(tframe3,tframe4,xn,yn);

      free_2d_farray(tframe1,xn);
      free_2d_farray(tframe2,xn);
      free_2d_farray(tframe3,xn);
      free_2d_farray(tframe4,xn);

      add_const_2d_farray(tframe  ,0,xn,0,yn,v0v); // Adjust if meanzero = 0
      add_const_2d_farray(tframe_r,0,xn,0,yn,v0v); // Adjust if meanzero = 0

      apply_aperture_scale(tframe  ,xn,yn,aptype,xc1,yc1,sz,sz,0.0,bgval,1.0);
      apply_aperture_scale(tframe_r,xn,yn,aptype,xc2,yc2,sz,sz,0.0,bgval,1.0);

      *rframl = tframe;
      *rframr = tframe_r;

    }else{

      // WYETH - NEW Sept 13, 2015
      //
      //  Must use '0.0' for background val, so backgrounds do not end up
      //  adding below.
      //
      //  Both frames have 0-mean range at this point
      //
      //  Note, I think the aperture is applied at this point in case the
      //  gratings are not spatially overlapping.
      //
      apply_aperture_scale(tframe1,xn,yn,aptype,xc1,yc1,sz,sz,0.0,0.0,1.0);
      apply_aperture_scale(tframe2,xn,yn,aptype,xc2,yc2,sz,sz,0.0,0.0,1.0);

      if (pldflag == 3) // `AND'
	tframe = and_2d_farrays(tframe1,tframe2,xn,yn,v1v,-v1v);
      else if (pldflag == 4) // `NAND'
	tframe = and_2d_farrays(tframe1,tframe2,xn,yn,-v1v,v1v);
      else{
	tframe = add_2d_farrays(tframe1,tframe2,xn,yn);
	//tframe = copy_2d_farray(tframe2,xn,yn);
	//print_min_max_2d_farray(tframe1,xn,yn,"tframe1");
	//print_min_max_2d_farray(tframe2,xn,yn,"tframe2");
      }

      free_2d_farray(tframe1,xn);
      free_2d_farray(tframe2,xn);

      add_const_2d_farray(tframe,0,xn,0,yn,v0v); // Adjust if meanzero = 0

      if ((xc1 == xc2) && (yc1 == yc2)){
	//
	//  WYETH - We generally cannot apply an aperture here if the
	//          the two patches are not perfectly overlaid.
	//
	//  Note, if "meanzero = 0", we're using 0..1 range, and so we need
	//        to add back the original 'bgval', because 'bgamp' will be 0.
	//
	apply_aperture_scale(tframe,xn,yn,aptype,xc1,yc1,sz,sz,0.0,bgval,1.0);
      }

      // Truncate values that are out of range
      if (meanzero == 0)
	stim_util_limit_2d(tframe,xn,yn,0.0,maxlum);
      else
	stim_util_limit_2d(tframe,xn,yn,-1.0,1.0);

      *rframl = tframe;  // Same in both eyes, or all in one eye
    }
  }else{
    exit_error("STIMX_PLAID","Unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              STIMX_GABOR_GRID                             */
/*                                                                           */
/*  Fields of Gabor patches in both eyes.                                    */
/*  Designed to capture stimulus of Rokers, Czuba, Cormack and Huk (2011)    */
/*                                                                           */
/*****************************************************************************/
void stimx_gabor_grid(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,
		      rframl,rframr)
     char *mylogf;
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // Time units per stim frame
     int ti;                       // Time index (tscale units)
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframl;              // Return 2D frame, left eye
     float ***rframr;              // Return 2D frame, right eye
{
  int i,j,k;
  int gxn,gyn,seed,regflag,*pdir;
  float boxw,boxh,cx,cy,*px,*py,*pph1,*pph2,dirbal;
  float ph1,dir,tdir,theta,thopp,theta1,theta2,t,amp,amp1,amp2,sdp,sdo;
  float **tframe1,**tframe2;
  char *gtype,*monoflag,*dichopt;
  static int meanzero,gpn_r,gpn_l,npatch,make_left,make_right;
  static float cx1,cy1,sf1,th1,th2,tf1;
  static float t0,tn,con1,con2,sz,v0v,v1v,bgval,bgamp,maxlum;
  static char *pchtype;
  static struct stmh_gpatch **gp_r,**gp_l;

  if (task == -1){  // Clean up
    ;
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //

    t0       = param_getf_exit(ppl,"st0");
    tn       = param_getf_exit(ppl,"stn");
    cx1      = param_getf_dflt(ppl,"cx",0.0);
    cy1      = param_getf_dflt(ppl,"cy",0.0);
    dir      = param_getf_exit(ppl,"direction");
    param_getf_2(ppl,"theta","theta1","theta2",&theta1,&theta2);
    dirbal   = param_getf_dflt(ppl,"dir_balance",0.5);
    sf1      = param_getf_exit(ppl,"sf");
    tf1      = param_getf_exit(ppl,"tf");
    pchtype  = param_getc_dflt(ppl,"patch_type","gabor");
    if (strcmp(pchtype,"sine")==0){
      sdo      = param_getf_exit(ppl,"size");
      sdp      = param_getf_exit(ppl,"taper_len");
    }else{
      sdo      = param_getf_exit(ppl,"sd_orth");
      sdp      = param_getf_exit(ppl,"sd_par");
    }
    gtype    = param_getc_exit(ppl,"grid_type");
    dichopt  = param_getc_dflt(ppl,"dichoptic","opposite");
    gxn      = param_geti_exit(ppl,"grid_w");
    gyn      = param_geti_exit(ppl,"grid_h");
    boxw     = param_getf_exit(ppl,"box_w");
    boxh     = param_getf_exit(ppl,"box_h");
    seed     = param_geti_exit(ppl,"seed");
    if (strcmp(gtype,"regular")==0){
      npatch = gxn*gyn;
      regflag = 1;
    }else
      npatch = param_geti_exit(ppl,"n_patch");  // -1 indicates full

    //con1     = param_getf_exit(ppl,"contrast");
    param_getf_2(ppl,"contrast","contrast1","contrast2",&con1,&con2);

    bgval    = param_getf_exit(ppl,"bgval");
    maxlum   = param_getf_exit(ppl,"maxlum");
    meanzero = param_geti_dflt(ppl,"stim_mean_zero",0);
    monoflag = param_getc_dflt(ppl,"stim_monocular","both"); // "left","right"

    // WYETH phase NOT NEEDED UNTIL WE HAVE A non-random phase setting.
    //param_getf_2(ppl,"phase","phase1","phase2",&phase1,&phase2);

    make_left = make_right = 1;
    if (strcmp(monoflag,"left")==0)
      make_right = 0;
    else if (strcmp(monoflag,"right")==0)
      make_left = 0;
    myfree(monoflag);

    if (npatch > (gxn * gyn))
      exit_error("STIMX_GABOR_GRID","Value for 'n_patch' is too high");
    else if (npatch < 0)
      npatch = gxn * gyn;

    th1 = dir + theta1;
    th2 = dir + theta2;

    if (meanzero == 0){
      bgamp = bgval * 2.0 - 1.0; // Convert 0..1 --> -1..1

      v0v = maxlum/2.0;  // The mean of the grating
      v1v = maxlum/2.0;  // The amplitude of a 100% contrast grating
    }else{
      // WYETH 2014 Nov - I'm trying these new def's for mean zero
      v0v = 0.0;      // The mean of the grating
      v1v = 1.0;      // The amplitude of a 100% contrast grating
      bgamp = bgval;  // 'bgval' should already be in -1..1
    }

    //
    //  Create g-patch lists
    //
    if ((strcmp(gtype,"simple")==0) || (strcmp(gtype,"regular")==0)){

      stm_gpatch_config_simp(gxn,gyn,boxw,boxh,seed,npatch,dirbal,regflag,
			     &px,&py,&pdir,&pph1,&pph2);
      gpn_l = gpn_r = npatch;
    }else{
      printf("  grid_type = %s\n",gtype);
      exit_error("STIMX_GABOR_GRID","Unknown grid type");
    }

    gp_l = (struct stmh_gpatch **)myalloc(gpn_l*sizeof(struct stmh_gpatch *));
    gp_r = (struct stmh_gpatch **)myalloc(gpn_r*sizeof(struct stmh_gpatch *));

    //
    //  Fill g-patch lists
    //
    amp1  = con1 * v1v;  // Amplitude depends on stimulus range (0..1 or -1..1)
    amp2  = con2 * v1v;  // Amplitude depends on stimulus range (0..1 or -1..1)

    for(i=0;i<npatch;i++){

      if (pdir[i] == 1){
	tdir = th1;
	amp = amp1;
      }else{
	tdir = th2;
	amp = amp2;
      }
      if (strcmp(dichopt,"opposite")==0)
	thopp = tdir + 180.0;
      else if (strcmp(dichopt,"same")==0)
	thopp = tdir;
      else
	exit_error("STIMX_GABOR_GRID",
		   "Unknown 'dichoptic' value: 'same' or 'opposite'");

      cx = px[i] + cx1;
      cy = py[i] + cy1;
      gp_l[i] = stm_gpatch_get(cx,cy,tdir ,sf1,tf1,sdp,sdo,pph1[i],amp);
      gp_r[i] = stm_gpatch_get(cx,cy,thopp,sf1,tf1,sdp,sdo,pph2[i],amp);
    }

    myfree(gtype);

  }else if (task == 1){
    //
    //  Get stimulus frame
    //

    tframe1 = get_zero_2d_farray(xn,yn);
    tframe2 = get_zero_2d_farray(xn,yn);

    t = (float)ti*tscale;

    // Was 'i' when i=1 on first frame, so now off by 1?
    ph1 = (float)ti*tscale*tf1 * 360.0;  // Phase offset

    if ((t >= t0) && (t < t0+tn)){

      if (make_left == 1){
	for(i=0;i<gpn_l;i++){
	  if (strcmp(pchtype,"sine")==0){
	    make_frame_add_gpatch_sine(tframe1,xn,yn,gp_l[i],sscale,ph1);
	  }else{
	    make_frame_add_gpatch_gabor(tframe1,xn,yn,gp_l[i],sscale,ph1);
	  }
	}
      }
      if (make_right == 1){
	for(i=0;i<gpn_r;i++)
	  if (strcmp(pchtype,"sine")==0){
	    make_frame_add_gpatch_sine(tframe2,xn,yn,gp_r[i],sscale,ph1);
	  }else{
	    make_frame_add_gpatch_gabor(tframe2,xn,yn,gp_r[i],sscale,ph1);
	  }
      }
    }

    add_const_2d_farray(tframe1,0,xn,0,yn,v0v); // Adjust if meanzero = 0
    add_const_2d_farray(tframe2,0,xn,0,yn,v0v); // Adjust if meanzero = 0

    *rframl = tframe1;
    *rframr = tframe2;
  }else{
    exit_error("STIMX_GABOR_GRID","Unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             STIMX_GLOBMO_PATCH                            */
/*                                                                           */
/*  Global and local motion of Gabor patches.                                */
/*                                                                           */
/*****************************************************************************/
void stimx_globmo_patch(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,
		      rframl,rframr)
     char *mylogf;
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // Time units per stim frame
     int ti;                       // Time index (tscale units)
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframl;              // Return 2D frame, left eye
     float ***rframr;              // Return 2D frame, right eye
{
  int i;
  int patch_env,seed,*seedlist,pixw;
  float cx,cy,phase,ph,dir,t,tt,amp,sdp,sdo,cx1,cy1,sf;
  float glob_dx,glob_dt,glob_dir,con,gdy,gdx,tdist,cx0,cy0;
  float **tframe1,**tframe2;
  static int meanzero,npatch,make_left,make_right,draw_left,draw_right;
  static float t0,tn,tf,v0v,v1v,bgval,bgamp,maxlum,patch_dur,*pt0;
  static char *ptype = NULL;
  static struct stmh_gpatch **gp_r,**gp_l;

  if (task == -1){  // Clean up
    if (ptype != NULL)
      myfree(ptype);
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //

    t0        = param_getf_exit(ppl,"st0");
    tn        = param_getf_exit(ppl,"stn");
    cx1       = param_getf_dflt(ppl,"cx",0.0);
    cy1       = param_getf_dflt(ppl,"cy",0.0);

    //  Global motion parameters
    npatch    = param_geti_exit(ppl,"patch_n");
    patch_dur = param_getf_exit(ppl,"patch_dur");  // (s)
    patch_env = param_geti_exit(ppl,"patch_env");  // 0-no envelope
    glob_dx   = param_getf_exit(ppl,"glob_dx");    // (deg)
    glob_dt   = param_getf_exit(ppl,"glob_dt");    // (s)
    glob_dir  = param_getf_exit(ppl,"glob_dir");   // (deg)

    //  Gabor and local motion parameters
    ptype     = param_getc_dflt(ppl,"patch_type","gabor");
    if (strcmp(ptype,"gabor")==0){
      sf      = param_getf_exit(ppl,"sf");
      tf      = param_getf_exit(ppl,"tf");
      phase   = param_getf_exit(ppl,"phase");
    }else if (strcmp(ptype,"noise")==0){
      seed    = param_geti_exit(ppl,"seed");
      pixw    = param_geti_exit(ppl,"noise_pix");
      seedlist = get_seeds(seed,1000000,npatch);
    }
    dir       = param_getf_exit(ppl,"direction");
    sdo       = param_getf_exit(ppl,"sd_orth");
    sdp       = param_getf_exit(ppl,"sd_par");
    con       = param_getf_exit(ppl,"contrast");
    bgval     = param_getf_exit(ppl,"bgval");
    maxlum    = param_getf_exit(ppl,"maxlum");
    meanzero  = param_geti_dflt(ppl,"stim_mean_zero",0);

    stimu_std_make_draw(ppl,&make_left,&make_right,&draw_left,&draw_right);

    if (meanzero == 0){
      bgamp = bgval * 2.0 - 1.0; // Convert 0..1 --> -1..1

      v0v = maxlum/2.0;  // The mean of the grating
      v1v = maxlum/2.0;  // The amplitude of a 100% contrast grating
    }else{
      // WYETH 2014 Nov - I'm trying these new def's for mean zero
      v0v = 0.0;      // The mean of the grating
      v1v = 1.0;      // The amplitude of a 100% contrast grating
      bgamp = bgval;  // 'bgval' should already be in -1..1
    }

    //
    //  Create g-patch lists
    //
    gp_l = (struct stmh_gpatch **)myalloc(npatch*sizeof(struct stmh_gpatch *));
    if (draw_right)
      gp_r = (struct stmh_gpatch **)myalloc(npatch*
					    sizeof(struct stmh_gpatch *));

    //
    //  Fill g-patch lists
    //
    amp  = con * v1v;  // Amplitude depends on stimulus range (0..1 or -1..1)
    ph = phase;

    tdist = (float)(npatch-1) * glob_dx/2.0;   // distance away from midpoint
    gdx = glob_dx * cos(glob_dir*M_PI/180.0);
    gdy = glob_dx * sin(glob_dir*M_PI/180.0);

    cx0 = cx1 - 0.5 * (float)(npatch-1) * gdx;
    cy0 = cy1 - 0.5 * (float)(npatch-1) * gdy;

    pt0 = get_farray(npatch);

    cx = cx0;
    cy = cy0;
    for(i=0;i<npatch;i++){

      gp_l[i] = stm_gpatch_get(cx,cy,dir,sf,tf,sdp,sdo,ph,amp);
      if (strcmp(ptype,"noise")==0){
	gp_l[i]->ptype = 2;
	gp_l[i]->seed = seedlist[i];
	gp_l[i]->pixw = pixw;
      }
      if (draw_right){
	gp_r[i] = stm_gpatch_get(cx,cy,dir,sf,tf,sdp,sdo,ph,amp);
	if (strcmp(ptype,"noise")==0){
	  gp_r[i]->ptype = 2;
	  gp_r[i]->seed = seedlist[i];
	  gp_r[i]->pixw = pixw;
	}
      }

      pt0[i] = t0 + (float)i * glob_dt;
      //printf("pt0[i] = %f\n",pt0[i]);

      cx += gdx;
      cy += gdy;
    }

  }else if (task == 1){
    //
    //  Get stimulus frame
    //

    tframe1 = get_zero_2d_farray(xn,yn);
    if (make_right)
      tframe2 = get_zero_2d_farray(xn,yn);

    t = (float)ti*tscale;  // Time in sec

    // Was 'i' when i=1 on first frame, so now off by 1?
    ph = (float)ti*tscale*tf * 360.0;  // Phase offset

    if ((t >= t0) && (t < t0+tn)){

      for(i=0;i<npatch;i++){

	tt = t - pt0[i];  // Time relative to this patch start
	if ((tt >= 0.0) && (tt < patch_dur)){

	  //ph = (float)ti*tscale*tf * 360.0;  // Phase offset
	  ph = tt*tf * 360.0;  // Phase offset

	  if (draw_left){
	    if (strcmp(ptype,"noise")==0){
	      make_frame_add_gpatch_noise(tframe1,xn,yn,gp_l[i],sscale,ph);
	    }else{
	      make_frame_add_gpatch_gabor(tframe1,xn,yn,gp_l[i],sscale,ph);
	    }
	  }
	  if (draw_right){
	    if (strcmp(ptype,"noise")==0){
	      make_frame_add_gpatch_noise(tframe1,xn,yn,gp_r[i],sscale,ph);
	    }else{
	      make_frame_add_gpatch_gabor(tframe2,xn,yn,gp_r[i],sscale,ph);
	    }
	  }
	}
      }
    }

    add_const_2d_farray(tframe1,0,xn,0,yn,v0v); // Adjust if meanzero = 0
    if (make_right)
      add_const_2d_farray(tframe2,0,xn,0,yn,v0v); // Adjust if meanzero = 0

    *rframl = tframe1;
    if (make_right)
      *rframr = tframe2;
  }else{
    exit_error("STIMX_GLOBMO_PATCH","Unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                STIMX_DOTS                                 */
/*                                                                           */
/*****************************************************************************/
void stimx_dots(ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,rframe)
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // Time units per stim frame
     int ti;                       // Time index (tscale units), frame number
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframe;              // Return 2D frame into this address
{
  int i,j,k;
  int dpf,dtf,aptype,seed1,seed2;
  float **tframe,x0,y0,fps,szx,szy,xoff,yoff;
  float t0,tn,theta,theta1,theta2,coh,speed1,speed2,ddens1,ddens2,xf,cx,cy;

  // Static
  static int t0i,nframes,dualflag;
  static int *dotn,*dotn2;
  static float ampl1,ampl2,bgval,minlum,maxlum,dsize,size;
  static float **dotx,**doty,**dotx2,**doty2;

  if (task == -1){  // Clean up
    if (dotx != NULL)
      free_2d_farray_ignore_null(dotx,nframes); // In case 0 dots
    if (doty != NULL)
      free_2d_farray_ignore_null(doty,nframes);
    if (dotn != NULL)
      myfree(dotn);

    if (dualflag == 1){
      if (dotx2 != NULL)
	free_2d_farray_ignore_null(dotx2,nframes);
      if (doty2 != NULL)
	free_2d_farray_ignore_null(doty2,nframes);
      if (dotn2 != NULL)
	myfree(dotn2);
    }
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //

    if (param_test(ppl,"theta2")){
      //
      //  There are two dot fields - Developed for "dot plaids"
      //
      param_getf_2(ppl,"theta","theta1","theta2",&theta1,&theta2);
      param_getf_2(ppl,"speed","speed1","speed2",&speed1,&speed2);
      param_getf_2(ppl,"seed","seed1","seed2",&seed1,&seed2);
      param_getf_2(ppl,"ampl","ampl1","ampl2",&ampl1,&ampl2);
      param_getf_2(ppl,"dotdens","dotdens1","dotdens2",&ddens1,&ddens2);

      //param_getf_2(ppl,"cx","cx1","cx2",&cx1,&cx2);

      dualflag = 1; // There are two dot fields
    }else{
      theta1 = theta2 = 0.0;
      speed1 = param_getf_exit(ppl,"speed");
      seed1  = param_geti_exit(ppl,"seed");
      ampl1  = param_getf_exit(ppl,"ampl");
      ddens1 = param_getf_exit(ppl,"ampl");
      dualflag = 0;  // There is only one dot field
    }

    t0     = param_getf_exit(ppl,"st0");
    tn     = param_getf_exit(ppl,"stn");
    cx     = param_getf_dflt(ppl,"cx",0.0);
    cy     = param_getf_dflt(ppl,"cy",0.0);
    //xp     = param_getf_dflt(ppl,"extra_param",0.0);  WYETH REMOVE - test
    theta  = param_getf_exit(ppl,"direction");
    coh    = param_getf_exit(ppl,"coherence");
    bgval  = param_getf_exit(ppl,"bgval");
    maxlum = param_getf_exit(ppl,"maxlum");
    dsize  = param_getf_exit(ppl,"dotsize");
    size   = param_getf_exit(ppl,"size");
    xf     = param_getf_exit(ppl,"exp_factor");
    dtf    = param_geti_exit(ppl,"dt_frames");
    aptype = param_geti_exit(ppl,"aptype");
    minlum = 0.0;

    fps = 1.0/(tsamp*tscale); // Frames/s, to compute dx, dy for motion

    // Work in entire 'xn' by 'yn' area, although wasteful if small aperture
    dpf = (int)((float)(xn*yn)*xf*xf*sscale*sscale * ddens1); // Area X density
    //if (dpf <= 0)
    //exit_error("STIMX_DOTS","Zero dots in field 1");
    szx = (float)xn * xf;
    szy = (float)yn * xf;
    //nframes = (int)(tn/(dwell*tscale*tsamp));
    nframes = (int)(tn/(tscale*tsamp));

    get_coherence_dot_movie(sscale,fps,coh,speed1,theta+theta1,szx,szy,nframes,
			    dpf,dtf,seed1,&dotx,&doty,&dotn);
    if (dualflag == 1){
      dpf = (int)((float)(xn*yn)*xf*xf*sscale*sscale * ddens2);
      //if (dpf <= 0)
      //exit_error("STIMX_DOTS","Zero dots in field 2");
      get_coherence_dot_movie(sscale,fps,coh,speed2,theta+theta2,szx,szy,
			      nframes,dpf,dtf,seed2,&dotx2,&doty2,&dotn2);
    }

    // Dots coords are pixels [0..], covering entire frame plus area from xf
    // WYETH HERE - check for xf < 1.0, and add const to dot coords.
    if (xf < 1.0){
      xoff = (float)xn * (1.0 - xf) / 2.0;
      yoff = (float)yn * (1.0 - xf) / 2.0;
      translate_dot_movie(dotx,doty,dotn,nframes,xoff,yoff);
      if (dualflag == 1)
	translate_dot_movie(dotx2,doty2,dotn2,nframes,xoff,yoff);
    }

    // Apply aperture to dots
    x0 = (float)((int)(xn-1)/2);  // Origin of dot frame (pixels)
    y0 = (float)((int)(yn-1)/2);  // E.g., 0..7 center is 3.0
    if ((aptype > 0) && (aptype != 3)){
      mask_dot_movie(dotx,doty,NULL,dotn,nframes,x0,y0,cx,cy,size,size,
		     sscale,aptype);
      if (dualflag == 1)
	mask_dot_movie(dotx2,doty2,NULL,dotn2,nframes,x0,y0,cx,cy,size,size,
		       sscale,aptype);
    }

    t0i = my_rint(t0 / tscale);  // Time index at which stim starts

  }else if (task == 1){
    //
    //  Get stimulus frame
    //

    if (dotx == NULL)
      exit_error("STIMX_DOTS","Not prepared to generate stim");

    //
    // Compute the index into the dot frame
    //
    k = (ti-t0i)/tsamp;  // Number of pattern frames from stim start

    // Make new flat frame
    tframe = get_2d_farray(xn,yn);
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	tframe[i][j] = bgval;

    if ((k >= 0) && (k < nframes)){
      make_frame_add_dots_0(tframe,xn,yn,dotx[k],doty[k],dotn[k],
			    ampl1,minlum,maxlum,dsize,NULL);

      if (dualflag == 1)
	make_frame_add_dots_0(tframe,xn,yn,dotx2[k],doty2[k],dotn2[k],
			      ampl2,minlum,maxlum,dsize,NULL);
      if (aptype == 3)
	multiply_frame_gaussian_0(tframe,xn,yn,0.0,0.0,size,sscale);
    }
    
    *rframe = tframe;
    
  }else{
    exit_error("STIMX_DOTS","unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              STIMX_DOT_PAIRED                             */
/*                                                                           */
/*****************************************************************************/
void stimx_dot_paired(ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,rframe)
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // Time units per stim frame
     int ti;                       // Time index (tscale units), frame number
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframe;              // Return 2D frame into this address
{
  int i,j,k;
  int dpf,aptype,seed,seed2,lifti;
  float **tframe,x0,y0,fps,szx,szy,lifet;
  float t0,tn,theta,theta2,dist,offo,offx,offy,ddens,xf,cx,cy;
  // Static
  static int t0i,nframes;
  static int *dotn,*dotn2;
  static float ampl1,ampl2,bgval,minlum,maxlum,dsize,size;
  static float **dotx,**doty,**dotx2,**doty2;

  if (task == -1){  // Clean up
    if (dotx != NULL)
      free_2d_farray(dotx,nframes);
    if (doty != NULL)
      free_2d_farray(doty,nframes);
    if (dotn != NULL)
      myfree(dotn);

    if (dotx2 != NULL)
      free_2d_farray(dotx2,nframes);
    if (doty2 != NULL)
      free_2d_farray(doty2,nframes);
    if (dotn2 != NULL)
      myfree(dotn2);
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //
    seed   = param_geti_exit(ppl,"seed");
    seed2  = param_geti_dflt(ppl,"seed_unpair",seed);
    if (seed2 == 0)
      seed2 = seed;
    dist   = param_getf_exit(ppl,"distance");
    offo   = param_getf_exit(ppl,"offset");
    //ampl   = param_getf_exit(ppl,"ampl");
    param_getf_2(ppl,"ampl","ampl1","ampl2",&ampl1,&ampl2);


    t0     = param_getf_exit(ppl,"st0");
    tn     = param_getf_exit(ppl,"stn");
    cx     = param_getf_dflt(ppl,"cx",0.0);
    cy     = param_getf_dflt(ppl,"cy",0.0);
    theta  = param_getf_exit(ppl,"direction");
    theta2 = param_getf_exit(ppl,"direction2");
    lifet  = param_getf_exit(ppl,"lifetime");    // (s)
    bgval  = param_getf_exit(ppl,"bgval");
    maxlum = param_getf_exit(ppl,"maxlum");
    dsize  = param_getf_exit(ppl,"dotsize");
    ddens  = param_getf_exit(ppl,"dotdens");
    size   = param_getf_exit(ppl,"size");
    xf     = param_getf_exit(ppl,"exp_factor");
    aptype = param_geti_exit(ppl,"aptype");
    minlum = 0.0;

    if (xf < 1.0)
      xf = 1.0;

    fps = 1.0/(tsamp*tscale); // Frames/s, to compute dx, dy for motion

    lifti = my_rint(lifet / tscale);

    // Work in entire 'xn' by 'yn' area, although wasteful if small aperture
    dpf = (int)((float)(xn*yn)*xf*xf*sscale*sscale * ddens); // Area X density
    szx = (float)xn * xf;
    szy = (float)yn * xf;
    //nframes = (int)(tn/(dwell*tscale*tsamp));
    nframes = (int)(tn/(tscale*tsamp));

    /*
      old way
    get_paired_dot_movie(sscale,fps,dist,offo,theta,theta2,szx,szy,lifti,
			 nframes,dpf,seed,seed2,&dotx,&doty,&dotn,
			 &dotx2,&doty2,&dotn2);
    */

    offx = offo * cos((theta + 90.0)/180.0 * M_PI);
    offy = offo * sin((theta + 90.0)/180.0 * M_PI);

    get_pairable_dot_movie(sscale,fps,dist,0.0,0.0,theta,szx,szy,
			   lifti,nframes,dpf,seed,&dotx,&doty,&dotn);
    get_pairable_dot_movie(sscale,fps,dist,offx,offy,theta2,szx,szy,
			   lifti,nframes,dpf,seed2,&dotx2,&doty2,&dotn2);


    // Dots coords are pixels [0..], covering entire frame plus area from xf

    // Apply aperture to dots
    x0 = (float)((int)(xn-1)/2);  // Origin of dot frame (pixels)
    y0 = (float)((int)(yn-1)/2);  // E.g., 0..7 center is 3.0
    if ((aptype > 0) && (aptype != 3)){
      mask_dot_movie(dotx,doty,NULL,dotn,nframes,x0,y0,cx,cy,size,size,
		     sscale,aptype);
      mask_dot_movie(dotx2,doty2,NULL,dotn2,nframes,x0,y0,cx,cy,size,size,
		     sscale,aptype);
    }

    t0i = my_rint(t0 / tscale);  // Time index at which stim starts

  }else if (task == 1){
    //
    //  Get stimulus frame
    //

    if (dotx == NULL)
      exit_error("STIMX_DOT_PAIRED","Not prepared to generate stim");

    //
    // Compute the index into the dot frame
    //
    k = (ti-t0i)/tsamp;  // Number of pattern frames from stim start

    // Make new flat frame
    tframe = get_2d_farray(xn,yn);
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	tframe[i][j] = bgval;

    if ((k >= 0) && (k < nframes)){
      make_frame_add_dots_0(tframe,xn,yn,dotx[k],doty[k],dotn[k],
			    ampl1,minlum,maxlum,dsize,NULL);
      make_frame_add_dots_0(tframe,xn,yn,dotx2[k],doty2[k],dotn2[k],
			    ampl2,minlum,maxlum,dsize,NULL);
      if (aptype == 3)
	multiply_frame_gaussian_0(tframe,xn,yn,0.0,0.0,size,sscale);
    }

    *rframe = tframe;

  }else{
    exit_error("STIMX_DOT_PAIRED","unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           STIMX_DOTS_STEREOGRAM                           */
/*                                                                           */
/*  Stimtypes:                                                               */
/*    randot_stereogram                                                      */
/*    rds_disp_mod        - Disparity modulation, eg. Sanada/DeAngelis       */
/*                                                                           */
/*****************************************************************************/
void stimx_dots_stereogram(ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,
			   rframl,rframr)
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // Time units per stim frame
     int ti;                       // Time index (tscale units), frame number
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframl;              // Return 2D frame for left stim.
     float ***rframr;              // Return 2D frame for right stim.
{
  int i,j,k;
  int dpf,dtf,aptype,seed,dotcorr,anticor,ptype,dispc,disps,nplanes;
  float **tframl,**tframr,x0,y0,fps,szx,szy,tf;
  float t0,tn,theta,coh,speed,ddens,xf,cx,cy,*vampl,*vampr;
  float px,py,pw,ph,pdx,pdy,px0_pix,py0_pix,pw_pix,ph_pix,dx_pix,dy_pix;
  float ***dl,***dr,pdx0,pdy0,dx0_pix,dy0_pix,disp0,dispr;

  // NEW for MidMod
  int distrib,dt,ord,seed_dot,seqn,denv;
  float v1p,v1a,*seq,maxshift;

  char *monoflag;
  char *pwf = NULL;
  struct stmh_dotmov *dml,*dmr;     // Left and right eye dot movies

  // Static
  static int t0i,nframes,meanzero,dwell,make_right,make_left;
  static int *dotn_l,*dotn_r,bipol;
  static float crx1,cry1,crx2,cry2;
  static float ampl,bgval,minlum,maxlum,dsize,size;
  static float **dotx_l,**doty_l,**dotx_r,**doty_r,**dampl,**dampr;
  static char *stimtype = NULL;

  if (task == -1){  // Clean up
    if (dotx_l != NULL)  free_2d_farray(dotx_l,nframes);
    if (doty_l != NULL)  free_2d_farray(doty_l,nframes);
    if (dotn_l != NULL)  myfree(dotn_l);
    if (dotx_r != NULL)  free_2d_farray(dotx_r,nframes);
    if (doty_r != NULL)  free_2d_farray(doty_r,nframes);
    if (dotn_r != NULL)  myfree(dotn_r);
    if (stimtype != NULL)  myfree(stimtype);
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //
    stimtype = param_getc_exit(ppl,"stim_type");

    t0      = param_getf_exit(ppl,"st0");
    tn      = param_getf_exit(ppl,"stn");
    cx      = param_getf_dflt(ppl,"cx",0.0);
    cy      = param_getf_dflt(ppl,"cy",0.0);
    theta   = param_getf_dflt(ppl,"direction",0.0);
    coh     = param_getf_dflt(ppl,"coherence",0.0);
    speed   = param_getf_dflt(ppl,"speed",0.0);
    v1p     = param_getf_dflt(ppl,"step_pref",1.0);
    v1a     = param_getf_dflt(ppl,"step_anti",-1.0);
    ampl    = param_getf_exit(ppl,"ampl");
    bgval   = param_getf_exit(ppl,"bgval");
    maxlum  = param_getf_exit(ppl,"maxlum");
    dsize   = param_getf_exit(ppl,"dotsize");
    ddens   = param_getf_exit(ppl,"dotdens");
    size    = param_getf_exit(ppl,"size");
    tf      = param_getf_dflt(ppl,"tf",0.0);
    xf      = param_getf_dflt(ppl,"exp_factor",1.0);
    bipol   = param_geti_dflt(ppl,"bipolar",1);
    dotcorr = param_geti_dflt(ppl,"dotcorr",1); // -1=opp.sign, 0=spatial uncor
    anticor = param_geti_dflt(ppl,"anticorr",1);
    dwell   = param_geti_dflt(ppl,"dwell",1);
    dtf     = param_geti_dflt(ppl,"dt_frames",1);
    aptype  = param_geti_exit(ppl,"aptype");
    seed    = param_geti_exit(ppl,"seed");
    meanzero= param_geti_dflt(ppl,"stim_mean_zero",0);
    monoflag= param_getc_dflt(ppl,"stim_monocular","both"); // "left","right"

    dampl = dampr = NULL;

    if (strcmp(stimtype,"rds_disp_mod")==0){
      //
      //  Change in Disparity
      //
      disp0 = param_getf_exit(ppl,"disp_mid");    // Starting disparity
      dispr = param_getf_exit(ppl,"disp_range");  // Range (deg) +/- this val
      disps = param_getf_exit(ppl,"disp_sign");   // 1-approach, -1 recede
      dispc = param_geti_dflt(ppl,"disp_mod_corr",1); // 0-uncorr. L,R eye
      if (tf == 0.0)
	exit_error("STIMX_DOTS_STEREOGRAM","TF should not be zero");
    }else if ((strcmp(stimtype,"rds_mid_mod")==0) ||
	      (strcmp(stimtype,"rds_cd_mod")==0)){
      //
      //  MidMod - MID random moduation
      //

      distrib  = param_geti_dflt(ppl,"distrib",7);
      dt       = param_geti_exit(ppl,"dt");
      if (distrib == 7){
	ord = param_geti_exit(ppl,"mseq_ord");
      }else
	exit_error("STIMX_DOTS_STEREOGRAM","Must be distrib = 7");
      seed_dot  = param_geti_exit(ppl,"dot_seed");

      maxshift = param_getf_exit(ppl,"max_shift");  // deg
      denv     = param_geti_dflt(ppl,"depth_envelope",0);

      nplanes  = param_geti_dflt(ppl,"n_planes",1);


    }else{
      ptype   = param_geti_dflt(ppl,"patch_type",0);  // 0-no patch
      if (ptype == 0){
	// Full field - there is no patch
	px = py = 0.0;
	pw = ph = -1.0;
      }else if (ptype != 0){
	// Pars for patch to be shifted
	px      = param_getf_exit(ppl,"patch_cx");
	py      = param_getf_exit(ppl,"patch_cy");
	pw      = param_getf_exit(ppl,"patch_w");
	ph      = param_getf_exit(ppl,"patch_h");

	if ((pw >= xn*sscale) || (ph >= yn*sscale))
	  exit_error("STIMX_DOTS_STEREOGRAM",
		     "Patch size exceeds stimulus field, use patch_type 0.");
      }
      pdx     = param_getf_exit(ppl,"patch_dx");
      pdy     = param_getf_exit(ppl,"patch_dy");
      pdx0    = param_getf_dflt(ppl,"patch_dx0",-pdx);
      pdy0    = param_getf_dflt(ppl,"patch_dy0",-pdy);
      pwf     = param_getc_dflt(ppl,"patch_waveform","const"); // "sine","ramp"
    }

    make_left = make_right = 1;
    if (strcmp(monoflag,"left")==0)
      make_right = 0;
    else if (strcmp(monoflag,"right")==0)
      make_left = 0;

    pw_pix  = pw/sscale;
    ph_pix  = ph/sscale;
    px0_pix = (float)xn/2.0 + px/sscale - pw_pix/2.0;
    py0_pix = (float)yn/2.0 + py/sscale - ph_pix/2.0;

    if ((pw < 0.0) || (ph < 0.0)){
      //
      //  There is no patch - the whole image is to be shifted
      //
      crx1 = 0.0;
      crx2 = (float)xn;
      cry1 = 0.0;
      cry2 = (float)yn;
    }else{
      crx1 = px0_pix;  // Clip region for patch that is shifted in depth
      crx2 = px0_pix + pw_pix - 1.0;
      cry1 = py0_pix;
      cry2 = py0_pix + ph_pix - 1.0;
    }

    dx_pix  = pdx/sscale;
    dy_pix  = pdy/sscale;

    dx0_pix  = pdx0/sscale;
    dy0_pix  = pdy0/sscale;

    if (meanzero == 1){
      minlum = -maxlum;
    }else
      minlum = 0.0;

    if (xf < 1.0)
      xf = 1.0;

    fps = 1.0/(dwell*tsamp*tscale); // Frames/s, to compute dx, dy for motion

    // Work in entire 'xn' by 'yn' area, although wasteful if small aperture
    dpf = (int)((float)(xn*yn)*xf*xf*sscale*sscale * ddens); // Area X density
    szx = (float)xn * xf;
    szy = (float)yn * xf;
    nframes = (int)(tn/((float)dwell*tscale*tsamp));

    //printf("nframes = %d\n",nframes);
    //printf("tn = %f\n",tn);
    //printf("dwell = %d\n",dwell);
    //printf("tscale = %f\n",tscale);
    //printf("tsamp = %f\n",(float)tsamp);

    //  Dots have ampl of +1, -1
    if ((strcmp(stimtype,"rds_mid_mod")==0) ||
	(strcmp(stimtype,"rds_cd_mod")==0)){
      //
      //  MidMod - MID random moduation
      //

      // Set float values for random sequences
      fps = 1.0/(tsamp*tscale);  // Frames/s - Pattern frames (not blanks)

      // Return a sequence of random values (no repeats for dt, minimal seq)
      //   Note, for order 11 mseq, 'seqn' will be 2048, regardless of 'tn'
      stm_get_ran_seq(NULL,distrib,ord,seed,fps,tn,dt,v1p,v1a,0.0,&seq,&seqn);

      if (strcmp(stimtype,"rds_mid_mod")==0)
	stm_dotmov_get_deep(sscale,szx,szy,fps,speed,maxshift,denv,
			    nframes,dpf,dt,seed,seq,seqn,bipol,dotcorr,
			    &dml,&dmr);
      else if (strcmp(stimtype,"rds_cd_mod")==0)
	stm_dotmov_get_planes(sscale,szx,szy,fps,speed,maxshift,denv,
			      nframes,dpf,nplanes,dt,seed,seq,seqn,bipol,
			      &dml,&dmr);
      else
	exit_error("STIMX_DOTS_STEREOGRAM","Unexpected stimtype");


      stm_dotmov_get_parts(dml,1,1,&dotx_l,&doty_l,&dampl,&dotn_l);
      stm_dotmov_get_parts(dmr,1,1,&dotx_r,&doty_r,&dampr,&dotn_r);

      // *** WYETH - NEED TO FREE 'dml' and 'dmr' here ***
      // ***   and ultimately free dotx_r etc, when new stim is created.

    }else if (strcmp(stimtype,"rds_disp_mod")==0){
      //
      //  For Sanada IOVD and Combined, but not for Sanada CD stim.
      //

      get_rds_disp_mod(sscale,fps,disp0,dispr,disps,theta,tf,szx,szy,
		       dx_pix,dy_pix,dx0_pix,dy0_pix,nframes,dpf,seed,dispc,
		       &dotx_l,&doty_l,&dotn_l,&dotx_r,&doty_r,&dotn_r,
		       &dampl,&dampr);
    }else{
      //
      //  stimtype "randot_stereogram"
      //
      //printf("__________________coh = %f\n",coh);
      get_dot_stereogram(ptype,sscale,fps,coh,speed,theta,tf,pwf,szx,szy,
			 crx1,crx2,cry1,cry2,dx_pix,dy_pix,dx0_pix,dy0_pix,
			 nframes,ddens,dtf,seed,dotcorr,
			 &dotx_l,&doty_l,&dotn_l,&dotx_r,&doty_r,&dotn_r,
			 &dampl,&dampr);
    }

    if (bipol == 1){
      if ((meanzero == 1) && (bgval != 0.0))
	exit_error("STIMX_DOTS_STEREOGRAM","Expecting bgval to be 0.0");
      else if (bgval == 0.5){
	for(i=0;i<nframes;i++){
	  multiply_farray(dampl[i],dotn_l[i],ampl/2.0);
	  multiply_farray(dampr[i],dotn_r[i],ampl/2.0);
	}
      }
    }

    if (dotcorr == -1){
      for(i=0;i<nframes;i++){
	multiply_farray(dampr[i],dotn_r[i],-1.0);
      }
    }

    // Dots coords are pixels [0..], covering entire frame plus area from xf

    // Apply aperture to dots
    x0 = (float)((int)(xn-1)/2);  // Origin of dot frame (pixels)
    y0 = (float)((int)(yn-1)/2);  // E.g., 0..7 center is 3.0
    if ((aptype > 0) && (aptype != 3)){
      mask_dot_movie(dotx_l,doty_l,dampl,dotn_l,nframes,x0,y0,cx,cy,size,size,
		     sscale,aptype);
      mask_dot_movie(dotx_r,doty_r,dampr,dotn_r,nframes,x0,y0,cx,cy,size,size,
		     sscale,aptype);
    }

    t0i = my_rint(t0 / tscale);  // Time index at which stim starts

    myfree(monoflag);
    if (pwf != NULL)  myfree(pwf);
  }else if (task == 1){
    //
    //  Get stimulus frame
    //

    if (dotx_l == NULL)
      exit_error("STIMX_DOTS_STEREOGRAM","Not prepared to generate stim");

    //
    // Compute the index into the dot frame
    //
    //k = (ti-t0i)/tsamp;
    k = (ti-t0i)/(dwell*tsamp);  // Number of pattern frames from stim start

    // Make new flat frame
    tframl = get_2d_farray(xn,yn);
    tframr = get_2d_farray(xn,yn);
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++){
	tframl[i][j] = bgval;
	tframr[i][j] = bgval;
      }

    if ((k >= 0) && (k < nframes)){

      //if (bipol == 1){
      if (dampl != NULL){
	vampl = dampl[k];
	vampr = dampr[k];
      }else{
	vampl = vampr = NULL;
      }

      if (make_left == 1){

	//
	//  *** WYETH:  the clip region is the 'mask' rectangle size.
	//      But should it not be the stimulus size ???
	//

	/*make_frame_add_dots_0(tframl,xn,yn,dotx_l[k],doty_l[k],dotn_l[k],
	  ampl,minlum,maxlum,dsize,vampl);*/
	make_frame_add_dots_clipr(tframl,xn,yn,dotx_l[k],doty_l[k],dotn_l[k],
				  ampl,minlum,maxlum,dsize,vampl,crx1,crx2,
				  cry1,cry2);


	if (aptype == 3)
	  multiply_frame_gaussian_0(tframl,xn,yn,0.0,0.0,size,sscale);
      }
      if (make_right == 1){
	/*make_frame_add_dots_0(tframr,xn,yn,dotx_r[k],doty_r[k],dotn_r[k],
	  ampl,minlum,maxlum,dsize,vampr);*/
	make_frame_add_dots_clipr(tframr,xn,yn,dotx_r[k],doty_r[k],dotn_r[k],
				  ampl,minlum,maxlum,dsize,vampr,crx1,crx2,
				  cry1,cry2);
	if (aptype == 3)
	  multiply_frame_gaussian_0(tframr,xn,yn,0.0,0.0,size,sscale);
      }
    }

    *rframl = tframl;
    *rframr = tframr;

  }else{
    exit_error("STIMX_DOTS_STEREOGRAM","unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                 STIMX_GP                                  */
/*                                                                           */
/*****************************************************************************/
void stimx_gp(ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,rframe)
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // Time units per stim frame
     int ti;                       // Time index (tscale units), frame number
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframe;              // Return 2D frame into this address
{
  int i,j,k;
  int seed;
  float **tframe,t0,tn,ddens,size,cx,cy;

  // Static
  static int t0i,nframes,aptype,dn;
  static int *seedlist;
  static float r,theta,amp1,amp2,bgval,ampmin,ampmax,dsize,meandot,xf,xc,yc,sz;

  if (task == -1){  // Clean up
    if (seedlist != NULL)
      myfree(seedlist);
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //

    t0     = param_getf_exit(ppl,"st0");
    tn     = param_getf_exit(ppl,"stn");
    cx     = param_getf_dflt(ppl,"cx",0.0);
    cy     = param_getf_dflt(ppl,"cy",0.0);
    r      = param_getf_exit(ppl,"gp_r");
    theta  = param_getf_exit(ppl,"gp_theta");
    amp1   = param_getf_exit(ppl,"dot1_amp");
    amp2   = param_getf_exit(ppl,"dot2_amp");
    bgval  = param_getf_exit(ppl,"bgval");
    ampmin = param_getf_exit(ppl,"ampmin");
    ampmax = param_getf_exit(ppl,"ampmax");
    dsize  = param_getf_exit(ppl,"dotsize");
    ddens  = param_getf_exit(ppl,"dotdens");
    xf     = param_getf_exit(ppl,"exp_factor");
    size   = param_getf_exit(ppl,"size");
    seed   = param_geti_exit(ppl,"seed");
    aptype = param_geti_exit(ppl,"aptype");

    xc = (float)(xn-1)/2.0 + cx/sscale;  // Center w.r.t. 0..xn-1 (pix)
    yc = (float)(yn-1)/2.0 + cy/sscale;  // Center w.r.t. 0..yn-1 (pix)

    meandot = (float)(xn*yn)*xf*xf*sscale*sscale * ddens; // Area times dens

    sz = size/sscale;

    nframes = (int)(tn/(tscale*tsamp));
    seedlist = get_seeds(seed,100000,nframes);

    t0i = my_rint(t0 / tscale);  // Time index at which stim starts

  }else if (task == 1){
    //
    //  Get stimulus frame
    //

    if (seedlist == NULL)
      exit_error("STIMX_GP","Not prepared to generate stim");

    //
    // Compute the index into the dot frame
    //
    k = (ti-t0i)/tsamp;  // Number of pattern frames from stim start

    // Make new flat frame
    tframe = get_2d_farray(xn,yn);
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	tframe[i][j] = bgval;

    if ((k >= 0) && (k < nframes)){

      dn = myrand_poisson_count(meandot,seedlist[k]);

      make_frame_gp_trans_0(tframe,xn,yn,sscale,r,theta,amp1,amp2,
			    bgval,ampmin,ampmax,dsize,xf,dn,seedlist[k]);
      apply_aperture_scale(tframe,xn,yn,aptype,xc,yc,sz,sz,0.0,bgval,1.0);

      //if (aptype == 3)  // FUTURE?  UNIFY APERTURES?
      //multiply_frame_gaussian_0(tframe,xn,yn,0.0,0.0,size,sscale);
    }

    *rframe = tframe;

  }else{
    exit_error("STIMX_GP","unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              STIMX_IMAGE_SET                              */
/*                                                                           */
/*  *** NOTE:  This routine re-uses the 'rframe' that is returned, thus the  */
/*             called MUST COPY data from 'rframe' and SHOULD NOT free it.   */
/*                                                                           */
/*****************************************************************************/
void stimx_image_set(ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,rframe)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;                    // Time units per stim frame (e.g., 1)
     int ti;                       // Time index (tscale units)
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframe;              // Return 2D frame into this address
{
  int k;
  int seed;
  float bgval,st0,stn,cx,cy,**gdata;

  // Static
  static int t0i,nframes,x0,y0,gxn,gyn,boxw,boxh,dwell;
  static int *seedlist = NULL;
  static float bgamp;
  static char *image_dir = NULL;
  static float **fprev = NULL;
  static int fprevk = -1;
    

  if (task == -1){  // Clean up
    if (seedlist != NULL)
      myfree(seedlist);
    if (image_dir != NULL)
      myfree(image_dir);
    if (fprev  != NULL){
      free_2d_farray(fprev,xn);
      fprev = NULL;
    }
    fprevk = -1;
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //

    st0     = param_getf_exit(ppl,"st0");
    stn     = param_getf_exit(ppl,"stn");
    cx      = param_getf_dflt(ppl,"cx",0.0);
    cy      = param_getf_dflt(ppl,"cy",0.0);
    gxn     = param_geti_exit(ppl,"grid_w");
    gyn     = param_geti_exit(ppl,"grid_h");
    boxw    = param_geti_exit(ppl,"box_w");
    boxh    = param_geti_exit(ppl,"box_h");
    dwell   = param_geti_exit(ppl,"dwell");
    seed    = param_geti_exit(ppl,"seed");
    bgval   = param_getf_exit(ppl,"bgval");
    image_dir = param_getc_exit(ppl,"set_dir");

    // WYETH - read in 'setname' param, which is "vanhateren"


    bgamp = bgval * 2.0 - 1.0; // Convert 0..1 --> -1..1

    nframes = stm_get_n(zn,tscale,st0,stn,dwell,tsamp);  // Unique frame count
    seedlist = get_seeds(seed,100000,nframes);  // One seed per uniqe frame

    t0i = stm_get_ti0(zn,tscale,st0,tsamp);  // Time index at which stim starts

    x0 = (int)((float)(xn - gxn*boxw)/2.0);  // Absolute center
    y0 = (int)((float)(yn - gyn*boxh)/2.0);

    x0 += my_rint(cx/sscale);  // User x-offset
    y0 += my_rint(cy/sscale);  // User y-offset

    fprev = NULL;  // If params have changed, we should recompute image
    fprevk = -1;

  }else if (task == 1){
    //
    //  Return stimulus frame
    //

    if (seedlist == NULL)
      exit_error("STIMX_IMAGE_SET","Not prepared to generate stim");

    k =  (int)((ti-t0i)/(dwell*tsamp)); // index in random seq. for this frame
    if ((k < 0) || (k >= nframes)){
      *rframe = NULL;  // Tell the caller the frame is background, not stim
      return;
    }

    if (k != fprevk){ // Must make new frame

      //  Fill 2D array 'gdata' with random values (will free these)
      gdata = stm_vanhat_get_patch(seedlist[k],gxn,gyn,image_dir,"-1 to 1");

      // Create and keep a pointer to this stimulus frame (caller must copy)
      if (fprev != NULL)
	free_2d_farray(fprev,xn);
      fprev = stm_noise_get_grid(xn,yn,gxn,gyn,boxw,boxh,x0,y0,gdata,bgamp);
      fprevk = k;

      free_2d_farray(gdata,gxn);
    }

    *rframe = fprev;

    //  Apply normalization here?

  }else{
    exit_error("STIMX_IMAGE_SET","unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                STIMX_NOISE                                */
/*                                                                           */
/*****************************************************************************/
void stimx_noise(ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,rframe)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;                    // Time units per stim frame (e.g., 1)
     int ti;                       // Time index (tscale units)
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframe;              // Return 2D frame into this address
{
  int i,j,k;
  int seed,nrpt;
  float bgval,st0,stn,**pow_data,ph;
  float **tframe,**tframe1,**tframe2;

  // Static
  static int t0i,nframes,x0,y0,gxn,gyn,boxw,boxh,dwell,ncol,binflag;
  static int *seedlist = NULL;
  static float bgamp,fpow;
  static float sf,tf,theta,phase,cx,cy;  // For 2nd order motion
  static float **gdata = NULL;
  static char *rgen = NULL;
  static char *motion;

  if (task == -1){  // Clean up
    if (seedlist != NULL)
      myfree(seedlist);
    if (gdata != NULL)
      free_2d_farray(gdata,gxn);
    if (rgen != NULL)
      myfree(rgen);
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //
    st0     = param_getf_exit(ppl,"st0");
    stn     = param_getf_exit(ppl,"stn");
    cx      = param_getf_dflt(ppl,"cx",0.0);
    cy      = param_getf_dflt(ppl,"cy",0.0);
    gxn     = param_geti_exit(ppl,"grid_w");
    gyn     = param_geti_exit(ppl,"grid_h");
    boxw    = param_geti_exit(ppl,"box_w");
    boxh    = param_geti_exit(ppl,"box_h");
    dwell   = param_geti_exit(ppl,"dwell");  // <= 0 for static noise
    ncol    = param_geti_dflt(ppl,"n_color",2);
    seed    = param_geti_exit(ppl,"seed");
    rgen    = param_getc_dflt(ppl,"randgen","default");
    bgval   = param_getf_exit(ppl,"bgval");
    fpow    = param_getf_dflt(ppl,"power_law_2d",0.0);

    motion  = param_getc_dflt(ppl,"motion","none");
    if (strcmp(motion,"sine_mult")==0){
      // A drifting sinusoid will multiply the noise
      
      sf      = param_getf_exit(ppl,"sf");
      tf      = param_getf_exit(ppl,"tf");
      phase   = param_getf_exit(ppl,"phase");
      theta   = param_getf_exit(ppl,"direction");
      binflag = param_geti_dflt(ppl,"binarize",0);
    }

    bgamp = bgval * 2.0 - 1.0; // Convert 0..1 --> -1..1

    // *** WYETH OLD WAY.  Use 'stm_get_n' INSTEAD
    //nframes = (int)(stn/(dwell*tscale*tsamp)); // Unique frame count
    if (dwell > 0)
      nframes = stm_get_n(zn,tscale,st0,stn,dwell,tsamp);  // Unique frame count
    else
      nframes = 1; // Static noise

    if (strcmp(rgen,"default")==0){
      seedlist = get_seeds(seed,100000,nframes);  // One seed per uniqe frame
    }else if (strcmp(rgen,"ejrand")==0){
      //myrand_ejrand_test(1000000,1);  // Test the random distribution
      seedlist = ejrand_get_seeds(seed,gxn*gyn,nframes);
      // *** USE THIS TO PRINT OUT SEED LIST FOR CHECKIN
      /***
      for(i=0;i<nframes;i++)
	printf("seed %d  %d\n",i,seedlist[i]);
      exit(0);
      ***/
    }else if (strcmp(rgen,"ejrand_genseq")==0){
      //
      //  This does not make a stimulus, it generates seed values to then
      //  be placed as a VAR parameter for the .stm file used for this 
      //  generation, and then re-run the simulation.
      //
      //  The sequence generated replicates the contiguous generation that
      //  is being done in the Horwitz lab, for Patrick's experiments.
      //
      
      nrpt = param_geti_exit(ppl,"stim_nrpt");
      printf("    Generating a seed list for %d repeats.\n",nrpt);
      printf("      %d x %d  grid\n",gxn,gyn);
      printf("      %d frames per stimulus\n",nframes);

      seedlist = ejrand_get_seeds(seed,gxn*gyn*nframes,nrpt);
      for(i=0;i<nrpt;i++)
	printf(" %d",seedlist[i]);
      printf("\n");

      exit_error("STIMX_NOISE","Generated sequence to use as VAR param");
    }else{
      printf("  *** randgen = %s\n",rgen);
      exit_error("STIMX_NOISE","Unknown 'randgen' value");
    }

    gdata = get_2d_farray(gxn,gyn);  // Storage for random pattern

    t0i = stm_get_ti0(zn,tscale,st0,tsamp);

    x0 = (int)((float)(xn - gxn*boxw)/2.0);  // Absolute center
    y0 = (int)((float)(yn - gyn*boxh)/2.0);

    x0 += my_rint(cx/sscale);  // User x-offset
    y0 += my_rint(cy/sscale);  // User y-offset

  }else if (task == 1){
    //
    //  Return stimulus frame
    //

    if (seedlist == NULL)
      exit_error("STIMX_NOISE","Not prepared to generate stim");

    if (dwell > 0)
      k =  (int)((ti-t0i)/(dwell*tsamp)); // index in random seq. for this frame
    else
      k = 0;  // Static noise

    if ((k < 0) || (k >= nframes)){
      *rframe = NULL;  // Tell the caller the frame is background, not stim
      return;
    }

    //  Fill 2D array 'gdata' with random values.
    if (ncol == 2)
      gdata = myrand_get_std_bin_float_2d(gxn,gyn,seedlist[k]);
    else if (ncol == 4){
      if (strcmp(rgen,"ejrand")==0){
	gdata = myrand_get_std_quad_ejrand_2d(gxn,gyn,seedlist[k]);
      }else{
	gdata = myrand_get_std_quad_float_2d(gxn,gyn,seedlist[k]);
      }
    }else
      exit_error("STIMX_NOISE","Number of colors, n_color, must be 2 or 4");


    //  Apply power law in frequency domain
    if (param_test(ppl,"power_law_2d"))
      power_law_2d_transform(gdata,gxn,gyn,fpow,"1 max 0.5 mid");

    tframe1 = stm_noise_get_grid(xn,yn,gxn,gyn,boxw,boxh,x0,y0,gdata,bgamp);
    //*rframe = stm_noise_get_grid(xn,yn,gxn,gyn,boxw,boxh,x0,y0,gdata,bgamp);

    if (strcmp(motion,"sine_mult")==0){  // Second-order motion
      // second order (second-order, non-fourier) motion  WYSO

      ph = phase + (float)ti*tscale*tf * 360.0;

      // Make sinusoid to multiply the noise
      tframe2 = get_2d_farray(xn,yn);
      make_frame_sinusoid_0(tframe2,xn,yn,sf,theta,sscale,ph,cx,cy,1.0);

      //make_frame_binarize_0(tframe2,xn,yn,0.0,1.0,0.0);
      if (binflag == 1)
	make_frame_binarize_0(tframe2,xn,yn,0.0,1.0,0.0);
      else if (binflag == 2)
	make_frame_binarize_0(tframe2,xn,yn,0.0,1.0,-1.0);


      *rframe = multiply_2d_farrays(tframe1,tframe2,xn,yn);

      free_2d_farray(tframe1,xn);
    }else
      *rframe = tframe1;
    
  }else{
    exit_error("STIMX_NOISE","unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            STIMX_PASUPATHY_SHAPE                          */
/*                                                                           */
/*****************************************************************************/
void stimx_pasupathy_shape(ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,rframe)
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // frames per stim pattern frame
     int ti;                       // Time index (tscale units), frame number
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframe;              // Return 2D frame into this address
{
  int i,j,k;
  int occ_seed,occ_gn;
  float **tframe,t0,tn,size,cx,cy,occ_diam,occ_gw,occ_jit,x0,y0;

  // Static
  static int t0i,nframes,shid,fill,normfl,occ_n;
  static float r,theta,fgval,bgval,xc,yc,sz,aa_sd;
  static char *shpset = NULL;
  static float spot_amp,spot_sd,spot_x,spot_y;
  static float *occ_x,*occ_y,occ_r,occ_amp;

  if (task == -1){  // Clean up
    if (shpset != NULL)  myfree(shpset);
    if (occ_x  != NULL)   myfree(occ_x);
    if (occ_y  != NULL)   myfree(occ_y);
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //

    t0     = param_getf_exit(ppl,"st0");
    tn     = param_getf_exit(ppl,"stn");
    cx     = param_getf_dflt(ppl,"cx",0.0);
    cy     = param_getf_dflt(ppl,"cy",0.0);
    shid   = param_geti_exit(ppl,"shape_id");
    theta  = param_getf_exit(ppl,"rotation");
    size   = param_getf_exit(ppl,"size");
    fgval  = param_getf_exit(ppl,"fgval");
    bgval  = param_getf_exit(ppl,"bgval");
    aa_sd  = param_getf_exit(ppl,"antialias");
    normfl = param_geti_dflt(ppl,"norm_flag",0);
    fill   = param_geti_dflt(ppl,"fill",1);
    shpset = param_getc_dflt(ppl,"shape_set","pc2001"); // "pc2001","yasm"

    spot_amp = param_getf_dflt(ppl,"spot_amp",0.0);
    if (spot_amp != 0.0){
      spot_sd = param_getf_exit(ppl,"spot_sd");
      spot_x  = param_getf_exit(ppl,"spot_x");
      spot_y  = param_getf_exit(ppl,"spot_y");
      printf("  spot_amp = %f\n",spot_amp);
    }

    //
    //  Optional grid of occluders
    //
    occ_n = param_geti_dflt(ppl,"occluder_count",0);   // 36
    if (occ_n > 0){
      occ_gn   = param_geti_exit(ppl,"occluder_grid_n");     // 9
      occ_gw   = param_getf_exit(ppl,"occluder_grid_width"); // (deg)
      occ_amp  = param_getf_exit(ppl,"occluder_ampl");       // 0 .. 1
      occ_diam = param_getf_exit(ppl,"occluder_diameter");   // (deg)
      occ_jit  = param_getf_exit(ppl,"occluder_jitter");     // (deg)
      occ_seed = param_geti_exit(ppl,"occluder_seed");       // points, noise

      // Compute lower left of grid in (deg) on screen
      //   ( Screen w - Grid w)
      x0 = (xn*sscale - occ_gw) / 2.0;   // (deg)
      y0 = (yn*sscale - occ_gw) / 2.0;   // (deg)

      stimu_rangrid(sscale,x0,y0,occ_gw,occ_gn,occ_n,occ_jit,occ_seed,
		    &occ_x,&occ_y);

      occ_r = 0.5 * occ_diam / sscale;  // Occluder radius (pix)
    }else{
      occ_x = occ_y = NULL;
    }


    xc = (float)(xn-1)/2.0 + cx/sscale;  // Center w.r.t. 0..xn-1 (pix)
    yc = (float)(yn-1)/2.0 + cy/sscale;  // Center w.r.t. 0..yn-1 (pix)

    nframes = (int)(tn/(tscale*tsamp));

    sz = size/sscale;

    t0i = my_rint(t0 / tscale);  // Time index at which stim starts

  }else if (task == 1){
    //
    //  Get stimulus frame
    //

    //
    // Compute the index into the frame
    //
    k = (ti-t0i)/tsamp;  // Number of pattern frames from stim start

    // Make new flat frame
    tframe = get_2d_farray(xn,yn);
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	tframe[i][j] = bgval;


    //
    //  If there is a spot, draw it
    //
    if (spot_amp != 0.0){
      make_frame_add_2d_gauss(tframe,xn,yn,spot_amp,spot_sd,spot_x,spot_y,
			      sscale);
    }


    //printf("wyeth stim_util.c  HERE 0.. - PASUUUUUUUU k = %d\n",k);
    //printf("HERE 0.. - PASUUUUUUUU k = %d\n",k);
    //printf("HERE 0.. - PASUUUUUUUU k = %d\n",k);

    if ((k >= 0) && (k < nframes)){
      //  WYETH - Inefficient to recompute shape if nothing is changing
      if (1){
	//if (strcmp(shpset,"pc2001")==0){
	stm_pasupathy_shape(NULL,tframe,xn,yn,xc,yc,shpset,shid,theta,sz,
			    fgval,bgval,aa_sd,fill,normfl);
      }else if (strcmp(shpset,"yasm")==0){
	;
	// *** WYETH HERE *** April 11, 2014 - DID I NOT FINISH/NEED THIS???
	// *** WYETH HERE *** April 11, 2014 - DID I NOT FINISH/NEED THIS???
	// *** WYETH HERE *** April 11, 2014 - DID I NOT FINISH/NEED THIS???
	// *** WYETH HERE *** April 11, 2014 - DID I NOT FINISH/NEED THIS???
	/*
	stm_pasupathy_shape_yasm(NULL,tframe,xn,yn,xc,yc,shid,theta,sz,fgval,
				 bgval,aa_sd,fill,normfl);
	*/
      }else{
	exit_error("STIMX_PASUPATHY_SHAPE","unknown 'shape_set'");
      }

      if (occ_n > 0){  // Draw occluders
	stimu_frame_disk_array(tframe,xn,yn,occ_x,occ_y,occ_r,occ_n,occ_amp);
      }
    }

    *rframe = tframe;

  }else{
    exit_error("STIMX_PASUPATHY_SHAPE","unknown task");
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                              STIMX_BLUR_EDGE                              */
/*                                                                           */
/*****************************************************************************/
void stimx_blur_edge(ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,rframe)
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // frames per stim pattern frame
     int ti;                       // Time index (tscale units), frame number
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     float ***rframe;              // Return 2D frame into this address
{
  int i,j,k;
  int seed,pxsz1,pxsz2;
  float **tframe,t0,tn,size,cx,cy,x0,y0;
  float sdpix,sd_lo,sd_hi,phase,sf;
  float **d1,**d2,tmin,tmax,tx,ty;

  // Static
  static int t0i,nframes,aptype;
  static float ori,bgval,xc,yc,sz;
  static float **dlo = NULL;
  static float **dhi = NULL;
  static float **dsin = NULL;

  if (task == -1){  // Clean up
    if (dlo  != NULL)  free_2d_farray(dlo,xn);
    if (dhi  != NULL)  free_2d_farray(dhi,xn);
    if (dsin != NULL)  free_2d_farray(dsin,xn);
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //

    t0     = param_getf_exit(ppl,"st0");
    tn     = param_getf_exit(ppl,"stn");
    cx     = param_getf_dflt(ppl,"cx",0.0);
    cy     = param_getf_dflt(ppl,"cy",0.0);
    ori    = param_getf_exit(ppl,"orientation");
    sf     = param_getf_exit(ppl,"sf");
    phase  = param_getf_exit(ppl,"phase");
    seed   = param_geti_exit(ppl,"seed");
    sd_lo  = param_getf_exit(ppl,"sd_1");
    sd_hi  = param_getf_exit(ppl,"sd_2");
    aptype = param_geti_exit(ppl,"aptype");
    pxsz1  = param_geti_dflt(ppl,"pix_size_1",1);
    pxsz2  = param_geti_dflt(ppl,"pix_size_2",1);
    size   = param_getf_exit(ppl,"size");
    bgval  = param_getf_exit(ppl,"bgval");

    //
    //  Create low and high SF images
    //
    d1 = get_2d_farray(xn,yn);
    d2 = get_2d_farray(xn,yn);

    make_frame_binary_noise_0(d1,xn,yn,pxsz1,0.0,1.0,seed);
    make_frame_binary_noise_0(d2,xn,yn,pxsz2,0.0,1.0,-1);

    /*
    genrand_init((unsigned long)seed);
    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){
	d1[i][j] = genrand_real1();
	d2[i][j] = genrand_real1();
      }
    }
    */
    sdpix = sd_lo / sscale;  // SD in pix
    dlo = smooth_2d_with_gaussian_circular(d1,xn,yn,sdpix,0.01);
    get_min_max_2d_farray(dlo,xn,yn,&tmin,&tmax);
    norm_01_2d_farray(dlo,xn,yn,0);

    sdpix = sd_hi / sscale;  // SD in pix
    dhi = smooth_2d_with_gaussian_circular(d2,xn,yn,sdpix,0.01);
    get_min_max_2d_farray(dhi,xn,yn,&tmin,&tmax);
    norm_01_2d_farray(dhi,xn,yn,0);

    //
    //  Get a sinusoidal array to use for mixing two images
    //
    dsin = get_2d_farray(xn,yn);
    tx = (float)xn/2.0 * sscale + cx;
    ty = (float)yn/2.0 * sscale + cy;
    make_frame_sinusoid_0(dsin,xn,yn,sf,ori,sscale,phase+90.0,tx,ty,1.0);
    
    xc = (float)(xn-1)/2.0 + cx/sscale;  // Center w.r.t. 0..xn-1 (pix)
    yc = (float)(yn-1)/2.0 + cy/sscale;  // Center w.r.t. 0..yn-1 (pix)

    nframes = (int)(tn/(tscale*tsamp));

    sz = size/sscale; // Pixels

    t0i = my_rint(t0 / tscale);  // Time index at which stim starts

  }else if (task == 1){
    //
    //  Get stimulus frame
    //

    //
    // Compute the index into the frame
    //
    k = (ti-t0i)/tsamp;  // Number of pattern frames from stim start

    // Make new flat frame
    tframe = get_2d_farray(xn,yn);
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	tframe[i][j] = bgval;

    if ((k >= 0) && (k < nframes)){
      for(i=0;i<xn;i++){
	for(j=0;j<yn;j++){
	  //tframe[i][j] = dsin[i][j];

	  // **** WYETH - some graded average of the 2 stimuli would help
	  //  ***         when the SF was relatively high, otherwise bands
	  //              will have variable thickness, jaggy edges (oblique)
	  if (dsin[i][j] < 0.0)
	    tframe[i][j] = dlo[i][j];
	  else
	    tframe[i][j] = dhi[i][j];
	}
      }
      apply_aperture_scale(tframe,xn,yn,aptype,xc,yc,sz,0.0,0.0,bgval,1.0);
    }

    *rframe = tframe;

  }else{
    exit_error("STIMX_BLUR_EDGE","unknown task");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              STIMX_RGB_SINE                               */
/*                                                                           */
/*****************************************************************************/
void stimx_rgb_sine(ppl,xn,yn,zn,sscale,tscale,tsamp,ti,task,rframe)
     struct param_pair_list *ppl;  // Stimulus parameter list
     int xn,yn,zn;                 // intended size of stimulus data
     float sscale,tscale;          // (deg/pix) and (s/frame)
     int tsamp;                    // frames per stim pattern frame
     int ti;                       // Time index (tscale units), frame number
     int task;                     // -1-cleanup 0-prep, 1-return stim at 'ti'
     int ***rframe;                // Return 2D frame into this address
{
  int i,j,k;
  int **tframe,tc;
  float t0,tn,tr,tg,tb,ph0;

  // Static
  static int t0i,nframes;
  static float tf,sv,r_off,r_amp,r_pha,g_off,g_amp,g_pha,b_off,b_amp,b_pha;


  if (task == -1){  // Clean up
    ;
  }else if (task == 0){
    //
    //  Prepare for stimulus creation
    //

    t0     = param_getf_exit(ppl,"st0");
    tn     = param_getf_exit(ppl,"stn");
    //cx     = param_getf_dflt(ppl,"cx",0.0);
    //cy     = param_getf_dflt(ppl,"cy",0.0);
    tf     = param_getf_exit(ppl,"tf");

    sv     = param_getf_exit(ppl,"c_scale");

    r_off  = param_getf_exit(ppl,"r_offset");
    r_amp  = param_getf_exit(ppl,"r_ampl");
    r_pha  = param_getf_exit(ppl,"r_phase");

    g_off  = param_getf_exit(ppl,"g_offset");
    g_amp  = param_getf_exit(ppl,"g_ampl");
    g_pha  = param_getf_exit(ppl,"g_phase");

    b_off  = param_getf_exit(ppl,"b_offset");
    b_amp  = param_getf_exit(ppl,"b_ampl");
    b_pha  = param_getf_exit(ppl,"b_phase");

    nframes = (int)(tn/(tscale*tsamp));

    t0i = my_rint(t0 / tscale);  // Time index at which stim starts

  }else if (task == 1){
    //
    //  Get stimulus frame
    //

    ph0 = (float)(ti-t0i)*tscale*tf * 360.0;
    tr = r_off + r_amp*sin(M_PI/180.0 * (r_pha + ph0));
    tg = g_off + g_amp*sin(M_PI/180.0 * (g_pha + ph0));
    tb = b_off + b_amp*sin(M_PI/180.0 * (b_pha + ph0));

    //printf("ti=%d  ph0 =  %f  tb = %f\n",ti,ph0,tb);

    tc = data_util_t11_int_for_rgb(sv*tr,sv*tg,sv*tb);

    // Make new flat frame
    tframe = get_2d_iarray(xn,yn);
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	tframe[i][j] = tc;

    *rframe = tframe;

  }else{
    exit_error("STIMX_RGB_SINE","unknown task");
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                         STIM_UTIL_3C_FROM_TERN_3D                         */
/*                                                                           */
/*  Convert a stimulus of form "3d" into a "3c" by replacing 0.0 and 1.0     */
/*  with the specified RGB colors.                                           */
/*                                                                           */
/*****************************************************************************/
int ***stim_util_3c_from_tern_3d(data,x0,xn,y0,yn,z0,zn,r0,g0,b0,r1,g1,b1,
				 r2,g2,b2)
     float ***data;           // Float data [xn][yn][zn]
     int x0,xn,y0,yn,z0,zn;   // start and duration
     float r0,g0,b0;          // Color to replace 0.0
     float r1,g1,b1;          // Color to replace 1.0
     float r2,g2,b2;          // Color to replace 0.5, Background
{
  int i,j,k;
  int ***dc,c0,c1,c2,*tc;
  float *t;

  c0 = data_util_t11_int_for_rgb(r0,g0,b0);
  c1 = data_util_t11_int_for_rgb(r1,g1,b1);
  c2 = data_util_t11_int_for_rgb(r2,g2,b2);

  dc = get_3d_iarray(xn,yn,zn);

  for(i=x0;i<(x0+xn);i++){
    for(j=y0;j<(y0+yn);j++){
      t = data[i][j];

      tc = dc[i-x0][j-y0];
      for(k=z0;k<(z0+zn);k++){

	if (t[k] < 0.001)
	  tc[k-z0] = c0;
	else if (t[k] > 0.999)
	  tc[k-z0] = c1;
	else
	  tc[k-z0] = c2;
      }
    }
  }
  
  return dc;
}
/**************************************-**************************************/
/*                                                                           */
/*                         STIM_UTIL_3C_2COL_FROM_3D                         */
/*                                                                           */
/*  Convert a stimulus of form "3d" into a "3c" by replacing 0-1 values      */
/*  with colors in the range from RGB_0 to RGB_1.                            */
/*                                                                           */
/*****************************************************************************/
int ***stim_util_3c_2col_from_3d(data,x0,xn,y0,yn,z0,zn,r0,g0,b0,r1,g1,b1)
     float ***data;           // Float data [xn][yn][zn]
     int x0,xn,y0,yn,z0,zn;   // start and duration
     float r0,g0,b0;          // Color 0
     float r1,g1,b1;          // Color 1
{
  int i,j,k;
  int ***dc,*tc;
  float *t,ro,go,bo,ra,ga,ba,tr,tg,tb,ts;

  ro = r0;  ra = r1 - r0;
  go = g0;  ga = g1 - g0;
  bo = b0;  ba = b1 - b0;

  dc = get_3d_iarray(xn,yn,zn);
  for(i=x0;i<(x0+xn);i++){
    for(j=y0;j<(y0+yn);j++){

      t = data[i][j];
      tc = dc[i-x0][j-y0];
      for(k=z0;k<(z0+zn);k++){

	ts = t[k];

	//  2011 Oct 8
	//  *** WYETH THIS IS A HACK, to keep values near mean gray
	//      at mean gray.  This was needed because adding the '..._blur'
	//      function caused small noise around 0.5000, and this ended up
	//      being magnified in the conversion to the 0..255 range in the
	//      t11 format.
	//
	//  *** WARNING: THIS will not help if the mean level is not 0.5
	//  *** WARNING: THIS will not help if the mean level is not 0.5
	//  *** WARNING: THIS will not help if the mean level is not 0.5
	//
	if ((ts > 0.4999) && (ts < 0.5001))
	  ts = 0.5;

	tr = ro + ts * ra;  // red_offset + ts * red_ampl
	tg = go + ts * ga;  // green ...
	tb = bo + ts * ba;  // blue ...

	tc[k-z0] = data_util_t11_int_for_rgb(tr,tg,tb);
      }
    }
  }

  return dc;
}
/**************************************-**************************************/
/*                                                                           */
/*                         STIM_UTIL_3C_2COL_FROM_3D                         */
/*                                                                           */
/*  Convert a stimulus of form "3d" into a "3c" by replacing 0-1 values      */
/*  with colors in the range from RGB_0 to RGB_1.                            */
/*                                                                           */
/*****************************************************************************/
int ***stim_util_3c_4col_from_3d(data,x0,xn,y0,yn,z0,zn,r0,g0,b0,r1,g1,b1,
				 r2,g2,b2,r3,g3,b3,r4,g4,b4)
     float ***data;           // Float data [xn][yn][zn]
     int x0,xn,y0,yn,z0,zn;   // start and duration
     float r0,g0,b0;          // Color 0
     float r1,g1,b1;          // Color 1
     float r2,g2,b2;          // Color 2
     float r3,g3,b3;          // Color 3
     float r4,g4,b4;          // Color 4 - Background color
{
  int i,j,k;
  int ***dc,*tc;
  float *t,tr,tg,tb,ts;

  dc = get_3d_iarray(xn,yn,zn);
  for(i=x0;i<(x0+xn);i++){
    for(j=y0;j<(y0+yn);j++){

      t = data[i][j];
      tc = dc[i-x0][j-y0];
      for(k=z0;k<(z0+zn);k++){

	ts = t[k];

	//  2011 Oct 8
	//  *** WYETH THIS IS A HACK, to keep values near mean gray
	//      at mean gray.  This was needed because adding the '..._blur'
	//      function caused small noise around 0.5000, and this ended up
	//      being magnified in the conversion to the 0..255 range in the
	//      t11 format.
	//
	//  *** WARNING: THIS will not help if the mean level is not 0.5
	//  *** WARNING: THIS will not help if the mean level is not 0.5
	//  *** WARNING: THIS will not help if the mean level is not 0.5
	//
	if ((ts > 0.4999) && (ts < 0.5001)){
	  //tr = tg = tb = 0.5;
	  tr = r4; tg = g4; tb = b4;
	}else if (ts < 0.15){
	  tr = r0; tg = g0; tb = b0;
	}else if (ts < 0.45){
	  tr = r1; tg = g1; tb = b1;
	}else if (ts < 0.85){
	  tr = r2; tg = g2; tb = b2;
	}else{
	  tr = r3; tg = g3; tb = b3;
	}

	tc[k-z0] = data_util_t11_int_for_rgb(tr,tg,tb);
      }
    }
  }

  return dc;
}
/**************************************-**************************************/
/*                                                                           */
/*                               STIM_UTIL_BLUR                              */
/*                                                                           */
/*  Convolve the stimulus in space and/or time.                              */
/*                                                                           */
/*****************************************************************************/
float ***stim_util_blur(mylogf,ppl,sdata,xn,yn,tn,sscale,tscale)
     char *mylogf;
     struct param_pair_list *ppl;
     float ***sdata;
     int xn,yn,tn;
     float sscale,tscale;
{
  int i,j,k;
  float sig,ts,tss,xc,yc;
  float ***d,***smf,**smf_s,***r,**tdata,**tsm;
  char ggstr[SLEN];

  //
  //  Spatial Blurring:   use 3D conv for blurring here
  //
  sig = paramfile_get_float_param_default(ppl,"stim_smooth",0.0);
  sig /= sscale;  // Convert from degr to pixels

  // 'stim_smooth_t' is time of peak of maxwell function, in sec
  ts  = paramfile_get_float_param_default(ppl,"stim_smooth_t",0.0);
  ts /= tscale;  // Convert from sec to time units

  if ((sig > 0.0) && (tn == 1)){
    // May 31, 2019 - to smooth the Pasupathy shapes for blur, for AlexNet

    //  Extract the single 2D frame from the 3D stimulus
    tdata = get_2d_farray(xn,yn);
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	tdata[i][j] = sdata[i+1][j+1][1];

    //printf(" AREA before: %f\n",sum_2d_farray(tdata,0,xn,0,yn));

    //printf(" sign = %f\n",sig);
    tsm = smooth_2d_with_gaussian(tdata,xn,yn,sig,0.01);

    sprintf(ggstr,"      Spatial smoothing, Gaussian SD %.2f pix\n",sig);
    mylog(mylogf,ggstr);

    //printf(" AREA smooth: %f\n",sum_2d_farray(tsm,0,xn,0,yn));

    // Store the smoothed frame back into a 3D stimulus array
    r = f3tensor(1,xn,1,yn,1,tn);
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	r[i+1][j+1][1] = tsm[i][j];

    free_2d_farray(tsm,xn);
    free_2d_farray(tdata,xn);

    //exit_error("WYETH HERE","Single frame smoothing: under development");

  }else if ((sig > 0.0) || (ts > 0.0)){

    mylog(mylogf,"  STIM_UTIL_BLUR\n");

    smf = get_stim_impulse(xn,yn,tn,1,1,1,0.0); // all zeros
    xc = (float)(xn/2 + 1);
    yc = (float)(yn/2 + 1);

    if ((sig > 0.0) && (ts == 0.0)){  // Spatial only
      sprintf(ggstr,"      Spatial smoothing, Gaussian SD %.2f pix\n",sig);
      mylog(mylogf,ggstr);
      for(i=1;i<=xn;i++)
	for(j=1;j<=yn;j++)
	  smf[i][j][1] = func_2d_gaussian((float)i,(float)j,xc,yc,sig,sig,1);
    }else if ((sig == 0.0) && (ts > 0.0)){  // Temporal only
      sprintf(ggstr,"      Temporal smooth, Maxwell peaks at %.2f pix\n",ts);
      mylog(mylogf,ggstr);
      //
      // WYETH - only using first half: assuming that second half is for -t
      // Spatial center of 1..n is (n/2 + 1)
      //
      tss = ts/sqrt(2.0);  // Convert time of peak to the maxwell parameter
      for(k=1;k<tn/2;k++)
	smf[xn/2+1][yn/2+1][k] = func_maxwell((float)(k-1),tss);
    }else{
      sprintf(ggstr,"      Spatial smooth, Gaussian SD %.2f pix\n",sig);
      mylog(mylogf,ggstr);
      sprintf(ggstr,"      Temporal smooth, Maxwell func peak %.2f pix\n",ts);
      mylog(mylogf,ggstr);
      mylog_exit(mylogf,"STIM_UTIL_BLUR  x and t smooth not impl'd yet!!!\n");
      /*
      for(i=1;i<=xn;i++)
	for(j=1;j<=yn;j++)
	  for(k=1;k<=tn;k++)
	    ;
      */
    }
    //write_3d_data_part("zz.SMF.3d",smf,1,xn,1,yn,1,tn,4,2,1); // 1=txy form

    //
    //  Compute FT of the smoothing function
    //
    smf_s = matrix(1,xn,1,2*yn);
    contort_real_3d_farray(smf,1,1,1,xn,yn,tn);
    rlft3(smf,smf_s,xn,yn,tn,1);
    //write_3d_data_part("zz.raw.3d",s->d,1,xn,1,yn,1,tn,4,2,1);  //1=txy form

    //
    //  Apply smoothing to 'sdata' and get smoothed result 'r'
    //
    compute_response_single_from_ft(sdata,smf,smf_s,xn,yn,tn,0,&r); // pflag 0

    // *** NOISE is introduced by smoothing, on the order of 0.000001
    // THIS CAN BECOME AMPLIFIED, if the stimulus is represented w/ 8 bits,
    // then the mean gray level can fluctuate between 127,128.

    //write_3d_data_part("zz.smoothed.3d",r,1,xn,1,yn,1,tn,4,2,1); //1=txy form
    //exit(0);

    free_f3tensor(smf,1,xn,1,yn,1,tn);
    free_matrix(smf_s,1,xn,1,2*yn);

    //write_3d_data_part("zz.sm.3d",s->d,1,xn,1,yn,1,tn,4,2,1);  //1=txy form
    //exit(0);

  }else{
    r = NULL;
  }

  return r;
}
/**************************************-**************************************/
/*                                                                           */
/*                            STIM_UTIL_ADD_MODSPOT                          */
/*                                                                           */
/*  Write a modulated spot into the stimulus array.                          */
/*                                                                           */
/*****************************************************************************/
void stim_util_add_modspot(ppl,data,x0,xn,y0,yn,z0,zn,tsamp,tscale,sscale)
     struct param_pair_list *ppl;
     float ***data;
     int x0,xn,y0,yn,z0,zn;
     int tsamp;      // frame increment
     float tscale;   // sec/frame
     float sscale;   // degr/pix
{
  int i,k;
  int xi,yi,x0pix,x1pix,y0pix,y1pix;
  int spot_type,varflag,seed,dt,distrib,ord,seqn,spotsz;
  float spotcon,stn,fps,v1,v2,*seqf,xc,yc,spotx,spoty,c;

  //printf("  STIM_UTIL_ADD_MODSPOT\n");
  //printf("    %d %d %d %d %d %d\n",x0,xn,y0,yn,z0,zn);

  if (tsamp != 1){
    printf("  tsamp = %d\n",tsamp);
    exit_error("STIM_UTIL_ADD_MODSPOT","tsamp > 1 not implemented");
  }

  spot_type = paramfile_get_int_param_or_exit(ppl,"spot_type");
  varflag   = paramfile_get_int_param_or_exit(ppl,"var");
  seed      = paramfile_get_int_param_or_exit(ppl,"seed");
  dt        = paramfile_get_int_param_or_exit(ppl,"dt");
  distrib   = param_getf_exit(ppl,"distrib");
  spotcon   = param_getf_exit(ppl,"spot_contrast");
  stn       = param_getf_exit(ppl,"stn");
  spotx     = param_getf_exit(ppl,"spot_x");
  spoty     = param_getf_exit(ppl,"spot_y");
  spotsz    = paramfile_get_int_param_or_exit(ppl,"spot_size_pix");

  xc = (float)(1+xn)/2.0 + spotx/sscale;
  yc = (float)(1+yn)/2.0 + spoty/sscale;

  // Pick x-limits for spot
  x0pix = my_rint(xc - (float)spotsz/2.0);
  x1pix = x0pix + spotsz - 1;
  if (x0pix < 1)
    x0pix = 1;
  if (x1pix > xn)
    x0pix = xn;

  // Pick y-limits for spot
  y0pix = my_rint(yc - (float)spotsz/2.0);
  y1pix = y0pix + spotsz - 1;
  if (y0pix < 1)
    y0pix = 1;
  if (y1pix > yn)
    y0pix = yn;

  //printf(" x:  %d %d\n",x0pix,x1pix);
  //printf(" y:  %d %d\n",y0pix,y1pix);
  
  v1 = 0.5 + 0.5*spotcon;
  if (v1 > 1.0)
    v1 = 1.0;
  v2 = 0.5 - 0.5*spotcon;
  if (v2 < 0.0)
    v2 = 0.0;

  if (distrib == 7)
    ord = param_geti_exit(ppl,"mseq_ord");
  else
    ord = 0;

  fps = 1.0/tscale;
  
  stm_get_ran_seq(NULL,distrib,ord,seed,fps,stn,dt,v1,v2,0.0,&seqf,&seqn);

  k = 0; 
  for(i=1;i<=zn;i++){
    if (((i-1)%dt) == 0){
      if (k >= seqn)
	exit_error("STIM_UTIL_ADD_MODSPOT","k too large");
      c = seqf[k];
      k += 1;
    }
    for(xi=x0pix;xi<=x1pix;xi++){    // Draw spot w/ value 'c'
      for(yi=y0pix;yi<=y1pix;yi++){
	data[xi][yi][i] = c;
      }
    }
  }

  myfree(seqf);
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_STIM_3D_NOISE                           */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_noise(xn,yn,zn,pixsize,pixdur,seed,amp,tensorflag)
     int xn,yn,zn;
     int pixsize,pixdur,seed;
     float amp;
     int tensorflag;
{
  int i,j,k;
  int x0,y0,z0,tseed;
  float ***data;

  if (tensorflag == 1){
    data = f3tensor(1,xn,1,yn,1,zn); // For 3d fft
    exit_error("GET_STIM_3D_NOISE","Not implemented yet");
    x0 = y0 = z0 = 1;
  }else{
    data = get_3d_farray(xn,yn,zn);
    x0 = y0 = z0 = 0;
  }

  tseed = seed;
  if (tseed > 0)
    tseed = -tseed;
  
  for(i=z0;i<(z0+zn);i++)
    if (((i-z0)%pixdur)==0){
      make_frame_binary_noise(data,x0,xn,y0,yn,i,pixsize,amp,&tseed);
    }else{ /* Copy last frame */
      for(j=0;j<xn;j++)
	for(k=0;k<yn;k++)
	  data[j+x0][k+y0][i] = data[j+x0][k+y0][i-1];
    }
  
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                         OFFSET_BLOCK_STIM_3D_NOISE                        */
/*                                                                           */
/*  Shift a block of pixels by some number of pixels, filling in with        */
/*  random pixels around the edge.                                           */
/*                                                                           */
/*  bx0,by0,bt0 - xyt coordinate of block (in pixsize units and frames)      */
/*  bxn,byn,btn - width, height, duration of block (in pixsize and frames)   */
/*  dx,dy       - x and y shift (in pixsize units)                           */
/*                                                                           */
/*****************************************************************************/
void offset_block_stim_3d_noise(data,xn,yn,zn,pixsize,pixdur,seed,amp,
				bx0,by0,bxn,byn,bt0,btn,dx,dy)
     float ***data;
     int xn,yn,zn;
     int pixsize,pixdur,seed;
     float amp;
     int bx0,by0,bxn,byn,bt0,btn,dx,dy;
{
  int i,j,k;
  int tseed,j1,j2,k1,k2,dj,dk;
  int mx0,mxn,my0,myn;
  float **tframe;

  tseed = seed;
  if (tseed > 0)
    tseed = -tseed;

  tframe = get_2d_farray(bxn*pixsize,byn*pixsize);
  j1 = bx0*pixsize;
  j2 = j1 + bxn*pixsize;
  k1 = by0*pixsize;
  k2 = k1 + byn*pixsize;
  dj = dx*pixsize;
  dk = dy*pixsize;


  for(i=bt0;i<(bt0+btn);i++){ /*** These are single frames ***/
    /*** Copy into the temp frame ***/
    for(j=j1;j<j2;j++)
      for(k=k1;k<k2;k++)
	tframe[j-j1][k-k1] = data[j][k][i];
    /*** Put it back shifted ***/
    for(j=j1;j<j2;j++)
      for(k=k1;k<k2;k++)
	data[j+dj][k+dk][i] = tframe[j-j1][k-k1];
    /*** Make new random vertical margin ***/
    if (dx != 0){
      if (dx > 0){
	mx0 = bx0*pixsize;
	mxn = dj;
	my0 = by0*pixsize;
	myn = byn*pixsize;
      }else{
	mx0 = bxn+dx - 1;
	mxn = dj;

	printf("****** WYETH what should my0 and myn be here? Using 0?\n");
	my0 = 0;  /* Error picked up by -Wall */
	myn = 0;

      }
      if (((i-bt0)%pixdur)==0)
	make_frame_binary_noise(data,mx0,mxn,my0,myn,i,pixsize,amp,&tseed);
      else{ /* Copy margin from last frame */
	for(j=mx0;j<(mx0+mxn);j++)
	  for(k=my0;k<(my0+myn);k++)
	    data[j][k][i] = data[j][k][i-1];
      }
    }
    if (dy != 0){
      exit_error("XXXXXXXXXXX","Not implemented yet");
    }
  }
  free_2d_farray(tframe,bxn*pixsize);
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_STIM_QUAD_SINE                            */
/*                                                                           */
/*  Return a stimulus that uses amplitude modulated sine and cosine          */
/*  waves to produce the appearance of motion of a single sine wave.         */
/*                                                                           */
/*    cos(ws*x - wt*t) = cos(ws*x)cos(wt*t) + sin(ws*x)sin(wt*t)             */
/*                                                                           */
/*  where wt = 2Pi*tf, ws = 2Pi*sf.                                          */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*    sf - spatial frequency (cyc/deg)                                       */
/*    tf - temporal frequency (cyc/sec)                                      */
/*    phase - phase offset (deg)                                             */
/*    theta - orientation of grating (deg)                                   */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_quad_sine(xn,yn,zn,sscale,tscale,tsamp,sf,tf,phase,theta,
			  contrast,maxlum)
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
     float sf,tf,phase,theta,contrast,maxlum;
{
  int i,k;
  float ***data,phi,mid,xc,yc;

  data = f3tensor(1,xn,1,yn,1,zn); /*** USE THIS FOR 3D FFT ***/

  xc = yc = 0.0;

  k = 0;
  for(i=1;i<=zn;i++){
    if (((i-1)%tsamp) == 0){
      if ((k%2)==0){ // Spatial SIN term
	make_frame_sinusoid(data,1,xn,1,yn,1,zn,i,sf,theta,sscale,0.0,xc,yc);
	phi = (float)i*tscale*tf * 2.0*M_PI;
	multiply_frame(data,1,xn,1,yn,1,zn,i,cos(phi));
      }else{ // Spatial COS term
	make_frame_sinusoid(data,1,xn,1,yn,1,zn,i,sf,theta,sscale,90.0,xc,yc);
	phi = (float)i*tscale*tf * 2.0*M_PI;
	multiply_frame(data,1,xn,1,yn,1,zn,i,sin(phi));
      }
      k += 1;
    }else
      make_frame_flat(data,1,xn,1,yn,1,zn,i,-1.0);
  }
  mid = (float)maxlum/2.0;
  multiply_3d_farray(data,1,xn,1,yn,1,zn,contrast*mid);
  add_const_3d_farray(data,1,xn,1,yn,1,zn,mid);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_STIM_3D_SINE                            */
/*                                                                           */
/*  Return a drifting grating which begins as a cosine centered in           */
/*  the middle of the spatial dimensions of the first frame (plus            */
/*  any specified phase offset.                                              */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*    sf - spatial frequency (cyc/deg)                                       */
/*    tf - temporal frequency (cyc/sec)                                      */
/*    phase - phase offset (deg)                                             */
/*    theta - orientation of grating (deg)                                   */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_sine(xn,yn,zn,sscale,tscale,tsamp,sf,tf,phase,theta,
			  contrast,maxlum)
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
     float sf,tf,phase,theta,contrast,maxlum;
{
  int i;
  float ***data,ph,mid,xc,yc;

  /***
    printf("%d %d %d %f %f %d %f %f %f %f %f %f\n",
    xn,yn,zn,sscale,tscale,tsamp,sf,tf,phase,theta,contrast,maxlum);***/

  data = f3tensor(1,xn,1,yn,1,zn); /*** USE THIS FOR 3D FFT ***/

  // WYETH - THIS IS NOT THE CONVENTIONAL CENTER!!
  printf("*** GET_STIM_3D_SINE - this is not the conventional center\n");
  xc = (float)xn/2.0 * sscale;
  yc = (float)yn/2.0 * sscale;

  for(i=1;i<=zn;i++){
    if (((i-1)%tsamp) == 0){
      ph = phase + (float)(i-1)*tscale*tf * 360.0;
      make_frame_sinusoid(data,1,xn,1,yn,1,zn,i,sf,theta,sscale,ph,xc,yc);
    }else
      make_frame_flat(data,1,xn,1,yn,1,zn,i,-1.0);
  }
  mid = (float)maxlum/2.0;
  multiply_3d_farray(data,1,xn,1,yn,1,zn,contrast*mid);
  add_const_3d_farray(data,1,xn,1,yn,1,zn,mid);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_STIM_SINE_INT_NOISE                         */
/*                                                                           */
/*  Return a drifting grating interleaved with binary noise of the           */
/*  specified "noise_contrast" and pixel size "npix".                        */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_sine_int_noise(xn,yn,zn,sscale,tscale,tsamp,sf,tf,phase,
				 theta,contrast,maxlum,noise_contrast,npix,
				 seed)
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
     float sf,tf,phase,theta,contrast,maxlum,noise_contrast;
     int npix;
     long seed;
{
  int i,j,k;
  int kk,tseed;
  float ***data,mid;

  printf("*** WARNING - WYETH changed usage of seed to send addr. ***\n");
  tseed = seed;
  if (tseed > 0)
    tseed = -tseed; /* Seed must be negative on first call */

  data = get_stim_3d_sine(xn,yn,zn,sscale,tscale,tsamp,sf,tf,phase,theta,
			  contrast,maxlum);
  mid = (float)maxlum/2.0;

  kk = 1+tsamp;
  if (zn >= 2){
    make_frame_binary_noise(data,1,xn,1,yn,kk,npix,noise_contrast*mid,
			    &tseed);
    for(i=1;i<=xn;i++)
      for(j=1;j<=yn;j++)
	data[i][j][kk] += mid;

    for(k=(kk+2*tsamp);k<=zn;k+=2*tsamp)
      for(i=1;i<=xn;i++)
	for(j=1;j<=yn;j++)
	  data[i][j][k] = data[i][j][kk];
  }
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                               GET_STIM_IMPULSE                            */
/*                                                                           */
/*  Return a single value in an otherwise 0.0 array.                         */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_impulse(xn,yn,zn,x,y,z,ampl)
     int xn,yn,zn,x,y,z;
     float ampl;
{
  int i,j,k;
  float ***data;

  data = f3tensor(1,xn,1,yn,1,zn); // USE THIS FOR 3D FFT
  
  for(i=1;i<=xn;i++)
    for(j=1;j<=yn;j++)
      for(k=1;k<=zn;k++)
	data[i][j][k] = 0.0;

  data[x][y][z] = ampl;

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_STIM_3D_PPL_DYNEDGE                         */
/*                                                                           */
/*  Return a dynamic edge stimulus.                                          */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_dynedge(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i;
  int seed;
  float ***data,t,t0,tn,ampl,bgval;
  float amp,sd,maxlum,contrast,mid;

  sd       = param_getf_exit(ppl,"sd");
  maxlum   = param_getf_exit(ppl,"maxlum");
  contrast = param_getf_exit(ppl,"contrast");
  seed     =   paramfile_get_int_param_or_exit(ppl,"seed");

  data = f3tensor(1,xn,1,yn,1,zn); // For 3D FFT

  if (seed > 0)
    seed = -seed;

  amp = 1.0;

  for(i=1;i<=zn;i++){
    make_frame_edge_01(data,1,xn,1,yn,1,zn,i,amp,sd,&seed);
  }
  mid = (float)maxlum/2.0;
  multiply_3d_farray(data,1,xn,1,yn,1,zn,contrast*mid);
  add_const_3d_farray(data,1,xn,1,yn,1,zn,mid);

  // WYETH - this prevents normalisation of 
  data[1][1][1] = -1.0;
  data[1][1][2] =  1.0;

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_STIM_3D_PPL_DISK                          */
/*                                                                           */
/*  Return a circular disk.                                                  */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_disk(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i;
  float ***data,t,t0,tn,size,ampl,bgval;

  size     = param_getf_exit(ppl,"size");
  ampl     = param_getf_exit(ppl,"ampl");
  bgval   = param_getf_exit(ppl,"bgval");
  t0       = param_getf_exit(ppl,"st0");
  tn       = param_getf_exit(ppl,"stn");

  data = f3tensor(1,xn,1,yn,1,zn); // USE THIS FOR 3D FFT

  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;
    if ((((i-1)%tsamp) == 0) && (t >= t0) && (t < t0+tn)){
      make_frame_disk_flat(data,1,xn,1,yn,1,zn,i,size,sscale,ampl,bgval);
    }else
      make_frame_flat(data,1,xn,1,yn,1,zn,i,bgval);
  }
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_STIM_3D_PPL_SINE                          */
/*                                                                           */
/*  Return a drifting grating which begins as a cosine centered in           */
/*  the middle of the spatial dimensions of the first frame (plus            */
/*  any specified phase offset.                                              */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*    sf - spatial frequency (cyc/deg)                                       */
/*    tf - temporal frequency (cyc/sec)                                      */
/*    phase - phase offset (deg)                                             */
/*    theta - orientation of grating (deg)                                   */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_sine(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i;
  int i0,binflag,aptype,spot_type,phjit_seed,*iseq,meanzero,cntrph,phorg;
  int customi;
  float ***data,ph,mid,sf,tf,phase,theta,contrast,maxlum,t0,tn,t,bgamp,size;
  float t0_2,tn_2,ph_off,phase2,tf2,sf2,theta2,size2,cx2,cy2,xc2,yc2;
  float binthr,sz,bgval,blankval,phase_t,pht,a,tfx,tfy;
  float phjit_sd,phjit_tsd,jsd,jsdt,*phj,blnk;
  float ff,ffrac,phtau,cx,cy,xc,yc;
  char *phjit;

  cx       = param_getf_dflt(ppl,"cx",0.0);
  cy       = param_getf_dflt(ppl,"cy",0.0);
  maxlum   = param_getf_exit(ppl,"maxlum");
  sf       = param_getf_exit(ppl,"sf");

  if (paramfile_test_param(ppl,"tfx")){
    //
    //  Allow the user to specify 'tfx' and 'tfy' instead of
    //     'direction' and 'tf'
    //
    tfx    = param_getf_exit(ppl,"tfx");
    tfy    = param_getf_exit(ppl,"tfy");
    tf     = sqrt(tfx*tfx + tfy*tfy);
    theta  = 180.0 * atan2(tfy,tfx)/M_PI;
    while(theta < 0.0)
      theta += 360.0;
    while(theta >= 360.0)
      theta -= 360.0;

    //printf(" tfx,y  %f  %f   dir %f  tf %f\n",tfx,tfy,theta,tf);
    
  }else{
    theta    = param_getf_exit(ppl,"direction");
    tf       = param_getf_exit(ppl,"tf");
  }
  contrast = param_getf_exit(ppl,"contrast");
  phase    = param_getf_exit(ppl,"phase");
  phorg    = param_geti_dflt(ppl,"phase_origin",0);
  cntrph   = param_geti_dflt(ppl,"counterphase",0);
  phase_t  = param_getf_dflt(ppl,"t_phase",0.0);
  t0       = param_getf_exit(ppl,"st0");
  tn       = param_getf_exit(ppl,"stn");
  bgval    = param_getf_exit(ppl,"bgval");
  size     = param_getf_exit(ppl,"size");
  aptype   = param_geti_exit(ppl,"aptype");
  phtau    = param_getf_dflt(ppl,"stim_phtau",0.006);
  binflag  = paramfile_test_param(ppl,"square_thresh");
  binthr   = param_getf_dflt(ppl,"square_thresh",0.0);
  meanzero = param_geti_dflt(ppl,"stim_mean_zero",0);
  blankval = param_getf_dflt(ppl,"stim_blank_val",0.0);

  customi  = param_geti_dflt(ppl,"custom",0);  // WYETH HACK for custom stimulus

  tn_2     = param_getf_dflt(ppl,"stn_2",-1.0);
  if (tn_2 > 0.0){  // WYETH - presumably must be greater than end of s1?
    t0_2   = param_getf_exit(ppl,"st0_2");
    sf2    = param_getf_dflt(ppl,"sf2",sf);
    tf2    = param_getf_dflt(ppl,"tf2",tf);
    theta2 = param_getf_dflt(ppl,"theta2",theta);
    size2  = param_getf_dflt(ppl,"size2",size);
    cx2    = param_getf_dflt(ppl,"cx2",cx);
    cy2    = param_getf_dflt(ppl,"cy2",cy);

    if (paramfile_test_param(ppl,"phase_offset")){
      ph_off = param_getf_exit(ppl,"phase_offset");
      phase2 = phase + ph_off;
    }else{
      phase2 = param_getf_dflt(ppl,"phase2",phase);
    }
  }else{
    t0_2 = -1.0;
  }

  // REMOVE THIS EVENTUALLY
  //  (WYETH maybe blankval should go from 0..maxlum ???)  2010 July
  if (blankval < 0.0)
    exit_error("*** WYETH ****  Blankval no longers goes neg, use [0..1]\n");

  blnk = 2.0*blankval - 1.0;

  // Phase jitter, _sd (deg vis ang) _tsd (s)
  phjit = paramfile_get_char_param_default(ppl,"phase_jit",NULL);

  if (phjit == NULL){
    phjit_sd = 0; // Do nothing
  }else if (strcmp(phjit,"gfgwn")==0){
    phjit_sd = paramfile_get_float_param_default(ppl,"phase_jit_sd",0.0);
    phjit_tsd = param_getf_exit(ppl,"phase_jit_tsd");
    phjit_seed = paramfile_get_int_param_or_exit(ppl,"phase_jit_seed");

    jsd  = phjit_sd * sf * 360.0;  // deg * cyc/deg * 360 = deg of Phase
    jsdt = phjit_tsd / tscale;     // s * frames/s = frames
  
    phj = gaussian_corr_noise(jsdt,0.0,jsd,zn,phjit_seed);
    //append_farray_plot("zz.phjit.pl","phjit",phj,zn,1);
  }else if (strcmp(phjit,"rwalk")==0){
    phjit_sd = paramfile_get_float_param_default(ppl,"phase_jit_step",0.0);
    phjit_seed = paramfile_get_int_param_or_exit(ppl,"phase_jit_seed");

    iseq = myrand_get_std_bin_seq(zn,phjit_seed);  // Integer sequence, 0,1
    phj = (float *)myalloc(zn*sizeof(float));
    phj[0] = 0.0;

    jsd  = phjit_sd * sf * 360.0;  // deg * cyc/deg * 360 = deg of Phase

    for(i=1;i<zn;i++){
      if (iseq[i] == 1)
	phj[i] = phj[i-1] + jsd;
      else
	phj[i] = phj[i-1] - jsd;
    }
    myfree(iseq);
    /*
    printf("  FINAL PHASE:  %.4f  = %.4f cycles\n",phj[zn-1],phj[zn-1]/360.0);
    printf("           SF:  %f cyc/deg\n",sf);
    printf("         DIST:  %f deg.v.a.\n",phj[zn-1]/360.0/sf);
    printf("           zn:  %d\n",zn);
    */
  }
  if (phjit != NULL)
    myfree(phjit);

  /*x0 = 1 + xn/2 + (cx - w/2.0 - barw/2.0)/sscale;*/
  /* Pix left edge of xrange */

  // Values of size and center in PIXELs w.r.t. 0..xn, 0..yn
  sz = size/sscale;
  xc = (float)(xn-1)/2.0 + cx/sscale;  // Center w.r.t. 0..xn-1 (pix)
  yc = (float)(yn-1)/2.0 + cy/sscale;  // Center w.r.t. 0..yn-1 (pix)

  if (tn_2 > 0.0){
    xc2 = (float)(xn-1)/2.0 + cx2/sscale;  // Center w.r.t. 0..xn-1 (pix)
    yc2 = (float)(yn-1)/2.0 + cy2/sscale;  // Center w.r.t. 0..yn-1 (pix)
  }

  if (phtau > 0.0)
    ffrac = exp(-tscale/phtau);
  else
    ffrac = 0.0;
  ff = ffrac;   // Ampl for next empty frame

  //printf("%d %d %d %f %f %d %f %f %f %f %f %f\n",
  //xn,yn,zn,sscale,tscale,tsamp,sf,tf,phase,theta,contrast,maxlum);

  data = f3tensor(1,xn,1,yn,1,zn); // USE THIS FOR 3D FFT

  //
  //  WYETH - NEED A BETTER SYSTEM FOR SPECIFYING STIMULUS RANGE and bgamp...
  //  bgval should be relative to maxlum?, whereas it is treated absolute
  //  between 0..1 here, regardless of maxlum
  //
  bgamp = bgval * 2.0 - 1.0; // Convert 0..1 --> -1..1

  i0 = 1 + my_rint(t0 / tscale);
  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;
    if (((i-1)%tsamp) == 0){
      if ((t >= t0) && (t < t0+tn)){

	if (cntrph == 0){
	  ph = phase + (float)(i-i0)*tscale*tf * 360.0;
	  a = 1.0;  // Amplitude scale factor for sinusoid
	}else{
	  // counterphase
	  pht = phase_t/360.0 + (float)(i-i0)*tscale*tf;  // mod. phase (cyc)
	  a = sin(2.0*M_PI*pht); // modulated contrast
	  ph = phase;
	}

	if (customi == 1){ // Change theta, 
	  theta += 0.35;  //  in 512 steps, change by 180 deg
	  printf("  (stim_util) theta = %f\n",theta);
	}else if (customi == 2){ // Change sf 
	  sf = 0.3 + 5.7 * (float)i/(float)zn;
	  printf("  (stim_util) sf = %f\n",sf);
	}else if (customi == 3){ // Change sf 
	  xc =  20.0 + 80.0 * (float)i/(float)zn;
	  printf("  (stim_util) xc = %f\n",xc);
	}

	if (phjit_sd > 0.0){
	  ph += phj[i-i0];
	}

	// WYETH CONFIX 2014 Nov
	make_frame_sinusoid_scale(data,1,xn,1,yn,1,zn,i,sf,theta,sscale,ph,
				  cx,cy,a*contrast,phorg);
	//cx,cy,a);

	if (binflag){
	  // WYETH CONFIX 2014 Nov, was ... 1.0,-1.0);
	  make_frame_binarize(data,1,xn,1,yn,1,zn,i,binthr,contrast,-contrast);
	}
	if (aptype == 3){
	  multiply_frame_gaussian(data,1,1,1,xn,yn,zn,i,xc,yc,size,sscale);
	}
      }else if ((t >= t0_2) && (t < t0_2+tn_2)){

	ph = phase2 + (float)(i-i0)*tscale*tf2 * 360.0;

	// WYETH CONFIX 2014 Nov
	//make_frame_sinusoid(data,1,xn,1,yn,1,zn,i,sf2,theta2,sscale,ph,
	//cx2,cy2);
	make_frame_sinusoid_scale(data,1,xn,1,yn,1,zn,i,sf2,theta2,sscale,ph,
				  cx2,cy2,contrast,phorg);

	if (binflag){
	  // WYETH CONFIX 2014 Nov, was ... 1.0,-1.0);
	  make_frame_binarize(data,1,xn,1,yn,1,zn,i,binthr,contrast,-contrast);
	}
	if (aptype == 3){
	  multiply_frame_gaussian(data,1,1,1,xn,yn,zn,i,xc2,yc2,size2,sscale);
	}

      }else
	make_frame_flat(data,1,xn,1,yn,1,zn,i,bgamp);
      ff = ffrac; // Reset phosphor decay ampl
    }else{
      if ((ff > 0.0)&&(i>1)){ // scale down for phosphor decay
	make_frame_scale_frame_zmid(data,1,xn,1,yn,1,zn,i,i-1,ff,blnk);
      }else
	make_frame_flat(data,1,xn,1,yn,1,zn,i,blnk); // typ. black, -1.0
      ff *= ffrac; // ampl for next time
    }

    if (aptype != 3){
      //printf("bgamp = %f\n",bgamp);
      frame_apply_aperture(data,1,xn,1,yn,1,zn,i,aptype,xc,yc,sz,sz,theta,
			   bgamp);
    }
  }
  if (meanzero == 0){
    mid = (float)maxlum/2.0;
    //multiply_3d_farray(data,1,xn,1,yn,1,zn,contrast*mid); // WYETH CONFIX
    multiply_3d_farray(data,1,xn,1,yn,1,zn,mid);
    add_const_3d_farray(data,1,xn,1,yn,1,zn,mid);
  }else{
    // WYETH CONFIX 2014 Nov
    //if (contrast*maxlum != 1.0)
    //multiply_3d_farray(data,1,xn,1,yn,1,zn,contrast*maxlum);
    if (maxlum != 1.0)
      multiply_3d_farray(data,1,xn,1,yn,1,zn,maxlum);
  }

  // Check for 'modspot'
  spot_type = paramfile_get_int_param_default(ppl,"spot_type",0);
  if (spot_type > 0){
    stim_util_add_modspot(ppl,data,1,xn,1,yn,1,zn,tsamp,tscale,sscale);
  }

  if (phjit_sd > 0.0){
    myfree(phj);
  }

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_STIM_3DB_PPL_SINE                          */
/*                                                                           */
/*  BINOCULAR                                                                */
/*                                                                           */
/*****************************************************************************/
void get_stim_3db_ppl_sine(ppl,xn,yn,zn,sscale,tscale,tsamp,rdl,rdr)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
     float ****rdl,****rdr;
{
  int i;
  int i0,binflag,aptype,res,make_right,make_left,meanzero,phdispfl;
  float ***dl,***dr,mid,maxlum,t0,tn,t,bgamp,size,dph;
  float phase_l,phase_r,phl,phr,cxl,cxr,cyl,cyr;
  float xcl,xcr,ycl,ycr,thetal,thetar,sfl,sfr,tfl,tfr,conl,conr;
  float binthr,sz,bgval,ff,ffrac,phtau,cx,cy,xc,yc;
  char *monoflag;

  param_getf_lrd("GET_STIM_3DB_PPL_SINE",ppl,"cx",&cxl,&cxr);
  param_getf_lrd("GET_STIM_3DB_PPL_SINE",ppl,"cy",&cyl,&cyr);
  param_getf_lrd("GET_STIM_3DB_PPL_SINE",ppl,"phase",&phase_l,&phase_r);
  param_getf_lrd("GET_STIM_3DB_PPL_SINE",ppl,"direction",&thetal,&thetar);
  param_getf_lrd("GET_STIM_3DB_PPL_SINE",ppl,"sf",&sfl,&sfr);
  param_getf_lrd("GET_STIM_3DB_PPL_SINE",ppl,"tf",&tfl,&tfr);
  param_getf_lrd("GET_STIM_3DB_PPL_SINE",ppl,"contrast",&conl,&conr);

  phdispfl = param_geti_dflt(ppl,"phase_disp_flag",1);
  maxlum   = param_getf_exit(ppl,"maxlum");
  t0       = param_getf_exit(ppl,"st0");
  tn       = param_getf_exit(ppl,"stn");
  bgval    = param_getf_exit(ppl,"bgval");
  size     = param_getf_exit(ppl,"size");
  aptype   = param_geti_exit(ppl,"aptype");
  phtau    = param_getf_dflt(ppl,"stim_phtau",0.006);
  monoflag = param_getc_dflt(ppl,"stim_monocular","both"); // "left","right"
  meanzero = param_geti_dflt(ppl,"stim_mean_zero",0);
  binflag  = paramfile_test_param(ppl,"square_thresh");
  binthr   = paramfile_get_float_param_default(ppl,"square_thresh",0.0);

  make_left = make_right = 1;
  if (strcmp(monoflag,"left")==0)
    make_right = 0;
  else if (strcmp(monoflag,"right")==0)
    make_left = 0;

  sz = size/sscale;

  xcl = (float)(xn-1)/2.0 + cxl/sscale;  // Center w.r.t. 0..xn-1 (pix)
  ycl = (float)(yn-1)/2.0 + cyl/sscale;  // Center w.r.t. 0..yn-1 (pix)
  xcr = (float)(xn-1)/2.0 + cxr/sscale;
  ycr = (float)(yn-1)/2.0 + cyr/sscale;

  if (phtau > 0.0)
    ffrac = exp(-tscale/phtau);
  else
    ffrac = 0.0;
  ff = ffrac;   // Ampl for next empty frame

  dl = f3tensor(1,xn,1,yn,1,zn); // For 3D FFT
  dr = f3tensor(1,xn,1,yn,1,zn); // For 3D FFT

  bgamp = bgval * 2.0 - 1.0; // Convert 0..1 --> -1..1
  i0 = 1 + my_rint(t0 / tscale);
  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;
    if (((i-1)%tsamp) == 0){
      if ((t >= t0) && (t < t0+tn)){

	//
	//  WYETH - new way - this allows the disparity to be interpreted
	//   as a maximum horizontal disparity that stays consistent regardless
	//   of the ori/direction.
	//
	if (phdispfl == 1){
	  dph = phase_l * cos(thetal * M_PI/180.0);
	  //printf("  theta_l %f  dph %f\n",thetal,dph);
	}else
	  dph = phase_l;
	phl = dph + (float)(i-i0)*tscale*tfl * 360.0;
	//phl = phase_l + (float)(i-i0)*tscale*tfl * 360.0;

	if (phdispfl == 1){
	  dph = phase_r * cos(thetar * M_PI/180.0);
	  //printf("    tht_r %f  dph %f\n",thetar,dph);
	}else
	  dph = phase_r;
	phr = dph + (float)(i-i0)*tscale*tfr * 360.0;
	//phr = phase_r + (float)(i-i0)*tscale*tfr * 360.0;

	if (make_left == 1)
	  make_frame_sinusoid(dl,1,xn,1,yn,1,zn,i,sfl,thetal,sscale,phl,
			      cxl,cyl);
	else
	  make_frame_flat(dl,1,xn,1,yn,1,zn,i,bgamp);

	if (make_right == 1)
	  make_frame_sinusoid(dr,1,xn,1,yn,1,zn,i,sfr,thetar,sscale,phr,
			      cxr,cyr);
	else
	  make_frame_flat(dr,1,xn,1,yn,1,zn,i,bgamp);

	if (binflag){
	  if (make_left == 1)
	    make_frame_binarize(dl,1,xn,1,yn,1,zn,i,binthr,1.0,-1.0);
	  if (make_right == 1)
	    make_frame_binarize(dr,1,xn,1,yn,1,zn,i,binthr,1.0,-1.0);
	}
	if (aptype == 3){
	  if (make_left == 1)
	    multiply_frame_gaussian(dl,1,1,1,xn,yn,zn,i,xcl,ycl,size,sscale);
	  if (make_right == 1)
	    multiply_frame_gaussian(dr,1,1,1,xn,yn,zn,i,xcr,ycr,size,sscale);
	}
      }else{
	make_frame_flat(dl,1,xn,1,yn,1,zn,i,bgamp);
	make_frame_flat(dr,1,xn,1,yn,1,zn,i,bgamp);
      }
      ff = ffrac; // Reset phosphor decay ampl
    }else{
      if ((ff > 0.0)&&(i>1)){ // scale down for phosphor decay
	make_frame_scale_frame_zmid(dl,1,xn,1,yn,1,zn,i,i-1,ff,-1.0);
	make_frame_scale_frame_zmid(dr,1,xn,1,yn,1,zn,i,i-1,ff,-1.0);
      }else{
	make_frame_flat(dl,1,xn,1,yn,1,zn,i,-1.0); // black
	make_frame_flat(dr,1,xn,1,yn,1,zn,i,-1.0); // black
      }
      ff *= ffrac; // ampl for next time
    }

    if (aptype != 3){
      if (make_left == 1)
	frame_apply_aperture(dl,1,xn,1,yn,1,zn,i,aptype,xcl,ycl,sz,sz,0.0,
			     bgamp);
      if (make_right == 1)
	frame_apply_aperture(dr,1,xn,1,yn,1,zn,i,aptype,xcr,ycr,sz,sz,0.0,
			     bgamp);
    }
  }

  if (meanzero == 0){
    mid = (float)maxlum/2.0;
    multiply_3d_farray(dl,1,xn,1,yn,1,zn,conl*mid);
    multiply_3d_farray(dr,1,xn,1,yn,1,zn,conr*mid);
    add_const_3d_farray(dl,1,xn,1,yn,1,zn,mid);
    add_const_3d_farray(dr,1,xn,1,yn,1,zn,mid);
  }else{
    //
    //  Using 'contrast' here allows control stimulus param 'contrast' 0 to
    //  still work, rather than changing to 'maxlum' 0.
    //
    if (conl*maxlum != 1.0)
      multiply_3d_farray(dl,1,xn,1,yn,1,zn,conl*maxlum);
    if (conr*maxlum != 1.0)
      multiply_3d_farray(dr,1,xn,1,yn,1,zn,conr*maxlum);
  }

  *rdl = dl;
  *rdr = dr;
}
/**************************************-**************************************/
/*                                                                           */
/*                      GET_STIM_3DB_PPL_RANDOT_STEREOGRAM                   */
/*                                                                           */
/*  BINOCULAR                                                                */
/*                                                                           */
/*****************************************************************************/
void get_stim_3db_ppl_randot_stereogram(ppl,xn,yn,zn,sscale,tscale,tsamp,
					rdl,rdr)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
     float ****rdl,****rdr;
{
  int i,j,k;
  float ***dl,***dr,t,t0,tn,bgval,**tframl,**tframr;

  //
  //  PREP
  //
  stimx_dots_stereogram(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL,NULL);

  t0        = param_getf_exit(ppl,"st0");
  tn        = param_getf_exit(ppl,"stn");
  bgval     = param_getf_exit(ppl,"bgval");

  dl = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT
  dr = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  for(i=1;i<=zn;i++){
    if (((i-1)%tsamp) == 0){

      stimx_dots_stereogram(ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,
			    &tframl,&tframr);

      //
      //  Copy this frame in and free it [WYETH - CHANGE THIS]
      //
      for(j=0;j<xn;j++){
	for(k=0;k<yn;k++){
	  dl[j+1][k+1][i] = tframl[j][k];
	  dr[j+1][k+1][i] = tframr[j][k];
	}
      }
      free_2d_farray(tframl,xn);
      free_2d_farray(tframr,xn);

    }else{
      make_frame_flat(dl,1,xn,1,yn,1,zn,i,bgval);
      make_frame_flat(dr,1,xn,1,yn,1,zn,i,bgval);
    }
  }
  // Clean up
  stimx_dots_stereogram(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL,NULL);

  *rdl = dl;
  *rdr = dr;
}
/**************************************-**************************************/
/*                                                                           */
/*                          STIM_UTIL_CHECK_SET_NSEQ                         */
/*                                                                           */
/*****************************************************************************/
void stim_util_check_set_nseq(rnseq,n)
     int *rnseq,n;
{
  if (*rnseq > 0){
    if (n > 0){
      if (n != *rnseq)
	exit_error("STIM_UTIL_CHECK_SET_NSEQ","conflicting number of values");
    }
    // Numbers match
  }else{
    if (n > 0){
      *rnseq = n;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_STIM_3D_PPL_SINE_SEQ                        */
/*                                                                           */
/*  Sequence of sinusoids                                                    */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_sine_seq(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i;
  int i0,si,binflag,aptype,dt,cnt,seqdone;
  float ***data,ph,mid,sf,tf,phase,theta,contrast,maxlum,t0,tn,t,bgamp,size;
  float binthr,sz,bgval,dir0;
  float ff,ffrac,phtau,cx,cy,xc,yc;
  
  int  nseq,seq_dir,seq_ph,seq_sf,seq_tf,seq_con,seq_sz;

  cx       = paramfile_get_float_param_default(ppl,"cx",0.0);
  cy       = paramfile_get_float_param_default(ppl,"cy",0.0);
  maxlum   = param_getf_exit(ppl,"maxlum");
  dir0     = param_getf_exit(ppl,"direction");
  sf       = param_getf_exit(ppl,"sf");
  tf       = param_getf_exit(ppl,"tf");
  contrast = param_getf_exit(ppl,"contrast");
  phase    = param_getf_exit(ppl,"phase");
  t0       = param_getf_exit(ppl,"st0");
  tn       = param_getf_exit(ppl,"stn");
  dt       =   paramfile_get_int_param_or_exit(ppl,"dt");
  bgval    = param_getf_exit(ppl,"bgval");
  size     = param_getf_exit(ppl,"size");
  aptype   = paramfile_get_int_param_or_exit(ppl,"aptype");
  phtau    = paramfile_get_float_param_default(ppl,"stim_phtau",0.006);
  binflag  = paramfile_test_param(ppl,"square_thresh");
  binthr   = paramfile_get_float_param_default(ppl,"square_thresh",0.0);

  // WYETH - should allow other params to vary here
  seq_dir = paramfile_count_values_param(ppl,"seq_dir"); // -1 if not found
  seq_ph  = paramfile_count_values_param(ppl,"seq_ph");  // -1 if not found
  seq_sf  = paramfile_count_values_param(ppl,"seq_sf");  // -1 if not found
  seq_tf  = paramfile_count_values_param(ppl,"seq_tf");  // -1 if not found
  seq_con = paramfile_count_values_param(ppl,"seq_con"); // -1 if not found
  seq_sz  = paramfile_count_values_param(ppl,"seq_sz");  // -1 if not found

  nseq = 0;
  stim_util_check_set_nseq(&nseq,seq_dir);
  stim_util_check_set_nseq(&nseq,seq_ph);
  stim_util_check_set_nseq(&nseq,seq_sf);
  stim_util_check_set_nseq(&nseq,seq_tf);
  stim_util_check_set_nseq(&nseq,seq_con);
  stim_util_check_set_nseq(&nseq,seq_sz);
  if (nseq <= 0)
    exit_error("GET_STIM_3D_PPL_SIN_SEQ","No sequence values");
  
  //printf("  Number of values: %d\n",nseq);

  sz = size/sscale;
  xc = (float)(xn-1)/2.0 + cx/sscale;  // Center w.r.t. 0..xn-1 (pix)
  yc = (float)(yn-1)/2.0 + cy/sscale;  // Center w.r.t. 0..yn-1 (pix)


  ffrac = exp(-tscale/phtau);
  ff = ffrac;   // Ampl for next empty frame

  //printf("%d %d %d %f %f %d %f %f %f %f %f %f\n",
  //xn,yn,zn,sscale,tscale,tsamp,sf,tf,phase,theta,contrast,maxlum);

  data = f3tensor(1,xn,1,yn,1,zn); // USE THIS FOR 3D FFT

  bgamp = bgval * 2.0 - 1.0; // Convert 0..1 --> -1..1
  mid = (float)maxlum/2.0;

  i0 = 1 + my_rint(t0 / tscale);
  si = 0; // Sequence index
  seqdone = 0;
  cnt = dt;
  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;
    if (((i-1)%tsamp) == 0){
      if ((t >= t0) && (t < t0+tn) && (!seqdone)){
	if (cnt == dt){ // Pick a new stimulus pattern

	  if (seq_dir > 0){
	    theta = paramfile_get_nth_float_param_or_exit(ppl,"seq_dir",si);
	    theta += dir0;
	  }
	  if (seq_ph > 0)
	    phase = paramfile_get_nth_float_param_or_exit(ppl,"seq_ph",si);
	  if (seq_sf > 0)
	    sf = paramfile_get_nth_float_param_or_exit(ppl,"seq_sf",si);
	  if (seq_tf > 0)
	    tf = paramfile_get_nth_float_param_or_exit(ppl,"seq_tf",si);
	  if (seq_con > 0)
	    contrast = paramfile_get_nth_float_param_or_exit(ppl,"seq_con",si);
	  if (seq_sz > 0){
	    size = paramfile_get_nth_float_param_or_exit(ppl,"seq_sz",si);
	    sz = size/sscale;
	  }

	  si += 1;
	  cnt = 0;
	}

	ph = phase + (float)(i-i0)*tscale*tf * 360.0;
	make_frame_sinusoid(data,1,xn,1,yn,1,zn,i,sf,theta,sscale,ph,cx,cy);
	if (binflag){
	  make_frame_binarize(data,1,xn,1,yn,1,zn,i,binthr,1.0,-1.0);
	}
	if (aptype == 3){
	  multiply_frame_gaussian(data,1,1,1,xn,yn,zn,i,0.0,0.0,size,sscale);
	}
	cnt += 1;
	if ((si >= nseq) && (cnt >= dt))
	  seqdone = 1;  // Done w/ sequence
      }else
	make_frame_flat(data,1,xn,1,yn,1,zn,i,bgamp);
      ff = ffrac; // Reset phosphor decay ampl
    }else{
      if ((ff > 0.0)&&(i>1)){ // scale down for phosphor decay
	make_frame_scale_frame_zmid(data,1,xn,1,yn,1,zn,i,i-1,ff,-1.0);
      }else
	make_frame_flat(data,1,xn,1,yn,1,zn,i,-1.0); // black
      ff *= ffrac; // ampl for next time
    }

    if (aptype != 3)
      frame_apply_aperture(data,1,xn,1,yn,1,zn,i,aptype,xc,yc,sz,sz,0.0,bgamp);
    
    multiply_frame(data,1,xn,1,yn,1,zn,i,contrast*mid);
  }
  add_const_3d_farray(data,1,xn,1,yn,1,zn,mid);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_STIM_3D_PPL_SOS                           */
/*                                                                           */
/*  Sum of sinusoids (VSK)                                                   */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_sos(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i;
  int i0,aptype,vskn,pcomp;
  float ***data,mid,theta,sf,phase,contrast,maxlum,t0,tn,t,bgampl,size;
  float sz,pratio,f,toff,bgval,*fot,cx,cy,xc,yc;

  cx       = paramfile_get_float_param_default(ppl,"cx",0.0);
  cy       = paramfile_get_float_param_default(ppl,"cy",0.0);

  maxlum   = param_getf_exit(ppl,"maxlum");
  theta    = param_getf_exit(ppl,"orientation");
  sf       = param_getf_exit(ppl,"sf");
  contrast = param_getf_exit(ppl,"contrast");
  phase    = param_getf_exit(ppl,"phase");
  t0       = param_getf_exit(ppl,"st0");
  tn       = param_getf_exit(ppl,"stn");
  bgval    = param_getf_exit(ppl,"bgval");
  size     = param_getf_exit(ppl,"size");
  aptype   =   paramfile_get_int_param_or_exit(ppl,"aptype");

  vskn     =   paramfile_get_int_param_or_exit(ppl,"vskn");
  pcomp    =   paramfile_get_int_param_or_exit(ppl,"pcomp",-1);
  pratio   = param_getf_exit(ppl,"pratio",1.0);
  toff     = param_getf_exit(ppl,"toff");

  /***printf("%d %d %d %f %f %d %f %f %f %f %f %f\n",
      xn,yn,zn,sscale,tscale,tsamp,sf,tf,phase,theta,contrast,maxlum);***/

  data = f3tensor(1,xn,1,yn,1,zn); /*** USE THIS FOR 3D FFT ***/

  xc = (float)(xn-1)/2.0 + cx/sscale;  // Center w.r.t. 0..xn-1 (pix)
  yc = (float)(yn-1)/2.0 + cy/sscale;  // Center w.r.t. 0..yn-1 (pix)

  bgampl = bgval * 2.0 - 1.0; /* Convert 0..1 --> -1..1 */
  i0 = 1 + my_rint(t0 / tscale);
  fot = (float *)myalloc(zn*sizeof(float));
  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;
    if (((i-1)%tsamp) == 0){
      if ((t >= t0) && (t < t0+tn)){
	make_frame_sinusoid(data,1,xn,1,yn,1,zn,i,sf,theta,sscale,phase,cx,cy);
	f = get_vsk_sum_of_sinusoids(vskn,pcomp,pratio,(t-t0+toff));
	fot[i] = f;
	multiply_frame(data,1,xn,1,yn,1,zn,i,f);
	if (aptype == 3)
	  multiply_frame_gaussian(data,1,1,1,xn,yn,zn,i,0.0,0.0,size,sscale);
      }else
	make_frame_flat(data,1,xn,1,yn,1,zn,i,bgampl);
    }else
      make_frame_flat(data,1,xn,1,yn,1,zn,i,-1.0);

    sz = size/sscale;
    if (aptype != 3)
      frame_apply_aperture(data,1,xn,1,yn,1,zn,i,aptype,xc,yc,
			   sz,sz,0.0,bgampl);
  }
  mid = (float)maxlum/2.0;
  multiply_3d_farray(data,1,xn,1,yn,1,zn,contrast*mid);
  add_const_3d_farray(data,1,xn,1,yn,1,zn,mid);

  printf("APPENDING time-varying modulation to 'zz.VSK.pl'\n");
  append_farray_plot("zz.VSK.pl","VSKt",fot,zn,1);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_STIM_3D_PPL_CPSINE                         */
/*                                                                           */
/*  Return a drifting grating which begins as a cosine centered in           */
/*  the middle of the spatial dimensions of the first frame (plus            */
/*  any specified phase offset.                                              */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*    sf - spatial frequency (cyc/deg)                                       */
/*    tf - temporal frequency (cyc/sec)                                      */
/*    phase - phase offset (deg)                                             */
/*    theta - orientation of grating (deg)                                   */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_cpsine(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i;
  int i0;
  float ***data,ph,con,mid;
  float sf,tf,phase,theta,contrast,maxlum,tphase,t0,tn,t,cx,cy;

  cx       = paramfile_get_float_param_default(ppl,"cx",0.0);
  cy       = paramfile_get_float_param_default(ppl,"cy",0.0);

  maxlum   = param_getf_exit(ppl,"maxlum");
  theta    = param_getf_exit(ppl,"orientation");
  sf       = param_getf_exit(ppl,"sf");
  tf       = param_getf_exit(ppl,"tf");
  contrast = param_getf_exit(ppl,"contrast");
  phase    = param_getf_exit(ppl,"phase");
  tphase   = param_getf_exit(ppl,"tphase");
  t0       = param_getf_exit(ppl,"st0");
  tn       = param_getf_exit(ppl,"stn");
  /*size     = param_getf_exit(ppl,"size");*/

  data = f3tensor(1,xn,1,yn,1,zn); /*** USE THIS FOR 3D FFT ***/

  i0 = 1 + my_rint(t0 / tscale);
  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;
    if ((((i-1)%tsamp) == 0) && (t >= t0) && (t < t0+tn)){
      make_frame_sinusoid(data,1,xn,1,yn,1,zn,i,sf,theta,sscale,phase,cx,cy);
      ph = tphase + (float)(i-i0)*tscale*tf * 360.0;
      con = cos(ph*M_PI/180.0);
      multiply_frame(data,1,xn,1,yn,1,zn,i,con);
    }else
      make_frame_flat(data,1,xn,1,yn,1,zn,i,-1.0);
  }
  mid = (float)maxlum/2.0;
  multiply_3d_farray(data,1,xn,1,yn,1,zn,contrast*mid);
  add_const_3d_farray(data,1,xn,1,yn,1,zn,mid);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_STIM_3D_PPL_PLAID                          */
/*                                                                           */
/*  Sum of two gratings.                                                     */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_plaid(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp)
     char *mylogf;
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;  // (deg/pix), (s/frame)
     int tsamp;
{
  int i,j,k;
  float ***data,t,maxlum,t0,tn,t0_2,tn_2,bgamp,bgval,**tframe,blackval;

  //
  //  PREP
  //
  stimx_plaid(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL,NULL);

  t0       = param_getf_exit(ppl,"st0");
  tn       = param_getf_exit(ppl,"stn");
  t0_2     = paramfile_get_float_param_default(ppl,"st0_2",t0);
  tn_2     = paramfile_get_float_param_default(ppl,"stn_2",tn);
  maxlum   = param_getf_exit(ppl,"maxlum");
  bgval    = param_getf_exit(ppl,"bgval");

  blackval = 2.0*bgval - maxlum;

  data = f3tensor(1,xn,1,yn,1,zn);

  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;  // Time (s)
    if (((i-1)%tsamp) == 0){  // This is a video frame
      if (((t >= t0)   && (t < t0+tn)) ||
	  ((t >= t0_2) && (t < t0_2+tn_2))){

	//  'tframe' grating values are centered in [0...maxlum]
	stimx_plaid(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,&tframe,
		    NULL);

	//
	//  Copy this frame in and free it [WYETH - CHANGE THIS]
	//
	for(j=0;j<xn;j++){
	  for(k=0;k<yn;k++){
	    data[j+1][k+1][i] = tframe[j][k];
	  }
	}
	free_2d_farray(tframe,xn);

      }else{
	make_frame_flat(data,1,xn,1,yn,1,zn,i,bgval);
      }
    }else{
      make_frame_flat(data,1,xn,1,yn,1,zn,i,blackval);  // should be blankval
    }
  }

  // Clean up
  stimx_plaid(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL,NULL);
  
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_STIM_3DB_PPL_PLAID                          */
/*                                                                           */
/*****************************************************************************/
void get_stim_3db_ppl_plaid(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,rdl,rdr)
     char *mylogf;
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;  // (deg/pix), (s/frame)
     int tsamp;
     float ****rdl,****rdr;
{
  int i,j,k;
  int dicho,make_left,make_right;
  float **tframl,**tframr;
  float ***dl,***dr,t,maxlum,t0,tn,t0_2,tn_2,bgamp,bgval,blackval;

  //
  //  PREP
  //
  stimx_plaid(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL,NULL);

  t0       = param_getf_exit(ppl,"st0");
  tn       = param_getf_exit(ppl,"stn");
  t0_2     = param_getf_dflt(ppl,"st0_2",t0);
  tn_2     = param_getf_dflt(ppl,"stn_2",tn);
  maxlum   = param_getf_exit(ppl,"maxlum");
  dicho    = param_geti_exit(ppl,"dichoptic");
  bgval    = param_getf_exit(ppl,"bgval");

  blackval = 2.0*bgval - maxlum;

  make_left = make_right = 0;
  if ((dicho == 0) || (dicho == 3))
    make_left = 1;
  if ((dicho == 1) || (dicho == 3))
    make_right = 1;

  dl = f3tensor(1,xn,1,yn,1,zn);
  dr = f3tensor(1,xn,1,yn,1,zn);

  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;  // Time (s)
    if (((i-1)%tsamp) == 0){  // This is a video frame
      if (((t >= t0)   && (t < t0+tn)) ||
	  ((t >= t0_2) && (t < t0_2+tn_2))){

	stimx_plaid(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,&tframl,
		    &tframr);

	//  Copy this frame in and free it
	if (dicho == 0){
	  for(j=0;j<xn;j++){
	    for(k=0;k<yn;k++){
	      dl[j+1][k+1][i] = tframl[j][k];  // Left eye only
	      dr[j+1][k+1][i] = bgval;
	    }
	  }
	}else if (dicho == 1){
	  for(j=0;j<xn;j++){
	    for(k=0;k<yn;k++){
	      dl[j+1][k+1][i] = bgval;
	      dr[j+1][k+1][i] = tframl[j][k];  // Right eye only
	    }
	  }
	}else if ((dicho == 2) || (dicho == 4) || (dicho == 5)){
	  for(j=0;j<xn;j++){
	    for(k=0;k<yn;k++){
	      dl[j+1][k+1][i] = tframl[j][k];  // Dichoptic split
	      dr[j+1][k+1][i] = tframr[j][k];
	    }
	  }
	  free_2d_farray(tframr,xn);
	}else if (dicho == 3){
	  for(j=0;j<xn;j++){
	    for(k=0;k<yn;k++){
	      dl[j+1][k+1][i] = tframl[j][k];  // Same in both eyes
	      dr[j+1][k+1][i] = tframl[j][k];  // Same in both eyes
	    }
	  }
	}else{
	  printf("  *** dicho = %d\n",dicho);
	  exit_error("GET_STIM_3DB_PPL_PLAID","Unknown value of 'dichoptic'");
	}
	free_2d_farray(tframl,xn);

      }else{
	make_frame_flat(dl,1,xn,1,yn,1,zn,i,bgval);
	make_frame_flat(dr,1,xn,1,yn,1,zn,i,bgval);
      }
    }else{
      make_frame_flat(dl,1,xn,1,yn,1,zn,i,blackval);  // should be blankval
      make_frame_flat(dr,1,xn,1,yn,1,zn,i,blackval);  // should be blankval
    }
  }
  
  // Clean up
  stimx_plaid(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL,NULL);

  *rdl = dl;
  *rdr = dr;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_STIM_3DB_PPL_GABOR_GRID                       */
/*                                                                           */
/*****************************************************************************/
void get_stim_3db_ppl_gabor_grid(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,
				 rdl,rdr)
     char *mylogf;
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;  // (deg/pix), (s/frame)
     int tsamp;
     float ****rdl,****rdr;
{
  int i,j,k;
  int make_left,make_right;
  float **tframl,**tframr;
  float ***dl,***dr,t,maxlum,t0,tn,bgval,blackval;

  //
  //  PREP
  //
  stimx_gabor_grid(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL,NULL);

  t0       = param_getf_exit(ppl,"st0");
  tn       = param_getf_exit(ppl,"stn");
  maxlum   = param_getf_exit(ppl,"maxlum");
  bgval    = param_getf_exit(ppl,"bgval");

  blackval = 2.0*bgval - maxlum;

  make_left = make_right = 1;

  dl = f3tensor(1,xn,1,yn,1,zn);
  dr = f3tensor(1,xn,1,yn,1,zn);

  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;  // Time (s)
    if (((i-1)%tsamp) == 0){  // This is a video frame
      if ((t >= t0)   && (t < t0+tn)){

	stimx_gabor_grid(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,&tframl,
			 &tframr);

	for(j=0;j<xn;j++){
	  for(k=0;k<yn;k++){
	    dl[j+1][k+1][i] = tframl[j][k];  // Dichoptic split
	    dr[j+1][k+1][i] = tframr[j][k];
	  }
	}
	free_2d_farray(tframr,xn);
	free_2d_farray(tframl,xn);

      }else{
	make_frame_flat(dl,1,xn,1,yn,1,zn,i,bgval);
	make_frame_flat(dr,1,xn,1,yn,1,zn,i,bgval);
      }
    }else{
      make_frame_flat(dl,1,xn,1,yn,1,zn,i,blackval);  // should be blankval
      make_frame_flat(dr,1,xn,1,yn,1,zn,i,blackval);  // should be blankval
    }
  }

  // Clean up
  stimx_gabor_grid(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL,NULL);

  *rdl = dl;
  *rdr = dr;
}
/**************************************-**************************************/
/*                                                                           */
/*                        GET_STIM_3DB_PPL_GLOBMO_PATCH                      */
/*                                                                           */
/*****************************************************************************/
void get_stim_3db_ppl_globmo_patch(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,
				 rdl,rdr)
     char *mylogf;
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;  // (deg/pix), (s/frame)
     int tsamp;
     float ****rdl,****rdr;
{
  int i,j,k;
  int make_left,make_right;
  float **tframl,**tframr;
  float ***dl,***dr,t,maxlum,t0,tn,bgval,blackval;
  char *stimform;

  //
  //  PREP
  //
  stimx_globmo_patch(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL,NULL);

  t0       = param_getf_exit(ppl,"st0");
  tn       = param_getf_exit(ppl,"stn");
  maxlum   = param_getf_exit(ppl,"maxlum");
  bgval    = param_getf_exit(ppl,"bgval");
  stimform = param_getc_dflt(ppl,"stim_form","3d");


  blackval = 2.0*bgval - maxlum;

  make_left = make_right = 1;

  dl = f3tensor(1,xn,1,yn,1,zn);
  if (strcmp(stimform,"3d_b")==0)
    dr = f3tensor(1,xn,1,yn,1,zn);

  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;  // Time (s)
    if (((i-1)%tsamp) == 0){  // This is a video frame
      if ((t >= t0)   && (t < t0+tn)){

	stimx_globmo_patch(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,
			   &tframl,&tframr);

	if (strcmp(stimform,"3d_b")==0){
	  for(j=0;j<xn;j++){
	    for(k=0;k<yn;k++){
	      dl[j+1][k+1][i] = tframl[j][k];
	      dr[j+1][k+1][i] = tframr[j][k];
	    }
	  }
	}else{
	  for(j=0;j<xn;j++){
	    for(k=0;k<yn;k++){
	      dl[j+1][k+1][i] = tframl[j][k];
	    }
	  }
	}

	free_2d_farray(tframl,xn);
	if (strcmp(stimform,"3d_b")==0)
	  free_2d_farray(tframr,xn);

      }else{
	make_frame_flat(dl,1,xn,1,yn,1,zn,i,bgval);
	if (strcmp(stimform,"3d_b")==0)
	  make_frame_flat(dr,1,xn,1,yn,1,zn,i,bgval);
      }
    }else{
      make_frame_flat(dl,1,xn,1,yn,1,zn,i,blackval);  // should be blankval
      if (strcmp(stimform,"3d_b")==0)
	make_frame_flat(dr,1,xn,1,yn,1,zn,i,blackval);  // should be blankval
    }
  }

  // Clean up
  stimx_globmo_patch(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL,NULL);

  *rdl = dl;
  if (strcmp(stimform,"3d_b")==0)
    *rdr = dr;

  myfree(stimform);
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_STIM_3D_PPL_DMMASK                          */
/*                                                                           */
/*  Sum of two randomly moving gratings.                                     */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_dmmask(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i;
  int seed1,seed2,distrib,ord,dt,seqn,cnt,si,aptype,overcount,nstim,dirflag;
  float ***data1,***data2,***data,ph1,ph2,mid,sz;
  float sf1,sf2,etf1p,etf1a,etf2p,etf2a,phase1,phase2,theta1,theta2,con1,con2;
  float size1,size2,maxlum,t0,tn,t,bgval,bgamp,cx,cy,xc,yc;
  float v1a,v1p,v2a,v2p,*seq1,*seq2,fps;

  // WYETH - looks like we need param for bg here

  nstim = 1;  // Assume one stimulus, unless we find params for two

  cx       = param_getf_dflt(ppl,"cx",0.0);
  cy       = param_getf_dflt(ppl,"cy",0.0);


  if (param_test(ppl,"sf2")){  // WYETH - OLD WAY - discourage

    printf("  *** OLD WAY - Please convert to new param names:\n");
    printf("  sf2    -->  s2_sf\n");
    printf("  etf1a  -->  s1_etfa\n");
    printf("         ...\n");

    etf1p    = param_getf_exit(ppl,"etf1p");
    etf1a    = param_getf_exit(ppl,"etf1a");
    etf2p    = param_getf_exit(ppl,"etf2p");
    etf2a    = param_getf_exit(ppl,"etf2a");
    sf1      = param_getf_exit(ppl,"sf");
    sf2      = param_getf_exit(ppl,"sf2");
    phase1   = param_getf_exit(ppl,"phase");
    phase2   = param_getf_dflt(ppl,"phase2",0.0);
    theta1   = param_getf_exit(ppl,"direction");
    con1     = param_getf_exit(ppl,"contrast");
    con2     = param_getf_exit(ppl,"contrast2");
    size1    = param_getf_exit(ppl,"size");
    size2    = param_getf_dflt(ppl,"size2",size1);
    seed1    = param_geti_exit(ppl,"seed1");
    seed2    = param_geti_exit(ppl,"seed2");

    nstim = 2;

  }else{
    if (param_test(ppl,"etf")){
      etf1p  = param_getf_exit(ppl,"etf");
      etf1a  = etf1p;
      etf2p  = etf1p;
      etf2a  = etf1a;
    }else if (param_test(ppl,"s1_etf")){
      param_getf_s2_exit(ppl,"etf",&etf1p,&etf2p,&nstim);
      etf1a = etf1p;
      etf2a = etf2p;
    }else{
      param_getf_s2_exit(ppl,"etfp",&etf1p,&etf2p,&nstim);
      param_getf_s2_exit(ppl,"etfa",&etf1a,&etf2a,&nstim);
    }
    param_getf_s2_exit(ppl,"sf",&sf1,&sf2,&nstim);
    param_getf_s2_exit(ppl,"phase",&phase1,&phase2,&nstim);
    dirflag = 1;  // Will be set to 2 in next call, if direction2 defined
    param_getf_s2_exit(ppl,"direction",&theta1,&theta2,&dirflag);
    //theta1   = param_getf_exit(ppl,"direction");
    param_getf_s2_exit(ppl,"contrast",&con1,&con2,&nstim);
    param_getf_s2_exit(ppl,"size",&size1,&size2,&nstim);
    param_geti_s2_exit(ppl,"seed",&seed1,&seed2,&nstim);
  }

  /*
    printf("sf1,2    = %f  %f\n",sf1,sf2);
    printf("etf1a,p  = %f  %f\n",etf1a,etf1p);
    printf("etf2a,p  = %f  %f\n",etf2a,etf2p);
    printf("con1,2   = %f  %f\n",con1,con2);
    printf("seed1,2  = %d  %d\n",seed1,seed2);
  */


  maxlum   = param_getf_dflt(ppl,"maxlum",1.0);
  aptype   = param_geti_exit(ppl,"aptype");
  bgval    = param_getf_exit(ppl,"bgval");
  bgamp = bgval * 2.0 - 1.0; // Convert 0..1 --> -1..1

  if (dirflag == 1)  // "direction2" was not specified, assume orthogonal
    theta2 = theta1 + 90.0;

  //
  // Random Sequences
  //
  t0       = param_getf_exit(ppl,"st0");
  tn       = param_getf_exit(ppl,"stn");
  distrib  = param_geti_exit(ppl,"distrib");
  dt       = param_geti_exit(ppl,"dt");

  if (distrib == 7){
    ord = param_geti_exit(ppl,"mseq_ord");
  }else
    ord = -1;

  // Set float values for random sequences
  fps = 1.0/(tsamp*tscale); // Frames/s
  v1p =  etf1p/fps * 360.0;  // step size in degrees of spatial phase
  v1a = -etf1a/fps * 360.0;
  v2p =  etf2p/fps * 360.0;
  v2a = -etf2a/fps * 360.0;

  stm_get_ran_seq(NULL,distrib,ord,seed1,fps,tn,dt,v1p,v1a,0.0,&seq1,&seqn);
  if (nstim == 2)
    stm_get_ran_seq(NULL,distrib,ord,seed2,fps,tn,dt,v2p,v2a,0.0,&seq2,&seqn);

  data1 = f3tensor(1,xn,1,yn,1,zn);
  if (nstim == 2)
    data2 = f3tensor(1,xn,1,yn,1,zn);

  xc = (float)(xn-1)/2.0 + cx/sscale;  // Center w.r.t. 0..xn-1 (pix)
  yc = (float)(yn-1)/2.0 + cy/sscale;  // Center w.r.t. 0..yn-1 (pix)

  ph1 = phase1;
  ph2 = phase2;
  cnt = 0;       // How many times the 'si' element was used
  si = 0;   // Index into the random sequences
  overcount = 0;
  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;

    if (((i-1)%tsamp) == 0){
      
      if ((t >= t0) && (t < t0+tn) && (si < seqn)){
	ph1 += seq1[si];
	make_frame_sinusoid(data1,1,xn,1,yn,1,zn,i,sf1,theta1,sscale,ph1,
			    cx,cy);
	if (nstim == 2){
	  ph2 += seq2[si];
	  make_frame_sinusoid(data2,1,xn,1,yn,1,zn,i,sf2,theta2,sscale,ph2,
			      cx,cy);
	}
	
	cnt += 1;
	if (cnt == dt){
	  si += 1;  // Move on to the next random element
	  cnt = 0;
	}
      }else{
	make_frame_flat(data1,1,xn,1,yn,1,zn,i,bgamp);
	if (nstim == 2)
	  make_frame_flat(data2,1,xn,1,yn,1,zn,i,bgamp);
	if (si >= seqn)
	  overcount += 1;
      }
    }else{
      make_frame_flat(data1,1,xn,1,yn,1,zn,i,-1.0);
      if (nstim == 2)
	make_frame_flat(data2,1,xn,1,yn,1,zn,i,-1.0);
    }
    sz = size1/sscale;
    frame_apply_aperture(data1,1,xn,1,yn,1,zn,i,aptype,xc,yc,sz,sz,0.0,bgamp);
    if (nstim == 2){
      sz = size2/sscale;
      frame_apply_aperture(data2,1,xn,1,yn,1,zn,i,aptype,xc,yc,sz,sz,0.0,bgamp);
    }
  }

  // WYETH - WE WOULD NEED mylogf TO PUT IN THESE LINES
  /*
  if (overcount > 0){
    printf("  *** NOTE: simulation 'tn' is longer than random sequence ");
    printf(" by %d time units\n",overcount);
  }
  if (si < seqn){
  printf("  *** NOTE: %d elements of the random stimulus sequence",seqn-si);
    printf(" were unused\n");
  }
  */


  // Stimulus values are now in -1...1

  mid = (float)maxlum/2.0;
  multiply_3d_farray(data1,1,xn,1,yn,1,zn,con1*mid);

  if (nstim == 2){
    multiply_3d_farray(data2,1,xn,1,yn,1,zn,con2*mid);
    data = add_3d_tensors(data1,data2,xn,yn,zn);
    // Values still centered on 0, but could exceed -mid...mid

    // Restrict stimulus values to -mid...mid
    stim_util_limit_3d(data,1,xn,1,yn,1,zn,-mid,mid);
  }

  if (nstim == 1)
    data = data1;

  add_const_3d_farray(data,1,xn,1,yn,1,zn,mid);
  // Final values should be in 0..maxlum

  if (nstim == 2){
    free_f3tensor(data1,1,xn,1,yn,1,zn);
    free_f3tensor(data2,1,xn,1,yn,1,zn);
  }
  
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_STIM_3D_PPL_DIRMOD                          */
/*                                                                           */
/*  Direction modulation - this includes the old 'dmmask'.                   */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_dirmod(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i,j,k;
  float ***data,t,t0,tn,bgval,**tframe;

  //
  //  PREP
  //
  stimx_dirmod(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL,NULL);

  t0    = param_getf_exit(ppl,"st0");
  tn    = param_getf_exit(ppl,"stn");
  bgval = param_getf_exit(ppl,"bgval");

  data = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;

    if (((i-1)%tsamp) == 0){

      if ((t >= t0) && (t < t0+tn)){

	stimx_dirmod(ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,&tframe,NULL);
	
	//
	//  Copy this frame in and free it [WYETH - CHANGE THIS]
	//
	for(j=0;j<xn;j++){
	  for(k=0;k<yn;k++){
	    data[j+1][k+1][i] = tframe[j][k];
	  }
	}
	free_2d_farray(tframe,xn);

      }else{
	make_frame_flat(data,1,xn,1,yn,1,zn,i,bgval);
      }
    }else{
      make_frame_flat(data,1,xn,1,yn,1,zn,i,0.0);  // should be blankval
    }
  }

  stimx_dirmod(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL,NULL);  // Clean up

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_STIM_3DB_PPL_DIRMOD                         */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
void get_stim_3db_ppl_dirmod(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,rdl,rdr)
     char *mylogf;
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
     float ****rdl,****rdr;
{
  int i,j,k;
  float ***dl,***dr,t,t0,tn,bgval,**tframl,**tframr;

  //
  //  PREP
  //
  stimx_dirmod(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL,NULL);

  t0    = param_getf_exit(ppl,"st0");
  tn    = param_getf_exit(ppl,"stn");
  bgval = param_getf_exit(ppl,"bgval");

  dl = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT
  dr = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;

    if (((i-1)%tsamp) == 0){

      if ((t >= t0) && (t < t0+tn)){

	stimx_dirmod(ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,&tframl,&tframr);

	//
	//  Copy this frame in and free it [WYETH - CHANGE THIS]
	//
	for(j=0;j<xn;j++){
	  for(k=0;k<yn;k++){
	    dl[j+1][k+1][i] = tframl[j][k];
	    dr[j+1][k+1][i] = tframr[j][k];
	  }
	}
	free_2d_farray(tframl,xn);
	free_2d_farray(tframr,xn);

      }else{
	make_frame_flat(dl,1,xn,1,yn,1,zn,i,bgval);
	make_frame_flat(dr,1,xn,1,yn,1,zn,i,bgval);
      }
    }else{
      make_frame_flat(dl,1,xn,1,yn,1,zn,i,0.0);  // should be blankval
      make_frame_flat(dr,1,xn,1,yn,1,zn,i,0.0);  // should be blankval
    }
  }

  stimx_dirmod(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL,NULL);  // Clean up

  *rdl = dl;
  *rdr = dr;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_STIM_3D_PPL_RANSINPHASE                       */
/*                                                                           */
/*  Also used for ransinori and ransindir.                                   */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*    sf - spatial frequency (cyc/deg)                                       */
/*    tf - temporal frequency (cyc/sec)                                      */
/*    phase - phase offset (deg) or phase jump for sflag=2                   */
/*    theta - orientation of grating (deg)                                   */
/*    sflag                                                                  */
/*      0 - ransinphase - random counterphase                                */
/*      1 - ransinori - random orthogonal                                    */
/*      2 - ransindir - random drift (start at ph=0, change by +-phase)      */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*  - 'ransindir' is probably the same as 'dirmod' stimtype.  Do we need     */
/*    it here?                                                               */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_ransinphase(ppl,xn,yn,zn,sscale,tscale,tsamp,sflag)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
     int sflag;  // 0 - ransinphase,  1 - ransinori,  2 - ransindir
{
  int i;
  int ri,cnt,dt,aptype,seed,*rseq,nrseq,nif,antiflag,aflag,distrib,dumpflag;
  int ms_ord,szn,epoch1,epoch2,ecnt,test_flag;
  float *frseq,bgval,gmu,con2,c2f,phase2,dmvar_fact,ffact;
  float ***data,ph,th,mid,sf,phase,theta,contrast,maxlum,size,bgampl,sz;
  float ff,ffrac,phtau,t0,tn,t,ph_anti,cx,cy,xc,yc,blankval,blnk;

  cx       = paramfile_get_float_param_default(ppl,"cx",0.0);
  cy       = paramfile_get_float_param_default(ppl,"cy",0.0);

  maxlum   = param_getf_exit(ppl,"maxlum");
  theta    = param_getf_exit(ppl,"direction");
  sf       = param_getf_exit(ppl,"sf");
  contrast = param_getf_exit(ppl,"contrast");
  size     = param_getf_exit(ppl,"size");
  phase    = param_getf_exit(ppl,"phase");
  ph_anti  = paramfile_get_float_param_default(ppl,"phase_anti",phase);
  t0       = param_getf_exit(ppl,"st0");
  tn       = param_getf_exit(ppl,"stn");
  distrib  =   paramfile_get_int_param_or_exit(ppl,"distrib");
  dt       =   paramfile_get_int_param_or_exit(ppl,"dt");
  seed     =   paramfile_get_int_param_or_exit(ppl,"seed");
  aptype   =   paramfile_get_int_param_or_exit(ppl,"aptype");
  aflag    =   paramfile_get_int_param_or_exit(ppl,"antiflag");
  bgval    = param_getf_exit(ppl,"bgval");
  phtau    = paramfile_get_float_param_default(ppl,"stim_phtau",0.006);
  blankval = paramfile_get_float_param_default(ppl,"stim_blank_val",0.0);

  blnk = 2.0*blankval - 1.0;

  /***
      etf3  = param_getf_exit(ppl,"etf3");
      etf4  = param_getf_exit(ppl,"etf4");
      // Convert to fractional shift of texture map grating
      frac3 = etf3/lss->fps;
      frac4 = etf4/lss->fps;
      ep1n = paramfile_get_int_param_or_exit(ppl,"epoch1");
      ep2n = paramfile_get_int_param_or_exit(ppl,"epoch2");
  ***/

  dmvar_fact = 1.0;  // Scale factor for phase step if DMVAR

  // Alternating epochs, e.g. contrast (for wyphasm)
  epoch1 = paramfile_get_int_param_default(ppl,"epoch1",0);
  if (epoch1 > 0){
    epoch2 = paramfile_get_int_param_or_exit(ppl,"epoch2");

    if (paramfile_test_param(ppl,"phase2")){ // DMVAR
      phase2 = param_getf_exit(ppl,"phase2");
      dmvar_fact = phase2 / phase;
      c2f = 1.0;  // Contrast scale factor for epoch2
    }else{
      con2 = paramfile_get_float_param_default(ppl,"contrast2",-1.0);
      if (con2 > contrast){
	exit_error("GET_STIM_3D_PPL_RANSINPHASE","contrast2 > contrast");
      }
      c2f = con2/contrast;  // Scale factor for epoch2
    }
  }


  dumpflag = paramfile_get_int_param_default(ppl,"dumpflag",0);

  if (phtau > 0.0)
    ffrac = exp(-tscale/phtau); // for phosphor decay w/ time const phtau
  else
    ffrac = 0.0;
  ff = ffrac;   // Ampl for next empty frame

  bgampl = bgval * 2.0 - 1.0; // Convert 0..1 --> -1..1

  if (0) printf("%f %d\n",size,aptype);

  data = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  xc = (float)(xn-1)/2.0 + cx/sscale;  // Center w.r.t. 0..xn-1 (pix)
  yc = (float)(yn-1)/2.0 + cy/sscale;  // Center w.r.t. 0..yn-1 (pix)

  if ((sflag == 0)||(sflag == 1))
    ph = phase;
  else{
    ph = paramfile_get_float_param_default(ppl,"ph0",0.0); // dm, initial phase
  }
  th = theta;

  /*** The calculation below was changed so that the number of random
       draws to be requested, 'nrseq', is determined by the stimulus
       'tn' rather than 'zn'.  If 'tn' is larger than 'zn', then use
       'zn'  ***/

  szn = (int)(tn/tscale);  /* Stim length in sampling units (tn is sec) */
  if (szn > zn)
    szn = zn;

  nif = (int)((szn + tsamp - 1)/tsamp); /* Number of image frames */
  nrseq = 1 + (int)((nif + dt -1)/ dt); /* Number of random draws */
                     /*  WYETH 1+ added because 'Exceeded rand...' below */ 

  /*printf("stim_util WYETH HERE   nrseq = %d\n",nrseq);*/

  // Usually, 0-uniform, 1-gauss, 2-bin, 3-tern ...
  rseq = NULL;
  frseq = NULL;
  if (distrib == 1){
    frseq = myrand_get_std_gauss_seq(nrseq,seed);
    gmu = param_getf_exit(ppl,"mu");
  }else if (distrib == 2){
    rseq = myrand_get_std_bin_seq(nrseq,seed);
  }else if (distrib == 7){
    ms_ord = paramfile_get_int_param_or_exit(ppl,"mu");
    rseq = myrand_get_std_mseq(nrseq,seed,ms_ord);
  }else
    exit_error("GET_STIM_3D_PPL_RANSINPHASE","Unknown distrib value");

  if ((sflag == 2) && (distrib != 1)){ // WYETH - new to allow test patterns
    frseq = i2farray(rseq,nrseq);
    for(i=0;i<nrseq;i++){
      if (rseq[i])
	frseq[i] = phase;
      else
	frseq[i] = -ph_anti;
    }
    myfree(rseq);
    rseq = NULL;
  }

  //
  // If 'dm_test', add in test sequences
  //
  test_flag = paramfile_get_int_param_default(ppl,"test_pat",-1);
  if (test_flag >= 0){
    int t_pat,n_anti,t_dur,t_int;
    float t_etf,fps;

    fps = 1.0/(tsamp*tscale); // Frames/s

    //printf("  Frames/sec = %f\n",fps);
    //printf("  nrseq = %d\n",nrseq);
    //printf("  phase = %f\n",phase);

    t_etf  = param_getf_exit(ppl,"test_etf");
    t_pat  = paramfile_get_int_param_or_exit(ppl,"test_pat");
    n_anti = paramfile_get_int_param_default(ppl,"test_nanti",3);
    t_dur  = paramfile_get_int_param_or_exit(ppl,"test_dur");
    t_int  = paramfile_get_int_param_or_exit(ppl,"test_int");

    // sflag = 1 to use proper sign for test pats, degflag = 1 for degrees 
    stm_dirmod_test_seq(NULL,t_dur,t_int,t_etf,t_pat,n_anti,fps,1,1,
			&frseq,&nrseq);
  }


  if (dumpflag == 1){
    printf("  Dumping, nrseq = %d\n",nrseq);
    if (rseq != NULL)
      append_iarray_plot("zzz.STIM.pl","rseq",rseq,nrseq,1);
    if (frseq != NULL)
      append_farray_plot("zzz.STIM.pl","fseq",frseq,nrseq,1);
    printf("Writing random time sequence to 'zzz.STIM.pl'\n");
  }
  
  antiflag = 0;

  ffact = 1.0;
  cnt = dt;
  ri = 0;
  ecnt = 0;   // Epoch timer
  for(i=1;i<=zn;i++){
    if (((i-1)%tsamp) == 0){
      t = (float)(i-1)*tscale;     // Time (s), float
      if ((t >= t0) && (t < (t0+tn))){
	if (cnt==dt){
	  if (ri >= nrseq){
	    printf("ri=%d  nrseq=%d  t,t0,tn= %f %f %f\n",ri,nrseq,t,t0,tn);
	    exit_error("GET_STIM_3D_PPL_RANSINPHASE","Exceeded rand sequence");
	  }
          if (rseq != NULL){
	    if (rseq[ri])
	      antiflag = 0;
	    else
	      antiflag = 1;
	  }
	  
	  if (sflag == 0){          // Random counterphase
	    if (rseq[ri])
	      ph = phase;
	    else
	      ph = phase + 180.0;
	  }else if (sflag == 1){    // Random orthogonal orientation
	    if (rseq[ri])
	      th = theta;
	    else
	      th = theta + 90.0;
	  }else if (sflag == 2){    // Random direction

	    if (dmvar_fact != 1.0){
	      if (ecnt < epoch1){
		ffact = 1.0;
		//printf("e1     %f\n",ffact);
	      }else{
		ffact = dmvar_fact;
		//printf("ffact  %f\n",ffact);
	      }
	    }
	      
	    if (distrib == 1){
	      ph += gmu + phase* ffact * frseq[ri]; // Gaussian
	    }else{
	      ph += ffact * frseq[ri];              // Binary
	    }
	  }else
	    exit_error("get_stim_3d_ransinphase","Unknown sflag value");
	  ri += 1;
	  cnt = 0;
	}else{
	  if (sflag == 2){ // Keep the grating moving
	    if (distrib == 1){
	      ph += (ffact * frseq[ri-1] * phase); // Gaussian
	    }else{
	      ph += ffact * frseq[ri-1];           // Binary
	    }
	  }
	}
	if ((aflag==0) && antiflag) // Use gray (Null) instead of antipref
	  make_frame_flat(data,1,xn,1,yn,1,zn,i,0.0);
	else{
	  if (epoch1 <= 0){ // WYETH, '=' added here to fix bug
	    make_frame_sinusoid(data,1,xn,1,yn,1,zn,i,sf,th,sscale,ph,cx,cy);
	  }else{
	    if (ecnt < epoch1){
	      make_frame_sinusoid(data,1,xn,1,yn,1,zn,i,sf,th,sscale,ph,cx,cy);
	    }else{
	      make_frame_sinusoid_scale(data,1,xn,1,yn,1,zn,i,sf,th,sscale,ph,
					cx,cy,c2f,0);
	    }
	  }
	}
	cnt += 1;
      }else{
	make_frame_flat(data,1,xn,1,yn,1,zn,i,bgampl); // Show BG value
      }

      ff = ffrac;
    }else{
      if ((ff > 0.0)&&(i>1))
	//make_frame_scale_frame_zmid(data,1,xn,1,yn,1,zn,i,i-1,ff,-1.0);
	make_frame_scale_frame_zmid(data,1,xn,1,yn,1,zn,i,i-1,ff,blnk);
      else
	//make_frame_flat(data,1,xn,1,yn,1,zn,i,-1.0);
	make_frame_flat(data,1,xn,1,yn,1,zn,i,blnk);
      ff *= ffrac;
    }

    sz = size/sscale;
    frame_apply_aperture(data,1,xn,1,yn,1,zn,i,aptype,xc,yc,sz,sz,0.0,bgampl);

    //printf("  frame %d  ri %d  cnt %d dt %d ph %f \n",i,ri,cnt,dt,ph);
    ecnt += 1;
    if (ecnt >= (epoch1+epoch2))
      ecnt = 0;
  }
  mid = (float)maxlum/2.0;
  multiply_3d_farray(data,1,xn,1,yn,1,zn,contrast*mid);
  add_const_3d_farray(data,1,xn,1,yn,1,zn,mid);

  /* WYETH - this doesn't work since I implemented st0,stn
  if (ri != nrseq)
  exit_error("GET_STIM_3D_PPL_RANSINPHASE","Excess random sequence");*/

  if (frseq != NULL)
    myfree(frseq);
  if (rseq != NULL)
    myfree(rseq);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_STIM_3D_PPL_RANSTEP                         */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_ranstep(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i;
  int ri,cnt,aptype,dumpflag,jn,jzero,dt,seed,seqn,meanzero;
  float sf,fps,*seqf,phase,jsz,st0,stn,bgval,t,cx,cy,xc,yc;
  float ***data,ph,th,mid,theta,contrast,maxlum,size,bgampl,sz;
  float blankval,blnk;

  cx       = paramfile_get_float_param_default(ppl,"cx",0.0);
  cy       = paramfile_get_float_param_default(ppl,"cy",0.0);

  //
  //  Params from 'ss_ranstep'
  //
  sf       = param_getf_exit(ppl,"sf");
  size     = param_getf_exit(ppl,"size");
  theta    = param_getf_exit(ppl,"direction");
  bgval    = param_getf_exit(ppl,"bgval");
  contrast = param_getf_exit(ppl,"contrast");

  // Prepare the random sequence of movements
  st0      = param_getf_exit(ppl,"st0");
  stn      = param_getf_exit(ppl,"stn");
  dt       = paramfile_get_int_param_or_exit(ppl,"dt");
  seed     = paramfile_get_int_param_or_exit(ppl,"seed");
  //varflag  = paramfile_get_int_param_or_exit(ppl,"var");
  jn       = paramfile_get_int_param_or_exit(ppl,"step_n");
  jsz      = param_getf_exit(ppl,"step_size");
  jzero    = paramfile_get_int_param_or_exit(ppl,"step_zero");
  //  phase is fraction of cycle in 'ranstep'
  phase    = param_getf_exit(ppl,"phase");

  maxlum   = param_getf_exit(ppl,"maxlum");
  aptype   =   paramfile_get_int_param_or_exit(ppl,"aptype");
  dumpflag =   paramfile_get_int_param_default(ppl,"dumpflag",0);

  meanzero = paramfile_get_int_param_default(ppl,"stim_mean_zero",0);
  blankval = paramfile_get_float_param_default(ppl,"stim_blank_val",0.0);

  blnk = 2.0*blankval - 1.0;

  fps = 1.0/(tsamp*tscale); // Frames/s

  // Jump sequence in terms of fraction of cycle to jump
  seqf = stm_ranstep_get_stim(NULL,sf,jsz,jzero,jn,seed,stn,fps,dt,0,&seqn);

  /***
  printf("seqn = %d\n",seqn);
  printf("stn = %.4f\n",stn);
  printf("fps = %.4f\n",fps);
  for(i=0;i<20;i++)
  printf("seqf[%d] = %8.4f degr\n",i,seqf[i]*360.0);***/


  bgampl = bgval * 2.0 - 1.0; // Convert 0..1 --> -1..1

  data = f3tensor(1,xn,1,yn,1,zn);  // Use this for 3D FFT

  xc = (float)(xn-1)/2.0 + cx/sscale;  // Center w.r.t. 0..xn-1 (pix)
  yc = (float)(yn-1)/2.0 + cy/sscale;  // Center w.r.t. 0..yn-1 (pix)

  ph = phase*360.0;  // Convert from fraction of cycle to degr
  th = theta;

  /*** The calculation below was changed so that the number of random
       draws to be requested, 'nrseq', is determined by the stimulus
       'tn' rather than 'zn'.  If 'tn' is larger than 'zn', then use
       'zn'  ***/

  if (dumpflag == 1){
    printf("  Dumping, nrseq = %d\n",seqn);
    append_farray_plot("zzz.STIM.pl","seqf",seqf,seqn,1);
    printf("Writing random time sequence to 'zzz.STIM.pl'\n");
  }

  //printf("WYETH HERE --------------- tsamp = %d\n",tsamp);
  
  cnt = dt;
  ri = 0;
  for(i=1;i<=zn;i++){
    if (((i-1)%tsamp) == 0){
      t = (float)(i-1)*tscale;     // Time (s), float
      //if ((t >= st0) && (t < (st0+stn))){
      if ((t >= st0) && (t < (st0+stn)) && (ri < seqn)){
	if (cnt==dt){
	  if (ri >= seqn){
	    printf("ri=%d  nrseq=%d  t,t0,tn= %f %f %f\n",ri,seqn,t,st0,stn);
	    exit_error("GET_STIM_3D_PPL_RANSINPHASE","Exceeded rand sequence");
	    // WYETH THIS CAN'T HAPPEN NOW (see 'if' above)
	  }
	  
	  ph -= seqf[ri]*360.0;  // Convert jump from frac to degrees

	  ri += 1;
	  cnt = 0;
	}else{
	  ph -= seqf[ri-1]*360.0;  // Use previous jump size 'dt' times
	}

	//if (ri < 20)
	//printf("ph(used) = %8.4f degr\n",ph);

	make_frame_sinusoid(data,1,xn,1,yn,1,zn,i,sf,th,sscale,ph,cx,cy);

	cnt += 1;
      }else{
	make_frame_flat(data,1,xn,1,yn,1,zn,i,bgampl); // Show BG value
      }

    }else{
      //make_frame_flat(data,1,xn,1,yn,1,zn,i,-1.0);
      make_frame_flat(data,1,xn,1,yn,1,zn,i,blnk);
    }

    sz = size/sscale;
    frame_apply_aperture(data,1,xn,1,yn,1,zn,i,aptype,xc,yc,sz,sz,0.0,bgampl);

    //printf("  frame %d  ri %d  cnt %d dt %d ph %f \n",i,ri,cnt,dt,ph);
  }

  if (meanzero == 0){
    mid = (float)maxlum/2.0;
    multiply_3d_farray(data,1,xn,1,yn,1,zn,contrast*mid);
    add_const_3d_farray(data,1,xn,1,yn,1,zn,mid);
  }else{
    if (contrast*maxlum != 1.0)
      multiply_3d_farray(data,1,xn,1,yn,1,zn,contrast*maxlum);
  }


  myfree(seqf);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_STIM_3D_PPL_TEST_SHAPE                        */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Returned array dimensions start at 1.  To free returned storage:       */
/*    free_f3tensor(data,1,xn,1,yn,1,zn).                                    */
/*  - Do not reseed here, the caller must initialize the seed                */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_test_shape(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i,j,k;
  float ***data,mid,t,contrast,maxlum,t0,tn,bgval,bgamp,**tframe;
  float sbot,sleft,sright,cx,cy,ampl;

  //
  //  WYETH - developed this to make a triangle to test phase randomizing
  //  for shapes like Anitha's
  //

  t0        = param_getf_exit(ppl,"st0");
  tn        = param_getf_exit(ppl,"stn");

  cx        = param_getf_dflt(ppl,"tr1_cx",0.0);
  cy        = param_getf_dflt(ppl,"tr1_cy",0.0);
  sbot      = param_getf_exit(ppl,"tr1_bottom");
  sleft     = param_getf_exit(ppl,"tr1_left");
  sright    = param_getf_exit(ppl,"tr1_right");
  ampl      = param_getf_exit(ppl,"tr1_ampl");

  bgval     = param_getf_exit(ppl,"bgval");
  contrast  = param_getf_exit(ppl,"contrast");
  maxlum    = param_getf_exit(ppl,"maxlum");

  data = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;  // Time (s)
    if (((i-1)%tsamp) == 0){  // Is this a video frame
      if ((t >= t0) && (t < t0+tn)){

	make_frame_flat(data,1,xn,1,yn,1,zn,i,bgval);

	make_frame_add_triangle(data,1,xn,1,yn,1,zn,i,sbot,sleft,sright,
				cx,cy,sscale,ampl,bgval);

      }else{
	make_frame_flat(data,1,xn,1,yn,1,zn,i,bgval);
      }
    }else{
      make_frame_flat(data,1,xn,1,yn,1,zn,i,-1.0); // typ. black, -1.0
    }
  }

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_STIM_3D_PPL_ADAMK                          */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_adamk(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i;
  int ri,cnt,dt,aptype,seed,nrseq,nif,aflag,dumpflag;
  int szn,tau,std1,std2,period,offset,floorp;
  float *frseq,bgval,scalep;
  float ***data,ph,th,mid,sf,phase,theta,contrast,maxlum,size,bgampl,sz;
  float ff,ffrac,phtau,t0,tn,t,cx,cy,xc,yc;

  cx       = paramfile_get_float_param_default(ppl,"cx",0.0);
  cy       = paramfile_get_float_param_default(ppl,"cy",0.0);

  maxlum   = param_getf_exit(ppl,"maxlum");
  theta    = param_getf_exit(ppl,"direction");
  sf       = param_getf_exit(ppl,"sf");
  contrast = param_getf_exit(ppl,"contrast");
  size     = param_getf_exit(ppl,"size");
  phase    = param_getf_exit(ppl,"phase0");
  t0       = param_getf_exit(ppl,"st0");
  tn       = param_getf_exit(ppl,"stn");
  dt       =   paramfile_get_int_param_or_exit(ppl,"dt");
  aptype   =   paramfile_get_int_param_or_exit(ppl,"aptype");
  bgval    = param_getf_exit(ppl,"bgval");
  phtau    = paramfile_get_float_param_default(ppl,"stim_phtau",0.0);

  /*** Params for random motion signal ***/
  tau      =   paramfile_get_int_param_or_exit(ppl,"tau");
  seed     =   paramfile_get_int_param_or_exit(ppl,"seed");
  std1     =   paramfile_get_int_param_or_exit(ppl,"std1");
  std2     =   paramfile_get_int_param_or_exit(ppl,"std1");
  period   =   paramfile_get_int_param_or_exit(ppl,"period");
  offset   =   paramfile_get_int_param_or_exit(ppl,"offset");
  floorp   =   paramfile_get_int_param_or_exit(ppl,"floorp");
  scalep   =   param_getf_exit(ppl,"scalep");

  dumpflag = paramfile_get_int_param_default(ppl,"dumpflag",0);

  if (phtau > 0.0)
    ffrac = exp(-tscale/phtau); /* for phosphor decay w/ time const phtau */
  else
    ffrac = 0.0;
  ff = ffrac;   /* Ampl for next empty frame */

  bgampl = bgval * 2.0 - 1.0; /* Convert 0..1 --> -1..1 */

  data = f3tensor(1,xn,1,yn,1,zn); /*** USE THIS FOR 3D FFT ***/

  xc = (float)(xn-1)/2.0 + cx/sscale;  // Center w.r.t. 0..xn-1 (pix)
  yc = (float)(yn-1)/2.0 + cy/sscale;  // Center w.r.t. 0..yn-1 (pix)

  ph = phase;
  th = theta;

  /*** The calculation below was changed so that the number of random
       draws to be requested, 'nrseq', is determined by the stimulus
       'tn' rather than 'zn'.  If 'tn' is larger than 'zn', then use
       'zn'  ***/

  szn = (int)(tn/tscale);  /* Stim length in sampling units (tn is sec) */
  if (szn > zn)
    szn = zn;

  nif = (int)((szn + tsamp - 1)/tsamp); /* Number of image frames */
  nrseq = 1 + (int)((nif + dt -1)/ dt); /* Number of random draws */
                     /*  WYETH 1+ added because 'Exceeded rand...' below */ 

  /*printf("nrseq = %d  szn %d  tn %f  nif %d  tsamp %d\n",
    nrseq,szn,tn,nif,tsamp);*/

  /*
    printf("tau %d seed %d std1,2, %d %d period %d\n",tau,seed,std1,std2);
    printf("period %d  offset %d   floorsca %f\n",period,offset,floorsca);*/

  frseq = akstim(tau,seed,std1,std2,period,offset,floorp,scalep,nrseq);

  if (dumpflag == 1){
    append_farray_plot("zzz.STIM.pl","fseq",frseq,nrseq,1);
    printf("Writing random time sequence to 'zzz.STIM.pl'\n");
  }

  cnt = dt;
  ri = 0;
  for(i=1;i<=zn;i++){
    if (((i-1)%tsamp) == 0){
      t = (float)(i-1)*tscale;     /* Time (s), float */
      /* if (1){ */
      if ((t >= t0) && (t < (t0+tn))){
	if (cnt==dt){
	  if (ri >= nrseq){
	    printf("ri=%d  nrseq=%d  t,t0,tn= %f %f %f\n",ri,nrseq,t,t0,tn);
	    exit_error("GET_STIM_3D_PPL_RANSINPHASE","Exceeded rand sequence");
	  }
	  
	  ph += frseq[ri]; /* Gaussian */
	  
	  ri += 1;
	  cnt = 0;
	}else{
	  /*** WYETH - I corrected this to be ri-1 here ***/
	  ph += frseq[ri-1]; /* Gaussian */
	}
	make_frame_sinusoid(data,1,xn,1,yn,1,zn,i,sf,th,sscale,ph,cx,cy);
	
	cnt += 1;
      }else{
	/* Use BG value before st0 or after st0+stn */
	make_frame_flat(data,1,xn,1,yn,1,zn,i,bgampl);
      }
      
      ff = ffrac;
    }else{
      if ((ff > 0.0)&&(i>1))
	make_frame_scale_frame_zmid(data,1,xn,1,yn,1,zn,i,i-1,ff,-1.0);
      else
	make_frame_flat(data,1,xn,1,yn,1,zn,i,-1.0);
      ff *= ffrac;
    }

    sz = size/sscale;
    frame_apply_aperture(data,1,xn,1,yn,1,zn,i,aptype,xc,yc,sz,sz,0.0,bgampl);
  }
  mid = (float)maxlum/2.0;
  multiply_3d_farray(data,1,xn,1,yn,1,zn,contrast*mid);
  add_const_3d_farray(data,1,xn,1,yn,1,zn,mid);

  myfree(frseq);
  
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_STIM_3D_PPL_RAN_SIN                       */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_ran_sin(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i,j,k,l;
  int m,ri,cnt,dt,aptype,seed,*rseq,nrseq,nif,szn,ori_n,ph_n,sf_n,blank_n;
  float bgval,t0,tn,t,sf_m,sf_low,ori,phase,cx,cy,xc,yc;
  float ***data,or,ph,mid,sf,contrast,maxlum,size,bgamp,sz;
  float *tab1_or,*tab2_ph,*tab3_sf;

  cx       = paramfile_get_float_param_default(ppl,"cx",0.0);
  cy       = paramfile_get_float_param_default(ppl,"cy",0.0);
  t0       = param_getf_exit(ppl,"st0");
  tn       = param_getf_exit(ppl,"stn");
  maxlum   = param_getf_exit(ppl,"maxlum");
  bgval    = param_getf_exit(ppl,"bgval");
  contrast = param_getf_exit(ppl,"contrast");
  blank_n  =   paramfile_get_int_param_or_exit(ppl,"blank_n");
  ori_n    =   paramfile_get_int_param_or_exit(ppl,"ori_n");
  ori      = paramfile_get_float_param_default(ppl,"ori",0.0);
  ph_n     =   paramfile_get_int_param_or_exit(ppl,"phase_n");
  phase    = paramfile_get_float_param_default(ppl,"phase",0.0);
  sf_n     =   paramfile_get_int_param_or_exit(ppl,"sf_n");
  sf_low   = param_getf_exit(ppl,"sf_low");
  sf_m     = param_getf_exit(ppl,"sf_m");
  seed     =   paramfile_get_int_param_or_exit(ppl,"seed");
  size     = param_getf_exit(ppl,"size");
  dt       =   paramfile_get_int_param_or_exit(ppl,"dt");
  aptype   =   paramfile_get_int_param_or_exit(ppl,"aptype");

  bgamp = bgval * 2.0 - 1.0; /* Convert 0..1 --> -1..1 */

  data = f3tensor(1,xn,1,yn,1,zn); /*** Use this for 3D FFT ***/

  /*** The calculation below was changed so that the number of random
       draws to be requested, 'nrseq', is determined by the stimulus
       'tn' rather than 'zn'.  If 'tn' is larger than 'zn', then use
       'zn'  ***/

  szn = (int)(tn/tscale);  /* Stim length in sampling units (tn is sec) */
  if (szn > zn)
    szn = zn;

  sz = size/sscale;
  xc = (float)(xn-1)/2.0 + cx/sscale;  // Center w.r.t. 0..xn-1 (pix)
  yc = (float)(yn-1)/2.0 + cy/sscale;  // Center w.r.t. 0..yn-1 (pix)


  nif = (int)((szn + tsamp - 1)/tsamp); /* Number of image frames */
  nrseq = 1 + (int)((nif + dt -1)/ dt); /* Number of random draws */
                     /*  WYETH 1+ added because 'Exceeded rand...' below */ 

  m = ori_n * ph_n * sf_n + blank_n;
  tab1_or = (float *)myalloc(m*sizeof(float));
  tab2_ph = (float *)myalloc(m*sizeof(float));
  tab3_sf = (float *)myalloc(m*sizeof(float));
  rseq = myrand_get_std_unif_int_seq(nrseq,seed,m);
  l = 0;
  for(i=0;i<ori_n;i++){
    for(j=0;j<ph_n;j++){
      for(k=0;k<sf_n;k++){
	tab1_or[l] = ori + (float)i * 180.0/(float)ori_n;
	tab2_ph[l] = phase + (float)j * 360.0/(float)ph_n;
	tab3_sf[l] = sf_low * pow(sf_m,(float)k);
	l += 1;
      }
    }
  }
  for(i=0;i<blank_n;i++){ /*** Add blank frames at end ***/
    tab1_or[l] = 0.0;
    tab2_ph[l] = 0.0;
    tab3_sf[l] = 0.0;  /* Indicates blank */
    l += 1;
  }

  cnt = dt;
  ri = 0;
  for(i=1;i<=zn;i++){
    if (((i-1)%tsamp) == 0){
      t = (float)(i-1)*tscale;     /* Time (s), float */
      if ((t >= t0) && (t < (t0+tn))){
	if (cnt == dt){ /* Pick a new stimulus pattern */
	  if (ri >= nrseq){
	    printf("ri=%d  nrseq=%d  t,t0,tn= %f %f %f\n",ri,nrseq,t,t0,tn);
	    exit_error("GET_STIM_3D_PPL_RAN_SIN","Exceeded rand sequence");
	  }
	  k = rseq[ri];     /* Index to a random stimulus from 0..m-1 */
	  ri += 1;
	  or = tab1_or[k];
	  ph = tab2_ph[k];
	  sf = tab3_sf[k];
	  cnt = 0;
	}
	if (sf == 0.0)
	  make_frame_flat(data,1,xn,1,yn,1,zn,i,bgamp);
	else
	  make_frame_sinusoid(data,1,xn,1,yn,1,zn,i,sf,or,sscale,ph,cx,cy);
	
	cnt += 1;
      }else{
	/* Use BG value before st0 or after st0+stn */
	make_frame_flat(data,1,xn,1,yn,1,zn,i,bgamp);
      }
    }else{
      make_frame_flat(data,1,xn,1,yn,1,zn,i,-1.0); /* Black */
    }
    frame_apply_aperture(data,1,xn,1,yn,1,zn,i,aptype,xc,yc,sz,sz,0.0,bgamp);
  }
  mid = (float)maxlum/2.0;
  multiply_3d_farray(data,1,xn,1,yn,1,zn,contrast*mid);
  add_const_3d_farray(data,1,xn,1,yn,1,zn,mid);

  /* WYETH - this doesn't work since I implemented st0,stn
  if (ri != nrseq)
  exit_error("GET_STIM_3D_PPL_RANSINPHASE","Excess random sequence");*/

  myfree(rseq);
  myfree(tab1_or);  myfree(tab2_ph);  myfree(tab3_sf);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_STIM_3D_PPL_BAR                           */
/*                                                                           */
/*  Return a drifting bar.                                                   */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*    sf - spatial frequency (cyc/deg)                                       */
/*    tf - temporal frequency (cyc/sec)                                      */
/*    phase - phase offset (deg)                                             */
/*    theta - orientation of grating (deg)                                   */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_bar(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i,j,k;
  float ***data,bgval,blankval,**tframe;

  //
  //  PREP
  //
  stimx_bar(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL,NULL);

  bgval    = param_getf_exit(ppl,"bgval");
  blankval = param_getf_dflt(ppl,"stim_blank_val",0.0);

  data = f3tensor(1,xn,1,yn,1,zn); // USE THIS FOR 3D FFT

  for(i=1;i<=zn;i++){
    if (((i-1)%tsamp) == 0){

      stimx_bar(ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,&tframe,NULL);

      //
      //  Copy this frame in and free it [WYETH - CHANGE THIS]
      //
      for(j=0;j<xn;j++){
	for(k=0;k<yn;k++){
	  data[j+1][k+1][i] = tframe[j][k];
	}
      }
      free_2d_farray(tframe,xn);

    }else{
      make_frame_flat(data,1,xn,1,yn,1,zn,i,blankval-bgval); // void value
      // This ends up as 'blankval' when 'bgval is added below
    }
  }
  add_const_3d_farray(data,1,xn,1,yn,1,zn,bgval);

  stimx_bar(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL,NULL);  // Clean up

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                        GET_STIM_3D_PPL_GABOR_BRITTEN                      */
/*                                                                           */
/*  Britten and Heuer (1999) stimulus.                                       */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_gabor_britten(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp)
     char *mylogf;                 // Log file
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i,j,k;
  float ***data,bgval,blankval,**tframe;

  //
  //  PREP
  //
  stimx_gabor_britten(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL);

  bgval    = param_getf_exit(ppl,"bgval");
  blankval = param_getf_dflt(ppl,"stim_blank_val",0.0);

  data = f3tensor(1,xn,1,yn,1,zn); // USE THIS FOR 3D FFT

  for(i=1;i<=zn;i++){
    if (((i-1)%tsamp) == 0){

      stimx_gabor_britten(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,
			  &tframe);

      //
      //  Copy this frame in and free it [WYETH - CHANGE THIS]
      //
      for(j=0;j<xn;j++){
	for(k=0;k<yn;k++){
	  data[j+1][k+1][i] = tframe[j][k];
	}
      }
      free_2d_farray(tframe,xn);

    }else{
      make_frame_flat(data,1,xn,1,yn,1,zn,i,blankval-bgval); // void value
      // This ends up as 'blankval' when 'bgval is added below
    }
  }
  multiply_3d_farray(data,1,xn,1,yn,1,zn,0.5);
  add_const_3d_farray(data,1,xn,1,yn,1,zn,bgval);

  // Clean up
  stimx_gabor_britten(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_STIM_3D_PPL_3BARLINK                        */
/*                                                                           */
/*  Three-bar linkage, like that of Nandy et al (2013).                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_3barlink(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i,j,k;
  float ***data,bgval,blankval,**tframe;

  //
  //  PREP
  //
  stimx_3barlink(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL);

  bgval    = param_getf_exit(ppl,"bgval");
  blankval = param_getf_dflt(ppl,"stim_blank_val",0.0);

  data = f3tensor(1,xn,1,yn,1,zn); // USE THIS FOR 3D FFT

  for(i=1;i<=zn;i++){
    if (((i-1)%tsamp) == 0){

      stimx_3barlink(ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,&tframe);

      //
      //  Copy this frame in and free it [WYETH - CHANGE THIS]
      //
      for(j=0;j<xn;j++){
	for(k=0;k<yn;k++){
	  data[j+1][k+1][i] = tframe[j][k];
	}
      }
      free_2d_farray(tframe,xn);

    }else{
      make_frame_flat(data,1,xn,1,yn,1,zn,i,blankval-bgval); // void value
      // This ends up as 'blankval' when 'bgval is added below
    }
  }
  add_const_3d_farray(data,1,xn,1,yn,1,zn,bgval);

  stimx_3barlink(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL);// Clean up

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_STIM_3D_PPL_PIXEL                           */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_pixel(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp)
     char *mylogf;
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i,j,k;
  float bgval,blankval,**tframl,***dl;

  //
  //  PREP
  //
  stimx_pixel(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL);

  bgval    = param_getf_exit(ppl,"bgval");
  blankval = param_getf_dflt(ppl,"stim_blank_val",0.0);

  dl = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  for(i=1;i<=zn;i++){
    if (((i-1)%tsamp) == 0){

      stimx_pixel(ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,&tframl);

      //
      //  Copy this frame in and free it [WYETH - CHANGE THIS]
      //
      for(j=0;j<xn;j++){
	for(k=0;k<yn;k++){
	  dl[j+1][k+1][i] = tframl[j][k];
	}
      }
      free_2d_farray(tframl,xn);

    }else{
      make_frame_flat(dl,1,xn,1,yn,1,zn,i,blankval-bgval); // void value
      // This ends up as 'blankval' when 'bgval is added below
    }
  }
  stimx_pixel(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL);  // Clean up

  return dl;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_STIM_3DB_PPL_BAR                           */
/*                                                                           */
/*****************************************************************************/
void get_stim_3db_ppl_bar(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,rdl,rdr)
     char *mylogf;
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
     float ****rdl,****rdr;
{
  int i,j,k;
  float bgval,blankval,**tframl,**tframr;
  float ***dl,***dr;

  //
  //  PREP
  //
  stimx_bar(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL,NULL);

  bgval    = param_getf_exit(ppl,"bgval");
  blankval = param_getf_dflt(ppl,"stim_blank_val",0.0);

  dl = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT
  dr = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  for(i=1;i<=zn;i++){
    if (((i-1)%tsamp) == 0){

      stimx_bar(ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,&tframl,&tframr);

      //
      //  Copy this frame in and free it [WYETH - CHANGE THIS]
      //
      for(j=0;j<xn;j++){
	for(k=0;k<yn;k++){
	  dl[j+1][k+1][i] = tframl[j][k];
	  dr[j+1][k+1][i] = tframr[j][k];
	}
      }
      free_2d_farray(tframl,xn);
      free_2d_farray(tframr,xn);

    }else{
      make_frame_flat(dl,1,xn,1,yn,1,zn,i,blankval-bgval); // void value
      make_frame_flat(dr,1,xn,1,yn,1,zn,i,blankval-bgval); // void value
      // This ends up as 'blankval' when 'bgval is added below
    }
  }
  add_const_3d_farray(dl,1,xn,1,yn,1,zn,bgval);
  add_const_3d_farray(dr,1,xn,1,yn,1,zn,bgval);

  stimx_bar(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL,NULL);  // Clean up

  *rdl = dl;
  *rdr = dr;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_STIM_3D_PPL_BAR2RAN                          */
/*                                                                           */
/*  Return a drifting bar.                                                   */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*    sf - spatial frequency (cyc/deg)                                       */
/*    tf - temporal frequency (cyc/sec)                                      */
/*    phase - phase offset (deg)                                             */
/*    theta - orientation of grating (deg)                                   */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_bar2ran(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp; // Rounded to integer, approximates "stim_samp"
{
  int i,j,k;
  float ***data,t,t0,tn,blankval,bgval,**tframe;

  //
  //  PREP
  //
  stimx_bar2ran(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL);

  t0       = param_getf_exit(ppl,"st0");
  tn       = param_getf_exit(ppl,"stn");
  bgval    = param_getf_exit(ppl,"bgval");
  blankval  = param_getf_dflt(ppl,"stim_blank_val",0.0);

  data = f3tensor(1,xn,1,yn,1,zn); // USE THIS FOR 3D FFT

  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;
    if (((i-1)%tsamp) == 0){
      if ((t >= t0) && (t < t0+tn)){

	stimx_bar2ran(ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,&tframe);

	//
	//  Copy this frame in and free it [WYETH - CHANGE THIS]
	//
	for(j=0;j<xn;j++){
	  for(k=0;k<yn;k++){
	    data[j+1][k+1][i] = tframe[j][k];
	  }
	}
	free_2d_farray(tframe,xn);

      }else
	make_frame_flat(data,1,xn,1,yn,1,zn,i,bgval);
    }else
      make_frame_flat(data,1,xn,1,yn,1,zn,i,blankval); // 0 is min value
  }

  stimx_bar2ran(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL);  // Clean up

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_STIM_3D_PPL_BARRAY_OLD                        */
/*                                                                           */
/*  *** WYETH - BEFORE 'stimx_...'  way                                      */
/*                                                                           */
/*  Return a drifting bar.                                                   */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*    sf - spatial frequency (cyc/deg)                                       */
/*    tf - temporal frequency (cyc/sec)                                      */
/*    phase - phase offset (deg)                                             */
/*    theta - orientation of grating (deg)                                   */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_barray_old(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp; // Rounded to integer, approximates "stim_samp"
{
  int i,j,k;
  int n,si,bi,cnt,seed,x0,y0,flag;
  int dt,dt_off,nbar,npos,seqtype,**seqi;
  float t0,tn,cx,cy,theta,barw,barl,contrast,bgval,maxlum,***data,t,x,y;
  float bar_gap,pos_gap,barval,size,**bar_cx,**bar_cy,sdaa;

  t0       = param_getf_exit(ppl,"st0");
  tn       = param_getf_exit(ppl,"stn");
  cx       = param_getf_exit(ppl,"cx");
  cy       = param_getf_exit(ppl,"cy");
  theta    = param_getf_exit(ppl,"direction");
  barw     = param_getf_exit(ppl,"bar_w");
  bar_gap  = param_getf_exit(ppl,"bar_gap");
  pos_gap  = param_getf_exit(ppl,"pos_gap");
  nbar     = paramfile_get_int_param_or_exit(ppl,"bar_n");
  npos     = paramfile_get_int_param_or_exit(ppl,"pos_n");
  dt       = paramfile_get_int_param_or_exit(ppl,"dt");
  dt_off   = paramfile_get_int_param_or_exit(ppl,"dt_off");
  seed     = paramfile_get_int_param_or_exit(ppl,"seed");
  contrast = param_getf_exit(ppl,"contrast");
  barval   = param_getf_exit(ppl,"bar_amp");
  bgval    = param_getf_exit(ppl,"bgval");
  maxlum   = param_getf_exit(ppl,"maxlum");
  seqtype  = paramfile_get_int_param_or_exit(ppl,"seq_type");
  size     = param_getf_exit(ppl,"size");
  sdaa     = paramfile_get_float_param_default(ppl,"sdaa",0.0);

  barl = (size - (float)(nbar-1)*bar_gap)/(float)nbar; // Bar length (degr)
  n = (int)(tn/tscale / (float)(dt+dt_off)); // Bar patterns in random sequ.

  // Get the random sequence [n][nbar]
  //seqi = get_zero_2d_iarray(n,nbar);  DEBUG
  seqi = stm_barray_get_seq(NULL,n,nbar,npos,seqtype,seed);

  /**
  for(i=0;i<10;i++)
    for(j=0;j<nbar;j++)
      printf(" %d bar %d   %d\n",i,j,seqi[i][j]);
      exit(0);***/

  // Get center coords for each possible bar in the grid [nbar][npos]
  // in degrees, w.r.t. 0,0  (Thus x,y, offsets must be added)
  stm_barray_get_bar_centers(NULL,nbar,npos,theta,barw,barl,bar_gap,
			     pos_gap,&bar_cx,&bar_cy);

  /***
  for(i=0;i<nbar;i++)
    for(j=0;j<npos;j++)
      printf(" bar pos %d %d   xy  %f %f\n",i,j,bar_cx[i][j],bar_cy[i][j]);
  ***/
  

  data = f3tensor(1,xn,1,yn,1,zn); // for 3D FFT

  si = 0;  // Sequence index
  cnt = 0;
  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;
    //if (0 == 1){
    if (((i-1)%tsamp) == 0){
      if ((t >= t0) && (t < t0+tn)){

	if (cnt >= dt){
	  make_frame_flat(data,1,xn,1,yn,1,zn,i,bgval);
	}else{

	  make_frame_flat(data,1,xn,1,yn,1,zn,i,bgval);

	  for(j=0;j<nbar;j++){ // Draw each bar (or gray if no bar)
	    if (si >= n){
	      printf("  si = %d\n",si);
	      exit_error("GET_STIM_3D_PPL_BARRAY","Sequence too short");
	    }
	    
	    bi = seqi[si][j] - 1;
	    //printf("  %d bar %d    %d\n",i,j,bi);
	    if (bi >= 0){ // This bar is drawn

	      if (seqtype == 99){ // Draw all bars - for demo
		for(k=0;k<npos;k++){
		  x = bar_cx[j][k] + cx;
		  y = bar_cy[j][k] + cy;
		  make_frame_bar(data,1,xn,1,yn,1,zn,i,x,y,theta,barl,barw,
				 sdaa,sscale,barval,bgval,1);
		}
	      }else{
		x = bar_cx[j][bi] + cx;
		y = bar_cy[j][bi] + cy;
		
		//printf("i=%d bar=%d  pos=%d  x,y  %f %f\n",i,j,bi,x,y);
		
		make_frame_bar(data,1,xn,1,yn,1,zn,i,x,y,theta,barl,barw,
			       sdaa,sscale,barval,bgval,1);
	      }
	    }
	  }
	}
	cnt += 1;
	if (cnt >= (dt+dt_off)){
	  cnt = 0;
	  si += 1;  // Sequence index
	}
      }else
	make_frame_flat(data,1,xn,1,yn,1,zn,i,bgval);
    }else
      make_frame_flat(data,1,xn,1,yn,1,zn,i,0.0); // 0 is min value
  }

  free_2d_farray(bar_cx,nbar);
  free_2d_farray(bar_cy,nbar);
  free_2d_iarray(seqi,n);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_STIM_3D_PPL_BARRAY                          */
/*                                                                           */
/*  Return a drifting bar.                                                   */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*    sf - spatial frequency (cyc/deg)                                       */
/*    tf - temporal frequency (cyc/sec)                                      */
/*    phase - phase offset (deg)                                             */
/*    theta - orientation of grating (deg)                                   */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_barray(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp; // Rounded to integer, approximates "stim_samp"
{
  int i,j,k;
  float ***data,t,t0,tn,bgval,**tframe;

  //
  //  PREP
  //
  stimx_barray(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL);

  t0    = param_getf_exit(ppl,"st0");
  tn    = param_getf_exit(ppl,"stn");
  bgval = param_getf_exit(ppl,"bgval");

  data = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;  // Time (s)
    if (((i-1)%tsamp) == 0){  // Is this a video frame
      if ((t >= t0) && (t < t0+tn)){

	stimx_barray(ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,&tframe);

	//
	//  Copy this frame in and free it [WYETH - CHANGE THIS]
	//
	for(j=0;j<xn;j++){
	  for(k=0;k<yn;k++){
	    data[j+1][k+1][i] = tframe[j][k];
	  }
	}
	free_2d_farray(tframe,xn);

      }else{
	make_frame_flat(data,1,xn,1,yn,1,zn,i,bgval);
      }
    }else{
      //make_frame_flat(data,1,xn,1,yn,1,zn,i,-1.0); // typ. black, -1.0
      make_frame_flat(data,1,xn,1,yn,1,zn,i,0.0);
    }
  }

  stimx_barray(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL);  // Clean up

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_STIM_3D_PPL_SURRT                         */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*    sf - spatial frequency (cyc/deg)                                       */
/*    tf - temporal frequency (cyc/sec)                                      */
/*    phase - phase offset (deg)                                             */
/*    theta - orientation of grating (deg)                                   */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned storage:   */
/*        free_f3tensor(data,1,xn,1,yn,1,zn).                                */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_surrt(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i;
  int i0,cs_flag;
  float t0,tn,t20,t2n,cx,cy,t,cn1,cn2;
  float ***data,ph1,ph2,mid;
  float size,idiam,odiam,sf1,sf2,tf1,tf2,phase1,phase2,theta1,theta2;
  float con1,con2,maxlum,bgamp,bgval;

  t0       = param_getf_exit(ppl,"st0");
  tn       = param_getf_exit(ppl,"stn");
  cx       = param_getf_exit(ppl,"cx");
  cy       = param_getf_exit(ppl,"cy");
  t20      = param_getf_exit(ppl,"s2t0");
  t2n      = param_getf_exit(ppl,"s2tn");
  size     = param_getf_exit(ppl,"size");
  idiam    = param_getf_exit(ppl,"idiam");
  odiam    = param_getf_exit(ppl,"odiam");
  sf1      = param_getf_exit(ppl,"sf1");
  sf2      = param_getf_exit(ppl,"sf2");
  tf1      = param_getf_exit(ppl,"tf1");
  tf2      = param_getf_exit(ppl,"tf2");
  phase1   = param_getf_exit(ppl,"phase1");
  phase2   = param_getf_exit(ppl,"phase2");
  theta1   = param_getf_exit(ppl,"direction1");
  theta2   = param_getf_exit(ppl,"direction2");
  con1     = param_getf_exit(ppl,"contrast1");
  con2     = param_getf_exit(ppl,"contrast2");
  maxlum   = param_getf_exit(ppl,"maxlum");
  bgval    = param_getf_exit(ppl,"bgval");

  //   S C  flag
  //   0 0   0   neither
  //   0 1   1   center
  //   1 0   2   surround
  //   1 1   3   both
  cs_flag  = paramfile_get_int_param_default(ppl,"cs_flag",3);


  data = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  bgamp = bgval * 2.0 - 1.0; // Convert 0..1 --> -1..1
  i0 = 1 + my_rint(t0 / tscale);
  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;
    if (((i-1)%tsamp) == 0){

      if ((t >= t0) && (t < t0+tn)) // Center is on or off
	cn1 = con1;
      else
	cn1 = 0.0;
      
      if ((t >= t20) && (t < t20+t2n)) // Surround is on or off
	cn2 = con2;
      else
	cn2 = 0.0;

      if ((cn1 != 0.0) || (cn2 != 0.0)){
	ph1 = phase1 + (float)(i-i0)*tscale*tf1 * 360.0;
	ph2 = phase2 + (float)(i-i0)*tscale*tf2 * 360.0;
	make_frame_circ_cs_sine(data,1,xn,1,yn,1,zn,i,size,idiam,odiam,sf1,sf2,
				theta1,theta2,ph1,ph2,sscale,cn1,cn2,bgamp,
				cx,cy,cs_flag);
      }else
	make_frame_flat(data,1,xn,1,yn,1,zn,i,bgamp);
    }else{
      make_frame_flat(data,1,xn,1,yn,1,zn,i,-1.0); // black
    }
  }
  mid = (float)maxlum/2.0;
  multiply_3d_farray(data,1,xn,1,yn,1,zn,mid);
  add_const_3d_farray(data,1,xn,1,yn,1,zn,mid);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_STIM_3D_PPL_WY4ST                         */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*    sf - spatial frequency (cyc/deg)                                       */
/*    tf - temporal frequency (cyc/sec)                                      */
/*    phase - phase offset (deg)                                             */
/*    theta - orientation of grating (deg)                                   */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned storage:   */
/*        free_f3tensor(data,1,xn,1,yn,1,zn).                                */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_wy4st(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i,k;
  int dwell,n,*seqi,seed,state,cnt,nst;
  float ***data,ph1,ph2,mid,ttc,tts,t0,tn,t;
  float size,idiam,odiam,sf1,sf2,tf1,tf2,phase1,phase2,theta1,theta2;
  float con1,con2,maxlum,bgamp,bgval,cx,cy;
  char *ststr;

  cx       = param_getf_exit(ppl,"cx");
  cy       = param_getf_exit(ppl,"cy");

  size     = param_getf_exit(ppl,"size");
  idiam    = param_getf_exit(ppl,"idiam");
  odiam    = param_getf_exit(ppl,"odiam");
  sf1      = param_getf_exit(ppl,"sf1");
  sf2      = param_getf_exit(ppl,"sf2");
  tf1      = param_getf_exit(ppl,"tf1");
  tf2      = param_getf_exit(ppl,"tf2");
  phase1   = param_getf_exit(ppl,"phase1");
  phase2   = param_getf_exit(ppl,"phase2");
  theta1   = param_getf_exit(ppl,"direction1");
  theta2   = param_getf_exit(ppl,"direction2");
  con1     = param_getf_exit(ppl,"contrast1");
  con2     = param_getf_exit(ppl,"contrast2");
  maxlum   = param_getf_exit(ppl,"maxlum");
  bgval    = param_getf_exit(ppl,"bgval");
  seed     = paramfile_get_int_param_or_exit(ppl,"seed");
  ststr    = paramfile_get_char_param_default(ppl,"state_seq",NULL);

  t0       = param_getf_exit(ppl,"st0");
  nst      = paramfile_get_float_param_default(ppl,"stn_state",-1);
  if (nst == -1)
    tn     = param_getf_exit(ppl,"stn");
  else{
    tn = tscale * zn;
    //printf("    SETTING tn to %f sec\n",tn);
  }

  bgamp = bgval * 2.0 - 1.0; // Convert 0..1 --> -1..1

  data = f3tensor(1,xn,1,yn,1,zn); // USE THIS FOR 3D FFT

  //
  //  WYETH - should be careful about dwell here.
  //
  dwell = my_rint(1.0/(tf1*tscale));  // Number of frames per cycle
  //printf("    TF1      = %f\n",tf1);
  //printf("    Dwell    = %d\n",dwell);
  //printf("    TF ideal = %f\n",1.0/(tscale * (float)dwell));

  // Determine duration 'n' in number of states, get sequence
  if (nst >= 0){
    n = nst;
  }else{
    n = zn/dwell;
  }
  stm_wy4st_get_seq(NULL,n,seed,ststr,&seqi); // Get state sequence

  if (seed > 0)
    seed = -seed;

  k = 0;        // Current state index
  cnt = dwell;  // Dwell counter
  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;
    if (((i-1)%tsamp) == 0){
      if ((t >= t0) && (t < t0+tn) &&              // In time window
	  ((k < n) || ((k==n) && (cnt < dwell)))){ // Not done w/ final state

	if (cnt == dwell){
	  state = seqi[k];
	  if ((state == 0)||(state == 1))
	    ttc = theta1 + 90.0;
	  else
	    ttc = theta1;
	  if ((state == 0)||(state == 2))
	    tts = theta2 + 90.0;
	  else
	    tts = theta2;

	  cnt = 1; // This will be first frame shown of this state
	  k += 1;  // Use next state next time
	}else{
	  cnt += 1; // How many times the state has been used
	}
	
	ph1 = phase1 + (float)(i-1)*tscale*tf1 * 360.0;
	ph2 = phase2 + (float)(i-1)*tscale*tf2 * 360.0;
	make_frame_circ_cs_sine(data,1,xn,1,yn,1,zn,i,size,idiam,odiam,sf1,sf2,
				ttc,tts,ph1,ph2,sscale,con1,con2,bgamp,cx,cy,
				3);
      }else{
	make_frame_flat(data,1,xn,1,yn,1,zn,i,bgamp);  // Plain background
      }
    }else{
      make_frame_flat(data,1,xn,1,yn,1,zn,i,-1.0); // Black
    }
  }
  mid = (float)maxlum/2.0;
  multiply_3d_farray(data,1,xn,1,yn,1,zn,mid);
  add_const_3d_farray(data,1,xn,1,yn,1,zn,mid);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_STIM_3D_PPL_X4ST                          */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*    sf - spatial frequency (cyc/deg)                                       */
/*    tf - temporal frequency (cyc/sec)                                      */
/*    phase - phase offset (deg)                                             */
/*    theta - orientation of grating (deg)                                   */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned storage:   */
/*        free_f3tensor(data,1,xn,1,yn,1,zn).                                */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_x4st(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i,k;
  int dwell,n,*seqi,seed,state,cnt,nst,aptype;
  float ***data,***data1,***data2,ph1,ph2,mid,t0,tn,t,sz;
  float size,sf1,sf2,tf1,tf2,phase1,phase2,theta1,theta2;
  float con1,con2,c1,c2,maxlum,bgamp,bgval,cx,cy,xc,yc;
  char *ststr;

  cx       = paramfile_get_float_param_default(ppl,"cx",0.0);
  cy       = paramfile_get_float_param_default(ppl,"cy",0.0);

  size     = param_getf_exit(ppl,"size");
  sf1      = param_getf_exit(ppl,"sf1");
  sf2      = param_getf_exit(ppl,"sf2");
  tf1      = param_getf_exit(ppl,"tf1");
  tf2      = param_getf_exit(ppl,"tf2");
  phase1   = param_getf_exit(ppl,"phase1");
  phase2   = param_getf_exit(ppl,"phase2");
  theta1   = param_getf_exit(ppl,"direction1");
  theta2   = param_getf_exit(ppl,"direction2");
  con1     = param_getf_exit(ppl,"contrast1");
  con2     = param_getf_exit(ppl,"contrast2");
  maxlum   = param_getf_exit(ppl,"maxlum");
  bgval    = param_getf_exit(ppl,"bgval");
  seed     = paramfile_get_int_param_or_exit(ppl,"seed");
  ststr    = paramfile_get_char_param_default(ppl,"state_seq",NULL);
  aptype   = paramfile_get_int_param_or_exit(ppl,"aptype");

  t0       = param_getf_exit(ppl,"st0");
  nst      = paramfile_get_float_param_default(ppl,"stn_state",-1);
  if (nst == -1)
    tn     = param_getf_exit(ppl,"stn");
  else{
    tn = tscale * zn;
    //printf("    SETTING tn to %f sec\n",tn);
  }

  bgamp = bgval * 2.0 - 1.0; // Convert 0..1 --> -1..1
  sz = size/sscale;

  data1 = f3tensor(1,xn,1,yn,1,zn); // USE THIS FOR 3D FFT
  data2 = f3tensor(1,xn,1,yn,1,zn);

  xc = (float)(xn-1)/2.0 + cx/sscale;  // Center w.r.t. 0..xn-1 (pix)
  yc = (float)(yn-1)/2.0 + cy/sscale;  // Center w.r.t. 0..yn-1 (pix)

  //
  //  WYETH - should be careful about dwell here.
  //
  dwell = my_rint(1.0/(tf1*tscale));  // Number of frames per cycle
  //printf("    TF1      = %f\n",tf1);
  //printf("    Dwell    = %d\n",dwell);
  //printf("    TF ideal = %f\n",1.0/(tscale * (float)dwell));

  // Determine duration 'n' in number of states, get sequence
  if (nst >= 0){
    n = nst;
  }else{
    n = zn/dwell;
  }
  stm_wy4st_get_seq(NULL,n,seed,ststr,&seqi); // Get state sequence

  if (seed > 0)
    seed = -seed;

  k = 0;        // Current state index
  cnt = dwell;  // Dwell counter
  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;
    if (((i-1)%tsamp) == 0){
      if ((t >= t0) && (t < t0+tn) &&              // In time window
	  ((k < n) || ((k==n) && (cnt < dwell)))){ // Not done w/ final state

	if (cnt == dwell){
	  state = seqi[k];

	  c1 = con1;
	  if ((state == 0)||(state==1))
	    c1 = 0.0;
	  c2 = con2;
	  if ((state == 0)||(state==2))
	    c2 = 0.0;

	  cnt = 1; // This will be first frame shown of this state
	  k += 1;  // Use next state next time
	}else{
	  cnt += 1; // How many times the state has been used
	}
	
	ph1 = phase1 + (float)(i-1)*tscale*tf1 * 360.0;
	ph2 = phase2 + (float)(i-1)*tscale*tf2 * 360.0;

	make_frame_sinusoid(data1,1,xn,1,yn,1,zn,i,sf1,theta1,sscale,ph1,
			    cx,cy);
	make_frame_sinusoid(data2,1,xn,1,yn,1,zn,i,sf2,theta2,sscale,ph2,
			    cx,cy);
	multiply_frame(data1,1,xn,1,yn,1,zn,i,c1);
	multiply_frame(data2,1,xn,1,yn,1,zn,i,c2);
	frame_apply_aperture(data1,1,xn,1,yn,1,zn,i,aptype,xc,yc,
			     sz,sz,0.0,bgamp);
	frame_apply_aperture(data2,1,xn,1,yn,1,zn,i,aptype,xc,yc,
			     sz,sz,0.0,bgamp);
      }else{
	make_frame_flat(data1,1,xn,1,yn,1,zn,i,bgamp);  // Plain background
	make_frame_flat(data2,1,xn,1,yn,1,zn,i,bgamp);  // Plain background
      }
    }else{
      make_frame_flat(data1,1,xn,1,yn,1,zn,i,-1.0); // Black
      make_frame_flat(data2,1,xn,1,yn,1,zn,i,-1.0); // Black
    }
  }

  mid = (float)maxlum/2.0;
  multiply_3d_farray(data1,1,xn,1,yn,1,zn,mid);
  multiply_3d_farray(data2,1,xn,1,yn,1,zn,mid);
  // At c=1, values now go from -mid to mid

  data = add_3d_tensors(data1,data2,xn,yn,zn);
  free_f3tensor(data1,1,xn,1,yn,1,zn);
  free_f3tensor(data2,1,xn,1,yn,1,zn);

  // At c=1, values go from -maxlum to maxlum
  // Thus, c=0.5 should be largest contrast used to avoid truncation next

  stim_util_limit_3d(data,1,xn,1,yn,1,zn,-mid,mid);
  add_const_3d_farray(data,1,xn,1,yn,1,zn,mid); // Values should be 0..maxlum

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*    WYETH REMOVE               GET_STIM_3D_PPL_X4ST                        */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*    sf - spatial frequency (cyc/deg)                                       */
/*    tf - temporal frequency (cyc/sec)                                      */
/*    phase - phase offset (deg)                                             */
/*    theta - orientation of grating (deg)                                   */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned storage:   */
/*        free_f3tensor(data,1,xn,1,yn,1,zn).                                */
/*                                                                           */
/*  - Use contrast 1 for maximum.                                            */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_x4st_OLD(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i;
  int state,dwell,aptype,seed,i0;
  float ***data,***data1,***data2,mid,rv,con1,con2,c1,c2,maxlum,bgampl,bgval;
  float size,sz,sf1,sf2,tf1,tf2,phase1,phase2,ph1,ph2,theta1,theta2,tt,t0,tn,t;
  float cx,cy,xc,yc;

  printf("  GET_STIM_3D_PPL_X4ST\n");

  cx       = paramfile_get_float_param_default(ppl,"cx",0.0);
  cy       = paramfile_get_float_param_default(ppl,"cy",0.0);
  size     = param_getf_exit(ppl,"size");
  sf1      = param_getf_exit(ppl,"sf1");
  sf2      = param_getf_exit(ppl,"sf2");
  tf1      = param_getf_exit(ppl,"tf1");
  tf2      = param_getf_exit(ppl,"tf2");
  phase1   = param_getf_exit(ppl,"phase1");
  phase2   = param_getf_exit(ppl,"phase2");
  theta1   = param_getf_exit(ppl,"direction1");
  theta2   = param_getf_exit(ppl,"direction2");
  con1     = param_getf_exit(ppl,"contrast1");
  con2     = param_getf_exit(ppl,"contrast2");
  maxlum   = param_getf_exit(ppl,"maxlum");
  bgval    = param_getf_exit(ppl,"bgval");
  aptype   =   paramfile_get_int_param_or_exit(ppl,"aptype");
  seed     =   paramfile_get_int_param_or_exit(ppl,"seed");
  t0       = param_getf_exit(ppl,"st0");
  tn       = param_getf_exit(ppl,"stn");

  bgampl = bgval * 2.0 - 1.0; /* Convert 0..1 --> -1..1 */

  if (seed > 0)
    seed = -seed;

  sz = size/sscale;
  dwell = my_rint(1.0/(tf1*tscale));
  tt = dwell*tf1*tscale;
  printf("    Dwell = %d (1/tf*tscale = %f)\n",dwell,1.0/(tf1*tscale));

  if (((tt-1.0) > 0.001)||((tt-1.0) < -0.001))
    exit_error("GET_STIM_3D_PPL_X4ST","Dwell, TF and tscale are incompatible");

  data1 = f3tensor(1,xn,1,yn,1,zn);
  data2 = f3tensor(1,xn,1,yn,1,zn);

  xc = (float)(xn-1)/2.0 + cx/sscale;  // Center w.r.t. 0..xn-1 (pix)
  yc = (float)(yn-1)/2.0 + cy/sscale;  // Center w.r.t. 0..yn-1 (pix)
  
  state = 0;
  i0 = 1 + my_rint(t0 / tscale);
  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;
    if (((i-i0)%dwell) == 0){
      rv = myrand_util_ran2(&seed);
      if (rv < 0.25)
	state = 3;
      else if (rv < 0.50)
	state = 2;
      else if (rv < 0.75)
	state = 1;
      else
	state = 0;
    }

    printf("WYETH - Does this routine work?, not last time I checked\n");
    printf("State = %d\n",state);
    c1 = con1;
    if ((state == 0)||(state==1))
      c1 = 0.0;
    c2 = con2;
    if ((state == 0)||(state==2))
      c2 = 0.0;

    printf("con1,2  = %f %f\n",con1,con2);

    if (((i-1)%tsamp) == 0){
      if ((t >= t0) && (t < t0+tn)){
	ph1 = phase1 + (float)(i-i0)*tscale*tf1 * 360.0;
	ph2 = phase2 + (float)(i-i0)*tscale*tf2 * 360.0;
	
	make_frame_sinusoid(data1,1,xn,1,yn,1,zn,i,sf1,theta1,sscale,ph1,
			    cx,cy);
	make_frame_sinusoid(data2,1,xn,1,yn,1,zn,i,sf2,theta2,sscale,ph2,
			    cx,cy);
	
	multiply_frame(data1,1,xn,1,yn,1,zn,i,c1);
	multiply_frame(data2,1,xn,1,yn,1,zn,i,c2);
	
	frame_apply_aperture(data1,1,xn,1,yn,1,zn,i,aptype,xc,yc,
			     sz,sz,0.0,bgampl);
	frame_apply_aperture(data2,1,xn,1,yn,1,zn,i,aptype,xc,yc,
			     sz,sz,0.0,bgampl);
      }else{
	make_frame_flat(data1,1,xn,1,yn,1,zn,i,bgampl);
	make_frame_flat(data2,1,xn,1,yn,1,zn,i,bgampl);
      }
    }else{
      make_frame_flat(data1,1,xn,1,yn,1,zn,i,-1.0);
      make_frame_flat(data2,1,xn,1,yn,1,zn,i,-1.0);
    }
  }
  
  mid = (float)maxlum/2.0;
  multiply_3d_farray(data1,1,xn,1,yn,1,zn,mid);
  multiply_3d_farray(data2,1,xn,1,yn,1,zn,mid);
  // At c=1, values now go from -mid to mid

  data = add_3d_tensors(data1,data2,xn,yn,zn);
  free_f3tensor(data1,1,xn,1,yn,1,zn);
  free_f3tensor(data2,1,xn,1,yn,1,zn);

  /* At c=1, values go from -maxlum to maxlum */
  /* Thus, c=0.5 should be largest contrast used to avoid truncation next */

  stim_util_limit_3d(data,1,xn,1,yn,1,zn,-mid,mid);
  add_const_3d_farray(data,1,xn,1,yn,1,zn,mid);
  /* Final values should be in 0..maxlum */
  
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_STIM_3D_PPL_WY2PA                           */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*    sf - spatial frequency (cyc/deg)                                       */
/*    tf - temporal frequency (cyc/sec)                                      */
/*    phase - phase offset (deg)                                             */
/*    theta - orientation of grating (deg)                                   */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned storage:   */
/*        free_f3tensor(data,1,xn,1,yn,1,zn).                                */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_wy2pa(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i;
  int count,full,tstatic,pflag,maxj1,minj1,maxj2,minj2;
  int toverlap,prefdur,antidur,antimax,p0,p1,a0,a1,tanti,apat,ppat,jump1,jump2;
  float ***data,ph1,ph2,mid,cx,cy,tcx,tcy,r,theta;
  float size,pdist,ptheta,sf1,sf2,tf1,tf2,phase1,phase2,theta1,theta2;
  float maxlum,bgampl,con1,con2,bgval;
  int seed;

  size    = param_getf_exit(ppl,"size");
  pdist   = param_getf_exit(ppl,"pdist");
  ptheta  = param_getf_exit(ppl,"ptheta");
  sf1     = param_getf_exit(ppl,"sf1");
  sf2     = param_getf_exit(ppl,"sf2");
  tf1     = param_getf_exit(ppl,"tf1");
  tf2     = param_getf_exit(ppl,"tf2");
  phase1     = param_getf_exit(ppl,"phase1");
  phase2     = param_getf_exit(ppl,"phase2");
  theta1     = param_getf_exit(ppl,"direction1");
  theta2     = param_getf_exit(ppl,"direction2");
  con1     = param_getf_exit(ppl,"contrast1");
  con2     = param_getf_exit(ppl,"contrast2");
  maxlum   = param_getf_exit(ppl,"maxlum");
  bgval   = param_getf_exit(ppl,"bgval");
  seed     = paramfile_get_int_param_or_exit(ppl,"seed");

  bgampl = bgval * 2.0 - 1.0; /* Convert 0..1 --> -1..1 */


  printf("  GET_STIM_3D_PPL_WY2PA\n");
  printf("    Only partially implemented --- many params ignored.\n");

  printf("  SCALE = %f\n",sscale);
  printf("  size = %f\n",size);
  printf("  phase1,phase2 = %f %f\n",phase1,phase2);

  data = f3tensor(1,xn,1,yn,1,zn); /*** USE THIS FOR 3D FFT ***/

  if (seed > 0)
    seed = -seed;

  prefdur = 4;   /* Should be sent in */
  antidur = 0;   /* Should be sent in */
  toverlap = 0;  /* Should be sent in */
  cx = cy = 0.0; /* These could be sent in */
  
  maxj1 =  tf1*tscale * 360.0;  /* Degrees of phase per time unit */
  minj1 = -tf1*tscale * 360.0;  /* Degrees of phase per time unit */
  maxj2 =  tf2*tscale * 360.0;  /* Degrees of phase per time unit */
  minj2 = -tf2*tscale * 360.0;  /* Degrees of phase per time unit */
  
  tstatic = 4;
  pflag = 0;
  if (antidur < 0){ /* Use this for simultaneous pref stimuli */
    pflag = 1;
    tstatic = -antidur;
    antidur = prefdur;
    antimax = 0;
  }else if (antidur == 0) /* Use this to implement other anti sequences */
    antimax = 10;
  else
    antimax = antidur;
  
  if (pflag){
    p0 = tstatic;
    p1 = p0+prefdur-1;
    a1 = p1;
    full = tstatic + prefdur;
    toverlap = 0;
  }else{
    p0 = antimax + tstatic - toverlap;    /* First frame of pref */
    p1 = p0+prefdur-1;                    /* Last frame of pref */
    full = p0+prefdur + 2;
    a1 = p0-1 + toverlap;                 /* Last frame of antipref */
  }

  a0 = apat = ppat = 0; /* Avoid -Wall warning */
  
  ph1 = ph2 = 0.0;
  count = full;
  for(i=0;i<zn;i++){
    if (count == full){
      count = 0;
      if (antidur == 0){
	/* tanti is duration of ANTI pulse */
	tanti = (int)(myrand_util_ran2(&seed)*8.0);
	if (tanti == 5)
	  tanti = 6;
	else if (tanti == 6)
	  tanti = 8;
	else if (tanti == 7)
	  tanti = 10;
      }else{
	if (myrand_util_ran2(&seed) < 0.5)
	  tanti = 0;
	else
	  tanti = antidur;
      }
      a0 = a1 - tanti + 1;
      if (myrand_util_ran2(&seed) < 0.5)   /* 'apat' patch does anti motion */
	apat = 0;
      else
	apat = 1;

      if (myrand_util_ran2(&seed) < 0.5)   /* 'ppat' patch does pref motion */
	ppat = 0;
      else
	ppat = 1;
    }
    
    jump1 = jump2 = 0;
    if ((count >= a0)&&(count <= a1)){  /* Show anti */
      if (apat)
	jump1 = minj1;
      else
	jump2 = minj2;
    }
    if ((count >= p0)&&(count <= p1)){  /* Show pref (may overwrite anti) */
      if (ppat)
	jump1 = maxj1;
      else
	jump2 = maxj2;
    }
    ph1 -= jump1;
    ph2 -= jump2;
    
    if (((i-1)%tsamp) == 0){
      make_frame_flat(data,1,xn,1,yn,1,zn,i,bgampl);
      
      r = (float)pdist/2.0;
      theta = (float)ptheta*M_PI/180.0;
      
      tcx = cx + r*cos(theta);
      tcy = cy + r*sin(theta); // Negate y-coords
      draw_circ_sine_in_frame(data,1,xn,1,yn,1,zn,i,tcx,tcy,size,sf1,theta1,
			      ph1,sscale,con1);
      tcx = cx - r*cos(theta);
      tcy = cy - r*sin(theta); // Negate y-coords
      draw_circ_sine_in_frame(data,1,xn,1,yn,1,zn,i,tcx,tcy,size,sf2,theta2,
			      ph2,sscale,con2);
    }else
      make_frame_flat(data,1,xn,1,yn,1,zn,i,-1.0);

    count += 1;
  }
  mid = (float)maxlum/2.0;
  multiply_3d_farray(data,1,xn,1,yn,1,zn,mid);
  add_const_3d_farray(data,1,xn,1,yn,1,zn,mid);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_STIM_3D_PPL_GP                            */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*    r - distance between dots (degr)                                       */
/*    theta - orientation of dot pair                                        */
/*    amp1 - luminance of dot 1                                              */
/*    amp2 - luminance of dot 2                                              */
/*    ampbg - background luminance                                           */
/*    ampdbg3 - background luminance                                         */
/*    dsize - dotsize (SD of Gaussian dot)                                   */
/*    ddens - dotpairs per square degree per frame                           */
/*    aptype - 0-none, 1-circ, 2-rect, 3-tiltrec                             */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Returned array dimensions start at 1.  To free returned storage:       */
/*    free_f3tensor(data,1,xn,1,yn,1,zn).                                    */
/*  - Do not reseed here, the caller must initialize the seed                */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_gp(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i,j,k;
  float ***data,t,t0,tn,bgval,**tframe;

  //
  //  PREP
  //
  stimx_gp(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL);

  t0        = param_getf_exit(ppl,"st0");
  tn        = param_getf_exit(ppl,"stn");
  bgval     = param_getf_exit(ppl,"bgval");

  data = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  for(i=1;i<=zn;i++){
    if (((i-1)%tsamp) == 0){

      stimx_gp(ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,&tframe);

      //
      //  Copy this frame in and free it [WYETH - CHANGE THIS]
      //
      for(j=0;j<xn;j++){
	for(k=0;k<yn;k++){
	  data[j+1][k+1][i] = tframe[j][k];
	}
      }
      free_2d_farray(tframe,xn);

    }else
      make_frame_flat(data,1,xn,1,yn,1,zn,i,bgval);
  }

  stimx_gp(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL);  // Clean up
  
  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_STIM_3D_PPL_NOISE                          */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Returned array dimensions start at 1.  To free returned storage:       */
/*    free_f3tensor(data,1,xn,1,yn,1,zn).                                    */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_noise(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;  // (deg/pix), (s/frame)
     int tsamp;            // Time units per stim frame
{
  int i,j,k;
  float ***data,mid,contrast,maxlum,bgval,bgamp,**tframe;

  //  For a class of stimuli that is like this one,
  //  For a class of stimuli that is like this one,
  // ******* WYETH This is the model stimulus, OTHERS SHOULD BE LIKE THIS
  // ******* WYETH This is the model stimulus, OTHERS SHOULD BE LIKE THIS
  // ******* WYETH This is the model stimulus, OTHERS SHOULD BE LIKE THIS
  // ******* WYETH This is the model stimulus, OTHERS SHOULD BE LIKE THIS
  //   Noted on 2011, May 14th


  //
  //  PREP
  //
  stimx_noise(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL);

  contrast  = param_getf_exit(ppl,"contrast");
  maxlum    = param_getf_exit(ppl,"maxlum");
  bgval     = param_getf_exit(ppl,"bgval");

  bgamp = 2.0*bgval - 1.0; // Convert 0..1 --> -1..1

  data = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  for(i=0;i<zn;i++){
    if ((i % tsamp) == 0){  // Is this a pattern frame
      stimx_noise(ppl,xn,yn,zn,sscale,tscale,tsamp,i,1,&tframe);
      if (tframe != NULL){
	for(j=0;j<xn;j++)  // Copy this frame and free [WYETH - CHANGE THIS]
	  for(k=0;k<yn;k++)
	    data[j+1][k+1][i+1] = tframe[j][k];
	free_2d_farray(tframe,xn);
      }else{
	make_frame_flat(data,1,xn,1,yn,1,zn,i+1,bgamp);
      }
    }else{
      make_frame_flat(data,1,xn,1,yn,1,zn,i+1,-1.0); // typ. black, -1.0
    }
  }

  if (param_test(ppl,"power_law_3d")){
    float ***pow_data;
    float fpow,min,max;

    // This should be done to the stimulus random sequence, all computed
    // in the PREP part of stimx_noise...
    exit_error("GET_STIM_3D_PPL_NOISE","WYETH - CHANGE IMPLEMENTATION");

    fpow   = param_getf_exit(ppl,"power_law_3d");

    pow_data = apply_power_law_3d(data,xn,yn,zn,fpow,1);  // 1-tensor_flag

    free_f3tensor(data,1,xn,1,yn,1,zn);
    data = pow_data;

    get_min_max_3d_farray(data,1,xn,1,yn,1,zn,&min,&max);
    if (-min > max)
      multiply_3d_farray(data,1,xn,1,yn,1,zn,1.0/-min);
    else
      multiply_3d_farray(data,1,xn,1,yn,1,zn,1.0/max);
  }


  // Assume all values in 'data' are in [-1,1]
  mid = (float)maxlum/2.0;
  multiply_3d_farray(data,1,xn,1,yn,1,zn,contrast*mid);
  add_const_3d_farray(data,1,xn,1,yn,1,zn,mid);

  stimx_noise(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL);  // Clean up

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_STIM_3D_PPL_IMAGE_SET                        */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Returned array dimensions start at 1.  To free returned storage:       */
/*    free_f3tensor(data,1,xn,1,yn,1,zn).                                    */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_image_set(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;  // (deg/pix), (s/frame)
     int tsamp;            // Time units per stim frame
{
  int i,j,k;
  float ***data,mid,contrast,maxlum,bgval,bgamp,**tframe;

  //
  //  PREP
  //
  stimx_image_set(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL);

  contrast  = param_getf_exit(ppl,"contrast");
  maxlum    = param_getf_exit(ppl,"maxlum");
  bgval     = param_getf_exit(ppl,"bgval");

  bgamp = 2.0*bgval - 1.0; // Convert 0..1 --> -1..1

  data = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  for(i=0;i<zn;i++){
    if ((i % tsamp) == 0){  // Is this a pattern frame
      stimx_image_set(ppl,xn,yn,zn,sscale,tscale,tsamp,i,1,&tframe);
      //
      // *** 'stimx_image_set' RE-USES 'tframe', *** MUST COPY *** from it
      // *** and DO NOT FREE 'tframe'
      //
      if (tframe != NULL){
	for(j=0;j<xn;j++)
	  for(k=0;k<yn;k++)
	    data[j+1][k+1][i+1] = tframe[j][k];
      }else{
	make_frame_flat(data,1,xn,1,yn,1,zn,i+1,bgamp);
      }
    }else{
      make_frame_flat(data,1,xn,1,yn,1,zn,i+1,-1.0); // typ. black, -1.0
    }
  }

  // Assume all values in 'data' are in [-1,1]
  mid = (float)maxlum/2.0;
  multiply_3d_farray(data,1,xn,1,yn,1,zn,contrast*mid);
  add_const_3d_farray(data,1,xn,1,yn,1,zn,mid);

  stimx_image_set(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL);  // Clean up

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_STIM_3D_PPL_BARNOISE                         */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Returned array dimensions start at 1.  To free returned storage:       */
/*    free_f3tensor(data,1,xn,1,yn,1,zn).                                    */
/*  - Do not reseed here, the caller must initialize the seed                */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_barnoise(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp)
     char *mylogf;
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;  // (deg/pix), (s/frame)
     int tsamp;
{
  int i,j,k;
  float ***data,t,maxlum,t0,tn,bgval,**tframe,blackval;

  //
  //  PREP
  //
  stimx_barnoise(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL,NULL);

  t0        = param_getf_exit(ppl,"st0");
  tn        = param_getf_exit(ppl,"stn");
  //contrast  = param_getf_exit(ppl,"contrast");
  maxlum    = param_getf_exit(ppl,"maxlum");
  bgval     = param_getf_exit(ppl,"bgval");

  // *** NOTE:  meanzero not needed here - bgval controls 0..1 vs. -1..1

  //bgamp = 2.0*bgval - 1.0; // Convert 0..1 --> -1..1
  //bgamp = bgval; // Convert 0..1 --> -1..1
  blackval = 2.0*bgval - maxlum;

  data = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;  // Time (s)
    if (((i-1)%tsamp) == 0){  // Is this a video frame
      if ((t >= t0) && (t < t0+tn)){
	
	stimx_barnoise(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,&tframe,
		       NULL);

	//
	//  Copy this frame in and free it [WYETH - CHANGE THIS]
	//
	for(j=0;j<xn;j++){
	  for(k=0;k<yn;k++){
	    data[j+1][k+1][i] = tframe[j][k];
	  }
	}
	free_2d_farray(tframe,xn);

      }else{
	make_frame_flat(data,1,xn,1,yn,1,zn,i,bgval);
      }
    }else{
      make_frame_flat(data,1,xn,1,yn,1,zn,i,blackval); // black
    }
  }

  // Assume all values in 'data' are in [-1,1]
  //mid = (float)maxlum/2.0;
  //multiply_3d_farray(data,1,xn,1,yn,1,zn,contrast*mid);
  //add_const_3d_farray(data,1,xn,1,yn,1,zn,mid);

  // Clean up
  stimx_barnoise(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL,NULL);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_STIM_3DB_PPL_BARNOISE                         */
/*                                                                           */
/*****************************************************************************/
void get_stim_3db_ppl_barnoise(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,rdl,rdr)
     char *mylogf;
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;  // (deg/pix), (s/frame)
     int tsamp;
     float ****rdl,****rdr;
{
  int i,j,k;
  float ***dl,***dr,t,maxlum,t0,tn,bgval,**tframl,**tframr;
  float blackval;

  //
  //  PREP
  //
  stimx_barnoise(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL,NULL);

  t0        = param_getf_exit(ppl,"st0");
  tn        = param_getf_exit(ppl,"stn");
  //contrast  = param_getf_exit(ppl,"contrast");
  maxlum    = param_getf_exit(ppl,"maxlum");
  bgval     = param_getf_exit(ppl,"bgval");

  // *** NOTE:  meanzero not needed here - bgval controls 0..1 vs. -1..1

  //bgamp = 2.0*bgval - 1.0; // Convert 0..1 --> -1..1
  //bgamp = bgval; // Convert 0..1 --> -1..1
  blackval = 2.0*bgval - maxlum;

  dl = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT
  dr = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;  // Time (s)
    if (((i-1)%tsamp) == 0){  // Is this a video frame
      if ((t >= t0) && (t < t0+tn)){

	stimx_barnoise(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,
		       &tframl,&tframr);

	//
	//  Copy this frame in and free it [WYETH - CHANGE THIS]
	//
	for(j=0;j<xn;j++){
	  for(k=0;k<yn;k++){
	    dl[j+1][k+1][i] = tframl[j][k];
	    dr[j+1][k+1][i] = tframr[j][k];
	  }
	}
	free_2d_farray(tframl,xn);
	free_2d_farray(tframr,xn);

      }else{
	make_frame_flat(dl,1,xn,1,yn,1,zn,i,bgval);
	make_frame_flat(dr,1,xn,1,yn,1,zn,i,bgval);
      }
    }else{
      make_frame_flat(dl,1,xn,1,yn,1,zn,i,blackval); // typ. black, -1.0
      make_frame_flat(dr,1,xn,1,yn,1,zn,i,blackval); // typ. black, -1.0
    }
  }

  // Assume all values in 'data' are in [-1,1]
  /*
  mid = (float)maxlum/2.0;
  multiply_3d_farray(dl,1,xn,1,yn,1,zn,contrast*mid);
  add_const_3d_farray(dl,1,xn,1,yn,1,zn,mid);
  multiply_3d_farray(dr,1,xn,1,yn,1,zn,contrast*mid);
  add_const_3d_farray(dr,1,xn,1,yn,1,zn,mid);
  */

  stimx_barnoise(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL,NULL);
  // Clean up

  *rdl = dl;
  *rdr = dr;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_STIM_3DB_PPL_BAR1RAN                         */
/*                                                                           */
/*****************************************************************************/
void get_stim_3db_ppl_bar1ran(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,rdl,rdr)
     char *mylogf;
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;  // (deg/pix), (s/frame)
     int tsamp;
     float ****rdl,****rdr;
{
  int i,j,k;
  float ***dl,***dr,t,t0,tn,bgval,maxlum,blackval,**tframl,**tframr;

  //
  //  PREP
  //
  stimx_bar1ran(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL,NULL);

  t0        = param_getf_exit(ppl,"st0");
  tn        = param_getf_exit(ppl,"stn");
  bgval     = param_getf_exit(ppl,"bgval");
  maxlum    = param_getf_exit(ppl,"maxlum");

  // *** NOTE:  meanzero not needed here - bgval controls 0..1 vs. -1..1

  //blackval = 0.0;   // Assume range is 0..1
  blackval = 2.0*bgval - maxlum;  // bg,ml = 0,1   ==> -1 
                                  // bg,ml = 0.5,1 ==>  0

  dl = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT
  dr = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;  // Time (s)
    if (((i-1)%tsamp) == 0){  // Is this a video frame
      if ((t >= t0) && (t < t0+tn)){

	stimx_bar1ran(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,
		      &tframl,&tframr);

	//
	//  Copy this frame in and free it [WYETH - CHANGE THIS]
	//
	for(j=0;j<xn;j++){
	  for(k=0;k<yn;k++){
	    dl[j+1][k+1][i] = tframl[j][k];
	    dr[j+1][k+1][i] = tframr[j][k];
	  }
	}
	free_2d_farray(tframl,xn);
	free_2d_farray(tframr,xn);

      }else{
	make_frame_flat(dl,1,xn,1,yn,1,zn,i,bgval);
	make_frame_flat(dr,1,xn,1,yn,1,zn,i,bgval);
      }
    }else{
      make_frame_flat(dl,1,xn,1,yn,1,zn,i,blackval);
      make_frame_flat(dr,1,xn,1,yn,1,zn,i,blackval);
    }
  }

  stimx_bar1ran(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL,NULL);
  // Clean up

  *rdl = dl;
  *rdr = dr;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_STIM_3D_PPL_BARPATMO                        */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Returned array dimensions start at 1.  To free returned storage:       */
/*    free_f3tensor(data,1,xn,1,yn,1,zn).                                    */
/*  - Do not reseed here, the caller must initialize the seed                */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_barpatmo(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i,j,k;
  float ***data,t,t0,tn,bgval,**tframe;

  //
  //  PREP
  //
  stimx_barpatmo(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL);

  t0        = param_getf_exit(ppl,"st0");
  tn        = param_getf_exit(ppl,"stn");
  bgval     = param_getf_exit(ppl,"bgval");

  data = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;  // Time (s)
    if (((i-1)%tsamp) == 0){  // Is this a video frame
      if ((t >= t0) && (t < t0+tn)){

	stimx_barpatmo(ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,&tframe);

	//
	//  Copy this frame in and free it [WYETH - CHANGE THIS]
	//
	for(j=0;j<xn;j++){
	  for(k=0;k<yn;k++){
	    data[j+1][k+1][i] = tframe[j][k];
	  }
	}
	free_2d_farray(tframe,xn);

      }else{
	make_frame_flat(data,1,xn,1,yn,1,zn,i,bgval);
      }
    }else{
      make_frame_flat(data,1,xn,1,yn,1,zn,i,0.0); // typ. black, 0.0
    }
  }

  stimx_barpatmo(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL);  // Clean up

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_STIM_3D_PPL_DOTS                          */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*    ddens - dotpairs per square degree per frame                           */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Returned array dimensions start at 1.  To free returned storage:       */
/*    free_f3tensor(data,1,xn,1,yn,1,zn).                                    */
/*  - Returned array has values within [minlum,maxlum].                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_dots(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp; // Time units/frame
{
  int i,j,k;
  float ***data,t,t0,tn,bgval,**tframe;

  //
  //  PREP
  //
  stimx_dots(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL);

  t0        = param_getf_exit(ppl,"st0");
  tn        = param_getf_exit(ppl,"stn");
  bgval     = param_getf_exit(ppl,"bgval");

  data = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  for(i=1;i<=zn;i++){
    if (((i-1)%tsamp) == 0){

      stimx_dots(ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,&tframe);

      //
      //  Copy this frame in and free it [WYETH - CHANGE THIS]
      //
      for(j=0;j<xn;j++){
	for(k=0;k<yn;k++){
	  data[j+1][k+1][i] = tframe[j][k];
	}
      }
      free_2d_farray(tframe,xn);

    }else
      make_frame_flat(data,1,xn,1,yn,1,zn,i,bgval);
  }

  stimx_dots(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL);  // Clean up

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_STIM_3D_PPL_DOT_PAIRED                        */
/*                                                                           */
/*  To implement the stimuli of Qian et al. (1994), where moving dots can    */
/*  be paired across two dot fields moving in independent directions.        */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Returned array dimensions start at 1.  To free returned storage:       */
/*    free_f3tensor(data,1,xn,1,yn,1,zn).                                    */
/*  - Returned array has values within [minlum,maxlum].                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_dot_paired(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp; // Time units/frame
{
  int i,j,k;
  float ***data,t,t0,tn,bgval,**tframe;

  //
  //  PREP
  //
  stimx_dot_paired(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL);

  t0        = param_getf_exit(ppl,"st0");
  tn        = param_getf_exit(ppl,"stn");
  bgval     = param_getf_exit(ppl,"bgval");

  data = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  for(i=1;i<=zn;i++){
    if (((i-1)%tsamp) == 0){

      stimx_dot_paired(ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,&tframe);

      //
      //  Copy this frame in and free it [WYETH - CHANGE THIS]
      //
      for(j=0;j<xn;j++){
	for(k=0;k<yn;k++){
	  data[j+1][k+1][i] = tframe[j][k];
	}
      }
      free_2d_farray(tframe,xn);

    }else
      make_frame_flat(data,1,xn,1,yn,1,zn,i,bgval);
  }

  stimx_dot_paired(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL);  // Clean up

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                       GET_STIM_3D_PPL_PASUPATHY_SHAPE                     */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Returned array dimensions start at 1.  To free returned storage:       */
/*    free_f3tensor(data,1,xn,1,yn,1,zn).                                    */
/*  - Returned array has values within [minlum,maxlum].                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_pasupathy_shape(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;               // Time units/frame
{
  int i,j,k;
  //float t0,tn;
  float ***data,t,bgval,**tframe;

  //
  //  PREP
  //
  stimx_pasupathy_shape(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL);

  //t0        = param_getf_exit(ppl,"st0");
  //tn        = param_getf_exit(ppl,"stn");
  bgval     = param_getf_exit(ppl,"bgval");

  data = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  for(i=1;i<=zn;i++){
    if (((i-1)%tsamp) == 0){
      stimx_pasupathy_shape(ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,&tframe);
      //
      //  Copy this frame in and free it [WYETH - CHANGE THIS]
      //
      for(j=0;j<xn;j++){
	for(k=0;k<yn;k++){
	  data[j+1][k+1][i] = tframe[j][k];
	}
      }
      free_2d_farray(tframe,xn);

    }else
      make_frame_flat(data,1,xn,1,yn,1,zn,i,bgval);
  }

  // Clean up
  stimx_pasupathy_shape(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_STIM_3D_PPL_BLUR_EDGE                        */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Returned array dimensions start at 1.  To free returned storage:       */
/*    free_f3tensor(data,1,xn,1,yn,1,zn).                                    */
/*  - Returned array has values within [minlum,maxlum].                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_blur_edge(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;               // Time units/frame
{
  int i,j,k;
  //float t0,tn;
  float ***data,t,bgval,**tframe;

//  exit_error("WYETH HERE","Under development (stim_util) blur_edge");

  //
  //  PREP
  //
  stimx_blur_edge(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL);

  //t0        = param_getf_exit(ppl,"st0");
  //tn        = param_getf_exit(ppl,"stn");
  bgval     = param_getf_exit(ppl,"bgval");

  data = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT

  for(i=1;i<=zn;i++){
    if (((i-1)%tsamp) == 0){
      stimx_blur_edge(ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,&tframe);
      //
      //  Copy this frame in and free it [WYETH - CHANGE THIS]
      //
      for(j=0;j<xn;j++){
	for(k=0;k<yn;k++){
	  data[j+1][k+1][i] = tframe[j][k];
	}
      }
      free_2d_farray(tframe,xn);
    }else
      make_frame_flat(data,1,xn,1,yn,1,zn,i,bgval);
  }

  // Clean up
  stimx_blur_edge(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_STIM_3D_PPL_FRAMESET                        */
/*                                                                           */
/*  Read stimulus movie from a data file.                                    */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_frameset(mylogf,s,ppl,xn,yn,zn,sscale,tscale,tsamp)
     char *mylogf;
     struct stim_struct *s;
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i,j,k;
  int txn,tyn,ttn,drange,fcode,binoc,colcode,i0,j0,fst_offx,fst_offy;
  int txn2,tyn2,ttn2,fcode2,binoc2,colcode2,fst_yrev;
  float ***data,***data2,fst_bgval;
  char *stim_file,*stim_file2,*path,*fst_op,*fst_form;
  char ggstr[SLEN];

  mylog(mylogf,"  GET_STIM_3D_PPL_FRAMESET\n");

  //
  //  Get up to 2 file names
  //
  stim_file  = param_getc_exit(ppl,"fst_file_1");
  stim_file2 = param_getc_dflt(ppl,"fst_file_2",NULL);
  fst_op     = param_getc_dflt(ppl,"fst_op","none");
  fst_bgval  = param_getf_dflt(ppl,"fst_bgval",0.0);
  fst_offx   = param_geti_dflt(ppl,"fst_pix_offset_x",-1);
  fst_offy   = param_geti_dflt(ppl,"fst_pix_offset_y",-1);
  fst_form   = param_getc_dflt(ppl,"fst_format","fst");
  fst_yrev   = param_geti_dflt(ppl,"fst_reverse_y",0);

  //
  //  Path to add to filename, or NULL
  //
  path = get_path_without_name(s->paramfile);

  drange = 0;

  if (strcmp(fst_form,"text")==0){

    // WYETH - 2019 - I am wondering how this routine compares to this
    //  other method of getting text files:
    // get_stim_3rgb_ppl_frameset(...)  <=== This allows for rotating RGB
    //  this '3rgb' format relates to the images that Dean supplied for
    //   AlexNet.


    // WYTEH - text files start with top line in image
    data = data_util_fst_read_txt(mylogf,path,stim_file,drange,fst_yrev,
				  &txn,&tyn,&ttn,&fcode,&binoc,&colcode);
  }else{
    // WYTEH - binary files start with bottom line in image
    data = data_util_fst_read(mylogf,path,stim_file,drange,&txn,&tyn,
			      &ttn,&fcode,&binoc,&colcode);
  }

  // *** TENSOR FORMAT  indexing starts at 1...

  if (stim_file2 != NULL){
    data2 = data_util_fst_read(mylogf,path,stim_file2,drange,&txn2,&tyn2,&ttn2,
			       &fcode2,&binoc2,&colcode2);
    // *** TENSOR FORMAT  indexing starts at 1...

    if (strcmp(fst_op,"mask")==0){
      stimu_fst_mask(data,txn,tyn,ttn,data2,txn2,tyn2,ttn2,fst_bgval);
    }else{
      exit_error("GET_STIM_3D_PPL_FRAMESET","Unknown 'fst_op' operation");
    }

    free_f3tensor(data2,1,txn2,1,tyn2,1,ttn2);
  }

  if ((txn > xn) || (tyn > yn) || (zn != ttn)){
    printf("  *** Model expects input xn,yn,zn = %d %d %d\n",xn,yn,zn);
    exit_error("GET_STIM_3D_PPL_FRAMESET","Stimulus size mismatch");
  }

  if (((txn < xn) && (tyn <= yn)) || ((tyn < yn) && (txn <= xn))){
    mylog(mylogf,"    Placing smaller stimulus within larger model field.\n");

    if (fst_offx == -1){ // Center
      i0 = xn - txn;
      if (i0 > 0)
	i0 /= 2;
    }else{
      i0 = fst_offx;
    }

    if (fst_offy == -1){ // Center
      j0 = yn - tyn;
      if (j0 > 0)
	j0 /= 2;
    }else{
      j0 = fst_offy;
    }

    //printf(" i0 = %d  j0 = %d\n",i0,j0);

    data2 = f3tensor(1,xn,1,yn,1,zn); // Use this for 3D FFT
    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){
	for(k=0;k<zn;k++){
	  data2[i+1][j+1][k+1] = fst_bgval;
	}
      }
    }
    for(i=0;i<txn;i++){
      for(j=0;j<tyn;j++){
	for(k=0;k<ttn;k++){
	  data2[i+1+i0][j+1+j0][k+1] = data[i+1][j+1][k+1];
	}
      }
    }

    free_f3tensor(data,1,txn,1,tyn,1,ttn);
    data = data2;
    
  }else{
    if ((fst_offx != -1) || (fst_offy != -1)){
      exit_error("GET_STIM_3D_PPL_FRAMESET",
		 "FST pixel offset not valid for this stimulus size.");
    }
  }

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_STIM_3D_PPL_UPLOAD                          */
/*                                                                           */
/*  Upload the stimulus data from a file.                                    */
/*                                                                           */
/*    sscale - spatial scale (deg/pixel)                                     */
/*    tscale - temporal scale (sec/frame)                                    */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl_upload(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i,k;
  int dxn,dyn,dtn,nbytes,tcode,winx0,winy0,winxn,winyn,destx0,desty0;
  float ***data,***d,bgval,maxlum,t0,stim_samp,mod_samp,t,stim_max;
  char *data_format,*stim_path,*stim_file,infile[SLEN];

  data_format  =  paramfile_get_char_param_or_exit(ppl,"data_format");
  stim_path    =  paramfile_get_char_param_or_exit(ppl,"stim_path");
  stim_file    =  paramfile_get_char_param_or_exit(ppl,"stim_file");
  stim_max     = param_getf_exit(ppl,"stim_max");
  bgval        = param_getf_exit(ppl,"bgval");
  t0           = param_getf_exit(ppl,"t0");
  stim_samp    = param_getf_exit(ppl,"stim_samp");

  mod_samp = my_rint(1.0/tscale);
  if (stim_samp > mod_samp){
    printf("    stim_samp = %f   >   tscale = %f\n",stim_samp,tscale);
    exit_error("GET_STIM_3D_PPL_UPLOAD","stim_samp is too high");
  }

  sprintf(infile,"%s/%s",stim_path,stim_file);
  if (strcmp(data_format,"3d")==0){
    read_3d_data(infile,&d,&dxn,&dyn,&dtn,&nbytes,&tcode,1);
    if ((nbytes != 4) || (tcode != 2))
      exit_error("GET_STIM_3D_PPL_UPLOAD","Expecting 4 byte float format");
  }else
    exit_error("GET_STIM_3D_PPL_UPLOAD","Unknown data format");
  
  myfree(data_format); myfree(stim_path); myfree(stim_file);

  if ((dxn != xn) || (dyn != yn)){
    printf("  Data:  %d X %d\n",dxn,dyn);
    exit_error("GET_STIM_3D_PPL_UPLOAD","Spatial size does not match model");
  }

  /* Scale the raw data */
  multiply_3d_farray(d,0,dxn,0,dyn,0,dtn,1.0/stim_max);

  /* Be sure the data is in [0..1] */
  /*stim_util_limit_3d(d,0,dxn,0,dyn,0,dtn,0.0,1.0);*/

  data = f3tensor(1,xn,1,yn,1,zn); /*** USE THIS FOR 3D FFT ***/

  maxlum = 1.0;           /* For now, fixed */
  destx0 = desty0 = 1;
  winx0 = winy0 = 0;
  winxn = xn;
  winyn = yn;

  k = 0; // Count frames in source data
  for(i=1;i<=zn;i++){
    t = (float)(i-1)*tscale;
    if (((i-1)%tsamp) == 0){
      if ((t >= t0) && (k < dtn)){  // Ends if source runs out, tn not used
	make_frame_from_3d(data,1,xn,1,yn,1,zn,i,d,dxn,dyn,dtn,k,winx0,winy0,
			   winxn,winyn,destx0,desty0,bgval,maxlum);
	k += 1;
      }else
	make_frame_flat(data,1,xn,1,yn,1,zn,i,bgval);
    }else
      make_frame_flat(data,1,xn,1,yn,1,zn,i,-1.0);
  }

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_STIM_3C_PPL_SINE                           */
/*                                                                           */
/*  Returned array dimensions [0..xn-1][0..yn-1][0..zn-1]                    */
/*                                                                           */
/*****************************************************************************/
int ***get_stim_3c_ppl_sine(mylogf,ppl,xn,yn,zn,sscale,tscale,tsamp)
     char *mylogf;
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int ***data;
  float ***fdata,***smooth;
  float r0,g0,b0,r1,g1,b1;

  fdata = get_stim_3d_ppl_sine(ppl,xn,yn,zn,sscale,tscale,tsamp);

  smooth = stim_util_blur(mylogf,ppl,fdata,xn,yn,zn,sscale,tscale);
  if (smooth != NULL){
    free_f3tensor(fdata,1,xn,1,yn,1,zn);
    fdata = smooth;
  }

  r0 = param_getf_exit(ppl,"color0_r");
  g0 = param_getf_exit(ppl,"color0_g");
  b0 = param_getf_exit(ppl,"color0_b");

  r1 = param_getf_exit(ppl,"color1_r");
  g1 = param_getf_exit(ppl,"color1_g");
  b1 = param_getf_exit(ppl,"color1_b");

  //printf("    color 0 rgb  %.4f %.4f %.4f\n",r0,g0,b0);
  //printf("    color 1 rgb  %.4f %.4f %.4f\n",r1,g1,b1);
  //printf("    color 2 rgb  %.4f %.4f %.4f\n",r2,g2,b2);

  //data = stim_util_3c_from_tern_3d(fdata,1,xn,1,yn,1,zn,r0,g0,b0,r1,g1,b1,
  //r2,g2,b2);

  data = stim_util_3c_2col_from_3d(fdata,1,xn,1,yn,1,zn,r0,g0,b0,r1,g1,b1);

  free_f3tensor(fdata,1,xn,1,yn,1,zn);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_STIM_3C_PPL_NOISE                          */
/*                                                                           */
/*  Returned array dimensions [0..xn-1][0..yn-1][0..zn-1]                    */
/*                                                                           */
/*****************************************************************************/
int ***get_stim_3c_ppl_noise(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int ***data,ncol;
  float ***fdata;
  float r0,g0,b0,r1,g1,b1,r2,g2,b2,r3,g3,b3;
  float bg_r,bg_g,bg_b,bgval;

  fdata = get_stim_3d_ppl_noise(ppl,xn,yn,zn,sscale,tscale,tsamp);

  r0 = param_getf_exit(ppl,"color0_r");
  g0 = param_getf_exit(ppl,"color0_g");
  b0 = param_getf_exit(ppl,"color0_b");

  r1 = param_getf_exit(ppl,"color1_r");
  g1 = param_getf_exit(ppl,"color1_g");
  b1 = param_getf_exit(ppl,"color1_b");


  if (paramfile_test_param(ppl,"bg_r")){ // DMVAR
    bg_r = param_getf_exit(ppl,"bg_r");
    bg_g = param_getf_exit(ppl,"bg_g");
    bg_b = param_getf_exit(ppl,"bg_b");
  }else{
    bgval = param_getf_exit(ppl,"bgval"); // 0..1
    bg_r = bg_g = bg_b = bgval;
  }

  ncol = param_geti_dflt(ppl,"n_color",2);
  if (ncol == 4){
    r2 = param_getf_exit(ppl,"color2_r");
    g2 = param_getf_exit(ppl,"color2_g");
    b2 = param_getf_exit(ppl,"color2_b");

    r3 = param_getf_exit(ppl,"color3_r");
    g3 = param_getf_exit(ppl,"color3_g");
    b3 = param_getf_exit(ppl,"color3_b");

    data = stim_util_3c_4col_from_3d(fdata,1,xn,1,yn,1,zn,r0,g0,b0,r1,g1,b1,
				     r2,g2,b2,r3,g3,b3,bg_r,bg_g,bg_b);
  }else{
    //
    //  2 color values
    //
    data = stim_util_3c_2col_from_3d(fdata,1,xn,1,yn,1,zn,r0,g0,b0,r1,g1,b1);
  }

  //printf("    color 0 rgb  %.4f %.4f %.4f\n",r0,g0,b0);
  //printf("    color 1 rgb  %.4f %.4f %.4f\n",r1,g1,b1);


  free_f3tensor(fdata,1,xn,1,yn,1,zn);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_STIM_3C_PPL_SURRT                          */
/*                                                                           */
/*  Returned array dimensions [0..xn-1][0..yn-1][0..zn-1]                    */
/*                                                                           */
/*****************************************************************************/
int ***get_stim_3c_ppl_surrt(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int ***data;
  float ***fdata;
  float r0,g0,b0,r1,g1,b1;

  fdata = get_stim_3d_ppl_surrt(ppl,xn,yn,zn,sscale,tscale,tsamp);

  r0 = param_getf_exit(ppl,"color0_r");
  g0 = param_getf_exit(ppl,"color0_g");
  b0 = param_getf_exit(ppl,"color0_b");

  r1 = param_getf_exit(ppl,"color1_r");
  g1 = param_getf_exit(ppl,"color1_g");
  b1 = param_getf_exit(ppl,"color1_b");

  //printf("    color 0 rgb  %.4f %.4f %.4f\n",r0,g0,b0);
  //printf("    color 1 rgb  %.4f %.4f %.4f\n",r1,g1,b1);

  data = stim_util_3c_2col_from_3d(fdata,1,xn,1,yn,1,zn,r0,g0,b0,r1,g1,b1);

  free_f3tensor(fdata,1,xn,1,yn,1,zn);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_STIM_3C_PPL_RGB_SINE                         */
/*                                                                           */
/*  Returned array dimensions [0..xn-1][0..yn-1][0..zn-1]                    */
/*                                                                           */
/*****************************************************************************/
int ***get_stim_3c_ppl_rgb_sine(ppl,xn,yn,zn,sscale,tscale,tsamp)
     struct param_pair_list *ppl;
     int xn,yn,zn;
     float sscale,tscale;
     int tsamp;
{
  int i,j,k;
  int ***data,**tframe,bgval;
  float r0,g0,b0,sv;


  //
  //  WYETH - Jul 2014 - Make the stimulus of 
  //  PDF #60007:  Sun, Smithson, Zaidi, Lee (2006)
  //

  r0        = param_getf_exit(ppl,"r_offset");
  g0        = param_getf_exit(ppl,"g_offset");
  b0        = param_getf_exit(ppl,"b_offset");
  sv        = param_getf_exit(ppl,"c_scale");

  bgval = data_util_t11_int_for_rgb(sv*r0,sv*g0,sv*b0);

  //
  //  PREP
  //
  stimx_rgb_sine(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,0,NULL);


  data = get_3d_iarray(xn,yn,zn);

  for(i=0;i<zn;i++){
    if ((i%tsamp) == 0){

      stimx_rgb_sine(ppl,xn,yn,zn,sscale,tscale,tsamp,i-1,1,&tframe);

      //
      //  Copy this frame in and free it [WYETH - CHANGE THIS]
      //
      for(j=0;j<xn;j++){
	for(k=0;k<yn;k++){
	  data[j][k][i] = tframe[j][k];
	}
      }
      free_2d_iarray(tframe,xn);

    }else{
      //
      //  make frame flat for integer
      //
      for(j=0;j<xn;j++){
	for(k=0;k<yn;k++){
	  data[j][k][i] = bgval;
	}
      }
    }
  }

  // Clean up
  stimx_rgb_sine(ppl,xn,yn,zn,sscale,tscale,tsamp,-1,-1,NULL);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_STIM_3RGB_PPL_FRAMESET                       */
/*                                                                           */
/*  Returned array dimensions [0,1,2][0..xn-1][0..yn-1][0..tn-1]             */
/*                                                                           */
/*****************************************************************************/
float ****get_stim_3rgb_ppl_frameset(mylogf,s,ppl,xn,yn,tn,sscale,tscale,tsamp)
     char *mylogf;
     struct stim_struct *s;
     struct param_pair_list *ppl;
     int xn,yn,tn;
     float sscale,tscale;
     int tsamp;
{
  FILE *fopen(),*fin;
  int i,j,k,it;
  int ns,t,dxn,dyn,dzn,dtn,rgb_ord,ii[3];
  float ****data;
  char *dformat,*infile,ts1[SLEN],ts2[SLEN];

  dformat    = param_getc_exit(ppl,"fst_format");

  if (strcmp(dformat,"text")!=0){
    exit_error("GET_STIM_3RGB_PPL_FRAMESET",
	       "'frameset_format' must be 'text'");
  }
  infile  = param_getc_exit(ppl,"fst_file_1");  // Name of image data file

  rgb_ord = param_geti_dflt(ppl,"fst_rgb_order",0);


  //
  //  WYETH - This should probably be a utility to be called by others...
  //
  if ((fin = fopen(infile,"r"))==NULL){
    printf("  File:  %s\n",infile);
    exit_error("GET_STIM_3RGB_PPL_FRAMESET","Cannot open file");
  }

  ns = fscanf(fin,"%s %s",ts1,ts2);

  if ((strcmp(ts1,"3")==0) && (strcmp(ts2,"227")==0)){
    //
    //  Wyeth - this is a hack for Dean's simpler file format
    //
    dtn = 1;
    dzn = 3;
    dxn = 227;
    ns = fscanf(fin,"%d",&dyn);  // Must be 227
  }else{
    if (strcmp(ts1,"frameset")!=0){
      exit_error("GET_STIM_3RGB_PPL_FRAMESET",
		 "First string in file must be 'frameset'");
    }
    if (strcmp(ts2,"text")!=0){
      exit_error("GET_STIM_3RGB_PPL_FRAMESET",
		 "Second string in file must be 'text'");
    }
    ns = fscanf(fin,"%d %d %d %d",&dtn,&dzn,&dxn,&dyn);
  }
  //printf("    Image size:  %d %d %d\n",dxn,dyn,dzn);
  //printf("    Number of images:  %d\n",dtn);

  if ((dxn != xn) || (dyn != yn)){
    printf("  Expecting %d x %d images, but found %d x %d\n",xn,yn,dxn,dyn);
    exit_error("GET_STIM_3RGB_PPL_FRAMESET","Image size mis-match");
  }

  if (dtn != tn){
    printf("  tn = %d  dtn = %d\n",tn,dtn);
    exit_error("GET_STIM_3RGB_PPL_FRAMESET","Number of images mis-match");
  }

  if (dzn != 3)
    exit_error("GET_STIM_3RGB_PPL_FRAMESET","Image depth mis-match");

  data = get_4d_farray(dzn,xn,yn,tn);

  if       (rgb_ord == 0){  ii[0] = 0;  ii[1] = 1;  ii[2] = 2;
  }else if (rgb_ord == 1){  ii[0] = 0;  ii[1] = 2;  ii[2] = 1;
  }else if (rgb_ord == 2){  ii[0] = 1;  ii[1] = 0;  ii[2] = 2;
  }else if (rgb_ord == 3){  ii[0] = 1;  ii[1] = 2;  ii[2] = 0;
  }else if (rgb_ord == 4){  ii[0] = 2;  ii[1] = 0;  ii[2] = 1;
  }else if (rgb_ord == 5){  ii[0] = 2;  ii[1] = 1;  ii[2] = 0;
  }else{
    exit_error("GET_STIM_3RGB_PPL_FRAMESET","Invalid 'rgb_order' value");
  }

  // *** NOTE ***
  //     The *BLUE* image data comes first, starting from upper left pixel
  //     then the *GREEN* image data comes next
  //     then the *RED* image data comes last.
  //     This is the case in Dean's simple file format.

  for(t=0;t<tn;t++){            // For each image
    for(it=0;it<dzn;it++){      //   For R,G,B planes
      i = ii[it];               //     (allow rotation of RGB)
      for(k=(yn-1);k>=0;k--){   //     For each row (starting with top row),
	for(j=0;j<xn;j++){      //       read from left to right.
	  ns = fscanf(fin,"%f",&(data[i][j][k][t]));
	}
      }
    }
  }
  fclose(fin);

  myfree(infile);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                              GET_STIM_3D_PPL                              */
/*                                                                           */
/*  NOTE: Returned array dimensions start at 1.  To free returned            */
/*        storage:  free_f3tensor(data,1,xn,1,yn,1,zn).                      */
/*                                                                           */
/*****************************************************************************/
float ***get_stim_3d_ppl(mylogf,s,xn,yn,tn,sscale,tscale,tsamp)
     char *mylogf;
     struct stim_struct *s;  // Stimulus params
     int xn,yn,tn;
     float sscale,tscale;
     int tsamp;
{
  float ***d;
  char ggstr[SLEN];
  struct param_pair_list *sppl;
  char *stimtype;

  sppl = s->ppl;
  stimtype = param_getc_exit(sppl,"stim_type");

  if (strcmp(stimtype,"adamk")==0){
    d = get_stim_3d_ppl_adamk(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"disk")==0){
    d = get_stim_3d_ppl_disk(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"sine")==0){
    d = get_stim_3d_ppl_sine(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"sine_seq")==0){
    d = get_stim_3d_ppl_sine_seq(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"sos")==0){
    d = get_stim_3d_ppl_sos(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"pixel")==0){
    d = get_stim_3d_ppl_pixel(mylogf,sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"bar")==0){
    d = get_stim_3d_ppl_bar(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"3barlink")==0){
    d = get_stim_3d_ppl_3barlink(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"barnoise")==0){
    d = get_stim_3d_ppl_barnoise(mylogf,sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"barpatmo")==0){
    d = get_stim_3d_ppl_barpatmo(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"barray")==0){
    d = get_stim_3d_ppl_barray(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"barray_old")==0){
    d = get_stim_3d_ppl_barray_old(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"bar2ran")==0){
    d = get_stim_3d_ppl_bar2ran(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"cpsine")==0){
    d = get_stim_3d_ppl_cpsine(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"plaid")==0){
    d = get_stim_3d_ppl_plaid(mylogf,sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"dmmask")==0){
    d = get_stim_3d_ppl_dmmask(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"dirmod")==0){
    d = get_stim_3d_ppl_dirmod(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"surrt")==0){
    d = get_stim_3d_ppl_surrt(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"wy4st")==0){
    d = get_stim_3d_ppl_wy4st(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"x4st")==0){
    d = get_stim_3d_ppl_x4st(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"wy2pa")==0){
    d = get_stim_3d_ppl_wy2pa(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"grid")==0){
    exit_error("GET_STIM_3D_PPL","'grid' has been changed to 'noise'");
  }else if (strcmp(stimtype,"noise")==0){
    d = get_stim_3d_ppl_noise(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"gp")==0){
    d = get_stim_3d_ppl_gp(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"ransinphase")==0){
    d = get_stim_3d_ppl_ransinphase(sppl,xn,yn,tn,sscale,tscale,tsamp,0);
  }else if (strcmp(stimtype,"ransinori")==0){
    d = get_stim_3d_ppl_ransinphase(sppl,xn,yn,tn,sscale,tscale,tsamp,1);
  }else if (strcmp(stimtype,"ransindir")==0){
    d = get_stim_3d_ppl_ransinphase(sppl,xn,yn,tn,sscale,tscale,tsamp,2);
  }else if (strcmp(stimtype,"ran_sin")==0){
    d = get_stim_3d_ppl_ran_sin(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"ranstep")==0){
    d = get_stim_3d_ppl_ranstep(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"test_shape")==0){
    d = get_stim_3d_ppl_test_shape(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"dots")==0){
    d = get_stim_3d_ppl_dots(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"dot_paired")==0){
    d = get_stim_3d_ppl_dot_paired(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"dynedge")==0){
    d = get_stim_3d_ppl_dynedge(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"pasupathy_shape")==0){
    d = get_stim_3d_ppl_pasupathy_shape(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"blur_edge")==0){
    d = get_stim_3d_ppl_blur_edge(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"gabor_britt")==0){
    d = get_stim_3d_ppl_gabor_britten(mylogf,sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"globmo_patch")==0){
    get_stim_3db_ppl_globmo_patch(mylogf,sppl,xn,yn,tn,sscale,tscale,tsamp,
				  &d,NULL);
  }else if (strcmp(stimtype,"upload")==0){
    // WYETH - this was an older upload method - do we still need it?
    // WYETH - this was an older upload method - do we still need it?
    d = get_stim_3d_ppl_upload(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"frameset")==0){
    d = get_stim_3d_ppl_frameset(mylogf,s,sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"image_set")==0){
    //system("date");
    mylog(mylogf,"    Constructing stimulus from image set ...\n");
    d = get_stim_3d_ppl_image_set(sppl,xn,yn,tn,sscale,tscale,tsamp);
    //system("date");
    mylog(mylogf,"    ... done.\n");
  }else if (strcmp(stimtype,"null")==0){
    d = NULL;
  }else{
    sprintf(ggstr,"stimtype %s\n",stimtype);
    mylog(mylogf,ggstr);
    mylogx(mylogf,"GET_STIM_3D_PPL","Unknown stimulus type");
  }

  myfree(stimtype);

  return d;
}
/**************************************-**************************************/
/*                                                                           */
/*                              GET_STIM_3C_PPL                              */
/*                                                                           */
/*  Returned arrays [0..xn-1][0..yn-1][0..tn-1]                              */
/*                                                                           */
/*****************************************************************************/
int ***get_stim_3c_ppl(mylogf,stimtype,sppl,xn,yn,tn,sscale,tscale,tsamp)
     char *mylogf;
     char *stimtype;
     struct param_pair_list *sppl;
     int xn,yn,tn;
     float sscale,tscale;
     int tsamp;
{
  int ***d;

  if (strcmp(stimtype,"sine")==0){
    d = get_stim_3c_ppl_sine(mylogf,sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"grid")==0){
    exit_error("GET_STIM_3C_PPL","'grid' has been changed to 'noise'");
  }else if (strcmp(stimtype,"noise")==0){
    d = get_stim_3c_ppl_noise(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"surrt")==0){
    d = get_stim_3c_ppl_surrt(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else if (strcmp(stimtype,"rgb_sine")==0){
    d = get_stim_3c_ppl_rgb_sine(sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else{
    printf("stimtype %s\n",stimtype);
    exit_error("GET_STIM_3C_PPL","Unknown stimulus type");
  }

  return d;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_STIM_3RGB_PPL                             */
/*                                                                           */
/*  Returned arrays [0..xn-1][0..yn-1][0..tn-1][0,1,2]                       */
/*                                                                           */
/*****************************************************************************/
float ****get_stim_3rgb_ppl(mylogf,s,stimtype,sppl,xn,yn,tn,sscale,tscale,
			    tsamp)
     char *mylogf;
     struct stim_struct *s;  // Stimulus params
     char *stimtype;
     struct param_pair_list *sppl;
     int xn,yn,tn;
     float sscale,tscale;
     int tsamp;
{
  float ****d;

  //  WYETH NOTE, this is called only by this:
  //    void mod_util_get_3rgb_stim(s,k,xn,yn,tn,sscale,tscale,tsamp,mylogf)
  //  which gets and can smooth the stimulus, and which is called only by
  //    wm_util.c in 'get_stimulus_data_tuning_curve' ...


  if (strcmp(stimtype,"frameset")==0){
    d = get_stim_3rgb_ppl_frameset(mylogf,s,sppl,xn,yn,tn,sscale,tscale,tsamp);
  }else{
    printf("stimtype %s\n",stimtype);
    exit_error("GET_STIM_3RGB_PPL","Unknown stimulus type");
  }

  return d;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_STIM_3D_B_PPL                             */
/*                                                                           */
/*  Binocular stimuli.                                                       */
/*                                                                           */
/*  NOTES                                                                    */
/*  - Returned array dimensions start at 1                                   */
/*  - To free returned storage:  free_f3tensor(data,1,xn,1,yn,1,zn)          */
/*                                                                           */
/*****************************************************************************/
void get_stim_3d_b_ppl(mylogf,stimtype,sppl,xn,yn,tn,sscale,tscale,tsamp,
		       rdl,rdr)  // Return left and right stimuli
     char *mylogf;
     char *stimtype;
     struct param_pair_list *sppl;
     int xn,yn,tn;
     float sscale,tscale;
     int tsamp;
     float ****rdl,****rdr;
{
  float ***dr,***dl;

  if (strcmp(stimtype,"dot_stereo")==0){
    dl = (float ***)NULL;
    dr = (float ***)NULL;
  }else if (strcmp(stimtype,"sine")==0){
    get_stim_3db_ppl_sine(sppl,xn,yn,tn,sscale,tscale,tsamp,&dl,&dr);
  }else if (strcmp(stimtype,"randot_stereogram")==0){
    get_stim_3db_ppl_randot_stereogram(sppl,xn,yn,tn,sscale,tscale,tsamp,
				       &dl,&dr);
  }else if (strcmp(stimtype,"rds_disp_mod")==0){
    get_stim_3db_ppl_randot_stereogram(sppl,xn,yn,tn,sscale,tscale,tsamp,
				       &dl,&dr);
  }else if (strcmp(stimtype,"rds_mid_mod")==0){
    get_stim_3db_ppl_randot_stereogram(sppl,xn,yn,tn,sscale,tscale,tsamp,
				       &dl,&dr);
  }else if (strcmp(stimtype,"rds_cd_mod")==0){
    get_stim_3db_ppl_randot_stereogram(sppl,xn,yn,tn,sscale,tscale,tsamp,
				       &dl,&dr);
  }else if (strcmp(stimtype,"barnoise")==0){
    get_stim_3db_ppl_barnoise(mylogf,sppl,xn,yn,tn,sscale,tscale,tsamp,
			      &dl,&dr);
  }else if (strcmp(stimtype,"bar1ran")==0){
    get_stim_3db_ppl_bar1ran(mylogf,sppl,xn,yn,tn,sscale,tscale,tsamp,
			       &dl,&dr);
  }else if (strcmp(stimtype,"bar")==0){
    get_stim_3db_ppl_bar(mylogf,sppl,xn,yn,tn,sscale,tscale,tsamp,
			 &dl,&dr);
  }else if (strcmp(stimtype,"plaid")==0){
    get_stim_3db_ppl_plaid(mylogf,sppl,xn,yn,tn,sscale,tscale,tsamp,&dl,&dr);
  }else if (strcmp(stimtype,"gabor_grid")==0){
    get_stim_3db_ppl_gabor_grid(mylogf,sppl,xn,yn,tn,sscale,tscale,tsamp,
				&dl,&dr);
  }else if (strcmp(stimtype,"globmo_patch")==0){
    get_stim_3db_ppl_globmo_patch(mylogf,sppl,xn,yn,tn,sscale,tscale,tsamp,
				&dl,&dr);
  }else if (strcmp(stimtype,"dirmod")==0){
    get_stim_3db_ppl_dirmod(mylogf,sppl,xn,yn,tn,sscale,tscale,tsamp,&dl,&dr);
  }else{
    printf("stimtype %s\n",stimtype);
    exit_error("GET_STIM_3D_B_PPL","Unknown stimulus type");
  }

  *rdl = dl;
  *rdr = dr;
}
/**************************************-**************************************/
/*                                                                           */
/*                             STIM_GET_2D_IMPULSE                           */
/*                                                                           */
/*****************************************************************************/
float **stim_get_2d_impulse(xn,yn,x0,y0,ampl0,ampl)
     int xn,yn,x0,y0;
     float ampl0,ampl;
{
  float **s;

  s = get_zero_2d_farray(xn,yn);
  add_const_2d_farray(s,0,xn,0,yn,ampl0);
  s[x0][y0] = ampl;

  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                              STIM_MAKE_XT_IMPULSE                         */
/*                                                                           */
/*****************************************************************************/
float **stim_make_xt_impulse(p)
     struct param_pair_list *p;
{
  int xn,tn,x0,t0;
  float **s,stim_null,ampl;

  xn = paramfile_get_int_param_or_exit(p,"xn");
  tn = paramfile_get_int_param_or_exit(p,"tn");
  x0 = paramfile_get_int_param_or_exit(p,"stim_x0");
  t0 = paramfile_get_int_param_or_exit(p,"stim_t0");
  stim_null = param_getf_exit(p,"stim_null");
  ampl = param_getf_exit(p,"stim_ampl");

  s = stim_get_2d_impulse(xn,tn,x0,t0,stim_null,ampl);

  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                        STIM_MAKE_XT_DRIFTING_GRATING                      */
/*                                                                           */
/*****************************************************************************/
float **stim_make_xt_drifting_grating(p)
     struct param_pair_list *p;
{
  int i,j;
  int xn,tn,x0,stim_xn,t0,stim_tn,dt,dir;
  float tt,**s,sscale,tscale,sf,tf,phase,stim_null,ampl,ampl0;

  xn = paramfile_get_int_param_or_exit(p,"xn");
  tn = paramfile_get_int_param_or_exit(p,"tn");
  sscale = param_getf_exit(p,"sscale");
  tscale = param_getf_exit(p,"tscale");
  x0 = paramfile_get_int_param_or_exit(p,"stim_x0");
  stim_xn = paramfile_get_int_param_or_exit(p,"stim_xn");
  t0 = paramfile_get_int_param_or_exit(p,"stim_t0");
  stim_tn = paramfile_get_int_param_or_exit(p,"stim_tn");
  dt = paramfile_get_int_param_or_exit(p,"stim_dt");
  sf = param_getf_exit(p,"stim_sf");       /* cyc/deg */
  tf = param_getf_exit(p,"stim_tf");       /* Hz */
  dir = paramfile_get_int_param_or_exit(p,"stim_dir");     /* +/- 1 */
  phase = param_getf_exit(p,"stim_phase"); /* Degrees */
  stim_null = param_getf_exit(p,"stim_null");
  ampl0 = param_getf_exit(p,"stim_ampl0");
  ampl = param_getf_exit(p,"stim_ampl");

  phase /= 360.0; /* Change phase to cycles */
  s = (float **)myalloc(tn*sizeof(float *));
  tf *= (float)dir;
  for(i=0;i<tn;i++){
    s[i] = (float *)myalloc(xn*sizeof(float));
    if ((i%dt) == 0){
      if ((i > t0)&&(i < t0+stim_tn)){
	tt = (float)(i-t0)*tscale*tf;
	for(j=0;j<xn;j++){
	  if ((j >= x0)&&(j < x0+stim_xn))
	    s[i][j] = ampl0 + ampl*sin(((float)(j-x0)*sscale*sf - phase - tt) *
				       2.0*M_PI);
	  else
	    s[i][j] = ampl0;
	}
      }else{
	for(j=0;j<xn;j++)
	  s[i][j] = ampl0;
      }
    }else{
      for(j=0;j<xn;j++)
	s[i][j] = stim_null; /* Fill in between-frame value. */
    }
  }
  transpose_2d_farray(&s,tn,xn);
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                         STIM_MAKE_XT_JUMPING_GRATING                      */
/*                                                                           */
/*****************************************************************************/
float **stim_make_xt_jumping_grating(p)
     struct param_pair_list *p;
{
  int i,j;
  int xn,tn,x0,stim_xn,t0,stim_tn,dt,dir,seed;
  float **s,sscale,tscale,sf,tf,phase,stim_null,ampl,ampl0,delta;

  xn = paramfile_get_int_param_or_exit(p,"xn");
  tn = paramfile_get_int_param_or_exit(p,"tn");
  sscale = param_getf_exit(p,"sscale");
  tscale = param_getf_exit(p,"tscale");
  x0 = paramfile_get_int_param_or_exit(p,"stim_x0");
  stim_xn = paramfile_get_int_param_or_exit(p,"stim_xn");
  t0 = paramfile_get_int_param_or_exit(p,"stim_t0");
  stim_tn = paramfile_get_int_param_or_exit(p,"stim_tn");
  dt = paramfile_get_int_param_or_exit(p,"stim_dt");
  sf = param_getf_exit(p,"stim_sf");       /* cyc/deg */
  tf = param_getf_exit(p,"stim_tf");       /* Hz */
  dir = paramfile_get_int_param_or_exit(p,"stim_dir");     /* +/- 1 */
  phase = param_getf_exit(p,"stim_phase"); /* Degrees */
  stim_null = param_getf_exit(p,"stim_null");
  ampl0 = param_getf_exit(p,"stim_ampl0");
  ampl = param_getf_exit(p,"stim_ampl");
  seed = paramfile_get_int_param_or_exit(p,"stim_jumpseed");
  if (seed > 0)
    seed *= -1;

  delta = tf * tscale * (float)dt;
  printf("delta = %f\n",delta);
  phase /= 360.0; /* Change phase to cycles */
  s = (float **)myalloc(tn*sizeof(float *));
  tf *= (float)dir;
  for(i=0;i<tn;i++){
    s[i] = (float *)myalloc(xn*sizeof(float));
    if ((i%dt) == 0)
      if ((i > t0)&&(i < t0+stim_tn)){
	if (myrand_util_ran2(&seed) > 0.5)
	  phase -= delta;
	else
	  phase += delta;
	for(j=0;j<xn;j++)
	  if ((j >= x0)&&(j < x0+stim_xn))
	    s[i][j] = ampl0 + ampl*sin(((float)(j-x0)*sscale*sf - phase)*
				       2.0*M_PI);
	  else
	    s[i][j] = ampl0;
      }else
	for(j=0;j<xn;j++)
	  s[i][j] = ampl0;
    else
      for(j=0;j<xn;j++)
	s[i][j] = stim_null; /* Fill in between-frame value. */
  }
  transpose_2d_farray(&s,tn,xn);
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                           STIM_MAKE_XT_AMC_GRATING                        */
/*                                                                           */
/*  Modulate the contrast of the amplitude according to a random walk.       */
/*                                                                           */
/*****************************************************************************/
float **stim_make_xt_amc_grating(p)
     struct param_pair_list *p;
{
  int i,j;
  int xn,tn,x0,stim_xn,t0,stim_tn,dt,dir,seed,curramp,step,maxamp;
  float tt,**s,sscale,tscale,sf,tf,phase,stim_null,ampl,ampl0;

  xn = paramfile_get_int_param_or_exit(p,"xn");
  tn = paramfile_get_int_param_or_exit(p,"tn");
  sscale = param_getf_exit(p,"sscale");
  tscale = param_getf_exit(p,"tscale");
  x0 = paramfile_get_int_param_or_exit(p,"stim_x0");
  stim_xn = paramfile_get_int_param_or_exit(p,"stim_xn");
  t0 = paramfile_get_int_param_or_exit(p,"stim_t0");
  stim_tn = paramfile_get_int_param_or_exit(p,"stim_tn");
  dt = paramfile_get_int_param_or_exit(p,"stim_dt");
  sf = param_getf_exit(p,"stim_sf");       /* cyc/deg */
  tf = param_getf_exit(p,"stim_tf");       /* Hz */
  dir = paramfile_get_int_param_or_exit(p,"stim_dir");     /* +/- 1 */
  phase = param_getf_exit(p,"stim_phase"); /* Degrees */
  stim_null = param_getf_exit(p,"stim_null");
  ampl0 = param_getf_exit(p,"stim_ampl0");
  seed = paramfile_get_int_param_or_exit(p,"stim_jumpseed");
  step = paramfile_get_int_param_or_exit(p,"stim_step");
  maxamp = paramfile_get_int_param_or_exit(p,"stim_maxamp");
  if (seed > 0)
    seed *= -1;

  phase /= 360.0; /* Change phase to cycles */
  s = (float **)myalloc(tn*sizeof(float *));
  tf *= (float)dir;
  curramp = 0;
  for(i=0;i<tn;i++){
    s[i] = (float *)myalloc(xn*sizeof(float));
    if ((i%dt) == 0)
      if ((i > t0)&&(i < t0+stim_tn)){
	/*** Compute the amplitude for this frame ***/
	if (myrand_util_ran2(&seed) > 0.5)
	  curramp -= step;
	else
	  curramp += step;
	while((curramp < -maxamp)||(curramp > maxamp)){
	  if (curramp > maxamp)
	    curramp = maxamp-(curramp-maxamp); /* Reflect the jump. */
	  if (curramp < -maxamp)
	    curramp = -maxamp-(curramp+maxamp); /* Reflect the jump. */
	}
	ampl = (float)curramp/maxamp;
	tt = (float)(i-t0)*tscale*tf;
	for(j=0;j<xn;j++)
	  if ((j >= x0)&&(j < x0+stim_xn))
	    s[i][j] = ampl0 + ampl*sin(((float)(j-x0)*sscale*sf - phase - tt)*
				       2.0*M_PI);
	  else
	    s[i][j] = ampl0;
      }else
	for(j=0;j<xn;j++)
	  s[i][j] = ampl0;
    else
      for(j=0;j<xn;j++)
	s[i][j] = stim_null; /* Fill in between-frame value. */
  }
  transpose_2d_farray(&s,tn,xn);
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                        STIM_MAKE_XT_SUMOSIN_GRATING                       */
/*                                                                           */
/*  Make a sum of sinusoid stimulus (Victor, Shapley, Knight, 1977).         */
/*                                                                           */
/*****************************************************************************/
float **stim_make_xt_sumosin_grating(p)
     struct param_pair_list *p;
{
  int i,j;
  int xn,tn,x0,stim_xn,t0,stim_tn,dt,vskn,pcomp;
  float **s,sscale,tscale,sf,phase,stim_null,ampl,ampl0,sos_ampl,pratio;
  float *tmp;

  vskn = paramfile_get_int_param_or_exit(p,"stim_vsk_n");
  xn = paramfile_get_int_param_or_exit(p,"xn");
  tn = paramfile_get_int_param_or_exit(p,"tn");
  sscale = param_getf_exit(p,"sscale");
  tscale = param_getf_exit(p,"tscale");
  x0 = paramfile_get_int_param_or_exit(p,"stim_x0");
  stim_xn = paramfile_get_int_param_or_exit(p,"stim_xn");
  t0 = paramfile_get_int_param_or_exit(p,"stim_t0");
  stim_tn = paramfile_get_int_param_or_exit(p,"stim_tn");
  dt = paramfile_get_int_param_or_exit(p,"stim_dt");
  sf = param_getf_exit(p,"stim_sf");       /* cyc/deg */
  phase = param_getf_exit(p,"stim_phase"); /* Degrees */
  stim_null = param_getf_exit(p,"stim_null");
  ampl0 = param_getf_exit(p,"stim_ampl0");
  ampl = param_getf_exit(p,"stim_ampl");
  pcomp = paramfile_get_int_param_default(p,"stim_vsk_perturb_comp",-1);
  pratio = paramfile_get_float_param_default(p,"stim_vsk_perturb_ratio",1.0);

  // WYETH HERE - TEMP
  tmp = get_zero_farray(tn);

  phase /= 360.0; /* Change phase to cycles */
  s = (float **)myalloc(tn*sizeof(float *));
  for(i=0;i<tn;i++){
    s[i] = (float *)myalloc(xn*sizeof(float));
    if ((i%dt) == 0)
      if ((i > t0)&&(i < t0+stim_tn)){
	tmp[i]   =        get_vsk_sum_of_sinusoids(vskn,pcomp,pratio,
						   (float)(i-t0)*tscale);
	sos_ampl = ampl * get_vsk_sum_of_sinusoids(vskn,pcomp,pratio,
						   (float)(i-t0)*tscale);
	for(j=0;j<xn;j++)
	  if ((j >= x0)&&(j < x0+stim_xn))
	    s[i][j] = ampl0 + sos_ampl*sin(((float)(j-x0)*sscale*sf - phase)*
					   2.0*M_PI);
	  else
	    s[i][j] = ampl0;
      }else
	for(j=0;j<xn;j++)
	  s[i][j] = ampl0;
    else
      for(j=0;j<xn;j++)
	s[i][j] = stim_null; /* Fill in between-frame value. */
  }
  transpose_2d_farray(&s,tn,xn);

  // WYETH HERE - TEMP
  append_farray_plot("zzt.VSK.pl","VSKt",tmp,tn,1);
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 GET_XT_STIM                               */
/*                                                                           */
/*  Return 's[space][time]' where s[0][0] is thought of as the lower left    */
/*  corner with space running to the right and time running up.              */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - phase should be relative to 'x0' and 't0' of the stimulus patch w/i    */
/*    the background.                                                        */
/*                                                                           */
/*****************************************************************************/
float **get_xt_stim(p)
     struct param_pair_list *p;
{
  char *stimtype;
  float **s;

  stimtype = paramfile_get_char_param_or_exit(p,"stim_type");

  s = NULL;
  if (strcmp(stimtype,"xt_drifting_grating")==0)
    s = stim_make_xt_drifting_grating(p);
  else if (strcmp(stimtype,"xt_jumping_grating")==0)
    s = stim_make_xt_jumping_grating(p);
  else if (strcmp(stimtype,"xt_amc_grating")==0)
    s = stim_make_xt_amc_grating(p);
  else if (strcmp(stimtype,"xt_sumosin_grating")==0)
    s = stim_make_xt_sumosin_grating(p);
  else if (strcmp(stimtype,"xt_impulse")==0)
    s = stim_make_xt_impulse(p);
  else
    exit_error("GET_XT_STIM","Unknown stim_type");
  myfree(stimtype);

  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                             COMPUTE_FPROD_ARRAY                           */
/*                                                                           */
/*  Create an array that will speed the computation of the convolution of    */
/*  the positional noise stimulus with the filter.                           */
/*                                                                           */
/*  Returns:  fprod[zn][fzn]                                                 */
/*                                                                           */
/*  Notes:                                                                   */
/*  - Filter is time reversed.                                               */
/*                                                                           */
/*****************************************************************************/
float **compute_fprod_array(stim,x0,xn,y0,yn,z0,zn,filter,fx0,fxn,fy0,fyn,
			    fz0,fzn)
     float ***stim;
     int x0,xn,y0,yn,z0,zn;
     float ***filter;
     int fx0,fxn,fy0,fyn,fz0,fzn;
{
  int i,j,k,l,jrev;
  float **fprod,sum,norm,**pf,**ps;
  int xoff,yoff,xc,yc;

  /*printf("  COMPUTE_FPROD_ARRAY\n");*/
  /*
    printf("  Frame product array is %d x %d (stim frames by filter frames)\n",
    zn,fzn);
    printf("    Each frame is %d x %d pixels\n",fxn,fyn);*/

  /* Compute offsets (in stimulus) to center filter within stimulus. */
  xoff = x0 + (xn-fxn)/2;
  yoff = y0 + (yn-fyn)/2;
  xc = xoff-fx0;
  yc = yoff-fy0;

  norm = (float)(fxn*fyn);

  /* Compute "frame product" array:  multiply filter frames by stim frames */
  fprod = get_2d_farray(zn,fzn);
  for(i=z0;i<(z0+zn);i++){ /*** For each stimulus frame. ***/
    for(j=fz0;j<(fz0+fzn);j++){
      jrev = 2*fz0+fzn - j - 1; /* Filter is time-reversed. */
      sum = 0.0;
      for(k=fx0;k<(fx0+fxn);k++){
	pf = filter[k];
	ps = stim[k+xc];
	for(l=fy0;l<(fy0+fyn);l++)
	  sum += pf[l][jrev] * ps[l+yc][i];
      }
      fprod[i-z0][j-fz0] = sum/norm;
    }
  }
  /*printf("    Done %d multiplies.\n",zn*fzn*fxn*fyn);*/
  return fprod;
}
/**************************************-**************************************/
/*                                                                           */
/*                          COMPUTE_FPROD_ARRAY_PERIODIC                     */
/*                                                                           */
/*  This version of "compute_fprod_array" takes advantage of period stimuli  */
/*  and computes "fprod" for each unique pattern, using pointers to fill     */
/*  out the rest of the array.                                               */
/*                                                                           */
/*****************************************************************************/
float **compute_fprod_array_periodic(stim,x0,xn,y0,yn,z0,zn,period,
				     filter,fx0,fxn,fy0,fyn,fz0,fzn)
     float ***stim;
     int x0,xn,y0,yn,z0,zn,period;
     float ***filter;
     int fx0,fxn,fy0,fyn,fz0,fzn;
{
  int i,j,k,l,jrev;
  float **fprod,sum,norm;
  int xoff,yoff,xc,yc;

  /*printf("  COMPUTE_FPROD_ARRAY_PERIODIC\n");*/

  /* Compute offsets (in stimulus) to center filter within stimulus. */
  xoff = x0 + (xn-fxn)/2;
  yoff = y0 + (yn-fyn)/2;
  xc = xoff-fx0;
  yc = yoff-fy0;

  norm = (float)(fxn*fyn);

  /* Compute "frame product" array:  multiply filter frames by stim frames */
  /* fprod = get_2d_farray(zn,fzn);*/
  fprod = (float **)myalloc(zn*sizeof(float *));
  for(i=z0;i<(z0+period);i++){ /*** For each stimulus frame. ***/
    fprod[i-z0] = (float *)myalloc(fzn*sizeof(float));
    /*printf("  stim frame %d\n",i);*/
    for(j=fz0;j<(fz0+fzn);j++){
      /*jrev = 2*fz0+fzn - j - 1;*/ /*** FILTER IS TIME-REVERSED  ***/
      jrev = j; /*** FILTER IS NOT TIME-REVERSED. ***/
      sum = 0.0;
      for(k=fx0;k<(fx0+fxn);k++)
	for(l=fy0;l<(fy0+fyn);l++)
	  sum += filter[k][l][jrev]*stim[k+xc][l+yc][i];
      fprod[i-z0][j-fz0] = sum/norm;
    }
  }
  for(i=(z0+period);i<(z0+zn);i++){ /*** For each stimulus frame. ***/
    k = (i-z0)%period + z0;
    fprod[i-z0] = fprod[k];
    /**
      for(j=0;j<fzn;j++){
      k = (i-z0)%period + z0;
      fprod[i-z0][j] = fprod[k][j];
      }**/
  }

  /*printf("    zn= %d  fzn= %d  period= %d\n",zn,fzn,period);*/
  /*
    for(i=0;i<period;i++)
    append_farray_plot("zz.fprod","fprod",fprod[i],fzn,1);*/

  /*printf("    Done.\n");*/
  return fprod;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_RESPONSE_FROM_FPROD                         */
/*                                                                           */
/*  Compute the response to the sequence of frames "fseq" where the product  */
/*  of each stimulus frame with each time slice of the filter is given by    */
/*  "fprod".  "fseq" is a list of integers between 0 and nstim-1.  The       */
/*  response from frame 0 to frame "nfr"-1 is returned.  "mag" indicates     */
/*  the number of filter time units between each stimulus frame.  The        */
/*  returned array has "mag"*"nfr" entries.                                  */
/*                                                                           */
/*  NOTES                                                                    */
/*  - Assuming that filter was reversed in time to make fprod.               */
/*  - fprod[stimframe][filterframe]                                          */
/*  - Assuming values between stimulus frames contribute 0 to response.      */
/*  - If 't0flag'=0, take the first element of the filter to be the origin,  */
/*    otherwise, take the middle of the filter to be zero.                   */
/*  - 'tdelay' is the delay of response relative to stimulus.                */
/*                                                                           */
/*****************************************************************************/
float *get_response_from_fprod(fprod,nstim,nfilt,fseq,nfr,mag,t0flag,tdelay)
     float **fprod;
     int nstim,nfilt,*fseq,nfr,mag,t0flag,tdelay;
{
  int i,j,k,jstim;
  int n,dtmi,a,b,t1,t2,t,dt;
  float *response,r;

  if (t0flag == 0)
    dt = nfilt-1 + tdelay;
  else
    dt = nfilt/2 + tdelay;

  n = mag*nfr;
  response = get_zero_farray(n);
  for(i=0;i<n;i++){
    a = i-dt;  /*** OLD LINE:  a = i-mid;, where mid=nfilt/2  ***/
    b = a + nfilt-1;
    if (a<0)
      t1 = 0;
    else
      t1 = a;
    if (b >= n)
      t2 = n-1;
    else
      t2 = b;
    t = t1+mag-1 - ((t1+mag-1)%mag);
    k = t/mag; /* To speed up inner loop. */
    dtmi = dt-i; /* To speed up inner loop. */
    r = 0.0;
    for(j=t;j<=t2;j+=mag){
      jstim = fseq[k++]; /*** was [j/mag] ***/
      r += fprod[jstim][j + dtmi];
      /*
	if (i==0)
	printf("fprod[%d][%d] = %.4f\n",jstim,j+dtmi,fprod[jstim][j + dtmi]);
      */
    }
    /*printf("  r = %f\n",r);*/
    response[i] = r;
  }

  return response;
}
/**************************************-**************************************/
/*                                                                           */
/*                             TIME_DOMAIN_RESPONSE                          */
/*                                                                           */
/*  This used to be called from 'mod_me_util', before I removed 'fseq'       */
/*  related processing, on 28 Jan 2009.                                     */
/*                                                                           */
/*  Convolve in the time domain.  Each stimulus frame accounts for "k"       */
/*  filter frames.                                                           */
/*                                                                           */
/*  Notes:                                                                   */
/*  - if period > 0, use faster method to compute fprod.                     */
/*                                                                           */
/*****************************************************************************/
void time_domain_response(mtype,stim,x0,xn,y0,yn,z0,zn,period,filter,fx0,fxn,
			  fy0,fyn,fz0,fzn,stim_samp,filt_samp,fseq,nseq,rresp,
			  rn)
     int mtype;
     float ***stim;
     int x0,xn,y0,yn,z0,zn,period;
     float ***filter;
     int fx0,fxn,fy0,fyn,fz0,fzn;
     float stim_samp,filt_samp;  /* sampling rates for stimulus and filter */
     int *fseq,nseq;
     float **rresp;
     int *rn;
{
  int i;
  float *response,**fprod;
  int mag,n,*tseq;

  /*printf("  TIME_DOMAIN_RESPONSE\n");*/

  if ((period > 0) && (mtype == 0))
    fprod = compute_fprod_array_periodic(stim,x0,xn,y0,yn,z0,zn,period,
					 filter,fx0,fxn,fy0,fyn,fz0,fzn);
  else
    fprod = compute_fprod_array(stim,x0,xn,y0,yn,z0,zn,filter,fx0,fxn,fy0,fyn,
				fz0,fzn);

  /* Compute response from fprod array. */
  mag = (int)(0.5 + filt_samp/stim_samp);
  /*printf("    Magnification = %d\n",mag);*/
  n = mag * nseq;

  tseq = (int *)myalloc(nseq*sizeof(int));

  if (mtype == 0){
    for(i=0;i<nseq;i++){
      tseq[i] = fseq[i]%NPAT;
      if (tseq[i] != fseq[i])
	printf("WYETH *\n");
    }
  }else if (mtype == 1){
    for(i=0;i<nseq;i++){
      tseq[i] = fseq[i];
      if (tseq[i] >= NAMP)
	printf("WYETH * ERROR *\n");
    }
  }
  response = get_response_from_fprod(fprod,zn,fzn,tseq,nseq,mag,1,0);

  if (period > 0){
    for(i=0;i<period;i++)
      myfree(fprod[i]);
    myfree(fprod);
  }else
    free_2d_farray(fprod,zn);
  myfree(tseq);

  *rresp = response;
  *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                                SINE_PATTERN                               */
/*                                                                           */
/*  ******* COPIED FROM "pattern.c" (for PEP) and modified. ********         */
/*                                                                           */
/*  Fill the global array "pattern" with a sinusoid.                         */
/*  The values in "pattern" will be centered at 512, ranging from            */
/*  0 to 1024 at most, depending on "max_amp".                               */
/*                                                                           */
/*****************************************************************************/
int *sine_pattern(n,sper,max_amp)
     int n,sper,max_amp;
{
  int i;
  int *p;
  float cf;
  
  if (max_amp > 512)
    max_amp = 512;
  cf = 2.0*M_PI/(float)sper;
  
  p = (int *)myalloc(n*sizeof(int));
  
  for(i=0;i<n;i++){
    p[i] = (int)((float)max_amp*sin(cf*(float)i)) + 512;
    if (p[i] > 1024)
      p[i] = 1024;
    else if (p[i] < 0)
      p[i] = 0;
  }
  return p;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MAKE_FRAME_SET_POS_NOISE                        */
/*                                                                           */
/*  Fill a 3D array with the spatial noise at all possible offsets, where    */
/*  the offset is in units of bar widths.                                    */
/*                                                                           */
/*  "hpos" specifies the horizontal start position of the window in number   */
/*  of bars.  NOTE:  we imagine that the bars are vertical, and the pattern  */
/*  extends horizontally.                                                    */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - 11-06-00 - now the stimulus has zero mean.                             */
/*  - *** WYETH - this is a big hack, should update to arb. sinusoid.        */
/*                                                                           */
/*****************************************************************************/
void make_frame_set_pos_noise(min_wavlen,wh,rstim,x0,xn,y0,yn,z0,zn,hpos)
     int min_wavlen,wh; // wh - is width and height
     float ****rstim;
     int *x0,*xn,*y0,*yn,*z0,*zn,hpos;
{
  int i,j,k,l;
  int *pattern,win_w,win_h,bar_width,patamp,nframes,pati;
  float ***stim; // [x][y][t]

  bar_width = 1;
  patamp = 512;
  nframes = NPAT; /* WYETH - Probably not nec. to always make 256 posistions */
  win_w = win_h = wh;

  if ((hpos < 0)||((hpos+win_w/bar_width) > NPAT)){
    printf("hpos = %d\n",hpos);
    exit_error("MAKE_FRAME_SET_POS_NOISE","Invalid hpos");
  }

  pattern = sine_pattern(NPAT,min_wavlen,patamp);
  /*write_iarray("zz.pattern",pattern,NPAT);*/

  stim = f3tensor(1,win_w,1,win_h,1,nframes);
  for(i=0;i<nframes;i++){
    for(j=0;j<(win_w/bar_width);j++){ /* j is bar number */
      pati = (i+hpos+j)%NPAT;
      for(l=0;l<bar_width;l++) /*** Fill in one bar. ***/
	for(k=1;k<=win_h;k++){
	  stim[bar_width*j+l+1][k][i+1] = (float)pattern[pati]/4 - 128.0;
	}
    }
  }

  *rstim = stim;
  *x0 = *y0 = *z0 = 1;
  *xn = win_w; *yn = win_h;
  *zn = nframes;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MAKE_FRAME_SET_WYPHASM_NOISE                     */
/*                                                                           */
/*  Fill a 3D array with two frames, preferred and counterphase, for each    */
/*  sper.  For each 'sper' value, sper, sper*2, sper*4, ..., make a pair of  */
/*  frames.                                                                  */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Should keep phase constant in center of stimulus.                      */
/*  - Looks like order is:  black center then white center.                  */
/*                                                                           */
/*****************************************************************************/
void make_frame_set_wyphasm_noise(sper,nsper,rstim,x0,xn,y0,yn,z0,zn,hpos)
     int sper,nsper;
     float ****rstim;
     int *x0,*xn,*y0,*yn,*z0,*zn,hpos;
{
  int i,j,k;
  int *pattern,win_w,win_h,patamp,nframes,pati,tsper,tpos;
  float ***stim; // [x][y][t]

  patamp = 512;
  nframes = 2*nsper;
  win_w = win_h = 64;

  stim = f3tensor(1,win_w,1,win_h,1,nframes);

  if ((hpos < 0)||((hpos+win_w) > NPAT))
    exit_error("MAKE_FRAME_SET_WYPHASM_NOISE","Invalid hpos");

  tsper = sper;
  for(i=0;i<nsper;i++){
    pattern = sine_pattern(NPAT,tsper,patamp);
    /*append_iarray_plot("zzz.sinetab","tt",pattern,NPAT,1);*/
    tpos = NPAT + tsper/4 - win_w/2; /* To keep phase constant in center */
    for(j=0;j<win_w;j++){ /* j is bar number */
      pati = (tpos+hpos+j + tsper/2)%NPAT;
      for(k=1;k<=win_h;k++)
	stim[j+1][k][i*2+1] = (float)pattern[pati]/4 - 128.0;
    }
    for(j=0;j<win_w;j++){ /* j is bar number */
      pati = (tpos+hpos+j)%NPAT;
      for(k=1;k<=win_h;k++)
	stim[j+1][k][i*2+2] = (float)pattern[pati]/4 - 128.0;
    }

    myfree(pattern);
    tsper *= 2;
  }
  
  *rstim = stim;
  *x0 = *y0 = *z0 = 1;
  *xn = win_w; *yn = win_h;
  *zn = nframes;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MAKE_FRAME_SET_PHASE_CON_NOISE                    */
/*                                                                           */
/*  Fill a 3D array with two frames, preferred and counterphase, at several  */
/*  contrast values (ncon).                                                  */
/*                                                                           */
/*****************************************************************************/
void make_frame_set_phase_con_noise(sper,ncon,rstim,x0,xn,y0,yn,z0,zn,hpos)
     int sper,ncon;
     float ****rstim;
     int *x0,*xn,*y0,*yn,*z0,*zn,hpos;
{
  int i,j,k;
  int *pattern,win_w,win_h,patamp,nframes,pati,tpos;
  float ***stim; /* [x][y][t] */

  patamp = 512;
  nframes = 2*ncon;
  win_w = win_h = 64;

  stim = f3tensor(1,win_w,1,win_h,1,nframes);

  if ((hpos < 0)||((hpos+win_w) > NPAT))
    exit_error("MAKE_FRAME_SET_PHASE_CON_NOISE","Invalid hpos");

  pattern = sine_pattern(NPAT,sper,patamp);

  /**** WYETH FIX ****/

  for(i=0;i<ncon;i++){
    /*append_iarray_plot("zzz.sinetab","tt",pattern,NPAT,1);*/
    tpos = NPAT + sper/4 - win_w/2; /* To keep phase constant in center */
    for(j=0;j<win_w;j++){ /* j is bar number */
      pati = (tpos+hpos+j + sper/2)%NPAT;
      for(k=1;k<=win_h;k++)
	stim[j+1][k][i*2+1] = (float)pattern[pati]/4 - 128.0;
    }
    for(j=0;j<win_w;j++){ /* j is bar number */
      pati = (tpos+hpos+j)%NPAT;
      for(k=1;k<=win_h;k++)
	stim[j+1][k][i*2+2] = (float)pattern[pati]/4 - 128.0;
    }
  }
  myfree(pattern);
  
  *rstim = stim;
  *x0 = *y0 = *z0 = 1;
  *xn = win_w; *yn = win_h;
  *zn = nframes;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MAKE_FRAME_SET_AMP_NOISE                         */
/*                                                                           */
/*  Fill a 3D array with the spatial noise at all possible amplitudes.       */
/*                                                                           */
/*  "hpos" specifies the horizontal start position of the window in number   */
/*  of bars.  NOTE:  we imagine that the bars are vertical, and the pattern  */
/*  extends horizontally.                                                    */
/*                                                                           */
/*****************************************************************************/
void make_frame_set_amp_noise(min_wavlen,rstim,x0,xn,y0,yn,z0,zn,hpos)
     int min_wavlen;
     float ****rstim;
     int *x0,*xn,*y0,*yn,*z0,*zn,hpos;
{
  int i,j,k,l;
  int *pattern,win_w,win_h,bar_width,patamp,nframes,pati;
  float ***stim,tamp,tp; /* [x][y][t] */

  bar_width = 1;
  /*patseed = 1777;*/
  patamp = 512;

  nframes = NAMP;  /* e.g. 257 */

  win_w = win_h = 64;

  if ((hpos < 0)||((hpos+win_w/bar_width) > NPAT))
    exit_error("MAKE_FRAME_SET_AMP_NOISE","Invalid hpos");

  pattern = sine_pattern(NPAT,min_wavlen,patamp);
  /*write_iarray("zz.pattern",pattern,NPAT);*/

  stim = f3tensor(1,win_w,1,win_h,1,nframes);
  for(i=0;i<nframes;i++){
    tamp = (-(float)((NAMP-1)/2) + (float)i) / (float)((NAMP-1)/2);
    /*printf("i=%d  tamp = %.2f\n",i,tamp);*/
    for(j=0;j<(win_w/bar_width);j++){ /* j is bar number */
      pati = (hpos+j)%NPAT;  /*** DONT LET THE PATTERN MOVE ***/
      for(l=0;l<bar_width;l++) /*** Fill in one bar. ***/
	for(k=1;k<=win_h;k++){
	  tp = (float)patamp + tamp * (float)(pattern[pati] - patamp);
	  stim[bar_width*j+l+1][k][i+1] = tp;  /* [i] changed to [i+1] */
	}
    }
  }

  *rstim = stim;
  *x0 = *y0 = *z0 = 1;
  *xn = win_w; *yn = win_h;
  *zn = nframes;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_FRAME_SEQUENCE_DIRMOD                       */
/*                                                                           */
/*  Get the sequence of frame offsets for a positional noise stimulus.       */
/*                                                                           */
/*****************************************************************************/
void get_frame_sequence_dirmod(framerate,jumpseed,dtframe,period,distrib,
			       minjump,maxjump,rpt,mu,sigma,rfseq,rn,rvel)
     int framerate,jumpseed,dtframe;
     int period,distrib,minjump,maxjump,rpt,mu,sigma,**rfseq,*rn,**rvel;
{
  int i;
  int n,offset,*fseq,*vel;
  float *xdata,*ydata;

  get_wn_stim_for_params(framerate,jumpseed,dtframe,period,distrib,mu,sigma,
			 minjump,maxjump,rpt,1.0,&xdata,&ydata,&n);
  offset = 0;
  fseq = (int *)myalloc(n*sizeof(int));
  vel = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++){
    vel[i] = my_rint(ydata[i]);

    /*offset -= my_rint(ydata[i]);*/
    offset += my_rint(ydata[i]); /*** WYETH - try plus here ***/

    if (offset >= NPAT) /*** keep offset in range 0..NPAT-1 ***/
      offset -= NPAT;
    else if (offset < 0)
      offset += NPAT;
    fseq[i] = offset;
  }

  myfree(xdata); myfree(ydata);
  *rfseq = fseq; *rn = n; *rvel = vel;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_FRAME_SEQUENCE_AMP                         */
/*                                                                           */
/*  Get the sequence of frame indices for contrast modulation.               */
/*                                                                           */
/*****************************************************************************/
void get_frame_sequence_amp(framerate,jumpseed,dtframe,period,distrib,minjump,
			    maxjump,mu,sigma,walkflag,steptype,amp0,dt,
			    rfseq,rn,rvel)
     int framerate,jumpseed,dtframe,period,distrib,minjump,maxjump,mu,sigma;
     int walkflag,steptype,amp0;
     float dt;
     int **rfseq,*rn,**rvel;
{
  int i;
  int n,offset,*fseq,*vel,maxampv;
  float *xdata,*ydata;

  get_amp_mod_stim_for_params(framerate,jumpseed,dtframe,period,distrib,mu,
			      sigma,minjump,maxjump,walkflag,steptype,amp0,
			      dt,&xdata,&ydata,&n);
  maxampv = (NAMP-1);
  fseq = (int *)myalloc(n*sizeof(int));
  vel = (int *)myalloc(n*sizeof(int));
  offset = (NAMP-1)/2; /* Middle of table is 0 contrast. */
  for(i=0;i<n;i++){
    vel[i] = my_rint(ydata[i]); /* This isn't velocity */
    offset += my_rint(ydata[i]); /* Changed to += (from -=) */
    if ((offset < 0) || (offset > maxampv)){
      printf("offset = %d\n",offset);
      exit_error("GET_FRAME_SEQUENCE_AMP","Offset out of bounds");
    }
    fseq[i] = offset;
  }
  myfree(xdata); myfree(ydata);

  *rfseq = fseq; *rn = n; *rvel = vel;
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_FRAME_SEQUENCE_WYPHASM                       */
/*                                                                           */
/*  Get the sequence of frame offsets for a positional noise stimulus.       */
/*                                                                           */
/*  NOTES                                                                    */
/*  fseq - holds 0,1                                                         */
/*  vel  - holds -1,+1                                                       */
/*                                                                           */
/*****************************************************************************/
void get_frame_sequence_wyphasm(framerate,jumpseed,dtframe,period,distrib,
				minjump,maxjump,rpt,mu,sigma,rfseq,rn,rvel)
     int framerate,jumpseed,dtframe;
     int period,distrib,minjump,maxjump,rpt,mu,sigma,**rfseq,*rn,**rvel;
{
  int i;
  int n,offset,*fseq,*vel;
  float *xdata,*ydata;

  /*printf("%d %d %d %d %d\n",framerate,jumpseed,dtframe,period,distrib);
    printf("%d %d %d %d\n",mu,sigma,minjump,maxjump);*/
  get_wn_stim_for_params(framerate,jumpseed,dtframe,period,distrib,mu,sigma,
			 minjump,maxjump,rpt,1.0,&xdata,&ydata,&n);

  offset = 0; /* Wyeth, this line added on April 24th, 2000 */

  fseq = (int *)myalloc(n*sizeof(int));
  vel = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++){
    vel[i] = my_rint(ydata[i]);
    if (vel[i] == 0)
      exit_error("GET_FRAME_SEQUENCE_WYPHASM","vel is 0\n");

    offset += my_rint(ydata[i]); /*** WYETH - try plus here ***/

    if (vel[i] > 0)
      fseq[i] = 1;
    else
      fseq[i] = 0;
  }

  myfree(xdata); myfree(ydata);
  *rfseq = fseq; *rn = n; *rvel = vel;
}
