/*****************************************************************************/
/*                                                                           */
/*  wm_util.c                                                                */
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
#include "myrand_util.h"
#include "spike_util.h"
#include "ndata_util.h"
#include "paramfile_util.h"
#include "fft_util.h"
#include "stim_util.h"
#include "stm_util.h"
#include "mod_util.h"
#include "pop_util.h"
#include "ifc_util.h"
#include "mod_dog_util.h"
#include "mod_lin_util.h"
#include "mod_me_util.h"
#include "mod_wav_util.h"
#include "mod_mesh_util.h"
#include "mod_pop_util.h"
#include "mod_vhf_util.h"
#include "mod_dcn_util.h"
#include "mod_test_util.h"
#include "mod_x_util.h"
#include "mod_srf_util.h"
#include "mod.h"
#include "ndata.h"
#include "paramfile.h"

/**************************************-**************************************/
/*                                                                           */
/*                               WM_MOO_VAR_PREP                             */
/*                                                                           */
/*****************************************************************************/
void wm_moo_var_prep(m,s,mylogf)
     struct model_struct *m; // Model params
     struct stim_struct *s;  // Stimulus params
     char *mylogf;
{
  struct onode *to,*vpo;
  char *infile;

  mylog(mylogf,"    WM_MOO_VAR_PREP\n");

  vpo = onode_child_get_unique(m->o,"var_param");
  if (vpo == NULL){
    if (onode_item(m->o,"var_param_file")==1){
      mylog(mylogf,"      Reading moo var params from file\n");
      //
      //  The var params are stored in a separate file
      //
      infile = onode_getpar_chr_exit(m->o,"var_param_file");
      to = paramfile_onode_file_read(infile,1);   // Read from the file
      myfree(infile);

      vpo = onode_child_get_unique(to,"var_param");
      if (vpo == NULL)
	exit_error("WM_MOO_VAR_PREP","<var_param> not found in file.");

      // Add this as a child of m->o

      onode_insert_child_at_end(m->o,vpo);
    }else{
      return; // There are no moo var params, we're done here
    }
  }else{
    mylog(mylogf,"      Found 'var_param' child\n");
  }

  //
  //  Extract the model var params from onodes, build var par records
  //
  m->vps = mod_util_var_par(mylogf,m,vpo);

  m->nrun = m->vps->m_n;  // The number of different model configs to run
}
/**************************************-**************************************/
/*                                                                           */
/*                            WM_PREP_VAR_CONST_GEN                          */
/*                                                                           */
/*****************************************************************************/
void wm_prep_var_const_gen(m,s,mylogf,arglist,narg,parflag)
     struct model_struct *m; // Model params
     struct stim_struct *s;  // Stimulus params
     char *mylogf;
     char **arglist;  // [narg]
     int narg;
     int *parflag;    // [narg] 1's indicate which params used by .moo file
{
  int i,k;
  int n,seed,nmoo,flag_moo_tn,flag_moo_tscale,flag_moo_sscale;
  char ggstr[LONG_SLEN],tval[SLEN];
  char **moo_name,**moo_val;

  // Model variable params
  if (m != NULL)  // If there is a .moo file
    if (m->ppl == NULL)  // Added 2013 Oct, so that .mod files work
      wm_moo_var_prep(m,s,mylogf);

  // Get variable param names, values, number, and number of trials w/ rpt
  prep_var_params(mylogf,s->ppl,1,&(s->vval),&(s->vname),&(s->vtype),
		  &(s->nvar),&(s->ntr),&(s->val),&(s->vcnt),&(s->nvl),
		  &(s->nvp));

  //
  //  (1) Any .moo params that were changed on the command line will be stored
  //      as 'const' params.
  //  (2) 'tn' and 'tscale' will ALWAYS be stored.
  //
  // Get model params that were changed on the param line

  if (m == NULL)
    nmoo = 0;
  else{

    flag_moo_tn = 0;      // Assume 'tn' has not yet been included
    flag_moo_tscale = 0;  // Assume 'tscale' has not yet been included
    flag_moo_sscale = 0;  // Assume 'sscale' has not yet been included
    nmoo = 3;             // There are at least 3 .moo params to be stored

    if ((narg > 0) && (parflag != NULL))
      nmoo += sum_iarray(parflag,narg,0,narg);

    moo_name = (char **)myalloc(nmoo*sizeof(char *));
    moo_val  = (char **)myalloc(nmoo*sizeof(char *));

    k = 0;
    for(i=0;i<narg;i++){  // Include any .moo params named on the command line
      if (parflag[i] == 1){
	moo_name[k] = strdup(arglist[i]);
	moo_val[k]  = strdup(arglist[i+1]);
	
	if (strcmp(moo_name[k],"tn")==0)
	  flag_moo_tn = 1;
	else if (strcmp(moo_name[k],"tscale")==0)
	  flag_moo_tscale = 1;
	else if (strcmp(moo_name[k],"sscale")==0)
	  flag_moo_sscale = 1;

	//printf("moo_name  _val  %s  %s\n",moo_name[k],moo_val[k]);
	k += 1;
      }
    }

    // Now add 'tn' & 'tscale' if not done already
    if (flag_moo_tn == 0){
      moo_name[k] = strdup("tn");
      if (m->ppl != NULL)
	moo_val[k] = param_getc_exit(m->ppl,"tn");
      else{
	moo_val[k] = onode_getpar_chr_exit(m->o,"tn");
      }
      k += 1;
    }
    if (flag_moo_tscale == 0){
      moo_name[k] = strdup("tscale");
      if (m->ppl != NULL)
	moo_val[k] = param_getc_exit(m->ppl,"tscale");
      else{
	moo_val[k] = onode_getpar_chr_exit(m->o,"tscale");
      }
      k += 1;
    }
    if (flag_moo_sscale == 0){
      moo_name[k] = strdup("sscale");
      if (m->ppl != NULL)
	moo_val[k] = param_getc_exit(m->ppl,"sscale");
      else{
	moo_val[k] = onode_getpar_chr_exit(m->o,"sscale");
      }
      k += 1;
    }
    nmoo = k;  // There might be fewer than original 'nmoo' value.
  }

  prep_const_params(mylogf,s,moo_name,moo_val,nmoo,&(s->cval),&(s->cname),
		    &(s->ctype),&(s->ncon));

  mylog(mylogf,"    Const params\n");
  for(i=0;i<s->ncon;i++){
    if (strlen(s->cval[i]) > 64){
      //
      //  This added to protect against very long const param values, which
      //  occurs for long lists of values, such as filenames for .fst stimuli.
      //
      memcpy(tval,s->cval[i],60);
      tval[60] = ' ';
      tval[61] = '.';
      tval[62] = '.';
      tval[63] = '.';
      tval[64] = '\0';
      sprintf(ggstr,"    %2d %12s %c %s\n",i,s->cname[i],s->ctype[i],tval);
    }else{
      sprintf(ggstr,"    %2d %12s %c %s\n",i,s->cname[i],s->ctype[i],
	      s->cval[i]);
    }
    mylog(mylogf,ggstr);
  }
  if (s->nvar > 0){
    mylog(mylogf,"    Var params\n");
    for(i=0;i<s->nvar;i++){
      sprintf(ggstr,"    %2d %12s %c\n",i,s->vname[i],s->vtype[i]);
      mylog(mylogf,ggstr);
    }
  }

  // Set up seeds for each trial for model noise
  n = s->ntr;

  if (m == NULL){
    mylog(mylogf,"  WM_PREP_VAR_CONST_GEN - No model, model seeds not set.\n");
  }else{
    if (m->ppl != NULL)
      seed = paramfile_get_int_param_default(m->ppl,"model_noise_seed",1777);
    else{
      seed = onode_getpar_int_dflt(m->o,"model_noise_seed",1777);
    }

    m->nseed = n;
    m->mseed = get_seeds(seed,100000,n);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              WM_RESPONSE_COUNT                            */
/*                                                                           */
/*  Compute total number of 'ndata' responses.  Modify 'r->n' only.          */
/*                                                                           */
/*****************************************************************************/
void wm_response_count(m,r,mylogf)
     struct model_struct *m;
     struct response_struct *r;
     char *mylogf;
{
  int nx,ny,nz;
  char *pop_name;
  struct param_pair *pp;

  pp = paramfile_get_first_prefix_pointer(r->ppl,"resp_save_");
  if (pp != NULL)
    mylog_exit(mylogf,"WM_RESPONSE_COUNT CHANGE resp_save_ to save_as_  !!\n");


  r->n = 0;
  r->n += paramfile_count_prefix_param(r->ppl,"save_as");
  r->n += paramfile_count_prefix_param(r->ppl,"save_pop_unit_as");
  r->n += paramfile_count_prefix_param(r->ppl,"save_pop_syn_as");

  pp = paramfile_get_first_prefix_pointer(r->ppl,"save_pop_grid_as_");
  while(pp != NULL){
    nx = paramfile_get_nth_int_param_ptr_or_exit(pp,6);
    ny = paramfile_get_nth_int_param_ptr_or_exit(pp,7);
    nz = paramfile_get_nth_int_param_ptr_or_exit(pp,8);
    r->n += (nx * ny * nz);

    if ((nx == 0) || (ny == 0) || (nz == 0)){
      printf("  *** In 'save_pop_grid_as_... %d %d %d ...'",nx,ny,nz);
      printf("   none of these lengths should be < 1\n");
      exit_error("WM_RESPONSE_COUNT","Instruction to record 0 data");
    }

    pp = paramfile_get_next_prefix_pointer(pp,"save_pop_grid_as_");
  }

  pp = paramfile_get_first_prefix_pointer(r->ppl,"save_pop_layer_as_");
  while(pp != NULL){
    pop_name = paramfile_get_nth_char_param_pp_or_exit(pp,2);
    pop_util_onode_get_pop_xn_yn_zn(m->o,pop_name,&nx,&ny,&nz);
    if (nx == -1){ // All -1's are returned if pop not found
      mylog(mylogf,"    Using model xn,yn for layer dimensions");
      nx = onode_getpar_int_exit(m->o,"xn");
      ny = onode_getpar_int_exit(m->o,"yn");
      nz = 1;
    }

    r->n += (nx * ny * nz);

    pp = paramfile_get_next_prefix_pointer(pp,"save_pop_layer_as_");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              WM_RESPONSE_CONFIG                           */
/*                                                                           */
/*  Fill in the responses description information in 'r' (see mod.h).        */
/*                                                                           */
/*****************************************************************************/
void wm_response_config(m,r,mylogf)
     struct model_struct *m;
     struct response_struct *r;
     char *mylogf;
{
  int i,j,k,l;
  int nval,gx0,gy0,gz0,gxn,gyn,gzn;
  float gsamp;
  char *savtype,*gname,*gnd_name,tstr[SLEN],*gdid,*tdup;
  struct param_pair *pp;

  r->rformat =   (int *)myalloc(r->n*sizeof(int));
  r->datid   = (char **)myalloc(r->n*sizeof(char *));

  //r->rflag   = get_zero_iarray(r->n);  // Nothing stored yet, all flags to 0

  r->plname  = (char **)myalloc(r->n*sizeof(char *));
  r->xi      =   (int *)myalloc(r->n*sizeof(int));
  r->yi      =   (int *)myalloc(r->n*sizeof(int));
  r->zi      =   (int *)myalloc(r->n*sizeof(int));

  r->plname1 = (char **)myalloc(r->n*sizeof(char *));
  r->xi1     =   (int *)myalloc(r->n*sizeof(int));
  r->yi1     =   (int *)myalloc(r->n*sizeof(int));
  r->zi1     =   (int *)myalloc(r->n*sizeof(int));

  r->nd_name = (char **)myalloc(r->n*sizeof(char *));
  r->rtype   = (char **)myalloc(r->n*sizeof(char *));
  r->samp    = (float *)myalloc(r->n*sizeof(float));
  r->ri      =   (int *)myalloc(r->n*sizeof(int));

  k = 0;
  r->ns = r->nf = 0;

  //  SAVE_AS_
  pp = paramfile_get_first_prefix_pointer(r->ppl,"save_as_");
  while(pp != NULL){
    nval = paramfile_count_values_pointer(pp);
    if (nval != 3)
      exit_error("WM_RESPONSE_CONFIG","save_as_ bad value list");

    r->rformat[k] = 0; // save_as_
    r->nd_name[k] = strdup(pp->name+8);

    savtype = paramfile_get_nth_char_param_pp_or_exit(pp,0);
    r->samp[k] = paramfile_get_nth_float_param_ptr_or_exit(pp,1);
    r->plname[k] = strdup("NULL");
    r->datid[k] = paramfile_get_nth_char_param_pp_or_exit(pp,2);
    r->xi[k] = r->yi[k] = r->zi[k] = 0; // Not used for rformat 0

    if (strcmp(savtype,"s")==0){
      r->rtype[k] = strdup("s");
      r->ri[k] = r->ns;
      r->ns += 1;
    }else if (strcmp(savtype,"f")==0){
      r->rtype[k] = strdup("f");
      r->ri[k] = r->nf;
      r->nf += 1;
    }else{
      printf("*** Error in .rsp file:\n");
      printf("  *** %s  %s\n",pp->name,pp->value);
      printf("  *** 'savtype' = %s\n",savtype);
      exit_error("WM_RESPONSE_CONFIG","save_as_ bad type");
    }

    myfree(savtype);
    k += 1;
    pp = paramfile_get_next_prefix_pointer(pp,"save_as_");
  }

  //  SAVE_POP_UNIT_AS_
  pp = paramfile_get_first_prefix_pointer(r->ppl,"save_pop_unit_as_");
  while(pp != NULL){
    nval = paramfile_count_values_pointer(pp);
    if (nval != 7)
      exit_error("WM_RESPONSE_CONFIG","save_pop_unit_as bad value list");

    r->rformat[k] = 1; // save_pop_unit_as_
    r->nd_name[k] = strdup(pp->name+17);

    savtype      = paramfile_get_nth_char_param_pp_or_exit(pp,0);
    r->samp[k]   = paramfile_get_nth_float_param_ptr_or_exit(pp,1);
    r->plname[k] = paramfile_get_nth_char_param_pp_or_exit(pp,2);
    r->xi[k]     = paramfile_get_nth_int_param_ptr_or_exit(pp,3);
    r->yi[k]     = paramfile_get_nth_int_param_ptr_or_exit(pp,4);
    r->zi[k]     = paramfile_get_nth_int_param_ptr_or_exit(pp,5);
    r->datid[k]  = paramfile_get_nth_char_param_pp_or_exit(pp,6);

    if (strcmp(savtype,"s")==0){
      r->rtype[k] = strdup("s");
      r->ri[k] = r->ns;
      r->ns += 1;
    }else if (strcmp(savtype,"f")==0){
      r->rtype[k] = strdup("f");
      r->ri[k] = r->nf;
      r->nf += 1;
    }else{
      printf("*** Error in .rsp file:\n");
      printf("  *** %s  %s\n",pp->name,pp->value);
      printf("  *** 'savtype' = %s\n",savtype);
      exit_error("WM_RESPONSE_CONFIG","save_pop_unit_as_ bad type");
    }

    myfree(savtype);
    k += 1;
    pp = paramfile_get_next_prefix_pointer(pp,"save_pop_unit_as_");
  }

  //  SAVE_POP_SYN_AS_
  pp = paramfile_get_first_prefix_pointer(r->ppl,"save_pop_syn_as_");
  while(pp != NULL){
    nval = paramfile_count_values_pointer(pp);
    if (nval != 11)
      exit_error("WM_RESPONSE_CONFIG","save_pop_syn_as bad value list");

    r->rformat[k] = 10; // save_pop_syn_as_
    r->nd_name[k] = strdup(pp->name+16);

    savtype      = paramfile_get_nth_char_param_pp_or_exit(pp,0);
    r->samp[k]   = paramfile_get_nth_float_param_ptr_or_exit(pp,1);

    r->plname[k] = paramfile_get_nth_char_param_pp_or_exit(pp,2);
    r->xi[k]     = paramfile_get_nth_int_param_ptr_or_exit(pp,3);
    r->yi[k]     = paramfile_get_nth_int_param_ptr_or_exit(pp,4);
    r->zi[k]     = paramfile_get_nth_int_param_ptr_or_exit(pp,5);

    r->plname1[k] = paramfile_get_nth_char_param_pp_or_exit(pp,6);
    r->xi1[k]     = paramfile_get_nth_int_param_ptr_or_exit(pp,7);
    r->yi1[k]     = paramfile_get_nth_int_param_ptr_or_exit(pp,8);
    r->zi1[k]     = paramfile_get_nth_int_param_ptr_or_exit(pp,9);

    r->datid[k]  = paramfile_get_nth_char_param_pp_or_exit(pp,10);

    if (strcmp(savtype,"s")==0){
      r->rtype[k] = strdup("s");
      r->ri[k] = r->ns;
      r->ns += 1;
    }else if (strcmp(savtype,"f")==0){
      r->rtype[k] = strdup("f");
      r->ri[k] = r->nf;
      r->nf += 1;
    }else{
      printf("*** Error in .rsp file:\n");
      printf("  *** %s  %s\n",pp->name,pp->value);
      printf("  *** 'savtype' = %s\n",savtype);
      exit_error("WM_RESPONSE_CONFIG","save_pop_syn_as_ bad type");
    }

    myfree(savtype);
    k += 1;
    pp = paramfile_get_next_prefix_pointer(pp,"save_pop_syn_as_");
  }

  //  SAVE_POP_GRID_AS_
  pp = paramfile_get_first_prefix_pointer(r->ppl,"save_pop_grid_as_");
  while(pp != NULL){
    nval = paramfile_count_values_pointer(pp);
    if (nval != 10)
      exit_error("WM_RESPONSE_CONFIG","save_pop_grid_as_ bad value list");
    
    gnd_name = strdup(pp->name+17);
    
    savtype = paramfile_get_nth_char_param_pp_or_exit(pp,0);
    gsamp = paramfile_get_nth_float_param_ptr_or_exit(pp,1);
    gname = paramfile_get_nth_char_param_pp_or_exit(pp,2);
    gx0   = paramfile_get_nth_int_param_ptr_or_exit(pp,3);
    gy0   = paramfile_get_nth_int_param_ptr_or_exit(pp,4);
    gz0   = paramfile_get_nth_int_param_ptr_or_exit(pp,5);
    gxn   = paramfile_get_nth_int_param_ptr_or_exit(pp,6);
    gyn   = paramfile_get_nth_int_param_ptr_or_exit(pp,7);
    gzn   = paramfile_get_nth_int_param_ptr_or_exit(pp,8);
    gdid  = paramfile_get_nth_char_param_pp_or_exit(pp,9);
    
    for(i=0;i<gxn;i++){
      for(j=0;j<gyn;j++){
	for(l=0;l<gzn;l++){
	  r->rformat[k] = 2; // save_pop_grid_as_

	  r->samp[k] = gsamp;
	  if (gzn == 1)
	    sprintf(tstr,"%s_%d_%d",gnd_name,gx0+i,gy0+j);
	  else
	    sprintf(tstr,"%s_%d_%d_%d",gnd_name,gx0+i,gy0+j,gz0+l);
	  
	  r->nd_name[k] = strdup(tstr);
	  r->plname[k] = strdup(gname);
	  r->xi[k]   = gx0+i;
	  r->yi[k]   = gy0+j;
	  r->zi[k]   = gz0+l;
	  r->datid[k]  = strdup(gdid);
	  
	  if (strcmp(savtype,"s")==0){
	    r->rtype[k] = strdup("s");
	    r->ri[k] = r->ns;
	    r->ns += 1;
	  }else if (strcmp(savtype,"f")==0){
	    r->rtype[k] = strdup("f");
	    r->ri[k] = r->nf;
	    r->nf += 1;
	  }else{
	    printf("*** Error in .rsp file:\n");
	    printf("  *** %s\n",tstr);
	    printf("  *** 'savtype' = %s\n",savtype);
	    exit_error("WM_RESPONSE_CONFIG","save_pop_grid_as_ bad type");
	  }
	  
	  k += 1;
	}
      }
    }
    myfree(savtype); myfree(gname); myfree(gnd_name); myfree(gdid);
    pp = paramfile_get_next_prefix_pointer(pp,"save_pop_grid_as_");
  }

    
  //  SAVE_POP_LAYER_AS_
  pp = paramfile_get_first_prefix_pointer(r->ppl,"save_pop_layer_as_");
  while(pp != NULL){
    nval = paramfile_count_values_pointer(pp);
    if (nval != 4)
      exit_error("WM_RESPONSE_CONFIG","save_pop_layer_as_ bad value list");

    gnd_name = strdup(pp->name+18);

    savtype = paramfile_get_nth_char_param_pp_or_exit(pp,0);
    gsamp = paramfile_get_nth_float_param_ptr_or_exit(pp,1);
    gname = paramfile_get_nth_char_param_pp_or_exit(pp,2);
    gdid  = paramfile_get_nth_char_param_pp_or_exit(pp,3);

    pop_util_onode_get_pop_xn_yn_zn(m->o,gname,&gxn,&gyn,&gzn);
    if (gxn == -1){ // All -1's are returned if pop not found
      gxn = onode_getpar_int_exit(m->o,"xn");
      gyn = onode_getpar_int_exit(m->o,"yn");
      gzn = 1;
    }

    for(i=0;i<gxn;i++){
      for(j=0;j<gyn;j++){
	for(l=0;l<gzn;l++){
	  r->rformat[k] = 3; // save_pop_layer_as_

	  r->samp[k] = gsamp;
	  if (gzn == 1)
	    sprintf(tstr,"%s_%d_%d",gnd_name,i,j);
	  else
	    sprintf(tstr,"%s_%d_%d_%d",gnd_name,i,j,l);
	  
	  r->nd_name[k] = strdup(tstr);
	  r->plname[k] = strdup(gname);
	  r->xi[k]   = i;
	  r->yi[k]   = j;
	  r->zi[k]   = l;
	  r->datid[k]  = strdup(gdid);

	  if (strcmp(savtype,"s")==0){
	    r->rtype[k] = strdup("s");
	    r->ri[k] = r->ns;
	    r->ns += 1;
	  }else if (strcmp(savtype,"f")==0){
	    r->rtype[k] = strdup("f");
	    r->ri[k] = r->nf;
	    r->nf += 1;
	  }else{
	    printf("*** Error in .rsp file:\n");
	    printf("  *** %s\n",tstr);
	    printf("  *** 'savtype' = %s\n",savtype);
	    exit_error("WM_RESPONSE_CONFIG","save_pop_layer_as_ bad type");
	  }

	  k += 1;
	}
      }
    }
    myfree(savtype); myfree(gname); myfree(gnd_name); myfree(gdid);
    pp = paramfile_get_next_prefix_pointer(pp,"save_pop_layer_as_");
  }

  //
  //  Check that there are no redundant 'nd_names'
  //
  tdup = carray_get_duplicate_ptr(r->nd_name,r->n);
  if (tdup != NULL){
    sprintf(tstr,"Duplicate name for saved data: %s",tdup);
    mylogx(mylogf,"WM_RESPONSE_CONFIG",tstr);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                        WM_RESPONSE_CHECK_TRIAL_FULL                       */
/*                                                                           */
/*  Verify that all requested responses have been filled for trial 'ti'.     */
/*                                                                           */
/*****************************************************************************/
void wm_response_check_trial_full(r,ti,mylogf,proc_name)
     struct response_struct *r;
     int ti;                     // trial index
     char *mylogf;
     char *proc_name;            // computer name for mm, or NULL
{
  int i;
  int di,n;
  char tstr[SLEN];

  for(i=0;i<r->n;i++){  // For each response requested
    di = r->ri[i];  // Index into either 'cnt' or 'fcnt'

    // Get the number 'n' of data values stored for this response.
    if (strcmp(r->rtype[i],"s")==0){
      n = r->cnt[di][ti];
    }else if (strcmp(r->rtype[i],"f")==0){
      n = r->fcnt[di][ti];
    }else{
      exit_error("WM_RESPONSE_CHECK_TRIAL_FULL","Unknown response type");
    }

    // -1 is the signal that the response was never filled.
    if (n < 0){
      if (proc_name != NULL)
	sprintf(tstr,"  DataID %s not filled on trial %d for processor %s",
		r->datid[i],ti,proc_name);
      else
	sprintf(tstr,"  DataID %s not filled on trial %d",r->datid[i],ti);
      printf("  *** WM_RESPONSE_CHECK_TRIAL_FULL  %s\n  *** Exiting.\n",tstr);
      mylogx(mylogf,"WM_RESPONSE_CHECK_TRIAL_FULL",tstr);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            WM_GET_SFORM_XN_YN_TN                          */
/*                                                                           */
/*****************************************************************************/
void wm_get_sform_xn_yn_tn(m,s,rsform,rxn,ryn,rtn)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     char **rsform;
     int *rxn;
     int *ryn;
     int *rtn;
{
  int xn,yn,tn;
  char *sform;

  sform = paramfile_get_char_param_default(s->ppl,"stim_form","3d");

  if (m->ppl != NULL)
    tn = paramfile_get_int_param_or_exit(m->ppl,"tn");
  else
    tn = onode_getpar_int_exit(m->o,"tn");

  if ((strcmp(sform,"3d")==0) || (strcmp(sform,"3d_b")==0) ||
      (strcmp(sform,"3rgb")==0) ||
      (strcmp(sform,"3c")==0) || (strcmp(sform,"3c_b")==0)){
    if (m->ppl != NULL){
      xn = paramfile_get_int_param_or_exit(m->ppl,"xn");
      yn = paramfile_get_int_param_or_exit(m->ppl,"yn");
    }else{
      xn = onode_getpar_int_exit(m->o,"xn");
      yn = onode_getpar_int_exit(m->o,"yn");
    }
  }else{
    xn = yn = 0;
  }

  *rsform = sform;
  *rxn = xn;
  *ryn = yn;
  *rtn = tn;
}
/**************************************-**************************************/
/*                                                                           */
/*                               WM_RESPONSE_INIT                            */
/*                                                                           */
/*****************************************************************************/
void wm_response_init(m,s,r,mylogf)
     struct model_struct *m;
     struct stim_struct *s;
     struct response_struct *r;
     char *mylogf;
{
  int i,j;
  int n,ntr,tn;
  float tscale;
  char ggstr[LONG_SLEN];

  // Check for old way, and exit
  if (m->ppl != NULL){
    r->ns = paramfile_count_values_param(m->ppl,"save_spikes") / 2;
    r->nf = paramfile_count_values_param(m->ppl,"save_data") / 2;
    if ((r->ns > 0) || (r->nf > 0)){
      mylog_exit(mylogf,"WM_RESPONSE_INIT - OLD FORMAT\n");
    }
  }

  wm_response_count(m,r,mylogf);   // set 'r->n'
  wm_response_config(m,r,mylogf);  // fill in details of each response

  mylog(mylogf,"  Response Configuration\n");
  for(i=0;i<r->n;i++){
    if ((r->n < 40) || (i < 5) || (i >= (r->n - 6))){
      if  (r->rformat[i] == 10){
	sprintf(ggstr,
	      "%4d %2s %8.1f %2d %8s %2d %2d %2d %8s %2d %2d %2d %14s %14s\n",
	      r->rformat[i],r->rtype[i],r->samp[i],r->ri[i],r->plname[i],
	      r->xi[i],r->yi[i],r->zi[i],r->plname1[i],
	      r->xi1[i],r->yi1[i],r->zi1[i],r->datid[i],r->nd_name[i]);
      }else{
	sprintf(ggstr,"%4d %2s %8.1f %2d %8s %2d %2d %2d %17s %14s %14s\n",
	      r->rformat[i],r->rtype[i],r->samp[i],r->ri[i],r->plname[i],
	      r->xi[i],r->yi[i],r->zi[i],"  ",r->datid[i],r->nd_name[i]);
      }
      mylog(mylogf,ggstr);
    }else if (i == 6){
      //strcpy(ggstr,"    .\n    .\n    .\n");
      strcpy(ggstr,"   ...     ...\n");
      mylog(mylogf,ggstr);
    }
  }

  // Time in seconds, to be used later to compute durations for records
  r->ts0 = 0.0;
  r->ttref = 0;
  if (m->ppl != NULL){
    tscale = paramfile_get_float_param_or_exit(m->ppl,"tscale");
    tn = paramfile_get_int_param_or_exit(m->ppl,"tn");
  }else{
    tscale = onode_getpar_flt_exit(m->o,"tscale");
    tn     = onode_getpar_int_exit(m->o,"tn");
    //  Allow the param 'tn_response" to override the value in the .nd file
    tn     = onode_getpar_int_dflt(m->o,"tn_response",tn);
  }
  r->tsn = (float)tn*tscale;   // frames * sec/frame ==> sec

  r->tsi  = 0;  // Set trial sequence pointer to zero
  r->gtsi = 0;  // Set global trial sequence pointer to zero
  r->tri  = 0;  // Set trial index pointer to zero

  //
  //  WYETH Aug 8, 2014:
  //  When this is running as an 'mm' slave, do we need to make storage for
  //  all combos?  Wouldn't we want 'ntr' to be set to 1, and simply not
  //  bother to have the slaves hold the data across multiple runs?
  //
  ntr = s->ntr * m->nrun;  // Multiply by the number of model configs to run

  r->dflag = get_zero_iarray(ntr);      // No trials done, used by master
  s->tribyx = get_const_iarray(ntr,-1); // Will hold 'tri' indices

  // Spikes
  n = r->ns;
  if (n > 0){
    r->s = (float ***)myalloc(n*sizeof(float **));
    r->cnt = (int **)myalloc(n*sizeof(int *));
    for(i=0;i<n;i++){
      r->s[i] = (float **)myalloc(ntr*sizeof(float *));
      r->cnt[i] = (int *)myalloc(ntr*sizeof(int));
      for(j=0;j<ntr;j++){
	r->s[i][j] = NULL;   // Marker, so we know if anything already stored
	r->cnt[i][j] = -1;   // So we know nothing has been stored
      }
    }
  }

  // Continuous
  n = r->nf;
  if (n > 0){
    r->f = (float ***)myalloc(n*sizeof(float **));
    r->fcnt = (int **)myalloc(n*sizeof(int *));
    for(i=0;i<n;i++){
      r->f[i] = (float **)myalloc(ntr*sizeof(float *));
      r->fcnt[i] = (int *)myalloc(ntr*sizeof(int));
      for(j=0;j<ntr;j++){
	r->f[i][j] = NULL;  // Marker, so we know if anything already stored
	r->fcnt[i][j] = -1; // so we know nothing has been stored
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             WM_GET_ND_T1_OUTFILE                          */
/*                                                                           */
/*  Return the outfile prefix, format, and names for .nd and .t1 files.      */
/*                                                                           */
/*****************************************************************************/
void wm_get_nd_t1_outfile(r,mylogf,myid,rpre,rform,rout,rout2)
     struct response_struct *r; // Response params
     char *mylogf;
     int myid;
     char **rpre;      // prefix for output file name
     char **rform;     // format for output file(s)
     char **rout;      // name of 1st outfile (.nd) or NULL
     char **rout2;     // name of 2nd outfile (.t1) or NULL
{
  char *outfile,*outfile2,*prefix,*format;

  prefix = paramfile_get_char_param_default(r->ppl,"outfile","zz");
  format = paramfile_get_char_param_default(r->ppl,"ndata_format","nd");

  if (strcmp(format,"nd")==0)
    outfile = paramfile_make_fname(r->ppl,"write_ndata",prefix,".nd");
  else if (strcmp(format,"nd,t1")==0){
    outfile = paramfile_make_fname(r->ppl,"write_ndata",prefix,".nd");
    outfile2 = paramfile_make_fname(r->ppl,"write_ndata",prefix,".t1");
  }else
    outfile = paramfile_make_fname(r->ppl,"write_ndata",prefix,".t1");

  *rpre  = prefix;
  *rform = format;
  *rout  = outfile;
  *rout2 = outfile2;
}
/**************************************-**************************************/
/*                                                                           */
/*                               WM_MAKE_NDTRIAL                             */
/*                                                                           */
/*****************************************************************************/
void wm_make_ndtrial(m,s,r,mylogf,myid,ti,nchan,nd,t)
     struct model_struct *m;
     struct stim_struct *s;
     struct response_struct *r;
     char mylogf[];
     int myid;
     int ti;
     int nchan;
     struct ndata_struct *nd;
     struct ndtrial_struct *t;  // Pointer to trial to fill
{
  int j,k;
  int ri,vi,mi,li;
  char *name,tstr[SLEN],*pval,*tval;

  //mylog(mylogf,"  WM_MAKE_NDTRIAL\n");

  t->tcode = ti;
  t->tref = my_rint_double(r->ttref * 1000.0);  // Convert to msec

  t->nparam = nd->nvar;
  t->pname = (char **)myalloc(t->nparam*sizeof(char *));
  t->pval  = (char **)myalloc(t->nparam*sizeof(char *));
  vi = s->tribyx[ti];

  for(j=0;j<s->nvar;j++){
    t->pname[j] = strdup(nd->vname[j]);
    t->pval[j]  = strdup(s->vval[vi][j]); // [ntrial][nvar]
  }
  if (m->vps != NULL){  // Moo var params
    mi = ti / s->ntr;  // Model index
    for(j=0;j<m->vps->m_nvar;j++){
      t->pname[s->nvar + j] = strdup(m->vps->m_name[j]);
      pval = m->vps->m_vval[mi][j];  // Convenient pointer [ntrial][nvar]

      //
      //  If there are varpar_lookup objects, check if we need to replace
      //  the actual value with the lookup key.  For example, a white-space
      //  separated list of entries could be replaced by a single key.
      //
      if (m->vps->vplk_n > 0){
	tval = mod_util_var_par_lookup_reverse(m->vps->m_name[j],pval,m->vps);
	if (tval != NULL){
	  //printf("  tval: %s   pval: %s\n",tval,pval);
	  pval = tval;  // Replace the value with the lookup key
	}
      }

      t->pval[s->nvar + j] = strdup(pval);
    }
    t->pname[s->nvar + m->vps->m_nvar] = strdup("MOO_var_i");
    sprintf(tstr,"%d",mi);
    t->pval[ s->nvar + m->vps->m_nvar] = strdup(tstr);
  }

  t->nrec = nchan;
  t->r = (struct ndrec_struct *)myalloc(t->nrec*sizeof(struct ndrec_struct));
  k = 0;
  for(j=0;j<r->n;j++){
    // Next line creates storage, otherwise all trials point to one name,
    //   and this creates problems when freeing.
    t->r[k].name = strdup(r->nd_name[j]);
    t->r[k].rcode = 0;
    t->r[k].sampling = r->samp[j];
    t->r[k].t0 = (int)(r->ts0 * r->samp[j]);
    t->r[k].tn = (int)(r->tsn * r->samp[j]);
    
    ri = r->ri[j]; // Index into data storage
    if (strcmp(r->rtype[j],"s")==0){
      if (r->cnt[ri][ti] < 0){
	sprintf(tstr,"No response stored for %s",r->datid[j]);
	mylogx(mylogf,"WM_MAKE_NDTRIAL",tstr);
      }
      t->r[k].rtype = 0;
      t->r[k].n = r->cnt[ri][ti];
      t->r[k].p = f2iarray_round(r->s[ri][ti],r->cnt[ri][ti]);
    }else if (strcmp(r->rtype[j],"f")==0){
      if (r->fcnt[ri][ti] < 0){
	sprintf(tstr,"No response stored for %s",r->datid[j]);
	mylogx(mylogf,"WM_MAKE_NDTRIAL",tstr);
      }
      t->r[k].rtype = 1;
      t->r[k].tn = r->fcnt[ri][ti];
      t->r[k].n = r->fcnt[ri][ti];
      t->r[k].x = copy_farray(r->f[ri][ti],r->fcnt[ri][ti]);
    }
    k += 1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               WM_MAKE_NDATA                               */
/*                                                                           */
/*****************************************************************************/
struct ndata_struct *wm_make_ndata(m,s,r,mylogf,myid)
     struct model_struct *m;
     struct stim_struct *s;
     struct response_struct *r;
     char mylogf[];
     int myid;
{
  int i,j,k,ti;
  int n,nchan,ri,vi,tot_ntr;
  float samp;
  char *name,tstr[SLEN];
  struct ndata_struct *nd;

  mylog(mylogf,"  WM_MAKE_NDATA\n");

  // Create and fill 'nd'
  nd = get_ndata();
  nd->class = strdup("MODEL");

  // Add spike channels and continuous channels
  nchan = r->ns + r->nf;
  if (nchan <= 0)
    mylog_exit(mylogf,
	       "WM_MAKE_NDATA  No response data to write ndata file\n");
  nd->nconst = s->ncon; // Const params
  nd->cname = (char **)myalloc(nd->nconst*sizeof(char *));
  nd->ctype = (char *)myalloc(nd->nconst*sizeof(char));
  nd->cval = (char **)myalloc(nd->nconst*sizeof(char *));
  for(i=0;i<nd->nconst;i++){
    nd->ctype[i] = s->ctype[i];
    nd->cname[i] = strdup(s->cname[i]);
    nd->cval[i] = strdup(s->cval[i]);
  }

  nd->nvar = s->nvar; // Stim var params
  if (m->vps != NULL){
    nd->nvar += m->vps->m_nvar; // Moo var params
    nd->nvar += 1;  // Add the generated param MOO_var_i
  }

  nd->vname = (char **)myalloc(nd->nvar*sizeof(char *));
  nd->vtype = (char *)myalloc(nd->nvar*sizeof(char));
  //for(i=0;i<nd->nvar;i++){
  for(i=0;i<s->nvar;i++){
    nd->vtype[i] = s->vtype[i];
    nd->vname[i] = strdup(s->vname[i]);
  }
  if (m->vps != NULL){
    for(i=0;i<m->vps->m_nvar;i++){
      nd->vtype[s->nvar + i] = 'c';  // WYETH FIX - Get type from var par recs
      nd->vname[s->nvar + i] = strdup(m->vps->m_name[i]);
    }
    nd->vtype[s->nvar + m->vps->m_nvar] = 'i';
    nd->vname[s->nvar + m->vps->m_nvar] = strdup("MOO_var_i");
  }

  nd->ntable = 0;

  tot_ntr = s->ntr * m->nrun;

  if (myid == 0){ // Master should count results in case of partial save
    n = 0;
    //for(i=0;i<s->ntr;i++){
    for(i=0;i<tot_ntr;i++){
      if (r->dflag[i] == 1)
	n += 1;
    }
  }else{
    //n = s->ntr;  // Assume all trials are complete
    n = tot_ntr;  // Assume all trials are complete
  }

  nd->ntrial = n;
  nd->t = (struct ndtrial_struct *)myalloc(n*sizeof(struct ndtrial_struct));
  ti = 0; // Trial save index will differ from 'i' for partial save
  //for(i=0;i<s->ntr;i++){
  for(i=0;i<tot_ntr;i++){
    if (r->dflag[i] == 1){
      wm_make_ndtrial(m,s,r,mylogf,myid,i,nchan,nd,&(nd->t[ti]));
      ti += 1;
    }
  }

  /*** Table of var values across all stimulus trials
  printf("s->ntr =%d\n",s->ntr);
  for(i=0;i<s->ntr;i++){
    for(j=0;j<s->nvar;j++){
      printf("vval[%d][%d] = %s\n",i,j,s->vval[i][j]);
    }
    }***/

  return nd;
}
/**************************************-**************************************/
/*                                                                           */
/*                                WM_STAT1_CUT                               */
/*                                                                           */
/*  Take an nd struct for a 'stat1d' stimulus and cut it into trials.        */
/*                                                                           */
/*****************************************************************************/
struct ndata_struct *wm_stat1_cut(mylogf,nd)
     char *mylogf;
     struct ndata_struct *nd;
{
  int i,j,k;
  int n,di,seqflag,seed,nblk,ndur,*ton,*toff,**ndx,nvar,ts,dt,ti,*t0,*tn;
  int *sdata,ns,dur,**s,*cnt,period,pflag;
  float tscale,tsep,sampling;
  char **vname,*vtype,***vval,tstr[SLEN],ggstr[LONG_SLEN];
  struct ndata_struct *ndcut;

  mylog(mylogf,"  WM_STAT1_CUT\n");

  pflag = 0;

  if (nd->ntrial != 1)
    exit_error("WM_STAT1_CUT","File must have only 1 trial to cut.");
  if (nd->t[0].nrec != 1)
    exit_error("WM_STAT1_CUT","Trial must have only 1 record to cut.");

  dur =    nd->t[0].r[0].tn;
  sdata =  nd->t[0].r[0].p;
  ns =     nd->t[0].r[0].n;
  sampling = nd->t[0].r[0].sampling;
  
  sprintf(ggstr,"    Single trial, duration %d, spikes %d\n",dur,ns);
  mylog(mylogf,ggstr);
  

  seed    = ndata_get_const_param_int(nd,"seed");
  seqflag = ndata_get_const_param_int(nd,"seq_type");
  nblk    = ndata_get_const_param_int(nd,"nblock");
  tsep    = ndata_get_const_param_float(nd,"tsep");

  tscale  = 0.001;  // Use msec

  ts = my_rint(tsep/tscale);  // Trial separation, msec
  
  stm_stat1d_get_seq(NULL,nblk,seqflag,seed,tscale,&ndur,&ton,&toff,&ndx);

  if (pflag){
    printf("  Stimuli:\n");
    for(i=0;i<ndur;i++){
      printf("    %2d %6d %6d\n",i,ton[i],toff[i]);
    }
  }

  n = nblk * ndur; // Maximum number of trials
  t0 = (int *)myalloc(n*sizeof(int));  // Trial start time
  tn = (int *)myalloc(n*sizeof(int));  // Trial duration

  s = (int **)myalloc(n*sizeof(int *));  // sarray
  cnt = (int *)myalloc(n*sizeof(int));

  // Make new ndata structure, use old const params
  /***
  ndcut = get_ndata();
  init_set_const_param(ndcut,nd->nconst,nd->cname,nd->ctype,nd->cval);
  ***/

  // Set up var param info
  nvar = 2;
  vname = (char **)myalloc(nvar*sizeof(char *));
  vtype = (char *)myalloc(nvar*sizeof(char));
  vname[0] = strdup("stn");   vtype[0] = 'f';
  vname[1] = strdup("stn2");  vtype[1] = 'f';
  vval = get_2d_pointer_carray(n,nvar);
  k = 0; // Trial counter
  ti = 0;  // Start time of curren stimulus
  period = 0;  // Find the maximum trial length
  for(i=0;i<nblk;i++){
    for(j=0;j<ndur;j++){
      di = ndx[i][j];
      sprintf(tstr,"%.2f",(float)(ton[di]) * tscale);
      vval[k][0] = strdup(tstr);
      sprintf(tstr,"%.2f",(float)(toff[di]) * tscale);
      vval[k][1] = strdup(tstr);

      t0[k] = ti;                    // Start of this stimulus
      dt = ts + ton[di] + toff[di];  // length of this stimulus
      tn[k] = dt;                    // Stimulus duration
      ti += dt;                      // Start of next stimulus

      if (ti > dur){ // WYETH - exit if not all stimuli have been used
	mylog_exit(mylogf,"");
	sprintf(ggstr,"*** WM_STAT1_CUT  Stimuli truncated");
      }

      // Get spikes and count for kth trial
      s[k] = extract_spikes(sdata,ns,t0[k],tn[k],&(cnt[k]),1);
      if (tn[k] > period)
	period = tn[k];

      k += 1;
    }
  }

  if (pflag){
    printf("  Var Params:\n");
    for(i=0;i<n;i++){
      printf("    trial %2d %8d %6d  ",i,t0[i],tn[i]);
      for(j=0;j<nvar;j++){
	printf("  %s",vval[i][j]);
      }
      printf("\n");
    }
  }

  ndcut = make_ndata_sarray(s,cnt,n,period,nd->nconst,nvar,nd->cname,nd->ctype,
			    nd->cval,vname,vtype,vval,nd->class,sampling);

  myfree(t0); myfree(tn);
  myfree(ton); myfree(toff); free_2d_iarray(ndx,nblk);
  free_sarray(s,cnt,n);

  myfree(vtype);
  free_2d_pointer_carray(vval,n,nvar);
  free_2d_carray(vname,nvar);

  return ndcut;
}
/**************************************-**************************************/
/*                                                                           */
/*                              WM_WRITE_RESPONSE                            */
/*                                                                           */
/*  Write model response to output file.                                     */
/*                                                                           */
/*****************************************************************************/
void wm_write_response(m,s,r,mylogf,myid)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     char *mylogf;
     int myid;
{
  int cutflag;
  char *outfile,*prefix,*format,tstr[SLEN],ggstr[LONG_SLEN],*stim_type;
  char *outfile2;
  struct ndata_struct *nd,*ndcut;
  struct param_pair *pp;

  //
  //  Get the outfile 'prefix', 'format', and outfile names.
  //
  wm_get_nd_t1_outfile(r,mylogf,myid,&prefix,&format,&outfile,&outfile2);

  /*
  prefix = paramfile_get_char_param_default(r->ppl,"outfile","zz");
  format = paramfile_get_char_param_default(r->ppl,"ndata_format","nd");

  if (strcmp(format,"nd")==0)
    outfile = paramfile_make_fname(r->ppl,"write_ndata",prefix,".nd");
  else if (strcmp(format,"nd,t1")==0){
    outfile = paramfile_make_fname(r->ppl,"write_ndata",prefix,".nd");
    outfile2 = paramfile_make_fname(r->ppl,"write_ndata",prefix,".t1");
  }else
    outfile = paramfile_make_fname(r->ppl,"write_ndata",prefix,".t1");
  */

  if (outfile != (char *)NULL){
    nd = wm_make_ndata(m,s,r,mylogf,myid);
    if (nd->ntrial < s->ntr){
      sprintf(tstr,"%s.partial",outfile);
      myfree(outfile);
      outfile = strdup(tstr);
      printf("    Partial results:  %d of %d trials\n",nd->ntrial,s->ntr);
    }

    //
    //  Check if we need to cut a long trial into multiple trials
    //
    cutflag = paramfile_get_int_param_default(r->ppl,"ndata_cut",0);
    if (cutflag == 1){
      if (!ndata_get_const_param(nd,"stim_type",&stim_type)){
	mylog_exit(mylogf,"WM_WRITE_RESPONSE  Cannot find stim_type");
      }else{
	ndcut = NULL;
	if (strcmp(stim_type,"stat1")==0){
	  ndcut = wm_stat1_cut(mylogf,nd);
	  //ndcut = nda_util_stat1_cut(mylogf,nd);
	}else{
	  sprintf(ggstr,"  stim_type = %s\n",stim_type);
	  mylog(mylogf,ggstr);
	  mylog_exit(mylogf,"WM_WRITE_RESPONSE  Do not know how to cut");
	}
      }
      free_ndata(nd);  // Free original file
      nd = ndcut;      // Replace with file cut into trials
      myfree(stim_type);
    }

    //
    //  Check for decision computation
    //
    //                    1        2       3    4      5  6
    //    decision_gen choice prefix_diff ds_ dsopp_  t0 tn
    //

    pp = paramfile_get_param_pointer(r->ppl,"decision_gen");
    if (pp != NULL){
      int nval,t0,tn;
      char *ch_name,*ctype,*name_add,*name_sub;
      float *dval,*dvadd,*dvsub;

      mylog(mylogf,"  WM_WRITE_RESPONSE  Writing decision file\n");

      nval = paramfile_count_values_pointer(pp);
      if (nval != 6)
	exit_error("WM_WRITE_RESPONSE","decision_gen has bad value list");
    
      ch_name  = paramfile_get_nth_char_param_pp_or_exit(pp,0);
      ctype    = paramfile_get_nth_char_param_pp_or_exit(pp,1);
      name_add = paramfile_get_nth_char_param_pp_or_exit(pp,2);
      name_sub = paramfile_get_nth_char_param_pp_or_exit(pp,3);
      t0       = paramfile_get_nth_int_param_ptr_or_exit(pp,4);
      tn       = paramfile_get_nth_int_param_ptr_or_exit(pp,5);

      /*
	printf("  ch_name   %s\n",ch_name);
	printf("  ctype     %s\n",ctype);
	printf("  name_add  %s\n",name_add);
	printf("  name_sub  %s\n",name_sub);
	printf("  t0        %d\n",t0);
	printf("  tn        %d\n",tn);
      */
      
      ndata_util_choice_01(nd,ctype,t0,tn,name_add,name_sub,ch_name,
			   &dval,&dvadd,&dvsub);
      myfree(dval); myfree(dvadd); myfree (dvsub);
    }
    

    // Write appropriate file format
    if (strcmp(format,"nd")==0)
      write_ndata(outfile,nd);
    else if (strcmp(format,"t1")==0)
      write_ndata_t1(outfile,"iModel",nd,4,0); // 4-digits, 0-all records
    else if (strcmp(format,"nd,t1")==0){
      write_ndata(outfile,nd);
      write_ndata_t1(outfile2,"iModel",nd,4,0); // 4-digits, 0-all records
    }else
      mylog_exit(mylogf,"WM_WRITE_RESPONSE  Unknown ndata format\n");

    free_ndata(nd);
    myfree(outfile);
    myfree(format);
  }
  myfree(prefix);
}
/**************************************-**************************************/
/*                                                                           */
/*                            WMU_CALIBRATION_INIT                           */
/*                                                                           */
/*  Initialize data for calibration.                                         */
/*  Return '1' if the mode is "generate".                                    */
/*                                                                           */
/*****************************************************************************/
int wmu_calibration_init(m,r,mylogf)
     struct model_struct *m;
     struct response_struct *r;
     char *mylogf;
{
  FILE *fopen(),*fin;
  int i,j;
  int nm,nv,si,mvpn,nr;
  char *calmode,*calfile,*calconp,tstr[LONG_SLEN],**calval,**mvpl;
  struct onode *co;

  mylog(mylogf,"  WMU_CALIBRATION_INIT\n");

  co = onode_child_get_unique(m->o,"calibrate");
  if (co == NULL)
    return 0;  // because there is no <calibrate> object

  //
  //  Check the <calibrate> "mode" and act accordingly
  //
  calmode = onode_getpar_chr_ptr_exit(co,"mode");
  mylog_cval(mylogf,"    mode =",calmode);
  
  if (strcmp(calmode,"off")==0)
    return 0;  // All calibration processing is turned off
  else if (strcmp(calmode,"generate")==0){
    if (r->n != 1){
      printf(" r->n = %d\n",r->n);
      exit_error("WMU_CALIBRATION_INIT",
		 "The .rsp file must request exactly one response");
    }
    return 1;  // No reading to do here; signal the need to write "cal_file"
  }else if (strcmp(calmode,"run")!=0)
    exit_error("WMU_CALIBRATION_INIT","Unknown <calibrate> 'mode'");


  calfile = onode_getpar_chr_ptr_exit(co,"cal_file");
  calconp = onode_getpar_chr_ptr_exit(co,"control_param");

  if ((fin = fopen(calfile,"r")) == NULL){
    printf("  *** File %s\n",calfile);
    exit_error("WMU_CALIBRATION_INIT","Cannot open file");
  }
  nr = fscanf(fin,"%*s %d",&nm);  // Number of model configurations
  nr = fscanf(fin,"%*s %d",&nv);  // Number of moo variables
  for(i=0;i<nv;i++){
    nr = fscanf(fin,"%s",tstr);  // Number of model configurations
    //printf("tstr = %s\n",tstr);
    // *** CHECK IF THESE MATCH THE MOO_VAR ???
  }
  nr = fscanf(fin,"%*s");  // Read the constant, "PARAMETERS"
  for(i=0;i<nm;i++){
    for(j=0;j<nv;j++){
      nr = fscanf(fin,"%s",tstr);  // Number of model configurations
      //printf("tstr = %s\n",tstr);
      // *** CHECK IF THESE MATCH THE MOO_VAR ???
    }
  }
  nr = fscanf(fin,"%s",tstr);  // Read the constant, "CALIBRATION"
  if (strcmp(tstr,"CALIBRATION")!=0){
    printf("  string is:  %s\n",tstr);
    exit_error("WMU_CALIBRATION_INIT","Expecting 'CALIBRATION'");
  }

  //
  //  Read the calibrated parameter values
  //
  calval = (char **)myalloc(nm*sizeof(char *));
  for(i=0;i<nm;i++){
    nr = fscanf(fin,"%s",tstr);  // Number of model configurations
    calval[i] = strdup(tstr);
  }
  fclose(fin);


  m->calib_pval = calval;  // Save in main model strcuture

  //
  //  Verify that the calibration parameter is not a MOO VAR param
  //
  mvpn = m->vps->m_nvar;  // Number of MOO VAR param names
  mvpl = m->vps->m_name;  // List of MOO VAR param names
  si = search_2d_carray(mvpl,calconp,mvpn);
  if (si != -1)
    exit_error("WMU_CALIBRATION_INIT",
	       "Calibration 'control_param' cannot be a <var_param>");

  //printf("calconp = %s\n",calconp);

  //exit_error("WM_CALIBRATION_INIT","No <calibrate> object in .moo file");
  //exit_error("WM_CALIBRATION_INIT","Under development");

  return 0;  // Return value indicates not to generate a new "cal_file"
}
/**************************************-**************************************/
/*                                                                           */
/*                            WMU_CALIBRATION_WRITE                          */
/*                                                                           */
/*  Write output for calibration file.                                       */
/*                                                                           */
/*****************************************************************************/
void wmu_calibration_write(m,s,r,mylogf)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     char *mylogf;
{
  FILE *fopen(),*fout;
  int i,j;
  int nm,nv;
  char *calfile;
  struct onode *co;

  mylog(mylogf,"  WMU_CALIBRATION_WRITE\n");

  nm = m->nrun;
  if (nm > 1){
    if (m->vps->m_n != nm)
      exit_error("WMU_CALIBRATION_WRITE","Mis-match number of models");
    nv = m->vps->m_nvar;
    printf("    for %d model configurations\n",nm);
  }

  co = onode_child_get_unique(m->o,"calibrate");
  if (co == NULL)
    exit_error("WMU_CALIBRATION_WRITE","No <calibrate> object in .moo file");

  calfile = onode_getpar_chr_ptr_exit(co,"cal_file");

  if (file_exists(calfile)){
    printf("  WARNING:  Over-writing existing cal_file.\n");
  }

  if ((fout = fopen(calfile,"w")) == NULL){
    printf("  *** File %s\n",calfile);
    exit_error("WMU_CALIBRATION_WRITE","Cannot open file");
  }
  fprintf(fout,"MODELS %d\n",nm);
  if (nm > 1){
    fprintf(fout,"VARIABLES %d\n",nv);

    for(j=0;j<nv;j++)  // for each model variation
      fprintf(fout,"%s\n",m->vps->m_name[j]);

    fprintf(fout,"PARAMETERS\n");

    for(i=0;i<nm;i++){  // for each model variation
      for(j=0;j<nv;j++){  // for each parameter
	fprintf(fout," %s",m->vps->m_vval[i][j]);
      }
      fprintf(fout,"\n");
    }
  }

  //
  //  Write Calibration values
  //
  fprintf(fout,"CALIBRATION\n");
  {
    int ti,t0,tn;
    int ns;     // Number of stimuli per model
    float *d;   // [ns]
    float *td;  // Data array
    float dmax,targv;

    t0 = onode_getpar_int_exit(co,"resp_t0");
    tn = onode_getpar_int_exit(co,"resp_tn");

    targv = onode_getpar_flt_exit(co,"target_value");

    printf("  t0,n = %d %d   targv = %f\n",t0,tn,targv);

    if (t0+tn > r->fcnt[0][0])  // WYETH - using first float response
      exit_error("WMU_CALIBRATION_WRITE","Response window too long");

    ns = s->ntr;  // Number of stimuli for each model

    d = (float *)myalloc(ns*sizeof(float));

    ti = 0;  // trial index
    for(i=0;i<nm;i++){   // for each model variation
      for(j=0;j<ns;j++){  // for each stimulus
	td = r->f[0][ti];  // Pointer to response data
	d[j] = mean_farray(&(td[t0]),tn);
	ti += 1;
      }
      dmax = max_of_farray(d,ns);
      fprintf(fout,"%f\n",targv/dmax);
    }
    myfree(d);
  }

  fclose(fout);

  //exit_error("WMU_CALIBRATION_WRITE","Calibrate is under development");
}
/**************************************-**************************************/
/*                                                                           */
/*                               WMU_WMPI_STIM                               */
/*                                                                           */
/*****************************************************************************/
void wmu_wmpi_stim(s,k,ktr,xn,yn,tn,sscale,tscale,tsamp,mylogf)
     struct stim_struct *s;  // Stimulus params
     int k;                  // Stimulus number among unique stimuli
     int ktr;                // Trial number
     int xn,yn,tn;
     float sscale,tscale;
     int tsamp;
     char *mylogf;
{
  FILE *fopen(),*fin;
  int i,j,t;
  int nr;
  float ***d;
  char *wmpi_path,*wmpi_pref,fname[SLEN],ggstr[LONG_SLEN];

  wmpi_path = param_getc_exit(s->ppl,"stim_wmpi_path");
  wmpi_pref = param_getc_exit(s->ppl,"stim_wmpi_prefix");

  sprintf(fname,"%s/%s%d",wmpi_path,wmpi_pref,k);
  sprintf(ggstr,"  [WMPI reading from file: %s]\n",fname);
  mylog(mylogf,ggstr);

  d = f3tensor(1,xn,1,yn,1,tn);

  if ((fin = fopen(fname,"r")) == NULL){
    printf("  *** File %s\n",fname);
    exit_error("WMU_WMPI_STIM","Cannot open file");
  }

  for(t=1;t<=tn;t++) // For each frame
    for(j=yn;j>=1;j--)  // Read lines (top line comes first)
      for(i=1;i<=xn;i++){
	nr = fread((char *)&(d[i][j][t]),sizeof(float),1,fin);
	if (nr != 1){
	  printf("(t=%d,y=%d,x=%d) nr = %d\n",t,j,i,nr);
	  exit_error("WMU_WMPI_STIM","Failed to read float data");
	}
      }

  fclose(fin);

  s->d = d;

  myfree(wmpi_path);
  myfree(wmpi_pref);
}
/**************************************-**************************************/
/*                                                                           */
/*                        GET_STIMULUS_DATA_TUNING_CURVE                     */
/*                                                                           */
/*  NOTES                                                                    */
/*  - Data is returned in 's' struct.                                        */
/*  - Special case for "dirmod" removed, see ./archive/wm.c.002.             */
/*                                                                           */
/*****************************************************************************/
void get_stimulus_data_tuning_curve(m,s,r,k,ktr,mylogf)
     struct model_struct *m;     // Model params.
     struct stim_struct *s;      // Stimulus params.
     struct response_struct *r;  // Response params.
     int k;                      // Stimulus number among unique stimuli
     int ktr;                    // Trial number
     char *mylogf;
{
  int i;
  int tsamp,xn,yn,tn,write_nd,ndstim_flag,wmpi_flag;
  float sscale,tscale,stim_samp;
  char *sform,tname[2048];

  s->stimno = k;

  //
  //  Fill 's->name' with a descriptive name for this stimulus, or 'stim0'
  //                 if no var params.
  //
  tname[0] = '\0';
  if (s->nvar > 0){
    for(i=0;i<s->nvar;i++){
      strcat(tname,s->vname[i]);
      strcat(tname,"_");
      strcat(tname,s->vval[ktr][i]);
      if (i < (s->nvar-1))
	strcat(tname,"_");
    }
  }else
    strcpy(tname,"stim0");  // Default name when no var params
  if (s->name != NULL)
    myfree(s->name);
  s->name = strdup(tname);


  //
  //  Get 'sscale' and 'tscale'
  //
  if (m->ppl != NULL)
    tscale = paramfile_get_float_param_or_exit(m->ppl,"tscale");
  else
    tscale = onode_getpar_flt_exit(m->o,"tscale");

  wm_get_sform_xn_yn_tn(m,s,&sform,&xn,&yn,&tn);
  if ((strcmp(sform,"3d")==0) || (strcmp(sform,"3d_b")==0) ||
      (strcmp(sform,"3c")==0) || (strcmp(sform,"3c_b")==0)){
    if (m->ppl != NULL)
      sscale = paramfile_get_float_param_or_exit(m->ppl,"sscale");
    else
      sscale = onode_getpar_flt_exit(m->o,"sscale");
  }else
    sscale = 0.0;

  stim_samp = paramfile_get_float_param_or_exit(s->ppl,"stim_samp");
  tsamp = stm_get_tsamp(tscale,stim_samp);

  wmpi_flag = param_geti_dflt(s->ppl,"stim_wmpi",0);  // WM Python interface

  if (strcmp(sform,"3d")==0){
    // WYETH - where is s->d freed?
    if (wmpi_flag == 1)
      wmu_wmpi_stim(s,k,ktr,xn,yn,tn,sscale,tscale,tsamp,mylogf);
    else
      mod_util_get_3d_stim(s,ktr,xn,yn,tn,sscale,tscale,tsamp,mylogf);
    // Now s->d has stimulus, s->d_r is NULL
    //   Or, s->d could be null
  }else if (strcmp(sform,"3c")==0){
    mod_util_get_3c_stim(s,ktr,xn,yn,tn,sscale,tscale,tsamp,mylogf);
  }else if (strcmp(sform,"3rgb")==0){
    mod_util_get_3rgb_stim(s,ktr,xn,yn,tn,sscale,tscale,tsamp,mylogf);
  }else if (strcmp(sform,"3d_b")==0){
    mod_util_get_3d_stim_binoc(s,ktr,xn,yn,tn,sscale,tscale,tsamp,mylogf);
    // Now both s->d and s->d_r should be set
  }else if (strcmp(sform,"1d")==0){
    s->d1 = mod_util_get_1d_stim(s,ktr,tn,tscale,tsamp,mylogf);
  }else
    mylog_exit(mylogf,"GET_STIMULUS_DATA_TUNING_CURVE stim_form invalid\n");

  // WYETH - for writing stimulus into ndata file???
  write_nd = paramfile_get_int_param_default(r->ppl,"write_nd",0);
  if (write_nd)
    ndstim_flag = paramfile_get_int_param_default(r->ppl,"write_nd_stim",0);

  if (write_nd && ndstim_flag){
    if (0){ // WYETH OLD, FIX
      float *tseq;
      int ntseq;

      mylog(mylogf,"***** WYETH - writing 10 zeros instead of random sequ.\n");
      ntseq = 10;
      tseq = get_zero_farray(ntseq);
      //wm_ndata_store_trial_data(tseq,(int *)NULL,ntseq,1); // 1-stim
      myfree(tseq);
    }else
      mylog_exit(mylogf,
		 "GET_STIMULUS_DATA_TUNING_CURVE  Expecting tseq data\n");
  }

  myfree(sform);
}
/**************************************-**************************************/
/*                                                                           */
/*                             WM_MOO_VAR_PAR_SET                            */
/*                                                                           */
/*  Update the parameters of the model.                                      */
/*                                                                           */
/*****************************************************************************/
void wm_moo_var_par_set(mylogf,m,k)
     char *mylogf;
     struct model_struct *m;    // Model params
     int k;                     // index for model config
{
  int i;
  int n,flag;
  char ggstr[SLEN],*pval,*pname;

  mylog(mylogf,"  WM_MOO_VAR_PAR_SET\n");

  n = m->vps->m_nvar;  // Number of moo params varying
  for(i=0;i<n;i++){
    pval  = m->vps->m_vval[k][i];  // Value of param
    pname = m->vps->m_name[i];     // Name (ostr) of param

    //printf("================= %s =============== %s =================\n",
    //pname,pval);

    flag = onode_update_value(m->o,pname,pval);
    if (flag == 0){
      sprintf(ggstr,"*** WM_MOO_VAR_PAR_SET Cannot update: %s",pname);
      mylog_exit(mylogf,ggstr);
    }
  }

  //  Check for <calibrate> parameter values

  if (m->calib_pval != NULL){

    // Get the name of the calibrated parameter to be set
    pname = onode_getpar_child_chr_ptr_exit(m->o,"calibrate","control_param");

    // Set the parameter
    flag = onode_update_value(m->o,pname,m->calib_pval[k]);
    if (flag == 0){
      sprintf(ggstr,"*** WM_MOO_VAR_PAR_SET Cannot update calib: %s",pname);
      mylog_exit(mylogf,ggstr);
    }

    sprintf(ggstr,"    Calibrate:  %s <== %s\n",pname,m->calib_pval[k]);
    mylog(mylogf,ggstr);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               WM_SLAVE_ACTION                             */
/*                                                                           */
/*****************************************************************************/
void wm_slave_action(m,s,r,actflag,mylogf)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     int actflag;               // 0-prep, 1-run, 2-fit,
                                //  -1-cleanup, 10-mar, 20-rspgen
     char *mylogf;
{
  int k;
  char *modtype;

  m->action = actflag;  // Will be used for MAR and RSPGEN

  if (actflag == 1)
    k = s->repno;  // WYETH - phase out the use of 'k' here
  else
    k = -1;

  if (m->ppl != NULL)
    modtype = paramfile_get_char_param_or_exit(m->ppl,"mod_type");
  else
    modtype = onode_getpar_chr_exit(m->o,"mod_type");

  if (strcmp(modtype,"lin_01")==0)
    mod_lin_run_lin_01(m,s,r,k,actflag);
  else if (strcmp(modtype,"dog_01")==0)
    mod_dog_run_dog_01(m,s,r,k,actflag);
  else if (strcmp(modtype,"simp_01")==0)
    mod_dog_run_simp_01(m,s,r,k,actflag);
  else if (strcmp(modtype,"comp_01")==0)
    mod_dog_run_comp_01(m,s,r,k,actflag);
  else if (strcmp(modtype,"ds_01")==0)
    mod_dog_run_ds_01(m,s,r,k,actflag);
  else if (strcmp(modtype,"ds_02")==0)
    mod_dog_run_ds_02(m,s,r,k,actflag);
  else if (strcmp(modtype,"ifc_01")==0)
    ifc_util_run_ifc_01(m,s,r,k,actflag);
  else if (strcmp(modtype,"ifc_02")==0)
    ifc_util_run_ifc_02(m,s,r,k,actflag);
  else if (strcmp(modtype,"ifc_02_02")==0)
    ifc_util_run_ifc_02_02(m,s,r,k,actflag);
  else if (strcmp(modtype,"ifc_03")==0)
    ifc_util_run_ifc_03(m,s,r,k,actflag);

  //
  //  mod_me_util
  //
  else if (strcmp(modtype,"me_gabor_01")==0)
    mod_me_run_me_gabor_01(m,s,r,k,actflag);
  else if (strcmp(modtype,"me_ab_01")==0)
    mod_me_run_me_gabor_01(m,s,r,k,actflag);
  else if (strcmp(modtype,"me_v5")==0)
    mod_me_run_me_v5(m,s,r,k,actflag);
  else if (strcmp(modtype,"random_filter")==0)
    mod_me_run_me_gabor_01(m,s,r,k,actflag);
  else if (strcmp(modtype,"binoc_filter")==0)
    mod_me_run_me_gabor_01(m,s,r,k,actflag);
  else if (strcmp(modtype,"filter_xyt")==0)
    mod_me_run_filter_xyt(m,s,r,k,actflag);
  else if (strcmp(modtype,"gabor_comp")==0)
    mod_me_run_gabor_comp(m,s,r,k,actflag);
  else if (strcmp(modtype,"rd_2gabor")==0)
    mod_me_run_gabor_comp(m,s,r,k,actflag);
  else if (strcmp(modtype,"ds_simp_01")==0)
    mod_me_run_ds_simp_01(m,s,r,k,actflag);
  else if (strcmp(modtype,"ds_liv3")==0)
    mod_me_run_ds_liv3(m,s,r,k,actflag);
  else if (strcmp(modtype,"mem")==0)
    mod_me_run_mem(m,s,r,k,actflag);
  else if (strcmp(modtype,"reichardt")==0)
    mod_me_run_reich(m,s,r,k,actflag);

  else if (strcmp(modtype,"mesh_rgc_01")==0)
    mod_mesh_run_mesh_rgc_01(m,s,r,k,actflag);
  else if (strcmp(modtype,"mesh_mag_01")==0)
    mod_mesh_run_mesh_mag_01(m,s,r,k,actflag);
  else if (strcmp(modtype,"pop")==0)
    mod_pop_run(m,s,r,actflag);
  else if (strcmp(modtype,"vhf")==0){
    //exit_error("WM_SLAVE_ACTION","MOD_VHF_RUN commented out to avoid MPI");
    mod_vhf_run_01(m,s,r,actflag);
  }else if (strcmp(modtype,"caffe_net")==0)
    mod_dcn_run_01(m,s,r,actflag);
  else if (strcmp(modtype,"spop_01")==0)
    mylog_exit(mylogf,"WM_SLAVE_ACTION  'spop_01' deactivated\n");
  else if (strcmp(modtype,"spop_02")==0)
    mylog_exit(mylogf,"WM_SLAVE_ACTION  'spop_02' deactivated\n");
  else if (strcmp(modtype,"xytcorr")==0)
    mod_test_run_xytcorr01(m,s,r,actflag);
  else if (strcmp(modtype,"x")==0)
    mod_x_run_01(m,s,r,actflag);
  else if (strcmp(modtype,"srf")==0)
    mod_srf_run_01(m,s,r,actflag);
  else if (strcmp(modtype,"wave_fbars97")==0)
    mod_wav_run_fbars(m,s,r,actflag);
  else
    mylog_exit(mylogf,"WM_SLAVE_ACTION  Unknown model type\n");

  myfree(modtype);
  //mylog(mylogf,"  WM_SLAVE_ACTION_done\n");
}
