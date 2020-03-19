/*****************************************************************************/
/*                                                                           */
/*  wm.c                                                                     */
/*  wyeth bair                                                               */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*  Hardware issues                                                          */
/*  (1) I found that                                                         */
/*      s[i] = (int)(tf * (float)scale);                                     */
/*      converts to double in the multiplication on 'koch' but not on my     */
/*      MacBook Pro, thus given a different value.  I have changed this to   */
/*      s[i] = (int)(tf * (double)scale);                                    */
/*      to try to achieve consistency.  This occurred in "get_seeds" in      */
/*      myrand_util.c, where 'tf' was the float value from ran2.             */
/*                                                                           */
/*****************************************************************************

  mpicc -O2 -o ../bin/wm wm.c -L../lib -I../libc \
  -lwm_util -lmod_pop_util -lmod_dog_util -lmod_wav_util -lmod_vhf_util \
  -lmod_dcn_util -lmod_mesh_util -lmod_me_util -lifc_util -lpop_util \
  -lmod_conn_util -lmod_lin_util -lmod_test_util \
  -lmod_x_util -lmod_srf_util -lpop_cell_util -lmm_util \
  -lmod_util -lpop_low_util -lstim_util -lstm_util -lretina_util \
  -lmctrl_util -lkernel_util -lnoise_util -lsig_util -lfft_util \
  -lparamfile_util -lndata_util -lspike_util \
  -lmin_util -ldata_util -lmyrand_util -lfunc_util \
  -lplot_util -lfarray_util -liarray_util -lcarray_util -lmisc_util \
  -lnr_util -lmy_util -lm

******************************************************************************/
//  -lnr_util -lmy_util -ldmalloc -lm -lGL -lX11   // USE for dmalloc debuging
//  -lnr_util -lmy_util -lm -lGL -lX11             // REGULAR LINE
//   -L/usr/X11R6/lib -I/usr/X11R6/include \
//
//  2019 WYETH - removed these to eliminate OpenGL, X11 dependence
//   -lmod_gui_util -lmod_mon_util -lmod_gr_util
//   -lwygit_util -lglplot_util
//   -L/usr/X11/lib -I/usr/X11/include \
//   -lGL -lX11


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>


#include "my_util.h"
#include "nr_util.h"
#include "carray_util.h"
#include "iarray_util.h"
#include "farray_util.h"
#include "myrand_util.h"
#include "data_util.h"
#include "paramfile_util.h"
#include "ndata_util.h"
#include "mctrl_util.h"
#include "mod_conn_util.h"
#include "mod_util.h"
#include "wm_util.h"
#include "paramfile.h"
#include "ndata.h"
#include "mod.h"
#include "ifc.h"

char *mylogf = NULL;

/**************************************-**************************************/
/*                                                                           */
/*                           WM_MAKE_NAME_FOR_TRIAL                          */
/*                                                                           */
/*****************************************************************************/
char *wm_make_name_for_trial(s,r)
     struct stim_struct *s;
     struct response_struct *r;
{
  int i,j,k;
  int n;
  char *tname,*tval;
  char tstr[LONG_SLEN];

  n = s->nvar;

  // For stim_type "null", Start name with "Null"
  if ((s->d == NULL) && (s->dc == NULL) && (s->drgb_l == NULL)){
    tstr[0] = 'N';
    tstr[1] = 'u';
    tstr[2] = 'l';
    tstr[3] = 'l';
    k = 4;
  }else{
    k = 0;
  }

  for(i=0;i<n;i++){ // For each var param
    tname = s->vname[i];
    tval  = s->vval[r->tri][i];
    for(j=0;j<strlen(tname);j++){
      tstr[k] = tname[j];
      k += 1;
    }
    tstr[k] = '=';
    k += 1;
    for(j=0;j<strlen(tval);j++){
      tstr[k] = tval[j];
      k += 1;
    }
    if (i < (n-1)){
      tstr[k] = '_';
      k += 1;
    }
  }
  tstr[k] = '\0';

  return strdup(tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                              WM_WRITE_3D_STIM                             */
/*                                                                           */
/*  This routine should be used to write all stimuli.                        */
/*                                                                           */
/*****************************************************************************/
void wm_write_3d_stim(stimfile,outfile,arglist,narg,stim_num,frame_num,fstflag)
     char stimfile[],outfile[];
     char **arglist;
     int narg;
     int stim_num;    // Stimulus number,  0...
     int frame_num;   // Frame number, 0...
     int fstflag;     // 0-write_3d, 1-fst
{
  int i,j,k;
  int xn,yn,tn,tsamp,dxn,***idata,fst_binoc,fst_tn,wmpi_flag;
  float sscale,tscale,stim_samp,***data,**t,dxa;
  char *sform;
  struct stim_struct *s;  // Stimulus params

  printf("  WM_WRITE_3D_STIM\n");

  // Initialize stimulus struct
  s = (struct stim_struct *)myalloc(sizeof(struct stim_struct));
  s->ppl = paramfile_get_initial_params_update_only(stimfile,arglist,narg,1);
  strcpy(s->paramfile,stimfile);
  s->name = (char *)NULL;
  s->stimno = -1; // WYETH - added 2007 Jan

  if (!paramfile_list_update(s->ppl,"stim_nrpt","1"))
    exit_error("WM_WRITE_3D_STIM","Could not set 'stim_nrpt' to 1\n");

  wm_prep_var_const_gen(NULL,s,mylogf,NULL,0,NULL); // Will set model seed 1777

  // Check if 'stim_num' is valid
  if ((stim_num < 0) || (stim_num > (s->ntr-1))){
    printf("  Valid stimulus range is 0..%d\n",s->ntr-1);
    exit_error("WM_WRITE_3D_STIM","Stimulus number out of range");
  }

  xn = paramfile_get_int_param_or_exit(s->ppl,"stim_frame_xn");
  yn = paramfile_get_int_param_or_exit(s->ppl,"stim_frame_yn");
  tn = paramfile_get_int_param_or_exit(s->ppl,"stim_frame_tn");
  sscale = paramfile_get_float_param_or_exit(s->ppl,"stim_frame_sscale");
  tscale = paramfile_get_float_param_or_exit(s->ppl,"stim_frame_tscale");
  stim_samp = paramfile_get_float_param_or_exit(s->ppl,"stim_samp");
  tsamp = my_rint(1.0/(tscale*stim_samp)); // Time units per stimulus frame
  if (tsamp <= 0)
    tsamp = 1; // Minimum value for tsamp

  fst_binoc = 0;  // Assume not binocular, may be changed below
  fst_tn = tn;    // Assume, but change below as needed

  printf("  Will write data for stimulus %d (0..%d)\n",stim_num,s->ntr-1);

  wmpi_flag = param_geti_dflt(s->ppl,"stim_wmpi",0);  // WM Python interface

  sform = paramfile_get_char_param_default(s->ppl,"stim_form","3d");
  if (strcmp(sform,"3d")==0){
    // *** WYETH - THIS CODE DUPLICATED IN 'get_stimulus_data_tuning_curve'
    //    *** we should eliminate this duplication

    if (wmpi_flag == 1)
      // Send trial number as -1
      wmu_wmpi_stim(s,stim_num,-1,xn,yn,tn,sscale,tscale,tsamp,mylogf);
    else
      mod_util_get_3d_stim(s,stim_num,xn,yn,tn,sscale,tscale,tsamp,mylogf);
    data = s->d;
  }else if (strcmp(sform,"3d_b")==0){
    dxn = paramfile_get_int_param_default(s->ppl,"stim_frame_binoc_dxn",10);
    dxa = param_getf_dflt(s->ppl,"stim_frame_binoc_amp",0.5);
    // Sets s->d and s->d_r
    printf("    Getting binoc data...\n");
    mod_util_get_3d_stim_binoc(s,stim_num,xn,yn,tn,sscale,tscale,tsamp,mylogf);
    printf("    done.\n");

    if (fstflag == 1){
      fst_binoc = 2;
      data = f3tensor(1,xn,1,yn,1,2*tn); // For 3d fft
      for(i=1;i<=xn;i++)
	for(j=1;j<=yn;j++)
	  for(k=1;k<=tn;k++){
	    data[i][j][2*k-1] = s->d[i][j][k];
	    data[i][j][2*k]   = s->d_r[i][j][k];
	  }
      fst_tn = 2*tn;
    }else{
      data = f3tensor(1,2*xn+dxn,1,yn,1,tn); // For 3d fft
      for(i=1;i<=xn;i++)
	for(j=1;j<=yn;j++)
	  for(k=1;k<=tn;k++){
	    data[i       ][j][k] = s->d[i][j][k];
	    data[i+xn+dxn][j][k] = s->d_r[i][j][k];
	  }

      // Set 'dxn' columns of gray between images
      for(j=1;j<=yn;j++){
	for(k=1;k<=tn;k++){
	  for(i=1;i<=dxn;i++){
	    data[xn+i][j][k] = dxa;
	  }
	}
      }
      xn = 2*xn + dxn;
    }

  }else if (strcmp(sform,"3c")==0){
    mod_util_get_3c_stim(s,stim_num,xn,yn,tn,sscale,tscale,tsamp,mylogf);
    idata = s->dc;
    data = NULL;
  }else if (strcmp(sform,"3rgb")==0){
    mod_util_get_3rgb_stim(s,stim_num,xn,yn,tn,sscale,tscale,tsamp,mylogf);

    //exit_error("WYETH HERE"," Have to handle this RGB format somehow????");
    //data = s->drgb_l[0];  // WYETH HACK Send the R (or RGB) for now.

    data = data_util_3rgb_to_3d(s->drgb_l,xn,yn,tn);
    tn *= 3;  // The time dimension is now 3 times as long
    //
    //  *** WYETH - currently, this write the R,G,B frames separately in 
    //        in sequence.  We should make up a new tcode to have this get
    //        put back as a real color movie.
    //

    idata = NULL;
  }else{
    printf("  *** stim_form  %s\n",sform);
    exit_error("WM_WRITE_3D_STIM","Unknown 'stim_form'");
  }

  if (frame_num < 0){
    if (strcmp(sform,"3c")==0){
      write_3d_data_part(outfile,idata,0,xn,0,yn,0,tn,4,11,1); // 1=txy format
    }else if (strcmp(sform,"3rgb")==0){
      // Note, 'data' is starts at "0" (not 1)
      //write_3d_data_part(outfile,data,0,xn,0,yn,0,tn,4,2,1);
      write_3d_data_part(outfile,data,0,xn,0,yn,0,tn,4,3,1);
    }else{ // assume "3d"
      if (fstflag == 1){
	//
	// .fst format for web
	//
	printf(" *********** WYETH ASSUMING 0..1 value range MUST CHECK\n");
	data_util_fst_write(outfile,data,1,xn,1,yn,1,fst_tn,0,2,fst_binoc);
	//write_3d_data_part(outfile,data,1,xn,1,yn,1,tn,4,2,1); //1=txy format
      }else
	//
	// Standard 3D write
	//
	write_3d_data_part(outfile,data,1,xn,1,yn,1,tn,4,2,1); // 1=txy format
    }
  }else{
    if (frame_num > (tn-1)){
      printf("  stim_frame_tn  %d\n",tn);
      printf("  stim_frame     %d\n",frame_num);
      exit_error("WM_WRITE_3D_STIM","Frame number too large");
    }

    t = get_2d_farray(xn,yn);
    for(i=0;i<xn;i++)
      for(j=0;j<yn;j++)
	t[i][j] = data[1+i][1+j][1+frame_num];

    printf("  Will write frame %d to 2d data file.\n",frame_num);
    write_2d_data(outfile,t,0,0,xn,yn,4,2,0,1); // No transp for consistency
    free_2d_farray(t,xn);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           WM_WRITE_STIM_PPM_FILES                         */
/*                                                                           */
/*  Write a bunch of PPM files that can be converted to GIF.                 */
/*                                                                           */
/*  NOTES:                                                                   */
/*  use the script 'ppm2gif'.  See mail from Stuart Gilson.                  */
/*  - Assumes image values range from 0 to 1.                                */
/*                                                                           */
/*****************************************************************************/
void wm_write_stim_ppm_files(stimfile,outfile,arglist,narg,frame_num,stim_num)
     char stimfile[],outfile[];
     char **arglist;
     int narg,frame_num;
     int stim_num;
{
  int xn,yn,tn,tsamp;
  float sscale,tscale,stim_samp,***data;
  struct stim_struct *s; // Stimulus params

  printf("  WM_WRITE_STIM_PPM_FILES\n");

  //  Initialize stimulus struct
  s = (struct stim_struct *)myalloc(sizeof(struct stim_struct));
  s->ppl = paramfile_get_initial_params_update_only(stimfile,arglist,narg,1);
  strcpy(s->paramfile,stimfile);

  xn = paramfile_get_int_param_or_exit(s->ppl,"stim_frame_xn");
  yn = paramfile_get_int_param_or_exit(s->ppl,"stim_frame_yn");
  tn = paramfile_get_int_param_or_exit(s->ppl,"stim_frame_tn");
  sscale = paramfile_get_float_param_or_exit(s->ppl,"stim_frame_sscale");
  tscale = paramfile_get_float_param_or_exit(s->ppl,"stim_frame_tscale");
  stim_samp = paramfile_get_float_param_or_exit(s->ppl,"stim_samp");
  tsamp = my_rint(1.0/(tscale*stim_samp)); // Time units per stimulus frame
  if (tsamp <= 0)
    tsamp = 1; // Minimum value for tsamp

  // WYETH - I added this (quickly) Dec 20, 2018
  wm_prep_var_const_gen(NULL,s,mylogf,NULL,0,NULL); // Will set model seed 1777

  mod_util_get_3d_stim(s,stim_num,xn,yn,tn,sscale,tscale,tsamp,mylogf);

  write_ppm_6_3d_gray_frames(outfile,s->d,xn,yn,tn,1); // 1=dt, >1 skip frms

  // FREE data
}
/**************************************-**************************************/
/*                                                                           */
/*                              RUN_INDEFINITE                               */
/*                                                                           */
/*  Keep running the same model, appending returned data as next trial.      */
/*                                                                           */
/*****************************************************************************/
void run_indefinite(m,s,r)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response params
{
  int n,nchan;
  long nd_ntr_addr;
  struct ndtrial_struct *t;
  struct ndata_struct *nd;
  char *outfile,*outfile2,*prefix,*format;

  printf("  RUN_INDEFINITE\n");

  // Get 'outfile' name (contains .nd).  Optional 'outfile2' would be .t1
  wm_get_nd_t1_outfile(r,mylogf,-1,&prefix,&format,&outfile,&outfile2);

  //printf("prefix = %s\n",prefix);
  //printf("format = %s\n",format);
  //printf("outfile = %s\n",outfile);
  //printf("outfile2 = %s\n",outfile2);


  // ***** WYETH TO DO ***********
  // ***** WYETH TO DO ***********
  //  (2) check how the record 'tn' is set, is it in the 'r' struct???
  //  (3) freeing old trial storage???
  // ************************************************
  // ************************************************

  //
  //  Prepare ndata file that will be appended
  //
  nd = get_ndata();                  // Create an ndata structure
  nd->class = strdup("WM");          // Set the class
  write_ndata(outfile,nd);   // Write the file with no trials
  nd_ntr_addr = ndata_get_file_pos_ntr(outfile);  // Offset to overwrite ntr

  t = (struct ndtrial_struct *)myalloc(sizeof(struct ndtrial_struct));

  nchan = r->ns + r->nf;

  n = 0;  // Number of trials completed
  while(m->run_mode == 8){

    wm_slave_action(m,s,r,1,mylogf); // THIS MODIFIES s->d, s->d_r

    wm_response_check_trial_full(r,0,mylogf,NULL);  // Re-use 0th trial

    wm_make_ndtrial(m,s,r,mylogf,-1,0,nchan,nd,t);  // 0-resp trial

    //ndata_trial_print_trial(t,2);  // WYETH DEBUG

    ndata_append_trial(outfile,t);
    n += 1;
    ndata_overwrite_ntr(outfile,nd_ntr_addr,n);
    printf("    Appended %s, trials %d\n",outfile,n);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            RUN_GEN_TUNING_CURVE                           */
/*                                                                           */
/*  General tuning curve.                                                    */
/*                                                                           */
/*****************************************************************************/
void run_gen_tuning_curve(m,s,r)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response params
{
  int i,j;
  int xn,yn,tn,nstim,nrpt,write_stim_num; //OLD REMOVE: ,gti;
  char *sform;

  printf("  RUN_GEN_TUNING_CURVE\n");

  wm_get_sform_xn_yn_tn(m,s,&sform,&xn,&yn,&tn);

  write_stim_num = param_geti_dflt(s->ppl,"write_stim_num",-1);
  nrpt           = param_geti_exit(s->ppl,"stim_nrpt");

  nstim = s->ntr / nrpt;

  for(i=0;i<nstim;i++){

    //      i = stim number (among unique stimuli)
    // r->tri = trial index number
    get_stimulus_data_tuning_curve(m,s,r,i,r->tri,mylogf);
    // Stim data now in 's->d', or 's->d1' for 1d data
    //   Or, s->d is null if stim_type is "null"

    if (write_stim_num == i){
      char tname[SLEN];

      if ((strcmp(sform,"3d")==0) && (s->d != NULL)){
	sprintf(tname,"zzz.stim_%d.3d",i);
	write_3d_data_part(tname,s->d,1,xn,1,yn,1,tn,4,2,1);
      }else if ((strcmp(sform,"3d_b")==0) && (s->d != NULL)){
	sprintf(tname,"zzz.stim_%d.3d",i);
	write_3d_data_part(tname,s->d,1,xn,1,yn,1,tn,4,2,1);
	sprintf(tname,"zzz.stim_%d_R.3d",i);
	write_3d_data_part(tname,s->d_r,1,xn,1,yn,1,tn,4,2,1);
      }else if (strcmp(sform,"1d")==0){
	sprintf(tname,"zzz.stim_%d.pl",i);
	append_farray_plot(tname,tname,s->d1,tn,1);
      }
    }

    if (r->tname != NULL)
      myfree(r->tname);
    r->tname = wm_make_name_for_trial(s,r);

    for(j=0;j<nrpt;j++){
      s->repno = j; // Store repeat number
      printf("  _________________________\n");
      if (r->tname != NULL)
	printf("  Stim %d (%s) Rpt %d\n",i+1,r->tname,j+1);
      else
	printf("  Stim %d  Rpt %d\n",i+1,j+1);

      // Pointer for response storage, factoring in moo var
      r->gtsi = m->m_i * s->ntr + r->tsi; // global trial index

      wm_slave_action(m,s,r,1,mylogf); // THIS MODIFIES s->d, s->d_r


      //wm_response_check_trial_full(r,r->tri,mylogf,NULL);
      wm_response_check_trial_full(r,r->gtsi,mylogf,NULL);

      //r->dflag[r->tri] = 1; // Results done for this trial
      r->dflag[r->gtsi] = 1; // Results done for this trial
      //s->tribyx[r->tsi] = r->tri;  // Save 'tri' to match storage to vvals
      s->tribyx[r->gtsi] = r->tri;  // Save 'tri' to match storage to vvals
      r->tsi += 1;          // Trial stimulus index
      r->tri = r->tsi;      // Lock 'tri' and 'tsi' together
    }

    if (strcmp(sform,"3d")==0){
      if (s->d != NULL)
	free_f3tensor(s->d,1,xn,1,yn,1,tn);
    }else if (strcmp(sform,"3d_b")==0){
      if (s->d != NULL)
	free_f3tensor(s->d,1,xn,1,yn,1,tn);
      if (s->d_r != NULL)
	free_f3tensor(s->d_r,1,xn,1,yn,1,tn);
    }else if (strcmp(sform,"3c")==0){
      if (s->dc != NULL)
	free_3d_iarray(s->dc,xn,yn,tn);
    }else if (strcmp(sform,"3rgb")==0){
      if (s->drgb_l != NULL)
	free_4d_farray(s->drgb_l,3,xn,yn,tn);
    }else if (strcmp(sform,"1d")==0){
      myfree(s->d1);
    }else{
      exit_error("WM_RUN_GEN_TUNING_CURVE","Unknown stimulus form");
    }
  }

  if (r->tname != NULL)  myfree(r->tname);
  if (sform    != NULL)  myfree(sform);
}
/**************************************-**************************************/
/*                                                                           */
/*                         RUN_GEN_TUNING_CURVE_RAND                         */
/*                                                                           */
/*  Like 'run_gen_tuning_curve' but with blockwise randomisation and some    */
/*  features omitted:                                                        */
/*  - No writing of stimuli                                                  */
/*  - No 'fseq' condition.                                                   */
/*                                                                           */
/*****************************************************************************/
void run_gen_tuning_curve_rand(m,s,r,shuff_seed)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response params
     int shuff_seed;             // Shuffle seed, for blockwise randomisation
{
  int i,j;
  int xn,yn,tn,nstim,nrpt,si,**randi;
  char *sform;

  printf("  RUN_GEN_TUNING_CURVE_RAND\n");

  wm_get_sform_xn_yn_tn(m,s,&sform,&xn,&yn,&tn);

  nrpt = paramfile_get_int_param_or_exit(s->ppl,"stim_nrpt");
  nstim = s->ntr / nrpt;

  randi = (int **)myalloc(nrpt*sizeof(int *));
  for(i=0;i<nrpt;i++)
    randi[i] = get_shuffle_index_return_seed(nstim,3,&shuff_seed);

  for(j=0;j<nrpt;j++){
    for(i=0;i<nstim;i++){
    
      si = randi[j][i];     // si = stimulus ID number, among unique stimuli
      r->tri = si*nrpt + j; // Index in ordered list, e.g. for 'vval'
      get_stimulus_data_tuning_curve(m,s,r,si,r->tri,mylogf);
      // Stim data now in 's->d', or 's->d1' for 1d data
      
      if (r->tname != NULL)
	myfree(r->tname);

      r->tname = wm_make_name_for_trial(s,r);


      s->repno = j; // Store repeat number
      if (r->tname != NULL)
	printf("  Stim %d (%s) Rpt %d\n",si+1,r->tname,j+1);
      else
	printf("  Stim %d  Rpt %d\n",si+1,j+1);
      wm_slave_action(m,s,r,1,mylogf); // THIS MODIFIES s->d, s->d_r

      wm_response_check_trial_full(r,r->tri,mylogf,NULL);

      r->dflag[r->tri] = 1;  // Results done for this trial
      s->tribyx[r->tsi] = r->tri;  // Save 'tri' to match storage to vvals
      r->tsi += 1;           // Trial Sequential Index, for response storage
      
      if (strcmp(sform,"3d")==0){
	free_f3tensor(s->d,1,xn,1,yn,1,tn);
      }else if (strcmp(sform,"3c")==0){
	free_3d_iarray(s->dc,xn,yn,tn);
	//free_f3tensor(s->dc,1,xn,1,yn,1,tn);
      }else if (strcmp(sform,"3d_b")==0){
	free_f3tensor(s->d,1,xn,1,yn,1,tn);
	free_f3tensor(s->d_r,1,xn,1,yn,1,tn);
      }else if (strcmp(sform,"1d")==0){
	myfree(s->d1);
      }
    }
  }
  free_2d_iarray(randi,nrpt);
}
/**************************************-**************************************/
/*                                                                           */
/*                              RUN_INTERACTIVE                              */
/*                                                                           */
/*****************************************************************************/
void run_interactive(m,s,r)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response params
{
  int k;
  int xn,yn,tn,done,told,sold,nrpt;
  char *sform;

  printf("  RUN_INTERACTIVE\n");

  wm_get_sform_xn_yn_tn(m,s,&sform,&xn,&yn,&tn);

  nrpt = paramfile_get_int_param_or_exit(s->ppl,"stim_nrpt");
  
  done = 0;
  while(!done){
    // Note, s->stimno gets set in 'mod_gui_pop_main'
    if (s->stimno < 0){
      //exit_error("RUN_INTERACTIVE","stimno < 0");
      printf("  *** WARNING:  setting stimno to zero\n");
      s->stimno = 0;  // WYETH 2008 Sept, used to exit
    }

    r->tri = s->stimno * nrpt; // Use first trial num
    get_stimulus_data_tuning_curve(m,s,r,s->stimno,r->tri,mylogf);
    sold = s->stimno;
    // Stim data now in 's->d', 's->d_r', or 's->d1' for 1d data
    told = r->tri;
    s->repno = 0; // Set repeat number to 0 for interactive

    printf("    ");
    for(k=0;k<s->nvar;k++)
      printf("  %s=%s",s->vname[k],s->vval[r->tri][k]);
    printf("\n");

    while(sold == s->stimno){
      wm_slave_action(m,s,r,1,mylogf); // THIS MODIFIES s->d, s->d_r
    }

    if (strcmp(sform,"3d")==0){
      free_f3tensor(s->d,1,xn,1,yn,1,tn);
    }else if (strcmp(sform,"3c")==0){
      free_3d_iarray(s->dc,xn,yn,tn);
      //free_f3tensor(s->dc,1,xn,1,yn,1,tn);
    }else if (strcmp(sform,"3d_b")==0){
      free_f3tensor(s->d,1,xn,1,yn,1,tn);
      free_f3tensor(s->d_r,1,xn,1,yn,1,tn);
    }else if (strcmp(sform,"1d")==0){
      myfree(s->d1);
    }else{
      exit_error("RUN_INTERACTIVE","Bad sform");
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                MODEL_READ_M                               */
/*                                                                           */
/*****************************************************************************/
void model_read_m(paramfile,arglist,narg,rm)
     char paramfile[];
     char **arglist;
     int narg;
     struct model_struct **rm;    // Model params.
{
  int i;
  int *ulist;                // Flags for cmd line args used in moo file
  struct model_struct *m;    // Model params.

  // Initialize model struct
  m = (struct model_struct *)myalloc(sizeof(struct model_struct));
  m->marfile = NULL;
  m->action = -999;   // Arbitrary value
  m->calib_pval = NULL;
  m->mds = NULL;
  m->mds_n = 0;

  if (is_suffix_string(paramfile,".mod")){
    m->ppl = paramfile_get_initial_params_update_only(paramfile,arglist,narg,
						      1);
    m->o = NULL;
  }else{
    m->ppl = NULL;
    m->o = paramfile_onode_file_read(paramfile,1);
    paramfile_onode_resolve(m->o,mylogf); // Resolve 'COPY' assignments
    ulist = onode_update_from_slist(m->o,mylogf,arglist,narg);
    // 'ulist' not used here.
    myfree(ulist);
  }

  strcpy(m->paramfile,paramfile);
  m->process_id = -1; // Single process
  m->run_mode = 0;    // Normal run mode, can change to interactive

  // Check if all command line params. occur in at least one param file
  if (m->ppl != NULL){
    for(i=0;i<narg;i+=2)
      if (!paramfile_test_param(m->ppl,arglist[i])){
	printf("  *** Parameter name:  %s\n",arglist[i]);
	exit_error("MODEL_READ_M","Command line param not found");
      }
  }else{
    for(i=0;i<narg;i+=2)
      if (!onode_test_ostr(m->o,arglist[i])){
	printf("  *** Parameter name:  %s\n",arglist[i]);
	exit_error("MODEL_READ_M","Command line param not found");
      }
  }
  *rm = m;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MODEL_READ_DESCRIPTION                          */
/*                                                                           */
/*****************************************************************************/
void model_read_description(paramfile,stimfile,respfile,arglist,narg,rm,rs,rr,
			    rulist)
     char paramfile[],stimfile[],respfile[];
     char **arglist;
     int narg;
     struct model_struct **rm;    // Model params.
     struct stim_struct **rs;     // Stimulus params.
     struct response_struct **rr; // Response params.
     int **rulist;   // [narg] flag array, 1-this cmd line par was moo.
{
  int i;
  int *ulist;                // Flags for cmd line args used in moo file
  struct model_struct *m;    // Model params.
  struct stim_struct *s;     // Stimulus params.
  struct response_struct *r; // Response params.

  // Initialize stimulus struct
  //
  //  WYETH - could make stim struct init routine in e.g. 'wm_util' to set
  //          NULL val
  //
  //        - Use '.../model/wm/dev/a/null.stm' 
  //
  if (strcmp(stimfile,"NULL")==0){
    s = NULL;
  }else{
    s = (struct stim_struct *)myalloc(sizeof(struct stim_struct));
    s->ppl = paramfile_get_initial_params_update_only(stimfile,arglist,narg,1);
    strcpy(s->paramfile,stimfile);
    s->name = (char *)NULL;
    s->stimno = -1;
  }
  //paramfile_print_ppl(s->ppl);   // WYETH DEBUG
  //exit(0);


  // Initialize model struct
  m = (struct model_struct *)myalloc(sizeof(struct model_struct));
  m->marfile = NULL;
  m->action = -999;   // Arbitrary value
  m->calib_pval = NULL;
  m->mds = NULL;
  m->mds_n = 0;

  if (is_suffix_string(paramfile,".mod")){
    m->ppl = paramfile_get_initial_params_update_only(paramfile,arglist,narg,
						      1);
    m->o = NULL;
    ulist = get_zero_iarray(narg); // WYETH Sept 2010
  }else{
    m->ppl = NULL;
    m->o = paramfile_onode_file_read(paramfile,1);
    paramfile_onode_resolve(m->o,mylogf); // Resolve 'COPY' assignments
    ulist = onode_update_from_slist(m->o,mylogf,arglist,narg);
  }

  strcpy(m->paramfile,paramfile);
  m->process_id = -1; // Single process
  m->run_mode = 0;    // Normal run mode, can change to interactive
  m->nrun = 1;        // Assume one model, change below if moo var pars
  m->vps = NULL;

  // Initialize response struct
  r = (struct response_struct *)myalloc(sizeof(struct response_struct));
  r->ppl = paramfile_get_initial_params_update_only(respfile,arglist,narg,1);
  strcpy(r->paramfile,respfile);
  r->ns = -1;      // Number of spike responses per trial
  r->nf = -1;      // Number of continuous responses per trial
  r->tname = NULL; // trial name, constructed from var params

  // Check if all command line params. occur in at least one param file
  if (m->ppl != NULL){
    for(i=0;i<narg;i+=2)
      if (!(paramfile_test_param(m->ppl,arglist[i])||
	    paramfile_test_param(s->ppl,arglist[i])||
	    paramfile_test_param(r->ppl,arglist[i]))){
	printf("  *** Parameter name:  %s\n",arglist[i]);
	exit_error("MODEL_READ_DESCRIPTION","Command line param not found");
      }
  }else{
    if (s != NULL){
      for(i=0;i<narg;i+=2)
	if (!(onode_test_ostr(m->o,arglist[i])||
	      paramfile_test_param(s->ppl,arglist[i])||
	      paramfile_test_param(r->ppl,arglist[i]))){
	  printf("  *** Parameter name:  %s\n",arglist[i]);
	  exit_error("MODEL_READ_DESCRIPTION","Command line param not found");
	}
    }else{
      for(i=0;i<narg;i+=2)
	if (!(onode_test_ostr(m->o,arglist[i])||
	      paramfile_test_param(r->ppl,arglist[i]))){
	  printf("  *** Parameter name:  %s\n",arglist[i]);
	  exit_error("MODEL_READ_DESCRIPTION","Command line param not found");
	}
    }
  }

  *rm = m; *rs = s; *rr = r;
  *rulist = ulist;
}
/**************************************-**************************************/
/*                                                                           */
/*                                 MODEL_MAIN                                */
/*                                                                           */
/*****************************************************************************/
void model_main(paramfile,stimfile,respfile,arglist,narg)
     char paramfile[],stimfile[],respfile[];
     char **arglist;
     int narg;
{
  int i,k;
  int *parflag;               // flag of which cmd line arg used
  int calgen;                 // 1-calibration generate mode, 0-otherwise
  int sseed;                  // Shuffle seed
  struct model_struct *m;     // Model params
  struct stim_struct *s;      // Stimulus params
  struct response_struct *r;  // Response params

  /*
  int debug_wait;
  debug_wait = 1;
  while(debug_wait) ;
  printf("HERE HERE__________\n");
  */

  model_read_description(paramfile,stimfile,respfile,arglist,narg,&m,&s,&r,
			 &parflag);

  // 'parflag' uses 1's to indicate which params were used by the .moo file
  wm_prep_var_const_gen(m,s,mylogf,arglist,narg,parflag);

  wm_response_init(m,s,r,mylogf);

  // Must occur BEFORE 'wm_moo_var_par_set' and AFTER 'wm_response_init'
  calgen = wmu_calibration_init(m,r,mylogf);  // Return 1 if "generate" mode

  sseed = param_geti_dflt(s->ppl,"stim_shuff_seed",0);

  if (m->vps != NULL){ // If MOO params vary, set up the first values
    wm_moo_var_par_set(mylogf,m,0);
  }

  // Here, stim can be changed, e.g. "stim_auto_center_mi" in mod_mesh_util
  m->m_i = 0;
  wm_slave_action(m,s,r,0,mylogf);   // Prep

  if (m->run_mode == 0){
    
    for(i=0;i<m->nrun;i++){ // For each model config

      if (i > 0){  // Any model var params were already set for i=0
	printf("  ________________________________________________________\n");
	printf("\n");
	wm_moo_var_par_set(mylogf,m,i); // Update the params of the model
	m->m_i = i;
	wm_slave_action(m,s,r,0,mylogf);   // Re-prep

	r->tri = 0;  // Stimulus trial pointers
	r->tsi = 0;
	r->tname = NULL;
      }

      if (sseed == 0)
	run_gen_tuning_curve(m,s,r);             // Sequential stimuli
      else
	run_gen_tuning_curve_rand(m,s,r,sseed);  // Blockwise randomise

      wm_slave_action(m,s,r,-1,mylogf);   // Clean up
    }

    wm_write_response(m,s,r,mylogf,-1); // Write response

    if (calgen == 1) //  Write a calibration file *** CALIBRATE 2018 Nov
      wmu_calibration_write(m,s,r,mylogf);

  }else if (m->run_mode == 8){
    run_indefinite(m,s,r);        // Run until run_mode is changed to '9'
  }else
    run_interactive(m,s,r);             // Run interactive
}
/**************************************-**************************************/
/*                                                                           */
/*                                 MODEL_MAR                                 */
/*                                                                           */
/*****************************************************************************/
void model_mar(paramfile,marfile,arglist,narg)
     char paramfile[];
     char marfile[];
     char **arglist;
     int narg;
{
  struct model_struct *m;     // Model params
  struct stim_struct *s;      // Stimulus params
  struct response_struct *r;  // Response params

  s = NULL;
  r = NULL;

  model_read_m(paramfile,arglist,narg,&m);

  m->marfile = marfile;

  wm_slave_action(m,s,r,10,mylogf);   // Prep model and write .mar file
}
/**************************************-**************************************/
/*                                                                           */
/*                                 MODEL_MDS                                 */
/*                                                                           */
/*  Write out a list of the model data available for response request.       */
/*                                                                           */
/*****************************************************************************/
void model_mds(paramfile,arglist,narg)
     char paramfile[];
     char **arglist;
     int narg;
{
  struct model_struct *m;     // Model params
  struct stim_struct *s;      // Stimulus params
  struct response_struct *r;  // Response params

  s = NULL;
  r = NULL;

  model_read_m(paramfile,arglist,narg,&m);

  wm_slave_action(m,s,r,20,mylogf);   // Prep model and write .mar file
}
/**************************************-**************************************/
/*                                                                           */
/*                               MODEL_RSP_GEN                               */
/*                                                                           */
/*  Response file generation.                                                */
/*                                                                           */
/*****************************************************************************/
void model_rsp_gen(paramfile,respfile,outfile,arglist,narg)
     char paramfile[];
     char respfile[];    // Commands to generate an response file
     char outfile[];     // output file - resulting .rsp file
     char **arglist;
     int narg;
{
  struct model_struct *m;     // Model params
  struct stim_struct *s;      // Stimulus params
  struct response_struct *r;  // Response params

  s = NULL;

  model_read_m(paramfile,arglist,narg,&m);

  // Initialize response struct, just to store 'respfile' name
  r = (struct response_struct *)myalloc(sizeof(struct response_struct));
  strcpy(r->paramfile,respfile);
  r->ppl = NULL;
  r->tname = NULL;
  r->ns = r->nf = -1;

  m->marfile = outfile; // WYRESP

  wm_slave_action(m,s,r,20,mylogf);   // Prep model and generate .rsp file
}
/**************************************-**************************************/
/*                                                                           */
/*                                    MAIN                                   */
/*                                                                           */
/*****************************************************************************/
int main(argc,argv)
     int argc;
     char *argv[];
{
  int fnum,snum;
  char infile[SLEN],outfile[SLEN],option[SLEN],marfile[SLEN];
  char paramfile[SLEN],paramfile2[SLEN],stimfile[SLEN],respfile[SLEN];

  if (argc == 1){
    printf("\n  wm <option>\n");
    printf("    where <option> is\n");
    printf("      mod <model_file> <stimfile> <respfile>");
    printf(" [<pname> <pval> ...]\n");
    printf("      modm <model_file> <stimfile> <respfile>");
    printf(" [<pname> <pval> ...]\n");
    printf("      diff <model_file_1> <model_file_2>\n");
    printf("      write_3d_stim <stim_file> <outfile> <stim_num> ");
    printf("[<pname> <pval> ...]\n");
    printf("      write_3d_stim_ppm <stimfile> <outfile> <stim_num>");
    printf("[<pname> <pval>..]\n");
    printf("      write_fst <stim_file> <outfile> <stim_num> ");
    printf("[<pname> <pval> ...]\n");
    printf("      write_2d_stim_frame <stim_file> <outfile> <stim_num> ");
    printf(" <frame_num> [<pname> <pval> ...]\n");
    printf("      conn2txt <infile> <outfile>\n");
    printf("      mar <model_file> <outfile>");
    printf(" [<pname> <pval> ...]\n");
    printf("      rsp_gen <model_file> <rsp_gen_file> <outfile>");
    printf(" [<pname> <pval> ...]\n");
    printf("      list_resp_data <model_file>");
    printf(" [<pname> <pval> ...]\n");
    printf("\n");

    printf("  Size of pop_cell struct  %ld bytes\n",sizeof(struct pop_cell));
    printf("  Size of pop_syn  struct  %ld bytes\n",sizeof(struct pop_syn));
    printf("\n");

    printf("  Version 3.0\n");
    printf("\n");

    exit(0);
  }
  printf("\nWM.C\n");

  strcpy(option,argv[1]);

  if (strcmp(option,"mod")==0){
    if (argc < 5)
      exit_error("MAIN","Too few parameters for option 'mod'");
    strcpy(paramfile,argv[2]);
    strcpy(stimfile,argv[3]);
    strcpy(respfile,argv[4]);
    model_main(paramfile,stimfile,respfile,&(argv[5]),argc-5);
  }else if (strcmp(option,"mar")==0){
    strcpy(paramfile,argv[2]);
    strcpy(marfile,argv[3]);
    model_mar(paramfile,marfile,&(argv[4]),argc-4);
  }else if (strcmp(option,"rsp_gen")==0){
    strcpy(paramfile,argv[2]);
    strcpy(respfile,argv[3]);
    strcpy(marfile,argv[4]);
    model_rsp_gen(paramfile,respfile,marfile,&(argv[4]),argc-5);
  }else if (strcmp(option,"modm")==0){
    printf("<br>\n<pre>\n");
    mctrl_util_determine_dir_structure();
    mctrl_util_report_pid_run();
    strcpy(paramfile,argv[2]);
    strcpy(stimfile,argv[3]);
    strcpy(respfile,argv[4]);
    model_main(paramfile,stimfile,respfile,&(argv[5]),argc-5);

    // *** WYETH - added Sept 5th.
    // *** WYETH ADDED FOR TEMP CHANGE TO KAMBLIPOOCHI
    append_string_to_file("JOB_DONE","Done\n");

  }else if (strcmp(option,"diff")==0){
    strcpy(paramfile,argv[2]);
    strcpy(paramfile2,argv[3]);
    paramfile_read_diff(paramfile,paramfile2);
  }else if (strcmp(option,"write_2d_stim_frame")==0){
    strcpy(infile,argv[2]);
    strcpy(outfile,argv[3]);
    snum = atoi(argv[4]);
    fnum = atoi(argv[5]);
    wm_write_3d_stim(infile,outfile,&(argv[6]),argc-6,snum,fnum,0);
  }else if (strcmp(option,"write_3d_stim")==0){
    strcpy(infile,argv[2]);
    strcpy(outfile,argv[3]);
    snum = atoi(argv[4]);
    wm_write_3d_stim(infile,outfile,&(argv[5]),argc-5,snum,-1,0);
  }else if (strcmp(option,"write_fst")==0){
    strcpy(infile,argv[2]);
    strcpy(outfile,argv[3]);
    snum = atoi(argv[4]);
    wm_write_3d_stim(infile,outfile,&(argv[5]),argc-5,snum,-1,1);
  }else if (strcmp(option,"write_3d_stim_ppm")==0){
    // WYETH - I added (quickly) the snum param here, Dec 20, 2018
    strcpy(infile,argv[2]);
    strcpy(outfile,argv[3]);
    snum = atoi(argv[4]);
    wm_write_stim_ppm_files(infile,outfile,&(argv[5]),argc-5,-1,snum);
  }else if (strcmp(option,"conn2txt")==0){
    strcpy(infile,argv[2]);
    strcpy(outfile,argv[3]);
    mod_conn_conn_to_text(mylogf,infile,outfile);
  }else if (strcmp(option,"list_resp_data")==0){
    strcpy(paramfile,argv[2]);
    model_mds(paramfile,&(argv[3]),argc-3);
  }else
    printf("  *** Unknown option.\n");

  printf("\n");
  return 0;
}
