/*****************************************************************************/
/*                                                                           */
/*   mod_dog_util.c                                                          */
/*   wyeth bair                                                              */
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
#include "carray_util.h"
#include "iarray_util.h"
#include "farray_util.h"
#include "plot_util.h"
#include "func_util.h"
#include "myrand_util.h"
#include "data_util.h"
#include "paramfile_util.h"
#include "fft_util.h"
#include "kernel_util.h"
#include "ifc_util.h"
#include "mod_mesh_util.h"
#include "mod_util.h"
#include "pop_util.h"
#include "mod.h" // Data structures
#include "paramfile.h"

// For cluster computation  (static needed w/ initialization
static int myid = -1;         // Process ID, -1-Single, 0-Master, >0-Slave
static char *mylogf = NULL;   // Logfile name
char ggstr[LONG_SLEN];        // Global log string

struct mod_dog_filt{ // Added Jan 10, 2009
  char *name;        // name of pop/layer
  float sig1;        // center SD
  float sig2;        // surround SD
  float amp1;        // center amplitude
  float amp2;        // surround amplitude
  float cs_x;        // Allow center to be offset w.r.t. surround
  float cs_y;
  float cs_delay;    // surround delay w.r.t. center
  float tsig;
  float tcent;
  float tdelay;
  float tau_ad;

  float power;

  float ***f;        // DOG model filter
  float  **s;        // for DOG model filter FT
  float ***r;        // response, kept for efficient repeats
  float ***rr;       // response for right eye, binocular stimuli

  float ***rg;       // response gain, e.g., r^2, kept for efficient repeats
  float ***rrg;      // response gain, right eye

  int oldstimno;     // Old stim number, to avoid convolving same stim twice

  float fsum;        // Sum of filter, used for OFF response, (psum + nsum)
  float psum;        // Sum of negative values in filter; negative (or 0)
  float nsum;        // Sum of positive values in filter; positive (or 0)
};

int mod_dog_fsn;                   // Number of filter structs
struct mod_dog_filt **mod_dog_fs;  // pointers to filter structs [...fsn]

int     mod_dog_stim_si = -1;      // For index assoc'd w/ ..._stim_s
float **mod_dog_stim_s  = NULL;    // For FT of stimulus

float ***mod_dog_surr   = NULL;    // Filter for non-classical surround
float  **mod_dog_surr_s = NULL;    // For FT of filter
float ***mod_dog_sup    = NULL;    // Suppressive signal
float ***mod_dog_sup_r  = NULL;    // Suppressive signal, right eye, binoc

float mod_dog_sscale;
float mod_dog_tscale;
int   mod_dog_xn;
int   mod_dog_yn;
int   mod_dog_tn;

float mod_dog_surr_sig;    // For extra-classical surround
float mod_dog_surr_tau;    // For extra-classical surround

int      mod_dog_simp_n;    // Number of simple cells
int     *mod_dog_simp_nin;  // [n] Number of inputs per cell
int    **mod_dog_simp_x;    // [n][nin] x-coord of DOG input
int    **mod_dog_simp_y;    // [n][nin] y-coord of DOG input
float  **mod_dog_simp_w;    // [n][nin] weight of DOG input
float ***mod_dog_simp_rf;   // [n][xn][yn] rf profile
int   ***mod_dog_simp_flag; // [n][xn][yn] set to 1 if input

int    mod_dog_cg_flag;     // Contrast Gain

int    mod_dog_ds_nsub;     // Number of ds subunits
int   *mod_dog_ds_subx1;    // x coord of subunit member 1
int   *mod_dog_ds_suby1;    // y coord of subunit member 1
int   *mod_dog_ds_subx2;    // x coord of subunit member 2
int   *mod_dog_ds_suby2;    // y coord of subunit member 2
float *mod_dog_ds_subd;     // delay

// Squash control globals, Blake
int mod_dog_sq_flag;        // Type of squash, 0 for none
int mod_dog_sq_dump;        // Flag to dump out intermediate values

// SQUASH TABLE
int mod_dog_squash_flag = 0;  // Use squash table
int    sqtab_n;           // Length of table
float  sqtab_max;         // maximum output value, asymptote of squashing func
float  sqtab_off_pre;     // subtract before squashing, add back after
float  sqtab_dx;          // minimum value, less than this y = x
float *sqtab = NULL;      // y-values [sqtab_n]

// Gain control globals, Blake
int     mod_dog_nl_flag;      // Type of nonlinearity, 0 for none
int     mod_dog_nl_num;       // Number of filters
int     mod_dog_nl_length;    // Length of filter
int     mod_dog_nl_sigdump;   // Flag to dump out intermediate values
int     mod_dog_nl_fildump;   // Flag to dump out filter values
float   mod_dog_nl_gstep;     // The spacing of the gain for the filters
float   mod_dog_nl_tres;      // The spacing of the gain for the filters
float   mod_dog_nl_mag;       // Controls the total amplitude of the filter
float   mod_dog_nl_off;       // Controls the offset of the filter
float **mod_dog_nl_filt = NULL;    // Filter data [nfilt][filt_len]
float   mod_dog_nl_tau;       // Time const for gain signal integration


char *mod_dog_logistic_name;     // "logistic", or NULL if non
float mod_dog_logistic_slope;
float mod_dog_logistic_asymp; 

/**************************************-**************************************/
/*                                                                           */
/*                             MOD_DOG_FS_BY_NAME                            */
/*                                                                           */
/*  Return a pointer to the filter structures for the given pop 'name'.      */
/*                                                                           */
/*****************************************************************************/
struct mod_dog_filt *mod_dog_fs_by_name(name)
     char *name;              // population name, or NULL to use fs0
{
  int i;
  int n;
  struct mod_dog_filt *mdfs;

  mdfs = NULL;
  n = 0;
  for(i=0;i<mod_dog_fsn;i++){
    if (strcmp(name,mod_dog_fs[i]->name)==0){
      mdfs = mod_dog_fs[i];
      n += 1;
    }
  }

  if (n == 0){
    sprintf(ggstr,"  name:  %s\n",name);
    mylog(mylogf,ggstr);
    mylog_exit(mylogf,"MOD_DOG_FS_BY_NAME  population name not found.\n");
  }else if (n > 1){
    sprintf(ggstr,"  found %d times, name:  %s\n",n,name);
    mylog(mylogf,ggstr);
    mylog_exit(mylogf,"MOD_DOG_FS_BY_NAME  population name not unique.\n");
  }

  return mdfs;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MODEL_DOG_GET_SAMPLING_POINTS                     */
/*                                                                           */
/*  Pick integer sampling points under a Gaussian RF profile of SD 'ssd'.    */
/*  Set values of RF profile to 0.0 that are less than 'eps'.                */
/*  Each point chosen contributes a Gaussian weight of SD 'wsd' to a weight  */
/*  mask, which is never allowed to exceed 'rfs' times the RF profile.       */
/*                                                                           */
/*****************************************************************************/
void model_dog_get_sampling_points(xn,yn,n,seed,ssd,rfs,wsd,eps,rsx,rsy,rn)
     int xn,yn,n,seed;
     float ssd,rfs,wsd,eps;
     int **rsx,**rsy,*rn;
{
  int i,j,k;
  int xx,yy,xi,yi,wxn,wyn,wxc,wyc,wflag,*sx,*sy,nfail,failmax;
  float **rf,min,max,xc,yc;
  float **wm,**wcurr;
  char tstr[SLEN];

  mylog(mylogf,"  MODEL_DOG_GET_SAMPLING_POINTS\n");
  failmax = 400;

  /*** Compute a weight mask ***/
  wxn = wyn = 1 + 2.0*my_rint(wsd * 3.0);
  sprintf(ggstr,"    Weight mask is %d x %d\n",wxn,wyn);
  mylog(mylogf,ggstr);
  wxc = (wxn-1)/2;
  wyc = (wyn-1)/2;
  wm = gaussian_2d(wxn,wyn,(float)wxc,(float)wyc,wsd,1.0);
  wcurr = get_zero_2d_farray(xn,yn);
  write_2d_data("zz.rf.wmask.2d",wm,0,0,wxn,wyn,4,2,1,0);

  if (seed > 0)
    seed = -seed;
  
  /*** Make RF profile, where peak value is 1.0 ***/
  xc = (float)(xn - 1)/2.0;
  yc = (float)(yn - 1)/2.0;
  rf = gaussian_2d(xn,yn,xc,yc,ssd,1.0);
  make_max_const_2d_farray(rf,xn,yn,1.0);

  /*** Constrain RF to be zero where it is less than 'eps' ***/
  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++)
      if (rf[i][j] < eps)
	rf[i][j] = 0.0;

  max = max_of_2d_farray(rf,xn,yn);
  min = min_of_2d_farray(rf,xn,yn);
  sprintf(ggstr,"    RF min, max =  %f %f\n",min,max);
  mylog(mylogf,ggstr);

  write_2d_data("zz.rf.2d",rf,0,0,xn,yn,4,2,1,0);

  /*** Pick 'n' points in the Gaussian RF profile ***/
  sx = (int *)myalloc(n*sizeof(int));
  sy = (int *)myalloc(n*sizeof(int));

  nfail = 0;
  k = 0;
  while((k<n) && (nfail < failmax)){
    xi = (int)((float)xn * myrand_util_ran2(&seed));
    yi = (int)((float)yn * myrand_util_ran2(&seed));
    if (myrand_util_ran2(&seed) < rf[xi][yi]){
      wflag = 1; /* It is OK to add this sampling point */
      for(i=0;i<wxn;i++){
	xx = xi - wxc + i;
	if ((xx >= 0)&&(xx < xn)){
	  for(j=0;j<wyn;j++){
	    yy = yi - wyc + j;
	    if ((yy >= 0)&&(yy < yn)){
	      if (rf[xx][yy] > 0.0){
		if ((wcurr[xx][yy] + wm[i][j]) > (rfs * rf[xx][yy])){
		  wflag = 0;
		}
	      }
	    }
	  }
	}
      } 
      if (wflag == 1){ /*** OK to add this point ***/
	for(i=0;i<wxn;i++){
	  xx = xi - wxc + i;
	  if ((xx >= 0)&&(xx < xn)){
	    for(j=0;j<wyn;j++){
	      yy = yi - wyc + j;
	      if ((yy >= 0)&&(yy < yn))
		wcurr[xx][yy] += wm[i][j];
	    }
	  }
	}
	sx[k] = xi;
	sy[k] = yi;
	k += 1;
	nfail = 0;
      }else{
	nfail += 1;
      }
    }
  }

  if (nfail < failmax){
    sprintf(ggstr,
	    "    Chose %d sampling points (failure limit not exceeded)\n",k);
    mylog(mylogf,ggstr);
  }else{
    sprintf(ggstr,
	    "    Chose %d sampling points (stopped after %d failures)\n",k,
	    nfail);
    mylog(mylogf,ggstr);
  }

  write_2d_data("zz.rf.weight.2d",wcurr,0,0,xn,yn,4,2,1,0);

  if (myid == -1){
    remove_file("zz.rf.points.pl");
    for(i=0;i<k;i++){
      sprintf(tstr,"%d %d\n",sx[i],sy[i]);
      append_string_to_file("zz.rf.points.pl",tstr);
    }
    remove_file("zz.rf.slice.pl");
    for(i=0;i<yn;i++){
      sprintf(tstr,"%d",i);
      append_farray_plot("zz.rf.slice.pl",tstr,wcurr[i],xn,1);
    }
  }

  free_2d_farray(rf,xn); free_2d_farray(wm,wxn); free_2d_farray(wcurr,xn);

  *rsx = sx; *rsy = sy; *rn = k;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_DOG_SQUASH_TABLE_PREP                        */
/*                                                                           */
/*  Make a look-up table to apply a squashing function to the filter         */
/*  output.                                                                  */
/*                                                                           */
/*****************************************************************************/
void mod_dog_squash_table_prep(m)
     struct model_struct *m; // Model params
{
  int i;
  int dump_flag;
  float x,*sqtabx;
  struct onode *o;

  mylog(mylogf,"  MOD_DOG_SQUASH_TABLE_PREP\n");

  if (m->o != NULL)
    o = onode_get_node_type_item_val(m->o,"pop","name","lgn_on");
  else
    o = NULL;

  if (o == NULL){
    //mylog_exit(mylogf,"MOD_DOG_SQUASH_TABLE_PREP  pop lgn not found.\n");
    sqtab_max = paramfile_get_float_param_or_exit(m->ppl,"squash_out_max");
    sqtab_off_pre = paramfile_get_float_param_default(m->ppl,
						      "squash_offset_pre",0.0);
  }else{
    sqtab_max = onode_getpar_flt_exit(o,"squash_out_max");
    sqtab_off_pre = onode_getpar_flt_dflt(o,"squash_offset_pre",0.0);
  }
  
  sqtab_n = 100;         // Length of table
  sqtab_dx = 4.0* sqtab_max / (float)sqtab_n;      // increment along x-axis

  if (sqtab != NULL)
    myfree(sqtab);
  sqtab = get_farray(sqtab_n);

  sqtabx = get_farray(sqtab_n);

  x = 0.0;
  for(i=0;i<sqtab_n;i++){
    // TANH has slope 1/2, so mult by 2.
    sqtab[i] = sqtab_max * func_sigmoid_tanh(x*2.0/sqtab_max,1.0,0.0);
    sqtabx[i] = x;
    x += sqtab_dx;
  }
  
  if (myid == -1){
    if (o == NULL)
      dump_flag = paramfile_get_float_param_default(m->ppl,"squash_dump",0);
    else
      dump_flag = onode_getpar_int_dflt(o,"squash_dump",0);
    if (dump_flag == 1)
      append_farray_xy_plot("zz.sqtab.pl",sqtabx,sqtab,sqtab_n,"Squash_Table");
  }
  myfree(sqtabx);
}
/**************************************-**************************************/
/*                                                                           */
/*                                  Blake                                    */
/*                            MOD_DOG_NONLIN1_PREP                           */
/*                                                                           */
/*  Set up a family of filters that will be used depending on the current    */
/*  level of the gain control signal.                                        */
/*                                                                           */
/*****************************************************************************/
void mod_dog_nonlin1_prep(m)
     struct model_struct *m; /* Model params. */
{
  int i,j; // a counter
  float gain,posamp,negamp,negadj,posadj,minratio; // for controlling filter
  float alpha,beta,a,scale,magadj,mag,negscale,posscale; // for filter
  struct onode *o; // the object node
  
  // this is a debugging statement
  sprintf(ggstr,"  MOD_DOG_NONLIN1_PREP\n");
  mylog(mylogf,ggstr);
  
  // check file type
  if (m->ppl != NULL){
    
    // .mod file
    
    // get all of the parameter values
    mod_dog_nl_num = paramfile_get_int_param_default(m->ppl,"nl_num",20);
    mod_dog_nl_length = paramfile_get_int_param_default(m->ppl,"nl_length",40);
    mod_dog_nl_sigdump= paramfile_get_int_param_default(m->ppl,"nl_sigdump",0);
    mod_dog_nl_fildump= paramfile_get_int_param_default(m->ppl,"nl_fildump",0);
    mod_dog_nl_gstep = paramfile_get_float_param_default(m->ppl,"nl_gstep",
							 0.1);
    mod_dog_nl_tres = paramfile_get_float_param_default(m->ppl,"nl_tres",
							 0.1);
    mod_dog_nl_mag = paramfile_get_float_param_default(m->ppl,"nl_mag",
							0.001);
    mod_dog_nl_off = paramfile_get_float_param_default(m->ppl,"nl_offset",
							0.0);
    mod_dog_nl_tau = paramfile_get_float_param_default(m->ppl,"nl_tau",20);
    posamp = paramfile_get_float_param_default(m->ppl,"nl_posamp",1.0);
    negamp = paramfile_get_float_param_default(m->ppl,"nl_negamp",1.0);
    minratio = paramfile_get_float_param_default(m->ppl,"nl_minratio",1.0);
    magadj = paramfile_get_float_param_default(m->ppl,"nl_magadj",1.0);
    alpha = paramfile_get_float_param_default(m->ppl,"nl_alpha",2.0);
    beta = paramfile_get_float_param_default(m->ppl,"nl_beta",0.5);
  }else{
    
    // .moo file
    
    // get the object node pointer for lgn_on, note this means only lgn_on
    // needs to have these parameteres specified
    o = onode_get_node_type_item_val(m->o,"pop","name","lgn_on");
    
    // get all of the parameter values
    mod_dog_nl_num = onode_getpar_int_dflt(o,"nl_num",20);
    mod_dog_nl_length = onode_getpar_int_dflt(o,"nl_length",40);
    mod_dog_nl_sigdump = onode_getpar_int_dflt(o,"nl_sigdump",0);
    mod_dog_nl_fildump = onode_getpar_int_dflt(o,"nl_fildump",0);
    mod_dog_nl_gstep = onode_getpar_flt_dflt(o,"nl_gstep",0.1);
    mod_dog_nl_tres = onode_getpar_flt_dflt(o,"nl_tres",0.1);
    mod_dog_nl_mag = onode_getpar_flt_dflt(o,"nl_mag",0.001);
    mod_dog_nl_off = onode_getpar_flt_dflt(o,"nl_offset",0.0);
    mod_dog_nl_tau = onode_getpar_flt_dflt(o,"nl_tau",20);
    posamp = onode_getpar_flt_dflt(o,"nl_posamp",1.0);
    negamp = onode_getpar_flt_dflt(o,"nl_negamp",1.0);
    minratio = onode_getpar_flt_dflt(o,"nl_minratio",1.0);
    magadj = onode_getpar_flt_dflt(o,"nl_magadj",1.0);
    alpha = onode_getpar_flt_dflt(o,"nl_alpha",2.0);
    beta = onode_getpar_flt_dflt(o,"nl_beta",0.5);
  }

  // print some info
  sprintf(ggstr,"   num %i, length %i, gstep %f, gmod %f, tau %f, offset %f\n",
	  mod_dog_nl_num,mod_dog_nl_length,mod_dog_nl_gstep,mod_dog_nl_mag,
	  mod_dog_nl_tau,mod_dog_nl_off);
  mylog(mylogf,ggstr);
  sprintf(ggstr,"    posamp: %f, negamp: %f, alpha: %f, beta: %f\n",posamp,
	  negamp,alpha,beta);
  mylog(mylogf,ggstr);
  
  // create all of the mod_dog_nl_filt functions
  if (mod_dog_nl_filt != NULL){ /* Assuming that sizes haven't changed */
    free_2d_farray(mod_dog_nl_filt,mod_dog_nl_num);
  }

  mod_dog_nl_filt = (float **)myalloc(mod_dog_nl_num*sizeof(float *));
  for(i = 0; i < mod_dog_nl_num; i++){
    
    gain = mod_dog_nl_gstep*i; // set shape controlling variable, a, with gain
    a = alpha ;//+ gain;
    
    // scale the amplitudes
    if (mod_dog_nl_num == 1) {
      negscale = 1;
      posscale = 1;
    } else {
      negscale = (minratio +(1-minratio)*((float)i/(float)(mod_dog_nl_num-1)));
      // EDIT negscale = (minratio + 
      //   (1-minratio)*((float)(mod_dog_nl_num-i)/(float)(mod_dog_nl_num-1)));
      posscale = (magadj + (1-magadj)*((float)i/(float)(mod_dog_nl_num-1)));
      // EDIT posscale = (magadj + 
      // (1-magadj)*((float)(mod_dog_nl_num-i)/(float)(mod_dog_nl_num-1)));
    }
    negadj = negscale*posscale*negamp;
    posadj = posscale*posamp;
    
    // build the filter
    mod_dog_nl_filt[i] = 
      (float *)diff_inv_gamma_farray(mod_dog_tscale,a,beta,posadj,
				     negadj,mod_dog_nl_length);
    
    // if the dump flag is set, store the filter
    if(mod_dog_nl_fildump){
      sprintf(ggstr,"NL_Filter_%i",i);
      append_farray_plot("zz.nl_fildump.pl",ggstr,mod_dog_nl_filt[i],
			 mod_dog_nl_length,1);
    }
    
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                  Blake                                    */
/*                              MOD_DOG_NONLIN1                              */
/*                                                                           */
/*  Apply the non-linearity to the floating signal 't'.                      */
/*                                                                           */
/*****************************************************************************/
void mod_dog_nonlin1(m,t,tn)
     struct model_struct *m; /* Model params. */
     float *t;               /* Response data [tn] */
     int tn;                 /* Length of response data */
{
  int i,j, filtnum, inum;
  int *indices, *ival, *icount;
  float sum, meannl, sdevnl, meanl, sdevl, meang, sdevg;
  float *tempt,*gainsums;
  
  // A degubbing printout
  sprintf(ggstr,"    MOD_DOG_NONLIN1\n");
  mylog(mylogf,ggstr);
  
  // make an array to store t temporarily for calculations
  // and an array to store which filters to use
  tempt = (float *)myalloc(tn*sizeof(float));
  gainsums = (float *)myalloc(tn*sizeof(float));
  indices = (int *)myalloc(tn*sizeof(int));
  ival = (int *)myalloc(mod_dog_nl_num*sizeof(int));
  icount = (int *)myalloc(mod_dog_nl_num*sizeof(int));
  
  // figure out which filters to use
  find_nl_bin_index(indices,mod_dog_nl_num,tn,t,mod_dog_nl_tres,
		    mod_dog_nl_tau,gainsums);
  
  // report the mean filter index
  get_unique_count_iarray(indices,tn,&ival,&icount,&inum);
  sprintf(ggstr,"      Index count: ");
  mylog(mylogf,ggstr);
  for(i=0;i<inum;i++) {
    sprintf(ggstr,"%i:%i ",ival[i],icount[i]);
    mylog(mylogf,ggstr);
  }
  sprintf(ggstr,"\n");
  mylog(mylogf,ggstr);
  
  // if the dump flag is set, store the input to the 
  // non-linearity and the gainsums
  if(mod_dog_nl_sigdump){
    append_farray_plot("zz.nl_sigdump.pl","Linear_Drive",t,tn,1);
    append_farray_plot("zz.nl_sigdump.pl","Gain_sums",gainsums,tn,1);
  }
  
  // report the linear signal
  mean_sdev_farray(t,tn,&meanl,&sdevl);
  sprintf(ggstr,"      Linear signal: m: %f, sd: %f\n",meanl,sdevl);
  mylog(mylogf,ggstr);
  
  // report the gain signal
  mean_sdev_farray(gainsums,tn,&meang,&sdevg);
  sprintf(ggstr,"      Gain signal: m: %f, sd: %f\n",meang,sdevg);
  mylog(mylogf,ggstr);
  
  // modify t using one of the temporal filters
  for(i=0; i<tn; i++){
    
    // store the initial t value (mean centered)
    tempt[i] = t[i];
    
    // integrate over the past and modify t
    sum = 0.0;
    for(j=0; j < mod_dog_nl_length && j<=i; j++){

      sum += tempt[i-j]*mod_dog_nl_mag*mod_dog_nl_filt[indices[i]][j];
    }
    t[i] = sum + mod_dog_nl_off;
  }
  
  // report the mean non-linear results
  mean_sdev_farray(t,tn,&meannl,&sdevnl);
  sprintf(ggstr,"      Non-Linear signal: m: %f, sd: %f\n",meannl,sdevnl);
  mylog(mylogf,ggstr);
  
  // if the dump flag is set, store the output from the non-linearity
  if(mod_dog_nl_sigdump){
    append_farray_plot("zz.nl_sigdump.pl","NonLinear_Output",t,tn,1);
  }
  
  // free the temporary storage arrays
  myfree(tempt);
  myfree(gainsums);
  myfree(indices);
  myfree(ival);
  myfree(icount);
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_DOG_SQUASH_BLAKE                           */
/*                                                                           */
/*  Contributed by Blake Richards for his MSc project.                       */
/*                                                                           */
/*  A function to send a signal through a chosen squash function.            */
/*                                                                           */
/*****************************************************************************/
void mod_dog_squash_blake(m,t,tn,pflag)
     struct model_struct *m; /* Model params. */
     float *t;               /* Response data [tn] */
     int tn;                 /* Length of response data */
     int pflag; // 1-print
{
  int i;
  int squashflag; // a flag for the chosen function
  float severity,midpoint,scale,shift; // severity, scale, shift and midpoint
  float inm, insdev, sqm, sqsdev; // vars for reporting stats
  struct onode *o; // the object node
  
  // check file type
  if (m->ppl != NULL){
    
    // .mod file
    mod_dog_sq_dump = paramfile_get_int_param_default(m->ppl,"sq_dump",0);
    squashflag = paramfile_get_int_param_default(m->ppl,"sq_type",0);
    midpoint = paramfile_get_float_param_default(m->ppl,"sq_midpoint",0.0);
    severity = paramfile_get_float_param_default(m->ppl,"sq_severity",1.0);
    scale = paramfile_get_float_param_default(m->ppl,"sq_scale",1.0);
    shift = paramfile_get_float_param_default(m->ppl,"sq_shift",0.0);
  }else{
    
    // .moo file
    o = onode_get_node_type_item_val(m->o,"pop","name","lgn_on");
    mod_dog_sq_dump = onode_getpar_int_dflt(o,"sq_dump",0);
    squashflag = onode_getpar_int_dflt(o,"sq_type",0);
    midpoint = onode_getpar_flt_dflt(o,"sq_midpoint",0.0);
    severity = onode_getpar_flt_dflt(o,"sq_severity",1.0);
    scale = onode_getpar_flt_dflt(o,"sq_scale",1.0);
    shift = onode_getpar_flt_dflt(o,"sq_shift",0.0);
  }

  // A printout
  if (pflag){
    sprintf(ggstr,"    MOD_DOG_SQUASH\n");
    mylog(mylogf,ggstr);
    sprintf(ggstr,"      Midpoint: %f, severity: %f, scale: %f, shift: %f\n",
	    midpoint,severity,scale,shift);
    mylog(mylogf,ggstr);
  }
  
  // dump info if requested
  if(mod_dog_sq_dump)
    append_farray_plot("zz.sqdump.pl","Input_to_Squash",t,tn,1);
  
  // report the mean and sdev of the input
  mean_sdev_farray(t,tn,&inm,&insdev);
  if (pflag){
    sprintf(ggstr,"      Input: m: %f, sd: %f\n",inm,insdev);
    mylog(mylogf,ggstr);
  }
  
  // step through each time
  for(i=0; i<tn; i++){
    
    if(squashflag == 0) // 0 == tanh
      t[i] = scale*(midpoint - (float)func_sigmoid_tanh((t[i]+shift),
							severity,0.0));
    else if (squashflag == 1) // 1 == logistic
      t[i] = (float)func_sigmoid_log(t[i],scale,midpoint,severity,shift); 
  }
  
  // dump info if requested
  if(mod_dog_sq_dump)
    append_farray_plot("zz.sqdump.pl","Output_from_Squash",t,tn,1);
  
  // report the mean and sdev of the input
  mean_sdev_farray(t,tn,&sqm,&sqsdev);
  if (pflag){
    sprintf(ggstr,"      Squash: m: %f, sd: %f\n",sqm,sqsdev);
    mylog(mylogf,ggstr);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               MOD_DOG_SQUASH                              */
/*                                                                           */
/*  This version uses a squash table, and linear interpolation.              */
/*                                                                           */
/*****************************************************************************/
void mod_dog_squash(d,n)
     float *d;  // Data to squash
     int n;     // length of data
{
  int i,k;
  float x,dx;
  
  // step through each time
  for(i=0;i<n;i++){
    x = d[i] + sqtab_off_pre;
    if (x > 0.0){  // Only squashing positive values
      k = x / sqtab_dx;
      if (k >= (sqtab_n - 1))
	x = sqtab_max;  // Asymptotic value
      else{ // Linear interp.
	dx = (x - (float)k*sqtab_dx) / sqtab_dx;
	x = sqtab[k] + (sqtab[k+1] - sqtab[k]) * dx;
//printf(" k = %d   sqtab[k],[k+1] = %f %f   dx = %f   d[i]= %f  x %f\n",
//     k,sqtab[k],sqtab[k+1],dx,d[i],x);
      }
      d[i] = x - sqtab_off_pre;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MODEL_DOG_01_PREP                            */
/*                                                                           */
/*  1.  Set parameters for the DOG filter.                                   */
/*  2.  Create the filter at 'mod_dog_f'.                                    */
/*                                                                           */
/*****************************************************************************/
void model_dog_01_prep(m,r)
     struct model_struct *m;    // Model params
     struct response_struct *r; // Response params
{
  int xn,yn,tn,ix,iy,iz,pflag;
  float ***dog,tscale,sscale,surr_sig,surr_tau;
  char *outfile,*prefix;
  struct param_pair_list *mppl;
  struct mod_dog_filt *f;

  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  mylog(mylogf,"  MODEL_DOG_01_PREP\n");

  pflag = 0;
  if (myid == -1)
    pflag = 1;

  // Get outfile prefix
  prefix = paramfile_get_char_param_default(r->ppl,"outfile","zz");

  mppl = m->ppl;
  sscale = paramfile_get_float_param_or_exit(mppl,"sscale");
  tscale = paramfile_get_float_param_or_exit(mppl,"tscale");
  xn = paramfile_get_int_param_or_exit(mppl,"xn");
  yn = paramfile_get_int_param_or_exit(mppl,"yn");
  tn = paramfile_get_int_param_or_exit(mppl,"tn");

  mod_dog_sscale = sscale;
  mod_dog_tscale = tscale;
  mod_dog_xn = xn;
  mod_dog_yn = yn;
  mod_dog_tn = tn;

  f = (struct mod_dog_filt *)myalloc(sizeof(struct mod_dog_filt));

  f->sig1 = paramfile_get_float_param_or_exit(mppl,"sig1");
  f->sig2 = paramfile_get_float_param_or_exit(mppl,"sig2");
  f->amp1 = paramfile_get_float_param_or_exit(mppl,"amp1");
  f->amp2 = paramfile_get_float_param_or_exit(mppl,"amp2");
  f->cs_x = paramfile_get_float_param_default(mppl,"cs_x",0.0);
  f->cs_y = paramfile_get_float_param_default(mppl,"cs_y",0.0);

  // WYETH - OLD WAY, new uses "cs_delay_s" in *SECONDS*
  f->cs_delay = paramfile_get_float_param_or_exit(mppl,"cs_delay");

  f->tsig = paramfile_get_float_param_or_exit(mppl,"tsig");
  f->tcent = paramfile_get_float_param_default(mppl,"tcent",0.0);
  f->tdelay = paramfile_get_float_param_default(mppl,"tdelay",0.0);
  f->tau_ad = paramfile_get_float_param_default(mppl,"tau_ad",0.0);

  f->power = paramfile_get_float_param_default(mppl,"power",1.0);
  

  mod_dog_nl_flag = paramfile_get_int_param_default(mppl,"nl_flag",0);
  mod_dog_sq_flag = paramfile_get_int_param_default(mppl,"sq_flag",0);
  mod_dog_cg_flag = paramfile_get_int_param_default(mppl,"contrast_gain",0);


  mylog(mylogf,"    Computing 3D DOG filter.\n");

  outfile = paramfile_make_fname(mppl,"write_dog_temporal",prefix,".dogt.pl");

  f->f = dog_space_time_tensor(xn,yn,tn,sscale,tscale,f->sig1,f->sig2,f->amp1,
			       f->amp2,f->cs_x,f->cs_y,f->cs_delay,f->tsig,
			       f->tcent,outfile,mppl,m->o,pflag,0,0,0.0,
			       f->tau_ad);
  
  kernel_util_norm_xyt_filter(f->f,1,xn,1,yn,1,tn,&(f->psum),&(f->nsum));
  printf("  Filter 1, sum pos,neg  %f %f\n",f->psum,f->nsum);

  dog = f->f;

  mod_dog_fsn = 1;
  mod_dog_fs = (struct mod_dog_filt **)myalloc(mod_dog_fsn*
					       sizeof(struct mod_dog_filt *));
  mod_dog_fs[0] = f;


  if (outfile != (char *)NULL)
    myfree(outfile);

  f->fsum = sum_3d_farray(dog,1,xn,1,yn,1,tn);
  sprintf(ggstr,"    sum over 3d filter = %f",f->fsum);
  mylog(mylogf,ggstr);
  get_max_coord_3d_farray(dog,1,xn,1,yn,1,tn,&ix,&iy,&iz);
  sprintf(ggstr,"    Max at (%d,%d,%d) in [1,%d][1,%d][1,%d]\n",
	  ix,iy,iz,xn,yn,tn);
  mylog(mylogf,ggstr);

  // For writing a 1D plot through the filter
  outfile = paramfile_make_fname(mppl,"write_trace_filter",prefix,".trf.pl");
  if (outfile != (char *)NULL){
    int wx,wy,wz,tdn;
    float *td;

    wx = paramfile_get_nth_int_param_or_exit(mppl,"write_trace_filter",1);
    wy = paramfile_get_nth_int_param_or_exit(mppl,"write_trace_filter",2);
    wz = paramfile_get_nth_int_param_or_exit(mppl,"write_trace_filter",3);

    if (wx == -1){
      td = get_1d_from_3d_farray(dog,1,xn,wy,1,wz,1);
      tdn = xn;
    }else if (wy == -1){
      td = get_1d_from_3d_farray(dog,wx,1,1,yn,wz,1);
      tdn = yn;
    }else if (wz == -1){
      td = get_1d_from_3d_farray(dog,wx,1,wy,1,1,tn);
      tdn = tn;
    }else{
      td = NULL;
      tdn = 0;
      mylog_exit(mylogf,"MODEL_DOG_01_PREP  Must set one index to -1\n");
    }
    write_farray_plot(outfile,td,tdn);
    myfree(td);
    myfree(outfile);
  }
//mod_dog_f = dog;

  // Make filter for surround
  surr_sig = paramfile_get_float_param_default(mppl,"surr_sd",0.0);
  if (surr_sig > 0.0){
    surr_tau = paramfile_get_float_param_or_exit(mppl,"surr_tau");
    mod_dog_surr = gauss_space_exp_time_tensor(mylogf,xn,yn,tn,sscale,tscale,
					       surr_sig,surr_tau);
    mod_dog_surr_sig = surr_sig;
    mod_dog_surr_tau = surr_tau;

    mylog(mylogf,"    Computing FFT of Surround filter.\n");
    mod_dog_surr_s = matrix(1,xn,1,2*yn);
    contort_real_3d_farray(mod_dog_surr,1,1,1,xn,yn,tn);
    rlft3(mod_dog_surr,mod_dog_surr_s,xn,yn,tn,1);  // In nr_util.c
  }

  outfile = paramfile_make_fname(mppl,"write_3d_filter",prefix,".dog.3d");
  if (outfile != (char *)NULL){
    write_3d_data_part(outfile,dog,1,xn,1,yn,1,tn,4,2,1);
    myfree(outfile);
  }

  // Check for squashing function table
  //   NOTE - 'lgn_on' params control LGN ON and OFF cells.
  mod_dog_squash_flag = paramfile_get_int_param_default(mppl,"squash_flag",0);
  if (mod_dog_squash_flag)
    mod_dog_squash_table_prep(m);

  myfree(prefix);
}
/**************************************-**************************************/
/*                                                                           */
/*                     MOD_DOG_FILT_STRUCT_FREE_CONTENTS                     */
/*                                                                           */
/*****************************************************************************/
void mod_dog_filt_struct_free_contents(f)
     struct mod_dog_filt *f;
{
  int xn,yn,tn;

  if (f->name != NULL)  myfree(f->name);

  xn = mod_dog_xn;
  yn = mod_dog_yn;
  tn = mod_dog_tn;

  if (f->f   != NULL)  free_f3tensor(f->f,1,xn,1,yn,1,tn);
  if (f->s   != NULL)  free_matrix(f->s,1,xn,1,2*yn);
  if (f->r   != NULL)  free_f3tensor(f->r,1,xn,1,yn,1,tn);
  if (f->rr  != NULL)  free_f3tensor(f->rr,1,xn,1,yn,1,tn);
  if (f->rg  != NULL)  free_3d_farray(f->rg,xn,yn,tn);
  if (f->rrg != NULL)  free_3d_farray(f->rrg,xn,yn,tn);

}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_DOG_GET_FILTER_STRUCT                        */
/*                                                                           */
/*****************************************************************************/
struct mod_dog_filt *mod_dog_get_filter_struct(po,pflag)
     struct onode *po;        // Population onode for an LGN pop
     int pflag;               // For printing output
{
  int xn,yn,tn,ix,iy,iz;
  float sscale,tscale;
  char *dumpfile;
  struct mod_dog_filt *f;

  f = (struct mod_dog_filt *)myalloc(sizeof(struct mod_dog_filt));

  f->name     = onode_getpar_chr_exit(po,"name");  // Free this storage
  f->sig1     = onode_getpar_flt_exit(po,"sig1");
  f->sig2     = onode_getpar_flt_exit(po,"sig2");
  f->amp1     = onode_getpar_flt_exit(po,"amp1");
  f->amp2     = onode_getpar_flt_exit(po,"amp2");
  f->cs_x     = onode_getpar_flt_dflt(po,"cs_x",0.0);
  f->cs_y     = onode_getpar_flt_dflt(po,"cs_y",0.0);
  if (onode_test_ostr(po,"cs_delay") == 1)
    //
    //  OLD WAY - do not use any more
    //
    f->cs_delay = onode_getpar_flt_exit(po,"cs_delay");
  else{
    f->cs_delay = onode_getpar_flt_exit(po,"cs_delay_s");
    f->cs_delay /= mod_dog_tscale;
  }
  f->tsig     = onode_getpar_flt_exit(po,"tsig");
  f->tcent    = onode_getpar_flt_dflt(po,"tcent",0.0);
  f->tdelay   = onode_getpar_flt_dflt(po,"tdelay",0.0);
  f->power    = onode_getpar_flt_dflt(po,"power",1.0);
  f->tau_ad   = onode_getpar_flt_dflt(po,"tau_ad",0.0);

  if (myid == -1){
    dumpfile = onode_getpar_chr_dflt(po,"write_dog_temporal",NULL);
  if (dumpfile != NULL){
    if (strcmp(dumpfile,"NULL")==0) // Interpret string "NULL" as NULL
      dumpfile = NULL;
  }
  }else
    dumpfile = NULL;

  // WYETH - THIS WAS an older and more general way to get dumpfile names:
  //outfile =paramfile_make_fname(mppl,"write_dog_temporal",prefix,".dogt.pl");

  xn = mod_dog_xn;
  yn = mod_dog_yn;
  tn = mod_dog_tn;
  sscale = mod_dog_sscale;
  tscale = mod_dog_tscale;


  //  
  //  Create filter
  //  
  f->f = dog_space_time_tensor(xn,yn,tn,sscale,tscale,f->sig1,f->sig2,
			       f->amp1,f->amp2,f->cs_x,f->cs_y,f->cs_delay,
			       f->tsig,f->tcent,
			       dumpfile,NULL,po,pflag,0,0,0.0,f->tau_ad);
  //printf("WYETH Debug:  getting DELTA instead of DOG\n");
  //f->f = delta_tensor(mylogf,xn,yn,tn,sscale,tscale);

  //  
  //  Normalize filter
  //  
  kernel_util_norm_xyt_filter(f->f,1,xn,1,yn,1,tn,&(f->psum),&(f->nsum));

  //  
  //  Compute sums
  //  
  f->fsum = sum_3d_farray(f->f,1,xn,1,yn,1,tn);
  sprintf(ggstr,"    sum over 3d filter = %f",f->fsum);
  mylog(mylogf,ggstr);
  get_max_coord_3d_farray(f->f,1,xn,1,yn,1,tn,&ix,&iy,&iz);
  sprintf(ggstr,"    Max at (%d,%d,%d) in [1,%d][1,%d][1,%d]\n",
	  ix,iy,iz,xn,yn,tn);
  mylog(mylogf,ggstr);

  f->s  = NULL;
  f->r  = NULL;
  f->rr = NULL;

  f->rg  = NULL; // Gain
  f->rrg = NULL;

  f->oldstimno = -1;

  if (dumpfile != NULL)
    myfree(dumpfile);

  return f;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MODEL_DOG_01_PREP_O                            */
/*                                                                           */
/*  1.  Set parameters for the DOG filter.                                   */
/*  2.  Create the filter at 'mod_dog_f'.                                    */
/*                                                                           */
/*****************************************************************************/
void model_dog_01_prep_o(m,r)
     struct model_struct *m;    // Model params
     struct response_struct *r; // Response params
{
  int k;
  int xn,yn,tn,ix,iy,iz,pflag;
  float tscale,sscale,surr_sig,surr_tau;
  char *ltype;
  struct onode *o,*t,*nlo;

  if (mylogf == NULL)  // Do not do this if this is a "re-prep"
    myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  mylog(mylogf,"  MODEL_DOG_01_PREP_O\n");

  mod_dog_fsn = 0;  // Assume for now, in case there are no LGN pops

  pflag = 0;
  if (myid == -1)
    pflag = 1;

  mod_dog_sscale = sscale = onode_getpar_flt_exit(m->o,"sscale");
  mod_dog_tscale = tscale = onode_getpar_flt_exit(m->o,"tscale");
  mod_dog_xn     = xn     = onode_getpar_int_exit(m->o,"xn");
  mod_dog_yn     = yn     = onode_getpar_int_exit(m->o,"yn");
  mod_dog_tn     = tn     = onode_getpar_int_exit(m->o,"tn");

  //
  //  WYETH - THIS CAN BE REMOVED, OLD FILES probably won't have 'type' vals
  //
  o = onode_get_node_type_item_val(m->o,"pop","name","lgn");
  if (o == NULL){
    o = onode_get_node_type_item_val(m->o,"pop","name","lgn_on");
    if (o == NULL){
      mylog(mylogf,"  *** Found neither pop 'lgn' nor 'lgn_on'\n");
      return;  // WYETH - now allowing there to be no LGN
      //mylog_exit(mylogf,"MODEL_DOG_01_PREP_O  pop 'lgn' not found.\n");
    }else{
      exit_error("MODEL_DOG_01_PREP_O","Old lgn pop format, please update");
    }
  }

  
  //  Count the number of LGN populations
  mod_dog_fsn = 0;
  t = onode_get_next_type(m->o->o,"pop");
  while(t != NULL){
    ltype = onode_getpar_chr_ptr(t,"type");
    if (ltype != NULL){
      if (pop_util_check_lclass(ltype,"lgn0"))
	mod_dog_fsn += 1;
    }
    t = onode_get_next_type(t->next,"pop");
  }
  mod_dog_fs = (struct mod_dog_filt **)myalloc(mod_dog_fsn*
					       sizeof(struct mod_dog_filt *));
  sprintf(ggstr,"    %d LGN populations\n",mod_dog_fsn);
  mylog(mylogf,ggstr);
  if (mod_dog_fsn == 0)
    mylog_exit("MODEL_DOG_01_PREP_O  No LGN populations found.\n");

  
  //  Create filter structures for each LGN population
  k = 0;
  t = onode_get_next_type(m->o->o,"pop");
  while(t != NULL){
    ltype = onode_getpar_chr_ptr(t,"type");
    if (ltype != NULL){
      if (pop_util_check_lclass(ltype,"lgn0")){
	mod_dog_fs[k] = mod_dog_get_filter_struct(t,pflag);
	k += 1;
      }
    }
    t = onode_get_next_type(t->next,"pop");
  }


  mod_dog_nl_flag = onode_getpar_int_dflt(o,"nl_flag",0);
  mod_dog_sq_flag = onode_getpar_int_dflt(o,"sq_flag",0);
  mod_dog_cg_flag = onode_getpar_int_dflt(o,"contrast_gain",0);


  // Make filter for surround
  surr_sig = onode_getpar_flt_dflt(o,"surr_sd",0.0);
  if (surr_sig > 0.0){
    //exit_error("(mod_dog_util)  WYETH REMOVE THIS EXIT","Testing");
    surr_tau = onode_getpar_flt_exit(o,"surr_tau");
    mod_dog_surr = gauss_space_exp_time_tensor(mylogf,xn,yn,tn,sscale,tscale,
					       surr_sig,surr_tau);
    mod_dog_surr_sig = surr_sig;
    mod_dog_surr_tau = surr_tau;

    mylog(mylogf,"    Computing FFT of Surround filter.\n");
    mod_dog_surr_s = matrix(1,xn,1,2*yn);
    contort_real_3d_farray(mod_dog_surr,1,1,1,xn,yn,tn);
    rlft3(mod_dog_surr,mod_dog_surr_s,xn,yn,tn,1);  // In nr_util.c
  }

  // Check for squashing function table
  //   NOTE - 'lgn_on' params control LGN ON and OFF cells.
  mod_dog_squash_flag = onode_getpar_int_dflt(o,"squash_flag",0);
  if (mod_dog_squash_flag){
    //exit_error("(mod_dog_util)  WYETH REMOVE THIS EXIT","Testing");
    mod_dog_squash_table_prep(m);
  }

  nlo = onode_child_get_unique(o,"nonlin");
  if (nlo != NULL){
    //exit_error("(mod_dog_util)  WYETH REMOVE THIS EXIT","Testing");
    mod_dog_logistic_name  = onode_getpar_chr_exit(nlo,"name");
    mod_dog_logistic_slope = onode_getpar_flt_exit(nlo,"slope");
    mod_dog_logistic_asymp = onode_getpar_flt_exit(nlo,"asymptote");

    mylog(mylogf,"  LOGISTIC\n");
    sprintf(ggstr,"    name  %s\n",mod_dog_logistic_name);
    mylog(mylogf,ggstr);
    sprintf(ggstr,"    slope  %f\n",mod_dog_logistic_slope);
    mylog(mylogf,ggstr);
    sprintf(ggstr,"    asymp  %f\n",mod_dog_logistic_asymp);
    mylog(mylogf,ggstr);
  }else{
    mod_dog_logistic_name = NULL;
    mod_dog_logistic_slope = -1.0;
    mod_dog_logistic_asymp = 0.0;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_DOG_CG01_RESPONSE                         */
/*                                                                           */
/*  Implement a form of contrast gain control.                               */
/*                                                                           */
/*****************************************************************************/
void mod_dog_cg01_response()
{
  int i,j,k;
  int xn,yn,tn;
  int nm;
  float ***r,*t,*lowpass,*mask,*v,mean,var,*lowpv,*rnorm,*rnew;

  mylog(mylogf,"******* MOD_DOG_CG01_RESPONSE\n");
  mylog(mylogf,"******* WYETH - this attempt should be done w/ quad filts?\n");

  xn = mod_dog_xn;
  yn = mod_dog_yn;
  tn = mod_dog_tn;
  r = mod_dog_fs[0]->r;

  mean = mod_dog_fs[0]->fsum/2.0;
  
  mylog_exit(mylogf,"WYETH ---- THIS psum was calc'd *BEFORE* norm in past\n");
  var  = mod_dog_fs[0]->psum - mean;

  nm = tn/2;
  if (nm > 256)
    nm = 256;
  //mask = maxwell_farray(0.0,20.0,1.0,nm);
  mask = alpha_farray(0.0,0.05,1.0,nm);
  norm_area_farray(mask,nm,1.0);
  
  if (myid == -1)
    append_farray_plot("zz.cg.pl","mask",mask,nm,1);
  
  /*
    for(i=1;i<=xn;i++){
    for(j=1;j<=yn;j++){*/
  
  i = j = xn/2;
  t = r[i][j];
  
  rnorm = (float *)myalloc(tn*sizeof(float));
  v = (float *)myalloc(tn*sizeof(float));

  for(k=0;k<tn;k++){
    rnorm[k] = (t[k+1] - mean);

    if (rnorm[k] < 0.0)                 /* Use ABS instead of SQR */
      v[k] = -rnorm[k]/var;
    else
      v[k] = rnorm[k]/var;
  }

  lowpass = convolve_with_mask_causal(rnorm,tn,mask,nm);

  lowpv = convolve_with_mask_causal(v,tn,mask,nm);

  rnew = (float *)myalloc(tn*sizeof(float));
  for(k=0;k<tn;k++){
    /*** WYETH - instead of lowpv, use spatial integral of v ***/
    rnew[k] = rnorm[k] - 2.0*lowpv[k] * lowpass[k];
    rnew[k] += mean; /* Add mean back */
  }


  if (myid == -1){
    append_farray_plot("zz.cg.pl","v",v,tn,1);
    append_farray_plot("zz.cg.pl","lowpv",lowpv,tn,1);
    append_farray_plot("zz.cg.pl","rnew",rnew,tn,1);
    append_farray_plot("zz.cg.pl","lowpass",lowpass,tn,1);
    append_farray_plot("zz.cg.pl","rnorm",rnorm,tn,1);
    append_farray_plot("zz.cg.pl","r",t+1,tn,1);
  }

  exit(0);
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_DOG_SURR_RESPONSE                          */
/*                                                                           */
/*  Implement a form of surround suppression.                                */
/*                                                                           */
/*****************************************************************************/
void mod_dog_surr_response(eye_flag)
     int eye_flag;  // 0-left, 1-right
{
  int i,j,k;
  int xn,yn,tn,xi,yi;
  float ***r,***rsq,*t,***surr,***f,**f_s,mu,y;

  mylog(mylogf,"  MOD_DOG_SURR_RESPONSE\n");

  xn = mod_dog_xn;
  yn = mod_dog_yn;
  tn = mod_dog_tn;
  if (eye_flag == 0)
    r = mod_dog_fs[0]->r;
  else
    r = mod_dog_fs[0]->rr;
  f = mod_dog_surr;
  f_s = mod_dog_surr_s;

  xi = 1 + xn/2;
  yi = 1 + yn/2;
  
  // Make a new array that is the square the response array
  rsq = f3tensor(1,xn,1,yn,1,tn);
  for(i=1;i<=xn;i++){
    for(j=1;j<=yn;j++){
      t = r[i][j];

      mu = t[1];  // Start running mean at first value
      for(k=1;k<=tn;k++){
	y = t[k] - mu;      // Subtract the running mean
	rsq[i][j][k] = y*y; // Square the zero-mean signal
	mu = 0.98*mu + 0.02*t[k];  // Update running mean
      }
    }
  }

  /***
  if (myid == -1){
    append_farray_plot("zzz.mut.pl","mu",mm,tn,1);
    append_farray_plot("zzz.mut.pl","t",yy,tn,1);
    }
  myfree(mm);
  myfree(yy);***/


  // Convolve the array with a Gaussian in space, decay exp in time
  compute_response_single_from_ft(rsq,f,f_s,xn,yn,tn,0,&surr);

  // Keep the suppressive signal
  if (eye_flag == 0){
    if (mod_dog_sup != NULL)
      free_f3tensor(mod_dog_sup,1,xn,1,yn,1,tn);
    mod_dog_sup = surr;
  }else{
    if (mod_dog_sup_r != NULL)
      free_f3tensor(mod_dog_sup_r,1,xn,1,yn,1,tn);
    mod_dog_sup_r = surr;
  }

  free_f3tensor(rsq,1,xn,1,yn,1,tn);
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_DOG_GAIN_PTR                             */
/*                                                                           */
/*  Return a pointer to the gain response array.  0 indexing (not tensor).   */
/*                                                                           */
/*****************************************************************************/
float ***mod_dog_gain_ptr(name,eyeflag)
     char *name;
     int eyeflag; // 0-left 1-right
{
  struct mod_dog_filt *mdfs;
  
  if (name != NULL)
    mdfs = mod_dog_fs_by_name(name);
  else{
    mylog(mylogf,"MOD_DOG_GAIN_PTR  'name' is NULL, using _fs[0].\n");
    mdfs = mod_dog_fs[0];
  }

  if (eyeflag == 0)
    return mdfs->rg;
  else
    return mdfs->rrg;
}
/**************************************-**************************************/
/*                                                                           */
/*                                MOD_DOG_GAIN                               */
/*                                                                           */
/*  Compute the square of the response, and keep this for computing a        */
/*  gain control signal used by V1 cells that get LGN input.                 */
/*                                                                           */
/*****************************************************************************/
float ***mod_dog_gain(data,xn,yn,tn,mdfs)
     float ***data;
     int xn,yn,tn;
     struct mod_dog_filt *mdfs;
{
  int i,j,k;
  float ***g,*t,*tg,fhalf,fmax;

  fhalf = 0.5 * mdfs->fsum;  // Middle of range of possible values
  fmax = (mdfs->psum - mdfs->nsum)/2.0;  // 'nsum' is stored as negative value

  g = get_3d_farray(xn,yn,tn);
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      t = &(data[i+1][j+1][1]);
      tg = g[i][j];
      for(k=0;k<tn;k++){
	
	//  Absolute Value, expressed as fraction of max possible
	tg[k] = (t[k]-fhalf)/fmax;
	if (tg[k] < 0.0)
	  tg[k] *= -1.0;
	
	// Square
	// tg[k] = (t[k]-fhalf) * (t[k]-fhalf);  // Use gray reference
	
      }
    }
  }

  return g;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_DOG_GET_RESPONSE                           */
/*                                                                           */
/*  Return a pointer to 'mod_dog_r'.                                         */
/*                                                                           */
/*****************************************************************************/
float ***mod_dog_get_response(s,k,name,gain_flag)
     struct stim_struct *s;   // Stimulus data and param pair list
     int k;                   // repeat number
     char *name;              // population name, or NULL to use fs[0]
     int gain_flag;           // 1-compute a gain signal
{
  int xn,yn,tn,pflag;
  float ***dog,**dog_s,***rdog;
  //static int oldstimno = -1;
  struct mod_dog_filt *mdfs;

  pflag = 0;  // print flag
  if (myid == -1)
    pflag = 1;

  if (name != NULL)
    mdfs = mod_dog_fs_by_name(name);
  else{
    mylog(mylogf,"MOD_DOG_GET_RESPONSE  'name' is NULL, using _fs[0].\n");
    mdfs = mod_dog_fs[0];
  }

  //printf("===================> _fs name = %s\n",mdfs->name);
  //printf("              =====> _fs fsum = %f\n",mdfs->fsum);

  xn = mod_dog_xn;
  yn = mod_dog_yn;
  tn = mod_dog_tn;
  dog   = mdfs->f;
  dog_s = mdfs->s;

  if (mdfs->oldstimno != s->stimno){ // First time for this stimulus
    if (mdfs->r != NULL)
      free_f3tensor(mdfs->r,1,xn,1,yn,1,tn);

    // Check if some other LGN population has already used the stimulus,
    //   thus it is already in the Freq Domain...
    if (mod_dog_stim_si != s->stimno){

      // Free old storage
      if (mod_dog_stim_s != NULL)
	free_matrix(mod_dog_stim_s,1,xn,1,2*yn);

      // Next call modifies s->d
      mod_dog_stim_s = compute_response_single_from_ft_ret_s(s->d,dog,dog_s,xn,
							     yn,tn,pflag,
							     &rdog);

      // Gain signal, 23 July 2009
      if (gain_flag == 1){
	if (mdfs->rg != NULL){
	  free_3d_farray(mdfs->rg,xn,yn,tn);
	}
	mdfs->rg = mod_dog_gain(rdog,xn,yn,tn,mdfs); // 0 indexing
      }

      mod_dog_stim_si = s->stimno;
    }else{
      // Use 's->d' and 'mod_dog_stim_s', which are already in Freq. Domain
      compute_response_single_from_ft_2(s->d,mod_dog_stim_s,dog,dog_s,
					xn,yn,tn,pflag,&rdog);
    }

    mdfs->r = rdog;

    if (mod_dog_surr_sig > 0.0){
      mod_dog_surr_response(0); // 'mod_dog_surr' will be kept
    }

    if (mod_dog_cg_flag){ // Apply non-linearity for gain control
      mod_dog_cg01_response();
    }
    mdfs->oldstimno = s->stimno;
  }else{
    rdog = mdfs->r; // Use stored response to this stimulus
  }

  return rdog;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_DOG_GET_RESPONSE_BINOC                       */
/*                                                                           */
/*  Return a pointer to 'mod_dog_r' and 'mod_dog_rr'                         */
/*                                                                           */
/*****************************************************************************/
float ***mod_dog_get_response_binoc(s,k,name,rright)
     struct stim_struct *s;   // Stimulus data and param pair list
     int k;                   // repeat number
     char *name;              // population name, or NULL to use fs[0]
     float ****rright;
{
  int xn,yn,tn,pflag;
  float ***dog,**dog_s,***rdog,***rdogr;
  //static int oldstimno = -1;
  struct mod_dog_filt *mdfs;

  pflag = 0;  // print flag
  if (myid == -1)
    pflag = 1;

  if (name != NULL)
    mdfs = mod_dog_fs_by_name(name);
  else{
    mylog(mylogf,"MOD_DOG_GET_RESPONSE  'name' is NULL, using _fs[0].\n");
    mdfs = mod_dog_fs[0];
  }

  xn = mod_dog_xn;
  yn = mod_dog_yn;
  tn = mod_dog_tn;
  dog   = mdfs->f;
  dog_s = mdfs->s;
  
  if (mdfs->oldstimno != s->stimno){ // First time for this stimulus
    if (mdfs->r != NULL){
      free_f3tensor(mdfs->r,1,xn,1,yn,1,tn);
      free_f3tensor(mdfs->rr,1,xn,1,yn,1,tn);
    }

    // Next call modifies s->d
    compute_response_single_from_ft(s->d,dog,dog_s,xn,yn,tn,pflag,&rdog);
    mdfs->r = rdog;

    compute_response_single_from_ft(s->d_r,dog,dog_s,xn,yn,tn,pflag,&rdogr);
    mdfs->rr = rdogr;

    if (mod_dog_surr_sig > 0.0){
      mod_dog_surr_response(0); // 'mod_dog_surr' will be kept
      mod_dog_surr_response(1);
    }

    // Apply non-linearity for gain control
    if (mod_dog_cg_flag)
      exit_error("MOD_DOG_GET_RESPONSE_BINOC","CG not imp'd for binoc");

    mdfs->oldstimno = s->stimno;
  }else{
    rdog  = mdfs->r;  // Use stored response to this stimulus
    rdogr = mdfs->rr;
  }

  *rright = rdogr;
  return rdog;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MODEL_DOG_01_PREP_FT                           */
/*                                                                           */
/*****************************************************************************/
void model_dog_01_prep_ft()
{
  int i;
  int xn,yn,tn;
  struct mod_dog_filt *mdfs;

  mylog(mylogf,"  MODEL_DOG_01_PREP_FT\n");

  xn = mod_dog_xn;
  yn = mod_dog_yn;
  tn = mod_dog_tn;

  for(i=0;i<mod_dog_fsn;i++){

    mdfs = mod_dog_fs[i];

    sprintf(ggstr,"    Computing FFT of DOG filter, %s\n",mdfs->name);
    mylog(mylogf,ggstr);
    mdfs->s = matrix(1,xn,1,2*yn);
    contort_real_3d_farray(mdfs->f,1,1,1,xn,yn,tn);
    rlft3(mdfs->f,mdfs->s,xn,yn,tn,1);  // In nr_util.c
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MODEL_DOG_01_DONE                            */
/*                                                                           */
/*****************************************************************************/
void model_dog_01_done(m)
     struct model_struct *m; // Model params
{
  int i;
  int xn,yn,tn;
  struct mod_dog_filt *mdfs;

  mylog(mylogf,"  MODEL_DOG_01_DONE\n");

  xn = mod_dog_xn;
  yn = mod_dog_yn;
  tn = mod_dog_tn;

  if (1){ // New way

    //
    //  Free storage
    //
    for(i=0;i<mod_dog_fsn;i++){
      mod_dog_filt_struct_free_contents(mod_dog_fs[i]);
    }
    myfree(mod_dog_fs);

    //
    //  Set values to unused state
    //
    mod_dog_fs = NULL;
    mod_dog_fsn = 0;

  }else{ // Old way

    for(i=0;i<mod_dog_fsn;i++){

      mdfs = mod_dog_fs[i];

      free_f3tensor(mdfs->f,1,xn,1,yn,1,tn);
      free_matrix(mdfs->s,1,xn,1,2*yn);

      mdfs->f = (float ***)NULL;
      mdfs->s = (float **)NULL;

      if (mdfs->r != NULL)
	free_f3tensor(mdfs->r,1,xn,1,yn,1,tn);
      mdfs->r = (float ***)NULL;
    }
  }


  if (mod_dog_surr != NULL){
    free_f3tensor(mod_dog_surr,1,xn,1,yn,1,tn);
    free_matrix(mod_dog_surr_s,1,xn,1,2*yn);
  }
  if (mod_dog_sup != NULL){
    free_f3tensor(mod_dog_sup,1,xn,1,yn,1,tn);
  }

  mod_dog_surr = (float ***)NULL;
  mod_dog_surr_s = (float **)NULL;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_DOG_NL_LOGISTIC                            */
/*                                                                           */
/*****************************************************************************/
void mod_dog_nl_logistic(t,n,fsum)
     float *t;
     int n;
     float fsum;  // This is the mean reference value
{
  int i;
  float a,a2,m,fs2;

  a = mod_dog_logistic_asymp;
  a2 = 2.0 * a;
  //m = mod_dog_logistic_slope;
  fs2 = fsum / 2.0;

  //printf("______________________fsum = %f\n",fs2);

  m = 2.0 * mod_dog_logistic_slope / mod_dog_logistic_asymp;


  for(i=0;i<n;i++){
    t[i] = a2/(1.0 + exp(-m * (t[i]-fs2))) - a;
    t[i] += fs2;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_DOG_GET_RESP_TRACE                          */
/*                                                                           */
/*****************************************************************************/
float *mod_dog_get_resp_trace(r,xi,yi,tn,onflag,layname)
     float ***r;     // Raw response [1..xn][1..yn][1..tn]
     int xi,yi;      // Coordinates, 0..xn-1,  0..yn-1
     int tn;         // Time duration of returned response
     int onflag;     // 1-on 0-off
     char *layname;  // LGN layer name
{
  int i;
  float *t,*tp,c;
  struct mod_dog_filt *mdfs;

  if (layname != NULL)
    mdfs = mod_dog_fs_by_name(layname);
  else{
    mylog(mylogf,"MOD_DOG_GET_RESP_TRACE  'name' is NULL, using _fs[0].\n");
    mdfs = mod_dog_fs[0];
  }

  tp = r[xi+1][yi+1];  // This matches 0 origin in coords to 1 origin in 'r'

  t = (float *)myalloc(tn*sizeof(float));

  if (onflag == 1){ // ON
    for(i=0;i<tn;i++)
      t[i] = tp[i+1];
  }else{
    for(i=0;i<tn;i++) // 1 - stim (stim should be [0,1] bounded)
      //t[i] = mod_dog_fs[0]->fsum - tp[i+1];
      t[i] = mdfs->fsum - tp[i+1];
  }


  /***
  if ((xi == mod_dog_xn/2) && (yi == mod_dog_yn/2)){
    printf("*************** dumping here...........\n");
    append_farray_plot("zz.raw.pl","raw",t,tn,1);
    }***/

  
  
  //  Logistic non-linearity for Michael
  if (mod_dog_logistic_name != NULL){
    mod_dog_nl_logistic(t,tn,mdfs->fsum);
  }


  if (mdfs->power != 1.0){
    printf(" POWER HERE %f\n",mdfs->power);
    pow_farray(t,tn,mdfs->power);
  }
      

  //
  //  Apply suppression
  //
  if (mod_dog_surr_sig > 0.0){

    //printf("WYETH NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO  777\n");

    tp = mod_dog_sup[xi+1][yi+1];
    for(i=0;i<tn;i++){
      c =  8.0 * tp[i+1]; // WYETH 
      if (c < 0.0)
	c = 0.0;
      t[i] *= 1.0 - c;
    }
    /***
    if (myid == -1){
      float *ts;
      
      ts = &(mod_dog_sup[xi][yi][1]);
      append_farray_plot("zzz.surr.pl","surr_xi_yi",ts,tn,1);
      }***/
  }

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MODEL_DOG_01_GET_RESPONSE                         */
/*                                                                           */
/*  Calling sequence:                                                        */
/*    1.  model_dog_01_prep                                                  */
/*    2.  model_dog_01_prep_ft                                               */
/*    3.  this routine                                                       */
/*        ...                                                                */
/*    N.  model_dog_01_done                                                  */
/*                                                                           */
/*  NOTES:                                                                   */
/*    s->d, the stimulus data, is modified here.                             */
/*                                                                           */
/*****************************************************************************/
void model_dog_01_get_response(m,s,r,k)
     struct model_struct *m;        /* Model parameter pair list. */
     struct stim_struct *s;         /* Stimulus data and param pair list */
     struct response_struct *r;     /* Response data */
     int k;                         /* Repeat number */
{
  int i,xi,yi,ri;
  int xn,yn,tn,ns,vflag,*seedlist,rn; // OLD REMOVE: ,gtsi;
  float ***rdog,*t,*fs,sampling;
  char *outfile,tname[LONG_SLEN],*prefix,tstr[LONG_SLEN4],*sgen,*lname;
  char *dump_name;
  static int oldstimno = -1;
  float *gti,*gtx,samp,*rv;
  struct onode *po,*spko;

  // Get outfile prefix
  prefix = paramfile_get_char_param_default(r->ppl,"outfile","zz");

  // Get name for plots for this trial
  if (s->nvar > 0)
    strcpy(tname,s->name);  // Name specifies var param values
  else
    strcpy(tname,"raw");    // If no var params, name is ""

  xn = mod_dog_xn;  // Convenient values
  yn = mod_dog_yn;
  tn = mod_dog_tn;

  // Determine layer name, e.g., "lgn", set onode pointers
  if (m->ppl == NULL){
    lname = strdup("lgn"); // *** WYETH - assuming pop name is LGN
    po = onode_get_node_type_item_val(m->o,"pop","name",lname);
    spko = onode_child_get_unique(po,"spike_gen");

    if (onode_test_int(po,"ifc_dump_center",1)){
      dump_name = tname;
    }else
      dump_name = NULL;

  }else{
    lname = NULL;
    po = NULL;
    spko = NULL;
  }

  /***
  if (myid == -1){
    if (oldstimno != s->stimno){ // First time for this stimulus
      outfile =paramfile_make_fname(m->ppl,"write_raw_stim",prefix,".stim.pl");
      if (outfile != (char *)NULL){
	exit_error("MODEL_DOG_01_GET_RESPONSE","UNDEFINED xi, yi WYETH");
	t = get_1d_from_3d_farray(s->d,xi,1,yi,1,1,tn);
	append_farray_plot(outfile,tname,t,tn,1);
	myfree(t);
	myfree(outfile);
      }
      oldstimno = s->stimno;
    }
  }
  ***/

  // Next call modifies s->d
  rdog = mod_dog_get_response(s,k,lname,0);// Set 'mod_dog_r' t'which rdog pnts

  if (m->ppl == NULL)
    outfile = onode_make_fname(po,"write_3d_response",prefix,".resp.3d");
  else
    outfile = paramfile_make_fname(m->ppl,"write_3d_response",prefix,
				   ".resp.3d");

  if (outfile != (char *)NULL){
    write_3d_data_part(outfile,rdog,1,xn,1,yn,1,tn,4,2,1);
    myfree(outfile);
  }

  samp = 1.0/mod_dog_tscale;

  // Pick randomization seed for each response
  seedlist = get_seeds(m->mseed[r->tsi],1000000,r->n);

  for(ri=0;ri<r->n;ri++){
    if ((strcmp(r->plname[ri],"lgn")==0)||(strcmp(r->plname[ri],"lgn_on")==0)){
      xi = r->xi[ri];
      yi = r->yi[ri];

      if ((xi < 0) || (xi >= xn) || (yi < 0) || (yi >= yn))
	mylog_exit(mylogf,
		   "MODEL_DOG_01_GET_RESPONSE  response x,y out of bounds\n");
      else{
	sprintf(ggstr,"   LGN %d %d\n",xi,yi);
	mylog(mylogf,ggstr);
      }

      // May do gain control and/or suppression ?
      t = mod_dog_get_resp_trace(rdog,xi,yi,tn,1,lname);  // 1 = onflag

      
      //for(i=0;i<tn;i++) // Get response at coord
      //t[i] = rdog[xi+1][yi+1][1+i];


      // pump the drive through the non-linear filter
      if (mod_dog_nl_flag == 1){
	if (mod_dog_nl_filt == NULL) // Only do this once
	  mod_dog_nonlin1_prep(m); // Blake
	mod_dog_nonlin1(m,t,tn); // Blake
      }

      // squash the signal
      if (mod_dog_squash_flag)
	mod_dog_squash(t,tn);
      else if (mod_dog_sq_flag)  // Probably just using this for now
	mod_dog_squash_blake(m,t,tn,1);

      /*** WYETH GAIN CONTROL HERE ??? - only for mod_dog run.  ***/

      if (myid == -1){
	if (m->ppl == NULL)
	  outfile = onode_make_fname(po,"write_raw",prefix,".raw.pl");
	else
	  outfile = paramfile_make_fname(m->ppl,"write_raw",prefix,".raw.pl");
	if (outfile != (char *)NULL){
	  sprintf(tstr,"%s_%d_%d",tname,xi,yi);
	  append_farray_plot(outfile,tstr,t,tn,1);
	  myfree(outfile);
	}
      }

      /*** SPIKE GEN *********************************************
	   INPUT
	   t[tn] - array holding raw input
	   samp  - samples per second for 't'
	   OUTPUT
	   fs[ns] - Array of float spike times
	   rv[rn] - Array of float data - Vm or Poisson Prob
      ***/

      vflag = 0;
      if (r->nf > 0)
	vflag = 1;

      if (m->ppl == NULL)
	sgen = onode_getpar_chr_exit(spko,"type");
      else
	sgen = paramfile_get_char_param_default(m->ppl,"spike_gen","poisson");
      
      if (strcmp(sgen,"ifc")==0){

	// Determine whether to return - WYETH does this work anymore?
	gtx = t;

	gti = get_zero_farray(tn);

	if (m->ppl != NULL){
	  // WYETH - should arrange to use 'dump_name' here?
	  ifc_test(myid,m,"ifc",gtx,gti,samp,tn,seedlist[ri],1,tname,
		   &fs,&ns,vflag,&rv,&rn);
	}else{
	  ifc_test(myid,m,"lgn",gtx,gti,samp,tn,seedlist[ri],1,dump_name,
		   &fs,&ns,vflag,&rv,&rn);
	}

	myfree(gti);
      }else if (strcmp(sgen,"poisson")==0){
	ifc_util_poisson(mylogf,m,spko,t,tn,samp,seedlist[ri],0,vflag,&fs,&ns,
			 &rv,&rn);
      }else
	mylog_exit(mylogf,"MODEL_DOG_01_GET_RESPONSE  Unknown spike gen\n");
      myfree(sgen);
      /********************************************/
	
      sampling = r->samp[ri]; // Adjust from ms to 'save_...' sampling
      multiply_farray(fs,ns,sampling/1000.0);

      add_const_farray(fs,ns,mod_dog_fs[0]->tdelay*sampling);// Retinal latency

      //gtsi = m->m_i * s->ntr + r->tsi; // Global trial sequence index
      //printf("m->m_i = %d\n",m->m_i);

      //sprintf(ggstr,"(mod_dog_u) m->m_i %d  s->ntr %d  r->tsi %d  gtsi %d\n",
      //m->m_i,s->ntr,r->tsi,gtsi);
      //mylog(mylogf,ggstr);

      //r->s[ri][r->tsi] = fs;
      //r->cnt[ri][r->tsi] = ns;
      r->s[ri][r->gtsi] = fs;
      r->cnt[ri][r->gtsi] = ns;
      
      if (r->nf > 0)
	mylog_exit(mylogf,"MODEL_DOG_01_GET_RESPONSE  - not re-impl'd yet\n");
      if (r->nf > 0){
	mylog(mylogf,"    Storing Vm\n");
	//r->f[0][r->tsi] = rv;
	//r->fcnt[0][r->tsi] = rn;
	r->f[0][r->gtsi] = rv;
	r->fcnt[0][r->gtsi] = rn;
      }
      myfree(t);
    }else if (strcmp(r->plname[ri],"lgn_off")==0){
      mylog_exit(mylogf,
		 "MODEL_DOG_01_GET_RESPONSE  'lgn_off' not yet implemented\n");
    }
  }

  myfree(seedlist);
  myfree(prefix); // Prefix for outfile names
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_DOG_GET_SPIKES_2D_FLAG                       */
/*                                                                           */
/*  Get spike trains for the DOG response at any element where flag = 1.     */
/*                                                                           */
/*  NOTES                                                                    */
/*  - Set 'ss' to 1.0 for O cells, -1.0 for OFF cells.                       */
/*  - Spikes are returned in 'msec'                                          */
/*                                                                           */
/*****************************************************************************/
void mod_dog_get_spikes_2d_flag(m,r,rdog,lname,onflag,xn,yn,tn,flag,samp,seed,
				rs2,rcnt2,sav_flag,rx,ri,rv,rvn)
     struct model_struct *m;        // Model parameter pair list
     struct response_struct *r;     // Response data
     float ***rdog;                 // [1..xn][1..yn][1..tn]
     char *lname;                   // LGN pop layer name
     int onflag;                    // 1=ON, 0=OFF
     int xn,yn,tn,**flag;           // 
     float samp;                    // Sampling (1/tscale), for 'rdog'
     int seed;                      // Randomization seed
     float ****rs2;                 // [xn][yn][rcnt2] spikes at flagged coords
     int ***rcnt2;                  // [xn][yn] num. spikes at flagged coords
     int sav_flag;                  // Number of saves in this layer
     float ****rx,****ri,****rv;    // Returned sav values [xn][yn][tn]
     int *rvn;                      // Length of 'rv' trace (ms)
{
  int i,j,k;
  int **cnt2,flag_count,spike_count,vn,savit,sav_vn;
  int tseed,dumpx,dumpy,*seedlist;
  float ***s2,*t,*gti,*gtx,*v;
  float ***sav_vm,***sav_gx,***sav_gi;
  char pname[SLEN],prefix[SLEN],tstr[SLEN],*dumpname,*dumpname_ptr,*sgen;
  struct onode *po,*spko;
  struct mod_dog_filt *mdfs;

  // Get a pointer to the filter structure for layer 'lname'
  if (lname != NULL)
    mdfs = mod_dog_fs_by_name(lname);
  else{
    //mylog(mylogf,"MOD_DOG_GET_RESP_TRACE  'name' is NULL, using _fs[0].\n");
    mdfs = mod_dog_fs[0];
  }


  // Get 'sgen' - spike generation algorithm
  if (m->ppl != NULL){
    strcpy(prefix,"ifc");
    sgen = paramfile_get_char_param_default(m->ppl,"spike_gen","poisson");
  }else{
    strcpy(prefix,lname);
    po = onode_get_node_type_item_val(m->o,"pop","name",lname);
    if (po == NULL)
      mylog_exit(mylogf,"MOD_DOG_GET_SPIKES_2D_FLAG  No lgn layer.");
    spko = onode_child_get_unique(po,"spike_gen");
    sgen = onode_getpar_chr_exit(spko,"type");
  }

  // Determine whether to dump output for an LGN unit
  dumpname = (char *)NULL;
  dumpx = dumpy = -1;
  if (myid == -1){
    strcpy(pname,"dog_ifc_dump_coord");
    if (m->ppl != NULL){
      if (paramfile_test_param(m->ppl,pname)){
	dumpx = paramfile_get_nth_int_param_or_exit(m->ppl,pname,0);
	dumpy = paramfile_get_nth_int_param_or_exit(m->ppl,pname,1);
	sprintf(ggstr,"    Dumping LGN results for unit %d %d\n",dumpx,dumpy);
	mylog(mylogf,ggstr);
	if (onflag)
	  sprintf(tstr,"%s_ON_%d_%d",lname,dumpx,dumpy);
	else
	  sprintf(tstr,"%s_OFF_%d_%d",lname,dumpx,dumpy);
	dumpname = strdup(tstr);
      }
    }else{
      // WYETH - HACK, pname is set above
      if ((onode_test_ostr(po,pname))||
	  onode_test_int(po,"ifc_dump_center",1)){
	//(onode_test_ostr(po,"ifc_dump_coord"))){
	//
	//  WYETH - THIS IS A HACK TO FORCE MIDDLE unit to be dumped
	//
	dumpx = xn/2;
	dumpy = yn/2;

	sprintf(ggstr,"    Dumping LGN results for unit %d %d\n",dumpx,dumpy);
	mylog(mylogf,ggstr);
	if (onflag)
	  sprintf(tstr,"%s_ON_%d_%d",lname,dumpx,dumpy);
	else
	  sprintf(tstr,"%s_OFF_%d_%d",lname,dumpx,dumpy);
	dumpname = strdup(tstr);
      }
    }
  }


  cnt2 = get_zero_2d_iarray(xn,yn);
  s2 = (float ***)myalloc(xn*sizeof(float **));
  for(i=0;i<xn;i++){
    s2[i] = (float **)myalloc(yn*sizeof(float *));
  }

  gti = get_zero_farray(tn);  // Set inhibitory input to zero

  // WYETH NEW JAN 13 2009, seed can change w/i trial for multiple LGN pops
  if (onflag)
    seedlist = get_seeds(seed,1000000,xn*yn);
  else
    seedlist = get_seeds(17*seed,1000000,xn*yn);

  if (mod_dog_nl_flag == 1){
    if (mod_dog_nl_filt == NULL) // Only do this once
      mod_dog_nonlin1_prep(m);   // Blake
  }
  
  // Save traces for possible response storage
  if (sav_flag > 0){
    sav_vm = get_3d_farray(xn,yn,0);
    sav_gx = get_3d_farray(xn,yn,0);
    sav_gi = get_3d_farray(xn,yn,0);
  }else{
    sav_vm = sav_gx = sav_gi = NULL;
    savit = 0;
  }

  flag_count = 0;
  spike_count = 0;
  t = (float *)myalloc(tn*sizeof(float));
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      if (flag[i][j]){

	if (sav_flag > 0){
	  if (flag[i][j] > 1)
	    savit = 1;
	  else
	    savit = 0;
	}

	// Will apply suppression, if applicable, and match 0 vs. 1 origin
	t = mod_dog_get_resp_trace(rdog,i,j,tn,onflag,lname);



	  /*
	  {
	    int k;
	    int zn;
	    float *zz,*zzc;
	    char tstr[SLEN];

	    zz = get_zero_farray(tn);
	    for(k=0;k<tn;k++){
	      zz[k] = 1.0/(float)tn *   ((float)k - tn/2.0);
	    }
	    zzc = copy_farray(zz,tn);

	    mod_dog_nl_logistic(zz,tn);

	    sprintf(tstr,"Logistic_slope_%f__asymp_%f",mod_dog_logistic_slope,
		    mod_dog_logistic_asymp);
	    append_farray_xy_plot("zz.log.pl",zzc,zz,tn,tstr);

	    exit(0);
	    }
	  */

	/*
	if ((i == xn/2) && (j == yn/2)){
	  append_farray_plot("zz.raw.pl","log",t,tn,1);
	  }*/


	// Non-linear processing
	if (mod_dog_nl_flag == 1){
	  mod_dog_nonlin1(m,t,tn); // Blake
	}


	// squash the signal
	if (mod_dog_squash_flag){
	  if ((i==dumpx) && (j==dumpy))   // Dumping impl'd for IFC
	    append_farray_plot("zz.squash.dump.pl","before",t,tn,1);
	  mod_dog_squash(t,tn);

	  if ((i==dumpx) && (j==dumpy))   // Dumping impl'd for IFC
	    append_farray_plot("zz.squash.dump.pl","after",t,tn,1);
	  
	}else if (mod_dog_sq_flag){
	  mod_dog_squash_blake(m,t,tn,0);
	}

	gtx = t;

	tseed = seedlist[i*yn + j];

	if (strcmp(sgen,"poisson")==0){
	  ifc_util_poisson(mylogf,m,spko,gtx,tn,samp,tseed,0,0,&(s2[i][j]),
			   &(cnt2[i][j]),&v,&vn);

	}else if (strcmp(sgen,"ifc")==0){
	  if ((i==dumpx) && (j==dumpy)){   // Dumping impl'd for IFC
	    dumpname_ptr = dumpname;
	  }else
	    dumpname_ptr = (char *)NULL;

	  ifc_test(myid,m,prefix,gtx,gti,samp,tn,tseed,0,dumpname_ptr,
		   &(s2[i][j]),&(cnt2[i][j]),savit,&v,&vn);

	  if (savit){
	    sav_gx[i][j] = copy_farray(gtx,tn);
	    sav_gi[i][j] = copy_farray(gti,tn);
	    sav_vm[i][j] = v;
	    sav_vn = vn;
	  }


	}else{
	  mylog_exit(mylogf,"MOD_DOG_GET_SPIKES_2D_FLAG  Unknown spike gen\n");
	}

	// Add a constant delay to all spike times WYETH - spikes in msec
	add_const_farray(s2[i][j],cnt2[i][j],mdfs->tdelay*1000.0);

	flag_count += 1;
	spike_count += cnt2[i][j];

	myfree(t);  // WYETH - Jun 30, 2009, moved here from end of routine

      }else{ // Flag was not set
	s2[i][j] = NULL;  // Important for freeing s2 later
	cnt2[i][j] = -1;  // -1 indicates spikes never computed

	if ((i==dumpx) && (j==dumpy)){
	  mylog(mylogf,"  *** LGN cell to be dumped not selected as input\n");
	}
      }
    }
  }

  sprintf(ggstr,
	  "      Spikes generated for %d of %d positions, %.4f spk/cell\n",
	  flag_count,xn*yn,(float)spike_count/(float)flag_count);
  mylog(mylogf,ggstr);

  if (dumpname != (char *)NULL)
    myfree(dumpname);
  myfree(seedlist);
  myfree(gti);
  myfree(sgen);

  if (sav_flag > 0){
    *rx = sav_gx;
    *ri = sav_gi;
    *rv = sav_vm;
    *rvn = sav_vn;
  }else{
    *rv = *rx = *ri = NULL;
    *rvn = -1;
  }

  *rs2 = s2; *rcnt2 = cnt2;
}
/**************************************-**************************************/
/*                                                                           */
/*                       MOD_DOG_GET_GT_PROD_PAIR_SPIKES                     */
/*                                                                           */
/*  Multiply two spike trains after smoothing both and delaying 's2'.        */
/*                                                                           */
/*  - Take two spikes trains 's1' and 's2'                                   */
/*  - Add 'delay' to spike times of 's2'                                     */
/*  - Convolve spike trains with 'mask'                                      */
/*  - Multiply the smoothed spike trains                                     */
/*                                                                           */
/*****************************************************************************/
void mod_dog_get_gt_prod_pair_sarray(s1,n1,s2,n2,dur,delay,mask,nm,rgt)
     float *s1;
     int n1;
     float *s2;
     int n2;
     int dur,delay;
     float *mask;
     int nm;
     float **rgt;               /* Conductance vs. time */
{
  int i;
  int t;
  float *sx1,*sx2,*g1,*g2,*gt;

  sx1 = get_zero_farray(dur);
  sx2 = get_zero_farray(dur);
  for(i=0;i<n1;i++){
    t = my_rint(s1[i]);
    if ((t >= 0)&&(t < dur)){
      sx1[t] += 1.0;
    }
  }
  for(i=0;i<n2;i++){
    t = delay + my_rint(s2[i]);  /* Delay the second spike train */
    if ((t >= 0)&&(t < dur)){
      sx2[t] += 1.0;
    }
  }

  g1 = convolve_with_mask_causal(sx1,dur,mask,nm);
  g2 = convolve_with_mask_causal(sx2,dur,mask,nm);
  gt = get_farray(dur);
  for(i=0;i<dur;i++)
    gt[i] = g1[i] * g2[i];

  myfree(g1); myfree(g2);
  myfree(sx1); myfree(sx2);

  *rgt = gt;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_DOG_GET_GT_SARRAY                          */
/*                                                                           */
/*****************************************************************************/
void mod_dog_get_gt_sarray(m,s,cnt,n,t1,t2,amp,rgt,rnms)
     struct model_struct *m;    /* Model parameter pair list. */
     float **s;
     int *cnt,n;
     float t1,t2,amp;
     float **rgt;               /* Conductance vs. time */
     int *rnms;
{
  int i,j;
  int nx,t,tn,nms;
  float *xshape,*gt,*g;

  nx = (int)(t1 * 7.0);        /* Length of mask is 7 times decay time */
  xshape = diff_exp_farray(0.0,amp,t1,t2,nx);

  tn = mod_dog_tn;
  nms = (int)((float)tn * 1000.0*mod_dog_tscale); /* Duration in ms */
  gt = get_zero_farray(nms);
  for(i=0;i<n;i++){
    for(j=0;j<cnt[i];j++){
      t = my_rint(s[i][j]);
      if ((t >= 0)&&(t < nms)){
	gt[t] += 1.0;
      }
    }
  }
  if (myid == -1)
    append_farray_plot("zz.simp_out.pl","spikes_out",gt,nms,1);

  g = convolve_with_mask_causal(gt,nms,xshape,nx);
  myfree(gt); myfree(xshape);
  gt = g;

  if (myid == -1)
    append_farray_plot("zz.simp_out.pl","gtx_to_complex",g,nms,1);

  *rgt = gt; *rnms = nms;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_DOG_GET_GTX_GTI                          */
/*                                                                           */
/*  Combine some spike trains from s0 and s1 (depending on 'inw'), and       */
/*  convolve with EPSC shape.                                                */
/*                                                                           */
/*  - Returned arrays are in time units of milliseconds.                     */
/*                                                                           */
/*****************************************************************************/
void mod_dog_get_gtx_gti(m,s0,s1,cnt0,cnt1,nin,inx,iny,inw,rgtx,rgti,rnms,cn)
     struct model_struct *m;    /* Model parameter pair list. */
     float ***s0,***s1;         /* Spike times in ms */
     int **cnt0,**cnt1;
     int nin,*inx,*iny;
     float *inw;
     float **rgtx,**rgti;
     int *rnms;
     int cn;
{
  int i,j;
  int nx,ni,xi,yi,t,cnt,tn,nms;
  float *xshape,*ishape,*s,t1x,t2x,ampx,t1i,t2i,ampi,*gtx,*gti,*g;
  char tstr[SLEN];

  /*** Time constants in msec ***/
  t1x   = paramfile_get_float_param_or_exit(m->ppl,"simp_ex_tau_f");
  t2x   = paramfile_get_float_param_or_exit(m->ppl,"simp_ex_tau_r");
  ampx  = paramfile_get_float_param_or_exit(m->ppl,"simp_ex_amp");
  t1i   = paramfile_get_float_param_or_exit(m->ppl,"simp_in_tau_f");
  t2i   = paramfile_get_float_param_or_exit(m->ppl,"simp_in_tau_r");
  ampi  = paramfile_get_float_param_or_exit(m->ppl,"simp_in_amp");

  nx = (int)(t1x * 7.0);        /* Length of mask is 7 times decay time */
  xshape = diff_exp_farray(0.0,ampx,t1x,t2x,nx);
  ni = (int)(t1i * 7.0);        /* Length of mask is 7 times decay time */
  ishape = diff_exp_farray(0.0,ampi,t1i,t2i,ni);

  tn = mod_dog_tn;
  nms = (int)((float)tn * 1000.0*mod_dog_tscale); /* Duration in ms */
  
  gtx = get_zero_farray(nms);
  gti = get_zero_farray(nms);

  for(i=0;i<nin;i++){
    xi = inx[i];
    yi = iny[i];
    if (inw[i] < 0.0){  /*** WYETH - only sign of weight is used !!! ***/
      s = s0[xi][yi];
      cnt = cnt0[xi][yi];
    }else{
      s = s1[xi][yi];
      cnt = cnt1[xi][yi];
    }
    for(j=0;j<cnt;j++){
      t = my_rint(s[j]);
      if ((t >= 0)&&(t < nms)){
	gtx[t] += 1.0;
      }
    }
  }
  if (myid == -1){
    append_farray_plot("zz.simp.pl","gtx_spikes",gtx,nms,1);
    if (cn == 0)
      append_farray_plot("zz.simp.pl","xshape",xshape,nx,1);
  }

  g = convolve_with_mask_causal(gtx,nms,xshape,nx);
  myfree(gtx);
  gtx = g;

  if (myid == -1){
    sprintf(tstr,"gtx_%d",cn);
    append_farray_plot("zz.simp.pl",tstr,gtx,nms,1);
  }

  myfree(xshape); myfree(ishape);

  *rgtx = gtx; *rgti = gti; *rnms = nms;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MODEL_DOG_SIMP_01_PREP                         */
/*                                                                           */
/*****************************************************************************/
void model_dog_simp_01_prep(mppl,r)
     struct param_pair_list *mppl; /* Model parameter pair list. */
     struct response_struct *r;    /* Response params. */
{
  int i,j,k;
  int xn,yn,nsamp,seed,xi,yi;
  int **flag;
  float **rf,sscale,theta,ssd,sf,ph,min,max,fsamp;

  mylog(mylogf,"  MODEL_DOG_SIMP_01_PREP\n");

  sscale = mod_dog_sscale;
  xn = mod_dog_xn;
  yn = mod_dog_yn;
  flag = get_zero_2d_iarray(xn,yn);   /* Has this cell been chosen yet? */

  ssd   = paramfile_get_float_param_or_exit(mppl,"simp_ssd");
  sf    = paramfile_get_float_param_or_exit(mppl,"simp_sf");
  ph    = paramfile_get_float_param_or_exit(mppl,"simp_phase");
  theta = paramfile_get_float_param_or_exit(mppl,"simp_ori");
  fsamp = paramfile_get_float_param_or_exit(mppl,"simp_rf_fsamp");
  seed  = paramfile_get_int_param_or_exit(mppl,"simp_rf_seed");
  nsamp = paramfile_get_int_param_or_exit(mppl,"simp_rf_nsamp");
  if (seed > 0)
    seed = -seed;

  mod_dog_simp_n = 1;
  mod_dog_simp_nin = get_zero_iarray(1);
  mod_dog_simp_nin[0] = nsamp;
  mod_dog_simp_x = get_2d_iarray(1,nsamp);
  mod_dog_simp_y = get_2d_iarray(1,nsamp);
  mod_dog_simp_w = get_2d_farray(1,nsamp);
  mod_dog_simp_rf = (float ***)myalloc(sizeof(float **));
  mod_dog_simp_flag = (int ***)myalloc(sizeof(int **));
  mod_dog_simp_flag[0] = flag;

  rf = gabor_2d_space(xn,yn,sscale,ssd,sf,theta,ph); /* in kernel_util.c */
  mod_dog_simp_rf[0] = rf;
  max = max_of_2d_farray(rf,xn,yn);
  min = min_of_2d_farray(rf,xn,yn);
  sprintf(ggstr,"  min, max =  %f %f\n",min,max);
  mylog(mylogf,ggstr);

  /*** Select inputs from 2D spatial array of DOG responses ***/
  k = 0;
  while(k < nsamp){
    xi = (int)((float)xn * myrand_util_ran2(&seed));
    yi = (int)((float)yn * myrand_util_ran2(&seed));
    if (flag[xi][yi] == 0){
      if ((rf[xi][yi] > fsamp*max)|| (rf[xi][yi] < fsamp*min)){
	flag[xi][yi] = 1;
	mod_dog_simp_x[0][k] = xi;
	mod_dog_simp_y[0][k] = yi;
	mod_dog_simp_w[0][k] = rf[xi][yi];
	k += 1;
      }
    }
  }

  for(i=0;i<mod_dog_simp_nin[0];i++){
    sprintf(ggstr,"%d (%d,%d) %.2f\n",i,mod_dog_simp_x[0][i],
	    mod_dog_simp_y[0][i],mod_dog_simp_w[0][i]);
    mylog(mylogf,ggstr);
  }

  for(i=0;i<yn;i++){
    for(j=0;j<xn;j++)
      if (flag[i][j])
	if (rf[i][j] > 0.0)
	  mylog(mylogf," O");
	else
	  mylog(mylogf," X");
      else
	mylog(mylogf," .");
    mylog(mylogf,"\n");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                       MODEL_DOG_SIMP_01_GET_RESPONSE                      */
/*                                                                           */
/*  Calling sequence:                                                        */
/*    1.   model_dog_01_prep                                                 */
/*    2.   model_dog_01_prep_ft                                              */
/*    3.   model_dog_simp_01_prep                                            */
/*    4.   this routine                                                      */
/*        ...                                                                */
/*    N-1. model_dog_01_done                                                 */
/*    N.   model_dog_simp_01_done                                            */
/*                                                                           */
/*  NOTES:                                                                   */
/*    s->d, the stimulus data, is modified here.                             */
/*                                                                           */
/*****************************************************************************/
void model_dog_simp_01_get_response(m,s,r,k)
     struct model_struct *m;        /* Model parameter pair list. */
     struct stim_struct *s;         /* Stimulus data and param pair list */
     struct response_struct *r;     /* Response data */
     int k; /* Repeat number */
{
  int xn,yn,tn,ns,vn;
  int **cnt0,**cnt1,*inx,*iny,nin,**flag,gnms,seed,tvn;
  float ***rdog,*fs,sampling,***s0,***s1,*inw,*gtx,*gti,sampms,*v,samp;
  float ***tgx,***tgi,***tgv;
  char *outfile;

  xn = mod_dog_xn;
  yn = mod_dog_yn;
  tn = mod_dog_tn;

  // Next call modifies s->d
  rdog = mod_dog_get_response(s,k,NULL,0); // Set 'mod_dog_r' t'which rdog pnts

  outfile = paramfile_get_char_param_default(m->ppl,"write_3d_response","");
  if (strcmp(outfile,"")!=0)
    write_3d_data_part(outfile,rdog,1,xn,1,yn,1,tn,4,2,1);
  myfree(outfile);

  /*** Get spike trains for all relevant DOG grid positions ***/
  samp = 1.0/mod_dog_tscale;
  flag = mod_dog_simp_flag[0];

  seed = m->mseed[r->tsi];
  mod_dog_get_spikes_2d_flag(m,r,rdog,NULL,1,xn,yn,tn,flag,samp,seed,
			     &s1,&cnt1,0,&tgx,&tgi,&tgv,&tvn); //on
  mod_dog_get_spikes_2d_flag(m,r,rdog,NULL,0,xn,yn,tn,flag,samp,seed,
			     &s0,&cnt0,0,&tgx,&tgi,&tgv,&tvn); //off

  /*** Use spike trains to get excitatory and inhibitory conductances ***/
  nin = mod_dog_simp_nin[0];
  inx = mod_dog_simp_x[0];
  iny = mod_dog_simp_y[0];
  inw = mod_dog_simp_w[0];
  mod_dog_get_gtx_gti(m,s0,s1,cnt0,cnt1,nin,inx,iny,inw,&gtx,&gti,&gnms,0);
  /*** Data returned in msec time units ***/

  /*** Get spike trains for this simple cell ***/
  sampms = 1000.0;
  ifc_test(myid,m,"ifc_simp",gtx,gti,sampms,gnms,m->mseed[r->tsi],1,
	   (char *)NULL,&fs,&ns,0,&v,&vn);
  myfree(gti); myfree(gtx);

  /* Spikes come in ms, adjust based on the 'save_spikes' sampling */
  /*sampling = paramfile_get_nth_float_param_or_exit(m->ppl,"save_spikes",1);*/
  sampling = mod_util_resp_get_samp(r,"spikes");
  multiply_farray(fs,ns,sampling/1000.0);
  
  /*** Store response, link to storage, copyflag = 0 ***/
  (void)mod_util_resp_check_store_s(r,"spikes",mylogf,fs,ns,0);

  /**if (r->ns > 0){
     r->s[0][r->tsi] = fs;
     r->cnt[0][r->tsi] = ns; * 'ti' is updated in 'run_gen_tuning_curve' *
     }**/
}
/**************************************-**************************************/
/*                                                                           */
/*                            MODEL_DOG_COMP_01_PREP                         */
/*                                                                           */
/*****************************************************************************/
void model_dog_comp_01_prep(mppl,r)
     struct param_pair_list *mppl; /* Model parameter pair list. */
     struct response_struct *r;    /* Response params. */
{
  int i,j,k;
  int xn,yn,nsamp,seed,xi,yi,n;
  int **flag;
  float **rf,sscale,theta,ssd,sf,ph,min,max,fsamp,dtheta;

  mylog(mylogf,"  MODEL_DOG_COMP_01_PREP\n");

  sscale = mod_dog_sscale;
  xn = mod_dog_xn;
  yn = mod_dog_yn;

  ssd   = paramfile_get_float_param_or_exit(mppl,"simp_ssd");
  sf    = paramfile_get_float_param_or_exit(mppl,"simp_sf");
  ph    = paramfile_get_float_param_or_exit(mppl,"simp_phase");
  theta = paramfile_get_float_param_or_exit(mppl,"simp_ori");
  fsamp = paramfile_get_float_param_or_exit(mppl,"simp_rf_fsamp");
  seed  = paramfile_get_int_param_or_exit(mppl,"simp_rf_seed");
  nsamp = paramfile_get_int_param_or_exit(mppl,"simp_rf_nsamp");
  if (seed > 0)
    seed = -seed;

  n = paramfile_get_int_param_or_exit(mppl,"comp_nsimp");
  mod_dog_simp_n = n;
  mod_dog_simp_nin = get_const_iarray(n,nsamp); /* Each cell has 'nsamp' inp */
  mod_dog_simp_x = get_2d_iarray(n,nsamp);
  mod_dog_simp_y = get_2d_iarray(n,nsamp);
  mod_dog_simp_w = get_2d_farray(n,nsamp);
  mod_dog_simp_rf = (float ***)myalloc(n*sizeof(float **));
  mod_dog_simp_flag = get_zero_3d_iarray(n,xn,yn);

  dtheta = 360.0/(float)n;

  for(i=0;i<n;i++){
    rf = gabor_2d_space(xn,yn,sscale,ssd,sf,theta,ph+(float)i*dtheta);
    mod_dog_simp_rf[i] = rf;
    
    max = max_of_2d_farray(rf,xn,yn);
    min = min_of_2d_farray(rf,xn,yn);
    sprintf(ggstr,"  min, max =  %f %f\n",min,max);
    mylog(mylogf,ggstr);
    
    /*** Select inputs from 2D spatial array of DOG responses ***/
    flag = mod_dog_simp_flag[i];
    k = 0;
    while(k < nsamp){
      xi = (int)((float)xn * myrand_util_ran2(&seed));
      yi = (int)((float)yn * myrand_util_ran2(&seed));
      if (flag[xi][yi] == 0){
	if ((rf[xi][yi] > fsamp*max)|| (rf[xi][yi] < fsamp*min)){
	  flag[xi][yi] = 1;
	  mod_dog_simp_x[i][k] = xi;
	  mod_dog_simp_y[i][k] = yi;
	  mod_dog_simp_w[i][k] = rf[xi][yi];
	  k += 1;
	}
      }
    }
    
    for(j=0;j<mod_dog_simp_nin[i];j++){
      sprintf(ggstr,"%d (%d,%d) %.2f\n",j,mod_dog_simp_x[i][j],
	      mod_dog_simp_y[i][j],mod_dog_simp_w[i][j]);
      mylog(mylogf,ggstr);
    }
    
    for(j=0;j<yn;j++){
      for(k=0;k<xn;k++)
	if (flag[j][k])
	  if (rf[j][k] > 0.0)
	    mylog(mylogf," O");
	  else
	    mylog(mylogf," X");
	else
	  mylog(mylogf," .");
      mylog(mylogf,"\n");
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                       MODEL_DOG_COMP_01_GET_RESPONSE                      */
/*                                                                           */
/*  Calling sequence:                                                        */
/*    1.   model_dog_01_prep                                                 */
/*    2.   model_dog_01_prep_ft                                              */
/*    3.   model_dog_simp_01_prep                                            */
/*    4.   this routine                                                      */
/*        ...                                                                */
/*    N-1. model_dog_01_done                                                 */
/*    N.   model_dog_simp_01_done                                            */
/*                                                                           */
/*  NOTES:                                                                   */
/*    s->d, the stimulus data, is modified here.                             */
/*                                                                           */
/*****************************************************************************/
void model_dog_comp_01_get_response(m,s,r,k)
     struct model_struct *m;        /* Model parameter pair list. */
     struct stim_struct *s;         /* Stimulus data and param pair list */
     struct response_struct *r;     /* Response data */
     int k; /* Repeat number */
{
  int i,xi,yi;
  int xn,yn,tn,ns,vn,seed,tvn;
  int *cnts,**cnt0,**cnt1,*inx,*iny,nin,**flag,gnms;
  float ***rdog,*fs,sampling,amp,t1,t2,*v;
  float **ss,***s0,***s1,*inw,*gtx,*gti,sampms,samp;
  float ***tgx,***tgi,***tgv;
  char *outfile;

  sampms = 1000.0;

  xn = mod_dog_xn;
  yn = mod_dog_yn;
  tn = mod_dog_tn;

  // Next call modifies s->d
  rdog = mod_dog_get_response(s,k,NULL,0); // Set 'mod_dog_r' t'which rdog pnts

  outfile = paramfile_get_char_param_default(m->ppl,"write_3d_response","");
  if (strcmp(outfile,"")!=0)
    write_3d_data_part(outfile,rdog,1,xn,1,yn,1,tn,4,2,1);
  myfree(outfile);
  
  /*** Get spike trains for all relevant DOG grid positions ***/
  flag = get_zero_2d_iarray(xn,yn);
  for(i=0;i<mod_dog_simp_n;i++){
    for(xi=0;xi<xn;xi++)
      for(yi=0;yi<yn;yi++)
	if (mod_dog_simp_flag[i][xi][yi])
	  flag[xi][yi] = 1;
  }
  samp = 1.0/mod_dog_tscale;
  seed = m->mseed[r->tsi];
  mod_dog_get_spikes_2d_flag(m,r,rdog,NULL,1,xn,yn,tn,flag,samp,seed,
			     &s1,&cnt1,0,&tgx,&tgi,&tgv,&tvn); //on
  mod_dog_get_spikes_2d_flag(m,r,rdog,NULL,0,xn,yn,tn,flag,samp,seed,
			     &s0,&cnt0,0,&tgx,&tgi,&tgv,&tvn); //off

  /*** Use spike trains to get exc and inh conductances for simple cells ***/
  cnts = (int *)myalloc(mod_dog_simp_n*sizeof(int));
  ss = (float **)myalloc(mod_dog_simp_n*sizeof(float *));
  for(i=0;i<mod_dog_simp_n;i++){
    nin = mod_dog_simp_nin[i];
    inx = mod_dog_simp_x[i];
    iny = mod_dog_simp_y[i];
    inw = mod_dog_simp_w[i];
    mod_dog_get_gtx_gti(m,s0,s1,cnt0,cnt1,nin,inx,iny,inw,&gtx,&gti,&gnms,i);
    /*** Data returned in msec time units ***/

    /*** Get spike trains for this simple cell ***/
    ifc_test(myid,m,"ifc_simp",gtx,gti,sampms,gnms,m->mseed[r->tsi],1,
	     (char *)NULL,&ss[i],&cnts[i],0,&v,&vn);
    myfree(gti); myfree(gtx);
  }

  mylog(mylogf,"    Summing conductance from 'mod_dog_simp_n' simple cells.\n");
  t1   = paramfile_get_float_param_or_exit(m->ppl,"comp_ex_tau_f");
  t2   = paramfile_get_float_param_or_exit(m->ppl,"comp_ex_tau_r");
  amp  = paramfile_get_float_param_or_exit(m->ppl,"comp_ex_amp");
  mod_dog_get_gt_sarray(m,ss,cnts,mod_dog_simp_n,t1,t2,amp,&gtx,&gnms);

  mylog(mylogf,"    Computing complex cell spike train.\n");
  gti = get_zero_farray(gnms);
  ifc_test(myid,m,"ifc_comp",gtx,gti,sampms,gnms,m->mseed[r->tsi],1,
	   (char *)NULL,&fs,&ns,0,&v,&vn);
  myfree(gtx); myfree(gti);

  /* Spikes come in ms, adjust based on the 'save_spikes' sampling */
  /*sampling = paramfile_get_nth_float_param_or_exit(m->ppl,"save_spikes",1);*/
  sampling = mod_util_resp_get_samp(r,"spikes");
  multiply_farray(fs,ns,sampling/1000.0);

  /*** Store response, link to storage, copyflag = 0 ***/
  (void)mod_util_resp_check_store_s(r,"spikes",mylogf,fs,ns,0);
  
  /**if (r->ns > 0){
     r->s[0][r->tsi] = fs;
     r->cnt[0][r->tsi] = ns; * 'ti' is updated in 'run_gen_tuning_curve' *
     }**/
}
/**************************************-**************************************/
/*                                                                           */
/*                           MODEL_DOG_DS_GET_POOL_G                         */
/*                                                                           */
/*****************************************************************************/
void mod_dog_ds_get_pool_g(s1,s0,cnt1,cnt0,tf,ds_tmu,ds_tsd,t1,t2,amp,rg,rn)
     float ***s1,***s0;
     int **cnt1,**cnt0;
     float tf,ds_tmu,ds_tsd,t1,t2,amp;
     float **rg;
     int *rn;
{
  int i,j;
  int xi,yi,n,nx,t,tn,nms;
  static int trn = 0;
  float *ts,*xshape,*st0,*st1,*g0,*g1,*gt,*gtot,ampl,theta,period;
  char tstr[SLEN];

  nx = (int)(t1 * 7.0);        /* Length of mask is 7 times decay time */
  xshape = diff_exp_farray(0.0,amp,t1,t2,nx);

  tn = mod_dog_tn;
  nms = (int)((float)tn * 1000.0*mod_dog_tscale); /* Duration in ms */
  st0 = get_zero_farray(nms);
  st1 = get_zero_farray(nms);
  gt  = get_zero_farray(nms);
  gtot = get_zero_farray(nms);

  if (myid == -1){
    if (tf > 0.0){
      append_string_to_file("zz.ampthet.pl","/newplot\n");
      sprintf(tstr,"/plotname p_%d\n",trn);
      append_string_to_file("zz.ampthet.pl",tstr);
    }
  }
  n = mod_dog_ds_nsub;
  for(i=0;i<n;i++){ /*** For each subunit ***/
    xi = mod_dog_ds_subx1[i];
    yi = mod_dog_ds_suby1[i];

    /*** OFF units ***/
    ts = s0[xi][yi];
    for(j=0;j<cnt0[xi][yi];j++){ /* Expand spike train for DOG 0 */
      t = my_rint(ts[j]);
      if ((t >= 0)&&(t < nms)){
	st0[t] += 1.0;
      }
    }
    xi = mod_dog_ds_subx2[i];
    yi = mod_dog_ds_suby2[i];
    ts = s0[xi][yi];
    for(j=0;j<cnt0[xi][yi];j++){ /* Expand spike train for DOG 1 */
      t = my_rint(ts[j] + ds_tmu);
      if ((t >= 0)&&(t < nms)){
	st1[t] += 1.0;
      }
    }
    g0 = convolve_with_mask_causal(st0,nms,xshape,nx);
    g1 = convolve_with_mask_causal(st1,nms,xshape,nx);
    if (myid == -1){
      if (i==0){
	append_farray_plot("zz.subunit.pl","s0_subunit_0",st0,nms,1);
	append_farray_plot("zz.subunit.pl","s1_subunit_0",st1,nms,1);
	append_farray_plot("zz.subunit.pl","g0_subunit_0",g0,nms,1);
	append_farray_plot("zz.subunit.pl","g1_subunit_0",g1,nms,1);
      }
    }

    /*** PRODUCT ***/
    for(j=0;j<nms;j++){
      gt[j] = g0[j]*g1[j];
      gtot[j] += gt[j];
      st0[j] = st1[j] = 0.0;  /* Reset to zero */
    }
    myfree(g0); myfree(g1);
    if (myid == -1){
      if (i==1)
	append_farray_plot("zz.subunit.pl","Prod_subunit_0",gt,nms,1);
    }
    
    /*** Compute F1 ampl and phase ***/
    if (myid == -1){
      if (tf > 0.0){
	period = 1000.0/tf;
	get_fourier_harmonic_farray(gt,nms,1,period,&ampl,&theta);
	sprintf(tstr,"%.4f %.4f\n",ampl,theta);
	append_string_to_file("zz.ampthet.pl",tstr);
      }
    }


    /*** ON units ***/
    ts = s1[xi][yi];
    for(j=0;j<cnt1[xi][yi];j++){ /* Expand spike train for DOG 0 */
      t = my_rint(ts[j]);
      if ((t >= 0)&&(t < nms)){
	st0[t] += 1.0;
      }
    }
    xi = mod_dog_ds_subx2[i];
    yi = mod_dog_ds_suby2[i];
    ts = s1[xi][yi];
    for(j=0;j<cnt1[xi][yi];j++){ /* Expand spike train for DOG 1 */
      t = my_rint(ts[j] + ds_tmu);
      if ((t >= 0)&&(t < nms)){
	st1[t] += 1.0;
      }
    }
    g0 = convolve_with_mask_causal(st0,nms,xshape,nx);
    g1 = convolve_with_mask_causal(st1,nms,xshape,nx);
    if (myid == -1){
      if (i==0){
	append_farray_plot("zz.subunit.pl","s0_subunit_0",st0,nms,1);
	append_farray_plot("zz.subunit.pl","s1_subunit_0",st1,nms,1);
	append_farray_plot("zz.subunit.pl","g0_subunit_0",g0,nms,1);
	append_farray_plot("zz.subunit.pl","g1_subunit_0",g1,nms,1);
      }
    }

    /*** PRODUCT ***/
    for(j=0;j<nms;j++){
      gt[j] = g0[j]*g1[j];
      gtot[j] += gt[j];
      st0[j] = st1[j] = 0.0;  /* Reset to zero */
    }
    myfree(g0); myfree(g1);
    if (myid == -1){
      if (i==1)
	append_farray_plot("zz.subunit.pl","Prod_subunit_0",gt,nms,1);
    }
    /*** Compute F1 ampl and phase ***/
    if (myid == -1){
      if (tf > 0.0){
	period = 1000.0/tf;
	get_fourier_harmonic_farray(gt,nms,1,period,&ampl,&theta);
	sprintf(tstr,"%.4f %.4f\n",ampl,theta);
	append_string_to_file("zz.ampthet.pl",tstr);
      }
    }
  }
  myfree(st0); myfree(st1); myfree(gt); myfree(xshape);

  /*append_farray_plot("zz.subunit.pl","Tot_of_Prods",gtot,nms,1);*/
  gt = smooth_with_gaussian(gtot,nms,1.0,0.01);
  myfree(gtot);
  gtot = gt;
  if (myid == -1)
    append_farray_plot("zz.subunit.pl","Tot_of_Prods__smoothsig=1",gtot,nms,1);

  trn += 1;
  *rg = gtot; *rn = nms;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MODEL_DOG_DS_01_PREP                           */
/*                                                                           */
/*  Pick pairs of coordinates to be used for DS subunits.                    */
/*                                                                           */
/*****************************************************************************/
void model_dog_ds_01_prep(mppl,r)
     struct param_pair_list *mppl; /* Model parameter pair list. */
     struct response_struct *r;    /* Response params. */
{
  int i;
  int xn,yn,seed,gseed,n,*sx,*sy,sn,done;
  float ssd,xmu,xsd,ymu,ysd,xmuo2,ymuo2,xir,yir,rfs,eps,wsd;
  char tstr[SLEN];

  mylog(mylogf,"  MODEL_DOG_DS_01_PREP\n");

  /*sscale = mod_dog_sscale;*/
  xn = mod_dog_xn;
  yn = mod_dog_yn;

  ssd   = paramfile_get_float_param_or_exit(mppl,"ds_ssd");
  xmu = paramfile_get_float_param_or_exit(mppl,"ds_xmu");
  xsd = paramfile_get_float_param_or_exit(mppl,"ds_xsd");
  ymu = paramfile_get_float_param_or_exit(mppl,"ds_ymu");
  ysd = paramfile_get_float_param_or_exit(mppl,"ds_ysd");
  seed  = paramfile_get_int_param_or_exit(mppl,"ds_seed");
  if (seed > 0)
    seed = -seed;
  gseed = seed;

  /* Number of DS subunits */
  n = paramfile_get_int_param_or_exit(mppl,"ds_nsub");

  /***  rfs = 20.0   wsd =  2.0    eps =  0.2 ***/
  rfs = paramfile_get_float_param_or_exit(mppl,"ds_rfs");
  wsd = paramfile_get_float_param_or_exit(mppl,"ds_wsd");
  eps = paramfile_get_float_param_or_exit(mppl,"ds_eps");
  model_dog_get_sampling_points(xn,yn,n,seed,ssd,rfs,wsd,eps,&sx,&sy,&sn);
  n = sn;
  mod_dog_ds_nsub = n;

  /*** Pick 'n' points in the Gaussian RF profile ***/
  mod_dog_ds_subx1 = (int *)myalloc(n*sizeof(int));
  mod_dog_ds_suby1 = (int *)myalloc(n*sizeof(int));
  mod_dog_ds_subx2 = (int *)myalloc(n*sizeof(int));
  mod_dog_ds_suby2 = (int *)myalloc(n*sizeof(int));
  xmuo2 = xmu/2.0;
  ymuo2 = ymu/2.0;

  for(i=0;i<n;i++){
    xir = (float)sx[i] + 0.5*nr_util_gasdev(&gseed);
    yir = (float)sy[i] + 0.5*nr_util_gasdev(&gseed);
    done = 0;
    while(!done){
      mod_dog_ds_subx1[i] = my_rint(xir - xmuo2 + xsd*nr_util_gasdev(&gseed));
      mod_dog_ds_suby1[i] = my_rint(yir - ymuo2 + ysd*nr_util_gasdev(&gseed));
      mod_dog_ds_subx2[i] = my_rint(xir + xmuo2 + xsd*nr_util_gasdev(&gseed));
      mod_dog_ds_suby2[i] = my_rint(yir + ymuo2 + ysd*nr_util_gasdev(&gseed));
      if ((mod_dog_ds_subx1[i] >= 0)&&(mod_dog_ds_subx1[i] < xn)&&
	  (mod_dog_ds_subx2[i] >= 0)&&(mod_dog_ds_subx2[i] < xn)&&
	  (mod_dog_ds_suby1[i] >= 0)&&(mod_dog_ds_suby1[i] < yn)&&
	  (mod_dog_ds_suby2[i] >= 0)&&(mod_dog_ds_suby2[i] < yn)){
	done = 1;
      }
    }
  }

  if (myid == -1){
    mylog(mylogf,"    Writing x,y coords for subunits to:  zzz.DS.RF.pl\n");
    remove_file("zzz.DS.RF.pl");
    append_string_to_file("zzz.DS.RF.pl","/newplot\n");
    sprintf(tstr,"0 0\n0 %d\n%d %d\n%d 0\n0 0\n",yn-1,xn-1,yn-1,xn-1);
    append_string_to_file("zzz.DS.RF.pl",tstr);
    for(i=0;i<n;i++){
      append_string_to_file("zzz.DS.RF.pl","/newplot\n");
      sprintf(tstr,"%4d %4d\n%4d %4d\n",mod_dog_ds_subx1[i],
	      mod_dog_ds_suby1[i],mod_dog_ds_subx2[i],mod_dog_ds_suby2[i]);
      append_string_to_file("zzz.DS.RF.pl",tstr);
    }
  }

  myfree(sx); myfree(sy);
}
/**************************************-**************************************/
/*                                                                           */
/*                          MODEL_DOG_DS_01_GET_RESPONSE                     */
/*                                                                           */
/*  Calling sequence:                                                        */
/*    1.   model_dog_01_prep                                                 */
/*    2.   model_dog_01_prep_ft                                              */
/*    3.   model_dog_simp_01_prep                                            */
/*    4.   this routine                                                      */
/*        ...                                                                */
/*    N-1. model_dog_01_done                                                 */
/*    N.   model_dog_simp_01_done                                            */
/*                                                                           */
/*  NOTES:                                                                   */
/*    s->d, the stimulus data, is modified here.                             */
/*                                                                           */
/*****************************************************************************/
void model_dog_ds_01_get_response(m,s,r,k)
     struct model_struct *m;        /* Model parameter pair list. */
     struct stim_struct *s;         /* Stimulus data and param pair list */
     struct response_struct *r;     /* Response data */
     int k;                         /* Repeat number */
{
  int i,j;
  int xn,yn,tn,ns,**cnt0,**cnt1,**flag,gnms,vn,seed,tvn;
  float ***rdog,*fs,sampling,amp,t1,t2,ds_tmu,ds_tsd,*v;
  float ***s0,***s1,*gtx,*gti,sampms,tf,samp;
  float ***tgx,***tgi,***tgv;
  char *outfile;

  xn = mod_dog_xn;
  yn = mod_dog_yn;
  tn = mod_dog_tn;

  if (paramfile_test_param(s->ppl,"tf"))
    tf = paramfile_get_float_param_or_exit(s->ppl,"tf","");
  else
    tf = -1.0;

  // Next call modifies s->d
  rdog = mod_dog_get_response(s,k,NULL,0); // Set 'mod_dog_r' t'which rdog pnts

  outfile = paramfile_get_char_param_default(m->ppl,"write_3d_response","");
  if (strcmp(outfile,"")!=0)
    write_3d_data_part(outfile,rdog,1,xn,1,yn,1,tn,4,2,1);
  myfree(outfile);

  /*** Get spike trains for all relevant DOG grid positions ***/
  flag = get_zero_2d_iarray(xn,yn);
  for(i=0;i<mod_dog_ds_nsub;i++){ /* For each subunit */
    flag[mod_dog_ds_subx1[i]][mod_dog_ds_suby1[i]] = 1;
    flag[mod_dog_ds_subx2[i]][mod_dog_ds_suby2[i]] = 1;
  }
  /*** Print out the 2D 'flag' array ***/
  for(i=0;i<yn;i++){
    for(j=0;j<xn;j++){
      if (flag[j][i])
	mylog(mylogf," *");
      else
	mylog(mylogf," .");
    }
    mylog(mylogf,"\n");
  }
  samp = 1.0/mod_dog_tscale;
  seed = m->mseed[r->tsi];
  mod_dog_get_spikes_2d_flag(m,r,rdog,NULL,1,xn,yn,tn,flag,samp,seed,
			     &s1,&cnt1,0,&tgx,&tgi,&tgv,&tvn);//on
  mod_dog_get_spikes_2d_flag(m,r,rdog,NULL,0,xn,yn,tn,flag,samp,seed,
			     &s0,&cnt0,0,&tgx,&tgi,&tgv,&tvn);//ff

  /*** Use paired spike trains to get excitatory conductance ***/
  t1  = paramfile_get_float_param_or_exit(m->ppl,"ds_ex_tau_f");
  t2  = paramfile_get_float_param_or_exit(m->ppl,"ds_ex_tau_r");
  amp = paramfile_get_float_param_or_exit(m->ppl,"ds_ex_amp");
  ds_tmu  = paramfile_get_float_param_or_exit(m->ppl,"ds_tmu");
  ds_tsd  = paramfile_get_float_param_or_exit(m->ppl,"ds_tsd");
  mod_dog_ds_get_pool_g(s1,s0,cnt1,cnt0,tf,ds_tmu,ds_tsd,t1,t2,amp,&gtx,&gnms);

  mylog(mylogf,"    Computing complex cell spike train.\n");
  gti = get_zero_farray(gnms);
  sampms = 1000.0;
  ifc_test(myid,m,"ifc_comp",gtx,gti,sampms,gnms,m->mseed[r->tsi],1,
	   (char *)NULL,&fs,&ns,0,&v,&vn);
  myfree(gtx); myfree(gti);

  /* Spikes come in ms, adjust based on the 'save_spikes' sampling */
  /*sampling = paramfile_get_nth_float_param_or_exit(m->ppl,"save_spikes",1);*/
  sampling = mod_util_resp_get_samp(r,"spikes");
  multiply_farray(fs,ns,sampling/1000.0);

  /*** Store response, link to storage, copyflag = 0 ***/
  (void)mod_util_resp_check_store_s(r,"spikes",mylogf,fs,ns,0);
  
  /**if (r->ns > 0){
     r->s[0][r->tsi] = fs;
     r->cnt[0][r->tsi] = ns; * 'ti' is updated in 'run_gen_tuning_curve' *
     }**/
}
/**************************************-**************************************/
/*                                                                           */
/*                             MODEL_DOG_DS_02_PREP                          */
/*                                                                           */
/*  Pick simple cell params to be used in DS computation.                    */
/*                                                                           */
/*****************************************************************************/
void model_dog_ds_02_prep(mppl,r)
     struct param_pair_list *mppl; /* Model parameter pair list. */
     struct response_struct *r;    /* Response params. */
{
  int i,j,k;
  int xn,yn,seed,n;
  int **flag,nsimp,nsamp,xi,yi;
  float **rf,sscale,theta,ssd,sf,ph,min,max,fsamp;
  float tph;

  mylog(mylogf,"  MODEL_DOG_DS_02_PREP\n");

  sscale = mod_dog_sscale;
  xn = mod_dog_xn;
  yn = mod_dog_yn;

  /* Number of DS subunits */
  n = paramfile_get_int_param_or_exit(mppl,"ds_nsub");
  nsimp = 2*n;

  ssd   = paramfile_get_float_param_or_exit(mppl,"simp_ssd");
  sf    = paramfile_get_float_param_or_exit(mppl,"simp_sf");
  ph    = paramfile_get_float_param_or_exit(mppl,"simp_phase");
  theta = paramfile_get_float_param_or_exit(mppl,"simp_ori");
  fsamp = paramfile_get_float_param_or_exit(mppl,"simp_rf_fsamp");
  seed  = paramfile_get_int_param_or_exit(mppl,"simp_rf_seed");
  nsamp = paramfile_get_int_param_or_exit(mppl,"simp_rf_nsamp");
  if (seed > 0)
    seed = -seed;

  mod_dog_simp_n = nsimp;
  mod_dog_simp_nin = get_const_iarray(nsimp,nsamp);
  mod_dog_simp_x = get_2d_iarray(nsimp,nsamp);
  mod_dog_simp_y = get_2d_iarray(nsimp,nsamp);
  mod_dog_simp_w = get_2d_farray(nsimp,nsamp);
  mod_dog_simp_rf = (float ***)myalloc(nsimp*sizeof(float **));
  mod_dog_simp_flag = get_zero_3d_iarray(nsimp,xn,yn);


  tph = ph;
  for(i=0;i<nsimp;i++){
    if (i%2 == 0){
      rf = gabor_2d_space(xn,yn,sscale,ssd,sf,theta,tph);
      
      /*** WYETH - randomly change phase ***/
    }else{
      rf = gabor_2d_space(xn,yn,sscale,ssd,sf,theta,tph+90.0);
      tph += 90.0; /* Next pair will have 90 degr shift */
    }
    mod_dog_simp_rf[i] = rf;
    
    max = max_of_2d_farray(rf,xn,yn);
    min = min_of_2d_farray(rf,xn,yn);
    sprintf(ggstr,"  min, max =  %f %f\n",min,max);
    mylog(mylogf,ggstr);
    
    /*** Select inputs from 2D spatial array of DOG responses ***/
    flag = mod_dog_simp_flag[i];
    k = 0;
    while(k < nsamp){
      xi = (int)((float)xn * myrand_util_ran2(&seed));
      yi = (int)((float)yn * myrand_util_ran2(&seed));
      if (flag[xi][yi] == 0){
	if ((rf[xi][yi] > fsamp*max)|| (rf[xi][yi] < fsamp*min)){
	  flag[xi][yi] = 1;
	  mod_dog_simp_x[i][k] = xi;
	  mod_dog_simp_y[i][k] = yi;
	  mod_dog_simp_w[i][k] = rf[xi][yi];
	  k += 1;
	}
      }
    }
    
    for(j=0;j<mod_dog_simp_nin[i];j++){
      sprintf(ggstr,"%d (%d,%d) %.2f\n",j,mod_dog_simp_x[i][j],
	      mod_dog_simp_y[i][j],mod_dog_simp_w[i][j]);
      mylog(mylogf,ggstr);
    }
    
    for(j=0;j<yn;j++){
      for(k=0;k<xn;k++)
	if (flag[j][k])
	  if (rf[j][k] > 0.0)
	    mylog(mylogf," O");
	  else
	    mylog(mylogf," X");
	else
	  mylog(mylogf," .");
      mylog(mylogf,"\n");
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_DOG_GET_GT_DS_02                          */
/*                                                                           */
/*  Combine pairs of simple cell spike trains to compute total conductance   */
/*  to DS unit.                                                              */
/*                                                                           */
/*****************************************************************************/
void mod_dog_get_gt_ds_02(m,s,cnt,n,t1,t2,amp,delay,rgt,rnms)
     struct model_struct *m;   /* Model parameter pair list. */
     float **s;
     int *cnt,n;
     float t1,t2,amp;
     int *delay;               /* [n/2] Delay for each pair of spike trains */
     float **rgt;              /* Conductance vs. time */
     int *rnms;                /* Return duration in msec */
{
  int i,j;
  int nx,tn,nms;
  float *xshape,*gt,*g,xarea;

  /*** Determine duration in msec ***/
  tn = mod_dog_tn;
  nms = (int)((float)tn * 1000.0*mod_dog_tscale); /* Duration in ms */

  /*** For each consecutive pair of spikes, delay the second and multiply ***/
  gt = get_zero_farray(nms);

  nx = (int)(t1 * 7.0);        /* Length of mask is 7 times decay time */
  xshape = diff_exp_farray(0.0,amp,t1,t2,nx);
  xarea = (float)nx * mean_farray(xshape,nx);
  if (myid == -1)
    append_farray_plot("zz.xshape.pl","xshape_0",xshape,nx,1);
  for(i=0;i<n;i+=2){
    mod_dog_get_gt_prod_pair_sarray(s[i],cnt[i],s[i+1],cnt[i+1],nms,
				    delay[i/2],xshape,nx,&g);

    for(j=0;j<nms;j++)
      gt[j] += g[j];
    myfree(g);
  }
  myfree(xshape);

  t1 *= 1.414;
  t2 *= 1.414;
  nx = (int)(t1 * 7.0);        /* Length of mask is 7 times decay time */
  xshape = diff_exp_farray(0.0,amp,t1,t2,nx);
  norm_area_farray(xshape,nx,xarea*1.0);  /*** WYETH lower weight ***/
  if (myid == -1)
    append_farray_plot("zz.xshape.pl","xshape_0",xshape,nx,1);
  for(i=0;i<n;i+=2){
    mod_dog_get_gt_prod_pair_sarray(s[i],cnt[i],s[i+1],cnt[i+1],nms,
				    delay[i/2]*2,xshape,nx,&g);
    for(j=0;j<nms;j++)
      gt[j] += g[j];
    myfree(g);
  }
  myfree(xshape);

  t1 *= 1.414;
  t2 *= 1.414;
  nx = (int)(t1 * 7.0);        /* Length of mask is 7 times decay time */
  xshape = diff_exp_farray(0.0,amp,t1,t2,nx);
  norm_area_farray(xshape,nx,xarea*1.0); /*** WYETH - lower weight ***/
  if (myid == -1)
    append_farray_plot("zz.xshape.pl","xshape_0",xshape,nx,1);
  for(i=0;i<n;i+=2){
    mod_dog_get_gt_prod_pair_sarray(s[i],cnt[i],s[i+1],cnt[i+1],nms,
				    delay[i/2]*4,xshape,nx,&g);
    for(j=0;j<nms;j++)
      gt[j] += g[j];
    myfree(g);
  }
  myfree(xshape);

  t1 *= 2.0;
  t2 *= 2.0;
  nx = (int)(t1 * 7.0);        /* Length of mask is 7 times decay time */
  xshape = diff_exp_farray(0.0,amp,t1,t2,nx);
  norm_area_farray(xshape,nx,xarea*0.2); /*** WYETH - more area ***/
  if (myid == -1)
    append_farray_plot("zz.xshape.pl","xshape_0",xshape,nx,1);
  for(i=0;i<n;i+=2){
    mod_dog_get_gt_prod_pair_sarray(s[i],cnt[i],s[i+1],cnt[i+1],nms,
				    delay[i/2]*8,xshape,nx,&g);
    for(j=0;j<nms;j++)
      gt[j] += g[j];
    myfree(g);
  }
  myfree(xshape);

  if (myid == -1)
    append_farray_plot("zz.ds02.gtot.out.pl","gtx_ds_tot",gt,nms,1);

  *rgt = gt; *rnms = nms;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MODEL_DOG_DS_02_GET_RESPONSE                     */
/*                                                                           */
/*  Calling sequence:                                                        */
/*    1.   model_dog_01_prep                                                 */
/*    2.   model_dog_01_prep_ft                                              */
/*    3.   model_dog_ds_02_prep                                              */
/*    4.   this routine                                                      */
/*        ...                                                                */
/*    N-1. model_dog_01_done                                                 */
/*    N.   model_dog_ds_02_done                                              */
/*                                                                           */
/*  NOTES:                                                                   */
/*    s->d, the stimulus data, is modified here.                             */
/*                                                                           */
/*****************************************************************************/
void model_dog_ds_02_get_response(m,s,r,k)
     struct model_struct *m;        // Model parameter pair list
     struct stim_struct *s;         // Stimulus data and param pair list
     struct response_struct *r;     // Response data
     int k; // Repeat number
{
  int i,xi,yi;
  int xn,yn,tn,ns,vn,seed,tvn;
  int *cnts,**cnt0,**cnt1,*inx,*iny,nin,**flag,gnms,*delay,ds_delay;
  int meshflag,sc_simp;
  float ***rdog,*fs,sampling,amp,t1,t2,*v,fscale;
  float **ss,***s0,***s1,*inw,*gtx,*gti,samp,sampms;
  float ***tgx,***tgi,***tgv;
  char *outfile;

  sampms = 1000.0;

  xn = mod_dog_xn;
  yn = mod_dog_yn;
  tn = mod_dog_tn;

  meshflag = paramfile_get_int_param_default(m->ppl,"mesh_flag",0);
  if (meshflag){
    fscale = 22.0;
    sprintf(ggstr,"      Scaling 'rdog' by %f, to match non-mesh ds02.\n",
	    fscale);
    mylog(mylogf,ggstr);
    printf("  ***  See .../libc/archive/mod_mesh_util.* for old code\n");
    exit_error("MODEL_DOG_DS_02_GET_RESPONSE","Old way no longer supported");
    //rdog = model_mesh_rgc_01_get_3d_response(m,s,r,k,fscale);
  }else{
    // Next call modifies s->d
    rdog = mod_dog_get_response(s,k,NULL,0); // Set mod_dog_r t'which rdog pnts
  }


  /*** CHECK MIN and MAX of RESPONSE ***/
  { 
    float min,max;

    get_min_max_3d_farray(rdog,1,xn,1,yn,1,tn,&min,&max);
    sprintf(ggstr,"  3D retina response min, max =  %f %f\n",min,max);
    mylog(mylogf,ggstr);

    if (myid == -1)
      append_farray_plot("zz.rdog_mid_v_t.pl","rdog_mid_v_t",
			 rdog[xn/2][yn/2]+1,tn,1);
  }

  outfile = paramfile_get_char_param_default(m->ppl,"write_3d_response","");
  if (strcmp(outfile,"")!=0)
    write_3d_data_part(outfile,rdog,1,xn,1,yn,1,tn,4,2,1);
  myfree(outfile);

  /*** Get spike trains for all relevant DOG grid positions ***/
  flag = get_zero_2d_iarray(xn,yn);
  for(i=0;i<mod_dog_simp_n;i++){
    for(xi=0;xi<xn;xi++)
      for(yi=0;yi<yn;yi++)
	if (mod_dog_simp_flag[i][xi][yi])
	  flag[xi][yi] = 1;
  }
  samp = 1.0/mod_dog_tscale;

  seed = m->mseed[r->tsi];
  mod_dog_get_spikes_2d_flag(m,r,rdog,NULL,1,xn,yn,tn,flag,samp,seed,
			     &s1,&cnt1,0,&tgx,&tgi,&tgv,&tvn);//on
  mod_dog_get_spikes_2d_flag(m,r,rdog,NULL,0,xn,yn,tn,flag,samp,seed,
			     &s0,&cnt0,0,&tgx,&tgi,&tgv,&tvn);//ff

  /*** Use spike trains to get exc and inh conductances for simple cells ***/
  cnts = (int *)myalloc(mod_dog_simp_n*sizeof(int));
  ss = (float **)myalloc(mod_dog_simp_n*sizeof(float *));
  sc_simp = 0; /* Spike count for simple cells */
  for(i=0;i<mod_dog_simp_n;i++){
    nin = mod_dog_simp_nin[i];
    inx = mod_dog_simp_x[i];
    iny = mod_dog_simp_y[i];
    inw = mod_dog_simp_w[i];
    mod_dog_get_gtx_gti(m,s0,s1,cnt0,cnt1,nin,inx,iny,inw,&gtx,&gti,&gnms,i);
    /*** Data returned in msec time units ***/

    /*** Get spike trains for this simple cell ***/
    ifc_test(myid,m,"ifc_simp",gtx,gti,sampms,gnms,m->mseed[r->tsi],0,
	     (char *)NULL,&ss[i],&cnts[i],0,&v,&vn);
    sc_simp += cnts[i];
    myfree(gti); myfree(gtx);
  }
  mylog(mylogf,"  Simple cells\n");
  sprintf(ggstr,"    %d cells, %.4f spikes/cell\n",mod_dog_simp_n,
	  (float)sc_simp/(float)mod_dog_simp_n);
  mylog(mylogf,ggstr);

  /*** Sum DS conductances from paired simple cell spike trains ***/
  t1   = paramfile_get_float_param_or_exit(m->ppl,"comp_ex_tau_f");
  t2   = paramfile_get_float_param_or_exit(m->ppl,"comp_ex_tau_r");
  amp  = paramfile_get_float_param_or_exit(m->ppl,"comp_ex_amp");
  ds_delay  = paramfile_get_int_param_or_exit(m->ppl,"ds_delay");

  delay = get_zero_iarray(mod_dog_simp_n/2);
  for(i=0;i<mod_dog_simp_n/2;i++)
    delay[i] = ds_delay;
  mod_dog_get_gt_ds_02(m,ss,cnts,mod_dog_simp_n,t1,t2,amp,delay,&gtx,&gnms);
  myfree(delay);

  sprintf(ggstr,"  Complex DS cell (%d DS subunits)\n",mod_dog_simp_n/2);
  mylog(mylogf,ggstr);
  gti = get_zero_farray(gnms);

  /*** WYETH - this smoothing in attempt to avoid 'kmax exceeded' ***/
  /*** But it didn't seem to help ***/
  if (myid == -1)
    append_farray_plot("zz.ifc_comp.pl","gtx",gtx,gnms,1);

  {
    float *gsm;
    gsm = smooth_with_gaussian(gtx,gnms,1.0,0.01);
    myfree(gtx);
    gtx = gsm;
    if (myid == -1)
      append_farray_plot("zz.ifc_comp.pl","gtx_smooth",gtx,gnms,1);
  }

  ifc_test(myid,m,"ifc_comp",gtx,gti,sampms,gnms,m->mseed[r->tsi],1,
	   (char *)NULL,&fs,&ns,0,&v,&vn);
  myfree(gtx); myfree(gti);
  
  // Spikes come in ms, adjust based on the 'save_spikes' sampling
  //sampling = paramfile_get_nth_float_param_or_exit(m->ppl,"save_spikes",1);
  sampling = mod_util_resp_get_samp(r,"spikes");
  multiply_farray(fs,ns,sampling/1000.0);

  /*** Store response, link to storage, copyflag = 0 ***/
  (void)mod_util_resp_check_store_s(r,"spikes",mylogf,fs,ns,0);
  
  /**if (r->ns > 0){
     r->s[0][r->tsi] = fs;
     r->cnt[0][r->tsi] = ns; * 'ti' is updated in 'run_gen_tuning_curve' *
     }**/
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_DOG_RUN_DOG_01                            */
/*                                                                           */
/*****************************************************************************/
void mod_dog_run_dog_01(m,s,r,k,action)
     struct model_struct *m;    // Model params.
     struct stim_struct *s;     // Stimulus params.
     struct response_struct *r; // Response params.
     int k;                     // Trial number, when action = 1
     int action;                // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  if (action == 0){

    if (m->ppl == NULL){
      model_dog_01_prep_o(m,r);
    }else{
      model_dog_01_prep(m,r);
    }
    model_dog_01_prep_ft();
  }else if (action == 1){
    model_dog_01_get_response(m,s,r,k);
  }else if (action == -1){
    model_dog_01_done(m);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_DOG_RUN_SIMP_01                           */
/*                                                                           */
/*****************************************************************************/
void mod_dog_run_simp_01(m,s,r,k,action)
     struct model_struct *m;    // Model params.
     struct stim_struct *s;     // Stimulus params.
     struct response_struct *r; // Response params.
     int k;                     // Trial number, when action = 1
     int action;                // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // DO THIS FIRST

  if (action == 0){
    model_dog_01_prep(m,r);
    model_dog_01_prep_ft();
    model_dog_simp_01_prep(m->ppl,r);
  }else if (action == 1){
    model_dog_simp_01_get_response(m,s,r,k);
  }else if (action == -1){
    model_dog_01_done(m);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_DOG_RUN_COMP_01                           */
/*                                                                           */
/*****************************************************************************/
void mod_dog_run_comp_01(m,s,r,k,action)
     struct model_struct *m;    /* Model params. */
     struct stim_struct *s;     /* Stimulus params. */
     struct response_struct *r; /* Response params. */
     int k;                     /* Trial number, when action = 1 */
     int action;                /* -1-cleanup, 0-prep, 1-run */
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  /*** DO THIS FIRST ***/

  if (action == 0){
    model_dog_01_prep(m,r);
    model_dog_01_prep_ft();
    model_dog_comp_01_prep(m->ppl,r);
  }else if (action == 1){
    model_dog_comp_01_get_response(m,s,r,k);
  }else if (action == -1){
    model_dog_01_done(m);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_DOG_RUN_DS_01                             */
/*                                                                           */
/*****************************************************************************/
void mod_dog_run_ds_01(m,s,r,k,action)
     struct model_struct *m;    /* Model params. */
     struct stim_struct *s;     /* Stimulus params. */
     struct response_struct *r; /* Response params. */
     int k;                     /* Trial number, when action = 1 */
     int action;                /* -1-cleanup, 0-prep, 1-run */
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  /*** DO THIS FIRST ***/

  if (action == 0){
    model_dog_01_prep(m,r);
    model_dog_01_prep_ft();
    model_dog_ds_01_prep(m->ppl,r);
  }else if (action == 1){
    model_dog_ds_01_get_response(m,s,r,k);
  }else if (action == -1){
    model_dog_01_done(m);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_DOG_RUN_DS_02                             */
/*                                                                           */
/*****************************************************************************/
void mod_dog_run_ds_02(m,s,r,k,action)
     struct model_struct *m;    /* Model params. */
     struct stim_struct *s;     /* Stimulus params. */
     struct response_struct *r; /* Response params. */
     int k;                     /* Trial number, when action = 1 */
     int action;                /* -1-cleanup, 0-prep, 1-run */
{
  int meshflag;

  myid = mylog_set_id_log(m->process_id,&mylogf);  /*** DO THIS FIRST ***/

  meshflag = paramfile_get_int_param_default(m->ppl,"mesh_flag",0);

  if (action == 0){
    if (meshflag)
      model_mesh_rgc_01_prep(m,s,r);
    model_dog_01_prep(m,r);
    model_dog_01_prep_ft();
    model_dog_ds_02_prep(m->ppl,r);
  }else if (action == 1){
    model_dog_ds_02_get_response(m,s,r,k);
  }else if (action == -1){
    model_dog_01_done(m);
    if (meshflag)
      model_mesh_rgc_01_done(m->ppl);
  }
}
