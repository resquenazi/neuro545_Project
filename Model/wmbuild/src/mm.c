/*****************************************************************************/
/*                                                                           */
/*  mm.c                                                                     */
/*  wyeth bair                                                               */
/*                                                                           */
/*  NOTES                                                                    */
/*  - On slaves, files get written to directory containing executable.       */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************

  mpicc -O2 -o ../bin/mm mm.c -L../lib -I../libc \
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
//  -L/usr/X11/lib -I/usr/X11/include \
//  -lmod_gui_util -lmod_mon_util -lmod_gr_util 
//  -lwygit_util -lglplot_util 
//  -lGL -lX11

#include "mpi.h"
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
//
//
#include "paramfile_util.h"
#include "ndata_util.h"
#include "mctrl_util.h"
//
#include "mod_util.h"
#include "wm_util.h"
#include "paramfile.h"
#include "ndata.h"
#include "mod.h"
#include "mm_util.h"
#include "mm.h"

#define BUFFSIZE  100000
#define MAX_BUF_RESULT 1000000  // Size of buffer to received float result

// Global state for MPI (static to avoid name confict arising from init val
static int myid = -1;        // Process ID, -1-Single, 0-Master, >0-Slave
static char *mylogf = NULL;  // Log file name
static char *tlogf = NULL;   // ja file, job assignment file
char ggstr[LONG_SLEN];       // Global log string

int numprocs;        // Total number of processes, INCLUDING master
char processor_name[MPI_MAX_PROCESSOR_NAME];

int mm_gui_flag = 1;         // 0-no gui windows (to run in background)

/**************************************-**************************************/
/*                                                                           */
/*                              MM_MAKE_VAR_NAME                             */
/*                                                                           */
/*****************************************************************************/
char *mm_make_var_name(s,ti)
     struct stim_struct *s;      // Stimulus params
     int ti;                     // Trial index
{
  int i;
  char tstr[LONG_SLEN];

  if (s->nvar == 0){
    sprintf(tstr,"No parameters vary.");
  }else{
    tstr[0] = '\0';
    for(i=0;i<s->nvar;i++){
      strcat(tstr,s->vname[i]);
      strcat(tstr,"=");
      strcat(tstr,s->vval[ti][i]);
      strcat(tstr,"  ");
    }
  }
  return strdup(tstr);
}

/**************************************-**************************************/
/*                                                                           */
/*                            MM_SLAVE_WAIT_TO_SEND                          */
/*                                                                           */
/*****************************************************************************/
void mm_slave_wait_to_send(r,gti)
     struct response_struct *r; // Response params
     int gti;                   // Global trial index (over stim & models)
{
  int i;
  int done,tsig;
  MPI_Status status;

  //
  //  Tell the master that we want to send a result
  //
  MPI_Send(&gti,1,MPI_INT,0,tag_req2send,MPI_COMM_WORLD); // Request to send
  //mylog(mylogf," ==> Want to send results, waiting for OK to send ...\n");

  //
  //  Wait to receive 1 integer from proc 0
  //
  MPI_Recv(&tsig,1,MPI_INT,0,tag_ok2send,MPI_COMM_WORLD,&status);
  //mylog(mylogf,"     Received OK to send.\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                           MM_SLAVE_SEND_RESULTS                           */
/*                                                                           */
/*****************************************************************************/
void mm_slave_send_results(r,gti)
     struct response_struct *r; // Response params
     int gti;                   // Global trial index (over stim & models)
{
  int i;
  int n,ns,nf;
  float *pf;

  ns = r->ns;
  nf = r->nf;

  // First check if all results have been stored, exit if not
  wm_response_check_trial_full(r,gti,mylogf,processor_name);

  //
  //  WYETH NEW 2013 Nov
  //
  mm_slave_wait_to_send(r,gti);

  //sprintf(ggstr," *** Sending %d spike trains (%d float)\n",ns,nf);
  //mylog(mylogf,ggstr);

  MPI_Send(&gti,1,MPI_INT,0,tag_result,MPI_COMM_WORLD);  // Global trial index
  MPI_Send(&ns,1,MPI_INT,0,0,MPI_COMM_WORLD);  // No. of spike responses
  for(i=0;i<ns;i++){
    pf = r->s[i][gti];
    n  = r->cnt[i][gti];

    if (n > MAX_BUF_RESULT){
      sprintf(ggstr,"gti %d  i %d  n %d  MAX %d\n",gti,i,n,MAX_BUF_RESULT);
      mylog(mylogf,ggstr);
      mylog_exit(mylogf,"MM_SLAVE_SEND_RESULTS  cnt > MAX_BUF_RESULT\n");
    }else if (n == -1){
      sprintf(ggstr,"gti %d  i %d  n %d  %s\n",gti,i,n,r->datid[gti]);
      mylog(mylogf,ggstr);
      mylog_exit(mylogf,"MM_SLAVE_SEND_RESULTS  cnt > MAX_BUF_RESULT\n");
    }

    MPI_Send(pf,n,MPI_FLOAT,0,0,MPI_COMM_WORLD);
  }
  MPI_Send(&nf,1,MPI_INT,0,0,MPI_COMM_WORLD);  // No. of float responses
  for(i=0;i<nf;i++){
    pf = r->f[i][gti];
    n  = r->fcnt[i][gti];
    if (n > MAX_BUF_RESULT)
      mylog_exit(mylogf,"MM_SLAVE_SEND_RESULTS  fcnt > MAX_BUF_RESULT\n");
    MPI_Send(pf,n,MPI_FLOAT,0,0,MPI_COMM_WORLD);
  }

  sprintf(ggstr,"    xxxxxx DONE slave_send_results.\n");
  mylog(mylogf,ggstr);

}
/**************************************-**************************************/
/*                                                                           */
/*                           MM_MASTER_RECEIVE_RESULTS                       */
/*                                                                           */
/*  NOTES                                                                    */
/*  - Do not store results from a suspended (state 2) processor.             */
/*                                                                           */
/*****************************************************************************/
void mm_master_receive_results(r,ti,source)
     struct response_struct *r;  // Response params
     int ti;      // Trial index for results
     int source;  // process id of sender
{
  int i;
  int n,ns,nf,state,tsig;
  float fbuf[MAX_BUF_RESULT];
  MPI_Status status;

  //printf("_____MM_MASTER_RECIEVE_RESULTS ...\n");

  state = mm_get_proc_state(source); // Should be 1 (2=suspended)


  // **** WYETH HERE  ******** Under development Dec 2013 - Feb 2014
  // **** WYETH HERE  ******** Under development Dec 2013 - Feb 2014
  // **** WYETH HERE  ******** Under development Dec 2013 - Feb 2014
  //not = done = here;
  //
  //  Who to send to??
  //  set pnum
  tsig = 1;  // Arbitrary integer, value does not matter?
  MPI_Send(&tsig,1,MPI_INT,source,tag_ok2send,MPI_COMM_WORLD); // OK to send


  //printf("    ___ PID %d   State = %d\n",source,state);

  // Get number of point responses
  //
  //  For lars' 200 x 100 neurons, this ns = 20,000
  //
  MPI_Recv(&ns,1,MPI_INT,source,0,MPI_COMM_WORLD,&status);

  //printf("MASTER: ______***____must receive %d spike trains\n",ns);

  if (ns != r->ns)
    mylog_exit(mylogf,"MM_MASTER_RECEIVE_RESULTS  ns != r->ns\n");
  for(i=0;i<ns;i++){
    MPI_Recv(fbuf,MAX_BUF_RESULT,MPI_FLOAT,source,0,MPI_COMM_WORLD,&status);
    MPI_Get_count(&status,MPI_FLOAT,&n);
    if (state == 1){
      r->cnt[i][ti] = n;
      r->s[i][ti] = copy_farray(fbuf,n);
    }
  }

  //printf("         DONE___ ns = %d\n",ns);

  // Get number of float responses
  MPI_Recv(&nf,1,MPI_INT,source,0,MPI_COMM_WORLD,&status);

  //printf("MASTER: ______***____must receive %d float signals\n",nf);

  if (nf != r->nf)
    mylog_exit(mylogf,"MM_MASTER_RECEIVE_RESULTS  nf != r->nf\n");
  for(i=0;i<nf;i++){
    MPI_Recv(fbuf,MAX_BUF_RESULT,MPI_FLOAT,source,0,MPI_COMM_WORLD,&status);
    MPI_Get_count(&status,MPI_FLOAT,&n);
    if (state == 1){
      r->fcnt[i][ti] = n;
      r->f[i][ti] = copy_farray(fbuf,n);
    }
  }

  if (state == 1)
    r->dflag[ti] = 1; // Results received for this trial

  if (state == 2) // Results received from a suspended proc
    mm_proc_allow_reactivate(source); // will change state to 3

  //printf("_____MM_MASTER_RECIEVE_RESULTS_____done\n");

}
/**************************************-**************************************/
/*                                                                           */
/*                                MM_SLAVE_RUN                               */
/*                                                                           */
/*****************************************************************************/
int mm_slave_run(m,s,r)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response params
{
  int si,tnum,snum,oldsnum,oldmi,done,calgen;
  int xn,yn,tn,nrpt,mi,jcnt,newmoo;
  char *sform,*tname;
  MPI_Status status;

  mylog(mylogf,"  MM_SLAVE_RUN\n");

  jcnt = 0; // Count the number of jobs we do

  wm_get_sform_xn_yn_tn(m,s,&sform,&xn,&yn,&tn);

  nrpt = paramfile_get_int_param_or_exit(s->ppl,"stim_nrpt");

  // Must occur BEFORE 'wm_moo_var_par_set' and AFTER 'wm_response_init'
  // Returned value ignored, because SLAVE does not write response or cal file
  calgen = wmu_calibration_init(m,r,mylogf);  

  oldsnum = -1;   // Should not match new 'snum'
  oldmi   = -1;   // Should not match new 'mi'
  done = 0;
  while (!done){
    mylog(mylogf,"  Waiting to receive a job ...");

    MPI_Recv(&tnum,1,MPI_INT,0,tag_runtrial,MPI_COMM_WORLD,&status);

    if (tnum == -1){
      done = 1;
      mylog(mylogf,"  Instructed to terminate.\n");
    }else{

      //
      //  *** WYETH - SHOULD r->gtsi BE SET TO 'tnum' ???
      //  *** This value is used by all the 'mod_util_resp_store...' routines
      //  It appears to me that each slave makes storage for all jobs, thus
      //  all repeats, all stim, all model configs, but doesn't need this.
      //  Nevertheless, using 'tnum' here means we will be writing sparsely
      //  into our own copy of that global response array.
      //
      r->gtsi = tnum;  // *** WYETH Aug 8th

      // **** OR *******
      // ********* MOVE BELOW AFTER IT IS CALC'd
      // ********* MOVE BELOW AFTER IT IS CALC'd
      // NO - I DON"T THINK SO, r->gtsi = si; 

      // Compute model, stimulus and repeat from trial number
      mi = (int)(tnum / s->ntr);   // Model index
      si = tnum - mi*s->ntr;       // index among Stim x NRpt
      snum = (int)(si / nrpt);
      s->repno = si - snum*nrpt;
      //snum = (int)(tnum / nrpt);
      //s->repno = tnum - snum*nrpt;

      m->m_i = mi;  // Used by models to determine global trial index (gti)

      //sprintf(ggstr,"  Got job:  Stim %d  Rpt %d\n",snum,s->repno);
      sprintf(ggstr,"  Got job:  Mod %d  Stim %d  Rpt %d\n",mi,snum,s->repno);
      mylog(mylogf,ggstr);

      if (mi != oldmi){  // MOOVAR
	// (1) Clean up any current model
	if (oldmi > -1)  // If this is not the first run
	  wm_slave_action(m,s,r,-1,mylogf); // -1 = Done

	// (2) Substitute the model var param values
	// *** WYETH BUG FIX ??? July 17, 2018
	//if (m->nrun > 1)  THIS CREATED TROUBLE if there was only one moo-var
	//  configuration - of course, why would you have just one?  But, it
	//  might be convenient to allow just one, for testing, and so we made
	//  this change in the "if" condition below:
	//
	if (m->vps != NULL) // If MOO params vary
	  wm_moo_var_par_set(mylogf,m,mi); // Update the params of the model

	// Moved up here in case 'prep' needs this information
	r->tri = si;
	r->tsi = si;  // WYETH - Lock these together, sequence means nothing?

	// (3) Call 'prep' for new model
	wm_slave_action(m,s,r,0,mylogf);  //  0 = Prep

	newmoo = 1;
	oldmi = mi;
      }else{
	r->tri = si;  // Be sure to do this
	r->tsi = si;
	newmoo = 0;
      }


      //
      //  Fill s->d... to contain stimulus data, if a new stim is needed
      //
      if ((snum != oldsnum) || (newmoo == 1)){ // This requires new stimulus data  WYETH 2018

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
	    free_3d_iarray(s->dc,xn,yn,tn);  // WYETH - IS THIS RIGHT?
	}else if (strcmp(sform,"3rgb")==0){
	  if (s->drgb_l != NULL)
	    free_4d_farray(s->drgb_l,3,xn,yn,tn);
	}else if (strcmp(sform,"1d")==0){
	  if (s->d1 != NULL)
	    myfree(s->d1);
	}else if (strcmp(sform,"fseq")==0){
	  mylog_exit(mylogf,"MM_SLAVE_RUN  WYETH - must clean up fseq\n");
	}

	get_stimulus_data_tuning_curve(m,s,r,snum,si,mylogf);
	// Stim data now in 's->d', 's->d_r', or 's->d1' for 1d data

	oldsnum = snum;
      }

      //tname = mm_make_var_name(s,tnum);
      tname = mm_make_var_name(s,si);
      sprintf(ggstr,"    Stimulus:  %s\n",tname);
      mylog(mylogf,ggstr);

      // Next call modifies s->d, s->d_r
      wm_slave_action(m,s,r,1,mylogf);  // actionflag 1 means run

      mm_slave_send_results(r,tnum);

      jcnt += 1;  // Count the number of jobs done.

      myfree(tname);
    }
  }
  myfree(sform);

  return jcnt;
}
/**************************************-**************************************/
/*                                                                           */
/*                                MM_MASTER_RUN                              */
/*                                                                           */
/*****************************************************************************/
void mm_master_run(m,s,r,arglist,narg,parflag)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response params
     char **arglist;             // [narg] command line param args
     int narg;                   // number of param args
     int *parflag;               // [narg] flag indicate which in .moo file
{
  int i,k;
  int n,done,tnum,rpti,stimi;
  int nmod,nstim,nrpt,jnum,pnum;
  int flag,aflag;
  int tmin_curr,tmin_last_save;
  MPI_Request *req;
  MPI_Request req1;
  MPI_Status status;

  mylog(mylogf,"  MM_MASTER_RUN\n");

  wm_response_init(m,s,r,mylogf);

  nrpt = paramfile_get_int_param_or_exit(s->ppl,"stim_nrpt");
  nstim = s->ntr / nrpt;

  nmod = m->nrun;

  // Initialize job control, including get names of slaves
  mm_init_job_control(numprocs,nmod,nstim,nrpt,mylogf,processor_name,NULL);

  mm_ja_log_proc_names(tlogf);  // Write to ja file

  if (mm_gui_flag == 1)
    mm_gui_init(r);

  n = s->ntr * nmod;   // Total jobs to complete

  req = NULL; // Avoid -Wall warning

  tmin_curr = tmin_last_save = 0; // Time in minutes, for saving

  // Initial job assignment
  aflag = mm_suggest_next_assignment(&pnum,&jnum,tlogf);
  while(aflag == 1){
    // (job,proc) Master is proc 0
    mm_assign_job(jnum,pnum,mylogf,mm_gui_flag,tlogf);
    aflag = mm_suggest_next_assignment(&pnum,&jnum,tlogf);
  }

  while (mm_get_ndone() < n){  // While not all trials are completed

    //  This is a non-blocking Receive
    MPI_Irecv(&tnum,1,MPI_INT,MPI_ANY_SOURCE,tag_req2send,MPI_COMM_WORLD,
	      &req1);

    //MPI_Irecv(&tnum,1,MPI_INT,MPI_ANY_SOURCE,tag_result,MPI_COMM_WORLD,&req1);
    done = 0;

    while(!done){
      MPI_Test(&req1,&flag,&status);  // Test for incoming results
      if (flag){
	done = 1;
	k = status.MPI_SOURCE; // processor ID that just finished

	/*** WYETH DEBUG ***/
	//mylog(mylogf,"MM_MASTER_RUN  before receive results\n");

	// I believe that the following routine stores data according to
	// the stim index rather than in sequential (received) order.  Thus,
	// 'tribyx' below is set accordingly.
	mm_master_receive_results(r,tnum,k);
	stimi = mm_job_get_stimi(tnum);  // stim index w/i NStim x NRpt
	s->tribyx[tnum] = stimi;

	if (mm_get_proc_state(k) == 1){ // If state was valid (i.e. busy)
	  mm_job_done(r,tnum,k,mm_gui_flag,tlogf); // Log that job completed
	                                     // and send response for monitor
	  if (mm_gui_flag == 0){
	    sprintf(ggstr,"    Job %d done by processor %d\n",tnum,k);
	    mylog(mylogf,ggstr);
	  }
	  
	  // Make new assingment (presumably the loop runs once)
	  aflag = mm_suggest_next_assignment(&pnum,&jnum,tlogf);
	  while(aflag == 1){
	    // (job,proc) Master is proc 0
	    mm_assign_job(jnum,pnum,mylogf,mm_gui_flag,tlogf);
	    aflag = mm_suggest_next_assignment(&pnum,&jnum,tlogf);
	  }
	  
	  // Save partial results, if time interval long enough
	  tmin_curr = mm_time_elapsed_min();
	  if (tmin_curr >= (tmin_last_save + 5)){ // 5 minutes
	    wm_write_response(m,s,r,mylogf,myid);
	    tmin_last_save = tmin_curr;
	  }
	}
      }else{

	// Make new assingment (presumably the loop runs once)
	aflag = mm_suggest_next_assignment(&pnum,&jnum,tlogf);
	while(aflag == 1){
	  // (job,proc) Master is proc 0
	  mm_assign_job(jnum,pnum,mylogf,mm_gui_flag,tlogf);
	  aflag = mm_suggest_next_assignment(&pnum,&jnum,tlogf);
	}

	if (mm_gui_flag == 1){
	  mm_gui_update();
	  mm_gui_event();
	  mm_mon_event(r); // Check for events in response monitor window
	}
	usleep(10000); // So we don't run the processor all the time
	//sleep(1);
      }
    }
  }
  if (mm_gui_flag == 1)
    mm_gui_hold(r,-1);

  mm_print_stat();
}
/**************************************-**************************************/
/*                                                                           */
/*                                   MM_FIT                                  */
/*                                                                           */
/*****************************************************************************/
void mm_fit(m,s,r,arglist,narg,parflag)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response params
     char **arglist;             // [narg] command line param args
     int narg;                   // number of param args
     int *parflag;               // [narg] flag indicate which in .moo file
{
  int i,k;
  int n,done,tnum,rpti,stimi;
  int nmod,nstim,nrpt,jnum,pnum;
  int flag,aflag;
  int tmin_curr,tmin_last_save;
  MPI_Request *req;
  MPI_Request req1;
  MPI_Status status;

  mylog(mylogf,"  MM_FIT\n");

  if (myid == 0){

    wm_slave_action(m,s,r,2,mylogf);  //  2 = Fit

    mylog_exit(mylogf," *** MM_FIT --- Master - Under development\n");
  }else{

    wm_slave_action(m,s,r,2,mylogf);  //  2 = Fit

    mylog_exit(mylogf," *** MM_FIT --- Slave - Under development\n");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                 MM_MAKE_MSR                               */
/*                                                                           */
/*****************************************************************************/
void mm_make_msr(paramfile,stimfile,respfile,mo,sppl,rppl,rm,rs,rr)
     char paramfile[],stimfile[],respfile[];
     struct onode *mo;
     struct param_pair_list *sppl;
     struct param_pair_list *rppl;
     struct model_struct **rm;    // Model params
     struct stim_struct **rs;     // Stimulus params
     struct response_struct **rr; // Response params
{
  struct model_struct *m;    // Model params
  struct stim_struct *s;     // Stimulus params
  struct response_struct *r; // Response params

  // Initialize stimulus struct
  s = (struct stim_struct *)myalloc(sizeof(struct stim_struct));
  s->ppl = sppl;
  strcpy(s->paramfile,stimfile);
  s->name = (char *)NULL;
  s->d = (float ***)NULL;
  s->d_r = (float ***)NULL;
  s->d1 = (float *)NULL;
  s->dc = (int ***)NULL;
  s->dc_r = (int ***)NULL;

  // Initialize model struct
  m = (struct model_struct *)myalloc(sizeof(struct model_struct));
  m->ppl = NULL;
  m->o = mo;
  strcpy(m->paramfile,paramfile);
  m->process_id = myid;
  m->process_n = numprocs;
  m->run_mode = 0;    // Normal run mode.  No interactive mode for mm
  m->marfile = NULL;
  m->action = -999;   // Arbitrary value
  m->calib_pval = NULL;
  m->mds = NULL;
  m->mds_n = 0;

  m->nrun = 1;        // Assume one model, change below if moo var pars
  m->vps = NULL;


  // Initialize response struct
  r = (struct response_struct *)myalloc(sizeof(struct response_struct));
  r->ppl = rppl;
  strcpy(r->paramfile,respfile);
  r->ns = -1; // Number of spike responses per trial
  r->nf = -1; // Number of continuous responses per trial

  *rm = m; *rs = s; *rr = r;
}
/**************************************-**************************************/
/*                                                                           */
/*                                  MM_SLAVE                                 */
/*                                                                           */
/*****************************************************************************/
void mm_slave(mmode,paramfile,stimfile,respfile)
     char mmode[];
     char paramfile[],stimfile[],respfile[];
{
  int cn,jcnt;
  char *c,*ct;
  struct onode *mo;
  struct param_pair_list *sppl,*rppl;
  struct model_struct *m;    // Model params
  struct stim_struct *s;     // Stimulus params
  struct response_struct *r; // Response params
  MPI_Status status;

  // Initialize model struct
  mm_mpi_receive_carray(0,tag_mppl,100000,mylogf,&c,&cn);
  mylog(mylogf,"  MM_SLAVE received carray for mppl\n");
  ct = c;  // The 'ct' param will be changed by the next call
  mo = paramfile_get_onode_from_carray(&ct);
  myfree(c);

  // Initialize stim struct
  mm_mpi_receive_carray(0,tag_sppl,100000,mylogf,&c,&cn);
  mylog(mylogf,"  MM_SLAVE received carray for sppl\n");
  sppl = paramfile_carray_to_ppl(c,cn);
  myfree(c);

  // Initialize response struct
  if (strcmp(mmode,"fit")==0){
    rppl = NULL;
  }else{
    mm_mpi_receive_carray(0,tag_rppl,100000,mylogf,&c,&cn);
    mylog(mylogf,"  MM_SLAVE received carray for rppl\n");
    rppl = paramfile_carray_to_ppl(c,cn);
    myfree(c);
  }

  // Send the name of the slave processor (+1 to send '\0')
  MPI_Send(processor_name,strlen(processor_name)+1,MPI_CHAR,0,tag_procname,
	   MPI_COMM_WORLD);

  mm_make_msr(paramfile,stimfile,respfile,mo,sppl,rppl,&m,&s,&r);

  wm_prep_var_const_gen(m,s,mylogf,NULL,0,NULL);



  if (strcmp(mmode,"fit")==0){
    mm_fit(m,s,r,NULL,0,NULL);  // For fitting
    jcnt = 1; // WYETH - Assuming this is OK
  }else{
    wm_response_init(m,s,r,mylogf);

    jcnt = mm_slave_run(m,s,r);  // Returns the number of jobs done
  }

  //
  // Send 'done' signal to slave
  //
  //  Wyeth - currently it is possible that this gets called even when
  //  no jobs were run by this slave, and thus no "PREP" was ever done.
  //  In this case, we don't want to carry out any clean-up.
  //
  if (jcnt > 0){
    if (strcmp(mmode,"fit")!=0){
      sprintf(ggstr,"    Did %d jobs\n",jcnt);
      mylog(mylogf,ggstr);
    }
    wm_slave_action(m,s,r,-1,mylogf); // -1 = Done
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                 MM_MASTER                                 */
/*                                                                           */
/*****************************************************************************/
void mm_master(mmode,paramfile,stimfile,respfile,arglist,narg)
     char mmode[];
     char paramfile[],stimfile[],respfile[];
     char **arglist;
     int narg;
{
  int i;
  int cn;
  int *parflag;              // Flags for cmd line args used in moo file
  char *c,*calmode;
  struct param_pair_list *sppl,*rppl;
  struct model_struct *m;    // Model params
  struct stim_struct *s;     // Stimulus params
  struct response_struct *r; // Response params
  struct onode *mo,*co;

  if (is_suffix_string(paramfile,".mod"))
    mylog_exit(mylogf,"MM_MASTER  'mm' no longer accepts .mod files.\n");

  mo = paramfile_onode_file_read(paramfile,1);
  paramfile_onode_resolve(mo,mylogf); // Resolve 'COPY' assignments
  parflag = onode_update_from_slist(mo,mylogf,arglist,narg);

  sppl = paramfile_get_initial_params_update_only(stimfile,arglist,narg,1);

  if (strcmp(mmode,"fit")==0){
    // THERE IS NO .rsp FILE FOR "fit"
    rppl = NULL;
    // Check if all command line params. occur in at least one param file
    for(i=0;i<narg;i+=2){
      if (!(onode_test_ostr(mo,arglist[i])||
	    paramfile_test_param(sppl,arglist[i]))){
	sprintf(ggstr,"  *** Parameter name:  %s\n",arglist[i]);
	mylog(mylogf,ggstr);
	mylog_exit(mylogf,"MM_MASTER  Command line param not found\n");
      }
    }
  }else{
    rppl = paramfile_get_initial_params_update_only(respfile,arglist,narg,1);

    // Check if all command line params. occur in at least one param file
    for(i=0;i<narg;i+=2){
      if (!(onode_test_ostr(mo,arglist[i])||
	    paramfile_test_param(sppl,arglist[i])||
	    paramfile_test_param(rppl,arglist[i]))){
	sprintf(ggstr,"  *** Parameter name:  %s\n",arglist[i]);
	mylog(mylogf,ggstr);
	mylog_exit(mylogf,"MM_MASTER  Command line param not found\n");
      }
    }
  }

  mm_make_msr(paramfile,stimfile,respfile,mo,sppl,rppl,&m,&s,&r);

  //  This may add a child to m->o for <var_param>, so do before sending 'm'
  wm_prep_var_const_gen(m,s,mylogf,arglist,narg,parflag);

  //
  //  Send the param lists for Model, Stim, and Resp to all slaves
  //
  c = paramfile_onode_carray(mo,&cn);
  for(i=1;i<numprocs;i++)
    MPI_Send(c,cn,MPI_CHAR,i,tag_mppl,MPI_COMM_WORLD);
  myfree(c);

  paramfile_ppl_to_carray(sppl,&c,&cn);
  for(i=1;i<numprocs;i++)
    MPI_Send(c,cn,MPI_CHAR,i,tag_sppl,MPI_COMM_WORLD);
  myfree(c);

  if (strcmp(mmode,"fit")!=0){
    paramfile_ppl_to_carray(rppl,&c,&cn);
    for(i=1;i<numprocs;i++)
      MPI_Send(c,cn,MPI_CHAR,i,tag_rppl,MPI_COMM_WORLD);
    myfree(c);
  }


  //
  //  Run all trials, accumulate results in 'r'
  //
  if (strcmp(mmode,"fit")==0){
    mm_fit(m,s,r,arglist,narg,parflag);  // For fitting
    exit(0);
  }else
    mm_master_run(m,s,r,arglist,narg,parflag);

  for(i=1;i<numprocs;i++)
    mm_assign_job(-1,i,mylogf,mm_gui_flag,tlogf); // Terminate slaves

  wm_write_response(m,s,r,mylogf,myid);

  //
  //  This exits because child node not found.
  //

  co = onode_child_get_unique(m->o,"calibrate");
  if (co != NULL){
    calmode = onode_getpar_child_chr_ptr_exit(m->o,"calibrate","mode");
    if (strcmp(calmode,"generate")==0) //  Write calibration file WYETH 2018 Dec
      wmu_calibration_write(m,s,r,mylogf);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                 MM_REPLAY                                 */
/*                                                                           */
/*****************************************************************************/
void mm_replay(jafile,ndfile)
     char jafile[];  // Job assignment file
     char ndfile[];  // ndata file
{
  int i;
  FILE *fopen(),*fin;
  int nproc,nstim,nrpt,pnum,jnum,tmsec,done,ntr,dt,pflag,delay_10th_sec,ns;
  int nmod;
  unsigned long tnow,t0ms;
  char **pname,tstr[LONG_SLEN];
  struct response_struct *r;
  struct ndata_struct *nd;
  struct ndtrial_struct *t0;

  printf("  MM_REPLAY\n");

  pflag = 0; // Print flag, for debugging

  //  Read nd file
  read_ndata(ndfile,&nd);


  //  Open the Job Assignment file
  if ((fin = fopen(jafile,"r")) == NULL){
    printf("  *** Cannot open file %s.\n",jafile);
    exit(0);
  }else
    printf("    Reading job log file:  %s\n",jafile);

  // Read list of processor names
  ns = fscanf(fin,"%*s %d",&nproc);
  printf("    NProc %d\n",nproc);
  pname = (char **)myalloc(nproc*sizeof(char *));
  for(i=0;i<nproc;i++){
    ns = fscanf(fin,"%*s %d %s",&pnum,tstr);
    pname[i] = strdup(tstr);
  }
  for(i=0;i<nproc;i++){
    printf("   %4d  %s\n",i,pname[i]);
  }

  // Read NSTIM NRPT
  nmod = 1; // *** WYETH FIX for moovar ??? ****
  ns = fscanf(fin,"%*s %d",&nstim);
  ns = fscanf(fin,"%*s %d",&nrpt);
  printf("    NStim %d\n",nstim);
  printf("    NRpt %d\n",nrpt);
  printf("    *** Assuming:  NMod %d\n",nmod);

  ntr = nmod*nstim*nrpt;  // Total number of trials

  mm_init_job_control(nproc,nmod,nstim,nrpt,NULL,pname[0],pname);


  //
  //  Create Response Structure 'r'
  //
  r = (struct response_struct *)myalloc(sizeof(struct response_struct));
  r->paramfile[0] = '\0';
  r->ppl = NULL;

  //
  //  WYETH - COULD GET FROM NDATA FILE
  //
  r->n = 1;  // Number of responses per trial
  r->rtype = (char **)myalloc(r->n * sizeof(char *));
  r->rtype[0] = strdup("s");
  r->ns = 1;  // Number of "s" responses

  r->nd_name = (char **)myalloc(r->n * sizeof(char *));
  r->nd_name[0] = strdup("spikes");

  t0 = &(nd->t[0]);
  r->tsn = (float)t0->r[0].tn / t0->r[0].sampling;

  // float ***s;                  // Spike train responses [ns][ti][cnt]
  //   int **cnt;                   // Spike count [ns][ti]  -1 MEANS UNFILLED

  r->s   = get_3d_farray(r->ns,ntr,0);
  r->cnt = get_2d_iarray(r->ns,ntr);
  for(i=0;i<ntr;i++){
    t0 = &(nd->t[i]);
    r->s[0][i] = i2farray(t0->r[0].p,t0->r[0].n);
    r->cnt[0][i] = t0->r[0].n;
  }
  if (ntr != nd->ntrial){
    exit_error("MM_REPLAY","ndata trials does not match nstim x nrpt.");
  }


  printf("    Trial duration %f s\n",r->tsn);

  mm_gui_init(r);

  printf("\n    *** Press 'q' to begin ***\n\n");    
  mm_gui_hold(r,1);  // Wait for 'q' to be pressed
  printf("    *** Press 'q' to exit ***\n\n");    

  done = 0;
  while(!done){

    //fscanf(fin,"%s",&tstr);
    ns = fscanf(fin,"%s",tstr);
    if (strcmp(tstr,"Assign")==0){

      ns = fscanf(fin,"%*s %d %*s %d",&jnum,&pnum);
      if (pflag) printf("  Assign job %d proc %d\n",jnum,pnum);

      mm_assign_job(jnum,pnum,NULL,mm_gui_flag,NULL);


      if (jnum == -1)
	done = 1;

    }else if (strcmp(tstr,"jnum")==0){

      ns = fscanf(fin,"%d %*s %d %*s %d",&jnum,&pnum,&tmsec);
      if (pflag) printf("  Done job %d proc %d in %d msec\n",jnum,pnum,tmsec);

      t0ms = mm_get_start_time_for_job(jnum);  // Start time (msec)
      tnow = mm_current_time_ms();
      dt = (int)(tnow - t0ms);

      // Wait if necessary for the real time the job would have taken
      if (tmsec > dt){
	
	if (pflag) printf("   Only %d of %d ms elapsed, Sleep it off\n",
			  dt,tmsec);

	// Every second, update the gui so the time changes in the lower right
	delay_10th_sec = (tmsec-dt)/100; // Number of 1/10 second delays
	for(i=0;i<delay_10th_sec;i++){
	  usleep(100000); // Sleep for 0.1 sec
	  mm_gui_update();
	  mm_mon_event(r);
	}
      }

      mm_job_done(r,jnum,pnum,mm_gui_flag,NULL); // Log that job completed


    }else{
      printf("  *** Found  %s\n",tstr);
      exit_error("MM_REPLAY","Unexpected string.");
    }

  }
  fclose(fin);

  mm_gui_hold(r,1);  // Wait for 'q' to be pressed

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
  int namelen,sysres;
  char option[SLEN],paramfile[SLEN],stimfile[SLEN],respfile[SLEN],tstr[SLEN];
  char jafile[SLEN],ndfile[SLEN];

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid); // Set global value
  MPI_Get_processor_name(processor_name,&namelen);
  MPI_Buffer_attach(malloc(BUFFSIZE),BUFFSIZE); // WYETH - necessary?
  
  if (myid == 0){
    mylogf = NULL;

    // Create a 'ja' (job assignment) file to log the job control
    tlogf = mylog_get_global_ja_fname();
    remove_file(tlogf);

    sprintf(tstr,"NPROC %d\n",numprocs);  // Append number of processors
    append_string_to_file(tlogf,tstr);

  }else{
    mylogf = mylog_get_global_fname(myid);
    remove_file(mylogf);
    mylog(mylogf,"Begin\n");
    sprintf(tstr,"chmod ugo+rwx %s",mylogf);
    sysres = system(tstr);
  }

  if (argc == 1){
    if (myid==0){
      mylog(mylogf,"\n  mm <option>\n");
      mylog(mylogf,"    where <option> is\n");
      mylog(mylogf,"      mod <paramfile> <stimfile> <respfile>");
      mylog(mylogf," [<pname> <pval> ...]\n");
      mylog(mylogf,"      modm <paramfile> <stimfile> <respfile>");
      mylog(mylogf," [<pname> <pval> ...]\n");
      mylog(mylogf,"      modq <paramfile> <stimfile> <respfile>");
      mylog(mylogf," [<pname> <pval> ...]\n");
      mylog(mylogf,"        The quiet option - no graphics.\n");
      mylog(mylogf,"      replay <ja_file> <ndata_file>\n");
      mylog(mylogf,"      fit <paramfile> <stimfile>\n");
      mylog(mylogf,"\n");
    }
    fflush(stdout);
    MPI_Finalize();
    mylog(mylogf,"  Version 3.0\n\n");
    return 0;
  }else{
    strcpy(option,argv[1]);
  }

  if ((strcmp(option,"mod")==0)||(strcmp(option,"modm")==0)||
      (strcmp(option,"modq")==0)||(strcmp(option,"fit")==0)){

    if (strcmp(option,"modm")==0){ 
      mctrl_util_determine_dir_structure();
      if (myid == 0){
	mylog(mylogf,"<br>\n<pre>\n");
	mctrl_util_report_pid_run();
      }
      mm_gui_flag = 0;
    }

    if (strcmp(option,"modq")==0){ 
      mm_gui_flag = 0;
    }

    strcpy(paramfile,argv[2]);
    strcpy(stimfile, argv[3]);
    if (strcmp(option,"fit")==0)
      respfile[0] = '\0';   // No respfile for fitting
    else
      strcpy(respfile, argv[4]);

    if (myid == 0){
      mm_master(option,paramfile,stimfile,respfile,&(argv[5]),argc-5);

      // *** WYETH ADDED FOR TEMP CHANGE TO KAMBLIPOOCHI
      append_string_to_file("JOB_DONE","Done\n");

    }else{
      mm_slave(option,paramfile,stimfile,respfile);
    }

  }else if (strcmp(option,"replay")==0){
    strcpy(jafile,argv[2]);
    strcpy(ndfile,argv[3]);
    mm_replay(jafile,ndfile);
  }else
    mylog(mylogf,"  *** Unknown option.\n");

  mylog(mylogf,"\n");
  
  fflush(stdout);
  MPI_Finalize();
  return 0;
}
