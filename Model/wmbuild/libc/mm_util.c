/*****************************************************************************/
/*                                                                           */
/*  mm_util.c                                                                */
/*  wyeth bair                                                               */
/*                                                                           */
/*  ISSUES for moo var                                                       */
/*   - Response monitor, needs to use 'mm_mod_mi'                            */
/*                                                                           */
/*****************************************************************************

  rm ../lib/libmm_util.a
  mpicc -c mm_util.c
  ar ruv ../lib/libmm_util.a mm_util.o  ;  ranlib ../lib/libmm_util.a

******************************************************************************/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>  // For usleep()
#include <string.h>
//#include <X11/Xlib.h>
//#include <GL/glx.h>
//#include <GL/gl.h>
#include <math.h>
#include <sys/types.h> // For times()
#include <sys/times.h> // For times()
#include <time.h>      // For time(), localtime(), time_t, tm

#include "my_util.h"
#include "carray_util.h"
#include "farray_util.h"
//#include "glplot_util.h"
#include "mod_util.h"
#include "mod.h"
#include "mm.h"

#define MAX_DATA_N_INT   1000000  // Arbitrary max size of data array
#define MAX_DATA_N_FLOAT 1000000  // Arbitrary max size of data array

struct mm_job_struct ***mm_jtable = NULL; // [nmod][nstim][nrpt]
int mm_nmod;                              // Number of models
int mm_nstim;                             // Number of stimuli
int mm_nrpt;                              // Number of repeats
int mm_ndone;                             // Number of jobs done
int mm_ndoa;                              // Number of jobs done or assigned
unsigned long mm_t0;                      // Start time (ms)
unsigned long mm_t1;                      // Time last job finished (ms)

struct mm_proc_struct *mm_ptable = NULL;
int mm_nproc;

char ggstr[LONG_SLEN];

//  GUI
int mm_winw = 500;
int mm_winh = 500;
float mm_zoom = 1.0;
//Display *mm_dpy;    *** WYETH 2019 comment out GL/X dependency
//Window mm_win;
//GLXContext mm_cx;
int mm_dbflag;
float mm_pixel_proc_top;     // Pixel location of top processor in list
float mm_pixel_proc_height;  // Pixel height of a processor name
int mm_proc_current;         // Index of current processor, for highlight

int mm_gui_hold_flag = 0;    // 1-wait at end before closing windows

//  MM Response Monitor
int mm_monflag = 0;
int mm_moni;            // Index of response in 'r' to monitor
int mm_mon_mi = 0;      // Model index of response in 'r' to monitor
int mm_sort = 0;        // How to sort trials:  0-stimulus, 1-processor
int mm_winw2 = 1024;
int mm_winh2 = 600;
//Display *mm_dpy2;   *** WYETH 2019 comment out GL/X dependency
//Window mm_win2;
//GLXContext mm_cx2;
int mm_dbflag2;

//  Constants for Drawing
int   mm_nlines;      // Number of lines
float mm_line_h;      // line height
float mm_spike_h;     // spike height


/*****************************************************************************/
/*
Bool mm_true_test()
{
  return 1;
}
*/
/**************************************-**************************************/
/*                                                                           */
/*                             MM_DATA_SEND_IARRAY                           */
/*                                                                           */
/*  Send one or more ints.                                                   */
/*                                                                           */
/*****************************************************************************/
void mm_data_send_iarray(p0,pn,d,n,mylogf)
     int p0;        // process number (or starting number) to send to
     int pn;        // Number of processes to send to
     int *d;        // [n] Data to send
     int n;         // Length of data array
     char *mylogf;  // log file
{
  int i;

  if (p0<0)
    mylogx(mylogf,"MM_DATA_SEND_IARRAY","Processor ID too low");

  if (pn<1)
    mylogx(mylogf,"MM_DATA_SEND_IARRAY","Processor count too low");

  if (n<0)
    mylogx(mylogf,"MM_DATA_SEND_IARRAY","Invalid length");

  if (n > MAX_DATA_N_INT)
    mylogx(mylogf,"MM_DATA_SEND_IARRAY","Length exceeds arbitrary maximum");

  for(i=p0;i<(p0+pn);i++)
    MPI_Send(d,n,MPI_INT,i,tag_data_int,MPI_COMM_WORLD);
}
/**************************************-**************************************/
/*                                                                           */
/*                             MM_DATA_RECV_IARRAY                           */
/*                                                                           */
/*  Receive one or more ints.                                                */
/*                                                                           */
/*****************************************************************************/
void mm_data_recv_iarray(source,mylogf,d,dmax,rn)
     int source;    // process number to receive from, -1 for MPI_ANY_SOURCE
     char *mylogf;  // log file
     int *d;        // [dmax] array to be filled
     int dmax;      // Maximum capacity of 'd'
     int *rn;       // Number of ints returned in 'd'
{
  MPI_Status stat;

  if (source == -1)
    source = MPI_ANY_SOURCE;
  else if (source < -1)
    mylogx(mylogf,"MM_DATA_RECV_IARRAY","Source ID too low");

  if (dmax<0)
    mylogx(mylogf,"MM_DATA_RECV_IARRAY","Invalid length");

  MPI_Recv(d,MAX_DATA_N_INT,MPI_INT,source,tag_data_int,MPI_COMM_WORLD,&stat);

  MPI_Get_count(&stat,MPI_INT,rn);

  if (*rn > dmax)
    mylogx(mylogf,"MM_DATA_RECV_IARRAY","Too many ints");
}
/**************************************-**************************************/
/*                                                                           */
/*                             MM_DATA_SEND_FARRAY                           */
/*                                                                           */
/*  Send one or more floats.                                                 */
/*                                                                           */
/*****************************************************************************/
void mm_data_send_farray(p0,pn,d,n,mylogf)
     int p0;        // process number (or starting number) to send to
     int pn;        // Number of processes to send to
     float *d;      // Pointer to value(s) to send
     int n;         // Number of values to send
     char *mylogf;  // log file
{
  int i;

  if (p0<0)
    mylogx(mylogf,"MM_DATA_SEND_FLOAT","Processor ID too low");

  if (pn<1)
    mylogx(mylogf,"MM_DATA_SEND_FLOAT","Processor count too low");

  if (n > MAX_DATA_N_FLOAT)
    mylogx(mylogf,"MM_DATA_SEND_FLOAT","Length exceeds arbitrary maximum");

  for(i=p0;i<(p0+pn);i++)
    MPI_Send(d,n,MPI_FLOAT,i,tag_data_flt,MPI_COMM_WORLD);
}
/**************************************-**************************************/
/*                                                                           */
/*                             MM_DATA_RECV_FARRAY                           */
/*                                                                           */
/*  Receive one or more floats.                                              */
/*                                                                           */
/*****************************************************************************/
void mm_data_recv_farray(source,mylogf,d,dmax,rn,rsrc)
     int source;    // process number to receive from, -1 for MPI_ANY_SOURCE
     char *mylogf;  // log file
     float *d;      // [dmax] array to be filled
     int dmax;      // Maximum capacity of 'd'
     int *rn;       // Number of float returned in 'd'
     int *rsrc;     // Return the source, if 'source' was -1
{
  MPI_Status stat;

  if (source == -1)
    source = MPI_ANY_SOURCE;
  else if (source < -1)
    mylogx(mylogf,"MM_DATA_RECV_FARRAY","Source ID too low");

  if (dmax < 0)
    mylogx(mylogf,"MM_DATA_RECV_FARRAY","Invalid length");

  MPI_Recv(d,MAX_DATA_N_FLOAT,MPI_FLOAT,source,tag_data_flt,MPI_COMM_WORLD,
	   &stat);

  if (source == MPI_ANY_SOURCE)
    *rsrc = stat.MPI_SOURCE; // processor ID that just finished

  MPI_Get_count(&stat,MPI_FLOAT,rn);

  if (*rn > dmax)
    mylogx(mylogf,"MM_DATA_RECV_FARRAY","Too many floats");
}
/**************************************-**************************************/
/*                                                                           */
/*                           MM_MPI_RECEIVE_CARRAY                           */
/*                                                                           */
/*  The character data received may have '\0' characters, thus is not        */
/*  limited to null terminated strings.                                      */
/*                                                                           */
/*****************************************************************************/
void mm_mpi_receive_carray(source,tag,nmax,mylogf,rc,rn)
     int source;    // process to receive from
     int tag;       // mpi tag
     int nmax;      // maximum number of characters
     char mylogf[]; // log file
     char **rc;     // return the string
     int *rn;       // number of elements in carray
{
  int i;
  int n;
  char *pstr,*c;
  MPI_Status status;

  //sprintf(ggstr,"(mm_util.c) mm_mpi_receive_carray  tag=%d nmax=%d\n",
  //tag,nmax);
  //mylog(mylogf,ggstr);

  pstr = (char *)myalloc(nmax*sizeof(char));
  MPI_Recv(pstr,nmax,MPI_CHAR,source,tag,MPI_COMM_WORLD,&status);
  MPI_Get_count(&status,MPI_CHAR,&n);
  /*sprintf(ggstr,"    %d chars recieved\n",n);
    mylog(mylogf,ggstr);*/

  c = (char *)myalloc(n*sizeof(char));
  for(i=0;i<n;i++)
    c[i] = pstr[i];

  myfree(pstr);

  *rc = c; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                                MM_CMD_SEND                                */
/*                                                                           */
/*  Send a string command.  Commands should be no longer than 256 chars,     */
/*  including a terminating '\0'                                             */
/*                                                                           */
/*****************************************************************************/
void mm_cmd_send(p0,pn,str,mylogf)
     int p0;        // process number (or starting number) to send to
     int pn;        // Number of processes to send to
     char *str;     // string to send
     char *mylogf;  // log file
{
  int i;

  if (p0<0)
    mylogx(mylogf,"MM_CMD_SEND","Processor ID too low");

  if (pn<1)
    mylogx(mylogf,"MM_CMD_SEND","Processor count too low");

  if (strlen(str) >= SLEN-1)
    mylogx(mylogf,"MM_CMD_SEND","Command string too long");

  for(i=p0;i<(p0+pn);i++){
    //printf("MASTER SEND:  i=%d  strlen = %d\n",i,(int)strlen(str));
    MPI_Send(str,strlen(str)+1,MPI_CHAR,i,tag_cmd_str,MPI_COMM_WORLD);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                MM_CMD_RECV                                */
/*                                                                           */
/*  Receive a string command.  The result should be \0 terminated.           */
/*                                                                           */
/*****************************************************************************/
char *mm_cmd_recv(p0,mylogf)
     int p0;        // process number
     char *mylogf;  // log file
{
  int n;
  char *tstr;

  mm_mpi_receive_carray(p0,tag_cmd_str,SLEN,mylogf,&tstr,&n);

  return tstr;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MM_JOB_GET_STIMI                             */
/*                                                                           */
/*****************************************************************************/
int mm_job_get_stimi(jnum)
     int jnum;  // Job number
{
  int k;
  int mi;

  mi = (int)(jnum/(mm_nstim * mm_nrpt));  // Model number
  k = jnum - mi*(mm_nstim * mm_nrpt);

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MM_JOB_GET_INDICES                            */
/*                                                                           */
/*****************************************************************************/
void mm_job_get_indices(jnum,rmi,rsi,rri)
     int jnum;  // Job number
     int *rmi;  // return model index
     int *rsi;  // return stim index
     int *rri;  // return rpt index
{
  int k;
  int mi,si,ri;

  mi = (int)(jnum/(mm_nstim * mm_nrpt));  // Model number
  k = jnum - mi*(mm_nstim * mm_nrpt);
  si = (int)(k/mm_nrpt);  // Stimulus number
  ri = k - si*mm_nrpt;    // Repeate number

  *rmi = mi;
  *rsi = si;
  *rri = ri;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MM_JOB_GET_PTR                              */
/*                                                                           */
/*****************************************************************************/
struct mm_job_struct *mm_job_get_ptr(jnum)
     int jnum;  // Job number
{
  int mi,si,ri;
  struct mm_job_struct *jp;

  mm_job_get_indices(jnum,&mi,&si,&ri);
  jp = &(mm_jtable[mi][si][ri]);

  return jp;
}
/**************************************-**************************************/
/*                                                                           */
/*                                MM_MON_HELP                                */
/*                                                                           */
/*****************************************************************************/
void mm_mon_help()
{
  printf("\n");
  printf("  MM_MON_HELP\n");
  printf("    h - help (print this list)\n");
  printf("    p - sort trials by processor number\n");
  printf("    s - sort trials by stimulus number\n");
  printf("    > - increase model index for monitor display\n");
  printf("    < - decrease model index for monitor display\n");
  printf("    ? - help (print this list)\n");
  printf("\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                           MM_MON_SET_LINES_HEIGHT                         */
/*                                                                           */
/*****************************************************************************/
void mm_mon_set_lines_height()
{
  if (mm_sort == 1){ // by processor
    mm_nlines = mm_nstim * mm_nrpt + mm_nproc - 1;
    mm_line_h = 1.8 / (float)mm_nlines;  // line height, allowing 10% extra
  }else{ // by stimulus
    mm_nlines = mm_nstim * (mm_nrpt + 1) - 1;
    mm_line_h = 1.8 / (float)mm_nlines;  // line height, allowing 10% extra
  }

  mm_spike_h = 0.9 * mm_line_h;  // Spike height

}
/**************************************-**************************************/
/*                                                                           */
/*                                 MM_MON_INIT                               */
/*                                                                           */
/*  Before calling this, it makes sense to check for 'config' events.        */
/*                                                                           */
/*****************************************************************************/
void mm_mon_init(r)
     struct response_struct *r;  // Response params, can be NULL
{
  int i;
  int n,li,dur_msec;
  float x,y,w,h,yoff,xtext,colr,colg,colb,tr,tg,tb;
  char tstr[SLEN3],chan_name[SLEN];
  struct mm_proc_struct *p;

  exit_error("MM_MON_INIT","Disabled");

  /*

  // Extract values from 'r'
  strcpy(chan_name,r->nd_name[mm_moni]);
  dur_msec = my_rint(r->tsn*1000.0);

  x = -0.90;
  xtext = -0.95;
  w = 1.8;

  colr = 0.1;
  colg = 0.2;
  colb = 0.2;

  glplot_clear_double(mm_dpy2,mm_win2);

  //
  //  Response channel name, top center
  //
  if (mm_nmod > 1)
    sprintf(tstr,"%s   (model %d)",chan_name,mm_mon_mi);
  else
    sprintf(tstr,"%s",chan_name);
  glplot_draw_text(mm_dpy2,mm_win2,mm_cx2,tstr,0.0,0.95,0.0,1.0,0.0);

  //
  //  Draw colored boxes for response grouping
  //
  if (mm_sort == 1){ // By processor
    tr = 1.0;
    tg = 1.0;
    tb = 0.0;
    y = 0.90;
    for(i=0;i<mm_nproc;i++){
      p = &(mm_ptable[i]);
      if (p->ndone > 0){
	h = (float)p->ndone * mm_line_h;  // height of box relates to ndone
	y -= h;                           // y-coord for this box
	yoff = (h - mm_line_h) / 2.0;
	glplot_draw_box(colr,colg,colb,x,y,w,h);
	sprintf(tstr,"%3d",i);
	glplot_draw_text(mm_dpy2,mm_win2,mm_cx2,tstr,xtext,y+yoff,tr,tg,tb);
	y -= mm_line_h;            // leave a blank space for one box
      }
    }
  }else{  // By stimulus
    tr = 0.4;
    tg = 0.8;
    tb = 0.8;

    n = mm_nstim;
    h = mm_nrpt * mm_line_h;
    yoff = (h - mm_line_h) / 2.0;
    for(i=0;i<n;i++){
      y = -0.90 + (float)i * (h + mm_line_h); // y-coord for this line
      glplot_draw_box(colr,colg,colb,x,y,w,h);
      sprintf(tstr,"%3d",n-1-i);
      glplot_draw_text(mm_dpy2,mm_win2,mm_cx2,tstr,xtext,y+yoff,tr,tg,tb);
    }
  }

  //
  //  Time label for x-axis, bottom center
  //
  sprintf(tstr,"%d ms",dur_msec);
  glplot_draw_text(mm_dpy2,mm_win2,mm_cx2,tstr,0.0,-0.97,tr,tg,tb);

  glClear(GL_DEPTH_BUFFER_BIT);    // Later draws will be on top, Dec 2009

  */
}
/**************************************-**************************************/
/*                                                                           */
/*                            MM_MON_UPDATE_DRAW_ONE                         */
/*                                                                           */
/*****************************************************************************/
void mm_mon_update_draw_one(r,pnum,ti,si,ri)
     struct response_struct *r;  // Response params
     int pnum;                   // Processor number
     int ti;                     // Trial number
     int si;                     // Stimulus number
     int ri;                     // Repeat number
{
  int i;
  int n,li;
  float *s,cfx,x,y,x0;
  char tstr[SLEN];

  exit_error("MM_MON_UPDATE_DRAW_ONE","Disabled");

  /*

  // Compute conversion factor from trial duration
  cfx = 1.8 / (r->tsn * 1000.0);  // GL_Coord / msec
  x0 = -0.90;

  // Determine y-value for current line to draw
  li = ri + si * (mm_nrpt + 1);   // Line index
  li = mm_nlines - li;            // Write from top down (not bottom up)
  y = -0.90 + (float)(li-1) * mm_line_h; // y-coord for this line

  // Get spikes on 0th channel
  s = r->s[0][ti];
  n = r->cnt[0][ti];

  // Draw spikes
  for(i=0;i<n;i++){
    x = x0 + cfx * s[i];
    glBegin(GL_LINE_STRIP);
    glColor3f(1.0,1.0,1.0);   // Spikes are white
    glVertex2f(x,y);
    glVertex2f(x,y+mm_spike_h);
    glEnd();
  }

  // Draw stimulus or processor number
  if (mm_sort == 1)
    sprintf(tstr,"%d",si);
  else
    sprintf(tstr,"%d",pnum);
  glplot_draw_text(mm_dpy2,mm_win2,mm_cx2,tstr,0.93,y,1.0,1.0,0.0);
  */
}
/**************************************-**************************************/
/*                                                                           */
/*                            MM_MON_UPDATE_BY_PROC                          */
/*                                                                           */
/*****************************************************************************/
void mm_mon_update_by_proc(r)
     struct response_struct *r;  // Response params
{
  int i,j,k;
  int n,ti,si;
  float *s,cfx,x,y,x0;
  char tstr[SLEN];
  struct mm_proc_struct *p;
  struct mm_jel_struct *jel;

  exit_error("MM_MON_UPDATE_BY_PROC","Disabled");

  /*
  // Compute conversion factor from trial duration
  cfx = 1.8 / (r->tsn * 1000.0);  // GL_Coord / msec
  x0 = -0.90;

  y = 0.90;
  for(i=0;i<mm_nproc;i++){
    p = &(mm_ptable[i]);
    jel = p->done;
    for(j=0;j<p->ndone;j++){

      // Get trial number
      ti = jel->j->trial;
      si = jel->j->stim;
      s = r->s[0][ti];
      n = r->cnt[0][ti];

      y -= mm_line_h;                    // y-coord for this box
      for(k=0;k<n;k++){
	x = x0 + cfx * s[k];
	glBegin(GL_LINE_STRIP);
	glColor3f(1.0,1.0,1.0);  // Spikes are white
	glVertex2f(x,y);
	glVertex2f(x,y+mm_spike_h);
	glEnd();
      }
      sprintf(tstr,"%d",si);  // Stimulus number
      glplot_draw_text(mm_dpy2,mm_win2,mm_cx2,tstr,0.93,y,0.4,0.8,0.8);

      jel = jel->next;
    }
    if (p->done > 0)
      y -= mm_line_h;                 // leave a blank space for one box
  }
  */
}
/**************************************-**************************************/
/*                                                                           */
/*                             MM_MON_UPDATE_ALL                             */
/*                                                                           */
/*****************************************************************************/
void mm_mon_update_all(r)
     struct response_struct *r;  // Response params
{
  int i,j,k;
  struct mm_job_struct *js;

  if (mm_sort == 1){
    mm_mon_update_by_proc(r);
  }else{
    for(i=0;i<mm_nstim;i++){
      for(j=0;j<mm_nrpt;j++){
	js = &(mm_jtable[mm_mon_mi][i][j]);  // For currently displayed model

	if (js->tn > -1){  // Job has finished
	  mm_mon_update_draw_one(r,js->proc->id,js->trial,js->stim,js->rpt);
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MM_MON_HANDLE_ONE_EVENT                         */
/*                                                                           */
/*****************************************************************************/
/*
int mm_mon_handle_one_event(event)
     XEvent event;
{
  int redrawflag,n;
  XEvent evt;
  XKeyEvent *kevent;
  char kbuf[10];

  redrawflag = 0;
  switch(event.type){
  case ButtonPress:
    switch(event.xbutton.button){
    case 1:
      printf("  left\n");
      break;
    case 2:
      printf("  middle\n");
      break;
    case 3:
      printf("  right\n");
      break;
    }
    break;
  case MotionNotify:
    break;
  case KeyPress:
    kevent = (XKeyEvent *)&event;
    XLookupString(kevent,kbuf,sizeof(kbuf),NULL,NULL);
    
    if ((kbuf[0]=='h')||(kbuf[0]=='?')){
      mm_mon_help();
    }else if (kbuf[0]=='p'){
      mm_sort = 1; // Sort by processor
      mm_mon_set_lines_height();  // Update global drawing params
      redrawflag = 1;
    }else if (kbuf[0]=='s'){
      mm_sort = 0; // Sort by stimulus
      mm_mon_set_lines_height();  // Update global drawing params
      redrawflag = 1;
    }else if (kbuf[0]=='>'){
      mm_mon_mi += 1;
      if (mm_mon_mi >= mm_nmod)
	mm_mon_mi = 0;
      //printf("mm_mon_mi = %d\n",mm_mon_mi);
      redrawflag = 1;
    }else if (kbuf[0]=='<'){
      mm_mon_mi -= 1;
      if (mm_mon_mi < 0)
	mm_mon_mi = mm_nmod-1;
      //printf("mm_mon_mi = %d\n",mm_mon_mi);
      redrawflag = 1;
    }
    break;
  case ConfigureNotify:
    // Delaying here does *not* help us catch more events.
    
    // Clear accumulated events off of queue
    // If more events, ignore all but the last ConfigureNotify
    n = XQLength(mm_dpy2);
    while (n > 0){
      XPeekEvent(mm_dpy2,&evt); // Blocks until event received
      if (evt.type == ConfigureNotify){
	XNextEvent(mm_dpy2,&event); // Ignore all but last of these
      }else
	XNextEvent(mm_dpy2,&evt);   // Ignore these
      n -= 1;
    }

    // This avoids reconfig for Apple move events, which don't change size
    if ((event.xconfigure.width != mm_winw2) ||
	(event.xconfigure.height != mm_winh2)){
      
      glViewport(0,0,event.xconfigure.width,event.xconfigure.height);
      
      mm_winw2 = event.xconfigure.width;
      mm_winh2 = event.xconfigure.height;
      
      redrawflag = 1;
    }
  case Expose:  // On Apple, these occur as window is moved
    break;
  }
  
  return redrawflag;
}
*/
/**************************************-**************************************/
/*                                                                           */
/*                                MM_MON_EVENT                               */
/*                                                                           */
/*****************************************************************************/
void mm_mon_event(r)
     struct response_struct *r;  // Response params
{
  int redrawflag,tflag;
//  XEvent event;  WYETH 2019 comment GL/X dependency

  exit_error("MM_MON_EVENT","Disabled");

  /*

  glplot_util_make_current(mm_dpy2,mm_win2,mm_cx2);
  glDrawBuffer(GL_FRONT_AND_BACK); // Draw both front and back buffers
  //
  //  BEGIN MON WINDOW (i.e., drawing in 'mon' window now OK)
  //
  
  redrawflag = 0;
  while (XCheckIfEvent(mm_dpy2,&event,mm_true_test,(char *)NULL)){
    tflag = mm_mon_handle_one_event(event);
    if (tflag == 1){  // We got a CONFIG, must redraw, but wait first
      glplot_clear_double(mm_dpy2,mm_win2);
      usleep(50000);  // wait 50 msec, for other config events
      redrawflag = 1; // In the end, we must redraw for CONFIG event
    }
  }
  
  if (redrawflag){
    mm_mon_init(r);
    mm_mon_update_all(r);
  }
  
  //
  //  END MON WINDOW
  //
  glXSwapBuffers(mm_dpy2,mm_win2);
  glplot_util_make_current(mm_dpy,mm_win,mm_cx); // Switch back to GUI window

  */
}
/**************************************-**************************************/
/*                                                                           */
/*                              MM_HMS_FOR_TIME                              */
/*                                                                           */
/*  Get the string of the form "hr:mn:sc" for the time 't' in milliseconds.  */
/*                                                                           */
/*****************************************************************************/
char *mm_hms_for_time(t)
     int t; /* Time in ms */
{
  int hr,min,sec,s,m;
  char tstr[SLEN],*pstr;

  s = t/1000;     // Total seconds
  m = s/60;       // Total minutes
  hr = m/60;      // Total hours

  min = m - hr*60;
  sec = s - hr*3600 - min*60;

  sprintf(tstr,"%2d:%2d:%2d",hr,min,sec);

  pstr = strdup(tstr);

  return pstr;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MM_SUSPEND_PROCESSOR                           */
/*                                                                           */
/*  Suspend the processor, free any current job, set state.                  */
/*                                                                           */
/*****************************************************************************/
void mm_suspend_processor(pnum,targ_state)
     int pnum;   /* Proc number */
     int targ_state;  /* State to be assigned to processor (should be 2) */
{
  struct mm_proc_struct *pp;

  if ((pnum <= 0) || (pnum >= mm_nproc))
    exit_error("MM_SUSPEND_PROCESSOR","Invalid processor number");

  if (targ_state == 2){
    pp = &(mm_ptable[pnum]);
    if ((pp->state == 1) || (pp->state == 0)){
      pp->state = 2; /* user suspended */
      printf("    Processor %d suspended by user.\n",pnum);
      
      if (pp->job != NULL){
	pp->job->proc = NULL;  /* Open up this job for reassignment */
	pp->job = NULL;        /* No more job pending for this processor */
	mm_ndoa -= 1;          /* Update number of jobs assigned or done */
      }
    }else{
      printf("    Cannot suspend processor %d, not currently active.\n",pnum);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               MM_GUI_UPDATE                               */
/*                                                                           */
/*****************************************************************************/
void mm_gui_update()
{
  int i;
  unsigned long tcurr;
  int np,ntot;
  float xname,xndone,xndonet,xntot,xtlast,xcurr,xtime,xtarg,xsr;
  float y,yinc,ytit0,ytit1,ybot0,ybot1,ybot2,ybot3,ytarg;
  float c1r,c1g,c1b,c2r,c2g,c2b;
  float ftop,cf,cfh;
  char tstr[SLEN],*pstr;
  struct mm_proc_struct *pp;
  struct tms buffer;

  exit_error("MM_GUI_UPDATE","Disabled");

  /*
  glplot_util_make_current(mm_dpy,mm_win,mm_cx);

  ftop = (float)mm_winh;
  cf = 2.0 / ftop; // GL-Coord / pixel

  cfh = 2.0 / (float)mm_winw; // GL-Coord / pixel

  np = mm_nproc;

  xname   = -0.80;  // x-coord for name
  xsr     = -0.55;
  xndonet = -0.20;
  xtarg   = -0.11;
  xndone  =  0.10;
  xntot   =  0.12;
  xtlast  =  0.40;
  xcurr   =  0.60;  // x-coord for current job info
  xtime   =  0.65;

  ytit0 = -1.0 + cf * (ftop - 33.0);  // 33 pixels down from top
  ytit1 = -1.0 + cf * (ftop - 50.0);  // 50 pixels down from top
  y     = -1.0 + cf * (ftop - 75.0);  // 75 pixels down, 0th processor
  yinc  = cf * mm_pixel_proc_height;  // pixels per process


  ybot0 = -0.40;
  ybot1 = -0.50;
  ytarg = -0.60;

  ybot0 = -1.0 + cf * 110.0;  // 30 pixels up from bottom
  ybot1 = -1.0 + cf *  90.0;  // 30 pixels up from bottom
  ytarg = -1.0 + cf *  70.0;  // 30 pixels up from bottom
  ybot2 = -1.0 + cf *  50.0;  // 30 pixels up from bottom
  ybot3 = -1.0 + cf *  30.0;  // 30 pixels up from bottom


  c1r = 0.70;  c1g = 0.70;  c1b = 1.00;   // Titles, light blue
  c2r = 1.00;  c2g = 1.00;  c2b = 0.00;   // Table,  yellow

  glplot_clear(mm_dpy,mm_win);

  // Processor name
  sprintf(tstr,"Processors");
  glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xname,ytit1,c1r,c1g,c1b);

  sprintf(tstr,"Jobs");
  glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xndone,ytit0,c1r,c1g,c1b);
  sprintf(tstr,"done");
  glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xndone,ytit1,c1r,c1g,c1b);

  sprintf(tstr,"Last");
  glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xtlast,ytit0,c1r,c1g,c1b);
  sprintf(tstr,"time");
  glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xtlast,ytit1,c1r,c1g,c1b);

  sprintf(tstr,"Current");
  glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xcurr,ytit0,c1r,c1g,c1b);
  sprintf(tstr,"Mod  Stim  Rpt");
  glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xcurr,ytit1,c1r,c1g,c1b);

  ntot = 0;
  for(i=0;i<np;i++){
    pp = &(mm_ptable[i]);

    // Processor name
    sprintf(tstr,"%2d.   %s",i,pp->name);
    if (mm_ptable[i].state == 2) // Suspended
      glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xname,y,1.0,0.0,0.0);
    else if (mm_ptable[i].state == 3) // Suspended, but could reactivate
      glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xname,y,0.5,0.5,0.5);
    else
      glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xname,y,c2r,c2g,c2b);

    if (i == mm_proc_current) // Current
      glplot_draw_disk(0.7,0.7,1.0,xname-cf*10.0,y,cfh*10.0,cf*10.0);


    // Number of jobs done
    sprintf(tstr,"%3d",pp->ndone);
    glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xndone,y,c2r,c2g,c2b);
    ntot += pp->ndone;

    // Time taken on last job
    if (pp->ndone > 0){
      pstr = mm_hms_for_time(pp->donelast->j->tn);
      glplot_draw_text(mm_dpy,mm_win,mm_cx,pstr,xtlast,y,c2r,c2g,c2b);
      myfree(pstr);
    }else
      glplot_draw_text(mm_dpy,mm_win,mm_cx,"  :  :  ",xtlast,y,c2r,c2g,c2b);

    // Write current job
    if (pp->job == NULL){
      sprintf(tstr,"    -     -     -");
    }else
      sprintf(tstr,"%4d  %4d  %4d",pp->job->mi,pp->job->stim,pp->job->rpt);
    glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xcurr,y,c2r,c2g,c2b);

    y -= yinc;
  }

  sprintf(tstr,"-----");
  glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xndone,ybot0,c1r,c1g,c1b);

  sprintf(tstr,"Total done:");
  glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xndonet,ybot1,c1r,c1g,c1b);
  sprintf(tstr,"%d",ntot);
  glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xntot,ybot1,0.2,1.0,0.2);

  sprintf(tstr,"Target:");
  glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xtarg,ytarg,c1r,c1g,c1b);
  sprintf(tstr,"%d",mm_nstim * mm_nrpt);
  glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xntot,ytarg,c2r,c2g,c2b);

  tcurr = 10*times(&buffer);  // Get start time in ms
  pstr = mm_hms_for_time(tcurr - mm_t0);
  glplot_draw_text(mm_dpy,mm_win,mm_cx,pstr,xtime,ybot3,1.0,0.2,0.2);
  myfree(pstr);

  sprintf(tstr,"Stimuli");
  glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xname,ybot2,c1r,c1g,c1b);
  sprintf(tstr,"%d",mm_nstim);
  glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xsr,ybot2,1.0,1.0,1.0);

  sprintf(tstr,"Repeats");
  glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xname,ybot3,c1r,c1g,c1b);
  sprintf(tstr,"%d",mm_nrpt);
  glplot_draw_text(mm_dpy,mm_win,mm_cx,tstr,xsr,ybot3,1.0,1.0,1.0);

  glXSwapBuffers(mm_dpy,mm_win);
  */
}
/**************************************-**************************************/
/*                                                                           */
/*                                MM_GUI_INIT                                */
/*                                                                           */
/*****************************************************************************/
void mm_gui_init(r)
     struct response_struct *r;  // Response params
{
  int i;
  int ns;

  exit_error("MM_GUI_INIT","MM GUI disabled\n");

  /*

  mm_pixel_proc_top = 75.0;     // Store global value (pixels)
  mm_pixel_proc_height = 10.0;  // Store global value (pixels)

  ns = r->ns;  // Number of spike responses to save

  mm_winh = mm_pixel_proc_height * mm_nproc; // Win height det'd by # of procs
  mm_winh += 60;           // Top margin (pixels)
  mm_winh += 120;          // Bottom margin (pixels)

  if (mm_winh < 200)
    mm_winh = 200;
  
  glplot_init(mm_winw,mm_winh,mm_zoom,"mm",&mm_dpy,&mm_win,&mm_cx,&mm_dbflag,
	      0); // 0 - don't look at motion events
  glplot_clear_double(mm_dpy,mm_win);
  //glplot_set_xfont_generic(mm_dpy,mm_win,mm_cx,"medium");
  glplot_set_xfont_generic(mm_dpy,mm_win,mm_cx,"small");

  mm_proc_current = -1;  // No currently selected processor

  // mm_sort is already set
  mm_mon_set_lines_height();  // Update global drawing params

  // Second window to monitor spike output
  if (ns > 0){
    mm_monflag = 1; // Spike monitor
    glplot_init(mm_winw2,mm_winh2,mm_zoom,"mm_mon",&mm_dpy2,&mm_win2,&mm_cx2,
		&mm_dbflag2,0); // 0 - don't look at motion events
    glplot_clear_double(mm_dpy2,mm_win2);
    glplot_set_xfont_generic(mm_dpy2,mm_win2,mm_cx2,"medium");
    glplot_util_make_current(mm_dpy,mm_win,mm_cx);

    // Find index of response to monitor
    mm_moni = -1;
    i = 0;
    while((mm_moni < 0) && (i < r->n)){
      if (strcmp(r->rtype[i],"s")==0){
	mm_moni = i;
      }else{
	i += 1;
      }
    }
    if (mm_moni < 0)
      exit_error("MM_GUI_INIT","Failed to find response index");

    //  Initialize the monitor window
    glplot_util_make_current(mm_dpy2,mm_win2,mm_cx2);
    glDrawBuffer(GL_FRONT_AND_BACK); // Draw both front and back buffers
    //
    //  BEGIN MON WINDOW
    //

    mm_mon_init(r);

    //
    //  END MON WINDOW
    //
    glXSwapBuffers(mm_dpy2,mm_win2);
    glplot_util_make_current(mm_dpy,mm_win,mm_cx); // Switch back to GUI window
    
  }else{
    exit_error("MM_GUI_INIT","No spike responses requested");
    // WYETH - this creates a problem for updating the window
  }


  mm_gui_update();
  */
}
/**************************************-**************************************/
/*                                                                           */
/*                                MM_GUI_HELP                                */
/*                                                                           */
/*****************************************************************************/
void mm_gui_help()
{
  printf("\n");
  printf("  MM_GUI_HELP\n");
  printf("\n");
  printf("  COMMANDS:\n");
  printf("    h - help (print this list)\n");
  printf("    H - hold display after all responses are done\n");
  printf("    LEFT CLICK - select processor as current\n");
  printf("    s - suspend current processor (reassign any current job)\n");
  printf("    S - reactivate current processor (may accept new jobs)\n");
  printf("    q - exit, writing any response\n");
  printf("    ? - help (print this list)\n");
  printf("\n");
  printf("  COLORS:\n");
  printf("    yellow - processor is active or idle\n");
  printf("    red    - processor has been suspended by user.\n");
  printf("    gray   - user-suspended processor can be reactivated.\n");
  printf("\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                                MM_GUI_EVENT                               */
/*                                                                           */
/*****************************************************************************/
void mm_gui_event()
{
  int pi;
//  XEvent event;             *** WYETH 2019 comment out GL/X dependency
//  XKeyEvent *kevent;
//  XWindowAttributes xwa;
  char kbuf[10];

  exit_error("MM_GUI_EVENT","Disabled");

  /*
  glplot_util_make_current(mm_dpy,mm_win,mm_cx);

  if (XCheckIfEvent(mm_dpy,&event,mm_true_test,(char *)NULL)){
    switch(event.type){
    case ButtonPress:
      switch(event.xbutton.button){
      case 1:
	// Calculate processor index
	pi = my_rint(((float)event.xbutton.y - mm_pixel_proc_top)/
		     mm_pixel_proc_height);
	if ((pi == mm_proc_current) || (pi < 0) || (pi > mm_nproc))
	  mm_proc_current = -1;
	else
	  mm_proc_current = pi;
	break;
      case 2:
	break;
      case 3:
	break;
      }
      break;
    case MotionNotify:
      break;
    case KeyPress:
      kevent = (XKeyEvent *)&event;
      XLookupString(kevent,kbuf,sizeof(kbuf),NULL,NULL);
      
      if ((kbuf[0]=='h')||(kbuf[0]=='?')){
	mm_gui_help();
      }else if (kbuf[0]=='H'){
	mm_gui_hold_flag = 1;  // Hold at end
      }else if (kbuf[0]=='q'){
	mm_gui_hold_flag = 0;  // Stop holding
      }else if (kbuf[0]=='s'){
	if ((mm_proc_current > 0) && (mm_proc_current < mm_nproc))
	  mm_suspend_processor(mm_proc_current,2);
      }else if (kbuf[0]=='S'){
	if ((mm_proc_current > 0) && (mm_proc_current < mm_nproc)){
	  if (mm_ptable[mm_proc_current].state == 3){
	    mm_ptable[mm_proc_current].state = 0; // re-activate
	    printf("    Processor %d re-activated by user.\n",mm_proc_current);
	  }else if (mm_ptable[mm_proc_current].state == 2){
	    printf("    Cannot reactivate processor %d (no response yet).\n",
		   mm_proc_current);
	  }else{
	    printf("    Processor %d is already active.\n",mm_proc_current);
	  }
	}
      }
      break;
    case ConfigureNotify:
      break;
    case Expose:  // On Apple, these occur as window is moved
      break;
    }
  }
  */
}
/**************************************-**************************************/
/*                                                                           */
/*                                MM_GUI_HOLD                                */
/*                                                                           */
/*****************************************************************************/
void mm_gui_hold(r,hflag)
     struct response_struct *r;  // Response params
     int hflag;
{
  if (hflag != -1)
    mm_gui_hold_flag = hflag;

  //printf("MM_GUI_HOLD\n");

  // Wait for user to press 'q' to start, or 'q' to finish

  while (mm_gui_hold_flag == 1){
    mm_gui_update();
    mm_gui_event();
    mm_mon_event(r); // Check for events in response monitor window
    usleep(20000);   // So we don't run the processor all the time (20ms)
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MM_TIME_ELAPSED_MIN                           */
/*                                                                           */
/*  Return the number of minutes elapsed from t0 until 't' in minutes.       */
/*                                                                           */
/*****************************************************************************/
int mm_time_elapsed_min()
{
  unsigned long t;
  int s,m;
  struct tms buffer;

  t = 10*times(&buffer);         // Get current time in ms.
  s = ((int)(t - mm_t0))/1000;   // Total seconds elapsed since start
  m = s/60;                      // Total minutes elapsed since start

  return m;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MM_GET_PROC_STATE                             */
/*                                                                           */
/*****************************************************************************/
int mm_get_proc_state(pi)
     int pi; 
{
  if ((pi < 0) || (pi >= mm_nproc))
    exit_error("MM_GET_PROC_STATE","Bad processor ID");
  
  return mm_ptable[pi].state;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MM_PROC_ALLOW_REACTIVATE                        */
/*                                                                           */
/*****************************************************************************/
void mm_proc_allow_reactivate(pi)
     int pi; 
{
  int state;

  state = mm_get_proc_state(pi);
  if (state == 2)
    mm_ptable[pi].state = 3;
}
/**************************************-**************************************/
/*                                                                           */
/*                                MM_GET_NDONE                               */
/*                                                                           */
/*****************************************************************************/
int mm_get_ndone()
{
  return mm_ndone;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MM_PRINT_STAT                               */
/*                                                                           */
/*****************************************************************************/
void mm_print_stat()
{
  int i,j;
  int np;
  struct mm_proc_struct *pp;
  struct mm_jel_struct *jel;
  struct mm_job_struct *jp;

  np = mm_nproc;
  for(i=0;i<np;i++){
    pp = &(mm_ptable[i]);
    printf("  Proc %d did %d jobs (%s)\n",pp->id,pp->ndone,pp->name);
    struct mm_jel_struct *done;     // Head of list of completed jobs
    jel = pp->done;
    for(j=0;j<pp->ndone;j++){
      jp = jel->j;
      printf("    Mod Stm Rpt %4d %4d %4d   Time %.1f s\n",jp->mi,jp->stim,
	     jp->rpt,(float)jp->tn/1000.0);
      jel = jel->next;
    }
  }
  printf("  Total time = %.1f s\n",(float)(mm_t1-mm_t0)/1000.0);
}
/**************************************-**************************************/
/*                                                                           */
/*                                MM_JOB_DONE                                */
/*                                                                           */
/*  Update the job control data to reflect the completion of this job.       */
/*                                                                           */
/*****************************************************************************/
void mm_job_done(r,jnum,pnum,gui_flag,tlogf)
     struct response_struct *r;  // Response params
     int jnum;                   // Job number
     int pnum;                   // Proc number
     int gui_flag;
     char *tlogf;     //  If not NULL, keep track of assigned jobs
{
  int mi,si,ri;
  char tstr[SLEN];
  struct mm_job_struct *jp;
  struct mm_proc_struct *pp;
  struct mm_jel_struct *jel;
  struct tms buffer;

  mm_job_get_indices(jnum,&mi,&si,&ri);
  jp = mm_job_get_ptr(jnum);

  pp = &(mm_ptable[pnum]);

  jp->t2 = 10*times(&buffer);        // Get stop time in ms.
  mm_t1 = jp->t2;                    // Global finish time of last job
  jp->tn = (int)(jp->t2 - jp->t1);   // Compute run time for this job
  jp->trial = jnum;

  // Make a new job element
  jel = (struct mm_jel_struct *)myalloc(sizeof(struct mm_jel_struct));
  jel->j = jp;
  jel->next = NULL;

  // Add to list of done jobs for this proc
  if (pp->ndone == 0){
    pp->done = jel;
  }else{
    pp->donelast->next = jel;
  }
  pp->donelast = jel;
  pp->ndone += 1;    // Add to total number of jobs done

  pp->job = NULL;    // No job currently running
  pp->state = 0;     // Idle

  mm_ndone += 1;

  if (tlogf != NULL){
    sprintf(tstr,"  jnum %d proc %d took %d\n",jnum,pnum,jp->tn);
    append_string_to_file(tlogf,tstr);
  }

  if (gui_flag == 1){
    exit_error("MM_JOB_DONE","GUI Flag disabled");
    /*
    mm_gui_update();

    if (r->ns > 0){
      glplot_util_make_current(mm_dpy2,mm_win2,mm_cx2);
      glDrawBuffer(GL_FRONT_AND_BACK); // Draw both front and back buffers
      //
      //  BEGIN MON WINDOW
      //

      if (mm_sort == 1){
	mm_mon_init(r);
	mm_mon_update_by_proc(r);
      }else{
	mm_mon_update_draw_one(r,pnum,jnum,si,ri);
      }

      //
      //  END MON WINDOW
      //
      glXSwapBuffers(mm_dpy2,mm_win2);
      glplot_util_make_current(mm_dpy,mm_win,mm_cx); // Back to GUI window
    }
    */
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MM_CURRENT_TIME_MS                           */
/*                                                                           */
/*****************************************************************************/
unsigned long mm_current_time_ms()
{
  unsigned long ut;
  struct tms buffer;
  
  ut = 10*times(&buffer);  // Get start time in ms

  //printf("======>lu   %lu\n",ut);  // Get start time in ms

  return ut;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MM_ASSIGN_JOB                               */
/*                                                                           */
/*****************************************************************************/
void mm_assign_job(jnum,dest,mylogf,gui_flag,tlogf)
     int jnum,dest;
     char mylogf[];
     int gui_flag;
     char *tlogf;     //  If not NULL, keep track of assigned jobs
{
  //int m,s,r;
  struct mm_job_struct *jp;
  struct mm_proc_struct *pp;
  struct tms buffer;
  char tstr[SLEN];

  /*
    if (jnum > 0){
    sprintf(ggstr,"  Assigning P_%d job %d\n",dest,jnum);
    mylog(mylogf,ggstr);
    }*/


  if (tlogf != NULL){

    /*** MOOVAR - FIRST we might have to "re-prep" ***/
    /*** MOOVAR - FIRST we might have to "re-prep" ***/
    /*** MOOVAR Add 'prep' flag as param to enclosing routine ??  ***/

    MPI_Send(&jnum,1,MPI_INT,dest,tag_runtrial,MPI_COMM_WORLD);
  }

  // Update job control for this assignment
  if (jnum > -1){ // Assigning jnum -1 is a signal to terminate

    /*
    m = (int)(jnum/(mm_nstim * mm_nrpt));
    k = jnum - m*(mm_nstim * mm_nrpt);
    s = (int)(k/mm_nrpt);
    r = k - s*mm_nrpt;
    jp = &(mm_jtable[m][s][r]);
    */

    jp = mm_job_get_ptr(jnum);

    pp = &(mm_ptable[dest]);

    jp->proc = pp;
    jp->t1 = 10*times(&buffer);  // Get start time in ms

    pp->state = 1;  // Indicate that processor is now busy
    pp->job = jp;   // Current job

    mm_ndoa += 1;
  }
  if (tlogf != NULL){
    sprintf(tstr,"Assign job %d proc %d\n",jnum,dest);
    append_string_to_file(tlogf,tstr);
  }

  if (gui_flag == 1)
    mm_gui_update();
}
/**************************************-**************************************/
/*                                                                           */
/*                          MM_GET_START_TIME_FOR_JOB                        */
/*                                                                           */
/*****************************************************************************/
unsigned long mm_get_start_time_for_job(jnum)
     int jnum;
{
  unsigned long tms;
  struct mm_job_struct *jp;

  tms = -1;  // Time in msec to return

  if (jnum > -1){ // Assigning jnum -1 is a signal to terminate
    jp = mm_job_get_ptr(jnum);
    //printf("jp->t1 = %lu\n",jp->t1);
    tms = jp->t1;
  }

  return tms;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MM_JA_LOG_PROC_NAMES                          */
/*                                                                           */
/*****************************************************************************/
void mm_ja_log_proc_names(tlogf)
     char *tlogf;
{
  int i;
  char tstr[SLEN];

  if (tlogf != NULL){

    for(i=0;i<mm_nproc;i++){
      sprintf(tstr,"PROC %4d %s\n",i,mm_ptable[i].name);
      append_string_to_file(tlogf,tstr);
    }

    sprintf(tstr,"NSTIM %d\n",mm_nstim);
    append_string_to_file(tlogf,tstr);

    sprintf(tstr,"NRPT %d\n",mm_nrpt);
    append_string_to_file(tlogf,tstr);
    
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MM_INIT_JOB_CONTROL                           */
/*                                                                           */
/*****************************************************************************/
void mm_init_job_control(nproc,nmod,nstim,nrpt,mylogf,mastername,pnames)
     int nproc;     // Number of processors
     int nmod;      // Number of models
     int nstim;     // Number of stimuli
     int nrpt;      // Number of repeats
     char mylogf[];      // Logfile
     char mastername[];  // Name of master process
     char **pnames;      // NULL, or list of proc names, for replay
{
  int i,j,k;
  int nstr;
  char *pstr;
  struct tms buffer;

  mylog(mylogf,"  MM_INIT_JOB_CONTROL\n");

  // Initialize job table
  mm_ndone = mm_ndoa = 0;
  mm_nmod = nmod;
  mm_nstim = nstim;
  mm_nrpt = nrpt;
  mm_jtable = (struct mm_job_struct ***)myalloc(nmod*
				        sizeof(struct mm_job_struct **));
  for(k=0;k<nmod;k++){
    mm_jtable[k] = (struct mm_job_struct **)myalloc(nstim*
					    sizeof(struct mm_job_struct *));
    for(i=0;i<nstim;i++){
      mm_jtable[k][i] = (struct mm_job_struct *)myalloc(nrpt*
						sizeof(struct mm_job_struct));
      for(j=0;j<nrpt;j++){
	mm_jtable[k][i][j].mi = k;
	mm_jtable[k][i][j].stim = i;
	mm_jtable[k][i][j].rpt = j;
	mm_jtable[k][i][j].trial = -1;
	mm_jtable[k][i][j].proc = NULL;
	mm_jtable[k][i][j].t1 = -1;
	mm_jtable[k][i][j].t2 = -1;
	mm_jtable[k][i][j].tn = -1;
      }
    }
  }

  // Initialize processor table
  mm_nproc = nproc;
  mm_ptable = (struct mm_proc_struct *)myalloc(nproc*
					       sizeof(struct mm_proc_struct));
  for(i=0;i<nproc;i++){
    mm_ptable[i].name = NULL;
    mm_ptable[i].id = i;
    mm_ptable[i].state = -1;
    mm_ptable[i].job = NULL;
    mm_ptable[i].ndone = 0;
    mm_ptable[i].done = NULL;
    mm_ptable[i].donelast = NULL;
  }

  mm_ptable[0].name = strdup(mastername);

  // Get slave names
  if (pnames == NULL){
    for(i=1;i<nproc;i++){
      mm_mpi_receive_carray(i,tag_procname,10000,mylogf,&pstr,&nstr);
      mm_ptable[i].name = strdup(pstr);
      mm_ptable[i].state = 0; // Indicate that the processor is available
      myfree(pstr);
    }
  }else{
    for(i=1;i<nproc;i++){
      mm_ptable[i].name = strdup(pnames[i]);
      mm_ptable[i].state = 0; // Indicate that the processor is available
    }
  }

  mm_t0 = 10*times(&buffer);  // Get start time in ms.
}
/**************************************-**************************************/
/*                                                                           */
/*                          MM_GET_FASTEST_FREE_PROC                         */
/*                                                                           */
/*  Return the number of the fastest free processor, or -1 if none free.     */
/*                                                                           */
/*****************************************************************************/
int mm_get_fastest_free_proc()
{
  int i;
  int pnum;
  float tmin;
  struct mm_proc_struct *p;

  pnum = -1;
  for(i=1;i<mm_nproc;i++){ // 0th processor is master
    p = &(mm_ptable[i]);  // Processor i
    if (p->state == 0){
      if (pnum < 0){ // This is the first available processor found
	pnum = i;
	if (p->ndone > 0) // If this proc has done a job in the past
	  tmin = p->donelast->j->tn; // time taken on last job by this proc
	else
	  tmin = -1;      // No history for this proc
      }else{
	if (p->ndone > 0){
	  if (tmin > p->donelast->j->tn){ // This proc was faster
	    pnum = i;
	    tmin = p->donelast->j->tn;
	  }
	}else{
	  if (tmin > -1){
	    pnum = i;      // Use this proc, which is untested
	    tmin = -1;
	  }
	}
      }
    }
  }
  return pnum;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MM_GET_NEXT_AVAIL_RPT                           */
/*                                                                           */
/*  Return the next repeat number to be run for this 'stim', or -1 if none.  */
/*                                                                           */
/*****************************************************************************/
int mm_get_next_avail_rpt(mi,stim)
     int mi;
     int stim;
{
  int i;
  int r,done;
  struct mm_job_struct *jp;

  r = -1;

  done = 0;
  i = 0;
  while(!done && (i < mm_nrpt)){
    jp = &(mm_jtable[mi][stim][i]);
    if (jp->proc == NULL){
      done = 1;
      r = i;
    }else
      i += 1;
  }

  return r;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MM_GET_NEXT_JOB_FOR_PROC                         */
/*                                                                           */
/*  Return the pointer to the next job to be done on processor 'pnum'.       */
/*                                                                           */
/*****************************************************************************/
struct mm_job_struct *mm_get_next_job_for_proc(pnum,tlogf)
     int pnum;
     char *tlogf;
{
  int m,s,r;
  int done;
  char tstr[SLEN];
  struct mm_proc_struct *p;
  struct mm_job_struct *jp;

  p = &(mm_ptable[pnum]);   // Process pointer
  jp = NULL;                // Job pointer
  
  // If this proc is working on a series, assign next available repeat
  if (p->donelast != NULL){

    // *** WYETH DEBUG REMOVE
    if (tlogf != NULL){
      sprintf(tstr,"    pnum %d  trying method 1\n",pnum);
      append_string_to_file(tlogf,tstr);
    }

    m = p->donelast->j->mi;
    s = p->donelast->j->stim;
    r = mm_get_next_avail_rpt(m,s);
    if (r >= 0){
      jp = &(mm_jtable[m][s][r]);
    }
  }

  // Otherwise, try to find a new stimulus series
  if (jp == NULL){

    // *** WYETH DEBUG REMOVE
    if (tlogf != NULL){
      sprintf(tstr,"    pnum %d  trying method 2\n",pnum);
      append_string_to_file(tlogf,tstr);
    }

// NSTIM 19
// NRPT 20
//    pnum 1  trying method 2
// Assign job 0 proc 1
//    pnum 2  trying method 2
//    pnum 2  trying method 3
// Assign job 1 proc 2
//    pnum 3  trying method 2
//    pnum 3  trying method 3
//
//  *** WYETH HERE - does the next line make sense???  March 2018
//  *** WYETH HERE - does the next line make sense???
//  *** WYETH HERE - does the next line make sense???
//  *** WYETH HERE - does the next line make sense???
//  *** WYETH HERE - does the next line make sense???
//  Even if we have not done one before, we want to find one in a new
//     series ???

    if (p->donelast != NULL){  // Try to use a stimulus from the same model
      m = p->donelast->j->mi;

      done = 0;
      s = 0;
      while(!done && (s < mm_nstim)){
	if (mm_jtable[m][s][0].proc == NULL){  // Found a new series
	  jp = &(mm_jtable[m][s][0]);
	  done = 1;
	}else{
	  s += 1;
	}
      }
    }else{  // Try to start a new model config

      done = 0;
      m = 0;
      while(!done && (m < mm_nmod)){
	if (mm_jtable[m][0][0].proc == NULL){  // Found a new model
	  jp = &(mm_jtable[m][0][0]);
	  done = 1;
	}else{
	  m += 1;
	}
      }
    }
  }

  // Otherwise, take the next available stimulus repeat
  if (jp == NULL){

    // *** WYETH DEBUG REMOVE
    if (tlogf != NULL){
      sprintf(tstr,"    pnum %d  trying method 3\n",pnum);
      append_string_to_file(tlogf,tstr);
    }

    done = 0;
    m = s = r = 0;
    while(!done){
      if (mm_jtable[m][s][r].proc == NULL){  // Found a job
	jp = &(mm_jtable[m][s][r]);
	done = 1;
      }else{
	r += 1; // Increment repeat number
	if (r >= mm_nrpt){
	  r = 0;
	  s += 1; // Go to next stimulus
	  if (s >= mm_nstim){
	    s = 0;
	    m += 1;
	    if (m >= mm_nmod){
	      printf("  mm_ndone = %d\n",mm_ndone);
	      printf("  mm_ndoa = %d\n",mm_ndoa);
	      exit_error("MM_GET_NEXT_JOB_FOR_PROC","This should not happen");
	    }
	  }
	}
      }
    }
  }

  /*** WYETH MOOVAR ****/
  // Check if the 'mi' for 'jp' is not equal to last done, and if so, set
  // a flag to say 'new PREP needed'
  //


  return jp;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MM_SUGGEST_NEXT_ASSIGNMENT                        */
/*                                                                           */
/*  Return  -2 if all jobs are done                                          */
/*          -1 if there are jobs running, but none left to be assigned       */
/*           0 if there are jobs to assign, but all processors are busy      */
/*           1 if a valid assignment is suggested                            */
/* MOOVAR ?? 2 if a valid assignment for a NEW MODEL is suggested            */
/*                                                                           */
/*****************************************************************************/
int mm_suggest_next_assignment(rp,rj,tlogf)
     int *rp;   // Processor to assign
     int *rj;   // Job index to assign, over all mod, stim, rpts.
     char *tlogf;
{
  int i;
  int p,j,flag;
  struct mm_job_struct *jp;
  struct mm_proc_struct *pp;
  time_t rawtime; 
  struct tm *timeinfo;

  // Check for suspended processors, and reassign jobs
  for(i=0;i<mm_nproc;i++){
    pp = &(mm_ptable[i]);
    if ((pp->state == 2) && (pp->job != NULL)){
      pp->job->proc = NULL;  // Open up this job for reassignment
      pp->job = NULL;        // No more job pending for this processor
      mm_ndoa -= 1;          // Update number of jobs assigned or done
    }
  }

  time(&rawtime);
  timeinfo = localtime(&rawtime);
  //printf("===> Current date and time are: %s",asctime(timeinfo));

  if (mm_ndone == (mm_nmod * mm_nstim * mm_nrpt))
    flag = -2;
  else if (mm_ndoa == (mm_nmod * mm_nstim * mm_nrpt))
    flag = -1;
  else{
    p = mm_get_fastest_free_proc();
    if (p <= 0)
      flag = 0;
    else{
      flag = 1;
      // Find best job to run on this processor
      jp = mm_get_next_job_for_proc(p,tlogf);
      if (jp == NULL){
	printf("MM_SUGGEST_NEXT_ASSIGNMENT:  THIS SHOULD NOT HAPPEN\n");
	exit(0);
      }else{
	//j = jp->stim * mm_nrpt + jp->rpt; // JobNo = stim * nrpt + rpt
	//
	//  WYETH - is this equivalent to 'r->gtsi' in 'wm' ?
	//  Yes: mi * s->ntr              +   r->tsi (inc's over stim, rpt)
	j = jp->mi * (mm_nstim * mm_nrpt) + jp->stim * mm_nrpt + jp->rpt;
      }
    }
  }
  
  *rp = p; *rj = j;

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MM_FIT_ASSIGN_JOB                            */
/*                                                                           */
/*****************************************************************************/
void mm_fit_assign_job(jnum,dest,mylogf,gui_flag,tlogf)
     int jnum,dest;
     char mylogf[];
     int gui_flag;
     char *tlogf;     //  If not NULL, keep track of assigned jobs
{
  int ccode;
  //int m,s,r;
  struct mm_job_struct *jp;
  struct mm_proc_struct *pp;
  struct tms buffer;
  char tstr[SLEN];

  /*

  //
  //  Send a command value to say "receive data"
  //
  //  *** WYETH - Why not send a string, that is a comman?
  //
  ccode = 2;  // Code to transmit a array of ints
  MPI_Send(&ccode,1,MPI_INT,dest,tag_fit_cmd,MPI_COMM_WORLD);

  //
  //  Send a string that describes the data,
  //  Then send the array of ints
  //
  strcpy(tstr,"ci");
  MPI_Send(tstr,strlen(tstr)+1,MPI_CHAR,dest,tag_fit_str,MPI_COMM_WORLD);
  MPI_Send(ci,cin,MPI_INT,dest,tag_fit_dvi,MPI_COMM_WORLD);

  */
}
