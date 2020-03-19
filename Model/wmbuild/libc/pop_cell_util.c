/*****************************************************************************/
/*                                                                           */
/*  pop_cell_util.c                                                          */
/*  wyeth bair                                                               */
/*                                                                           */
/*  Low level utilities to operate on 'pop_cell' units.                      */
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
#include "carray_util.h"
#include "iarray_util.h"
#include "farray_util.h"
#include "data_util.h"
#include "spike_util.h"
#include "paramfile_util.h"
#include "paramfile.h"
#include "mod_util.h"
#include "ifc.h" // Data structures

/**************************************-**************************************/
/*                                                                           */
/*                          POP_CELL_PARSE_UNIT_NAME                         */
/*                                                                           */
/*  'name' is of the form:  layer_x_y_z                                      */
/*                                                                           */
/*****************************************************************************/
void pop_cell_parse_unit_name(name,rlay,rx,ry,rz)
     char *name;
     char **rlay;
     int *rx,*ry,*rz;
{
  int i,k;
  int x,y,z,n,cnt;
  char lay[SLEN],tstr[SLEN];

  strcpy(tstr,name);
  n = strlen(name);

  cnt = 0;
  for(i=0;i<n;i++){
    if (tstr[i] == '_')
      cnt += 1;
  }

  // WYETH.EDIT.2010.Aug.26
  if (cnt < 3){
    printf("  Offending name:  %s\n",name);
    exit_error("POP_CELL_PARSE_UNIT_NAME","expecting 3 '_'");
  }else{
    k = 0;
    for(i=0;i<n;i++){
      if (tstr[i] == '_'){
	if (k >= (cnt-3)){
	  tstr[i] = ' ';
	}
	k += 1;
      }
    }
  }
      
  sscanf(tstr,"%s %d %d %d",lay,&x,&y,&z);

  *rx = x;
  *ry = y;
  *rz = z;
  *rlay = strdup(lay);
}
/**************************************-**************************************/
/*                                                                           */
/*                          POP_CELL_GET_LAYER_INDEX                         */
/*                                                                           */
/*  Return the index of the named layer from the laylist.                    */
/*                                                                           */
/*****************************************************************************/
int pop_cell_get_layer_index(laylist,n,name)
     struct pop_layer **laylist;
     int n;
     char *name;
{
  int i,k;

  k = -1;
  for(i=0;i<n;i++){
    if (strcmp(name,laylist[i]->name)==0)
      k = i;
  }
  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                         POP_CELL_GET_LAYER_POINTER                        */
/*                                                                           */
/*  Return the pointer to the named layer from the laylist.                  */
/*                                                                           */
/*****************************************************************************/
struct pop_layer *pop_cell_get_layer_pointer(laylist,n,name)
     struct pop_layer **laylist;
     int n;
     char *name;
{
  int k;

  k = pop_cell_get_layer_index(laylist,n,name);
  if (k == -1)
    return NULL;
  else
    return laylist[k];
}
/**************************************-**************************************/
/*                                                                           */
/*                     POP_CELL_LAYLIST_GET_CELL_POINTER                     */
/*                                                                           */
/*  Return the pointer to the cell from the named layer for coords.          */
/*                                                                           */
/*****************************************************************************/
struct pop_cell *pop_cell_laylist_get_cell_pointer(laylist,n,name,xi,yi,zi)
     struct pop_layer **laylist;
     int n;
     char *name;
     int xi,yi,zi;
{
  int k;
  struct pop_cell *c;
  struct pop_layer *pl;

  k = pop_cell_get_layer_index(laylist,n,name);
  if (k == -1){
    printf("name = %s\n",name);
    exit_error("POP_CELL_LAYLIST_GET_CELL_POINTER","No such layer");
  }

  pl = laylist[k];

  if ((xi<0)||(xi>=pl->xn)||(yi<0)||(yi>=pl->yn)||(zi<0)||(zi>=pl->zn)){
    printf("  *** Requested:  %d %d %d\n",xi,yi,zi);
    printf("  *** Limits:  0..%d 0..%d 0..%d\n",pl->xn-1,pl->yn-1,pl->zn-1);
    exit_error("POP_CELL_LAYLIST_GET_CELL_POINTER","Cell index error");
  }

  c = &(pl->c[xi][yi][zi]);

  return c;
}
/**************************************-**************************************/
/*                                                                           */
/*                          POP_CELL_GET_INIT_CIRCBUF                        */
/*                                                                           */
/*****************************************************************************/
struct pop_circbuf *pop_cell_get_init_circbuf(n,t,dt,zflag)
     int n;
     double t,dt;
     int zflag;
{
  struct pop_circbuf *cb;

  cb = (struct pop_circbuf *)myalloc(sizeof(struct pop_circbuf));
  cb->n = n;
  if (zflag)
    cb->d = get_zero_farray(n);
  else
    cb->d = (float *)myalloc(n*sizeof(float));
  cb->i = 0;
  cb->t = t;
  cb->dt = dt;
  cb->maskn = -1;

  return cb;
}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_CELL_ZERO_CIRCBUF                          */
/*                                                                           */
/*****************************************************************************/
void pop_cell_zero_circbuf(cb)
     struct pop_circbuf *cb;
{
  int i;
  int n;
  float *d;

  n = cb->n;
  d = cb->d;
  for(i=0;i<n;i++)
    d[i] = 0.0;         // WYETH BUG FOUND 23 Apr 2009, was 'd[0]'
  cb->i = 0;
  cb->t = 0.0;    // See also, ..._TM1 below
}
/**************************************-**************************************/
/*                                                                           */
/*                          POP_CELL_ZERO_CIRCBUF_TM1                        */
/*                                                                           */
/*****************************************************************************/
void pop_cell_zero_circbuf_tm1(cb,t0)
     struct pop_circbuf *cb;
     double t0;
{
  int i;
  int n;
  float *d;

  n = cb->n;
  d = cb->d;
  for(i=0;i<n;i++)
    d[i] = 0.0;
  cb->i = 0;
  cb->t = t0;   // "tmi":  t = -1
}
/**************************************-**************************************/
/*                                                                           */
/*                              POP_CELL_SI_FREE                             */
/*                                                                           */
/*****************************************************************************/
void pop_cell_si_free(si)
     struct pop_si *si;
{

  // Currently don't free it, because the other part of the SI 
  // (other SI index synapses) might need it, if they haven't been freed
  
  //  WYETH - should write 
  //

  //printf("*** POP_CELL_SI_FREE  Warning, 'siu' not freed\n");

  myfree(si);  // Free the 'si' but not the 'siu' in case it is used.

  return;



  //
  //  WYETH - not finished...not tested  17 May 2009
  //

  if (si->siui == 1){
    struct pop_si001 *si001;

    si001 = si->siu->s001;

    // FREE THESE
    // si001->wt
    // si001->sisv

    // The 'mech' is probably shared, don't free?

  }else if (si->siui == 2){
    struct pop_si002 *si002;

    si002 = si->siu->s002;

    // FREE THESE
    // si001->wt1
    // si001->wt2
    // si001->sisv
 
    // The 'mech' is probably shared, don't free?

  }else{
    exit_error("POP_CELL_SI_FREE","Unknown SI index");
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_CELL_SI_SAV_CLEAR                          */
/*                                                                           */
/*****************************************************************************/
void pop_cell_si_sav_clear(sv)
     struct pop_si_sav *sv;
{
  int i,j;

  if (sv != NULL){
    for(i=0;i<sv->n;i++){  // For each stored trace
      for(j=0;j<sv->cnt[i];j++)
	sv->d[i][j] = 0.0;       // zero all data elements
      sv->cnt[i] = 0;            // zero count
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_CELL_SI_MAKE_001                           */
/*                                                                           */
/*****************************************************************************/
void pop_cell_si_make_001(s0,s1,msi,circn,circdt)
     struct pop_syn *s0,*s1; // synapse 0 (controls storage) and 1
     struct pop_mech *msi;   // Pointer to SI definition
     int circn;              // Size of circular buffer
     float circdt;           // (sec) time resolution of circ buffer
{
  struct pop_circbuf *cb;
  struct pop_si001 *psi001;

  // Synapse 0 should control 'si' storage
  s0->si = (struct pop_si *)myalloc(sizeof(struct pop_si));
  s0->si->ci = 0;
  s0->si->syn_code = 0;
  s0->si->siui = 1; // SI union type
  s0->si->siu = (union si_union *)myalloc(sizeof(union si_union));
  psi001 = (struct pop_si001 *)myalloc(sizeof(struct pop_si001));
  s0->si->siu->s001 = psi001;
  s0->si->siu->s001->msi = msi; // Index into layer list
  cb = (struct pop_circbuf *)myalloc(sizeof(struct pop_circbuf));
  s0->si->siu->s001->wt = cb;
  s0->si->siu->s001->sisv = NULL; // Indicates that nothing is to be saved

  cb->n = circn;
  cb->d = (float *)myalloc(circn*sizeof(float));
  cb->i = 0; // Index of origin
  cb->t = -1.0;
  cb->dt = circdt;
  cb->maskn = -1;

  //
  //  's1' will have it's own 'si', but share the 'siu'
  //
  s1->si = (struct pop_si *)myalloc(sizeof(struct pop_si));
  s1->si->ci = 1;
  s0->si->syn_code = 0;
  s1->si->siui = 1; // SI union type
  s1->si->siu = s0->si->siu;
}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_CELL_SI_RESET_001                          */
/*                                                                           */
/*****************************************************************************/
void pop_cell_si_reset_001(s)
     struct pop_si001 *s;
{
  int i,j;
  int n;
  struct pop_circbuf *cb;

  if (s->wt == NULL)
    exit_error("POP_CELL_SI_RESET_001","Weight v. time circbuf is NULL");

  // Clear the circbuf, and reset the index and time (time to -1.0)
  pop_cell_zero_circbuf_tm1(s->wt,-1.0);

  // Reset any response data storage
  if (s->sisv != NULL)
    pop_cell_si_sav_clear(s->sisv);

}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_CELL_SI_MAKE_002                           */
/*                                                                           */
/*  Two weight v. time buffers.                                              */
/*                                                                           */
/*****************************************************************************/
void pop_cell_si_make_002(s0,s1,msi,circn,circdt,t0,zflag)
     struct pop_syn *s0,*s1;
     struct pop_mech *msi;  // Pointer to SI definition
     int circn;             // Size of circular buffer
     float circdt;          // (sec) time resolution of circ buffer
     float t0;              // start time (or -1.0)
     int zflag;             // 1-zero the circbuf
{
  struct pop_circbuf *cb1,*cb2;
  struct pop_si002 *psi002;

  // Synapse 0 should control 'si' storage
  s0->si = (struct pop_si *)myalloc(sizeof(struct pop_si));
  s0->si->ci = 0;
  s0->si->syn_code = 0;
  s0->si->siui = 2; // SI union type
  s0->si->siu = (union si_union *)myalloc(sizeof(union si_union));
  psi002 = (struct pop_si002 *)myalloc(sizeof(struct pop_si002));
  s0->si->siu->s002 = psi002;
  s0->si->siu->s002->msi = msi; // Index into layer list

  //
  //  Create 2 circular buffers, one for each cell
  //
  cb1 = pop_cell_get_init_circbuf(circn,t0,circdt,zflag);
  s0->si->siu->s002->wt1 = cb1;
  cb1->maskn = msi->mask1n;
  cb2 = pop_cell_get_init_circbuf(circn,t0,circdt,zflag);
  s0->si->siu->s002->wt2 = cb2;
  cb2->maskn = msi->mask1n;

  s0->si->siu->s002->sisv = NULL; // Indicates that nothing is to be saved

  s1->si = (struct pop_si *)myalloc(sizeof(struct pop_si));
  s1->si->ci = 1;
  s0->si->syn_code = 0;
  s1->si->siui = 2; // SI union type
  s1->si->siu = s0->si->siu;
}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_CELL_SI_RESET_002                          */
/*                                                                           */
/*****************************************************************************/
void pop_cell_si_reset_002(s,t0)
     struct pop_si002 *s;
     double t0;
{
  int i,j;

  if (s->wt1 == NULL)
    exit_error("POP_CELL_SI_RESET_002","Weight v. time circbuf is NULL");
  if (s->wt2 == NULL)
    exit_error("POP_CELL_SI_RESET_002","Weight v. time circbuf is NULL");

  // Clear both circbuf's, and reset the index and time (time to -1.0)
  pop_cell_zero_circbuf_tm1(s->wt1,t0);
  pop_cell_zero_circbuf_tm1(s->wt2,t0);

  // Reset any response data storage
  if (s->sisv != NULL)
    pop_cell_si_sav_clear(s->sisv);
}
/**************************************-**************************************/
/*                                                                           */
/*                              POP_CELL_SI_RESET                            */
/*                                                                           */
/*****************************************************************************/
void pop_cell_si_reset(syn)
     struct pop_syn *syn;
{
  if (syn->si->siui == 1){
    pop_cell_si_reset_001(syn->si->siu->s001);
  }else if (syn->si->siui == 2){

    //printf("pre/post   %s  %s\n",syn->pre->name,syn->post->name);

    if (syn->si->syn_code == 100002){
      pop_cell_si_reset_002(syn->si->siu->s002,0.0);
    }else{
      pop_cell_si_reset_002(syn->si->siu->s002,-1.0);
    }

  }else{
    exit_error("POP_CELL_SI_RESET","Unknown SI Union Index");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_CELL_SI_SAV_SET_N                          */
/*                                                                           */
/*****************************************************************************/
void pop_cell_si_sav_set_n(mylogf,s,datid,n)
     char *mylogf;
     struct pop_syn *s;
     char *datid;
     int n;
{
  int i;
  struct pop_si_sav *sisv;

  /*mylog(mylogf,"  POP_CELL_SI_SAV_SET_N");*/

  if (s->si->siui == 1)
    sisv = s->si->siu->s001->sisv;
  else if (s->si->siui == 2)
    sisv = s->si->siu->s002->sisv;
  else
    mylog_exit(mylogf,"*** POP_CELL_SI_SAV_SET_N  Unknown SI index\n");

  if (sisv == NULL)
    mylog_exit(mylogf,"*** POP_CELL_SI_SAV_SET_N  sisv is NULL\n");

  i = search_2d_carray(sisv->name,datid,sisv->n);

  if (i < 0)
    mylog_exit(mylogf,"*** POP_CELL_SI_SAV_SET_N  Unknown SI name\n");

  sisv->cnt[i] = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                          POP_CELL_SI_GET_DATA_PTR                         */
/*                                                                           */
/*****************************************************************************/
float *pop_cell_si_get_data_ptr(mylogf,s,datid,rn)
     char *mylogf;
     struct pop_syn *s;
     char *datid;
     int *rn;
{
  int i;
  struct pop_si_sav *sisv;

  /*mylog(mylogf,"  POP_CELL_SI_GET_DATA_PTR");*/


  if (s->si->siui == 1)
    sisv = s->si->siu->s001->sisv;
  else if (s->si->siui == 2)
    sisv = s->si->siu->s002->sisv;
  else
    mylog_exit(mylogf,"*** POP_CELL_SI_GET_DATA_PTR  Unknown SI index\n");

  if (sisv == NULL)
    mylog_exit(mylogf,"*** POP_CELL_SI_GET_DATA_PTR  sisv is NULL\n");

  i = search_2d_carray(sisv->name,datid,sisv->n);

  if (i < 0){
    printf("  VALID CHOICES ARE:\n");
    for(i=0;i<sisv->n;i++)
      printf("    %2d  %s\n",i,sisv->name[i]);

    printf("  s->pre   %s\n",s->pre->name);
    printf("  s->post  %s\n",s->post->name);
    printf("  datid    %s\n",datid);
    mylog_exit(mylogf,"*** POP_CELL_SI_GET_DATA_PTR  Unknown SI name\n");
  }
  
  *rn = sisv->cnt[i];
  return sisv->d[i];
}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_CELL_SI_ADD_SAV                            */
/*                                                                           */
/*****************************************************************************/
void pop_cell_si_add_sav(mylogf,s,datid)
     char *mylogf;
     struct pop_syn *s;
     char *datid;
{
  int i;
  int n,*tcnt,*tnmax;
  float *tsamp,**td;
  char **tname;
  struct pop_si_sav *sisv;

  mylog(mylogf,"  POP_CELL_SI_ADD_SAVE\n");

  if (s->si->siui == 1)
    sisv = s->si->siu->s001->sisv;
  else if (s->si->siui == 2)
    sisv = s->si->siu->s002->sisv;
  else
    mylog_exit(mylogf,"*** POP_CELL_SI_ADD_SAVE  Unknown SI index\n");

  // Create new storage, set index 'i'
  if (sisv == NULL){
    sisv = (struct pop_si_sav *)myalloc(sizeof(struct pop_si_sav));
    n = 1;
    sisv->name = (char **)myalloc(sizeof(char *));
    sisv->nmax = (int *)myalloc(sizeof(int));
    sisv->cnt  = (int *)myalloc(sizeof(int));
    sisv->samp = (float *)myalloc(sizeof(float));
    sisv->d    = (float **)myalloc(sizeof(float *));
    i = 0;

    if (s->si->siui == 1)
      s->si->siu->s001->sisv = sisv;
    else if (s->si->siui == 2)
      s->si->siu->s002->sisv = sisv;
    
  }else{
    n = sisv->n + 1;
    tname = (char **)myalloc(n*sizeof(char *));
    tcnt  = (int *)myalloc(n*sizeof(int));
    tnmax  = (int *)myalloc(n*sizeof(int));
    tsamp = (float *)myalloc(n*sizeof(float));
    td    = (float **)myalloc(n*sizeof(float *));
    for(i=0;i<sisv->n;i++){
      tname[i] = sisv->name[i];
      tnmax[i] = sisv->nmax[i];
      tcnt[i]  = sisv->cnt[i];
      tsamp[i] = sisv->samp[i];
      td[i]    = sisv->d[i];
    }
    myfree(sisv->name);
    myfree(sisv->nmax);
    myfree(sisv->cnt);
    myfree(sisv->samp);
    myfree(sisv->d);
    sisv->name = tname;
    sisv->nmax = tnmax;
    sisv->cnt = tcnt;
    sisv->samp = tsamp;
    sisv->d = td;
    i = sisv->n;
  }

  sisv->name[i] = strdup(datid);
  sisv->nmax[i] = 2000000;  // WYETH - also see gtex1max below
  sisv->cnt[i] = 0;
  sisv->samp[i] = 1000.0;   // WYETH - warning, setting all sampling to msec
  sisv->d[i] = get_zero_farray(sisv->nmax[i]);  // Must have zero's

  sisv->n = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                              POPC_CSAV_GET_CSI                            */
/*                                                                           */
/*  Return the csav index if the 'datid' exists in the layer for the cell,   */
/*  or -1 otherwise.                                                         */
/*                                                                           */
/*****************************************************************************/
int popc_csav_get_csi(c,datid)
     struct pop_cell *c;
     char *datid;
{
  int csi;

  csi = -1;

  if (c->pl != NULL){
    if (c->pl->csav_n > 0){
      csi = search_2d_carray(c->pl->csav_datid,datid,c->pl->csav_n);
      // -1 if not found
    }
  }

  return csi;
}
/**************************************-**************************************/
/*                                                                           */
/*                               POPC_CSAV_SET_F                             */
/*                                                                           */
/*  Store the data trace in the appropriate csav location for the cell.      */
/*                                                                           */
/*****************************************************************************/
void popc_csav_set_f(c,datid,d,n,samp)
     struct pop_cell *c;
     char *datid;
     float *d;
     int n;
     float samp;
{
  int csi;
  struct pop_layer *pl;

  pl = c->pl;

  if (pl == NULL)
    exit_error("POPC_CSAV_SET_F","cell layer is NULL");
  else if (pl->csav_n == 0)
    exit_error("POPC_CSAV_SET_F","layer has no csav");
    
  csi = search_2d_carray(pl->csav_datid,datid,pl->csav_n);
  if (csi < 0){
    exit_error("POPC_CSAV_SET_F","datid not found in layer csav list");
  }else if ((c->csav_cnt == NULL) || (c->csav_f == NULL))
    exit_error("POPC_CSAV_SET_F","no storage for cell csav");

  if (c->csav_f[csi] != NULL){
    myfree(c->csav_f[csi]);    // Assume this was filled on a previous trial
  }

  c->csav_f[csi] = d;
  c->csav_cnt[csi] = n;

  if (mod_util_sampling_equal(samp,pl->csav_samp[csi]) == 0){
    printf("  data samp %f   csav_samp %f\n",samp,pl->csav_samp[csi]);
    exit_error("POPC_CSAV_SET_F","Sampling mismatch for csav");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                      POP_CELL_SYN_GET_PTR_DATA_FLOAT                      */
/*                                                                           */
/*****************************************************************************/
void pop_cell_syn_get_ptr_data_float(syn,name,rx,ry,rn,rxflag,rtscale)
  struct pop_syn *syn;
  char *name;
  float **rx,**ry;
  int *rn;
  int *rxflag;           /*  How to handle x-coords:
			     0 - x-data is included
			     1 - lgn time base
			     2 - spike data
			 */
  float *rtscale;  /* Multiply by this value to get sec */
{
  int n,xflag;
  float *x,*y,tscale;

  n = 0;
  x = (float *)NULL;
  y = (float *)NULL;
  xflag = 0;
  tscale = 1.0;

  if ((strcmp(name,"si_mask")==0) ||
      (strcmp(name,"si_mask1")==0) ||
      (strcmp(name,"si_mask2")==0)){
    if (syn->si != NULL){
      y = pop_cell_si_get_data_ptr(NULL,syn,name,&n);
      xflag = 1;       // WYETH - what to do w/ time base?
      tscale = 0.001;  // WYETH should match to 'samp' in si_sav
    }
  }else{
    printf("  *** name = %s\n",name);
    exit_error("POP_CELL_SYN_GET_PTR_DATA_FLOAT","Unknown name");
  }
  
  *rx = x;
  *ry = y;
  *rn = n;
  *rxflag = xflag;
  *rtscale = tscale;
}
/**************************************-**************************************/
/*                                                                           */
/*                         POP_CELL_GET_PTR_DATA_FLOAT                       */
/*                                                                           */
/*****************************************************************************/
void pop_cell_get_ptr_data_float(c,name,rx,ry,rn,rxflag,rtscale)
  struct pop_cell *c;
  char *name;
  float **rx,**ry;
  int *rn;
  int *rxflag;           /*  How to handle x-coords:
			     0 - x-data is included
			     1 - lgn time base
			     2 - spike data
			 */
  float *rtscale;  // Multiply by this value to get sec
{
  int n,xflag,csi;
  float *x,*y,tscale;

  n = 0;
  x = (float *)NULL;
  y = (float *)NULL;
  xflag = 0;
  tscale = 1.0;

  if (strcmp(name,"spikes")==0){
    if (c->savs){
      xflag = 2;
      y = c->s;
      n = c->ns;
      tscale = 0.001;
    }
  }else if (strcmp(name,"lgn_gx")==0){
    if (c->savg){
      exit_error("POP_CELL_GET_PTR_DATA_FLOAT","Not implemented.");
    }
    /*    
	  if (c->pl->nmda_p != NULL){
	  gtot = add_farrays(c->gtx0,c->gtn0,c->n0);
	  mod_util_resp_store_f(r,i,gtot,c->n0,0,mylogf);
	  }else{
	  mod_util_resp_store_f(r,i,c->gtx0,c->n0,1,mylogf);
	  }*/
    
  }else if (strcmp(name,"lgn_ga")==0){
    if (c->savg){
      xflag = 1;
      y = c->gtx0;
      n = c->n0;
    }
  }else if (strcmp(name,"lgn_gn")==0){
    if (c->savg){
      if (c->pl->nmda_p != NULL){
	xflag = 1;
	y = c->gtn0;
	n = c->n0;
      }
    }
  }else if (strcmp(name,"lgn_gain")==0){
    if ((c->savg) && (c->ggain != NULL)){
      xflag = 1;
      y = c->ggain;
      n = c->n0;
    }
  }else if (strcmp(name,"ex_gx")==0){
    if (c->savg){
      x = c->gtex1t;
      y = c->gtex1;
      n = c->gtex1n;
      tscale = 0.001;
    }
  }else if (strcmp(name,"in_gi")==0){
    if (c->savg){
      x = c->gtex1t;
      y = c->gtin1;
      n = c->gtex1n;
      tscale = 0.001;
    }
  }else if (strcmp(name,"gad")==0){
    if (c->savg){
      //x = c->gtex1t;
      xflag = 1;
      y = c->gta0;
      n = c->x * c->samp0;  // Number of values filled in so far
      //n = c->n0;
    }
  }else if ((strcmp(name,"vm")==0)||(strcmp(name,"rate")==0)){
    if (c->savv){
      x = c->vmt;
      y = c->vm;
      n = c->vn;
    }
  }else if (strcmp(name,"vd")==0){
    if (c->savd){
      x = c->vmt;
      y = c->vmd;
      n = c->vn;
    }
  }else{
    //
    //  Check for csav
    //
    csi = popc_csav_get_csi(c,name);
    if (csi >= 0){
      xflag = 1;  // Use LGN time base
      y = c->csav_f[csi];
      n = c->csav_cnt[csi];
    }else{
      printf("  *** name = %s\n",name);
      exit_error("POP_CELL_GET_PTR_DATA_FLOAT","Unknown name");
    }
  }
  
  *rx = x;
  *ry = y;
  *rn = n;
  *rxflag = xflag;
  *rtscale = tscale;
}
/**************************************-**************************************/
/*                                                                           */
/*                           POP_CELL_GET_DATA_FLOAT                         */
/*                                                                           */
/*****************************************************************************/
void pop_cell_get_data_float(c,name,rx,ry,rn,rxflag,rtscale,rptrflag)
  struct pop_cell *c;
  char *name;
  float **rx,**ry;
  int *rn;
  int *rxflag;     /* See comments in '...GET_PTR_DATA_FLOAT' */
  float *rtscale;  /* Multiply by this value to get sec */
  int *rptrflag;   /* 1 - ptr to data, DO NOT FREE
		      0 - y-data needs to be freed */
{
  int yn;
  float *y;

  // Handle derived quantities first
  if (strcmp(name,"lgn_ia")==0){
    // Most return values can come from 'lgn_ga'
    pop_cell_get_ptr_data_float(c,"lgn_ga",rx,&y,rn,rxflag,rtscale);

    // WYETH - i think c->ggain gets applied before here? so it is OK,
    // to derive current based on gtx0?

    y = mod_util_resp_derive_current(c->gtx0,c->samp0,(float *)NULL,c->n0,
				     *rtscale,c->vm,c->vmt,c->vn,
				     c->cifcp->v_th_x,
				     c->cifcp->v_ex,
				     (char *)NULL,rn); // Get acutal 'rn'
    *ry = y;
    *rptrflag = 0; // Not a pointer, this is malloc'd storage
  }else if (strcmp(name,"in_ii")==0){
    // Most return values can come from 'lgn_ga'
    pop_cell_get_ptr_data_float(c,"in_gi",rx,&y,rn,rxflag,rtscale);

    y = mod_util_resp_derive_current(c->gtin1,0.0,c->gtex1t,c->gtex1n,
				     *rtscale,c->vm,c->vmt,c->vn,
				     c->cifcp->v_th_x,
				     c->cifcp->v_in,
				     (char *)NULL,rn); /* Get acutal 'rn' */
    *ry = y;
    *rptrflag = 0; // Not a pointer, this is malloc'd storage
  }else if (strcmp(name,"ex_ix")==0){
    // Most return values can come from 'lgn_ga'
    pop_cell_get_ptr_data_float(c,"ex_gx",rx,&y,rn,rxflag,rtscale);

    y = mod_util_resp_derive_current(c->gtex1,0.0,c->gtex1t,c->gtex1n,
				     *rtscale,c->vm,c->vmt,c->vn,
				     c->cifcp->v_th_x,
				     c->cifcp->v_ex,
				     (char *)NULL,rn); /* Get acutal 'rn' */
    *ry = y;
    *rptrflag = 0; // Not a pointer, this is malloc'd storage
  }else{  // This is not a derived value, just get the pointer
    //printf("POP_CELL_GET_DATA_FLOAT,  name %s, NOT DERIVED VALUE\n",name);
    pop_cell_get_ptr_data_float(c,name,rx,ry,rn,rxflag,rtscale);
    //printf("    *rn = %d\n",*rn);
    *rptrflag = 1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         POP_CELL_SYN_GET_DATA_FLOAT                       */
/*                                                                           */
/*  Like POP_CELL_GET_DATA_FLOAT, but for syn instead of cell.               */
/*                                                                           */
/*****************************************************************************/
void pop_cell_syn_get_data_float(syn,name,rx,ry,rn,rxflag,rtscale,rptrflag)
  struct pop_syn *syn;
  char *name;
  float **rx,**ry;
  int *rn;
  int *rxflag;     // See comments in '...GET_PTR_DATA_FLOAT'
  float *rtscale;  // Multiply by this value to get sec
  int *rptrflag;   // 1 - ptr to data, DO NOT FREE
		   // 0 - y-data needs to be freed
{
  // This is not a derived value, just get the pointer
  pop_cell_syn_get_ptr_data_float(syn,name,rx,ry,rn,rxflag,rtscale);

  *rptrflag = 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                        POP_CELL_WRITE_SAV_LAYER_LIST                      */
/*                                                                           */
/*  Write output for:                                                        */
/*    savg                                                                   */
/*    savv                                                                   */
/*                                                                           */
/*****************************************************************************/
void pop_cell_write_sav_layer_list(laylist,nlay,outfile,sname,rpti)
     struct pop_layer **laylist;
     int nlay;
     char outfile[];   // Outfile
     char sname[];     // Stimulus name
     int rpti;         // Repeat number
{
  int i,j,k,l;
  float *x;
  struct pop_cell *c;
  char tname[SLEN];

  // WYETH - anything that has 'savv' or 'savg' will be dumped here,
  //   even if it was intended as 'ndata' response

  for(l=0;l<nlay;l++){ // For each layer
    for(i=0;i<laylist[l]->xn;i++){
      for(j=0;j<laylist[l]->yn;j++){

	if ((laylist[l]->c[i][j] != NULL) &&  // LGN layers, maybe no storage?
	    (strcmp(laylist[l]->laytype,"lgn")!=0)){
	  // Don't do LGN layers, because there are different storage
	  // conventions for sav's.  Condition added Aug 2009
	  
	  for(k=0;k<laylist[l]->zn;k++){ // For each cell
	    c = &(laylist[l]->c[i][j][k]);
	    if (c != NULL){
	      if (c->savv){
		x = copy_farray(c->vmt,c->vn);
		multiply_farray(x,c->vn,1000.0);
		
		sprintf(tname,"%s_Vm_rpt_%d_%s",c->name,rpti,sname);
		append_farray_xy_plot(outfile,x,c->vm,c->vn,tname);
		
		if (c->savd){
		  sprintf(tname,"%s_VmD_rpt_%d_%s",c->name,rpti,sname);
		  append_farray_xy_plot(outfile,x,c->vmd,c->vn,tname);
		}
	      }
	      if (c->savg){
		sprintf(tname,"%s_Dyn_gx_rpt_%d_%s",c->name,rpti,sname);
		append_farray_xy_plot(outfile,c->gtex1t,c->gtex1,c->gtex1n,
				      tname);
		
		sprintf(tname,"%s_Dyn_gi_rpt_%d_%s",c->name,rpti,sname);
		append_farray_xy_plot(outfile,c->gtex1t,c->gtin1,c->gtex1n,
				      tname);
		
		sprintf(tname,"%s_LGN_gx_rpt_%d_%s",c->name,rpti,sname);
		append_farray_plot(outfile,tname,c->gtx0,c->n0,1);

		if (c->ggain != NULL){
		  sprintf(tname,"%s_LGN_gGain_rpt_%d_%s",c->name,rpti,sname);
		  append_farray_plot(outfile,tname,c->ggain,c->n0,1);
		}
		
		if (laylist[l]->nmda_p != NULL){
		  sprintf(tname,"%s_LGN_gNMDA_rpt_%d_%s",c->name,rpti,sname);
		  append_farray_plot(outfile,tname,c->gtn0,c->n0,1);
		}
		
		sprintf(tname,"%s_LGN_gi_rpt_%d_%s",c->name,rpti,sname);
		append_farray_plot(outfile,tname,c->gti0,c->n0,1);

		if (c->gta0 != NULL){
		  sprintf(tname,"%s_Adapt_g_rpt_%d_%s",c->name,rpti,sname);
		  append_farray_plot(outfile,tname,c->gta0,c->n0,1);
		}
		
		if (c->nad > 0)
		  append_farray_plot(outfile,"gta_pulse",c->gta_pulse,c->nad,
				     1);
	      }
	    }
	  }


	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             POPC_PRECOMP_RESET                            */
/*                                                                           */
/*  Reset the pre-computed data (creating storage if needed).                */
/*                                                                           */
/*****************************************************************************/
void popc_precomp_reset(c,lay,tn)
     struct pop_cell *c;
     struct pop_layer *lay;
     int tn;
{

  //  (1)  Clear 'gta0'
  //  (2)  Clear 'gtx0', 'gti0', and 'gtn0'
  //  (3)  Clear 'ggain'

  if (c->gta0 != NULL){
    zero_farray(c->gta0,tn);
  }

  if (c->gtx0 == NULL){
    //
    //  Create storage
    //
    c->gtx0 = get_zero_farray(tn);  // Initialize x, i, and a,
    c->gti0 = get_zero_farray(tn);  //   assuming that all 3 linked
    if (lay->nmda_p != NULL)
      c->gtn0 = get_zero_farray(tn);

    if (lay->lgn_gain == 1){
      c->ggain = get_zero_farray(tn);
    }
    c->n0 = tn;
    // c->samp0 should already be set

  }else{
    //
    //  Clear (zero) existing storage.
    //
    if (tn != c->n0)
      exit_error("POPC_PRECOMP_RESET","tn changed");

    zero_farray(c->gtx0,tn);
    zero_farray(c->gti0,tn);
    if (lay->nmda_p != NULL)
      zero_farray(c->gtn0,tn);

    if (lay->lgn_gain == 1){
      if (c->ggain == NULL){
	printf("===================== c->ggain is NULL\n");
      }else{
	zero_farray(c->ggain,tn);
      }
    }
  }
  
}
/**************************************-**************************************/
/*                                                                           */
/*                      POPC_PRECOMP_GET_EXPANDED_SPIKES                     */
/*                                                                           */
/*****************************************************************************/
float *popc_precomp_get_expanded_spikes(c,tn,sdf,sdtau)
     struct pop_cell *c;
     int tn;
     float sdf,sdtau;
{
  int xi,yi,zi,n;
  float *s,*x,w;
  struct pop_syn *t;
  struct pop_cell *cpre;

  //printf("POPC_PRECOMP_GET_EXPANDED_SPIKES  cell %s  %d\n",c->name,c->nin);

  x = get_zero_farray(tn);

  t = c->in;  // pointer to first pre-syn input synapse
  while (t != NULL){
    if (t->stype == 11){  // This syn comes from a pre-compute layer

      cpre = t->pre;  // Pre-syn cell

      xi = cpre->layx;  // Coords of pre-syn cell
      yi = cpre->layy;
      zi = cpre->layz;
      n = cpre->pl->cnt[zi][xi][yi];   // Number of spikes
      s = cpre->pl->s[zi][xi][yi];     // Spike train

      //printf("  popc---       %d spikes  (%d,%d,%d)\n",n,xi,yi,zi);

      spikeu_expand_spikes_wsd_accum(s,n,x,tn,t->w,sdf,sdtau);

    }
    //printf("    next syn  (cell %s)\n",t->pre->name);
    t = t->post_next;
  }

  return x;
}
/**************************************-**************************************/
/*                                                                           */
/*                             POP_CELL_INIT_CELL                            */
/*                                                                           */
/*  Initialize cell fields to safe values.                                   */
/*  No storage is created here!                                              */
/*                                                                           */
/*****************************************************************************/
void pop_cell_init_cell(c)
     struct pop_cell *c;
{
  c->name = (char *)NULL;
  c->subclass = (char *)NULL;
  c->pl = (struct pop_layer *)NULL;
  c->cx = c->cy = 0.0;
  c->rfx = 0.0;
  c->rfy = 0.0;

  c->cifcp = NULL;
  c->cpoissp = NULL;  // 2012 Feb 6

  c->attrib_f = NULL;

  c->nin = c->nout = 0;
  c->in = c->out = NULL;

  c->syn_proc_code = 0;

  // Precomputed conductances
  c->n0 = -1;
  c->samp0 = 0.0;
  c->gtx0 = c->gtn0 = c->gti0 = c->gta0 = (float *)NULL;
  c->ggain = NULL;

  c->bg_x_rate = c->bg_i_rate = 0.0;
  c->bg_x_s = c->bg_i_s = NULL;
  c->bg_x_corr_c = c->bg_i_corr_c = NULL;
  c->bg_x_n = c->bg_i_n = -1;               // No spikes generated yet
  c->bg_x_k = c->bg_i_k = -1;  // Needed for test in pop_util: pop_input_..

  c->savv = c->savg = c->savd = c->savs = c->sava = 0;

  c->lgn_mflag = 0; // WYETH - TEMPORARY

  c->lgnin = NULL;
  c->lgn_n = 0;   // Start with no LGN inputs

  c->s = NULL;
  c->ns = 0;

  c->csav_cnt = NULL;
  c->csav_f = NULL;
}
/**************************************-**************************************/
/*                                                                           */
/*                             POP_CELL_PRINT_CELL                           */
/*                                                                           */
/*****************************************************************************/
void pop_cell_print_cell(c)
     struct pop_cell *c;
{
  int i;

  printf("  POP_CELL_PRINT_CELL\n");
  if (c->name != (char *)NULL)
    printf("    Name:  %s\n",c->name);
  if (c->subclass != (char *)NULL)
    printf("    Subclass:  %s\n",c->subclass);
  else
    printf("    Subclass:  NULL\n");
  if (c->pl != NULL)
    printf("    Layer:  %s\n",c->pl->name);

  printf("    Location (cx,cy):  %f, %f\n",c->cx,c->cy);

  printf("    Attribs n = %d\n",c->pl->attrib_fn);
  for(i=0;i<c->pl->attrib_fn;i++)
    printf("      %12s  %f\n",c->pl->attrib_fname[i],c->attrib_f[i]);

  printf("    Synaptic Connections:\n");
  printf("      nin, nout  %d %d\n",c->nin,c->nout);
  if (c->out == NULL)
    printf("    c->out is NULL\n");
  if (c->in == NULL)
    printf("    c->in is NULL\n");

  printf("    n0    = %d\n",c->n0);
  printf("    samp0 = %f\n",c->samp0);

  printf("    Adaptation pulse array has %d elements.\n",c->nad);
  
  printf("    x = %f\n",c->x);
  printf("    y = %f\n",c->y);
  printf("    ystart = %f\n",c->ystart);
  printf("    yscal = %f\n",c->yscal);

  printf("    h    = %f\n",c->h);
  printf("    h1   = %f\n",c->h1);
  printf("    hmin = %f\n",c->hmin);
  printf("    hmax = %f\n",c->hmax);
  
  printf("    nbad = %d\n",c->nbad);
  printf("    nok  = %d\n",c->nok);

  printf("    eps  = %f\n",c->eps);

  printf("    savv  = %d\n",c->savv);
  if (c->savv){
    printf("      dxsav = %f\n",c->dxsav);
    printf("      xsav = %f\n",c->xsav);
    printf("      vnmax = %d\n",c->vnmax);
    printf("      vn    = %d\n",c->vn);
  }
  printf("    savd  = %d\n",c->savd);

  printf("    savs  = %d\n",c->savs);
  if (c->savs){
    printf("      maxscnt = %d\n",c->maxscnt);
    printf("      ns      = %d\n",c->ns);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          POP_CELL_ATTRIB_INIT_F                           */
/*                                                                           */
/*****************************************************************************/
void pop_cell_attrib_init_f(pl,n,aname,cellflag)
     struct pop_layer *pl;   // Layer post-syn for input
     int n;                  // Number of attribs
     char **aname;           // [n] Attribute names
     int cellflag;           // 1-initialize cell storage also, 0-not
{
  int i,j,k;
  struct pop_cell *c;

  pl->attrib_fn    = n;
  pl->attrib_fname = copy_2d_carray(aname,n);

  if (cellflag == 1){
    //
    // Allocate the storage for the attributes in the cells
    //
    for(i=0;i<pl->xn;i++)
      for(j=0;j<pl->yn;j++)
	for(k=0;k<pl->zn;k++){
	  c = &(pl->c[i][j][k]);
	  c->attrib_f = (float *)myalloc(n*sizeof(float));
	}
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           POP_CELL_ATTRIB_INDEX_F                         */
/*                                                                           */
/*  Return the attrib index, or -1 if not found.                             */
/*                                                                           */
/*****************************************************************************/
int pop_cell_attrib_index_f(pl,tname)
     struct pop_layer *pl;
     char *tname;
{
  int i,k;
  int n;
  char **alist;

  n     = pl->attrib_fn;
  alist = pl->attrib_fname;

  k = -1;
  i = 0;
  while(i<n){
    if (strcmp(tname,alist[i])==0){
      k = i;
      i = n; // end loop
    }else
      i += 1;
  }

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                             POP_CELL_ATTRIB_SET                           */
/*                                                                           */
/*****************************************************************************/
void pop_cell_attrib_set(c,aname,fval,caller)
     struct pop_cell *c;  // cell
     char *aname;         // attrib name
     float fval;          // attrib value
     char *caller;        // name of caller, exit on failure if this not NULL
{
  int ai;
  struct pop_layer *pl;

  pl = c->pl;
  ai = pop_cell_attrib_index_f(pl,aname);
  if (ai != -1)
    c->attrib_f[ai] = fval;
  else if (caller != NULL){
    printf("  *** Called by %s\n",caller);
    exit_error("POP_CELL_ATTRIB_SET","Attrib name not found");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_CELL_ATTRIB_GET_F                          */
/*                                                                           */
/*****************************************************************************/
float pop_cell_attrib_get_f(c,tname)
     struct pop_cell *c;
     char *tname;
{
  int i,k;
  int n;
  char **alist;

  n     = c->pl->attrib_fn;
  alist = c->pl->attrib_fname;

  k = -1;
  i = 0;
  while(i<n){
    if (strcmp(tname,alist[i])==0){
      k = i;
      i = n; // end loop
    }else
      i += 1;
  }

  if (k == -1){
    printf("  *** attrib_fn  %d\n",n);
    printf("  *** attrib name  %s\n",tname);
    printf("  *** cell name  %s\n",c->name);
    exit_error("POP_CELL_ATTRIB_GET_F","Attrib name not found");
  }

  return c->attrib_f[k];
}
/**************************************-**************************************/
/*                                                                           */
/*                           POP_CELL_ATTRIB_PRINT_F                         */
/*                                                                           */
/*****************************************************************************/
void pop_cell_attrib_print_f(c)
     struct pop_cell *c;
{
  int i;
  int n;
  char **alist;

  n     = c->pl->attrib_fn;
  alist = c->pl->attrib_fname;

  for(i=0;i<n;i++)
    printf("  %16s %f\n",alist[i],c->attrib_f[i]);

}
/**************************************-**************************************/
/*                                                                           */
/*                         POP_CELL_ATTRIB_DEFAULT_SET                       */
/*                                                                           */
/*****************************************************************************/
void pop_cell_attrib_default_set(c)
     struct pop_cell *c;   // cell
{
  int ai;
  struct pop_layer *pl;

  //
  //  Default values for otherwise uninitialized attributes.
  //
  //  Note, these will attempt to set attrib values, but will not fail
  //  if the attrib names are not found.
  //

  pl = c->pl;

  //
  //  Alphabetical Order.
  //
  pop_cell_attrib_set(c,"conetype",-1.0,NULL);  // Unassigned
  pop_cell_attrib_set(c,"ocdom",-1.0,NULL);     // Assume left eye
  pop_cell_attrib_set(c,"ori",0.0,NULL);
  pop_cell_attrib_set(c,"ori_cv",0.0,NULL);
  pop_cell_attrib_set(c,"phase",0.0,NULL);
  pop_cell_attrib_set(c,"sz1",0.0,NULL);
}
/**************************************-**************************************/
/*                                                                           */
/*                             POPC_ATTRIB_VARY_XY                           */
/*                                                                           */
/*  Return 1 if the attrib at index 'ai' varies in the (x,y) plane at 'zi',  */
/*  return 0 otherwise.                                                      */
/*                                                                           */
/*****************************************************************************/
int popc_attrib_vary_xy(pl,ai,zi)
     struct pop_layer *pl;    // Layer
     int ai;                  // Attrib index in layer
     int zi;                  // z-index to check
{
  int i,j;
  int xn,yn;
  float aval;

  if ((zi < 0) || (zi >= pl->zn))
    exit_error("POPC_ATTRIB_VARY_XY","Bad 'zi' value");
  if ((ai < 0) || (ai >= pl->attrib_fn))
    exit_error("POPC_ATTRIB_VARY_XY","Bad 'ai' value");

  xn = pl->xn;
  yn = pl->yn;

  aval = pl->c[0][0][zi].attrib_f[ai];  // Attrib value of first unit
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      if (aval != pl->c[i][j][zi].attrib_f[ai])
	return 1;       // Varies across (x,y), no z variation possible
    }
  }
  return 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                         POPC_ATTRIB_VARIATION_CODE                        */
/*                                                                           */
/*  Examine the attrib at index 'ai' in this layer; report how it varies:    */
/*    0  - does not vary across cells                                        */
/*    1  - varies across cells, but there is no variation in 'z'             */
/*    2  - varies with 'z' but there is no variation across (x,y)            */
/*    3  - varies across (x,y) and z.                                        */
/*                                                                           */
/*****************************************************************************/
int popc_attrib_variation_code(pl,ai)
     struct pop_layer *pl;    // Layer
     int ai;                  // Attrib index in layer
{
  int i,j,k;
  int xn,yn,zn,varcode,done;
  float aval,tval;

  xn = pl->xn;
  yn = pl->yn;
  zn = pl->zn;

  if (zn == 1){  // No z-variation
    varcode = popc_attrib_vary_xy(pl,ai,0);  // Return either 0 or 1
  }else{
    varcode = popc_attrib_vary_xy(pl,ai,0);
    if (varcode == 0){  // No variation in layer zi=0
      for(i=1;i<zn;i++)
	varcode += popc_attrib_vary_xy(pl,ai,i);
      if (varcode == 0){  // No variation in any z-sheet

	aval = pl->c[0][0][0].attrib_f[ai];  // Attrib value of first unit
	for(i=1;i<zn;i++){
	  if (aval != pl->c[0][0][i].attrib_f[ai]){
	    varcode = 2;  // Varies in z, but not in (x,y)
	  }
	}
	//
	//  Varcode is now either 0 or 2.
	//
      }else{
	varcode = 3;  // flag that some (x,y) variation found
      }
    }else{
      varcode = 3;  // flag that some (x,y) variation found
    }

    if (varcode == 3){
      //
      //  There must have been variation found across (x,y), now we check z.
      //  Thus, 1 or 3 are the only possibilities now.
      //

      varcode = 1;  // Assume there will be no z-variation
      i = 0;
      j = 0;
      done = 0;
      while(done == 0){
	aval = pl->c[i][j][0].attrib_f[ai];  // Val for first unit at (i,j)
	k = 1;
	while((k<zn) && (done == 0)){
	  if (aval != pl->c[i][j][k].attrib_f[ai]){
	    varcode = 3;  // There is z-variation also
	    done = 1;
	  }else{
	    k += 1;
	  }
	}
	j += 1;
	if (j == yn){
	  j = 0;
	  i += 1;
	  if (i == xn)
	    done = 1;
	}
      }
    }
  }

  return varcode;
}
/**************************************-**************************************/
/*                                                                           */
/*                         POPC_ATTRIB_LAYER_VARIATION                       */
/*                                                                           */
/*  For each attrib in this layer, return an int array that encodes how      */
/*  the attrib varies:                                                       */
/*    0  - does not vary across cells                                        */
/*    1  - varies across cells, but there is no variation in 'z'             */
/*    2  - varies with 'z' but there is no variation across (x,y)            */
/*    3  - varies across (x,y) and z.                                        */
/*                                                                           */
/*  Returned array has length 'pl->attrib_fn'                                */
/*                                                                           */
/*****************************************************************************/
int *popc_attrib_layer_variation(pl)
     struct pop_layer *pl;
{
  int i;
  int n,*v;

  n = pl->attrib_fn;  // Number of floating point attribs
  v = (int *)myalloc(n*sizeof(int));

  for(i=0;i<n;i++){
    v[i] = popc_attrib_variation_code(pl,i);
  }

  return v;
}
/**************************************-**************************************/
/*                                                                           */
/*                            POPC_ATTRIB_DUMP_LAYER                         */
/*                                                                           */
/*  Write the attrib value for each unit in the layer.                       */
/*                                                                           */
/*****************************************************************************/
void popc_attrib_dump_layer(pl,outfile,aname)
     struct pop_layer *pl;
     char *outfile;
     char *aname;
{
  int i,j,k;
  int xn,yn,zn;
  float fval;
  char tstr[SLEN];
  struct pop_cell *c;

  xn = pl->xn;
  yn = pl->yn;
  zn = pl->zn;

  remove_file(outfile);

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	c = &(pl->c[i][j][k]);
	fval = pop_cell_attrib_get_f(c,aname);
	sprintf(tstr,"%s %d %d %d %s %f\n",pl->name,i,j,k,aname,fval);
	append_string_to_file(outfile,tstr);
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                        POPC_LAY_COUNT_INPUT_FROM_LAY                      */
/*                                                                           */
/*  Return the number of inputs in the input list of 'lpost' that come from  */
/*  'lpre'.                                                                  */
/*                                                                           */
/*****************************************************************************/
int popc_lay_count_input_from_lay(lpre,lpost)
     struct pop_layer *lpre,*lpost;
{
  int i;
  int n,cnt;
  char *oname;
  struct onode *o;

  n = lpost->ninlist;

  cnt = 0;
  for(i=0;i<n;i++){  // For each input
    o = lpost->inlist[i];   // Get the onode for this input
    oname = onode_getpar_chr_ptr(o,"pop_origin");
    if (oname != NULL){
      if (strcmp(oname,lpre->name)==0){
	cnt += 1;
      }
    }
  }

  return cnt;
}
/**************************************-**************************************/
/*                                                                           */
/*                              POP_CELL_SYN_FREE                            */
/*                                                                           */
/*  This should be the ONLY routine to free a synapse.                       */
/*  Other search terms:  synapse_free free_synpase                           */
/*                                                                           */
/*****************************************************************************/
void pop_cell_syn_free(s)
     struct pop_syn *s;
{

  // If there is a synaptic interaction
  if (s->si != NULL){
    if (s->si->ci == 0){  // SI index 0 controls 'siu' storage
      pop_cell_si_free(s->si);
    }
  }

  myfree(s);
}
/**************************************-**************************************/
/*                                                                           */
/*                           POP_CELL_PRINT_SYNINT                           */
/*                                                                           */
/*****************************************************************************/
void pop_cell_print_synint(s)
     struct pop_syn *s;
{
  struct pop_si *si;
  
  si = s->si;
  
  if (si->siui == 1){
    printf("       SI type %d, cell index %d\n",si->siui,si->ci);
  }else if (si->siui == 2){
    printf("       SI type %d, cell index %d\n",si->siui,si->ci);
  }else{
    exit_error("POP_CELL_PRINT_SYNINT","Unknown SI Union Index");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                      POP_CELL_PRINT_PRE_SYNAPTIC_LIST                     */
/*                                                                           */
/*****************************************************************************/
void pop_cell_print_pre_synaptic_list(c)
     struct pop_cell *c;
{
  int i;
  struct pop_syn *t;

  printf("  Pre-Synaptic Inputs\n");
  t = c->in;
  i = 1;
  while(t != NULL){
    printf("  %6d %s %.3f %.3f  %s   type %d  weight %f  tlast %f  alast %f tdelay %f\n",
	   i,t->pre->name,
	   t->pre->cx,t->pre->cy,
	   t->post->name,t->stype,t->w,t->tlast,t->alast,
	   t->tdelay);
    if (t->si != NULL)
      pop_cell_print_synint(t);

    t = t->post_next;
    i += 1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                      POP_CELL_PRINT_POST_SYNAPTIC_LIST                    */
/*                                                                           */
/*****************************************************************************/
void pop_cell_print_post_synaptic_list(c)
     struct pop_cell *c;
{
  int i;
  struct pop_syn *t;

  printf("  Post-Synaptic Targets\n");
  t = c->out;
  i = 1;
  while(t != NULL){
    printf("  %6d %s %s   type %d  weight %f  tlast %f  alast %f tdelay %f\n",
	   i,t->pre->name,t->post->name,t->stype,t->w,t->tlast,t->alast,
	   t->tdelay);
    if (t->si != NULL)
      pop_cell_print_synint(t);

    t = t->pre_next; // Traverse synapses associate with same pre-syn cell
    i += 1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         POP_CELL_GET_NEXT_SYN_CELL                        */
/*                                                                           */
/*  Return a pointer to the next synapse in the list to/from the named cell. */
/*                                                                           */
/*****************************************************************************/
struct pop_syn *pop_cell_get_next_syn_cell(s,pflag,layname,xi,yi,zi)
     struct pop_syn *s;
     int pflag;  // 0-presyn, 1-postsyn
     char layname[];
     int xi,yi,zi;
{
  struct pop_syn *syn,*t;

  syn = NULL;

  t = s;
  if (pflag == 0){
    while ((t != NULL) && (syn == NULL)){
      if ((strcmp(t->post->pl->name,layname)==0) &&
	  (t->post->layx == xi) &&
	  (t->post->layy == yi) &&
	  (t->post->layz == zi))
	syn = t;
      else{
	/*printf("  layer name:  %s\n",t->post->pl->name);*/
	t = t->pre_next;
      }
    }
  }else{
    while ((t != NULL) && (syn == NULL)){
      if ((strcmp(t->pre->pl->name,layname)==0) && 
	  (t->pre->layx == xi) &&
	  (t->pre->layy == yi) &&
	  (t->pre->layz == zi))
	syn = t;
      else{
	/*printf("  layer name:  %s\n",t->pre->pl->name);*/
	t = t->post_next;
      }
    }
  }

  return syn;
}
/**************************************-**************************************/
/*                                                                           */
/*                      POP_CELL_GET_NEXT_SYN_LAYER_NAME                     */
/*                                                                           */
/*  Return a pointer to the next synapse in the list to a cell in the named  */
/*  layer.                                                                   */
/*                                                                           */
/*****************************************************************************/
struct pop_syn *pop_cell_get_next_syn_layer_name(s,pflag,layname)
     struct pop_syn *s;
     int pflag;  // 0-presyn, 1-postsyn
     char layname[];
{
  struct pop_syn *syn,*t;

  syn = NULL;

  t = s;
  if (pflag == 0){ // Find synapses going to this layer

    while ((t != NULL) && (syn == NULL)){

      if (t->post->pl == NULL){  // *** WYETH ADDED Sep 5, 2015
	exit_error("POP_CELL_GET_NEXT_SYN_LAYER_NAME","Postsyn layer is NULL");
      }
      
      if (strcmp(t->post->pl->name,layname)==0)
	syn = t;
      else{
	//printf("  layer name:  %s\n",t->post->pl->name);
	t = t->pre_next;
      }
    }
  }else{ // Find synapses coming from this layer

    while ((t != NULL) && (syn == NULL)){

      if (t->pre->pl == NULL){  // *** WYETH ADDED Sep 5, 2015
	exit_error("POP_CELL_GET_NEXT_SYN_LAYER_NAME","Presyn layer is NULL");
      }

      if (strcmp(t->pre->pl->name,layname)==0){
	syn = t;
      }else{
	//printf("  layer name:  %s\n",t->pre->pl->name);
	t = t->post_next;
      }
    }
  }

  return syn;
}
/**************************************-**************************************/
/*                                                                           */
/*                         POP_CELL_GET_NEXT_SYN_INDEX                       */
/*                                                                           */
/*  Return a pointer to the next synapse in the list to a cell in the named  */
/*  layer and with the input index value 'inindex'.                          */
/*                                                                           */
/*****************************************************************************/
struct pop_syn *pop_cell_get_next_syn_index(s,pflag,inindex,layname)
     struct pop_syn *s;
     int pflag;          // 0-presyn, 1-postsyn
     short inindex;
     char layname[];     // NULL if pflag = 0, inindex is unique among presyn
{
  struct pop_syn *syn,*t;

  syn = NULL;

  t = s;
  if (pflag == 0){ // Find synapses going to this layer
    while ((t != NULL) && (syn == NULL)){
      if ((strcmp(t->post->pl->name,layname)==0)&&
	  (t->inindx == inindex))
	syn = t;
      else{
	t = t->pre_next;
      }
    }
  }else{ // Find synapses coming from this layer
    while ((t != NULL) && (syn == NULL)){
      if (t->inindx == inindex)
	syn = t;
      else{
	t = t->post_next;
      }
    }
  }

  return syn;
}
/**************************************-**************************************/
/*                                                                           */
/*                      POP_CELL_LAYLIST_GET_SYN_POINTER                     */
/*                                                                           */
/*  Return the pointer to the 'k'th synapse between the units.               */
/*                                                                           */
/*****************************************************************************/
struct pop_syn *pop_cell_laylist_get_syn_pointer(laylist,n,name,xi,yi,zi,
						 name1,xi1,yi1,zi1,k)
     struct pop_layer **laylist;
     int n;
     char *name;       /* Pre-syn unit */
     int xi,yi,zi;
     char *name1;      /* Post-syn unit */
     int xi1,yi1,zi1;
     int k;            /* get the kth synapse, 0,1,.. */
{
  int i;
  struct pop_cell *c,*c1;
  struct pop_syn *s,*t;

  /*
    printf("name xyz  %s  %d %d %d\n",name,xi,yi,zi);
    printf("name1 xyz1  %s  %d %d %d\n",name1,xi1,yi1,zi1);
  */

  /* Make sure both cells exist */
  c = pop_cell_laylist_get_cell_pointer(laylist,n,name,xi,yi,zi);
  if (c == NULL)
    exit_error("pop_cell_laylist_get_syn_pointer","Can't find pre-syn cell");

  c1 = pop_cell_laylist_get_cell_pointer(laylist,n,name1,xi1,yi1,zi1);
  if (c1 == NULL)
    exit_error("pop_cell_laylist_get_syn_pointer","Can't find post-syn cell");

  /*** Get next 'k'th synapse from c to c1 */
  s = pop_cell_get_next_syn_cell(c->out,0,name1,xi1,yi1,zi1);
  if (k > 0){
    i = 0;
    while((s != NULL) && (i < k)){
      i += 1;
      t = s->pre_next;
      s = pop_cell_get_next_syn_cell(t,0,name1,xi1,yi1,zi1);
    }
  }
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                    POP_CELL_COUNT_SYN_FROM_LAYER_TO_CELL                  */
/*                                                                           */
/*****************************************************************************/
int pop_cell_count_syn_from_layer_to_cell(layname,c)
     char layname[];
     struct pop_cell *c;
{
  int n;
  struct pop_syn *s,*t;

  // Scan pre-syn list for synapses from 'lay'
  s = pop_cell_get_next_syn_layer_name(c->in,1,layname);
  n = 0;
  while(s != NULL){
    n += 1;
    t = s->post_next; // The post-syn cell of an input syn is this cell
    s = pop_cell_get_next_syn_layer_name(t,1,layname);
  }
  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                    POP_CELL_MAX_WEIGHT_FROM_LAYER_TO_CELL                 */
/*                                                                           */
/*****************************************************************************/
float pop_cell_max_weight_from_layer_to_cell(layname,c)
     char layname[];
     struct pop_cell *c;
{
  int n;
  float wmax;
  struct pop_syn *s,*t;

  wmax = 0.0;

  // Scan pre-syn list for synapses from 'lay'
  s = pop_cell_get_next_syn_layer_name(c->in,1,layname);
  while(s != NULL){
    if (s->w > wmax)
      wmax = s->w;
    t = s->post_next; // The post-syn cell of an input syn is this cell
    s = pop_cell_get_next_syn_layer_name(t,1,layname);
  }
  return wmax;
}
/**************************************-**************************************/
/*                                                                           */
/*                    POP_CELL_TOT_WEIGHT_FROM_LAYER_TO_CELL                 */
/*                                                                           */
/*****************************************************************************/
float pop_cell_tot_weight_from_layer_to_cell(layname,c)
     char layname[];
     struct pop_cell *c;
{
  float wtot;
  struct pop_syn *s,*t;

  wtot = 0.0;

  // Scan pre-syn list for synapses from 'lay'
  s = pop_cell_get_next_syn_layer_name(c->in,1,layname);
  while(s != NULL){
    wtot += s->w;
    t = s->post_next; // The post-syn cell of an input syn is this cell
    s = pop_cell_get_next_syn_layer_name(t,1,layname);
  }
  return wtot;
}
/**************************************-**************************************/
/*                                                                           */
/*                     POP_CELL_TOT_WEIGHT_ININDEX_TO_CELL                   */
/*                                                                           */
/*****************************************************************************/
float pop_cell_tot_weight_inindex_to_cell(c,inindex)
     struct pop_cell *c;
     short inindex;
{
  float wtot;
  struct pop_syn *s,*t;

  wtot = 0.0;

  // Scan pre-syn list for synapses from 'lay'
  s = pop_cell_get_next_syn_index(c->in,1,inindex,NULL);
  //s = pop_cell_get_next_syn_layer_name(c->in,1,layname);
  while(s != NULL){
    wtot += s->w;
    t = s->post_next; // The post-syn cell of an input syn is this cell
    s = pop_cell_get_next_syn_index(t,1,inindex,NULL);
  //s = pop_cell_get_next_syn_layer_name(t,1,layname);
  }

  return wtot;
}
/**************************************-**************************************/
/*                                                                           */
/*                    POP_CELL_COUNT_SYN_FROM_CELL_TO_LAYER                  */
/*                                                                           */
/*****************************************************************************/
int pop_cell_count_syn_from_cell_to_layer(layname,c)
     char layname[];
     struct pop_cell *c;
{
  int n;
  struct pop_syn *s,*t;

  s = pop_cell_get_next_syn_layer_name(c->out,0,layname);
  n = 0;
  while(s != NULL){
    n += 1;
    t = s->pre_next;
    s = pop_cell_get_next_syn_layer_name(t,0,layname);
  }
  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                    POP_CELL_MAX_WEIGHT_FROM_CELL_TO_LAYER                 */
/*                                                                           */
/*****************************************************************************/
float pop_cell_max_weight_from_cell_to_layer(layname,c)
     char layname[];
     struct pop_cell *c;
{
  int n;
  float wmax;
  struct pop_syn *s,*t;

  wmax = 0.0;

  s = pop_cell_get_next_syn_layer_name(c->out,0,layname);
  while(s != NULL){
    if (s->w > wmax)
      wmax = s->w;
    t = s->pre_next;
    s = pop_cell_get_next_syn_layer_name(t,0,layname);
  }
  return wmax;
}
/**************************************-**************************************/
/*                                                                           */
/*                    POP_CELL_TOT_WEIGHT_FROM_CELL_TO_LAYER                 */
/*                                                                           */
/*****************************************************************************/
float pop_cell_tot_weight_from_cell_to_layer(layname,c)
     char layname[];
     struct pop_cell *c;
{
  float wtot;
  struct pop_syn *s,*t;

  wtot = 0.0;

  s = pop_cell_get_next_syn_layer_name(c->out,0,layname);
  while(s != NULL){
    wtot += s->w;
    t = s->pre_next;
    s = pop_cell_get_next_syn_layer_name(t,0,layname);
  }
  return wtot;
}
/**************************************-**************************************/
/*                                                                           */
/*                          POPC_WF_ACCUM_LGN_TO_UNIT                        */
/*                                                                           */
/*  Mark all units in the flag array 'cmap' that send connections to the     */
/*  post-syn unit.                                                           */
/*                                                                           */
/*  Return the number of *NEWLY* marked units.                               */
/*                                                                           */
/*****************************************************************************/
int popc_wf_accum_lgn_to_unit(lpre,lpost,ci,cj,ck,wp,flag,fmap,optype,fcrit)
     struct pop_layer *lpre;     // Pre-syn layer
     struct pop_layer *lpost;    // Post-syn layer
     int ci,cj,ck;               // Post-syn coords
     float wp;                   // prior weight for this unit 
     int   ***flag;              // [lpre->xn][->yn][->zn]
     float ***fmap;              // [lpre->xn][->yn][->zn]

     // *** THESE ARE IGNORED, GIVEN THAT LGN as w = 1
     char *optype;               // "w_min", "w_max"
     float fcrit;                // Critical value, e.g., min weight
{
  int i,j;
  int nsyn,nlgn,n,xi,yi,*x,*y;
  struct pop_cell *c;

  //printf("    wp %f  otype %s  fcrit %f\n",wp,optype,fcrit);

  nsyn = 0;  // Number of newly marked units

  if ((ci == -1)||(cj == -1)||(ck == -1))
    return nsyn;

  c = &(lpost->c[ci][cj][ck]);    // Pointer to the post-syn cell

  nlgn = c->lgn_n;
  for(i=0;i<nlgn;i++){
    if (c->lgnin[i]->lay == lpre){
      n = c->lgnin[i]->cn0;
      x = c->lgnin[i]->cx0;
      y = c->lgnin[i]->cy0;
      for(j=0;j<n;j++){
	xi = x[j];
	yi = y[j];
	fmap[xi][yi][0] += wp;    // Total weight for this pre-syn unit
	flag[xi][yi][0] += 1;     // Total count for this pre-syn unit
	nsyn += 1;
      }
      n = c->lgnin[i]->cn1;
      x = c->lgnin[i]->cx1;
      y = c->lgnin[i]->cy1;
      for(j=0;j<n;j++){
	xi = x[j];
	yi = y[j];
	fmap[xi][yi][1] += wp;    // Total weight for this pre-syn unit
	flag[xi][yi][1] += 1;     // Total count for this pre-syn unit
	nsyn += 1;
      }
    }
  }

  return nsyn;  // Newly marked units
}
/**************************************-**************************************/
/*                                                                           */
/*                          POPC_WF_ACCUM_LAY_TO_UNIT                        */
/*                                                                           */
/*  Mark all units in the flag array 'cmap' that send connections to the     */
/*  post-syn unit.                                                           */
/*                                                                           */
/*  Return the number of *NEWLY* marked units.                               */
/*                                                                           */
/*****************************************************************************/
int popc_wf_accum_lay_to_unit(lpre,lpost,ci,cj,ck,wp,flag,fmap,optype,fcrit)
     struct pop_layer *lpre;     // Pre-syn layer
     struct pop_layer *lpost;    // Post-syn layer
     int ci,cj,ck;               // Post-syn coords
     float wp;                   // prior weight for this unit 
     int   ***flag;              // [lpre->xn][->yn][->zn]
     float ***fmap;              // [lpre->xn][->yn][->zn]
     char *optype;               // "w_min", "w_max"
     float fcrit;                // Critical value, e.g., min weight
{
  int i,j,k;
  int nsyn,accept;
  struct pop_cell *c;
  struct pop_syn *tsyn;

  //printf("    wp %f  otype %s  fcrit %f\n",wp,optype,fcrit);

  nsyn = 0;  // Number of newly marked units

  if ((ci == -1)||(cj == -1)||(ck == -1))
    return nsyn;

  c = &(lpost->c[ci][cj][ck]); // Pointer to the post-syn cell

  // Get first synapse coming from lpre to 'c'
  tsyn = pop_cell_get_next_syn_layer_name(c->in,1,lpre->name); // 1-inputs
  while (tsyn != NULL){

    accept = 0;
    if (strcmp(optype,"w_min")==0){
      if (tsyn->w >= fcrit)
	accept = 1;
    }else if (strcmp(optype,"w_max")==0){
      if (tsyn->w <= fcrit)
	accept = 1;
    }else if (strcmp(optype,"all")==0){
      accept = 1;
    }else{
      printf("  *** optype:  %s\n",optype);
      exit_error("POPC_WF_ACCUM_LAY_TO_UNIT","Unknown optype");
    }

    if (accept == 1){
      i = tsyn->pre->layx;
      j = tsyn->pre->layy;
      k = tsyn->pre->layz;

      fmap[i][j][k] += wp * tsyn->w;  // Total weight for this pre-syn unit
      flag[i][j][k] += 1;             // Total count for this pre-syn unit

      nsyn += 1;
    }

    tsyn = tsyn->post_next;
    tsyn = pop_cell_get_next_syn_layer_name(tsyn,1,lpre->name); // postsyn
  }
  return nsyn;  // Newly marked units
}
/**************************************-**************************************/
/*                                                                           */
/*                          POPC_WF_ACCUM_LAY_TO_LAY                         */
/*                                                                           */
/*****************************************************************************/
float ***popc_wf_accum_lay_to_lay(lpre,lpost,xi,yi,zi,wprior,optype,fcrit)
     struct pop_layer *lpre;     // Pre-syn layer
     struct pop_layer *lpost;    // Post-syn layer
     int xi,yi,zi;               // Unit, or -1 for entire layer
     float ***wprior;            // prior weights for post-syn
     char *optype;               // Type of operation, criterion
     float fcrit;                // critical value
{
  int i,j,k;
  int nsyn,xn,yn,zn,***flag,lgnflag;
  float ***fmap,wp;

  printf("  POPC_WF_ACCUM_LAY_TO_LAY\n");
  printf("    %s   %s -> %s\n",optype,lpre->name,lpost->name);

  nsyn = 0;  // Number of newly marked units

  xn = lpost->xn;
  yn = lpost->yn;
  zn = lpost->zn;

  //  Maps are for the pre-syn layer
  flag = get_zero_3d_iarray(lpre->xn,lpre->yn,lpre->zn);
  fmap = get_zero_3d_farray(lpre->xn,lpre->yn,lpre->zn);
  //printf("    fmap is  %d,%d,%d\n",lpre->xn,lpre->yn,lpre->zn);

  //  Check if pre-syn is in post-syn LGN list
  if (strcmp(lpre->laytype,"lgn")==0){
    lgnflag = 1;
  }else
    lgnflag = 0;
  

  wp = 1.0;  // Use this value for all, if wprior is NULL

  if ((xi < 0)||(yi < 0)||(zi < 0)){
    //
    //  Entire layer
    //
    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){
	for(k=0;k<zn;k++){  //  For each post-syn unit (i,j,k)
	
	  if (wprior != NULL)
	    wp = wprior[i][j][k];

	  if (lgnflag == 1)
	    nsyn += popc_wf_accum_lgn_to_unit(lpre,lpost,i,j,k,wp,flag,fmap,
					      optype,fcrit);
	  else
	    nsyn += popc_wf_accum_lay_to_unit(lpre,lpost,i,j,k,wp,flag,fmap,
					      optype,fcrit);
	}
      }
    }
  }else{
    //
    //  Single unit
    //
    if (wprior != NULL)
      wp = wprior[xi][yi][zi];

    if (lgnflag == 1)
      nsyn += popc_wf_accum_lgn_to_unit(lpre,lpost,xi,yi,zi,wp,flag,fmap,
					optype,fcrit);
    else
      nsyn += popc_wf_accum_lay_to_unit(lpre,lpost,xi,yi,zi,wp,flag,fmap,
					optype,fcrit);
  }

  printf("    %d synapses total\n",nsyn);
  //printf("    done  POPC_WF_ACCUM_LAY_TO_LAY\n");

  free_3d_iarray(flag,xn,yn,zn);

  return fmap;  // Return the pre-syn weight map
}
/**************************************-**************************************/
/*                                                                           */
/*                           POPC_WF_ACCUM_MULTI_LAY                         */
/*                                                                           */
/*****************************************************************************/
float ***popc_wf_accum_multi_lay(laylist,nlay,xi,yi,zi,optype,fcrit)
     struct pop_layer **laylist;   // Sequence of layers [nlay]
     int nlay;                     // Number of layers in 'laylist'
     int xi,yi,zi;                 // Unit coord, or -1 for entire layer
     char *optype;                 // Type of operation, criterion
                                   // "aw" - aggregate weight
                                   // "an" - aggregate number of connections
     float fcrit;                  // critical value, min_w, or min_n
{
  float ***wprior,***w,***wp;
  char *opt;

  // *********** WYETH WYETH HERE **************  2011 April 19
  // *********** WYETH WYETH HERE **************
  //  what is this supposed to do?  Apply 'fcrit', or not?  It seems better
  //  to apply no condition, then let the caller apply the condition to the
  //  final aggregate?  WHAT IS intended here????
  // *********** WYETH WYETH HERE **************
  // *********** WYETH WYETH HERE **************

  printf("  POPC_WF_ACCUM_MULTI_LAY\n");
  printf("    %d %s\n",nlay,optype);

  if (nlay < 2)
    exit_error("POPC_W_ACCUM_MULTI_LAY","Too few layers");

  if (nlay > 2){
    //
    //  This is a multi-layer operation.  If the condition is "w_agg_", then
    //  we must accumulate the weights using "all" until the final stage.
    //
    if ((strcmp(optype,"w_agg_min")==0)||
	(strcmp(optype,"w_agg_max")==0)){
      opt = strdup("w_agg_all");
    }else{
      opt = strdup(optype);
    }

    wprior = popc_wf_accum_multi_lay(&(laylist[1]),nlay-1,xi,yi,zi,opt,fcrit);
    myfree(opt);

    if (strcmp(optype,"w_agg_min")==0){
      //opt = strdup("w_min");
      //  *** WYETH basically, the 'fcrit' is on the overall weight, so we
      //  *** WYETH basically, the 'fcrit' is on the overall weight, so we
      //  *** WYETH basically, the 'fcrit' is on the overall weight, so we
      //      should use 'all' and let the caller apply the criterion????
      opt = strdup("all");
      wp = wprior;
    }else if (strcmp(optype,"w_agg_all")==0){
      opt = strdup("all");
      wp = wprior;
    }else{
      opt = strdup(optype);  // e.g., "all"
      wp = NULL;
    }
    //printf("   _____1_____ fcrit = %f  (opt = %s)\n",fcrit,opt);
    w = popc_wf_accum_lay_to_lay(laylist[0],laylist[1],-1,-1,-1,wp,opt,fcrit);
    myfree(opt);

  }else{
    //
    //  This is a two-layer operation.
    //
    wp = NULL;  // No prior weights

    if (strcmp(optype,"w_agg_min")==0){
      opt = strdup("w_min");
    }else if (strcmp(optype,"w_agg_all")==0){
      opt = strdup("all");
    }else{
      printf("  *** optype:  %s\n",optype);
      exit_error("POPC_W_ACCUM_MULTI_LAY","Unexpected operation type");
    }

    //printf("   _____2_____ fcrit = %f  (opt = %s)\n",fcrit,opt);
    w = popc_wf_accum_lay_to_lay(laylist[0],laylist[1],xi,yi,zi,wp,opt,fcrit);
    myfree(opt);
  }

  //printf("  done  POPC_WF_ACCUM_MULTI_LAY\n");
  {
    char outfile[SLEN];

    sprintf(outfile,"zz.w.%d_%s_to_%s.3d",nlay,laylist[0]->name,
	    laylist[1]->name);

    //printf("    3D writing -  %s\n",outfile);

    write_3d_data_part(outfile,w,0,laylist[0]->xn,0,laylist[0]->yn,0,
		       laylist[0]->zn,4,2,1); // 1 write as txy, for show3d
  }

  return w;
}
/**************************************-**************************************/
/*                                                                           */
/*                         POPC_FLAG_ACCUM_LAY_TO_UNIT                       */
/*                                                                           */
/*  Mark all units in the flag array 'cmap' that send connections to the     */
/*  post-syn unit.                                                           */
/*                                                                           */
/*  Return the number of *NEWLY* marked units.                               */
/*                                                                           */
/*****************************************************************************/
int popc_flag_accum_lay_to_unit(lpre,lpost,ci,cj,ck,cmap,minw)
     struct pop_layer *lpre;     // Pre-syn layer
     struct pop_layer *lpost;    // Post-syn layer
     int ci,cj,ck;               // Post-syn coords
     int ***cmap;                // [lpre->xn][->yn][->zn]
     float minw;                 // Only flag weights >= this amount
{
  int i,j,k;
  int nsyn;
  struct pop_cell *c;
  struct pop_syn *tsyn;

  nsyn = 0;  // Number of newly marked units

  if ((ci == -1)||(cj == -1)||(ck == -1))
    return nsyn;

  c = &(lpost->c[ci][cj][ck]); // Pointer to the post-syn cell

  // Get first synapse coming from lpre to 'c'
  tsyn = pop_cell_get_next_syn_layer_name(c->in,1,lpre->name); // 1-inputs
  while (tsyn != NULL){

    if ((minw < 0.0) || (tsyn->w >= minw)){
      i = tsyn->pre->layx;
      j = tsyn->pre->layy;
      k = tsyn->pre->layz;

      if (cmap[i][j][k] == 0)  // Only count newly marked units
	nsyn += 1;

      cmap[i][j][k] += 1;
    }

    tsyn = tsyn->post_next;
    tsyn = pop_cell_get_next_syn_layer_name(tsyn,1,lpre->name); // postsyn
  }
  
  return nsyn;  // Newly marked units
}
/**************************************-**************************************/
/*                                                                           */
/*                         POPC_FLAG_ACCUM_LAY_TO_LAY                        */
/*                                                                           */
/*  Mark all units in the flag array 'cmap' that send connections to the     */
/*  post-syn layer.                                                          */
/*                                                                           */
/*  Return the number of *NEWLY* marked units.                               */
/*                                                                           */
/*****************************************************************************/
int popc_flag_accum_lay_to_lay(lpre,lpost,cmap,minw)
     struct pop_layer *lpre;     // Pre-syn layer
     struct pop_layer *lpost;    // Post-syn layer
     int ***cmap;                // [lpre->xn][->yn][->zn]
     float minw;                 // Only flag weights >= this amount
{
  int i,j,k;
  int nsyn,xn,yn,zn;

  nsyn = 0;  // Number of newly marked units

  xn = lpost->xn;
  yn = lpost->yn;
  zn = lpost->zn;

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	nsyn += popc_flag_accum_lay_to_unit(lpre,lpost,i,j,k,cmap,minw);
      }
    }
  }

  return nsyn;  // Newly marked units
}
/**************************************-**************************************/
/*                                                                           */
/*                          POPC_FLAG_GET_LAY_TO_ALL                         */
/*                                                                           */
/*  Return a flag array indicating the synapses from 'lpre' to all other     */
/*  layers.                                                                  */
/*                                                                           */
/*  Return the number of marked units.                                       */
/*                                                                           */
/*****************************************************************************/
int popc_flag_get_lay_to_all(lpre,laylist,nlay,rcmap,minw)
     struct pop_layer *lpre;       // Pre-syn layer
     struct pop_layer **laylist;   // All layers, might include 'lpre' [nlay]
     int nlay;                     // Number of layers in 'laylist'
     int ****rcmap;                // [lpre->xn][->yn][->zn]
     float minw;                   // Only flag weights >= this amount
{
  int i;
  int nsyn;
  int ***cmap;
  struct pop_layer *lpost;       // Post-syn layer

  cmap = get_zero_3d_iarray(lpre->xn,lpre->yn,lpre->zn);

  nsyn = 0;  // Number of newly marked units

  for(i=0;i<nlay;i++){
    lpost = laylist[i];  // for each layer
    if (lpost != lpre){
      if (popc_lay_count_input_from_lay(lpre,lpost) > 0){  // if relevant input
	nsyn += popc_flag_accum_lay_to_lay(lpre,lpost,cmap,minw);
      }
    }
  }

  *rcmap = cmap;

  return nsyn;  // Newly marked units
}
/**************************************-**************************************/
/*                                                                           */
/*                           POPC_FLAG_ACCUM_LAY_ORI                         */
/*                                                                           */
/*  Mark all units in the flag array 'cmap' that have acceptable 'ori'       */
/*  attribute values.                                                        */
/*                                                                           */
/*  Return the number of *NEWLY* marked units.                               */
/*                                                                           */
/*****************************************************************************/
int popc_flag_accum_lay_ori(lay,cmap,optype,ori_0,ori_1)
     struct pop_layer *lay;      // Pre-syn layer
     int ***cmap;                // [lpre->xn][->yn][->zn]
     char *optype;               // "ori180_ctr_dev"
     float ori_0;                // Center of ori range (deg)
     float ori_1;                // Maximum deviation from 'ori_0' (deg)
{
  int i,j,k;
  int xn,yn,zn,cnt,accept;
  float ori;
  struct pop_cell *c;

  xn = lay->xn;
  yn = lay->yn;
  zn = lay->zn;

  cnt = 0;  // Number of newly marked units
  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){  //  For each post-syn unit (i,j,k)

	c = &(lay->c[i][j][k]); // Pointer to unit
	ori = pop_cell_attrib_get_f(c,"ori");

	accept = 0;
	if (strcmp(optype,"ori180_ctr_dev")==0){
	  if (get_circular_diff(ori,ori_0,180.0) <= ori_1){
	    accept = 1; // This unit should be marked for saving
	  }
	}

	if (accept == 1){
	  if (cmap[i][j][k] == 0)  // Only count newly marked units
	    cnt += 1;
	  cmap[i][j][k] += 1;
	}
      }
    }
  }

  return cnt;  // Newly marked units
}
/**************************************-**************************************/
/*                                                                           */
/*                        POP_CELL_GET_LGN_INPUT_BY_NAME                     */
/*                                                                           */
/*  Return the pointer to the LGN inputs arising from the LGN layer named    */
/*  'name', or NULL if no inputs come from 'name'.                           */
/*                                                                           */
/*****************************************************************************/
struct pop_lgn_in *pop_cell_get_lgn_input_by_name(c,name)
     struct pop_cell *c;     // Cell to be scanned for LGN inputs
     char *name;             // LGN name to search for in LGN input list
{
  int i;
  int n;
  struct pop_lgn_in *ln;

  ln = NULL;

  n = c->lgn_n;
  for(i=0;i<n;i++){
    if (strcmp(name,c->lgnin[i]->lay->name)==0){
      if (ln == NULL)
	ln = c->lgnin[i];
      else
	exit_error("POP_CELL_GET_LGN_INPUT_BY_NAME","Redundant LGN input.");
    }
  }
  return ln;
}
/**************************************-**************************************/
/*                                                                           */
/*                        POP_CELL_GET_LGN_INPUT_BY_LAY                      */
/*                                                                           */
/*  Return the pointer to the LGN inputs arising from the 'lay', or NULL     */
/*  if no inputs come from 'lay'.                                            */
/*                                                                           */
/*****************************************************************************/
struct pop_lgn_in *pop_cell_get_lgn_input_by_lay(c,lay)
     struct pop_cell *c;        // Cell to be scanned for LGN inputs  
     struct pop_layer *lay;     // layer to search for in LGN input list   
{
  return pop_cell_get_lgn_input_by_name(c,lay->name);
}
/**************************************-**************************************/
/*                                                                           */
/*                     POP_CELL_GET_LGN_INPUT_PAIR_BY_LAY                    */
/*                                                                           */
/*  Return the pointer to the LGN inputs arising from the 'lay' that are     */
/*  PAIRED with the regular LGN inputs for this cell.                        */
/*                                                                           */
/*****************************************************************************/
struct pop_lgn_in *pop_cell_get_lgn_input_pair_by_lay(c,lay)
     struct pop_cell *c;        // Cell to be scanned for LGN pair inputs
     struct pop_layer *lay;     // layer to search for in LGN input list   
{
  int i;
  int n;
  struct pop_lgn_in *ln;

  ln = NULL;

  n = c->lgn_n;
  for(i=0;i<n;i++){
    if (c->lgnin[i]->cpair != NULL){  // If there are paired inputs
      //printf("pairlayername ====== %s    layname %s\n",
      //c->lgnin[i]->cpair->lay->name,lay->name);
      if (strcmp(lay->name,c->lgnin[i]->cpair->lay->name)==0){
	if (ln == NULL)
	  ln = c->lgnin[i]->cpair;
	else
	  exit_error("POP_CELL_GET_LGN_INPUT_PAIR_BY_LAY","Redundant input");
      }
    }
  }
    
  return ln;
}
/**************************************-**************************************/
/*                                                                           */
/*                           POP_CELL_ADD_LGN_INPUT                          */
/*                                                                           */
/*  NOTES                                                                    */
/*  - Coordinate arrays are attached by pointer, caller DO NOT FREE.         */
/*                                                                           */
/*****************************************************************************/
void pop_cell_add_lgn_input(c,lay,mech,cx0,cy0,cn0,cx1,cy1,cn1,nposs,weight)
     struct pop_cell *c;        // Cell to which connections are added
     struct pop_layer *lay;     // LGN layer from which connections arise
     struct pop_mech *mech;     // post-syn receptor mechanism
     int *cx0,*cy0;             // OFF coords
     int cn0;                   // number of OFF connections
     int *cx1,*cy1;             // ON coords
     int cn1;                   // number of ON connections
     int nposs;                 // Possible LGN inputs, given some template
     float weight;              // Constant weight for each synapse
{
  int i;
  int n;
  struct pop_lgn_in **ln;

  n = c->lgn_n;  // Number of LGN inputs before adding this one

  // Make new list of LGN inputs, and copy existing LGN inputs to new list
  ln = (struct pop_lgn_in **)myalloc((n+1)*sizeof(struct pop_lgn_in *));
  for(i=0;i<n;i++)
    ln[i] = c->lgnin[i];

  ln[n] = (struct pop_lgn_in *)myalloc(sizeof(struct pop_lgn_in));

  ln[n]->lay = lay;
  ln[n]->mech = mech;
  ln[n]->cn0 = cn0;  ln[n]->cx0 = cx0;  ln[n]->cy0 = cy0;
  ln[n]->cn1 = cn1;  ln[n]->cx1 = cx1;  ln[n]->cy1 = cy1;
  ln[n]->nposs = nposs;
  ln[n]->gw = weight;  // Default global weight  WYETH LGNW

  ln[n]->cpair = NULL;

  if ((c->lgnin != NULL) && (c->lgn_n > 0))
    myfree(c->lgnin);  // Free the old list

  c->lgnin = ln;    // Attach the new list
  c->lgn_n = n+1;  // Update number of LGN inputs
}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_CELL_ADD_LGN_PAIR                          */
/*                                                                           */
/*  NOTES                                                                    */
/*  - Coordinate arrays are attached by pointer, caller DO NOT FREE.         */
/*                                                                           */
/*****************************************************************************/
void pop_cell_add_lgn_pair(lgnin,lay,mech,cx0,cy0,cn0,cx1,cy1,cn1,nposs)
     struct pop_lgn_in *lgnin;  // Existing LGN inputs to be paired with
     struct pop_layer *lay;     // LGN layer from which connections arise
     struct pop_mech *mech;     // post-syn receptor mechanism
     int *cx0,*cy0;             // OFF coords
     int cn0;                   // number of OFF connections
     int *cx1,*cy1;             // ON coords
     int cn1;                   // number of ON connections
     int nposs;                 // Possible LGN inputs, given some template
{
  struct pop_lgn_in *ln;

  /*** WYETH - do we need to add 'weight' as an input (as above) ? ***/

  ln = (struct pop_lgn_in *)myalloc(sizeof(struct pop_lgn_in));

  ln->lay = lay;
  ln->mech = mech;

  ln->cn0 = cn0;
  ln->cx0 = cx0;
  ln->cy0 = cy0;

  ln->cn1 = cn1;
  ln->cx1 = cx1;
  ln->cy1 = cy1;

  ln->nposs = nposs;
  ln->gw = 1.0;  // Default global weight WYETH LGNW
  ln->cpair = NULL;

  lgnin->cpair = ln;    // Attach this to the existing LGN input
}
/**************************************-**************************************/
/*                                                                           */
/*                          POP_CELL_GET_MECH_BY_NAME                        */
/*                                                                           */
/*****************************************************************************/
struct pop_mech *pop_cell_get_mech_by_name(lay,name)
     struct pop_layer *lay;
     char *name;
{
  int i;
  struct pop_mech *pm;
  char *tname;
  
  pm = NULL;
  for(i=0;i<lay->nmech;i++){
    tname = onode_getpar_chr_ptr_exit(lay->mech[i]->o,"name");
    if (strcmp(tname,name)==0){
      if (pm == NULL)
	pm = lay->mech[i];
      else
	exit_error("POP_CELL_GET_MECH_BY_NAME","Duplicate mech names");
    }
  }
  return pm;
}
/**************************************-**************************************/
/*                                                                           */
/*                           POP_CELL_LAYER_ADD_MECH                         */
/*                                                                           */
/*****************************************************************************/
void pop_cell_layer_add_mech(mylogf,lay,mech)
     char *mylogf;
     struct pop_layer *lay;
     struct pop_mech *mech;
{
  int i;
  int n;
  struct pop_mech *tm,**tmlist;
  char tstr[SLEN];

  tm = pop_cell_get_mech_by_name(lay,mech->name);
  if (tm != NULL){
    sprintf(tstr,"  Mechanism name  %s\n",mech->name);
    mylog(mylogf,tstr);
    mylog_exit(mylogf,"POP_CELL_LAYER_ADD_MECH  Mechanism already exists\n");
  }

  /* Copy any existing SIDs for this layer */
  n = 1 + lay->nmech;
  tmlist = (struct pop_mech **)myalloc(n*sizeof(struct pop_mech *));
  for(i=0;i<lay->nmech;i++)
    tmlist[i] = lay->mech[i];
  tmlist[n-1] = mech; /* Add the new one to the end of the list */

  myfree(lay->mech);
  lay->mech = tmlist;
  lay->nmech = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                   POP_CELL_CELLIST_GET_PRESYN_FROM_LAYER                  */
/*                                                                           */
/*  Return a list of pointers to all cells in layer 'layname' that synpapse  */
/*  on 'c'.                                                                  */
/*                                                                           */
/*****************************************************************************/
struct pop_cell **pop_cell_cellist_get_presyn_from_layer(layname,c,rn)
     char layname[];
     struct pop_cell *c;
     int *rn;
{
  int k;
  int n;
  struct pop_cell **cl;
  struct pop_syn *s,*t;

  n = pop_cell_count_syn_from_layer_to_cell(layname,c);

  cl = (struct pop_cell **)myalloc(n*sizeof(struct pop_cell *));

  // Scan pre-syn list for synapses from 'lay'
  s = pop_cell_get_next_syn_layer_name(c->in,1,layname);
  k = 0;
  while(s != NULL){
    cl[k] = s->pre;  // The pre-syn cell is the one sending to 'c'
    k += 1;
    t = s->post_next; // The post-syn cell of an input syn is this cell
    s = pop_cell_get_next_syn_layer_name(t,1,layname);
  }

  *rn = n;
  return cl;
}
/**************************************-**************************************/
/*                                                                           */
/*                      POP_CELL_RESET_ALL_OUT_SYNAPSES                      */
/*                                                                           */
/*  Reset synaptic depression params.                                        */
/*                                                                           */
/*****************************************************************************/
void pop_cell_reset_all_out_synapses(c)
     struct pop_cell *c;
{
  int i;
  struct pop_syn *t;

  t = c->out;
  i = 1;
  while(t != NULL){

    //printf("Pre/PostCell    %s   %s\n",t->pre->name,t->post->name);

    t->tlast = -1000000.0; // Previous spike was a long time ago
    t->alast = 1.0;        // Full efficacy

    if (t->si != NULL)
      if (t->si->ci == 0)  // Only the 0-cell does the reset
	pop_cell_si_reset(t);

    t = t->pre_next;       // Traverse synapses of same pre-syn cell
    i += 1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_CELL_SYN_SI_FIND                           */
/*                                                                           */
/*  Informally, this routine finds the "other half" of an SI.                */
/*  Find the synapse that uses 'siu' and has SI index 'sii'.                 */
/*                                                                           */
/*****************************************************************************/
struct pop_syn *pop_cell_syn_si_find(mylogf,s,siu,sii)
     char *mylogf;
     struct pop_syn *s;        // pointer to pre-syn list
     union si_union *siu;      // Pointer to union
     int sii;                  // SI index
{
  struct pop_syn *t,*syn;

  syn = NULL;

  // Scan all synapses in this input list
  t = s;
  while(t != NULL){

    //printf("pre,post = %s %s\n",t->pre->name,t->post->name);
    
    if (t->si != NULL){         // If there is an SI
      if (t->si->siu == siu){   //   If the union is the one sought
	if (t->si->ci == sii){  //      If the index is the one sought
	  syn = t;   // save this to return
	  t = NULL;  // signal that we are done
	  //printf("FOUND SI match!!!!!!!!!!!!!!\n");
	}
      }
    }
    
    if (t != NULL)
      t = t->post_next;  // Next syn that has same POST-syn cell
  }

  return syn;
}
/**************************************-**************************************/
/*                                                                           */
/*                             POP_CELL_ADD_SYNAPSE                          */
/*                                                                           */
/*  The returned pointer can be ignored, because the synapse is connected    */
/*  into a list.                                                             */
/*                                                                           */
/*****************************************************************************/
struct pop_syn *pop_cell_add_synapse(pre,post,stype,w,tdelay,inindex)
     struct pop_cell *pre,*post;
     int stype;      // type, specified in 'ifc.h'
     float w;        // weight, 0-1
     float tdelay;
     short inindex;  // Index into population 'inlist'
{
  struct pop_syn *s;

  //printf("  FIRST:  inindx = %d\n",inindex);
  //printf("          stype = %d\n",stype);
  //printf("          w = %f\n",w);
  //printf("          tdelay = %f\n",tdelay);

  s = (struct pop_syn *)myalloc(sizeof(struct pop_syn));
  s->pre = pre;
  s->post = post;

  s->stype = stype;
  s->w = w;
  // Used for synaptic depression
  s->tlast = -1000000.0;  // Time of last spike (far in the past)
  s->alast = 1.0;         // Ampl at time of last spike, 1 is full
  s->tdelay = tdelay;

  s->si = NULL;  // NULL indicates there is no synaptic interaction

  s->inindx = inindex;

  pre->nout += 1;
  post->nin += 1;

  // Add to pre-synaptic list
  s->pre_next = pre->out;      // New synapse points to head of old 'out' list
  pre->out = s;                // Presyn out list now starts with new synapse
  if (s->pre_next != NULL)
    s->pre_next->pre_prev = s; // Next in list points back at new synapse
  s->pre_prev = NULL;          // There is no synapse before this one

  // Add to post-synaptic list
  s->post_next = post->in;     // New synapse points to head of old 'in' list
  post->in = s;                // Postsyn in list now starts with new synapse
  if (s->post_next != NULL){
    s->post_next->post_prev = s; // Next in list points back at new synapse
  }
  s->post_prev = NULL;         // There is no synapse before this one

  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                           POP_CELL_REMOVE_SYNAPSE                         */
/*                                                                           */
/*****************************************************************************/
void pop_cell_remove_synapse(s)
     struct pop_syn *s;
{
  struct pop_syn *prev,*next;
  
  // Remove from pre-synaptic list
  prev = s->pre_prev;
  next = s->pre_next;
  if (prev != NULL)
    prev->pre_next = next; // Make previous point to next
  else
    s->pre->out = next; // 'next' is the new head of 'out' list

  if (next != NULL)
    next->pre_prev = prev; // Make next to back to prev

  s->pre->nout -= 1;  // Adjust output count for pre-syn cell


  // Remove from post-synaptic list
  prev = s->post_prev;
  next = s->post_next;
  if (prev != NULL)
    prev->post_next = next; // Make previous point to next
  else
    s->post->in = next;  // 'next' is the new head of 'in' list

  if (next != NULL)
    next->post_prev = prev; // Make next to back to prev

  s->post->nin -= 1;  // Adjust input count for post-syn cell

  pop_cell_syn_free(s);
}
/**************************************-**************************************/
/*                                                                           */
/*                      POP_CELL_REMOVE_ALL_SYN_FROM_CELL                    */
/*                                                                           */
/*  Remove all synapses listed for this cell.                                */
/*                                                                           */
/*****************************************************************************/
void pop_cell_remove_all_syn_from_cell(mylogf,c)
     char *mylogf;
     struct pop_cell *c;
{
  int cnt;
  struct pop_syn *s,*t;

  cnt = 0;
  while(c->out != NULL){
    pop_cell_remove_synapse(c->out);
    cnt += 1;
  }
  c->nout = 0;
  //printf("OUT cnt = %d\n",cnt);

  cnt = 0;
  while(c->in != NULL){
    pop_cell_remove_synapse(c->in);
    cnt += 1;
  }
  c->nin = 0;
  //printf("IN cnt = %d\n",cnt);
}
/**************************************-**************************************/
/*                                                                           */
/*                 POP_CELL_REMOVE_ALL_SYN_FROM_LAYER_TO_CELL                */
/*                                                                           */
/*****************************************************************************/
void pop_cell_remove_all_syn_from_layer_to_cell(mylogf,lay,c)
     char *mylogf;
     struct pop_layer *lay;
     struct pop_cell *c;
{
  int n;
  struct pop_syn *s,*t;
  char tstr[SLEN];

  sprintf(tstr,"    * Removing all synapses from layer '%s' to cell %s\n",
	  lay->name,c->name);
  mylog(mylogf,tstr);

  if (strcmp(lay->name,"lgn_on")==0){
    mylog_exit(mylogf,"***WYETH - MUST REWRITE REMOVE FOR LGN");
    /*
      n = c->cn1;
      if (n > 0){
      c->cn1 = 0;
      myfree(c->cx1);
      myfree(c->cy1);
      }*/
  }else if (strcmp(lay->name,"lgn_off")==0){
    mylog_exit(mylogf,"***WYETH - MUST REWRITE REMOVE FOR LGN");
    /*
      n = c->cn0;
      if (n > 0){
      c->cn0 = 0;
      myfree(c->cx0);
      myfree(c->cy0);
      }
    */
  }else{
    // Scan pre-syn list and remove all synapses from 'lay'
    s = pop_cell_get_next_syn_layer_name(c->in,1,lay->name);
    n = 0;
    while(s != NULL){
      t = s->post_next; // The post-syn cell of an input syn is this cell
      pop_cell_remove_synapse(s);
      n += 1;
      s = pop_cell_get_next_syn_layer_name(t,1,lay->name);
    }
  }
  sprintf(tstr,"      Removed %d synapses\n",n);
  mylog(mylogf,tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                 POP_CELL_REMOVE_ALL_SYN_FROM_CELL_TO_LAYER                */
/*                                                                           */
/*****************************************************************************/
void pop_cell_remove_all_syn_from_cell_to_layer(mylogf,c,lay)
     char *mylogf;
     struct pop_cell *c;
     struct pop_layer *lay;
{
  int n,sii_other;
  struct pop_syn *s,*t,*so;
  char tstr[SLEN];

  sprintf(tstr,"    * Removing all synapses from cell '%s' to layer %s\n",
	  c->name,lay->name);
  mylog(mylogf,tstr);

  // Scan pre-syn list and remove all synapses from 'lay'
  s = pop_cell_get_next_syn_layer_name(c->out,0,lay->name);
  n = 0;
  while(s != NULL){
    t = s->pre_next; // The post-syn cell of an input syn is this cell

    // REMOVE RELATED SI's  *** WYETH SHOULD MAKE A ROUTINE FOR THIS
    if (s->si != NULL){        // There is an SI

      // Find the other synapse involved in this interaction
      if (s->si->siui == 1){ // This SI type only has 2 index values, 0 and 1
	sii_other = 1 - s->si->ci;
	so = pop_cell_syn_si_find(mylogf,s->post->in,s->si->siu,sii_other);
	if (so == NULL)
	  mylog_exit(mylogf,"POP_CELL_REMOVE_ALL_SYN_FROM_CELL_TO_LAYER",
		     "Could not find paired synapse for SI");
	pop_cell_remove_synapse(so);
      }else{
	mylog_exit(mylogf,"POP_CELL_REMOVE_ALL_SYN_FROM_CELL_TO_LAYER",
		   "What about paired synapses for this SI type?");
      }
    }

    pop_cell_remove_synapse(s);
    n += 1;
    s = pop_cell_get_next_syn_layer_name(t,0,lay->name);
  }
  sprintf(tstr,"      Removed %d synapses\n",n);
  mylog(mylogf,tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*               POP_CELL_REMOVE_ALL_SYN_FROM_CELL_TO_LAYER_SI               */
/*                                                                           */
/*****************************************************************************/
void pop_cell_remove_all_syn_from_cell_to_layer_si(mylogf,c,lay,sii)
     char *mylogf;
     struct pop_cell *c;
     struct pop_layer *lay;
     int sii;                  // SI index
{
  int n,sii_other;
  struct pop_syn *s,*t,*so;
  char tstr[SLEN];

  sprintf(tstr,"    * Removing all SI %d syn. from cell '%s' to layer %s\n",
	  sii,c->name,lay->name);
  mylog(mylogf,tstr);

  // Scan pre-syn list and remove all synapses from 'lay'
  s = pop_cell_get_next_syn_layer_name(c->out,0,lay->name);
  n = 0;
  while(s != NULL){
    t = s->pre_next; // The pre-syn cell of an output syn is this cell

    if (s->si != NULL){        // There is an SI
      if (s->si->ci == sii){      // The SI index matches

	// Find the other synapse involved in this interaction
	if (s->si->siui == 1){ // This SI type only has 2 index values, 0 and 1
	  sii_other = 1 - sii;
	  so = pop_cell_syn_si_find(mylogf,s->post->in,s->si->siu,sii_other);
	  if (so == NULL)
	    mylog_exit(mylogf,"POP_CELL_REMOVE_ALL_SYN_FROM_CELL_TO_LAYER_SI",
		       "Could not find paired synapse for SI");
	  pop_cell_remove_synapse(so);
	}else{
	    mylog_exit(mylogf,"POP_CELL_REMOVE_ALL_SYN_FROM_CELL_TO_LAYER_SI",
		       "What about paired synapses for this SI type?");
	}

	// Now remove this synapse
	pop_cell_remove_synapse(s);

	n += 1;
      }
    }
    s = pop_cell_get_next_syn_layer_name(t,0,lay->name);
  }
  sprintf(tstr,"      Removed %d synapses\n",n);
  mylog(mylogf,tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                 POP_CELL_ADJUST_ALL_SYN_FROM_CELL_TO_LAYER                */
/*                                                                           */
/*****************************************************************************/
void pop_cell_adjust_all_syn_from_cell_to_layer(mylogf,c,lay,synw)
     char *mylogf;
     struct pop_cell *c;
     struct pop_layer *lay;
     float synw;
{
  int n,sii_other;
  struct pop_syn *s,*t,*so;
  char tstr[SLEN];

  sprintf(tstr,"    * Adjusting all synapses from cell '%s' to layer %s\n",
	  c->name,lay->name);
  mylog(mylogf,tstr);

  // Scan pre-syn list and remove all synapses from 'lay'
  s = pop_cell_get_next_syn_layer_name(c->out,0,lay->name);
  n = 0;
  while(s != NULL){
    t = s->pre_next; // The post-syn cell of an input syn is this cell

    // REMOVE RELATED SI's  *** WYETH SHOULD MAKE A ROUTINE FOR THIS
    if (s->si == NULL){        // No SI
      s->w = synw;
    }else{       // There is an SI

      // Find the other synapse involved in this interaction
      if (s->si->siui == 1){ // This SI type only has 2 index values, 0 and 1
	if (s->si->ci == 0){
	  s->w = synw;
	  printf("WYETH___________  ci = 0\n");
	}else{
	  ;
	  /**** DONT ADJUST?
	  printf("WYETH___________  ci = 1\n");
	  sii_other = 0;
	  so = pop_cell_syn_si_find(mylogf,s->post->in,s->si->siu,sii_other);
	  if (so == NULL)
	    mylog_exit(mylogf,"POP_CELL_ADJUST_ALL_SYN_FROM_CELL_TO_LAYER",
		       "Could not find paired synapse for SI");
	  so->w = synw;
	  */
	}
      }else{
	mylog_exit(mylogf,"POP_CELL_ADJUST_ALL_SYN_FROM_CELL_TO_LAYER",
		   "What about paired synapses for this SI type?");
      }
    }

    n += 1;
    s = pop_cell_get_next_syn_layer_name(t,0,lay->name);
  }
  sprintf(tstr,"      Adjusted %d synapses\n",n);
  mylog(mylogf,tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                  POP_CELL_SAVE_ALL_SYN_FROM_LAYER_TO_CELL                 */
/*                                                                           */
/*  Modify what data is saved for all cells from 'lay' to 'c'.               */
/*                                                                           */
/*****************************************************************************/
void pop_cell_save_all_syn_from_layer_to_cell(mylogf,layname,c)
     char *mylogf;
     char layname[];
     struct pop_cell *c;
{
  int n;
  struct pop_syn *s,*t;
  struct pop_cell *pre;

  // Scan pre-syn list and remove all synapses from 'lay'
  s = pop_cell_get_next_syn_layer_name(c->in,1,layname);
  while(s != NULL){
    pre = s->pre; // The pre-synaptic cell
    pre->savs = 1;        // Save spikes
    pre->maxscnt = 1000;  // Will be doubled as needed on the fly

    t = s->post_next; // The post-syn cell of an input syn is 'c'
    s = pop_cell_get_next_syn_layer_name(t,1,layname);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                 POP_CELL_COPY_SYN_LAYER_FROM_CELL_TO_CELL                 */
/*                                                                           */
/*  Copy certain synapses between layer 'lay' and cell 'c1' to cell 'c2'.    */
/*    pflag 0:  copy lay->c1 synapses (input synapses)                       */
/*    pflag 1:  copy c1->lay synapses (output synapses)                      */
/*                                                                           */
/*****************************************************************************/
void pop_cell_copy_syn_layer_from_cell_to_cell(mylogf,lay,c1,c2,pflag)
     char *mylogf;
     struct pop_layer *lay;
     struct pop_cell *c1;
     struct pop_cell *c2;
     int pflag;               // 0: lay->c1, 1: c1->lay
{
  int n;
  struct pop_syn *s,*t;
  char tstr[SLEN];

  n = 0;
  if (pflag == 0){
    sprintf(tstr,"    * Copying all synapses '%s' -> %s to %s\n",
	    lay->name,c1->name,c2->name);
    mylog(mylogf,tstr);

    // Scan input list for synapses from 'lay'
    s = pop_cell_get_next_syn_layer_name(c1->in,1,lay->name);
    while(s != NULL){
      (void)pop_cell_add_synapse(s->pre,c2,s->stype,s->w,0.0,-1);
      n += 1;
      t = s->post_next; // The post-syn cell of an input syn is this cell
      s = pop_cell_get_next_syn_layer_name(t,1,lay->name);
    }
    

  }else{
    sprintf(tstr,"    * Copying all synapses '%s' -> %s to %s\n",
	    c1->name,lay->name,c2->name);
    mylog(mylogf,tstr);
    
    // Scan output list for synapses to 'lay'
    s = pop_cell_get_next_syn_layer_name(c1->out,0,lay->name);
    while(s != NULL){
      (void)pop_cell_add_synapse(c2,s->post,s->stype,s->w,0.0,-1);
      n += 1;
      t = s->pre_next; // The pre-syn cell of an output syn is this cell
      s = pop_cell_get_next_syn_layer_name(t,0,lay->name);
    }
  }
  sprintf(tstr,"      Copied %d synapses\n",n);
  mylog(mylogf,tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                        POP_CELL_NORM_WEIGHT_FROM_LAYER                    */
/*                                                                           */
/*  Normalize the total weight of all connections from the layer to the      */
/*  cell to be the given total weight.                                       */
/*                                                                           */
/*****************************************************************************/
float pop_cell_norm_weight_from_layer(mylogf,lay,c,wf)
     char *mylogf;
     struct pop_layer *lay; // layer
     struct pop_cell *c;    // cell
     float wf;              // desired total weight, 1.0 is unitary strength
{
  float wtot,ww;
  struct pop_syn *syn;

  wtot = pop_cell_tot_weight_from_layer_to_cell(lay->name,c);
  ww = wf/wtot;

  syn = pop_cell_get_next_syn_layer_name(c->in,1,lay->name); // inputs
  while(syn != NULL){
    syn->w *= ww;
    syn = syn->post_next;
    syn = pop_cell_get_next_syn_layer_name(syn,1,lay->name);
  }

  return wtot;
}
/**************************************-**************************************/
/*                                                                           */
/*                     POP_CELL_NORM_WEIGHT_FROM_ININDEX                     */
/*                                                                           */
/*  Normalize the total weight of all inputs having 'inindex' to the cell    */
/*  to be 'wf'.                                                              */
/*                                                                           */
/*****************************************************************************/
float pop_cell_norm_weight_from_inindex(mylogf,c,inindex,wf)
     char *mylogf;
     //struct pop_layer *lay; // layer
     struct pop_cell *c;    // cell
     short inindex;
     float wf;              // desired total weight, 1.0 is unitary strength
{
  float wtot,ww;
  struct pop_syn *syn;

  //wtot = pop_cell_tot_weight_from_layer_to_cell(lay->name,c);
  wtot = pop_cell_tot_weight_inindex_to_cell(c,inindex);
  ww = wf/wtot;

  syn = pop_cell_get_next_syn_index(c->in,1,inindex,NULL);
  //syn = pop_cell_get_next_syn_layer_name(c->in,1,lay->name); // inputs
  while(syn != NULL){
    syn->w *= ww;
    syn = syn->post_next;
    syn = pop_cell_get_next_syn_index(syn,1,inindex,NULL);
    //syn = pop_cell_get_next_syn_layer_name(syn,1,lay->name);
  }

  return wtot;
}
/**************************************-**************************************/
/*                                                                           */
/*                       POP_CELL_NORM_WEIGHT_LAYER_LAYER                    */
/*                                                                           */
/*  For each cell in lay2, normalize total weight from lay1 to target value. */
/*                                                                           */
/*****************************************************************************/
void pop_cell_norm_weight_layer_layer(mylogf,lay1,lay2,wf,inindex)
     char *mylogf;
     struct pop_layer *lay1; // pre-syn layer
     struct pop_layer *lay2; // post-syn layer
     float wf;               // desired total weight, 1.0 is unitary strength
     short inindex;          // Index into population 'inlist'
{
  int i,j,k;
  int xn,yn,zn;
  float old_totw;
  struct pop_cell *c;

  xn = lay2->xn;
  yn = lay2->yn;
  zn = lay2->zn;

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	c = &(lay2->c[i][j][k]);
	old_totw = pop_cell_norm_weight_from_layer(mylogf,lay1,c,wf);
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                     POP_CELL_NORM_WEIGHT_LAYER_ININDEX                    */
/*                                                                           */
/*  For each cell in lay2, normalize total weight for 'inindex' inputs.      */
/*                                                                           */
/*****************************************************************************/
void pop_cell_norm_weight_layer_inindex(mylogf,lay,wf,inindex)
     char *mylogf;
     struct pop_layer *lay;  // post-syn layer
     float wf;               // desired total weight, 1.0 is unitary strength
     short inindex;          // Index into population 'inlist'
{
  int i,j,k;
  int xn,yn,zn;
  float old_totw;
  struct pop_cell *c;

  xn = lay->xn;
  yn = lay->yn;
  zn = lay->zn;

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      for(k=0;k<zn;k++){
	c = &(lay->c[i][j][k]);
	old_totw = pop_cell_norm_weight_from_inindex(mylogf,c,inindex,wf);
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                      POPC_LAYIRR_GET_CONDITIONAL_COORDS                   */
/*                                                                           */
/*  Return the subset of the irregular coordinates that satisfy the          */
/*  specified condition.                                                     */
/*                                                                           */
/*****************************************************************************/
void popc_layirr_get_conditional_coords(lay,cond,rx,ry,rn)
     struct pop_layer *lay;  // Layer w/ irregular geometry
     char *cond;             // condition to apply
     float **rx;             // [*rn] x-coords
     float **ry;             // [*rn] y-coords
     int *rn;                // Number of qualified points
{
  int i,k;
  float *x,*y;
  int tid;

  //
  //  *** WYETH - THIS MIGHT NOT BE USED?  I ended up using the routine below.
  //

  if (strcmp(cond,"all")==0)
    tid = -1;
  else if (strcmp(cond,"0")==0)
    tid = 0;
  else if (strcmp(cond,"1")==0)
    tid = 1;
  else if (strcmp(cond,"2")==0)
    tid = 2;
  else{
    printf("  Condition:  %s\n",cond);
    exit_error("POPC_LAYIRR_GET_CONDITIONAL_COORDS","Unknown condition");
  }

  x = get_farray(lay->irr_n);
  y = get_farray(lay->irr_n);

  k = 0;
  for(i=0;i<lay->irr_n;i++){
    if ((tid == -1) || (lay->irr_id[i] == tid)){
      x[k] = lay->irr_x[i];
      y[k] = lay->irr_y[i];
      k += 1;
    }
  }

  *rx = x;
  *ry = y;
  *rn = k;
}
/**************************************-**************************************/
/*                                                                           */
/*                      POPC_LAYIRR_GET_CONDITIONAL_TFLAG                    */
/*                                                                           */
/*  Return the subset of the irregular coordinates that satisfy the          */
/*  specified condition.                                                     */
/*                                                                           */
/*****************************************************************************/
int **popc_layirr_get_conditional_tflag(lay,cond)
     struct pop_layer *lay;  // Layer w/ irregular geometry
     char *cond;             // condition to apply
{
  int i;
  int n,tid,**tflag;

  if (strcmp(cond,"all")==0)
    tid = -1;
  else if (strcmp(cond,"0")==0)
    tid = 0;
  else if (strcmp(cond,"1")==0)
    tid = 1;
  else if (strcmp(cond,"2")==0)
    tid = 2;
  else{
    printf("  Condition:  %s\n",cond);
    exit_error("POPC_LAYIRR_GET_CONDITIONAL_TFLAG","Unknown condition");
  }

  n = lay->irr_n;
  tflag = get_const_2d_iarray(n,1,1);  // Set to all 1's

  if (tid != -1){  // If this is not "all"

    for(i=0;i<n;i++){
      if (lay->irr_id[i] != tid){
	tflag[i][0] = 0;
      }
    }

  }

  return tflag;
}
/**************************************-**************************************/
/*                                                                           */
/*                        POP_CELL_ADD_ON_OFF_INPUT_IRR                      */
/*                                                                           */
/*****************************************************************************/
void pop_cell_add_on_off_input_irr(c,lay,cx0,cy0,cn0,cx1,cy1,cn1,nposs,
				   inindex)
     struct pop_cell *c;        // Cell to which connections are added
     struct pop_layer *lay;     // LGN layer from which connections arise
     int *cx0,*cy0;             // OFF coords
     int cn0;                   // number of OFF connections
     int *cx1,*cy1;             // ON coords
     int cn1;                   // number of ON connections
     int nposs;                 // Possible LGN inputs, given some template
     short inindex;             // Input index number
{
  int i;
  int xi,syn_type;
  float w,tdelay;
  struct pop_cell *cpre;

  w = 1.0;        // Weight is constant
  tdelay = 0.0;
  syn_type = 11;  // stype 11 indicates a pre-computed input, lgn-like

  for(i=0;i<cn0;i++){
    xi = cx0[i];
    cpre = &(lay->c[xi][0][0]);  //  0 - OFF
    if (cpre == NULL)
      exit_error("POP_CELL_ADD_ON_OFF_INPUT_IRR","OFF cell is NULL");

    (void)pop_cell_add_synapse(cpre,c,syn_type,w,tdelay,inindex);
  }

  for(i=0;i<cn1;i++){
    xi = cx1[i];
    cpre = &(lay->c[xi][0][1]);  //  1 - ON
    if (cpre == NULL)
      exit_error("POP_CELL_ADD_ON_OFF_INPUT_IRR","ON cell is NULL");
    
    (void)pop_cell_add_synapse(cpre,c,syn_type,w,tdelay,inindex);
  }

  //
  //  WYETH HERE - I needed to set this above '0' so responses could be
  //     req'd for 'lgn_ga' for a V1 cell getting input
  //
  // NOTE, 'c->lgnin' will remain NULL
  c->lgn_n += 1;  // Update number of LGN inputs WYETH HERE

}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_CELL_SYN_SET_SAVE                          */
/*                                                                           */
/*****************************************************************************/
void pop_cell_syn_set_save(mylogf,s,datid)
     char *mylogf;
     struct pop_syn *s;
     char *datid;
{
  mylog(mylogf,"  POP_CELL_SYN_SET_SAVE\n");

  if ((strcmp(datid,"si_mask")==0) ||
      (strcmp(datid,"si_mask1")==0) ||
      (strcmp(datid,"si_mask2")==0)){

    pop_cell_si_add_sav(mylogf,s,datid);

  }else{
    mylog_exit(mylogf,"*** POP_CELL_SYN_SET_SAVE  Invalid data ID\n");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              POP_CELL_SET_SAVE                            */
/*                                                                           */
/*****************************************************************************/
void pop_cell_set_save(mylogf,c,datid,nmdaflag)
     char *mylogf;
     struct pop_cell *c;
     char *datid;
     int nmdaflag;
{
  int csi;
  char tstr[SLEN];

  if (compare_prefix_string_order("all_spk_",datid)){
    exit_error("POP_CELL_SET_SAVE","This doesn't handle 'all_spk_'");
  }else if (strcmp(datid,"spikes")==0){
    c->savs = 1;        // Save spikes
    c->maxscnt = 1000;  // Will be doubled as needed on the fly
  }else if ((strcmp(datid,"vm")==0)||(strcmp(datid,"rate")==0)){
    c->savv = 1;        // Save Vm values
    c->dxsav = 0.0004;  // WYETH- Save values no more often than this (s)
  }else if (strcmp(datid,"vd")==0){
    c->savv = 1;        // Save Vm values
    c->savd = 1;        // Save Vm Dendr values
    c->dxsav = 0.0004;  // WYETH- Save values no more often than this (s)
  }else if (strcmp(datid,"gad")==0){
    c->sava = 1;        // Save g_ad values
  }else if ((strcmp(datid,"lgn_gx")==0)||  // AMPA + NMDA
	    (strcmp(datid,"lgn_gain")==0)||
	    (strcmp(datid,"lgn_ga")==0)||
	    (strcmp(datid,"lgn_ia")==0)||
	    (strcmp(datid,"lgn_gn")==0)||
	    (strcmp(datid,"gad")==0)||
	    (strcmp(datid,"ex_ix")==0)||
	    (strcmp(datid,"ex_gx")==0)||
	    (strcmp(datid,"in_ii")==0)||
	    (strcmp(datid,"in_gi")==0)){

    // Prevent requests for LGN-related data if no LGN connections
    if ((strcmp(datid,"lgn_gx")==0)||
	(strcmp(datid,"lgn_gain")==0)||
	(strcmp(datid,"lgn_ga")==0)||
	(strcmp(datid,"lgn_ia")==0)||
	(strcmp(datid,"lgn_gn")==0)){
      if (c->lgn_n == 0){
	sprintf(tstr,"*** Cell name:  %s\n",c->name);
	mylog(mylogf,tstr);
	mylog(mylogf,"*** Cannot save LGN-related data\n");
	mylog_exit(mylogf,"POP_CELL_SET_SAVE  No LGN input to cell\n");
      }
    }

    c->savg = 1;
    c->gtex1max = 2000000; // WYETH - trouble at 40k, 50k w/ pop02
    c->gtex1  = (float *)myalloc(c->gtex1max * sizeof(float));
    c->gtin1  = (float *)myalloc(c->gtex1max * sizeof(float));
    c->gtex1t = (float *)myalloc(c->gtex1max * sizeof(float));
    if (nmdaflag)
      c->gtex1_nmda = (float *)myalloc(c->gtex1max * sizeof(float));
    else
      c->gtex1_nmda = (float *)NULL;

    c->gtex1[0]  = c->gex1;  // Store initial value
    c->gtin1[0]  = c->gin1;  // Store initial value
    c->gtex1t[0] = c->gex1t;
    if (nmdaflag)
      c->gtex1_nmda[0] = c->gex1;  // Store initial value

    c->gtex1n = 1;

    // For currents, we also need Vm
    if ((strcmp(datid,"lgn_ia")==0)||
	(strcmp(datid,"ex_ix")==0)||
	(strcmp(datid,"in_ii")==0)){
      c->savv = 1;       // Save Vm values
      c->dxsav = 0.0004; // WYETH- Save values no more often than this (s)
    }

  }else if ((strcmp(datid,"gx")==0)||(strcmp(datid,"gi")==0)){
    ; // These are LGN layers, no prep needed.
  }else{

    csi = popc_csav_get_csi(c,datid);  // This name might be a 'csav'
    //
    //  WYETH - THIS won't get recorded unless we add to response list ??
    //

    if (csi < 0){
      sprintf(tstr,"datid = %s\n",datid);
      mylog(mylogf,tstr);
      mylog_exit(mylogf,"*** POP_CELL_SET_SAVE  Invalid data ID\n");
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            POP_CELL_UNSET_SAVE                            */
/*                                                                           */
/*****************************************************************************/
void pop_cell_unset_save(mylogf,c,datid)
     char *mylogf;
     struct pop_cell *c;
     char *datid;
{
  int csi;

  if (compare_prefix_string_order("all_spk_",datid)){
    exit_error("POP_CELL_UNSET_SAVE","This doesn't handle 'all_spk_'");
  }else if (strcmp(datid,"spikes")==0){
    c->savs = 0;
    if (c->s != NULL){
      myfree(c->s);
      c->s = NULL;
    }
  }else if ((strcmp(datid,"vm")==0)||(strcmp(datid,"rate")==0)){
    c->savv = 0;

    if (c->vm != NULL){
      myfree(c->vm); c->vm = NULL;
      myfree(c->vmt); c->vmt = NULL;
      if (c->savd){
	c->savd = 0;        // Save Vm Dendr values
	myfree(c->vmd); c->vmd = NULL;
      }
    }
  }else if (strcmp(datid,"vd")==0){
    c->savd = 0;
    if (c->vmd != NULL){
      myfree(c->vmd);
      c->vmd = NULL;
    }
  }else if ((strcmp(datid,"lgn_gx")==0)||  // AMPA + NMDA
	    (strcmp(datid,"lgn_gain")==0)||
	    (strcmp(datid,"lgn_ga")==0)||
	    (strcmp(datid,"lgn_ia")==0)||
	    (strcmp(datid,"lgn_gn")==0)||
	    (strcmp(datid,"gad")==0)||
	    (strcmp(datid,"ex_ix")==0)||
	    (strcmp(datid,"ex_gx")==0)||
	    (strcmp(datid,"in_ii")==0)||
	    (strcmp(datid,"in_gi")==0)){

    /**** WYEH - Big problem here - these are not independent, cannot
	  free individually, must know if cell has other savs in effect.  ***/

    c->savg = 0;

    if (c->gtex1 != NULL){
      myfree(c->gtex1);
      c->gtex1 = NULL;
    }

    if (c->gtex1_nmda != NULL){
      myfree(c->gtex1_nmda);
      c->gtex1_nmda = NULL;
    }
    if (c->gtin1 != NULL){
      myfree(c->gtin1);
      c->gtin1 = NULL;
    }
    if (c->gtex1t != NULL){
      myfree(c->gtex1t);
      c->gtex1t = NULL;
    }
  }else if ((strcmp(datid,"si_mask2")==0)||
	    (strcmp(datid,"si_mask1")==0)){

    printf("*** SI_MASK - not impl'd yet\n");
    mylog_exit(mylogf,"*** POP_CELL_UNSET_SAVE  Not implemented yet\n");
    
  }else{

    //
    //  Check for csav
    //
    csi = popc_csav_get_csi(c,datid);
    if (csi >= 0){
      ; // WYETH HERE - so what do we do?
    }else{
      printf("  Data ID:  %s\n",datid);
      mylog_exit(mylogf,"*** POP_CELL_UNSET_SAVE  Invalid data ID\n");
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          POP_CELL_UNSET_SAVE_SYN                          */
/*                                                                           */
/*****************************************************************************/
void pop_cell_unset_save_syn(mylogf,s,datid)
     char *mylogf;
     struct pop_syn *s;
     char *datid;
{
  int k;
  struct pop_si001 *si001;
  struct pop_si002 *si002;
  struct pop_si_sav *sisv;   // Saved data STRUCTURE to return for responses

  printf("  POP_CELL_UNSET_SAVE_SYN\n");
  printf("    pre-syn cell:  %s\n",s->pre->name);
  printf("    post-syn cell:  %s\n",s->post->name);
  printf("    data ID:  %s\n",datid);

  // WYETH VERY CRUDE - JUST DROPPING ALL STORAGE
  // WYETH VERY CRUDE - JUST DROPPING ALL STORAGE
  // WYETH VERY CRUDE - JUST DROPPING ALL STORAGE

  if (s->si->siui == 1){
    s->si->siu->s001->sisv = NULL;
  }else if (s->si->siui == 2){
    s->si->siu->s002->sisv = NULL;
  }
  return;

  //*****
  //*****
  //*****  WYETH - BIG HACK, DOESN"T TRY TO FIND THE DATA ELEMENT
  //*****
  //*****


  if (strcmp(datid,"si_mask1")==0){
    if (s->si->siui == 1){
      sisv = s->si->siu->s001->sisv;
    }else if (s->si->siui == 2){
      sisv = s->si->siu->s002->sisv;
    }

    k = search_2d_carray(sisv->name,datid,sisv->n);
    if (sisv->d[k] != NULL){
      myfree(sisv->d[k]);
      sisv->d[k] = NULL;
      sisv->cnt[k] = -1;
      if (sisv->name[k] != NULL){
	myfree(sisv->name[k]);
	sisv->name[k] = strdup("DELETED");
      }
    }
  }else if (strcmp(datid,"si_mask2")==0){
    if (s->si->siui == 1){
      exit_error("POP_CELL_UNSET_SAVE_SYN","should not happen, no mask2");
    }else if (s->si->siui == 2){
      sisv = s->si->siu->s002->sisv;
    }

    k = search_2d_carray(sisv->name,datid,sisv->n);
    if (sisv->d[k] != NULL){
      myfree(sisv->d[k]);
      sisv->d[k] = NULL;
      sisv->cnt[k] = -1;
      if (sisv->name[k] != NULL){
	myfree(sisv->name[k]);
	sisv->name[k] = strdup("DELETED");
      }
    }
  }else{
    printf("  Data ID:  %s\n",datid);
    mylog_exit(mylogf,"*** POP_CELL_UNSET_SAVE_SYN  Invalid data ID\n");
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                        POP_CELL_DERIVE_ORI_FROM_LAYER                     */
/*                                                                           */
/*  NOTES                                                                    */
/*  - If the synapse has a synaptic interaction (SI), only use if ci=0.      */
/*  -                                                                        */
/*                                                                           */
/*****************************************************************************/
void pop_cell_derive_ori_from_layer(mylogf,lay,c,dirflag)
     char *mylogf;
     struct pop_layer *lay;  // From this layer
     struct pop_cell *c;     // to this cell
     int dirflag;            // 0-assume 0..180,  1-assume 0..360
{
  int i;
  int n,flag,dir,ai_ori,ai_ori_cv;
  float *w,*ori,cr,ctheta,cv,ph_shift;
  char tstr[SLEN];
  struct pop_syn *s,*t;
  struct onode *sio;

  // Make arrays long enough to hold all relevant inputs
  ori = (float *)myalloc(c->nin*sizeof(float));
  w   = (float *)myalloc(c->nin*sizeof(float));

  // Scan pre-syn list for all synapses from 'lay'
  s = pop_cell_get_next_syn_layer_name(c->in,1,lay->name);
  n = 0;
  while(s != NULL){
    dir = 1;  // Assume default direction
    flag = 1;
    if (s->si != NULL){
      if (s->si->ci > 0)
	flag = 0;  // Not the primary player in this SI
      else if (s->si->ci == 0){
	if (s->si->siui == 1){
	  if (strcmp(s->si->siu->s001->msi->type,"si_ds02")==0){
	    sio = s->si->siu->s001->msi->o;  // Onode w/ description
	    ph_shift = onode_getpar_flt_exit(sio,"ph_shift");
	    if ((ph_shift > 0.0) && (ph_shift < 180.0)){
	      dir = -1;
	    }
	  }
	}
      }
    }

    if (flag){
      if (dir == 0){
	// WYETH ATTRIB
	ori[n] = pop_cell_attrib_get_f(s->pre,"ori");
	//ori[n] = s->pre->ori;
      }else{
	// WYETH ATTRIB
	ori[n] = pop_cell_attrib_get_f(s->pre,"ori") + 180.0;
	//ori[n] = s->pre->ori + 180.0;
      }
      w[n]   = s->w;
      n += 1;
    }

    t = s->post_next; // The post-syn cell of an input syn is this cell
    s = pop_cell_get_next_syn_layer_name(t,1,lay->name);
  }

  //for(i=0;i<n;i++)
  //printf("  w %f  ori %f\n",w[i],ori[i]);


  //  cr,ctheta  - centroid r,theta
  //  cv         - cicular variance
  if (dirflag == 1)
    centroid_circular_variance_polar_farray(w,ori,n,360.0,&cr,&ctheta,&cv);
  else
    centroid_circular_variance_polar_farray(w,ori,n,180.0,&cr,&ctheta,&cv);

  //printf("  Centroid r,theta   %f %f\n",cr,ctheta);
  //printf("  Circular variance  %f\n",cv);

  myfree(w); myfree(ori);

  ai_ori    = pop_cell_attrib_index_f(c->pl,"ori");
  ai_ori_cv = pop_cell_attrib_index_f(c->pl,"ori_cv");

  // WYETH ATTRIB
  c->attrib_f[ai_ori] = ctheta;
  c->attrib_f[ai_ori_cv] = cv;
  //c->ori = ctheta;
  //c->ori_cv = cv;
}
