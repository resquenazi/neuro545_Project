/*****************************************************************************/
/*                                                                           */
/*  ndata_util.c                                                             */
/*  Wyeth Bair                                                               */
/*                                                                           */
/*  This file contains routines to read and write the multichannel           */
/*  data format.                                                             */
/*                                                                           */
/*  Currently, "channels" is assumed to be constant over all trials.         */
/*                                                                           */
/*  Naming Conventions                                                       */
/*  - Routines accepting *nd should be named "NDATA_..."                     */
/*  - Routines accepting ndtrial_struct * should be named "NDATA_TRIAL_..."  */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h> // Added to make "atof" work in Linux
#include <math.h>
#include <string.h>

#include "my_util.h"
#include "misc_util.h"
#include "carray_util.h"
#include "iarray_util.h"
#include "farray_util.h"
#include "data_util.h"
#include "spike_util.h" // Had to add this for mean_spike_count_sarray
#include "ndata.h"

/**************************************-**************************************/
/*                                                                           */
/*                               NDATA_READ_NCHAR                            */
/*                                                                           */
/*  *** REPLACE w/ CARRAY_NCHAR_READ                                         */
/*                                                                           */
/*****************************************************************************/
void ndata_read_nchar(p,fin,revflag)
     char **p;
     FILE *fin;
     int revflag;
{
  int n;
  int nread;

  nread = revflag_fread((char *)&n,sizeof(int),1,fin,revflag);
  if (nread == 0)
    exit_error("NDATA_READ_NCHAR","Error");
  *p = (char *)myalloc((n+1)*sizeof(char));
  nread = fread((char *)(*p),sizeof(char),n,fin);
  (*p)[n] = '\0';
}
/**************************************-**************************************/
/*                                                                           */
/*                              NDATA_WRITE_NCHAR                            */
/*                                                                           */
/*  *** REPLACE w/ CARRAY_NCHAR_WRITE                                        */
/*                                                                           */
/*****************************************************************************/
void ndata_write_nchar(p,fout)
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
/*                               INIT_CONST_PARAM                            */
/*                                                                           */
/*****************************************************************************/
void init_const_param(nd,n)
     struct ndata_struct *nd;
     int n;
{
  nd->nconst = n;
  nd->cname = (char **)myalloc(n*sizeof(char *));
  nd->ctype = (char *)myalloc(n*sizeof(char));
  nd->cval = (char **)myalloc(n*sizeof(char *));
}
/**************************************-**************************************/
/*                                                                           */
/*                             INIT_SET_CONST_PARAM                          */
/*                                                                           */
/*****************************************************************************/
void init_set_const_param(nd,n,cname,ctype,cval)
     struct ndata_struct *nd;
     int n;
     char **cname,*ctype,**cval;
{
  int i;

  init_const_param(nd,n);
  for(i=0;i<n;i++){
    nd->cname[i] = strdup(cname[i]);
    nd->ctype[i] = ctype[i];
    nd->cval[i] = strdup(cval[i]);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              INIT_SET_VAR_PARAM                           */
/*                                                                           */
/*****************************************************************************/
void init_set_var_param(nd,n,nvar,vname,vtype,vval)
     struct ndata_struct *nd;
     int n,nvar;
     char **vname,*vtype,***vval; /* vval[ntrial][nvar] */
{
  int i,j;

  nd->nvar = nvar;
  nd->vname = (char **)myalloc(nvar*sizeof(char *));
  nd->vtype = (char *)myalloc(nvar*sizeof(char));

  for(i=0;i<nvar;i++){
    //printf("%s %c\n",vname[i],vtype[i]);
    nd->vname[i] = strdup(vname[i]);
    nd->vtype[i] = vtype[i];
  }

  for(i=0;i<n;i++){
    nd->t[i].nparam = nvar;
    nd->t[i].pname = (char **)myalloc(nvar*sizeof(char *));
    nd->t[i].pval = (char **)myalloc(nvar*sizeof(char *));
    for(j=0;j<nvar;j++){
      nd->t[i].pname[j] = strdup(vname[j]);
      nd->t[i].pval[j] = strdup(vval[i][j]);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                      NDATA_ADD_UNINITIALIZED_VAR_PARAM                    */
/*                                                                           */
/*****************************************************************************/
void ndata_add_uninitialized_var_param(nd,vname,vtype)
     struct ndata_struct *nd;
     char *vname,vtype;
{
  int i,j;
  int nvar,n;
  char **tname,*ttype,**tval;

  /*** Add the name for the varying parameter to the header. ***/
  nvar = nd->nvar + 1;
  tname = (char **)myalloc(nvar*sizeof(char *));
  ttype = (char *)myalloc(nvar*sizeof(char));
  for(i=0;i<(nvar-1);i++){
    tname[i] = nd->vname[i];
    ttype[i] = nd->vtype[i];
  }
  tname[nvar-1] = strdup(vname);
  ttype[nvar-1] = vtype;
  nd->nvar = nvar;
  myfree(nd->vname);
  nd->vname = tname;
  myfree(nd->vtype);
  nd->vtype = ttype;
  
  /*** Add a pointer for the value of the var param in each trial ***/
  n = nd->ntrial;
  for(i=0;i<n;i++){
    nd->t[i].nparam = nvar;

    tname = (char **)myalloc(nvar*sizeof(char *));
    tval = (char **)myalloc(nvar*sizeof(char *));
    for(j=0;j<(nvar-1);j++){
      tname[j] = nd->t[i].pname[j];
      tval[j] = nd->t[i].pval[j];
    }
    tname[nvar-1] = strdup(vname);
    tval[nvar-1] = NULL; /*** Add a null pointer. ***/
    myfree(nd->t[i].pname);
    nd->t[i].pname = tname;
    myfree(nd->t[i].pval);
    nd->t[i].pval = tval;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            NDATA_OPEN_SET_REVFLAG                         */
/*                                                                           */
/*****************************************************************************/
FILE *ndata_open_set_revflag(infile,rrevflag)
     char infile[];
     int *rrevflag;
{
  FILE *fopen(),*fin;
  int i,j,k;
  struct ndata_struct *nd;
  int temp,revflag,nread;

  printf("  NDATA_OPEN_SET_REVFLAG\n");

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** Cannot open file %s.\n",infile);
    return NULL;
  }
  printf("    Reading %s\n",infile);

  // Endian: If the first INT is not 1, then set "revflag"
  revflag = 0;
  nread = revflag_fread((char *)&temp,sizeof(int),1,fin,revflag);
  if (temp != 1)
    revflag = 1;

  *rrevflag = revflag;
  
  return fin;
}
/**************************************-**************************************/
/*                                                                           */
/*                       NDATA_READ_HEAD_AFTER_REVFLAG                       */
/*                                                                           */
/*  Assume the file is already open and the 'revflag' has been read.         */
/*  Read the rest of the header, stopping before the first trial.            */
/*                                                                           */
/*****************************************************************************/
struct ndata_struct *ndata_read_head_after_revflag(fin,revflag)
     FILE *fin;
     int revflag;
{
  int i,j,k;
  struct ndata_struct *nd;
  int temp,nread;

  //
  //  WYETH - SHOULD USE THIS BELOW TO AVOID REDUNDANT CODE ***** !!!!
  //  WYETH - SHOULD USE THIS BELOW TO AVOID REDUNDANT CODE ***** !!!!
  //  WYETH - SHOULD USE THIS BELOW TO AVOID REDUNDANT CODE ***** !!!!
  //

  printf("  NDATA_READ_HEAD_AFTER_REVFLAG\n");

  nd = (struct ndata_struct *)myalloc(sizeof(struct ndata_struct));

  ndata_read_nchar(&(nd->class),fin,revflag);
  // Read number, names, and values for constant parameters
  nread = revflag_fread((char *)&(nd->nconst),sizeof(int),1,fin,revflag);
  nd->cname = (char **)myalloc(nd->nconst*sizeof(char *));
  nd->ctype = (char *)myalloc(nd->nconst*sizeof(char));
  nd->cval = (char **)myalloc(nd->nconst*sizeof(char *));
  for(i=0;i<nd->nconst;i++){
    ndata_read_nchar(&(nd->cname[i]),fin,revflag);
    nread = fread((char *)&(nd->ctype[i]),sizeof(char),1,fin);
    ndata_read_nchar(&(nd->cval[i]),fin,revflag);
  }

  // Read number and names and types for variable parameters
  nread = revflag_fread(&(nd->nvar),sizeof(int),1,fin,revflag);
  nd->vname = (char **)myalloc(nd->nvar*sizeof(char *));
  nd->vtype = (char *)myalloc(nd->nvar*sizeof(char));
  for(i=0;i<nd->nvar;i++){
    ndata_read_nchar(&(nd->vname[i]),fin,revflag);
    nread = fread((char *)&(nd->vtype[i]),sizeof(char),1,fin);
  }

  // Read event code tables
  nread = revflag_fread(&(nd->ntable),sizeof(int),1,fin,revflag);
  nd->table = (struct ect_struct *)myalloc(nd->ntable*
					   sizeof(struct ect_struct));
  for(i=0;i<nd->ntable;i++){
    ndata_read_nchar(&(nd->table[i].tname),fin,revflag);
    nread = revflag_fread((char *)&(nd->table[i].n),sizeof(int),1,fin,revflag);
    nd->table[i].name = (char **)myalloc(nd->table[i].n*sizeof(char *));
    nd->table[i].tnum = i;
    nd->table[i].num = (int *)myalloc(nd->table[i].n*sizeof(int));
    for(j=0;j<nd->table[i].n;j++){
      ndata_read_nchar(&(nd->table[i].name[j]),fin,revflag);
      nread = revflag_fread((char *)&(nd->table[i].num[j]),sizeof(int),1,fin,
			    revflag);
    }
  }

  // Read number of trials
  nread = revflag_fread(&(nd->ntrial),sizeof(int),1,fin,revflag);

  nd->t = NULL;  // No trials read yet

  return nd;
}
/**************************************-**************************************/
/*                                                                           */
/*                          NDATA_READ_GET_NEXT_TRIAL                        */
/*                                                                           */
/*  Read the next trial from the open ndata file.                            */
/*                                                                           */
/*****************************************************************************/
struct ndtrial_struct *ndata_read_get_next_trial(fin,nd,revflag)
     FILE *fin;
     struct ndata_struct *nd;
     int revflag;
{
  int i,j,k;
  struct ndtrial_struct *t;
  int temp,nread;

  t = (struct ndtrial_struct *)myalloc(sizeof(struct ndtrial_struct));

  nread = revflag_fread(&(t->tcode),sizeof(int),1,fin,revflag);
  nread = revflag_fread(&(t->tref),sizeof(int),1,fin,revflag);
  nread = revflag_fread(&(t->nparam),sizeof(int),1,fin,revflag);
  t->pname = (char **)myalloc(t->nparam*sizeof(char *));
  t->pval = (char **)myalloc(t->nparam*sizeof(char *));
  for(j=0;j<t->nparam;j++){
    ndata_read_nchar(&(t->pname[j]),fin,revflag);
    ndata_read_nchar(&(t->pval[j]),fin,revflag);
  }
  // Read all records for this trial
  nread = revflag_fread(&(t->nrec),sizeof(int),1,fin,revflag);
  t->r = (struct ndrec_struct *)myalloc(t->nrec*sizeof(struct ndrec_struct));
  for(j=0;j<t->nrec;j++){
    nread = revflag_fread(&(t->r[j].rtype),sizeof(int),1,fin,revflag);
    ndata_read_nchar(&(t->r[j].name),fin,revflag);
    nread = revflag_fread(&(t->r[j].rcode),sizeof(int),1,fin,revflag);
    nread = revflag_fread(&(t->r[j].sampling),sizeof(float),1,fin,revflag);
    nread = revflag_fread(&(t->r[j].t0),sizeof(int),1,fin,revflag);
    nread = revflag_fread(&(t->r[j].tn),sizeof(int),1,fin,revflag);
    nread = revflag_fread(&(t->r[j].n),sizeof(int),1,fin,revflag);
    if (t->r[j].rtype == 0){
      t->r[j].p = (int *)myalloc(t->r[j].n*sizeof(int));
      nread = revflag_fread(t->r[j].p,sizeof(int),t->r[j].n,fin,revflag);
    }else if (t->r[j].rtype == 1){ // Use 'tn'
      t->r[j].x = (float *)myalloc(t->r[j].tn*sizeof(float));
      nread = revflag_fread(t->r[j].x,sizeof(float),t->r[j].tn,fin,revflag);
    }else if (t->r[j].rtype == 2){ // Use 'n'
      t->r[j].p = (int *)myalloc(t->r[j].n*sizeof(int));
      t->r[j].x = (float *)myalloc(t->r[j].n*sizeof(float));
      for(k=0;k<t->r[j].n;k++){
	nread = revflag_fread(&(t->r[j].p[k]),sizeof(int),1,fin,revflag);
	nread = revflag_fread(&(t->r[j].x[k]),sizeof(float),1,fin,revflag);
      }
    }else if (t->r[j].rtype == 3){ // Use 'n'
      nread = revflag_fread(&temp,sizeof(int),1,fin,revflag);
      if (temp >= 0)
	t->r[j].ect = &(nd->table[temp]);
      else
	t->r[j].ect = NULL;  // WYETH added to allow no EC Table, Jan 5, 2017
      t->r[j].p = (int *)myalloc(t->r[j].n*sizeof(int));
      t->r[j].ec = (int *)myalloc(t->r[j].n*sizeof(int));
      for(k=0;k<t->r[j].n;k++){
	nread = revflag_fread(&(t->r[j].p[k]),sizeof(int),1,fin,revflag);
	nread = revflag_fread(&(t->r[j].ec[k]),sizeof(int),1,fin,revflag);
      }
    }else
      exit_error("READ_NDATA","Unknown record type");
  }

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                           NDATA_GET_FILE_POS_NTR                          */
/*                                                                           */
/*  Find the position in the file where 'ntr' is stored.  This should be     */
/*  useful for updating the file.                                            */
/*                                                                           */
/*****************************************************************************/
long ndata_get_file_pos_ntr(infile)
     char infile[];
{
  FILE *fopen(),*fin;
  int i,j;
  int nc,nv,nt,n1,n2,temp,revflag,nread;
  //char s1[SLEN],c1;
  char *s1,c1;
  long plong;

  /*printf("  NDATA_GET_FILE_POS_NTR\n");*/

  //
  //  WYETH - can use 'ndata_open_set_revflag()' above ********************
  //  WYETH - can use 'ndata_open_set_revflag()' above ********************
  //  WYETH - can use 'ndata_open_set_revflag()' above ********************
  //

  s1 = (char *)myalloc(SLEN*sizeof(char)); // WYETH to satisfy compiler warning

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** Cannot open file %s.\n",infile);
    return (long)(-1);
  }
  /*printf("    Reading %s\n",infile);*/

  /*** We want this to work on Suns and DECs, which have endian conflict. ***/
  /*** If the first INT is not 1, then set "revflag". ***/
  revflag = 0;
  nread = revflag_fread((char *)&temp,sizeof(int),1,fin,revflag);
  if (temp != 1)
    revflag = 1;

  ndata_read_nchar(&s1,fin,revflag);
  /*** Read number, names, and values for constant parameters. ***/
  nread = revflag_fread((char *)&nc,sizeof(int),1,fin,revflag);
  for(i=0;i<nc;i++){
    ndata_read_nchar(&s1,fin,revflag);
    nread = fread((char *)&c1,sizeof(char),1,fin);
    ndata_read_nchar(&s1,fin,revflag);
  }

  /*** Read number and names and types for variable parameters. ***/
  nread = revflag_fread(&nv,sizeof(int),1,fin,revflag);
  for(i=0;i<nv;i++){
    ndata_read_nchar(&s1,fin,revflag);
    nread = fread((char *)&c1,sizeof(char),1,fin);
  }

  /*** Read event code tables ***/
  nread = revflag_fread(&nt,sizeof(int),1,fin,revflag);
  for(i=0;i<nt;i++){
    ndata_read_nchar(&s1,fin,revflag); /* Table name */
    nread = revflag_fread((char *)&n1,sizeof(int),1,fin,revflag);
    for(j=0;j<n1;j++){
      ndata_read_nchar(&s1,fin,revflag);
      nread = revflag_fread((char *)&n2,sizeof(int),1,fin,revflag);
    }
  }

  /*** GET POSITION OF NTRIALS ***/
  plong = ftell(fin);

  /*** Read number of trials. ***/
  nread = revflag_fread(&nt,sizeof(int),1,fin,revflag);

  /*printf("*******NTRIALS:  %d\n",nt);*/

  return plong;
}
/**************************************-**************************************/
/*                                                                           */
/*                            NDATA_OVERWRITE_NTR                            */
/*                                                                           */
/*  Go to position 'pos' and write the int 'ntr'.                            */
/*                                                                           */
/*****************************************************************************/
void ndata_overwrite_ntr(infile,pos,ntr)
     char infile[];
     long pos;
     int ntr;
{
  FILE *fopen(),*fin;
  int i,j;
  int temp,revflag,nread;

  /*printf("  NDATA_OVERWRITE_NTR\n");*/

  if ((fin = fopen(infile,"r+")) == NULL){
    printf("  *** Cannot open file %s.\n",infile);
    return;
  }

  /*** We want this to work on Suns and DECs, which have endian conflict. ***/
  /*** If the first INT is not 1, then set "revflag". ***/
  revflag = 0;
  nread = revflag_fread((char *)&temp,sizeof(int),1,fin,revflag);
  if (temp != 1)
    revflag = 1;


  /*printf("pos, ntr  %ld  %d\n",pos,ntr);*/
  
  fseek(fin,pos,0); /* 0 indicates position from *beginning* of file. */

  /*fwrite(&ntr,sizeof(int),1,fin);*/

  nread = revflag_fwrite(&ntr,sizeof(int),1,fin,revflag);
  if (nread != 1)
    exit_error("NDATA_OVERWRITE_NTR","Write failed");

  fclose(fin);
}
/**************************************-**************************************/
/*                                                                           */
/*                                 READ_NDATA                                */
/*                                                                           */
/*****************************************************************************/
void read_ndata(infile,rnd)
     char infile[];
     struct ndata_struct **rnd;
{
  FILE *fopen(),*fin;
  int i,j,k;
  struct ndata_struct *nd;
  int temp,revflag,nread;

  //
  //  WYETH - can use 'ndata_open_set_revflag()' above ********************
  //  WYETH - can use 'ndata_open_set_revflag()' above ********************
  //  WYETH - can use 'ndata_open_set_revflag()' above ********************
  //

  printf("  READ_NDATA\n");

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** Cannot open file %s.\n",infile);
    *rnd = NULL;
    return;
  }
  printf("    Reading %s\n",infile);
  nd = (struct ndata_struct *)myalloc(sizeof(struct ndata_struct));

  // We want this to work on Suns and DECs, which have endian conflict
  // If the first INT is not 1, then set "revflag"
  revflag = 0;
  nread = revflag_fread((char *)&temp,sizeof(int),1,fin,revflag);
  if (temp != 1)
    revflag = 1;

  ndata_read_nchar(&(nd->class),fin,revflag);
  // Read number, names, and values for constant parameters
  nread = revflag_fread((char *)&(nd->nconst),sizeof(int),1,fin,revflag);
  nd->cname = (char **)myalloc(nd->nconst*sizeof(char *));
  nd->ctype = (char *)myalloc(nd->nconst*sizeof(char));
  nd->cval = (char **)myalloc(nd->nconst*sizeof(char *));
  for(i=0;i<nd->nconst;i++){
    ndata_read_nchar(&(nd->cname[i]),fin,revflag);
    nread = fread((char *)&(nd->ctype[i]),sizeof(char),1,fin);
    ndata_read_nchar(&(nd->cval[i]),fin,revflag);
  }

  // Read number and names and types for variable parameters
  nread = revflag_fread(&(nd->nvar),sizeof(int),1,fin,revflag);
  nd->vname = (char **)myalloc(nd->nvar*sizeof(char *));
  nd->vtype = (char *)myalloc(nd->nvar*sizeof(char));
  for(i=0;i<nd->nvar;i++){
    ndata_read_nchar(&(nd->vname[i]),fin,revflag);
    nread = fread((char *)&(nd->vtype[i]),sizeof(char),1,fin);
  }

  // Read event code tables
  nread = revflag_fread(&(nd->ntable),sizeof(int),1,fin,revflag);
  nd->table = (struct ect_struct *)myalloc(nd->ntable*
					   sizeof(struct ect_struct));
  for(i=0;i<nd->ntable;i++){
    ndata_read_nchar(&(nd->table[i].tname),fin,revflag);
    nread = revflag_fread((char *)&(nd->table[i].n),sizeof(int),1,fin,revflag);
    nd->table[i].name = (char **)myalloc(nd->table[i].n*sizeof(char *));
    nd->table[i].tnum = i;
    nd->table[i].num = (int *)myalloc(nd->table[i].n*sizeof(int));
    for(j=0;j<nd->table[i].n;j++){
      ndata_read_nchar(&(nd->table[i].name[j]),fin,revflag);
      nread = revflag_fread((char *)&(nd->table[i].num[j]),sizeof(int),1,fin,
			    revflag);
    }
  }

  // Read number of trials
  nread = revflag_fread(&(nd->ntrial),sizeof(int),1,fin,revflag);
  nd->t = (struct ndtrial_struct *)myalloc(nd->ntrial*
					 sizeof(struct ndtrial_struct));
  for(i=0;i<nd->ntrial;i++){
    nread = revflag_fread(&(nd->t[i].tcode),sizeof(int),1,fin,revflag);
    nread = revflag_fread(&(nd->t[i].tref),sizeof(int),1,fin,revflag);
    nread = revflag_fread(&(nd->t[i].nparam),sizeof(int),1,fin,revflag);
    nd->t[i].pname = (char **)myalloc(nd->t[i].nparam*sizeof(char *));
    nd->t[i].pval = (char **)myalloc(nd->t[i].nparam*sizeof(char *));
    for(j=0;j<nd->t[i].nparam;j++){
      ndata_read_nchar(&(nd->t[i].pname[j]),fin,revflag);
      ndata_read_nchar(&(nd->t[i].pval[j]),fin,revflag);
    }
    // Read all records for this trial
    nread = revflag_fread(&(nd->t[i].nrec),sizeof(int),1,fin,revflag);
    nd->t[i].r = (struct ndrec_struct *)myalloc(nd->t[i].nrec*
					      sizeof(struct ndrec_struct));
    for(j=0;j<nd->t[i].nrec;j++){
      nread = revflag_fread(&(nd->t[i].r[j].rtype),sizeof(int),1,fin,revflag);
      ndata_read_nchar(&(nd->t[i].r[j].name),fin,revflag);
      nread = revflag_fread(&(nd->t[i].r[j].rcode),sizeof(int),1,fin,revflag);
      nread = revflag_fread(&(nd->t[i].r[j].sampling),sizeof(float),1,fin,
			    revflag);
      nread = revflag_fread(&(nd->t[i].r[j].t0),sizeof(int),1,fin,revflag);
      nread = revflag_fread(&(nd->t[i].r[j].tn),sizeof(int),1,fin,revflag);
      nread = revflag_fread(&(nd->t[i].r[j].n),sizeof(int),1,fin,revflag);
      if (nd->t[i].r[j].rtype == 0){
	nd->t[i].r[j].p = (int *)myalloc(nd->t[i].r[j].n*sizeof(int));
	nread = revflag_fread(nd->t[i].r[j].p,sizeof(int),nd->t[i].r[j].n,
			      fin,revflag);
      }else if (nd->t[i].r[j].rtype == 1){ // USE "tn"
	nd->t[i].r[j].x = (float *)myalloc(nd->t[i].r[j].tn*sizeof(float));
	nread = revflag_fread(nd->t[i].r[j].x,sizeof(float),nd->t[i].r[j].tn,
			      fin,revflag);
      }else if (nd->t[i].r[j].rtype == 2){ // USE "n"
	nd->t[i].r[j].p = (int *)myalloc(nd->t[i].r[j].n*sizeof(int));
	nd->t[i].r[j].x = (float *)myalloc(nd->t[i].r[j].n*sizeof(float));
	for(k=0;k<nd->t[i].r[j].n;k++){
	  nread = revflag_fread(&(nd->t[i].r[j].p[k]),sizeof(int),1,
				fin,revflag);
	  nread = revflag_fread(&(nd->t[i].r[j].x[k]),sizeof(float),1,
				fin,revflag);
	}
      }else if (nd->t[i].r[j].rtype == 3){ // USE "n"
	nread = revflag_fread(&temp,sizeof(int),1,fin,revflag);
	if (temp >= 0) // WYETH - added this to allow no ECTable, Jan 5, 2017
	  nd->t[i].r[j].ect = &(nd->table[temp]);
	else
	  nd->t[i].r[j].ect = NULL;
	nd->t[i].r[j].p = (int *)myalloc(nd->t[i].r[j].n*sizeof(int));
	nd->t[i].r[j].ec = (int *)myalloc(nd->t[i].r[j].n*sizeof(int));
	for(k=0;k<nd->t[i].r[j].n;k++){
	  nread = revflag_fread(&(nd->t[i].r[j].p[k]),sizeof(int),1,
				fin,revflag);
	  nread = revflag_fread(&(nd->t[i].r[j].ec[k]),sizeof(int),1,
				fin,revflag);
	}
      }else
	exit_error("READ_NDATA","Unknown record type");
    }
  }
  fclose(fin);
  printf("    %d trials found.\n",nd->ntrial);
  *rnd = nd;
  i = nread; // Keep lint from complaining
}
/**************************************-**************************************/
/*                                                                           */
/*                              READ_NDATA_TEXT                              */
/*                                                                           */
/*****************************************************************************/
void read_ndata_text(infile,rnd)
     char infile[];
     struct ndata_struct **rnd;
{
  FILE *fopen(),*fin;
  int i,j,k;
  struct ndata_struct *nd;
  int ni,temp,ns;
  char t[SLEN],tc[SLEN],ts[SLEN],**ilist,*tg;

  printf("  READ_NDATA_TEXT\n");

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** Cannot open file %s.\n",infile);
    *rnd = NULL;
    return;
  }
  printf("    Reading %s\n",infile);
  nd = (struct ndata_struct *)myalloc(sizeof(struct ndata_struct));

  ns = fscanf(fin,"%s %s",t,ts); check_string_exit_carray(t,"CLASS");
  nd->class = strdup(ts);
  /*** Read number, names, and values for constant parameters. ***/
  ns = fscanf(fin,"%s %d",t,&nd->nconst); check_string_exit_carray(t,"NCONST");
  nd->cname = (char **)myalloc(nd->nconst*sizeof(char *));
  nd->ctype = (char *)myalloc(nd->nconst*sizeof(char));
  nd->cval = (char **)myalloc(nd->nconst*sizeof(char *));
  for(i=0;i<nd->nconst;i++){
    ns = fscanf(fin,"%s %s",tc,t);
    nd->ctype[i] = tc[0];
    nd->cname[i] = strdup(t);
    /*** WYETH - below is a back hack to read the rest of the line ***/
    tg = fgets(ts,SLEN,fin); // Read until the carriage return.
    ts[strlen(ts)-1] = '\0'; // Remove final CR
    get_items_from_string(ts,&ilist,&ni);
    nd->cval[i] = make_string_from_items(ilist,ni);
    free_2d_carray(ilist,ni);
  }

  /*** Read number and names and types for variable parameters. ***/
  ns = fscanf(fin,"%s %d",t,&nd->nvar); check_string_exit_carray(t,"NVAR");
  nd->vname = (char **)myalloc(nd->nvar*sizeof(char *));
  nd->vtype = (char *)myalloc(nd->nvar*sizeof(char));
  for(i=0;i<nd->nvar;i++){
    ns = fscanf(fin,"%s %s",tc,t);
    nd->vtype[i] = tc[0];
    nd->vname[i] = strdup(t);
  }

  /*** Read event code tables ***/
  ns = fscanf(fin,"%s %d",t,&nd->ntable); check_string_exit_carray(t,"NTABLE");
  nd->table = (struct ect_struct *)myalloc(nd->ntable*
					   sizeof(struct ect_struct));
  for(i=0;i<nd->ntable;i++){
    ns = fscanf(fin,"%s %d",t,&nd->table[i].n);
    nd->table[i].tname = strdup(t);
    nd->table[i].name = (char **)myalloc(nd->table[i].n*sizeof(char *));
    nd->table[i].num = (int *)myalloc(nd->table[i].n*sizeof(int));
    for(j=0;j<nd->table[i].n;j++){
      ns = fscanf(fin,"%s %d",t,&nd->table[i].num[j]);
      nd->table[i].name[j] = strdup(t);
    }
  }

  /*** Read number of trials. ***/
  ns = fscanf(fin,"%s %d",t,&nd->ntrial); check_string_exit_carray(t,"NTRIALS"); 
  nd->t = (struct ndtrial_struct *)myalloc(nd->ntrial*
					 sizeof(struct ndtrial_struct));
  for(i=0;i<nd->ntrial;i++){
    ns = fscanf(fin,"%s %d",t,&nd->t[i].tcode);
    check_string_exit_carray(t,"TCODE");
    ns = fscanf(fin,"%s %d",t,&nd->t[i].tref);
    check_string_exit_carray(t,"TREF");
    ns = fscanf(fin,"%s %d",t,&nd->t[i].nparam);
    check_string_exit_carray(t,"NPARAM");
    nd->t[i].pname = (char **)myalloc(nd->t[i].nparam*sizeof(char *));
    nd->t[i].pval = (char **)myalloc(nd->t[i].nparam*sizeof(char *));
    for(j=0;j<nd->t[i].nparam;j++){
      ns = fscanf(fin,"%s %s",t,ts);
      nd->t[i].pname[j] = strdup(t);
      nd->t[i].pval[j] = strdup(ts);
    }
    /*** Read all records for this trial. ***/
    ns = fscanf(fin,"%s %d",t,&nd->t[i].nrec);
    check_string_exit_carray(t,"NREC");
    nd->t[i].r = (struct ndrec_struct *)myalloc(nd->t[i].nrec*
						sizeof(struct ndrec_struct));
    for(j=0;j<nd->t[i].nrec;j++){
      ns = fscanf(fin,"%s %d",t,&nd->t[i].r[j].rtype);
      check_string_exit_carray(t,"RTYPE");
      ns = fscanf(fin,"%s %s",t,ts); check_string_exit_carray(t,"NAME");
      nd->t[i].r[j].name = strdup(ts);
      ns = fscanf(fin,"%s %d",t,&nd->t[i].r[j].rcode);
      check_string_exit_carray(t,"RCODE");
      ns = fscanf(fin,"%s %f",t,&nd->t[i].r[j].sampling);
      check_string_exit_carray(t,"SAMPLING");
      ns = fscanf(fin,"%s %d",t,&nd->t[i].r[j].t0);
      check_string_exit_carray(t,"t0");
      ns = fscanf(fin,"%s %d",t,&nd->t[i].r[j].tn);
      check_string_exit_carray(t,"tN");
      ns = fscanf(fin,"%s %d",t,&nd->t[i].r[j].n);
      check_string_exit_carray(t,"N");
      if (nd->t[i].r[j].rtype == 0){
	nd->t[i].r[j].p = (int *)myalloc(nd->t[i].r[j].n*sizeof(int));
	for(k=0;k<nd->t[i].r[j].n;k++)
	  ns = fscanf(fin,"%d",&nd->t[i].r[j].p[k]);
      }else if (nd->t[i].r[j].rtype == 1){ /*** USE "tn" ***/
	nd->t[i].r[j].x = (float *)myalloc(nd->t[i].r[j].tn*sizeof(float));
	for(k=0;k<nd->t[i].r[j].n;k++)
	  ns = fscanf(fin,"%f",&nd->t[i].r[j].x[k]);
      }else if (nd->t[i].r[j].rtype == 2){ /*** USE "n" ***/
	nd->t[i].r[j].p = (int *)myalloc(nd->t[i].r[j].n*sizeof(int));
	nd->t[i].r[j].x = (float *)myalloc(nd->t[i].r[j].n*sizeof(float));
	for(k=0;k<nd->t[i].r[j].n;k++)
	  ns = fscanf(fin,"%d %f",&nd->t[i].r[j].p[k],&nd->t[i].r[j].x[k]);
      }else if (nd->t[i].r[j].rtype == 3){ /*** USE "n" ***/
	ns = fscanf(fin,"%s %d",t,&temp);
	check_string_exit_carray(t,"TABLENUM");
	if (temp >= 0)
	  nd->t[i].r[j].ect = &(nd->table[temp]);
	else
	  nd->t[i].r[j].ect = NULL;
	nd->t[i].r[j].p = (int *)myalloc(nd->t[i].r[j].n*sizeof(int));
	nd->t[i].r[j].ec = (int *)myalloc(nd->t[i].r[j].n*sizeof(int));
	for(k=0;k<nd->t[i].r[j].n;k++)
	  ns = fscanf(fin,"%d %d",&nd->t[i].r[j].p[k],&nd->t[i].r[j].ec[k]);
      }else
	exit_error("READ_NDATA","Unknown record type");
    }
  }
  fclose(fin);
  printf("    %d trials found.\n",nd->ntrial);
  *rnd = nd;
}
/**************************************-**************************************/
/*                                                                           */
/*                              NDATA_WRITE_TRIAL                            */
/*                                                                           */
/*****************************************************************************/
void ndata_write_trial(fout,t,debugflag)
     FILE *fout;
     struct ndtrial_struct *t;
     int debugflag;
{
  int j,k;
  int ectnum;

  if (debugflag) printf("  NDATA_WRITE_TRIAL\n");

  fwrite(&(t->tcode),sizeof(int),1,fout);
  fwrite(&(t->tref),sizeof(int),1,fout);
  fwrite(&(t->nparam),sizeof(int),1,fout);
  for(j=0;j<t->nparam;j++){
    ndata_write_nchar(t->pname[j],fout);
    ndata_write_nchar(t->pval[j],fout);
  }
  if (debugflag) printf("  lower\n");
  
  /*** Write all records for this trial. ***/
  fwrite(&(t->nrec),sizeof(int),1,fout);
  if (debugflag) printf("  wrote nrec (%d)\n",t->nrec);
  for(j=0;j<t->nrec;j++){
    if (debugflag) printf("  j=%d\n",j);
    fwrite(&(t->r[j].rtype),sizeof(int),1,fout);
    if (debugflag) printf("  wrote rtype (%d)\n",t->r[j].rtype);
    ndata_write_nchar(t->r[j].name,fout);
    if (debugflag) printf("  wrote name (%s)\n",t->r[j].name);
    fwrite(&(t->r[j].rcode),sizeof(int),1,fout);
    if (debugflag) printf("  wrote rcode (%d)\n",t->r[j].rcode);
    fwrite(&(t->r[j].sampling),sizeof(float),1,fout);
    if (debugflag) printf("  wrote sampling (%f)\n",t->r[j].sampling);
    fwrite(&(t->r[j].t0),sizeof(int),1,fout);
    if (debugflag) printf("  wrote t0 (%d)\n",t->r[j].t0);
    fwrite(&(t->r[j].tn),sizeof(int),1,fout);
    if (debugflag) printf("  wrote tn (%d)\n",t->r[j].tn);
    fwrite(&(t->r[j].n),sizeof(int),1,fout);
    if (debugflag) printf("  wrote n\n");
    
    if (t->r[j].rtype == 0){
      if (debugflag) printf("    rtype 0\n");
      if (debugflag) printf("      n=%d\n",t->r[j].n);
      fwrite(t->r[j].p,sizeof(int),t->r[j].n,fout);
    }else if (t->r[j].rtype == 1){ /*** USE "tn" ***/
      if (debugflag) printf("    rtype 1\n");
      if (debugflag) printf("      tn = %d\n",t->r[j].tn);
      fwrite(t->r[j].x,sizeof(float),t->r[j].tn,fout);
      if (debugflag) printf("      wrote 1\n");
    }else if (t->r[j].rtype == 2){ /*** USE "n" ***/
      if (debugflag) printf("    rtype 2\n");
      for(k=0;k<t->r[j].n;k++){
	fwrite(&(t->r[j].p[k]),sizeof(int),1,fout);
	fwrite(&(t->r[j].x[k]),sizeof(float),1,fout);
      }
    }else if (t->r[j].rtype == 3){ /*** USE "n" ***/
      if (debugflag) printf("    rtype 3\n");

      // WYETH - added condition to allow NULL EC table, Jan 5, 2017
      if (t->r[j].ect != NULL)
	ectnum = t->r[j].ect->tnum;
      else
	ectnum = -1;  // Indicate that there is no table
      fwrite(&ectnum,sizeof(int),1,fout);
      //fwrite(&(t->r[j].ect->tnum),sizeof(int),1,fout);


      for(k=0;k<t->r[j].n;k++){
	fwrite(&(t->r[j].p[k]),sizeof(int),1,fout);
	fwrite(&(t->r[j].ec[k]),sizeof(int),1,fout);
      }
    }else
      exit_error("NDATA_WRITE_TRIAL","Unknown record type");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                WRITE_NDATA                                */
/*                                                                           */
/*****************************************************************************/
void write_ndata(outfile,nd)
     char outfile[];
     struct ndata_struct *nd;
{
  FILE *fopen(),*fout;
  int i,j,k;
  int temp,debugflag;

  printf("  WRITE_NDATA\n");

  debugflag = 0;

  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** File %s\n",outfile);
    exit_error("WRITE_NDATA","Cannot open file");
  }
  printf("    Writing %d trials to %s\n",nd->ntrial,outfile);

  temp = 1;
  fwrite((char *)&temp,sizeof(int),1,fout);

  if (debugflag) printf("WRITE_NDATA debug ***  00\n");

  ndata_write_nchar(nd->class,fout);
  /*** Write number, names, types, and values for constant parameters. ***/
  fwrite(&(nd->nconst),sizeof(int),1,fout);
  for(i=0;i<nd->nconst;i++){
    ndata_write_nchar(nd->cname[i],fout);
    fwrite(&(nd->ctype[i]),sizeof(char),1,fout);
    ndata_write_nchar(nd->cval[i],fout);
  }
  if (debugflag) printf("WRITE_NDATA debug ***  10\n");

  /*** Write number and names and types for variable parameters. ***/
  fwrite(&(nd->nvar),sizeof(int),1,fout);
  for(i=0;i<nd->nvar;i++){
    ndata_write_nchar(nd->vname[i],fout);
    fwrite(&(nd->vtype[i]),sizeof(char),1,fout);
  }
  if (debugflag) printf("WRITE_NDATA debug ***  20\n");

  /*** Write event code tables ***/
  fwrite(&(nd->ntable),sizeof(int),1,fout);
  if (debugflag) printf("  wrote number of tables:  %d\n",nd->ntable);
  for(i=0;i<nd->ntable;i++){
    ndata_write_nchar(nd->table[i].tname,fout);
    fwrite(&(nd->table[i].n),sizeof(int),1,fout);
    for(j=0;j<nd->table[i].n;j++){
      ndata_write_nchar(nd->table[i].name[j],fout);
      fwrite(&(nd->table[i].num[j]),sizeof(int),1,fout);
      if (debugflag) printf("  wrote table entry:  %3d  %s\n",
			    nd->table[i].num[j],nd->table[i].name[j]);
    }
  }
  if (debugflag) printf("WRITE_NDATA debug ***  30\n");

  // Write number of trials, followed by trial data
  fwrite(&(nd->ntrial),sizeof(int),1,fout);
  for(i=0;i<nd->ntrial;i++){
    if (debugflag)
      printf("i=%d\n",i);

    ndata_write_trial(fout,&(nd->t[i]),debugflag);
  }
  if (debugflag) printf("WRITE_NDATA debug ***  40\n");

  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                             NDATA_APPEND_TRIAL                            */
/*                                                                           */
/*****************************************************************************/
void ndata_append_trial(outfile,t)
     char outfile[];
     struct ndtrial_struct *t;
{
  FILE *fopen(),*fout;
  int debugflag;

  /*printf("  NDATA_APPEND_TRIAL\n");*/

  debugflag = 0;

  if ((fout = fopen(outfile,"a")) == NULL){
    printf("  *** File %s\n",outfile);
    exit_error("NDATA_APPEND_TRIAL","Cannot open file");
  }
  /*printf("    Appending one trial to %s\n",outfile);*/

  ndata_write_trial(fout,t,debugflag);
  
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                              WRITE_NDATA_LINK                             */
/*                                                                           */
/*  Use the links between trials.                                            */
/*                                                                           */
/*****************************************************************************/
void write_ndata_link(outfile,nd)
     char outfile[];
     struct ndata_struct *nd;
{
  FILE *fopen(),*fout;
  int i,j,k;
  int temp,debugflag,ectnum;
  struct ndtrial_struct *tp;

  printf("  WRITE_NDATA_LINK\n");

  debugflag = 0;

  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** File %s\n",outfile);
    exit_error("WRITE_NDATA_LINK","Cannot open file");
  }
  printf("    Writing %d trials to %s\n",nd->ntrial,outfile);

  temp = 1;
  fwrite((char *)&temp,sizeof(int),1,fout);

  if (debugflag) printf("WRITE_NDATA_LINK debug ***  00\n");

  ndata_write_nchar(nd->class,fout);
  /*** Write number, names, types, and values for constant parameters. ***/
  fwrite(&(nd->nconst),sizeof(int),1,fout);
  for(i=0;i<nd->nconst;i++){
    ndata_write_nchar(nd->cname[i],fout);
    fwrite(&(nd->ctype[i]),sizeof(char),1,fout);
    ndata_write_nchar(nd->cval[i],fout);
  }
  if (debugflag) printf("WRITE_NDATA_LINK debug ***  10\n");

  /*** Write number and names and types for variable parameters. ***/
  fwrite(&(nd->nvar),sizeof(int),1,fout);
  for(i=0;i<nd->nvar;i++){
    ndata_write_nchar(nd->vname[i],fout);
    fwrite(&(nd->vtype[i]),sizeof(char),1,fout);
  }
  if (debugflag) printf("WRITE_NDATA_LINK debug ***  20\n");

  /*** Write event code tables ***/
  fwrite(&(nd->ntable),sizeof(int),1,fout);
  for(i=0;i<nd->ntable;i++){
    ndata_write_nchar(nd->table[i].tname,fout);
    fwrite(&(nd->table[i].n),sizeof(int),1,fout);
    for(j=0;j<nd->table[i].n;j++){
      ndata_write_nchar(nd->table[i].name[j],fout);
      fwrite(&(nd->table[i].num[j]),sizeof(int),1,fout);
    }
  }
  if (debugflag) printf("WRITE_NDATA_LINK debug ***  30\n");

  // Write number of trials
  fwrite(&(nd->ntrial),sizeof(int),1,fout);
  tp = nd->t;
  for(i=0;i<nd->ntrial;i++){
    if (debugflag) printf("i=%d\n",i);
    fwrite(&(tp->tcode),sizeof(int),1,fout);
    fwrite(&(tp->tref),sizeof(int),1,fout);
    fwrite(&(tp->nparam),sizeof(int),1,fout);
    for(j=0;j<tp->nparam;j++){
      ndata_write_nchar(tp->pname[j],fout);
      ndata_write_nchar(tp->pval[j],fout);
    }
    if (debugflag) printf("  lower\n");

    /*** Write all records for this trial. ***/
    fwrite(&(tp->nrec),sizeof(int),1,fout);
    if (debugflag) printf("  wrote nrec (%d)\n",tp->nrec);
    for(j=0;j<tp->nrec;j++){
      if (debugflag) printf("  j=%d\n",j);
      fwrite(&(tp->r[j].rtype),sizeof(int),1,fout);
      if (debugflag) printf("  wrote rtype (%d)\n",tp->r[j].rtype);
      ndata_write_nchar(tp->r[j].name,fout);
      if (debugflag) printf("  wrote name (%s)\n",tp->r[j].name);
      fwrite(&(tp->r[j].rcode),sizeof(int),1,fout);
      if (debugflag) printf("  wrote rcode (%d)\n",tp->r[j].rcode);
      fwrite(&(tp->r[j].sampling),sizeof(float),1,fout);
      if (debugflag) printf("  wrote sampling (%f)\n",tp->r[j].sampling);
      fwrite(&(tp->r[j].t0),sizeof(int),1,fout);
      if (debugflag) printf("  wrote t0 (%d)\n",tp->r[j].t0);
      fwrite(&(tp->r[j].tn),sizeof(int),1,fout);
      if (debugflag) printf("  wrote tn (%d)\n",tp->r[j].tn);
      fwrite(&(tp->r[j].n),sizeof(int),1,fout);
      if (debugflag) printf("  wrote n\n");

      if (tp->r[j].rtype == 0){
	if (debugflag) printf("    rtype 0\n");
	if (debugflag) printf("      n=%d\n",tp->r[j].n);
	fwrite(tp->r[j].p,sizeof(int),tp->r[j].n,fout);
      }else if (tp->r[j].rtype == 1){ /*** USE "tn" ***/
	if (debugflag) printf("    rtype 1\n");
	if (debugflag) printf("      tn = %d\n",tp->r[j].tn);
	fwrite(tp->r[j].x,sizeof(float),tp->r[j].tn,fout);
	if (debugflag) printf("      wrote 1\n");
      }else if (tp->r[j].rtype == 2){ /*** USE "n" ***/
	if (debugflag) printf("    rtype 2\n");
	for(k=0;k<tp->r[j].n;k++){
	  fwrite(&(tp->r[j].p[k]),sizeof(int),1,fout);
	  fwrite(&(tp->r[j].x[k]),sizeof(float),1,fout);
	}
      }else if (tp->r[j].rtype == 3){ /*** USE "n" ***/
	if (debugflag) printf("    rtype 3\n");

	// WYETH - added condition to allow NULL EC table, Jan 5, 2017
	if (tp->r[j].ect != NULL)
	  ectnum = tp->r[j].ect->tnum;
	else
	  ectnum = -1;  // Indicate that there is no table
	fwrite(&ectnum,sizeof(int),1,fout);

	if (debugflag) printf("      wrote table num\n");
	for(k=0;k<tp->r[j].n;k++){
	  fwrite(&(tp->r[j].p[k]),sizeof(int),1,fout);
	  fwrite(&(tp->r[j].ec[k]),sizeof(int),1,fout);
	}
	if (debugflag) printf("      wrote event list\n");
      }else
	exit_error("WRITE_NDATA_LINK","Unknown record type");
    }
    tp = tp->next;
  }
  if (debugflag) printf("WRITE_NDATA_LINK debug ***  40\n");

  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                              WRITE_NDATA_TEXT                             */
/*                                                                           */
/*****************************************************************************/
void write_ndata_text(outfile,nd)
     char outfile[];
     struct ndata_struct *nd;
{
  FILE *fopen(),*fout;
  int i,j,k;
  int ectnum;

  printf("  WRITE_NDATA_TEXT\n");

  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** File %s\n",outfile);
    exit_error("WRITE_NDATA_TEXT","Cannot open file");
  }

  fprintf(fout,"CLASS %s\n",nd->class);
  /*** Write number, names, types, and values for constant parameters. ***/
  fprintf(fout,"NCONST %d\n",nd->nconst);
  for(i=0;i<nd->nconst;i++)
    fprintf(fout,"  %c %s %s\n",nd->ctype[i],nd->cname[i],nd->cval[i]);
  /*** Write number and names and types for variable parameters. ***/
  fprintf(fout,"NVAR %d\n",nd->nvar);
  for(i=0;i<nd->nvar;i++)
    fprintf(fout,"  %c %s\n",nd->vtype[i],nd->vname[i]);
  /*** Write event code tables ***/
  fprintf(fout,"NTABLE %d\n",nd->ntable);
  for(i=0;i<nd->ntable;i++){
    fprintf(fout,"%s %d\n",nd->table[i].tname,nd->table[i].n);
    for(j=0;j<nd->table[i].n;j++)
      fprintf(fout,"  %20s %6d\n",nd->table[i].name[j],nd->table[i].num[j]);
  }
  /*** Write number of trials. ***/
  fprintf(fout,"NTRIALS %d\n",nd->ntrial);
  for(i=0;i<nd->ntrial;i++){
    fprintf(fout,"TCODE %d  ",nd->t[i].tcode);
    fprintf(fout,"TREF %d  ",nd->t[i].tref);
    fprintf(fout,"NPARAM %d\n",nd->t[i].nparam);
    for(j=0;j<nd->t[i].nparam;j++)
      fprintf(fout,"  %s %s\n",nd->t[i].pname[j],nd->t[i].pval[j]);
    /*** Write all records for this trial. ***/
    fprintf(fout,"NREC %d\n",nd->t[i].nrec);
    for(j=0;j<nd->t[i].nrec;j++){
      fprintf(fout,"RTYPE %d  ",nd->t[i].r[j].rtype);
      fprintf(fout,"NAME %s  ",nd->t[i].r[j].name);
      fprintf(fout,"RCODE %d  ",nd->t[i].r[j].rcode);
      fprintf(fout,"SAMPLING %f  ",nd->t[i].r[j].sampling);
      fprintf(fout,"t0 %d  ",nd->t[i].r[j].t0);
      fprintf(fout,"tN %d  ",nd->t[i].r[j].tn);
      fprintf(fout,"N %d\n",nd->t[i].r[j].n);
      if (nd->t[i].r[j].rtype == 0){
	print_iarray_80_column(fout,nd->t[i].r[j].p,nd->t[i].r[j].n);
      }else if (nd->t[i].r[j].rtype == 1){ /*** USE "tn" ***/
	for(k=0;k<nd->t[i].r[j].tn;k++)
	  fprintf(fout,"%f ",nd->t[i].r[j].x[k]);
	fprintf(fout,"\n");
      }else if (nd->t[i].r[j].rtype == 2){ /*** USE "n" ***/
	for(k=0;k<nd->t[i].r[j].n;k++)
	  fprintf(fout,"%d %f ",nd->t[i].r[j].p[k],nd->t[i].r[j].x[k]);
	fprintf(fout,"\n");
      }else if (nd->t[i].r[j].rtype == 3){ /*** USE "n" ***/

	// WYETH - added condition to allow NULL EC table, Jan 5, 2017
	if (nd->t[i].r[j].ect != NULL)
	  ectnum = nd->t[i].r[j].ect->tnum;
	else
	  ectnum = -1;  // Indicate that there is no table
	fprintf(fout,"TABLENUM %d\n",ectnum);

	for(k=0;k<nd->t[i].r[j].n;k++)
	  fprintf(fout,"%d %d ",nd->t[i].r[j].p[k],nd->t[i].r[j].ec[k]);
	fprintf(fout,"\n");
      }else
	exit_error("WRITE_NDATA_TEXT","Unknown record type");
    }
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                               WRITE_NDATA_TXT                             */
/*                                                                           */
/*  Write a simplified text format file.                                     */
/*                                                                           */
/*****************************************************************************/
void write_ndata_txt(outfile,nd)
     char outfile[];
     struct ndata_struct *nd;
{
  FILE *fopen(),*fout;
  int i,j,k;
  int ectnum;

  printf("  WRITE_NDATA_TXT\n");

  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** File %s\n",outfile);
    exit_error("WRITE_NDATA_TXT","Cannot open file");
  }

  /*fprintf(fout,"CLASS %s\n",nd->class);*/
  /*** Write number, names, types, and values for constant parameters. ***/
  fprintf(fout,"NCONST %d\n",nd->nconst);
  for(i=0;i<nd->nconst;i++)
    fprintf(fout,"  %s %s\n",nd->cname[i],nd->cval[i]);
  /*** Write number and names and types for variable parameters. ***/
  fprintf(fout,"NVAR %d\n",nd->nvar);
  for(i=0;i<nd->nvar;i++)
    fprintf(fout,"  %s\n",nd->vname[i]);
  /*** Write event code tables ***/
  if (nd->ntable > 0){
    fprintf(fout,"NTABLE %d\n",nd->ntable);
    for(i=0;i<nd->ntable;i++){
      fprintf(fout,"%s %d\n",nd->table[i].tname,nd->table[i].n);
      for(j=0;j<nd->table[i].n;j++)
	fprintf(fout,"  %20s %6d\n",nd->table[i].name[j],nd->table[i].num[j]);
    }
  }
  /*** Write number of trials. ***/
  fprintf(fout,"NTRIAL %d\n",nd->ntrial);
  /*** Make big assumption here that all trials have same number of recs. ***/
  fprintf(fout,"NREC %d\n",nd->t[0].nrec);
  for(i=0;i<nd->ntrial;i++){
    fprintf(fout,"TRIAL %d\n",i);
    for(j=0;j<nd->t[i].nparam;j++)
      fprintf(fout,"  %s %s\n",nd->t[i].pname[j],nd->t[i].pval[j]);
    /*** Write all records for this trial. ***/
    for(j=0;j<nd->t[i].nrec;j++){
      fprintf(fout,"%s ",nd->t[i].r[j].name);
      fprintf(fout,"%.2f ",nd->t[i].r[j].sampling);
      fprintf(fout,"%d ",nd->t[i].r[j].tn);
      fprintf(fout,"%d\n",nd->t[i].r[j].n);
      if (nd->t[i].r[j].rtype == 0){
	print_iarray_80_column(fout,nd->t[i].r[j].p,nd->t[i].r[j].n);
      }else if (nd->t[i].r[j].rtype == 1){ /*** USE "tn" ***/
	for(k=0;k<nd->t[i].r[j].n;k++)
	  fprintf(fout,"%f ",nd->t[i].r[j].x[k]);
	fprintf(fout,"\n");
      }else if (nd->t[i].r[j].rtype == 2){ /*** USE "n" ***/
	for(k=0;k<nd->t[i].r[j].n;k++)
	  fprintf(fout,"%d %f ",nd->t[i].r[j].p[k],nd->t[i].r[j].x[k]);
	fprintf(fout,"\n");
      }else if (nd->t[i].r[j].rtype == 3){ /*** USE "n" ***/

	// WYETH - added condition to allow NULL EC table, Jan 5, 2017
	if (nd->t[i].r[j].ect != NULL)
	  ectnum = nd->t[i].r[j].ect->tnum;
	else
	  ectnum = -1;  // Indicate that there is no table
	fprintf(fout,"TABLENUM %d\n",ectnum);
	//fprintf(fout,"TABLENUM %d\n",nd->t[i].r[j].ect->tnum);

	for(k=0;k<nd->t[i].r[j].n;k++)
	  fprintf(fout,"%d %d ",nd->t[i].r[j].p[k],nd->t[i].r[j].ec[k]);
	fprintf(fout,"\n");
      }else
	exit_error("WRITE_NDATA_TXT","Unknown record type");
    }
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                               WRITE_NDATA_T1                              */
/*                                                                           */
/*****************************************************************************/
void write_ndata_t1(outfile,name,nd,fdig,nrecwrite)
     char outfile[],name[];
     struct ndata_struct *nd;
     int fdig;                  // Number of digits if float
     int nrecwrite;             // Number of records to write, all if <1
{
  FILE *fopen(),*fout;
  int i,j,k;
  int n,nrec;

  printf("  WRITE_NDATA_T1\n");

  nrec = nd->t[0].nrec;
  if ((nrecwrite < 1) || (nrecwrite > nrec))
    nrecwrite = nrec;
    
  printf("    %d records in first trial.\n",nrec);

  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** File %s\n",outfile);
    exit_error("WRITE_NDATA_T1","Cannot open file");
  }

  /*** Write number and names and types for variable parameters. ***/
  fprintf(fout,"Name %s\n",name);

  fprintf(fout,"Start %d\n",nd->t[0].r[0].t0);
  fprintf(fout,"Duration %d\n",nd->t[0].r[0].tn);
  
  /* WYETH - for TF data */
  /***
      fprintf(fout,"Start 0\n");
      fprintf(fout,"Duration 5138\n");
  ***/
  
  fprintf(fout,"Sampling %.1f\n",nd->t[0].r[0].sampling);
  if (nrec > 1)
    fprintf(fout,"Channels %d\n",nrec);
  fprintf(fout,"Params");
  for(i=0;i<nd->nvar;i++)
    fprintf(fout," %s",nd->vname[i]);
  fprintf(fout,"\n");
  fprintf(fout,"Trials %d\n",nd->ntrial);

  /*** WYETH - should check that all trials have start/dur ***/

  for(i=0;i<nd->ntrial;i++){
    /*** Print var param values all on one line ***/
    fprintf(fout,"T %d ",i+1);

    /* 
       if (nd->t[i].nparam > 0)
       fprintf(fout,"%s",nd->t[i].pval[0]);
       for(j=1;j<nd->t[i].nparam;j++)
    */
    for(j=0;j<nd->t[i].nparam;j++)
      fprintf(fout," %s",nd->t[i].pval[j]);

    fprintf(fout,"\n");

    /*** Write all records for this trial. ***/
    for(j=0;j<nrecwrite;j++){
      if (nd->t[i].r[j].rtype == 0){
	n = nd->t[i].r[j].n;
	fprintf(fout,"R %d",n);
	for(k=0;k<n;k++)
	  fprintf(fout," %d",nd->t[i].r[j].p[k]);
	fprintf(fout,"\n");
      }else if (nd->t[i].r[j].rtype == 1){
	n = nd->t[i].r[j].n;
	fprintf(fout,"R %d",n);
	for(k=0;k<n;k++)
	  if (fdig == 1)
	    fprintf(fout," %.1f",nd->t[i].r[j].x[k]);
	  else if (fdig == 2)
	    fprintf(fout," %.2f",nd->t[i].r[j].x[k]);
	  else if (fdig == 3)
	    fprintf(fout," %.3f",nd->t[i].r[j].x[k]);
	  else if (fdig == 4)
	    fprintf(fout," %.4f",nd->t[i].r[j].x[k]);
	  else
	    fprintf(fout," %f",nd->t[i].r[j].x[k]);
	fprintf(fout,"\n");
      }else{
	printf("  Unknown record type:   %d\n",nd->t[i].r[j].rtype);
	printf("  *** Ignored.\n");
      }
    }
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                             NDATA_UTIL_READ_T1                            */
/*                                                                           */
/*****************************************************************************/
void ndata_util_read_t1(infile,rnd)
     char infile[];
     struct ndata_struct **rnd;
{
  FILE *fopen(),*fin;
  int i,j,k;
  int n,nrec,ns,nconst,nvar,start,duration,*s,done,nn;
  float sampling;
  char t[SLEN],ts[SLEN],**cname,**cval,*ctype,*vtype,**vname,***vval,*tc;
  struct ndata_struct *nd;

  printf("  NDATA_UTIL_READ_T1\n");
  printf("    Assumes all records are spike times\n");

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** Cannot open file %s.\n",infile);
    exit_error("NDATA_UTIL_READ_T1","Cannot open file");
  }
  printf("    Reading %s\n",infile);
  nd = (struct ndata_struct *)myalloc(sizeof(struct ndata_struct));
  nd->class = strdup("T1");
  nd->ntable = 0;
  nd->table = NULL;

  /*** Write number and names and types for variable parameters. ***/
  nconst = 0;
  cname = (char **)myalloc(1000*sizeof(char *));
  cval  = (char **)myalloc(1000*sizeof(char *));
  ctype = (char  *)myalloc(1000*sizeof(char));
  vtype = NULL; /* to suppress compiler warning */
  nn = fscanf(fin,"%s",t);
  duration = -1;
  start = -1;
  sampling = 0.0;
  nvar = -1;
  while(strcmp(t,"Trials")!=0){
    if (strcmp(t,"Params")==0){
      tc = fgets(ts,SLEN,fin); /* Read until the carriage return. */
      ts[strlen(ts)-1] = '\0'; /* Remove final CR */
      get_items_from_string(ts,&vname,&nvar);

      if (nvar > 0) /*** Print var param names ***/
	printf("    Var params: ");
      for(j=0;j<nvar;j++)
	printf(" %s",vname[j]);
      if (nvar > 0)
	printf("\n");

      vtype = (char *)myalloc(nvar*sizeof(char));
      for(i=0;i<nvar;i++)
	vtype[i] = 'c';
    }else if (strcmp(t,"Duration")==0){
      nn = fscanf(fin,"%d",&duration);
    }else if (strcmp(t,"Start")==0){
      nn = fscanf(fin,"%d",&start);
    }else if (strcmp(t,"Sampling")==0){
      nn = fscanf(fin,"%f",&sampling);
    }else{
      cname[nconst] = strdup(t);
      nn = fscanf(fin,"%s",t);
      cval[nconst] = strdup(t);
      ctype[nconst] = 'c';
      nconst += 1;
      if (nconst > 1000)
	exit_error("NDATA_UTIL_READ_T1","Too many const params");
    }
    nn = fscanf(fin,"%s",t);
    /*printf("==>%s<==\n",t);*/
  }
  nn = fscanf(fin,"%d",&n);
  printf("    %d trials found.\n",n);

  if (start == -1)
    exit_error("NDATA_UTIL_READ_T1","Start not found, or is -1");
  if (duration == -1)
    exit_error("NDATA_UTIL_READ_T1","Duration not found, or is -1");
  if (sampling == 0.0)
    exit_error("NDATA_UTIL_READ_T1","Sampling not found, or is 0");
  if (nvar == -1)
    exit_error("NDATA_UTIL_READ_T1","Params not found");

  init_set_const_param(nd,nconst,cname,ctype,cval);

  /*** Check all trials, and determine number of records ***/
  nrec = 0; /* to suppress compiler warning */
  nn = fscanf(fin,"%s",t);
  for(i=0;i<n;i++){
    if (strcmp(t,"T")!=0){
      printf("  at trial %d, found %s\n",i+1,t);
      exit_error("NDATA_UTIL_READ_T1","Expecting T");
    }
    nn = fscanf(fin,"%d",&k);
    if (k != i+1)
      printf("  *** WARNING, trial %d is numbered %d\n",i+1,k);
    for(j=0;j<nvar;j++) /* Read through var params */
      nn = fscanf(fin,"%s",t);
    nn = fscanf(fin,"%s",t);
    if (strcmp(t,"R")!=0){
      printf("  at trial %d\n",i);
      exit_error("NDATA_UTIL_READ_T1","Expecting R");
    }
    nn = fscanf(fin,"%d",&ns);
    /*printf("  Trial %d:  %d data points\n",i,ns);*/
    for(j=0;j<ns;j++) /* Read through data */
      nn = fscanf(fin,"%s",t);
    
    k = 1; /* Count number of records for this trial */
    if (fscanf(fin,"%s",t) != EOF){
      done = 0;
      while ((strcmp(t,"R")==0) && (!done)){
	nn = fscanf(fin,"%d",&ns);
	for(j=0;j<ns;j++){ /* Read through data */
	  nn = fscanf(fin,"%s",t);
	}
	k += 1;
	if (fscanf(fin,"%s",t) == EOF)
	  done = 1;
      }
    }
    if (i==0)
      nrec = k;
    else
      if (k != nrec){
	printf("  at trial %d\n",i);
	exit_error("NDATA_UTIL_READ_T1","Number of records varies");
      }
  }
  fclose(fin);

  printf("    %d records per trial\n",nrec);

  vval = (char ***)myalloc(n*sizeof(char **));
  nd->ntrial = n;
  nd->t = (struct ndtrial_struct *)myalloc(n*sizeof(struct ndtrial_struct));

  /*** Reopen the file, and skip over the header ***/
  fin = fopen(infile,"r");
  nn = fscanf(fin,"%s",t);
  while(strcmp(t,"Trials")!=0){
    tc = fgets(ts,SLEN,fin); /* Read until the carriage return. */
    nn = fscanf(fin,"%s",t);
  }
  nn = fscanf(fin,"%s",t); /* Skip the number of trials */

  for(i=0;i<n;i++){
    nd->t[i].tcode = i;
    nn = fscanf(fin,"%s",t);  /* T */

    nn = fscanf(fin,"%d",&k); /* trial number */
    nd->t[i].tref = k;

    vval[i] = (char **)myalloc(nvar*sizeof(char *));
    for(j=0;j<nvar;j++){ /* Read var param values */
      nn = fscanf(fin,"%s",t);
      vval[i][j] = strdup(t);
    }
    nd->t[i].nrec = nrec;
    nd->t[i].r = (struct ndrec_struct *)myalloc(nrec*
						sizeof(struct ndrec_struct));
    for(j=0;j<nrec;j++){
      nd->t[i].r[j].rtype = 0;
      sprintf(ts,"unit%d",j);
      nd->t[i].r[j].name = strdup(ts);
      nd->t[i].r[j].rcode = 0;
      nd->t[i].r[j].sampling = sampling;
      nd->t[i].r[j].t0 = start;
      nd->t[i].r[j].tn = duration;
      nn = fscanf(fin,"%s",t);  /* R */
      nn = fscanf(fin,"%d",&ns); /* number of datum */
      nd->t[i].r[j].n = ns;
      nd->t[i].r[j].p = (int *)myalloc(ns*sizeof(int));
      s = nd->t[i].r[j].p; /* Pointer */
      for(k=0;k<ns;k++){
	nn = fscanf(fin,"%d",&(s[k]));
	/*printf("%d ",s[k]);*/
      }
    }
  }

  fclose(fin);
  init_set_var_param(nd,n,nvar,vname,vtype,vval);
  
  *rnd = nd;
}
/**************************************-**************************************/
/*                                                                           */
/*                              NDU_READ_T2_CONST                            */
/*                                                                           */
/*****************************************************************************/
void ndu_read_t2_const(fin,nd,name)
     FILE *fin;
     struct ndata_struct *nd;
     char *name;               // This name was already read in.
{
  int k;
  int nmax,ns;
  char **cname,**cval,*ctype,ts[SLEN],tval[SLEN];

  printf("  NDU_READ_T2_CONST\n");

  nmax = 10000;  // Arbitrary limit
  cname = (char **)myalloc(nmax*sizeof(char *));
  cval  = (char **)myalloc(nmax*sizeof(char *));
  ctype = (char  *)myalloc(nmax*sizeof(char));

  k = 0;
  if (strlen(name) > 0){
    cname[0] = strdup("Name");
    cval[0]  = strdup(name);
    ctype[0] = 'c';
    k += 1;
  }

  ns = fscanf(fin,"%s",ts);  // name of first const param
  while(strcmp(ts,"EndConst")!=0){
    ns = fscanf(fin,"%s",tval);  // value of const param

    cname[k] = strdup(ts);
    cval[k]  = strdup(tval);
    if (is_number_string(tval)){
      if (is_int_string(tval)){
	ctype[k] = 'i';  // Integer
      }else{
	ctype[k] = 'f';  // Float
      }
    }else{
	ctype[k] = 'c';  // Char string
    }

    k += 1;
    if (k >= nmax)
      exit_error("NDU_READ_T2_CONST","Too many 'const' parameters");

    ns = fscanf(fin,"%s",ts);  // Read next name, or "EndConst"
  }

  init_set_const_param(nd,k,cname,ctype,cval);

  printf("    %d constant parameters, including 'Name'\n",k);

  myfree(cname);
  myfree(cval);
  myfree(ctype);
}
/**************************************-**************************************/
/*                                                                           */
/*                              NDU_READ_T2_TABLE                            */
/*                                                                           */
/*****************************************************************************/
void ndu_read_t2_table(fin,nd)
     FILE *fin;
     struct ndata_struct *nd;
{
  int i,k;
  int nmax,*tcode,ns;
  char **tstr,ts[SLEN];

  printf("  NDU_READ_T2_TABLE\n");

  nmax = 10000;  // Arbitrary limit
  tcode = (int   *)myalloc(nmax*sizeof(int));
  tstr  = (char **)myalloc(nmax*sizeof(char *));

  k = 0;
  ns = fscanf(fin,"%s",ts);  // Code number, read as string in case "EndTable"
  while(strcmp(ts,"EndTable")!=0){

    tcode[k] = atoi(ts);
    ns = fscanf(fin,"%s",ts);  // value of const param
    tstr[k] = strdup(ts);
    k += 1;
    if (k >= nmax)
      exit_error("NDU_READ_T2_TABLE","Too many events in table");

    ns = fscanf(fin,"%s",ts);  // Read next code, or "EndTable"
  }

  printf("    %d codes in table\n",k);

  nd->table = (struct ect_struct *)myalloc(sizeof(struct ect_struct));
  nd->table->tname = strdup("Default");
  nd->table->tnum = 0;
  nd->table->n = k;
  nd->table->name = (char **)myalloc(k*sizeof(char *));
  nd->table->num  = (int   *)myalloc(k*sizeof(int));
  for(i=0;i<k;i++){
    nd->table->name[i] = tstr[i];
    nd->table->num[i] = tcode[i];
  }
  nd->ntable = 1;

  myfree(tstr);
  myfree(tcode);
}
/**************************************-**************************************/
/*                                                                           */
/*                             NDATA_UTIL_READ_T2                            */
/*                                                                           */
/*  *** WYETH - THIS SHOULD PROBABLY SHARE CODE w/ READ_T1 ****              */
/*                                                                           */
/*****************************************************************************/
void ndata_util_read_t2(infile,rnd)
     char infile[];
     struct ndata_struct **rnd;
{
  FILE *fopen(),*fin;
  int i,j,k;
  int n,nrec,ns,nconst,nvar,start,duration,*s,*e,done,tref,si,ei,nn;
  float sampling;
  char t[SLEN],ts[SLEN],**cname,**cval,*ctype,*vtype,**vname,***vval,*tc;
  char name[SLEN];
  struct ndata_struct *nd;

  printf("  NDATA_UTIL_READ_T2\n");

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** Cannot open file %s.\n",infile);
    exit_error("NDATA_UTIL_READ_T2","Cannot open file");
  }
  printf("    Reading %s\n",infile);
  nd = (struct ndata_struct *)myalloc(sizeof(struct ndata_struct));
  nd->class = strdup("T2");
  nd->ntable = 0;
  nd->table = NULL;

  /*** Write number and names and types for variable parameters. ***/

  /*
  nconst = 0;
  cname = (char **)myalloc(1000*sizeof(char *));
  cval = (char **)myalloc(1000*sizeof(char *));
  ctype = (char *)myalloc(1000*sizeof(char));
  */


  vtype = NULL; // to suppress compiler warning
  nn = fscanf(fin,"%s",t);
  duration = -1;
  start = -1;
  sampling = 0.0;
  name[0] = '\0';  // Is this NULL ?
  nvar = -1;
  while(strcmp(t,"Trials")!=0){
    if (strcmp(t,"Params")==0){
      tc = fgets(ts,SLEN,fin); /* Read until the carriage return. */
      ts[strlen(ts)-1] = '\0'; /* Remove final CR */
      get_items_from_string(ts,&vname,&nvar);

      if (nvar > 0) // Print var param names
	printf("  Var params:\n    ");
      for(j=0;j<nvar;j++)
	printf(" %s",vname[j]);
      if (nvar > 0)
	printf("\n");

      vtype = (char *)myalloc(nvar*sizeof(char));
      for(i=0;i<nvar;i++)
	vtype[i] = 'c';
    }else if (strcmp(t,"Duration")==0){
      nn = fscanf(fin,"%d",&duration);
    }else if (strcmp(t,"Name")==0){
      nn = fscanf(fin,"%s",name);
    }else if (strcmp(t,"Start")==0){
      nn = fscanf(fin,"%d",&start);
    }else if (strcmp(t,"Sampling")==0){
      nn = fscanf(fin,"%f",&sampling);
    }else if (strcmp(t,"BeginConst")==0){
      ndu_read_t2_const(fin,nd,name);
    }else if (strcmp(t,"BeginTable")==0){
      ndu_read_t2_table(fin,nd);
    }else{
      printf("  Found: %s\n",t);
      exit_error("NDATA_UTIL_READ_T2","Unknown item");
    }
    nn = fscanf(fin,"%s",t);
    /*printf("==>%s<==\n",t);*/
  }
  nn = fscanf(fin,"%d",&n);
  printf("  %d trials found.\n",n);

  if (start == -1)
    exit_error("NDATA_UTIL_READ_T2","Start not found, or is -1");
  if (duration == -1)
    exit_error("NDATA_UTIL_READ_T2","Duration not found, or is -1");
  if (sampling == 0.0)
    exit_error("NDATA_UTIL_READ_T2","Sampling not found, or is 0");
  if (nvar == -1)
    exit_error("NDATA_UTIL_READ_T2","Params not found");

  //init_set_const_param(nd,nconst,cname,ctype,cval);

  /*** Check all trials, and determine number of records ***/
  nrec = 0; // to suppress compiler warning
  nn = fscanf(fin,"%s",t);
  for(i=0;i<n;i++){
    if (strcmp(t,"T")!=0){
      printf("  at trial %d, found %s\n",i+1,t);
      exit_error("NDATA_UTIL_READ_T2","Expecting T");
    }
    nn = fscanf(fin,"%d",&k);

    /*  THIS IS NO LONGER DESIRABLE TO ENFORCE
    if (k != i+1)
      printf("  *** WARNING, trial %d is numbered %d\n",i+1,k);
    */

    nn = fscanf(fin,"%d",&tref);

    for(j=0;j<nvar;j++) // Read through var params
      nn = fscanf(fin,"%s",t);

    //
    //  WYETH - handles R0 and R3, can extend to R1 and R2 later
    //
    done = 0;
    k = 0;  // Count records for this trial
    while(!done){
      if (fscanf(fin,"%s",t) == EOF){
	done = 1;
      }else if (t[0] == 'R'){
	nn = fscanf(fin,"%d",&ns);
	//printf("  Trial %d:  %d data points\n",i,ns);
	if (strcmp(t,"R3")==0)
	  ns *= 2;  // Two values for each event
	for(j=0;j<ns;j++) // Skip over the record data
	  nn = fscanf(fin,"%s",t);
	k += 1;
      }else if (t[0] == 'T'){
	done = 1;  // ti += 1;  // Next trial
      }else{
	printf("  at trial %d, found:  %s\n",i,t);
	exit_error("NDATA_UTIL_READ_T2","Expecting R...");
      }
    }

    if (i==0)
      nrec = k;
    else
      if (k != nrec){
	printf("  at trial %d\n",i);
	exit_error("NDATA_UTIL_READ_T2","Number of records varies");
      }
  }
  fclose(fin);

  printf("    %d records per trial\n",nrec);

  vval = (char ***)myalloc(n*sizeof(char **));
  nd->ntrial = n;
  nd->t = (struct ndtrial_struct *)myalloc(n*sizeof(struct ndtrial_struct));

  // Reopen the file, and skip over the header
  fin = fopen(infile,"r");
  nn = fscanf(fin,"%s",t);
  while(strcmp(t,"Trials")!=0){
    tc = fgets(ts,SLEN,fin); // Read until the carriage return
    nn = fscanf(fin,"%s",t);
  }
  nn = fscanf(fin,"%s",t); // Skip the number of trials

  for(i=0;i<n;i++){
    nn = fscanf(fin,"%s",t);  // T

    nn = fscanf(fin,"%d",&k); // trial number
    nd->t[i].tcode = k;

    nn = fscanf(fin,"%d",&k); // time reference
    nd->t[i].tref = k;

    vval[i] = (char **)myalloc(nvar*sizeof(char *));
    for(j=0;j<nvar;j++){ // Read var param values
      nn = fscanf(fin,"%s",t);
      vval[i][j] = strdup(t);
    }
    nd->t[i].nrec = nrec;
    nd->t[i].r = (struct ndrec_struct *)myalloc(nrec*
						sizeof(struct ndrec_struct));
    si = ei = 0;  // Index of spike and event records
    for(j=0;j<nrec;j++){

      nn = fscanf(fin,"%s",t);  // "R..."
      if (strcmp(t,"R0")==0){
	nd->t[i].r[j].rtype = 0;
	sprintf(ts,"unit%d",si);
	si += 1;
	nd->t[i].r[j].name = strdup(ts);
	nd->t[i].r[j].rcode = 0;
	nd->t[i].r[j].sampling = sampling;
	nd->t[i].r[j].t0 = start;
	nd->t[i].r[j].tn = duration;
	nn = fscanf(fin,"%d",&ns); /* number of datum */
	nd->t[i].r[j].n = ns;
	nd->t[i].r[j].p = (int *)myalloc(ns*sizeof(int));
	s = nd->t[i].r[j].p; /* Pointer */
	for(k=0;k<ns;k++){
	  nn = fscanf(fin,"%d",&(s[k]));
	  /*printf("%d ",s[k]);*/
	}
      }else if (strcmp(t,"R3")==0){

	nd->t[i].r[j].rtype = 3;
	sprintf(ts,"event%d",ei);
	ei += 1;
	nd->t[i].r[j].name = strdup(ts);
	nd->t[i].r[j].rcode = 0;
	nd->t[i].r[j].sampling = sampling;
	nd->t[i].r[j].t0 = start;
	nd->t[i].r[j].tn = duration;
	nn = fscanf(fin,"%d",&ns); // number of datum
	nd->t[i].r[j].n = ns;
	nd->t[i].r[j].ect = nd->table;  // Pointer to SINGLE TABLE
	nd->t[i].r[j].p = (int *)myalloc(ns*sizeof(int));
	nd->t[i].r[j].ec = (int *)myalloc(ns*sizeof(int));
	s = nd->t[i].r[j].p; // Pointer
	e = nd->t[i].r[j].ec; // Pointer
	for(k=0;k<ns;k++){
	  nn = fscanf(fin,"%d %d",&(e[k]),&(s[k]));
	  //printf("%d %d  ",e[k],s[k]);
	}
	//printf("\n");

      }else{
	printf("  Record type:  %s\n",t);
	exit_error("NDATA_UTIL_READ_T2","Bad record type");
      }
    }
  }

  fclose(fin);
  init_set_var_param(nd,n,nvar,vname,vtype,vval);
  
  *rnd = nd;
}
/**************************************-**************************************/
/*                                                                           */
/*                             READ_PSYCH_TEXT_NDATA                         */
/*                                                                           */
/*  Read a text file written by the OpenGL psych program 'psych.c'.          */
/*                                                                           */
/*****************************************************************************/
void read_psych_text_ndata(infile,qtflag,rnd)
     char infile[];
     int qtflag; /* Query for var param types */
     struct ndata_struct **rnd;
{
  FILE *fopen(),*fin;
  int i,j,k;
  struct ndata_struct *nd;
  int ni,ti,nresp,nvar0,ns;
  char t[SLEN],ts[SLEN],**ilist,c0,*tc;

  printf("  READ_PSYCH_TEXT_NDATA\n");
  
  c0 = '\0';

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** Cannot open file %s.\n",infile);
    *rnd = NULL;
    return;
  }
  printf("    Reading %s\n",infile);
  nd = (struct ndata_struct *)myalloc(sizeof(struct ndata_struct));
  nd->class = strdup("PSYCH");
  nd->ntable = 0;

  /*** Find the number and names of the variable parameters ***/
  ns = fscanf(fin,"%s",t);
  while(strcmp(t,"UNIQUE_TRIAL_LIST")!=0){
    if (strcmp(t,"nresponse")==0)
      ns = fscanf(fin,"%d",&nresp);
    ns = fscanf(fin,"%s",t);
  }

  tc = fgets(ts,SLEN,fin); // Read until the carriage return
  tc = fgets(ts,SLEN,fin); // Read until the carriage return
  get_items_from_string(ts,&ilist,&ni);
  nvar0 = ni;
  nd->nvar = nvar0 + nresp*2 + 1; /* 2=time and response, 1=stimulus time */
  nd->vname = (char **)myalloc(nd->nvar*sizeof(char *));
  nd->vtype = (char *)myalloc(nd->nvar*sizeof(char));
  k = 0; /* Counter of var params */
  nd->vname[k] = strdup("psych_display_time");
  nd->vtype[k] = 'i';

  k += 1;
  for(i=0;i<nresp;i++){
    //sprintf(t,"response%d_time\0",i+1);
    sprintf(t,"response%d_time%c",i+1,c0);
    nd->vname[k] = strdup(t);
    nd->vtype[k] = 'i';
    k += 1;
    //sprintf(t,"response%d\0",i+1);
    sprintf(t,"response%d%c",i+1,c0);
    nd->vname[k] = strdup(t);
    if (qtflag){
      printf("  Enter var param type (i,f,c) %s: ",nd->vname[k]);
      fflush(stdout);
      ns = scanf("%s",t);
      nd->vtype[k] = t[0];
    }else
      nd->vtype[k] = 'c';
    k += 1;
  }
  for(i=0;i<ni;i++){
    nd->vname[k] = strdup(ilist[i]);
    if (qtflag){
      printf("  Enter var param type (i,f,c) %s: ",nd->vname[k]);
      fflush(stdout);
      ns = scanf("%s",t);
      nd->vtype[k] = t[0];
    }else
      nd->vtype[k] = 'c';
    k += 1;
  }
  fclose(fin);
  printf("    %d var params found\n",nd->nvar);
  for(i=0;i<nd->nvar;i++)
    printf("      %d  %s\n",i,nd->vname[i]);

  free_2d_carray(ilist,ni);
  if (k != nd->nvar)
    exit_error("*********","SHOULD NOT HAPPEN\n");

  /*** Count const params ***/
  fin = fopen(infile,"r");

  nd->nconst = 0;
  tc = fgets(ts,SLEN,fin); // Read until the carriage return
  get_items_from_string(ts,&ilist,&ni);
  while((strcmp(ilist[0],"ptable")!=0)&&
	(strcmp(ilist[0],"UNIQUE_TRIAL_LIST")!=0)){
    /*printf("ilist[0] = %s\n",ilist[0]);*/
    if (search_2d_carray(nd->vname,ilist[0],nd->nvar) == -1)
      nd->nconst += 1;
    free_2d_carray(ilist,ni);
    tc = fgets(ts,SLEN,fin); // Read until the carriage return
    get_items_from_string(ts,&ilist,&ni);
  }
  free_2d_carray(ilist,ni);

  /*** Find and count trial data ***/
  ns = fscanf(fin,"%s",t);
  while(strcmp(t,"TRIAL_DATA")!=0)
    ns = fscanf(fin,"%s",t);
  nd->ntrial = 0;
  tc = fgets(ts,SLEN,fin); // Carriage return
  while(fgets(ts,SLEN,fin)!=NULL)
    nd->ntrial += 1;

  fclose(fin);
  printf("    %d const params found\n",nd->nconst);
  printf("    %d trials found\n",nd->ntrial);

  /*** Fill in const param values and trial data ***/
  nd->cname = (char **)myalloc(nd->nconst*sizeof(char *));
  nd->ctype = (char *)myalloc(nd->nconst*sizeof(char));
  nd->cval = (char **)myalloc(nd->nconst*sizeof(char *));

  fin = fopen(infile,"r");
  k = 0;
  tc = fgets(ts,SLEN,fin); /* Read until the carriage return. */
  get_items_from_string(ts,&ilist,&ni);
  while((strcmp(ilist[0],"ptable")!=0)&&
	(strcmp(ilist[0],"UNIQUE_TRIAL_LIST")!=0)){
    if (search_2d_carray(nd->vname,ilist[0],nd->nvar) == -1){
      nd->cname[k] = strdup(ilist[0]);
      nd->ctype[k] = 'c';
      nd->cval[k] = make_string_from_items(ilist+1,ni-1);
      k += 1;
    }
    free_2d_carray(ilist,ni);
    tc = fgets(ts,SLEN,fin); /* Read until the carriage return. */
    get_items_from_string(ts,&ilist,&ni);
  }
  free_2d_carray(ilist,ni);

  /*** Find and count trial data ***/
  nd->t = (struct ndtrial_struct *)myalloc(nd->ntrial*
					 sizeof(struct ndtrial_struct));
  ns = fscanf(fin,"%s",t);
  while(strcmp(t,"TRIAL_DATA")!=0)
    ns = fscanf(fin,"%s",t);
  tc = fgets(ts,SLEN,fin); /* Carriage return */
  for(i=0;i<nd->ntrial;i++){
    nd->t[i].nrec = 0;
    nd->t[i].tcode = i;
    nd->t[i].nparam = nd->nvar;
    nd->t[i].pname = (char **)myalloc(nd->t[i].nparam*sizeof(char *));
    nd->t[i].pval = (char **)myalloc(nd->t[i].nparam*sizeof(char *));
    for(j=0;j<nd->t[i].nparam;j++)
      nd->t[i].pname[j] = strdup(nd->vname[j]);
    
    ns = fscanf(fin,"%d",&ti); /* Trial number */
    nd->t[i].tref = ti;

    ns = fscanf(fin,"%s",t); /* stimulus display time, measured */
    k = 0;
    nd->t[i].pval[k] = strdup(t);
    k += 1;

    k += 2*nresp; /*** Jump over responses/times ***/
    for(j=0;j<nvar0;j++){
      ns = fscanf(fin,"%s",t);
      nd->t[i].pval[k] = strdup(t);
      k += 1;
    }
    k = 1;
    for(j=0;j<nresp;j++){
      ns = fscanf(fin,"%s",t); /* response time */
      nd->t[i].pval[k] = strdup(t);
      k += 1;
      ns = fscanf(fin,"%s",t); /* response */
      nd->t[i].pval[k] = strdup(t);
      k += 1;
    }
  }
  fclose(fin);

  *rnd = nd;
}
/**************************************-**************************************/
/*                                                                           */
/*                                NDATA_FREE_ECT                             */
/*                                                                           */
/*****************************************************************************/
void ndata_free_ect(ect)
     struct ect_struct *ect;
{
  int i;

  myfree(ect->tname);
  myfree(ect->num);

  for(i=0;i<ect->n;i++){
    if (ect->name[i] != NULL)
      myfree(ect->name[i]);
  }
  myfree(ect->name);


  /*
ndata_util.c:1298: warning: cast to pointer from integer of different size
ndata_util.c: In function `free_ndata':
  */


  free_2d_carray(ect->name,ect->n);
  myfree(ect);
}
/**************************************-**************************************/
/*                                                                           */
/*                                 GET_NDATA                                 */
/*                                                                           */
/*  This 'nd' struct can be safely freed using 'free_ndata'.                 */
/*                                                                           */
/*****************************************************************************/
struct ndata_struct *get_ndata()
{
  struct ndata_struct *nd;
  
  nd = (struct ndata_struct *)myalloc(sizeof(struct ndata_struct));

  nd->class = (char *)NULL;
  nd->cname = (char **)NULL;
  nd->nconst = 0;
  nd->ctype = (char *)NULL;
  nd->cval = (char **)NULL;

  nd->nvar = 0;
  nd->vname = (char **)NULL;
  nd->vtype = (char *)NULL;

  nd->ntable = 0;
  nd->table = (struct ect_struct *)NULL;

  nd->ntrial = 0;
  nd->t = (struct ndtrial_struct *)NULL;

  return nd;
}
/**************************************-**************************************/
/*                                                                           */
/*                               FREE_NDATA_REC                              */
/*                                                                           */
/*****************************************************************************/
void free_ndata_rec(rp)
     struct ndrec_struct *rp;   // Pointer to trial
{
  if (rp->name != NULL)
    myfree(rp->name);

  if (rp->rtype == 0){
    if (rp->p != NULL)
      myfree(rp->p);
  }else if (rp->rtype == 1){
    if (rp->x != NULL)
      myfree(rp->x);
  }else if (rp->rtype == 2){
    if (rp->p != NULL)
      myfree(rp->p);
    if (rp->x != NULL)
      myfree(rp->x);
  }else if (rp->rtype == 3){
    if (rp->p != NULL)
      myfree(rp->p);
    if (rp->ec != NULL)
      myfree(rp->ec);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              FREE_NDATA_TRIAL                             */
/*                                                                           */
/*  NOTE:  This does *NOT* traverse or free the 'prev' and 'next' links.     */
/*                                                                           */
/*****************************************************************************/
void free_ndata_trial(tp)
     struct ndtrial_struct *tp;   // Pointer to trial
{
  int i;

  free_2d_carray(tp->pname,tp->nparam);
  free_2d_carray(tp->pval,tp->nparam);

  for(i=0;i<tp->nrec;i++)
    free_ndata_rec(&(tp->r[i]));

  myfree(tp->r);
}
/**************************************-**************************************/
/*                                                                           */
/*                                 FREE_NDATA                                */
/*                                                                           */
/*****************************************************************************/
void free_ndata(nd)
     struct ndata_struct *nd;
{
  int i,j;
  int debugflag;
  
  debugflag = 0;

  if (debugflag) printf("  FREE_NDATA\n");
  
  myfree(nd->class);
  if (debugflag) printf("  Here 01\n");
  free_2d_carray(nd->cname,nd->nconst);
  if (debugflag) printf("  Here 02\n");
  myfree(nd->ctype);
  if (debugflag) printf("  Here 03\n");
  free_2d_carray(nd->cval,nd->nconst);
  if (debugflag) printf("  Here 04\n");

  free_2d_carray(nd->vname,nd->nvar);
  if (debugflag) printf("  Here 05\n");
  myfree(nd->vtype);
  if (debugflag) printf("  Here 06\n");

  for(i=0;i<nd->ntable;i++){ // This is an array, *not* an array of pointers
    myfree(nd->table[i].tname);
    myfree(nd->table[i].num);

    /*
  ndata_util.c:1355: warning: cast to pointer from integer of different size
  ar: creating ../lib/libndata_util.a
    */

    free_2d_carray(nd->table[i].name,nd->table[i].n);
  }
  if (debugflag) printf("  Here 07\n");
  if (nd->table != NULL)
    myfree(nd->table);

  if (debugflag) printf("  Here 08\n");

  for(i=0;i<nd->ntrial;i++){

    free_ndata_trial(&(nd->t[i]));  // Free contents of trial

    /*  WYETH - THIS REPLACED BY LINE ABOVE on 2011 Aug 01
    free_2d_carray(nd->t[i].pname,nd->t[i].nparam);
    free_2d_carray(nd->t[i].pval,nd->t[i].nparam);
    for(j=0;j<nd->t[i].nrec;j++){
      myfree(nd->t[i].r[j].name);
      if (nd->t[i].r[j].rtype == 0){
	if (nd->t[i].r[j].p != NULL)
	  myfree(nd->t[i].r[j].p);
      }else if (nd->t[i].r[j].rtype == 1){
	myfree(nd->t[i].r[j].x);
      }else if (nd->t[i].r[j].rtype == 2){
	myfree(nd->t[i].r[j].p);
	myfree(nd->t[i].r[j].x);
      }else if (nd->t[i].r[j].rtype == 3){
	myfree(nd->t[i].r[j].p);
	myfree(nd->t[i].r[j].ec);
      }
    }
    myfree(nd->t[i].r);
    */
  }
  if (debugflag) printf("  Here 09\n");

  myfree(nd->t);
  if (debugflag) printf("  Here 10\n");

  myfree(nd);
  if (debugflag) printf("  Here 11\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                             FREE_NDATA_LINK                               */
/*                                                                           */
/*****************************************************************************/
void free_ndata_link(nd)
     struct ndata_struct *nd;
{
  int i,j;
  struct ndtrial_struct *tp,*tnext;
  
  /*printf("  FREE_NDATA_LINK\n");*/

  myfree(nd->class);
  free_2d_carray(nd->cname,nd->nconst);
  myfree(nd->ctype);
  free_2d_carray(nd->cval,nd->nconst);

  free_2d_carray(nd->vname,nd->nvar);
  myfree(nd->vtype);

  for(i=0;i<nd->ntable;i++){ /* This is an array, *not* an array of pointers */
    myfree(nd->table[i].tname);
    myfree(nd->table[i].num);
    free_2d_carray(nd->table[i].name,nd->table[i].n);
  }
  if (nd->table != NULL)
    myfree(nd->table);

  tp = nd->t;
  for(i=0;i<nd->ntrial;i++){
    free_2d_carray(tp->pname,tp->nparam);
    free_2d_carray(tp->pval,tp->nparam);
    for(j=0;j<tp->nrec;j++){
      myfree(tp->r[j].name);
      if (tp->r[j].rtype == 0){
	if (tp->r[j].p != NULL)
	  myfree(tp->r[j].p);
      }else if (tp->r[j].rtype == 1){
	myfree(tp->r[j].x);
      }else if (tp->r[j].rtype == 2){
	myfree(tp->r[j].p);
	myfree(tp->r[j].x);
      }else if (tp->r[j].rtype == 3){
	myfree(tp->r[j].p);
	myfree(tp->r[j].ec);
      }
    }
    myfree(tp->r);
    tnext = tp->next;
    myfree(tp);
    tp = tnext;
  }
  myfree(nd);
}
/**************************************-**************************************/
/*                                                                           */
/*                               NDATA_TRIAL_PRINT                           */
/*                                                                           */
/*****************************************************************************/
void ndata_trial_print_trial(t,maxr)
     struct ndtrial_struct *t;
     int maxr;  // Maximum number of records to show, -1 for all
{
  int i,k;
  int n;

  printf("tcode = %d\n",t->tcode);
  printf("tref = %d\n",t->tref);
  printf("nparam = %d\n",t->nparam);
  for(i=0;i<t->nparam;i++)
    printf("    %s = %s\n",t->pname[i],t->pval[i]);
  printf("  nrec = %d\n",t->nrec);

  // Number of records to print
  n = t->nrec;
  if ((maxr >= 1) && (maxr < t->nrec))
    n = maxr;
      
  for(i=0;i<n;i++){
    printf("%d %s %d %.2f %d %d %d\n",t->r[i].rtype,t->r[i].name,
	   t->r[i].rcode,t->r[i].sampling,t->r[i].t0,
	   t->r[i].tn,t->r[i].n);
    if (t->r[i].rtype == 0)
      for(k=0;k<t->r[i].n;k++)
	printf("%d ",t->r[i].p[k]);
    else if (t->r[i].rtype == 1) /*** USE "tn" ***/
      for(k=0;k<t->r[i].n;k++)
	printf("%.4f ",t->r[i].x[k]);
    else if (t->r[i].rtype == 2) /*** USE "n" ***/
      for(k=0;k<t->r[i].n;k++)
	printf("(%d %.4f) ",t->r[i].p[k],t->r[i].x[k]);
    else if (t->r[i].rtype == 3){ /*** USE "n" ***/
      if (t->r[i].ect != NULL)
	printf("  TABLENUM %d\n",t->r[i].ect->tnum);
      else
	printf("  TABLENUM %d\n",-1);
      for(k=0;k<t->r[i].n;k++)
	printf("(%d %d) ",t->r[i].p[k],t->r[i].ec[k]);
    }else
      exit_error("PRINT_NDATA","Unknown record type");
    printf("\n");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            NDATA_TRIAL_GET_TCODE                          */
/*                                                                           */
/*****************************************************************************/
int ndata_trial_get_tcode(t)
     struct ndtrial_struct *t;
{
  return t->tcode;
}
/**************************************-**************************************/
/*                                                                           */
/*                       NDATA_TRIAL_GET_VAR_PARAM_INDEX                     */
/*                                                                           */
/*  Return the index of the record "name" in the "k"th trial.  If not        */
/*  found, return -1.                                                        */
/*                                                                           */
/*****************************************************************************/
int ndata_trial_get_var_param_index(t,name)
     struct ndtrial_struct *t;
     char name[];
{
  return search_2d_carray(t->pname,name,t->nparam);
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_TRIAL_GET_RECORD_INDEX                      */
/*                                                                           */
/*  Return the index of the record "name" in the "k"th trial.  If not        */
/*  found, return -1.                                                        */
/*                                                                           */
/*****************************************************************************/
int ndata_trial_get_record_index(t,name)
     struct ndtrial_struct *t;
     char name[];
{
  int i;
  int n;

  i = 0;
  n = -1;
  while((i < t->nrec) && (n == -1))
    if (strcmp(t->r[i].name,name)==0)
      n = i;
    else
      i+= 1;

  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_TRIAL_GET_RECORD_PTR                        */
/*                                                                           */
/*****************************************************************************/
struct ndrec_struct *ndata_trial_get_record_ptr(t,name)
     struct ndtrial_struct *t;
     char name[];
{
  int i;
  struct ndrec_struct *r;

  r = NULL;

  i = ndata_trial_get_record_index(t,name);
  if (i >= 0)
    r = &(t->r[i]);

  return r;
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_TRIAL_GET_RECORD_TYPE                       */
/*                                                                           */
/*  Return the type of the record "name" in the "k"th trial.  If not        */
/*  found, return -1.                                                        */
/*                                                                           */
/*****************************************************************************/
int ndata_trial_get_record_type(t,name)
     struct ndtrial_struct *t;
     char name[];
{
  int i,k;

  i = ndata_trial_get_record_index(t,name);
  if (i >= 0)
    k = t->r[i].rtype;
  else
    k = -1;

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_TRIAL_GET_SAMPLING_CHAN                     */
/*                                                                           */
/*  Return the sampling for this channel.                                    */
/*                                                                           */
/*****************************************************************************/
float ndata_trial_get_sampling_chan(t,name)
     struct ndtrial_struct *t;
     char name[];
{
  int k;

  k = ndata_trial_get_record_index(t,name);
  if (k < 0)
    exit_error("NDATA_TRIAL_GET_SAMPLING_CHAN","Channel not found");

  return t->r[k].sampling;
}
/**************************************-**************************************/
/*                                                                           */
/*                          NDATA_TRIAL_GET_RCODE_CHAN                       */
/*                                                                           */
/*  Return the rcode for this channel on this trial.                         */
/*                                                                           */
/*****************************************************************************/
int ndata_trial_get_rcode_chan(t,name)
     struct ndtrial_struct *t;
     char name[];
{
  int k;

  k = ndata_trial_get_record_index(t,name);
  if (k < 0)
    exit_error("NDATA_TRIAL_GET_RCODE_CHAN","Channel not found");

  return t->r[k].rcode;
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_TRIAL_GET_START_CHAN                        */
/*                                                                           */
/*  Return the start time of this channel in raw sampling units.             */
/*                                                                           */
/*****************************************************************************/
int ndata_trial_get_start_chan(t,name)
     struct ndtrial_struct *t;
     char name[];
{
  int k;

  k = ndata_trial_get_record_index(t,name);
  if (k < 0)
    exit_error("NDATA_TRIAL_GET_START_CHAN","Channel not found");

  return t->r[k].t0;
}
/**************************************-**************************************/
/*                                                                           */
/*                        NDATA_TRIAL_GET_START_SEC_CHAN                     */
/*                                                                           */
/*  Return the start time of this channel in seconds.                        */
/*                                                                           */
/*****************************************************************************/
float ndata_trial_get_start_sec_chan(t,name)
     struct ndtrial_struct *t;
     char name[];
{
  int k;
  float start;

  k = ndata_trial_get_record_index(t,name);
  if (k < 0)
    exit_error("NDATA_TRIAL_GET_START_SEC_CHAN","Channel not found");
  start = (float)t->r[k].t0 / t->r[k].sampling;

  return start;
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_TRIAL_GET_DUR_SEC_CHAN                      */
/*                                                                           */
/*  Return the duration of this channel in seconds.                          */
/*                                                                           */
/*****************************************************************************/
float ndata_trial_get_dur_sec_chan(t,name)
     struct ndtrial_struct *t;
     char name[];
{
  int k;
  float dur;

  k = ndata_trial_get_record_index(t,name);
  if (k < 0)
    exit_error("NDATA_TRIAL_GET_DUR_SEC_CHAN","Channel not found");
  dur = (float)t->r[k].tn / t->r[k].sampling;

  return dur;
}
/**************************************-**************************************/
/*                                                                           */
/*                      NDATA_RECORD_GET_NTH_OCCUR_EVCO_TIME                 */
/*                                                                           */
/*  Return the time of the "occur"th occurrence of "ecode" in the specified  */
/*  record.  Return -1 if not found.                                         */
/*                                                                           */
/*****************************************************************************/
int ndata_record_get_nth_occur_evco_time(r,ecode,occur,rtime)
     struct ndrec_struct *r;
     int ecode,occur,*rtime;
{
  int k;

  k = get_index_search_nth_iarray(r->ec,r->n,ecode,occur);
  if (k >= 0)
    *rtime = r->p[k];
  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                      NDATA_TRIAL_GET_START_SEC_SCHAN                      */
/*                                                                           */
/*  Return an array containing start times for each trial in the index.      */
/*  Each start time is determined by "schan" and "k".                        */
/*                                                                           */
/*****************************************************************************/
float ndata_trial_get_start_sec_schan(t,schan,k)
     struct ndtrial_struct *t;
     char *schan;
     int k;
{
  int j;
  int rtype,tt,flag;
  float tsec;

  j = ndata_trial_get_record_index(t,schan);
  if (j < 0)
    exit_error("NDATA_TRIAL_GET_START_SEC_SCHAN","Channel not found");

  rtype = t->r[j].rtype;
  if (rtype == 0){ /*** This is a point channel, use the kth point ***/
    tsec = (float)(t->r[j].p[k] + t->r[j].t0)/t->r[j].sampling;
  }else if (rtype == 3){
    flag = ndata_record_get_nth_occur_evco_time(&(t->r[j]),k,1,&tt);
    if (flag < 0)
      exit_error("NDATA_TRIAL_GET_START_SEC_SCHAN","ecode not found");
    tsec = (float)(tt + t->r[j].t0)/t->r[j].sampling;
  }else{
    exit_error("NDATA_TRIAL_GET_START_SEC_SCHAN","Record type error");
    tsec = 0.0;
  }

  return tsec;
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_TRIAL_GETP_ADATA_CHAN                       */
/*                                                                           */
/*  Return a pointer to the analog data for the named channel.               */
/*                                                                           */
/*****************************************************************************/
float *ndata_trial_getp_adata_chan(t,chan,rn)
     struct ndtrial_struct *t;
     char *chan;
     int *rn;
{
  int k;

  k = ndata_trial_get_record_index(t,chan);
  if (k < 0)
    exit_error("NDATA_TRIAL_GETP_ADATA_CHAN","Cannot find channel");
  if (t->r[k].rtype != 1)
    return NULL; // This is not a float data record
  else{
    *rn = t->r[k].tn;
    return t->r[k].x;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_TRIAL_GETP_POINT_CHAN                       */
/*                                                                           */
/*  Return a pointer to the point (typically spike) data for "chan" of the   */
/*  trial pointed to by "t".  Also, return the number of points.             */
/*                                                                           */
/*****************************************************************************/
int *ndata_trial_getp_point_chan(t,chan,rn)
     struct ndtrial_struct *t;
     char *chan;
     int *rn;
{
  int k;

  k = ndata_trial_get_record_index(t,chan);
  //if (t->r[k].rtype != 0){  // WYETH CHANGE Jan 5, 2017
  if (t->r[k].rtype == 1){  // All other record types have point data
    *rn = -1;
    return NULL; // This is not a point record
  }else{
    *rn = t->r[k].n;
    return t->r[k].p;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                       NDATA_TRIAL_GETP_EVCO_LIST_CHAN                     */
/*                                                                           */
/*  Return a pointer to the event code list for the named channel of the     */
/*  trial pointed to by "t".  Also, return the number of events.             */
/*                                                                           */
/*****************************************************************************/
int *ndata_trial_getp_evco_list_chan(t,chan,rn)
     struct ndtrial_struct *t;
     char *chan;
     int *rn;
{
  int k;

  k = ndata_trial_get_record_index(t,chan);
  if (t->r[k].rtype != 3)
    return NULL; // This is not an event code record
  else{
    *rn = t->r[k].n;
    return t->r[k].ec;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_TRIAL_CHAN_EVCO_TEST                        */
/*                                                                           */
/*  Return 1 if the channel exists and contains the 'ecode', 0 otherwise.    */
/*                                                                           */
/*****************************************************************************/
int ndata_trial_chan_evco_test(t,chan,ecode)
     struct ndtrial_struct *t;
     char *chan;
     int ecode;
{
  int i;
  int *elist,n,flag;

  flag = 0;  // Assume it will not be found

  elist = ndata_trial_getp_evco_list_chan(t,chan,&n); //  Get event code list

  if (elist != NULL){
    i = get_index_search_iarray(elist,n,ecode);  // Find code in list
    if (i >= 0)
      flag = 1;
  }

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                        NDATA_TRIAL_GET_N_POINTS_CHAN                      */
/*                                                                           */
/*  Return the number of points on the named channel of this trial.          */
/*                                                                           */
/*  Error return values:                                                     */
/*    -1 Channel name not found.                                             */
/*    -2 Channel found, but not of point type.                               */
/*                                                                           */
/*****************************************************************************/
int ndata_trial_get_n_points_chan(t,chan)
     struct ndtrial_struct *t;
     char *chan;
{
  int k;

  k = ndata_trial_get_record_index(t,chan);
  if (k < 0)
    return -1;
  else if (t->r[k].rtype != 0)
    return -2;
  else
    return t->r[k].n; /* Number of points in this record. */
}
/**************************************-**************************************/
/*                                                                           */
/*                       NDATA_TRIAL_GET_NTH_POINT_CHAN                      */
/*                                                                           */
/*  Get the value of the nth point in the specified record name in the       */
/*  specified trial number.  If the trial doesn't have an nth point, return  */
/*  0, else return 1.                                                        */
/*                                                                           */
/*****************************************************************************/
int ndata_trial_get_nth_point_chan(t,chan,n,rp)
     struct ndtrial_struct *t;
     char *chan;
     int n,*rp;
{
  int k;
  int np;

  if (n < 0)
    exit_error("NDATA_TRIAL_GET_NTH_POINT_CHAN","n < 0");

  k = ndata_trial_get_record_index(t,chan);
  if (k < 0)
    exit_error("NDATA_TRIAL_GET_NTH_POINT_CHAN","Record not found");
  else if (t->r[k].rtype != 0)
    exit_error("NDATA_TRIAL_GET_NTH_POINT_CHAN","Record type is not point");
  else{
    np = t->r[k].n; /* Number of points in this record. */
    if (n >= np){
      *rp = 0;
      return 0;
    }else{
      *rp = t->r[k].p[n];
      return 1;
    }
  }
  return 0; /* This should never happen, here for lint. */
}
/**************************************-**************************************/
/*                                                                           */
/*                       NDATA_TRIAL_GET_EVENT_TIME_CHAN                     */
/*                                                                           */
/*  This is meant to be a general time-finding routine for events and        */
/*  point data.  Time is returned in sampling units.                         */
/*                                                                           */
/*****************************************************************************/
int ndata_trial_get_event_time_chan(t,chan,ecode,occur,rtime)
     struct ndtrial_struct *t;
     char *chan;
     int ecode,occur,*rtime;
{
  int i,k;
  int flag;

  i = ndata_trial_get_record_index(t,chan);
  if (i < 0)
    exit_error("NDATA_TRIAL_GET_EVENT_TIME_CHAN","Channel not found");
  if (t->r[i].rtype == 0){ // This is a point channel, use the kth point
    if (occur < t->r[i].n){
      k = t->r[i].p[occur-1];
    }else
      exit_error("NDATA_TRIAL_GET_EVENT_TIME_CHAN","Occurrence not found");
  }else if (t->r[i].rtype == 3){ // Event channel
    flag = ndata_record_get_nth_occur_evco_time(&(t->r[i]),ecode,occur,&k);
    if (flag==0)
      exit_error("NDATA_TRIAL_GET_EVENT_TIME_CHAN","Event code not found");
  }else
    exit_error("NDATA_TRIAL_GET_EVENT_TIME_CHAN","Record type error");

  *rtime = k;
  return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_FREE_RECORD_CONTENTS                        */
/*                                                                           */
/*****************************************************************************/
void ndata_free_record_contents(r)
     struct ndrec_struct *r;
{
  int rt;

  myfree(r->name);
  
  rt = r->rtype;
  if ((rt==0)||(rt==2)||(rt==3))
    myfree(r->p);
  if ((rt==1)||(rt==2))
    myfree(r->x);
  if (rt==3)
    myfree(r->ec);
  if ((rt<0)||(rt>3))
    exit_error("NDATA_FREE_RECORD","Unknown record type");
}
/**************************************-**************************************/
/*                                                                           */
/*                           NDATA_RECORD_CREATE_TYPE_1                      */
/*                                                                           */
/*****************************************************************************/
struct ndrec_struct *ndata_record_create_type_1(name,sampling,t0,tn,data)
     char name[];
     float sampling;
     int t0,tn;
     float *data;
{
  struct ndrec_struct *r;

  r = (struct ndrec_struct *)myalloc(sizeof(struct ndrec_struct));
  
  r->rtype = 1; // float
  r->name = strdup(name);
  r->rcode = 0;
  r->sampling = sampling;
  r->t0 = t0;
  r->tn = r->n = tn;  // Should pick one or the other
  r->x = copy_farray(data,tn);

  // Unused
  r->p = NULL;
  r->ec = NULL;
  r->ect = NULL;

  return r;
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_RECORD_COPY_CONTENTS                        */
/*                                                                           */
/*  Copy the contents of record "r2" to record "r1".  Presumably, "r1" is    */
/*  empty.                                                                   */
/*                                                                           */
/*****************************************************************************/
void ndata_record_copy_contents(r1,r2)
     struct ndrec_struct *r1,*r2;
{
  int rt;

  rt = r2->rtype;
  r1->rtype = rt;
  r1->name = strdup(r2->name);
  r1->rcode = r2->rcode;
  r1->sampling = r2->sampling;
  r1->t0 = r2->t0;
  r1->tn = r2->tn;
  r1->n = r2->n;

  if ((rt==0)||(rt==2)||(rt==3))
    r1->p = copy_iarray(r2->p,r2->n);
  if ((rt==1)||(rt==2))
    r1->x = copy_farray(r2->x,r2->tn);
  if (rt==3)
    r1->ec = copy_iarray(r2->ec,r2->n);
  if ((rt<0)||(rt>3))
    exit_error("NDATA_RECORD_COPY_CONTENTS","Unknown record type");
  r1->ect = r2->ect;
}
/**************************************-**************************************/
/*                                                                           */
/*                             NDATA_TRIAL_CREATE                            */
/*                                                                           */
/*****************************************************************************/
struct ndtrial_struct *ndata_trial_create(nrec)
     int nrec;
{
  struct ndtrial_struct *tr;

  tr = (struct ndtrial_struct *)myalloc(sizeof(struct ndtrial_struct));
  tr->tcode = 0;
  tr->tref = 0;
  tr->nparam = 0;
  tr->pname = NULL;
  tr->pval = NULL;
  tr->nrec = nrec;
  tr->next = NULL;
  tr->prev = NULL;

  if (nrec > 0){
    tr->r = (struct ndrec_struct *)myalloc(nrec*sizeof(struct ndrec_struct));
  }else
    tr->r = NULL;

  return tr;
}
/**************************************-**************************************/
/*                                                                           */
/*                          NDATA_TRIAL_INSERT_RECORD                        */
/*                                                                           */
/*  Insert the record 'r' into the trial 't'.                                */
/*                                                                           */
/*****************************************************************************/
void ndata_trial_insert_record(t,r)
     struct ndtrial_struct *t;
     struct ndrec_struct *r;
{
  int i;
  int n;
  struct ndrec_struct *rlist;

  /* Make storage for new record list (too bad no pointers!) */
  n = t->nrec + 1;
  rlist = (struct ndrec_struct *)myalloc(n*sizeof(struct ndrec_struct));

  /* Copy existing records, then add new record at end */
  for(i=0;i<(n-1);i++)
    ndata_record_copy_contents(&(rlist[i]),&(t->r[i]));
  i = n-1;
  ndata_record_copy_contents(&(rlist[i]),r);

  /* Free old records */
  for(i=0;i<t->nrec;i++)
    ndata_free_record_contents(&(t->r[i]));
  myfree(t->r);

  t->nrec = n;
  t->r = rlist;
}
/**************************************-**************************************/
/*                                                                           */
/*                       NDATA_TRIAL_DELETE_RECORD_FLAG                      */
/*                                                                           */
/*  Delete the records in 't' that have 'flag' values != 1.                  */
/*                                                                           */
/*****************************************************************************/
void ndata_trial_delete_record_flag(t,flag)
     struct ndtrial_struct *t;
     int *flag;  /* [t->nrec] 0-delete, 1-keep */
{
  int i,k;
  int n;
  struct ndrec_struct *rlist;
  
  /* Make storage for new record list */
  n = t->nrec;
  for(i=0;i<t->nrec;i++)
    if (flag[i] != 1)
      n -= 1;
  if (n<=0)
    exit_error("NDATA_TRIAL_DETELE_RECORD_FLAG","No records would remain");
  rlist = (struct ndrec_struct *)myalloc(n*sizeof(struct ndrec_struct));

  /* Copy existing records, then add new record at end */
  k = 0;
  for(i=0;i<t->nrec;i++){
    if (flag[i] == 1){
      ndata_record_copy_contents(&(rlist[k]),&(t->r[i]));
      k += 1;
    }
  }

  /* Free old records */
  for(i=0;i<t->nrec;i++)
    ndata_free_record_contents(&(t->r[i]));
  myfree(t->r);

  t->nrec = n;
  t->r = rlist;
}
/**************************************-**************************************/
/*                                                                           */
/*                           NDATA_TRIAL_COPY_CONTENTS                       */
/*                                                                           */
/*  Copy the contents of trial "t2" to trial "t1".  Presumably, "t1" is      */
/*  empty.                                                                   */
/*                                                                           */
/*****************************************************************************/
void ndata_trial_copy_contents(t1,t2)
     struct ndtrial_struct *t1,*t2;
{
  int i;

  t1->tcode = t2->tcode;
  t1->tref = t2->tref;
  t1->nparam = t2->nparam;
  t1->pname = copy_2d_carray(t2->pname,t2->nparam);
  t1->pval = copy_2d_carray(t2->pval,t2->nparam);
  t1->nrec = t2->nrec;
  t1->r = (struct ndrec_struct *)myalloc(t1->nrec*sizeof(struct ndrec_struct));
  for(i=0;i<t1->nrec;i++)
    ndata_record_copy_contents(&(t1->r[i]),&(t2->r[i]));
}
/**************************************-**************************************/
/*                                                                           */
/*                                 NDATA_COPY                                */
/*                                                                           */
/*  Copy the contents of 'nd' to a new nd struct.                            */
/*                                                                           */
/*****************************************************************************/
struct ndata_struct *ndata_copy(nd)
     struct ndata_struct *nd;
{
  int i;
  struct ndata_struct *nd2;

  nd2 = (struct ndata_struct *)myalloc(sizeof(struct ndata_struct));
  nd2->class = strdup(nd->class);

  nd2->ntrial = nd->ntrial;
  nd2->nconst = nd->nconst;
  nd2->nvar   = nd->nvar;
  nd2->cname = copy_2d_carray(nd->cname,nd->nconst);
  nd2->cval  = copy_2d_carray(nd->cval,nd->nconst);
  nd2->vname = copy_2d_carray(nd->vname,nd->nvar);

  nd2->ctype = (char *)myalloc(nd->nconst*sizeof(char));
  for(i=0;i<nd->nconst;i++)
    nd2->ctype[i] = nd->ctype[i];

  nd2->vtype = (char *)myalloc(nd->nvar*sizeof(char));
  for(i=0;i<nd->nvar;i++)
    nd2->vtype[i] = nd->vtype[i];

  nd2->ntable = nd->ntable;
  if (nd2->ntable > 0){
    nd2->table = (struct ect_struct *)myalloc(nd->ntable*
					      sizeof(struct ect_struct));
    for(i=0;i<nd2->ntable;i++){
      nd2->table[i].tname = strdup(nd->table[i].tname);
      nd2->table[i].tnum  = nd->table[i].tnum;
      nd2->table[i].n     = nd->table[i].n;
      nd2->table[i].name  = copy_2d_carray(nd->table[i].name,nd->table[i].n);
      nd2->table[i].num   = copy_iarray(nd->table[i].num,nd->table[i].n);
    }
  }else{
    nd2->table = NULL;
  }

  nd2->t = (struct ndtrial_struct *)myalloc(nd->ntrial*
					    sizeof(struct ndtrial_struct));
  for(i=0;i<nd->ntrial;i++)
    ndata_trial_copy_contents(&(nd2->t[i]),&(nd->t[i]));

  return nd2;
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_INSERT_EVENT_RECORDS                        */
/*                                                                           */
/*  Insert the event records into the existing nd structure.                 */
/*                                                                           */
/*****************************************************************************/
void ndata_insert_event_records(nd,s,cnt,nchan,duration,sampling,chname)
  struct ndata_struct *nd;
  int ***s,**cnt,nchan,duration;
  float sampling;
  char **chname;
{
  int i,j,k;
  struct ndrec_struct *tr;
  int totchan,nold,t0;

  printf("  NDATA_INSERT_EVENT_RECORDS\n");

  for(i=0;i<nd->ntrial;i++){
    k = ndata_trial_get_record_index(&(nd->t[i]),"sync0");
    if (k < 0)
      exit_error("NDATA_INSERT_EVENT_RECORDS","Cannot find record name sync0");
    /* Use sync channel t0 value (converted for sampling) for new records. */
    t0 = (int)(0.5 + (float)nd->t[i].r[k].t0 / nd->t[i].r[k].sampling
	       * sampling);
    nold = nd->t[i].nrec; /* Original number of channels. */
    totchan = nchan + nold; /* New plus original channels. */
    tr = nd->t[i].r; /* Point to old record list. */
    nd->t[i].r = (struct ndrec_struct *)myalloc(totchan*
						sizeof(struct ndrec_struct));
    nd->t[i].nrec = totchan;
    for(j=0;j<nold;j++){
      ndata_record_copy_contents(&(nd->t[i].r[j]),&(tr[j]));
      ndata_free_record_contents(&(tr[j]));
    }

    for(j=0;j<nchan;j++){
      nd->t[i].r[nold+j].rtype = 0;
      nd->t[i].r[nold+j].name = strdup(chname[j]);
      nd->t[i].r[nold+j].rcode = 0; /* Unused. */
      nd->t[i].r[nold+j].sampling = sampling;
      nd->t[i].r[nold+j].t0 = t0;
      nd->t[i].r[nold+j].tn = duration;
      nd->t[i].r[nold+j].n = cnt[i][j];
      nd->t[i].r[nold+j].p = copy_iarray(s[i][j],cnt[i][j]);
      nd->t[i].r[nold+j].x = NULL;
      nd->t[i].r[nold+j].ec = NULL;
      nd->t[i].r[nold+j].ect = NULL;
    }
    myfree(tr); /* Free the list of pointers to old records. */
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_TIMETRIAL_GET_MERGE_I                       */
/*                                                                           */
/*  Return the merged spike train for channel 'i', assuming this file        */
/*  has timetrial format.                                                    */
/*                                                                           */
/*****************************************************************************/
int *ndata_timetrial_get_merge_i(nd,ci,rn)
     struct ndata_struct *nd;
     int ci;      // Channel index
     int *rn;     // Return number of spikes in merged trial
{
  int i,j,k;
  int n,sn,*s,tref,toff;
  struct ndrec_struct *tr;

  n = nd->ntrial;
  sn = 0;
  for(i=0;i<n;i++){
    tr = &(nd->t[i].r[ci]);  // pointer to record
    if (tr->rtype != 0)
      exit_error("NDATA_TIMETRIAL_GET_MERGE_I","Record type not zero");
    sn += tr->n;            // Number of spikes
  }

  printf("sn = %d\n",sn);
  s = (int *)myalloc(sn*sizeof(int));

  k = 0;
  for(i=0;i<n;i++){
    tref = nd->t[i].tref;

    tr = &(nd->t[i].r[ci]);  // pointer to record

    toff = my_rint(tref * tr->sampling/1000.0);  // ASSUMING 'tref' is msec

    for(j=0;j<tr->n;j++){
      s[k] = toff + tr->p[j];       // Number of spikes
      k += 1;
    }
  }

  *rn = sn;
  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                           NDATA_LOOKUP_ECODE_NAME                         */
/*                                                                           */
/*  Return the string value for the ecode in the ecode table.                */
/*                                                                           */
/*****************************************************************************/
void ndata_lookup_ecode_name(t,ecode,rname)
     struct ect_struct *t;
     int ecode;
     char **rname;
{
  int i;

  if (t != NULL)
    i = get_index_search_iarray(t->num,t->n,ecode);
  else
    i = -1;

  //printf("HERE i = %d\n",i);

  if (i != -1)
    *rname = strdup(t->name[i]);
  else
    *rname = NULL;
}
/**************************************-**************************************/
/*                                                                           */
/*                           NDATA_LOOKUP_ECODE_NUM                          */
/*                                                                           */
/*  Return the ecode value given the name.                                   */
/*                                                                           */
/*****************************************************************************/
int ndata_lookup_ecode_num(t,ev_name,rnum)
     struct ect_struct *t;
     char *ev_name;
     int *rnum;
{
  int i;

  i = search_2d_carray(t->name,ev_name,t->n);
  if (i >= 0){
    *rnum = t->num[i];
    i = 1;  // Success
  }

  return i;
}
/**************************************-**************************************/
/*                                                                           */
/*                             NDATA_ECODE_GET_INT                           */
/*                                                                           */
/*  Return the ecode value an integer, where 'estr' is either the string     */
/*  of the int value or the string of the ecode name.                        */
/*                                                                           */
/*****************************************************************************/
int ndata_ecode_get_int(tab,estr)
     struct ect_struct *tab;   // Ecode table
     char *estr;               // This could be the name or integer
{
  int ecode,flag;

  if (is_int_string(estr) == 1)
    ecode = atoi(estr);
  else{
    flag = ndata_lookup_ecode_num(tab,estr,&ecode);
    if (flag == 0){
      printf("  Event:  %s\n",estr);
      exit_error("NDATA_ECODE_GET_INT","Event name not found in table");
    }
  }
  return ecode;
}
/**************************************-**************************************/
/*                                                                           */
/*                        NDU_TRIAL_CHAN_ECODE_STR_DIFF                      */
/*                                                                           */
/*  Return the time difference for the event stringss on the channel.        */
/*                                                                           */
/*****************************************************************************/
int ndu_trial_chan_ecode_str_diff(tr,chan,e1,e2,rdiff)
     struct ndtrial_struct *tr;
     char *chan;
     char *e1,*e2;  // Can be either a number or name (e.g., "5" or "Start")
     int  *rdiff;
{
  int flag,evt1,evt2,t1,t2;
  struct ndrec_struct *rp;

  rp = ndata_trial_get_record_ptr(tr,chan);

  //  Convert event codes to integers
  evt1 = ndata_ecode_get_int(rp->ect,e1);  // Will exit if not found
  evt2 = ndata_ecode_get_int(rp->ect,e2);

  //  Set 'flag' to 1 if both codes occur on this trial/chan
  flag = ndata_trial_chan_evco_test(tr,chan,evt1);
  if (flag == 1)  // If first code was found
    flag = ndata_trial_chan_evco_test(tr,chan,evt2);

  if (flag){  // Both codes exist
    flag = ndata_record_get_nth_occur_evco_time(rp,evt1,1,&t1);
    flag = ndata_record_get_nth_occur_evco_time(rp,evt2,1,&t2);

    *rdiff = t2 - t1;
    flag = 1;
  }else{
    *rdiff = 0.0;
    flag = -1;  // Failure
  }

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                            NDATA_GET_DEFAULT_CHAN                         */
/*                                                                           */
/*  If a "chan" is unspecified for an analysis, then the default will be     */
/*  assumed: the first record name in the first trial.                       */
/*                                                                           */
/*****************************************************************************/
void ndata_get_default_chan(nd,rname)
     struct ndata_struct *nd;
     char **rname;
{
  *rname = NULL;
  if (nd->ntrial > 0)
    if (nd->t[0].nrec > 0)
      *rname = strdup(nd->t[0].r[0].name);
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_GET_FIRST_UNIT_CHAN                         */
/*                                                                           */
/*  Return the name of the first channel whose name begins with "unit".      */
/*  Look in the *first* trial.                                               */
/*                                                                           */
/*****************************************************************************/
void ndata_get_first_unit_chan(nd,rname)
     struct ndata_struct *nd;
     char **rname;
{
  int i;
  char *c;

  *rname = NULL;
  if (nd->ntrial > 0)
    for(i=0;i<nd->t[0].nrec;i++){
      c = nd->t[0].r[i].name;
      if ((c[0]=='u')&&(c[1]=='n')&&(c[2]=='i')&&(c[3]=='t')&&(*rname==NULL))
	*rname = strdup(nd->t[0].r[i].name);
    }
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_GET_FIRST_PREFIX_CHAN                       */
/*                                                                           */
/*  Return the name of the first channel whose name begins with "prefix".    */
/*  Look in the *first* trial.                                               */
/*                                                                           */
/*****************************************************************************/
void ndata_get_first_prefix_chan(nd,prefix,rname)
     struct ndata_struct *nd;
     char *prefix;
     char **rname;
{
  int i;
  char *c;

  *rname = NULL;
  if (nd->ntrial > 0){
    for(i=0;i<nd->t[0].nrec;i++){
      c = nd->t[0].r[i].name;
      if (compare_prefix_string_order(prefix,c)){
	*rname = strdup(c);
	return;
      }
    }
  }else{
    exit_error("NDATA_GET_FIRST_PREFIX_CHAN","There are no trials");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_GET_FIRST_EVENT_CHAN                        */
/*                                                                           */
/*  Return the name of the first channel whose name begins with "event".     */
/*  Look in the *first* trial.                                               */
/*                                                                           */
/*****************************************************************************/
void ndata_get_first_event_chan(nd,rname)
     struct ndata_struct *nd;
     char **rname;
{
  int i;
  char *c;

  *rname = NULL;
  if (nd->ntrial > 0)
    for(i=0;i<nd->t[0].nrec;i++){
      c = nd->t[0].r[i].name;
      if ((c[0]=='e')&&(c[1]=='v')&&(c[2]=='e')&&(c[3]=='n')&&(c[4]=='t')&&
	  (*rname==NULL))
	*rname = strdup(nd->t[0].r[i].name);
    }
}
/**************************************-**************************************/
/*                                                                           */
/*                            NDATA_GET_RECORD_INDEX                         */
/*                                                                           */
/*  An "ndrx_struct" is intended to point to storage defined elsewhere,      */
/*  and the records pointed to by the "r" pointers should not be freed       */
/*  with the structure.                                                      */
/*                                                                           */
/*****************************************************************************/
void ndata_get_record_index(nd,name,rrx)
     struct ndata_struct *nd;
     char name[];
     struct ndrx_struct **rrx;
{
  int i,j,k;
  struct ndrx_struct *rx;
  int count;

  // Count the number of times the record appears in "nd"
  count = 0;
  for(i=0;i<nd->ntrial;i++){
    j = ndata_trial_get_record_index(&(nd->t[i]),name);
    if (j >= 0)
      count += 1;
  }
  if (count < nd->ntrial)
    printf("    *** NDATA_GET_RECORD_INDEX:  Note, %s in %d of %d trials.\n",
	   name,count,nd->ntrial);
  
  rx = (struct ndrx_struct *)myalloc(sizeof(struct ndrx_struct));
  rx->name = strdup(name);
  rx->n = count;
  rx->tnum = (int *)myalloc(rx->n*sizeof(int));
  rx->r = (struct ndrec_struct **)myalloc(rx->n*sizeof(struct ndrec_struct *));

  k = 0;
  for(i=0;i<nd->ntrial;i++){
    j = ndata_trial_get_record_index(&(nd->t[i]),name);
    if (j >= 0){
      rx->tnum[k] = i;
      rx->r[k] = &(nd->t[i].r[j]);
      k += 1;
    }
  }
  *rrx = rx;
}
/**************************************-**************************************/
/*                                                                           */
/*                               NDATA_FREE_NDRX                             */
/*                                                                           */
/*  An "ndrx_struct" is intended to point to storage defined else where,     */
/*  and the records pointed to by the "r" pointers should not be freed       */
/*  with the structure.                                                      */
/*                                                                           */
/*****************************************************************************/
void ndata_free_ndrx(rx)
     struct ndrx_struct *rx;
{
  myfree(rx->name);
  myfree(rx->tnum);
  myfree(rx->r);
  myfree(rx);
}
/**************************************-**************************************/
/*                                                                           */
/*                       NDATA_GET_ADATA_RECORD_POINTER                      */
/*                                                                           */
/*  Return an array of pointers to the analog data for "name" record.        */
/*                                                                           */
/*****************************************************************************/
void ndata_get_adata_record_pointer(nd,name,rdata,rcnt,rn)
     struct ndata_struct *nd;
     char name[];
     float ***rdata;
     int **rcnt,*rn;
{
  int i;
  struct ndrx_struct *rx;
  int *cnt,n;
  float **data;

  ndata_get_record_index(nd,name,&rx);

  n = rx->n;
  if (n != nd->ntrial)
    exit_error("NDATA_GET_ADATA_RECORD_POINTER","Trial count violation");

  data = (float **)myalloc(n*sizeof(float *));
  cnt = (int *)myalloc(n*sizeof(int));
  for(i=0;i<n;i++){
    if (rx->tnum[i] != i){
      printf("  *** i=%d  rx->tnum[i]=%d\n",i,rx->tnum[i]);
      exit_error("NDATA_GET_ADATA_RECORD_POINTER","Trial order violation");
    }
    data[i] = rx->r[i]->x;
    cnt[i] = rx->r[i]->tn;
  }
  ndata_free_ndrx(rx);
  
  *rdata = data; *rcnt = cnt; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                        NDATA_GET_SAMPLING_RECORD_INDEX                    */
/*                                                                           */
/*  Values returned:                                                         */
/*    -1.0  the sampling varied across trials for this rec.                  */
/*     0.0  there were no records.                                           */
/*  1000.0  the sampling was always 1000.0 for this record.                  */
/*                                                                           */
/*****************************************************************************/
float ndata_get_sampling_record_index(nd,ndrx)
     struct ndata_struct *nd;
     struct ndrx_struct *ndrx;
{
  int i;
  float sampling;

  if (ndrx->n <= 0){
    return 0.0;
  }else{
    sampling = ndrx->r[0]->sampling;
    for(i=1;i<ndrx->n;i++)
      if (sampling != ndrx->r[i]->sampling)
	return -1.0;
    return sampling;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          NDATA_GET_TYPE_RECORD_INDEX                      */
/*                                                                           */
/*  Values returned:                                                         */
/*    -2 there were no records.                                              */
/*    -1 the type varied across trials for this record name.                 */
/*     n the type was 0,1,2,3, ...                                           */
/*                                                                           */
/*****************************************************************************/
int ndata_get_type_record_index(nd,ndrx)
     struct ndata_struct *nd;
     struct ndrx_struct *ndrx;
{
  int i;
  int rtype;

  if (ndrx->n <= 0){
    return -2;
  }else{
    rtype = ndrx->r[0]->rtype;
    for(i=1;i<ndrx->n;i++)
      if (rtype != ndrx->r[i]->rtype)
	return -1;
    return rtype;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                       NDATA_GET_SAMPLING_RECORD_NAME                      */
/*                                                                           */
/*  Values returned:                                                         */
/*    -1.0  the sampling varied across trials for this rec.                  */
/*     0.0  there were no records.                                           */
/*  1000.0  the sampling was always 1000.0 for this record.                  */
/*                                                                           */
/*****************************************************************************/
float ndata_get_sampling_record_name(nd,name)
     struct ndata_struct *nd;
     char name[];
{
  struct ndrx_struct *rx;
  float sampling;

  ndata_get_record_index(nd,name,&rx);
  sampling = ndata_get_sampling_record_index(nd,rx);
  ndata_free_ndrx(rx);

  return sampling;
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_GET_TYPE_RECORD_NAME                        */
/*                                                                           */
/*  Values returned:                                                         */
/*    -2 there were no records.                                              */
/*    -1 the type varied across trials for this record name.                 */
/*     n the type was 0,1,2,3, ...                                           */
/*                                                                           */
/*****************************************************************************/
int ndata_get_type_record_name(nd,name)
     struct ndata_struct *nd;
     char name[];
{
  struct ndrx_struct *rx;
  int type;

  ndata_get_record_index(nd,name,&rx);
  type = ndata_get_type_record_index(nd,rx);
  ndata_free_ndrx(rx);

  if (type == -2){
    printf("  *** record name = %s\n",name);
    exit_error("NDATA_GET_TYPE_RECORD_NAME","No records in index");
  }else if (type == -1){
    printf("  *** record name = %s\n",name);
    exit_error("NDATA_GET_TYPE_RECORD_NAME","Type varies across trials");
  }

  return type;
}
/**************************************-**************************************/
/*                                                                           */
/*                            NDATA_GET_CONST_TYPE                           */
/*                                                                           */
/*  Returns 1 if param "cname" found in list of constant parameters, other-  */
/*  wise return 0.  Return the character type of this parameter.             */
/*                                                                           */
/*****************************************************************************/
int ndata_get_const_type(nd,cname,rctype)
     struct ndata_struct *nd;
     char cname[];
     char *rctype;
{
  int i;

  i = search_2d_carray(nd->cname,cname,nd->nconst);
  if (i == -1){
    *rctype = '\0';
    return 0;
  }else{
    *rctype = nd->ctype[i];
    return 1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            NDATA_GET_VAR_TYPE                             */
/*                                                                           */
/*  Returns 1 if param "vname" found in list of parameters, 0 other-         */
/*  wise.  Return the character type of this parameter.                      */
/*                                                                           */
/*****************************************************************************/
int ndata_get_var_type(nd,vname,rvtype)
     struct ndata_struct *nd;
     char vname[];
     char *rvtype;
{
  int i;

  i = search_2d_carray(nd->vname,vname,nd->nvar);
  if (i == -1){
    *rvtype = '\0';
    return 0;
  }else{
    *rvtype = nd->vtype[i];
    return 1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                        NDATA_GET_VAR_OR_CONST_TYPE                        */
/*                                                                           */
/*  Get the parameter data type.  Check first for a variable parameter,      */
/*  then check the constant param list.                                      */
/*                                                                           */
/*****************************************************************************/
int ndata_get_var_or_const_type(nd,name,rtype)
     struct ndata_struct *nd;
     char name[];
     char *rtype;
{
  int flag;

  flag = ndata_get_var_type(nd,name,rtype);
  if (!flag)
    flag = ndata_get_const_type(nd,name,rtype);
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                            NDATA_TEST_NDA_PARAM                           */
/*                                                                           */
/*  Returns 1 if param "pname" found in list of parameters, 0 other-         */
/*  wise.  Does *not* malloc storage.                                        */
/*                                                                           */
/*****************************************************************************/
int ndata_test_nda_param(nda,pname)
     struct nda_struct *nda;
     char pname[];
{
  int i;

  i = search_2d_carray(nda->pname,pname,nda->nparam);
  if (i == -1)
    return 0;
  else
    return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                            NDATA_GET_NDA_PARAM                            */
/*                                                                           */
/*  Returns 1 if param "pname" found in list of parameters, 0 otherwise.     */
/*  Return the character value of this parameter.                            */
/*                                                                           */
/*****************************************************************************/
int ndata_get_nda_param(nda,pname,rpval)
     struct nda_struct *nda;
     char pname[];
     char **rpval;
{
  int i;

  i = search_2d_carray(nda->pname,pname,nda->nparam);
  if (i == -1){
    *rpval = (char *)NULL;
    return 0;
  }else{
    *rpval = strdup(nda->pval[i]);
    return 1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            NDATA_SET_NDA_PARAM                            */
/*                                                                           */
/*  Returns 1 if param "pname" found in list of parameters, 0 otherwise.     */
/*  Set a new value for this param.                                          */
/*                                                                           */
/*  NOTE: this subroutine name uses 'set' but in a different way than        */
/*  the other 'ndata_set_...' routines that actually 'get' values.           */
/*                                                                           */
/*****************************************************************************/
int ndata_set_nda_param(nda,pname,pval)
     struct nda_struct *nda;
     char pname[];
     char *pval;
{
  int i;

  i = search_2d_carray(nda->pname,pname,nda->nparam);
  if (i == -1){
    return 0;
  }else{
    myfree(nda->pval[i]);         // Free old storage
    nda->pval[i] = strdup(pval);  // Point to new value
    return 1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                       NDATA_GET_NDA_PARAM_CHAR_OR_EXIT                    */
/*                                                                           */
/*  Returns 1 if param "pname" found in list of parameters, 0 otherwise.     */
/*  Return the character value of this parameter.                            */
/*                                                                           */
/*****************************************************************************/
char *ndata_get_nda_param_char_or_exit(nda,pname)
     struct nda_struct *nda;
     char pname[];
{
  int i;

  i = search_2d_carray(nda->pname,pname,nda->nparam);
  if (i == -1){
    printf("  %s\n",pname);
    exit_error("NDATA_GET_NDA_PARAM_CHAR_OR_EXIT","Param not found");
  }

  return strdup(nda->pval[i]);
}
/**************************************-**************************************/
/*                                                                           */
/*                        NDATA_GET_NDA_PARAM_INT_OR_EXIT                    */
/*                                                                           */
/*****************************************************************************/
int ndata_get_nda_param_int_or_exit(nda,pname)
     struct nda_struct *nda;
     char pname[];
{
  int flag;
  int k;
  char *tstr;
  
  flag = ndata_get_nda_param(nda,pname,&tstr);
  if (flag){
    k = atoi(tstr);
    myfree(tstr);
  }else{
    printf("  %s\n",pname);
    exit_error("NDATA_GET_NDA_PARAM_INT_OR_EXIT","Param not found");
  }
  
  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                      NDATA_GET_NDA_PARAM_FLOAT_OR_EXIT                    */
/*                                                                           */
/*****************************************************************************/
float ndata_get_nda_param_float_or_exit(nda,pname)
     struct nda_struct *nda;
     char pname[];
{
  int flag;
  float x;
  char *tstr;

  flag = ndata_get_nda_param(nda,pname,&tstr);
  if (flag){
    x = atof(tstr);
    myfree(tstr);
  }else{
    printf("  %s\n",pname);
    exit_error("NDATA_GET_NDA_PARAM_FLOAT_OR_EXIT","Param not found");
  }

  return x;
}
/**************************************-**************************************/
/*                                                                           */
/*                               NDA_GETC_EXIT                               */
/*                                                                           */
/*  Short name for "ndata_get_nda_param_char_or_exit"                        */
/*                                                                           */
/*****************************************************************************/
char *nda_getc_exit(nda,pname)
     struct nda_struct *nda;
     char pname[];
{
  return ndata_get_nda_param_char_or_exit(nda,pname);
}
/**************************************-**************************************/
/*                                                                           */
/*                               NDA_GETI_EXIT                               */
/*                                                                           */
/*  Short name for "ndata_get_nda_param_int_or_exit"                         */
/*                                                                           */
/*****************************************************************************/
int nda_geti_exit(nda,pname)
     struct nda_struct *nda;
     char pname[];
{
  return ndata_get_nda_param_int_or_exit(nda,pname);
}
/**************************************-**************************************/
/*                                                                           */
/*                               NDA_GETF_EXIT                               */
/*                                                                           */
/*  Short name for "ndata_get_nda_param_float_or_exit"                       */
/*                                                                           */
/*****************************************************************************/
float nda_getf_exit(nda,pname)
     struct nda_struct *nda;
     char pname[];
{
  return ndata_get_nda_param_float_or_exit(nda,pname);
}
/**************************************-**************************************/
/*                                                                           */
/*                       NDATA_SET_NDA_PARAM_DEFAULT_INT                     */
/*                                                                           */
/*  If "pname" is found, use the associated integer value, otherwise, use    */
/*  the specified default value.                                             */
/*                                                                           */
/*****************************************************************************/
void ndata_set_nda_param_default_int(nda,pname,rparam,defval)
     struct nda_struct *nda;
     char pname[];
     int *rparam,defval;
{
  int flag;
  char *tstr;

  flag = ndata_get_nda_param(nda,pname,&tstr);
  if (flag){
    *rparam = atoi(tstr);
    myfree(tstr);
  }else{
    *rparam = defval;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                      NDATA_SET_NDA_PARAM_DEFAULT_FLOAT                    */
/*                                                                           */
/*  If "pname" is found, use the associated float value, otherwise, use the  */
/*  specified default value.                                                 */
/*                                                                           */
/*****************************************************************************/
void ndata_set_nda_param_default_float(nda,pname,rparam,defval)
     struct nda_struct *nda;
     char pname[];
     float *rparam,defval;
{
  int flag;
  char *tstr;

  flag = ndata_get_nda_param(nda,pname,&tstr);
  if (flag){
    *rparam = atof(tstr);
    myfree(tstr);
  }else{
    *rparam = defval;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                      NDATA_SET_NDA_PARAM_DEFAULT_CHAR                     */
/*                                                                           */
/*  If "pname" is found, use the associated string value, otherwise, use     */
/*  the specified default value.                                             */
/*                                                                           */
/*****************************************************************************/
void ndata_set_nda_param_default_char(nda,pname,rparam,defval)
     struct nda_struct *nda;
     char pname[];
     char **rparam,*defval;
{
  int flag;

  flag = ndata_get_nda_param(nda,pname,rparam);
  if (flag==0){
    if (defval == NULL)
      *rparam = NULL;
    else
      *rparam = strdup(defval);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                NDA_GETI_DFLT                              */
/*                                                                           */
/*  Short name for 'ndata_set_nda_param_default_int'.                        */
/*                                                                           */
/*****************************************************************************/
int nda_geti_dflt(nda,pname,defval)
     struct nda_struct *nda;
     char pname[];
     int defval;
{
  int flag;
  int ival;
  char *tstr;

  flag = ndata_get_nda_param(nda,pname,&tstr);
  if (flag){
    ival = atoi(tstr);
    myfree(tstr);
  }else{
    ival = defval;
  }

  return ival;
}
/**************************************-**************************************/
/*                                                                           */
/*                                NDA_GETF_DFLT                              */
/*                                                                           */
/*  Short name for 'ndata_set_nda_param_default_float'.                      */
/*                                                                           */
/*****************************************************************************/
float nda_getf_dflt(nda,pname,defval)
     struct nda_struct *nda;
     char pname[];
     float defval;
{
  int flag;
  float fval;
  char *tstr;

  flag = ndata_get_nda_param(nda,pname,&tstr);
  if (flag){
    fval = atof(tstr);
    myfree(tstr);
  }else{
    fval = defval;
  }

  return fval;
}
/**************************************-**************************************/
/*                                                                           */
/*                                NDA_GETC_DFLT                              */
/*                                                                           */
/*  Short name for 'ndata_set_nda_param_default_char'.                       */
/*                                                                           */
/*****************************************************************************/
char *nda_getc_dflt(nda,pname,defval)
     struct nda_struct *nda;
     char pname[];
     char *defval;
{
  int flag;
  char *tstr;

  flag = ndata_get_nda_param(nda,pname,&tstr);
  if (flag==0){
    if (defval == NULL)
      tstr = NULL;
    else
      tstr = strdup(defval);
  }

  return tstr;
}
/**************************************-**************************************/
/*                                                                           */
/*                           NDATA_GET_CONST_PARAM                           */
/*                                                                           */
/*  Returns 1 if param "pname" found, 0 otherwise.  "rstring" returns the    */
/*  string value associated with the parameter.                              */
/*                                                                           */
/*****************************************************************************/
int ndata_get_const_param(nd,pname,rstring)
     struct ndata_struct *nd;
     char pname[];
     char **rstring;
{
  int i;

  i = search_2d_carray(nd->cname,pname,nd->nconst);
  if (i == -1)
    return 0;
  else{
    *rstring = strdup(nd->cval[i]);
    return 1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                        NDATA_GET_CONST_PARAM_INDEX                        */
/*                                                                           */
/*  Returns the index of 'pname' in the const param list.                    */
/*                                                                           */
/*****************************************************************************/
int ndata_get_const_param_index(nd,pname)
     struct ndata_struct *nd;
     char pname[];
{
  int i;

  i = search_2d_carray(nd->cname,pname,nd->nconst);
  return i;
}
/**************************************-**************************************/
/*                                                                           */
/*                          NDATA_GET_CONST_PARAM_INT                        */
/*                                                                           */
/*****************************************************************************/
int ndata_get_const_param_int(nd,pname)
     struct ndata_struct *nd;
     char pname[];
{
  int i;

  i = ndata_get_const_param_index(nd,pname);
  if (i < 0){
    printf("pname = %s\n",pname);
    exit_error("NDATA_GET_CONST_PARAM_INT","Const param not found");
  }
  return atoi(nd->cval[i]);
}
/**************************************-**************************************/
/*                                                                           */
/*                        NDATA_GET_CONST_PARAM_FLOAT                        */
/*                                                                           */
/*****************************************************************************/
float ndata_get_const_param_float(nd,pname)
     struct ndata_struct *nd;
     char pname[];
{
  int i;

  i = ndata_get_const_param_index(nd,pname);
  if (i < 0){
    printf("pname = %s\n",pname);
    exit_error("NDATA_GET_CONST_PARAM_FLOAT","Const param not found");
  }
  return atof(nd->cval[i]);
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_GET_CONST_PARAM_CHAR                        */
/*                                                                           */
/*****************************************************************************/
char *ndata_get_const_param_char(nd,pname)
     struct ndata_struct *nd;
     char pname[];
{
  int i;

  i = ndata_get_const_param_index(nd,pname);
  if (i < 0){
    printf("pname = %s\n",pname);
    exit_error("NDATA_GET_CONST_PARAM_CHAR","Const param not found");
  }
  return strdup(nd->cval[i]);
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_GET_VAR_PARAM_INDEX                         */
/*                                                                           */
/*  Returns the index of 'pname' in the var param list.                      */
/*                                                                           */
/*****************************************************************************/
int ndata_get_var_param_index(nd,pname)
     struct ndata_struct *nd;
     char pname[];
{
  int i;

  i = search_2d_carray(nd->vname,pname,nd->nvar);
  return i;
}
/**************************************-**************************************/
/*                                                                           */
/*                          NDATA_GET_VAR_VALUE_TRIAL                        */
/*                                                                           */
/*  Returns 1 if param "vname" found in trial "k", 0 otherwise.              */
/*  "rstring" returns the string value associated with the parameter         */
/*  from trial "k".                                                          */
/*                                                                           */
/*****************************************************************************/
int ndata_get_var_value_trial(nd,k,vname,rstring)
     struct ndata_struct *nd;
     int k;
     char vname[];
     char **rstring;
{
  int i;

  i = search_2d_carray(nd->t[k].pname,vname,nd->t[k].nparam);
  if (i == -1)
    return 0;
  else{
    *rstring = strdup(nd->t[k].pval[i]);
    return 1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                        NDATA_GET_VAR_VALUE_TRIAL_INT                      */
/*                                                                           */
/*  Returns 1 if param "vname" found in trial "k", 0 otherwise.  "rval"      */
/*  returns the integer value of the parameter.                              */
/*                                                                           */
/*****************************************************************************/
int ndata_get_var_value_trial_int(nd,k,vname,rval)
     struct ndata_struct *nd;
     int k;
     char vname[];
     int *rval;
{
  int flag;
  char *rstring;

  flag = ndata_get_var_value_trial(nd,k,vname,&rstring);
  if (flag){
    *rval = atoi(rstring);
    myfree(rstring);
  }
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                       NDATA_GET_VAR_VALUE_TRIAL_FLOAT                     */
/*                                                                           */
/*  Returns 1 if param "vname" found in trial "k", 0 otherwise.  "rval"      */
/*  returns the float value of the parameter.                                */
/*                                                                           */
/*****************************************************************************/
int ndata_get_var_value_trial_float(nd,k,vname,rval)
     struct ndata_struct *nd;
     int k;
     char vname[];
     float *rval;
{
  int flag;
  char *rstring;

  flag = ndata_get_var_value_trial(nd,k,vname,&rstring);
  if (flag){
    *rval = atof(rstring);
    myfree(rstring);
  }
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                        NDATA_GET_VAR_OR_CONST_PARAM                       */
/*                                                                           */
/*  If "pname" is a varying parameter, get its value from trial "k", if it   */
/*  is not a varying paramters, get its value from the header, if it is not  */
/*  in the header, return 0.                                                 */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Return 1 on success, 0 on failure.                                     */
/*  - 'rstring' returns new storage - must be freed by caller.               */
/*                                                                           */
/*****************************************************************************/
int ndata_get_var_or_const_param(nd,k,pname,rstring)
     struct ndata_struct *nd;
     int k;
     char pname[],**rstring;
{
  int flag;

  flag = ndata_get_var_value_trial(nd,k,pname,rstring); // New storage
  if (flag==0)
    flag = ndata_get_const_param(nd,pname,rstring); // New storage

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                           NDATA_GET_VC_PARAM_INT                          */
/*                                                                           */
/*  Get the integer value of the parameter, var or const.                    */
/*                                                                           */
/*****************************************************************************/
int ndata_get_vc_param_int(nd,k,pname,rval)
     struct ndata_struct *nd;
     int k;
     char pname[];
     int *rval;
{
  int flag;
  char *sval;

  flag =  ndata_get_var_or_const_param(nd,k,pname,&sval);
  if (flag){
    *rval = atoi(sval);
    myfree(sval);
  }
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                        NDATA_GET_VC_PARAM_INT_OR_EXIT                     */
/*                                                                           */
/*  Get the integer value of the parameter, var or const.                    */
/*                                                                           */
/*  SHORT NAME SEE BELOW!                                                    */
/*                                                                           */
/*****************************************************************************/
int ndata_get_vc_param_int_or_exit(nd,k,pname,callername)
     struct ndata_struct *nd;
     int k;
     char pname[],callername[];
{
  int flag,pval;
  char *sval,temp[SLEN];

  flag =  ndata_get_var_or_const_param(nd,k,pname,&sval);
  if (flag){
    pval = atoi(sval);
    myfree(sval);
  }else{
    sprintf(temp,"Param %s not found",pname);
    printf("*** NDATA_GET_VC_PARAM_INT_OR_EXIT\n");
    exit_error(callername,temp);
    pval = 0.0;
  }

  return pval;
}
/**************************************-**************************************/
/*                                                                           */
/*                             NDATA_GET_VC_INT_EXIT                         */
/*                                                                           */
/*  SHORT NAME FOR ABOVE.                                                    */
/*                                                                           */
/*****************************************************************************/
int ndata_get_vc_int_exit(nd,k,pname,callername)
     struct ndata_struct *nd;
     int k;
     char callername[],pname[];
{
  return ndata_get_vc_param_int_or_exit(nd,k,pname,callername);
}
/**************************************-**************************************/
/*                                                                           */
/*                            NDATA_GET_VC_INT_DFLT                          */
/*                                                                           */
/*****************************************************************************/
int ndata_get_vc_int_dflt(nd,k,pname,dval)
     struct ndata_struct *nd;
     int k;
     char pname[];
     int dval;
{
  int ival,flag;

  flag = ndata_get_vc_param_int(nd,k,pname,&ival);
  if (flag == 0)
    ival = dval;

  return ival;
}
/**************************************-**************************************/
/*                                                                           */
/*                          NDATA_GET_VC_PARAM_FLOAT                         */
/*                                                                           */
/*  Get the float value of the parameter, var or const.                      */
/*                                                                           */
/*****************************************************************************/
int ndata_get_vc_param_float(nd,k,pname,rval)
     struct ndata_struct *nd;
     int k;
     char pname[];
     float *rval;
{
  int flag;
  char *sval;

  flag =  ndata_get_var_or_const_param(nd,k,pname,&sval);
  if (flag){
    *rval = atof(sval);
    myfree(sval);
  }
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                      NDATA_GET_VC_PARAM_FLOAT_OR_EXIT                     */
/*                                                                           */
/*  Get the float value of the parameter, var or const.                      */
/*                                                                           */
/*****************************************************************************/
float ndata_get_vc_param_float_or_exit(nd,k,pname,callername)
     struct ndata_struct *nd;
     int k;
     char pname[],callername[];
{
  int flag;
  float pval;
  char *sval,temp[SLEN];

  flag =  ndata_get_var_or_const_param(nd,k,pname,&sval);
  if (flag){
    pval = atof(sval);
    myfree(sval);
  }else{
    sprintf(temp,"Param %s not found",pname);
    printf("*** NDATA_GET_VC_PARAM_FLOAT_OR_EXIT\n");
    exit_error(callername,temp);
    pval = 0.0;
  }
  return pval;
}
/**************************************-**************************************/
/*                                                                           */
/*                           NDATA_GET_VC_FLT_EXIT                           */
/*                                                                           */
/*  Short name for above.                                                    */
/*                                                                           */
/*****************************************************************************/
float ndata_get_vc_flt_exit(nd,k,pname,callername)
     struct ndata_struct *nd;
     int k;
     char pname[],callername[];
{
  return ndata_get_vc_param_float_or_exit(nd,k,pname,callername);
}
/**************************************-**************************************/
/*                                                                           */
/*                            NDATA_GET_VC_INT_DFLT                          */
/*                                                                           */
/*****************************************************************************/
float ndata_get_vc_flt_dflt(nd,k,pname,dval)
     struct ndata_struct *nd;
     int k;
     char pname[];
     float dval;    // default value
{
  int flag;
  float fval;

  flag = ndata_get_vc_param_float(nd,k,pname,&fval);
  if (flag == 0)
    fval = dval;

  return fval;
}
/**************************************-**************************************/
/*                                                                           */
/*                        NDATA_MATCH_TRIALS_VAR_PARAMS                      */
/*                                                                           */
/*  Test consistency of the numbers and values of the variable paramters.    */
/*  Return values:                                                           */
/*                                                                           */
/*  -1 - trials have explicit mismatch in at least one value.                */
/*   0 - trials have no explicit mismatch, but param names not identical.    */
/*   1 - trials match in number and values of all variable parameters.       */
/*                                                                           */
/*****************************************************************************/
int ndata_match_trials_var_params(t1,t2)
     struct ndtrial_struct *t1,*t2;
{
  int i,k;
  int *tx,tn,flag;

  flag = 1;
  get_exclusive_index_compare_2d_carrays(t1->pname,t1->nparam,
					 t2->pname,t2->nparam,&tx,&tn);
  if (tn > 0){ /* If one list has elements not in the other. */
    flag = 0;
    myfree(tx);
  }
  for(i=0;i<t1->nparam;i++){
    k = search_2d_carray(t2->pname,t1->pname[i],t2->nparam);
    if (k >= 0)
      if (strcmp(t1->pval[i],t2->pval[k])!=0){
	flag = -1;
	return flag;
      }
  }
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                      NDATA_GET_TF_HZ_PEP_CYCLE_TRIAL                      */
/*                                                                           */
/*  Return the temporal frequency for the trial assuming PEP conventions:    */
/*  - A sync pulse is given before each cycle of a grating, and after the    */
/*    last cycle.                                                            */
/*  Return -1.0 if the number of cycles and syncs are incompatible.          */
/*                                                                           */
/*****************************************************************************/
float ndata_get_tf_hz_pep_cycle_trial(nd,k)
     struct ndata_struct *nd;
     int k;
{
  int flag,ncyc,np,s0,s1;
  float sampling,tf;
  struct ndtrial_struct *t;

  flag =  ndata_get_vc_param_int(nd,k,"cycles",&ncyc);
  if (flag != 1)
    exit_error("NDATA_GET_TF_HZ_PEP_CYCLE_TRIAL","Parameter cycles not found");

  t = &(nd->t[k]);
  np = ndata_trial_get_n_points_chan(t,"sync0");
  if (np <= 0)
    exit_error("NDATA_GET_TF_HZ_PEP_CYCLE_TRIAL","np <= 0");
  flag = ndata_trial_get_nth_point_chan(t,"sync0",0,&s0);
  flag = ndata_trial_get_nth_point_chan(t,"sync0",np-1,&s1);

  if (np != (ncyc+1)){
    printf("  *** WARNING NDATA_GET_TF_HZ_PEP_CYCLE_TRIAL Nsync != Ncyc+1\n");
    tf = -1.0;
  }else{
    sampling = ndata_trial_get_sampling_chan(t,"sync0");
    tf = (float)ncyc / ((float)(s1-s0)/sampling); /* Cycles per second */
  }

  return tf;
}
/**************************************-**************************************/
/*                                                                           */
/*                      NDATA_GET_DURATION_PEP_CYCLE_TRIAL                   */
/*                                                                           */
/*  Return the duration for the trial assuming PEP periodic conventions:     */
/*  - A sync pulse is given before each cycle of a grating, and after the    */
/*    last cycle.                                                            */
/*                                                                           */
/*  If "secflag" is 0, the value is returned in sampling units, if it is 1,  */
/*  the value is returned in seconds.                                        */
/*                                                                           */
/*****************************************************************************/
float ndata_get_duration_pep_cycle_trial(nd,k,secflag)
     struct ndata_struct *nd;
     int k,secflag;
{
  int flag,np,s0,s1;
  float sampling,duration;
  struct ndtrial_struct *t;

  t = &(nd->t[k]);
  np = ndata_trial_get_n_points_chan(t,"sync0");
  if (np <= 0){
    /*exit_error("NDATA_GET_DURATION_PEP_CYCLE_TRIAL","np <= 0");*/
    printf("    Note:  No sync pulses, using duration of record 0\n");
    if (secflag)
      duration = (float)t->r[0].tn/t->r[0].sampling;
    else
      duration = (float)t->r[0].tn;
    return duration;
  }
  flag = ndata_trial_get_nth_point_chan(t,"sync0",0,&s0);
  flag = ndata_trial_get_nth_point_chan(t,"sync0",np-1,&s1);
  if (!flag)
    exit_error("NDATA_GET_DURATION_PEP_CYCLE","Should not happen"); // Lint

  sampling = ndata_trial_get_sampling_chan(t,"sync0");
  if (secflag==1)
    duration = (float)(s1-s0)/sampling; // In seconds
  else
    duration = (float)(s1-s0); // In sampling units

  return duration;
}
/**************************************-**************************************/
/*                                                                           */
/*                             NDATA_ADD_CONST_PARAM                         */
/*                                                                           */
/*  Add a const param to the param list, or update the value if 'cname' is   */
/*  already in the list of constant params.                                  */
/*                                                                           */
/*****************************************************************************/
void ndata_add_const_param(nd,cname,ctype,cval)
     struct ndata_struct *nd;
     char *cname,ctype,*cval;
{
  int i,k;
  int tn;
  char **tname,**tval,*ttype;
  int ndata_get_const_param_index();

  /*** WYETH - SHOULD WE CHECK THAT THIS IS NOT A vparam NAME ??? ***/

  k = ndata_get_const_param_index(nd,cname);
  if (k >= 0){
    printf("    Parameter %s already exists, updating value\n",cname);
    myfree(nd->cval[k]);
    nd->cval[k] = strdup(cval);
    nd->ctype[k] = ctype;
  }else{
    tn = nd->nconst;
    tn += 1;
    tname = (char **)myalloc(tn*sizeof(char *));
    tval = (char **)myalloc(tn*sizeof(char *));
    ttype = (char *)myalloc(tn*sizeof(char)); 
    for(i=0;i<(tn-1);i++){
      tname[i] = strdup(nd->cname[i]);
      tval[i] = strdup(nd->cval[i]);
      ttype[i] = nd->ctype[i];
    }
    i = tn-1;
    ttype[i] = ctype;
    tname[i] = strdup(cname);
    tval[i] = strdup(cval);

    free_2d_carray(nd->cname,nd->nconst);
    free_2d_carray(nd->cval,nd->nconst);
    myfree(nd->ctype);

    nd->cname = tname;
    nd->cval = tval;
    nd->ctype = ttype;
    nd->nconst = tn;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              NDATA_ADD_VAR_PARAM                          */
/*                                                                           */
/*  Add a var param to the param list, and to each trial.                    */
/*                                                                           */
/*  Note: the value 'vval' is a constant, so the caller must change the      */
/*  value one a trial-by-trial basis if needed.                              */
/*                                                                           */
/*****************************************************************************/
void ndata_add_var_param(nd,vname,vtype,vval)
     struct ndata_struct *nd;
     char *vname,vtype,*vval;
{
  int i,j,k;
  int tn;
  char **tname,**tval,*ttype;
  int ndata_get_const_param_index();
  struct ndtrial_struct *tr;

  k = ndata_get_var_param_index(nd,vname);
  if (k >= 0){
    printf("    Parameter %s already exists.\n",vname);
    exit_error("NDATA_ADD_VAR_PARAM","Parameter already exists");
  }else{
    // Update var params in file header
    tn = nd->nvar;
    tn += 1;
    tname = (char **)myalloc(tn*sizeof(char *));
    ttype = (char *)myalloc(tn*sizeof(char)); 
    for(i=0;i<(tn-1);i++){
      tname[i] = strdup(nd->vname[i]);
      ttype[i] = nd->vtype[i];
    }
    i = tn-1;
    ttype[i] = vtype;
    tname[i] = strdup(vname);

    free_2d_carray(nd->vname,nd->nvar);
    myfree(nd->vtype);

    nd->vname = tname;
    nd->vtype = ttype;
    nd->nvar = tn;

    // Update var params in trials
    for(j=0;j<nd->ntrial;j++){
      tr = &(nd->t[j]);

      tn = tr->nparam;
      tn += 1;
      tname = (char **)myalloc(tn*sizeof(char *));
      tval = (char **)myalloc(tn*sizeof(char *));

      for(i=0;i<(tn-1);i++){
	tname[i] = strdup(tr->pname[i]);
	tval[i] = strdup(tr->pval[i]);
      }
      i = tn-1;
      tname[i] = strdup(vname);
      tval[i] = strdup(vval);

      free_2d_carray(tr->pname,tr->nparam);
      free_2d_carray(tr->pval,tr->nparam);

      tr->pname = tname;
      tr->pval = tval;
      tr->nparam = tn;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_GET_UNIQUE_VAR_VALUES                       */
/*                                                                           */
/*  Find the number of different values associated with the given parameter  */
/*  name.                                                                    */
/*                                                                           */
/*****************************************************************************/
void ndata_get_unique_var_values(nd,vname,ruval,rnuval)
     struct ndata_struct *nd;
     char vname[],***ruval;
     int *rnuval;
{
  int i,j,k;
  int flag,result;
  char *vstring,**uval; // Unique values
  
  uval = (char **)myalloc(nd->ntrial*sizeof(char *));
  k = 0;
  for(i=0;i<nd->ntrial;i++){
    result = ndata_get_var_value_trial(nd,i,vname,&vstring);
    if (result==0)
      exit_error("NDATA_GET_UNIQUE_VAR_VALUES","Variable not found");
    flag = 1; // Assume this value is new
    for(j=0;j<k;j++)
      if (strcmp(vstring,uval[j])==0)
	flag = 0;
    if (flag == 1){
      uval[k] = vstring;
      k += 1;
    }
  }
  *ruval = uval; *rnuval = k;
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_GET_VAR_PARAM_FREQENCY                      */
/*                                                                           */
/*  Determine the number of unique values and the frequencies of usage       */
/*  across trials for each value of the var param "vname".                   */
/*                                                                           */
/*****************************************************************************/
void ndata_get_var_param_frequency(nd,vname,ruval,rnuval,rvfrac)
     struct ndata_struct *nd;
     char vname[],***ruval;
     int *rnuval;
     float **rvfrac;
{
  int i,j;
  int nuval,result;
  char *vstring,**uval; /* Unique values. */
  float *vfrac; /* fraction of trials using the ith value */

  ndata_get_unique_var_values(nd,vname,&uval,&nuval);
  vfrac = get_zero_farray(nuval);
  
  for(i=0;i<nd->ntrial;i++){
    result = ndata_get_var_value_trial(nd,i,vname,&vstring);
    if (result==0)
      exit_error("NDATA_GET_VAR_PARAM_INFO","Variable not found");
    for(j=0;j<nuval;j++)
      if (strcmp(vstring,uval[j])==0)
	vfrac[j] += 1.0;
    myfree(vstring);
  }
  multiply_farray(vfrac,nuval,1.0/(float)nd->ntrial);

  *ruval = uval; *rnuval = nuval; *rvfrac = vfrac;
}
/**************************************-**************************************/
/*                                                                           */
/*                          NDATA_GET_VAR_PARAM_INFO                         */
/*                                                                           */
/*  For each variable parameter, determine the following:                    */
/*  - the number of unique values which it assumes across trials.            */
/*  - the most common value of that parameter.                               */
/*  - the fraction of trials containing the most common value.               */
/*                                                                           */
/*****************************************************************************/
void ndata_get_var_param_info(nd,rnuval,rfrac,rcval)
     struct ndata_struct *nd;
     int **rnuval;
     float **rfrac;
     char ***rcval;
{
  int i,k;
  int *nuval;
  float *frac,*tfrac;
  char **cval,**uval;

  nuval = (int *)myalloc(nd->nvar*sizeof(int));
  frac = (float *)myalloc(nd->nvar*sizeof(float));
  cval = (char **)myalloc(nd->nvar*sizeof(char *));
  for(i=0;i<nd->nvar;i++){
    ndata_get_var_param_frequency(nd,nd->vname[i],&uval,&nuval[i],&tfrac);
    k = max_coord_farray(tfrac,nuval[i]); /* Find (first) most common value. */
    cval[i] = strdup(uval[k]); /* Store the (or a) most common value. */
    frac[i] = tfrac[k]; /* Store the fraction of trials on which it occurs. */
    free_2d_carray(uval,nuval[i]);
    myfree(tfrac);
  }

  *rnuval = nuval; *rfrac = frac; *rcval = cval;
}
/**************************************-**************************************/
/*                                                                           */
/*                          NDATA_SET_VAR_VALUE_TRIAL                        */
/*                                                                           */
/*  Set the value of the var param.                                          */
/*                                                                           */
/*****************************************************************************/
void ndata_set_var_value_trial(nd,k,vname,vstr,freeflag)
     struct ndata_struct *nd;
     int k;
     char vname[],*vstr;
     int freeflag;
{
  int i;

  i = search_2d_carray(nd->t[k].pname,vname,nd->t[k].nparam);
  if (i == -1)
    exit_error("NDATA_SET_VAR_VALUE_TRIAL","Cannot find var param");
  else{
    if (freeflag)
      myfree(nd->t[k].pval[i]);
    nd->t[k].pval[i] = strdup(vstr);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                      NDATA_CONVERT_VAR_TYPE_INT_FLOAT                     */
/*                                                                           */
/*  Return values:                                                           */
/*    0 - not converted                                                      */
/*    1 - converted to int 'i'                                               */
/*    2 - converted to float 'f'                                             */
/*                                                                           */
/*****************************************************************************/
int ndata_convert_var_type_int_float(nd,vname)
     struct ndata_struct *nd;
     char *vname;
{
  int i;
  int nuval,iflag,fflag,nflag,retval;
  char **uval;

  iflag = 1;  // All are ints
  fflag = 1;  // All are floats
  nflag = 1;  // All are numbers (could be a mix of float and int)

  ndata_get_unique_var_values(nd,vname,&uval,&nuval);
  for(i=0;i<nuval;i++){
    if (is_int_string(uval[i]) == 0)
      iflag = 0;
  }
  for(i=0;i<nuval;i++){
    if (is_float_string(uval[i]) == 0)
      fflag = 0;
  }
  for(i=0;i<nuval;i++){
    if (is_number_string(uval[i]) == 0)
      nflag = 0;
  }
  free_2d_carray(uval,nuval);

  i = search_2d_carray(nd->vname,vname,nd->nvar);
  if (i == -1)
    exit_error("NDATA_CONVERT_VAR_TYPE_INT_FLOAT","var name not found");
    
  if (iflag == 1){
    nd->vtype[i] = 'i';
    retval = 1;
  }else if ((fflag == 1) || (nflag == 1)){  // All float, or mix of float/int
    nd->vtype[i] = 'f';
    retval = 2;
  }else
    retval = 0;

  return retval;
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_CONVERT_VAR_NUMERIC                         */
/*                                                                           */
/*  Attempt to convert all non-numeric ('c') type var params to int 'i' or   */
/*  to float 'f', if possible.                                               */
/*                                                                           */
/*****************************************************************************/
void ndata_convert_var_numeric(nd)
     struct ndata_struct *nd;
{
  int i;
  int flag,cflag;
  char vtype;

  printf("  NDATA_CONVERT_VAR_NUMERIC\n");

  for(i=0;i<nd->nvar;i++){
    flag = ndata_get_var_type(nd,nd->vname[i],&vtype);
    if (flag == 0)
      exit_error("NDATA_CONVERT_VAR_NUMERIC","vname not found");
    if (vtype == 'c'){
      cflag = ndata_convert_var_type_int_float(nd,nd->vname[i]);
      if (cflag == 1)
	printf("    %s converted to type 'i'\n",nd->vname[i]);
      else if (cflag == 2)
	printf("    %s converted to type 'f'\n",nd->vname[i]);
      else
	printf("    %s not converted.\n",nd->vname[i]);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               NDATA_PRINT_ECT                             */
/*                                                                           */
/*****************************************************************************/
void ndata_print_ect(nd,k)
     struct ndata_struct *nd;
     int k;
{
  int i;
  
  printf("  Table name:  %s\n",nd->table[k].tname);
  printf("  Table number:  %d\n",nd->table[k].tnum);
  printf("  Number of codes:  %d\n",nd->table[k].n);
  for(i=0;i<nd->table[k].n;i++)
    printf("  %16d %s\n",nd->table[k].num[i],nd->table[k].name[i]);
}
/**************************************-**************************************/
/*                                                                           */
/*                                PRINT_NDATA                                */
/*                                                                           */
/*****************************************************************************/
void print_ndata(nd)
     struct ndata_struct *nd;
{
  int i,j,k;
  
  printf("  PRINT_NDATA\n");
  
  printf("CLASS %s\n",nd->class);
  printf("NConst %d\n",nd->nconst);
  for(i=0;i<nd->nconst;i++)
    printf("  %s %s\n",nd->cname[i],nd->cval[i]);
  printf("NVar %d\n",nd->nvar);
  for(i=0;i<nd->nvar;i++)
    printf("  %s\n",nd->vname[i]);
  for(i=0;i<nd->ntable;i++){
    printf("TABLE %s %d\n",nd->table[i].tname,nd->table[i].n);
    for(j=0;j<nd->table[i].n;j++){
      printf("  %5d %s\n",nd->table[i].num[j],nd->table[i].name[j]);
    }
  }
  printf("NTrial %d\n",nd->ntrial);
  printf("--------------------------------------------------------------\n");
  for(i=0;i<nd->ntrial;i++){
    printf("  Tcode %d\n",nd->t[i].tcode);
    printf("  Tref %d\n",nd->t[i].tref);
    printf("  Nparam %d\n",nd->t[i].nparam);
    for(j=0;j<nd->t[i].nparam;j++)
      printf("    %s %s\n",nd->t[i].pname[j],nd->t[i].pval[j]);
    printf("  Nrec %d\n",nd->t[i].nrec);

    for(j=0;j<nd->t[i].nrec;j++){
      printf("%d %s %d %.2f %d %d %d\n",nd->t[i].r[j].rtype,nd->t[i].r[j].name,
	     nd->t[i].r[j].rcode,nd->t[i].r[j].sampling,nd->t[i].r[j].t0,
	     nd->t[i].r[j].tn,nd->t[i].r[j].n);
      if (nd->t[i].r[j].rtype == 0)
	for(k=0;k<nd->t[i].r[j].n;k++)
	  printf("%d ",nd->t[i].r[j].p[k]);
      else if (nd->t[i].r[j].rtype == 1) /*** USE "tn" ***/
	for(k=0;k<nd->t[i].r[j].n;k++)
	  printf("%.4f ",nd->t[i].r[j].x[k]);
      else if (nd->t[i].r[j].rtype == 2) /*** USE "n" ***/
	for(k=0;k<nd->t[i].r[j].n;k++)
	  printf("(%d %.4f) ",nd->t[i].r[j].p[k],nd->t[i].r[j].x[k]);
      else if (nd->t[i].r[j].rtype == 3){ /*** USE "n" ***/
	if (nd->t[i].r[j].ect != NULL)
	  printf("  TABLENUM %d\n",nd->t[i].r[j].ect->tnum);
	else
	  printf("  TABLENUM %d\n",-1);
	for(k=0;k<nd->t[i].r[j].n;k++)
	  printf("(%d %d) ",nd->t[i].r[j].p[k],nd->t[i].r[j].ec[k]);
      }else
	exit_error("PRINT_NDATA","Unknown record type");
      printf("\n");
    }
    printf("--------------------------------------------------------------\n");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            PRINT_PARTIAL_NDATA                            */
/*                                                                           */
/*****************************************************************************/
void print_partial_ndata(nd,option,k)
     struct ndata_struct *nd;
     char option[];
     int k;
{
  int i,j;
  int nuval,nn,c;
  char **uval,*name;

  if (strcmp(option,"all")==0)
    print_ndata(nd);
  else if (strcmp(option,"tables")==0){
    /*printf("NConst %d\n",nd->nconst);*/
    for(i=0;i<nd->ntable;i++)
      ndata_print_ect(nd,i);
  }else if (strcmp(option,"const")==0){
    /*printf("NConst %d\n",nd->nconst);*/
    for(i=0;i<nd->nconst;i++)
      printf("  %s = %s\n",nd->cname[i],nd->cval[i]);
  }else if (strcmp(option,"var")==0){
    /*printf("NVar %d\n",nd->nvar);*/
    for(i=0;i<nd->nvar;i++){
      ndata_get_unique_var_values(nd,nd->vname[i],&uval,&nuval);
      printf("  %s:",nd->vname[i]);
      for(j=0;j<nuval;j++)
	printf(" %s",uval[j]);
      printf("\n");
      free_2d_carray(uval,nuval);
    }
  }else if (strcmp(option,"header")==0){
    printf("  Class = %s\n",nd->class);
    printf("  Trials = %d\n",nd->ntrial);
  }else if (strcmp(option,"trial")==0){
    if ((k<0)||(k>=nd->ntrial))
      printf("  No such trial.\n");
    else{
      printf("  Tcode %d\n",nd->t[k].tcode);
      printf("  Tref %d\n",nd->t[k].tref);
      printf("  Nparam %d\n",nd->t[k].nparam);
      for(i=0;i<nd->t[k].nparam;i++)
	printf("    %s %s\n",nd->t[k].pname[i],nd->t[k].pval[i]);
      printf("  Nrec %d\n",nd->t[k].nrec);
      
      for(i=0;i<nd->t[k].nrec;i++){
	printf("%d %s %d %.2f %d %d %d\n",nd->t[k].r[i].rtype,
	       nd->t[k].r[i].name,nd->t[k].r[i].rcode,nd->t[k].r[i].sampling,
	       nd->t[k].r[i].t0,nd->t[k].r[i].tn,nd->t[k].r[i].n);
	if (nd->t[k].r[i].rtype == 0)
	  for(j=0;j<nd->t[k].r[i].n;j++)
	    printf("%d ",nd->t[k].r[i].p[j]);
	else if (nd->t[k].r[i].rtype == 1){ /*** USE "tn" ***/
	  nn = nd->t[k].r[i].tn;
	  if (nn > 20){
	    for(j=0;j<10;j++)
	      printf("%.4f ",nd->t[k].r[i].x[j]);
	    printf("\n...\n");
	    for(j=nn-10;j<nn;j++)
	      printf("%.4f ",nd->t[k].r[i].x[j]);
	    printf("\n");
	  }else{
	    for(j=0;j<nn;j++)
	      printf("%.4f ",nd->t[k].r[i].x[j]);
	    printf("\n");
	  }
	  /***for(j=0;j<nd->t[k].r[i].tn;j++)
	    printf("%.4f ",nd->t[k].r[i].x[j]);***/
	}else if (nd->t[k].r[i].rtype == 2) /*** USE "n" ***/
	  for(j=0;j<nd->t[k].r[i].n;j++)
	    printf("(%d %.4f) ",nd->t[k].r[i].p[j],nd->t[k].r[i].x[j]);
	else if (nd->t[k].r[i].rtype == 3){ /*** USE "n" ***/
	  for(j=0;j<nd->t[k].r[i].n;j++){
	    c = nd->t[k].r[i].ec[j];
	    ndata_lookup_ecode_name(nd->t[k].r[i].ect,c,&name);
	    if (name == NULL)
	      printf("  %10d %6d *** Not in table.\n",nd->t[k].r[i].p[j],c);
	    else
	      printf("  %10d %6d %s\n",nd->t[k].r[i].p[j],c,name);
	    myfree(name);
	  }
	}
	printf("\n");
      }
    }
  }else if (strcmp(option,"trialsum")==0){  // Summary of trial
    if ((k<0)||(k>=nd->ntrial))
      printf("  No such trial.\n");
    else{
      printf("  Tcode %d\n",nd->t[k].tcode);
      printf("  Tref %d\n",nd->t[k].tref);
      printf("  Nparam %d\n",nd->t[k].nparam);
      for(i=0;i<nd->t[k].nparam;i++)
	printf("    %s %s\n",nd->t[k].pname[i],nd->t[k].pval[i]);
      printf("  Nrec %d\n",nd->t[k].nrec);
      for(i=0;i<nd->t[k].nrec;i++)
	printf("    %d %s %d %.2f %d %d %d\n",nd->t[k].r[i].rtype,
	       nd->t[k].r[i].name,nd->t[k].r[i].rcode,nd->t[k].r[i].sampling,
	       nd->t[k].r[i].t0,nd->t[k].r[i].tn,nd->t[k].r[i].n);
    }
  }else if (strcmp(option,"trialhead")==0){
    if ((k<0)||(k>=nd->ntrial))
      printf("  No such trial.\n");
    else{
      printf("  Tcode %d\n",nd->t[k].tcode);
      printf("  Tref %d\n",nd->t[k].tref);
      printf("  Nparam %d\n",nd->t[k].nparam);
      for(i=0;i<nd->t[k].nparam;i++)
	printf("    %s %s\n",nd->t[k].pname[i],nd->t[k].pval[i]);
      printf("  Nrec %d\n",nd->t[k].nrec);
      /*
	for(i=0;i<nd->t[k].nrec;i++)
	printf("    %d %s %d %.2f %d %d %d\n",nd->t[k].r[i].rtype,
	nd->t[k].r[i].name,nd->t[k].r[i].rcode,nd->t[k].r[i].sampling,
	nd->t[k].r[i].t0,nd->t[k].r[i].tn,nd->t[k].r[i].n);
      */
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           NDU_GET_CHECK_ECT_CHAN                          */
/*                                                                           */
/*  Return a pointer to the event code table associated with a channel       */
/*  name, or NULL if the same ECT is not associated with each instance of    */
/*  the channel.                                                             */
/*                                                                           */
/*****************************************************************************/
struct ect_struct *ndu_get_check_ect_chan(nd,chan)
     struct ndata_struct *nd;
     char *chan;
{
  int i,k;
  int rtype,err_flag;
  struct ndrec_struct *rp;
  struct ect_struct *tt;

  //
  //  Check that 'chan' is always an event code channel (rec type 3)
  //
  rtype = ndata_get_type_record_name(nd,chan);
  if (rtype != 3)
    return NULL;

  tt = NULL;     // Table pointer
  err_flag = 0;  // Everything is OK

  for(i=0;i<nd->ntrial;i++){
    k = ndata_trial_get_record_index(&(nd->t[i]),chan);
    if (k >= 0){
      rp = &(nd->t[i].r[k]);
      if (tt == NULL){
	//
	//  Save pointer to first table found
	//
	if (rp->ect != NULL)
	  tt = rp->ect;
	else
	  exit_error("NDU_GET_CHECK_ECT_CHAN","No table for type 3 rec");
      }else if (tt != rp->ect)
	err_flag = 1;
    }
  }

  if (err_flag == 1)  // If there was an error, return NULL
    tt = NULL;

  return tt;
}
/**************************************-**************************************/
/*                                                                           */
/*                               NDATA_DIFF_ECT                              */
/*                                                                           */
/*  Return 0 if the two event code tables are the same, otherwise, return:   */
/*    1 - table names differ                                                 */
/*    2 - table numbers differ                                               */
/*    3 - number of codes differ                                             */
/*    4 - at least one code name or code number differs.                     */
/*                                                                           */
/*****************************************************************************/
int ndata_diff_ect(ect1,ect2)
     struct ect_struct *ect1,*ect2;
{
  int i;

  if (strcmp(ect1->tname,ect2->tname)!=0)
    return 1;
  else if (ect1->tnum != ect2->tnum)
    return 2;
  else if (ect1->n != ect2->n)
    return 3;
  else
    for(i=0;i<ect1->n;i++){
      if ((strcmp(ect1->name[i],ect2->name[i])!=0)||
	  (ect1->num[i] != ect2->num[i]))
	return 4;
    }
  return 0;
}
/**************************************-**************************************/
/*                                                                           */
/*                               NDATA_MERGE_ECT                             */
/*                                                                           */
/*****************************************************************************/
struct ect_struct *ndata_merge_ect(ect1,ect2,tname,tnum)
     struct ect_struct *ect1,*ect2;
     char *tname;
     int tnum;
{
  int i,j,k;
  int *t,*ti,n,*flag,count;
  char **tc;
  struct ect_struct *ect;

  /* Make a list 't' of all codes, get an index 'ti' into the sorted list */
  n = ect1->n + ect2->n;
  t = (int *)myalloc(n*sizeof(int));
  tc = (char **)myalloc(n*sizeof(char *));
  for(i=0;i<ect1->n;i++){
    t[i] = ect1->num[i];
    tc[i] = ect1->name[i];
  }
  for(i=0;i<ect2->n;i++){
    t[ect1->n+i] = ect2->num[i];
    tc[ect1->n+i] = ect2->name[i];
  }
  ti = index_bubble_sort_iarray(t,n);

  /* Look for inconsistencies in the list.  Reuse 't' for flag array. */
  flag = get_zero_iarray(n);
  count = 0;
  for(i=0;i<(n-1);i++){
    j = ti[i];
    k = ti[i+1];
    if (t[j] != t[k]){
      flag[i] = 1;
      count += 1;
    }else{
      if (strcmp(tc[j],tc[k])!=0){
	printf("  *** Code %d called %s and %s\n",t[j],tc[j],tc[k]);
	exit_error("NDATA_MERGE_ECT","Name conflict for ecode");
      }
    }
  }
  flag[n-1] = 1; /* Last code must be valid. */
  count += 1;

  ect = (struct ect_struct *)myalloc(sizeof(struct ect_struct));
  ect->tname = strdup(tname);
  ect->tnum = tnum;
  ect->n = count;
  ect->name = (char **)myalloc(ect->n*sizeof(char *));
  ect->num = (int *)myalloc(ect->n*sizeof(int));
  k = 0;
  for(i=0;i<n;i++)
    if (flag[i] == 1){
      ect->name[k] = strdup(tc[ti[i]]);
      ect->num[k] = t[ti[i]];
      k += 1;
    }

  myfree(t); myfree(ti); myfree(flag); myfree(tc);

  return ect;
}
/**************************************-**************************************/
/*                                                                           */
/*                                NDATA_COPY_ECT                             */
/*                                                                           */
/*****************************************************************************/
struct ect_struct *ndata_copy_ect(ect)
     struct ect_struct *ect;
{
  struct ect_struct *t;

  t = (struct ect_struct *)myalloc(sizeof(struct ect_struct));
  t->tname = strdup(ect->tname);
  t->tnum = ect->tnum;
  t->n = ect->n;
  t->name = copy_2d_carray(ect->name,ect->n);
  t->num = copy_iarray(ect->num,ect->n);
  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                             GET_EMPTY_CONDITION                           */
/*                                                                           */
/*  Return a condition structure with storage for "n" conditions.  Called    */
/*  with n=0, this will allow all trials.                                    */
/*                                                                           */
/*****************************************************************************/
struct ndcond_struct *get_empty_condition(n)
     int n;
{
  struct ndcond_struct *ndc;

  ndc = (struct ndcond_struct *)myalloc(sizeof(struct ndcond_struct));
  ndc->n = n;
  if (n > 0){
    ndc->ns = (int *)myalloc(n*sizeof(int));
    ndc->slist = (char ***)myalloc(n*sizeof(char **));
  }else{
    ndc->ns = NULL;
    ndc->slist = NULL;
  }
  
  return ndc;
}
/**************************************-**************************************/
/*                                                                           */
/*                             NDCOND_FREE_CONDITION                         */
/*                                                                           */
/*  Free a condition structure.                                              */
/*                                                                           */
/*****************************************************************************/
void ndcond_free_condition(ndc)
  struct ndcond_struct *ndc;
{
  int i;

  for(i=0;i<ndc->n;i++)
    free_2d_carray(ndc->slist[i],ndc->ns[i]);
  myfree(ndc->slist);
  myfree(ndc->ns);
  myfree(ndc);
}
/**************************************-**************************************/
/*                                                                           */
/*                            PRINT_CONDITION_LIST                           */
/*                                                                           */
/*****************************************************************************/
void print_condition_list(ndc)
     struct ndcond_struct *ndc;
{
  int i,j;

  if (ndc->n == 0)
    printf("  No conditions.\n");
  for(i=0;i<ndc->n;i++){
    printf("  (%d) %s",i,ndc->slist[i][0]);
    for(j=1;j<ndc->ns[i];j++)
      printf(" %s",ndc->slist[i][j]);
    printf("\n");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             TEST_CONDITION_LIST                           */
/*                                                                           */
/*****************************************************************************/
void test_condition_list(nd,ndc)
     struct ndata_struct *nd;
     struct ndcond_struct *ndc;
{
  int count_trials_ndata_condition(); /* Defined below */

  if (nd != NULL)
    printf("  %d of %d trials match conditions.\n",
	   count_trials_ndata_condition(nd,ndc),nd->ntrial);
  else
    printf("  No ndata file.\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                                NDCOND_REMOVE                              */
/*                                                                           */
/*  Remove the "k"th condition from the list.                                */
/*                                                                           */
/*****************************************************************************/
void ndcond_remove(ndc,k)
     struct ndcond_struct *ndc;
     int k;
{
  int i,j;
  int n,*ns;
  char ***slist;

  n = ndc->n;
  if ((k<0)||(k>=n))
    printf("  No such condition.\n");
  else{
    ns = (int *)myalloc((n-1)*sizeof(int)); /*** BUILD NEW LIST ***/
    slist = (char ***)myalloc((n-1)*sizeof(char **));
    j = 0;
    for(i=0;i<n;i++){
      if (i!=k){
	ns[j] = ndc->ns[i];
	slist[j] = copy_2d_carray(ndc->slist[i],ndc->ns[i]);
	j += 1;
      }
    }
    for(i=0;i<n;i++) /*** FREE OLD LIST ***/
      free_2d_carray(ndc->slist[i],ndc->ns[i]);
    myfree(ndc->slist);
    myfree(ndc->ns);

    ndc->n = n-1; /*** SET POINTERS TO NEW LIST ***/
    ndc->ns = ns;
    ndc->slist = slist;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                 NDCOND_ADD                                */
/*                                                                           */
/*  Add a condition to the list.                                             */
/*                                                                           */
/*****************************************************************************/
void ndcond_add(ndc,tlist,tn)
     struct ndcond_struct *ndc;
     char **tlist;
     int tn;
{
  int i;
  int n,*ns;
  char ***slist;

  n = ndc->n;
  if (tn<1)
    printf("*** NDCOND_ADD  Warning, empty condition ignored.\n");
  else{
    ns = (int *)myalloc((n+1)*sizeof(int)); // BUILD NEW LIST
    slist = (char ***)myalloc((n+1)*sizeof(char **));
    for(i=0;i<n;i++){
      ns[i] = ndc->ns[i];
      slist[i] = copy_2d_carray(ndc->slist[i],ndc->ns[i]);
    }
    ns[n] = tn; // ADD NEW CONDITION
    slist[n] = copy_2d_carray(tlist,tn);
    
    for(i=0;i<n;i++) // FREE OLD LIST
      free_2d_carray(ndc->slist[i],ndc->ns[i]);
    myfree(ndc->slist);
    myfree(ndc->ns);
    
    ndc->n = n+1; // SET POINTERS TO NEW LIST
    ndc->ns = ns;
    ndc->slist = slist;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           NDCOND_CHECK_POINT_COUNT                        */
/*                                                                           */
/*****************************************************************************/
int ndcond_check_point_count(nd,k,slist,ns)
     struct ndata_struct *nd;
     int k;
     char **slist;
     int ns;
{
  int i;
  int *s,cnt,flag,start,period,min,max,count,pctype;

  if (strcmp(slist[0],"point_count_min")==0)
    pctype = 1;
  else if (strcmp(slist[0],"point_count_max")==0)
    pctype = 2;
  else if (strcmp(slist[0],"point_count")==0)
    pctype = 3;
  else{
    pctype = -1;
    exit_error("NDCOND_CHECK_POINT_COUNT","Unknown point cound condition");
  }

  flag = 1;
  if ((pctype == 1) || (pctype == 2)){  // min or max
    if (ns != 5){
      printf("  Expecting point_count_min/max <chan_name> <start> <period> ");
      printf("<min/max_points>\n");
    }else{
      start = atoi(slist[2]);
      period = atoi(slist[3]);
      min = atoi(slist[4]);  // might be max value
      i = ndata_trial_get_record_index(&(nd->t[k]),slist[1]);
      if (i < 0)
	exit_error("NDCOND_CHECK_POINT_COUNT","Channel not found");
      s = nd->t[k].r[i].p;
      cnt = nd->t[k].r[i].n;
      count = count_spikes(s,cnt,start,period);
      if ((pctype == 1) && (count < min))
	flag = 0;
      else if ((pctype == 2) && (count > min))  // 'min' holds max
	flag = 0;
    }
  }else if (pctype == 3){
    if (ns != 6){
      printf("  Expecting:  point_count <chan_name> <start> <period> ");
      printf("<min_points> <max_points>\n");
    }else{
      start = atoi(slist[2]);
      period = atoi(slist[3]);
      min = atoi(slist[4]);
      max = atoi(slist[5]);
      i = ndata_trial_get_record_index(&(nd->t[k]),slist[1]);
      if (i < 0)
	exit_error("NDCOND_CHECK_POINT_COUNT","Channel not found");
      s = nd->t[k].r[i].p;
      cnt = nd->t[k].r[i].n;
      count = count_spikes(s,cnt,start,period);
      if ((count < min)  || (count > max))
	flag = 0;
    }
  }else
    printf(" *** Unknown condition is ignored.\n");

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDCOND_CHECK_EVCO_SEQUENCE                        */
/*                                                                           */
/*****************************************************************************/
int ndcond_check_evco_sequence(nd,k,slist,ns)
     struct ndata_struct *nd;
     int k;
     char **slist;
     int ns;
{
  int flag;

  printf(" *** WYETH - ACCEPTING ALL TRIALS - condition not implemented\n");

  flag = 1;

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                          NDCOND_CHECK_EVCO_VALUE                          */
/*                                                                           */
/*****************************************************************************/
int ndcond_check_evco_value(nd,k,slist,ns)
     struct ndata_struct *nd;
     int k;
     char **slist;
     int ns;
{
  int i;
  int flag,min,max,evt1,evt2,t1,t2,tdiff;
  int *clist,cn;
  struct ndtrial_struct *t;
  struct ndrec_struct *rp;

  t = &(nd->t[k]);  // Trial pointer

  flag = 1;  //  Assume the condition is met

  if (strcmp(slist[0],"ecode_range")==0){
    if (ns != 4){
      exit_error("NDCOND_CHECK_EVCO_VALUE",
		 "Expecting:  ecode_range <chan_name> <min> <max>\n");
    }else{
      min = atoi(slist[2]);
      max = atoi(slist[3]);
      clist = ndata_trial_getp_evco_list_chan(t,slist[1],&cn);
      // Note: don't free clist, it is a pointer to other storage
      i = get_index_search_range_iarray(clist,cn,min,max);
      if (i < 0)
	flag = 0;
    }
  }else if (strcmp(slist[0],"ecode_name")==0){
    if (ns != 3){
      exit_error("NDCOND_CHECK_EVCO_VALUE",
		 "Expecting:  ecode_name <chan_name> <event_name_string>\n");
    }else{
      rp = ndata_trial_get_record_ptr(t,slist[1]);  // Record pointer
      flag = ndata_lookup_ecode_num(rp->ect,slist[2],&evt1);
      if (flag == 0){
	printf("  Event:  %s\n",slist[2]);
	exit_error("NDCOND_CHECK_EVCO_VALUE","Event name not found in table");
      }

      flag = ndata_trial_chan_evco_test(t,slist[1],evt1);
    }
  }else if (strcmp(slist[0],"ecode_time_diff")==0){
    //
    //   0     1     2     3        4         5         6
    //  ... <chan> <id1> <id2> <condition> <value1> [<value2>]
    //
    //    where <condition> is one of the following strings:
    //
    //      >  <  >=  <=  =  range
    //
    if (ns < 6){
      exit_error("NDCOND_CHECK_EVCO_VALUE",
		 "Expecting:  ecode_time_diff <chan_name> <event_id_1>...\n");
    }else{

      //
      //  Set 'tdiff' to the time difference
      //
      flag = ndu_trial_chan_ecode_str_diff(t,slist[1],slist[2],slist[3],
					   &tdiff);

      if (flag){  // Both codes exist

	//printf("tdiff = %d\n",tdiff);

	flag = 0; // Assume failure
	if (strcmp("<",slist[4])==0){
	  if (tdiff < atoi(slist[5]))
	    flag = 1;
	}else if (strcmp("<=",slist[4])==0){
	  if (tdiff <= atoi(slist[5]))
	    flag = 1;
	}else if (strcmp(">",slist[4])==0){
	  if (tdiff > atoi(slist[5]))
	    flag = 1;
	}else if (strcmp(">=",slist[4])==0){
	  if (tdiff >= atoi(slist[5]))
	    flag = 1;
	}else if (strcmp("=",slist[4])==0){
	  if (tdiff == atoi(slist[5]))
	    flag = 1;
	}else if (strcmp("range",slist[4])==0){
	  if ((tdiff >= atoi(slist[5])) &&
	      (tdiff <= atoi(slist[6])))
	    flag = 1;
	}
      }
    }
  }else
    exit_error("NDCOND_CHECK_EVCO_VALUE","Unknown condition\n");

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                          NDCOND_CHECK_VAR_PARAM                           */
/*                                                                           */
/*****************************************************************************/
int ndcond_check_var_param(nd,k,slist,ns)
     struct ndata_struct *nd;
     int k;
     char **slist;
     int ns;
{
  int flag,t;
  char ptype,*pval;

  // ********** WYETH - IS THIS UNUSED ???? *************
  // ********** WYETH - IS THIS UNUSED ???? *************
  // ********** WYETH - IS THIS UNUSED ???? *************

  flag = 1;
  if (strcmp(slist[0],"param_range")==0){
    if (ns != 4){
      printf("  Expecting param_range <param_name> <min> <max>\n");
    }else{
      t = ndata_get_var_type(nd,slist[1],&ptype);
      if (t<=0){
	printf("*** WARNING  Param name (%s) not found, condition ignored.\n",
	       slist[1]);
	return flag;
      }
      t = ndata_get_var_value_trial(nd,k,slist[1],&pval);
      if (ptype=='i'){
	if ((atoi(pval) < atoi(slist[2]))||(atoi(pval) > atoi(slist[3])))
	  flag = 0;
      }else if (ptype=='f'){
	if ((atof(pval) < atof(slist[2]))||(atof(pval) > atof(slist[3])))
	  flag = 0;
      }else{
	printf("*** WARNING  Param type is char, condition ignored.\n");
	return flag;
      }
    }
  }else
    printf(" *** Unknown condition is ignored.\n");

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                           NDCOND_CHECK_VC_PARAM                           */
/*                                                                           */
/*****************************************************************************/
int ndcond_check_vc_param(nd,k,slist,ns)
     struct ndata_struct *nd;
     int k;  /* Trial number begins at zero */
     char **slist;
     int ns;
{
  int i;
  int flag,t;
  char ptype,*pval;

  flag = 1;
  if ((strcmp(slist[0],"param_range")==0) ||
      (strcmp(slist[0],"param_range_prev")==0) ||
      (strcmp(slist[0],"param_range_prev1")==0)){
    if (ns != 4){
      printf("  Expecting %s <param_name> <min> <max>\n",slist[0]);
      exit_error("NDCOND_CHECK_VC_PARAM","Format error");
    }else{
      t = ndata_get_var_or_const_type(nd,slist[1],&ptype);
      if (t<=0){
	printf("*** WARNING  Param name (%s) not found, condition ignored.\n",
	       slist[1]);
	return flag;
      }
      if (strcmp(slist[0],"param_range")==0){
	t = ndata_get_var_or_const_param(nd,k,slist[1],&pval);
      }else{
	if (k > 0)
	  t = ndata_get_var_or_const_param(nd,k-1,slist[1],&pval);
	else{
	  if (strcmp(slist[0],"param_range_prev")==0){
	    flag = 0;     /* No previous trial, must fail */
	    return flag;
	  }else{
	    flag = 1;

	    /*** WYETH - test/debug ***/
	    printf("  Trial %d succeeded\n",k);

	    return flag;  /* No previous trial, must succeed */
	  }
	}
      }
      
      if (ptype=='i'){
	if ((atoi(pval) < atoi(slist[2]))||(atoi(pval) > atoi(slist[3])))
	  flag = 0;
      }else if (ptype=='f'){
	if ((atof(pval) < atof(slist[2]))||(atof(pval) > atof(slist[3])))
	  flag = 0;
      }else{
	printf("*** WARNING  Param type is char, condition ignored.\n");
	return flag;
      }
    }
  }else if (strcmp(slist[0],"param_prev_range")==0){

    /**** WEYTH ****/
    if (ns != 4){
      printf("  Expecting param_prev_range <param_name> <min> <max>\n");
      exit_error("NDCOND_CHECK_VC_PARAM","Format error");
    }else{
      t = ndata_get_var_or_const_type(nd,slist[1],&ptype);
      if (t<=0){
	printf("*** WARNING  Param name (%s) not found, condition ignored.\n",
	       slist[1]);
	return flag;
      }
      t = ndata_get_var_or_const_param(nd,k,slist[1],&pval);
      if (ptype=='i'){
	if ((atoi(pval) < atoi(slist[2]))||(atoi(pval) > atoi(slist[3])))
	  flag = 0;
      }else if (ptype=='f'){
	if ((atof(pval) < atof(slist[2]))||(atof(pval) > atof(slist[3])))
	  flag = 0;
      }else{
	printf("*** WARNING  Param type is char, condition ignored.\n");
	return flag;
      }
    }
  }else if (strcmp(slist[0],"param_val")==0){
    if (ns != 3){
      printf("  Expecting param_val <param_name> <value>\n");
      exit_error("NDCOND_CHECK_VC_PARAM","Format error");
    }else{
      t = ndata_get_var_or_const_type(nd,slist[1],&ptype);
      if (t<=0){
	printf("*** WARNING  Param name (%s) not found, condition ignored.\n",
	       slist[1]);
	return flag;
      }
      t = ndata_get_var_or_const_param(nd,k,slist[1],&pval);
      if (ptype=='i'){
	if (atoi(pval) != atoi(slist[2]))
	  flag = 0;
      }else if (ptype=='f'){
	if (atof(pval) != atof(slist[2]))
	  flag = 0;
      }else{
	if (strcmp(pval,slist[2])!=0)
	  flag = 0;
      }
    }
  }else if (strcmp(slist[0],"param_list")==0){
    if (ns < 3){
      printf("  Expecting param_list <param_name> <val1> ... <valn>\n");
      exit_error("NDCOND_CHECK_VC_PARAM","Format error");
    }else{
      t = ndata_get_var_or_const_type(nd,slist[1],&ptype);
      if (t<=0){
	printf("*** WARNING  Param name (%s) not found, condition ignored.\n",
	       slist[1]);
	return flag;
      }
      t = ndata_get_var_or_const_param(nd,k,slist[1],&pval);
      if (ptype=='i'){
	flag = 0;
	for(i=0;i<(ns-2);i++)
	  if (atoi(pval) == atoi(slist[i+2]))
	    flag = 1;
      }else if (ptype=='f'){
	flag = 0;
	for(i=0;i<(ns-2);i++)
	  if (atof(pval) == atof(slist[i+2]))
	    flag = 1;
      }else{
	flag = 0;
	for(i=0;i<(ns-2);i++)
	  if (strcmp(pval,slist[i+2])==0)
	    flag = 1;
      }
    }
  }else
    printf(" *** Unknown condition is ignored.\n");


  /*** WYETH - test/debug ***/
  if ((strcmp(slist[0],"param_range_prev")==0) ||
      (strcmp(slist[0],"param_range_prev1")==0)){
    if (flag == 1)
      printf("  Trial %d succeeded\n",k);
  }

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDCOND_CHECK_SPECIAL_DECISION                     */
/*                                                                           */
/*  WYETH - this is a temporaray hack.  This condition looks for pref'd or   */
/*  null decisions (rather than correct/incorrect) using the two variables   */
/*  "response" and "direction" (in old Newsome data).                        */
/*                                                                           */
/*****************************************************************************/
int ndcond_check_special_decision(nd,k,slist,ns)
     struct ndata_struct *nd;
     int k;
     char **slist;
     int ns;
{
  int flag,t,dir,resp;

  flag = 1;
  if (strcmp(slist[0],"special_decision")==0){
    if (ns != 2){
      printf("  Expecting special_decision [pref|null]\n");
    }else{
      t = ndata_get_vc_param_int(nd,k,"direction",&dir);
      if (t != 1) exit_error("NDCOND_CHECK_SPECIAL_DECISION","direction");
      t = ndata_get_vc_param_int(nd,k,"response",&resp);
      if (t != 1) exit_error("NDCOND_CHECK_SPECIAL_DECISION","response");
      if ((strcmp(slist[1],"null")==0) &&
	  (((dir==1)&&(resp==1)) || ((dir==0)&&(resp==0))))
	flag = 0;
      else if ((strcmp(slist[1],"pref")==0) &&
	       (((dir==1)&&(resp==0)) || ((dir==0)&&(resp==1))))
	flag = 0;
    }
  }else
    printf(" *** Unknown condition is ignored.\n");

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                      CHECK_TRIAL_NDATA_CONDITION_ONE                      */
/*                                                                           */
/*  Check the "c"th condition for the "k"th trial.                           */
/*                                                                           */
/*  If "chan" == NULL, condition is not checked, and a 1 is returned.        */
/*                                                                           */
/*****************************************************************************/
int check_trial_ndata_condition_one(nd,k,ndc,c)
     struct ndata_struct *nd;
     int k;
     struct ndcond_struct *ndc;
     int c;
{
  int flag;

  flag = 1;
  if (ndc->ns[c] <= 0)
    printf("*** CHECK_TRIAL_NDATA_CONDITION_ONE  Ignoring empty condition.\n");
  else if ((strcmp(ndc->slist[c][0],"point_count_min")==0)||
	   (strcmp(ndc->slist[c][0],"point_count_max")==0)||
	   (strcmp(ndc->slist[c][0],"point_count")==0))
    flag = ndcond_check_point_count(nd,k,ndc->slist[c],ndc->ns[c]);
  else if (strcmp(ndc->slist[c][0],"ecode_range")==0)
    flag = ndcond_check_evco_value(nd,k,ndc->slist[c],ndc->ns[c]);
  else if (strcmp(ndc->slist[c][0],"ecode_name")==0)
    flag = ndcond_check_evco_value(nd,k,ndc->slist[c],ndc->ns[c]);
  else if (strcmp(ndc->slist[c][0],"ecode_time_diff")==0)
    flag = ndcond_check_evco_value(nd,k,ndc->slist[c],ndc->ns[c]);
  else if (strcmp(ndc->slist[c][0],"ecode_sequence")==0)
    flag = ndcond_check_evco_sequence(nd,k,ndc->slist[c],ndc->ns[c]);
  else if (strcmp(ndc->slist[c][0],"param_val")==0)
    flag = ndcond_check_vc_param(nd,k,ndc->slist[c],ndc->ns[c]);
  else if (strcmp(ndc->slist[c][0],"param_list")==0)
    flag = ndcond_check_vc_param(nd,k,ndc->slist[c],ndc->ns[c]);
  else if ((strcmp(ndc->slist[c][0],"param_range")==0) ||
	   (strcmp(ndc->slist[c][0],"param_range_prev")==0) ||
	   (strcmp(ndc->slist[c][0],"param_range_prev1")==0))
    flag = ndcond_check_vc_param(nd,k,ndc->slist[c],ndc->ns[c]);
  else if (strcmp(ndc->slist[c][0],"special_decision")==0)
    flag = ndcond_check_special_decision(nd,k,ndc->slist[c],ndc->ns[c]);
  else if (strcmp(ndc->slist[c][0],"trial_number")==0){
    if (atoi(ndc->slist[c][1]) != k)
      flag = 0;
  }else if (strcmp(ndc->slist[c][0],"trial_last")==0){
    if (k != (nd->ntrial-1))
      flag = 0;
  }else if (strcmp(ndc->slist[c][0],"trial_even")==0){
    if (k%2 == 1)
      flag = 0;
  }else if (strcmp(ndc->slist[c][0],"trial_odd")==0){
    if (k%2 == 0)
      flag = 0;
  }else if (strcmp(ndc->slist[c][0],"trial_range")==0){
    if ((k < atoi(ndc->slist[c][1])) || (k > atoi(ndc->slist[c][2])))
      flag = 0;
  }else{
    printf("*** CHECK_TRIAL_NDATA_CONDITION_ONE Unknown condition ignored.\n");
    printf("%s\n",ndc->slist[c][0]);
  }

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                          NDCOND_POINT_COUNT_MAX_DEV                       */
/*                                                                           */
/*  Check the condition specified in "slist" for *all* trials, set flags     */
/*  to zero for trials that fail to meet criterion.                          */
/*                                                                           */
/*****************************************************************************/
void ndcond_point_count_max_dev(nd,slist,ns,flag)
     struct ndata_struct *nd;
     char **slist;
     int ns,*flag;
{
  int i;
  int start,period;
  float nsd,*data,mean,sdev;
  struct ndrx_struct *rx;

  if (strcmp(slist[0],"point_count_max_dev")==0){
    if (ns != 5){
      printf("  Expecting  point_count_max_dev <chan_name> <start> <period> ");
      printf("<nSD>\n");
    }else{
      start = atoi(slist[2]);
      period = atoi(slist[3]);
      nsd = atof(slist[4]);
      ndata_get_record_index(nd,slist[1],&rx);
      if (rx->n != nd->ntrial)
	exit_error("NDCOND_POINT_COUNT_MAX_DEV","Index misses some trials");
      data = (float *)myalloc(rx->n*sizeof(float));
      for(i=0;i<rx->n;i++){
	data[i] = (float)count_spikes(rx->r[i]->p,rx->r[i]->n,start,period);
      }
      mean_sdev_farray(data,rx->n,&mean,&sdev);
      for(i=0;i<rx->n;i++)
	if ((((data[i] - mean)/sdev) > nsd)||
	    ((-(data[i] - mean)/sdev) > nsd)){
	  printf("  Excluding trial %d:  sc=%d mean %.2f SD %.2f\n",i,
		 (int)data[i],mean,sdev);
	  flag[i] = 0;
	}
      myfree(data);
      ndata_free_ndrx(rx);
    }
  }else
    exit_error("NDCOND_POINT_COUNT_MAX_DEV","Unexpected failure");
}
/**************************************-**************************************/
/*                                                                           */
/*                    NDCOND_POINT_COUNT_MAX_DEV_TWO_SIDED                   */
/*                                                                           */
/*  Check the condition specified in "slist" for *all* trials, set flags     */
/*  to zero for trials that fail to meet criterion.                          */
/*                                                                           */
/*****************************************************************************/
void ndcond_point_count_max_dev_two_sided(nd,slist,ns,flag)
     struct ndata_struct *nd;
     char **slist;
     int ns,*flag;
{
  int i;
  int start,period;
  float nsdl,nsdr,*data,mean,sdl,sdr,dev;
  struct ndrx_struct *rx;

  if (strcmp(slist[0],"point_count_max_dev_two_sided")==0){
    if (ns != 6){
      printf("  Expecting  point_count_max_dev_two_sided <chan_name> <start>");
      printf(" <period> <nSD_left> <nSD_right>\n");
    }else{
      start = atoi(slist[2]);
      period = atoi(slist[3]);
      nsdl = atof(slist[4]);
      nsdr = atof(slist[5]);
      ndata_get_record_index(nd,slist[1],&rx);
      if (rx->n != nd->ntrial)
	exit_error("NDCOND_POINT_COUNT_MAX_DEV_TWO_SIDE","Indices missing");
      data = (float *)myalloc(rx->n*sizeof(float));
      for(i=0;i<rx->n;i++){
	data[i] = (float)count_spikes(rx->r[i]->p,rx->r[i]->n,start,period);
      }
      mean_two_sided_sdev_farray(data,rx->n,&mean,&sdl,&sdr);
      /*mean_sdev_farray(data,rx->n,&mean,&sdev);*/
      for(i=0;i<rx->n;i++){
	dev = data[i] - mean;
	if (((dev > 0.0) && (dev/sdr > nsdr)) ||
	    ((dev < 0.0) && (-dev/sdl > nsdl))){
	  printf("  Excluding trial %d:  SC=%d  mean %.2f lSD %.2f rSD %.2f\n",
		 i,(int)data[i],mean,sdl,sdr);
	  flag[i] = 0;
	}
      }
      myfree(data);
      ndata_free_ndrx(rx);
    }
  }else
    exit_error("NDCOND_POINT_COUNT_MAX_DEV_TWO_SIDED","Unexpected failure");
}
/**************************************-**************************************/
/*                                                                           */
/*                          NDATA_CONDITION_GET_FLAGS                        */
/*                                                                           */
/*  Return a list of ints that flags trials as accepted, 1, or unaccepted,   */
/*  0.  This routine allows implementation of conditions that require        */
/*  measures to be made across all trials, for example, mean and SD of       */
/*  firing rate.                                                             */
/*                                                                           */
/*****************************************************************************/
int *ndata_condition_get_flags(nd,ndc)
     struct ndata_struct *nd;
     struct ndcond_struct *ndc;
{
  int i,j;
  int *flag;

  flag = get_const_iarray(nd->ntrial,1);
  for(i=0;i<ndc->n;i++){
    if (strcmp(ndc->slist[i][0],"point_count_max_dev")==0)
      ndcond_point_count_max_dev(nd,ndc->slist[i],ndc->ns[i],flag);
    else if (strcmp(ndc->slist[i][0],"point_count_max_dev_two_sided")==0)
      ndcond_point_count_max_dev_two_sided(nd,ndc->slist[i],ndc->ns[i],flag);
    else{
      for(j=0;j<nd->ntrial;j++){
	if (check_trial_ndata_condition_one(nd,j,ndc,i) == 0)
	  flag[j] = 0;
      }
    }
  }
  
  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                         COUNT_TRIALS_NDATA_CONDITION                      */
/*                                                                           */
/*  Return the number of trials which satisfy the conditions.                */
/*                                                                           */
/*****************************************************************************/
int count_trials_ndata_condition(nd,ndc)
     struct ndata_struct *nd;
     struct ndcond_struct *ndc;
{
  int count,*cflag;

  cflag = ndata_condition_get_flags(nd,ndc);
  count = count_non_zero_iarray(cflag,nd->ntrial);
  myfree(cflag);

  return count;
}
/**************************************-**************************************/
/*                                                                           */
/*                     COUNT_TRIALS_NDATA_CONDITION_CHAN                     */
/*                                                                           */
/*  Return the number of trials which satisfy the conditions and contain     */
/*  the channel name "chan".                                                 */
/*                                                                           */
/*****************************************************************************/
int count_trials_ndata_condition_chan(nd,ndc,chan)
     struct ndata_struct *nd;
     struct ndcond_struct *ndc;
     char *chan;
{
  int i,k;
  int count,*cflag;

  cflag = ndata_condition_get_flags(nd,ndc);

  count = 0;
  for(i=0;i<nd->ntrial;i++){
    k = ndata_trial_get_record_index(&(nd->t[i]),chan);
    if (k >= 0)
      if (cflag[i])
	count += 1;
  }
  myfree(cflag);
  return count;
}
/**************************************-**************************************/
/*                                                                           */
/*                         GET_SARRAY_NDATA_CONDITION                        */
/*                                                                           */
/*  This mallocs a list of pointers into the "nd" data structure, thus,      */
/*  only the list of pointers, not what they point to, should be freed.      */
/*                                                                           */
/*  Note:  The "chan" specified here is the one from which the sarray is     */
/*         taken, which may be different than the one that is used in the    */
/*         condition data structure, "ndc".                                  */
/*                                                                           */
/*****************************************************************************/
void get_sarray_ndata_condition(nd,ndc,chan,rs,rcnt,rn)
     struct ndata_struct *nd;
     struct ndcond_struct *ndc;
     char *chan;
     int ***rs,**rcnt,*rn;
{
  int i,j,k;
  int count,**s,*cnt,*cflag;

  count = count_trials_ndata_condition_chan(nd,ndc,chan);
  s = (int **)myalloc(count*sizeof(int *));
  cnt = (int *)myalloc(count*sizeof(int));

  cflag = ndata_condition_get_flags(nd,ndc);

  k = 0;
  for(i=0;i<nd->ntrial;i++){
    j = ndata_trial_get_record_index(&(nd->t[i]),chan);
    if (j >= 0)
      if (cflag[i]){
	s[k] = nd->t[i].r[j].p;
	cnt[k] = nd->t[i].r[j].n;
	k += 1;
      }
  }
  myfree(cflag);
  *rs = s; *rcnt = cnt; *rn = count;
}
/**************************************-**************************************/
/*                                                                           */
/*                       GET_DOUBLE_SARRAY_NDATA_CONDITION                   */
/*                                                                           */
/*  This mallocs a list of pointers into the "nd" data structure, thus,      */
/*  only the list of pointers, not what they point to, should be freed.      */
/*                                                                           */
/*  Note:  "chan1" and "chan2" here are the ones from which the sarrays are  */
/*         taken, which may be different than the channel(s) used in the     */
/*         condition data structure, "ndc".                                  */
/*                                                                           */
/*  **** THIS IS A HACK - WYETH                                              */
/*  **** WYETH - NEVER USED.                                                 */
/*                                                                           */
/*****************************************************************************/
void get_double_sarray_ndata_condition(nd,ndc,chan1,chan2,rs1,rcnt1,rs2,rcnt2,
				       rn)
     struct ndata_struct *nd;
     struct ndcond_struct *ndc;
     char *chan1,*chan2;
     int ***rs1,**rcnt1,***rs2,**rcnt2,*rn;
{
  int i,j,k;
  int count,**s1,*cnt1,**s2,*cnt2,*cflag;

  printf("  *** WARNING: Wyeth, conditions applied only to chan1\n");
  exit_error("GET_DOUBLE_SARRAY_NDATA_CONDITION","Never used");

  count = count_trials_ndata_condition_chan(nd,ndc,chan1);
  s1 = (int **)myalloc(count*sizeof(int *));
  cnt1 = (int *)myalloc(count*sizeof(int));
  s2 = (int **)myalloc(count*sizeof(int *));
  cnt2 = (int *)myalloc(count*sizeof(int));

  cflag = ndata_condition_get_flags(nd,ndc);

  k = 0;
  for(i=0;i<nd->ntrial;i++){
    j = ndata_trial_get_record_index(&(nd->t[i]),chan1);
    if (j >= 0){
      if (cflag[i]){
	s1[k] = nd->t[i].r[j].p;
	cnt1[k] = nd->t[i].r[j].n;
	
	j = ndata_trial_get_record_index(&(nd->t[i]),chan1);
	if (j >= 0){
	  if (cflag[i]){ /************ WYETH, IS THIS CORRECT???? ******/
	    s2[k] = nd->t[i].r[j].p;
	    cnt2[k] = nd->t[i].r[j].n;
	    k += 1;
	  }else
	    exit_error("GET_DOUBLE_SARRAY_NDATA_CONDITION","Channel error 1");
	}else
	  exit_error("GET_DOUBLE_SARRAY_NDATA_CONDITION","Channel error 2");
      }
    }
  }
  myfree(cflag);
  *rs1 = s1; *rcnt1 = cnt1; *rs2 = s2; *rcnt2 = cnt2; *rn = count;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_MAX_PERIOD_RECORD                           */
/*                                                                           */
/*  Return the maximum period for the specified record.                      */
/*                                                                           */
/*****************************************************************************/
int get_max_period_record(ndrx)
     struct ndrx_struct *ndrx;
{
  int i;
  int max;

  max = -1;
  for(i=0;i<ndrx->n;i++)
    if (max < ndrx->r[i]->tn)
      max = ndrx->r[i]->tn;

  return max;
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_GROUP_IS_PARAM_CONST                        */
/*                                                                           */
/*  Return 1 or 2 if the parameter is constant across all trials.            */
/*  Return 0 if the parameter varies across trials.                          */
/*  Return -1 if there is any other type of exception.                       */
/*                                                                           */
/*  Return values                                                            */
/*    2 - parameter is 'const' in this nd file                               */
/*    1 - parameter is 'var' but value is const in this group                */
/*    0 - parameter is 'var' and value varies in this group                  */
/*   -1 - parameter could not be found, ERROR                                */
/*   -3 - there are no trials in this group                                  */
/*                                                                           */
/*****************************************************************************/
int ndata_group_is_param_const(nd,ndg,k,pname,xflag)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k;
     char *pname;
     int xflag;     // 1-exit if not constant, 0-do not exit
{
  int i,j;
  int n,flag,cflag;
  char *pstr,*tval;

  cflag = 1;  // Assume param will be constant

  n = ndg->cnt[k];  // Number of trials in group
  if (n <= 0){
    flag = -3;
  }else{

    flag = ndata_get_const_param(nd,pname,&tval);
    if (flag == 1)
      cflag = 2;  // This is a constant parameter
    else{

      // Get target value of parameter from first trial in group
      j = ndg->tnum[k][0];  // index of 1st trial in group
      flag = ndata_get_var_or_const_param(nd,j,pname,&tval);
      if (flag != 1)
	cflag = -1;

      for(i=1;i<n;i++){
	j = ndg->tnum[k][i]; //   get trial number

	flag = ndata_get_var_or_const_param(nd,j,pname,&pstr);
	if (flag != 1)
	  cflag = -1;
	else{
	  if (strcmp(pstr,tval)!=0)
	    cflag = 0;  // Parameter varies across trials
	  myfree(pstr);
	}
      }
    }
  }

  if ((xflag == 1) && (cflag < 1)){
    printf("  *** cflag = %d\n",cflag);
    exit_error("NDATA_GROUP_IS_PARAM_CONST","Parameter varies");
  }

  return cflag;
}
/**************************************-**************************************/
/*                                                                           */
/*                           NDATA_GROUP_GET_PAR_INT                         */
/*                                                                           */
/*  Return the unique value of the param for this group, or exit if a        */
/*  unique value does not exist.                                             */
/*                                                                           */
/*****************************************************************************/
int ndata_group_get_par_int(nd,ndg,k,pname)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k;
     char *pname;
{
  int i;
  int ival,cflag;

  cflag = ndata_group_is_param_const(nd,ndg,k,pname,1);  // 1-exit
  // returned flag value must be 1 or 2, otherwise it would have exited

  i = ndg->tnum[k][0];  // index of 1st trial in group
  ival = ndata_get_vc_int_exit(nd,i,pname,"NDATA_GROUP_GET_PAR_INT");

  return ival;
}
/**************************************-**************************************/
/*                                                                           */
/*                          NDATA_COUNT_RECORD_IN_GROUP                      */
/*                                                                           */
/*  Count the number of times the record "chan" appears in the group.        */
/*                                                                           */
/*****************************************************************************/
int ndata_count_record_in_group(nd,ndg,k,chan)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k;
     char *chan;
{
  int i,j;
  int count;

  count = 0;
  for(i=0;i<ndg->cnt[k];i++){
    j = ndata_trial_get_record_index(&(nd->t[ndg->tnum[k][i]]),chan);
    if (j >= 0)
      count += 1;
  }
  return count;
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_SEARCH_GROUPS_FOR_TRIAL                     */
/*                                                                           */
/*  Return the index of the group which contains the 'k'th trial.  If the    */
/*  trial does not occur in any group, return -1.                            */
/*                                                                           */
/*****************************************************************************/
int ndata_search_groups_for_trial(ndg,k)
     struct ndgroup_struct *ndg;
     int k;
{
  int i,j;

  for(i=0;i<ndg->n;i++)
    for(j=0;j<ndg->cnt[i];j++)
      if (ndg->tnum[i][j] == k)
	return i;

  return -1;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_FARRAY_NDATA_GROUP                          */
/*                                                                           */
/*  This mallocs a list of pointers into the "nd" data structure, thus,      */
/*  only the list of pointers, not what they point to, should be freed.      */
/*                                                                           */
/*****************************************************************************/
void get_farray_ndata_group(nd,ndg,k,chan,rsampling,rdata,rm,rn,rnames)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k;
     char *chan;
     float *rsampling,***rdata;
     int *rm;                       // Number of trials in group
     int *rn;                       // Length of data
     char ***rnames;
{
  int i,j;
  float **d,sampling;
  int m,n,tn,cnum,n0;
  char **names,temp[SLEN];

  int dflag = 0;

  // WYETH - currently, all records must have same data length, otherwise exit

  //
  //  WYETH DEBUG Aug 18, 2017 **** I realized that the next line is
  //                             *VERY SLOW* and replaced it.
  //
  //sampling = ndata_get_sampling_record_name(nd,chan);
  //if (sampling == -1)
  //exit_error("GET_FARRAY_NDATA_GROUP","Sampling varies across group");

  n0 = -1;
  n = -1; // to avoid compiler warning

/*
  if (dflag == 1){
    m = 1;
  }else{
*/
    m = ndata_count_record_in_group(nd,ndg,k,chan);
    if (m != ndg->cnt[k])
      exit_error("GET_FARRAY_NDATA_GROUP","Channel not found on all trials");
/*  } */

  sampling = 0.0;

  d = (float **)myalloc(m*sizeof(float *));
  names = (char **)myalloc(m*sizeof(char *));
  tn = 0;
  for(i=0;i<ndg->cnt[k];i++){
    j = ndg->tnum[k][i];
    cnum = ndata_trial_get_record_index(&(nd->t[j]),chan);
    if (cnum >= 0){

      if (sampling == 0.0){
	sampling = nd->t[j].r[cnum].sampling;
      }else{
	if (sampling != nd->t[j].r[cnum].sampling)
	  exit_error("GET_FARRAY_NDATA_GROUP","Sampling varies across group");
      }

/*      if (dflag == 1){
	d[tn] = (float *)myalloc(1*sizeof(float));
	d[tn][0] = 1.0;
	names[0] = strdup("test3");
	n = 1;
      }else{
*/
	d[tn] = nd->t[j].r[cnum].x;
	//printf("n = %d  tn = %d\n",nd->t[j].r[cnum].n,nd->t[j].r[cnum].tn);
	n = nd->t[j].r[cnum].n; // or tn?
	if (n0 == -1) // If this is the first trial, save the length 'n'
	  n0 = n;
	if (n != n0){
	  exit_error("GET_FARRAY_NDATA_GROUP","Variable data length, not imp'd");
	}
	sprintf(temp,"%s_%d",ndg->name[k],j);
	names[tn] = strdup(temp);
/*      } */
      tn += 1;
    }
  }

  *rsampling = sampling; *rdata = d; *rm = m; *rn = n; *rnames = names;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_SARRAY_NDATA_GROUP                          */
/*                                                                           */
/*  This mallocs a list of pointers into the "nd" data structure, thus,      */
/*  only the list of pointers, not what they point to, should be freed.      */
/*                                                                           */
/*****************************************************************************/
void get_sarray_ndata_group(nd,ndg,k,chan,rs,rcnt,rn,rnames)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k;
     char *chan;
     int ***rs,**rcnt,*rn;
     char ***rnames;
{
  int i,j;
  int **s,*cnt,n,tn,cnum;
  char **names,temp[SLEN];

  n = ndata_count_record_in_group(nd,ndg,k,chan);
  s = (int **)myalloc(n*sizeof(int *));
  cnt = (int *)myalloc(n*sizeof(int));
  names = (char **)myalloc(n*sizeof(char *));
  tn = 0;
  for(i=0;i<ndg->cnt[k];i++){
    j = ndg->tnum[k][i];
    cnum = ndata_trial_get_record_index(&(nd->t[j]),chan);
    if (cnum >= 0){
      if (nd->t[j].r[cnum].rtype != 0)
	exit_error("GET_SARRAY_NDATA_GROUP","Record type is not spikes");
      s[tn] = nd->t[j].r[cnum].p;
      cnt[tn] = nd->t[j].r[cnum].n;
      sprintf(temp,"%s_%d",ndg->name[k],j);
      names[tn] = strdup(temp);
      tn += 1;
    }
  }
  *rs = s; *rcnt = cnt; *rn = n; *rnames = names;
}
/**************************************-**************************************/
/*                                                                           */
/*                        GET_SARRAY_EM_RASTER_NDATA_GROUP                   */
/*                                                                           */
/*  Create and return an sarray that corresponds to eye movements for the    */
/*  given group.  Storage is created here, and should be freed.              */
/*                                                                           */
/*  This routine assumes the channel names "heye0" and "veye0" for the       */
/*  horizontal and vertical eye position records.                            */
/*                                                                           */
/*****************************************************************************/
void get_sarray_em_raster_ndata_group(nd,ndg,k,chan,sampling,sigma,crit,
				      min_gap,schan,sval,rs,rcnt,rn,rnames)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k;
     char *chan;
     float sampling,sigma,crit;
     int min_gap;
     char *schan;
     int sval,***rs,**rcnt,*rn;
     char ***rnames;
{
  int i,j,l;
  int **s,*cnt,n,nn,nv,*idata,offset;
  float *deriv,*tf,*xdata,*ydata,analog_sampling,tsec,ts_eye;
  char **names,temp[SLEN];
  struct ndtrial_struct *t;

  n = ndg->cnt[k];
  s = (int **)myalloc(n*sizeof(int *));
  cnt = (int *)myalloc(n*sizeof(int));
  names = (char **)myalloc(n*sizeof(char *));
  for(i=0;i<n;i++){
    j = ndg->tnum[k][i];
    t = &(nd->t[j]);
    xdata = ndata_trial_getp_adata_chan(t,"heye0",&nn);
    ydata = ndata_trial_getp_adata_chan(t,"veye0",&nv);
    if (nn != nv)
      exit_error("GET_SARRAY_EM_RASTER_NDATA_GROUP","nHoriz != nVert");

    deriv = get_smooth_abs_deriv_xy_data(xdata,ydata,nn,sigma);
    tf = get_binary_threshold_farray(deriv,nn,crit);
    idata = get_round_to_int_farray(tf,nn);
    myfree(tf);
    replace_min_runs_middle_iarray(idata,nn,0,1,min_gap); /* Fill short gaps */
    get_times_from_expanded_int_spikes(idata,nn,0,&(s[i]),&(cnt[i]));

    analog_sampling =  ndata_trial_get_sampling_chan(t,"heye0");
    for(l=0;l<cnt[i];l++)
      s[i][l] = my_rint((float)s[i][l] * (sampling/analog_sampling));

    if (schan != NULL){ /*** Adjust start for spike times. ***/
      tsec = ndata_trial_get_start_sec_schan(t,schan,sval);
      ts_eye = ndata_trial_get_start_sec_chan(t,"heye0");
      offset = my_rint((ts_eye - tsec)*sampling);
      for(l=0;l<cnt[i];l++)
	s[i][l] += offset;
    }

    sprintf(temp,"%s_%d",ndg->name[k],j);
    names[i] = strdup(temp);
  }
  *rs = s; *rcnt = cnt; *rn = n; *rnames = names;
}
/**************************************-**************************************/
/*                                                                           */
/*                            NDATA_GET_ADATA_GROUP                          */
/*                                                                           */
/*  Create and return an sarray that corresponds to eye movements for the    */
/*  given group.  Storage is created here, and should be freed.              */
/*                                                                           */
/*  This routine assumes the channel names "heye0" and "veye0" for the       */
/*  horizontal and vertical eye position records.                            */
/*                                                                           */
/*****************************************************************************/
void ndata_get_adata_group(nd,ndg,k,chan,schan,sval,sigma,rsampling,rdata,rcnt,
			   rn,rnames)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k;
     char *chan,*schan;
     int sval;
     float sigma,*rsampling,***rdata;
     int **rcnt,*rn;
     char ***rnames;
{
  int i,j;
  int *cnt,n,nn,offset;
  float *tdata,*tsdata,**data,analog_sampling,tsec,tsec_a;
  char **names,temp[SLEN];
  struct ndtrial_struct *t;

  analog_sampling = -1;
  n = ndg->cnt[k];
  data = (float **)myalloc(n*sizeof(float *));
  cnt = (int *)myalloc(n*sizeof(int));
  names = (char **)myalloc(n*sizeof(char *));
  for(i=0;i<n;i++){
    j = ndg->tnum[k][i];
    t = &(nd->t[j]);
    tdata = ndata_trial_getp_adata_chan(t,chan,&nn);
    tsdata = smooth_with_gaussian(tdata,nn,sigma,0.01);
    analog_sampling =  ndata_trial_get_sampling_chan(t,chan);

    if (schan != NULL){ /*** Adjust start times. ***/
      tsec = ndata_trial_get_start_sec_schan(t,schan,sval);
      tsec_a = ndata_trial_get_start_sec_chan(t,chan);
      offset = my_rint((tsec - tsec_a)*analog_sampling);
      if ((offset < 0)||(offset > nn))
	exit_error("NDATA_GET_ADATA_GROUP","Start is beyond analog data");
      data[i] = copy_farray(tsdata+offset,nn-offset);
      cnt[i] = nn-offset;
      myfree(tsdata);
    }else{
      data[i] = tsdata;
      cnt[i] = nn;
    }

    sprintf(temp,"%s_%d",ndg->name[k],j);
    names[i] = strdup(temp);
  }
  *rsampling = analog_sampling; *rdata = data; *rcnt = cnt; *rn = n;
  *rnames = names;
}
/**************************************-**************************************/
/*                                                                           */
/*                              PRINT_NDATA_INDEX                            */
/*                                                                           */
/*  Mainly for debugging.                                                    */
/*                                                                           */
/*****************************************************************************/
void print_ndata_index(nd,ndc,ndx)
     struct ndata_struct *nd;
     struct ndcond_struct *ndc;
     struct ndindex_struct *ndx;
{
  int i,j;
  int n,nrec;

  printf("  PRINT_NDATA_INDEX\n");

  if (ndx==NULL)
    printf("  Null index.\n");
  else{
    nrec = ndx->nrec;
    n = ndx->n;
    printf("    nrec = %d\n",nrec);
    printf("    nparam = %d\n",n);
    for(i=0;i<n;i++)
      if (ndx->ptype[i] == 'c')
	printf("      %s %c (%d unique values)\n",ndx->pname[i],ndx->ptype[i],
	       ndx->nuval[i]);
      else
	printf("      %s %c\n",ndx->pname[i],ndx->ptype[i]);
    for(i=0;i<n;i++)
      if (ndx->ptype[i] == 'c'){
	printf("    Unique values for 'c' param:\n     ");
	for(j=0;j<ndx->nuval[i];j++)
	  printf(" %s",ndx->uval[i][ndx->nuval[i]]); /*** WYETH -untested? ***/
	printf("\n");
      }
    if (ndx->ndxval != NULL){
      printf("    Index values:\n");
      for(i=0;i<nrec;i++){
	if (n > 0){
	  printf("(%.2f",ndx->ndxval[i][0]);
	  for(j=1;j<n;j++)
	    printf(",%.2f",ndx->ndxval[i][j]);
	  printf(") ");
	}
      }
      printf("\n");
    }
    if (ndx->tnum != NULL){
      printf("    Trial numbers:\n");
      for(i=0;i<nrec;i++)
	printf(" %d",ndx->tnum[i]);
      printf("\n");
    }
    if (ndx->tnum != NULL){
      printf("    Index pointers:\n");
      for(i=0;i<nrec;i++)
	printf(" %d",ndx->ndx[i]);
      printf("\n");
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MAKE_NDATA_INDEX                             */
/*                                                                           */
/*  The following must be assigned before calling:                           */
/*    ndx->n                                                                 */
/*    ndx->pname                                                             */
/*    ndx->ptype                                                             */
/*                                                                           */
/*****************************************************************************/
void make_ndata_index(nd,ndc,ndx)
     struct ndata_struct *nd;
     struct ndcond_struct *ndc;
     struct ndindex_struct *ndx;
{
  int i,j,k;
  int t,flag,done,bottom,*cflag;
  char *rstring;

  printf("  MAKE_NDATA_INDEX\n");

  cflag = ndata_condition_get_flags(nd,ndc);
  ndx->nrec = count_non_zero_iarray(cflag,nd->ntrial);
  //ndx->nrec = count_trials_ndata_condition(nd,ndc);
  
  if (ndc->n > 0)
    printf("    %d trials passed condition tests\n",ndx->nrec);

  ndx->ndxval = get_2d_farray(ndx->nrec,ndx->n);
  ndx->ndx = (int *)myalloc(ndx->nrec*sizeof(int));
  ndx->tnum = (int *)myalloc(ndx->nrec*sizeof(int));

  // Get unique value lists for string (category) parameters
  ndx->nuval = (int *)myalloc(ndx->n*sizeof(int));
  ndx->uval = (char ***)myalloc(ndx->n*sizeof(char **));
  for(i=0;i<ndx->n;i++){ // For each index parameter
    if (ndx->ptype[i] == 'c'){
      ndata_get_unique_var_values(nd,ndx->pname[i],&(ndx->uval[i]),
				  &(ndx->nuval[i]));
    }
  }

  // Compute float values for all index params for all trials
  k = 0;
  for(i=0;i<nd->ntrial;i++){ // For each trial
    if (cflag[i]){
      ndx->ndx[k] = k; // Index pointer
      ndx->tnum[k] = i; // Store this trial number in the index list
      for(j=0;j<ndx->n;j++){ // For each index parameter
	flag = ndata_get_var_value_trial(nd,i,ndx->pname[j],&rstring);
	if (flag == 0)
	  exit_error("MAKE_NDATA_INDEX","Parameter value not found");
	if ((ndx->ptype[j] == 'f')||(ndx->ptype[j] == 'i')){
	  //printf("k=%d i=%d rstring=%s %.2f\n",k,i,rstring,atof(rstring));
	  ndx->ndxval[k][j] = atof(rstring);
	}else if (ndx->ptype[j] == 'c')
	  ndx->ndxval[k][j] = (float)search_2d_carray(ndx->uval[j],rstring,
						      ndx->nuval[j]);
	else
	  exit_error("MAKE_NDATA_INDEX","Unknown parameter type");
	myfree(rstring);
      }
      k += 1;
    }
  }
  myfree(cflag);

  /***
    exit(0);
    for(i=0;i<k;i++){
    for(j=0;j<ndx->n;j++){
    printf("%.2f ",ndx->ndxval[i][j]);
    }
    printf("\n");
    }
    exit(0);***/

  /*** THE FOLLOWING IS A BUBBLE SORT ***/
  /*** The index is "ndx->ndx" and it had "ndx->nrec" entries in it. ***/
  /*** Each entry points to a trial in the "nd" structure. ***/
  printf("    Bubble sorting.\n");
  bottom = ndx->nrec - 1;
  done = 0;
  while (!done){
    done = 1;
    for (i=0;i<bottom;i++){

      // Decide whether we need to swap.
      flag = 0; // No movement needed
      k = 0;
      while((k < ndx->n)&&(flag==0)){
	if (ndx->ndxval[ndx->ndx[i]][k] > ndx->ndxval[ndx->ndx[i+1]][k]){
	  flag = 1; // Must swap this record
	}else if (ndx->ndxval[ndx->ndx[i]][k]== ndx->ndxval[ndx->ndx[i+1]][k]){
	  k += 1; // Compare the next parameters
	}else{
	  k = ndx->n; // This will force us out of the loop
	}
      }

      // Perform swap if necessary
      if (flag){
	done = 0;
	t = ndx->ndx[i];
	ndx->ndx[i] = ndx->ndx[i+1];
	ndx->ndx[i+1] = t;
      }
    }
    bottom -= 1;
  }
  //print_ndata_index(nd,ndc,ndx);  // For debugging
}
/**************************************-**************************************/
/*                                                                           */
/*                              FREE_NDATA_INDEX                             */
/*                                                                           */
/*****************************************************************************/
void free_ndata_index(ndx)
     struct ndindex_struct *ndx;
{
  int i;

  myfree(ndx->ndx);
  myfree(ndx->tnum);
  free_2d_farray(ndx->ndxval,ndx->nrec);
  for(i=0;i<ndx->n;i++)
    if (ndx->ptype[i] == 'c')
      free_2d_carray(ndx->uval[i],ndx->nuval[i]);
  myfree(ndx->uval);
  myfree(ndx->nuval);
  myfree(ndx->ptype);
  free_2d_carray(ndx->pname,ndx->n);
}
/**************************************-**************************************/
/*                                                                           */
/*                          COUNT_RECORD_IN_NDATA_INDEX                      */
/*                                                                           */
/*  Count the number of times the record "chan" appears in the group.        */
/*                                                                           */
/*****************************************************************************/
int count_record_in_ndata_index(nd,ndx,chan)
     struct ndata_struct *nd;
     struct ndindex_struct *ndx;
     char *chan;
{
  int i,j,k;
  int count;

  count = 0;
  for(i=0;i<ndx->nrec;i++){
    j = ndx->tnum[ndx->ndx[i]];
    k = ndata_trial_get_record_index(&(nd->t[j]),chan);
    if (k >= 0)
      count += 1;
  }
  return count;
}
/**************************************-**************************************/
/*                                                                           */
/*                           GET_SARRAY_NDATA_INDEX                          */
/*                                                                           */
/*  This mallocs a list of pointers into the "nd" data structure, thus,      */
/*  only the list of pointers, not what they point to, should be freed.      */
/*                                                                           */
/*****************************************************************************/
void get_sarray_ndata_index(nd,ndx,chan,rs,rcnt,rn)
     struct ndata_struct *nd;
     struct ndindex_struct *ndx;
     char *chan;
     int ***rs,**rcnt,*rn;
{
  int i,j,k;
  int cnum,count,**s,*cnt;

  count = count_record_in_ndata_index(nd,ndx,chan);
  s = (int **)myalloc(count*sizeof(int *));
  cnt = (int *)myalloc(count*sizeof(int));
  k = 0;
  for(i=0;i<ndx->nrec;i++){
    j = ndx->tnum[ndx->ndx[i]];
    cnum = ndata_trial_get_record_index(&(nd->t[j]),chan);
    if (cnum >= 0){
      s[k] = nd->t[j].r[cnum].p;
      cnt[k] = nd->t[j].r[cnum].n;
      k += 1;
    }
  }
  *rs = s; *rcnt = cnt; *rn = count;
}
/**************************************-**************************************/
/*                                                                           */
/*                      NDATA_GET_ADATA_POINTER_NDATA_INDEX                  */
/*                                                                           */
/*  Return *pointers* to the analog (float) data for the given channel       */
/*  "chan" according to the given index "ndx".                               */
/*  This mallocs a list of pointers into the "nd" data structure, thus,      */
/*  only the list of pointers, not what they point to, should be freed.      */
/*                                                                           */
/*****************************************************************************/
void ndata_get_adata_pointer_ndata_index(nd,ndx,chan,rdata,rcnt,rn)
     struct ndata_struct *nd;
     struct ndindex_struct *ndx;
     char *chan;
     float ***rdata;
     int **rcnt,*rn;
{
  int i,j,k;
  int cnum,count,*cnt;
  float **data;

  count = count_record_in_ndata_index(nd,ndx,chan);
  data = (float **)myalloc(count*sizeof(float *));
  cnt = (int *)myalloc(count*sizeof(int));
  k = 0;
  for(i=0;i<ndx->nrec;i++){
    j = ndx->tnum[ndx->ndx[i]];
    cnum = ndata_trial_get_record_index(&(nd->t[j]),chan);
    if (cnum >= 0){
      data[k] = nd->t[j].r[cnum].x;
      cnt[k] = nd->t[j].r[cnum].tn;
      k += 1;
    }
  }
  *rdata = data; *rcnt = cnt; *rn = count;
}
/**************************************-**************************************/
/*                                                                           */
/*                     NDATA_GET_ADATA_ADJUSTED_NDATA_INDEX                  */
/*                                                                           */
/*  Return analog data that has been adjusted for the start time and the     */
/*  period (in seconds).                                                     */
/*                                                                           */
/*****************************************************************************/
void ndata_get_adata_adjusted_ndata_index(nd,ndx,chan,schan,sval,start_sec,
					  period_sec,rdata,rm,rn)
     struct ndata_struct *nd;
     struct ndindex_struct *ndx;
     char *chan,*schan;
     int sval;
     float start_sec,period_sec,***rdata;
     int *rm,*rn;
{
  int i,j;
  int count,nn,n,offset;
  float **data,*tdata,analog_sampling,tsec,tsec_a;
  struct ndtrial_struct *t;

  /*** WYETH - this routine modified on Wed Apr 28 23:42:35 EDT 1999 ***/
  /*** "start_sec" was added, and changes were made to the "offset" ***/
  /*** computation for the case when "schan == NULL".  ***/

  n = -1;
  count = count_record_in_ndata_index(nd,ndx,chan);
  data = (float **)myalloc(count*sizeof(float *));
  for(i=0;i<ndx->nrec;i++){
    j = ndx->tnum[ndx->ndx[i]];
    t = &(nd->t[j]);
    tdata = ndata_trial_getp_adata_chan(t,chan,&nn); /* Exits if no "chan" */
    analog_sampling =  ndata_trial_get_sampling_chan(t,chan);
    if (i==0)
      n = my_rint(period_sec * analog_sampling);
    else if (n != my_rint(period_sec * analog_sampling))
      exit_error("NDATA_GET_ADATA_ADJUSTED_NDATA_INDEX","Sampling error");

    if (schan != NULL) /*** Adjust start times. ***/
      tsec = ndata_trial_get_start_sec_schan(t,schan,sval);
    else
      tsec = 0;
    tsec_a = ndata_trial_get_start_sec_chan(t,chan);
    offset = my_rint((tsec - tsec_a + start_sec)*analog_sampling);
    if ((offset < 0)||(offset > nn))
      exit_error("NDATA_GET_ADATA_ADJUSTED_NDATA_INDEX",
		 "Start is beyond analog data");
    if ((offset + n) > nn)
      exit_error("NDATA_GET_ADATA_ADJUSTED_NDATA_INDEX",
		 "End is beyond analog data");
    data[i] = copy_farray(tdata+offset,n);
  }
  *rdata = data; *rm = count; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_RECORD_ADJUST_START_TIME                    */
/*                                                                           */
/*  Adjust the record to have a start time "ts" seconds (following the       */
/*  trial TRef, which is not needed here.)                                   */
/*                                                                           */
/*****************************************************************************/
void ndata_record_adjust_start_time(pr,ts)
     struct ndrec_struct *pr;
     float ts;
{
  int i;
  int dt;
  float ts1;

  // WYETH - Old assumptions about start used here?

  ts1 = (float)(pr->t0)/pr->sampling;
  dt = my_rint((ts - ts1) * pr->sampling);
  if ((pr->rtype == 0)||(pr->rtype == 2)||(pr->rtype == 3)){
    pr->t0 += dt;
    for(i=0;i<pr->n;i++) /* Subtract from each point */
      pr->p[i] -= dt;
  }else{
    printf("  *** NDATA_RECORD_ADJUST_START_TIME  Warning, no adjustment to");
    printf("record of type 1.\n");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                      NDATA_TRIAL_GET_ADJUST_TIME_SEC                      */
/*                                                                           */
/*  Compute the appropriate 'start' time in seconds relative to 'tref' for   */
/*  this trial.  Return 1 if successful, 0 if failure.                       */
/*                                                                           */
/*****************************************************************************/
int ndata_trial_get_adjust_time_sec(t,schan,sval,rft)
     struct ndtrial_struct *t;
     char *schan;
     int sval;
     float *rft;
{
  int i,k;
  int rtype,et,flag;
  float ft;

  ft = 0;
  k = 0;
  i = ndata_trial_get_record_index(t,schan);
  if (i >= 0){
    rtype = t->r[i].rtype;
    if (rtype == 0){
      ft = (float)(t->r[i].p[sval-1] + t->r[i].t0)/t->r[i].sampling;
      k = 1;
    }else if (rtype == 3){
      flag = ndata_record_get_nth_occur_evco_time(&(t->r[i]),sval,1,&et);
      if (flag >= 0){
	ft = (float)(et + t->r[i].t0)/t->r[i].sampling;
	k = 1;
      }else{
	//printf("  sval = %d\n",sval);
	//exit_error("NDATA_TRIAL_GET_ADJUST_TIME_SEC","Wyeth ?");
	k = -1;  // No event found on trial
      }
    }else{
      exit_error("NDATA_TRIAL_GET_ADJUST_TIME_SEC","Invalid record type");
    }
  }
  *rft = ft;
  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                     NDATA_TRIAL_ADJUST_START_WITH_POINT                   */
/*                                                                           */
/*  Adjust the time of point record "rec0" by making the time of the "k"th   */
/*  point in record "rec1" be time zero.                                     */
/*                                                                           */
/*  One application of this routine is to set the start time to the time of  */
/*  the second sync pulse.                                                   */
/*                                                                           */
/*  NOTE:  The sum, t0+p[i] must remain constant, since it is the time from  */
/*         tref to the point.                                                */
/*                                                                           */
/*****************************************************************************/
void ndata_trial_adjust_start_with_point(pt,rec0,rec1,k)
     struct ndtrial_struct *pt;
     char rec0[],rec1[];
     int k;
{
  int r0,r1;
  float tstref;

  r0 = ndata_trial_get_record_index(pt,rec0);
  r1 = ndata_trial_get_record_index(pt,rec1);
  if ((r0>=0)&&(r1>=0)){
    tstref = (float)(pt->r[r1].p[k] + pt->r[r1].t0)/pt->r[r1].sampling;
    ndata_record_adjust_start_time(&(pt->r[r0]),tstref);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                    NDATA_TRIAL_ADJUST_START_WITH_EVENT                    */
/*                                                                           */
/*  Adjust the time of point record "rec0" so that the time of event         */
/*  "ecode" in the EVENT record "rec1" is time zero.                         */
/*                                                                           */
/*  One application of this routine is to set the start time to the time of  */
/*  a particular ecode, such as "stimon" etc.                                */
/*                                                                           */
/*  NOTE:  The sum, t0+p[i] must remain constant, since it is the time from  */
/*         tref to the point.                                                */
/*                                                                           */
/*****************************************************************************/
void ndata_trial_adjust_start_with_event(pt,rec0,rec1,ecode)
     struct ndtrial_struct *pt;
     char rec0[],rec1[];
     int ecode;
{
  int t,r0,r1,flag;
  float tstref;  /* Time in seconds from tref. */

  r0 = ndata_trial_get_record_index(pt,rec0);
  r1 = ndata_trial_get_record_index(pt,rec1);
  if ((r0>=0)&&(r1>=0)){
    /* Get time of event code. */
    flag = ndata_record_get_nth_occur_evco_time(&(pt->r[r1]),ecode,1,&t);
    if (flag < 0)
      exit_error("NDATA_ADJUST_START_EVENT_POINT_RECORD","ecode not found");
    tstref = (float)(t + pt->r[r1].t0)/pt->r[r1].sampling;
    ndata_record_adjust_start_time(&(pt->r[r0]),tstref);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          NDATA_ADJUST_START_TIME                          */
/*                                                                           */
/*  Adjust all of the times of "chan" so that some time value from another   */
/*  channel "schan" is taken to be time zero.  "k" may hold an ordinal       */
/*  number of a point, or the value of an event code that has an             */
/*  associated time.                                                         */
/*                                                                           */
/*  NOTE:  (Jan 2000) this is used by STA analyses to adjust spike and       */
/*  sync times.                                                              */
/*                                                                           */
/*****************************************************************************/
void ndata_adjust_start_time(nd,chan,schan,k)
     struct ndata_struct *nd;
     char *chan,*schan;
     int k;
{
  int i;
  int rtype1,rtype2;

  rtype1 = ndata_get_type_record_name(nd,schan);
  rtype2 = ndata_get_type_record_name(nd,chan);
  if (rtype1 == 0){ /*** This is a point channel, use the kth point ***/
    if (rtype2 == 0)
      for(i=0;i<nd->ntrial;i++)
	ndata_trial_adjust_start_with_point(&(nd->t[i]),chan,schan,k);
    else
      exit_error("NDATA_ADJUST_START_TIME","Not implemented yet.");
  }else if (rtype1 == 3){
    if (rtype2 == 0)
      for(i=0;i<nd->ntrial;i++)
	ndata_trial_adjust_start_with_event(&(nd->t[i]),chan,schan,k);
    else
      exit_error("NDATA_ADJUST_START_TIME","Not implemented yet.");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                     NDATA_ADJUST_START_TIME_NDATA_INDEX                   */
/*                                                                           */
/*  Adjust all of the times of "chan" so that some time value from another   */
/*  channel "schan" is taken to be time zero.  "k" may hold an ordinal       */
/*  number of a point, or the value of an event code that has an             */
/*  associated time.                                                         */
/*                                                                           */
/*  Adjust only those records in the given index.                            */
/*                                                                           */
/*  Jan 01 2000 - Most calls to this commented out - use global adjust.      */
/*                                                                           */
/*****************************************************************************/
void ndata_adjust_start_time_ndata_index(nd,chan,schan,k,ndx)
     struct ndata_struct *nd;
     char *chan,*schan;
     int k;
     struct ndindex_struct *ndx;
{
  int i;
  int rtype1,rtype2;
  struct ndtrial_struct *pt;


  rtype1 = ndata_get_type_record_name(nd,schan);
  rtype2 = ndata_get_type_record_name(nd,chan);
  if (rtype1 == 0){ /*** This is a point channel, use the kth point ***/
    if (rtype2 == 0)
      for(i=0;i<ndx->nrec;i++){
	pt = &(nd->t[ndx->tnum[i]]);
	ndata_trial_adjust_start_with_point(pt,chan,schan,k-1);
      }
    else
      exit_error("NDATA_ADJUST_START_TIME_NDATA_INDEX","Not implemented yet.");
  }else if (rtype1 == 3){
    if (rtype2 == 0){
      for(i=0;i<ndx->nrec;i++){
	pt = &(nd->t[ndx->tnum[i]]);
	ndata_trial_adjust_start_with_event(pt,chan,schan,k);
      }
    }else
      exit_error("NDATA_ADJUST_START_TIME_NDATA_INDEX","Not implemented yet.");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                        NDATA_ADJUST_START_TIME_GLOBAL                     */
/*                                                                           */
/*  Adjust start times of all channels (that can be suitably adjusted)       */
/*  for all trials.                                                          */
/*                                                                           */
/*  Jan 01 2000 - This routine is used in 'prepare_nda' instead of the old   */
/*                method with 'ndata_adjust_start_time_ndata_index'.         */
/*                                                                           */
/*****************************************************************************/
void ndata_adjust_start_time_global(nd,nda)
     struct ndata_struct *nd;
     struct nda_struct *nda;
{
  int i,j;
  int sval,no_event,flag,rtype;
  char *schan,*temp;
  float ft;
  struct ndtrial_struct *t;

  if (ndata_get_nda_param(nda,"start_chan",&schan)){
    if (ndata_get_nda_param(nda,"start_value",&temp)){

      if (is_number_string(temp) == 1){
	sval = atoi(temp);
	myfree(temp);
      }else{  // If 'sval' is a string, it should be an Event name

	rtype = ndata_get_type_record_name(nd,schan);
	if (rtype == 3){
	  if (nd->ntable > 1)
	    exit_error("NDATA_ADJUST_START_TIME_GLOBAL",
		       "NOT IMP'D YET for mutliple Event Tables\n");
	  //
	  //  *** WYETH THIS IS A HACK - assuming only 1 table
	  //
	  flag = ndata_lookup_ecode_num(&(nd->table[0]),temp,&sval);
	  if (flag != 1)
	    exit_error("NDATA_ADJUST_START_TIME_GLOBAL",
		       "Event name not found in table");
	  printf("    %s is code %d\n",temp,sval);
	}else{
	  printf("  rtype:  %d\n",rtype);
	  exit_error("NDATA_ADJUST_START_TIME_GLOBAL","Expecting Event chan");
	}
      }

      printf("  NDATA_ADJUST_START_TIME_GLOBAL\n");
      printf("    Note, only adjusts 'point' records, rtype = 0\n");

      no_event = 0;
      for(i=0;i<nd->ntrial;i++){ // For every trial
	//printf("trial %d\n",i);
	t = &(nd->t[i]);

	flag = ndata_trial_get_adjust_time_sec(t,schan,sval,&ft);
	if (flag == 1){  // Success
	  for(j=0;j<t->nrec;j++){ // Adjust all suitable records
	    if (t->r[j].rtype == 0)
	      ndata_record_adjust_start_time(&(t->r[j]),ft);
	  }
	}else if (flag == -1){  // Failed to find event in trial
	  no_event += 1;
	}else{
	  exit_error("NDATA_ADJUST_START_TIME_GLOBAL","Start time error");
	}
      }
      if (no_event > 0){
	printf("    *** %d of %d trials were not aligned: event not found.\n",
	       no_event,nd->ntrial);
      }
    }else
      exit_error("NDATA_ADJUST_START_TIME_GLOBAL",
		 "start_chan specified without start_value");
    myfree(schan);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                       NDATA_GET_START_TIME_NDATA_INDEX                    */
/*                                                                           */
/*  Return an array containing start times for each trial in the index.      */
/*  Each start time is determined by "schan" and "k".                        */
/*                                                                           */
/*  This could be rewritten as a simple loop using a function like           */
/*  "ndata_trial_get_start_sec_schan" that does not convert to seconds.      */
/*                                                                           */
/*****************************************************************************/
void ndata_get_start_time_ndata_index(nd,ndx,schan,k,rt)
     struct ndata_struct *nd;
     struct ndindex_struct *ndx;
     char *schan;
     int k,**rt;
{
  int i,j;
  int rtype1,*t,tt,flag;
  struct ndtrial_struct *pt;

  t = (int *)myalloc(ndx->nrec*sizeof(int));

  rtype1 = ndata_get_type_record_name(nd,schan);
  if (rtype1 == 0){ /*** This is a point channel, use the kth point ***/

    printf("*** WYETH k changed to k-1, to change from 0,1,... to 1,2...\n");
    printf("*** WYETH k changed to k-1, to change from 0,1,... to 1,2...\n");

    for(i=0;i<ndx->nrec;i++){
      pt = &(nd->t[ndx->tnum[i]]);
      j = ndata_trial_get_record_index(pt,schan);
      if (j>=0)
	t[i] = pt->r[j].p[k-1] + pt->r[j].t0;
      else
	exit_error("NDATA_GET_START_TIME_NDATA_INDEX","channel not found");
    }
  }else if (rtype1 == 3){
    for(i=0;i<ndx->nrec;i++){
      pt = &(nd->t[ndx->tnum[i]]);
      j = ndata_trial_get_record_index(pt,schan);
      if (j>=0){
	flag = ndata_record_get_nth_occur_evco_time(&(pt->r[j]),k,1,&tt);
	if (flag < 0)
	  exit_error("NDATA_GET_START_TIME_NDATA_INDEX","ecode not found");
	t[i] = tt + pt->r[j].t0;
      }else
	exit_error("NDATA_GET_START_TIME_NDATA_INDEX","channel not found");
    }
  }else
    exit_error("NDATA_GET_START_TIME_NDATA_INDEX","Record type error");

  *rt = t;
}
/**************************************-**************************************/
/*                                                                           */
/*                        NDATA_VAR_START_GET_OFFSET_LIST                    */
/*                                                                           */
/*  Look for specific offsets for the variable parameters.  The names have   */
/*  the following format:                                                    */
/*                                                                           */
/*    var_start_off_[var_param_value]  <value>                               */
/*                                                                           */
/*****************************************************************************/
void ndata_var_start_get_offset_list(nda,rname,roff,rn)
     struct nda_struct *nda;
     char ***rname;     // Return the list of 'names' (var vals as strings)
     float **roff;      // Return the list of offsets (seconds)
     int *rn;           // Lengths of lists returned
{
  int i,k;
  int n;
  float *voff;
  char *tc,**vname;
  struct param_pair *pp;

  n = 0;
  for(i=0;i<nda->nparam;i++){
    if (compare_prefix_string_order("var_start_off_",nda->pname[i]))
      n += 1;
  }

  vname = (char **)myalloc(n*sizeof(char *));
  voff = (float *)myalloc(n*sizeof(float));

  k = 0;
  for(i=0;i<nda->nparam;i++){
    tc = nda->pname[i];
    if (compare_prefix_string_order("var_start_off_",tc)){

      vname[k] = strdup(&tc[14]);   // Name starts at 14 position
      voff[k] = atof(nda->pval[i]); // offset is in seconds

      printf("    param val %s    offset %f (s)\n",vname[k],voff[k]);

      k += 1;
    }
  }

  *rname = vname;
  *roff  = voff;
  *rn    = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                              NDATA_VAR_START                              */
/*                                                                           */
/*  Adjust the start time of trials in the ndata file based on the value     */
/*  of a variable parameter.                                                 */
/*                                                                           */
/*****************************************************************************/
void ndata_var_start(nda,nd,vsp,vsadd,vsmult)
     struct nda_struct *nda;
     struct ndata_struct *nd;
     char *vsp;
     float vsadd;
     float vsmult;
{
  int i,j,k;
  int flag,dt;
  float fval,ts,ts1,toff;
  char *sval;
  struct ndtrial_struct *t;
  struct ndrec_struct *r;
  int voffn;        // Length of offset list
  char **voffname;  // List of values for which there is a special offset
  float *voffval;   // Offset values for variable param values in '

  ndata_var_start_get_offset_list(nda,&voffname,&voffval,&voffn);

  for(i=0;i<nd->ntrial;i++){

    flag = ndata_get_var_or_const_param(nd,i,vsp,&sval);
    if (flag == 0){
      //printf("sval = %s\n",sval);
      //printf("vsp = %s\n",sval);
      exit_error("PSTH_NDA","Start param not found");
    }
    k = search_2d_carray(voffname,sval,voffn);  // Search for offset
    if (k >= 0)
      toff = voffval[k];
    else
      toff = 0.0;

    flag = ndata_get_vc_param_float(nd,i,vsp,&fval);
    if (flag == 0)
      exit_error("PSTH_NDA","Start param not found");
    //printf("    start param %s has value %f (trial %d)\n",vsp,fval,i);
    ts = fval*vsmult + vsadd + toff; // New start in seconds
    //printf("    ndata_util.c:  new start value is %f (sec)\n",ts);

    t = &(nd->t[i]);
    t->tref += my_rint(ts * 1000.0);  // Update TREF (msec)

    for(j=0;j<t->nrec;j++){
      r = &(t->r[j]);

      dt = my_rint(ts * r->sampling); // Change (sampling units)
      if ((r->rtype == 0)||(r->rtype == 2)||(r->rtype == 3)){
	r->t0 -= dt;
	for(k=0;k<r->n;k++) /* Subtract from each point */
	  r->p[k] -= dt;
      }else{
	printf("  NDATA_VAR_START  Warning:  rec of type %d ignored\n",
	       r->rtype);
      }
    }
  }

  if (voffn > 0){
    free_2d_carray(voffname,voffn);
    myfree(voffval);
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                              NDU_VARDUR_TRIAL                             */
/*                                                                           */
/*  Get the duration for this trial.                                         */
/*                                                                           */
/*****************************************************************************/
int ndu_vardur_trial(tr,nda)
     struct ndtrial_struct *tr;
     struct nda_struct *nda;
{
  int dur,flag,ns,toff,tdiff;
  char *pstr,**slist,*chan;
  struct ndrec_struct *rp;
  
  flag = ndata_get_nda_param(nda,"vardur",&pstr);
  if (flag != 1)
    exit_error("NDU_VARDUR_TRIAL","Cannot find 'vardur' in nda script.");
  get_items_from_string(pstr,&slist,&ns);
  myfree(pstr);

  if (strcmp(slist[0],"events")==0){

    chan = slist[1];
    flag = ndu_trial_chan_ecode_str_diff(tr,chan,slist[2],slist[3],&tdiff);
    toff = atoi(slist[4]);

    if (flag != 1){
      printf("flag = %d\n",flag);
      exit_error("NDU_VARDUR_TRIAL","Error computing time difference");
    }

    dur = tdiff + toff;

    /*
    printf("  Channel:  %s\n",chan);
    printf("  Event 1:  %s\n",slist[2]);
    printf("  Event 2:  %s\n",slist[3]);
    printf("  Offset:   %d\n",toff);
    printf("  T Diff:   %d\n",tdiff);
    */

    free_2d_carray(slist,ns);
  }else{
    printf("  VarDur type:  %s\n",slist[0]);
    exit_error("NDU_VARDUR_TRIAL","Unknown vardur type");
  }

  return dur;
}
/**************************************-**************************************/
/*                                                                           */
/*                              NDU_VARDUR_GROUP                             */
/*                                                                           */
/*  Get the array of durations for the group.                                */
/*                                                                           */
/*****************************************************************************/
int *ndu_vardur_group(nd,ndg,k,nda)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     int k;                          // Group index
     struct nda_struct *nda;
{
  int i;
  int n,*vdur;
  struct ndtrial_struct *tr;

  n = ndg->cnt[k];  // Number of trials in group

  vdur = (int *)myalloc(n*sizeof(int));

  for(i=0;i<n;i++){
    tr = &(nd->t[ndg->tnum[k][i]]);
    vdur[i] = ndu_vardur_trial(tr,nda);
  }

  return vdur;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MAKE_NDATA_ANALYSIS                           */
/*                                                                           */
/*  Note, "ndx->pname" and "ndx->ptype" are not created for the "group 0"    */
/*  condition because they require information from the "nd" file.  The      */
/*  routine "set_ndata_analysis_index" is used for this purpose.             */
/*                                                                           */
/*****************************************************************************/
void make_ndata_analysis(atype,grouping,outtype,ngroup,pname,ptype,rnda)
     char atype[],grouping[],outtype[];
     int ngroup;
     char **pname,*ptype;
     struct nda_struct **rnda;
{
  int i;
  struct nda_struct *nda;

  printf("  MAKE_NDATA_ANALYSIS\n");

  nda = (struct nda_struct *)myalloc(sizeof(struct nda_struct));
  nda->atype = strdup(atype);
  nda->grouping = strdup(grouping);
  if (strcmp(nda->grouping,"group")==0){
    nda->ndx = (struct ndindex_struct *)myalloc(sizeof(struct ndindex_struct));
    nda->ndx->n = ngroup; // Number of grouping params
    if (nda->ndx->n > 0){
      nda->ndx->pname = (char **)myalloc(nda->ndx->n*sizeof(char *));
      nda->ndx->ptype = (char *)myalloc(nda->ndx->n*sizeof(char));
      for(i=0;i<nda->ndx->n;i++)
	nda->ndx->pname[i] = strdup(pname[i]);
    }
  }else
    nda->ndx = NULL;
  nda->outtype = strdup(outtype);
  nda->nparam = 0; // Parameters will be handled elsewhere

  nda->ndc = get_empty_condition(0); // No conditions

  *rnda = nda;
}
/**************************************-**************************************/
/*                                                                           */
/*                             NDU_LOOKUP_RESOLVE                            */
/*                                                                           */
/*  Format example:                                                          */
/*                                                                           */
/*    lookup <filename> <condition> <param> [<field>] <param> [<field>]...   */
/*                                                                           */
/*  Must look up each <field> name in the <filename> and replace it.         */
/*                                                                           */
/*****************************************************************************/
void ndu_lookup_resolve(cond,k,keyname)
     struct ndcond_struct *cond;   // Condition structure
     int k;                        // Index of condition to process
     char *keyname;                // Name of data file, key to lookup table
{
  int i,j;
  int n,ii,lcnt,tlen;
  char **tc,*filename,fieldname[SLEN],*pstr;

  printf("  NDU_LOOKUP_RESOLVE\n");

  n  = cond->ns[k];        // Number of strings in this condition line
  tc = cond->slist[k];     // Array os strings [n]

  if (strcmp(tc[0],"lookup")!=0)
    exit_error("NDU_LOOKUP_RESOLVE","Expecting 'lookup'");

  filename = tc[1];

  lcnt = 0;
  for(i=2;i<n;i++){  // Skip over "lookup <filename>"
    if (tc[i][0] == '['){
      tlen = strlen(tc[i]);
      if (tc[i][tlen-1] != ']')
	exit_error("NDU_LOOKUP_RESOLVE","No matching ']' in lookup");
      for(j=1;j<tlen-1;j++)
	fieldname[j-1] = tc[i][j];
      fieldname[tlen-2] = '\0';     // Terminating NULL char

      //printf("Fieldname = %s\n",fieldname);

      ii = lookup_key_field_file(filename,keyname,fieldname,&pstr);
      if (ii <= 0){
	printf("  look-up file:  %s\n",filename);
	printf("  look-up key:  %s\n",keyname);
	printf("  look-up field:  %s\n",fieldname);
	exit_error("NDU_LOOKUP_RESOLVE","Lookup failed");
      }

      myfree(tc[i]);  // Replace the string with the looked-up value
      tc[i] = pstr;

      lcnt += 1;
    }
  }
  if (lcnt == 0){
    printf("  *** No fields found in '[...]' notation on look-up line\n");
    exit_error("NDU_LOOKUP_RESOLVE","Please use '[...]' format");
  }

  //print_condition_list(cond);
  printf("    %d values replaced by lookup from table %s\n",lcnt,filename);

  //  Compress the first two strings out of the list
  myfree(tc[0]);
  myfree(tc[1]);
  cond->ns[k] = n-2;
  printf("    ");
  for(i=0;i<n-2;i++){
    tc[i] = tc[i+2];
    printf("%s ",tc[i]);
  }
  printf("\n");
  //cond->slist[i] = cond->slist[i+2];


  //print_condition_list(cond);
}
/**************************************-**************************************/
/*                                                                           */
/*                             READ_NDATA_ANALYSIS                           */
/*                                                                           */
/*  Note, "ndx->pname" and "ndx->ptype" are not created for the "group 0"    */
/*  condition because they require information from the "nd" file.  The      */
/*  routine "set_ndata_analysis_index" is used for this purpose.             */
/*                                                                           */
/*  Notes:                                                                   */
/*  - Parameter values can be a list of values, but are stored as 1 string.  */
/*  - Lookup params added 9/21/00.                                           */
/*                                                                           */
/*****************************************************************************/
void read_ndata_analysis(infile,name,rnda)
     char infile[];      // .nda file name
     char name[];        // Name (w/o path) of the data file
     struct nda_struct **rnda;
{
  FILE *fopen(),*fin;
  int i,k;
  int ns,pcount,ccount,ival,nnc,nn;
  float fval;
  struct nda_struct *nda;
  struct ndcond_struct *tc;
  char temp[SLEN],temp2[SLEN],**slist,*pstr,*tg;

  printf("  READ_NDATA_ANALYSIS\n");

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** Cannot open file %s.\n",infile);
    *rnda = NULL;
    return;
  }else
    printf("    Reading %s\n",infile);

  nda = (struct nda_struct *)myalloc(sizeof(struct nda_struct));
  nn = fscanf(fin,"%s",temp);
  while(temp[0] == '#'){
    tg = fgets(temp,SLEN,fin);
    nn = fscanf(fin,"%s",temp);
  }
  nda->atype = strdup(temp);
  nn = fscanf(fin,"%s",temp); nda->grouping = strdup(temp);
  if (strcmp(nda->grouping,"group")==0){
    nda->ndx = (struct ndindex_struct *)myalloc(sizeof(struct ndindex_struct));

    // Orig REMOVE
    //nn = fscanf(fin,"%d",&(nda->ndx->n)); // Number of grouping params

    nn = fscanf(fin,"%s",temp2);
    if (strcmp(temp2,"all")==0){
      printf("  *** 'group all' is an invalid grouping, please use one of:\n");
      printf("      'all' - to put all trials together in one group, or\n");
      printf("      'group 0' - to automatically create all unique groups\n");
      printf("      (see nData documentation for more options)\n");
      exit_error("READ_NDATA_ANALYSIS","Please correct the .nda file");
    }else if (is_int_string(temp2) == 0){
      printf("  *** Offending string:  %s\n",temp2);
      exit_error("READ_NDATA_ANALYSIS",
		 "The value after 'group' is not an integer (in .nda file)");
    }else
      nda->ndx->n = atoi(temp2);

    
    //nn = fscanf(fin,"%d",&(nda->ndx->n)); // Number of grouping params

    if (nda->ndx->n > 0){
      nda->ndx->pname = (char **)myalloc(nda->ndx->n*sizeof(char *));
      nda->ndx->ptype = (char *)myalloc(nda->ndx->n*sizeof(char));
      for(i=0;i<nda->ndx->n;i++){
	nn = fscanf(fin,"%s",temp);
	nda->ndx->pname[i] = strdup(temp);
      }
      nda->ndx->nuval = NULL; nda->ndx->uval = NULL; nda->ndx->ndxval = NULL;
      nda->ndx->tnum = NULL; nda->ndx->ndx = NULL;
    }
  }else{
    nda->ndx = (struct ndindex_struct *)myalloc(sizeof(struct ndindex_struct));
    nda->ndx->n = 0;
  }
  nn = fscanf(fin,"%s",temp); nda->outtype = strdup(temp);

  /***************************************************************/
  /*** WYETH - new code to get rid of the number of parameters ***/
  /*** and the number after 'condition'                        ***/
  /*** Thu Mar 30 15:16:52                                     ***/

  pcount = 0;
  ccount = -1;
  while(fscanf(fin,"%s",temp) != EOF){
    /*printf("temp ==>%s<==\n",temp);*/
    tg = fgets(temp2,SLEN,fin);
    get_items_from_string(temp2,&slist,&ns);

    if (temp[0] != '#'){ /* Ignore comments, beginning with '#' */
      if (strcmp(temp,"condition")==0)
	ccount = 0;
      else if (ns == 0){
	if (pcount == 0)
	  printf("  ignoring item %s\n",temp);
	else{
	  if (ccount != -1) /* Count conditions that are only one item */
	    ccount += 1;
	  else
	    exit_error("READ_NDATA_ANALYSIS","No values given for param");
	}
      }else{
	if (ccount == -1)
	  pcount += 1;
	else
	  ccount += 1;
      }
    }
    /*printf("-->%s<--\n",temp2);*/
  }
  printf("    %d parameters.\n",pcount);
  if (ccount < 0)
    printf("    No conditions.\n");
  else
    printf("    Conditions %d\n",ccount);
  
  fclose(fin);

  nda->nparam = pcount;
  nda->pname = (char **)myalloc(nda->nparam*sizeof(char *));
  nda->pval = (char **)myalloc(nda->nparam*sizeof(char *));

  if (ccount < 0)
    ccount = 0;

  nda->ndc = get_empty_condition(0);
  tc = nda->ndc;
  tc->n = ccount; // Number of conditions
  tc->ns = (int *)myalloc(tc->n*sizeof(int));
  tc->slist = (char ***)myalloc(tc->n*sizeof(char **));

  fin = fopen(infile,"r");

  /*** Read through the header, and ignore. ***/
  nn = fscanf(fin,"%s",temp);
  while(temp[0] == '#'){
    tg = fgets(temp,SLEN,fin);
    nn = fscanf(fin,"%s",temp);
  } // Last read was atype
  nn = fscanf(fin,"%s",temp); // grouping type
  if (strcmp(temp,"group")==0){
    nn = fscanf(fin,"%d",&k); // Number of grouping params
    if ((k > 0) && (k < 30)){  // WYETH - 30 is arbitrary limit
      for(i=0;i<k;i++)
	nn = fscanf(fin,"%s",temp);
    }else if (k != 0){
      printf("k = %d\n",k);
      exit_error("HERE","invalid number following 'group'");
    }
  }
  nn = fscanf(fin,"%s",temp); // outtype

  //  READ PARAMETER NAMES AND VALUES
  pcount = 0;
  while(pcount < nda->nparam){
    nn = fscanf(fin,"%s",temp); // Get param name
    tg = fgets(temp2,SLEN,fin); // Get param values
    if (temp[0] != '#'){   // To ignore lines beginning with '#'
      get_items_from_string(temp2,&slist,&ns);
      if (ns == 0){
	if (pcount != 0)
	  exit_error("READ_NDATA_ANALYSIS","No values given for param");
      }else{
	nda->pname[pcount] = strdup(temp);
	if (ns==1)
	  nda->pval[pcount] = strdup(slist[0]); // Store single value
	else{
	  // Store multiple values
	  nnc = ns;  // Number of non-comment values
	  i = 1;
	  while((i<ns)&&(nnc == ns)){
	    if (slist[i][0] == '#')
	      nnc = i;
	    else
	      i += 1;
	  }
	  nda->pval[pcount] = make_string_from_items(slist,nnc);
	}
	free_2d_carray(slist,ns);
	pcount += 1;
      }
    }
  }

  //  READ CONDITIONS
  ccount = 0;
  while(ccount < tc->n){
    if (fgets(temp,SLEN,fin)==NULL)
      exit_error("READ_NDATA_ANALYSIS","Cannot find condition string");
    get_items_from_string(temp,&(tc->slist[ccount]),&(tc->ns[ccount]));
    if (tc->ns[ccount] > 0){
      if (strcmp(tc->slist[ccount][0],"condition")==0)
	; /*printf("condition\n");*/
      else if (tc->slist[ccount][0][0] == '#')
	; /*printf("condition COMMENT\n");*/
      else{
	if (tc->ns[ccount] > 0)
	  ccount += 1;
      }
    }
  }

  /***
  printf("PARAMETERS\n");
  for(i=0;i<nda->nparam;i++){
    printf(" %s %s\n",nda->pname[i],nda->pval[i]);
  }  

  printf("CONDITIONS\n");
  for(i=0;i<tc->n;i++){
    for(k=0;k<tc->ns[i];k++)
      printf(" %s",tc->slist[i][k]);
    printf("\n");
  }
  ***/

  // Handle look-up parameters
  for(i=0;i<pcount;i++){
    if ((strcmp(nda->pname[i],"lookup")==0)||
	(strcmp(nda->pname[i],"lookup_float_op")==0)||
	(strcmp(nda->pname[i],"lookup_int_op")==0)){
      get_items_from_string(nda->pval[i],&slist,&ns);
      /*printf("LOOKUP %s in file %s as %s\n",slist[0],slist[1],slist[2]);*/
      k = lookup_key_field_file(slist[1],name,slist[2],&pstr);
      if (strcmp(nda->pname[i],"lookup_int_op")==0){
	ival = atoi(pstr);
	if (strcmp(slist[3],"+")==0)
	  ival += atoi(slist[4]);
	else if (strcmp(slist[3],"-")==0)
	  ival -= atoi(slist[4]);
	else if (strcmp(slist[3],"*")==0)
	  ival *= atoi(slist[4]);
	else if (strcmp(slist[3],"/")==0)
	  ival /= atoi(slist[4]);
	else
	  exit_error("READ_NDATA_ANALYSIS","Unknown operation");
	myfree(pstr);
	sprintf(temp,"%d",ival);
	pstr = strdup(temp);
      }else if (strcmp(nda->pname[i],"lookup_float_op")==0){
	fval = atof(pstr);
	if (strcmp(slist[3],"+")==0)
	  fval += atof(slist[4]);
	else if (strcmp(slist[3],"-")==0)
	  fval -= atof(slist[4]);
	else if (strcmp(slist[3],"*")==0)
	  fval *= atof(slist[4]);
	else if (strcmp(slist[3],"/")==0)
	  fval /= atof(slist[4]);
	else if (strcmp(slist[3],"pow")==0)
	  fval = pow(fval,atof(slist[4]));
	else
	  exit_error("READ_NDATA_ANALYSIS","Unknown operation");
	myfree(pstr);
	sprintf(temp,"%f",fval);
	pstr = strdup(temp);
      }
      
      myfree(nda->pname[i]);
      nda->pname[i] = strdup(slist[0]);
      myfree(nda->pval[i]);
      nda->pval[i] = pstr;
      printf("    Look-up:  %s is %s\n",slist[0],pstr);
      free_2d_carray(slist,ns);
      
      /*printf("k=%d\n",k);*/
      if (k <= 0)
	exit_error("READ_NDATA_ANALYSIS","Lookup value not found");
    }
  }

  //  Handle look-up conditions
  //                  0        1          2        3       4       5
  //  OLD EXAMPLE:  lookup <filename> param_range sper <pname1> <pname2>
  for(i=0;i<tc->n;i++){

    if (strcmp(tc->slist[i][0],"lookup")==0)
      ndu_lookup_resolve(tc,i,name);
    
    /*
    if (strcmp(tc->slist[i][0],"lookup")==0){
      if (strcmp(tc->slist[i][2],"param_range")==0){
	k = lookup_key_field_file(tc->slist[i][1],name,tc->slist[i][4],&pstr);
	if (k <= 0) exit_error("READ_NDATA_ANALYSIS",
			       "First lookup value not found");
	myfree(tc->slist[i][0]); // Free 'lookup'
	tc->slist[i][0] = tc->slist[i][2]; // Move condit. type to beginning
	tc->slist[i][2] = pstr;
	k = lookup_key_field_file(tc->slist[i][1],name,tc->slist[i][5],&pstr);
	if (k <= 0) exit_error("READ_NDATA_ANALYSIS",
			       "2nd lookup value not found");
	myfree(tc->slist[i][1]); // Free 'lookup'
	tc->slist[i][1] = tc->slist[i][3]; // Move condit. type to beginning
	tc->slist[i][3] = pstr;
	myfree(tc->slist[i][4]); // Free lookup field names
	myfree(tc->slist[i][5]);
	tc->ns[i] = 4;
      }else
	exit_error("READ_NDATA_ANALYSIS","Condit. type not impl'd for lookup");
    }
    */
    /*else printf("==>%s\n",tc->slist[i][0]);*/
  }

  fclose(fin);
  *rnda = nda;
}
/**************************************-**************************************/
/*                                                                           */
/*                          SET_NDATA_ANALYSIS_INDEX                         */
/*                                                                           */
/*  The analysis description in "nda" may specify grouping parameters        */
/*  explicitly by name or implicity with "group 0".  The types of these      */
/*  parameters, and their names, have to be placed in "ndx" struct.          */
/*                                                                           */
/*****************************************************************************/
void set_ndata_analysis_index(nd,nda)
     struct ndata_struct *nd;
     struct nda_struct *nda;
{
  int i;
  int flag;

  if (nda->ndx != NULL){
    if (strcmp(nda->grouping,"group")==0){
      if (nda->ndx->n == 0){ // Default condition: "group 0"
	nda->ndx->n = nd->nvar; // Use all var params
	nda->ndx->pname = (char **)myalloc(nda->ndx->n*sizeof(char *));
	nda->ndx->ptype = (char *)myalloc(nda->ndx->n*sizeof(char));
	for(i=0;i<nda->ndx->n;i++){
	  nda->ndx->pname[i] = nd->vname[i];
	  nda->ndx->ptype[i] = nd->vtype[i];
	}
      }else
	for(i=0;i<nda->ndx->n;i++){ // Copy type fields from "nd" header
	  flag = ndata_get_var_type(nd,nda->ndx->pname[i],
				    &(nda->ndx->ptype[i]));
	  if (flag == 0){
	    printf("Parameter name %s\n",nda->ndx->pname[i]);
	    exit_error("SET_NDATA_ANALYSIS_INDEX",
		       "Param name not found, or value does not vary");
	  }
	}
    }
    // Conditions are applied here
    make_ndata_index(nd,nda->ndc,nda->ndx); // for "group" or "all"
 }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MAKE_NAME_FOR_GROUP                           */
/*                                                                           */
/*  Note that "name" must already point to memory location.                  */
/*                                                                           */
/*  flag                                                                     */
/*     1 - use '=' between names and values                                  */
/*     2 - use '_' between names and values                                  */
/*                                                                           */
/*****************************************************************************/
void make_name_for_group(ndg,i,flag,infile,name)
     struct ndgroup_struct *ndg;
     int i,flag;
     char infile[],name[];
{
  int j;
  int tlen;

  if (flag == 1){
    if ((strcmp(ndg->name[i],"all")==0)&&(ndg->n == 1))
      get_name_without_path(infile,name);
    else
      sprintf(name,"%s",ndg->name[i]);
  }else if (flag == 2){
    if ((strcmp(ndg->name[i],"all")==0)&&(ndg->n == 1))
      get_name_without_path(infile,name);
    else{
      sprintf(name,"%s",ndg->name[i]);
      tlen = strlen(name);
      for(j=0;j<tlen;j++){
	if (name[j] == '=')
	  name[j] = '_';
      }
    }
  }else{
    exit_error("MAKE_NAME_FOR_GROUP","Unknown flag value");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          MAKE_SPACE_NAME_FOR_GROUP                        */
/*                                                                           */
/*  Note that "name" must already point to memory location.                  */
/*                                                                           */
/*  Name string has format:                                                  */
/*    "name1 10.0 name2 20.0 name3 30.0"                                     */
/*                                                                           */
/*****************************************************************************/
void make_space_name_for_group(ndg,i,flag,infile,name)
     struct ndgroup_struct *ndg;
     int i,flag;
     char infile[],name[];
{    
  if (flag == 1){
    if ((strcmp(ndg->name[i],"all")==0)&&(ndg->n == 1))
      get_name_without_path(infile,name);
    else{
      sprintf(name,"%s",ndg->name[i]);
      for(i=0;i<strlen(name);i++)
	if ((name[i]=='=')||(name[i]==','))
	  name[i] = ' ';
    }
  }else{
    exit_error("MAKE_SPACE_NAME_FOR_GROUP","Unknown flag value");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           NDATA_MAKE_NAME_VAR_TRIAL                       */
/*                                                                           */
/*  Create a name for the kth trial in the ndata file.  Create the name by   */
/*  using the names and values of the var parameters.                        */
/*                                                                           */
/*****************************************************************************/
char *ndata_make_name_var_trial(nd,k)
     struct ndata_struct *nd;
     int k;
{
  int i;
  struct ndtrial_struct *t;
  char temp[SLEN],tname[SLEN];

  t = &(nd->t[k]);
  tname[0] = '\0';
  for(i=0;i<t->nparam;i++){
    sprintf(temp,"%s=%s,",t->pname[i],t->pval[i]);
    strcat(tname,temp);
  }
  tname[strlen(tname)-1] = '\0'; // Remove the final comma

  return strdup(tname);
}
/**************************************-**************************************/
/*                                                                           */
/*                          NDATA_MAKE_NAME_INDEX_TRIAL                      */
/*                                                                           */
/*  Create a name for the kth trial in the index.  The name is created by    */
/*  concatenating all the parameter names and values in the index for that   */
/*  trial.                                                                   */
/*                                                                           */
/*****************************************************************************/
char *ndata_make_name_index_trial(nd,ndx,k)
     struct ndata_struct *nd;
     struct ndindex_struct *ndx;
     int k;
{
  int i;
  int flag,tk;
  char temp[SLEN],tname[SLEN],*sval;

  tk = ndx->tnum[ndx->ndx[k]]; // Get trial no. for kth entry in index.
  tname[0] = '\0';
  for(i=0;i<ndx->n;i++){
    flag = ndata_get_var_value_trial(nd,tk,ndx->pname[i],&sval);
    if (flag == 0)
      exit_error("NDATA_MAKE_NAME_INDEX_TRIAL","Parameter name not found");
    sprintf(temp,"%s=%s,",ndx->pname[i],sval);
    myfree(sval);
    strcat(tname,temp);
  }
  if (ndx->n > 0)
    tname[strlen(tname)-1] = '\0'; // Remove the final comma
  else
    strcpy(tname,"all");

  return strdup(tname);
}
/**************************************-**************************************/
/*                                                                           */
/*                      NDATA_MAKE_NAME_INDEX_TRIAL_SPACE                    */
/*                                                                           */
/*  Create a name for the kth trial in the index.  The name is created by    */
/*  concatenating all the parameter names and values in the index for that   */
/*  trial.                                                                   */
/*                                                                           */
/*****************************************************************************/
char *ndata_make_name_index_trial_space(nd,ndx,k)
     struct ndata_struct *nd;
     struct ndindex_struct *ndx;
     int k;
{
  int i;
  int flag,tk;
  char temp[SLEN],tname[SLEN],*sval;

  tk = ndx->tnum[ndx->ndx[k]]; // Get trial no. for kth entry in index.
  tname[0] = '\0';
  for(i=0;i<ndx->n;i++){
    flag = ndata_get_var_value_trial(nd,tk,ndx->pname[i],&sval);
    if (flag == 0)
      exit_error("NDATA_MAKE_NAME_INDEX_TRIAL","Parameter name not found");
    sprintf(temp,"%s %s ",ndx->pname[i],sval);
    myfree(sval);
    strcat(tname,temp);
  }
  if (ndx->n > 0)
    tname[strlen(tname)-1] = '\0'; // Remove the final space
  else
    strcpy(tname,"all");

  return strdup(tname);
}
/**************************************-**************************************/
/*                                                                           */
/*                             PRINT_GROUP_STRUCT                            */
/*                                                                           */
/*****************************************************************************/
void print_group_struct(ndg)
     struct ndgroup_struct *ndg;
{
  int i,j;

  printf("  PRINT_GROUP_STRUCT\n");

  for(i=0;i<ndg->n;i++){
    printf("    %s ",ndg->name[i]);
    if (ndg->value==NULL)
      printf("(null) ");
    else
      printf("%s ",ndg->value[i]);
    printf("%d - ",ndg->cnt[i]);
    for(j=0;j<ndg->cnt[i];j++)
      printf(" %d",ndg->tnum[i][j]);
    printf("\n");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                               NDATA_FREE_GROUP                            */
/*                                                                           */
/*****************************************************************************/
void ndata_free_group(ndg)
     struct ndgroup_struct *ndg;
{
  free_2d_carray(ndg->name,ndg->n);
  free_2d_carray(ndg->value,ndg->n);
  myfree(ndg->cnt);
  free_2d_iarray(ndg->tnum,ndg->n);
  myfree(ndg);
}
/**************************************-**************************************/
/*                                                                           */
/*                                GET_GROUP_ALL                              */
/*                                                                           */
/*  Make one group containing all trials satisfying the conditions.          */
/*                                                                           */
/*****************************************************************************/
void get_group_all(nd,ndc,rndg)
     struct ndata_struct *nd;
     struct ndcond_struct *ndc;
     struct ndgroup_struct **rndg;
{
  int i,k;
  int *cflag;
  struct ndgroup_struct *ndg;

  printf("  GET_GROUP_ALL\n");

  cflag = ndata_condition_get_flags(nd,ndc);

  /*** Create group structure. ***/
  ndg = (struct ndgroup_struct *)myalloc(sizeof(struct ndgroup_struct));
  ndg->n = 1;
  ndg->name = (char **)myalloc(ndg->n*sizeof(char *));
  ndg->name[0] = strdup("all");
  ndg->value = (char **)myalloc(ndg->n*sizeof(char *));
  ndg->value[0] = strdup("0");
  ndg->cnt = (int *)myalloc(ndg->n*sizeof(int));
  ndg->cnt[0] = count_non_zero_iarray(cflag,nd->ntrial);
  /*count_trials_ndata_condition(nd,ndc);*/

  ndg->tnum = (int **)myalloc(ndg->n*sizeof(int *));
  ndg->tnum[0] = (int *)myalloc(ndg->cnt[0]*sizeof(int));

  /*** Fill the trial list with all trials satisfying the condition "ndc" ***/
  k = 0;
  for(i=0;i<nd->ntrial;i++)
    if (cflag[i]){
      ndg->tnum[0][k] = i;
      k += 1;
    }
  *rndg = ndg;
  myfree(cflag);
}
/**************************************-**************************************/
/*                                                                           */
/*                          GET_GROUP_ALL_UNCONDITIONAL                      */
/*                                                                           */
/*  Make one group containing all trials.                                    */
/*                                                                           */
/*****************************************************************************/
void get_group_all_unconditional(nd,ndc,rndg)
     struct ndata_struct *nd;
     struct ndcond_struct *ndc;
     struct ndgroup_struct **rndg;
{
  int i,k;
  int *cflag;
  struct ndgroup_struct *ndg;

  printf("  GET_GROUP_ALL_UNCONDITIONAL\n");

  /*** Create group structure. ***/
  ndg = (struct ndgroup_struct *)myalloc(sizeof(struct ndgroup_struct));
  ndg->n = 1;
  ndg->name = (char **)myalloc(ndg->n*sizeof(char *));
  ndg->name[0] = strdup("all");
  ndg->value = (char **)myalloc(ndg->n*sizeof(char *));
  ndg->value[0] = strdup("0");
  ndg->cnt = (int *)myalloc(ndg->n*sizeof(int));
  ndg->cnt[0] = count_trials_ndata_condition(nd,ndc);
  ndg->tnum = (int **)myalloc(ndg->n*sizeof(int *));
  ndg->tnum[0] = (int *)myalloc(ndg->cnt[0]*sizeof(int));

  /*** Fill the trial list with all trials satisfying the condition "ndc" ***/
  cflag = ndata_condition_get_flags(nd,ndc);
  k = 0;
  for(i=0;i<nd->ntrial;i++)
    if (cflag[i]){
      ndg->tnum[0][k] = i;
      k += 1;
    }
  myfree(cflag);
  *rndg = ndg;
}
/**************************************-**************************************/
/*                                                                           */
/*                            GET_GROUP_INDIVIDUAL                           */
/*                                                                           */
/*  Make a group of size 1 for each trial satisfying the conditions.         */
/*                                                                           */
/*****************************************************************************/
void get_group_individual(nd,nda,rndg)
     struct ndata_struct *nd;
     struct nda_struct *nda;
     struct ndgroup_struct **rndg;
{
  int i,k;
  int *cflag;
  struct ndgroup_struct *ndg;
  char temp[SLEN];

  printf("  GET_GROUP_INDIVIDUAL\n");

  /*** Create group structure. ***/
  ndg = (struct ndgroup_struct *)myalloc(sizeof(struct ndgroup_struct));
  ndg->n = count_trials_ndata_condition(nd,nda->ndc);
  printf("    %d groups.\n",ndg->n);
  ndg->name = (char **)myalloc(ndg->n*sizeof(char *));
  ndg->value = (char **)myalloc(ndg->n*sizeof(char *));
  ndg->cnt = (int *)myalloc(ndg->n*sizeof(int));
  ndg->tnum = (int **)myalloc(ndg->n*sizeof(int *));

  cflag = ndata_condition_get_flags(nd,nda->ndc);
  k = 0;
  for(i=0;i<nd->ntrial;i++)
    if (cflag[i]){
      sprintf(temp,"trial_%d",i);
      ndg->name[k] = strdup(temp);
      sprintf(temp,"%d",i);
      ndg->value[k] = strdup(temp);
      ndg->cnt[k] = 1;
      ndg->tnum[k] = (int *)myalloc(sizeof(int));
      ndg->tnum[k][0] = i;
      k += 1;
    }
  myfree(cflag);
  *rndg = ndg;
}
/**************************************-**************************************/
/*                                                                           */
/*                          NDATA_GET_GROUPDEF_GROUP                         */
/*                                                                           */
/*  Get a group containing trials satisfying the "groupdef" string:          */
/*                                                                           */
/*    groupdef <pname1> <pval1> ... <pnamen> <pvaln>                         */
/*                                                                           */
/*****************************************************************************/
void ndata_get_groupdef_group(gstr,nd,rndg)
     char *gstr;
     struct ndata_struct *nd;
     struct ndgroup_struct **rndg;
{
  int i;
  int nc,ns;
  char **slist,tstr[SLEN],*pname,*pval,gname[SLEN],temp[SLEN];
  struct ndcond_struct *ndc;

  get_items_from_string(gstr,&slist,&ns);
  nc = (ns - 1)/2;
  if (nc == 0)
    exit_error("NDATA_GET_GROUPDEF_GROUP","No parameters in groupdef");
  else if (((ns-1) % 2) != 0)
    exit_error("NDATA_GET_GROUPDEF_GROUP","Name to value pairing violated");

  gname[0] = '\0';
  ndc = get_empty_condition(nc);
  for(i=0;i<nc;i++){
    pname = slist[2*i+1];
    pval = slist[2*i+2];
    sprintf(tstr,"param_range %s %s %s",pname,pval,pval);
    get_items_from_string(tstr,&(ndc->slist[i]),&(ndc->ns[i]));
    sprintf(temp,"%s=%s_",pname,pval);
    strcat(gname,temp);
  }
  get_group_all(nd,ndc,rndg);
  gname[strlen(gname)-1] = '\0'; // Remove the final '_'
  (*rndg)->name[0] = strdup(gname);

  ndcond_free_condition(ndc);
  free_2d_carray(slist,ns);
}
/**************************************-**************************************/
/*                                                                           */
/*                              GET_GROUPED_DATA                             */
/*                                                                           */
/*  Make a group for each unique set of values for the parameters specified  */
/*  which satisfy the conditions.                                            */
/*                                                                           */
/*****************************************************************************/
void get_grouped_data(nd,nda,rndg)
     struct ndata_struct *nd;
     struct nda_struct *nda;
     struct ndgroup_struct **rndg;
{
  int i,j,k,l;
  int ntrial,nparam,flag,gcount,tcount;
  struct ndgroup_struct *ndg;
  char *gname;

  printf("  GET_GROUPED_DATA\n");

  // Cases "all" and "individual" are handled by simpler routines
  if (strcmp(nda->grouping,"all")==0){
    get_group_all(nd,nda->ndc,rndg);
    return;
  }else if (strcmp(nda->grouping,"individual")==0){
    get_group_individual(nd,nda,rndg);
    return;
  }else if (strcmp(nda->grouping,"group")!=0){
    printf("  grouing type:  %s\n",nda->grouping);
    exit_error("GET_GROUPED_DATA","Unknown grouping type"); 
  }

  // Count groups
  nparam = nda->ndx->n;
  ntrial = nda->ndx->nrec;
  if (ntrial <= 0)
    exit_error("GET_GROUPED_DATA","No trials to anlyze");
  gcount = 1; // Count the first group
  for(i=1;i<ntrial;i++){
    k = nda->ndx->ndx[i];
    l = nda->ndx->ndx[i-1];
    flag = 0;
    for(j=0;j<nparam;j++)
      if (nda->ndx->ndxval[k][j] != nda->ndx->ndxval[l][j])
	flag = 1;
    gcount += flag; // Increment count each time the param values change
  }
  printf("    %d groups.\n",gcount);

  // Create group structure
  ndg = (struct ndgroup_struct *)myalloc(sizeof(struct ndgroup_struct));
  ndg->n = gcount;
  ndg->name = (char **)myalloc(ndg->n*sizeof(char *));
  ndg->value = NULL;
  ndg->cnt = (int *)myalloc(ndg->n*sizeof(int));
  ndg->tnum = (int **)myalloc(ndg->n*sizeof(int *));

  gname = ndata_make_name_index_trial(nd,nda->ndx,0);
  gcount = 1; // Count the first group
  tcount = 1; // Count trials in current group
  for(i=1;i<ntrial;i++){
    k = nda->ndx->ndx[i];
    l = nda->ndx->ndx[i-1];
    flag = 0; // flag will be set if param value changes
    for(j=0;j<nparam;j++)
      if (nda->ndx->ndxval[k][j] != nda->ndx->ndxval[l][j])
	flag = 1;
    if (flag){
      ndg->cnt[gcount-1] = tcount;
      ndg->tnum[gcount-1] = (int *)myalloc(tcount*sizeof(int));
      for(j=0;j<tcount;j++)
	ndg->tnum[gcount-1][j] = nda->ndx->tnum[nda->ndx->ndx[i-tcount+j]];
      tcount = 0;
      ndg->name[gcount-1] = strdup(gname);
      gcount += 1;
      gname = ndata_make_name_index_trial(nd,nda->ndx,i);
    }
    tcount += 1;
  }

  // Set values for last group
  ndg->cnt[gcount-1] = tcount;
  ndg->tnum[gcount-1] = (int *)myalloc(tcount*sizeof(int));
  for(j=0;j<tcount;j++)
    ndg->tnum[gcount-1][j] = nda->ndx->tnum[nda->ndx->ndx[i-tcount+j]];
  ndg->name[gcount-1] = strdup(gname);

  *rndg = ndg;
}
/**************************************-**************************************/
/*                                                                           */
/*                       GET_PRIMARY_LIST_FROM_GROUPED_DATA                  */
/*                                                                           */
/*  Determine the list of groups which consist of the primary points on a    */
/*  tuning curve.  This is designed to separate regular stimuli from         */
/*  control stimuli.                                                         */
/*                                                                           */
/*****************************************************************************/
void get_primary_list_from_grouped_data(nd,nda,ndg,pname,rpindex,rpn)
     struct ndata_struct *nd;
     struct nda_struct *nda;
     struct ndgroup_struct *ndg;
     char *pname;
     int **rpindex,*rpn;
{
  int i,j,k;
  int *nuval,*pindex,flag,result,t;
  float *frac,pfrac;
  char **comval,*vstring;

  printf("  GET_PRIMARY_LIST_FROM_GROUPED_DATA\n");

  ndata_get_var_param_info(nd,&nuval,&frac,&comval);
  pfrac = -1.0;
  for(i=0;i<nd->nvar;i++)
    if (strcmp(nd->vname[i],pname)==0)
      pfrac = frac[i];
  if (pfrac < 0.0)
    exit_error("GET_PRIMARY_LIST_FROM_GROUPED_DATA","pfrac < 0.0");

  /*** Check each group to see whether it is "primary" or a "control". ***/
  pindex = (int *)myalloc(ndg->n*sizeof(int)); /* Index to primary groups. */
  k = 0; /* Current number of primary groups. */
  for(i=0;i<ndg->n;i++){
    t = ndg->tnum[i][0]; /* Use first trial of each group. */
    flag = 1; /* Assume this will be a primary group. */
    for(j=0;j<nd->nvar;j++){ /*** Check each variable parameter. ***/
      if (frac[j] > pfrac){ /* This varies less than the primary parameter. */
	result = ndata_get_var_value_trial(nd,t,nd->vname[j],&vstring);
	if (result == 0)
	  exit_error("GET_PRIMARY_LIST_FROM_GROUPED_DATA","vname not found");
	if (strcmp(vstring,comval[j])!=0) /* Not the most common value. */
	  flag = 0;
	myfree(vstring);
      }
    }
    if (flag){
      pindex[k] = i;
      k += 1;
    }
  }
  myfree(frac); myfree(nuval); free_2d_carray(comval,nd->nvar);
  *rpindex = pindex; *rpn = k;
}
/**************************************-**************************************/
/*                                                                           */
/*                              PREPARE_NDA_ARG                              */
/*                                                                           */
/*  This is the beginning of all analyses (in 'nda.c' and 'nda_corr.c'.      */
/*  If the start time is redefined, that should happen here so that trial    */
/*  conditions involving times are aligned consistently.                     */
/*                                                                           */
/*****************************************************************************/
void prepare_nda_arg(afile,infile,outfile,slist,ns,rnd,rnda,rndg)
     char afile[],infile[],outfile[];
     char **slist;   // Command line arguments [ns], name val name val ...
     int ns;         // There are ns/2 name-value pairs
     struct ndata_struct **rnd;
     struct nda_struct **rnda;
     struct ndgroup_struct **rndg;
{
  int i,k;
  int tn;
  char name[SLEN],**tname,**tval,tstr[SLEN];

  get_name_without_path(infile,name);

  read_ndata_analysis(afile,name,rnda);
  if (*rnda == NULL)
    exit_error("PREPARE_NDA","Failed to read nda file");

  read_ndata(infile,rnd);
  if (*rnd == NULL)
    exit_error("PREPARE_NDA","Failed to read ndata file");

  // Add or update parameter definitions from command line
  tname = (*rnda)->pname;
  tval  = (*rnda)->pval;
  tn    = (*rnda)->nparam;
  for(i=0;i<ns;i+=2){
    if ((i+1) < ns){
      k = search_2d_carray(tname,slist[i],tn);
      if (k < 0){
	sprintf(tstr,"Name %s does not appear in .nda file.",slist[i]);
	exit_error("PREPARE_NDA_ARG",tstr);
      }
      myfree((*rnda)->pval[k]);
      (*rnda)->pval[k] = strdup(slist[i+1]);
    }else{
      printf("*** WARNING last parameter name ignored---no value.\n");
    }
  }

  ndata_adjust_start_time_global(*rnd,*rnda);

  // Make args float here
  ndata_convert_var_numeric(*rnd);

  set_ndata_analysis_index(*rnd,*rnda); // Conditions set in here
  get_grouped_data(*rnd,*rnda,rndg);
}
/**************************************-**************************************/
/*                                                                           */
/*                                PREPARE_NDA                                */
/*                                                                           */
/*  This is the beginning of all analyses (in 'nda.c' and 'nda_corr.c'.      */
/*  If the start time is redefined, that should happen here so that trial    */
/*  conditions involving times are aligned consistently.                     */
/*                                                                           */
/*****************************************************************************/
void prepare_nda(afile,infile,outfile,rnd,rnda,rndg)
     char afile[],infile[],outfile[];
     struct ndata_struct **rnd;
     struct nda_struct **rnda;
     struct ndgroup_struct **rndg;
{
  char name[SLEN];

  get_name_without_path(infile,name);

  read_ndata_analysis(afile,name,rnda);
  if (*rnda == NULL)
    exit_error("PREPARE_NDA","Failed to read nda file");
  read_ndata(infile,rnd);
  if (*rnd == NULL)
    exit_error("PREPARE_NDA","Failed to read ndata file");

  ndata_adjust_start_time_global(*rnd,*rnda);

  set_ndata_analysis_index(*rnd,*rnda);
  get_grouped_data(*rnd,*rnda,rndg);
}
/**************************************-**************************************/
/*                                                                           */
/*                      NDATA_GET_MAX_PERIOD_SAMPLING_CHAN                   */
/*                                                                           */
/*****************************************************************************/
void ndata_get_max_period_sampling_chan(nd,chan,rperiod,rsampling)
     struct ndata_struct *nd;
     char chan[];
     int *rperiod;
     float *rsampling;
{
  struct ndrx_struct *ndrx;

  /*** Determine start and period. ***/
  ndata_get_record_index(nd,chan,&ndrx);

  *rperiod = get_max_period_record(ndrx);
  /*** Get sampling rate for this channel. ***/

  *rsampling = ndata_get_sampling_record_index(nd,ndrx);
  if (*rsampling == -1.0)
    exit_error("NDATA_GET_MAX_PERIOD_SAMPLING_CHAN",
	       "Sampling varies across trials");

  ndata_free_ndrx(ndrx);
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_ESTABLISH_BASIC_PARAMS                      */
/*                                                                           */
/*  If parameter values are specified in "nda", they are used, otherwise,    */
/*    "chan" is the name of the first record in the first group              */
/*    "start" is 0                                                           */
/*    "period" is maximum period for "chan"                                  */
/*    "sampling" is the sampling value associated with "chan".               */
/*                                                                           */
/*  WYETH - WARNING - 'start' and 'period' may be in data file sampling      */
/*  units (if taken from nd file) or in other units (if taken from the       */
/*  nda file).                                                               */
/*                                                                           */
/*****************************************************************************/
void ndata_establish_basic_params(nd,nda,flag,rchan,rstart,rperiod,rsampling)
     struct ndata_struct *nd;
     struct nda_struct *nda;
     int flag;      // 0 - Find the first channel name starting with "unit..."
                    // 1 - Find the first channel name starting w/ '*rchan'
                    // 3 - Find the first event channel
     char **rchan;
     int *rstart,*rperiod;
     float *rsampling;
{
  struct ndrx_struct *ndrx;
  int max_period;
  float st0,stn,ssamp;
  char *defchan,*sparam;

  // Determine the channel name
  defchan = NULL;
  if (flag == 0){
    ndata_get_first_unit_chan(nd,&defchan);
  }else if (flag == 1){
    ndata_get_first_prefix_chan(nd,*rchan,&defchan);
  }else if (flag == 3){
    ndata_get_first_event_chan(nd,&defchan);
  }else
    exit_error("NDATA_ESTABLISH_BASIC_PARAMS","Unknown flag value");

  if (defchan == NULL){
    if (ndata_test_nda_param(nda,"chan") == 0)
      exit_error("NDATA_ESTABLISH_BASIC_PARAMS","Cannot find unit record");
  }

  ndata_set_nda_param_default_char(nda,"chan",rchan,defchan);
  myfree(defchan);

  // Determine start and period
  ndata_get_record_index(nd,*rchan,&ndrx);
  max_period = get_max_period_record(ndrx);

  if (ndata_test_nda_param(nda,"start_param")){
    sparam = ndata_get_nda_param_char_or_exit(nda,"start_param");
    st0 = ndata_get_const_param_float(nd,sparam);
    ndata_set_nda_param_default_float(nda,"start_param_scale",&ssamp,1000.0);
    st0 *= ssamp;
    *rstart = st0;
    myfree(sparam);
  }else{
    ndata_set_nda_param_default_int(nda,"start",rstart,0);
  }

  if (ndata_test_nda_param(nda,"period_param")){
    sparam = ndata_get_nda_param_char_or_exit(nda,"period_param");
    stn = ndata_get_const_param_float(nd,sparam);
    stn *= ssamp;
    *rperiod = stn;
    myfree(sparam);
  }else{
    ndata_set_nda_param_default_int(nda,"period",rperiod,max_period);
  }
  
  // Get sampling rate for this channel
  *rsampling = ndata_get_sampling_record_index(nd,ndrx);
  if (*rsampling == -1.0)
    exit_error("SPIKE_PLOT_NDA","Sampling varies across trials");


  ndata_free_ndrx(ndrx);
}
/**************************************-**************************************/
/*                                                                           */
/*                    NDATA_ESTABLISH_START_PERIOD_SAMPLING                  */
/*                                                                           */
/*  'start' and 'period' are returned in 'asampling' units.  'sampling'      */
/*  is the value for the data on channel 'chan'.                             */
/*                                                                           */
/*****************************************************************************/
void ndata_establish_start_period_sampling(nd,nda,chan,rstart,rperiod,
					   rasampling,rsampling)
     struct ndata_struct *nd;
     struct nda_struct *nda;
     char *chan;
     int *rstart,*rperiod;
     float *rasampling;       // specified by .nda analysis file
     float *rsampling;        // specified for 'chan' in .nd file
{
  struct ndrx_struct *ndrx;
  int max_period;

  ndata_get_record_index(nd,chan,&ndrx);

  // Get sampling rate for this channel
  *rsampling = ndata_get_sampling_record_index(nd,ndrx);
  if (*rsampling == -1.0)
    exit_error("NDATA_ESTABLISH_START_PERIOD_SAMPLING",
	       "Sampling varies across trials");

  // Get sampling rate for analysis (for 'start' and 'period')
  ndata_set_nda_param_default_float(nda,"sampling",rasampling,*rsampling);

  // Get start time from analysis file, or use time 0
  ndata_set_nda_param_default_int(nda,"start",rstart,0);

  // Determine period
  ndata_set_nda_param_default_int(nda,"period",rperiod,-1);
  if (*rperiod == -1){
    max_period = get_max_period_record(ndrx);
    *rperiod = (float)max_period/(*rsampling) * (*rasampling);
  }
  ndata_free_ndrx(ndrx);
}
/**************************************-**************************************/
/*                                                                           */
/*                              MAKE_NDATA_SARRAY                            */
/*                                                                           */
/*  Make an ndata file from the sarray.                                      */
/*                                                                           */
/*  Notes:                                                                   */
/*  - vval[i][j] holds the value of the jth var param on the ith trial.      */
/*                                                                           */
/*****************************************************************************/
struct ndata_struct *make_ndata_sarray(s,cnt,n,period,ncon,nvar,cname,ctype,
				       cval,vname,vtype,vval,class,sampling)
     int **s,*cnt,n,period,ncon,nvar;
     char **cname,*ctype,**cval,**vname,*vtype,***vval,*class;
     float sampling;
{
  int i;
  struct ndata_struct *nd;

  nd = (struct ndata_struct *)myalloc(sizeof(struct ndata_struct));
  nd->class = strdup(class);
  init_set_const_param(nd,ncon,cname,ctype,cval);
  nd->ntable = 0;
  nd->table = NULL;
  nd->ntrial = n;
  nd->t = (struct ndtrial_struct *)myalloc(n*sizeof(struct ndtrial_struct));
  init_set_var_param(nd,n,nvar,vname,vtype,vval);

  for(i=0;i<n;i++){
    nd->t[i].tcode = i;
    nd->t[i].tref = i;
    nd->t[i].nrec = 1;
    nd->t[i].r = (struct ndrec_struct *)myalloc(sizeof(struct ndrec_struct));
    nd->t[i].r[0].rtype = 0;
    nd->t[i].r[0].name = strdup("unit0");
    nd->t[i].r[0].rcode = 0;
    nd->t[i].r[0].sampling = sampling;
    nd->t[i].r[0].t0 = 0;
    nd->t[i].r[0].tn = period;
    nd->t[i].r[0].n = cnt[i];
    nd->t[i].r[0].p = copy_iarray(s[i],cnt[i]);
  }
  return nd;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MAKE_NDATA_SARRAY_PAIR                          */
/*                                                                           */
/*  Make an ndata file with 2 spike channels.                                */
/*                                                                           */
/*  Notes:                                                                   */
/*  - vval[i][j] holds the value of the jth var param on the ith trial.      */
/*                                                                           */
/*****************************************************************************/
struct ndata_struct *make_ndata_sarray_pair(s1,cnt1,s2,cnt2,n,period,ncon,nvar,
					    cname,ctype,cval,vname,vtype,vval,
					    class,sampling)
     int **s1,*cnt1,**s2,*cnt2,n,period,ncon,nvar;
     char **cname,*ctype,**cval,**vname,*vtype,***vval,*class;
     float sampling;
{
  int i;
  struct ndata_struct *nd;

  nd = (struct ndata_struct *)myalloc(sizeof(struct ndata_struct));
  nd->class = strdup(class);

  init_set_const_param(nd,ncon,cname,ctype,cval);
  nd->ntable = 0;
  nd->table = NULL;
  nd->ntrial = n;
  nd->t = (struct ndtrial_struct *)myalloc(n*sizeof(struct ndtrial_struct));
  init_set_var_param(nd,n,nvar,vname,vtype,vval);

  for(i=0;i<n;i++){
    nd->t[i].tcode = i;
    nd->t[i].tref = i;
    nd->t[i].nrec = 2;
    nd->t[i].r = (struct ndrec_struct *)myalloc(2*sizeof(struct ndrec_struct));

    nd->t[i].r[0].rtype = 0;
    nd->t[i].r[0].name = strdup("unit0");
    nd->t[i].r[0].rcode = 0;
    nd->t[i].r[0].sampling = sampling;
    nd->t[i].r[0].t0 = 0;
    nd->t[i].r[0].tn = period;
    nd->t[i].r[0].n = cnt1[i];
    nd->t[i].r[0].p = copy_iarray(s1[i],cnt1[i]);

    nd->t[i].r[1].rtype = 0;
    nd->t[i].r[1].name = strdup("unit1");
    nd->t[i].r[1].rcode = 0;
    nd->t[i].r[1].sampling = sampling;
    nd->t[i].r[1].t0 = 0;
    nd->t[i].r[1].tn = period;
    nd->t[i].r[1].n = cnt2[i];
    nd->t[i].r[1].p = copy_iarray(s2[i],cnt2[i]);
  }
  return nd;
}
/**************************************-**************************************/
/*                                                                           */
/*                       NDATA_GET_NORMD_RESP_CHAN_VDUR                      */
/*                                                                           */
/*  Compute the group-wise z-score normalized response for all trials for    */
/*  the specified channel.                                                   */
/*                                                                           */
/*  *** VARIABLE TRIAL DURATION ***                                          */
/*                                                                           */
/*****************************************************************************/
float *ndata_get_normd_resp_chan_vdur(nd,nda,ndg,chan,start,period,
				      sampling,rrn)
     struct ndata_struct *nd;
     struct nda_struct *nda;
     struct ndgroup_struct *ndg;
     char chan[];
     int start,period;
     float sampling;
     int *rrn;
{
  int i,j;
  int **s,*cnt,n,rn,tnum,*flag,*vdurg,*psthn,fn,vdur,spkcnt,flagcnt;
  char **names;
  float *r,mean,fbin,*fpsth;

  printf("  NDATA_GET_NORMD_RESP_CHAN_VDUR\n");
  
  rn = nd->ntrial;
  r = (float *)myalloc(rn*sizeof(float));
  flag = get_zero_iarray(rn);

  for(i=0;i<ndg->n;i++){ // Compute firing rates for each group
    get_sarray_ndata_group(nd,ndg,i,chan,&s,&cnt,&n,&names);

    // 'vdurg' holds variable duration for each trial in group
    vdurg = ndu_vardur_group(nd,ndg,i,nda);  // Variable duration

    // Get PSTH from variable duration trials
    fbin = 1.0;  // binsize
    su_psth_float_bin_sarray_vdur(s,cnt,n,(float)start,vdurg,fbin,
				  &fpsth,&fn,&psthn);
    // 'psthn' [fn] contains the number of trials that covered each bin

    // Integrate the PSTH to get the mean spike count vs. duration
    for(j=0;j<fn;j++){
      fpsth[j] /= (float)psthn[j];  // Conver to spikes/bin
      if (j > 0)
	fpsth[j] += fpsth[j-1];  // Integrate
    }

    for(j=0;j<ndg->cnt[i];j++){
      tnum = ndg->tnum[i][j];

      vdur = vdurg[j];
      mean = fpsth[vdur-1];  // Integrated spike count for duration
      spkcnt = count_spikes(s[j],cnt[j],start,vdur);

      if (mean == 0.0){
	if (spkcnt == 0.0)
	  r[tnum] = 0.0;
	else{
	  r[tnum] = (float)spkcnt - mean;  // - mean / SD
	  printf("  *** WARNING:  mean is zero - assuming SD = 1\n");
	  printf("    %3d vdur: %5d   mean: %f   spkcnt: %d\n",j,vdur,mean,
		 spkcnt);
	}
      }else{
	//  NOTE, If Poisson, mean = var
	r[tnum] = ((float)spkcnt - mean) / sqrt(mean);  // - mean / SD
      }
      //printf("%3d vdur: %5d   mean: %f   spkcnt: %d\n",j,vdur,mean,spkcnt);
  
      flag[tnum] = 1;
    }
    myfree(s); myfree(cnt); free_2d_carray(names,n);
  }
  flagcnt = sum_iarray(flag,rn,0,rn);
  if (flagcnt != rn){
    printf(" *** WARNING:  %d trials failed condition, response set to 0\n",
	   rn - flagcnt);
    //exit_error("NDATA_GET_NORMALIZED_RESPONSE_CHAN","Some flags are 0");
  }
  myfree(flag);

  *rrn = rn;
  return r;
}
/**************************************-**************************************/
/*                                                                           */
/*                      NDATA_GET_NORMALIZED_RESPONSE_CHAN                   */
/*                                                                           */
/*  Compute the group-wise z-score normalized response for all trials for    */
/*  the specified channel.                                                   */
/*                                                                           */
/*****************************************************************************/
float *ndata_get_normalized_response_chan(nd,nda,ndg,chan,start,period,
					  sampling,rrn)
     struct ndata_struct *nd;
     struct nda_struct *nda;
     struct ndgroup_struct *ndg;
     char chan[];
     int start,period;
     float sampling;
     int *rrn;
{
  int i,j;
  int **s,*cnt,n,rn,tnum,*flag;
  char **names;
  float *r,mean,sdev;
  /*float *sc;*/

  rn = nd->ntrial;
  r = (float *)myalloc(rn*sizeof(float));
  /*sc = (float *)myalloc(rn*sizeof(float));*/
  flag = get_zero_iarray(rn);

  for(i=0;i<ndg->n;i++){ // Compute firing rates for each group
    get_sarray_ndata_group(nd,ndg,i,chan,&s,&cnt,&n,&names);

    mean = mean_spike_count_sarray(s,cnt,n,start,period,sampling,&sdev);
    for(j=0;j<ndg->cnt[i];j++){
      tnum = ndg->tnum[i][j];
      /*sc[tnum] = (float)count_spikes(s[j],cnt[j],start,period);*/
      if (sdev > 0.0)
	r[tnum] = ((float)count_spikes(s[j],cnt[j],start,period) - mean)/sdev;
      else
	r[tnum] = 0.0;
      flag[tnum] = 1;

      /*printf("name = %s\n",ndg->name[i]);*/
      /*** Used for printing points for TCC figure in Udi paper
	if (strcmp(ndg->name[i],"coherence=999")==0)
	printf("%d %.4f\n",tnum,r[tnum]);***/
    }
    myfree(s); myfree(cnt); free_2d_carray(names,n);
  }
  if (sum_iarray(flag,rn,0,rn) != rn){
    printf(" *** WARNING: some flags are zero - WYETH - is this needed?\n");
    //exit_error("NDATA_GET_NORMALIZED_RESPONSE_CHAN","Some flags are 0");
  }
  myfree(flag);

  /***append_farray_plot("zz.zz","raw_sc",sc,rn,1);
    myfree(sc);***/

  *rrn = rn;
  return r;
}
/**************************************-**************************************/
/*                                                                           */
/*                        NDATA_GET_ADJUSTED_SPIKE_COUNT                     */
/*                                                                           */
/*  Adjust each spike count by taking into account the value of the          */
/*  normalized response function.                                            */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - "cflag" 1 to include, 0 to omit, center value when smoothing.          */
/*  - "adjflag" determines the adjustment algorithm                          */
/*                                                                           */
/*****************************************************************************/
float *ndata_get_adjusted_spike_count(nd,ndg,chan,start,period,sampling,sigma,
				      r,rn,cflag,adjflag)
     struct ndata_struct *nd;
     struct ndgroup_struct *ndg;
     char chan[];
     int start,period;
     float sampling; /**** ,*rmean,*rsdev; ****/
     float sigma,*r; /* Normalized responses for all trials. */
     int rn,cflag,adjflag;
{
  int i,j;
  int **s,*cnt,n,tnum,*flag;
  char **names;
  float *radj,*rsm,sc,mean,sdev;

  radj = (float *)myalloc(rn*sizeof(float));
  flag = get_zero_iarray(rn); /* To check if groups account for all trials. */

  if (cflag == 1)
    rsm = smooth_with_gaussian(r,rn,sigma,0.01);
  else{
    rsm = smooth_with_gaussian_omit_center(r,rn,sigma,0.01);
    /* write_farray_xy_plot("zz.debug.pl","rsm_v_r",r,rsm,rn);*/
    /*append_farray_xy_plot("zz.debug.pl",r,rsm,rn,"rsm_v_r");*/
    /*rsm = gaussian_corr_noise(0.0,0.0,0.2,rn,1777); TO GENERATE NOISE */
    /*printf("NOISE for rsm\n");*/
  }

  /***append_farray_plot("zz.zz","r",r,rn,1);
    append_farray_plot("zz.zz","r_smooth",rsm,rn,1);***/

  for(i=0;i<ndg->n;i++){ /*** Compute firing rates for each group. ***/
    get_sarray_ndata_group(nd,ndg,i,chan,&s,&cnt,&n,&names);
    mean = mean_spike_count_sarray(s,cnt,n,start,period,sampling,&sdev);
    for(j=0;j<ndg->cnt[i];j++){
      tnum = ndg->tnum[i][j];
      if (sdev > 0.0){
	sc = (float)count_spikes(s[j],cnt[j],start,period);
	if (adjflag == 1){
	  radj[tnum] = (sc - rsm[tnum]*sdev);
	  /*printf("SDEV\n");*/
	}else if (adjflag == 2)
	  radj[tnum] = (sc - rsm[tnum]*sqrt(mean));
	else
	  exit_error("NDATA_GET_ADJUSTED_SPIKE_COUNT","Unknown adjflag value");
      }else
	radj[tnum] = 0.0;
      flag[tnum] = 1;
    }
    myfree(s); myfree(cnt); free_2d_carray(names,n);
  }
  if (sum_iarray(flag,rn,0,rn) != rn)
    exit_error("NDATA_GET_ADJUSTED_SPIKE_COUNT","Some flags are 0");

  myfree(flag); myfree(rsm);

  /***append_farray_plot("zz.zz","r_adj",radj,rn,1);
    exit(0);***/

  return radj;
}
/**************************************-**************************************/
/*                                                                           */
/*                       NDATA_GET_ADJUSTED_MEAN_SDEV_GROUP                  */
/*                                                                           */
/*  For the "k"th group, compute the mean and SD of the specified            */
/*  (probably adjusted) response measure "r".                                */
/*                                                                           */
/*****************************************************************************/
void ndata_get_adjusted_mean_sdev_group(ndg,k,r,rn,rmean,rsdev)
     struct ndgroup_struct *ndg; 
     int k;                  // index of group in 'ndg'
     float *r;               // [rn] response metric
     int rn;                 // Number of trials overall
     float *rmean,*rsdev;    // Return mean and SD
{
  int i;
  int n,tnum;
  float *data;

  //
  //  Fill 'data' with the 'r' values for the trials in the 'k' group
  //
  n = ndg->cnt[k];
  data = (float *)myalloc(n*sizeof(float));
  for(i=0;i<n;i++){
    tnum = ndg->tnum[k][i];
    data[i] = r[tnum];
  }

  //  Analyze the selected 'r' values
  mean_sdev_farray(data,n,rmean,rsdev);

  myfree(data);
}
/**************************************-**************************************/
/*                                                                           */
/*                               READ_NDATA_GEN                              */
/*                                                                           */
/*  Read into the "ndgen_struct" which is used to hold a description for     */
/*  generating spike trains.  The file format is specified in "ndata.h".     */
/*                                                                           */
/*  Notes:                                                                   */
/*  - Parameter values can be a list of values, but are stored as 1 string.  */
/*                                                                           */
/*****************************************************************************/
void read_ndata_gen(infile,rg)
     char infile[];
     struct ndgen_struct **rg;
{
  FILE *fopen(),*fin;
  int i,k;
  int ns,done,nn;
  struct ndgen_struct *g;
  char temp[SLEN],cstr[SLEN],**slist,*cptr,*tc;

  printf("  READ_NDATA_GEN\n");

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** Cannot open file %s.\n",infile);
    *rg = NULL;
    return;
  }else
    printf("    Reading %s\n",infile);

  g = (struct ndgen_struct *)myalloc(sizeof(struct ndgen_struct));
  nn = fscanf(fin,"%s",temp);
  while(temp[0] == '#'){
    tc = fgets(temp,SLEN,fin);
    nn = fscanf(fin,"%s",temp);
  }
  g->gtype = strdup(temp);

  g->nparam = 0;
  done = 0;
  while(!done){
    if (fgets(temp,SLEN,fin)==NULL){
      done = 1;
    }else{
      get_items_from_string(temp,&slist,&ns);
      if (ns > 0){
	if ((strcmp(slist[0],"nvar")==0)||(strcmp(slist[0],"npop")==0)||
	    (strcmp(slist[0],"nvar_link")==0))
	  done = 1;
	else if (strcmp(slist[0],"nfiles")!=0) /* Don't count "nfile" */
	  g->nparam += 1;
	free_2d_carray(slist,ns);
      }
    }
  }
  fclose(fin);
  /*printf("    %d params found\n",g->nparam);*/
  g->pname = (char **)myalloc(g->nparam*sizeof(char *));
  g->ptype = (char *)myalloc(g->nparam*sizeof(char));
  g->pval = (char **)myalloc(g->nparam*sizeof(char *));

  fin = fopen(infile,"r");
  nn = fscanf(fin,"%s",temp);
  while(temp[0] == '#'){
    tc = fgets(temp,SLEN,fin);
    nn = fscanf(fin,"%s",temp);
  }

  done = 0;
  k = 0;
  while(!done){
    if (fscanf(fin,"%s",temp)==EOF)
      done = 1;
    else{
      if (strcmp(temp,"nfiles")==0)
	nn = fscanf(fin,"%d",&g->nfiles);
      else if ((strcmp(temp,"nvar")==0)||(strcmp(temp,"npop")==0)||
	       (strcmp(temp,"nvar_link")==0))
	done = 1;
      else{
	g->pname[k] = strdup(temp);
	if (fgets(temp,SLEN,fin)==NULL)
	  exit_error("READ_NDATA_GEN","Cannot find parameter values");
	get_items_from_string(temp,&slist,&ns);
	if (ns==1){
	  g->pval[k] = strdup(slist[0]);
	  g->ptype[k] = 'c'; /* Default, since no type given */
	}else if (ns==2){ /* Assume ns > 1 */
	  g->pval[k] = strdup(slist[1]);
	  if (!((strcmp(slist[0],"c")==0)||
		(strcmp(slist[0],"i")==0)||
		(strcmp(slist[0],"f")==0))){
	    printf("*** For pname = %s\n",g->pname[k]);
	    exit_error("READ_NDATA_GEN","Type specifier must be c, i, f");
	  }
	  g->ptype[k] = slist[0][0];
	}else{ /* assume ns > 2 */
	  if (!((strcmp(slist[0],"c")==0)||
		(strcmp(slist[0],"i")==0)||
		(strcmp(slist[0],"f")==0))){
	    printf("*** For pname = %s\n",g->pname[k]);
	    exit_error("READ_NDATA_GEN","Type specifier must be c, i, f");
	  }
	  g->ptype[k] = slist[0][0];
	  if (temp[strlen(temp)-1]=='\n')
	    temp[strlen(temp)-1]='\0';
	  cptr = (char *)temp + 2;
	  while(cptr[0] == ' ') /* Pass leading spaces before type char */
	    cptr += 1;
	  while(cptr[0] != ' ') /* Pass type char */
	    cptr += 1;
	  while(cptr[0] == ' ') /* Pass spaces after type char */
	    cptr += 1;
	  g->pval[k] = strdup(cptr);
	}
	free_2d_carray(slist,ns);
	k += 1;
      }
    }
  }

  g->nvar = 0;
  if ((strcmp(temp,"nvar")==0)||(strcmp(temp,"nvar_link")==0)){
    nn = fscanf(fin,"%d",&g->nvar);
    if (strcmp(temp,"nvar_link")==0)
      g->nvar_link = 1;
    else
      g->nvar_link = 0;
    g->vname = (char **)myalloc(g->nvar*sizeof(char *));
    g->vtype = (char *)myalloc(g->nvar*sizeof(char));
    g->vvalstr = (char **)myalloc(g->nvar*sizeof(char *));
    for(i=0;i<g->nvar;i++){
      nn = fscanf(fin,"%s %s",temp,cstr);
      if (search_2d_carray(g->pname,temp,g->nparam) != -1){
	printf("    name: %s\n",temp);
	exit_error("READ_NDATA_GEN","Variable parameter name conflict");
      }
      g->vname[i] = strdup(temp);
      g->vtype[i] = cstr[0];
      if (fgets(temp,SLEN,fin)==NULL)
	exit_error("READ_NDATA_GEN","Cannot find var parameter values");
      get_items_from_string(temp,&slist,&ns);
      g->vvalstr[i] = make_string_from_items(slist,ns);
      free_2d_carray(slist,ns);
      /*printf("g->vname[i] = %s  %s\n",g->vname[i],g->vvalstr[i]);*/
    }
    nn = fscanf(fin,"%s",temp); // This should get "npop" or nothing
  }

  g->npop = 0;
  if (strcmp(temp,"npop")==0){
    nn = fscanf(fin,"%d",&g->npop);
    g->popname = (char **)myalloc(g->npop*sizeof(char *));
    g->popdd = (char **)myalloc(g->npop*sizeof(char *));
    g->popval = (char **)myalloc(g->npop*sizeof(char *)); /* Filled in later */
    for(i=0;i<g->npop;i++){
      nn = fscanf(fin,"%s",temp);
      if (search_2d_carray(g->pname,temp,g->nparam) != -1){
	printf("    name: %s\n",temp);
	exit_error("READ_NDATA_GEN","Population parameter name conflict");
      }
      if (search_2d_carray(g->vname,temp,g->nvar) != -1){
	printf("    name: %s\n",temp);
	exit_error("READ_NDATA_GEN","Population parameter name conflict");
      }
      g->popname[i] = strdup(temp);
      if (fgets(temp,SLEN,fin)==NULL)
	exit_error("READ_NDATA_GEN","Cannot find pop. distr. description");
      get_items_from_string(temp,&slist,&ns);
      g->popdd[i] = make_string_from_items(slist,ns);
      free_2d_carray(slist,ns);
      /*printf("g->popname[i] = %s  %s\n",g->popname[i],g->popdd[i]);*/
    }
  }
  fclose(fin);

  *rg = g;
}
/**************************************-**************************************/
/*                                                                           */
/*                           NDGEN_GET_PARAM_CONST_POP                       */
/*                                                                           */
/*  Search the lists of constant and population parameters for value.        */
/*  Returns 1 if param "pname" found in list of parameters, 0 other-         */
/*  wise.  Return the character value of this parameter.                     */
/*                                                                           */
/*****************************************************************************/
int ndgen_get_param_const_pop(g,pname,rpval)
     struct ndgen_struct *g;
     char pname[];
     char **rpval;
{
  int i;
  char *t;

  t = NULL;
  i = search_2d_carray(g->pname,pname,g->nparam);
  if (i >= 0)
    t = strdup(g->pval[i]);
  else{
    i = search_2d_carray(g->popname,pname,g->npop);
    if (i >= 0)
      t = strdup(g->popval[i]);
  }
  if (i == -1){
    *rpval = NULL;
    return 0;
  }else{
    *rpval = t;
    return 1;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                      NDGEN_GET_PARAM_CONST_POP_DEF_INT                    */
/*                                                                           */
/*  Search the lists of constant and population parameters for value.        */
/*  If not found, return the default value, otherwise return the found       */
/*  value as an int.                                                         */
/*                                                                           */
/*****************************************************************************/
void ndgen_get_param_const_pop_def_int(g,pname,rparam,defval)
     struct ndgen_struct *g;
     char pname[];
     int *rparam,defval;
{
  int flag;
  char *tstr;

  flag = ndgen_get_param_const_pop(g,pname,&tstr);
  if (flag){
    *rparam = atoi(tstr);
    myfree(tstr);
  }else{
    *rparam = defval;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                     NDGEN_GET_PARAM_CONST_POP_DEF_FLOAT                   */
/*                                                                           */
/*  Search the lists of constant and population parameters for value.        */
/*  If not found, return the default value, otherwise return the found       */
/*  value as a float.                                                        */
/*                                                                           */
/*****************************************************************************/
void ndgen_get_param_const_pop_def_float(g,pname,rparam,defval)
     struct ndgen_struct *g;
     char pname[];
     float *rparam,defval;
{
  int flag;
  char *tstr;

  flag = ndgen_get_param_const_pop(g,pname,&tstr);
  if (flag){
    *rparam = atof(tstr);
    myfree(tstr);
  }else{
    *rparam = defval;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                     NDGEN_GET_PARAM_CONST_POP_DEF_CHAR                    */
/*                                                                           */
/*  Search the lists of constant and population parameters for value.        */
/*  If not found, return the default value, otherwise return the found       */
/*  value as a float.                                                        */
/*                                                                           */
/*****************************************************************************/
void ndgen_get_param_const_pop_def_char(g,pname,rparam,defval)
     struct ndgen_struct *g;
     char pname[],**rparam,*defval;
{
  int flag;
  char *tstr;

  flag = ndgen_get_param_const_pop(g,pname,&tstr);
  if (flag)
    *rparam = tstr;
  else
    *rparam = strdup(defval);
}
/**************************************-**************************************/
/*                                                                           */
/*                            NDGEN_GET_PARAM_TRIAL                          */
/*                                                                           */
/*  Returns 1 if param "pname" found in list of parameters, 0 other-         */
/*  wise.  Return the character value of this parameter.                     */
/*                                                                           */
/*  "vval" contains the variable parameter values for the current trial.     */
/*                                                                           */
/*****************************************************************************/
int ndgen_get_param_trial(g,pname,vval,rpval)
     struct ndgen_struct *g;
     char pname[],**vval,**rpval;
{
  int i;
  char *t;
  
  i = search_2d_carray(g->pname,pname,g->nparam);   /* Check const. params. */
  if (i >= 0)
    t = strdup(g->pval[i]);
  else{
    i = search_2d_carray(g->popname,pname,g->npop); /* Check pop. params. */
    if (i >= 0)
      t = strdup(g->popval[i]);
    else if (vval != NULL){
      i = search_2d_carray(g->vname,pname,g->nvar); /* Check var. params. */
      if (i >= 0)
	t = strdup(vval[i]);
      else{
	i = -1;
	t = NULL;
      }
    }else{
      i = -1;
      t = NULL;
    }
  }

  *rpval = t;

  if (i == -1)
    return 0;
  else
    return 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                          NDGEN_GET_PARAM_TRIAL_INT                        */
/*                                                                           */
/*  Return integer value of parameter if found.                              */
/*  "vval" contains the variable parameter values for the current trial.     */
/*                                                                           */
/*****************************************************************************/
int ndgen_get_param_trial_int(g,pname,vval)
     struct ndgen_struct *g;
     char pname[],**vval;
{
  int i,k;
  char *t;

  i = ndgen_get_param_trial(g,pname,vval,&t);
  if (i == 0){
    printf("  pname:  %s\n",pname);
    exit_error("NDGEN_GET_PARAM_TRIAL_INT","Cannot find parameter");
  }
  k = atoi(t);
  myfree(t);
  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                        NDGEN_GET_PARAM_TRIAL_FLOAT                        */
/*                                                                           */
/*  Return integer value of parameter if found.                              */
/*  "vval" contains the variable parameter values for the current trial.     */
/*                                                                           */
/*****************************************************************************/
float ndgen_get_param_trial_float(g,pname,vval)
     struct ndgen_struct *g;
     char pname[],**vval;
{
  int i;
  char *t;
  float x;
  
  i = ndgen_get_param_trial(g,pname,vval,&t);
  if (i == 0){
    printf("  pname:  %s\n",pname);
    exit_error("NDGEN_GET_PARAM_TRIAL_FLOAT","Cannot find parameter");
  }
  x = atof(t);
  myfree(t);
  return x;
}
/**************************************-**************************************/
/*                                                                           */
/*                        NDGEN_GET_PARAM_TRIAL_CHAR                         */
/*                                                                           */
/*  Return integer value of parameter if found.                              */
/*  "vval" contains the variable parameter values for the current trial.     */
/*                                                                           */
/*****************************************************************************/
char *ndgen_get_param_trial_char(g,pname,vval)
     struct ndgen_struct *g;
     char pname[],**vval;
{
  int i;
  char *t;
  
  i = ndgen_get_param_trial(g,pname,vval,&t);
  if (i == 0){
    printf("  pname:  %s\n",pname);
    exit_error("NDGEN_GET_PARAM_TRIAL_CHAR","Cannot find parameter");
  }
  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                            NDATA_GLIST_GET_INIT                           */
/*                                                                           */
/*****************************************************************************/
struct ndglist *ndata_glist_get_init()
{
  struct ndglist *g;

  g = (struct ndglist *)myalloc(sizeof(struct ndglist));
  g->n = 0;
  g->head = g->tail = NULL;

  return g;
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_GLIST_FREE_POINTERS                         */
/*                                                                           */
/*****************************************************************************/
void ndata_glist_free_pointers(g)
     struct ndglist *g;
{
  int i;
  struct ndgnode *curr,*next;

  curr = g->head;
  for(i=0;i<g->n;i++){
    next = curr->next;
    myfree(curr);
    curr = next;
  }
  myfree(g);
}
/**************************************-**************************************/
/*                                                                           */
/*                              NDATA_GLIST_PRINT                            */
/*                                                                           */
/*****************************************************************************/
void ndata_glist_print(g)
     struct ndglist *g;
{
  int tcnt,gcnt;
  struct ndgnode *p;
  struct ndtrial_struct *t;

  gcnt = 0;
  p = g->head;
  while(p!=NULL){
    printf("  Group %d has %d trials\n",gcnt,p->n);
    tcnt = 0;
    t = p->head;
    while(t!=NULL){
      printf("    trial %d tref = %d\n",tcnt,t->tref);
      tcnt += 1;
      t = t->next;
    }
    if (tcnt != p->n)
      printf("***WARNING found %d trials, not equal to stored value %d\n",
	     tcnt,p->n);
    gcnt += 1;
    p = p->next;
  }
  if (gcnt != g->n)
    printf("***WARNING found %d groups, not equal to stored value %d\n",
	   gcnt,g->n);
}
/**************************************-**************************************/
/*                                                                           */
/*                       NDATA_GLIST_FIND_GROUP_FOR_TRIAL                    */
/*                                                                           */
/*  Return a pointer to the appropriate group for this trial, or NULL if no  */
/*  group matches this trial.  Match is based on the trial "t" having the    */
/*  same variable parameters as the first trial in a group node's list.      */
/*                                                                           */
/*****************************************************************************/
struct ndgnode *ndata_glist_find_group_for_trial(g,t)
     struct ndglist *g;
     struct ndtrial_struct *t;
{
  struct ndgnode *p;

  p = g->head;
  while(p != NULL){
    if (p->n == 0)
      exit_error("NDATA_GLIST_FIND_GROUP_FOR_TRIAL","Empty gnode");
    if (ndata_match_trials_var_params(p->head,t) == 1) /* found a match */
      return p;
    else
      p = p->next; /* Point to the next gnode */
  }
  return p;
}
/**************************************-**************************************/
/*                                                                           */
/*                             NDATA_GLIST_ADD_GNODE                         */
/*                                                                           */
/*  Add a new group node to the list that has one trial "t" in its list.     */
/*  Make this the last group in the list.                                    */
/*                                                                           */
/*  Return a pointer to the newly added gnode.                               */
/*                                                                           */
/*****************************************************************************/
struct ndgnode *ndata_glist_add_gnode(g,t)
     struct ndglist *g;
     struct ndtrial_struct *t;
{
  struct ndgnode *p;

  p = (struct ndgnode *)myalloc(sizeof(struct ndgnode));
  p->n = 1;
  p->head = p->tail = t;
  t->next = t->prev = NULL;
  p->next = NULL; /* New last element has nothing after it. */
  p->prev = g->tail; /* New last element has old tail before it. */
  if (g->tail != NULL) /* If list is not empty */
    g->tail->next = p; /* Old last element now points to new last element. */
  g->tail = p; /* Tail points to new last element. */
  if (g->n == 0) /* This is first node in list */
    g->head = p;
  g->n += 1; /* Update number of groups */

  return p;
}
/**************************************-**************************************/
/*                                                                           */
/*                            NDATA_GNODE_ADD_TRIAL                          */
/*                                                                           */
/*  Add a new trial to the list pointed to by the gnode "gn".  Add to the    */
/*  tail of the list.                                                        */
/*                                                                           */
/*****************************************************************************/
void ndata_gnode_add_trial(gn,t)
     struct ndgnode *gn;
     struct ndtrial_struct *t;
{
  t->next = NULL; /* New last trial has nothing after it. */
  t->prev = gn->tail; /* New last trial has old tail before it. */
  if (gn->tail != NULL) /* If list not empty */
    gn->tail->next = t; /* Old last trial now points to new last trial. */
  gn->tail = t; /* Tail points to new last trial. */
  if (gn->n == 0) /* This is the first trial */
    gn->head = t;
  gn->n += 1; /* Update number of trials in group. */
}
/**************************************-**************************************/
/*                                                                           */
/*                            NDATA_GLIST_ADD_TRIAL                          */
/*                                                                           */
/*  Add the trial to an existing list for an appropriate gnode, or add a     */
/*  new gnode and start a new list if no appropriate group already exists.   */
/*                                                                           */
/*  Return a pointer to the group to which the trial was just added.         */
/*                                                                           */
/*****************************************************************************/
struct ndgnode *ndata_glist_add_trial(g,t)
     struct ndglist *g;
     struct ndtrial_struct *t;
{
  struct ndgnode *p;

  p = ndata_glist_find_group_for_trial(g,t);
  if (p==NULL)
    p = ndata_glist_add_gnode(g,t);
  else
    ndata_gnode_add_trial(p,t);

  return p;
}
/**************************************-**************************************/
/*                                                                           */
/*                     NDATA_GNODE_COUNT_TRIALS_WITH_CHAN                    */
/*                                                                           */
/*  Count the number of times the record "chan" appears in the gnode list.   */
/*                                                                           */
/*****************************************************************************/
int ndata_gnode_count_trials_with_chan(gn,chan)
     struct ndgnode *gn;
     char *chan;
{
  int count;
  struct ndtrial_struct *t;

  count = 0;
  t = gn->head;
  while(t!=NULL){
    if (ndata_trial_get_record_index(t,chan) >= 0)
      count += 1;
    t = t->next;
  }
  return count;
}
/**************************************-**************************************/
/*                                                                           */
/*                           NDATA_GNODE_GET_SARRAY                          */
/*                                                                           */
/*  Return lists of pointers to records for "chan" in the trials pointed     */
/*  to by the gnode.                                                         */
/*                                                                           */
/*****************************************************************************/
void ndata_gnode_get_sarray(gn,chan,rs,rcnt,rn)
     struct ndgnode *gn;
     char *chan;
     int ***rs,**rcnt,*rn;
{
  int i,k;
  int **s,*cnt,n;
  struct ndtrial_struct *t;

  n = ndata_gnode_count_trials_with_chan(gn,chan); /* number of valid trials */
  s = (int **)myalloc(n*sizeof(int *));
  cnt = (int *)myalloc(n*sizeof(int));
  k = 0;
  t = gn->head;
  while(t!=NULL){
    i = ndata_trial_get_record_index(t,chan);
    if (i >= 0){
      s[k] = t->r[i].p;
      cnt[k] = t->r[i].n;
      k += 1;
    }
    t = t->next;
  }
  *rs = s; *rcnt = cnt; *rn = n;
}
/**************************************-**************************************/
/*                                                                           */
/*                         NDATA_GET_MEAN_SDEV_2D_FDATA                      */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Assume fdata starts at time 0.                                         */
/*                                                                           */
/*****************************************************************************/
void ndata_get_mean_sdev_2d_fdata(fdata,m,n,start,period,dsamp,tsamp,rmean,rsd)
     float **fdata;
     int m,n,start,period;
     float dsamp,tsamp;    // dsamp applies to data, tsamp applies to times
     float *rmean,*rsd;
{
  int i;
  int t0,tn;
  float *t;

  //printf("m,n = %d %d  start,period = %d %d  d/tsamp = %f %f\n",
  //m,n,start,period,dsamp,tsamp);

  // Determine start and duration for averaging
  if (dsamp == tsamp){
    t0 = start;
    tn = period;
  }else{
    t0 = (int) ((float)start  * dsamp/tsamp);
    tn = (int) ((float)period * dsamp/tsamp);
  }
  if ((t0 < 0) || (t0+tn > n))
    exit_error("NDATA_GET_MEAN_SDEV_2D_FDATA","Bad start or period");
  
  t = (float *)myalloc(m*sizeof(float));
  for(i=0;i<m;i++)
    t[i] =  mean_farray(fdata[i]+t0,tn);

  mean_sdev_farray(t,m,rmean,rsd);


 // 2019 Aug 26 ***** WYETH, WHERE DOES 't' get freed ?????????
 // ************ WYETH, WHERE DOES 't' get freed ?????????
 // ************ WYETH, WHERE DOES 't' get freed ?????????
 // ************ WYETH, WHERE DOES 't' get freed ?????????
  myfree(t);  // WYETH added 2019 Aug 26
}
/**************************************-**************************************/
/*                                                                           */
/*                           NDATA_GET_STAT_2D_FDATA                         */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Assume fdata starts at time 0.                                         */
/*  - There are 'm' data traces, each of time length 'n'                     */
/*  - First, compute the mean of each trace.                                 */
/*  - Second, compute stats across the 'm' mean values.                      */
/*                                                                           */
/*****************************************************************************/
void ndata_get_stat_2d_fdata(fdata,m,n,start,period,dsamp,tsamp,rmu,rsd,rmin,
			     rmax)
     float **fdata;
     int m,n,start,period;
     float dsamp,tsamp;    // dsamp applies to data, tsamp applies to times
     float *rmu,*rsd;
     float *rmin,*rmax;
{
  int i;
  int t0,tn;
  float *t;

  // Determine start and duration for averaging
  if (dsamp == tsamp){
    t0 = start;
    tn = period;
  }else{
    t0 = (int) ((float)start  * dsamp/tsamp);
    tn = (int) ((float)period * dsamp/tsamp);
  }
  if ((t0 < 0) || (t0+tn > n))
    exit_error("NDATA_GET_STAT_2D_FDATA","Bad start or period");
  
  t = (float *)myalloc(m*sizeof(float));
  for(i=0;i<m;i++)
    t[i] =  mean_farray(fdata[i]+t0,tn);

  mean_sdev_farray(t,m,rmu,rsd);
  get_min_max_farray(t,m,rmin,rmax);

  myfree(t);
}
/**************************************-**************************************/
/*                                                                           */
/*                           NDATA_UTIL_CHOICE_GEN                           */
/*                                                                           */
/*  Generate a choice for this trial.                                        */
/*  Return                                                                   */
/*    0 - test value <= 0                                                    */
/*    1 - test value  > 0                                                    */
/*   -1 - no decision                                                        */
/*                                                                           */
/*  'ctype'                                                                  */
/*    chan_diff     difference of two channels                               */
/*    prefix_diff   difference of two sets of channels, names are prefixes   */
/*                                                                           */
/*****************************************************************************/
int ndata_util_choice_gen(t,ctype,t0,tn,name_add,name_sub,rval,rvadd,rvsub)
     struct ndtrial_struct *t;  // Pointer to trial
     char *ctype;               // Calculation type
     int t0,tn;                 // Start and duration for response window
     char *name_add;            // Add response for this channel/pop
     char *name_sub;            // Subtract response for this channel/pop
     float *rval;               // Return the decision value
     float *rvadd;              // Return the "add" value
     float *rvsub;              // Return the "sub" value
{
  int i;
  int *s,n,d,cnum,nc,nadd,nsub;
  float tval,tadd,tsub,*tlist;
  char *tname;

  d = -1;  // No decision

  if (strcmp(ctype,"chan_diff")==0){
    cnum = ndata_trial_get_record_index(t,name_add);
    if (cnum < 0){
      printf("  *** Channel name: %s\n",name_add);
      exit_error("NDATA_UTIL_CHOICE_GEN","Could not find channel in trial");
    }else{
      s = t->r[cnum].p;
      n = t->r[cnum].n;
      tadd = (float)count_spikes(s,n,t0,tn);
    }

    cnum = ndata_trial_get_record_index(t,name_sub);
    if (cnum < 0){
      printf("  *** Channel name: %s\n",name_sub);
      exit_error("NDATA_UTIL_CHOICE_GEN","Could not find channel in trial");
    }else{
      s = t->r[cnum].p;
      n = t->r[cnum].n;
      tsub = (float)count_spikes(s,n,t0,tn);
    }

  }else if (strcmp(ctype,"prefix_diff")==0){

    if (compare_prefix_string(name_add,name_sub) == 1)
      exit_error("NDATA_UTIL_CHOICE_GEN",
		 "'choice_add' or 'choice_sub' is a prefix of the other");

    tadd = tsub = 0.0;
    nadd = nsub = 0;

    nc = t->nrec;
    for(i=0;i<nc;i++){
      //
      //  Check the name of each channel in this trial
      //
      tname = t->r[i].name;
      if (compare_prefix_string_order(name_add,tname) == 1){
	tadd += (float)count_spikes(t->r[i].p,t->r[i].n,t0,tn);
	nadd += 1;
	//printf("ADDing %s\n",tname);
      }else if (compare_prefix_string_order(name_sub,tname) == 1){
	tsub += (float)count_spikes(t->r[i].p,t->r[i].n,t0,tn);
	nsub += 1;
	//printf("Subtracting   %s\n",tname);
      }
    }

    if (nadd == 0){
      exit_error("NDATA_UTIL_CHOICE_GEN",
		 "Prefix for 'choice_add' has no matches");
    }
    if (nsub == 0){
      exit_error("NDATA_UTIL_CHOICE_GEN",
		 "Prefix for 'choice_sub' has no matches");
    }

  }else{
    exit_error("NDATA_UTIL_CHOICE_GEN","Unknown choice calculation");
  }

  tval = tadd - tsub;
  if (tval > 0.0)
    d = 1;
  else
    d = 0;

  *rval  = tval;
  *rvadd = tadd;
  *rvsub = tsub;

  return d;
}

/**************************************-**************************************/
/*                                                                           */
/*                            NDATA_UTIL_CHOICE_01                           */
/*                                                                           */
/*  Add the choice variable name 'ch_name' to the 'nd' structure.            */
/*                                                                           */
/*  'ctype'                                                                  */
/*    chan_diff     difference of two channels                               */
/*    prefix_diff   difference of two sets of channels, names are prefixes   */
/*                                                                           */
/*****************************************************************************/
void ndata_util_choice_01(nd,ctype,t0,tn,name_add,name_sub,ch_name,
			  rval,rvadd,rvsub)
     struct ndata_struct *nd;
     char *ctype;               // Calculation type
     int t0,tn;                 // Start and duration for response window
     char *name_add;            // Add response for this channel/pop
     char *name_sub;            // Subtract response for this channel/pop
     char *ch_name;             // Choice variable name
     float **rval;              // [ntr] return decision values
     float **rvadd;             // [ntr] return "add" value
     float **rvsub;             // [ntr] return "sub" value
{
  int i;
  int ntr,*d;
  float *dval,*dvadd,*dvsub,fval,vadd,vsub;
  char tstr[SLEN];
  struct ndtrial_struct *tp;  // Pointer to trial

  ndata_add_var_param(nd,ch_name,'i',"-1");  // Create a new var param

  ntr = nd->ntrial;
  d = (int *)myalloc(ntr*sizeof(int));
  dval  = get_farray(ntr);
  dvadd = get_farray(ntr);
  dvsub = get_farray(ntr);
  for(i=0;i<ntr;i++){ // Set value of var param on each trial

    tp = &(nd->t[i]);  // Trial pointer

    // Get the decision: 0,1,-1
    d[i] = ndata_util_choice_gen(tp,ctype,t0,tn,name_add,name_sub,
				 &fval,&vadd,&vsub);
    dval[i]  = fval;
    dvadd[i] = vadd;
    dvsub[i] = vsub;
    sprintf(tstr,"%i",d[i]);
    ndata_set_var_value_trial(nd,i,ch_name,tstr,1); // 1-freeflag
  }

  *rval  = dval;
  *rvadd = dvadd;
  *rvsub = dvsub;
}
