/*****************************************************************************/
/*                                                                           */
/*  mctrl_util.c                                                             */
/*  wyeth bair                                                               */
/*                                                                           */
/*  For use with mcontrol.c, mfilter.c, etc.                                 */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#include "my_util.h"
#include "carray_util.h"
#include "paramfile_util.h"
#include "mctrl_util.h"

int DIR_FLAG = 0;

//
//  NOTE: These are EXTERNed in 'mctrl_util.h'
//
char MCTRL_DIR_ADMIN[SLEN];
char MCTRL_DIR_USR[SLEN];
char MCTRL_DIR_WORK[SLEN];       // For running wm/mmm jobs.
char MCTRL_DIR_WORK_STIM[SLEN];
char MCTRL_DIR_MM2[SLEN];
char MCTRL_DIR_M[SLEN];
//REMOVE char MCTRL_DIR_M_MOD[SLEN];
char MCTRL_DIR_MM[SLEN];
char MCTRL_DIR_M_RSP[SLEN];
char MCTRL_DIR_M_STM[SLEN];
char MCTRL_DIR_M_UPT[SLEN];

/**************************************-**************************************/
/*                                                                           */
/*                       MCTRL_UTIL_DETERMINE_DIR_STRUCTURE                  */
/*                                                                           */
/*  Write the process ID to the file.                                        */
/*                                                                           */
/*****************************************************************************/
void mctrl_util_determine_dir_structure()
{
  //
  //  WYETH - NOTE, if this is called by 'imodgen', it might not be user
  //  'wyeth', and thus the user would need permission on MAIN_DIR_ADMIN
  //  So, perhaps we want to think about whether some other file should be
  //  tested?  But probably OK, to have the permissions on the directory,
  //  but keep the files unreadable to others.  Just don't want them to be
  //  open to the web.
  //
  if (file_exists(MAIN_DIR_WORK)){
    //printf("  Using Main Linux directory structure.\n");
    strcpy(MCTRL_DIR_ADMIN,MAIN_DIR_ADMIN);
    strcpy(MCTRL_DIR_USR  ,MAIN_DIR_USR);
    strcpy(MCTRL_DIR_WORK ,MAIN_DIR_WORK);
    strcpy(MCTRL_DIR_WORK_STIM,MAIN_DIR_WORK_STIM);
    strcpy(MCTRL_DIR_MM2  ,MAIN_DIR_MM2);
    strcpy(MCTRL_DIR_M,    MAIN_DIR_M);
    //REMOVE strcpy(MCTRL_DIR_M_MOD,MAIN_DIR_M_MOD);
    strcpy(MCTRL_DIR_MM,   MAIN_DIR_MM);
    strcpy(MCTRL_DIR_M_RSP,MAIN_DIR_M_RSP);
    strcpy(MCTRL_DIR_M_STM,MAIN_DIR_M_STM);
    strcpy(MCTRL_DIR_M_UPT,MAIN_DIR_M_UPT);
    //mcontrol_log("mcontrol started.  Using main directory structure\n",0);
  }else{
    printf("  Using Apple directory structure.\n");
    strcpy(MCTRL_DIR_ADMIN,MAC_DIR_ADMIN);
    strcpy(MCTRL_DIR_USR  ,MAC_DIR_USR);
    strcpy(MCTRL_DIR_WORK ,MAC_DIR_WORK);
    strcpy(MCTRL_DIR_WORK_STIM,MAC_DIR_WORK_STIM);
    strcpy(MCTRL_DIR_MM2  ,MAC_DIR_MM2);
    strcpy(MCTRL_DIR_M,    MAC_DIR_M);
    //REMOVE strcpy(MCTRL_DIR_M_MOD,MAC_DIR_M_MOD);
    strcpy(MCTRL_DIR_MM,   MAC_DIR_MM);
    strcpy(MCTRL_DIR_M_RSP,MAC_DIR_M_RSP);
    strcpy(MCTRL_DIR_M_STM,MAC_DIR_M_STM);
    strcpy(MCTRL_DIR_M_UPT,MAC_DIR_M_UPT);
    //mcontrol_log("mcontrol started.  Using Apple directory structure\n",0);
  }
  DIR_FLAG = 1;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MCTRL_UTIL_REPORT_PID_RUN                         */
/*                                                                           */
/*  Write the process ID and the parent process ID to the file.              */
/*                                                                           */
/*  Called by 'mm' or 'wm', for use w/ 'mcontrol.c'                          */
/*                                                                           */
/*****************************************************************************/
void mctrl_util_report_pid_run()
{
  int pid,ppid;
  char tstr[SLEN],outfile[SLEN];

  pid = getpid();
  ppid = getppid();  // Parent is probably a shell, not who called system

  //
  //  *** WYETH - MAKE THIS WRITE LOCALLY - AND CHANGE THE DONE CONDITION
  //  *** THIS BEING DONE FOR CHANGE - OVER TO KAMBPLIPOOCHI.
  //
  //sprintf(outfile,"%s/%s",MCTRL_DIR_ADMIN,FILE_ADM_PID_JOB);
  sprintf(outfile,"%s",FILE_ADM_PID_JOB);
  sprintf(tstr,"%d %d\n",pid,ppid);
  //append_string_to_file(outfile,tstr);
  append_string_to_file(outfile,tstr);
}
/**************************************-**************************************/
/*                                                                           */
/*                         MCTRL_UTIL_MODLIST_GET_ONODE                      */
/*                                                                           */
/*  Write the process ID and the parent process ID to the file.              */
/*                                                                           */
/*  For use w/ 'mcontrol.c'                                                  */
/*                                                                           */
/*****************************************************************************/
struct onode *mctrl_util_modlist_get_onode(modname,errstr)
     char *modname;
     char **errstr;  // Return an error string, or NULL
{
  char modlist_file[SLEN3],terr[LONG_SLEN];
  struct onode *omods,*mo,*o;

  //
  //  Attempt to read the main ModelList.txt onode file.
  //
  sprintf(modlist_file,"%s/%s",MCTRL_DIR_MM,FILE_MOD_MODLIST);
  omods = paramfile_onode_file_read(modlist_file,0);
  if (omods == NULL){
    sprintf(terr,"MCTRL_UTIL_MODLIST_GET_ONODE  Model List Not Found:  %s",
	    modlist_file);
    *errstr = strdup(terr);
    return NULL;
  }

  //
  //  Find the node for the specified model 
  //
  mo = onode_get_node_type_item_val(omods,"model","name",modname);
  if (mo == NULL){
    sprintf(terr,"MCTRL_UTIL_MODLIST_GET_ONODE  Model node not found:  %s",
	    modname);
    *errstr = strdup(terr);
    return NULL;
  }

  o = onode_copy(mo);  // Keep a copy of the onode for this model

  paramfile_onode_free(omods);  // Free the onode list

  *errstr = NULL;
  return o;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MCTRL_UTIL_CMD_LINE_PARS                         */
/*                                                                           */
/*  Return a string with command line parameters, or the empty string "".    */
/*                                                                           */
/*  *** For now, this only handles 'tn'                                      */
/*                                                                           */
/*****************************************************************************/
char *mctrl_util_cmd_line_pars(tn,pstr)
     int tn;
     char *pstr;
{
  char sstr[LONG_SLEN],fstr[LONG_SLEN*2],*tstr;


  if (tn > 0){
    sprintf(sstr,"tn %d",tn);
  }else
    strcpy(sstr,"");  // Empty string


  if (strcmp(pstr,"none")!=0){
    if (strlen(sstr) > 0)
      sprintf(fstr,"%s %s",pstr,sstr);
    else
      sprintf(fstr,"%s",pstr);
  }else{
    sprintf(fstr,"%s",sstr);
  }

  tstr = strdup(fstr);

  return tstr;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MCTRL_UTIL_UF_FI_READ                           */
/*                                                                           */
/*  Read the contents of a 'file_list' file.                                 */
/*                                                                           */
/*****************************************************************************/
void mctrl_util_uf_fi_read(infile,rn,ra,ru)
     char *infile;  // 'file_list' to read
     int *rn;       // Number of entries
     char ***ra;    // List of archive file names [*rn]
     char ***ru;    // List of user file names [*rn]
{
  FILE *fopen(),*fin;
  int i,n,ns;
  char **sa,**su,ta[SLEN],tu[SLEN];

  if ((fin = fopen(infile,"r"))==NULL){
    printf("  *** File %s\n",infile);
    exit_error("MCTRL_UTIL_UF_FI_READ","Can't open file");
  }

  ns = fscanf(fin,"%d",&n);  // Read number of entries in table

  if (n <= 0){
    sa = su = (char **)NULL;
  }else{
    sa = (char **)myalloc(n*sizeof(char *));
    su = (char **)myalloc(n*sizeof(char *));

    for(i=0;i<n;i++){
      ns = fscanf(fin,"%s %s",ta,tu);
      sa[i] = strdup(ta);
      su[i] = strdup(tu);
    }
  }

  fclose(fin);

  *rn = n;
  *ra = sa;
  *ru = su;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MCTRL_UTIL_UF_FMAN_GET                          */
/*                                                                           */
/*  Similar to iModel PHP .../m/include/imod_uf_fman_get.php.                */
/*                                                                           */
/*  Return information from the user file manager file.                      */
/*                                                                           */
/*****************************************************************************/
char *mctrl_util_uf_fman_get(uname,ftype,field)
     char *uname;  // User name
     char *ftype;  // file type, e.g., model, stim, resp, frameset
     char *field;  // Field name
{
  int i,k;
  int dn,ns,fi,done;
  char *infile,**data,**slist,*sval;

  //printf("  MCTRL_UTIL_UF_FMAN_GET  uname %s  ftype %s  field %s\n",
  //uname,ftype,field);

  infile = mctrl_util_uf_fname("fman",uname,NULL,NULL);

  //printf("INFILE:  %s\n",infile);

  //  Read the file_manager as a list of strings
  read_2d_carray(infile,1,0,&data,&dn);

  //
  //  Parse the field names in the first line, and get the field index
  //
  get_items_from_string(data[0],&slist,&ns);
  fi = search_2d_carray(slist,field,ns);
  free_2d_carray(slist,ns);
  if (k < 0){
    printf("  *** field:  %s\n",field);
    exit_error("MCTRL_UTIL_UF_FMAN_GET","Could not find field name");
  }

  //printf("fi = %d\n",fi);

  //
  //  Read the data lines until first string matches 'ftype'
  //
  k = -1;
  done = 0;
  i = 1;
  sval = NULL;
  while(done == 0){

    get_items_from_string(data[i],&slist,&ns);
    if (strcmp(slist[0],ftype)==0){
      sval = strdup(slist[fi]);
      done = 1;
    }else{
      i += 1;
      if (i >= dn)
	done = 1;
    }
    free_2d_carray(slist,ns);
  }

  /*
  if (sval != NULL)
    printf("SVAL:....%s\n",sval);
  else
    printf("SVAL is NULL *****\n");
  */
    

  myfree(infile);
  free_2d_carray(data,dn);

  return sval;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MCTRL_UTIL_UF_FNAME                            */
/*                                                                           */
/*  Similar to iModel PHP .../m/include/imod_uf_fname.php.                   */
/*                                                                           */
/*  Return the full path filename relatd to the UF (User filesystem).        */
/*                                                                           */
/*****************************************************************************/
char *mctrl_util_uf_fname(mode,uname,ftype,tname)
     char *mode;   // E.g., fi, fman, arch, aname
     char *uname;  // E.g., lab
     char *ftype;  // E.g., model, stim, resp, frameset
     char *tname;  // User defined file name
{
  int k;
  int n;
  char *locpath,*rstr,fname[LONG_SLEN],**alist,**ulist;

  //printf("  MCTRL_UTIL_UF_FNAME  mode %s  uname %s  ftype %s  tname %s\n",
  //mode,uname,ftype,tname);

  rstr = NULL;  // String to return

  if ((strcmp(mode,"fi")==0) || (strcmp(mode,"aname")==0)){
    //
    //  Get full name of 'file_index' for this user, this type of file
    //
    locpath = mctrl_util_uf_fman_get(uname,ftype,"Path");
    sprintf(fname,"%s/%s/uf/%s/file_index",MCTRL_DIR_USR,uname,locpath);

    if (strcmp(mode,"fi")==0){
      rstr = strdup(fname);
    }else if (strcmp(mode,"aname")==0){
      //
      //  Get FULL ARCHIVE name of the USER file from the archive
      //   e.g. .../f0001
      //
      mctrl_util_uf_fi_read(fname,&n,&alist,&ulist);
      //printf("  (from...uf_fi_read())  n = %d\n",n);


      k = search_2d_carray(ulist,tname,n);
      if (k >= 0){

	sprintf(fname,"%s/%s/uf/%s/%s",MCTRL_DIR_USR,uname,locpath,alist[k]);
	rstr = strdup(fname);
      }
      //printf("  mctrl_util_uf_fname:  k = %d\n",k);

      free_2d_carray(alist,n);
      free_2d_carray(ulist,n);

      if (locpath != NULL)
	myfree(locpath);
    }

  }else if (strcmp(mode,"fman")==0){
    //
    //  Get full name of 'file_manager' for this user
    //
    sprintf(fname,"%s/%s/uf/file_manager",MCTRL_DIR_USR,uname);
    rstr = strdup(fname);

  }else if (strcmp(mode,"arch")==0){
    //
    //  Get full name of a file from the archive
    //
    locpath = mctrl_util_uf_fman_get(uname,ftype,"Path");
    //printf("_____________locpath = %s\n",locpath);
    sprintf(fname,"%s/%s/uf/%s/%s",MCTRL_DIR_USR,uname,locpath,tname);
    if (locpath != NULL)
      myfree(locpath);

    rstr = strdup(fname);

  }else{
    printf("  *** mode = %s\n",mode);
    exit_error("MCTRL_UTIL_UF_FNAME","Unknown 'mode' value");
  }

  return rstr;
}
