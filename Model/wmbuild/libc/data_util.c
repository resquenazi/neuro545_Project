/*****************************************************************************/
/*                                                                           */
/*  data_util.c                                                              */
/*  Wyeth Bair                                                               */
/*                                                                           */
/*  This file contains routines related to processing data files and         */
/*  performing database types of operations.                                 */
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
#include "wytif.h"

/**************************************-**************************************/
/*                                                                           */
/*                           LOOKUP_KEY_FIELD_FILE                           */
/*                                                                           */
/*  Make a lot of assumptions, this is quick and dirty.  Search for the      */
/*  record that has "key" in the first field.  Return the value of the       */
/*  field named "field" for that record.                                     */
/*                                                                           */
/*  File structure:                                                          */
/*  <nrec> <nfield>                                                          */
/*  <field_name_1> ... <field_name_n>                                        */
/*  ...             ...                                                      */
/*  ... DATA MATRIX ...                                                      */
/*  ...             ...                                                      */
/*                                                                           */
/*****************************************************************************/
int lookup_key_field_file(infile,key,field,rstring)
     char infile[],key[],field[],**rstring;
{
  FILE *fopen(),*fin;
  int i,j;
  int nfield,nrec,xfield,xrec,flag,ns;
  char ***data,**fname,temp[SLEN],*tc;

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** File %s\n",infile);
    exit_error("LOOKUP_KEY_FIELD_FILE","Cannot open file\n");
  }

  /*** WYETH - added on 03/11/01 to allow comments at beginning of file ***/
  ns = fscanf(fin,"%s",temp);
  while(temp[0] == '#'){
    tc = fgets(temp,SLEN,fin);
    ns = fscanf(fin,"%s",temp);
  }
  
  /*** Read the number of records and fields. ***/
  sscanf(temp,"%d",&nrec);
  ns = fscanf(fin,"%d",&nfield);
  //printf("nrec= %d  nfield= %d\n",nrec,nfield);
  data = get_2d_pointer_carray(nfield,nrec);
  fname = (char **)myalloc(nfield*sizeof(char *));

  for(j=0;j<nfield;j++){
    ns = fscanf(fin,"%s",temp);
    fname[j] = strdup(temp);
  }

  for(i=0;i<nrec;i++){
    for(j=0;j<nfield;j++){
      ns = fscanf(fin,"%s",temp);
      data[j][i] = strdup(temp);
    }
    //printf("%s\n",data[0][i]);
  }
  fclose(fin);

  flag = 1;
  xrec = -1;
  xfield = search_2d_carray(fname,field,nfield);
  if (xfield < 0){
    printf("  *** WARNING LOOKUP_KEY_FIELD_FILE Field name %s not found\n",
	   field);
    flag = -1;
  }else{
    xrec = search_2d_carray(data[0],key,nrec);
    if (xrec < 0){
      printf(" *** WARNING LOOKUP_KEY_FIELD_FILE Key %s not found\n",key);
      flag = 0;
    }
  }
  if (flag > 0)
    *rstring = strdup(data[xfield][xrec]);

  free_2d_pointer_carray(data,nfield,nrec);
  free_2d_carray(fname,nfield);

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                         DATA_UTIL_T11_INT_FOR_RGB                         */
/*                                                                           */
/*****************************************************************************/
int data_util_t11_int_for_rgb(r,g,b)
     float r,g,b;  // color values [0..1]
{
  int ri,gi,bi,c;

  ri = my_rint(r * 255.0);
  gi = my_rint(g * 255.0);
  bi = my_rint(b * 255.0);

  c = ri;
  c = c << 8;
  c = c | gi;
  c = c << 8;
  c = c | bi;

  return c;
}
/**************************************-**************************************/
/*                                                                           */
/*                            DATA_UTIL_T11_GET_RGB                          */
/*                                                                           */
/*  Return the R, G, and B components as ints from the tcode 11 data.        */
/*                                                                           */
/*****************************************************************************/
void data_util_t11_get_rgb(d,rr,rg,rb)
     int d;
     int *rr,*rg,*rb;
{
  *rb = d & 0377;
  *rg = (d>>8) & 0377;
  *rr = (d>>16) & 0377;

  //printf("  RGB = %d %d %d\n",*rr,*rg,*rb);
  //exit(0);
}
/**************************************-**************************************/
/*                                                                           */
/*                         DATA_UTIL_T11_GET_RGB_FLOAT                       */
/*                                                                           */
/*  Return the R, G, and B components as floats from the tcode 11 data.      */
/*                                                                           */
/*****************************************************************************/
void data_util_t11_get_rgb_float(d,rr,rg,rb)
     int d;
     float *rr,*rg,*rb;
{
  *rb = (float)(d & 0377);
  *rg = (float)((d>>8) & 0377);
  *rr = (float)((d>>16) & 0377);
}
/**************************************-**************************************/
/*                                                                           */
/*                          DATA_UTIL_T11_GET_RGB_01                         */
/*                                                                           */
/*  Return the R, G, and B components in [0,1].                              */
/*                                                                           */
/*****************************************************************************/
void data_util_t11_get_rgb_01(d,rr,rg,rb)
     int d;
     float *rr,*rg,*rb;
{
  *rb = (float)(d       & 0377) / 255.0;
  *rg = (float)((d>>8)  & 0377) / 255.0;
  *rr = (float)((d>>16) & 0377) / 255.0;
}
/**************************************-**************************************/
/*                                                                           */
/*                      DATA_UTIL_T11_GET_MIN_MAX_RGB_ALL                    */
/*                                                                           */
/*  Return the minimum and maximum values among the R, G, B components.      */
/*                                                                           */
/*****************************************************************************/
void data_util_t11_get_min_max_rgb_all(d,rmin,rmax)
     int d;
     int *rmin,*rmax;
{
  int r,g,b,min,max;

  data_util_t11_get_rgb(d,&r,&g,&b);

  min = r;
  max = r;

  if (g < min)
    min = g;
  else if (g > max)
    max = g;

  if (b < min)
    min = b;
  else if (b > max)
    max = b;

  *rmin = min;
  *rmax = max;
}
/**************************************-**************************************/
/*                                                                           */
/*                          DATA_UTIL_T11_2D_FROM_RGB                        */
/*                                                                           */
/*  Convert R,G,B 2D arrays into T11 color image.                            */
/*                                                                           */
/*****************************************************************************/
int **data_util_t11_2d_from_rgb(dr,dg,db,xn,yn)
     float **dr,**dg,**db;   // [xn][yn] red, green, blue values
                             //    MUST BE IN RANGE [0..1]
     int xn,yn;              // image size
{
  int i,j;
  int **d;

  d = get_2d_iarray(xn,yn);

  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++)
      d[i][j] = data_util_t11_int_for_rgb(dr[i][j],dg[i][j],db[i][j]);

  return d;
}
/**************************************-**************************************/
/*                                                                           */
/*                            DATA_UTIL_3RGB_TO_3D                           */
/*                                                                           */
/*  Convert a 4D array that has a set of RGB float images into a 3D array    */
/*  where the 3rd (time) dimension interleaves the R,G,B planes for each     */
/*  image in sequence.                                                       */
/*                                                                           */
/*****************************************************************************/
float ***data_util_3rgb_to_3d(data,xn,yn,tn)
     float ****data;   // [3][xn][yn][tn]  RGB is first index
     int xn,yn,tn;
{
  int i,j,k;
  int tt;
  float ***s;

  s = get_3d_farray(xn,yn,3*tn);

  for(i=0;i<xn;i++){
    for(j=0;j<yn;j++){
      tt = 0;
      for(k=0;k<tn;k++){
	s[i][j][tt] = data[0][i][j][k];
	tt += 1;
	s[i][j][tt] = data[1][i][j][k];
	tt += 1;
	s[i][j][tt] = data[2][i][j][k];
	tt += 1;
      }
    }
  }

  return s;
}
/**************************************-**************************************/
/*                                                                           */
/*                          DATA_UTIL_3RGB_RAW1_READ                         */
/*                                                                           */
/*  Read a single image using Dean's simple text format for AlexNet          */
/*  images.                                                                  */
/*                                                                           */
/*****************************************************************************/
float ***data_util_3rgb_raw1_read(infile,rzn,rxn,ryn)
     char *infile;
     int *rzn;      // Note, returned array is [3][x][y]
     int *rxn;
     int *ryn;
{
  FILE *fopen(),*fin;
  int i,j,k;
  int ns,xn,yn,zn;
  float ***data;

  if ((fin = fopen(infile,"r"))==NULL){
    printf("  File:  %s\n",infile);
    exit_error("DATA_UTIL_3RGB_RAW1_READ","Cannot open file");
  }

  ns = fscanf(fin,"%d %d %d",&zn,&xn,&yn); // e.g.:  "3 227 227"

  //printf("    Image size:  %d %d %d\n",xn,yn,zn);

  if (zn != 3)
    exit_error("DATA_UTIL_3RGB_RAW1_READ","Image depth mis-match");

  data = get_3d_farray(zn,xn,yn);

  // *** NOTE ***
  //     The *BLUE* image data comes first, starting from upper left pixel
  //     then the *GREEN* image data comes next
  //     then the *RED* image data comes last.
  //     This is the case in Dean's simple file format.

  for(i=(zn-1);i>=0;i--){    // For R,G,B planes (reverse for R,G,B order)
    for(k=(yn-1);k>=0;k--){   //   For each row (starting with top row),
      for(j=0;j<xn;j++){      //     read from left to right.
	ns = fscanf(fin,"%f",&(data[i][j][k]));
      }
    }
  }
  fclose(fin);

  *rxn = xn;
  *ryn = yn;
  *rzn = zn;

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                             DATA_UTIL_FST_WRITE                           */
/*                                                                           */
/*  Write a frame set file, to be used for visual stimuli.                   */
/*                                                                           */
/*****************************************************************************/
void data_util_fst_write(outfile,data,x0,xn,y0,yn,t0,tn,drange,fcode,binoc)
     char outfile[];   // Output file name, should be .fst
     float ***data;    // [w][h][t]
     int x0,xn,y0,yn,t0,tn;    // width, height, time duration
     int drange;       // 0-[0..1], 1-[-1..1]
     int fcode;        // 1-unsigned char, 2-unsigned short, 3-int, 4-float
     int binoc;        // 0-monoc, 2-binoc L/R frame interleave, 3-L before R
{
  int i,j,k;
  FILE *fopen(),*fout;
  int temp,nbytes,vcode,color_code;
  unsigned short int ds;

  printf("  DATA_UTIL_FST_WRITE\n");
  printf("    xn,yn,tn = %d %d %d\n",xn,yn,tn);
  printf("    writing %s\n",outfile);

  if (fcode != 2)
    exit_error("DATA_UTIL_FST_WRITE","fcode expected to be 2");

  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** File %s\n",outfile);
    exit_error("DATA_UTIL_FST_WRITE","Cannot open file");
  }

  temp = 16909061; // This is hex 01 02 03 05, for checking byte order
  vcode = 1;       // Version 1
  color_code = 0;  // Achromatic

  fwrite((char *)&temp,sizeof(int),1,fout);
  fwrite((char *)&vcode,sizeof(int),1,fout);
  fwrite((char *)&xn,sizeof(int),1,fout);
  fwrite((char *)&yn,sizeof(int),1,fout);
  fwrite((char *)&tn,sizeof(int),1,fout);
  fwrite((char *)&fcode,sizeof(int),1,fout);
  fwrite((char *)&color_code,sizeof(int),1,fout);
  fwrite((char *)&binoc,sizeof(int),1,fout);

  if (fcode == 1)
    nbytes = 1;
  else if (fcode == 2)
    nbytes = 2;
  else if ((fcode == 3) || (fcode == 4))
    nbytes = 4;
  else
    exit_error("DATA_UTIL_FST_WRITE","Invalid fcode");

  for(k=0;k<tn;k++){
    for(j=0;j<yn;j++){
      for(i=0;i<xn;i++){
	// 65535 is used (rather than 65536) to avoid turning 1 into 0
	ds = (unsigned short int)(data[i+x0][j+y0][k+t0] * 65535);
	fwrite(&ds,nbytes,1,fout);
      }
    }
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                             DATA_UTIL_FST_READ                            */
/*                                                                           */
/*  Read a frameset (.fst) file, to be used for visual stimuli.              */
/*                                                                           */
/*  *** Data returned in NumRec tensor format.                               */
/*                                                                           */
/*****************************************************************************/
float ***data_util_fst_read(logfile,path,infile,drange,rxn,ryn,rtn,
			    rfcode,rbinoc,rccode)
     char *logfile;      // file name to print output, or NULL to print
     char *path;         // extension to add to 'infile', or NULL
     char *infile;       // Input file name, typically .fst
     int drange;         // 0-[0..1], 1-[-1..1]
     int *rxn,*ryn,*rtn; // width, height, time duration
     int *rfcode;        // 1-unsigned char, 2-unsigned short, 3-int, 4-float
     int *rbinoc;        // 0-monoc, 2-binoc L/R frame interleave, 3-L before R
     int *rccode;        // color code
{
  FILE *fopen(),*fin;
  int i,j,k;
  int nr,temp,nbytes,vcode,revflag;
  unsigned short int ds;
  float ***data;
  char ggstr[SLEN3],fname[SLEN];

  if (path == NULL){
    //printf("  Path is null\n");
    sprintf(fname,"%s",infile);
  }else{
    //printf("  Path is **NOT** null\n");
    sprintf(fname,"%s%s",path,infile);
  }

  mylog(logfile,"  DATA_UTIL_FST_READ\n");
  sprintf(ggstr,"    Reading  %s\n",fname);
  mylog(logfile,ggstr);

  if ((fin = fopen(fname,"r")) == NULL){
    printf("  *** File %s\n",fname);
    exit_error("DATA_UTIL_FST_READ","Cannot open file");
  }

  nr = fread((char *)&temp,sizeof(int),1,fin);
  if (temp == 16909061)
    revflag = 0;
  else
    revflag = 1;

  revflag_fread((char *)(&vcode),sizeof(int),1,fin,revflag);
  revflag_fread((char *)rxn,sizeof(int),1,fin,revflag);
  revflag_fread((char *)ryn,sizeof(int),1,fin,revflag);
  revflag_fread((char *)rtn,sizeof(int),1,fin,revflag);
  revflag_fread((char *)rfcode,sizeof(int),1,fin,revflag);
  revflag_fread((char *)rccode,sizeof(int),1,fin,revflag);
  revflag_fread((char *)rbinoc,sizeof(int),1,fin,revflag);

  sprintf(ggstr,"    %d x %d x %d (binoc %d, fcode %d, color %d, range %d)\n",
	  *rxn,*ryn,*rtn,*rbinoc,*rfcode,*rccode,drange);
  mylog(logfile,ggstr);

  if (*rfcode != 2)
    exit_error("DATA_UTIL_FST_READ","fcode expected to be 2");
  if (*rccode != 0)
    exit_error("DATA_UTIL_FST_READ","color_code expected to be 0");

  if (*rfcode == 1)
    nbytes = 1;
  else if (*rfcode == 2)
    nbytes = 2;
  else if ((*rfcode == 3) || (*rfcode == 4))
    nbytes = 4;
  else
    exit_error("DATA_UTIL_FST_READ","Invalid fcode");

  data = f3tensor(1,*rxn,1,*ryn,1,*rtn);

  for(k=0;k<*rtn;k++){
    for(j=0;j<*ryn;j++){
      for(i=0;i<*rxn;i++){
	revflag_fread((char *)&ds,nbytes,1,fin,revflag);

	if (drange == 0)
	  data[i+1][j+1][k+1] = (float)ds / 65535.0;
	else
	  data[i+1][j+1][k+1] = 2.0*((float)ds / 65535.0) - 1.0;

	// 65535 is used (rather than 65536) to avoid turning 1 into 0
	//ds = (unsigned short int)(data[i+x0][j+y0][k+t0] * 65535);
	//fwrite(&ds,nbytes,1,fout);
      }
    }
  }
  fclose(fin);

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                           DATA_UTIL_FST_READ_TXT                          */
/*                                                                           */
/*  Read a frameset text file, to be used for visual stimuli.                */
/*                                                                           */
/*  *** Data returned in NumRec tensor format.                               */
/*                                                                           */
/*****************************************************************************/
float ***data_util_fst_read_txt(logfile,path,infile,drange,yrev,rxn,ryn,rtn,
				rfcode,rbinoc,rccode)
     char *logfile;      // file name to print output, or NULL to print
     char *path;         // extension to add to 'infile', or NULL
     char *infile;       // Input file name, typically .fst
     int drange;         // 0-[0..1], 1-[-1..1]
     int yrev;           // 1-reverse the y-axis
     int *rxn,*ryn,*rtn; // width, height, time duration
     int *rfcode;        // 1-unsigned char, 2-unsigned short, 3-int, 4-float
     int *rbinoc;        // 0-monoc, 2-binoc L/R frame interleave, 3-L before R
     int *rccode;        // color code
{
  FILE *fopen(),*fin;
  int i,j,k;
  int ns,temp,dzn,xn,yn,tn;
  unsigned short int ds;
  float ***data;
  char ggstr[SLEN3],fname[SLEN],ts1[SLEN],ts2[SLEN];

  if (path == NULL){
    //printf("  Path is null\n");
    sprintf(fname,"%s",infile);
  }else{
    //printf("  Path is **NOT** null\n");
    sprintf(fname,"%s%s",path,infile);
  }

  mylog(logfile,"  DATA_UTIL_FST_READ_TXT\n");
  sprintf(ggstr,"    Reading  %s\n",fname);
  mylog(logfile,ggstr);

  if ((fin = fopen(fname,"r")) == NULL){
    printf("  *** File %s\n",fname);
    exit_error("DATA_UTIL_FST_READ_TXT","Cannot open file");
  }

  //
  //  Read from the open file
  //
  ns = fscanf(fin,"%s %s",ts1,ts2);

  if (strcmp(ts1,"frameset")!=0){
    exit_error("DATA_UTIL_FST_READ_TXT",
	       "First string in file must be 'frameset'");
  }
  if (strcmp(ts2,"text")!=0){
    exit_error("DATA_UTIL_FST_READ_TXT",
	       "Second string in file must be 'text'");
  }
  ns = fscanf(fin,"%d %d %d %d",&tn,&dzn,&xn,&yn);

  if (dzn != 1)
    exit_error("DATA_UTIL_FST_READ_TXT","Depth of image must be 1'");

  //printf("    Image size:  %d x %d\n",xn,yn);
  //printf("    Number of images:  %d\n",tn);

  *rbinoc = 0;  // Assume monocular for now
  *rfcode = 4;  // Text input is read as floats
  *rccode = 0;  // Assume no color for now

  sprintf(ggstr,"    %d x %d x %d (binoc %d, fcode %d, color %d, range %d)\n",
	  xn,yn,tn,*rbinoc,*rfcode,*rccode,drange);
  mylog(logfile,ggstr);


  data = f3tensor(1,xn,1,yn,1,tn);

  if (yrev == 0){
    for(k=0;k<tn;k++){          // For each image
      for(j=(yn-1);j>=0;j--){   //     For each row (starting with top row),
	for(i=0;i<xn;i++){      //       read from left to right.
	  ns = fscanf(fin,"%f",&(data[i+1][j+1][k+1]));
	}
      }
    }
  }else{
    for(k=0;k<tn;k++){          // For each image
      for(j=0;j<yn;j++){        //   For each row (starting with bottom row),
	for(i=0;i<xn;i++){      //       read from left to right.
	  ns = fscanf(fin,"%f",&(data[i+1][j+1][k+1]));
	}
      }
    }
  }
  fclose(fin);

  *rxn = xn;
  *ryn = yn;
  *rtn = tn;

  return data;
}
/**************************************-**************************************/
/*                                                                           */
/*                               WRITE_2D_DATA                               */
/*                                                                           */
/*   tcode 1 - INT                                                           */
/*   tcode 2 - FLOAT                                                         */
/*                                                                           */
/*****************************************************************************/
void write_2d_data(outfile,data,s1,s2,n1,n2,nbytes,tcode,transpose,pflag)
     char outfile[];
     void **data;
     int s1,s2,n1,n2,nbytes,tcode,transpose,pflag;
{
  int i,j;
  FILE *fopen(),*fout;
  int temp;

  if (pflag) printf("  WRITE_2D_DATA\n");

  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** File %s\n",outfile);
    exit_error("WRITE_2D_DATA","Cannot open file");
  }
  if (pflag) printf("    Writing file %s  (%d x %d)\n",outfile,n1,n2);

  temp = 16909060; // This is hex 01 02 03 04, for checking byte order
  fwrite((char *)&temp,sizeof(int),1,fout);
  if (transpose==0){
    fwrite((char *)&n1,sizeof(int),1,fout);
    fwrite((char *)&n2,sizeof(int),1,fout);
  }else{
    fwrite((char *)&n2,sizeof(int),1,fout);
    fwrite((char *)&n1,sizeof(int),1,fout);
  }
  fwrite((char *)&nbytes,sizeof(int),1,fout);
  fwrite((char *)&tcode,sizeof(int),1,fout);
  if (transpose==0){
    for(i=s1;i<(s1+n1);i++)
      fwrite((char *)data[i]+s2*nbytes,nbytes,n2,fout);
  }else{
    for(i=s2;i<(s2+n2);i++)
      for(j=s1;j<(s1+n1);j++)
	fwrite((char *)data[j]+i*nbytes,nbytes,1,fout);
  }

  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                                READ_2D_DATA                               */
/*                                                                           */
/*****************************************************************************/
void read_2d_data(infile,rdata,rn1,rn2,rnbytes,rtcode)
     char infile[];
     void ***rdata;
     int *rn1,*rn2,*rnbytes,*rtcode;
{
  int i;
  FILE *fopen(),*fin;
  int temp,revflag,nr;
  void **data;

  printf("  READ_2D_DATA\n");

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** File %s\n",infile);
    exit_error("READ_2D_DATA","Cannot open file");
  }

  nr = fread((char *)&temp,sizeof(int),1,fin);
  if (temp == 16909060)
    revflag = 0;
  else
    revflag = 1;

  revflag_fread((char *)rn1,sizeof(int),1,fin,revflag);
  revflag_fread((char *)rn2,sizeof(int),1,fin,revflag);
  revflag_fread((char *)rnbytes,sizeof(int),1,fin,revflag);
  revflag_fread((char *)rtcode,sizeof(int),1,fin,revflag);

  data = (void **)myalloc(*rn1 * sizeof(void *));
  for(i=0;i<*rn1;i++){
    data[i] = (void *)myalloc(*rn2 * *rnbytes);
    revflag_fread((char *)data[i],*rnbytes,*rn2,fin,revflag);
  }
  fclose(fin);
  *rdata = data;
}
/**************************************-**************************************/
/*                                                                           */
/*                               WRITE_3D_DATA                               */
/*                                                                           */
/*  NOTES:                                                                   */
/*  1. 'zn' is thought of as the number of frames of                         */
/*  2.  Set tcode = 1 for integer 4 byte data, =2 for float 4 byte data.     */
/*                                                                           */
/*  'xytflag' is set to 1 to cause all of the first frame to be written      */
/*  before starting the next frame.                                          */
/*                                                                           */
/*****************************************************************************/
void write_3d_data(outfile,data,n1,n2,n3,nbytes,tcode,xytflag)
     char outfile[];
     void ***data;
     int n1,n2,n3,nbytes,tcode,xytflag;
{
  int i,j,k;
  FILE *fopen(),*fout;
  int temp,in;
  void *t;

  printf("  WRITE_3D_DATA\n");

  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** File %s\n",outfile);
    exit_error("WRITE_3D_DATA","Cannot open file");
  }

  temp = 16909060; // This is hex 01 02 03 04, for checking byte order
  fwrite((char *)&temp,sizeof(int),1,fout);

  fwrite((char *)&n1,sizeof(int),1,fout);
  fwrite((char *)&n2,sizeof(int),1,fout);
  fwrite((char *)&n3,sizeof(int),1,fout);
  fwrite((char *)&nbytes,sizeof(int),1,fout);
  fwrite((char *)&tcode,sizeof(int),1,fout);

  if (xytflag == 0){
    printf("    nbytes = %d,  n1,n2,n3 = %d,%d,%d\n",nbytes,n1,n2,n3);
    for(i=0;i<n1;i++){
      for(j=0;j<n2;j++){
	fwrite((char *)(data[i][j]),nbytes,n3,fout);
      }
    }
  }else{
    for(i=0;i<n3;i++){
      in = i*nbytes;
      for(j=0;j<n1;j++)
	for(k=0;k<n2;k++){
	  t = (void *)((long)data[j][k] + in);
	  fwrite((char *)t,nbytes,1,fout);
	}
    }
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                             WRITE_3D_DATA_PART                            */
/*                                                                           */
/*  NOTES:                                                                   */
/*  1. 'zn' is thought of as the number of frames of                         */
/*  2.  Set tcode = 1 for integer 4 byte data, =2 for float 4 byte data.     */
/*                                                                           */
/*  TCODE                                                                    */
/*         1 - integer 4-byte                                                */
/*         2 - float 4-byte                                                  */
/*        11 - integer packed color data, bytes are:  0RGB                   */
/*                                                                           */
/*****************************************************************************/
void write_3d_data_part(outfile,data,i1,n1,i2,n2,i3,n3,nbytes,tcode,xytflag)
     char outfile[];
     void ***data;
     int i1,n1,i2,n2,i3,n3,nbytes,tcode;
     int xytflag;   // 1-write x,y,z as t,x,y, thus to be shown as movie in t
{
  int i,j,k;
  FILE *fopen(),*fout;
  int temp,in;
  void *t;

  printf("  WRITE_3D_DATA_PART\n");

  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** File %s\n",outfile);
    exit_error("WRITE_3D_DATA_PART","Cannot open file");
  }

  temp = 16909060; // This is hex 01020304, for checking byte order
  fwrite((char *)&temp,sizeof(int),1,fout);

  fwrite((char *)&n1,sizeof(int),1,fout);
  fwrite((char *)&n2,sizeof(int),1,fout);
  fwrite((char *)&n3,sizeof(int),1,fout);
  fwrite((char *)&nbytes,sizeof(int),1,fout);
  fwrite((char *)&tcode,sizeof(int),1,fout);

  if (xytflag == 0){
    printf("nbytes = %d,  n1,n2,n3 = %d,%d,%d\n",nbytes,n1,n2,n3);
    printf("i1,i2,i3 = %d,%d,%d\n",i1,i2,i3);
    for(i=i1;i<(i1+n1);i++){
      for(j=i2;j<(i2+n2);j++){
	t = (void *)((long)data[i][j] + i3*nbytes);
	//t = (void *)((long)&(data[i][j][0]) + i3*nbytes);
	fwrite((char *)t,nbytes,n3,fout);
	/*fwrite((char *)(&((char *)data[i][j][i3])),nbytes,n3,fout);*/
      }
    }
  }else{
    for(i=i3;i<(i3+n3);i++){
      in = i*nbytes;
      for(j=i1;j<(i1+n1);j++)
	for(k=i2;k<(i2+n2);k++){
	  t = (void *)((long)data[j][k] + in);
	  fwrite((char *)t,nbytes,1,fout);
	}
    }
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                              READ_3D_DATA_OLD                             */
/*                                                                           */
/*****************************************************************************/
void read_3d_data_old(infile,rdata,rn1,rn2,rn3,rnbytes,rtcode,xytflag)
     char infile[];
     void ****rdata;
     int *rn1,*rn2,*rn3,*rnbytes,*rtcode,xytflag;
{
  int i,j,k;
  FILE *fopen(),*fin;
  int temp,revflag,nr;
  void ***data,*t;
  /*short sd[1000];*/

  printf("  READ_3D_DATA_OLD\n");
  printf("  **** WYETH - this is old routine, conventions changed in new\n");
  printf("  **** WYETH - this is old routine, conventions changed in new\n");

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** File %s\n",infile);
    exit_error("READ_3D_DATA","Cannot open file");
  }

  nr = fread((char *)&temp,sizeof(int),1,fin);
  if (temp == 16909060)
    revflag = 0;
  else
    revflag = 1;

  revflag_fread((char *)rn1,sizeof(int),1,fin,revflag);
  revflag_fread((char *)rn2,sizeof(int),1,fin,revflag);
  revflag_fread((char *)rn3,sizeof(int),1,fin,revflag);
  revflag_fread((char *)rnbytes,sizeof(int),1,fin,revflag);
  revflag_fread((char *)rtcode,sizeof(int),1,fin,revflag);

  /*** WYETH - note, the 'revflag' param was missing in the two calls
    below - could this have affected the design of the code? ***/

  printf("    %d x %d x %d\n",*rn1,*rn2,*rn3);

  if (xytflag == 0){
    data = (void ***)myalloc(*rn3 * sizeof(void **));
    for(i=0;i<*rn3;i++){
      data[i] = (void **)myalloc(*rn1 * sizeof(void *));
      for(j=0;j<*rn1;j++){
	data[i][j] = (void *)myalloc(*rn2 * *rnbytes);
	revflag_fread((char *)data[i][j],*rnbytes,*rn2,fin,revflag);
      }
    }
  }else{
    data = (void ***)myalloc(*rn1 * sizeof(void **));
    for(i=0;i<*rn1;i++){
      data[i] = (void **)myalloc(*rn2 * sizeof(void *));
      for(j=0;j<*rn2;j++)
	data[i][j] = (void *)myalloc(*rn3 * *rnbytes);
    }
    
    for(i=0;i<*rn3;i++)
      for(j=0;j<*rn1;j++)
	for(k=0;k<*rn2;k++){
	  t = (void *)((long)data[j][k] + i* *rnbytes);
	  revflag_fread((char *)t,*rnbytes,1,fin,revflag);
	}
  }
  fclose(fin);
  *rdata = data;
}
/**************************************-**************************************/
/*                                                                           */
/*                               READ_3D_DATA                                */
/*                                                                           */
/*  WYETH -the old one looks wrong (Oct 9, 2004)                             */
/*                                                                           */
/*****************************************************************************/
void read_3d_data(infile,rdata,rn1,rn2,rn3,rnbytes,rtcode,xytflag)
     char infile[];
     void ****rdata;      // Pointer to 3D array to create and return
     int *rn1,*rn2,*rn3;  // Dimensions of data
     int *rnbytes;        // number of bytes per entry
     int *rtcode;         // 1-int, 2-float
     int xytflag;         //
{
  int i,j,k;
  FILE *fopen(),*fin;
  int temp,revflag,nr;
  void ***data,*t;
  /*short sd[1000];*/

  printf("  READ_3D_DATA\n");

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** File %s\n",infile);
    exit_error("READ_3D_DATA","Cannot open file");
  }

  nr = fread((char *)&temp,sizeof(int),1,fin);
  if (temp == 16909060)
    revflag = 0;
  else
    revflag = 1;

  revflag_fread((char *)rn1,sizeof(int),1,fin,revflag);
  revflag_fread((char *)rn2,sizeof(int),1,fin,revflag);
  revflag_fread((char *)rn3,sizeof(int),1,fin,revflag);
  revflag_fread((char *)rnbytes,sizeof(int),1,fin,revflag);
  revflag_fread((char *)rtcode,sizeof(int),1,fin,revflag);

  data = (void ***)myalloc(*rn1 * sizeof(void **));
  for(i=0;i<*rn1;i++){
    data[i] = (void **)myalloc(*rn2 * sizeof(void *));
    for(j=0;j<*rn2;j++)
      data[i][j] = (void *)myalloc(*rn3 * *rnbytes);
  }
  
  if (xytflag == 0){
    for(i=0;i<*rn1;i++){
      for(j=0;j<*rn2;j++){
	revflag_fread((char *)data[i][j],*rnbytes,*rn3,fin,revflag);
      }
    }
  }else{
    for(i=0;i<*rn3;i++)
      for(j=0;j<*rn1;j++)
	for(k=0;k<*rn2;k++){
	  t = (void *)((long)data[j][k] + i* *rnbytes);
	  revflag_fread((char *)t,*rnbytes,1,fin,revflag);
	}
  }
  fclose(fin);
  *rdata = data;
}
/**************************************-**************************************/
/*                                                                           */
/*                          COUNT_COLUMNS_COLUMN_DATA                        */
/*                                                                           */
/*****************************************************************************/
int count_columns_column_data(fin)
     FILE *fin;
{
  char c,oldc;
  int cols,done,ic;

  cols = 0;
  oldc = ' ';
  done = 0;
  if ((ic = getc(fin))==EOF) // Must assign to int for EOF check.
    done = 1;
  c = (char)ic;
  while ((!done)&&(c!='\n')){
    if (((c!=' ')&&(c!='\t'))&&((oldc==' ')||(oldc=='\t')))
      cols += 1;
    oldc = c;
    if ((ic=getc(fin))==EOF)
      done = 1;
    c = (char)ic;
  }
  if (done && (cols==0))
    return -1;
  else
    return cols;
}
/**************************************-**************************************/
/*                                                                           */
/*                              READ_COLUMN_DATA                             */
/*                                                                           */
/*   rcdata[rnrec][rncol]                                                    */
/*                                                                           */
/*****************************************************************************/
void read_column_data(infile,pflag,rcdata,rnrec,rncol) // WYETH.2010.Aug
     char infile[];
     int pflag;            // 1-print, 0-quiet
     char ****rcdata;
     int *rnrec,*rncol;
{
  FILE *fopen(),*fin;
  int i,j;
  int cols,columns,records,flag,ns;
  char ***cdata,tstr[1024];

  if (pflag) printf("  READ_COLUMN_DATA\n");

  // read data, counting records and columns
  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** Cannot open file %s.\n",infile);
    exit_error("READ_COLUMN_DATA","Cannot open file");
  }
  if (pflag) printf("    Reading %s\n",infile);


  records = 0;
  columns = count_columns_column_data(fin);
  cols = columns;
  flag = 0;
  while ((cols!=-1)&&(flag==0)){
    if (cols>0)
      records+=1;
    if ((cols!=columns)&&(cols!=0))
      flag = records;
    cols = count_columns_column_data(fin);
  }
  fclose(fin);

  /* report results, read and store data */
  if (flag > 0){
    printf("    at record %d\n",flag);
    exit_error("READ_COLUMN_DATA","Unequal column count");
    cdata = NULL;
  }else{
    /* print results */
    if (pflag){
      printf("    %d columns found\n",columns);
      printf("    %d records found\n",records);
    }
    /* allocate storage */
    cdata = (char ***)malloc(records*sizeof(char **));
    for (i=0;i<records;i++)
      cdata[i] = (char **)malloc(columns*sizeof(char *));
    /* read all data */
    fin = fopen(infile,"r");
    for(i=0;i<records;i++)
      for(j=0;j<columns;j++){
	ns = fscanf(fin,"%s",tstr);
	cdata[i][j] = (char *)malloc((strlen(tstr)+1)*sizeof(char));
	strcpy(cdata[i][j],tstr);
      }
    fclose(fin);
  }
  *rcdata = cdata; *rncol = columns; *rnrec = records;
}
/**************************************-**************************************/
/*                                                                           */
/*                            READ_VAN_HATEREN_OLD                           */
/*                                                                           */
/*  Read the van hateren image with the given index 'k'.                     */
/*                                                                           */
/*****************************************************************************/
int **read_vanhateren_old(imagedir,k,x0,y0,xn,yn)
     char imagedir[];    // directory name
     int k;              // image index
     int x0,y0;          // start coords
     int xn,yn;          // size of patch
{
  int i,j;
  int revflag,w,h,**data,**patch,offset,offset_line;
  short ts;
  FILE *fopen(),*fin;
  char fullname[SLEN],*numstr;

  // WYETH - need a way to check for endian.  We could read the first
  // value of image 0 (or 1) and compare to a known value.
  revflag = 1;

  w = 1536;  // Image width (pix)
  h = 1024;  // Image height (pix)

  numstr = get_leading_zero_string_from_int(k,5);
  sprintf(fullname,"%s/imk%s.iml",imagedir,numstr);
  myfree(numstr);

  //printf("READ_VAN_HATEREN\n");
  //printf("  %s\n",fullname);

  if ((fin = fopen(fullname,"r")) == NULL){
    printf("  *** File %s\n",fullname);
    exit_error("READ_VAN_HATEREN","Cannot open file");
  }
  data = get_2d_iarray(w,h);

  for(i=0;i<h;i++){
    for(j=0;j<w;j++){
      revflag_fread(&ts,2,1,fin,revflag);
      data[j][h-i-1] = (int)ts;  // First line read is top of image
    }
  }
  fclose(fin);

  //
  // data[0][0] holds lower left corner, and [x][y] is order
  //

  write_2d_data("zz.vh.22.2d",(void **)data,0,0,w,h,4,1,1,1);  // 1-transpose
  exit(0);

  patch = get_2d_iarray(xn,yn);
  for(i=0;i<xn;i++)
    for(j=0;j<yn;j++)
      patch[i][j] = data[x0+i][y0+j];

  free_2d_farray(data,h);

  //write_2d_data("zz.patch.2d",patch,0,0,xn,yn,4,1,1,1);  // 1-transpose
  //exit(0);
  
  return patch;
}
/**************************************-**************************************/
/*                                                                           */
/*                              READ_VAN_HATEREN                             */
/*                                                                           */
/*  Read the van hateren image with the given index 'k'.                     */
/*                                                                           */
/*****************************************************************************/
int **read_vanhateren(imagedir,k,x0,y0,xn,yn)
     char imagedir[];    // directory name
     int k;              // image index
     int x0,y0;          // start coords
     int xn,yn;          // size of patch
{
  int i,j;
  int revflag,w,h,**patch,offset,offset_line,yi;
  short ts,*tsa;
  FILE *fopen(),*fin;
  char fullname[SLEN],*numstr;

  // WYETH - need a way to check for endian.  We could read the first
  // value of image 0 (or 1) and compare to a known value.
  revflag = 1;

  w = 1536;  // Image width (pix)
  h = 1024;  // Image height (pix)

  numstr = get_leading_zero_string_from_int(k,5);
  sprintf(fullname,"%s/imk%s.iml",imagedir,numstr);
  myfree(numstr);

  //printf("READ_VAN_HATEREN\n");
  //printf("  %s\n",fullname);

  if ((fin = fopen(fullname,"r")) == NULL){
    printf("  *** File %s\n",fullname);
    exit_error("READ_VAN_HATEREN","Cannot open file");
  }

  offset = ((h - (y0+yn)) * w + x0) * 2;  // Offset from start of file
  offset_line = (w - xn) * 2;   // Offset to start of patch in next line

  /*
  printf("x0,y0 = %d %d\n",x0,y0);
  printf("xn,yn = %d %d\n",xn,yn);
  printf("offset = %d (bytes)\n",offset);
  printf("offset_line = %d (bytes)\n",offset_line);
  */

  tsa = (short *)myalloc(xn*sizeof(short));
  patch = get_2d_iarray(xn,yn);

  fseek(fin,(long)offset,SEEK_SET);  // Seek rel. to begin
  for(i=0;i<yn;i++){
    revflag_fread(tsa,2,xn,fin,revflag);
    yi = yn-i-1;
    for(j=0;j<xn;j++){
      patch[j][yi] = (int)tsa[j];
    }
    fseek(fin,(long)offset_line,SEEK_CUR);  // Seek rel. to Current position
  }
  fclose(fin);

  myfree(tsa);
  //
  // patch[0][0] holds lower left corner, and [x-horiz][y-vert] applies
  //

  //write_2d_data("zz.patch.2d",patch,0,0,xn,yn,4,1,1,1);  // 1-transpose
  //exit(0);
  
  return patch;
}
/**************************************-**************************************/
/*                                                                           */
/*                            READ_PPM_IMAGE_TYPE                            */
/*                                                                           */
/*****************************************************************************/
int read_ppm_image_type(infile)
     char infile[];
{
  int k,ns;
  FILE *fopen(),*fin;
  char temp[SLEN];

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** File %s\n",infile);
    exit_error("READ_PPM_IMAGE_TYPE","Cannot open file");
  }

  ns = fscanf(fin,"%s",temp);
  k = (int)temp[1] - (int)'0';
  fclose(fin);

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                              READ_PPM_HEADER                              */
/*                                                                           */
/*  For more information:  man ppm                                           */
/*                                                                           */
/*  type 1 - bitmap (PBM) ascii                                              */
/*  type 2 - grayscale (PGM) ascii                                           */
/*  type 3 - color (PPM) ascii                                               */
/*  type 4 -                                                                 */
/*  type 5 - binary form of P2                                               */
/*  type 6 - byte format, 1 byte per color component                         */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/
void read_ppm_header(infile,rt,rw,rh,rm)
     char infile[];
     int *rt,*rw,*rh,*rm; /* type, width, height, max value */
{
  FILE *fopen(),*fin;
  char temp[SLEN],*tc;
  int t,w,h,m,ns;

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** File %s\n",infile);
    exit_error("READ_PPM_IMAGE_TYPE","Cannot open file");
  }

  t = w = h = m = -1;
  while (m == -1){
    ns = fscanf(fin,"%s",temp);
    if (temp[0] == '#'){
      tc = fgets(temp,SLEN,fin); // read rest of line
    }else{
      if (t == -1)
	t = (int)temp[1] - (int)'0';
      else if (w == -1)
	w = atoi(temp);
      else if (h == -1)
	h = atoi(temp);
      else if (m == -1)
	m = atoi(temp);
    }
  }
  fclose(fin);

  *rt = t; *rw = w; *rh = h; *rm = m;
}
/**************************************-**************************************/
/*                                                                           */
/*                               READ_PGM_DATA                               */
/*                                                                           */
/*  ASCII gray scale file.                                                   */
/*                                                                           */
/*  Format:                                                                  */
/*    P2                                                                     */
/*    # comment                                                              */
/*    w h m # comment                                                        */
/*                                                                           */
/*****************************************************************************/
void read_pgm_data(infile,rt,rw,rh,rm,rdata)
     char infile[];
     int *rt,*rw,*rh,*rm; // type, width, height, max value
     int ***rdata;
{
  int i,j;
  FILE *fopen(),*fin;
  char temp[SLEN],*tc;
  int t,w,h,m,ns;
  int **data;

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** File %s\n",infile);
    exit_error("READ_PGM_DATA","Cannot open file");
  }

  t = w = h = m = -1;
  while (m == -1){
    ns = fscanf(fin,"%s",temp);
    if (temp[0] == '#'){
      tc = fgets(temp,SLEN,fin); // read rest of line
    }else{
      if (t == -1)
	t = (int)temp[1] - (int)'0';
      else if (w == -1)
	w = atoi(temp);
      else if (h == -1)
	h = atoi(temp);
      else if (m == -1)
	m = atoi(temp);
    }
  }

  if (t==2){
    data = (int **)myalloc(h*sizeof(int *));
    for(i=0;i<h;i++){
      data[i] = (int *)myalloc(w*sizeof(int));
      for(j=0;j<w;j++)
	ns = fscanf(fin,"%d",&(data[i][j]));
    }
  }else{
    printf("  Type = %d\n",t);
    exit_error("READ_PGM_DATA","PPM type is not 2");
  }
  
  *rt = t; *rw = w; *rh = h; *rm = m;
  *rdata = data;
}
/**************************************-**************************************/
/*                                                                           */
/*                               READ_PGM_P5_DATA                            */
/*                                                                           */
/*  BINARY data, 1 or 2 bytes, depending on 'm'                              */
/*                                                                           */
/*  Format:                                                                  */
/*    P5                                                                     */
/*    # comment                                                              */
/*    w h m # comment                                                        */
/*                                                                           */
/*  *** Returned data is [h][w]                                              */
/*                                                                           */
/*****************************************************************************/
void read_pgm_p5_data(infile,rt,rw,rh,rm,rdata)
     char infile[];
     int *rt,*rw,*rh,*rm; // type, width, height, max value
     int ***rdata;
{
  int i,j;
  FILE *fopen(),*fin;
  char temp[SLEN],*tc;
  int t,w,h,m,ns;
  int **data,revflag;
  unsigned short ts;
  char c;

  // WYETH I needed to set revflag to 1 for Bartlett's PGM images:
  revflag = 1;
  printf("  *** READ_PGM_P5_DATA:  revflag is set to 1 for 2byte endian\n");

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** File %s\n",infile);
    exit_error("READ_PGM_P5_DATA","Cannot open file");
  }

  t = w = h = m = -1;
  while (m == -1){
    ns = fscanf(fin,"%s",temp);
    if (temp[0] == '#'){
      tc = fgets(temp,SLEN,fin); // read rest of line
    }else{
      if (t == -1)
	t = (int)temp[1] - (int)'0';
      else if (w == -1)
	w = atoi(temp);
      else if (h == -1)
	h = atoi(temp);
      else if (m == -1)
	m = atoi(temp);
    }
  }

  ns = fscanf(fin,"%c",&c); // Remove the last CR
  if ((int)c != 10)
    exit_error("READ_PGM_P5_DATA","Expecting ASCII 10 (newline)");

  if (t!=5)
    exit_error("READ_PGM_P5_DATA","PPM type is not 5");

  data = (int **)myalloc(h*sizeof(int *));
  for(i=0;i<h;i++){
    data[i] = (int *)myalloc(w*sizeof(int));
    for(j=0;j<w;j++){
      if (m <= 255)
	ns = fread(&(data[i][j]),1,1,fin);
      else{
	//fread(&ts,2,1,fin);
	revflag_fread(&ts,2,1,fin,revflag);
	data[i][j] = (int)ts;
	if ((i == 177) && (j > 100) && (j < 120))
	  printf("  ts = %d\n",(int)ts);
      }
    }
  }
  
  *rt = t; *rw = w; *rh = h; *rm = m;
  *rdata = data;
}
/**************************************-**************************************/
/*                                                                           */
/*                            WRITE_PGM_DATA_FLOAT                           */
/*                                                                           */
/*  ASCII gray scale file.                                                   */
/*                                                                           */
/*  't' should be 2.                                                         */
/*                                                                           */
/*****************************************************************************/
void write_pgm_data_float(outfile,t,w,h,m,w0,h0,data)
     char outfile[];
     int t,w,h,m; /* type, width, height, max value */
     int w0,h0;   // WYETH ADDED THIS LINE Jun 2010, had been omitted
                  // just seen at random, was not testing this code ***
     float **data;
{
  FILE *fopen(),*fout;
  int i,j,k;
  int x,nhi,nlo;

  if (t!=2)
    exit_error("WRITE_PGM_DATA_FLOAT","t is not 2");

  if ((fout = fopen(outfile,"w")) == NULL){
    printf("  *** File %s\n",outfile);
    exit_error("WRITE_PGM_DATA_FLOAT","Cannot open file");
  }

  fprintf(fout,"P%d\n",t);
  fprintf(fout,"# CREATOR:  WRITE_PGM_DATA_FLOAT by Wyeth Bair\n");
  fprintf(fout,"%d %d %d\n",w,h,m);

  k = 0;
  nlo = nhi = 0;
  for(j=0;j<h;j++){
    for(i=0;i<w;i++){
      x = (int)(0.5 + data[h0+j][w0+i]); /* Round to integer */
      if (x < 0){
	x = 0;
	nlo += 1;
      }else if (x > m){
	x = m;
	nhi += 1;
      }
      k += 1;
      if (m < 1000){
	fprintf(fout,"%3d ",x);
	if (k == 17){
	  fprintf(fout,"\n");
	  k = 0;
	}
      }else if (m < 10000){
	fprintf(fout,"%4d ",x);
	if (k == 13){
	  fprintf(fout,"\n");
	  k = 0;
	}
      }else
	exit_error("WRITE_PGM_DATA_FLOAT","m too large"); /*** hack ***/
    }
  }
  fclose(fout);

  if (nlo > 0)
    printf("   *** WARNING:  %d values < 0 set to 0.\n",nlo);
  if (nhi > 0)
    printf("   *** WARNING:  %d values > %d set to %d.\n",nhi,m,m);
}
/**************************************-**************************************/
/*                                                                           */
/*                              WRITE_PPM_6_DATA                             */
/*                                                                           */
/*  Write 24 bit color data.                                                 */
/*                                                                           */
/*****************************************************************************/
void write_ppm_6_data(outfile,data,w,h)
     char outfile[];
     unsigned char *data;
     int w,h;
{
  FILE *fout,*fopen();

//  printf("DATA[0] = %d , %d , ...\n",(int)data[0],(int)data[1]);
  
  fout = fopen(outfile,"w");
  fprintf(fout,"P6\n%d %d\n255\n",w,h);
  fwrite(data,sizeof(char),w*h*3,fout);

  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                             WRITE_PPM_6_RGB_DATA                          */
/*                                                                           */
/*  Assume that float data is in range from [0..1].                          */
/*                                                                           */
/*****************************************************************************/
void write_ppm_6_rgb_data(outfile,w,h,r,g,b,yflag)
     char outfile[];
     int w,h;
     float **r,**g,**b;     // [w][h] image
     int yflag;             // 1-reverse y-axis, 0-do not reverse
{
  int i,j,l;
  int yi;
  unsigned char *u;

  u = (unsigned char *)myalloc(w*h*3*sizeof(unsigned char));
  l = 0;
  for(i=0;i<h;i++){

    if (yflag == 0)
      yi = i;
    else
      yi = h-1 - i;  // Reverse the y-axis order

    for(j=0;j<w;j++){
      u[l]   = (unsigned char)my_rint(255.0 * r[j][yi]);
      u[l+1] = (unsigned char)my_rint(255.0 * g[j][yi]);
      u[l+2] = (unsigned char)my_rint(255.0 * b[j][yi]);
      l += 3;
    }
  }
  write_ppm_6_data(outfile,u,w,h);

  myfree(u);
}
/**************************************-**************************************/
/*                                                                           */
/*                          WRITE_PPM_6_3D_GRAY_FRAMES                       */
/*                                                                           */
/*  Assume that float data is [1..xn][1..yn][1..tn].                         */
/*                                                                           */
/*  Write filters every 'dt' steps.                                          */
/*                                                                           */
/*****************************************************************************/
void write_ppm_6_3d_gray_frames(outfile,data,xn,yn,tn,dt)
     char outfile[];
     float ***data;
     int xn,yn,tn;
     int dt;            // NOTE - 'dt' should be 1, unless skipping frames
{
  int i,j,k,l;
  int i0,rangeflag;
  unsigned char *u,uu;
  float dmin,dmax;
  char name[SLEN];

  //
  //  WYETH 2018Apr20 - fixed this to check range, and to be able to handle
  //                    ranges from 0-1 and 0-255.
  //

  get_min_max_3d_farray(data,1,xn,1,yn,1,tn,&dmin,&dmax);
  printf("    Stimulus values range from %f to %f\n",dmin,dmax);
  if (dmin < 0.0)
    exit_error("WRITE_PPM_6_3D_GRAY_FRAMES","Min value is negative");
  if (dmax >= 255.5)
    exit_error("WRITE_PPM_6_3D_GRAY_FRAMES","Max value is > 255");

  if (dmax > 1.0001)
    rangeflag = 255;
  else
    rangeflag = 1;

  u = (unsigned char *)myalloc(xn*yn*3*sizeof(unsigned char));
  for(i=1;i<=tn;i+=dt){
    i0 = i-1;
    if (i0 < 10)
      sprintf(name,"%s.000%d.ppm",outfile,i0);
    else if (i0 < 100)
      sprintf(name,"%s.00%d.ppm",outfile,i0);
    else if (i0 < 1000)
      sprintf(name,"%s.0%d.ppm",outfile,i0);
    else
      sprintf(name,"%s.%d.ppm",outfile,i0);
    l = 0;

    for(k=yn;k>=1;k--){
      for(j=1;j<=xn;j++){
	if (rangeflag == 1)
	  uu = (unsigned char)my_rint(255.0 * data[j][k][i]);
	else
	  uu = (unsigned char)my_rint(data[j][k][i]);
	u[l  ] = uu; // R, G, B
	u[l+1] = uu;
	u[l+2] = uu;
	l += 3;
      }
    }
    write_ppm_6_data(name,u,xn,yn);
  }
  myfree(u);
}
/**************************************-**************************************/
/*                                                                           */
/*                           WYTIFF_GET_EMPTY_STRUCT                         */
/*                                                                           */
/*****************************************************************************/
struct wytif_struct *wytiff_get_empty_struct()
{
  struct wytif_struct *w;

  w = (struct wytif_struct *)myalloc(sizeof(struct wytif_struct));

  w->w = 0;
  w->h = 0;
  w->bitpsamp = -1;
  w->compression = -1;
  w->photo_interp = -1;
  w->orientation = -1;
  w->sppix = -1;
  w->planarconfig = -1;
  w->descrip = (char *)NULL;
  w->nstrip = 0;
  w->soff = NULL;
  w->row_p_strip = 0;
  w->sbyte = NULL;
  w->xres = 0.0;
  w->yres = 0.0;
  w->resunit = 0;
  w->software = (char *)NULL;
  w->datetime = (char *)NULL;
  w->predictor = -1;
  w->newsubfiletype = -1;

  w->ifd_next = 0;
  w->next = NULL;

  return w;
}
/**************************************-**************************************/
/*                                                                           */
/*                                WYTIFF_FREE                                */
/*                                                                           */
/*****************************************************************************/
void wytiff_free(w,data_flag)
     struct wytif_struct *w;
     int data_flag;            // 1-free image data, 0-do not free image data
{
  int i;

  if (w->descrip != NULL)
    myfree(w->descrip);

  if (w->soff != NULL)
    myfree(w->soff);

  if (w->sbyte != NULL)
    myfree(w->sbyte);

  if (w->software != NULL)
    myfree(w->software);

  if (w->datetime != NULL)
    myfree(w->datetime);

  if (data_flag == 1){
    for(i=0;i<w->w;i++)
      myfree(w->data[i]);
    myfree(w->data);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            WYTIFF_PRINT_STRUCT                            */
/*                                                                           */
/*****************************************************************************/
void wytiff_print_struct(w)
     struct wytif_struct *w;
{
  int i;

  printf("  WYTIFF_PRINT_STRUCT\n");

  printf("    %d x %d  (%d bits/pix)\n",w->w,w->h,w->bitpsamp);
  if (w->compression == 1)
    printf("    No compression.\n");
  else
    printf("    Compression:  %d\n",w->compression);
  printf("    Photo Interpretation:  %d\n",w->photo_interp);
  if (w->descrip != NULL)
    printf("    Description:  ==>%s<==\n",w->descrip);
  printf("    Strips:  %d\n",w->nstrip);
  printf("    Rows/strip:  %d\n",w->row_p_strip);
  printf("    Resolution unit:  %d\n",w->resunit);
  printf("    Orientation flag:  %d\n",w->orientation);
  printf("    Planar Configuration:  %d\n",w->planarconfig);
  printf("    Samp per pixel:  %d\n",w->sppix);
  printf("    Next IFD Offset:  %d\n",(int)w->ifd_next);

  if (w->software != NULL)
    printf("    Software:  %s\n",w->software);
  if (w->datetime != NULL)
    printf("    Date Time:  %s\n",w->datetime);
  if (w->predictor != -1)
    printf("    Predictor: %d\n",w->predictor);
  if (w->newsubfiletype != -1)
    printf("    NewSubFileType  %d\n",w->newsubfiletype);

  printf("    Strip Offsets\n");
  if (w->nstrip > 0){
    for(i=0;i<w->nstrip;i++)
      printf("      %d\n",(int)(w->soff[i]));
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                             WYTIFF_READ_STRING                            */
/*                                                                           */
/*****************************************************************************/
char *wytiff_read_string(fin,offset,cnt)
     FILE *fin;
     unsigned int offset;
     int cnt;
{
  int nr;
  char *tstr;

  if (cnt > 0){
    tstr = (char *)myalloc(cnt*sizeof(char));

    fseek(fin,(long)offset,SEEK_SET);  // Seek rel. to begin
    nr = fread(tstr,1,cnt,fin);
    if (tstr[cnt-1] != '\0'){
      //printf("  *** READ_TIFF_STRING:  adding NULL to end of string.\n");
      tstr[cnt-1] = '\0';
    }
    
  }else{
    tstr = NULL;
  }

  return tstr;
}
/**************************************-**************************************/
/*                                                                           */
/*                            WYTIFF_READ_IFD_ENTRY                           */
/*                                                                           */
/*  Read a TIF "image file directory" (IFD) entry, which provides info in a  */
/*  TIF file.                                                                */
/*                                                                           */
/*****************************************************************************/
void wytiff_read_ifd_entry(fin,revflag,rtid,rfty,rcnt,roff)
     FILE *fin;
     int revflag;
     unsigned short *rtid,*rfty;
     unsigned int *rcnt,*roff;
{
  unsigned int cnt,offset;
  unsigned short tagid,fieldtype;


  // Bytes 0-1:  Tag that IDs the field  (IFDs are sorted by tag number)
  revflag_fread(&tagid,2,1,fin,revflag);

  // Bytes 2-3:  Field type
  revflag_fread(&fieldtype,2,1,fin,revflag);
  // Field types and sizes:
  // 1-Byte
  // 2-ASCII, 8-bits, last byte may be NULL
  // 3-SHORT, 2-byte unsigned int
  // 4-LONG,  4-byte unsigned int
  // 5-RATIONAL,  Two LONGs, 1st is numerator, 2nd denominator
  //
  //   new in TIFF 6.0
  //
  // 6-SBYTE  8-bit signed (2s-complement) int
  // 7-UNDEFINED  byte may contain anything
  // 8-SSHORT 2-byte signed (2s-complement) int
  // 9-SLONG  4-byte signed (2s-complement) int
  // 10-SRATIONAL  two SLONGs, representing fraction
  // 11-FLOAT single precision (4-byte) IEEE format
  // 12-DOUBLE double precision (8-byte) IEEE format
  //
  //  Readers should skip any unknown field types
  //

  // Bytes 4-7:  The number of values of the indicated type
  revflag_fread(&cnt,4,1,fin,revflag);

  // Bytes 8-11:  The Value Offset (an offset w/i the file)
  //              May contain the value, if it fits in 4 bytes
  revflag_fread(&offset,4,1,fin,revflag);

  *rtid = tagid;
  *rfty = fieldtype;
  *rcnt = cnt;
  *roff = offset;
}
/**************************************-**************************************/
/*                                                                           */
/*                              WYTIFF_READ_IFD                              */
/*                                                                           */
/*  Read a TIF "image file directory" (IFD).                                 */
/*                                                                           */
/*****************************************************************************/
void wytiff_read_ifd(fin,revflag,w,pflag)
     FILE *fin;
     int revflag;
     struct wytif_struct *w;
     int pflag;  // Print flag
{
  int i;
  unsigned int cnt,offset;
  unsigned short n,tagid,fldtyp;
  int cnt_desc,cnt_date,cnt_soft;
  unsigned int off_desc,off_date,off_soft;
  unsigned int off_strip;
  unsigned int off_byte;

  off_desc  = cnt_desc = 0;
  off_date  = cnt_date = 0;
  off_soft  = cnt_soft = 0;

  // Bytes 0-1:  Tag that IDs the field  (IFDs are sorted by tag number)
  revflag_fread(&n,2,1,fin,revflag);
  if (pflag) printf("    n_ifd_entries  %d\n",(int)n);

  for(i=0;i<n;i++){
    wytiff_read_ifd_entry(fin,revflag,&tagid,&fldtyp,&cnt,&offset);

    if (pflag)
      printf("  Entry %2d   TagID  %6d   FldTyp  %2d   Cnt  %4d   Offset  %d\n",
	     i,(int)tagid,(int)fldtyp,(int)cnt,(int)offset);

    if ((tagid == 254) && (cnt == 1)){
      w->newsubfiletype = (int)offset;
    }else if ((tagid == 256) && (cnt == 1)){
      w->w = (int)offset;
    }else if ((tagid == 257) && (cnt == 1)){
      w->h = (int)offset;
    }else if ((tagid == 258) && (cnt == 1)){
      w->bitpsamp = (int)offset;
    }else if ((tagid == 259) && (cnt == 1)){
      w->compression = (int)offset;
    }else if ((tagid == 262) && (cnt == 1)){
      w->photo_interp = (int)offset;
    }else if (tagid == 270){
      off_desc = offset;
      cnt_desc = (int)cnt;
      // cnt is string length, including NULL
    }else if (tagid == 273){
      off_strip = offset;
      w->nstrip = (int)cnt;
    }else if (tagid == 274){
      w->orientation = (int)offset;
    }else if (tagid == 277){
      w->sppix = (int)offset;
    }else if ((tagid == 278) && (cnt == 1)){
      w->row_p_strip = (int)offset;
    }else if (tagid == 279){
      if (w->nstrip != (int)cnt)
	exit_error("WYTIFF_READ_IFD","Strip byte counts mis-match.");
      off_byte = offset;
    }else if ((tagid == 282) && (cnt == 1)){
      ; // x-resolution
    }else if ((tagid == 283) && (cnt == 1)){
      ; // y-resolution
    }else if (tagid == 284){
      w->planarconfig = (int)offset;
    }else if ((tagid == 296) && (cnt == 1)){
      w->resunit = (int)offset;
    }else if (tagid == 305){
      off_soft = offset;
      cnt_soft = (int)cnt;
    }else if (tagid == 306){
      off_date = offset;
      cnt_date = (int)cnt;
    }else if ((tagid == 317) && (cnt == 1)){
      w->predictor = (int)offset;
    }else if (tagid > 1000){
      ; // Some high-value tag; Pairie uses, e.g., 33628
    }else{
      printf("tagid = %d\n",tagid);
      exit_error("WYTIFF_READ_IFD","Unknown condition.");
    }
  }

  //
  //  Read offset of next IFD (for multiple images)
  //
  revflag_fread(&offset,4,1,fin,revflag);
  w->ifd_next = offset;

  //
  //  Read some strings
  //
  if (cnt_desc != 0){
    w->descrip = wytiff_read_string(fin,off_desc,cnt_desc);
  }
  if (cnt_date != 0){
    w->datetime = wytiff_read_string(fin,off_date,cnt_date);
  }
  if (cnt_soft != 0){
    w->software = wytiff_read_string(fin,off_soft,cnt_soft);
  }

  //
  //  Read strip offsets
  //
  w->soff = (unsigned int *)myalloc(w->nstrip*sizeof(unsigned int));
  fseek(fin,(long)off_strip,SEEK_SET);  // Seek rel. to begin
  revflag_fread(w->soff,4,w->nstrip,fin,revflag);

}
/**************************************-**************************************/
/*                                                                           */
/*                            WYTIFF_READ_DATA_FLOAT                         */
/*                                                                           */
/*  Read data from strips and return as a 2D farray.                         */
/*                                                                           */
/*****************************************************************************/
float **wytiff_read_data_float(fin,revflag,w)
     FILE *fin;
     int revflag;
     struct wytif_struct *w;
{
  int i,j;
  int xi,yi,nlast,nstrip,nrows;
  float **fdata;
  unsigned short *tdata;

  //
  //  The data are stored in strips, but the last strip might not have the
  //  same number of rows as the other strips.
  //

  //  Assuming 'h' rows, then there are h % row_p_strip rows in last strip
  nlast = w->h % w->row_p_strip;
  if (nlast == 0)
    nlast = w->row_p_strip;
  /*
    printf("=========>  w->h = %d\n",w->h);
    printf("=========>  row_p_strip = %d\n",w->row_p_strip);
    printf("=========>  nlast = %d\n",nlast);*/


  //
  //  Read data from strips
  //
  fdata = get_2d_farray(w->w,w->h);  // [xi][yi]
  yi = w->h-1;  // Run yi backwards, assuming top of image is scanned first

  tdata = (unsigned short *)myalloc(w->w * sizeof(unsigned short));

  nstrip = w->nstrip;

  for(i=0;i<nstrip;i++){
    fseek(fin,(long)(w->soff[i]),SEEK_SET);
    //printf("done seek %d\n",i);

    if (i == (nstrip-1))  // If this is the last strip, read the remaining rows
      nrows = nlast;
    else
      nrows = w->row_p_strip;
            
    for(j=0;j<nrows;j++){
      revflag_fread(tdata,2,w->w,fin,revflag);
      for(xi=0;xi<w->w;xi++){
	fdata[xi][yi] = (float)tdata[xi];
	//printf("xi yi %d  %d    i j   %d %d\n",xi,yi,i,j);
      }
      yi -= 1;
    }
  }
  myfree(tdata);

  return fdata;
}
/**************************************-**************************************/
/*                                                                           */
/*                             WYTIFF_IMAGE_COUNT                            */
/*                                                                           */
/*****************************************************************************/
int wytiff_image_count(w)
     struct wytif_struct *w;
{
  int n,done;
  struct wytif_struct *t;

  //printf("  WYTIFF_IMAGE_COUNT\n");

  t = w;
  n = 1;
  done = 0;
  while(done == 0){
    if (t->next == NULL)
      done = 1;
    else{
      n += 1;
      t = t->next;
    }
  }

  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                            WYTIFF_CHECK_HEADER                            */
/*                                                                           */
/*  Check the consistency of header data across images.                      */
/*                                                                           */
/*****************************************************************************/
int wytiff_check_header(w,flag)
     struct wytif_struct *w;
     int flag;  // 0-default, 1-exit on error
{
  int cflag,done;
  struct wytif_struct *t;

  //printf("  WYTIFF_CHECK_HEADER\n");

  t = w;
  cflag = 1; // Consistent
  done = 0;
  while(done == 0){

    if (t->next == NULL)
      done = 1;
    else{
      t = t->next;

      if (t->w != w->w){
	printf("    Width != %d %d\n",t->w,w->w);
	cflag = 0;
      }else if (t->h != w->h){
	printf("    Height != %d %d\n",t->h,w->h);
	cflag = 0;
      }else if (t->bitpsamp != w->bitpsamp){
 	printf("    Bits/sample != %d %d\n",t->bitpsamp,w->bitpsamp);
	cflag = 0;
      }else if (t->compression != w->compression){
 	printf("    Compression != %d %d\n",t->compression,w->compression);
	cflag = 0;
      }else if (t->photo_interp != w->photo_interp){
 	printf("    Photo_interp != %d %d\n",t->photo_interp,w->photo_interp);
	cflag = 0;
      }else if (t->orientation != w->orientation){
 	printf("    Orientation != %d %d\n",t->orientation,w->orientation);
	cflag = 0;
      }else if (t->sppix != w->sppix){
 	printf("    Samp/pix != %d %d\n",t->sppix,w->sppix);
	cflag = 0;
      }else if (t->planarconfig != w->planarconfig){
 	printf("    Planar Config != %d %d\n",t->planarconfig,w->planarconfig);
	cflag = 0;
      }
      //  ... there are other checks that could be done ...
      //  ... there are other checks that could be done ...
      //  ... there are other checks that could be done ...

      
      if (cflag == 0)
	done = 0;
    }
  }

  if ((flag == 1) && (cflag == 0))
    exit_error("WYTIFF_CHECK_HEADER","Mismatch found across images");

  return cflag;
}
/**************************************-**************************************/
/*                                                                           */
/*                                WYTIFF_READ                                */
/*                                                                           */
/*  Read TIF data from Prarie 2-photon software.                             */
/*                                                                           */
/*****************************************************************************/
struct wytif_struct *wytiff_read(infile,mflag,pflag) // ,rt,rw,rh,rm,rdata)
     char infile[];
     int mflag;     // Mode:  0-Read all headers only
                    //        1-read one image
                    //        2-read all images, including data
     int pflag;
{
  FILE *fopen(),*fin;
  int i,j;
  int t,w,h,m,nret,nr,tn;
  int revflag;
  unsigned short ts;
  //unsigned long ifd1;
  unsigned int ifd1;
  char temp[SLEN];
  char c1,c2;
  struct wytif_struct *tiffs,*twt,*wt;
  float **fdata;

  if ((fin = fopen(infile,"r")) == NULL){
    printf("  *** File %s\n",infile);
    exit_error("WYTIFF_READ","Cannot open file");
  }

  //
  //  Endian check
  //
  nr = fread(&c1,1,1,fin);
  nr = fread(&c2,1,1,fin);

  if ((c1 == 'I') && (c2 == 'I')){
    if (pflag) printf("  First two bytes are 'II' - Little Endian (Intel).\n");
    revflag = 0;
  }else if ((c1 == 'M') && (c2 == 'M')){
    if (pflag) printf("  First two bytes are 'MM' - Big Endian (Motorola).\n");
    revflag = 1;
  }else{
    printf("  First two bytes are '%c%c'\n",c1,c2);
    exit_error("WYTIFF_READ","Unexpected first two bytes, neither II nor MM");
  }
  nret = revflag_fread(&ts,2,1,fin,0);
  if (pflag) printf("N-returned = %d\n",nret);
  //revflag_fread(&ts,2,1,fin,revflag);
  if (pflag) printf("This should be 42:  %d\n",(int)ts);
  if (ts != 42){
    printf(" *** found %d\n",(int)ts);
    exit_error("WYTIFF_READ","Expecting 42");
  }

  //
  //  Read offset of first IFD (Image File Directory)
  //
  nret = revflag_fread(&ifd1,4,1,fin,revflag);
  if (pflag) printf("N-returned = %d\n",nret);
  //printf("IFD 1:   %ld\n",ifd1);
  if (pflag) printf("IFD 1 offset:   %d\n",ifd1);

  fseek(fin,(long)ifd1,SEEK_SET);  // Seek rel. to begin
  tiffs = wytiff_get_empty_struct();
  wytiff_read_ifd(fin,revflag,tiffs,pflag);
  if (mflag > 0)
    tiffs->data = wytiff_read_data_float(fin,revflag,tiffs);

  tn = 1;
  if ((mflag == 0)||(mflag == 2)){
    //
    //  Read all IFD header data
    //
    twt = tiffs;
    while(twt->ifd_next > 0){
      fseek(fin,(long)twt->ifd_next,SEEK_SET);  // Seek rel. to begin
      wt = wytiff_get_empty_struct();
      wytiff_read_ifd(fin,revflag,wt,pflag);
      if (mflag == 2)
	wt->data = wytiff_read_data_float(fin,revflag,wt);

      twt->next = wt;
      twt = wt;
      tn += 1;
    }
  }

  if (pflag) // Debug
    wytiff_print_struct(tiffs);

  //printf("  Total images (IFDs):  %d\n",tn);

  //
  //  Read the data from the strips.
  //
  if (mflag > 0){
    tiffs->data = wytiff_read_data_float(fin,revflag,tiffs);
  }

  fclose(fin);


  /*
  {
    int transpose,printflag,tcode,w,h;

    w = tiffs->w;
    h = tiffs->h;
    transpose = 1;  // This shows the data as 'ImageJ' does.
    printflag = 1;
    tcode = 2;      // float
    write_2d_data("zzz.2d",fdata,0,0,w,h,4,tcode,transpose,printflag);
  }
  */

  return tiffs;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MYXML_GET_TAG_FROM_STRING                       */
/*                                                                           */
/*  's' should be of the form "...<xyz ..." and "xyz" will be returned.      */
/*                         or "...</xyz .." and "/xyz" will be returned.     */
/*                                                                           */
/*****************************************************************************/
char *myxml_get_tag_from_string(s)
     char *s;
{
  char *cp,*tag;

  cp = strchr(s,'<');
  tag = get_first_item_from_string(cp+1); // This creates storage

  cp = strchr(tag,'>');  // If there is a '>', then end the string there
  if (cp != NULL)
    *cp = '\0';    
    
  string_make_lower_case(tag);

  return tag;
}
/**************************************-**************************************/
/*                                                                           */
/*                     MYXML_GET_TAG_FROM_STRING_ADV_INDEX                   */
/*                                                                           */
/*  's' should be of the form "...<xyz ..." and "xyz" will be returned.      */
/*                         or "...</xyz .." and "/xyz" will be returned.     */
/*                                                                           */
/*  Advance the index along the string.                                      */
/*                                                                           */
/*****************************************************************************/
char *myxml_get_tag_from_string_adv_index(s,ri)
     char *s;  // Character string containing tag
     int *ri;  // Index currently points to '<' to be advanced to '>'
{
  int k;
  int done,nmax;
  char *c0,*ce,*cp,*tag;

  //
  //  Verify that we are at '<'
  //
  if (s[*ri] != '<'){
    printf("*** Index %d in String:  %s\n",*ri,s);
    exit_error("MYXML_GET_TAG_FROM_STRING_ADV_INDEX","Index not at '<'");
  }
  c0 = &(s[*ri+1]);

  //
  //  Find the matching '>'
  //
  ce = strchr(c0,'>');  // Include <string.h>
  if (ce == NULL){
    printf("*** Index %d in String:  %s\n",*ri,s);
    exit_error("MYXML_GET_TAG_FROM_STRING_ADV_INDEX","No matching '>'");
  }

  //
  //  Create storage for tag string
  //
  nmax = ce - c0 + 1;
  tag = (char *)myalloc(nmax*sizeof(char));

  cp = c0;
  k = 0;
  done = 0;
  while(done == 0){
    tag[k] = *cp;
    k += 1;
    cp += 1;
    if ((*cp == '>')||(*cp ==' ')||(*cp == '/')){
      done = 1;
    }
    // *** WYETH - NO END OF STRING CONDITION CHECKED HERE, CAN OVERRUN
  }
  tag[k] = '\0';

  *ri += (ce - c0) + 1;  // Adjust pointer

  return tag;
}
/**************************************-**************************************/
/*                                                                           */
/*                   MYXML_GET_TAGEND_FROM_STRING_ADV_INDEX                  */
/*                                                                           */
/*  's' should be of the form "...<xyz ..." and "xyz" will be returned.      */
/*                         or "...</xyz .." and "/xyz" will be returned.     */
/*                                                                           */
/*  Advance the index along the string.                                      */
/*                                                                           */
/*****************************************************************************/
char *myxml_get_tagend_from_string_adv_index(s,ri)
     char *s;  // Character string containing tag
     int *ri;  // Index currently points to '<' to be advanced to '>'
{
  int k;
  int done,nmax;
  char *c0,*ce,*cp,*tag;

  //
  //  Verify that we are at '/'
  //
  if (s[*ri] != '/'){
    printf("*** Index %d in String:  %s\n",*ri,s);
    exit_error("MYXML_GET_TAGEND_FROM_STRING_ADV_INDEX","Index not at '/'");
  }
  c0 = &(s[*ri+1]);

  //
  //  Find the matching '>'
  //
  ce = strchr(c0,'>');  // Include <string.h>
  if (ce == NULL){
    printf("*** Index %d in String:  %s\n",*ri,s);
    exit_error("MYXML_GET_TAGEND_FROM_STRING_ADV_INDEX","No matching '>'");
  }

  //
  //  Create storage for tag string
  //
  nmax = ce - c0 + 1;
  if (nmax > 0)
    tag = (char *)myalloc(nmax*sizeof(char));
  else{
    return NULL;
  }

  //printf("nmax = %d  ri=%d  ==>%s<== \n",nmax,*ri,s);

  cp = c0;
  k = 0;
  done = 0;
  while(done == 0){
    tag[k] = *cp;
    k += 1;
    cp += 1;
    if ((*cp == '>')||(*cp ==' ')||(*cp == '/')){
      done = 1;
    }
    // *** WYETH - NO END OF STRING CONDITION CHECKED HERE, CAN OVERRUN
  }
  tag[k] = '\0';

  *ri += (ce - c0) + 1;  // Adjust pointer

  return tag;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MYXML_GET_TAG_ASSIGNMENTS                        */
/*                                                                           */
/*  's' should be of the form '...name1=val1 name2=val2 name3="val 3"..."    */
/*                                                                           */
/*  A list of names and values are returned for all assignments, and the     */
/*  number of assignments is returned.  If there are zero, then the          */
/*  lists are set to NULL.                                                   */
/*                                                                           */
/*****************************************************************************/
int myxml_get_tag_assignments(s,rname,rval)
     char *s;               // string containing tag
     char ***rname,***rval; // returned lists of names and values, or NULL
{
  int i,j,k;
  int n,sn,done;
  char **name,**val,tstr[SLEN_MAX];

  n = 0;  // Number of assignments found in tag
  sn = strlen(s);

  // Count the number of '=' that are outside of double quotes
  i = 0;
  while(i < sn){
    if (s[i] == '='){        // Found an assignment
      n += 1;
      i += 1;
    }else if (s[i] == '"'){  // Found opening quotes, skip to closing quotes
      i += 1;
      done = 0;
      while(!done){
	if (i >= sn){
	  exit_error("MYXML_GET_TAG_ASSIGNMENTS","No closing quotes");
	}
	if (s[i] == '"')
	  done = 1;          // Found the closing quote
	i += 1;
      }
    }else{
      i += 1;                // Go to next character
    }
  }

  if (n > 0){
    name = (char **)myalloc(n*sizeof(char *));
    val  = (char **)myalloc(n*sizeof(char *));
  }else{
    name = NULL;
    val = NULL;
  }

  k = 0;
  while(s[k] != '<')  // Find the opening bracket
    k += 1;

  for(i=0;i<n;i++){
    while(s[k] != ' ')  // Move to the next space
      k += 1;
    while(s[k] == ' ')  // Move to the next non-space
      k += 1;
    j = 0;
    while(s[k] != '='){
      if (s[k] != ' '){  // Ignore spaces before '='
	tstr[j] = s[k];
	j += 1;
	if (j >= SLEN_MAX-1)
	  exit_error("MYXML_GET_TAG_ASSIGNMENTS","String too long");
      }
      k += 1;
    }
    tstr[j] = '\0';         // End the string
    name[i] = strdup(tstr); // Save the name

    // s[k] should now be '='
    k += 1;
    while(s[k] == ' ')  // skip any spaces
      k += 1;

    j = 0;
    if (s[k] == '"'){
      k += 1;
      while(s[k] != '"'){
	tstr[j] = s[k];
	k += 1;
	j += 1;
	if (j >= SLEN_MAX-1)
	  exit_error("MYXML_GET_TAG_ASSIGNMENTS","String too long");
      }
    }else{
      done = 0;
      while(!done){
	if (s[k] == ' ')
	  done = 1;
	else if (s[k] == '>')
	  done = 1;
	else{
	  tstr[j] = s[k];
	  j += 1;
	  k += 1;
	  if (j >= SLEN_MAX-1)
	    exit_error("MYXML_GET_TAG_ASSIGNMENTS","String too long");
	}
      }
    }

    tstr[j] = '\0';        // End the string
    val[i] = strdup(tstr); // Save the value
  }

  *rname = name;
  *rval = val;
  return n;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MYXML_GET_VALUE_FROM_STRING                      */
/*                                                                           */
/*  's' should be of the form '...name="xYz"...' and "xYz" will be returned. */
/*                                                                           */
/*  ***** WYETH - this is a HACK                                             */
/*                                                                           */
/*****************************************************************************/
char *myxml_get_value_from_string(s,name)
     char *s;
     char name[];
{
  char *cp,*ct,*valstr;

  cp = strstr(s,name);
  if (cp != NULL){
    ct = cp + strlen(name);  // Skip over name
    cp = strchr(ct,'=');     // Verify '='
    if (cp == NULL)
      exit_error("MYXML_GET_TAG_FROM_STRING","No '=' found");
    cp = strchr(ct,'"');     // Find opening quote
    valstr = strdup(cp+1);
    cp = strchr(valstr,'"');
    *cp = '\0';              // Terminate string at next quote
  }else
    valstr = (char *)NULL;

  return valstr;
}
