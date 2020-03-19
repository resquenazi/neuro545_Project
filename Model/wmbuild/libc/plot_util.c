/*****************************************************************************/
/*                                                                           */
/*  plot_util.c                                                              */
/*  wyeth bair                                                               */
/*  caltech                                                                  */
/*  04/15/95                                                                 */
/*                                                                           */
/*  This file contains utilities to write plots.  Other utilities            */
/*  already exist in other places:  (1) myplot.c (2) farray_util.c           */
/*  (3) iarray_util.c.                                                       */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "my_util.h"
#include "farray_util.h"

/**************************************-**************************************/
/*                                                                           */
/*                       WRITE_XPLOT_SARRAY_SAMPLING                         */
/*                                                                           */
/*  Write a spike "raster plot" for Matt Wilson's xplot.                     */
/*                                                                           */
/*****************************************************************************/
void write_xplot_sarray_sampling(outfile,s,cnt,n,offset,period,sampling)
     char outfile[];
     int **s,*cnt,n,offset,period;
     float sampling;
{
  int i,j;
  int t;
  FILE *fout;

  fout = my_fopen(outfile,"w","WRITE_XPLOT_SARRAY_SAMPLING");
  fprintf(fout,"/graphtitle %s\n",outfile);
  for(i=0;i<n;i++){
    fprintf(fout,"/newplot\n");
    fprintf(fout,"/plotname %d\n",i);
    fprintf(fout,"%d %d\n",offset,i);
    for(j=0;j<cnt[i];j++){
      if ((s[i][j] >= offset)&&(s[i][j] < offset+period)){
	t = (int)(0.5 + (float)s[i][j]*1000.0/sampling);
	fprintf(fout,"%d %d\n",t,i);
	fprintf(fout,"%d %.1f\n",t,(float)i+0.8);
	fprintf(fout,"%d %d\n",t,i);
      }    
    }
    fprintf(fout,"%d %d\n",(int)((float)(offset+period)*1000.0/sampling),i);
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                     WRITE_XPLOT_SARRAY_SAMPLING_COLOR                     */
/*                                                                           */
/*  Write a spike "raster plot" for Matt Wilson's xplot.                     */
/*                                                                           */
/*****************************************************************************/
void write_xplot_sarray_sampling_color(outfile,s,cnt,n,offset,period,sampling,
				       color)
     char outfile[];
     int **s,*cnt,n,offset,period;
     float sampling;
     int color;
{
  int i,j;
  int t;
  FILE *fout;

  fout = my_fopen(outfile,"w","WRITE_XPLOT_SARRAY_SAMPLING_COLOR");
  fprintf(fout,"/graphtitle %s\n",outfile);
  for(i=0;i<n;i++){
    fprintf(fout,"/newplot\n");
    fprintf(fout,"/plotname %d\n",i);
    fprintf(fout,"/color %d\n",color);
    fprintf(fout,"%d %d\n",offset,i);
    for(j=0;j<cnt[i];j++){
      if ((s[i][j] >= offset)&&(s[i][j] < offset+period)){
	t = (int)(0.5 + (float)s[i][j]*1000.0/sampling);
	fprintf(fout,"%d %d\n",t,i);
	fprintf(fout,"%d %.1f\n",t,(float)i+0.8);
	fprintf(fout,"%d %d\n",t,i);
      }    
    }
    fprintf(fout,"%d %d\n",(int)((float)(offset+period)*1000.0/sampling),i);
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                             WRITE_XPLOT_SARRAY                            */
/*                                                                           */
/*  Write a spike "raster plot" for Matt Wilson's xplot.                     */
/*                                                                           */
/*****************************************************************************/
void write_xplot_sarray(outfile,s,cnt,n,offset,period)
     char outfile[];
     int **s,*cnt,n,offset,period;
{
  int i,j;
  FILE *fout;

  fout = my_fopen(outfile,"w","WRITE_XPLOT_SARRAY");
  fprintf(fout,"/graphtitle %s\n",outfile);
  for(i=0;i<n;i++){
    fprintf(fout,"/newplot\n");
    fprintf(fout,"/plotname %d\n",i);
    fprintf(fout,"%d %d\n",offset,i);
    for(j=0;j<cnt[i];j++){
      if ((s[i][j] >= offset)&&(s[i][j] < offset+period)){
	fprintf(fout,"%d %d\n",s[i][j],i);
	fprintf(fout,"%d %.1f\n",s[i][j],(float)i+0.8);
	fprintf(fout,"%d %d\n",s[i][j],i);
      }    
    }
    fprintf(fout,"%d %d\n",offset+period,i);
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                        APPEND_XPLOT_SPIKE_SAMPLING                        */
/*                                                                           */
/*  Append one spike "raster plot" for Matt Wilson's xplot.                  */
/*                                                                           */
/*****************************************************************************/
void append_xplot_spike_sampling(outfile,name,s,n,offset,period,yoff,sampling,
				 t0)
     char outfile[],name[];
     int *s,n,offset,period,yoff;
     float sampling;
     int t0;
{
  int i;
  FILE *fout;

  fout = my_fopen(outfile,"a","APPEND_XPLOT_SPIKE_SAMPLING");
  fprintf(fout,"/newplot\n");
  fprintf(fout,"/plotname %s\n",name);
  fprintf(fout,"%d %d\n",offset-t0,yoff);
  for(i=0;i<n;i++){
    if ((s[i] >= offset)&&(s[i] < offset+period)){
      fprintf(fout,"%d %d\n",(int)((float)(s[i]-t0)*1000.0/sampling),yoff);
      fprintf(fout,"%d %.1f\n",(int)((float)(s[i]-t0)*1000.0/sampling),
	      (float)yoff+0.8);
      fprintf(fout,"%d %d\n",(int)((float)(s[i]-t0)*1000.0/sampling),yoff);
    }
  }    
  fprintf(fout,"%d %d\n",(int)((float)(offset-t0+period)*1000.0/sampling),
	  yoff);

  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                        APPEND_XPLOT_SARRAY_SAMPLING                       */
/*                                                                           */
/*  Append a spike "raster plot" for Matt Wilson's xplot.                    */
/*                                                                           */
/*****************************************************************************/
void append_xplot_sarray_sampling(outfile,names,s,cnt,n,offset,period,yoff,
				  sampling)
     char outfile[],*names[];
     int **s,*cnt,n,offset,period,yoff;
     float sampling;
{
  int i,j;
  FILE *fout;

  fout = my_fopen(outfile,"a","APPEND_XPLOT_SARRAY_SAMPLING");
  for(i=0;i<n;i++){
    fprintf(fout,"/newplot\n");
    if (names != NULL)
      fprintf(fout,"/plotname %s\n",names[i]);
    else
      fprintf(fout,"/plotname %d\n",i);
    fprintf(fout,"%d %d\n",offset,i+yoff);
    for(j=0;j<cnt[i];j++){
      if ((s[i][j] >= offset)&&(s[i][j] < offset+period)){
	fprintf(fout,"%d %d\n",(int)((float)s[i][j]*1000.0/sampling),i+yoff);
	fprintf(fout,"%d %.1f\n",(int)((float)s[i][j]*1000.0/sampling),
		(float)(i+yoff)+0.8);
	fprintf(fout,"%d %d\n",(int)((float)s[i][j]*1000.0/sampling),i+yoff);
      }
    }    
    fprintf(fout,"%d %d\n",(int)((float)(offset+period)*1000.0/sampling),
				 i+yoff);
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                     APPEND_XPLOT_SARRAY_SAMPLING_COLOR                    */
/*                                                                           */
/*  Append a spike "raster plot" for Matt Wilson's xplot.                    */
/*                                                                           */
/*****************************************************************************/
void append_xplot_sarray_sampling_color(outfile,names,s,cnt,n,offset,period,
					yoff,sampling,color)
     char outfile[],*names[];
     int **s,*cnt,n,offset,period,yoff;
     float sampling;
     int color;
{
  int i,j;
  FILE *fout;

  fout = my_fopen(outfile,"a","APPEND_XPLOT_SARRAY_SAMPLING_COLOR");
  for(i=0;i<n;i++){
    fprintf(fout,"/newplot\n");
    fprintf(fout,"/color %d\n",color);
    if (names != NULL)
      fprintf(fout,"/plotname %s\n",names[i]);
    else
      fprintf(fout,"/plotname %d\n",i);
    fprintf(fout,"%d %d\n",offset,i+yoff);
    for(j=0;j<cnt[i];j++){
      if ((s[i][j] >= offset)&&(s[i][j] < offset+period)){
	fprintf(fout,"%d %d\n",(int)((float)s[i][j]*1000.0/sampling),i+yoff);
	fprintf(fout,"%d %.1f\n",(int)((float)s[i][j]*1000.0/sampling),
		(float)(i+yoff)+0.8);
	fprintf(fout,"%d %d\n",(int)((float)s[i][j]*1000.0/sampling),i+yoff);
      }
    }    
    fprintf(fout,"%d %d\n",(int)((float)(offset+period)*1000.0/sampling),
				 i+yoff);
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                            APPEND_XPLOT_SARRAY                            */
/*                                                                           */
/*  Append a spike "raster plot" for Matt Wilson's xplot.                    */
/*                                                                           */
/*****************************************************************************/
void append_xplot_sarray(outfile,s,cnt,n,offset,period,yoff)
     char outfile[];
     int **s,*cnt,n,offset,period,yoff;
{
  int i,j;
  FILE *fout;
  
  fout = my_fopen(outfile,"a","APPEND_XPLOT_SARRAY");
  for(i=0;i<n;i++){
    fprintf(fout,"/newplot\n");
    fprintf(fout,"/plotname %d\n",i);
    fprintf(fout,"%d %d\n",offset,i+yoff);
    for(j=0;j<cnt[i];j++){
      if ((s[i][j] >= offset)&&(s[i][j] < offset+period)){
	fprintf(fout,"%d %d\n",s[i][j],i+yoff);
	fprintf(fout,"%d %.1f\n",s[i][j],(float)(i+yoff)+0.8);
	fprintf(fout,"%d %d\n",s[i][j],i+yoff);
      }
    }    
    fprintf(fout,"%d %d\n",offset+period,i+yoff);
  }
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                        APPEND_XPLOT_FARRAY_HISTOGAM                       */
/*                                                                           */
/*****************************************************************************/
void append_xplot_farray_histogram(outfile,name,data,n,binsize,start,yoff)
     char outfile[],name[];
     float *data;
     int n;
     float binsize,start,yoff;
{
  FILE *fout;
  int i;
  float x;

  fout = my_fopen(outfile,"a","APPEND_XPLOT_FARRAY_HISTOGRAM");
  fprintf(fout,"/newplot\n");
  fprintf(fout,"/plotname %s\n",name);
  fprintf(fout,"%.4e %.4e\n",start,yoff);
  for(i=0;i<n;i++){
    x = start + (float)i*binsize;
    fprintf(fout,"%.4e %.4e\n",x,data[i]+yoff);
    fprintf(fout,"%.4e %.4e\n",x+binsize,data[i]+yoff);
    fprintf(fout,"%.4e 0.0\n",x+binsize+yoff);
  }
  fprintf(fout,"%.4e %.4e\n",start,yoff);
  
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                           APPEND_XPLOT_FARRAY_DX                          */
/*                                                                           */
/*****************************************************************************/
void append_xplot_farray_dx(outfile,name,data,n,dx,start,yoff)
     char outfile[],name[];
     float *data;
     int n;
     double dx;
     float start,yoff;
{
  FILE *fout;
  int i;

  fout = my_fopen(outfile,"a","APPEND_XPLOT_FARRAY_DX");
  fprintf(fout,"/newplot\n");
  fprintf(fout,"/plotname %s\n",name);
  for(i=0;i<n;i++)
    fprintf(fout,"%.6e %.4e\n",start+(float)((double)i*dx),data[i]+yoff);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                        APPEND_XPLOT_FARRAY_DX_INT                         */
/*                                                                           */
/*****************************************************************************/
void append_xplot_farray_dx_int(outfile,name,data,n,dx,start,yoff)
     char outfile[],name[];
     float *data;
     int n;
     int dx;
     int start;
     float yoff;
{
  FILE *fout;
  int i;

  fout = my_fopen(outfile,"a","APPEND_XPLOT_FARRAY_DX_INT");
  fprintf(fout,"/newplot\n");
  fprintf(fout,"/plotname %s\n",name);
  for(i=0;i<n;i++)
    fprintf(fout,"%d %.4e\n",start+i*dx,data[i]+yoff);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                           APPEND_XPLOT_CONST_PLOT                         */
/*                                                                           */
/*  Plot a horizontal line of the specified value.                           */
/*                                                                           */
/*****************************************************************************/
void append_xplot_const_plot(outfile,name,x0,x1,y0)
     char *outfile,*name;
     float x0,x1,y0;
{
  FILE *fout;

  fout = my_fopen(outfile,"a","APPEND_XPLOT_CONST_PLOT");
  fprintf(fout,"/newplot\n");
  if (name != NULL)
    fprintf(fout,"/plotname %s\n",name);
  fprintf(fout,"%.6e %.6e\n",x0,y0);
  fprintf(fout,"%.6e %.6e\n",x1,y0);
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                           APPEND_XPLOT_CONST_PLOT                         */
/*                                                                           */
/*  Append automatically generated histograms for two data arrays to an      */
/*  output file.                                                             */
/*                                                                           */
/*****************************************************************************/
void append_xplot_auto_hist_pair(outfile,name1,name2,d1,d2,n1,n2,nbin)
     char *outfile;       // File to append
     char *name1,*name2;  // Names of two plots
     float *d1,*d2;       // Data to be histogrammed [n1], [n2]
     int n1,n2;           // Length of data
     int nbin;            // Number of bins
{
  int npob;
  float dmin1,dmin2,dmax1,dmax2,dmin,dmax,drange,*hh,binsize;

  get_min_max_farray(d1,n1,&dmin1,&dmax1);
  get_min_max_farray(d2,n2,&dmin2,&dmax2);
  dmin = myminf(dmin1,dmin2);
  dmax = mymaxf(dmax1,dmax2);
  drange = dmax - dmin;
  dmax = dmax + 0.1*drange;
  dmin = dmin - 0.1*drange;

  hh = get_histogram_farray(d1,n1,dmin,dmax,nbin,&binsize,&npob);
  if (npob > 0) printf("  H_add:  %d points out of bounds\n",npob);
  append_xplot_farray_histogram(outfile,name1,hh,nbin,binsize,dmin,0.0);
  myfree(hh);

  hh = get_histogram_farray(d2,n2,dmin,dmax,nbin,&binsize,&npob);
  if (npob > 0) printf("  H_add:  %d points out of bounds\n",npob);
  append_xplot_farray_histogram(outfile,name2,hh,nbin,binsize,dmin,0.0);
  myfree(hh);
}
