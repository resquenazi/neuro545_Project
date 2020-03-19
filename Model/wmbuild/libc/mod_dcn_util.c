/*****************************************************************************/
/*                                                                           */
/*   mod_dcn_util.c                                                          */
/*   wyeth bair                                                              */
/*                                                                           */
/*   Deep Convolutional Network                                              */
/*   Designed to implement AlexNet to match Caffe model.                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>  // For sleep
#include <math.h>
#include <string.h>

// user defined functions
#include "my_util.h"
#include "nr_util.h"
#include "iarray_util.h"
#include "carray_util.h"
#include "farray_util.h"
#include "myrand_util.h"
#include "data_util.h"
#include "min_util.h"
#include "paramfile_util.h"
#include "mod_util.h"
#include "ifc.h"           // For .mar file
#include "mod.h"           // Data structures
#include "paramfile.h"

// For cluster computation  (static needed w/ initialization)
static int myid = -1;         // Process ID, -1-Single, 0-Master, >0-Slave
static int numproc = -1;      // Number of processors (MM)
static char *mylogf = NULL;   // Logfile name
char ggstr[LONG_SLEN];        // Global log string

int mod_dcn_action;           // Action value from caller

struct dcn_layer_struct{
  char *type;   // 'conv', 'relu', ...
  char *name;   // name of layer
  char *top;    // name of layer
  char *bot;    // name of layer

  char *pool;   // Type of pooling, 'MAX',...

  int  nkern;   // Number of kernels in conv layer; No. units in FC layer
  int  kx;      // width of kernel (pix)
  int  ky;      // height of kernel (pix)
  int  kz;      // depth of kernel (planes, e.g., RGB is 3)
  int  stride;  // stride (default 1)
  int  pad;     // padding (default 0)
  int  group;   // number of groups (default 1)

  // This gets computed from the layer parameters above
  int  out_xn;  // width of layer output
  int  out_yn;  // height of layer output
  int  out_zn;  // Depth of output

  float ****kern;  // [nkern][kz][kx][ky]  for Convolution
                   // [1][1][kx][ky]      for InnerProduct (fully connected)
  float *bias;     // [nkern], or [kx] for InnerProduct

  float ***kern1;  // [nkern][kx][ky] If input is 'data', then create a
                   //    set of 2D kernels that sum over the z (RGB) axis.

  float ***r;      // Response of conv layer [out_zn][out_xn][out_yn]
                   //             FC layer   [1][1][out_yn]

  float ****rt;    // Like 'r', but store for multiple time points. ...[t]

  int loc_size; // LRN param
  float alpha;  // LRN param
  float beta;   // LRN param

  //
  //  Visualization
  //
  int vxi,vyi,vzi;    // Index of unit to visualize, or -1
  int ***vis_pool_x;  // [out_zn][out_xn][out_yn] x-coord (in input) of max
  int ***vis_pool_y;  // [out_zn][out_xn][out_yn] y-coord (in input) of max
};


//
//  Global variables - model specific
//
int         mod_dcn_xn;      // Width of stimulus (pix)
int         mod_dcn_yn;      // Height of stimulus (pix)
int         mod_dcn_zn;      // Stimulus depth (e.g., 3 for RGB)
int         mod_dcn_tn;      // Duration of stimulus (frames)
float       mod_dcn_tscale;  // (s/frame)
float       mod_dcn_sscale;  // (deg/pix)

char       *mod_dcn_stimform;

char       *mod_dcn_vistype;    // Type of visualization "ZF14" [NULL]
char       *mod_dcn_vismode;    // "max" - find maximum response in layer
                                // "unit" - use specified unit
char       *mod_dcn_vislay;     // Target layer for visualization
int         mod_dcn_visu_xi;    // Target unit coordinates
int         mod_dcn_visu_yi;    // Target unit coordinates
int         mod_dcn_visu_zi;    // Target unit coordinates
char       *mod_dcn_visopre;    // Output file prefix
int         mod_dcn_viso3d;     // 1-write .3D file output
float    ***mod_dcn_vis_rec;    // Reconstructed stimulus

int mod_dcn_lay_n;  // Current number of built layers
struct dcn_layer_struct **mod_dcn_lay;    // [mod_dcn_lay_n] all model layers

int mod_dcn_maxli;  // Maximum layer to compute responses

float ***mod_dcn_tstim;  // Test stimulus [zn][xn][yn]

int mod_dcn_rsize;  // Size of response for all layers

char *mod_dcn_resp_3d_fname;

char *mod_dcn_wmover;           // WM override mode, 'max' or NULL
char *mod_dcn_wmover_flist;     // Name of file containing stimulus list
char *mod_dcn_wmover_fpath;     // Path to image files
char *mod_dcn_wmover_ofile;     // Name of output file to append stats

/**************************************-**************************************/
/*                                                                           */
/*                            MOD_DCN_LAYI_FOR_NAME                          */
/*                                                                           */
/*  Return the index of the layer w/i 'mod_dcn_lay' given the 'name'.        */
/*  Or, return -1 if not found.                                              */
/*                                                                           */
/*****************************************************************************/
int mod_dcn_layi_for_name(name)
     char *name;
{
  int i,k;

  k = -1;
  i = 0;
  while((i < mod_dcn_lay_n) && (k == -1)){
    if (strcmp(mod_dcn_lay[i]->name,name)==0)
      k = i;
    i += 1;
  }

  return k;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_DCN_LAY_FOR_NAME                           */
/*                                                                           */
/*  Return a pointer to the layer structure for the given name.              */
/*  Or, return NULL if not found.                                            */
/*                                                                           */
/*****************************************************************************/
struct dcn_layer_struct *mod_dcn_lay_for_name(name)
     char *name;
{
  int i;
  struct dcn_layer_struct *t;

  i = mod_dcn_layi_for_name(name);

  if (i < 0)
    t = NULL;
  else
    t = mod_dcn_lay[i];

  return t;
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_DCN_STAT_KERN                             */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_stat_kern(dl)
     struct dcn_layer_struct *dl;
{
  int i;
  int n,xn,yn,zn;
  float *kmu,mu,sd;

  printf("  MOD_DCN_STAT_KERN\n");

  n  = dl->nkern;
  xn = dl->kx;
  yn = dl->ky;
  zn = dl->kz;

  kmu = get_farray(n);

  for(i=0;i<n;i++)
    single_mean_3d_farray(dl->kern[i],0,zn,0,xn,0,yn,&(kmu[i]));

  mean_sdev_farray(kmu,n,&mu,&sd);
  printf("    %s   Mean, SD of kernel means:  %f %f\n",dl->name,mu,sd);

  mean_sdev_farray(dl->bias,n,&mu,&sd);
  printf("            Bias:  %f %f\n",mu,sd);

  myfree(kmu);
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_DCN_WRITE_RGB_KERNELS                        */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_write_rgb_kernels(outfile)
     char *outfile;
{
  int i,j,k;
  int n,xn,yn,***d,**t;
  float **fr,**fg,**fb,r_min,r_max,g_min,g_max,b_min,b_max,min,max;
  struct dcn_layer_struct *dl;

  dl = mod_dcn_lay[0];

  xn = dl->kx;
  yn = dl->ky;
  n  = dl->nkern;

  d = get_3d_iarray(xn,yn,n);

  for(i=0;i<n;i++){
    fr = dl->kern[i][0];  // [nkern][z][x][y]
    fg = dl->kern[i][1];  // [nkern][z][x][y]
    fb = dl->kern[i][2];  // [nkern][z][x][y]

    //
    //  Scale the RGB values to be in [0..1] based on the single
    //  min and max across all three color planes.
    //
    get_min_max_2d_farray(fr,xn,yn,&r_min,&r_max);
    get_min_max_2d_farray(fg,xn,yn,&g_min,&g_max);
    get_min_max_2d_farray(fb,xn,yn,&b_min,&b_max);
    min = min_of_three(r_min,g_min,b_min);
    max = max_of_three(r_max,g_max,b_max);

    add_const_2d_farray(fr,0,xn,0,yn,-min);
    add_const_2d_farray(fg,0,xn,0,yn,-min);
    add_const_2d_farray(fb,0,xn,0,yn,-min);

    multiply_2d_farray(fr,xn,yn,1.0/(max-min));
    multiply_2d_farray(fg,xn,yn,1.0/(max-min));
    multiply_2d_farray(fb,xn,yn,1.0/(max-min));

    // Add this red pixel to verify that the first dimension is 'x' and the
    // second is 'y' (horiz and vert, respectively)
    //fr[7][3] = 1.0;
    //fg[7][3] = 0.0;
    //fb[7][3] = 0.0;

    t = data_util_t11_2d_from_rgb(fr,fg,fb,xn,yn);

    //  Copy into x,y,t order
    for(j=0;j<xn;j++)
      for(k=0;k<yn;k++)
	d[j][k][i] = t[j][k];

    free_2d_iarray(t,xn);
  }

  write_3d_data_part(outfile,d,0,xn,0,yn,0,n,4,11,1); // 1=txy format

  free_3d_iarray(d,xn,yn,n);
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_DCN_WRITE_KERNEL                           */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_write_kernel(outfile,lname,ki)
     char *outfile;  // outfile name
     char *lname;    // layer name, e.g., conv2
     int ki;         // index of kernel in layer
{
  int i,j,k;
  int n,xn,yn,zn,cnt,fpr,mar,nw,nh,w,h,xi,yi;
  float **d,tot;
  struct dcn_layer_struct *dl;

  mylog(mylogf,"  MOD_DCN_WRITE_KERNEL\n");

  // Constants
  fpr = 10;  // At most 10 frames per row
  mar = 2;   // Pixel margin

  dl = mod_dcn_lay_for_name(lname);

  xn = dl->kx;
  yn = dl->ky;
  zn = dl->kz;
  n  = dl->nkern;

  if ((ki < 0) || (ki > dl->nkern)){
    printf("  Index:  %d is not in 0..%d\n",ki,dl->nkern-1);
    exit_error("MOD_DCN_WRITE__KERNEL","Bad kernel index");
  }

  printf("    Layer %s, kernel %d  (%d,%d,%d)\n",lname,ki,xn,yn,zn);

  nw = mymini(zn,fpr);  // Number of frames horizontally
  nh = (int)(zn / nw);  // Number of frames vertically
  if (nw*nh < zn)
    nh += 1;

  w = xn*nw + (nw-1)*mar; // Image width, 2 pix border between frames
  h = yn*nh + (nh-1)*mar; // Image height, 2 pix border between frames

  d = get_zero_2d_farray(w,h);

  printf("    %d x %d frames\n",nw,nh);
  printf("    %d x %d pix image\n",w,h);

  xi = 0;
  yi = h-yn;
  cnt = 0;
  tot = 0.0;
  for(k=0;k<zn;k++){

    for(i=0;i<xn;i++){  // Draw the 'k'th spatial frame into the image
      for(j=0;j<yn;j++){
	d[xi+i][yi+j] = dl->kern[ki][k][i][j];
	tot += dl->kern[ki][k][i][j];
      }
    }
    xi += xn+mar;

    cnt += 1;
    if (cnt >= fpr){
      yi -= (yn+mar);
      xi = 0;
      cnt = 0;
    }
  }

  printf("    Total:  %f    Avg:  %f\n",tot,tot/(float)(zn*xn*yn));

  write_2d_data(outfile,d,0,0,w,h,4,2,1,1); // No transp for consistency
  free_2d_farray(d,w);
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_DCN_MAP_FC                               */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_map_fc(outfile,lname,ki)
     char *outfile;  // outfile name
     char *lname;    // layer name, e.g., conv2
     int ki;         // index of unit to map in layer
{
  int i,j,k;
  int li;
  struct dcn_layer_struct *dl,*dl2;

  mylog(mylogf,"  MOD_DCN_MAP_FC\n");

  li = mod_dcn_layi_for_name(lname);

  printf("    Layer %s  lay_i %d  unit %d\n",lname,li,ki);

  dl2 = mod_dcn_lay[li+2];    // [mod_dcn_lay_n] all model layers
  printf("    Layer %s  lay_i %d\n",dl2->name,li+2);

  // For each unit
  //   Make 64x64 array of weights
  //   Compute sum of abs of weights
  //   scale each entry to be % of abs of weights
  //   

  /***  WYETH UNDER CONSTRUCTION - - NEVRER COMPLETED ***/
  /***  WYETH UNDER CONSTRUCTION - - NEVRER COMPLETED ***/


  dl = mod_dcn_lay_for_name(lname);

  exit(0);
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_DCN_WRITE_3D_RESP                          */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_write_3d_resp(outfile,dl)
     char *outfile;
     struct dcn_layer_struct *dl;
{
  int i,j,k;
  int zn,xn,yn;
  float ***d,**t;

  xn = dl->out_xn;
  yn = dl->out_yn;
  zn = dl->out_zn;

  if ((xn > 1) && (yn > 1)){

    d = get_3d_farray(xn,yn,zn);

    //  Copy into x,y,t order
    for(i=0;i<zn;i++){
      t = dl->r[i];
      for(j=0;j<xn;j++)
	for(k=0;k<yn;k++)
	  d[j][k][i] = t[j][k];
    }

    write_3d_data_part(outfile,d,0,xn,0,yn,0,zn,4,2,1); // 1=txy format

    free_3d_farray(d,xn,yn,zn);
  }
/*
  // print out the lower left corner of 0th channel
  printf("_______________ %s _________________\n",dl->name);


  if ((xn > 1) && (yn > 1)){
    //
    //  Print a 5x5 patch
    //
    t = dl->r[0];
    for(j=yn-1;j>(yn-1-5);j--){
      for(i=0;i<5;i++){
	printf(" %11.6f",t[i][j]);
      }
      printf("\n");
    }
  }else{
    //
    //  Print the first 20 - this is 1D data.
    //
    t = dl->r[0];
    for(i=0;i<20;i++){
      if (strcmp(dl->type,"Softmax")==0){
	printf(" %5e",t[0][i]);
      }else{
	printf(" %11.6f",t[0][i]);
      }
      if ((i+1)%5 == 0)
	printf("\n");
    }
  }
*/
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_DCN_READ_WEIGHTS                            */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_read_weights(infile)
     char *infile;  // Name of the input file
{
  FILE *fopen(),*fin;
  int i1,i2,i3,i4,n1,n2,n3,n4,sn,ndim,layi,ns;
  float ***kp,**kp1,t,*td;
  char temp[SLEN],*tc,**slist,*lname,*dtype;
  struct dcn_layer_struct *lay;

  mylog(mylogf,"  MOD_DCN_READ_WEIGHTS\n");

  if ((fin = fopen(infile,"r"))==NULL){
    printf("filename: %s\n",infile);
    exit_error("MOD_DCN_READ_WEIGHTS","Cannot open file");
  }

  tc = fgets(temp,SLEN,fin); // read one full line, including carriage return
  while(tc != NULL){
    printf("    %s",tc);  // Print the 'LAYER ... ' line

    get_items_from_string(tc,&slist,&sn);
    ndim = sn - 3;

    lname = slist[1];  // Pointer to layer name
    dtype = slist[2];  // Pointer to data type

    // Search for 'lname' in the list of layers
    layi = mod_dcn_layi_for_name(lname);
    if (layi > mod_dcn_maxli){
      // Stop reading from the file, the higher layers provide no responses.
      tc = NULL;
    }else{ // not done

      lay = mod_dcn_lay_for_name(lname);
      if (lay == NULL){
	printf("  Layer name:  %s\n",lname);
	exit_error("MOD_DCN_READ_WEIGHTS","Cannot find layer");
      }

      if (strcmp(dtype,"weights")==0){
	//
	//  Read weights
	//

	if (strcmp(lay->type,"Convolution")==0){
	  if (ndim != 4)
	    exit_error("MOD_DCN_READ_WEIGHTS","Convolution needs 4 values");
	  n1 = atoi(slist[3]);
	  n2 = atoi(slist[4]);
	  n3 = atoi(slist[5]);
	  n4 = atoi(slist[6]);

	  if (n1 != lay->nkern){
	    printf("  nkern = %d  n1 = %d\n",lay->nkern,n1);
	    exit_error("MOD_DCN_READ_WEIGHTS","nkernel does not match");
	  }
	  if (n2 != lay->kz){
	    printf("  kz = %d  n2 = %d\n",lay->kz,n2);
	    exit_error("MOD_DCN_READ_WEIGHTS","kz does not match");
	  }
	  if (n3 != lay->kx)
	    exit_error("MOD_DCN_READ_WEIGHTS","kx does not match");
	  if (n4 != lay->ky)
	    exit_error("MOD_DCN_READ_WEIGHTS","ky does not match");

	}else if (strcmp(lay->type,"InnerProduct")==0){
	  if (ndim != 2)
	    exit_error("MOD_DCN_READ_WEIGHTS","InnerProduct needs 2 values");
	  n1 = 1;
	  n2 = 1;
	  n3 = atoi(slist[3]);
	  n4 = atoi(slist[4]);

	  if (n3 != lay->kx)
	    exit_error("MOD_DCN_READ_WEIGHTS","InnProd: kx does not match");
	  if (n4 != lay->ky){
	    printf("n4 = %d   lay->ky = %d\n",n4,lay->ky);
	    exit_error("MOD_DCN_READ_WEIGHTS","InnProd: ky does not match");
	  }
	}

	printf("      %d dimensions\n",ndim);
	printf("      %d values to read\n",n1*n2*n3*n4);

	if (strcmp(lay->type,"InnerProduct")==0){

	  for(i3=0;i3<n3;i3++){  // For each unit
	    td = lay->kern[0][0][i3];
	    for(i4=0;i4<n4;i4++){  // read all of its weights
	      ns = fscanf(fin,"%f",&(td[i4]));
	    }
	  }

	}else{
	  for(i1=0;i1<n1;i1++)
	    for(i2=0;i2<n2;i2++)
	      for(i4=(n4-1);i4>=0;i4--){
		for(i3=0;i3<n3;i3++)
		  ns = fscanf(fin,"%f",&(lay->kern[i1][i2][i3][i4]));
	      }
	}

	if ((strcmp(lay->bot,"data")==0)&&(lay->kz == 3)){
	  //
	  //  If this convolutional layer works on RGB inputs, prepare flattened
	  //     kernels 'kern1' to quickly handle grayscale stimuli
	  //
	  printf("      ---> Preparing flattened RGB kernels.\n");

	  // Make a 2D kernel 'kern1' that combines the RGB layers
	  for(i1=0;i1<n1;i1++){
	    kp = lay->kern[i1];
	    kp1 = lay->kern1[i1];
	    for(i3=0;i3<n3;i3++)
	      for(i4=0;i4<n4;i4++)
		kp1[i3][i4] = kp[0][i3][i4] + kp[1][i3][i4] + kp[2][i3][i4];
	  }
	}

      }else if (strcmp(dtype,"biases")==0){
	//
	//  Read biases
	//

	n1 = atoi(slist[3]);

	if (strcmp(lay->type,"Convolution")==0){
	  if (n1 != lay->nkern)
	    exit_error("MOD_DCN_READ_WEIGHTS","Conv: bad bias count");
	}else if (strcmp(lay->type,"InnerProduct")==0){
	  if (n1 != lay->kx)
	    exit_error("MOD_DCN_READ_WEIGHTS","InnerProd: bad bias count");
	}

	printf("      %d values to read\n",n1);
	for(i1=0;i1<n1;i1++)
	  ns = fscanf(fin,"%f",&(lay->bias[i1]));

      }else{
	printf("  Data type = %s\n",dtype);
	exit_error("MOD_DCN_READ_WEIGHTS","Unknown data type");
      }

      tc = fgets(temp,SLEN,fin); // Discard the carriage return
      tc = fgets(temp,SLEN,fin); // read the next 'LAYER' line

    }
    free_2d_carray(slist,sn);
  }

  fclose(fin);
  mylog(mylogf,"    Done reading weight file.\n");
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_DCN_LAY_SET_SIZE                           */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_lay_set_size(dl,layo)
     struct dcn_layer_struct *dl;
     struct onode *layo;
{
  int zn,xn_in,yn_in;
  struct dcn_layer_struct *layin;

  //
  //  Compute the depth of the layer input 'zn'
  //  Compute the width of layer input 'xn_in' (needed below)
  //
  if (strcmp(dl->bot,"data")==0){
    zn = mod_dcn_zn;
    xn_in = mod_dcn_xn;
    yn_in = mod_dcn_yn;
  }else{
    layin = mod_dcn_lay_for_name(dl->bot);
    if (layin == NULL){
      printf("  Bottom (input) name:  %s\n",dl->bot);
      exit_error("MOD_DCN_LAY_SET_SIZE","Cannot find layer");
    }
    zn = layin->out_zn / dl->group;
    xn_in = layin->out_xn;
    yn_in = layin->out_yn;
  }

  dl->kz = zn;

  if (strcmp(dl->type,"Convolution")==0){
    dl->out_xn = (int)((xn_in + 2*dl->pad - dl->kx) / dl->stride) + 1;
    dl->out_yn = dl->out_xn;
    dl->out_zn = dl->nkern;  // Number of kernels is output depth
  }else if (strcmp(dl->type,"Pooling")==0){
    dl->out_xn = (int)((xn_in + 2*dl->pad - dl->kx) / dl->stride) + 1;
    dl->out_yn = dl->out_xn;
    dl->out_zn = zn;     // Default is that output matches input
  }else if (strcmp(dl->type,"InnerProduct")==0){
    dl->out_zn = 1;
    dl->out_xn = 1;
    dl->out_yn = dl->nkern;  // Number of units in full-connect layer

    // The y-size of the 'kernel' is the total number of outputs (x*y*z)
    //   from the previous layer.
    // The x-size is the number of units in this layer
    dl->kz = 1;
    dl->kx = dl->nkern;
    dl->ky = layin->out_xn * layin->out_yn * layin->out_zn;

  }else{
    dl->out_xn = xn_in;  // Default is that output matches input
    dl->out_yn = yn_in;
    dl->out_zn = zn;     // Default is that output matches input
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_DCN_LAY_COMMON                            */
/*                                                                           */
/*  Read and set parameters in the layer structure that are common to all    */
/*  (or most?) layers.                                                       */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_lay_common(dl,layo)
     struct dcn_layer_struct *dl;
     struct onode *layo;
{
  dl->name  = onode_getpar_chr_exit(layo,"name");
  dl->type  = onode_getpar_chr_exit(layo,"type");
  dl->bot   = onode_getpar_chr_exit(layo,"bottom");
  dl->top   = onode_getpar_chr_exit(layo,"top");

  // These may be over-written for some layers, e.g., "Convolution"
  dl->group  = 1;  // Default value
  dl->pad    = 0;  // Default value
  dl->stride = 1;  // Default value

  dl->kx     = 1;  // Default value
  dl->ky     = 1;  // Default value

  dl->vxi    = -1;
  dl->vyi    = -1;
  dl->vzi    = -1;

  dl->vis_pool_x = NULL;
  dl->vis_pool_y = NULL;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_DCN_PREP_RESPONSE                           */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_prep_response()
{
  int i;
  int xn,yn,zn;
  struct dcn_layer_struct *dl;

  mylog(mylogf,"  MOD_DCN_PREP_RESPONSE\n");

  for(i=0;i<=mod_dcn_maxli;i++){

    dl = mod_dcn_lay[i];

    if (strcmp(dl->type,"InnerProduct")==0){
      xn = 1;
      yn = dl->kx;
      zn = 1;
    }else{
      xn = dl->out_xn;
      yn = dl->out_yn;
      zn = dl->out_zn;
    }

    //
    //printf("  x,y,z  %d %d %d\n",xn,yn,zn);
    // This (above) printf reveals the following:
    //
    //           xn   yn   zn
    //   conv1   55   55   96
    //   ...
    //   pool5    6    6  256
    //   fc6      1 4096    1
    //   ...       
    //   relu7    1 4096    1
    //   prob     1 1000    1
    //

    dl->r  = get_3d_farray(zn,xn,yn);
    if (mod_dcn_tn > 1)
      dl->rt = get_4d_farray(zn,xn,yn,-1); // Don't allocate last dimension yet
    mod_dcn_rsize += zn * xn * yn;
    printf("    %6s %8d responses\n",dl->name,zn*xn*yn);
  }
  printf("    %6s %8d responses\n","Total",mod_dcn_rsize);
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_DCN_PREP_LAYTYPE_CONV                        */
/*                                                                           */
/*****************************************************************************/
struct dcn_layer_struct *mod_dcn_prep_laytype_conv(m,lay)
     struct model_struct *m;   // Model params
     struct onode *lay;
{
  struct onode *opar;
  struct dcn_layer_struct *dl;

  mylog(mylogf,"  MOD_DCN_PREP_LAYTYPE_CONV\n");

  dl = (struct dcn_layer_struct *)myalloc(sizeof(struct dcn_layer_struct));
  mod_dcn_lay_common(dl,lay);  // name, type, bottom, top, group, pad

  opar = onode_child_get_unique(lay,"param");
  if (opar == NULL)
    exit_error("MOD_DCN_PREP_LAYTYPE_CONV","Cannot find <param>.");

  dl->nkern  = onode_getpar_int_exit(opar,"num_output");
  dl->kx     = onode_getpar_int_exit(opar,"kernel_size");
  dl->ky     = dl->kx;
  dl->stride = onode_getpar_int_dflt(opar,"stride",1);
  dl->pad    = onode_getpar_int_dflt(opar,"pad",0);
  dl->group  = onode_getpar_int_dflt(opar,"group",1);

  mod_dcn_lay_set_size(dl,lay);  // Set sizes of inputs and outputs

  printf("    %s (layer %d)\n",dl->name,mod_dcn_lay_n);
  printf("    Kernels:  %d (%d,%d,%d)\n",dl->nkern,dl->kx,dl->ky,dl->kz);
  printf("    Stride:   %d\n",dl->stride);
  printf("    Padding:  %d\n",dl->pad);
  printf("    Output:   %d (%d,%d)\n",dl->out_zn,dl->out_xn,dl->out_yn);

  dl->kern = get_4d_farray(dl->nkern,dl->kz,dl->kx,dl->ky);
  dl->bias = get_farray(dl->nkern);

  if ((strcmp(dl->bot,"data")==0)&&(dl->kz == 3))
    dl->kern1 = get_3d_farray(dl->nkern,dl->kx,dl->ky);

  return dl;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_DCN_PREP_LAYTYPE_RELU                        */
/*                                                                           */
/*****************************************************************************/
struct dcn_layer_struct *mod_dcn_prep_laytype_relu(m,lay)
     struct model_struct *m;    // Model params
     struct onode *lay;
{
  struct dcn_layer_struct *dl;

  mylog(mylogf,"  MOD_DCN_PREP_LAYTYPE_RELU\n");

  dl = (struct dcn_layer_struct *)myalloc(sizeof(struct dcn_layer_struct));
  mod_dcn_lay_common(dl,lay);  // Set name, type, bottom and top

  mod_dcn_lay_set_size(dl,lay);  // Set sizes of inputs and outputs

  printf("    %s (layer %d)\n",dl->name,mod_dcn_lay_n);
  printf("    Output:   %d (%d,%d)\n",dl->out_zn,dl->out_xn,dl->out_yn);

  return dl;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_DCN_PREP_LAYTYPE_POOL                        */
/*                                                                           */
/*****************************************************************************/
struct dcn_layer_struct *mod_dcn_prep_laytype_pool(m,lay)
     struct model_struct *m;    // Model params
     struct onode *lay;
{
  struct onode *opar;
  struct dcn_layer_struct *dl;

  mylog(mylogf,"  MOD_DCN_PREP_LAYTYPE_POOL\n");

  dl = (struct dcn_layer_struct *)myalloc(sizeof(struct dcn_layer_struct));
  mod_dcn_lay_common(dl,lay);  // Set name, type, bottom and top

  opar = onode_child_get_unique(lay,"param");
  if (opar == NULL)
    exit_error("MOD_DCN_PREP_LAYTYPE_CONV","Cannot find <param>.");

  dl->kx     = onode_getpar_int_exit(opar,"kernel_size");
  // WYETH - NOTE 'dl->ky' is not set
  dl->stride = onode_getpar_int_dflt(opar,"stride",1);
  dl->pool   = onode_getpar_chr_exit(opar,"pool");

  mod_dcn_lay_set_size(dl,lay);  // Set sizes of inputs and outputs

  printf("    %s (layer %d)\n",dl->name,mod_dcn_lay_n);
  printf("    Stride:   %d\n",dl->stride);
  printf("    Padding:  %d\n",dl->pad);
  printf("    Output:   %d (%d,%d)\n",dl->out_zn,dl->out_xn,dl->out_yn);

  if (mod_dcn_vistype != NULL){
    //
    //  Create 'switch' arrays to save coordinates of max
    //
    dl->vis_pool_x = get_3d_iarray(dl->out_zn,dl->out_xn,dl->out_yn);
    dl->vis_pool_y = get_3d_iarray(dl->out_zn,dl->out_xn,dl->out_yn);
  }

  return dl;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_DCN_PREP_LAYTYPE_IPROD                        */
/*                                                                           */
/*****************************************************************************/
struct dcn_layer_struct *mod_dcn_prep_laytype_iprod(m,lay)
     struct model_struct *m;    // Model params
     struct onode *lay;
{
  struct onode *opar;
  struct dcn_layer_struct *dl;

  mylog(mylogf,"  MOD_DCN_PREP_LAYTYPE_IPROD\n");

  dl = (struct dcn_layer_struct *)myalloc(sizeof(struct dcn_layer_struct));
  mod_dcn_lay_common(dl,lay);  // Set name, type, bottom and top

  opar = onode_child_get_unique(lay,"param");
  if (opar == NULL)
    exit_error("MOD_DCN_PREP_LAYTYPE_IPROD","Cannot find <param>.");

  dl->nkern     = onode_getpar_int_exit(opar,"num_output");

  mod_dcn_lay_set_size(dl,lay);  // Set sizes of inputs and outputs

  printf("    %s (layer %d)\n",dl->name,mod_dcn_lay_n);
  printf("    Kernels:  %d (%d,%d,%d)\n",dl->nkern,dl->kx,dl->ky,dl->kz);
  printf("    Output:   %d (%d,%d)\n",dl->out_zn,dl->out_xn,dl->out_yn);

  dl->kern = get_4d_farray(1,dl->kz,dl->kx,dl->ky); // [1][1][kx = nkern][ky]
  dl->bias = get_farray(dl->kx);

  return dl;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_DCN_PREP_LAYTYPE_LRN                        */
/*                                                                           */
/*****************************************************************************/
struct dcn_layer_struct *mod_dcn_prep_laytype_lrn(m,lay)
     struct model_struct *m;    // Model params
     struct onode *lay;
{
  struct onode *opar;
  struct dcn_layer_struct *dl;

  mylog(mylogf,"  MOD_DCN_PREP_LAYTYPE_LRN\n");

  dl = (struct dcn_layer_struct *)myalloc(sizeof(struct dcn_layer_struct));
  mod_dcn_lay_common(dl,lay);  // Set name, type, bottom and top

  opar = onode_child_get_unique(lay,"param");
  if (opar == NULL)
    exit_error("MOD_DCN_PREP_LAYTYPE_LRN","Cannot find <param>.");

  dl->loc_size  = onode_getpar_int_exit(opar,"local_size");
  dl->alpha     = onode_getpar_flt_exit(opar,"alpha");
  dl->beta      = onode_getpar_flt_exit(opar,"beta");

  if ((dl->loc_size % 2) != 1)
    exit_error("MOD_DCN_PREP_LAYTYPE_LRN","Expecting 'loc_size' to be odd");

  mod_dcn_lay_set_size(dl,lay);  // Set sizes of inputs and outputs

  printf("    %s (layer %d)\n",dl->name,mod_dcn_lay_n);
  printf("    Output:   %d (%d,%d)\n",dl->out_zn,dl->out_xn,dl->out_yn);

  return dl;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_DCN_PREP_LAYTYPE_SFTMX                      */
/*                                                                           */
/*****************************************************************************/
struct dcn_layer_struct *mod_dcn_prep_laytype_sftmx(m,lay)
     struct model_struct *m;    // Model params
     struct onode *lay;
{
  struct dcn_layer_struct *dl;

  mylog(mylogf,"  MOD_DCN_PREP_LAYTYPE_SFTMX\n");

  dl = (struct dcn_layer_struct *)myalloc(sizeof(struct dcn_layer_struct));
  mod_dcn_lay_common(dl,lay);  // Set name, type, bottom and top

  mod_dcn_lay_set_size(dl,lay);  // Set sizes of inputs and outputs

  printf("    %s (layer %d)\n",dl->name,mod_dcn_lay_n);
  printf("    Output:   %d (%d,%d)\n",dl->out_zn,dl->out_xn,dl->out_yn);

  return dl;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_DCN_01_PREP                              */
/*                                                                           */
/*  Prepare constructs that will be needed when running trials.  Do things   */
/*  here that only need to be done once, to avoid repeating them for each    */
/*  trial.                                                                   */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_01_prep(m,s,r)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
{
  int i,k;
  int n,xn,yn,tn,kmax;
  float tscale,sscale;
  char *laytype,*wfname,*stimfile;
  struct onode *fito,*s1o,*c1o,*s2o,*c2o,*t,*wmovr;
  void mod_dcn_01_get_response_override();

  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  mylog(mylogf,"  MOD_DCN_01_PREP\n");

  xn     = mod_dcn_xn     = onode_getpar_int_exit(m->o,"xn");
  yn     = mod_dcn_yn     = onode_getpar_int_exit(m->o,"yn");
  tn     = mod_dcn_tn     = onode_getpar_int_exit(m->o,"tn");
  tscale = mod_dcn_tscale = onode_getpar_flt_exit(m->o,"tscale");
  sscale = mod_dcn_sscale = onode_getpar_flt_exit(m->o,"sscale");
  mod_dcn_zn              = onode_getpar_int_exit(m->o,"stim_depth");
  mod_dcn_vistype         = onode_getpar_chr_dflt(m->o,"visualization",NULL);

  //
  //  Configure visualization.
  //
  if (mod_dcn_vistype != NULL){
    if ((strcmp(mod_dcn_vistype,"NULL")==0)||strcmp(mod_dcn_vistype,"null")==0){
      myfree(mod_dcn_vistype);
      mod_dcn_vistype = NULL;  // This indicates no visualization
    }else{
      mod_dcn_vismode       = onode_getpar_chr_exit(m->o,"vis_mode");
      mod_dcn_vislay        = onode_getpar_chr_exit(m->o,"vis_layer");
      mod_dcn_visopre       = onode_getpar_chr_exit(m->o,"vis_outfile");
      mod_dcn_viso3d        = onode_getpar_int_dflt(m->o,"vis_out_3d",0);

      if (strcmp(mod_dcn_vismode,"unit")==0){
	mod_dcn_visu_xi = onode_getpar_int_dflt(m->o,"vis_unit_xi");
	mod_dcn_visu_yi = onode_getpar_int_dflt(m->o,"vis_unit_yi");
	mod_dcn_visu_zi = onode_getpar_int_dflt(m->o,"vis_unit_zi");
      }
    
      sprintf(ggstr,"    Visualization turned on:  %s\n",mod_dcn_vistype);
      mylog(mylogf,ggstr);
      sprintf(ggstr,"    Visualization mode:  %s\n",mod_dcn_vismode);
      mylog(mylogf,ggstr);
      sprintf(ggstr,"    Visualization layer:  %s\n",mod_dcn_vislay);
      mylog(mylogf,ggstr);
      if (strcmp(mod_dcn_vismode,"unit")==0){
	sprintf(ggstr,"    Visualization unit:  %d %d %d\n",mod_dcn_visu_zi,
		mod_dcn_visu_xi,mod_dcn_visu_yi);
	mylog(mylogf,ggstr);
      }

      // Storage for reconstructed visualization
      mod_dcn_vis_rec = get_3d_farray(3,xn,yn);
    }
  }

  //
  //  PROCESS EACH LAYER
  //
  n = onode_count_otype(m->o,"layer");
  printf("    %d layers\n",n);

  mod_dcn_lay = (struct dcn_layer_struct **)
    myalloc(n*sizeof(struct dcn_layer_struct *));

  mod_dcn_rsize = 0;  // Count number of floats for respones storage
  mod_dcn_lay_n = 0;  // This must increase as layers are built
  k = 0;
  t = onode_get_next_type(m->o->o,"layer");
  while(t != NULL){
    //
    //  Loop across all layers, preparing each one in sequence
    //
    laytype = onode_getpar_chr_exit(t,"type");

    if (strcmp(laytype,"Convolution")==0){
      mod_dcn_lay[k] = mod_dcn_prep_laytype_conv(m,t);
    }else if (strcmp(laytype,"ReLU")==0){
      mod_dcn_lay[k] = mod_dcn_prep_laytype_relu(m,t);
    }else if (strcmp(laytype,"Pooling")==0){
      mod_dcn_lay[k] = mod_dcn_prep_laytype_pool(m,t);
    }else if (strcmp(laytype,"InnerProduct")==0){
      mod_dcn_lay[k] = mod_dcn_prep_laytype_iprod(m,t);
    }else if (strcmp(laytype,"LRN")==0){
      mod_dcn_lay[k] = mod_dcn_prep_laytype_lrn(m,t);
    }else if (strcmp(laytype,"Softmax")==0){
      mod_dcn_lay[k] = mod_dcn_prep_laytype_sftmx(m,t);
    }else{
      printf("  laytype = %s\n",laytype);
      exit_error("MOD_DCN_01_PREP","Unknown layer type.");
    }

    k += 1;
    mod_dcn_lay_n = k;  // This must increase as layers are built
    t = onode_get_next_type(t->next,"layer");
  }

  //
  //  Determine deepest level of response requested
  //
  kmax = -1;
  for(i=0;i<r->n;i++){ // For each response requested
    if ((r->rformat[i] == 1) || // 'save_pop_unit_as_'
	(r->rformat[i] == 2)){  // 'save_pop_grid_as_'
      k = mod_dcn_layi_for_name(r->plname[i]);
      if (k < 0){
	printf("  Layer name:  %s\n",r->plname[i]);
	exit_error("MOD_DCN_01_PREP","Layer name not found");
      }
      if (k > kmax)
	kmax = k;
    }
  }
  // Check visualization layer
  if (mod_dcn_vistype != NULL){
    k = mod_dcn_layi_for_name(mod_dcn_vislay);
    if (k < 0){
      printf("  Layer name:  %s\n",r->plname[i]);
      exit_error("MOD_DCN_01_PREP","Visualizaation layer not found");
    }
    if (k > kmax)
      kmax = k;
  }

  mod_dcn_maxli = kmax;
  printf("  Deepest level required for response request:  %d\n",kmax);

  mod_dcn_prep_response();

  //
  //  Read in the weights and biases for all layers.
  //
  wfname = onode_getpar_chr_exit(m->o,"weight_filename");
  mod_dcn_read_weights(wfname);


  // Compute statistics for kernels
  /*  WYETH - THIS SHOULD NOT RUN EXCEPT IN DEBUG OR 'wm' MODE, not in 'mm'
  for(i=0;i<mod_dcn_maxli;i++){
    if (strcmp(mod_dcn_lay[i]->type,"Convolution")==0){
      mod_dcn_stat_kern(mod_dcn_lay[i]);
    }
  }
  */


  //
  //  If kernel images are requested for checking, dump kernels and exit.
  //
  {
    int xflag;
    char *outfile;
    char *lname;
    char *tstr;

    xflag = 0;  // Default - do not exit

    outfile  = onode_getpar_chr_dflt(m->o,"write_kernels_rgb",NULL);
    if (outfile != NULL){
      mod_dcn_write_rgb_kernels(outfile);
      xflag = 1;
    }

    t = onode_get_next_type(m->o->o,"item");
    while(t != NULL){
      if (strcmp(t->name,"write_kernel_conv")==0){
	printf("  Name:  %s\n",t->name);

	lname   = onode_get_nth_val(t,0);
	tstr    = onode_get_nth_val(t,1);
	outfile = onode_get_nth_val(t,2);

	mod_dcn_write_kernel(outfile,lname,atoi(tstr));

	myfree(lname);
	myfree(tstr);
	myfree(outfile);
	xflag = 1;
      }
      t = onode_get_next_type(t->next,"item");
    }

    if (xflag == 1){
      printf("  Exiting after writing output files for checking.\n\n");
      exit(0);
    }
  }

  //
  //  Determine stimulus type.  If this is "3rgb" we have to do special
  //  processing.
  //
  mod_dcn_stimform = param_getc_dflt(s->ppl,"stim_form","3d");
  if (strcmp(mod_dcn_stimform,"3rgb")==0){
    mod_dcn_tstim = get_3d_farray(3,xn,yn); // Allocate storage for RGB image
  }else{
    mod_dcn_tstim = NULL;
  }

  mod_dcn_resp_3d_fname = onode_getpar_chr_dflt(m->o,"write_3d_resp",NULL);


  //
  //  Check for 'wm_override' object
  //
  wmovr = onode_child_get_unique(m->o,"wm_override");
  if (wmovr != NULL){
    mod_dcn_wmover        = onode_getpar_chr_exit(wmovr,"mode");
    if ((strcmp(mod_dcn_wmover,"NULL")==0)||
	(strcmp(mod_dcn_wmover,"null")==0)){
      myfree(mod_dcn_wmover);
      mod_dcn_wmover = NULL;
    }else{
      mod_dcn_wmover_flist  = onode_getpar_chr_exit(wmovr,"stim_list");
      mod_dcn_wmover_fpath  = onode_getpar_chr_exit(wmovr,"stim_path");
      mod_dcn_wmover_ofile  = onode_getpar_chr_exit(wmovr,"stat_file");
      mod_dcn_01_get_response_override(m,s,r);
    }
  }else{
    mod_dcn_wmover = NULL;
  }

  // WYETH HERE TESTING a way to visualize drive from FC6 units to FC7
  //   WYETH - never finished this
  // mod_dcn_map_fc("zz.zz","fc6",23);

}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_DCN_01_DONE                              */
/*                                                                           */
/*  Free any storage here when done.                                         */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_01_done(m)
     struct model_struct *m;        // Model parameters
{
  mylog(mylogf,"  MOD_DCN_01_DONE\n");
  
  // Nothing to do here.
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_DCN_GET_PADDED_R                           */
/*                                                                           */
/*****************************************************************************/
float ***mod_dcn_get_padded_r(layin,pad)
     struct dcn_layer_struct *layin;   // Pointer to layer
     int pad;
{
  int i,j,k;
  int xn,zn,xnp,xn2p;
  float ***d,**t,**tin;

  zn   = layin->out_zn;
  xn   = layin->out_xn;
  if (layin->out_xn != layin->out_yn)
    exit_error("MOD_DCN_GET_PADDED_R","xn != yn");
  xnp  = xn + pad;
  xn2p = xn + 2*pad;

  d = get_3d_farray(zn,xn2p,xn2p);
  //printf("Making 3d Farray %d %d %d\n",zn,xn2p,xn2p);

  // Fill padded margins with 0.0
  for(i=0;i<zn;i++){
    t = d[i];
    for(j=0;j<pad;j++){
      for(k=0;k<xnp;k++){
	t[k][j]         = 0.0;  // bottom
	t[j][k+pad]     = 0.0;  // left
	t[k+pad][j+xnp] = 0.0;  // top
	t[j+xnp][k]     = 0.0;  // right
      }
    }
  }

  // Fill middle with data
  for(i=0;i<zn;i++){

    t = d[i];
    tin = layin->r[i];

    for(j=0;j<xn;j++){
      for(k=0;k<xn;k++){
	t[pad+j][pad+k] = tin[j][k];
      }
    }
  }

  return d;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_DCN_PUT_PADDED_R                           */
/*                                                                           */
/*  Write a padded response array back into the smaller original array.      */
/*                                                                           */
/*  This is used for VISUALIZATION processing only.                          */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_put_padded_r(lay,pad,d)
     struct dcn_layer_struct *lay;   // Pointer to layer
     int pad;
     float ***d;     // Data with padding
{
  int i,j,k;
  int xn,zn;
  float **t,**tin;

  zn   = lay->out_zn;
  xn   = lay->out_xn;
  if (lay->out_xn != lay->out_yn)
    exit_error("MOD_DCN_PUT_PADDED_R","xn != yn");

  // Fill target response array with middle of padded data 'd'
  for(i=0;i<zn;i++){

    t = d[i];
    tin = lay->r[i];

    for(j=0;j<xn;j++){
      for(k=0;k<xn;k++){
	tin[j][k] = t[pad+j][pad+k];
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_DCN_VIS_PATCH_GET                           */
/*                                                                           */
/*  Extract a patch from an image in [3][xn][yn] format, padding out of      */
/*  bounds regions.                                                          */
/*                                                                           */
/*****************************************************************************/
float ***mod_dcn_vis_patch_get(data,xn,yn,x0,y0,w,h,padval)
     float ***data;  // [3][xn][yn] image data
     int xn,yn;      // size of image
     int x0,y0;      // lower left corner of patch.
     int w,h;        // width and height of patch
     float padval;   // 0..1
{
  int i,j,k;
  int xi,yi;
  float ***td;

  td = get_3d_farray(3,w,h);
  for(i=0;i<w;i++){
    for(j=0;j<h;j++){
      xi = x0 + i;
      yi = y0 + j;
      if ((xi >= 0) && (xi < xn) && (yi >= 0) && (yi < yn)){
	for(k=0;k<3;k++)
	  td[k][i][j] = data[k][xi][yi];
      }else{
	for(k=0;k<3;k++)
	  td[k][i][j] = padval;  // Assume data is already in 0..1
      }
    }
  }

  return td;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_DCN_VIS_PATCH_BOX                           */
/*                                                                           */
/*  Draw a box on the image.                                                 */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_vis_patch_box(d,xn,yn,x0,y0,w,h,r,g,b)
     float ***d;  // [3][xn][yn]
     int xn,yn;   // image size
     int x0,y0;   // box origin
     int w,h;     // box width, height
     float r,g,b; // color of line
{
  int i,j;

  //
  //  Draw a red box around the patch
  //
  for(i=x0;i<(x0+w);i++){
    j = y0;
    if ((i>=0)&&(i<xn)&&(j>=0)&&(j<yn)){
      d[0][i][j] = b;
      d[1][i][j] = g;
      d[2][i][j] = r;
    }
    j = y0+h-1;
    if ((i>=0)&&(i<xn)&&(j>=0)&&(j<yn)){
      d[0][i][j] = b;
      d[1][i][j] = g;
      d[2][i][j] = r;
    }
  }
  for(j=y0;j<(y0+h);j++){
    i = x0;
    if ((i>=0)&&(i<xn)&&(j>=0)&&(j<yn)){
      d[0][i][j] = b;
      d[1][i][j] = g;
      d[2][i][j] = r;
    }
    i = x0+w-1;
    if ((i>=0)&&(i<xn)&&(j>=0)&&(j<yn)){
      d[0][i][j] = b;
      d[1][i][j] = g;
      d[2][i][j] = r;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_DCN_VIS_PATCH_WRITE                         */
/*                                                                           */
/*  Write the patch of the recon image corresponding to the unit RF.         */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_vis_patch_write(lay,xi,yi)
     struct dcn_layer_struct *lay;     // Layer
     int xi,yi;                        // x,y coordinates of unit
{
  int i,j,k;
  int sz,x0,dx,x,y,xx,yy,xn,yn;
  float ***td;
  char fname[SLEN];

  printf("  MOD_DCN_VIS_PATCH_WRITE\n");

  xn = mod_dcn_xn;
  yn = mod_dcn_yn;

  if (strcmp(lay->name,"conv1")==0){
    //
    //  11 + (55-1)*4 = 227 pix   Input
    //  19 + (27-1)*8 = 227 pix   Output covers input image
    //
    x0 = 0;
    dx = 4;
    sz = 11;
  }else if (strcmp(lay->name,"conv2")==0){  // pad 2
    //
    //  51 + (27-1)*8 = 259 [Thus, 32 pixels wider than the image]
    //  Thus, 16 pix on either end
    //
    //  67 + (13-1)*16 = 259       Output matches input in pix
    //
    x0 = -16;
    dx = 8;
    sz = 51;
  }else if (strcmp(lay->name,"conv3")==0){  // pad 1
    //
    //  99 + (13-1)*16 = 291 pix  [Thus, 64 pixels wider than image]
    //
    x0 = -32;
    dx = 16;
    sz = 99;
  }else if (strcmp(lay->name,"conv4")==0){  // pad 1
    //
    //  131 + (13-1)*16 = 323 pix  [Thus, 96 pixels wider than image]
    //
    x0 = -48;
    dx = 16;
    sz = 131;
  }else if (strcmp(lay->name,"conv5")==0){  // pad 1
    //
    //  163 + (13-1)*16 = 355 pix  [Thus, 128 pixels wider than image]
    //
    x0 = -64;
    dx = 16;
    sz = 163;
  }else
    exit_error("MOD_DCN_VIS_PATCH_WRITE","Layer name error");

  x = x0 + xi*dx;
  y = x0 + yi*dx;  // Use 'dx' for 'dy' because everything is square

  printf("    Layer:   %s\n",lay->name);
  printf("    Origin:  %d\n",x0);
  printf("    Patch size:  %d\n",sz);
  printf("    Patch x,y:   %d %d\n",x,y);

  td =  mod_dcn_vis_patch_get(mod_dcn_vis_rec,xn,yn,x,y,sz,sz,1.0);
  sprintf(fname,"%s.e.sub.ppm",mod_dcn_visopre);
  write_ppm_6_rgb_data(fname,sz,sz,td[2],td[1],td[0],1);
  free_3d_farray(td,3,sz,sz);

  if (mod_dcn_tstim != NULL){
    td =  mod_dcn_vis_patch_get(mod_dcn_tstim,xn,yn,x,y,sz,sz,1.0);
    sprintf(fname,"%s.d.so.ppm",mod_dcn_visopre);
    write_ppm_6_rgb_data(fname,sz,sz,td[2],td[1],td[0],1);
    free_3d_farray(td,3,sz,sz);

    //x -= 1;
    //y -= 1;
    //sz += 2;
    //
    //  Draw a double line red box around the patch
    //
    mod_dcn_vis_patch_box(mod_dcn_tstim,xn,yn,x-1,y-1,sz+2,sz+2,1.0,0.0,0.0);
    mod_dcn_vis_patch_box(mod_dcn_tstim,xn,yn,x-2,y-2,sz+4,sz+4,1.0,0.0,0.0);

    sprintf(fname,"%s.b.box.ppm",mod_dcn_visopre);
    write_ppm_6_rgb_data(fname,xn,yn,mod_dcn_tstim[2],
			 mod_dcn_tstim[1],mod_dcn_tstim[0],1);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_DCN_VIS_CONV                             */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_vis_conv(dl)
     struct dcn_layer_struct *dl;     // Pointer to layer
{
  int i,j;
  int ki,kx,ky,kz,nk,nx,ny,nz,ndp_x,ndp_y,stride,x0,y0,z0,pad,grp_sz,inz,cnt;
  int ngr;
  float ***kern,**kp1,***d,tot,**s,bias,rr;
  struct dcn_layer_struct *layin;

  //mylog(mylogf,"  MOD_DCN_VIS_CONV\n");

  ndp_x  = dl->out_xn;  // Number of dot products, x-axis
  ndp_y  = dl->out_yn;  // Number of dot products, y-axis
  nx     = dl->kx;      // Size of kernel, x-axis
  ny     = dl->ky;      // Size of kernel, y-axis
  nz     = dl->kz;      // Size of kernel, z-axis
  nk     = dl->nkern;   // Number of kernels
  stride = dl->stride;  // Stride (pixels between dot products)
  pad    = dl->pad;     // Padding

  //
  //  Make 'd' point to the "deconv" target, and zero it
  //
  if (strcmp(dl->bot,"data")==0){  // Input is the global test image
    d = mod_dcn_vis_rec;
    inz = mod_dcn_zn;
    zero_3d_farray(d,inz,mod_dcn_xn,mod_dcn_yn);
  }else{

    layin = mod_dcn_lay_for_name(dl->bot);

    if (pad == 0){
      d = layin->r;
      zero_3d_farray(d,layin->out_zn,layin->out_xn,layin->out_yn);
    }else{
      d = get_zero_3d_farray(layin->out_zn,layin->out_xn+2*pad,
			                   layin->out_yn+2*pad);
      //d = mod_dcn_get_padded_r(layin,pad);
    }
    inz = layin->out_zn;
  }

  //printf("NUMBER______OF_____KERN:  %d   nz= %d\n",nk,nz);
  ngr = inz/nz;  // Number of goups
  if (ngr != dl->group)
    exit_error("MOD_DCN_GETRESP_CONV","Group size error");

  grp_sz = dl->nkern/ngr;

  cnt = 0;  // Counter for determining if/when to start new group
  z0 = 0;   // Offset only changes for subsequent groups, if any

  for(ki=0;ki<nk;ki++){  // For each unique kernel

    kern = dl->kern[ki];  // Get a pointer to this kernel
    bias = dl->bias[ki];  // Bias value

    if (cnt >= grp_sz){  // If we have finished the first group
      z0 += nz;  // increment the kernel offset
      cnt = 0;
    }
    //printf("ki %d  z0 %d\n",ki,z0);

    x0 = 0;  // Initialize x- and y-offsets into the input data 'd'
    y0 = 0;

    for(i=0;i<ndp_x;i++){    // For each location to compute dot product,
      for(j=0;j<ndp_y;j++){  //   multiply kernel 'kern' by input data 'd'.

	// We already have the response
	rr = dl->r[ki][i][j];  // WYETH, ?? Subtract the bias ???

	for(kx=0;kx<nx;kx++){
	  for(ky=0;ky<ny;ky++){
	    for(kz=0;kz<nz;kz++){
	      // Add the kernel weight back into the input
	      // WYETH WYETH - MUST ZERO THE INPUT FIRST
	      //tot += kern[kz][kx][ky] * d[z0+kz][x0+kx][y0+ky];
	      d[z0+kz][x0+kx][y0+ky] += rr * kern[kz][kx][ky];
	    }
	  }
	}
	//dl->r[ki][i][j] = tot + bias;  //  Store response in layer struct

	y0 += stride;  // Step the y-offset within 'd'
      }
      x0 += stride;  // Step the x-offset within 'd'
      y0 = 0;        // Reset the y-offset
    }

    cnt += 1;
  }

  if (pad > 0){
    //
    //  Write the valid data from the padded array back into the input
    //
    mod_dcn_put_padded_r(layin,pad,d);

    free_3d_farray(d,nz,nx+2*pad,nx+2*pad);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_DCN_VIS_IPROD                             */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_vis_iprod(dl)
     struct dcn_layer_struct *dl;     // Pointer to layer
{
  int i,j,k;
  int ki,inx,iny,inz,kk,nkern,tn;
  float *kdata,*t,rr;
  struct dcn_layer_struct *layin;

  nkern = dl->kx;

  layin = mod_dcn_lay_for_name(dl->bot);
  inx = layin->out_xn;
  iny = layin->out_yn;
  inz = layin->out_zn;

  tn = inx * iny * inz;

  if (dl->ky != tn)
    exit_error("MOD_DCN_VIS_IPROD","Kernel size mismatch");

  //printf("    Input layer is:  %d x %d x %d\n",inx,iny,inz);
  //printf("    Number of kernels:  %d\n",nkern);
  //printf("      Length of kernels:  %d\n",dl->ky);

  if ((inz == 1) && (inx == 1) && (iny == dl->ky)){
    //
    //  Previous layer is also FC ???
    //
    t = layin->r[0][0];
    zero_farray(t,tn);
  }else{
    //
    //  Arrange responses from previous layer in a single array
    //

    t = get_zero_farray(tn);
    /*
    kk = 0;
    for(k=0;k<inz;k++){
      for(j=(iny-1);j>=0;j--){
	for(i=0;i<inx;i++){
	  t[kk] = layin->r[k][i][j];
	  kk += 1;
	}
      }
    }*/
  }

  for(ki=0;ki<nkern;ki++){  // For each kernel

    kdata = dl->kern[0][0][ki];  // [1][1][kx][ky]

    rr = dl->r[0][0][ki];  // Response for this kernel

    if (rr != 0.0){

      for(i=0;i<tn;i++)
	t[i] += rr * kdata[i];
        //tot += t[i] * kdata[i];
    }

    //dl->r[0][0][ki] = tot + dl->bias[ki];
  }

  if ((inz == 1) && (inx == 1) && (iny == dl->ky)){
    ; // No-op
  }else{
    //
    //  Put responses back into previous layer
    //
    kk = 0;
    for(k=0;k<inz;k++){
      for(j=(iny-1);j>=0;j--){
	for(i=0;i<inx;i++){
	  //t[kk] = layin->r[k][i][j];
	  layin->r[k][i][j] = t[kk];
	  kk += 1;
	}
      }
    }
    myfree(t);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_DCN_VIS_LRN                              */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_vis_lrn(dl)
     struct dcn_layer_struct *dl;   // Pointer to layer
{
  int i,j,k;
  int xn,yn,zn;
  float ***d,***r;
  struct dcn_layer_struct *layin;

  xn = dl->out_xn;
  yn = dl->out_yn;
  zn = dl->out_zn;

  layin = mod_dcn_lay_for_name(dl->bot);

  d = layin->r;
  r = dl->r;

  for(i=0;i<zn;i++){
    for(j=0;j<xn;j++){
      for(k=0;k<yn;k++){
	d[i][j][k] = r[i][j][k];  // Move data to the output of previous layer
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_DCN_VIS_POOL                             */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_vis_pool(dl)
     struct dcn_layer_struct *dl;   // Pointer to layer
{
  int i,j,k;
  int xn,yn,zn,mxi,myi;
  float ***d,***r,**t;
  struct dcn_layer_struct *layin;

  if (strcmp(dl->pool,"MAX")!=0)
    exit_error("MOD_DCN_VIS_POOL","Pooling type must be 'MAX'");

  xn = dl->out_xn;
  yn = dl->out_yn;
  zn = dl->out_zn;

  layin = mod_dcn_lay_for_name(dl->bot);

  // Zero all responses, because only some will be over-written below
  zero_3d_farray(layin->r,layin->out_zn,layin->out_xn,layin->out_yn);

  d = layin->r;
  r = dl->r;

  for(i=0;i<zn;i++){
    t = d[i];  // Pointer to 2D response plane

    for(j=0;j<xn;j++){
      for(k=0;k<yn;k++){  // For each output value

	mxi = dl->vis_pool_x[i][j][k];
	myi = dl->vis_pool_y[i][j][k];

	//if (t[mxi][myi] != 0.0){

	if (1){
	  //
	  //  New way - do not overwrite larger values
	  //
	  if (r[i][j][k] > t[mxi][myi])
	    t[mxi][myi] = r[i][j][k];  // Put response back at location of max.
	  //else if (r[i][j][k] < t[mxi][myi])
	  //printf("Not overwriting %f  with  %f\n",t[mxi][myi],r[i][j][k]);
	  //
	  //  NOTE, no need to write anything negative, because the previous
	  //    relu layer will set these to zero anyway.
	}else{
	  //
	  //  OLD WAY - can overwrite values
	  //
	  t[mxi][myi] = r[i][j][k];  // Put response back at location of max.
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_DCN_VIS_RELU                             */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_vis_relu(dl)
     struct dcn_layer_struct *dl;   // Pointer to layer
{
  int i,j,k;
  int xn,yn,zn;
  float ***d,***r;
  struct dcn_layer_struct *layin;

  xn = dl->out_xn;
  yn = dl->out_yn;
  zn = dl->out_zn;

  layin = mod_dcn_lay_for_name(dl->bot);

  if ((layin->out_xn != xn) || (layin->out_yn != yn) || (layin->out_zn != zn))
    exit_error("MOD_DCN_VIS_RELU","Input / Output size mismatch");

  d = layin->r;
  r = dl->r;

  for(i=0;i<zn;i++){
    for(j=0;j<xn;j++){
      for(k=0;k<yn;k++){
	if (r[i][j][k] < 0.0)
	  d[i][j][k] = 0.0;
	else
	  d[i][j][k] = r[i][j][k];  // Put response back in previous layer
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_DCN_VISUALIZE_01                           */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_visualize_01(vl)
     struct dcn_layer_struct *vl;   // Visualization target layer
{
  int i,j,k;
  int vli,zn,xn,yn,zi,xi,yi,txn,tyn,tzn;
  float mu,sd,rmax,***d;
  char fname[SLEN];
  struct dcn_layer_struct *dl;

  mylog(mylogf,"  MOD_DCN_VISUALIZE_01\n");

  //
  // *** WYETH - KNOWN ISSUES
  //  (1) I have only tested this with "conv" layers.  It may not work with
  //      other layer types.
  //

  xn = vl->out_xn;
  yn = vl->out_yn;
  zn = vl->out_zn;

  if (strcmp(mod_dcn_vismode,"max")==0){
    //
    //  Find the max response in vl->r[zn][xn][yn]
    //
    get_max_coord_3d_farray(vl->r,0,zn,0,xn,0,yn,&zi,&xi,&yi);
    rmax = vl->r[zi][xi][yi];
    printf("    Found max at:  %d %d %d  (z,x,y) within %d %d %d\n",zi,xi,yi,
	   zn,xn,yn);
    printf("          max =  %f\n",rmax);

    single_mean_sdev_3d_farray(vl->r,0,zn,0,xn,0,yn,&mu,&sd);
    printf("          mean = %f (SD %f)\n",mu,sd);
  }else if (strcmp(mod_dcn_vismode,"unit")==0){
    //
    //  Use the specified unit
    //
    zi = mod_dcn_visu_zi;
    xi = mod_dcn_visu_xi;
    yi = mod_dcn_visu_yi;
    rmax = vl->r[zi][xi][yi];
    printf("    Target unit is: %s %d %d %d (z,x,y) within %d %d %d\n",
	   mod_dcn_vislay,zi,xi,yi,zn,xn,yn);
    printf("      response value:  %f\n",rmax);
  }else{
    exit_error("MOD_DCN_VISUALIZE_01","Unknown 'vis_mode' value.");
  }

  //
  //  Zero all responses, except the max
  //
  zero_3d_farray(vl->r,zn,xn,yn);
  vl->r[zi][xi][yi] = rmax;

  vl->vxi = xi;  // Save the index of the unit to visualize
  vl->vyi = yi;
  vl->vzi = zi;

  //
  //  Determine the layer index
  //
  vli = mod_dcn_layi_for_name(vl->name);

  //
  //  Work back down to the pixel level.
  //
  for(i=vli;i>=0;i--){  // Only compute layers needed for response

    dl = mod_dcn_lay[i];

    //
    //  DEBUGGING:  Write out the 3D response
    //
    //sprintf(fname,"zzz.%d.%s.3d",i,dl->name);
    //mod_dcn_write_3d_resp(fname,dl);

    sprintf(ggstr,"    %2d %s\n",i,dl->name);
    mylog(mylogf,ggstr);

    if (strcmp(dl->type,"Convolution")==0){
      mod_dcn_vis_conv(dl);
    }else if (strcmp(dl->type,"LRN")==0){
      mod_dcn_vis_lrn(dl);
    }else if (strcmp(dl->type,"Pooling")==0){
      mod_dcn_vis_pool(dl);
    }else if (strcmp(dl->type,"ReLU")==0){
      mod_dcn_vis_relu(dl);
    }else if (strcmp(dl->type,"InnerProduct")==0){
      mod_dcn_vis_iprod(dl);
    }else{
      printf("  Layer type:  %s\n",dl->type);
      exit_error("MOD_DCN_VISUALIZE_01","Unknown layer type.");
    }

    /*
    dl = mod_dcn_lay[i-1];
    printf("===> name: %s\n",dl->name);
    printf("===> new val:  %f\n",dl->r[149][13][10]);
    exit(0);
    */
  }

  //
  //  Write the reconstructed image to an output file.
  //
  txn = mod_dcn_xn;
  tyn = mod_dcn_yn;
  tzn = mod_dcn_zn;

  if (mod_dcn_viso3d == 1){
    //
    //  Write r,g,b planes separately to a .3d file.
    //
    //  Copy into x,y,t order
    d = get_3d_farray(txn,tyn,tzn);
    for(i=0;i<tzn;i++)
      for(j=0;j<txn;j++)
	for(k=0;k<tyn;k++)
	  d[j][k][i] = mod_dcn_vis_rec[i][j][k];

    sprintf(fname,"%s.3d",mod_dcn_visopre);
    write_3d_data_part(fname,d,0,txn,0,tyn,0,tzn,4,2,1); // 1=txy format
    printf("    Wrote r,g,b planes separately to:  %s\n",fname);
  }

  //
  //  Write a color .ppm image
  //
  norm_01_3d_farray(mod_dcn_vis_rec,tzn,txn,tyn,1);

  sprintf(fname,"%s.c.ppm",mod_dcn_visopre);
  //write_ppm_6_rgb_data(fname,txn,tyn,mod_dcn_vis_rec[0],
  //		       mod_dcn_vis_rec[1],mod_dcn_vis_rec[2],1);
  write_ppm_6_rgb_data(fname,txn,tyn,mod_dcn_vis_rec[2],
		       mod_dcn_vis_rec[1],mod_dcn_vis_rec[0],1);
  printf("    Wrote RGB reconstruction image to:  %s\n",fname);


  //
  //  Write original stimulus
  //
  if (mod_dcn_tstim != NULL){
    norm_01_3d_farray(mod_dcn_tstim,tzn,txn,tyn,1);

    sprintf(fname,"%s.a.orig.ppm",mod_dcn_visopre);

    write_ppm_6_rgb_data(fname,txn,tyn,mod_dcn_tstim[2],
			 mod_dcn_tstim[1],mod_dcn_tstim[0],1);

    printf("    Wrote RGB reconstruction image to:  %s\n",fname);
  }

  //
  //  Write the segment corresponding to the patch
  //
  mod_dcn_vis_patch_write(vl,xi,yi);

}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_DCN_RESP_TSAVE                            */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_resp_tsave(r,ti)
     struct response_struct *r;  // Response params
     int ti;                     // Time index [0...tn-1]
{
  int i;
  int xi,yi,zi,tn;
  float fval;
  struct dcn_layer_struct *lay;

  tn = mod_dcn_tn;

  for(i=0;i<r->n;i++){ // For each response requested

    if ((r->rformat[i] == 1) ||  // 'save_pop_unit_as_'
	(r->rformat[i] == 2)){   // 'save_pop_grid_as_'

      //
      // Example:                                      x y z
      // save_pop_unit_as_conv1_0_0_0  f 1000  conv1   0 0 0  raw
      //

      lay = mod_dcn_lay_for_name(r->plname[i]);
      if (lay == NULL){
	printf("  Layer name:  %s\n",r->plname[i]);
	exit_error("MOD_DCN_RESP_TSAVE","Layer name not found");
      }

      xi = r->xi[i];
      yi = r->yi[i];
      zi = r->zi[i];

      fval = lay->r[zi][xi][yi];

    }else{
      printf("  Response format:  %d\n",r->rformat[i]);
      exit_error("MOD_DCN_RESP_TSAVE","Format not implemented");
    }

    if (lay->rt[zi][xi][yi] == NULL) // Allocate memory only the first time.
      lay->rt[zi][xi][yi] = get_farray(tn);

    lay->rt[zi][xi][yi][ti] = fval;  // Save this value within time sequence
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_DCN_HANDOVER                             */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_handover(r)
     struct response_struct *r; // Response params
{
  int i;
  int xi,yi,zi,tn;
  float tscale,fval,*t;
  struct dcn_layer_struct *lay;

  tscale = mod_dcn_tscale;
  tn     = mod_dcn_tn;

  for(i=0;i<r->n;i++){ // For each response requested

    if (strcmp(r->datid[i],"raw")!=0){
      printf("  Data ID:  %s\n",r->datid[i]);
      exit_error("MOD_DCN_HANDOVER","Only 'raw' response is available");
    }

    t = NULL;

    if ((r->rformat[i] == 1) ||  // 'save_pop_unit_as_'
	(r->rformat[i] == 2)){   // 'save_pop_grid_as_'

      //
      // Example:                                      x y z
      // save_pop_unit_as_conv1_0_0_0  f 1000  conv1   0 0 0  raw
      //

      lay = mod_dcn_lay_for_name(r->plname[i]);
      if (lay == NULL){
	printf("  Layer name:  %s\n",r->plname[i]);
	exit_error("MOD_DCN_HANDOVER","Layer name not found");
      }

      xi = r->xi[i];
      yi = r->yi[i];
      zi = r->zi[i];

      //printf("%s  xi,yi,zi = %d %d %d\n",r->plname[i],xi,yi,zi);

      if (tn == 1)
	fval = lay->r[zi][xi][yi];  // There is a single value
      else
	t = lay->rt[zi][xi][yi];     // There is a time sequence

    }else{
      printf("  Response format:  %d\n",r->rformat[i]);
      exit_error("MOD_DCN_HANDOVER","Format not implemented");
    }

    if (tn == 1){
      t = get_const_farray(tn,fval);
      mod_util_resp_store_f_samp(r,i,t,tn,1.0/tscale,mylogf);
      myfree(t);
    }else{
      mod_util_resp_store_f_samp(r,i,t,tn,1.0/tscale,mylogf);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_DCN_GETRESP_CONV                           */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_getresp_conv(dl,stim)
     struct dcn_layer_struct *dl;   // Pointer to layer
     float ***stim;                 // Stimulus data, or NULL
{
  int i,j;
  int ki,kx,ky,kz,nk,nx,ny,nz,ndp_x,ndp_y,stride,x0,y0,z0,pad,grp_sz,inx,inz;
  int ngr,cnt;
  float ***kern,**kp1,***d,tot,**s,bias;
  struct dcn_layer_struct *layin;

  //mylog(mylogf,"  MOD_DCN_GETRESP_CONV\n");

  // *** UNRESOLVED ISSUES ***
  //  (1) deal with padding
  //  (3) deal with ordering of responses along y-axis.

  ndp_x  = dl->out_xn;  // Number of dot products, x-axis
  ndp_y  = dl->out_yn;  // Number of dot products, y-axis
  nx     = dl->kx;      // Size of kernel, x-axis
  ny     = dl->ky;      // Size of kernel, y-axis
  nz     = dl->kz;      // Size of kernel, z-axis
  nk     = dl->nkern;   // Number of kernels
  stride = dl->stride;  // Stride (pixels between dot products)
  pad    = dl->pad;     // Padding

  if ((strcmp(dl->bot,"data")==0) && (mod_dcn_tstim == NULL)){
    //
    //  The input is the visual stimulus, and not the test image.
    //
    if (pad != 0)
      exit_error("MOD_DCN_GETRESP_CONV","Padding not implemented on stimulus");

    s = stim[0];

    for(ki=0;ki<nk;ki++){  // For each unique kernel

      kp1 = dl->kern1[ki];  // Get a pointer to the summed RGB kernel

      x0 = 0;  // Initialize the x-offset into the input data 'd'
      y0 = 0;  // Initialize the y-offset into the input data 'd'
      for(i=0;i<ndp_x;i++){
	for(j=0;j<ndp_y;j++){  // For each location to compute dot product

	  tot = 0.0;
	  for(kx=0;kx<nx;kx++){
	    for(ky=0;ky<ny;ky++){
	      tot += kp1[kx][ky] * s[x0+kx][y0+ky];
	    }
	  }
	  dl->r[ki][i][j] = tot;  // Store response in layer struct

	  y0 += stride;  // Step the y-offset within 'd'
	}
	x0 += stride;  // Step the x-offset within 'd'
	y0 = 0;  // Reset the y-offset
      }
    }

  }else{  // The input is the output of the previous layer

    if (strcmp(dl->bot,"data")==0){  // Input is the global test image
      d = mod_dcn_tstim;
      inz = mod_dcn_zn;
    }else{

      layin = mod_dcn_lay_for_name(dl->bot);

      if (pad == 0){
	d = layin->r;
      }else{
	//
	//  *** WYETH - ultimately we could do this once at 'prep' time,
	//              and use offsets for the response computation.
	//
	d = mod_dcn_get_padded_r(layin,pad);
      }
      inz = layin->out_zn;
      inx = layin->out_xn;
    }

    //printf("NUMBER______OF_____KERN:  %d   nz= %d\n",nk,nz);
    ngr = inz/nz;  // Number of goups
    if (ngr != dl->group)
      exit_error("MOD_DCN_GETRESP_CONV","Group size error");

    grp_sz = dl->nkern/ngr;

    cnt = 0;  // Counter for determining if/when to start new group
    z0 = 0;   // Offset only changes for subsequent groups, if any

    for(ki=0;ki<nk;ki++){  // For each unique kernel

      kern = dl->kern[ki];  // Get a pointer to this kernel
      bias = dl->bias[ki];  // Bias value

      if (cnt >= grp_sz){  // If we have finished the first group
	z0 += nz;  // increment the kernel offset
	cnt = 0;
      }
      //printf("ki %d  z0 %d\n",ki,z0);

      x0 = 0;  // Initialize x- and y-offsets into the input data 'd'
      y0 = 0;

      for(i=0;i<ndp_x;i++){    // For each location to compute dot product,
	for(j=0;j<ndp_y;j++){  //   multiply kernel 'kern' by input data 'd'.

	  tot = 0.0;
	  for(kx=0;kx<nx;kx++){
	    for(ky=0;ky<ny;ky++){
	      for(kz=0;kz<nz;kz++){
		tot += kern[kz][kx][ky] * d[z0+kz][x0+kx][y0+ky];
	      }
	    }
	  }
	  dl->r[ki][i][j] = tot + bias;  //  Store response in layer struct

	  y0 += stride;  // Step the y-offset within 'd'
	}
	x0 += stride;  // Step the x-offset within 'd'
	y0 = 0;        // Reset the y-offset
      }

      cnt += 1;
    }

    if (pad > 0){
      free_3d_farray(d,inz,inx+2*pad,inx+2*pad);
      //printf("____Freeing 3d Farray %d %d %d\n",inz,inx+2*pad,inx+2*pad);
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_DCN_GETRESP_RELU                           */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_getresp_relu(dl,stim)
     struct dcn_layer_struct *dl;   // Pointer to layer
     float ***stim;                 // Stimulus data, or NULL
{
  int i,j,k;
  int xn,yn,zn;
  float ***d,***r;
  struct dcn_layer_struct *layin;

  xn = dl->out_xn;
  yn = dl->out_yn;
  zn = dl->out_zn;

  layin = mod_dcn_lay_for_name(dl->bot);

  if ((layin->out_xn != xn) || (layin->out_yn != yn) || (layin->out_zn != zn))
    exit_error("MOD_DCN_GETRESP_RELU","Input / Output size mismatch");

  d = layin->r;
  r = dl->r;

  for(i=0;i<zn;i++){
    for(j=0;j<xn;j++){
      for(k=0;k<yn;k++){
	if (d[i][j][k] < 0.0)
	  r[i][j][k] = 0.0;
	else
	  r[i][j][k] = d[i][j][k];
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_DCN_GETRESP_POOL                           */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_getresp_pool(dl,stim)
     struct dcn_layer_struct *dl;   // Pointer to layer
     float ***stim;                 // Stimulus data, or NULL
{
  int i,j,k;
  int pi,pj,xn,yn,zn,kx,ky,xi,yi,stride,mxi,myi;
  float ***d,***r,**t,max;
  struct dcn_layer_struct *layin;

  kx = dl->kx;
  ky = dl->kx;
  stride = dl->stride;

  if (strcmp(dl->pool,"MAX")!=0)
    exit_error("MOD_DCN_GETRESP_POOL","Pooling type must be 'MAX'");

  xn = dl->out_xn;
  yn = dl->out_yn;
  zn = dl->out_zn;

  layin = mod_dcn_lay_for_name(dl->bot);

  d = layin->r;
  r = dl->r;

  for(i=0;i<zn;i++){
    t = d[i];  // Pointer to 2D response plane

    xi = 0;
    yi = 0;
    for(j=0;j<xn;j++){
      for(k=0;k<yn;k++){

	max = t[xi][yi];  // Start with 1st value
	if (mod_dcn_vistype == NULL){
	  //
	  //  Simply find the max
	  //
	  for(pi=0;pi<kx;pi++){
	    for(pj=0;pj<ky;pj++){
	      if (t[xi+pi][yi+pj] > max)
		max = t[xi+pi][yi+pj];
	    }
	  }
	}else{  // Track for Visualization
	  //
	  //  Track and save (x,y) coord of max
	  //
	  mxi = xi;
	  myi = yi;
	  for(pi=0;pi<kx;pi++){
	    for(pj=0;pj<ky;pj++){
	      if (t[xi+pi][yi+pj] > max){
		max = t[xi+pi][yi+pj];
		mxi = xi+pi;
		myi = yi+pj;
	      }
	    }
	  }
	  dl->vis_pool_x[i][j][k] = mxi;
	  dl->vis_pool_y[i][j][k] = myi;
	}

	r[i][j][k] = max;

	yi += stride;
      }
      xi += stride;
      yi = 0;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_DCN_GETRESP_LRN                            */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_getresp_lrn(dl,stim)
     struct dcn_layer_struct *dl;   // Pointer to layer
     float ***stim;                 // Stimulus data, or NULL
{
  int i,j,k;
  int zi,z0,z1,xn,yn,zn,lsize,zoff;
  float ***d,***r,alpha,beta,ac,tot;
  struct dcn_layer_struct *layin;

  xn = dl->out_xn;
  yn = dl->out_yn;
  zn = dl->out_zn;

  lsize = dl->loc_size;
  alpha = dl->alpha;
  beta  = dl->beta;
  //printf("  loc size %d   alpha %f   beta %f\n",lsize,alpha,beta);

  ac = alpha / (float)lsize;  // Constant for normalization equation

  zoff = lsize/2;

  layin = mod_dcn_lay_for_name(dl->bot);

  d = layin->r;
  r = dl->r;

  for(i=0;i<zn;i++){

    z0 = i-zoff;
    if (z0 < 0)
      z0 = 0;

    z1 = i+zoff;
    if (z1 >= zn)
      z1 = (zn-1);

    for(j=0;j<xn;j++){
      for(k=0;k<yn;k++){

	tot = 0.0;
	for(zi=z0;zi<=z1;zi++){
	  tot += d[zi][j][k] * d[zi][j][k];
	}

	r[i][j][k] = d[i][j][k] / (float)pow(1.0 + ac*tot,beta);

	if (0){
	  // Note, values in the cat image are being remapped from
	  //   things like -137 to 281 down to 1.1 t 3.0 ish.
	  if ((i==0) && (j == 0) && (k == yn-1)){
	    printf("   Tot:  %f   beta,alpha  %f %f\n",tot,beta,alpha);
	    printf("   Norm:  %f --> %f\n",d[i][j][k],r[i][j][k]);
	  }
	}

      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_DCN_GETRESP_INPROD                          */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_getresp_inprod(dl,stim)
     struct dcn_layer_struct *dl;   // Pointer to layer
     float ***stim;                 // Stimulus data, or NULL
{
  int i,j,k;
  int ki,inx,iny,inz,kk,nkern,tn;
  float tot,*kdata,*t;
  struct dcn_layer_struct *layin;

  nkern = dl->kx;

  layin = mod_dcn_lay_for_name(dl->bot);
  inx = layin->out_xn;
  iny = layin->out_yn;
  inz = layin->out_zn;

  tn = inx * iny * inz;

  if (dl->ky != tn)
    exit_error("MOD_DCN_GETRESP_INPROD","Kernel size mismatch");

  //printf("    Input layer is:  %d x %d x %d\n",inx,iny,inz);
  //printf("    Number of kernels:  %d\n",nkern);
  //printf("      Length of kernels:  %d\n",dl->ky);

  if ((inz == 1) && (inx == 1) && (iny == dl->ky)){
    t = layin->r[0][0];
  }else{
    //
    //  Arrange responses from previous layer in a single array
    //

    t = get_farray(tn);
    kk = 0;
    for(k=0;k<inz;k++){
      for(j=(iny-1);j>=0;j--){
	for(i=0;i<inx;i++){
	  t[kk] = layin->r[k][i][j];
	  kk += 1;
	}
      }
    }
  }

  for(ki=0;ki<nkern;ki++){  // For each kernel

    kdata = dl->kern[0][0][ki];  // [1][1][kx][ky]

    tot = 0.0;
    for(i=0;i<tn;i++)
      tot += t[i] * kdata[i];

    dl->r[0][0][ki] = tot + dl->bias[ki];
  }

  if ((inz == 1) && (inx == 1) && (iny == dl->ky))
    ;
  else
    myfree(t);
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_DCN_GETRESP_SFTMAX                          */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_getresp_sftmax(dl,stim)
     struct dcn_layer_struct *dl;   // Pointer to layer
     float ***stim;                 // Stimulus data, or NULL
{
  int i;
  int n,maxi;
  float *d,*tr,tot,vmax;
  struct dcn_layer_struct *layin;

  layin = mod_dcn_lay_for_name(dl->bot);

  if ((layin->out_zn != 1) || (layin->out_xn != 1))
    exit_error("MOD_DCN_GETRESP_SFTMAX","Expecting 1D input");

  n = layin->out_yn;
  d = layin->r[0][0];
  tr = dl->r[0][0];

  vmax = max_of_farray(d,n);

  tot = 0.0;
  for(i=0;i<n;i++)
    tot += exp(d[i] - vmax);

  for(i=0;i<n;i++)
    tr[i] = exp(d[i] - vmax) / tot;

  maxi = max_index_farray(tr,n);

  //printf("    max_index = %d\n",maxi);
  //printf("    max_value = %f\n",tr[maxi]);
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_DCN_01_GET_RESPONSE                          */
/*                                                                           */
/*  This routine computes the response of the DCN model to the stimlus       */
/*  stored in 's'.  The response is stored in 'r'.  The model has already    */
/*  been configured by a call to 'mod_dcn_01_prep()'.                        */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_01_get_response(m,s,r)
     struct model_struct *m;        // Model parameters
     struct stim_struct *s;         // Stimulus parameters
     struct response_struct *r;     // Response parameters
{
  int i,j;
  int xn,yn,tn,ti;
  float ***stim;
  char outfile[SLEN];
  struct dcn_layer_struct *dl;

  xn = mod_dcn_xn;
  yn = mod_dcn_yn;
  tn = mod_dcn_tn;

  stim = get_3d_farray(1,xn,yn);

  for(ti=0;ti<tn;ti++){  // For each frame in the stimulus
    //
    //  Extract frame 'ti' from the visual stimulus movie.
    //
    if (strcmp(mod_dcn_stimform,"3rgb")==0){
      for(i=0;i<xn;i++){
	for(j=0;j<yn;j++){
	  mod_dcn_tstim[0][i][j] = s->drgb_l[0][i][j][ti];
	  mod_dcn_tstim[1][i][j] = s->drgb_l[1][i][j][ti];
	  mod_dcn_tstim[2][i][j] = s->drgb_l[2][i][j][ti];
	}
      }
    }else{
      for(i=0;i<xn;i++){
	for(j=0;j<yn;j++){
	  stim[0][i][j] = s->d[i+1][j+1][ti+1];  // Use ti^th stimulus frame
	}
      }
    }

    //printf("   "); fflush(stdout);
    for(i=0;i<=mod_dcn_maxli;i++){  // Only compute layers needed for response
      dl = mod_dcn_lay[i];
      if (strcmp(dl->type,"Convolution")==0){
	mod_dcn_getresp_conv(dl,stim);
      }else if (strcmp(dl->type,"ReLU")==0){
	mod_dcn_getresp_relu(dl,stim);
      }else if (strcmp(dl->type,"Pooling")==0){
	mod_dcn_getresp_pool(dl,stim);
      }else if (strcmp(dl->type,"LRN")==0){
	mod_dcn_getresp_lrn(dl,stim);
      }else if (strcmp(dl->type,"InnerProduct")==0){
	mod_dcn_getresp_inprod(dl,stim);
      }else if (strcmp(dl->type,"Softmax")==0){
	mod_dcn_getresp_sftmax(dl,stim);
      }else{
	printf("  Layer type:  %s\n",dl->type);
	exit_error("MOD_DCN_01_GET_RESPONSE","Unknown layer type.");
      }
      //printf(" %s",dl->name); fflush(stdout);


      //if (noise_flag == 1){
      //; // Call routine to add Gaussian noise here
      //}

      if (mod_dcn_resp_3d_fname != NULL){
	if (strcmp(mod_dcn_resp_3d_fname,"NULL")!=0){
	  sprintf(outfile,"%s.%s.3d",mod_dcn_resp_3d_fname,dl->name);
	  mod_dcn_write_3d_resp(outfile,dl);
	}
      }

      if (mod_dcn_vistype != NULL){
	//
	//  VISUALIZATION
	//
	if (strcmp(mod_dcn_vislay,dl->name)==0){
	  mod_dcn_visualize_01(dl);
	  mylog(mylogf,"  Exiting after writing visualization output.\n\n");
	  exit(0);
	}
      }
    }

    if (tn > 1){
      //  If there are multiple stimulus frames, save the results in 'rt'
      mod_dcn_resp_tsave(r,ti);
    }

    //printf("\n");
  }

  free_3d_farray(stim,1,xn,yn);

  mod_dcn_handover(r);
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_DCN_STIM_RAW_TEXT_READ                       */
/*                                                                           */
/*  Read a stimulus from a raw text file.                                    */
/*                                                                           */
/*  *** WYETH - this is a hack because Dean's text images do not match       */
/*      exactly the frameset text image, which is read in                    */
/*      'stim_util.c' by 'get_stim_3rgb_ppl_frameset()'                      */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_stim_raw_text_read(infile)
     char *infile;
{
  FILE *fopen(),*fin;
  int i,j,k;
  int ns,t,xn,yn,zn,dxn,dyn,dzn;

  xn = mod_dcn_xn;
  yn = mod_dcn_yn;
  zn = mod_dcn_zn;

  // ******
  // WYETH 2019 Sep: NOTE, this code could instead call the following to
  //   read in the image (and corret the RGB order):
  // float ***data_util_3rgb_raw1_read(infile,rxn,ryn,rzn)
  // ******

  if ((fin = fopen(infile,"r"))==NULL){
    printf("  File:  %s\n",infile);
    exit_error("MOD_DCN_STIM_RAW_TEXT_READ","Cannot open file");
  }

  ns = fscanf(fin,"%d %d %d",&dzn,&dxn,&dyn);
  //printf("    Image size:  %d %d %d\n",dxn,dyn,dzn);

  if ((dxn != xn) || (dyn != yn)){
    printf("  Expecting %d x %d images, but found %d x %d\n",xn,yn,dxn,dyn);
    exit_error("MOD_DCN_STIM_RAW_TEXT_READ","Image size mis-match");
  }

  if (dzn != 3)
    exit_error("MOD_DCN_STIM_RAW_TEXT_READ","Image depth mis-match");

  for(i=0;i<dzn;i++){         //   For R,G,B planes
    for(k=(yn-1);k>=0;k--){   //     For each row (starting with top row),
      for(j=0;j<xn;j++){      //       read from left to right.
	ns = fscanf(fin,"%f",&(mod_dcn_tstim[i][j][k]));
      }
    }
  }

/*** TESTING AND DEBUGGING
  {
    float rmin,rmax;
    get_min_max_2d_farray(mod_dcn_tstim[0],xn,yn,&rmin,&rmax);
    printf("0 min,max = %f %f\n",rmin,rmax);
    get_min_max_2d_farray(mod_dcn_tstim[1],xn,yn,&rmin,&rmax);
    printf("1 min,max = %f %f\n",rmin,rmax);
    get_min_max_2d_farray(mod_dcn_tstim[2],xn,yn,&rmin,&rmax);
    printf("2 min,max = %f %f\n",rmin,rmax);

    // -104.00 to  150.99 - Blue
    // -116.66 to  138.33 - Green
    // -122.67 to  132.32 - Red
    // 
  }
***/


  fclose(fin);
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_DCN_OVERRIDE_RESP_CHECK                       */
/*                                                                           */
/*  Check the response format once before we start, so we do not have to do  */
/*  this repeatedly.                                                         */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_override_resp_check(r)
     struct response_struct *r; // Response params
{
  int i;
  int xi,yi,zi;
  struct dcn_layer_struct *lay;

  printf("  MOD_DCN_OVERRIDE_RESP_CHECK\n");

  if (mod_dcn_tn != 1)
    exit_error("MOD_DCN_OVERRIDE_RESP_CHECK","Value of 'tn' must be 1.");

  for(i=0;i<r->n;i++){ // For each response requested

    if ((r->rformat[i] == 1) ||  // 'save_pop_unit_as_'
	(r->rformat[i] == 2)){   // 'save_pop_grid_as_'

      // Example:                                      x y z
      // save_pop_unit_as_conv1_0_0_0  f 1000  conv1   0 0 0  raw

      lay = mod_dcn_lay_for_name(r->plname[i]);
      if (lay == NULL){
	printf("  Layer name:  %s\n",r->plname[i]);
	exit_error("MOD_DCN_OVERRIDE_RESP_CHECK","Layer name not found");
      }

      xi = r->xi[i];
      yi = r->yi[i];
      zi = r->zi[i];

      if ((xi < 0)||(xi >= lay->out_xn))
	exit_error("MOD_DCN_OVERRIDE_RESP_CHECK","'xi' is out of bounds");
      if ((yi < 0)||(yi >= lay->out_yn))
	exit_error("MOD_DCN_OVERRIDE_RESP_CHECK","'yi' is out of bounds");
      if ((zi < 0)||(zi >= lay->out_zn))
	exit_error("MOD_DCN_OVERRIDE_RESP_CHECK","'zi' is out of bounds");

      //printf("%s  xi,yi,zi = %d %d %d\n",r->plname[i],xi,yi,zi);

    }else{
      printf("  Response format:  %d\n",r->rformat[i]);
      exit_error("MOD_DCN_HANDOVER","Format not implemented");
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_DCN_OVERRIDE_RESP                           */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_override_resp(r,img_name,outfile)
     struct response_struct *r;  // Response params
     char *img_name;             // Name of image file
     char *outfile;
{
  FILE *fopen(),*fout;
  int i;
  int xi,yi,zi,xmax,ymax,xmin,ymin;
  float vmax,vmin;
  char tstr[SLEN];
  struct dcn_layer_struct *lay;

  //
  //  Open file
  //
  if ((fout = fopen(outfile,"a"))==NULL){
    //
    // WYETH - this sleep not needed any more - it was a bug with too many
    //         open files.
    sleep(5);
    if ((fout = fopen(outfile,"a"))==NULL){
      printf("  FILENAME:  %s\n",outfile);
      exit_error("MOD_DCN_OVERRIDE_RESP","Cannot open file");
    }
    printf("  Opened after sleeping for 5 s\n");
  }

  for(i=0;i<r->n;i++){ // For each response requested

    // Assume r->rformat[i] == 1 or 2, as checked by '..._check()' routine

    lay = mod_dcn_lay_for_name(r->plname[i]);
    xi = r->xi[i];
    yi = r->yi[i];
    zi = r->zi[i];

    //
    //  Report response for exact unit
    //
    //fval = lay->r[zi][xi][yi];  // There is a single value
    //sprintf(tstr,"%s %d %d %d %.4f %s\n",r->plname[i],zi,xi,yi,fval,img_name);

    //
    //  Find max response for the feature 'zi' plane
    //
    /*
    max_coord_2d_farray(lay->r[zi],lay->out_xn,lay->out_yn,&xmax,&ymax);
    vmax = lay->r[zi][xmax][ymax];
    fprintf(fout,"%s %d %d %d %.4f %s\n",r->plname[i],zi,xmax,ymax,
	    vmax,img_name);
    */

    if (strcmp(mod_dcn_wmover,"unit_fc")==0){
      vmax = lay->r[zi][xi][yi];
      fprintf(fout,"%s %d %.4f %s\n",r->plname[i],yi,vmax,img_name);
    }else if (strcmp(mod_dcn_wmover,"max")==0){
      min_max_coord_2d_farray(lay->r[zi],lay->out_xn,lay->out_yn,&xmin,&ymin,
			      &xmax,&ymax);
      vmax = lay->r[zi][xmax][ymax];
      vmin = lay->r[zi][xmin][ymin];
      fprintf(fout,"%s %d %d %d %.4f %d %d %.4f %s\n",r->plname[i],zi,xmax,ymax,
	      vmax,xmin,ymin,vmin,img_name);
    }else{
      exit_error("MOD_DCN_OVERRIDE_RESP","Unknown override 'mode'");
    }
  }

  //
  //  Close file
  //
  fclose(fout);
}
/**************************************-**************************************/
/*                                                                           */
/*                      MOD_DCN_01_GET_RESPONSE_OVERRIDE                     */
/*                                                                           */
/*  Process many images from a list in sequence, without control of 'wm'.    */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_01_get_response_override(m,s,r)
     struct model_struct *m;        // Model parameters
     struct stim_struct *s;         // Stimulus parameters
     struct response_struct *r;     // Response parameters
{
  int i,j;
  int xn,yn,tn,nf;
  char fullname[SLEN],**flist;
  struct dcn_layer_struct *dl;

  printf("  MOD_DCN_01_GET_RESPONSE_OVERRIDE\n");

  //
  //  Check validity of responses requested
  //
  mod_dcn_override_resp_check(r);

  xn = mod_dcn_xn;
  yn = mod_dcn_yn;
  tn = mod_dcn_tn;

  read_2d_carray(mod_dcn_wmover_flist,1,0,&flist,&nf);
  printf("    Found %d file names in %s\n",nf,mod_dcn_wmover_flist);


  if (mod_dcn_tstim == NULL)
    mod_dcn_tstim = get_3d_farray(3,xn,yn);

  for(i=0;i<nf;i++){  // For each file in list
    //
    //  Read stimulus from file
    //
    sprintf(fullname,"%s%s",mod_dcn_wmover_fpath,flist[i]);
    mod_dcn_stim_raw_text_read(fullname);

    for(j=0;j<=mod_dcn_maxli;j++){  // Only compute layers needed for response
      dl = mod_dcn_lay[j];
      if (strcmp(dl->type,"Convolution")==0){
	mod_dcn_getresp_conv(dl,NULL);
      }else if (strcmp(dl->type,"ReLU")==0){
	mod_dcn_getresp_relu(dl,NULL);
      }else if (strcmp(dl->type,"Pooling")==0){
	mod_dcn_getresp_pool(dl,NULL);
      }else if (strcmp(dl->type,"LRN")==0){
	mod_dcn_getresp_lrn(dl,NULL);
      }else if (strcmp(dl->type,"InnerProduct")==0){
	mod_dcn_getresp_inprod(dl,NULL);
      }else if (strcmp(dl->type,"Softmax")==0){
	mod_dcn_getresp_sftmax(dl,NULL);
      }else{
	printf("  Layer type:  %s\n",dl->type);
	exit_error("MOD_DCN_01_GET_RESPONSE_OVERRIDE","Unknown layer type.");
      }
    }

    //
    //  Find maximum response coord for each unit / layer in .rsp file
    //
    mod_dcn_override_resp(r,flist[i],mod_dcn_wmover_ofile);
    printf("  %6d  %s\n",i,flist[i]);
  }

  printf("    Exiting after processing all stimuli.\n");
  exit(0);
}
/**************************************-**************************************/
/*                                                                           */
/*                               MOD_DCN_RUN_01                              */
/*                                                                           */
/*****************************************************************************/
void mod_dcn_run_01(m,s,r,action)
     struct model_struct *m;    // Model params
     struct stim_struct *s;     // Stimulus params
     struct response_struct *r; // Response params
     int action;                // -1-cleanup, 0-prep, 1-run
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // DO THIS FIRST
  numproc = m->process_n;

  mod_dcn_action = action;

  if ((action == 0) || (action == 2) || (action == 10)){
    mod_dcn_01_prep(m,s,r);
  }else if (action == 1){
    mod_dcn_01_get_response(m,s,r);
  }else if (action == -1){
    mod_dcn_01_done(m);
  }
}
