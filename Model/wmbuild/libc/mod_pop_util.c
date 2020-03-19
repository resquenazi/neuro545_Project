/*****************************************************************************/
/*                                                                           */
/*  mod_pop_util.c                                                           */
/*  wyeth bair                                                               */
/*                                                                           */
/*  Population model.                                                        */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "my_util.h"
#include "carray_util.h"
#include "iarray_util.h"
#include "farray_util.h"
#include "func_util.h"
#include "myrand_util.h"
#include "paramfile_util.h"
#include "mod_util.h"
#include "ifc_util.h"
#include "pop_cell_util.h"
#include "pop_util.h"
#include "pop_low_util.h"
#include "mod_conn_util.h"
#include "mod_dog_util.h"
#include "mod_mesh_util.h"
//#include "mod_gui_util.h"  WYETH 2019 to remove OpenGL/X11 dependency
//#include "mod_mon_util.h"  WYETH 2019 to remove OpenGL/X11 dependency
#include "mod.h"
#include "ifc.h"
#include "paramfile.h"

// For cluster computation  (static needed w/ initialization
static int myid = -1;         // Process ID, -1-Single, 0-Master, >0-Slave
static char *mylogf = NULL;   // Logfile name
char ggstr[LONG_SLEN];        // Global log string

// Global for pop model
struct pop_top *mpt;          // Global model structure

/**************************************-**************************************/
/*                                                                           */
/*                         MOD_POP_INIT_LGN_CELLS_RESP                       */
/*                                                                           */
/*  Create 'cell' structs for LGN cells that have been requested for output  */
/*  using 'save_spikes' or 'save_data'.                                      */
/*                                                                           */
/*****************************************************************************/
void mod_pop_init_lgn_cells_resp(laylist,nlay,r)
     struct pop_layer **laylist;
     int nlay;
     struct response_struct *r;    // Response params
{
  int i,j;
  int xi,yi,zi;
  char tname[SLEN],*datid;
  struct pop_cell *cl,*tc;
  struct pop_layer *lay;

  //  WYETH - original routine (before NEWLGN) is below this one

  for(i=0;i<r->n;i++){
    if ((strcmp(r->plname[i],"lgn_off")==0)||
	(strcmp(r->plname[i],"lgn_on")==0)){
      mylog_exit(mylogf,
	"MOD_POP_INIT_LGN_CELLS_RESP  Old LGN pop names in .rsp file");
    }
    
    // Get the layer pointer for the requested 'plname'
    lay = pop_cell_get_layer_pointer(laylist,nlay,r->plname[i]);
    if (lay == NULL){
      sprintf(ggstr,"*** Layer name for requested response:  %s\n",
	      r->plname[i]);
      mylog(mylogf,ggstr);
      mylog_exit(mylogf,"MOD_POP_INIT_LGN_CELLS_RESP  Layer name not found.");
    }

    if (pop_util_lclass_is_lgn0(lay)){  // E.g., "lgn" (OLD"retina0_gc0")
      xi = r->xi[i];
      yi = r->yi[i];
      zi = r->zi[i];

      if ((xi < 0) || (xi >= lay->xn))
	exit_error("MOD_POP_INIT_LGN_CELLS_RESP","Response xi out of bounds");
      if ((yi < 0) || (yi >= lay->yn))
	exit_error("MOD_POP_INIT_LGN_CELLS_RESP","Response yi out of bounds");
      if ((zi < 0) || (zi >= lay->zn))
	exit_error("MOD_POP_INIT_LGN_CELLS_RESP","Response zi out of bounds");
      
      sprintf(tname,"%s_%d_%d_%d",r->plname[i],xi,yi,zi);

      // Create storage for the entire list of cells at [xi][yi]
      // Assuming this was called from ..._prep, it may never be freed?
      if (lay->c[xi][yi] == NULL){  // No storage has been created here yet
	cl = (struct pop_cell *)myalloc(lay->zn*sizeof(struct pop_cell));
	for(j=0;j<lay->zn;j++){
	  pop_cell_init_cell(&(cl[j])); // Initialize all cells
	}
	lay->c[xi][yi] = cl; // Attach column to x,y coord
      }else
	cl = lay->c[xi][yi];  // Storage already exists, use pointer to it

      tc = &(cl[zi]); // A pointer to the relevant cell

      if (tc->name == NULL){
	tc->layx = xi; // Store cell coordinates
	tc->layy = yi;
	tc->layz = zi;
	tc->name = strdup(tname);
	tc->subclass = NULL;
	tc->pl = lay;
      }

      datid = r->datid[i];

      if (strcmp(datid,"spikes")==0){
	tc->savs = 1;  // Will use c->s, ns
      }else if (strcmp(datid,"vm")==0){
	tc->savv = 1;  // Will use c->vn, vm, vmt
      }else if (strcmp(datid,"rate")==0){
	tc->savv = 1;  // Will use c->vn, vm, vmt
      }else if (strcmp(datid,"gx")==0){
	tc->savg = 1;  // Will use c->gtx0
      }else if (strcmp(datid,"gi")==0){
	tc->savg = 2;  // Will use c->gti0
      }else if (strcmp(datid,"raw")==0){
	tc->savg = 3;  // Will use c->gtn0
      }else{
	printf("***  %s\n",datid);
	exit_error("MOD_POP_INIT_LGN_CELLS_RESP","Unknown Data ID");
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_POP_STORE_LGN                            */
/*                                                                           */
/*****************************************************************************/
void mod_pop_store_lgn(laylist,nlay,r)
     struct pop_layer **laylist;
     int nlay;
     struct response_struct *r;    // Response params
{
  int i;
  int xi,yi,zi,tcnt,n;
  float *ts,samp;
  char *datid;
  struct pop_cell *c;
  struct pop_layer *lay;

  // For each response
  for(i=0;i<r->n;i++){

    // Find the layer from which this response will arise
    lay = pop_low_get_layer_for_name(laylist,nlay,r->plname[i]);
    if (lay != NULL){
      
      // Check if this layer is associated with LGN response storage
      if (pop_util_lclass_is_lgn0(lay) ||  // E.g., "lgn" (OLD"retina0_gc0")
	  pop_util_lclass_is_rgci(lay)){   // E.g., mesh rgc irreg
	
	// Identify the cell for this response request
	xi = r->xi[i];
	yi = r->yi[i];
	zi = r->zi[i];
	// Get cell pointer after doing bounds checking on cell indices
	c = pop_cell_laylist_get_cell_pointer(laylist,nlay,r->plname[i],
					      xi,yi,zi);

	datid = r->datid[i];
	if (strcmp(datid,"spikes")==0){

	  // 'zi' FIRST DIMENSION, decides ON/OFF, Left/Right
	  tcnt = lay->cnt[zi][xi][yi];
	  ts   = lay->s[zi][xi][yi];

	  // Store the spike data
	  c->ns = tcnt;
	  if (c->ns < 0){
	    sprintf(ggstr,"*** LGN %d %d %d  NOT CONNECTED, NO SPIKES MADE\n",
		    xi,yi,zi);
	    mylog(mylogf,ggstr);
	    c->ns = 0;
	    c->s = (float *)NULL;
	  }else{
	    c->s = copy_farray(ts,tcnt);
	  }
	}else if (strcmp(datid,"vm")==0){
	  ;
	}else if (strcmp(datid,"rate")==0){
	  ;
	}else if (strcmp(datid,"gx")==0){
	  ;
	}else if (strcmp(datid,"gi")==0){
	  ;
	}else{
	  exit_error("MOD_POP_STORE_LGN","Unknown data ID.");
	}
      }else if (strcmp(lay->laytype,"virtual")==0){

	xi = r->xi[i];  // Identify the cell for this response request
	yi = r->yi[i];
	zi = r->zi[i];
	// Get cell pointer after doing bounds checking on cell indices
	c = pop_cell_laylist_get_cell_pointer(laylist,nlay,r->plname[i],
					      xi,yi,zi);
	
	datid = r->datid[i];
	if (strcmp(datid,"photo")==0){
	  ts = mod_mesh_get_mm_resp(xi,yi,zi,"photo",&n,&samp);
	  popc_csav_set_f(c,"photo",ts,n,samp);
	}else if (strcmp(datid,"horiz")==0){
	  ts = mod_mesh_get_mm_resp(xi,yi,zi,"horiz",&n,&samp);
	  popc_csav_set_f(c,"horiz",ts,n,samp);
	}else if (strcmp(datid,"h2")==0){
	  ts = mod_mesh_get_mm_resp(xi,yi,zi,"h2",&n,&samp);
	  popc_csav_set_f(c,"h2",ts,n,samp);
	}else if (strcmp(datid,"diff")==0){
	  ts = mod_mesh_get_mm_resp(xi,yi,zi,"diff",&n,&samp);
	  popc_csav_set_f(c,"diff",ts,n,samp);
	}else{
	  exit_error("MOD_POP_STORE_LGN","Unknown data ID for virtual.");
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_POP_STORE_LGN_SAV                          */
/*                                                                           */
/*  Save gtx, gti, and vm for LGN cells.                                     */
/*  *** NOTE: Spikes are currently handled differently.                      */
/*                                                                           */
/*****************************************************************************/
void mod_pop_store_lgn_sav(mpt,r,lgn_lay_name,lzi,gx,gi,vm,vn,tn)
     struct pop_top *mpt;
     struct response_struct *r;    // Response params
     char *lgn_lay_name;           // LGN layer name
     int lzi;                      // layer z-index
     float ***gx,***gi;            // ex and in conductance [tn]
     float ***vm;                  // Vm [vn] msec
     int vn;                       // length of 'vm'
     int tn;                       // length of 'gx' and 'gi'
{
  int i;
  int xi,yi,zi;
  char *datid;
  struct pop_cell *c;

  // For each response
  for(i=0;i<r->n;i++){                         // For each response,
    if (strcmp(lgn_lay_name,r->plname[i])==0){ // from this LGN layer,
      if (strcmp(r->datid[i],"spikes")!=0){    // which is not spikes

	zi = r->zi[i];
	if (r->zi[i] == lzi){  // If the ON/OFF index matches

	  xi = r->xi[i];
	  yi = r->yi[i];

	  // Get cell pointer after doing bounds checking on cell indices
	  c = pop_cell_laylist_get_cell_pointer(mpt->lay,mpt->nlay,
						r->plname[i],xi,yi,zi);
	  datid = r->datid[i];
	  if ((strcmp(datid,"vm")==0)||(strcmp(datid,"rate")==0)){
	    c->vm = copy_farray(vm[xi][yi],vn);  // Will use c->vn, vm, vmt
	    c->vn = vn;
	    c->dxsav = 1000.0;  // Sampling stored here
	    c->vmt = NULL;      // Signal value for response handover
	  }else if (strcmp(datid,"gx")==0){
	    c->gtx0 = copy_farray(gx[xi][yi],tn);  // Will use c->gtx0
	    c->n0 = tn;
	  }else if (strcmp(datid,"gi")==0){
	    c->gti0 = copy_farray(gi[xi][yi],tn);  // Will use c->gti0
	    c->n0 = tn;
	  }/* else if (strcmp(datid,"raw")==0){
	    tc->savg = 3;  // Will use c->gtn0 */

	  // ****** WYETH ** Should deal with spikes and raw rdog right HERE
	  // ****** WYETH ** Should deal with spikes and raw rdog right HERE
	  // ****** WYETH ** Should deal with spikes and raw rdog right HERE

	  /*
	  // Store the spike data
	  c->ns = tcnt;
	  if (c->ns < 0){
	    sprintf(ggstr,"*** LGN %d %d %d  NOT CONNECTED, NO SPIKES MADE\n",
	    xi,yi,zi);
	    mylog(mylogf,ggstr);
	    c->ns = 0;
	    c->s = (float *)NULL;
	    }else{
	    c->s = copy_farray(ts,tcnt);
	    }
	    
	  */
	}
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_POP_LGN_SET_SAV_FLAG                        */
/*                                                                           */
/*****************************************************************************/
int mod_pop_lgn_set_sav_flag(laylist,nlay,r,flag,xn,yn,lgn_lay_name,lzi)
     struct pop_layer **laylist;
     int nlay;
     struct response_struct *r;    // Response params
     int **flag;      // [xn][yn] Set to '2' if response is requested
     int xn,yn;
     char *lgn_lay_name;
     int lzi;            // Z-layer index
{
  int i;
  int xi,yi,zi,sav_cnt;

  sav_cnt = 0;  // How many saves are flagged in this layer

  for(i=0;i<r->n;i++){                         // For each response,
    if (strcmp(lgn_lay_name,r->plname[i])==0){ // from this LGN layer,
      if (strcmp(r->datid[i],"spikes")!=0){    // which is not spikes

	zi = r->zi[i];      
	if (r->zi[i] == lzi){  // If the ON/OFF index matches

	  xi = r->xi[i];  // Get coords
	  yi = r->yi[i];
	  
	  flag[xi][yi] = 2;  // set flag to save traces
	  sav_cnt += 1;
	
	}
      }
    }
  }
  return sav_cnt;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_POP_GET_MECH                             */
/*                                                                           */
/*  Return a pointer to the named "mech".                                    */
/*                                                                           */
/*****************************************************************************/
struct pop_mech *mod_pop_get_mech(pl,name)
     struct pop_layer *pl;
     char *name;
{
  int k;
  int done;
  struct pop_mech *pm;

  pm = NULL;

  k = 0;
  done = 0;
  while(!done){
    if (k < pl->nmech){
      if (strcmp(pl->mech[k]->name,name)==0){
	pm = pl->mech[k];
	done = 1;
      }else
	k += 1;
    }else
      done = 1;
  }

  if (pm == NULL){
    printf("  *** mechanism name:  %s\n",name);
    exit_error("MOD_POP_GET_MECH","Mechanism not found");
  }

  return pm;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_POP_CONFIG_CELL                            */
/*                                                                           */
/*  WYETH - this routine and POP_CELL_INIT_CELL  can set initial             */
/*     properties for cells.  There should probably be only one.             */
/*  OR, 'pop_cell_init_cell', which is called right before this routine,     */
/*    should be moved as the first thing here???                             */
/*                                                                           */
/*****************************************************************************/
void mod_pop_config_cell(o,po,pl,c,r,xi,yi,zi)
     struct onode *o;              // Head node
     struct onode *po;             // Population node
     struct pop_layer *pl;         // Layer
     struct pop_cell *c;           // Cell to config
     struct response_struct *r;    // Response params
     int xi,yi,zi;
{
  int nmdaflag,ri,nbx,ai_ocdom;
  char tstr[SLEN3],xyzstr[SLEN];
  float tau_r_ad,tau_f_ad;
  struct onode *bgo;
  struct pop_mech *mech;

  // Wyeth - this is very inefficient, if it reads the params from the
  // lists for each cell.  Should read layer params once only

  nmdaflag = 0;
  if (pl->ifcp != NULL)
    if (pl->ifcp->g_tran > 0.0)
      nmdaflag = 1;

  sprintf(xyzstr,"%d_%d_%d",xi,yi,zi);
  sprintf(tstr,"%s_%s",pl->name,xyzstr);
  c->name = strdup(tstr);
  c->layx = xi; // Position in layer
  c->layy = yi;
  c->layz = zi;
  c->pl   = pl; // Make cell point to its layer

  // Attributes
  if (pl->attrib_fn > 0){
    c->attrib_f = (float *)myalloc(pl->attrib_fn*sizeof(float)); // storage
    pop_cell_attrib_default_set(c); //  Use default value for 'ocdom'
  }else
    c->attrib_f = NULL;

  // Synaptic Connections
  c->nout = c->nin = 0;
  c->out = (struct pop_syn *)NULL;
  c->in = (struct pop_syn *)NULL;

  // csav
  if (pl->csav_n > 0){
    c->csav_cnt = get_zero_iarray(pl->csav_n);
    c->csav_f = get_2d_farray(pl->csav_n,0);  // 2nd dim is NULL
  }


  if ((strcmp(c->pl->laytype,"virtual")==0) ||
      (strcmp(c->pl->laytype,"retina0_gc0")==0)){
    return;
  }

  // Buffers to hold upcoming input spikes, and recent g
  //   100 x 0.5 ms = 50 ms DEFAULT
  nbx = onode_getpar_int_dflt(po,"circbuf_ex_n",100);
  c->ex1s = pop_cell_get_init_circbuf(nbx,0.0,0.0005,1); // To hold
  c->in1s = pop_cell_get_init_circbuf(nbx,0.0,0.0005,1); // Recent spikes

  // Background noise firing rates, added to upcoming spike buffers
  bgo = onode_get_node_type_item_val(po,"input","type","bg");
  if (bgo != NULL){
    //
    //  WYETH - THIS SEEMS TO BE REDUNDANT WITH  mod_pop_input_config_bg(io,pl)
    //  WYETH - THIS SEEMS TO BE REDUNDANT WITH  mod_pop_input_config_bg(io,pl)
    //  WYETH - THIS SEEMS TO BE REDUNDANT WITH  mod_pop_input_config_bg(io,pl)
    //  WYETH - THIS SEEMS TO BE REDUNDANT WITH  mod_pop_input_config_bg(io,pl)
    //
    //  The other routine gets called for inputs, which happens after this
    //  so these values get over-written?  THUS, this should be removed?
    //
    c->bg_x_rate = onode_getpar_flt_dflt(bgo,"ex_rate",0.0);
    c->bg_x_amp  = onode_getpar_flt_dflt(bgo,"ex_amp",1.0);
    c->bg_i_rate = onode_getpar_flt_dflt(bgo,"in_rate",0.0);
    c->bg_i_amp  = onode_getpar_flt_dflt(bgo,"in_amp",1.0);

    // WYETH - added Dec 8, 2008  WOULD WE EVER WANT TO USE bg_x_corr_f ????
    // WYETH - added Dec 8, 2008  WOULD WE EVER WANT TO USE bg_x_corr_f ????

    c->bg_x_corr_f = onode_getpar_flt_dflt(bgo,"ex_corr_frac",0.0);
    c->bg_i_corr_f = onode_getpar_flt_dflt(bgo,"in_corr_frac",0.0);
    //c->bg_x_corr_f = 0.0;
    //c->bg_i_corr_f = 0.0;
  }else{
    c->bg_x_rate = 0.0;
    c->bg_x_amp  = 0.0;
    c->bg_i_rate = 0.0;
    c->bg_i_amp  = 0.0;
    c->bg_x_corr_f = 0.0; // WYETH - added Dec 8, 2008
    c->bg_i_corr_f = 0.0;
  }

  c->bg_x_s = c->bg_i_s = NULL;
  c->bg_x_corr_c = c->bg_i_corr_c = NULL;
  c->bg_x_n = c->bg_i_n = -1; // No spikes generated yet

  // Hold recent dynamic conductance input
  c->ex1g = pop_cell_get_init_circbuf(100,0.0,0.00005,1); // 5 ms tot dur
  if (nmdaflag)
    c->ex1g_nmda = pop_cell_get_init_circbuf(100,0.0,0.00005,1);
  c->in1g = pop_cell_get_init_circbuf(100,0.0,0.00005,1); // 5 ms tot dur

  c->gex1t = 0.0; // Overall time for computing g, shared for e/i

  c->gex1 = 0.0;  // Set starting conductance and deriv to zero
  c->gdex1 = 0.0;

  // WYETH - getting for each cell seems inefficient
  mech = mod_pop_get_mech(pl,"ex");
  c->ex1tau = mech->tau;
  c->ex1amp = mech->amp;

  c->input_ex_sdf   = onode_getpar_flt_dflt(po,"ctx2ctx_ex_sdf",1.0);
  // WYETH - BUG FIX, was reading param name '_sdf' for tau, 7 Oct 2008
  c->input_ex_sdtau = onode_getpar_flt_dflt(po,"ctx2ctx_ex_sdtau",0.0);
  c->input_ex_sdtau /= 1000.0;  // Convert to sec

  c->gin1 = 0.0;  // Set starting conductance and deriv to zero
  c->gdin1 = 0.0;

  // WYETH - at the moment, this is required, even if not used.
  mech = mod_pop_get_mech(pl,"in");
  c->in1tau = mech->tau;
  c->in1amp = mech->amp;

  c->gex1_nmda1 = c->gex1_nmda2 = c->gex1_nmda3 = 0.0;
  c->gdex1_nmda1 = c->gdex1_nmda2 = c->gdex1_nmda3 = 0.0;
  // Analogous NMDA params for PSC will be stored for layer, not for cell

  // All-time storage, WYETH - should make optional

  c->savg = 0;

  // FIXED (non-dynamic) INPUT CONDUCTANCES
  c->n0 = -1;
  c->samp0 = 1000.0;  // Sampling for pre-computed conductances
  c->gtx0 = c->gtn0 = c->gti0 = c->gta0 = (float *)NULL;
  c->ggain = NULL;

  // CAN USE POINTERS, instead of getting for each cell
  if (pl->ifcp != NULL){
    c->cifcp = pl->ifcp;
    c->cpoissp = pl->poissp;
  }else
    exit_error("MOD_POP_CONFIG_CELL","No ifc param pointer given");

  // Get adaptation shapes (Should this be done dynamically w/ alphafunc?)
  if (c->cifcp->gbar_ad > 0.0){
    tau_r_ad = c->cifcp->tau_r_ad * c->samp0/1000.0;
    tau_f_ad = c->cifcp->tau_f_ad * c->samp0/1000.0;
    c->nad = (int)(c->cifcp->tau_f_ad*4.0 * c->samp0/1000.0);
    // WYETH - save storage and have one per layer
    c->gta_pulse = diff_exp_farray(0.0,c->cifcp->gbar_ad,tau_f_ad,tau_r_ad,
				   c->nad);

    // WYETH NEW (used to be gta0)
    // Time resolution is c->samp0, duration is 5 longer than gta pulse
    // (5 was chosen arbitrarily)
    c->ad1g = pop_cell_get_init_circbuf(c->nad+5,0.0,1.0/c->samp0,1);
    c->ad1g->maskn = c->nad;  // Circbuf needs to know the mask length

  }else{
    c->nad = 0;
    c->gta_pulse = (float *)NULL;
    c->ad1g = NULL;
  }

  // Vm storage params
  c->vm = NULL;
  c->vmd = NULL; // 2-comp
  c->vmt = NULL;
  c->vnmax = 0;  // THIS safe, 0 init should be don in pop_cell_init_cell
  c->vn = 0;
  c->savv = 0;  // Unnecessary, handled in 'pop_cell_init_cell'

  //
  //  If 'cpoissp' != NULL, we do not need to do the below.
  //
  c->h1 = 0.0005;      // 1/2 msec
  c->eps = 0.0001;     // Accuracy
  c->reft = -1.0;      // Make sure we do not start in refractory state
  c->hmin = 1.0e-16;   // Minimal allowed step size.
  c->nok = c->nbad = 0;
  c->ystart = c->cifcp->v_leak_x; // Start at leakage reversal
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_POP_INIT_POP_NULL                          */
/*                                                                           */
/*  Init a population structure with all NULL cells.                         */
/*                                                                           */
/*****************************************************************************/
void mod_pop_init_pop_null(pl)
     struct pop_layer *pl;
{
  int i,j,k;
  struct pop_cell ***c;

  pl->runflag = 0;
  //pl->laytype = NULL;   WYETH - REMOVED FOR NEWLGN

  pl->x0 = pl->y0 = 0.0;
  pl->xf = pl->yf = 1.0;
  pl->oddxoff = 0.0;
  pl->umx = pl->umy = -1.0;  // Old way; new see 'area'

  pl->nmda_p = NULL;

  pl->ifcp = NULL;
  pl->poissp = NULL;

  c = (struct pop_cell ***)myalloc(pl->xn*sizeof(struct pop_cell **));
  for(i=0;i<pl->xn;i++){
    c[i] = (struct pop_cell **)myalloc(pl->yn*sizeof(struct pop_cell *));
    for(j=0;j<pl->yn;j++){
      c[i][j] = NULL;
    }
  }
  pl->c = c;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_POP_INIT_POP_LGN                           */
/*                                                                           */
/*  Special initialization for LGN populations.                              */
/*                                                                           */
/*****************************************************************************/
void mod_pop_init_pop_lgn(o,po,pl)
     struct onode *o;       // Head onode
     struct onode *po;      // Population onode
     struct pop_layer *pl;  // Layer
{
  int i,j,k;
  int ci,cj,xn,yn,zn;
  int mblock,mblock_n,mblock_seed;
  char *mblock_alg;
  struct pop_cell ***c,*cp;
  static char ***celltype = NULL;

  xn = pl->xn;
  yn = pl->yn;
  zn = pl->zn;

  //
  //  Make a population of all NULL cells
  //
  mod_pop_init_pop_null(pl);


  //
  //  Set up storage for spike trains
  //  NOTE:  'zn' is FIRST DIMENSION
  //  Typically, zn will be 2 for OFF/ON, or 4 for R/L x OFF/ON
  //
  pl->cnt = (int ***)myalloc(zn*sizeof(int **));
  pl->s   = (float ****)myalloc(zn*sizeof(float ***));
  for(i=0;i<zn;i++){
    pl->cnt[i] = NULL;
    pl->s[i]   = NULL;
  }
  pl->rl = NULL;
  pl->rr = NULL;


  //
  //  The following is used to assign 'subclass' or 'p' or 'm' to cells
  //  Presumably, was used for .../proj/size/  for limiting to m-cells
  //  Might be able to get rid of this, and use multiple LGN pops instead
  //
  mblock = onode_getpar_int_dflt(po,"mblock",-1);
  if (mblock > 0){

    // WYETH - NEWLGN
    exit_error("MOD_POP_INIT_POP_LGN","Is this needed w/ newlgn?");
    // WYETH - NEWLGN

    mblock_alg = onode_getpar_chr_exit(po,"mblock_alg");
    mblock_seed = onode_getpar_int_exit(po,"mblock_seed");
    if (mblock_seed > 0)
      mblock_seed = -mblock_seed;

    // Create a single, static cell type array, to use for ON and OFF
    if (celltype == NULL){
      celltype = get_2d_pointer_carray(xn,yn);
      for(i=0;i<xn;i++)
	for(j=0;j<yn;j++)
	  celltype[i][j] = strdup("p");  // Default type
      
      for(i=0;i<xn;i+=mblock)
	for(j=0;j<yn;j+=mblock){
	  if (strcmp(mblock_alg,"center")==0){
	    ci = i + mblock/2;
	    cj = j + mblock/2;
	  }else{ // Assume 'random'
	    ci = i + (int)((float)mblock * myrand_util_ran2(&mblock_seed));
	    cj = j + (int)((float)mblock * myrand_util_ran2(&mblock_seed));
	    //printf("  %d %d\n",ci,cj);
	  }
	  if ((ci < xn)&&(cj < yn)){
	    myfree(celltype[ci][cj]);
	    celltype[ci][cj] = strdup("m");
	  }
	}
    }
    
    // Create storage for all cells
    c = pl->c;
    for(i=0;i<xn;i++){
      for(j=0;j<yn;j++){
	c[i][j] = (struct pop_cell *)myalloc(sizeof(struct pop_cell));
	cp = &(c[i][j][0]);      // Pointer to the cell
	pop_cell_init_cell(cp);  // Set safe values
	//cp->subclass = strdup("p");
	cp->subclass = strdup(celltype[i][j]);
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_POP_MAKE_MECH_PSG_ALPHA                      */
/*                                                                           */
/*****************************************************************************/
struct pop_mech *mod_pop_make_mech_psg_alpha(o)
     struct onode *o;     // Population onode
{
  struct pop_mech *pm;

  pm = (struct pop_mech *)myalloc(sizeof(struct pop_mech));

  pm->o  = o;

  // WYETH these may not be needed, since we point to the onode
  pm->name = onode_getpar_chr_exit(o,"name");
  pm->type = onode_getpar_chr_exit(o,"type");
  pm->tau  = onode_getpar_flt_exit(o,"tau");
  pm->amp  = onode_getpar_flt_exit(o,"amp");

  return pm;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_POP_MAKE_MECH_PSG_DOE                       */
/*                                                                           */
/*****************************************************************************/
struct pop_mech *mod_pop_make_mech_psg_doe(o)
     struct onode *o;     // Mechanism onode
{
  struct pop_mech *pm;

  pm = (struct pop_mech *)myalloc(sizeof(struct pop_mech));

  pm->o  = o;

  // WYETH these may not be needed, since we point to the onode
  pm->name   = onode_getpar_chr_exit(o,"name");
  pm->type   = onode_getpar_chr_exit(o,"type");
  pm->tau_r  = onode_getpar_flt_exit(o,"tau_r");
  pm->tau_f  = onode_getpar_flt_exit(o,"tau_f");
  pm->amp    = onode_getpar_flt_exit(o,"amp");
  pm->sdf    = onode_getpar_flt_exit(o,"sdf");
  pm->sdtau  = onode_getpar_flt_exit(o,"sdtau");

  return pm;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_POP_MAKE_MECH_SI_DS02                        */
/*                                                                           */
/*****************************************************************************/
struct pop_mech *mod_pop_make_mech_si_ds02(o)
     struct onode *o;     // Mechanism onode
{
  struct pop_mech *pm;

  pm = (struct pop_mech *)myalloc(sizeof(struct pop_mech));
  pm->name   = onode_getpar_chr_exit(o,"name");
  pm->type   = onode_getpar_chr_exit(o,"type");
  pm->o = o;

  pop_util_mech_config_mask(mylogf,pm);

  return pm;
}
/**************************************-**************************************/
/*                                                                           */
/*                       MOD_POP_MAKE_MECH_SI_PAIR_RFDIST                    */
/*                                                                           */
/*****************************************************************************/
struct pop_mech *mod_pop_make_mech_si_pair_rfdist(o)
     struct onode *o;     // Mechanism onode
{
  char *nonlin;
  struct pop_mech *pm;

  pm = (struct pop_mech *)myalloc(sizeof(struct pop_mech));
  pm->name   = onode_getpar_chr_exit(o,"name");
  pm->type   = onode_getpar_chr_exit(o,"type");
  pm->o = o;


  nonlin    = onode_getpar_chr_exit(o,"nonlin");  // Type of non-linearity
  if (strcmp(nonlin,"mask")==0){
    pop_util_mech_config_mask(mylogf,pm);
  }else if (strcmp(nonlin,"symmask")==0){
    pop_util_mech_config_mask(mylogf,pm);
  }else if ((strcmp(nonlin,"mult")==0)||(strcmp(nonlin,"divide")==0)){
    pop_util_mech_config_mask(mylogf,pm);
  }else
    exit_error("MOD_POP_MAKE_MECH_SI_PAIR_RFDIST","Unknown mask type");
  myfree(nonlin);


  return pm;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_POP_MAKE_MECH_MULTIPLY                        */
/*                                                                           */
/*****************************************************************************/
struct pop_mech *mod_pop_make_mech_multiply(o)
     struct onode *o;  // Mechanism onode
{
  struct pop_mech *pm;

  pm = (struct pop_mech *)myalloc(sizeof(struct pop_mech));

  pm->o = o;

  pm->normv = 1.0;  // Default value for normalization constant; Ignored
  pm->a = 1.0;      // not used for "multiply"
  pm->b = 1.0;      // not used for "multiply"

  pm->type = strdup("multiply");

  return pm;
}
/**************************************-**************************************/
/*                                                                           */
/*                          MOD_POP_MAKE_MECH_DIVIDE                         */
/*                                                                           */
/*****************************************************************************/
struct pop_mech *mod_pop_make_mech_divide(o)
     struct onode *o;  // Mechanism onode
{
  struct pop_mech *pm;

  pm = (struct pop_mech *)myalloc(sizeof(struct pop_mech));

  pm->o = o;

  pm->normv = 1.0;  // Not used for "divide"
  pm->a = onode_getpar_flt_exit(o,"div_a",1.0);
  pm->b = onode_getpar_flt_exit(o,"div_b",1.0);

  pm->type = strdup("divide");

  return pm;
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_POP_MAKE_MECH                            */
/*                                                                           */
/*****************************************************************************/
void mod_pop_make_mech(po,pl)
     struct onode *po;     // Population onode
     struct pop_layer *pl; // layer
{
  int k;
  struct onode *t;
  char *mtype;

  pl->nmech = onode_count_otype(po,"mech");
  pl->mech = (struct pop_mech **)myalloc(pl->nmech*sizeof(struct pop_mech *));

  // "mech"
  k = 0;
  t = po->o;
  while(t != NULL){
    if (strcmp(t->otype,"mech")==0){

      mtype = onode_getpar_chr_exit(t,"type");
      //printf("    Making mechanism  %s\n",mtype);

      if (strcmp(mtype,"psg_alpha")==0)
	pl->mech[k] = mod_pop_make_mech_psg_alpha(t);
      else if (strcmp(mtype,"psg_doe")==0)
	pl->mech[k] = mod_pop_make_mech_psg_doe(t);
      else if (strcmp(mtype,"si_ds02")==0)
	pl->mech[k] = mod_pop_make_mech_si_ds02(t);
      else if (strcmp(mtype,"si_pair_rfdist")==0)
	// WYETH - GET RID OF THIS NAME???
	pl->mech[k] = mod_pop_make_mech_si_pair_rfdist(t);
      else if (strcmp(mtype,"si_nonlin")==0)
	pl->mech[k] = mod_pop_make_mech_si_pair_rfdist(t);
      else if (strcmp(mtype,"multiply")==0)
	pl->mech[k] = mod_pop_make_mech_multiply(t);
      else if (strcmp(mtype,"divide")==0)
	pl->mech[k] = mod_pop_make_mech_divide(t);
      else{
	printf("  *** mech %s\n",mtype);
	exit_error("MOD_POP_MAKE_MECH","Unknown mechanism");
      }
      myfree(mtype);

      k += 1;
    }
    t = t->next;
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_POP_INIT_GEOM                            */
/*                                                                           */
/*  Initialize the geometry informatin.                                      */
/*                                                                           */
/*****************************************************************************/
void mod_pop_init_geom(o,po,pl)
     struct onode *o;    // Head onode
     struct onode *po;   // Population onode
     struct pop_layer *pl;
{
  char *gsource;
  struct onode *geo;

  mylog(mylogf,"  MOD_POP_INIT_GEOM\n");

  popu_init_geom(mylogf,o,po,pl,mpt->retflag);  // Basic initialization

  geo = onode_child_get_unique(po,"geometry");

  if (strcmp(pl->geomt,"irregular")==0){
    gsource = onode_getpar_chr_dflt(geo,"source",NULL);
    if (strcmp(gsource,"retina0")==0){
      //
      //  Get cone mosaic information from 'mod_mesh_util'
      //
      // WHAT DO WE NEED:
      //  1.  n - number of units
      //  2.  x[],y[] - x,y coords in stim space
      //  3.  conetype - 0,1,2 = S,M,L

      pl->irr_n = mod_mesh_get_mm_mosaic_cn();
      pl->irr_x  = copy_farray(mod_mesh_ptr_mm_mosaic_cx(),pl->irr_n);
      pl->irr_y  = copy_farray(mod_mesh_ptr_mm_mosaic_cy(),pl->irr_n);
      pl->irr_id = copy_iarray(mod_mesh_ptr_mm_mosaic_cid(),pl->irr_n);
      pl->irr_dens = mod_mesh_get_mm_mosaic_dens_deg();  // cones / deg^2

      // get 'area' and set 'umx' for microns per stim grid unit on retina
      if (pl->area != NULL){
	pl->area->umx = 1000.0 / mod_mesh_get_mm_mosaic_degpmm() * mpt->sscale;
	//printf("************WYETH HERE  umx  %f  um / stim grid on retina\n",
	//pl->area->umx);
	pl->area->umy = pl->area->umx;
      }

      pl->xn = pl->irr_n;
      pl->yn = 1;

      //printf("    %s over-writing xn,yn for irreg.\n",pl->name);

      if (mpt->retflag == 100)
	mpt->retflag = 101;
    }else{
      mylogx(mylogf,"MOD_POP_INIT_GEOM","Unknown 'source' for irreg geom");
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_POP_INIT_POP                             */
/*                                                                           */
/*  Initialize the population layer, and make storage for the cells.         */
/*                                                                           */
/*****************************************************************************/
void mod_pop_init_pop(o,po,r,pl)
     struct onode *o;    // Head onode
     struct onode *po;   // Population onode
     struct response_struct *r;
     struct pop_layer *pl;
{
  int i,j,k;
  int nattrib,cidflag;
  float tx,ty;
  char tstr[SLEN],*name_ptr,**tattrib,*maptype;
  struct pop_cell ***c,*cp;
  struct onode *spko,*t;
  struct pop_map *orimap,*odmap,*rfxmap,*sfmap;

  mylog(mylogf,"  MOD_POP_INIT_POP\n");

  pl->nmda_p = NULL;

  pl->ifcp = NULL;
  pl->poissp = NULL;

  spko = onode_child_get_unique(po,"spike_gen");
  if (spko != NULL){
    if (onode_test_chr(spko,"type","ifc") == 1){

      // WYETH - this does noise params, which probably shouldn't be part
      //  of IFC struct???
      pl->ifcp = ifc_util_get_param_o(spko);  // g_tran read here for IFC
      // This is used below in 'mod_pop_config_cell'
    }else if (onode_test_chr(spko,"type","bernoulli") == 1){
      pl->ifcp = ifc_util_get_param_poiss_o(spko);  // Params shared IFC+Poiss
      pl->poissp = ifc_util_get_param_poisson(spko); // Params specific to Pois

      dll_string_insert(pl->save_name_list,"rate");  // "vm" was added before
      // Should we remove "vm"?
    }else if (onode_test_chr(spko,"type","poisson") == 1){
      //  WYETH - Bringing this back for RGC type units
      //  WYETH - Bringing this back for RGC type units  ?? TEST THIS ??
      //  WYETH - Bringing this back for RGC type units
      pl->ifcp = NULL;
      pl->poissp = NULL;
    }else{
      //  Give a meaningful error message for bad <spike_gen> type
      sprintf(ggstr,"*** <pop> = %s, <spike_gen> type = %s\n",
	      onode_getpar_chr_ptr(po,"name"),
	      onode_getpar_chr_ptr(spko,"type"));
      mylog(mylogf,ggstr);
      mylog_exit(mylogf,"MOD_POP_INIT_POP  Unknown type in <spike_gen>\n");
    }
  }

  // WYETH REMOVE below is redundant w/ previous setting to null
  //else{
  //pl->ifcp = NULL;
  //}

  //
  //  Attributes
  //
  // WYETH - could check for an 'attrib_f' object, otherwise use this as
  //         default.

  nattrib = 7;
  tattrib = (char **)myalloc(nattrib*sizeof(char *));
  tattrib[0] = strdup("ori");
  tattrib[1] = strdup("ori_cv");
  tattrib[2] = strdup("phase");
  tattrib[3] = strdup("ocdom");
  tattrib[4] = strdup("sz1");
  tattrib[5] = strdup("conetype");
  tattrib[6] = strdup("sf");
  pop_cell_attrib_init_f(pl,nattrib,tattrib,0); // 0-cells don't exist yet!
  free_2d_carray(tattrib,nattrib);

  if (onode_test_chr(po,"attrib_conetype","012SML") == 1)
    cidflag = 1;
  else
    cidflag = 0;


  // Process all "mech" objects, e.g., synpatic PSG's, SIs
  mod_pop_make_mech(po,pl);

  c = (struct pop_cell ***)myalloc(pl->xn*sizeof(struct pop_cell **));
  for(i=0;i<pl->xn;i++){
    c[i] = (struct pop_cell **)myalloc(pl->yn*sizeof(struct pop_cell *));
    for(j=0;j<pl->yn;j++){

      if (j%2 == 1)
	tx = (pl->x0 + (float)i * pl->xf + pl->oddxoff) * pl->area->umx;
      else
	tx = (pl->x0 + (float)i * pl->xf) * pl->area->umx;
      ty = (pl->y0 + (float)j * pl->yf) * pl->area->umy;

      c[i][j] = (struct pop_cell *)myalloc(pl->zn*sizeof(struct pop_cell));
      for(k=0;k<pl->zn;k++){    // For each cell
	cp = &(c[i][j][k]);     // Pointer to the cell

	pop_cell_init_cell(cp); // Set safe blank values

	mod_pop_config_cell(o,po,pl,cp,r,i,j,k);  // Makes attrib storage

	// Assign cortical locations, in um; (0,0) is Area lower left

	if (strcmp(pl->geomt,"irregular")==0){
	  cp->cx = pl->irr_x[i];
	  cp->cy = pl->irr_y[i];
	  if (cidflag == 1){
	    pop_cell_attrib_set(cp,"conetype",(float)pl->irr_id[i],
				"MOD_POP_INIT_POP");
	  }
	}else{
	  cp->cx = tx;
	  cp->cy = ty;
	}

	// NEW, w.r.t. stim grid
	cp->rfx = pl->area->x0 + (pl->x0 + (float)i*pl->xf) * pl->area->xf;
	cp->rfy = pl->area->y0 + (pl->y0 + (float)j*pl->yf) * pl->area->yf;
      }
    }
  }
  pl->c = c;


  //
  //  Apply maps
  //
  //  *** WYETH HERE:  Maps were originally built here, but should be built
  //  *** WYETH HERE:    earlier.  Thus, eventually, we should clean this up
  //  *** WYETH HERE:    and exit here if the map does not already exist.
  //  *** WYETH HERE:  AND, we should check that 'area' is consitent in the
  //                      pop and in the map
  //
  name_ptr = onode_getpar_chr_ptr(po,"map_ori");
  if (name_ptr != NULL){

    // Get the map struct
    orimap = pop_util_get_map_by_name(mpt,name_ptr);
    if (orimap == NULL)
      mylog_exit(mylogf,"MOD_POP_INIT_POP  orimap not found\n");

    if (orimap->map2 == NULL){ // Make map if not already made
      mod_conn_onode_make_map_ori(mylogf,orimap,pl->area->xn,pl->area->yn);
    }

    maptype = onode_getpar_chr_ptr(orimap->o,"type");
    //printf("maptype = %s\n",maptype);

    //****** WYETH - 1 or 0. .....
    if (strcmp(maptype,"dir_1")==0)
      mod_conn_set_layer_ori_map(mylogf,pl,orimap,0);  // No div.by 2
    else
      mod_conn_set_layer_ori_map(mylogf,pl,orimap,1);
  }


  // OD map
  name_ptr = onode_getpar_chr_ptr(po,"map_od");
  if (name_ptr != NULL){

    // Get the map struct
    odmap = pop_util_get_map_by_name(mpt,name_ptr);
    if (odmap == NULL)
      mylog_exit(mylogf,"MOD_POP_INIT_POP  odmap not found\n");

    if (odmap->map2 == NULL){ // Make map if not already made
      mod_conn_onode_make_map_od(mylogf,odmap,pl->area->xn,pl->area->yn);
    }

    mod_conn_set_layer_od_map(mylogf,pl,odmap);
  }


  // rfx map
  name_ptr = onode_getpar_chr_ptr(po,"map_rfx");
  if (name_ptr != NULL){

    // Get the map struct
    rfxmap = pop_util_get_map_by_name(mpt,name_ptr);
    if (rfxmap == NULL)
      mylog_exit(mylogf,"MOD_POP_INIT_POP  rfxmap not found\n");

    if (rfxmap->map2 == NULL){ // Make map if not already made
      mod_conn_onode_make_map_rfx(mylogf,rfxmap,pl->area->xn,pl->area->yn,
				  pl->area->xf);
    }

    mod_conn_set_layer_rfx_map(mylogf,pl,rfxmap);
  }


  // dir map
  //
  // Check for a flag string that tells how to set direction based on ori
  // If such a param does not exist, the 'ori' value will remain as it was,
  // which is the same as if the param had the value 'ori'
  name_ptr = onode_getpar_chr_ptr(po,"map_dir_flag");
  if (name_ptr != NULL){
    mod_conn_set_layer_dir_map_flag(mylogf,pl,name_ptr);
  }


  //  SF
  name_ptr = onode_getpar_chr_ptr(po,"map_sf");
  if (name_ptr != NULL){

    // Get the map struct
    sfmap = pop_util_get_map_by_name(mpt,name_ptr);
    if (sfmap == NULL)
      mylog_exit(mylogf,"MOD_POP_INIT_POP  sfmap not found\n");

    if (sfmap->map2 == NULL){ // Make map if not already made
      mod_conn_onode_make_map_sf(mylogf,sfmap,pl->area->xn,pl->area->yn);
    }
    // WYETH - MUST SET MAP
    mod_conn_set_layer_sf_map(mylogf,pl,sfmap);  // No div.by 2
    //exit_error("HERE WYETH - mod_pop_util","SF Map under development");
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                              MOD_POP_INIT_NMDA                            */
/*                                                                           */
/*****************************************************************************/
void mod_pop_init_nmda(pl,o,po,m)
     struct pop_layer *pl;       // Layer
     struct onode *o;            // Head node
     struct onode *po;           // Population node
     struct model_struct *m;
{
  float vt;
  struct pop_layer *lay;
  struct nmda_param *tp;

  // Create NMDA g(V) table

  tp = (struct nmda_param *)myalloc(sizeof(struct nmda_param));
  pl->nmda_p = tp;  // A non-NULL value here acts as the NMDA flag

  tp->nmda_ppmv = 10;
  tp->nmda_v0 = -100;
  tp->nmda_v0f = (float)tp->nmda_v0;
  pop_util_compute_nmda_table(tp->nmda_ppmv,tp->nmda_v0,100,
			      &(tp->nmda_gv),&(tp->nmda_n));
  
  tp->nmda_frac      = onode_getpar_flt_exit(po,"nmda_frac");
  tp->nmda_alpha_1   = onode_getpar_flt_exit(po,"nmda_alpha_1");
  tp->nmda_amp_1     = 1.0; // Will be adjusted automatically
  tp->nmda_alpha_2   = onode_getpar_flt_exit(po,"nmda_alpha_2");
  tp->nmda_amp_2     = onode_getpar_flt_exit(po,"nmda_amp_2");
  tp->nmda_alpha_3   = onode_getpar_flt_exit(po,"nmda_alpha_3");
  tp->nmda_amp_3     = onode_getpar_flt_exit(po,"nmda_amp_3");
  tp->lgn_nmda_frac  = onode_getpar_flt_exit(po,"lgn_nmda_frac");
  tp->lgn_nmda_amp_1 = tp->nmda_amp_1;
  tp->lgn_nmda_amp_2 = tp->nmda_amp_2;
  tp->lgn_nmda_amp_3 = tp->nmda_amp_3;
  
  vt = pl->c[0][0][0].cifcp->v_th_x;
  pop_util_compute_intcurr_ampa_nmda_layer(myid,mylogf,pl,vt,mpt->out_prefix);
  pop_util_nmda_prep_constants_layer(pl);
  //exit_error("MOD_POP_INIT_NMDA","WYETH - next call sends m, but needs o");
  pop_util_compute_intcurr_ampa_nmda_lgn(myid,mylogf,pl,m,vt,mpt->out_prefix);
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_POP_MAKE_SAVE_NAME_LIST                       */
/*                                                                           */
/*****************************************************************************/
struct dll_slist *mod_pop_make_save_name_list()
{
  struct dll_slist *slist;

  slist = dll_string_get_init();

  dll_string_insert(slist,"vm");
  dll_string_insert(slist,"spikes");

  return slist;
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_POP_LAY_RESP_CSAV                          */
/*                                                                           */
/*  Find any <response> objects in the population onode, and use these to    */
/*  create 'csav' values for the layer, and to add to the 'save_name_list'.  */
/*                                                                           */
/*****************************************************************************/
void mod_pop_lay_resp_csav(pl,po)
     struct pop_layer *pl;     // Pop layer
     struct onode *po;         // Population node
{
  int k;
  int nresp;
  struct onode *t;

  mylog(mylogf,"  MOD_POP_LAY_RESP_CSAV\n");

  nresp = onode_count_otype(po,"response");        // Count number of itemsa
  if (nresp > 0){

    if (pl->csav_n != 0)
      exit_error("MOD_POP_LAY_RESP_CSAV","csav_n is not zero");

    sprintf(ggstr,"    %d response definitions\n",nresp);
    mylog(mylogf,ggstr);
    pl->csav_n = nresp;
    pl->csav_datid = (char **)myalloc(pl->csav_n * sizeof(char *));
    pl->csav_source = (char **)myalloc(pl->csav_n * sizeof(char *));
    pl->csav_samp = (float *)myalloc(pl->csav_n * sizeof(float));

    k = 0;
    t = onode_get_next_type(po->o,"response");
    while(t != NULL){

      pl->csav_datid[k]  = onode_getpar_chr_exit(t,"data_id");
      pl->csav_source[k] = onode_getpar_chr_exit(t,"source");
      pl->csav_samp[k]   = onode_getpar_flt_exit(t,"sampling");

      dll_string_insert(pl->save_name_list,pl->csav_datid);

      k += 1;
      t = onode_get_next_type(t->next,"response");
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_POP_MAKE_POP_LAYER                          */
/*                                                                           */
/*****************************************************************************/
struct pop_layer *mod_pop_make_pop_layer(o,po,m,r)
     struct onode *o;          // Head node
     struct onode *po;         // Population node
     struct model_struct *m;
     struct response_struct *r;
{
  int biflag;
  struct pop_layer *pl;
  struct onode *t,*ico;
  char *aname,tstr[SLEN];

  mylog(mylogf,"  MOD_POP_MAKE_POP_LAYER\n");

  t = onode_get_item_child(po,"name");
  sprintf(ggstr,"    Making layer %s\n",t->val);
  mylog(mylogf,ggstr);


  pl = popu_make_pop_layer(NULL);  // Constructor, set safe values

  pl->name = onode_getpar_chr_exit(po,"name");
  aname = onode_getpar_chr_ptr_exit(po,"area");
  pl->area = pop_util_get_area_by_name(mpt,aname); // NULL if not found

  pl->laytype = onode_getpar_chr_dflt(po,"type","default");

  pl->icon = (struct pop_icon *)myalloc(sizeof(struct pop_icon));
  pl->icon->r1 = pl->icon->g1 = pl->icon->b1 = -1.0;
  ico = onode_child_get_unique(po,"icon");
  if (ico == NULL){
    pl->icon->shape = strdup("circle");
    pl->icon->r = pl->icon->g = pl->icon->b = 0.5;
    pl->icon->nside = 0;
    pl->icon->o = NULL;
  }else{
    pl->icon->shape = onode_getpar_chr_dflt(ico,"shape","circle");
    pl->icon->r     = onode_getpar_flt_dflt(ico,"r",0.5);
    pl->icon->g     = onode_getpar_flt_dflt(ico,"g",0.5);
    pl->icon->b     = onode_getpar_flt_dflt(ico,"b",0.5);
    if (strcmp(pl->icon->shape,"cs")==0){
      pl->icon->r1     = onode_getpar_flt_dflt(ico,"r1",0.5);
      pl->icon->g1     = onode_getpar_flt_dflt(ico,"g1",0.5);
      pl->icon->b1     = onode_getpar_flt_dflt(ico,"b1",0.5);
    }else if (strcmp(pl->icon->shape,"polygon")==0){
      pl->icon->nside = onode_getpar_flt_exit(ico,"nside");
    }else if (strcmp(pl->icon->shape,"stellate")==0){
      pl->icon->nside = onode_getpar_flt_exit(ico,"npoint");
    }
    pl->icon->o     = ico;
  }

  mod_pop_init_geom(o,po,pl);  // Configure the layer geometry


  pl->save_name_list = mod_pop_make_save_name_list(); // 'vm' + 'spikes'
  dll_string_insert(pl->save_name_list,"gad");
  dll_string_insert(pl->save_name_list,"in_gi");
  dll_string_insert(pl->save_name_list,"ex_gx");


  //
  //  Process <response> objects for this layer
  //
  mod_pop_lay_resp_csav(pl,po);


  //
  //  layer type-dependent initialization
  //
  if (strcmp(pl->laytype,"default") == 0){
    sprintf(tstr,"    Default pop initialization for %s\n",pl->name);
    mylog(mylogf,tstr);
    
    pl->runflag = onode_getpar_int_dflt(po,"run_flag",1);
    
    mod_pop_init_pop(o,po,r,pl);  // pl->ifcp is created here

    if (pl->ifcp->g_tran > 0.0){
      mod_pop_init_nmda(pl,o,po,m);
    }

  }else if ((strcmp(pl->laytype,"virtual")==0) || pop_util_lclass_is_rgci(pl)){

    sprintf(tstr,"    Virtual pop initialization for %s\n",pl->name);
    mylog(mylogf,tstr);
    //
    //  WYETH - for irregular 'rgc' layer using mesh 
    //
    pl->runflag = 0;
    mod_pop_init_pop(o,po,r,pl);
  }else if (pop_util_lclass_is_lgn0(pl) ||   // "lgn" or square "retina0_gc0"
	    ((strcmp(pl->laytype,"retina0_gc0") == 0) && (mpt->retflag==100))){
    sprintf(tstr,"    LGN DOG pop initialization for %s\n",pl->name);
    mylog(mylogf,tstr);
    //
    //  WYETH - for DOG lgn, and for REGULAR 'rgc' layer using mesh
    //
    
    // Check for binocular 
    biflag = onode_getpar_int_exit(po,"binocular",0);
    if (biflag == 1){
      mpt->binoc = 1; // Binocular model
      if (pl->zn != 4)
	mylog_exit(mylogf,"MOD_POP_MAKE_POP_LAYER  Use zn 4 for binoc lgn");
    }
    pl->runflag = 0;  // Do not run

    mod_pop_init_pop_lgn(o,po,pl);

  }else{
    mylog_exit(mylogf,"MOD_POP_MAKE_POP_LAYER  No initialization");
  }

  return pl;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_POP_SI_CONFIG_DS02                          */
/*                                                                           */
/*  Direction selective interaction.                                         */
/*                                                                           */
/*****************************************************************************/
void mod_pop_si_config_ds02(o,po,io,msi,inindex)
     struct onode *o;       // Top node
     struct onode *po;      // Population node
     struct onode *io;      // Input node
     struct pop_mech *msi;  // Synaptic Interaction Mechanism
     short inindex;         // Index into population 'inlist'
{
  struct pop_layer *lay1,*lay2;

  mylog(mylogf,"    MOD_POP_SI_CONFIG_DS02\n");

  lay2 = popu_get_lay_named(mpt,po,"name","MOD_POP_SI_CONFIG_DS02 lay2");
  lay1 = popu_get_lay_named(mpt,io,"pop_origin","MOD_POP_SI_CONFIG_DS02 lay1");

  mod_conn_syn_convert_layer(mylogf,lay1,lay2,msi,inindex);
}
/**************************************-**************************************/
/*                                                                           */
/*                       MOD_POP_SI_CONFIG_PAIR_RFDIST                       */
/*                                                                           */
/*  Pair synapses based on RF distance.                                      */
/*                                                                           */
/*****************************************************************************/
void mod_pop_si_config_pair_rfdist(o,po,io,msi,inindex)
     struct onode *o;       // Top node
     struct onode *po;      // Population node
     struct onode *io;      // Input node
     struct pop_mech *msi;  // Synaptic Interaction Mechanism
     short inindex;         // Index into population 'inlist'
{
  struct pop_layer *lay1,*lay2;

  //
  //  WYETH - this is just like ...ds02 above
  //

  mylog(mylogf,"    MOD_POP_SI_CONFIG_PAIR_RFDIST\n");

  lay2 = popu_get_lay_named(mpt,po,"name",
			    "MOD_POP_SI_CONFIG_PAIR_RFDIST lay2");
  lay1 = popu_get_lay_named(mpt,io,"pop_origin",
			    "MOD_POP_SI_CONFIG_PAIR_RFDIST lay1");

  mod_conn_syn_convert_layer(mylogf,lay1,lay2,msi,inindex);
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_POP_SI_INPUT_MAIN                           */
/*                                                                           */
/*  Process synaptic interactions specified in input objects.                */
/*                                                                           */
/*****************************************************************************/
void mod_pop_si_input_main(o,po,io,msi,inindex)
     struct onode *o;       // Top node
     struct onode *po;      // Population node
     struct onode *io;      // Input node
     struct pop_mech *msi;  // Synaptic Interaction Mechanism
     short inindex;         // Index into population 'inlist'
{
  mylog(mylogf,"  MOD_POP_SI_INPUT_MAIN\n");

  if (strcmp(msi->type,"si_ds02")==0){
    mod_pop_si_config_ds02(o,po,io,msi,inindex);
  }else if (strcmp(msi->type,"si_pair_rfdist")==0){
    mod_pop_si_config_pair_rfdist(o,po,io,msi,inindex);
  }else{
    sprintf(ggstr,"*** %s\n",msi->type);
    mylog(mylogf,ggstr);
    mylog_exit(mylogf,"MOD_POP_SI_INPUT_MAIN  Unknown SI type\n");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_POP_INPUT_CONFIG_BG                         */
/*                                                                           */
/*****************************************************************************/
void mod_pop_input_config_bg(io,pl)
     struct onode *io;       // Input node
     struct pop_layer *pl;   // Layer
{
  int i,j,k;
  float xrate,irate,xamp,iamp,xcfr,icfr;
  struct pop_cell *c;

  xrate = onode_getpar_flt_dflt(io,"ex_rate",0.0);
  irate = onode_getpar_flt_dflt(io,"in_rate",0.0);
  xamp  = onode_getpar_flt_dflt(io,"ex_amp",1.0);
  iamp  = onode_getpar_flt_dflt(io,"in_amp",1.0);
  xcfr  = onode_getpar_flt_dflt(io,"ex_corr_frac",0.0);
  icfr  = onode_getpar_flt_dflt(io,"in_corr_frac",0.0);

  for(i=0;i<pl->xn;i++){
    for(j=0;j<pl->yn;j++){
      for(k=0;k<pl->zn;k++){    // For each cell
	c = &(pl->c[i][j][k]);  // Pointer to the cell

	c->bg_x_rate = xrate;
	c->bg_i_rate = irate;
	c->bg_x_amp  = xamp;
	c->bg_i_amp  = iamp;

	c->bg_x_s = c->bg_i_s = NULL;
	c->bg_x_corr_c = c->bg_i_corr_c = NULL;
	c->bg_x_n = c->bg_i_n = -1; // No spikes generated yet

	c->bg_x_corr_f = xcfr;
	c->bg_i_corr_f = icfr;
      }
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                        MOD_POP_GET_CONN_FNAME_READ                        */
/*                                                                           */
/*  Return the filename for reading the input's connections, or NULL.        */
/*                                                                           */
/*  *** Should free any non-null returned value.                             */
/*                                                                           */
/*****************************************************************************/
char *mod_pop_get_conn_fname_read(mpt,io)
     struct pop_top *mpt;    // Model top level
     struct onode *io;       // input onode
{
  char *fname,*rname;
  char tname[LONG_SLEN];

  rname = NULL;

  if (strcmp(mpt->conn_rw,"read")==0){  // If allowed to read from file

    fname = onode_getpar_chr_ptr(io,"file_read"); // ptr, do not free
    if (fname != NULL){ // If there is a request to read this input from file
  
      // Build the full input file name
      if (mpt->conn_read_dir != NULL){ // Use this dir to build full path
	sprintf(tname,"%s/%s",mpt->conn_read_dir,fname);
      }else{
	sprintf(tname,"%s",fname);
      }
      rname = strdup(tname);  // A returnable copy
    }
  }

  return rname;
}
/**************************************-**************************************/
/*                                                                           */
/*                        MOD_POP_GET_CONN_FNAME_WRITE                       */
/*                                                                           */
/*  Return the filename for writing the connections, or NULL.                */
/*                                                                           */
/*  *** Should free any non-null returned value.                             */
/*                                                                           */
/*****************************************************************************/
char *mod_pop_get_conn_fname_write(mpt,io)
     struct pop_top *mpt;    // Model top level
     struct onode *io;       // input onode
{
  char *fname,*rname;
  char tname[LONG_SLEN];

  rname = NULL;

  if (strcmp(mpt->conn_rw,"write")==0){  // If allowed to write to file
    fname = onode_getpar_chr_ptr(io,"file_write"); // ptr, do not free
    if (fname != NULL){ // If there is a request to read this input from file

      // Build the full input file name
      if (mpt->conn_write_dir != NULL){ // Use this dir to build full path
	sprintf(tname,"%s/%s",mpt->conn_write_dir,fname);
      }else{
	sprintf(tname,"%s",fname);
      }
      rname = strdup(tname);  // A returnable copy
    }
  }

  return rname;
}
/**************************************-**************************************/
/*                                                                           */
/*                               MOD_POP_INPUT                               */
/*                                                                           */
/*****************************************************************************/
void mod_pop_input(o,po,io,inindex)
     struct onode *o;    // Top node
     struct onode *po;   // Population node
     struct onode *io;   // Input node
     short inindex;      // Input index number
{
  int i,j,k;
  int syntype,lii;
  float sig,sscale,axonv,axondt,normw;
  char *intype,*tname,*dtype,*wfile,*si_name,*pre_pop_name,*fname;
  char *ori_pop,tstr[SLEN];
  struct pop_layer *lay1,*lay2,*lay,*layo1,*layo2;
  struct pop_mech *mechr,*msi,*mechp;
  struct pop_cell *c;
  struct onode *disto,*t;

  mylog(mylogf,"  MOD_POP_INPUT\n");

  // Get pointer to Post-synaptic layer
  lay2 = popu_get_lay_named(mpt,po,"name","MOD_POP_INPUT post-syn");

  intype = onode_getpar_chr_dflt(io,"type","regular");
  if (strcmp(intype,"regular")==0){

    // Get pointer to pre-syn layer
    lay1 = popu_get_lay_named(mpt,io,"pop_origin","MOD_POP_INPUT regular");

    mechr = popu_get_mech_named(mpt,lay2,io,"receptor","MOD_POP_INPUT");

    syntype = onode_getpar_int_exit(mechr->o,"syntype");

    // Distribution
    disto = onode_child_get_unique(io,"distrib");
    if (disto == NULL){
      mylog_exit(mylogf,"  MOD_POP_INPUT  No distrib for input.\n");
    }


    // WYETH- Should incorporate axon velocity into other connection routines!

    // Read from file, or create
    fname = mod_pop_get_conn_fname_read(mpt,io);
    if (fname != NULL){
      axondt = onode_getpar_flt_dflt(io,"axon_dt",0.0); // (s)
      axonv = onode_getpar_flt_dflt(io,"axon_vel",0.0); // (m/s)
      mod_conn_read_conn_t1(mylogf,lay1,lay2,fname,axonv,axondt,inindex);
      myfree(fname);
    }else{
      dtype = onode_getpar_chr_ptr_exit(disto,"type");
      if (strcmp(dtype,"ori_dist")==0){
	mod_conn_onode_ori_dist_01(mylogf,io,lay1,lay2,syntype,inindex);
      }else if (strcmp(dtype,"corr_lgn")==0){

	//  Note, "lgn" is default
	tname = onode_getpar_chr_dflt(disto,"lgn_pname","lgn");
	t = onode_get_node_type_item_val(o,"pop","name",tname); // Get LGN pop
	myfree(tname);
	if (t == NULL)
	  mylog_exit(mylogf,"MOD_POP_INPUT  Cannot find 'lgn_pname'.");

	sig = onode_getpar_flt_exit(t,"sig1"); // LGN center SD
	sscale = mpt->sscale;
	// lay1 is origin
	// lay2 is target
	mod_conn_onode_corr_on_off(mylogf,io,lay2,lay1,syntype,sig,sscale,
				   mpt->xn,mpt->yn,inindex); // 2nd -> 1st arg
      }else if (strcmp(dtype,"binoc_mask")==0){

	// WYETH NEWLGN - allow other names?  WILL NEED TO CHANGE FOR retina0
	mylog(mylogf,"********* WARNING - ASSUMING 'lgn' name for layer.\n");
	t = onode_get_node_type_item_val(o,"pop","name","lgn");
	if (t == NULL)
	  mylog_exit(mylogf,"MOD_POP_INPUT  Cannot find pop 'lgn'.");

	sig = onode_getpar_flt_exit(t,"sig1");
	sscale = mpt->sscale;
	mod_conn_onode_binoc_mask(mylogf,io,lay2,lay1,syntype,sig,sscale,
				  mpt->xn,mpt->yn,1,inindex); // 2nd -> 1st arg
      }else if (strcmp(dtype,"gabor_mask")==0){

	// Should be something like 'lgn' or 'rgc'
	// WYETH - should update 'onode_binoc_mask' to be like this:
	pre_pop_name = onode_getpar_chr_ptr_exit(disto,"mask_pop_pre");
	t = onode_get_node_type_item_val(o,"pop","name",pre_pop_name);
	if (t == NULL)
	  mylog_exit(mylogf,"MOD_POP_INPUT  Cannot find 'mask_pop_pre'.");

	sig = onode_getpar_flt_exit(t,"sig1");
	sscale = mpt->sscale;
	mod_conn_onode_binoc_mask(mylogf,io,lay2,lay1,syntype,sig,sscale,
				  mpt->xn,mpt->yn,0,inindex); // 2nd -> 1st arg
      }else{
	sprintf(ggstr,"*** %s\n",dtype);
	mylog(mylogf,ggstr);
	mylog_exit(mylogf,"MOD_POP_INPUT  Unknown distrib type\n");
      }
    }

    //wfile = onode_getpar_chr_ptr(io,"write_file"); // ptr, do not free
    wfile = mod_pop_get_conn_fname_write(mpt,io);
    if (wfile != NULL){
      mod_conn_write_conn_t1(mylogf,lay1,lay2,syntype,wfile);
      myfree(wfile);
    }


    //
    //  <input_pair>
    //
    if (onode_get_next_type(io->o,"input_pair") != NULL){
      struct input_pair *ip;

      ip = pop_util_input_pair(mpt,lay2,io,inindex);

      mod_conn_input_pair_main(mpt,lay2,io,ip);

      if (ip->z1 != NULL)	myfree(ip->z1);
      if (ip->distrib != NULL)	myfree(ip->distrib);
      myfree(ip);
    }


    // Synaptic Interaction
    si_name = onode_getpar_chr_ptr(io,"syn_int");
    if (si_name != NULL){
      msi = pop_cell_get_mech_by_name(lay2,si_name); // SI Mechanism

      mylog(mylogf," *** Using old method for SYN_INT PAIRING\n");

      if (msi == NULL){
	sprintf(ggstr,"*** Mechanism %s\n",si_name);
	mylog(mylogf,ggstr);
	sprintf(ggstr,"*** Layer %s\n",lay2->name);
	mylog(mylogf,ggstr);
	mylog_exit(mylogf,"MOD_POP_INPUT  Mechanism not found in layer.");
      }
      mod_pop_si_input_main(o,po,io,msi,inindex);
    }

    // Add this input as an "all_spk_" output for the gui/monitor
    // Name it after the pre-syn layer (lay1), which gets extracted by mon
    sprintf(tstr,"all_spk_%s",lay1->name);
    dll_string_insert(lay2->save_name_list,tstr);

  }else if (strcmp(intype,"lgn_on_off")==0){
    lay2->lgn_in_flag = 1; // Used for preparing sub-cortical inputs

    // Get pointer to origin layer and receptor mechanism
    lay1 = popu_get_lay_named(mpt,io,"pop_origin","MOD_POP_INPUT lgn_on_off");
    mechr = popu_get_mech_named(mpt,lay2,io,"receptor","MOD_POP_INPUT");

    lay2->lgn_gain     = onode_getpar_int_dflt(io,"gainflag",0);
    lay2->lgn_gain_tau = onode_getpar_flt_dflt(io,"gain_tau",0.010);
    lay2->lgn_gain_ca  = onode_getpar_flt_dflt(io,"gain_ca",1.0);
    lay2->lgn_gain_cb  = onode_getpar_flt_dflt(io,"gain_cb",1.0);
    if (lay2->lgn_gain == 1)
      lay1->lgn_gain = 1;  // Must mark LGN layer to compute gain signal

    // Read from file, or create
    fname = mod_pop_get_conn_fname_read(mpt,io);
    if (fname != NULL){

      if (strcmp(lay1->geomt,"irregular")==0){
	axondt = onode_getpar_flt_dflt(io,"axon_dt",0.0); // (s)
	axonv = onode_getpar_flt_dflt(io,"axon_vel",0.0); // (m/s)
	mod_conn_read_conn_t1(mylogf,lay1,lay2,fname,axonv,axondt,inindex);
      }else{
	normw = onode_getpar_flt_dflt(io,"normw",-1.0);
	if (normw > 0.0)
	  mylog(mylogf,"    NOTE: applying 'normw' to LGN inputs\n");
	mod_conn_read_conn_t0(mylogf,lay1,lay2,mechr,fname,normw);
      }
      // But now call a routine that sets up phase values.
      // ...

      mod_conn_onode_layer_gabor_map_set_phase(mylogf,io,lay2,mpt->sscale);

    }else{
      // WYETH REMOVE RECENT inindex added to call below
      //printf("(mod_pop_util) lay1Name %s   Lay2Name %s\n",lay1->name,lay2->name);
      //printf("(mod_pop_util) lay2 ninlist: %d\n",lay2->ninlist);

      // lay1 is pre, lay2 is post


      mod_conn_onode_layer_gabor_map(mylogf,io,lay1,lay2,mechr,
				     mpt->xn,mpt->yn,mpt->sscale,inindex);
    }


    // WYETH - LGN inputs get added in order, so this will be most recent
    //wfile = onode_getpar_chr_ptr(io,"write_file"); // ptr, do not free
    wfile = mod_pop_get_conn_fname_write(mpt,io);
    if (wfile != NULL){

      // WYETH 2013 Mar
      if (strcmp(lay1->geomt,"irregular")==0){
	mod_conn_write_conn_t1(mylogf,lay1,lay2,11,wfile);  // 11-syntype
      }else{
	lii = lay2->lgn_in_n;  // LGN input index
	mod_conn_write_conn_t0(mylogf,lay2,lii,1,wfile); // 1-syntype
      }
      myfree(wfile);
    }

    // WYETH - could be combined w/ next clause for any lgn_in_flag layer??
    dll_string_insert(lay2->save_name_list,"lgn_ga");
    if (lay2->nmda_p != NULL){
      dll_string_insert(lay2->save_name_list,"lgn_gn");
    }
    if (lay2->lgn_gain == 1){
      dll_string_insert(lay2->save_name_list,"lgn_gain");
    }

    //  This is an important indicator of whether the "lgn" inputs are of
    //  the original DOG type, or of the new type.
    if (strcmp(lay1->geomt,"irregular")==0){
      lay2->lgn_in_n = 0;  // Irregular inputs are stored in syn list.
    }else
      lay2->lgn_in_n += 1; // Keep track of total number of LGN inputs

  }else if (strcmp(intype,"lgn_pair")==0){

    lay2->lgn_in_flag = 1; // Used for preparing sub-cortical inputs

    lay2->lgn_gain = onode_getpar_int_dflt(io,"gainflag",0);

    // Get pointers to the two origin layers, and the receptor mechanism
    layo1 = popu_get_lay_named(mpt,io,"pop_origin_1","MOD_POP_INPUT lgn_pair");
    layo2 = popu_get_lay_named(mpt,io,"pop_origin_2","MOD_POP_INPUT lgn_pair");
    mechr = popu_get_mech_named(mpt,lay2,io,"receptor","MOD_POP_INPUT");
    mechp = popu_get_mech_named(mpt,lay2,io,"pair_mech","MOD_POP_INPUT");

    fname = mod_pop_get_conn_fname_read(mpt,io);
    if (fname != NULL){
      myfree(fname);
      exit_error("MOD_POP_INPUT","FILE not imp'd yet for 'lgn_pair'");
      // See above use of 'lii' for lgn input index
    }else
      mod_conn_lgn_pair(mpt,io,lay2,layo1,layo2,mechr,mechp);

    //wfile = onode_getpar_chr_ptr(io,"write_file"); // ptr, do not free
    wfile = mod_pop_get_conn_fname_write(mpt,io);
    if (wfile != NULL){
      exit_error("MOD_POP_INPUT","FILE_WRITE not imp'd yet for 'lgn_pair'");
      myfree(wfile);
    }

    dll_string_insert(lay2->save_name_list,"lgn_ga");
    if (lay2->nmda_p != NULL){
      dll_string_insert(lay2->save_name_list,"lgn_gn");
      // dll_string_insert(lay2->save_name_list,"lgn_gx");  NOT Imp'd
    }
    if (lay2->lgn_gain == 1)
      dll_string_insert(lay2->save_name_list,"lgn_gain");

    lay2->lgn_in_n += 1; // Keep track of total number of LGN inputs

  }else if (strcmp(intype,"mesh_bp_to_rgc")==0){

    //  Get pre-syn 'lay1  ('lay2' is post-syn, from above)
    lay1 = popu_get_lay_named(mpt,io,"pop_origin",
			      "MOD_POP_INPUT mesh_bp_to_rgc");

    {
      int i,j,k,l;
      int xn,yn,zn,xi,zi,cn,*ci,nsyn;
      float *cw;
      struct pop_cell *cpre,*cpost;

      xn = lay2->xn;
      yn = lay2->yn;
      zn = lay2->zn;
    
      nsyn = 0;
      for(i=0;i<xn;i++){  // Number of units
	for(j=0;j<yn;j++){  // yn = 1
	  for(k=0;k<zn;k++){  // zn = 2, 0-OFF, 1-ON

	    cpost = &(lay2->c[i][j][k]);  // Pointer to post-syn cell

	    if (cpost == NULL){
	      printf("******  cpost is null\n");
	      printf("******  %d %d %d\n",i,j,k);
	      exit(0);
	    }
	    
	    //  Match layers, 0-0, 1-1 for OFF and ON
	    zi = k;
	    mod_mesh_get_rgc_conn(i,zi,&cn,&ci,&cw);  // *** zi ignored here
	    for(l=0;l<cn;l++){ // For each synapse
	      
	      xi = ci[l];
	      cpre = &(lay1->c[xi][0][zi]);  // Pointer to pre-syn cell

	      //  1-Ex, 0.0-tdelay
	      (void)pop_cell_add_synapse(cpre,cpost,1,cw[l],0.0,inindex);
	      nsyn += 1;
	    }
	    if (cn > 0){
	      myfree(ci);
	      myfree(cw);
	    }

	  }
	}
      }
      sprintf(ggstr,"    %d synapses for 'mesh_bp_to_rgc'\n",nsyn);
      mylog(mylogf,ggstr);
    }

  }else if (strcmp(intype,"bg")==0){
    mod_pop_input_config_bg(io,lay2);
  }else{
    sprintf(ggstr,"  *** %s\n",intype);
    mylog(mylogf,ggstr);
    mylog_exit(mylogf,"MOD_POP_INPUT  Unknown input type\n");
  }



  //
  //  WYETH - THIS SHOULD BE A SEPARATE ROUTINE
  //

  // Derive or rederive any attributes following the change
  ori_pop = onode_getpar_chr_dflt(io,"attr_ori_derive_from_pop",NULL);
  if (ori_pop != NULL){
    lay = pop_util_get_layer_by_name(mpt,ori_pop);  // Pointer to pop layer
    for(i=0;i<lay2->xn;i++){
      for(j=0;j<lay2->yn;j++){
	for(k=0;k<lay2->zn;k++){ // For each cell
	  c = &(lay2->c[i][j][k]);  // Pointer to the cell
	  pop_cell_derive_ori_from_layer(mylogf,lay,c,0);
	}
      }
    }
    myfree(ori_pop);
  }
  ori_pop = onode_getpar_chr_dflt(io,"attr_dir_derive_from_pop",NULL);
  if (ori_pop != NULL){
    lay = pop_util_get_layer_by_name(mpt,ori_pop);  // Pointer to pop layer
    for(i=0;i<lay2->xn;i++){
      for(j=0;j<lay2->yn;j++){
	for(k=0;k<lay2->zn;k++){ // For each cell
	  c = &(lay2->c[i][j][k]);  // Pointer to the cell
	  pop_cell_derive_ori_from_layer(mylogf,lay,c,1);
	}
      }
    }
    myfree(ori_pop);
  }

  myfree(intype);
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_POP_SET_CONN_RW                           */
/*                                                                           */
/*  Set values for:                                                          */
/*    mpt->conn_rw                                                           */
/*    mpt->conn_read_dir                                                     */
/*                                                                           */
/*****************************************************************************/
void mod_pop_set_conn_rw(mpt,o)
     struct pop_top *mpt;    // Model top level
     struct onode *o;        // Top level onode
{
  char *tdir;
  char ggstr[LONG_SLEN];

  mylog(mylogf,"    MOD_POP_SET_CONN_RW\n");

  // State for read/writing connections from a file: "none", "read", "write"
  mpt->conn_rw = onode_getpar_chr_dflt(o,"mod_conn_rw","none");

  if (strcmp(mpt->conn_rw,"none")==0){
    mpt->conn_read_dir = NULL;  // There will be no reading
    mpt->conn_write_dir = NULL;  // There will be no writing
    return;
  }else if (strcmp(mpt->conn_rw,"read")==0){

    tdir = onode_getpar_chr_ptr(o,"mod_conn_dir_r");
    if (tdir != NULL){
      if (is_dir(tdir) == 0){
	tdir = onode_getpar_chr_ptr(o,"mod_conn_dir_r2");
	if (tdir != NULL){
	  if (is_dir(tdir) == 0){
	    sprintf(ggstr,"  mod_conn_dir_r2 = %s\n",tdir);
	    mylog(mylogf,ggstr);
	    mylogx(mylogf,"MOD_POP_SET_CONN_RW","mod_conn_dir_r2 not found");
	  }
	}else{
	  tdir = onode_getpar_chr_ptr(o,"mod_conn_dir_r");
	  sprintf(ggstr,"  mod_conn_dir_r = %s\n",tdir);
	  mylog(mylogf,ggstr);
	  mylogx(mylogf,"MOD_POP_SET_CONN_RW","mod_conn_dir_r not found");
	}
      }
    }

    if (tdir == NULL)
      mylog(mylogf,"      No default connections directory for reading.\n");
    else{
      sprintf(ggstr,"      Connections directory for reading:  %s\n",tdir);
      mylog(mylogf,ggstr);
    }
    
    mpt->conn_read_dir = tdir;
    mpt->conn_write_dir = NULL;

  }else if (strcmp(mpt->conn_rw,"write")==0){

    tdir = onode_getpar_chr_ptr(o,"mod_conn_dir_w");
    if (tdir != NULL){
      if (is_dir(tdir) == 0){
	sprintf(ggstr,"  mod_conn_dir_w = %s\n",tdir);
	mylog(mylogf,ggstr);
	mylogx(mylogf,"MOD_POP_SET_CONN_RW","mod_conn_dir_w not found");
      }
    }

    if (tdir == NULL)
      mylog(mylogf,"      No default connections directory for writing.\n");
    else{
      sprintf(ggstr,"      Connections directory for writing:  %s\n",tdir);
      mylog(mylogf,ggstr);
    }
    
    mpt->conn_read_dir = NULL;
    mpt->conn_write_dir = tdir;

  }else{
    mylogx(mylogf,"MOD_POP_SET_CONN_RW","Unknown 'mod_conn_rw' value");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                                MOD_POP_PREP                               */
/*                                                                           */
/*****************************************************************************/
void mod_pop_prep(m,r,s)
     struct model_struct *m;     // Model params
     struct response_struct *r;  // Response params
     struct stim_struct *s;      // Stimulus params
{
  int i,k;
  int gui_flag,n,xflag,lgn_flag;
  char *tstr;
  struct onode *t,*tin,*wo;
  struct pop_layer *lay;

  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first
  mylog(mylogf,"  MOD_POP_PREP\n");


  //
  //  WYETH - SHOULD USE THIS LINE, BUT PLEASE TEST FIRST
  //  WYETH - SHOULD USE THIS LINE, BUT PLEASE TEST FIRST
  //  WYETH - SHOULD USE THIS LINE, BUT PLEASE TEST FIRST
  //    instead of the code below
  //
  //mpt = popu_get_init_pop_top(mylogf,sscale,tscale,xn,yn,tn);
  //


  // Basic global params, determine duration of sim in seconds
  mpt = (struct pop_top *)myalloc(sizeof(struct pop_top));

 if  (mylogf == NULL){
   mpt->logf = NULL;
 }else
   mpt->logf = strdup(mylogf);

  mpt->tscale = onode_getpar_flt_exit(m->o,"tscale");
  mpt->sscale = onode_getpar_flt_exit(m->o,"sscale");
  mpt->xn = onode_getpar_int_exit(m->o,"xn");
  mpt->yn = onode_getpar_int_exit(m->o,"yn");
  mpt->tn = onode_getpar_int_exit(m->o,"tn");
  mpt->x1 = 0.0;
  mpt->x2 = mpt->tscale * (float)mpt->tn;
  mpt->out_prefix = NULL;
  if (r != NULL){
    if (r->ppl != NULL)
      mpt->out_prefix = param_getc_dflt(r->ppl,"outfile","zz");
  }

  mpt->binoc = 0;   // Assume monocular left eye, unless 'lgn_on_r' below
  if (onode_child_get_unique(m->o,"retina0") != NULL)
    mpt->retflag = 100; // 'retina0' indicator, might change to 101 later
  else
    mpt->retflag = 0;   // Assume DOG filter

  mpt->conn_rw = NULL;
  mpt->conn_read_dir = NULL;


  //
  //  Global settings
  //   - usage of connections files
  //
  mod_pop_set_conn_rw(mpt,m->o);


  //
  // AREAS
  //
  n = onode_count_otype(m->o,"area");
  mpt->area = (struct pop_area **)myalloc(n*sizeof(struct pop_area *));
  mpt->narea = n;
  sprintf(ggstr,"    %d areas\n",mpt->narea);
  mylog(mylogf,ggstr);

  k = 0;
  t = onode_get_next_type(m->o->o,"area");
  while(t != NULL){
    mpt->area[k] = popu_make_area(mylogf,t,mpt->sscale); // WYETH Aug 2015

    /*
    mpt->area[k] = (struct pop_area *)myalloc(sizeof(struct pop_area));
    mpt->area[k]->name = onode_getpar_chr_exit(t,"name");
    mpt->area[k]->x0   = onode_getpar_flt_exit(t,"x0");
    mpt->area[k]->y0   = onode_getpar_flt_exit(t,"y0");

    mpt->area[k]->xf   = onode_getpar_flt_dflt(t,"xf",1.0);
    mpt->area[k]->yf   = onode_getpar_flt_dflt(t,"yf",1.0);

    mpt->area[k]->xn   = onode_getpar_int_exit(t,"xn");
    mpt->area[k]->yn   = onode_getpar_int_exit(t,"yn");
    mpt->area[k]->umx  = onode_getpar_flt_dflt(t,"umx",0.0);
    mpt->area[k]->umy  = onode_getpar_flt_dflt(t,"umy",mpt->area[k]->umx);
    mpt->area[k]->sscale = mpt->sscale;
    */

    k += 1;
    t = onode_get_next_type(t->next,"area");
  }



  

  //
  // MAPS
  //
  mpt->nmap = onode_count_otype(m->o,"map");
  mpt->map = (struct pop_map **)myalloc(mpt->nmap*sizeof(struct pop_map *));
  sprintf(ggstr,"    %d maps\n",mpt->nmap);
  mylog(mylogf,ggstr);

  k = 0;
  t = onode_get_next_type(m->o->o,"map");
  while(t != NULL){
    mpt->map[k] = (struct pop_map *)myalloc(sizeof(struct pop_map));
    mpt->map[k]->o  = t; 
    mpt->map[k]->xn = 0;
    mpt->map[k]->yn = 0;
    mpt->map[k]->zn = 0;
    mpt->map[k]->map2 = NULL;

    pop_util_map_make(mpt,k);  // WYETH 2014 Apr: try to make maps here

    k += 1;
    t = onode_get_next_type(t->next,"map");
  }


  //
  //  WARP - Apply a warp to the map (2014 Apr, for Lars' project)
  //
  wo = onode_child_get_unique(m->o,"warp");
  if (wo != NULL){
    popu_warp_map(mpt,wo);
  }


  //
  // LAYERS i.e., <POP>
  //
  mpt->nlay = onode_count_otype(m->o,"pop");
  mpt->lay = (struct pop_layer **)myalloc(mpt->nlay*
					  sizeof(struct pop_layer *));
  sprintf(ggstr,"    %d layers\n",mpt->nlay);
  mylog(mylogf,ggstr);
  k = 0;
  t = onode_get_next_type(m->o->o,"pop");
  while(t != NULL){
    mpt->lay[k] = mod_pop_make_pop_layer(m->o,t,m,r); // "MECH" made in here
    k += 1;
    t = onode_get_next_type(t->next,"pop");
  }

  //
  // INPUTS
  //
  for(i=0;i<mpt->nlay;i++){
    lay = mpt->lay[i];
    t = onode_get_node_type_item_val(m->o,"pop","name",lay->name);
    if (t == NULL){
      printf("*** i=%d  %s\n",i,lay->name);
      exit_error("MOD_POP_PREP","Layer name not found");
    }

    // Create list of pointers to onodes (the input onodes)
    lay->ninlist = onode_count_otype(t,"input");
    lay->inlist = (struct onode **)myalloc(lay->ninlist*
					   sizeof(struct onode *));
    k = 0;
    tin = onode_get_next_type(t->o,"input");
    while(tin != NULL){
      mod_pop_input(m->o,t,tin,(short)k);  // Process this input
      lay->inlist[k] = tin;         // Store a pointer to the onode
      k += 1;
      tin = onode_get_next_type(tin->next,"input");
    }
  }

  /*  TO EXAMINE THE CONTENTS OF A CELL
  {
    struct pop_layer *pl;

    pl = pop_util_get_layer_by_name(mpt,"bipolar");  // Pointer to pop layer
    //pl = pop_util_get_layer_by_name(mpt,"ex");  // Pointer to pop layer
    pop_cell_print_cell(&(pl->c[0][0][0]));
    exit_error("WYETH","Devel___");
  }
  */

  //
  //  CUSTOMIZE
  //
  pop_util_customize(m,mpt);

  //
  //   WRITE .mar FILE and RETURN
  //
  if ((m->action == 10) && (m->marfile != NULL)){
    //popu_mar_write(mpt,m->marfile,2);  // Version 2
    //popu_mar_write(mpt,m->marfile,3);  // Version 3
    popu_mar_write(mpt,m->marfile,4);  // Version 4
    return;
  }else if ((m->action == 20) && (m->marfile != NULL)){
    popu_rsp_gen(mpt,m,r);
    return;
  }

  //
  //  CHECK FOR LGN LAYERS
  //
  lgn_flag = 0;
  for(i=0;i<mpt->nlay;i++){
    lay = mpt->lay[i];
    if (pop_util_lclass_is_lgn0(lay)){
      lgn_flag = 1;
    }
  }
  if ((mpt->retflag != 100) && (mpt->retflag != 101) && (lgn_flag == 0))
    mpt->retflag = -1;  // No retina

  if (mpt->retflag == -1){
    // Set bg spike units for random prob. comparison
    pop_util_init_bg_rate(mpt->lay,mpt->nlay);
  }


  //
  // ANATOMY STATS
  //
  if (myid == -1){
    if (paramfile_test_param(r->ppl,"astat_lgn_proj_layer")){
      //
      //  For Valeria's project
      //
      mod_conn_astat_lgn_postsyn(mpt,r);
      xflag = paramfile_get_int_param_default(r->ppl,"astat_exit",1);
      printf("    Exiting after computing statistics.\n\n");
      exit(0);
    }else if (paramfile_test_param(r->ppl,"astat_lms_input_post")){
      //
      //  For Ruben's project
      //
      mod_conn_astat_lms_input(mpt,r);
      xflag = paramfile_get_int_param_default(r->ppl,"astat_exit",1);
      printf("    Exiting after computing statistics.\n\n");
      exit(0);
    }else if (paramfile_test_param(r->ppl,"astat_unit_attrib_layer")){
      char *popname,*a_name,*outfile;
      struct param_pair *pp;

      // WYETH - ONLY ALLOWS ONE OF THESE, CURRENTLY
      pp = paramfile_get_first_prefix_pointer(r->ppl,"astat_unit_attrib_layer");

      popname = paramfile_get_nth_char_param_pp_or_exit(pp,0);  // Layer name
      a_name  = paramfile_get_nth_char_param_pp_or_exit(pp,1);  // Attribute
      outfile = paramfile_get_nth_char_param_pp_or_exit(pp,2);  // Outfile name

      lay = pop_cell_get_layer_pointer(mpt->lay,mpt->nlay,popname);

      popc_attrib_dump_layer(lay,outfile,a_name);

      myfree(popname);
      myfree(a_name);
      myfree(outfile);
    }
  }

  // RESPONSE LINKING
  mod_pop_init_lgn_cells_resp(mpt->lay,mpt->nlay,r);
  pop_util_response_link(r,mpt->lay,mpt->nlay,mylogf);


  // GUI
  if (myid == -1){
    gui_flag = paramfile_get_int_param_default(r->ppl,"gui_flag",0);
    if (gui_flag){
      exit_error("MOD_POP_PREP","Parameter 'gui_flag' disabled");
      // WYETH 2019 commented out to remove OpenGL, X11 dependence
      //mod_gui_pop_main(mpt->lay,mpt->nlay,mpt,m,s,r);
    }
  }

  if (m->run_mode == 0) // Normal, not interactive
    mpt->monflag = paramfile_get_int_param_default(r->ppl,"mon_flag",0);
  else
    mpt->monflag = 0;

  // MONITOR
  if ((myid == -1) && (mpt->monflag)){
    printf("*** MOD_POP_PREP Warning: 'mon_flag' disabled.\n");
    //mod_mon_init(m,r,s,mpt->lay,mpt->nlay,"monitor");
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                             MOD_POP_GET_LGN_FLAG                          */
/*                                                                           */
/*  Return 2D array of ints 'flag0' and 'flag1' that specify using 0,1       */
/*  which elements of the particular LGN population are required for         */
/*  synaptic inputs.                                                         */
/*                                                                           */
/*****************************************************************************/
int **mod_pop_get_lgn_flag(m,r,pname,onflag)
     struct model_struct *m;        // Model parameter pair list
     struct response_struct *r;     // Response data
     char *pname;                   // Name of LGN population
     int onflag;                    // 1-ON, 0-OFF
{
  int i;
  int xn,yn,nlay,**flag;
  float corr_d;
  struct pop_layer *lay;

  mylog(mylogf,"  MOD_POP_GET_LGN_FLAG\n");

  nlay = mpt->nlay;

  corr_d = onode_getpar_int_dflt(m->o,"lgn_corr_distance",1);

  xn = mpt->xn;
  yn = mpt->yn;

  // Get flag array, marking each subcortical point to contribute spikes
  if (corr_d <= 1){
    flag = get_zero_2d_iarray(xn,yn);
    for(i=0;i<nlay;i++){ // Scan all layers for LGN inputs
      lay = mpt->lay[i];
      if ((lay->runflag == 1) && (lay->lgn_in_flag == 1)){

	printf("    LayerName %s,  lgn_in_n %d \n",lay->name,lay->lgn_in_n);

	// WYETH SEND 'pname' and only count relevant LGN inputs for lgn_mflag
	// WYETH SEND 'pname' and only count relevant LGN inputs for lgn_mflag
	//   WYETH - 'pname' added, 31 March, 2009
	pop_util_set_flag_grid_input_layer(mylogf,flag,xn,yn,0,lay,pname);
      }
    }
  }else{ // Make all spike trains.  May be needed depending on corr_d
    flag = get_const_2d_iarray(xn,yn,1);
  }

  return flag;
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_POP_GET_LGN_SPIKES                          */
/*                                                                           */
/*****************************************************************************/
void mod_pop_get_lgn_spikes(m,r,rdog,lgn_name,seed0,seed1,seed2,
			    rs0,rs1,rcnt0,rcnt1)
     struct model_struct *m;        // Model parameter pair list
     struct response_struct *r;     // Response data
     float ***rdog;
     char lgn_name[];               // LGN layer name
     int seed0,seed1,seed2;
     float ****rs0,****rs1;
     int ***rcnt0,***rcnt1;
{
  int i;
  int xn,yn,tn,nlay,sv_vn,sav0,sav1;
  int **flag0,**flag1,**cnt0,**cnt1,corr_d;
  float samp,***s0,***s1,corr_p,corr_tsd,periodms;
  float ***sv0gx,***sv0gi,***sv0vm,***sv1gx,***sv1gi,***sv1vm;
  struct pop_layer *lay;

  nlay = mpt->nlay;
  
  mylog(mylogf,"  MOD_POP_GET_LGN_SPIKES\n");

  corr_d = onode_getpar_int_dflt(m->o,"lgn_corr_distance",1);
  
  xn = mpt->xn;
  yn = mpt->yn;
  tn = mpt->tn;
  samp = 1.0/mpt->tscale;

  // Get flag array, marking each subcortical point to contribute spikes
  if (corr_d <= 1){
    flag0 = get_zero_2d_iarray(xn,yn); // OFF
    flag1 = get_zero_2d_iarray(xn,yn); // ON
    for(i=0;i<nlay;i++){ // Scan all layers for LGN inputs
      lay = mpt->lay[i];
      if ((lay->runflag == 1) && (lay->lgn_in_flag == 1)){
	pop_util_set_flag_grid_input_layer(mylogf,flag0,xn,yn,0,lay,lgn_name);
	pop_util_set_flag_grid_input_layer(mylogf,flag1,xn,yn,1,lay,lgn_name);
      }
    }
  }else{ // Make all spike trains.  May be needed depending on corr_d
    flag0 = get_const_2d_iarray(xn,yn,1); // OFF
    flag1 = get_const_2d_iarray(xn,yn,1); // ON
  }

  // WYETH - the part of this for spikes is done below, in '...get_response('
  // WYETH -  by this command:  mod_pop_store_lgn(laylist,nlay,r);

  // Mark flag arrays with '2' where non-spike saves are requested
  sav0 = mod_pop_lgn_set_sav_flag(mpt->lay,mpt->nlay,r,flag0,xn,yn,lgn_name,0);
  sav1 = mod_pop_lgn_set_sav_flag(mpt->lay,mpt->nlay,r,flag1,xn,yn,lgn_name,1);

  //printf("===================== sav0,1  %d  %d\n",sav0,sav1);

  // get ON and OFF spikes, returned in msec
  mod_dog_get_spikes_2d_flag(m,r,rdog,lgn_name,1,xn,yn,tn,flag1,samp,seed2,
			     &s1,&cnt1,sav1,&sv1gx,&sv1gi,&sv1vm,&sv_vn);
  mod_dog_get_spikes_2d_flag(m,r,rdog,lgn_name,0,xn,yn,tn,flag0,samp,seed2,
			     &s0,&cnt0,sav0,&sv0gx,&sv0gi,&sv0vm,&sv_vn);
  free_2d_iarray(flag0,xn);
  free_2d_iarray(flag1,xn);

  if (corr_d > 1){
    periodms = (float)tn*mpt->tscale * 1000.0; // Assume spikes in ms
    corr_p = onode_getpar_flt_dflt(m->o,"lgn_corr_ppick",1.0);
    corr_tsd = onode_getpar_flt_dflt(m->o,"lgn_corr_tsd",0.0);
    mylog(mylogf,"    Adjusting LGN spikes for correlation.\n");
    pop_util_corr_nearby_spikes(xn,yn,s0,cnt0,periodms,corr_d,corr_p,
				corr_tsd,seed0);
    pop_util_corr_nearby_spikes(xn,yn,s1,cnt1,periodms,corr_d,corr_p,
				corr_tsd,seed1);
  }


  // Store any continuous traces into LGN cells, for response handover
  if (sav0 > 0){
    mod_pop_store_lgn_sav(mpt,r,lgn_name,0,sv0gx,sv0gi,sv0vm,sv_vn,tn);
    free_3d_farray_null(sv0gx,xn,yn,0);
    free_3d_farray_null(sv0gi,xn,yn,0);
    free_3d_farray_null(sv0vm,xn,yn,0);
  }
  if (sav1 > 0){
    mod_pop_store_lgn_sav(mpt,r,lgn_name,1,sv1gx,sv1gi,sv1vm,sv_vn,tn);
    free_3d_farray_null(sv1gx,xn,yn,0);
    free_3d_farray_null(sv1gi,xn,yn,0);
    free_3d_farray_null(sv1vm,xn,yn,0);
  }

  *rs0 = s0;  *rcnt0 = cnt0;
  *rs1 = s1;  *rcnt1 = cnt1;
}
/**************************************-**************************************/
/*                                                                           */
/*                         MOD_POP_CELL_DUMP_LGN_INPUT                       */
/*                                                                           */
/*  For debugging, dump out all LGN input coming to a particular cell.       */
/*                                                                           */
/*****************************************************************************/
void mod_pop_cell_dump_lgn_input(popname,xi,yi,zi,outfile,rdog)
     char *popname;
     int xi,yi,zi;
     char *outfile;
     float ***rdog;   // response from filters
{
  int i;
  int xn,yn,tn,n1,n0,tx,ty;
  float *tr;
  char tname[SLEN];
  struct pop_layer *lay;
  struct pop_cell *c;

  printf("  MOD_POP_CELL_DUMP_LGN_INPUT\n");

  xn = mpt->xn;
  yn = mpt->yn;
  tn = mpt->tn;

  lay = pop_cell_get_layer_pointer(mpt->lay,mpt->nlay,popname);
  
  c = &(lay->c[xi][yi][zi]);
  n1 = c->lgnin[0]->cn1;
  n0 = c->lgnin[0]->cn0;

  printf("    %d ON  %d OFF inputs\n",n1,n0);

  for(i=0;i<n1;i++){
    tx = c->lgnin[0]->cx1[i];
    ty = c->lgnin[0]->cy1[i];
    tr = &(rdog[tx+1][ty+1][1]);
    sprintf(tname,"%s_%d_%d_%d_ON_in_%d_%d",popname,xi,yi,zi,tx,ty);
    append_farray_plot(outfile,tname,tr,tn,1);
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                           MOD_POP_GET_RESP_LGN_LAY                        */
/*                                                                           */
/*  Save spike train responses into layers.                                  */
/*  Save 'rdog' responses (presumably for output?)                           */
/*  Save 'rdogg' gain signals, for later computation.                        */
/*                                                                           */
/*****************************************************************************/
void mod_pop_get_resp_lgn_lay(m,s,r,lay)
     struct model_struct *m;        // Model parameter pair list
     struct stim_struct *s;         // Stimulus data and param pair list
     struct response_struct *r;     // Response data
     struct pop_layer *lay;         // LGN layer
{
  int i;
  int xn,yn,zn,seed0,seed1,seed0r,seed1r,seed2,seed2r;
  int **cnt0,**cnt1,**cnt0r,**cnt1r,**sflag,nrsp_flagged;
  float ***s0,***s1,***s0r,***s1r;
  float ***rdog,***rdogr;
  float ***rdogg,***rdoggr;
  float ****s3;
  int ***cnt3;

  mylog(mylogf,"  MOD_POP_GET_RESP_LGN_LAY\n");

  sprintf(ggstr,"    %s\n",lay->name);
  mylog(mylogf,ggstr);

  xn = lay->xn;  // WYETH HERE changed Feb 5, 2011, from xn = mpt->xn
  yn = lay->yn;
  zn = lay->zn;

  seed0  = lay->layseed[2];  // Use layer seeds for corr spikes
  seed1  = lay->layseed[3];
  seed0r = lay->layseed[4];
  seed1r = lay->layseed[5];

  seed2  = lay->layseed[6];  // WYETH - this makes a 2nd seed for OFF, x17
                             //   but probably should be a separate seed
  seed2r = lay->layseed[7];  // WYETH - this makes a 2nd seed for OFF, x17
                             //   but probably should be a separate seed

  // Start all spike and response values at NULL
  cnt0 = cnt1 = cnt0r = cnt1r = NULL;
  s0 = s1 = s0r = s1r = NULL;
  rdog = rdogr = NULL;
  s3 = NULL;
  cnt3 = NULL;

  if (mpt->binoc == 0){
    if ((mpt->retflag >= 100) && (mpt->retflag <= 101)){   // 'mod_mesh'
      
      mod_mesh_compute_diff_all(m,s,r); // Compute PH difference everywhere

      // WYETH - HACK, hard-coded "rgc"
      // WYETH - note, seeds for correlated LGN spikes aren't used here.
      // WYETH - THESE FLAG ARRAYS COULD BE COMPUTED IN PREP???

      if (mpt->retflag == 101){
	int nflagged;
	int ***cmap;

	nflagged = popc_flag_get_lay_to_all(lay,mpt->lay,mpt->nlay,&cmap,-1.0);

        // printf("    nflagged = %d  of %d x %d x %d\n",nflagged,xn,yn,zn);

	//
	//  Set any additional flags to 1 for cells in 'lay' that are 
	//  requested by the .rsp file
	//
	nrsp_flagged = popu_flag_resp(r,lay,cmap);

	mod_mesh_get_spikes_3d_flag(m,s,r,"rgc",cmap,xn,yn,zn,&s3,&cnt3);
	free_3d_iarray(cmap,xn,yn,zn);

      }else{

	//  **** WYETH HERE - WHEN DOES THIS OPTION GET RUN ?????????
	//  **** WYETH HERE - WHEN DOES THIS OPTION GET RUN ?????????
	//  **** WYETH HERE - WHEN DOES THIS OPTION GET RUN ?????????
	//  **** WYETH HERE - WHEN DOES THIS OPTION GET RUN ?????????
	//  **** WYETH HERE - WHEN DOES THIS OPTION GET RUN ?????????
	//  **** WYETH HERE - WHEN DOES THIS OPTION GET RUN ?????????
	//  Jun 2013
	//

	// *** WYETH NOTE we are not flagging responses here as above
	// *** WYETH NOTE we are not flagging responses here as above
	// *** WYETH NOTE we are not flagging responses here as above

	sflag = mod_pop_get_lgn_flag(m,r,"rgc",0); // 0-OFF, *** "rgc" IGNORED
	mod_mesh_get_spikes_2d_flag(m,s,r,"rgc",0,xn,yn,sflag,&s0,&cnt0);
	free_2d_iarray(sflag,xn);

	sflag = mod_pop_get_lgn_flag(m,r,"rgc",1); // 0-ON, *** "rgc" IGNORED
	mod_mesh_get_spikes_2d_flag(m,s,r,"rgc",1,xn,yn,sflag,&s1,&cnt1);
	free_2d_iarray(sflag,xn);
      }

    }else{ // Default is 'mod_dog'

      // Modifies s->d, returned is [1..xn][1..yn][1..tn], DO NOT FREE
      rdog = mod_dog_get_response(s,s->repno,lay->name,lay->lgn_gain);

      if (lay->lgn_gain == 1){
	rdogg = mod_dog_gain_ptr(lay->name,0);  // 0-left eye
	rdoggr = NULL;
      }else{
	rdogg = rdoggr = NULL;
      }

      //
      //  WYETH - next call can probably be streamlined by sending the
      //  layer pointers, since it contains the layseeds and the spike storage
      //
      // LGN SPIKES get float spikes: s[xn][yn][cnt] and cnt[xn][yn]

      mod_pop_get_lgn_spikes(m,r,rdog,lay->name,seed0,seed1,seed2,&s0,&s1,
			     &cnt0,&cnt1);

      //mod_pop_get_lgn_spikes(m,r,rdog,"lgn",seed0,seed1,seed2,&s0,&s1,&cnt0,
      //&cnt1);
    }
    //s0r = s1r = NULL;     // No right eye spikes
    //cnt0r = cnt1r = NULL; // No right eye spikes
  }else{
    if (s->d_r == NULL){
      //mylog_exit(mylogf,"MOD_POP_GET_RESPONSE  Right-eye stimulus is NULL");
      mylog(mylogf,"  *** MOD_POP_GET_RESPONSE - Using 0 spikes for R-eye.\n");
      s0r = get_3d_farray(xn,yn,0);
      s1r = get_3d_farray(xn,yn,0);
      cnt0r = get_zero_2d_iarray(xn,yn);
      cnt1r = get_zero_2d_iarray(xn,yn);

      // Modifies s->d, returned is [1..xn][1..yn][1..tn], DO NOT FREE
      rdog = mod_dog_get_response(s,s->repno,lay->name,lay->lgn_gain);

      //mod_pop_get_lgn_spikes(m,r,rdog,"lgn",seed0,seed1,seed2,&s0,&s1,&cnt0,
      mod_pop_get_lgn_spikes(m,r,rdog,lay->name,seed0,seed1,seed2,&s0,&s1,
			     &cnt0,&cnt1);
    }else{
      // Modifies s->d, returned is [1..xn][1..yn][1..tn], DO NOT FREE

      rdog = mod_dog_get_response_binoc(s,s->repno,lay->name,&rdogr);
      //mod_pop_get_lgn_spikes(m,r,rdog,"lgn",seed0,seed1,seed2,&s0,&s1,&cnt0,
      mod_pop_get_lgn_spikes(m,r,rdog,lay->name,seed0,seed1,seed2,&s0,&s1,
			     &cnt0,&cnt1);
      //mod_pop_get_lgn_spikes(m,r,rdogr,"lgn",seed0r,seed1r,seed2r,
      mod_pop_get_lgn_spikes(m,r,rdogr,lay->name,seed0r,seed1r,seed2r,
			     &s0r,&s1r,&cnt0r,&cnt1r);
    }
  }


  //
  // WYETH - looking for mem leak - it comes in this 'if' but does not free
  //  anything on the first time through
  //
  if (s3 != NULL){   // Used for irreg mesh

    if (lay->s != NULL){

      //
      //  To avoid freeing non-allocated pointers.  2012 May 23
      //
      //free_4d_farray(lay->s,zn,xn,yn,0);
      {
	int i1,i2,i3;
	int ***dcnt;
	float ****dd;

	dd = lay->s;
	dcnt = lay->cnt;

	for(i1=0;i1<zn;i1++){
	  for(i2=0;i2<xn;i2++){
	    for(i3=0;i3<yn;i3++)
	      if (lay->cnt[i1][i2][i3] > 0)
		myfree(dd[i1][i2][i3]);
	    myfree(dd[i1][i2]);
	  }
	  myfree(dd[i1]);
	}
	myfree(dd);
      }

    }

    if (lay->cnt != NULL){
      free_3d_iarray(lay->cnt,zn,xn,yn);
    }

    lay->cnt = cnt3;
    lay->s = s3;
  }else{

    // Free old spike storage
    for(i=0;i<zn;i++){
      if (lay->cnt[i] != NULL)
	free_2d_farray(lay->cnt[i],xn);
      if (lay->s[i] != NULL){
	free_3d_farray(lay->s[i],xn,yn,0);  // last val ignored
      }
    }

    // Store new spikes
    lay->cnt[0] = cnt0;  lay->s[0] = s0;    //  Left OFF
    lay->cnt[1] = cnt1;  lay->s[1] = s1;    //  Left ON
    lay->rl = rdog;
    lay->rgl = rdogg;
    if (lay->zn >= 4){
      lay->cnt[2] = cnt0r;  lay->s[2] = s0r;    //  Right OFF
      lay->cnt[3] = cnt1r;  lay->s[3] = s1r;    //  Right ON
      lay->rr = rdogr;
      lay->rgr = rdoggr;
    }
  }
}
/**************************************-**************************************/
/*                                                                           */
/*                            MOD_POP_GET_RESPONSE                           */
/*                                                                           */
/*****************************************************************************/
void mod_pop_get_response(m,s,r)
     struct model_struct *m;        // Model parameter pair list
     struct stim_struct *s;         // Stimulus data and param pair list
     struct response_struct *r;     // Response data
{
  int i;
  int xn,yn,tn,nlay,**sflag;
  int seed,seed1,seed2,*seedlist;
  float ***rdog,***rdogr,***s0,***s1,tnsec;
  char outfile[SLEN],tname[SLEN],lname[SLEN];
  struct pop_layer **laylist;
  struct pop_layer *lay;

  mylog(mylogf,"  MOD_POP_GET_RESPONSE\n");


  // Set convenient pointers
  laylist = mpt->lay;
  nlay = mpt->nlay;
  xn = mpt->xn;
  yn = mpt->yn;
  tn = mpt->tn;


  // Choose overall seeds here
  seed = m->mseed[r->tsi];
  sprintf(ggstr,"    trial seed = %d\n",seed);
  mylog(mylogf,ggstr);
  if (seed > 0)
    seed = -seed;
  seed1 = (int)( 57913.0 * myrand_util_ran2(&seed)); // layseeds
  seed2 = (int)( 73119.0 * myrand_util_ran2(&seed)); // run_layer_list


  // Choose layer seeds
  seedlist = get_seeds(seed1,99999,nlay);
  for(i=0;i<nlay;i++){
    lay = mpt->lay[i];
    lay->layseed_n = 10; // WYETH arbitrary constant, should be set in prep
    if (lay->layseed != NULL)
      myfree(lay->layseed);
    lay->layseed = get_seeds(seedlist[i],77777,lay->layseed_n);
    //
    //       Regular     LGN
    //
    // [0] - gt input    unused
    // [1] - bg spks     unused
    // [2] - unused      lgn corr spikes left OFF
    // [3] - unused      lgn corr spikes left ON
    // [4] - unused      lgn corr spikes right OFF
    // [5] - unused      lgn corr spikes right ON
    // [6] - unused      lgn g noise
    // [7] - unused      lgn g noise, right
    // [8] - unused      unused
    // [9] - unused      unused
  }
  myfree(seedlist);



  //  RUN SUBCORTICAL MODULE on stimulus for each relevant layer
  if (mpt->retflag != -1){  // If there is an lgn/retina level
    mylog(mylogf,"    Running subcortical module\n");
    for(i=0;i<nlay;i++){
      lay = mpt->lay[i];
      if (pop_util_lclass_is_lgn0(lay) ||
	  (strcmp(lay->laytype,"retina0_gc0")==0)){
	mod_pop_get_resp_lgn_lay(m,s,r,lay);
      }
    }
  }

  // Preparation of LGN layers for later response handover
  //   WYETH - this seems to store spikes only, into the cells where
  //           spikes are requested.
  mod_pop_store_lgn(laylist,nlay,r);


  // INIT CORTICAL INPUT FROM SUBCORT RESPONSES
  if (mpt->retflag != -1){ // If not on-the-fly mode

    mylog(mylogf,"    Pre-computing subcortical g's and bg spikes.\n");
    tnsec = (float)tn * mpt->tscale;
    sprintf(tname,"%s.thalamocort.pl",mpt->out_prefix);
    tname[0] = '\0'; // WYETH - when changed, think about myid

    for(i=0;i<nlay;i++){  // EACH LAYER
      lay = mpt->lay[i];

      if (lay->runflag == 1){

	// Set sub-cort input and reset dynamic cortical storage
	if (lay->lgn_in_flag == 1){

	  //  Gain control, if present, is applied to LGN signals here

	  pop_util_set_gt_on_off_input_layer(mylogf,lay,m,xn,yn,tnsec,
					     lay->layseed[0],tname,mpt->binoc);
	}else{
	  pop_util_set_gad_layer(lay,xn,yn,tnsec);
	}

	// Background spikes
	pop_util_init_bg_spikes(lay,&(lay->layseed[1]),mpt->tn,mpt->tscale);
      }
    }
  }

  // RUN CORTICAL MODEL
  if ((myid == -1) && (mpt->monflag || (m->run_mode == 1))){
    printf("*** MOD_POP_GET_RESPONSE Warning: 'mon_flag' disabled.\n");
    // WYETH 2019 commented to remove OpenGL, X11 dependency
    //mod_mon_reset(m,r,s);
    // This will handle configure events, so LGN gx is not erased
    //mod_mon_monitor_check(mpt->lay,mpt->nlay);
  }

  //printf("####........  before run layer list ...### (mod_pop_util)\n");

  pop_util_run_layer_list(myid,mylogf,m,r,s,laylist,nlay,mpt->x1,mpt->x2,
			  seed2,mpt->retflag);

  // Response Dumping  - write plots for savv savg
  if (myid == -1){
    sprintf(outfile,"%s.sav.pl",mpt->out_prefix);
    pop_cell_write_sav_layer_list(laylist,nlay,outfile,s->name,s->repno);
  }

  // Response Handover
  if (m->run_mode == 0){  // 0=normal
    pop_util_response_handover(m,r,laylist,nlay,mylogf);
  }

  // Interactive mode
  if (m->run_mode == 1){
    printf("*** MOD_POP_GET_RESPONSE Warning: Interactive Mode disabled\n");
    // WYETH 2019 commented out to remove OpenGL, X11 dependence
    //mod_gui_pop_main(laylist,nlay,mpt,m,s,r);
  }

}
/**************************************-**************************************/
/*                                                                           */
/*                                MOD_POP_RUN                                */
/*                                                                           */
/*****************************************************************************/
void mod_pop_run(m,s,r,action)
     struct model_struct *m;     // Model params
     struct stim_struct *s;      // Stimulus params
     struct response_struct *r;  // Response params
     int action;                 // -1-cleanup, 0-prep, 1-run, 10-mar, 20-rsgen
{
  myid = mylog_set_id_log(m->process_id,&mylogf);  // Do this first

  mylog(mylogf,"  MOD_POP_RUN\n");

  if (action == 0){

    // WYETH - TRY MOVING
    // WYETH - TRY MOVING  mod_pop_prep first, and then use 'retflag'
    // WYETH - TRY MOVING

    if (onode_child_get_unique(m->o,"retina0") != NULL){
      model_mesh_rgc_01_prep(m,s,r);
    }else{
      model_dog_01_prep_o(m,r);
      model_dog_01_prep_ft();
    }
    mod_pop_prep(m,r,s); // Send 's' for GUI
  }else if (action == 1){
    mod_pop_get_response(m,s,r);
  }else if (action == -1){
    if (onode_child_get_unique(m->o,"retina0") != NULL){
      model_mesh_rgc_01_done();
    }else{
      model_dog_01_done(m->ppl); // Doesn't use argument
    }
  }else if (action == 10){
    if (onode_child_get_unique(m->o,"retina0") != NULL){
      model_mesh_rgc_01_prep(m,s,r);
    }
    mod_pop_prep(m,r,s); // r,s are NULL, write MAR file
  }else if (action == 20){
    mod_pop_prep(m,r,s); // s is NULL, write .RSP file
  }
}
