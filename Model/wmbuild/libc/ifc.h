/*****************************************************************************/
/*                                                                           */
/*  ifc.h                                                                    */
/*  Wyeth Bair                                                               */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*  input_pair                                                               */
/*                                                                           */
/*****************************************************************************/
struct input_pair{

  struct pop_mech *msi;    // Mechanism of synaptic interaction, "NULL" if none

  int in_index;            // Input index of original input

  struct pop_layer *pop0;  // First part of pair
  struct pop_layer *pop1;  // Second part of pair

  // WYETH - z0 doesn't belong here
  //char *z0;                // z-layer specification, first pop
  char *z1;                // z-layer specification, second pop

  char *distrib;           // Type of pairing rule

  // For jittering the x,y RF location of the paired input
  float rf_offset_dir;          // direction of offset
  float rf_offset_dist;         // distance of offset
  float rf_offset_sd_par;       //
  float rf_offset_sd_orth;      //

  float rf_ph_shift;            // phase shift (deg)
  float rf_ph_dev;              // phase deviation (deg)

  int   pair_seed;          // Randomization seed
  int   tseed;              // Temporary seed, will vary as matches are made
};
/*****************************************************************************/
/*                                                                           */
/*  nmda_param                                                               */
/*                                                                           */
/*****************************************************************************/
struct nmda_param{
  int    nmda_ppmv;  // points per mV in NMDA table
  int    nmda_v0;    // Starting value of table in mV
  float  nmda_v0f;   // Starting value of table in mV (float)
  int    nmda_n;     // Number of points in table
  float *nmda_gv;    // Relative conductance vs mV

  float nmda_frac;   // Relative conductance vs mV

  float nmda_alpha_1;  //
  float nmda_amp_1;    // This is set automatically, NOT in paramfile
  float nmda_alpha_2;  //
  float nmda_amp_2;    //
  float nmda_alpha_3;  //
  float nmda_amp_3;    //

  float nmda1_aetau;   // Precomputed constants for use in 'pop_input_02'
  float nmda1_tau2inv; //
  float nmda1_totau;   //

  float nmda2_aetau;   // Same as above
  float nmda2_tau2inv; //
  float nmda2_totau;   //

  float nmda3_aetau;   // Same as above
  float nmda3_tau2inv; //
  float nmda3_totau;   //

  // Allow separate scaling for LGN
  float lgn_nmda_frac; // Relative conductance vs mV
  float lgn_nmda_amp_1;  // This is set automatically, NOT in paramfile
  float lgn_nmda_amp_2;  //
  float lgn_nmda_amp_3;  //
};
/*****************************************************************************/
/*                                                                           */
/*  IFC = Integrate and Fire, Conductance                                    */
/*                                                                           */
/*****************************************************************************/
struct ifc_param{
  float v_spike;
  float tau_r_ad;      // (ms)
  float tau_f_ad;      // (ms)
  float gbar_ad;
  float v_reset_x;
  float v_ex;
  float v_in;
  float v_ad;
  float v_th_x;
  float v_leak_x;
  float g_leak_x;
  float c_x;
  float trefr_x;        // (ms)
  float trefr_x_s;      // (s) trefr_x/1000.0
  float trefr_x_sd;     // (ms) for rectified Gaussian noise

  float g_tran;         // Transfer resistance, = 1/rt

  float gx_scale;
  float gx_bias;
  float gi_scale;
  float gi_bias;
  char *gx_noise;       // Type of noise
  float gx_noise_mu;    // mean
  float gx_noise_sd;    // SD
  float gx_noise_tsd;   // SD for temporal filtering
  char *gi_noise;       // Type of noise
  float gi_noise_mu;    // mean
  float gi_noise_sd;    // SD
  float gi_noise_tsd;   // SD for temporal filtering
  int grect;            // 1/2 wave rectify gx, gi after adding noise, def 0
  int dump;
  int clip_spike;
};
/*****************************************************************************/
/*                                                                           */
/*  Poisson paramters                                                        */
/*                                                                           */
/*****************************************************************************/
struct poiss_param{
  char *style;        // "default" subtract scaled g_in & g_ad from g_ex
  float scale_in;     // Weight applied to g_in
  float scale_ad;     // Weight applied to g_ad
  float offset0;      // Added before scaling
  float scale;        // scale the combined conductance
  float offset1;      // added after scaling
};
/*****************************************************************************/
/*                                                                           */
/*  List of names.                                                           */
/*                                                                           */
/*****************************************************************************/
struct pop_name_list{
  int n;
  char **name;
};
/*****************************************************************************/
/*                                                                           */
/*  POP_CIRCBUF - circular buffer, i.e., storage that wraps around.          */
/*                                                                           */
/*****************************************************************************/
struct pop_circbuf{
  int n;          // length of buffer
  float *d;       // data storage [n]
  int i;          // index of origin  (storage wraps around)
  double t;       // time at origin
  double dt;      // time change between consecutive bins

  int maskn;      // Mask duration associated w/ buffer
};
/*****************************************************************************/
/*                                                                           */
/*  pop_mech                                                                 */
/*                                                                           */
/*****************************************************************************/
struct pop_mech{
  char *name;
  char *type;

  // PSG_ALPHA
  float tau;             // Time constant
  float amp;             // Amplitude

  // PSG_DOE
  float tau_r;           // Time constant rise
  float tau_f;           // Time constant fall

  // Synaptic depression
  float sdf;             // Time constant rise
  float sdtau;           // Time constant fall

  // Synaptic interaction
  int wtn;        // Length of circ buffer
  float wdt;      // Time resolution of buffer (msec)
  int mask1n;     // Length of mask
  float *mask1;   // [mask1n] mask for synaptic interactions
  int si_comp;    // 1='prod', 2='prod1m'  *** THIS CAN BE DONE w/ CONSTS
                  // For "mult" or "divide" 1=mult, 2=divide

  // "mult", "divide", and LGN pair, "multiply" or "divide"
  float normv;    // normalization value
  float a,b;      // additive and multiplicative constants for "divide"

  struct onode *o;  // The onode from which this mech was created
};
/*****************************************************************************/
/*                                                                           */
/*  POP_SI_SAV - data to save and return for responses                       */
/*                                                                           */
/*****************************************************************************/
struct pop_si_sav{
  int n;         // number of saved traces
  char **name;   // names [n]
  int *nmax;     // Max length of saved data [n]
  int *cnt;      // Length of saved data [n]
  float *samp;   // Sampling of saved data [n]
  float **d;     // [n][cnt[]] data values
};
/*****************************************************************************/
/*                                                                           */
/*  POP_SI_01 - synaptic interaction                                         */
/*                                                                           */
/*****************************************************************************/
struct pop_si001{
  struct pop_mech *msi;     // Pointer to SI definition
  struct pop_circbuf *wt;   // Weight vs. time
  struct pop_si_sav *sisv;  // Saved data to return for responses
};
/*****************************************************************************/
/*                                                                           */
/*  POP_SI_02 - synaptic interaction DUMMY                                   */
/*                                                                           */
/*****************************************************************************/
struct pop_si002{
  struct pop_mech *msi;      // Pointer to SI definition
  struct pop_circbuf *wt1;   // Weight vs. time
  struct pop_circbuf *wt2;   // Weight vs. time
  struct pop_si_sav *sisv;   // Saved data STRUCTURE to return for responses
};
/*****************************************************************************/
/*                                                                           */
/*  SI_UNION - one of many possible synaptic interaction types.              */
/*                                                                           */
/*****************************************************************************/
union si_union{
  struct pop_si001 *s001;
  struct pop_si002 *s002;
};
/*****************************************************************************/
/*                                                                           */
/*  POP_SI - synaptic interaction                                            */
/*                                                                           */
/*****************************************************************************/
struct pop_si{
  int ci;               // Cell index, for multi-presyn-cell interactions.

  // WYETH-NEW-MULT
  int syn_code;         // Identifies the type of synaptic interaction
                        // for the purpose of fast post-syn readout
                        //        0 - no post-syn readout
                        //   100002 - multiply two buffers


  int siui;             // SI union index
  union si_union *siu;  // One of several synaptic interactions,
			// shared by all cells involved in the
			//  interaction. Only the cell with 'ci' = 0
			//  should create or free the 'siu' storage.
};
/*****************************************************************************/
/*                                                                           */
/*  POP_SYN - synapse parameters                                             */
/*                                                                           */
/*****************************************************************************/
struct pop_syn{
  // CONNECTIONS
  struct pop_cell *post;     // Post-synaptic cell
  struct pop_cell *pre;      // Pre-synaptic cell

  // WYETH - Should probably put SD synaptic depression in synapse
  // rather than in cell.

  int stype;                 // Type of synapse
			     //  1 - ex1 - fast excitatory, e.g. non-NMDA
			     //  2 - in1 - fast inhibitory, e.g. GABA_A
			     //  11 - pre-computed (mesh input, like lgn)
			     //  ...
  float w;                   // Weight, from 0.0 to 1.0
  float tlast;               // Most recent spike time
  float alast;               // Amp at last spike time (synaptic depression)
  float tdelay;              // e.g., distance-related delay (sec)


  struct pop_si *si;         // Synaptic interaction, NULL-none

  // This index defines the input that created this synapse
  short inindx;              // Input Index, into pop_layer 'inlist'
                             // -1 - created by some other means
                             // -2 - created by customize

  struct pop_syn *pre_prev;  // Previous in PRE-synaptic list
  struct pop_syn *pre_next;  // Next in PRE-synaptic list
  struct pop_syn *post_prev; // Previous in POST-synaptic list
  struct pop_syn *post_next; // Next in POST-synaptic list
};
/*****************************************************************************/
/*                                                                           */
/*  POP_LGN_IN - special for LGN inputs                                      */
/*                                                                           */
/*****************************************************************************/
struct pop_lgn_in{
  struct pop_layer *lay;    // pointer to LGN pop layer (has name and xn,yn)
  struct pop_mech *mech;    // pointer to 'mech' for post-syn receptor
  int cn0;                  // Number of OFF inputs
  int *cx0,*cy0;            // OFF input coordinates [cn0]
  int cn1;                  // Number of ON inputs
  int *cx1,*cy1;            // ON input coordinates [cn1]


  /*********************/
  /*********************/
  // WYETH LGNW *** NOTE - Must check how this is handled when synapses are
  // read from a file ****** 2013 Jul
  float gw;                 // Global weight [1.0 by default]


  // TO BE ADDED
  //int cn0r;                 // RIGHT EYE Number of OFF inputs
  //int *cx0r,*cy0r;          // RIGHT EYE OFF input coordinates [cn0]
  //int cn1r;                 // RIGHT EYE Number of ON inputs
  //int *cx1r,*cy1r;          // RIGHT EYE ON input coordinates [cn1]

  int nposs;                // Total possible inputs

  struct pop_lgn_in *cpair; // Paired connections

};
/*****************************************************************************/
/*                                                                           */
/*  Parameters for evolving the state of one cell in time                    */
/*                                                                           */
/*  Initialized by:  'pop_cell_init_cell'                                    */
/*                                                                           */
/*****************************************************************************/
struct pop_cell{
  char *name;            // Cell ID - EXTRANEOUS, can be computed?
  char *subclass;        // Cell sub-class, e.g., "ds01"
  struct pop_layer *pl;  // Layer for this cell
  int layx,layy,layz;    // position in layer
  float cx,cy;           // Cortical location, um, (0,0) is stimulus center
                         // or, irregular position (0,0) is lower left
  float rfx,rfy;         // Center of RF on LGN/stim grid (for ocdom)

  // Cell attributes, e.g.
  //
  //  'ori'           Ori/dir assigned to cell, or derived for cell
  //  'ori_cv'        Orientation Circular Variance for cell
  //  'phase'         Phase assigned to cell, 0..360
  //  'sz1'           Size ?
  //  'ocdom'         Ocular dominance -1.0 to 1.0 = L-R
  //  'conetype'      0.0 1.0 2.0  (also stored in layer 'irr_id')
  //
  float *attrib_f;      // [pl->attrib_fn]  Attribute values

  // SYNAPTIC CONNECTIONS
  int nout;                  // Number of cells targeted
  struct pop_syn *out;       // List of "nout" synapses
  int nin;                   // Number of cells sending input
  struct pop_syn *in;        // List of "nin" synapses

  // Indicate any special processing to advance the cell's conductance input
  int syn_proc_code;         // 0-no special processing
                             // 1-yes special processing

  float input_ex_sdf;        // Synaptic depression for excitatory input
  float input_ex_sdtau;      // Synaptic depression for excitatory input

  // LGN inputs
  int lgn_n;                 // Number of LGN populations giving input
  struct pop_lgn_in **lgnin; // LGN inputs, array of pointers [n_lgn_in]
                             // *** THIS LIST remains NULL if irreg LGN input

  int lgn_mflag;             // 0-p, 1-m;  WYETH TEMPORARY


  struct   ifc_param *cifcp;    // Parameters for IFC model (old name was p)
  struct poiss_param *cpoissp;  // Poisson params, for <spike_gen> poisson


  // EXCITATORY INPUT SPIKES FOR NEAR FUTURE
  struct pop_circbuf *ex1s;       // Upcoming spikes, like a PSTH
  struct pop_circbuf *in1s;       // Upcoming spikes, like a PSTH
  struct pop_circbuf *ex1g;       // Recently computed g
  struct pop_circbuf *ex1g_nmda;  // Recently computed gNMDA
  struct pop_circbuf *in1g;       // Recently computed g

  // WYETH - NEW - 2010 Jun 19
  struct pop_circbuf *ad1g;       // Recently computed gad

  double gex1t;     // time of "gex1" value
  float gex1;       // conductance value
  float gdex1;      // derivative of "gex1"
  float ex1tau;     // WYETH - this should come from IFC params!!!
  float ex1amp;     // WYETH - this should come from IFC params!!!

  float gex1_nmda1,gex1_nmda2,gex1_nmda3;     // NMDA
  float gdex1_nmda1,gdex1_nmda2,gdex1_nmda3;

  int savg;         // Whether to store dynamic conductance input
  float *gtex1;     // Save "gex1"
  float *gtex1_nmda;// Save g NMDA
  float *gtin1;     // Save "gin1"
  float *gtex1t;    // Save time assoc'd w/ "gex1"
  int gtex1n;       // number currently stored
  int gtex1max;     // max storage

  int    sava;      // Whether to save the adaptation conductance
  float *gta0;      // Adapting input for full time [gta_max]
  int    gta_n;     // Current value stored.
  int    gta_max;   // Length of gta0

  // INHIBITORY INPUT
  float gin1;       // conductance value
  float gdin1;      // derivative of "gin1"
  float in1tau;     // WYETH - this should come from IFC params!!!
  float in1amp;     // WYETH - this should come from IFC params!!!

  // PRE-COMPUTED CONDUCTANCE INPUTS
  int n0;           // length of "gtx0" and "gti0"
  float samp0;      // Samples/sec for "gtx0" and "gti0" and "ggain"
  float *gtx0;      // Known excitatory input for full time
  float *gtn0;      // Known NMDA excitatory input (to be scaled by V later)
  float *gti0;      // Known inhibitory input for full time
  float *ggain;     // Gain to multiply 'gtx0' from LGN

  // ADAPTATION PULSE STORAGE
  int nad;            // length of "gta_pulse"
  float *gta_pulse;   // Single adaptation g pulse

  //
  //  WYETH - could use sub-structures to modularize storage size
  //

  // BACKGROUND SPIKES
  float bg_x_rate;    // Background rate to add spikes to 'ex1s'
  float bg_i_rate;    // Background rate to add spikes to 'in1s'
  float bg_x_amp;     // Synaptic weight
  float bg_i_amp;     // Synaptic weight
  float *bg_x_s;      // spike times (use 'samp0' for units)
  float *bg_i_s;      // spike times (use 'samp0' for units)
  int bg_x_n;         // Number of spikes (-1 for none generated yet)
  int bg_i_n;         // Number of spikes (-1 for none generated yet)
  int bg_x_k;         // Index of first un-used spike
  int bg_i_k;         // Index of first un-used spike
  struct pop_cell *bg_x_corr_c; // Pointer to cell for correlated BG
  float            bg_x_corr_f; // Fraction of spikes to share
  struct pop_cell *bg_i_corr_c; // Pointer to cell for correlated BG
  float            bg_i_corr_f; // Fraction of spikes to share

  //
  //  WYETH - could use sub-structures to modularize storage size
  //

  // RUNGE-KUTTA PARAMETERS
  double x;        // Current time (sec)
  float y;         // Current value of Vm (mV) [for soma]
  float yd;        // Current value of Vm (mV) [for dendr] 2-comp 
  float ystart;    //
  float yscal;     // Scaling used to monitor accuracy...
  float h;         // Step size
  float h1;        // initial step size
  float hmin;      // Smallest step size
  float hmax;      // Largest allowable step size, to limit V_thr overshoot
  int nbad,nok;    // Number of good and bad step sizes
  float eps;       // Desired accuracy
  float reft;      // Time at which current refractory period ends

  // LONG-TERM STORAGE OF Vm
  int savv;        // Whether to store Vm
  char *savvname;  // Name to associate with output
  float dxsav;     // Save values no more often that this, storage interval
  float xsav;      // x-value of last point saved
  int vnmax;       // length of "vm" and "vmt"
  int vn;          // Number of values stored so far 
  float *vm;       // ith Membrane voltage (mV)
  float *vmt;      // Time at which ith voltage was recorded

  // same as above, but for DENDRITE 2-comp
  int savd;        // Whether to store Vm, dendr, must also have savv=1
  char *savdname;  // Name to associate with output
  float xsavd;     // x-value of last point saved
  int vnd;         // Number of values stored so far 
  float *vmd;      // ith Membrane voltage (mV)

  // CURRENT PRE-SYN SPIKES
  int sdt;     // 1/sampling of "sx" and "si", perhaps 0.2 to 0.5 msec
  int sn;      // length of "sx" and "si", 100 x 0.5 "sdt" allows up to 5
	       //   delay for a spike to travel from pre-syn to post-syn
  int sxp,sip; // pointers to current time in "sx" and "si"
  float *sx;   // current excitatory pre-syn spikes [time]
  float *si;   // current inhibitory pre-syn spikes [time]

  // HISTORY OF POST-SYNAPTIC SPIKES
  int savs;        // Flag to save spikes
  char *savsname;  // Name to associate with output
  int maxscnt;     // Max spikes per trial, can increase on fly
  float *s;        // Spike times
  int ns;          // Number of spikes


  // SAVE ARBITRARY DATA
  int    *csav_cnt;     // [pl->csav_n] Number of data values in each trace
  float **csav_f;       // [pl->csav_n][csav_cnt[]] Data ID:  pl->csav_datid

};
/*****************************************************************************/
/*                                                                           */
/*  pop_icon                                                                 */
/*                                                                           */
/*****************************************************************************/
struct pop_icon{
  char *shape;      // Shape
  int nside;        // number of sides, for shape = "polygon"
  float r,g,b;      // color 0
  float r1,g1,b1;   // color 1
  struct onode *o;  // Node of type 'icon', Added for .mar writing, Jun 2009
};
/*****************************************************************************/
/*                                                                           */
/*  pop_layer                                                                */
/*                                                                           */
/*****************************************************************************/
struct pop_layer{
  char *name;       // Layer ID
  int runflag;      // 0-don't run, 1-run

  int lgn_in_flag;  // 0-no LGN input, 1-LGN input.  Essentially, this means
                    // that spikes must be retrieved to compute pre-computed
                    // inputs.  Could be DOG or MESH models.
  int lgn_in_n;     // Number of LGN inputs, 0..

  int lgn_gain;       // 0-no gain signal, 1-gain signal,
                      //   If LGN layer, compute gain sign, else use gain
  float lgn_gain_tau; // Time const for smoothing gain signal (s)
  float lgn_gain_ca;  // Additive constant
  float lgn_gain_cb;  // Multiplicative constant

  char *laytype;    // Layer type, "default" if not spec'd
                    //   "default"      - none specified
                    //   "lgn"          - LGN for mod_dog
                    //   "retina0_gc0"  - mod_mesh
                    //   "irregular"    - irregular cell placement
                    //   "virtual"      - No spike generation

  struct pop_icon *icon;  // GUI icon configuration

  struct pop_area *area;  // Area

  int ninlist;            // Number of inputs in 'inlist'
  struct onode **inlist;  // [nin] index values are referenced from
                          // pop_syn 'inindx' input index
                          // *** WARNING, if this list is reordered,
                          // 'inindx' values would need to be changed

  //
  //  GEOMETRY.  Conceptually, there is an underlying grid that has
  //  umx,umy spacing on the cortex. Currently, this is the
  //  stimulus/lgn grid
  //
  char *geomt;     // "irregular" or "default"
  int xn,yn,zn;    // dimensions of layer; 'xn' is 'irr_n' if "irregular"
  float x0,y0;     // Origin of layer wrt overall 'xn', 'yn' stimgrid
  float xf,yf;     // Scaling of coords wrt stim xn yn grid (stimgrid/cell)
  // WYETH WARNING:  These values are not used?? - See area->umx
  float umx,umy;   // Microns on cortex between cells in (cortical?) grid.
		   //   Ultimately, this should be shared among all layers in
		   //   the same "area"
  float oddxoff;   // Odd row, x-offset, for checkerboard or hexagonal,
                   //   in grid units
  int irr_n;       // Number of irregular coords, = 'xn' if "irregular"
  float *irr_x;    // [irr_n] x-coords
  float *irr_y;    // [irr_n] y-coords
  int *irr_id;     // [irr_n] identifier
  float irr_dens;  // units per deg^2

  struct dll_slist *save_name_list;  // Names of values to save in monitor

  int  layseed_n;    // Number of seeds for this layer
  int *layseed;      // Seed values [layseed_n]

  int nmech;               // Number of "mech" objects
  struct pop_mech **mech;  // List of pointers [nmech]

  struct pop_cell ***c;       // Cells in layer [xn][yn][zn]

  int attrib_fn;         // Number of float attribs in each cell
  char **attrib_fname;   // Attribute names


  struct ifc_param *ifcp;     // IFC params, may be shared by all cells
  struct poiss_param *poissp; // Poisson params, may be shared by all cells
  struct nmda_param *nmda_p;  // NMDA parameters, NULL if no NMDA

  //  Saving arbitrary data
  int    csav_n;              // Number of saved data traces
  char **csav_datid;          // Data ID (datid) for each trace, the index
                              // here matches that in c->csav_i
  char **csav_source;         // Source of data, e.g., 'retina0'
  float *csav_samp;           // Samples per second

  // LGN layers
  // Pre-computed LGN spikes, NOTE, 'zn' IS FIRST DIMENSION
  // zn: 0-off 1-on 2-off_r 3-on-r
  int   ***cnt;               // spike count ***[zn]***[xn][yn]
  float ****s;                // spike times ***[zn]***[xn][yn][cnt]
  float ***rl;                // LGN response, left
  float ***rr;                // LGN response, right
  float ***rgl;               // Gain signals (not tensors, 0 indexing)
  float ***rgr;               // Gain signals (not tensors, 0 indexing)

};
/*****************************************************************************/
/*                                                                           */
/*  pop_map                                                                  */
/*                                                                           */
/*****************************************************************************/
struct pop_map{
  struct onode *o;    // Node of type 'map'
  int xn,yn,zn;       // Size of map
  float **map2;       // Map data [xn][yn]
  // float ***map3d;  FUTURE
};
/*****************************************************************************/
/*                                                                           */
/*  pop_area                                                                 */
/*                                                                           */
/*****************************************************************************/
struct pop_area{
  char *name;         // lgn, v1, etc
  float x0,y0;        // origin wrt stim grid
  float xf,yf;        // expansion factor wrt stim grid
  int xn,yn;          // width, height of area
  float umx,umy;      // um per area grid unit (for now Area = Stim grid)
  float sscale;       // stim grid sscale (for access from cells via layer)
};
/*****************************************************************************/
/*                                                                           */
/*  pop_top                                                                  */
/*                                                                           */
/*****************************************************************************/
struct pop_top{
  char *name;               // Name of model
  char *logf;               // Log file name, typically copy of mylogf

  int xn,yn;                // Grid size (currently = LGN = STIM size
  int tn;                   // Stimulus duration
  float tscale;             // Stimulus time scale
  float sscale;             // Stimulus spatial scale
  float x1,x2;              // start and stop time (s) to run simulation
  char *out_prefix;         // Output file prefix
  int monflag;              // Monitor

  int retflag;              // 0-DOG, 100-retina0, 101-retina0 irreg, -1-none
  int binoc;                // 0-monocular left eye, 1-binocular l+r

  int narea;                // Number of areas
  struct pop_area **area;   // Pointers to areas [narea]

  int nlay;                 // Number of layers
  struct pop_layer **lay;   // Pointers to layers [nlay]

  int nmap;                 // Number of maps
  struct pop_map **map;     // Pointers to maps [nmap]

  char *conn_rw;            // "none", "read", "write"
  char *conn_read_dir;      // Read connections from this directory, or NULL
  char *conn_write_dir;     // Read connections from this directory, or NULL
};

/***
 ***  SETS OF STORED CONNECTIONS
 ***/
/*****************************************************************************/
/*                                                                           */
/*  pop_conn_el                                                              */
/*                                                                           */
/*****************************************************************************/
struct pconn_el{
  int nsyn;  /* Number of synapses */
  int *xi;   /* Coords and weight */
  int *yi;
  int *zi;
  float *w;
};
/*****************************************************************************/
/*                                                                           */
/*  pop_conn                                                                 */
/*                                                                           */
/*****************************************************************************/
struct pconn{
  int ftype;                 // File type
  char *name1;               // Pre-syn name
  char *name2;               // Post-syn name
  char *syntype;             // Synapse type
  int xn,yn,zn;              // Size of pre-syn layer
  struct pconn_el ***c;      // [xn][yn][zn]
};
