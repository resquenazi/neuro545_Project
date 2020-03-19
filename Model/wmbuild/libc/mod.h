/*****************************************************************************/
/*                                                                           */
/*  mod.h                                                                    */
/*  Wyeth Bair                                                               */
/*                                                                           */
/*  Data structures used by 'wm.c' and by 'mod_xxx_util.c' functions.        */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/

/**************************************-**************************************/
/*                                                                           */
/*                             STIM_VLINK_STRUCT                             */
/*                                                                           */
/*  WYETH - for the future, to allow vertical tables in stim file, but       */
/*  this will move away from the param-pair-list model.                      */
/*                                                                           */
/*****************************************************************************/
struct stim_vlink_struct{ // For visual stimulus
  int npar;     // Number of parameters
  char *pname; 
  int nval;     // Number of values
  char **pval;  // Values [npar][nval]
};
/**************************************-**************************************/
/*                                                                           */
/*                                VAR_PAR_REC                                */
/*                                                                           */
/*  Uses:                                                                    */
/*   -one var par: n=1, uval contains all values                             */
/*   -set of n varlinks                                                      */
/*   -set of specific values for a single control stimulus, nval=1           */
/*                                                                           */
/*****************************************************************************/
struct var_par_rec{
  int n;              // Number of linked params
  char **name;        // [n] Names of linked params
  char *ptype;        // [n] Type of param, 'f', 'i', 'c'
  int nval;           // Number of (unique) values in list
  char ***val;        // [n][nval] Parameter values
};
/**************************************-**************************************/
/*                                                                           */
/*                               VAR_PAR_STRUCT                              */
/*                                                                           */
/*****************************************************************************/
struct var_par_struct{

  //
  //  Model params that vary
  //
  int m_var_n;
  int m_one_n;
  struct var_par_rec **m_var;  // [m_var_n] Vary all against each other
  struct var_par_rec **m_one;  // [m_one_n] Additional one-off stimuli

  int     m_n;       // Total number of model configurations to run
  int     m_nvar;    // Number of model params to vary
  char  **m_name;    // [m_nvar] Pointer to param name (in var_par_rec)
  char ***m_vval;    // [m_n][m_nvar] table of parameter values

  //
  //  Lookup for varpars, <varpar_lookup>
  //
  int     vplk_n;        // Number of <varpar_lookup> objects
  char  **vplk_name;     // [vplk_n] name
  int    *vplk_nval;     // [vplk_n] number of variable values
  char ***vplk_key;      // [vplk_n][vplk_nval[]] Key names for values
  char ***vplk_val;      // [vplk_n][vplk_nval[]] Values for this variable


  // 2018 March
  // - Adding 'varpar_lookup' to allow var params with a list of values
  //   for example, for varying the MT weights

  /*** WYETH NEXT STEP ***/
  // - Write code to make a list (table) of all the var mod param assignments
  //     and count up the total models to run
  // 
  // - Create a command in 'mm' to tell the model to clean-up / free.
  // - Create a command in 'mm' to tell the model to re-prep itself.
  // - Make sure that the model can clean up all of its storage.
  //

  //
  //  Stim params that vary
  //
  /*** FUTURE
  int s_var_n;
  int s_one_n;
  struct var_par_rec **s_var;  // [m_var_n] Vary all against each othear
  struct var_par_rec **s_one;  // [m_one_n] Additional one-off stimuli
  ***/
};
/**************************************-**************************************/
/*                                                                           */
/*                                 STIM_STRUCT                               */
/*                                                                           */
/*****************************************************************************/
struct stim_struct{
  char paramfile[SLEN];        // Parameter file name
  struct param_pair_list *ppl; // List of parameter names and values

  //
  //  2017 Aug 11 - WYETH: I want to allow the AlexNet model to work with
  //                an RGB stimulus with floating point accuracy, and I am
  //                considering the following additional fields here:
  float ****drgb_l;  // [3][x][y][t] RGB image, left eye
  float ****drgb_r;  // [3][x][y][t] RGB image, right eye

  int ***dc;                   // 3c stimulus (packed RGB) color
  int ***dc_r;                 // 3c stimulus, right eye, color

  float ***d;                  // 3d stimulus data for current stim
  float ***d_r;                // right eye stimulus
  float *d1;                   // 1d stimulus data for current stim

  int repno;                   // Repeat number
  int stimno;                  // Stimulus number, among UNIQUE STIMULI
			       //   when this changes, new stim data is needed
  char *name;                  // Stimulus name

  // NOTE: 'vval' is stored in a fixed sequential order, which used to match
  // the execution order until blockwise randomization (BR) was allowed.
  // In the new case of BR, response storage is in execution order, and so
  // we needed the new index 'tribyx' to relate each executed trial to its
  // 'tri' index in 'vval'.

  int nvar,ncon,ntr;           // number of params, number of trials
  char ***vval,**vname,*vtype; // Var params [ntr][nvar] or [nvar]
  char ***val;                 // Unique values of var params [nvar][vcnt[]]
  int *vcnt;                   // Unique values of var params [nvar]
  int nvl;                     // Number of linked values
  int nvp;                     // Number of linked pairs
                               //   This should be more general

  char **cval,**cname,*ctype;  // Const params [ncon]

  int *tribyx;              // Trial index, stored by execution order [ntr]
                            //   This is the index into 'vval', it stores
                            //   'r->tri' in exec. order, for blockwise rand
};
/*****************************************************************************/
/*                                                                           */
/*                             MODEL_DATA_STRUCT                             */
/*                                                                           */
/*  2018 June.  This is meant to describe data that is available to request  */
/*  in the .rsp file.                                                        */
/*                                                                           */
/*****************************************************************************/
struct model_data_struct{
  char *pname;   // population name
  char *dataid;  // Data ID
  char  dtype;   // 'f' (float), 's' (spike)
  int   ndim;    // dimension of population, e.g., 3
  int  *dlist;   // [dim] length of each dimension
  void *dptr;    // Pointer to multi-dimensional data, or NULL
  char *davail;  // 'full' - all positions in 'dptr' are available
                 // 'part_nullsig' - elements in 'dptr' are NULL if no data
                 // 'part_unknow' - it is unknown which are available
                 // 'other' - ???
  char *stimdep; // 'none', 'binoc' - data available if stimulus is binocular
  char *ospec;   // Other specification, or NULL
                 //  'stimdep' - stimulus dependent as to whether it will exist
                 //  

  struct model_data_struct *prev;
  struct model_data_struct *next;
};
/*****************************************************************************/
/*                                                                           */
/*                                MODEL_STRUCT                               */
/*                                                                           */
/*****************************************************************************/
struct model_struct{
  char paramfile[SLEN];        // Parameter file name

  // FOR .mod
  struct param_pair_list *ppl; // List of parameter names and values

  // FOR .moo
  struct onode *o;             // o-node description of model, tree

  // For variable params
  int nrun;                    // Number of variable param configs
  int m_i;                     // Index of currently running config
  struct var_par_struct *vps;  // Variable parameters

  // WYETH added Dec 2018
  char **calib_pval;           // [nrun] Calibration parameter values, or NULL

  int process_id;              // -1-Single, 0-Master, >0-Slave#
  int process_n;               // Number of processes

  int run_mode;                // 0-normal
                               // 1-interactive
                               // 2-replay

                               // 8-indefinite continuing
                               // 9-indefinite finished

  int nseed;                   // Number of seeds in mseed list
  int *mseed;                  // [nseed] Seeds change model noise w/ repeats
                               //    Note: subscribed by 'r->tsi'

  float ***fpe,***fpo,***fne,***fno;     // ME filters, quadrature opponent
  float **fpe_s,**fpo_s,**fne_s,**fno_s; // For FT of filters

  float ***dog;                // For DOG model filter
  float **dog_s;               // For DOG model FT

  int action;                  // wm_slave_action flag

  char *marfile;               // Used for writing .mar file

  int mds_n;                      // Current number of records
  struct model_data_struct *mds;  // Doubly-linked list, or NULL
};
/*****************************************************************************/
/*                                                                           */
/*                               RESPONSE_STRUCT                             */
/*                                                                           */
/*****************************************************************************/
struct response_struct{
  char paramfile[SLEN];        // Parameter file name
  struct param_pair_list *ppl; // List of parameter names and values

  int tri;                     // Index for current trial, maybe random order
                               //   e.g. as index into 'vval'
  int tsi;                     // Trial sequence index, always sequential 0,1..
  int gtsi;                    // Global trial sequence index, uses m->m_mi
                               //   The trial index into 'f' or 's', which 
                               //   factors in m->m_mi.
                               //   WYETH - JULY 29, 2014
  char *tname;                 // Name for current trial, NULL if no var pars

  // Only used by Master Process in 'mm'
  int *dflag;                  // Done flag [ntr] ('ntr' in 'stim_struct')

  double ts0,tsn;               // Start time and duration (s)
  double ttref;                 // tref for trial (s)
                                //   NOTE - only used by 'mod_wav_util'?

  //
  //  Response description
  //
  //  Notes
  //   - Unfilled responses are indicated by -1 in 'cnt' or 'fcnt'
  //
  int n;                       // Number of response traces (per trial)
  int *rformat;                // [n] 0 - save_as_
			       //     1 - save_pop_unit_as_ 
			       //     2 - save_pop_grid_as_ 
			       //     3 - save_pop_layer_as_
			       //     4 - save_pop_all_spikes_  - FUTURE
			       //    10 - save_pop_syn_as_
  char **plname;               // [n] Name within model, or pop layer name
  char **datid;                // [n] Data ID, e.g., 'lgn_gx', 'spikes'
  int   *xi,*yi,*zi;           // [n] coordinates within pop (for format 1)
  char **nd_name;              // [n] Name to write to ndata
  float *samp;                 // [n] Response sampling
  char **rtype;                // [n] Data type "s"-spike, "f"-float
  int   *ri;                   // [n] Index within type, e.g. 0..ns-1, 0..nf-1

  char **plname1;              // post-syn unit layer
  int *xi1,*yi1,*zi1;          // post-syn unit coords

  // The above 'ri' could be used as a flag for completion

  // "s" - Spike times (float)
  int ns;                      // Number of spike responses/trial (channels)
  float ***s;                  // Spike train responses [ns][ti][cnt]
  int **cnt;                   // Spike count [ns][ti]  -1 MEANS UNFILLED

  // "f" - Continuous float
  int nf;                      // Number of float responses/trial (channels)
  float ***f;                  // Float responses [chan][ntr][cnt]
  int **fcnt;                  // Float count [chan][ntr], -1 MEANS UNFILLED
};
