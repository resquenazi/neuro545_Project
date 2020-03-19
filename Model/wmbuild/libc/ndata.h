/*****************************************************************************/
/*                                                                           */
/*  ndata.h                                                                  */
/*  Wyeth Bair                                                               */
/*                                                                           */
/*  This file contains data structures used for spike train processing.      */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*                                                                           */
/*  MULTI CHANNEL SPIKE DATA                                                 */
/*                                                                           */
/*  NOTES:                                                                   */
/*  - Channels, or records, should be referred to by name only.  They are    */
/*    not constrained to be in the same order in all trials, or to be        */
/*    present in all trials.                                                 */
/*  - A record name cannot be repeated within a trial block.                 */
/*                                                                           */
/*  RESERVED NAMES                                                           */
/*  - params:  MOO_<...> for command line model params as const params       */
/*  - params:  VAR_<...> for VAR param lists                                 */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/
/***

WYETH - SHOULD CONSIDER USING FIRST <int> AS VERSION, AND MAKING SOME CHANGES
  (*) 'tref' - could be long


FILE
  (int)1 (the number 1 written as an INT)
  File_class - typically experiment class, such as "cseries" or
              or "model01"
  N_Constant_Params - STIMULUS parameters constant for all Trials.
    nchars name nchars value --- Names and values stored as characters.
       ...
    nchars name nchars value --- Names and values stored as characters.
  N_Varied_Params - STIMULUS parameters not constant for all Trials.
    nchars name --- All stored as characters.
     ...
    nchars name --- All stored as characters.
  N_Trials - Number of trials in this file (typically number of trials).
TRIAL
  Trial_code - could be stimulus number (somewhat redundant with stim_code)
  Ref_time - could be real world time.
  N_Params - number of parameters to specify for this trial, including
             STIMULUS params, and other things such as time codes, response
             codes, abort codes, etc.
    nchars name nchars value --- All stored as characters.
     ...
    nchars name nchars value
  N_Records - number of records in this trial.  Could be 2 spike channels
              and one LFP chan.

REC_CODE nchar name stim_code sampling time_zero period n
         offset, offset, offset

  name - identifies the source of the points.  For example, "cell01",
         "cell02", "trigger01", "mua01".
  stim_code - identifies the class of stimulus, or experimental condition
              to which this belongs.

analog data:

  name stim

TYPES OF RECORDS:

Code Description
 0   point - point process, spikes, or any type of timing pulse.
 1   analog - regular sampling, floating point, for LFP data.
 2   point value - time and value representation, float values at int times.

WHEN A FILE IS READ:

(1) Check for common sampling on all point channels
***/
/*****************************************************************************/
/*                                                                           */
/*                                   NDATA                                   */
/*                                                                           */
/*****************************************************************************/
struct ndata_struct {
  char  *class;             // File class, i.e., "cseries" or "model".
  int    nconst;            // Number of params constant across all trials.
  int    nvar;              // Number of stimulus params varied.
  char **cname;             // Param names [npconst].
  char  *ctype;             // Constant parameter types: 'i', 'f', 'c'
  char **cval;              // Param values
  char **vname;             // Variable parameter names.
  char  *vtype;             // Variable parameter types: 'i', 'f', 'c'
  int     ntable;           // Number of event code tables
  struct ect_struct *table; // list of tables.
  int    ntrial;            // Number of trials.
  struct ndtrial_struct *t; // list of trials.
};
/*****************************************************************************/
/*                                                                           */
/*                                ECT_STRUCT                                 */
/*                                                                           */
/*  Event code table.  Useful for REX and EXP file formats.                  */
/*                                                                           */
/*****************************************************************************/
struct ect_struct {
  char  *tname; // Name of table
  int    tnum;  // Ordinal number of table, important for file storage.
  int    n;     // Number of codes in table.
  char **name;  // Code description.
  int   *num;   // Code number
};
/*****************************************************************************/
/*                                                                           */
/*                                NDTRIAL                                    */
/*                                                                           */
/*****************************************************************************/
struct ndtrial_struct {
  int    tcode;         // Used for stimulus ID, or stimno, see 'sr_nda'

  int    tref;          // Reference time, could be stim order, or real time
                        // Labr convention is TREF in MSEC
                        // wm convention, TREF in MSEC

  int    nparam;        // Number of params specified for this trial.
  char **pname;         // Param names [nparams].
  char **pval;          // Param values
  int    nrec;          // Number of records
  struct ndrec_struct *r; // pointer to records
  struct ndtrial_struct *next; // For linked lists of trials
  struct ndtrial_struct *prev; // For linked lists of trials
};
/*****************************************************************************/
/*                                                                           */
/*                                  RECORD                                   */
/*                                                                           */
/*  The first "p" value in a point record is the time from "t0" to the       */
/*  time of the first point.                                                 */
/*                                                                           */
/*  Types of records (rtype):                                                */
/*   0 - point, a list of event (spike) times in "p", for spikes, syncs.     */
/*   1 - float, a list of float values in "x", for eye position, LFP.        */
/*   2 - point-value, a list of float values "x" at times "p".               */
/*   3 - event-code, a list of codes "ec" at times "p".                      */
/*                                                                           */
/*  WYETH - perhaps add a 'table' record that holds a 2D table of string     */
/*  values and that points to a table format description stored in the       */
/*  header (like event code tables).  This would allow storage of arbitrary  */
/*  values in each trial, whose names and types were stored only once.       */
/*  Such could be used for Andrew Parker's saccade parameters.               */
/*                                                                           */
/*****************************************************************************/
struct ndrec_struct {
  int    rtype;        // Record type:  0--3
  char  *name;         // channel name
  int    rcode;        // USED BY OLD PROGRAMS, for stimulus ID?
  float  sampling;     // samples per second

  // LABR
  // t0 - start time of recording window
  // tn - duration of rec window (typically longer than stimulus on period)
  // trial 'tref' is global time of stimulus start (msec)
  // Thus, a spike time of 0 occurred at stimulus onset.

  int    t0;           // time origin of data (in sampling units)
		       //  relative to trial "tref"
  int    tn;           // time period of data (in sampling units)
  int    n;            // number of data points
  int   *p;            // point data [n]
  float *x;            // analog data [tn]
  int   *ec;           // event code number [n]
  struct ect_struct *ect; // pointer to event code table, used with "p"
                          //   or NULL if no table is given
};
/*****************************************************************************/
/*                                                                           */
/*                                   GROUP                                   */
/*                                                                           */
/*  This structure is used to group trials according to parameter values.    */
/*  Used when the number of trials is knowable ahead of time.                */
/*                                                                           */
/*****************************************************************************/
struct ndgroup_struct{ // For grouping trials
  int n;          // Number of groups
  char **name;    // Group names, could be parameter names [n]
  char **value;   // Group values, could be parameter values [n]
  int *cnt;       // Number of trials in each group [n]
  int **tnum;     // Trial numbers for each group [n][cnt[n]]
};
/*****************************************************************************/
/*                                                                           */
/*                                GROUP_NODE                                 */
/*                                                                           */
/*  A linked list of pointers to linked lists of trials.  Used when trials   */
/*  become available sequentially over time, and the total number of trials  */
/*  and groups may not be knowable at the outset.                            */
/*                                                                           */
/*****************************************************************************/
struct ndgnode{
  int n;                       // number of trials in group at this node
  struct ndtrial_struct *head; // first trial in list
  struct ndtrial_struct *tail; // last trial in list
  struct ndgnode *next;        // next group
  struct ndgnode *prev;        // previous group
};
/*****************************************************************************/
/*                                                                           */
/*                                GROUP_NODE_LIST                            */
/*                                                                           */
/*  A linked list of pointers to linked lists of trials.  Used when trials   */
/*  become available sequentially over time, and the total number of trials  */
/*  and groups may not be knowable at the outset.                            */
/*                                                                           */
/*****************************************************************************/
struct ndglist{
  int n;                 /* Number of groups. */
  struct ndgnode *head;  /* Pointer to first group node */
  struct ndgnode *tail;  /* Pointer to last group node */
};
/*****************************************************************************/
/*                                                                           */
/*                               RECORD_INDEX                                */
/*                                                                           */
/*  This structure is needed because records are identified by name only,    */
/*  and may be ordered arbitrarily in a trial block.                         */
/*                                                                           */
/*****************************************************************************/
struct ndrx_struct{ /*** An index to all records with the same name. ***/
  char *name;     // Name of the record
  int n;          // Number of records in this index
  int *tnum;      // Trial number of each record
  struct ndrec_struct **r; // Pointers to each record
};
/*****************************************************************************/
/*                                                                           */
/*                                    OLD                                    */
/*                                CONDITIONAL                                */
/*                                                                           */
/*  If "chan" < 0, no conditions are checked, all trials are passed.         */
/*  Set "min_spikes" = "max_spikes" = -1 for no spike count condition.       */
/*                                                                           */
/*****************************************************************************/
struct ndcond_struct_OLD{ /*** For conditions on spike trains. ***/
  char *chan;        /* Channel name for spike train condition. */
  int start;         /* Start time of window for spike count condition,
			relative to rec "start"	given in sampling units. */
  int period;        /* Period of window for spike count condition. */
  int min_spikes;    /* Minimum number of spikes for inclusion of window */
  int max_spikes;    /* Maximum number of spikes for inclusion of window */
  int ncond;         /* Number of conditions used */
  int *flag;         /* 0=unused, 1=min, 2=max, 3=min and max, 4=pval list */
  char **pname;      /* Parameter name for condition [ncond] */
  float *min,*max;   /* Minimum and maximum parameter values [ncond] */
  int *nval;         /* Number of values for condition [ncond] */
  char ***pval;      /* Set of values [ncond][nval] */
};
/*****************************************************************************/
/*                                                                           */
/*                                   INDEX                                   */
/*                                                                           */
/*  This is freed by "free_ndata_index" in "ndata_util".  Be sure to update  */
/*  that routine if this structure is modified.                              */
/*                                                                           */
/*  After sorting, tnum[ndx[0]] ... tnum[ndx[nrec]] is the ordering.         */
/*                                                                           */
/*****************************************************************************/
struct ndindex_struct{ // For sorting trials
  int nrec;       // Number of TRIALS to index.
  int n;          // Number of index parameters.
  char **pname;   // Names of index parameters.
  char *ptype;    // 'i', 'f', 'c' for int, float, category (character).
  int *nuval;     // Number of unique param values, only for 'c' params, [n]
  char ***uval;   // Unique values, used only for 'c' params, [n][nuval[n]]
  float **ndxval; // Index values [nrec][n]
  int *tnum;      // [nrec] Trial no. in "nd", monotonic, but may skip trials
  int *ndx;       // [nrec] Index pointer.
};
/*****************************************************************************/
/*                                                                           */
/*                                CONDITIONAL                                */
/*                                                                           */
/*  This condition structure holds a list of conditions which are used to    */
/*  determine which trials will be processed.                                */
/*                                                                           */
/*****************************************************************************/
struct ndcond_struct{ // For conditions on spike trains
  int n;         // Number of conditions
  int *ns;       // Number of strings in each condition [n]
  char ***slist; // Strings specifying conditions [n][ns]
};
/*****************************************************************************/
/*                                                                           */
/*                                 ANALYSIS                                  */
/*                                                                           */
/*  This data structure attempts to capture the attributes of an "analysis"  */
/*  of neuronal data.                                                        */
/*                                                                           */
/*  Analysis File:                                                           */
/*                                                                           */
/*  <type>                                                                   */
/*  <grouping> ...                                                           */
/*  <output_format>                                                          */
/*  <nparam>                                                                 */
/*  <name> <value>                                                           */
/*  ...                                                                      */
/*  condition                                                                */
/*  ...UNSPECIFIED...                                                        */
/*                                                                           */
/*  For example:                                                             */
/*                                                                           */
/*  psth                                                                     */
/*  group 0                                                                  */
/*  histogram                                                                */
/*  2                                                                        */
/*  binsize 10                                                               */
/*  sampling 1000.0                                                          */
/*                                                                           */
/*  IN GENERAL                                                               */
/*    <grouping> is one of the following:                                    */
/*       all                                                                 */
/*       individual                                                          */
/*       group <k> <name_1> ... <name_k>                                     */
/*    and "group 0" indicates the default grouping.                          */
/*                                                                           */
/*  STANDARD ANALYSIS PARAMETER NAMES (with example values):                 */
/*    sampling 1000.0                                                        */
/*    start 0                                                                */
/*    period 2000                                                            */
/*    chan                                                                   */
/*                                                                           */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/
struct nda_struct{ // For analysis
  char *atype;     // Type of analysis, e.g., "psth_simple", "isi", "power"
  char *grouping;  // Maybe, "all", "individual", "group" ...
  char *outtype;   // Type of plot to make, e.g., "simple", "histogram"
  int nparam;      // Number of parameters for analysis
  char **pname;    // Param names, e.g., "binsize", "chan1", "m", "sampling"
  char **pval;     // Param values
  struct ndindex_struct *ndx;  // Index on trials
  struct ndcond_struct *ndc;  // Pointer to condition structure
};
/*****************************************************************************/
/*                                                                           */
/*                                GENERATING                                 */
/*                                                                           */
/*  This data structure contains the specification for making spike trains.  */
/*                                                                           */
/*  Format:                                                                  */
/*                                                                           */
/*  <mtype>                                                                  */
/*  nfiles 3                                                                 */
/*  <name> [i,c,f] <value>                                                   */
/*  <name> [i,c,f] <value>                                                   */
/*  ...                                                                      */
/*  [nvar|nvar_link] <nvparam>                                               */
/*  <vname1> <vtype1> <vval1_1> ... <vval1_n>                                */
/*  ...                                                                      */
/*  <vnameX> <vtypeX> <vvalX_1> ... <vvalX_n>                                */
/*  population <npop>                                                        */
/*  <popname1> <distrib_name> <distr_val> ... <distrib_val>                  */
/*  ...                                                                      */
/*  <popnameX> <distrib_name> <distr_val> ... <distrib_val>                  */
/*                                                                           */
/*  For example:                                                             */
/*  -----------------------------------------------------------------------  */
/*  poisson                                                                  */
/*  nfiles 20                                                                */
/*  ntrials 100                                                              */
/*  period 2000                                                              */
/*  sampling 1000.0                                                          */
/*  lambda1 50.0                                                             */
/*  sigma 2.0                                                                */
/*  ...                                                                      */
/*                                                                           */
/*  nvar 2   [OR, use 'nvar_link 2' for linked params]                       */
/*  lambda2 f 10.0 30.0 50.0 70.0 90.0                                       */
/*  njumps  i 5 6 7 8 9                                                      */
/*                                                                           */
/*  npop 1                                                                   */
/*  seed ran2 1777                                                           */
/*  -----------------------------------------------------------------------  */
/*                                                                           */
/*  STANDARD PARAMETER NAMES (with example values):                          */
/*    sampling 1000.0                                                        */
/*    period 2000                                                            */
/*    ntrials 100                                                            */
/*    nrepeats 2  (used when nvar > 0, and ntrials is not specified)         */
/*                                                                           */
/*****************************************************************************/
struct ndgen_struct{ /*** For generating spikes. ***/
  char *gtype;     /* Name of algorithm */
  int nfiles;      /* No. of files to make (could be for cells or sites) */
  int nparam;      /* Number of constant parameters. */
  char *ptype;     /* Const param types ('c' by default) */
  char **pname;    /* Param names */
  char **pval;     /* Param values */
  int nvar_link;   /* Are variable params linked? 0, 1 */
  int nvar;        /* Number of variable ("stimulus") parameters. */
  char **vname;    /* Var. param names */
  char *vtype;     /* Var. param types: i,c,f */
  char **vvalstr;  /* Var. param values (in one string) */
  int npop;        /* Number of population variables */
  char **popname;  /* Population variable names */
  char **popdd;    /* Population distribution descriptions */
  char **popval;   /* Population variable values (reused, passed to routine) */
};
