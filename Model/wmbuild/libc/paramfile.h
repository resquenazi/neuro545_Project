/*****************************************************************************/
/*                                                                           */
/*  paramfile.h                                                              */
/*  Wyeth Bair                                                               */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/

struct param_pair{
  char *name;
  char *value;
  struct param_pair *next;
};

struct param_pair_list{
  int n;
  struct param_pair *p;    // This list checked for unique 'name'


  // "CUSTOMIZE"  OLD - ONLY USED BY mod_spop_util, and pop_util ???
  int nc;
  struct param_pair *c;    // This list not checked for unique 'name'

  int ns;     // Number of string data lines
  char **s;   // [ns] String data
};


struct onode {

  //  ***
  //  ***
  //  ***  ANY CHANGES here must be reflected in 'paramfile_util.c' routines
  //  ***  that convert this structure to/from a carray
  //  ***
  //  ***

  char *otype;      // Node type, e.g., 'pop', 'geometry'
                    //   Reserved values:  'item'

  int resolve;      // 0-nothing to resolve
                    // 1-something to resolve, e.g., 'COPY'


  int tv_n;         // Number of tag values
  char **tv_name;   // Names of tag values [tv_n]
  char **tv;        // Tag values [tv_n]

  char *otag;       // Label for this node - WYETH - currently unused.
  int n;            // Number of children
  struct onode *o;  // head node in child list

  struct onode *prev;  // prev in this list
  struct onode *next;  // next in this list

  char *name;       // Item name
  char *val;        //      value, e.g., "2.0" or "1.0 0.5 0.0 0.5 1.0"
  char *unit;       //      units
  char *comm;       //      comment

  int nval;         // Number of values, added 20 Jan 09 to allow multiple vals
                    // Note, even if nval > 1, a single string is stored in 
                    // 'val', as shown above.


  //  ***
  //  ***
  //  ***  ANY CHANGES here must be reflected in 'paramfile_util.c' routines
  //  ***  that convert this structure to/from a carray
  //  ***
  //  ***
};
