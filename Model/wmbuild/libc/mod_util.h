// Cluster computing, logging to files
extern char *mylog_get_global_fname();
extern char *mylog_get_global_ja_fname();
extern int   mylog_set_id_log();

extern int mod_util_sampling_equal();

extern void mod_util_substitute_var_params();

extern int mod_util_set_var_tn();

// Variable parameters
extern                  void  mod_util_var_par_rec_free();
extern                  void  mod_util_var_par_lookup_read();
extern                  void  mod_util_var_par_lookup_replace();
extern                  char *mod_util_var_par_lookup_reverse();
extern    struct var_par_rec *mod_util_var_par_rec_get_one();
extern    struct var_par_rec *mod_util_var_par_rec_from_onode();
extern                   int  mod_util_var_par_rec_check_moo();
extern                  void  mod_util_var_par_free();
extern struct var_par_struct *mod_util_var_par_get_init();
//                   void     mod_util_var_par_cross_table();
extern struct var_par_struct *mod_util_var_par();


// MDS - model data struct
extern void modu_mds_write();
extern void modu_mds_create();
extern void modu_mds_insert();
extern void modu_mds_add_3d();


// Putting results into the response structure
extern int    mod_util_resp_unfilled_count();
extern void   mod_util_resp_unfilled_exit();
extern char  *mod_util_resp_get_rtype_ptr();
extern int    mod_util_resp_get_ri();
extern float  mod_util_resp_get_samp();
extern void   mod_util_resp_set_samp();
extern int    mod_util_resp_check_name();
extern int    mod_util_resp_get_check_name();
extern int    mod_util_resp_check_store_f();
extern int    mod_util_resp_check_store_s();
extern float *mod_util_resp_derive_current();
extern void   mod_util_resp_store_s();           // NEW WAY
extern void   mod_util_resp_store_f();
extern void   mod_util_resp_store_f_samp();
extern void   mod_util_resp_store_f_reg();

// Response:  save_pop_unit_as
extern int    mod_util_resp_pop_unit();


// retina0 cone mosaic utilities
extern float mod_util_cone_weight();


//**** WYETH useful stuff was moved to KERNEL_UTIL ****/
//  2D filters
//
//    float **mod_util_o_get_filt2d_gabor();   // OLD STYLE
//
extern float **mod_util_o_get_filter_2d();    // OLD

//
//  Stimulus generation
//
extern void mod_util_stim_smooth();
extern void mod_util_get_3d_stim();
extern void mod_util_get_3d_stim_binoc();
extern void mod_util_get_3c_stim();
extern void mod_util_get_3rgb_stim();
extern float *mod_util_get_1d_stim();
