/*extern void ndata_read_nchar();*/
/*extern void ndata_write_nchar();*/
extern void init_const_param();
extern void init_set_const_param();
extern void init_set_var_param();
extern void ndata_add_uninitialized_var_param();
extern FILE *ndata_open_set_revflag();
extern struct ndata_struct *ndata_read_head_after_revflag();
extern struct ndtrial_struct *ndata_read_get_next_trial();
extern long ndata_get_file_pos_ntr();
extern void ndata_overwrite_ntr();
extern void read_ndata();
extern void read_ndata_text();

/*extern void ndata_write_trial();*/
extern void write_ndata();
extern void ndata_append_trial();
extern void write_ndata_link();

extern void write_ndata_text();
extern void write_ndata_txt();
extern void write_ndata_t1();
extern void ndata_util_read_t1();
extern void ndata_util_read_t2();
extern void read_psych_text_ndata();
extern void ndata_free_ect();
extern struct ndata_struct *get_ndata();
extern void free_ndata_rec();
extern void free_ndata_trial();
extern void free_ndata();
extern void free_ndata_link();


// Getting information from a trial
extern void   ndata_trial_print_trial();
extern int    ndata_trial_get_tcode();
extern int    ndata_trial_get_var_param_index();
extern int    ndata_trial_get_record_index();
extern struct ndrec_struct *ndata_trial_get_record_ptr();
extern int    ndata_trial_get_record_type();
extern float  ndata_trial_get_sampling_chan();
extern int    ndata_trial_get_rcode_chan();
extern int    ndata_trial_get_start_chan();
extern float  ndata_trial_get_start_sec_chan();
extern float  ndata_trial_get_dur_sec_chan();
extern int    ndata_record_get_nth_occur_evco_time();
extern float  ndata_trial_get_start_sec_schan();
extern float *ndata_trial_getp_adata_chan();
extern int   *ndata_trial_getp_point_chan();
extern int   *ndata_trial_getp_evco_list_chan();
extern int    ndata_trial_chan_evco_test();
extern int    ndata_trial_get_n_points_chan();
extern int    ndata_trial_get_nth_point_chan();
extern int    ndata_trial_get_event_time_chan();

// Creating, Coping, Merging records and trials
extern void ndata_free_record_contents();
extern struct ndrec_struct *ndata_record_create_type_1();  // float
extern void ndata_record_copy_contents();
extern struct ndtrial_struct *ndata_trial_create();
extern void ndata_trial_insert_record();
extern void ndata_trial_delete_record_flag();
extern void ndata_trial_copy_contents();
extern struct ndata_struct *ndata_copy();
extern void ndata_insert_event_records();  // not general

//  Timetrial (channels are broken across trials into smaller time segments)
extern int *ndata_timetrial_get_merge_i();


//  Event ecode table
extern void ndata_lookup_ecode_name();
extern int  ndata_lookup_ecode_num();
extern int  ndata_ecode_get_int();      // Integer ecode from either string/int

//  Event queries
extern int ndu_trial_chan_ecode_str_diff();


extern void ndata_get_default_chan();
extern void ndata_get_first_unit_chan();
extern void ndata_get_first_prefix_chan();
extern void ndata_get_first_event_chan();
extern void ndata_get_record_index();
extern void ndata_free_ndrx();
extern void ndata_get_adata_record_pointer();
extern float ndata_get_sampling_record_index();
extern int ndata_get_type_record_index();
extern float ndata_get_sampling_record_name();
extern int ndata_get_type_record_name(); // Get a `channel' type

// Parameter query functions
extern int ndata_get_const_type();
extern int ndata_get_var_type();
extern int ndata_get_var_or_const_type();
extern int ndata_convert_var_type_int_float();
extern void ndata_convert_var_numeric();

// nda script params
extern int    ndata_test_nda_param();
extern int    ndata_get_nda_param();
extern int    ndata_set_nda_param();   // To overwrite values in nda struct

extern char  *ndata_get_nda_param_char_or_exit();   // Long names (original)
extern int    ndata_get_nda_param_int_or_exit();
extern float  ndata_get_nda_param_float_or_exit();

extern char  *nda_getc_exit();      // Shorter names (new)
extern int    nda_geti_exit();      // Shorter names
extern float  nda_getf_exit();      // Shorter names

extern void   ndata_set_nda_param_default_int();
extern void   ndata_set_nda_param_default_float();
extern void   ndata_set_nda_param_default_char();

extern int    nda_geti_dflt();      // Shorter names (new)
extern float  nda_getf_dflt();      // Shorter names
extern char  *nda_getc_dflt();      // Shorter names



extern int    ndata_get_const_param();
extern int    ndata_get_const_param_index();
extern int    ndata_get_const_param_int();
extern float  ndata_get_const_param_float();
extern char  *ndata_get_const_param_char();
extern int    ndata_get_var_param_index();
extern int    ndata_get_var_value_trial();
extern int    ndata_get_var_value_trial_int();
extern int    ndata_get_var_value_trial_float();
extern int    ndata_get_var_or_const_param();
extern int    ndata_get_vc_param_int();
extern int    ndata_get_vc_param_int_or_exit();
extern int    ndata_get_vc_int_exit();
extern int    ndata_get_vc_int_dflt();
extern int    ndata_get_vc_param_float();
extern float  ndata_get_vc_param_float_or_exit();
extern float  ndata_get_vc_flt_exit();
extern float  ndata_get_vc_flt_dflt();

extern int ndata_match_trials_var_params();

extern float ndata_get_tf_hz_pep_cycle_trial();
extern float ndata_get_duration_pep_cycle_trial();

extern void ndata_add_const_param();
// See also:  "ndata_add_uninitialized_var_param"
extern void ndata_add_var_param();


/*** Var parameters ***/
extern void ndata_get_unique_var_values();
extern void ndata_get_var_param_frequency();
extern void ndata_get_var_param_info();

/*** Parameter set functions. ***/
extern void ndata_set_var_value_trial();


extern void ndata_print_ect();
extern void print_ndata();
extern void print_partial_ndata();

// Event code tables
extern struct ect_struct *ndu_get_check_ect_chan();
extern int ndata_diff_ect();
extern struct ect_struct *ndata_merge_ect();
extern struct ect_struct *ndata_copy_ect();

// Condition utilities
extern struct ndcond_struct *get_empty_condition();
extern void ndcond_free_condition();
extern void print_condition_list();
extern void test_condition_list();
extern void ndcond_remove();
extern void ndcond_add();
extern int ndcond_check_point_count();
extern int ndcond_check_evco_sequence();
extern int ndcond_check_evco_value();
extern int ndcond_check_var_param();
extern int ndcond_check_vc_param();
extern int ndcond_check_special_decision();
extern int check_trial_ndata_condition();
extern int check_trial_ndata_condition_one();
extern void ndcond_point_count_max_dev();
extern void ndcond_point_count_max_dev_two_sided();
extern int *ndata_condition_get_flags();
extern int count_trials_ndata_condition();
extern int count_trials_ndata_condition_chan();
extern void get_sarray_ndata_condition();
extern void get_double_sarray_ndata_condition();

extern int get_max_period_record();
extern int ndata_group_is_param_const();
extern int ndata_group_get_par_int();
extern int ndata_count_record_in_group();
extern int ndata_search_groups_for_trial();
extern void get_farray_ndata_group();
extern void get_sarray_ndata_group();
extern void get_sarray_em_raster_ndata_group();
extern void ndata_get_adata_group();

// Ndata index utilities
extern void print_ndata_index();
extern void make_ndata_index();
extern void free_ndata_index();
extern int count_record_in_ndata_index();
extern void get_sarray_ndata_index();
extern void ndata_get_adata_pointer_ndata_index();
extern void ndata_get_adata_adjusted_ndata_index();

// Start time adjustment
extern void ndata_record_adjust_start_time();
extern int ndata_trial_get_adjust_time_sec();
extern void ndata_trial_adjust_start_with_point();
extern void ndata_trial_adjust_start_with_event();
extern void ndata_adjust_start_time();
extern void ndata_adjust_start_time_ndata_index();
extern void ndata_adjust_start_time_global();
extern void ndata_get_start_time_ndata_index();
extern void ndata_var_start();

// Variable duration
extern int  ndu_vardur_trial();
extern int *ndu_vardur_group();


// Group and Index utilities
extern void make_ndata_analysis();
extern void ndu_lookup_resolve();
extern void read_ndata_analysis();
extern void set_ndata_analysis_index();
extern void make_name_for_group();
extern void make_space_name_for_group();
extern char *ndata_make_name_var_trial();
extern char *ndata_make_name_index_trial();
extern char *ndata_make_name_index_trial_space();

extern void print_group_struct();
extern void ndata_free_group();
extern void get_group_all();
extern void get_group_all_unconditional();
extern void get_group_individual();
extern void ndata_get_groupdef_group();
extern void get_grouped_data();
extern void get_primary_list_from_grouped_data();
extern void prepare_nda_arg();
extern void prepare_nda();

extern void ndata_get_max_period_sampling_chan();
extern void ndata_establish_basic_params();
extern void ndata_establish_start_period_sampling();

// Making ndata from sarrays
extern struct ndata_struct *make_ndata_sarray();
extern struct ndata_struct *make_ndata_sarray_pair();

// Analysis-specific routines
extern float *ndata_get_normd_resp_chan_vdur();
extern float *ndata_get_normalized_response_chan();
extern float *ndata_get_adjusted_spike_count();
extern void ndata_get_adjusted_mean_sdev_group();

// Generating Spike Trains
extern void read_ndata_gen();
extern int ndgen_get_param_const_pop();
extern void ndgen_get_param_const_pop_def_int();
extern void ndgen_get_param_const_pop_def_float();
extern void ndgen_get_param_const_pop_def_char();
extern int ndgen_get_param_trial();
extern int ndgen_get_param_trial_int();
extern float ndgen_get_param_trial_float();
extern char *ndgen_get_param_trial_char();

// Group Lists and Group Nodes
extern struct ndglist *ndata_glist_get_init();
extern void ndata_glist_free_pointers();
extern void ndata_glist_print();
extern struct ndgnode *ndata_glist_find_group_for_trial();
extern struct ndgnode *ndata_glist_add_gnode();
extern void ndata_gnode_add_trial();
extern struct ndgnode *ndata_glist_add_trial();
extern int ndata_gnode_count_trials_with_chan();
extern void ndata_gnode_get_sarray();

// Custom data processing, called from 'nda.c'
extern void ndata_get_mean_sdev_2d_fdata();
extern void ndata_get_stat_2d_fdata();

// Choice prob analysis
extern int ndata_util_choice_gen();
extern void ndata_util_choice_01();
