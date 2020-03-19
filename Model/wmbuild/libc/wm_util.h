// Stimulus Structure
extern void wm_prep_var_const_gen();
extern void wm_response_count();
extern void wm_response_config();
extern void wm_response_check_trial_full();

extern void wm_get_sform_xn_yn_tn();
extern void wm_response_init();
extern void wm_get_nd_t1_outfile();
extern void wm_make_ndtrial();
extern struct ndata_struct *wm_make_ndata();
extern void wm_write_response();

extern  int wmu_calibration_init();
extern void wmu_calibration_write();

extern void wmu_wmpi_stim();
extern void get_stimulus_data_tuning_curve();

// Model param variation
extern void wm_moo_var_par_set(); // Update the params of the model

extern void wm_slave_action();
