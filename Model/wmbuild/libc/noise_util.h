extern void get_mseq_tap_array();
extern void noise_util_get_seed_array();

extern float noise_util_ran2();
extern float noise_util_gasdev();

extern void get_stimulus_word_index_format();

extern void convert_jump_to_position();
extern void convert_position_to_toward_away();

extern int get_stim_screen_y_const();
extern void get_wn_stim_params();
extern void get_wy2pat_stim_params();
extern void get_wy2pat_stim_params_02();
extern void get_wy2pat_stim_params_03();
extern void get_wy2pat_rel_stim_params();
extern void get_wy2pat_simple_stim_params();
extern void get_wy2flat_rel_stim_params();
extern void get_wy2pa_stim_params();
extern void get_wytfv_stim_params();
extern void get_addot_stim_params();
extern float get_wn_actual_dt();
extern float get_sgi_actual_dt();
extern void get_actual_screen_dimensions();
extern int get_actual_latency();
extern int get_sync_code_int();

extern void get_wy2pat_stim_for_params();
extern void get_wn_stim_for_params();
extern void get_amp_mod_stim_for_params();
extern void get_4state_rel_stim_for_params();
extern void get_4state_simple_stim_for_params();
extern void get_wy2pa_antidur_list();
extern void get_wy2pa_stim_for_params();
extern void get_wytfv_stim_for_params();

extern void get_wn_stim_for_trial();
extern void get_wy2pat_stim_for_trial_01();
extern void get_wy2pat_stim_for_trial_02();
extern void get_wy2pat_stim_for_trial_03();
extern void get_wn_md_stim_for_trial();
extern void get_wy2pa_stim_for_trial();
extern void get_wytfv_stim_for_trial();

extern void get_stim_at_resolution_simple();
extern void get_stim_at_resolution_piecewise();
extern void get_stim_at_resolution_impulse();
extern void get_stim_at_resolution();
extern void get_new_noise_simple_stim_group();
extern void get_noise_simple_stim_group();
extern void get_noise_simple_md_stim_group();
extern void get_noise_simple_stim_vel_group();
/* INTERNAL: get_special_wy4st_group */
extern void get_two_noise_stim_group();
/* INTERNAL: get_new_noise_stim_group(); */
extern void get_noise_stim_group();
extern void get_model_grid_stim_group();
extern void get_labr_grid_stim_group();
extern void get_ndata_stim_group();
extern void get_ndata_expanded_stim_group();
extern void get_std_stim_group();
extern void get_std_stim_group_expanded();


extern void get_4state_stim_for_trial();
extern void get_4state_stim_lim_trial();
extern void get_four_state_stim_group();
extern void get_trig_for_4state_stim();
extern void get_4state_rel_trig_trial();
extern void get_trig_4state_rel_group();
extern void get_4state_simple_trig_trial();
extern void get_trig_4state_simple_group();

extern void get_wy2pa_stim_group();
extern void get_wytfv_stim_group();
extern void write_wytfv_stim();

extern void get_addot_stim_seq_group();

extern int get_seed40_index();

extern void get_map_stim_params();
extern void get_map_stim_for_trial();
extern void get_map_stim_group();
extern void get_map_psth();

// star
extern void noise_util_star_get_trig();

