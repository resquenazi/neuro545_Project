/*** pop_util.h ***/

// HID float pop_util_ran2();
// HID float pop_util_gasdev();
extern float *pop_util_get_mask_diff_exp();
extern float *pop_util_get_mask_diff_exp_mech();
extern float *popu_mech_onode_get_mask_alpha();

// Getting values from onodes
extern void pop_util_onode_get_pop_xn_yn_zn();

// Layer, low level
extern struct pop_area  *pop_util_get_area_by_name();    // NEW mpt
extern struct pop_layer *pop_util_get_layer_for_name(); // old
extern struct pop_layer *pop_util_get_layer_by_name();  // NEW mpt
extern               int popu_get_layer_index_by_name();
extern struct pop_layer *popu_get_lay_named();
extern struct pop_mech  *popu_get_mech_named();
extern struct pop_cell  *popu_get_c_named();

extern void pop_util_get_xn_yn_zn_by_name();
extern struct pop_map *pop_util_get_map_by_name();      // NEW mpt


extern void pop_util_map_make();  // Main map maker, calls mod_conn_...


//  For warping the maps (created for Lars' project)
extern float **popu_warp_apply();
extern float   popu_warp_remap();
extern void    popu_warp_map_coords();
extern float **popu_warp_get_map_arrays();
extern float **popu_warp_gen_params();
extern float **popu_warp_map_simple();
extern void    popu_warp_map();


// Queries related to "laytype"
extern int pop_util_check_lclass();
extern int pop_util_lclass_is_lgn0();
extern int pop_util_lclass_is_rgci();


// structure creation, initialization
extern              void popu_free_pop_cell();
extern              void popu_free_pop_layer();
extern              void popu_free_pop_area();
extern              void popu_free_pop_map();
extern              void popu_free_pop_top();

extern struct pop_cell ***popu_make_init_cells();
extern struct pop_layer *popu_make_pop_layer();
extern struct pop_layer *popu_make_layer_irreg();
extern struct pop_layer *popu_make_layer_cart();
extern struct pop_layer *popu_make_layer_template();

extern struct pop_top   *popu_make_pop_top();

extern struct pop_area  *popu_make_area_default();
extern struct pop_area  *popu_make_area();


// Geometry
extern              void popu_init_geom();


// Inputs
extern struct pop_layer *popu_input_get_origin();
extern int popu_input_is_paired();


// NMDA related
extern float *pop_util_get_nmda_shape();
extern void pop_util_compute_nmda_table();
extern void pop_util_nmda_match_intcurr_waveform();
extern void pop_util_compute_intcurr_ampa_nmda_layer();
extern void pop_util_compute_intcurr_ampa_nmda_lgn();
extern void pop_util_nmda_prep_constants_layer();


extern int popu_flag_resp();


// Writing and returning data from one trial run
extern void pop_util_response_link();
extern void pop_util_response_handover();

// HID void pop_util_cell_spike_inc();
// HID void pop_util_add_spike_to_cb();

//    float pop_util_circbuf_val_at_t();
//    float pop_util_circbuf_val_at_t_interp();
//     void pop_util_circbuf_mask_update_t();

// HID void pop_util_process_sid_ds02();
// HID void pop_util_process_sid_symmask();
// HID void pop_util_process_si001();
// HID void pop_util_process_si_spike();
// HID void pop_util_process_spike_postsyn();


// Mechanisms
extern void pop_util_mech_config_mask();


// Input Pair
extern struct input_pair *pop_util_input_pair();


// Cell styles
extern void pop_util_make_cell_style_ds01();
extern void pop_util_make_cell_style_ds02();



// General customization
extern void pop_util_customize();


// Accumulating inputs for synapses
extern void   pop_util_add_noise();
extern void   pop_util_set_flag_grid_input_layer();
extern float *pop_util_expand_spikes();
extern float *pop_util_expand_grid_spikes();
extern void   pop_util_add_gt_sum_grid_spikes();
extern void   pop_util_gt_sum_grid_pair();
extern void   pop_util_set_gad_layer();
extern void   pop_util_set_gt_on_off_input_layer();

// Inducing correlation between spike trains
extern void pop_util_corr_nearby_spikes();


// Background firing
extern void pop_util_init_bg_spikes();
extern void pop_util_init_bg_rate();


// Run Simulation

// HID void pop_input_01();
// HID void pop_input_02();
// HID void pop_util_deriv_01();
// HID void pop_util_deriv_02();
// HID void pop_rkck_01();
// HID void pop_rkqs_01();
// HID void pop_process_spike_01();
// HID void pop_step_01();
// HID void pop_util_init_cell_trial();

extern void pop_util_run_layer_list();
extern void pop_util_run_cell_list();


// Writing .mar files  [Note, perhaps this belongs in 'mod_util' ??]
extern void popu_mar_write_onode_items();
extern void popu_mar_write_area_name();
extern void popu_mar_write_area();
extern void popu_mar_write_lay();
extern void popu_mar_write_syn_reg();
extern void popu_mar_write_syn_lgn();
extern void popu_mar_write_input();
extern void popu_mar_write();

//
// Response file generation
//
// struct pop_layer *popu_rsp_gen_get_lay_ptr();
//            int ***popu_rsp_gen_cmap_from_w();
//              void popu_rsp_gen_append_by_map();
//               int popu_rsp_gen_expand()
extern void popu_rsp_gen();
