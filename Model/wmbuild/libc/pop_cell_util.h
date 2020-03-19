

// General utilities
extern void pop_cell_parse_unit_name();
extern int pop_cell_get_layer_index();
extern struct pop_layer *pop_cell_get_layer_pointer();
extern struct pop_cell *pop_cell_laylist_get_cell_pointer();

// Circular buffer
extern struct pop_circbuf *pop_cell_get_init_circbuf();
extern void pop_cell_zero_circbuf();
extern void pop_cell_zero_circbuf_tm1();

// SI utilities, making and resetting
extern void   pop_cell_si_free();
extern void   pop_cell_si_sav_clear();
extern void   pop_cell_si_make_001();
extern void   pop_cell_si_reset_001();
extern void   pop_cell_si_make_002();
extern void   pop_cell_si_reset_002();
extern void   pop_cell_si_reset();

extern void   pop_cell_si_sav_set_n();
extern float *pop_cell_si_get_data_ptr();
extern void   pop_cell_si_add_sav();

// csav
extern   int popc_csav_get_csi();
extern  void popc_csav_set_f();


extern void pop_cell_syn_get_ptr_data_float();
extern void pop_cell_get_ptr_data_float();
extern void pop_cell_get_data_float();
extern void pop_cell_syn_get_data_float();

extern void pop_cell_write_sav_layer_list();


// Cell utilities
extern void popc_precomp_reset();
extern float *popc_precomp_get_expanded_spikes();
extern void pop_cell_init_cell();
extern void pop_cell_print_cell();

// attribs
extern void  pop_cell_attrib_init_f();
extern int   pop_cell_attrib_index_f();
extern void  pop_cell_attrib_set();
extern float pop_cell_attrib_get_f();
extern void  pop_cell_attrib_print_f();
extern void  pop_cell_attrib_default_set();
extern int   popc_attrib_vary_xy();
extern int   popc_attrib_variation_code();
extern int  *popc_attrib_layer_variation();
extern void  popc_attrib_dump_layer();




// Querying inputs
extern int popc_lay_count_input_from_lay();


// Synapse Utilities - printing, querying
extern void pop_cell_syn_free();
extern void pop_cell_print_synint();
extern void pop_cell_print_pre_synaptic_list();
extern void pop_cell_print_post_synaptic_list();
extern struct pop_syn *pop_cell_get_next_syn_cell();
extern struct pop_syn *pop_cell_get_next_syn_layer_name();
extern struct pop_syn *pop_cell_get_next_syn_index();

extern struct pop_syn *pop_cell_laylist_get_syn_pointer();

extern int   pop_cell_count_syn_from_layer_to_cell();
extern float pop_cell_max_weight_from_layer_to_cell();
extern float pop_cell_tot_weight_from_layer_to_cell();
extern float pop_cell_tot_weight_inindex_to_cell();

extern int   pop_cell_count_syn_from_cell_to_layer();
extern float pop_cell_max_weight_from_cell_to_layer();
extern float pop_cell_tot_weight_from_cell_to_layer();

// 3D maps of synapses
extern int      popc_wf_accum_lgn_to_unit();
extern int      popc_wf_accum_lay_to_unit();
extern float ***popc_wf_accum_lay_to_lay();
extern float ***popc_wf_accum_multi_lay();

extern int popc_flag_accum_lay_to_unit();
extern int popc_flag_accum_lay_to_lay();
extern int popc_flag_get_lay_to_all();

extern int popc_flag_accum_lay_ori();


// LGN connections
extern struct pop_lgn_in *pop_cell_get_lgn_input_by_name();
extern struct pop_lgn_in *pop_cell_get_lgn_input_by_lay();
extern struct pop_lgn_in *pop_cell_get_lgn_input_pair_by_lay();
extern void pop_cell_add_lgn_input();
extern void pop_cell_add_lgn_pair();

// Mechanism queries
extern struct pop_mech *pop_cell_get_mech_by_name();

// Mechanism creation
extern void pop_cell_layer_add_mech();




// Cellist utilities
extern struct pop_cell **pop_cell_cellist_get_presyn_from_layer();


// More synapse utilities

extern void pop_cell_reset_all_out_synapses();

extern struct pop_syn *pop_cell_syn_si_find();  // Find paired synapse

extern struct pop_syn *pop_cell_add_synapse();
extern void pop_cell_remove_synapse();
extern void pop_cell_remove_all_syn_from_cell();
extern void pop_cell_remove_all_syn_from_layer_to_cell();
extern void pop_cell_remove_all_syn_from_cell_to_layer();
extern void pop_cell_remove_all_syn_from_cell_to_layer_si();

extern void pop_cell_adjust_all_syn_from_cell_to_layer();

// Modifying cell properties
extern void pop_cell_save_all_syn_from_layer_to_cell();
extern void pop_cell_copy_syn_layer_from_cell_to_cell();
extern float pop_cell_norm_weight_from_layer();
extern float pop_cell_norm_weight_from_inindex();
extern void pop_cell_norm_weight_layer_layer();
extern void pop_cell_norm_weight_layer_inindex();

// Irregular syn connections
extern void popc_layirr_get_conditional_coords();
extern int **popc_layirr_get_conditional_tflag();
extern void pop_cell_add_on_off_input_irr();

// Setting save parameters
extern void pop_cell_syn_set_save();
extern void pop_cell_set_save();
extern void pop_cell_unset_save();
extern void pop_cell_unset_save_syn();

// Derive attributes for cells
extern void pop_cell_derive_ori_from_layer();

