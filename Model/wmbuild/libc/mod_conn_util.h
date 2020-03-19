// zstring conditions, z-flag
extern int  mod_conn_zstring_get_unique_z();
extern int *mod_conn_zstring_get_zflag();


// Comparing ori differences
extern float mod_conn_get_ori_diff_cells();
extern float mod_conn_get_ori_diff();
extern float mod_conn_get_ori_diff_one();

// For color map
extern    void mod_conn_color_angle();
extern    void mod_conn_color_sf();
extern    void mod_conn_test_color_map();
extern    void mod_conn_write_map_text();
extern    void mod_conn_write_ori_map();
extern float **mod_conn_get_ori_map_old();
extern float **mod_conn_get_ori_map();
extern    void mod_conn_onode_make_map_ori();
extern    void mod_conn_set_layer_ori_map();

extern void mod_conn_onode_make_map_sf();      // SF
extern void mod_conn_set_layer_sf_map();
extern void mod_conn_onode_make_map_od();      // OD
extern void mod_conn_set_layer_od_map();

extern void mod_conn_onode_make_map_rfx();     // RFX
extern void mod_conn_set_layer_rfx_map();

extern void mod_conn_set_layer_dir_map_flag(); // Dir


extern void mod_conn_pick_paired_coord();


// Reading/Writing connections to/from file
//extern void mod_conn_write_conn_txt();
extern void mod_conn_write_conn_t1();
extern void mod_conn_read_conn_t1();
extern void mod_conn_write_conn_t0();  // LGN
extern void mod_conn_read_conn_t0();  // LGN
extern void mod_conn_read_conn();
extern void mod_conn_conn_to_text();

// Querying an LGN connection
extern int mod_conn_is_lgn_connected();
// Anatomy Stats
extern void mod_conn_astat_lgn_postsyn();
extern void mod_conn_astat_lms_input();

// Picking paired connections
extern void mod_conn_pick_paired_connections();
extern int mod_conn_convert_syn_to_phase_pair();
extern void mod_conn_syn_convert_cell();
extern void mod_conn_syn_convert_layer();

extern void mod_conn_input_pair_main();

// Picking connections from probability matrix
extern void mod_conn_meth_01();
extern void mod_conn_meth_02();

// Connecting a layer to a grid according to some rule
extern void mod_conn_gabor_01();
extern int mod_conn_gabor_02();
extern int mod_conn_gabor_03();
extern int mod_conn_gauss_02();
extern void mod_conn_layer_gabor_map();
extern void mod_conn_onode_layer_gabor_map_set_phase();
extern void mod_conn_onode_layer_gabor_map();
extern void mod_conn_lgn_pair();

// clist - connection lists
extern int mod_conn_clist_count_shared();

//
// Inter-layer connections
//
extern void mod_conn_input_plain();

extern void mod_conn_corr_on_off();
//  float **mod_conn_composite_rf_irr();
//  float **mod_conn_composite_rf();
extern void mod_conn_corr_on_off_01();
extern void mod_conn_onode_corr_on_off();

extern void mod_conn_onode_binoc_mask();
extern void mod_conn_onode_gabor_mask();

extern void mod_conn_layer_to_cell_ori_dist_01();
extern void mod_conn_ori_dist_01();
extern void mod_conn_onode_ori_dist_01();

extern void mod_conn_cell_to_layer_gauss();
extern void mod_conn_layer_to_cell_gauss();
extern void mod_conn_layers_gauss();

// Define specific SI (synaptic interaction) types
//REMOVED extern struct pop_mech *mod_conn_def_si_ds02();

// Writing connection patterns for graphical output
extern float **mod_conn_image_2d_lgn_grid();

// Derive cell attributes


