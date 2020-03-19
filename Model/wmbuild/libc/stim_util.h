
//     void get_cross_product_params();
//      int stim_util_get_vargen_list();
//      int stim_util_get_vargenpair_list();
//      int stimu_varseq_get_list();
//      int stimu_varseq_get_n();

extern void free_var_params();
extern void prep_var_params();
extern void stimu_var_param_check();
extern void prep_const_params();
extern void update_const_param();

extern  int stim_util_count_new_frames();

//     void get_coh_dot_movie_ampl();
//     void get_coherence_dot_movie();
//     void get_pairable_dot_movie();
//     void mask_dot_movie();
//     void translate_dot_movie();
extern void get_dot_stereogram();

extern float get_vsk_sum_of_sinusoids();

extern void stimu_rangrid();                // Coords w/i a random grid


extern float *get_stim_1d_ppl_pulse();
extern float *get_stim_1d_ppl_sine();
extern float *get_stim_1d_ppl_vsksos();
extern float *get_stim_1d_ppl_noise();
extern float *get_stim_1d_ppl_stat1();


extern void make_frame_flat();
extern void make_frame_scale_frame_zmid();
extern void make_frame_from_3d();
extern void frame_apply_aperture();
extern void apply_aperture_scale();
extern void multiply_frame();
extern void multiply_frame_gaussian();
extern void multiply_frame_gaussian_0();
extern void multiply_frame_sinusoid();
extern void make_frame_circ_cs_sine();
extern void make_frame_disk_flat();
extern void make_frame_circ_cs_flat();
extern void make_frame_add_triangle();
extern void make_frame_edge_01();
extern void make_frame_sinusoid();
extern void make_frame_sinusoid_0();  // Like prev, but not for tensor
extern void make_frame_sinusoid_scale();
extern void make_frame_bar_add_0();
extern void make_frame_bar();
extern void make_frame_bar_0();       // Like prev, but not for tensor
extern void make_frame_bars();
extern float **make_frame_bars_tmplt();
extern float **make_frame_barpatmo();
extern void make_frame_add_gpatch_sine();
extern void make_frame_add_gpatch_gabor();
extern void make_frame_add_2d_gauss();
extern void make_frame_add_rect_pix();
extern void make_frame_add_rect_pix_0();
extern void make_frame_add_sine_bump();
extern void make_frame_add_sine_bump_0();
extern void make_frame_binarize();
extern void make_frame_binarize_0();
extern void make_frame_add_dots();
extern void make_frame_add_dots_0();  // Like prev, but not for tensor
extern void make_frame_gp_trans();
extern void make_frame_gp_trans_0();
extern void make_frame_binary_noise();
extern void stimu_frame_disk_array();  // Write into an existing 2D array

extern void stim_util_limit_3d();
extern void stim_util_limit_2d();

// extern void stimu_std_make_draw();

//
//  New way to write stimuli
//
extern void stimx_bar();
extern void stimx_gabor_britten();
extern void stimx_3barlink();
extern void stimx_barnoise();
extern void stimx_barpatmo();
extern void stimx_barray();
extern void stimx_bar2ran();
extern void stimx_dirmod();
extern void stimx_plaid();
extern void stimx_gabor_grid();
extern void stimx_globmo_patch();
extern void stimx_dots();
extern void stimx_dot_paired();
extern void stimx_dots_stereogram();
extern void stimx_gp();
extern void stimx_image_set();
extern void stimx_noise();
extern void stimx_pasupathy_shape();
extern void stimx_rgb_sine();


// 3C - for 3-D color stimuli.  See also 'data_util' t11 format utilities
extern int ***stim_util_3c_from_tern_3d();
extern int ***stim_util_3c_2col_from_3d();
extern int ***stim_util_3c_4col_from_3d();


// Operations on 3D stimulus
extern float ***stim_util_blur();
extern void stim_util_add_modspot();

extern float ***get_stim_3d_noise();
extern void offset_block_stim_3d_noise();

extern float ***get_stim_quad_sine();
extern float ***get_stim_3d_sine();
extern float ***get_stim_sine_int_noise();
extern float ***get_stim_impulse();

extern float ***get_stim_3d_ppl_dynedge();
extern float ***get_stim_3d_ppl_disk();
extern float ***get_stim_3d_ppl_sine();
extern void     get_stim_3db_ppl_sine();

extern void get_stim_3db_ppl_randot_stereogram();
//     void     stim_util_check_set_nseq();
extern float ***get_stim_3d_ppl_sine_seq();
extern float ***get_stim_3d_ppl_sos();
extern float ***get_stim_3d_ppl_cpsine();
extern float ***get_stim_3d_ppl_plaid();
extern     void get_stim_3db_ppl_plaid();
extern     void get_stim_3db_ppl_gabor_grid();
extern float ***get_stim_3d_ppl_dmmask();
extern float ***get_stim_3d_ppl_dirmod();
extern     void get_stim_3db_ppl_dirmod();
extern float ***get_stim_3d_ppl_ransinphase();
extern float ***get_stim_3d_ppl_adamk();
extern float ***get_stim_3d_ppl_ran_sin();
extern float ***get_stim_3d_ppl_bar();
extern float ***get_stim_3d_ppl_3barlink();
extern float ***get_stim_3d_ppl_pixel();
extern     void get_stim_3db_ppl_bar();
extern float ***get_stim_3d_ppl_bar2ran();
extern float ***get_stim_3d_ppl_barray();
extern float ***get_stim_3d_ppl_barray_old();  // WYETH - OLD - REMOVE
extern float ***get_stim_3d_ppl_wy4st();
extern float ***get_stim_3d_ppl_x4st();
extern float ***get_stim_3d_ppl_wy2pa();
extern float ***get_stim_3d_ppl_gp();
extern float ***get_stim_3d_ppl_noise();
extern float ***get_stim_3d_ppl_image_set();
extern float ***get_stim_3d_ppl_barnoise();
extern     void get_stim_3db_ppl_barnoise();
extern     void get_stim_3db_ppl_bar1ran();
extern float ***get_stim_3d_ppl_barpatmo();
extern float ***get_stim_3d_ppl_dots();
extern float ***get_stim_3d_ppl_dot_paired();
extern float ***get_stim_3d_ppl_pasupathy_shape();
extern float ***get_stim_3d_ppl_frameset();
extern float ***get_stim_3d_ppl_upload();

extern   int ***get_stim_3c_ppl_sine();
extern   int ***get_stim_3c_ppl_noise();
extern   int ***get_stim_3c_ppl_surrt();
extern   int ***get_stim_3c_ppl_rgb_sine();


extern float  ***get_stim_3d_ppl();    // Monocular stimuli
extern   int  ***get_stim_3c_ppl();    // Monocular color, packec INT
extern float ****get_stim_3rgb_ppl();  // Moncoular full RGB color
extern  void     get_stim_3d_b_ppl();  // Binocular stimuli

// General 2D routines
extern float **stim_get_2d_impulse();

// 2D xt stimuli - made from 'param_pair_list'
extern float **stim_make_xt_impulse();
extern float **stim_make_xt_drifting_grating();
extern float **stim_make_xt_jumping_grating();
extern float **stim_make_xt_amc_grating();
extern float **stim_make_xt_sumosin_grating();
extern float **get_xt_stim();

// Frame sequence computation
extern float **compute_fprod_array();
extern float **compute_fprod_array_periodic();
extern float *get_response_from_fprod();
extern void time_domain_response();

// Make frame sets for various stimuli
extern int *sine_pattern(); // Perhaps should be internal
extern void make_frame_set_pos_noise();
extern void make_frame_set_wyphasm_noise();
extern void make_frame_set_phase_con_noise();
extern void make_frame_set_amp_noise();

// Frame sequence for various stimuli
extern void get_frame_sequence_dirmod();
extern void get_frame_sequence_amp();
extern void get_frame_sequence_wyphasm();
