
// fundamental
extern int stm_get_tsamp();
extern int stm_get_ti0();
extern int stm_get_t1i();
extern int stm_get_n();
extern int stm_util_count_new_frames();


// general
extern void stm_get_ran_seq();
extern void stm_get_ran_seq_2d();

extern void stm_util_grid_coords();
extern void stm_util_grid_coords_1d();

// shape
extern void stm_shape_fill_recursive();
extern void stm_shape_draw_contour();
extern void stm_shape_draw_contour_antialias();
extern void stm_shape_draw_fill();

// noise
extern float **stm_noise_get_grid();

// vanhateren
extern float **stm_vanhat_get_patch();

// ranstep
extern float *stm_ranstep_get_stim();

// barnoise
extern float **stm_barnoise_get_template();

// barray
extern void stm_barray_get_bar_centers();
extern int **stm_barray_get_seq();

// stat1d
extern void stm_stat1d_get_seq();

// dirmod dm_test
extern float **stm_dirmod_get_test_pat();
extern void stm_dirmod_test_seq();

// wy4st
extern void stm_wy4st_get_seq();

// gabor patch
extern struct stmh_gpatch *stm_gpatch_get();
extern void stm_gpatch_config_simp();

// dots
extern void stm_get_coherence_dot_movie();
extern void stm_mask_dot_movie();
extern void stm_dot_movie_shift();
extern void stm_dot_movie_copy();

extern struct stmh_dotmov *stm_dotmov_init();
extern               void  stm_dotmov_get_parts();
extern struct stmh_dotmov *stm_dotmov_get_deep_f0();
extern               void  stm_dotmov_get_deep();
extern                void stm_dotmov_get_planes();

// star
extern void stm_star_get_trig();

// pasupathy_shape
extern   void stm_pasupathy_shape_idrot_list();
extern   void stm_pasupathy_shape_get_control();
extern   void stm_pasupathy_yasm_get_control();
extern   void stm_pasupathy_yosh_get_control();
extern   void stm_pasupathy_shape_get_inside_point();
extern double stm_pasupathy_spline_eqn();
extern   void stm_pasupathy_shape_get_coords();
extern   void stm_pasupathy_shape();

// gabor_britten
extern void stm_gabor_britt();
extern void stm_gabor_britt_timing();
