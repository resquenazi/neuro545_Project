
extern void kernu_norm_2d_filter();
extern void kernel_util_norm_xyt_filter();

extern float *rectangle();
extern float *decay_exp();
extern float *doe();
extern float *alpha();
extern float *doa();

//  1D filters
//
//     float  *kernu_o_filter_boxcar();
//     float  *kernu_o_filter_exp();
//     float  *kernu_o_filter_alpha();
//     float  *kernu_o_filter_gaussian();
//     float  *kernu_o_filter_maxwell();
//
//     float **kernu_o_filt2d_gaussian();    // NEW STYLE
//
extern float  *kernu_o_filter();
extern float **kernu_o_filt2d();

extern float *wa_temporal();
extern float *ab_temporal();
extern float *cwq_temporal(int,float,int,float,float,float,float);
//     float *dd_exp_temporal();
//     float *filt_01_temporal();

extern float ***random_filter();
extern float  **gaussian_2d();
extern float  **gauss_2d_noise_raw();
extern float  **gabor_2d_space_raw();
extern float   *gabor_2d_space_raw_irreg();
extern float  **gauss_2d_space_raw();
extern float  **disk_2d_space_raw();
extern float  **dog_2d_space_raw();
extern float  **gabor_2d_space();
extern float ***gabor_space_time();
extern float ***gabor_space_time_tensor();
extern float ***gabor_space_time_tensor_curve();

//       float *kernu_tshift_pow();
//       float *kernu_tshift_sig();
//       float *kernu_tshift();
//     float ***gabor_space_time_tshift();
//         void rotate_space_time();

extern float ***gabor_tilt_tensor();
extern float ***gabor_tilt2_tensor();
extern float ***gabor_tilt3_tensor();
extern float ***gabor_space_time_sep_tensor();
extern float ***gauss_space_exp_time_tensor();
extern float ***exp_time_tensor();
extern float ***delta_tensor();

extern void get_opponent_random();
extern void get_opponent_random_x();
extern void get_quad_opponent_gabor();
extern void get_quad_opponent_gabor_curve();
extern void get_quad_opponent_gabor_tshift();
extern void get_arb_opponent_gabor();
extern void get_quad_opp_tilt_gabor();
extern void get_quad_gabor();
extern float ***get_ppl_gabor();

extern float ***dog_space_time_tensor();
extern float **separable_space();
extern float ***stf3d_01();
extern float ***stf3d_02();
extern void quad_opponent_causal_stf3d_01();
extern void quad_opponent_causal_stf3d_01_tensor();
extern void bde_cwq_stf3d_tensor();
//extern void bde_cwq_stf3d_tensor(struct onode *,int,int,int,float,float,
//float ****,float ****,float ****,float ****);
extern float **get_xt_dog();
