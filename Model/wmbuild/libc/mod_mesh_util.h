extern int    mod_mesh_get_mm_mosaic_cn();
extern int   *mod_mesh_ptr_mm_mosaic_cid();
extern float *mod_mesh_ptr_mm_mosaic_cx();
extern float *mod_mesh_ptr_mm_mosaic_cy();
extern float  mod_mesh_get_mm_mosaic_dens_mm();
extern float  mod_mesh_get_mm_mosaic_dens_deg();
extern float  mod_mesh_get_mm_mosaic_degpmm();
extern float *mod_mesh_get_mm_resp();
extern void   mod_mesh_get_rgc_conn();

//     void model_mesh_cone_mosaic_prep();
//     void mod_mesh_ad1_reset();
extern void model_mesh_rgc_01_prep();
extern void model_mesh_rgc_01_done();

//     void mod_mesh_spike_gen();

extern float ***model_mesh_rgc_01_compute_3d_response();
extern float ***model_mesh_rgc_01_get_3d_response();

extern void mod_mesh_compute_diff_all();
extern void mod_mesh_get_spikes_2d_flag();
extern void mod_mesh_get_spikes_3d_flag();

extern void model_mesh_rgc_01_get_response();
extern void model_mesh_mag_01_get_response();

extern void mod_mesh_run_mesh_rgc_01();
extern void mod_mesh_run_mesh_mag_01();
