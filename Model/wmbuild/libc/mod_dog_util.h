extern void mod_dog_normalize_xyt_filter();
extern void model_dog_01_prep();
extern void model_dog_01_prep_o();

extern float ***mod_dog_gain_ptr();
//extern float ***mod_dog_gain();

extern float ***mod_dog_get_response();
extern float ***mod_dog_get_response_binoc();

extern void model_dog_01_prep_ft();
extern void model_dog_01_done();
extern void model_dog_01_get_response();

extern void mod_dog_get_spikes_2d_flag();

extern void model_dog_simp_01_prep();
extern void model_dog_simp_01_get_response();
extern void model_dog_comp_01_prep();
extern void model_dog_comp_01_get_response();
extern void model_dog_ds_01_prep();
extern void model_dog_ds_01_get_response();

extern void model_dog_ds_02_prep();
extern void mod_dog_get_gt_ds_02();              /* Internal use only */
extern void model_dog_ds_02_get_response();


extern void mod_dog_run_dog_01();
extern void mod_dog_run_simp_01();
extern void mod_dog_run_comp_01();
extern void mod_dog_run_ds_01();
extern void mod_dog_run_ds_02();
