extern                void ifc_util_rm_dump_file();
extern                void ifc_util_free_ifc_param();
extern   struct ifc_param *ifc_util_get_params();
extern   struct ifc_param *ifc_util_get_param_o();
extern   struct ifc_param *ifc_util_get_param_poiss_o();
extern struct poiss_param *ifc_util_get_param_poisson();


// void ifc_rkck();
// void ifc_rkqs();

// void ifc_odeint();
// void test_derivs();
// void ifc_add_noise();
// void ifc_add_noise_o();

extern void ifc_test();

//   float *get_gt_from_sarray();
//     void mod_ifc_01_prep();


extern void model_ifc_01_get_response();
extern void model_ifc_02_get_response();
extern void model_ifc_02_02_get_response();
extern void model_ifc_03_get_response();

// Poisson spikes
extern void ifc_util_poisson();


extern void ifc_util_run_ifc_01();
extern void ifc_util_run_ifc_02();
extern void ifc_util_run_ifc_02_02();
extern void ifc_util_run_ifc_03();
