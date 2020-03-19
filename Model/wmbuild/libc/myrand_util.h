extern int yarlrand();
extern void yarlsrand();
extern float myrand_util_ran2();
extern int myrand_get_poisson_count();
extern float myrand_get_poisson_count_mean(float, int, int *);
extern int myrand_poisson_count();
extern void hpsort();
extern void shuffle_iarray();
extern void shuffle_iarray_nyu();
extern void shuffle_iarray_return_seed();
extern int   *get_shuffle_index();
extern int   *get_shuffle_index_return_seed();
extern int   *get_seeds();
extern int  **get_seeds_2d();
extern int ***get_seeds_3d();
extern float *get_random_floats();

extern   int  *myrand_get_std_unif_int_seq();
extern float  *myrand_get_std_unif_float_seq();
extern   int  *myrand_get_std_bin_seq();
extern   int  *myrand_get_std_tern_seq();
extern   int  *myrand_get_std_quad_seq();
extern   int  *myrand_get_std_mseq();
extern float  *myrand_get_std_gauss_seq();
extern float **myrand_get_std_bin_float_2d();
extern float **myrand_get_std_quad_float_2d();

extern char **get_deviate_list_from_description();

// From Greg H's lab
//       int   StimRandShort();
extern   int   ejrand_getnums();
extern float **myrand_get_std_quad_ejrand_2d();
extern   int  *ejrand_get_seeds();


// Mersenne twister download
extern void genrand_init(unsigned long s);
extern void genrand_init_by_array(unsigned long init_key[], int key_length);
extern unsigned long genrand_int32(void);
extern long genrand_int31(void);
extern double genrand_real1(void);
extern double genrand_real2(void);
extern double genrand_real3(void);
extern double genrand_res53(void);

// Routines based on 'genrand' above
void gr_shuffle_iarray();
int *gr_get_shuffle_index();


// WYETH - could try Xorshift in future for speed-up???
