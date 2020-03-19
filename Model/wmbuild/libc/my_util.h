/***  MY_UTIL.H  ***/

#define SLEN 256
#define SLEN3 384
#define SLEN5 512
#define LONG_SLEN 2048
#define LONG_SLEN4 4096
#define SLEN_MAX 16384
#define SLINE_MAX 131072

#define myfree(ptr) free((char *)(ptr))

extern int mymini();
extern int mymaxi();
extern float myminf();
extern float mymaxf();
extern float min_of_three();
extern float max_of_three();
extern void exit_error();
extern void ee();
extern void *myalloc();
extern float *get_farray();
extern double *get_darray();

extern FILE *my_fopen();

//
//  Log file
//
extern int mylog();
extern void mylog_exit();
extern void mylogx();
extern void mylog_xci();
extern void mylog_ival();
extern void mylog_fval();
extern void mylog_cval();


extern int my_rint();
extern int my_rint_double();
extern double my_log2();
extern int power_of_two();
extern int my_sign();
extern float n_choose_m();
extern int count_ones_binary_int();
extern float get_binary_entropy();

extern void reverse_bytes();
extern int rev_fread();
extern int revflag_fread();
extern int revflag_fwrite();

extern char **get_2d_carray();
extern void free_2d_carray();

extern float *get_zero_farray();
extern float *get_const_farray();

extern int *get_zero_iarray();
extern int *get_const_iarray();
extern void check_int_range();

extern void remove_file();
extern int is_dir();
extern int is_file();
extern char **filelist();
extern char **make_filelist();

extern int file_exists();
extern int file_not_empty();
extern void create_file();

extern int process_exists();

extern int my_true_test();

//
//  For plotting
//
extern double smart_div_size();

//
//  Uncomment this for debugging with 'dmalloc', and see 'wm.c'
//
//#include "dmalloc.h"
