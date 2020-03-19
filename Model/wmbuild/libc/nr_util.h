/***  NR_UTIL.H  ***/

extern void nrerror();
extern void free_f3tensor();
extern float ***f3tensor();
extern void get_f3tensor_from_data();
extern void free_matrix();
extern float **matrix();
extern float *vector();
extern void free_vector();
extern int *ivector();
extern void free_ivector();

extern void avevar();
extern float erfcc();
extern float gammln();
extern float betacf();
extern float betai();

extern float mybetacf();
extern float mybetai();

extern void crank();
extern float ran1();
extern float nr_util_ran2();
extern float nr_util_gasdev();

extern void sort();
extern void sort2();
extern void hpsort_int();
extern float gasdev();
extern void pearsn();
extern void simple_pearsn();
extern void simple_pearsn0();
extern void spear();
extern void gser();
extern void gcf();
extern float gammq();
extern void cntab1();
extern void cntab2();
extern void fourn();
extern void my_fourn();
extern void rlft3();
extern void my_rlft3();

extern void fit();

extern float pythag();
extern void tridag();
extern void eigsrt();
extern void jacobi();
extern void tred2();
extern void tqli();

extern void mrqcof();
extern void gaussj();
extern void covsrt();
extern void mrqmin();
extern void fgauss();
extern void fgauss4();
extern void fsin4();
extern void fexp3();
extern void sigpsych();
extern void flat_line();

extern void my_ttest();
extern void ttest();
extern void tutest();

extern void spline();
extern void splint();
