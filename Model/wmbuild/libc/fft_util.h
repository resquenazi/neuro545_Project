// Manipulation of complex numbers
extern void ab_to_rtheta();
extern void rtheta_to_ab();
extern void rotate_ab();

extern float ***three_d_power();
extern void three_d_fft_prod();
extern void get_power_law_fft();

extern void print_complex_array();
extern float *make_complex_array();
extern float *make_complex_array_from_real_imag();
extern float *real();
extern float *imag();
extern int get_min_power_of_two();
extern void contort_real_farray();
extern void contort_real_2d_farray();
extern void contort_real_3d_farray();
extern void contort_real_3d_farray_space();

extern float **power_2d();
extern float **apply_power_law_2d();
extern float ***apply_power_law_3d();
extern void power_law_2d_transform();

extern void fft();
// void four1();
extern void spctrm();
extern void twofft();
extern void realft();
extern void correl();
extern float *correlate();
extern void convlv();
extern void deconvolve();
extern float *fft_convolve();
extern float *fft_hilbert();
extern float *fft_interpolate();

extern void get_farray_at_resolution();
extern void get_2d_farray_at_resolution();
extern float *fft_smooth_with_gaussian();
extern float **fft_smooth_2d_with_gaussian();
extern float *auto_corr();
extern float *autocorr_farray();
extern float *bandpass_farray();
extern float *hipass_farray();
extern float *gaussian_bandpass_farray();
extern float *dft();
extern float *power_of_complex();
extern float *modulus_of_complex();
extern float *phase_of_complex();
extern float *phase_unwrap();
extern float *power_of_spikes();
extern float *periodogram_normalization();
extern float *periodogram_of_spikes();
extern float *pwrwin_trial_segment();
extern float window_spike_rate();
extern void normalize_pwrwin_by_spike_rate();
extern void normalize_autocorr();
extern void write_pwrwin_record();
extern float *get_padded_spikes();
extern float *auto_spikes();
extern float *auto_convolve_spikes();
extern float *cross_spikes();
extern void write_autocorr_record();
extern void data_trans_inv();
extern void power_farray();
extern void get_mod_phase_farray();
extern void power();
extern void plain_pwrwin();
extern float *pwrwin_data();
extern void avg_pwrwin_rows_2d_farray();

extern void compute_response_single_from_ft();
extern void compute_response_single_from_ft_2();
extern float **compute_response_single_from_ft_ret_s();

extern void optimal_linear_filter_spike_train();

// Adam Kohn's stimulus sequence
extern float *akstim();
