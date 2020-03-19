/*
 *  See bottom for STRUCTURES
 */
extern float su_gasdev();
//extern int gaussian_int();
extern float su_expdev();
extern float su_gaussian_float();

extern void free_sarray();
extern void write_spike_plot();
extern void write_float_spike_plot();
extern float spike_util_ran2();
extern int count_spikes();
extern int count_spikes_float();
extern int count_spikes_dbl();
extern float mean_spike_rate();
extern float mean_spike_rate_float();
extern int total_spikes_sarray();
extern int get_max_time_sarray();
extern int *convert_sampling_units_trial();
extern int **get_sarray_convert_sampling();
extern void sarray_convert_sampling();
extern int get_first_spike();
extern int get_first_spike_in_window();
extern int get_last_spike();
extern int *extract_spikes();
extern int *extract_int_spikes_from_float();
extern int *extract_int_spikes_from_dbl();
extern void add_const_sarray();

extern int count_unique_spikes();
extern int *extract_unique_spikes();
extern void get_unique_sarray();

extern void spikeu_expand_spikes_wsd_accum(); // synaptic depression
extern float *expand_spike_array();
extern int *expand_spike_int_array();
extern float *expand_spike_farray();
extern void get_times_from_expanded_int_spikes();
extern void get_times_from_expanded_float_spikes();
extern int *histogram_expand_sarray();

extern float *merge_spike_arrays_float();
extern int *merge_spike_arrays();
extern int *merge_sarray();
extern int *float_to_int_spike();

extern char *su_compare_trains();

extern void print_sarray();
extern void print_sarray_append();
extern void align_to_first_spike_in_window_sarray();
extern void spike_count_sarrays();
extern void spike_count_r_sarrays();

extern void spike_trig_stim_sarray();
extern void spike_trig_avg_sarray();
extern void spike_trig_avg_sarray_varcond();
extern void spike_trig_avg_condition_sarray();
extern void spike_trig_distrib_sarray();
extern void spike_trig_2d_avg_trial();
extern void second_order_wiener_kernel_sarray();
extern void get_pattern_trig_sarray_from_sarray();
extern int pattern_trig_avg_sarray();
extern int pattern_trig_avg_self();
extern void multi_pattern_trig_avg_sarray();
extern void multi_pattern_trig_avg_farray();
extern void avg_spikes_at_trig_sarray();
extern void spikeu_patt_trig_avg_evcode_sarray();

extern void get_coincidence_spikes();
extern void get_coincidence_sarrays();
extern void get_suppress_spikes_sarrays();
extern void transmission_stat_sarrays();

extern void shift_sarray();
extern void get_filtered_sarray();
extern void get_sc_filtered_double_sarray();
extern int check_trial_condition_su();

extern void get_sarray_from_events();

// vardur - variable duration analysis
extern float *su_vardur_rate_sarray();

extern float mean_spike_count_sarray();
extern float mean_spike_count_sarray_window();
extern float mean_spike_rate_sarray();
extern float *get_norm_spike_count_sarray();
extern void get_cycle_trig_avg_trial();
extern void get_fourier_harmonic_sarray();
extern void get_fourier_harmonic_sarray_old();
extern void increment_psth_float_bin_trial();
extern void psth_float_bin_trial();
extern void psth_float_bin_sarray();
extern void su_psth_float_bin_sarray_vdur();
extern int *simple_pst_sarray();
extern float *psth_prob_sarray();
extern int *simple_psth_bin_sarray();
extern float *adaptive_pst_sarray();

extern void spikeu_sarray_condition_sarray();

extern float *get_first_spike_distrib_sarray();

extern int *simple_isi_sarray();
extern void get_longest_isi_sarray();
extern int get_fraction_isi_sarray();
extern int get_time_fraction_isi_sarray();

extern void power_burst_metric();

extern int get_tuning_time();
extern int get_response_onset();
extern float *shifted_difference_floats();
extern float *distance_trial();

extern float log_prob_spikes();
extern float *log_prob_sarray();
extern float *compare_spike_to_pst();

extern float avg_isi_trial();
extern float *get_trial_intervals();
extern void interval_threshold();
extern void get_intervals_xy_plot();
extern void isi_mean_sdev_skew_trial();
extern float *avg_rate_intervals_trial();
extern float *avg_rate_intervals_sarray();

extern float *interval_cv_sarray();
extern void avg_trial_cv_sarray();
extern void avg_trial_cv2_sarray();

extern void get_gap_sarray();

extern void first_last_count_events_trial();
extern struct event_trial *get_events_sarray();
extern int count_events();
extern float count_weighted_events();

extern void make_poisson_float_spikes();
extern void make_poisson_float_spikes_new();

extern int make_poisson_spikes_int_refract();
extern int make_gaussian_spikes();
extern int make_uniform_spikes();

extern int spike_to_gaussian_burst();
extern int spike_to_uniform_burst();

extern void spike_util_random_select_spikes();

extern void make_random_sarray_from_isi();
extern void make_random_spikes_from_prob();
extern void make_random_spikes_from_prob_refract();
extern void make_random_sarray_from_prob();
extern void make_random_sarray_from_sarray();
extern void make_poisson_float_spikes_from_prob();
extern void make_poisson_f_spikes_from_pr_refr();
extern void make_poisson_spikes_from_prob();
extern void make_poisson_sarray_from_prob();
extern void make_integrate_fire_spikes_from_farray();

extern float *spikeu_poisson_g();
extern float *spikeu_poisson_g_prob();

extern float roc_sarrays();

extern void pick_cross_spikes();
extern void get_stationarity_sarray();

extern void spike_word_recode_sarray();
extern void spike_word_distrib_sarray();
extern void stim_condition_word_distrib_sarray();

extern void spike_util_sta_r_v_stim_sarray();
extern void spike_util_sta_r_v_stim_bayes_sarray();

/*****************************************************************************/
/*                                                                           */
/*                          SPIKE_UTIL_COND_STRUCT                           */
/*                                                                           */
/*****************************************************************************/
struct spike_util_cond_struct {
  char  *ctype;             /* Condition type (XXX = "min" or "max")
			       - abs_val_avg_XXX
			       -     sqr_avg_XXX
			       -         var_XXX
			       -         val_XXX
			    */
  char  *data;              /* Type of data, e.g., "conv" or "stim" */
  int    win0;              /* start time for applying condition */
  int    winn;              /* duration of window */
  float  vcrit;             /* critical value */
  float  vcrit1;            /* critical value */
  int    nsuc;              /* Number of successful matches */
  int    nfail;             /* Number of failures */
  float  fsuc;              /* Fraction of successes */
  float  ftarg;             /* Target success fraction */
};
