/***  FARRAY_UTIL.H  ***/

extern void print_farray();

extern  float **get_2d_farray();
extern double **get_2d_darray();
extern float **get_const_2d_farray();
extern float **get_zero_2d_farray();
extern void free_2d_farray();
extern void free_2d_farray_ignore_null();
extern float *copy_farray();
extern double *copy_darray();
extern float *unsigned_char_to_farray();
extern float *i2farray();
extern float *d2farray();
extern float **i2farray_2d();
extern void reverse_farray();
extern float *get_regular_sample_farray();
extern float **copy_2d_farray();
extern float **get_transpose_2d_farray();
extern void transpose_2d_farray();
extern float ***copy_3d_farray();
extern float *get_column_2d_farray();
extern float *get_column_2d_double_farray();
extern void get_sub_2d_farray();
extern float ***get_3d_farray();
extern float ***get_zero_3d_farray();
extern void zero_3d_farray();
extern void free_3d_farray();
extern void free_3d_farray_null();
extern float ***get_xyt_from_txy_3d_farray();
extern float ***get_txy_from_xyt_3d_farray();
extern void get_sub_3d_farray();
extern float **get_2d_from_3d_farray();
extern float *get_1d_from_3d_farray();
extern float ****get_const_4d_farray();
extern float ****get_4d_farray();
extern float ****read_4d_farray();
extern float ***read_3d_variable_farray();
extern void free_4d_farray();
extern float *****get_5d_farray();
extern void free_5d_farray();
extern int equal_farrays();
extern void zero_farray();
extern void round_to_int_farray();
extern int *get_round_to_int_farray();
extern float **rotation_matrix_3d();
extern void farray_rotate_coords();
extern void rotate_four_pairs();
extern void rotate_farray();
extern void shift_farray();

extern void get_y_at_x_farray();
//     float interp_mask_circular_w();
extern float *interp_mask_circular_farray();
extern float lin_interp_farray();
extern float interpolate_2d_farray();
extern float *resample_farray();
extern float *over_sample_farray();
extern float *over_sample_farray_interp();
extern void subsample_farray();
extern void subsample_2d_farray();

extern float *get_binary_threshold_farray();
extern void binary_threshold_farray();
extern void get_threshold_crossing_indices_farray();

extern void discretize_farray();
extern void add_const_farray();
extern void add_const_2d_farray();
extern float *add_scale_farrays();
extern float *add_squared_farrays();
extern float *add_farrays();
extern double *add_darrays();
extern void add_to_farray();
extern void accumulate_2d_farray();
extern float *subtract_farrays();
extern double *subtract_darrays();
extern float subtract_integrate_farrays();
extern float **add_2d_farrays();
extern float **subtract_2d_farrays();
extern float ***subtract_3d_farrays();
extern float mean_square_error_farrays();
extern float ***add_3d_farrays();
extern float ***add_3d_tensors();
extern float ***and_3d_tensors();
extern float  **and_2d_farrays();
extern void add_const_3d_farray();
extern void multiply_farray();
extern void multiply_darray();
extern float *get_scale_farray();
extern void multiply_2d_farray();
extern void multiply_3d_farray();
extern float *multiply_farrays();
extern double *multiply_darrays();
extern void multiply_farrays_in_place();
extern float multiply_integrate_farrays();
extern double dot_prod_2d_farrays();
extern float *divide_farrays();
extern float **multiply_2d_farrays();
extern void half_wave_rectify_farray();
extern void absolute_value_farray();
extern void square_farray();
extern void square_darray();
extern void pow_farray();
extern void pow_2d_farray();
extern float *get_exp_farray();
extern float *get_log_farray();
extern float sum_farray();
extern double sum_square_farray();
extern double sum_square_2d_farray();
extern double sum_2d_farray();
extern double sum_3d_farray();
extern double sum_abs_3d_farray();
extern double sum_square_3d_farray();
extern double sum_abs_farray();
extern float *get_sequence_farray();

extern float max_of_farray();
extern float min_of_farray();
extern void get_min_max_farray();
extern int max_index_farray();
extern int get_index_min_level_right_farray();
extern int get_index_max_level_right_farray();
extern int get_index_local_min_right_farray();
extern int get_index_local_min_left_farray();
extern float get_magnitude_range_farray();
extern int get_index_nearest_value_farray();
extern int max_coord_farray();
extern int min_coord_farray();

extern float max_of_2d_farray();
extern float min_of_2d_farray();
extern float get_min_max_2d_farray();
extern void absolute_value_2d_farray();
extern void half_wave_rectify_2d_farray();
extern void max_coord_2d_farray();
extern void min_coord_2d_farray();
extern void min_max_coord_2d_farray();

extern void get_min_max_3d_farray();
extern void get_max_coord_3d_farray();
extern void get_min_max_coord_3d_farray();

extern int *find_extrema_coord_farray();
extern int strictly_increasing_farray();
extern int strictly_decreasing_farray();

extern float *spline_interp_farray();
extern void   spline_max_farray();

extern float *pad_end_farray();
extern float *diff_farrays();
extern double *diff_darrays();
extern float roc_farrays();

extern void write_xview();
extern void write_xview_transpose();
extern void read_xview();
extern void write_xview_movie();
extern void write_xview_movie_tensor();
extern void write_hips_movie();
extern float ***convert_2d_to_tensor();
extern float ***convert_3d_to_tensor();

extern void get_fractional_change_xy_farray();
extern void write_farray_sdev_plot();
extern void append_farray_sdev_plot();
extern void append_farray_sdev_plot_offset();
extern void write_farray_xy_plot();
extern void write_farray_xy_plot_d7();
extern void write_farray_xy_seq_plot();
extern void write_farray_xy_sdev_plot();
extern void append_farray_xy_sdev_plot();
extern void append_farray_xy_plot();
extern void append_farray_xy_plot_nonew();
extern void append_farray_xy_plot_offset();
extern void append_farray_xy_seq_plot();
extern void write_farray_plot();
extern void write_farray_plot_name();
extern void write_farray_scale_plot();
extern void write_farray_plot_offset();
extern void append_farray_plot_offset();
extern float *get_histogram_farray();
extern float *get_histogram_auto_farray();
extern void write_farray_hist_plot();
extern void append_farray_plot_xscale();
extern void append_farray_plot();
extern void append_farray_plot_color();

extern void free_farray_plot_file_old();
extern void read_farray_plot_file_old();
extern void read_farray_ptr_until_string();
extern int read_farray();
extern int read_2d_farray();
extern int read_2col_farray();
extern void write_farray();
extern void append_farray();
extern void append_line_farray();
extern void write_farray_format();
extern void append_farray_format();
extern void write_2col_farray();
extern void append_2col_farray();
extern void write_3col_farray();
extern void append_2d_farray();
extern void write_2d_farray();
extern void add_farray_to_file();

extern float mean_abs_farray();
extern float mean_farray();
extern void mean_adev_farray();
extern void mean_sdev_farray();
extern void mean_sdev_darray();
extern void mean_sdev_skew_farray();
extern float skew_farray();
extern void weighted_mean_sdev_farray();
extern void mean_two_sided_sdev_farray();
extern void subtract_mean_farray();
extern float *get_zero_mean_farray();
extern void make_max_const_farray();
extern void scale_mean_to_value_farray();
extern void z_score_farray();
extern void z_score_exact_farray();
extern void mean_sdev_freq_hist_farray();

extern void partial_correlations_farrays();
extern void vector_average_farray();

extern void centroid_circular_variance_polar_farray();

// Two-dimensional farrays
extern void mean_2d_farray();
extern void mean_3d_farray();
extern void single_mean_sdev_2d_farray();
extern void single_mean_3d_farray();
extern void single_mean_sdev_3d_farray();
extern void subtract_mean_2d_farray();
extern void mean_sdev_2d_farray();
extern void mean_sd_trial_2d_farray();
extern float *avg_1st_dim_2d_farray();
extern void mean_sdev_col_2d_farray();
extern void mean_sdev_2d_flag_farray();
extern void mean_sdev_2d_segment_farray();
extern void weighted_mean_sdev_2d_farray();
extern void mean_sdev_diff_2d_farray();
extern void weighted_mean_sdev_diff_2d_farray();
extern void make_max_const_2d_farray();
extern void norm_01_2d_farray();
extern void norm_01_3d_farray();

extern float **rect_to_polar_2d_farray();
extern float **polar_to_rect_2d_farray();

extern void norm_area_farray();
extern void norm_area_2d_farray();
extern void norm_positive_area_farray();
extern void norm_max_farray();
extern void norm_variance_farray();
extern void norm_variance_quiet_farray();
extern float *get_cumulative_farray();

extern void get_fourier_harmonic_farray();
extern void get_fourier_harmonic_trial_2d_farray();

extern void sort_farray();
extern int *index_bubble_sort_farray();
extern void farray_sort_carray();
extern void sort_three_farray();
extern int *index_sort_farray();

extern float gaussian();
extern void discrete_gaussian();
extern void fit_to_gaussian();
extern void fit_to_gaussian_4();
extern void fit_to_flat_line();
extern void fit_to_sigmoid();
extern void fit_to_sine();
extern void fit_to_decay_exp();

extern float *convolve_with_mask_renorm();
extern float *convolve_with_mask();
extern float *convolve_with_mask_causal();
extern float *convolve_with_mask_causal_x0();
extern float *convolve_circular_with_mask();
extern float *convolve_sparse_with_mask();
extern float *smooth_with_gaussian();
extern float *smooth_with_gaussian_circular();
extern float *smooth_with_gaussian_omit_center();
extern float *smooth_circular_with_gaussian();
extern float *smooth_sparse_with_gaussian();
extern float **smooth_2d_with_gaussian();
extern float **smooth_2d_with_gaussian_circular();
extern float **smooth_sparse_2d_with_gaussian();
extern float ***smooth_3d_2d_with_gaussian();

extern float *derivative_simple_farray();
extern float *derivative_farray();
extern float *derivative_smooth_farray();
extern float **derivative_smooth_2d_farray();
extern float *second_derivative_farray();

extern float *get_smooth_abs_deriv_xy_data();

extern float *get_paradiagonal_2d_float();
extern float *get_band_diagonal_2d_float();
extern float *get_paradiag_projection_2d_float();
extern float **project32();
extern float ***project43();

extern void width_at_frac_height();
extern void width_at_frac_height_interp();
extern int find_tallest_middle_peak();
extern void find_nl_bin_index();
extern void width_at_frac_height_check_double();
extern void first_left_frac_rise();
extern void first_left_frac_rise_interp();
extern void last_left_frac_rise_interp();
extern void width_at_frac_height_circ_farray();
extern void width_at_frac_area();
extern void estimate_area_peak();
extern int first_left_min_below_crit_farray();
extern int find_peak_start();
extern int first_zero_before_peak();
extern float first_frac_before_peak_interp();

extern float get_entropy_farray();

/*** for PCA principle component analysis SVD ***/
extern float **covar_matrix_2d_farray();
extern void get_eigenvectors_covar_farray();

// Matrix

extern void    matrix_print();
extern float **matrix_transpose();
extern float  *matrix_times_vector();
extern float **matrix_multiply();
extern float **matrix_invert_3x3();

// For Debugging
extern void print_min_max_farray();
extern void print_min_max_2d_farray();
extern void print_min_max_3d_farray();

// Specific and pecuilar
extern float *distrib_box_match();
