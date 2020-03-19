extern void zero_iarray();
extern int **get_const_2d_iarray();
extern int **get_zero_2d_iarray();
extern int **get_2d_iarray();
extern int ***get_2d_pointer_iarray();
extern void free_2d_iarray();
extern int ***get_3d_iarray();
extern int ***get_zero_3d_iarray();
extern void free_3d_iarray();

extern int  *copy_iarray();
extern int **copy_2d_iarray();
extern void append_iarray();
extern int sum_iarray();

extern int is_const_iarray();
extern int is_const_col_2d_iarray();
extern int is_nondecreasing_iarray();
extern int is_nonnegative_iarray();

extern int *convert_2dv_to_iarray();

extern void reverse_iarray();
extern void rotate_iarray();
extern char *get_string_from_iarray();
extern int count_occurrences_iarray();
extern void iarray_pattern_search();
extern int get_index_min_search_iarray();
extern int get_index_search_iarray();
extern int get_index_search_nth_iarray();
extern int get_index_search_range_iarray();
extern int diff_iarrays();
extern int get_index_max_diff_iarrays();
extern int get_max_diff_iarrays();
extern int *add_iarrays();
extern void add_const_iarray();
extern int *multiply_iarray();
extern void subsample_iarray();
extern int count_min_runs_middle_iarray();
extern void replace_min_runs_middle_iarray();
extern int *get_column_2d_iarray();
extern int get_decimal_number_from_iarray();
extern void get_bit_iarray_from_int();
extern int count_digits();
extern void print_iarray_80_column();
extern void heap_sort_iarray();
extern void bubble_sort_iarray();
extern int *index_bubble_sort_iarray_double();
extern int *index_bubble_sort_iarray();
extern int *get_permutation_iarray();
extern void permute_iarray();
extern void get_all_permutations_iarray();
extern int convert_multi_index_to_1d();
extern void get_nd_cross_product_list_iarray();
extern void get_nd_cross_product_list_iarray_rev();

extern int *get_unique_iarray();
extern int *get_unique_sort_iarray();
extern void get_unique_count_iarray();

extern void mean_sdev_freq_hist_iarray();
extern void mean_sdev_iarray();
extern int count_non_zero_iarray();
extern int sum_column_2d_iarray();
extern void make_mean_zero_iarray();
extern int *f2iarray();
extern int *f2iarray_round();

extern int min_of_iarray();
extern int max_of_iarray();
extern int max_coord_iarray();
extern int min_coord_iarray();
extern int max_of_2d_iarray();
extern void min_max_2d_iarray();
extern void print_min_max_2d_iarray();  // For debugging

extern void print_2d_iarray();
extern void read_iarray_ptr_until_string();
extern void read_iarray();
extern void write_iarray();
extern void write_iarray_plot();
extern void write_iarray_xy_plot();
extern void append_iarray_xy_plot();
extern void append_iarray_plot();
extern void append_iarray_plot_offset();
extern float lin_interp_iarray();

extern int *pascal_triangle();


