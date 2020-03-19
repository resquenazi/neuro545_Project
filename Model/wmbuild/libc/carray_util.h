/*****************************************************************************/
/*                                                                           */
/*  All names should have "carray" or "string" in them.                      */
/*                                                                           */
/*****************************************************************************/

struct strel{ // String Element
  char *s;
  struct strel *prev,*next;
};

struct dll_slist{ // Doubly linked string list
  int n;
  struct strel *first,*last;
};

struct math_strel{ /*** Math string element ***/
  char *etype;      /* Type of element
		       int
		       float
		       op1 */
  float fval;
  int ival;
  char *sval;
  struct math_strel *prev,*next;
};

extern struct dll_slist *dll_string_get_init();
extern void dll_string_free();
extern void dll_string_print();
extern void dll_string_write_file();
extern void dll_string_insert();
extern struct strel *dll_string_search();
extern void dll_string_append();
extern char **dll_string_get_slist();
extern char *dll_string_get_nth_string_ptr();

extern int  string_count_char();
extern void string_make_lower_case();

// Numbers
extern int is_digit_char();
extern int is_numeric_char();
extern int is_name_char();
extern int is_number_string();

extern char ***get_2d_pointer_carray();
extern void free_2d_pointer_carray();
extern void free_2d_pointer_var_count_carray();
extern void free_string_list_carray();

extern void append_null_to_carray();
extern char **copy_2d_carray();
extern char **concat_2d_carrays();
extern int is_const_2d_carray();
extern int is_suffix_string();
extern int compare_prefix_string();
extern int compare_prefix_string_order();
extern void carray_replace_char();
extern char  *replace_html_gtlt();
extern char **replace_html_gtlt_2d_carray();
extern char  *replace_string_segment();
extern char  *replace_string();
extern char **replace_string_2d_carray();

extern char *get_string_from_int();
extern char *get_leading_zero_string_from_int();
extern char *get_bit_string_from_int();
extern void print_bit_string_for_int();

// Used by nda.c, pta_nda
extern void convert_01_string_to_iarray();
extern void convert_n0p_string_to_iarray();
extern void convert_digit_string_to_iarray();
extern void convert_multi_string_to_iarray();

extern void append_string_to_file();
extern int is_int_string();
extern int is_float_string();
extern int count_decimals();
extern int count_items_in_string();
extern void get_items_from_string();
extern char *get_first_item_from_string();
extern void get_iarray_from_string();
extern void get_farray_from_string();
extern float *get_farray_from_carray();
extern int *get_iarray_from_carray();
extern char *make_string_from_items();
extern void check_string_exit_carray();
extern void bubble_sort_2d_carray();
extern int is_unique_sorted_2d_carray();
extern int is_unique_2d_carray();
extern char *carray_get_duplicate_ptr();
extern int search_2d_carray();
extern int search_2d_carray_exit();
extern int search_2d_carray_flag();
extern void update_key_value_pair_carray();

extern void print_carray();
extern void print_2d_pointer_carray();
extern void print_string_list_carray();
extern  int carray_string_read_count();
extern void read_carray();
extern void read_2d_carray();
extern int search_file_2d_carray();
extern void write_carray();
extern void write_carray_segment();

extern void encode_slist_as_carray();

extern void carray_nchar_read();
extern void carray_nchar_write();

extern void get_exclusive_index_compare_2d_carrays();
extern int *get_shared_index_compare_2d_carrays();
extern char **get_column_3d_carray();
extern char **get_unique_slist_from_3d_carray();


extern int *carray_util_parse_underscore();
extern int  carray_util_underscore_int_final();

/*** Process strings of the 'underscore equal' format ***/
extern int carray_get_float_val_underscore_eq();

// Process item lists
extern char *carray_item_list_remove_comment();



/*** Math string processing ***/
/* */
/* struct math_strel *mstrel_get_init(); */
/* void mstrel_free(); */
/* void mstrel_push(); */
/* struct math_strel *mstrel_pop(); */
/* void mstrel_free_stack(); */
/* int mstrel_stack_size(); */
/* int mstrel_etype_n(); */
/* void mstrel_stack_print(); */
/* int mstrel_paren_match_remove(); */
/* void mstrel_calc(); */
/* struct math_strel *get_next_math_string_element(); */
extern float mstr_evaluate();
