extern struct param_pair_list *paramfile_ppl_get_init();
extern void paramfile_free_pair();
extern void paramfile_free_ppl();

extern void paramfile_ppl_to_carray();
extern struct param_pair_list *paramfile_carray_to_ppl();

extern void paramfile_print_ppl();

extern int   paramfile_test_param();
extern int   param_test();              // Shorter name for prev.

extern struct param_pair *paramfile_get_next_prefix_pointer();
extern struct param_pair *paramfile_get_first_prefix_pointer();
extern int paramfile_count_prefix_param();

extern struct param_pair *paramfile_get_param_pointer();
extern struct param_pair *paramfile_get_param_pointer_or_exit();
extern struct param_pair *paramfile_get_param_pointer_prev();

extern int  paramfile_list_delete();
extern int  paramfile_list_update();
extern int  paramfile_list_update_index();
extern void paramfile_update_pair_in_list();
extern void paramfile_update_add_pair_to_list();
extern void paramfile_add_pair_to_list();
extern void paramfile_ppl_add_name_val();
extern void *paramfile_include();
extern struct param_pair_list *paramfile_read();
extern void paramfile_update_from_list();
extern struct param_pair_list *paramfile_init_from_list();
extern struct param_pair_list *paramfile_get_initial_params();
extern struct param_pair_list *paramfile_get_initial_params_update_only();

extern int   paramfile_prefix_delete();

extern int   paramfile_get_nth_int_param_ptr_or_exit();
extern int   paramfile_get_nth_int_param_or_exit();
extern int   paramfile_get_int_param_or_exit();

extern float paramfile_get_nth_float_param_ptr_or_exit();
extern float paramfile_get_nth_float_param_or_exit();
extern float paramfile_get_float_param_or_exit();

extern char *paramfile_get_nth_char_param_pp_or_exit();
//extern char *paramfile_get_nth_char_param_ptr_or_exit_DEBUG();
extern char *paramfile_get_nth_char_param_or_exit();
extern char *paramfile_get_char_param_or_exit();

extern int   paramfile_get_int_param_default();
extern float paramfile_get_float_param_default();
extern char *paramfile_get_char_param_default();

extern int   param_geti_dflt();  // Shorter names
extern float param_getf_dflt();
extern char *param_getc_dflt();

extern int   param_geti_exit();  // Shorter names
extern float param_getf_exit();
extern char *param_getc_exit();

extern void  param_geti_s2_exit();  // Params for two stimuli: s1_, s2_
extern void  param_getf_s2_exit();  // Params for two stimuli

extern void  param_getf_2();  // Params for two stimuli, arbitrary names
extern   int paramfile_get_flt_2of3();
extern  void param_getf_lrd();  // Binoc:  foo foo_L foo_R foo_disp

extern char *paramfile_get_comment_for_param();
extern char  paramfile_get_type_from_comment();

extern int   paramfile_count_values_pointer();
extern int   paramfile_count_values_param();
extern int   paramfile_get_slist_pointer();
extern int   paramfile_get_slist_param();

extern char *paramfile_make_fname();  /*** Build outfile names ***/

extern void paramfile_read_diff();


//
// ONODE
//
extern         void  paramfile_onode_free();
extern         void  paramfile_onode_write_file();
extern         void  paramfile_onode_print();
extern struct onode *paramfile_onode_create_onode();
extern         char *paramfile_onode_carray();
extern struct onode *paramfile_get_onode_from_carray();
extern         void  paramfile_onode_insert_next_onode();
extern         void  onode_insert_child_at_end();
extern struct onode *paramfile_onode_get_init_onode();

extern struct onode *onode_copy();

// Searching through the onode tree
extern          int   onode_tree_count();
extern         void   onode_tree_fill_list();
extern struct onode **onode_tree_get_list();

// Reading from file
extern struct onode *paramfile_onode_extract_item();
extern struct onode *paramfile_onode_read();
extern struct onode *paramfile_onode_file_read();
extern struct onode *paramfile_onode_xml_read();
extern struct onode *paramfile_onode_xml_file_read();

// Tags
extern char *onode_get_tag_val_ptr();


// Finding onode type in a list
extern struct onode *onode_get_next_type();

// Get child nodes
extern struct onode *onode_get_item_child();    // Use to check existence
extern          int  onode_item();              // Use to check existence, or
extern struct onode *onode_child_get_unique();
// See "onode_get_node_type_item_val()" below

// Get item params
extern  char *onode_get_nth_val();
extern    int onode_getpar_nvals();
extern    int onode_getpar_int_dflt();
extern    int onode_getpar_int_exit();
extern   int *onode_getpar_int_list();
extern  float onode_getpar_flt_dflt();
extern  float onode_getpar_flt_exit();
extern  float onode_getpar_nth_flt_exit();
extern float *onode_getpar_flt_list();
extern double onode_getpar_dbl_exit();
extern  char *onode_getpar_chr_dflt();
extern  char *onode_getpar_chr_null();
extern  char *onode_getpar_chr_exit();
extern  char *onode_getpar_chr_ptr();
// NOTE: A function like "onode_getpar_chr_ptr_dflt" should not exist
//   becuase you should not get a ptr to a default value
extern  char *onode_getpar_chr_ptr_exit();
extern  char *onode_getpar_child_chr_ptr_exit();
extern char **onode_getpar_chr_list();
extern  char *onode_getpar_unit();

extern int onode_count_otype();

extern struct onode *onode_get_node_type_item_val();  // Should move above?

extern struct onode *onode_get_item_for_ostr();

// Check for the existence of items, or item-value pairs
extern int onode_test_ostr();  // Whether the COMMAND LINE STRING item exists
extern int onode_test_int();   // Whether the item has the INT value
extern int onode_test_chr();   // Whether the item has the CHR value

extern char *onode_make_fname();  // Uses default '*' prefix for wm
extern int onode_update_value();  // Change a value for an ostr
extern int  *onode_update_from_slist();

// Resolving ID and reference tags
extern int onode_conflict_with_children();
extern void onode_copy_contents_novel();
extern void paramfile_onode_resolve();

// XML
extern struct onode *myxml_get_node_tag_with_child_val();
extern char *myxml_get_tag_val();
