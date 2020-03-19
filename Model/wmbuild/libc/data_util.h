extern int lookup_key_field_file();

// Packed color data [tcode = 11]
extern   int data_util_t11_int_for_rgb();
extern  void data_util_t11_get_rgb();
extern  void data_util_t11_get_rgb_float();
extern  void data_util_t11_get_rgb_01();
extern  void data_util_t11_get_min_max_rgb_all();
extern int **data_util_t11_2d_from_rgb();

// 3rgb to 3d
extern float ***data_util_3rgb_to_3d();
extern float ***data_util_3rgb_raw1_read();

//  .fst - frame set file
extern void     data_util_fst_write();
extern float ***data_util_fst_read();
extern float ***data_util_fst_read_txt();

extern void write_2d_data();
extern void read_2d_data();
extern void write_3d_data();
extern void write_3d_data_part();
extern void read_3d_data_old();
extern void read_3d_data();
extern  int count_columns_column_data();
extern void read_column_data();

// Van Hateren image database
extern int **read_vanhateren();


// PPM/PGM/PBM image format routines
extern  int read_ppm_image_type();
extern void read_ppm_header();
extern void read_pgm_data();
extern void read_pgm_p5_data();
extern void write_pgm_data_float();
extern void write_ppm_6_data();
extern void write_ppm_6_rgb_data();
extern void write_ppm_6_3d_gray_frames();

// TIFF

// extern struct wytif_struct  *wytiff_get_empty_struct();
extern void wytiff_free();
//extern                void   wytiff_print_struct();
//extern                char  *wytiff_read_string();
//extern                void   wytiff_read_ifd_entry();
//extern                void   wytiff_read_ifd();
//extern              float  **wytiff_read_data_float();

extern struct wytif_struct *wytiff_read();  // Read tiff data from file
extern int wytiff_image_count();
extern int wytiff_check_header();


// XML
extern char *myxml_get_tag_from_string();
extern char *myxml_get_tag_from_string_adv_index();
extern char *myxml_get_tagend_from_string_adv_index();
extern int   myxml_get_tag_assignments();
extern char *myxml_get_value_from_string();
