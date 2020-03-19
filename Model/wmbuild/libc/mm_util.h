/*** MM_UTIL.H ***/

//
//  DATA TRANSMISSION
//
extern void  mm_data_send_iarray();
extern void  mm_data_recv_iarray();
extern void  mm_data_send_farray();
extern void  mm_data_recv_farray();
extern void  mm_mpi_receive_carray();
extern void  mm_cmd_send();
extern char *mm_cmd_recv();

extern  int  mm_job_get_stimi();
extern void  mm_job_get_indices();

/* struct mm_job_struct *mm_job_get_ptr(); */

//
//  MONITOR and GUI
//

//     void mm_mon_help();
//     void mm_mon_get_lines_height();
//     void mm_mon_init();
//     void mm_mon_update();
//     void mm_mon_update_by_proc();
//     void mm_mon_update_all();
extern void mm_mon_event();

//    char *mm_hms_for_time();
//     void mm_suspend_processor();
extern void mm_gui_update();
extern void mm_gui_init();
//     void mm_gui_help()
extern void mm_gui_event();
extern void mm_gui_hold();

//
//  JOB CONTROL
//
extern  int mm_time_elapsed_min();
extern  int mm_get_proc_state();
extern  int mm_proc_allow_reactivate();
extern  int mm_get_ndone();
extern void mm_print_stat();
extern void mm_job_done();
extern unsigned long mm_current_time_ms();
extern void mm_assign_job();
extern unsigned long mm_get_start_time_for_job();
extern void mm_ja_log_proc_names();
extern void mm_init_job_control();
//      int mm_get_fastest_free_proc();
//      int mm_get_next_avail_rpt();
//   struct mm_job_struct *mm_get_next_job_for_proc()
extern int mm_suggest_next_assignment();

