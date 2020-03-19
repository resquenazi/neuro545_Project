/*****************************************************************************/
/*                                                                           */
/*  mm.h                                                                     */
/*  Wyeth Bair                                                               */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/

// Tags for MPI communication
#define tag_mform      99
#define tag_mppl      100
#define tag_sppl      101
#define tag_rppl      102
#define tag_procname  110
#define tag_runtrial  200
#define tag_result    201  // The output of simulation is sent from slave

// ADDED 2013, Nov, to prevent delays for very large response data transfers
#define tag_req2send  204  // Master says it is OK for slave to send results
#define tag_ok2send   205  // Master says it is OK for slave to send results


// General command
#define tag_cmd_str  500  // Send a string command
#define tag_data_int 501  // Send one or more ints
#define tag_data_flt 502  // Send one or more floats


// Fitting
#define tag_fit_cmd  300  // Send a command code (int) to slave
                          //   1-run, 2-get ints, 3-get floats, 4-get chars
#define tag_fit_str  301  // Send a name string
#define tag_fit_dvi  311  // Data vector int
#define tag_fit_dvf  312  // Data vector float
#define tag_fit_dvc  313  // Data vector char

/**************************************-**************************************/
/*                                                                           */
/*                               MM_JOB_STRUCT                               */
/*                                                                           */
/*****************************************************************************/
struct mm_job_struct {
  int mi;                        // Model index
  int stim;                      // stimulus number
  int rpt;                       // repeat number
  int trial;                     // trial, 0...stim*npt-1; -1 until finished
  struct mm_proc_struct *proc;   // Processor if assigned, or NULL
  unsigned long t1;              // Start time, -1 until started
  unsigned long t2;              // Finish time, -1 until finished
  int tn;                        // t2 - t1, -1 until finished
};
/**************************************-**************************************/
/*                                                                           */
/*                               MM_JEL_STRUCT                               */
/*                                                                           */
/*****************************************************************************/
struct mm_jel_struct { // Job list element
  struct mm_job_struct *j;       // Pointer to job in job table
  struct mm_jel_struct *next;    // Next job list element
};
/**************************************-**************************************/
/*                                                                           */
/*                               MM_PROC_STRUCT                              */
/*                                                                           */
/*****************************************************************************/
struct mm_proc_struct {
  char *name;   // Processor name
  int id;       // Processor ID number
  int state;    // State:  -1 unknown
                //   0 idle, available
                //   1 busy
                //   2 deactivated by user
                //   3 deactivated by user, but results recieved

  struct mm_job_struct *job;       // Current job
  
  int ndone;                      // Number of completed jobs
  struct mm_jel_struct *done;     // Head of list of completed jobs
  struct mm_jel_struct *donelast; // Tail of list of completed jobs
};
