#define  MAC_DIR_ADMIN "/Users/wyeth/arch/VAR_WWW_HTML/emod/admin"
#define  MAC_DIR_USR   "/Users/wyeth/arch/VAR_WWW_HTML/emod/usr"
#define  MAC_DIR_WORK  "/Users/wyeth/arch/VAR_WWW_HTML/emod/work"
#define  MAC_DIR_WORK_STIM  "/Users/wyeth/arch/VAR_WWW_HTML/emod/work_stim"
#define  MAC_DIR_MM2   "/Users/wyeth/arch/VAR_WWW_HTML/emod/t"
#define  MAC_DIR_M     "/Users/wyeth/arch/WILLIS/public_html/m"
//#define  MAC_DIR_M_MOD "/Users/wyeth/arch/WILLIS/public_html/m/mod/mod"
#define  MAC_DIR_MM    "/Users/wyeth/arch/WILLIS/public_html/m/m"
#define  MAC_DIR_M_RSP "/Users/wyeth/arch/WILLIS/public_html/m/mod/rsp"
#define  MAC_DIR_M_STM "/Users/wyeth/arch/WILLIS/public_html/m/stim"
#define  MAC_DIR_M_UPT "/Users/wyeth/arch/WILLIS/public_html/m/uptab"

// Was:  /var/www/html/emod/...
#define  MAIN_DIR_ADMIN      "/var/www/m/admin"
#define  MAIN_DIR_USR        "/var/www/m/usr/u/00"
#define  MAIN_DIR_WORK       "/home/imod/work"
#define  MAIN_DIR_WORK_STIM  "/home/imod/work_stim"
#define  MAIN_DIR_MM2        "/var/www/m/m"

// Was:  /home/willis/wyeth/public_html/m/...
#define  MAIN_DIR_M       "/var/www/m"
//#define  MAIN_DIR_M_MOD   "/var/www/m/mod/mod"
#define  MAIN_DIR_MM      "/var/www/m/m"
#define  MAIN_DIR_M_RSP   "/var/www/m/mod/rsp"
#define  MAIN_DIR_M_STM   "/var/www/m/stim"
#define  MAIN_DIR_M_UPT   "/var/www/m/uptab"

// For 'mmm'
#define  MAIN_DIR_MF      "/home/wyeth/m/machfile/imod"
#define  MAIN_BIN_MMM     "/home/wyeth/bin/mm_imodel"
#define  MAIN_DIR_LOC_STM "/home/wyeth/m/stm"

// WYETH The first line below (commented out) doesn't work, not sure why,
// perhaps 'www.imodel.org' cannot be decoded by the Java calls that I'm using?
// 
//#define  WWW_M_USR    "http://www.imodel.org/usr/u/00"
//#define  WWW_M_USR    "http://tito.biostr.washington.edu/m/usr/u/00"
#define  WWW_M_USR    "http://wartburg.biostr.washington.edu/m/usr/u/00"

#define  PASSWORD_FILE    "/var/www/m/admin/data/userpwd.txt"

#define  FILE_ADM_USR_RUN    "usr_run.list"
#define  FILE_ADM_JOB_RUN    "job_run.list"
#define  FILE_ADM_JOB_CURR   "job_run.curr"
#define  FILE_ADM_PID_MC     "pid_mc.num"
#define  FILE_ADM_PID_JOB    "pid_job.num"
#define  FILE_ADM_MCTRL_LOG  "mctrl.log"
#define  FILE_ADM_MCTRL_TMP  "mctrl.temp.log"
#define  FILE_ADM_MCTRL_DATE "mctrl.date"
#define  FILE_MOD_MODLIST    "ModelList.txt"
#define  LOGSTR              "MCTRL-LOG-DATE "
#define  FILE_OBS_PARAM      "obspar.txt"

// *** MUST create these at top of "mctr_util.c" ALSO ***
// *** MUST create these at top of "mctr_util.c" ALSO ***
// *** MUST create these at top of "mctr_util.c" ALSO ***
extern char MCTRL_DIR_ADMIN[];
extern char MCTRL_DIR_USR[];
extern char MCTRL_DIR_WORK[];
extern char MCTRL_DIR_WORK_STIM[];
extern char MCTRL_DIR_MM2[];
extern char MCTRL_DIR_M[];
//extern char MCTRL_DIR_M_MOD[];
extern char MCTRL_DIR_MM[];
extern char MCTRL_DIR_M_RSP[];
extern char MCTRL_DIR_M_STM[];
extern char MCTRL_DIR_M_UPT[];

extern         void  mctrl_util_determine_dir_structure();
extern         void  mctrl_util_report_pid_run();
extern struct onode *mctrl_util_modlist_get_onode();
extern         char *mctrl_util_cmd_line_pars();

//
//  UF - User file system
//
extern void  mctrl_util_uf_fi_read();
extern char *mctrl_util_uf_fman_get();
extern char *mctrl_util_uf_fname();
