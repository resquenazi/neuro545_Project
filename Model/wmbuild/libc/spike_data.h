/*****************************************************************************/
/*                                                                           */
/*  wyeth bair                                                               */
/*  caltech                                                                  */
/*  02/27/93                                                                 */
/*                                                                           */
/*  This file contains data structures used for spike train processing.      */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/

/**********************************************************************/
/*                                                                    */
/*  MULTI CHANNEL SPIKE DATA                                          */
/*                                                                    */
/**********************************************************************/
struct mchan_struct {
  char    name[SLEN];   /* channel name */
  int     num;          /* USED BY OLD PROGRAMS */
  int     n;            /* number of spikes recorded on this channel */
  int    *s;            /* spike time, msec */
};

struct multi_trial {
  int     number;       /* trial number in original file */
  int     channels;     /* channels recorded on this trial */
  int     coherency;    /* stimulus coherency, tenths of percent */
  int     direction;    /* stimulus direction, if binary: 0=null, 1=pref'd */
  int     response;     /* animal response, 0=Incorrect, 1=Correct */
  struct mchan_struct *chan;  /* number of spikes per channel */
};

struct multi_channel {
  char    name[SLEN];   /* cell name, usually */
  int     ncomments;    /* number of comment lines */
  char  **comments;     /* comments about this recording site */
  char    mode[SLEN];   /* "dseries", "cseries", "single" added 1/95 */
  int     trials;       /* number of trials */
  int     channels;     /* number of channels recorded per trial */
  float   sampling;     /* sampling rate in Hz */
  int     start;        /* trial start time, msec */
  int     duration;     /* trial duration, msec */
  struct multi_trial *trial; /* array of trial data */
};

/**********************************************************************/
/*                                                                    */
/*  PST:  Post-stimulut Time Histogram                                */
/*                                                                    */
/**********************************************************************/
struct pst_struct {
  char   name[SLEN];
  float  sampling;           /* samples per second */
  int    trials;             /* trials accumulated */
  int    total_counts;       /* total counts */
  int    binsize;            /* in sampling units */
  int    rstart,rduration;   /* response period, msec */
  int    period;             /* duration in sampling units */
  int    n;                  /* number of bins */
  float *bin;                /* total spikes per bin */
};

/**********************************************************************/
/*                                                                    */
/*  ISI:  Inter-spike Interval Histogram                              */
/*                                                                    */
/**********************************************************************/
struct isi_struct {
  char   name[SLEN];
  int    coherency;
  int    direction;
  float  sampling;       /* samples per second */
  int    trials;         /* number of trials averaged */
  int    interval_total; /* total number of intervals */
  int    time_total;     /* total number of msec */
  int    n;              /* number of bins */
  int   *bin;            /* counts for each interval length */
  float *fraction;       /* density distribution of intervals */
};

/**********************************************************************/
/*                                                                    */
/*  PSERIES:  for power spectra, auto and cross correlations          */
/*                                                                    */
/**********************************************************************/
struct p_struct {
  char   name[SLEN];
  int    coherency;
  int    direction;
  int    weight;      /* related to quantity of spikes in average */
  int    n;           /* length of data */
  int    count;       /* number of records in average */
  float *mean;        /* mean pwrwin value */
  float *sdev;        /* standard deviation of pwrwin value */
};

struct pseries_struct {
  char   name[SLEN];
  float  sampling;        /* samples per second */
  int    max_lag;         /* maximum lag time for display */
  int    start;           /* start time */
  int    duration;        /* maximum valid duration */
  int    m;               /* half-width of data window */
  int    k;               /* use 2*k windows */
  int    n;               /* number of values in series */
  int   *svalue;          /* parameter value */
  struct p_struct **rec;  /* spectrum data for this svalue */
};

/**********************************************************************/
/*                                                                    */
/*  AUTOCORRELATION                                                   */
/*                                                                    */
/**********************************************************************/
struct ac_struct {
  char   name[SLEN];
  int    time;            /* start time of window */
  int    count;           /* number of records averaged */
  float  max_freq;        /* max in power spectrum of this acg */
  int    spk_squared;     /* average of squares of spike rates */
  int    n;               /* number of data points in acg */
  float *sdev;            /* standard deviations */
  float *data;
};

/**********************************************************************/
/*                                                                    */
/*  JPSTH:  Joint Post-Stimulus Time Histogram                        */
/*                                                                    */
/**********************************************************************/
struct jpsth_struct {
  char   name[SLEN];
  float  sampling;           /* samples per second */
  int    trials;             /* trials accumulated */
  int    total_counts;       /* total counts */
  int    n;                  /* number of bins */
  struct pst_struct *pst1;   /* */
  struct pst_struct *pst2;   /* */
  short int **bin;           /* coincidences [-32768,32767] */
};

/**********************************************************************/
/*                                                                    */
/*  EVENTS:  This structure holds data about "events" (bursts).       */
/*                                                                    */
/**********************************************************************/
struct event_struct {
  int n;   /* number of spikes in this event */
  int *t;  /* times of spikes in this event */
};

struct event_trial {
  int nev;                  /* number of events in this trial */
  struct event_struct *ev;  /* array of events */
};
