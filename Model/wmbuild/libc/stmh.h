/*****************************************************************************/
/*                                                                           */
/*  stmh.h                                                                   */
/*  Wyeth Bair                                                               */
/*                                                                           */
/*  Data structures stimuli.                                                 */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 2016-2017                                     */
/*                                                                           */
/*****************************************************************************/

/**************************************-**************************************/
/*                                                                           */
/*                                STMH_DOTMOV                                */
/*                                                                           */
/*  Dot movie.                                                               */
/*                                                                           */
/*****************************************************************************/
struct stmh_dotmov{
  int     nframes; // Number of frames
  int    *cnt;     // [nframes] Number of dots per frame
  int   **flag;    // [nframes][cnt[...]] type or condition of dot
  float **x;       // [nframes][cnt[...]] x-coords
  float **y;       // [nframes][cnt[...]] y-coords
  float **z;       // [nframes][cnt[...]] z-coords, for depth
  float **a;       // [nframes][cnt[...]] amplitude, e.g. for black / white
};

/**************************************-**************************************/
/*                                                                           */
/*                                STMH_GPATCH                                */
/*                                                                           */
/*  Gabor patch.                                                             */
/*                                                                           */
/*****************************************************************************/
struct stmh_gpatch{
  int   ptype;     // patch type: 1-gabor, 2-noise
  float x0,y0;     // center
  float dir;       // direction (deg)
  float sf;        // spatial frequency (cyc/deg)
  float tf;        // temporal frequency (Hz)
  float sd_par;    // SD parallel (deg)
  float sd_orth;   // SD orthogonal (deg)
  float ph0;       // Initial spatial phase (deg)
  float amp;       // Amplitude
  int   seed;      // Randomization seed
  int   pixw;      // Noise pixel size (pix)
};
