/**************************************-**************************************/
/*                                                                           */
/*   wytif.h                                                                 */
/*   Wyeth Bair                                                              */
/*                                                                           */
/*   For TIFF file format utilities.                                         */
/*                                                                           */
/*C  (C) Copyright  Wyeth Bair 1993-2010                                     */
/*                                                                           */
/*****************************************************************************/

//
//  For TIFF Revision 6.0 file format, see:
//    http://partners.adobe.com/public/developer/en/tiff/TIFF6.pdf
//

struct wytif_struct {
  int w;               // Image width, Number of columns (Tag 256)
  int h;               // Image height (or length), Number of rows (Tag 257)
  int bitpsamp;        // Bits per sample
  int compression;     // 1-none, 5-LZW, 6-JPEG
  int photo_interp;    // 0-white is zero, 1-black is zero
  int orientation;     // Ori of image w.r.t. rows/cols, 1..8
                       //   1-0th row is visual top, 0th col is visual left
                       //   2-0th row is visual top, 0th col is visual right
                       //   ...
  int sppix;           // Samples per pixel (typ. 1-binary, 2-graylev, 3-RGB)
  int planarconfig;    // Planar configuration, 1-book interleaved, 2-not
  char *descrip;       // Image description
  int nstrip;          // Number of strips
  unsigned int *soff;  // Strip offsets [nstrip]
  int row_p_strip;     // Rows per strip (except possibly last strip)
  // Assuming 'h' rows, then there are h % row_p_strip rows in last strip
  int *sbyte;          // Strip byte counts, after any compression [nstrip]
  double xres;         // X resolution; pix/ResUnit
  double yres;         // Y resolution; pix/ResUnit
  int resunit;         // 1-none, 2-inch, 3-cm (default 2)
  char *software;      // Software name/version
  char *datetime;      // DateTime, "YYYY:MM:DD HH:MM:SS" + NULL
  int predictor;       // Predictor, Used w/ LZW Compression
  int newsubfiletype;  // New Subfile Type, the general kind of data
  float **data;        // Image data [w][h]

  unsigned int ifd_next;      // Offset of next IFD, or 0 if none
  struct wytif_struct *next;  // Pointer to next IFD, or NULL
};
