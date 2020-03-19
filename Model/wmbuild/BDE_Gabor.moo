#
#  BDE_Gabor.moo
#

mod_type binoc_filter   # Binocular energy filter

sscale         0.1      # degr/pix
tscale         0.002    # sec/frame
xn            32        # filter x size
yn            32        # filter y size
tn           128        # filter t size

phase_shift   0   (deg) # Phase shift of right eye RF
phase_1       0   (deg) # Phase of cosine for first left filter
binoc_nonlin  halfsq    # Type of nonlinearity
right_sign  ++++        # Sign of R.Eye signal (Eq 18 Read,Parker,Cumming 2002)
                        # A 4-char combination of + and -, e.g., ++-- [++++]
simp_rect     0         # 1-rectify simple cells, 0-do not rectify [0]
simp_thresh   0.0       # Subtracted from all simple monoc responses
                        #   before rect. (Fraction of mean filter area) [0.0]

#
#  Oriented linear filter
#
<filter>
  type           Gabor             # type of filter
  direction      0       (deg)     # preferred direction
  sf             1.25    (cyc/deg) # Spatial frequency
  ssd            0.18    (deg)     # SD space
  tf 		 0.0	 (cyc/s)   # Temporal frequency
  tsd            0.015   (s)       # SD time
  tilt_type	 none	 	   # ["none"] "rotate" "shift"
  tilt_frac	 0.0		   # Amount of rotation, 0-none, 1-vel.

  write_filter   0                 # 0-none, 1,2,3,4 or 5=all
  opponent	 no
</filter>

#
#  Spike Generation
#
<spike_gen>
  type          poisson     # Spike generation algorithm, dflt "poisson"; "ifc"
  offset0       0.0         # Add to filter output, BEFORE scaling (0.0)
  scale         0.0003      # Multiply filter output (1.0)
  offset       -0.1         # Add to filter output, AFTER scaling (0.0)
  toffset       0.060   (s) # Time delay added to spike times
  spike_dump    0           # Dump firing probability to 'zz.dump.pl' if 1
</spike_gen>


#
#  Responses that can be requested in the .rsp file for 'binoc_filter' are
#
#    spikes
#    poisson_prob
#
#    filter_left_1
#    filter_left_2
#    filter_right_1
#    filter_right_2
#
#    binoc_1_pos
#    binoc_1_neg
#    binoc_2_pos
#    binoc_2_neg
#
#    binoc_1_pos_sq
#    binoc_1_neg_sq
#    binoc_2_pos_sq
#    binoc_2_neg_sq
#
#    binoc_1
#    binoc_2
#
#    binoc_energy
#
