#
#  MEO_Gabor
#

mod_type me_gabor_01    # Motion energy

sscale         0.1      # degr/pix
tscale         0.002    # sec/frame
xn            32        # filter x size
yn            32        # filter y size
tn          4096        # filter t size

#
#  Oriented linear filter
#
<filter>
  type           Gabor             # type of filter
  direction      0       (deg)     # preferred direction
  sf             1.25    (cyc/deg) # Spatial frequency
  ssd            0.18    (deg)     # SD space
  tf            10.0     (cyc/s)   # Temporal frequency
  tsd            0.025   (s)       # SD time
  tilt_type      none              # ["none"] "rotate", "shift"
  tilt_frac      0.0               # Amount of rotation, 0.0-none, 1.0-vel.
  write_filter   0                 # 0-none, 1,2,3,4 or 5=all
  opponent       yes               # ["yes"] "no"
</filter>

#
#  Spike Generation
#
<spike_gen>
  type          poisson     # Spike generation algorithm, dflt "poisson"; "ifc"
  offset0       0.0         # Add to filter output, BEFORE scaling (0.0)
  scale         0.0004      # Multiply filter output (1.0)
  offset        0.0         # Add to filter output, AFTER scaling (0.0)
  toffset       0.060   (s) # Time delay added to spike times
  spike_dump    0           # Dump firing probability to 'zz.dump.pl' if 1
</spike_gen>
