First, you need to compile the libraries. 

Compiling and Running WM

You need 'cc' and 'mpicc' for compiling.
Run the following three commands in this directory.
Note, the 1st line may take ~30 sec to compile the libraries.

Also note, the library mm_util may not compile. If this happens,
you will need to install it manually.

  source ./s/compile_lib
  source ./s/compile_main
  source ./s/compile_main_mm

Now, run a test to verify that WM is working:

  ./bin/wm mod ./demo/MEO_Gabor.moo ./demo/sine_dir.stm ./demo/d.rsp


to run the model: 

cd wmbuild
./bin/wm mod BDE_Gabor.moo sine_phase_disp.stm t_all.rsp tn 256 stim_nrpt 4

this will write an output file 'zz.nd'. Then you would use the nData browser (the nd.jar)
file to view the model responses. 

NOTE: if you're on a PC, you will need to download a graphical linux application
to view the model responses. 

find instructions here: https://seanthegeek.net/234/graphical-linux-applications-bash-ubuntu-windows/

Once you are viewing the responses, the model is interactive so you can play around with 
different parameters. The most informative result can be found by adjusting the following 
parameters: 

Set Channel 0 to 'prob'
Set Analysis to 'Tuning_Curve'
Set phase_disp to 'all'

This displays the average firing probability for neurons at phase disparieties ranging from 
0-360. 

