# Binocular Model for Contrast Reversed Stimuli

This model will produce probabalistic neuronal responses in V1 to an sine wave stimulus in one eye that is contrast reversed in the other eye. For example, the stimulus is an edge going from black to white in one eye, and going from white to black in the other eye. In other words, the stimuli are offset by a phase of 180 degrees in each eye. This matches the central idea of this project (described in the main README file of this repository). However, rather than a complex grating being applie to an image, we're simply looking at the neuronal responses to dichoptic stimuli that are contrast reversed.

## To compile the libraries

download WM file from www.iModel.org

`wget http://www.imodel.org/wm/download.html`


You need 'cc' and 'mpicc' for compiling.
Run the following three commands in this directory.
Note, the 1st line may take ~30 sec to compile the libraries.

Also note, the library mm_util may not compile. If this happens,
you will need to install it manually.

  `source ./s/compile_lib`
  
  `source ./s/compile_main`
  
  `source ./s/compile_main_mm`

Now, run a test to verify that WM is working:

  `./bin/wm mod ./demo/MEO_Gabor.moo ./demo/sine_dir.stm ./demo/d.rsp`


## To run the model: 
### NOTE: the model has already been run and the output file is called `zz.nd`. 
#### This information is just for future runs of different models. Skip down below to view model responses.

`cd wmbuild`

`./bin/wm mod BDE_Gabor.moo sine_phase_disp.stm t_all.rsp tn 256 stim_nrpt 4`

this will write an output file `'zz.nd'`. Then you would use the nData browser (the nd.jar)
file to view the model responses. 

## View Model Responses

You'll need to download the nD browser on the iModel.org website

`wget http://www.imodel.org/nd/down/v1/nd.jar`

* NOTE: if you're on a PC, you will need to download a graphical linux application. To view the model responses using a PC, find instructions here: https://seanthegeek.net/234/graphical-linux-applications-bash-ubuntu-windows/

To view the responses using a Mac or Linux:

`java -jar nd.jar zz.nd`

Once you are viewing the responses, the model is interactive so you can play around with 
different parameters. The most informative result can be found by adjusting the following 
parameters: 

* Set Channel 0 to 'prob'

* Set Analysis to 'Tuning_Curve'

* Set phase_disp to 'all'


This displays the average firing probability for neurons at phase disparieties ranging from 
0-360. As you can see, for binocular neurons that prefer stimuli that have zero disparity, the firing rate is minimal. This indicates that these binocular neurons are compromised for contrast reversed stimuli. However, this model is only useful for showing the response of neurons that prefer a zero disparity stimulus. As we know, there are other neurons that do not have this same preference, so these neurons would likely display normal firing patterns for contrast reversed stimuli in each eye.   

