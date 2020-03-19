# Neuro545 Project

## Background
Degenerative eye diseases are a leading cause of blindness in adults. While only limited treatment options for diseases such as age related macular degeneration and retinitis pigmentosa exist, there have been a variety of technologies developed that attempt to restore some visual functioning. These include electronic prostheses (either cortical or retinal), optogenetics, and gene therapy. However to date only retinal prosthetic devices have been implanted in patients, with cortical implant and optogenetic clinical trials still under way.  These technologies are not likely to completely recreate normal vision, but rather provide essential visual input that can improve everyday functioning in patients. 

In normally sighted individuals, visual processing, from the retina to the early cortical areas, includes complementary firing of on and off-center cells in response to visual stimulation. If an on-center cell fires, the corresponding off-center cell is suppressed, and vice versa. While this on-off pattern of reciprocal firing originates in the retina, it remains a fundamental characteristic of the majority of cells throughout the early visual areas. Unless electrodes are small enough to target a single cell, electrical stimulation will result in on- and off-center cells being stimulated simultaneously: leading to unnatural cell population responses in early visual areas. 


In my research, we proposed a unique way to mimic this distortion. While we cannot recreate these distortions in the retina, it is possible to simultaneously stimulate on- and off-center populations in early visual cortex, at the level of V1 layer 4. 


![](https://github.com/resquenazi/neuro545_Project/blob/master/figures/matrix.png)



The image above shows an example of the filtering process. The upper two panels show an example scene (I), and the contrast-reversed version of that scene (Icr). The leftmost panels show the two radial checkerboard Fourier filters we will use: F and F′. Filters are shown in the Fourier domain, with spatial frequency increasing with distance from the center of the image and orientation changing along the polar angle dimension. Each filter is a complement of the other, so the full spatial frequency and orientation content of the scenes is divided equally across the two filters. Here, we are assuming that original and contrast-reversed images can act as a proxy for inappropriate on and off-cell population responses, because regions that would, in the original image, produce strong on-responses will produce strong off-responses in the contrast reversed image, and vice versa.

The original (I) and the contrast-reversed scene (Icr= 1-I) were each converted into the Fourier domain, multiplied by one of the two Fourier filters, and then converted back to image space using the inverse Fourier transform. The middle panels show the 4 examples of possible filtering: I * F, I * F′, Icr*F, and Icr* F′ (where * denotes convolution). Note that the sum, [I * F] + [I * F′], equals the original image I. The sum, [Icr * F] + [Icr * F′], equals the original contrast reversed image Icr. In one eye, we will present the sum of two filtered images, [I * F′] + [Icr * F], such that half the spatial frequency and orientation content is based on the original image and the other half is based on the contrast reversed image. In the other eye we will present the sum of the complementary filtered images, [I * F] + [Icr * F′]. The sum of the distorted images (rightmost panels, Figure 1) in each eye results in a blank image, indicating that all the spatial frequency and orientation information of both the original and contrast reversed image is preserved, but reorganized to create visual input that is impossible in the real world

## The Project

### Matlab portion

For this portion of the project, I decided to explore the differences between the two images using different versions of the filter. As you can see above, the radial checkerboard filter has hard cutoffs in the transition between parts of the image that are passed through, and parts of the image that are not. As we discussed a bit in class during the fourier lesson, this can have harsh effects on the resulting image when filtering. This is because just a simple edge is actually made up of several different sine waves to create that perfect line. This has consequences when filtering images using these edges because it creates ringing in the processed image. However, because of the nature of the project described above, it is imperative that the two images do not share 'on' and 'off' information. 

I decided to explore what the differences in the images are between using this filter with hard edges, to the same filter with blurred edges. It turns out that the resulting images using a filter with smooth edges results in a pair of images that is much smoother. However, as is discovered in the matlab code, the resulting fourier spectrums of these images indicate that a lot of information is shared between the two eyes. 

If you go to the 'Matlab Scripts' folder, you can run through the code I created in the 'Project.m' file. 

### Model portion

For this portion of the project, I worked with Wyeth Bair's binocular energy model. 


![](https://github.com/resquenazi/neuro545_Project/blob/master/figures/model_description.png)

This model was chosen mainly because of the need to present 180 degree phase offset stimuli to two eyes rather than 1 in the way that was done in the stimulus description above. When presenting two sine wave gratings that are contrast reversed in each eye, the model indicates that the firing probability for binocular neurons that prefer zero disparity is minimal. 

This indicates that these binocular neurons are compromised for contrast reversed stimuli. However, this model is only useful for showing the response of neurons that prefer a zero disparity stimulus. As we know, there are other neurons that do not have this same preference, so these neurons would likely display normal firing patterns for contrast reversed stimuli in each eye.

to see how the model was created, and to run/interact with it, go to the 'Model' section of this repository. 
