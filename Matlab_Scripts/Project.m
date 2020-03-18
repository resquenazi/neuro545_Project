% Neuro 545 Project written by Rebecca Esquenazi

% The following code uses concepts from Fourier transforms covered in
% class. This code was adapted and extended from an existing project I'm
% working on which aims to explore perceptual learning for retinal
% prosthesis users (see more in the README section of this repository). In
% the following sections, a series of fourier filters will be created and
% convolved with an 1) image of a natural scene, and 2) its contrast reversed
% complement. The main goal of this portion of the project is to examine
% (i.e. 'clipping') the filter, or blurring the edges of the filter.
% As discussed in class, fourier transforms of an edge versus a gaussian
% result in different results that are of interest in image processing, as
% well as the current research question.

%% Parameters for Radial Checkerboard Filter
% fpar is a structure that stores all the necessary parameters for the
% filter

fpar.nPix = 1024; %size of the filter 
fpar.spokes =4; % number of spokes (i.e. rings)
fpar.freq = 13; % number of cycles per degree
fpar.wDeg = 5; % size of filter in degrees of visual angle
fpar.clip = 1; % this option binarizes the filter (1), or causes the filter to have a large sigmar that blurs the edges
fpar.plot = 1; % turn plotting on (1) or off (0)
fpar.n = 20; % radial scaling factor
fpar.rect = [1 1 1024 1024]; % for plotting purposes
fpar.radPhase = pi/2; % phase of the radius
fpar.spokePhase = pi/4; % phase of the spokes
filterName = 'RadialCheckerboard'; %name of the filter 

%% Make filter

% This function takes in the name of a filter that you want to create, and
% parameters for that filter. If there are no parameters given, then
% function automatically creates filter parameters in a structure called
% 'fpar'. The function outputs a bank of two filters that are
% contrast reversed complements of each other. These filters are already in
% the fourier domain and are ready to be convolved with images later on.
% There are two options - either create two radial checkerboard filters, or
% two checkerboard filters.

% First, we will start with a radial checkerboard filter (as indicated by
% the string "filterName").

[filterbankClip, fpar] = MakeFilter(filterName, fpar);

%% Define weights to convert filtered inputs into eye inputs

% the idea behind these stimuli is that the resulting images will have half
% the spatial frequency and orientation information 'flipped' or 'contrast
% reversed' in each eye. In order to do that, weights have been assigned to
% each eye that rearrange the filters in the filter bank.

% F = Original Filter 
% F' = Contrast Reversed Filter
% I = Original Image
% I' = Contrast Reversed Image

eyeWeight(1).w = zeros(length(filterbankClip), 2); % left eye
eyeWeight(1).w(1, 1) = 1; eyeWeight(1).w(2, 2) = 1; % F * I + F' * I'
eyeWeight(2).w = zeros(length(filterbankClip), 2); % right eye
eyeWeight(2).w(2, 1) = 1; eyeWeight(2).w(1, 2) = 1; % F * I' + F' * I

%% Fourier the image and it's inverse

% This is where we get to actually apply the filters to each image using
% convolution. First, the image is read in and converted into proper
% formatting for the function 'FilterScene'.

Img = imread('Scene10_EXSEMSYN.png');
Img = mean(Img,3); % take average of R, G, B matrices
Img = imscale(Img./255);% scale the image between 0 and 1
Img = imresize(Img,max(size(Img))/fpar.nPix); % makes the long edge of image fill fpar.nPix
Img = insertImg(ones(fpar.nPix)*mean(Img(:)), Img); 

% Running this function will output 4 figures. 
[filtImgClip, eyeImgClip] = FilterScene(Img, filterbankClip, fpar, eyeWeight);

% EXPLANATION OF FIGURES BELOW: 

% 2) This is a plot of the original image (I), and contrast and reversed
% image (I'). These two images are derived from the file
% 'Scene10_EXSEMSYN.png'. These are the two images that will be convolved
% with the two filters F and F'.

% 3) This is a plot of a matrix of possible image combinations (see more in
% the README section of this repository). Having two images, (I and I') and
% 2 filters (F and F') make it possible to have 4 image convolution combinations. The
% following combinations are possible: 

% 1. I * F - the original image (I) convolved with the original filter (F)

% 2. I * F' - the original image (I) convolved with the contrast reversed
%           filter (F')

% 3. I' * F - the contrast reversed image (I') convolved with the original
%           filter (F)

% 4. I' * F' - the contrast reversed image (I') convolved with the contrast
%            reversed filter (F)

% 4) This plot shows the resulting left and right eye images that would be
% presented to a subject. The left eye consists of: [I * F?] + [I’* F] and
% the right eye consists of [I * F] + [I’ * F?]. This rearranges the
% spatial frequency and orientation information so that it is half reversed
% in each eye, creating input that is impossible in the real world. This is
% meant to mimic the distortion caused by electronic prosthesis in which
% on/off cells are simultaneously stimulated. However, in these stimuli,
% simultaneous stimulation of on and off cells are mimiced cortically
% (rather than retinally) because regions producing on-responses in one eye
% produce off-responses in the other eye. 

% 5) Figure 5 is a sanity check. This is the result of summing [I * F + I *
% F'] and checking that we get the original image back. Getting the
% original image back tells us that we are not losing any information in
% the image manipulation. 

%% Binarization (i.e. clipping)

% fpar.clip = 1 is what the original fpar.clip parameter is set to. Notice in figure
% 1 when the filter is created that the edges are hard cutoffs rather than
% blurred. This filter is crucial for my research because there is no
% shared spatial frequency or orientation information information across
% the two eyes. However, if we don't clip the filter, it looks different.

fpar.clip = 0; 
figure (1) 
clf

[filterbankSmooth, fpar] = MakeFilter(filterName, fpar);

% now, what we have are regions of the image (i.e. the blurred parts) that
% actually permit the sharing of spatial frequency and orientation
% information across the two eyes. Beucase I am trying to cortically
% stimulate on and off-cells simultaneously, this filter does not allow me
% to achieve that goal. However, a gaussian ramp up/down of an edge is
% what's often used in image processing to decrease the amount of ringing
% in a filtered image. For example, using this filter causes the resulting
% left and right eye images to be different in each eye:

close all

[filtImgSmooth, eyeImgSmooth] = FilterScene(Img, filterbankSmooth, fpar, eyeWeight);

% although figure 5 ensures that we are still preserving all the
% information in images I and I', what we see in the resulting image
% combination is a much clearer representation of the scene. This is
% because there is less ringing artifacts in the scenes since there is
% shared information across the two eyes. 

% as a further sanity check, the phases of each image can be plotted next
% to each other: 

[X,Y] = meshgrid(linspace(-fpar.wDeg/2,fpar.wDeg/2,fpar.nPix+1)); %generating matrix that goes from -wDeg/2 (-2.5) to wdeg/2 (2.5) in steps of the size of "testImg" (n=1080)
X = X(1:end-1,1:end-1); %putting 0 in the center
Y = Y(1:end-1,1:end-1); % putting 0 in the center

% put the final images in the fourier domain
Img_Fourier = complex2real2(fft2(eyeImgClip(1).Img, fpar.nPix,fpar.nPix), X,Y);
ImgCR_Fourier = complex2real2(fft2(eyeImgClip(2).Img,fpar.nPix,fpar.nPix),X,Y);

% plot the first 10 frequencies (low frequencies) along with their phase.
% This demostrates that at any given point, the each image is always 180
% degrees out of phase with each other

subplot(1,2,1)
plot(Img_Fourier.freq(512:522),abs(Img_Fourier.ph(512:522)));
title('Phase of Left Eye Image');
xlabel('Frequency');
ylabel('Phase (in degrees)');
subplot(1,2,2)
plot(ImgCR_Fourier.freq(512:522),abs(ImgCR_Fourier.ph(512:522)));
title('Phase of Right Eye Image');
xlabel('Frequency');
ylabel('Phase (in degrees)');

% this would also be the same for the smooth image, but for brevity, I'm
% only showing the phase behavior of the clipped images. 

%% Exploring the differences between clipped and smooth images

% the code above takes the fft of the image and creates a nice structure of
% bits of information about the image such as the phases and amplitudes.
% However, we can't plot the magnitude spectrum since all of the
% information has been converted from complex to real numbers. Therefore I
% worked with the standard fft function to create magnitude spectrums of
% the smooth and clipped resulting images.

close all
figure
subplot(2,3,1); 
imagesc(filterbankClip(1).filt);
title('Clipped Filter'); axis square;

subplot(2,3,2)
imagesc(eyeImgClip(1).Img); colormap gray; axis square;
title('Clip-filtered Image');

% fft of clipped left eye image
smoothFiltFFT = fft2(double(eyeImgClip(1).Img));
clipFiltFFT_shift = fftshift(smoothFiltFFT);
imFFT2_clip = log2(clipFiltFFT_shift);
imFFT3_clip = abs(imFFT2_clip);
subplot(2,3,3)
showImage(imFFT3_clip); title('Clip-filtered image power spectrum');

subplot(2,3,4)
imagesc(filterbankSmooth(1).filt);
title('Smooth Filter'); axis square;

subplot(2,3,5)
imagesc(eyeImgSmooth(1).Img); colormap gray; axis square;
title('Smooth-filtered Image'); axis square;

% fft of smooth left eye image
imFFT_smooth = fft2(double(eyeImgSmooth(1).Img));
imFFT1_smooth = fftshift(imFFT_smooth);
imFFT2_smooth = log2(imFFT1_smooth);
imFFT3_smooth = abs(imFFT2_smooth);
subplot(2,3,6)
showImage(imFFT3_smooth); title('Smooth-filtered image power spectrum');

% The fourier spectra of the resulting images are really interesting to me.
% I actually thought the spectra would be reversed, and that the image that
% went through the clipped filter would look much more ringy. The
% clip-filtered image fourier spectrum indicates that there is a lot of
% energy at the cardinal orientations, with a bit of energy at some other
% orientations. Then I realized that this is because the filter restricted
% more information coming from multiple other orientations, while the
% smooth-filtered image contains much more information since the blurred
% edges let more spatial frequency and orientation information into the
% image.
