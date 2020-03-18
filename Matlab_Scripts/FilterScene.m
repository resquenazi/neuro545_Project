function [filtImg, eyeImg] = FilterScene(Img, filterbank, fpar, eyeWeight)
% This function requires the following functions
%       - complex2real2
%       - real2complex2
%       - myifft2

%% Fourier the image and it's inverse

[X,Y] = meshgrid(linspace(-fpar.wDeg/2,fpar.wDeg/2,fpar.nPix+1)); %generating matrix that goes from -wDeg/2 (-2.5) to wdeg/2 (2.5) in steps of the size of "testImg" (n=1080)
X = X(1:end-1,1:end-1); %putting 0 in the center
Y = Y(1:end-1,1:end-1);
Img_Fourier = complex2real2(fft2(Img, fpar.nPix,fpar.nPix), X,Y);
ImgCR_Fourier = complex2real2(fft2(1-Img,fpar.nPix,fpar.nPix),X,Y);

%% original image, original Fourier filter

for f=1:length(filterbank)
    tmp = Img_Fourier;
    tmp.amp = [filterbank(f).filt.*tmp.amp]; % multiply the filter times the amplitudes of the image in the fourier domain
    filtImg(f, 1).img = myifft2(tmp); % apply inverse FFT to take image out of fourier domain
    
    tmp = ImgCR_Fourier;
    tmp.amp = [filterbank(f).filt.*tmp.amp]; % same for the contrast reversed image
    filtImg(f, 2).img = myifft2(tmp); % apply inverse FFT to take image out of fourier domain
    
end

for e = 1:2 % for each eye
    eyeImg(e).Img = zeros(size(filtImg(f, 1).img));
    for f=1:length(filterbank)
        for c = 1:2
            eyeImg(e).Img = eyeImg(e).Img + [filtImg(f,c).img .*eyeWeight(e).w(f, c)];
        end
    end
end

% Plotting
if fpar.plot
    
    cstr={' Image I', ' Image I'''}; % I = original image, I' = Contrast reversed image
    estr={'Left Eye', 'Right Eye'};
    figure(2); set(gcf, 'Name', 'Original & contrast rev. image')

    % Figure 2: Plot original contrast and reversed image
    for c=1:2 
        subplot(1, 2, c);
        if c==1
            img = Img(fpar.rect(2):fpar.rect(4), fpar.rect(1):fpar.rect(3)); % scaling image to be 1024x1024
        else
            img = 1-Img(fpar.rect(2):fpar.rect(4), fpar.rect(1):fpar.rect(3)); % doing the same for other eye
        end    
    [min(img(:)), max(img(:))];
    image(uint8(255.*img)); axis equal; axis off;  colormap(gray(256));
    title(cstr{c});
    sgtitle('Input Images');
    end

    % Figure 3: plot matrix of [I * F], [I * F'], [I' * F] and [I' * F']
    figure(3)
    for f=1:length(filterbank)
        for c=1:2 % original contrast and reversed
            subplot(length(filterbank),2,[[(f-1)*length(filterbank)]+c]);
            img = filtImg(f, c).img(fpar.rect(2):fpar.rect(4), fpar.rect(1):fpar.rect(3));
            image(uint8(255.*img)); axis equal; axis off;  colormap(gray(256));
            title(['Filter ', num2str(f), ' *', cstr{c}]);
            sgtitle('Convolution of Filters with Images');
        end
    end
    
    % Figure 4: Plot resulting Left and Right Eye Images
    figure(4); set(gcf, 'Name', 'Eye Input')
    for e = 1:2
        subplot(1, 2, e);
        img = eyeImg(e).Img(fpar.rect(2):fpar.rect(4), fpar.rect(1):fpar.rect(3));
        imagesc(img./sum(eyeWeight(e).w(:))); axis equal; axis off;  colormap(gray(256));
        title(estr{e})
        sgtitle('Resulting Images');
    end
    
    % Figure 5: Sanity Check. Make sure we get the original images back by
    % summing up the columns of figure 3 to catch whether all information
    % is preserved during the process of reversing half the spatial
    % frequency and orientaiton information in each eye. 
    figure(5); set(gcf, 'Name', 'Sanity Check: Sum of Figure 3 Columns')
    for i = 1:2     
        tmp = [2, 1];
        subplot(1, 2, i);
        img = [filtImg(i,i).img + filtImg(tmp(i),i).img] ./ 2.0;
        imagesc(uint8(255.*img)); axis equal; axis off;  colormap(gray(256));
        title(cstr{i})
        sgtitle('Sum of Figure 3 Columns');
    end
end
