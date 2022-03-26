function smoothedBetas = smooth_betas(betas,Fs,fwhm,whichdim)

% smoothedBetas = smooth_betas(betas,Fs,fwhm,whichdim)
%
% function to smooth n-dimensional matrix of betas (or, indeed, any timeseries) using gaussian smoothing
% kernel
%
% betas is 1-dimensional vector or n-dimensional matrix, with time along one of the axes
% Fs is sampling rate in Hz (default 100)
% fwhm is full-width half maximum of smoothing kernel in ms (default 30)
% whichdim is the dimension of betas along which time runs (and smoothing is to occur) (default 3rd dimension)
%
% [the default settings are based on what is typically used in Maria's experiment, YMMV for other experiments]
%
% LH 17/03/2022

%% handle input arguments

if nargin<2
    Fs = 100; %sampling rate, Hz
end

if nargin<3
    fwhm = 30; % FWHM of smoothing kernel, ms
end

if ~isvector(betas)
    if nargin<4
        whichdim = 3;
    end
    
    if whichdim>ndims(betas) %sanity check
        error('dimension along which smoothing is to occur does not exist');
    end
end

%% set up smoothing kernel
fwhm_samples = Fs*fwhm/1000; % convert fwhm from ms to samples
kernel = gauss(fwhm_samples/2.355,1,ceil(fwhm_samples*5))'; % why 2.355? see https://en.wikipedia.org/wiki/Full_width_at_half_maximum

%% run convolution
if isvector(betas) %1-D convolution, easy peasy
    smoothedBetas = conv(betas,kernel,'same');
else %N-D matrix, convolution along whichdim
    
    %make whichdim the first dimension
    dimSort = [whichdim 1:(whichdim-1) whichdim+1:ndims(betas)];       
    betasPermuted = permute(betas,dimSort);
    
    %do smoothing on the first dimension using conv2, and reshape back into N-D rather than 2-D
    smoothedBetasPermuted = conv2(kernel,1,betasPermuted(:,:),'same');
    smoothedBetasPermuted = reshape(smoothedBetasPermuted,size(betasPermuted));
    
    %permute back into original dimensionality
    tmp = 1:length(dimSort); dimUnsort(dimSort) = tmp; clear tmp
    smoothedBetas = permute(smoothedBetasPermuted,dimUnsort);
end