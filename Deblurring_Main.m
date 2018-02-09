
clc;
clear;

im_name  = 'Starfish.tif';
       
fn                 =    fullfile('Data\Deblurring_test_images',im_name); 
blur_type          =    1 ;                                            % 1: uniform blur kernel;  2: Gaussian blur kernel; 
if blur_type == 1                                                      % When blur_type = 1, blur_par denotes the kernel size; When blur_type = 2, blur_par denotes the standard variance of Gaussian kernel
    blur_par       =    9;                                             % the default blur kernel size is 9 for uniform blur;
else
    blur_par       =   1.6;                                            % the default standard deviation of Gaussian blur kernel is 1.6
end
nSig                    =    sqrt(2);                                  % The standard variance of the additive Gaussian noise;
par                     =    INSR_Deblurring_Par( nSig, blur_type );
par.I                   =    double( imread( fn ) );
[par.bim, par.fft_h]    =    Generate_blur_image(par.I, blur_type, blur_par, nSig);

[im, PSNR, SSIM, FSIM]  =   INSR_Deblurring( par );

imwrite(im./255, fullfile('Results\Deblurring_results',im_name)); 
fprintf('%s: PSNR = %3.2f   SSIM = %f   FSIM = %f\n\n', fn, PSNR, SSIM, FSIM);

