
clc;
clear;
im_name  = 'Raccoon.tif';
fn                 =    fullfile('Data\SR_test_images',im_name);     
psf                =    fspecial('gauss', 7, 1.6);
scale              =    3;
nSig               =    5;                                     % 0;

par                =    INSR_SR_Par( nSig, scale, psf );
par.I              =    double( imread( fn ) );
LR                 =    Blur('fwd', par.I, psf);               % The blur operator
LR                 =    LR(1:par.scale:end,1:par.scale:end,:); % Dowmsampling  
par.LR             =    Add_noise(LR, nSig);                   % Add noise 
par.B              =    Set_blur_matrix( par );                
    
[im, PSNR, SSIM, FSIM]     =    INSR_Superresolution( par );     
fprintf('%s: PSNR = %3.2f   SSIM = %f   FSIM = %f\n\n', fn, PSNR, SSIM, FSIM);

if nSig == 0       
   imwrite(im./255, fullfile('Results\SR_results\Noiseless',im_name));
 else
    imwrite(im./255, fullfile('Results\SR_results\Noisy',im_name));
end




