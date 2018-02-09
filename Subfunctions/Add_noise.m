function im_o   =  Add_noise( im, v )

seed  =  0;
randn( 'state', seed );
noise     =  randn( size(im) );
im_o      =  double(im) + v*noise;