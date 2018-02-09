function z = Blur (mode,x,psf)
% The blur operator.

ws   =  size(psf);
t    =  (ws-1)/2;

if mode == 'fwd'
    s  = x;
    se = [s(:,t:-1:1,:), s, s(:,end:-1:end-t+1,:)];
    se = [se(t:-1:1,:,:); se; se(end:-1:end-t+1,:,:)];
    
    if size(x,3)==3
        z(:,:,1) = conv2(se(:,:,1),psf,'valid');
        z(:,:,2) = conv2(se(:,:,2),psf,'valid');
        z(:,:,3) = conv2(se(:,:,3),psf,'valid');
    else
        z = conv2(se,psf,'valid');
    end

elseif mode == 'trn'
    y  = x;
    ye = [y(:,t:-1:1,:), y, y(:,end:-1:end-t+1,:)];
    ye = [ye(t:-1:1,:,:); ye; ye(end:-1:end-t+1,:,:)];
    
    if size(ye,3)==3
        z(:,:,1) = conv2(ye(:,:,1),psf,'valid');
        z(:,:,2) = conv2(ye(:,:,2),psf,'valid');
        z(:,:,3) = conv2(ye(:,:,3),psf,'valid');
    else
        z = conv2(ye,psf,'valid');
    end
end
