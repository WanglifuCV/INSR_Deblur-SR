function  X_m   =  Update_NLM( im, opts, blk_arr, wei_arr )
[h, w, ch] =   size(im);
b          =   opts.win;
b2         =   b*b*ch;
k          =   0;
s          =   opts.step;
N     =  h-b+1;
M     =  w-b+1;
L     =  N*M;
r     =  1:s:N;
r     =  [r r(end)+1:N];
c     =  1:s:M;
c     =  [c c(end)+1:M];
X     =  zeros(b*b,L,'single');
for i  = 1:b
    for j  = 1:b
        k    =  k+1;
        blk  =  im(i:h-b+i,j:w-b+j);
        blk  =  blk(:);
        X(k,:) =  blk';            
    end
end
X_m     =   zeros(length(r)*length(c),b*b,'single');
X0      =   X';
for i = 1:opts.nblk
   v        =  wei_arr(:,i);
   X_m      =  X_m(:,:) + X0(blk_arr(:,i),:) .*v(:, ones(1,b2));
end
X_m    =   X_m';
return;

