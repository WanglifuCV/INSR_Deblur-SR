function  Dict   =  KMeans_PCA( im, par, cls_num )
b         =   par.win;
X   =   Get_patches(im, b, 1, 1);
for i=2:2:6
   Xs =   Get_patches(im, b, 1, 0.8^i);
   X          =   [X Xs];
end
itn         =  14;     
rand('seed',0);

[cls_idx, vec, cls_num]   =  Clustering(X, cls_num, itn);
vec   =  vec';
[s_idx, seg]   =  Proc_cls_idx( cls_idx );
PCA_D          =  zeros(b^4, cls_num);
for  i  =  1 : length(seg)-1
    idx    =   s_idx(seg(i)+1:seg(i+1));    
    cls    =   cls_idx(idx(1));    
    X1     =   X(:, idx);
    P      =  getpca(X1, par.nSig); 
    PCA_D(:,cls)    =  P(:);
end


X      =   Get_patches(im, b,par.step, 1);
cls_idx     =   zeros(size(X, 2), 1);
L           =   size(X,2);
vec         =   vec';

SSIM=dis_SSIM(vec,X');
[val,ind]    =   max(SSIM, [], 2);
 for k=1:L
    cls_idx(k)   =   ind(k);
 end
[s_idx, seg]    =   Proc_cls_idx( cls_idx );
Dict.PCA_D      =   PCA_D;
Dict.cls_idx    =   cls_idx;
Dict.s_idx      =   s_idx;
Dict.seg        =   seg;



function  Px =  Get_patches( im, b, s, scale )
im        =  imresize( im, scale, 'bilinear' );
[h, w, ch]   =  size(im);
if  ch==3
    lrim      =  rgb2ycbcr( uint8(im) );
    im        =  double( lrim(:,:,1));    
end

N         =  h-b+1;
M         =  w-b+1;
r         =  [1:s:N];
r         =  [r r(end)+1:N];
c         =  [1:s:M];
c         =  [c c(end)+1:M];
L         =  length(r)*length(c);
Px        =  zeros(b*b, L, 'single');

k    =  0;
for i  = 1:b
    for j  = 1:b
        k       =  k+1;
        blk     =  im(r-1+i,c-1+j);
        Px(k,:) =  blk(:)';        
    end
end


function   [cls_idx,vec,cls_num]  =  Clustering(X, cls_num, itn)
X         =   X';
[L,b2]    =   size(X);
P         =   randperm(L);
P2        =   P(1:cls_num);
vec       =   X(P2(1:end), :);
m_num     =   250;    

for i = 1 : itn
    cnt       =  zeros(1, cls_num);    
    
    v_dis    =   zeros(L, cls_num);

[pc1,score1,latent1] = pca(X);
[pc2,score2,latent2] = pca(vec);
tem=find(latent1<100);
dimen=tem(1);
pca_vec=score2(:,1:dimen);
pca_X=score1(:,1:dimen);

for  k = 1 : cls_num
        v_dis(:, k) = (pca_X(:,1) - pca_vec(k,1)).^2;
        for c = 2:dimen
            v_dis(:,k) =  v_dis(:,k) + (pca_X(:,c) - pca_vec(k,c)).^2;
        end
end 
[val, cls_idx]     =   min(v_dis, [], 2); 

    [s_idx, seg]   =  Proc_cls_idx( cls_idx );
    for  k  =  1 : length(seg)-1
        idx    =   s_idx(seg(k)+1:seg(k+1));    
        cls    =   cls_idx(idx(1));
        vec(cls,:)    =   mean(X(idx, :));
        cnt(cls)      =   length(idx);
    end        
    
    if (i==itn-2)
        [val, ind]  =  min( cnt );      
        while (val<m_num) && (cls_num>=40)
            vec(ind, :)    =  [];
            cls_num       =  cls_num - 1;
            cnt(ind)      =  [];
            [val, ind]    =  min(cnt);
        end        
    end
end