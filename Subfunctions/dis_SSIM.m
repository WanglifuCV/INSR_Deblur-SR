function [ SSIM ] = dis_SSIM( m1,m2)

m1=double(m1);
m2=double(m2);
[r1,c1]=size(m1);
[r2,c2]=size(m2);

u1=mean(m1,2);                            
u2=mean(m2,2);
u1_2=u1.^2;                               
u2_2=u2.^2;
tempo1=bsxfun(@plus,u1_2',u2_2);          
              
m1_u=bsxfun(@minus,m1,u1);                
m2_u=bsxfun(@minus,m2,u2);

sigma1_2=mean(m1_u.^2,2);                 
sigma2_2=mean(m2_u.^2,2);
tempo2=bsxfun(@plus,sigma1_2',sigma2_2);

k(1)=0.01;
k(2)=0.03;
L=255;
C1=(k(1)*L)^2;
C2=(k(2)*L)^2;
u1_u2=u2*u1';
sigma12=(m2_u*m1_u')/c1;
SSIM=((2*u1_u2+C1).*(2*sigma12+C2))./((tempo1+C1).*(tempo2+C2));


