clear all;
addpath 'imgs'
%%%%%%%%%%%%%%构建水印矩阵%%%%%%%%%%%%%
filename = '1';
format = '.bmp';
%%%
I=imread([filename  format]);
I = imresize(I,[512,512]);
%I=rgb2gray(I);
im = rgb2ycbcr(I);
I1 = im(:,:,1); %嵌入Y通道
I1=double(I1);
figure,imshow(I1);
[dx,dy]=size(I1);
d=min(dy,dx);

M=64;
N=360;

Lc=64;
Np=72;
key=0;  %密钥
rng('default');
rng(key);
m=zeros(1,Lc);
for i=1:Lc
    m(1,i)=randi([0,1]);
end

rng('default');
rng(key);
p=zeros(1,Np);
for i=1:Np
    p(1,i)=randi([0,1]);
    if p(1,i)==0
        p(1,i)=-1;
    end
end
Wi=zeros(Lc,Np);   %%生成水印的信息部分W
for i=1:Lc
    for j=1:Np
        if m(1,i)==1
            Wi(i,j)=p(1,j);
        elseif m(1,i)==0
            Wi(i,j)=(-1)*p(1,j);
        end
    end
end
W=reshape(Wi,1,Lc*Np);
W1=W(1,1:Lc*Np/2);
W2=W(1,Lc*Np/2+1:Lc*Np);

NT=M*N/4-Lc*Np;   %%生成模板T
rng('default');
rng(key);
T=zeros(1,NT);
for i=1:NT
    T(1,i)=randi([0,1]);
    if T(1,i)==0
        T(1,i)=-1;
    end
end
T1=T(1,1:NT/3);
T2=T(1,NT/3+1:NT*2/3);
T3=T(1,NT*2/3+1:NT);

WT1=[T1,W1,T2,W2,T3];             %%生成水印矩阵
WT=reshape(WT1,M/2,N/2);
WT2=zeros(M,N/2);      %%上采样
for i=1:M/2
    for j=1:N/2
        WT2(2*i-1,j)=WT(i,j);
        WT2(2*i,j)=WT(i,j);
    end
end
WT3=zeros(M,N);
for i=1:M
    for j=1:N/2
        WT3(i,2*j-1)=WT2(i,j);
        WT3(i,2*j)=WT2(i,j);
    end
end

Ws1=zeros(1,Lc*Np/2);          %%生成水印模板
Ws2=zeros(1,Lc*Np/2);
WTm=[T1,Ws1,T2,Ws2,T3];
Tm1=reshape(WTm,M/2,N/2);
WT5=zeros(M,N/2);      %%上采样
for i=1:M/2
    for j=1:N/2
        WT5(2*i-1,j)=Tm1(i,j);
        WT5(2*i,j)=Tm1(i,j);
    end
end
Tm=zeros(M,N);
for i=1:M
    for j=1:N/2
        Tm(i,2*j-1)=WT5(i,j);
        Tm(i,2*j)=WT5(i,j);
    end
end
%%%%%%%%%%%图像的Fourier变换与LPM映射%%%%%%%%%

a=power(2,1/M);

F1=fft2(I1);             %%傅里叶变换，建立直角坐标系
F2=fftshift(F1);
amp=abs(F2);
x=floor(dx/2)+1;
y=floor(dy/2)+1;
%imshow(log(amp+1),[]);

x0=x;     %%%%选取嵌入信息的位置，为半圆弧
y0=y;
R=0.25*(d-1);
r1max=R*power(a,M/2);
r1min=R*power(a,-M/2);
amp1=zeros(dx,dy);
for alpha=-pi:pi/2000:0           %傅里叶谱的半圆
    for x1=floor(x0+r1max*sin(alpha)):x0
        for y1=floor(y0+r1max*cos(alpha)):y0
            amp1(x1,y1)=255;
        end
    end
    for x1=floor(x0+r1max*sin(alpha)):x0
        for y1=y0:floor(y0+r1max*cos(alpha))
            amp1(x1,y1)=255;
        end
    end
end
for alpha=-pi:pi/2000:0           
    for x1=floor(x0+r1min*sin(alpha)):x0
        for y1=floor(y0+r1min*cos(alpha)):y0
            amp1(x1,y1)=0;
        end
    end
    for x1=floor(x0+r1min*sin(alpha)):x0
        for y1=y0:floor(y0+r1min*cos(alpha))
            amp1(x1,y1)=0;
        end
    end
end
% figure;
% imshow(amp1);


[x2,y2]=find(amp1==255);     %%嵌入位置的傅里叶直角坐标系坐标
A1=[x2,y2];

A11=A1(:,1)-x;
A12=A1(:,2)-y;
A=[A11,A12];
[sizeA,~]=size(A);

amp2=zeros(dx,dy);        %%进行ULPM变换
B=zeros(sizeA,2);
for i=1:sizeA
    x3=A(i,1);
    y3=A(i,2);
    r=sqrt(x3^2+y3^2);
    theta=atan(y3/x3);
    l1=floor(log(r/R)/log(a))+M/2;
    l2=floor(N*theta/pi);
    amp2(l1+x,l2+y)=255;
    B(i,1)=l1;
    B(i,2)=l2+N/2;
end
% figure;
% imshow(amp2);   

for i=1:sizeA
    if B(i,1)<=0
        B(i,1)=1;
    elseif B(i,1)>M
        B(i,1)=M;
    end
    if B(i,2)<=0
        B(i,2)=1;
    elseif B(i,2)>N
        B(i,2)=N;
    end
end

B1=unique(B,'rows');
% B2=B;
% for i=1:sizeA
%     if B2(i,1)<=0 || B2(i,2)<=0
%         B2(i,1)=NaN;
%     end
% end
% [m1,n1]=find(isnan(B2));
% B2(m1,:)=[];
% [sizeB,~]=size(B2);

FF2=F2;                   %%对应的水印信息嵌入
alpha=0.2;             %%水印嵌入强度
for i=1:sizeA
    if B(i,1)<=0 || B(i,2)<=0 || B(i,1)>M || B(i,2)>N
        continue
    end
    WT4=WT3(B(i,1),B(i,2)); 
    FF2(A(i,1)+x,A(i,2)+y)=(1+alpha*WT4)*F2(A(i,1)+x,A(i,2)+y);%%%傅里叶对称性保证反傅里叶变换后为实数
    FF2(2*x-(A(i,1)+x),2*y-(A(i,2)+y))=(1+alpha*WT4)*F2 (2*x-(A(i,1)+x),2*y-(A(i,2)+y));
end
II1=ifft2(fftshift(FF2));

figure;
imshow(II1);
savename = [filename '_embed' format];
im(:,:,1) = II1;
nim = ycbcr2rgb(im);
figure;
imshow(nim);
imwrite(nim,savename);

%JPEG
for quality=20:10:100
    imwrite(nim,[filename '_' num2str(quality) '.jpeg'],'Quality',quality);
    disp(['Writing ' filename '_' num2str(quality) '.jpeg']);
end

%Rotate
for angle=10:10:90
    nim1 = imrotate(nim,angle,'crop');
    imwrite(nim1,[filename '_Angle' num2str(angle) format]);
    disp(['Writing ' filename '_Angle' num2str(angle) format]);
end
disp('Done.');
close all;
%J1=II1;           %%没有攻击

% rate=0.25;             %%剪切
% x1=512*rate/2;
% len1=512-512*rate;
% J1=imcrop(II1,[x1 x1 len1 len1]);  
% %J1=imresize(J11,[512 512]);
% figure;
% imshow(J1);

% theta=90;
% J1=imrotate(II1,theta);             %%旋转
% figure;
% imshow(J1);

% theta=0;
% J11=imrotate(II1,theta);             %%旋转剪切
% [size11,~]=size(J11);
% xx=512*tan(theta*pi/180)/(1+tan(theta*pi/180));
% yy=floor(xx/sin(theta*pi/180));
% aa=floor((size11-yy)/2);
% J1=J11(aa:aa+yy,aa:aa+yy);
% %J1=imresize(J1,[512 512]);
% figure;
% imshow(J1);

% sigma=0.71;           %%缩放 
% J1=imresize(II1,sigma);
% figure;
% imshow(J1);

% J1=II1;          %%jpeg压缩
% imwrite(J1,'J1.jpg','quality',90);

% J1=imnoise(II1,'gaussian',0.005);  %%加噪
% figure;
% imshow(J1);