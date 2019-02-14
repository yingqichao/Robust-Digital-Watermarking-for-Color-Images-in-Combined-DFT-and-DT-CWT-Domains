clear all;
addpath 'imgs'
errors_rotate = zeros(5,1);ssim = 0;
errors_quanti = zeros(5,1);psps = 0;
format = '.bmp';
for filename=1:120
    %%%%%%%%%%%%%%构建水印矩阵%%%%%%%%%%%%%
        %%%
        disp(['File: ' num2str(filename)]);
        I=imread([num2str(filename)  format]);
        I = imresize(I,[512,512]);
        %I=rgb2gray(I);
        im = rgb2ycbcr(I);
        I1 = im(:,:,1); %嵌入Y通道
        I1=double(I1);
%         figure,imshow(I1);
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
        alpha=0.35;             %%水印嵌入强度
        for i=1:sizeA
            if B(i,1)<=0 || B(i,2)<=0 || B(i,1)>M || B(i,2)>N
                continue
            end
            WT4=WT3(B(i,1),B(i,2)); 
            FF2(A(i,1)+x,A(i,2)+y)=(1+alpha*WT4)*F2(A(i,1)+x,A(i,2)+y);%%%傅里叶对称性保证反傅里叶变换后为实数
            FF2(2*x-(A(i,1)+x),2*y-(A(i,2)+y))=(1+alpha*WT4)*F2 (2*x-(A(i,1)+x),2*y-(A(i,2)+y));
        end
        II1=ifft2(fftshift(FF2));

%         figure;
%         imshow(II1);
%         savename = [filename '_embed' format];
%         [PSNR,~] = psnr(im(:,:,1),II1);
%       
%         psps= ((psps*(filename-1))+PSNR)/filename;

        im(:,:,1) = II1;
        nim = ycbcr2rgb(im);
%         [mssim, ~] = fxssim(I,nim);
%         ssim = ((ssim*(filename-1))+mssim)/filename;
%                 disp(['SSIM: ' num2str( ssim) ' PSNR: ' num2str( psps)]);
%         figure;
%         imshow(nim);
%         imwrite(nim,savename);

        for counter=1:10
            if counter<=5
%                 %JPEG
%                     quality=10+counter*10;
%                     imwrite(nim,'temp.jpeg','Quality',quality);
%                     J = imread('temp.jpeg');
%                     disp(['JPEG level: ' num2str(quality)]);
                    %crop
                    J=nim(ceil((counter)/40*512):ceil((1-(counter)/40)*512),ceil((counter)/40*512):ceil((1-(counter)/40)*512),:);
            else
%                 %Rotate
%                 	angle=(counter-9)*10;
%                     J = imrotate(nim,angle);
%                     disp(['Angle Rotation: ' num2str(angle)]);
                    %scale
                    J=imresize(nim,1-(counter-5)/20);
            end
            im = rgb2ycbcr(J);
            J1 = im(:,:,1); %嵌入Y通道
            J1=double(J1);
            F1a=fft2(J1);             %%傅里叶变换，建立直角坐标系
            F2a=fftshift(F1a);
            ampa=abs(F2a);
            [dxa,dya]=size(J1);
            xa=floor(dxa/2+1);
            ya=floor(dya/2+1);
            da=min(dxa,dya);

            fn=0.25;             %%%%选取嵌入信息的位置，为半圆弧
            x0=xa;
            y0=ya;
            Ra=0.24*(da-1);
            Ma=2*M;   
            if mod(floor(N*4/3),2)==1
                Na=floor(N*4/3)+1;
            else
                Na=floor(N*4/3);
            end
            r1maxa=Ra*power(a,Ma/2);
            r1mina=Ra*power(a,-Ma/2);
            amp1a=zeros(dxa,dya);
            for alpha=-pi:pi/2000:0           %傅里叶谱的半圆
                for x1=ceil(x0+r1maxa*sin(alpha)):x0
                    for y1=ceil(y0+r1maxa*cos(alpha)):y0
                        amp1a(x1,y1)=255;
                    end
                end
                for x1=ceil(x0+r1maxa*sin(alpha)):x0
                    for y1=y0:ceil(y0+r1maxa*cos(alpha))
                        amp1a(x1,y1)=255;
                    end
                end
            end
            for alpha=-pi:pi/2000:0           
                for x1=ceil(x0+r1mina*sin(alpha)):x0
                    for y1=ceil(y0+r1mina*cos(alpha)):y0
                        amp1a(x1,y1)=0;
                    end
                end
                for x1=ceil(x0+r1mina*sin(alpha)):x0
                    for y1=y0:ceil(y0+r1mina*cos(alpha))
                        amp1a(x1,y1)=0;
                    end
                end
            end
          

            [x2,y2]=find(amp1a==255);     %%嵌入位置的傅里叶直角坐标系坐标
            A1a=[x2,y2];
            A11a=A1a(:,1)-xa;
            A12a=A1a(:,2)-ya;
            Aa=[A11a,A12a];
            [sizeAa,~]=size(Aa);

            amp2a=zeros(dxa,dya);        %%进行ULPM变换
            Ba=zeros(sizeAa,2);
            for i=1:sizeAa
                x3=Aa(i,1);
                y3=Aa(i,2);
                r=sqrt(x3^2+y3^2);
                theta=atan(y3/x3);
                l1a=floor(log(r/Ra)/log(a))+Ma/2;
                l2a=floor(N*theta/pi);
                amp2a(l1a+xa,l2a+ya)=255;
                Ba(i,1)=l1a;
                Ba(i,2)=l2a+Na/2;
            end

            B1a=unique(Ba,'rows');

            %%%%%%%%%%%%%%%%%%%做水印同步得到同步幅值矩阵%%%%%%%%%%%%%%%%
            amp3=zeros(Ma,Na);
            for i=1:sizeAa
                if Ba(i,1)<=0 || Ba(i,2)<=0
                    continue
                end
                amp3(Ba(i,1),Ba(i,2))=ampa(A1a(i,1),A1a(i,2));
            end

            rk=zeros(Ma-M+1,Na-N+1);      %%%求模板与待测的相关性确定水印位置
            for i=1:Ma-M+1
                for j=1:Na-N+1
                amp4a=amp3(i:i+M-1,j:j+N-1);
                RK=corr2(amp4a,Tm);
                rk(i,j)=RK;
                end
            end

            [rk1max,rk2max]=find(rk==max(max(rk)));
            amp4=amp3(rk1max:rk1max+M-1,rk2max:rk2max+N-1);

            %%%%%%%%%%%%%%%%%%%%对同步矩阵进行解扩频调制%%%%%%%%%%%%%%%%
            key=0;                  %密钥
            rng('default');
            rng(key);
            pa=zeros(1,Np);
            for i=1:Np
                pa(1,i)=randi([0,1]);
                if pa(1,i)==0
                    pa(1,i)=-1;
                end
            end
            amp5=zeros(M/2,N/2);
            for i=1:M/2
                for j=1:N/2
                    amp5(i,j)=(amp4(2*i-1,2*j-1)+amp4(2*i-1,2*j)+amp4(2*i,2*j-1)+amp4(2*i,2*j))/4;
                end
            end
            amp5a=reshape(amp5,1,M*N/4);
            amp5b1=amp5a(:,NT/3+1:NT/3+Lc*Np/2);
            amp5b2=amp5a(:,NT/3*2+Lc*Np/2+1:NT/3*2+Lc*Np);
            amp5b=[amp5b1,amp5b2];
            amp5c=reshape(amp5b,Lc,Np);
            ma=-2*zeros(1,Lc);
            for i=1:Lc                    %%与pa做相关，相关系数大于0为1，否则为0
                amp6=amp5c(i,:);
                xg=corrcoef(amp6,pa);
                XG=xg(1,2);
                if XG>0
                    ma(1,i)=1;
                else
                    ma(1,i)=0;
                end                          %%ma即为提取到的水印信息
            end

            error1=0;                   %%显示误码数
            for i=1:Lc
                if m(1,i)~=ma(1,i)
                    error1=error1+1;
                end
            end
            if counter<=5
                errors_quanti(counter) = ((errors_quanti(counter)*(filename-1))+error1)/filename;
                disp(['Cropping :' num2str(counter/10) ': '  num2str(errors_quanti(counter))]);
            else
                errors_rotate(counter-5) = ((errors_rotate(counter-5)*(filename-1))+error1)/filename;
                disp(['Scaling : ' num2str((counter-5)/10) ': '  num2str(errors_rotate(counter-5))]);
            end
        end
      
        
        
end
disp('Done.');
save crop.mat;
close all;
