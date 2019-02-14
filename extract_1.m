%%%%%%%%%%% Extracting Watermark %%%%%%%%
addpath 'imgs'
filename = '1';
format = '.bmp';
%%%
I=imread([filename '_embed' format]);

M=64;N=360;Lc=64;
Np=72;
key=0;
a=power(2,1/M);
Tm=zeros(M,N);

im = rgb2ycbcr(I);
J1 = im(:,:,1); %Ƕ��Yͨ��
J1=double(J1);
F1a=fft2(J1);             %%����Ҷ�任������ֱ������ϵ
F2a=fftshift(F1a);
ampa=abs(F2a);
[dxa,dya]=size(J1);
xa=floor(dxa/2+1);
ya=floor(dya/2+1);
da=min(dxa,dya);

fn=0.25;             %%%%ѡȡǶ����Ϣ��λ�ã�Ϊ��Բ��
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
for alpha=-pi:pi/2000:0           %����Ҷ�׵İ�Բ
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
figure;
imshow(amp1a);


[x2,y2]=find(amp1a==255);     %%Ƕ��λ�õĸ���Ҷֱ������ϵ����
A1a=[x2,y2];
A11a=A1a(:,1)-xa;
A12a=A1a(:,2)-ya;
Aa=[A11a,A12a];
[sizeAa,~]=size(Aa);

amp2a=zeros(dxa,dya);        %%����ULPM�任
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
figure;
imshow(amp2a);      %%%����������ϵ�е�����
B1a=unique(Ba,'rows');

%%%%%%%%%%%%%%%%%%%��ˮӡͬ���õ�ͬ����ֵ����%%%%%%%%%%%%%%%%
amp3=zeros(Ma,Na);
for i=1:sizeAa
    if Ba(i,1)<=0 || Ba(i,2)<=0
        continue
    end
    amp3(Ba(i,1),Ba(i,2))=ampa(A1a(i,1),A1a(i,2));
end

% phase3=zeros(Ma,Na);
% [sizeB1a,~]=size(B1a);
% for i=1:sizeB1a
%     if B1a(i,1)<=0 || B1a(i,2)<=0
%         continue
%     end
%     aa=0;
%     C=[];
%     for j=1:sizeAa
%         if B1a(i,1)==Ba(j,1) && B1a(i,2)==Ba(j,2)
%             aa=aa+1;
%             C=[C;j];
%         end
%     end
%     sum1=0;
%     for k1=1:aa
%         sum1=sum1+phasea(A1a(C(k1,1),1),A1a(C(k1,1),2));
%     end
%     mean1=sum1/aa;
%     phase3(B1a(i,1),B1a(i,2))=mean1;
% end

rk=zeros(Ma-M+1,Na-N+1);      %%%��ģ�������������ȷ��ˮӡλ��
for i=1:Ma-M+1
    for j=1:Na-N+1
    amp4a=amp3(i:i+M-1,j:j+N-1);
    RK=corr2(amp4a,Tm);
    rk(i,j)=RK;
    end
end

[rk1max,rk2max]=find(rk==max(max(rk)));
amp4=amp3(rk1max:rk1max+M-1,rk2max:rk2max+N-1);

%%%%%%%%%%%%%%%%%%%%��ͬ��������н���Ƶ����%%%%%%%%%%%%%%%%
key=0;                  %��Կ
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
for i=1:Lc                    %%��pa����أ����ϵ������0Ϊ1������Ϊ0
    amp6=amp5c(i,:);
    xg=corrcoef(amp6,pa);
    XG=xg(1,2);
    if XG>0
        ma(1,i)=1;
    else
        ma(1,i)=0;
    end                          %%ma��Ϊ��ȡ����ˮӡ��Ϣ
end

error1=0;                   %%��ʾ������
for i=1:Lc
    if m(1,i)~=ma(1,i)
        error1=error1+1;
    end
end
disp(['error=',num2str(error1)]);