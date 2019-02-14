clear;
clc;

disp('请选择水印图像：');
[filename, pathname] = uigetfile('*.jpg', '读取图片文件');
pathfile=fullfile(pathname, filename);
markbefore=imread(pathfile); 
disp('请选择载体图像：');
[filename2, pathname2] = uigetfile('*.jpg', '读取图片文件');
pathfile2=fullfile(pathname2, filename2);
image=imread(pathfile2); 

markbefore2=rgb2gray(markbefore);
mark=im2bw(markbefore2);    %使水印图像变为二值图
figure(1);      %打开窗口
subplot(2,3,1);    %该窗口内的图像可以有两行三列
imshow(mark),title('水印图像');   %显示水印图像
marksize=size(mark);   %计算水印图像的长宽
rm=marksize(1);      %rm为水印图像的行数
cm=marksize(2);     %cm为水印图像的列数

I=mark;
alpha=30;     %尺度因子,控制水印添加的强度,决定了频域系数被修改的幅度
k1=randn(1,8);  %产生两个不同的随机序列
k2=randn(1,8);
subplot(2,3,2),imshow(image,[]),title('载体图像'); %[]表示显示时灰度范围为image上的灰度最小值到最大值
yuv=rgb2ycbcr(image);   %将RGB模式的原图变成YUV模式
Y=yuv(:,:,1);    %分别获取三层，该层为灰度层
U=yuv(:,:,2);      %因为人对亮度的敏感度大于对色彩的敏感度，因此水印嵌在色彩层上
V=yuv(:,:,3);
[rm2,cm2]=size(U);   %新建一个和载体图像色彩层大小相同的矩阵
before=blkproc(U,[8 8],'dct2');   %将载体图像的灰度层分为8×8的小块，每一块内做二维DCT变换，结果记入矩阵before

after=before;   %初始化载入水印的结果矩阵
for i=1:rm          %在中频段嵌入水印
    for j=1:cm
        x=(i-1)*8;
        y=(j-1)*8;
        if mark(i,j)==1
            k=k1;
        else
            k=k2;
        end;
        after(x+1,y+8)=before(x+1,y+8)+alpha*k(1);
        after(x+2,y+7)=before(x+2,y+7)+alpha*k(2);
        after(x+3,y+6)=before(x+3,y+6)+alpha*k(3);
        after(x+4,y+5)=before(x+4,y+5)+alpha*k(4);
        after(x+5,y+4)=before(x+5,y+4)+alpha*k(5);
        after(x+6,y+3)=before(x+6,y+3)+alpha*k(6);
        after(x+7,y+2)=before(x+7,y+2)+alpha*k(7);
        after(x+8,y+1)=before(x+8,y+1)+alpha*k(8);
    end;
end;
result=blkproc(after,[8 8],'idct2');    %将经处理的图像分为8×8的小块，每一块内做二维DCT逆变换
yuv_after=cat(3,Y,result,V);      %将经处理的色彩层和两个未处理的层合成
rgb=ycbcr2rgb(yuv_after);    %使YUV图像变回RGB图像
imwrite(rgb,'markresule.jpg','jpg');      %存储添加水印后的图像
subplot(2,3,3),imshow(rgb,[]),title('嵌入水印的图像');    %显示添加水印后的图像

%攻击图像，测试其鲁棒性
disp('请选择对图像的攻击方式：');
disp('1.添加白噪声');
disp('2.对图像进行部分剪切');
disp('3.将图像旋转十度');
disp('4.将图像压缩处理');
disp('5.不处理图像，直接显示提取水印');
disp('输入其它数字则直接显示提取水印');
choice=input('请输入选择：');
figure(1);
switch choice        %读入输入的选择  withmark为等待提取水印的图像
case 1
result_1=rgb;
noise=10*randn(size(result_1));    %生成随机白噪声
result_1=double(result_1)+noise;        %添加白噪声
withmark=uint8(result_1);
subplot(2,3,4);
imshow(withmark,[]);
title('加入白噪声后的图像');     %显示加了白噪声的图像
case 2
result_2=rgb;
A=result_2(:,:,1);
B=result_2(:,:,2);
C=result_2(:,:,3);
A(1:64,1:400)=512;   %使图像上方被剪裁
B(1:64,1:400)=512;   %分别对三个图层操作
C(1:64,1:400)=512; 
result_2=cat(3,A,B,C);
subplot(2,3,4);
imshow(result_2);
title('上方剪切后图像');
figure(1);
withmark=result_2;
case 3
result_3=rgb;
result_3=imrotate(rgb,10,'bilinear','crop');   %最邻近线性插值算法旋转10度
subplot(2,3,4);
imshow(result_3);
title('旋转10度后图像');
withmark=result_3;
case 4
[cA1,cH1,cV1,cD1]=dwt2(rgb,'Haar');    %通过小波变换对图像进行压缩
cA1=compress(cA1);
cH1=compress(cH1);
cV1=compress(cV1);
cD1=compress(cD1);
result_4=idwt2(cA1,cH1,cV1,cD1,'Haar');
result_4=uint8(result_4);
subplot(2,3,4);
imshow(result_4);
title('经小波压缩后的图像');
figure(1);
withmark=result_4;
case 5
subplot(2,3,4);
imshow(rgb,[]);
title('未受攻击的水印图像');
withmark=rgb;
otherwise
disp('选择无效，图像未受攻击，直接提取水印');
subplot(2,3,4);
imshow(rgb,[]);
title('未受攻击的水印图像');
withmark=rgb;
end

% ↓ 这里应该是要先变回YUV模式，我大意了_(:з」∠)_
U_2=withmark(:,:,2);         %取出withmark图像的灰度层
after_2=blkproc(U_2,[8,8],'dct2');   %此步开始提取水印，将灰度层分块进行DCT变换
p=zeros(1,8);        %初始化提取数值用的矩阵
for i=1:marksize(1)
for j=1:marksize(2)
x=(i-1)*8;y=(j-1)*8;
p(1)=after_2(x+1,y+8);         %将之前改变过数值的点的数值提取出来
p(2)=after_2(x+2,y+7);
p(3)=after_2(x+3,y+6);
p(4)=after_2(x+4,y+5);
p(5)=after_2(x+5,y+4);
p(6)=after_2(x+6,y+3);
p(7)=after_2(x+7,y+2);
p(8)=after_2(x+8,y+1);
if corr2(p,k1)>corr2(p,k2)  %corr2计算两个矩阵的相似度，越接近1相似度越大
mark_2(i,j)=0;              %比较提取出来的数值与随机频率k1和k2的相似度，还原水印图样
else
mark_2(i,j)=1;
end
end
end
subplot(2,3,5);
imshow(mark_2,[]),title('提取出的水印');
subplot(2,3,6);
imshow(mark),title('原嵌入水印');