load('crop.mat')
Roc = (50-Roc)/50;
% Roc(21:30)=1;
Roc2(1:4)=1-2/500;Roc2(5:27)=1-1/500;Roc2(28:30)=1;
ft=zeros(30,1);ft2=zeros(30,1);
ft(1)=1/power(2,64);ft2(1)=1/power(2,100);
for i=2:30
   ft(i)=ft(i-1)+nchoosek(64,i-1)/power(2,64); 
end
for i=2:30
   ft2(i)=ft2(i-1)+nchoosek(100,i-1)/power(2,100); 
end
% roc3=[-0.25,-0.32,-0.40,-0.50,-0.66,-0.88,-1,-1.5,-2,-2.67,-3.5,-4.22,-5.45,-6.25,-7.92,-10];%//rotate30
roc3 = [-0.23 -0.24 -0.25 -0.26 -0.5 -0.75 -1 -1.45 -1.95 -2.55 -3.1 -3.97 -4.92 -5.87 -7 -8.1];%jpeg60
% roc3 = [-0.30,-0.40,-0.50,-0.60,-0.80,-1.0,-1.5,-2,-2.9,-3.8,-4.5,-5.5,-6.4,-7.25,-8.3,-10];%//crop70
% roc3 =   [-0.40,-0.50,-0.60,-0.70 -1.00 -1.45 -1.95 -2.25 -3.25 -4.10
% -5.15 -6.00 -7.10 -8.05 -9 -10.5];//scale
rate=-15:1:0;ft3 = power(10,rate);
Roc3 = 1-power(10,roc3);
figure,
plot(ft2,Roc2,'-sr');
hold on,plot(ft,Roc,'-b.');
plot(ft3,Roc3,'-cx');
legend('Proposed Scheme','Method in [9]','Method in [12]','Location','SouthEast');
xlabel('True Negative Rate'),ylabel('True Positive Rate');
set(gcf,'position',[100,100,500,500]);
set(gca, 'Fontname', 'Times New Roman','FontSize',12);
 saveas(gca,'jpeg.emf')
saveas(gca,'jpeg.fig')