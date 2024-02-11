
close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% 8 el. mesh %%%%%%%%%%%%%%%%%%%%%%%%%
% I used 8.th GP
n=25;
num=(n-1)*35+2;
% For the Gauss point, closest to the bottom surface of the 8.th element
% IP=27


format short e
sim8Lc1=load('el_column8Lc1_8.txt');
time8Lc1=sim8Lc1(:,1);
trs8Lc1=sim8Lc1(:,num+21);
s_inv8Lc1=sim8Lc1(:,num+22);
trrel8Lc1=sim8Lc1(:,num+23);
rel_inv8Lc1=sim8Lc1(:,num+24);
trm8Lc1=sim8Lc1(:,num+25); %
hos_inv8Lc1=sim8Lc1(:,num+26);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% 16 el. mesh %%%%%%%%%%%%%%%%%%%%%%%%%
sim16Lc1=load('el_column16Lc1_15.txt');
n=25;
num=(n-1)*35+2;n=25;

% for the element 16, closest IP point to the bottom surface IP=27
%  n=27;
%  num=(n-1)*35+2;n=25;
 
%sim16Lc1=load('el_column16Lc1_16.txt');
time16Lc1=sim16Lc1(:,1);
trs16Lc1=sim16Lc1(:,num+21);
s_inv16Lc1=sim16Lc1(:,num+22);
trrel16Lc1=sim16Lc1(:,num+23);
rel_inv16Lc1=sim16Lc1(:,num+24);
trm16Lc1=sim16Lc1(:,num+25); %||trm||
hos_inv16Lc1=sim16Lc1(:,num+26);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% 32 el. mesh %%%%%%%%%%%%%%%%%%%%%%%%%
n=26; % other IP points can 26, and 27
num=(n-1)*35+2;
% for the element 32, the closest IP point to the other location is IP=26

sim32Lc1=load('el_column32Lc1_29.txt');
% 
% n=26; % other IP points can 26, and 27
% num=(n-1)*35+2;
% sim32Lc1=load('el_column32Lc1_32.txt');
time32Lc1=sim32Lc1(:,1);
trs32Lc1=sim32Lc1(:,num+21);
s_inv32Lc1=sim32Lc1(:,num+22);
trrel32Lc1=sim32Lc1(:,num+23);
rel_inv32Lc1=sim32Lc1(:,num+24);
trm32Lc1=sim32Lc1(:,num+25); %||trm||
hos_inv32Lc1=sim32Lc1(:,num+26);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% 64 el. mesh %%%%%%%%%%%%%%%%%%%%%%%%%
n=27; % other IP points can 26, and 27
num=(n-1)*35+2;

% for the element 64, the closest IP point to the other location is IP=25
% n=25; %
% num=(n-1)*35+2;

sim64Lc1=load('el_column64Lc1_57.txt');
%sim64Lc1=load('el_column64Lc1_64.txt');
time64Lc1=sim64Lc1(:,1);
trs64Lc1=sim64Lc1(:,num+21);
s_inv64Lc1=sim64Lc1(:,num+22);
trrel64Lc1=sim64Lc1(:,num+23);
rel_inv64Lc1=sim64Lc1(:,num+24);
trm64Lc1=sim64Lc1(:,num+25); %||trm||
hos_inv64Lc1=sim64Lc1(:,num+26);


% 
pos1=[10,1700,650,500];

figure('Position',pos1) 
figure(1)
%figure('Name','Results at the closest Gauss point (GP-1) under the corner point ','NumberTitle','off')
%subplot(3,3,1);
plot(time8Lc1,s_inv8Lc1,'-sk','LineWidth',1, 'MarkerSize',8)
hold on
plot(time16Lc1,s_inv16Lc1,'-xk','LineWidth',1, 'MarkerSize',8)
plot(time32Lc1,s_inv32Lc1,'-ok','LineWidth',1, 'MarkerSize',8)
plot(time64Lc1,s_inv64Lc1,'-+k','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('||devsigma|| MPa    ','FontSize',16)
h=legend('8 el','16 el','32 el','64 el','Location','NorthWest');
legend('boxoff')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)
% 
pos2=[670,1700,650,500];

figure('Position',pos2)
figure(2)
%figure('Name','Results at the closest Gauss point (GP-1) under the corner point ','NumberTitle','off')
plot(time8Lc1,rel_inv8Lc1,'-sk','LineWidth',1, 'MarkerSize',8)
%subplot(3,3,2);
hold on
plot(time16Lc1,rel_inv16Lc1,'-xk','LineWidth',1, 'MarkerSize',8)
plot(time32Lc1,rel_inv32Lc1,'-ok','LineWidth',1, 'MarkerSize',8)
plot(time64Lc1,rel_inv64Lc1,'-+k','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('||dev(s-sigma)|| MPa    ','FontSize',16)
h=legend('8 el','16 el','32 el','64 el','Location','NorthWest');
legend('boxoff')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)
% 



pos3=[1340,1700,650,500];
figure('Position',pos3)
figure(3)
%figure('Name','Results at the closest Gauss point (GP-1) under the corner point ','NumberTitle','off')
plot(time8Lc1,hos_inv8Lc1,'-sk','LineWidth',1, 'MarkerSize',8)
%subplot(3,3,2);
hold on
plot(time16Lc1,hos_inv16Lc1,'-xk','LineWidth',1, 'MarkerSize',8)
plot(time32Lc1,hos_inv32Lc1,'-ok','LineWidth',1, 'MarkerSize',8)
plot(time64Lc1,hos_inv64Lc1,'-+k','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('||devm|| MPa    ','FontSize',16)
h=legend('8 el','16 el','32 el','64 el','Location','NorthWest');
legend('boxoff')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)
% 
pos4=[10,10,650,500];
figure('Position',pos4)
figure(4)
%figure('Name','Results at the closest Gauss point (GP-1) under the corner point ','NumberTitle','off')
plot(time8Lc1,trs8Lc1,'-sk','LineWidth',1, 'MarkerSize',8)
%subplot(3,3,2);
hold on
plot(time16Lc1,trs16Lc1,'-xk','LineWidth',1, 'MarkerSize',8)
plot(time32Lc1,trs32Lc1,'-ok','LineWidth',1, 'MarkerSize',8)
plot(time64Lc1,trs64Lc1,'-+k','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('trsigma MPa    ','FontSize',16)
h=legend('8 el','16 el','32 el','64 el','Location','NorthWest');
legend('boxoff')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)
% 
pos5=[670,10,650,500];
figure('Position',pos5)
figure(5)
%figure('Name','Results at the closest Gauss point (GP-1) under the corner point ','NumberTitle','off')
plot(time8Lc1,trrel8Lc1,'-sk','LineWidth',1, 'MarkerSize',8)
%subplot(3,3,2);
hold on
plot(time16Lc1,trrel16Lc1,'-xk','LineWidth',1, 'MarkerSize',8)
plot(time32Lc1,trrel32Lc1,'-ok','LineWidth',1, 'MarkerSize',8)
plot(time64Lc1,trrel64Lc1,'-+k','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('trrel MPa    ','FontSize',16)
h=legend('8 el','16 el','32 el','64 el','Location','NorthWest');
legend('boxoff')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)
% 
pos6=[1340,10,650,500];
figure('Position',pos6)
figure(6)
%figure('Name','Results at the closest Gauss point (GP-1) under the corner point ','NumberTitle','off')
plot(time8Lc1,trm8Lc1,'-sk','LineWidth',1, 'MarkerSize',8)
%subplot(3,3,2);
hold on
plot(time16Lc1,trm16Lc1,'-xk','LineWidth',1, 'MarkerSize',8)
plot(time32Lc1,trm32Lc1,'-ok','LineWidth',1, 'MarkerSize',8)
plot(time64Lc1,trm64Lc1,'-+k','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('trrel MPa    ','FontSize',16)
h=legend('8 el','16 el','32 el','64 el','Location','NorthWest');
legend('boxoff')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)

