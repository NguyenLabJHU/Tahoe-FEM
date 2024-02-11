
close all
clear all

% I used 8.th GP
n=25;
num=(n-1)*35+2;

format short e
sim88Lc1x2=load('el_punch2x2x2_7.txt');
time88Lc1x2=sim88Lc1x2(:,1);
trs88Lc1x2=sim88Lc1x2(:,num+21);
s_inv88Lc1x2=sim88Lc1x2(:,num+22);
trrel88Lc1x2=sim88Lc1x2(:,num+23);
rel_inv88Lc1x2=sim88Lc1x2(:,num+24);
trm88Lc1x2=sim88Lc1x2(:,num+25); %||trm||
hos_inv88Lc1x2=sim88Lc1x2(:,num+26);

n=19;
num=(n-1)*35+2;
sim88Lc1x4=load('el_punch4x4x2_11.txt');
time88Lc1x4=sim88Lc1x4(:,1);
trs88Lc1x4=sim88Lc1x4(:,num+21);
s_inv88Lc1x4=sim88Lc1x4(:,num+22);
trrel88Lc1x4=sim88Lc1x4(:,num+23);
rel_inv88Lc1x4=sim88Lc1x4(:,num+24);
trm88Lc1x4=sim88Lc1x4(:,num+25); %||trm||
hos_inv88Lc1x4=sim88Lc1x4(:,num+26);






% sim88Lc10new=load('el_columnLc10newbuild_1.txt');
% time88Lc10new=sim88Lc10new(:,1);
% s_inv88Lc10new=sim88Lc10new(:,1460);
% rel_inv88Lc10new=sim88Lc10new(:,1461);
% hos_inv88Lc10new=sim88Lc10new(:,1462);
% trs88Lc10new=sim88Lc10new(:,1463);
% trrel88Lc10new=sim88Lc10new(:,1464);
% trm88Lc10new=sim88Lc10new(:,1465);


% 

% sim88Lc1new=load('el_columnLc1newbuild_1.txt');
% time88Lc1new=sim88Lc1new(:,1);
% s_inv88Lc1new=sim88Lc1new(:,1460);
% rel_inv88Lc1new=sim88Lc1new(:,1461);
% hos_inv88Lc1new=sim88Lc1new(:,1462);
% trs88Lc1new=sim88Lc1new(:,1463);
% trrel88Lc1new=sim88Lc1new(:,1464);
% trm88Lc1new=sim88Lc1new(:,1465);
% 
% sim88Lc100=load('el_columnLc100_1.txt');
% time88Lc100=sim88Lc100(:,1);
% s_inv88Lc100=sim88Lc100(:,1460);
% rel_inv88Lc100=sim88Lc100(:,1461);
% hos_inv88Lc100=sim88Lc100(:,1462);
% trs88Lc100=sim88Lc100(:,1463);
% trrel88Lc100=sim88Lc100(:,1464);
% trm88Lc100=sim88Lc100(:,1465);
% 
%  sim88Lc100new=load('el_columnLc100newbuild_1.txt');
%  time88Lc100new=sim88Lc100new(:,1);
%  s_inv88Lc100new=sim88Lc100new(:,1460);
%  rel_inv88Lc100new=sim88Lc100new(:,1461);
%  hos_inv88Lc100new=sim88Lc100new(:,1462);
%  trs88Lc100new=sim88Lc100new(:,1463);
%  trrel88Lc100new=sim88Lc100new(:,1464);
% trm88Lc100new=sim88Lc100new(:,1465);

% bdwidth = 5;
% topbdwidth = 30;
% 
% set(0,'Units','pixels') 
% scnsize = get(0,'ScreenSize');
% 
% pos1  = [bdwidth,... 
%  2/3*scnsize(4) + bdwidth,...
%  scnsize(3)/2 - 2*bdwidth,...
%  scnsize(4)/3 - (topbdwidth + bdwidth)];
% pos2 = [pos1(1) + scnsize(3)/2,...
%  pos1(2),...
%  pos1(3),...
%  pos1(4)];
% 
pos1=[10,1700,650,500];

figure('Position',pos1) 
figure(1)
%figure('Name','Results at the closest Gauss point (GP-1) under the corner point ','NumberTitle','off')
%subplot(3,3,1);
plot(time88Lc1x2,s_inv88Lc1x2*1e6,'-sk','LineWidth',1, 'MarkerSize',8)
hold on
plot(time88Lc1x4,s_inv88Lc1x4*1e6,'-+k','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100new,s_inv88Lc100new*1e6,'or','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100,s_inv88Lc100*1e6,'xb','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('||devsigma|| MPa    ','FontSize',16)
h=legend('2x2x2','4x4x2','Location','NorthWest');
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
plot(time88Lc1x2,rel_inv88Lc1x2*1e6,'-sk','LineWidth',1, 'MarkerSize',8)
%subplot(3,3,2);
hold on
plot(time88Lc1x4,rel_inv88Lc1x4*1e6,'-+k','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100new,rel_inv88Lc100new*1e6,'or','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100,rel_inv88Lc100*1e6,'xb','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('||dev(s-sigma)|| MPa    ','FontSize',16)
h=legend('2x2x2','4x4x2','Location','NorthWest');
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
plot(time88Lc1x2,hos_inv88Lc1x2*1e6,'-sk','LineWidth',1, 'MarkerSize',8)
%subplot(3,3,3);
hold on
plot(time88Lc1x4,hos_inv88Lc1x4*1e6,'-+k','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100new,hos_inv88Lc100new*1e6,'or','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100,hos_inv88Lc100*1e6,'xb','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('||devm|| MPa    ','FontSize',16)
h=legend('2x2x2','4x4x2','Location','NorthWest');
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
plot(time88Lc1x2,trs88Lc1x2*1e6,'-ok','LineWidth',1, 'MarkerSize',8)
hold on
plot(time88Lc1x4,trs88Lc1x4*1e6,'-+k','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100new,trs88Lc100new*1e6,'or','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100,trs88Lc100*1e6,'xb','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('tr(sigma) MPa    ','FontSize',16)
h=legend('2x2x2','4x4x2','Location','SouthWest');
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
plot(time88Lc1x2,trrel88Lc1x2*1e6,'-sk','LineWidth',1, 'MarkerSize',8)
hold on
plot(time88Lc1x4,trrel88Lc1x4*1e6,'-+k','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('tr(s-sigma) MPa   ','FontSize',16)
h=legend('2x2x2','4x4x2','Location','NorthWest');
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
plot(time88Lc1x2,trm88Lc1x2*1e6,'-sk','LineWidth',1, 'MarkerSize',8)
hold on
plot(time88Lc1x4,trm88Lc1x4*1e6,'-+k','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('||tr(m)|| MPa    ','FontSize',16)
h=legend('2x2x2','4x4x2','Location','NorthWest');
legend('boxoff')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)


