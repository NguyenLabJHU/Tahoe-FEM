
close all
clear all

% I used 8.th GP

format short e
sim88Lc01=load('el_columnLc01_1.txt');
time88Lc01=sim88Lc01(:,1);
s_inv88Lc01=sim88Lc01(:,213);
rel_inv88Lc01=sim88Lc01(:,215);
hos_inv88Lc01=sim88Lc01(:,217);
trs88Lc01=sim88Lc01(:,212);
trrel88Lc01=sim88Lc01(:,214);
trm88Lc01=sim88Lc01(:,216);
% 
sim88Lc1=load('el_columnLc1_1.txt');
time88Lc1=sim88Lc1(:,1);
s_inv88Lc1=sim88Lc1(:,213);
rel_inv88Lc1=sim88Lc1(:,215);
hos_inv88Lc1=sim88Lc1(:,217);
trs88Lc1=sim88Lc1(:,212);
trrel88Lc1=sim88Lc1(:,214);
trm88Lc1=sim88Lc1(:,216);
% 
sim88Lc10=load('el_columnLc10_1.txt');
time88Lc10=sim88Lc10(:,1);
s_inv88Lc10=sim88Lc10(:,213);
rel_inv88Lc10=sim88Lc10(:,215);
hos_inv88Lc10=sim88Lc10(:,217);
trs88Lc10=sim88Lc10(:,212);
trrel88Lc10=sim88Lc10(:,214);
trm88Lc10=sim88Lc10(:,216);
%
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
plot(time88Lc10,s_inv88Lc10*1e6,'-sk','LineWidth',1, 'MarkerSize',8)
hold on
plot(time88Lc1,s_inv88Lc1*1e6,'-+k','LineWidth',1, 'MarkerSize',8)
plot(time88Lc01,s_inv88Lc01*1e6,'-ok','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100new,s_inv88Lc100new*1e6,'or','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100,s_inv88Lc100*1e6,'xb','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('||devsigma|| MPa    ','FontSize',16)
h=legend('Lcten','Lcone','Lcptone','Location','NorthWest');
legend('boxoff')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)

pos2=[670,1700,650,500];

figure('Position',pos2)
figure(2)
%figure('Name','Results at the closest Gauss point (GP-1) under the corner point ','NumberTitle','off')
plot(time88Lc10,rel_inv88Lc10*1e6,'-sk','LineWidth',1, 'MarkerSize',8)
%subplot(3,3,2);
hold on
plot(time88Lc1,rel_inv88Lc1*1e6,'-+k','LineWidth',1, 'MarkerSize',8)
plot(time88Lc01,rel_inv88Lc01*1e6,'-ok','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100new,rel_inv88Lc100new*1e6,'or','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100,rel_inv88Lc100*1e6,'xb','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('||dev(s-sigma)|| MPa    ','FontSize',16)
h=legend('Lcten','Lcone','Lcptone','Location','NorthWest');
legend('boxoff')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)

pos3=[1340,1700,650,500];
figure('Position',pos3)
figure(3)
%figure('Name','Results at the closest Gauss point (GP-1) under the corner point ','NumberTitle','off')
plot(time88Lc10,hos_inv88Lc10*1e6,'-sk','LineWidth',1, 'MarkerSize',8)
%subplot(3,3,3);
hold on
plot(time88Lc1,hos_inv88Lc1*1e6,'-+k','LineWidth',1, 'MarkerSize',8)
plot(time88Lc01,hos_inv88Lc01*1e6,'-ok','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100new,hos_inv88Lc100new*1e6,'or','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100,hos_inv88Lc100*1e6,'xb','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('||m|| MPa    ','FontSize',16)
h=legend('Lcten','Lcone','Lcptone','Location','NorthWest');
legend('boxoff')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)

pos4=[10,10,650,500];
figure('Position',pos4)
figure(4)
%figure('Name','Results at the closest Gauss point (GP-1) under the corner point ','NumberTitle','off')
plot(time88Lc01,trs88Lc01*1e6,'-ok','LineWidth',1, 'MarkerSize',8)
hold on
plot(time88Lc1,trs88Lc1*1e6,'-+k','LineWidth',1, 'MarkerSize',8)
plot(time88Lc10,trs88Lc10*1e6,'-sk','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100new,trs88Lc100new*1e6,'or','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100,trs88Lc100*1e6,'xb','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('tr(sigma) MPa    ','FontSize',16)
h=legend('Lcptone','Lcone','Lcten','Location','SouthWest');
legend('boxoff')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)

pos5=[670,10,650,500];
figure('Position',pos5)
figure(5)
%figure('Name','Results at the closest Gauss point (GP-1) under the corner point ','NumberTitle','off')
plot(time88Lc10,trrel88Lc10*1e6,'-sk','LineWidth',1, 'MarkerSize',8)
hold on
plot(time88Lc1,trrel88Lc1*1e6,'-+k','LineWidth',1, 'MarkerSize',8)
plot(time88Lc01,trrel88Lc01*1e6,'-ok','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100new,trrel88Lc100new*1e6,'or','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100,trrel88Lc100*1e6,'xb','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('tr(s-sigma) MPa   ','FontSize',16)
h=legend('Lcten','Lcone','Lcptone','Lchundred','Location','NorthWest');
legend('boxoff')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)

pos6=[1340,10,650,500];
figure('Position',pos6)
figure(6)
%figure('Name','Results at the closest Gauss point (GP-1) under the corner point ','NumberTitle','off')
plot(time88Lc10,trm88Lc10*1e6,'-sk','LineWidth',1, 'MarkerSize',8)
hold on
plot(time88Lc1,trm88Lc1*1e6,'-+k','LineWidth',1, 'MarkerSize',8)
plot(time88Lc01,trm88Lc01*1e6,'-ok','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100new,trm88Lc100new*1e6,'or','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100,trm88Lc100*1e6,'xb','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('||tr(m)|| MPa    ','FontSize',16)
h=legend('Lcten','Lcone','Lcptone','Location','NorthWest');
legend('boxoff')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)


