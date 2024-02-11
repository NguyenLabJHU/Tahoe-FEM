
close all
clear all

% I used 8.th GP

format short e




% from top to bottom, nodes: 67,121,139
% col8=load('el_column8Lc1nodal_67.txt');
% time8=col8(:,1);
% trs8=col8(:,38);
% devsigma8=col8(:,39);
% trrel8=col8(:,40);
% devrel8=col8(:,41);
% trm8=col8(:,42);
% devm8=col8(:,43);
% 
% %from top to bottom, nodes:139,247,283
% col16=load('el_column16Lc1nodal_139.txt');
% time16=col16(:,1);
% trs16=col16(:,38);
% devsigma16=col16(:,39);
% trrel16=col16(:,40);
% devrel16=col16(:,41);
% trm16=col16(:,42);
% devm16=col16(:,43);


% New stress points for 32 elements: 571, 553, 535
% for 64 elements: 1147, 1111, 1075
% for 128 elements: 2299, 2227, 2155

%from top to bottom, nodes:283,499,571
col32=load('el_column32Lc1nodal_499.txt');
time32=col32(:,1);
trs32=col32(:,38);
devsigma32=col32(:,39);
trrel32=col32(:,40);
devrel32=col32(:,41);
trm32=col32(:,42);
devm32=col32(:,43);

%from top to bottom, nodes:571,1003,1147
col64=load('el_column64Lc1nodal_1003.txt');
time64=col64(:,1);
trs64=col64(:,38);
devsigma64=col64(:,39);
trrel64=col64(:,40);
devrel64=col64(:,41);
trm64=col64(:,42);
devm64=col64(:,43);
%from top to bottom, nodes:1147,2011,2299
col128=load('el_column128Lc1nodal_2011.txt');
time128=col128(:,1);
trs128=col128(:,38);
devsigma128=col128(:,39);
trrel128=col128(:,40);
devrel128=col128(:,41);
trm128=col128(:,42);
devm128=col128(:,43);



pos1=[10,1700,650,500];

figure('Position',pos1) 
figure(1)
%figure('Name','Results at the closest Gauss point (GP-1) under the corner point ','NumberTitle','off')
%subplot(3,3,1);
%%plot(time8,devsigma8*1e6,'-sk','LineWidth',1, 'MarkerSize',8)

%plot(time16,devsigma16*1e6,'-xk','LineWidth',1, 'MarkerSize',8)
plot(time32,devsigma32*1e6,'-ok','LineWidth',1, 'MarkerSize',8)
hold on
plot(time64,devsigma64*1e6,'-^k','LineWidth',1, 'MarkerSize',8)
plot(time128,devsigma128*1e6,'-vk','LineWidth',1, 'MarkerSize',8)

%plot(time88Lc1x4,s_inv88Lc1x4*1e6,'-+k','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100new,s_inv88Lc100new*1e6,'or','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100,s_inv88Lc100*1e6,'xb','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('||devsigma|| MPa    ','FontSize',16)
h=legend('col32','col64','col128','Location','NorthWest');
legend('boxoff')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)
% 
pos2=[670,1700,650,500];
% 
figure('Position',pos2)
 figure(2)
%figure('Name','Results at the closest Gauss point (GP-1) under the corner point ','NumberTitle','off')
%subplot(3,3,1);
%plot(time8,devrel8*1e6,'-sk','LineWidth',1, 'MarkerSize',8)

%plot(time16,devrel16*1e6,'-xk','LineWidth',1, 'MarkerSize',8)
plot(time32,devrel32*1e6,'-ok','LineWidth',1, 'MarkerSize',8)
hold on
plot(time64,devrel64*1e6,'-^k','LineWidth',1, 'MarkerSize',8)
plot(time128,devrel128*1e6,'-vk','LineWidth',1, 'MarkerSize',8)

%plot(time88Lc1x4,s_inv88Lc1x4*1e6,'-+k','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100new,s_inv88Lc100new*1e6,'or','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100,s_inv88Lc100*1e6,'xb','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('||devrel|| MPa    ','FontSize',16)
h=legend('col32','col64','col128','Location','NorthWest');
legend('boxoff')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)


pos3=[1340,1700,650,500];
 figure('Position',pos3)
 figure(3)
%plot(time8,devm8*1e6,'-sk','LineWidth',1, 'MarkerSize',8)

%plot(time16,devm16*1e6,'-xk','LineWidth',1, 'MarkerSize',8)
plot(time32,devm32*1e6,'-ok','LineWidth',1, 'MarkerSize',8)
hold on
plot(time64,devm64*1e6,'-^k','LineWidth',1, 'MarkerSize',8)
plot(time128,devm128*1e6,'-vk','LineWidth',1, 'MarkerSize',8)
%plot(time88Lc1x4,s_inv88Lc1x4*1e6,'-+k','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100new,s_inv88Lc100new*1e6,'or','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100,s_inv88Lc100*1e6,'xb','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('||devm|| MPa    ','FontSize',16)
h=legend('col32','col64','col128','Location','NorthWest');
legend('boxoff')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)






 pos4=[10,10,650,500];
 figure('Position',pos4)
 figure(4)
%plot(time8,trs8*1e6,'-sk','LineWidth',1, 'MarkerSize',8)
%hold on
%plot(time16,trs16*1e6,'-xk','LineWidth',1, 'MarkerSize',8)
plot(time32,trs32*1e6,'-ok','LineWidth',1, 'MarkerSize',8)
hold on
plot(time64,trs64*1e6,'-^k','LineWidth',1, 'MarkerSize',8)
plot(time128,trs128*1e6,'-vk','LineWidth',1, 'MarkerSize',8)

%plot(time88Lc1x4,s_inv88Lc1x4*1e6,'-+k','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100new,s_inv88Lc100new*1e6,'or','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100,s_inv88Lc100*1e6,'xb','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('trsigma MPa    ','FontSize',16)
h=legend('col32','col64','col128','Location','NorthWest');
legend('boxoff')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)

% % 
 pos5=[670,10,650,500];
 figure('Position',pos5)
 figure(5)
%plot(time8,trrel8*1e6,'-sk','LineWidth',1, 'MarkerSize',8)
%hold on
%plot(time16,trrel16*1e6,'-xk','LineWidth',1, 'MarkerSize',8)

plot(time32,trrel32*1e6,'-ok','LineWidth',1, 'MarkerSize',8)
hold on
plot(time64,trrel64*1e6,'-^k','LineWidth',1, 'MarkerSize',8)
plot(time128,trrel128*1e6,'-vk','LineWidth',1, 'MarkerSize',8)

%plot(time88Lc1x4,s_inv88Lc1x4*1e6,'-+k','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100new,s_inv88Lc100new*1e6,'or','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100,s_inv88Lc100*1e6,'xb','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('trrel MPa    ','FontSize',16)
h=legend('col32','col64','col128','Location','NorthWest');
legend('boxoff')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)

% % 
pos6=[1340,10,650,500];
 figure('Position',pos6)
 figure(6)
%plot(time8,trm8*1e6,'-sk','LineWidth',1, 'MarkerSize',8)
%hold on
%plot(time16,trm16*1e6,'-xk','LineWidth',1, 'MarkerSize',8)
plot(time32,trm32*1e6,'-ok','LineWidth',1, 'MarkerSize',8)
hold on
plot(time64,trm64*1e6,'-^k','LineWidth',1, 'MarkerSize',8)
plot(time128,trm128*1e6,'-vk','LineWidth',1, 'MarkerSize',8)

%plot(time88Lc1x4,s_inv88Lc1x4*1e6,'-+k','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100new,s_inv88Lc100new*1e6,'or','LineWidth',1, 'MarkerSize',8)
% plot(time88Lc100,s_inv88Lc100*1e6,'xb','LineWidth',1, 'MarkerSize',8)
hold off
xlabel('time','FontSize',16)
ylabel('trm MPa    ','FontSize',16)
h=legend('col32','col64','col128','Location','NorthWest');
legend('boxoff')
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)


