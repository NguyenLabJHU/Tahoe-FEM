
close all
clear all

format short e
%Different elements on 9 elements punch area
%Element on the edge is 541 ---1
%Element in the middle of betwwen two edge elements ( in the middle of
% the symmetric part of punch area) is 631---2
%Element which ist eh middle od the whole punch area ( not only symmetric
% punch area) is 721---3

% Nodes on the edge of the element 1 on the diagonal are
% 4390+1=4391
% 5194+1=5195
% 5188+1=5189


% First element middle of the punch 5993
% First element vertex 5987
% Second element  vertex 5189
% third element vertex 4391

% Enter GP number "n" 
n=25;
num=(n-1)*35+2;

% 
simset8=load('el_cube10_set8all_991.txt');
timeset8=simset8(:,1);
%ip1
s11set8=simset8(:,num);
s22set8=simset8(:,num+1);
s33set8=simset8(:,num+2);
s12set8=simset8(:,num+3);
s13set8=simset8(:,num+4);
s21set8=simset8(:,num+5);
s23set8=simset8(:,num+6);
s31set8=simset8(:,num+7);
s32set8=simset8(:,num+8);
e11set8=simset8(:,num+9);
e22set8=simset8(:,num+10);
e33set8=simset8(:,num+11);
e12set8=simset8(:,num+12);
e13set8=simset8(:,num+13);
e21set8=simset8(:,num+14);
e23set8=simset8(:,num+15);
e31set8=simset8(:,num+16);
e32set8=simset8(:,num+17);

%final
cset8=(simset8(length(timeset8),num+18));
%Zc=(sim(length(time),21));
hcset8=(simset8(length(timeset8),num+19));
Delgammaset8=(simset8(length(timeset8),num+20));
%initial
c0set8=(simset8(1,num+18));
%Zc0=(sim(1,21));
hc0set8=(simset8(1,num+19));
Delgamma0set8=(simset8(1,num+20));

trsigmaset8=simset8(:,num+21);
devsigmaset8=simset8(:,num+22);
trrelset8=simset8(:,num+23);
devrelset8=simset8(:,num+24);
trmset8=simset8(:,num+25); %||trm||
devmset8=simset8(:,num+26);



simset9=load('el_cube10_set9_991.txt');
timeset9=simset9(:,1);
%ip1
s11set9=simset9(:,num);
s22set9=simset9(:,num+1);
s33set9=simset9(:,num+2);
s12set9=simset9(:,num+3);
s13set9=simset9(:,num+4);
s21set9=simset9(:,num+5);
s23set9=simset9(:,num+6);
s31set9=simset9(:,num+7);
s32set9=simset9(:,num+8);
e11set9=simset9(:,num+9);
e22set9=simset9(:,num+10);
e33set9=simset9(:,num+11);
e12set9=simset9(:,num+12);
e13set9=simset9(:,num+13);
e21set9=simset9(:,num+14);
e23set9=simset9(:,num+15);
e31set9=simset9(:,num+16);
e32set9=simset9(:,num+17);

%final
cset9=(simset9(length(timeset9),num+18));
%Zc=(sim(length(time),21));
hcset9=(simset9(length(timeset9),num+19));
Delgammaset9=(simset9(length(timeset9),num+20));
%initial
c0set9=(simset9(1,num+18));
%Zc0=(sim(1,21));
hc0set9=(simset9(1,num+19));
Delgamma0set9=(simset9(1,num+20));

trsigmaset9=simset9(:,num+21);
devsigmaset9=simset9(:,num+22);
trrelset9=simset9(:,num+23);
devrelset9=simset9(:,num+24);
trmset9=simset9(:,num+25); %||trm||
devmset9=simset9(:,num+26);


simset10=load('el_cube10_set10_991.txt');
timeset10=simset10(:,1);

%ip1
s11set10=simset10(:,num);
s22set10=simset10(:,num+1);
s33set10=simset10(:,num+2);
s12set10=simset10(:,num+3);
s13set10=simset10(:,num+4);
s21set10=simset10(:,num+5);
s23set10=simset10(:,num+6);
s31set10=simset10(:,num+7);
s32set10=simset10(:,num+8);
e11set10=simset10(:,num+9);
e22set10=simset10(:,num+10);
e33set10=simset10(:,num+11);
e12set10=simset10(:,num+12);
e13set10=simset10(:,num+13);
e21set10=simset10(:,num+14);
e23set10=simset10(:,num+15);
e31set10=simset10(:,num+16);
e32set10=simset10(:,num+17);

%final
cset10=(simset10(length(timeset10),num+18));
%Zc=(sim(length(time),21));
hcset10=(simset10(length(timeset10),num+19));
Delgammaset10=(simset10(length(timeset10),num+20));
%initial
c0set10=(simset10(1,num+18));
%Zc0=(sim(1,21));
hc0set10=(simset10(1,num+19));
Delgamma0set10=(simset10(1,num+20));

trsigmaset10=simset10(:,num+21);
devsigmaset10=simset10(:,num+22);
trrelset10=simset10(:,num+23);
devrelset10=simset10(:,num+24);
trmset10=simset10(:,num+25); %||trm||
devmset10=simset10(:,num+26);


simset11=load('el_cube10_set11_991.txt');
timeset11=simset11(:,1);

%ip1
s11set11=simset11(:,num);
s22set11=simset11(:,num+1);
s33set11=simset11(:,num+2);
s12set11=simset11(:,num+3);
s13set11=simset11(:,num+4);
s21set11=simset11(:,num+5);
s23set11=simset11(:,num+6);
s31set11=simset11(:,num+7);
s32set11=simset11(:,num+8);
e11set11=simset11(:,num+9);
e22set11=simset11(:,num+10);
e33set11=simset11(:,num+11);
e12set11=simset11(:,num+12);
e13set11=simset11(:,num+13);
e21set11=simset11(:,num+14);
e23set11=simset11(:,num+15);
e31set11=simset11(:,num+16);
e32set11=simset11(:,num+17);

%final
cset11=(simset11(length(timeset11),num+18));
%Zc=(sim(length(time),21));
hcset11=(simset11(length(timeset11),num+19));
Delgammaset11=(simset11(length(timeset11),num+20));
%initial
c0set11=(simset11(1,num+18));
%Zc0=(sim(1,21));
hc0set11=(simset11(1,num+19));
Delgamma0set11=(simset11(1,num+20));

trsigmaset11=simset11(:,num+21);
devsigmaset11=simset11(:,num+22);
trrelset11=simset11(:,num+23);
devrelset11=simset11(:,num+24);
trmset11=simset11(:,num+25); %||trm||
devmset11=simset11(:,num+26);



simset12=load('el_cube10_set12_991.txt');
timeset12=simset12(:,1);

%ip1
s11set12=simset12(:,num);
s22set12=simset12(:,num+1);
s33set12=simset12(:,num+2);
s12set12=simset12(:,num+3);
s13set12=simset12(:,num+4);
s21set12=simset12(:,num+5);
s23set12=simset12(:,num+6);
s31set12=simset12(:,num+7);
s32set12=simset12(:,num+8);
e11set12=simset12(:,num+9);
e22set12=simset12(:,num+10);
e33set12=simset12(:,num+11);
e12set12=simset12(:,num+12);
e13set12=simset12(:,num+13);
e21set12=simset12(:,num+14);
e23set12=simset12(:,num+15);
e31set12=simset12(:,num+16);
e32set12=simset12(:,num+17);

%final
cset12=(simset12(length(timeset12),num+18));
%Zc=(sim(length(time),21));
hcset12=(simset12(length(timeset12),num+19));
Delgammaset12=(simset12(length(timeset12),num+20));
%initial
c0set12=(simset12(1,num+18));
%Zc0=(sim(1,21));
hc0set12=(simset12(1,num+19));
Delgamma0set12=(simset12(1,num+20));

trsigmaset12=simset12(:,num+21);
devsigmaset12=simset12(:,num+22);
trrelset12=simset12(:,num+23);
devrelset12=simset12(:,num+24);
trmset12=simset12(:,num+25); %||trm||
devmset12=simset12(:,num+26);



% sim33=load('el_cube33M_25.txt');
% time33=sim33(:,1);
% trmall33=sim33(:,601);
% 
% figure(1)
% plot(timeFSE,s_invFSE*1e-6,'-k','LineWidth',1, 'MarkerSize',14)
% hold on
% plot(timekappanu,s_invkappanu*1e-6,'-ok','LineWidth',1, 'MarkerSize',14)
% plot(timeetakappanu,s_invetakappanu*1e-6,'-vk','LineWidth',1, 'MarkerSize',14)
% plot(timenoTau7,s_invnoTau7*1e-6,'-^k','LineWidth',1, 'MarkerSize',14)
% plot(timeall,s_invall*1e-6,'-+k','LineWidth',1, 'MarkerSize',14)
% hold off
% 
%  xlabel('time','FontSize',18)
%  ylabel('||devsigma||  (MPa)','FontSize',18)
%  h=legend('Finite strain elasticity (FSE)','Micromorphic FSE with \kappa,\nu','Micromorphic FSE with \eta,\kappa,\nu','Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma','Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma and \tau_{7}');
% %axis([0 2.5 0 30])
% %set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
% %set(gca,'YTickLabel',{'0','2','4','6','8','10'})
% set(gca,'FontName','Helvetica','FontSize',22)


figure(1)
plot(timeset12,devmset12*1e-6,'-^k','LineWidth',1, 'MarkerSize',14)
hold on
plot(timeset11,devmset11*1e-6,'-xk','LineWidth',1, 'MarkerSize',14)
plot(timeset8,devmset8*1e-6,'-vk','LineWidth',1, 'MarkerSize',14)
plot(timeset9,devmset9*1e-6,'-ok','LineWidth',1, 'MarkerSize',14)
plot(timeset10,devmset10*1e-6,'-+k','LineWidth',1, 'MarkerSize',14)
 xlabel('time','FontSize',18)
 ylabel('||trm|| (MPa)','FontSize',18)
h=legend('Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma,\tau_{7}, punch area s_{5 } ',...
'Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma,\tau_{7}, punch area s_{4 } ',...
'Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma,\tau_{7}, punch area s_{1 } ',...
'Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma,\tau_{7}, punch area s_{2 } ',...
'Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma,\tau_{7}, punch area s_{3 } ','Location','NorthWest');
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',22)

hold off


figure(2)
plot(timeset12,trmset12*1e-6,'-^k','LineWidth',1, 'MarkerSize',14)
hold on
plot(timeset11,trmset11*1e-6,'-xk','LineWidth',1, 'MarkerSize',14)
plot(timeset8,trmset8*1e-6,'-vk','LineWidth',1, 'MarkerSize',14)
plot(timeset9,trmset9*1e-6,'-ok','LineWidth',1, 'MarkerSize',14)
plot(timeset10,trmset10*1e-6,'-+k','LineWidth',1, 'MarkerSize',14)
 xlabel('time','FontSize',18)
 ylabel('||devm|| (MPa)','FontSize',18)
h=legend('Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma,\tau_{7}, punch area s_{5 } ',...
'Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma,\tau_{7}, punch area s_{4 } ',...
'Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma,\tau_{7}, punch area s_{1 } ',...
'Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma,\tau_{7}, punch area s_{2 } ',...
'Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma,\tau_{7}, punch area s_{3 } ','Location','NorthWest');
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',22)

hold off