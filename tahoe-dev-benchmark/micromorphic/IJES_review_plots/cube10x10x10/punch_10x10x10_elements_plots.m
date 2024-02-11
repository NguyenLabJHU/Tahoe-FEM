
close all
clear all

format short e
% Enter GP number "n" 
n=25;
num=(n-1)*35+2;


simFSE=load('el_cube10_FSE_991.txt');
timeFSE=simFSE(:,1);
s11FSE=simFSE(:,578);
s22FSE=simFSE(:,579);
s33FSE=simFSE(:,580);
s12FSE=simFSE(:,581);
s13FSE=simFSE(:,582);
s21FSE=simFSE(:,583);
s23FSE=simFSE(:,584);
s31FSE=simFSE(:,585);
s32FSE=simFSE(:,586);
e11FSE=simFSE(:,587);
e22FSE=simFSE(:,588);
e33FSE=simFSE(:,589);
e12FSE=simFSE(:,590);
e13FSE=simFSE(:,591);
e21FSE=simFSE(:,592);
e23FSE=simFSE(:,593);
e31FSE=simFSE(:,594);
e32FSE=simFSE(:,595);
s_invFSE=simFSE(:,596);
rel_invFSE=simFSE(:,597);
hos_invFSE=simFSE(:,598);
trsFSE=simFSE(:,599);
trrelFSE=simFSE(:,600);
trmFSE=simFSE(:,601);

simkappanu=load('el_cube10_kappanu_991.txt');
timekappanu=simkappanu(:,1);
s11kappanu=simkappanu(:,num);
s22kappanu=simkappanu(:,num+1);
s33kappanu=simkappanu(:,num+2);
s12kappanu=simkappanu(:,num+3);
s13kappanu=simkappanu(:,num+4);
s21kappanu=simkappanu(:,num+5);
s23kappanu=simkappanu(:,num+6);
s31kappanu=simkappanu(:,num+7);
s32kappanu=simkappanu(:,num+8);
e11kappanu=simkappanu(:,num+9);
e22kappanu=simkappanu(:,num+10);
e33kappanu=simkappanu(:,num+11);
e12kappanu=simkappanu(:,num+12);
e13kappanu=simkappanu(:,num+13);
e21kappanu=simkappanu(:,num+14);
e23kappanu=simkappanu(:,num+15);
e31kappanu=simkappanu(:,num+16);
e32kappanu=simkappanu(:,num+17);

%final
ckappanu=(simkappanu(length(timekappanu),num+18));
%Zc=(sim(length(time),21));
hckappanu=(simkappanu(length(timekappanu),num+19));
Delgammakappanu=(simkappanu(length(timekappanu),num+20));
%initial
c0kappanu=(simkappanu(1,num+18));
%Zc0=(sim(1,21));
hc0kappanu=(simkappanu(1,num+19));
Delgamma0kappanu=(simkappanu(1,num+20));

trsigmakappanu=simkappanu(:,num+21);
devsigmakappanu=simkappanu(:,num+22);
trrelkappanu=simkappanu(:,num+23);
devrelkappanu=simkappanu(:,num+24);
trmkappanu=simkappanu(:,num+25); %||trm||
devmkappanu=simkappanu(:,num+26);



simetakappanu=load('el_cube10_etakappanu_991.txt');
timeetakappanu=simetakappanu(:,1);
s11etakappanu=simetakappanu(:,num);
s22etakappanu=simetakappanu(:,num+1);
s33etakappanu=simetakappanu(:,num+2);
s12etakappanu=simetakappanu(:,num+3);
s13etakappanu=simetakappanu(:,num+4);
s21etakappanu=simetakappanu(:,num+5);
s23etakappanu=simetakappanu(:,num+6);
s31etakappanu=simetakappanu(:,num+7);
s32etakappanu=simetakappanu(:,num+8);
e11etakappanu=simetakappanu(:,num+9);
e22etakappanu=simetakappanu(:,num+10);
e33etakappanu=simetakappanu(:,num+11);
e12etakappanu=simetakappanu(:,num+12);
e13etakappanu=simetakappanu(:,num+13);
e21etakappanu=simetakappanu(:,num+14);
e23etakappanu=simetakappanu(:,num+15);
e31etakappanu=simetakappanu(:,num+16);
e32etakappanu=simetakappanu(:,num+17);

%final
cetakappanu=(simetakappanu(length(timeetakappanu),num+18));
%Zc=(sim(length(time),21));
hcetakappanu=(simetakappanu(length(timeetakappanu),num+19));
Delgammaetakappanu=(simetakappanu(length(timeetakappanu),num+20));
%initial
c0etakappanu=(simetakappanu(1,num+18));
%Zc0=(sim(1,21));
hc0etakappanu=(simetakappanu(1,num+19));
Delgamma0etakappanu=(simetakappanu(1,num+20));

trsigmaetakappanu=simetakappanu(:,num+21);
devsigmaetakappanu=simetakappanu(:,num+22);
trreletakappanu=simetakappanu(:,num+23);
devreletakappanu=simetakappanu(:,num+24);
trmetakappanu=simetakappanu(:,num+25); %||trm||
devmetakappanu=simetakappanu(:,num+26);




simnoTau7=load('el_cube10_noTau7_991.txt');
timenoTau7=simnoTau7(:,1);
s11noTau7=simnoTau7(:,num);
s22noTau7=simnoTau7(:,num+1);
s33noTau7=simnoTau7(:,num+2);
s12noTau7=simnoTau7(:,num+3);
s13noTau7=simnoTau7(:,num+4);
s21noTau7=simnoTau7(:,num+5);
s23noTau7=simnoTau7(:,num+6);
s31noTau7=simnoTau7(:,num+7);
s32noTau7=simnoTau7(:,num+8);
e11noTau7=simnoTau7(:,num+9);
e22noTau7=simnoTau7(:,num+10);
e33noTau7=simnoTau7(:,num+11);
e12noTau7=simnoTau7(:,num+12);
e13noTau7=simnoTau7(:,num+13);
e21noTau7=simnoTau7(:,num+14);
e23noTau7=simnoTau7(:,num+15);
e31noTau7=simnoTau7(:,num+16);
e32noTau7=simnoTau7(:,num+17);
s_invnoTau7=simnoTau7(:,num+18);

%final
cnoTau7=(simnoTau7(length(timenoTau7),num+18));
%Zc=(sim(length(time),21));
hcnoTau7=(simnoTau7(length(timenoTau7),num+19));
DelgammanoTau7=(simnoTau7(length(timenoTau7),num+20));
%initial
c0noTau7=(simnoTau7(1,num+18));
%Zc0=(sim(1,21));
hc0noTau7=(simnoTau7(1,num+19));
Delgamma0noTau7=(simnoTau7(1,num+20));

trsigmanoTau7=simnoTau7(:,num+21);
devsigmanoTau7=simnoTau7(:,num+22);
trrelnoTau7=simnoTau7(:,num+23);
devrelnoTau7=simnoTau7(:,num+24);
trmnoTau7=simnoTau7(:,num+25); %||trm||
devmnoTau7=simnoTau7(:,num+26);


simall=load('el_cube10_set8all_991.txt');
timeall=simall(:,1);
s11all=simall(:,num);
s22all=simall(:,num+1);
s33all=simall(:,num+2);
s12all=simall(:,num+3);
s13all=simall(:,num+4);
s21all=simall(:,num+5);
s23all=simall(:,num+6);
s31all=simall(:,num+7);
s32all=simall(:,num+8);
e11all=simall(:,num+9);
e22all=simall(:,num+10);
e33all=simall(:,num+11);
e12all=simall(:,num+12);
e13all=simall(:,num+13);
e21all=simall(:,num+14);
e23all=simall(:,num+15);
e31all=simall(:,num+16);
e32all=simall(:,num+17);

%final
call=(simall(length(timeall),num+18));
%Zc=(sim(length(time),21));
hcall=(simall(length(timeall),num+19));
Delgammaall=(simall(length(timeall),num+20));
%initial
c0all=(simall(1,num+18));
%Zc0=(sim(1,21));
hc0all=(simall(1,num+19));
Delgamma0all=(simall(1,num+20));

trsigmaall=simall(:,num+21);
devsigmaall=simall(:,num+22);
trrelall=simall(:,num+23);
devrelall=simall(:,num+24);
trmall=simall(:,num+25); %||trm||
devmall=simall(:,num+26);

% sim33=load('el_cube33M_25.txt');
% time33=sim33(:,1);
% trmall33=sim33(:,601);

figure(1)
plot(timeFSE,s_invFSE*1e-6,'-k','LineWidth',1, 'MarkerSize',14)
hold on
plot(timekappanu,devsigmakappanu*1e-6,'-ok','LineWidth',1, 'MarkerSize',14)
%hold on
plot(timeetakappanu,devsigmaetakappanu*1e-6,'-vk','LineWidth',1, 'MarkerSize',14)
plot(timenoTau7,devsigmanoTau7*1e-6,'-^k','LineWidth',1, 'MarkerSize',14)
plot(timeall,devsigmaall*1e-6,'-+k','LineWidth',1, 'MarkerSize',14)
hold off

 xlabel('time','FontSize',18)
 ylabel('||devsigma||  (MPa)','FontSize',18)
 h=legend('Finite strain elasticity (FSE)','Micromorphic FSE with \kappa,\nu','Micromorphic FSE with \eta,\kappa,\nu','Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma'...
     ,'Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma and \tau_{7}','Location','NorthWest');
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',22)

figure(2)
plot(timekappanu,devrelkappanu*1e-6,'-ok','LineWidth',1, 'MarkerSize',14)
hold on
plot(timeetakappanu,devreletakappanu*1e-6,'-vk','LineWidth',1, 'MarkerSize',14)
%hold on
plot(timenoTau7,devrelnoTau7*1e-6,'-^k','LineWidth',1, 'MarkerSize',14)
plot(timeall,devrelall*1e-6,'-+k','LineWidth',1, 'MarkerSize',14)
hold off
xlabel('time','FontSize',18)
ylabel('||dev(s-sigma)||    (MPa)','FontSize',18)
h=legend('Micromorphic FSE with \kappa,\nu','Micromorphic FSE with \eta,\kappa,\nu','Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma'...
    ,'Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma and \tau_{7}','Location','NorthWest');
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',22)


% figure(3)
% plot(timeall,devmall*1e-6,'-+k','LineWidth',1, 'MarkerSize',14)
% xlabel('time','FontSize',18)
% ylabel('||dev(m)||    (MPa)','FontSize',18)
%  h=legend('Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma and \tau_{7}');
% %axis([0 2.5 0 30])
% %set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
% %set(gca,'YTickLabel',{'0','2','4','6','8','10'})
% set(gca,'FontName','Helvetica','FontSize',22)


figure(4)
plot(timeFSE,trsFSE*1e-6,'-k','LineWidth',1, 'MarkerSize',14)
 hold on
plot(timekappanu,trsigmakappanu*1e-6,'-ok','LineWidth',1, 'MarkerSize',14)
%hold on
plot(timeetakappanu,trsigmaetakappanu*1e-6,'-vk','LineWidth',1, 'MarkerSize',14)
plot(timenoTau7,trsigmanoTau7*1e-6,'-^k','LineWidth',1, 'MarkerSize',14)
plot(timeall,trsigmaall*1e-6,'-+k','LineWidth',1, 'MarkerSize',14)
hold off

 xlabel('time','FontSize',18)
 ylabel('trsigma (MPa)','FontSize',18)
 h=legend('Finite strain elasticity (FSE)','Micromorphic FSE with \kappa,\nu','Micromorphic FSE with \eta,\kappa,\nu','Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma'...
     ,'Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma and \tau_{7}','Location','SouthWest');
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',22)


figure(5)
plot(timekappanu,trrelkappanu*1e-6,'-ok','LineWidth',1, 'MarkerSize',14)
hold on
plot(timeetakappanu,trreletakappanu*1e-6,'-vk','LineWidth',1, 'MarkerSize',14)
plot(timenoTau7,trrelnoTau7*1e-6,'-^k','LineWidth',1, 'MarkerSize',14)
plot(timeall,trrelall*1e-6,'-+k','LineWidth',1, 'MarkerSize',14)
hold off

 xlabel('time','FontSize',18)
 ylabel('trrel (MPa)','FontSize',18)
 h=legend('Micromorphic FSE with \kappa,\nu','Micromorphic FSE with \eta,\kappa,\nu','Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma'...
     ,'Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma and \tau_{7}','Location','NorthWest');
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',22)

% figure(6)
% plot(timeall,trmall*1e-6,'-+k','LineWidth',1, 'MarkerSize',14)
%  xlabel('time','FontSize',18)
%  ylabel('trm (MPa)','FontSize',18)
%  h=legend('Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma and \tau_{7}');
% %axis([0 2.5 0 30])
% %set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
% %set(gca,'YTickLabel',{'0','2','4','6','8','10'})
% set(gca,'FontName','Helvetica','FontSize',22)

