
close all
clear all

format short e

simFSE=load('el_cube99FSE_721.txt');
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

simkappanu=load('el_kappannu_721.txt');
timekappanu=simkappanu(:,1);
s11kappanu=simkappanu(:,578);
s22kappanu=simkappanu(:,579);
s33kappanu=simkappanu(:,580);
s12kappanu=simkappanu(:,581);
s13kappanu=simkappanu(:,582);
s21kappanu=simkappanu(:,583);
s23kappanu=simkappanu(:,584);
s31kappanu=simkappanu(:,585);
s32kappanu=simkappanu(:,586);
e11kappanu=simkappanu(:,587);
e22kappanu=simkappanu(:,588);
e33kappanu=simkappanu(:,589);
e12kappanu=simkappanu(:,590);
e13kappanu=simkappanu(:,591);
e21kappanu=simkappanu(:,592);
e23kappanu=simkappanu(:,593);
e31kappanu=simkappanu(:,594);
e32kappanu=simkappanu(:,595);
s_invkappanu=simkappanu(:,596);
rel_invkappanu=simkappanu(:,597);
hos_invkappanu=simkappanu(:,598);
trskappanu=simkappanu(:,599);
trrelkappanu=simkappanu(:,600);
trmkappanu=simkappanu(:,601);


simetakappanu=load('el_cube99etakappanu_721.txt');
timeetakappanu=simetakappanu(:,1);
s11etakappanu=simetakappanu(:,578);
s22etakappanu=simetakappanu(:,579);
s33etakappanu=simetakappanu(:,580);
s12etakappanu=simetakappanu(:,581);
s13etakappanu=simetakappanu(:,582);
s21etakappanu=simetakappanu(:,583);
s23etakappanu=simetakappanu(:,584);
s31etakappanu=simetakappanu(:,585);
s32etakappanu=simetakappanu(:,586);
e11etakappanu=simetakappanu(:,587);
e22etakappanu=simetakappanu(:,588);
e33etakappanu=simetakappanu(:,589);
e12etakappanu=simetakappanu(:,590);
e13etakappanu=simetakappanu(:,591);
e21etakappanu=simetakappanu(:,592);
e23etakappanu=simetakappanu(:,593);
e31etakappanu=simetakappanu(:,594);
e32etakappanu=simetakappanu(:,595);
s_invetakappanu=simetakappanu(:,596);
rel_invetakappanu=simetakappanu(:,597);
hos_invetakappanu=simetakappanu(:,598);
trsetakappanu=simetakappanu(:,599);
trreletakappanu=simetakappanu(:,600);
trmetakappanu=simetakappanu(:,601);

simnoTau7=load('el_cube99noTau7_721.txt');
timenoTau7=simnoTau7(:,1);
s11noTau7=simnoTau7(:,578);
s22noTau7=simnoTau7(:,579);
s33noTau7=simnoTau7(:,580);
s12noTau7=simnoTau7(:,581);
s13noTau7=simnoTau7(:,582);
s21noTau7=simnoTau7(:,583);
s23noTau7=simnoTau7(:,584);
s31noTau7=simnoTau7(:,585);
s32noTau7=simnoTau7(:,586);
e11noTau7=simnoTau7(:,587);
e22noTau7=simnoTau7(:,588);
e33noTau7=simnoTau7(:,589);
e12noTau7=simnoTau7(:,590);
e13noTau7=simnoTau7(:,591);
e21noTau7=simnoTau7(:,592);
e23noTau7=simnoTau7(:,593);
e31noTau7=simnoTau7(:,594);
e32noTau7=simnoTau7(:,595);
s_invnoTau7=simnoTau7(:,596);
rel_invnoTau7=simnoTau7(:,597);
hos_invnoTau7=simnoTau7(:,598);
trsnoTau7=simnoTau7(:,599);
trrelnoTau7=simnoTau7(:,600);
trmnoTau7=simnoTau7(:,601);

simall=load('el_cube99M_721.txt');
timeall=simall(:,1);
s11all=simall(:,578);
s22all=simall(:,579);
s33all=simall(:,580);
s12all=simall(:,581);
s13all=simall(:,582);
s21all=simall(:,583);
s23all=simall(:,584);
s31all=simall(:,585);
s32all=simall(:,586);
e11all=simall(:,587);
e22all=simall(:,588);
e33all=simall(:,589);
e12all=simall(:,590);
e13all=simall(:,591);
e21all=simall(:,592);
e23all=simall(:,593);
e31all=simall(:,594);
e32all=simall(:,595);
s_invall=simall(:,596);
rel_invall=simall(:,597);
hos_invall=simall(:,598);
trsall=simall(:,599);
trrelall=simall(:,600);
trmall=simall(:,601);


% sim33=load('el_cube33M_25.txt');
% time33=sim33(:,1);
% trmall33=sim33(:,601);

figure(1)
plot(timeFSE,s_invFSE*1e-6,'-k','LineWidth',1, 'MarkerSize',14)
hold on
plot(timekappanu,s_invkappanu*1e-6,'-ok','LineWidth',1, 'MarkerSize',14)
plot(timeetakappanu,s_invetakappanu*1e-6,'-vk','LineWidth',1, 'MarkerSize',14)
plot(timenoTau7,s_invnoTau7*1e-6,'-^k','LineWidth',1, 'MarkerSize',14)
plot(timeall,s_invall*1e-6,'-+k','LineWidth',1, 'MarkerSize',14)
hold off

 xlabel('time','FontSize',18)
 ylabel('||devsigma||  (MPa)','FontSize',18)
 h=legend('Finite strain elasticity (FSE)','Micromorphic FSE with \kappa,\nu','Micromorphic FSE with \eta,\kappa,\nu','Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma','Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma and \tau_{7}');
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',22)

figure(2)
plot(timekappanu,rel_invkappanu*1e-6,'-ok','LineWidth',1, 'MarkerSize',14)
hold on
plot(timeetakappanu,rel_invetakappanu*1e-6,'-vk','LineWidth',1, 'MarkerSize',14)
plot(timenoTau7,rel_invnoTau7*1e-6,'-^k','LineWidth',1, 'MarkerSize',14)
plot(timeall,rel_invall*1e-6,'-+k','LineWidth',1, 'MarkerSize',14)
hold off
xlabel('time','FontSize',18)
ylabel('||dev(s-sigma)||    (MPa)','FontSize',18)
h=legend('Micromorphic FSE with \kappa,\nu','Micromorphic FSE with \eta,\kappa,\nu','Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma','Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma and \tau_{7}');
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',22)


figure(3)
plot(timeall,hos_invall*1e-6,'-+k','LineWidth',1, 'MarkerSize',14)
xlabel('time','FontSize',18)
ylabel('||dev(m)||    (MPa)','FontSize',18)
 h=legend('Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma and \tau_{7}');
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',22)


figure(4)
plot(timeFSE,trsFSE*1e-6,'-k','LineWidth',1, 'MarkerSize',14)
hold on
plot(timekappanu,trskappanu*1e-6,'-ok','LineWidth',1, 'MarkerSize',14)
plot(timeetakappanu,trsetakappanu*1e-6,'-vk','LineWidth',1, 'MarkerSize',14)
plot(timenoTau7,trsnoTau7*1e-6,'-^k','LineWidth',1, 'MarkerSize',14)
plot(timeall,trsall*1e-6,'-+k','LineWidth',1, 'MarkerSize',14)
hold off

 xlabel('time','FontSize',18)
 ylabel('trsigma (MPa)','FontSize',18)
 h=legend('Finite strain elasticity (FSE)','Micromorphic FSE with \kappa,\nu','Micromorphic FSE with \eta,\kappa,\nu','Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma','Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma and \tau_{7}');
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
 h=legend('Micromorphic FSE with \kappa,\nu','Micromorphic FSE with \eta,\kappa,\nu','Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma','Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma and \tau_{7}');
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',22)

figure(6)
plot(timeall,trmall*1e-6,'-+k','LineWidth',1, 'MarkerSize',14)
 xlabel('time','FontSize',18)
 ylabel('trm (MPa)','FontSize',18)
 h=legend('Micromorphic FSE with \eta,\kappa,\nu,\tau,\sigma and \tau_{7}');
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',22)

