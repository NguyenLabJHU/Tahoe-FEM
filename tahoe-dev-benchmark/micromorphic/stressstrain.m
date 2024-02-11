
close all
clear all

format short e

%sim=load('el_data_1.txt');
%sim=load('el_data_file_1.txt');
%sim=load('el_data_rot_1.txt');
%sim=load('el_data_newBC_1.txt');
%sim=load('el_iso_data_1.txt');
%sim=load('el_data_etaandkappa_1.txt');
%sim=load('el_data_eta_1.txt');
%sim=load('el_datanu_1.txt');
sim=load('el_FSE_data_1.txt');
time=sim(:,1);
s11=(sim(:,2));
s22=(sim(:,3));
s33=(sim(:,4));
s23=(sim(:,5));
s13=(sim(:,6));
s12=(sim(:,7));
e11=sim(:,8);
e22=sim(:,9);
e33=sim(:,10);
e23=sim(:,11);
e13=sim(:,12);
e12=sim(:,13);
% thing1=sim(:,14);
% thing2=sim(:,15);
% J=sim(:,16);

% Principal Stresses
for i=1:length(time)
    sig = [s11(i) s12(i) s13(i) ;
           s12(i) s22(i) s23(i) ;
           s13(i) s23(i) s33(i) ];
    lambda = eig(sig);
    s1(i,1) = max(lambda);
    s2(i,1) = lambda(2);
    s3(i,1) = min(lambda);
end

I1=s1+s2+s3;
meanstress=I1./3;

% Principal Strains
for i=1:length(time)
    eps = [e11(i) e12(i) e13(i) ;
           e12(i) e22(i) e23(i) ;
           e13(i) e23(i) e33(i) ];
    lambda = eig(eps);
    e1(i,1) = max(lambda);
    e2(i,1) = lambda(2);
    e3(i,1) = min(lambda);
end

% stress strain
figure(1)
plot(e33,s33,'--rs','LineWidth',2)
%plot(abs(e33),abs(s33),'-k','LineWidth',2)
%plot(abs(e11),abs(s11),'--k',abs(e33),abs(s33),'-k','LineWidth',2)
xlabel('e33')
ylabel('s33')
%ylabel('s33')
%legend('s11-e11','s33-e33')
set(gca,'FontName','Helvetica','FontSize',16)
%set(gca,'XTickLabel',{'0','1','2','3','4','5'})
% 
%principal stress strain
% figure(2)
% %plot(abs(e3),abs(s3),'LineWidth',2)
% plot(e3,s3,'--rs','LineWidth',2)
% xlabel('e3')
% ylabel('s3')
% % set(gca,'FontName','Helvetica','FontSize',16)
% 
% % principal stress strain
% figure(3)
% %plot(abs(e3),abs(s3),'LineWidth',2)
% plot(e2,s2,'--rs','LineWidth',2)
% xlabel('e2')
% ylabel('s2')
% set(gca,'FontName','Helvetica','FontSize',16)
% %principal stress strain
% 
% figure(4)
% %plot(abs(e3),abs(s3),'LineWidth',2)
% plot(e22,s22,'--rs','LineWidth',2)
% xlabel('e22')
% ylabel('s22')
% 
% figure(5)
% %plot(abs(e3),abs(s3),'LineWidth',2)
% plot(e11,s11,'--rs','LineWidth',2)
% xlabel('e11')
% ylabel('s11')
% 
% 
% figure(6)
% %plot(abs(e3),abs(s3),'LineWidth',2)
% plot(e2,s1,'--rs','LineWidth',2)
% xlabel('e2')
% ylabel('s1')

%set(gca,'FontName','Helvetica','FontSize',16)
% 
% 
% % 
% % 
% % % principal stress strain
% % figure(4)
% % %plot(abs(e3),abs(s3),'LineWidth',2)
% % plot(e11,s11,'--rs','LineWidth',2)
% % xlabel('e11')
% % ylabel('s11')
% % set(gca,'FontName','Helvetica','FontSize',16)
% % 
% % % principal stress strain
% % figure(5)
% % %plot(abs(e3),abs(s3),'LineWidth',2)
% % plot(e22,s22,'--rs','LineWidth',2)
% % xlabel('e22')
% % ylabel('s22')
% % set(gca,'FontName','Helvetica','FontSize',16)
% % % principal stress strain
% 
% 
% 
% % % principal stress strain
% % figure(6)
% % %plot(abs(e3),abs(s3),'LineWidth',2)
% % plot(e13,s13,'--rs','LineWidth',2)
% % xlabel('e13')
% % ylabel('s13')
% % set(gca,'FontName','Helvetica','FontSize',16)
% % % principal stress strain




