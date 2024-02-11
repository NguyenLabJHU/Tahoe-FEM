clear all
%
%sim=load('eldataloc_1_cnstrsh_el1.txt');
sim=load('eldataloc_1.txt');
sig12_weak_disc=(sim(:,7));
eps12_weak_disc=(sim(:,13));
hardmod=(sim(:,25));
chihardmod=(sim(:,26));
%
figure(1)
plot(eps12_weak_disc,sig12_weak_disc,'-k','LineWidth',2)
xlabel('SHEAR STRAIN')
ylabel('SHEAR STRESS (MPa)')
%legend('constrained shear')
set(gca,'FontName','Helvetica','FontSize',16)
set(gca,'XTickLabel',{'0','0.0005','0.001','0.0015','0.002','0.0025'})
%
figure(2)
plot(eps12_weak_disc,hardmod,'-k','LineWidth',2)
xlabel('SHEAR STRAIN')
ylabel('HARDENING MODULUS (MPa^3)')
%legend('constrained shear')
set(gca,'FontName','Helvetica','FontSize',16)
set(gca,'XTickLabel',{'0','0.0005','0.001','0.0015','0.002','0.0025'})
%
figure(3)
plot(eps12_weak_disc,chihardmod,'-k','LineWidth',2)
xlabel('SHEAR STRAIN')
ylabel('Chi Hardening Modulus (MPa^3)')
%legend('constrained shear')
set(gca,'FontName','Helvetica','FontSize',16)
set(gca,'XTickLabel',{'0','0.0005','0.001','0.0015','0.002','0.0025'})
%
figure(4)
[AX,H1,H2] = plotyy(eps12_weak_disc,sig12_weak_disc,eps12_weak_disc,hardmod)
ax1=AX(1);
ax2=AX(2);
set(H1,'LineStyle','-','LineWidth',1.5,'Color','black')
set(H2,'LineStyle','--','LineWidth',1.5,'Color','black')
xlabel('SHEAR STRAIN')
set(get(AX(1),'Ylabel'),'String','SHEAR STRESS (MPa)','Color','black','FontName','Helvetica','FontSize',16)
set(ax1,'XColor','k','YColor','k')
set(ax1,'YTickLabel',{'0','10','20'})
set(get(AX(2),'Ylabel'),'String','HARDENING MODULUS (MPa^3)','Color','black','FontName','Helvetica','FontSize',16)
set(ax2,'XColor','k','YColor','k')
set(ax2,'YTickLabel',{'0','2e5','4e5'})
set(ax1,'FontName','Helvetica','FontSize',16)
set(ax2,'FontName','Helvetica','FontSize',16)
set(ax1,'XTickLabel',{'0','0.0005','0.001','0.0015','0.002','0.0025'})
set(ax2,'XTickLabel',{'0','0.0005','0.001','0.0015','0.002','0.0025'})
%
figure(5)
[AX,H1,H2] = plotyy(eps12_weak_disc,sig12_weak_disc,eps12_weak_disc,chihardmod)
ax1=AX(1);
ax2=AX(2);
set(H1,'LineStyle','-','LineWidth',1.5,'Color','black')
set(H2,'LineStyle','--','LineWidth',1.5,'Color','black')
xlabel('SHEAR STRAIN')
set(get(AX(1),'Ylabel'),'String','SHEAR STRESS (MPa)','Color','black','FontName','Helvetica','FontSize',16)
set(ax1,'XColor','k','YColor','k')
%set(ax1,'YTickLabel',{'0','10','20'})
set(get(AX(2),'Ylabel'),'String','CHI HARDENING MODULUS, (MPa^3)','Color','black','FontName','Helvetica','FontSize',16)
set(ax2,'XColor','k','YColor','k')
set(ax2,'YTickLabel',{'0','1e6','2e6','3e6','4e6','5e6','6e6'})
set(ax1,'FontName','Helvetica','FontSize',16)
set(ax2,'FontName','Helvetica','FontSize',16)
set(ax1,'XTickLabel',{'0','0.0005','0.001','0.0015','0.002','0.0025'})
set(ax2,'XTickLabel',{'0','0.0005','0.001','0.0015','0.002','0.0025'})
%