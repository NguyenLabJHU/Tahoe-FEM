clear all
%
sim=load('nddata_0699.txt');
%sim=load('nddata_1036.txt');
%sim=load('nddata_2470.txt');
%sim=load('nddata_2479.txt');
time=(sim(:,1));
%nd
Tt_1=(sim(:,2));
Tn_1=(sim(:,3));
tension_1=(sim(:,4));
cohesion_1=(sim(:,5));
friction_1=(sim(:,6));
dilation_1=(sim(:,7));
%
%plot yield surface
Tn = -30:0.01:0.55;
%initial
c=4; %MPa
chi=0.5; %MPa
phi=0.4; %radians 
Fyieldfunc0 = (c - Tn*tan(phi)).^2 - (c - chi*tan(phi))^2;
Fyield0 = 0.5*(Fyieldfunc0+abs(Fyieldfunc0));
Tt0 = sqrt(Fyield0);
Tt0neg = -sqrt(Fyield0);
%final at ip 1
c=cohesion_1(length(cohesion_1));
chi=tension_1(length(cohesion_1));
tanphi=friction_1(length(cohesion_1)); 
Fyieldfunc = (c - Tn*tan(phi)).^2 - (c - chi*tan(phi))^2;
Fyield = 0.5*(Fyieldfunc+abs(Fyieldfunc));
Tt = sqrt(Fyield);
Ttneg = -sqrt(Fyield);
%
figure(1)
plot(Tn_1,Tt_1,'-.ok',Tn,Tt0,'--k',Tn,Tt,'-k',Tn,Tt0neg,'--k',Tn,Ttneg,'-k','LineWidth',1)
xlabel('Tn (MPa)')
ylabel('Tt (MPa)')
legend('traction path','initial failure surface','final failure surface')
grid on
%axis([0 2.5 0 30])
%set(gca,'XTickLabel',{'0','0.5','1','1.5','2','2.5','3'})
%set(gca,'YTickLabel',{'0','2','4','6','8','10'})
set(gca,'FontName','Helvetica','FontSize',16)

