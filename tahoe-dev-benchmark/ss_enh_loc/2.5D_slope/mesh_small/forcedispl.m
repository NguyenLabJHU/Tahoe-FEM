%
sim=load('slope.io2_R.txt');
d=abs(sim(:,1));
f=abs(sim(:,2));
%
figure(1)
plot(d,f,'LineWidth',1)
xlabel('DISPLACEMENT (mm)')
%set(get(gco,'XLabel'),'DISPLACEMENT (mm)','FontName','Helvetica','FontSize',16)
ylabel('FORCE (MN)')
%legend('node 27')
set(gca,'FontName','Helvetica','FontSize',16)
%