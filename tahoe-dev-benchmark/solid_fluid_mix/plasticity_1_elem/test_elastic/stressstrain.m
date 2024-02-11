
%close all
clear all

%"test" is the Simo elastic implementation
%"elastic" is the ln elastic implementation in solid fluid mix
%"elasticlin" is the linear elastic implementation in solid fluid mix

format short e

sim=load('nddata_test_27.txt');
time=sim(:,1);
s11test=(sim(:,2))/1000;
s22test=(sim(:,3))/1000;
s33test=(sim(:,4))/1000;
s23test=(sim(:,5))/1000;
s13test=(sim(:,6))/1000;
s12test=(sim(:,7))/1000;
e11test=sim(:,8);
e22test=sim(:,9);
e33test=sim(:,10);
e23test=sim(:,11);
e13test=sim(:,12);
e12test=sim(:,13);

sim=load('nddata_elastic_27.txt');
time=sim(:,1);
s11elas=(sim(:,2))/1000;
s22elas=(sim(:,3))/1000;
s33elas=(sim(:,4))/1000;
s23elas=(sim(:,5))/1000;
s13elas=(sim(:,6))/1000;
s12elas=(sim(:,7))/1000;
e11elas=sim(:,8);
e22elas=sim(:,9);
e33elas=sim(:,10);
e23elas=sim(:,11);
e13elas=sim(:,12);
e12elas=sim(:,13);

sim=load('nddata_elasticlin_1.txt');
time=sim(:,1);
s11elaslin=(sim(:,2))/1000;
s22elaslin=(sim(:,3))/1000;
s33elaslin=(sim(:,4))/1000;
s23elaslin=(sim(:,5))/1000;
s13elaslin=(sim(:,6))/1000;
s12elaslin=(sim(:,7))/1000;
e11elaslin=sim(:,8);
e22elaslin=sim(:,9);
e33elaslin=sim(:,10);
e23elaslin=sim(:,11);
e13elaslin=sim(:,12);
e12elaslin=sim(:,13);


% stress strain
figure(1)
plot(abs(e33test),abs(s33test),'-k',abs(e33elas),abs(s33elas),'--k',abs(e33elaslin),abs(s33elaslin),'ok','LineWidth',2)
xlabel('e33')
ylabel('s33')
legend('test - Simo elastic','elastic - ln','elastic - lin')
set(gca,'FontName','Helvetica','FontSize',16)
%set(gca,'XTickLabel',{'0','1','2','3','4','5'})

