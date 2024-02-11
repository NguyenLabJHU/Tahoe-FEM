
%close all
clear all

format short e

sim33=load('cube3M_force.txt');
sim33FSE=load('cube3FSE_force.txt');
sim66=load('cube6M_force.txt');
sim66FSE=load('cube6FSE_force.txt');
sim99=load('cube9M_force.txt');
sim99FSE=load('cube9FSE_force.txt');

sim12M=load('cube12M_force.txt');
sim12FSE=load('cube12FSE_force.txt');

Displ333=sim33(:,1);
Force333=sim33(:,2);

Displ3FSE=sim33FSE(:,1);
Force3FSE=sim33FSE(:,2);


Displ666=sim66(:,1);
Force666=sim66(:,2);


Displ6FSE=sim66FSE(:,1);
Force6FSE=sim66FSE(:,2);


Displ999=sim99(:,1);
Force999=sim99(:,2);

Displ9FSE=sim99FSE(:,1);
Force9FSE=sim99FSE(:,2);


Displ12M=sim12M(:,1);
Force12M=sim12M(:,2);

Displ12FSE=sim12FSE(:,1);
Force12FSE=sim12FSE(:,2);

%plot(-Displ333,Force333,'-^k','LineWidth',1, 'MarkerSize',8)
%hold on
%plot(-Displ666,Force666,'-ok','LineWidth',1, 'MarkerSize',8)
%plot(-Displ999,-Force999,'-+k','LineWidth',1, 'MarkerSize',8)
%hold off
%h=legend('MFSE 3x3x3','MFSE 6x6x6','Simo material 9x9x9');
% 
% sim111=load('cube111force.txt');
% Displ111=sim111(:,1);
% Force111=sim111(:,2);
% 
% sim222=load('cube222force.txt');
% Displ222=sim222(:,1); 
% Force222=sim222(:,2);
figure(1)
plot(Displ333,Force333,'-ok',Displ666,Force666,'-+k',Displ999,Force999,'-^k',Displ12M,Force12M,'-sk')
legend('3x3x3 Microm.','6x6x6 Microm.','9x9x9 Microm.','12x12x12 Microm.');
set(gca,'FontName','Helvetica','FontSize',16)


figure(2)

 plot(Displ3FSE,Force3FSE,'-^k','LineWidth',1, 'MarkerSize',8)
 hold on
 plot(Displ6FSE,Force6FSE,'-ok','LineWidth',1, 'MarkerSize',8)
 plot(Displ9FSE,Force9FSE,'-+k','LineWidth',1, 'MarkerSize',8)
 plot(Displ12FSE,Force12FSE,'-sk','LineWidth',1, 'MarkerSize',8)
 hold off
h=legend('3x3x3 FSE',' 6x6x6 FSE','9x9x9 FSE','12x12x12 FSE');
