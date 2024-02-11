close all

        
n1=8;
for i=1:n1+1
    x8(i)=n1*12.5-12.5*(i-1);
end
%  n2=16;
% for i=1:n2+1
%     x16(i)=n2*6.25-6.25*(i-1);
% end
% n3=32;
% for i=1:n3+1
%     x32(i)=n3*3.125-3.125*(i-1);
% end
% n4=64;
% for i=1:n4+1
%     x64(i)=n4*1.5625-1.5625*(i-1);
% end
%  n5=128;
% for i=1:n5+1
%     x128(i)=n5*0.78125-0.78125*(i-1);
% end


phi33_fixed_Lc1_load=load('cube8_Lc1_phi33_fixed.txt');
phi33_Lc1_phi33_fixed(1)=phi33_fixed_Lc1_load(4847,11);
phi33_Lc1_phi33_fixed(2)=phi33_fixed_Lc1_load(4846,11);
k=4858;
for i=3:n1+1
    phi33_Lc1_phi33_fixed(i)=phi33_fixed_Lc1_load(k,11);
    k=k+8;
    
end


phi33_fixed_Lc01_load=load('cube8_Lc01_phi33_fixed.txt');
phi33_Lc01_phi33_fixed(1)=phi33_fixed_Lc01_load(4847,11);
phi33_Lc01_phi33_fixed(2)=phi33_fixed_Lc01_load(4846,11);
k=4858;
for i=3:n1+1
    phi33_Lc01_phi33_fixed(i)=phi33_fixed_Lc01_load(k,11);
    k=k+8;
    
end


phi33_fixed_Lc10_load=load('cube8_Lc10_phi33_fixed.txt');
phi33_Lc10_phi33_fixed(1)=phi33_fixed_Lc10_load(4847,11);
phi33_Lc10_phi33_fixed(2)=phi33_fixed_Lc10_load(4846,11);
k=4858;
for i=3:n1+1
    phi33_Lc10_phi33_fixed(i)=phi33_fixed_Lc10_load(k,11);
    k=k+8;
    
end



% % 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 

figure (1)
plot(phi33_Lc01_phi33_fixed,x8,'-sk')
hold on
plot(phi33_Lc1_phi33_fixed,x8,'-xk')
plot(phi33_Lc10_phi33_fixed,x8,'-+k')
hold off
xlabel('phi33')
ylabel('height')
h=legend('Lcptone','Lcone','Lcten','Location','SouthEast');
set(gca,'FontName','Helvetica','FontSize',16)
%set(gca,'YTickLabel',{'0','12.5','25','37.5','50','62.5','75','87.5','100'})
grid on

% 
% figure (2)
% plot(phi33_Lc01_phi33_fixed,x8,'-sr')
% hold off
% xlabel('phi33')
% ylabel('height')
% h=legend('Lc01','Location','SouthWest');
% set(gca,'FontName','Helvetica','FontSize',16)
% %set(gca,'YTickLabel',{'0','12.5','25','37.5','50','62.5','75','87.5','100'})
% grid on

% figure (3)
% plot(phi33_Lc10_phi33_fixed,x8,'-sr')
% hold off
% xlabel('phi33')
% ylabel('height')
% h=legend('Lc10','Location','SouthWest');
% set(gca,'FontName','Helvetica','FontSize',16)
% %set(gca,'YTickLabel',{'0','12.5','25','37.5','50','62.5','75','87.5','100'})
% grid on
