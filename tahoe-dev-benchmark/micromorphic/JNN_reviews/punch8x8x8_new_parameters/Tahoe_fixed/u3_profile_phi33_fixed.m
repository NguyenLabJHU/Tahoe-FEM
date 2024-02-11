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


u3_fixed_Lc1_load=load('cube8_Lc1_phi33_fixed.txt');
u3_Lc1_phi33_fixed(1)=phi33_fixed_Lc1_load(4847,14);
u3_Lc1_phi33_fixed(2)=phi33_fixed_Lc1_load(4846,14);
k=4858;
for i=3:n1+1
    u3_Lc1_phi33_fixed(i)=phi33_fixed_Lc1_load(k,14);
    k=k+8;
    
end


u3_fixed_Lc01_load=load('cube8_Lc01_phi33_fixed.txt');
u3_Lc01_phi33_fixed(1)=phi33_fixed_Lc01_load(4847,14);
u3_Lc01_phi33_fixed(2)=phi33_fixed_Lc01_load(4846,14);
k=4858;
for i=3:n1+1
    u3_Lc01_phi33_fixed(i)=phi33_fixed_Lc01_load(k,14);
    k=k+8;
    
end


u3_fixed_Lc10_load=load('cube8_Lc10_phi33_fixed.txt');
u3_Lc10_phi33_fixed(1)=phi33_fixed_Lc10_load(4847,14);
u3_Lc10_phi33_fixed(2)=phi33_fixed_Lc10_load(4846,14);
k=4858;
for i=3:n1+1
    u3_Lc10_phi33_fixed(i)=phi33_fixed_Lc10_load(k,14);
    k=k+8;
    
end



% % 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 

figure (1)
plot(u3_Lc01_phi33_fixed,x8,'-sk')
hold on
plot(u3_Lc1_phi33_fixed,x8,'-xk')
plot(u3_Lc10_phi33_fixed,x8,'-+k')
hold off
xlabel('phi33')
ylabel('height')
h=legend('Lcptone','Lcone','Lcten','Location','SouthWest');
set(gca,'FontName','Helvetica','FontSize',16)
%set(gca,'YTickLabel',{'0','12.5','25','37.5','50','62.5','75','87.5','100'})
grid on

