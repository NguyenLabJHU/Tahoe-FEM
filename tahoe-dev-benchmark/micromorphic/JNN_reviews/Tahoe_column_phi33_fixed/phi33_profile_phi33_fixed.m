        
n1=8;
for i=1:n1+1
    x8(i)=n1*12.5-12.5*(i-1);
end
 n2=16;
for i=1:n2+1
    x16(i)=n2*6.25-6.25*(i-1);
end
n3=32;
for i=1:n3+1
    x32(i)=n3*3.125-3.125*(i-1);
end
n4=64;
for i=1:n4+1
    x64(i)=n4*1.5625-1.5625*(i-1);
end
 n5=128;
for i=1:n5+1
    x128(i)=n5*0.78125-0.78125*(i-1);
end


simcol8_phi33_fixed=load('column8Lc1_phi33_fixed.txt');
phi33_col8_phi33_fixed(1)=simcol8_phi33_fixed(8,11);
phi33_col8_phi33_fixed(2)=simcol8_phi33_fixed(7,11);
k=31;
for i=3:n1+1
    phi33_col8_phi33_fixed(i)=simcol8_phi33_fixed(k,11);
    k=k+18;
    
end
% 
simcol16_phi33_fixed=load('column16Lc1_phi33_fixed.txt');
phi33_col16_phi33_fixed(1)=simcol16_phi33_fixed(8,11);
phi33_col16_phi33_fixed(2)=simcol16_phi33_fixed(7,11);

k=31;
for i=3:n2+1
    phi33_col16_phi33_fixed(i)=simcol16_phi33_fixed(k,11);
    k=k+18;
    
end
% % 
simcol32_phi33_fixed=load('column32Lc1_phi33_fixed.txt');
phi33_col32_phi33_fixed(1)=simcol32_phi33_fixed(8,11);
phi33_col32_phi33_fixed(2)=simcol32_phi33_fixed(7,11);

k=31;
for i=3:n3+1
    phi33_col32_phi33_fixed(i)=simcol32_phi33_fixed(k,11);
    k=k+18;
    
end
% % 

simcol64_phi33_fixed=load('column64Lc1_phi33_fixed.txt');
phi33_col64_phi33_fixed(1)=simcol64_phi33_fixed(8,11);
phi33_col64_phi33_fixed(2)=simcol64_phi33_fixed(7,11);

k=31;
for i=3:n4+1
    phi33_col64_phi33_fixed(i)=simcol64_phi33_fixed(k,11);
    k=k+18;
    
end

simcol128_phi33_fixed=load('column128Lc1_phi33_fixed.txt');
phi33_col128_phi33_fixed(1)=simcol128_phi33_fixed(8,11);
phi33_col128_phi33_fixed(2)=simcol128_phi33_fixed(7,11);

k=31;
for i=3:n5+1
    phi33_col128_phi33_fixed(i)=simcol128_phi33_fixed(k,11);
    k=k+18;
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


simcol8=load('column8Lc1.txt');
phi33_col8(1)=simcol8(8,11);
phi33_col8(2)=simcol8(7,11);

k=31;
for i=3:n1+1
   phi33_col8(i)=simcol8(k,11);
    k=k+18;
    
end

simcol16=load('column16Lc1.txt');
phi33_col16(1)=simcol16(8,11);
phi33_col16(2)=simcol16(7,11);

k=31;
for i=3:n2+1
    phi33_col16(i)=simcol16(k,11);
    k=k+18;
    
end

simcol32=load('column32Lc1.txt');
phi33_col32(1)=simcol32(8,11);
phi33_col32(2)=simcol32(7,11);
k=31;
for i=3:n3+1
    phi33_col32(i)=simcol32(k,11);
    k=k+18;
    
end


simcol64=load('column64Lc1.txt');
phi33_col64(1)=simcol64(8,11);
phi33_col64(2)=simcol64(7,11);
k=31;
for i=3:n4+1
    phi33_col64(i)=simcol64(k,11);
    k=k+18;
    
end

simcol128=load('column128Lc1.txt');
phi33_col128(1)=simcol128(8,11);
phi33_col128(2)=simcol128(7,11);
k=31;%2154;%2190;
for i=3:n5+1
    phi33_col128(i)=simcol128(k,11);
    k=k+18;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 

figure (1)
plot(phi33_col8_phi33_fixed,x8,'-+k')
hold on
plot(phi33_col16_phi33_fixed,x16,'-ok')
plot(phi33_col32_phi33_fixed,x32,'-s')
plot(phi33_col64_phi33_fixed,x64,'-ob')
plot(phi33_col128_phi33_fixed,x128,'-g')
hold off
xlabel('phi33')
ylabel('heigth')
h=legend('8 el','16 el','32 el','64 el','128 el','Location','SouthWest');
grid on


% 
figure (2)
plot(phi33_col128,x128,'-k')
hold on
plot(phi33_col128_phi33_fixed,x128,'-r')
hold off
xlabel('phi33')
ylabel('heigth')
h=legend('128 elelements','128 elements \phi_{33} fixed at the top surface','Location','SouthWest');
grid on

figure (3)
plot(phi33_col8,x8,'-+k')
hold on
plot(phi33_col16,x16,'-ok')
plot(phi33_col32,x32,'-s')
plot(phi33_col64,x64,'-ob')
plot(phi33_col128,x128,'-g')
hold off
xlabel('phi33')
ylabel('heigth')
h=legend('8 el','16 el','32 el','64 el','128 el','Location','SouthWest');
grid on
