        
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
u3_col8_phi33_fixed(1)=simcol8_phi33_fixed(8,14);
u3_col8_phi33_fixed(2)=simcol8_phi33_fixed(7,14);
k=31;
for i=3:n1+1
    u3_col8_phi33_fixed(i)=simcol8_phi33_fixed(k,14);
    k=k+18;
    
end
% 
simcol16_phi33_fixed=load('column16Lc1_phi33_fixed.txt');
u3_col16_phi33_fixed(1)=simcol16_phi33_fixed(8,14);
u3_col16_phi33_fixed(2)=simcol16_phi33_fixed(7,14);

k=31;
for i=3:n2+1
    u3_col16_phi33_fixed(i)=simcol16_phi33_fixed(k,14);
    k=k+18;
    
end
% % 
simcol32_phi33_fixed=load('column32Lc1_phi33_fixed.txt');
u3_col32_phi33_fixed(1)=simcol32_phi33_fixed(8,14);
u3_col32_phi33_fixed(2)=simcol32_phi33_fixed(7,14);

k=31;
for i=3:n3+1
    u3_col32_phi33_fixed(i)=simcol32_phi33_fixed(k,14);
    k=k+18;
    
end
% % 

simcol64_phi33_fixed=load('column64Lc1_phi33_fixed.txt');
u3_col64_phi33_fixed(1)=simcol64_phi33_fixed(8,14);
u3_col64_phi33_fixed(2)=simcol64_phi33_fixed(7,14);

k=31;
for i=3:n4+1
    u3_col64_phi33_fixed(i)=simcol64_phi33_fixed(k,14);
    k=k+18;
    
end

simcol128_phi33_fixed=load('column128Lc1_phi33_fixed.txt');
u3_col128_phi33_fixed(1)=simcol128_phi33_fixed(8,14);
u3_col128_phi33_fixed(2)=simcol128_phi33_fixed(7,14);

k=31;
for i=3:n5+1
    u3_col128_phi33_fixed(i)=simcol128_phi33_fixed(k,14);
    k=k+18;
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


simcol8=load('column8Lc1.txt');
u3_col8(1)=simcol8(8,14);
u3_col8(2)=simcol8(7,14);

k=31;
for i=3:n1+1
   u3_col8(i)=simcol8(k,14);
    k=k+18;
    
end

simcol16=load('column16Lc1.txt');
u3_col16(1)=simcol16(8,14);
u3_col16(2)=simcol16(7,14);

k=31;
for i=3:n2+1
    u3_col16(i)=simcol16(k,14);
    k=k+18;
    
end

simcol32=load('column32Lc1.txt');
u3_col32(1)=simcol32(8,14);
u3_col32(2)=simcol32(7,14);
k=31;
for i=3:n3+1
    u3_col32(i)=simcol32(k,14);
    k=k+18;
    
end


simcol64=load('column64Lc1.txt');
u3_col64(1)=simcol64(8,14);
u3_col64(2)=simcol64(7,14);
k=31;
for i=3:n4+1
    u3_col64(i)=simcol64(k,14);
    k=k+18;
    
end

simcol128=load('column128Lc1.txt');
u3_col128(1)=simcol128(8,14);
u3_col128(2)=simcol128(7,14);
k=31;%2154;%2190;
for i=3:n5+1
    u3_col128(i)=simcol128(k,14);
    k=k+18;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 

figure (1)
plot(u3_col8_phi33_fixed,x8,'-+k')
hold on
plot(u3_col16_phi33_fixed,x16,'-ok')
plot(u3_col32_phi33_fixed,x32,'-s')
plot(u3_col64_phi33_fixed,x64,'-ob')
plot(u3_col128_phi33_fixed,x128,'-g')
hold off
xlabel('u3')
ylabel('heigth')
h=legend('8 el','16 el','32 el','64 el','128 el','Location','SouthWest');
grid on


% 
figure (2)
plot(u3_col128,x128,'-k')
hold on
plot(u3_col128_phi33_fixed,x128,'-r')
hold off
xlabel('u3')
ylabel('heigth')
h=legend('128 elelements','128 elements \phi_{33} fixed at the top surface','Location','SouthWest');
grid on

figure (3)
plot(u3_col8,x8,'-+k')
hold on
plot(u3_col16,x16,'-ok')
plot(u3_col32,x32,'-s')
plot(u3_col64,x64,'-ob')
plot(u3_col128,x128,'-g')
hold off
xlabel('u3')
ylabel('heigth')
h=legend('8 el','16 el','32 el','64 el','128 el','Location','SouthWest');
grid on
