% This M-file calculates the total system energy as a function of time
% Takes Tahoe EnSight output and reads it in
% DEF 19 Mar 2005
% DEF 21 Mar 2005 - modified to take 2 sets, plot them together
% DEF 29 Mar 2005 - modified to do 2D
% DEF 27 Aug 2005 - modified to do 2d benchmark comparison


clear;
close all;

% now begin time loop (requires input)
maxnumstep1 = 99;       % highest output number to load (0-N)
maxnumstep2 = 99;       % highest output number to load (0-N)

% initialize some stuff
total_E_1 = zeros(maxnumstep1+1,1);  
total_E_2 = zeros(maxnumstep2+1,1);

% begin time loop
for i=1:maxnumstep1+1
    % set up the filenames - assumes form: atom.3D.0.gp1.ps%04i.KE/.PE
    KE_file1 = sprintf('atom.2D.0.gp1.ps%04i.KE',i-1);
    PE_file1 = sprintf('atom.2D.0.gp1.ps%04i.PE',i-1);
    
    % now read in the data
    KE1 = dlmread(char(KE_file1),' ',[352,0,539,0]);
    PE1 = dlmread(char(PE_file1),' ',[352,0,539,0]);
    
    % System kinetic and potential energy
    total_PE1(i) = sum(PE1);
    total_KE1(i) = sum(KE1);
    
    % Now determine total system energy as sum of the above
    total_E_1(i) = sum(KE1+PE1);
end

% begin time loop
for i=1:maxnumstep2+1
    % set up the filenames - assumes form: atom.3D.0.gp1.ps%04i.KE/.PE
    KE_file2 = sprintf('/Users/davidfarrell/Code/Tahoe/Tahoe_David_5/development_benchmark_XML/bridging/2D_simple_wave/atom.2D.0.gp1.ps%04i.KE',i-1);
    PE_file2 = sprintf('/Users/davidfarrell/Code/Tahoe/Tahoe_David_5/development_benchmark_XML/bridging/2D_simple_wave/atom.2D.0.gp1.ps%04i.PE',i-1);
    
    % now read in the data
    KE2 = dlmread(char(KE_file2),' ',[19,0,206,0]);
    PE2 = dlmread(char(PE_file2),' ',[19,0,206,0]);
    
    % System kinetic and potential energy
    total_PE2(i) = sum(PE2);
    total_KE2(i) = sum(KE2);
    
    % Now determine total system energy as sum of the above
    total_E_2(i) = sum(KE2+PE2);
end

% figure out dissipation percentage:
dissp_percent = 1 - abs((total_E_2(maxnumstep2+1) - total_E_1(maxnumstep1+1))/(total_E_1(1)-total_E_1(maxnumstep1+1)))

% now plot some stuff
figure
    hold on;
    plot([1:maxnumstep1+1],total_E_1,'k-',[1:maxnumstep2+1],total_E_2,'r-')
    title('1NN 2D Simple Wave LJ Hexagonal Model');  
    ylabel('System Energy + Applied Work');
    xlabel('Output Step');
    hold off;
    
figure
    hold on;
    plot([1:maxnumstep1+1],total_PE1,'k-',[1:maxnumstep2+1],total_PE2,'r-');
    title('1NN 2D Simple Wave LJ Hexagonal Model'); 
    ylabel('System Potential Energy');
    xlabel('Output Step');
    hold off;    
 
figure
    hold on;
    plot([1:maxnumstep1+1],total_KE1,'k-',[1:maxnumstep2+1],total_KE2,'r-');
    title('1NN 2D Simple Wave LJ Hexagonal Model'); 
    ylabel('System Kinetic Energy');
    xlabel('Output Step');
    hold off;     