% used to plot average velocity of the poly result and its corresponding
% ellip result together for comparison

clc
clear
close all
timestep = 1.0e-07;

result_same_mineral = load('dep_progress_eccentric_same_mineral');   
result_diff_mineral = load('dep_progress_eccentric_different_mineral');  

%interval_print = 10;
interval_plot= 500;



timestep_same = result_same_mineral(:,1);
timestep_diff = result_diff_mineral(:,1);

kinetic_same = result_same_mineral(:,13);
kinetic_diff = result_diff_mineral(:,13);

total_same = result_same_mineral(:,15);
total_diff = result_diff_mineral(:,15);


figure(1)
hold on
plot(timestep_same*timestep, kinetic_same, '-b','LineWidth',2)
plot(timestep_diff*timestep, kinetic_diff, '-r','LineWidth',2)

xlabel('time (second)')
ylabel('kinetic energy (Pa)')
legend('two particles with same mineral', 'two particles with different mineral')


figure(2)
hold on
plot(timestep_same*timestep, total_same, '-b','LineWidth',2)
plot(timestep_diff*timestep, total_diff, '-r','LineWidth',2)

xlabel('time (second)')
ylabel('total energy (Pa)')
legend('two particles with same mineral', 'two particles with different mineral')
