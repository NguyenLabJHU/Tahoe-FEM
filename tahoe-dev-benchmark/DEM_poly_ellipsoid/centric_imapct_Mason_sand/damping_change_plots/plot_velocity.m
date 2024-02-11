% used to plot average velocity of the poly result and its corresponding
% ellip result together for comparison

clc
clear
close all


no_back_no_contact_result = load('no_contact_background_damping_dep_progress');
contact_05_result = load('contact_damping_0.3_dep_progress');
contact_1_result = load('contact_damping_0.5_dep_progress');
contact_3_result = load('contact_damping_0.7_dep_progress');
contact_5_result = load('contact_damping_0.9_dep_progress');

%interval_print = 10;
interval_plot= 1;
timestep = 1.0e-7;

% no background no cantact damping
no_back_no_contact_result = no_back_no_contact_result(1:interval_plot:end,:);

no_back_no_contact_timestep = no_back_no_contact_result(:,1);
no_back_no_contact_velocity = no_back_no_contact_result(:,7);
no_back_no_contact_kinetic = no_back_no_contact_result(:,13);   % kinetic energy
no_back_no_contact_potential = no_back_no_contact_result(:,13); % potential energy
no_back_no_contact_total = no_back_no_contact_result(:,15);     % total energy

no_back_no_contact_time = no_back_no_contact_timestep*timestep;


% 0.05 contact damping
contact_05_result = contact_05_result(1:interval_plot:end,:);

contact_05_timestep = contact_05_result(:,1);
contact_05_velocity = contact_05_result(:,7);
contact_05_kinetic = contact_05_result(:,13);   % kinetic energy
contact_05_potential = contact_05_result(:,13); % potential energy
contact_05_total = contact_05_result(:,15);     % total energy

contact_05_time = contact_05_timestep*timestep;


% 0.1 contact damping
contact_1_result = contact_1_result(1:interval_plot:end,:);

contact_1_timestep = contact_1_result(:,1);
contact_1_velocity = contact_1_result(:,7);
contact_1_kinetic = contact_1_result(:,13);   % kinetic energy
contact_1_potential = contact_1_result(:,13); % potential energy
contact_1_total = contact_1_result(:,15);     % total energy

contact_1_time = contact_1_timestep*timestep;


% 0.3 contact damping
contact_3_result = contact_3_result(1:interval_plot:end,:);

contact_3_timestep = contact_3_result(:,1);
contact_3_velocity = contact_3_result(:,7);
contact_3_kinetic = contact_3_result(:,13);   % kinetic energy
contact_3_potential = contact_3_result(:,13); % potential energy
contact_3_total = contact_3_result(:,15);     % total energy

contact_3_time = contact_3_timestep*timestep;


% 0.5 contact damping
contact_5_result = contact_5_result(1:interval_plot:end,:);

contact_5_timestep = contact_5_result(:,1);
contact_5_velocity = contact_5_result(:,7);
contact_5_kinetic = contact_5_result(:,13);   % kinetic energy
contact_5_potential = contact_5_result(:,13); % potential energy
contact_5_total = contact_5_result(:,15);     % total energy

contact_5_time = contact_5_timestep*timestep;

figure(1)
hold on
plot(no_back_no_contact_time(1:end/2), no_back_no_contact_kinetic(1:end/2), 'r-', 'LineWidth', 2)
plot(contact_05_time(1:end/2), contact_05_kinetic(1:end/2), 'g-', 'LineWidth', 2)
plot(contact_1_time(1:end/2), contact_1_kinetic(1:end/2), 'b-', 'LineWidth', 2)
plot(contact_3_time(1:end/2), contact_3_kinetic(1:end/2), 'm-', 'LineWidth', 2)
plot(contact_5_time(1:end/2), contact_5_kinetic(1:end/2), 'k-', 'LineWidth', 2)
%axis([0 0.25 0 3.5e-6])
xlabel('time (s)')
ylabel('kinetic energy (J)')
legend('simulation without background and contact damping', 'simulation with 0.3 contact damping',...
       'simulation with 0.5 contact damping', 'simulation with 0.7 contact damping',...
       'simulation with 0.9 contact damping')