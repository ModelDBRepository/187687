% parameter_studies.m
%
% This script gives an example on how the code can be used to perform
% parameter studies. This code will run the main file five times for
% different values on the degredation rate g and reproduce Fig. 9 of 
% the original paper.
%
%
% Note: For faster simulations for the extreme values of g
% the preallocation guess is increased in the main file:
%
%   axon_growth_simulations_with_time_and_space_scaling.m
%   line 129 is selected: where N_guess is round(2000/k);
%
% The code will take several minutes to run since, for big values of g, the
% axon never growths very long and therefore really big time steps (in
% original time) are never taken. A solution might be to run these 
% simulations with a bigger time step (in scaled time), denoted by k in the
% main file.
%
%
% Erik Henningsson
% May, 2016
% Lund University, Sweden
% erikh@maths.lth.se

clear
close all

gs = [1/2 1 2 4 8]*5e-7;

figure
hold on
colors = 'bcrgm';
fig_handle=figure();
for iii = 1:length(gs)
    g = gs(iii);
    disp(['Starting simulation with g = ' num2str(g) '.'])
    hold_off = true; % prevents clearing and closing in axon_...
    axon_growth_simulation_with_time_and_space_scaling;
    figure(fig_handle);
    hold on
    plot(t/24/3600, l*1000, colors(iii))
end

title('Fig 9')
xlabel('Time [days]')
ylabel('Axon length [mm]')
