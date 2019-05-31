function concentration = cs(t, cs0)
%CS     Soma concentration as a function of time.
%   
% Input function to axon_growth_simulation_with_time_and_space_scaling.m.
% CS defines the soma concentration as a function of original time t. This
% file describes an example soma concentration with jumps between two
% different values. The user should modify this function to fit their
% simulations.
%
%
% INPUT:
%
% t:    The time point (in original time) to evaluate the soma 
%       concentration.
%
% cs0:  An optional parameter describing the initial soma concentration.
%       Can be used for easy scaling of the soma profile from the main
%       program.
%
%
% OUTPUT:
%
% concentration:    The soma concentration at time point t (original time).
%
%
% Erik Henningsson
% March, 2016
% Lund University, Sweden
% erikh@maths.lth.se


a = cs0;
b = cs0/4;
period = 2e8;

count = 0;
while t >= 0
    t = t - period;
    count = count + 1;
end

if mod(count,2) == 1 || count >= 4
    concentration = a;
else
    concentration = b;
end