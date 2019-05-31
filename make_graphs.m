% make_graphs.m
% these were not the original commands to make paper figures but
% make figures that resemble three in the paper:
% Fig 2a, 2b, 3b


%%% Axon length, l %%%

figure
plot(t/24/3600, l*1000)
title('Fig 2a')
xlabel('Time [days]')
ylabel('Axon length [mm]')


%%% Concentration along the axon, c %%%

% Down-sampling of the solution for cheaper 3D-plot.
ndisp = 10; % In time by a factor 10.
ndispy = 10; % In space by a factor 10.

y_down = [0; y(ndispy:ndispy:end-ndispy+1); 1];
x3D = 1000*[l0; l(ndisp:ndisp:end)]*y_down';
t_down = [0; t(ndisp:ndisp:end)];
css = zeros(1,size(x3D,1));
for kk = 1:length(css)
    css(kk) = cs(t_down(kk),cs0);
end
c_down = [css; ...
    cinit(y_down(2:end-1)) c(ndispy:ndispy:end-ndispy+1,ndisp:ndisp:end); ...
    [cc(1); cc(ndisp:ndisp:end)]'];

figure
mesh(x3D, t_down/24/3600, c_down')
title('Fig 2b')
xlabel('x [mm]')
ylabel('Time [days]')
zlabel('Concentration [mol/m^3]')


%%% Growth cone concentration, cc %%%

figure
plot(t,cc)
axis([0 3500 0.01 0.024])
grid on
title('Fig 3b')
xlabel('Time [s]')
ylabel('Growth cone concentration [mol/m^3]')