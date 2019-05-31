% MAIN
%
% Script for performing axon growth simulations. The script implements the
% Peaceman--Rachford and the explicit Euler time discretizations of the
% fully scaled axon growth model. All described in 
%
% S. Diehl, E. Henningsson, A. Heyden:
% "Efficient simulations of tubulin-driven axonal growth"
% Journal of Computational Neuroscience (2016)
%
% The files were created with and run with MATLAB version 8.3 (R2014a).
%
%
% INPUT: 
%
% To run the code the user needs to supply a file cs.m which describes the
% soma concentration as a function of time. (An example file can be
% downloaded together with this main file.)
%
% Physical, biological, and numerical parameters can be changed as desired 
% below. Nominal values are set by default, cf. the article referenced 
% above. Input should be given using SI units.
%
%
% OUTPUT: (All output is given in SI units.)
%
% tau:  Vector containing the scaled time points.
%
% t:    Vector containing the original time as a function of the scaled
%       time tau.
%
% cc:   Vector containing the growth cone tubulin concentration as a 
%       function of tau.
%
% l:	Vector containing the axon length as a function of tau.
%
% y:    Vector containing the scaled spatial points. (I.e., a 
%       discretization of the interval (0,1).)
%
% c:    Matrix containing the tubulin concentration along the axon as a
%       function of the scaled spatial variable y and scaled temporal
%       variable tau.
%
%
% EXAMPLE: A plot of the axon length vs. original time is given by first 
%          running this script and then calling:
%          plot(t, l)
%
% EXAMPLE: The original space variable can be obtained at a given time
%          point t(i) by the formula:
%          x = y*l(i)
%
%
% Erik Henningsson
% March, 2016
% Lund University, Sweden
% erikh@maths.lth.se

if ~exist('hold_off') % if hold_off does not exists then clear and close
    clear
    close all
end

%%% MODIFY AFTER YOUR DESIRE THE METHOD CHOICE AND PARAMETER VALUES BELOW %%%


%%% Choose time stepping method %%%

% 0: Peaceman--Rachford, 1: Explicit Euler, 2: Peaceman--Rachford with
% extra time scaling.
%
% For methods 0 and 1 the time is scaled with the advection rate. All 
% parameter studies in the above mentioned article are performed using this
% scaling. For method 2 the time is instead scaled with the diffusion rate. 
% This gives a numerical scheme that uses far more time steps at short axon
% lengths. This option can be used for studying transient behavoiur but it
% is not recommended to use it when studying the convergence to steady
% state.
%
% ATTENTION: The numerical parameters below are chosen for
% Peaceman--Rachford. If explicit Euler is to be used either k or M has to
% be chosen substantially smaller. We strongly suggest to use explicit
% Euler only if you are interested in comparing numerical methods.

method = 0;


%%% Choose parameters %%%

% Biological and physical parameters

a = 1e-8;
D = 10e-12;
if ~exist('g')
    g = 5e-7;
end
lc = 4e-6;
rg = 1.783e-5;
cinf = 11.9e-3;
rgt = 0.053;


% Initial and boundary data

l0 = 1e-6; % Initial axon length.
cs0 = 2*cinf; % Initial soma concentration, see also the file cs.m!
cc0 = 2*cinf; % Initial cone concentration.
cinit = @(y) cc0*ones(size(y)); % Initial concentration profile along the
    % axon. Denoted by c^0 in article.
T = 6e8; % End time.


% Numerical parameters

M = 1000 - 1; % Number spatial points. 
k = 5e-4; % Time step size, denoted \Delta\tau in article.


% Storage parameters

% For big values of M and small values of k the code will downsample the
% output to limit the memory useage.
M_max = 10000 - 1; % This controls the amount of spatial data that is 
    % stored when M is chosen large.
N_max = 10000; % This controls the amount of temporal data that is stored
    % when k is chosen small. (Actually when N_guess is chosen large.)

if exist('hold_off') % if hold_off exists then set N_guess as
    N_guess = round(2000/k); % appropriate for use by parameter_studies.m (see below)
else
    N_guess = round(200/k); % The user is asked to supply a guess of how many
    % steps the scheme will take. A too large value of N_guess means that
    % more memory than required is used during simulation. A too small
    % value means slow simulation as reallocation is performed each step
    % after the preallocation is filled up. The value set here is for when
    % axon_growth_... is run by itself (before parameter_studies.m)
end

include_alls = false; % If true, all but concentration along the axon will 
    % be stored without downsampling to the variables tall, lall and ccall.

   
%%% NO MORE USER OPTIONS AFTER THIS LINE %%%
    

%%% Initialize %%%

% Define spatial variable
yext = linspace(0,1,M+2)';
y = yext(2:end-1);
h = y(2)-y(1); % Spatial step, denoted \Delta y in article.

N = N_guess;
    
% Allocate memory
len = ceil(N/round(N/N_max));
if M <= M_max,  c = zeros(M, min(len+1, N+1));
else            c = zeros(M_max, min(len+1, N+1)); end
j = 1;

tau = zeros(min(len+1, N+1),1);
t = zeros(min(len+1, N+1),1);
cc = zeros(min(len+1, N+1),1);
l = zeros(min(len+1, N+1),1);

% Set initial data
n = 1;
taun1 = 0;
ccn1 = cc0;
ln1 = l0;
tn1 = 0;
cn1 = cinit(y);

t(1) = 0;
if M <= M_max,  c(:,1) = cinit(y);
else            step = (M+1)/(M_max+1); c(:,1) = cinit(y(step:step:end-step+1)); end
cc(1) = cc0;
l(1) = l0;

% Allocate memory and set initial data for "all"-variables
if include_alls
    lall = zeros(N+1,1);
    ccall = zeros(N+1,1);
    tall = zeros(N+1,1);
    
    lall(1) = l0;
    ccall(1) = cc0;
    tall(1) = 0;
end

% Construct finite difference matrices

I = speye(M,M);
ett = ones(M,1);
TL = spdiags([ett -2*ett ett]/h^2, [-1,0,1], M, M);
TN = spdiags([-ett ett]/2/h, [-1, 1], M, M);
Y = spdiags(y, 0, M, M);
    

%%% Time stepping schemes %%%

if method == 0

    %%% Time stepping -- Peaceman-Rachford %%%
    
    ms = zeros(1,N);
    
    tic;
    while tn1 < T
        
        % Output progress
        if mod(n,1e4) == 0
            disp(['Step: ' num2str(n) '. Part of time: ' num2str(tn1/T) '. Axon length: ' num2str(ln1) '.'])
        end
        if n == N_guess, warning('Preallocation filled up.'); end
        
        cold = cn1;
        ccold = ccn1;
        lold = ln1;
        tauold = taun1;
        told = tn1;
        
        % Explicit Euler for ODE:s
        coldM_y = (3*ccold - 4*cold(end) + cold(end-1))/2/h;
        ccnew = ccold + k/2*((a-g*lc)*ccold - D/lold*coldM_y - (rg*ccold + rgt*lc)*(ccold - cinf))*lold/a/lc;
        lnew = lold + k/2*rg*(ccold - cinf)*lold/a;
        tnew = told + k/2*lold/a;
        
        % Implicit Euler for PDE:s
        ccold = ccnew; lold = lnew; told = tnew;
        cstold = cs(told,cs0);
        alpha = (a - y*rg*(ccold - cinf))/lold;
        alpha_ = spdiags(alpha, 0, M, M); 
        rhs = cold;
        rhs(1) = rhs(1) + k/2*lold/a*(D/(h*lold)^2 + alpha(1)/2/h)*cstold;
        rhs(end) = rhs(end) + k/2*lold/a*(D/(h*lold)^2 - alpha(end)/2/h)*ccold;
        cnew = ((1+k/2*lold/a*g)*I - k/2*lold/a*D/lold^2*TL + k/2*lold/a*alpha_*TN)\rhs;
        
        % Explicit Euler for PDE:s
        cold = cnew;
        rhs = cold;
        rhs(1) = rhs(1) + k/2*lold/a*(D/(h*lold)^2 + alpha(1)/2/h)*cstold;
        rhs(end) = rhs(end) + k/2*lold/a*(D/(h*lold)^2 - alpha(end)/2/h)*ccold;
        cnew = rhs + k/2*lold/a*(-g*cold + D/lold^2*(TL*cold) - alpha.*(TN*cold));
        
        % Implicit Euler for ODE:s
        cold = cnew;
        u = [ccold; lold];
        % Newton iteration
        for m = 1:1e3
            c_y = (3*u(1) - 4*cold(end) + cold(end-1))/2/h;
            Fu = [(1 - k/2*u(2)/a/lc*(a-g*lc))*u(1) - ccold + k/2*u(2)/a/lc*((rg*u(1)+rgt*lc)*(u(1)-cinf)) + k/2/a/lc*D*c_y;
                u(2) - lold - k/2*u(2)/a*rg*(u(1) - cinf)];
            Ju = [(1 - k/2*u(2)/a/lc*(a-g*lc)) + k/2*u(2)/a/lc*(2*rg*u(1) - rg*cinf + rgt*lc) + k/2/a/lc*D*3/2/h, -k/2/a/lc*(a-g*lc)*u(1)+k/2/a/lc*((rg*u(1)+rgt*lc)*(u(1)-cinf));
                -k/2*u(2)/a*rg, 1 - k/2/a*rg*(u(1) - cinf)];
            uold = u;
            u = u - Ju\Fu;
            if all(abs(u - uold)./abs(u) < 1e-14)
                break;
            end
        end
        ms(n) = m;
        if m >= 1e3, error('Newton diverged.'); end
        tn1 = told + k/2*u(2)/a;
        
        cn1 = cold;
        ccn1 = u(1);
        ln1 = u(2);
        taun1 = tauold + k;
        
        % Store temporal slice of solution
        if N < N_max || mod(n, round(N/N_max)) == 0 || n == N
            if M <= M_max
                c(:,j+1) = cn1;
            else
                c(:,j+1) = cn1(step:step:end-step+1);
            end
            cc(j+1) = ccn1;
            l(j+1) = ln1;
            tau(j+1) = taun1;
            t(j+1) = tn1;
            j = j+1;
        end
        
        % Store each time step
        if include_alls
            lall(n+1) = ln1;
            ccall(n+1) = ccn1;
            tall(n+1) = tn1;
        end
        
        n = n+1;
    end
    time = toc;
    
    ms = ms(1:n-1);
end


if method == 1

    %%% Time stepping -- Explicit Euler %%%
    
    % In CFL an approximation of the left hand side of a CFL condition is
    % stored. At each time point it must be below 1/2 or numerical
    % stability may occur.
    initial_CFL = k*D/a/l0/h^2
    CFL = zeros(min(len+1, N+1),1);
    CFL(1) = initial_CFL;
    
    num = 1;
    
    tic;
    while tn1 < T
        
        % Output progress
        if tn1 >= num*T/10
            disp(tn1)
            num = num + 1;
        end
        if mod(n,1e5) == 0
            disp(['Step: ' num2str(n) '. Part of time: ' num2str(tn1/T) '. Axon length: ' num2str(ln1) '.'])
        end
        if n == N_guess, error('Allocate more. /Erik'); end
        
        cn = cn1;
        ccn = ccn1;
        ln = ln1;
        tn = tn1;
        
        taun = taun1;
        
        % Explicit Euler
        taun1 = taun + k;
        
        cMn_y = (3*ccn - 4*cn(end) + cn(end-1))/2/h;
        ccn1 = ccn + k*((a-g*lc)*ccn - D/ln*cMn_y - (rg*ccn + rgt*lc).*(ccn - cinf))*ln/a/lc;
        ln1 = ln + k*rg*(ccn - cinf)*ln/a;
        tn1 = tn + k/a*ln;
        
        alphan = (a - y*rg*(ccn - cinf))/ln;
        rhs = cn;
        rhs(1) = rhs(1) + k*ln/a*(D/(h*ln)^2 + alphan(1)/2/h)*cs(tn,cs0);
        rhs(end) = rhs(end) + k*ln/a*(D/(h*ln)^2 - alphan(end)/2/h)*ccn;
        cn1 = rhs - k*ln/a*g*cn + k*ln/a*D/ln^2*(TL*cn) - k*ln/a*alphan.*(TN*cn);
        
        % Store temporal slice of solution
        if N < N_max || mod(n, round(N/N_max)) == 0 || n == N
            c(:,j) = cn1;
            cc(j+1) = ccn1;
            l(j+1) = ln1;
            tau(j+1) = taun1;
            t(j+1) = tn1;
            CFL(j+1) = k*D/a/ln1/h^2;
            j = j+1;
        end
        
        % Store each time step
        if include_alls
            lall(n+1) = ln1;
            ccall(n+1) = ccn1;
            tall(n+1) = tn1;
        end
        
        n = n+1;
    end
    time = toc;
    
    % Compute CFL approximation.
    CFL(j+1) = k*D/a/ln1/h^2;
    CFL = CFL(1:j);
end


if method == 2
    
    %%% Time stepping -- Peaceman-Rachford -- extra time scaling %%%
    
    ms = zeros(1,N);
    cold = c(:,1);
    
    tic;
    while tn1 < T
        
        % Output progress
        if mod(n,1e4) == 0
            disp(['Step: ' num2str(n) '. Part of time: ' num2str(tn1/T) '. Axon length: ' num2str(ln1) '.'])
        end
        if n == N_guess, warning('Preallocation filled up.'); end
        
        cold = cn1;
        ccold = ccn1;
        lold = ln1;
        tauold = taun1;
        told = tn1;
        
        % Explicit Euler for ODE:s
        coldM_y = (3*ccold - 4*cold(end) + cold(end-1))/2/h;
        ccnew = ccold + k/2*((a-g*lc)*ccold - D/lold*coldM_y - (rg*ccold + rgt*lc)*(ccold - cinf))*lold^2/D/lc;
        lnew = lold + k/2*rg*(ccold - cinf)*lold^2/D;
        tnew = told + k/2*lold^2/D;
        
        % Implicit Euler for PDE:s
        ccold = ccnew; lold = lnew; told = tnew;
        cstold = cs(told,cs0);
        alpha = (a - y*rg*(ccold - cinf))/lold;
        alpha_ = spdiags(alpha, 0, M, M);
        rhs = cold;
        rhs(1) = rhs(1) + k/2*lold^2/D*(D/(h*lold)^2 + alpha(1)/2/h)*cstold;
        rhs(end) = rhs(end) + k/2*lold^2/D*(D/(h*lold)^2 - alpha(end)/2/h)*ccold;
        cnew = ((1+k/2*lold^2/D*g)*I - k/2*TL + k/2*lold^2/D*alpha_*TN)\rhs;
        
        % Explicit Euler for PDE:s
        cold = cnew;
        rhs = cold;
        rhs(1) = rhs(1) + k/2*lold^2/D*(D/(h*lold)^2 + alpha(1)/2/h)*cstold;
        rhs(end) = rhs(end) + k/2*lold^2/D*(D/(h*lold)^2 - alpha(end)/2/h)*ccold;
        cnew = rhs + k/2*lold^2/D*(-g*cold + D/lold^2*(TL*cold) - alpha.*(TN*cold));
        
        % Implicit Euler for ODE:s
        cold = cnew;
        u = [ccold; lold];
        for m = 1:1e3
            c_y = (3*u(1) - 4*cold(end) + cold(end-1))/2/h;
            Fu = [(1 - k/2*u(2)^2/D/lc*(a-g*lc))*u(1) - ccold + k/2*u(2)^2/D/lc*((rg*u(1)+rgt*lc)*(u(1)-cinf) + D*c_y/u(2));
                u(2) - lold - k/2*u(2)^2/D*rg*(u(1) - cinf)];
            Ju = [(1 - k/2*u(2)^2/D/lc*(a-g*lc)) + k/2*u(2)^2/D/lc*((2*rg*u(1) - rg*cinf + rgt*lc) + D*3/(2*h*u(2))), ...
                    -k*u(2)/D/lc*(a-g*lc)*u(1) + k*u(2)/D/lc*((rg*u(1)+rgt*lc)*(u(1)-cinf) + D*c_y/u(2)) - k/2*u(2)^2/lc*c_y/u(2)^2;
                -k/2*u(2)^2/D*rg, 1 - k*u(2)/D*rg*(u(1) - cinf)];
            uold = u;
            u = u - Ju\Fu;
            if all(abs(u - uold)./abs(u) < 1e-14)
                break;
            end
        end
        ms(n) = m;
        if m >= 1e3, error('Newton diverged.'); end
        
        cn1 = cold;
        ccn1 = u(1);
        ln1 = u(2);
        tn1 = told + k/2*ln1^2/D;
        taun1 = tauold + k;
        
        % Store temporal slice of solution
        if N < N_max || mod(n, round(N/N_max)) == 0 || n == N
            if M <= M_max
                c(:,j+1) = cn1;
            else
                c(:,j+1) = cn1(step:step:end-step+1);
            end
            cc(j+1) = ccn1;
            l(j+1) = ln1;
            tau(j+1) = taun1;
            t(j+1) = tn1;
            j = j+1;
        end
        
        % Store each time step
        if include_alls
            lall(n+1) = ln1;
            ccall(n+1) = ccn1;
            tall(n+1) = tn1;
        end
        
        n = n+1;
    end
    time = toc;
end


%%% Post processing %%%

% Store the last time step
if t(j) ~= tn1
    if M <= M_max
        c(:,j+1) = cn1;
    else
        c(:,j+1) = cn1(step:step:end-step+1);
    end
    cc(j+1) = ccn1;
    l(j+1) = ln1;
    tau(j+1) = taun1;
    t(j+1) = tn1;
    j = j+1;
end

% Shrink the solution vectors to correct size
tau = tau(1:j);
t = t(1:j);
cc = cc(1:j);
l = l(1:j);
c = c(:,1:j);

if include_alls
    tall = tall(1:n);
    lall = lall(1:n);
    ccall = ccall(1:n);
end
