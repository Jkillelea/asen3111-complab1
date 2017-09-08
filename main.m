% ASEN 3111 - Jacob Killelea [105510162] - Computational Lab 1

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIRST PROBLEM - Sphere in inviscid flow
clear all; clc;

radius  = 10;                        % meters
v_inf   = 25;                        % m/s
rho_inf = 0.9093;                    % kg/m^3
p_inf   = 7.012 * (10^4);            % Pa
q_inf   = 0.5 * rho_inf * (v_inf^2); % Pa

Cp = @(theta) 1 - 4.*(sin(theta).^2);

% We know the Cp distribution around the body and can work backwards to
% find the actual pressure distribution, given dynamic and static pressure
p = @(theta) Cp(theta) .* q_inf + p_inf;

up_force   = @(theta) -p(theta) .* sin(theta);
side_force = @(theta) -p(theta) .* cos(theta);

N_iters = 10000;
R       = radius; % meters
range   = [0, 2*pi];
h       = (2*pi)/N_iters;
t       = linspace(range(1), range(2), N_iters+1);

% Calculate lift from pressure distribution
accumulator = 0;
f = up_force;
for k = 1:(N_iters/2) % simpson's method for a line integral around a sphere
  t_1 = t(2*k-1);
  t_2 = t(2*k);
  t_3 = t(2*k+1);
  accumulator = accumulator + ( f(t_1) + 4*f(t_2) + f(t_3) );
end
lift = (h*R/3) * accumulator;

% Calculate drag from pressure distribution
accumulator = 0;
f = side_force;
for k = 1:(N_iters/2)
  t_1 = t(2*k-1);
  t_2 = t(2*k);
  t_3 = t(2*k+1);
  accumulator = accumulator + ( f(t_1) + 4*f(t_2) + f(t_3) );
end
drag = (h*R/3) * accumulator;

fprintf('Sphere of radius %d m: lift: %d N, drag: %d N.\n', radius, lift, drag);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECOND PROBLEM - Airfoil pressure distribution

clear all;
load Cp;

airfoil = @(x, c) (0.12 .* c ./ 0.2) .* (0.2969.*sqrt(x./c) - 0.126.*sqrt(x./c) - 0.3516.*(x./c).^2 + 0.2843.*(x./c).^3 - 0.1036.*(x./c).^4); % airfoil polynomial

alpha = 9; % degrees
chord = 0.5; % meters
airfoil = @(x) airfoil(x, chord);


truth_iterations = 50000; % fifty thousand iterations for 'truth' measurement
true_lift = estimate_lift(Cp_upper, Cp_lower, alpha, truth_iterations);

fprintf('Airfoil produces %.1f N of lift per unit span at %d degrees AoA.\n', true_lift, alpha);

% Progressively increase the number of iterations and get the error of that estimate
for n = 5:5:200
  lift = estimate_lift(Cp_upper, Cp_lower, alpha, n);
  percent_err = 100*abs((true_lift - lift)/true_lift);
  fprintf('%d panels: %f %% error \n', n, percent_err);
end