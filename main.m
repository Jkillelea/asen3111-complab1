% ASEN 3111 - Jacob Killelea [105510162] - Computational Lab 1

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIRST PROBLEM - Sphere in inviscid flow
clear all; clc;

radius  = 1;                        % meters
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

% Calculate lift and drag from pressure distribution
lift_accum = 0;
drag_accum = 0;
lif = up_force;
drg = side_force;
for k = 1:(N_iters/2) % simpson's method for a line integral around a sphere
  t_1 = t(2*k-1);
  t_2 = t(2*k);
  t_3 = t(2*k+1);
  lift_accum = lift_accum + ( lif(t_1) + 4*lif(t_2) + lif(t_3) );
  drag_accum = drag_accum + ( drg(t_1) + 4*drg(t_2) + drg(t_3) );
end
lift = (h*R/3) * lift_accum;
drag = (h*R/3) * drag_accum;

fprintf('Sphere of radius %d m: lift: %d N, drag: %d N.\n', radius, lift, drag);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECOND PROBLEM - Airfoil pressure distribution

clear all;
load Cp;

% airfoil polynomial
airfoil = @(x, c) (0.12 .* c ./ 0.2) .* (0.2969.*sqrt(x./c)  ...
                                       - 0.1260.*sqrt(x./c)   ...
                                       - 0.3516.*(x./c).^2   ...
                                       + 0.2843.*(x./c).^3   ...
                                       - 0.1036.*(x./c).^4);

alpha = 9; % degrees
chord = 0.5; % meters
airfoil = @(x) airfoil(x, chord);

disp('Calculating lift and drag. This may take a bit...');
truth_iterations = 50000; % fifty thousand iterations for 'truth' measurement
[true_lift, true_drag] = estimate_lift(Cp_upper, Cp_lower, alpha, truth_iterations);

fprintf('Airfoil produces %.1f N of lift and %0.1f N of drag per unit span at %d degrees AoA.\n', ...
         true_lift, true_drag, alpha);

% Progressively increase the number of iterations and get the error of that estimate
errs        = [];
n           = 15;
percent_err = 100;
while percent_err > 5 % 5% error
  [lift, drag]   = estimate_lift(Cp_upper, Cp_lower, alpha, n);
  percent_err    = 100*abs((true_lift - lift)/true_lift);
  errs(:, end+1) = [n; percent_err];
  n = n + 1;
end
fprintf('%d panels: %f %% error \n', n, percent_err);

while percent_err > 1 % 1% error
  [lift, drag]   = estimate_lift(Cp_upper, Cp_lower, alpha, n);
  percent_err    = 100*abs((true_lift - lift)/true_lift);
  errs(:, end+1) = [n; percent_err];
  n = n + 1;
end
fprintf('%d panels: %f %% error \n', n, percent_err);

while percent_err > 0.1 % 0.1% error
  [lift, drag]   = estimate_lift(Cp_upper, Cp_lower, alpha, n);
  percent_err    = 100*abs((true_lift - lift)/true_lift);
  errs(:, end+1) = [n; percent_err];
  n = n + 5;
end
fprintf('%d panels: %f %% error \n', n, percent_err);

plot(errs(1, :), errs(2, :));
title('Percent error versus number of points used in simulation');
xlabel('Number of points used');
ylabel('Percent error');
print('err_vs_num', '-dpng');
