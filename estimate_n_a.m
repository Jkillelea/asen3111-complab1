function [N, A] = estimate_n_a(cp_upper, cp_lower, n_iters)
  airfoil = @(x, c) (0.12 .* c ./ 0.2) .* (0.2969.*sqrt(x./c) - 0.126.*sqrt(x./c) - 0.3516.*(x./c).^2 + 0.2843.*(x./c).^3 - 0.1036.*(x./c).^4); % airfoil polynomial

  chord = 0.5; % meters
  airfoil = @(x) airfoil(x, chord);

  alpha   = 9;            % degrees
  v_inf   = 20;           % m/s
  rho_inf = 1.225;        % kg/m^3
  p_inf   = 10.13 * 10^4; % Pa
  q_inf   = 0.5 * rho_inf * (v_inf^2);

  % n_iters = 10000;
  x = linspace(0, chord, n_iters+1);
  y = airfoil(x);

  % normal and axial force on upper surface of airfoil
  n_acc = 0;
  a_acc = 0;
  for k = 1:n_iters
    x_1 = x(k);
    y_1 = y(k);
    x_2 = x(k+1);
    y_2 = y(k+1);

    dx  = x_2 - x_1;
    dy  = y_2 - y_1;

    Cp_1 = fnval(cp_upper, x_1/chord);
    Cp_2 = fnval(cp_upper, x_2/chord);

    P_1 = Cp_1 .* q_inf + p_inf;
    P_2 = Cp_2 .* q_inf + p_inf;
    P   = (P_1 + P_2) / 2;

    n_acc = n_acc - dx * P;
    a_acc = a_acc - dy * P; % when dy is negative this contributes to axial force
  end
  N_upper = n_acc;
  A_upper = a_acc;

  % normal and axial force on lower surface
  n_acc = 0;
  a_acc = 0;
  for k = 1:n_iters
    x_1 = x(k);
    y_1 = -y(k);
    x_2 = x(k+1);
    y_2 = -y(k+1);

    dx  = x_2 - x_1;
    dy  = -(y_2 - y_1);

    Cp_1 = fnval(cp_lower, x_1/chord);
    Cp_2 = fnval(cp_lower, x_2/chord);

    P_1 = Cp_1 .* q_inf + p_inf;
    P_2 = Cp_2 .* q_inf + p_inf;
    P   = (P_1 + P_2) / 2;

    n_acc = n_acc + dx * P;
    a_acc = a_acc + dy * P; % when dy is positive this contributes to axial force
  end
  N_lower = n_acc;
  A_lower = a_acc;

  N = N_lower + N_upper;
  A = A_lower + A_upper;
end
