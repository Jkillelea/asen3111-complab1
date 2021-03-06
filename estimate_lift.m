function [lift, drag] = estimate_lift(cp_upper, cp_lower, alpha, n_iters)
  [n, a] = estimate_n_a(cp_upper, cp_lower, n_iters);
  lift = n*cosd(alpha) - a*sind(alpha);
  drag = n*sind(alpha) + a*cosd(alpha);
end
