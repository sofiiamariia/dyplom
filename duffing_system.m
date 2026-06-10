function dxdt = duffing_system(t, x)
  dxdt = zeros(2,1);
  dxdt(1) = x(2);
  dxdt(2) = -0.2 * x(2) - x(1)^3 + 0.3 * cos(t);
end

