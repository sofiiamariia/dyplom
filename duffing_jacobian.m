function J = duffing_jacobian(t, x)
  J = [0, 1; -3*x(1)^2, -0.2];
end

