function J = vdp_forced_jacobian(~, x, mu)
    J = zeros(2,2);
    J(1,1) = 0;
    J(1,2) = 1;
    J(2,1) = -2 * mu * x(1) * x(2) - 1;
    J(2,2) = mu * (1 - x(1)^2);
end
