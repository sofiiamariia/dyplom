function dxdt = vdp_forced_system(t, x, mu)
     E = 1;
     T = 2*pi;
     omega = 2*pi / T;

    dxdt = zeros(2,1);
    dxdt(1) = x(2);
    dxdt(2) = mu * (1 - x(1)^2) * x(2) - x(1) + E * cos(omega * t);
end
