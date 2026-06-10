function [traj, xT, Phi] = loop_rk4(f, Jac_f, t_span, h, x0)
    N = length(t_span);
    m = length(x0);
    I = eye(m);

    traj = zeros(m, N);
    traj(:,1) = x0;
    Phi = I;

    for j = 1:N-1
        t_curr = t_span(j);
        t_next = t_span(j+1);
        u = traj(:,j);

        k1 = f(t_curr, u);
        k2 = f(t_curr + h/2, u + h/2 * k1);
        k3 = f(t_curr + h/2, u + h/2 * k2);
        k4 = f(t_next, u + h * k3);

        traj(:,j+1) = u + (h/6)*(k1 + 2*k2 + 2*k3 + k4);


        J1 = Jac_f(t_curr, u);
        J2 = Jac_f(t_curr + h/2, traj(:,j) + h/2 * k1);
        J3 = Jac_f(t_curr + h/2, traj(:,j) + h/2 * k2);
        J4 = Jac_f(t_next, traj(:,j) + h * k3);

        K1 = J1 * Phi;
        K2 = J2 * (Phi + h/2 * K1);
        K3 = J3 * (Phi + h/2 * K2);
        K4 = J4 * (Phi + h * K3);

        Phi = Phi + (h/6)*(K1 + 2*K2 + 2*K3 + K4);
    end

    xT = traj(:,N);
end
