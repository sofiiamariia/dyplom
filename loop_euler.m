function [traj, xT, Phi_inv] = loop_euler(f, Jac_f, t_span, h, x0, tol)
    N = length(t_span);
    m = length(x0);
    I = eye(m);

    traj = zeros(m, N);
    traj(:,1) = x0;

    Phi_inv = I;

    for j = 1:N-1
        t_next = t_span(j+1);
        u0     = traj(:,j);

        F = @(u) u - u0 - h * f(t_next, u);
        J = @(u) I  - h * Jac_f(t_next, u);
        [u_next, ~] = newton_solve(F, J, u0, tol);
        traj(:,j+1) = u_next;

        J_next  = Jac_f(t_next, u_next);
        Phi_inv = Phi_inv * (I - h * J_next);
    end

    xT = traj(:,N);
end
