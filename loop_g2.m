function [traj, xT, Phi] = loop_g2(f, Jac_f, t_span, h, x0, tol)
    N = length(t_span);
    m = length(x0);
    I = eye(m);

    traj = zeros(m, N);
    traj(:,1) = x0;

    t1 = t_span(2);
    u0 = x0;
    F_E = @(u) u - u0 - h*f(t1, u);
    J_E = @(u) I - h*Jac_f(t1, u);
    [u1, ~] = newton_solve(F_E, J_E, u0, tol);
    traj(:,2) = u1;

    J1   = Jac_f(t1, u1);
    Phi  = (I - h*J1) \ I;
    Phi_prev = I;

    u_prev = u0;
    u_curr = u1;

    for j = 3:N
        t_next = t_span(j);

        F_G = @(u) (3/2)*u - 2*u_curr + (1/2)*u_prev - h*f(t_next, u);
        J_G = @(u) (3/2)*I - h*Jac_f(t_next, u);
        [u_next, ~] = newton_solve(F_G, J_G, u_curr, tol);
        traj(:,j) = u_next;

        J_next  = Jac_f(t_next, u_next);
        LHS = (3/2)*I - h*J_next;
        RHS = 2*Phi - (1/2)*Phi_prev;
        Phi_new = LHS \ RHS;

        Phi_prev = Phi;
        Phi = Phi_new;

        u_prev = u_curr;
        u_curr = u_next;
    end

    xT = traj(:,N);
end
