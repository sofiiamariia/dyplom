function [traj, xT, Phi] = loop_g3(f, Jac_f, t_span, h, a, x0, tol)
    N = length(t_span);
    m = length(x0);
    I = eye(m);

    traj = zeros(m, N);

    X_start = get_start_points(f, t_span(1), x0, h, 3);
    traj(:,1) = X_start(:,1);
    traj(:,2) = X_start(:,2);
    traj(:,3) = X_start(:,3);

    Phi_prev2_m = I;
    J1 = Jac_f(t_span(1), traj(:,1));
    Phi_prev1_m = (I + h*J1) * Phi_prev2_m;

    J2 = Jac_f(t_span(2), traj(:,2));
    Phi = (I + h*J2) * Phi_prev1_m;

    u_prev2 = traj(:,1);
    u_prev1 = traj(:,2);
    u_curr  = traj(:,3);

    for j = 4:N
        t_next = t_span(j);

        F_G = @(u) (11/6)*u - 3*u_curr + (3/2)*u_prev1 - (1/3)*u_prev2 - h*f(t_next, u);
        J_G = @(u) (11/6)*I - h*Jac_f(t_next, u);
        [u_next, ~] = newton_solve(F_G, J_G, u_curr, tol);
        traj(:,j) = u_next;

        J_next = Jac_f(t_next, u_next);
        LHS = (11/6)*I - h*J_next;
        RHS = 3*Phi - (3/2)*Phi_prev1_m + (1/3)*Phi_prev2_m;
        Phi_new = LHS \ RHS;

        Phi_prev2_m = Phi_prev1_m;
        Phi_prev1_m = Phi;
        Phi = Phi_new;

        u_prev2 = u_prev1;
        u_prev1 = u_curr;
        u_curr  = u_next;
    end

    xT = traj(:,N);
end
