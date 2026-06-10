function [x0, xT, k, err_T_hist, err_x0_hist, traj] = gear_r3(f, Jac_f, a, b, h, x0, eps, tol, n)
    t_span = a:h:b;
    N = length(t_span);
    m = length(x0);
    I = eye(m);

    err   = inf;
    err_T = inf;
    max_k = 20;
    k = 0;
    xT = NaN(m,1);
    err_T_hist  = [];
    err_x0_hist = [];

    while (err >= tol || err_T >= eps) && k < max_k
        k = k + 1;
        x_prev = x0;

        [traj, xT, Phi] = loop_g3(f, Jac_f, t_span, h, a, x_prev, tol);

        Phi_inv = inv(Phi);
        A   = Phi_inv - I;
        rhs = Phi_inv * xT - x_prev;

        if rcond(A) < 1e-8
            x0 = pinv(A) * rhs;
        else
            x0 = A \ rhs;
        end

        err_T = norm(xT - x_prev);
        err   = norm(x0 - x_prev);

        err_T_hist(end+1)  = err_T;
        err_x0_hist(end+1) = err;

        fprintf('k = %d: ||x(T)-x0^k|| = %.3e, ||x0^{k+1}-x0^k|| = %.3e\n', k, err_T, err);
    end
    xT
##    [m1,m2]=size(traj);
##    m1=1;
##    for i=1:m1
##       plot(t_span, traj(i,:));hold on
##    endfor
##     plot(traj(1,:), traj(2,:));hold on

    if k == max_k
        warning('Досягнуто max ітерацій (%d). err_x0=%e, err_T=%e', max_k, err, err_T);
    end
end

