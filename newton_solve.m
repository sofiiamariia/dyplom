function [u, count] = newton_solve(Func, Jac, u_start, tol)
    u = u_start;
    max_iter = 1000;
    count = 0;

    for iter = 1:max_iter
        count = iter;
        F_val = Func(u);
        J_val = Jac(u);

        delta = -J_val \ F_val;
        u = u + delta;

        if norm(delta, inf) < tol
            return;
        end
    end
end


