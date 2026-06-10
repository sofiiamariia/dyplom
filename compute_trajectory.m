function traj = compute_trajectory(f, Jac_f, t_span, h, a, x0, tol)
    [traj, ~, ~] = outer_loop(f, Jac_f, t_span, h, a, x0, tol);
end
