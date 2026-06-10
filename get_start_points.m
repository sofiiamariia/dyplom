function X_start = get_start_points(f, t0, x0, h, r)
    b_start = t0 + (r-1)*h;
    [~, X_start, ~] = runge_kutta_system(f, t0, b_start, h, x0);
end
