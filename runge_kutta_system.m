function [t_values, x_values, rk_steps] = runge_kutta_system(f, a, b, h, x0)
    t_values = a:h:b;
    n_steps = length(t_values);
    n_eqs = length(x0);

    x_values = zeros(n_eqs, n_steps);
    x_values(:, 1) = x0;
    rk_steps = 0;

    for i = 1:n_steps-1
        t = t_values(i);
        x = x_values(:, i);

        k1 = f(t, x);
        k2 = f(t + h/2, x + h/2 * k1);
        k3 = f(t + h/2, x + h/2 * k2);
        k4 = f(t + h, x + h * k3);

        x_next = x + h/6 * (k1 + 2*k2 + 2*k3 + k4);

        if any(isnan(x_next)) || any(isinf(x_next))
            warning('NaN or Inf detected at step %d (t = %.5f)', i+1, t_values(i+1));
            x_values(:, i+1:end) = NaN;
            break;
        end

        x_values(:, i+1) = x_next;
        rk_steps = rk_steps + 1;
    end
end

