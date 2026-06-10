clc; clear; close all;

% Вибір системи: 1 - Дуффінг, 2 - Ван дер Поль
system_choice = 1;

eps_tol = 1e-5;
tol     = 1e-5;

if system_choice == 1
    fprintf('Система: Осцилятор Дуффінга\n');
    T = 2*pi;
    x0_guess = [-0.742; 0.729];
    f      = @duffing_system;
    Jac_f  = @duffing_jacobian;
elseif system_choice == 2
    fprintf('Система: Вимушений осцилятор Ван дер Поля \n');
    mu = 0.1;
    T = 2*pi;
    x0_guess = [0.1; 0.1];

    f      = @(t, x) vdp_forced_system(t, x, mu);
    Jac_f  = @(t, x) vdp_forced_jacobian(t, x, mu);

end

%% RK4, n = 2000
n_ref = 4000;
h_ref = T/n_ref;

[x0_ref, ~, ~, ~, ~, traj_ref] = rk4(f, Jac_f, 0, T, h_ref, x0_guess, eps_tol, tol, n_ref);
fprintf('"Еталонне" x0* = [%.10f; %.10f]\n\n', x0_ref(1), x0_ref(2));

%% n
if system_choice == 1
  n_values = [150];%, 200, 400, 600, 800];
%elseif system_choice == 3
%    n_values = [10000]
else
##  n_values = [150, 200, 400, 600, 800];
    n_values = [250];
end
num_n = length(n_values);
h_values = T ./ n_values;

errors = NaN(4, num_n);
times  = NaN(4, num_n);
iters  = NaN(4, num_n);

names = {'РK4', 'Неявний Ейлера', 'Гіра r=2', 'Гіра r=3'};
colors  = {'g', 'b', 'm', 'r'};
markers = {'o', 's', 'd', '^'};



for i = 1:num_n
    n = n_values(i);
    h = h_values(i);
    fprintf('--- Тестування для n = %d (h = %.4f) ---\n', n, h);

    % 1. RK4
    try
        tic;
        %[x0_1, ~, k_1, ~, ~] = rk4(f, Jac_f, 0, T, h, x0_guess, eps_tol, tol, n);
        [x0_1, xT_1, k_1, err_T_1, err_x0_1, traj_1] = rk4(f, Jac_f, 0, T, h, x0_guess, eps_tol, tol, n);
        times(1, i) = toc; errors(1, i) = norm(x0_1 - x0_ref); iters(1, i) = k_1;
    catch
        fprintf('  RK4 розійшовся!\n');
        times(1, i) = NaN; errors(1, i) = NaN; iters(1, i) = NaN;
    end

    % 2. Неявний м-д Ейлера
    try
        tic;
        %[x0_2, ~, k_2, ~, ~] = euler(f, Jac_f, 0, T, h, x0_guess, eps_tol, tol, n);
        [x0_2, xT_2, k_2, err_T_2, err_x0_2, traj_2] = euler(f, Jac_f, 0, T, h, x0_guess, eps_tol, tol, n);
        times(2, i) = toc; errors(2, i) = norm(x0_2 - x0_ref); iters(2, i) = k_2;
    catch
        fprintf('  Неявний Ейлер розійшовся!\n');
        times(2, i) = NaN; errors(2, i) = NaN; iters(2, i) = NaN;
    end

    % 3. Гіра r=2
    try
        tic;
        %[x0_3, ~, k_3, ~, ~] = gear_test(f, Jac_f, 0, T, h, x0_guess, eps_tol, tol, n);
        [x0_3, xT_3, k_3, err_T_3, err_x0_3, traj_3] = gear_test(f, Jac_f, 0, T, h, x0_guess, eps_tol, tol, n);
        times(3, i) = toc; errors(3, i) = norm(x0_3 - x0_ref); iters(3, i) = k_3;
    catch
        fprintf('  Гір r=2 розійшовся!\n');
        times(3, i) = NaN; errors(3, i) = NaN; iters(3, i) = NaN;
    end

    % 4. Гіра r=3
    try
        tic;
        %[x0_4, ~, k_4, ~, ~] = gear_r3(f, Jac_f, 0, T, h, x0_guess, eps_tol, tol, n);
        [x0_4, xT_4, k_4, err_T_4, err_x0_4, traj_4] = gear_r3(f, Jac_f, 0, T, h, x0_guess, eps_tol, tol, n);
        times(4, i) = toc; errors(4, i) = norm(x0_4 - x0_ref); iters(4, i) = k_4;
    catch ME
        fprintf('  Гір r=3 видав помилку: %s\n', ME.message);

        times(4, i) = NaN; errors(4, i) = NaN; iters(4, i) = NaN;
        x0_4 = [NaN; NaN];
        xT_4 = [NaN; NaN];
        traj_4 = NaN(2, length(0:h:T));
    end

    if n == 150
        fprintf('\n> Знайдені точки x0*:\n');
        fprintf('  RK4:           [%.8f; %.8f]\n', x0_1(1), x0_1(2));
        fprintf('  Неявний Ейлер: [%.8f; %.8f]\n', x0_2(1), x0_2(2));
        fprintf('  Гір r=2:       [%.8f; %.8f]\n', x0_3(1), x0_3(2));
        fprintf('  Гір r=3:       [%.8f; %.8f]\n', x0_4(1), x0_4(2));
        fprintf('  (Еталон:       [%.8f; %.8f])\n\n', x0_ref(1), x0_ref(2));
    end
end

%% Фазові портрети
fig1 = figure('Name', 'Фазові портрети', 'Position', [100, 100, 1000, 350]);

methods_names = {'РК4', 'Неявний Ейлера','Гіра r=2','Гірa r=3'};
trajs = {traj_1, traj_2, traj_3, traj_4};
x0_vals = {x0_1, x0_2, x0_3, x0_4};
xT_vals = {xT_1, xT_2, xT_3, xT_4};
line_colors = {'g', 'b', 'm', 'r'};

for m = 1:4
  subplot(1, 4, m);
  plot(trajs{m}(1,:), trajs{m}(2,:), 'Color', line_colors{m}, 'LineWidth', 1.5); hold on;
  plot(x0_vals{m}(1), x0_vals{m}(2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
  plot(xT_vals{m}(1), xT_vals{m}(2), 'rx', 'LineWidth', 2, 'MarkerSize', 8);

  err_norm = norm(xT_vals{m} - x0_vals{m});
  title(sprintf('%s\n||x(T)-x_0^*||=%.1e', methods_names{m}, err_norm), 'FontSize', 10);
  xlabel('x_1'); ylabel('x_2');
  grid on;%axis equal;

##  legend('орбіта', 'x_0^*', 'x(T)', 'Location', 'best');
##  xlim([-1.2 1.2]); ylim([-2.5 2.5]); %
  axis tight; axis square;
  margins = 0.05;
  xl = xlim(); yl = ylim();
  xlim([xl(1)-margins*diff(xl), xl(2)+margins*diff(xl)]);
  ylim([yl(1)-margins*diff(yl), yl(2)+margins*diff(yl)]);
  end

  %% x1(t), x2(t)
  fig2 = figure('Name', 'Часові залежності', 'Position', [150, 150, 800, 400]);
  t_span = 0:h:T;

  % x1(t)
  subplot(1, 2, 1);
  plot(t_span, traj_1(1,:), 'g-', 'LineWidth', 1); hold on;
  plot(t_span, traj_2(1,:), 'b--', 'LineWidth', 1.5);
  plot(t_span, traj_3(1,:), 'm:', 'LineWidth', 1.5);
  plot(t_span, traj_4(1,:), 'r-.', 'LineWidth', 1.5);
  title('x_1(t)'); xlabel('t'); grid on;
%  legend(methods_names, 'Location', 'northeast');

  % x2(t)
  subplot(1, 2, 2);
  plot(t_span, traj_1(2,:), 'g-', 'LineWidth', 1); hold on;
  plot(t_span, traj_2(2,:), 'b--', 'LineWidth', 1.5);
  plot(t_span, traj_3(2,:), 'm:', 'LineWidth', 1.5);
  plot(t_span, traj_4(2,:), 'r-.', 'LineWidth', 1.5);
  title('x_2(t)'); xlabel('t'); grid on;
##  legend(methods_names, 'Location', 'best', 'FontSize', 6);

  %% convergence
  fig3 = figure('Name', 'Збіжність Ньютона', 'Position', [200, 200, 1000, 300]);

  err_T_all = {err_T_1, err_T_2, err_T_3, err_T_4};
  err_x0_all = {err_x0_1, err_x0_2, err_x0_3, err_x0_4};

  for m = 1:4
    subplot(1, 4, m);
    iters_k = 1:length(err_T_all{m});
    semilogy(iters_k, err_T_all{m}, '-ob', 'LineWidth', 1, 'MarkerFaceColor', 'none'); hold on;
    semilogy(iters_k, err_x0_all{m}, '-sr', 'LineWidth', 1, 'MarkerFaceColor', 'none');

    title(methods_names{m}, 'FontSize', 10);
    xlabel('Ітерація'); ylabel('Похибка');
    grid on;

    set(gca, 'YMinorGrid', 'off');
    axis square;
%    legend('||x(T)-x_0^k||', '||x_0^{k+1}-x_0^k||', 'Location', 'southwest');
end


%% --- ГРАФІК 4:
n_coarse = 150;
h_coarse = T / n_coarse;

[~, ~, ~, ~, ~, traj_1_c] = rk4(f, Jac_f, 0, T, h_coarse, x0_guess, eps_tol, tol, n_coarse);
[~, ~, ~, ~, ~, traj_2_c] = euler(f, Jac_f, 0, T, h_coarse, x0_guess, eps_tol, tol, n_coarse);
[~, ~, ~, ~, ~, traj_3_c] = gear_test(f, Jac_f, 0, T, h_coarse, x0_guess, eps_tol, tol, n_coarse);
[~, ~, ~, ~, ~, traj_4_c] = gear_r3(f, Jac_f, 0, T, h_coarse, x0_guess, eps_tol, tol, n_coarse);

fig4 = figure('Position', [250, 250, 500, 500]);

plot(traj_ref(1,:), traj_ref(2,:), '-k', 'LineWidth', 3); hold on;
plot(traj_1_c(1,:), traj_1_c(2,:), '-g', 'LineWidth', 1.5);
plot(traj_2_c(1,:), traj_2_c(2,:), '--b', 'LineWidth', 2.5);
plot(traj_3_c(1,:), traj_3_c(2,:), '--m', 'LineWidth', 1.5);
plot(traj_4_c(1,:), traj_4_c(2,:), ':r', 'LineWidth', 2);

plot(x0_ref(1), x0_ref(2), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 8);

title(sprintf('Візуалізація на грубій сітці (n = %d)', n_coarse), 'FontSize', 12, 'FontWeight', 'bold');
xlabel('x_1(t)', 'FontSize', 11);
ylabel('x_2(t)', 'FontSize', 11);
%legend(sprintf('Еталонний розв''язок (n=%d)', n_ref), sprintf('РК4 (n=%d)', n_coarse), sprintf('Неявний Ейлер (n=%d)', n_coarse), sprintf('Гір r=2 (n=%d)', n_coarse), sprintf('Гір r=3 (n=%d)', n_coarse), 'Location', 'northeast', 'FontSize', 10);
axis equal;
grid on;


margins = 0.1; xl = xlim(); yl = ylim();
xlim([xl(1)-margins*diff(xl), xl(2)+margins*diff(xl)]);
ylim([yl(1)-margins*diff(yl), yl(2)+margins*diff(yl)]);


periods = 3;
T_multi = periods * T;
n_multi = periods * 500;
h_multi = T_multi / n_multi;
t_span_multi = 0:h_multi:T_multi;

[~, ~, ~, ~, ~, traj_multi] = rk4(f, Jac_f, 0, T_multi, h_multi, x0_ref, eps_tol, tol, n_multi);

fig5 = figure('Name', 'Періодичність розв''язку', 'Position', [300, 300, 800, 400]);

plot(t_span_multi, traj_multi(1,:), '-b', 'LineWidth', 1.5); hold on;

yl = ylim();
for p = 1:periods
    plot([p*T, p*T], yl, '--k', 'LineWidth', 1.2);
end

title('Періодичність розв''язку на інтервалі [0, 3T]', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Час t', 'FontSize', 11);
ylabel('Координата x_1(t)', 'FontSize', 11);
grid on;

xticks(0:T:T_multi);
xticklabels({'0', 'T', '2T', '3T'});
ylim(yl);

fprintf('\n\n=== ЗВЕДЕНА ТАБЛИЦЯ РЕЗУЛЬТАТІВ ===\n');
fprintf('%-15s | %-5s | %-7s | %-8s | %-12s | %-7s\n', 'Метод', 'n', 'h', 'Ітерацій', 'Похибка', 'Час (с)');
fprintf(repmat('-', 1, 65)); fprintf('\n');

methods_str = {'Явний РК4', 'Неявний Ейлер', 'Метод Гіра r=2', 'Метод Гіра r=3'};

for m = 1:4
    for i = 1:num_n
        if i == 1
            m_name = methods_str{m};
        else
            m_name = '';
        end

        fprintf('%-15s | %-5d | %-7.5f | %-8d | %-12.2e | %-7.2f\n', ...
            m_name, n_values(i), h_values(i), iters(m, i), errors(m, i), times(m, i));
    end
    if m < 4
        fprintf(repmat('-', 1, 65)); fprintf('\n');
    end
end
