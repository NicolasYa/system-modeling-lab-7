clear; clc; close all;

%% Параметры варианта 14
lambda = 10;        % интенсивность входного потока
mu = 6;             % производительность канала
eps_dop = 0.01;     % допустимая абсолютная погрешность
alpha_d = 3;        % коэффициент доверительного интервала 

%% Теоретические значения
p_theor = mu / (lambda + mu);
q_theor = lambda / (lambda + mu);
fprintf('Исходные данные:\n');
fprintf('  lambda = %.2f s^-1\n', lambda);
fprintf('  mu     = %.2f s^-1\n', mu);
fprintf('  eps_dop = %.4f\n', eps_dop);
fprintf('\nТеоретические значения:\n');
fprintf('  p_theor = %.6f\n', p_theor);
fprintf('  q_theor = %.6f\n\n', q_theor);

x_theor = p_theor;   % искомая характеристика - вероятность обслуживания

%% Начальная серия наблюдений (T = 100 с)
T_obs = 100;
[N, M, ~] = simulate_QS(T_obs, lambda, mu);
p_est = M / N;
eps_est = alpha_d * sqrt(p_est * (1 - p_est) / N);

fprintf('Начальная серия наблюдений (T = %.1f s)\n', T_obs);
fprintf('  N   = %d\n', N);
fprintf('  M   = %d\n', M);
fprintf('  p*  = %.6f\n', p_est);
fprintf('  eps = %.6f\n', eps_est);
if eps_est <= eps_dop
    fprintf('  Точность достигнута (eps <= %.4f)\n', eps_dop);
else
    fprintf('  Точность НЕ достигнута (eps > %.4f)\n', eps_dop);
end
fprintf('\n');

% История для графиков
N_history = N;
p_history = p_est;
eps_history = eps_est;

%% Итерационное уточнение
iteration = 0;
while eps_est > eps_dop
    iteration = iteration + 1;
    fprintf('Итерация %d\n', iteration);
     
    % Оценка требуемого количества опытов
    N_req = ceil(alpha_d^2 * p_est * (1 - p_est) / eps_dop^2);
    N_add = max(N_req - N, ceil(0.1 * N));
    fprintf('  Требуемое N = %d, добавка = %d\n', N_req, N_add);
    
    % Оценка интенсивности по текущей серии и пересчёт в дополнительное время
    lambda_est = N / T_obs;
    T_add = N_add / lambda_est;
    fprintf('  Оценка lambda = %.4f s^-1\n', lambda_est);
    fprintf('  Дополнительное T = %.4f s\n', T_add);
    
    % Дополнительная серия
    [N_add_actual, M_add_actual, ~] = simulate_QS(T_add, lambda, mu);
    N = N + N_add_actual;
    M = M + M_add_actual;
    T_obs = T_obs + T_add;
    
    % Обновление оценки
    p_est = M / N;
    eps_est = alpha_d * sqrt(p_est * (1 - p_est) / N);
    
    fprintf('  Новые N = %d, M = %d\n', N, M);
    fprintf('  p*  = %.6f\n', p_est);
    fprintf('  eps = %.6f\n\n', eps_est);
    
    N_history(end+1) = N;
    p_history(end+1) = p_est;
    eps_history(end+1) = eps_est;
    
    if iteration > 50
        fprintf('Превышено максимальное число итераций\n');
        break;
    end
end

%% Финальные результаты
fprintf('Итоговые результаты\n');
fprintf('  Общее N = %d\n', N);
fprintf('  Общее T = %.2f s\n', T_obs);
fprintf('  p* = %.6f\n', p_est);
fprintf('  Теоретическое p = %.6f\n', p_theor);
fprintf('  Абсолютная погрешность = %.6f\n', abs(p_est - p_theor));
fprintf('  Достигнутая точность eps = %.6f\n', eps_est);
fprintf('  Требуемая точность eps_dop = %.4f\n', eps_dop);

%% График сходимости
figure('Name','Сходимость оценки');
subplot(2,1,1);
semilogx(N_history, p_history, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
yline(p_theor, 'r--', 'LineWidth', 2);
xlabel('N (число опытов)');
ylabel('p*');
title('Зависимость оценки p* от числа опытов');
legend('Оценка p*', 'Теоретическое p', 'Location', 'best');
grid on;

subplot(2,1,2);
semilogx(N_history, eps_history, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
yline(eps_dop, 'r--', 'LineWidth', 2);
xlabel('N (число опытов)');
ylabel('eps');
title('Зависимость погрешности от числа опытов');
legend('Оценка eps', 'eps_{доп}', 'Location', 'best');
grid on;

%% Функция имитационного моделирования одноканальной СМО с отказами
function [N, M, T_actual] = simulate_QS(T_max, lambda, mu)
    t = 0;
    N = 0;
    M = 0;
    busy = false;
    remaining_service = 0;
    
    while t < T_max
        dt = -1/lambda * log(1 - rand());
        
        if ~busy
            t_new = t + dt;
            if t_new > T_max, break; end
            t = t_new;
            N = N + 1;
            M = M + 1;
            busy = true;
            remaining_service = -1/mu * log(1 - rand());
        else
            if dt < remaining_service
                t_new = t + dt;
                if t_new > T_max, break; end
                t = t_new;
                N = N + 1;
                remaining_service = remaining_service - dt;
            else
                t_new = t + dt;
                if t_new > T_max, break; end
                t = t_new;
                N = N + 1;
                M = M + 1;
                busy = true;
                remaining_service = -1/mu * log(1 - rand());
            end
        end
    end
    T_actual = t;
end