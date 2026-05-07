addpath("../");
clear all;
close all;

lb = [0.001, 5, 0.01, 10]; 
ub = [50.0,  1000, 50.0,  1000];  

nVars = 4;
options_pso = optimoptions('particleswarm', ...
    'SwarmSize', 40, ...    % Wielkość roju (liczba cząstek)
    'MaxIterations', 20, ... % Liczba epok
    'Display', 'iter', ...
    'PlotFcn', 'pswplotbestf');

% Uruchomienie PSO
[best_params, min_error] = particleswarm(@(x) objective_function(x), nVars, lb, ub, options_pso);

[~, y, y_zad, u, iterations] = objective_function(best_params);

figure('Name', 'Wyniki po optymalizacji PSO');
subplot(2,1,1);
stairs(iterations, y(:, 1), 'LineWidth', 1.5); hold on;
stairs(iterations, y_zad(:, 1), '--k');
stairs(iterations, y(:, 2), 'LineWidth', 1.5);
stairs(iterations, y_zad(:, 2), '--r');
grid on; legend("H", "H zad", "T", "T zad"); title('Wyjścia obiektu');

subplot(2,1,2);
stairs(iterations, u(:, 1)); hold on;
stairs(iterations, u(:, 2));
grid on; legend("FC", "FH"); title('Sygnały sterujące (Saturacja 0-100)');

function [total_error, y, y_zad, u, iterations] = objective_function(x)
    K_h = x(1); Ti_h = x(2); Td_h = 0;
    K_t = x(3); Ti_t = x(4); Td_t = 0;
    
    F_Cin=27; T_C=25; T_H=65; T_D=29; F_C=21; F_H=17; F_D=13;
    Tau_C=100; Tau=40; C=0.15; alpha=6;
    Hpp=72.25; Tpp=39.35; Tp = 1; tmax = 7900;
    iterations = 0:Tp:tmax;
    
    % Inicjalizacja
    y = ones(length(iterations), 2) .* [Hpp, Tpp];
    T_out_without_delay = ones(length(iterations), 1) * Tpp;
    e = zeros(length(iterations), 2);
    y_zad = ones(length(iterations), 2) .* [Hpp, Tpp];

    % Decoupling
    decoupling = true;
    delay_d = round(100 / Tp);
    delay_u1 = round(Tau_C / Tp);
    transfer_functions = get_transfer_func;
    K = [dcgain(transfer_functions(1,1)) dcgain(transfer_functions(1,2)); dcgain(transfer_functions(2,1)) dcgain(transfer_functions(2,2))];
    D12 = -K(1,2) / K(1,1);
    D21 = -K(2, 1) / K(2,2);
    M = [1, D12; D21, 1];
    
    u = zeros(length(iterations), 2);
    u(:, 1) = F_C;
    u(:, 2) = F_H;
    
    v = zeros(length(iterations), 2);
    if(decoupling)
        v_init = M \ [F_C; F_H]; 
        v(:, 1) = v_init(1);
        v(:, 2) = v_init(2);
    else
        v(:, 1) = F_C;
        v(:, 2) = F_H;
    end
   
   % Jumps - height
    y_zad(500:end, 1) = 1.4*Hpp;
    y_zad(3500:end, 1) = 0.7*Hpp;
    y_zad(6000:end, 1) = 1.0*Hpp;
    
    % Jumps - temperature
    y_zad(2000:end, 2) = 0.8*Tpp;
    y_zad(5000:end, 2) = 1.2*Tpp;

    % Obliczenie parametrów dyskretnych regulatora
    [r2_h, r1_h, r0_h] = discrete_pid_params(K_h, Ti_h, Td_h, Tp);
    [r2_t, r1_t, r0_t] = discrete_pid_params(K_t, Ti_t, Td_t, Tp);

    H0 = Hpp; T0 = Tpp;
    delay_steps = round(Tau / Tp);
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-7);

    try
        for k=3:length(iterations)
            % Pobranie sterowania z poprzedniego kroku
            if k > delay_u1
                FC_u = u(k - delay_u1, 1); 
            else
                FC_u = F_C;
            end
            FH_u = u(k-1, 2);
            
            % Symulacja modelu
            t_span = ((k-1)*Tp):(k*Tp);
            [~, S] = ode15s(@(t, x) nonlinear_model(x,t,FC_u,T_C,T_H,T_D,F_C,FH_u,F_D,Tau_C,C,alpha), t_span, [H0; T0], options);
            
            y(k, 1) = S(end, 1);
            T_out_without_delay(k) = S(end, 2);
            
            % Opóźnienie
            if k > delay_steps
                y(k, 2) = T_out_without_delay(k-delay_steps);
            else
                y(k, 2) = Tpp;
            end
            
            H0 = y(k, 1); T0 = T_out_without_delay(k);

            % Height regulator
            e(k, 1) = y_zad(k, 1) - y(k, 1);
            v(k, 1) = r2_h*e(k-2, 1) + r1_h*e(k-1, 1) + r0_h*e(k, 1) + v(k-1, 1);
            
            % Temperature regulator
            e(k, 2) = y_zad(k, 2) - y(k, 2);
            v(k, 2) = r2_t*e(k-2, 2) + r1_t*e(k-1, 2) + r0_t*e(k, 2) + v(k-1, 2);
        
            v(k, 1) = min(max(v(k, 1), -200), 300);
            v(k, 2) = min(max(v(k, 2), -200), 300);

            if k > delay_d
                v1_delayed = v(k - delay_d, 1);
            else
                v1_delayed = v(1, 1); 
            end
        
            if(decoupling)
                u(k, 1) = v(k, 1) + D12 * v(k, 2); 
                u(k, 2) = v(k, 2) + D21 * v1_delayed; 
            else
                u(k, :) = v(k, :);
            end
        
            u(k, 1) = min(max(u(k, 1), 0), 300);
            u(k, 2) = min(max(u(k, 2), 0), 300);
                    
            if y(k, 1) <= 0
                total_error = 1e12; return;
            end
        end
        
        error_h = sum((e(:, 1) / Hpp).^2);
        error_t = sum((e(:, 2) / Tpp).^2);
        total_error = 1.0*error_h + 1.0*error_t;

    catch
        total_error = 1e12;
    end
end