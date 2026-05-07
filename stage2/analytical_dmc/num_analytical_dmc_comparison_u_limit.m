addpath("../")
clear all;
close all;

%% Parametry obiektu i punktu pracy
F_Cin = 27; T_C = 25; T_H = 65; T_D = 29;
F_C = 21; F_H = 17; F_D0 = 13;
Tau_C = 100; Tau = 40; C = 0.15; alpha = 6;
Hpp = 72.25; Tpp = 39.35;

%% Parametry symulacji
tmax = 3500; 
Tp = 1;
start_k = 2;
N_steps = round(tmax/Tp) + 1;
t_span_total = 0:Tp:tmax;
delay_steps = round(Tau / Tp);
delay_u1 = round(Tau_C / Tp);

%% Nastawy regulatorów DMC
D = 750;      
N = 170;             
Nu = 3;             
lambda = [10.0; 5.0]; % [lambda_FC; lambda_FH]
psi = [1.0; 40.0];    % [psi_H; psi_T]

u_min = [0; 0];
u_max = [100; 100]; 
du_min = [-2; -2];
du_max = [2; 2];
y_min = [0; 0];      
y_max = [200; 100];  

%% Obliczenia Offline
[S_data, Sz_data] = get_s(D+1, F_C, F_H, F_D0, T_C, T_H, T_D, Tau_C, C, alpha, Hpp, Tpp, Tp, delay_u1, delay_steps);
[M, MP, MPz] = DMC_offline(S_data, Sz_data, D, N, Nu);

Psi_mat = diag(repmat(psi, N, 1));
Lambda_mat = diag(repmat(lambda, Nu, 1));

H_qp = 2 * (M' * Psi_mat * M + Lambda_mat);
H_qp = (H_qp + H_qp') / 2; % Symetryzacja
K_full = H_qp \ (2 * M' * Psi_mat);
K_ana = K_full(1:2, :); 

%% Wartości zadane
y_zad = ones(N_steps, 2) .* [Hpp, Tpp];
y_zad(500:end, 1) = 1.5*Hpp;
y_zad(500:end, 2) = 1.4*Tpp;

F_D_seq = F_D0 * ones(N_steps, 1);
F_D_seq(1750:end) = 1.2 * F_D0; 

ode_opt = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);


if isempty(gcp('nocreate'))
    parpool; 
end

%% Symulacja
y_results = cell(2, 1);
u_results = cell(2, 1);
err_results = zeros(2, 1);

parfor method = 1:2
    y_local = ones(N_steps, 2) .* [Hpp, Tpp];
    u_local = ones(N_steps, 2) .* [F_C, F_H];
    dup_local = zeros((D-1) * 2, 1);
    
    dzp_local = zeros(D-1, 1); 
    FD_prev = F_D0;
    
    T_nd_local = ones(N_steps, 1) * Tpp;
    H0_local = Hpp; T0_local = Tpp;
    
    for k = start_k:N_steps-1
        if mod(k, 10) == 0 
            metoda_nazwa = "Numeryczny"; if method == 2, metoda_nazwa = "Analityczny"; end
            fprintf('[%s] Iteracja %d / %d\n', metoda_nazwa, k, N_steps);
        end
        
        FC_u = F_C;
        if k > delay_u1, FC_u = u_local(k - delay_u1, 1); end
        
        t_span = [(k-1)*Tp, k*Tp];
        FD_curr = F_D_seq(k);
        
        [~, X_out] = ode15s(@(t, x) nonlinear_model(x, t, FC_u, T_C, T_H, T_D, F_C, u_local(k,2), FD_curr, Tau_C, C, alpha), ...
                        t_span, [H0_local; T0_local], ode_opt);
                    
        y_local(k, 1) = X_out(end, 1);
        T_nd_local(k) = X_out(end, 2);
        if k > delay_steps, y_local(k, 2) = T_nd_local(k - delay_steps); else, y_local(k, 2) = Tpp; end
        
        H0_local = y_local(k, 1); T0_local = T_nd_local(k);
        
        y_curr = [y_local(k, 1); y_local(k, 2)];
        y_zad_curr = [y_zad(k, 1); y_zad(k, 2)];
        u_prev = [u_local(k, 1); u_local(k, 2)];
        
        dz_curr = FD_curr - FD_prev;
        FD_prev = FD_curr;
        
        y_free = repmat(y_curr, N, 1) + MP * dup_local + MPz * dzp_local;
        
        % DMC numeryczny
        if method == 1
            [du1_n, du2_n] = numDMC_online(H_qp, M, Psi_mat, y_free, y_zad_curr, ...
                                            u_prev, u_min, u_max, y_min, y_max, du_min, du_max, N, Nu);
            
            u_local(k+1, 1) = u_local(k, 1) + du1_n;
            u_local(k+1, 2) = u_local(k, 2) + du2_n;
            dup_local = [du1_n; du2_n; dup_local(1:end-2)];
        % DMC analityczny
        else
            Y_zad_full = repmat(y_zad_curr, N, 1);
            
            du_step = K_ana * (Y_zad_full - y_free);
            
            du_step = min(max(du_step, du_min), du_max);
            u_temp = u_prev + du_step;
            u_temp = min(max(u_temp, u_min), u_max);
            
            du_realized = u_temp - u_prev;
            
            u_local(k+1, :) = u_temp';
            dup_local = [du_realized; dup_local(1:end-2)];
        end
        
        dzp_local = [dz_curr; dzp_local(1:end-1)];
    end
    
    y_local(N_steps, :) = y_local(N_steps-1, :);
    
    err_results(method) = sum(((y_zad - y_local) ./ [Hpp, Tpp]).^2, 'all');
    
    y_results{method} = y_local;
    u_results{method} = u_local;
end

y_num = y_results{1}; u_num = u_results{1}; err_num = err_results(1);
y_ana = y_results{2}; u_ana = u_results{2}; err_ana = err_results(2);
disp('Symulacja zakończona.');
fprintf('Błąd całkowity J - Numeryczny:  %.4f\n', err_num);
fprintf('Błąd całkowity J - Analityczny: %.4f\n', err_ana);

%% Wykresy

figure('Name', 'DMC: Wyjścia obiektu');
yyaxis left
stairs(t_span_total, y_num(:, 1), 'b-', 'LineWidth', 1.5); hold on;
stairs(t_span_total, y_ana(:, 1), 'r-.', 'LineWidth', 1.5);
stairs(t_span_total, y_zad(:, 1), 'k--', 'LineWidth', 1);
ylabel('Wysokość [cm]');
ylim([min(y_zad(:,1))*0.75, max(y_zad(:,1))*1.4]); 

yyaxis right
stairs(t_span_total, y_num(:, 2), 'c-', 'LineWidth', 1.5); hold on;
stairs(t_span_total, y_ana(:, 2), 'm-.', 'LineWidth', 1.5);
stairs(t_span_total, y_zad(:, 2), 'k--', 'LineWidth', 1);
ylabel('Temperatura [°C]');
ylim([min(y_zad(:,2))*0.75, max(y_zad(:,2))*1.4]);

xlabel('Iteracje [k]');
legend("H (Num)", "H (Ana)", "H_{zad}", "T_{out} (Num)", "T_{out} (Ana)", "T_{zad}", 'Location', 'best');
grid on;
set(gcf, 'position', [10, 10, 1000, 700]);
exportgraphics(gcf, "images/dmc_control_y_input_delta_limit.pdf", "ContentType", "vector");

figure('Name', 'DMC: Sygnały sterujące');
stairs(t_span_total, u_num(:, 1), 'b-', 'LineWidth', 1.2); hold on;
stairs(t_span_total, u_ana(:, 1), 'r-.', 'LineWidth', 1.2);
stairs(t_span_total, u_num(:, 2), 'c-', 'LineWidth', 1.2);
stairs(t_span_total, u_ana(:, 2), 'm-.', 'LineWidth', 1.2);
stairs(t_span_total, F_D_seq, 'k', 'LineWidth', 1.5);
xlabel('Iteracje [k]');
ylabel('u [cm^3/s]');
legend("FC (Num)", "FC (Ana)", "FH (Num)", "FH (Ana)", "FD", 'Location', 'best');
grid on;
set(gcf, 'position', [10, 10, 1000, 700]);
exportgraphics(gcf, "images/dmc_control_u_input_delta_limit.pdf", "ContentType", "vector");
