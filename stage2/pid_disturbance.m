addpath("../")
clear all;
close all;

F_Cin=27; T_C=25; T_H=65; T_D=29;
F_C=21; F_H=17; F_D=13;
Tau_C=100; Tau=40; C=0.15; alpha=6;
Hpp=72.25; Tpp=39.35;
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-7);
tmax = 1550; 
Tp = 1;
iterations = 0:Tp:tmax;

K1z = -1;
K2z = 0.4035;
decoupling = true;
transfer_functions = get_transfer_func;
K = [dcgain(transfer_functions(1,1)) dcgain(transfer_functions(1,2)); dcgain(transfer_functions(2,1)) dcgain(transfer_functions(2,2))];
D12 = -K(1,2) / K(1,1);
D21 = -K(2, 1) / K(2,2);
M = [1, D12; D21, 1];
delay_d = round(Tau_C / Tp);
delay_steps = round(Tau / Tp);
delay_u1 = round(Tau_C / Tp);

y_zad = ones(length(iterations), 2);
y_zad(:, 1) = Hpp;
y_zad(:, 2) = Tpp;

disturbance = zeros(length(iterations), 1);
disturbance(200:end) = 1.5*F_D;

if(decoupling)
    [r2_h, r1_h, r0_h] = discrete_pid_params(0.1491 , 109.0822 , 0, Tp);
    [r2_t, r1_t, r0_t] = discrete_pid_params(0.5272, 156.1887, 0, Tp);
else
    [r2_h, r1_h, r0_h] = discrete_pid_params(0.2726 , 137.0822 , 0, Tp);
    [r2_t, r1_t, r0_t] = discrete_pid_params(0.5272, 156.1887, 0, Tp);
end

y_all = cell(1,2);
u_all = cell(1,2);

for mode = 1:2
    disturbance_comp = (mode == 2);
    H0 = Hpp; T0 = Tpp;
    u = zeros(length(iterations), 2);
    u(:, 1) = F_C; u(:, 2) = F_H;
    v = zeros(length(iterations), 2);
    if(decoupling)
        v_init = M \ [F_C; F_H]; 
        v(:, 1) = v_init(1); v(:, 2) = v_init(2);
    else
        v(:, 1) = F_C; v(:, 2) = F_H;
    end
    y = zeros(length(iterations), 2);
    y(:, 1) = Hpp; y(:, 2) = Tpp;
    T_out_without_delay = ones(length(iterations), 1)*Tpp;
    e = zeros(length(iterations), 2);

    for k=3:length(iterations)
        if k > delay_u1
            FC_u = u(k - delay_u1, 1); 
        else
            FC_u = F_C;
        end
        FH_u = u(k-1, 2);
        FD_u = F_D + disturbance(k);
        
        t_span = ((k-1)*Tp):(k*Tp);
        [~, S]  = ode15s(@(t, x) nonlinear_model(x,t,FC_u,T_C,T_H,T_D,F_C,FH_u,FD_u,Tau_C,C,alpha), t_span, [H0; T0], options);
        y(k, 1)  = S(end, 1);
        T_out_without_delay(k)  = S(end, 2);
        if(k > delay_steps)
            y(k, 2) = T_out_without_delay(k-delay_steps);
        else 
            y(k, 2) = Tpp;
        end
        H0 = y(k, 1); T0 = T_out_without_delay(k);

        e(k, 1) = y_zad(k, 1) - y(k, 1);
        v(k, 1) = r2_h*e(k-2, 1) + r1_h*e(k-1, 1) + r0_h*e(k, 1) + v(k-1, 1);
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

        if(disturbance_comp)
            u(k, 1) = u(k, 1);
            u(k, 2) = u(k, 2) + K2z*disturbance(k);
        end

        u(k, 1) = min(max(u(k, 1), 0), 300);
        u(k, 2) = min(max(u(k, 2), 0), 300);
    end
    y_all{mode} = y;
    u_all{mode} = u;
end

figure('Name', 'Sygnały sterujące');
plot(iterations, u_all{1}(:, 1), 'b:', 'LineWidth', 1); hold on;
plot(iterations, u_all{2}(:, 1), 'b-', 'LineWidth', 1.2);
plot(iterations, u_all{1}(:, 2), 'r:', 'LineWidth', 1);
plot(iterations, u_all{2}(:, 2), 'r-', 'LineWidth', 1.2);
plot(iterations, disturbance + F_D, 'k--', 'LineWidth', 1.2);
xlabel('Iteracje [k]'); ylabel('Przepływ [cm^3/s]');
legend("FC (bez komp)", "FC (z komp)", "FH (bez komp)", "FH (z komp)", "FD");
grid on; 
set(gcf, 'position', [10, 10, 1000, 700]);
exportgraphics(gcf, "images/pid_control_u_disturbance.pdf", "ContentType", "vector");

figure('Name', 'Wyjście 1: Wysokość');
plot(iterations, y_all{1}(:, 1), 'LineWidth', 1.2); hold on;
plot(iterations, y_all{2}(:, 1), 'LineWidth', 1.5);
plot(iterations, y_zad(:, 1), 'k--', 'LineWidth', 1);
xlabel('Iteracje [k]'); ylabel('Wysokość [cm]');
legend("h", "h komp", "h_{zad}");
grid on; 
set(gcf, 'position', [10, 10, 1000, 700]);
exportgraphics(gcf, "images/pid_control_h_disturbance.pdf", "ContentType", "vector");

figure('Name', 'Wyjście 2: Temperatura');
plot(iterations, y_all{1}(:, 2), 'LineWidth', 1.2); hold on;
plot(iterations, y_all{2}(:, 2), 'LineWidth', 1.5);
plot(iterations, y_zad(:, 2), 'r--', 'LineWidth', 1);
xlabel('Iteracje [k]'); ylabel('Temperatura [°C]');
legend("T_{out}", "T_{out} komp", "T_{zad}");
grid on; 
set(gcf, 'position', [10, 10, 1000, 700]);
exportgraphics(gcf, "images/pid_control_t_disturbance.pdf", "ContentType", "vector");