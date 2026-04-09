clear all;
close all;

F_Cin=27; T_C=25; T_H=65; T_D=29;
F_C=21; F_H=17; F_D=13;
Tau_C=100; Tau=40; C=0.15; alpha=6;
H0=72.25; T0=39.35;

options = odeset('RelTol', 1e-9, 'AbsTol', 1e-7);
tmax = 1000; 
Tp = 1;
Ts = 0:Tp:tmax;

FH_jumps = [30, 20, 10, -10, -15]; 
FC_jumps = [15, 10, 5, -5, -10, -15]; 
FD_jumps = [10, 5, 2, -5, -10];

h_tf_FH = zeros(length(Ts), length(FH_jumps));
h_l_FH = zeros(length(Ts), length(FH_jumps));
T_out_tf_FH = zeros(length(Ts), length(FH_jumps));
T_out_l_FH = zeros(length(Ts), length(FH_jumps));

h_tf_FC = zeros(length(Ts), length(FC_jumps));
h_l_FC = zeros(length(Ts), length(FC_jumps));
T_out_tf_FC = zeros(length(Ts), length(FC_jumps));
T_out_l_FC = zeros(length(Ts), length(FC_jumps));

h_tf_FD = zeros(length(Ts), length(FD_jumps));
h_l_FD = zeros(length(Ts), length(FD_jumps));
T_out_tf_FD = zeros(length(Ts), length(FD_jumps));
T_out_l_FD = zeros(length(Ts), length(FD_jumps));

legend_str_FH = {};
legend_str_FC = {};
legend_str_FD = {};
colors = lines(6);
delay_steps = round(Tau / Tp);

G = get_transfer_func();

for k = 1:length(FH_jumps)
    FH_u = F_H + FH_jumps(k); 
    FC_u = F_C;                   
    FD_u = F_D;
    
    [~, S_l] = ode15s(@(t, x) linear_model(x,t,FC_u,T_C,T_H,T_D,F_C,FH_u,FD_u,Tau_C,C,alpha), Ts, [H0; T0], options);
    
    h_l_FH(:, k)  = S_l(:, 1);
    T_l  = S_l(:, 2);
    
    for i = 1:length(Ts)
        idx = max(1, i - delay_steps);
        T_out_l_FH(i, k)  = T_l(idx);
    end
    
    U = zeros(length(Ts), 3);
    U(:, 2) = FH_jumps(k);    
    
    Y_tf = lsim(G, U, Ts);    
    
    h_tf_FH(:, k) = Y_tf(:, 1) + H0;      
    T_out_tf_FH(:, k) = Y_tf(:, 2) + T0;  
    
    legend_str_FH{end+1} = sprintf('Skok %d, transmitancja', FH_jumps(k));
    legend_str_FH{end+1} = sprintf('Skok %d, liniowy ODE', FH_jumps(k));
end

for k = 1:length(FC_jumps)
    FH_u = F_H;
    FC_u = F_C + FC_jumps(k);                   
    FD_u = F_D;
    
    [~, S_l] = ode15s(@(t, x) linear_model(x,t,FC_u,T_C,T_H,T_D,F_C,FH_u,FD_u,Tau_C,C,alpha), Ts, [H0; T0], options);
    
    h_l_FC(:, k)  = S_l(:, 1);
    T_l  = S_l(:, 2);
    
    for i = 1:length(Ts)
        idx = max(1, i - delay_steps);
        T_out_l_FC(i, k)  = T_l(idx);
    end
    
    U = zeros(length(Ts), 3);
    U(:, 1) = FC_jumps(k); 
    
    Y_tf = lsim(G, U, Ts);
    
    h_tf_FC(:, k) = Y_tf(:, 1) + H0;
    T_out_tf_FC(:, k) = Y_tf(:, 2) + T0;
    
    legend_str_FC{end+1} = sprintf('Skok %d, transmitancja', FC_jumps(k));
    legend_str_FC{end+1} = sprintf('Skok %d, liniowy ODE', FC_jumps(k));
end

for k = 1:length(FD_jumps)
    FH_u = F_H;
    FC_u = F_C;                   
    FD_u = F_D + FD_jumps(k);
    
    [~, S_l] = ode15s(@(t, x) linear_model(x,t,FC_u,T_C,T_H,T_D,F_C,FH_u,FD_u,Tau_C,C,alpha), Ts, [H0; T0], options);
    
    h_l_FD(:, k)  = S_l(:, 1);
    T_l  = S_l(:, 2);
    
    for i = 1:length(Ts)
        idx = max(1, i - delay_steps);
        T_out_l_FD(i, k)  = T_l(idx);
    end
    
    U = zeros(length(Ts), 3);
    U(:, 3) = FD_jumps(k); 
    
    Y_tf = lsim(G, U, Ts);
    
    h_tf_FD(:, k) = Y_tf(:, 1) + H0;
    T_out_tf_FD(:, k) = Y_tf(:, 2) + T0;
    
    legend_str_FD{end+1} = sprintf('Skok %d, transmitancja', FD_jumps(k));
    legend_str_FD{end+1} = sprintf('Skok %d, liniowy ODE', FD_jumps(k));
end

figure(1); hold on; grid on;
for k = 1:length(FH_jumps)
    plot(Ts, h_tf_FH(:, k), 'Color', colors(k, :), 'LineStyle', '--', 'LineWidth', 2.0);
    plot(Ts, h_l_FH(:, k), 'Color', colors(k, :), 'LineStyle', '-', 'LineWidth', 0.5);
end
set(gcf, 'Color', 'w');
ax1 = gca;
set(ax1, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.8 0.8 0.8]);
xlabel('Czas [s]', 'Interpreter', 'latex');
ylabel('$h$ [cm]', 'Interpreter', 'latex');
lgd1 = legend(legend_str_FH, 'Location', 'eastoutside', 'TextColor', 'k');
title(lgd1, 'Skoki $F_H$', 'Interpreter', 'latex');
set(gcf, 'position', [10, 10, 1000, 700]);
exportgraphics(gcf, "images/multi_jump_FH_h_TF.pdf", "ContentType", "vector");

figure(2); hold on; grid on;
for k = 1:length(FH_jumps)
    plot(Ts, T_out_tf_FH(:, k), 'Color', colors(k, :), 'LineStyle', '--', 'LineWidth', 2.0);
    plot(Ts, T_out_l_FH(:, k), 'Color', colors(k, :), 'LineStyle', '-', 'LineWidth', 0.5);
end
set(gcf, 'Color', 'w');
ax2 = gca;
set(ax2, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.8 0.8 0.8]);
xlabel('Czas [s]', 'Interpreter', 'latex');
ylabel('$T_{out}$ [$^\circ$C]', 'Interpreter', 'latex');
lgd2 = legend(legend_str_FH, 'Location', 'eastoutside', 'TextColor', 'k');
title(lgd2, 'Skoki $F_H$', 'Interpreter', 'latex');
set(gcf, 'position', [10, 10, 1000, 700]);
exportgraphics(gcf, "images/multi_jump_FH_T_TF.pdf", "ContentType", "vector");

figure(3); hold on; grid on;
for k = 1:length(FC_jumps)
    plot(Ts, h_tf_FC(:, k), 'Color', colors(k, :), 'LineStyle', '--', 'LineWidth', 2.0);
    plot(Ts, h_l_FC(:, k), 'Color', colors(k, :), 'LineStyle', '-', 'LineWidth', 0.5);
end
set(gcf, 'Color', 'w');
ax3 = gca;
set(ax3, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.8 0.8 0.8]);
xlabel('Czas [s]', 'Interpreter', 'latex');
ylabel('$h$ [cm]', 'Interpreter', 'latex');
lgd3 = legend(legend_str_FC, 'Location', 'eastoutside', 'TextColor', 'k');
title(lgd3, 'Skoki $F_C$', 'Interpreter', 'latex');
set(gcf, 'position', [10, 10, 1000, 700]);
exportgraphics(gcf, "images/multi_jump_FC_h_TF.pdf", "ContentType", "vector");

figure(4); hold on; grid on;
for k = 1:length(FC_jumps)
    plot(Ts, T_out_tf_FC(:, k), 'Color', colors(k, :), 'LineStyle', '--', 'LineWidth', 2.0);
    plot(Ts, T_out_l_FC(:, k), 'Color', colors(k, :), 'LineStyle', '-', 'LineWidth', 0.5);
end
set(gcf, 'Color', 'w');
ax4 = gca;
set(ax4, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.8 0.8 0.8]);
xlabel('Czas [s]', 'Interpreter', 'latex');
ylabel('$T_{out}$ [$^\circ$C]', 'Interpreter', 'latex');
lgd4 = legend(legend_str_FC, 'Location', 'eastoutside', 'TextColor', 'k');
title(lgd4, 'Skoki $F_C$', 'Interpreter', 'latex');
set(gcf, 'position', [10, 10, 1000, 700]);
exportgraphics(gcf, "images/multi_jump_FC_T_TF.pdf", "ContentType", "vector");

figure(5); hold on; grid on;
for k = 1:length(FD_jumps)
    plot(Ts, h_tf_FD(:, k), 'Color', colors(k, :), 'LineStyle', '--', 'LineWidth', 2.0);
    plot(Ts, h_l_FD(:, k), 'Color', colors(k, :), 'LineStyle', '-', 'LineWidth', 0.5);
end
set(gcf, 'Color', 'w');
ax5 = gca;
set(ax5, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.8 0.8 0.8]);
xlabel('Czas [s]', 'Interpreter', 'latex');
ylabel('$h$ [cm]', 'Interpreter', 'latex');
lgd5 = legend(legend_str_FD, 'Location', 'eastoutside', 'TextColor', 'k');
title(lgd5, 'Skoki $F_D$', 'Interpreter', 'latex');
set(gcf, 'position', [10, 10, 1000, 700]);
exportgraphics(gcf, "images/multi_jump_FD_h_TF.pdf", "ContentType", "vector");

figure(6); hold on; grid on;
for k = 1:length(FD_jumps)
    plot(Ts, T_out_tf_FD(:, k), 'Color', colors(k, :), 'LineStyle', '--', 'LineWidth', 2.0);
    plot(Ts, T_out_l_FD(:, k), 'Color', colors(k, :), 'LineStyle', '-', 'LineWidth', 0.5);
end
set(gcf, 'Color', 'w');
ax6 = gca;
set(ax6, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.8 0.8 0.8]);
xlabel('Czas [s]', 'Interpreter', 'latex');
ylabel('$T_{out}$ [$^\circ$C]', 'Interpreter', 'latex');
lgd6 = legend(legend_str_FD, 'Location', 'eastoutside', 'TextColor', 'k');
title(lgd6, 'Skoki $F_D$', 'Interpreter', 'latex');
set(gcf, 'position', [10, 10, 1000, 700]);
exportgraphics(gcf, "images/multi_jump_FD_T_TF.pdf", "ContentType", "vector");
