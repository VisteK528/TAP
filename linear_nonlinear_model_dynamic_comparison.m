clear all;
close all;

F_Cin=27; T_C=25; T_H=65; T_D=29;
F_C=21; F_H=17; F_D=13;
Tau_C=100; Tau=40; C=0.15; alpha=6;
H0=72.25; V0=C*H0^2; T0=39.25;

tmax = 1000; 
Tp = 1;
Ts = 0:Tp:tmax;

FH_jumps = [30, 20, 10, -10, -20]; 
FC_jumps = [15, 10, 5, -5, -10]; 

colors = lines(max(length(FH_jumps), length(FC_jumps)));

legend_str_FH = {};
legend_str_FC = {};

figure(1); hold on; grid on;
figure(2); hold on; grid on;

for k = 1:length(FH_jumps)
    u_new = F_H + FH_jumps(k); 
    z = F_C;                   
    
    [~, S_nl] = ode15s(@(t, x) nonlinear_model(x,t,z,T_C,T_H,T_D,F_C,u_new,F_D,Tau_C,C,alpha), Ts, [V0; T0]);
    [~, S_l] = ode15s(@(t, x) linear_model(x,t,z,T_C,T_H,T_D,F_C,u_new,F_D,Tau_C,C,alpha), Ts, [H0; T0]);
    
    h_nl = sqrt(S_nl(:, 1) / C);
    h_l  = S_l(:, 1);
    T_nl = S_nl(:, 2);
    T_l  = S_l(:, 2);
    
    delay_steps = round(Tau / Tp);
    T_out_nl = zeros(size(T_nl));
    T_out_l  = zeros(size(T_l));
    for i = 1:length(Ts)
        idx = max(1, i - delay_steps);
        T_out_nl(i) = T_nl(idx);
        T_out_l(i)  = T_l(idx);
    end
    
    figure(1);
    plot(Ts, h_nl, 'Color', colors(k, :), 'LineWidth', 1.5);
    plot(Ts, h_l, 'Color', colors(k, :), 'LineStyle', '--', 'LineWidth', 1.2);
    
    figure(2);
    plot(Ts, T_out_nl, 'Color', colors(k, :), 'LineWidth', 1.5);
    plot(Ts, T_out_l, 'Color', colors(k, :), 'LineStyle', '--', 'LineWidth', 1.2);
    
    legend_str_FH{end+1} = sprintf('Skok %d, nieliniowy', FH_jumps(k));
    legend_str_FH{end+1} = sprintf('Skok %d, liniowy', FH_jumps(k));
end

figure(3); hold on; grid on;
figure(4); hold on; grid on;

for k = 1:length(FC_jumps)
    u_new = F_H; 
    z = F_C + FC_jumps(k);                   
    
    [~, S_nl] = ode15s(@(t, x) nonlinear_model(x,t,z,T_C,T_H,T_D,F_C,u_new,F_D,Tau_C,C,alpha), Ts, [V0; T0]);
    [~, S_l] = ode15s(@(t, x) linear_model(x,t,z,T_C,T_H,T_D,F_C,u_new,F_D,Tau_C,C,alpha), Ts, [H0; T0]);
    
    h_nl = sqrt(S_nl(:, 1) / C);
    h_l  = S_l(:, 1);
    T_nl = S_nl(:, 2);
    T_l  = S_l(:, 2);
    
    delay_steps = round(Tau / Tp);
    T_out_nl = zeros(size(T_nl));
    T_out_l  = zeros(size(T_l));
    for i = 1:length(Ts)
        idx = max(1, i - delay_steps);
        T_out_nl(i) = T_nl(idx);
        T_out_l(i)  = T_l(idx);
    end
    
    figure(3);
    plot(Ts, h_nl, 'Color', colors(k, :), 'LineWidth', 1.5);
    plot(Ts, h_l, 'Color', colors(k, :), 'LineStyle', '--', 'LineWidth', 1.2);
    
    figure(4);
    plot(Ts, T_out_nl, 'Color', colors(k, :), 'LineWidth', 1.5);
    plot(Ts, T_out_l, 'Color', colors(k, :), 'LineStyle', '--', 'LineWidth', 1.2);
    
    legend_str_FC{end+1} = sprintf('Skok %d, nieliniowy', FC_jumps(k));
    legend_str_FC{end+1} = sprintf('Skok %d, liniowy', FC_jumps(k));
end

figure(1);
set(gcf, 'Color', 'w');
ax1 = gca;
set(ax1, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.8 0.8 0.8]);
xlabel('Czas [s]', 'Interpreter', 'latex');
ylabel('$h$ [cm]', 'Interpreter', 'latex');
lgd1 = legend(legend_str_FH, 'Location', 'eastoutside', 'TextColor', 'k');
title(lgd1, 'Skoki $F_H$', 'Interpreter', 'latex');
set(gcf, 'position', [10, 10, 1000, 700]);
exportgraphics(gcf, "images/multi_jump_FH_h.pdf", "ContentType", "vector");

figure(2);
set(gcf, 'Color', 'w');
ax2 = gca;
set(ax2, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.8 0.8 0.8]);
xlabel('Czas [s]', 'Interpreter', 'latex');
ylabel('$T_{out}$ [$^\circ$C]', 'Interpreter', 'latex');
lgd2 = legend(legend_str_FH, 'Location', 'eastoutside', 'TextColor', 'k');
title(lgd2, 'Skoki $F_H$', 'Interpreter', 'latex');
set(gcf, 'position', [10, 10, 1000, 700]);
exportgraphics(gcf, "images/multi_jump_FH_T.pdf", "ContentType", "vector");

figure(3);
set(gcf, 'Color', 'w');
ax3 = gca;
set(ax3, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.8 0.8 0.8]);
xlabel('Czas [s]', 'Interpreter', 'latex');
ylabel('$h$ [cm]', 'Interpreter', 'latex');
lgd3 = legend(legend_str_FC, 'Location', 'eastoutside', 'TextColor', 'k');
title(lgd3, 'Skoki $F_C$', 'Interpreter', 'latex');
set(gcf, 'position', [10, 10, 1000, 700]);
exportgraphics(gcf, "images/multi_jump_FC_h.pdf", "ContentType", "vector");

figure(4);
set(gcf, 'Color', 'w');
ax4 = gca;
set(ax4, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.8 0.8 0.8]);
xlabel('Czas [s]', 'Interpreter', 'latex');
ylabel('$T_{out}$ [$^\circ$C]', 'Interpreter', 'latex');
lgd4 = legend(legend_str_FC, 'Location', 'eastoutside', 'TextColor', 'k');
title(lgd4, 'Skoki $F_C$', 'Interpreter', 'latex');
set(gcf, 'position', [10, 10, 1000, 700]);
exportgraphics(gcf, "images/multi_jump_FC_T.pdf", "ContentType", "vector");