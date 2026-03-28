clear all;
close all;

F_Cin=27; T_C=25; T_H=65; T_D=29;
F_C=21; F_H=17; F_D=13;
Tau_C=100; Tau=40; C=0.15; alpha=6;
H0=72.25; V0=C*H0^2; T0=39.25;

Ts = 1:1:5000; 
N = 10; 
U_jump_points = linspace(F_H*0.5, F_H*1.5, N);
Z_jump_points = linspace(F_C*0.5, F_C*1.5, N);

y_steady = zeros(N, N);
y_steady_nonlinear = zeros(N, N);
y_steady_2 = zeros(N, N);
y_steady_2_nonlinear = zeros(N, N);

for i = 1:N
    u = U_jump_points(i);
    for j = 1:N
        z = Z_jump_points(j);
        
        [~, S_nl] = ode15s(@(t, x) nonlinear_model(x,t,z,T_C,T_H,T_D,F_C,u,F_D,Tau_C,C,alpha), Ts, [V0; T0]);
        [~, S_l] = ode15s(@(t, x) linear_model(x,t,z,T_C,T_H,T_D,F_C,u,F_D,Tau_C,C,alpha), Ts, [H0; T0]);
        
        y_steady(i, j) = S_l(end, 1);
        y_steady_nonlinear(i, j) = sqrt(S_nl(end, 1) / C);
        
        y_steady_2(i, j) = S_l(end, 2);
        y_steady_2_nonlinear(i, j) = S_nl(end, 2);
    end
end

[U1, U2] = meshgrid(Z_jump_points, U_jump_points);

figure(1);
surf(U1, U2, y_steady_nonlinear, 'FaceAlpha', 0.6, 'EdgeColor', 'none'); 
hold on;
surf(U1, U2, y_steady, 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'EdgeAlpha', 0.2); 
plot3(F_C, F_H, H0, 'r.', 'MarkerSize', 35);
xlabel('$F_C$ [cm$^3$/s]', 'Interpreter', 'latex');
ylabel('$F_H$ [cm$^3$/s]', 'Interpreter', 'latex');
zlabel('$h$ [cm]', 'Interpreter', 'latex');
legend({'Model nieliniowy', 'Model liniowy', 'Punkt pracy'}, 'Location', 'best');
grid on;

x0 = 10; y0 = 10; width = 1000; height = 720;
set(gcf, 'position', [x0, y0, width, height]);
exportgraphics(gcf, "images/comparison_h_linear_nonlinear.pdf", "ContentType", "vector");

figure(2);
surf(U1, U2, y_steady_2_nonlinear, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
hold on;
surf(U1, U2, y_steady_2, 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r', 'EdgeAlpha', 0.2);
plot3(F_C, F_H, T0, 'r.', 'MarkerSize', 35);
xlabel('$F_C$ [cm$^3$/s]', 'Interpreter', 'latex');
ylabel('$F_H$ [cm$^3$/s]', 'Interpreter', 'latex');
zlabel('$T$ [$^\circ$C]', 'Interpreter', 'latex');
legend({'Model nieliniowy', 'Model liniowy', 'Punkt pracy'}, 'Location', 'best');
grid on;

set(gcf, 'position', [x0, y0, width, height]);
exportgraphics(gcf, "images/comparison_T_linear_nonlinear.pdf", "ContentType", "vector");