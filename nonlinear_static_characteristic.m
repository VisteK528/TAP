clear all;
close all;

F_C_op = 21; 
F_H_op = 17;
F_D = 13;
T_C = 25; T_H = 65; T_D = 29;
C = 0.15; alpha = 6; Tau_C = 100;
H0 = 72.25; V0 = C*H0^2; T0 = 39.25;

Tp = 1.0;
tmax = 5000; 
Ts = 0:Tp:tmax;
N = 15;

F_C_range = linspace(F_C_op*0.5, F_C_op*1.5, N);
F_H_range = linspace(F_H_op*0.5, F_H_op*1.5, N);

[FC_mesh, FH_mesh] = meshgrid(F_C_range, F_H_range);
H_steady = zeros(N, N);
T_steady = zeros(N, N);

for i = 1:N
    for j = 1:N
        f_c = FC_mesh(i,j);
        f_h = FH_mesh(i,j);
        
        [t, S] = ode15s(@(t, x) nonlinear_model(x,t,f_c,T_C,T_H,T_D,f_c,f_h,F_D,Tau_C,C,alpha), ...
                 [0, tmax], [V0; T0]);
        
        H_steady(i,j) = sqrt(S(end, 1) / C);
        T_steady(i,j) = S(end, 2);
    end
end

figure(1);
surf(FC_mesh, FH_mesh, H_steady, 'FaceAlpha', 0.8);
hold on;
plot3(F_C_op, F_H_op, H0, 'r.', 'MarkerSize', 30); 

xlabel('$F_C$ [cm$^3$/s]', 'Interpreter', 'latex');
ylabel('$F_H$ [cm$^3$/s]', 'Interpreter', 'latex');
zlabel('$h$ [cm]', 'Interpreter', 'latex');
grid on; colorbar;

x0 = 10; y0 = 10; width = 1000; height = 720;
set(gcf, 'position', [x0, y0, width, height]);
exportgraphics(gcf, "images/static_char_h_nonlinear.pdf", "Resolution", 200);

figure(2);
surf(FC_mesh, FH_mesh, T_steady, 'FaceAlpha', 0.8);
hold on;
plot3(F_C_op, F_H_op, T0, 'r.', 'MarkerSize', 30);

xlabel('$F_C$ [cm$^3$/s]', 'Interpreter', 'latex');
ylabel('$F_H$ [cm$^3$/s]', 'Interpreter', 'latex');
zlabel('$T$ [$^\circ$C]', 'Interpreter', 'latex');
grid on; colorbar;

set(gcf, 'position', [x0, y0, width, height]);
exportgraphics(gcf, "images/static_char_T_nonlinear.pdf", "Resolution", 200);