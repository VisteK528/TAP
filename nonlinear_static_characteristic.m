clear all;
close all;

model_config;
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-7);

Tp = 1.0;
tmax = 5000; 
Ts = 0:Tp:tmax;
N = 15;

F_C_range = linspace(F_C*0.5, F_C*1.5, N);
F_H_range = linspace(F_H*0.5, F_H*1.5, N);

[FC_mesh, FH_mesh] = meshgrid(F_C_range, F_H_range);
H_steady = zeros(N, N);
T_steady = zeros(N, N);

for i = 1:N
    for j = 1:N
        f_c = FC_mesh(i,j);
        f_h = FH_mesh(i,j);
        
        [t, S] = ode15s(@(t, x) nonlinear_model(x,t,f_c,T_C,T_H,T_D,f_c,f_h,F_D,Tau_C,C,alpha), ...
                 [0, tmax], [h0; T0], options);
        
        H_steady(i,j) = S(end, 1);
        T_steady(i,j) = S(end, 2);
    end
end

figure(1);
surf(FC_mesh, FH_mesh, H_steady, 'FaceAlpha', 0.8);
hold on;
plot3(F_C, F_H, h0, 'r.', 'MarkerSize', 30); 

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
plot3(F_C, F_H, T0, 'r.', 'MarkerSize', 30);

xlabel('$F_C$ [cm$^3$/s]', 'Interpreter', 'latex');
ylabel('$F_H$ [cm$^3$/s]', 'Interpreter', 'latex');
zlabel('$T$ [$^\circ$C]', 'Interpreter', 'latex');
grid on; colorbar;

set(gcf, 'position', [x0, y0, width, height]);
exportgraphics(gcf, "images/static_char_T_nonlinear.pdf", "Resolution", 200);