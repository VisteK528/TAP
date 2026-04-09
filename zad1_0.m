close all;
clear all;
model_config;

tmax=1000;
dt = 0.1;
Ts=0:dt:tmax;

options = odeset('RelTol', 1e-9, 'AbsTol', 1e-7);
[t, S]=ode15s(@(t, x) nonlinear_model(x,t,F_Cin,T_C,T_H,T_D,F_C,F_H,F_D,Tau_C,C,alpha), Ts, [h0; T0], options);

h = S(:,1);
T = S(:, 2);

T_out = zeros(length(Ts), 1);
delay_steps = round(Tau/dt);
for i=1:length(Ts)
    T_out(i) = T(max(1, i - delay_steps));
end

figure;

subplot(2, 1, 1);
plot(Ts, h, 'LineWidth', 1.5);
xlabel('Czas [s]');
ylabel('Wysokość [cm]');
grid on;

subplot(2, 1, 2);
plot(Ts, T_out, 'LineWidth', 1.5);
xlabel('Czas [s]');
ylabel('Temperatura [°C]');
legend({'T_{out}'});
grid on;

x0 = 10;
y0 = 10;
width = 1280;
height = 720;
set(gcf, 'position', [x0, y0, width, height]);
name = "images/operating_point_nonlinear_model.pdf";
exportgraphics(gcf, name, "Resolution", 200);