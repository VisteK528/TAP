F_Cin=21; T_C=25; T_H=65; T_D=29;
F_C=21; F_H=17; F_D=13;
Tau_C=100; Tau=40; C=0.15; alpha=6;
h0=72.25; V0=C*h0^2; T0=39.25;

tmax=1000;
dt = 0.1;
Ts=0:dt:tmax;

[t, S]=ode15s(@(t, x) nonlinear_model(x,t,F_Cin,T_C,T_H,T_D,F_C,F_H,F_D,Tau_C,C,alpha), Ts, [V0; T0]);

h = sqrt(S(:,1)/C);
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
name = "images/operating_point_nonlinear_model.png";
exportgraphics(gcf, name, "Resolution", 200);