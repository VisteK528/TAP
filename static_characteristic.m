clear all;
% --- Parametry ---
F_Cin=27; T_C=25; T_H=65; T_D=29;
F_C=21; F_H=17; F_D=13;
Tau_C=100; Tau=40; C=0.15; alpha=6;
H0=72.25; V0=C*H0^2; T0=39.25;

Ts = 1:1:5000; % Zwiększyłem krok dla szybkości, dostosuj wg potrzeb
Tp = 0.1;
N = 10; % Większe N da ładniejszą powierzchnię 3D

U_jump_points = linspace(F_H*0.5, F_H*1.5, N);
Z_jump_points = linspace(F_C*0.5, F_C*1.5, N);

y_steady = zeros(N, N);
y_steady_nonlinear = zeros(N, N);
y_steady_2 = zeros(N, N);
y_steady_2_nonlinear = zeros(N, N);

i = 1;
for u = U_jump_points
    j = 1;
    for z = Z_jump_points
        % Symulacja nieliniowa
        [t_nl, S_nl] = ode45(@(t, x) nonlinear_model(x,t,z,T_C,T_H,T_D,F_C,u,F_D,Tau_C,C,alpha), Ts, [V0; T0]);
        % Symulacja liniowa
        [t_l, S_l] = ode45(@(t, x) linear_model(x,t,z,T_C,T_H,T_D,F_C,u,F_D,Tau_C,C,alpha), Ts, [H0; T0]);
        
        % Wyniki h (dla nieliniowego przeliczamy z V na h)
        h_lin = S_l(:, 1);
        h_nonlin = sqrt(S_nl(:, 1) / C);
        
        % Wyniki T (zakładam, że Tout to stan T z opóźnieniem Tau)
        % W stanie ustalonym Tout = T, więc dla wykresu 3D bierzemy ostatnią wartość
        y_steady(i, j) = h_lin(end);
        y_steady_nonlinear(i, j) = h_nonlin(end);
        
        y_steady_2(i, j) = S_l(end, 2);
        y_steady_2_nonlinear(i, j) = S_nl(end, 2);
        
        j = j + 1;
    end
    fprintf("Postęp: %d%%\n", round(i/N*100));
    i = i + 1;
end

% Przygotowanie siatki do wykresu (transpozycja wymiarów dla zgodności z i,j)
[U1, U2] = meshgrid(Z_jump_points, U_jump_points);

% --- OKNO 1: Porównanie wysokości h ---
figure(1);
% Model nieliniowy jako kolorowa powierzchnia
surf(U1, U2, y_steady_nonlinear, 'FaceAlpha', 0.6, 'EdgeColor', 'none'); 
hold on;
% Model liniowy jako czarna siatka (mesh) dla kontrastu
mesh(U1, U2, y_steady, 'EdgeColor', 'k'); 
xlabel('$F_C$ (zakłócenie)', 'Interpreter', 'latex');
ylabel('$F_H$ (sterowanie)', 'Interpreter', 'latex');
zlabel('$h$', 'Interpreter', 'latex');
title('Charakterystyka statyczna h: Liniowa vs Nieliniowa');
legend('Nieliniowy (Sygnał rzeczywisty)', 'Liniowy (Aproksymacja)');
grid on;

% --- OKNO 2: Porównanie temperatury T ---
figure(2);
surf(U1, U2, y_steady_2_nonlinear, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
hold on;
mesh(U1, U2, y_steady_2, 'EdgeColor', 'k');
xlabel('$F_C$ (zakłócenie)', 'Interpreter', 'latex');
ylabel('$F_H$ (sterowanie)', 'Interpreter', 'latex');
zlabel('$T$', 'Interpreter', 'latex');
title('Charakterystyka statyczna T: Liniowa vs Nieliniowa');
legend('Nieliniowy (Sygnał rzeczywisty)', 'Liniowy (Aproksymacja)');
grid on;