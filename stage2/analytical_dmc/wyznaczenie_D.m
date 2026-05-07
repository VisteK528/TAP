close all; clear; clc;

% 1. Pobranie Twojego modelu ciągłego
Tp = 1;
sys = get_transfer_func(); 
metoda = 'zoh';

% 2. Dyskretyzacja
sysd = c2d(sys, Tp, metoda);

% 3. Symulacja skoku w otwartej pętli (bez regulatora) na długim czasie
czas_symulacji = 1000; % Zakładamy 1000 sekund, żeby na pewno zobaczyć wypłaszczenie
czas = 0:Tp:czas_symulacji;

% step generuje od razu odpowiedzi ze wszystkich wejść na wszystkie wyjścia
S_all = step(sysd, czas);

% Przypominam układ wejść i wyjść w sysd:
% Wejścia: (1) F_C, (2) F_H, (3) F_D
% Wyjścia: (1) h,   (2) T_out

%% 4. Rysowanie wykresów do odczytu D (dla sterowań F_H i F_C)
figure('Name', 'Wyznaczanie D (Sterowania)', 'Position', [100, 100, 900, 600]);

% Skok F_H na h i T
subplot(2,2,1);
plot(czas, S_all(:, 1, 2), 'r', 'LineWidth', 1.5);
title('Skok F_H \rightarrow Wysokość (h)');
grid on;

subplot(2,2,2);
plot(czas, S_all(:, 2, 2), 'r', 'LineWidth', 1.5);
title('Skok F_H \rightarrow Temperatura (T_{out})');
grid on;

% Skok F_C na h i T
subplot(2,2,3);
plot(czas, S_all(:, 1, 1), 'b', 'LineWidth', 1.5);
title('Skok F_C \rightarrow Wysokość (h)');
xlabel('Czas [s] (Liczba próbek D)');
grid on;

subplot(2,2,4);
plot(czas, S_all(:, 2, 1), 'b', 'LineWidth', 1.5);
title('Skok F_C \rightarrow Temperatura (T_{out})');
xlabel('Czas [s] (Liczba próbek D)');
grid on;

%% 5. Rysowanie wykresów do odczytu D_z (dla zakłócenia F_D)
figure('Name', 'Wyznaczanie D_z (Zakłócenie)', 'Position', [150, 150, 900, 300]);

subplot(1,2,1);
plot(czas, S_all(:, 1, 3), 'k', 'LineWidth', 1.5);
title('Skok zakłócenia F_D \rightarrow Wysokość (h)');
xlabel('Czas [s] (Liczba próbek D_z)');
grid on;

subplot(1,2,2);
plot(czas, S_all(:, 2, 3), 'k', 'LineWidth', 1.5);
title('Skok zakłócenia F_D \rightarrow Temperatura (T_{out})');
xlabel('Czas [s] (Liczba próbek D_z)');
grid on;