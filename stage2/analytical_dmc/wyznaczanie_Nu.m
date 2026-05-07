close all; clear; clc;

%% 1. Parametry wspólne
Tp = 1;              
t_p = 200;           
t_k = 2000;          

h0 = 72.25;          T0 = 39.35;          
F_H0 = 17;           F_C0 = 21;           F_D0 = 13;           
tau_C = 100;         tau   = 40;          

U_min = [0; 0];      U_max = [100; 100];              
dU_min = [-10; -10]; dU_max = [10; 10];               

h_zad = h0 * ones(1, t_k); h_zad(500:end) = 1.1 * h0; 
T_zad = T0 * ones(1, t_k); T_zad(1000:end) = 1.05 * T0; 
F_D = F_D0 * ones(1, t_k); 

% ZAMROŻONE PARAMETRY DMC (Testujemy Nu)
D = 450; 
D_z = 450;
N = 170;             % <--- TWÓJ FAWORYT Z POPRZEDNIEGO BADANIA           
lambda = [1 1];  
psi = [1 1];         

%% 2. Modele Offline
sys = get_transfer_func(); 
sysd = c2d(sys, Tp, 'zoh');

%% 3. PRZYGOTOWANIE DO EKSPERYMENTÓW
% BARDZIEJ SENSOWNE, ZAGĘSZCZONE WARTOŚCI Nu:
Nu_wartosci = [1, 10, 20, 30, 50, 80, 100, 150];
liczba_testow = length(Nu_wartosci);

h_wyniki = zeros(liczba_testow, t_k);
T_wyniki = zeros(liczba_testow, t_k);
FH_wyniki = zeros(liczba_testow, t_k);
FC_wyniki = zeros(liczba_testow, t_k);

disp('Rozpoczynam pętlę eksperymentów dla zagęszczonych wartości Nu...');

%% 4. PĘTLA PO RÓŻNYCH Nu
for iter = 1:liczba_testow
    Nu = Nu_wartosci(iter);
    fprintf('Liczenie symulacji nr %d/%d (N = %d, Nu = %d)...\n', iter, liczba_testow, N, Nu);
    
    czas_S = 0:Tp:(D + N + Nu)*Tp; 
    S_all = step(sysd, czas_S);
    S_all = S_all(2:end, :, :);
    S = S_all(:, :, [2, 1]); 
    S_z = S_all(:, :, 3);
    
    [Ke, Ku, Kz] = DMC_anal_offline(S, S_z, D, D_z, N, Nu, lambda, psi);
    dUp = zeros(2, D-1);
    dZ = zeros(1, D_z);
    
    F_Hin = F_H0 * ones(1, t_k); F_Cin = F_C0 * ones(1, t_k);
    h = h0 * ones(1, t_k);       T = T0 * ones(1, t_k);
    U = [F_H0; F_C0] * ones(1, round(t_p/Tp)); 
    z_zak = F_D0 * ones(1, round(t_p/Tp)); 
    
    for t = t_p:1:t_k
        indeks_FCin = max(1, t - tau_C);
        [h(t), T(t)] = obiekt(h(t-1), T(t-1), F_Hin(t), F_Cin(indeks_FCin), F_D(t)); 
        indeks_T = max(1, t - tau);
        Tout_measured = T(indeks_T);
        
        k = round(t / Tp);
        y = [h(t); Tout_measured]; 
        y_zad = [h_zad(t); T_zad(t)];
        
        z_zak(k) = F_D(t);
        dZ(2:end) = dZ(1:end-1);
        dZ(1) = z_zak(k) - z_zak(k-1);
        
        dU = DMC_anal_dU(y, y_zad, Ke, Ku, Kz, dUp, dZ);
        dU = min(max(dU, dU_min), dU_max);
        
        U_temp = U(:, k-1) + dU;
        U_temp = min(max(U_temp, U_min), U_max);
        dU = U_temp - U(:, k-1);
        
        U(:, k) = U(:, k-1) + dU;
        dUp(:, 2:end) = dUp(:, 1:end-1);
        dUp(:, 1) = dU;
        
        F_Hin(t+1) = U(1, k); 
        F_Cin(t+1) = U(2, k); 
    end
    
    h_wyniki(iter, :) = h;
    T_wyniki(iter, :) = T;
    FH_wyniki(iter, :) = F_Hin(1:t_k);
    FC_wyniki(iter, :) = F_Cin(1:t_k);
end
disp('Eksperymenty zakończone!');

%% 5. Rysowanie dynamicznych wykresów
figure('Name', 'Analiza wpływu Nu na jakość regulacji', 'Position', [100, 100, 1000, 900]);
kolory = lines(liczba_testow); 

% Wykres F_H (Sterowanie gorącą wodą)
subplot(3,1,1);
for i = 1:liczba_testow
    stairs(1:t_k, FH_wyniki(i, :), 'Color', kolory(i,:), 'LineWidth', 1); hold on;
end
title('Wpływ horyzontu Nu na sterowanie (F_H)');
ylabel('Przepływ [cm^3/s]');
legend_text = arrayfun(@(n) sprintf('Nu = %d', n), Nu_wartosci, 'UniformOutput', false);
legend(legend_text, 'Location', 'EastOutside');
grid on;

% Wykres h
subplot(3,1,2);
plot(1:t_k, h_zad, 'k--', 'LineWidth', 1.5); hold on;
for i = 1:liczba_testow
    plot(1:t_k, h_wyniki(i, :), 'Color', kolory(i,:), 'LineWidth', 1);
end
title('Wpływ horyzontu Nu na regulację wysokości (h)');
ylabel('Wysokość [cm]');
grid on;

% Wykres T
subplot(3,1,3);
plot(1:t_k, T_zad, 'k--', 'LineWidth', 1.5); hold on;
for i = 1:liczba_testow
    plot(1:t_k, T_wyniki(i, :), 'Color', kolory(i,:), 'LineWidth', 1);
end
title('Wpływ horyzontu Nu na regulację temperatury (T_{out})');
ylabel('Temperatura [°C]');
xlabel('Czas [s]');
grid on;