close all; clear; clc;

%% 1. Parametry symulacji i punkt pracy
Tp = 1;              
t_p = 200;
t_k = 2500;

h0 = 72.25;          
T0 = 39.35;          
F_H0 = 17;           
F_C0 = 21;           
F_D0 = 13;           

tau_C = 100;         
tau   = 40;          

%% 2. Przechowywanie wyników i wartości zadane
F_Hin = F_H0 * ones(1, t_k);
F_Cin = F_C0 * ones(1, t_k);

F_D = [F_D0 * ones(1, round(t_k/2)), 1.2 * F_D0 * ones(1, t_k - round(t_k/2))];

h = h0 * ones(1, t_k);
T = T0 * ones(1, t_k);

h_zad = h0 * ones(1, t_k); 
h_zad(500:end) = 1.5 * h0;

T_zad = T0 * ones(1, t_k); 
T_zad(1000:end) = 1.05 * T0;

%% 3. Nastawy regulatora DMC
D = 450;      %600 200 10? 10 5    1 70       
D_z = 450;           
N = 170;             
Nu = 3;             
lambda = [10 5];
psi = [1, 10];


%% 4. OBLICZENIA OFFLINE
sys = get_transfer_func(); 
metoda = 'zoh';
sysd = c2d(sys, Tp, metoda);

czas = 0:Tp:(D + N + Nu)*Tp;
S_all = step(sysd, czas);
S_all = S_all(2:end, :, :);


S = S_all(:, :, [2, 1]); 
S_z = S_all(:, :, 3);

[Ke, Ku, Kz] = DMC_anal_offline(S, S_z, D, D_z, N, Nu, lambda, psi);

dUp = zeros(2, D-1);
dZ = zeros(1, D_z);

U = [F_H0; F_C0] * ones(1, round(t_p/Tp)); 
z_zak = F_D0 * ones(1, round(t_p/Tp)); 

%% 5. Ograniczenia sterowania
U_min = [0; 0];                  
U_max = [100; 100];              
dU_min = [-10; -10];             
dU_max = [10; 10];               

%% 6. GŁÓWNA PĘTLA SYMULACYJNA ONLINE
disp('Rozpoczynam symulację...');
for t = t_p:1:t_k
    
    indeks_FCin = max(1, t - tau_C);
    [h(t), T(t)] = obiekt(h(t-1), T(t-1), F_Hin(t), F_Cin(indeks_FCin), F_D(t)); 
    
    indeks_T = max(1, t - tau);
    Tout_measured = T(indeks_T);

    if mod(t, Tp) == 0
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
    end
    
    F_Hin(t+1) = U(1, k); 
    F_Cin(t+1) = U(2, k); 
end
disp('Symulacja zakończona.');

%% 7. Wykresy
figure('Name', 'Wyniki symulacji DMC', 'Position', [100, 100, 1000, 800]);


subplot(4,1,1);
plot(1:t_k, h_zad, 'r--', 1:t_k, h, 'b-', 'LineWidth', 1.5);
title('Regulacja poziomu cieczy (h)');
ylabel('Wysokość [cm]');
legend('Wartość zadana', 'Pomiar');
grid on;

subplot(4,1,2);
plot(1:t_k, T_zad, 'r--', 1:t_k, T, 'b-', 'LineWidth', 1.5);
title('Regulacja temperatury (T_{out})');
ylabel('Temperatura [°C]');
legend('Wartość zadana', 'Pomiar');
grid on;



subplot(4,1,3);

stairs(1:t_k+1, F_Hin, 'r', 'LineWidth', 1.5); hold on; 
stairs(1:t_k+1, F_Cin, 'b', 'LineWidth', 1.5);          
title('Sygnały sterujące');
ylabel('Przepływ [cm^3/s]');
legend('F_H (Woda gorąca)', 'F_C (Woda zimna)');
grid on;

subplot(4,1,4);
stairs(1:t_k, F_D, 'k', 'LineWidth', 1.5);
title('Zakłócenie mierzone');
xlabel('Czas [s]');
ylabel('Przepływ F_D [cm^3/s]');
grid on;