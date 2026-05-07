%% 3. PRZYGOTOWANIE DO EKSPERYMENTÓW
% ZAMROŻONE PARAMETRY DMC (poza N i Nu)
D = 450; 
D_z = 450;           
lambda = [10 10];  % <--- ZWIĘKSZONO! Uspokaja to zawory do testu N
psi = [1 1];         

% BADANE WARTOŚCI N:
N_wartosci = [50, 100, 140, 150, 170, 200];
liczba_testow = length(N_wartosci);

% Macierze do przechowywania wyników
h_wyniki = zeros(liczba_testow, t_k);
T_wyniki = zeros(liczba_testow, t_k);

disp('Rozpoczynam pętlę eksperymentów dla wielu wartości N...');

%% 4. PĘTLA PO RÓŻNYCH N
for iter = 1:liczba_testow
    N = N_wartosci(iter);
    Nu = N; % <--- KLUCZOWE: Uwalniamy Nu! Równamy horyzont sterowania z N.
    
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
end
disp('Eksperymenty zakończone!');

%% 5. Rysowanie dynamicznych wykresów
figure('Name', 'Analiza wpływu N na jakość regulacji', 'Position', [100, 100, 1000, 800]);
kolory = lines(liczba_testow); % Wygenerowanie ładnej palety kolorów

% Wykres h
subplot(2,1,1);
% Zmniejszono grubość czerwonej przerywanej z 2 na 1.5
plot(1:t_k, h_zad, 'k--', 'LineWidth', 1.5); hold on; 
for i = 1:liczba_testow
    % Zmniejszono grubość linii z 1.2/1.5 na 1 (lub wpisz 0.5 dla bardzo cienkich)
    plot(1:t_k, h_wyniki(i, :), 'Color', kolory(i,:), 'LineWidth', 1);
end
title('Wpływ horyzontu N na regulację wysokości (h)');
ylabel('Wysokość [cm]');
legend_text = ['Wartość zadana', arrayfun(@(n) sprintf('N = %d', n), N_wartosci, 'UniformOutput', false)];
legend(legend_text, 'Location', 'SouthEast');
grid on;

% Wykres T
subplot(2,1,2);
% Zmniejszono grubość czerwonej przerywanej
plot(1:t_k, T_zad, 'k--', 'LineWidth', 1.5); hold on;
for i = 1:liczba_testow
    % Zmniejszono grubość linii na 1
    plot(1:t_k, T_wyniki(i, :), 'Color', kolory(i,:), 'LineWidth', 1);
end
title('Wpływ horyzontu N na regulację temperatury (T_{out})');
ylabel('Temperatura [°C]');
xlabel('Czas [s]');
legend(legend_text, 'Location', 'SouthEast');
grid on;