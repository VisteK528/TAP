clear; clc;

% 1. Parametry z pierwszego etapu
Tp = 1;     % Czas próbkowania (1 sekunda)
D = 450;    % Horyzont dynamiki (tyle próbek potrzebuje obiekt na ustabilizowanie)

% Deklaracja zmiennej 'z' żeby MATLAB wiedział, że robimy model dyskretny
z = tf('z', Tp);

% --- 2. ZBUDOWANIE MODELU Z LICZB Z WASZEGO SPRAWOZDANIA ---

% Modele od sterowań (F_H i F_C) na wyjścia (h i T_out)
G_FH_h = 0.04576 / (z - 0.9838);
G_FC_h = z^(-100) * 0.04576 / (z - 0.9838); % tu jest opóźnienie 100
G_FH_T = z^(-40) * 0.03171 / (z - 0.9369);  % tu jest opóźnienie 40
G_FC_T = z^(-140) * (-0.01774) / (z - 0.9369); % tu opóźnienie 140

% Złożenie ich w jeden system MIMO (2 wyjścia, 2 wejścia)
Model_U = [G_FH_h, G_FC_h;
           G_FH_T, G_FC_T];

% Modele od zakłócenia (F_D) na wyjścia (h i T_out)
G_FD_h = 0.04576 / (z - 0.9838);
G_FD_T = z^(-40) * (-0.0128) / (z - 0.9369);

% Złożenie ich w system zakłóceń (2 wyjścia, 1 wejście)
Model_Z = [G_FD_h;
           G_FD_T];

% --- 3. WŁAŚCIWE TWORZENIE MACIERZY S ---

% Funkcja step robi tu całą magię. Robi skoki jednostkowe na modelu
% i zapisuje je do zmiennych S i S_z przez D próbek.
[S, ~] = step(Model_U, D-1);
[S_z, ~] = step(Model_Z, D-1);

% Wyświetlmy wymiary, żeby potwierdzić, że się udało:
disp('Wymiary wygenerowanej macierzy S:');
disp(size(S));

disp('Wymiary wygenerowanej macierzy S_z:');
disp(size(S_z));