function [h_nowe, T_nowe] = obiekt(h_stare, T_stare, F_H, F_C, F_D)
    Tp = 1;
    C = 0.15;
    alpha = 6;
    Tau_C = 100;
    
    T_C = 25;
    T_H = 65;
    T_D = 29;
    
    % Stan początkowy to bezpośrednio [h; T], tak jak tego oczekuje Twój model
    x0 = [h_stare; T_stare];
    
    % Wywołujemy TWÓJ plik nonlinear_model.
    % Przekazujemy F_C dwa razy (dla F_Cin i F_C), bo opóźnienie 
    % ogarnęliśmy już na zewnątrz w pętli for.
    ode_fun = @(t, x) nonlinear_model(x, t, F_C, T_C, T_H, T_D, F_C, F_H, F_D, Tau_C, C, alpha);
    
    % Solver na dystansie jednej próbki
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    [~, X_out] = ode15s(ode_fun, [0 Tp], x0, options);
    
    % Zwracamy wyliczone wartości
    h_nowe = X_out(end, 1);
    T_nowe = X_out(end, 2);
end