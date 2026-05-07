function [S, Sz] = get_s(steps, F_C, F_H, F_D, T_C, T_H, T_D, Tau_C, C, alpha, Hpp, Tpp, Tp, delay_u1, delay_steps)
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-7);
    t_span = 0:Tp:(steps*Tp); H0 = Hpp; T0 = Tpp;
    
    dF_C = F_C * 0.1;
    [~, out1] = ode15s(@(t, x) linear_model(x, t, F_C+dF_C, T_C, T_H, T_D, F_C, F_H, F_D, Tau_C, C, alpha), t_span, [H0; T0], options);
    s11_r = (out1(:, 1) - Hpp) / dF_C;
    s21_r = (out1(:, 2) - Tpp) / dF_C;
    
    dF_H = F_H * 0.1;
    [~, out2] = ode15s(@(t, x) linear_model(x, t, F_C, T_C, T_H, T_D, F_C, F_H+dF_H, F_D, Tau_C, C, alpha), t_span, [H0; T0], options);
    s12_r = (out2(:, 1) - Hpp) / dF_H;
    s22_r = (out2(:, 2) - Tpp) / dF_H;
    
    dF_D = F_D * 0.1;
    [~, out3] = ode15s(@(t, x) linear_model(x, t, F_C, T_C, T_H, T_D, F_C, F_H, F_D+dF_D, Tau_C, C, alpha), t_span, [H0; T0], options);
    s13_r = (out3(:, 1) - Hpp) / dF_D;
    s23_r = (out3(:, 2) - Tpp) / dF_D;
    
    shift = @(v, d) [zeros(d, 1); v(1:end-d)];
    
    S_raw = [shift(s11_r, delay_u1), shift(s21_r, delay_u1 + delay_steps), s12_r, shift(s22_r, delay_steps)];
    S = S_raw(2:end, :);
    
    Sz_raw = [s13_r, shift(s23_r, delay_steps)];
    Sz = Sz_raw(2:end, :);
end