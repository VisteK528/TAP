function [du1, du2] = numDMC_online(H, M, Psi, y_free, y_zad_curr, u_prev, u_min, u_max, y_min, y_max, du_min, du_max, N, Nu)
    nu = 2; 
    Y_zad = repmat(y_zad_curr, N, 1);
    
    f = -2 * M' * Psi * (Y_zad - y_free);
    
    J = kron(tril(ones(Nu)), eye(nu));
    Umin = repmat(u_min, Nu, 1);
    Umax = repmat(u_max, Nu, 1);
    Uprev = repmat(u_prev, Nu, 1);    
    Ymin = repmat(y_min, N, 1);
    Ymax = repmat(y_max, N, 1);
    
    A = [-J; J; -M; M];
    b = [-(Umin - Uprev); (Umax - Uprev); -(Ymin - y_free); (Ymax - y_free)];
    
    lb = repmat(du_min, Nu, 1);
    ub = repmat(du_max, Nu, 1);
 
    options = optimoptions('quadprog', 'Display', 'off');
    du_full = quadprog(H, f, A, b, [], [], lb, ub, [], options);
    
    if isempty(du_full)
        du1 = 0; du2 = 0; 
    else
        du1 = du_full(1);
        du2 = du_full(2);
    end
end