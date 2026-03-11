function [G]=get_transfer_func
    T_C=25;
    T_H=65;
    T_D=29;
    F_C0=21;
    F_H0=17;
    F_D0=13;
    h0=72.25;
    T_0=39.25;
    C = 0.15;
    alpha = 6;

    A = [-alpha/(4*C*h0*sqrt(h0)) 0; 0 -(F_H0+F_C0+F_D0)/(C*h0^2)];
    B = [1/(2*C*h0) 1/(2*C*h0) 1/(2*C*h0); (T_C-T_0)/(C*h0^2) (T_H-T_0)/(C*h0^2) (T_D-T_0)/(C*h0^2)];
    C = [1 0; 0 1];
    D = [0 0 0; 0 0 0];
    
    sys = ss(A, B, C, D);

    G = tf(sys);
    G.InputDelay = [100, 0, 0];
    G.OutputDelay = [0, 40];
    ends