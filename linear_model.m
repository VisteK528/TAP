function [x]=linear_model(x,t,F_Cin,T_C,T_H,T_D,F_C,F_H,F_D, Tau_C,C,alpha)
    T_C0=25;
    T_H0=65;
    T_D0=29;
    F_C0=21;
    F_H0=17;
    F_D0=13;
    h0=72.25;
    T0=39.35;

    h=x(1);
    T=x(2);
    if t-Tau_C>1
        F_C=F_Cin;
    end

    x(1) = (-alpha*sqrt(h0))/(4*C*h0) - (alpha * h) / (4*C*h0*sqrt(h0)) + (F_H + F_C + F_D)/(2*C*h0);
    x(2) = (-(F_H0 + F_C0 + F_D0)*T + (T_H - T0) * F_H + (T_C - T0)*F_C + (T_D-T0)*F_D + T0*(F_H0 + F_C0 + F_D0)) / (C*h0^2);
    
end