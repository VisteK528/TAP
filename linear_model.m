function [x]=linear_model(x,t,F_Cin,T_C,T_H,T_D,F_C,F_H,F_D, Tau_C,C,alpha)
    T_C0=25;
    T_H0=65;
    T_D0=29;
    F_C0=21;
    F_H0=17;
    F_D0=13;
    %Tau_C=100;
    %Tau=40;
    %C=0.15;
    %alpha=6;
    h0=72.25;
    T0=39.25;

    h=x(1);
    T=x(2);
    if t-Tau_C>1
        F_C=F_Cin;
    end

    x(1)=(F_H0 + F_C0 + F_D0 - alpha * sqrt(h0))/(2*C*h0) + (1/(2*C*h0))*(F_H - F_H0 + F_C - F_C0 + F_D - F_D0) ...
            - (alpha/(4*C*h0*sqrt(h0)))*(h - h0);
    x(2)=(F_H0*T_H+F_C0*T_C+F_D0*T_D-(F_H0+F_C0+F_D0)*T0)/(C*h0^2)+(T_H0-T0)/(C*h0^2)*(F_H-F_H0)+(T_C0-T0)/(C*h0^2)*(F_C-F_C0)+(T_D0-T0)/(C*h0^2)*(F_D-F_D0)+(-F_H0-F_C0-F_D0)/(C*h0^2)*(T-T0);
    
end