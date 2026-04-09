function [x]=nonlinear_model(x,t,F_Cin,T_C,T_H,T_D,F_C,F_H,F_D,Tau_C,C,alpha)
    h=x(1);
    T=x(2);
    if t-Tau_C>1
        F_C=F_Cin;
    end

    x(1)=(F_H+F_C+F_D-alpha*sqrt(h))/(2*C*h);
    x(2)=(F_H*T_H+F_C*T_C+F_D*T_D-(F_H+F_C+F_D)*T)/(C*h^2);
    
end