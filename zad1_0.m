
F_Cin=50;
T_C=25;
T_H=65;
T_D=29;
F_C=21;
F_H=17;
F_D=13;
Tau_C=100;
Tau=40;
C=0.15;
alpha=6;
h=72.25;
V0=C*h^2;
T0=39.25;
tmax=1000;
Ts=0:0.1:tmax;
[t, S]=ode45(@(t, x) nonlinear_model(x,t,F_Cin,T_C,T_H,T_D,F_C,F_H,F_D,Tau_C,C,alpha), Ts, [V0; T0]);
h=sqrt(S(:,1)/C);
T=S(:, 2);
T_out = zeros(1, length(S(:, 2)));

for i=1:length(T_out)
    T_out(i) = T(max(1, floor(i-Tau/0.1)));
end

plot(h);
%plot(T);
hold on;
plot(T_out);
grid on;
legend({"h", "T_out"});
