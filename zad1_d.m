
Tp = 1;
metoda = 'zoh';
sys = get_transfer_func();
sysd = c2d(sys, Tp, metoda);
% sysd1 = c2d(sys, 1, metoda);
% sysd5 = c2d(sys, 5, metoda);
% sysd10 = c2d(sys, 10, metoda);
% sysd25 = c2d(sys, 25, metoda);

%% porównanie
% figure;
% 
% hold on;
% step(sys(2,3)); 
% step(sysd1(2,3));
% step(sysd5(2,3));
% step(sysd10(2,3));
% step(sysd25(2,3));
% title("")
% ylabel("T_{out} [°C]")
% xlabel("Czas ") %
% legend("ciągły","Tp=1","Tp=5","Tp=10","Tp=25")

%%
% symulacja w równaniach stanu
F_C=21;
F_H=17;
F_D=13;
Tau_C=100;
Tau=40;
h=72.25;
T=39.35;
kmax=400;

[Ad, Bd, Cd, Dd] = ssdata(sysd);

x0 = [h; T];
u0 = [F_C; F_H; F_D];
dx = zeros(2,kmax);
dy = zeros(2,kmax);

du = [0;0;0];
dx(:,1) = [0;0];
y=zeros(2,kmax);
% y(:,1)=x0;
delay_T = Tau/Tp;
y_buffer = zeros(1, delay_T);
for i=2:kmax

    if i > Tau_C/Tp
        du = [1;0;0]; 
    else
        du = [0;0;0];
    end

    [dy(:,i), dx(:,i)] = discrete_linear_model(dx(:,i-1),...
        du, Ad, Bd, Cd, Dd);
    y(1,i) = dy(1,i); %+ x0(1);

    y_buffer = [dy(2,i), y_buffer(1:end-1)];
    y(2,i) = y_buffer(end); %+ x0(2);
end


figure;
hold on;
stairs(y(1,2:end),'Color','red');
%step(sysd1(1,1));
step(sysd1(1,1));
title('')
legend('model w przestrzeni stanów','transmitancja')
ylabel("h [cm]")
%ylabel("T_{out} [°C]")
xlabel('Czas ');




