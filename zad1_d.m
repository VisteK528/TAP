
Tp = 1;
metoda = 'zoh';
sys = get_transfer_func();
sysd = c2d(sys, Tp, metoda);

% porównanie
step(sys); hold on;
step(sysd);
legend('ciągły','dyskretny');

%%
% symulacja w równaniach stanu
F_C=21;
F_H=17;
F_D=13;
Tau_C=100;
Tau=40;
h=72.25;
T=39.25;
kmax=400;

[Ad, Bd, Cd, Dd] = ssdata(sysd);

x0 = [h; T];
u0 = [F_C; F_H; F_D];
dx = zeros(2,kmax);
dy = zeros(2,kmax);

du = [0;0;0];
dx(:,1) = [0;0];
y=zeros(2,kmax);
delay_T = Tau/Tp;
y_buffer = zeros(1, delay_T);
for i=2:kmax
    if i > Tau_C/Tp
        du = [1;0;0]; 
    else
        du = [0;0;0];
    end

    [dy(:,i), dx(:,i)] = discrete_linear_model(dx(:,i-1), du, Ad, Bd, Cd, Dd);
    % do absolutnych trzeba dodac stan ustalony i 1szy indeks. 
    % do symulacji sie doda ale do stepa mozna zostawic jak jest bo lepiej widac 
    y(1,i) = dy(1,i);% + x0(1);

    y_buffer = [dy(2,i), y_buffer(1:end-1)];
    y(2,i) = y_buffer(end);% + x0(2);
end


figure;
stairs(y(1,:)); hold on;
stairs(y(2,:));
xlabel('chwila k');




