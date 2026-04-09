function [y,xk1] = discrete_linear_model(x,u,Ad,Bd,Cd,Dd)
xk1=Ad*x+Bd*u;
y=Cd*x+Dd*u;
end

