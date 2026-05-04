function [r2, r1, r0] = discrete_pid_params(Kp, Ti, Td, Tp)
r2 = Kp*Td/Tp;
r1 = Kp*(Tp/(2*Ti) - 2*(Td/Tp) - 1);
r0 = Kp*(1+(Tp/(2*Ti)+Td/Tp));
end