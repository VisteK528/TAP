
function [ Mpz ] = DMC_Mpz(S_z, D_z, N)
    [~, ny, ~] = size(S_z);
    S_temp = S_z';
    S_temp = S_temp(:);
    
    Mpz = [ S_temp(1 : N * ny, :) DMC_Mp(S_z, D_z, N) ];
end