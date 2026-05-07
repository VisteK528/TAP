function [ Mp ] = DMC_Mp(S, D, N)
    [~, ny, nu] = size(S);
    Mp = zeros(ny * N, nu * (D - 1));

    for i = 1:N
        for j = 1:(D - 1)
            Sp_diff = S(i + j, :, :) - S(j, :, :);
            Mp((i - 1) * ny + 1 : i * ny, (j - 1) * nu + 1 : j * nu) = Sp_diff;
        end
    end
end