function [ M ] = DMC_M(S, N, Nu)
    [~, ny, nu] = size(S);
    M = zeros(ny * N, nu * Nu);
    
    for i = 1:N
        for j = 1:Nu
            if i >= j
                Sp = S(i - j + 1, :, :);
            else
                continue;
            end
            M((i - 1) * ny + 1 : i * ny, (j - 1) * nu + 1 : j * nu) = Sp;
        end
    end
end 
