function [M, MP, MPz] = DMC_offline(S, Sz, D, N, Nu)
    ny = 2; nu = 2; nz = 1; S_cell = cell(1, D);
    for i = 1:D, S_cell{i} = [S(i, 1), S(i, 3); S(i, 2), S(i, 4)]; end
    
    M = zeros(N*ny, Nu*nu);
    for i = 1:N
        for j = 1:Nu
            if i >= j, M((i-1)*ny+1:i*ny, (j-1)*nu+1:j*nu) = S_cell{i-j+1}; end
        end
    end
    
    MP = zeros(N*ny, (D-1)*nu);
    for i = 1:N
        for j = 1:D-1
            if i+j <= D, MP((i-1)*ny+1:i*ny, (j-1)*nu+1:j*nu) = S_cell{i+j} - S_cell{j};
            else, MP((i-1)*ny+1:i*ny, (j-1)*nu+1:j*nu) = S_cell{D} - S_cell{j}; end
        end
    end
    
    MPz = zeros(N*ny, (D-1)*nz);
    for i = 1:N
        for j = 1:D-1
            if i+j <= D
                MPz((i-1)*ny+1 : i*ny, j) = Sz(i+j, :)' - Sz(j, :)';
            else
                MPz((i-1)*ny+1 : i*ny, j) = Sz(D, :)' - Sz(j, :)';
            end
        end
    end
end