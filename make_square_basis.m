function F = make_square_basis(kx_max, ky_max)
    %F = {};
    if mod(2*kx_max,2) == 0 && mod(2*ky_max,2) == 0
        [X,Y] = meshgrid(-2*kx_max+1:2:2*kx_max-1, ...
            2*ky_max-1:-2:-2*ky_max+1);
    else
        [X,Y] = meshgrid(-kx_max:kx_max, ky_max:-1:-ky_max);
    end

    %for kx = 1:size(X,2)
    %    for ky = 1:size(Y,1)
    %        F = [F {[X(ky,kx);Y(ky,kx)]}]; %#ok<AGROW> 
    %    end
    %end
    F = cat(3,X,Y);
end