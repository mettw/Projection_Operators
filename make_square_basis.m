function F = make_square_basis(kx_max, ky_max)
    %F = {};
    if mod(kx_max,2) == 0 && mod(ky_max,2) == 0
        [X,Y] = meshgrid(-kx_max+1:2:kx_max-1, ...
            ky_max-1:-2:-ky_max+1);
    else
        [X,Y] = meshgrid(-(kx_max-1)/2:(kx_max-1)/2, ...
            (ky_max-1)/2:-1:-(ky_max-1)/2);
    end

    %for kx = 1:size(X,2)
    %    for ky = 1:size(Y,1)
    %        F = [F {[X(ky,kx);Y(ky,kx)]}]; %#ok<AGROW> 
    %    end
    %end
    F = cat(3,X,Y);
end