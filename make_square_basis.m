function F = make_square_basis(kx_max, ky_max)
    %F = {};
    if mod(kx_max,2) == 0 
        kx_range = -kx_max+1:2:kx_max-1;
    else
        kx_range = -(kx_max-1)/2:(kx_max-1)/2;
    end
    if mod(ky_max,2) == 0
        ky_range = ky_max-1:-2:-ky_max+1;
    else
        ky_range = (ky_max-1)/2:-1:-(ky_max-1)/2;
    end
    
    [X,Y] = meshgrid(kx_range, ky_range);

    %for kx = 1:size(X,2)
    %    for ky = 1:size(Y,1)
    %        F = [F {[X(ky,kx);Y(ky,kx)]}]; %#ok<AGROW> 
    %    end
    %end
    F = cat(3,X,Y);
end