function F = make_square_basis3D(kx_max, ky_max, kz_max)
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
    if mod(kz_max,2) == 0 
        kz_range = kz_max-1:-2:-kz_max+1;
    else
        kz_range = (kz_max-1)/2:-1:-(kz_max-1)/2;
    end
    
    [X,Y,Z] = meshgrid(kx_range, ky_range, kz_range);

    %for kx = 1:size(X,2)
    %    for ky = 1:size(Y,1)
    %        F = [F {[X(ky,kx);Y(ky,kx)]}]; %#ok<AGROW> 
    %    end
    %end

    %F = cat(4, cat(3,X(:,:,1),Y(:,:,1),Z(:,:,1)), ...
    %            cat(3,X(:,:,2),Y(:,:,2),Z(:,:,2)));
    F = cat(4, X,Y,Z);
    %{
    F = [];
    for i=1:size(X,3)
        F = cat(4, F, cat(3,X(:,:,i),Y(:,:,i),Z(:,:,i)));
    end
    %}
end