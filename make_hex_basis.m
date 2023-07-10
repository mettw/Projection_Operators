function [F, Kx, Ky, B] = make_hex_basis(m_max, varargin)
   
    b1 = [sqrt(3)/2, 1/2];
    b2 = [sqrt(3)/2, -1/2];
    [M,N] = meshgrid(-m_max:m_max, m_max:-1:-m_max);

    Kx = (M(:)*b1(1)+N(:)*b2(1))/m_max;
    Ky = (M(:)*b1(2)+N(:)*b2(2))/m_max;
    %F = [Kx(sqrt(Kx.^2+Ky.^2)<=m_max),Ky(sqrt(Kx.^2+Ky.^2)<=m_max)];


    F = [Kx(-b1(1)-1/m_max^2<=Kx & Kx<=b1(1)+1/m_max^2 & ...
        -(b1(2)-b2(2))<=Ky & Ky<=(b1(2)-b2(2))), ...
        Ky(-b1(1)-1/m_max^2<=Kx & Kx<=b1(1)+1/m_max^2 & ...
        -(b1(2)-b2(2))<=Ky & Ky<=(b1(2)-b2(2)))];

    % Flat top or pointed top?
    if nargin == 2 && varargin{1} == true
        b1 = [-1/2, sqrt(3)/2];
        b2 = [1, 0];
        F = F*[cosd(30) -sind(30); sind(30) cosd(30)];
        K = [Kx Ky]*[cosd(30) -sind(30); sind(30) cosd(30)];
        Kx = K(:,1);
        Ky = K(:,2);
    end

    B = [b1; b2];
end