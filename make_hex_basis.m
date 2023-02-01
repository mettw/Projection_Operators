function [F, Kx, Ky, B] = make_hex_basis(m_max, varargin)
    % Flat top or pointed top?
    if nargin == 2 && varargin{1} == true
        b1 = [-1/2, sqrt(3)/2];
        b2 = [1, 0];
    else
        b1 = [sqrt(3)/2, 1/2];
        b2 = [sqrt(3)/2, -1/2];
    end

    [M,N] = meshgrid(-m_max:m_max, m_max:-1:-m_max);

    Kx = M(:)*b1(1)+N(:)*b2(1);
    Ky = M(:)*b1(2)+N(:)*b2(2);
    F = [Kx(sqrt(Kx.^2+Ky.^2)<=m_max),Ky(sqrt(Kx.^2+Ky.^2)<=m_max)];

    B = [b1.' b2.'];
end