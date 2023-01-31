classdef Harmonic
    %HARMONIC - vector representation of a planar harmonic
    %  
    % This object creates a vector describing the field in a planar
    % waveguide in terms of the Fourier coefficients of the field.
    %
    % USAGE:
    %
    % Create a Projector object for convenience
    %
    % > P=Projector("F1", "C4v")
    %
    % Create the harmonic using the projector to give us the basis for the
    % E irreducible representation of C4v.  The first two parameters are
    % the lattice vectors and the last one is the coefficients.  The
    % exp(1i*90) ensures that the phase is set so that the field is not
    % zero and the cosine and sines define the orientation of the mode
    % within the 2D eigenspace of the E irreducible representation.
    %
    % > hm=Harmonic([1 0],[0 1], [1 1], P.C4v.E_vec, [cosd(135)*exp(1i*90), sind(135)*exp(1i*90)])
    %
    % TODO:
    %
    % At the moment it only works with square lattices.
    
    properties
        % Hilbert space
        F;
        F_len;
        % lengths of the lattice vectors
        b_len;

        % Mode basis vectors for the irreducible representation of the mode.  It
        % is easiest to use the Projector object to define this.  This can
        % be two dimentional, in which case the coeffs property becomes
        % important.
        basis;
        % A scalling factor for the harmonic.  This is most important when
        % the mode is 2D, such as the E irreducible representation of C4v,
        % so that the coeffs will define the orientation of the mode within
        % the 2D space.
        coeffs;

        vec_type;

        % vector of Fourier coefficients
        harm;
    end

    properties (GetAccess='private', SetAccess='private')
        is_hexagonal = false;
        X;
        Y;
        K;
        X20;
        Y20;
        K20;
        vec_x;
        vec_y;
        vec_z;
        vec_x20;
        vec_y20;
        vec_z20;
        vec_norm;

        num_points;

        % General derivative operator that works on 3D matrix of fields.
        derivx = [];
        derivy = [];

    end
    
    methods
        function hObj = Harmonic(Hilbert_space, mode_basis, mode_coeffs, varargin)
            %HARMONIC Construct an instance of this class
            %   
            % Hilbert_space   Basis to use
            % mode_basis      Basis vectors of the mode.  Use the Projector
            %                 object for this.  Can be two dimentional.
            % mode_coeffs     Scaling factors for the mode
            % vec_type        "L", "M" or "N"
            % num_points      How fine a mesh to use
            % lattice_vecs    lengths of the lattice vectors

            switch nargin
                case 3
                    hObj.num_points = 100;
                    lattice_vecs = [1 1];
                    vec_type = "M";
                case 4
                    vec_type = varargin{1};
                    hObj.num_points = 100;
                    lattice_vecs = [1 1];
                case 5
                    vec_type = varargin{1};
                    hObj.num_points = varargin{2};
                    lattice_vecs = [1 1];
                otherwise
                    vec_type = varargin{1};
                    hObj.num_points = varargin{2};
                    lattice_vecs = varargin{3};
            end
            % if even make odd so that there is a zero frequency at the
            % center of the Fourier space
            if mod(hObj.num_points,2) == 0 
                hObj.num_points = hObj.num_points + 1;
            end

            % if Hilbert_space is a name of a standard F_i space
            if isstring(Hilbert_space)
                for F_i=Hilbert_space
                    switch lower(F_i)
                        case "f0"
                            hObj.F = [hObj.F; {[0; 0]}];
                        case "f1"
                            hObj.F = [hObj.F; {[1; 0]; [0; 1]; [-1; 0]; [0; -1]}];
                        case "f2"
                            hObj.F = [hObj.F; {[1; 1]; [-1; 1]; [-1; -1]; [1; -1]}];
                        case "f3"
                            hObj.F = [hObj.F; {[2; 0]; [0; 2]; [-2; 0]; [0; -2]}];
                        case "f4"
                            hObj.F = [hObj.F; {[2; 1]; [1; 2]; [-1; 2]; [-2; 1]; ...
                                [-2; -1]; [-1; -2]; [1; -2]; [2; -1]}];
                        case "f0_hex"
                            hObj.F = [hObj.F; {[0; 0]}];
                            hObj.is_hexagonal = true;
                        case "f1_hex"
                            hObj.F = [hObj.F; {[1/2; sqrt(3)/2]; [1; 0]; [1/2; -sqrt(3)/2]; ...
                                [-1/2; -sqrt(3)/2]; [-1; 0]; [-1/2; sqrt(3)/2]}];
                            %hObj.F = [hObj.F; {[1; 0]; [1/2; sqrt(3)/2]; [-1/2; sqrt(3)/2]; ...
                            %    [-1; 0]; [-1/2; -sqrt(3)/2]; [1/2; -sqrt(3)/2]}];
                            hObj.is_hexagonal = true;
                        case "f2_hex"
                            hObj.F = [hObj.F; {[3/2; sqrt(3)/2]; [0; sqrt(3)]; [-3/2; sqrt(3)/2]; ...
                                [-3/2; -sqrt(3)/2]; [0; -sqrt(3)]; [3/2; -sqrt(3)/2]}];
                            hObj.is_hexagonal = true;
                        case "f1_tri"
                            %hObj.F = [hObj.F; {[0; sqrt(3)/2]; [3/4; -sqrt(3)/4]; [-3/4; -sqrt(3)/4]}];
                            %hObj.F = [hObj.F; {[0; 1]; [sqrt(3)/2; -1/2]; [-sqrt(3)/2; -1/2]}];
                            hObj.F = [hObj.F; {[-1/2; sqrt(3)/2]; [1; 0]; [-1/2; -sqrt(3)/2]}];
                            %hObj.F = [hObj.F; {[-1/2; sqrt(3)/2]; [1/2; sqrt(3)/2]; [1;0]; [1/2; -sqrt(3)/2]; ...
                            %    [-1/2; -sqrt(3)/2]; [-1; 0]}];
                            hObj.is_hexagonal = true;
                        otherwise
                            error("Unknown Hibert space name");
                    end
                end
                hObj.F = single([hObj.F{:}].');
                hObj.F_len = length(hObj.F);
            elseif isnumeric(Hilbert_space) % array of values
                sz = size(Hilbert_space);
                if sz(2) == 2
                    hObj.F_len = sz(1);
                    hObj.F = Hilbert_space;
                else
                    hObj.F_len = sz(1)*sz(2);
                    hObj.F = reshape(int16(Hilbert_space),[hObj.F_len 2]);
                end
            else % Hilbert_space is a cell array of [m,n]
                if isrow(Hilbert_space)
                    hObj.F = int16([Hilbert_space{:}].');
                else
                    hObj.F = int16([Hilbert_space{:}]);
                end
                hObj.F_len = length(hObj.F);
            end


            hObj.b_len = 2*pi./lattice_vecs;
            hObj.harm = mode_basis.*mode_coeffs;
            hObj.basis = mode_basis;
            if isrow(mode_coeffs)
                hObj.coeffs = mode_coeffs.';
            else
                hObj.coeffs = mode_coeffs;
            end

            for i=1:length(hObj.F)
                hObj.derivx = cat(3, hObj.derivx, 1i*single(hObj.F(i,1)));
                hObj.derivy = cat(3, hObj.derivy, 1i*single(hObj.F(i,2)));
            end

            [hObj.X,hObj.Y,hObj.K,hObj.X20,hObj.Y20,hObj.K20] = hObj.plot_setup;

            hObj.vec_type = vec_type;
            switch lower(vec_type)
                case "l"
                    hObj.vec_y = sum(hObj.K.*hObj.derivy,3);
                    hObj.vec_x = sum(hObj.K.*hObj.derivx,3);
                    hObj.vec_y20 = sum(hObj.K20.*hObj.derivy,3);
                    hObj.vec_x20 = sum(hObj.K20.*hObj.derivx,3);
                    hObj.vec_z = zeros(size(hObj.vec_x));
                    hObj.vec_z20 = zeros(size(hObj.vec_x20));
                case "m"
                    hObj.vec_x = sum(hObj.K.*hObj.derivy,3);
                    hObj.vec_y = -sum(hObj.K.*hObj.derivx,3);
                    hObj.vec_x20 = sum(hObj.K20.*hObj.derivy,3);
                    hObj.vec_y20 = -sum(hObj.K20.*hObj.derivx,3);
                    hObj.vec_z = zeros(size(hObj.vec_x));
                    hObj.vec_z20 = zeros(size(hObj.vec_x20));
                case"n"
                    hObj.vec_z = sum(-(sum(hObj.K.*hObj.derivx.*hObj.derivx+hObj.K.*hObj.derivy.*hObj.derivy,3))...
                        /sqrt(hObj.b_len(1)^2+hObj.b_len(2)^2),3);
                    hObj.vec_z20 = -(sum(hObj.K20.*hObj.derivx.*hObj.derivx+hObj.K20.*hObj.derivy.*hObj.derivy,3))...
                        /sqrt(hObj.b_len(1)^2+hObj.b_len(2)^2);
                    hObj.vec_x = zeros(size(hObj.vec_z));
                    hObj.vec_x20 = zeros(size(hObj.vec_z20));
                    hObj.vec_y = zeros(size(hObj.vec_z));
                    hObj.vec_y20 = zeros(size(hObj.vec_z20));
                otherwise
                    error("Harmonic.m: Unknown type of vector planar harmonic");
            end
            hObj.vec_norm = sqrt(real(hObj.vec_x).^2+real(hObj.vec_y).^2+real(hObj.vec_z).^2);
        end
        
        % return the harmonic as a ket
        function out = ket_harm(hObj)
            out = sum(hObj.harm,2);
        end
        
        % return the harmonic as a bra (conj. transpose)
        function out = bra_harm(hObj)
            out = sum(hObj.harm,2)';
        end

        % return the x and y components of the vector planar harmonic as
        % kets
        function [x,y,z,xpoints,ypoints] = ket_vec(hObj)
            x = hObj.vec_x;
            y = hObj.vec_y;
            z = hObj.vec_z;
            xpoints = hObj.X;
            ypoints = hObj.Y;
        end

        % return the x and y components of the vector planar harmonic as
        % bras (conj. transpose)
        function [x,y,z,xpoints,ypoints] = bra_vec(hObj)
            x = hObj.vec_x';
            y = hObj.vec_y';
            z = hObj.vec_z';
            xpoints = hObj.X;
            ypoints = hObj.Y;
        end

        % return the x and y components of the vector planar harmonic as
        % kets
        function [x,y,xpoints,ypoints] = ket_vec20(hObj)
            x = hObj.vec_x20;
            y = hObj.vec_y20;
            xpoints = hObj.X20;
            ypoints = hObj.Y20;
        end

        % return the x and y components of the vector planar harmonic as
        % bras (conj. transpose)
        function [x,y,xpoints,ypoints] = bra_vec20(hObj)
            x = hObj.vec_x20';
            y = hObj.vec_y20';
            xpoints = hObj.X20;
            ypoints = hObj.Y20;
        end

        % The harmonic is in Fourier space, so just construct a square in
        % Fourier space with the harmonic values saved to the corresponding
        % point (m,n) as defined in hObj.F and then take the Fourier
        % transform of it to get the real space image.
        %
        % This is not usefull for the case of C_3 or C_6v etc since taking
        % a FFT of these is harder.  I am therefore no longer using this
        % method.
        function I = get_phi(hObj)
            A = zeros(hObj.num_points);
            H = hObj.basis*hObj.coeffs;
            center = (hObj.num_points-1)/2+1;
            for i=1:hObj.F_len
                A(center+hObj.F(i,1),center+hObj.F(i,2)) = ...
                    A(center+hObj.F(i,1),center+hObj.F(i,2))+H(i);
            end
            I=ifft2(ifftshift(-A*hObj.num_points^2));
            % The results if off center by 1 pixel.
            I = I(2:end,2:end);
        end

        % plot the scalar planar harmonic
        function I = plot_phi(hObj)

            figure;
            I = sum(hObj.K,3);%get_phi;
            imagesc([hObj.X(1) hObj.X(end)],[hObj.Y(1) hObj.Y(end)],real(I))
            axis tight;
            pbaspect([1 1 1])
            view([0 90])
            box on
            set(gca, 'FontSize', 16);
            xlabel('x/a', 'FontSize',18)
            ylabel('y/a', 'FontSize',18)
            set(gca, 'XTick', [-hObj.b_len(1)/2 0 hObj.b_len(1)/2])
            set(gca, 'XTickLabel', {'-a','0','a'})
            set(gca, 'YTick', [-hObj.b_len(2)/2 0 hObj.b_len(2)/2])
            set(gca, 'YTickLabel', {'-a','0','a'})
            if hObj.is_hexagonal == true
                xlim([-hObj.b_len(1) hObj.b_len(1)])
                ylim([-hObj.b_len(2) hObj.b_len(2)])
                plot([hObj.F(:,1); hObj.F(1,1)]*hObj.b_len(1)*2/3, ...
                    [hObj.F(:,2); hObj.F(1,2)]*hObj.b_len(2)*2/3, 'k', 'LineWidth', 2)
            else
                xlim([-hObj.b_len(1)/2 hObj.b_len(1)/2])
                ylim([-hObj.b_len(2)/2 hObj.b_len(2)/2])
            end
        end

        % Plot the vector planar harmonic
        % TODO: set pbaspect() via the lattice vectors
        function plot(hObj)
            figure;
            surf(hObj.X,hObj.Y,zeros(size(hObj.X)),hObj.vec_norm-max(hObj.vec_norm,[],'all'), 'EdgeColor','none')
            hold on;
            view([0 90])
            xlabel('$x_1$', 'FontSize', 18, 'Interpreter', 'latex')
            ylabel('$x_2$', 'FontSize', 18, 'Interpreter', 'latex')
            set(gca, 'XTick', [-hObj.b_len(1)/2 0 hObj.b_len(1)/2])
            set(gca, 'XTickLabel', {'-a','0','a'})
            set(gca, 'YTick', [-hObj.b_len(2)/2 0 hObj.b_len(2)/2])
            set(gca, 'YTickLabel', {'-a','0','a'})
            set(gca, 'LineWidth', 2)
            set(gca, 'FontSize', 16)
            grid off
            box on
            if hObj.is_hexagonal == true
                xlim([-hObj.b_len(1) hObj.b_len(1)])
                ylim([-hObj.b_len(2) hObj.b_len(2)])
                plot([1 1/2 -1/2 -1 -1/2 1/2 1]*hObj.b_len(1)*2/3, ...
                    [0 sqrt(3)/2 sqrt(3)/2 0 -sqrt(3)/2 -sqrt(3)/2 0]*hObj.b_len(2)*2/3, 'k', 'LineWidth', 2)
                %plot([hObj.F(:,1); hObj.F(1,1)]*hObj.b_len(1)*2/3, ...
                %    [hObj.F(:,2); hObj.F(1,2)]*hObj.b_len(2)*2/3, 'k', 'LineWidth', 2)
            else
                xlim([-hObj.b_len(1)/2 hObj.b_len(1)/2])
                ylim([-hObj.b_len(2)/2 hObj.b_len(2)/2])
            end

            quiver3(hObj.X20,hObj.Y20,zeros(size(hObj.X20)),real(hObj.vec_x20), ...
                real(hObj.vec_y20),real(hObj.vec_z20), 'k')
            pbaspect([1 1 1])
            if lower(hObj.vec_type) == "n"
                view([-15 30])
            else
                view([0 90])
            end
        end
    end

    methods (Access = private)
        % create the mesh grid and K vector
        function [X,Y,K,vX,vY,vK] = plot_setup(hObj)
            [X,Y,K,vX,vY,vK] = hObj.plot_setup_k([0 0]);
        end

        function [X,Y,K,vX,vY,vK] = plot_setup_k(hObj, k)
            if hObj.is_hexagonal == true
                [X,Y] =   meshgrid((-hObj.b_len(1):2*hObj.b_len(1)/(hObj.num_points-1):hObj.b_len(1)),...
                    (-hObj.b_len(2):2*hObj.b_len(2)/(hObj.num_points-1):hObj.b_len(2)));
                % lower resolution mesh for the vector plots
                [vX,vY] = meshgrid((-hObj.b_len(1):2*hObj.b_len(1)/(hObj.num_points/5-1):hObj.b_len(1)),...
                    (-hObj.b_len(2):2*hObj.b_len(2)/(hObj.num_points/5-1):hObj.b_len(2)));
            else
                [X,Y] =   meshgrid((-hObj.b_len(1)/2:hObj.b_len(1)/(hObj.num_points-1):hObj.b_len(1)/2),...
                    (-hObj.b_len(2)/2:hObj.b_len(2)/(hObj.num_points-1):hObj.b_len(2)/2));
                % lower resolution mesh for the vector plots
                [vX,vY] = meshgrid((-hObj.b_len(1)/2:hObj.b_len(1)/(hObj.num_points/5-1):hObj.b_len(1)/2),...
                    (-hObj.b_len(2)/2:hObj.b_len(2)/(hObj.num_points/5-1):hObj.b_len(2)/2));
            end


            K = [];
            vK = [];
            for i = 1:length(hObj.F)
                b = single(hObj.F(i,:));
                K =  cat(3,  K, sum(hObj.harm(i,:))*exp(1i*(b(1)+k(1))*X +1i*(b(2)+k(2))*Y));
                vK = cat(3, vK, sum(hObj.harm(i,:))*exp(1i*(b(1)+k(1))*vX+1i*(b(2)+k(2))*vY));
            end
        end
    end
end

