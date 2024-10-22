classdef Field3D < handle
    % Field
    % Handles the orientation etc of vector fields for use with projection
    % operators.

    properties
        x;
        y;
        z;
        Ex;
        Ey;
        Ez;
        size_x;
        size_y;
        size_z;
        % The Fourier space used to symmetrise a field
        F_x = [];
        F_y = [];
        F_z = [];
    end

    methods
        function obj = Field3D(x,y,z,Ex,Ey,Ez, size_x, size_y, size_z)
            obj.size_x = size_x;
            obj.size_y = size_y;
            obj.size_z = size_z;

            if iscolumn(x) || isrow(x)
                obj.x = reshape(x, [size_x size_y size_z]);
            else
                obj.x = x;
            end
            if iscolumn(y) || isrow(y)
                obj.y = reshape(y, [size_x size_y size_z]);
            else
                obj.y = y;
            end
            if iscolumn(z) || isrow(z)
                obj.z = reshape(z, [size_x size_y size_z]);
            else
                obj.z = z;
            end
            if iscolumn(Ex) || isrow(Ex)
                obj.Ex = reshape(Ex, [size_x size_y size_z]);
            else
                obj.Ex = Ex;
            end
            if iscolumn(Ey) || isrow(Ey)
                obj.Ey = reshape(Ey, [size_x size_y size_z]);
            else
                obj.Ey = Ey;
            end
            if iscolumn(Ez) || isrow(Ez)
                obj.Ez = reshape(Ez, [size_x size_y size_z]);
            else
                obj.Ez = Ez;
            end

            % correct the orientation of the array elements so that they
            % match Cartesian coordiantes.
            if obj.x(1,1,1) == obj.x(1,end,1)
                obj.x = rot90(obj.x);
                obj.y = rot90(obj.y);
                obj.z = rot90(obj.z);
                obj.Ex = rot90(obj.Ex);
                obj.Ey = rot90(obj.Ey);
                obj.Ez = rot90(obj.Ez);
            end
            if obj.y(1,1,1) < obj.y(end,1,1)
                obj.x = flipud(obj.x);
                obj.y = flipud(obj.y);
                obj.z = flipud(obj.z);
                obj.Ex = flipud(obj.Ex);
                obj.Ey = flipud(obj.Ey);
                obj.Ez = flipud(obj.Ez);
            end
            if obj.x(1,1,1) > obj.x(1,end,1)
                obj.x = fliplr(obj.x);
                obj.y = fliplr(obj.y);
                obj.z = fliplr(obj.z);
                obj.Ex = fliplr(obj.Ex);
                obj.Ey = fliplr(obj.Ey);
                obj.Ez = fliplr(obj.Ez);
            end
            if obj.z(1,1,1) < obj.z(1,1,end)
                obj.x = flip(obj.x,3);
                obj.y = flip(obj.y,3);
                obj.z = flip(obj.z,3);
                obj.Ex = flip(obj.Ex,3);
                obj.Ey = flip(obj.Ey,3);
                obj.Ez = flip(obj.Ez,3);
            end
        end

        % This fuction converts a vector field into a form that can be tested for
        % symmetry.  It changes the basis for each point so that the x-axis lies
        % along the line from the origin to that particular point.
        %
        % The results are returned in matrices of the same form as that passed
        % to the function, so you will have to use outx(:) and outy(:) when
        % multiplying by a projector.
        %
        % Since the x and y components are generated seperately the test for
        % symmetry is
        % (<x|P|x> + <y|P'|y>)/(<x|x>+<y|y>)
        % which is implemented by the function norm_field().  P is the projection
        % operator for the symmetry and P' is that with the characters of the
        % mirror symmetries negated.  That is, if P is C4v A1 then P' is C4v A2 and
        % if P is C4v B2 then P' is C4v B1.  If there are no mirror symmetries then
        % P=P'.  This comes about because the y axis is perpendicular the vectors
        % from the origin to the point, so mirror symmetries are reversed for y
        % components.

        % x   array of x components of the vector field.  They should be arranged
        %     in the same order as the cartesian coordiates so that the angle from
        %     the origin to each point can be determined by its location in the
        %     array.
        % y   array of y components of the vector field.
        %
        % outx, outy are the new values with the adjusted coordinates.
        % X and Y are returned so that you can check the orientation used in this
        % function.
        function out = symmetrise(obj)

            if isodd(obj.size_x)
                x_vals = -(obj.size_x-1)/2:(obj.size_x-1)/2;
            else
                x_vals = -obj.size_x+1:2:obj.size_x-1;
            end
            if isodd(obj.size_y)
                y_vals = (obj.size_y-1)/2:-1:-(obj.size_y-1)/2;
            else
                y_vals = obj.size_y-1:-2:-obj.size_y+1;
            end
            if isodd(obj.size_z)
                z_vals = (obj.size_z-1)/2:-1:-(obj.size_z-1)/2;
            else
                z_vals = obj.size_z-1:-2:-obj.size_z+1;
            end
            [X,Y,Z] = meshgrid(x_vals, y_vals, z_vals);

            % WrapToPi() makes the results easier to analyse.
            phi_ref = wrapToPi(angle(X+1i*Y));
            theta_ref = wrapToPi(angle(sqrt(X.^2+Y.^2)+1i*Z));

            % rotate about the z axes
            xsym = obj.Ex.*cos(phi_ref) + obj.Ey.*sin(phi_ref);
            ysym = -obj.Ex.*sin(phi_ref) + obj.Ey.*cos(phi_ref);
            zsym = obj.Ez;
            % now rotate about the y' axis
            xsym = xsym.*cos(theta_ref) + zsym.*sin(theta_ref);
            zsym = -xsym.*sin(theta_ref) + zsym.*cos(theta_ref);

            

            out = Field3D(obj.x, obj.y, obj.z, xsym, ysym, zsym, obj.size_x, obj.size_y, obj.size_z);
            out.F_x = X;
            out.F_y = Y;
            out.F_z = Z;
        end

        function out = norm(obj)
            out = sqrt(obj.Ex.*conj(obj.Ex)+obj.Ey.*conj(obj.Ey)+obj.Ez.*conj(obj.Ez));
        end

        function out = phase(obj)
            out = angle(real(obj.Ex)+1i*real(obj.Ey));
        end

        % This function is used to determine if a vector field has a particular
        % symmetry.
        %
        % This function returns
        %
        % <bra_x|oprtr|ket_x> + <bra_y|oprtr_complement|ket_y> + <bra_z|oprtr|ket_z>
        % --------------------------------------------------------------------------
        %              <bra_x|ket_x> + <bra_y|ket_y> + <bra_z|ket_z>
        %
        % where ket_x and ket_y are vectors, bra_x=ket_x' and both oprtr and
        % oprtr_complement are matrices of the correct size. oprtr is the
        % projection operator for the symmetry and oprtr_complement is that with
        % the characters of the mirror symmetries negated.  That is, if oprtr is
        % C4v A1 then oprtr_complement is C4v A2 and if oprtr is C4v B2 then
        % oprtr_complement is C4v B1.  If there are no mirror symmetries then
        % oprtr=oprtr_complement.  This comes about because the y axis is
        % perpendicular the vectors from the origin to the given point, so mirror
        % symmetries are reversed for y components.

        function out = projection(obj,oprtr, oprtr_complement)
            if isempty(obj.F_x) || isempty(obj.F_y)
                warning('Field is not symmetrised!');
            end

            if size(oprtr,2) ~= length(obj.Ex(:))
                error(['projection(): First argument does not have the same number ' ...
                    'of columns as the number of elements in the Ex field.']);
            end
            if size(oprtr_complement,2) ~= length(obj.Ey(:))
                error(['norm_vec(): Second argument does not have the same number ' ...
                    'of columns as the number of elements in the Ey field.']);
            end
            if size(oprtr_complement,2) ~= length(obj.Ez(:))
                error(['norm_vec(): Second argument does not have the same number ' ...
                    'of columns as the number of elements in the Ez field.']);
            end

            % need to avoid overflow errors with using single precision for memory
            % saving.
            out = (double(obj.Ex(:))'*double(oprtr)*double(obj.Ex(:))+ ...
                double(obj.Ey(:))'*double(oprtr_complement)*double(obj.Ey(:))+ ...
                double(obj.Ez(:))'*double(oprtr_complement)*double(obj.Ez(:)))/ ...
                (double(obj.Ex(:))'*double(obj.Ex(:))+double(obj.Ey(:))'*double(obj.Ey(:))+ ...
                +double(obj.Ez(:))'*double(obj.Ez(:)));
        end

        function plot_axis(obj, ax, varargin)

            if nargin == 2
                step_size = varargin{1};
            elseif nargin > 2
                step_size = varargin{1};
                quiver_varargin = varargin(2:end);
            else
                step_size = 1;
                quiver_varargin = {'k'};
            end

            norm2d = obj.norm;

            for z_plane_num = 1:length(unique(obj.z))
                surf(ax, squeeze(obj.x(:,:,z_plane_num))*1e9, squeeze(obj.y(:,:,z_plane_num))*1e9, ...
                    squeeze(obj.z(:,:,z_plane_num))*1e9, ...
                    squeeze(norm2d(:,:,z_plane_num)), 'EdgeColor','none');
                hold(ax, 'on');
            end

            set(ax, 'FontSize', 14)
            xlabel(ax, 'x (nm)', 'FontSize', 18)
            ylabel(ax, 'y (nm)', 'FontSize', 18)
            set(ax, 'XTick', fix([obj.x(1,1) 0 obj.x(1,end)]*1e9))
            set(ax, 'YTick', fix([obj.y(end,1) 0 obj.y(1,1)]*1e9))
            set(ax, 'LineWidth', 2)
            grid(ax, 'off');
            box(ax, 'on');
            xlim(ax, [obj.x(1,1) obj.x(1,end)]*1e9)
            ylim(ax, [obj.y(end,1) obj.y(1,1)]*1e9)
            title(ax, '');

            for z_plane_num = 1:length(unique(obj.z))
                Ex2D = squeeze(obj.Ex(:,:,z_plane_num));
                Ey2D = squeeze(obj.Ey(:,:,z_plane_num));
                Ez2D = squeeze(obj.Ez(:,:,z_plane_num));
                x2D = squeeze(obj.x(:,:,z_plane_num));
                y2D = squeeze(obj.y(:,:,z_plane_num));
                z2D = squeeze(obj.z(:,:,z_plane_num));
                quiver3(ax, x2D(1:step_size:end)*1e9, y2D(1:step_size:end)*1e9, z2D(1:step_size:end)*1e9, ...
                    real(Ex2D(1:step_size:end)), ...
                    real(Ey2D(1:step_size:end)), real(Ez2D(1:step_size:end)), ...
                    quiver_varargin{:})
            end

            pbaspect(ax, [1 1 1])
            %view(ax, [0 90])

            cb = colorbar(ax);
            ylabel(cb, '|E| (V/m)', 'FontSize', 14);

        end

        function plot(obj, varargin)
            % Create an input parser object
            p = inputParser;
        
            % Define parameters and their default values
            %addRequired(p, 'x', @isnumeric);  % x is a required numeric parameter
            addParameter(p, 'step_size', 1, @isnumeric);  % Optional parameter with default value
            addParameter(p, 'quiver_varargin', {'k'}, @iscell);  % Optional parameter with default value
            addParameter(p, 'z', 0, @isnumeric);  % Optional parameter with default value
        
            % Parse input arguments
            %parse(p, x, vararg{:});
            parse(p, varargin{:});
        
            % Access the parsed inputs
            %x = p.Results.x;
            step_size = p.Results.step_size;
            quiver_varargin = p.Results.quiver_varargin;
            z_plane = p.Results.z;
    %{
            if nargin == 2
                step_size = varargin{1};
            elseif nargin > 2
                step_size = varargin{1};
                quiver_varargin = varargin(2:end);
            else
                step_size = 1;
                quiver_varargin = {'k'};
            end
    %}
            norm2d = obj.norm;

            figure;
            if z_plane == 0
                for z_plane_num = 1:length(unique(obj.z))
                    surf(squeeze(obj.x(:,:,z_plane_num)), squeeze(obj.y(:,:,z_plane_num)), ...
                        squeeze(obj.z(:,:,z_plane_num)), ...
                        squeeze(norm2d(:,:,z_plane_num)), 'EdgeColor','none')
                end
            else
                surf(squeeze(obj.x(:,:,z_plane)), squeeze(obj.y(:,:,z_plane)), ...
                    squeeze(obj.z(:,:,z_plane)), ...
                    squeeze(norm2d(:,:,z_plane)), 'EdgeColor','none')
            end
            hold on;

            set(gca, 'FontSize', 20)
            %view([0 90])
            xlabel('x', 'FontSize', 24)
            ylabel('y', 'FontSize', 24)
            set(gca, 'XTick', [obj.x(1,1) 0 obj.x(1,end)])
            set(gca, 'YTick', [obj.y(end,1) 0 obj.y(1,1)])
            set(gca, 'LineWidth', 2)
            grid off
            box on
            xlim([obj.x(1,1) obj.x(1,end)])
            ylim([obj.y(end,1) obj.y(1,1)])

            if z_plane == 0
                for z_plane_num = 1:length(unique(obj.z))
                    Ex2D = squeeze(obj.Ex(:,:,z_plane_num));
                    Ey2D = squeeze(obj.Ey(:,:,z_plane_num));
                    Ez2D = squeeze(obj.Ez(:,:,z_plane_num));
                    x2D = squeeze(obj.x(:,:,z_plane_num));
                    y2D = squeeze(obj.y(:,:,z_plane_num));
                    z2D = squeeze(obj.z(:,:,z_plane_num));
                    quiver3(x2D(1:step_size:end), y2D(1:step_size:end), z2D(1:step_size:end), ...
                        real(Ex2D(1:step_size:end)), ...
                        real(Ey2D(1:step_size:end)), real(Ez2D(1:step_size:end)), ...
                        quiver_varargin{:})
                end
            else
                Ex2D = squeeze(obj.Ex(:,:,z_plane));
                Ey2D = squeeze(obj.Ey(:,:,z_plane));
                Ez2D = squeeze(obj.Ez(:,:,z_plane));
                x2D = squeeze(obj.x(:,:,z_plane));
                y2D = squeeze(obj.y(:,:,z_plane));
                z2D = squeeze(obj.z(:,:,z_plane));
                quiver3(x2D(1:step_size:end), y2D(1:step_size:end), z2D(1:step_size:end), ...
                    real(Ex2D(1:step_size:end)), ...
                    real(Ey2D(1:step_size:end)), real(Ez2D(1:step_size:end)), ...
                    quiver_varargin{:})
            end
            pbaspect([1 1 1])
            %view([0 90])
        end

        function plot_norm(obj)

            norm2d = obj.norm;

            figure;
            for z_plane_num = 1:length(unique(obj.z))
                surf(squeeze(obj.x(:,:,z_plane_num)), squeeze(obj.y(:,:,z_plane_num)), ...
                    squeeze(obj.z(:,:,z_plane_num)), squeeze(norm2d(:,:,z_plane_num)), 'EdgeColor','none')
                hold on;
            end
            %view([0 90])
            xlabel('$x_1$', 'FontSize', 18, 'Interpreter', 'latex')
            ylabel('$x_2$', 'FontSize', 18, 'Interpreter', 'latex')
            zlabel('$x_3$', 'FontSize', 18, 'Interpreter', 'latex')
            set(gca, 'XTick', [obj.x(1,1,1) 0 obj.x(1,end,1)])
            set(gca, 'YTick', [obj.y(end,1,1) 0 obj.y(1,1,1)])
            set(gca, 'ZTick', [obj.z(1,1,end) 0 obj.z(1,1,1)])
            set(gca, 'LineWidth', 2)
            set(gca, 'FontSize', 16)
            grid off
            box on
            xlim([obj.x(1,1,1) obj.x(1,end,1)])
            ylim([obj.y(end,1,1) obj.y(1,1,1)])
            zlim([obj.z(1,1,end) obj.z(1,1,1)])
            pbaspect([1 1 1])
            %view([0 90])
        end

        % plot only a single component of the electric field.
        function plot_component(obj, component)
            
            switch lower(component)
                case "x"
                    vec = obj.Ex;
                case "y"
                    vec = obj.Ey;
                case "z"
                    vec = obj.Ez;
                otherwise
                    error("Component must be either ""x"", ""y"" or ""z""");
            end

            figure;
            for z_plane_num = 1:length(unique(obj.z))
                surf(squeeze(obj.x(:,:,z_plane_num)), squeeze(obj.y(:,:,z_plane_num)), ...
                    squeeze(obj.z(:,:,z_plane_num)), ...
                    real(squeeze(vec(:,:,z_plane_num))), 'EdgeColor','none')
                hold on;
            end

            set(gca, 'FontSize', 18)
            %view([0 90])
            xlabel('x', 'FontSize', 20)
            ylabel('y', 'FontSize', 20)
            set(gca, 'XTick', [obj.x(1,1) 0 obj.x(1,end)])
            set(gca, 'YTick', [obj.y(end,1) 0 obj.y(1,1)])
            set(gca, 'LineWidth', 2)
            grid off
            box on
            xlim([obj.x(1,1) obj.x(1,end)])
            ylim([obj.y(end,1) obj.y(1,1)])
            pbaspect([1 1 1])
            %view([0 90])

        end

        function quiver(obj, varargin)

            if nargin == 2
                scale = varargin{1};
            else
                scale = 1; 
            end
            figure;

            for z_plane_num = 1:length(unique(obj.z))
                quiver3(squeeze(obj.x(:,:,z_plane_num)), squeeze(obj.y(:,:,z_plane_num)), ...
                    squeeze(obj.z(:,:,z_plane_num)), real(squeeze(obj.Ex(:,:,z_plane_num))), ...
                    real(squeeze(obj.Ey(:,:,z_plane_num))), real(squeeze(obj.Ez(:,:,z_plane_num))), ...
                    scale, 'k');
                hold on;
            end
            xlabel('$x_1$', 'FontSize', 18, 'Interpreter', 'latex')
            ylabel('$x_2$', 'FontSize', 18, 'Interpreter', 'latex')
            set(gca, 'XTick', [obj.x(1,1) 0 obj.x(1,end)])
            set(gca, 'YTick', [obj.y(end,1) 0 obj.y(1,1)])
            set(gca, 'LineWidth', 2)
            set(gca, 'FontSize', 16)
            grid off
            box on
            xlim([obj.x(1,1) obj.x(1,end)])
            ylim([obj.y(end,1) obj.y(1,1)])
            axis equal
        end

        function plot_phase(obj, varargin)

            %{
            if nargin == 2
                step_size = varargin{1};
            elseif nargin > 2
                step_size = varargin{1};
                quiver_varargin = varargin(2:end);
            else
                step_size = 1;
                quiver_varargin = {'k'};
            end
            %}

            [Fx,Fy] = gradient(obj.phase);
            phase2D = sqrt(Fx.^2+Fy.^2);

            figure;
            for z_plane_num = 1:length(unique(obj.z))
                surf(squeeze(obj.x(:,:,z_plane_num)), squeeze(obj.y(:,:,z_plane_num)), ...
                    squeeze(obj.z(:,:,z_plane_num)), squeeze(phase2D(:,:,z_plane_num)), ...
                    'EdgeColor','none')
                hold on;
            end

            hold on;
            set(gca, 'FontSize', 20)
            %view([0 90])
            xlabel('x', 'FontSize', 24)
            ylabel('y', 'FontSize', 24)
            set(gca, 'XTick', [obj.x(1,1) 0 obj.x(1,end)])
            set(gca, 'YTick', [obj.y(end,1) 0 obj.y(1,1)])
            set(gca, 'LineWidth', 2)
            grid off
            box on
            xlim([obj.x(1,1) obj.x(1,end)])
            ylim([obj.y(end,1) obj.y(1,1)])

            %quiver3(obj.x(1:step_size:end),obj.y(1:step_size:end),...
            %    zeros(size(obj.y(1:step_size:end))),real(obj.Ex(1:step_size:end)),...
            %    real(obj.Ey(1:step_size:end)),real(obj.Ez(1:step_size:end)), quiver_varargin{:})
            pbaspect([1 1 1])
            %view([0 90])
        end

        % Apply a projector to the field and return the resultant field as
        % a new Field object.  That is, return: |psi'> = P|psi>
        function out = get_projection(obj, oprtr, oprtr_complement)

            sym_field = obj.symmetrise;
            
            % Check that the projectors are of the right size
            if size(oprtr,2) ~= length(sym_field.Ex(:))
                error(['projection(): First argument does not have the same number ' ...
                    'of columns as the number of elements in the Ex field.']);
            end
            if size(oprtr_complement,2) ~= length(sym_field.Ey(:))
                error(['norm_vec(): Second argument does not have the same number ' ...
                    'of columns as the number of elements in the Ey field.']);
            end
            if size(oprtr_complement,2) ~= length(sym_field.Ez(:))
                error(['norm_vec(): Second argument does not have the same number ' ...
                    'of columns as the number of elements in the Ez field.']);
            end

            % apply the projectors to the field components.
            % need to avoid overflow errors with using single precision for memory
            % saving.
            newEx = double(oprtr)*double(sym_field.Ex(:));
            newEx = reshape(newEx, [obj.size_x obj.size_y obj.size_z]);
            newEy = double(oprtr_complement)*double(sym_field.Ey(:));
            newEy = reshape(newEy, [obj.size_x obj.size_y obj.size_z]);
            newEz = double(oprtr_complement)*double(sym_field.Ez(:));
            newEz = reshape(newEz, [obj.size_x obj.size_y obj.size_z]);

            % Desymmetrise the field
            if isodd(obj.size_x)
                x_vals = -(obj.size_x-1)/2:(obj.size_x-1)/2;
            else
                x_vals = -obj.size_x+1:2:obj.size_x-1;
            end
            if isodd(obj.size_y)
                y_vals = (obj.size_y-1)/2:-1:-(obj.size_y-1)/2;
            else
                y_vals = obj.size_y-1:-2:-obj.size_y+1;
            end
            if isodd(obj.size_z)
                z_vals = (obj.size_z-1)/2:-1:-(obj.size_z-1)/2;
            else
                z_vals = obj.size_z-1:-2:-obj.size_z+1;
            end
            [X,Y,Z] = meshgrid(x_vals, y_vals, z_vals);

            % WrapToPi() makes the results easier to analyse.
            phi_ref = -wrapToPi(angle(X+1i*Y));
            theta_ref = -wrapToPi(angle(sqrt(X.^2+Y.^2)+1i*Z));

            % rotate about z axis
            xsym = newEx.*cos(phi_ref) + newEy.*sin(phi_ref);
            ysym = -newEx.*sin(phi_ref) + newEy.*cos(phi_ref);
            zsym = newEz;
            % now rotate about the y' axis
            xsym = xsym.*cos(theta_ref) + zsym.*sin(theta_ref);
            zsym = -xsym.*sin(theta_ref) + zsym.*cos(theta_ref);

            % Create Field object
            out = Field3D(sym_field.x, sym_field.y, sym_field.z, xsym, ysym, zsym, ...
                sym_field.size_x, sym_field.size_y, sym_field.size_z);
        end

        % Shift the field to the left and down for cases where the highest
        % symmetry point is not at the center.  
        function out = shift(obj, varargin)

            switch nargin
                case 1
                    Ex_shifted = [obj.Ex(:,(obj.size_x+(obj.size_x-1)/2)/2+1:end-1,:) ...
                        obj.Ex(:,1:(obj.size_x+(obj.size_x-1)/2)/2+1,:)];
                    Ey_shifted = [obj.Ey(:,(obj.size_x+(obj.size_x-1)/2)/2+1:end-1,:) ...
                        obj.Ey(:,1:(obj.size_x+(obj.size_x-1)/2)/2+1,:)];
                    Ez_shifted = [obj.Ez(:,(obj.size_x+(obj.size_x-1)/2)/2+1:end-1,:) ...
                        obj.Ez(:,1:(obj.size_x+(obj.size_x-1)/2)/2+1,:)];
                    Ex_shifted = [Ex_shifted(((obj.size_y-1)/2-1)/2+1:end-1,:,:); ...
                        Ex_shifted(1:((obj.size_y-1)/2-1)/2+1,:,:)];
                    Ey_shifted = [Ey_shifted(((obj.size_y-1)/2-1)/2+1:end-1,:,:); ...
                        Ey_shifted(1:((obj.size_y-1)/2-1)/2+1,:,:)];
                    Ez_shifted = [Ez_shifted(((obj.size_y-1)/2-1)/2+1:end-1,:,:); ...
                        Ez_shifted(1:((obj.size_y-1)/2-1)/2+1,:,:)];
                otherwise
                    shift_x = varargin{1};
                    shift_y = varargin{2};
                    if shift_x < 0
                        shift_x = -shift_x;
                        Ex_shifted = [obj.Ex(:,end-shift_x+1:end-1,:) ...
                            obj.Ex(:,1:end-shift_x+1,:)];
                        Ey_shifted = [obj.Ey(:,end-shift_x+1:end-1,:) ...
                            obj.Ey(:,1:end-shift_x+1,:)];
                        Ez_shifted = [obj.Ez(:,end-shift_x+1:end-1,:) ...
                            obj.Ez(:,1:end-shift_x+1,:)];
                    elseif shift_x > 0
                        Ex_shifted = [obj.Ex(:,shift_x:end-1,:) ...
                            obj.Ex(:,1:shift_x,:)];
                        Ey_shifted = [obj.Ey(:,shift_x:end-1,:) ...
                            obj.Ey(:,1:shift_x,:)];
                        Ez_shifted = [obj.Ez(:,shift_x:end-1,:) ...
                            obj.Ez(:,1:shift_x,:)];
                    else
                        Ex_shifted = obj.Ex;
                        Ey_shifted = obj.Ey;
                        Ez_shifted = obj.Ez;
                    end
                    if shift_y < 0
                        shift_y = -shift_y;
                        Ex_shifted = [Ex_shifted(end-shift_y+1:end-1,:,:); ...
                            Ex_shifted(1:end-shift_y+1,:,:)];
                        Ey_shifted = [Ey_shifted(end-shift_y+1:end-1,:,:); ...
                            Ey_shifted(1:end-shift_y+1,:,:)];
                        Ez_shifted = [Ez_shifted(end-shift_y+1:end-1,:,:); ...
                            Ez_shifted(1:end-shift_y+1,:,:)];
                    elseif shift_y > 0
                        Ex_shifted = [Ex_shifted(shift_y:end-1,:,:); ...
                            Ex_shifted(1:shift_y,:,:)];
                        Ey_shifted = [Ey_shifted(shift_y:end-1,:,:); ...
                            Ey_shifted(1:shift_y,:,:)];
                        Ez_shifted = [Ez_shifted(shift_y:end-1,:,:); ...
                            Ez_shifted(1:shift_y,:,:)];
                    end
                    %Ex_shifted = [Ex_shifted(((obj.size_y-1)/2-1)/2+1:end-1,:); ...
                    %    Ex_shifted(1:((obj.size_y-1)/2-1)/2+1,:)];
                    %Ey_shifted = [Ey_shifted(((obj.size_y-1)/2-1)/2+1:end-1,:); ...
                    %    Ey_shifted(1:((obj.size_y-1)/2-1)/2+1,:)];
                    %Ez_shifted = [Ez_shifted(((obj.size_y-1)/2-1)/2+1:end-1,:); ...
                    %    Ez_shifted(1:((obj.size_y-1)/2-1)/2+1,:)];
            end


            out = Field3D(obj.x,obj.y,obj.z,Ex_shifted,Ey_shifted,Ez_shifted, ...
                obj.size_x, obj.size_y, obj.size_z);
        end

    end
end