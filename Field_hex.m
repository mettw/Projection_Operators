classdef Field_hex < Field
    % Field_hex
    % Handles the orientation etc of vector fields for use with projection
    % operators.
    %
    % This version does not store the data as a square array, so it
    % can be used with hexagonal and any general cut plane.  It should be
    % possible to merge this with the standard Field object, but I don't
    % have time to work on it at the moment.

    properties
        kx;
        ky;
    end

    methods
        function obj = Field_hex(x,y,Ex,Ey,Ez, kx, ky)
            obj = obj@Field(x,y,Ex,Ey,Ez, size(x,1), size(x,2));
            obj.kx = kx;
            obj.ky = ky;
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
            % use kx and ky instead
            %[X,Y] = meshgrid(-(obj.size_x-1)/2:(obj.size_x-1)/2, ...
            %    (obj.size_y-1)/2:-1:-(obj.size_y-1)/2);

            % WrapToPi() makes the results easier to analyse.
            %ang_ref = wrapToPi(angle(X+1i*Y));
            ang_ref = wrapToPi(angle(obj.x+1i*obj.y));


            xsym = obj.Ex.*cos(ang_ref) + obj.Ey.*sin(ang_ref);

            ysym = -obj.Ex.*sin(wrapToPi(ang_ref)) + obj.Ey.*cos(wrapToPi(ang_ref));

            out = Field_hex(obj.x, obj.y, xsym, ysym, obj.Ez, obj.kx, obj.ky);
            out.F_x = obj.x;
            out.F_y = obj.y;
        end

        function plot(obj, varargin)

            if nargin == 2
                step_size = varargin{1};
            elseif nargin > 2
                step_size = varargin{1};
                quiver_varargin = varargin(2:end);
            else
                step_size = 1;
                quiver_varargin = {'k'};
            end
            
            E_norm = obj.norm;
            I = ones(sqrt(length(obj.kx)), sqrt(length(obj.ky)))*NaN;
            for i = 1:length(obj.x)
                I((obj.kx - obj.x(i) < abs(obj.x(i)+1e-5)*1e-10 & ...
                    obj.kx - obj.x(i) > -abs(obj.x(i)+1e-5)*1e-10 & ...
                    obj.ky - obj.y(i) < abs(obj.y(i)+1e-5)*1e-10 & ...
                    obj.ky - obj.y(i) > -abs(obj.y(i)+1e-5)*1e-10)) ...
                    = E_norm(i);
            end

            figure;
            surf(reshape(obj.kx, size(I)), reshape(obj.ky, size(I)), ...
                zeros(size(I))+min(real(obj.Ez),[],'all'), I, ...
                'EdgeColor','none')
            hold on;
            set(gca, 'FontSize', 20)
            view([0 90])
            xlabel('x', 'FontSize', 24)
            ylabel('y', 'FontSize', 24)
            set(gca, 'XTick', [min(obj.x) 0 max(obj.x)])
            set(gca, 'YTick', [min(obj.y) 0 max(obj.y)])
            set(gca, 'LineWidth', 2)
            grid off
            box on
            xlim([min(obj.x) max(obj.x)])
            ylim([min(obj.y) max(obj.y)])

            quiver3(obj.x(1:step_size:end),obj.y(1:step_size:end),...
                zeros(size(obj.y(1:step_size:end))),real(obj.Ex(1:step_size:end)),...
                real(obj.Ey(1:step_size:end)),real(obj.Ez(1:step_size:end)), quiver_varargin{:})
            pbaspect([2 sqrt(3) 1])
            view([0 90])
        end

        function plot_norm(obj)

            E_norm = obj.norm;
            I = ones(sqrt(length(obj.kx)), sqrt(length(obj.ky)))*NaN;
            for i = 1:length(obj.x)
                I((obj.kx - obj.x(i) < abs(obj.x(i)+1e-5)*1e-10 & ...
                    obj.kx - obj.x(i) > -abs(obj.x(i)+1e-5)*1e-10 & ...
                    obj.ky - obj.y(i) < abs(obj.y(i)+1e-5)*1e-10 & ...
                    obj.ky - obj.y(i) > -abs(obj.y(i)+1e-5)*1e-10)) ...
                    = E_norm(i);
            end

            figure;
            surf(reshape(obj.kx, size(I)), reshape(obj.ky, size(I)), I, ...
                'EdgeColor','none')
            hold on;
            view([0 90])
            xlabel('$x_1$', 'FontSize', 18, 'Interpreter', 'latex')
            ylabel('$x_2$', 'FontSize', 18, 'Interpreter', 'latex')
            set(gca, 'XTick', [min(obj.x) 0 max(obj.x)])
            set(gca, 'YTick', [min(obj.y) 0 max(obj.y)])
            set(gca, 'LineWidth', 2)
            set(gca, 'FontSize', 16)
            grid off
            box on
            xlim([min(obj.x) max(obj.x)])
            ylim([min(obj.y) max(obj.y)])
            pbaspect([2 sqrt(3) 1])
            view([0 90])
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

            I = ones(sqrt(length(obj.kx)), sqrt(length(obj.ky)))*NaN;
            for i = 1:length(obj.x)
                I((obj.kx - obj.x(i) < abs(obj.x(i)+1e-5)*1e-10 & ...
                    obj.kx - obj.x(i) > -abs(obj.x(i)+1e-5)*1e-10 & ...
                    obj.ky - obj.y(i) < abs(obj.y(i)+1e-5)*1e-10 & ...
                    obj.ky - obj.y(i) > -abs(obj.y(i)+1e-5)*1e-10)) ...
                    = vec(i);
            end

            figure;
            surf(reshape(obj.kx, size(I)), reshape(obj.ky, size(I)), ...
                zeros(size(I))+min(real(obj.Ez),[],'all'), ...
                real(I), 'EdgeColor','none')
            hold on;
            set(gca, 'FontSize', 18)
            view([0 90])
            xlabel('x', 'FontSize', 20)
            ylabel('y', 'FontSize', 20)
            set(gca, 'XTick', [min(obj.x) 0 max(obj.x)])
            set(gca, 'YTick', [min(obj.y) 0 max(obj.y)])
            set(gca, 'LineWidth', 2)
            grid off
            box on
            xlim([min(obj.x) max(obj.x)])
            ylim([min(obj.y) max(obj.y)])
            pbaspect([2 sqrt(3) 1])
            view([0 90])

        end

        function quiver(obj, varargin)

            if nargin == 2
                scale = varargin{1};
            else
                scale = 1; 
            end
            figure;

            quiver3(obj.x,obj.y,zeros(size(obj.y)),real(obj.Ex),real(obj.Ey),real(obj.Ez), scale, 'k')
            xlabel('$x_1$', 'FontSize', 18, 'Interpreter', 'latex')
            ylabel('$x_2$', 'FontSize', 18, 'Interpreter', 'latex')
            set(gca, 'XTick', [min(obj.x) 0 max(obj.x)])
            set(gca, 'YTick', [min(obj.y) 0 max(obj.y)])
            set(gca, 'LineWidth', 2)
            set(gca, 'FontSize', 16)
            grid off
            box on
            xlim([min(obj.x) max(obj.x)])
            ylim([min(obj.y) max(obj.y)])
            pbaspect([2 sqrt(3) 2])
            axis equal
        end

        % Apply a projector to the field and return the resultant field as
        % a new Field object.  That is, return: |psi'> = P|psi>
        function out = get_projection(obj, oprtr, oprtr_complement)

            sym_field = obj.symmetrise;
            
            % Check that the projectors are of the right size
            if size(oprtr,2) ~= length(sym_field.Ex(:))
                error(['get_projection(): First argument does not have the same number ' ...
                    'of columns as the number of elements in the Ex field.']);
            end
            if size(oprtr,2) ~= length(sym_field.Ez(:))
                error(['get_projection(): First argument does not have the same number ' ...
                    'of columns as the number of elements in the Ez field.']);
            end
            if size(oprtr_complement,2) ~= length(sym_field.Ey(:))
                error(['get_projection(): Second argument does not have the same number ' ...
                    'of columns as the number of elements in the Ey field.']);
            end

            % apply the projectors to the field components.
            % need to avoid overflow errors with using single precision for memory
            % saving.
            newEx = double(oprtr)*double(sym_field.Ex(:));
            newEx = reshape(newEx, [sym_field.size_x sym_field.size_y]);
            newEy = double(oprtr_complement)*double(sym_field.Ey(:));
            newEy = reshape(newEy, [sym_field.size_x sym_field.size_y]);
            newEz = double(oprtr)*double(sym_field.Ez(:));
            newEz = reshape(newEz, [sym_field.size_x sym_field.size_y]);

            % Desymmetrise the field

            % WrapToPi() makes the results easier to analyse.
            ang_ref = -wrapToPi(angle(obj.x+1i*obj.y));


            xsym = newEx.*cos(ang_ref) + newEy.*sin(ang_ref);

            ysym = -newEx.*sin(wrapToPi(ang_ref)) + newEy.*cos(wrapToPi(ang_ref));

            % Create Field object
            out = Field_hex(sym_field.x, sym_field.y, xsym, ysym, newEz, ...
                sym_field.kx, sym_field.ky);
        end

    end
end