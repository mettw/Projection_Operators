classdef ScalarFieldU3D < handle
    % ScalarField
    % Handles the orientation etc of scalar fields for use with projection
    % operators.

    properties
        x; % x coordinates
        y; % y coordinates
        z; % z coordinates
        fn; % scalar function
        size_x;
        size_y;
        size_z;
    end

    methods
        function obj = ScalarFieldU3D(x, y, z, fn, size_x, size_y, size_z)
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
            if iscolumn(fn) || isrow(fn)
                obj.fn = reshape(fn, [size_x size_y size_z]);
            else
                obj.fn = fn;
            end

            % correct the orientation of the array elements so that they
            % match Cartesian coordiantes.
            if obj.x(1,1,1) == obj.x(1,end,1)
                obj.x = rot90(obj.x);
                obj.y = rot90(obj.y);
                obj.z = rot90(obj.z);
                obj.fn = rot90(obj.fn);
            end
            if obj.y(1,1,1) < obj.y(end,1,1)
                obj.x = flipud(obj.x);
                obj.y = flipud(obj.y);
                obj.z = flipud(obj.z);
                obj.fn = flipud(obj.fn);
            end
            if obj.x(1,1,1) > obj.x(1,end,1)
                obj.x = fliplr(obj.x);
                obj.y = fliplr(obj.y);
                obj.z = fliplr(obj.z);
                obj.fn = fliplr(obj.fn);
            end
            if obj.z(1,1,1) < obj.z(1,1,end)
                obj.x = flip(obj.x,3);
                obj.y = flip(obj.y,3);
                obj.z = flip(obj.z,3);
                obj.fn = flip(obj.fn,3);
            end
        end

        % Shift the field to the left and down for cases where the highest
        % symmetry point is not at the center.  This assumes a particular
        % case and is not generally usefull.
        %
        % Not modified for 3D fields yet!
        function out = shifti(obj, varargin)

            switch nargin
                case 1
                    fn_shifted = [obj.fn(:,(obj.size_x+(obj.size_x-1)/2)/2+1:end) ...
                        obj.fn(:,1:(obj.size_x+(obj.size_x-1)/2)/2)];
                    fn_shifted = [fn_shifted(((obj.size_y-1)/2-1)/2+1:end,:); ...
                        fn_shifted(1:((obj.size_y-1)/2-1)/2,:)];
                otherwise
                    shift_x = varargin{1};
                    shift_y = varargin{2};
                    if shift_x < 0
                        shift_x = -shift_x;
                        fn_shifted = [obj.fn(:,end-shift_x+1:end) ...
                            obj.fn(:,1:1:end-shift_x)];
                    elseif shift_x > 0
                        fn_shifted = [obj.fn(:,shift_x+1:end) ...
                            obj.fn(:,1:shift_x)];
                    else
                        fn_shifted = obj.fn;
                    end
                    if shift_y < 0
                        shift_y = -shift_y;
                        fn_shifted = [fn_shifted(end-shift_y+1:end,:); ...
                            fn_shifted(1:end-shift_y,:)];
                    elseif shift_y > 0
                        fn_shifted = [fn_shifted(shift_y+1:end,:); ...
                            fn_shifted(1:shift_y,:)];
                    end
            end

            out = ScalarFieldU(obj.x,obj.y,fn_shifted,obj.size_x,obj.size_y);
        end

        % This function is used to determine if a scalar field has a particular
        % symmetry.
        %
        % This function returns
        %
        % <bra_x|oprtr|ket_x> 
        % -------------------
        %    <bra_x|ket_x>
        %

        function out = projection(obj,oprtr)

            if size(oprtr,2) ~= length(obj.fn(:))
                error(['projection(): First argument does not have the same number ' ...
                    'of columns as the number of elements in the Ex field.']);
            end

            % need to avoid overflow errors with using single precision for memory
            % saving.
            out = (double(obj.fn(:))'*double(oprtr)*double(obj.fn(:)))/ ...
                (double(obj.fn(:))'*double(obj.fn(:)));
        end

        function out = projectionU(obj,oprtr_basis)

            if size(oprtr_basis,1) ~= length(obj.fn(:))
                error(['projection(): First argument does not have the same number ' ...
                    'of rows as the number of elements in the Ex field.']);
            end

            % need to avoid overflow errors with using single precision for memory
            % saving.
            out = (double(obj.fn(:))'*double(oprtr_basis)*double(oprtr_basis')*double(obj.fn(:)))/ ...
                (double(obj.fn(:))'*double(obj.fn(:)));
        end

        function out = projection2(obj,oprtr, normoprtr)

            if size(oprtr,2) ~= length(obj.fn(:))
                error(['projection(): First argument does not have the same number ' ...
                    'of columns as the number of elements in the Ex field.']);
            end
            if size(normoprtr,2) ~= length(obj.fn(:))
                error(['projection(): First argument does not have the same number ' ...
                    'of columns as the number of elements in the Ex field.']);
            end

            % need to avoid overflow errors with using single precision for memory
            % saving.
            out = (double(obj.fn(:))'*double(oprtr)*double(obj.fn(:)))/ ...
                (double(obj.fn(:))'*double(normoprtr)*double(obj.fn(:)));
        end

        function plot(obj, vararg)
            switch nargin
                case 1
                    cb_str = [];
                case 2
                    cb_str = vararg{1};
                otherwise
                    error("Too many arguments.");
            end
            obj.plot_fn(real(obj.fn), "Real part", cb_str);
            if ~isreal(obj.fn)
                obj.plot_fn(imag(obj.fn), "Imaginary part", cb_str);
            end
        end

        function plot_fn(obj, real_fn, title_str, cb_str)
            figure;
            for i = 1:obj.size_z
                surf(obj.x(:,:,i),obj.y(:,:,i),obj.z(:,:,i), real_fn(:,:,i), 'EdgeColor','none')
                %surf(obj.x,obj.y,real_fn, real_fn, 'EdgeColor','none')
                hold on;
            end

            xlabel('$x_1$', 'FontSize', 18, 'Interpreter', 'latex')
            ylabel('$x_2$', 'FontSize', 18, 'Interpreter', 'latex')
            zlabel('$x_3$', 'FontSize', 18, 'Interpreter', 'latex')
            title(title_str, 'FontSize', 18, 'Interpreter', 'latex')
            set(gca, 'XTick', [obj.x(1,1,1) 0 obj.x(1,end,1)])
            set(gca, 'XTickLabels', {num2str(round_exponential(obj.x(1,1,1))), ...
                '0', num2str(round_exponential(obj.x(1,end,1)))})
            set(gca, 'YTick', [obj.y(end,1,1) 0 obj.y(1,1,1)])
            set(gca, 'YTickLabels', {num2str(round_exponential(obj.y(end,1,1))), ...
                '0', num2str(round_exponential(obj.y(1,1,1)))})
            set(gca, 'ZTick', [obj.z(1,1,end) 0 obj.z(1,1,1)])
            set(gca, 'ZTickLabels', {num2str(round_exponential(obj.z(1,1,end))), ...
                '0', num2str(round_exponential(obj.z(1,1,1)))})
            set(gca, 'LineWidth', 2)
            set(gca, 'FontSize', 16)
            grid off
            box on
            xlim([obj.x(1,1,1) obj.x(1,end,1)])
            ylim([obj.y(end,1,1) obj.y(1,1,1)])
            zlim([obj.z(1,1,end) obj.z(1,1,1)])
            cb = colorbar;
            ylabel(cb, cb_str, 'FontSize', 18, 'Interpreter', 'latex');

            pbaspect([1 1 1])
            
            function rounded_num = round_exponential(num)
                % Calculate the logarithm base 10 of the absolute value of the number
                log_value = log10(abs(num));
                
                % Extract the exponent part by taking the floor of the log value
                exp = floor(log_value);
                
                % Calculate the factor to round the number to the nearest integer
                factor = 10^exp;
                
                % Round the number to the nearest integer
                rounded_num = round(num / factor) * factor;
            end
        end
    end
end