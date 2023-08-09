classdef ScalarField < handle
    % ScalarField
    % Handles the orientation etc of scalar fields for use with projection
    % operators.

    properties
        x; % x coordinates
        y; % y coordinates
        fn; % scalar function
        size_x;
        size_y;
        % The Hilbert space used to symmetrise a field
        F_x = [];
        F_y = [];
    end

    methods
        function obj = ScalarField(x, y, fn, size_x, size_y)
            obj.size_x = size_x;
            obj.size_y = size_y;

            if iscolumn(x) || isrow(x)
                obj.x = reshape(x, [size_x size_y]);
            else
                obj.x = x;
            end
            if iscolumn(y) || isrow(y)
                obj.y = reshape(y, [size_x size_y]);
            else
                obj.y = y;
            end
            if iscolumn(fn) || isrow(fn)
                obj.fn = reshape(fn, [size_x size_y]);
            else
                obj.fn = fn;
            end

            % correct the orientation of the array elements so that they
            % match Cartesian coordiantes.
            if obj.x(1,1) == obj.x(1,end)
                obj.x = rot90(obj.x);
                obj.y = rot90(obj.y);
                obj.fn = rot90(obj.fn);
            end
            if obj.y(1,1) < obj.y(end,1)
                obj.x = flipud(obj.x);
                obj.y = flipud(obj.y);
                obj.fn = flipud(obj.fn);
            end
            if obj.x(1,1) > obj.x(1,end)
                obj.x = fliplr(obj.x);
                obj.y = fliplr(obj.y);
                obj.fn = fliplr(obj.fn);
            end
        end

        % Shift the field to the left and down for cases where the highest
        % symmetry point is not at the center.  This assumes a particular
        % case and is not generally usefull.
        function out = shift(obj, varargin)

            switch nargin
                case 1
                    fn_shifted = [obj.fn(:,(obj.size_x+(obj.size_x-1)/2)/2+1:end-1) ...
                        obj.fn(:,1:(obj.size_x+(obj.size_x-1)/2)/2+1)];
                    fn_shifted = [fn_shifted(((obj.size_y-1)/2-1)/2+1:end-1,:); ...
                        fn_shifted(1:((obj.size_y-1)/2-1)/2+1,:)];
                otherwise
                    shift_x = varargin{1};
                    shift_y = varargin{2};
                    if shift_x < 0
                        shift_x = -shift_x;
                        fn_shifted = [obj.fn(:,end-shift_x+1:end-1) ...
                            obj.fn(:,1:1:end-shift_x+1)];
                    elseif shift_x > 0
                        fn_shifted = [obj.fn(:,shift_x:end-1) ...
                            obj.fn(:,1:shift_x)];
                    end
                    if shift_y < 0
                        shift_y = -shift_y;
                        fn_shifted = [fn_shifted(end-shift_y+1:end-1,:); ...
                            fn_shifted(1:end-shift_y+1,:)];
                    elseif shift_y > 0
                        fn_shifted = [fn_shifted(shift_y:end-1,:); ...
                            fn_shifted(1:shift_y,:)];
                    end
            end

            out = ScalarField(obj.x,obj.y,fn_shifted,obj.size_x,obj.size_y);
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

        function plot(obj)
            figure;
            %surf(obj.x,obj.y,zeros(size(obj.x)), obj.fn, 'EdgeColor','none')
            surf(obj.x,obj.y,obj.fn, obj.fn, 'EdgeColor','none')
            hold on;
            view([0 90])
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
            colorbar;

            pbaspect([1 1 1])
            view([0 90])
        end
    end
end