classdef Field < handle
    % Field
    % Handles the orientation etc of vector fields for use with projection
    % operators.

    properties
        x;
        y;
        Ex;
        Ey;
        Ez;
        size_x;
        size_y;
        % The Hilbert space used to symmetrise a field
        F_x = [];
        F_y = [];
    end

    methods
        function obj = Field(x,y,Ex,Ey,Ez, size_x, size_y)
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
            if iscolumn(Ex) || isrow(Ex)
                obj.Ex = reshape(Ex, [size_x size_y]);
            else
                obj.Ex = Ex;
            end
            if iscolumn(Ey) || isrow(Ey)
                obj.Ey = reshape(Ey, [size_x size_y]);
            else
                obj.Ey = Ey;
            end
            if iscolumn(Ez) || isrow(Ez)
                obj.Ez = reshape(Ez, [size_x size_y]);
            else
                obj.Ez = Ez;
            end

            % correct the orientation of the array elements so that they
            % match Cartesian coordiantes.
            if obj.x(1,1) == obj.x(1,end)
                obj.x = rot90(obj.x);
                obj.y = rot90(obj.y);
                obj.Ex = rot90(obj.Ex);
                obj.Ey = rot90(obj.Ey);
                obj.Ez = rot90(obj.Ez);
            end
            if obj.y(1,1) < obj.y(end,1)
                obj.x = flipud(obj.x);
                obj.y = flipud(obj.y);
                obj.Ex = flipud(obj.Ex);
                obj.Ey = flipud(obj.Ey);
                obj.Ez = flipud(obj.Ez);
            end
            if obj.x(1,1) > obj.x(1,end)
                obj.x = fliplr(obj.x);
                obj.y = fliplr(obj.y);
                obj.Ex = fliplr(obj.Ex);
                obj.Ey = fliplr(obj.Ey);
                obj.Ez = fliplr(obj.Ez);
            end
        end

        % Shift the field to the left and down for cases where the highest
        % symmetry point is not at the center.  This assumes a particular
        % case and is not generally usefull.
        function out = shift(obj)
            Ex_shifted = [obj.Ex(:,(obj.size_x+(obj.size_x-1)/2)/2+1:end-1) ...
                obj.Ex(:,1:(obj.size_x+(obj.size_x-1)/2)/2+1)];
            Ey_shifted = [obj.Ey(:,(obj.size_x+(obj.size_x-1)/2)/2+1:end-1) ...
                obj.Ey(:,1:(obj.size_x+(obj.size_x-1)/2)/2+1)];
            Ez_shifted = [obj.Ez(:,(obj.size_x+(obj.size_x-1)/2)/2+1:end-1) ...
                obj.Ez(:,1:(obj.size_x+(obj.size_x-1)/2)/2+1)];
            Ex_shifted = [Ex_shifted(((obj.size_y-1)/2-1)/2+1:end-1,:); ...
                Ex_shifted(1:((obj.size_y-1)/2-1)/2+1,:)];
            Ey_shifted = [Ey_shifted(((obj.size_y-1)/2-1)/2+1:end-1,:); ...
                Ey_shifted(1:((obj.size_y-1)/2-1)/2+1,:)];
            Ez_shifted = [Ez_shifted(((obj.size_y-1)/2-1)/2+1:end-1,:); ...
                Ez_shifted(1:((obj.size_y-1)/2-1)/2+1,:)];

            out = Field(obj.x,obj.y,Ex_shifted,Ey_shifted,Ez_shifted,obj.size_x,obj.size_y);
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

            [X,Y] = meshgrid(-(obj.size_x-1)/2:(obj.size_x-1)/2, ...
                (obj.size_y-1)/2:-1:-(obj.size_y-1)/2);

            % WrapToPi() makes the results easier to analyse.
            ang_ref = wrapToPi(angle(X+1i*Y));


            xsym = obj.Ex.*cos(ang_ref) + obj.Ey.*sin(ang_ref);

            ysym = -obj.Ex.*sin(wrapToPi(ang_ref)) + obj.Ey.*cos(wrapToPi(ang_ref));

            out = Field(obj.x, obj.y, xsym, ysym, obj.Ez, obj.size_x, obj.size_y);
            out.F_x = X;
            out.F_y = Y;
        end

        function out = norm(obj)
            out = sqrt(obj.Ex.*conj(obj.Ex)+obj.Ey.*conj(obj.Ey)+obj.Ez.*conj(obj.Ez));
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
            if size(oprtr,2) ~= length(obj.Ez(:))
                error(['norm_vec(): First argument does not have the same number ' ...
                    'of columns as the number of elements in the Ez field.']);
            end
            if size(oprtr_complement,2) ~= length(obj.Ey(:))
                error(['norm_vec(): Second argument does not have the same number ' ...
                    'of columns as the number of elements in the Ey field.']);
            end

            % need to avoid overflow errors with using single precision for memory
            % saving.
            out = (double(obj.Ex(:))'*double(oprtr)*double(obj.Ex(:))+ ...
                double(obj.Ey(:))'*double(oprtr_complement)*double(obj.Ey(:))+ ...
                double(obj.Ez(:))'*double(oprtr)*double(obj.Ez(:)))/ ...
                (double(obj.Ex(:))'*double(obj.Ex(:))+double(obj.Ey(:))'*double(obj.Ey(:))+ ...
                +double(obj.Ez(:))'*double(obj.Ez(:)));
        end


        function plot(obj)
            figure;
            surf(obj.x,obj.y,zeros(size(obj.x))+min(real(obj.Ez),[],'all'), obj.norm, ...
                'EdgeColor','none')
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

            quiver3(obj.x,obj.y,zeros(size(obj.y)),real(obj.Ex),real(obj.Ey),real(obj.Ez), 'k')
            pbaspect([1 1 1])
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
    end
end