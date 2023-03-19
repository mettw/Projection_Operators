classdef Field_shift < Field
    % FIELD_SHIFT Shift field to account for symmetry
    % 
    % Sometimes the center of symmetry does not correspond to the center of
    % the unit cell, so shift the Field data so that they correspond.  This
    % will allow us to apply the proper point group projectors to the
    % field.
    
    methods
        function obj = Field_shift(x,y,Ex,Ey,Ez, size_x, size_y)
            obj = obj@Field(x,y,Ex,Ey,Ez, size_x, size_y);
        end
        
        % Shift the field to the left and down for cases where the highest
        % symmetry point is not at the center.  This assumes a particular
        % case and is not generally usefull.
        function out = shift(obj, varargin)

            switch nargin
                case 1
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
                otherwise
                    shift_x = varargin{1};
                    shift_y = varargin{2};
                    if shift_x < 0
                        shift_x = -shift_x;
                        Ex_shifted = [obj.Ex(:,end-shift_x+1:end-1) ...
                            obj.Ex(:,1:end-shift_x+1)];
                        Ey_shifted = [obj.Ey(:,end-shift_x+1:end-1) ...
                            obj.Ey(:,1:end-shift_x+1)];
                        Ez_shifted = [obj.Ez(:,end-shift_x+1:end-1) ...
                            obj.Ez(:,1:end-shift_x+1)];
                    else
                        Ex_shifted = [obj.Ex(:,shift_x+1:end-1) ...
                            obj.Ex(:,1:shift_x+1)];
                        Ey_shifted = [obj.Ey(:,shift_x+1:end-1) ...
                            obj.Ey(:,1:shift_x+1)];
                        Ez_shifted = [obj.Ez(:,shift_x+1:end-1) ...
                            obj.Ez(:,1:shift_x+1)];
                    end
                    if shift_y < 0
                        shift_y = -shift_y;
                        Ex_shifted = [Ex_shifted(end-shift_y+1:end-1,:); ...
                            Ex_shifted(1:end-shift_y+1,:)];
                        Ey_shifted = [Ey_shifted(end-shift_y+1:end-1,:); ...
                            Ey_shifted(1:end-shift_y+1,:)];
                        Ez_shifted = [Ez_shifted(end-shift_y+1:end-1,:); ...
                            Ez_shifted(1:end-shift_y+1,:)];
                    else
                        Ex_shifted = [Ex_shifted(shift_y+1:end-1,:); ...
                            Ex_shifted(1:shift_y+1,:)];
                        Ey_shifted = [Ey_shifted(shift_y+1:end-1,:); ...
                            Ey_shifted(1:shift_y+1,:)];
                        Ez_shifted = [Ez_shifted(shift_y+1:end-1,:); ...
                            Ez_shifted(1:shift_y+1,:)];
                    end
                    %Ex_shifted = [Ex_shifted(((obj.size_y-1)/2-1)/2+1:end-1,:); ...
                    %    Ex_shifted(1:((obj.size_y-1)/2-1)/2+1,:)];
                    %Ey_shifted = [Ey_shifted(((obj.size_y-1)/2-1)/2+1:end-1,:); ...
                    %    Ey_shifted(1:((obj.size_y-1)/2-1)/2+1,:)];
                    %Ez_shifted = [Ez_shifted(((obj.size_y-1)/2-1)/2+1:end-1,:); ...
                    %    Ez_shifted(1:((obj.size_y-1)/2-1)/2+1,:)];
            end


            out = Field(obj.x,obj.y,Ex_shifted,Ey_shifted,Ez_shifted,obj.size_x,obj.size_y);
        end

    end
end

