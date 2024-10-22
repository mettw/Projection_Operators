classdef LinSym3D
    % LinSym - Implementation of linear symmetry caluclations and plotting

    properties
        U;  % The bases as a single matrix ie it is the direct sum of the E 
            % basis for each of the projectors passed to the constructor.
            % U = Exy (+) Eab (+) Elr

        % Size of each basis
        size_xy_E_11g;
        size_xy_E_22g;
        size_ab_E_11g;
        size_ab_E_22g;
        size_lr_E_11;
        size_lr_E_22;
    end

    methods
        function obj = LinSym3D(Pxy, Pab, Plr)
            obj.U = [Pxy.basis("E_11g") Pxy.basis("E_22g") ...
                Pab.basis("E_11g") Pab.basis("E_22g") ...
                Plr.basis("E_11") Plr.basis("E_22")];
            obj.U = sparse(obj.U);

            obj.size_xy_E_11g = size(Pxy.basis("E_11g"),2);
            obj.size_xy_E_22g = size(Pxy.basis("E_22g"),2);
            obj.size_ab_E_11g = size(Pab.basis("E_11g"),2);
            obj.size_ab_E_22g = size(Pab.basis("E_22g"),2);
            obj.size_lr_E_11 = size(Plr.basis("E_11"),2);
            obj.size_lr_E_22 = size(Plr.basis("E_22"),2);
        end

        function sigma = get_linear_symmetry(obj, vec)
            sym_vec = obj.U'*vec;

            base = 0;
            sigma_0 = sym_vec(base+1:base+obj.size_xy_E_11g)'*sym_vec(base+1:base+obj.size_xy_E_11g);
            sigma_1 = sigma_0;
            base = base+obj.size_xy_E_11g;
            sigma_0 = sigma_0 + sym_vec(base+1:base+obj.size_xy_E_22g)'*sym_vec(base+1:base+obj.size_xy_E_22g);
            sigma_1 = sigma_1 - sym_vec(base+1:base+obj.size_xy_E_22g)'*sym_vec(base+1:base+obj.size_xy_E_22g);

            base = base + obj.size_xy_E_22g;
            sigma_2 = sym_vec(base+1:base+obj.size_ab_E_11g)'*sym_vec(base+1:base+obj.size_ab_E_11g);
            base = base+obj.size_ab_E_11g;
            sigma_2 = sigma_2 - sym_vec(base+1:base+obj.size_ab_E_22g)'*sym_vec(base+1:base+obj.size_ab_E_22g);

            base = base+obj.size_ab_E_22g;
            sigma_3 = -sym_vec(base+1:base+obj.size_lr_E_11)'*sym_vec(base+1:base+obj.size_lr_E_11);
            base = base+obj.size_lr_E_11;
            sigma_3 = sigma_3 + sym_vec(base+1:base+obj.size_lr_E_22)'*sym_vec(base+1:base+obj.size_lr_E_22);

            sigma = [sigma_0 sigma_1 sigma_2 sigma_3]/(vec'*vec);
        end

        function sigma = get_linear_symmetry_field(obj, field)
            field_sym = field.symmetrise;
            vec_x = field_sym.Ex(:);
            vec_y = field_sym.Ey(:);
            vec_z = field_sym.Ez(:);
            sym_vec_x = obj.U'*vec_x;
            sym_vec_y = obj.U'*vec_y;
            sym_vec_z = obj.U'*vec_z;

            base = 0;
            sigma_0_x = sym_vec_x(base+1:base+obj.size_xy_E_11g)'*sym_vec_x(base+1:base+obj.size_xy_E_11g);
            sigma_1_x = sigma_0_x;

            sigma_0_z = sym_vec_z(base+1:base+obj.size_xy_E_11g)'*sym_vec_z(base+1:base+obj.size_xy_E_11g);
            sigma_1_z = sigma_0_z;
            sigma_0_y = sym_vec_y(base+1:base+obj.size_xy_E_11g)'*sym_vec_y(base+1:base+obj.size_xy_E_11g);
            sigma_1_y = sigma_0_y;
            base = base+obj.size_xy_E_11g;
            sigma_0_x = sigma_0_x + sym_vec_x(base+1:base+obj.size_xy_E_22g)'*sym_vec_x(base+1:base+obj.size_xy_E_22g);
            sigma_1_x = sigma_1_x - sym_vec_x(base+1:base+obj.size_xy_E_22g)'*sym_vec_x(base+1:base+obj.size_xy_E_22g);

            sigma_0_z = sigma_0_z + sym_vec_z(base+1:base+obj.size_xy_E_22g)'*sym_vec_z(base+1:base+obj.size_xy_E_22g);
            sigma_1_z = sigma_1_z - sym_vec_z(base+1:base+obj.size_xy_E_22g)'*sym_vec_z(base+1:base+obj.size_xy_E_22g);

            sigma_0_y = sigma_0_y + sym_vec_y(base+1:base+obj.size_xy_E_22g)'*sym_vec_y(base+1:base+obj.size_xy_E_22g);
            sigma_1_y = sigma_1_y - sym_vec_y(base+1:base+obj.size_xy_E_22g)'*sym_vec_y(base+1:base+obj.size_xy_E_22g);

            base = base + obj.size_xy_E_22g;
            sigma_2_x = sym_vec_x(base+1:base+obj.size_ab_E_11g)'*sym_vec_x(base+1:base+obj.size_ab_E_11g);
            sigma_2_z = sym_vec_z(base+1:base+obj.size_ab_E_11g)'*sym_vec_z(base+1:base+obj.size_ab_E_11g);
            sigma_2_y = sym_vec_y(base+1:base+obj.size_ab_E_11g)'*sym_vec_y(base+1:base+obj.size_ab_E_11g);
            base = base+obj.size_ab_E_11g;
            sigma_2_x = sigma_2_x - sym_vec_x(base+1:base+obj.size_ab_E_22g)'*sym_vec_x(base+1:base+obj.size_ab_E_22g);
            sigma_2_z = sigma_2_z - sym_vec_z(base+1:base+obj.size_ab_E_22g)'*sym_vec_z(base+1:base+obj.size_ab_E_22g);
            sigma_2_y = sigma_2_y - sym_vec_y(base+1:base+obj.size_ab_E_22g)'*sym_vec_y(base+1:base+obj.size_ab_E_22g);

            base = base+obj.size_ab_E_22g;
            sigma_3_x = -sym_vec_x(base:base+obj.size_lr_E_11)'*sym_vec_x(base:base+obj.size_lr_E_11);
            sigma_3_z = -sym_vec_z(base:base+obj.size_lr_E_11)'*sym_vec_z(base:base+obj.size_lr_E_11);
            sigma_3_y = -sym_vec_y(base:base+obj.size_lr_E_11)'*sym_vec_y(base:base+obj.size_lr_E_11);
            base = base+obj.size_lr_E_11;
            sigma_3_x = sigma_3_x + sym_vec_x(base+1:base+obj.size_lr_E_22)'*sym_vec_x(base+1:base+obj.size_lr_E_22);
            sigma_3_z = sigma_3_z + sym_vec_z(base+1:base+obj.size_lr_E_22)'*sym_vec_z(base+1:base+obj.size_lr_E_22);
            sigma_3_y = sigma_3_y + sym_vec_y(base+1:base+obj.size_lr_E_22)'*sym_vec_y(base+1:base+obj.size_lr_E_22);

            sigma = [(sigma_0_x+sigma_0_y+sigma_0_z) (sigma_1_x+sigma_1_y+sigma_1_z) ...
                (sigma_2_x+sigma_2_y+sigma_2_z) (sigma_3_x+sigma_3_y+sigma_3_z)]/ ...
                (vec_x'*vec_x+vec_y'*vec_y+vec_z'*vec_z);
        end

        function sigma = get_linear_symmetry_orig(obj, vec)
            base = 0;
            A = obj.U(:,base+1:base+obj.size_xy_E_11g);
            sym_vec = (vec'*A)*(A'*vec);
            sigma_0 = sym_vec;
            sigma_1 = sym_vec;
            base = base+obj.size_xy_E_11g;
            A = obj.U(:,base+1:base+obj.size_xy_E_22g);
            sym_vec = (vec'*A)*(A'*vec);
            sigma_0 = sigma_0 + sym_vec;
            sigma_1 = sigma_1 - sym_vec;

            base = base + obj.size_xy_E_22g;
            A = obj.U(:,base+1:base+obj.size_ab_E_11g);
            sigma_2 = (vec'*A)*(A'*vec);
            base = base+obj.size_ab_E_11g;
            A = obj.U(:,base+1:base+obj.size_ab_E_22g);
            sigma_2 = sigma_2 - (vec'*A)*(A'*vec);

            base = base+obj.size_ab_E_22g;
            A = obj.U(:,base:base+obj.size_lr_E_11);
            sigma_3 = -(vec'*A)*(A'*vec);
            base = base+obj.size_lr_E_11;
            A = obj.U(:,base+1:base+obj.size_lr_E_22);
            sigma_3 = sigma_3 + (vec'*A)*(A'*vec);

            sigma = [sigma_0 sigma_1 sigma_2 sigma_3]/(vec'*vec);
        end
        
        function cb = plot(obj, vec)
            sigma = obj.get_linear_symmetry(vec);
            S0 = sigma(1);
            S1 = sigma(2);
            S2 = sigma(3);
            S3 = sigma(4);
            % Calculate the x and y coordinates of the points to be plotted.
            r = 1./(1-log10(S0));
            phi = atan(sqrt(abs(S2./S1)));
            phi(isnan(phi)) = 0; % E==0 exactly gives NaN
            theta = atan(sqrt(abs(S3)./(abs(S1)+abs(S2))));
            theta(isnan(theta)) = 0; % E==0 exactly gives NaN
        
            x = r.*sign(S1).*cos(phi).*cos(theta);
            y = r.*sign(S2).*sin(phi).*cos(theta);
            z = r.*sign(S3).*sin(theta);
        
            %Plot the grid lines.
            circle(1/(1-log10(0.5)), 0.75);
            hold on;
            for i=1:20
                circle(1/(1-log10(10^(-i))), 1.5);
            end
        
            line([1,-1],[0,0], 'Color', [1 1 1]*0.8, 'LineWidth', 1.5);
            line([0,0],[1,-1], 'Color', [1 1 1]*0.8, 'LineWidth', 1.5);
            line(cosd(30)*[1,-1],sind(30)*[1,-1], 'Color', [1 1 1]*0.8, 'LineWidth', 1.5);
            line(cosd(30)*[-1,1],sind(30)*[1,-1], 'Color', [1 1 1]*0.8, 'LineWidth', 1.5);
            line(cosd(60)*[1,-1],sind(60)*[1,-1], 'Color', [1 1 1]*0.8, 'LineWidth', 1.5);
            line(cosd(60)*[-1,1],sind(60)*[1,-1], 'Color', [1 1 1]*0.8, 'LineWidth', 1.5);
            text(0,-0.95,0,'1','Color', [1 1 1]*0.3)
            text(0,-1/(1-log10(0.5)),0,'5 E-1','Color', [1 1 1]*0.3)
            text(0,-1/(1-log10(0.1)),0,'1 E-1','Color', [1 1 1]*0.3)
            text(0,-1/(1-log10(0.01)),0,'1 E-2','Color', [1 1 1]*0.3)
            text(0,-1/(1-log10(0.001)),0,'1 E-3','Color', [1 1 1]*0.3)
            text(0,0,0,'0','Color', [1 1 1]*0.3)
            % Plot the S axes.
            plotEquator;
    
            scatter3(x, y, z, 5, 'b', 'filled')
            cb = [];
        
            % sundry.
            ax_len = 1.1;
            text(ax_len-0.05,0,0, 'x', 'FontSize',18);
            text(-0.05,ax_len,0, 'a', 'FontSize',18);
            text(-ax_len-0.05,0,0, 'y', 'FontSize',18);
            text(-0.05,-ax_len,0, 'b', 'FontSize',18);
            text(0,0,ax_len+0.05, 'x+iy', 'FontSize',18);
            text(0,0,-0.05-ax_len, 'x-iy', 'FontSize',18);
        
            axis off;
            pbaspect([1 1 1]);
            set(gca, 'FontSize', 14)
            view([-15 15])
        
            %%%%%%%%%%%%%%%%%%
            %
            % EOF
            %
            %%%%%%%%%%%%%%%%%%
        
            function plotEquator
                center = [0,0,0];
                normal = [0,0,1];
                radius = 1;
                theta_eq=0:0.01:2*pi;
                v=null(normal);
                points=repmat(center',1,size(theta_eq,2))+radius*(v(:,1)*cos(theta_eq)+v(:,2)*sin(theta_eq));
                plot3(points(1,:),points(2,:),points(3,:),'k-', 'LineWidth', 2);
                plot3(points(2,:),points(3,:),points(1,:),'k-', 'LineWidth', 2);
                plot3(points(3,:),points(1,:),points(2,:),'k-', 'LineWidth', 2);
                scatter3(1,0,0, 40,'filled', 'MarkerFaceColor', 'k')
                scatter3(-1,0,0, 40,'filled', 'MarkerFaceColor', 'k')
                scatter3(0,1,0, 40,'filled', 'MarkerFaceColor', 'k')
                scatter3(0,-1,0, 40,'filled', 'MarkerFaceColor', 'k')
                scatter3(0,0,1, 40,'filled', 'MarkerFaceColor', 'k')
                scatter3(0,0,-1, 40,'filled', 'MarkerFaceColor', 'k')
            end
        
            function circle(r, lwidth)
                %x and y are the coordinates of the center of the circle
                %r is the radius of the circle
                %0.01 is the angle step, bigger values will draw the circle faster but
                %you might notice imperfections (not very smooth)
                ang=0:0.01:2*pi;
                xp=r*cos(ang);
                yp=r*sin(ang);
                plot(xp,yp, 'Color', [1 1 1]*0.8, 'LineWidth', lwidth);
            end
        end

    end
end