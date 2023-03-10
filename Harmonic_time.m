classdef Harmonic_time < Harmonic
    % Harmonic_time - example of subclass of Harmonic.
    %
    % If we take a C4v square lattice for example, the four wavevectors at
    % the lowest order are all related by the pi/2 rotation symmetry (as
    % well as the mirror symmetries).  That is, there are four coupled
    % waves in a waveguide layer travelling in the +/-x and +/-y
    % directions.  By symmetry they all have the same magnitude, but might
    % be in or out of phase with each other as determined by the characters
    % of the different irreducible representations.
    %
    % If we reduce the symmetry of the MS on waveguide layer from C4v to
    % C2v then the waves travelling in the +/-x directions are related by
    % symmetry, as are those in the +/-y directions, but there is no longer
    % any symmetry relation between waves travelling in the x and y
    % directions since there is only a pi rotation symmetry.  This opens up
    % the possibility of there being an arbitrary phase difference between
    % the x and y travelling waves.
    %
    % This is reflected in the basis for C2v with the F1 basis: The two
    % dimentional C4v E mode reduces to C2v B1+B2 and the B1 and A1 modes
    % both reduce to the C2v A1 mode, which is two dimentional.
    % Examination of the bases reveals that the two dimentions of A1 are
    % exclusively x for the one and exclusively y for the other.  Likewise
    % B1 is exclusively x and B2 is exclusively y.
    %
    % These phase differences can ofcourse introduce elliptical or circular
    % polarisation and so it is usefull to study the evolution of the field
    % with time.
    %
    % In this class the phase differences are introduced in the
    % `mode_coeffs' parameter to the constructor and some new plotting
    % functions are introduced to look at how the polarisation changes with
    % different mode coefficients.

    methods
        function hObj = Harmonic_time(Hilbert_space, mode_basis, mode_coeffs)
            hObj = hObj@Harmonic(Hilbert_space, mode_basis, mode_coeffs);
        end


        % plot_no_norm() method does not normalise the field vectors and 
        % uses a constant colour scaling for the norm of the field.
        function plot_no_norm(hObj)
            figure;

            surf(hObj.X,hObj.Y,zeros(size(hObj.X)),hObj.vec_norm, 'EdgeColor','none')
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
                real(hObj.vec_y20),real(hObj.vec_z20), 'off', 'k')
            pbaspect([1 1 1])
            if lower(hObj.vec_type) == "n"
                view([-15 30])
            else
                view([0 90])
            end
        end

        % plot_time() - plot the evolution of the field in time.
        % t - vector of phase angles in radians to use as the time 
        % parameter.
        function plot_time(hObj, time_vals)
            figure;
            
            for t=time_vals
                c = get_color_triplet(jet, t, time_vals);
                quiver3(hObj.X20,hObj.Y20,ones(size(hObj.X20))*t,...
                    real(hObj.vec_x20*exp(1i*t)), real(hObj.vec_y20*exp(1i*t)),...
                    real(hObj.vec_z20*exp(1i*t)), 'off', 'Color', c)
                hold on
            end
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

            pbaspect([1 1 1])

            function c = get_color_triplet(cmap, t, time_vals)
                c = time_vals - min(time_vals);
                c = (t-min(time_vals))/max(c);
                c = round(c*size(cmap,1));
                if c==0
                    c=1;
                end
                c = cmap(c,:);
            end
        end


        function plot_polarisation(hObj, time_vals, space_scale)
            figure;
            
            Xv = real(hObj.vec_x20(:)*exp(1i*time_vals));
            Yv = real(hObj.vec_y20(:)*exp(1i*time_vals));
            %for X=1:length(hObj.X20(:))
                x = hObj.X20(:)*ones(1,length(time_vals))*space_scale;
                x = x+real(hObj.vec_x20(:)*exp(1i*time_vals));
                y = hObj.Y20(:)*ones(1,length(time_vals))*space_scale;
                y = y+real(hObj.vec_y20(:)*exp(1i*time_vals));
             %   for t=1:length(time_vals)
              %      x(t) = x(t)+real(hObj.vec_x20(X)*exp(1i*time_vals(t)));
               %     y(t) = y(t)+real(hObj.vec_y20(X)*exp(1i*time_vals(t)));
                %end
                plot(x.',y.');
                hold on
            %end
            %{
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

            pbaspect([1 1 1])
            %}
        end
    end
end