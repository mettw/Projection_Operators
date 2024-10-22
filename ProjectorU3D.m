classdef ProjectorU3D < dynamicprops & matlab.mixin.CustomDisplay
    % Projector - Create a projection operator
    %
    % Creates a projector for the specified basis and point group.
    %
    % REQUIRED:
    %
    % if you want to create a C6v projector then you will need the gso 
    % package from the MATLAB file exchange 
    % https://au.mathworks.com/matlabcentral/fileexchange/55881-gram-schmidt-orthogonalization?s_tid=srchtitle
    %
    % USAGE:
    %
    % The first parameter is the basis where each basis vector |K> is represented
    % as [m,n] where K=m*b1+n*b2.  The basis is a cell array of these 1x2 arrays.
    % The second parameter is a list of strings of the point groups that you want
    % projectors for.  eg.
    %
    % > P=Projector({[1; 0] [0; 1] [-1; 0] [0; -1]}, "C4v")
    %
    % For C6v the basis is
    %
    % > P=Projector({[1; 0] [1/2; sqrt(3)/2] [-1/2; sqrt(3)/2] [-1; 0] [-1/2; -sqrt(3)/2] [1/2; -sqrt(3)/2]}, "C6v")
    %
    % Alternatively you can name a standard basis such as
    %
    % > P=Projector("F1", "C4v")
    %
    % or a set of standard bases (Note the use of square brackets in this
    % case)
    %
    % > P=Projector(["F0" "F1" "F2"], "C4v")
    %
    % and for C6v 
    %
    % > P=Projector("F1_hex", "C6v")
    %
    % Note that because of the long run time for large bases the
    % constructor does not automatically build any of the projectors,
    % instead you need to build them explicitly as
    %
    % > P.add_irrep("A1")
    %
    % or you can build all of them by telling it to build the change of
    % coordinates matrix
    %
    % > P.add_irrep("U")
    %
    % Get the A1 projector
    %
    % > P.A1
    %
    %   ans =
    %
    %       0.2500    0.2500    0.2500    0.2500
    %       0.2500    0.2500    0.2500    0.2500
    %       0.2500    0.2500    0.2500    0.2500
    %       0.2500    0.2500    0.2500    0.2500
    %
    % Get the change of coordinates matrix
    %
    % > P.U
    %
    %   ans =
    %
    %       0.5000    0.5000    0.7071         0
    %       0.5000   -0.5000         0    0.7071
    %       0.5000    0.5000   -0.7071         0
    %       0.5000   -0.5000         0   -0.7071
    %
    % Get the basis for the A1 irreducible representation
    %
    % > P.basis("A1")
    %
    %   ans =
    %
    %       0.5000
    %       0.5000
    %       0.5000
    %       0.5000
    %
    % TODO:
    %
    % * Add more point groups as needed
    

    % The properties and methods are all public so that you can easily 
    % create subclasses of this class.
    properties
        % Fourier space
        % These are the (p,q) values for the vectors
        % va{k} = p\va{b}_1+q\va{b}_2
        F;
        F_len;
        % reciprocal space basis vectors (square by default)
        b_vecs = [1 0 0;0 1 0; 0 0 1];

        projs; % Names of the projectors to create
        ops; % Names of the symmetry operations
        U = []; % Change of  coordinates matrix
        U_calculated = false;

        % Point group operators
        E_1 = [1 0 0;0 1 0; 0 0 1];
        C4_1 = [0 -1 0;1 0 0; 0 0 1];
        sigma_v_x = [1 0 0;0 -1 0; 0 0 1];
        i_op = [-1 0 0; 0 -1 0; 0 0 -1];

        C4_3 = [];
        C2_1 = [];
        sigma_v_y = [];
        sigma_d_d = [];
        sigma_d_a = [];
        S4_1 = [];
        S4_3 = [];
        sigma_h = [];
        C2p_x = [];
        C2p_y = [];
        C2pp_d = [];
        C2pp_a = [];

        point_group;
    end
    
    methods
        %
        % SETUP functions
        %
        
        function hObj = ProjectorU3D(Fourier_space, point_group, vararg)
            
            hObj.create_operator_matrices;
            
            if nargin == 3
                hObj.include_y_projs = vararg{1};
            end
            % if Fourier_space is a name of a standard F_i space
            if isstring(Fourier_space)
                for F_i=Fourier_space
                    switch lower(F_i)
                        case "g0"
                            hObj.F = [hObj.F; {[0; 0; 0]}];
                        case "g1"
                            hObj.F = [hObj.F; {[1; 0; 0]; [0; 1; 0]; [0; 0; 1]; ...
                                [-1; 0; 0]; [0; -1; 0]; [0; 0; -1]}];
                        case "g2"
                            hObj.F = [hObj.F; {[1; 0; 1]; [0; 1; 1]; [-1; 0; 1]; [0; -1; 1];...
                                [1; 1; 0]; [-1; 1; 0]; [-1; -1; 0]; [1; -1; 0];...
                                [1; 0; -1]; [0; 1; -1]; [-1; 0; -1];[0; -1; -1]}];
                   
                        case "g3"
                            hObj.F = [hObj.F; {[2; 0]; [0; 2]; [-2; 0]; [0; -2]}];
                        case "g4"
                            hObj.F = [hObj.F; {[2; 1]; [1; 2]; [-1; 2]; [-2; 1]; ...
                                [-2; -1]; [-1; -2]; [1; -2]; [2; -1]}];
                        case "g5"
                            hObj.F = [hObj.F; {[2; 2]; [-2; 2]; [-2; -2]; [2; -2]}];
                        case "g6"
                            hObj.F = [hObj.F; {[3; 0]; [0; 3]; [-3; 0]; [0; -3]}];
                        case "g7"
                            hObj.F = [hObj.F; {[3; 1]; [1; 3]; [-1; 3]; [-3; 1]; ...
                                [-3; -1]; [-1; -3]; [1; -3]; [3; -1]}];
                        case "g8"
                            hObj.F = [hObj.F; {[3; 2]; [2; 3]; [-2; 3]; [-3; 2]; ...
                                [-3; -2]; [-2; -3]; [2; -3]; [3; -2]}];
                        case "g9"
                            hObj.F = [hObj.F; {[3; 3]; [-3; 3]; [-3; -3]; [3; -3]}];
                        case "g0_hex"
                            hObj.b_vecs = [1/2 1; sqrt(3)/2 0];
                            hObj.F = [hObj.F; {[0; 0]}];
                        case "g1_hex"
                            hObj.b_vecs = [1/2 1; sqrt(3)/2 0];
                            hObj.F = [hObj.F; {[1;0]; [0; 1]; [-1; 1]; [-1; 0]; [0; -1]; [1; -1]}];
                                %{[1; 0]; [1/2; sqrt(3)/2]; [-1/2; sqrt(3)/2]; ...
                                %[-1; 0]; [-1/2; -sqrt(3)/2]; [1/2; -sqrt(3)/2]}];
                        case "g1p_hex"
                            hObj.b_vecs = [sqrt(3)/2 sqrt(3)/2; 1/2 -1/2]*sqrt(3)/2;
                            hObj.F = [hObj.F; {[1;0]; [0; 1]; [-1; 1]; [-1; 0]; [0; -1]; [1; -1]}];
                        case "g2_hex"
                            hObj.b_vecs = [1/2 1; sqrt(3)/2 0];
                            hObj.F = [hObj.F; {[1; 1]; [1; 2]; [-2; 1]; [-1; -1]; [1; -2]; [2; -1]}];
                                %{[3/2; sqrt(3)/2]; [0; sqrt(3)]; [-3/2; sqrt(3)/2]; ...
                                %[-3/2; -sqrt(3)/2]; [0; -sqrt(3)]; [3/2; -sqrt(3)/2]}];
                        case "g1_tri"
                            %hObj.b_vecs = [0 3/4;sqrt(3)/2 -sqrt(3)/4];
                            hObj.b_vecs = [-1/2 1;sqrt(3)/2 0];
                            hObj.F = [hObj.F; {[1;0]; [0; 1]; [-1; -1]}];
                            %hObj.F = [hObj.F; {[1;0]; [1; 1]; [0; 1]; [-1; 0]; [-1; -1]; [0; -1]}];
                        otherwise
                            error("Unknown Fourier space name");
                    end
                end
                %hObj.F = int16([hObj.F{:}].');
                hObj.F = ([hObj.F{:}].');
                hObj.F_len = length(hObj.F);
            elseif isfloat(Fourier_space) & (lower(point_group) == "c6v" ...
                    || lower(point_group) == "c6" || lower(point_group) == "c3v" ...
                    || lower(point_group) == "c3") 
                 % a hex array of values
                hObj.F = Fourier_space;
                hObj.F_len = length(hObj.F);
            elseif isfloat(Fourier_space) % a square array of values
                sz = size(Fourier_space);
                hObj.F_len = sz(1)*sz(2)*sz(3);
                %hObj.F = reshape(int16(Fourier_space),[hObj.F_len 2]);
                hObj.F = reshape((Fourier_space),[hObj.F_len 3]);
                %hObj.F = Fourier_space;
                %hObj.F_len = length(hObj.F);
            else % Fourier_space is a cell array of [m,n]
                if isrow(Fourier_space)
                    %hObj.F = int16([Fourier_space{:}].');
                    hObj.F = ([Fourier_space{:}].');
                else
                    %hObj.F = int16([Fourier_space{:}]);
                    hObj.F = ([Fourier_space{:}]);
                end
                hObj.F_len = length(hObj.F);
            end

            hObj.point_group = point_group;

            % List the symmetry operations in each group
            switch lower(point_group)
                case "c1"
                    hObj.ops = "E_1";
                    hObj.projs = "A";
                case "cs"
                    hObj.ops = ["E_1" "sigma_d_d"];
                    hObj.projs = ["Ap" "App"];
                case "cs_a"
                    hObj.ops = ["E_1" "sigma_d_a"];
                    hObj.projs = ["Ap" "App"];
                case "cs_hex"
                    hObj.ops = ["E_1" "sigma_d_1"];
                    hObj.projs = ["Ap" "App"];
                case "c1v"
                    hObj.ops = ["E_1" "sigma_v_x"];
                    hObj.projs = ["A1" "A2"];
                case "c1v_vert"
                    hObj.ops = ["E_1" "sigma_v_y"];
                    hObj.projs = ["A1" "A2"];
                case "c2"
                    hObj.ops = ["E_1" "C2_z"];
                    hObj.projs = ["A" "B"];
                case "c2v"
                    hObj.ops = ["E_1" "C2_z" "sigma_v_x" "sigma_v_y"];
                    hObj.projs = ["A1" "A2" "B1" "B2"];
                case "c2v_d"
                    hObj.ops = ["E_1" "C2_z" "sigma_d_d" "sigma_d_a"];
                    hObj.projs = ["A1" "A2" "B1" "B2"];
                case "c4"
                    hObj.ops = ["E_1" "C4_1" "C4_3" "C2_1"];
                    hObj.projs = ["A" "B" "E_11" "E_22" "E"];
                case "c4v"
                    hObj.ops = ["E_1" "C4_1" "C4_3" "C2_1" "sigma_v_x" "sigma_v_y" "sigma_d_d" "sigma_d_a"];
                    hObj.projs = ["A1" "A2" "B1" "B2" "E_11" "E_22" "E"];
                case "c4v_d"
                    hObj.ops = ["E_1" "C4_1" "C4_3" "C2_1" "sigma_d_d" "sigma_d_a" "sigma_v_x" "sigma_v_y"];
                    hObj.projs = ["A1" "A2" "B1" "B2" "E_11" "E_22" "E"];
                    case "d2h"
                        hObj.ops = ["E_1" "C2_1" "C2p_y" "C2p_x" ...
                            "i_op" "sigma_h" "sigma_v_x" "sigma_v_y"];
                        hObj.projs = ["Ag" "B1g" "B2g" "B3g" "Au" "B1u" "B2u" "B3u"];
                    case "d2h_d"
                        hObj.ops = ["E_1" "C2_1" "C2pp_a" "C2pp_d" ...
                            "i_op" "sigma_h" "sigma_d_d" "sigma_d_a"];
                        hObj.projs = ["Ag" "B1g" "B2g" "B3g" "Au" "B1u" "B2u" "B3u"];
                    case "d4h"
                        hObj.ops = ["E_1" "C4_1" "C4_3" "C2_1" "C2p_y" "C2p_x" "C2pp_a" "C2pp_d" ...
                            "i_op" "S4_1" "S4_3" "sigma_h" "sigma_v_x" "sigma_v_y" "sigma_d_d" "sigma_d_a"];
                        hObj.projs = ["A1g" "A2g" "B1g" "B2g" "E_11g" "E_22g" "Eg" ...
                            "A1u" "A2u" "B1u" "B2u" "E_11u" "E_22u" "Eu"];
                    case "d4h_d"
                        hObj.ops = ["E_1" "C4_1" "C4_3" "C2_1" "C2pp_a" "C2pp_d" "C2p_y" "C2p_x" ...
                            "i_op" "S4_1" "S4_3" "sigma_h" "sigma_d_d" "sigma_d_a" "sigma_v_x" "sigma_v_y"];
                        hObj.projs = ["A1g" "A2g" "B1g" "B2g" "E_11g" "E_22g" "Eg" ...
                            "A1u" "A2u" "B1u" "B2u" "E_11u" "E_22u" "Eu"];
                    case "s4"
                        hObj.ops = ["E_1" "S4_1" "S4_3" "C2_1"];
                        hObj.projs = ["A" "B" "E_11" "E_22" "E"];
                otherwise
                    error("Projector: Unkown group")
            end

            for p=hObj.projs
                addprop(hObj, p);
                addprop(hObj, strcat(p, "_calculated"));
                hObj.(strcat(p, "_calculated")) = false;
            end
        end

        function create_operator_matrices(hObj)

            hObj.C4_3 = hObj.C4_1^3;
            hObj.C2_1 = hObj.C4_1^2;
            hObj.sigma_v_y = hObj.sigma_v_x*hObj.C2_1;
            hObj.sigma_d_d = hObj.C4_1*hObj.sigma_v_x;
            hObj.sigma_d_a = hObj.C4_3*hObj.sigma_v_x;
            hObj.S4_1 = hObj.i_op*hObj.C4_1;
            hObj.S4_3 = hObj.i_op*hObj.C4_3;
            hObj.sigma_h = hObj.i_op*hObj.C2_1;
            hObj.C2p_y = hObj.i_op*hObj.sigma_v_x;
            hObj.C2p_x = hObj.i_op*hObj.sigma_v_y;
            hObj.C2pp_a = hObj.i_op*hObj.sigma_d_d;
            hObj.C2pp_d = hObj.i_op*hObj.sigma_d_a;
        end

        %%
        %
        % Sundry functions
        %

        function eta = symmetry_param(hObj, irrep, vec)
            basis_range = hObj.(irrep)(1):hObj.(irrep)(end);
            pvec = hObj.U(:,basis_range)'*vec;
            eta = pvec'*pvec/(vec'*vec);
        end

        function eta = symmetry_param_vec_field(hObj, irrep, irrep_complement, vec_x, vec_y, vec_z)
            basis_range = hObj.(irrep)(1):hObj.(irrep)(end);
            basis_range_complement = hObj.(irrep_complement)(1):hObj.(irrep_complement)(end);
            pvec_x = hObj.U(:,basis_range)'*vec_x;
            pvec_y = hObj.U(:,basis_range_complement)'*vec_y;
            pvec_z = hObj.U(:,basis_range)'*vec_z;
            eta = (pvec_x'*pvec_x+pvec_y'*pvec_y+pvec_z'*pvec_z)/(vec_x'*vec_x+vec_y'*vec_y+vec_z'*vec_z);
        end

        function proj = projector(hObj, irrep)
            if isempty(hObj.(irrep))
                proj = [];
            else
                basis_range = hObj.(irrep)(1):hObj.(irrep)(end);
                proj = hObj.U(:,basis_range)*hObj.U(:,basis_range)';
            end
        end

        % Produce a vector of a given symmetry
        function vec = vector(hObj, irrep)
            vec = zeros(size(hObj.U,1),1);
            basis_range = hObj.(irrep)(1):hObj.(irrep)(end);
            vec(basis_range) = 1;
        end

        function basis = basis(hObj, varargin)
            if nargin == 1
                basis = hObj.F;
            else
                basis_range = hObj.(varargin{1})(1):hObj.(varargin{1})(end);
                basis = hObj.U(:,basis_range);
            end
        end

        function dim_of_irrep = dim(hObj, irrep)
            if isscalar(hObj.(irrep))
                dim_of_irrep = 1;
            else
                dim_of_irrep = hObj.(irrep)(:,2) - hObj.(irrep)(:,1)+1;
            end
        end

        function space = fourier_space(hObj)
            space = hObj.F;
        end
        
        function out = group(hObj)
            out = hObj.point_group;
        end


        function out = get_char_table(hObj, group_name)
            % http://gernot-katzers-spice-pages.com/character_tables/index.html
            %
            % R. McWeeny, Symmetry; an introduction to group theory and 
            % its applications. (Pergamon Press; distributed in the 
            % Western Hemisphere by Macmillan, Oxford, New York,, 1963).

            switch lower(group_name)
                case "c1"
                    out = 1;
                case {"cs", "cs_a", "cs_hex", "c1v", "c1v_vert" "c2"}
                    out = [1  1; ...
                        1 -1];
                case {"c2v", "c2v_d"}
                    out = [1  1  1  1;...
                           1  1 -1 -1;...
                           1 -1  1 -1;...
                           1 -1 -1  1];
                case {"c4" "s4"}
                    out = [1  1  1  1;...
                           1 -1  -1   1;...
                           1  1i -1i -1;...
                           1 -1i  1i -1;...
                           2  0   0  -2];
                case {"c4v", "c4v_d"}
                    %      1  2  3  4  5  6  7  8
                    out = [1  1  1  1  1  1  1  1;...
                           1  1  1  1 -1 -1 -1 -1;...
                           1 -1 -1  1  1  1 -1 -1;...
                           1 -1 -1  1 -1 -1  1  1;...

                           % These rows are out of order since we need to
                           % do E_11 and E_22 first so that we can
                           % normalise E_12 and E_21
                           % E_11; E_22; E_12; E_21; E
                           1  0  0 -1  1 -1  0  0;...
                           1  0  0 -1 -1  1  0  0;...
                           %0 -1  1  0  0  0  1 -1;...
                           %0  1 -1  0  0  0  1 -1;...
                           2  0  0 -2  0  0  0  0];
                case {"d2h" "d2h_d"}
                    d2 = hObj.get_char_table("C2v");
                    out = [d2 d2;d2 -d2];
                case {"d4h" "d4h_d"}
                    d4 = hObj.get_char_table("C4v");
                    out = [d4 d4;d4 -d4];
                otherwise
                    error("Projector.get_char_table(): Unkown group")
            end
        end

        function add_irrep(hObj, irrep)
        
            % Only run if it hasn't already been calculated
            if ~hObj.(strcat(irrep, "_calculated"))
                char_table = hObj.get_char_table(hObj.point_group);

                switch irrep
                    case "E"
                        hObj.add_irrep("E_11");
                        hObj.add_irrep("E_22");
                        hObj.E = [hObj.E_11; hObj.E_22];
                    case "E1"
                        hObj.add_irrep("E1_11");
                        hObj.add_irrep("E1_22");
                        hObj.E1 = [hObj.E1_11; hObj.E1_22];
                    case "E2"
                        hObj.add_irrep("E2_11");
                        hObj.add_irrep("E2_22");
                        hObj.E2 = [hObj.E2_11; hObj.E2_22];
                    case "Eg"
                        hObj.add_irrep("E_11g");
                        hObj.add_irrep("E_22g");
                        hObj.Eg = [hObj.E_11g; hObj.E_22g];
                    case "Eu"
                        hObj.add_irrep("E_11u");
                        hObj.add_irrep("E_22u");
                        hObj.Eu = [hObj.E_11u; hObj.E_22u];
                    case "U"
                        for p = hObj.projs(1:end)
                            hObj.add_irrep(p);
                        end
                        hObj.U = sparse(hObj.U);
                        hObj.U_calculated = true;
                    otherwise
                        K_rep = sparse(hObj.F_len, 2*hObj.F_len);
                        F_vec = sparse(zeros(size(hObj.F)));
                        % F_tmp has F as columns repeated as many times as there are
                        % elements of F.
                        F_tmp = sparse(hObj.F.');
                        F_tmp = repmat(F_tmp(:).',hObj.F_len,1);
            
                        P_tmp = sparse(hObj.F_len, hObj.F_len);
        
                        irrep_num = find(hObj.projs == irrep);
    
                        %representation = sparse(zeros(hObj.F_len, hObj.F_len, 'logical'));
                        for sym_op = 1:length(hObj.ops)
                            get_representation;
                            %{
                            if ~isprop(hObj, strcat(hObj.ops(sym_op), "_representation"))
                                addprop(hObj, strcat(hObj.ops(sym_op), "_representation"));
                                hObj.(strcat(hObj.ops(sym_op), "_representation")) = representation;
                            else
                                representation = hObj.(strcat(hObj.ops(sym_op), "_representation"));
                            end
                            %}
                            %add_representation(char_table(:,sym_op));
                        end
        
                        [basis, pivots] = create_canonical_basis;
    
                        hObj.(irrep) = pivots+size(hObj.U,2);
                        hObj.U = [hObj.U basis];
                        hObj.U = sparse(hObj.U);
                end
                hObj.(strcat(irrep, "_calculated")) = true;
            end

            % Get a matrix representation of each symmetry operation in the
            % specified point group.
            function get_representation
                characters = char_table(:,sym_op);

                % K_out is F repeated so that there is one column of F for each
                % symmetry operation.  The corresponding symmetry operation is
                % then made on each column.

                % we are storing the Fourier space as the (p,q) values, but
                % we need to convert this to cartesian coordinates by
                % multiplying by the basis vectors - ie
                % p*b_1+q*b_2
                % and then convert back to (p,q) space afterwards.  This is
                % what the hObj.b_vecs in the code below does.
                F_vec = sparse((hObj.b_vecs\(hObj.(hObj.ops(sym_op))*...
                    (hObj.b_vecs*(hObj.F).'))).');


                % Create representations for each of the symmetry operations
                % ----------------------------------------------------------
                %
                % We are creating a matrix A such that A*F performs the relevant
                % symmetry operation on F as a vector.  To do this we take one
                % particular column of K_out, which represents a single symmetry
                % operation, and then repeat that column as, for example:
                %
                % [2;3;1] -> [2 2 2;3 3 3;1 1 1]
                %
                % and then compare this to F.'=[1 2 3] to get a matrix of the form
                %
                % [F.'; F.'; F.'] - [2 2 2;3 3 3;1 1 1] = [-1 0 1; -2 -1 0; 0 -1 -2]
                %
                % If we then take the negation of this result we get
                %
                % not([-1 0 1; -2 -1 0; 0 -1 -2]) = [0 1 0; 0 0 1; 1 0 0]
                %
                % which is the representation we were after since
                %
                % [0 1 0; 0 0 1; 1 0 0]*[1;2;3] = [2;3;1]
                %
                % Which matches the initial vector above.
                %
                % In reality it is a little more complicated since the actual form
                % of F has two columns, one for the x coord and one for the y coord
                % such as
                %
                % F = [1 0; 0 1; -1 0; 0 -1]
                %

                K_rep = repmat(F_vec, [1 hObj.F_len]);
                
                % Stop MATLAB from converting K_out into logical values,
                % which will ruin subsequent calcualtions.
                %K_rep = (F_tmp-K_rep)<1e-6&(F_tmp-K_rep)>-1e-6;
                K_rep = abs(F_tmp-K_rep)<1e-6;

                %representation(:,:) = K_rep(:,1:3:end)&K_rep(:,2:3:end)&K_rep(:,3:3:end);
                P_tmp = P_tmp + (K_rep(:,1:3:end)&K_rep(:,2:3:end)&K_rep(:,3:3:end))*characters(irrep_num);
            end
%{
            function add_representation(characters)
                P_tmp = P_tmp + representation*characters(irrep_num);
            end
%}

            function [basis, pivots] = create_canonical_basis

                % Reduce to row echelon form to get the matrix A such that
                % P = A*(A^T*A)^-1*A^T
                [basis, pivots] = frref(P_tmp);
                rnk = length(pivots);

                if rnk~= 0
                    basis=basis(1:rnk,:).';
                    % We need to apply Gram-Schmidt process to make it
                    % orthonormal so that we can get the relationship:
                    % Proj = basis*basis';
                    %[basis,~] = mgson(basis);
                    [basis, ~] = qr(basis,0);
                    if rnk > 1
                        pivots = [1, rnk];
                    else
                        pivots = 1;
                    end
                else
                    basis = [];
                end
            end
        end


    end

    methods (Access = protected)

        function propgrp = getPropertyGroups(hObj)
            propgrp = matlab.mixin.util.PropertyGroup([hObj.projs "U"]);
        end
    end
end

