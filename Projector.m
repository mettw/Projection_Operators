classdef Projector < dynamicprops & matlab.mixin.CustomDisplay
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
    % > P=Projector("G1", "C4v")
    %
    % (where "G1" refers to the \Gamma^{(1)} points) or a set of standard 
    % bases (Note the use of square brackets in this case)
    %
    % > P=Projector(["G0" "G1" "G2"], "C4v")
    %
    % and for C6v 
    %
    % > P=Projector("G1_hex", "C6v")
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
        % Hilbert space
        % These are the (p,q) values for the vectors
        % va{k} = p\va{b}_1+q\va{b}_2
        F;
        F_len;
        % reciprocal space basis vectors (square by default)
        b_vecs = [1 0;0 1];

        save_representations = false;

        projs; % Names of the projectors to create
        ops; % Names of the symmetry operations
        y_vec; % Used to create operator for y axis projector.
        include_y_projs = false;

        % Point group operators
        E_1 = ([1 0;0 1]);
        C6_1 = ([1/2 -sqrt(3)/2;sqrt(3)/2 1/2]);
        C6_5 = ([1/2 sqrt(3)/2;-sqrt(3)/2 1/2]);
        C4_1 = ([0 -1;1 0]);
        C4_3 = ([0 1;-1 0]);
        C3_1 = ([-1/2 -sqrt(3)/2;sqrt(3)/2 -1/2]);
        C3_2 = ([-1/2 sqrt(3)/2;-sqrt(3)/2 -1/2]);
        C2_1 = ([-1 0;0 -1]);
        sigma_v_x = ([1 0;0 -1]);
        sigma_v_y = ([-1 0;0 1]);
        sigma_d_d = ([0 1;1 0]);
        sigma_d_a = ([0 -1;-1 0]);
        % Hexagonal mirror symmetries
        sigma_v_1 = ([1 0;0 -1]);
        sigma_v_2 = ([1/4-3/4 2*1/2*sqrt(3)/2;2*1/2*sqrt(3)/2 3/4-1/4]);
        sigma_v_3 = ([1/4-3/4 -2*1/2*sqrt(3)/2;-2*1/2*sqrt(3)/2 3/4-1/4]);
        sigma_d_1 = ([3/4-1/4 2*1/2*sqrt(3)/2;2*1/2*sqrt(3)/2 1/4-3/4]);
        sigma_d_2 = ([-1 0;0 1]);
        sigma_d_3 = ([3/4-1/4 -2*1/2*sqrt(3)/2;-2*1/2*sqrt(3)/2 1/4-3/4]);
        %"mE_1" "mC4_1" "mC4_3" "mC2_1" "msigma_v_x" "msigma_v_y" "msigma_d_d" "msigma_d_a"
        % Negation of the operators for groups with inversion, since in the
        % 2D plane inversion is equivalent to a negation of the
        % coordinates.  ie, in 2D you multiply by [-1 0; 0 -1] == -I.
        mE_1 = -([1 0;0 1]);
        mC4_1 = -([0 -1;1 0]);
        mC4_3 = -([0 1;-1 0]);
        mC2_1 = -([-1 0;0 -1]);
        msigma_v_x = -([1 0;0 -1]);
        msigma_v_y = -([-1 0;0 1]);
        msigma_d_d = -([0 1;1 0]);
        msigma_d_a = -([0 -1;-1 0]);

        point_group;
        all_representations = [];
    end
    
    methods
        %
        % SETUP functions
        %
        
        function hObj = Projector(Hilbert_space, point_group, vararg)
            
            if nargin == 3
                hObj.include_y_projs = vararg{1};
            end
            % if Hilbert_space is a name of a standard F_i space
            if isstring(Hilbert_space)
                for F_i=Hilbert_space
                    switch lower(F_i)
                        case "g0"
                            hObj.F = [hObj.F; {[0; 0]}];
                        case "g1"
                            hObj.F = [hObj.F; {[1; 0]; [0; 1]; [-1; 0]; [0; -1]}];
                        case "g2"
                            hObj.F = [hObj.F; {[1; 1]; [-1; 1]; [-1; -1]; [1; -1]}];
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
                            error("Unknown Hibert space name");
                    end
                end
                %hObj.F = int16([hObj.F{:}].');
                hObj.F = ([hObj.F{:}].');
                hObj.F_len = length(hObj.F);
            elseif isfloat(Hilbert_space) & (lower(point_group) == "c6v" ...
                    || lower(point_group) == "c6" || lower(point_group) == "c3v" ...
                    || lower(point_group) == "c3") 
                 % a hex array of values
                hObj.F = Hilbert_space;
                hObj.F_len = length(hObj.F);
            elseif isfloat(Hilbert_space) % a square array of values
                sz = size(Hilbert_space);
                hObj.F_len = sz(1)*sz(2);
                %hObj.F = reshape(int16(Hilbert_space),[hObj.F_len 2]);
                hObj.F = reshape((Hilbert_space),[hObj.F_len 2]);
                %hObj.F = Hilbert_space;
                %hObj.F_len = length(hObj.F);
            else % Hilbert_space is a cell array of [m,n]
                if isrow(Hilbert_space)
                    %hObj.F = int16([Hilbert_space{:}].');
                    hObj.F = ([Hilbert_space{:}].');
                else
                    %hObj.F = int16([Hilbert_space{:}]);
                    hObj.F = ([Hilbert_space{:}]);
                end
                hObj.F_len = length(hObj.F);
            end

            hObj.point_group = point_group;

            % List the symmetry operations in each group
            switch lower(point_group)
                case "c1"
                    hObj.ops = "E_1";
                    hObj.projs = ["A" "U"];
                case "cs"
                    hObj.ops = ["E_1" "sigma_d_d"];
                    hObj.projs = ["Ap" "App" "U"];
                case "cs_a"
                    hObj.ops = ["E_1" "sigma_d_a"];
                    hObj.projs = ["Ap" "App" "U"];
                case "cs_hex"
                    hObj.ops = ["E_1" "sigma_d_1"];
                    hObj.projs = ["Ap" "App" "U"];
                case "c1v"
                    hObj.ops = ["E_1" "sigma_v_x"];
                    hObj.y_vec = [1 -1];
                    hObj.projs = ["A1" "A2" "U"];
                case "c1v_vert"
                    hObj.ops = ["E_1" "sigma_v_y"];
                    hObj.y_vec = [1 -1];
                    hObj.projs = ["A1" "A2" "U"];
                case "c2"
                    hObj.ops = ["E_1" "C2_1"];
                    hObj.projs = ["A" "B" "U"];
                case "c2v"
                    hObj.ops = ["E_1" "C2_1" "sigma_v_x" "sigma_v_y"];
                    hObj.y_vec = [1 1 -1 -1];
                    hObj.projs = ["A1" "A2" "B1" "B2" "U"];
                case "c2v_d"
                    hObj.ops = ["E_1" "C2_1" "sigma_d_d" "sigma_d_a"];
                    hObj.y_vec = [1 1 -1 -1];
                    hObj.projs = ["A1" "A2" "B1" "B2" "U"];
                case "c3v"
                    hObj.ops = ["E_1" "C3_1" "C3_2" "sigma_v_1" "sigma_v_2" "sigma_v_3"];
                    hObj.y_vec = [1 1 1 -1 -1 -1];
                    hObj.projs = ["A1" "A2" "E_11" "E_12" "E_21" "E_22" "E" "U"];
                case "c4"
                    hObj.ops = ["E_1" "C4_1" "C4_3" "C2_1"];
                    hObj.projs = ["A" "B" "E_11" "E_22" "E" "U"];
                case "c4v"
                    hObj.ops = ["E_1" "C4_1" "C4_3" "C2_1" "sigma_v_x" "sigma_v_y" "sigma_d_d" "sigma_d_a"];
                    hObj.y_vec = [1 1 1 1 -1 -1 -1 -1];
                    hObj.projs = ["A1" "A2" "B1" "B2" "E_11" "E_12" "E_21" "E_22" "E" "U"];
                case "c4v_d"
                    hObj.ops = ["E_1" "C4_1" "C4_3" "C2_1" "sigma_d_d" "sigma_d_a" "sigma_v_x" "sigma_v_y"];
                    hObj.y_vec = [1 1 1 1 -1 -1 -1 -1];
                    hObj.projs = ["A1" "A2" "B1" "B2" "E_11" "E_12" "E_21" "E_22" "E" "U"];
                case "c6v"
                    hObj.ops = ["E_1" "C6_1" "C6_5" "C3_1" "C3_2" "C2_1" "sigma_v_1" "sigma_v_2" "sigma_v_3" ...
                        "sigma_d_1" "sigma_d_2" "sigma_d_3"];
                    hObj.y_vec = [1 1 1 1 1 1 -1 -1 -1 -1 -1 -1];
                    hObj.projs = ["A1" "A2" "B1" "B2" ...
                        "E1_11" "E1_12" "E1_21" "E1_22" "E1" ...
                        "E2_11" "E2_12" "E2_21" "E2_22" "E2" "U"];
                case "d2h"
                    hObj.ops = ["E_1" "C2_1" "sigma_v_x" "sigma_v_y" ...
                        "mE_1" "mC2_1" "msigma_v_x" "msigma_v_y"];
                    hObj.projs = ["Ag" "B1g" "B2g" "B3g" "Au" "B1u" "B2u" "B3u" "U"];
                case "d4h"
                    hObj.ops = ["E_1" "C4_1" "C4_3" "C2_1" "sigma_v_x" "sigma_v_y" "sigma_d_d" "sigma_d_a" ...
                        "mE_1" "mC4_1" "mC4_3" "mC2_1" "msigma_v_x" "msigma_v_y" "msigma_d_d" "msigma_d_a"];
                    hObj.projs = ["A1g" "A2g" "B1g" "B2g" "E_11g" "E_12g" "E_21g" "E_22g" "Eg" ...
                        "A1u" "A2u" "B1u" "B2u" "E_11u" "E_12u" "E_21u" "E_22u" "Eu" "U"];
                otherwise
                    error("Projector: Unkown group")
            end

            for p=hObj.projs
                addprop(hObj, p);
            end
        end

        function add_irrep(hObj, irrep)
            % 'true' means normalise the matrix
            hObj.add_irrep_private(irrep, true);
        end

        % The method add_irrep() takes a lot of time because it uses
        % reduction to row echelon form (ffref()) to normalise the
        % matrices.  While this has the advantage of finding the basis for
        % the projectors we usually don't need to know what the basis is
        % and so it is a waste of resources to find the basis every time.
        % It is especially a problem for projectors with a large basis (say
        % 150x150) since it takes too long to create the projector.
        %
        % I am therefore adding this method that uses an alternative way to
        % normalise the projectors without having to find the basis.
        function add_irrep_fast(hObj, irrep)
            if irrep == "U"
                error("Can't add U with add_irrep_fast() as that requires finding the basis.");
            end

            % 'false' means don't normalise the matrix
            hObj.add_irrep_private(irrep, false);

            % Fast normalisation
            % ------------------
            %
            % If the normalisation factor is 'g' and the non-normalised
            % matrix is 'B' then, for any i
            % 
            % g = B(i,i)/(B(i,:)*B(:,i))
            %
            % but we need to make sure that the value of i chosen is one
            % for which the result is non-zero, so we get

            g = hObj.(irrep)./(hObj.(irrep)*hObj.(irrep));
            g = g(find(~isnan(g) & g ~= 0, 1)); % g becomes a scalar here

            hObj.(irrep) = g*hObj.(irrep);

        end

        %%
        %
        % Sundry functions
        %

        function basis = basis(hObj, varargin)
            if nargin == 1
                basis = hObj.F;
            else
                basis = hObj.(strcat(varargin{1}, "_basis"));
            end
        end

        function space = hilbert_space(hObj)
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
                case "c4"
                    out = [1  1  1  1;...
                           1 -1  -1   1;...
                           1  1i -1i -1;...
                           1 -1i  1i -1;...
                           2  0   0  -2];
                case "c3v"
                    c=1/2;
                    s=sqrt(3)/2;
                    %      1  2  3  4  5  6 
                    out = [1  1  1  1  1  1;...
                           1  1  1  -1 -1 -1 ;...

                           1 -c -c  1 -c -c;...
                           0 -s  s  0 -s  s;...
                           0  s -s  0 -s  s;...
                           1 -c -c -1  c  c;...
                           2 -1 -1  0  0  0 ];
                case {"c4v", "c4v_d"}
                    %      1  2  3  4  5  6  7  8
                    out = [1  1  1  1  1  1  1  1;...
                           1  1  1  1 -1 -1 -1 -1;...
                           1 -1 -1  1  1  1 -1 -1;...
                           1 -1 -1  1 -1 -1  1  1;...

                           1  0  0 -1  1 -1  0  0;...
                           0 -1  1  0  0  0  1 -1;...
                           0  1 -1  0  0  0  1 -1;...
                           1  0  0 -1 -1  1  0  0;...
                           2  0  0 -2  0  0  0  0];
                case "c6v"
                    % "E_1" "C6_1" "C6_5" "C3_1" "C3_2" "C2_1" "sigma_v_1" "sigma_v_2" "sigma_v_3" ...
                    % "sigma_d_1" "sigma_d_2" "sigma_d_3"
                    c=1/2;
                    s=sqrt(3)/2;
                    out = [1    1   1     1    1    1    1  1  1    1  1  1;...
                           1    1   1     1    1    1   -1 -1 -1   -1 -1 -1;...
                           1   -1  -1     1    1   -1    1  1  1   -1 -1 -1;...
                           1   -1  -1     1    1   -1   -1 -1 -1    1  1  1;...

                           1    c   c    -c   -c   -1    1 -c -c    c -1  c;...
                           0   -s   s    -s    s    0    0  s -s    s  0 -s;...
                           0    s  -s     s   -s    0    0  s -s    s  0  s;...
                           1    c   c    -c   -c   -1   -1  c  c   -c  1 -c;...
                           2  2*c 2*c  -2*c -2*c   -2    0  0  0    0  0  0;...

                           1   -c   -c   -c   -c    1    1 -c -c   -c  1 -c;...
                           0    s   -s   -s    s    0    0  s -s   -s  0  s;...
                           0   -s    s    s   -s    0    0  s -s   -s  0  s;...
                           1   -c   -c   -c   -c    1   -1  c  c    c -1  c;...
                           2 -2*c -2*c -2*c -2*c    2    0  0  0    0  0  0];
                case "d2h"
                    d2 = hObj.get_char_table("C2v");
                    out = [d2 d2;d2 -d2];
                case "d4h"
                    d4 = hObj.get_char_table("C4v");
                    out = [d4 d4;d4 -d4];
                otherwise
                    error("Projector.get_char_table(): Unkown group")
            end
        end


    end

    methods (Access = protected)

        function add_irrep_private(hObj, irrep, normalise)
        
            char_table = hObj.get_char_table(hObj.point_group);

            % Only run if it hasn't already been run
            if isempty(hObj.(irrep))
                %K_rep = ones(hObj.F_len, 2*hObj.F_len, 'int16');
                K_rep = ones(hObj.F_len, 2*hObj.F_len);
                F_vec = ones(size(hObj.F));
                % F_tmp has F as columns repeated as many times as there are
                % elements of F.
                F_tmp = (hObj.F.');
                F_tmp = repmat(F_tmp(:).',hObj.F_len,1);

                % Don't allocate any memory for U
                if irrep ~= "U"
                    irrep_num = find(hObj.projs == irrep);
                    %hObj.(irrep) = zeros(hObj.F_len, hObj.F_len, 'single');
                    hObj.(irrep) = zeros(hObj.F_len, hObj.F_len);
                    representation = zeros(hObj.F_len, hObj.F_len, 'logical');
                    for sym_op = 1:length(hObj.ops)
                        get_representation;
                        add_representation(char_table(:,sym_op));
                    end
                    if normalise
                        create_canonical_basis;
                    end
                    % I got the orientaion of the matrices wrong - which
                    % doesn't matter for a real projector, since they are
                    % symmetric.  But the Hermitian projectors of C4 E11 
                    % and E22 are oriented incorrectly, so I need to
                    % transpose the projector here.
                    hObj.(hObj.projs(irrep_num)) = hObj.(hObj.projs(irrep_num)).';
                else
                    % allocate memory for all not already calculated
                    for p = hObj.projs(1:end-1)
                        irrep_num = find(hObj.projs == p);
                        if isempty(hObj.(p))
                            %hObj.(p) = zeros(hObj.F_len, hObj.F_len, 'single');
                            hObj.(p) = zeros(hObj.F_len, hObj.F_len);
                            representation = zeros(hObj.F_len, hObj.F_len, 'logical');
                            for sym_op = 1:length(hObj.ops)
                                get_representation;
                                add_representation(char_table(:,sym_op));
                            end
                            if normalise
                                create_canonical_basis;
                            end
                            % I got the orientaion of the matrices wrong - which
                            % doesn't matter for a real projector, since they are
                            % symmetric.  But the Hermitian projectors of C4 E11 
                            % and E22 are oriented incorrectly, so I need to
                            % transpose the projector here.
                            hObj.(hObj.projs(irrep_num)) = hObj.(hObj.projs(irrep_num)).';
                        end
                    end
                    % Now add U
                    irrep_num = length(hObj.projs);
                    if normalise
                        create_canonical_basis;
                    end
                end

            end

            function create_canonical_basis

                % if the user is asking for "U" then get all irr. reps.
                if irrep_num == length(hObj.projs)
                    %don't pass "U" to get_basis()
                    for pp=hObj.projs(1:end-1)
                        % Only run if it hasn't already been run
                        if isempty(hObj.(pp))
                            get_basis(pp);
                        end
                    end

                    switch lower(hObj.point_group)
                        case "c1"
                            hObj.U = hObj.A_basis;
                        case {"cs", "cs_a"}
                            hObj.U = [hObj.Ap_basis hObj.App_basis];
                        case {"c1v", "c1v_vert"}
                            hObj.U = [hObj.A1_basis hObj.A2_basis];
                        case "c2"
                            hObj.U = [hObj.A_basis hObj.B_basis];
                        case {"c2v", "c2v_d"}
                            hObj.U = [hObj.A1_basis hObj.A2_basis ...
                                hObj.B1_basis hObj.B2_basis];
                        case "c4"
                            hObj.U = [hObj.A_basis hObj.B_basis ...
                                hObj.E_basis];
                        case "c3v"
                            hObj.U = [hObj.A1_basis hObj.A2_basis ...
                                hObj.E_basis];
                        case {"c4v", "c4v_d"}
                            hObj.U = [hObj.A1_basis hObj.A2_basis ...
                                hObj.B1_basis hObj.B2_basis ...
                                hObj.E_11_basis hObj.E_22_basis];
                        case "c6v"
                            % apply gram-schmidt process to make orthogonal
                            [hObj.E1_11_basis,~] = gsog(hObj.E1_11_basis);
                            [hObj.E1_12_basis,~] = gsog(hObj.E1_12_basis);
                            [hObj.E1_21_basis,~] = gsog(hObj.E1_21_basis);
                            [hObj.E1_22_basis,~] = gsog(hObj.E1_22_basis);
                            [hObj.E1_basis,~] = gsog(hObj.E1_basis);

                            [hObj.E2_11_basis,~] = gsog(hObj.E2_11_basis);
                            [hObj.E2_12_basis,~] = gsog(hObj.E2_12_basis);
                            [hObj.E2_21_basis,~] = gsog(hObj.E2_21_basis);
                            [hObj.E2_22_basis,~] = gsog(hObj.E2_22_basis);
                            [hObj.E2_basis,~] = gsog(hObj.E2_basis);

                            hObj.U = [hObj.A1_basis hObj.A2_basis ...
                                hObj.B1_basis hObj.B2_basis ...
                                hObj.E1_basis hObj.E2_basis];
                        case "d2h"
                            hObj.U = [hObj.Ag_basis hObj.B1g_basis ...
                                hObj.B2g_basis hObj.B3g_basis ...
                                hObj.Au_basis hObj.B1u_basis hObj.B2u_basis ...
                                hObj.B3u_basis];
                        case "d4h"
                            hObj.U = [hObj.A1g_basis hObj.A2g_basis ...
                                hObj.B1g_basis hObj.B2g_basis ...
                                hObj.Eg_basis hObj.A1u_basis hObj.A2u_basis ...
                                hObj.B1u_basis hObj.B2u_basis ...
                                hObj.Eu_basis];
                    end
                else
                    get_basis(hObj.projs(irrep_num));
                    % If we need to apply Gram-Schmidt process to make
                    % orthonormal
                    if lower(hObj.group) == "c6v" && irrep_num > 4
                        [hObj.(strcat(hObj.projs(irrep_num), "_basis")),~] = ...
                            gsog(hObj.(strcat(hObj.projs(irrep_num), "_basis")));
                    end
                end
                
            end

            % Get a matrix representation of each symmetry operation in the
            % specified point group.
            function get_representation

                % K_out is F repeated so that there is one column of F for each
                % symmetry operation.  The corresponding symmetry operation is
                % then made on each column.
                %F_vec = single(hObj.F);

                %for vec_num = 1:hObj.F_len
                %    F_vec(vec_num,:) = (hObj.(hObj.ops(sym_op))*F_vec(vec_num,:).').';
                %end

                % we are storing the Hilbert space as the (p,q) values, but
                % we need to convert this to cartesian coordinates by
                % multiplying by the basis vectors - ie
                % p*b_1+q*b_2
                % and then convert back to (p,q) space afterwards.  This is
                % what the hObj.b_vecs in the code below does.
                %F_vec = int16(hObj.b_vecs\(hObj.(hObj.ops(sym_op))*...
                %    (hObj.b_vecs*single(hObj.F).'))).';
                F_vec = (hObj.b_vecs\(hObj.(hObj.ops(sym_op))*...
                    (hObj.b_vecs*(hObj.F).'))).';


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

                %K_rep(:,1:2:end) = repmat(F_vec(:,1), [1 hObj.F_len]);
                %K_rep(:,2:2:end) = repmat(F_vec(:,2), [1 hObj.F_len]);
                K_rep = repmat(F_vec, [1 hObj.F_len]);
                % Stop MATLAB from converting K_out into logical values,
                % which will ruin subsequent calcualtions.
                %K_rep = ~(F_tmp-K_rep);%<1e-6&(F_tmp-K_rep)>-1e-6);
                K_rep = (F_tmp-K_rep)<1e-6&(F_tmp-K_rep)>-1e-6;

                representation(:,:) = K_rep(:,1:2:end)&K_rep(:,2:2:end);
                if hObj.save_representations
                    hObj.all_representations = cat(3,hObj.all_representations, representation);
                end

            end

            function add_representation(characters)
                hObj.(hObj.projs(irrep_num)) = hObj.(hObj.projs(irrep_num)) + ...
                    representation*characters(irrep_num);
            end


            function get_basis(irr_rep)
                vec_str = strcat(irr_rep, "_basis");

                addprop(hObj, vec_str);

                [hObj.(vec_str), rnk] = frref(hObj.(irr_rep));
                rnk = length(rnk);

                if rnk~= 0
                    % normalise the operator
                    %
                    % The E_12 and E_21 matrices are not projectors since
                    % they are no idempotent and hence require special
                    % treatment as the trace is zero.
                    if irr_rep =="E_12" || irr_rep == "E_21"
                        if hObj.include_y_projs
                            hObj.add_irrep(irr_rep);
                            p_y = strcat(irr_rep, "_y");
                            addprop(hObj, p_y);
                            hObj.(p_y) = hObj.(irr_rep).*hObj.y_vec;
                            if irr_rep == "E_12"
                                    hObj.(p_y) = hObj.(p_y) ...
                                        /(trace(hObj.E_22)/rank(hObj.E_22));
                            else
                                    hObj.(p_y) = hObj.(p_y) ...
                                        /(trace(hObj.E_11)/rank(hObj.E_11));
                            end
                        end
                        hObj.add_irrep(irr_rep);
                            if irr_rep == "E_12"
                                hObj.(irr_rep) = hObj.(irr_rep) ...
                                    /(trace(hObj.E_22)/rank(hObj.E_22));
                            else
                                hObj.(irr_rep) = hObj.(irr_rep) ...
                                    /(trace(hObj.E_11)/rank(hObj.E_11));
                            end
                    else
                        if hObj.include_y_projs
                            p_y = strcat(irr_rep, "_y");
                            addprop(hObj, p_y);
                            hObj.(p_y) = hObj.(irr_rep).*hObj.y_vec;
                            hObj.(p_y) = hObj.(p_y) ...
                                /(trace(hObj.(p_y))/rank(hObj.(p_y)));
                        end
                        hObj.(irr_rep) = hObj.(irr_rep) ...
                            /(trace(hObj.(irr_rep))/rnk);

                    end
                    hObj.(vec_str)=hObj.(vec_str)(1:rnk,:);
                else
                    %hObj.(irr_rep) = [];
                    hObj.(vec_str) = [];
                end

                if rnk~= 0
                    hObj.(vec_str) = hObj.(vec_str).'./norm(hObj.(vec_str).');
                end
            end
        end

        function propgrp = getPropertyGroups(hObj)
            propgrp = matlab.mixin.util.PropertyGroup(hObj.projs);
        end
    end
end

