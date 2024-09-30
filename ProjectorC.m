classdef ProjectorC < dynamicprops & matlab.mixin.CustomDisplay
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
        projs; % Names of the projectors to create
        ops; % Names of the symmetry operations

        % Point group operators
        num_ops = 8;    
        E_1 = eye(8);
        C4_1 = [0 0 1 0 0 0 0 0; ...
                1 0 0 0 0 0 0 0; ...
                0 0 0 1 0 0 0 0; ...
                0 1 0 0 0 0 0 0; ...
                0 0 0 0 0 0 0 1; ...
                0 0 0 0 0 0 1 0; ...
                0 0 0 0 1 0 0 0; ...
                0 0 0 0 0 1 0 0];
        C4_3 = [0 1 0 0 0 0 0 0; ...
                0 0 0 1 0 0 0 0; ...
                1 0 0 0 0 0 0 0; ...
                0 0 1 0 0 0 0 0; ...
                0 0 0 0 0 0 1 0; ...
                0 0 0 0 0 0 0 1; ...
                0 0 0 0 0 1 0 0; ...
                0 0 0 0 1 0 0 0];
        C2_1 = [0 0 0 1 0 0 0 0; ...
              0 0 1 0 0 0 0 0; ...
              0 1 0 0 0 0 0 0; ...
              1 0 0 0 0 0 0 0; ...
              0 0 0 0 0 1 0 0; ...
              0 0 0 0 1 0 0 0; ...
              0 0 0 0 0 0 0 1; ...
              0 0 0 0 0 0 1 0];
        sigma_v_x = [0 0 0 0 1 0 0 0; ...
                     0 0 0 0 0 0 0 1; ...
                     0 0 0 0 0 0 1 0; ...
                     0 0 0 0 0 1 0 0; ...
                     1 0 0 0 0 0 0 0; ...
                     0 0 0 1 0 0 0 0; ...
                     0 0 1 0 0 0 0 0; ...
                     0 1 0 0 0 0 0 0];
        sigma_v_y = [0 0 0 0 0 1 0 0; ...
                     0 0 0 0 0 0 1 0; ...
                     0 0 0 0 0 0 0 1; ...
                     0 0 0 0 1 0 0 0; ...
                     0 0 0 1 0 0 0 0; ...
                     1 0 0 0 0 0 0 0; ...
                     0 1 0 0 0 0 0 0; ...
                     0 0 1 0 0 0 0 0];
        sigma_d_d = [0 0 0 0 0 0 1 0; ...
                     0 0 0 0 1 0 0 0; ...
                     0 0 0 0 0 1 0 0; ...
                     0 0 0 0 0 0 0 1; ...
                     0 1 0 0 0 0 0 0; ...
                     0 0 1 0 0 0 0 0; ...
                     1 0 0 0 0 0 0 0; ...
                     0 0 0 1 0 0 0 0];
        sigma_d_a = [0 0 0 0 0 0 0 1; ...
                     0 0 0 0 0 1 0 0; ...
                     0 0 0 0 1 0 0 0; ...
                     0 0 0 0 0 0 1 0; ...
                     0 0 1 0 0 0 0 0; ...
                     0 1 0 0 0 0 0 0; ...
                     0 0 0 1 0 0 0 0; ...
                     1 0 0 0 0 0 0 0];
        
        point_group;
        all_representations = [];
    end
    
    methods
        %
        % SETUP functions
        %
        
        function hObj = ProjectorC(point_group)
            
            if isstring(point_group)
                hObj.point_group = point_group;
    
                hObj.ops = ["E_1" "C4_1" "C4_3" "C2_1" "sigma_v_x" "sigma_v_y" "sigma_d_d" "sigma_d_a"];

                % List the symmetry operations in each group
                switch lower(point_group)
                    case "c1"
                        %hObj.ops = "E_1";
                        hObj.projs = ["A" "U"];
                    case "cs"
                        %hObj.ops = ["E_1" "sigma_d_d"];
                        hObj.projs = ["Ap" "App" "U"];
                    case "cs_a"
                        %hObj.ops = ["E_1" "sigma_d_a"];
                        hObj.projs = ["Ap" "App" "U"];
                    case "cs_hex"
                        %hObj.ops = ["E_1" "sigma_d_1"];
                        hObj.projs = ["Ap" "App" "U"];
                    case "c1v"
                        %hObj.ops = ["E_1" "sigma_v_x"];
                        hObj.projs = ["A1" "A2" "U"];
                    case "c1v_vert"
                        %hObj.ops = ["E_1" "sigma_v_y"];
                        hObj.projs = ["A1" "A2" "U"];
                    case "c2"
                        %hObj.ops = ["E_1" "C2_1"];
                        hObj.projs = ["A" "B" "U"];
                    case "c2v"
                        %hObj.ops = ["E_1" "C2_1" "sigma_v_x" "sigma_v_y"];
                        hObj.projs = ["A1" "A2" "B1" "B2" "U"];
                    case "c2v_d"
                        %hObj.ops = ["E_1" "C2_1" "sigma_d_d" "sigma_d_a"];
                        hObj.projs = ["A1" "A2" "B1" "B2" "U"];
                    case "c3v"
                        hObj.ops = ["E_1" "C3_1" "C3_2" "sigma_v_1" "sigma_v_2" "sigma_v_3"];
                        hObj.projs = ["A1" "A2" "E_11" "E_12" "E_21" "E_22" "E" "U"];
                    case "c4"
                        %hObj.ops = ["E_1" "C4_1" "C4_3" "C2_1"];
                        hObj.projs = ["A" "B" "E_11" "E_22" "E" "U"];
                    case "c4v"
                        %hObj.ops = ["E_1" "C4_1" "C4_3" "C2_1" "sigma_v_x" "sigma_v_y" "sigma_d_d" "sigma_d_a"];
                        hObj.projs = ["A1" "A2" "B1" "B2" "E_11" "E_12" "E_21" "E_22" "E" "U"];
                    case "c4v_d"
                        %hObj.ops = ["E_1" "C4_1" "C4_3" "C2_1" "sigma_d_d" "sigma_d_a" "sigma_v_x" "sigma_v_y"];
                        hObj.projs = ["A1" "A2" "B1" "B2" "E_11" "E_12" "E_21" "E_22" "E" "U"];
                    case "c6v"
                        hObj.ops = ["E_1" "C6_1" "C6_5" "C3_1" "C3_2" "C2_1" "sigma_v_1" "sigma_v_2" "sigma_v_3" ...
                            "sigma_d_1" "sigma_d_2" "sigma_d_3"];
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
                    hObj.(p) = zeros(hObj.num_ops);
                    addprop(hObj, strcat("M_", p));
                    hObj.(strcat("M_", p)) = zeros(hObj.num_ops);
                end
    
                hObj.add_irreps;
            elseif isnumeric(point_group)
                if iscolumn(point_group)
                    [n,m] = size(point_group);
                    if n==hObj.num_ops && m==1
                        addprop(hObj, "Operator");
                        hObj.Operator = zeros(hObj.num_ops);
                        hObj.ops = ["E_1" "C4_1" "C4_3" "C2_1" "sigma_v_x" ...
                            "sigma_v_y" "sigma_d_d" "sigma_d_a"];
                        hObj.projs = "Operator";
                        for op_num = 1:hObj.num_ops
                            hObj.Operator = hObj.Operator + ...
                                point_group(op_num)*hObj.(hObj.ops(op_num));
                        end
                    else
                        error("point_group wrong size");
                    end
                else
                    error("point_group not column vector");
                end
            else
                error("point_group neither string nor numeric");
            end
        end

        %%
        %
        % Sundry functions
        %

        function basis = basis(hObj, irrep)
            basis = hObj.(strcat(irrep, "_basis"));
        end

        function out = group(hObj)
            out = hObj.point_group;
        end

        % Return the character space vector representation of irrep.
        function vec = vector(hObj, irrep)
            char_tbl = hObj.get_char_table(hObj.point_group);
            vec = char_tbl(hObj.projs(1:end-1) == irrep,:).';
        end

        function out = get_char_table(hObj, group_name)
            % http://gernot-katzers-spice-pages.com/character_tables/index.html
            %
            % R. McWeeny, Symmetry; an introduction to group theory and 
            % its applications. (Pergamon Press; distributed in the 
            % Western Hemisphere by Macmillan, Oxford, New York,, 1963).

            switch lower(group_name)
                case "c1"
                    out = [1 0 0 0 0 0 0 0];
                case {"cs", "cs_a", "cs_hex"}
                    out = [1  1; ...
                        1 -1];
                % "E_1" "C4_1" "C4_3" "C2_1" "sigma_v_x" "sigma_v_y" "sigma_d_d" "sigma_d_a"
                case "c1v"
                    out = [1 0 0 0  1 0 0 0; ...
                           1 0 0 0 -1 0 0 0];
                case "c1v_vert"
                    out = [1 0 0 0 0  1 0 0; ...
                           1 0 0 0 0 -1 0 0];
                case "c2"
                    out = [1 0 0  1 0 0 0 0; ...
                           1 0 0 -1 0 0 0 0];
                case "c2v"
                    out = [1 0 0  1  1  1 0 0;...
                           1 0 0  1 -1 -1 0 0;...
                           1 0 0 -1  1 -1 0 0;...
                           1 0 0 -1 -1  1 0 0];
                case "c2v_d"
                    out = [1 0 0  1 0 0  1  1;...
                           1 0 0  1 0 0 -1 -1;...
                           1 0 0 -1 0 0  1 -1;...
                           1 0 0 -1 0 0 -1  1];
                case "c4"
                    out = [1  1  1    1 0 0 0 0;...
                           1 -1  -1   1 0 0 0 0;...
                           1  1i -1i -1 0 0 0 0;...
                           1 -1i  1i -1 0 0 0 0;...
                           2  0   0  -2 0 0 0 0];
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
                case "c4v"
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
                case "c4v_d"
                    %      1  2  3  4  5  6  7  8
                    out = [1  1  1  1  1  1  1  1;...
                           1  1  1  1 -1 -1 -1 -1;...
                           1 -1 -1  1 -1 -1  1  1;...
                           1 -1 -1  1  1  1 -1 -1;...

                           1  0  0 -1  0  0  1 -1;...
                           0 -1  1  0  1 -1  0  0;...
                           0  1 -1  0  1 -1  0  0;...
                           1  0  0 -1  0  0 -1  1;...
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

        function add_irreps(hObj)
        
            char_table = hObj.get_char_table(hObj.point_group);

            % allocate memory for all not already calculated
            for p = hObj.projs(1:end-1)
                irrep_num = find(hObj.projs == p);
                for sym_op = hObj.ops(1:end)
                    sym_num = find(hObj.ops == sym_op);
                    hObj.(p) = hObj.(p) + ...
                        char_table(irrep_num, sym_num)*hObj.(sym_op); %#ok<FNDSB>
                end
                hObj.(strcat("M_", p)) = hObj.(p);
                get_basis(p);
            end
            create_canonical_basis;

            function create_canonical_basis

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
                            hObj.E_11_basis hObj.E_22_basis];
                    case "c3v"
                        hObj.U = [hObj.A1_basis hObj.A2_basis ...
                            hObj.E_11_basis hObj.E_22_basis];
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
                            hObj.E1_11_basis hObj.E1_22_basis ...
                            hObj.E2_11_basis hObj.E2_22_basis];
                    case "d2h"
                        hObj.U = [hObj.Ag_basis hObj.B1g_basis ...
                            hObj.B2g_basis hObj.B3g_basis ...
                            hObj.Au_basis hObj.B1u_basis hObj.B2u_basis ...
                            hObj.B3u_basis];
                    case "d4h"
                        hObj.U = [hObj.A1g_basis hObj.A2g_basis ...
                            hObj.B1g_basis hObj.B2g_basis ...
                            hObj.E_11g_basis hObj.E_22g_basis ...
                            hObj.A1u_basis hObj.A2u_basis ...
                            hObj.B1u_basis hObj.B2u_basis ...
                            hObj.E_11u_basis hObj.E_22u_basis];
                end
                
            end
                               
            function get_basis(irr_rep)
                vec_str = strcat(irr_rep, "_basis");
                rnk_str = strcat(irr_rep, "_rank");
                tr_str = strcat(irr_rep, "_trace");

                addprop(hObj, vec_str);
                addprop(hObj, rnk_str);
                addprop(hObj, tr_str);

                [hObj.(vec_str), rnk] = rref(hObj.(irr_rep));
                rnk = length(rnk);
                hObj.(rnk_str) = rnk;

                if rnk~= 0
                    % normalise the operator
                    %
                    % The E_12 and E_21 matrices are not projectors since
                    % they are no idempotent and hence require special
                    % treatment as the trace is zero.
                    if irr_rep =="E_12" || irr_rep == "E_21"
                        if irr_rep == "E_12"
                            hObj.(tr_str) = trace(hObj.E_22);
                            hObj.(irr_rep) = hObj.(irr_rep) ...
                                /(trace(hObj.E_22)/rank(hObj.E_22));
                        else
                            hObj.(tr_str) = trace(hObj.E_11);
                            hObj.(irr_rep) = hObj.(irr_rep) ...
                                /(trace(hObj.E_11)/rank(hObj.E_11));
                        end
                    else
                        hObj.(tr_str) = trace(hObj.(irr_rep));
                        hObj.(irr_rep) = hObj.(irr_rep) ...
                            /(trace(hObj.(irr_rep))/rnk);

                    end
                    hObj.(vec_str)=hObj.(vec_str)(1:rnk,:);
                else
                    hObj.(vec_str) = [];
                end

                if rnk~= 0
                    if irr_rep =="E_12" || irr_rep == "E_21"
                        hObj.(vec_str) = hObj.(vec_str).'./norm(hObj.(vec_str).');
                    else
                        hObj.(vec_str) = orth(hObj.(irr_rep));
                    end
                end
            end
        end

        function propgrp = getPropertyGroups(hObj)
            propgrp = matlab.mixin.util.PropertyGroup(hObj.projs);
        end
    end
end

