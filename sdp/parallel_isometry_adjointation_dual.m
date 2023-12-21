% -------------------------------------------------------------------------
%
% File: parallel_isometry_adjointation_dual.m
%
% Description :
% Dual SDP code to obtain the optimal approximation error or success
% probability of transforming parallel 'n' calls of any isometry operations
% with input dimension 'din' and output dimension 'dout' into its inverse
% map, corresponding POVM, or adjoint map.
%
% -------------------------------------------------------------------------

function opt_fom = parallel_isometry_adjointation_dual(din,dout,n,task,isComplex)

N = n+1;

run('young_diagrams.m')

cvx_begin SDP quiet

cvx_solver sedumi

% ------------------------------------------------------------------
%          List of young diagrams with depth <= din or dout
% ------------------------------------------------------------------
for i = 1:N
    list_young_diagrams_in{i}={};
    list_young_diagrams_out{i}={};
    for alpha = 1:num_young_diagrams{i}
        if depth{i}{alpha}<=din
            list_young_diagrams_in{i}{end+1} = alpha;
        end
        if depth{i}{alpha}<=dout
            list_young_diagrams_out{i}{end+1} = alpha;
        end
    end
end

% ------------------------------------------------------------------
%          Calculate young diagrams by adding or removing a box
% ------------------------------------------------------------------
for i = 1:N
    list_child_in{i}={};
    list_parent_in{i}={};
    for index_alpha = 1:length(list_young_diagrams_in{i})
        alpha = list_young_diagrams_in{i}{index_alpha};
        if i<N
            list_child_in{i}{alpha}={};
            for index_mu=1:length(list_young_diagrams_in{i+1})
                mu = list_young_diagrams_in{i+1}{index_mu};
                if sum(X{i+1}{alpha}{mu},'all') ~= 0
                    list_child_in{i}{alpha}{end+1} = mu;
                end
            end
        end
        if i>1
            list_parent_in{i}{alpha}={};
            for index_gamma=1:length(list_young_diagrams_in{i-1})
                gamma = list_young_diagrams_in{i-1}{index_gamma};
                if sum(X{i}{gamma}{alpha},'all') ~= 0
                    list_parent_in{i}{alpha}{end+1} = gamma;
                end
            end
        end
    end
end

for i = 1:N
    list_child_out{i}={};
    list_parent_out{i}={};
    for index_alpha = 1:length(list_young_diagrams_out{i})
        alpha = list_young_diagrams_out{i}{index_alpha};
        if i<N
            list_child_out{i}{alpha}={};
            for index_mu=1:length(list_young_diagrams_out{i+1})
                mu = list_young_diagrams_out{i+1}{index_mu};
                if sum(X{i+1}{alpha}{mu},'all') ~= 0
                    list_child_out{i}{alpha}{end+1} = mu;
                end
            end
        end
        if i>1
            list_parent_out{i}{alpha}={};
            for index_gamma=1:length(list_young_diagrams_out{i-1})
                gamma = list_young_diagrams_out{i-1}{index_gamma};
                if sum(X{i}{gamma}{alpha},'all') ~= 0
                    list_parent_out{i}{alpha}{end+1} = gamma;
                end
            end
        end
    end
end

% ------------------------------------------------------------------
%          Declare SDP variables 
% ------------------------------------------------------------------
variable p;
variable omega;
variable xi;

for index_alpha = 1:length(list_young_diagrams_in{N-1})
    alpha = list_young_diagrams_in{N-1}{index_alpha};
    for index_nu = 1:length(list_young_diagrams_out{N})
        nu = list_young_diagrams_out{N}{index_nu};
        W{alpha}{nu} = semidefinite(dim{N-1}{alpha}*dim{N}{nu},isComplex);
    end
end

% ------------------------------------------------------------------
%          Parallel comb conditions
% ------------------------------------------------------------------
for index_mu = 1:length(list_young_diagrams_in{N})
    mu = list_young_diagrams_in{N}{index_mu};
    for index_nu = 1:length(list_young_diagrams_out{N})
        nu = list_young_diagrams_out{N}{index_nu};
        C{mu}{nu} = 0;
        for index_alpha = 1:length(list_parent_in{N}{mu})
            alpha = list_parent_in{N}{mu}{index_alpha};
            XI = Tensor(X{N}{alpha}{mu}, eye(dim{N}{nu}));
            C{mu}{nu} = C{mu}{nu} + multin{N}{mu}/multin{N-1}{alpha} * XI'*W{alpha}{nu}*XI;
        end
    end
end

for index_alpha = 1:length(list_young_diagrams_in{N-1})
    alpha = list_young_diagrams_in{N-1}{index_alpha};
    norm = 0;
    for index_nu = 1:length(list_young_diagrams_out{N})
        nu = list_young_diagrams_out{N}{index_nu};
        norm = norm + PartialTrace(W{alpha}{nu}, [2], [dim{N-1}{alpha} dim{N}{nu}]);
    end
    norm == p*multin{N-1}{alpha}*eye(dim{N-1}{alpha});
end

% ------------------------------------------------------------------
%          Conditions for the figure of merit
% ------------------------------------------------------------------
for index_mu = 1:length(list_young_diagrams_in{N})
    mu = list_young_diagrams_in{N}{index_mu};
    CI = Tensor(cycle{N}{mu},eye(dim{N}{mu}));
    Omega{mu} = CI*pure_to_mixed(MaxEntangled(dim{N}{mu}))*CI'*dim{N}{mu}*multin{N}{mu}/(din^2);
end

for index_mu = 1:length(list_young_diagrams_in{N})
    mu = list_young_diagrams_in{N}{index_mu};
    for index_nu = 1:length(list_young_diagrams_out{N})
        nu = list_young_diagrams_out{N}{index_nu};
        Sigma{mu}{nu} = 0;
        Xi{mu}{nu} = 0;
    end
end
for index_alpha = 1:length(list_young_diagrams_in{N-1})
    alpha = list_young_diagrams_in{N-1}{index_alpha};
    for index_mu = 1:length(list_child_in{N-1}{alpha})
        mu = list_child_in{N-1}{alpha}{index_mu};
        for index_nu = 1:length(list_child_in{N-1}{alpha})
            nu = list_child_in{N-1}{alpha}{index_nu};
            IC = Tensor(eye(dim{N}{mu}),cycle{N}{nu})';
            XX = Tensor(X{N}{alpha}{mu},X{N}{alpha}{nu});
            Xi{mu}{nu} = Xi{mu}{nu} + IC*XX'*pure_to_mixed(MaxEntangled(dim{N-1}{alpha}))*XX*IC'*dim{N-1}{alpha}*multin{N}{mu}*multin{N}{nu}/(din*multin{N-1}{alpha});
            if task == 3 | task == 4
                Sigma{mu}{nu} = Sigma{mu}{nu} - IC*XX'*pure_to_mixed(MaxEntangled(dim{N-1}{alpha}))*XX*IC'*dim{N-1}{alpha}*multin{N}{mu}*multin{N}{nu}/((dout-din)*multin{N-1}{alpha});
            end
        end
        for index_nu = 1:length(list_child_out{N-1}{alpha})
            nu = list_child_out{N-1}{alpha}{index_nu};
            IC = Tensor(eye(dim{N}{mu}),cycle{N}{nu})';
            XX = Tensor(X{N}{alpha}{mu},X{N}{alpha}{nu});
            if task == 3 | task == 4
                Sigma{mu}{nu} = Sigma{mu}{nu} + IC*XX'*pure_to_mixed(MaxEntangled(dim{N-1}{alpha}))*XX*IC'*dim{N-1}{alpha}*multin{N}{mu}*multout{N}{nu}/((dout-din)*multout{N-1}{alpha});
            end
        end
    end
end

switch task
    case 1
        for index_mu = 1:length(list_young_diagrams_in{N})
            mu = list_young_diagrams_in{N}{index_mu};
            C{mu}{mu} >= (1+omega)*Omega{mu} - omega*Xi{mu}{mu};
        end
        minimize p
        cvx_end
        cvx_status
        opt_fom = p;
    case 2
        for index_mu = 1:length(list_young_diagrams_in{N})
            mu = list_young_diagrams_in{N}{index_mu};
            C{mu}{mu} >= Omega{mu};
        end
        minimize p
        cvx_end
        cvx_status
        opt_fom = p;
    case 3
        for index_mu = 1:length(list_young_diagrams_in{N})
            mu = list_young_diagrams_in{N}{index_mu};
            for index_nu = 1:length(list_young_diagrams_out{N})
                nu = list_young_diagrams_out{N}{index_nu};
                    C{mu}{nu} >= xi*Xi{mu}{nu} - Sigma{mu}{nu};
            end
        end
        fom = xi - p;
        maximize fom
        cvx_end
        cvx_status
        opt_fom = fom;
    case 4
        for index_mu = 1:length(list_young_diagrams_in{N})
            mu = list_young_diagrams_in{N}{index_mu};
            for index_nu = 1:length(list_young_diagrams_out{N})
                nu = list_young_diagrams_out{N}{index_nu};
                if mu == nu
                    C{mu}{nu} >= omega*Omega{mu}+xi*Xi{mu}{nu} - (1-omega)*Sigma{mu}{nu};
                else
                    C{mu}{nu} >= xi*Xi{mu}{nu} - (1-omega)*Sigma{mu}{nu};
                end
            end
        end
        omega <= 1;
        fom = omega+xi-p;
        maximize fom
        cvx_end
        cvx_status
        opt_fom = fom;
end