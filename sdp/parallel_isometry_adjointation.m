% -------------------------------------------------------------------------
%
% File: parallel_isometry_adjointation.m
%
% Description :
% Primal SDP code to obtain the optimal approximation error or success
% probability of transforming parallel 'n' calls of any isometry operations
% with input dimension 'din' and output dimension 'dout' into its inverse
% map, corresponding POVM, or adjoint map.
%
% -------------------------------------------------------------------------

function opt_fom = parallel_isometry_adjointation(din,dout,n,task,isComplex)

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
for index_mu = 1:length(list_young_diagrams_in{N})
    mu = list_young_diagrams_in{N}{index_mu};
    for index_nu = 1:length(list_young_diagrams_out{N})
        nu = list_young_diagrams_out{N}{index_nu};
        S{mu}{nu} = semidefinite(dim{N}{mu}*dim{N}{nu},isComplex);
        if task == 2
            C{mu}{nu} = S{mu}{nu};
        else
            D{mu}{nu} = semidefinite(dim{N}{mu}*dim{N}{nu},isComplex);
            C{mu}{nu} = S{mu}{nu} + D{mu}{nu};
        end
    end
end

% ------------------------------------------------------------------
%          Parallel comb conditions
% ------------------------------------------------------------------
for index_alpha = 1:length(list_young_diagrams_in{N-1})
    alpha = list_young_diagrams_in{N-1}{index_alpha};
    E=0;
    for index_mu = 1:length(list_child_in{N-1}{alpha})
        mu = list_child_in{N-1}{alpha}{index_mu};
        for index_nu = 1:length(list_young_diagrams_out{N})
            nu = list_young_diagrams_out{N}{index_nu};
            XI = Tensor(X{N}{alpha}{mu},eye(dim{N}{nu}));
            E = E + PartialTrace(XI*C{mu}{nu}*XI', [2], [dim{N-1}{alpha} dim{N}{nu}]);
        end
    end
    for index_nu = 1:length(list_young_diagrams_out{N})
        LHS = 0;
        nu = list_young_diagrams_out{N}{index_nu};
        for index_mu = 1:length(list_child_in{N-1}{alpha})
            mu = list_child_in{N-1}{alpha}{index_mu};
            XI = Tensor(X{N}{alpha}{mu},eye(dim{N}{nu}));
            LHS = LHS + XI*C{mu}{nu}*XI'/multout{N}{nu};
        end
        LHS == Tensor(E,eye(dim{N}{nu}))/(dout^(n+1));
    end
end
norm = 0;
for index_mu = 1:length(list_young_diagrams_in{N})
    mu = list_young_diagrams_in{N}{index_mu};
    for index_nu = 1:length(list_young_diagrams_out{N})
        nu = list_young_diagrams_out{N}{index_nu};
        norm = norm +trace(C{mu}{nu});
    end
end
norm == dout^(n+1);

% ------------------------------------------------------------------
%          Conditions for the figure of merit
% ------------------------------------------------------------------
for index_mu = 1:length(list_young_diagrams_in{N})
    mu = list_young_diagrams_in{N}{index_mu};
    CI = Tensor(cycle{N}{mu},eye(dim{N}{mu}));
    Omega{mu} = CI*pure_to_mixed(MaxEntangled(dim{N}{mu}))*CI'*dim{N}{mu}/(din^2*multout{N}{mu});
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
            Sigma{mu}{nu} = Sigma{mu}{nu} + IC*XX'*pure_to_mixed(MaxEntangled(dim{N-1}{alpha}))*XX*IC'*dim{N-1}{alpha}*multin{N}{nu}/(multin{N-1}{alpha}*multout{N}{nu});
            Xi{mu}{nu} = Xi{mu}{nu} - IC*XX'*pure_to_mixed(MaxEntangled(dim{N-1}{alpha}))*XX*IC'*dim{N-1}{alpha}*multin{N}{nu}/(multin{N-1}{alpha}*multout{N}{nu});
        end
        for index_nu = 1:length(list_child_out{N-1}{alpha})
            nu = list_child_out{N-1}{alpha}{index_nu};
            IC = Tensor(eye(dim{N}{mu}),cycle{N}{nu})';
            XX = Tensor(X{N}{alpha}{mu},X{N}{alpha}{nu});
            Xi{mu}{nu} = Xi{mu}{nu} + IC*XX'*pure_to_mixed(MaxEntangled(dim{N-1}{alpha}))*XX*IC'*dim{N-1}{alpha}/(multout{N-1}{alpha});
        end
    end
end

F = 0;
err = 0;
prob = 0;
for index_mu = 1:length(list_young_diagrams_in{N})
    mu = list_young_diagrams_in{N}{index_mu};
    F = F + real(trace(S{mu}{mu}*Omega{mu}));
end

for index_mu = 1:length(list_young_diagrams_in{N})
    mu = list_young_diagrams_in{N}{index_mu};
    for index_nu = 1:length(list_young_diagrams_in{N})
        nu = list_young_diagrams_in{N}{index_nu};
        if trace(Sigma{mu}{nu}) > 0
            prob = prob + real(trace(S{mu}{nu}*Sigma{mu}{nu}))/din;
        end
        if task == 3 | task == 4
            if trace(Xi{mu}{nu}) > 0
                err = err + real(trace(S{mu}{nu}*Xi{mu}{nu}))/(dout-din);
            end
        end
    end
end

switch task
    case 1
        prob == F;
        maximize prob
        cvx_end
        cvx_status
        opt_fom = prob;
    case 2
        prob == 1;
        maximize F
        cvx_end
        cvx_status
        opt_fom = F;
    case 3
        prob == 1;
        minimize err
        cvx_end
        cvx_status
        opt_fom = err;
    case 4
        prob == 1;
        diamond_norm = max([err, 1-F]);
        minimize diamond_norm
        cvx_end
        cvx_status
        opt_fom = diamond_norm;
end