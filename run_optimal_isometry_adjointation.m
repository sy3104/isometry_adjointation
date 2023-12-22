% -------------------------------------------------------------------------
%
% File : run_optimal_isometry_adjointation.m
%
% Discription : 
% Code to obtain the optimal approximation error or success probability of
% transforming 'n' calls of any isometry operation with input dimension
% 'din' and output dimension 'dout' into its inverse map, corresponding
% POVM, or adjoint operation.
% Please set the parameters 'din', 'dout', 'n', 'protocol', 'task',
% 'isComplex', and 'isDual'.
%
% -------------------------------------------------------------------------

clear

% ------------------------------------------------------------------
%                   Start of setting the parameters
% ------------------------------------------------------------------

din = 2                % Input dimension of the input isometry operation
dout = 3               % Output dimension of the input isometry operation
n = 1                  % Number of calls
protocol = 1           % 1 for parallel, 2 for sequential, 3 for indefinite
task = 1               % See below
isComplex = 0          % Set 0 for real Choi matrix
isDual = 0             % Set 0 for solve primal SDP (except for indefinite)

% Correspondence of 'task' and obtained optimal value
% 1: Probabilistic exact isometry inversion (maximal probability)
% 2: Deterministic isometry inversion (maximal fidelity)
% 3: Universal error detection (minimal one-sided error)
% 4: Isometry adjointation (minimal diamond-norm distance)

if task == 1 | task == 2
    if dout < din
        error('dout should be greater than or equal to din!')
    end
elseif task == 3 | task == 4
    if dout <= din
        error('dout should be greater than din!')
    end
else
    error('task should be 1, 2, 3 or 4!')
end

if isDual == 0 & protocol == 3
    error('Primal problem for indefinite causal order is not implemented. Please set isDual = 1!')
end

% ------------------------------------------------------------------
%                   End of setting the parameters
% ------------------------------------------------------------------

tic;

switch protocol
    case 1
        if isDual == 0
            opt_fom = parallel_isometry_adjointation(din,dout,n,task,isComplex);
        else
            opt_fom = parallel_isometry_adjointation_dual(din,dout,n,task,isComplex);
        end
    case 2
        if isDual == 0
            opt_fom = sequential_isometry_adjointation(din,dout,n,task,isComplex);
        else
            opt_fom = sequential_isometry_adjointation_dual(din,dout,n,task,isComplex);
        end
    case 3
        if n == 1
            opt_fom = parallel_isometry_adjointation_dual(din,dout,n,task,isComplex);
        else
            opt_fom = indefinite_isometry_adjointation_dual(din,dout,n,task,isComplex);
        end
end

optimal_task_of_merit = opt_fom

total_time_in_seconds = toc
