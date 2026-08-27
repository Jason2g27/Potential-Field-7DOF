function [x, result, status] = gurobi_lsqlin(C, d, A, b, Aeq, beq, lb, ub)
% GUROBI_LSQLIN Solves a constrained linear least-squares problem using Gurobi.
%
%   Matches the mathematical formulation of MATLAB's lsqlin:
%   min_x  ||C*x - d||^2   subject to:  A*x <= b,  Aeq*x == beq,  lb <= x <= ub
%
%   Outputs:
%     x      - Optimized decision variable vector
%     result - The full raw output structure from Gurobi
%     status - Solution status string (e.g., 'OPTIMAL')

    % --- 1. Initialize Model Structure ---
    model2 = struct();

    % --- 2. Convert Least-Squares to Quadratic Form ---
    % Objective: ||C*x - d||^2 = x'*(C'*C)*x - 2*(C'*d)*x + d'*d
    % Gurobi minimizes: 1/2 * x'*Q*x + obj'*x
    % Therefore: Q = 2 * C'*C  and  obj = -2 * C'*d
    model2.Q = sparse(C' * C); 
    linear_objective = -(C' * d);
    model2.obj = linear_objective(:); 
    
    epsilon = 1e-4; 
    model2.Q = model2.Q + epsilon * speye(size(model2.Q, 1));
    % --- 3. Stack and Align Constraints Vertically ---
    % Check which constraints were provided by the user
    has_ineq = (~isempty(A) && ~isempty(b));
    has_eq   = (~isempty(Aeq) && ~isempty(beq));

    if has_ineq && has_eq
        model2.A     = sparse([A; Aeq]);
        model2.rhs   = full([b(:); beq(:)]);
        model2.sense = [repmat('<', size(A, 1), 1); repmat('=', size(Aeq, 1), 1)];
    elseif has_ineq
        model2.A     = sparse(A);
        model2.rhs   = full(b(:));
        model2.sense = repmat('<', size(A, 1), 1);
    elseif has_eq
        model2.A     = sparse(Aeq);
        model2.rhs   = full(beq(:));
        model2.sense = repmat('=', size(Aeq, 1), 1);
    else
        % Unconstrained or bounded-only problem
        model2.A = sparse(0, size(C, 2));
        model2.rhs = [];
        model2.sense = [];
    end

    % --- 4. Handle Variable Bounds ---
    % If bounds are omitted or empty, Gurobi defaults lower bounds to 0.
    % We explicitly override this to match lsqlin behavior (-inf to inf).
    if nargin >= 7 && ~isempty(lb)
        model2.lb = lb(:);
    else
        model2.lb = -inf(size(C, 2), 1);
    end

    if nargin >= 8 && ~isempty(ub)
        model2.ub = ub(:);
    else
        model2.ub = inf(size(C, 2), 1);
    end
    % --- Explicitly Force Continuous Variables ---
    num_vars = size(C, 2);
    model2.vtype = repmat('C', num_vars, 1);
    
    % --- 5. Configure Solver Performance Parameters ---
    params = struct();
    params.Precheck = 1; 
    params.Presolve = 0;
    %params.TimeLimit = 5;       % Safety Net
    params.Threads = 4;         % System Stability
    
    params.outputflag = 1;   % Keep this at 1 for now to monitor progress
    
    % FIX 1: Tell Gurobi to treat the problem as convex even if matrix scaling is messy
    params.NonConvex = 1;    
    
    % FIX 2: Force Gurobi to use its fast Barrier (Interior-Point) algorithm
    params.Method = 2;       
    
    % FIX 3: Turn on aggressive scaling to handle the 1e-13 numbers
    params.ScaleFlag = 2;    
    
    % FIX 4: Add a numeric stability protector to prevent the solver from getting stuck
    params.NumericFocus = 3; 

    % --- 6. Execute Gurobi Optimization ---
    try
        result = gurobi(model2, params);
        status = result.status;
        
        % MODIFIED CODE: Accept both OPTIMAL and SUBOPTIMAL results
        if strcmp(status, 'OPTIMAL') || strcmp(status, 'SUBOPTIMAL')
            if isfield(result, 'x')
                x = result.x; 
            else
                x = [];
            end
        else
            x = [];
            warning('Gurobi finished with catastrophic status: %s', status);
        end
    catch ME
        fprintf('Error executing Gurobi solver.\n');
        rethrow(ME);
    end 
    % 1. Dynamically calculate the number of points
    %{
    num_coordinates = dimension; 
    num_points = length(x) / num_coordinates;

    % 2. Instantly reshape it into a 3-row matrix (3 by N)
    optimized_matrix = reshape(x, num_coordinates, num_points);
    x2=optimized_matrix;
    %}
end