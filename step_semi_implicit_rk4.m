function [q_next, v_next,acc] = step_semi_implicit_rk4(q, v, dt, forward_dynamics_fn)
    % STEP_SEMI_IMPLICIT_RK4 Executes one 4th-order symplectic integration step.
    %
    % Inputs:
    %   q                   - Current 7x1 joint position vector
    %   v                   - Current 7x1 joint velocity vector
    %   dt                  - Time step size (e.g., 0.001 for 1 kHz)
    %   forward_dynamics_fn - Handle to your function: acc = forward_dynamics(q, v)
    %
    % Outputs:
    %   q_next              - Stable 4th-order position vector at t + dt
    %   v_next              - Stable 4th-order velocity vector at t + dt

    % Forest-Ruth analytical coefficients for symplectic 4th-order accuracy
    w0 = -2^(1/3) / (2 - 2^(1/3));
    w1 = 1 / (2 - 2^(1/3));
    
    c = [w1/2, (w0+w1)/2, (w0+w1)/2, w1/2]; % Position steps
    d = [w1, w0, w1, 0];                    % Velocity steps
    
    % Initialize scratchpad variables
    q_temp = q;
    v_temp = v;
    
    % Stage 1
    q_temp = q_temp + c(1) * dt * v_temp;
    acc = forward_dynamics_fn(q_temp, v_temp);
    v_temp = v_temp + d(1) * dt * acc;
    
    % Stage 2
    q_temp = q_temp + c(2) * dt * v_temp;
    acc = forward_dynamics_fn(q_temp, v_temp);
    v_temp = v_temp + d(2) * dt * acc;
    
    % Stage 3
    q_temp = q_temp + c(3) * dt * v_temp;
    acc = forward_dynamics_fn(q_temp, v_temp);
    v_temp = v_temp + d(3) * dt * acc;
    
    % Stage 4 (Final update boundary)
    q_next = q_temp + c(4) * dt * v_temp;
    v_next = v_temp + d(4) * dt * acc; % d(4) is 0, evaluated for strict structure
end
