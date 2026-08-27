function [q_next, v_next] = rkn4_step(q, v, dt, calc_acceleration)
    % RKN4_STEP Advances a second-order system by one time step dt
    % Using the 4th-order Runge-Kutta-Nyström method.
    %
    % Inputs:
    %   q                 - Column vector of current positions (DoF x 1)
    %   v                 - Column vector of current velocities (DoF x 1)
    %   dt                - Time step size (scalar)
    %   calc_acceleration - Function handle to a(q, v) returning (DoF x 1)
    
    % --- Stage 1 ---
    k1 = calc_acceleration(q, v);
    
    % --- Stage 2 ---
    q2 = q + 0.5 * dt * v + 0.125 * dt^2 * k1;
    v2 = v + 0.5 * dt * k1;
    k2 = calc_acceleration(q2, v2);
    
    % --- Stage 3 ---
    q3 = q + 0.5 * dt * v + 0.125 * dt^2 * k1;
    v3 = v + 0.5 * dt * k2;
    k3 = calc_acceleration(q3, v3);
    
    % --- Stage 4 ---
    q4 = q + dt * v + 0.5 * dt^2 * k3;
    v4 = v + dt * k3;
    k4 = calc_acceleration(q4, v4);
    
    % --- Final Weight Assembly ---
    q_next = q + dt * v + (dt^2 / 6.0) * (k1 + k2 + k3);
    v_next = v + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
end