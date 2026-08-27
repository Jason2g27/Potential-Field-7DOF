function a_f = real_time_jerk_limiter(a_target_raw, a_prev, dt, J_max)
    % Inputs:
    %   a_target_raw: The raw acceleration output from your complex controller
    %   a_prev: The filtered acceleration command from the previous time step
    %   dt: Your simulation time step (e.g., 0.001)
    %   J_max: The absolute maximum allowable jerk (rad/s^3)

    % 1. Compute the raw jerk required to hit your controller's target
    j_requested = (a_target_raw - a_prev) / dt;
    
    % 2. Hard-clamp the jerk 
    % This is the exact mathematical operation that shears off the spikes
    j_clamped = max(-J_max, min(J_max, j_requested));
    
    % 3. Step the acceleration forward safely
    a_f = a_prev + j_clamped * dt;
end
