function [q_smooth, v_smooth, a_smooth] = quintic_spline_trajectory(q_recorded, dt)
    [MN,N] = size(q_recorded);
    t = (0:N-1)' * dt;
    
    % Estimate continuous velocity and acceleration from data using central differences
    v_rec = gradient(q_recorded) / dt;
    a_rec = gradient(v_rec) / dt;
    
    % Enforce boundary conditions to be perfectly resting (zero spikes at start/end)
    v_rec(:,1) = 0; v_rec(:,end) = 0;
    a_rec(:,1) = 0; a_rec(:,end) = 0;
    
    q_smooth = zeros(7, N);
    v_smooth = zeros(7, N);
    a_smooth = zeros(7, N);
    
    % Loop through each segment and fit a 5th-order (Quintic) polynomial
    for i = 1:N-1
        t0 = t(i);
        tf = t(i+1);
        T = tf - t0;
        
        % Segment boundary states
        q0 = q_recorded(:,i);   qf = q_recorded(:,i+1);
        v0 = v_rec(:,i);        vf = v_rec(:,i+1);
        a0 = a_rec(:,i);        af = a_rec(:,i+1);
        
        % Solve algebraic quintic coefficients for this specific segment
        % This forces boundary position, velocity, and acceleration to match perfectly.
        a0_coeff = q0;
        a1_coeff = v0;
        a2_coeff = a0 / 2;
        a3_coeff = (10*(qf - q0) - (4*vf + 6*v0)*T - (3*a0 - af)*(T^2)/2) / (T^3);
        a4_coeff = (-15*(qf - q0) + (7*vf + 8*v0)*T + (3*a0 - 2*af)*(T^2)/2) / (T^4);
        a5_coeff = (6*(qf - q0) - 3*(vf + v0)*T - (a0 - af)*(T^2)/2) / (T^5);
        
        % Evaluate the precise midpoint or segment steps (handling the dt interval)
        tau = 0; % Local time within segment
        q_smooth(:,i) = a0_coeff;
        v_smooth(:,i) = a1_coeff;
        a_smooth(:,i) = 2*a2_coeff;
    end
    
    % Lock the final endpoint exactly to your recorded point
    q_smooth(:,i) = q_recorded(end);
    v_smooth(:,i) = 0;
    a_smooth(:,i) = 0;
end
