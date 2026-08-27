function [q_ref, v_ref, a_ref] = smoothReferenceFilter(q_cmd, dt, omega_n)
% smoothReferenceFilter Smooths raw cubic spline commands to eliminate acceleration spikes.
%
% Inputs:
%   q_cmd   - Vector of raw cubic trajectory position commands (1 x N)
%   dt      - Sampling time period in seconds (scalar, e.g., 0.001)
%   omega_n - Cutoff frequency in rad/s (scalar, sweet spot: 30 to 150)
%
% Outputs:
%   q_ref   - Spikeless reference positions (1 x N)
%   v_ref   - Spikeless reference velocities (1 x N)
%   a_ref   - Spikeless reference accelerations (1 x N)

    % Get the number of trajectory points
    [MN,N] = size(q_cmd);
    
    % Pre-allocate output arrays for speed
    q_ref = zeros(7, N);
    v_ref = zeros(7, N);
    a_ref = zeros(7, N);
    
    % Initialize the filter state with the first command value
    q_ref(:,1) = q_cmd(:,1);
    v_ref(:,1) = 0.0; % Assuming starting from rest
    a_ref(:,1) = 0.0;
    
    % Loop through the trajectory step-by-step
    for k = 2:N
        % 1. Calculate tracking error relative to the PREVIOUS step
        e = q_cmd(:,k) - q_ref(:,k-1);
        
        % 2. Compute smooth reference acceleration for the CURRENT step
        a_ref(:,k) = (omega_n^2 * e) - (2 * omega_n * v_ref(:,k-1));
        
        % 3. Integrate to find current reference velocity
        v_ref(:,k) = v_ref(:,k-1) + a_ref(:,k) * dt;
        
        % 4. Integrate to find current reference position
        q_ref(:,k) = q_ref(:,k-1) + v_ref(:,k) * dt;
    end
end
