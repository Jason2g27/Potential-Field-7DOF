function [y_next, x_state_next] = butter3_step(u, x_state, fc, dt)
% BUTTER3_STEP Executes one real-time step of a 3rd-order Butterworth filter.
%
% Inputs:
%   u            - Current raw input scalar (e.g., target velocity)
%   x_state      - Current [3x1] internal state vector from the previous step
%   fc           - Cutoff frequency in Hz
%   dt           - Simulation time step in seconds
%
% Outputs:
%   y_next       - Current smooth filtered output scalar
%   x_state_next - Updated [3x1] internal state vector for the next step

    % 1. Calculate discrete filter coefficients via bilinear transform
    Wn = 2 * fc * dt; % Normalized frequency (Nyquist = 1)
    [b, a] = butter(3, Wn);

    % 2. Convert Transfer Function coefficients to State-Space Form
    % This prevents messy history array shifts and improves numerical stability
    [A, B, C, D] = tf2ss(b, a);

    % 3. Calculate current output
    y_next = C * x_state + D * u;

    % 4. Compute next internal states for the next time step
    x_state_next = A * x_state + B * u;
end
