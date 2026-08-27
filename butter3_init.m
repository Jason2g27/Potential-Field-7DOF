function [b, a] = butter3_init(fc, dt)
% BUTTER3_INIT Pre-calculates and normalizes 3rd-order Butterworth coefficients.
%
% Inputs:
%   fc - Desired cutoff frequency in Hz
%   dt - Simulation time step in seconds
%
% Outputs:
%   b  - Normalized feedforward coefficients [b0, b1, b2, b3]
%   a  - Normalized feedback coefficients [a0, a1, a2, a3]

    % 1. Calculate normalized frequency (Nyquist frequency = 1)
    Wn = 2 * fc * dt; 
    
    % 2. Generate raw 3rd-order Butterworth filter coefficients
    [b_raw, a_raw] = butter(3, Wn); 
    
    % 3. Normalize by a0 immediately so the runtime loop avoids costly division
    b = b_raw / a_raw(1); 
    a = a_raw / a_raw(1); 
end
