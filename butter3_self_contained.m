function [y_next, hist_next] = butter3_self_contained(u, hist, fc, dt)
% BUTTER3_SELF_CONTAINED Filter step that calculates coefficients internally.
%
% Inputs:
%   u         - Current raw input scalar
%   hist      - 6x1 history vector
%   fc        - Cutoff frequency in Hz
%   dt        - Time step in seconds

    % Define persistent variables to cache the coefficients across loop calls
    persistent b a current_fc current_dt
    
    % If this is the first run, or if fc/dt changed, compute the coefficients
    if isempty(b) || isempty(a) || (fc ~= current_fc) || (dt ~= current_dt)
        Wn = 2 * fc * dt;
        [b_raw, a_raw] = butter(3, Wn);
        
        % Normalize by a0 immediately to save execution speed
        b = b_raw / a_raw(1);
        a = a_raw / a_raw(1);
        
        % Save current settings to detect changes later
        current_fc = fc;
        current_dt = dt;
    end

    % 1. Extract historical data points
    x1 = hist(1); x2 = hist(2); x3 = hist(3);
    y1 = hist(4); y2 = hist(5); y3 = hist(6);

    % 2. Execute raw scalar arithmetic
    y_next = b(1)*u + b(2)*x1 + b(3)*x2 + b(4)*x3 ...
                    - a(2)*y1 - a(3)*y2 - a(4)*y3;

    % 3. Slide the history window forward
    hist_next = [u; x1; x2; y_next; y1; y2];
end
