function [y_next, hist_next] = butter3_traditional_step(u, hist, b, a)
% BUTTER3_TRADITIONAL_STEP Executes one real-time step of a 3rd-order 
% Butterworth filter using raw scalar arithmetic.
%
% Inputs:
%   u         - Current raw input scalar (e.g., target velocity)
%   hist      - 6x1 history vector: [x(n-1); x(n-2); x(n-3); y(n-1); y(n-2); y(n-3)]
%   b         - Pre-calculated feedforward coefficients [b0, b1, b2, b3]
%   a         - Pre-calculated feedback coefficients [a0, a1, a2, a3]
%
% Outputs:
%   y_next    - Current smooth filtered output scalar
%   hist_next - Updated 6x1 history vector for the next time step

    % 1. Extract historical data points from the vector
    x1 = hist(1); x2 = hist(2); x3 = hist(3);
    y1 = hist(4); y2 = hist(5); y3 = hist(6);

    % 2. Execute raw scalar arithmetic (the fastest possible execution format)
    % Note: Coefficients must be pre-normalized by a0 outside the loop
    y_next = b(1)*u + b(2)*x1 + b(3)*x2 + b(4)*x3 ...
                    - a(2)*y1 - a(3)*y2 - a(4)*y3;

    % 3. Slide the history window forward by one step to prepare for the next loop
    hist_next = [u; x1; x2; y_next; y1; y2];
end
