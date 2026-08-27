function [u_filtered, hist_next] = butter3_exo_step(u_raw, hist, b, a)
% BUTTER3_EXO_STEP Smooths a 7-DOF raw joint acceleration command vector.
%
% Inputs:
%   u_raw     - [7x1] current raw joint acceleration vector
%   hist      - [6x7] current history matrix tracking all 7 joints
%   b, a      - Pre-calculated and normalized filter coefficient vectors
%
% Outputs:
%   u_filtered - [7x1] smooth joint acceleration vector for your M matrix
%   hist_next  - [6x7] updated history matrix for the next time step

    % 1. Extract historical states for all 7 joints simultaneously
    x1 = hist(:,1); x2 = hist(:,2); x3 = hist(:,3);
    y1 = hist(:,4); y2 = hist(:,5); y3 = hist(:,6);

    % 2. Execute parallel vector arithmetic across all 7 joints
    % Transposing u_raw to a row vector [1x7] ensures matching dimensions
    %u_raw_row = u_raw';
    y_next = b(1)*u_raw + b(2)*x1 + b(3)*x2 + b(4)*x3 ...
                            - a(2)*y1 - a(3)*y2 - a(4)*y3;

    % 3. Package the updated rows back into the [6x7] matrix
    hist_next = [u_raw, x1, x2, y_next, y1, y2];

    % 4. Convert output back to a [7x1] column vector for your robot model
    u_filtered = y_next;
end
