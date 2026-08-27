function [N2, B2] = transportFrame(P1, T1, N1, P2, T2)
    % Inputs:
    % P1, P2 : 1x3 positions
    % T1, N1 : 1x3 previous tangent and normal
    % T2     : 1x3 target tangent
    
    % Step 1: First Reflection (R1)
    v1 = P2 - P1;
    c1 = dot(v1,v1);
    
    if c1 < 1e-15
        % P1 and P2 are essentially the same point
        TL = T1;
        NL = N1;
    else
        TL = T1 - (2 / c1) * (v1 * T1') * v1;
        NL = N1 - (2 / c1) * (v1 * N1') * v1;
    end
    
    % Step 2: Second Reflection (R2)
    v2 = T2 - TL;
    c2 = dot(v2,v2);
    
    if c2 < 1e-15
        % Tangent hasn't changed relative to TL
        N2 = NL;
    else
        N2 = NL - (2 / c2) * (v2 * NL') * v2;
    end
    
    % Step 3: Compute final vector
    B2 = cross(T2, N2);
    
    % Optional: Normalize to prevent numerical drift
    T2 = T2 / norm(T2);
    N2 = N2 / norm(N2);
    B2 = B2 / norm(B2);
end