function [Ni, Bi] = nextBishopFrame(Pi_prev, Ti_prev, Ni_prev, Pi, Ti)
    % Inputs:
    % Pi_prev, Pi : 1x3 positions (Previous and Current)
    % Ti_prev, Ti : 1x3 tangents (Previous and Current)
    % Ni_prev     : 1x3 normal vector (Previous)

    % 1. Reflection 1: Reflect onto the plane bisecting the segment P_prev to P
    v1 = Pi - Pi_prev;
    c1 = dot(v1,v1);
    
    if c1 < 1e-15
        Ni_mid = Ni_prev;
        Ti_mid = Ti_prev;
    else
        % Reflect both Tangent and Normal to an intermediate state
        Ti_mid = Ti_prev - (2/c1) * (v1 * Ti_prev') * v1;
        Ni_mid = Ni_prev - (2/c1) * (v1 * Ni_prev') * v1;
    end

    % 2. Reflection 2: Reflect from intermediate state to the current Tangent (Ti)
    v2 = Ti - Ti_mid;
    c2 = dot(v2,v2);
    
    
    if c2 < 1e-15
        Ni = Ni_mid;
    else
        % This second reflection aligns Ti_mid with Ti and adjusts Ni accordingly
        Ni = Ni_mid - (2/c2) * (v2 * Ni_mid') * v2;
    end

    % 3. Final Orthonormal Frame
    Bi = cross(Ti, Ni);
    Bi = Bi / norm(Bi); 
    Ni = cross(Bi, Ti); % Ensure Ni is perfectly orthogonal to Ti
end