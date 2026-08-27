function [b, a] = bessel3_init(fc, dt)
    % 1. Design a continuous-time 3rd-order Bessel filter
    w_rad = 2 * pi * fc;
    [num, den] = besself(3, w_rad);
    
    % 2. Convert to discrete-time via Bilinear Transform (matches your dt)
    [b_raw, a_raw] = bilinear(num, den, 1/dt);
    
    % 3. Normalize by a0 for your high-speed loop
    b = b_raw / a_raw(1);
    a = a_raw / a_raw(1);
end