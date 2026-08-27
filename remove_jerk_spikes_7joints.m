function [q_clean, v_smooth, a_smooth, j_smooth] = remove_jerk_spikes_7joints(q_recorded, dt)
    % q_recorded is assumed to be 7 rows (joints) x N columns (samples)
    [num_joints, N] = size(q_recorded);
    t = (0:N-1) * dt; % Time row vector

    % 1. Calculate raw derivatives (Single output gradient works across columns perfectly)
    v_raw = gradient(q_recorded) / dt;
    a_raw = gradient(v_raw) / dt;
    j_raw = gradient(a_raw) / dt; 

    % Allocate output matrices matching the 7 x N shape
    q_clean = q_recorded;

    % 2. Process each of the 7 joints individually
    for joint = 1:num_joints
        % Extract raw data for this specific joint (1 x N vector)
        jerk_single = j_raw(joint, :);
        q_single = q_recorded(joint, :);
        
        % Detect outliers on this specific joint's jerk profile
        [~, outlier_indices] = filloutliers(jerk_single, 'linear', 'movmedian', 15, 'ThresholdFactor', 3);
        
        % Bridge the gaps only if outliers were actually found
        if any(outlier_indices)
            good_indices = ~outlier_indices;
            % Note: transposing t and q_single inside ensures interp1 is happy
            q_clean(joint, outlier_indices) = interp1(t(good_indices), q_single(good_indices), t(outlier_indices), 'spline');
        end
    end

    % 3. Re-calculate the final clean profiles for all joints simultaneously
    v_smooth = gradient(q_clean) / dt;
    a_smooth = gradient(v_smooth) / dt;
    j_smooth = gradient(a_smooth) / dt;
end

