function state_dot = forwardDynamics(state,state_dot, torque,M,C,G)
   % Jw = arm.Jw;
    %Je = arm.Je;
    %state_dot = zeros (14,1);
    % [~, S, ~] = svd(M);
    % singular = min(diag(S));
    % if singular < 0.001
    M = M + 0.01*eye(size(M)); % regularize the mass to avoid near singularity

    % end
    state_dot(1:7) = state(8:14);
    state_dot(8:14) = M\(torque-C-G); %+ inv(M)*Jw'*f_wrist + inv(M)*[Je zeros(7,1)]'*f_elbow;
    %state_dot = state_dot';
    
end