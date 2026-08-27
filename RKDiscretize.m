function [state_new,state_dot] = RKDiscretize(ts, state_old,state_dot_old, u, f_elbow, f_wrist, model, arm,M,C,G)

    k1 = forwardDynamics(state_old,state_dot_old, u,M,C,G);
    k2 = forwardDynamics(state_old + k1*ts/2,state_dot_old, u,M,C,G);
    k3 = forwardDynamics(state_old + k2*ts/2,state_dot_old, u,M,C,G);
    k4 = forwardDynamics(state_old + k3*ts,state_dot_old, u,M,C,G);

    state_new = state_old + ts/6 * (k1 + 2*k2 + 2*k3 + k4);
    state_dot = (k1 + 2*k2 + 2*k3 + k4)/6;
    
end