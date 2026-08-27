function [tangentialTorque, normalTorque] = getInteractionInducedTorque(model, state, desiredTorque)
    % interactionForce [Ft, Fn] unit magnitude and direction
    q1d = state.state(8);
    q2d = state.state(9);
    q3d = state.state(10);
    q4d = state.state(11);
    q5d = state.state(12);
    q6d = state.state(13);
    q7d = state.state(14);
    T01 = state.T01;
    T02 = T01 * state.T12;
    T03 = T02 * state.T23;
    T0_p1= T03 * state.T3_p1;
    T04 = T0_p1 * state.Tp1_4;
    T05 = T04 * state.T45;
    T0_p2 = T05 * state.T5_p2;
    T06 = T0_p2 * state.TP2_6;
    T07 = T06 * state.T67;
    T08 = T07 * state.T78;
    Jve = calculateJacobian(state.T01,T02,T03,T04,T05,T06,T07,T08,4);    
    
    % for wrist
    % for tangential motion
    ve = Jve *[q1d;q2d;q3d;q4d;q5d;q6d;q7d];
    z4 = T04(1:3,3)/norm(T04(1:3,3));
    % z4 and ve will never align with each other since the arm will always
    % move and there will be motion component perpendicular to z4.
    ntw = cross(ve,z4)/norm(cross(ve,z4));
    nnw = z4;
    
    % get direction
    arm.Jw = Jve;
    arm.Je = calculateJacobian(T01,T02,T03,T04,T05,T06,T07,4);
    [state_new, ~] = RKDiscretize(model.ts, state.state(1:14), desiredTorque, zeros(3,1), zeros(3,1), model, arm);
    q1n = state_new(1);
    q2n = state_new(2);
    q3n = state_new(3);
    q4n = state_new(4);
    q5n = state_new(5);
    q6n = state_new(6);
    q7n = state_new(7);
    
    T01n = DHConvention(-pi/2, 0, 0, -pi/2 + q1n); %
    T12n = DHConvention(-pi/2, 0, 0, -pi/2 + q2n); %
    T02n = T01n*T12n;
    T23n = DHConvention(pi/2, 0, 0, q3n); %
    T03n = T02n*T23n;
    T3_p1n= DHConvention(-pi/2, 0, 0, -pi/2);
    T0_p1n= T03n * T3_p1n;
    Tp1_4n = DHConvention(-pi/2, lu, 0, pi/2+q4n); %
    T04n = T0_p1n * Tp1_4n;
    T45n = DHConvention(-pi/2, 0, 0, q5n); 
    T05n = T04n * T45n;
    T5_p2n = DHConvention(-pi/2, 0, 0, pi/2); 
    T0_p2n = T05n * T5_p2n;
    Tp2_6n = DHConvention(0, lf, 0, q6n); 
    T06n = T0_p2n * Tp2_6n;
    T67n = DHConvention(-pi/2, 0, 0, q7n); 
    T07n = T06n * T67n;
    T78n= DHConvention(0, lw, 0, 0); 
    T08n=T07n*T78n;
    


   


    pos_wrist= T06n(1:3,4);


    pos_wrist_ref = T06(1:3,4);


    displacement_wrist = (pos_wrist_ref - pos_wrist)/norm(pos_wrist_ref - pos_wrist);

    if dot(displacement_wrist, ntw) > 1e-3
        coeff_wt = 1;
    elseif dot(displacement_wrist, ntw) < -1e-3
        coeff_wt = -1;
    else
        coeff_wt = 0;
    end
    
    if dot(displacement_wrist, nnw) > 1e-3
        coeff_wn = 1;
    elseif dot(displacement_wrist, nnw) < -1e-3
        coeff_wn = -1;
    else
        coeff_wn = 0;
    end
    
    % for wrist
    Fnt = coeff_wt*ntw;  % force in the plane
    tangentialTorque = Jve'*Fnt;
    Fnn = coeff_wn*nnw;  % force normal to the plane
    normalTorque = Jve'*Fnn;
    

end