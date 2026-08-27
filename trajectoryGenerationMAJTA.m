function [trajectory, true_delay] = trajectoryGenerationMAJTA(arm_length, start_pos, end_pos, t_tot, ts)
    syms s d real;
    lf = arm_length(1);
    lu = arm_length(2);
    theta = sym('theta', [1, length(start_pos)]);

    if end_pos(end)-start_pos(end) == 0
        delay = 0;
    else
        theta(1:end-1) = start_pos(1:end-1) + (start_pos(1:end-1) - end_pos(1:end-1))*(15*s^4 - 6*s^5- 10*s^3);
        theta_elbow_s2 = start_pos(end) + (start_pos(end) - end_pos(end))*(15*((s-d)/(1-d))^4 - 6*((s-d)/(1-d))^5- 10*((s-d)/(1-d))^3);
        theta(end) = theta_elbow_s2;
        % forward kinematics
        T01 = DHConvention(-pi/2, 0, 0, -pi/2 + theta(1)); %
        T12 = DHConvention(-pi/2, 0, 0, -pi/2 + theta(2)); %
        T02 = T01*T12;
        T23 = DHConvention(pi/2, 0, 0, theta(3)); %
        T03 = T02*T23;
        T3_p1= DHConvention(-pi/2, 0, 0, -pi/2);
        T0_p1= T03 * T3_p1;
        Tp1_4 = DHConvention(-pi/2, lu, 0, pi/2+theta(4)); %
        T04 = T0_p1 * Tp1_4;
        T45 = DHConvention(-pi/2, 0, 0, theta(5)); 
        T05 = T04 * T45;
        T5_p2 = DHConvention(-pi/2, 0, 0, pi/2); 
        T0_p2 = T05 * T5_p2;
        Tp2_6 = DHConvention(0, lf, 0, theta(6)); 
        T06 = T0_p2 * Tp2_6;
        T67 = DHConvention(-pi/2, 0, 0, theta(7)); 
        T07 = T06 * T67;
        T78= DHConvention(0, lw, 0, 0); 
        T08=T07*T78;
        x = T05(1,4);
        y = T05(2,4);
        z = T05(3,4);
        vel  = norm([diff(x,s)^2,  diff(y,s)^2,  diff(z,s)^2]);
        acceleration_tip = diff(vel,s);
        acceleration_tip = subs(acceleration_tip, s, 0.5);
        equation = norm(acceleration_tip) == 0;
        solution = vpasolve(equation, d, [0, 1]);
        
        % % The Secant Method
        % x_n = 0.5;
        % x_nm1 = 0;
        % i = 0;
        % equation = subs(acceleration_tip, s, 0.5);
        % equation_with_value = norm(equation);
        % while 1
        %     disp("here")
        %     Q = (subs(equation_with_value, d, x_nm1)-subs(equation_with_value, d, x_n))/(x_nm1-x_n);
        %     disp("here")
        %     x_np1 = x_n - subs(equation_with_value, d, x_n)/Q;
        %     if abs(x_np1-x_n) <= 1e-3 && i <= 1000
        %         solution = x_np1;
        %         break
        %     elseif i > 1000
        %         solution = 0;
        %         break
        %     else
        %         x_n = x_np1;
        %         x_nm1 = x_n;
        %         i = i + 1;
        %     end
        % end


        if solution < 0.5 && solution > 0
            delay = double(solution);
        else 
            disp("no feasible delay")
            delay = 0;
        end
    end
    time = 0:ts:t_tot;
    ss = time/t_tot;
    theta = zeros(7,length(ss));
    theta_dot = zeros(7,length(ss));
    theta_ddot = zeros(7,length(ss));
    theta(1:6,:) = start_pos(1:6) + (start_pos(1:6) - end_pos(1:6)).*(15*ss.^4 - 6*ss.^5- 10*ss.^3);
    theta_dot(1:6,:) = (start_pos(1:6) - end_pos(1:6)).*(60*ss.^3*1/t_tot - 30*ss.^4*1/t_tot- 30*ss.^2*1/t_tot);
    theta_ddot(1:6,:) = (start_pos(1:6) - end_pos(1:6)).*(180*ss.^2*1/t_tot^2 - 120*ss.^3*1/t_tot^2- 60*ss*1/t_tot^2);

    if delay == 0
        theta(7,:) = start_pos(7) + (start_pos(7) - end_pos(7)).*(15*ss.^4 - 6*ss.^5- 10*ss.^3);
        theta_dot(7,:) = (start_pos(7) - end_pos(7)).*(60*ss.^3*1/t_tot - 30*ss.^4*1/t_tot- 30*ss.^2*1/t_tot);
        theta_ddot(7,:) = (start_pos(7) - end_pos(7)).*(180*ss.^2*1/t_tot^2 - 120*ss.^3*1/t_tot^2- 60*ss*1/t_tot^2);
    else
        t_s1 = floor(delay*t_tot/ts)+1;
        theta(7, :) = start_pos(end) + (start_pos(end) - end_pos(end)).*(15*((ss-delay)/(1-delay)).^4 - 6*((ss-delay)/(1-delay)).^5- 10*((ss-delay)/(1-delay)).^3);
        theta_dot(7, :) = (start_pos(end) - end_pos(end)).*(60*((ss-delay)/(1-delay)).^3*(1/t_tot)/(1-delay) - 30*((ss-delay)/(1-delay)).^4*(1/t_tot)/(1-delay)- 30*((ss-delay)/(1-delay)).^2*(1/t_tot)/(1-delay));
        theta_ddot(7, :) = (start_pos(end) - end_pos(end)).*(180*((ss-delay)/(1-delay)).^2*(1/t_tot)^2/(1-delay)^2 - 120*((ss-delay)/(1-delay)).^3*(1/t_tot)^2/(1-delay)^2- 60*((ss-delay)/(1-delay))*(1/t_tot)^2/(1-delay)^2);
        theta(7,1:t_s1) = 0;
        theta_dot(7,1:t_s1) = 0;
        theta_ddot(7,1:t_s1) = 0;
    end
    trajectory = [theta;theta_dot;theta_ddot];
    true_delay = delay*t_tot;
    disp(['Delay is: ', num2str(true_delay), ' sec.']);
end
