function [Torque_Wrist, U] = armSimulator(ss, t_traj, sim, trajectory, impairment, torquelimit, model, controller_type, plot_freq, debug)
    % impairment  2 - mild, 1 - severe, 0 - not care
    fs = 16;  % font size
    % initialize 
    t = sim.t_tot_sim;
    ts = sim.ts;
    t_now = 0;
    t_now_traj = 0;
    referencePos = getTrajectoryCartesian(model, trajectory);
    state = zeros(8,1);
    state(1:4) = trajectory(1:4,1);
    state_dot = trajectory(5:8,1);
    Contact_Elbow = [];
    Contact_Wrist = [];
    Control_Effort = [];
    Torque_Elbow = [];
    Torque_Wrist = [];
    Hand_Vel = [];
    Hand_Acc = [];
    hand_vel_old = 0;
    Endpoint_elbow = [];
    Endpoint_wrist = [];
    U = [];
    Time = [];
    r = model.r;
    k = model.k;
    c = model.c;
    Kd = model.Kd;
    Kp = model.Kp;
    DesiredTorque = [];
    % ind_t = [1, 2:plot_freq:length(ss)-1, length(ss)];
    % ss = ss(ind_t);
    % trajectory = trajectory(:,ind_t);
    for jj = 1:length(trajectory(1,:))
        M = getM(model, trajectory(:,jj));
        C = getC(model, trajectory(:,jj));
        G = getG(model, trajectory(:,jj));
        DesiredTorque = [DesiredTorque inverseDynamics(trajectory(:,jj), M, C, G)];
    end

    if debug
        arm = kinematicAnalysis(model, state);
        fig_simu = figure();
        figure(fig_simu);
        set(fig_simu, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); 
        pause(0.3);
        subplot(4,2,[1,3,5]);
        hold on;
        axis equal;
        view(135, 30); 
        xlabel("X [m]", 'FontSize', fs);
        ylabel("Y [m]", 'FontSize', fs);
        zlabel("Z [m]", 'FontSize', fs);
        % axis([-0.7 0.7 -0.1 0.7 -0.7 0.7]); 
        title_time = title(['Time: ', num2str(t_now), ' sec'], 'FontSize', fs); 
        createTunnel(fig_simu, referencePos.Pos_elbow_ref', "#9EC9E2", "Virtual tunnel for elbow", r, 40);
        createTunnel(fig_simu, referencePos.Pos_wrist_ref', "#9CCEA7", "Virtual tunnel for wrist", r, 40);
        config = [zeros(3,1),arm.pos_elbow', arm.pos_wrist'];
        armconfig = plot3(config(1,:),config(2,:),config(3,:),"-", "color", '#00CD6C', 'LineWidth',5, 'DisplayName',"Arm configuration");
        % set(get(get(armconfig, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off'); 

        ref_traj_elbow = plot3(referencePos.Pos_elbow_ref(:,1),referencePos.Pos_elbow_ref(:,2),referencePos.Pos_elbow_ref(:,3),":", "color", '#FFC61E', 'LineWidth',3, 'DisplayName',"Reference elbow trajectory");
        ref_traj_wrist = plot3(referencePos.Pos_wrist_ref(:,1),referencePos.Pos_wrist_ref(:,2),referencePos.Pos_wrist_ref(:,3),":", "color", '#009ADE', 'LineWidth',3, 'DisplayName',"Reference wrist trajectory");
        pause(0.5);
        ref_traj_elbow.Color(4) = 0.7; 
        ref_traj_wrist.Color(4) = 0.7;

        endpoint_elbow = plot3(arm.pos_elbow(1),arm.pos_elbow(2),arm.pos_elbow(3), "-", "color", '#F28522', 'LineWidth',3, 'DisplayName',"Actual elbow trajectory");
        endpoint_wrist = plot3(arm.pos_wrist(1),arm.pos_wrist(2),arm.pos_wrist(3), "-", "color", '#AF58BA', 'LineWidth',3, 'DisplayName',"Actual wrist trajectory");
        pause(0.5);
        endpoint_elbow.Color(4) = 0.7; 
        endpoint_wrist.Color(4) = 0.7;

        arr_elbow = plotContact(arm.pos_elbow, 0, zeros(3,1),"#8F0038","Contact force on elbow");
        arr_wrist = plotContact(arm.pos_wrist, 0, zeros(3,1),"#B10026","Contact force on wrist");

        legend('FontSize', fs, 'Location', 'northoutside', 'NumColumns',3 , 'AutoUpdate', 'off'); 
        hold off;
        
        subplot(4,2,7);
        hold on;
        yyaxis left;
        velplot = plot(0,0,"-", "color","#009ADE", 'LineWidth',2, 'DisplayName',"Wrist speed");
        ylabel("Wrist speed [m/s]", 'FontSize', fs, 'Color', [0 0 0]);
        set(gca, 'XColor', 'k', 'YColor', 'k'); 
        yyaxis right;
        accplot = plot(0,0,"-", "color","#FF1F5B", 'LineWidth',2, 'DisplayName',"Wrist acceleration");
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("Wrist acceleration [m/s^2]", 'FontSize', fs, 'Color',[0 0 0]);
        % xlim([0, t]); 
        grid on;
        hold off;
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside', 'AutoUpdate', 'off'); 
        set(gca, 'XColor', 'k', 'YColor', 'k'); 


        Ginit = getG(model, state);
        subplot(4,2,2);
        hold on;
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("Torque [Nm]", 'FontSize', fs);
        % xlim([0, t]); 

        uplot11 = plot(0, Ginit(1),"-", 'color', "#FF1F5B", 'LineWidth',2, 'DisplayName',"Human applied torque");
        uplot12 = plot(0, 0, "-", 'color', "#00CD6C", 'LineWidth',2, 'DisplayName',"Elbow contact induced torque");
        uplot13 = plot(0, 0, "-", 'color', "#009ADE", 'LineWidth',2, 'DisplayName',"Wrist contact induced torque");
        uplot14 = plot(ss, DesiredTorque(1,:), "-", 'color', "#AF58BA", 'LineWidth',2, 'DisplayName',"Desired torque");
        title("Joint 1", FontSize=fs)
        lgd = legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside','NumColumns',2 , 'AutoUpdate', 'off'); 
        grid on;
        hold off;

        subplot(4,2,4);
        hold on;
        hold on;
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("Torque [Nm]", 'FontSize', fs);
        % xlim([0, t]); 

        uplot21 = plot(0, Ginit(2),"-", 'color', "#FF1F5B", 'LineWidth',2, 'DisplayName',"Human applied torque");
        uplot22 = plot(0, 0, "-", 'color', "#00CD6C", 'LineWidth',2, 'DisplayName',"Elbow contact induced torque");
        uplot23 = plot(0, 0, "-", 'color', "#009ADE", 'LineWidth',2, 'DisplayName',"Wrist contact induced torque");
        uplot24 = plot(ss, DesiredTorque(2,:), "-", 'color', "#AF58BA", 'LineWidth',2,'DisplayName',"Desired torque");
        title("Joint 2", FontSize=fs)
        grid on;
        hold off;

        subplot(4,2,6);
        hold on;
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("Torque [Nm]", 'FontSize', fs);
        % xlim([0, t]); 

        uplot31 = plot(0, Ginit(3),"-", 'color', "#FF1F5B", 'LineWidth',2, 'DisplayName',"Human applied torque");
        uplot32 = plot(0, 0, "-", 'color', "#00CD6C", 'LineWidth',2, 'DisplayName',"Elbow contact induced torque");
        uplot33 = plot(0, 0, "-", 'color', "#009ADE", 'LineWidth',2, 'DisplayName',"Wrist contact induced torque");
        uplot34 = plot(ss, DesiredTorque(3,:), "-", 'color', "#AF58BA", 'LineWidth',2, 'DisplayName',"Desired torque");
        title("Joint 3", FontSize=fs)
        grid on;
        hold off;

        subplot(4,2,8);
        hold on;
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("Torque [Nm]", 'FontSize', fs);
        % xlim([0, t]); 

        uplot41 = plot(0, Ginit(4),"-", 'color', "#FF1F5B", 'LineWidth',2, 'DisplayName',"Human applied torque");
        % uplot42 = plot(0, 0, "-", 'color', "#00CD6C", 'LineWidth',2, 'DisplayName',"Elbow contact induced torque");
        uplot43 = plot(0, 0, "-", 'color', "#009ADE", 'LineWidth',2, 'DisplayName',"Wrist contact induced torque");
        uplot44 = plot(ss, DesiredTorque(4,:), "-", 'color', "#AF58BA", 'LineWidth',2, 'DisplayName',"Desired torque");
        title("Joint 4", FontSize=fs)
        grid on;
        hold off;

        set(lgd, 'Position', [0.73 0.98 0 0]);

    end
    [pos_desired,~,~] = getDesiredState(trajectory(1:4,1), trajectory(1:4,end), t_now_traj, t_traj);
    err_posk_1 = pos_desired - state(1:4);
    err_posk_2 = err_posk_1;
    eint = zeros(4,1);
    kk = 0;
    % ind_t_elbow_old = 0;
    % ind_n_elbow_old = 0;
    % ind_t_wrist_old = 0;
    % ind_n_wrist_old = 0;
    % syms Z00;
    % eqn = 1/(-Z00+exp(0.4*Z00)) ==model.cr;
    % Z0 = vpasolve(eqn,Z00);
    c_t_elbow = c;
    c_n_elbow = c;
    c_t_wrist = c;
    c_n_wrist = c;

    while 1
        arm = kinematicAnalysis(model, state);
    
        [ind_t_elbow, ind_n_elbow, dir_t_elbow, dir_n_elbow, dev_t_elbow, dist_n_elbow, hitendwall_elbow] = checkContact(arm.pos_elbow, referencePos.Pos_elbow_ref, r);
        [ind_t_wrist, ind_n_wrist, dir_t_wrist, dir_n_wrist, dev_t_wrist, dist_n_wrist, hitendwall_wrist] = checkContact(arm.pos_wrist, referencePos.Pos_wrist_ref, r);
        % if ~ind_t_elbow_old && ind_t_elbow
        %     h_dot_t_elbow0 = dot(arm.Je * arm.state(5:7), dir_t_elbow);
        %     if h_dot_t_elbow0 ==0
        %         h_dot_t_elbow0 = -1e-4;
        %     end
        %     c_t_elbow = Z0*k/h_dot_t_elbow0;
        % end
        % if ~ind_n_elbow_old && ind_n_elbow
        %     h_dot_n_elbow0 = dot(arm.Je * arm.state(5:7), dir_n_elbow);
        %     if h_dot_n_elbow0 ==0
        %         h_dot_n_elbow0 = -1e-4;
        %     end
        %     c_n_elbow = Z0*k/h_dot_n_elbow0;
        % end
        % if ~ind_t_wrist_old && ind_t_wrist
        %     h_dot_t_wrist0 = dot(arm.Jw * arm.state(5:8), dir_t_wrist);
        %     if h_dot_t_wrist0 ==0
        %         h_dot_t_wrist0 = -1e-4;
        %     end
        %     c_t_wrist = Z0*k/h_dot_t_wrist0;
        % end
        % if ~ind_n_wrist_old && ind_n_wrist
        %     h_dot_n_wrist0 = dot(arm.Jw * arm.state(5:8), dir_n_wrist);
        %     if h_dot_n_wrist0 ==0
        %         h_dot_n_wrist0 = -1e-4;
        %     end
        %     c_n_wrist = Z0*k/h_dot_n_wrist0;
        % end

        f_elbow = getContact(arm, r, k, c_t_elbow, c_n_elbow, ind_t_elbow, ind_n_elbow, dir_t_elbow, dir_n_elbow, dev_t_elbow, dist_n_elbow, "elbow");
        f_wrist = getContact(arm, r, k, c_t_wrist, c_n_wrist, ind_t_wrist, ind_n_wrist, dir_t_wrist, dir_n_wrist, dev_t_wrist, dist_n_wrist, "wrist");
        
        % if stateReached(pos_desired(find(impairment==2)), state(find(impairment==2)))
        [pos_desired,vel_desired,acc_desired] = getDesiredState(trajectory(1:4,1), trajectory(1:4,end), t_now, t_traj);
        % end

        err_pos = pos_desired - state(1:4);
        err_vel = vel_desired - state(5:8);

        M = getM(model,state);
        C = getC(model,state);
        G = getG(model,state);

        if controller_type =="pid"
            % PID
            [u, ek_1_new, ek_2_new, eint_new] = pidff(model,err_pos,err_posk_1,err_posk_2,eint,ts,G);
            err_posk_1 = ek_1_new;
            err_posk_2 = ek_2_new;
            eint = eint_new;
        elseif controller_type =="pd"
            % PD+gravity compensation
            u = Kp*err_pos + Kd*err_vel + G; % PD + gravity identification
        elseif controller_type =="impedance"
            u = M*acc_desired + C + G + M*inv(model.MD)*(model.KD* err_vel + model.KP* err_pos);
        end

        if abs(u(find(impairment==1)))> torquelimit(1)
            u(find(impairment==1)) = torquelimit(1);
        end

        if abs(u(find(impairment==2)))> torquelimit(2)
            u(find(impairment==2)) = torquelimit(2);
        end

        Control_Effort = [Control_Effort, u];
        hand_vel = norm(arm.Jw*state(5:8));
        Hand_Vel = [Hand_Vel, hand_vel];
        hand_acc = (hand_vel - hand_vel_old)/ts;
        Hand_Acc = [Hand_Acc, hand_acc];
        hand_vel_old = hand_vel;
        % [state_new, state_dot_new] = RKDiscretize(ts, state, u, f_elbow, f_wrist, model, arm)

        [state_new, state_dot_new] = RKDiscretize(ts, state, u, f_elbow, f_wrist, model, arm);
        [state_new, state_dot_new] = limitState(model,state_new, state_dot_new);

        state = state_new;
        state_dot = state_dot_new;
        Contact_Elbow = [Contact_Elbow, f_elbow];
        Contact_Wrist = [Contact_Wrist, f_wrist];
        Torque_Elbow = [Torque_Elbow, arm.Je'* f_elbow];
        Torque_Wrist = [Torque_Wrist, arm.Jw'* f_wrist];
        Endpoint_elbow = [Endpoint_elbow; arm.pos_elbow];
        Endpoint_wrist = [Endpoint_wrist; arm.pos_wrist];
        U = [U u];
        Time = [Time t_now];
        if debug
            if kk/plot_freq >= 1 || hitendwall_wrist || hitendwall_elbow || length(Time)<2 || abs(state(find(impairment==2))-trajectory(find(impairment==2),end)) <=pi/180*0.1
                figure(fig_simu); 
                config = [zeros(3,1),arm.pos_elbow', arm.pos_wrist'];
                set(armconfig,"XData", config(1,:),"YData",config(2,:),"ZData",config(3,:));
                set(endpoint_elbow,"XData", Endpoint_elbow(:,1),"YData",Endpoint_elbow(:,2),"ZData",Endpoint_elbow(:,3));
                set(endpoint_wrist,"XData", Endpoint_wrist(:,1),"YData",Endpoint_wrist(:,2),"ZData",Endpoint_wrist(:,3));
                if norm(f_elbow) > 1e-10
                    nor_f_elbow = 0.1*f_elbow/norm(f_elbow);
                else
                    nor_f_elbow = zeros(3,1);
                end
    
                if norm(f_wrist) > 1e-10
                    nor_f_wrist = 0.1*f_wrist/norm(f_wrist);
                else
                    nor_f_wrist = zeros(3,1);
                end
    
                set(arr_elbow, 'XData', arm.pos_elbow(1),'YData', arm.pos_elbow(2),'ZData', arm.pos_elbow(3),'UData', double(nor_f_elbow(1)), 'VData', double(nor_f_elbow(2)),'WData', double(nor_f_elbow(3)));
                set(arr_wrist, 'XData', arm.pos_wrist(1),'YData', arm.pos_wrist(2),'ZData', arm.pos_wrist(3),'UData', double(nor_f_wrist(1)), 'VData', double(nor_f_wrist(2)),'WData', double(nor_f_wrist(3)));
                set(title_time, 'String', ['Time: ', num2str(t_now), ' sec'])
    
                set(velplot,"XData", Time,"YData",Hand_Vel);
    
    
                set(accplot,"XData", Time,"YData",Hand_Acc);
    
                set(uplot11,"XData", Time,"YData",Control_Effort(1,:));
                set(uplot21,"XData", Time,"YData",Control_Effort(2,:));
                set(uplot31,"XData", Time,"YData",Control_Effort(3,:));
                set(uplot41,"XData", Time,"YData",Control_Effort(4,:));
    
                set(uplot12,"XData", Time,"YData",Torque_Elbow(1,:));
                set(uplot22,"XData", Time,"YData",Torque_Elbow(2,:));
                set(uplot32,"XData", Time,"YData",Torque_Elbow(3,:));
    
                set(uplot13,"XData", Time,"YData",Torque_Wrist(1,:));
                set(uplot23,"XData", Time,"YData",Torque_Wrist(2,:));
                set(uplot33,"XData", Time,"YData",Torque_Wrist(3,:));
                set(uplot43,"XData", Time,"YData",Torque_Wrist(4,:));
                pause(0.01);
                kk = 0;
            end
        end

        kk = kk+1;
        if t_now > t || hitendwall_elbow || hitendwall_wrist || abs(state(find(impairment==2))-trajectory(find(impairment==2),end)) <=pi/180*0.1
            break;
        end
        t_now = t_now + ts;
        t_now_traj = t_now_traj + ts;
        % if ind_t_elbow
        %     ind_t_elbow_old = 1;
        % else
        %     ind_t_elbow_old = 0;
        % end
        % if ind_n_elbow
        %     ind_n_elbow_old = 1;
        % else
        %     ind_n_elbow_old = 0;
        % end
        % if ind_t_wrist
        %     ind_t_wrist_old = 1;
        % else
        %     ind_t_wrist_old = 0;
        % end
        % if ind_n_wrist
        %     ind_n_wrist_old = 1;
        % else
        %     ind_n_wrist_old = 0;
        % end

    end


end

function [theta,theta_dot,theta_ddot] = getDesiredState(start_pos, end_pos, tnow, t)
    ss = tnow/t;
    if tnow <= t
        theta = start_pos + (start_pos - end_pos).*(15*ss.^4 - 6*ss.^5- 10*ss.^3);
        theta_dot = (start_pos - end_pos).*(60*ss.^3*1/t - 30*ss.^4*1/t- 30*ss.^2*1/t);
        theta_ddot = (start_pos - end_pos).*(180*ss.^2*1/t^2 - 120*ss.^3*1/t^2- 60*ss*1/t^2);
    else
        ss = 1;
        theta = start_pos + (start_pos - end_pos).*(15*ss.^4 - 6*ss.^5- 10*ss.^3);
        theta_dot = (start_pos - end_pos).*(60*ss.^3*1/t - 30*ss.^4*1/t- 30*ss.^2*1/t);
        theta_ddot = (start_pos - end_pos).*(180*ss.^2*1/t^2 - 120*ss.^3*1/t^2- 60*ss*1/t^2);
    end

end


function arm = kinematicAnalysis(model, state)

    lu = model.lu;
    lf = model.lf;
    t1 = state(1);
    t2 = state(2);
    t3 = state(3);
    t4 = state(4);

    T01 = DHConvention(-pi/2, 0, 0, -pi/2 + t1);
    T12 = DHConvention(-pi/2, 0, 0, -pi/2 + t2);
    T02 = T01 * T12;
    T23 = DHConvention(pi/2, 0, 0, t3);
    T03 = T02 * T23;
    T33p = DHConvention(-pi/2, 0, 0, -pi/2);
    T3p4 = DHConvention(0, lu, 0, t4);
    T34 = T33p*T3p4;
    T04 = T01 * T12 * T23 * T34;
    T45 = DHConvention(0, lf, 0, 0);
    T05 = T04 * T45;
    arm.pos_elbow = T04(1:3,4)';
    arm.pos_wrist = T05(1:3,4)';
    arm.state = state;
    arm.R05 = T05(1:3,1:3);
    arm.Jw = calculateJacobian(T01,T02,T03,T04,T05,5);
    arm.Je = calculateJacobian(T01,T02,T03,T04,T05,4);
    arm.Jew = [T01(1:3,3) T02(1:3,3) T03(1:3,3) T04(1:3,3)];

end

function referencePos = getTrajectoryCartesian(model, trajectory)

    lu = model.lu;
    lf = model.lf;
    PoseRef = [];
    PoswRef = [];
    for i = 1:length(trajectory(1,:))
        t1 = trajectory(1,i);
        t2 = trajectory(2,i);
        t3 = trajectory(3,i);
        t4 = trajectory(4,i);
    
        T01 = DHConvention(-pi/2, 0, 0, -pi/2 + t1);
        T12 = DHConvention(-pi/2, 0, 0, -pi/2 + t2);
        T23 = DHConvention(pi/2, 0, 0, t3);
        T33p = DHConvention(-pi/2, 0, 0, -pi/2);
        T3p4 = DHConvention(0, lu, 0, t4);
        T34 = T33p*T3p4;
        T04 = T01 * T12 * T23 * T34;
        T45 = DHConvention(0, lf, 0, 0);
        T05 = T04 * T45;
        PoseRef = [PoseRef;T04(1:3,4)'];
        PoswRef = [PoswRef;T05(1:3,4)'];
    end
    referencePos.Pos_elbow_ref = PoseRef;
    referencePos.Pos_wrist_ref = PoswRef;
end

function t = createTunnel(fig, path, RGB, name, r, n)
    figure(fig);
    num_segments = size(path, 2);
    tunnel_x = zeros(n, num_segments);
    tunnel_y = zeros(n, num_segments);
    tunnel_z = zeros(n, num_segments);
    
    theta = linspace(0, 2*pi, n);
    
    tangent_vectors = diff(path, 1, 2); 
    tangent_vectors = [tangent_vectors, tangent_vectors(:, end)]; 
    
    for i = 1:num_segments
        center = path(:, i);
        
        t = tangent_vectors(:, i);
        t = t / norm(t); 
        
        if abs(t(3)) < 1e-10 && abs(t(2)) < 1e-10 
            normal = cross(t, [0; 1; 0]);
        else
            normal = cross(t, [1; 0; 0]);
        end
        normal = normal / norm(normal); 
    
        binormal = cross(t, normal);
        
        circle_x = r * cos(theta);
        circle_y = r * sin(theta);
        circle_z = zeros(size(circle_x));
        
        for j = 1:n
            circle_point = normal * circle_x(j) + binormal * circle_y(j);
            tunnel_x(j, i) = circle_point(1) + center(1);
            tunnel_y(j, i) = circle_point(2) + center(2);
            tunnel_z(j, i) = circle_point(3) + center(3);
        end
    end
    
    t = surf(tunnel_x, tunnel_y, tunnel_z, 'FaceColor', RGB, ... 
         'FaceAlpha', 0.3, 'EdgeColor', 'none', DisplayName=name);
end

function arr = plotContact(start_point, f, direction, color,name)
    arr = quiver3(start_point(1), start_point(2), start_point(3), f*direction(1), f*direction(2), f*direction(3), "color", color, 'LineWidth', 3, 'MaxHeadSize', 2, DisplayName=name);
end

function [ind_t, ind_n, dir_t, dir_n, dev_t, dist_n, hitendwall] = checkContact(pos, Pos_ref, r)
    hitendwall = 0;
    k = dsearchn(Pos_ref, pos);

    if k == length(Pos_ref(:,1)) % the end point
        tan = (Pos_ref(k-1,:) - Pos_ref(k,:))/norm(Pos_ref(k-1,:) - Pos_ref(k,:));
    elseif k == 1 % the start point
        tan = (Pos_ref(k+1,:) - Pos_ref(k,:))/norm(Pos_ref(k+1,:) - Pos_ref(k,:));
    else
        ind_t = 0;
        dir_t = zeros(3,1);
        dir_n = (Pos_ref(k,:) - pos) / norm(Pos_ref(k,:) - pos);
        dir_n = dir_n';
        dev_t = 0;
        dist_n = norm(Pos_ref(k,:) - pos);
        if r-dist_n <=0
            ind_n = 1;
        else
            ind_n = 0;
        end
        return;
    end
    displacement = Pos_ref(k,:) - pos;
    displacement_t = dot(displacement, tan);
    displacement_n = displacement - displacement_t*tan;
    if displacement_t>=0 % outside the tunnel in tangential direction
        ind_t = 1;
        dir_t = tan';
        dev_t = -displacement_t;
        if k == length(Pos_ref(:,1))
            hitendwall = 1;
        end
    else
        ind_t = 0;
        dir_t = zeros(3,1);
        dev_t = 0;
    end
    dir_n = displacement_n/norm(displacement_n);
    dir_n = dir_n';
    dist_n = norm(displacement_n);
    if r-dist_n <=0
        ind_n = 1;
    else
        ind_n = 0;
    end
end

function [u, ek_1_new, ek_2_new, eint_new] = pidff(model,ek,ek_1,ek_2,eint,ts,grav)
    u = model.Kp*ek + model.Ki*(ts/2 * (ek+2*ek_1+ek_2) + eint) +model.Kd/ts*(ek-ek_1) + grav;
    ek_1_new = ek;
    ek_2_new = ek_1;
    eint_new = ts/2 * (ek+2*ek_1+ek_2) + eint;

end

function [new_state, new_state_dot] = limitState(model,state, state_dot)

    state(find((state(1:4)-model.jointlower)<0)) = model.jointlower(find((state(1:4)-model.jointlower)<0));
    state(find(state(4+find((state(1:4)-model.jointlower)<0)) <0)) = 0;
    state_dot(find(state_dot(4+find((state(1:4)-model.jointlower)<0)) <0)) = 0;

    state(find((state(1:4)-model.jointupper)>0)) = model.jointupper(find((state(1:4)-model.jointupper)>0));
    state(find(state(4+find((state(1:4)-model.jointupper)>0)) >0)) = 0;
    state_dot(find(state_dot(4+find((state(1:4)-model.jointupper)>0)) >0)) = 0;
    
    new_state = state;
    new_state_dot = zeros(8,1);
    new_state_dot(1:4) = state(5:8);
    new_state_dot(5:8) = state_dot(5:8);
end

function ind = stateReached(x_d,x)
    if norm(x_d-x)<0.1*pi/180
        ind = 1;
    else 
        ind = 0;
    end
end
