% xzcwss% main
close all
clear
clc
% 1. Point MATLAB to your specific Gurobi installation directory
addpath('C:\gurobi1302\win64\matlab')

% 2. Save this path so MATLAB remembers it every time it opens
savepath

% 3. Run Gurobi's internal setup verification script
gurobi_setup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% motion configuration
start_pos = [0;0;0;0;0;0;0]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
model.r = 0.02;
model.k = 50000;
model.c = 800000;
model.cr = 0.3;
% for pid or pd
model.Kp = 50;
model.Kd = 1;
model.Ki = 1;
% for impedance control
model.MD = diag([10 10 10 10 10 10 10 10]);
model.KD = diag([100 100 100 100 100 100 100]);
model.KP = diag([10000 10000 10000 10000 10000 10000 10000]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% joints of interest
group1 = 2;
group2 = 4; 
impairment = [0,1,0,2,0,1,0]; 
torquelimit = [0.5,20];  %
% for simulator
performance.overalpct = 0.95;
performance.truepct = 0.65;
performance.ratio = 1.5;
% for heuristic
percentage = [0.7 0.7];
threshold = [0.5 2 0.5 2 0.5 2 0.5]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation time step
sim.ts = 0.001;
sim.rt = 5;
ts_traj = sim.rt*sim.ts;
sim.t_tot_sim = 10;
n = 100; % dilute the plot
t_tot = 10;
plot_freq = 1/(10*sim.ts);
controller_type = "potential";

model.vm = 0.5;
model.vd = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fontsize = 16;
legendfontsize = 12;
q1_l = -deg2rad(180);
q1_u = deg2rad(30);
q2_l = -deg2rad(50);
q2_u = deg2rad(180);
q3_l = -deg2rad(100);
q3_u = deg2rad(80);
q4_l = 0;
q4_u = deg2rad(145); 
partition = 2;
feasibleWorkspace = [];
Results = [];
model.jointupper = pi/180*[30; 180; 80; 145;30; 180; 80]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% [-pi/2;0;0;pi/4];
model.jointlower = pi/180*[-180; -50; -100; 0;-180; -50; -100]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% [-pi/2;pi/2;0;pi/2];
model.ts = sim.ts;
model.lu = 0.3;
model.lf = 0.33;
model.lw = 0.05;
model.m1 = 0;
model.m2 = 0;
model.m4 = 0;
model.m6 = 0;
model.m3 = 1.4;
model.m5 = 1.1; 
model.m7 = 0.3;
model.Lcmu = 0.11; 
model.Lcmf = 0.16; 
model.Lcmw = 0.02;
model.lvua = 0.03/model.lu;
model.lvla = 0.03/model.lf;
model.lvw = 0.03/model.lw;
model.lvua21 = model.lvua/2;
model.lvua22 = model.lvua/2;
model.lvla2 = model.lvla/2;
model.lvw2 = model.lvw/2;
model.g = 9.81;    
Ixx_3 = 0.01;
Iyy_3 = 0.027;
Izz_3 = 0.027;
Ixx_5 = 0.015;
Iyy_5 = 0.045;
Izz_5 = 0.045; 
Ixx_7 = 0.003;
Iyy_7 = 0.003;
Izz_7 = 0.003;
I3c3 = [Ixx_3 0 0;0 Iyy_3 0;0 0 Izz_3];
I5c5 = [Ixx_5 0 0;0 Iyy_5 0;0 0 Izz_5];
I7c7 = [Ixx_7 0 0;0 Iyy_7 0;0 0 Izz_7];
model.Ic11 = zeros(3);
model.Ic22 = zeros(3);
model.Ic44 = zeros(3);
model.Ic66 = zeros(3);
model.Ic33 = [Ixx_3 0 0;0 Iyy_3 0;0 0 Izz_3];
model.Ic55 = [Ixx_5 0 0;0 Iyy_5 0;0 0 Izz_5];
model.Ic77 = [Ixx_7 0 0;0 Iyy_7 0;0 0 Izz_7];
%% Test
end_pos = (pi/180)*[-140;36.4;122.4;46.8;-25.2;-7.2;-3.6]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[trajectory, ss] = getTrajectoryProperty(model, start_pos, end_pos, ts_traj, t_tot, n, fontsize, legendfontsize, 1);
%% Trajectory smoothing
%{
% 1. Create a dummy cubic-like step command (causes an acceleration spike)
dt = ts_traj;                 % 1 kHz loop
%tsm = 0:dt:sim.t_tot_sim;                 % 2 seconds duration
q_recorded = trajectory (1:7,:);
%q_cmd(tsm >= 0.5) = 1.0;      % Sharp step change at 0.5 seconds

% 2. Filter the commands (Tuned to 50 rad/s)

%[q_smooth, v_smooth, a_smooth] = smoothReferenceFilter(q_cmd, dt, omega_n);

% 1. Define your parameters
            % Your controller sampling period (e.g., 1000 Hz)
           % Cutoff frequency (rad/s) - Adjust higher for tighter tracking

% 2. Create the continuous second-order critically damped filter
% G(s) = omega_n^2 / (s^2 + 2*omega_n*s + omega_n^2)
num_s = [omega_n^2];
den_s = [1, 2*omega_n, omega_n^2];
sys_s = tf(num_s, den_s);

% 3. Convert to discrete-time for filtfilt (Tustin/Bilinear method)
sys_d = c2d(sys_s, dt, 'tustin');
[b, a] = tfdata(sys_d, 'v');

% 4. Run zero-phase smoothing on your recorded data array
% Assuming 'q_recorded' is a vector of your practice run positions
q_smooth = filtfilt(b, a, q_recorded);

% 5. Extract perfectly continuous Velocity and Acceleration from the smooth path
v_smooth = gradient(q_recorded) / dt;
a_smooth = gradient(q_recorded) / dt;

trajectory2 = [q_smooth;v_smooth;a_smooth];
% 3. Plot the spikeless acceleration output
%}
% Replace Steps 2, 3, and 4 entirely with this single line:
% 'gaussian' or 'rlowess' smoothing preserves edge points beautifully
dt = ts_traj;  
omega_n = 30; 
q_recorded = trajectory (1:7,:);
%q_smooth = smoothdata(q_recorded, 2,'gaussian', 10); 

% Proceed to Step 5
%v_smooth = gradient(q_smooth) / dt;
%a_smooth = gradient(v_smooth) / dt;
%[q_clean, v_smooth, a_smooth, j_smooth] = remove_jerk_spikes_7joints(q_recorded, dt);

% Proceed to Step 5
%v_smooth = gradient(q_smooth) / dt;
%a_smooth = gradient(v_smooth) / dt;
[q_smooth, v_smooth, a_smooth] = quintic_spline_trajectory(q_recorded, dt);
trajectory2 = [q_smooth; v_smooth; a_smooth];

%{
figure;
plot(t, a_smooth, 'LineWidth', 2);
grid on;
title('Filtered Reference Acceleration (Completely Spikeless)');
xlabel('Time (s)');
ylabel('Acceleration (rad/s^2)');
%}
%% Simulator for arm
% [Torque_Wrist, U] = armSimulator(ss, t_tot, sim, trajectory, impairment, torquelimit, model,controller_type, plot_freq, 1);
% [ind1, ind2] = simulationEvaluator(U,Torque_Wrist,impairment, performance, 1);
%% Simulator for exo
[Torque_Wrist, U] = exoSimulator(ss, t_tot, sim, trajectory, trajectory2, impairment, torquelimit, model,controller_type, plot_freq, 1,start_pos, end_pos);


%% Heuristic

% [indicator, results] = trajectoryEvaluator(ss, model, group1, group2, trajectory, percentage, threshold, 1);


%% % Discretize the joint space
% for q1 = q1_l:(q1_u-q1_l)/partition:q1_u
%     for q2 = q2_l:(q2_u-q2_l)/partition:q2_u
%         for q3 = q3_l:(q3_u-q3_l)/partition:q3_u
%             for q4 = q4_l:(q4_u-q4_l)/partition:q4_u
% 
%                 end_pos = [q1;q2;q3;q4];
%                 if all(start_pos([group1 group2])-end_pos([group1 group2]) ~= 0)
%                     fprintf('========================================================\nEvaluating %s\n', num2str(end_pos'));
%                     [trajectory, ss] = getTrajectoryProperty(arm_length, start_pos, end_pos, ts, t_tot, fontsize, legendfontsize, 0);
%                     disp("Trajectory generated successfully");
%                     [indicator, results] = trajectoryEvaluator(ss, mass, inertia, cm_length, arm_length, group1, group2, trajectory, percentage, threshold, 0);
%                     if indicator
%                         q_feasible = [q1,q2,q3,q4];
%                         Results = [Results; results];
%                         feasibleWorkspace = [feasibleWorkspace; q_feasible];
%                         disp("Determined as a feasible end configuration!");
%                     else
%                         disp("Determined as an infeasible end configuration!");
%                     end
%                 end
% 
%             end
%         end
%     end
% end
% 
% columnNames = getColumnName(group1, group2);
% output = [columnNames; num2cell([feasibleWorkspace results])];
% writematrix(output, "FeasibleWorkspace.csv");
% 



%%
function columnNames = getColumnName(group1, group2)
    l1 = length(group1);
    l2 = length(group2);
    
    
    columnTan1 = [];
    columnNor1 = [];
    for i = 1:l1
        columnTan1 = [columnTan1 string(sprintf('Percentage_Tan_M_J%d', group1(i)))];
        columnNor1 = [columnNor1 string(sprintf('Percentage_Nor_M_J%d', group1(i)))];
    end
    columnTan2 = [];
    columnNor2 = [];
    for i = 1:l2
        columnTan2 = [columnTan2 string(sprintf('Percentage_Tan_S_J%d', group2(i)))];
        columnNor2 = [columnNor2 string(sprintf('Percentage_Nor_S_J%d', group2(i)))];
    end
    
    columnTan3 = [];
    columnNor3 = [];
    for i = 1:l2
        columnTan3 = [columnTan3 string(sprintf('Percentage_Tan_M_J%d', group2(i)))];
        columnNor3 = [columnNor3 string(sprintf('Percentage_Nor_M_J%d', group2(i)))];
    end
    columnTan4 = [];
    columnNor4 = [];
    for i = 1:l1
        columnTan4 = [columnTan4 string(sprintf('Percentage_Tan_S_J%d', group1(i)))];
        columnNor4 = [columnNor4 string(sprintf('Percentage_Nor_S_J%d', group1(i)))];
    end

    trans1 = [];
    trans2 = [];

    for i = 1:l1
        for j = 1:l2
            trans1 = [trans1 string(sprintf('TransmissionRatio_Tan_M_J%d_S_J%d ', group1(i),group2(j))) string(sprintf('TransmissionRatio_Nor_M_J%d_S_J%d ', group1(i),group2(j)))];
        end
    end
    for i = 1:l2
        for j = 1:l1
            trans2 = [trans2 string(sprintf('TransmissionRatio_Tan_M_J%d_S_J%d ', group2(i),group1(j))) string(sprintf('TransmissionRatio_Nor_M_J%d_S_J%d ', group2(i),group1(j)))];
        end
    end
    columnNames = ["q1 [deg]", "q2 [deg]", "q3 [deg]", "q4 [deg]", columnTan1, columnNor1, columnTan2, columnNor2, columnTan3, columnNor3,columnTan4, columnNor4, trans1, trans2];

end


function [trajectory, ss] = getTrajectoryProperty(model, start_pos, end_pos, ts, t_tot, n, fontsize, legendfontsize, debug)
    ss = 0:ts:t_tot;    
    lu = model.lu;
    lf = model.lf;
    lw = model.lw;

    trajectory = trajectoryGenerationMAJ(start_pos, end_pos, t_tot, ts);
    theta = trajectory(1:7,:);
    theta_dot = trajectory(8:14,:);
    
    
    idx = [1, 2:n:length(theta(1,:))-1, length(theta(1,:))];
    theta = theta(:,idx);
    theta_dot = theta_dot(:,idx);
    sst = ss(idx);

    Pos_elbow = zeros(3,length(idx));
    Pos_wrist = zeros(3,length(idx));
    Pos_tip = zeros(3,length(idx));
    Vel = zeros(1,length(idx));
    Acc = zeros(1,length(idx));

    for i = 1:length(sst)
        
        T01 = DHConvention(-pi/2, 0, 0, -pi/2 + theta(1,i)); %
        T12 = DHConvention(-pi/2, 0, 0, -pi/2 + theta(2,i)); %
        T02 = T01*T12;
        T23 = DHConvention(pi/2, 0, 0, theta(3,i)); %
        T03 = T02*T23;
        T3_p1= DHConvention(-pi/2, 0, 0, -pi/2);
        T0_p1= T03 * T3_p1;
        Tp1_4 = DHConvention(-pi/2, lu, 0, pi/2+theta(4,i)); %
        T04 = T0_p1 * Tp1_4;
        T45 = DHConvention(-pi/2, 0, 0, theta(5,i)); 
        T05 = T04 * T45;
        T5_p2 = DHConvention(-pi/2, 0, 0, pi/2); 
        T0_p2 = T05 * T5_p2;
        Tp2_6 = DHConvention(0, lf, 0, theta(6,i)); 
        T06 = T0_p2 * Tp2_6;
        T67 = DHConvention(-pi/2, 0, 0, theta(7,i)); 
        T07 = T06 * T67;
        T78= DHConvention(0, lw, 0, 0); 
        T08=T07*T78;
        


    
        pos_elbow = T04(1:3,4);
        pos_wrist= T06(1:3,4);
        pos_tip = T08(1:3,4);
        
        Pos_elbow(:,i) = pos_elbow;
        Pos_wrist(:,i) = pos_wrist;
        Pos_tip(:,i) = pos_tip;
    
        Jv = calculateJacobian(T01,T02,T03,T04,T05,T06,T07,T08,8);
        Vel(i) = norm(Jv * theta_dot(:,i));
    end
    
    
    for j = 2:length(Vel)-1
        Acc(j) = (Vel(j+1) - Vel(j-1)) / (sst(j+1)-sst(j-1));
    end
    
    Acc(1) = (Vel(2) - Vel(1)) / (sst(2)-sst(1));     
    Acc(end) = (Vel(end) - Vel(end-1)) / (sst(end)-sst(end-1));  
    
    if debug
        figure()
        
        h1 = plot(sst,rad2deg(theta(1,:)), 'LineWidth', 2.5);
        hold on
        h2 = plot(sst,rad2deg(theta(2,:)), 'LineWidth', 2.5);
        hold on
        h3 = plot(sst,rad2deg(theta(3,:)), 'LineWidth', 2.5);
        hold on
        h4 = plot(sst,rad2deg(theta(4,:)), 'LineWidth', 2.5);        
        hold on
        h5 = plot(sst,rad2deg(theta(5,:)), 'LineWidth', 2.5);
        hold on
        h6 = plot(sst,rad2deg(theta(6,:)), 'LineWidth', 2.5);
        hold on
        h7 = plot(sst,rad2deg(theta(7,:)), 'LineWidth', 2.5);
        hold off
        
        set(h1, 'Color', [h1.Color, 0.3]);  
        set(h2, 'Color', [h2.Color, 0.3]);  
        set(h3, 'Color', [h3.Color, 0.3]);  
        set(h4, 'Color', [h4.Color, 0.3]);  
        set(h5, 'Color', [h5.Color, 0.3]);  
        set(h6, 'Color', [h6.Color, 0.3]);  
        set(h7, 'Color', [h7.Color, 0.3]);  
        
        xlabel("Time [s]", 'FontSize', fontsize)
        ylabel("Angular position [deg]", 'FontSize', fontsize)
        legend(["\theta_1", "\theta_2", "\theta_3", "\theta_4", "\theta_5", "\theta_6", "\theta_7"], 'FontSize', legendfontsize)
        
        grid on
        
        figure()
        Pos_elbow0 = Pos_elbow(:,1);
        Pos_wrist0 = Pos_wrist(:,1);
        Pos_tip0 = Pos_tip(:,1);
        config0 = [[0;0;0] Pos_elbow0 Pos_wrist0 Pos_tip0];
        
        Pos_elbow25 = Pos_elbow(:, round(0.25*t_tot/ts/n));
        Pos_wrist25 = Pos_wrist(:, round(0.25*t_tot/ts/n));
        Pos_tip25 = Pos_tip(:, round(0.25*t_tot/ts/n));
        config25 = [[0;0;0] Pos_elbow25 Pos_wrist25 Pos_tip25];
        
        Pos_elbow50 = Pos_elbow(:, round(0.5*t_tot/ts/n));
        Pos_wrist50 = Pos_wrist(:, round(0.5*t_tot/ts/n));
        Pos_tip50 = Pos_tip(:, round(0.5*t_tot/ts/n));
        config50 = [[0;0;0] Pos_elbow50 Pos_wrist50 Pos_tip50];
        
        Pos_elbow75 = Pos_elbow(:, round(0.75*t_tot/ts/n));
        Pos_wrist75 = Pos_wrist(:, round(0.75*t_tot/ts/n));
        Pos_tip75 = Pos_tip(:, round(0.75*t_tot/ts/n));
        config75 = [[0;0;0] Pos_elbow75 Pos_wrist75 Pos_tip75];
        
        Pos_elbow100 = Pos_elbow(:, end);
        Pos_wrist100 = Pos_wrist(:, end);
        Pos_tip100 = Pos_tip(:, round(0.25*t_tot/ts/n));
        config100 = [[0;0;0] Pos_elbow100 Pos_wrist100 Pos_tip100];
        
        
        subplot(1,1,1)
        plot3(Pos_elbow(1,:),Pos_elbow(2,:),Pos_elbow(3,:),":", 'LineWidth', 2.5)
        hold on
        plot3(Pos_wrist(1,:),Pos_wrist(2,:),Pos_wrist(3,:),":", 'LineWidth', 2.5)
        hold on
        plot3(Pos_tip(1,:),Pos_tip(2,:),Pos_tip(3,:),":", 'LineWidth', 2.5)
        hold on
        plot3(config0(1,:),config0(2,:),config0(3,:), 'LineWidth', 1.5)
        hold on
        plot3(config25(1,:),config25(2,:),config25(3,:), 'LineWidth', 1.5)
        hold on
        plot3(config50(1,:),config50(2,:),config50(3,:), 'LineWidth', 1.5)
        hold on
        plot3(config75(1,:),config75(2,:),config75(3,:), 'LineWidth', 1.5)
        hold on
        plot3(config100(1,:),config100(2,:),config100(3,:), 'LineWidth', 1.5)
        hold off
        
        axis equal
        xlabel("x [m]", 'FontSize', fontsize)
        ylabel("y [m]", 'FontSize', fontsize)
        zlabel("z [m]", 'FontSize', fontsize)
        legend(["Hand position", "Elbow position","Tip position", "Configuration at 0% of progress", ...
            "Configuration at 25% of progress","Configuration at 50% of progress", ...
            "Configuration at 75% of progress","Configuration at 100% of progress"], 'FontSize', legendfontsize)
        view(135, 30); 
        %{
        subplot(2,1,2)
        yyaxis left;
        plot(sst, Vel, '-b', 'LineWidth', 1.5)
        ylabel("Hand speed [m/s]", 'FontSize', fontsize)
        yyaxis right;
        plot(sst, Acc, '-r', 'LineWidth', 1.5)
        ylabel("Hand acceleration [m/s^2]", 'FontSize', fontsize)
        xlabel("Time [s]", 'FontSize', fontsize)
        grid on;
        %}
    end
    
end
