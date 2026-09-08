start_pos = [0;0;0;0;0;0;0]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
end_pos = (pi/180)*[-140;36.4;122.4;46.8;-25.2;-7.2;-3.6];


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

Ixx_3 = 0.01;
Iyy_3 = 0.027;
Izz_3 = 0.027;
Ixx_5 = 0.015;
Iyy_5 = 0.045;
Izz_5 = 0.045; 
Ixx_7 = 0.003;
Iyy_7 = 0.003;
Izz_7 = 0.003;

% Parameters
    Kuamx1= 2.5; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Kuamx2= 2.5;
    Kuamn1= 1.5;
    Kuamn2= 1.5;
    Klamx= 2.0;
    Kwmx= 1.5;
    Klamn= 1.0;
    Kwmn= 0.5; 
    Kt=200;
    
    Sua = 3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Sla = 2.5; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Sw = 2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    60<=Sua
    45<=Sla
    30<=Sw
    
    J_max = 0.5 ;
    UUP = 0;
    sigma= 2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    winsize = 15;
    K_intent  = diag([15.0, 15.0, 12.0, 12.0,  6.0,  4.0,  3.0]);  
    
        % 1. Healthy Parameter Matrices (Relaxed State)
    K_healthy = diag([20.0, 18.0, 12.0, 16.0,  5.0,  2.0,  1.5]); % Nm/rad
    B_healthy = diag([ 3.0,  2.5,  2.0,  2.2,  0.6,  0.2,  0.15]); % Nms/rad
    q_rest_healthy = [0.0;  0.0;  0.0; 0.44;  0.0;  0.0;  0.0]; % Radians

    % 2. Pathological Parameter Matrices (Severe Chronic Stroke State)
    K_spastic = diag([110.0, 95.0, 65.0, 115.0, 30.0, 22.0, 14.0]); % Nm/rad
    B_spastic = diag([ 11.0,  9.0,  6.5,  10.0,  3.5,  1.8,  1.2]); % Nms/rad
    q_rest_stroke  = [0.26; 0.17; 0.78;  1.57; 0.52; -0.78; 0.0]; % Radians

    % 3. Controller Active Intent Gain Matrix (Healthy target tracking effort)
    % This models how aggressively a cooperative user fights to rejoin the path
    K_intent  = diag([15.0, 15.0, 12.0, 12.0,  6.0,  4.0,  3.0]); 