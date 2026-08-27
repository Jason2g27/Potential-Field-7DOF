function [Torque_Wrist, U] = exoSimulator(ss, t_traj, sim, trajectory, trajectory2, impairment2, torquelimit, model, controller_type, plot_freq, debug,start_pos, end_pos)
    
    hist = zeros (7,6);
    Torque_Wrist=0;
    brc=0;
    fc = 10;  
    [bbuf, abuf] = bessel3_init(fc, sim.ts);
    % impairment  2 - mild, 1 - severe, 0 - not care
    fs = 9;  % font size
    % initialize 
    t = sim.t_tot_sim;
    ts = sim.ts;
    t_now = 0;
    t_now_traj = 0;
    referencePos = getTrajectoryCartesian(model, trajectory);
    state = zeros(14,1);
    state_dot= zeros(14,1);
    state(1:7) = trajectory(1:7,1);
    state_dot (1:7) = trajectory(8:14,1);
    state(8:14) = trajectory(8:14,1);
    %state(15:21) = trajectory(15:21,1);
    %Contact_Elbow = [];
    %Contact_Wrist = [];
    Control_Effort = [];
    %Torque_Elbow = [];
    %Torque_Wrist = [];
    J_Vel = [];
    J_Pos=[];
    %Hand_Acc = [];
    %hand_vel_old = 0;
    Endpoint_elbow = [];
    Endpoint_wrist = [];
    Endpoint_tip = [];
    U = [];
    Time = [];
    r = model.r;
    k = model.k;
    c = model.c;
    Kd = model.Kd;
    Kp = model.Kp;
    J_max = 0.5 ;
    UUP = 0;
    sigma= 2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ind_t = [1, 2:plot_freq:length(ss)-1, length(ss)];
    % ss = ss(ind_t);
    % trajectory = trajectory(:,ind_t);
    tau_hu = zeros(7, 1);
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
    impairment.shoulder = 0;
    impairment.elbow = 0;
    impairment.wrist = 0;
    [DesiredTorque,DesiredTorque2, Kua, Kla, Kw, ualua, ualla, ualw, calua, calla, calw, BBBB, ThetauapdN1, ThetauapdN2 ,ThetalapdN, ThetawpdN] = Potentialfieldinitator(trajectory,trajectory2, model,sigma,[],[]);
   
    
    %Kua =Kua 
    
    mkkk=length(trajectory(1,:));
    winsize = 15;
    %{ 
    wuap = zeros(mkkk,1);
    wlap = zeros(mkkk,1);
    wwp = zeros(mkkk,1);
    uuap = zeros(mkkk);
    ulap = zeros(mkkk);
    uwp = zeros(mkkk);
    %}
    wuap = zeros(winsize+1);
    wlap = zeros(winsize+1);
    wwp = zeros(winsize+1);
    uuap = zeros(winsize+1);
    ulap = zeros(winsize+1);
    uwp = zeros(winsize+1);   
    %uuak = zeros(3,1);
    %ulak = zeros(2,1);
    %uwk = zeros(2,1);
    %mlknmn = 1;
    
    Thetauap = state (1:3);
    Thetalap = state (4:5);
    Thetawp = state (6:7);



    Thetauapp = trajectory(1:3,:);
    Thetalapp = trajectory (4:5,:);
    Thetawpp = trajectory (6:7,:);
    Thetauappdt = trajectory(8:10,:);
    Thetalappdt = trajectory(11:12,:);
    Thetawppdt = trajectory(13:14,:);
    Thetappdt = trajectory(8:14,:);
    
    Thetauapp2 = trajectory2(1:3,:);
    Thetalapp2 = trajectory2 (4:5,:);
    Thetawpp2 = trajectory2 (6:7,:);
    Thetauappdt2 = trajectory2(8:10,:);
    Thetalappdt2 = trajectory2(11:12,:);
    Thetawppdt2 = trajectory2(13:14,:);
    Thetappdt2 = trajectory2(8:14,:);    
    
     % Plotting the Potential field begining;
     Theta21 = 0:0.05:2*pi; 
     mmmL= length (Theta21);
     jjjL= mmmL;
     iiiL=mkkk; 
     wwp2 = zeros (iiiL,mmmL,mmmL);
     wlap2 = zeros (iiiL,mmmL,mmmL);
     ulap2 = zeros (iiiL,mmmL,mmmL);
     uwp2 = zeros (iiiL,mmmL,mmmL);
     UUla2 = zeros (mmmL,mmmL);
     UUw2= zeros (mmmL,mmmL);


     for mmm=1:mmmL
         for jjj=1:jjjL
             for iii=1:iiiL 
                wlap2 (iii,mmm,jjj)= exp(-(1/(sigma^2))*([Theta21(mmm);Theta21(jjj)]-Thetalapp(:,iii))'*([Theta21(mmm);Theta21(jjj)]-Thetalapp(:,iii)));
                wwp2 (iii,mmm,jjj)= exp(-(1/(sigma^2))*([Theta21(mmm);Theta21(jjj)]-Thetawpp(:,iii))'*([Theta21(mmm);Theta21(jjj)]-Thetawpp(:,iii))); 
                ulap2(iii,mmm,jjj)= ualla(iii)+0.5*([Theta21(mmm);Theta21(jjj)]-Thetalapp(:,iii))'*Kla(:,2*iii-1:2*iii)*([Theta21(mmm);Theta21(jjj)]-Thetalapp(:,iii));
                uwp2(iii,mmm,jjj)= ualw(iii)+0.5*([Theta21(mmm);Theta21(jjj)]-Thetawpp(:,iii))'*Kw(:,2*iii-1:2*iii)*([Theta21(mmm);Theta21(jjj)]-Thetawpp(:,iii));
             end
             wlap2(:,mmm,jjj)= wlap2(:,mmm,jjj)/sum(wlap2(:,mmm,jjj));
             wwp2(:,mmm,jjj)= wwp2(:,mmm,jjj)/sum(wwp2(:,mmm,jjj));
             for iii=1:iiiL
                UUla2(mmm,jjj)= UUla2(mmm,jjj) + wlap2(iii,mmm,jjj)* ulap2(iii,mmm,jjj);
                UUw2(mmm,jjj)=  UUw2(mmm,jjj) + wwp2(iii,mmm,jjj)* uwp2(iii,mmm,jjj);
             end
         end
     end


     fig_simu01 = figure();
     figure(fig_simu01);
     surf(Theta21,Theta21,UUla2);
     title('Elbow potential field');
     fig_simu02 = figure();
     figure(fig_simu02);    
     surf(Theta21,Theta21,UUw2);
     title('Wrist potential field');
     % Plotting the Potential field end;

    %diffsua = Thetauapp - Thetauap;
    %sqDistancesua = sum(diffsua.^2, 2);
    %[~, closestIdxua] = min(sqDistancesua);
    %prevIdxua = closestIdxua;
    
    %diffsla = Thetalapp - Thetalap;
    %sqDistancesla = sum(diffsla.^2, 2);
    %[~, closestIdxla] = min(sqDistancesla);
    %prevIdxla = closestIdxla;

     thetadrua = [1/sqrt(3);1/sqrt(3);1/sqrt(3)];
     thetadrla = [1/sqrt(2);1/sqrt(2)];
     thetadrw = [1/sqrt(2);1/sqrt(2)];
     %[prevIdxua,dist1ua ] = dsearchn(Thetauapp',Thetauap');
     %[prevIdxla,dist1la] = dsearchn(Thetalapp',Thetalap');
     %[prevIdxw, dist1w] = dsearchn(Thetawpp',Thetawp');
    prevIdxua=1;
    prevIdxla=1;
    prevIdxw=1;
    %{
    %prevIdxua=prevIdxua
    prevIdxua= int64 (prevIdxua);
    prevIdxla= int64 (prevIdxla);
    prevIdxw= int64 (prevIdxw);
    %}

    Sua = 3; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Sla = 2.5; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Sw = 2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    Ti_prevua= Thetauapp(:,2)-Thetauapp(:,1); 
    Ni_prevua= BBBB(:,1);
    Pi_prevua = Thetauapp(:,1);
    %P2 = Thetauap; %only for the first moment  + every turn

    %Dtt =[Dua, zeros(3,4);zeros (2,3), Dla, zeros(2,2); zeros (2,5), Dw];

    

    
    
    if debug
        arm = kinematicAnalysis(model, state);
        fig_simu = figure();
        figure(fig_simu);
        set(fig_simu, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); 
        pause(0.3);
        %subplot(6,2,[1,3,5]);
        hold on;
        axis equal;
        view(135, 30); 
        xlabel("X [m]", 'FontSize', fs);
        ylabel("Y [m]", 'FontSize', fs);
        zlabel("Z [m]", 'FontSize', fs);
        % axis([-0.7 0.7 -0.1 0.7 -0.7 0.7]); 
        title_time = title(['Time: ', num2str(t_now), ' sec'], 'FontSize', fs); 
        %createTunnel(fig_simu, referencePos.Pos_elbow_ref', "#9EC9E2", "Virtual tunnel for elbow", r, 40);
        %createTunnel(fig_simu, referencePos.Pos_wrist_ref', "#9CCEA7", "Virtual tunnel for wrist", r, 40);
        config = [zeros(3,1),arm.pos_elbow', arm.pos_wrist',arm.pos_tip'];
        armconfig = plot3(config(1,:),config(2,:),config(3,:),"-", "color", '#00CD6C', "LineWidth",5, "DisplayName", "Arm configuration");
        % set(get(get(armconfig, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off'); 

        ref_traj_elbow = plot3(referencePos.Pos_elbow_ref(:,1),referencePos.Pos_elbow_ref(:,2),referencePos.Pos_elbow_ref(:,3),":", "color", '#FFC61E', "LineWidth",3, "DisplayName", 'Reference elbow trajectory');
        ref_traj_wrist = plot3(referencePos.Pos_wrist_ref(:,1),referencePos.Pos_wrist_ref(:,2),referencePos.Pos_wrist_ref(:,3),":", "color", '#009ADE', "LineWidth",3, "DisplayName", 'Reference wrist trajectory');
        ref_traj_tip = plot3(referencePos.Pos_tip_ref(:,1),referencePos.Pos_tip_ref(:,2),referencePos.Pos_tip_ref(:,3),":", "color", '#FFC61E', "LineWidth",3, "DisplayName", 'Reference tip trajectory');
        pause(0.5);
        ref_traj_elbow.Color(4) = 0.7; 
        ref_traj_wrist.Color(4) = 0.7;
        ref_traj_tip.Color(4) = 0.7;

        endpoint_elbow = plot3(arm.pos_elbow(1),arm.pos_elbow(2),arm.pos_elbow(3), "-", "color", '#F28522', "LineWidth",3, "DisplayName",'Actual elbow trajectory');
        endpoint_wrist = plot3(arm.pos_wrist(1),arm.pos_wrist(2),arm.pos_wrist(3), "-", "color", '#AF58BA', "LineWidth",3, "DisplayName",'Actual wrist trajectory');
        endpoint_tip = plot3(arm.pos_tip(1),arm.pos_tip(2),arm.pos_tip(3), "-", "color", '#AF58BA', "LineWidth",3, "DisplayName",'Actual wrist trajectory');
        pause(0.5);
        endpoint_elbow.Color(4) = 0.7; 
        endpoint_wrist.Color(4) = 0.7;
        endpoint_tip.Color(4) = 0.7;

        %arr_elbow = plotContact(arm.pos_elbow, 0, zeros(3,1),"#8F0038","Contact force on elbow");
        %arr_wrist = plotContact(arm.pos_wrist, 0, zeros(3,1),"#B10026","Contact force on wrist");
        
        %function arr = plotContact(start_point, f, direction, color,name)
        %    arr = quiver3(start_point(1), start_point(2), start_point(3), f*direction(1), f*direction(2), f*direction(3), "color", color, 'LineWidth', 3, 'MaxHeadSize', 2, "DisplayName",name);
        %end
        fig_simu2 = figure();
        figure(fig_simu2);
        legend('FontSize', fs, 'Location', 'northoutside', 'NumColumns',3 , 'AutoUpdate', 'off'); 
        hold off;
        
        subplot(4,2,1);
        hold on;
        yyaxis left;
        velplot1 = plot(0,0,"-", "color","#FF1F5B", "LineWidth",2, "DisplayName",'DOF 1 angular velocity');
        velplot1ref = plot(ss, trajectory(8,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","DOF 1 reference angular velocity");
        ylabel("DOF 1 angular velocity [rad/s]", 'FontSize', fs, 'Color', [0 0 0]);
        %set(gca, 'XColor', 'k', 'YColor', 'k'); 
        yyaxis right;
        hold off;
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside', 'AutoUpdate', 'off'); 
        
        subplot(4,2,2);
        hold on;
        velplot2 = plot(0,0,"-", "color","#FF1F5B", "LineWidth",2, "DisplayName",'DOF 2 angular velocity');
        velplot2ref = plot(ss, trajectory(9,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","DOF 2 reference angular velocity");
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("DOF 2 angular velocity [rad/s]", 'FontSize', fs, 'Color',[0 0 0]);
        % xlim([0, t]); 
        grid on;
        hold off;
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside', 'AutoUpdate', 'off'); 
        %set(gca, 'XColor', 'k', 'YColor', 'k'); 
        
        subplot(4,2,3);
        hold on;
        velplot3 = plot(0,0,"-", "color","#FF1F5B", "LineWidth",2, "DisplayName",'DOF 3 angular velocity');
        velplot3ref = plot(ss, trajectory(10,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","DOF 3 reference angular velocity");
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("DOF 2 angular velocity [rad/s]", 'FontSize', fs, 'Color',[0 0 0]);
        % xlim([0, t]); 
        grid on;
        hold off;
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside', 'AutoUpdate', 'off'); 
        
        subplot(4,2,4);
        hold on;
        velplot4 = plot(0,0,"-", "color","#FF1F5B", "LineWidth",2, "DisplayName",'DOF 4 angular velocity');
        velplot4ref = plot(ss, trajectory(11,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","DOF 4 reference angular velocity");
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("DOF 4 angular velocity [rad/s]", 'FontSize', fs, 'Color',[0 0 0]);
        % xlim([0, t]); 
        grid on;
        hold off;
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside', 'AutoUpdate', 'off');         

        subplot(4,2,5);
        hold on;
        velplot5 = plot(0,0,"-", "color","#FF1F5B", "LineWidth",2, "DisplayName",'DOF 5 angular velocity');
        velplot5ref = plot(ss, trajectory(12,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","DOF 5 reference angular velocity");
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("DOF 5 angular velocity [rad/s]", 'FontSize', fs, 'Color',[0 0 0]);
        % xlim([0, t]); 
        grid on;
        hold off;
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside', 'AutoUpdate', 'off'); 
        
        subplot(4,2,6);
        hold on;
        velplot6 = plot(0,0,"-", "color","#FF1F5B", "LineWidth",2, "DisplayName",'DOF 6 angular velocity');
        velplot6ref = plot(ss, trajectory(13,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","DOF 6 reference angular velocity");
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("DOF 6 angular velocity [rad/s]", 'FontSize', fs, 'Color',[0 0 0]);
        % xlim([0, t]); 
        grid on;
        hold off;
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside', 'AutoUpdate', 'off');         
 
        subplot(4,2,7);
        hold on;
        velplot7 = plot(0,0,"-", "color","#FF1F5B", "LineWidth",2, "DisplayName",'DOF 7 angular velocity');
        velplot7ref = plot(ss, trajectory(14,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","DOF 7 reference angular velocity");
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("DOF 7 angular velocity [rad/s]", 'FontSize', fs, 'Color',[0 0 0]);
        % xlim([0, t]); 
        grid on;
        hold off;
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside', 'AutoUpdate', 'off');      
        
        fig_simu3 = figure();
        figure(fig_simu3);

        Ginit = getG(model, state);
        subplot(4,2,1);
        hold on;
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("Torque [Nm]", 'FontSize', fs);
        % xlim([0, t]); 

        uplot11 = plot(0, Ginit(1),"-", 'color', "#FF1F5B", "LineWidth",2, "DisplayName","Joint 1 Robot torque");
        %uplot12 = plot(0, 0, "-", 'color', "#00CD6C", "LineWidth",2, "DisplayName","Elbow contact induced torque");
        %uplot13 = plot(0, 0, "-", 'color', "#009ADE", "LineWidth",2, "DisplayName","Wrist contact induced torque");
        uplot14 = plot(ss, DesiredTorque(1,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","Joint 1 Desired torque");
        %title("Joint 1", "FontSize",fs);
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside','NumColumns',2 , 'AutoUpdate', 'off'); %lgd = 
        grid on;
        hold off;

        subplot(4,2,2);
        hold on;
        hold on;
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("Torque [Nm]", 'FontSize', fs);
        % xlim([0, t]); 

        uplot21 = plot(0, Ginit(2),"-", 'color', "#FF1F5B", "LineWidth",2, "DisplayName","Joint 2 Robot torque");
        %uplot22 = plot(0, 0, "-", 'color', "#00CD6C", "LineWidth",2, "DisplayName","Elbow contact induced torque");
        %uplot23 = plot(0, 0, "-", 'color', "#009ADE", "LineWidth",2, "DisplayName","Wrist contact induced torque");
        uplot24 = plot(ss, DesiredTorque(2,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","Joint 2 Desired torque");
        %title("Joint 2", "FontSize",fs);
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside','NumColumns',2 , 'AutoUpdate', 'off');
        grid on;
        hold off;

        subplot(4,2,3);
        hold on;
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("Torque [Nm]", 'FontSize', fs);
        % xlim([0, t]); 

        uplot31 = plot(0, Ginit(3),"-", 'color', "#FF1F5B", "LineWidth",2, "DisplayName","Joint 3 Robot torque");
        %uplot32 = plot(0, 0, "-", 'color', "#00CD6C", "LineWidth",2, "DisplayName","Elbow contact induced torque");
        %uplot33 = plot(0, 0, "-", 'color', "#009ADE", "LineWidth",2, "DisplayName","Wrist contact induced torque");
        uplot34 = plot(ss, DesiredTorque(3,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","Joint 3 Desired torque");        
        %title("Joint 3", "FontSize",fs);
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside','NumColumns',2 , 'AutoUpdate', 'off');
        grid on;
        hold off;

        subplot(4,2,4);
        hold on;
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("Torque [Nm]", 'FontSize', fs);
        % xlim([0, t]); 

        uplot41 = plot(0, Ginit(4),"-", 'color', "#FF1F5B", "LineWidth",2, "DisplayName","Joint 4 Robot torque");
        %uplot42 = plot(0, 0, "-", 'color', "#00CD6C", "LineWidth",2, "DisplayName","Elbow contact induced torque");
        %uplot43 = plot(0, 0, "-", 'color', "#009ADE", "LineWidth",2, "DisplayName","Wrist contact induced torque");
        uplot44 = plot(ss, DesiredTorque(4,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","Joint 4 Desired torque");
        %title("Joint 4", "FontSize",fs);
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside','NumColumns',2 , 'AutoUpdate', 'off');
        grid on;
        hold off;

        subplot(4,2,5);
        hold on;
        hold on;
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("Torque [Nm]", 'FontSize', fs);
        % xlim([0, t]); 

        uplot51 = plot(0, Ginit(5),"-", 'color', "#FF1F5B", "LineWidth",2, "DisplayName","Joint 5 Robot torque");
        %uplot52 = plot(0, 0, "-", 'color', "#00CD6C", "LineWidth",2, "DisplayName","Elbow contact induced torque");
        %uplot53 = plot(0, 0, "-", 'color', "#009ADE", "LineWidth",2, "DisplayName","Wrist contact induced torque");
        uplot54 = plot(ss, DesiredTorque(5,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","Joint 5 Desired torque");
        %title("Joint 5", "FontSize",fs);
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside','NumColumns',2 , 'AutoUpdate', 'off');
        grid on;
        hold off;

        subplot(4,2,6);
        hold on;
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("Torque [Nm]", 'FontSize', fs);
        % xlim([0, t]); 

        uplot61 = plot(0, Ginit(6),"-", 'color', "#FF1F5B", "LineWidth",2, "DisplayName","Joint 6 Robot torque");
        %uplot62 = plot(0, 0, "-", 'color', "#00CD6C", "LineWidth",2, "DisplayName","Elbow contact induced torque");
        %uplot63 = plot(0, 0, "-", 'color', "#009ADE", "LineWidth",2, "DisplayName","Wrist contact induced torque");
        uplot64 = plot(ss, DesiredTorque(6,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","Joint 6 Desired torque");
        %title("Joint 6", "FontSize",fs);
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside','NumColumns',2 , 'AutoUpdate', 'off');
        grid on;
        hold off;

        subplot(4,2,7);
        hold on;
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("Torque [Nm]", 'FontSize', fs);
        % xlim([0, t]); 

        uplot71 = plot(0, Ginit(7),"-", 'color', "#FF1F5B", "LineWidth",2, "DisplayName","Joint 7 Robot torque");
        %uplot72 = plot(0, 0, "-", 'color', "#00CD6C", "LineWidth",2, "DisplayName","Elbow contact induced torque");
        %uplot73 = plot(0, 0, "-", 'color', "#009ADE", "LineWidth",2, "DisplayName","Wrist contact induced torque");
        uplot74 = plot(ss, DesiredTorque(7,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","Joint 7 Desired torque");
        %title("Joint 7", "FontSize",fs);
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside','NumColumns',2 , 'AutoUpdate', 'off');
        grid on;
        hold off;
        %set(lgd, 'Position', [0.73 0.98 0 0]);
        
        fig_simu6 = figure();
        figure(fig_simu6);
        legend('FontSize', fs, 'Location', 'northoutside', 'NumColumns',3 , 'AutoUpdate', 'off'); 
        hold off;
        
        subplot(4,2,1);
        hold on;
        yyaxis left;
        xplot1 = plot(0,trajectory(1,1),"-", "color","#FF1F5B", "LineWidth",2, "DisplayName",'DOF 1 angular position');
        xplot1ref = plot(ss, trajectory(1,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","DOF 1 reference angular position");
        ylabel("DOF 1 angular position [rad]", 'FontSize', fs, 'Color', [0 0 0]);
        %set(gca, 'XColor', 'k', 'YColor', 'k'); 
        yyaxis right;
        hold off;
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside', 'AutoUpdate', 'off'); 
        
        subplot(4,2,2);
        hold on;
        xplot2 = plot(0,trajectory(2,1),"-", "color","#FF1F5B", "LineWidth",2, "DisplayName",'DOF 2 angular position');
        xplot2ref = plot(ss, trajectory(2,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","DOF 2 reference angular position");
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("DOF 2 angular position [rad]", 'FontSize', fs, 'Color',[0 0 0]);
        % xlim([0, t]); 
        grid on;
        hold off;
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside', 'AutoUpdate', 'off'); 
        %set(gca, 'XColor', 'k', 'YColor', 'k'); 
        
        subplot(4,2,3);
        hold on;
        xplot3 = plot(0,trajectory(3,1),"-", "color","#FF1F5B", "LineWidth",2, "DisplayName",'DOF 3 angular position');
        xplot3ref = plot(ss, trajectory(3,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","DOF 3 reference angular position");
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("DOF 2 angular position [rad]", 'FontSize', fs, 'Color',[0 0 0]);
        % xlim([0, t]); 
        grid on;
        hold off;
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside', 'AutoUpdate', 'off'); 
        
        subplot(4,2,4);
        hold on;
        xplot4 = plot(0,trajectory(4,1),"-", "color","#FF1F5B", "LineWidth",2, "DisplayName",'DOF 4 angular position');
        xplot4ref = plot(ss, trajectory(4,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","DOF 4 reference angular position");
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("DOF 4 angular position [rad]", 'FontSize', fs, 'Color',[0 0 0]);
        % xlim([0, t]); 
        grid on;
        hold off;
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside', 'AutoUpdate', 'off');         

        subplot(4,2,5);
        hold on;
        xplot5 = plot(0,trajectory(5,1),"-", "color","#FF1F5B", "LineWidth",2, "DisplayName",'DOF 5 angular position');
        xplot5ref = plot(ss, trajectory(5,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","DOF 5 reference angular position");
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("DOF 5 angular position [rad]", 'FontSize', fs, 'Color',[0 0 0]);
        % xlim([0, t]); 
        grid on;
        hold off;
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside', 'AutoUpdate', 'off'); 
        
        subplot(4,2,6);
        hold on;
        xplot6 = plot(0,trajectory(6,1),"-", "color","#FF1F5B", "LineWidth",2, "DisplayName",'DOF 6 angular position');
        xplot6ref = plot(ss, trajectory(6,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","DOF 6 reference angular position");
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("DOF 6 angular position [rad]", 'FontSize', fs, 'Color',[0 0 0]);
        % xlim([0, t]); 
        grid on;
        hold off;
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside', 'AutoUpdate', 'off');         
 
        subplot(4,2,7);
        hold on;
        xplot7 = plot(0,trajectory(7,1),"-", "color","#FF1F5B", "LineWidth",2, "DisplayName",'DOF 7 angular position');
        xplot7ref = plot(ss, trajectory(7,:), "-", 'color', "#AF58BA", "LineWidth",2, "DisplayName","DOF 7 reference angular position");
        xlabel("Time [sec]", 'FontSize', fs);
        ylabel("DOF 7 aposition [rad]", 'FontSize', fs, 'Color',[0 0 0]);
        % xlim([0, t]); 
        grid on;
        hold off;
        legend('FontSize', fs, 'Orientation', 'horizontal', 'Location', 'northoutside', 'AutoUpdate', 'off');      
                
    end
   
    %[pos_desired,~,~] = getDesiredState(trajectory(1:7,1), trajectory(1:7,end), t_now_traj, t_traj);
    %err_posk_1 = pos_desired - state(1:7);
    %err_posk_2 = err_posk_1;
    %eint = zeros(7,1);

    % ind_t_elbow_old = 0;
    % ind_n_elbow_old = 0;
    % ind_t_wrist_old = 0;
    % ind_n_wrist_old = 0;
    % syms Z00;
    % eqn = 1/(-Z00+exp(0.4*Z00)) ==model.cr;
    % Z0 = vpasolve(eqn,Z00);
    kk = 0;
    %c_t_elbow = c;
    %c_n_elbow = c;
    %c_t_wrist = c;
    %c_n_wrist = c;
    ggggg = 1;
    while 1
        %t_now
        arm = kinematicAnalysis(model, state);
    
        %[ind_t_elbow, ind_n_elbow, dir_t_elbow, dir_n_elbow, dev_t_elbow, dist_n_elbow, hitendwall_elbow,~] = checkContact(arm.pos_elbow, referencePos.Pos_elbow_ref, r);
        %[ind_t_wrist, ind_n_wrist, dir_t_wrist, dir_n_wrist, dev_t_wrist, dist_n_wrist, hitendwall_wrist,~] = checkContact(arm.pos_wrist, referencePos.Pos_wrist_ref, r);
        %arm.ind_t_elbow = ind_t_elbow;
        %arm.ind_n_elbow = ind_n_elbow;
        %arm.dir_t_elbow = dir_t_elbow;
        %arm.dir_n_elbow = dir_n_elbow;

        %arm.ind_t_wrist = ind_t_wrist;
        %arm.ind_n_wrist = ind_n_wrist;
        %arm.dir_t_wrist = dir_t_wrist;
        %arm.dir_n_wrist = dir_n_wrist;

        %f_elbow = getContact(arm, r, k, c_t_elbow, c_n_elbow, ind_t_elbow, ind_n_elbow, dir_t_elbow, dir_n_elbow, dev_t_elbow, dist_n_elbow, "elbow");
        %f_wrist = getContact(arm, r, k, c_t_wrist, c_n_wrist, ind_t_wrist, ind_n_wrist, dir_t_wrist, dir_n_wrist, dev_t_wrist, dist_n_wrist, "wrist");
        
        % if stateReached(pos_desired(find(impairment==2)), state(find(impairment==2)))
        %[pos_desired,vel_desired,acc_desired] = getDesiredState(trajectory(1:7,1), trajectory(1:7,end), t_now, t_traj);
        % end
        



        
        %err_pos = pos_desired - state(1:7);
        %err_vel = vel_desired - state(8:14);

        M = getM(model,state);
        C = getC(model,state);
        G = getG(model,state);
        
        Thetauap = state (1:3);
        Thetalap = state (4:5);
        Thetawp = state (6:7);
        Thetauapdt = state (8:10);
        Thetalapdt = state (11:12);
        Thetawpdt = state (13:14);
        %Thetapdt= state(8:14);
        
        [hist,UUP,diff221,diff222,diff223,disua,disla,disw,u,Sua,Sla,Sw,prevIdxua,prevIdxla,prevIdxw,wuap,uuap,wlap,ulap,wwp,uwp,Ni_prevua,Ti_prevua,Pi_prevua,thetadrua,thetadrla,thetadrw] = Potentialfield(bbuf, abuf,hist,J_max,UUP,DesiredTorque2,t_now,ThetauapdN1, ThetauapdN2 ,ThetalapdN, ThetawpdN,thetadrua,thetadrla,thetadrw,sim,model,winsize,Pi_prevua,Ti_prevua,Ni_prevua,wuap,uuap,wlap,ulap,wwp,uwp,calua,calla,calw,ualua,ualla,ualw,M,C,G,Sua,Sla,Sw,sigma,mkkk,Thetauap,Thetalap,Thetawp,Thetauapdt,Thetalapdt,Thetawpdt,Thetauapp,Thetalapp,Thetawpp,Thetauappdt,Thetalappdt,Thetawppdt,Kua,Kla,Kw,prevIdxua,prevIdxla,prevIdxw); 
       % brc=1;
%if controller_type =="pid"
            % PID
         %   [u, ek_1_new, ek_2_new, eint_new] = pidff(model,err_pos,err_posk_1,err_posk_2,eint,ts,G);
         %   err_posk_1 = ek_1_new;
         %   err_posk_2 = ek_2_new;
         %   eint = eint_new;
        %elseif controller_type =="pd"
            % PD+gravity compensation
        %    u = Kp*err_pos + Kd*err_vel + G; % PD + gravity identification
        %elseif controller_type =="impedance"
        %    u = M*acc_desired + C + G + (M\model.MD)*(model.KD* err_vel + model.KP* err_pos);
        %elseif controller_type =="Potential"
                %end
        %if abs(u(find(impairment==1)))> torquelimit(1)
            %u(find(impairment==1)) = torquelimit(1);
        %end

        %if abs(u(find(impairment==2)))> torquelimit(2)
            %u(find(impairment==2)) = torquelimit(2);
        %end
%        u(1:3) = 0;
%        u(1) = 15;
%        u(4) = 15;
        %arm_new = kinematicAnalysis(model, state_new);
        %v_elbow = arm_new.Je * state_new(8:14);
        %v_wrist = arm_new.Jw * state_new(8:14);
        %prevIdxua=prevIdxua
        %prevIdxla=prevIdxla
        %prevIdxw=prevIdxw
        q_d= [Thetauapp(:,prevIdxua);Thetalapp(:,prevIdxla);Thetawpp(:,prevIdxw)];
        if state(1:7)== q_d
            q_d= [Thetauapp(:,prevIdxua+1);Thetalapp(:,prevIdxla+1);Thetawpp(:,prevIdxw+1)];
        end
        hu= human_effort(state(1:7), state(8:14), q_d, 'active', tau_hu, K_healthy, B_healthy, q_rest_healthy, K_spastic, B_spastic, q_rest_stroke, K_intent,impairment);
        Control_Effort = [Control_Effort, u];
        J_Vel= [J_Vel , state(8:14)];
        J_Pos= [J_Pos , state(1:7)];
        %hand_vel = norm(arm.Jw*state(8:14));
        %Hand_Vel = [Hand_Vel, hand_vel];
        %hand_acc = (hand_vel - hand_vel_old)/ts;
        %Hand_Acc = [Hand_Acc, hand_acc];
        %hand_vel_old = hand_vel;
        state_old=state;
        state_dot_old=state_dot;
        %if t_now<0.5
            %u= u+ (M+0.01*eye(size(M)))*(start_pos-end_pos)*10;
        %end
        [state_new,state_dot_new] = RKDiscretize(ts, state_old,state_dot_old, u+hu*1, 0, 0, model, arm,M,C,G); %f_elbow, f_wrist
        
            %dynamics_wrapper = @(state_old, state_dot_old) forwardDynamics(state_old, state_dot_old, u+hu*0,M,C,G);
    
    % Pass the wrapper to the integrator
        %[state_new, state_dot_new,accc] = step_semi_implicit_rk4(state_old,state_dot_old,sim.ts,dynamics_wrapper);
        %state_dot_new(1:7) = state_new(8:14);
        %[state_new, state_dot_new] = limitState(model,state_new, state_dot_new);

        state = state_new;
        %state(8:14) = state_dot_new;
        state_dot = state_dot_new;
        %state_dot(8:14) = accc;
        %Contact_Elbow = [Contact_Elbow, f_elbow];
        %Contact_Wrist = [Contact_Wrist, f_wrist];
        %Torque_Elbow = [Torque_Elbow, arm.Je'* f_elbow];
        %Torque_Wrist = [Torque_Wrist, arm.Jw'* f_wrist];
        Endpoint_elbow = [Endpoint_elbow; arm.pos_elbow];
        Endpoint_wrist = [Endpoint_wrist; arm.pos_wrist];
        Endpoint_tip = [Endpoint_tip; arm.pos_tip];
        %U = [U u];
        Time = [Time t_now];
        if debug
            if kk/plot_freq >= 1 || ggggg<2%|| hitendwall_wrist || hitendwall_elbow  || abs(state(find(impairment==2))-trajectory(find(impairment==2),end)) <=pi/180*0.1
                figure(fig_simu); 
                config = [zeros(3,1),arm.pos_elbow', arm.pos_wrist',arm.pos_tip'];
                set(armconfig,"XData", config(1,:),"YData",config(2,:),"ZData",config(3,:));
                set(endpoint_elbow,"XData", Endpoint_elbow(:,1),"YData",Endpoint_elbow(:,2),"ZData",Endpoint_elbow(:,3));
                set(endpoint_wrist,"XData", Endpoint_wrist(:,1),"YData",Endpoint_wrist(:,2),"ZData",Endpoint_wrist(:,3));
                set(endpoint_tip,"XData", Endpoint_tip(:,1),"YData",Endpoint_tip(:,2),"ZData",Endpoint_tip(:,3));
                %if norm(f_elbow) > 1e-10
                %    nor_f_elbow = 0.1*f_elbow/norm(f_elbow);
                %else
                %    nor_f_elbow = zeros(3,1);
                %end
    
                %if norm(f_wrist) > 1e-10
                %    nor_f_wrist = 0.1*f_wrist/norm(f_wrist);
                %else
                %    nor_f_wrist = zeros(3,1);
                %end
    
                %set(arr_elbow, 'XData', arm.pos_elbow(1),'YData', arm.pos_elbow(2),'ZData', arm.pos_elbow(3));%'UData', double(nor_f_elbow(1)), 'VData', double(nor_f_elbow(2)),'WData', double(nor_f_elbow(3)));
                %set(arr_wrist, 'XData', arm.pos_wrist(1),'YData', arm.pos_wrist(2),'ZData', arm.pos_wrist(3));%'UData', double(nor_f_wrist(1)), 'VData', double(nor_f_wrist(2)),'WData', double(nor_f_wrist(3)));
                %set(arr_tip, 'XData', arm.pos_tip(1),'YData', arm.pos_tip(2),'ZData', arm.pos_tip(3));
                %{
                disua=disua
                disla= disla
                disw= disw
                diff221=diff221 
                diff222=diff222
                diff223= diff223
                %}
                figure(fig_simu2);
                set(title_time, 'String', ['Time: ', num2str(t_now), ' sec'])
                set(velplot1,"XData", Time,"YData",J_Vel(1,:));
                set(velplot2,"XData", Time,"YData",J_Vel(2,:));
                set(velplot3,"XData", Time,"YData",J_Vel(3,:));
                set(velplot4,"XData", Time,"YData",J_Vel(4,:));
                set(velplot5,"XData", Time,"YData",J_Vel(5,:));
                set(velplot6,"XData", Time,"YData",J_Vel(6,:));
                set(velplot7,"XData", Time,"YData",J_Vel(7,:));
               
                figure(fig_simu3);
                set(title_time, 'String', ['Time: ', num2str(t_now), ' sec'])
                set(uplot11,"XData", Time,"YData",Control_Effort(1,:));
                set(uplot21,"XData", Time,"YData",Control_Effort(2,:));
                set(uplot31,"XData", Time,"YData",Control_Effort(3,:));
                set(uplot41,"XData", Time,"YData",Control_Effort(4,:));
                set(uplot51,"XData", Time,"YData",Control_Effort(5,:));
                set(uplot61,"XData", Time,"YData",Control_Effort(6,:));
                set(uplot71,"XData", Time,"YData",Control_Effort(7,:));

                figure(fig_simu6);
                set(title_time, 'String', ['Time: ', num2str(t_now), ' sec'])
                set(xplot1,"XData", Time,"YData",J_Pos(1,:));
                set(xplot2,"XData", Time,"YData",J_Pos(2,:));
                set(xplot3,"XData", Time,"YData",J_Pos(3,:));
                set(xplot4,"XData", Time,"YData",J_Pos(4,:));
                set(xplot5,"XData", Time,"YData",J_Pos(5,:));
                set(xplot6,"XData", Time,"YData",J_Pos(6,:));
                set(xplot7,"XData", Time,"YData",J_Pos(7,:));                
                %set(uplot12,"XData", Time,"YData",Torque_Elbow(1,:));
                %set(uplot22,"XData", Time,"YData",Torque_Elbow(2,:));
                %set(uplot32,"XData", Time,"YData",Torque_Elbow(3,:));                
                %set(uplot42,"XData", Time,"YData",Torque_Elbow(4,:));
                %set(uplot52,"XData", Time,"YData",Torque_Elbow(5,:));
                %set(uplot62,"XData", Time,"YData",Torque_Elbow(6,:));
                %set(uplot72,"XData", Time,"YData",Torque_Elbow(7,:));
                
                
                %set(uplot13,"XData", Time,"YData",Torque_Wrist(1,:));
                %set(uplot23,"XData", Time,"YData",Torque_Wrist(2,:));
                %set(uplot33,"XData", Time,"YData",Torque_Wrist(3,:));
                %set(uplot43,"XData", Time,"YData",Torque_Wrist(4,:));
                %set(uplot53,"XData", Time,"YData",Torque_Wrist(5,:));
                %set(uplot63,"XData", Time,"YData",Torque_Wrist(6,:));
                %set(uplot73,"XData", Time,"YData",Torque_Wrist(7,:));
                pause(0.005);
                kk = 0;
                
                ggggg = ggggg+1;
            end
        end

        kk = kk+1;
        if t_now > t || brc==1 %|| hitendwall_elbow || hitendwall_wrist || abs(state(find(impairment==2))-trajectory(find(impairment==2),end)) <=pi/180*0.1
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
    lw = model.lw;
    t1 = state(1);
    t2 = state(2);
    t3 = state(3);
    t4 = state(4);
    t5 = state(5);
    t6 = state(6);
    t7 = state(7);
    


    T01 = DHConvention(-pi/2, 0, 0, -pi/2 + t1); %
    T12 = DHConvention(-pi/2, 0, 0, -pi/2 + t2); %
    T02 = T01*T12;
    T23 = DHConvention(pi/2, 0, 0, t3); %
    T03 = T02*T23;
    T3_p1= DHConvention(-pi/2, 0, 0, -pi/2);
    T0_p1= T03 * T3_p1;
    Tp1_4 = DHConvention(-pi/2, lu, 0, pi/2+t4); %
    T04 = T0_p1 * Tp1_4;
    T45 = DHConvention(-pi/2, 0, 0, t5); 
    T05 = T04 * T45;
    T5_p2 = DHConvention(-pi/2, 0, 0, pi/2); 
    T0_p2 = T05 * T5_p2;
    Tp2_6 = DHConvention(0, lf, 0, t6); 
    T06 = T0_p2 * Tp2_6;
    T67 = DHConvention(-pi/2, 0, 0, t7); 
    T07 = T06 * T67;
    T78= DHConvention(0, lw, 0, 0); 
    T08=T07*T78;
    
    arm.pos_elbow = T04(1:3,4)';
    arm.pos_wrist = T06(1:3,4)';
    arm.pos_tip = T08(1:3,4)';
    arm.state = state;
    arm.R05 = T05(1:3,1:3);
    arm.Jw = calculateJacobian(T01,T02,T03,T04,T05,T06,T07,T08,6);
    arm.Je = calculateJacobian(T01,T02,T03,T04,T05,T06,T07,T08,4);
    arm.Jew = [T01(1:3,3) T02(1:3,3) T03(1:3,3) T04(1:3,3) T05(1:3,3) zeros(3,1) zeros(3,1)];
    arm.Jww = [T01(1:3,3) T02(1:3,3) T03(1:3,3) T04(1:3,3) T05(1:3,3) T06(1:3,3) T07(1:3,3)];
end

function state_new = constrainedAdmittanceController(sim, arm, model, u,~,~) %constrainedAdmittanceController trajectoty  r
    
    expp = exp(-model.vd/model.vm*sim.ts);
    A = [1, (1-expp)*model.vm/model.vd;
        0, expp];
    B = [(model.vd*sim.ts-(1-expp)*model.vm)/model.vd^2;(1-expp)/model.vd];
    state_new([1,8]) = A*arm.state([1,8]) + B*u(1);
    state_new([2,9]) = A*arm.state([2,9]) + B*u(2);
    state_new([3,10]) = A*arm.state([3,10]) + B*u(3);
    state_new([4,11]) = A*arm.state([4,11]) + B*u(4);
    state_new([5,12]) = A*arm.state([5,12]) + B*u(5);
    state_new([6,13]) = A*arm.state([6,13]) + B*u(6);
    state_new([7,14]) = A*arm.state([7,14]) + B*u(7);
    state_new = state_new';
    arm_new = kinematicAnalysis(model, state_new);
    %[ind_t_elbow, ind_n_elbow, dir_t_elbow, dir_n_elbow, ~, ~, ~,~] = checkContact(arm_new.pos_elbow, referencePos.Pos_elbow_ref, r);
    %[ind_t_wrist, ind_n_wrist, dir_t_wrist, dir_n_wrist, ~, ~, ~,~] = checkContact(arm_new.pos_wrist, referencePos.Pos_wrist_ref, r);
    v_elbow = arm_new.Je * state_new(8:14);
    v_wrist = arm_new.Jw * state_new(8:14);
    double(norm(v_elbow))
    double(norm(v_wrist))
    %if (ind_t_elbow && dot(dir_t_elbow, v_elbow) < 0) || (ind_n_elbow && dot(dir_n_elbow, v_elbow) < 0) || (ind_t_wrist && dot(dir_t_wrist, v_wrist) < 0) || (ind_n_wrist && dot(dir_n_wrist, v_wrist) < 0)
       % if ind_t_elbow && dot(dir_t_elbow, v_elbow) < 0
        %    disp("elbow")
        %    v_elbow = v_elbow-dot(dir_t_elbow, v_elbow)*dir_t_elbow;
        %end
        %if ind_n_elbow && dot(dir_n_elbow, v_elbow) < 0
        %    disp("elbow")
       %     v_elbow = v_elbow-dot(dir_n_elbow, v_elbow)*dir_n_elbow;
        %end
        %if ind_t_wrist && dot(dir_t_wrist, v_wrist) < 0
        %    disp("wrist")
        %    v_wrist = v_wrist-dot(dir_t_wrist, v_wrist)*dir_t_wrist;
        %end
        %if ind_n_wrist && dot(dir_n_wrist, v_wrist) < 0
        %    disp("wrist")
        %    v_wrist = v_wrist-dot(dir_n_wrist, v_wrist)*dir_n_wrist;
        %end
        %double(norm(v_elbow))
        %double(norm(v_wrist))
        %Jc = [arm_new.Je ;arm_new.Jw];
        %vc = [v_elbow;v_wrist];
        %q_dot_new = inv(Jc'*Jc+0.001*eye(7))*Jc'*vc;
        %state_new = double([state_new(1:7);q_dot_new]);
    %end
end



function referencePos = getTrajectoryCartesian(model, trajectory)

    lu = model.lu;
    lf = model.lf;
    lw = model.lw;
    PoseRef = [];
    PoswRef = [];
    PostRef = [];
    for i = 1:length(trajectory(1,:))
        t1 = trajectory(1,i);
        t2 = trajectory(2,i);
        t3 = trajectory(3,i);
        t4 = trajectory(4,i);
        t5 = trajectory(5,i);
        t6 = trajectory(6,i);
        t7 = trajectory(7,i);
    
        T01 = DHConvention(-pi/2, 0, 0, -pi/2 + t1); %
        T12 = DHConvention(-pi/2, 0, 0, -pi/2 + t2); %
        T02 = T01*T12;
        T23 = DHConvention(pi/2, 0, 0, t3); %
        T03 = T02*T23;
        T3_p1= DHConvention(-pi/2, 0, 0, -pi/2);
        T0_p1= T03 * T3_p1;
        Tp1_4 = DHConvention(-pi/2, lu, 0, pi/2+t4); %
        T04 = T0_p1 * Tp1_4;
        T45 = DHConvention(-pi/2, 0, 0, t5); 
        T05 = T04 * T45;
        T5_p2 = DHConvention(-pi/2, 0, 0, pi/2); 
        T0_p2 = T05 * T5_p2;
        Tp2_6 = DHConvention(0, lf, 0, t6); 
        T06 = T0_p2 * Tp2_6;
        T67 = DHConvention(-pi/2, 0, 0, t7); 
        T07 = T06 * T67;
        T78= DHConvention(0, lw, 0, 0); 
        T08=T07*T78;
        
        PoseRef = [PoseRef;T04(1:3,4)'];
        PoswRef = [PoswRef;T06(1:3,4)'];
        PostRef = [PoswRef;T08(1:3,4)'];
        
        


    end
    referencePos.Pos_elbow_ref = PoseRef;
    referencePos.Pos_wrist_ref = PoswRef;
    referencePos.Pos_tip_ref = PostRef;
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
         'FaceAlpha', 0.3, 'EdgeColor', 'none', "DisplayName",name);
end

function arr = plotContact(start_point, f, direction, color,name)
    arr = quiver3(start_point(1), start_point(2), start_point(3), f*direction(1), f*direction(2), f*direction(3), "color", color, 'LineWidth', 3, 'MaxHeadSize', 2, "DisplayName",name);
end

function [ind_t, ind_n, dir_t, dir_n, dev_t, dist_n, hitendwall,k] = checkContact(pos, Pos_ref, r)
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

    state(find((state(1:7)-model.jointlower)<0)) = model.jointlower(find((state(1:7)-model.jointlower)<0));
    state(find(state(7+find((state(1:7)-model.jointlower)<0)) <0)) = 0;
    state_dot(find(state_dot(7+find((state(1:7)-model.jointlower)<0)) <0)) = 0;

    state(find((state(1:7)-model.jointupper)>0)) = model.jointupper(find((state(1:7)-model.jointupper)>0));
    state(find(state(7+find((state(1:7)-model.jointupper)>0)) >0)) = 0;
    state_dot(find(state_dot(7+find((state(1:7)-model.jointupper)>0)) >0)) = 0;
    
    new_state = state;
    new_state_dot = zeros(14,1);
    new_state_dot(1:7) = state(8:14);
    new_state_dot(8:14) = state_dot(8:14);
end

function ind = stateReached(x_d,x)
    if norm(x_d-x)<0.1*pi/180
        ind = 1;
    else 
        ind = 0;
    end
end

function tau_human = human_effort(q, qd, q_d, mode, tau_hu, K_healthy, B_healthy, q_rest_healthy, K_spastic, B_spastic, q_rest_stroke, K_intent,impairment)

    tau_human = tau_hu;
    tau_healthy_intent = K_intent * (q_d - q);
    tau_pathology = -K_spastic .* (q - q_rest_stroke) - B_spastic .* qd;
    tau_healthy_passive = -K_healthy .* (q - q_rest_healthy) - B_healthy .* qd;
    switch mode
                 case 'passive'
            % Model 1: Human is a passive relaxed spring-damper system
            % Exoskeleton must fight human stiffness (K_h) and damping (B_h)
            tau_healthy = tau_healthy_passive; 
            
        case 'active'
            % Model 2: Healthy user or recovering patient helping the robot
            % Human applies a force directed toward the target path
            tau_healthy = tau_healthy_intent;
            
    end
 

    % Shoulder (1:3)
    for j = 1:3
        tau_human(j) = (1 - impairment.shoulder) * tau_healthy(j) + ...
                             (impairment.shoulder) * tau_pathology(j);
    end
    
    % Elbow & Forearm (4:5)
    for j = 4:5
        tau_human(j) = (1 - impairment.elbow) * tau_healthy(j) + ...
                             (impairment.elbow) * tau_pathology(j);
    end
    
    % Wrist (6:7)
    for j = 6:7
        tau_human(j) = (1 - impairment.wrist) * tau_healthy(j) + ...
                             (impairment.wrist) * tau_pathology(j);
    end

end


                                       