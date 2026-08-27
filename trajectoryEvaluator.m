function [indicator, results] = trajectoryEvaluator(tt, model, group1, group2, trajectory, percentage, threshold, debug)
%  Input:
%  arm_length --- 2-dim [lu lf]
%  trajectoy --- 12xn first four rows: pos (rad); mid four rows: vel (rad/s);
%                last for rows: acc (rad/s^2)
%  debug --- debug mode

% Output:
% indicator --- whether it is a good trajectory
    impairmentLevel = [0,0,0,0]; % only must know the elbow impairment level
    impairmentLevel(group1) = 2;
    impairmentLevel(group2) = 1;
    % for the wrist
    [TangentialTorque1, NormalTorque1, DesiredTorque1] = parseTrajectory(model, trajectory, impairmentLevel);
    [indicator_torque1, percentage1] = checkInteractiveTorque(group1, group2, TangentialTorque1, NormalTorque1, DesiredTorque1, percentage, debug);
    [indicator_transmission1, ratio1] = checkTransmission(group1, group2, TangentialTorque1, NormalTorque1, threshold, debug);
    
    
    % for the wrist
    impairmentLevel(group1) = 1;
    impairmentLevel(group2) = 2;
    [TangentialTorque2, NormalTorque2, DesiredTorque2] = parseTrajectory(model, trajectory, impairmentLevel);
    [indicator_torque2, percentage2] = checkInteractiveTorque(group2, group1, TangentialTorque2, NormalTorque2, DesiredTorque2, percentage, debug);
    [indicator_transmission2, ratio2] = checkTransmission(group2, group1, TangentialTorque2, NormalTorque2, threshold, debug);
    results = [percentage1,percentage2,ratio1,ratio2];

    if debug
        plotTorque(tt, TangentialTorque1, NormalTorque1, DesiredTorque1, "Wrist contact: Group1 is mildly impaired")
        plotTorque(tt, TangentialTorque2, NormalTorque2, DesiredTorque2, "Wrist contact: Group2 is mildly impaired")
    end
    
    indicator = indicator_torque1 && indicator_torque2 && indicator_transmission1 && indicator_transmission2;
end


function [TangentialTorque, NormalTorque, DesiredTorque] = parseTrajectory(model, trajectory, impairmentLevel)

    TangentialTorque = zeros(7,length(trajectory(1,:)));
    NormalTorque = zeros(7,length(trajectory(1,:)));
    DesiredTorque = zeros(7,length(trajectory(1,:)));
    for i = 1:1:length(trajectory(1,:))
        % set up
        state.state = trajectory(:,i);

        

        
        state.T01 = DHConvention(-PI/2, 0, 0, -PI/2 + state.state(1)); %
        state.T12 = DHConvention(-PI/2, 0, 0, -PI/2 + state.state(2)); %
        state.T23 = DHConvention(PI/2, 0, 0, state.state(3)); %
        state.T3_p1 = DHConvention(-PI/2, 0, 0, -PI/2); %
        state.Tp1_4 = DHConvention(-PI/2, model.lu, 0, PI/2+state.state(4)); %
        state.T45 = DHConvention(-PI/2, 0, 0, state.state(5));
        state.T5_p2 = DHConvention(-PI/2, 0, 0, PI/2); %
        state.Tp2_6 = DHConvention(0, model.lf, 0, state.state(4));       
        state.T67 = DHConvention(-PI/2, 0, 0, state.state(7)); 
        state.T78= DHConvention(0, model.lw, 0, 0); 
        
        M = getM(model,state.state);
        C = getC(model,state.state);
        G = getG(model,state.state);
        DesiredTorque(:,i) = inverseDynamics(state.state, M, C, G);
        trueTorque = zeros(7,1);
        ind = find(impairmentLevel==2);
        trueTorque(ind) = DesiredTorque(ind,i);

        [tangentialTorque, normalTorque] = getInteractionInducedTorque(model, state, trueTorque);
        TangentialTorque(:,i) = tangentialTorque;
        NormalTorque(:,i) = normalTorque;
        
    end

end



function [pass, ratio] = checkInteractiveTorque(group1, group2, TangentialTorque, NormalTorque, DesiredTorque, percentage, debug)
    % group1 is mildly impaired
    % group2 is severely impaired
    pass = 1;
    productTan = TangentialTorque .* DesiredTorque;
    productNor = NormalTorque .* DesiredTorque;
    range = length(TangentialTorque(1,:));
    ratioTan1 = sum(productTan(group1, :)<-1e-6,2)/range;
    ratioTan2 = sum(productTan(group2, :)>1e-6,2)/range;

    ratioNor1 = sum(productNor(group1, :)<-1e-6,2)/range;
    ratioNor2 = sum(productNor(group2, :)>1e-6,2)/range;
    ratio = [ratioTan1', ratioNor1', ratioTan2', ratioNor2'];
    if debug
        disp("For mildly impaired joint ");
        disp(group1);
        disp("ratioTan1: ");
        disp(ratioTan1);
        disp(" ratioTan2: ");
        disp(ratioTan2);
        disp("ratioNor1: ");
        disp(ratioNor1);
        disp(" ratioNor2: ");
        disp(ratioNor2);
    end
    % check if the normal torque is always around zero

    if any([ratioTan1', ratioTan2'] < percentage(1))
        pass = 0;
    end

    for i = 1:length(group1)
        if any(ratioNor1' < percentage(2)) && max(abs(NormalTorque(group1(i),:))) > 1e-6
            pass = 0;
            break;
        end
    end

    for i = 1:length(group2)
        if any(ratioNor2' < percentage(2)) && max(abs(NormalTorque(group2(i),:))) > 1e-6
            pass = 0;
            break;
        end
    end
    
end


function [pass, Ratio] = checkTransmission(group1, group2, TangentialTorque, NormalTorque, threshold, debug)
    % group1 is mildly impaired
    % group2 is severely impaired
    pass = 1;
    Ratio = [];
    for i = 1:length(group1)
        for j = 1:length(group2)
            ratioTan = abs(TangentialTorque(group2(j), :))./abs(TangentialTorque(group1(i), :));
            ratioNor = abs(NormalTorque(group2(j), :))./abs(NormalTorque(group1(i), :));
            
            filteredRatioTan = ratioTan(~isnan(ratioTan) & ~isinf(ratioTan));
            filteredRatioNor = ratioNor(~isnan(ratioNor) & ~isinf(ratioNor));
            averageRatioTan = mean(filteredRatioTan);
            averageRatioNor = mean(filteredRatioNor);
            % not likely to be inf or nan
            if averageRatioTan < threshold(1) || averageRatioTan > threshold(2)
                pass = 0;
            end
            if max(abs(NormalTorque(group1(i), :))) > 1e-6
                if averageRatioNor < threshold(3) || averageRatioNor > threshold(4)
                    pass = 0;
                end
            end

            Ratio = [Ratio averageRatioTan averageRatioNor];
        end
    end
    if debug
        disp("the transmission ratio is ");
        disp(Ratio);
    end
end


function plotTorque(t, TangentialTorque, NormalTorque, DesiredTorque, text)
    fontsize = 16;
    legendfontsize = 12;

    fig = figure();
    subplot(7,1,1)
    yyaxis left;
    plot(t, TangentialTorque(1,:), "-", "color", "#FF1F5B", 'LineWidth', 2);
    hold on
    plot(t, NormalTorque(1,:), "-", "color", "#00CD6C", 'LineWidth', 2);
    hold off
    yyaxis right;
    plot(t, DesiredTorque(1,:), "-", "color", "#009ADE", 'LineWidth', 2);
    title("Joint 1", 'FontSize', fontsize)
    grid on

    subplot(7,1,2)
    yyaxis left;
    plot(t, TangentialTorque(2,:), "-", "color", "#FF1F5B", 'LineWidth', 2);
    hold on
    plot(t, NormalTorque(2,:), "-", "color", "#00CD6C", 'LineWidth', 2);
    hold off
    yyaxis right;
    plot(t, DesiredTorque(2,:), "-", "color", "#009ADE", 'LineWidth', 2);
    title("Joint 2", 'FontSize', fontsize)
    grid on

    subplot(7,1,3)
    yyaxis left;
    plot(t, TangentialTorque(3,:), "-", "color", "#FF1F5B", 'LineWidth', 2);
    hold on
    plot(t, NormalTorque(3,:), "-", "color", "#00CD6C", 'LineWidth', 2);
    hold off
    yyaxis right;
    plot(t, DesiredTorque(3,:), "-", "color", "#009ADE", 'LineWidth', 2);
    title("Joint 3", 'FontSize', fontsize)
    grid on

    subplot(7,1,4)
    yyaxis left;
    plot(t, TangentialTorque(4,:), "-", "color", "#FF1F5B", 'LineWidth', 2, "DisplayName", "Tangential torque");
    hold on
    plot(t, NormalTorque(4,:), "-", "color", "#00CD6C", 'LineWidth', 2, "DisplayName", "Normal torque");
    hold off
    yyaxis right;
    plot(t, DesiredTorque(4,:), "-", "color", "#009ADE", 'LineWidth', 2, "DisplayName", "Desired torque");
    title("Joint 4", 'FontSize', fontsize)
    grid on
    
    subplot(7,1,5)
    yyaxis left;
    plot(t, TangentialTorque(5,:), "-", "color", "#FF1F5B", 'LineWidth', 2);
    hold on
    plot(t, NormalTorque(5,:), "-", "color", "#00CD6C", 'LineWidth', 2);
    hold off
    yyaxis right;
    plot(t, DesiredTorque(5,:), "-", "color", "#009ADE", 'LineWidth', 2);
    title("Joint 2", 'FontSize', fontsize)
    grid on

    subplot(7,1,6)
    yyaxis left;
    plot(t, TangentialTorque(6,:), "-", "color", "#FF1F5B", 'LineWidth', 2);
    hold on
    plot(t, NormalTorque(6,:), "-", "color", "#00CD6C", 'LineWidth', 2);
    hold off
    yyaxis right;
    plot(t, DesiredTorque(6,:), "-", "color", "#009ADE", 'LineWidth', 2);
    title("Joint 3", 'FontSize', fontsize)
    grid on

    subplot(7,1,7)
    yyaxis left;
    plot(t, TangentialTorque(7,:), "-", "color", "#FF1F5B", 'LineWidth', 2, "DisplayName", "Tangential torque");
    hold on
    plot(t, NormalTorque(7,:), "-", "color", "#00CD6C", 'LineWidth', 2, "DisplayName", "Normal torque");
    hold off
    yyaxis right;
    plot(t, DesiredTorque(7,:), "-", "color", "#009ADE", 'LineWidth', 2, "DisplayName", "Desired torque");
    title("Joint 4", 'FontSize', fontsize)
    grid on

    % add legend
    Lgnd = legend('show', 'Orientation', 'horizontal', 'FontSize', legendfontsize);
    Lgnd.Position(1) = 0.1;
    Lgnd.Position(2) = 0.93;
    han=axes(fig,'visible','off');
    han.XLabel.Visible='on';
    yyaxis left;
    han.YLabel.Visible='on';
    ylabel(han,'Interaction-induced torque [Nm]', 'FontSize', fontsize);
    yyaxis right;
    han.YLabel.Visible='on';
    ylabel(han,'Desired torque [Nm]', 'FontSize', fontsize);
    xlabel(han,"Time [s]", 'FontSize', fontsize);
    set(fig, 'Name', text, 'NumberTitle', 'off');

    

end