function trajectory = trajectoryGenerationMAJ(start_pos, end_pos, t_tot, ts)
    time = 0:ts:t_tot;
    ss = time/t_tot;

    theta = start_pos + (start_pos - end_pos).*(15*ss.^4 - 6*ss.^5- 10*ss.^3);
    theta_dot = (start_pos - end_pos).*(60*ss.^3*1/t_tot - 30*ss.^4*1/t_tot- 30*ss.^2*1/t_tot);
    theta_ddot = (start_pos - end_pos).*(180*ss.^2*1/t_tot^2 - 120*ss.^3*1/t_tot^2- 60*ss*1/t_tot^2);
    trajectory = [theta;theta_dot;theta_ddot];
end
