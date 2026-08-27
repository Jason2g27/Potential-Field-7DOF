function f = getContact(arm, r, k, ct, cn, ind_t, ind_n, dir_t, dir_n, dev_t, dist_n, which)
    f = zeros(3,1);
    if ind_n
        x = r - dist_n;
        if which == "elbow"
            x_dot = dot(arm.Je * arm.state(8:14), dir_n);
        else
            x_dot = dot(arm.Jw * arm.state(8:14), dir_n);
        end
        f = (cn*x*x_dot - k*x) * dir_n;
    end

    if ind_t
        x_t = dev_t;
        if which == "elbow"
            x_t_dot = dot(arm.Je * arm.state(8:14), dir_t);
        else
            x_t_dot = dot(arm.Jw * arm.state(8:14), dir_t);
        end
        f = f + (ct*x_t*x_t_dot - k*x_t) * dir_t;
    end
    
end
