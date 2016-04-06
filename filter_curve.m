function [t,x,v] = filter_curve(t_unf,x_unf,x_max,v_unf)
    x(1,:) = x_unf(1,:);
    t(1) = t_unf(1);
    v(1) = v_unf(1);
    i = 2;
    while x(i-1,1) <= x_max
        x(i,:) = x_unf(i,:);
        t(i,1) = t_unf(i);
        v(i) = v_unf(i);
        i = i+1;
    end
end