function [a_mag, angle]= acceleration_linear(vx, vy, time)
for i = 2:length(time) - 1
    ax(i) = (vx(i+1) - vx(i-1))/(time(i+1) - time(i-1));
    ay(i) = (vy(i+1) - vy(i-1))/(time(i+1) - time(i-1));
    a_mag(i) = sqrt(ax(i)^2 + ay(i)^2);
    angle(i) = atan2(ax(i),ay(i));
end
a_mag(1) = 0;
a_mag(length(time)) = 0;
angle(1) = 0;
angle(length(time)) = 0;

end