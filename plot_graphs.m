function plot_graphs(x,y,vx, vy, omega, a_linear, a_angle, a_angular, impulse,impulse_angle, impulse_angular, time) 

for i = 1:length(vx)
    v_mag(i) = sqrt(vx(i)^2 + vy(i)^2);
    v_angle(i) = atan2(vx(i),vy(i));
end

figure
plot(x,y)
grid on

figure
subplot(2, 1, 1)
plot(time, v_mag)
title('Magnitude of Velocity')
subplot(2, 1, 2)
plot(time, v_angle*180/pi)
title('Angle of Velocity')

figure
subplot(2, 1, 1)
plot(time, a_linear)
title('Magnitude of Linear Acceleration')
subplot(2, 1, 2)
plot(time, a_angle*180/pi)
title('Angle of Linear Acceleration')

figure
subplot(2, 1, 1)
plot(time, omega)
title('Angular Velocity (rad/s)')
subplot(2, 1, 2)
plot(time, a_angular)
title('Angular Acceleration (rad/s)')

figure
subplot(2, 1, 1)
plot(time, impulse)
title('Magnitude of Impulse')
subplot(2, 1, 2)
plot(time, impulse_angle)
title('Angle of Impulse')

figure
plot(time, impulse_angular)
title('Angular Impulse')
end