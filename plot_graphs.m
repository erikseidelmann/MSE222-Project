function plot_graphs(x,y,vx, vy, omega, a_linear, a_angle, a_angular, impulse,impulse_angle, impulse_angular, time) 

for i = 1:length(vx)
    v_mag(i) = sqrt(vx(i)^2 + vy(i)^2);
    v_angle(i) = atan2(vx(i),vy(i));
end

figure
plot(x,y)
title('Position of Ball')
xlabel('X-position (m)')
ylabel('Y-position (m)')

figure
subplot(2, 1, 1)
plot(time, v_mag)
title('Magnitude of Linear Velocity')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
subplot(2, 1, 2)
plot(time, v_angle*180/pi)
title('Angle of Velocity')
xlabel('Time (s)')
ylabel('Velocity Direction (degrees from origin)')

figure
subplot(2, 1, 1)
plot(time, a_linear)
title('Magnitude of Linear Acceleration')
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')

subplot(2, 1, 2)
plot(time, a_angle*180/pi)
title('Angle of Linear Acceleration')
xlabel('Time (s)')
ylabel('Acceleration Direction (degrees from origin')

figure
subplot(2, 1, 1)
plot(time, omega)
title('Angular Velocity')
xlabel('Time (s)')
ylabel('Angular Velocity (rad/s)')
subplot(2, 1, 2)
plot(time, a_angular)
title('Angular Acceleration')
xlabel('Time (s)')
ylabel('Angular Acceleration (rad/s^2)')

figure
subplot(2, 1, 1)
plot(time, impulse)
title('Magnitude of Impulse')
xlabel('Time (s)')
ylabel('Impulse (Ns)')
subplot(2, 1, 2)
plot(time, impulse_angle)
title('Angle of Impulse')
xlabel('Time (s)')
ylabel('Impulse Direction (degrees from origin)')

figure
plot(time, impulse_angular)
title('Angular Impulse')
xlabel('Time (s)')
ylabel('Angular Impulse (Nm/s)')
end