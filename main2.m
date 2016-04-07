clear all
close all
clc
data = zeros(1,9);%[time, x_ball, y_ball, vx_ball, vy_ball, omega, ax_ball, ay_ball, alpha];

implulses = zeros(1,4);%[x_impulse, y_impulse, angular_impulse, time];

% Initialize variables
% Co-ordinate system (0,0) at bottom left corner
%Initial positions of where the ball is still in contact with the spring
%Initial height of ball with spring
%theta_bb = zeros(1,2920)
x_impulse(1) = 0;
y_impulse(1)=0;
angular_impulse(1) = 0;
y_top = 0.30;%height of the board
x_right = 0.30; %width of the board
m_ball= 0.006; % Mass of Ball
r_ball = 0.00815; %Radius of Ball
I_ball = 2/5*m_ball*r_ball^2; %Moment of inertia of ball
g = 9.81;  %gravity of the Earth at sea level
time(1) = 0;%:0.05:4; %Hardcoded time of which ball is still in contact with the spring
x_ball(1) = 0.016; %Initial x position of ball at t= t1
y_ball(1) = y_top - 0.02;%Initial x position of ball at t= t1
vx_ball(1) = 0;%Initial x position of ball at t= t1
vy_ball(1) = 0;%Initial y position of ball at t= t1
omega(1) = 0;%Initial angular v of ball at t= t1
ax_ball(1) = 0;%Initial x acceleration of ball at t= t1
ay_ball(1) = 0;%Initial y acceleration of ball at t= t1
alpha(1) = 0;%Initial alpha of ball at t= t1
%Matrix of data containing all of the ball values. When we call this matrix
%it updates the values at that certain time
data = [time, x_ball, y_ball, vx_ball, vy_ball, omega, ax_ball, ay_ball, alpha];

%% Spring Release (Ball slipping, not rolling)
% Ball is released through spring
%Solves the time that the spring and ball are in contact
[t_unf,x_unf] = solve_spring(m_ball);
%Filter allows us to solve the spring as a function of x given time
[time,x_spring] = filter_spring(t_unf,x_unf,0);
% x_spring = (displacement,velocity)
%Update ball parameters
a_spring = diff(x_spring(:,2));
for i=2:length(x_spring)
    %Position of the spring initial
    x_ball(i,1) = x_ball(1) + 0.024 + x_spring(i,1);
    y_ball(i,1) = y_ball(1);
    %vx is the velocity of the spring
    vx_ball(i,1) = x_spring(i,2);
    %vx is the velocity of the spring
    vy_ball(i,1) = 0;
    %angular velcoty is 0 while spring is applied to ball
    omega(i,1) = 0;
    ax_ball(i,1) = a_spring(i-1);
    ay_ball(i,1) = 0;
    alpha(i,1) = 0;
end
%% Ball leaves spring into track
% Assume sudden no-slip
%Matrix that stores all of the data of the ball at that certain time.
%Updates with new values as soon as 'data' is called. 
data = [time, x_ball, y_ball, vx_ball, vy_ball, omega, ax_ball, ay_ball, alpha];
%plot(time,x(:,1));
% Energy Equation
% T_1 = T_200
% 1/2mv_1^2 = 1/2mv_2^2 + 1/2I_gOmega^2
%solve for velcity and omega
A = [1 -r_ball^2;
    1/2*m_ball 1/2*I_ball];
b = [0; 1/2*m_ball*(vx_ball(i))^2];
%x(1) is velocity squared, x(2) is omega squared
x = A\b;
vx_ball(i+1) = sqrt(x(1));
vy_ball(i+1) = 0;
omega(i+1) = -sqrt(x(2));
%the time after the spring is 0.03s
time(i+1) = time(i) + 0.001;
%transition to next part of the track
x_ball(i+1) = x_ball(i);% + (time(i+1) - time(i))*vx_ball(i+1);
y_ball(i+1) = y_ball(i);
%ax_ball(i+1) = (vx_ball(i+1)-vx_ball(i))/(time(i+1) - time(i));
j = i+2;
%solve for straight track up until x-ball is = 0.169 which is the location
%of the halfpipe
while x_ball(j-1) <= 0.169
    time(j) = time(j-1) + 0.001;
    x_ball(j) = x_ball(j-1) + 0.001*vx_ball(i+1); %(time(j) - time(i+1))
    y_ball(j) = y_ball(i+1);
    vx_ball(j) = vx_ball(i+1);
    vy_ball(j) = vy_ball(i+1);
    omega(j) = omega(i+1);
    % alpha(j) = 0;
    % ax_ball(j) = 0;
    % ay_ball(j) = 0;
    j = j+1;
end
%Matrix to save the changing data of time, x, y position and velocity,and
%angular velocity
data = [time, x_ball, y_ball, vx_ball, vy_ball, omega];
index_curve_enter = j-1
% comet(x_ball,y_ball)
% xlim([0,x_right])
% ylim([0,y_top])

%% Ball Enters and Travels through curve
% Assume Omega changes direction instantaneously
omega(index_curve_enter) = -omega(index_curve_enter);
%diameter of the half pipe and radius
d_curve = 0.1114;
r_curve = d_curve/2;
%Solve curve ODE
[t_unf,theta_unf,omega_curve_unf] = solve_curve(r_ball,r_curve,omega(index_curve_enter));
%Filter curve ODE
[time_new,theta,omega_curve] = filter_curve(t_unf,theta_unf,pi,omega_curve_unf);
%Increase time array
time_new_merge = time(index_curve_enter) + time_new(2:length(time_new));
time = vertcat(time,time_new_merge);
%Update values
omega_ball_curve = (r_curve-r_ball)/r_ball*omega_curve;
for k = 2:length(theta)%index_curve_enter+1:length(time)
          %before the ball has reached 90 degrees ie. before it is perfectly
          %vertical
    if theta(k) < pi/2
          %after the ball has exceed 90 degrees ie. not perfectly vertical
        x_ball = vertcat(x_ball,x_ball(index_curve_enter)+(r_curve-r_ball)*sin(theta(k)));
    else
        x_ball = vertcat(x_ball,x_ball(index_curve_enter)-(r_curve-r_ball)*sin(-theta(k)));
    end
    y_ball = vertcat(y_ball,y_ball(index_curve_enter)-2*(r_curve-r_ball)+(r_curve-r_ball)*(1+cos(theta(k))));
    vx_ball = vertcat(vx_ball,omega_curve(k)*r_curve*sin(pi/2-theta(k)));
    vy_ball = vertcat(vy_ball,-omega_curve(k)*r_curve*cos(pi/2-theta(k)));
    omega = vertcat(omega, omega_ball_curve(k));
    ax_ball = vertcat(ax_ball,0);
    ay_ball = vertcat(ay_ball,0);
    alpha = vertcat(alpha,0);
end
[time,x_ball,y_ball];
%[vx_ball,vy_ball,omega]
%plot(time,vx_ball,time,vy_ball)

%% Ball Leaves curve (Yay!)
%Length of index of leave-curve is equal to the number of time intervals
%added into 'data'
index_leave_curve = length(time)
%first value of vx ball is the same as the one before 
vx_ball(index_leave_curve+1) = vx_ball(index_leave_curve);
%first value of omega ball is the same as the one before 
vy_ball(index_leave_curve+1) = 0;
%first value of the time of the ball is the same as the one before 
omega(index_leave_curve+1) = omega(index_leave_curve);
%first value of position of ball is the same as the one before 
time(index_leave_curve+1) = time(index_leave_curve) + 0.001;
x_ball(index_leave_curve+1) = x_ball(index_leave_curve) + (time(index_leave_curve+1) - time(index_leave_curve))*vx_ball(index_leave_curve);
y_ball(index_leave_curve+1) = y_ball(index_leave_curve);
l = index_leave_curve+2;
while x_ball(l-1) > 0.156 %0.17-0.018 = 0.152 is ideal
    time(l) = time(l-1) + 0.001;
    x_ball(l) = x_ball(l-1) + 0.001*vx_ball(index_leave_curve+1);
    y_ball(l) = y_ball(index_leave_curve+1);
    %Vx is the same for the entire flat section 
    vx_ball(l) = vx_ball(index_leave_curve+1);
    %Vy is the same for the entire flat section 
    vy_ball(l) = vy_ball(index_leave_curve+1);
    %Omega is the same for the entire flat section 
    omega(l) = omega(index_leave_curve+1);
    % alpha(l) = 0;
    % ax_ball(l) = 0;
    % ay_ball(l) = 0;
    l = l+1;
end
%Calling data to update matrix with new values of ball
data = [time, x_ball, y_ball, vx_ball, vy_ball, omega];
% figure
% comet(x_ball,y_ball)
% xlim([0,x_right])
% ylim([0,y_top])

%% Ball Enters Ramp
%Index is the same as the number of time values that have been inputted
%into 'data'
index_enters_ramp = length(time);
%Since we have no friction we are declaring a drag factor of 0.8
%decrease magnitude of the velocity upon enting the ramp by a factor of 0.8
%to account for drag/friction
vx_ball(index_enters_ramp+1) = 0.8*vx_ball(index_enters_ramp)*cosd(19.1);
vy_ball(index_enters_ramp+1) = -0.8*vx_ball(index_enters_ramp)*sind(19.1);
%omega remains the same 


v = sqrt((vx_ball(index_enters_ramp+1))^2 + (vy_ball(index_enters_ramp+1))^2);

omega(index_enters_ramp+1) = v*r_ball;
%time is increased by 1 ms
time(index_enters_ramp+1) = time(index_enters_ramp) + 0.001;


x_impulse (index_enters_ramp+1,1) = m_ball*(vx_ball(index_enters_ramp+1)-vx_ball(index_enters_ramp));
y_impulse(index_enters_ramp+1,1)= m_ball*(vy_ball(index_enters_ramp+1)-vy_ball(index_enters_ramp));
angular_impulse(index_enters_ramp+1,1) = I_ball*(omega(index_enters_ramp+1) - omega(index_enters_ramp));

%implulses = [x_impulse, y_impulse, angular_impulse, time];




implulses = [x_impulse, y_impulse, angular_impulse, time];
%position of the ball remains the same for transitino between flat track
%and ramp
x_ball(index_enters_ramp+1) = x_ball(index_enters_ramp);% + 0.001*vx_ball(index_enters_ramp);
y_ball(index_enters_ramp+1) = y_ball(index_enters_ramp);
l = index_enters_ramp+2;
%calculate acceleration while on the ramp
a_x = g*sind(19.1)*cosd(19.1);
a_y = -g*sind(19.1)*sind(19.1);
t_step = 0.001;
while x_ball(l-1) > 0.0985 %0.17-0.018-56.61cos19.1 is ideal
    time(l) = time(l-1) + t_step;
    %use x = x0 + v0t + 0.5at^2 and v = v0 + at to calculate position and
    %velocity
    vx_ball(l) = vx_ball(l-1) + t_step*a_x;
    x_ball(l) = x_ball(l-1) + t_step*vx_ball(l-1) + 1/2*t_step^2*a_x;
    vy_ball(l) = vy_ball(l-1) + t_step*a_y;
    y_ball(l) = y_ball(l-1) + t_step*vy_ball(l-1) + 1/2*t_step^2*a_y;
    %omega = v/r ie. NO slip condition
    omega(l) = sqrt(vx_ball(l)^2 + vy_ball(l)^2)/r_ball; % Non-slip condition
    l = l+1;
end
%Calling data to update values of ball
data = [time, x_ball, y_ball, vx_ball, vy_ball, omega]
% comet(x_ball,y_ball)
% xlim([0,x_right])
% ylim([0,y_top])
%theta
%time = time + time_new
%plot(time_new,theta)
%% Ball embarks on magnificent flight
index_projectile_start = length(time);
%Since we have no friction we are declaring a drag factor of 0.6
%decrease magnitude of the velocity upon enting the ramp by a factor of 0.6
%to account for drag/friction
vx_ball(index_projectile_start+1) = 0.6*vx_ball(index_projectile_start);
vy_ball(index_projectile_start+1) = 0.6*vy_ball(index_projectile_start);

v = sqrt((vx_ball(index_projectile_start+1))^2 + (vy_ball(index_projectile_start+1))^2);

omega(index_projectile_start+1) = v*r_ball;
time(index_projectile_start+1) = time(index_projectile_start) + 0.001;


x_impulse (index_projectile_start+1,1) = m_ball*(vx_ball(index_projectile_start+1)-vx_ball(index_projectile_start));
y_impulse(index_projectile_start+1,1)= m_ball*(vy_ball(index_projectile_start+1)-vy_ball(index_projectile_start));;
angular_impulse(index_projectile_start+1,1) = I_ball*(omega(index_projectile_start+1) - omega(index_projectile_start));
time(index_projectile_start+1) = time(index_projectile_start) + 0.001;



x_ball(index_projectile_start+1) = x_ball(index_projectile_start);% + 0.001*vx_ball(index_projectile_start);
y_ball(index_projectile_start+1) = y_ball(index_projectile_start);
l = index_projectile_start+2;
%Zero acceleration in the x direction for projectile motion
a_x = 0;
%gravity is -9.8
a_y = -g;
t_step = 0.001;
% syms t
% t_max_projectile = vpasolve(-1/2*g*t^2+(vx_ball(index_projectile_start)*28/53+vy_ball(index_projectile_start))*t + 0.0680377,t)
%ramp. Solution, 0.093426921...used in to end the while loop and exit
%projectile motion. 
while time(l-1)-time(index_projectile_start) <  1.15*0.093426921662130394418769034451829 % t_max_projectile(2) %0.17-0.018-56.61cos19.1 is ideal
    time(l) = time(l-1) + t_step;
    %use x = x0 + v0t + 0.5at^2 and v = v0 + at to calculate position and
    %velocity
    vx_ball(l) = vx_ball(l-1) + t_step*a_x;
    x_ball(l) = x_ball(l-1) + t_step*vx_ball(l-1) + 1/2*t_step^2*a_x;
    vy_ball(l) = vy_ball(l-1) + t_step*a_y;
    %omega = v/r
    y_ball(l) = y_ball(l-1) + t_step*vy_ball(l-1) + 1/2*t_step^2*a_y;
    omega(l) = omega(index_projectile_start); % Non-slip condition
    l = l+1;
end

%% Ball meets penultimate doom
% Perfectly inelastic collision ie. Stuck together for the duration of the
% rotation
index_metal_start = length(time);
vx_ball(index_metal_start+1) = 0;
vy_ball(index_metal_start+1) = 0;
time(index_metal_start+1) = time(index_metal_start) + 0.001;

v = sqrt((vx_ball(index_metal_start+1))^2 + (vy_ball(index_metal_start+1))^2);

omega(index_metal_start+1) = v*r_ball;

x_impulse (index_metal_start+1,1) = m_ball*(vx_ball(index_metal_start+1)-vx_ball(index_metal_start));
y_impulse(index_metal_start+1,1)= m_ball*(vy_ball(index_metal_start+1)-vy_ball(index_metal_start));
angular_impulse(index_metal_start+1,1) = I_ball*(omega(index_metal_start+1)-omega(index_metal_start));
time(index_metal_start+1) = time(index_metal_start) + 0.001;


x_ball(index_metal_start+1) = x_ball(index_metal_start);
y_ball(index_metal_start+1) = y_ball(index_metal_start);
l = index_metal_start+2;
%0.8 is a drag factor
a_x = 0.8*g*sind(25)*cosd(25);
a_y = -0.8*g*sind(25)*sind(25);
t_step = 0.001;
while x_ball(l-1) < 0.0791
    time(l) = time(l-1) + t_step;
    vx_ball(l) = vx_ball(l-1) + t_step*a_x;
    x_ball(l) = x_ball(l-1) + t_step*vx_ball(l-1) + 1/2*t_step^2*a_x;
    vy_ball(l) = vy_ball(l-1) + t_step*a_y;
    y_ball(l) = y_ball(l-1) + t_step*vy_ball(l-1) + 1/2*t_step^2*a_y;
    %solving for omega no slip conditon
    omega(l) = -sqrt(vx_ball(l)^2 + vy_ball(l)^2)/r_ball; 
    l = l+1;
end
data = [time, x_ball, y_ball, vx_ball, vy_ball, omega];

%% Ball finds hope in darkness
% Straight slope
index_pipe_top_start = length(time);
omega(index_pipe_top_start+1) = 5/7*omega(index_pipe_top_start);
vx_ball(index_pipe_top_start+1) = -omega(index_pipe_top_start+1)*r_ball;
vy_ball(index_pipe_top_start+1) = 0;
time(index_pipe_top_start+1) = time(index_pipe_top_start) + 0.001;
x_ball(index_pipe_top_start+1) = x_ball(index_pipe_top_start);
y_ball(index_pipe_top_start+1) = y_ball(index_pipe_top_start);
l = index_pipe_top_start+2;
a_x = -0.3;
a_y = 0;
t_step = 0.001;
while x_ball(l-1) < 0.0791+0.195
    time(l) = time(l-1) + t_step;
    vx_ball(l) = vx_ball(l-1) + t_step*a_x;
    x_ball(l) = x_ball(l-1) + t_step*vx_ball(l-1) + 1/2*t_step^2*a_x;
    vy_ball(l) = 0;
    y_ball(l) = y_ball(l-1);
    omega(l) = -sqrt(vx_ball(l)^2 + vy_ball(l)^2)/r_ball; % Non-slip condition
    l = l+1;
end
data = [time, x_ball, y_ball, vx_ball, vy_ball, omega];

%% Another failed attempt at flight
% Assume projectile motion
index_proj_start = length(time);
vx_ball(index_proj_start+1) = vx_ball(index_proj_start);
vy_ball(index_proj_start+1) = 0;
time(index_proj_start+1) = time(index_proj_start) + 0.001;

v = sqrt((vx_ball(index_proj_start+1))^2 + (vy_ball(index_proj_start+1))^2);
omega(index_proj_start+1) = r_ball*v;

x_impulse (index_proj_start+1,1) = m_ball*(vx_ball(index_proj_start+1)-vx_ball(index_proj_start));
y_impulse(index_proj_start+1,1)= m_ball*(vy_ball(index_proj_start+1)-vy_ball(index_proj_start));
angular_impulse(index_proj_start+1,1) = I_ball*(omega(index_proj_start+1)-omega(index_proj_start));
time(index_proj_start+1) = time(index_proj_start) + 0.001;


x_ball(index_proj_start+1) = x_ball(index_proj_start);
y_ball(index_proj_start+1) = y_ball(index_proj_start);
l = index_proj_start+2;
a_x = 0;
a_y = -g;
t_step = 0.001;
% syms t2
% t_max_projectile_2 = vpasolve(-1/2*g*t2^2-255/441*vx_ball(index_proj_start)*t2+0.0255,t2)

while time(l-1)-time(index_proj_start) <  1.2*0.051427348524224903787683960221078 %t_max_projectile_2(2)
    time(l) = time(l-1) + t_step;
    vx_ball(l) = vx_ball(l-1) + t_step*a_x;
    x_ball(l) = x_ball(l-1) + t_step*vx_ball(l-1) + 1/2*t_step^2*a_x;
    vy_ball(l) = vy_ball(l-1) + t_step*a_y;
    y_ball(l) = y_ball(l-1) + t_step*vy_ball(l-1) + 1/2*t_step^2*a_y;
    omega(l) = omega(index_proj_start);
    % Non-slip condition
    l = l+1;
end

%% Ramp into the jaws of death
% Assume Straight slope
index_slope_death_start = length(time);
vx_ball(index_slope_death_start+1) = 0; %omega(index_slope_death_start)*r_ball;
vy_ball(index_slope_death_start+1) = 0;
omega(index_slope_death_start+1) = 0; %omega(index_slope_death_start);
time(index_slope_death_start+1) = time(index_slope_death_start) + 0.001;
x_ball(index_slope_death_start+1) = x_ball(index_slope_death_start);
y_ball(index_slope_death_start+1) = y_ball(index_slope_death_start);
l = index_slope_death_start+2;
a_x = -0.8*g*sind(30)*cosd(30); %0.8 is a drag factor
a_y = -0.8*g*sind(30)*sind(30);
t_step = 0.001;
while x_ball(l-1) > 0.0791+0.195
    time(l) = time(l-1) + t_step;
    vx_ball(l) = vx_ball(l-1) + t_step*a_x;
    x_ball(l) = x_ball(l-1) + t_step*vx_ball(l-1) + 1/2*t_step^2*a_x;
    vy_ball(l) = vy_ball(l-1) + t_step*a_y;
    y_ball(l) = y_ball(l-1) + t_step*vy_ball(l-1) + 1/2*t_step^2*a_y;
    omega(l) = sqrt(vx_ball(l)^2 + vy_ball(l)^2)/r_ball; % Non-slip condition
    l = l+1;
end

%% Ball valliantly slaughters cerberus while crossing the river styx
% Straight slope
index_pipe_bottom_start = length(time)
vx_ball(index_pipe_bottom_start+1) = -omega(index_pipe_bottom_start)*r_ball;
vy_ball(index_pipe_bottom_start+1) = 0;
time(index_pipe_bottom_start+1) = time(index_pipe_bottom_start) + 0.001;

v = sqrt((vx_ball(index_pipe_bottom_start+1))^2 + (vy_ball(index_pipe_bottom_start+1))^2);
omega(index_pipe_bottom_start+1) = v*r_ball;

x_impulse (index_pipe_bottom_start+1,1) = m_ball*(vx_ball(index_pipe_bottom_start+1)-vx_ball(index_pipe_bottom_start));
y_impulse(index_pipe_bottom_start+1,1)= m_ball*(vy_ball(index_pipe_bottom_start+1)-vy_ball(index_pipe_bottom_start));
angular_impulse(index_pipe_bottom_start+1,1) = I_ball*(omega(index_pipe_bottom_start+1) - omega(index_pipe_bottom_start));


x_ball(index_pipe_bottom_start+1) = x_ball(index_pipe_bottom_start);
y_ball(index_pipe_bottom_start+1) = y_ball(index_pipe_bottom_start);
l = index_pipe_bottom_start+2;
a_x = 0.1;
a_y = 0;
t_step = 0.001;
while x_ball(l-1) > 0.028
    time(l) = time(l-1) + t_step;
    vx_ball(l) = vx_ball(l-1) + t_step*a_x;
    x_ball(l) = x_ball(l-1) + t_step*vx_ball(l-1) + 1/2*t_step^2*a_x;
    vy_ball(l) = 0;
    y_ball(l) = y_ball(l-1);
    omega(l) = sqrt(vx_ball(l)^2 + vy_ball(l)^2)/r_ball; % Non-slip condition
    l = l+1;
end

%% Ball takes break from bounty hunting and visits local amusement park
% Introducing 'Big Boy'
omega_bb = zeros(length(time)+1,1);
theta_bb = zeros(length(time)+1,1);
alpha_bb = -45.900;

index_bb_start = length(time);
vx_ball(index_bb_start+1) = vx_ball(index_bb_start); %omega(index_bb_start)*r_ball;
vy_ball(index_bb_start+1) = 0;
omega(index_bb_start+1) = 0; %omega(index_bb_start);
time(index_bb_start+1) = time(index_bb_start) + 0.001;
x_ball(index_bb_start+1) = x_ball(index_bb_start);
y_ball(index_bb_start+1) = y_ball(index_bb_start);

l = index_bb_start+2;
r=0
t_step = 0.001;
tol_angle = 36;
while theta_bb(l-1) > -tol_angle*pi/180
    time(l) = time(l-1) + t_step;
    omega_bb(l) = omega_bb(l-1) + alpha_bb * t_step;
    theta_bb(l) = theta_bb(l-1) + omega_bb(l-1)*t_step + 1/2*alpha_bb*t_step^2;
    [r_x, r_y] = get_r(x_ball(l-1), y_ball(l-1));
    r = sqrt(r_x^2 + r_y^2);
    a_x = alpha_bb * r_y;
    a_y = alpha_bb * r_x;
    vx_ball(l) = -omega_bb(l-1)*r_y;
    x_ball(l) = x_ball(index_bb_start+1) - theta_bb(l-1)*r_y;
    vy_ball(l) = omega_bb(l-1)*r_x;
    y_ball(l) = y_ball(index_bb_start+1) + theta_bb(l-1)*r_x;
    omega(l) = 0;
    data = [time(l), x_ball(l), y_ball(l), vx_ball(l), vy_ball(l), omega(l), theta_bb(l), omega_bb(l)]
    l = l+1;
end

%% Ramp BB
%Update data for ball values
index_ramp_bb = length(time);
theta_bb(index_ramp_bb+1) = theta_bb(index_ramp_bb)
vx_ball(index_ramp_bb+1) = 0;
vy_ball(index_ramp_bb+1) = vy_ball(index_ramp_bb);
omega(index_ramp_bb+1) = 0;
time(index_ramp_bb+1) = time(index_ramp_bb) + 0.001;
x_ball(index_ramp_bb+1) = x_ball(index_ramp_bb);
y_ball(index_ramp_bb+1) = y_ball(index_ramp_bb);
l = index_ramp_bb+2;
a_x = 0.8*g*sind(tol_angle)*cosd(tol_angle); %0.8 is a drag factor
a_y = -0.8*g*sind(tol_angle)*sind(tol_angle);
t_step = 0.001;
while r < 0.0275
    [r_x, r_y] = get_r(x_ball(l-1), y_ball(l-1));
    r = sqrt(r_x^2 + r_y^2);
    time(l) = time(l-1) + t_step;
    vx_ball(l) = vx_ball(l-1) + t_step*a_x;
    x_ball(l) = x_ball(l-1) + t_step*vx_ball(l-1) + 1/2*t_step^2*a_x;
    vy_ball(l) = vy_ball(l-1) + t_step*a_y;
    y_ball(l) = y_ball(l-1) + t_step*vy_ball(l-1) + 1/2*t_step^2*a_y;
    omega(l) = -(sqrt(vx_ball(l)^2 + vy_ball(l)^2)/r_ball); % Non-slip condition
    theta_bb(l) = theta_bb(index_ramp_bb);
    l = l+1;
end
data = [time, x_ball, y_ball, vx_ball, vy_ball, omega];

%% Another another failed attempt at flight
% Assume projectile motion
index_proj_bb = length(time);
theta_bb(index_proj_bb) = 0;
vx_ball(index_proj_bb+1) = vx_ball(index_proj_bb);
vy_ball(index_proj_bb+1) = vy_ball(index_proj_bb);


v = sqrt((vx_ball(index_proj_bb+1))^2 + (vy_ball(index_proj_bb+1))^2);

omega(index_proj_bb+1) = r_ball*v;
time(index_proj_bb+1) = time(index_proj_bb) + 0.001;


x_impulse (index_proj_bb+1,1) = m_ball*(vx_ball(index_proj_bb+1)-vx_ball(index_proj_bb));
y_impulse(index_proj_bb+1,1)= m_ball*(vy_ball(index_proj_bb+1)-vy_ball(index_proj_bb));
angular_impulse(index_proj_bb+1,1) = I_ball*(omega(index_proj_bb+1)-omega(index_proj_bb));
time(index_proj_bb+1) = time(index_proj_bb) + 0.001;



x_ball(index_proj_bb+1) = x_ball(index_proj_bb);
y_ball(index_proj_bb+1) = y_ball(index_proj_bb);
l = index_proj_bb+2;
a_x = 0;
a_y = -g;
t_step = 0.001;
% syms t2
% t_max_projectile_2 = vpasolve(-1/2*g*t2^2-255/441*vx_ball(index_proj_bb)*t2+0.0255,t2)

while y_ball(l-1) > 0.088 %
    time(l) = time(l-1) + t_step;
    vx_ball(l) = vx_ball(l-1) + t_step*a_x;
    x_ball(l) = x_ball(l-1) + t_step*vx_ball(l-1) + 1/2*t_step^2*a_x;
    vy_ball(l) = vy_ball(l-1) + t_step*a_y;
    y_ball(l) = y_ball(l-1) + t_step*vy_ball(l-1) + 1/2*t_step^2*a_y;
    omega(l) = omega(index_proj_bb+1); % Non-slip condition
    l = l+1;
end


angle_ball(k) = 0;
for k = 1:length(omega)-1
    angle_ball(k+1) = angle_ball(k)+(omega(k+1) +omega(k))*(time(k+1)-time(k))/2;
end

%animate_ball(theta_bb, 0, time, vx_ball,vy_ball, r_ball, x_ball, y_ball, angle_ball,0,0)
%comet1(1000*x_ball,1000*y_ball,0)

%% Ball takes break from bounty hunting and visits local amusement park
%Introducing 'Little Guy'
omega_lg = zeros(length(time)+1,1);
theta_lg = zeros(length(time)+1,1);
alpha_lg = 0.1;

index_lg_start = length(time);
omega_lg(index_lg_start+1) = -m_ball*0.0117*vy_ball(index_lg_start)/15387.62*10^9.5
vx_ball(index_lg_start+1) = 0; %omega(index_lg_start)*r_ball;
vy_ball(index_lg_start+1) = 0;

v = sqrt((vx_ball(index_lg_start+1))^2 + (vy_ball(index_lg_start+1))^2);

omega(index_lg_start+1) = v*r_ball;
time(index_lg_start+1) = time(index_lg_start) + 0.001;


x_impulse (index_lg_start+1,1) = m_ball*(vx_ball(index_lg_start+1)-vx_ball(index_lg_start));
y_impulse(index_lg_start+1,1)= m_ball*(vy_ball(index_lg_start+1)-vy_ball(index_lg_start));
angular_impulse(index_lg_start+1,1) = I_ball*(omega(index_lg_start+1)-omega(index_lg_start));
time(index_lg_start+1) = time(index_lg_start) + 0.001;



x_ball(index_lg_start+1) = x_ball(index_lg_start);
y_ball(index_lg_start+1) = y_ball(index_lg_start);
l = index_lg_start+2;
r=0
t_step = 0.001;
tol_angle = 18.5;
while theta_lg(l-1) < tol_angle*pi/180
    time(l) = time(l-1) + t_step;
    omega_lg(l) = omega_lg(l-1);
    theta_lg(l) = theta_lg(l-1) + omega_lg(l)*t_step;
    [r_x, r_y] = get_r2(x_ball(l-1), y_ball(l-1));
    vx_ball(l) = omega_lg(l-1)*r_y;
    x_ball(l) = x_ball(index_lg_start+1) - theta_lg(l-1)*r_y;
    vy_ball(l) = omega_lg(l-1)*r_x;
    y_ball(l) = y_ball(index_lg_start+1) + theta_lg(l-1)*r_x;
    omega(l) = 0;
    data = [time(l), x_ball(l), y_ball(l), vx_ball(l), vy_ball(l), omega(l), theta_lg(l), omega_lg(l)]
    l = l+1;
end
%% Ramp into the train
% Assume Straight slope
index_slope_lg_start = length(time);
theta_lg(index_slope_lg_start+1) = theta_lg(index_slope_lg_start)
vx_ball(index_slope_lg_start+1) = 0; %omega(index_slope_lg_start)*r_ball;
vy_ball(index_slope_lg_start+1) = 0;%vy_ball(index_slope_lg_start);
omega(index_slope_lg_start+1) = 0; %omega(index_slope_lg_start);
time(index_slope_lg_start+1) = time(index_slope_lg_start) + 0.001;
x_ball(index_slope_lg_start+1) = x_ball(index_slope_lg_start);
y_ball(index_slope_lg_start+1) = y_ball(index_slope_lg_start);
l = index_slope_lg_start+2;
a_x = -g*sind(tol_angle)*cosd(tol_angle); %0.8 is a drag factor
a_y = -g*sind(tol_angle)*sind(tol_angle);
t_step = 0.001;

while y_ball(l-1) - y_ball(index_slope_lg_start+1) > -0.0111
    time(l) = time(l-1) + t_step;
    vx_ball(l) = vx_ball(l-1) + t_step*a_x;
    x_ball(l) = x_ball(l-1) + t_step*vx_ball(l) + 1/2*t_step^2*a_x;
    vy_ball(l) = vy_ball(l-1) + t_step*a_y;
    y_ball(l) = y_ball(l-1) + t_step*vy_ball(l) + 1/2*t_step^2*a_y;
    omega(l) = sqrt(vx_ball(l)^2 + vy_ball(l)^2)/r_ball; % Non-slip condition
    theta_lg(l) = theta_lg(index_slope_lg_start);
    l = l+1;
end

%% Ball embarks on averageness
index_projectile_lg_start = length(time);
vx_ball(index_projectile_lg_start+1) = vx_ball(index_projectile_lg_start);
vy_ball(index_projectile_lg_start+1) = vy_ball(index_projectile_lg_start);
omega(index_projectile_lg_start+1) = omega(index_projectile_lg_start);
time(index_projectile_lg_start+1) = time(index_projectile_lg_start) + 0.001;
x_ball(index_projectile_lg_start+1) = x_ball(index_projectile_lg_start);% + 0.001*vx_ball(index_projectile_lg_start);
y_ball(index_projectile_lg_start+1) = y_ball(index_projectile_lg_start);
l = index_projectile_lg_start+2;
%Zero acceleration in the x direction for projectile motion
a_x = 0;
a_y = -g;
t_step = 0.001;
% syms t
% t_max_projectile = vpasolve(-1/2*g*t^2+(vx_ball(index_projectile_lg_start)*28/53+vy_ball(index_projectile_lg_start))*t + 0.0680377,t)
while y_ball(l-1) > 0.058
time(l) = time(l-1) + t_step;
vx_ball(l) = vx_ball(l-1) + t_step*a_x;
x_ball(l) = x_ball(l-1) + t_step*vx_ball(l-1) + 1/2*t_step^2*a_x;
vy_ball(l) = vy_ball(l-1) + t_step*a_y;
y_ball(l) = y_ball(l-1) + t_step*vy_ball(l-1) + 1/2*t_step^2*a_y;
omega(l) = omega(index_projectile_lg_start); % Non-slip condition
l = l+1;
end

%% Ball smashes cart
% Perfectly inelastic collision ie. Stuck together for the duration of the
% rotation
index_incart_start = length(time);
vx_ball(index_incart_start+1) = 0;
vy_ball(index_incart_start+1) = 0;
omega(index_incart_start+1) = 0;
time(index_incart_start+1) = time(index_incart_start) + 0.001;


x_impulse (index_incart_start+1,1) = m_ball*(vx_ball(index_incart_start+1)-vx_ball(index_incart_start));
y_impulse(index_incart_start+1,1)= m_ball*(vy_ball(index_incart_start+1)-vy_ball(index_incart_start));
angular_impulse(index_incart_start+1,1) = I_ball *omega(index_incart_start);
time(index_incart_start+1) = time(index_incart_start) + 0.001;


x_ball(index_incart_start+1) = x_ball(index_incart_start);
y_ball(index_incart_start+1) = y_ball(index_incart_start);
l = index_incart_start+2;
angle_cart = 9.83;
a_x = g*sind(angle_cart)*cosd(angle_cart); %0.8 is a drag factor
a_y = -g*sind(angle_cart)*sind(angle_cart);
t_step = 0.001;
while x_ball(l-1) < 0.03
    time(l) = time(l-1) + t_step;
    vx_ball(l) = vx_ball(l-1) + t_step*a_x;
    x_ball(l) = x_ball(l-1) + t_step*vx_ball(l-1) + 1/2*t_step^2*a_x;
    vy_ball(l) = vy_ball(l-1) + t_step*a_y;
    y_ball(l) = y_ball(l-1) + t_step*vy_ball(l-1) + 1/2*t_step^2*a_y;
    omega(l) = -sqrt(vx_ball(l)^2 + vy_ball(l)^2)/r_ball; % Non-slip condition
    l = l+1;
end

%% Cart
for i = 1:length(time)+1
    x_cart(i,1) = 44;
    y_cart(i,1) = 57;
end

%% Ball rides cart
% Perfectly inelastic collision ie. Stuck together for the duration of the
% rotation
index_withcart_start = length(time);
vx_ball(index_withcart_start+1) = 0;
vy_ball(index_withcart_start+1) = 0;
omega(index_withcart_start+1) = 0;
time(index_withcart_start+1) = time(index_withcart_start) + 0.001;


x_ball(index_withcart_start+1) = x_ball(index_withcart_start);
y_ball(index_withcart_start+1) = y_ball(index_withcart_start);
l = index_withcart_start+2;
angle_cart = 9.83;
a_x = g*sind(angle_cart)*cosd(angle_cart); %0.8 is a drag factor
a_y = -g*sind(angle_cart)*sind(angle_cart);
t_step = 0.001;
while x_ball(l-1) < 0.18
    time(l) = time(l-1) + t_step;
    vx_ball(l) = vx_ball(l-1) + t_step*a_x;
    x_ball(l) = x_ball(l-1) + t_step*vx_ball(l-1) + 1/2*t_step^2*a_x;
    x_cart(l) = x_cart(l-1) +  1000*(x_ball(l) -x_ball(l-1));
    vy_ball(l) = vy_ball(l-1) + t_step*a_y;
    y_ball(l) = y_ball(l-1) + t_step*vy_ball(l-1) + 1/2*t_step^2*a_y;
    y_cart(l) = y_cart(l-1) + 1000*(y_ball(l)- y_ball(l-1)) ;
    omega(l) = 0; % Non-slip condition
    l = l+1;
    index_rotate = l;
end

%% Ball rides cart (flat)
% Perfectly inelastic collision ie. Stuck together for the duration of the
% rotation
index_withcart_start = length(time);
vx_ball(index_withcart_start+1) = vx_ball(index_withcart_start);
vy_ball(index_withcart_start+1) = 0;
omega(index_withcart_start+1) = 0;
time(index_withcart_start+1) = time(index_withcart_start) + 0.001;
x_ball(index_withcart_start+1) = x_ball(index_withcart_start);
y_ball(index_withcart_start+1) = y_ball(index_withcart_start);
x_cart(index_withcart_start+1) = x_cart(index_withcart_start);
y_cart(index_withcart_start+1) = y_cart(index_withcart_start);
l = index_withcart_start+2;
t_step = 0.001;
while x_ball(l-1) < 0.28
    time(l) = time(l-1) + t_step;
    vx_ball(l) = vx_ball(l-1);
    x_ball(l) = x_ball(l-1) + t_step*vx_ball(l-1);
    x_cart(l) = x_cart(l-1) + 1000*(x_ball(l) - x_ball(l-1));
    vy_ball(l) = vy_ball(l-1);
    y_ball(l) = y_ball(l-1) + t_step*vy_ball(l-1);
    y_cart(l) = y_cart(l-1) + 1000*(y_ball(l) - y_ball(l-1));
    omega(l) = 0; % Non-slip condition
    l = l+1;
end
data = [time, x_ball, y_ball, vx_ball, vy_ball, omega];
    

x_impulse (l-1) = m_ball*vx_ball(l-1);
y_impulse(l-1)= m_ball*(vy_ball(l-1));
angular_impulse(l-1) = 0;
time(l-1) = time(l-1);


%% Cart Rotation 

%start rotating at index 2350

for i = 1:length(time)+1
    theta_cart(i,1) = -9.81;
    if (i > index_rotate) && (theta_cart(i,1) < 0)
        
        theta_cart(i,1) = theta_cart(i-1,1) + .25;
    end
    
    if (theta_cart(i,1) > 0)
        
        theta_cart(i,1) = 0;
        
    end
    
        
end

%--------------------------------------------------------


%% rotator big return to original position

%start rotating at index 2350
% 
theta_bb(index_proj_bb:length(time)) = 0;
for i = 1:length(time)
%     if (i < 1646)
%         
%     theta_bb(i,1) = 0;
%     end
%     
    if (i > 1870)
        
        theta_bb(i,1) = theta_bb(i-1,1) + 0.001;
        
    end
    
    if (i > 2919) 
        
       theta_bb(i,1) = 0;
        
    end
    
        
end

%--------------------------------------------------------


%% Rotating Figures
data = [time, x_ball, y_ball, vx_ball, vy_ball, omega];
theta_bb(index_proj_bb:length(time)) = 0;
theta_lg(index_projectile_lg_start:length(time)) = 0;
%% Ball Angle
angle_ball(k) = 0;
for k = 1:length(omega)-1
    angle_ball(k+1) = angle_ball(k)+(omega(k+1) +omega(k))*(time(k+1)-time(k))/2;
end
%% Draw Board and Animate
figure
window_x = [0;0;x_right];
window_y = [y_top;0;0];
draw_board1();
start_val  = 1;
animate_ball(theta_bb(start_val:length(theta_bb)), theta_lg(start_val:length(x_ball)), time(start_val:length(x_ball)), vx_ball(start_val:length(x_ball)),vy_ball(start_val:length(x_ball)), r_ball, x_ball(start_val:length(x_ball)), y_ball(start_val:length(x_ball)), angle_ball(start_val:length(x_ball)),x_cart(start_val:length(x_ball)),y_cart(start_val:length(x_ball)), theta_cart(start_val:length(x_ball)))
%comet(1000*x_ball(start_val:length(x_ball)),1000*y_ball(start_val:length(y_ball)),0)
data = [time, x_ball, y_ball, vx_ball, vy_ball, omega]

%% Draw Plots
[impulse, impulse_angle]= impulse_linear(x_impulse, y_impulse, time);
[a_linear, a_linear_angle] = acceleration_linear(vx_ball, vy_ball, time);
a_angular = acceleration_anglular(omega, time);
plot_graphs(x_ball,y_ball,vx_ball, vy_ball, omega, a_linear, a_linear_angle, a_angular, impulse,impulse_angle, angular_impulse, time)

