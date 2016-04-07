function [t,x] = solve_spring(m_ball)
tspan = 0:0.001:0.5;
x0 = [-0.024;0];
[t,x] = ode45(@spring,tspan,x0);

function dxdt = spring(t,x)
k=400; %Experimentally measured
m= m_ball + 0.15; %Mass of Ball + Mass of Mass of rod + spring
%Ffr = -0.5;
dxdt = [x(2); -k/m*x(1)]; %/m
end
end