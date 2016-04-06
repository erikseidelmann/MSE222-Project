function [t,theta,omega_curve] = solve_curve(r_ball,r_curve,omega)
tspan = 0:0.001:0.14;
x0 = 0;%omega*r_ball/(r_curve-r_ball)];
[t,theta] = ode113(@curve,tspan,x0);
for j = 1:length(t)
omega_curve(j) = curve(1,theta(j));
end

function dthetadt = curve(t,theta)
dthetadt = sqrt(((7/5)*r_ball^2*omega^2 + 2*9.81*(r_curve-r_ball) - 9.81*(1+cos(theta))*(r_curve-r_ball))/((7/5)*(r_curve-r_ball)^2)); %/m
end
end