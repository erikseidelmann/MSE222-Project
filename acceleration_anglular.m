function [a_ang]= acceleration_anglular(omega, time)
for i = 2:length(time)-1
    a_ang(i) = (omega(i+1) - omega(i-1))/(time(i+1) - time(i-1));
end
a_ang(1) = 0;
a_ang(length(time)) = 0;
end