function draw_board %(x_t)
clc
close all
clear all

hold all

%draw frame
line([0, 300],[0, 0], 'LineWidth', 1, 'Color', 'black');
line([-15.4, 300],[300, 300], 'LineWidth', 1, 'Color', 'black');
line([300, 300],[0, 300], 'LineWidth', 1, 'Color', 'black');
line([0, 0],[0, 300-8.5-20.4-8.5], 'LineWidth', 1, 'Color', 'black');
%draw L shape 
line([-15.4, -15.4],[300-136.7, 300], 'LineWidth', 1, 'Color', 'black');
line([0, 300],[300-8.5, 300-8.5], 'LineWidth', 1, 'Color', 'black');
line([0, 0],[300-8.5-20.4, 300-8.5], 'LineWidth', 1, 'Color', 'black');
line([0, 170],[300-8.5-20.4-8.5, 300-8.5-20.4-8.5], 'LineWidth', 1, 'Color', 'black');
line([0, 170],[300-8.5-20.4, 300-8.5-20.4], 'LineWidth', 1, 'Color', 'black');
line([170, 170],[300-8.5-20.4-8.5, 300-8.5-20.4], 'LineWidth', 1, 'Color', 'black');
line([-15.4, 0],[300-136.7, 300-136.7], 'LineWidth', 1, 'Color', 'black');
%draw track with ramp
line([80, 170],[178.5, 178.5], 'LineWidth', 1, 'Color', 'black');
line([80, 170],[170, 170], 'LineWidth', 1, 'Color', 'black');
line([80, 80],[170, 170], 'LineWidth', 1, 'Color', 'black');
line([170, 170],[178.5, 178.5], 'LineWidth', 1, 'Color', 'black');

line([170, 150.5, 97, 92.06866, 125.54727, 80, 80, 170, 170], [178.5, 178.5, 197, 190.07672, 178.5, 178.5, 170, 170, 176.7], 'LineWidth', 1, 'Color', 'black');

%draw pipe
line([79.1, 274.1],[114, 114], 'LineWidth', 1, 'Color', 'black')
line([79.1, 274.1],[114 + 25.5, 114 + 25.5], 'LineWidth', 1, 'Color', 'black')
line([79.1, 79.1],[114, 114 + 25.5], 'LineWidth', 1, 'Color', 'black')
line([274.1, 274.1],[114, 114 + 25.5], 'LineWidth', 1, 'Color', 'black')
line([274.1+14.5,300+10],[100+0.965*14.5,100+25.5+0.965*10],'LineWidth', 1, 'Color', 'black')
line([274.1,274.1+14.5],[114,100+0.965*14.5],'LineWidth', 1, 'Color', 'black')

%draw track under pipe
line([44.3, 199.3],[114, 114], 'LineWidth', 1, 'Color', 'black');
line([44.3, 199.3],[114-8.5, 114-8.5], 'LineWidth', 1, 'Color', 'black');
line([44.3, 44.3],[114-8.5, 114], 'LineWidth', 1, 'Color', 'black');
line([199.3, 199.3],[114-8.5, 114], 'LineWidth', 1, 'Color', 'black');

%bottom portion of track right corner
line([300-8.5, 300-8.5],[8.5, 41.5], 'LineWidth', 1, 'Color', 'black');
line([300-8.5, 300],[41.5, 41.5], 'LineWidth', 1, 'Color', 'black');
line([300-8.5-111.5, 300-8.5],[8.5, 8.5], 'LineWidth', 1, 'Color', 'black');
%plot angled portion of track
line([180, -21.9901, -17.32489, -32.30172, -38.39684, 180],[8.5, 43.5, 70.42361, 73.01873, 37.84289, 0],'LineWidth', 1, 'Color', 'black');


% ang = 9.83*pi/180;
% line([300-8.5-111.5, 300-8.5-111.5 - 221.65*cos(ang)],[0, 221.65*sin(ang)], 'LineWidth', 1, 'Color', 'black');
% line([300-8.5-111.5, 300-8.5-111.5 - (221.65*cos(ang)) + (15.2*cos(ang)+8.5*sin(ang))],[8.5/cos(ang), 221.65*sin(ang) - 15.2*sin(ang) + 8.5*cos(ang)], 'LineWidth', 1, 'Color', 'black');
% 
% line([300-8.5-111.5 - 221.65*cos(ang),300-8.5-111.5 - 221.65*cos(ang)+35.7*sin(ang) ],[221.65*sin(ang), 221.65*sin(ang)+35.7*cos(ang)], 'LineWidth', 1, 'Color', 'black');
% line([300-8.5-111.5 - (221.65*cos(ang)) + (15.2*cos(ang)+8.5*sin(ang)),300-8.5-111.5 - (221.65*cos(ang)) + (15.2*cos(ang)+8.5*sin(ang))+(35.7-8.5)*sin(ang) ],[221.65*sin(ang) - 15.2*sin(ang) + 8.5*cos(ang), 221.65*sin(ang) - 15.2*sin(ang) + 8.5*cos(ang) + (35.7-8.5)*cos(ang)], 'LineWidth', 1, 'Color', 'black');
% 
% 
% line([-38.395877 + (35.7*sin(ang)), -38.395877 + (35.7*sin(ang))+15.2*cos(ang)],[221.65*sin(ang) + 35.7*cos(ang),(221.65-15.2)*sin(ang) + 35.7*cos(ang) ], 'LineWidth', 1, 'Color', 'black');
%line([(xangle(1)-14.97684336), (xangle(length(xangle))-14.97684336)],[(-tan(9.83*pi/360)*(xangle(1)-14.97684336)) + (tan(9.83*pi/360)*(300-111.5-8.5))+8.5,(-tan(9.83*pi/360)*xangle(length(xangle)) + (tan(9.83*pi/360)*(300-111.5-8.5)) + 8.5)], 'black')

%first rotator
% line([13.8, 37.8],[114-9.5, 114-9.5], 'LineWidth', 1, 'Color', 'black');
% line([13.8, 37.8],[114-9.5+14.9, 114-9.5+14.9], 'LineWidth', 1, 'Color', 'black');
% line([0, 13.8],[114-9.5, 114-9.5], 'LineWidth', 1, 'Color', 'black');

% xc = [170:0.1:(170+55.7)];
% for i = 1:length(xc)
%     yc(i) = sqrt(55.7^2 - (xc(i)-170)^2) + 234.2;
%     y_c(i) = -sqrt(55.7^2 - (xc(i)-170)^2) + 234.2;
%     yout(i) = sqrt(57.5^2 - (xc(i)-170)^2) + 234.2;
%     y_out(i) = -sqrt(57.5^2 - (xc(i)-170)^2) + 234.2;
% end
% plot(xc, yc, 'black')
% plot(xc, y_c, 'black')
% plot(xc, yout, 'black')
% plot(xc, y_out, 'black')
 rc = 55.7;
 c_angle = [(pi/2):-0.01:(-pi/2)];
 for i = 1:length(c_angle)
     xc(i) = rc*cos(c_angle(i)) + 170;
     yc(i) = rc*sin(c_angle(i)) + 234.2;
 end
 
 rc = 57.5;
 for i = 1:length(c_angle)
     xc1(i) = rc*cos(c_angle(i)) + 170;
     yc1(i) = rc*sin(c_angle(i)) + 234.2;
 end
 
 plot(xc, yc, 'LineWidth', 1, 'Color', 'black');
 plot(xc1, yc1, 'LineWidth', 1, 'Color', 'black');
 
 line([170,170],[289.9, 291.7], 'LineWidth', 1, 'Color', 'black');
%plot second ramp 
%line([0, 79.1],[88.5*2, 139.5], 'LineWidth', 1, 'Color', 'black');

line([77.80304226, 0.26416+(13*sin(0.7))],[138.1990043, 173.6663501], 'LineWidth', 1, 'Color', 'black');

axis equal
axis([-60 330 -30 330])
end