function [im_mag, angle]= impulse_linear(imx, imy, time)
for i = 1:length(time)

    im_mag(i) = sqrt(imx(i)^2 + imy(i)^2);
    angle(i) = atan2(imx(i),imy(i));
end
im_mag(1) = 0;
im_mag(length(time)) = 0;
angle(1) = 0;
angle(length(time)) = 0;

end