function l = lambda(cCoordnb1,cCoordnb2,fCoord)
%LAMBDA Summary of this function goes here
%   Detailed explanation goes here
    x1 = cCoordnb1(1);
    x2 = cCoordnb2(1);
    y1 = cCoordnb1(2);
    y2 = cCoordnb2(2);
    xf = fCoord(1);
    yf = fCoord(2);
    d = sqrt((x2-x1)^2+(y2-y1)^2);
    d1 = sqrt((xf-x1)^2+(yf-y1)^2);
    l = d1/d;
end

