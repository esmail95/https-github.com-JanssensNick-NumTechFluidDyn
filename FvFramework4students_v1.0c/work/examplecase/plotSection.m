function plotSection(cCoord,result,indices,xy,LineWidth)
%PLOTSECTION Summary of this function goes here
%   Detailed explanation goes here
    x = 0*indices;
    y = 0*indices;
    cntr = 1;
    if xy == 'x'
        for k = indices
            x(cntr) = cCoord(1,k);
            y(cntr) = result(k);
            cntr = cntr + 1;
        end
    else
        for k = indices
            x(cntr) = cCoord(2,k);
            y(cntr) = result(k);
            cntr = cntr + 1;
        end
    end 
    %scatter(x,y)
    plot(x,y,'-','LineWidth',LineWidth)
    
end

