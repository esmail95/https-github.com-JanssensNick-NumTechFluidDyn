function fArea = calc_fArea(nF,fNbVLoc,fNbV,vCoord)
    fArea = zeros(1,nF);
    for i = 1:nF
        v1 = fNbV(2*i-1);
        v2 = fNbV(2*i);
        fArea(1,i) = sqrt((vCoord(1,v1)-vCoord(1,v2))^2+(vCoord(2,v1)-vCoord(2,v2))^2);
    end 
end



