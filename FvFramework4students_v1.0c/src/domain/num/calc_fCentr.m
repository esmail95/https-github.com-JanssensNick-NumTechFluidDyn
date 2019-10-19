function fCentr = calc_fCentr(nF,fNbVLoc,fNbV,vCoord)
    %loop over faces, for each face find the neighbor vertices
    fCentr = zeros(fNbVLoc,nF);
    for i = 1:nF
        v1 = fNbV(2*i-1);
        v2 = fNbV(2*i);
        xCentr = (vCoord(1,v1) + vCoord(1,v2))/2;
        yCentr = (vCoord(2,v1) + vCoord(2,v2))/2;
        fCentr(1,i) = xCentr;
        fCentr(2,i) = yCentr;
    end

end



