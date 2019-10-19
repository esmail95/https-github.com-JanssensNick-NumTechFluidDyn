function fNormal = calc_fNormal(fNbVLoc,fNbV,vCoord,fArea)
   nN = length(fArea);
   fNormal = zeros(fNbVLoc,nN);
   for i = 1:nN
        v1 = fNbV(2*i-1);
        v2 = fNbV(2*i);
        %vh = max(v1,v2);
        %vl = min(v1,v2);
        tx = vCoord(1,v2) - vCoord(1,v1);
        ty = vCoord(2,v2) - vCoord(2,v1);
        fNormal(1,i) = -ty/fArea(i);
        fNormal(2,i) = tx/fArea(i);
   end
end