function fNormal = calc_fNormal(fNbVLoc,fNbV,vCoord,fArea)
   nF = length(fArea);
   fNormal = zeros(fNbVLoc,nF);
   for i = 1:nF
        v_1 = fNbV(2*i-1);
        v_2 = fNbV(2*i);
        x_dist = vCoord(1,v_2) - vCoord(1,v_1);
        y_dist = vCoord(2,v_2) - vCoord(2,v_1);
        fNormal(1,i) = -y_dist/fArea(i);
        fNormal(2,i) = x_dist/fArea(i);
   end
    
end

