function fXi = calc_fXi(nF,cCoord,fNbCLoc,fNbC)
   fXi = zeros(fNbCLoc,nF);
   for i = 1:nF
        c1 = fNbC(2*i-1);
        c2 = fNbC(2*i);
        ch = max(c1,c2);
        cl = min(c1,c2);
        fXi(1,i) = cCoord(1,ch) - cCoord(1,cl);
        fXi(2,i) = cCoord(2,ch) - cCoord(2,cl);
   end
end



