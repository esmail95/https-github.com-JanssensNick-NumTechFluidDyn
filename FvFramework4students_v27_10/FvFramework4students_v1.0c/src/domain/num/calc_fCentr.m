function fCentr = calc_fCentr(nF,fNbVLoc,fNbV,vCoord)
    fCentr=[];
    a=nF;
    b=fNbVLoc;
    c=fNbV;
    d=vCoord;
    
for i=1:2:2*nF   
    k=[fNbV(i),fNbV(i+1)];
    fCentr=[fCentr,0.5*(vCoord(:,k(1))+vCoord(:,k(2)))];    
end

end



