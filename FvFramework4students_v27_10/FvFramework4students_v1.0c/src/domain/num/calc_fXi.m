function fXi = calc_fXi(nF,cCoord,fNbCLoc,fNbC)
g=cCoord;
h=fNbCLoc;
psi=fNbC;
    fXi=[];
    a=size(cCoord);
for i=1:2:2*nF   
    k=[fNbC(i),fNbC(i+1)];
    fXi=[fXi,cCoord(:,k(2))-cCoord(:,k(1))];    
end

end




