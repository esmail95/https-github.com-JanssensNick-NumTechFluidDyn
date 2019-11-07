function fArea = calc_fArea(nF,fNbVLoc,fNbV,vCoord)
   
  fArea=[];
    
for i=1:2:2*nF   
    k=[fNbV(i),fNbV(i+1)];
    fArea_vec=vCoord(:,k(1))-vCoord(:,k(2)); 
    fArea=[fArea,norm(fArea_vec)];
end
  

end



