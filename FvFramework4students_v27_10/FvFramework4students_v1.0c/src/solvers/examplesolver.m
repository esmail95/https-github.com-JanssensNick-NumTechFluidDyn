%==========================================================================
%
% Example solver using the FVMLab framework 4 students
%
% Purpose: Provides code structure for solving a scalar conservation
%          equation using data structures provided by the framework.
%
% by Frederik Rogiers
%
%==========================================================================
function result = examplesolver(casedef)

dom = casedef.dom;
material=casedef.material;
% Create field objects
T = Field(dom.allCells,0);      % Temperature [K] (scalar); empty field
reset(T,0);                     % Reset with all zeros

% ... Create all other required data structures
nC=dom.nC;
nPc=dom.nPc;
nGc=dom.nGc;
fNbC=dom.fNbC;
fNbCLoc=dom.fNbCLoc;
nIf=dom.nIf;
nBf=dom.nBf;
nF=dom.nF;
cNbF=dom.cNbF;
vCoord=dom.vCoord;
cCoord=dom.cCoord;
fNbV=dom.fNbV;
fNbVLoc=dom.fNbVLoc;

U=casedef.vars.U.data; %put to zero to have just diffusion
k=material.k;
fArea = calc_fArea(nF,dom.fNbVLoc,dom.fNbV,dom.vCoord);
fXi = calc_fXi(nF,dom.cCoord,dom.fNbCLoc,fNbC);
fNormal = calc_fNormal(fNbVLoc,fNbV,vCoord,fArea);
fTangent = calc_fTangent(nF,fNormal);
fXiLambda = calc_fXiLambda(...
      nF,...
      fXi,fTangent,cCoord,vCoord,...
      fNbCLoc,fNbC,fNbVLoc,fNbV);
adata=[];
adiag=zeros(nC,1);
aNb_int=zeros(2*nIf,1);
aNb_bound=zeros(2*nBf,1);
b=zeros(nC,1);
BC = casedef.BC;
nbZones = length(BC);
ranges = zeros(2,nbZones);
% Create an equation object for holding a scalar conservation equation
eqn = ScalarFvEqn2(dom);

iterate = true;
niter = 0;
while iterate   
   niter = niter+1;
   
   % Set all terms to zero
reset(eqn); 
adiag=zeros(nC,1);
aNb_int=zeros(2*nIf,1);
aNb_bound=zeros(2*nBf,1);
b=zeros(nC,1);
   
for i = 1:nbZones
    zone = BC{i}.zoneID;
    boundary = dom.getzone(zone);
    ranges(:,i) = boundary.range;
end

% Create an equation object for holding a scalar conservation equation
eqn = ScalarFvEqn2(dom);

iterate = true;
niter = 0;
while iterate   
   niter;
   niter = niter+1;
   
   % Set all terms to zero
   reset(eqn); 
   
   % Compute coefficients for physical cell eqns and add them to eqn object
   for i = 1:nIf
      nb1 = fNbC(2*i-1);
      nb2 = fNbC(2*i);
      Af = fArea(i);
      Xi = norm(fXi(:,i));
      Normal=fNormal(:,i);
      %fXiLambda weighted average between main cell and neighbour to find phi_f  
      XiLambda=fXiLambda(i);
      %advective term
      U_f=XiLambda*U(:,nb1)+(1-XiLambda)*U(:,nb2);
      U_n=dot(U_f,Normal); 
      
      % Diagonal elements
      adiag(nb1) = adiag(nb1) - k*Af/Xi - Af*U_n*XiLambda;
      adiag(nb2) = adiag(nb2) - k*Af/Xi - Af*U_n*XiLambda;
      % Off-diagonal elements for internal faces
      aNb_int(2*i) = aNb_int(2*i) + k*Af/Xi + Af*U_n*XiLambda; % Benedendiagonaal
      aNb_int(2*i-1) = aNb_int(2*i-1) + k*Af/Xi + Af*U_n*XiLambda; % Bovendiagonaal
   end
   % Compute coefficients for ghost cell eqns and add them to eqn object
   for i = 1:nBf
      nb1 = fNbC(2*i-1 + 2*nIf); % Physicl cell
      nb2 = fNbC(2*i + 2*nIf); % Ghost Cell
      Af = fArea(i + nIf);
      Xi = norm(fXi(:,i + nIf));
      Normal=fNormal(:,i+nIf);
      %fXiLambda weighted average between main cell and neighbour to find phi_f  
      XiLambda=fXiLambda(i+nIf);
      %advective term
      U_f=XiLambda*U(:,nb1)+(1-XiLambda)*U(:,nb2);
      U_n=(dot(U_f,Normal));
      
      boundaryFound = 0;
      for j = 1:nbZones
        % Determine to which boundary the face belongs
        if ismember(i + nIf,ranges(1,j):ranges(2,j)) && boundaryFound == 0
            
            switch BC{j}.kind
                case 'Dirichlet' 
                    %l*phi(GC) + (1-l)*phi(PC) = phi*
                    l = 0.5;
                    % Diagonal elements from BC: ghost cell
                    adiag(nb2) = l;
                    % Diagonal elements: physical cells
                    adiag(nb1) = adiag(nb1) - k*Af/Xi - Af*U_n*XiLambda;
                    % Off-diagonal elements: physical cell
                    aNb_bound(2*i) = 1-l;
                    % Off-diagonal elements: ghost cells
                    aNb_bound(2*i-1) = aNb_bound(2*i-1) + k*Af/Xi + Af*U_n*XiLambda;
                    
                    % Forcing terms
                    phi_star = BC{j}.data.bcval;
                    %b(i+nPc) = phi_star;  
                    b(nb2) = phi_star;
                case 'Neumann'
                    % (phi(PC)-phi(GC))/Ef=phi*
                    % Diagonal elements from BC: ghost cell
                    adiag(nb2) = -1/Xi;
                    % Diagonal elements: physical cells
                    adiag(nb1) = adiag(nb1) - k*Af/Xi - Af*U_n*XiLambda;
                    
                    % Off-diagonal elements: physical cell
                    aNb_bound(2*i) = 1/Xi;
                    % Off-diagonal elements: ghost cells
                    aNb_bound(2*i-1) = aNb_bound(2*i-1) + k*Af/Xi + Af*U_n*XiLambda;
                    
                    % Forcing terms
                    phi_star = BC{j}.data.bcval;
                    b(i+nPc) = phi_star;  
            end
        end
      end
   end
   % eqn.adata = rand(eqn.nnz,1); % just a meaningless example adata
   % eqn.bdata = rand(eqn.n,1);   % just a meaningless example bdata
   eqn.adata = [adiag ; aNb_int ; aNb_bound];
   eqn.bdata = b;
     
   % Create a matlab sparse linear system from the eqn object
   [A,b] = to_msparse(eqn);
   x = get(T);
   x = x';
   
   % Check tolerance and iteration count
   TRes = b-A*x;
   TResnorm = norm(TRes);         
   if TResnorm < casedef.iteration.TTol
      Tconverged = true;
      iterate = false;
   elseif niter > casedef.iteration.maxniter
      Tconverged = false;
      iterate = false;
%    elseif checkstoprequest(stopmon)
%       Tconverged = false;
%       iterate = false;
   else
      x = A\b; % Direct sparse solver.
               % Alternatives: gmres, bicgstabb, ...
      set(T,x'); % Put algebraic solution in the Field
   end
      
   
end % iterate

result.endtime = now; % call datestr(now) for displaying this time 
result.Tconverged = Tconverged;
result.niter = niter;
result.TResnorm = TResnorm;
result.TRes = Field(dom.allCells,0);
   set(result.TRes,TRes');
result.T = T;

  size(adiag)
   size(aNb_int)
   size(aNb_bound)
end



