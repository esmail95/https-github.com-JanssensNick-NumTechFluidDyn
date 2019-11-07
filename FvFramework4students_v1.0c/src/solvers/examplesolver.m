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

% Create field objects
T = Field(dom.allCells,0);      % Temperature [K] (scalar); empty field
reset(T,0);                     % Reset with all zeros
% Import velocities
U = casedef.vars.U.data;             

% Import variables, create data structures
nC = dom.nC;                    % Number of cells
nIf = dom.nIf;                  % Number of interior cells
cCoord = dom.cCoord;            % Cell coordinates

fCoord = dom.fCoord;            % Face coordinates
nF = dom.nF;                    % Number of faces
nBf = dom.nBf;                  % Number of boundary faces
fArea = dom.fArea;              % Face area
fNbC = dom.fNbC;                % Face neighbouring cells
fNormal = dom.fNormal;          % Normal on faces
Xi = dom.fXi;                   % Xi-vector on faces

k = casedef.material.k;         % Material property

BC = casedef.BC;                % Boundary conditions
nbZones = length(BC);           % Number of region boundaries
ranges = zeros(2,nbZones);      % Faces per boundary region
for i = 1:nbZones
    zone = BC{i}.zoneID;
    boundary = dom.getzone(zone);
    ranges(:,i) = boundary.range;
end

% Lambda
fXiLamba = dom.fXiLambda;

% Create an equation object for holding a scalar conservation equation
eqn = ScalarFvEqn2(dom);
%2 vgl u en v
%time = -dt;


iterate = true;
niter = 0;
while iterate   
   % Initialise matrics:
   % - A = [ap ; aNbInt ; aNbBound]
   ap = zeros(nC,1);
   aNbInt = zeros(2*nIf,1);
   aNbBound = zeros(2*nBf,1); 
   b = zeros(nC,1);
   
   % Update iteration counter
   niter = niter+1;
   %time = time + dt;
   
   % Set all terms to zero
   reset(eqn); 
   
   % Compute coefficients for physical cell eqns and add them to eqn object
   for i = 1:nIf
      nb1 = fNbC(2*i-1);
      nb2 = fNbC(2*i);
      %c1 = cCoord(:,nb1);
      %c2 = cCoord(:,nb2);
      %fc = fCoord(:,i);
      %l = lambda(c1,c2,fc); 
      l = fXiLamba(i);
      Af = fArea(i);
      Xif = norm(Xi(:,i));
      Uf = l*U(:,nb1) + (1-l)*U(:,nb2);
      Unf = dot(Uf,fNormal(:,i)); 
      %%% Compute sign of the normal
      %diff = cCoord(:,nb2)-cCoord(:,nb1);
      %outw = sign(dot(diff,fNormal(:,i)));
      outw = 1;
      
      %%% With diffusion in reversed direction and advection in the right
      %%% direction (very accurate result)
      ap(nb1) = ap(nb1) + k*Af/Xif + outw*l*Unf*Af; % advection + diffusion term
      ap(nb2) = ap(nb2) + k*Af/Xif - outw*l*Unf*Af;
      aNbInt(2*i-1) = aNbInt(2*i-1) - k*Af/Xif + outw*(1-l)*Unf*Af; % Upper diag, corresponds to nb1 contribution
      aNbInt(2*i) = aNbInt(2*i) - k*Af/Xif - outw*(1-l)*Unf*Af;
      
      % With diffusion in original direction and advection in the wrong
      % direction (okay result)
      %ap(nb1) = ap(nb1) - k*Af/Xif + outw*l*Unf*Af; % advection + diffusion term
      %ap(nb2) = ap(nb2) - k*Af/Xif - outw*l*Unf*Af;
      %aNbInt(2*i) = aNbInt(2*i) + k*Af/Xif + outw*(1-l)*Unf*Af; %
      %aNbInt(2*i-1) = aNbInt(2*i-1) + k*Af/Xif - outw*(1-l)*Unf*Af; %

      % With diffusion in original direction and advection in the right
      % direction (erronous result)
      %ap(nb1) = ap(nb1) - k*Af/Xif + outw*l*Unf*Af; % advection + diffusion term
      %ap(nb2) = ap(nb2) - k*Af/Xif - outw*l*Unf*Af;      
      %aNbInt(2*i-1) = aNbInt(2*i-1) + k*Af/Xif + outw*(1-l)*Unf*Af; % Upper diag, corresponds to nb1 contribution
      %aNbInt(2*i) = aNbInt(2*i) + k*Af/Xif - outw*(1-l)*Unf*Af; % Lower diag, corresponds to nb2 contribution
   end
   % Compute coefficients for ghost cell eqns and add them to eqn object
   for i = 1:nBf
      nb1 = fNbC(2*i-1 + 2*nIf); % Physicl cell
      nb2 = fNbC(2*i + 2*nIf); % Ghost Cell
      l = fXiLamba(i + nIf);
      Af = fArea(i + nIf);
      Xif = norm(Xi(:,i + nIf));
      Uf = l*U(:,nb1) + (1-l)*U(:,nb2);
      Unf = dot(Uf,fNormal(:,i + nIf));
      %%%
      %diff = cCoord(:,nb2)-cCoord(:,nb1);
      %outw = sign(dot(diff,fNormal(:,i + nIf))); %%%%%% CHANGED 19/10
      outw = 1;
      %%%
      boundaryFound = 0;
      for j = 1:nbZones
        % Determine to which boundary the face belongs
        if ismember(i + nIf,ranges(1,j):ranges(2,j)) && boundaryFound == 0
            boundaryFound = 1;
            switch BC{j}.kind
                case 'Dirichlet' 
                    %l*phi(GC) + (1-l)*phi(PC) = phi*
                    % Diagonal elements from BC: ghost cell
                    ap(nb2) = l;
                    aNbBound(2*i) = 1-l;
                    % Diagonal elements: physical cells
                    
                    % With inversed diffusion
                    ap(nb1) = ap(nb1) + k*Af/Xif + outw*l*Unf*Af;
                    aNbBound(2*i-1) = aNbBound(2*i-1) - k*Af/Xif + outw*(1-l)*Unf*Af;
                    
                    % With diffusion in original direction
                    %ap(nb1) = ap(nb1) - k*Af/Xif + outw*l*Unf*Af;
                    %aNbBound(2*i-1) = aNbBound(2*i-1) + k*Af/Xif + outw*(1-l)*Unf*Af;
                    
                    % Forcing terms
                    phi_star = BC{j}.data.bcval;
                    %b(i+nPc) = phi_star;  
                    b(nb2) = phi_star;
                case 'Neumann'
                    % (phi(PC)-phi(GC))/Ef=phi*
                    % Diagonal elements from BC: ghost cell
                    ap(nb2) = -1/Xif;
                    aNbBound(2*i) = 1/Xif;
                    
                    % Diagonal elements: physical cells
                    % With inversed diffusion
                    ap(nb1) = ap(nb1) + k*Af/Xif + outw*(1-l)*Unf*Af;
                    aNbBound(2*i-1) = aNbBound(2*i-1) - k*Af/Xif + outw*(1-l)*Unf*Af;
                    
                    % With original diffusion            
                    %ap(nb1) = ap(nb1) - k*Af/Xif + outw*(1-l)*Unf*Af;
                    %aNbBound(2*i-1) = aNbBound(2*i-1) + k*Af/Xif + outw*(1-l)*Unf*Af;
                    
                    % Forcing terms
                    phi_star = BC{j}.data.bcval;
                    b(nb2) = phi_star;  
            end
        end
      end
   end
   % eqn.adata = rand(eqn.nnz,1); % just a meaningless example adata
   % eqn.bdata = rand(eqn.n,1);   % just a meaningless example bdata
   eqn.adata = [ap ; aNbInt ; aNbBound];
   eqn.bdata = b;
   
   % Create a matlab sparse linear system from the eqn object
   [A,b] = to_msparse(eqn);
   %AA = full(A);
   %save('AA')
   %save('b')
   %spy(A)
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


end



