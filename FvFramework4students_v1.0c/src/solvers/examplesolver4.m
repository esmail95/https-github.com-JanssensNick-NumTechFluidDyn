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
function result = examplesolver4(casedef)

dom = casedef.dom;

% Create field objects
U = Field(dom.allCells,1);      % Velocity field [m/s]
U_old = Field(dom.allCells,1); 
set(U,casedef.vars.U.data);
set(U_old,casedef.vars.U.data);
%U = Field(dom.allCells,0);      % Velocity in x [m/s]
%V = Field(dom.allCells,0);      % Velocity in y [m/s]
%reset(U,0);
%reset(V,0);     

% Import variables, create data structures
nC = dom.nC;                    % Number of cells
nPc = dom.nPc;                  % Number of physical cells
nIf = dom.nIf;                  % Number of interior cells
cCoord = dom.cCoord;            % Cell coordinates
cVol = dom.cVol(1:dom.nPc);     % Cell volumes for physical cells

fCoord = dom.fCoord;            % Face coordinates
nF = dom.nF;                    % Number of faces
nBf = dom.nBf;                  % Number of boundary faces
fArea = dom.fArea;              % Face area
fNbC = dom.fNbC;                % Face neighbouring cells
fNormal = dom.fNormal;          % Normal on faces
Xi = dom.fXi;                   % Xi-vector on faces

k = casedef.material.k;         % Material property
rho = casedef.material.rho;     % Material property

BC = casedef.BC;                % Boundary conditions
nbZones = length(BC);           % Number of region boundaries
ranges = zeros(2,nbZones);      % Faces per boundary region
for i = 1:nbZones
    zone = BC{i}.zoneID;
    boundary = dom.getzone(zone);
    ranges(:,i) = boundary.range;
end

dX = casedef.vars.dX;
dY = casedef.vars.dY;
%cos = casedef.vars.cosTheta;
%sin = casedef.vars.sinTheta;
pIn = casedef.vars.pIn;
pOut = casedef.vars.pOut;
dPx = -(pIn-pOut)/dX;
dPy = -(pIn-pOut)/dY;
dPy = 0;

dt = casedef.iteration.dt;

% Lambda
fXiLamba = dom.fXiLambda;

% Create an equation object for holding a scalar conservation equation
eqn_u = ScalarFvEqn2(dom);
eqn_v = ScalarFvEqn2(dom);
time = -dt;
stepping = true;

niter = 0;

while stepping
    time = time + dt;
    iterate = true;

while iterate  
   % Initialise matrics:
   % - Au = [apu ; aNbIntu ; aNbBoundu]
   % - Av = [apv ; aNbIntv ; aNbBoundv]
   apu = zeros(nC,1);
   aNbIntu = zeros(2*nIf,1);
   aNbBoundu = zeros(2*nBf,1); 
   bu = zeros(nC,1);
   apv = zeros(nC,1);
   aNbIntv = zeros(2*nIf,1);
   aNbBoundv = zeros(2*nBf,1); 
   bv = zeros(nC,1);
   
   % Update iteration counter
   niter = niter+1;
   
   % Set all terms to zero
   reset(eqn_u);
   reset(eqn_v);
   % Compute coefficients for physical cell eqns and add them to eqn object
   for i = 1:nPc
      Ui = U_old.data(1,i); % Current u
      Vi = U_old.data(2,i); % Current v
      apu(i) = apu(i) + cVol(i)/dt;
      apv(i) = apv(i) + cVol(i)/dt;
      bu(i) = bu(i) - cVol(i)*dPx/rho + cVol(i)*Ui/dt;
      bv(i) = bv(i) - cVol(i)*dPy/rho + cVol(i)*Vi/dt;
   end
   for i = 1:nIf
      nb1 = fNbC(2*i-1);
      nb2 = fNbC(2*i);
      l = fXiLamba(i);
      Af = fArea(i);
      Xif = norm(Xi(:,i));
      U1 = U_old.data(:,nb1);
      U2 = U_old.data(:,nb2);
      Uf = l*U1 + (1-l)*U2;
      Unf = dot(Uf,fNormal(:,i));
      %%% Compute sign of the normal
      %diff = cCoord(:,nb2)-cCoord(:,nb1);
      %outw = sign(dot(diff,fNormal(:,i)));
      outw = 1;
      %%%
      
      %%% With diffusion in reversed direction and advection in the right
      %%% direction (very accurate result)
      apu(nb1) = apu(nb1) + outw*l*Unf*Af + k*Af/Xif; 
      apu(nb2) = apu(nb2) - outw*l*Unf*Af + k*Af/Xif;
      apv(nb1) = apv(nb1) + outw*l*Unf*Af + k*Af/Xif; 
      apv(nb2) = apv(nb2) - outw*l*Unf*Af + k*Af/Xif;
      aNbIntu(2*i-1) = aNbIntu(2*i-1) + outw*(1-l)*Unf*Af  - k*Af/Xif; % Bovendiagonaal,nb1
      aNbIntu(2*i) = aNbIntu(2*i) - outw*(1-l)*Unf*Af - k*Af/Xif; % Benedendiagonaal, nb2
      aNbIntv(2*i-1) = aNbIntv(2*i-1) + outw*(1-l)*Unf*Af  - k*Af/Xif; % Bovendiagonaal, nb1
      aNbIntv(2*i) = aNbIntv(2*i) - outw*(1-l)*Unf*Af - k*Af/Xif; % Benedendiagonaal, nb2
      
      % With diffusion in original direction and advection in the wrong
      % direction (okay result)
%       apu(nb1) = apu(nb1) + outw*l*Unf*Af - k*Af/Xif; 
%       apu(nb2) = apu(nb2) - outw*l*Unf*Af - k*Af/Xif;
%       apv(nb1) = apv(nb1) + outw*l*Unf*Af - k*Af/Xif; 
%       apv(nb2) = apv(nb2) - outw*l*Unf*Af - k*Af/Xif;
%       aNbIntu(2*i-1) = aNbIntu(2*i-1) - outw*(1-l)*Unf*Af  + k*Af/Xif; % Bovendiagonaal,nb1
%       aNbIntu(2*i) = aNbIntu(2*i) + outw*(1-l)*Unf*Af + k*Af/Xif; % Benedendiagonaal, nb2
%       aNbIntv(2*i-1) = aNbIntv(2*i-1) - outw*(1-l)*Unf*Af  + k*Af/Xif; % Bovendiagonaal, nb1
%       aNbIntv(2*i) = aNbIntv(2*i) + outw*(1-l)*Unf*Af + k*Af/Xif; % Benedendiagonaal, nb2
      
      % With diffusion in original direction and advection in the right
      % direction (erronous result)
      %apu(nb1) = apu(nb1) + outw*l*Unf*Af - k*Af/Xif; 
      %apu(nb2) = apu(nb2) - outw*l*Unf*Af - k*Af/Xif;
      %apv(nb1) = apv(nb1) + outw*l*Unf*Af - k*Af/Xif; 
      %apv(nb2) = apv(nb2) - outw*l*Unf*Af - k*Af/Xif;
      %aNbIntu(2*i-1) = aNbIntu(2*i-1) + outw*(1-l)*Unf*Af  + k*Af/Xif; % Bovendiagonaal,nb1
      %aNbIntu(2*i) = aNbIntu(2*i) - outw*(1-l)*Unf*Af + k*Af/Xif; % Benedendiagonaal, nb2
      %aNbIntv(2*i-1) = aNbIntv(2*i-1) + outw*(1-l)*Unf*Af  + k*Af/Xif; % Bovendiagonaal, nb1
      %aNbIntv(2*i) = aNbIntv(2*i) - outw*(1-l)*Unf*Af + k*Af/Xif; % Benedendiagonaal, nb2
   end
   % Compute coefficients for ghost cell eqns and add them to eqn object
   for i = 1:nBf
      nb1 = fNbC(2*i-1 + 2*nIf); % Physicl cell
      nb2 = fNbC(2*i + 2*nIf); % Ghost Cell
      l = fXiLamba(i + nIf);
      Af = fArea(i + nIf);
      Xif = norm(Xi(:,i + nIf));
      U1 = U_old.data(:,nb1);
      U2 = U_old.data(:,nb2);
      Uf = l*U1 + (1-l)*U2;
      Unf = dot(Uf,fNormal(:,i + nIf));
      %%%
      diff = cCoord(:,nb2)-cCoord(:,nb1);
      outw = sign(dot(diff,fNormal(:,i + nIf))); %%% CHANGED 19/10 +nIf
      %%%
      boundaryFound = 0;
      for j = 1:nbZones
        % Determine to which boundary the face belongs
        if ismember(i + nIf,ranges(1,j):ranges(2,j)) && boundaryFound == 0
            boundaryFound = 1;
            switch BC{j}.kind_u
                case 'Dirichlet' 
                    % Diagonal elements from BC: ghost cell
                    apu(nb2) = l;
                    % Off-diagonal elements: physical cell
                    aNbBoundu(2*i) = 1-l;
                    
                    % With inversed diffusion
                    apu(nb1) = apu(nb1) + outw*l*Unf*Af + k*Af/Xif;
                    aNbBoundu(2*i-1) = aNbBoundu(2*i-1) + outw*(1-l)*Unf*Af - k*Af/Xif;
                    
                    % With original diffusion
                    % Diagonal elements: physical cells
                    %apu(nb1) = apu(nb1) + outw*l*Unf*Af - k*Af/Xif;
                    % Off-diagonal elements: ghost cells
                    %aNbBoundu(2*i-1) = aNbBoundu(2*i-1) + outw*(1-l)*Unf*Af + k*Af/Xif;
                    
                    % Forcing terms
                    phi_star = BC{j}.data.bcval_u;
                    bu(nb2) = phi_star;
                    
                    % Diagonal elements from BC: ghost cell
                    apv(nb2) = l;
                    % Off-diagonal elements: physical cell
                    aNbBoundv(2*i) = 1-l;
                    
                    % With inversed diffusion
                    apv(nb1) = apv(nb1) + outw*l*Unf*Af + k*Af/Xif;
                    aNbBoundv(2*i-1) = aNbBoundv(2*i-1) + outw*(1-l)*Unf*Af - k*Af/Xif;
                    
                    % With original diffusion
                    % Diagonal elements: physical cells
                    %apv(nb1) = apv(nb1) + outw*l*Unf*Af - k*Af/Xif;
                    % Off-diagonal elements: ghost cells
                    %aNbBoundv(2*i-1) = aNbBoundv(2*i-1) + outw*(1-l)*Unf*Af + k*Af/Xif;
                    
                    % Forcing terms
                    phi_star = BC{j}.data.bcval_v; 
                    bv(nb2) = phi_star;
                case 'Neumann'
                    % (phi(PC)-phi(GC))/Ef=phi*
                    % Diagonal elements from BC: ghost cell
                    apu(nb2) = -1/Xif;
                    % Off-diagonal elements: physical cell
                    aNbBoundu(2*i) = 1/Xif;
                    
                    % With inversed diffusion
                    apu(nb1) = apu(nb1) + outw*(1-l)*Unf*Af + k*Af/Xif;
                    aNbBoundu(2*i-1) = aNbBoundu(2*i-1) - k*Af/Xif + outw*(1-l)*Unf*Af;
                    
                    % With original diffusion
                    % Diagonal elements: physical cells
                    %apu(nb1) = apu(nb1) + outw*(1-l)*Unf*Af - k*Af/Xif;
                    % Off-diagonal elements: ghost cells
                    %aNbBoundu(2*i-1) = aNbBoundu(2*i-1) + k*Af/Xif + outw*(1-l)*Unf*Af;
                    
                    % Forcing terms
                    phi_star = BC{j}.data.bcval_u;
                    bu(nb2) = phi_star;  
                    
                    % Diagonal elements from BC: ghost cell
                    apv(nb2) = -1/Xif;
                    % Off-diagonal elements: physical cell
                    aNbBoundv(2*i) = 1/Xif;
                    
                    % With reversed diffusion
                    apv(nb1) = apv(nb1) + outw*(1-l)*Unf*Af + k*Af/Xif;
                    aNbBoundv(2*i-1) = aNbBoundv(2*i-1) - k*Af/Xif + outw*(1-l)*Unf*Af;
                    
                    % With original diffusion
                    % Diagonal elements: physical cells
                    %apv(nb1) = apv(nb1) + outw*(1-l)*Unf*Af - k*Af/Xif;
                    % Off-diagonal elements: ghost cells
                    %aNbBoundv(2*i-1) = aNbBoundv(2*i-1) + k*Af/Xif + outw*(1-l)*Unf*Af;
                    
                    % Forcing terms
                    phi_star = BC{j}.data.bcval_v;
                    bv(nb2) = phi_star; 
            end
        end
      end
   end
   eqn_u.adata = [apu ; aNbIntu ; aNbBoundu];
   eqn_u.bdata = bu;
   eqn_v.adata = [apv ; aNbIntv ; aNbBoundv];
   eqn_v.bdata = bv; 
   
   % Create a matlab sparse linear system from the eqn object
   [Au,bu] = to_msparse(eqn_u);
   [Av,bv] = to_msparse(eqn_v);
   x = get(U);
   xu = x(1,:);
   xv = x(2,:);
   xu = xu';
   xv = xv';
    %[Au,bu] = to_msparse(eqn_u);
%    AA = full(Au);
%    save('AA')
%    save('bu')
%    spy(AA)
   %%%
%    axis off; axis equal; colormap(jet(50));
%    colorbar
%    scale = 'lin'; lw = 1;
%    fvmplotfield(U,scale,lw,1); % 1 plots u, 2 plots v
%    fvmplotmesh(casedef.dom,lw);
   %%%
   % Check tolerance and iteration count
   URes = bu-Au*xu;
   VRes = bv-Av*xv;
   UResnorm = norm(URes);
   VResnorm = norm(VRes);
   Resnorm = max(UResnorm,VResnorm);      
   if Resnorm < casedef.iteration.tol
      %converged = true;
      iterate = false;
   elseif niter > casedef.iteration.maxniter
      %converged = false;
      iterate = false;
%    elseif checkstoprequest(stopmon)
%       converged = false;
%       iterate = false;
   else
      %xu = Au\bu; % Direct sparse solver.
               % Alternatives: gmres, bicgstabb, ...
      %xv = Av\bv;
      x = [Au\bu,Av\bv];
      set(U,x'); % Put algebraic solution in the Field
   end
end % iterate
if max(norm(U.data(1,:)-U_old.data(1,:)),norm(U.data(2,:)-U_old.data(2,:))) < casedef.iteration.tol
    stepping = false;
else
    stepping = true;
    set(U_old,U.data);
end


end % stepping

result.endtime = now; % call datestr(now) for displaying this time 
result.converged = 1-stepping;
%result.niter = niter;
result.steps = time/50;
result.UResnorm = UResnorm;
result.URes = Field(dom.allCells,0);
   set(result.URes,URes');
result.VResnorm = VResnorm;
result.VRes = Field(dom.allCells,0);
   set(result.VRes,VRes');
result.U = U;
result.time = time;
end



