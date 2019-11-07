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
function result_u = examplesolver(casedef)

dom = casedef.dom;
material=casedef.material;
% Create field objects
%T = Field(dom.allCells,0);      % Temperature [K] (scalar); empty field
%reset(T,0);                     % Reset with all zeros

U = casedef.vars.U.data;
Up_s =Field(dom.allCells,0);      
reset(Up_s,0); 

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
volPc=dom.cVol;
dt=casedef.iteration.dt;
Nx=casedef.vars.Nx;
Ny=casedef.vars.Ny;
Lx=casedef.vars.Lx;
Ly=casedef.vars.Ly;
P_in=casedef.vars.P_in;
P_out=casedef.vars.P_out;


k=material.k;
rho = casedef.material.rho;
fArea = calc_fArea(nF,dom.fNbVLoc,dom.fNbV,dom.vCoord);
fXi = calc_fXi(nF,dom.cCoord,dom.fNbCLoc,fNbC);
fNormal = calc_fNormal(fNbVLoc,fNbV,vCoord,fArea);
fTangent = calc_fTangent(nF,fNormal);
fXiLambda = calc_fXiLambda(...
      nF,...
      fXi,fTangent,cCoord,vCoord,...
      fNbCLoc,fNbC,fNbVLoc,fNbV);
  

adiag_u=zeros(nC,1);
aNb_int_u=zeros(2*nIf,1);
aNb_bound_u=zeros(2*nBf,1);
b_u=zeros(nC,1);


Up=[zeros(1,nC);zeros(1,nC)];

t=0;
dx=Lx/Nx;
dy=Ly/Ny;
dPx=(P_out-P_in)/Nx;
dPy=0;
resnorm_u_t=5;
%FOLLOWING CYCLE IS FALSE TIME STEPPING
% while resnorm_u_t>= casedef.iteration.TTol
% adiag_u=zeros(nC,1);
% aNb_int_u=zeros(2*nIf,1);
% aNb_bound_u=zeros(2*nBf,1);
% b_u=zeros(nC,1);

    U_old=U;    
    t = t + 1;
    
BC = casedef.BC;                % Boundary conditions
nbZones = length(BC);           % Number of region boundaries
ranges = zeros(2,nbZones);      % Faces per boundary region
for i = 1:nbZones
    zone = BC{i}.zoneID;
    boundary = dom.getzone(zone);
    ranges(:,i) = boundary.range;
end

% Create an equation object for holding a scalar conservation equation
%eqn = ScalarFvEqn2(dom);
eqn_u = ScalarFvEqn2(dom);
iterate = true;
niter = 0;

%FOLLOWING CYCLE IS SOLVING AT EACH FALSE TIME STEPPING
while iterate   
   
   adiag_u=zeros(nC,1);
   aNb_int_u=zeros(2*nIf,1);
   aNb_bound_u=zeros(2*nBf,1);
   b_u=zeros(nC,1);

   niter = niter+1;
   
% Set all terms to zero
%reset(eqn); 
reset(eqn_u);


% dx=seedI.displ.x;
% dy=seedJ.displ.y;

  %U COMPONENT 
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
      U_f=XiLambda*U_old(:,nb1)+(1-XiLambda)*U_old(:,nb2);
      U_n=dot(U_f,Normal); 
      %%% Compute sign of the normal
      diff = cCoord(:,nb2)-cCoord(:,nb1);
      outw = sign(dot(diff,fNormal(:,i)));
      % Diagonal elements
      adiag_u(nb1) = adiag_u(nb1) - k*Af/Xi + outw*Af*U_n*XiLambda + volPc(nb1)/dt;
      adiag_u(nb2) = adiag_u(nb2) - k*Af/Xi - outw*Af*U_n*(1-XiLambda) + volPc(nb2)/dt;
      % Off-diagonal elements for internal faces
      aNb_int_u(2*i) = aNb_int_u(2*i) + k*Af/Xi + outw*Af*U_n*XiLambda; % Benedendiagonaal
      aNb_int_u(2*i-1) = aNb_int_u(2*i-1) + k*Af/Xi - outw*Af*U_n*(1-XiLambda); % Bovendiagonaal
      
      %bVector
      b(nb1) = - dPx/dx * volPc(nb1) / rho + U_old(1,nb1) * volPc(nb1)/dt;
      b(nb2) = - dPx/dx * volPc(nb2) / rho + U_old(1,nb2) * volPc(nb2)/dt;
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
      U_f=XiLambda*U_old(:,nb1)+(1-XiLambda)*U_old(:,nb2);
      U_n=(dot(U_f,Normal));
      
      %Outward normal
      diff = cCoord(:,nb2)-cCoord(:,nb1);
      outw = sign(dot(diff,fNormal(:,i)));
      
      boundaryFound = 0;
      for j = 1:nbZones
        % Determine to which boundary the face belongs
        if ismember(i + nIf,ranges(1,j):ranges(2,j)) && boundaryFound == 0
            
            switch BC{j}.kind
                case 'Dirichlet' 
                    %l*phi(GC) + (1-l)*phi(PC) = phi*
                    % Diagonal elements from BC: ghost cell
                    adiag_u(nb2) = XiLambda;
                    % Diagonal elements: physical cells
                    adiag_u(nb1) = adiag_u(nb1) - k*Af/Xi + outw*Af*U_n*(1-XiLambda);
                    % Off-diagonal elements: physical cell
                    aNb_bound_u(2*i) = 1-XiLambda;
                    % Off-diagonal elements: ghost cells
                    aNb_bound_u(2*i-1) = aNb_bound_u(2*i-1) + k*Af/Xi + outw*Af*U_n*(1-XiLambda);
                    
                    % Forcing terms
                    bcell_forcing_u = - dPx/dx * volPc(nb1) / rho + U_old(1,nb1) * volPc(nb1)/dt;
                    phi_star_u = BC{j}.data.bcval_u;
                    %b(i+nPc) = phi_star;
                    b(nb1) = bcell_forcing_u;
                    b(nb2) = phi_star_u;
                    
                case 'Neumann'
                    % (phi(PC)-phi(GC))/Ef=phi*
                    % Diagonal elements from BC: ghost cell
                    adiag_u(nb2) = -1/Xi;
                    % Diagonal elements: physical cells
                    adiag_u(nb1) = adiag_u(nb1) - k*Af/Xi + outw*Af*U_n*(1-XiLambda);
                    
                    % Off-diagonal elements: physical cell
                    aNb_bound_u(2*i) = 1/Xi;
                    % Off-diagonal elements: ghost cells
                    aNb_bound_u(2*i-1) = aNb_bound_u(2*i-1) + k*Af/Xi + outw*Af*U_n*(1-XiLambda);
                    
                    % Forcing terms
                    bcell_forcing_u = - dPx/dx * volPc(nb1) / rho + U_old(1,nb1) * volPc(nb1)/dt;
                    phi_star_u = BC{j}.data.bcval_u;
                    b(nb1) = bcell_forcing_u;
                    b(nb2) = phi_star_u;  
            end
        end
      end
   end
   
   
  % eqn.adata = rand(eqn.nnz,1); % just a meaningless example adata
   % eqn.bdata = rand(eqn.n,1);   % just a meaningless example bdata
   eqn_u.adata = [adiag_u ; aNb_int_u ; aNb_bound_u];
   eqn_u.bdata = b_u
   
       
   % Create a matlab sparse linear system from the eqn object
   [A_u,b_u] = to_msparse(eqn_u);
   x = get(Up_s);
   x = x';
   
     
   % Check tolerance and iteration count
   res_u = b_u-A_u*x;
   norm_u = norm(res_u);
   
   
   if norm_u < casedef.iteration.TTol
      u_converged = true;
      iterate = false;
   elseif niter > casedef.iteration.maxniter
      u_converged = false;
      iterate = false;
      
%    elseif checkstoprequest(stopmon)
%       Tconverged = false;
%       iterate = false;
   else
      x = A_u\b_u; % Direct sparse solver.
               % Alternatives: gmres, bicgstabb, ...
      set(Up_s,x'); % Put algebraic solution in the Field
      
    end
      
   
end % iterate

%NOW RESIDUALS IN TIME
 

   

 %end

% result_u_v.endtime = now; % call datestr(now) for displaying this time 
% result_u_v.Tconverged_u_v = u_v_converged_t;
% result_u_v.niter = niter;
% result_u_v.resnorm_u_t = resnorm_u_t;
% result_u_v.resnorm_v_t = resnorm_v_t;
% result_u_v.res_u_t = Field(dom.allCells,0);
%    set(result_u_v.res_u_t,res_u_t');
% result_u_v.res_v_t = Field(dom.allCells,0);
%    set(result_u_v.res_v_t,res_v_t');
result_u.Up_s = Up_s;





%   size(adiag_u)
%    size(aNb_int_u)
%    size(aNb_bound_u)
end



