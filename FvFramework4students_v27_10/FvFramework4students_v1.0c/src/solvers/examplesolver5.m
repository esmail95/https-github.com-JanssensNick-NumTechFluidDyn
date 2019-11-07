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
function result = examplesolver5(casedef)

dom = casedef.dom;
material=casedef.material;

U = Field(casedef.dom.allCells,1);     % Velocity [m/s] (vector);
set(U,[0*ones(1,U.elcountzone);0*ones(1,U.elcountzone)]);

Up=casedef.vars.U.data;
Vp=Up;

Up_s =Field(dom.allCells,0);      
reset(Up_s,0); 
Vp_s = Field(dom.allCells,0);     
reset(Vp_s,0);



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
fCoord=dom.fCoord;
fNbV=dom.fNbV;
fNbVLoc=dom.fNbVLoc;
volPc=dom.cVol;
dt=casedef.iteration.dt;
Nx=casedef.vars.Nx;
Ny=casedef.vars.Ny;
Lx=casedef.vars.Lx;
Ly=casedef.vars.Ly;
% P_in=casedef.vars.P_in;
% P_out=casedef.vars.P_out;


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
  
 Source = casedef.vars.Source;
dx=Lx/Nx;
dy=Ly/Ny;
% dPx=(P_out-P_in)/Nx;
% dPy=0;

x_source=cCoord(1,:);
y_source=cCoord(2,:);
S = Source(x_source,y_source);

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
eqn_v = ScalarFvEqn2(dom);
iterate = true;
time = true;
niter = 0;
t = 0;
    U_old=Up;
    V_old=Vp;
    
% % Import variables, create data structures
% nC = dom.nC;                    % Number of cells
% nPc = dom.nPc;                  % Number of physical cells
% nIf = dom.nIf;                  % Number of interior cells
% cCoord = dom.cCoord;            % Cell coordinates
% cVol = dom.cVol(1:dom.nPc);     % Cell volumes for physical cells
% 
% fCoord = dom.fCoord;            % Face coordinates
% nF = dom.nF;                    % Number of faces
% nBf = dom.nBf;                  % Number of boundary faces
% fArea = dom.fArea;              % Face area
% fNbC = dom.fNbC;                % Face neighbouring cells
% fNormal = dom.fNormal;          % Normal on faces
% Xi = dom.fXi;                   % Xi-vector on faces
% 
% k = casedef.material.k;         % Material property
% rho = casedef.material.rho;     % Material property
% 
% BC = casedef.BC;                % Boundary conditions
% nbZones = length(BC);           % Number of region boundaries
% ranges = zeros(2,nbZones);      % Faces per boundary region
% for i = 1:nbZones
%     zone = BC{i}.zoneID;
%     boundary = dom.getzone(zone);
%     ranges(:,i) = boundary.range;
% end
% 
% %dX = casedef.vars.dX;
% %dY = casedef.vars.dY;
% %pIn = casedef.vars.pIn;
% %pOut = casedef.vars.pOut;
% %dPx = (pIn-pOut)/dX;
% %dPy = (pIn-pOut)/dY;
% %dPy = 0;
% Source = casedef.vars.Source;
% 
% 
% dt = casedef.iteration.dt;
% 
% % Lambda
% fXiLamba = dom.fXiLambda;
% 
% % Create an equation object for holding a scalar conservation equation
% eqn_u = ScalarFvEqn2(dom);
% eqn_v = ScalarFvEqn2(dom);
% time = -dt;
% stepping = true;
% 
% iterate = true;
% niter = 0;

while time
    t = t + dt;
    
    %FOLLOWING CYCLE IS SOLVING AT EACH FALSE TIME STEPPING
while iterate   
   % Set all terms to zero
   adiag_u=zeros(nC,1);
   aNb_int_u=zeros(2*nIf,1);
   aNb_bound_u=zeros(2*nBf,1);
   b_u=zeros(nC,1);

   adiag_v=zeros(nC,1);
   aNb_int_v=zeros(2*nIf,1);
   aNb_bound_v=zeros(2*nBf,1);
   b_v=zeros(nC,1); 
    
   niter = niter+1;
   
   reset(eqn_u);
   reset(eqn_v);
   
  
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
      adiag_u(nb1) = adiag_u(nb1) - k*Af/Xi + outw*Af*U_n*XiLambda + volPc(nb1)/(4*dt);
      adiag_u(nb2) = adiag_u(nb2) - k*Af/Xi - outw*Af*U_n*(1-XiLambda) + volPc(nb2)/(4*dt);
      % Off-diagonal elements for internal faces
      aNb_int_u(2*i) = aNb_int_u(2*i) + k*Af/Xi + outw*Af*U_n*XiLambda; % Benedendiagonaal
      aNb_int_u(2*i-1) = aNb_int_u(2*i-1) + k*Af/Xi - outw*Af*U_n*(1-XiLambda); % Bovendiagonaal
      
      %bVector
      b_u(nb1) =  S(1,nb1) * volPc(nb1) + U_old(1,nb1) * volPc(nb1)/dt;
      b_u(nb2) =  S(1,nb2) * volPc(nb2) + U_old(1,nb2) * volPc(nb2)/dt;
   end
   % Compute coefficients for ghost cell eqns and add them to eqn object
   for i = 1:nBf
      fCx = fCoord(1,i + nIf);
      fCy = fCoord(2,i + nIf); 
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
                    adiag_u(nb1) = adiag_u(nb1) - k*Af/Xi + outw*Af*U_n*(1-XiLambda)+volPc(nb1)/(4*dt);
                    % Off-diagonal elements: physical cell
                    aNb_bound_u(2*i) = 1-XiLambda;
                    % Off-diagonal elements: ghost cells
                    aNb_bound_u(2*i-1) = aNb_bound_u(2*i-1) + k*Af/Xi + outw*Af*U_n*(1-XiLambda);
                    
                    % Forcing terms
                    %for Pc already accessed
                    %by Pc
                    %bcell_forcing_u = - dPx/dx * volPc(nb1) / rho + U_old(1,nb1) * volPc(nb1)/dt;
                    phi_star_u = BC{j}.data.bcval(fCx,fCy);
                    %b(i+nPc) = phi_star;
                    %b_u(nb1) = bcell_forcing_u;
                    b_u(nb2) = phi_star_u(1);
                    
                case 'Neumann'
                    % (phi(PC)-phi(GC))/Ef=phi*
                    % Diagonal elements from BC: ghost cell
                    adiag_u(nb2) = -1/Xi;
                    % Diagonal elements: physical cells
                    adiag_u(nb1) = adiag_u(nb1) - k*Af/Xi + outw*Af*U_n*(1-XiLambda)+volPc(nb1)/(4*dt);
                    
                    % Off-diagonal elements: physical cell
                    aNb_bound_u(2*i) = 1/Xi;
                    % Off-diagonal elements: ghost cells
                    aNb_bound_u(2*i-1) = aNb_bound_u(2*i-1) + k*Af/Xi + outw*Af*U_n*(1-XiLambda);
                    
                    % Forcing terms
                    %bcell_forcing_u = - dPx/dx * volPc(nb1) / rho +
                    %U_old(1,nb1) * volPc(nb1)/dt; for Pc already accessed
                    %by If
                    phi_star_u = BC{j}.data.bcval;
                    %b_u(nb1) = bcell_forcing_u;
                    b_u(nb2) = phi_star_u;  
            end
        end
      end
   end
   
   
 %V COMPONENT  
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
      V_f=XiLambda*V_old(:,nb1)+(1-XiLambda)*V_old(:,nb2);
      V_n=dot(V_f,Normal); 
      
      %%% Compute sign of the normal
      diff = cCoord(:,nb2)-cCoord(:,nb1);
      outw = sign(dot(diff,fNormal(:,i)));
      
      % Diagonal elements
      adiag_v(nb1) = adiag_v(nb1) - k*Af/Xi + outw*Af*V_n*XiLambda + volPc(nb1)/(4*dt);
      adiag_v(nb2) = adiag_v(nb2) - k*Af/Xi - outw*Af*V_n*(1-XiLambda) + volPc(nb2)/(4*dt);
      % Off-diagonal elements for internal faces
      aNb_int_v(2*i) = aNb_int_v(2*i) + k*Af/Xi + outw*Af*V_n*XiLambda; % Benedendiagonaal
      aNb_int_v(2*i-1) = aNb_int_v(2*i-1) + k*Af/Xi - outw*Af*V_n*(1-XiLambda); % Bovendiagonaal
      
      %bVector
      b_v(nb1) =  S(2,nb1) * volPc(nb1) + V_old(2,nb1) * volPc(nb1)/dt;
      b_v(nb2) =  S(2,nb2) * volPc(nb2) + V_old(2,nb2) * volPc(nb2)/dt;
   end
   % Compute coefficients for ghost cell eqns and add them to eqn object
   for i = 1:nBf
      fCx = fCoord(1,i + nIf);
      fCy = fCoord(2,i + nIf); 
      nb1 = fNbC(2*i-1 + 2*nIf); % Physicl cell
      nb2 = fNbC(2*i + 2*nIf); % Ghost Cell
      Af = fArea(i + nIf);
      Xi = norm(fXi(:,i + nIf));
      Normal=fNormal(:,i+nIf);
      %fXiLambda weighted average between main cell and neighbour to find phi_f  
      XiLambda=fXiLambda(i+nIf);
      %advective term
      V_f=XiLambda*V_old(:,nb1)+(1-XiLambda)*V_old(:,nb2);
      V_n=(dot(V_f,Normal));
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
                    adiag_v(nb2) = XiLambda;
                    % Diagonal elements: physical cells
                    adiag_v(nb1) = adiag_v(nb1) - k*Af/Xi + outw*Af*V_n*(1-XiLambda)+ V_old(2,nb1) * volPc(nb1)/(4*dt);
                    % Off-diagonal elements: physical cell
                    aNb_bound_v(2*i) = 1-XiLambda;
                    % Off-diagonal elements: ghost cells
                    aNb_bound_v(2*i-1) = aNb_bound_v(2*i-1) + k*Af/Xi + outw*Af*V_n*(1-XiLambda);
                    
                    % Forcing terms
                    %bcell_forcing_v = - dPy/dy * volPc(nb1) / rho + V_old(2,nb1) * volPc(nb1)/dt;
                    phi_star_v = BC{j}.data.bcval(fCx,fCy);
                    %b(i+nPc) = phi_star;
                    %b_v(nb1) = bcell_forcing_v;
                    b_v(nb2) = phi_star_v(2);
                    
                case 'Neumann'
                    % (phi(PC)-phi(GC))/Ef=phi*
                    % Diagonal elements from BC: ghost cell
                    adiag_v(nb2) = -1/Xi;
                    % Diagonal elements: physical cells
                    adiag_v(nb1) = adiag_v(nb1) - k*Af/Xi + outw*Af*V_n*(1-XiLambda)+ V_old(2,nb1) * volPc(nb1)/(4*dt);
                    
                    % Off-diagonal elements: physical cell
                    aNb_bound_v(2*i) = 1/Xi;
                    % Off-diagonal elements: ghost cells
                    aNb_bound_v(2*i-1) = aNb_bound_v(2*i-1) + k*Af/Xi + outw*Af*V_n*(1-XiLambda);
                    
                    % Forcing terms
                    %bcell_forcing_v = - dPy/dy * volPc(nb1) / rho + V_old(2,nb1) * volPc(nb1)/dt;
                    phi_star_v = BC{j}.data.bcval_v;
                    %b_v(nb1) = bcell_forcing_v;
                    b_v(nb2) = phi_star_v;  
            end
        end
      end
   end
   % eqn.adata = rand(eqn.nnz,1); % just a meaningless example adata
   % eqn.bdata = rand(eqn.n,1);   % just a meaningless example bdata
   eqn_u.adata = [adiag_u ; aNb_int_u ; aNb_bound_u];
   eqn_u.bdata = b_u;
   
   eqn_v.adata = [adiag_v ; aNb_int_v ; aNb_bound_v];
   eqn_v.bdata = b_v;
     
   % Create a matlab sparse linear system from the eqn object
   [A_u,b_u] = to_msparse(eqn_u);
   x = get(Up_s);
   x = x';
   
   [A_v,b_v] = to_msparse(eqn_v);
   y = get(Vp_s);
   y = y';
   
   % Check tolerance and iteration count
   res_u = b_u-A_u*x;
   resnorm_u = norm(res_u);
   
   res_v = b_v-A_v*y;
   resnorm_v = norm(res_v); 
   
   if max(resnorm_v,resnorm_u) < casedef.iteration.TTol
      u_v_converged = true;
      iterate = false;
   elseif niter > casedef.iteration.maxniter
      u_v_converged = false;
      iterate = false;
      
%    elseif checkstoprequest(stopmon)
%       Tconverged = false;
%       iterate = false;
   else
      x = A_u\b_u; % Direct sparse solver.
               % Alternatives: gmres, bicgstabb, ...
      set(Up_s,x'); % Put algebraic solution in the Field
      
      y = A_v\b_v; % Direct sparse solver.
               % Alternatives: gmres, bicgstabb, ...
      set(Vp_s,y');
   end
   
      
end % iterate

if max(norm([x';zeros(1,nC)]-U_old),norm([zeros(1,nC);y']-V_old))< casedef.iteration.TTol
      u_v_converged_t = true;
      time = false;
else   
       
       U_old=[x';zeros(1,nC)];
       V_old=[zeros(1,nC);y'];
       time=true;
       
end

end %time iteration
 t

set(U,[x';y']);
result.endtime = now; % call datestr(now) for displaying this time 
result.Tconverged_u_v = u_v_converged;
result.niter = niter;
result.resnorm_u = resnorm_u;
result.resnorm_v = resnorm_v;
result.res_u = Field(dom.allCells,0);
   set(result.res_u,res_u');
result.res_v = Field(dom.allCells,0);
   set(result.res_v,res_v');
   
result.U = U;



%   size(adiag_u)
%    size(aNb_int_u)
%    size(aNb_bound_u)
end
