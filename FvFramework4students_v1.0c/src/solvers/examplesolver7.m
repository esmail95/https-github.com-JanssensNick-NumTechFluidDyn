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
function result = examplesolver6bis(casedef)

dom = casedef.dom;

% Create field objects
U = Field(dom.allCells,1);      % Velocity field [m/s]
U_old = Field(dom.allCells,1);  % Previous iteration velocity field [m/s]
set(U,casedef.vars.U.data);
set(U_old,casedef.vars.U.data);     
P = Field(dom.allCells,0);      % Pressure field [Pa]
P_prime = Field(dom.allCells,0);% Pressure field corrections [Pa]
reset(P_prime,0);
set(P,casedef.vars.P.data);

% Import variables, create data structures
nC = dom.nC;                    % Number of cells
nPc = dom.nPc;                  % Number of physical cells
nIf = dom.nIf;                  % Number of interior cells
cCoord = dom.cCoord;            % Cell coordinates
cVol = dom.cVol(1:dom.nPc);     % Cell volumes for physical cells

fCoord = dom.fCoord;            % Face coordinates
fXiLamba = dom.fXiLambda;       % Lambda parameter
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

dt = casedef.iteration.dt;
alpha = casedef.iteration.alpha;
maxNbTimeSteps = casedef.iteration.maxniter_stepping;

% Create an equation object for holding a scalar conservation equation
eqn_u = ScalarFvEqn2(dom);
eqn_v = ScalarFvEqn2(dom);
eqn_p = ScalarFvEqn2(dom);

time = -dt;
stepping = true;

while stepping % Loop for (false) time stepping
    time = time + dt;
    disp("Timestep: ")
    disp(time/dt)
    
    %%%%% PART I: MOMENTUM EQUATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Given the current iterate of the pressure field, compute the 
    % solution of the momentum equation for u* and v*
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    apu = zeros(nC,1); aNbIntu = zeros(2*nIf,1); aNbBoundu = zeros(2*nBf,1); 
    bu = zeros(nC,1);
    apv = zeros(nC,1); aNbIntv = zeros(2*nIf,1); aNbBoundv = zeros(2*nBf,1); 
    bv = zeros(nC,1);
    % Set all terms to zero
    reset(eqn_u); reset(eqn_v);
    
    % Coefficients for the momentum equation, given the pressure field
    for i = 1:nPc % Loop over physical cells
        Ui = U_old.data(1,i); Vi = U_old.data(2,i); % Current velocity fields
        apu(i) = apu(i) + cVol(i)/dt;
        apv(i) = apv(i) + cVol(i)/dt;
        bu(i) = bu(i) + cVol(i)*Ui/dt; % Pressure term is added while looping over the faces
        bv(i) = bv(i) + cVol(i)*Vi/dt; % Pressure term is added while looping over the faces
    end
       
    for i = 1:nIf % Loop over internal faces
       nb1 = fNbC(2*i-1); nb2 = fNbC(2*i);
       l = fXiLamba(i); Af = fArea(i); Xif = norm(Xi(:,i));
       U1 = U_old.data(:,nb1); U2 = U_old.data(:,nb2);
       Uf = l*U1 + (1-l)*U2;
       Unf = dot(Uf,fNormal(:,i));
          
       apu(nb1) = apu(nb1) + l*Unf*Af + k*Af/Xif; 
       apu(nb2) = apu(nb2) - l*Unf*Af + k*Af/Xif;
       apv(nb1) = apv(nb1) + l*Unf*Af + k*Af/Xif; 
       apv(nb2) = apv(nb2) - l*Unf*Af + k*Af/Xif;
       aNbIntu(2*i-1) = aNbIntu(2*i-1) + (1-l)*Unf*Af  - k*Af/Xif; % Upper diagonal (contributes to nb1 equation)
       aNbIntu(2*i) = aNbIntu(2*i) - (1-l)*Unf*Af - k*Af/Xif; % Lower diagonal (contributes to nb2 equation)
       aNbIntv(2*i-1) = aNbIntv(2*i-1) + (1-l)*Unf*Af  - k*Af/Xif; % Upper diagonal (contributes to nb1 equation)
       aNbIntv(2*i) = aNbIntv(2*i) - (1-l)*Unf*Af - k*Af/Xif; % Lower diagonal (contributes to nb2 equation)
          
       % Pressure terms: see Moukalled, p. 566
       % Note: if the face lies along x (y), the pressure contributes to bv
       % (bu). The first nIf/2 faces lie along y.
       Pf = l*P.data(nb1) + (1-l)*P.data(nb2);
       if i <= ceil(nIf/2) % Face along y
          dX = cCoord(1,nb2)-cCoord(1,nb1);
          bu(nb1) = bu(nb1) - cVol(nb1)*(Pf/dX)/rho;
          bu(nb2) = bu(nb2) + cVol(nb2)*(Pf/dX)/rho;
       else % Face along x
          dY = cCoord(2,nb2)-cCoord(2,nb1);
          bv(nb1) = bv(nb1) - cVol(nb1)*(Pf/dY)/rho;
          bv(nb2) = bv(nb2) + cVol(nb2)*(Pf/dY)/rho;
       end
    end

    for i = 1:nBf % Loop over boundary faces
       nb1 = fNbC(2*i-1 + 2*nIf); nb2 = fNbC(2*i + 2*nIf);
       l = fXiLamba(i + nIf); Af = fArea(i + nIf); Xif = norm(Xi(:,i + nIf));
       U1 = U_old.data(:,nb1); U2 = U_old.data(:,nb2);
       Uf = l*U1 + (1-l)*U2;
       Unf = dot(Uf,fNormal(:,i + nIf));
          
       % Pressure terms: see Moukalled p. 566
       % Note: if the face lies along x (y), the pressure contributes to bv
       % (bu). The first nBf/2 faces lie along y. Note that the pressure 
       % now only contributes to nb1 (the physical cell). 
       % Note that the - in dX or dY fixes the sign of dP.
       Pf = l*P.data(nb1) + (1-l)*P.data(nb2);
       if i <= ceil(nBf/2) % Face along y
          dX = cCoord(1,nb2)-cCoord(1,nb1);
          bu(nb1) = bu(nb1) - cVol(nb1)*(Pf/dX)/rho;
       else % Face along x
          dY = cCoord(2,nb2)-cCoord(2,nb1);
          bv(nb1) = bv(nb1) - cVol(nb1)*(Pf/dY)/rho; 
       end
          
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
                     apu(nb1) = apu(nb1) + l*Unf*Af + k*Af/Xif;
                     aNbBoundu(2*i-1) = aNbBoundu(2*i-1) + (1-l)*Unf*Af - k*Af/Xif;

                     % Forcing terms
                     phi_star = BC{j}.data.bcval_u;
                     bu(nb2) = phi_star;

                     % Diagonal elements from BC: ghost cell
                     apv(nb2) = l;
                     % Off-diagonal elements: physical cell
                     aNbBoundv(2*i) = 1-l;

                     % With inversed diffusion
                     apv(nb1) = apv(nb1) + l*Unf*Af + k*Af/Xif;
                     aNbBoundv(2*i-1) = aNbBoundv(2*i-1) + (1-l)*Unf*Af - k*Af/Xif;

                     % Forcing terms
                     phi_star = BC{j}.data.bcval_v; 
                     bv(nb2) = phi_star;
                 case 'Neumann'
                     % Diagonal elements from BC: ghost cell
                     apu(nb2) = -1/Xif;
                     % Off-diagonal elements: physical cell
                     aNbBoundu(2*i) = 1/Xif;

                     % With inversed diffusion
                     apu(nb1) = apu(nb1) + (1-l)*Unf*Af + k*Af/Xif;
                     aNbBoundu(2*i-1) = aNbBoundu(2*i-1) - k*Af/Xif + (1-l)*Unf*Af;

                     % Forcing terms
                     phi_star = BC{j}.data.bcval_u;
                     bu(nb2) = phi_star;  

                     % Diagonal elements from BC: ghost cell
                     apv(nb2) = -1/Xif;
                     % Off-diagonal elements: physical cell
                     aNbBoundv(2*i) = 1/Xif;

                     % With reversed diffusion
                     apv(nb1) = apv(nb1) + (1-l)*Unf*Af + k*Af/Xif;
                     aNbBoundv(2*i-1) = aNbBoundv(2*i-1) - k*Af/Xif + (1-l)*Unf*Af;

                     % Forcing terms
                     phi_star = BC{j}.data.bcval_v;
                     bv(nb2) = phi_star; 
             end
         end
       end
    end
       
    % Constitute matrices for the sparse solvers
    eqn_u.adata = [apu ; aNbIntu ; aNbBoundu]; eqn_u.bdata = bu;
    eqn_v.adata = [apv ; aNbIntv ; aNbBoundv]; eqn_v.bdata = bv; 
    [Au,bu] = to_msparse(eqn_u);
    [Av,bv] = to_msparse(eqn_v);
    x = get(U); 
    xu = x(1,:); xv = x(2,:);
    xu = xu'; xv = xv';
       
    % Check tolerance and iteration count (for the system solver)
    URes = bu-Au*xu; VRes = bv-Av*xv;
    UResnorm = norm(URes); VResnorm = norm(VRes);
    Resnorm = max(UResnorm,VResnorm);
    disp(Resnorm)

    % Update
    x = [Au\bu,Av\bv];
    set(U,x'); 
    
%     figure(1)
%     subplot(2,2,1)
%     title("U_{1,start}")
%     colorbar
%     fvmplotfield(U,'linear',1,1); 
%     subplot(2,2,2)
%     title("U_{2,start}")
%     colorbar
%     fvmplotfield(U,'linear',1,2); 
      
    %%%%% PART 2: PRESSURE EQUATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Given the current iterate of the velocities field (u*,v*), compute
    % the solution of the pressure equation for p'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    reset(P_prime,0); 
    
    app = zeros(nC,1); aNbIntp = zeros(2*nIf,1); aNbBoundp = zeros(2*nBf,1); 
    bp = zeros(nC,1);
    reset(eqn_p);
       
    % Coefficients for the pressure equation, given the pressure field     
    for i = 1:nIf % Loop over internal faces
       nb1 = fNbC(2*i-1); nb2 = fNbC(2*i);
       l = fXiLamba(i); Af = fArea(i);
       fVol = l*cVol(nb1)+(1-l)*cVol(nb2); % Volume of the shifted control cell
       apuNb1 = apu(nb1); apvNb1 = apv(nb1); 
       apuNb2 = apu(nb2); apvNb2 = apv(nb2); 
       apuf = l*apuNb1 + (1-l)*apuNb2;
       apvf = l*apvNb1 + (1-l)*apvNb2;
       P1 = P.data(:,nb1); P2 = P.data(:,nb2);
       Pf = l*P1 + (1-l)*P2;
       % Note: the first nIf/2 faces lie along y, watch out for the signs.
       if i <= ceil(nIf/2) % Face along y
          dX = cCoord(1,nb2)-cCoord(1,nb1);
          dPx = (P2-P1)/dX;
          dPx_ = Pf/dX;
          U1 = U.data(:,nb1); U2 = U.data(:,nb2);
          Uf_ = l*U1(1) + (1-l)*U2(1); % Linear interpolation
          Uf = Uf_ - fVol*(1/apuf)*(dPx-dPx_) ; %Rie-Chow interpolation
          %Unf = dot(Uf,fNormal(:,i));
          Unf = Uf;
          bp(nb1) = bp(nb1) - Unf*Af; % - because term brought to RHS
          bp(nb2) = bp(nb2) + Unf*Af;
           
          app(nb1) = app(nb1) + fVol*Af/(rho*dX*apuf); 
          app(nb2) = app(nb2) + fVol*Af/(rho*dX*apuf);
          aNbIntp(2*i-1) = aNbIntp(2*i-1) - fVol*Af/(rho*dX*apuf); % Upper diagonal
          aNbIntp(2*i) = aNbIntp(2*i) - fVol*Af/(rho*dX*apuf); % Lower diagonal
       else % Face along x
          dY = cCoord(2,nb2)-cCoord(2,nb1);
          dPy = (P2-P1)/dX;
          dPy_ = Pf/dX;
          U1 = U.data(:,nb1); U2 = U.data(:,nb2);
          Uf_ = l*U1(2) + (1-l)*U2(2); % Linear interpolation
          Uf = Uf_ - fVol*(1/apvf)*(dPy-dPy_) ; %Rie-Chow interpolation
          %Unf = dot(Uf,fNormal(:,i));
          Unf = Uf;
          bp(nb1) = bp(nb1) - Unf*Af; % - because term brought to RHS
          bp(nb2) = bp(nb2) + Unf*Af;
          app(nb1) = app(nb1) + fVol*Af/(rho*dY*apvf); 
          app(nb2) = app(nb2) + fVol*Af/(rho*dY*apvf); % sign flip
          aNbIntp(2*i-1) = aNbIntp(2*i-1) - fVol*Af/(rho*dY*apvf);
          aNbIntp(2*i) = aNbIntp(2*i) - fVol*Af/(rho*dY*apvf); % sign flip
       end
    end  

    for i = 1:nBf % Loop over boundary faces
       nb1 = fNbC(2*i-1 + 2*nIf); nb2 = fNbC(2*i + 2*nIf);
       l = fXiLamba(i + nIf); Af = fArea(i + nIf); Xif = norm(Xi(:,i + nIf));
       %U1 = U.data(:,nb1); U2 = U.data(:,nb2);
       %Uf = l*U1 + (1-l)*U2;
       %Unf = dot(Uf,fNormal(:,i + nIf));
       %fVol = l*cVol(nb1)+(1-l)*cVol(nb2); % cVol is not included for ghost cells, do not use fvol here!
       fVol = cVol(nb1); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       apuNb1 = apu(nb1); apvNb1 = apv(nb1);
          
       %bp(nb1) = bp(nb1) - Unf*Af; % - because term brought to RHS
       
       P1 = P.data(:,nb1); P2 = P.data(:,nb2);
       Pf = l*P1 + (1-l)*P2;
          
       boundaryFound = 0;
       for j = 1:nbZones
         % Determine to which boundary the face belongs
         if ismember(i + nIf,ranges(1,j):ranges(2,j)) && boundaryFound == 0
             boundaryFound = 1;
             switch BC{j}.kind_p
                 case 'Dirichlet' 
                     % Note: the first nBf/2 faces lie along y. Note that 
                     % dX and dY may become negative depending on the
                     % boundary.                       
                     if i <= ceil(nBf/2) % Face along y
                        dX = cCoord(1,nb2)-cCoord(1,nb1);
                        dPx = (P2-P1)/dX;
                        dPx_ = Pf/dX;
                        U1 = U.data(:,nb1); U2 = U.data(:,nb2);
                        Uf_ = l*U1(1) + (1-l)*U2(1); % Linear interpolation
                        Uf = Uf_ - fVol*(1/apuNb1)*(dPx-dPx_) ; %Rie-Chow interpolation
                        %Unf = dot(Uf,fNormal(:,i));
                        Unf = Uf;
                        bp(nb1) = bp(nb1) - Unf*Af; % - because term brought to RHS
                         
                        app(nb1) = app(nb1) + Af*Af/(rho*apuNb1); % sign flip
                        aNbBoundp(2*i-1) = aNbBoundp(2*i-1) - Af*Af/(rho*apuNb1); % sign flip
                     else % Face along x
                        dY = cCoord(2,nb2)-cCoord(2,nb1);
                        
                        dPy = (P2-P1)/dY;
                        dPy_ = Pf/dY;
                        U1 = U.data(:,nb1); U2 = U.data(:,nb2);
                        Uf_ = l*U1(2) + (1-l)*U2(2); % Linear interpolation
                        Uf = Uf_ - fVol*(1/apvNb1)*(dPy-dPy_) ; %Rie-Chow interpolation
                        %Unf = dot(Uf,fNormal(:,i));
                        Unf = Uf;
                        bp(nb1) = bp(nb1) - Unf*Af; % - because term brought to RHS
                        
                        app(nb1) = app(nb1) + Af*Af/(rho*apvNb1); 
                        aNbBoundp(2*i-1) = aNbBoundp(2*i-1) - Af*Af/(rho*apvNb1);
                     end
                     app(nb2) = l;
                     aNbBoundp(2*i) = 1-l; % Lower diagonal
                      
                     % Forcing terms
                     % Note: the boundary condition must be imposed for p',
                     % and p = p' + p*, so the value for the boundary
                     % condition is actually p' = p - p*
                     p_star = l*P.data(nb1) + (1-l)*P.data(nb2);
                     phi_star = BC{j}.data.bcval_p - p_star;
                     bp(nb2) = phi_star;
                      
                 case 'Neumann'
                     % Note: the first nBf/2 faces lie along y. Note that 
                     % dX and dY may become negative depending on the 
                     % boundary.
                     if i <= ceil(nBf/2) % Face along y
                        dX = abs(cCoord(1,nb2)-cCoord(1,nb1));
                        app(nb1) = app(nb1) + Af*Af/(rho*apuNb1); 
                        aNbBoundp(2*i-1) = aNbBoundp(2*i-1) - Af*Af/(rho*apuNb1);
                     else % Face along x
                        dY = abs(cCoord(2,nb2)-cCoord(2,nb1));
                        app(nb1) = app(nb1) + Af*Af/(rho*apvNb1); 
                        aNbBoundp(2*i-1) = aNbBoundp(2*i-1) - Af*Af/(rho*apvNb1);
                     end
                     app(nb2) = -1/Xif;
                     aNbBoundp(2*i) = 1/Xif;   
                      
                     % Forcing terms
                     phi_star = BC{j}.data.bcval_p - (P.data(nb1) - P.data(nb2))/Xif;
                     bp(nb2) = phi_star;                                              
             end
         end
       end
    end
       
    % Constitute matrices for the sparse solvers
    eqn_p.adata = [app ; aNbIntp ; aNbBoundp]; eqn_p.bdata = bp;
    [Ap,bp] = to_msparse(eqn_p);
    xp = get(P_prime); 
    xp = xp';
       
    % Check tolerance and iteration count (for the system solver)
    PRes = bp-Ap*xp;
    PResnorm = norm(PRes);   
    disp(PResnorm)

    % Update
    xp = Ap\bp;
    set(P_prime,xp'); % Put algebraic solution in the Field
    
    %%%%% PART 3: Corrections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Given the pressure correction, adjust the current pressure and 
    % velocities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure(2)
%     subplot(1,2,1)
%     title("P_{start}")
%     colorbar
%     fvmplotfield(P,'linear',1);  
    
    set(P,P.data + alpha*P_prime.data)
    Ucorr = 0*U.data;
    
%     figure(2)
%     subplot(1,2,2)
%     title("newP")
%     colorbar
%     fvmplotfield(P,'linear',1); 
    
    for i = 1:nIf
        nb1 = fNbC(2*i-1); nb2 = fNbC(2*i);
        l = fXiLamba(i); Af = fArea(i);
        %fVol = l*cVol(nb1)+(1-l)*cVol(nb2); % Volume of the shifted control cell
        apuNb1 = apu(nb1); apvNb1 = apv(nb1); 
        apuNb2 = apu(nb2); apvNb2 = apv(nb2); 
        %apuf = l*apuNb1 + (1-l)*apuNb2; % These are not used, see paper
        %apvf = l*apvNb1 + (1-l)*apvNb2;
        P1 = P_prime.data(:,nb1); P2 = P_prime.data(:,nb2);
        Pf = l*P1 + (1-l)*P2;
        % Note: the first nIf/2 faces lie along y.
        if i <= ceil(nIf/2) % Face along y
           Ucorr(1,nb1) = Ucorr(1,nb1) - Pf*Af/(rho*apuNb1);
           Ucorr(1,nb2) = Ucorr(1,nb2) + Pf*Af/(rho*apuNb2);
        else % Face along x
           %dY = cCoord(2,nb2)-cCoord(2,nb1);
           Ucorr(2,nb1) = Ucorr(2,nb1) - Pf*Af/(rho*apvNb1);
           Ucorr(2,nb2) = Ucorr(2,nb2) + Pf*Af/(rho*apvNb2);
        end 
    end
    
    for i = 1:nBf
        nb1 = fNbC(2*i-1 + 2*nIf); nb2 = fNbC(2*i + 2*nIf);
        l = fXiLamba(i + nIf); Af = fArea(i + nIf); % Xif = norm(Xi(:,i + nIf));
        %fVol = l*cVol(nb1)+(1-l)*cVol(nb2); % Volume of the shifted control cell
        apuNb1 = apu(nb1); apvNb1 = apv(nb1); 
        P1 = P_prime.data(:,nb1); P2 = P_prime.data(:,nb2);
        Pf = l*P1 + (1-l)*P2;
        % Note: the first nIf/2 faces lie along y.
        % See paper collocated, question what about Af???
        if i <= ceil(nBf/2) % Face along y
           dX = cCoord(1,nb2)-cCoord(1,nb1);           
           Ucorr(1,nb1) = Ucorr(1,nb1) - sign(dX)*Pf*Af/(rho*apuNb1);
           %Ucorr(1,nb2) = Ucorr(1,nb2) + Pf*Af/(rho*apuNb2);
        else % Face along x
           dY = cCoord(2,nb2)-cCoord(2,nb1);
           Ucorr(2,nb1) = Ucorr(2,nb1) - sign(dY)*Pf*Af/(rho*apvNb1);
           %Ucorr(2,nb2) = Ucorr(2,nb2) + Pf*Af/(rho*apvNb2);
        end 
        boundaryFound = 0;
        for j = 1:nbZones
            % Determine to which boundary the face belongs
            if ismember(i + nIf,ranges(1,j):ranges(2,j)) && boundaryFound == 0
                boundaryFound = 1;
                switch BC{j}.kind_u
                    case 'Dirichlet'
                        %phi_star_u = (BC{j}.data.bcval_u - l*Ucorr(1,nb1))/(1-l);
                        %phi_star_u = (- l*Ucorr(1,nb1))/(1-l);
                        %Ucorr(1,nb2) = phi_star_u;
                        Ucorr(1,nb2) = l*Ucorr(1,nb1)/(1-l);
                        %phi_star_v = (BC{j}.data.bcval_v - l*Ucorr(2,nb1))/(1-l);
                        %phi_star_v = (- l*Ucorr(2,nb1))/(1-l);
                        %Ucorr(2,nb2) = phi_star_v;       
                        Ucorr(2,nb2) = l*Ucorr(2,nb1)/(1-l);
                    case 'Neumann'
                        %phi_star_u = -(BC{j}.data.bcval_u*Xif - Ucorr(1,nb1));
                        %phi_star_u = Ucorr(1,nb1);
                        %Ucorr(1,nb2) = phi_star_u;
                        Ucorr(1,nb2) = Ucorr(1,nb1);
                        %phi_star_v = -(BC{j}.data.bcval_v*Xif - Ucorr(2,nb1));
                        %phi_star_v = Ucorr(2,nb1);
                        %Ucorr(2,nb2) = phi_star_v;
                        Ucorr(2,nb2) = Ucorr(2,nb1);
                end
            end
        end     
    end
    
    set(U,U.data + Ucorr);
    
%     figure(1)
%     subplot(2,2,3)
%     title("U_{1,corrected}")
%     colorbar
%     fvmplotfield(U,'linear',1,1); 
%     subplot(2,2,4)
%     title("U_{2,corrected}")
%     colorbar
%     fvmplotfield(U,'linear',1,2); 
     
    if stoppingCriterion(Ucorr(1,:),Ucorr(2,:),P_prime.data,casedef.iteration.tol) || (time/dt >= maxNbTimeSteps)
        stepping = false;
        if time/dt >= maxNbTimeSteps
            disp("maximum number of time steps reached")
        end
    else
        stepping = true;
        set(U_old,U.data);
    end
end % stepping

result.endtime = now; % call datestr(now) for displaying this time 
result.converged = 1-stepping;
%result.niter = niter;
result.steps = time/dt;
result.UResnorm = UResnorm;
result.URes = Field(dom.allCells,0);
   set(result.URes,URes');
result.VResnorm = VResnorm;
result.VRes = Field(dom.allCells,0);
   set(result.VRes,VRes');
result.U = U;
result.time = time;
result.P = P;
end



