%==========================================================================
%
% Case setup script for the FVMLab framework 4 students
%
% Purpose: Creating a mesh, defining materials, defining
%          boundary conditions, defining iteration parameters, and finally
%          calling the solver for different applications: diffusion,
%          advection-diffusion, advection-diffusion with fixed pressure
%          gradient, verification with method of manufactured solutions...
%
% by Esmail Mostatira, Nick Janssens (source code: Frederik Rogiers)
%
%==========================================================================

clear variables
close all
clc

%% User parameters 

% Imlemented cases are: 
%   - advection-diffusion
%   - advection-diffusion with fixed pressure gradient
%   - verification with method of manufactured solutions
%   - Pressuree-velocity with odd-even decoupling
case_def = 'rie-chow'; % 'advection-diffusion', 'advection-diffusion-fixed-pressure', 'MMS', 'odd-even-decoupling', 'rie-chow'
MMS_case = 'example-slides';
boundary_cond = 'standard'; %'standard'
p_init = 'checkerboard'; %'checkerboard'
analytical = 1;

% Mesh parameters (are ignored if mesh_def is set)
nCx = 50; nCy = 50; % Number of cells
origin = [0, 0]; % Origin of the mesh
xAxis = [1, 0]; yAxis = [0, 1]; % x-axis and y-axis direction
skewX = 1.00; skewY = 1.00; % Mesh skewness

% Initial fields (are ignored if init_velocity_field is set)
uValue = 0;
vValue = 0;

% Pressure difference 
pIn = 50; pOut = 0;
pValue = 0; % Initialise constant pressure field

% Material properties
% Define material properties
nu = 0.1; % Diffusion contant
k = 16; % Thermal conductivity [W/(m K)]
rho = 1; % Density [kg/m^3]

% MMS chosen solution (is ignored if MMS_case is set)
U_chosen = @(x,y) [x.^2;-y.^2];  % Imposed velocity field
Source = @(x,y) [4*x.^3 - 2*x.^2.*y - 2*nu ; -2*x.*y.^2 + 4*y.^3 + 2*nu]; % Analytically computed source term

% Boundary conditions (are ignored if boundary_cond is set)
% Order: west, east, south, north
uBC_names = {'Neumann','Neumann','Dirichlet','Dirichlet'};
vBC_names = {'Neumann','Neumann','Dirichlet','Dirichlet'};
TBC_names = {'Dirichlet','Dirichlet','Neumann','Neumann'};
pBC_names = {'Dirichlet','Dirichlet','Neumann','Neumann'};
uBC_values = [0, 0, 0, 0];
vBC_values = [0, 0, 0, 0];
TBC_values = [0, 1, 0, 0];
pBC_values = [pIn, pOut, 0, 0];

% Iteration parameters
maxniter_solver = 1000;
maxniter_stepping = 5000;
tol      = 1e-6;
dt       = 50;
alpha    = 0.1; % Relaxation, the more cells, the lower this value should be

%% Case setup
% Creation of the mesh   
seedI = LineSeed.lineSeedOneWayBias(origin,xAxis,nCx,skewX,'o'); 
seedJ = LineSeed.lineSeedOneWayBias(origin,yAxis,nCy,skewY,'o'); 
casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
casedef.dom = newdomain(mesh,'MyDomain');
dX = seedI.displX; dY = seedJ.displY;
casedef.vars.dX = dX; casedef.vars.dY = dY; % Length of the domain 

% Setting up initial fields
% Velocity field
U = Field(casedef.dom.allCells,1);     % Velocity [m/s] 
set(U,[uValue*ones(1,U.elcountzone);vValue*ones(1,U.elcountzone)]);     
casedef.vars.U = U;                    % Save U in the casedef
% Temperature field
T = Field(casedef.dom.allCells,0);     % Temperature [K] 
randomdata = rand(T.elsize,T.elcountzone)-0.5;
set(T,randomdata);                     
% MMS solution and error fields
U_sol = Field(casedef.dom.allCells,1);     % Velocity [m/s]
U_diff = Field(casedef.dom.allCells,1);    % Error velocity [m/s] 
% Pressure field
P = Field(casedef.dom.allCells,0);     % Pressure [Pa] 
set(P,pValue*ones(1,P.elcountzone));   
% Constant pressure gradient
% pVals = [100*ones(1,10) 90*ones(1,10) 80*ones(1,10) 70*ones(1,10) 60*ones(1,10) 50*ones(1,10) 40*ones(1,10) 30*ones(1,10) 20*ones(1,10) 10*ones(1,10) 110*ones(1,10) 0*ones(1,10) 100 90 80 70 60 50 40 30 20 10 100 90 80 70 60 50 40 30 20 10];
% Checkerboard:
if strcmp(p_init,'checkerboard')
    pVals = 0*P.data;
    for k = 1:(nCx*nCy/2)
        pVals(2*k) = pValue;
    end
    for k = 1:(nCx+nCy)
        pVals(2*k-1) = pValue;
    end
    set(P,pVals);
end
casedef.vars.P = P;                    % Save P in the casedef

% Define material properties
switch case_def
    case 'advection-diffusion'
        casedef.material.k = k; % Thermal conductivity [W/(m K)]
        casedef.material.rho = rho; % density [kg/m^3]
    case 'advection-diffusion-fixed-pressure' 
        casedef.material.k = nu; % Thermal conductivity [W/(m K)]
        casedef.material.rho = rho; % density [kg/m^3]
    case 'odd-even-decoupling'
        casedef.material.k = nu; % Thermal conductivity [W/(m K)]
        casedef.material.rho = rho; % density [kg/m^3]    
    case 'rie-chow'
        casedef.material.k = nu; % Thermal conductivity [W/(m K)]
        casedef.material.rho = rho; % density [kg/m^3]     
    otherwise
        casedef.material.k = nu; % Thermal conductivity [W/(m K)]
        casedef.material.rho = rho; % density [kg/m^3]
end

% Pressure definition
casedef.vars.pIn = pIn;
casedef.vars.pOut = pOut;

% MMS source terms
switch MMS_case
    case 'example-slides'
        U_chosen = @(x,y) [x.^2;-y.^2];  % Imposed pressure field
        set(U_sol,U_chosen(casedef.dom.cCoord(1,:),casedef.dom.cCoord(2,:)));
        set(U_diff,U_chosen(casedef.dom.cCoord(1,:),casedef.dom.cCoord(2,:))); 
        Source = @(x,y) [4*x.^3 - 2*x.^2.*y - 2*nu ; -2*x.*y.^2 + 4*y.^3 + 2*nu]; % Analytically computed source term
        casedef.vars.Source = Source;
    otherwise
        set(U_sol,U_chosen(casedef.dom.cCoord(1,:),casedef.dom.cCoord(2,:)));
        set(U_diff,U_chosen(casedef.dom.cCoord(1,:),casedef.dom.cCoord(2,:))); 
        casedef.vars.Source = Source;
end

% Boundary conditions
switch boundary_cond
    case 'standard' % Boundary conditions corresponding to the ones in the slides
        switch case_def
            case 'advection-diffusion'
                jBC = 0; jBC = jBC+1; casedef.BC{jBC}.zoneID = 'WESTRAND'; 
                casedef.BC{jBC}.kind = 'Dirichlet'; casedef.BC{jBC}.data.bcval = 0;
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'OOSTRAND'; 
                casedef.BC{jBC}.kind = 'Dirichlet'; casedef.BC{jBC}.data.bcval = 1;
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'ZUIDRAND'; 
                casedef.BC{jBC}.kind = 'Neumann'; casedef.BC{jBC}.data.bcval = 0;
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'NOORDRAND'; 
                casedef.BC{jBC}.kind = 'Neumann'; casedef.BC{jBC}.data.bcval = 0;
            case 'advection-diffusion-fixed-pressure'
                jBC = 0; jBC = jBC+1; casedef.BC{jBC}.zoneID = 'WESTRAND'; 
                casedef.BC{jBC}.kind_u = 'Neumann'; casedef.BC{jBC}.data.bcval_u = 0;
                casedef.BC{jBC}.kind_v = 'Neumann'; casedef.BC{jBC}.data.bcval_v = 0;
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'OOSTRAND';
                casedef.BC{jBC}.kind_u = 'Neumann'; casedef.BC{jBC}.data.bcval_u = 0;
                casedef.BC{jBC}.kind_v = 'Neumann'; casedef.BC{jBC}.data.bcval_v = 0;
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'ZUIDRAND';
                casedef.BC{jBC}.kind_u = 'Dirichlet'; casedef.BC{jBC}.data.bcval_u = 0;
                casedef.BC{jBC}.kind_v = 'Dirichlet'; casedef.BC{jBC}.data.bcval_v = 0;
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'NOORDRAND';
                casedef.BC{jBC}.kind_u = 'Dirichlet'; casedef.BC{jBC}.data.bcval_u = 0;
                casedef.BC{jBC}.kind_v = 'Dirichlet'; casedef.BC{jBC}.data.bcval_v = 0;
            case 'MMS'
                jBC = 0; jBC = jBC+1; casedef.BC{jBC}.zoneID = 'WESTRAND';
                casedef.BC{jBC}.kind_u = 'Dirichlet'; casedef.BC{jBC}.data.bcval_u = U_chosen;
                casedef.BC{jBC}.kind_v = 'Dirichlet'; casedef.BC{jBC}.data.bcval_v = U_chosen;
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'OOSTRAND';
                casedef.BC{jBC}.kind_u = 'Dirichlet'; casedef.BC{jBC}.data.bcval_u = U_chosen;
                casedef.BC{jBC}.kind_v = 'Dirichlet'; casedef.BC{jBC}.data.bcval_v = U_chosen;
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'ZUIDRAND';
                casedef.BC{jBC}.kind_u = 'Dirichlet'; casedef.BC{jBC}.data.bcval_u = U_chosen;
                casedef.BC{jBC}.kind_v = 'Dirichlet'; casedef.BC{jBC}.data.bcval_v = U_chosen;
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'NOORDRAND';
                casedef.BC{jBC}.kind_u = 'Dirichlet'; casedef.BC{jBC}.data.bcval_u = U_chosen;
                casedef.BC{jBC}.kind_v = 'Dirichlet'; casedef.BC{jBC}.data.bcval_v = U_chosen; 
            case 'odd-even-decoupling'
                jBC = 0; jBC = jBC+1; casedef.BC{jBC}.zoneID = 'WESTRAND'; 
                casedef.BC{jBC}.kind_u = 'Neumann'; casedef.BC{jBC}.data.bcval_u = 0;
                casedef.BC{jBC}.kind_v = 'Neumann'; casedef.BC{jBC}.data.bcval_v = 0;
                casedef.BC{jBC}.kind_p = 'Dirichlet'; casedef.BC{jBC}.data.bcval_p = pIn;
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'OOSTRAND';
                casedef.BC{jBC}.kind_u = 'Neumann'; casedef.BC{jBC}.data.bcval_u = 0;
                casedef.BC{jBC}.kind_v = 'Neumann'; casedef.BC{jBC}.data.bcval_v = 0;
                casedef.BC{jBC}.kind_p = 'Dirichlet'; casedef.BC{jBC}.data.bcval_p = pOut;
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'ZUIDRAND';
                casedef.BC{jBC}.kind_u = 'Dirichlet'; casedef.BC{jBC}.data.bcval_u = 0;
                casedef.BC{jBC}.kind_v = 'Dirichlet'; casedef.BC{jBC}.data.bcval_v = 0;
                casedef.BC{jBC}.kind_p = 'Neumann'; casedef.BC{jBC}.data.bcval_p = 0;
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'NOORDRAND';
                casedef.BC{jBC}.kind_u = 'Dirichlet'; casedef.BC{jBC}.data.bcval_u = 0;
                casedef.BC{jBC}.kind_v = 'Dirichlet'; casedef.BC{jBC}.data.bcval_v = 0;
                casedef.BC{jBC}.kind_p = 'Neumann'; casedef.BC{jBC}.data.bcval_p = 0;
            case 'rie-chow'
                jBC = 0; jBC = jBC+1; casedef.BC{jBC}.zoneID = 'WESTRAND'; 
                casedef.BC{jBC}.kind_u = 'Neumann'; casedef.BC{jBC}.data.bcval_u = 0;
                casedef.BC{jBC}.kind_v = 'Neumann'; casedef.BC{jBC}.data.bcval_v = 0;
                casedef.BC{jBC}.kind_p = 'Dirichlet'; casedef.BC{jBC}.data.bcval_p = pIn;
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'OOSTRAND';
                casedef.BC{jBC}.kind_u = 'Neumann'; casedef.BC{jBC}.data.bcval_u = 0;
                casedef.BC{jBC}.kind_v = 'Neumann'; casedef.BC{jBC}.data.bcval_v = 0;
                casedef.BC{jBC}.kind_p = 'Dirichlet'; casedef.BC{jBC}.data.bcval_p = pOut;
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'ZUIDRAND';
                casedef.BC{jBC}.kind_u = 'Dirichlet'; casedef.BC{jBC}.data.bcval_u = 0;
                casedef.BC{jBC}.kind_v = 'Dirichlet'; casedef.BC{jBC}.data.bcval_v = 0;
                casedef.BC{jBC}.kind_p = 'Neumann'; casedef.BC{jBC}.data.bcval_p = 0;
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'NOORDRAND';
                casedef.BC{jBC}.kind_u = 'Dirichlet'; casedef.BC{jBC}.data.bcval_u = 0;
                casedef.BC{jBC}.kind_v = 'Dirichlet'; casedef.BC{jBC}.data.bcval_v = 0;
                casedef.BC{jBC}.kind_p = 'Neumann'; casedef.BC{jBC}.data.bcval_p = 0;
        end
    otherwise
        switch case_def
            case 'advection-diffusion'
                jBC = 0; jBC = jBC+1; casedef.BC{jBC}.zoneID = 'WESTRAND'; 
                casedef.BC{jBC}.kind = TBC_names(jBC); casedef.BC{jBC}.data.bcval = TBC_values(jBC);
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'OOSTRAND'; 
                casedef.BC{jBC}.kind = TBC_names(jBC); casedef.BC{jBC}.data.bcval = TBC_values(jBC);
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'ZUIDRAND'; 
                casedef.BC{jBC}.kind = TBC_names(jBC); casedef.BC{jBC}.data.bcval = TBC_values(jBC);
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'NOORDRAND'; 
                casedef.BC{jBC}.kind = TBC_names(jBC); casedef.BC{jBC}.data.bcval = TBC_values(jBC);
            case 'advection-diffusion-fixed-pressure'
                jBC = 0; jBC = jBC+1; casedef.BC{jBC}.zoneID = 'WESTRAND'; 
                casedef.BC{jBC}.kind_u = uBC_names(jBC); casedef.BC{jBC}.data.bcval_u = uBC_values(jBC);
                casedef.BC{jBC}.kind_v = vBC_names(jBC); casedef.BC{jBC}.data.bcval_v = vBC_values(jBC);
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'OOSTRAND';
                casedef.BC{jBC}.kind_u = uBC_names(jBC); casedef.BC{jBC}.data.bcval_u = uBC_values(jBC);
                casedef.BC{jBC}.kind_v = vBC_names(jBC); casedef.BC{jBC}.data.bcval_v = vBC_values(jBC);
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'ZUIDRAND';
                casedef.BC{jBC}.kind_u = uBC_names(jBC); casedef.BC{jBC}.data.bcval_u = uBC_values(jBC);
                casedef.BC{jBC}.kind_v = vBC_names(jBC); casedef.BC{jBC}.data.bcval_v = vBC_values(jBC);
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'NOORDRAND';
                casedef.BC{jBC}.kind_u = uBC_names(jBC); casedef.BC{jBC}.data.bcval_u = uBC_values(jBC);
                casedef.BC{jBC}.kind_v = vBC_names(jBC); casedef.BC{jBC}.data.bcval_v = vBC_values(jBC);
            case 'MMS'
                jBC = 0; jBC = jBC+1; casedef.BC{jBC}.zoneID = 'WESTRAND';
                casedef.BC{jBC}.kind_u = 'Dirichlet'; casedef.BC{jBC}.data.bcval_u = U_chosen;
                casedef.BC{jBC}.kind_v = 'Dirichlet'; casedef.BC{jBC}.data.bcval_v = U_chosen;
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'OOSTRAND';
                casedef.BC{jBC}.kind_u = 'Dirichlet'; casedef.BC{jBC}.data.bcval_u = U_chosen;
                casedef.BC{jBC}.kind_v = 'Dirichlet'; casedef.BC{jBC}.data.bcval_v = U_chosen;
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'ZUIDRAND';
                casedef.BC{jBC}.kind_u = 'Dirichlet'; casedef.BC{jBC}.data.bcval_u = U_chosen;
                casedef.BC{jBC}.kind_v = 'Dirichlet'; casedef.BC{jBC}.data.bcval_v = U_chosen;
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'NOORDRAND';
                casedef.BC{jBC}.kind_u = 'Dirichlet'; casedef.BC{jBC}.data.bcval_u = U_chosen;
                casedef.BC{jBC}.kind_v = 'Dirichlet'; casedef.BC{jBC}.data.bcval_v = U_chosen; 
            case 'odd-even-decoupling'
                jBC = 0; jBC = jBC+1; casedef.BC{jBC}.zoneID = 'WESTRAND'; 
                casedef.BC{jBC}.kind_u = uBC_names(jBC); casedef.BC{jBC}.data.bcval_u = uBC_values(jBC);
                casedef.BC{jBC}.kind_v = vBC_names(jBC); casedef.BC{jBC}.data.bcval_v = vBC_values(jBC);
                casedef.BC{jBC}.kind_p = pBC_names(jBC); casedef.BC{jBC}.data.bcval_p = pBC_values(jBC);
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'OOSTRAND';
                casedef.BC{jBC}.kind_u = uBC_names(jBC); casedef.BC{jBC}.data.bcval_u = uBC_values(jBC);
                casedef.BC{jBC}.kind_v = vBC_names(jBC); casedef.BC{jBC}.data.bcval_v = vBC_values(jBC);
                casedef.BC{jBC}.kind_p = pBC_names(jBC); casedef.BC{jBC}.data.bcval_p = pBC_values(jBC);
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'ZUIDRAND';
                casedef.BC{jBC}.kind_u = uBC_names(jBC); casedef.BC{jBC}.data.bcval_u = uBC_values(jBC);
                casedef.BC{jBC}.kind_v = vBC_names(jBC); casedef.BC{jBC}.data.bcval_v = vBC_values(jBC);
                casedef.BC{jBC}.kind_p = pBC_names(jBC); casedef.BC{jBC}.data.bcval_p = pBC_values(jBC);
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'NOORDRAND';
                casedef.BC{jBC}.kind_u = uBC_names(jBC); casedef.BC{jBC}.data.bcval_u = uBC_values(jBC);
                casedef.BC{jBC}.kind_v = vBC_names(jBC); casedef.BC{jBC}.data.bcval_v = vBC_values(jBC);
                casedef.BC{jBC}.kind_p = pBC_names(jBC); casedef.BC{jBC}.data.bcval_p = pBC_values(jBC);
            case 'rie-chow'
                jBC = 0; jBC = jBC+1; casedef.BC{jBC}.zoneID = 'WESTRAND'; 
                casedef.BC{jBC}.kind_u = uBC_names(jBC); casedef.BC{jBC}.data.bcval_u = uBC_values(jBC);
                casedef.BC{jBC}.kind_v = vBC_names(jBC); casedef.BC{jBC}.data.bcval_v = vBC_values(jBC);
                casedef.BC{jBC}.kind_p = pBC_names(jBC); casedef.BC{jBC}.data.bcval_p = pBC_values(jBC);
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'OOSTRAND';
                casedef.BC{jBC}.kind_u = uBC_names(jBC); casedef.BC{jBC}.data.bcval_u = uBC_values(jBC);
                casedef.BC{jBC}.kind_v = vBC_names(jBC); casedef.BC{jBC}.data.bcval_v = vBC_values(jBC);
                casedef.BC{jBC}.kind_p = pBC_names(jBC); casedef.BC{jBC}.data.bcval_p = pBC_values(jBC);
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'ZUIDRAND';
                casedef.BC{jBC}.kind_u = uBC_names(jBC); casedef.BC{jBC}.data.bcval_u = uBC_values(jBC);
                casedef.BC{jBC}.kind_v = vBC_names(jBC); casedef.BC{jBC}.data.bcval_v = vBC_values(jBC);
                casedef.BC{jBC}.kind_p = pBC_names(jBC); casedef.BC{jBC}.data.bcval_p = pBC_values(jBC);
                jBC = jBC+1; casedef.BC{jBC}.zoneID = 'NOORDRAND';
                casedef.BC{jBC}.kind_u = uBC_names(jBC); casedef.BC{jBC}.data.bcval_u = uBC_values(jBC);
                casedef.BC{jBC}.kind_v = vBC_names(jBC); casedef.BC{jBC}.data.bcval_v = vBC_values(jBC);
                casedef.BC{jBC}.kind_p = pBC_names(jBC); casedef.BC{jBC}.data.bcval_p = pBC_values(jBC);
        end
end

% Set up iteration parameters
casedef.iteration.maxniter_solver   = maxniter_solver; 
casedef.iteration.maxniter_stepping = maxniter_stepping; 
casedef.iteration.tol               = tol;
casedef.iteration.dt                = dt;
casedef.iteration.alpha             = alpha;

%% Call solver
switch case_def
    case 'advection-diffusion'
        result = examplesolver(casedef);
        disp("Max norm error: ")
        display(result.TResnorm)
    case 'advection-diffusion-fixed-pressure'
        result = examplesolver4(casedef);
        disp("Number of false time steps: ")
        display(result.steps)
        disp("End time: ")
        display(result.time)
        disp("Max norm error: ")
        display(result.UResnorm)
    case 'MMS'
        result = examplesolver5(casedef);
        set(U_diff,abs(U_diff.data-result.U.data));
        disp("Number of false time steps: ")
        display(result.steps)
        disp("End time: ")
        display(result.time)
        disp("Max norm error: ")
        display(result.UResnorm)
    case 'odd-even-decoupling'
        result = examplesolver6bis(casedef);
        disp("Number of false time steps: ")
        display(result.steps)
        disp("End time: ")
        display(result.time)
        disp("Max norm error: ")
        display(result.UResnorm) 
    case 'rie-chow'
        result = examplesolver7(casedef);
        disp("Number of false time steps: ")
        display(result.steps)
        disp("End time: ")
        display(result.time)
        disp("Max norm error: ")
        display(result.UResnorm)
end

normal = Field(casedef.dom.allFaces,1);
tangent = Field(casedef.dom.allFaces,1);
xi = Field(casedef.dom.allFaces,1);
set(normal,[casedef.dom.fNormal]);
set(tangent,[casedef.dom.fTangent]);
set(xi,[casedef.dom.fXi]);

%% Plot results: Solution plot
figure; hold on; axis off; axis equal; colormap(jet(50)); colorbar; scale = 'lin'; lw = 1;
switch case_def
    case 'advection-diffusion'
        title("Temperature distribution")
        fvmplotfield(result.T,scale,lw);
    case 'advection-diffusion-fixed-pressure'
        subplot(1,2,1); colorbar
        title("U_x")
        fvmplotfield(result.U,scale,lw,1); 
        subplot(1,2,2); colorbar
        title("U_y")
        fvmplotfield(result.U,scale,lw,2); 
    case 'MMS'
        subplot(1,2,1); colorbar
        title("U_x")
        fvmplotfield(result.U,scale,lw,1); 
        subplot(1,2,2); colorbar
        title("U_y")
        fvmplotfield(result.U,scale,lw,2);
    case 'odd-even-decoupling'
        subplot(1,3,1); colorbar
        title("U_x")
        fvmplotfield(result.U,scale,lw,1); 
        subplot(1,3,2); colorbar
        title("U_y")
        fvmplotfield(result.U,scale,lw,2);
        subplot(1,3,3); colorbar
        title("P")
        fvmplotfield(result.P,scale,lw);
    case 'rie-chow'
        subplot(1,3,1); colorbar
        title("U_x")
        fvmplotfield(result.U,scale,lw,1); 
        subplot(1,3,2); colorbar
        title("U_y")
        fvmplotfield(result.U,scale,lw,2);
        subplot(1,3,3); colorbar
        title("P")
        fvmplotfield(result.P,scale,lw);
end
fvmplotmesh(casedef.dom,lw);
%fvmplotvectorfield(normal,lw);
%fvmplotvectorfield(tangent,lw);
%fvmplotvectorfield(xi,lw);
%fvmplotcellnumbers(casedef.dom,8);
%fvmplotfacenumbers(casedef.dom,8);
%fvmplotvertexnumbers(casedef.dom,8);

%% Plot results: MMS comparison plots
if strcmp(case_def,'MMS')
    figure; hold on; axis off; axis equal; colormap(jet(50)); colorbar; scale = 'lin'; lw = 1;

    subplot(2,3,1); colorbar; % Computed solution X
    fvmplotfield(result.U,scale,lw,1);
    %fvmplotmesh(casedef.dom,lw);
    title('U_x code')

    subplot(2,3,4); colorbar; % Computed solution Y
    fvmplotfield(result.U,scale,lw,2);
    %fvmplotmesh(casedef.dom,lw);
    title('U_y code')

    subplot(2,3,2); colorbar; % Analytical solution X
    fvmplotfield(U_sol,scale,lw,1);
    %fvmplotmesh(casedef.dom,lw);
    title('U_x analytical')

    subplot(2,3,5); colorbar; % Analytical solution Y
    fvmplotfield(U_sol,scale,lw,2);
    %fvmplotmesh(casedef.dom,lw);
    title('U_y analytical')

    subplot(2,3,3); colorbar; % Error X
    fvmplotfield(U_diff,scale,lw,1);
    %fvmplotmesh(casedef.dom,lw);
    title('U_x error')

    subplot(2,3,6); colorbar; % Error Y
    fvmplotfield(U_diff,scale,lw,2);
    %fvmplotmesh(casedef.dom,lw);
    title('U_y error')
end

%% Plot results: plot section
switch case_def
    case 'advection-diffusion'
        figure
        Lx = xAxis(1);
        Uvalues = uValue;
        Pe = Lx*Uvalues*1/(casedef.material.k);
        ind0 = ceil(nCy/2); step = nCy*(0:1:(nCx-1)); indices = ind0 + step;
        leg = {}; leg2 = {};
        [XX,T] = plotSection(casedef.dom.cCoord,result.T.data,indices,'x',2);
        leg{1} = strcat('U_x = ', num2str(Uvalues), ', Pe = ',num2str(Pe));
        lgd = legend(leg); lgd.Location = 'northwest';
        grid on  
        xlabel("x")
        ylabel("T")
    case 'advection-diffusion-fixed-pressure'
        figure
        pValues = pIn-pOut;
        ind0 = nCy*ceil(nCx/2)+1;
        indices = ind0:(ind0+nCy-1);
        leg = {};
        plotSection(casedef.dom.cCoord,result.U.data(1,:),indices,'y',2);
        leg{1} = strcat('p_{in} = ', num2str(pValues));
        lgd = legend(leg);
        lgd.Location = 'northwest';
        grid on;
        xlabel("x")
        ylabel("U_x")
    case 'odd-even-decoupling'
        figure
        subplot(1,2,1)
        ind0 = nCy*ceil(nCx/2)+1;
        indices = ind0:(ind0+nCy-1);
        leg = {};
        plotSection(casedef.dom.cCoord,result.U.data(1,:),indices,'y',2);
        leg{1} = strcat('p_{in} = ', num2str(pIn), ', p_{out} = ', num2str(pOut));
        lgd = legend(leg);
        lgd.Location = 'northwest';
        grid on;
        xlabel("x")
        ylabel("U_x")
        subplot(1,2,2)
        ind0 = ceil(nCy/2);
        indices = ind0 + nCy*(0:(nCx-1));
        leg = {};
        plotSection(casedef.dom.cCoord,result.P.data,indices,'x',2);
        leg{1} = strcat('p_{in} = ', num2str(pIn), 'p_{out} = ', num2str(pOut));
        lgd = legend(leg);
        lgd.Location = 'northwest';
        grid on;
        xlabel("x")
        ylabel("p")
    case 'rie-chow'
        figure
        subplot(1,2,1)
        ind0 = nCy*ceil(nCx/2)+1;
        indices = ind0:(ind0+nCy-1);
        leg = {};
        plotSection(casedef.dom.cCoord,result.U.data(1,:),indices,'y',2);
        leg{1} = strcat('p_{in} = ', num2str(pIn), ', p_{out} = ', num2str(pOut));
        lgd = legend(leg);
        lgd.Location = 'northwest';
        grid on;
        xlabel("x")
        ylabel("U_x")
        subplot(1,2,2)
        ind0 = ceil(nCy/2);
        indices = ind0 + nCy*(0:(nCx-1));
        leg = {};
        plotSection(casedef.dom.cCoord,result.P.data,indices,'x',2);
        leg{1} = strcat('p_{in} = ', num2str(pIn), 'p_{out} = ', num2str(pOut));
        lgd = legend(leg);
        lgd.Location = 'northwest';
        grid on;
        xlabel("x")
        ylabel("p")
end
        
%% Plot results: error on analytical soltion
if analytical
    switch case_def
        case 'advection-diffusion'
            if ~strcmp(boundary_cond,'standard')
                error("Boundary conditions must be set to standard for analytical solution")
            else
                T_anal = @(x,Pe) (exp(Pe*x) - 1 ) ./ (exp(Pe)-1);
                Lx = xAxis(1); Ly = yAxis(2);
                x0 = Lx/nCx/2; dx = Lx/nCx; xL = x0 + (nCx-1)*dx;
                x = x0:dx:xL;

                figure
                Uvalues = [20 60 100];
                Pe = Lx*Uvalues*1/(casedef.material.k);
                ind0 = ceil(nCy/2); step = nCy*(0:1:(nCx-1)); indices = ind0 + step;
                leg = {}; leg2 = {};
                cntr = 1;
                for Ux = Uvalues
                    set(U,[Ux*ones(1,U.elcountzone);0*ones(1,U.elcountzone)]);
                    casedef.vars.U = U;                    
                    result = examplesolver(casedef);
                    subplot(1,2,1); hold on;
                    [XX,T] = plotSection(casedef.dom.cCoord,result.T.data,indices,'x',2);
                    leg{2*cntr-1} = strcat('U_x = ', num2str(Ux), ', Pe = ',num2str(Pe(cntr)));
                    T_anal_x = T_anal(x,Pe(cntr));
                    plot(x,T_anal_x,'--','LineWidth',1)
                    leg{2*cntr} = strcat('Analytic: U_x = ', num2str(Ux), ', Pe = ',num2str(Pe(cntr)));
                    subplot(1,2,2); hold on;
                    stem(x,T-T_anal_x)
                    leg2{cntr} = strcat('Analytic: U_x = ', num2str(Ux), ', Pe = ',num2str(Pe(cntr)));
                    cntr = cntr + 1;
                end
                subplot(1,2,1)
                lgd = legend(leg); lgd.Location = 'northwest';
                grid on
                subplot(1,2,2)
                xlabel('x'); ylabel('error');
                lgd = legend(leg2); lgd.Location = 'southwest';
                grid on
            end
        case 'odd-even-decoupling'
            if ~strcmp(boundary_cond,'standard')
                error("Boundary conditions must be set to standard for analytical solution")
            else
                Lx = xAxis(1); Ly = yAxis(2);
                y0 = Ly/nCy/2; dy = Ly/nCy; yL = y0 + (nCy-1)*dy;
                y = y0:dy:yL;
                U_anal = @(y,dPx) -1/2*(Ly/2)^2/(nu*rho)*dPx*(1-(y/(Ly/2)).^2); %http://www.nar-associates.com/nar-publishing/lfachap2.pdf

                figure
                title(strcat('p_{in} = ', num2str(pIn), 'p_{out} = ', num2str(pOut)))
                subplot(1,2,1)
                hold on
                ind0 = nCy*ceil(nCx/2)+1;
                indices = ind0:(ind0+nCy-1);
                leg = {};
                [U,XX] = plotSection(casedef.dom.cCoord,result.U.data(1,:),indices,'y',2);
                U_anal_y = U_anal(y-Ly/2,-pIn);
                plot(U_anal_y,y,'--','LineWidth',1)
                leg{1} = strcat('p_{in} = ', num2str(pIn), 'p_{out} = ', num2str(pOut));
                leg{2} = "Analytic";
                lgd = legend(leg);
                lgd.Location = 'northeast';
                grid on;
                xlabel("U_x")
                ylabel("y")
                
                subplot(1,2,2)
                stem(y,U-U_anal_y)
                leg2 = strcat('error');
                lgd = legend(leg2);
                grid on
                xlabel("y")
                ylabel("error")
            end
        case 'rie-chow'
            if ~strcmp(boundary_cond,'standard')
                error("Boundary conditions must be set to standard for analytical solution")
            else
                Lx = xAxis(1); Ly = yAxis(2);
                y0 = Ly/nCy/2; dy = Ly/nCy; yL = y0 + (nCy-1)*dy;
                y = y0:dy:yL;
                U_anal = @(y,dPx) -1/2*(Ly/2)^2/(nu*rho)*dPx*(1-(y/(Ly/2)).^2); %http://www.nar-associates.com/nar-publishing/lfachap2.pdf

                figure
                title(strcat('p_{in} = ', num2str(pIn), 'p_{out} = ', num2str(pOut)))
                subplot(1,2,1)
                hold on
                ind0 = nCy*ceil(nCx/2)+1;
                indices = ind0:(ind0+nCy-1);
                leg = {};
                [U,XX] = plotSection(casedef.dom.cCoord,result.U.data(1,:),indices,'y',2);
                U_anal_y = U_anal(y-Ly/2,-pIn);
                plot(U_anal_y,y,'--','LineWidth',1)
                leg{1} = strcat('p_{in} = ', num2str(pIn), 'p_{out} = ', num2str(pOut));
                leg{2} = "Analytic";
                lgd = legend(leg);
                lgd.Location = 'northeast';
                grid on;
                xlabel("U_x")
                ylabel("y")
                
                subplot(1,2,2)
                stem(y,U-U_anal_y)
                leg2 = strcat('error');
                lgd = legend(leg2);
                grid on
                xlabel("y")
                ylabel("error")
            end
    end
end











