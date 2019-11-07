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
case_def = 'advection-diffusion'; % 'advection-diffusion-fixed-pressure', 'MMS'
mesh_def = 'standard-orthogonal';
init_velocity_field = 'zero-field';
MMS_case = 'example-slides';
boundary_cond = 'standard';

% Mesh parameters (are ignored if mesh_def is set)
nCx = 100; nCy = 100; % Number of cells
origin = [0, 0]; % Origin of the mesh
xAxis = [1, 0]; yAxis = [0, 1]; % x-axis and y-axis direction
skewX = 1.00; skewY = 1.00; % Mesh skewness

% Initial fields (are ignored if init_velocity_field is set)
uValue = 0;
vValue = 0;

% Pressure difference 
pIn = 50; pOut = 0;

% Material properties
% Define material properties
nu = 0.1; % Diffusion contant
k = 16; % Thermal conductivity [W/(m K)]
rho = 1; % Density [kg/m^3]

% MMS chosen solution (is ignored if MMS_case is set)
U_chosen = @(x,y) [x.^2;-y.^2];  % Imposed velocity field
Source = @(x,y) [4*x.^3 - 2*x.^2.*y - 2*nu ; -2*x.*y.^2 + 4*y.^3 + 2*nu]; % Analytically computed source term

% Boundary conditions (are ignored if boundary_cond is set)
    
%% Creation of the mesh
switch mesh_def
    case 'standard-orthogonal'
        nCx = 10; nCy = 10; 
        origin = [0, 0]; xAxis = [1, 0]; yAxis = [0, 1]; 
        skewX = 1.00; skewY = 1.00;
end
        
seedI = LineSeed.lineSeedOneWayBias(origin,xAxis,nCx,skewX,'o'); 
seedJ = LineSeed.lineSeedOneWayBias(origin,yAxis,nCy,skewY,'o'); 
casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
casedef.dom = newdomain(mesh,'MyDomain');
dX = seedI.displX; dY = seedJ.displY;
casedef.vars.dX = dX; casedef.vars.dY = dY; % Length of the domain 

%% Setting up initial fields
% Velocity field
U = Field(casedef.dom.allCells,1);     % Velocity [m/s] 
switch init_velocity_field
    case 'zero-field'
        set(U,[zeros(1,U.elcountzone);zeros(1,U.elcountzone)]);
    otherwise
        set(U,[uValueX*ones(1,U.elcountzone);uValueY*ones(1,U.elcountzone)]);     
end
casedef.vars.U = U;                    % Save U in the casedef
switch case_def
    case 'advection-diffusion' % Temperature field
        T = Field(casedef.dom.allCells,0);     % Temperature [K] 
        randomdata = rand(T.elsize,T.elcountzone)-0.5;
        set(T,randomdata);                     
    case 'MMS'
        U_sol = Field(casedef.dom.allCells,1);     % Velocity [m/s]
        U_diff = Field(casedef.dom.allCells,1);    % Error velocity [m/s] 
end

%% Define material properties
switch case_def
    case 'advection-diffusion'
        casedef.material.k = k; % Thermal conductivity [W/(m K)]
        casedef.material.rho = rho; % density [kg/m^3]
    case 'advection-diffusion-fixed-pressure'
        casedef.material.k = nu; % Thermal conductivity [W/(m K)]
        casedef.material.rho = rho; % density [kg/m^3]
    otherwise
        casedef.material.k = nu; % Thermal conductivity [W/(m K)]
        casedef.material.rho = rho; % density [kg/m^3]
end

%% Pressure definition
casedef.vars.pIn = pIn;
casedef.vars.pOut = pOut;

%% MMS source terms
switch MMS_case
    case 'example-slides'
        U_chosen = @(x,y) [x.^2;-y.^2];  % Imposed pressure field
        Source = @(x,y) [4*x.^3 - 2*x.^2.*y - 2*nu ; -2*x.*y.^2 + 4*y.^3 + 2*nu]; % Analytically computed source term
        casedef.vars.Source = Source;
    otherwise
        casedef.vars.Source = Source;
end

%% Boundary conditions
switch boundary_cond
    case 'standard' % Boundary conditions corresponding to the ones in the slides
        switch case_def
            case 'advection-diffusion'
                jBC = 0; jBC = jBC+1;
                casedef.BC{jBC}.zoneID = 'WESTRAND'; casedef.BC{jBC}.kind = 'Dirichlet'; casedef.BC{jBC}.data.bcval = 0;
                jBC = jBC+1;
                casedef.BC{jBC}.zoneID = 'OOSTRAND'; casedef.BC{jBC}.kind = 'Dirichlet'; casedef.BC{jBC}.data.bcval = 1;
                jBC = jBC+1;
                casedef.BC{jBC}.zoneID = 'ZUIDRAND'; casedef.BC{jBC}.kind = 'Neumann'; casedef.BC{jBC}.data.bcval = 0;
                jBC = jBC+1;
                casedef.BC{jBC}.zoneID = 'NOORDRAND'; casedef.BC{jBC}.kind = 'Neumann'; casedef.BC{jBC}.data.bcval = 0;
            case 'advection-diffusion-fixed-pressure'
        end
end




% Define boundary conditions: for u and v
switch MMS_case
   case 'example'
        jBC = 0;
        jBC = jBC+1;
        casedef.BC{jBC}.zoneID = 'WESTRAND';
        casedef.BC{jBC}.kind_u   = 'Dirichlet';
        casedef.BC{jBC}.data.bcval_u = U_chosen;
        casedef.BC{jBC}.kind_v   = 'Dirichlet';
        casedef.BC{jBC}.data.bcval_v = U_chosen;
        jBC = jBC+1;
        casedef.BC{jBC}.zoneID = 'OOSTRAND';
        casedef.BC{jBC}.kind_u   = 'Dirichlet';
        casedef.BC{jBC}.data.bcval_u = U_chosen;
        casedef.BC{jBC}.kind_v   = 'Dirichlet';
        casedef.BC{jBC}.data.bcval_v = U_chosen;
        jBC = jBC+1;
        casedef.BC{jBC}.zoneID = 'ZUIDRAND';
        casedef.BC{jBC}.kind_u   = 'Dirichlet';
        casedef.BC{jBC}.data.bcval_u = U_chosen;
        casedef.BC{jBC}.kind_v   = 'Dirichlet';
        casedef.BC{jBC}.data.bcval_v = U_chosen;
        jBC = jBC+1;
        casedef.BC{jBC}.zoneID = 'NOORDRAND';
        casedef.BC{jBC}.kind_u = 'Dirichlet';
        casedef.BC{jBC}.data.bcval_u = U_chosen;
        casedef.BC{jBC}.kind_v  = 'Dirichlet';
        casedef.BC{jBC}.data.bcval_v = U_chosen;
    case 'degree1'
        jBC = 0;
        jBC = jBC+1;
        casedef.BC{jBC}.zoneID = 'WESTRAND';
        casedef.BC{jBC}.kind_u   = 'Dirichlet';
        casedef.BC{jBC}.data.bcval_u = U_chosen;
        casedef.BC{jBC}.kind_v   = 'Dirichlet';
        casedef.BC{jBC}.data.bcval_v = U_chosen;
        jBC = jBC+1;
        casedef.BC{jBC}.zoneID = 'OOSTRAND';
        casedef.BC{jBC}.kind_u   = 'Dirichlet';
        casedef.BC{jBC}.data.bcval_u = U_chosen;
        casedef.BC{jBC}.kind_v   = 'Dirichlet';
        casedef.BC{jBC}.data.bcval_v = U_chosen;
        jBC = jBC+1;
        casedef.BC{jBC}.zoneID = 'ZUIDRAND';
        casedef.BC{jBC}.kind_u   = 'Dirichlet';
        casedef.BC{jBC}.data.bcval_u = U_chosen;
        casedef.BC{jBC}.kind_v   = 'Dirichlet';
        casedef.BC{jBC}.data.bcval_v = U_chosen;
        jBC = jBC+1;
        casedef.BC{jBC}.zoneID = 'NOORDRAND';
        casedef.BC{jBC}.kind_u = 'Dirichlet';
        casedef.BC{jBC}.data.bcval_u = U_chosen;
        casedef.BC{jBC}.kind_v  = 'Dirichlet';
        casedef.BC{jBC}.data.bcval_v = U_chosen;
end


% Set up iteration parameters
casedef.iteration.maxniter = 10; 
casedef.iteration.TTol     = 1e-6;
casedef.iteration.dt       = 50;

% Call solver
result = examplesolver5(casedef);
set(U_diff,abs(U_diff.data-result.U.data));

normal = Field(casedef.dom.allFaces,1);
tangent = Field(casedef.dom.allFaces,1);
xi = Field(casedef.dom.allFaces,1);
set(normal,[casedef.dom.fNormal]);
set(tangent,[casedef.dom.fTangent]);
set(xi,[casedef.dom.fXi]);

%% Plot result
%figure; hold on; axis off; axis equal; colormap(jet(50));
%colorbar
%scale = 'lin'; lw = 1;
%fvmplotfield(result.U,scale,lw,1); % 1 plots u, 2 plots v
%fvmplotmesh(casedef.dom,lw);
display("Number of false time steps: ")
display(result.steps)
display("End time: ")
display(result.time)
display("Max norm error: ")
display(result.UResnorm)
%Uoost = restrictto(U,getzone(casedef.dom,'OOSTRAND'));
%fvmplotvectorfield(Uoost,lw);
%fvmplotmesh(casedef.dom,lw);
%fvmplotvectorfield(normal,lw);
%fvmplotvectorfield(tangent,lw);
%fvmplotvectorfield(xi,lw);
%fvmplotcellnumbers(casedef.dom,8);
%fvmplotfacenumbers(casedef.dom,8);
%fvmplotvertexnumbers(casedef.dom,8);

%% Compare analtical
figure; hold on; axis off; axis equal; colormap(jet(50));
colorbar
scale = 'lin'; lw = 1;

subplot(2,3,1) % Computed solution X
colorbar
fvmplotfield(result.U,scale,lw,1);
%fvmplotmesh(casedef.dom,lw);
title('U_x code')

subplot(2,3,4) % Computed solution Y
colorbar
fvmplotfield(result.U,scale,lw,2);
%fvmplotmesh(casedef.dom,lw);
title('U_y code')

subplot(2,3,2) % Analytical solution X
colorbar
fvmplotfield(U_sol,scale,lw,1);
%fvmplotmesh(casedef.dom,lw);
title('U_x analytical')

subplot(2,3,5) % Analytical solution Y
colorbar
fvmplotfield(U_sol,scale,lw,2);
%fvmplotmesh(casedef.dom,lw);
title('U_y analytical')

subplot(2,3,3) % Error X
colorbar
fvmplotfield(U_diff,scale,lw,1);
%fvmplotmesh(casedef.dom,lw);
title('U_x error')

subplot(2,3,6) % Error Y
colorbar
fvmplotfield(U_diff,scale,lw,2);
%fvmplotmesh(casedef.dom,lw);
title('U_y error')











