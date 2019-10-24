%==========================================================================
%
% Example case using the FVMLab framework 4 students
%
% Purpose: Provides an example for setting up a case and calling a solver.
%          This involves creating a mesh, defining materials, defining
%          boundary conditions, defining iteration parameters, and finally
%          calling the solver.
%
% by Frederik Rogiers
%
%==========================================================================

% TIP: use "clear variables" instead of "clear all" to clear variables
%      use "clear classes" when the interface of a class has changes
%      use "close all" to close figures
%      use "clc" to clear the command window
% 
% TIP: pressing CTRL+D while the cursor is on a function opens that function
%      in the m-editor. This is the most convenient way of browsing through
%      your source code.

clear variables
clc


% Create a mesh
% Orthogonal mesh
nCx = 10;
nCy = 10;
Lx = 1;
Ly = 1;
seedI = LineSeed.lineSeedOneWayBias([0 0],[Lx 0],nCx,1.00,'o'); 
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 Ly],nCy,1.00,'o'); 
casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
% Create domain from mesh
casedef.dom = newdomain(mesh,'MyDomain');

dX = seedI.displX;
dY = seedJ.displY;
casedef.vars.dX = dX;
casedef.vars.dY = dY;

% Set up initial fields
uValueX = 0;
uValueY = 0;
U = Field(casedef.dom.allCells,1);     % Velocity [m/s] (vector);
set(U,[uValueX*ones(1,U.elcountzone);uValueY*ones(1,U.elcountzone)]);
%reset(U,[1;0.2]);

% Define material properties
k = 0.1;
rho = 1;
casedef.material.k = k; % Thermal conductivity [W/(m K)]
casedef.material.rho = rho; % density [kg/m^3]

casedef.vars.U = U;
casedef.vars.pIn = 100;
casedef.vars.pOut = 0;

% MMS solution
U_sol = Field(casedef.dom.allCells,1);     % Velocity [m/s] (vector);
U_diff = Field(casedef.dom.allCells,1);
MMS_case = 'example';
switch MMS_case
   case 'example'
      U_chosen = @(x,y) [x.^2;-y.^2];  % Imposed pressure field
      set(U_sol,U_chosen(casedef.dom.cCoord(1,:),casedef.dom.cCoord(2,:)));
      set(U_diff,U_chosen(casedef.dom.cCoord(1,:),casedef.dom.cCoord(2,:))); 
      % Derivation in slides is wrong, k must be included
      Source =@(x,y) [4*x.^3 - 2*x.^2.*y - 2*k ; -2*x.*y.^2 + 4*y.^3 + 2*k]; % Analytically computed source term
      casedef.vars.Source = Source;
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
end


% Set up iteration parameters
casedef.iteration.maxniter = 10; 
casedef.iteration.TTol     = 1e-6;
casedef.iteration.dt       = 50;

% Call solver
result = examplesolver5(casedef);
set(U_diff,U_diff.data-result.U.data);

normal = Field(casedef.dom.allFaces,1);
tangent = Field(casedef.dom.allFaces,1);
xi = Field(casedef.dom.allFaces,1);
set(normal,[casedef.dom.fNormal]);
set(tangent,[casedef.dom.fTangent]);
set(xi,[casedef.dom.fXi]);

%% Plot result
% figure; hold on; axis off; axis equal; colormap(jet(50));
% colorbar
% scale = 'lin'; lw = 1;
% fvmplotfield(result.U,scale,lw,1); % 1 plots u, 2 plots v
% fvmplotmesh(casedef.dom,lw);
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
fvmplotmesh(casedef.dom,lw);
title('U_x code')

subplot(2,3,4) % Computed solution Y
colorbar
fvmplotfield(result.U,scale,lw,2);
fvmplotmesh(casedef.dom,lw);
title('U_y code')

subplot(2,3,2) % Analytical solution X
colorbar
fvmplotfield(U_sol,scale,lw,1);
fvmplotmesh(casedef.dom,lw);
title('U_x analytical')

subplot(2,3,5) % Analytical solution Y
colorbar
fvmplotfield(U_sol,scale,lw,2);
fvmplotmesh(casedef.dom,lw);
title('U_y analytical')

subplot(2,3,3) % Error X
colorbar
fvmplotfield(U_diff,scale,lw,1);
fvmplotmesh(casedef.dom,lw);
title('U_x error')

subplot(2,3,6) % Error Y
colorbar
fvmplotfield(U_diff,scale,lw,2);
fvmplotmesh(casedef.dom,lw);
title('U_y error')











