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
nCx = 100;
nCy = 100;
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
casedef.material.k = 0.1; % Thermal conductivity [W/(m K)]
casedef.material.rho = 1; % density [kg/m^3]

casedef.vars.U = U;
casedef.vars.pIn = 100;
casedef.vars.pOut = 0;

% Define boundary conditions: for u and v
jBC = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'WESTRAND';
casedef.BC{jBC}.kind_u   = 'Neumann';
casedef.BC{jBC}.data.bcval_u = 0;
casedef.BC{jBC}.kind_v   = 'Neumann';
casedef.BC{jBC}.data.bcval_v = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.kind_u   = 'Neumann';
casedef.BC{jBC}.data.bcval_u = 0;
casedef.BC{jBC}.kind_v   = 'Neumann';
casedef.BC{jBC}.data.bcval_v = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'ZUIDRAND';
casedef.BC{jBC}.kind_u   = 'Dirichlet';
casedef.BC{jBC}.data.bcval_u = 0;
casedef.BC{jBC}.kind_v   = 'Dirichlet';
casedef.BC{jBC}.data.bcval_v = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.kind_u = 'Dirichlet';
casedef.BC{jBC}.data.bcval_u = 0;
casedef.BC{jBC}.kind_v  = 'Dirichlet';
casedef.BC{jBC}.data.bcval_v = 0;

% Set up iteration parameters
casedef.iteration.maxniter = 10; 
casedef.iteration.TTol     = 1e-6;
casedef.iteration.dt       = 50;%50;

% Call solver
result = examplesolver4(casedef);

normal = Field(casedef.dom.allFaces,1);
tangent = Field(casedef.dom.allFaces,1);
xi = Field(casedef.dom.allFaces,1);
set(normal,[casedef.dom.fNormal]);
set(tangent,[casedef.dom.fTangent]);
set(xi,[casedef.dom.fXi]);

%% Plot result
figure; hold on; axis off; axis equal; colormap(jet(50));
colorbar
scale = 'lin'; lw = 1;
fvmplotfield(result.U,scale,lw,1); % 1 plots u, 2 plots v
fvmplotmesh(casedef.dom,lw);
display("Number of false time steps: ")
display(result.steps)
display("End time: ")
display(result.time)
display("Max norm error: ")
display(result.UResnorm)
%%%
%figure; hold on; axis off; axis equal; colormap(jet(50));
%fvmplotvectorfield(result.U,lw,1);
%fvmplotmesh(casedef.dom,lw);
%%%
%Uoost = restrictto(U,getzone(casedef.dom,'OOSTRAND'));
%fvmplotvectorfield(Uoost,lw);
%fvmplotmesh(casedef.dom,lw);
%fvmplotvectorfield(normal,lw);
%fvmplotvectorfield(tangent,lw);
%fvmplotvectorfield(xi,lw);
%fvmplotcellnumbers(casedef.dom,8);
%fvmplotfacenumbers(casedef.dom,8);
%fvmplotvertexnumbers(casedef.dom,8);

%% Session 4: 
% North, South: Dirichlet, 0
% East, West: Neumann, 0
figure
pValues = [0 20 40 60 80 100];
ind0 = nCy*ceil(nCx/2)+1;
indices = ind0:(ind0+nCy-1);
leg = {};
cntr = 1;
for px = pValues
    casedef.vars.pIn = px;
    result = examplesolver4(casedef);
    plotSection(casedef.dom.cCoord,result.U.data(1,:),indices,'y',2)
    hold on
    leg{cntr} = strcat('p_{in} = ', num2str(px));
    cntr = cntr + 1;
end
lgd = legend(leg);
lgd.Location = 'northwest';

