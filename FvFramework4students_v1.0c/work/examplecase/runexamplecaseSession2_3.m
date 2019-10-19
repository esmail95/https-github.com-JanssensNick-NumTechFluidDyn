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


% Set up initial fields
T = Field(casedef.dom.allCells,0);     % Temperature [K] (scalar); empty field
% reset(T,0);                          % Reset with all zeros
randomdata = rand(T.elsize,T.elcountzone)-0.5;
set(T,randomdata);                     % Set with random numbers

uValueX = 0;
uValueY = 0;
U = Field(casedef.dom.allCells,1);     % Velocity [m/s] (vector);
%set(U,[rand(1,U.elcountzone);rand(1,U.elcountzone)]);
set(U,[uValueX*ones(1,U.elcountzone);uValueY*ones(1,U.elcountzone)]);
%reset(U,[1;0.2]);


% Define material properties
casedef.material.k = 16;  % Thermal conductivity [W/(m K)]
casedef.vars.U = U;

% Define boundary conditions
jBC = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'WESTRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 10;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'ZUIDRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = 0;

jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = 0;


% Set up iteration parameters
casedef.iteration.maxniter = 1000; 
casedef.iteration.TTol     = 1e-6;


% Call solver
result = examplesolver(casedef);

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
fvmplotfield(result.T,scale,lw);
%Uoost = restrictto(U,getzone(casedef.dom,'OOSTRAND'));
%fvmplotvectorfield(Uoost,lw);
fvmplotmesh(casedef.dom,lw);
%fvmplotvectorfield(normal,lw);
%fvmplotvectorfield(tangent,lw);
%fvmplotvectorfield(xi,lw);
%fvmplotcellnumbers(casedef.dom,8);
%fvmplotfacenumbers(casedef.dom,8);
%fvmplotvertexnumbers(casedef.dom,8);

%% Session 3: 
% North, South: Neumann, 0
% East: Dirichlet, 10
% West; Dirichlet, 0
figure
Uvalues = [0 25 50 75 100 125];
Pe = Lx*Uvalues*1000*4.18/(casedef.material.k);
ind0 = ceil(nCy/2);
step = nCy*(0:1:(nCx-1));
indices = ind0 + step;
leg = {};
cntr = 1;
for Ux = Uvalues
    set(U,[Ux*ones(1,U.elcountzone);0*ones(1,U.elcountzone)]);
    result = examplesolver(casedef);
    plotSection(casedef.dom.cCoord,result.T.data,indices,'x',2)
    hold on
    leg{cntr} = strcat('U_x = ', num2str(Ux), ', Pe = ',num2str(Pe(cntr)));
    cntr = cntr + 1;
end
lgd = legend(leg);
lgd.Location = 'northwest';
