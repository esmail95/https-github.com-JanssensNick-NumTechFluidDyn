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
Nx=100;Lx=1;
Ny=100;Ly=1;
casedef.vars.Nx=Nx;
casedef.vars.Ny=Ny;
casedef.vars.Lx=Lx;
casedef.vars.Ly=Ly;
% Create a mesh
seedI = LineSeed.lineSeedOneWayBias([0 0],[1 0],Nx,Lx,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 1],Ny,Ly,'o');
casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
% Create domain from mesh
casedef.dom = newdomain(mesh,'MyDomain');
normal = Field(casedef.dom.allFaces,1);
tangent = Field(casedef.dom.allFaces,1);
xi = Field(casedef.dom.allFaces,1);
set(normal,[casedef.dom.fNormal]);
set(tangent,[casedef.dom.fTangent]);
set(xi,[casedef.dom.fXi]);


% Set up initial fields
T = Field(casedef.dom.allCells,0);     % Temperature [K] (scalar); empty field
% reset(T,0);                          % Reset with all zeros
randomdata = rand(T.elsize,T.elcountzone)-0.5;
set(T,randomdata);                     % Set with random numbers

Up_s =Field(casedef.dom.allCells,0);      
randomdata = rand(Up_s.elsize,Up_s.elcountzone)-0.5;
set(Up_s,randomdata);

Vp_s = Field(casedef.dom.allCells,0);
randomdata = rand(Up_s.elsize,Up_s.elcountzone)-0.5;
set(Vp_s,randomdata);

uValueX = 50;
uValueY = 0;
U = Field(casedef.dom.allCells,1);     % Velocity [m/s] (vector);
%set(U,[rand(1,U.elcountzone);rand(1,U.elcountzone)]);
set(U,[0*ones(1,U.elcountzone);uValueY*ones(1,U.elcountzone)]);
%reset(U,[1;0.2]);
% dx=seedI.displ.x;
% dy=seedJ.displ.y;

P_in=100;
P_out=0;
dt=20;
% Define material properties
casedef.material.k = 0.1;  %16 % Thermal conductivity [W/(m K)]
casedef.material.rho = 1;
casedef.vars.U = U;
casedef.vars.P_in=P_in;
casedef.vars.P_out=P_out;


% Define boundary conditions
jBC = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'WESTRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'ZUIDRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 0;


% Define boundary conditions_u
jBC = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'WESTRAND';
casedef.BC{jBC}.kind_u   = 'Neumann';
casedef.BC{jBC}.data.bcval_u = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.kind_u   = 'Neumann';
casedef.BC{jBC}.data.bcval_u = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'ZUIDRAND';
casedef.BC{jBC}.kind_u   = 'Dirichlet';
casedef.BC{jBC}.data.bcval_u = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.kind_u   = 'Dirichlet';
casedef.BC{jBC}.data.bcval_u = 0;

% Define boundary conditions_v
jBC = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'WESTRAND';
casedef.BC{jBC}.kind_v   = 'Neumann';
casedef.BC{jBC}.data.bcval_v = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.kind_v   = 'Neumann';
casedef.BC{jBC}.data.bcval_v = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'ZUIDRAND';
casedef.BC{jBC}.kind_v   = 'Dirichlet';
casedef.BC{jBC}.data.bcval_v = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.kind_v   = 'Dirichlet';
casedef.BC{jBC}.data.bcval_v = 0;


% Set up iteration parameters
casedef.iteration.maxniter = 1000;
casedef.iteration.TTol     = 1e-6;
casedef.iteration.dt       = dt;


% Call solver
%result = examplesolver(casedef);
result_u_v=examplesolver_u_v(casedef);



% Plot result
figure; hold on; axis off; axis equal; colormap(jet(50));
colorbar
scale = 'lin'; lw = 1;
%fvmplotfield(T,scale,lw);
fvmplotfield(result_u_v.U,scale,lw,1);
% % Uoost = restrictto(U,getzone(casedef.dom,'OOSTRAND'));
% % fvmplotvectorfield(Uoost,lw);
fvmplotmesh(casedef.dom,lw);
%fvmplotcellnumbers(casedef.dom,8);
%fvmplotvectorfield(normal,lw);
%fvmplotvectorfield(tangent,lw);
%fvmplotvectorfield(xi,lw);
%fvmplotfacenumbers(casedef.dom,8);
% fvmplotvertexnumbers(casedef.dom,8);

%Plot cross section species
% t = result.T.data;
% 
% xmesh = linspace(0,Lx,Nx);
% t_plot = t(1:Ny:casedef.dom.nPc);
% 
% figure
% plot(xmesh,t_plot,'linewidth',2)
% colorbar;
% set(gca,'Fontsize',12);xlabel('x','Fontsize',14);ylabel('\phi','FontSize',14)
% title('Cross section of species value')



