%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Memristor FEM simulation - script only
% Reizinger Patrik (2018)
%
% Functionality:
% this function creates a simplified model for a memristor, specifies
% and solves a PDE for it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = memristor_pde(semi_a, semi_b, offset, plot_flag)
%% Model creation
model = createpde;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define geometry
% Each geometric object is a column vector in the geometry definition
% matrix, the following types are valid

% Circle: [1, x0, y0, r]'
% Polygon [2, #segments, x0,...,x(n-1), y0,...,y(n-1)]'
% Rectangle: [3, 4, x0, x1, x2, x3, y0, y1, y2, y3]'
% Ellipse: [4, x0, y0, semi_a, semi_b, angle_in_rad]

%%%%%%%%%%%%%%%%%%%%
% Bounding rectangle
% A - along x axis
% B - along y axis
rect_x = 4e-08;
rect_y = 2.5e-08;
A = 2*rect_x;
B = 2*rect_y;
R1 = [3, 4, -rect_x, -rect_x, rect_x, rect_x, -rect_y, rect_y, rect_y, -rect_y]';

%%%%%%%%%%%%%%%%%%%%
% Simulated dendrite-like part
% A - along x axis
% B - along y axis
% semi_a = 20e-9;
% semi_b = 10e-9;
E1 = [4, -rect_x, 0, semi_a, semi_b, 0]';
E1 = [E1; zeros(length(R1)-length(E1),1)];  % append zeros for dimension matching

gm = [R1, E1];  % geometry description matrix
sf = 'R1-E1';   % set formula
ns = char('R1', 'E1');  %namespace matrix
ns = ns';

%%%%%%%%%%%%%%%%%%%%
% Create the geometry
g = decsg(gm, sf, ns);

% Include the gemoetry in the model and plot it
geometryFromEdges(model, g);

%%%%%%%%%%%%%%%%%%%%
% plot the geometry
% notice the G in the command name, it's pdeGplot
pdegplot(model, 'EdgeLabels', 'on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Boundary conditions

%%%%%%%%%%%%%%%%%%%%
% Dirichlet
% the scheme is:
% h*u = r

% homogeneous for edges 4-7
applyBoundaryCondition(model, 'dirichlet', 'Edge', 4:7, 'u', 0);

% apply the excitation on the other side - edge 2
u_excitation = 0.3;
applyBoundaryCondition(model, 'dirichlet', 'Edge', 2, 'h', 1, 'r', u_excitation);

%%%%%%%%%%%%%%%%%%%%
% Neumann (homogeneous)
% the scheme is:
% n*c*grad(u) + q*u = g
applyBoundaryCondition(model, 'neumann', 'Edge', [1, 3], 'g', 0, 'q', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

generateMesh(model, 'Hmax', 3e-9); % default: mean jiggle, 10 max jiggle iterations, preR2013a mehing algorithm
pdeplot(model); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PDE coefficient specification
% the scheme is:
% -div(c*grad(u)) + a*u = f

eps0 = 8.85e-12;
eps_r = 6;
% offset = 3.5e-8;
omega_area = A*B - pi/2*semi_a*semi_b;
n = 1.2e27; %[m^-3] average ionic concentration

c = @(region, state) eps0*eps_r*(1-1./(1+exp(18e7*(region.x-offset))));

specifyCoefficients(model, 'm', 0,...
                           'd', 0,...
                           'c', c,...
                           'a', 0,...
                           'f', n*omega_area*1.6e-19);
                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PDE solving

results = solvepde(model);

% export the structure elements into variables with more specific names
u = results.NodalSolution;
ux = results.XGradients;
uy = results.YGradients;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualization
if plot_flag
    pdeplot(model, 'XYData', u, ... specify what to plot on the XY-plane with colors
                   'ZData', u, ... specify value showed on the Z-axis
                   'FaceAlpha', 0.5, ... set the opacity of the face of the plot
                   'FlowData', [ux, uy], ... include a Quiver plot
                   'ColorMap', 'jet', ... colormap
                   'Contour', 'on' ... plot contours
                   );
end

end