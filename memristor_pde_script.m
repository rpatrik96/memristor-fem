%% Memristor FEM simulation - script only
% Reizinger Patrik
% W5PDBR

%% Model creation
model = createpde;

%% Define geometry
% Each geometric object is a column vector in the geometry definition
% matrix, the following types are valid

% Circle: [1, x0, y0, r]'
% Polygon [2, #segments, x0,...,x(n-1), y0,...,y(n-1)]'
% Rectangle: [3, 4, x0, x1, x2, x3, y0, y1, y2, y3]'
% Ellipse: [4, x0, y0, semi_a, semi_b, angle_in_rad]

% Bounding rectangle
rect_x = 4e-08;
rect_y = 2.5e-08;
R1 = [3, 4, -rect_x, -rect_x, rect_x, rect_x, -rect_y, rect_y, rect_y, -rect_y]';

% Simulated dendrite-like part
% A - along x axis
% B - along y axis
semi_a = 20e-9;
semi_b = 10e-9;
E1 = [4, -rect_x, 0, semi_a, semi_b, 0]';
E1 = [E1; zeros(length(R1)-length(E1),1)];  % append zeros for dimension matching

gm = [R1, E1];  % geometry description matrix
sf = 'R1-E1';   % set formula
ns = char('R1', 'E1');  %namespace matrix
ns = ns';

% Create the geometry
g = decsg(gm, sf, ns);

% Include the gemoetry in the model and plot it
geometryFromEdges(model, g);

% plot the geometry
% notice the G in the command name, it's pdeGplot
pdegplot(model, 'EdgeLabels', 'on');




