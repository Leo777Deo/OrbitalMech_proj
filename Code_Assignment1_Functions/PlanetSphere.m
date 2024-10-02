function PlanetSphere(R,X,Y,Z,C)
% 
% PLANETSPHERE   The function gives the Planet's texture loaded in a plot.
% 
% PlanetSphere(R,X,Y,Z,C)  
%  
% Input arguments:
% ------------------------------------------------------------------------
%   R          [1x1]       Planet mean radius                          [km]
%   X          [1x1]       x coordinate of the center of the planet wrt the
%                          origin of the reference frame               [km]
%   Y          [1x1]       y coordinate of the center of the planet wrt the
%                          origin of the reference frame               [km]
%   Z          [1x1]       z coordinate of the center of the planet wrt the
%                          origin of the reference frame               [km]
%   C          [1x1 struct]It cointains the picture of the Planet
% Output arguments:
% ------------------------------------------------------------------------
%   []          [figure]    Figure open with the Planet picture loaded
%
% AUTHOR:
%   Andrea Barbiera
%
%   Leo De Luca
%   
%   Gianluca Perusini
%
%   Viola Poverini
%
%  CHANGELOG:
%   25/11/23, Andrea Barbiera: Cartesian axes 2 times as long as the radius
%                              of the planet

% Link to the Planet image:
Planet_Image = C.foto;   

% Background colour:
background_plot = 'w';

hold on;
grid on;
axis equal;

% Initial view
view(120,30);

%% Create the Earth surface as a wireframe:

% Number of panels necessary to model the sphere:
N_panels = 180; 

% 3D meshgrid of the sphere points through the 'ellipsoid' function:
[x, y, z] = ellipsoid(X, Y, Z, R, R, R, N_panels);

% Create the globe with the 'surf' function
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 'none','handlevisibility','off');


%% Texturemap the globe:

% Load the Earth image for texturemap:
cdata = imread(Planet_Image);


% Transparency of the globe: 
% if alpha = 1: opaque; if alpha = 0: invisible
alpha = 1;

% Set the 'FaceColor' to 'texturemap' to apply an image on the globe, and
% specify the image data using the 'CData' property with the data loaded 
% from the image. Finally, set the transparency and remove the edges.
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');

hold on
grid on

% I,J,K axis in cartesian state:
quiver3(0,0,0,2*R,0,0,'k','linewidth', 1, 'handlevisibility', 'off')
quiver3(0,0,0,0,2*R,0,'k', 'linewidth', 1, 'handlevisibility', 'off')
quiver3(0,0,0,0,0,2*R,'k', 'linewidth', 1 , 'handlevisibility', 'off')

end