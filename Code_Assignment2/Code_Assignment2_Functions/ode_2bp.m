function dy = ode_2bp(~, y)

%ode_2bp ODE system for the two-body problem (Keplerian motion)
%----------------------------------------------------------------------------------------
% PROTOTYPE
% dy = ode_2bp( t, y)
%----------------------------------------------------------------------------------------
% INPUT:
% t          [1]      Time (can be omitted, as the system is autonomous)    [s]
% y          [6x1]    State of the body (rx, ry, rz, vx, vy, vz)            [km][km/s]
%----------------------------------------------------------------------------------------
% OUTPUT:
% dy         [6x1]    Derivative of the state                               [km/s][km/s^2]
%----------------------------------------------------------------------------------------
% CONTRIBUTORS:
% Viola Poverini
% Leo De Luca 
%----------------------------------------------------------------------------------------
% VERSIONS
% 2023-10-05: First version
%----------------------------------------------------------------------------------------

% Parameters 
mu_E = astroConstants(13);

% Position and velocity
r = y(1:3);
v = y(4:6);

% Distance from the primary
rnorm = norm(r);

% Set the derivatives of the state
dy = [v; (-mu_E/rnorm^3)*r];

end 