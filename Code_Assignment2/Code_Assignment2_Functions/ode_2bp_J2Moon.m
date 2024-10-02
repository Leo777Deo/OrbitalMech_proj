function dy = ode_2bp_J2Moon(t, y, fun_a3B)

% Evaluation of the acceleration due to J2 and Moon perturbations in ECI frame 
% ---------------------------------------------------------------------------------------
% PROTOTYPE:
% dy = ode_2bp_J2Moon(t, y, fun_a3B)
% ---------------------------------------------------------------------------------------
% INPUT:
% t               [1x1]       Time                                             [s]
% y               [6x1]       State of the body (rx, ry, rz, vx, vy, vz)       [km][km/s]
% fun_a3B:        Function to evaluate the acceleration due to the Moon in ECI
% ---------------------------------------------------------------------------------------
% OUTPUT:
% dy              [6x1]       Derivative of the state                          [km/s][km/s^2]
% ---------------------------------------------------------------------------------------
% CONTRIBUTORS:
% Viola Poverini
% Andrea Barbiera
% ---------------------------------------------------------------------------------------
% VERSIONS:
% 2023-12-13: First Version
% 2024-01-01: Last version 
% ---------------------------------------------------------------------------------------

% Parameters 
mu_E = astroConstants(13); %Earth 
R_E = astroConstants(23);
mu_M = astroConstants(20); %Moon 

% Position and velocity of the S/C
r = y(1:3);
v = y(4:6);

% Distance of the S/C from the primary
rnorm = norm(r);

% J2 effect 
J2 = astroConstants(9);
aj2i = (3/2)*((J2*mu_E*R_E^2)/rnorm^4)*(y(1)/rnorm)*((5*y(3)^2/rnorm^2)-1);
aj2j = (3/2)*((J2*mu_E*R_E^2)/rnorm^4)*(y(2)/rnorm)*((5*y(3)^2/rnorm^2)-1);
aj2k = (3/2)*((J2*mu_E*R_E^2)/rnorm^4)*(y(3)/rnorm)*((5*y(3)^2/rnorm^2)-3);
aj2 = [aj2i aj2j aj2k];
aj2norm = norm(aj2);

% Perturbation of the Moon
modr = rnorm; 
[a3B] = fun_a3B(t,y);

% Set the derivatives of the state
dy = [v; (-mu_E/rnorm^3)*r + aj2' + a3B'];

end 

