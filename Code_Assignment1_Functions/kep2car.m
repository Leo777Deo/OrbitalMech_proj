function [r, v] = kep2car(a, e, i, OM, om, theta, mu)

% Trasformation from orbital (Keplerian) elements to Cartesian coordinates.
% ------------------------------------------------------------------------
% PROTOTYPE:
% [r, v] = kep2car (a, e, i, OM, om, theta, mu)
% ------------------------------------------------------------------------
% INPUT:
% a         [1x1]       semi-major axis                    [km]
% e         [1x1]       eccentricity                       [-]
% i         [1x1]       inclination                        [rad]
% OM        [1x1]       RAAN                               [rad]
% om        [1x1]       argument of pericentre             [rad]
% theta     [1x1]       true anomaly                       [rad]
% mu        [1x1]       gravitational parameter            [km^3/s^2]
% ------------------------------------------------------------------------
% OUTPUT:
% r        [3x1]       position vector                 [km]
% v        [3x1]       velocity vector                 [km/s]
% ------------------------------------------------------------------------
% CONTRIBUTORS:
% Andrea Barbiera 
% Leo De Luca
% ------------------------------------------------------------------------
%VERSIONS:
% 2023-12-01

% Directions of the perifocal reference frame
pp = [1; 0; 0]; % radial direction
qq = [0; 1; 0]; % transversal direction
ww = [0; 0; 1]; % out of plane direction

p = a * (1 - e^2);
h = sqrt(mu * p);

rr_pf = (h^2 / mu) * (1 / (1 + e * cos(theta))) * (cos(theta) * pp + sin(theta) * qq);

vv_pf = (mu / h) * (-sin(theta) * pp + (e + cos(theta)) * qq);

% Definition of the Direction Cosine Matrices
R_OM = [cos(OM) sin(OM) 0;
        -sin(OM) cos(OM) 0;
                     0 0 1];

R_i = [1 0 0;
       0 cos(i) sin(i);
       0 -sin(i) cos(i)];

R_w = [cos(om) sin(om) 0;
        -sin(om) cos(om) 0;
                     0 0 1];

% Passage from perfical r.f. (pf) to Earth centred equatorial r.f. (Ece)
A_Ece_pf = R_w * R_i * R_OM; 

A_pf_Ece = A_Ece_pf';
                                                 

r = A_pf_Ece * rr_pf;
v = A_pf_Ece * vv_pf;

end
