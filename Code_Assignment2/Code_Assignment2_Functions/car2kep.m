function [a, e, i, OM, om, th] = car2kep (rr, vv, mu)

% Transformation from Cartesian state to orbital elements 
% ------------------------------------------------------------------------
% PROTOTYPE:
% [a, e , i, OM, om, theta] = car2kep (rr, vv, mu)
% ------------------------------------------------------------------------
% INPUT:
% rr        [3x1]       Position vector                        [km]
% vv        [3x1]       Velocity vector                        [km/s]
% mu        [1x1]       Gravitational parameter of the primary [km^3/s^2]
% ------------------------------------------------------------------------
% OUTPUT:
% a         [1x1]       Semi-major axis                        [km]
% e         [1x1]       Eccentricity                           [-]
% i         [1x1]       Inclination                            [rad]
% OM        [1x1]       RAAN                                   [rad]
% om        [1x1]       Pericenter anomaly                     [rad]
% theta     [1x1]       True anomaly                           [rad]
% ------------------------------------------------------------------------
% CONTRIBUTORS:
% Andrea Barbiera
% Gianluca Perusini
% ------------------------------------------------------------------------
% VERSIONS:
% 2023-11-25: First version
% 2023-12-18: Last version

% Directions:
% I = [1; 0; 0];
% J = [0; 1; 0];
K = [0; 0; 1];

% Distance:
r = norm(rr);

% Velocity:
v = norm(vv);

% Radial velocity;
v_r = dot(rr, vv) / r;

% Semi-major axis:
a = - mu * r / (v^2 * r - 2 * mu);

% Angular momentum:
hh = cross(rr, vv);
h = norm(hh);

% Inclination:
i1 = acos(hh(3) / h);

if i1<1e-6
    i=0;
else
    i=i1;
end

% Line of nodes:
NN = cross(K, hh);
N = norm(NN);

% RAAN;

if i == 0
    OM = 0;
else
        if NN(2) >= 0
            OM = acos(NN(1)/ N);
        else
            OM = 2*pi - acos(NN(1) / N);
        end
end



% Eccentricity;
ee = (1 / mu) * ((v^2 - mu / r) * rr - r * v_r * vv);
e1 = norm(ee);

if e1 <= 1e-4
    e = 0;
else
   e = e1;
end

% Pericentre argument:

if e==0
            om=0;
end

if e~=0
    if ee(3) >= 0
            om = acos(dot(NN, ee) / (N * e));
    else
        om = 2 * pi - acos(dot(NN, ee) / (N * e));
    end

end

% True anomaly;

if e == 0
    if v_r >=0
            th= acos(dot(NN, rr) / (N * r));
    else  
            th = 2 * pi - acos(dot(NN, rr) / (N * r));
    end
end

if e ~= 0
        if v_r >=0
            th = acos(dot(ee, rr) / (e * r));
        else
            th = 2*pi - acos(dot(ee, rr) / (e * r));
        end
end
 
end

