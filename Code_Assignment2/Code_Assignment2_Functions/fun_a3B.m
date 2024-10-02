function [a3B] = fun_a3B (t,y)

% Evaluation of the acceleration due to the third body perturbation of the Moon in ECI frame 
% ---------------------------------------------------------------------------------------
% PROTOTYPE:
% [a3B] = fun_a3B (t,y)
% ---------------------------------------------------------------------------------------
% INPUT:
% t               [1x1]       Time                                             [s]
% y               [6x1]       State of the body (rx, ry, rz, vx, vy, vz)       [km][km/s]
% ---------------------------------------------------------------------------------------
% OUTPUT:
% a3B             [3x1]       Acceleration caused by the presence of the Moon [km/s^2]
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
mu_M = astroConstants(20); %Moon 

% Position of the S/C 
r = y(1:3); 

% Distance of the S/C from the primary
rnorm = norm(r);
modr = rnorm; 

% Position of the Moon 
mjd2000i = date2mjd2000([2015,08,03,00,00,00]); %initial date
date = mjd2000i + t/86400; 
[r_EM, ~] = ephMoon(date);
r_EM = r_EM'; %distance of the Moon wrt the Earth (radius from Earth to the Moon) 

% Distance of the Moon from the primary 
modr_EM = norm(r_EM); 

% Acceleration in ECI frame
r_SCM = r_EM-r; %distance of the Moon wrt the spacecraft (radius from S/C to the Moon) 
modr_SCM = norm(r_SCM);

a3Bi = mu_M*((r_SCM(1)/(modr_SCM)^3)-(r_EM(1)/(modr_EM)^3));
a3Bj = mu_M*((r_SCM(2)/(modr_SCM)^3)-(r_EM(2)/(modr_EM)^3));
a3Bk = mu_M*((r_SCM(3)/(modr_SCM)^3)-(r_EM(3)/(modr_EM)^3));

a3B = [a3Bi a3Bj a3Bk]; 

end 