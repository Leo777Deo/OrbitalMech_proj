function acc_pert_vec = fun_a3_J2Moon (t,s)

% Evaluation of the acceleration due to J2+Moon perturbations in RSW frame 
% -----------------------------------------------------------------------------------------
% PROTOTYPE:
% acc_pert_vec = fun_a3_J2Moon (t,s)
% -----------------------------------------------------------------------------------------
% INPUT:
% t               [1x1]       Time                                             [s]
% s               [1x6]       Vector of the orbital elements (a,e,i,OM,om,th)  [km][-][rad]
% -----------------------------------------------------------------------------------------
% OUTPUT:
% acc_pert_vec    [3x1]       Acceleration J2 and Moon                         [km/s^2]
% -----------------------------------------------------------------------------------------
% CONTRIBUTORS:
% Viola Poverini
% Gianluca Perusini
% -----------------------------------------------------------------------------------------
% VERSIONS:
% 2023-12-18: First Version
% 2024-01-01: Last version 
% -----------------------------------------------------------------------------------------

%Parameters 
mu_E = astroConstants(13); %Earth's gravitational parameter 
mu_M = astroConstants(20); %Moon's gravitational parameter
mu = mu_E; %Gravitational parameter of the primary 
R_E = astroConstants(23); %Earth's radius
J2 = astroConstants(9); 

%Orbital elements
a = s(1,:);
e = s(2,:);
i = s(3,:);
OM = s(4,:);
om = s(5,:);
th = s(6,:);

%Conversion of keplerian elements into Cartesian coordinates 
[r, ~] = kep2car(a,e,i,OM,om,th,mu); 

% Distance of the S/C from the primary
rnorm = norm(r);
modr = rnorm; 

% Position of the Moon 
mjd2000i = date2mjd2000([2015,08,03,00,00,00]); %initial day
date = mjd2000i + t/86400; %initial date 2015,08,03,00:00:00
[r_EM, ~] = ephMoon(date);
r_EM = r_EM'; 

% Distance of the Moon from the primary 
modr_EM = norm(r_EM); 

% Acceleration's evaluation 

r_SCM = r_EM-r; 
modr_SCM = norm(r_SCM);

a3Bi = mu_M*((r_SCM(1)/(modr_SCM)^3)-(r_EM(1)/(modr_EM)^3));
a3Bj = mu_M*((r_SCM(2)/(modr_SCM)^3)-(r_EM(2)/(modr_EM)^3));
a3Bk = mu_M*((r_SCM(3)/(modr_SCM)^3)-(r_EM(3)/(modr_EM)^3));


[arM, asM, awM] = rot_ECI_RSW(a3Bi, a3Bj, a3Bk,s); 


arJ2 = (1-3*(sin(i))^2*(sin(th+om))^2); 
asJ2 = ((sin(i))^2*(sin(2*(th+om))));
awJ2 = (sin(2*i)*sin(th+om)); 
p = a*(1-e^2);
h = sqrt(mu*p);
r = p/(1+e*cos(th));
acc_pert_vec = ((-3/2)*(J2*mu*R_E^2)/(r^4))*[arJ2 asJ2 awJ2] + [arM asM awM]; 

end 