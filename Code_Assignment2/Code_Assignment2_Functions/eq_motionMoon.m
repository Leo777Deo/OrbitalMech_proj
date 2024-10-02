function ds = eq_motionMoon( t, s, fun_a3_J2Moon)

% Evaluation of the acceleration due to J2+Moon perturbations in ECI frame 
% ----------------------------------------------------------------------------------------------
% PROTOTYPE:
% ds = eq_motionMoon( t, s, fun_a3_J2Moon)
% ----------------------------------------------------------------------------------------------
% INPUT:
% t               [1x1]       Time                                              [s]
% s               [6x1]       Vector of the orbital elements (a,e,i,OM,om,th)   [km][-][rad]
% fun_a3_J2Moon   Function to evaluate the J2+Moon acceleration in RSW
% ----------------------------------------------------------------------------------------------
% OUTPUT:
% ds              [6x1]       Derivative of the state                           [km/s][1/s][rad/s]
% ----------------------------------------------------------------------------------------------
% CONTRIBUTORS:
% Viola Poverini
% Leo De Luca
% ----------------------------------------------------------------------------------------------
% VERSIONS:
% 2023-12-13: First Version
% 2024-01-02: Last version 
% ----------------------------------------------------------------------------------------------

% Parameters
mu = astroConstants(13); %gravitational parameter of the Earth 
J2 = astroConstants(9);
R_E = astroConstants(23); %radius of the Earth 

%Orbital elements
a = s(1,:);
e = s(2,:);
i = s(3,:);
OM = s(4,:);
om = s(5,:);
th = s(6,:);

% Evaluation the perturbing accelerations
acc_pert_vec = fun_a3_J2Moon(t,s);

% Evaluation of the Gauss's planetary equations in RSW frame
p = a*(1-e^2);
h = sqrt(mu*p);
r = p/(1+e*cos(th)); 

da = ((2*a^2)/h)*(e*sin(th)*acc_pert_vec(1)+(p/r)*acc_pert_vec(2)); 
dh = r*acc_pert_vec(2);
de = (1/h)*(p*sin(th)*acc_pert_vec(1)+((p+r)*cos(th)+r*e)*acc_pert_vec(2));
di = ((r*cos(th+om))/h)*acc_pert_vec(3); 
dOM = ((r*sin(th+om))/(h*sin(i)))*acc_pert_vec(3); 
dom = (1/(h*e))*(-p*cos(th)*acc_pert_vec(1)+(p+r)*sin(th)*acc_pert_vec(2))-((r*sin(th+om)*cos(i))/(h*sin(i)))*acc_pert_vec(3);
dth = (h/r^2)+(1/(h*e))*(p*cos(th)*acc_pert_vec(1)-(p+r)*sin(th)*acc_pert_vec(2));

ds = [da; de; di; dOM; dom; dth];

end