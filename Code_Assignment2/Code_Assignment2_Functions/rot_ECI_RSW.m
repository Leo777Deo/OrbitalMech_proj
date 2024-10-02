function [ar, as, aw] = rot_ECI_RSW(ai, aj, ak, s)

% Rotation from ECI frame to RSW frame 
% -----------------------------------------------------------------------------------------
% PROTOTYPE:
% [ar, as, aw] = rot_ECI_RSW(ai, aj, ak, s)
% ------------------------------------------------------------------------------------------------
% INPUT:
% ai               [1x1]       Acceleration due to the Moon on the i direction             [km/s^2]                                   
% aj               [1x1]       Acceleration due to the Moon on the j direction             [km/s^2]
% ak               [1x1]       Acceleration due to the Moon on the k direction             [km/s^2]
% s                [6x1]       Vector of the orbital elements (a,e,i,OM,om,th)             [km][-][rad]
% ------------------------------------------------------------------------------------------------
% OUTPUT:
% ar               [1x1]       Acceleration due to the Moon on the radial direction        [km/s^2]
% as               [1x1]       Acceleration due to the Moon on the transversal direction   [km/s^2]
% aw               [1x1]       Acceleration due to the Moon on the out-of-plane direction  [km/s^2]
% ------------------------------------------------------------------------------------------------
% CONTRIBUTORS:
% Viola Poverini
% Leo De Luca
% ------------------------------------------------------------------------------------------------
% VERSIONS:
% 2023-10-20: First Version 
% 2024-02-01: Last Version
% ------------------------------------------------------------------------------------------------

%Parameters 
i = s(3,:);
OM = s(4,:);
om = s(5,:);
th = s(6,:);

%Rotation matrices 
R3_OM = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
R1_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
R3_omth = [cos(om+th) sin(om+th) 0; -sin(om+th) cos(om+th) 0; 0 0 1]; 

%Product between the three rotation matrices 
R_ECI_orbit = R3_omth*R1_i*R3_OM; 

%Acceleration in the ECI frame
aECI = [ai;aj;ak];

%Acceleration in the orbit frame
a3BRSW = R_ECI_orbit*aECI; 

%Components of the acceleration in the orbit frame
ar = a3BRSW(1,:); 
as = a3BRSW(2,:); 
aw = a3BRSW(3,:);

end 
