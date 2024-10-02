function [vrot,deltav] = Rodrigues(v,delta,u)
%
% RODRIGUES The function gives velocity obtained from an initial velocity 
%           rotated by the counterclockwise delta angle around the versor u.
% 
% [vrot,deltav] = Rodrigues(v,delta,u)
%  
% Input arguments:
% ------------------------------------------------------------------------
%   v          [3x1]       Initial Velocity  [km/s]
%  delta       [1x1]       Turn Angle        [rad]
%   u          [3x1]       Axis of rotation
% Output arguments:
% ------------------------------------------------------------------------
%   vrot       [3x1]       Rotated Velocity  [km/s]
%  deltav      [3x1]       Difference between Final and Initial velocity [km/s]
%
%
% AUTHOR:
%   Andrea Barbiera


vrot=v*cos(delta)+cross(u,v)*sin(delta)+u*dot(u,v)*(1-cos(delta));
deltav=vrot-v;


end

