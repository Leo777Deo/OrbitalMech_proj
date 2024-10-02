function [deltaVtot] = computedeltaVtot(t,ibody1,ibody2,ibody3,mu_primary,mu_planet2,r_planet2)
%
% COMPUTEDELTAVTOT   Function that calculates the total delta v to be provided 
%                    assuming that the initial heliocentric orbit is equal 
%                    to that of the starting planet, that we have a fly-by of 
%                    planet 2, and that we enter the final orbit assuming 
%                    that the final heliocentric orbit is equal to that of 
%                    the arrival asteroid.
%  
% [deltaVtot] = computedeltaVtot(t,ibody1,ibody2,ibody3,mu_S,mu_planet2,r_planet2)
% 
% Input arguments:
% ------------------------------------------------------------------------------------------
% t           [3x1]       Time vector with departure,flyby,arrival (mjd2000)     [s]
% ibody1      [1x1]       Integer number identifying the celestial body (< 11)   [-]
% ibody2      [1x1]       Integer number identifying the celestial body (< 11)   [-]
% ibody3      [1x1]       Integer number identifying the celestial body (< 11)   [-]           
%                             1:   Mercury
%                             2:   Venus
%                             3:   Earth
%                             4:   Mars
%                             5:   Jupiter
%                             6:   Saturn
%                             7:   Uranus
%                             8:   Neptune
%                             9:   Pluto
%                             10:  Sun  
% mu_primary  [1x1]       Gravitational parameter of the primary                 [km^3/s^2]
% mu_planet2  [1x1]       Gravitational parameter of the planet of the flyby     [km^3/s^2]
% r_planet2   [1x1]       Radius of the planet of the flyby                      [km]
%
%
% Output arguments:
% ------------------------------------------------------------------------------------------
% deltaVtot  [1x1]        Magnitude of the manoeuvres required                   [km/s]
%
%
%
% AUTHOR
%   Andrea Barbiera
%
%   Gianluca Perusini


ti=t(1);
tfb=t(2);
tf=t(3);
orbitType=0;
Nrev=0;
Ncase=0;
optionsLMR=0;

%Transfer from 1 to 2
[kep1,~] = uplanet(ti, ibody1);
[R1, V1] = kep2car (kep1(1),kep1(2),kep1(3),kep1(4),kep1(5),kep1(6), mu_primary);
[kep2,~] = uplanet(tfb, ibody2);
[R2, V2] = kep2car (kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6), mu_primary);
TOFt1=(tfb-ti)*24*60*60;
[~,~,~,~,VI,VFB1,~,~] = lambertMR(R1,R2,TOFt1,mu_primary,orbitType,Nrev,Ncase,optionsLMR);  
VI=VI';  
VFB1=VFB1';
deltaV11=VI-V1;       %IN REALTA' DOVREI TROVARE IL DELTAVP ALL INJECTION
deltaVt1=norm(deltaV11);

%Transfer from 2 to 3
[kep_A,~,~] = ephNEO(tf,ibody3);
[R_A, V_A] = kep2car (kep_A(1),kep_A(2),kep_A(3),kep_A(4),kep_A(5),kep_A(6), mu_primary);
TOFt2=(tf-tfb)*24*60*60;
[~,~,~,~,VFB2,VF,~,~] = lambertMR(R2,R_A,TOFt2,mu_primary,orbitType,Nrev,Ncase,optionsLMR);  
VFB2=VFB2';  
VF=VF';
deltaV22=V_A-VF;
deltaVt2=norm(deltaV22);


%Fly-by of 2
vinfa=VFB1-V2;
vinfan=norm(vinfa);
vinfd=VFB2-V2;
vinfdn=norm(vinfd);
delta=acos(dot(vinfa,vinfd)/(vinfan*vinfdn));

delta_minus=@(rp)   2*asin(1/(1+rp*vinfan^2/mu_planet2));
delta_plus=@(rp)    2*asin(1/(1+rp*vinfdn^2/mu_planet2));
rp_solve=@(rp)      delta_minus(rp)/2+ delta_plus(rp)/2-delta;
rp0=1*r_planet2;
[rp,~,~,~] = fzero(rp_solve, rp0,optimset('Display','none'));

v1=sqrt(vinfan^2+2*mu_planet2/rp);     
v2=sqrt(vinfdn^2+2*mu_planet2/rp);    
deltavp2=v2-v1;    
deltavp2n=norm(deltavp2);
deltaVfb=deltavp2n;



deltaVtot=deltaVt1+deltaVfb+deltaVt2;



end

