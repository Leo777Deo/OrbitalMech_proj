%% --------------PLANETARY EXPLORER MISSION---------------------

%% ---------------------Assignment 2---------------------------
%
% By running this script the whole assignment 2 is analysed.
%
%
clear all; 
clc;

%########################################################################
%########################################################################

%Parameters 
mu_E = astroConstants(13); %Earth's gravitational parameter [km^3/s^2]
mu = mu_E; 
mu_M = astroConstants(20); %Moon's gravitational parameter [km^3/s^2]
R_E = astroConstants(23); %Earth's radius [km]
J2 = astroConstants(9);

% Initial conditions
t0 = 0; %[s] 
tetaG0 = 0; %[rad]  
wE = (15.04*pi/180)/3600; %[rad/s]
a0 = 39689; %[km] 
e0 = 0.8264; 
i0 = deg2rad(21.8645); %[rad] 
OM0 = deg2rad(180); %[rad] 
om0 = deg2rad (90); %[rad] 
th0 = deg2rad(0); %[rad] 

%########################################################################
%########################################################################

%% UNPERTURBED 2BP NOMINAL ORBIT 

%------------------------------------------------------------------------
%%ORBIT PROPAGATION OVER SIX MONTHS

% Initial conditions
[r0, v0] = kep2car(a0, e0, i0, OM0, om0, th0, mu); % Startng from Orbital Elements in km or rad 
y0 = [r0; v0];
Tp = 2*pi*sqrt(a0^3/mu); % Orbital period [1/s]
z = 202.029; %Period's repetitions
tspan = [0:3:z*Tp]; %6 months (184 days)

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y), tspan, y0, options );

% Values of r and v along the orbit
r = Y(:,1:3);
v = Y(:,4:6);
modr = vecnorm(r,2,2); 

%------------------------------------------------------------------------
%%ORBIT PLOT
h = cross(r0,v0); %Angular momentum 
B = ((cross(v0,h)) - (mu_E/(norm(r0)))*r0)./10; %Scaled eccentricity vector 
fig1 = figure;
scatter3( Y(:,1), Y(:,2), Y(:,3),3,modr) %plot3D
hold on 
quiver3(0,0,0,cos(OM0)*40000,sin(OM0)*40000,0,'y','LineWidth',1.5); %Nodal line
hold on
quiver3(0,0,0,40000,0,0,'r','LineWidth',1.5); %Gamma line
hold on
quiver3(0,0,0,0,40000,0,'g','LineWidth',1.5); %Y axis
hold on
quiver3(0,0,0,0,0,40000,'b','LineWidth',1.5); %Rotation axis of E
hold on
quiver3(0,0,0,B(1,1),B(2,1),B(3,1),'k','LineWidth',1.5); %Apse Line
hold on
quiver3(0,0,0,h(1,1)/2,h(2,1)/2,h(3,1)/2,'m','LineWidth',1.5); %Scaled h 
hold on 
[xs,ys,zs] = sphere;
xs = R_E*xs;
ys = R_E*ys;
zs = R_E*zs; 
s = surf(xs,ys,zs);
set(s,'edgecolor','none')
sph1 = findobj('Type','surface'); 
im1 = imread('earth.png'); % Plot the Earth on the sphere 
set(sph1,'CData',im1,'FaceColor','texturemap','FaceAlpha',1);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit without perturbations');
axis equal;
grid on
cb = colorbar('eastoutside','Ticks',[round(min(modr)):5000:round(max(modr))]);
cb.Title.String = 'Spacecraft position [km]';
cb.Ruler.Exponent = 0;

%------------------------------------------------------------------------
%% GROUND TRACK
% Declination and Right ascension 
Declination = asin(r(:,3)./modr); %rad
RA = atan2(r(:,2),r(:,1));
tetaG = tetaG0 + wE.*tspan; %Longitude of Greenwich meridian  
Lon = RA - tetaG'; %rad
Longitude = wrapToPi(Lon);
Londeg = Longitude*180/pi; %deg 
Lat = Declination; %rad
Latitude = wrapToPi(Lat);
Latdeg = Latitude*180/pi; %deg 

%Plot of the Ground Track 
fig2 = figure;
axis([-180 180 -90 90]);
im2 = imread('earth.png');
map_image = image(-180:180, -90:90, im2); 
hold on
scatter(Londeg,Latdeg,10,modr) %100 periods
hold on 
plot(Londeg(1,1),Latdeg(1,1),'Color','g','Marker','o','MarkerSize',16); %Start point 
hold on
plot(Londeg(end,1),Latdeg(end,1),'Color','g','Marker','square','MarkerSize',16); %End point
ax = gca;
ax.XTick = -180:30:180;
ax.YDir = "normal";
ax.YTick = -90:30:90;
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
legend('','Start','End')
title ('Ground Track of the nominal unperturbed orbit')
cb = colorbar('eastoutside','Ticks',[round(min(modr)):5000:round(max(modr))]);
cb.Title.String = 'Spacecraft position [km]';
cb.Ruler.Exponent = 0;

%########################################################################

%% PERTURBED J2+MOON 2BP NOMINAL ORBIT

%------------------------------------------------------------------------
%%ORBIT PROPAGATION

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[ TP, YP ] = ode113( @(tp,yp) ode_2bp_J2Moon(tp,yp, @(tp,yp) fun_a3B(tp,yp)), tspan, y0, options );

% Values of r and v along the orbit
rp = YP(:,1:3);
vp = YP(:,4:6);
modrp = vecnorm(rp,2,2);

%------------------------------------------------------------------------
%%ORBIT PLOT
fig3 = figure;
plot3( YP(:,1), YP(:,2), YP(:,3), 'm' ) %plot3D
hold on 
plot3 (Y(:,1), Y(:,2), Y(:,3),'b', 'Linewidth', 3) %First unperturbed orbit
hold on
[xs,ys,zs] = sphere;
xs = R_E*xs;
ys = R_E*ys;
zs = R_E*zs; 
s = surf(xs,ys,zs);
set(s,'edgecolor','none')
sph1 = findobj('Type','surface'); 
im1 = imread('earth.png'); %P lot the Earth on the sphere 
set(sph1,'CData',im1,'FaceColor','texturemap','FaceAlpha',1);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem nominal orbit with perturbations');
axis equal;
grid on

%------------------------------------------------------------------------
%%GROUND TRACK
% Declination and Right ascension 
Declinationp = asin(rp(:,3)./modrp); % [rad]
RAp = atan2(rp(:,2),rp(:,1));
tetaG = tetaG0 + wE.*tspan; % Longitude of Greenwich meridian at initial time and omegaE Earth's rotation speed 
Lonp = RAp - tetaG'; % [rad]
Longitudep = wrapToPi(Lonp);
Londegp = Longitudep*180/pi; % [deg] 
Latp = Declinationp; % [rad]
Latitudep = wrapToPi(Latp);
Latdegp = Latitudep*180/pi; % [deg]

%Plot of the Ground Track 
fig4 = figure;
axis([-180 180 -90 90]);
im2 = imread('earth.png');
map_image = image(-180:180, -90:90, im2); 
hold on
plot(Londeg,Latdeg,'.b','LineStyle','none'); %Nominal Unperturbed 
hold on
plot(Londegp,Latdegp,'.m','LineStyle','none'); %Nominal Perturbed 
hold on  
plot(Londeg(1,1),Latdeg(1,1),'Color','g','Marker','o','MarkerSize',16,'LineWidth',3); %Start point 
hold on
plot(Londeg(end,1),Latdeg(end,1),'Color','g','Marker','square','MarkerSize',16,'LineWidth',3); %End point
hold on
plot(Londegp(1,1),Latdegp(1,1),'Color','y','Marker','o','MarkerSize',16,'LineWidth',3);
hold on
plot(Londegp(end,1),Latdegp(end,1),'Color','y','Marker','square','MarkerSize',16,'LineWidth',3);
ax = gca;
ax.XTick = -180:30:180;
ax.YDir = "normal";
ax.YTick = -90:30:90;
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
title ('Ground Track of the nominal orbit with J2+Moon perturbations')
legend('Unperturbed Nominal orbit', 'Perturbed Nominal orbit','Start for the Unperturbed','End for the Unperturbed','Start for the Perturbed','End for the Perturbed');

%########################################################################
%########################################################################

%% UNPERTURBED 2BP MODIFIED ORBIT 
%Semi-major axis for repeating ground track
m = 5;
k = 6; 
Tr = (2*pi/wE)*(m/k);
ar = (mu*(Tr/(2*pi))^2)^(1/3); %Semi-major axis of the repeating orbit

%------------------------------------------------------------------------
%%ORBIT PROPAGATION
[r0r, v0r] = kep2car(ar, e0, i0, OM0, om0, th0, mu); %in rad 
y0r = [r0r; v0r];
p=100; 
tspanr = [0:3:p*Tr];

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[ TR, YR ] = ode113( @(tr,yr) ode_2bp(tr,yr), tspanr, y0r, options );

% Values of r and v along the orbit
rr = YR(:,1:3);
vr = YR(:,4:6);
modrr = vecnorm(rr,2,2);

%------------------------------------------------------------------------
%%GROUND TRACK 
% Declination and Right ascension 
Declinationr = asin(rr(:,3)./modrr); %rad
RAr = atan2(rr(:,2),rr(:,1));
tetaGr = tetaG0 + wE.*tspanr; %longitude of Greenwich meridian at initial time and omegaE Earth's rotation speed 
Lonr = RAr - tetaGr'; %rad
Longituder = wrapToPi(Lonr);
Londegr = Longituder*180/pi; %deg 
Latr = Declinationr; %rad
Latituder = wrapToPi(Latr);
Latdegr = Latituder*180/pi; %deg 

%Plot of the Ground Track 
fig5 = figure;
axis([-180 180 -90 90]);
im2 = imread('earth.png');
map_image = image(-180:180, -90:90, im2); 
hold on
scatter(Londegr,Latdegr,7,modrr);
plot(Londegr(1,1),Latdegr(1,1),'Color','g','Marker','o','MarkerSize',16); %Start point
hold on
plot(Londegr(end,1),Latdegr(end,1),'Color','g','Marker','square','MarkerSize',16); %End point
ax = gca;
ax.XTick = -180:30:180;
ax.YDir = "normal";
ax.YTick = -90:30:90;
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
legend('','Start','End')
title ('Repeating Ground Track for the unperturbed orbit')
cb = colorbar('eastoutside','Ticks',[round(min(modrr)):5000:round(max(modrr))]);
cb.Title.String = 'Spacecraft position [km]';
cb.Ruler.Exponent = 0;

%########################################################################

%% PERTURBED J2+MOON 2BP MODIFIED ORBIT

%------------------------------------------------------------------------
%%ORBIT PROPAGATION
% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[ TRP, YRP ] = ode113( @(trp,yrp) ode_2bp_J2Moon(trp,yrp, @(trp,yrp) fun_a3B(trp,yrp)), tspanr, y0r, options );

% Values of r and v along the orbit
rrp = YRP(:,1:3);
vrp = YRP(:,4:6);
modrrp = vecnorm(rrp,2,2);

%------------------------------------------------------------------------
%% GROUND TRACK 
% Declination and Right ascension 
Declinationrp = asin(rrp(:,3)./modrrp); %rad
RArp = atan2(rrp(:,2),rrp(:,1));
tetaGr = tetaG0 + wE.*tspanr; %longitude of Greenwich meridian at initial time and omegaE Earth's rotation speed 
Lonrp = RArp - tetaGr'; %rad
Longituderp = wrapToPi(Lonrp);
Londegrp = Longituderp*180/pi; %deg 
Latrp = Declinationrp; %rad
Latituderp = wrapToPi(Latrp);
Latdegrp = Latituderp*180/pi; %deg 

%Plot of the Ground Track 
fig6 = figure;
axis([-180 180 -90 90]);
im2 = imread('earth.png');
map_image = image(-180:180, -90:90, im2); 
hold on
plot(Londegrp,Latdegrp,'.r','LineStyle','none');
hold on
plot(Londegr,Latdegr,'.c','LineStyle','none');
hold on  
plot(Londegr(1,1),Latdegr(1,1),'Color','g','Marker','o','MarkerSize',16,'LineWidth',3); %Start point 
hold on
plot(Londegr(end,1),Latdegr(end,1),'Color','g','Marker','square','MarkerSize',16,'LineWidth',3); %End point
hold on
plot(Londegrp(1,1),Latdegrp(1,1),'Color','y','Marker','o','MarkerSize',16,'LineWidth',3);
hold on
plot(Londegrp(end,1),Latdegrp(end,1),'Color','y','Marker','square','MarkerSize',16,'LineWidth',3);
ax = gca;
ax.XTick = -180:30:180;
ax.YDir = "normal";
ax.YTick = -90:30:90;
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
title ('Repeating Ground Track with J2+Moon perturbations')
legend('Perturbed Repeating orbit','Unperturbed Repeating orbit','Start for the Unperturbed','End for the Unperturbed','Start for the Perturbed','End for the Perturbed');

%########################################################################
%########################################################################

%% GAUSS'S PLANETARY EQUATIONS

% Initial conditions
s0 = [a0; e0; i0; OM0; om0; th0]; 

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration in Keplerian elements 
[ TS, S ] = ode113( @(ts,s) eq_motionMoon( ts, s, @(ts,s) fun_a3_J2Moon(ts,s) ), tspan, s0, options );


%Conversion of the values from Cartesian coordinates to Keplerian elements 
Jr = YP(:,1:3)'; 
Jv = YP(:,4:6)';
for k=1:length(TS)
[ac(k),ec(k),ic(k),OMc(k),omc(k),thc(k)] = car2kep(Jr(:,k), Jv(:,k), mu);
end 

%% HISTORY OF THE KEPLERIAN ELEMENTS AND RELATIVE ERRORS OVER SIX MONTH
%------------------------------------------------------------------------
%% Semi-major axis 
fig7 = figure;
sm = 27*24*3600+7*3600+43*60+12; %Sidereal month 
subplot(1,2,1)
plot(tspan,ac(:),'m') %Cartesian coordinates
hold on 
plot(tspan,S(:,1),'Color',[0.4940 0.1840 0.5560]); %Gauss's planetary equations 
hold on 
plot(tspan,movmean(ac(:),sm+Tp),'g','LineWidth',1.5);
ylabel('a [km]')
xlabel('time [Tp]')
legend('Cartesian','Gauss equations','Secular (filtered)');
set(gca,'XTickLabel',[0,25,51,76,102,127,152,178,203]);
title('Semi-major axis')
grid on 

%Error on the semi-major axis
subplot(1,2,2) 
erra = abs(ac(:)-S(:,1))/a0; 
semilogy(tspan,erra,'Color',[0.9290 0.6940 0.1250])
ylabel('|aCar-aGauss|/a0 [-]')
xlabel('time [Tp]')
set(gca,'XTickLabel',[0,25,51,76,102,127,152,178,203]);
title('Error over the semi-major axis')
grid on 

%------------------------------------------------------------------------
%% Eccentricity 
fig8 = figure; 
subplot(1,2,1)
plot(tspan,ec(:),'m'); %Cartesian coordinates
hold on
plot(tspan,S(:,2),'Color',[0.4940 0.1840 0.5560]); %Gauss's planetary equations 
hold on 
plot(tspan,movmean(ec(:),sm),'g','LineWidth',1.5);
ylabel('e [-]')
xlabel('time [Tp]')
legend('Cartesian','Gauss equations','Secular (filtered)')
set(gca,'XTickLabel',[0,25,51,76,102,127,152,178,203]);
title('Eccentricity')
grid on 

%Error on the eccentricity 
subplot(1,2,2) 
erre = abs(ec(:)-S(:,2)); 
semilogy(tspan,erre,'Color',[0.9290 0.6940 0.1250])
ylabel('|eCar-eGauss| [-]')
xlabel('time [Tp]')
set(gca,'XTickLabel',[0,25,51,76,102,127,152,178,203]);
title('Error over the eccentricity')
grid on 

%------------------------------------------------------------------------
%% Inclination 
fig9 = figure;
subplot(1,2,1)
plot(tspan,ic(:)*180/pi,'m'); %Cartesian coordinates
hold on
plot(tspan,S(:,3)*180/pi,'Color',[0.4940 0.1840 0.5560]); %Gauss's planetary equations 
hold on
plot(tspan,movmean(ic(:)*180/pi,sm),'g','LineWidth',1.5); %Secular perturbation of the Moon
ylabel('i [deg]')
xlabel('time [Tp]')
legend('Cartesian','Gauss equations','Secular (filtered)','Location','northwest')
set(gca,'XTickLabel',[0,25,51,76,102,127,152,178,203]);
title('Inclination')
grid on 

%Error on the inclination
subplot(1,2,2)
erri = abs(ic(:)*180/pi-S(:,3)*180/pi)/(2*180); 
semilogy(tspan,erri,'Color',[0.9290 0.6940 0.1250])
ylabel('|iCar-iKep|/|2π| [-]')
xlabel('time [Tp]')
set(gca,'XTickLabel',[0,25,51,76,102,127,152,178,203]);
title('Error over the inclination')
grid on 

%------------------------------------------------------------------------
%% Right ascention of the ascending node 
fig10 = figure; 
subplot(1,2,1)
plot(tspan,OMc(:)*180/pi,'m'); %Cartesian coordinates
hold on
plot(tspan,S(:,4)*180/pi,'Color',[0.4940 0.1840 0.5560]); %Gauss's planetary equations 
hold on 
plot(tspan,movmean(OMc(:)*180/pi,Tp),'g','LineWidth',1.5);
ylabel('Ω [deg]')
xlabel('time [Tp]')
legend('Cartesian','Gauss equations','Secular (filtered)')
set(gca,'XTickLabel',[0,25,51,76,102,127,152,178,203]);
title('Right ascention of the ascending node')
grid on 

%Error on the RAAN
subplot(1,2,2) 
errOM = abs(OMc(:)*180/pi-S(:,4)*180/pi)/(2*180); 
semilogy(tspan,errOM,'Color',[0.9290 0.6940 0.1250])
ylabel('|ΩCar-ΩGauss|/|2π| [-]')
xlabel('time [Tp]')
set(gca,'XTickLabel',[0,25,51,76,102,127,152,178,203]);
title('Error over the RAAN')
grid on 

%------------------------------------------------------------------------
%% Argument of pericentre 
fig11 = figure; 
subplot(1,2,1)
plot(tspan,omc(:)*180/pi,'m'); %Cartesian coordinates
hold on
plot(tspan,S(:,5)*180/pi,'Color',[0.4940 0.1840 0.5560]); %Gauss's planetary equations 
hold on
plot(tspan,movmean(omc(:)*180/pi,sm+Tp),'g','LineWidth',1);
ylabel('ω [deg]')
xlabel('time [Tp]')
legend('Cartesian','Gauss equations','Secular (filtered)','Location','best')
set(gca,'XTickLabel',[0,25,51,76,102,127,152,178,203]);
title('Argument of pericentre')
grid on 

%Error on the Argument of pericentre  
subplot(1,2,2) 
errom = abs(omc(:)*180/pi-S(:,5)*180/pi)/(2*180); 
semilogy(tspan,errom,'Color',[0.9290 0.6940 0.1250])
ylabel('|ωCar-ωGauss|/|2π| [-]')
xlabel('time [Tp]')
set(gca,'XTickLabel',[0,25,51,76,102,127,152,178,203]);
title('Error over the argument of pericentre')
grid on 

%------------------------------------------------------------------------
%% True anomaly 
fig12 = figure;
subplot(1,2,1)
plot(tspan(1,65:end),(unwrap(thc(1,65:end)))*180/pi,'m') %Cartesian coordinates
hold on
plot(tspan(1,65:end),(unwrap(S(65:end,6)))*180/pi,'Color',[0.4940 0.1840 0.5560]); %Gauss's planetary equations 
hold on
plot(tspan(1,65:end),movmean(unwrap(thc(1,65:end))*180/pi,Tp),'g','LineWidth',1);
ylabel('th [deg]')
xlabel('time [Tp]')
legend('Cartesian','Gauss equations','Secular (filtered)','Location','best')
set(gca,'XTickLabel',[0,25,51,76,102,127,152,178,203]);
title('True anomaly')
grid on 

%Error on the True Anomaly  
subplot(1,2,2)
errth = abs((unwrap(thc(1,65:end))')*180/pi-(S(65:end,6))*180/pi)./abs((S(65:end,6))*180/pi); 
semilogy(tspan(1,65:end),errth,'Color',[0.9290 0.6940 0.1250])
ylabel('|thCar-thGauss|/|thGauss| [-]')
xlabel('time [Tp]')
set(gca,'XTickLabel',[0,25,51,76,102,127,152,178,203]);
title('Error over the true anomaly')
grid on 

%########################################################################
%########################################################################

%% PERTURBED REAL ORBIT 

%------------------------------------------------------------------------
% Initial conditions
t0 = 0; %[s] 
mjd2000i = date2mjd2000([2015,08,03,00,00,00]); %Initial day
mjd2000f = date2mjd2000([2021,08,03,00,00,00]); %Final day
r02 = [-1.990983933520086E+04; 3.456648088166495E+04; 4.928636246930589E+03]; %Initial position delta 4 R\B
v02 = [-2.728512058041876E+00; 1.482057828630650E+00; -3.871912720354164E-01]; %Initial velocity
y02 = [r02;v02];
[a02,e02,i02,OM02,om02,th02] = car2kep(r02, v02, mu);
s02 = [a02; e02; i02; OM02; om02; th02]; 

%------------------------------------------------------------------------
%%ORBIT PROPAGATION IN CARTESIAN COORDINATES
L = 2190; %days between 2015,08,03-2021,08,03
tspanp2 = [0:3600:L*86400]; % sec

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% Perform the integration
[ TP2, YP2 ] = ode113( @(tp2,yp2) ode_2bp_J2Moon(tp2,yp2, @(tp2,yp2) fun_a3B(tp2,yp2)), tspanp2, y02, options );

% Values of r and v along the orbit
rp2 = YP2(:,1:3);
vp2 = YP2(:,4:6);
modrp2 = vecnorm(rp2,2,2); 

Jrp2 = YP2(:,1:3)'; 
Jvp2 = YP2(:,4:6)';

for t=1:length(TP2)
[ac2(t),ec2(t),ic2(t),OMc2(t),omc2(t),thc2(t)] = car2kep(Jrp2(:,t), Jvp2(:,t), mu);
end 

%------------------------------------------------------------------------
%%INTEGRATION WITH GAUSS'S PLANETARY EQUATIONS 
[ TS2, S2 ] = ode113( @(ts2,s2) eq_motionMoon( ts2, s2, @(ts2,s2) fun_a3_J2Moon(ts2,s2) ), tspanp2, s02, options );
load('KEPEPH.mat'); % loads the Ephemerides data from Horizon

%------------------------------------------------------------------------
%% ORBIT PLOT 
fig13 = figure;
plot3(YP2(:,1), YP2(:,2), YP2(:,3),'m') % plot of real orbit
hold on
[xs,ys,zs] = sphere;
xs = R_E*xs;
ys = R_E*ys;
zs = R_E*zs; 
s = surf(xs,ys,zs);
set(s,'edgecolor','none')
sph1 = findobj('Type','surface'); 
im1 = imread('earth.png'); % Plot the Earth on the sphere 
set(sph1,'CData',im1,'FaceColor','texturemap','FaceAlpha',1);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit without perturbations');
axis equal;
grid on


%% HISTORY OF THE KEPLERIAN ELEMENTS OVER TWO YEARS

%-----------------------------------------------------------------------
%% Semi-major axis 
fig14 = figure;
plot(tspanp2,S2(:,1),'Color',[0.9290 0.6940 0.1250]); %Gauss's planetary equations
hold on 
plot(tspanp2,EPH(1:52561,5),'Color',[0.4940 0.1840 0.5560])
ylabel('a [km]')
xlabel('time [MJD2000]')
set(gca,'XTickLabel',[5692.5,5923.9,6155.4,6386.9,6618.4,6849.9,7081.3,7312.8,7544.3,7775.8,8007.3]);
title('Semi-major axis')
legend('Propagated','Delta 4 R/B','Location','best')
grid on 

%------------------------------------------------------------------------
%% Eccentricity 
fig15= figure; 
subplot(1,2,1)
plot(tspanp2,S2(:,2),'Color',[0.9290 0.6940 0.1250]); %Gauss's planetary equations
hold on
plot(tspanp2,EPH(1:52561,1),'Color',[0.4940 0.1840 0.5560])
ylabel('e [-]')
xlabel('time [MJD2000]')
set(gca,'XTickLabel',[5692.5,6271.2,6849.9,7428.6,8007.3]);
title('Eccentricity')
legend('Propagated','Delta 4 R/B','Location','northeast')
grid on 

% Inclination 
subplot(1,2,2)
plot(tspanp2,S2(:,3)*180/pi,'Color',[0.9290 0.6940 0.1250]); %Gauss's planetary equations
hold on
plot(tspanp2,EPH(1:52561,2),'Color',[0.4940 0.1840 0.5560])
ylabel('i [deg]')
xlabel('time [MJD2000]')
set(gca,'XTickLabel',[5692.5,6271.2,6849.9,7428.6,8007.3]);
title('Inclination')
legend('Propagated','Delta 4 R/B','Location','best')
grid on 

%------------------------------------------------------------------------
%% RAAN 
fig16 = figure; 
subplot(1,2,1)
plot(tspanp2,S2(:,4)*180/pi,'Color',[0.9290 0.6940 0.1250]); %Gauss's planetary equations
hold on
plot(tspanp2,EPH(1:52561,3),'Color',[0.4940 0.1840 0.5560])
ylabel('e [-]')
xlabel('time [MJD2000]')
set(gca,'XTickLabel',[5692.5,6271.2,6849.9,7428.6,8007.3]);xlabel('time [Tp]')
title('RAAN')
legend('Propagated','Delta 4 R/B','Location','best')
grid on 

% Argument of pericentre
subplot(1,2,2)
plot(tspanp2,S2(:,5)*180/pi,'Color',[0.9290 0.6940 0.1250]); %Gauss's planetary equations
hold on
plot(tspanp2,EPH(1:52561,4),'Color',[0.4940 0.1840 0.5560])
ylabel('AOP [deg]')
xlabel('time [MJD2000]')
set(gca,'XTickLabel',[5692.5,6271.2,6849.9,7428.6,8007.3]);
title('AOP')
legend('Propagated','Delta 4 R/B','Location','best')
grid on 
