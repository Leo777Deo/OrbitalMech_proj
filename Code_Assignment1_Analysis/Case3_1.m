clear;
close all;
clc;
%%

%% -------------------------3.1--------------------
% Initial values in mjd2000
%          Departure: 11927
%          FLYBY:     15942
%          Arrival:   20221

nd=18; % dummy variables
md=60;
nfb=2;
mfb=60;
na=2;
ma=60;

mjd2000di=11927;
datedi=mjd20002date(mjd2000di);
mjd2000df=mjd2000di+md*nd;                
datedf=mjd20002date(mjd2000df);

mjd2000fbi= 15942-60;
datefbi=mjd20002date(mjd2000fbi);
mjd2000fbf=mjd2000fbi+mfb*nfb;               
datefbf=mjd20002date(mjd2000fbf);

mjd2000ai= 20221-60;  
dateai=mjd20002date(mjd2000ai);
mjd2000af= mjd2000ai+ma*na;
dateaf=mjd20002date(mjd2000af);

%Sun
mu_S=astroConstants(4);
r_S=astroConstants(3);

%Saturn
ibody_S=6;
mu_Sat=astroConstants(16);
r_Sat=astroConstants(26);
[kep_Sat,~] = uplanet(mjd2000di, ibody_S);
[R_Sat, ~] = kep2car (kep_Sat(1),kep_Sat(2),kep_Sat(3),kep_Sat(4), ...
    kep_Sat(5),kep_Sat(6), mu_S);
T_Sat=2*pi*sqrt(kep_Sat(1)^3/mu_S);

%Jupiter
ibody_J=5;
mu_J=astroConstants(15);
r_J=astroConstants(25);
[kep_J,~] = uplanet(mjd2000fbi, ibody_J);
T_J=2*pi*sqrt(kep_J(1)^3/mu_S);
[R_J, ~] = kep2car (kep_J(1),kep_J(2),kep_J(3),kep_J(4),kep_J(5),kep_J(6), mu_S);



%Asteroid 1979XB
id_A=59;
r_A=r_Sat/100;    
[kep_A,~,~] = ephNEO(mjd2000ai,id_A);
T_A=2*pi*sqrt(kep_A(1)^3/mu_S);
[R_A, ~] = kep2car (kep_A(1),kep_A(2),kep_A(3),kep_A(4),kep_A(5),kep_A(6), mu_S);

% Computation of parabolic time as a stop condition for the loop

orbitType=0;
Nrev=0;
Ncase=0;
optionsLMR=0;

TOF1=1;
RI=R_Sat;
RF=R_J;
[~,~,~,~,~,~,TPAR1,~] = lambertMR(RI,RF,TOF1,mu_S,orbitType,Nrev,Ncase,optionsLMR);

TOF2=1;
RI=R_J;
RF=R_A;
[~,~,~,~,~,~,TPAR2,~] = lambertMR(RI,RF,TOF2,mu_S,orbitType,Nrev,Ncase,optionsLMR);

a1=(norm(R_Sat)+norm(R_J))/2;
TOFH1=pi*sqrt(a1^3/mu_S);  


a2=(norm(R_J)+norm(R_A))/2;
TOFH2=pi*sqrt(a2^3/mu_S);


%% Lambert
MU=mu_S;
orbitType=0;
Nrev=0;
Ncase=0;
optionsLMR=0;


vettdep=mjd2000di:nd:mjd2000df;
vettfb=mjd2000fbi:nfb:mjd2000fbf;
vettarr=mjd2000ai:na:mjd2000af;

y=1; % dummy variables
u=1;
p=1;

TOF1Matrix=NaN(md+1,mfb+1);
TOF2Matrix=NaN(mfb+1,ma+1);
deltaVMatrix1=NaN(md+1,mfb+1);
deltaVMatrix2=NaN(mfb+1,ma+1);
deltaVMatrixfb=NaN(md+1,mfb+1,ma+1);
deltaVTOT=NaN(md+1,mfb+1,ma+1);
FLAG=NaN(md+1,mfb+1,ma+1);


for q=mjd2000di:nd:mjd2000df
    
    [kep_Sat,~] = uplanet(q, ibody_S);   %Kep è riga
    [R_Sat, V_Sat] = kep2car (kep_Sat(1),kep_Sat(2),kep_Sat(3),kep_Sat(4), ...
    kep_Sat(5),kep_Sat(6), mu_S);  %R è colonna
    RI=R_Sat;




    for w=mjd2000fbi:nfb:mjd2000fbf
  
        
        TOFt1=(w-q)*24*60*60;
        
        if (TOFt1> 2920*(3600*24)  && TOFt1<=4015*(3600*24))                             %(TOFt1> 2920*(3600*24)  && TOFt1<=4015*(3600*24))   %TOFH1+-1.5 ANNI                                                
            [kep_J,~] = uplanet(w, ibody_J);
            [R_J, V_J] = kep2car (kep_J(1),kep_J(2),kep_J(3),kep_J(4),kep_J(5),kep_J(6), mu_S);
            RFB=R_J;
            [~,~,Et1,~,VI,VFB1,TPAR1,~] = lambertMR(RI,RFB,TOFt1,MU,orbitType,Nrev,Ncase,optionsLMR);
            % TPAR1Matrix(y,u)=TPAR1;
            VI=VI';
            VFB1=VFB1';
            deltaV1=VI-V_Sat;      %delta velocità che devo dare per entrare in t1
            deltaV1n=norm(deltaV1);
            if deltaV1n<4
            deltaVMatrix1(y,u)=deltaV1n;
            TOF1Matrix(y,u)=TOFt1/(3600*24); %così è in giorni

                for t=mjd2000ai:na:mjd2000af
                    
                    TOFt2=(t-w)*24*60*60;
                    
                    if (TOFt2>TPAR2 && TOFt2<=1.1*T_J)  % minimum is the parabolic and the maximum is through energy reasoning
                        [kep_A,~,~] = ephNEO(t,id_A);
                        [R_A, V_A] = kep2car (kep_A(1),kep_A(2),kep_A(3),kep_A(4),kep_A(5),kep_A(6), mu_S);
                        
                        RF=R_A;
                        [~,~,~,~,VFB2,VF,~,~] = lambertMR(RFB,RF,TOFt2,MU,orbitType,Nrev,Ncase,optionsLMR);
                        
                            VFB2=VFB2';
                            VF=VF';
                            deltaV2=V_A-VF;        %delta velocità che devo dare per entrare in orbita asteroide
                            deltaV2n=norm(deltaV2);
                            TOF2Matrix(u,p)=TOFt2/(3600*24); %così è in giorni
                            vinfa=VFB1-V_J;
                            vinfan=norm(vinfa);
                            vinfd=VFB2-V_J;
                            vinfdn=norm(vinfd);
                            delta=acos(dot(vinfa,vinfd)/(vinfan*vinfdn));
                            delta_minus=@(rp)   2*asin(1/(1+rp*vinfan^2/mu_J));
                            delta_plus=@(rp)    2*asin(1/(1+rp*vinfdn^2/mu_J));
                            rp_solve=@(rp)      delta_minus(rp)/2+ delta_plus(rp)/2-delta;
                            rp0=1*r_J;
                           
                            [rp,~,flag,~] = fzero(rp_solve, rp0,optimset('Display','none'));
                            FLAG(y,u,p)=flag;
                            altitude=rp-r_J;
                            if (altitude>5000 && flag==1)   
                                v1=sqrt(vinfan^2+2*mu_J/rp);     
                                v2=sqrt(vinfdn^2+2*mu_J/rp);    
            
            
                                deltavp=v2-v1;  
                                deltafbn=norm(deltavp);
                                deltaVMatrix2(u,p)=deltaV2n;
                                deltaVMatrixfb(y,u,p)=deltafbn;
                                deltaVTOT(y,u,p)=deltaV1n+deltaV2n+deltafbn;
                            end
                        
                    end
                    p=p+1;
                end
                p=1;
            
            end
        
        end
        u=u+1;

    end
    u=1;
    y=y+1;
end


somma=sum(sum(sum(FLAG==1)));

%%

deltaVTOTmin=min(min(min(deltaVTOT)));

[row,column,plane]=findND(deltaVTOT==deltaVTOTmin);

deltaVdott=deltaVMatrix1(row,column);
deltaVfbott=deltaVMatrixfb(row,column,plane);
deltaVaott=deltaVMatrix2(column,plane);

deltaVOTT=deltaVdott+deltaVfbott+deltaVaott;

mjd2000d=mjd2000di+nd*(row-1);
mjd2000fb=mjd2000fbi+nfb*(column-1);
mjd2000a=mjd2000ai+na*(plane-1);
dated=mjd20002date(mjd2000d);
datefb=mjd20002date(mjd2000fb);
datea=mjd20002date(mjd2000a);

mu_planet2=mu_J;
r_planet2=r_J;
ibody1=ibody_S;
ibody2=ibody_J;
ibody3=id_A;

t=[mjd2000d;mjd2000fb;mjd2000a];

[deltaVtotcompute] = computedeltaVtot(t,ibody1,ibody2,ibody3,mu_S,mu_planet2,r_planet2);


deltav1min=min(min(deltaVMatrix1));
[rowmind,columnmind]=find(deltaVMatrix1==deltav1min);
mjd2000dmind=mjd2000di+nd*(rowmind-1);
mjd2000fbmind=mjd2000fbi+nfb*(columnmind-1);

deltavfbmin=min(min(min(deltaVMatrixfb)));

deltav2min=min(min(deltaVMatrix2));
[rowmina,columnmina]=find(deltaVMatrix2==deltav2min);
mjd2000fbmina=mjd2000fbi+nfb*(rowmina-1);
mjd2000amina=mjd2000ai+na*(columnmina-1);

deltavfbmin=min(min(min(deltaVMatrixfb)));

% Solution: 8.3642 km/s

%deltaVdott=1.7964   km/s
%deltaVfbott=1.5314  km/s
%deltaVaott=5.0364   km/s

%          Departure: 13-2-2035
%          FLYBY:    26-6-2043
%          Arrival:   27-4-2055
%%
clear
close all
clc
%% Parameters definition

dated=[2035,2,13,12,0,0];
mjd2000d= date2mjd2000(dated);                

datefb=[2043,6,26,12,0,0];   
mjd2000fb= date2mjd2000(datefb);

datea=[2055,4,27,12,0,0];       
mjd2000a= date2mjd2000(datea);   

%Parametrs
G=astroConstants(1);
AU=astroConstants(2);


%Sun
mu_S=astroConstants(4);
r_S=astroConstants(3);
mass_S=mu_S/G;


%Saturn
ibody_S=6;
mu_Sat=astroConstants(16);
r_Sat=astroConstants(26);
[kep_Sat,~] = uplanet(mjd2000d, ibody_S);
[R_Sat, V_Sat] = kep2car (kep_Sat(1),kep_Sat(2),kep_Sat(3),kep_Sat(4), ...
    kep_Sat(5),kep_Sat(6), mu_S);
T_Sat=2*pi*sqrt(kep_Sat(1)^3/mu_S);


%Jupiter
ibody_J=5;
mu_J=astroConstants(15);
mass_J=mu_J/G;
r_J=astroConstants(25);
[kep_J,~] = uplanet(mjd2000fb, ibody_J);
[R_J, V_J] = kep2car (kep_J(1),kep_J(2),kep_J(3),kep_J(4),kep_J(5),kep_J(6), mu_S);
meanR_J=778e6;  %from Internet
T_J=2*pi*sqrt(kep_J(1)^3/mu_S);
R_JSOI=meanR_J*(mass_J/mass_S)^(2/5);
ratioRSOI_rJ=R_JSOI/r_J;

%Asteroid 1979XB
id_A=59;
r_A=r_Sat/100;    %The constant is computeted in order to have a good plot in the Transfer section
[kep_A,~,~] = ephNEO(mjd2000a,id_A);
[R_A, V_A] = kep2car (kep_A(1),kep_A(2),kep_A(3),kep_A(4),kep_A(5),kep_A(6), mu_S);
T_A=2*pi*sqrt(kep_A(1)^3/mu_S);



%% Lambert and flyby
MU=mu_S;
orbitType=0;
Nrev=0;
Ncase=0;
optionsLMR=0;
RI=R_Sat;
RFB=R_J;
RF=R_A;
TOFt1=(mjd2000fb-mjd2000d)*24*60*60;
TOFt2=(mjd2000a-mjd2000fb)*24*60*60;
[At1,Pt1,Et1,ERROR1,VI,VFB1,TPAR1,THETA1] = lambertMR(RI,RFB,TOFt1,MU,orbitType,Nrev,Ncase,optionsLMR);
[At2,Pt2,Et2,ERROR2,VFB2,VF,TPAR2,THETA2] = lambertMR(RFB,RF,TOFt2,MU,orbitType,Nrev,Ncase,optionsLMR);
VI=VI';
VFB1=VFB1';
VFB2=VFB2';
VF=VF';

deltaV1opt=VI-V_Sat;     %delta speed that we have to provide to get into the First Transfer Leg 
deltaV1nopt=norm(deltaV1opt);
deltaV2opt=V_A-VF;        %delta speed that we have to provide to get into the Asteroid Orbit
deltaV2nopt=norm(deltaV2opt);

%Fly-By
vinfa=VFB1-V_J;
vinfan=norm(vinfa);
vinfd=VFB2-V_J;
vinfdn=norm(vinfd);
delta=acos(dot(vinfa,vinfd)/(vinfan*vinfdn));
delta_minus=@(rp)   2*asin(1/(1+rp*vinfan^2/mu_J));
delta_plus=@(rp)    2*asin(1/(1+rp*vinfdn^2/mu_J));
rp_solve=@(rp)      delta_minus(rp)/2+ delta_plus(rp)/2-delta;
rp0=1*r_J;
[rp,~,flag,~] = fzero(rp_solve, rp0);
altitude=rp-r_J;


v1=sqrt(vinfan^2+2*mu_J/rp);     %velocity at the pericenter for the incoming hyperbola
e1=1+(rp*vinfan^2)/mu_J;
a1=rp/(1-e1);
p1=a1*(1-e1^2);
thetaSOI1=acos(1/e1*(p1/R_JSOI-1));
F1=2*atanh(sqrt((e1-1)/(e1+1))*tan(thetaSOI1/2));
n1=sqrt(mu_J/abs(a1)^3);
deltaT1=1/n1*(e1*sinh(F1)-F1);
deltaT1days=deltaT1/(3600*24);
v2=sqrt(vinfdn^2+2*mu_J/rp);     %velocity at the pericenter for the outgoing hyperbola
e2=1+(rp*vinfdn^2)/mu_J;
a2=rp/(1-e2);
p2=a2*(1-e2^2);
thetaSOI2=acos(1/e2*(p2/R_JSOI-1));
F2=2*atanh(sqrt((e2-1)/(e2+1))*tan(thetaSOI2/2));
n2=sqrt(mu_J/abs(a2)^3);
deltaT2=1/n2*(e2*sinh(F2)-F2);
deltaT2days=deltaT2/(3600*24);

deltaT=deltaT1+deltaT2;

deltaTdays=deltaT/(3600*24);

deltavpopt=v2-v1;    %to be provided at the pericentre
deltavfbnopt=norm(deltavpopt);

deltaVfb=vinfd-vinfa;


%Rotation of -beta1 and other steps to get a good graph in the flyby section

u=cross(vinfa,vinfd)/norm(cross(vinfa,vinfd));
DELTA1=rp*sqrt(1+(2*mu_J)/(rp*vinfan^2));
delta1=2*atan2(-a1,DELTA1);
beta1=(pi-delta1)/2;
[edir,~] = Rodrigues(vinfa,-beta1,u); %computation of the e direction

evers=edir/norm(edir);

vp1=v1;
vp2=v2;
rpvect=rp*evers;
vp1vect=vp1*cross(u,evers);
vp2vect=vp2*cross(u,evers);
centro1=(rp+abs(a1))*evers;
point1aux=centro1-44000000*vinfa;     %The constants are computeted in order to have a good plot in the flyby section(40000000 if we use -9000000 s instead of -deltaT1 in tspanpast)
centro2=(rp+abs(a2))*evers;
point2aux=centro2+2000000*vinfd;
point3aux=centro2-30*rpvect;
vettore1=[centro1 point1aux]';
vettore2=[centro2 point2aux]';
vettore3=[centro1 point3aux]';



deltaVtotopt=deltaV1nopt+deltaV2nopt+deltavfbnopt;


%Natural Flyby
%Assumption:The versor u is the same as for rp and the impact parameter. In
%this way we can have a meaningful comparison
[vinfdnat,deltaVfbnat]= Rodrigues(vinfa,delta1,u);

VFB2nat=V_J+vinfdnat;

ratio_deltaVinfd=norm(vinfdnat)/norm(vinfd);      %It's less than 1 as expected

ratio_deltaVfb=norm(deltaVfbnat)/norm(deltaVfb);  %It's less than 1 as expected

ratio_Vinfnat=norm(vinfdnat)/norm(vinfa);         %It's 1 as expected

%% Transfer and Orbit Propagation

[~, ~, inclit1, OMt1, omt1, thetat1] = car2kep (RI, VI, MU);
[~, ~, inclit2, OMt2, omt2, thetat2] = car2kep (RFB, VFB2, MU);
[a2nat, e2nat, inclit2nat, OMt2nat, omt2nat, thetat2nat] = car2kep (RFB, VFB2nat, MU);
Np=100;
tspan1=linspace( 0, TOFt1, Np);
tspan2=linspace( 0, TOFt2, Np);
tspan2nat=linspace(0,TOFt2, Np);  %What we have in the case of a natural flyby
tspanS=linspace( 0, T_Sat, Np);
tspanJ=linspace( 0, T_J, Np);
tspanA=linspace( 0, T_A, Np);
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
y0t1=[RI;VI];
y0t2=[RFB;VFB2];
y0t2nat=[RFB;VFB2nat];
y0S=[R_Sat;V_Sat];
y0J=[R_J;V_J];
y0A=[R_A;V_A];
[ ~, Yt1 ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan1, y0t1, options ); 
[ ~, Yt2 ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan2, y0t2, options );
[ ~, Yt2nat ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan2nat, y0t2nat, options );
[ ~, YS ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspanS, y0S, options ); 
[ ~, YJ ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspanJ, y0J, options ); 
[ ~, YA ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspanA, y0A, options ); 

figure(1)
plot3( YS(:,1)/AU, YS(:,2)/AU, YS(:,3)/AU, '-' );
hold on 
plot3( YJ(:,1)/AU, YJ(:,2)/AU, YJ(:,3)/AU, '-' );
hold on 
plot3( YA(:,1)/AU, YA(:,2)/AU, YA(:,3)/AU, '-' );
hold on 
plot3( Yt1(:,1)/AU, Yt1(:,2)/AU, Yt1(:,3)/AU, '-' );
hold on 
plot3( Yt2(:,1)/AU, Yt2(:,2)/AU, Yt2(:,3)/AU, '-' );
hold on 
plot3( Yt2nat(:,1)/AU, Yt2nat(:,2)/AU, Yt2nat(:,3)/AU, '-' );
S.foto='Sun.jpg';
PlanetSphere(100*r_S/AU,0,0,0,S);
Sat.foto='Saturn.jpg';
PlanetSphere(1000*r_Sat/AU,R_Sat(1)/AU,R_Sat(2)/AU,R_Sat(3)/AU,Sat);
J.foto='Jupiter.jpg';
PlanetSphere(1000*r_J/AU,Yt2(1,1)/AU,Yt2(1,2)/AU,Yt2(1,3)/AU,J);
A.foto='Asteroid.jpg';
PlanetSphere(20000*r_A/AU,Yt2(end,1)/AU,Yt2(end,2)/AU,Yt2(end,3)/AU,A);
legend('$Saturn$','$Jupiter$','$Asteroid$','$First$ $Leg$','$Second$ $Leg$', ...
    '$Second$ $Leg$ $with$ $Natural$ $Flyby$','Interpreter','latex');
title('$Interplanetary$ $Mission$','Interpreter','latex');
xlabel('$x$  $[AU]$','Interpreter','latex');
ylabel('$y$  $[AU]$','Interpreter','latex');
zlabel('$z$  $[AU]$','Interpreter','latex');
