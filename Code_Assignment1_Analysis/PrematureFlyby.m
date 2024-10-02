clear;
close all;
clc;
%
%
% Script for data in the report
%
%PrematureFlyby:
%          Departure: [10227  13165]     
%          FLYBY:    13165+-365         17-1-2036
%          Arrival:   14045+-365         15-6-2038


nd=30;
md=98;
nfb=8;
mfb=98;
na=8;
ma=98;
   
mjd2000di= 10227;
mjd2000df= mjd2000di+md*nd;

mjd2000fbi= 13165-365;
mjd2000fbf= mjd2000fbi+mfb*nfb;

mjd2000ai= 14045-365;       
mjd2000af= mjd2000ai+ma*na;


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
T_Satdays=T_Sat/(3600*24);


%Jupiter
ibody_J=5;
mu_J=astroConstants(15);
r_J=astroConstants(25);
[kep_J,~] = uplanet(mjd2000fbi, ibody_J);
T_J=2*pi*sqrt(kep_J(1)^3/mu_S);
T_Jdays=T_J/(3600*24);

[R_J, ~] = kep2car (kep_J(1),kep_J(2),kep_J(3),kep_J(4),kep_J(5),kep_J(6), mu_S);



%Asteroid 1979XB
id_A=59;
r_A=r_Sat/100;    %così si vede
[kep_A,~,~] = ephNEO(mjd2000ai,id_A);
T_A=2*pi*sqrt(kep_A(1)^3/mu_S);
T_Adays=T_A/(3600*24);

[R_A, ~] = kep2car (kep_A(1),kep_A(2),kep_A(3),kep_A(4),kep_A(5),kep_A(6), mu_S);


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



%% Lambert

MU=mu_S;
orbitType=0;
Nrev=0;
Ncase=0;
optionsLMR=0;


vettdep=mjd2000di:nd:mjd2000df;
vettfb=mjd2000fbi:nfb:mjd2000fbf;
vettarr=mjd2000ai:na:mjd2000af;



%%
y=1;
u=1;
p=1;

TOF1Matrix=NaN(md+1,mfb+1);
TOF2Matrix=NaN(mfb+1,ma+1);
deltaVMatrix1=NaN(md+1,mfb+1);
deltaVMatrix2=NaN(mfb+1,ma+1);
deltaVMatrixfb=NaN(md+1,mfb+1,ma+1);
deltaVTOT=NaN(md+1,mfb+1,ma+1);

TPAR1Matrix=NaN(md+1,mfb+1);

% TMIN=NaN(m+1,m+1);
FLAG=NaN(md+1,mfb+1,ma+1);


for q=mjd2000di:nd:mjd2000df
    
    [kep_Sat,~] = uplanet(q, ibody_S);  
    [R_Sat, V_Sat] = kep2car (kep_Sat(1),kep_Sat(2),kep_Sat(3),kep_Sat(4), ...
    kep_Sat(5),kep_Sat(6), mu_S);  
    RI=R_Sat;




    for w=mjd2000fbi:nfb:mjd2000fbf
  
        
        TOFt1=(w-q)*24*60*60;
        
        if (TOFt1> TPAR1  && TOFt1<T_Sat)                         
            [kep_J,~] = uplanet(w, ibody_J);
            [R_J, V_J] = kep2car (kep_J(1),kep_J(2),kep_J(3),kep_J(4),kep_J(5),kep_J(6), mu_S);
            RFB=R_J;
            [~,~,~,~,VI,VFB1,TPAR1,~] = lambertMR(RI,RFB,TOFt1,MU,orbitType,Nrev,Ncase,optionsLMR);
            TPAR1Matrix(y,u)=TPAR1;
            VI=VI';
            VFB1=VFB1';
            deltaV1=VI-V_Sat;      
            deltaV1n=norm(deltaV1);
            deltaVMatrix1(y,u)=deltaV1n;
            TOF1Matrix(y,u)=TOFt1;
           
                for t=mjd2000ai:na:mjd2000af
                    
                    TOFt2=(t-w)*24*60*60;
                    
                    if (TOFt2>TPAR2 && TOFt2<=1.1*T_J) 
                        [kep_A,~,~] = ephNEO(t,id_A);
                        [R_A, V_A] = kep2car (kep_A(1),kep_A(2),kep_A(3),kep_A(4),kep_A(5),kep_A(6), mu_S);
                        
                        RF=R_A;
                        [~,~,~,~,VFB2,VF,~,~] = lambertMR(RFB,RF,TOFt2,MU,orbitType,Nrev,Ncase,optionsLMR);
                        
                            VFB2=VFB2';
                            VF=VF';
                            deltaV2=V_A-VF;        
                            deltaV2n=norm(deltaV2);
                            TOF2Matrix(u,p)=TOFt2;
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
        u=u+1;

    end
    u=1;
    y=y+1;
end


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

deltav2min=min(min(deltaVMatrix2));
[rowmina,columnmina]=find(deltaVMatrix2==deltav2min);
mjd2000fbmina=mjd2000fbi+nfb*(rowmina-1);
mjd2000amina=mjd2000ai+na*(columnmina-1);

deltavfbmin=min(min(min(deltaVMatrixfb)));

%% Contour

figure(1)
contourf(vettdep,vettfb,deltaVMatrix1',linspace(deltav1min,deltav1min+5,20),'ShowText','on');
xlabel('Departure');
ylabel('flyby');
title('T1');

figure(2)
contourf(vettfb,vettarr,deltaVMatrix2',linspace(deltav2min,deltav2min+5,20),'ShowText','on'); 
xlabel('flyby');
ylabel('Arrival');
title('T2');


figure(3)
contour(vettdep,vettfb,deltaVMatrix1',linspace(deltav1min,deltav1min+5,20),'ShowText','on');
c=colorbar;
c.Label.String= 'deltaV [km/s]';
xlabel('Departure');
ylabel('flyby');
title('T1f');
hold on
plot(mjd2000d,mjd2000fb,'*');


figure(5)
contour(vettfb,vettarr,deltaVMatrix2',linspace(deltaVaott,deltaVaott+5,20),'ShowText','on');
c=colorbar;
c.Label.String= 'deltaV [km/s]';
xlabel('flyby');
ylabel('Arrival');
title('T2f with optimum');
hold on
plot(mjd2000fb,mjd2000a,'*');
hold on
plot(mjd2000fbmina,mjd2000amina,'o');

figure(6)
contour(vettfb,vettarr,deltaVMatrix2',linspace(deltav2min,deltav2min+5,20),'ShowText','on');
c=colorbar;
c.Label.String= 'deltaV [km/s]';
xlabel('flyby');
ylabel('Arrival');
title('T2f with minimum of 2');
hold on
plot(mjd2000fb,mjd2000a,'*');
hold on
plot(mjd2000fbmina,mjd2000amina,'o');

%%  The best solution for this new region is:
clear
close all
clc

%PrematureFlyby:
%          Departure: 14-06-2031    11487
%          FLYBY:    17-4-2036     13256
%          Arrival:   15-4-2038     13984

% Daily discretization in the neighbourhood of +-45 from departure
% and +- 12 from fly by

nd=1;
md=90;
nfb=1;
mfb=24;
na=1;
ma=24;



mjd2000di= 11487-45;
mjd2000df= mjd2000di+md*nd;

mjd2000fbi= 13256-12;
mjd2000fbf= mjd2000fbi+mfb*nfb;

mjd2000ai= 13984-12;       
mjd2000af= mjd2000ai+ma*na;


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
T_Satdays=T_Sat/(3600*24);


%Jupiter
ibody_J=5;
mu_J=astroConstants(15);
r_J=astroConstants(25);
[kep_J,~] = uplanet(mjd2000fbi, ibody_J);
T_J=2*pi*sqrt(kep_J(1)^3/mu_S);
T_Jdays=T_J/(3600*24);

[R_J, ~] = kep2car (kep_J(1),kep_J(2),kep_J(3),kep_J(4),kep_J(5),kep_J(6), mu_S);



%Asteroid 1979XB
id_A=59;
r_A=r_Sat/100;    %così si vede
[kep_A,~,~] = ephNEO(mjd2000ai,id_A);
T_A=2*pi*sqrt(kep_A(1)^3/mu_S);
T_Adays=T_A/(3600*24);

[R_A, ~] = kep2car (kep_A(1),kep_A(2),kep_A(3),kep_A(4),kep_A(5),kep_A(6), mu_S);

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



%% Lambert

MU=mu_S;
orbitType=0;
Nrev=0;
Ncase=0;
optionsLMR=0;


vettdep=mjd2000di:nd:mjd2000df;
vettfb=mjd2000fbi:nfb:mjd2000fbf;
vettarr=mjd2000ai:na:mjd2000af;

y=1;
u=1;
p=1;

TOF1Matrix=NaN(md+1,mfb+1);
TOF2Matrix=NaN(mfb+1,ma+1);
deltaVMatrix1=NaN(md+1,mfb+1);
deltaVMatrix2=NaN(mfb+1,ma+1);
deltaVMatrixfb=NaN(md+1,mfb+1,ma+1);
deltaVTOT=NaN(md+1,mfb+1,ma+1);
e1matrix=NaN(md+1,mfb+1);
TPAR1Matrix=NaN(md+1,mfb+1);
TPAR2Matrix=NaN(mfb+1,ma+1);
% TMIN=NaN(m+1,m+1);
FLAG=NaN(md+1,mfb+1,ma+1);


for q=mjd2000di:nd:mjd2000df
    
    [kep_Sat,~] = uplanet(q, ibody_S);   %Kep è riga
    [R_Sat, V_Sat] = kep2car (kep_Sat(1),kep_Sat(2),kep_Sat(3),kep_Sat(4), ...
    kep_Sat(5),kep_Sat(6), mu_S);  %R è colonna
    RI=R_Sat;




    for w=mjd2000fbi:nfb:mjd2000fbf
  
        
        TOFt1=(w-q)*24*60*60;
        
        if (TOFt1> TPAR1  && TOFt1<T_Sat)                      
            [kep_J,~] = uplanet(w, ibody_J);
            [R_J, V_J] = kep2car (kep_J(1),kep_J(2),kep_J(3),kep_J(4),kep_J(5),kep_J(6), mu_S);
            RFB=R_J;
            [~,~,~,~,VI,VFB1,TPAR1,~] = lambertMR(RI,RFB,TOFt1,MU,orbitType,Nrev,Ncase,optionsLMR);
            TPAR1Matrix(y,u)=TPAR1;
            VI=VI';
            VFB1=VFB1';
            deltaV1=VI-V_Sat;     
            deltaV1n=norm(deltaV1);
            
            deltaVMatrix1(y,u)=deltaV1n;
            TOF1Matrix(y,u)=TOFt1;
            
                for t=mjd2000ai:na:mjd2000af
                    
                    TOFt2=(t-w)*24*60*60;
                    
                    if (TOFt2>TPAR2 && TOFt2<=1.1*T_J) 
                        [kep_A,~,~] = ephNEO(t,id_A);
                        [R_A, V_A] = kep2car (kep_A(1),kep_A(2),kep_A(3),kep_A(4),kep_A(5),kep_A(6), mu_S);
                        
                        RF=R_A;
                        [~,~,~,~,VFB2,VF,~,~] = lambertMR(RFB,RF,TOFt2,MU,orbitType,Nrev,Ncase,optionsLMR);
                        
                            VFB2=VFB2';
                            VF=VF';
                            deltaV2=V_A-VF;        
                            deltaV2n=norm(deltaV2);
                            TOF2Matrix(u,p)=TOFt2;
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
        u=u+1;

    end
    u=1;
    y=y+1;
end

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

deltav2min=min(min(deltaVMatrix2));
[rowmina,columnmina]=find(deltaVMatrix2==deltav2min);
mjd2000fbmina=mjd2000fbi+nfb*(rowmina-1);
mjd2000amina=mjd2000ai+na*(columnmina-1);

deltavfbmin=min(min(min(deltaVMatrixfb)));

%% Contour

figure(1)
contourf(vettdep,vettfb,deltaVMatrix1',linspace(deltav1min,deltav1min+5,20),'ShowText','on');
xlabel('Departure');
ylabel('flyby');
title('T1');

figure(2)
contourf(vettfb,vettarr,deltaVMatrix2',linspace(deltav2min,deltav2min+5,20),'ShowText','on'); 
xlabel('flyby');
ylabel('Arrival');
title('T2');


figure(3)
contour(vettdep,vettfb,deltaVMatrix1',linspace(deltav1min,deltav1min+5,20),'ShowText','on');
c=colorbar;
c.Label.String= 'deltaV [km/s]';
xlabel('Departure');
ylabel('flyby');
title('T1f');
hold on
plot(mjd2000d,mjd2000fb,'*');


figure(5)
contour(vettfb,vettarr,deltaVMatrix2',linspace(deltaVaott,deltaVaott+5,20),'ShowText','on');
c=colorbar;
c.Label.String= 'deltaV [km/s]';
xlabel('flyby');
ylabel('Arrival');
title('T2f with optimum');
hold on
plot(mjd2000fb,mjd2000a,'*');
hold on
plot(mjd2000fbmina,mjd2000amina,'o');

figure(6)
contour(vettfb,vettarr,deltaVMatrix2',linspace(deltav2min,deltav2min+5,20),'ShowText','on');
c=colorbar;
c.Label.String= 'deltaV [km/s]';
xlabel('flyby');
ylabel('Arrival');
title('T2f showing the minimum of 2');
hold on
plot(mjd2000fb,mjd2000a,'*');
hold on
plot(mjd2000fbmina,mjd2000amina,'o');

%% The solution for pre-mature fly by is:


%PrematureFlyby
%          Departure: 17-06-2031 alle 12   
%          FLYBY:    13-4-2036  alle 12   
%          Arrival:   20-4-2038  alle 12   


clear
close all
clc
dated=[2031,6,17,12,0,0];
mjd2000d= date2mjd2000(dated);                

datefb=[2036,4,13,12,0,0];   
mjd2000fb= date2mjd2000(datefb);

datea=[2038,4,20,12,0,0];       
mjd2000a= date2mjd2000(datea);   

%Sun
mu_S=astroConstants(4);
r_S=astroConstants(3);


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
r_J=astroConstants(25);
[kep_J,~] = uplanet(mjd2000fb, ibody_J);
[R_J, V_J] = kep2car (kep_J(1),kep_J(2),kep_J(3),kep_J(4),kep_J(5),kep_J(6), mu_S);
T_J=2*pi*sqrt(kep_J(1)^3/mu_S);


%Asteroid 1979XB
id_A=59;
r_A=r_Sat/100;    %così si vede
[kep_A,~,~] = ephNEO(mjd2000a,id_A);
[R_A, V_A] = kep2car (kep_A(1),kep_A(2),kep_A(3),kep_A(4),kep_A(5),kep_A(6), mu_S);
T_A=2*pi*sqrt(kep_A(1)^3/mu_S);



%% Lambert
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

deltaV1opt=VI-V_Sat;     
deltaV1nopt=norm(deltaV1opt);
deltaV2opt=V_A-VF;       
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

v1=sqrt(vinfan^2+2*mu_J/rp);    
e1=1+(rp*vinfan^2)/mu_J;
a1=rp/(1-e1);
v2=sqrt(vinfdn^2+2*mu_J/rp);    
e2=1+(rp*vinfdn^2)/mu_J;
a2=rp/(1-e2);

deltavpopt=v2-v1;    
deltavfbnopt=norm(deltavpopt);

% Rotation of -beta1
u=cross(vinfa,vinfd)/norm(cross(vinfa,vinfd));
DELTA1=rp*sqrt(1+(2*mu_J)/(rp*vinfan^2));
delta1=2*atan2(-a1,DELTA1);
beta1=(pi-delta1)/2;
[edir,~] = Rodrigues(vinfa,-beta1,u); % computation of e

evers=edir/norm(edir);

vp1=v1;
vp2=v2;
rpvect=rp*evers;
vp1vect=vp1*cross(u,evers);
vp2vect=vp2*cross(u,evers);
centro1=(rp+abs(a1))*evers;
point1aux=centro1-40000000*vinfa;
centro2=(rp+abs(a2))*evers;
point2aux=centro2+2000000*vinfd;
point3aux=centro2-30*rpvect;
vettore1=[centro1 point1aux]';
vettore2=[centro2 point2aux]';
vettore3=[centro1 point3aux]';



deltaVtotopt=deltaV1nopt+deltaV2nopt+deltavfbnopt;


%% Transfer and Orbit Propagation
[~, ~, inclit1, OMt1, omt1, thetat1] = car2kep (RI, VI, MU);
[~, ~, inclit2, OMt2, omt2, thetat2] = car2kep (RFB, VFB2, MU);
Np=10000;
tspan1=linspace( 0, TOFt1, Np);
tspan2=linspace( 0, TOFt2, Np);
tspanS=linspace( 0, T_Sat, Np);
tspanJ=linspace( 0, T_J, Np);
tspanA=linspace( 0, T_A, Np);
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
y0t1=[RI;VI];
y0t2=[RFB;VFB2];
y0S=[R_Sat;V_Sat];
y0J=[R_J;V_J];
y0A=[R_A;V_A];
[ ~, Yt1 ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan1, y0t1, options ); 
[ ~, Yt2 ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan2, y0t2, options );
[ ~, YS ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspanS, y0S, options ); 
[ ~, YJ ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspanJ, y0J, options ); 
[ ~, YA ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspanA, y0A, options ); 

figure(5)
plot3( YS(:,1), YS(:,2), YS(:,3), '-' );
hold on 
plot3( YJ(:,1), YJ(:,2), YJ(:,3), '-' );
hold on 
plot3( YA(:,1), YA(:,2), YA(:,3), '-' );
hold on 
plot3( Yt1(:,1), Yt1(:,2), Yt1(:,3), '-' );
hold on 
plot3( Yt2(:,1), Yt2(:,2), Yt2(:,3), '-' );
S.foto='Sun.jpg';
PlanetSphere(10*r_S,0,0,0,S);
Sat.foto='Saturn.jpg';
PlanetSphere(1000*r_Sat,R_Sat(1),R_Sat(2),R_Sat(3),Sat);
J.foto='Jupiter.jpg';
PlanetSphere(1000*r_J,Yt2(1,1),Yt2(1,2),Yt2(1,3),J);
A.foto='Asteroid.jpg';
PlanetSphere(20000*r_A,Yt2(end,1),Yt2(end,2),Yt2(end,3),A);

legend('Saturn','Jupiter','Asteroid','T1','T2');
title('Interplanetary Mission');

