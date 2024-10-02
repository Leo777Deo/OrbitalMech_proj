%% ----------------------Starbusters_Group_2341------------------
%
% Barbiera Andrea
%
% De Luca Leo
%
% Perusini Gianluca
%
% Poverini Viola
%
%% ------------------ Assignment 1: Interplanetary Mission-------------------
%
% The final solution is analysed in this script.
%
% The evaluation of the different cases that are not the optimal one 
% are analysed in the folder Code_Assignment1_Analysis.
%
% In the following script the final optimal solution is presented.
%


clear;
close all;
clc;


%% Optimization of the second leg of the interplanetary mission


%% Parameters definition

n=61; % dummy variables to allocate the right space for vectors
m=120;

datedi=[2028,1,1,0,0,0];

datefbi=[2036,1,1,0,0,0];           % We use a start time of about 8 years 
%                                     (two years less than a TOFH1) subsequent to datedi     

dateai=[2038,1,1,0,0,0];            % We have a window that starts about 10 years 
%                                     after launch window start(see datedf)                              

mjd2000di= date2mjd2000(datedi);
mjd2000df=mjd2000di+m*n;            % We use an end time of departures of about 10 years 
%                                     before the end of the window (a little less than TOFH1+TPAR2 
%                                     because then we are comfortable on first section which is almost planar
%                                     planar and circular while the second section I guard with TPAR2)


mjd2000fbi= date2mjd2000(datefbi);
mjd2000fbf=mjd2000fbi+m*n;          % We can use an end flyby time of about 2 years (TPAR2 is about 1.2 years
%                                    but we know that around TPAR we don't have cheap solutions)  before the 
%                                    end of arrivals

datefbf=mjd20002date(mjd2000fbf);        
mjd2000ai= date2mjd2000(dateai);        
mjd2000af= mjd2000ai+m*n;
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
atmosphheight=5000;   %[km]


%Asteroid 1979XB
id_A=59;
[kep_A,~,~] = ephNEO(mjd2000ai,id_A);
T_A=2*pi*sqrt(kep_A(1)^3/mu_S);
[R_A, ~] = kep2car (kep_A(1),kep_A(2),kep_A(3),kep_A(4),kep_A(5),kep_A(6), mu_S);

%SYNODIC PERIOD between Jupiter and Asteroid
TsynJA=(T_A*T_J)/abs(T_A-T_J);          %it's about 5 years
%SYNODIC PERIOD between Jupiter and Saturn
TsynSJ=(T_Sat*T_J)/abs(T_Sat-T_J);      %it's about 20 years


% By calculating TPAR2,we can use conditions on the various TOFs to not  
% calculate the values when we do not meet the requirements

orbitType=0;
Nrev=0;
Ncase=0;
optionsLMR=0;


TOF2=1;
RI=R_J;
RF=R_A;
[~,~,~,~,~,~,TPAR2,~] = lambertMR(RI,RF,TOF2,mu_S,orbitType,Nrev,Ncase,optionsLMR);

a1=(norm(R_Sat)+norm(R_J))/2;
TOFH1=pi*sqrt(a1^3/mu_S);  %about 10 years


a2=(norm(R_J)+norm(R_A))/2;
TOFH2=pi*sqrt(a2^3/mu_S);   %about 3 years  

%% Lambert
MU=mu_S;
orbitType=0;
Nrev=0;
Ncase=0;
optionsLMR=0;

vettfb=mjd2000fbi:n:mjd2000fbf;
vettarr=mjd2000ai:n:mjd2000af;
y=1;
p=1;
u=1;
deltaVMatrix2=NaN(m+1,m+1);
TOF2Matrix=NaN(m+1,m+1);

for w=mjd2000fbi:n:mjd2000fbf
            [kep_J,~] = uplanet(w, ibody_J);
            [R_J, ~] = kep2car (kep_J(1),kep_J(2),kep_J(3),kep_J(4),kep_J(5),kep_J(6), mu_S);
            RFB=R_J;
    for t=mjd2000ai:n:mjd2000af
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
                        deltaVMatrix2(u,p)=deltaV2n;   
                        TOF2Matrix(u,p)=TOFt2;
                     end
     p=p+1;                   
    end
    p=1;
    u=u+1;
end
deltav2min=min(min(deltaVMatrix2));


%% Contour Plot

figure(1)
startDate1 = datenum('01-01-2036', 'dd mmm yyyy');
startDate2 = datenum('16-01-2056', 'dd mmm yyyy'); %end of flyby
endDate1 = datenum('01-01-2038',  'dd mmm yyyy');
endDate2 = datenum('16-01-2058',  'dd mmm yyyy');
xData = linspace(startDate1, startDate2, length(vettfb));
yData = linspace(endDate1, endDate2, length(vettarr));
[T1, T2] = meshgrid(xData, yData);

contour(T1,T2,deltaVMatrix2',linspace(deltav2min,deltav2min+6,20),'ShowText','off');
c=colorbar;
c.Label.String= "\Deltav_{arrival} [km/s]";     
ha=gca;
clim(ha,deltav2min+[0 5.5]);
clim(ha,'manual');
hold on
contour(T1,T2,TOF2Matrix',800:500:6000,'k','ShowText','on');
title('$Arrival$ $Cost$ $Plot$','Interpreter','latex')
xlabel('$Flyby$ $Date$ $of$ $Jupiter$','Interpreter','latex')
ylabel('$Arrival$ $Date$ $to$ $the$ $Asteroid$','Interpreter','latex')
datetick('x', 'dd-mmm-yyyy', 'keeplimits', 'keepticks');
datetick('y', 'dd-mmm-yyyy', 'keeplimits', 'keepticks');
grid on

%% We can see from the Contour plot which regions are the best: 

dateIFB=mjd20002date(1.54e4);
dateFFB=mjd20002date(1.73e4);
dateAI=mjd20002date(1.7e4);
dateAF=mjd20002date(2.05e4);

%% Optimization of the interplanetary mission

nd=60; % dummy variables
md=55;
nfb=60;
mfb=37;
na=60;
ma=61;

datedi=[2031,1,1,0,0,0];
datefbi=[2042,1,1,0,0,0];     
dateai=[2046,1,1,0,0,0];       
mjd2000di= date2mjd2000(datedi);
mjd2000df=mjd2000di+md*nd;                
datedf=mjd20002date(mjd2000df);

mjd2000fbi= date2mjd2000(datefbi);
mjd2000fbf=mjd2000fbi+mfb*nfb;               
datefbf=mjd20002date(mjd2000fbf);

mjd2000ai= date2mjd2000(dateai);        
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
[kep_A,~,~] = ephNEO(mjd2000ai,id_A);
T_A=2*pi*sqrt(kep_A(1)^3/mu_S);
[R_A, ~] = kep2car (kep_A(1),kep_A(2),kep_A(3),kep_A(4),kep_A(5),kep_A(6), mu_S);

%SYNODIC PERIOD between Jupiter and Asteroid
TsynJA=(T_A*T_J)/abs(T_A-T_J);          %it's about 5 years
%SYNODIC PERIOD between Jupiter and Saturn
TsynSJ=(T_Sat*T_J)/abs(T_Sat-T_J);      %it's about 20 years


% By calculating TPAR2, we can use the conditions on the various TOFs 
% to not  calculate the values when we do not meet the requirements

orbitType=0;
Nrev=0;
Ncase=0;
optionsLMR=0;


TOF2=1;
RI=R_J;
RF=R_A;
[~,~,~,~,~,~,TPAR2,~] = lambertMR(RI,RF,TOF2,mu_S,orbitType,Nrev,Ncase,optionsLMR);

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
    
    [kep_Sat,~] = uplanet(q, ibody_S);   % Kep is row
    [R_Sat, V_Sat] = kep2car (kep_Sat(1),kep_Sat(2),kep_Sat(3),kep_Sat(4), ...
    kep_Sat(5),kep_Sat(6), mu_S);  % R is column
    RI=R_Sat;




    for w=mjd2000fbi:nfb:mjd2000fbf
  
        
        TOFt1=(w-q)*24*60*60;
        
        if (TOFt1> 2920*(3600*24)  && TOFt1<=4015*(3600*24))                             %(TOFt1> 2920*(3600*24)  && TOFt1<=4015*(3600*24))   %TOFH1+-1.5 Years                                             
            [kep_J,~] = uplanet(w, ibody_J);
            [R_J, V_J] = kep2car (kep_J(1),kep_J(2),kep_J(3),kep_J(4),kep_J(5),kep_J(6), mu_S);
            RFB=R_J;
            [~,~,~,~,VI,VFB1,~,~] = lambertMR(RI,RFB,TOFt1,MU,orbitType,Nrev,Ncase,optionsLMR);
            VI=VI';
            VFB1=VFB1';
            deltaV1=VI-V_Sat;      
            deltaV1n=norm(deltaV1);
            deltaVMatrix1(y,u)=deltaV1n;
            TOF1Matrix(y,u)=TOFt1/(3600*24);    % [days]

                for t=mjd2000ai:na:mjd2000af
                    
                    TOFt2=(t-w)*24*60*60;
                    
                    if (TOFt2>TPAR2 && TOFt2<=1.1*T_J)  % minimum is the parabolic and the maximum is through energy reasoning
                        [kep_A,~,~] = ephNEO(t,id_A);
                        [R_A, V_A] = kep2car (kep_A(1),kep_A(2),kep_A(3),kep_A(4),kep_A(5),kep_A(6), mu_S);
                        
                        RF=R_A;
                        [~,~,~,~,VFB2,VF,~,~] = lambertMR(RFB,RF,TOFt2,MU,orbitType,Nrev,Ncase,optionsLMR);
                        
                            VFB2=VFB2';
                            VF=VF';
                            deltaV2=V_A-VF;        
                            deltaV2n=norm(deltaV2);
                            TOF2Matrix(u,p)=TOFt2/(3600*24);  % [days]
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
                            if (altitude>atmosphheight && flag==1)        % atmosphere condition and condition on FLAG 
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

%% Analysis of the results:

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



deltav2min=min(min(deltaVMatrix2));
[rowmina,columnmina]=find(deltaVMatrix2==deltav2min);
mjd2000fbmina=mjd2000fbi+nfb*(rowmina-1);
mjd2000amina=mjd2000ai+na*(columnmina-1);

deltavfbmin=min(min(min(deltaVMatrixfb)));

%% Contour Plots

figure(2)
startDate1 = datenum('01-01-2042', 'dd mmm yyyy');
startDate2 = datenum('30-01-2048', 'dd mmm yyyy'); %fine flyby
endDate1 = datenum('01-01-2046',  'dd mmm yyyy');
endDate2 = datenum('09-01-2056',  'dd mmm yyyy');
xData = linspace(startDate1, startDate2, length(vettfb));
yData = linspace(endDate1, endDate2, length(vettarr));
[T1, T2] = meshgrid(xData, yData);

contour(T1,T2,deltaVMatrix2',linspace(deltav2min,deltav2min+6,20),'ShowText','off');
c=colorbar;
c.Label.String= "\Deltav_{arrival} [km/s]";     
ha=gca;
clim(ha,deltav2min+[0 5.5]);
clim(ha,'manual');
hold on
contour(T1,T2,TOF2Matrix',800:500:6000,'k','ShowText','on');
title('$Arrival$ $Cost$ $Plot$','Interpreter','latex')
xlabel('$Flyby$ $Date$ $of$ $Jupiter$','Interpreter','latex')
ylabel('$Arrival$ $Date$ $to$ $the$ $Asteroid$','Interpreter','latex')
datetick('x', 'dd-mmm-yyyy', 'keeplimits', 'keepticks');
datetick('y', 'dd-mmm-yyyy', 'keeplimits', 'keepticks');
grid on


% 1.1 The one presented as optimal is:
%          DEPARTURE: 18-8-2033
%          FLYBY:     25-6-2043
%          ARRIVAL:   24-08-2047
%
%
%
% 2.1 The one presented with smaller deltaV on arrival:
%          DEPARTURE: 
%          FLYBY:     15941   [mjd2000]
%          ARRIVAL:   18842   [mjd2000]
%
%
%
% 3.1 Potentially interesting region:
%          DEPARTURE: 
%          FLYBY:     15942   [mjd2000]
%          ARRIVAL:   20221   [mjd2000]
%
%
% 1.3 Potentially interesting region:
%          DEPARTURE: 
%          FLYBY:     17143   [mjd2000]
%          ARRIVAL:   17994   [mjd2000]
%
% The cases 1.3, 2.1 and 3.1 were analysed in the Code_Assignment1_Analysis
% 

%% 1.1
% We take +-2month wide windows around the optimum.
% Now the discretization is 2 days

%1.1 The one presented as optimal is:
%          DEPARTURE: 18-8-2033
%          FLYBY:     25-6-2043
%          ARRIVAL:   24-08-2047

nd=2; % dummy variables to allocate right space for vectors
md=60;
nfb=2;
mfb=60;
na=2;
ma=60;

datedi=[2033,6,18,0,0,0];
datefbi=[2043,4,25,0,0,0];     
dateai=[2047,6,24,0,0,0];       
mjd2000di= date2mjd2000(datedi);
mjd2000df=mjd2000di+md*nd;                
datedf=mjd20002date(mjd2000df);

mjd2000fbi= date2mjd2000(datefbi);
mjd2000fbf=mjd2000fbi+mfb*nfb;               
datefbf=mjd20002date(mjd2000fbf);

mjd2000ai= date2mjd2000(dateai);        
mjd2000af= mjd2000ai+ma*na;
dateaf=mjd20002date(mjd2000af);

%Sun
mu_S=astroConstants(4);
r_S=astroConstants(3);

%Saturn
ibody_S=6;
r_Sat=astroConstants(26);


%Jupiter
ibody_J=5;
mu_J=astroConstants(15);
r_J=astroConstants(25);
[kep_J,~] = uplanet(mjd2000fbi, ibody_J);
T_J=2*pi*sqrt(kep_J(1)^3/mu_S);
[R_J, ~] = kep2car (kep_J(1),kep_J(2),kep_J(3),kep_J(4),kep_J(5),kep_J(6), mu_S);


%Asteroid 1979XB
id_A=59;
[kep_A,~,~] = ephNEO(mjd2000ai,id_A);
[R_A, ~] = kep2car (kep_A(1),kep_A(2),kep_A(3),kep_A(4),kep_A(5),kep_A(6), mu_S);


TOF2=1;
RI=R_J;
RF=R_A;
[~,~,~,~,~,~,TPAR2,~] = lambertMR(RI,RF,TOF2,mu_S,orbitType,Nrev,Ncase,optionsLMR);


%% Lambert
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
    
    [kep_Sat,~] = uplanet(q, ibody_S);   %Kep is row
    [R_Sat, V_Sat] = kep2car (kep_Sat(1),kep_Sat(2),kep_Sat(3),kep_Sat(4), ...
    kep_Sat(5),kep_Sat(6), mu_S);  %R is column
    RI=R_Sat;




    for w=mjd2000fbi:nfb:mjd2000fbf
  
        
        TOFt1=(w-q)*24*60*60;
        
        if (TOFt1> 2920*(3600*24)  && TOFt1<=4015*(3600*24))                             %(TOFt1> 2920*(3600*24)  && TOFt1<=4015*(3600*24))   %TOFH1+-1.5 ANNI                                                
            [kep_J,~] = uplanet(w, ibody_J);
            [R_J, V_J] = kep2car (kep_J(1),kep_J(2),kep_J(3),kep_J(4),kep_J(5),kep_J(6), mu_S);
            RFB=R_J;
            [~,~,~,~,VI,VFB1,~,~] = lambertMR(RI,RFB,TOFt1,MU,orbitType,Nrev,Ncase,optionsLMR);
            VI=VI';
            VFB1=VFB1';
            deltaV1=VI-V_Sat;     
            deltaV1n=norm(deltaV1);
            deltaVMatrix1(y,u)=deltaV1n;
            TOF1Matrix(y,u)=TOFt1/(3600*24);         %[days]

                for t=mjd2000ai:na:mjd2000af
                    
                    TOFt2=(t-w)*24*60*60;
                    
                    if (TOFt2>TPAR2 && TOFt2<=1.1*T_J)  % minimum is the parabolic and the maximum is through energy reasoning
                        [kep_A,~,~] = ephNEO(t,id_A);
                        [R_A, V_A] = kep2car (kep_A(1),kep_A(2),kep_A(3),kep_A(4),kep_A(5),kep_A(6), mu_S);
                        
                        RF=R_A;
                        [~,~,~,~,VFB2,VF,~,~] = lambertMR(RFB,RF,TOFt2,MU,orbitType,Nrev,Ncase,optionsLMR);
                        
                            VFB2=VFB2';
                            VF=VF';
                            deltaV2=V_A-VF;        
                            deltaV2n=norm(deltaV2);
                            TOF2Matrix(u,p)=TOFt2/(3600*24);   %[days]
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
                            if (altitude>atmosphheight && flag==1)   
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

%% Analysis of the results

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


%The total solution for deltaV is 7.8658 km/s    

%deltaVdott=1.9382 km/s 
%deltaVfbott=0.5875 km/s 
%deltaVaott=5.34 km/s 

%          DEPARTURE: 7-8-2033
%          FLYBY:     30-7-2043
%          ARRIVAL:   1-8-2047


%% Use of fminunc 


fun=@(t) computedeltaVtot(t,ibody1,ibody2,ibody3,mu_S,mu_planet2,r_planet2);
t0 = [mjd2000d mjd2000fb mjd2000a];
[T] = fminunc(@(t) fun(t),t0);
dateDfminunc=mjd20002date(T(1));
dateFBfminunc=mjd20002date(T(2));
dateAfminunc=mjd20002date(T(3));
deltaVtotfminunc=computedeltaVtot(T,ibody1,ibody2,ibody3,mu_S,mu_planet2,r_planet2);

  
%          DEPARTURE: 8-8-2033 AT 3:35
%          FLYBY:     29-7-2043 AT 6:11
%          ARRIVAL:   31-7-2047  AT 00:04

%We stop at a choice with minute discretization 

%% -------------INTERPLANETARY MISSION--------------

%% Parameters definition

dated=[2033,8,8,3,35,0];
mjd2000d= date2mjd2000(dated);                

datefb=[2043,7,29,6,11,0];   
mjd2000fb= date2mjd2000(datefb);

datea=[2047,7,31,0,4,0];       
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
meanR_J=778e6;   %Mean Heliocentric Radius of Jupiter [km]               
T_J=2*pi*sqrt(kep_J(1)^3/mu_S);
R_JSOI=meanR_J*(mass_J/mass_S)^(2/5);
ratioRSOI_rJ=R_JSOI/r_J;

%Asteroid 1979XB
id_A=59;
r_A=r_Sat/100;    % The constant is computeted in order to have a good plot in the Transfer section
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

deltaV1opt=VI-V_Sat;      % delta speed that we have to provide to get into the First Transfer Leg 
deltaV1nopt=norm(deltaV1opt);
deltaV2opt=V_A-VF;        % delta speed that we have to provide to get into the Asteroid's Orbit
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


v1=sqrt(vinfan^2+2*mu_J/rp);     %velocity at the pericentre for the incoming hyperbola
e1=1+(rp*vinfan^2)/mu_J;
a1=rp/(1-e1);
p1=a1*(1-e1^2);
thetaSOI1=acos(1/e1*(p1/R_JSOI-1));
F1=2*atanh(sqrt((e1-1)/(e1+1))*tan(thetaSOI1/2));
n1=sqrt(mu_J/abs(a1)^3);
deltaT1=1/n1*(e1*sinh(F1)-F1);
deltaT1days=deltaT1/(3600*24);
v2=sqrt(vinfdn^2+2*mu_J/rp);     %velocity at the pericentre for the outgoing hyperbola
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

ratioTIMESOI_TIMEMISSION=deltaTdays/(mjd2000a-mjd2000d);
ratioRSOI_RJ=R_JSOI/norm(R_J);

deltavpopt=v2-v1;    % Delta V to be provided at the pericentre
deltavfbnopt=norm(deltavpopt);

deltaVfb=vinfd-vinfa;
ratioFB=norm(deltaVfb)/deltavfbnopt;

%Rotation of -beta1 with additional steps to get a clear graph in the flyby section

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
point1aux=centro1-44000000*vinfa;     %The constants are computeted in order to have a clear plot in the flyby section(40000000 if we use -9000000 s instead of -deltaT1 in tspanpast)
centro2=(rp+abs(a2))*evers;
point2aux=centro2+7000000*vinfd;
point3aux=centro2-30*rpvect;
vettore1=[centro1 point1aux]';
vettore2=[centro2 point2aux]';
vettore3=[centro1 point3aux]';



deltaVtotopt=deltaV1nopt+deltaV2nopt+deltavfbnopt;


%Natural Flyby

% Assumption:The versor u is the same as for rp and the impact parameter.
% In this way we can have a meaningful comparison

[vinfdnat,deltaVfbnat]= Rodrigues(vinfa,delta1,u);

VFB2nat=V_J+vinfdnat;

ratio_deltaVinfd=norm(vinfdnat)/norm(vinfd);      %It's less than 1 as expected (In the natural case we do not have the additional boost)

ratio_deltaVfb=norm(deltaVfbnat)/norm(deltaVfb);  %It's less than 1 as expected

ratio_Vinfnat=norm(vinfdnat)/norm(vinfa);         %It's 1 as expected

ratio_Vinf=norm(vinfd)/norm(vinfa);               %It's greater than 1  as expected (we have a powered flyby)

deltaVsavedFB=norm(VFB2-VFB1)-deltavfbnopt;       %It's the amount of deltaV saved by doing the flyby
%                                                  instead of an instantaneous manoeuvre with the same result of
%                                                  the one done thanks to Jupiter gravitational field

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

figure(3)
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


%% Fly-By Plot  
pointSun=-R_J/norm(R_J);
centre=[0;0;0];
vettSun=[centre 100*pointSun];
% tspanpast= linspace( 0, -9000000, 36000);    % Propagation from about 100 days before reaching the pericentre
% tspanfuture= linspace( 0, 1000000, 36000);   % Propagation till about 11 days after reaching the pericentre
tspanpast= linspace( 0, -deltaT1, 36000); 
tspanfuture= linspace( 0, deltaT2, 36000);
tspanfuturenat= tspanfuture;
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% Perform the integration
y0past=[rpvect;vp1vect];
y0future=[rpvect;vp2vect];
y0futurenat=[rpvect;vp1vect];
[ ~, ypast ] = ode113( @(t,y) ode_2bp(t,y,mu_J), tspanpast, y0past, options );
[ ~, yfuture ] = ode113( @(t,y) ode_2bp(t,y,mu_J), tspanfuture, y0future, options);
[ ~, yfuturenat ] = ode113( @(t,y) ode_2bp(t,y,mu_J), tspanfuturenat, y0futurenat, options);

figure(4)
plot3( ypast(:,1)/r_J, ypast(:,2)/r_J, ypast(:,3)/r_J, 'b-' );
hold on
plot3( yfuture(:,1)/r_J, yfuture(:,2)/r_J, yfuture(:,3)/r_J, 'r-' );
hold on
plot3(yfuture(1,1)/r_J, yfuture(1,2)/r_J, yfuture(1,3)/r_J, 'bo' );   
hold on
plot3(vettore1(:,1)/r_J,vettore1(:,2)/r_J,vettore1(:,3)/r_J);
hold on
plot3(vettore2(:,1)/r_J,vettore2(:,2)/r_J,vettore2(:,3)/r_J);
hold on
plot3(vettore3(:,1)/r_J,vettore3(:,2)/r_J,vettore3(:,3)/r_J,'--');
hold on
plot3(vettSun(1,:),vettSun(2,:),vettSun(3,:),'--');
hold on
plot3( yfuturenat(:,1)/r_J, yfuturenat(:,2)/r_J, yfuturenat(:,3)/r_J, '-' );
J.foto='Jupiter.jpg';
PlanetSphere(1,0,0,0,J);
xlabel('$x$ $[r_J]$','Interpreter','latex');
ylabel('$y$  $[r_J]$','Interpreter','latex');
zlabel('$z$  $[r_J]$','Interpreter','latex');
legend('$Incoming$ $hyperbola$','$Outgoing$ $hyperbola$','','$Asymptote$ $for$ $incoming$ $hyperbola$', ...
    '$Asymptote$ $for$ $outgoing$ $hyperbola$','$Apse$ $Line$','$Sun$ $Direction$', ...
    '$Outgoing$ $hyperbola$ $with$ $Natural$ $Flyby$','Interpreter','latex');
title('$Planetocentric$ $Flyby$','Interpreter','latex');

R1=ypast(end,1:3);
R2=yfuture(end,1:3);

ratioSOI1=norm(R1)/R_JSOI;         %It's 1 as expected 
ratioSOI2=norm(R2)/R_JSOI;         %It's 1 as expected (when we use deltaT2 instead of 1000000 s in tspanfuture)

%% Comparison of the first leg with a Hohmann

aHoh1=(norm(R_Sat)+norm(R_J))/2;
TOFH1=pi*sqrt(aHoh1^3/mu_S);  

ratio_comparison1=TOFt1/TOFH1;
comparison_days1=(TOFt1-TOFH1)/(3600*24);   % About 158 days,i.e. the 4% of the TOFt1 (the hypothesis of an Hohmann Transfer                                                         is reasonable and verified)

%% Comparison of the second leg with a Hohmann

aHoh2=(norm(R_J)+norm(R_A))/2;
TOFH2=pi*sqrt(aHoh2^3/mu_S);  

ratio_comparison2=TOFt2/TOFH2;
comparison_days2=(TOFH2-TOFt2)/(3600*24);   % About 171 days,i.e. the 12% of the TOFt2 (the hypothesis of an Hohmann Transfer                                                        isn't reasonable and it's not verified)



%% ------------------ RESULTS--------------------

%The total deltaV is 7.8658 km/s and we have:
%          deltaVd=1.9358 km/s
%          deltaVfb=0.5872 km/s
%          deltaVa=5.3427 km/s
%          DEPARTURE: 8-8-2033 AT 3:35
%          FLYBY:     29-7-2043 AT 6:11
%          ARRIVAL:   31-7-2047  AT 00:04

% Time
Departure=dated                  % [year, month, day, hour, minute, second]
Flyby=datefb                     % [year, month, day, hour, minute, second]
Arrival=datea                    % [year, month, day, hour, minute, second]


% First Leg
At1                              % [km]
Et1                              % [-]
inclit1=rad2deg(inclit1)         % [deg]
OMt1=rad2deg(OMt1)               % [deg]
omt1=rad2deg(omt1)               % [deg]
thetat1=rad2deg(thetat1)         % [deg]
deltaDeparture=deltaV1nopt       % [km/s]

% Second Leg
At2                              % [km]
Et2                              % [-]
inclit2=rad2deg(inclit2)         % [deg]
OMt2=rad2deg(OMt2)               % [deg]
omt2=rad2deg(omt2)               % [deg]
thetat2=rad2deg(thetat2)         % [deg]
deltaArrival=deltaV2nopt         % [km/s]

%Flyby
altitude                         % [km]
deltaTdays                       % [days]
ratioFB                          % [-]
deltaPowered=deltavpopt          % [km/s]
