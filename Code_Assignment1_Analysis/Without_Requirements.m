clear;
close all;
clc;
%
% Script for data in the report
%
%% Parameters definition
datedi=[2028,1,1,0,0,0];
datefbi=[2028,1,1,0,0,0];     
dateai=[2028,1,1,0,0,0];       
mjd2000di= date2mjd2000(datedi);


mjd2000fbi= date2mjd2000(datefbi);

mjd2000ai= date2mjd2000(dateai);        

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
nd=T_Satdays/96;

%Jupiter
ibody_J=5;
mu_J=astroConstants(15);
r_J=astroConstants(25);
[kep_J,~] = uplanet(mjd2000fbi, ibody_J);
T_J=2*pi*sqrt(kep_J(1)^3/mu_S);
T_Jdays=T_J/(3600*24);
nfb=T_Jdays/96;
%nfb=T_Jdays/48;  %nfb=T_Jdays/96;
[R_J, ~] = kep2car (kep_J(1),kep_J(2),kep_J(3),kep_J(4),kep_J(5),kep_J(6), mu_S);
atmosphheight=5000;   %[km]


%Asteroid 1979XB
id_A=59;
r_A=r_Sat/100;      % We decided this value only for aesthetic reasons
[kep_A,~,~] = ephNEO(mjd2000ai,id_A);
T_A=2*pi*sqrt(kep_A(1)^3/mu_S);
T_Adays=T_A/(3600*24);
na=T_Adays/96;
[R_A, ~] = kep2car (kep_A(1),kep_A(2),kep_A(3),kep_A(4),kep_A(5),kep_A(6), mu_S);

%% Lambert

% md=52;
% mfb=124;
% ma=404;
md=98;
mfb=243;    %parameters taken to optimally cover all windows
ma=800;

datedi=[2028,1,1,0,0,0];
mjd2000di= date2mjd2000(datedi);
mjd2000df=mjd2000di+md*nd;                
datedf=mjd20002date(mjd2000df);

datefbi=[2028,1,1,0,0,0];   
mjd2000fbi= date2mjd2000(datefbi);
mjd2000fbf=mjd2000fbi+mfb*nfb;                
datefbf=mjd20002date(mjd2000fbf);


dateai=[2028,1,1,0,0,0]; 
mjd2000ai= date2mjd2000(dateai);
mjd2000af=mjd2000ai+ma*na;                
dateaf=mjd20002date(mjd2000af);
%%
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
FLAG=NaN(md+1,mfb+1,ma+1);


for q=mjd2000di:nd:mjd2000df
    
    [kep_Sat,~] = uplanet(q, ibody_S);   
    [R_Sat, V_Sat] = kep2car (kep_Sat(1),kep_Sat(2),kep_Sat(3),kep_Sat(4), ...
    kep_Sat(5),kep_Sat(6), mu_S);  
    RI=R_Sat;




    for w=mjd2000fbi:nfb:mjd2000fbf
  
        
        TOFt1=(w-q)*24*60*60;
        
        if (TOFt1> 0  && TOFt1<T_Sat)                    % We have only an energy constraint    
            [kep_J,~] = uplanet(w, ibody_J);
            [R_J, V_J] = kep2car (kep_J(1),kep_J(2),kep_J(3),kep_J(4),kep_J(5),kep_J(6), mu_S);
            RFB=R_J;
            [~,~,Et1,~,VI,VFB1,TPAR1,~] = lambertMR(RI,RFB,TOFt1,MU,orbitType,Nrev,Ncase,optionsLMR);
            VI=VI';
            VFB1=VFB1';
            deltaV1=VI-V_Sat;      
            deltaV1n=norm(deltaV1);
            if TOFt1>TPAR1
            deltaVMatrix1(y,u)=deltaV1n;
            TOF1Matrix(y,u)=TOFt1/(3600*24);    % [days]
          
                for t=mjd2000ai:na:mjd2000af
                    
                    TOFt2=(t-w)*24*60*60;
                    
                    if (TOFt2>0 && TOFt2<=1.1*T_J)   % We have only an energy constraint
                        [kep_A,~,~] = ephNEO(t,id_A);
                        [R_A, V_A] = kep2car (kep_A(1),kep_A(2),kep_A(3),kep_A(4),kep_A(5),kep_A(6), mu_S);
                        
                        RF=R_A;
                        [~,~,~,~,VFB2,VF,TPAR2,~] = lambertMR(RFB,RF,TOFt2,MU,orbitType,Nrev,Ncase,optionsLMR);
                     if TOFt2>TPAR2    
                            VFB2=VFB2';
                            VF=VF';
                            deltaV2=V_A-VF;        
                            deltaV2n=norm(deltaV2);
                            TOF2Matrix(u,p)=TOFt2/(3600*24);    % [days]
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
                            if (altitude>atmosphheight && flag==1)   % 5000km of atmosphere and condition on FLAG 
                                v1=sqrt(vinfan^2+2*mu_J/rp);     
                                v2=sqrt(vinfdn^2+2*mu_J/rp);     
                                deltavp=v2-v1;    
                                deltafbn=norm(deltavp);
                                deltaVMatrix2(u,p)=deltaV2n;
                                deltaVMatrixfb(y,u,p)=deltafbn;
                                deltaVTOT(y,u,p)=deltaV1n+deltaV2n+deltafbn;  
                            end
                        
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


deltav2min=min(min(deltaVMatrix2));
[rowmina,columnmina]=find(deltaVMatrix2==deltav2min);
mjd2000fbmina=mjd2000fbi+nfb*(rowmina-1);
mjd2000amina=mjd2000ai+na*(columnmina-1);

deltavfbmin=min(min(min(deltaVMatrixfb)));

%% Contour Plots with all the informations

figure(1)
contour(vettdep,vettfb,deltaVMatrix1',linspace(deltav1min,deltav1min+2,20),'ShowText','off');
c=colorbar;
c.Label.String= "\Deltav_{departure} [km/s]";     
ha=gca;
clim(ha,deltav1min+[0 1.5]);
clim(ha,'manual');
hold on
contour(vettdep,vettfb,TOF1Matrix',1000:1000:10000,'k','ShowText','on');
hold on
plot(mjd2000d,mjd2000fb,'*');
hold on
plot(mjd2000dmind,mjd2000fbmind,'o');
title('$Departure$ $Cost$ $Plot$','Interpreter','latex')
xlabel('$Departure$ $Date$ $from$ $Saturn$','Interpreter','latex')
ylabel('$Flyby$ $Date$ $of$ $Jupiter$','Interpreter','latex')
legend('$deltaV$','$TOF$','$deltaVott$','$deltaVmin$','Interpreter','latex');
grid on

figure(2)
contour(vettfb,vettarr,deltaVMatrix2',linspace(deltav2min,deltav2min+6,20),'ShowText','off');
c=colorbar;
c.Label.String= "\Deltav_{arrival} [km/s]";     
ha=gca;
clim(ha,deltav2min+[0 5.5]);
clim(ha,'manual');
hold on
contour(vettfb,vettarr,TOF2Matrix',800:500:6000,'k','ShowText','on');
hold on
plot(mjd2000fb,mjd2000a,'*r');
hold on
plot(mjd2000fbmina,mjd2000amina,'om');
title('$Arrival$ $Cost$ $Plot$','Interpreter','latex')
xlabel('$Flyby$ $Date$ $of$ $Jupiter$','Interpreter','latex')
ylabel('$Arrival$ $Date$ $to$ $the$ $Asteroid$','Interpreter','latex')
legend('$deltaV$','$TOF$','$deltaVott$','$deltaVmin$','Interpreter','latex');
grid on

%%
figure(3)
% For the cost plot:
% Use meshgrid to create 2-D grid coordinates based on the
% time discretizations for both departure and arrival.
%   [X,Y] = meshgrid(x,y)
startDate1 = datenum('01-01-2028', 'dd mmm yyyy');
startDate2 = datenum('24-02-2058', 'dd mmm yyyy'); %fine dep
endDate1 = datenum('01-01-2028',  'dd mmm yyyy');
endDate2 = datenum('14-01-2058',  'dd mmm yyyy');
xData = linspace(startDate1, startDate2, length(vettdep));
yData = linspace(endDate1, endDate2, length(vettfb));
[T1, T2] = meshgrid(xData, yData);
contour(T1,T2,deltaVMatrix1',linspace(deltav1min,deltav1min+6,20),'ShowText','off');
c=colorbar;
c.Label.String= "\Deltav_{departure} [km/s]";     
ha=gca;
clim(ha,deltav1min+[0 5.5]);
clim(ha,'manual');
hold on
contour(T1,T2,TOF1Matrix',1000:1000:10000,'k','ShowText','on');
title('$Departure$ $Cost$ $Plot$','Interpreter','latex');
xlabel('$Departure$ $Date$ $from$ $Saturn$','Interpreter','latex');
ylabel('$Flyby$ $Date$ $of$ $Jupiter$','Interpreter','latex');
legend('$deltaV$','$TOF$','Interpreter','latex');
datetick('x', 'dd-mmm-yyyy', 'keeplimits', 'keepticks');
datetick('y', 'dd-mmm-yyyy', 'keeplimits', 'keepticks');
grid on


figure(4)
startDate1 = datenum('01-01-2028', 'dd mmm yyyy');
startDate2 = datenum('14-01-2058', 'dd mmm yyyy'); %fine flyby
endDate1 = datenum('01-01-2028',  'dd mmm yyyy');
endDate2 = datenum('08-01-2058',  'dd mmm yyyy');
xData = linspace(startDate1, startDate2, length(vettfb));
yData = linspace(endDate1, endDate2, length(vettarr));
[T1, T2] = meshgrid(xData, yData);

contour(T1,T2,deltaVMatrix2',linspace(deltav2min,deltav2min+10,20),'ShowText','off');
c=colorbar;
c.Label.String= "\Deltav_{arrival} [km/s]";     
ha=gca;
clim(ha,deltav2min+[0 9.5]);
clim(ha,'manual');
hold on
contour(T1,T2,TOF2Matrix',800:800:6000,'k','ShowText','on');
title('$Arrival$ $Cost$ $Plot$','Interpreter','latex')
xlabel('$Flyby$ $Date$ $of$ $Jupiter$','Interpreter','latex')
ylabel('$Arrival$ $Date$ $to$ $the$ $Asteroid$','Interpreter','latex')
legend('$deltaV$','$TOF$','Interpreter','latex');
datetick('x', 'dd-mmm-yyyy', 'keeplimits', 'keepticks');
datetick('y', 'dd-mmm-yyyy', 'keeplimits', 'keepticks');
grid on

% New potentially interesting region characterized by a flyby relatively close 
% to the beginning of the windows
% PrematureFlyby:         
%          DEPARTURE: 
%          FLYBY:     13165   [mjd2000]
%          ARRIVAL:   14045   [mjd2000]
