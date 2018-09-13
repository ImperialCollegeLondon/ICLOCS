function Weight_dot = calcWeight(V,H,throttle,FFModel)

%% constants

% coefficients for maximum throttle setting
am1 = -1.5374*10^-13;
am2 = 2.3875*10^-11;
am3 = -5.7989*10^-10;
bm1 = 2.1483*10^-9;
bm2 = -3.2731*10^-7;
bm3 = -9.1353*10^-7;
cm1 = -5.8602*10^-6;
cm2 = 9.7605*10^-4;
cm3 = 0.11389; 

% coefficients for minimum throttle setting (i for idle)
ai1 = 4.2242*10^-15;         %                      [-]
ai2 = -9.0868*10^-13;        %                      [-]
ai3 = 1.6801*10^-10;         %                      [-]   
bi1 = -4.8122*10^-11;        %                      [-]
bi2 = 1.0818*10^-8;          %                      [-]
bi3 = -2.7199*10^-6;         %                      [-] 
ci1 = 2.5788*10^-8;          %                      [-]
ci2 = -6.8894*10^-6;         %                      [-]
ci3 = 0.025804;              %                      [-]   
Vcas=V;
  
%% Fuel flow calculations
% Maximum fuel flow

am=am1.*Vcas.^2+am2.*Vcas+am3;  %                      [-]
bm=bm1.*Vcas.^2+bm2.*Vcas+bm3;  %                      [-]
cm=cm1.*Vcas.^2+cm2.*Vcas+cm3;  %                      [-]
Fmax=am.*H.^2+bm.*H+cm;         % maximum fuel flow    [kg/s]   

% Minimum fuel flow (idle throttle setting)

ai=ai1.*Vcas.^2+ai2.*Vcas+ai3;  %                      [-]
bi=bi1.*Vcas.^2+bi2.*Vcas+bi3;  %                      [-]
ci=ci1.*Vcas.^2+ci2.*Vcas+ci3;  %                      [-]
Fidle=ai.*H.^2+bi.*H+ci;        % minimum fuel flow    [kg/s]

% Interpolation
% FF= Fidle+(Fmax-Fidle).*((throttle-0)/(1-0)); % Current fuel flow [kg/s]
FF= Fidle+(Fmax-Fidle).*ppval(FFModel,throttle); % Current fuel flow [kg/s]

% Calculation of weight
Weight_dot=-2*FF*9.81;

