function [Isp,wf] = RESTM12int(M, Alpha, t_ratio, Efficiency, scattered, SPARTAN_SCALE)
% Engine Interpolator for RESTM12 Data
% Reading from RESTM12DATA 


% T150kPa = scattered.temp(M,Alpha); % from communicator, get temp behind shock at 50kpa from communicator (T1,50kPa)
% T1 = T150kPa.*t_ratio;
% P1 = scattered.pres(M,Alpha).*Efficiency; % note this is at 50kPa, modified by efficiency
% M1 = scattered.M1(M, Alpha);

T150kPa = scattered.tempgridded(M,Alpha); % from communicator, get temp behind shock at 50kpa from communicator (T1,50kPa)
T1 = T150kPa.*t_ratio;
P1 = scattered.presgridded(M,Alpha).*Efficiency; % note this is at 50kPa, modified by efficiency
M1 = scattered.M1gridded(M, Alpha);


% Isp = scattered.IspScattered(M1,T1);

Isp = scattered.IspGridded(M1,T1);

% 
% phi = scattered.phi(M1,T1)

phi = scattered.equivalence(M1,T1);
% data = scattered.data;

% Isp = griddata(data(:,1),data(:,2),data(:,6),M1,T1,'cubic');
% 
% phi = griddata(data(:,1),data(:,2),data(:,7),M1,T1,'cubic');

wcap = 0.65*SPARTAN_SCALE^(2/3);

wcapstandard = 0.2156; %meters

R0 = wcap^2/wcapstandard^2; % meters

Acap = 0.0470*R0; %                   m^2          

gam0 = 1.4000000;
r = 287.035;

a1 = 0.5730261;
a2 = 1.65047e-2;
a3 = 8.889289e-3;
a4 = -6.0995262e-4;

mc     = a1 + a2*M1 + a3*M1.^2 + a4*M1.^3;

w = mc.*Acap.*P1.*M1.*sqrt(gam0./r./T1).*4.0; % 4 engine modules in Full capture engine 

fst = 0.0291;
wf = fst.*w.*phi; %(kg/s)
end


