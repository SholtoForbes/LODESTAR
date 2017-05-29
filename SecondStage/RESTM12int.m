function [Isp,wf,eq] = RESTM12int(M, Alpha, scattered, SPARTAN_SCALE,T0,P0)
% Engine Interpolator for RESTM12 Data
% Reading from RESTM12DATA 

% THIS IS ACTUALLY FOR THE CREST M10

% NEW SPARTAN MODEL

T1 = scattered.tempgridded(M,Alpha).*T0;
P1 = scattered.presgridded(M,Alpha).*P0; % note this is at 50kPa, modified by efficiency
M1 = scattered.M1gridded(M, Alpha);

Isp = scattered.IspGridded(M1,T1);

% eq = scattered.equivalence(M1,T1);
eq = scattered.eqGridded(M1,T1);
for i = 1: length(eq)
    if eq(i) > 1
        eq(i) = 1;
    end
end


wcap = 0.65;

wcapstandard = 0.2156; %meters

R0 = wcap^2/wcapstandard^2; % meters

Acap = 0.0470*R0*SPARTAN_SCALE^(2/3); %                   m^2          

gam0 = 1.4000000;
r = 287.035;

a1 = 0.5730261;
a2 = 1.65047e-2;
a3 = 8.889289e-3;
a4 = -6.0995262e-4;

mc     = a1 + a2*M1 + a3*M1.^2 + a4*M1.^3;

w = mc.*Acap.*P1.*M1.*sqrt(gam0./r./T1).*4.0; % 4 engine modules in Full capture engine 

fst = 0.0291;
wf = fst.*w.*eq; %(kg/s)
end


