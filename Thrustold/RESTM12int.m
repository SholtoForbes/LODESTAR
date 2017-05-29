function [Isp,wf,Pnom] = RESTM12int(M1, T1, P1, scattered)
% Engine Interpolator for RESTM12 Data
% Reading from RESTM12DATA 


Isp = scattered.IspScattered(M1,T1);

Pnom = scattered.PnomScattered(M1,T1);

wcap = 0.65;

wcapstandard = 0.2156; %meters

R0 = wcap^2/wcapstandard^2; % meters

Acap = 0.0470*R0; %                   m^2          

gam0 = 1.4000000;
r = 287.035;

a1 = 0.5730261;
a2 = 1.65047e-2;
a3 = 8.889289e-3;
a4 = -6.0995262e-4;

mc     = a1 + a2*M1 + a3*M1^2 + a4*M1^3;

w = mc*Acap*P1*M1*sqrt(gam0/r/T1)*4.0; % 4 engine modules in Full capture engine 

phi = 1; 
fst = 0.0291;
wf = fst*w*phi; %(kg/s)
end


