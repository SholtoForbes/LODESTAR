function [Isp,wf,eq,q1,w] = RESTint(M, Alpha, auxdata,T0,P0)
% Engine Interpolator for engine data


T1 = auxdata.interp.tempgridded(M,Alpha).*T0;
P1 = auxdata.interp.presgridded(M,Alpha).*P0; % note this is at 50kPa, modified by efficiency
M1 = auxdata.interp.M1gridded(M, Alpha);

q1 = 0.5.*1.4.*P1.*M1.^2;

Isp = auxdata.interp.IspGridded(M1,T1);

% eq = scattered.equivalence(M1,T1);
eq = auxdata.interp.eqGridded(M1,T1);
for i = 1: length(eq)
    if eq(i) > 1
        eq(i) = 1;
    end
end

wcap = 0.65;

% From C-REST thrust calculator

wcapstandard = 0.5395; % meters

Acapstandard = 0.3057 + 0.06891*(5 + Alpha);

Acap = (wcap./wcapstandard).^2.*Acapstandard;

Mrat = M./M1;

MM = 0.175097746892639.*M1.^(0.776790959520025).*(Mrat).^(0.57952831191643) - 0.121263698193072;

a1 = -0.08084;
a2 = 0.9422;
a3 = 0.7429;
a4 = -0.6744;

mc = a1 + a2*MM + a3*MM.^2 + a4*MM.^3;

gam0 = 1.4000000;
RR = 287.035;

w = 0.9.*mc.*Acap.*P0.*M.*sqrt(gam0./RR./T0)*4.0 ;

fst = 0.0291;
wf = fst.*w.*eq; %(kg/s)
end


