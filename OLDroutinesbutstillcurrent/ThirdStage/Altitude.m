function [c, Constraint] = Altitude(x)
mScale = .01;
AoA_list = x(1:end-1)/AScale ;
mfuel_burn = x(end)/mScale;

[AltF, vF, Alt, v, t, mpayload, Alpha] = ThirdStageSim(AoA_list, mfuel_burn);
Constraint = AltF - 566890;
c = [];
% Constraint = []
end