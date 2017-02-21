function Target = Payload(x,k,j,u, phi0, zeta0)
mScale = 1;
mfuel_burn = x/mScale;



[AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma] = ThirdStageSim(mfuel_burn,k,j,u, phi0, zeta0);

Target = -mpayload;
% Target = (Alt(end) - 160000)^2 -mpayload*10000;
end