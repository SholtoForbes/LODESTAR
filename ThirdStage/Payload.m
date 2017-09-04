function Target = Payload(x,k,j,u, phi0, zeta0, lb,num_div,plotflag)
mScale = 1;
mfuel_burn = x/mScale;



[AltF_actual, vF, Alt, v, t, mpayload, Alpha, m,AoA_init,q,gamma,D,AoA_max,zeta,phi, inc,Vec_angle,T,CL,L] = ThirdStageSim(mfuel_burn,k,j,u, phi0, zeta0, lb,num_div,plotflag);

Target = -mpayload ;

% Target = -vF ;

end