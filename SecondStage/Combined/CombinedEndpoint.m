function output = CombinedEndpoint(input)

auxdata = input.auxdata;

returnflag = auxdata.returnflag;

t01 = input.phase(1).initialtime;
tf1 = input.phase(1).finaltime;
x01 = input.phase(1).initialstate;
xf1 = input.phase(1).finalstate;

t021 = input.phase(2).initialtime;
tf21 = input.phase(2).finaltime;
x021 = input.phase(2).initialstate;
xf21 = input.phase(2).finalstate;

if returnflag
t022 = input.phase(3).initialtime;
tf22 = input.phase(3).finaltime;
x022 = input.phase(3).initialstate;
xf22 = input.phase(3).finalstate;
end

%%

output.eventgroup(1).event = [x021(1)-xf1(1) x021(2)-xf1(9) x021(3)-xf1(8) x021(4)-xf1(2) x021(5)-xf1(4) x021(6)-xf1(6) x021(7)-xf1(5)];


%%
if returnflag
output.eventgroup(2).event = [x022(1:9)-xf21(1:9) t021-tf1];

output.eventgroup(3).event = [t022-tf21];

else

output.eventgroup(2).event = [t021-tf1];
    
end

%%
phase_no = auxdata.phase_no;

altf3 = input.phase(phase_no).finalstate(1);
vf3 = input.phase(phase_no).finalstate(2);
gammaf3 = input.phase(phase_no).finalstate(3);
mf3 = input.phase(phase_no).finalstate(4);
phif3 = input.phase(phase_no).finalstate(6);
zetaf3 = input.phase(phase_no).finalstate(7);
%%

t03 = input.phase(phase_no).initialtime;
tf3 = input.phase(phase_no).finaltime;
x03 = input.phase(phase_no).initialstate;
xf3 = input.phase(phase_no).finalstate;


[AltF_actual, vF, Alt, v, t, mpayload, Alpha, m,q,gamma,D,zeta,phi, inc,T,CL,L,inc_diff] = ThirdStageSim(altf3,gammaf3,vf3, phif3, 0,zetaf3, mf3, input.auxdata);

output.eventgroup(phase_no).event =[xf21(1)-x03(1) xf21(4)-x03(2) xf21(5)-x03(3) xf21(6)-x03(7) xf21(3)-x03(6) t03-tf21 tf3-t03 AltF_actual inc_diff];



%% cost
output.objective = -mpayload;


end
