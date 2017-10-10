function [adigatorFunInfo, adigatorOutputs] = adigatortempfunc5(adigatorFunInfo,adigatorInputs)
[flag, adigatorFunInfo, adigatorInputs] = adigatorFunctionInitialize(5,adigatorFunInfo,adigatorInputs);
if flag; adigatorOutputs = adigatorInputs; return; end;
M = adigatorInputs{1};
Alpha = adigatorInputs{2};
auxdata = adigatorInputs{3};
T0 = adigatorInputs{4};
P0 = adigatorInputs{5};
nargin = 5; nargout = 3;adigatorVarAnalyzer('% Engine Interpolator for engine data');
T1 = auxdata.tempgridded(M,Alpha).*T0;
T1 = adigatorVarAnalyzer('T1 = auxdata.tempgridded(M,Alpha).*T0;',T1,'T1',0);
P1 = auxdata.presgridded(M,Alpha).*P0;
P1 = adigatorVarAnalyzer('P1 = auxdata.presgridded(M,Alpha).*P0;',P1,'P1',0);
M1 = auxdata.M1gridded(M, Alpha);
M1 = adigatorVarAnalyzer('M1 = auxdata.M1gridded(M, Alpha);',M1,'M1',0);
Isp = auxdata.IspGridded(M1,T1);
Isp = adigatorVarAnalyzer('Isp = auxdata.IspGridded(M1,T1);',Isp,'Isp',0);
adigatorVarAnalyzer('% eq = scattered.equivalence(M1,T1);');
eq = auxdata.eqGridded(M1,T1);
eq = adigatorVarAnalyzer('eq = auxdata.eqGridded(M1,T1);',eq,'eq',0);
% ADiGator FOR Statement #2: START
cadaforvar2 = 1: length(eq);
cadaforvar2 = adigatorVarAnalyzer('cadaforvar2 = 1: length(eq);',cadaforvar2,'cadaforvar2',0);
[adigatorForVariable2, adigatorForEvalStr, adigatorForEvalVar] = adigatorForInitialize(2,cadaforvar2,0);%#ok<NASGU>
if ~isempty(adigatorForEvalStr)
    adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorForEvalStr); adigatorSetCellEvalFlag(0);
end
for adigatorForVariable2i = adigatorForVariable2
cadaforcount2 = adigatorForIterStart(2,adigatorForVariable2i);
i = cadaforvar2(:,cadaforcount2);
i = adigatorVarAnalyzer('i = cadaforvar2(:,cadaforcount2);',i,'i',0);
    % ADiGator IF Statement #1: START
    cadaconditional1 = eq(i) > 1;
    cadaconditional1 = adigatorVarAnalyzer('cadaconditional1 = eq(i) > 1;',cadaconditional1,'cadaconditional1',0);
    adigatorIfInitialize(1,cadaconditional1);
    adigatorIfIterStart(1,1);
        if ~exist('eq','var'); eq = cadastruct([],'eq',[],0); end
        eq(i) = 1;
        eq = adigatorVarAnalyzer('eq(i) = 1;',eq,'eq',1);
    [adigatorIfEvalStr, adigatorIfEvalVar] = adigatorIfIterEnd(1,1);%#ok<NASGU>
    if ~isempty(adigatorIfEvalStr)
adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorIfEvalStr); adigatorSetCellEvalFlag(0);
end
    % ADiGator IF Statement #1: END
[adigatorForEvalStr, adigatorForEvalVar]= adigatorForIterEnd(2,adigatorForVariable2i);
if ~isempty(adigatorForEvalStr)
    adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorForEvalStr); adigatorSetCellEvalFlag(0);
end
end
% ADiGator FOR Statement #2: END
wcap = 0.65;
wcap = adigatorVarAnalyzer('wcap = 0.65;',wcap,'wcap',0);
adigatorVarAnalyzer('% From C-REST thrust calculator');
wcapstandard = 0.5395;
wcapstandard = adigatorVarAnalyzer('wcapstandard = 0.5395;',wcapstandard,'wcapstandard',0);
Acapstandard = 0.3057 + 0.06891*(5 + Alpha);
Acapstandard = adigatorVarAnalyzer('Acapstandard = 0.3057 + 0.06891*(5 + Alpha);',Acapstandard,'Acapstandard',0);
Acap = (wcap./wcapstandard).^2.*Acapstandard;
Acap = adigatorVarAnalyzer('Acap = (wcap./wcapstandard).^2.*Acapstandard;',Acap,'Acap',0);
Mrat = M./M1;
Mrat = adigatorVarAnalyzer('Mrat = M./M1;',Mrat,'Mrat',0);
MM = 0.175097746892639.*M1.^(0.776790959520025).*(Mrat).^(0.57952831191643) - 0.121263698193072;
MM = adigatorVarAnalyzer('MM = 0.175097746892639.*M1.^(0.776790959520025).*(Mrat).^(0.57952831191643) - 0.121263698193072;',MM,'MM',0);
a1 = -0.08084;
a1 = adigatorVarAnalyzer('a1 = -0.08084;',a1,'a1',0);
a2 = 0.9422;
a2 = adigatorVarAnalyzer('a2 = 0.9422;',a2,'a2',0);
a3 = 0.7429;
a3 = adigatorVarAnalyzer('a3 = 0.7429;',a3,'a3',0);
a4 = -0.6744;
a4 = adigatorVarAnalyzer('a4 = -0.6744;',a4,'a4',0);
mc = a1 + a2*MM + a3*MM.^2 + a4*MM.^3;
mc = adigatorVarAnalyzer('mc = a1 + a2*MM + a3*MM.^2 + a4*MM.^3;',mc,'mc',0);
gam0 = 1.4000000;
gam0 = adigatorVarAnalyzer('gam0 = 1.4000000;',gam0,'gam0',0);
RR = 287.035;
RR = adigatorVarAnalyzer('RR = 287.035;',RR,'RR',0);
w = 0.9.*mc.*Acap.*P0.*M.*sqrt(gam0./RR./T0)*4.0 ;
w = adigatorVarAnalyzer('w = 0.9.*mc.*Acap.*P0.*M.*sqrt(gam0./RR./T0)*4.0 ;',w,'w',0);
fst = 0.0291;
fst = adigatorVarAnalyzer('fst = 0.0291;',fst,'fst',0);
wf = fst.*w.*eq;
wf = adigatorVarAnalyzer('wf = fst.*w.*eq;',wf,'wf',0);
adigatorOutputs = {Isp;wf;eq};
[adigatorFunInfo, adigatorOutputs] = adigatorFunctionEnd(5,adigatorFunInfo,adigatorOutputs);
