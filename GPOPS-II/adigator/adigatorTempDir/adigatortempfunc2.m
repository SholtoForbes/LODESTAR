function [adigatorFunInfo, adigatorOutputs] = adigatortempfunc2(adigatorFunInfo,adigatorInputs)
[flag, adigatorFunInfo, adigatorInputs] = adigatorFunctionInitialize(2,adigatorFunInfo,adigatorInputs);
if flag; adigatorOutputs = adigatorInputs; return; end;
x = adigatorInputs{1};
params = adigatorInputs{2};
nargin = 2; nargout = 1;adigatorVarAnalyzer('%GAUSSMF Gaussian curve membership function.');
adigatorVarAnalyzer('%   GAUSSMF(X, PARAMS) returns a matrix which is the Gaussian');
adigatorVarAnalyzer('%   membership function evaluated at X. PARAMS is a 2-element vector');
adigatorVarAnalyzer('%   that determines the shape and position of this membership function.');
adigatorVarAnalyzer('%   Specifically, the formula for this membership function is:');
adigatorVarAnalyzer('%');
adigatorVarAnalyzer('%   GAUSSMF(X, [SIGMA, C]) = EXP(-(X - C).^2/(2*SIGMA^2));');
adigatorVarAnalyzer('%');
adigatorVarAnalyzer('%   For example:');
adigatorVarAnalyzer('%');
adigatorVarAnalyzer('%       x = (0:0.1:10)'';');
adigatorVarAnalyzer('%       y1 = gaussmf(x, [0.5 5]);');
adigatorVarAnalyzer('%       y2 = gaussmf(x, [1 5]);');
adigatorVarAnalyzer('%       y3 = gaussmf(x, [2 5]);');
adigatorVarAnalyzer('%       y4 = gaussmf(x, [3 5]);');
adigatorVarAnalyzer('%       subplot(211); plot(x, [y1 y2 y3 y4]);');
adigatorVarAnalyzer('%       y1 = gaussmf(x, [1 2]);');
adigatorVarAnalyzer('%       y2 = gaussmf(x, [1 4]);');
adigatorVarAnalyzer('%       y3 = gaussmf(x, [1 6]);');
adigatorVarAnalyzer('%       y4 = gaussmf(x, [1 8]);');
adigatorVarAnalyzer('%       subplot(212); plot(x, [y1 y2 y3 y4]);');
adigatorVarAnalyzer('%       set(gcf, ''name'', ''gaussmf'', ''numbertitle'', ''off'');');
adigatorVarAnalyzer('%');
adigatorVarAnalyzer('%   See also DSIGMF,EVALMF, GAUSS2MF, GBELLMF, MF2MF, PIMF, PSIGMF, SIGMF,');
adigatorVarAnalyzer('%   SMF, TRAPMF, TRIMF, ZMF.');
adigatorVarAnalyzer('%   Roger Jang, 6-29-93, 10-5-93.');
adigatorVarAnalyzer('%   Copyright 1994-2002 The MathWorks, Inc.');
% ADiGator IF Statement #1: START
cadaconditional1 = nargin ~= 2;
cadaconditional1 = adigatorVarAnalyzer('cadaconditional1 = nargin ~= 2;',cadaconditional1,'cadaconditional1',0);
cadaconditional2 = length(params) < 2;
cadaconditional2 = adigatorVarAnalyzer('cadaconditional2 = length(params) < 2;',cadaconditional2,'cadaconditional2',0);
cadaconditional3 = params(1) == 0;
cadaconditional3 = adigatorVarAnalyzer('cadaconditional3 = params(1) == 0;',cadaconditional3,'cadaconditional3',0);
adigatorIfInitialize(1,cadaconditional1,cadaconditional2,cadaconditional3);
adigatorIfIterStart(1,1);
    adigatorError(1,1,'error(''Two arguments are required by the Gaussian MF.'');');
adigatorIfIterEnd(1,1);
[adigatorIfEvalStr, adigatorIfEvalVar] = adigatorIfIterStart(1,2);%#ok<NASGU>
if ~isempty(adigatorIfEvalStr)
adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorIfEvalStr); adigatorSetCellEvalFlag(0);
end
    adigatorError(1,2,'error(''The Gaussian MF needs at least two parameters.'');');
adigatorIfIterEnd(1,2);
[adigatorIfEvalStr, adigatorIfEvalVar] = adigatorIfIterStart(1,3);%#ok<NASGU>
if ~isempty(adigatorIfEvalStr)
adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorIfEvalStr); adigatorSetCellEvalFlag(0);
end
    adigatorError(1,3,'error(''The Gaussian MF needs a non-zero sigma.'');');
[adigatorIfEvalStr, adigatorIfEvalVar] = adigatorIfIterEnd(1,3);%#ok<NASGU>
if ~isempty(adigatorIfEvalStr)
adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorIfEvalStr); adigatorSetCellEvalFlag(0);
end
% ADiGator IF Statement #1: END
sigma = params(1);
sigma = adigatorVarAnalyzer('sigma = params(1);',sigma,'sigma',0);
c = params(2);
c = adigatorVarAnalyzer('c = params(2);',c,'c',0);
y = exp(-(x - c).^2/(2*sigma^2));
y = adigatorVarAnalyzer('y = exp(-(x - c).^2/(2*sigma^2));',y,'y',0);
adigatorOutputs = {y};
[adigatorFunInfo, adigatorOutputs] = adigatorFunctionEnd(2,adigatorFunInfo,adigatorOutputs);
