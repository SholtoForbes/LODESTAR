function [adigatorFunInfo, adigatorOutputs] = adigatortempfunc3(adigatorFunInfo,adigatorInputs)
[flag, adigatorFunInfo, adigatorInputs] = adigatorFunctionInitialize(3,adigatorFunInfo,adigatorInputs);
if flag; adigatorOutputs = adigatorInputs; return; end;
angleInDegrees = adigatorInputs{1};
nargin = 1; nargout = 1;adigatorVarAnalyzer('% DEG2RAD Convert angles from degrees to radians');
adigatorVarAnalyzer('%');
adigatorVarAnalyzer('%   DEG2RAD has been replaced by DEGTORAD.');
adigatorVarAnalyzer('%');
adigatorVarAnalyzer('%   angleInRadians = DEG2RAD(angleInDegrees) converts angle units from');
adigatorVarAnalyzer('%   degrees to radians.');
adigatorVarAnalyzer('% Copyright 2007-2009 The MathWorks, Inc.');
angleInRadians = (pi/180) * angleInDegrees;
angleInRadians = adigatorVarAnalyzer('angleInRadians = (pi/180) * angleInDegrees;',angleInRadians,'angleInRadians',0);
adigatorOutputs = {angleInRadians};
[adigatorFunInfo, adigatorOutputs] = adigatorFunctionEnd(3,adigatorFunInfo,adigatorOutputs);
