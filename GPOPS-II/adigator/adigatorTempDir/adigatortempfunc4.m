function [adigatorFunInfo, adigatorOutputs] = adigatortempfunc4(adigatorFunInfo,adigatorInputs)
[flag, adigatorFunInfo, adigatorInputs] = adigatorFunctionInitialize(4,adigatorFunInfo,adigatorInputs);
if flag; adigatorOutputs = adigatorInputs; return; end;
angleInRadians = adigatorInputs{1};
nargin = 1; nargout = 1;adigatorVarAnalyzer('% RAD2DEG Convert angles from radians to degrees');
adigatorVarAnalyzer('%');
adigatorVarAnalyzer('%   RAD2DEG has been replaced by RADTODEG.');
adigatorVarAnalyzer('%');
adigatorVarAnalyzer('%   angleInDegrees = RAD2DEG(angleInRadians) converts angle units from');
adigatorVarAnalyzer('%   radians to degrees.');
adigatorVarAnalyzer('% Copyright 2007-2008 The MathWorks, Inc.');
angleInDegrees = (180/pi) * angleInRadians;
angleInDegrees = adigatorVarAnalyzer('angleInDegrees = (180/pi) * angleInRadians;',angleInDegrees,'angleInDegrees',0);
adigatorOutputs = {angleInDegrees};
[adigatorFunInfo, adigatorOutputs] = adigatorFunctionEnd(4,adigatorFunInfo,adigatorOutputs);
