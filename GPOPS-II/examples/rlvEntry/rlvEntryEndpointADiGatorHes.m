% This code was generated using ADiGator version 1.4
% ©2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function output = rlvEntryEndpointADiGatorHes(input)
global ADiGator_rlvEntryEndpointADiGatorHes
if isempty(ADiGator_rlvEntryEndpointADiGatorHes); ADiGator_LoadData(); end
Gator1Data = ADiGator_rlvEntryEndpointADiGatorHes.rlvEntryEndpointADiGatorHes.Gator1Data;
Gator2Data = ADiGator_rlvEntryEndpointADiGatorHes.rlvEntryEndpointADiGatorHes.Gator2Data;
% ADiGator Start Derivative Computations
%User Line: % Inputs
%User Line: % input.phase(phasenumber).initialstate -- row
%User Line: % input.phase(phasenumber).finalstate -- row
%User Line: % input.phase(phasenumber).initialtime -- scalar
%User Line: % input.phase(phasenumber).finaltime -- scalar
%User Line: % input.phase(phasenumber).integral -- row
%User Line: %
%User Line: % input.parameter -- row
%User Line: % input.auxdata = auxiliary information
%User Line: % Output
%User Line: % output.objective -- scalar
%User Line: % output.eventgroup(eventnumber).event -- row
latf.dv = input.phase.finalstate.dv(3);
latf.f = input.phase.finalstate.f(3);
%User Line: latf = input.phase.finalstate(3);
%User Line: % cost
output.objective.dv = uminus(latf.dv);
output.objective.f = uminus(latf.f);
%User Line: output.objective = -latf;
output.objective.dv_size = 14;
output.objective.dv_location = Gator1Data.Index1;
end


function ADiGator_LoadData()
global ADiGator_rlvEntryEndpointADiGatorHes
ADiGator_rlvEntryEndpointADiGatorHes = load('rlvEntryEndpointADiGatorHes.mat');
return
end