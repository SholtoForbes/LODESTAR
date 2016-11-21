function [Drag, Alpha, flapdeflection,lift_search] = OutForce(theta,M,q,m,scattered,v,V,thetadot,time, lift_search)
%THIS IS THE SLOWEST PART OF THE ROUTINE

% v_H = v.*cos(theta);
% % 
% % % find aerodynamics using only gravity of vehicle
% gravity = m.*(- 6.674e-11.*5.97e24./(V + 6371e3).^2 + v_H.^2./(V + 6371e3)); %Includes Gravity Variation and Centripetal Force 
% 
% lift_search = -gravity.*cos(theta);


%use LiftForceInterp splines
Alpha = scattered.AoA(v,V,lift_search);
flapdeflection = scattered.flapdeflection(v,V,lift_search);
% flapdeflection = 1;
Drag = scattered.drag(v,V,lift_search); 
% Flap_pitchingmoment = scattered.flap_pm(v,V,lift_search);

% omegadot = diff(thetadot)./diff(time);
% I = 150000; % from hand calculation
% extramoment = [0 -omegadot*I];
% 
% flapdeflection = scattered.flap_def(M,Alpha,Flap_pitchingmoment + extramoment);
end






