function [Drag, Alpha, flapdeflection,lift_search] = OutForce(theta,M,q,m,scattered,v,V,thetadot,time, lift_search,c,rho,Atmosphere)
%THIS IS THE SLOWEST PART OF THE ROUTINE

% v_H = v.*cos(theta);
% % 
% % % find aerodynamics using only gravity of vehicle
% gravity = m.*(- 6.674e-11.*5.97e24./(V + 6371e3).^2 + v_H.^2./(V + 6371e3)); %Includes Gravity Variation and Centripetal Force 
% 
% lift_search = -gravity.*cos(theta);


%use LiftForceInterp splines

Alpha = scattered.AoA(v,V,lift_search);
% flapdeflection = scattered.flapdeflection(v,V,lift_search);
Drag = scattered.drag(v,V,lift_search); 


Flap_pitchingmoment = scattered.flap_pm(v,V,lift_search);

omegadot = diff(thetadot)./diff(time);
alp0hadot = diff(Alpha)./diff(time);
I = 267400; % from Creo
extramoment = [0 (omegadot+alphadot)*I];

flapdeflection = scattered.flap_def(M,Alpha,Flap_pitchingmoment + extramoment);

% flapdeflection = scattered.flap_def(M,Alpha,Flap_pitchingmoment);

%========================================================
% 
%             SPARTAN_SCALE = 1;
%             A = 62.77;
% %             c = spline( Atmosphere(:,1),  Atmosphere(:,5), alt); % Calculate speed of sound using atmospheric data
% % 
% %             rho = spline( Atmosphere(:,1),  Atmosphere(:,4), alt); % Calculate density using atmospheric data
% 
%             q = 0.5 * rho .* (v .^2); % Calculating Dynamic Pressure
% 
%             M = v./c; % Calculating Mach No (Descaled)
%             
%             %% Calculate Thrust Component ==================================
%             
%             %% Determine AoA ==============================================================
%             T0 = spline( Atmosphere(:,1),  Atmosphere(:,2), V); 
%             P0 = spline( Atmosphere(:,1),  Atmosphere(:,3), V); 
% 
% %             alphaguess = 5;
% %             error(i) = 2;
% %            while error(i) > 1 && alphaguess < 7
% for i = 1:length(V)
% %     i
% %     M(i), Alpha(i), scattered, SPARTAN_SCALE,scattered.pitchingmoment_spline,scattered.flaplift_spline,scattered.Cl_spline,q(i),A,lift_search(i),T0(i),P0(i)
%             Alpha(i) = fminsearch(@(Alpha)LiftError(M(i), Alpha, scattered, SPARTAN_SCALE,scattered.pitchingmoment_spline,scattered.flaplift_spline,scattered.Cl_spline,q(i),A,lift_search(i),T0(i),P0(i)),5);
% end
%             %Fuel Cost ===========================================================================
% %             M, Alpha, scattered, SPARTAN_SCALE,T0,P0
%             [Isp,Fueldt] = RESTM12int(M, Alpha, scattered, SPARTAN_SCALE,T0,P0);
% 
%             Thrust = Isp.*Fueldt*9.81;
% 
%             %======================================================================
% 
%             Cl1 = scattered.Cl_spline(M,Alpha);
% 
%             body_pitchingmoment = scattered.pitchingmoment_spline(M, Alpha);% first approximation of pitchingmoment using only body lift
% 
%             Flap_lift = q./50000.*scattered.flaplift_spline(M,Alpha,-body_pitchingmoment).*SPARTAN_SCALE^(2/3);% first approximation of flap lift, scale is only applied here as it will cancel for pitchingmoments
% 
%             total_lift = Cl1.*A.*q + Flap_lift + Thrust.*sin(deg2rad(Alpha)); %first total lift force, with normalised dynamic pressure, this needs to iterate to equal the original liftq
% 
% %             error(i) = abs(total_lift - Lift(i));
% %             alphaguess = alphaguess + 0.2;
% %            end
%             
% 
%             flapdeflection = scattered.flapdeflection_spline(M,Alpha,-body_pitchingmoment);
% 
%             Drag = scattered.Cd_spline(M,Alpha).*A.*q +  q/50000.*scattered.flapdrag_spline(M,Alpha,-body_pitchingmoment).*SPARTAN_SCALE^(2/3);
% 

end






