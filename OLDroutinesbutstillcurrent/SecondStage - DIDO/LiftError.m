function error = LiftError(M, x, scattered, SPARTAN_SCALE,pitchingmoment_spline,flaplift_spline,Cl_spline,q,A,Lift,T0,P0)
Alpha = x(1);
flapdeflection = x(2);

            %Fuel Cost ===========================================================================
            
            [Isp,Fueldt] = RESTM12int(M, Alpha, scattered, SPARTAN_SCALE,T0,P0);

            Thrust = Isp.*Fueldt*9.81;

            %======================================================================

            Cl1 = Cl_spline(M,Alpha);

            body_pitchingmoment = pitchingmoment_spline(M, Alpha);% first approximation of pitchingmoment using only body lift
flap_pitchingmoment = scattered.flappm_spline(M, Alpha,flapdeflection);
%             Flap_lift = q./50000*flaplift_spline(M,Alpha,-body_pitchingmoment)*SPARTAN_SCALE^(2/3);% first approximation of flap lift, scale is only applied here as it will cancel for pitchingmoments
Flap_lift = q./50000*flaplift_spline(M,Alpha,flapdeflection);
            total_lift = Cl1*A*q + Flap_lift + Thrust*sin(deg2rad(Alpha)); %first total lift force, with normalised dynamic pressure, this needs to iterate to equal the original liftq

            error = (total_lift - Lift)^2 + (body_pitchingmoment + flap_pitchingmoment)^2;
% error = abs(total_lift - Lift)
            
end