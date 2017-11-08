function [Drag, Alpha, flapdeflection] = OutForce2(theta,M,q,m,scattered,v,V,thetadot,time,t_ratio, Efficiency, SPARTAN_SCALE,A,Lift)
for i = 1:length(time)
% parfor i = 1:length(time)
% Alpha(i) = fminsearch(@(Alpha)LiftError(M(i), Alpha, t_ratio(i), Efficiency(i), scattered, SPARTAN_SCALE,scattered.pitchingmoment_spline,scattered.flaplift_spline,scattered.Cl_spline,q(i),A,Lift(i)),5);
Alpha(i) = fminbnd(@(Alpha)LiftError(M(i), Alpha, t_ratio(i), Efficiency(i), scattered, SPARTAN_SCALE,scattered.pitchingmoment_spline,scattered.flaplift_spline,scattered.Cl_spline,q(i),A,Lift(i)),3,7);
end
LiftError(M(1), Alpha(1), t_ratio(1), Efficiency(1), scattered, SPARTAN_SCALE,scattered.pitchingmoment_spline,scattered.flaplift_spline,scattered.Cl_spline,q(1),A,Lift(1))
[Isp,Fueldt] = RESTM12int(M, Alpha, t_ratio, Efficiency, scattered, SPARTAN_SCALE);

Thrust = Isp.*Fueldt*9.81;

%======================================================================

Cl1 = scattered.Cl_spline(M,Alpha);

body_pitchingmoment = scattered.pitchingmoment_spline(M, Alpha);% first approximation of pitchingmoment using only body lift

Flap_lift = q./50000.*scattered.flaplift_spline(M,Alpha,-body_pitchingmoment)*SPARTAN_SCALE^(2/3);% first approximation of flap lift, scale is only applied here as it will cancel for pitchingmoments

total_lift = Cl1*A.*q + Flap_lift + Thrust.*sin(deg2rad(Alpha)); %first total lift force, with normalised dynamic pressure, this needs to iterate to equal the original liftq

error = abs(total_lift - Lift);

flapdeflection = scattered.flapdeflection_spline(M,Alpha,-body_pitchingmoment);

Drag = scattered.Cd_spline(M,Alpha)*A.*q +  q/50000.*scattered.flapdrag_spline(M,Alpha,-body_pitchingmoment)*SPARTAN_SCALE^(2/3);
end