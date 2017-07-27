function [Vec_angle_search] = AlphaFind(Alpha,v,rho,Drag_interp,Lift_interp,T,M,A)

CD = Drag_interp(M,rad2deg(Alpha));
    
    CL = Lift_interp(M,rad2deg(Alpha));

%     CA(i) = 0.346 + 0.183 - 0.058*M(i)^2 + 0.00382*M(i)^3;
%     
%     CN(i) = (5.006 - 0.519*M(i) + 0.031*M(i)^2)*rad2deg(Alpha(i));
    
    

    D = 1/2*rho*(v^2)*A*CD;
    L = 1/2*rho*(v^2)*A*CL; % Aerodynamic lift
    
    
    %% Thrust vectoring
        Vec_angle = asin(2.5287/2.9713*L/T); % calculate the thrust vector angle necessary to resist the lift force moment.

   
Vec_angle_search = (Vec_angle - deg2rad(25))^2;
end

