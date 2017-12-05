function path = ThirdStagePath(primal)

global Vec_angle


Vec_angle(end+1) = 0;
Vec_angle_constraint = Vec_angle;


path = Vec_angle_constraint;

end
% path = DynamicPressure ;




