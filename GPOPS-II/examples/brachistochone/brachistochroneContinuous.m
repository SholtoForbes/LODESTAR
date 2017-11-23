%---------------------------------------------%
% BEGIN: function brachistochroneContinuous.m %
%---------------------------------------------%
function phaseout = brachistochroneContinuous(input)

g                 = input.auxdata.g;

t                 = input.phase.time;
x                 = input.phase.state(:,1);
y                 = input.phase.state(:,2);
v                 = input.phase.state(:,3);
u                 = input.phase.control;
xdot              = v.*cos(u);
ydot              = v.*sin(u);
vdot              = -g*sin(u);
phaseout.dynamics = [xdot, ydot, vdot];

%---------------------------------------------%
% END: function brachistochroneContinuous.m   %
%---------------------------------------------%
