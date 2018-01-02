function xdot = FirstStageDynamics(primal)

global iteration
iteration = iteration + 1;

global Throttle
global interp
global q_forward
global xi
global iterative_V
global iterative_t
global iterative_V_f
global Vehicle
global Atmosphere

h = primal.states(1,:);		
v = primal.states(2,:);		
m = primal.states(3,:);	
gamma = primal.states(4,:);	
alpha = primal.states(5,:);	
zeta = primal.states(6,:);	

alphadot = primal.states(7,:);

phi = primal.states(8,:);	

t = primal.nodes;

x = [h; v ;m ;gamma; alpha; zeta; alphadot; phi];

u = primal.controls;
phase = 'postpitch';

[dz,q,xi] = rocketDynamics(x,u,t,phase,interp,Throttle,Vehicle,Atmosphere); % Pass primal variables into dynamics simulator

q_forward = q;

if rem(iteration,5000) == 0
    iterative_V_f(end+1,:) = cumtrapz(t,v.*sin(gamma));
    cla
    hold on
    todisp = ['Total Iterations Done:',num2str(iteration)];
    disp(todisp);
    figure(1010)
    iterative_V(end+1,:) = h;
    iterative_t(end+1,:) = t;
    if length(iterative_V(:,1)) > 4
        plot(iterative_t(end-4,:),iterative_V(end-4,:),'Color',[0.8 0.8 0.8])
    end
    
    if length(iterative_V(:,1)) > 2
        plot(iterative_t(end-2,:),iterative_V(end-2,:),'Color',[0.4 0.4 0.4])
    end
    
    if length(iterative_V(:,1)) > 3
        plot(iterative_t(end-3,:),iterative_V(end-3,:),'Color',[0.6 0.6 0.6])
    end
    plot(iterative_t(end-1,:),iterative_V(end-1,:),'Color',[0.2 0.2 0.2])
    plot(iterative_t(end,:),iterative_V(end,:),'Color',[0 0 0])
%     plot(iterative_t(end,:),iterative_V_f(end,:),'Color','r')

    filename = 'testnew51.gif';
    frame = getframe(1010);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,filename,'gif','WriteMode','append');

    drawnow
    
end

xdot = dz;

