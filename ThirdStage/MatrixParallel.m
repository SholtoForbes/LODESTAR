clear all
mat = [];




% u = [2700:25:2900 3000];
u = 2900


% options.Display = 'iter';
options.Display = 'final';
options.TolFun = 1;
options.TolX = 1;
% for phi0 = [-0.1271-0.005 -0.1271 -0.1271+0.005]
    phi0 = -0.13 % this has very minimal effect
% for zeta0 = [1.70 1.7040 1.7080]
zeta0 = 1.69 % this is the phi to reach close to 1.704 rad heading angle (SSO)
    
%     phi0 = -0.1271
%     zeta0 = 1.7011
%     
%     guess = [1500; deg2rad(13); deg2rad(13)]*ones(1,length(u));
%     guess = [1500; deg2rad(20); deg2rad(20)]*ones(1,length(u));
% phi0 = -0.13154;
% zeta0 = deg2rad(96.9);

% for k = [30000:1000:35000 35000:250:38000 38500:500:40000]
for k = [32000:1000:40000]
    for j = [0:0.02:0.1]
        
        
        
        phi0
        zeta0
        k
        j
        AltF = [];
        vF = [];
        Alt = [];
        t = [];
        mpayload = [];
        Alpha = [];
        m = [];
        AoA = [];
        
        parfor i = 1:length(u)
        [AltF_temp, vF_temp, Alt_temp, v_temp, t_temp, mpayload_temp, Alpha_temp, m_temp,AoA_temp,q_temp,gamma_temp,D_temp,AoA_max] = ThirdStageSim([0 0 0],k,j,u(i), phi0, zeta0);
        
%         guess = [1700 AoA_max-0.01 AoA_max/2+AoA_max/2*(0.1-j)/0.1-0.01];
%         guess = [1600 AoA_max(i)-0.01];
guess = [2200  deg2rad(10) deg2rad(10)];
%         x0 = guess(:,i); 
        x0 = guess; 
        x = fminsearch(@(x)Payload(x,k,j,u(i),phi0,zeta0),x0,options);

        mfuel_burn = x(1);

        [AltF(i), vF(i), Alt, v, t, mpayload(i), Alpha, m,AoA] = ThirdStageSim(x,k,j,u(i),phi0,zeta0);
        
        temp(:,i) = x;
        
        end
        
%         guess = temp;
        
%         mat = [mat;[phi0*ones(length(u),1),zeta0*ones(length(u),1),k*ones(length(u),1),j*ones(length(u),1),u.',mpayload.',temp(1,:).',temp(2,:).',temp(3,:).']];
mat = [mat;[phi0*ones(length(u),1),zeta0*ones(length(u),1),k*ones(length(u),1),j*ones(length(u),1),u.',mpayload.',temp(1,:).',temp(2,:).',]];
    end
end
% end
% end
dlmwrite('thirdstagenew.dat', mat,'delimiter','\t')