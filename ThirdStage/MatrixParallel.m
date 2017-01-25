clear all
mat = [];



% u = [2500:100:2800 2800:20:2900 3000];
u = [2700:25:2900 3000];
% u = [2880:20:2940];


% options.Display = 'iter';
options.Display = 'final';
options.TolFun = 1;
options.TolX = 10;
for phi0 = [-deg2rad(7.5) -deg2rad(7)]
for zeta0 = [deg2rad(96) deg2rad(98)]
    
%     phi0 = -deg2rad(7)
%     zeta0 = deg2rad(96)
%     
%     guess = [1500; deg2rad(13); deg2rad(13)]*ones(1,length(u));
%     guess = [1500; deg2rad(20); deg2rad(20)]*ones(1,length(u));
% phi0 = -0.13154;
% zeta0 = deg2rad(96.9);

% for k = [30000:1000:35000 35000:250:38000 38500:500:40000]
for k = [30000:500:40000]
    for j = [0:0.01:0.05]
        
        
        
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
        
        guess = [1500 AoA_max-0.01 AoA_max-0.01];
        
%         x0 = guess(:,i); 
        x0 = guess; 
        x = fminsearch(@(x)Payload(x,k,j,u(i),phi0,zeta0),x0,options);

        mfuel_burn = x(1);

        [AltF(i), vF(i), Alt, v, t, mpayload(i), Alpha, m,AoA(i)] = ThirdStageSim(x,k,j,u(i),phi0,zeta0);
        
        temp(:,i) = x;
        
        end
        
%         guess = temp;
        
%         mat = [mat;[phi0*ones(length(u),1),zeta0*ones(length(u),1),k*ones(length(u),1),j*ones(length(u),1),u.',mpayload.',temp(1,:).',temp(2,:).',temp(3,:).']];
mat = [mat;[phi0*ones(length(u),1),zeta0*ones(length(u),1),k*ones(length(u),1),j*ones(length(u),1),u.',mpayload.',temp(1,:).',temp(2,:).',]];
    end
end
end
end
dlmwrite('thirdstagenew.dat', mat,'delimiter','\t')