clear all
mat = [];



u = [2500:100:2800 2810:10:2900 3000];



% options.Display = 'iter';
options.Display = 'final';
options.TolFun = 0.1;
options.TolX = 0.1;
for phi0 = [-deg2rad(8.5):deg2rad(1):-deg2rad(6.5)]
for zeta0 = [deg2rad(96):deg2rad(1):deg2rad(98)]
    
    guess = [1500; deg2rad(13); deg2rad(13)]*ones(1,length(u));
    
% phi0 = -0.13154;
% zeta0 = deg2rad(96.9);
% for k = [30000:1000:34000 34100:100:36000 37000:1000:40000]
for k = [36000:100:37000]
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

        x0 = guess(:,i); 


        x = fminsearch(@(x)Payload(x,k,j,u(i),phi0,zeta0),x0,options);

        mfuel_burn = x(1);

        [AltF(i), vF(i), Alt, v, t, mpayload(i), Alpha, m,AoA(i)] = ThirdStageSim(x,k,j,u(i),phi0,zeta0);
        
        temp(:,i) = x;
        

        end
        
        guess = temp;
        
        mat = [mat;[phi0*ones(length(u),1),zeta0*ones(length(u),1),k*ones(length(u),1),j*ones(length(u),1),u.',mpayload.',AoA.']];
    end
end
end
end
dlmwrite('thirdstagenew.dat', mat,'delimiter','\t')