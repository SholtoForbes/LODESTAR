clear all
mat = [];




% u = [2850:25:2925];
% u = 2900
u = [2900:25:2975];

% options.Display = 'iter';
options.Display = 'final';
options.Algorithm = 'sqp';
options.TolFun = 1e-4;
options.TolX = 1e-4;
% for phi0 = [-0.1271-0.005 -0.1271 -0.1271+0.005]
    phi0 = -0.13 % this has very minimal effect
% for zeta0 = [1.70 1.7040 1.7080]
zeta0 = 1.69 % this is the phi to reach close to 1.704 rad heading angle (SSO)


%     guess = [2500/10000  deg2rad(20) deg2rad(20) deg2rad(20) deg2rad(20) 200/1000;
%         2500/10000  deg2rad(20) deg2rad(20) deg2rad(20) deg2rad(20) 200/1000;
%         2500/10000  deg2rad(20) deg2rad(20) deg2rad(20) deg2rad(20) 200/1000;
%         2500/10000  deg2rad(20) deg2rad(20) deg2rad(20) deg2rad(20) 200/1000];
%    guess = []; 
    
%     phi0 = -0.1271
%     zeta0 = 1.7011
%     
%     guess = [1500; deg2rad(13); deg2rad(13)]*ones(1,length(u));
%     guess = [1500; deg2rad(20); deg2rad(20)]*ones(1,length(u));
% phi0 = -0.13154;
% zeta0 = deg2rad(96.9);


% options.TypicalX = [2600 0.2 0.2 0.2 0.2 250];
% for k = [30000:1000:35000 35000:250:38000 38500:500:40000]
for k = [33000:1000:38000]
    for j = [0:0.025:0.05]
        
        temp_guess_no = 1;
        
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
     [AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA,q,gamma,D,AoA_max] = ThirdStageSim([0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 10],k,j,u(i), phi0, zeta0);
%         guess = [1700 AoA_max-0.01 AoA_max/2+AoA_max/2*(0.1-j)/0.1-0.01];
%         guess = [1600 AoA_max(i)-0.01];
% guess = [2700  deg2rad(7.5) deg2rad(10)];
% guess = [2750  (deg2rad(14)+(deg2rad(6)-deg2rad(14))*j/0.05) (deg2rad(10)+((AoA_max-0.01)-deg2rad(10))*j/0.05) (deg2rad(14)+((AoA_max-0.01)-deg2rad(14))*j/0.05) (deg2rad(12)+(deg2rad(8)-deg2rad(12))*j/0.05)]
% if k >= 34000 % variable guess
% x0 = [2500/10000  AoA_max AoA_max AoA_max AoA_max 200/1000]
% else
%     x0 = [2500/10000  0 AoA_max AoA_max AoA_max 200/1000]
% end

% if temp_guess_no == 1;
% x0 = [2500/10000  AoA_max AoA_max AoA_max AoA_max 200/1000];
% elseif temp_guess_no > 1
%    x0 =  guess(i,:);
% end
% x0 = [2590/10000  AoA_max*ones(1,16) 250/1000];
nodesalt = [33000; 33000; 34000; 36000 ;36000];
nodesgam = [0;0.05; 0; 0; 0.05];
vals =  [deg2rad(2);deg2rad(.5); deg2rad(2); 0 ;0];
interp = scatteredInterpolant(nodesalt,nodesgam,vals);
x0 = [2590/10000  AoA_max*ones(1,16)-interp(k,j) 250/1000] % this problem is extremely sensitive to initital guess! mostly at low altitude low gamma


%         x = fminsearch(@(x)Payload(x,k,j,u(i),phi0,zeta0),x0,options);
% x = fmincon(@(x)Payload(x,k,j,u(i), phi0, zeta0),x0,[],[],[],[],[1400 0 0 0 0 50],[3200 AoA_max AoA_max AoA_max AoA_max 250],@(x)Constraint(x,k,j,u(i), phi0, zeta0),options);
x0;
% x = fmincon(@(x)Payload(x,k,j,u(i), phi0, zeta0),x0,[],[],[],[],[2200/10000 0 0 0 0 50/1000],[3000/10000 AoA_max AoA_max AoA_max AoA_max 300/1000],@(x)Constraint(x,k,j,u(i), phi0, zeta0),options);
x = fmincon(@(x)Payload(x,k,j,u(i), phi0, zeta0),x0,[],[],[],[],[2200/10000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 200/1000],[3000/10000 AoA_max AoA_max AoA_max AoA_max AoA_max AoA_max AoA_max AoA_max AoA_max AoA_max AoA_max AoA_max AoA_max AoA_max AoA_max AoA_max 270/1000],@(x)Constraint(x,k,j,u(i), phi0, zeta0),options);

        mfuel_burn = x(1);

        [AltF(i), vF(i), Alt, v, t, mpayload(i), Alpha, m,AoA] = ThirdStageSim(x,k,j,u(i),phi0,zeta0);
        
        temp(i,:) = x;
        temp_payload(i) = mpayload(i);
        u(i)
        mpayload(i)
        end
temp_payload
%         guess = temp;
        
%         mat = [mat;[phi0*ones(length(u),1),zeta0*ones(length(u),1),k*ones(length(u),1),j*ones(length(u),1),u.',mpayload.',temp(1,:).',temp(2,:).',temp(3,:).']];
mat = [mat;[phi0*ones(length(u),1),zeta0*ones(length(u),1),k*ones(length(u),1),j*ones(length(u),1),u.',temp_payload.']];

temp_guess_no = temp_guess_no + 1;
    end
end
% end
% end
dlmwrite('thirdstagenew.dat', mat,'delimiter','\t')