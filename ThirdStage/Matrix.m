clear all
temp = [];

guess = [1500 deg2rad(13)];

for k = [30000:1000:34000 34100:100:36000 37000:1000:40000]
    for j = 0.0:0.01:.05
        for u = [2500:100:2800 2810:10:2900 3000]

        mScale = 1; % This needs to be manually changed in altitude and velocity files as well
        x0 = guess; % NOTE, NOT USING AoA AT THE MOMENT
        options.Display = 'iter';
        options.TolFun = 0.1;
        options.TolX = 0.1;

        x = fminsearch(@(x)Payload(x,k,j,u),x0,options);

        mfuel_burn = x(1)

        [AltF, vF, Alt, v, t, mpayload, Alpha, m,AoA] = ThirdStageSim(x,k,j,u);

        temp(end+1,:) = [k,j,u,mpayload,AoA];
        

        end
    end
end

dlmwrite('thirdstage.dat', temp,'delimiter','\t')