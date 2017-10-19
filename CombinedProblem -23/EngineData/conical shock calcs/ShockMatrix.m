% Calculates conditions after shock for a variety of flight conditions
% Saves as file

hangle = 5; % cone half angle (deg)

mat = [];
for i = 5:.1:10
    for j = 0:.1:8
        delete('input')
        delete('output')
        
        dlmwrite('input',[i hangle j], 'delimiter',' ');
        
        system('Shock.exe')
        pause(5)
        
        temp = dlmread('output');
        
        mat(end+1,:) = [i j temp];
    end
end

dlmwrite('ShockMat',mat, 'delimiter',' ')