%2ndStageEngineInterpolator
clear all
% this runs the file EngineData.exe which takes input of Mach, T and p (behind initial shock, given in communicator matrix) and
% produces output of thrust force (N) and fuel usage (kg/s).

delete input
delete output
delete engineoutput_matrixnew
delete grid.Isp_eng

wcap = 0.65;

data = dlmread('ENGINEDATA.txt');

M_englist = unique(sort(data(:,1))); % create unique list of Mach numbers from engine data
M_eng_interp = unique(sort(data(:,1)));

T_englist = unique(sort(data(:,2))); % create unique list of angle of attack numbers from engine data 
T_eng_interp = unique(sort(data(:,2)));

[grid.Mgrid_eng,grid.T_eng] =  ndgrid(M_eng_interp,T_eng_interp);

% scat = scatteredInterpolant(data(:,1),data(:,2),data(:,3))
for i = 1:30
    i
    for j= 1:30
delete output



    Mshock = grid.Mgrid_eng(i,j) ;


            Tshock=grid.T_eng(i,j);
   
            

        pshock = 1; % doesnt matter


        %write input file
        input = fopen('input','w');

        inputstring = [num2str(Mshock,'%10.4e') ' ' num2str(Tshock,'%10.4e') ' ' num2str(pshock,'%10.4e') ' ' num2str(wcap,'%10.4e')];

        fprintf(input,inputstring);

        fclose('all');

        %run program, this exe was created in fortran silverfrost, takes
        %'input', gives 'output'
        system(['engineinterpolator.exe']);
        % !EngineData.exe

        
        %this reads the output for thrust and fuel use 
        if exist('output')
        output_temp = dlmread('output');
        else
            output_temp = 0;
        end
            

        %write to file
        engineoutput_matrix = fopen('engineoutput_matrixnew','a+');

        outputstring = [num2str(Mshock,'%10.4e') ' ' num2str(Tshock,'%10.4e') ' ' num2str(output_temp(1),'%10.4e') '\r\n'];

        gridIsp_eng(i,j)  = output_temp(1);
        fprintf(engineoutput_matrix,outputstring);

        fclose('all');

    end
end
save('gridIsp_eng')