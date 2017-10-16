%2ndStageEngineInterpolator
clear all
% this runs the file EngineData.exe which takes input of Mach, T and p (behind initial shock, given in communicator matrix) and
% produces output of thrust force (N) and fuel usage (kg/s).

delete input
delete output
delete engineoutput_matrix

com = dlmread('communicator_extrapolate_MpT.txt'); % import data from atmosphere matrix

wcap = 0.65;

for i = 5:.2:6
    Mshock = i ;

    for j = 250.:10:300
            Tshock=j;

        pshock = 1; % doesnt matter


        %write input file
        input = fopen('input','w');

        inputstring = [num2str(Mshock,'%10.4e') ' ' num2str(Tshock,'%10.4e') ' ' num2str(pshock,'%10.4e') ' ' num2str(wcap,'%10.4e')];

        fprintf(input,inputstring);

        fclose('all');

        %run program, this exe was created in fortran silverfrost, takes
        %'input', gives 'output'
        system(['EngineData.exe']);
        % !EngineData.exe

        
        %this reads the output for thrust and fuel use 
        output_temp = dlmread('output');

        %write to file
        engineoutput_matrix = fopen('engineoutput_matrixnew','a+');

        outputstring = [num2str(M,'%10.4e') ' ' num2str(AoA,'%10.4e') ' ' num2str(t_ratio,'%10.4e') ' ' num2str(output_temp(1),'%10.4e') ' ' num2str(output_temp(2),'%10.4e') '\r\n'];

        fprintf(engineoutput_matrix,outputstring);

        fclose('all');
        end
    end
end
