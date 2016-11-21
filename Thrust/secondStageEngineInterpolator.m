%2ndStageEngineInterpolator

% this runs the file EngineData.exe which takes input of Mach, T and p (behind initial shock, given in communicator matrix) and
% produces output of thrust force (N) and fuel usage (kg/s).

delete input
delete output
delete engineoutput_matrix

com = dlmread('communicator_extrapolate_MpT.txt'); % import data from atmosphere matrix

wcap = 0.65;

for i = 6:1.:13
    M = i ;

    for j = 0.:1.:6
        AoA = j;

        Mshock = griddata(com(:,1), com(:,2), com(:,4), M, AoA); %calculates M after shock from extrapolated hypaero communicator

        pshock = griddata(com(:,1), com(:,2), com(:,5), M, AoA);

        Tshock = griddata(com(:,1), com(:,2), com(:,6), M, AoA);
        
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
        engineoutput_matrix = fopen('engineoutput_matrix','a+');

        outputstring = [num2str(M,'%10.4e') ' ' num2str(AoA,'%10.4e') ' ' num2str(output_temp(1),'%10.4e') ' ' num2str(output_temp(2),'%10.4e') '\r\n'];

        fprintf(engineoutput_matrix,outputstring);

        fclose('all');
    end
end
