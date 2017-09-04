% LODESTAR Runner
clear all

str = {'Create Payload Matrix' 'Create LIFTARRAY' 'Max Payload Optimisation' 'Constant q' '1st Stage Optimisation' '3rd Stage Optimisation'};
[Selection,ok] = listdlg('PromptString','What Would You Like To Do?','SelectionMode','single','ListString',str);


if Selection == 1
    %% Payload Matrix Generation
    addpath('.\ThirdStage');
    run('PayloadMatrix.m')
elseif Selection == 2
    %% LIFTARRAY Generation
    addpath('.\SecondStage'); 
    run('LiftForceMAT.m') % CURRENTLY CRASHES WHEN RUN FROM HERE
elseif Selection == 3
    %% Max Payload Optimisation
    % Includes min fuel first stage optimisation
    addpath('.\SecondStage'); 
    [configfile,path] = uigetfile('*.m','Select Config File');
    run(configfile);
  
elseif Selection == 4
    %% Constant Dynamic Pressure Simulation
    % Includes min fuel first stage optimisation
    
elseif Selection == 5
    %% Separate First Stage Optimisation
    
elseif Selection == 6
    %% Separate Third Stage Optimisation
    
end