% Startup file for DIDO
% 
%---------------------------------------------------------------------------------------------
% This file is executed only when matlab is started from the same folder as where this file is
% located.
%---------------------------------------------------------------------------------------------
%


%================
% Initialize DIDO
%================
current_dir = pwd;
                                    %**************************
warning off MATLAB:oldPfileVersion  % DO NOT DELETE THIS LINE
InitializeDido;                     % DO NOT DELETE THIS LINE
                                    %**************************
%================


% You may alter or delete the following lines:

% make thicker lines for plotting:
set(0, 'defaultaxesfontsize', 14, 'defaultaxeslinewidth', .7, ...
    'defaultlinelinewidth', .9, 'defaultpatchlinewidth', .7, ...
    'DefaultTextFontSize', 14);