% Optimisation Config


%% First Stage

% Initial and End Conditions

% Bounds

% Number of Nodes

% Initial Guesses

%% Second Stage %%=========================================================

% Number of Nodes %--------------------------------------------------------


% Bounds %--------------------------------------------------------
bounds.lower.controls = [-0.00003];
bounds.upper.controls = [0.0002];

% Initial Guesses %--------------------------------------------------------
guess.states(1,:) = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000/v0^2)-100,33500];
guess.states(2,:) = [v0, 2900];
guess.states(3,:) = [0.0,deg2rad(4.5)];
guess.states(4,:) = [mfuelU, 0];
guess.states(5,:) = [0,0];
guess.states(6,:) = [1.682,1.699];
guess.controls(1,:) = [0,0]; 
guess.time = [t0 ,360];

%% Third Stage


%% Second Stage Return