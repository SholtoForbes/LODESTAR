% Optimisation Config


%% First Stage
% DIDO Inputs
% Initial Conditions %-----------------------------------------------------

% End Conditions %---------------------------------------------------------

% Number of Nodes %--------------------------------------------------------

% Bounds %-----------------------------------------------------------------
% Define the search space of the optimiser. 

% Variable Constraints %---------------------------------------------------
% Define variable Bounds.

% Initial Guesses %--------------------------------------------------------
% Kickstart the optimiser.

%% Second Stage %%=========================================================
% DIDO Inputs
% Initial Conditions %-----------------------------------------------------
Stage2.v0               = 1520; % Initial Velocity, m/s
Stage2.mFuel0           = 1562; % Initial Fuel Mass, kg
Stage2.phi0             = -0.264; % Initial Latitude, rad
Stage2.xi0              = 0; % Initial Longitude, rad

% End Conditions %---------------------------------------------------------
Stage2.mFuelf           = 0; % End Fuel Mass, kg
Stage2.Zetaf            = 1.78; % End Heading Angle, rad
Stage2.gammaf           = [deg2rad(2) deg2rad(4)]; % End Trajectory Angle

% Number of Nodes %--------------------------------------------------------
Stage2.Nodes = 104;

% Bounds %-----------------------------------------------------------------
Stage2.Bounds.Alt       = [20000 50000]; % Altitude Bounds, m
Stage2.Bounds.v         = [1500 3100]; % Velocity Bounds, m/s
Stage2.Bounds.gamma     = [-0.1 deg2rad(5)]; % Trajectory Angl Bounds, rad
Stage2.Bounds.mFuel     = [0 Stage2.mFuel0]; % Fuel Mass Bounds, kg
Stage2.Bounds.gammadot  = [-0.001 0.002]; % Trajectory Angle Derivative Bounds, rad/s
Stage2.Bounds.zeta      = [1 2]; % Heading Angle Bounds, rad
Stage2.Bounds.control   = [-0.00003 0.0002]; % Trajectory Angle Double-Derivative Bounds, rad/s^2
Stage2.Bounds.time      = [100 800]; % Time Bounds, s

% Variable Constraints %---------------------------------------------------
Stage2.Bounds.q         = [0 50000]; % Dynamic Pressure Bounds, pa 
Stage2.Bounds.AoA       = [0 9]; % Angle of Attack Bounds, deg

% Initial Guesses %--------------------------------------------------------
Stage2.Guess.Alt        = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000/v0^2)-100,33500]; % Altitude Guess, m
Stage2.Guess.v          = [v0, 2900]; % Velocity Guess, m/s
Stage2.Guess.gamma      = [0.0,deg2rad(4.5)]; % Trajectory Angle guess, rad
Stage2.Guess.mFuel      = [mFuelU, 0]; % Fuel Mass Guess, kg
Stage2.Guess.gammadot   = [0,0]; % Trajectory Angle Derivative Guess, rad/s
Stage2.Guess.zeta       = [1.682,1.699]; % Heading Angle Guess, rad
Stage2.Guess.control    = [0,0]; % Trajectory Angle Double-Derivative Guess, rad/s^2
Stage2.Guess.time       = [t0 ,360]; % Time Guess, s

%% Third Stage


%% Second Stage Return
% DIDO Inputs
% Initial Conditions %-----------------------------------------------------

% End Conditions %---------------------------------------------------------

% Number of Nodes %--------------------------------------------------------

% Bounds %-----------------------------------------------------------------

% Variable Constraints %---------------------------------------------------

% Initial Guesses %--------------------------------------------------------