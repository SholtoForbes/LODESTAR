% Optimisation Config

% Bounding the state and control variables creates a 'search space'. The
% optimal solution is found within this search space.These bounds should be sufficiently wide to allow a wide
% search space, but tight enough that a solution can be found efficiently.  These bounds must be
% chosen carefully, with the physical model in mind. 

%% First Stage %%==========================================================
% DIDO Inputs
% Initial Conditions %-----------------------------------------------------
Stage1.Initial.Alt      = 90; % Initial Altitude, m (note that this is at pitchover)
Stage1.Initial.v        = 30; % Initial Velocity, m/s (note that this is at pitchover)
Stage1.Initial.gamma    = deg2rad(89.9); % Initial trajectory Angle, rad (note that this is at pitchover) 
Stage1.Initial.AoA      = 0; % Initial Angle of Attack, rad 

% End Conditions %---------------------------------------------------------
Stage1.End.mTot         = Stage1.m+Stage2.mStruct+Stage2.mFuel+Stage3.mTot-(Stage1.FMF*(Stage1.m-Stage1.mEngine)+Stage1.mEngine); % End Mass, kg (note that this includes the SPARTAN and third stage)

% Number of Nodes %--------------------------------------------------------
Stage1.Nodes            = 91;

% Bounds %-----------------------------------------------------------------
Stage1.Bounds.Alt       = [0, 30000]; % Altitude Bounds, m
Stage1.Bounds.v         = [0, 3000]; % Velocity Bounds, m/s
Stage1.Bounds.gamma     = [deg2rad(-.1), Stage1.Initial.gamma]; % Trajectory Angl Bounds, rad
Stage1.Bounds.mTot      = [Stage1.m+Stage2.mStruct+Stage2.mFuel+Stage3.mTot,...
                           Stage1.m+Stage2.mStruct+Stage2.mFuel+Stage3.mTot-(Stage1.FMF*(Stage1.m-Stage1.mEngine)+Stage1.mEngine)]; % Total Mass Bounds, kg (note that this includes the SPARTAN and third stage)
Stage1.Bounds.AoA       = [-deg2rad(5), deg2rad(2)]; % Angle of Attack Bounds, rad
Stage1.Bounds.zeta      = [0, 2*pi]; % Heading Angle Bounds, rad
Stage1.Bounds.AoAdot    = [-0.1, 0.1]; % Angle of Attack Derivative Bounds, rad/s
Stage1.Bounds.phi       = [-0.5, -0.2]; % Latitude Bounds, rad
Stage1.Bounds.control   = [-.0005, .0005]; % (Control) Angle of Attack Double-Derivative Bounds, rad/s^2
Stage1.Bounds.time      = [0 300]; % Time Bounds, s 

% Guesses are adaptively set 

%% Second Stage %%=========================================================
% Initial Conditions %-----------------------------------------------------
Stage2.Initial.v        = 1520; % Initial Velocity, m/s
Stage2.Initial.mFuel    = Stage2.mFuel; % Initial Fuel Mass, kg
Stage2.Initial.phi      = -0.264; % Initial Latitude, rad
Stage2.Initial.xi       = 0; % Initial Longitude, rad

% End Conditions %---------------------------------------------------------
% These can be either a set number, or a bound.
Stage2.End.mFuel        = 0; % End Fuel Mass, kg
Stage2.End.Zeta         = 1.78; % End Heading Angle, rad

Stage2.End.gammaOpt     = [deg2rad(1), deg2rad(6)]; % End Trajectory Angle Bounds, rad (for max payload trajectories)
Stage2.End.gammaConst   = deg2rad(1.5); % End Trajectory Angle Bound, rad (for constant dynamic pressure trajectories)

% Number of Nodes %--------------------------------------------------------
Stage2.Nodes            = 104;

% Bounds %-----------------------------------------------------------------
Stage2.Bounds.Alt       = [20000, 50000]; % Altitude Bounds, m
Stage2.Bounds.v         = [1400, 3100]; % Velocity Bounds, m/s
Stage2.Bounds.gamma     = [-0.5, deg2rad(6)]; % Trajectory Angle Bounds, rad
Stage2.Bounds.mFuel     = [0, Stage2.mFuel]; % Fuel Mass Bounds, kg
Stage2.Bounds.gammadot  = [-0.01, 0.02]; % Trajectory Angle Derivative Bounds, rad/s
Stage2.Bounds.zeta      = [-4/3*pi, 2*pi]; % Heading Angle Bounds, rad
Stage2.Bounds.control   = [-0.00003, 0.0002]; % (Control) Trajectory Angle Double-Derivative Bounds, rad/s^2
Stage2.Bounds.time      = [100, 800]; % Time Bounds, s

% Variable Constraints %---------------------------------------------------
Stage2.Bounds.q         = [0, 50000]; % Dynamic Pressure Bounds, pa 
Stage2.Bounds.AoA       = [0, 9]; % Angle of Attack Bounds, deg

% Initial Guesses %--------------------------------------------------------
Stage2.Guess.Alt        = [interp1(Atmosphere(:,4),Atmosphere(:,1),2*50000/Stage2.Initial.v^2)-100, 35100]; % Altitude Guess, m
Stage2.Guess.v          = [Stage2.Initial.v, 2900]; % Velocity Guess, m/s
Stage2.Guess.gamma      = [0.0, deg2rad(4.3)]; % Trajectory Angle guess, rad
Stage2.Guess.mFuel      = [Stage2.mFuel, 0]; % Fuel Mass Guess, kg
Stage2.Guess.gammadot   = [0, 0]; % Trajectory Angle Derivative Guess, rad/s
Stage2.Guess.zeta       = [1.682, 1.699]; % Heading Angle Guess, rad
Stage2.Guess.control    = [0, 0]; % Trajectory Angle Double-Derivative Guess, rad/s^2
Stage2.Guess.time       = [0, 430]; % Time Guess, s

%% Third Stage %%==========================================================
% Fmincon Inputs

% Bounds %-----------------------------------------------------------------
Stage3.Bounds.AoA       = 20; % Maximum Angle of Attack, deg

% Adaptive Constraints %---------------------------------------------------
Stage3.Bounds.Alt       = [100000 400000]; % End Altitude Bounds, m
Stage3.Bounds.VecAng    = 8; % Maximum Vector Angle, deg

%% Second Stage Return
% Initial Conditions %-----------------------------------------------------

% End Conditions %---------------------------------------------------------

% Number of Nodes %--------------------------------------------------------

% Bounds %-----------------------------------------------------------------

% Variable Constraints %---------------------------------------------------

% Initial Guesses %--------------------------------------------------------