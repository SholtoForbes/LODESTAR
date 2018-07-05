% Vehicle Config File For LODESTAR

%% First Stage

% Sea Level Thrust
Stage1.T = 555900; % N

% Sea Level Specific Impulse
Stage1.Isp = 275; % s

% Nozzle Exit Area
Stage1.Anoz = 0.5518; % m^2

% Throttle
Stage1.Throttle = 0.85; % To throttle down the rocket so that it can pitch over more easily.

% Scaled Mass 
Stage1.m = 21816; % kg. This is the mass that the Falcon 1e would be if scaled down. Note that this gives the maximum possible fuel mass, which will be minimised.

% Fuel Mass Fraction
Stage1.FMF = 0.939; 

% Engine Mass
Stage1.mEngine = 470; % kg. Mass of Merlin 1C


%% Second Stage

% Structural Mass
Stage2.mStruct = 4910.5 - 132.8 + 179.41;% kg. Mass of everything but fuel from dawids work.

% Fuel Mass
Stage2.mFuel = 1562; % kg. Fuel tank mass scaled by surface area to hold 1562kg fuel, see Fuel Tank Sizing.txt

% Reference Area
Stage2.Aref = 62.77; % m^2


%% Third Stage

% Total Mass
Stage3.mTot = 3300;

% Heat Shield Mass
Stage3.mHS = 130.9; %kg 

% Engine Mass
Stage3.mEng = 52; %kg. Kestrel

% Structural Mass Fraction (No Heat Shield)
Stage3.SMF = 0.09;

% Specific Impulse
Stage3.Isp = 317; %Kestrel, from Falcon 1 users guide

% Mass Flow Rate
Stage3.mdot = 9.86977*1.5; %Kestrel Modified

% Efficiency
Stage3.Eff = 0.98; % efficiency reduction due to mass flow rate increase. 

% Reference Area
Stage3.Aref = 0.95; % diameter of 1.1m