clear all

addpath('..\')
Atmosphere = dlmread('atmosphere.txt');
run VehicleConfig.m
run TrajectoryConfig50kPa.m

CG_SPARTAN = 14.5; % m from nose. Specified for good aerodynamics.
CG_ThirdStage = 18.04;% m from nose. calculated using Creo.
CG_ConicalFuel = 9.33;% m from nose. calculated using Creo.
CG_CylindricalFuel = 16.55; % m from nose. calculated using Creo.

Mass_SPARTAN = Stage2.mStruct;
Mass_ThirdStage = Stage3.mTot;
Mass_ConicalFuel = 852;
Mass_CylindricalFuel = 353*2;

mTot = Mass_SPARTAN + Mass_ThirdStage + Mass_ConicalFuel + Mass_CylindricalFuel;

CG = CG_SPARTAN*Mass_SPARTAN/mTot + CG_ThirdStage*Mass_ThirdStage/mTot + CG_ConicalFuel*Mass_ConicalFuel/mTot + CG_CylindricalFuel*Mass_CylindricalFuel/mTot;

