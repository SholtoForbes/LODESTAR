% This code was generated using ADiGator version 1.4
% �2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function phaseout = rlvEntryContinuousADiGatorHes(input)
global ADiGator_rlvEntryContinuousADiGatorHes
if isempty(ADiGator_rlvEntryContinuousADiGatorHes); ADiGator_LoadData(); end
Gator1Data = ADiGator_rlvEntryContinuousADiGatorHes.rlvEntryContinuousADiGatorHes.Gator1Data;
Gator2Data = ADiGator_rlvEntryContinuousADiGatorHes.rlvEntryContinuousADiGatorHes.Gator2Data;
% ADiGator Start Derivative Computations
%User Line: % ---------------------------------------------------%
%User Line: % ------ Extract Each Component of the State ------- %
%User Line: % ---------------------------------------------------%
rad.dV = input.phase.state.dV(:,1);
rad.f = input.phase.state.f(:,1);
%User Line: rad  = input.phase.state(:,1);
lon.dV = input.phase.state.dV(:,2);
lon.f = input.phase.state.f(:,2);
%User Line: lon  = input.phase.state(:,2);
lat.dV = input.phase.state.dV(:,3);
lat.f = input.phase.state.f(:,3);
%User Line: lat  = input.phase.state(:,3);
v.dV = input.phase.state.dV(:,4);
v.f = input.phase.state.f(:,4);
%User Line: v    = input.phase.state(:,4);
fpa.dV = input.phase.state.dV(:,5);
fpa.f = input.phase.state.f(:,5);
%User Line: fpa  = input.phase.state(:,5);
azi.dV = input.phase.state.dV(:,6);
azi.f = input.phase.state.f(:,6);
%User Line: azi  = input.phase.state(:,6);
aoa.dV = input.phase.control.dV(:,1);
aoa.f = input.phase.control.f(:,1);
%User Line: aoa  = input.phase.control(:,1);
bank.dV = input.phase.control.dV(:,2);
bank.f = input.phase.control.f(:,2);
%User Line: bank = input.phase.control(:,2);
%User Line: % ---------------------------------------------------%
%User Line: % ------- Compute the Aerodynamic Quantities --------%
%User Line: % ---------------------------------------------------%
cd0.f = input.auxdata.cd(1);
%User Line: cd0      = input.auxdata.cd(1);
cd1.f = input.auxdata.cd(2);
%User Line: cd1      = input.auxdata.cd(2);
cd2.f = input.auxdata.cd(3);
%User Line: cd2      = input.auxdata.cd(3);
cl0.f = input.auxdata.cl(1);
%User Line: cl0      = input.auxdata.cl(1);
cl1.f = input.auxdata.cl(2);
%User Line: cl1      = input.auxdata.cl(2);
mu = input.auxdata.mu;
%User Line: mu       = input.auxdata.mu;
rho0 = input.auxdata.rho0;
%User Line: rho0     = input.auxdata.rho0;
H = input.auxdata.H;
%User Line: H        = input.auxdata.H;
S = input.auxdata.S;
%User Line: S        = input.auxdata.S;
mass = input.auxdata.mass;
%User Line: mass     = input.auxdata.mass;
altitude.dV = rad.dV;
altitude.f = rad.f - input.auxdata.Re;
%User Line: altitude = rad - input.auxdata.Re;
cada1f1dV = cd1.f*aoa.dV;
cada1f1 = cd1.f*aoa.f;
cada1f2dV = cada1f1dV;
cada1f2 = cd0.f + cada1f1;
cada2f1dV = 1.*aoa.f.^(1-1).*aoa.dV;
cada2f1 = aoa.f.^1;
cada2f2dV = 2.*cada2f1dV;
cada2f2 = 2*cada2f1;
cada1f3dVdV = aoa.dV.*cada2f2dV;
cada1f3dV = cada2f2.*aoa.dV;
cada1f3 = aoa.f.^2;
cada1f4dVdV = cd2.f.*cada1f3dVdV;
cada1f4dV = cd2.f*cada1f3dV;
cada1f4 = cd2.f*cada1f3;
cada1td1 = cada1f2dV;
cada1td1dV = cada1f4dVdV;
cada1td1 = cada1td1 + cada1f4dV;
CD.dVdV = cada1td1dV; CD.dV = cada1td1;
CD.f = cada1f2 + cada1f4;
%User Line: CD       = cd0+cd1*aoa+cd2*aoa.^2;
cada1f1dV = uminus(altitude.dV);
cada1f1 = uminus(altitude.f);
cada1f2dV = cada1f1dV/H;
cada1f2 = cada1f1/H;
cada2f1dV = exp(cada1f2).*cada1f2dV;
cada2f1 = exp(cada1f2);
cada1f3dVdV = cada1f2dV.*cada2f1dV;
cada1f3dV = cada2f1.*cada1f2dV;
cada1f3 = exp(cada1f2);
rho.dVdV = rho0.*cada1f3dVdV;
rho.dV = rho0*cada1f3dV;
rho.f = rho0*cada1f3;
%User Line: rho      = rho0*exp(-altitude/H);
cada1f1dV = cl1.f*aoa.dV;
cada1f1 = cl1.f*aoa.f;
CL.dV = cada1f1dV;
CL.f = cl0.f + cada1f1;
%User Line: CL       = cl0+cl1*aoa;
%User Line: % speedOfSound = spline(auxdata.Atmosphere(:,1),auxdata.Atmosphere(:,5),altitude);
%User Line: % mach = v./speedOfSound;
%User Line: % density = spline(auxdata.Atmosphere(:,1),auxdata.Atmosphere(:,4),altitude);
%User Line: % % interpolate coefficients
%User Line: % CD = auxdata.interp.Cd_spline(mach,rad2deg(aoa));
%User Line: % CL = auxdata.interp.Cl_spline(mach,rad2deg(aoa));
%User Line: %
%User Line: %
cada1f1dVdV = 0.5.*rho.dVdV;
cada1f1dV = 0.5*rho.dV;
cada1f1 = 0.5*rho.f;
cada2f1dV = 1.*v.f.^(1-1).*v.dV;
cada2f1 = v.f.^1;
cada2f2dV = 2.*cada2f1dV;
cada2f2 = 2*cada2f1;
cada1f2dVdV = v.dV.*cada2f2dV;
cada1f2dV = cada2f2.*v.dV;
cada1f2 = v.f.^2;
cada2f1 = size(cada1f1dV,1);
cada1td1 = zeros(cada2f1,2);
cada2td1 = zeros(size(cada1f2dV,1),2);
cada2td1(:,2) = cada1f1dV.*cada1f2dV;
cada2td1(:,1) = cada2td1(:,1) + cada1f2.*cada1f1dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1f2.*cada1f1dV;
cada1td1dV = cada2f1dV;
cada1td1(:,1) = cada2f1;
cada2f1 = cada1td1(:,2);
cada2td1 = zeros(size(cada1f1dV,1),2);
cada2td1(:,1) = cada1f2dV.*cada1f1dV;
cada2td1(:,2) = cada2td1(:,2) + cada1f1.*cada1f2dVdV;
cada2f2dV = cada2td1;
cada2f2 = cada1f1.*cada1f2dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),4);
cada2td1(:,Gator2Data.Index1) = cada2f3dV;
cada2td1(:,Gator2Data.Index2) = cada1td1dV(:,Gator2Data.Index3);
cada1td1dV = cada2td1;
cada1td1(:,2) = cada2f3;
q.dVdV = cada1td1dV; q.dV = cada1td1;
q.f = cada1f1.*cada1f2;
%User Line: q        = 0.5*rho.*v.^2;
cada1f1dVdV = S.*q.dVdV;
cada1f1dV = S*q.dV;
cada1f1 = q.f*S;
cada1tf1dV = CD.dV(:,Gator2Data.Index4);
cada1tf1 = CD.f(:,Gator1Data.Index1);
cada2f1 = size(cada1f1dV,1);
cada1td1 = zeros(cada2f1,3);
cada2td1 = zeros(size(cada1tf1dV,1),6);
cada2td1(:,Gator2Data.Index5) = cada1f1dV.*cada1tf1dV;
cada2tf1 = cada1tf1(:,Gator2Data.Index6);
cada2td1(:,Gator2Data.Index7) = cada2td1(:,Gator2Data.Index7) + cada2tf1.*cada1f1dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1tf1.*cada1f1dV;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index2) = cada2f1;
cada2f1 = cada1td1(:,3);
cada2tf1 = CD.dV(:,Gator2Data.Index8);
cada2td1 = zeros(size(cada1f1dV,1),3);
cada2td1(:,Gator2Data.Index9) = cada2tf1.*cada1f1dV;
cada2td1(:,3) = cada2td1(:,3) + cada1f1.*CD.dVdV;
cada2f2dV = cada2td1;
cada2f2 = cada1f1.*CD.dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),9);
cada2td1(:,Gator2Data.Index10) = cada2f3dV;
cada2td1(:,Gator2Data.Index11) = cada1td1dV(:,Gator2Data.Index12);
cada1td1dV = cada2td1;
cada1td1(:,3) = cada2f3;
cada1f2dVdV = cada1td1dV; cada1f2dV = cada1td1;
cada1f2 = cada1f1.*CD.f;
D.dVdV = cada1f2dVdV./mass;
D.dV = cada1f2dV/mass;
D.f = cada1f2/mass;
%User Line: D        = q.*S.*CD./mass;
cada1f1dVdV = S.*q.dVdV;
cada1f1dV = S*q.dV;
cada1f1 = q.f*S;
cada1tf1dV = CL.dV(:,Gator2Data.Index13);
cada1tf1 = CL.f(:,Gator1Data.Index3);
cada2f1 = size(cada1f1dV,1);
cada1td1 = zeros(cada2f1,3);
cada2td1 = zeros(size(cada1tf1dV,1),6);
cada2td1(:,Gator2Data.Index14) = cada1f1dV.*cada1tf1dV;
cada2tf1 = cada1tf1(:,Gator2Data.Index15);
cada2td1(:,Gator2Data.Index16) = cada2td1(:,Gator2Data.Index16) + cada2tf1.*cada1f1dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1tf1.*cada1f1dV;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index4) = cada2f1;
cada2f1 = cada1td1(:,3);
cada2tf1 = CL.dV(:,Gator2Data.Index17);
cada2f2dV = cada2tf1.*cada1f1dV;
cada2f2 = cada1f1.*CL.dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),8);
cada2td1(:,Gator2Data.Index18) = cada2f3dV;
cada2td1(:,Gator2Data.Index19) = cada1td1dV(:,Gator2Data.Index20);
cada1td1dV = cada2td1;
cada1td1(:,3) = cada2f3;
cada1f2dVdV = cada1td1dV; cada1f2dV = cada1td1;
cada1f2 = cada1f1.*CL.f;
L.dVdV = cada1f2dVdV./mass;
L.dV = cada1f2dV/mass;
L.f = cada1f2/mass;
%User Line: L        = q.*S.*CL./mass;
cada2f1dV = 1.*rad.f.^(1-1).*rad.dV;
cada2f1 = rad.f.^1;
cada2f2dV = 2.*cada2f1dV;
cada2f2 = 2*cada2f1;
cada1f1dVdV = rad.dV.*cada2f2dV;
cada1f1dV = cada2f2.*rad.dV;
cada1f1 = rad.f.^2;
cada2f1 = uminus(mu);
cada2f2dV = 2.*cada1f1.^(2-1).*cada1f1dV;
cada2f2 = cada1f1.^2;
cada2f3dV = -cada2f1./cada2f2.^2.*cada2f2dV;
cada2f3 = cada2f1./cada2f2;
cada2td1 = cada1f1dV.*cada2f3dV;
cada2td1 = cada2td1 + cada2f3.*cada1f1dVdV;
gravity.dVdV = cada2td1;
gravity.dV = cada2f3.*cada1f1dV;
gravity.f = mu./cada1f1;
%User Line: gravity  = mu./rad.^2;
%User Line: % ---------------------------------------------------%
%User Line: % ---- Evaluate Right-Hand Side of the Dynamics ---- %
%User Line: % ---------------------------------------------------%
cada2f1dV = -sin(fpa.f).*fpa.dV;
cada2f1 = cos(fpa.f);
cada1f1dVdV = fpa.dV.*cada2f1dV;
cada1f1dV = cada2f1.*fpa.dV;
cada1f1 = sin(fpa.f);
cada2f1 = size(v.dV,1);
cada1td1 = zeros(cada2f1,2);
cada2f1dV = v.dV.*cada1f1dV;
cada2f1 = cada1f1.*v.dV;
cada1td1dV = cada2f1dV;
cada1td1(:,1) = cada2f1;
cada2f1 = cada1td1(:,2);
cada2td1 = zeros(size(v.dV,1),2);
cada2td1(:,1) = cada1f1dV.*v.dV;
cada2td1(:,2) = cada2td1(:,2) + v.f.*cada1f1dVdV;
cada2f2dV = cada2td1;
cada2f2 = v.f.*cada1f1dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),3);
cada2td1(:,Gator2Data.Index21) = cada2f3dV;
cada2td1(:,2) = cada1td1dV(:,1);
cada1td1dV = cada2td1;
cada1td1(:,2) = cada2f3;
raddot.dVdV = cada1td1dV; raddot.dV = cada1td1;
raddot.f = v.f.*cada1f1;
%User Line: raddot = v.*sin(fpa);
cada2f1dV = cos(fpa.f).*fpa.dV;
cada2f1 = sin(fpa.f);
cada2f2dV = -cada2f1dV;
cada2f2 = uminus(cada2f1);
cada1f1dVdV = fpa.dV.*cada2f2dV;
cada1f1dV = cada2f2.*fpa.dV;
cada1f1 = cos(fpa.f);
cada2f1 = size(v.dV,1);
cada1td1 = zeros(cada2f1,2);
cada2f1dV = v.dV.*cada1f1dV;
cada2f1 = cada1f1.*v.dV;
cada1td1dV = cada2f1dV;
cada1td1(:,1) = cada2f1;
cada2f1 = cada1td1(:,2);
cada2td1 = zeros(size(v.dV,1),2);
cada2td1(:,1) = cada1f1dV.*v.dV;
cada2td1(:,2) = cada2td1(:,2) + v.f.*cada1f1dVdV;
cada2f2dV = cada2td1;
cada2f2 = v.f.*cada1f1dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),3);
cada2td1(:,Gator2Data.Index22) = cada2f3dV;
cada2td1(:,2) = cada1td1dV(:,1);
cada1td1dV = cada2td1;
cada1td1(:,2) = cada2f3;
cada1f2dVdV = cada1td1dV; cada1f2dV = cada1td1;
cada1f2 = v.f.*cada1f1;
cada2f1dV = -sin(azi.f).*azi.dV;
cada2f1 = cos(azi.f);
cada1f3dVdV = azi.dV.*cada2f1dV;
cada1f3dV = cada2f1.*azi.dV;
cada1f3 = sin(azi.f);
cada1tf1dV = cada1f3dV(:,Gator2Data.Index23);
cada1tf1 = cada1f3(:,Gator1Data.Index5);
cada2f1 = size(cada1f2dV,1);
cada1td1 = zeros(cada2f1,3);
cada2td1 = zeros(size(cada1tf1dV,1),5);
cada2td1(:,Gator2Data.Index24) = cada1f2dV.*cada1tf1dV;
cada2tf1 = cada1tf1(:,Gator2Data.Index25);
cada2td1(:,Gator2Data.Index26) = cada2td1(:,Gator2Data.Index26) + cada2tf1.*cada1f2dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1tf1.*cada1f2dV;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index6) = cada2f1;
cada2f1 = cada1td1(:,3);
cada2tf1 = cada1f3dV(:,Gator2Data.Index27);
cada2td1 = zeros(size(cada1f2dV,1),3);
cada2td1(:,Gator2Data.Index28) = cada2tf1.*cada1f2dV;
cada2td1(:,3) = cada2td1(:,3) + cada1f2.*cada1f3dVdV;
cada2f2dV = cada2td1;
cada2f2 = cada1f2.*cada1f3dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),8);
cada2td1(:,Gator2Data.Index29) = cada2f3dV;
cada2td1(:,Gator2Data.Index30) = cada1td1dV(:,Gator2Data.Index31);
cada1td1dV = cada2td1;
cada1td1(:,3) = cada2f3;
cada1f4dVdV = cada1td1dV; cada1f4dV = cada1td1;
cada1f4 = cada1f2.*cada1f3;
cada2f1dV = cos(lat.f).*lat.dV;
cada2f1 = sin(lat.f);
cada2f2dV = -cada2f1dV;
cada2f2 = uminus(cada2f1);
cada1f5dVdV = lat.dV.*cada2f2dV;
cada1f5dV = cada2f2.*lat.dV;
cada1f5 = cos(lat.f);
cada2f1 = size(rad.dV,1);
cada1td1 = zeros(cada2f1,2);
cada2f1dV = rad.dV.*cada1f5dV;
cada2f1 = cada1f5.*rad.dV;
cada1td1dV = cada2f1dV;
cada1td1(:,1) = cada2f1;
cada2f1 = cada1td1(:,2);
cada2td1 = zeros(size(rad.dV,1),2);
cada2td1(:,1) = cada1f5dV.*rad.dV;
cada2td1(:,2) = cada2td1(:,2) + rad.f.*cada1f5dVdV;
cada2f2dV = cada2td1;
cada2f2 = rad.f.*cada1f5dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),3);
cada2td1(:,Gator2Data.Index32) = cada2f3dV;
cada2td1(:,2) = cada1td1dV(:,1);
cada1td1dV = cada2td1;
cada1td1(:,2) = cada2f3;
cada1f6dVdV = cada1td1dV; cada1f6dV = cada1td1;
cada1f6 = rad.f.*cada1f5;
cada1tf1dV = cada1f6dV(:,Gator2Data.Index33);
cada1tf1 = cada1f6(:,Gator1Data.Index7);
cada2f1 = size(cada1f4dV,1);
cada1td1 = zeros(cada2f1,5);
cada2tf1 = cada1tf1(:,Gator2Data.Index34);
cada2td1 = zeros(size(cada1f4dVdV,1),14);
cada2td1(:,Gator2Data.Index35) = cada1f4dVdV./cada2tf1;
cada2tf1 = cada1f4dV(:,Gator2Data.Index36);
cada2tf2 = cada1tf1(:,Gator2Data.Index37);
cada2td1(:,Gator2Data.Index38) = cada2td1(:,Gator2Data.Index38) + -cada2tf1./cada2tf2.^2.*cada1tf1dV;
cada2f1dV = cada2td1;
cada2f1 = cada1f4dV./cada1tf1;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index8) = cada2f1;
cada1tf1dV = cada1f4dV(:,Gator2Data.Index39);
cada1tf1 = cada1f4(:,Gator1Data.Index9);
cada1tf2dV = cada1f6dV(:,Gator2Data.Index40);
cada1tf2 = cada1f6(:,Gator1Data.Index10);
cada2f1 = cada1td1(:,Gator1Data.Index11);
cada2f2dV = -cada1tf1dV;
cada2f2 = uminus(cada1tf1);
cada2tf2 = cada1tf2(:,Gator2Data.Index41);
cada2f3dV = 2.*cada2tf2.^(2-1).*cada1tf2dV;
cada2f3 = cada1tf2.^2;
cada2tf1 = cada2f3(:,Gator2Data.Index42);
cada2td1 = zeros(size(cada2f2dV,1),10);
cada2td1(:,Gator2Data.Index43) = cada2f2dV./cada2tf1;
cada2tf1 = cada2f2(:,Gator2Data.Index44);
cada2tf2 = cada2f3(:,Gator2Data.Index45);
cada2td1(:,Gator2Data.Index46) = cada2td1(:,Gator2Data.Index46) + -cada2tf1./cada2tf2.^2.*cada2f3dV;
cada2f4dV = cada2td1;
cada2f4 = cada2f2./cada2f3;
cada2tf1 = cada1f6dV(:,Gator2Data.Index47);
cada2td1 = cada2tf1.*cada2f4dV;
cada2tf1 = cada2f4(:,Gator2Data.Index48);
cada2td1(:,Gator2Data.Index49) = cada2td1(:,Gator2Data.Index49) + cada2tf1.*cada1f6dVdV;
cada2f5dV = cada2td1;
cada2f5 = cada2f4.*cada1f6dV;
cada2f6dV = cada2f5dV;
cada2f6 = cada2f1 + cada2f5;
cada2td1 = zeros(size(cada1td1,1),24);
cada2td1(:,Gator2Data.Index50) = cada2f6dV;
cada2td1(:,Gator2Data.Index51) = cada1td1dV(:,Gator2Data.Index52);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index11) = cada2f6;
londot.dVdV = cada1td1dV; londot.dV = cada1td1;
londot.f = cada1f4./cada1f6;
%User Line: londot = v.*cos(fpa).*sin(azi)./(rad.*cos(lat));
cada2f1dV = cos(fpa.f).*fpa.dV;
cada2f1 = sin(fpa.f);
cada2f2dV = -cada2f1dV;
cada2f2 = uminus(cada2f1);
cada1f1dVdV = fpa.dV.*cada2f2dV;
cada1f1dV = cada2f2.*fpa.dV;
cada1f1 = cos(fpa.f);
cada2f1 = size(v.dV,1);
cada1td1 = zeros(cada2f1,2);
cada2f1dV = v.dV.*cada1f1dV;
cada2f1 = cada1f1.*v.dV;
cada1td1dV = cada2f1dV;
cada1td1(:,1) = cada2f1;
cada2f1 = cada1td1(:,2);
cada2td1 = zeros(size(v.dV,1),2);
cada2td1(:,1) = cada1f1dV.*v.dV;
cada2td1(:,2) = cada2td1(:,2) + v.f.*cada1f1dVdV;
cada2f2dV = cada2td1;
cada2f2 = v.f.*cada1f1dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),3);
cada2td1(:,Gator2Data.Index53) = cada2f3dV;
cada2td1(:,2) = cada1td1dV(:,1);
cada1td1dV = cada2td1;
cada1td1(:,2) = cada2f3;
cada1f2dVdV = cada1td1dV; cada1f2dV = cada1td1;
cada1f2 = v.f.*cada1f1;
cada2f1dV = cos(azi.f).*azi.dV;
cada2f1 = sin(azi.f);
cada2f2dV = -cada2f1dV;
cada2f2 = uminus(cada2f1);
cada1f3dVdV = azi.dV.*cada2f2dV;
cada1f3dV = cada2f2.*azi.dV;
cada1f3 = cos(azi.f);
cada1tf1dV = cada1f3dV(:,Gator2Data.Index54);
cada1tf1 = cada1f3(:,Gator1Data.Index12);
cada2f1 = size(cada1f2dV,1);
cada1td1 = zeros(cada2f1,3);
cada2td1 = zeros(size(cada1tf1dV,1),5);
cada2td1(:,Gator2Data.Index55) = cada1f2dV.*cada1tf1dV;
cada2tf1 = cada1tf1(:,Gator2Data.Index56);
cada2td1(:,Gator2Data.Index57) = cada2td1(:,Gator2Data.Index57) + cada2tf1.*cada1f2dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1tf1.*cada1f2dV;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index13) = cada2f1;
cada2f1 = cada1td1(:,3);
cada2tf1 = cada1f3dV(:,Gator2Data.Index58);
cada2td1 = zeros(size(cada1f2dV,1),3);
cada2td1(:,Gator2Data.Index59) = cada2tf1.*cada1f2dV;
cada2td1(:,3) = cada2td1(:,3) + cada1f2.*cada1f3dVdV;
cada2f2dV = cada2td1;
cada2f2 = cada1f2.*cada1f3dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),8);
cada2td1(:,Gator2Data.Index60) = cada2f3dV;
cada2td1(:,Gator2Data.Index61) = cada1td1dV(:,Gator2Data.Index62);
cada1td1dV = cada2td1;
cada1td1(:,3) = cada2f3;
cada1f4dVdV = cada1td1dV; cada1f4dV = cada1td1;
cada1f4 = cada1f2.*cada1f3;
cada1tf1dV = rad.dV(:,Gator2Data.Index63);
cada1tf1 = rad.f(:,Gator1Data.Index14);
cada2f1 = size(cada1f4dV,1);
cada1td1 = zeros(cada2f1,4);
cada2tf1 = cada1tf1(:,Gator2Data.Index64);
cada2td1 = zeros(size(cada1f4dVdV,1),11);
cada2td1(:,Gator2Data.Index65) = cada1f4dVdV./cada2tf1;
cada2td1(:,Gator2Data.Index66) = cada2td1(:,Gator2Data.Index66) + -cada1f4dV./cada1tf1.^2.*cada1tf1dV;
cada2f1dV = cada2td1;
cada2f1 = cada1f4dV./cada1tf1;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index15) = cada2f1;
cada2f1 = cada1td1(:,1);
cada2f2dV = -cada1f4dV;
cada2f2 = uminus(cada1f4);
cada2f3dV = 2.*rad.f.^(2-1).*rad.dV;
cada2f3 = rad.f.^2;
cada2tf1 = cada2f3(:,Gator2Data.Index67);
cada2td1 = zeros(size(cada2f2dV,1),4);
cada2td1(:,Gator2Data.Index68) = cada2f2dV./cada2tf1;
cada2td1(:,1) = cada2td1(:,1) + -cada2f2./cada2f3.^2.*cada2f3dV;
cada2f4dV = cada2td1;
cada2f4 = cada2f2./cada2f3;
cada2tf1 = rad.dV(:,Gator2Data.Index69);
cada2f5dV = cada2tf1.*cada2f4dV;
cada2f5 = cada2f4.*rad.dV;
cada2f6dV = cada2f5dV;
cada2f6 = cada2f1 + cada2f5;
cada2td1 = zeros(size(cada1td1,1),15);
cada2td1(:,Gator2Data.Index70) = cada2f6dV;
cada2td1(:,Gator2Data.Index71) = cada1td1dV(:,Gator2Data.Index72);
cada1td1dV = cada2td1;
cada1td1(:,1) = cada2f6;
latdot.dVdV = cada1td1dV; latdot.dV = cada1td1;
latdot.f = cada1f4./rad.f;
%User Line: latdot = v.*cos(fpa).*cos(azi)./rad;
cada1f1dVdV = -D.dVdV;
cada1f1dV = uminus(D.dV);
cada1f1 = uminus(D.f);
cada2f1dV = -sin(fpa.f).*fpa.dV;
cada2f1 = cos(fpa.f);
cada1f2dVdV = fpa.dV.*cada2f1dV;
cada1f2dV = cada2f1.*fpa.dV;
cada1f2 = sin(fpa.f);
cada2f1 = size(gravity.dV,1);
cada1td1 = zeros(cada2f1,2);
cada2td1 = zeros(size(cada1f2dV,1),2);
cada2td1(:,2) = gravity.dV.*cada1f2dV;
cada2td1(:,1) = cada2td1(:,1) + cada1f2.*gravity.dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1f2.*gravity.dV;
cada1td1dV = cada2f1dV;
cada1td1(:,1) = cada2f1;
cada2f1 = cada1td1(:,2);
cada2td1 = zeros(size(gravity.dV,1),2);
cada2td1(:,1) = cada1f2dV.*gravity.dV;
cada2td1(:,2) = cada2td1(:,2) + gravity.f.*cada1f2dVdV;
cada2f2dV = cada2td1;
cada2f2 = gravity.f.*cada1f2dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),4);
cada2td1(:,Gator2Data.Index73) = cada2f3dV;
cada2td1(:,Gator2Data.Index74) = cada1td1dV(:,Gator2Data.Index75);
cada1td1dV = cada2td1;
cada1td1(:,2) = cada2f3;
cada1f3dVdV = cada1td1dV; cada1f3dV = cada1td1;
cada1f3 = gravity.f.*cada1f2;
cada2f1 = size(cada1f1dV,1);
cada1td1 = zeros(cada2f1,4);
cada1td1dV = cada1f1dVdV;
cada1td1(:,Gator1Data.Index16) = cada1f1dV;
cada2f1dV = cada1td1dV(:,Gator2Data.Index76);
cada2f1 = cada1td1(:,Gator1Data.Index17);
cada2f2dV = -cada1f3dVdV;
cada2f2 = uminus(cada1f3dV);
cada2td1 = zeros(size(cada2f1dV,1),6);
cada2td1(:,Gator2Data.Index77) = cada2f1dV;
cada2td1(:,Gator2Data.Index78) = cada2td1(:,Gator2Data.Index78) + cada2f2dV;
cada2f3dV = cada2td1;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),12);
cada2td1(:,Gator2Data.Index79) = cada2f3dV;
cada2td1(:,Gator2Data.Index80) = cada1td1dV(:,Gator2Data.Index81);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index17) = cada2f3;
vdot.dVdV = cada1td1dV; vdot.dV = cada1td1;
vdot.f = cada1f1 - cada1f3;
%User Line: vdot   = -D-gravity.*sin(fpa);
cada2f1dV = cos(bank.f).*bank.dV;
cada2f1 = sin(bank.f);
cada2f2dV = -cada2f1dV;
cada2f2 = uminus(cada2f1);
cada1f1dVdV = bank.dV.*cada2f2dV;
cada1f1dV = cada2f2.*bank.dV;
cada1f1 = cos(bank.f);
cada1tf1dV = cada1f1dV(:,Gator2Data.Index82);
cada1tf1 = cada1f1(:,Gator1Data.Index18);
cada2f1 = size(L.dV,1);
cada1td1 = zeros(cada2f1,4);
cada2td1 = zeros(size(cada1tf1dV,1),11);
cada2td1(:,Gator2Data.Index83) = L.dV.*cada1tf1dV;
cada2tf1 = cada1tf1(:,Gator2Data.Index84);
cada2td1(:,Gator2Data.Index85) = cada2td1(:,Gator2Data.Index85) + cada2tf1.*L.dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1tf1.*L.dV;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index19) = cada2f1;
cada2f1 = cada1td1(:,4);
cada2tf1 = cada1f1dV(:,Gator2Data.Index86);
cada2td1 = zeros(size(L.dV,1),4);
cada2td1(:,Gator2Data.Index87) = cada2tf1.*L.dV;
cada2td1(:,4) = cada2td1(:,4) + L.f.*cada1f1dVdV;
cada2f2dV = cada2td1;
cada2f2 = L.f.*cada1f1dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),15);
cada2td1(:,Gator2Data.Index88) = cada2f3dV;
cada2td1(:,Gator2Data.Index89) = cada1td1dV(:,Gator2Data.Index90);
cada1td1dV = cada2td1;
cada1td1(:,4) = cada2f3;
cada1f2dVdV = cada1td1dV; cada1f2dV = cada1td1;
cada1f2 = L.f.*cada1f1;
cada2f1dV = cos(fpa.f).*fpa.dV;
cada2f1 = sin(fpa.f);
cada2f2dV = -cada2f1dV;
cada2f2 = uminus(cada2f1);
cada1f3dVdV = fpa.dV.*cada2f2dV;
cada1f3dV = cada2f2.*fpa.dV;
cada1f3 = cos(fpa.f);
cada2f1dV = 1.*v.f.^(1-1).*v.dV;
cada2f1 = v.f.^1;
cada2f2dV = 2.*cada2f1dV;
cada2f2 = 2*cada2f1;
cada1f4dVdV = v.dV.*cada2f2dV;
cada1f4dV = cada2f2.*v.dV;
cada1f4 = v.f.^2;
cada2f1 = size(cada1f4dV,1);
cada1td1 = zeros(cada2f1,2);
cada2td1 = zeros(size(cada1f4dVdV,1),2);
cada2td1(:,2) = cada1f4dVdV./rad.f;
cada2td1(:,1) = cada2td1(:,1) + -cada1f4dV./rad.f.^2.*rad.dV;
cada2f1dV = cada2td1;
cada2f1 = cada1f4dV./rad.f;
cada1td1dV = cada2f1dV;
cada1td1(:,2) = cada2f1;
cada2f1 = cada1td1(:,1);
cada2f2dV = -cada1f4dV;
cada2f2 = uminus(cada1f4);
cada2f3dV = 2.*rad.f.^(2-1).*rad.dV;
cada2f3 = rad.f.^2;
cada2td1 = zeros(size(cada2f2dV,1),2);
cada2td1(:,2) = cada2f2dV./cada2f3;
cada2td1(:,1) = cada2td1(:,1) + -cada2f2./cada2f3.^2.*cada2f3dV;
cada2f4dV = cada2td1;
cada2f4 = cada2f2./cada2f3;
cada2tf1 = rad.dV(:,Gator2Data.Index91);
cada2f5dV = cada2tf1.*cada2f4dV;
cada2f5 = cada2f4.*rad.dV;
cada2f6dV = cada2f5dV;
cada2f6 = cada2f1 + cada2f5;
cada2td1 = zeros(size(cada1td1,1),4);
cada2td1(:,Gator2Data.Index92) = cada2f6dV;
cada2td1(:,Gator2Data.Index93) = cada1td1dV(:,Gator2Data.Index94);
cada1td1dV = cada2td1;
cada1td1(:,1) = cada2f6;
cada1f5dVdV = cada1td1dV; cada1f5dV = cada1td1;
cada1f5 = cada1f4./rad.f;
cada2f1 = size(gravity.dV,1);
cada1td1 = zeros(cada2f1,2);
cada1td1dV = gravity.dVdV;
cada1td1(:,1) = gravity.dV;
cada2f1dV = -cada1f5dVdV;
cada2f1 = uminus(cada1f5dV);
cada2td1 = zeros(size(cada1td1dV,1),4);
cada2td1(:,1) = cada1td1dV;
cada2td1 = cada2td1 + cada2f1dV;
cada1td1dV = cada2td1;
cada1td1 = cada1td1 + cada2f1;
cada1f6dVdV = cada1td1dV; cada1f6dV = cada1td1;
cada1f6 = gravity.f - cada1f5;
cada2f1 = size(cada1f3dV,1);
cada1td1 = zeros(cada2f1,3);
cada2tf1 = cada1f3dV(:,Gator2Data.Index95);
cada2td1 = zeros(size(cada1f6dV,1),3);
cada2td1(:,Gator2Data.Index96) = cada2tf1.*cada1f6dV;
cada2td1(:,3) = cada2td1(:,3) + cada1f6.*cada1f3dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1f6.*cada1f3dV;
cada1td1dV = cada2f1dV;
cada1td1(:,3) = cada2f1;
cada1tf1dV = cada1f3dV(:,Gator2Data.Index97);
cada1tf1 = cada1f3(:,Gator1Data.Index20);
cada2f1 = cada1td1(:,Gator1Data.Index21);
cada2td1 = zeros(size(cada1tf1dV,1),6);
cada2td1(:,Gator2Data.Index98) = cada1f6dV.*cada1tf1dV;
cada2tf1 = cada1tf1(:,Gator2Data.Index99);
cada2td1(:,Gator2Data.Index100) = cada2td1(:,Gator2Data.Index100) + cada2tf1.*cada1f6dVdV;
cada2f2dV = cada2td1;
cada2f2 = cada1tf1.*cada1f6dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),9);
cada2td1(:,Gator2Data.Index101) = cada2f3dV;
cada2td1(:,Gator2Data.Index102) = cada1td1dV(:,Gator2Data.Index103);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index21) = cada2f3;
cada1f7dVdV = cada1td1dV; cada1f7dV = cada1td1;
cada1f7 = cada1f3.*cada1f6;
cada2f1 = size(cada1f2dV,1);
cada1td1 = zeros(cada2f1,5);
cada1td1dV = cada1f2dVdV;
cada1td1(:,Gator1Data.Index22) = cada1f2dV;
cada2f1dV = cada1td1dV(:,Gator2Data.Index104);
cada2f1 = cada1td1(:,Gator1Data.Index23);
cada2f2dV = -cada1f7dVdV;
cada2f2 = uminus(cada1f7dV);
cada2td1 = zeros(size(cada2f1dV,1),13);
cada2td1(:,Gator2Data.Index105) = cada2f1dV;
cada2td1(:,Gator2Data.Index106) = cada2td1(:,Gator2Data.Index106) + cada2f2dV;
cada2f3dV = cada2td1;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),20);
cada2td1(:,Gator2Data.Index107) = cada2f3dV;
cada2td1(:,Gator2Data.Index108) = cada1td1dV(:,Gator2Data.Index109);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index23) = cada2f3;
cada1f8dVdV = cada1td1dV; cada1f8dV = cada1td1;
cada1f8 = cada1f2 - cada1f7;
cada1tf1dV = v.dV(:,Gator2Data.Index110);
cada1tf1 = v.f(:,Gator1Data.Index24);
cada2tf1 = cada1tf1(:,Gator2Data.Index111);
cada2td1 = cada1f8dVdV./cada2tf1;
cada2td1(:,Gator2Data.Index112) = cada2td1(:,Gator2Data.Index112) + -cada1f8dV./cada1tf1.^2.*cada1tf1dV;
cada1td1dV = cada2td1;
cada1td1 = cada1f8dV./cada1tf1;
cada2f1dV = cada1td1dV(:,Gator2Data.Index113);
cada2f1 = cada1td1(:,2);
cada2f2dV = -cada1f8dV;
cada2f2 = uminus(cada1f8);
cada2f3dV = 2.*v.f.^(2-1).*v.dV;
cada2f3 = v.f.^2;
cada2tf1 = cada2f3(:,Gator2Data.Index114);
cada2td1 = cada2f2dV./cada2tf1;
cada2td1(:,2) = cada2td1(:,2) + -cada2f2./cada2f3.^2.*cada2f3dV;
cada2f4dV = cada2td1;
cada2f4 = cada2f2./cada2f3;
cada2tf1 = v.dV(:,Gator2Data.Index115);
cada2f5dV = cada2tf1.*cada2f4dV;
cada2f5 = cada2f4.*v.dV;
cada2td1 = cada2f1dV;
cada2td1 = cada2td1 + cada2f5dV;
cada2f6dV = cada2td1;
cada2f6 = cada2f1 + cada2f5;
cada2td1 = zeros(size(cada1td1,1),20);
cada2td1(:,Gator2Data.Index116) = cada2f6dV;
cada2td1(:,Gator2Data.Index117) = cada1td1dV(:,Gator2Data.Index118);
cada1td1dV = cada2td1;
cada1td1(:,2) = cada2f6;
fpadot.dVdV = cada1td1dV; fpadot.dV = cada1td1;
fpadot.f = cada1f8./v.f;
%User Line: fpadot = (L.*cos(bank)-cos(fpa).*(gravity-v.^2./rad))./v;
cada2f1dV = -sin(bank.f).*bank.dV;
cada2f1 = cos(bank.f);
cada1f1dVdV = bank.dV.*cada2f1dV;
cada1f1dV = cada2f1.*bank.dV;
cada1f1 = sin(bank.f);
cada1tf1dV = cada1f1dV(:,Gator2Data.Index119);
cada1tf1 = cada1f1(:,Gator1Data.Index25);
cada2f1 = size(L.dV,1);
cada1td1 = zeros(cada2f1,4);
cada2td1 = zeros(size(cada1tf1dV,1),11);
cada2td1(:,Gator2Data.Index120) = L.dV.*cada1tf1dV;
cada2tf1 = cada1tf1(:,Gator2Data.Index121);
cada2td1(:,Gator2Data.Index122) = cada2td1(:,Gator2Data.Index122) + cada2tf1.*L.dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1tf1.*L.dV;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index26) = cada2f1;
cada2f1 = cada1td1(:,4);
cada2tf1 = cada1f1dV(:,Gator2Data.Index123);
cada2td1 = zeros(size(L.dV,1),4);
cada2td1(:,Gator2Data.Index124) = cada2tf1.*L.dV;
cada2td1(:,4) = cada2td1(:,4) + L.f.*cada1f1dVdV;
cada2f2dV = cada2td1;
cada2f2 = L.f.*cada1f1dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),15);
cada2td1(:,Gator2Data.Index125) = cada2f3dV;
cada2td1(:,Gator2Data.Index126) = cada1td1dV(:,Gator2Data.Index127);
cada1td1dV = cada2td1;
cada1td1(:,4) = cada2f3;
cada1f2dVdV = cada1td1dV; cada1f2dV = cada1td1;
cada1f2 = L.f.*cada1f1;
cada2f1dV = cos(fpa.f).*fpa.dV;
cada2f1 = sin(fpa.f);
cada2f2dV = -cada2f1dV;
cada2f2 = uminus(cada2f1);
cada1f3dVdV = fpa.dV.*cada2f2dV;
cada1f3dV = cada2f2.*fpa.dV;
cada1f3 = cos(fpa.f);
cada1tf1dV = cada1f3dV(:,Gator2Data.Index128);
cada1tf1 = cada1f3(:,Gator1Data.Index27);
cada2f1 = size(cada1f2dV,1);
cada1td1 = zeros(cada2f1,5);
cada2tf1 = cada1tf1(:,Gator2Data.Index129);
cada2td1 = zeros(size(cada1f2dVdV,1),19);
cada2td1(:,Gator2Data.Index130) = cada1f2dVdV./cada2tf1;
cada2td1(:,Gator2Data.Index131) = cada2td1(:,Gator2Data.Index131) + -cada1f2dV./cada1tf1.^2.*cada1tf1dV;
cada2f1dV = cada2td1;
cada2f1 = cada1f2dV./cada1tf1;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index28) = cada2f1;
cada2f1 = cada1td1(:,3);
cada2f2dV = -cada1f2dV;
cada2f2 = uminus(cada1f2);
cada2f3dV = 2.*cada1f3.^(2-1).*cada1f3dV;
cada2f3 = cada1f3.^2;
cada2tf1 = cada2f3(:,Gator2Data.Index132);
cada2td1 = zeros(size(cada2f2dV,1),5);
cada2td1(:,Gator2Data.Index133) = cada2f2dV./cada2tf1;
cada2td1(:,3) = cada2td1(:,3) + -cada2f2./cada2f3.^2.*cada2f3dV;
cada2f4dV = cada2td1;
cada2f4 = cada2f2./cada2f3;
cada2tf1 = cada1f3dV(:,Gator2Data.Index134);
cada2td1 = cada2tf1.*cada2f4dV;
cada2td1(:,3) = cada2td1(:,3) + cada2f4.*cada1f3dVdV;
cada2f5dV = cada2td1;
cada2f5 = cada2f4.*cada1f3dV;
cada2f6dV = cada2f5dV;
cada2f6 = cada2f1 + cada2f5;
cada2td1 = zeros(size(cada1td1,1),24);
cada2td1(:,Gator2Data.Index135) = cada2f6dV;
cada2td1(:,Gator2Data.Index136) = cada1td1dV(:,Gator2Data.Index137);
cada1td1dV = cada2td1;
cada1td1(:,3) = cada2f6;
cada1f4dVdV = cada1td1dV; cada1f4dV = cada1td1;
cada1f4 = cada1f2./cada1f3;
cada2f1dV = 1.*v.f.^(1-1).*v.dV;
cada2f1 = v.f.^1;
cada2f2dV = 2.*cada2f1dV;
cada2f2 = 2*cada2f1;
cada1f5dVdV = v.dV.*cada2f2dV;
cada1f5dV = cada2f2.*v.dV;
cada1f5 = v.f.^2;
cada2f1dV = cos(fpa.f).*fpa.dV;
cada2f1 = sin(fpa.f);
cada2f2dV = -cada2f1dV;
cada2f2 = uminus(cada2f1);
cada1f6dVdV = fpa.dV.*cada2f2dV;
cada1f6dV = cada2f2.*fpa.dV;
cada1f6 = cos(fpa.f);
cada2f1 = size(cada1f5dV,1);
cada1td1 = zeros(cada2f1,2);
cada2td1 = zeros(size(cada1f6dV,1),2);
cada2td1(:,2) = cada1f5dV.*cada1f6dV;
cada2td1(:,1) = cada2td1(:,1) + cada1f6.*cada1f5dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1f6.*cada1f5dV;
cada1td1dV = cada2f1dV;
cada1td1(:,1) = cada2f1;
cada2f1 = cada1td1(:,2);
cada2td1 = zeros(size(cada1f5dV,1),2);
cada2td1(:,1) = cada1f6dV.*cada1f5dV;
cada2td1(:,2) = cada2td1(:,2) + cada1f5.*cada1f6dVdV;
cada2f2dV = cada2td1;
cada2f2 = cada1f5.*cada1f6dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),4);
cada2td1(:,Gator2Data.Index138) = cada2f3dV;
cada2td1(:,Gator2Data.Index139) = cada1td1dV(:,Gator2Data.Index140);
cada1td1dV = cada2td1;
cada1td1(:,2) = cada2f3;
cada1f7dVdV = cada1td1dV; cada1f7dV = cada1td1;
cada1f7 = cada1f5.*cada1f6;
cada2f1dV = -sin(azi.f).*azi.dV;
cada2f1 = cos(azi.f);
cada1f8dVdV = azi.dV.*cada2f1dV;
cada1f8dV = cada2f1.*azi.dV;
cada1f8 = sin(azi.f);
cada1tf1dV = cada1f8dV(:,Gator2Data.Index141);
cada1tf1 = cada1f8(:,Gator1Data.Index29);
cada2f1 = size(cada1f7dV,1);
cada1td1 = zeros(cada2f1,3);
cada2td1 = zeros(size(cada1tf1dV,1),6);
cada2td1(:,Gator2Data.Index142) = cada1f7dV.*cada1tf1dV;
cada2tf1 = cada1tf1(:,Gator2Data.Index143);
cada2td1(:,Gator2Data.Index144) = cada2td1(:,Gator2Data.Index144) + cada2tf1.*cada1f7dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1tf1.*cada1f7dV;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index30) = cada2f1;
cada2f1 = cada1td1(:,3);
cada2tf1 = cada1f8dV(:,Gator2Data.Index145);
cada2td1 = zeros(size(cada1f7dV,1),3);
cada2td1(:,Gator2Data.Index146) = cada2tf1.*cada1f7dV;
cada2td1(:,3) = cada2td1(:,3) + cada1f7.*cada1f8dVdV;
cada2f2dV = cada2td1;
cada2f2 = cada1f7.*cada1f8dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),9);
cada2td1(:,Gator2Data.Index147) = cada2f3dV;
cada2td1(:,Gator2Data.Index148) = cada1td1dV(:,Gator2Data.Index149);
cada1td1dV = cada2td1;
cada1td1(:,3) = cada2f3;
cada1f9dVdV = cada1td1dV; cada1f9dV = cada1td1;
cada1f9 = cada1f7.*cada1f8;
cada2f1dV = sec(lat.f).*tan(lat.f).*lat.dV;
cada2f1 = sec(lat.f);
cada2f2dV = 2.*cada2f1.^(2-1).*cada2f1dV;
cada2f2 = cada2f1.^2;
cada1f10dVdV = lat.dV.*cada2f2dV;
cada1f10dV = cada2f2.*lat.dV;
cada1f10 = tan(lat.f);
cada1tf1dV = cada1f10dV(:,Gator2Data.Index150);
cada1tf1 = cada1f10(:,Gator1Data.Index31);
cada2f1 = size(cada1f9dV,1);
cada1td1 = zeros(cada2f1,4);
cada2td1 = zeros(size(cada1tf1dV,1),12);
cada2td1(:,Gator2Data.Index151) = cada1f9dV.*cada1tf1dV;
cada2tf1 = cada1tf1(:,Gator2Data.Index152);
cada2td1(:,Gator2Data.Index153) = cada2td1(:,Gator2Data.Index153) + cada2tf1.*cada1f9dVdV;
cada2f1dV = cada2td1;
cada2f1 = cada1tf1.*cada1f9dV;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index32) = cada2f1;
cada2f1 = cada1td1(:,1);
cada2tf1 = cada1f10dV(:,Gator2Data.Index154);
cada2td1 = zeros(size(cada1f9dV,1),4);
cada2td1(:,Gator2Data.Index155) = cada2tf1.*cada1f9dV;
cada2td1(:,1) = cada2td1(:,1) + cada1f9.*cada1f10dVdV;
cada2f2dV = cada2td1;
cada2f2 = cada1f9.*cada1f10dV;
cada2f3dV = cada2f2dV;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(size(cada1td1,1),16);
cada2td1(:,Gator2Data.Index156) = cada2f3dV;
cada2td1(:,Gator2Data.Index157) = cada1td1dV(:,Gator2Data.Index158);
cada1td1dV = cada2td1;
cada1td1(:,1) = cada2f3;
cada1f11dVdV = cada1td1dV; cada1f11dV = cada1td1;
cada1f11 = cada1f9.*cada1f10;
cada1tf1dV = rad.dV(:,Gator2Data.Index159);
cada1tf1 = rad.f(:,Gator1Data.Index33);
cada2f1 = size(cada1f11dV,1);
cada1td1 = zeros(cada2f1,5);
cada2tf1 = cada1tf1(:,Gator2Data.Index160);
cada2td1 = zeros(size(cada1f11dVdV,1),20);
cada2td1(:,Gator2Data.Index161) = cada1f11dVdV./cada2tf1;
cada2td1(:,Gator2Data.Index162) = cada2td1(:,Gator2Data.Index162) + -cada1f11dV./cada1tf1.^2.*cada1tf1dV;
cada2f1dV = cada2td1;
cada2f1 = cada1f11dV./cada1tf1;
cada1td1dV = cada2f1dV;
cada1td1(:,Gator1Data.Index34) = cada2f1;
cada2f1 = cada1td1(:,1);
cada2f2dV = -cada1f11dV;
cada2f2 = uminus(cada1f11);
cada2f3dV = 2.*rad.f.^(2-1).*rad.dV;
cada2f3 = rad.f.^2;
cada2tf1 = cada2f3(:,Gator2Data.Index163);
cada2td1 = zeros(size(cada2f2dV,1),5);
cada2td1(:,Gator2Data.Index164) = cada2f2dV./cada2tf1;
cada2td1(:,1) = cada2td1(:,1) + -cada2f2./cada2f3.^2.*cada2f3dV;
cada2f4dV = cada2td1;
cada2f4 = cada2f2./cada2f3;
cada2tf1 = rad.dV(:,Gator2Data.Index165);
cada2f5dV = cada2tf1.*cada2f4dV;
cada2f5 = cada2f4.*rad.dV;
cada2f6dV = cada2f5dV;
cada2f6 = cada2f1 + cada2f5;
cada2td1 = zeros(size(cada1td1,1),25);
cada2td1(:,Gator2Data.Index166) = cada2f6dV;
cada2td1(:,Gator2Data.Index167) = cada1td1dV(:,Gator2Data.Index168);
cada1td1dV = cada2td1;
cada1td1(:,1) = cada2f6;
cada1f12dVdV = cada1td1dV; cada1f12dV = cada1td1;
cada1f12 = cada1f11./rad.f;
cada2f1 = size(cada1f4dV,1);
cada1td1 = zeros(cada2f1,7);
cada1td1dV = cada1f4dVdV;
cada1td1(:,Gator1Data.Index35) = cada1f4dV;
cada2f1dV = cada1td1dV(:,Gator2Data.Index169);
cada2f1 = cada1td1(:,Gator1Data.Index36);
cada2td1 = zeros(size(cada2f1dV,1),31);
cada2td1(:,Gator2Data.Index170) = cada2f1dV;
cada2td1(:,Gator2Data.Index171) = cada2td1(:,Gator2Data.Index171) + cada1f12dVdV;
cada2f2dV = cada2td1;
cada2f2 = cada2f1 + cada1f12dV;
cada2td1 = zeros(size(cada1td1,1),40);
cada2td1(:,Gator2Data.Index172) = cada2f2dV;
cada2td1(:,Gator2Data.Index173) = cada1td1dV(:,Gator2Data.Index174);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index36) = cada2f2;
cada1f13dVdV = cada1td1dV; cada1f13dV = cada1td1;
cada1f13 = cada1f4 + cada1f12;
cada1tf1dV = v.dV(:,Gator2Data.Index175);
cada1tf1 = v.f(:,Gator1Data.Index37);
cada2tf1 = cada1tf1(:,Gator2Data.Index176);
cada2td1 = cada1f13dVdV./cada2tf1;
cada2td1(:,Gator2Data.Index177) = cada2td1(:,Gator2Data.Index177) + -cada1f13dV./cada1tf1.^2.*cada1tf1dV;
cada1td1dV = cada2td1;
cada1td1 = cada1f13dV./cada1tf1;
cada2f1dV = cada1td1dV(:,Gator2Data.Index178);
cada2f1 = cada1td1(:,3);
cada2f2dV = -cada1f13dV;
cada2f2 = uminus(cada1f13);
cada2f3dV = 2.*v.f.^(2-1).*v.dV;
cada2f3 = v.f.^2;
cada2tf1 = cada2f3(:,Gator2Data.Index179);
cada2td1 = cada2f2dV./cada2tf1;
cada2td1(:,3) = cada2td1(:,3) + -cada2f2./cada2f3.^2.*cada2f3dV;
cada2f4dV = cada2td1;
cada2f4 = cada2f2./cada2f3;
cada2tf1 = v.dV(:,Gator2Data.Index180);
cada2f5dV = cada2tf1.*cada2f4dV;
cada2f5 = cada2f4.*v.dV;
cada2td1 = cada2f1dV;
cada2td1 = cada2td1 + cada2f5dV;
cada2f6dV = cada2td1;
cada2f6 = cada2f1 + cada2f5;
cada2td1 = zeros(size(cada1td1,1),40);
cada2td1(:,Gator2Data.Index181) = cada2f6dV;
cada2td1(:,Gator2Data.Index182) = cada1td1dV(:,Gator2Data.Index183);
cada1td1dV = cada2td1;
cada1td1(:,3) = cada2f6;
azidot.dVdV = cada1td1dV; azidot.dV = cada1td1;
azidot.f = cada1f13./v.f;
%User Line: azidot = (L.*sin(bank)./cos(fpa)+v.^2.*cos(fpa).*sin(azi).*tan(lat)./rad)./v;
cada2f1 = size(raddot.f,1);
cada1td1 = zeros(cada2f1,27);
cada1td1dV = raddot.dVdV;
cada1td1(:,Gator1Data.Index38) = raddot.dV;
cada2td1 = zeros(size(cada1td1,1),27);
cada2td1(:,Gator2Data.Index184) = londot.dVdV;
cada2td1(:,Gator2Data.Index185) = cada1td1dV(:,Gator2Data.Index186);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index39) = londot.dV;
cada2td1 = zeros(size(cada1td1,1),42);
cada2td1(:,Gator2Data.Index187) = latdot.dVdV;
cada2td1(:,Gator2Data.Index188) = cada1td1dV(:,Gator2Data.Index189);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index40) = latdot.dV;
cada2td1 = zeros(size(cada1td1,1),54);
cada2td1(:,Gator2Data.Index190) = vdot.dVdV;
cada2td1(:,Gator2Data.Index191) = cada1td1dV(:,Gator2Data.Index192);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index41) = vdot.dV;
cada2td1 = zeros(size(cada1td1,1),74);
cada2td1(:,Gator2Data.Index193) = fpadot.dVdV;
cada2td1(:,Gator2Data.Index194) = cada1td1dV(:,Gator2Data.Index195);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index42) = fpadot.dV;
cada2td1 = zeros(size(cada1td1,1),114);
cada2td1(:,Gator2Data.Index196) = azidot.dVdV;
cada2td1(:,Gator2Data.Index197) = cada1td1dV(:,Gator2Data.Index198);
cada1td1dV = cada2td1;
cada1td1(:,Gator1Data.Index43) = azidot.dV;
phaseout.dynamics.dVdV = cada1td1dV; phaseout.dynamics.dV = cada1td1;
phaseout.dynamics.f = [raddot.f londot.f latdot.f vdot.f fpadot.f azidot.f];
%User Line: phaseout.dynamics  = [raddot, londot, latdot, vdot, fpadot, azidot];
phaseout.dynamics.dV_size = [6 9];
phaseout.dynamics.dV_location = Gator1Data.Index44;
phaseout.dynamics.dVdV_size = [phaseout.dynamics.dV_size,9];
phaseout.dynamics.dVdV_location = [phaseout.dynamics.dV_location(Gator2Data.Index199,:), Gator2Data.Index200];
end


function ADiGator_LoadData()
global ADiGator_rlvEntryContinuousADiGatorHes
ADiGator_rlvEntryContinuousADiGatorHes = load('rlvEntryContinuousADiGatorHes.mat');
return
end